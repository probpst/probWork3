!=======================================================================
! patch_flux.F90  -  VOLUME FLUX THROUGH THE PHYSICAL BOUNDARY PATCHES
!                    (heart CFD extension, update H1)
!
!   Every non-wall patch of the .bcs file (inlets, outlets, pressure
!   boundaries) is grouped by its LABEL (the comment of the direction
!   line of its record, e.g. "! outlet_descending_aorta (type 12)";
!   patches without a comment are named b<block>_<face>).  Every
!   Data_Intvl steps the corrected volume flux J U^q / fvw of the faces
!   of each group is summed (positive = OUT of the fluid domain) and
!   written to Output_data/<case>_patch_flux.dat [m3/s], one column per
!   group, then the sum of all groups.  For a closed domain with moving
!   walls the sum equals -dV/dt of the fluid volume (ALE continuity).
!=======================================================================
module patch_flux
  use mpi
  use global_vars
  use mpi_vars
  use mom_ops, only: face_range
  use sixdof, only: wall_patch_range, is_wall_patch
  implicit none
  private
  public :: patch_flux_init, patch_flux_export, patch_flux_finalize

  character(len=*), parameter :: out_dir = 'Output_data/'
  integer :: ngrp = 0, npat = 0, iu = -1
  character(len=48), allocatable :: gname(:)
  integer, allocatable :: pb(:), pn(:), pg(:)        ! patch: block, patch index, group
  logical :: on = .false.

contains

  subroutine patch_flux_init
    integer :: b, n, g, nmax
    character(len=48) :: lab
    character(len=8) :: s1, s2
    logical :: ex
    nmax = 0
    do b = 1, nblk
       nmax = nmax + blk(b)%nbc
    end do
    if (nmax == 0) return
    allocate(gname(nmax), pb(nmax), pn(nmax), pg(nmax))
    do b = 1, nblk
       do n = 1, blk(b)%nbc
          if (is_wall_patch(b, n)) cycle
          if (blk(b)%bc_flag(1,n) == 3) cycle                 ! periodic
          lab = blk(b)%bc_label(n)
          if (len_trim(lab) == 0) then
             write(s1,'(I0)') b;  write(s2,'(I0)') blk(b)%bc_dir(n)
             lab = 'b'//trim(s1)//'_d'//trim(s2)
          end if
          g = find_group(lab)
          if (g == 0) then
             ngrp = ngrp + 1;  gname(ngrp) = lab;  g = ngrp
          end if
          npat = npat + 1;  pb(npat) = b;  pn(npat) = n;  pg(npat) = g
       end do
    end do
    if (npat == 0) return
    on = .true.
    if (myrank == 0) then
       inquire(file=out_dir//'.', exist=ex)
       if (.not. ex) call execute_command_line('mkdir -p '//out_dir, wait=.true.)
       open(newunit=iu, file=out_dir//trim(casename)//'_patch_flux.dat', status=merge('unknown', 'replace', run_mode == 1), &
            position='append', action='write')
       if (run_mode /= 1) then
          write(iu,'(A)') '# probCFD volume flux through the non-wall boundary patches [m3/s], positive = out of the fluid domain'
          write(iu,'(A,I0,A,I0,A)') '# ', ngrp, ' groups (patches with the same .bcs label summed) from ', npat, ' patches; last column = sum'
          write(iu,'(A)',advance='no') '#   step      time[s]'
          do g = 1, ngrp
             write(iu,'(A)',advance='no') '  '//adjustr(gname(g)(1:max(16,len_trim(gname(g)))))
          end do
          write(iu,'(A)') '               sum'
       end if
       write(*,'(A)') ' '
       write(*,'(A,I0,A,I0,A)') '  patch fluxes     : ', ngrp, ' boundary groups from ', npat, ' non-wall patches -> '//out_dir// &
          trim(casename)//'_patch_flux.dat (positive = outflow)'
       do g = 1, ngrp
          write(*,'(A,I0,A)') '    '//trim(gname(g))//' (', count(pg(1:npat) == g), ' patches)'
       end do
    end if
  contains
    integer function find_group(l)
      character(len=*), intent(in) :: l
      integer :: i
      find_group = 0
      do i = 1, ngrp
         if (gname(i) == l) find_group = i
      end do
    end function find_group
  end subroutine patch_flux_init

  subroutine patch_flux_export(istep, t)
    integer,  intent(in) :: istep
    real(dp), intent(in) :: t
    integer :: ib, d, b, ip, n, ax, sd, t1, t2, rr(2,3), c1s, c1e, c2s, c2e, m1, m2, idx(3), H(3), ierr, g
    real(dp) :: sgn, fl
    real(dp), allocatable :: q(:)
    real(dp), pointer :: Uq(:,:,:)
    if (.not. on) return
    if (mod(istep, data_interval) /= 0 .and. istep /= 0) return
    allocate(q(ngrp));  q = 0.0_dp
    do ib = 1, nown
       d = own_list(ib);  b = dom(d)%b
       call set_block(d)
       do ip = 1, npat
          if (pb(ip) /= b) cycle
          n = pn(ip)
          if (.not. wall_patch_range(d, b, n, ax, sd, t1, t2, rr, c1s, c1e, c2s, c2e)) cycle
          select case (ax)
          case (1);  Uq => U1
          case (2);  Uq => U2
          case default;  Uq => U3
          end select
          sgn = merge(-1.0_dp, 1.0_dp, sd == 1)                  ! outward: -U at a min face, +U at a max face
          do m2 = c2s, c2e
             do m1 = c1s, c1e
                idx(t1) = m1;  idx(t2) = m2
                idx(ax) = merge(0, ncell(b, ax), sd == 1)       ! face index of the family (0 .. ncell)
                if (fstat(ax, idx(1), idx(2), idx(3)) /= 1) cycle
                H = 2*idx;  H(ax) = H(ax) + 1
                fl = sgn*Jh(H(1),H(2),H(3))/fvw(H(1),H(2),H(3))*Uq(idx(1), idx(2), idx(3))
                q(pg(ip)) = q(pg(ip)) + fl
             end do
          end do
       end do
    end do
    call MPI_Allreduce(MPI_IN_PLACE, q, ngrp, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    q = q*Uch*Lch*Lch                                             ! non-dim flux -> m3/s
    if (myrank == 0) then
       write(iu,'(I8,ES13.5)',advance='no') istep, t
       do g = 1, ngrp
          write(iu,'(ES18.8)',advance='no') q(g)
       end do
       write(iu,'(ES18.8)') sum(q)
       flush(iu)
    end if
    deallocate(q)
  end subroutine patch_flux_export

  subroutine patch_flux_finalize
    if (on .and. myrank == 0 .and. iu > 0) close(iu)
  end subroutine patch_flux_finalize

end module patch_flux
