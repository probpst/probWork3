!=======================================================================
! mod_mpi.F90
!   MPI data and helpers for probCFD (multi-block + intra-block split).
!
!   Every grid block b is split into blk(b)%sp(1) x sp(2) x sp(3)
!   sub-domains (from Input_polyMesh_probCFD.csv).  The total CPU
!   count prescribed by the csv is  ntot = sum over blocks of
!   blk(b)%ncpu.  When nprocs == ntot every rank owns exactly one
!   sub-domain, allocated block by block in csv column order (block 1
!   gets ranks 0..ncpu(1)-1, block 2 the next ncpu(2) ranks, ...).
!   With nprocs < ntot the sub-domains are distributed round-robin;
!   nprocs > ntot is an error.
!=======================================================================
module mpi_vars
  use mpi
  use global_vars
  implicit none
  save

  integer :: nprocs, myrank

contains

  !---------------------------------------------------------------------
  subroutine mpi_setup
    integer :: ierr
    call MPI_Init(ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)
  end subroutine mpi_setup

  !---------------------------------------------------------------------
  ! Build the sub-domain table from the per-block split counts and
  ! assign the CPUs.  Called after read_grid + read_split_csv.
  !---------------------------------------------------------------------
  subroutine mpi_assign_domains
    integer :: b, d, ierr, ntot, pi, pj, pk, a
    integer :: nc(3), q1(3), q2(3)

    !----- count and allocate -------------------------------------------
    ntot = 0
    do b = 1, nblk
       do a = 1, 3
          if (blk(b)%sp(a) < 1 .or. blk(b)%sp(a) > ncell(b,a)/2) then
             if (myrank == 0) write(*,'(A,I0,A,I0,A,I0,A)')              &
                ' ERROR: block ', b, ' split ', blk(b)%sp(a),            &
                ' in direction ', a, ' is invalid (each part needs'//    &
                ' at least 2 cells).'
             call MPI_Abort(MPI_COMM_WORLD, 30, ierr)
          end if
       end do
       blk(b)%dom0 = ntot + 1
       ntot = ntot + blk(b)%ncpu
    end do
    ndom = ntot
    allocate(dom(ndom))

    !----- sub-domain cell ranges (balanced split) ----------------------
    d = 0
    do b = 1, nblk
       nc = [blk(b)%nx, blk(b)%ny, blk(b)%nz]
       do pk = 0, blk(b)%sp(3)-1
        do pj = 0, blk(b)%sp(2)-1
         do pi = 0, blk(b)%sp(1)-1
            d = d + 1
            dom(d)%b  = b
            dom(d)%pc = [pi, pj, pk]
            do a = 1, 3
               call split_range(nc(a), blk(b)%sp(a), dom(d)%pc(a),       &
                                q1(a), q2(a))
            end do
            dom(d)%r1 = q1;  dom(d)%r2 = q2
         end do
        end do
       end do
    end do

    !----- CPU assignment ------------------------------------------------
    if (nprocs > ndom) then
       if (myrank == 0) write(*,'(A,I0,A,I0,A)')                          &
          ' ERROR: nprocs = ', nprocs, ' exceeds the number of sub-'//    &
          'domains prescribed by the csv file (', ndom, ').'
       call MPI_Abort(MPI_COMM_WORLD, 31, ierr)
    end if
    do d = 1, ndom
       if (nprocs == ndom) then
          dom(d)%owner = d - 1              ! one CPU per sub-domain,
       else                                 ! block by block (csv order)
          dom(d)%owner = mod(d-1, nprocs)   ! fallback: round-robin
       end if
    end do

    nown = count(dom(1:ndom)%owner == myrank)
    allocate(own_list(nown))
    nown = 0
    do d = 1, ndom
       if (dom(d)%owner == myrank) then
          nown = nown + 1
          own_list(nown) = d
       end if
    end do

    !----- locate the sub-domain holding the reference pressure cell ----
    if (p_ref_switch) then
       do d = 1, ndom
          if (dom(d)%b /= p_ref_blk) cycle
          if (p_ref_i < dom(d)%r1(1) .or. p_ref_i > dom(d)%r2(1)) cycle
          if (p_ref_j < dom(d)%r1(2) .or. p_ref_j > dom(d)%r2(2)) cycle
          if (p_ref_k < dom(d)%r1(3) .or. p_ref_k > dom(d)%r2(3)) cycle
          p_ref_dom = d
          exit
       end do
    end if

    if (myrank == 0) then
       write(*,'(A)') ' '
       write(*,'(A)') ' ---------------- block splitting / CPU allocation ----------------'
       do d = 1, ndom
          write(*,'(A,I4,A,I3,A,I3,A,I3,A,I3,A,I4,A,I4,A,I4,A,I4,A,I4,A,I4)')  &
             ' PC_#:', dom(d)%owner, '     BLK_#:', dom(d)%b,             &
             '     i_th:', dom(d)%pc(1)+1, ' j_th:', dom(d)%pc(2)+1,      &
             ' k_th:', dom(d)%pc(3)+1,                                     &
             '     i-th:', dom(d)%r1(1), '~', dom(d)%r2(1)+1,              &
             ' j-th:', dom(d)%r1(2), '~', dom(d)%r2(2)+1,                  &
             ' k-th:', dom(d)%r1(3), '~', dom(d)%r2(3)+1
       end do
       write(*,'(A,I0,A,I0,A)') '  MPI             : ', nprocs,          &
            ' rank(s) for ', ndom, ' sub-domains'
       if (nprocs == ndom) then
          write(*,'(A)') '  CPU allocation  : one CPU per sub-domain'//  &
               ' (as prescribed by the csv file)'
       else
          write(*,'(A,I0,A)') '  NOTE: the csv file prescribes ', ndom,  &
               ' CPUs; running with fewer ranks -> sub-domains are'//    &
               ' distributed round-robin'
       end if
    end if
  end subroutine mpi_assign_domains

  !---------------------------------------------------------------------
  ! Balanced 1-D split of n cells into np parts; part p (0-based) gets
  ! range [q1, q2].
  !---------------------------------------------------------------------
  subroutine split_range(n, np, p, q1, q2)
    integer, intent(in)  :: n, np, p
    integer, intent(out) :: q1, q2
    integer :: base, rem
    base = n / np
    rem  = mod(n, np)
    q1 = p*base + min(p, rem) + 1
    q2 = q1 + base - 1
    if (p < rem) q2 = q2 + 1
  end subroutine split_range

  !---------------------------------------------------------------------
  real(dp) function global_max(x)
    real(dp), intent(in) :: x
    integer :: ierr
    call MPI_Allreduce(x, global_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX,  &
                       MPI_COMM_WORLD, ierr)
  end function global_max

  !---------------------------------------------------------------------
  real(dp) function global_sum(x)
    real(dp), intent(in) :: x
    integer :: ierr
    call MPI_Allreduce(x, global_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM,  &
                       MPI_COMM_WORLD, ierr)
  end function global_sum

end module mpi_vars
