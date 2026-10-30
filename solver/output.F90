!=======================================================================
! output.F90
!   PLOT3D BINARY output (Fortran unformatted), MULTI-block with
!   intra-block MPI sub-domains.
!
!   write_grid_plot3d : writes 'casename.x' once (node coords, meters)
!   write_solution    : writes 'casename_NNNNNNNN.q', a PLOT3D function
!                       file with nvar variables interpolated to the
!                       grid nodes of every block and re-dimensionalized:
!                         1: u [m/s]   2: v [m/s]   3: w [m/s]   4: p [Pa]
!                         5: alpha (liquid volume fraction, 0..1) -
!                            written only for Num_of_Phase = 2 (nvar = 5)
!                       Velocity components are GLOBAL Cartesian on
!                       every block.
!
!   Each rank interpolates the node sub-range of every sub-domain it
!   owns (the ghost values carry the boundary conditions, the halo and
!   the b2b donor data, so wall, sub-domain and interface nodes come
!   out right) and sends the node piece to rank 0, which assembles the
!   per-block node arrays.  Nodes shared by neighboring sub-domains
!   receive consistent values from either side.
!
!   SINGLE-VALUED NODES ACROSS BLOCKS: a node on a block-to-block
!   interface exists in every block sharing it.  Away from junction
!   lines the blocks interpolate it from the same (exchanged) data, but
!   where three or more blocks meet (O/H junction lines of the cylinder
!   grid, window edges of partial-face interfaces) the ghost cells of a
!   block are geometrically ambiguous and the per-block node values can
!   differ.  Rank 0 therefore averages, once per output, every field
!   over all block copies of each coincident boundary node (coincidence
!   groups built once from the node coordinates with the b2b tolerance,
!   module node_groups), so every node of the PLOT3D file carries one
!   value regardless of the block it is read from.
!=======================================================================
module node_groups
  use global_vars
  implicit none
  integer :: ngrp = 0                          ! coincidence groups (size > 1)
  integer, allocatable :: gstart(:)            ! group g = members gstart(g):gstart(g+1)-1
  integer, allocatable :: gblk(:), gi(:), gj(:), gk(:)   ! member node (block, i, j, k)
  logical :: groups_built = .false.
contains
  !--------------------------------------------------------------------
  ! Build the coincidence groups of the block-boundary nodes (rank 0).
  !--------------------------------------------------------------------
  subroutine build_node_groups
    integer :: b, i, j, k, n, nb, m, a, c, ra, rc, g
    integer, allocatable :: bb(:), ii(:), jj(:), kk(:), ord(:), root(:), cnt(:), gid(:), fill(:)
    real(dp), allocatable :: xx(:), yy(:), zz(:)
    real(dp) :: tol
    logical :: onb
    tol = max(b2b_tol, 1.0e-12_dp)
    ! count boundary nodes
    nb = 0
    do b = 1, nblk
       do k = 1, blk(b)%nk
        do j = 1, blk(b)%nj
         do i = 1, blk(b)%ni
            onb = (i == 1 .or. i == blk(b)%ni .or. j == 1 .or. j == blk(b)%nj .or. &
                   k == 1 .or. k == blk(b)%nk)
            if (onb) nb = nb + 1
         end do
        end do
       end do
    end do
    allocate(bb(nb), ii(nb), jj(nb), kk(nb), xx(nb), yy(nb), zz(nb), ord(nb), root(nb))
    n = 0
    do b = 1, nblk
       do k = 1, blk(b)%nk
        do j = 1, blk(b)%nj
         do i = 1, blk(b)%ni
            onb = (i == 1 .or. i == blk(b)%ni .or. j == 1 .or. j == blk(b)%nj .or. &
                   k == 1 .or. k == blk(b)%nk)
            if (.not. onb) cycle
            n = n + 1
            bb(n) = b;  ii(n) = i;  jj(n) = j;  kk(n) = k
            xx(n) = blk(b)%xn(i,j,k);  yy(n) = blk(b)%yn(i,j,k);  zz(n) = blk(b)%zn(i,j,k)
         end do
        end do
       end do
    end do
    ! sort by x (heap sort of the index array)
    do n = 1, nb
       ord(n) = n
    end do
    call heapsort(nb, xx, ord)
    ! union-find over x-windows of width tol
    do n = 1, nb
       root(n) = n
    end do
    do m = 1, nb
       a = ord(m)
       do c = m+1, nb
          if (xx(ord(c)) - xx(a) > tol) exit
          if (abs(yy(ord(c)) - yy(a)) <= tol .and. abs(zz(ord(c)) - zz(a)) <= tol) then
             ra = find(a);  rc = find(ord(c))
             if (ra /= rc) root(rc) = ra
          end if
       end do
    end do
    ! groups with more than one member
    allocate(cnt(nb), gid(nb));  cnt = 0
    do n = 1, nb
       ra = find(n);  cnt(ra) = cnt(ra) + 1
    end do
    ngrp = 0;  gid = 0
    do n = 1, nb
       if (cnt(n) > 1) then
          ngrp = ngrp + 1;  gid(n) = ngrp
       end if
    end do
    allocate(gstart(ngrp+1));  gstart = 0
    do n = 1, nb
       if (cnt(n) > 1) gstart(gid(n)+1) = cnt(n)
    end do
    gstart(1) = 1
    do g = 1, ngrp
       gstart(g+1) = gstart(g) + gstart(g+1)
    end do
    m = gstart(ngrp+1) - 1
    allocate(gblk(m), gi(m), gj(m), gk(m), fill(ngrp))
    fill = 0
    do n = 1, nb
       ra = find(n)
       g = gid(ra)
       if (g == 0) cycle
       c = gstart(g) + fill(g);  fill(g) = fill(g) + 1
       gblk(c) = bb(n);  gi(c) = ii(n);  gj(c) = jj(n);  gk(c) = kk(n)
    end do
    groups_built = .true.
    write(*,'(A,I0,A,I0,A)') '  node averaging  : ', ngrp,               &
       ' coincident block-boundary node groups (', m, ' node copies)'
  contains
    recursive integer function find(x) result(r)
      integer, intent(in) :: x
      r = x
      do while (root(r) /= r)
         root(r) = root(root(r))
         r = root(r)
      end do
    end function find
  end subroutine build_node_groups

  subroutine heapsort(n, key, idx)
    integer,  intent(in)    :: n
    real(dp), intent(in)    :: key(n)
    integer,  intent(inout) :: idx(n)
    integer :: i, t
    do i = n/2, 1, -1
       call sift(i, n)
    end do
    do i = n, 2, -1
       t = idx(1);  idx(1) = idx(i);  idx(i) = t
       call sift(1, i-1)
    end do
  contains
    subroutine sift(s, e)
      integer, intent(in) :: s, e
      integer :: r, c, t
      r = s
      do
         c = 2*r
         if (c > e) exit
         if (c < e) then
            if (key(idx(c+1)) > key(idx(c))) c = c + 1
         end if
         if (key(idx(c)) > key(idx(r))) then
            t = idx(r);  idx(r) = idx(c);  idx(c) = t
            r = c
         else
            exit
         end if
      end do
    end subroutine sift
  end subroutine heapsort
end module node_groups


!-----------------------------------------------------------------------
subroutine write_grid_plot3d(istep)
  use global_vars
  implicit none
  integer, intent(in) :: istep
  ! Volume_solutions_grids/<case>_NNNNNNNN.grd (update 13: step number also
  ! on the initial grid, extension .grd, all volume files in the sub-folder)
  call write_grid_to(istep, vol_dir)
end subroutine write_grid_plot3d

! the grid of step istep into the sub-folder dir (rank 0; Mode = 2 uses
! Check_grid_movement/)
subroutine write_grid_to(istep, dir)
  use global_vars
  use mpi_vars
  implicit none
  integer, intent(in) :: istep
  character(len=*), intent(in) :: dir
  integer :: iu, i, j, k, b
  character(len=256) :: fname
  if (myrank /= 0) return
  call make_dir(dir)
  write(fname,'(A,A,A,I8.8,A)') trim(dir), trim(casename), '_', istep, '.grd'
  open(newunit=iu, file=trim(fname), form='unformatted', status='replace', action='write')
  write(iu) nblk
  write(iu) (blk(b)%ni, blk(b)%nj, blk(b)%nk, b=1,nblk)
  ! single precision (real*4): the format expected by the PLOT3D loaders
  ! of Tecplot / ParaView; the solver itself computes in double precision
  do b = 1, nblk
     associate(bb=>blk(b))
       write(iu) (((real(bb%xn(i,j,k),4), i=1,bb%ni), j=1,bb%nj), k=1,bb%nk), &
                 (((real(bb%yn(i,j,k),4), i=1,bb%ni), j=1,bb%nj), k=1,bb%nk), &
                 (((real(bb%zn(i,j,k),4), i=1,bb%ni), j=1,bb%nj), k=1,bb%nk)
     end associate
  end do
  close(iu)
  if (istep == 0) write(*,'(A)') '  wrote binary PLOT3D grid: '//trim(fname)
end subroutine write_grid_to

!-----------------------------------------------------------------------
! the output sub-folder of the volume files (rank 0, once)
!-----------------------------------------------------------------------
subroutine make_output_dir
  use global_vars
  use mpi_vars
  implicit none
  logical, save :: done = .false.
  if (done .or. myrank /= 0) return
  call make_dir(vol_dir)
  done = .true.
end subroutine make_output_dir

! create a sub-folder (rank 0) when it does not exist
subroutine make_dir(path)
  use mpi_vars
  implicit none
  character(len=*), intent(in) :: path
  logical :: ex
  if (myrank /= 0) return
  inquire(file=trim(path)//'.', exist=ex)
  if (.not. ex) call execute_command_line('mkdir -p '//trim(path), wait=.true.)
end subroutine make_dir

!-----------------------------------------------------------------------
subroutine write_solution(istep)
  use global_vars
  use mpi_vars
  use node_groups
  use overset, only: ovs_node_iblank
  use ovs_geom, only: plic_normal27, plic_plane_const, plic_box_frac
  use surface_output, only: wall_node_data, cell_shear, wall_node_shear
  implicit none
  integer, intent(in) :: istep

  type :: nodearr_t
     real(dp), allocatable :: f(:,:,:,:)
  end type nodearr_t
  type(nodearr_t), allocatable :: fb(:)
  real(dp), allocatable :: piece(:,:,:,:)
  real(dp) :: avg(12)
  integer :: iu, b, d, ib, i, j, k, l, ncnt, ierr, nvar, g, c, nvf, ivd, ivb, ivv, ivt, ivs, ord(12), nord
  real(dp) :: tauw
  integer, allocatable :: ibl(:,:,:)
  integer :: n1(3), n2(3)
  integer :: stat(MPI_STATUS_SIZE)
  character(len=256) :: fname

  nvar = merge(5, 4, Num_of_Phase == 2)
  ivv  = 0
  nvar = nvar + 1;  ivs = nvar             ! + shear stress magnitude [Pa] (update 22), cell field at the nodes
  if (cav_on) then                         ! + vapour volume fraction (update 15)
     nvar = nvar + 1;  ivv = nvar
  end if
  ivt  = 0
  if (energy_on) then                      ! + temperature [K] (update 16)
     nvar = nvar + 1;  ivt = nvar
  end if
  nvf  = nvar                              ! variables in the file
  ivd  = 0
  if (Num_of_Phase == 2) then              ! + dynamic pressure p - p_hyd (update 13)
     nvf = nvf + 1;  ivd = nvf
  end if
  ivb = 0
  if (ovst_on) then                        ! + iblank (1 active, 0 fringe, -1 hole)
     nvf = nvf + 1;  ivb = nvf
  end if
  call wall_node_data(istep)               ! wall-face stresses -> wall-node |tau_w| (update 22; all ranks take part)
  if (myrank == 0) then
     if (.not. groups_built) call build_node_groups
     call make_output_dir
     allocate(fb(nblk))
     do b = 1, nblk
        allocate(fb(b)%f(blk(b)%ni, blk(b)%nj, blk(b)%nk, nvf))
        do d = blk(b)%dom0, blk(b)%dom0 + blk(b)%ncpu - 1
           call node_range(d, n1, n2)
           allocate(piece(n1(1):n2(1), n1(2):n2(2), n1(3):n2(3), nvar))
           if (dom(d)%owner == 0) then
              call node_interp(d, n1, n2, piece)
           else
              ncnt = (n2(1)-n1(1)+1)*(n2(2)-n1(2)+1)*(n2(3)-n1(3)+1)*nvar
              call MPI_Recv(piece, ncnt, MPI_DOUBLE_PRECISION,           &
                            dom(d)%owner, 400+d, MPI_COMM_WORLD, stat, ierr)
           end if
           fb(b)%f(n1(1):n2(1), n1(2):n2(2), n1(3):n2(3), 1:nvar) = piece
           deallocate(piece)
        end do
     end do
     ! solid-wall nodes: the friction shear stress magnitude of the wall faces (update 22), the same
     ! node average that integrates to the friction force of the force files
     do b = 1, nblk
        do k = 1, blk(b)%nk
         do j = 1, blk(b)%nj
          do i = 1, blk(b)%ni
             if (wall_node_shear(b, i, j, k, tauw)) fb(b)%f(i,j,k,ivs) = tauw
          end do
         end do
        end do
     end do
     ! single-valued nodes: average every coincidence group over its
     ! block copies and write the average back to all of them
     do g = 1, ngrp
        avg = 0.0_dp
        do c = gstart(g), gstart(g+1)-1
           avg(1:nvar) = avg(1:nvar) + fb(gblk(c))%f(gi(c), gj(c), gk(c), 1:nvar)
        end do
        avg = avg / real(gstart(g+1) - gstart(g), dp)
        do c = gstart(g), gstart(g+1)-1
           fb(gblk(c))%f(gi(c), gj(c), gk(c), 1:nvar) = avg(1:nvar)
        end do
     end do
     ! two-phase: dynamic pressure = p - p_hyd, p_hyd(z) = rho_ref(z) |g| (z_fs - z)
     ! (rho_ref = water below the still free surface, air above; the same
     ! still-water hydrostatic part that the overset interpolation subtracts)
     if (ivd > 0) then
        do b = 1, nblk
           do k = 1, blk(b)%nk
            do j = 1, blk(b)%nj
             do i = 1, blk(b)%ni
                fb(b)%f(i,j,k,ivd) = fb(b)%f(i,j,k,4) - phyd_dim(blk(b)%zn(i,j,k)) - merge(p_level, 0.0_dp, cav_on)
             end do
            end do
           end do
        end do
     end if
     ! overset: node iblank as the last variable (not averaged)
     if (ivb > 0) then
        do b = 1, nblk
           allocate(ibl(blk(b)%ni, blk(b)%nj, blk(b)%nk))
           call ovs_node_iblank(b, ibl)
           fb(b)%f(:,:,:,ivb) = real(ibl, dp)
           deallocate(ibl)
        end do
     end if
     write(fname,'(A,A,A,I8.8,A)') trim(vol_dir), trim(casename), '_', istep, '.p3d'
     open(newunit=iu, file=trim(fname), form='unformatted',              &
          status='replace', action='write')
     write(iu) nblk
     write(iu) (blk(b)%ni, blk(b)%nj, blk(b)%nk, nvf, b=1,nblk)
     ! variable order in the file (update 13b): u, v, w, TTL_Prs, Shear_Stress (update 22), [Dym_Prs, VoF], [Blank]
     nord = 0
     do l = 1, 4
        nord = nord + 1;  ord(nord) = l
     end do
     nord = nord + 1;  ord(nord) = ivs               ! Shear_Stress [Pa] (update 22), directly after the pressure
     if (ivd > 0) then
        nord = nord + 1;  ord(nord) = ivd            ! Dym_Prs
        nord = nord + 1;  ord(nord) = 5              ! VoF
     end if
     if (ivv > 0) then
        nord = nord + 1;  ord(nord) = ivv            ! Vapor (update 15)
     end if
     if (ivt > 0) then
        nord = nord + 1;  ord(nord) = ivt            ! Temperature (update 16)
     end if
     if (ivb > 0) then
        nord = nord + 1;  ord(nord) = ivb            ! Blank
     end if
     do b = 1, nblk
        write(iu) ((((real(fb(b)%f(i,j,k,ord(l)),4), i=1,blk(b)%ni),          &
                     j=1,blk(b)%nj), k=1,blk(b)%nk), l=1,nvf)
        deallocate(fb(b)%f)
     end do
     deallocate(fb)
     close(iu)
     write(*,'(A)') '  wrote solution: '//trim(fname)
     if (istep == 0 .or. run_mode == 1) then
        if (Num_of_Phase == 2) then
           write(*,'(A)') '  variables: u, v, w [m/s], TTL_Prs [Pa], Shear_Stress [Pa], Dym_Prs = p - p_hyd [Pa], VoF'// &
              merge(', Vapor', '       ', cav_on)// &
              merge(', Temperature [K]', '                 ', energy_on)//merge(', Blank', '       ', ovst_on)
        else
           write(*,'(A)') '  variables: u, v, w [m/s], TTL_Prs [Pa], Shear_Stress [Pa]'//merge(', Vapor', '       ', cav_on)// &
              merge(', Temperature [K]', '                 ', energy_on)//merge(', Blank', '       ', ovst_on)
        end if
        write(*,'(A)') '  Shear_Stress (update 22): mu_eff sqrt(2 S_ij S_ij) at the cells (laminar + eddy viscosity), '// &
           'the friction stress |tau_w| of the wall faces at the solid-wall nodes'
        ! Tecplot PLOT3D name file: the variable names of the function file
        ! (load it with the grid and the solution: "Name file" in the
        ! PLOT3D loader; without it Tecplot shows F1V1, F1V2, ...)
        open(newunit=iu, file=trim(vol_dir)//trim(casename)//'.nam', status='replace', action='write')
        ! (update 23: the SI unit follows the name, e.g. "u [m/s]")
        write(iu,'(A)') 'u [m/s]'
        write(iu,'(A)') 'v [m/s]'
        write(iu,'(A)') 'w [m/s]'
        write(iu,'(A)') 'TTL_Prs [Pa]'
        write(iu,'(A)') 'Shear_Stress [Pa]'
        if (Num_of_Phase == 2) then
           write(iu,'(A)') 'Dym_Prs [Pa]'
           write(iu,'(A)') 'VoF [-]'
        end if
        if (cav_on) write(iu,'(A)') 'Vapor [-]'
        if (energy_on) write(iu,'(A)') 'Temperature [K]'
        if (ovst_on) write(iu,'(A)') 'Blank [-]'
        close(iu)
        write(*,'(A)') '  variable names for Tecplot: '//trim(vol_dir)//trim(casename)//'.nam (PLOT3D name file)'
     end if
  else
     do ib = 1, nown
        d = own_list(ib)
        call node_range(d, n1, n2)
        allocate(piece(n1(1):n2(1), n1(2):n2(2), n1(3):n2(3), nvar))
        call node_interp(d, n1, n2, piece)
        ncnt = (n2(1)-n1(1)+1)*(n2(2)-n1(2)+1)*(n2(3)-n1(3)+1)*nvar
        call MPI_Send(piece, ncnt, MPI_DOUBLE_PRECISION, 0, 400+d,       &
                      MPI_COMM_WORLD, ierr)
        deallocate(piece)
     end do
  end if

contains
  !--------------------------------------------------------------------
  ! Node index range covered by sub-domain d: cells r1..r2 span the
  ! nodes r1 .. r2+1 in every direction.
  !--------------------------------------------------------------------
  subroutine node_range(d_, n1_, n2_)
    integer, intent(in)  :: d_
    integer, intent(out) :: n1_(3), n2_(3)
    n1_ = dom(d_)%r1
    n2_ = dom(d_)%r2 + 1
  end subroutine node_range

  !--------------------------------------------------------------------
  ! Interpolate the staggered fields of (owned) sub-domain d to its
  ! nodes and re-dimensionalize.  Ghost values (indices r1-1 / r2+1)
  ! carry the BCs, the halo and the b2b donor data.
  !--------------------------------------------------------------------
  ! still-water hydrostatic pressure [Pa] at the dimensional height z
  real(dp) function phyd_dim(z)
    real(dp), intent(in) :: z
    real(dp) :: zn
    zn = z/Lch
    phyd_dim = -gnd(3)*(zfs - zn)*merge(1.0_dp, rho_ratio, zn <= zfs) * (rho1*Uch*Uch)
  end function phyd_dim

  subroutine node_interp(d_, n1_, n2_, fn)
    integer,  intent(in)  :: d_, n1_(3), n2_(3)
    real(dp), intent(out) :: fn(n1_(1):n2_(1), n1_(2):n2_(2),            &
                                n1_(3):n2_(3), nvar)
    integer :: im, jm, km
    call set_block(d_)
    do k = n1_(3), n2_(3)
     do j = n1_(2), n2_(2)
      do i = n1_(1), n2_(1)
         im = i-1;  jm = j-1;  km = k-1        ! face index of node (i,j,k)
         if (ovst_on) then
            ! overset: averages over the non-hole faces / cells only
            fn(i,j,k,1) = favg(1, im, jm, km) * Uch
            fn(i,j,k,2) = favg(2, im, jm, km) * Uch
            fn(i,j,k,3) = favg(3, im, jm, km) * Uch
            fn(i,j,k,4) = cavg(p, im, jm, km) * (rho1*Uch*Uch)
            if (Num_of_Phase == 2) fn(i,j,k,5) = merge(anode_plic(im, jm, km), cavg(alpha, im, jm, km), vof_scheme == 3)
            if (ivv > 0) fn(i,j,k,ivv) = cavg(vap, im, jm, km)
            if (ivt > 0) fn(i,j,k,ivt) = cavg(temp, im, jm, km)
            fn(i,j,k,ivs) = shear_node(im, jm, km)
            cycle
         end if
         fn(i,j,k,ivs) = shear_node(im, jm, km)
         ! u lives at (xf, yc, zc): average over 4 neighbors in j,k
         fn(i,j,k,1) = 0.25_dp*( u(im,jm,km) + u(im,jm+1,km)             &
                               + u(im,jm,km+1) + u(im,jm+1,km+1) ) * Uch
         ! v lives at (xc, yf, zc)
         fn(i,j,k,2) = 0.25_dp*( v(im,jm,km) + v(im+1,jm,km)             &
                               + v(im,jm,km+1) + v(im+1,jm,km+1) ) * Uch
         ! w lives at (xc, yc, zf)
         fn(i,j,k,3) = 0.25_dp*( w(im,jm,km) + w(im+1,jm,km)             &
                               + w(im,jm+1,km) + w(im+1,jm+1,km) ) * Uch
         ! p lives at (xc, yc, zc): average over 8 neighbors
         fn(i,j,k,4) = 0.125_dp*( p(im,jm,km)     + p(im+1,jm,km)        &
                                + p(im,jm+1,km)   + p(im+1,jm+1,km)      &
                                + p(im,jm,km+1)   + p(im+1,jm,km+1)      &
                                + p(im,jm+1,km+1) + p(im+1,jm+1,km+1) )  &
                       * (rho1*Uch*Uch)
         if (ivv > 0) fn(i,j,k,ivv) = 0.125_dp*( vap(im,jm,km) + vap(im+1,jm,km) + vap(im,jm+1,km) + vap(im+1,jm+1,km) &
                + vap(im,jm,km+1) + vap(im+1,jm,km+1) + vap(im,jm+1,km+1) + vap(im+1,jm+1,km+1) )
         if (ivt > 0) fn(i,j,k,ivt) = 0.125_dp*( temp(im,jm,km) + temp(im+1,jm,km) + temp(im,jm+1,km) + temp(im+1,jm+1,km) &
                + temp(im,jm,km+1) + temp(im+1,jm,km+1) + temp(im,jm+1,km+1) + temp(im+1,jm+1,km+1) )
         ! alpha lives at (xc, yc, zc): average over 8 neighbors
         if (Num_of_Phase == 2) then
            if (vof_scheme == 3) then
               fn(i,j,k,5) = anode_plic(im, jm, km)
            else
               fn(i,j,k,5) = 0.125_dp*( alpha(im,jm,km)                   &
                + alpha(im+1,jm,km)   + alpha(im,jm+1,km)                 &
                + alpha(im+1,jm+1,km) + alpha(im,jm,km+1)                 &
                + alpha(im+1,jm,km+1) + alpha(im,jm+1,km+1)               &
                + alpha(im+1,jm+1,km+1) )
            end if
         end if
      end do
     end do
    end do
  end subroutine node_interp

  ! overset: Cartesian component q of the 4 faces of family q around
  ! node (im+1,jm+1,km+1), hole faces (status -1) excluded
  real(dp) function favg(q, im, jm, km)
    integer, intent(in) :: q, im, jm, km
    integer :: a, c, id(3), e1(3), e2(3), t1, t2, n
    real(dp) :: s
    real(dp), pointer :: f(:,:,:)
    select case (q)
    case (1);  f => u;  t1 = 2;  t2 = 3
    case (2);  f => v;  t1 = 1;  t2 = 3
    case default;  f => w;  t1 = 1;  t2 = 2
    end select
    e1 = 0;  e1(t1) = 1;  e2 = 0;  e2(t2) = 1
    s = 0.0_dp;  n = 0
    do a = 0, 1
       do c = 0, 1
          id = [im, jm, km] + a*e1 + c*e2
          if (fstat(q, id(1), id(2), id(3)) /= -1) then
             s = s + f(id(1),id(2),id(3));  n = n + 1
          end if
       end do
    end do
    favg = 0.0_dp
    if (n > 0) favg = s / real(n, dp)
  end function favg

  ! overset: cell field over the 8 cells around the node, holes excluded
  !--------------------------------------------------------------------
  ! GEOMETRIC VoF: node value of alpha = liquid fraction of the dual cell
  ! around the node (the 8 octants of the surrounding cells adjacent to
  ! the node), each octant evaluated on its cell's PLIC plane.  The
  ! sharp reconstructed field is then exported with a one-node band
  ! instead of the two-layer smear of the 8-cell average.
  !--------------------------------------------------------------------
  real(dp) function anode_plic(im, jm, km) result(an)
    integer, intent(in) :: im, jm, km
    integer  :: a, c, e, ci, cj, ck, d1, d2, d3, n
    real(dp) :: st(-1:1,-1:1,-1:1), m(3), acst, x0(3), w(3), ac, s, vc, vs
    logical  :: mir(3)
    s = 0.0_dp;  vs = 0.0_dp;  n = 0
    do e = 0, 1
     do c = 0, 1
      do a = 0, 1
         ci = im + a;  cj = jm + c;  ck = km + e
         if (ovst_on) then
            if (cstat(ci,cj,ck) < 0) cycle
         end if
         vc = Jh(2*ci, 2*cj, 2*ck)
         ac = alpha(ci,cj,ck)
         if (ac > 0.0_dp .and. ac < 1.0_dp .and. ci-1 >= is-ng .and. cj-1 >= js-ng .and. ck-1 >= ks-ng &
             .and. ci+1 <= ie+ng .and. cj+1 <= je+ng .and. ck+1 <= ke+ng) then
            do d3 = -1, 1
             do d2 = -1, 1
              do d1 = -1, 1
                 st(d1,d2,d3) = alpha(ci+d1, cj+d2, ck+d3)
              end do
             end do
            end do
            call plic_normal27(st, m, mir)
            acst = plic_plane_const(m, ac)
            ! octant of the cell that touches the node: x in [0,1/2] when
            ! the node is at the cell's low side (a = 1: node index = cell
            ! index -> low side), [1/2,1] otherwise; mirrored axes flip
            x0 = [merge(0.0_dp, 0.5_dp, a == 1), merge(0.0_dp, 0.5_dp, c == 1), merge(0.0_dp, 0.5_dp, e == 1)]
            do d1 = 1, 3
               if (mir(d1)) x0(d1) = 0.5_dp - x0(d1)
            end do
            w = 0.5_dp
            ac = plic_box_frac(m, acst, x0, w)
         end if
         s = s + ac*vc;  vs = vs + vc;  n = n + 1
      end do
     end do
    end do
    if (vs > 0.0_dp) then
       an = s/vs
    else
       an = 0.0_dp
    end if
    an = min(max(an, 0.0_dp), 1.0_dp)
  end function anode_plic

  ! shear stress magnitude at the node: average of the 8 cells around it whose gradient
  ! stencil lies inside the domain's ghost range (holes excluded)
  real(dp) function shear_node(im, jm, km)
    integer, intent(in) :: im, jm, km
    integer :: a, c, e, n
    real(dp) :: s
    s = 0.0_dp;  n = 0
    do e = 0, 1
     do c = 0, 1
      do a = 0, 1
         if (cstat(im+a, jm+c, km+e) < 0) cycle
         if (im+a-1 < is-ng .or. im+a+1 > ie+ng .or. jm+c-1 < js-ng .or. jm+c+1 > je+ng .or. &
             km+e-1 < ks-ng .or. km+e+1 > ke+ng) cycle
         s = s + cell_shear(im+a, jm+c, km+e);  n = n + 1
      end do
     end do
    end do
    shear_node = 0.0_dp
    if (n > 0) shear_node = s/real(n, dp)
  end function shear_node

  real(dp) function cavg(f, im, jm, km)
    real(dp), intent(in) :: f(is-ng:, js-ng:, ks-ng:)
    integer,  intent(in) :: im, jm, km
    integer :: a, c, e, n
    real(dp) :: s
    s = 0.0_dp;  n = 0
    do e = 0, 1
     do c = 0, 1
      do a = 0, 1
         if (cstat(im+a, jm+c, km+e) >= 0) then
            s = s + f(im+a, jm+c, km+e);  n = n + 1
         end if
      end do
     end do
    end do
    cavg = 0.0_dp
    if (n > 0) cavg = s / real(n, dp)
  end function cavg
end subroutine write_solution

!=======================================================================
! Overset grid-check mode: node iblank (1 active, 0 fringe, -1 hole) of
! every block as a one-variable binary PLOT3D function file
! <case>_iblank.q (rank 0), matching <case>.x.
!=======================================================================
subroutine write_iblank_plot3d
  use global_vars
  use mpi_vars
  use overset, only: ovs_node_iblank
  implicit none
  integer :: iu, b, i, j, k, n1, n2, n3
  integer, allocatable :: ibl(:,:,:)
  real(4), allocatable :: bv(:,:,:), bv2(:,:,:)
  if (myrank /= 0) return
  call make_output_dir
  ! (overset off: every node is active -> iblank = 1 everywhere)
  open(newunit=iu, file=trim(vol_dir)//trim(casename)//'_00000000_iblank.p3d', form='unformatted', status='replace', action='write')
  write(iu) nblk
  write(iu) (blk(b)%ni, blk(b)%nj, blk(b)%nk, 1, b=1,nblk)
  do b = 1, nblk
     allocate(ibl(blk(b)%ni, blk(b)%nj, blk(b)%nk))
     if (ovst_on) then
        call ovs_node_iblank(b, ibl)
     else
        ibl = 1
     end if
     ! blank value: 1 active, -1 hole, 0 fringe - except that a fringe
     ! node DIRECTLY ADJACENT (6 index neighbours) to an active node is
     ! written as 0.5 (the interpolation boundary that the active region
     ! touches)
     n1 = blk(b)%ni;  n2 = blk(b)%nj;  n3 = blk(b)%nk
     allocate(bv(n1, n2, n3))
     do k = 1, n3
      do j = 1, n2
       do i = 1, n1
          bv(i,j,k) = real(ibl(i,j,k), 4)
          if (ibl(i,j,k) /= 0) cycle
          if (i > 1) then
             if (ibl(i-1,j,k) == 1) bv(i,j,k) = 0.5
          end if
          if (i < n1) then
             if (ibl(i+1,j,k) == 1) bv(i,j,k) = 0.5
          end if
          if (j > 1) then
             if (ibl(i,j-1,k) == 1) bv(i,j,k) = 0.5
          end if
          if (j < n2) then
             if (ibl(i,j+1,k) == 1) bv(i,j,k) = 0.5
          end if
          if (k > 1) then
             if (ibl(i,j,k-1) == 1) bv(i,j,k) = 0.5
          end if
          if (k < n3) then
             if (ibl(i,j,k+1) == 1) bv(i,j,k) = 0.5
          end if
       end do
      end do
     end do
     ! second layer: a fringe node (0) directly adjacent to a 0.5 node is
     ! written as 0.25 (uses the 0.5 values set above, so the pass is
     ! made on a copy)
     allocate(bv2(n1, n2, n3));  bv2 = bv
     do k = 1, n3
      do j = 1, n2
       do i = 1, n1
          if (bv(i,j,k) /= 0.0) cycle
          if (i > 1) then
             if (bv(i-1,j,k) == 0.5) bv2(i,j,k) = 0.25
          end if
          if (i < n1) then
             if (bv(i+1,j,k) == 0.5) bv2(i,j,k) = 0.25
          end if
          if (j > 1) then
             if (bv(i,j-1,k) == 0.5) bv2(i,j,k) = 0.25
          end if
          if (j < n2) then
             if (bv(i,j+1,k) == 0.5) bv2(i,j,k) = 0.25
          end if
          if (k > 1) then
             if (bv(i,j,k-1) == 0.5) bv2(i,j,k) = 0.25
          end if
          if (k < n3) then
             if (bv(i,j,k+1) == 0.5) bv2(i,j,k) = 0.25
          end if
       end do
      end do
     end do
     ! third layer: a fringe node (0) directly adjacent to a 0.25 node is
     ! written as 0.125 (again on a copy, one node deep)
     bv = bv2
     do k = 1, n3
      do j = 1, n2
       do i = 1, n1
          if (bv2(i,j,k) /= 0.0) cycle
          if (i > 1) then
             if (bv2(i-1,j,k) == 0.25) bv(i,j,k) = 0.125
          end if
          if (i < n1) then
             if (bv2(i+1,j,k) == 0.25) bv(i,j,k) = 0.125
          end if
          if (j > 1) then
             if (bv2(i,j-1,k) == 0.25) bv(i,j,k) = 0.125
          end if
          if (j < n2) then
             if (bv2(i,j+1,k) == 0.25) bv(i,j,k) = 0.125
          end if
          if (k > 1) then
             if (bv2(i,j,k-1) == 0.25) bv(i,j,k) = 0.125
          end if
          if (k < n3) then
             if (bv2(i,j,k+1) == 0.25) bv(i,j,k) = 0.125
          end if
       end do
      end do
     end do
     write(iu) (((bv(i,j,k), i=1,n1), j=1,n2), k=1,n3)
     deallocate(ibl, bv, bv2)
  end do
  close(iu)
  write(*,'(A)') '  wrote node iblank: '//trim(vol_dir)//trim(casename)// &
     '_00000000_iblank.p3d (1 active, 0.5 / 0.25 / 0.125 fringe layers from the active region, 0 fringe, -1 hole)'
end subroutine write_iblank_plot3d
