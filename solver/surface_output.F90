!=======================================================================
! surface_output.F90 - solid-wall surface data (update 22)
!
!   * wall_node_data(istep): on every output step, the pressure and
!     friction stress of every solid-wall face (is_wall_patch: no-slip
!     flag-1 patches with zero values and wall-function flag-11 patches
!     of every block, with or without a 6-DoF body) from
!     sixdof.wall_face_stress - the SAME per-face values the force
!     integration sums - gathered on rank 0 and averaged to the patch
!     nodes with the dual-area rule
!        v_node = sum_{faces f at the node} w_f v_f A_f / 4
!               / sum_{faces f at the node} A_f / 4
!     (w_f = overset wall weight of the face, A_f its area), so that the
!     trapezoidal integration over the patch cells (mean of the four node
!     values times the cell area) reproduces sum_f w_f v_f A_f exactly:
!     the integral of the friction vector over the surface file equals
!     the friction force of the force files, that of the pressure vector
!     the pressure force (hole faces of an overset wall count with zero
!     stress and full area, as in the force integration).
!   * surf_write(istep): SW_Prs_Frc_surf = 1 - the surface files
!     Volume_solutions_grids/<case>_surf_NNNNNNNN.grd / .p3d (one PLOT3D
!     block of ni x nj x 1 nodes per wall patch) with the six variables
!     Prs_x Prs_y Prs_z Frc_x Frc_y Frc_z [Pa]: the pressure stress
!     vector -p n (n into the fluid: the pressure force per unit area ON
!     the wall) and the friction stress vector on the wall, in the
!     earth-fixed frame (Frame_Prs_Frc_surf = 0: the moved grid
!     coordinates) or in the background frame (1: components and node
!     coordinates turned back by the background pose; identical when the
!     background does not move).  <case>_surf.nam (variable names) and
!     <case>_surf_patches.txt (zone -> block, patch, body) once.
!   * cell_shear(i,j,k): shear stress magnitude of a cell [Pa] for the
!     volume files and the probes: mu_eff sqrt(2 S_ij S_ij) with the
!     strain rate from the cell-centred velocity gradient and mu_eff =
!     laminar mixture viscosity + eddy viscosity (turbulence / IDDES
!     model); for simple shear du/dy this is mu_eff du/dy.
!   * wall_node_shear(b, i, j, k, tau): the wall-node value |tau_w| of
!     the Shear_Stress variable of the volume files (the same node
!     average of the face friction magnitudes).
!=======================================================================
module surface_output
  use global_vars
  use mpi_vars
  use mpi
  use sixdof, only: sixdof_on_flag, wall_patch_range, wall_face_stress, is_wall_patch, sixdof_part_name, tang_axes
  use overset, only: ovs_wall_weights
  use turbulence, only: velocity_gradient, strain_vorticity, mu_lam
  use fsi, only: fsi_wall_node_value
  implicit none
  private
  public :: wall_node_data, surf_write, cell_shear, wall_node_shear, surf_npatch, surf_patch_of_node

  type :: spatch_t
     integer :: b = 0, n = 0, ax = 0, sd = 0, t1 = 0, t2 = 0, rr(2,3) = 0
     integer :: n1 = 0, n2 = 0                       ! face cells along t1, t2
     real(dp), allocatable :: fv(:,:,:)               ! (8, n1, n2): Pvec(3) Fvec(3) area weight  (face values)
     real(dp), allocatable :: nv(:,:,:)               ! (7, n1+1, n2+1): Pvec Fvec |tau|  (node values, rank 0)
  end type spatch_t
  type(spatch_t), allocatable :: sp(:)
  integer :: nsp = 0
  integer :: data_step = -1
  logical :: listed = .false.

contains

  integer function surf_npatch()
    surf_npatch = nsp
  end function surf_npatch

  !---------------------------------------------------------------------
  ! the wall patches of all blocks (every rank: the list is replicated)
  !---------------------------------------------------------------------
  subroutine build_list
    integer :: b, n, ax, sd, t1, t2, rr(2,3), cnt, pass
    do pass = 1, 2
       cnt = 0
       do b = 1, nblk
          do n = 1, blk(b)%nbc
             if (.not. is_wall_patch(b, n)) cycle
             ax = abs(blk(b)%bc_dir(n));  sd = merge(1, 2, blk(b)%bc_dir(n) > 0)
             rr(:,1) = blk(b)%bc_irng(:,n);  rr(:,2) = blk(b)%bc_jrng(:,n);  rr(:,3) = blk(b)%bc_krng(:,n)
             call tang_axes(ax, t1, t2)
             if (rr(2,t1) - rr(1,t1) < 1 .or. rr(2,t2) - rr(1,t2) < 1) cycle
             cnt = cnt + 1
             if (pass == 2) then
                sp(cnt)%b = b;  sp(cnt)%n = n;  sp(cnt)%ax = ax;  sp(cnt)%sd = sd
                sp(cnt)%t1 = t1;  sp(cnt)%t2 = t2;  sp(cnt)%rr = rr
                sp(cnt)%n1 = rr(2,t1) - rr(1,t1);  sp(cnt)%n2 = rr(2,t2) - rr(1,t2)
                allocate(sp(cnt)%fv(8, sp(cnt)%n1, sp(cnt)%n2))
                allocate(sp(cnt)%nv(7, sp(cnt)%n1 + 1, sp(cnt)%n2 + 1))
             end if
          end do
       end do
       if (pass == 1) then
          nsp = cnt
          allocate(sp(max(nsp, 1)))
       end if
    end do
    listed = .true.
  end subroutine build_list

  !---------------------------------------------------------------------
  ! face values of every wall patch (all ranks fill their part, the sum
  ! lands on rank 0), then the node values on rank 0.  Once per step.
  !---------------------------------------------------------------------
  subroutine wall_node_data(istep)
    integer, intent(in) :: istep
    integer :: ip, ib, d, b, n, ax, sd, t1, t2, rr(2,3), c1s, c1e, c2s, c2e, m1, m2, i1, i2, ierr, a, c
    real(dp), allocatable :: wgt(:,:)
    real(dp) :: Av(3), amag, nrm(3), xf(3), pw, tauw(3), wt, num(7), den
    logical :: ok
    if (.not. listed) call build_list
    if (istep == data_step) return
    data_step = istep
    do ip = 1, nsp
       sp(ip)%fv = 0.0_dp
    end do
    do ib = 1, nown
       d = own_list(ib);  b = dom(d)%b
       call set_block(d)
       do ip = 1, nsp
          if (sp(ip)%b /= b) cycle
          n = sp(ip)%n
          if (.not. wall_patch_range(d, b, n, ax, sd, t1, t2, rr, c1s, c1e, c2s, c2e)) cycle
          allocate(wgt(sp(ip)%n1, sp(ip)%n2))
          call ovs_wall_weights(b, n, wgt)
          do m2 = c2s, c2e
           do m1 = c1s, c1e
              i1 = m1 - rr(1,t1) + 1;  i2 = m2 - rr(1,t2) + 1
              wt = wgt(i1, i2)
              call wall_face_stress(b, n, ax, sd, t1, t2, m1, m2, ok, Av, amag, nrm, xf, pw, tauw)
              if (.not. ok) then
                 ! hole / degenerate face: no stress; the face keeps its geometric area so that the
                 ! node average integrates to the force (the force integration skips it as well)
                 sp(ip)%fv(7, i1, i2) = face_area(b, ax, sd, t1, t2, m1, m2)
                 sp(ip)%fv(8, i1, i2) = 0.0_dp
                 cycle
              end if
              sp(ip)%fv(1:3, i1, i2) = -pw*nrm
              sp(ip)%fv(4:6, i1, i2) = tauw
              sp(ip)%fv(7, i1, i2)   = amag
              sp(ip)%fv(8, i1, i2)   = max(wt, 0.0_dp)
           end do
          end do
          deallocate(wgt)
       end do
    end do
    do ip = 1, nsp
       if (myrank == 0) then
          call MPI_Reduce(MPI_IN_PLACE, sp(ip)%fv, 8*sp(ip)%n1*sp(ip)%n2, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
       else
          call MPI_Reduce(sp(ip)%fv, sp(ip)%fv, 8*sp(ip)%n1*sp(ip)%n2, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
       end if
    end do
    if (myrank /= 0) return
    ! node values: dual-area weighted average of the adjacent faces
    do ip = 1, nsp
       do i2 = 1, sp(ip)%n2 + 1
          do i1 = 1, sp(ip)%n1 + 1
             num = 0.0_dp;  den = 0.0_dp
             do c = max(i2-1, 1), min(i2, sp(ip)%n2)
                do a = max(i1-1, 1), min(i1, sp(ip)%n1)
                   wt = sp(ip)%fv(8, a, c)*sp(ip)%fv(7, a, c)*0.25_dp
                   num(1:6) = num(1:6) + wt*sp(ip)%fv(1:6, a, c)
                   num(7)   = num(7)   + wt*sqrt(sum(sp(ip)%fv(4:6, a, c)**2))
                   den = den + sp(ip)%fv(7, a, c)*0.25_dp
                end do
             end do
             if (den > 0.0_dp) then
                sp(ip)%nv(:, i1, i2) = num/den
             else
                sp(ip)%nv(:, i1, i2) = 0.0_dp
             end if
          end do
       end do
    end do
  end subroutine wall_node_data

  ! geometric area of a wall face [m2] (hole faces, for the node weights)
  real(dp) function face_area(b, ax, sd, t1, t2, m1, m2) result(a)
    integer, intent(in) :: b, ax, sd, t1, t2, m1, m2
    integer  :: nn(4,3), idx(3), q
    real(dp) :: x(3,4), Av(3)
    idx(t1) = m1;  idx(t2) = m2;  idx(ax) = merge(1, ncell(b, ax) + 1, sd == 1)
    nn(1,:) = idx
    nn(2,:) = nn(1,:);  nn(2,t1) = nn(2,t1) + 1
    nn(3,:) = nn(1,:);  nn(3,t1) = nn(3,t1) + 1;  nn(3,t2) = nn(3,t2) + 1
    nn(4,:) = nn(1,:);  nn(4,t2) = nn(4,t2) + 1
    do q = 1, 4
       x(:,q) = [blk(b)%xn(nn(q,1),nn(q,2),nn(q,3)), blk(b)%yn(nn(q,1),nn(q,2),nn(q,3)), blk(b)%zn(nn(q,1),nn(q,2),nn(q,3))]
    end do
    Av = 0.5_dp*cross(x(:,3) - x(:,1), x(:,4) - x(:,2))
    a = sqrt(sum(Av*Av))
  end function face_area

  pure function cross(a, b) result(c)
    real(dp), intent(in) :: a(3), b(3)
    real(dp) :: c(3)
    c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
  end function cross

  !---------------------------------------------------------------------
  ! patch index and node indices (i1, i2) of block node (i, j, k) of block b on a wall patch
  ! (0 when the node is on no wall patch); the first patch found wins.  Rank 0, after wall_node_data.
  !---------------------------------------------------------------------
  integer function surf_patch_of_node(b, i, j, k, i1, i2) result(ip)
    integer, intent(in)  :: b, i, j, k
    integer, intent(out) :: i1, i2
    integer :: q, idx(3), nf
    idx = [i, j, k];  i1 = 0;  i2 = 0;  ip = 0
    do q = 1, nsp
       if (sp(q)%b /= b) cycle
       nf = merge(1, ncell(b, sp(q)%ax) + 1, sp(q)%sd == 1)
       if (idx(sp(q)%ax) /= nf) cycle
       if (idx(sp(q)%t1) < sp(q)%rr(1,sp(q)%t1) .or. idx(sp(q)%t1) > sp(q)%rr(2,sp(q)%t1)) cycle
       if (idx(sp(q)%t2) < sp(q)%rr(1,sp(q)%t2) .or. idx(sp(q)%t2) > sp(q)%rr(2,sp(q)%t2)) cycle
       i1 = idx(sp(q)%t1) - sp(q)%rr(1,sp(q)%t1) + 1;  i2 = idx(sp(q)%t2) - sp(q)%rr(1,sp(q)%t2) + 1
       ip = q;  return
    end do
  end function surf_patch_of_node

  ! |tau_w| at a wall node of block b (rank 0, after wall_node_data); .false. when not a wall node
  logical function wall_node_shear(b, i, j, k, tau) result(is)
    integer,  intent(in)  :: b, i, j, k
    real(dp), intent(out) :: tau
    integer :: ip, i1, i2
    tau = 0.0_dp;  is = .false.
    if (.not. listed) return
    ip = surf_patch_of_node(b, i, j, k, i1, i2)
    if (ip == 0) return
    tau = sp(ip)%nv(7, i1, i2);  is = .true.
  end function wall_node_shear

  !---------------------------------------------------------------------
  ! shear stress magnitude of cell (i, j, k) of the current block [Pa]
  !---------------------------------------------------------------------
  real(dp) function cell_shear(i, j, k) result(tau)
    integer, intent(in) :: i, j, k
    real(dp) :: g(3,3), smag, wmag, mu_eff
    call velocity_gradient(i, j, k, g)
    call strain_vorticity(g, smag, wmag)
    mu_eff = mu_lam(i, j, k)
    if (turb_on) mu_eff = mu_eff + mut(i,j,k)
    tau = mu_eff*mu1*smag*Uch/Lch
  end function cell_shear

  !---------------------------------------------------------------------
  ! the surface files of step istep (rank 0, after wall_node_data)
  !---------------------------------------------------------------------
  subroutine surf_write(istep)
    integer, intent(in) :: istep
    integer :: iu, ip, i1, i2, l, idx(3), b, q, nv
    real(dp) :: x(3), Rt(3,3), dsp(3), vmis
    real(dp), allocatable :: xs(:,:,:), vs(:,:,:)
    character(len=64) :: pname
    character(len=256) :: fname
    logical, save :: named = .false.
    if (sw_prs_frc_surf /= 1) return
    call wall_node_data(istep)
    if (myrank /= 0) return
    if (nsp == 0) then
       if (.not. named) write(*,'(A)') '  SW_Prs_Frc_surf = 1 but the grid has no solid-wall patch: no surface files'
       named = .true.
       return
    end if
    Rt = reshape([1.0_dp,0.0_dp,0.0_dp, 0.0_dp,1.0_dp,0.0_dp, 0.0_dp,0.0_dp,1.0_dp], [3,3])
    if (frame_prs_frc_surf == 1 .and. bgr_on) Rt = transpose(bgr_Rg)
    ! grid
    write(fname,'(A,A,A,I8.8,A)') trim(vol_dir), trim(casename), '_surf_', istep, '.grd'
    open(newunit=iu, file=trim(fname), form='unformatted', status='replace', action='write')
    write(iu) nsp
    write(iu) (sp(ip)%n1 + 1, sp(ip)%n2 + 1, 1, ip=1,nsp)
    do ip = 1, nsp
       b = sp(ip)%b
       allocate(xs(3, sp(ip)%n1 + 1, sp(ip)%n2 + 1))
       do i2 = 1, sp(ip)%n2 + 1
          do i1 = 1, sp(ip)%n1 + 1
             idx(sp(ip)%t1) = sp(ip)%rr(1,sp(ip)%t1) + i1 - 1
             idx(sp(ip)%t2) = sp(ip)%rr(1,sp(ip)%t2) + i2 - 1
             idx(sp(ip)%ax) = merge(1, ncell(b, sp(ip)%ax) + 1, sp(ip)%sd == 1)
             x = [blk(b)%xn(idx(1),idx(2),idx(3)), blk(b)%yn(idx(1),idx(2),idx(3)), blk(b)%zn(idx(1),idx(2),idx(3))]
             if (frame_prs_frc_surf == 1 .and. bgr_on) x = matmul(Rt, x - bgr_O0 - bgr_dg) + bgr_O0
             xs(:, i1, i2) = x
          end do
       end do
       write(iu) ((real(xs(1,i1,i2), 4), i1=1,sp(ip)%n1+1), i2=1,sp(ip)%n2+1), &
                 ((real(xs(2,i1,i2), 4), i1=1,sp(ip)%n1+1), i2=1,sp(ip)%n2+1), &
                 ((real(xs(3,i1,i2), 4), i1=1,sp(ip)%n1+1), i2=1,sp(ip)%n2+1)
       deallocate(xs)
    end do
    close(iu)
    ! solution
    write(fname,'(A,A,A,I8.8,A)') trim(vol_dir), trim(casename), '_surf_', istep, '.p3d'
    open(newunit=iu, file=trim(fname), form='unformatted', status='replace', action='write')
    nv = 6
    if (fsi_on) nv = 10                           ! + Disp_x Disp_y Disp_z vonMises (update 24)
    write(iu) nsp
    write(iu) (sp(ip)%n1 + 1, sp(ip)%n2 + 1, 1, nv, ip=1,nsp)
    do ip = 1, nsp
       b = sp(ip)%b
       allocate(vs(nv, sp(ip)%n1 + 1, sp(ip)%n2 + 1))
       do i2 = 1, sp(ip)%n2 + 1
          do i1 = 1, sp(ip)%n1 + 1
             vs(1:3, i1, i2) = matmul(Rt, sp(ip)%nv(1:3, i1, i2))
             vs(4:6, i1, i2) = matmul(Rt, sp(ip)%nv(4:6, i1, i2))
             if (nv == 10) then
                idx(sp(ip)%t1) = sp(ip)%rr(1,sp(ip)%t1) + i1 - 1
                idx(sp(ip)%t2) = sp(ip)%rr(1,sp(ip)%t2) + i2 - 1
                idx(sp(ip)%ax) = merge(1, ncell(b, sp(ip)%ax) + 1, sp(ip)%sd == 1)
                vs(7:10, i1, i2) = 0.0_dp
                if (fsi_wall_node_value(b, idx(1), idx(2), idx(3), dsp, vmis)) then
                   vs(7:9, i1, i2) = matmul(Rt, dsp);  vs(10, i1, i2) = vmis
                end if
             end if
          end do
       end do
       write(iu) (((real(vs(l,i1,i2), 4), i1=1,sp(ip)%n1+1), i2=1,sp(ip)%n2+1), l=1,nv)
       deallocate(vs)
    end do
    close(iu)
    write(*,'(A)') '  wrote wall surface files: '//trim(fname)//' (+ .grd)'
    if (.not. named) then
       open(newunit=iu, file=trim(vol_dir)//trim(casename)//'_surf.nam', status='replace', action='write')
       write(iu,'(A)') 'Prs_x [Pa]';  write(iu,'(A)') 'Prs_y [Pa]';  write(iu,'(A)') 'Prs_z [Pa]'
       write(iu,'(A)') 'Frc_x [Pa]';  write(iu,'(A)') 'Frc_y [Pa]';  write(iu,'(A)') 'Frc_z [Pa]'
       if (fsi_on) then
          write(iu,'(A)') 'Disp_x [m]';  write(iu,'(A)') 'Disp_y [m]';  write(iu,'(A)') 'Disp_z [m]';  write(iu,'(A)') 'vonMises [Pa]'
       end if
       close(iu)
       open(newunit=iu, file=trim(vol_dir)//trim(casename)//'_surf_patches.txt', status='replace', action='write')
       write(iu,'(A)') '# wall surface files <case>_surf_NNNNNNNN.grd / .p3d (update 22): one PLOT3D block (ni x nj x 1) per solid-wall patch'
       write(iu,'(A)') '# variables [Pa]: Prs_x Prs_y Prs_z = pressure stress vector -p n on the wall (n into the fluid),'
       write(iu,'(A)') '#                 Frc_x Frc_y Frc_z = friction stress vector on the wall;'
       if (fsi_on) write(iu,'(A)') '#                 Disp_x Disp_y Disp_z [m] = deformation of the wall node from its rigid-body position,'// &
          ' vonMises [Pa] = shell stress (flexible parts, update 24; 0 on rigid walls)'
       if (frame_prs_frc_surf == 1) then
          write(iu,'(A)') '#   frame: background grid (Frame_Prs_Frc_surf = 1)'
       else
          write(iu,'(A)') '#   frame: earth-fixed (Frame_Prs_Frc_surf = 0)'
       end if
       write(iu,'(A)') '# node values = dual-area weighted averages of the face values (overset wall weights included), so that'
       write(iu,'(A)') '#   sum over the patch cells of (mean of the 4 node values x cell area) = the force of the force files exactly'
       write(iu,'(A)') '# zone  block  patch  face(+/-i,j,k)  i-range  j-range  k-range  body(0 = none)'
       do q = 1, nsp
          b = sp(q)%b
          pname = ' '
          if (blk_obj(b) > 0 .and. sixdof_on_flag()) pname = '  '//sixdof_part_name(blk_obj(b))
          write(iu,'(I5,I7,I7,I8,3(I6,A,I0),I6,A)') q, b, sp(q)%n, blk(b)%bc_dir(sp(q)%n), &
             sp(q)%rr(1,1), '-', sp(q)%rr(2,1), sp(q)%rr(1,2), '-', sp(q)%rr(2,2), sp(q)%rr(1,3), '-', sp(q)%rr(2,3), &
             blk_obj(b), trim(pname)
       end do
       close(iu)
       write(*,'(A,I0,A)') '  wall surface zones: ', nsp, ' patches, listed in '//trim(vol_dir)//trim(casename)//'_surf_patches.txt'
    end if
    named = .true.
  end subroutine surf_write
end module surface_output
