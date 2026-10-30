!=======================================================================
! fsi.F90 - FLUID-STRUCTURE INTERACTION of flexible parts (update 24)
!
!   A Body / Child_body / Sub_child_body with SW_FSI = 1 is a flexible
!   part: its closed wall surface (the no-slip patches of its blocks) is
!   turned into a shell MID-SURFACE whose thickness is measured from the
!   surface grid itself (ray pairing), the shell is solved by fsi_shell
!   (co-rotational MITC4 shell, generalized-alpha), the fluid grid of the
!   part follows the deformation (mesh morphing of the part's blocks,
!   re-assembled overset), and fluid and structure are coupled strongly
!   (sub-iterations with Aitken relaxation inside every time step).
!
!   MID-SURFACE FROM THE WALL PATCH.  Every wall patch (block face) of
!   the part is a structured sheet of n1 x n2 nodes.  The AROUND
!   direction is the periodic one (O-grid around a section) or the one
!   along which the wall normals of the two ends are opposed (C-grid);
!   the other is the SPAN.  At every span station the section's chord
!   is the principal axis of its nodes, n_up = chord x e_span, and a
!   node is UPPER (n . n_up > 0.5), LOWER (< -0.5) or CAP (leading /
!   trailing edge, |n . n_up| <= 0.5).  The upper nodes of the middle
!   station give the around-index range used at every station (nup
!   nodes).  From every upper node a ray along -n hits the lower side:
!   mid-surface node = midpoint, thickness = distance (floored by
!   FSI_Min_Thickness x max thickness; a miss falls back to the lower
!   node nearest to the ray).  Shell elements = the quads (m, s) of the
!   nup x nspan lattice; coincident nodes of different sheets of one
!   part (block cuts of a blade) are merged.  Every wall node (upper,
!   lower, cap) is attached to the mid-surface by the bilinear weights of
!   its projection (nearest node for the cap nodes beyond the edge) and
!   the rigid offset r0 = X_wall - X_node of each attached node.
!   NoFSI groups (block index ranges) clamp the mid-surface nodes whose
!   upper or lower source node lies in the range; Submat regions give
!   the elements whose nodes lie in the range their own rho, E, nu.
!   The enclosed volume V (divergence theorem over the wall faces) gives
!   Mass = Density V, or Density = Mass / V.
!
!   LOADS AND MOTION TRANSFER (energy-conserving): the wall node forces
!   (-p n + tau) A/4 of the faces around a node (overset wall weights
!   included, the same values the force files integrate) go to the
!   attached nodes a as  f_a = w_a f,  m_a = w_a (R_a r_a0) x f; the wall
!   node moves as  x = sum_a w_a (X_a + u_a + R_a r_a0).  The structure
!   lives in the part's t = 0 frame: fluid forces are turned by Rg^T,
!   displacements by Rg; gravity and the inertial loads of a moving
!   frame (-m [a_O + alpha x r + omega x (omega x r) + 2 omega x v]) are
!   added to the nodal loads.
!
!   MESH MORPHING of a part's block: along every grid line normal to the
!   wall (the column of a wall node),
!     x = X0 + d_w + (R(w(s) theta_w) - I)(X0 - X_w0)
!   with d_w, theta_w the wall node's displacement and rotation and w(s)
!   = 1 near the wall (s < 0.25 of the arc length of the line), falling
!   as a cosine to 0 at the outer boundary: the boundary-layer cells
!   move rigidly with the wall, the outer cells absorb the distortion;
!   the outer boundary (overset) follows the translation.  Then the
!   rigid pose of the part is applied as for every moving block.  The
!   metrics of the morphed sub-domains are rebuilt (curvi.metrics_domain),
!   the discrete grid velocity (x^(n+1) - x^n)/dt of every half-index
!   position feeds the ALE fluxes and the wall velocity, the stored face
!   fluxes are re-projected with the new metrics, and the momentum
!   equation carries the geometric-conservation term phi div(u_grid) on
!   the morphing blocks (convection.F90).
!
!   COUPLING (every time step, Fabbri et al. 2023 Fig. 7): predictor
!   d0 = d_n + dt (1.5 v_n - 0.5 v_n-1) -> morph -> for every inner
!   iteration of the flow (momentum + Poisson + correction from the level
!   n field): wall forces -> structural step from t_n with the new loads
!   (generalized-alpha, Newton) -> Aitken relaxation of the interface
!   displacement -> convergence  ||d~_k - d_k-1|| / ||d~_k - d_n|| <
!   FSI_Tol -> morph + overset update -> next inner iteration; at most
!   FSI_Max_Iter iterations.  FSI_Start_T holds the structure before the
!   flow has developed.  FSI_Mode = 0: quasi-static equilibrium instead
!   of the dynamic step.
!
!   FILES (Output_data/): <case>_fsi_<part>_node<n>.dat (FSI_Monitor_node:
!   displacement in the part frame, velocity, von Mises stress, every
!   Data_Intvl steps), <case>_fsi_<part>_shell.grd/.p3d (the mid-surface
!   mesh with thickness, material id, clamped flag, once), the wall
!   surface files gain Disp_x Disp_y Disp_z [m] and vonMises [Pa];
!   structural states in the restart file.
!=======================================================================
module fsi
  use mpi
  use global_vars
  use mpi_vars
  use fsi_shell, only: shell_t, shell_setup, shell_step, shell_accept, shell_static, shell_stress, shell_energy, &
                       shell_set_state_from_trial, rot_exp, rot_log, shell_mode
  use sixdof, only: sixdof_nbody, sixdof_on_flag, sixdof_fsi_def, sixdof_mass_data, sixdof_frame_kin, fsidef_t, &
                    wall_patch_range, wall_face_stress, is_wall_patch, tang_axes, sixdof_part_name, sixdof_part_stem, &
                    sixdof_remesh, MAXFG
  use overset, only: ovs_wall_weights
  use curvi, only: fill_ghost_coords, metrics_domain, cart_to_contra
  use b2b_runtime, only: refresh_uvw
  use mom_ops, only: sync_face_jacobian
  use morph_frames, only: frames_part, frames_blk, frames_morph, frames_hvel, frames_commit, frames_info
  implicit none
  private
  public :: fsi_init, fsi_predict, fsi_iterate, fsi_finish_step, fsi_export, fsi_restart_write, fsi_restart_read, &
            fsi_wall_node_value, fsi_finalize, fsi_summary, fsi_nparts, fsi_monitor, fsi_max_iterations

  character(len=*), parameter :: out_dir = 'Output_data/'

  type :: sheet_t
     integer :: b = 0, n = 0, ax = 0, sd = 0, t1 = 0, t2 = 0, rr(2,3) = 0
     integer :: n1 = 0, n2 = 0                     ! patch nodes along t1, t2
     integer :: ar = 0                             ! around direction (1 = t1, 2 = t2)
     logical :: periodic = .false.
     logical :: passive = .false.                  ! no structural nodes: the wall nodes are attached to the mid-surface of the other sheets
     integer :: nup = 0, nsp = 0                   ! upper nodes, span stations
     integer, allocatable :: iup(:)                ! (nup) around index of the upper side
     integer, allocatable :: snode(:,:)            ! (nup, nsp) structural node
     integer, allocatable :: att(:,:,:)            ! (4, n1, n2) attached structural nodes (0 = none)
     real(dp), allocatable :: wgt(:,:,:)           ! (4, n1, n2)
     real(dp), allocatable :: off(:,:,:,:)         ! (3, 4, n1, n2) offsets X_wall - X_a (t = 0 frame) [m]
     integer, allocatable :: side(:,:)             ! (n1, n2) 1 upper, -1 lower, 0 cap
     real(dp), allocatable :: dw(:,:,:), thw(:,:,:)   ! (3, n1, n2) wall node displacement [m] / rotation vector
     real(dp), allocatable :: fw(:,:,:)            ! (3, n1, n2) wall node force (grid frame) [N]
  end type sheet_t

  type :: part_t
     integer :: k = 0
     type(fsidef_t) :: def
     type(shell_t) :: sh
     integer :: nsheet = 0
     type(sheet_t), allocatable :: sheet(:)
     real(dp), allocatable :: fext(:,:)            ! (6, nn) current load (part frame) [N], [N m]
     real(dp), allocatable :: drel(:,:)            ! (6, nn) relaxed interface displacement (translation, rotation increment)
     real(dp), allocatable :: rprev(:,:)           ! (6, nn) previous Aitken residual
     real(dp), allocatable :: rhist(:,:), dhist(:,:)   ! (6 nn, MAXQN) residual / structural-output history of the step (IQN-ILS)
     real(dp), allocatable :: vprev(:,:)           ! (6, nn) discrete velocity of the step before the last (predictor)
     real(dp), allocatable :: vd(:,:)              ! (6, nn) discrete velocity of the last step: displacement / rotation increment / dt (24c)
     real(dp), allocatable :: dasm(:,:)            ! (6, nn) interface displacement at the last overset assembly
     real(dp), allocatable :: vm(:)                ! (nn) von Mises stress [Pa]
     real(dp), allocatable :: thick(:)             ! (nn) nodal thickness [m]
     integer,  allocatable :: matid(:)             ! (ne) 0 = body material, n = Submat region n
     real(dp) :: omega = 0.5_dp
     integer  :: iter = 0, nit_total = 0, nsteps = 0, nnewton = 0, nfail = 0
     real(dp) :: resid = 0.0_dp, dmax = 0.0_dp, resid0 = 0.0_dp
     logical  :: active = .false., converged = .false.
     real(dp) :: volume = 0.0_dp, mass = 0.0_dp, area = 0.0_dp
     integer  :: iu_mon(MAXFG) = 0                   ! 0 = no file (newunit gives negative units)
     integer  :: mon_sheet(MAXFG) = 0, mon_i1(MAXFG) = 0, mon_i2(MAXFG) = 0
     integer  :: nclamp = 0
     ! FSI_Test_Motion (24c)
     real(dp), allocatable :: txi(:)              ! normalized distance from the clamp of every node
     real(dp) :: tL = 0.0_dp, tes(3) = 0.0_dp, ten(3) = 0.0_dp   ! free length, span axis, bending direction
     real(dp) :: txc(3) = 0.0_dp                   ! the clamp centroid (origin of the chord line)
     real(dp), allocatable :: tphu(:,:), tpht(:)  ! mode 2: the pluck shape (translation, rotation about e_s) per unit tip deflection
     real(dp), allocatable :: tsg(:), tthg(:), th0(:) ! its chord stations, rotation per station, node offsets from the chord line
     integer,  allocatable :: tix(:)               ! station of every node
     integer  :: tns = 0
     integer, allocatable :: pgrp(:)               ! FSI_Plane_Strain: span group of every node (the nodes that differ only along the span axis)
     integer  :: ngrp = 0
  end type part_t

  type :: morph_t                                   ! per morphing block
     logical :: on = .false.
     integer :: part = 0, sheet = 0
     real(dp), allocatable :: X0(:,:,:,:)           ! (3, ni, nj, nk) original coordinates [m]
     real(dp), allocatable :: sarc(:,:,:)           ! (ni, nj, nk) arc-length fraction from the wall
     real(dp), allocatable :: xprev(:,:,:,:)        ! (3, ni, nj, nk) coordinates at the end of the previous step
     real(dp), allocatable :: vh(:,:,:,:)           ! (3, 2ni-1, 2nj-1, 2nk-1) discrete velocity, non-dim
  end type morph_t

  type :: oldm_t
     real(dp), allocatable :: Jh(:,:,:), mh(:,:,:,:)
  end type oldm_t

  type :: savn_t                                    ! level-n fields of a sub-domain (24b)
     real(dp), allocatable :: u(:,:,:), v(:,:,:), w(:,:,:), U1(:,:,:), U2(:,:,:), U3(:,:,:), p(:,:,:)
  end type savn_t

  type(part_t), allocatable, target :: pt(:)
  type(morph_t), allocatable :: mb(:)
  type(oldm_t), allocatable :: om(:)
  type(savn_t), allocatable, target :: sv(:)
  logical :: step_start = .false.                    ! the first morph of a time step (the old metrics are those of level n)
  integer :: npart = 0
  logical :: inited = .false., has_prev_vel = .false.
  real(dp) :: geo_tol = 1.0e-9_dp               ! coincidence tolerance of the current part [m]
  ! sub-iteration restart (24b): 1 = the level-n state - Cartesian field, face fluxes and (Prs_form = 2)
  ! the pressure - of every owned sub-domain (default); 2 = the level-n fluxes only (tried in 24c: the
  ! sub-iterations converged less well); 0 = keep the iterate
  integer :: fsi_restore_mode = 1
  real(dp), parameter :: s_rigid = 1.0_dp            ! rigid part of the morph blend (24b: 1 = every column moves rigidly with its wall node)
  integer, parameter :: MAXQN = 12                   ! columns kept by the interface quasi-Newton
  real(dp), parameter :: pi = 3.14159265358979323846_dp

contains

  integer function fsi_nparts()
    fsi_nparts = npart
  end function fsi_nparts

  ! the largest FSI_Max_Iter of the flexible parts (the inner-iteration count of Main)
  integer function fsi_max_iterations()
    integer :: ip
    fsi_max_iterations = 1
    do ip = 1, npart
       fsi_max_iterations = max(fsi_max_iterations, pt(ip)%def%max_iter)
    end do
  end function fsi_max_iterations

  pure function cross3(a, b) result(c)
    real(dp), intent(in) :: a(3), b(3)
    real(dp) :: c(3)
    c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
  end function cross3

  pure function unit3(a) result(c)
    real(dp), intent(in) :: a(3)
    real(dp) :: c(3), l
    l = sqrt(sum(a*a))
    c = a/max(l, 1.0e-300_dp)
  end function unit3

  subroutine fsi_abort(msg)
    character(len=*), intent(in) :: msg
    integer :: ierr
    if (myrank == 0) write(*,'(A)') ' ERROR (FSI): '//msg
    flush(6)
    call MPI_Abort(MPI_COMM_WORLD, 80, ierr)
  end subroutine fsi_abort

  ! block node indices of patch node (i1, i2) of sheet s
  pure function node_idx(s, i1, i2) result(idx)
    type(sheet_t), intent(in) :: s
    integer, intent(in) :: i1, i2
    integer :: idx(3)
    idx(s%t1) = s%rr(1,s%t1) + i1 - 1
    idx(s%t2) = s%rr(1,s%t2) + i2 - 1
    idx(s%ax) = merge(1, ncell(s%b, s%ax) + 1, s%sd == 1)
  end function node_idx

  function node_x(b, idx) result(x)
    integer, intent(in) :: b, idx(3)
    real(dp) :: x(3)
    x = [blk(b)%xn(idx(1),idx(2),idx(3)), blk(b)%yn(idx(1),idx(2),idx(3)), blk(b)%zn(idx(1),idx(2),idx(3))]
  end function node_x

  !=====================================================================
  ! Set-up (before the first grid motion, grids at t = 0)
  !=====================================================================
  subroutine fsi_init
    integer :: nb, k, ip, b, n, ax, sd, t1, t2, rr(2,3), cnt, is, ierr
    type(fsidef_t) :: def
    if (.not. sixdof_on_flag()) return
    nb = sixdof_nbody()
    npart = 0
    do k = 1, nb
       def = sixdof_fsi_def(k)
       if (def%sw == 1) npart = npart + 1
    end do
    if (npart == 0) return
    if (Time_scheme == 2) call fsi_abort('SW_FSI = 1 needs Time_scheme = 1 (the coupling loop uses the inner iterations of the step)')
    allocate(pt(npart), mb(nblk), fsi_blk(nblk), om(ndom), sv(ndom))
    fsi_blk = .false.
    ip = 0
    do k = 1, nb
       def = sixdof_fsi_def(k)
       if (def%sw /= 1) cycle
       ip = ip + 1
       pt(ip)%k = k;  pt(ip)%def = def
       if (frames_part(k) > 0) then
          ! prescribed deformation from a frame file (heart extension H1, morph_frames.F90): no structure,
          ! no coupling; every block of the part morphs (fsi_blk) and the moving-grid machinery is shared
          pt(ip)%def%start_t = huge(1.0_dp)
          pt(ip)%nsheet = 0;  allocate(pt(ip)%sheet(0))
          do b = 1, nblk
             if (blk_obj(b) == k) fsi_blk(b) = .true.
          end do
          if (any(def%mon(1,1:max(def%nmon,1)) > 0) .and. def%nmon > 0) call fsi_abort(sixdof_part_name(k)// &
             ': FSI_Monitor_node is not available with FSI_Frames')
          cycle
       end if
       if (def%E <= 0.0_dp) call fsi_abort(sixdof_part_name(k)//': SW_FSI = 1 needs Youngs_modulus > 0')
       if (def%nu < 0.0_dp .or. def%nu >= 0.5_dp) call fsi_abort(sixdof_part_name(k)//': Poisson_ratio must be 0 <= nu < 0.5')
       ! sheets: the wall patches of the part's blocks
       cnt = 0
       do b = 1, nblk
          if (blk_obj(b) /= k) cycle
          do n = 1, blk(b)%nbc
             if (.not. is_wall_patch(b, n)) cycle
             ax = abs(blk(b)%bc_dir(n));  sd = merge(1, 2, blk(b)%bc_dir(n) > 0)
             rr(:,1) = blk(b)%bc_irng(:,n);  rr(:,2) = blk(b)%bc_jrng(:,n);  rr(:,3) = blk(b)%bc_krng(:,n)
             call tang_axes(ax, t1, t2)
             if (rr(2,t1) - rr(1,t1) < 2 .or. rr(2,t2) - rr(1,t2) < 1) cycle
             cnt = cnt + 1
          end do
       end do
       if (cnt == 0) call fsi_abort(sixdof_part_name(k)//': SW_FSI = 1 but none of its blocks has a no-slip wall patch')
       allocate(pt(ip)%sheet(cnt))
       cnt = 0
       do b = 1, nblk
          if (blk_obj(b) /= k) cycle
          do n = 1, blk(b)%nbc
             if (.not. is_wall_patch(b, n)) cycle
             ax = abs(blk(b)%bc_dir(n));  sd = merge(1, 2, blk(b)%bc_dir(n) > 0)
             rr(:,1) = blk(b)%bc_irng(:,n);  rr(:,2) = blk(b)%bc_jrng(:,n);  rr(:,3) = blk(b)%bc_krng(:,n)
             call tang_axes(ax, t1, t2)
             if (rr(2,t1) - rr(1,t1) < 2 .or. rr(2,t2) - rr(1,t2) < 1) cycle
             cnt = cnt + 1
             associate(s => pt(ip)%sheet(cnt))
               s%b = b;  s%n = n;  s%ax = ax;  s%sd = sd;  s%t1 = t1;  s%t2 = t2;  s%rr = rr
               s%n1 = rr(2,t1) - rr(1,t1) + 1;  s%n2 = rr(2,t2) - rr(1,t2) + 1
             end associate
          end do
       end do
       pt(ip)%nsheet = cnt
       call build_structure(pt(ip))
       call setup_morph(pt(ip))
       if (myrank == 0) call write_shell_mesh(pt(ip))
       if (myrank == 0) call open_monitors(pt(ip))
    end do
    fsi_morph_fn => fsi_morph
    fsi_hvel_fn  => fsi_hvel
    fsi_flux_fn  => fsi_flux_reproject
    inited = .true.
    if (myrank == 0) then
       write(*,'(A)') ' '
       write(*,'(A)') ' ---------------- fluid-structure interaction (update 24) ----------------'
       do ip = 1, npart
          if (frames_part(pt(ip)%k) > 0) then
             call frames_info(pt(ip)%k)
          else
             call print_part(pt(ip))
          end if
       end do
       if (any([(frames_part(pt(ip)%k) == 0, ip = 1, npart)])) then
          write(*,'(A)') '  coupling: strong (sub-iterations with Aitken relaxation inside every time step), the structure'
          write(*,'(A)') '            in the part frame, the fluid grid of the part morphed and the overset re-assembled'
       end if
    end if
    is = 0;  ierr = 0
  end subroutine fsi_init

  subroutine print_part(p)
    type(part_t), intent(in), target :: p
    integer :: i
    write(*,'(A)') '  '//sixdof_part_name(p%k)//': flexible (SW_FSI = 1)'
    write(*,'(A,ES11.4,A,F6.3,A,ES11.4,A)') '    E = ', p%def%E, ' Pa, nu = ', p%def%nu, ', rho = ', p%def%rho, ' kg/m3'
    write(*,'(A,ES11.4,A,ES11.4,A,ES11.4,A)') '    enclosed volume ', p%volume, ' m3, wall area ', p%area, ' m2, mass ', p%mass, ' kg'
    write(*,'(A,I0,A,I0,A,I0,A)') '    mid-surface: ', p%sh%nn, ' nodes, ', p%sh%ne, ' shell elements from ', p%nsheet, ' wall patch(es)'
    do i = 1, p%nsheet
       if (p%sheet(i)%passive) then
          write(*,'(A,I0,A,I0,A,I0,A,I0,A)') '      block ', p%sheet(i)%b, ' patch ', p%sheet(i)%n, ': ', p%sheet(i)%n1, ' x ', &
             p%sheet(i)%n2, ' wall nodes, PASSIVE (overlapping or not sheet-like: attached to the mid-surface of the other sheets)'
       else
          write(*,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,L1)') '      block ', p%sheet(i)%b, ' patch ', p%sheet(i)%n, ': ', p%sheet(i)%n1, ' x ', &
             p%sheet(i)%n2, ' wall nodes, ', p%sheet(i)%nup, ' upper nodes x ', p%sheet(i)%nsp, ' span stations, periodic ', p%sheet(i)%periodic
       end if
    end do
    write(*,'(A,ES10.3,A,ES10.3,A,ES10.3,A)') '    thickness from the surface grid: ', p%sh%tmin, ' .. ', p%sh%tmax, &
       ' m (floor ', p%def%tmin_frac, ' x max)'
    write(*,'(A,I0,A,I0)') '    clamped nodes (NoFSI groups): ', p%nclamp, ' of ', p%sh%nn
    if (p%nclamp == 0) write(*,'(A)') '    WARNING: no clamped node (Num_NoFSI_Group = 0) - a free flexible part'
    if (p%def%nsub > 0) then
       do i = 1, p%def%nsub
          write(*,'(A,I0,A,I0,A,ES10.3,A,ES10.3,A,F6.3)') '    Submat region ', i, ': ', count(p%matid == i), ' elements, rho ', &
             p%def%sub_rho(i), ' E ', p%def%sub_E(i), ' nu ', p%def%sub_nu(i)
       end do
    end if
    write(*,'(A,ES11.4,A,F9.4,A)') '    shell mass ', p%sh%mtot, ' kg, first natural frequency (in vacuum) ', p%sh%omega1/(2.0_dp*pi), ' Hz'
    if (p%def%zeta > 0.0_dp) write(*,'(A,F6.3,A,ES10.3,A)') '    Rayleigh damping ratio ', p%def%zeta, ' at the first mode (a_K = ', p%sh%aK, ' s)'
    write(*,'(A,ES10.3,A,I0,A,ES9.2,A,F5.2)') '    FSI_Start_T ', p%def%start_t, ' s, FSI_Max_Iter ', p%def%max_iter, ', FSI_Tol ', &
       p%def%tol, ', FSI_Relax ', p%def%relax
    if (p%def%mode == 1) then
       write(*,'(A)') '    dynamic (generalized-alpha, rho_inf = 0.8, Newton iterations)'
    else
       write(*,'(A)') '    quasi-static equilibrium at every sub-iteration (FSI_Mode = 0)'
    end if
    do i = 1, p%def%nmon
       if (p%mon_sheet(i) > 0) write(*,'(A,I0,A,4I5,A)') '    monitor node ', i, ': block i j k =', p%def%mon(:,i), &
          ' -> '//out_dir//trim(casename)//'_fsi_'//sixdof_part_stem(p%k)//'_node'//itoa(i)//'.dat'
    end do
  end subroutine print_part

  function itoa(i) result(s)
    integer, intent(in) :: i
    character(len=:), allocatable :: s
    character(len=16) :: a
    write(a,'(I0)') i
    s = trim(a)
  end function itoa

  !=====================================================================
  ! The structural mesh of a part from its sheets (replicated on all ranks)
  !=====================================================================
  subroutine build_structure(p)
    type(part_t), intent(inout), target :: p
    integer :: is, i1, i2, m, sidx, nn, ne, e, a, c, nmax, i, j, nbuf, idx(3), jdx(3), q, best, cntq
    real(dp), allocatable :: bufx(:,:), buft(:)
    integer, allocatable :: bufsrc(:,:), gid(:), srcs(:,:)
    real(dp) :: tol, bmin(3), bmax(3), vol3(3), area, t, tmed, x(3), xh(3), nh(3), d(3), w4(4), r(3), thk
    logical :: hit
    type(sheet_t), pointer :: s
    ! ---- geometry of every sheet: node positions, normals, sides, upper range ----
    bmin = huge(1.0_dp);  bmax = -huge(1.0_dp)
    vol3 = 0.0_dp;  area = 0.0_dp
    do is = 1, p%nsheet
       s => p%sheet(is)
       call sheet_geometry(s, vol3, area)
       do i2 = 1, s%n2
          do i1 = 1, s%n1
             x = node_x(s%b, node_idx(s, i1, i2))
             bmin = min(bmin, x);  bmax = max(bmax, x)
          end do
       end do
    end do
    ! overlapping wall patches (overset blocks of one blade: a tip cap over the root block): a
    ! sheet most of whose wall nodes lie on another, larger sheet is passive too - its nodes
    ! follow the other sheet's mid-surface instead of forming a second, unconnected structure
    call mark_overlapped_sheets(p, vol3, area)
    if (count(.not. p%sheet(:)%passive) == 0) call fsi_abort(sixdof_part_name(p%k)//': none of its wall patches is '// &
       'sheet-like (closed or folded around direction) - check the grid, or give the blocks to a rigid part')
    ! enclosed volume: the three divergence-theorem estimates int x_c n_c dA agree for a closed
    ! surface; where the surface is open (the span ends of a quasi-2-D plate, a blade root) the
    ! estimate along the opening's normal is wrong, so the median is taken
    vol3 = abs(vol3)
    p%volume = max(min(vol3(1), vol3(2)), min(max(vol3(1), vol3(2)), vol3(3)))
    p%area = area
    tol = 1.0e-6_dp*maxval(bmax - bmin);  geo_tol = tol
    ! ---- mid nodes of every sheet (ray pairing), collected in a buffer ----
    nmax = 0
    do is = 1, p%nsheet
       nmax = nmax + p%sheet(is)%nup*p%sheet(is)%nsp
    end do
    allocate(bufx(3, nmax), buft(nmax), bufsrc(4, nmax), gid(nmax))
    nbuf = 0
    do is = 1, p%nsheet
       s => p%sheet(is)
       if (s%passive) cycle
       allocate(s%snode(s%nup, s%nsp))
       do sidx = 1, s%nsp
          do m = 1, s%nup
             call pair_node(s, s%iup(m), span_index(s, sidx), xh, nh, thk, hit, idx)
             x = node_x(s%b, node_idx(s, ar1(s, s%iup(m), sidx), ar2(s, s%iup(m), sidx)))
             nbuf = nbuf + 1
             bufx(:, nbuf) = 0.5_dp*(x + xh);  buft(nbuf) = thk
             bufsrc(1:2, nbuf) = [ar1(s, s%iup(m), sidx), ar2(s, s%iup(m), sidx)]   ! upper node (patch indices)
             bufsrc(3:4, nbuf) = idx(1:2)                                           ! lower node (patch indices, 0 = none)
             s%snode(m, sidx) = nbuf
          end do
       end do
    end do
    ! ---- merge coincident nodes of different sheets ----
    nn = 0
    do i = 1, nbuf
       gid(i) = 0
       do j = 1, i - 1
          if (gid(j) == 0) cycle
          if (sum((bufx(:,i) - bufx(:,j))**2) < tol*tol) then
             gid(i) = gid(j);  exit
          end if
       end do
       if (gid(i) == 0) then
          nn = nn + 1;  gid(i) = nn
       end if
    end do
    ! thickness floor and the structural node list
    tmed = 0.0_dp
    allocate(p%sh%X0(3, nn), p%thick(nn), srcs(5, nn))
    p%thick = 0.0_dp;  srcs = 0
    do i = 1, nbuf
       p%sh%X0(:, gid(i)) = bufx(:, i)
       p%thick(gid(i)) = max(p%thick(gid(i)), buft(i))
       if (srcs(1, gid(i)) == 0) then
          srcs(1, gid(i)) = sheet_of_buf(p, i);  srcs(2:5, gid(i)) = bufsrc(:, i)
       end if
    end do
    do is = 1, p%nsheet
       s => p%sheet(is)
       if (s%passive) cycle
       do sidx = 1, s%nsp
          do m = 1, s%nup
             s%snode(m, sidx) = gid(s%snode(m, sidx))
          end do
       end do
    end do
    ! floor of the thickness: FSI_Min_Thickness x the largest thickness of the nodes outside the
    ! NoFSI groups (a blade root buried in a hub is thick but clamped)
    t = 0.0_dp;  tmed = 0.0_dp
    do i = 1, nn
       hit = .false.
       do q = 1, p%def%nnof
          if (in_region(p, srcs(:, i), p%def%nof_blk(q), p%def%nof_rng(:,:,q))) hit = .true.
       end do
       if (.not. hit) t = max(t, p%thick(i))
       tmed = max(tmed, p%thick(i))
    end do
    if (t <= 0.0_dp) t = tmed
    p%thick = max(p%thick, p%def%tmin_frac*t)
    p%sh%nn = nn
    ! ---- elements ----
    ne = 0
    do is = 1, p%nsheet
       if (p%sheet(is)%passive) cycle
       ne = ne + (p%sheet(is)%nup - 1)*(p%sheet(is)%nsp - 1)
    end do
    p%sh%ne = ne
    allocate(p%sh%conn(4, ne), p%sh%te(ne), p%sh%Ee(ne), p%sh%nue(ne), p%sh%rhoe(ne), p%sh%fixed(6, nn), p%matid(ne))
    e = 0
    do is = 1, p%nsheet
       s => p%sheet(is)
       if (s%passive) cycle
       do sidx = 1, s%nsp - 1
          do m = 1, s%nup - 1
             e = e + 1
             p%sh%conn(:, e) = [s%snode(m, sidx), s%snode(m+1, sidx), s%snode(m+1, sidx+1), s%snode(m, sidx+1)]
          end do
       end do
    end do
    ! ---- density / mass ----
    if (p%def%rho > 0.0_dp) then
       p%mass = p%def%rho*p%volume
       call sixdof_mass_data(p%k, p%mass, x, xh, .true.)
    else
       p%mass = 0.0_dp
       call sixdof_mass_data(p%k, p%mass, x, xh, .false.)
       if (p%mass <= 0.0_dp) call fsi_abort(sixdof_part_name(p%k)//': SW_FSI = 1 needs Density or Mass > 0')
       if (p%volume <= 0.0_dp) call fsi_abort(sixdof_part_name(p%k)//': the enclosed volume of the wall surface is 0 - give Density')
       p%def%rho = p%mass/p%volume
    end if
    ! ---- element thickness and material ----
    do e = 1, ne
       p%sh%te(e) = 0.25_dp*sum(p%thick(p%sh%conn(:, e)))
       p%sh%Ee(e) = p%def%E;  p%sh%nue(e) = p%def%nu;  p%sh%rhoe(e) = p%def%rho
       p%matid(e) = 0
       best = 0;  cntq = 0
       do q = 1, p%def%nsub
          c = 0
          do a = 1, 4
             if (in_region(p, srcs(:, p%sh%conn(a, e)), p%def%sub_blk(q), p%def%sub_rng(:,:,q))) c = c + 1
          end do
          if (c > cntq) then
             cntq = c;  best = q
          end if
       end do
       if (best > 0 .and. cntq >= 2) then
          p%matid(e) = best
          if (p%def%sub_E(best) > 0.0_dp) p%sh%Ee(e) = p%def%sub_E(best)
          if (p%def%sub_nu(best) >= 0.0_dp) p%sh%nue(e) = p%def%sub_nu(best)
          if (p%def%sub_rho(best) > 0.0_dp) p%sh%rhoe(e) = p%def%sub_rho(best)
       end if
    end do
    ! ---- clamped nodes ----
    p%sh%fixed = .false.;  p%nclamp = 0
    do i = 1, nn
       do q = 1, p%def%nnof
          if (in_region(p, srcs(:, i), p%def%nof_blk(q), p%def%nof_rng(:,:,q))) p%sh%fixed(:, i) = .true.
       end do
       if (all(p%sh%fixed(:, i))) p%nclamp = p%nclamp + 1
       ! FSI_Plane_Strain (24c): a quasi-2-D case bends cylindrically - no displacement along the span
       ! axis, no rotation about the two other axes (without it the plate twisted along the span and
       ! bent the span walls of the grid)
       if (p%def%plane > 0) then
          p%sh%fixed(p%def%plane, i) = .true.
          do q = 1, 3
             if (q /= p%def%plane) p%sh%fixed(3 + q, i) = .true.
          end do
       end if
    end do
    ! ---- shell set-up ----
    p%sh%zeta = p%def%zeta
    p%sh%rho_inf = 0.8_dp
    p%sh%max_newton = 25
    call shell_setup(p%sh)
    if (p%def%plane > 0) call plane_groups(p)
    allocate(p%fext(6, nn), p%drel(6, nn), p%rprev(6, nn), p%vprev(6, nn), p%dasm(6, nn), p%vm(nn), &
             p%rhist(6*nn, MAXQN), p%dhist(6*nn, MAXQN))
    p%fext = 0.0_dp;  p%drel = 0.0_dp;  p%rprev = 0.0_dp;  p%vprev = 0.0_dp;  p%dasm = 0.0_dp;  p%vm = 0.0_dp
    allocate(p%vd(6, nn));  p%vd = 0.0_dp
    ! ---- attachment of every wall node ----
    do is = 1, p%nsheet
       s => p%sheet(is)
       allocate(s%att(4, s%n1, s%n2), s%wgt(4, s%n1, s%n2), s%off(3, 4, s%n1, s%n2), s%dw(3, s%n1, s%n2), &
                s%thw(3, s%n1, s%n2), s%fw(3, s%n1, s%n2))
       s%att = 0;  s%wgt = 0.0_dp;  s%off = 0.0_dp;  s%dw = 0.0_dp;  s%thw = 0.0_dp;  s%fw = 0.0_dp
       call attach_nodes(p, s)
    end do
    ! ---- monitor nodes ----
    do q = 1, p%def%nmon
       p%mon_sheet(q) = 0
       do is = 1, p%nsheet
          s => p%sheet(is)
          if (s%b /= p%def%mon(1, q)) cycle
          jdx = p%def%mon(2:4, q)
          if (jdx(s%ax) /= merge(1, ncell(s%b, s%ax) + 1, s%sd == 1)) cycle
          if (jdx(s%t1) < s%rr(1,s%t1) .or. jdx(s%t1) > s%rr(2,s%t1)) cycle
          if (jdx(s%t2) < s%rr(1,s%t2) .or. jdx(s%t2) > s%rr(2,s%t2)) cycle
          p%mon_sheet(q) = is;  p%mon_i1(q) = jdx(s%t1) - s%rr(1,s%t1) + 1;  p%mon_i2(q) = jdx(s%t2) - s%rr(1,s%t2) + 1
       end do
       if (p%mon_sheet(q) == 0) call fsi_abort(sixdof_part_name(p%k)//': FSI_Monitor_node'//itoa(q)// &
          ' is not a wall node of the part (block i j k must lie on a no-slip patch of one of its blocks)')
    end do
    deallocate(bufx, buft, bufsrc, gid, srcs)
    d = 0.0_dp;  w4 = 0.0_dp;  r = 0.0_dp;  tmed = 0.0_dp
  end subroutine build_structure

  ! sheet number of buffer entry i (the buffer is filled sheet by sheet)
  integer function sheet_of_buf(p, i)
    type(part_t), intent(in) :: p
    integer, intent(in) :: i
    integer :: is, c
    c = 0
    do is = 1, p%nsheet
       if (p%sheet(is)%passive) cycle
       c = c + p%sheet(is)%nup*p%sheet(is)%nsp
       if (i <= c) then
          sheet_of_buf = is;  return
       end if
    end do
    sheet_of_buf = p%nsheet
  end function sheet_of_buf

  ! sheets overlapped by a larger sheet of the same part become passive (the volume and area
  ! of the part are recomputed without them)
  subroutine mark_overlapped_sheets(p, vol3, area)
    type(part_t), intent(inout), target :: p
    real(dp), intent(inout) :: vol3(3), area
    integer :: is, js, i1, i2, j1, j2, nover, ntot, order(64), nact, q, tmp
    real(dp) :: x(3), dmin, dd, h, xw1(3), xw2(3)
    logical :: changed
    type(sheet_t), pointer :: s, o
    nact = 0
    do is = 1, p%nsheet
       if (p%sheet(is)%passive) cycle
       nact = nact + 1;  order(nact) = is
    end do
    if (nact < 2) return
    ! largest sheet first
    do is = 1, nact - 1
       do js = is + 1, nact
          if (p%sheet(order(js))%n1*p%sheet(order(js))%n2 > p%sheet(order(is))%n1*p%sheet(order(is))%n2) then
             tmp = order(is);  order(is) = order(js);  order(js) = tmp
          end if
       end do
    end do
    changed = .false.
    do q = 2, nact
       s => p%sheet(order(q))
       nover = 0;  ntot = s%n1*s%n2
       do i2 = 1, s%n2
          do i1 = 1, s%n1
             x = node_x(s%b, node_idx(s, i1, i2))
             ! nearest wall node of the larger, still active sheets ahead in the order
             dmin = huge(1.0_dp);  h = 0.0_dp
             do js = 1, q - 1
                o => p%sheet(order(js))
                if (o%passive) cycle
                do j2 = 1, o%n2
                   do j1 = 1, o%n1
                      xw1 = node_x(o%b, node_idx(o, j1, j2))
                      dd = sum((xw1 - x)**2)
                      if (dd < dmin) then
                         dmin = dd
                         ! local spacing of the larger sheet: the longer of its two node distances
                         xw2 = node_x(o%b, node_idx(o, min(j1+1, o%n1), j2))
                         if (j1 == o%n1) xw2 = node_x(o%b, node_idx(o, j1-1, j2))
                         h = sqrt(sum((xw2 - xw1)**2))
                         xw2 = node_x(o%b, node_idx(o, j1, min(j2+1, o%n2)))
                         if (j2 == o%n2) xw2 = node_x(o%b, node_idx(o, j1, j2-1))
                         h = max(h, sqrt(sum((xw2 - xw1)**2)))
                      end if
                   end do
                end do
             end do
             if (sqrt(dmin) < h) nover = nover + 1
          end do
       end do
       ! more than 40 % of the nodes lie on the larger sheet: an overlapping patch
       if (5*nover >= 2*ntot) then
          s%passive = .true.;  s%nup = 0;  s%nsp = 0;  s%side = 0;  changed = .true.
          if (allocated(s%iup)) deallocate(s%iup)
       end if
    end do
    if (changed) then
       vol3 = 0.0_dp;  area = 0.0_dp
       do is = 1, p%nsheet
          s => p%sheet(is)
          if (s%passive) cycle
          call sheet_area_volume(s, vol3, area)
       end do
    end if
  end subroutine mark_overlapped_sheets

  ! volume / area contribution of a sheet (as in sheet_geometry)
  subroutine sheet_area_volume(s, vol3, area)
    type(sheet_t), intent(in) :: s
    real(dp), intent(inout) :: vol3(3), area
    real(dp) :: x1(3), x2(3), x3(3), x4(3), Av(3), xf(3), xin(3)
    integer :: i1, i2, idx(3)
    do i2 = 1, s%n2 - 1
       do i1 = 1, s%n1 - 1
          x1 = node_x(s%b, node_idx(s, i1, i2));  x2 = node_x(s%b, node_idx(s, i1+1, i2))
          x3 = node_x(s%b, node_idx(s, i1+1, i2+1));  x4 = node_x(s%b, node_idx(s, i1, i2+1))
          Av = 0.5_dp*cross3(x3 - x1, x4 - x2)
          xf = 0.25_dp*(x1 + x2 + x3 + x4)
          idx = node_idx(s, i1, i2);  idx(s%ax) = idx(s%ax) + merge(1, -1, s%sd == 1)
          xin = node_x(s%b, idx)
          if (sum(Av*(xin - xf)) < 0.0_dp) Av = -Av
          vol3 = vol3 + xf*Av
          area = area + sqrt(sum(Av*Av))
       end do
    end do
  end subroutine sheet_area_volume

  ! is the (upper or lower) source node of a structural node inside the block index range?
  logical function in_region(p, src, b, rng)
    type(part_t), intent(in) :: p
    integer, intent(in) :: src(5), b, rng(2,3)
    integer :: is, idx(3), q, lo(3), hi(3)
    in_region = .false.
    is = src(1)
    if (is == 0) return
    if (p%sheet(is)%b /= b) return
    do q = 1, 3
       lo(q) = rng(1,q);  hi(q) = rng(2,q)
       if (lo(q) == 0) then
          lo(q) = 1;  hi(q) = huge(1)
       end if
       if (hi(q) == 0) hi(q) = lo(q)
    end do
    idx = node_idx(p%sheet(is), src(2), src(3))
    if (all(idx >= lo) .and. all(idx <= hi)) in_region = .true.
    if (src(4) > 0) then
       idx = node_idx(p%sheet(is), src(4), src(5))
       if (all(idx >= lo) .and. all(idx <= hi)) in_region = .true.
    end if
  end function in_region

  ! patch indices (i1, i2) of around index ia and span station sidx
  pure integer function ar1(s, ia, sidx)
    type(sheet_t), intent(in) :: s
    integer, intent(in) :: ia, sidx
    ar1 = merge(ia, sidx, s%ar == 1)
  end function ar1
  pure integer function ar2(s, ia, sidx)
    type(sheet_t), intent(in) :: s
    integer, intent(in) :: ia, sidx
    ar2 = merge(sidx, ia, s%ar == 1)
  end function ar2
  pure integer function span_index(s, sidx)
    type(sheet_t), intent(in) :: s
    integer, intent(in) :: sidx
    span_index = sidx
  end function span_index

  !---------------------------------------------------------------------
  ! Geometry of a sheet: normals, periodicity, around direction, sides,
  ! upper range; adds the patch's volume / area contributions
  !---------------------------------------------------------------------
  subroutine sheet_geometry(s, vol3, area)
    type(sheet_t), intent(inout) :: s
    real(dp), intent(inout) :: vol3(3), area
    real(dp), allocatable :: xw(:,:,:), nw(:,:,:)
    real(dp) :: x1(3), x2(3), x3(3), x4(3), Av(3), xf(3), xin(3), c(3), chord(3), esp(3), nup(3), nup0(3)
    real(dp) :: cov(3,3), ev(3), w, dn
    integer :: i1, i2, na, ns, ia, sidx, idx(3), nupc, best, bestlen, start, len_, j, jj, mid, q
    logical, allocatable :: up(:)
    integer, allocatable :: run(:)
    allocate(xw(3, s%n1, s%n2), nw(3, s%n1, s%n2))
    nw = 0.0_dp
    do i2 = 1, s%n2
       do i1 = 1, s%n1
          xw(:, i1, i2) = node_x(s%b, node_idx(s, i1, i2))
       end do
    end do
    ! face normals (into the fluid) accumulated at the nodes; volume by the divergence theorem (outward = into the fluid)
    do i2 = 1, s%n2 - 1
       do i1 = 1, s%n1 - 1
          x1 = xw(:, i1, i2);  x2 = xw(:, i1+1, i2);  x3 = xw(:, i1+1, i2+1);  x4 = xw(:, i1, i2+1)
          Av = 0.5_dp*cross3(x3 - x1, x4 - x2)
          xf = 0.25_dp*(x1 + x2 + x3 + x4)
          idx = node_idx(s, i1, i2);  idx(s%ax) = idx(s%ax) + merge(1, -1, s%sd == 1)
          xin = node_x(s%b, idx)                                    ! first node inside the fluid
          if (sum(Av*(xin - xf)) < 0.0_dp) Av = -Av
          nw(:, i1, i2) = nw(:, i1, i2) + Av;  nw(:, i1+1, i2) = nw(:, i1+1, i2) + Av
          nw(:, i1+1, i2+1) = nw(:, i1+1, i2+1) + Av;  nw(:, i1, i2+1) = nw(:, i1, i2+1) + Av
          vol3 = vol3 + xf*Av                                      ! the three estimates int x_c n_c dA
          area = area + sqrt(sum(Av*Av))
       end do
    end do
    do i2 = 1, s%n2
       do i1 = 1, s%n1
          nw(:, i1, i2) = unit3(nw(:, i1, i2))
       end do
    end do
    ! periodicity: first and last node line coincide
    dn = 1.0e-8_dp*maxval(abs(xw))
    s%periodic = .false.;  s%ar = 0
    if (maxval(abs(xw(:, 1, :) - xw(:, s%n1, :))) < max(dn, 1.0e-12_dp)) then
       s%ar = 1;  s%periodic = .true.
    else if (maxval(abs(xw(:, :, 1) - xw(:, :, s%n2))) < max(dn, 1.0e-12_dp)) then
       s%ar = 2;  s%periodic = .true.
    else
       mid = (s%n2 + 1)/2
       if (sum(nw(:, 1, mid)*nw(:, s%n1, mid)) < 0.0_dp) then
          s%ar = 1
       else
          mid = (s%n1 + 1)/2
          if (sum(nw(:, mid, 1)*nw(:, mid, s%n2)) < 0.0_dp) s%ar = 2
       end if
    end if
    if (s%ar == 0) then
       ! not sheet-like (no closed or folded direction): a PASSIVE patch - its wall nodes are
       ! attached to the mid-surface of the other sheets of the part (a tip cap, a fillet block)
       s%passive = .true.;  s%nup = 0;  s%nsp = 0
       allocate(s%side(s%n1, s%n2));  s%side = 0
       deallocate(xw, nw)
       return
    end if
    na = merge(s%n1, s%n2, s%ar == 1);  ns = merge(s%n2, s%n1, s%ar == 1)
    s%nsp = ns
    ! sides per span station
    allocate(s%side(s%n1, s%n2), up(na))
    s%side = 0
    nup0 = 0.0_dp
    do sidx = 1, ns
       ! chord = principal axis of the section, e_span = direction of the neighbouring station
       c = 0.0_dp
       do ia = 1, na
          c = c + xw(:, ar1(s, ia, sidx), ar2(s, ia, sidx))
       end do
       c = c/na
       cov = 0.0_dp
       do ia = 1, na
          xf = xw(:, ar1(s, ia, sidx), ar2(s, ia, sidx)) - c
          do q = 1, 3
             cov(:, q) = cov(:, q) + xf*xf(q)
          end do
       end do
       call principal_axis(cov, chord)
       esp = 0.0_dp
       do ia = 1, na
          if (sidx < ns) esp = esp + xw(:, ar1(s, ia, sidx+1), ar2(s, ia, sidx+1)) - xw(:, ar1(s, ia, sidx), ar2(s, ia, sidx))
          if (sidx > 1)  esp = esp + xw(:, ar1(s, ia, sidx), ar2(s, ia, sidx)) - xw(:, ar1(s, ia, sidx-1), ar2(s, ia, sidx-1))
       end do
       esp = unit3(esp)
       nup = unit3(cross3(chord, esp))
       if (sidx > 1) then
          if (sum(nup*nup0) < 0.0_dp) nup = -nup
       end if
       nup0 = nup
       do ia = 1, na
          w = sum(nw(:, ar1(s, ia, sidx), ar2(s, ia, sidx))*nup)
          if (w > 0.5_dp) then
             s%side(ar1(s, ia, sidx), ar2(s, ia, sidx)) = 1
          else if (w < -0.5_dp) then
             s%side(ar1(s, ia, sidx), ar2(s, ia, sidx)) = -1
          end if
       end do
    end do
    ! upper range: the longest (cyclic) run of upper nodes at the middle station
    mid = (ns + 1)/2
    do ia = 1, na
       up(ia) = (s%side(ar1(s, ia, mid), ar2(s, ia, mid)) == 1)
    end do
    nupc = merge(na - 1, na, s%periodic)          ! distinct nodes
    bestlen = 0;  best = 0
    do start = 1, nupc
       if (.not. up(start)) cycle
       len_ = 0
       do j = 0, nupc - 1
          jj = start + j
          if (jj > nupc) then
             if (.not. s%periodic) exit
             jj = jj - nupc
          end if
          if (.not. up(jj)) exit
          len_ = len_ + 1
       end do
       if (len_ > bestlen) then
          bestlen = len_;  best = start
       end if
    end do
    if (bestlen < 2) then
       s%passive = .true.;  s%nup = 0;  s%nsp = 0;  s%side = 0
       deallocate(xw, nw, up)
       return
    end if
    s%nup = bestlen
    allocate(s%iup(s%nup))
    do j = 1, s%nup
       jj = best + j - 1
       if (jj > nupc) jj = jj - nupc
       s%iup(j) = jj
    end do
    ! the upper range is used at every station: mark it
    do sidx = 1, ns
       do j = 1, s%nup
          s%side(ar1(s, s%iup(j), sidx), ar2(s, s%iup(j), sidx)) = 1
       end do
    end do
    deallocate(xw, nw, up)
    allocate(run(1));  deallocate(run)
    ev = 0.0_dp
  end subroutine sheet_geometry

  ! eigenvector of the largest eigenvalue of a symmetric 3 x 3 matrix (power iteration)
  subroutine principal_axis(cov, v)
    real(dp), intent(in)  :: cov(3,3)
    real(dp), intent(out) :: v(3)
    real(dp) :: w(3)
    integer :: it
    v = [1.0_dp, 0.7_dp, 0.3_dp]
    do it = 1, 60
       w = matmul(cov, v)
       if (sum(w*w) <= 0.0_dp) exit
       v = unit3(w)
    end do
  end subroutine principal_axis

  !---------------------------------------------------------------------
  ! Ray from the wall node (patch indices ia along the around direction,
  ! sidx station) along the inward normal to the opposite side: hit point
  ! xh, hit normal nh, thickness; idx(1:2) = patch indices of the nearest
  ! lower node of the hit face (0 = none).  A miss falls back to the
  ! lower node nearest to the ray.
  !---------------------------------------------------------------------
  subroutine pair_node(s, ia, sidx, xh, nh, thk, hit, idx)
    type(sheet_t), intent(in) :: s
    integer, intent(in) :: ia, sidx
    real(dp), intent(out) :: xh(3), nh(3), thk
    logical, intent(out) :: hit
    integer, intent(out) :: idx(3)
    real(dp) :: xo(3), nrm(3), dir(3), x(3,4), tmin, tt, u, v, best(3), dmin, dd, xl(3), xq(3)
    integer :: i1, i2, j1, j2, c, q, fi1, fi2, a
    i1 = ar1(s, ia, sidx);  i2 = ar2(s, ia, sidx)
    xo = node_x(s%b, node_idx(s, i1, i2))
    nrm = wall_normal(s, i1, i2)
    dir = -nrm
    hit = .false.;  tmin = huge(1.0_dp);  idx = 0
    do j2 = 1, s%n2 - 1
       do j1 = 1, s%n1 - 1
          ! faces with a lower / cap node only, and not the ones touching the origin node
          if (all(s%side(j1:j1+1, j2:j2+1) == 1)) cycle
          x(:,1) = node_x(s%b, node_idx(s, j1, j2));    x(:,2) = node_x(s%b, node_idx(s, j1+1, j2))
          x(:,3) = node_x(s%b, node_idx(s, j1+1, j2+1)); x(:,4) = node_x(s%b, node_idx(s, j1, j2+1))
          if (any([sum((x(:,1)-xo)**2), sum((x(:,2)-xo)**2), sum((x(:,3)-xo)**2), sum((x(:,4)-xo)**2)] < geo_tol**2)) cycle
          do c = 1, 2
             if (c == 1) then
                call ray_tri(xo, dir, x(:,1), x(:,2), x(:,3), tt, u, v, q)
             else
                call ray_tri(xo, dir, x(:,1), x(:,3), x(:,4), tt, u, v, q)
             end if
             if (q == 1 .and. tt > 10.0_dp*geo_tol .and. tt < tmin) then
                tmin = tt;  hit = .true.;  fi1 = j1;  fi2 = j2
             end if
          end do
       end do
    end do
    if (hit) then
       xh = xo + tmin*dir;  thk = tmin
       ! normal of the hit face (into the fluid) and its nearest non-upper node
       x(:,1) = node_x(s%b, node_idx(s, fi1, fi2));    x(:,2) = node_x(s%b, node_idx(s, fi1+1, fi2))
       x(:,3) = node_x(s%b, node_idx(s, fi1+1, fi2+1)); x(:,4) = node_x(s%b, node_idx(s, fi1, fi2+1))
       nh = unit3(cross3(x(:,3) - x(:,1), x(:,4) - x(:,2)))
       if (sum(nh*wall_normal(s, fi1, fi2)) < 0.0_dp) nh = -nh
       dmin = huge(1.0_dp)
       do a = 0, 3
          j1 = fi1 + merge(1, 0, a == 1 .or. a == 2);  j2 = fi2 + merge(1, 0, a >= 2)
          if (s%side(j1, j2) == 1) cycle
          dd = sum((node_x(s%b, node_idx(s, j1, j2)) - xh)**2)
          if (dd < dmin) then
             dmin = dd;  idx(1) = j1;  idx(2) = j2
          end if
       end do
       if (sum(nh*nrm) > -0.2_dp) hit = .false.       ! not an opposed surface: fall back
    end if
    if (.not. hit) then
       ! lower node nearest to the ray line, ahead of the origin
       dmin = huge(1.0_dp);  best = xo
       do j2 = 1, s%n2
          do j1 = 1, s%n1
             if (s%side(j1, j2) /= -1) cycle
             xl = node_x(s%b, node_idx(s, j1, j2))
             tt = sum((xl - xo)*dir)
             if (tt <= 0.0_dp) cycle
             xq = xo + tt*dir
             dd = sum((xl - xq)**2) + 0.01_dp*tt*tt
             if (dd < dmin) then
                dmin = dd;  best = xl;  idx(1) = j1;  idx(2) = j2
             end if
          end do
       end do
       xh = best;  thk = max(sum((best - xo)*dir), 0.0_dp)
       nh = -nrm
    end if
  end subroutine pair_node

  ! unit wall normal (into the fluid) at patch node (i1, i2) from the adjacent faces
  function wall_normal(s, i1, i2) result(nrm)
    type(sheet_t), intent(in) :: s
    integer, intent(in) :: i1, i2
    real(dp) :: nrm(3), x1(3), x2(3), x3(3), x4(3), Av(3), xf(3), xin(3)
    integer :: a, c, idx(3)
    nrm = 0.0_dp
    do c = max(i2-1, 1), min(i2, s%n2-1)
       do a = max(i1-1, 1), min(i1, s%n1-1)
          x1 = node_x(s%b, node_idx(s, a, c));  x2 = node_x(s%b, node_idx(s, a+1, c))
          x3 = node_x(s%b, node_idx(s, a+1, c+1));  x4 = node_x(s%b, node_idx(s, a, c+1))
          Av = 0.5_dp*cross3(x3 - x1, x4 - x2)
          xf = 0.25_dp*(x1 + x2 + x3 + x4)
          idx = node_idx(s, a, c);  idx(s%ax) = idx(s%ax) + merge(1, -1, s%sd == 1)
          xin = node_x(s%b, idx)
          if (sum(Av*(xin - xf)) < 0.0_dp) Av = -Av
          nrm = nrm + Av
       end do
    end do
    nrm = unit3(nrm)
  end function wall_normal

  ! Moeller-Trumbore ray / triangle intersection: q = 1 hit at t with barycentric u, v
  pure subroutine ray_tri(o, d, v0, v1, v2, t, u, v, q)
    real(dp), intent(in)  :: o(3), d(3), v0(3), v1(3), v2(3)
    real(dp), intent(out) :: t, u, v
    integer, intent(out)  :: q
    real(dp) :: e1(3), e2(3), p(3), det, inv, tv(3), qv(3)
    q = 0;  t = 0.0_dp;  u = 0.0_dp;  v = 0.0_dp
    e1 = v1 - v0;  e2 = v2 - v0
    p = cross3(d, e2)
    det = sum(e1*p)
    if (abs(det) < 1.0e-30_dp) return
    inv = 1.0_dp/det
    tv = o - v0
    u = sum(tv*p)*inv
    if (u < -1.0e-6_dp .or. u > 1.0_dp + 1.0e-6_dp) return
    qv = cross3(tv, e1)
    v = sum(d*qv)*inv
    if (v < -1.0e-6_dp .or. u + v > 1.0_dp + 1.0e-6_dp) return
    t = sum(e2*qv)*inv
    q = 1
  end subroutine ray_tri

  !---------------------------------------------------------------------
  ! Attachment of the wall nodes of a sheet to the mid-surface
  !---------------------------------------------------------------------
  subroutine attach_nodes(p, s)
    type(part_t), intent(inout), target :: p
    type(sheet_t), intent(inout) :: s
    integer :: i1, i2, m, sidx, a, nd, nb, idx(3), m0, s0, e1, e2, mm, ss, c, q, js, ks
    real(dp) :: x(3), pm(3), xh(3), nh(3), thk, dmin, dd, xi, eta, w(4), xc(3,4), best(2), dbest
    logical :: hit, found
    type(sheet_t), pointer :: o
    do i2 = 1, s%n2
       do i1 = 1, s%n1
          x = node_x(s%b, node_idx(s, i1, i2))
          ! the point on the mid-surface this wall node belongs to
          if (s%passive) then
             pm = x
          else if (s%side(i1, i2) == 1) then
             ! upper node: its own mid node
             m = 0
             do mm = 1, s%nup
                if (s%iup(mm) == merge(i1, i2, s%ar == 1)) m = mm
             end do
             sidx = merge(i2, i1, s%ar == 1)
             if (m > 0) then
                nd = s%snode(m, sidx)
                s%att(1, i1, i2) = nd;  s%wgt(1, i1, i2) = 1.0_dp;  s%off(:, 1, i1, i2) = x - p%sh%X0(:, nd)
                cycle
             end if
             pm = x
          else if (s%side(i1, i2) == -1) then
             call pair_node_any(s, i1, i2, xh, hit)
             pm = 0.5_dp*(x + xh)
          else
             pm = x
          end if
          ! nearest mid node of this sheet (a passive sheet: of every active sheet of the part)
          dmin = huge(1.0_dp);  m0 = 1;  s0 = 1;  js = 0
          do ks = 1, p%nsheet
             if (p%sheet(ks)%passive) cycle
             if (.not. s%passive .and. ks /= own_sheet(p, s)) cycle
             do sidx = 1, p%sheet(ks)%nsp
                do mm = 1, p%sheet(ks)%nup
                   dd = sum((p%sh%X0(:, p%sheet(ks)%snode(mm, sidx)) - pm)**2)
                   if (dd < dmin) then
                      dmin = dd;  m0 = mm;  s0 = sidx;  js = ks
                   end if
                end do
             end do
          end do
          if (js == 0) call fsi_abort('attach_nodes: no active sheet')
          o => p%sheet(js)
          ! the element around it that contains the projection
          found = .false.;  dbest = huge(1.0_dp)
          do e2 = max(s0-1, 1), min(s0, o%nsp-1)
             do e1 = max(m0-1, 1), min(m0, o%nup-1)
                xc(:,1) = p%sh%X0(:, o%snode(e1, e2));    xc(:,2) = p%sh%X0(:, o%snode(e1+1, e2))
                xc(:,3) = p%sh%X0(:, o%snode(e1+1, e2+1)); xc(:,4) = p%sh%X0(:, o%snode(e1, e2+1))
                call inverse_bilinear(xc, pm, xi, eta, dd)
                if (abs(xi) <= 1.05_dp .and. abs(eta) <= 1.05_dp) then
                   if (dd < dbest) then
                      dbest = dd;  found = .true.;  best = [xi, eta];  mm = e1;  ss = e2
                   end if
                end if
             end do
          end do
          if (found) then
             xi = max(-1.0_dp, min(1.0_dp, best(1)));  eta = max(-1.0_dp, min(1.0_dp, best(2)))
             w = 0.25_dp*[(1-xi)*(1-eta), (1+xi)*(1-eta), (1+xi)*(1+eta), (1-xi)*(1+eta)]
             s%att(:, i1, i2) = [o%snode(mm, ss), o%snode(mm+1, ss), o%snode(mm+1, ss+1), o%snode(mm, ss+1)]
             do a = 1, 4
                s%wgt(a, i1, i2) = w(a)
                s%off(:, a, i1, i2) = x - p%sh%X0(:, s%att(a, i1, i2))
             end do
          else
             nd = o%snode(m0, s0)
             s%att(1, i1, i2) = nd;  s%wgt(1, i1, i2) = 1.0_dp;  s%off(:, 1, i1, i2) = x - p%sh%X0(:, nd)
          end if
       end do
    end do
    idx = 0;  nb = 0;  c = 0;  q = 0;  nh = 0.0_dp;  thk = 0.0_dp
  end subroutine attach_nodes

  ! index of sheet s in part p
  integer function own_sheet(p, s)
    type(part_t), intent(in), target :: p
    type(sheet_t), intent(in), target :: s
    integer :: is
    own_sheet = 0
    do is = 1, p%nsheet
       if (p%sheet(is)%b == s%b .and. p%sheet(is)%n == s%n) own_sheet = is
    end do
  end function own_sheet

  ! ray pairing from any wall node (patch indices) towards the opposite side
  subroutine pair_node_any(s, i1, i2, xh, hit)
    type(sheet_t), intent(in) :: s
    integer, intent(in) :: i1, i2
    real(dp), intent(out) :: xh(3)
    logical, intent(out) :: hit
    real(dp) :: xo(3), dir(3), x(3,4), tmin, tt, u, v
    integer :: j1, j2, c, q
    xo = node_x(s%b, node_idx(s, i1, i2))
    dir = -wall_normal(s, i1, i2)
    hit = .false.;  tmin = huge(1.0_dp);  xh = xo
    do j2 = 1, s%n2 - 1
       do j1 = 1, s%n1 - 1
          if (all(s%side(j1:j1+1, j2:j2+1) == -1)) cycle          ! own side
          x(:,1) = node_x(s%b, node_idx(s, j1, j2));    x(:,2) = node_x(s%b, node_idx(s, j1+1, j2))
          x(:,3) = node_x(s%b, node_idx(s, j1+1, j2+1)); x(:,4) = node_x(s%b, node_idx(s, j1, j2+1))
          if (any([sum((x(:,1)-xo)**2), sum((x(:,2)-xo)**2), sum((x(:,3)-xo)**2), sum((x(:,4)-xo)**2)] < geo_tol**2)) cycle
          do c = 1, 2
             if (c == 1) then
                call ray_tri(xo, dir, x(:,1), x(:,2), x(:,3), tt, u, v, q)
             else
                call ray_tri(xo, dir, x(:,1), x(:,3), x(:,4), tt, u, v, q)
             end if
             if (q == 1 .and. tt > 10.0_dp*geo_tol .and. tt < tmin) then
                tmin = tt;  hit = .true.
             end if
          end do
       end do
    end do
    if (hit) xh = xo + tmin*dir
  end subroutine pair_node_any

  ! local coordinates (xi, eta) of the projection of point pm onto the bilinear quad xc; dist = distance
  subroutine inverse_bilinear(xc, pm, xi, eta, dist)
    real(dp), intent(in)  :: xc(3,4), pm(3)
    real(dp), intent(out) :: xi, eta, dist
    real(dp) :: N(4), dNx(4), dNe(4), xp(3), r(3), gx(3), ge(3), A(2,2), bb(2), det, dxi, deta
    integer :: it
    xi = 0.0_dp;  eta = 0.0_dp
    do it = 1, 20
       N   = 0.25_dp*[(1-xi)*(1-eta), (1+xi)*(1-eta), (1+xi)*(1+eta), (1-xi)*(1+eta)]
       dNx = 0.25_dp*[-(1-eta), (1-eta), (1+eta), -(1+eta)]
       dNe = 0.25_dp*[-(1-xi), -(1+xi), (1+xi), (1-xi)]
       xp = matmul(xc, N);  gx = matmul(xc, dNx);  ge = matmul(xc, dNe)
       r = pm - xp
       A(1,1) = sum(gx*gx);  A(1,2) = sum(gx*ge);  A(2,1) = A(1,2);  A(2,2) = sum(ge*ge)
       bb = [sum(gx*r), sum(ge*r)]
       det = A(1,1)*A(2,2) - A(1,2)*A(2,1)
       if (abs(det) < 1.0e-300_dp) exit
       dxi = (A(2,2)*bb(1) - A(1,2)*bb(2))/det;  deta = (A(1,1)*bb(2) - A(2,1)*bb(1))/det
       xi = xi + dxi;  eta = eta + deta
       if (abs(dxi) + abs(deta) < 1.0e-10_dp) exit
       if (abs(xi) > 3.0_dp .or. abs(eta) > 3.0_dp) exit
    end do
    N = 0.25_dp*[(1-xi)*(1-eta), (1+xi)*(1-eta), (1+xi)*(1+eta), (1-xi)*(1+eta)]
    xp = matmul(xc, N)
    dist = sqrt(sum((pm - xp)**2))
  end subroutine inverse_bilinear

  !---------------------------------------------------------------------
  ! morph data of the part's blocks
  !---------------------------------------------------------------------
  subroutine setup_morph(p)
    type(part_t), intent(inout), target :: p
    integer :: is, b, i, j, k, ni, nj, nk, ax, sd, c, q, idx(3), jdx(3), lo, hi, step
    real(dp) :: len_, s_
    type(sheet_t), pointer :: s
    do is = 1, p%nsheet
       s => p%sheet(is)
       b = s%b
       if (mb(b)%on) call fsi_abort('block '//itoa(b)//' has two no-slip wall patches - one flexible sheet per block is supported')
       mb(b)%on = .true.;  mb(b)%part = p%k;  mb(b)%sheet = is
       fsi_blk(b) = .true.
       ni = blk(b)%ni;  nj = blk(b)%nj;  nk = blk(b)%nk
       allocate(mb(b)%X0(3, ni, nj, nk), mb(b)%sarc(ni, nj, nk), mb(b)%xprev(3, ni, nj, nk), &
                mb(b)%vh(3, 2*ni-1, 2*nj-1, 2*nk-1))
       do k = 1, nk
          do j = 1, nj
             do i = 1, ni
                mb(b)%X0(:, i, j, k) = [blk(b)%xn(i,j,k), blk(b)%yn(i,j,k), blk(b)%zn(i,j,k)]
             end do
          end do
       end do
       mb(b)%xprev = mb(b)%X0;  mb(b)%vh = 0.0_dp
       ! arc-length fraction along the columns (from the wall face)
       ax = s%ax;  sd = s%sd
       mb(b)%sarc = 1.0_dp
       do k = 1, nk
          do j = 1, nj
             do i = 1, ni
                idx = [i, j, k]
                if (idx(ax) /= 1) cycle
                ! this column: nodes idx with idx(ax) = 1 .. n(ax), from the wall side
                if (sd == 1) then
                   lo = 1;  hi = nnode(b, ax);  step = 1
                else
                   lo = nnode(b, ax);  hi = 1;  step = -1
                end if
                len_ = 0.0_dp
                jdx = idx;  jdx(ax) = lo
                mb(b)%sarc(jdx(1), jdx(2), jdx(3)) = 0.0_dp
                do q = lo + step, hi, step
                   jdx(ax) = q
                   c = q - step
                   len_ = len_ + sqrt(sum((mb(b)%X0(:, jdx(1), jdx(2), jdx(3)) - &
                          mb(b)%X0(:, merge(c, jdx(1), ax == 1), merge(c, jdx(2), ax == 2), merge(c, jdx(3), ax == 3)))**2))
                   mb(b)%sarc(jdx(1), jdx(2), jdx(3)) = len_
                end do
                if (len_ > 0.0_dp) then
                   do q = lo, hi, step
                      jdx(ax) = q
                      s_ = mb(b)%sarc(jdx(1), jdx(2), jdx(3))/len_
                      mb(b)%sarc(jdx(1), jdx(2), jdx(3)) = s_
                   end do
                end if
             end do
          end do
       end do
    end do
  end subroutine setup_morph

  pure integer function nnode(b, a)
    integer, intent(in) :: b, a
    select case (a)
    case (1);  nnode = blk(b)%ni
    case (2);  nnode = blk(b)%nj
    case default;  nnode = blk(b)%nk
    end select
  end function nnode

  ! blend of the wall rotation along a column.  24b: s_rigid = 1, every column of the wall node
  ! moves rigidly (translation and rotation) - the block bends with the wall like a thick shell and
  ! its overset outer boundary follows; the first version decayed the rotation with a cosine from
  ! 0.25 of the arc length to the outer boundary, which tangled the fan of columns around the
  ! plate's rounded end once the tip rotation reached 0.3 rad (the columns of the cap crossed those
  ! of the flat sides)
  pure real(dp) function blend(s_)
    real(dp), intent(in) :: s_
    if (s_ <= s_rigid) then
       blend = 1.0_dp
    else if (s_ >= 1.0_dp) then
       blend = 0.0_dp
    else
       blend = 0.5_dp*(1.0_dp + cos(pi*(s_ - s_rigid)/(1.0_dp - s_rigid)))
    end if
  end function blend

  !=====================================================================
  ! wall node displacement / rotation of every sheet from the relaxed
  ! interface displacement of its part (part frame)
  !=====================================================================
  subroutine wall_from_structure(p)
    type(part_t), intent(inout), target :: p
    integer :: is, i1, i2, a, nd
    real(dp) :: x(3), th(3), Ra(3,3), w
    type(sheet_t), pointer :: s
    do is = 1, p%nsheet
       s => p%sheet(is)
       do i2 = 1, s%n2
          do i1 = 1, s%n1
             x = 0.0_dp;  th = 0.0_dp
             do a = 1, 4
                nd = s%att(a, i1, i2)
                if (nd == 0) cycle
                w = s%wgt(a, i1, i2)
                Ra = matmul(rot_exp(p%drel(4:6, nd)), p%sh%R(:,:,nd))
                x = x + w*(p%drel(1:3, nd) + matmul(Ra, s%off(:, a, i1, i2)) - s%off(:, a, i1, i2))
                th = th + w*rot_log(Ra)
             end do
             s%dw(:, i1, i2) = x;  s%thw(:, i1, i2) = th
          end do
       end do
    end do
  end subroutine wall_from_structure

  !=====================================================================
  ! HOOK (sixdof.move_grids): morph the blocks of the flexible parts
  ! after the rigid motion, rebuild the metrics of the owned morphing
  ! sub-domains, discrete grid velocities, ALE fluxes
  !=====================================================================
  subroutine fsi_morph(restoring)
    logical, intent(in) :: restoring
    integer :: b, ip, is, i, j, k, ni, nj, nk, ax, idx(3), ib, d, ndeg, h1, h2, h3, q, hx, hy, hz, n1, n2, n3
    real(dp) :: Rg(3,3), O0(3), dg(3), vG(3), aG(3), omG(3), alG(3), x(3), dm(3), xw(3), th(3), Rw(3,3), s_
    real(dp) :: dtt, vn(3), sbar
    real(dp), allocatable :: ths(:,:,:,:), tt(:,:,:), dmb(:,:,:,:), dtmp(:,:,:,:), sbl(:)
    integer :: nl, l, m, it, i1, i2, c1, c2
    type(sheet_t), pointer :: s
    if (.not. inited) return
    do ip = 1, npart
       if (frames_part(pt(ip)%k) > 0) then              ! prescribed deformation (H1): coordinates and grid velocity from the frames
          call frames_morph(pt(ip)%k, restoring)
          cycle
       end if
       call sixdof_frame_kin(pt(ip)%k, Rg, O0, dg, vG, aG, omG, alG)
       do is = 1, pt(ip)%nsheet
          s => pt(ip)%sheet(is)
          b = s%b;  ax = s%ax
          ni = blk(b)%ni;  nj = blk(b)%nj;  nk = blk(b)%nk
          ! the rotation of the columns, level by level away from the wall (24c): the wall rotation
          ! diffused along the wall (around / span) by nint(100 sbar^2) Jacobi sweeps, sbar = the mean
          ! arc-length fraction of the level - the near-wall cells turn rigidly with their wall node,
          ! the outer ones with a smoothed rotation.  Rigid columns with the wall rotation itself
          ! sheared the outer cells between the clamped columns (no motion) and the first rotating
          ! ones of a 20-mm long fan at the plate root: at a 22 mm tip deflection the cell area there
          ! fell to 6 % and the imposed-motion test diverged at full amplitude (with the smoothing
          ! 80 %; 64 % at 30 mm)
          nl = nnode(b, ax)
          allocate(ths(3, s%n1, s%n2, nl), tt(3, s%n1, s%n2), dmb(3, ni, nj, nk), dtmp(3, ni, nj, nk), sbl(nl))
          do l = 1, nl
             idx = 1;  idx(ax) = l
             sbar = 0.0_dp;  q = 0
             do k = 1, merge(1, nk, ax == 3)
                do j = 1, merge(1, nj, ax == 2)
                   do i = 1, merge(1, ni, ax == 1)
                      idx = [i, j, k];  idx(ax) = l
                      sbar = sbar + mb(b)%sarc(idx(1), idx(2), idx(3));  q = q + 1
                   end do
                end do
             end do
             sbar = sbar/max(q, 1);  sbl(l) = sbar
             m = nint(100.0_dp*sbar**2)
             ths(:, :, :, l) = s%thw
             do it = 1, m
                tt = ths(:, :, :, l)
                do i2 = 1, s%n2
                   do i1 = 1, s%n1
                      ! per direction: a field uniform along one direction stays uniform (the span of a
                      ! quasi-2-D case: the k planes of the grid stay identical)
                      th = 0.0_dp
                      c1 = 0;  vn = 0.0_dp
                      call nbr(i1 - 1, i2, 1);  call nbr(i1 + 1, i2, 1)
                      if (c1 > 0) th = th + 0.5_dp*(vn/c1 - tt(:, i1, i2))
                      c1 = 0;  vn = 0.0_dp
                      call nbr(i1, i2 - 1, 2);  call nbr(i1, i2 + 1, 2)
                      if (c1 > 0) th = th + 0.5_dp*(vn/c1 - tt(:, i1, i2))
                      ths(:, i1, i2, l) = tt(:, i1, i2) + th
                   end do
                end do
             end do
          end do
          do k = 1, nk
             do j = 1, nj
                do i = 1, ni
                   idx = [i, j, k]
                   ! the wall node of this column
                   idx(ax) = merge(1, nnode(b, ax), s%sd == 1)
                   n1 = idx(s%t1) - s%rr(1,s%t1) + 1;  n2 = idx(s%t2) - s%rr(1,s%t2) + 1
                   if (n1 < 1 .or. n1 > s%n1 .or. n2 < 1 .or. n2 > s%n2) then
                      dm = 0.0_dp
                   else
                      xw = mb(b)%X0(:, idx(1), idx(2), idx(3))
                      s_ = mb(b)%sarc(i, j, k)
                      l = merge(i, merge(j, k, ax == 2), ax == 1)
                      th = blend(s_)*ths(:, n1, n2, l)
                      Rw = rot_exp(th)
                      dm = s%dw(:, n1, n2) + matmul(Rw, mb(b)%X0(:, i, j, k) - xw) - (mb(b)%X0(:, i, j, k) - xw)
                   end if
                   dmb(:, i, j, k) = dm
                end do
             end do
          end do
          ! the displacement of the outer levels diffused along the wall (24c): nint(200 sbar^2)
          ! Jacobi sweeps along the block's tangential index directions - on the concave side of a
          ! strongly bent plate the columns converge and crushed the outer cells (the free benchmark
          ! run diverged at t = 0.65 s with 6 % of the area left there; with the diffusion 37 %)
          do l = 1, nl
             m = nint(200.0_dp*sbl(l)**2)
             do it = 1, m
                call diffuse_level(l)
             end do
          end do
          do k = 1, nk
             do j = 1, nj
                do i = 1, ni
                   x = matmul(Rg, mb(b)%X0(:, i, j, k) + dmb(:, i, j, k) - O0) + O0 + dg
                   blk(b)%xn(i,j,k) = x(1);  blk(b)%yn(i,j,k) = x(2);  blk(b)%zn(i,j,k) = x(3)
                end do
             end do
          end do
          deallocate(ths, tt, dmb, dtmp, sbl)
          ! discrete velocity at the nodes -> half lattice
          dtt = max(dt_dim, 1.0e-300_dp)
          if (restoring) then
             do k = 1, nk
                do j = 1, nj
                   do i = 1, ni
                      mb(b)%xprev(:, i, j, k) = [blk(b)%xn(i,j,k), blk(b)%yn(i,j,k), blk(b)%zn(i,j,k)]
                   end do
                end do
             end do
             mb(b)%vh = 0.0_dp
          else
             mb(b)%vh = 0.0_dp
             do h3 = 1, 2*nk - 1
                do h2 = 1, 2*nj - 1
                   do h1 = 1, 2*ni - 1
                      vn = 0.0_dp;  q = 0
                      do hz = (h3+1)/2, min((h3+2)/2, nk)
                         do hy = (h2+1)/2, min((h2+2)/2, nj)
                            do hx = (h1+1)/2, min((h1+2)/2, ni)
                               vn = vn + ([blk(b)%xn(hx,hy,hz), blk(b)%yn(hx,hy,hz), blk(b)%zn(hx,hy,hz)] - mb(b)%xprev(:, hx, hy, hz))
                               q = q + 1
                            end do
                         end do
                      end do
                      mb(b)%vh(:, h1, h2, h3) = vn/q/dtt/Uch
                   end do
                end do
             end do
          end if
       end do
    end do
    ! ghost coordinates of every block (replicated), metrics of the owned morphing sub-domains
    call fill_ghost_coords
    do ib = 1, nown
       d = own_list(ib);  b = dom(d)%b
       if (.not. fsi_blk(b)) cycle
       ! the old metrics of the flux re-projection: those of the level-n grid - kept through the
       ! sub-iterations of a step (24b), each of which restores the level-n fluxes and re-projects
       ! them onto its own grid
       if (.not. allocated(om(d)%Jh)) then
          allocate(om(d)%Jh, source=dom(d)%Jh);  allocate(om(d)%mh, source=dom(d)%mh)
       else if (step_start .or. restoring) then
          om(d)%Jh = dom(d)%Jh;  om(d)%mh = dom(d)%mh
       end if
       ndeg = 0
       call metrics_domain(d, ndeg)
       call set_block(d)
       ! ALE grid fluxes from the half-lattice velocity
       do k = lbound(U1,3), ubound(U1,3)
          do j = lbound(U1,2), ubound(U1,2)
             do i = lbound(U1,1), ubound(U1,1)
                do q = 1, 3
                   hx = 2*i;  hy = 2*j;  hz = 2*k
                   select case (q)
                   case (1);  hx = hx + 1
                   case (2);  hy = hy + 1
                   case default;  hz = hz + 1
                   end select
                   hx = min(max(hx, hlo(1)+1), hhi(1)-1)
                   hy = min(max(hy, hlo(2)+1), hhi(2)-1)
                   hz = min(max(hz, hlo(3)+1), hhi(3)-1)
                   x = fsi_hvel(b, [hx, hy, hz])
                   select case (q)
                   case (1);  Ug1(i,j,k) = fvw(hx,hy,hz)*(mq(1,1,hx,hy,hz)*x(1) + mq(2,1,hx,hy,hz)*x(2) + mq(3,1,hx,hy,hz)*x(3))
                   case (2);  Ug2(i,j,k) = fvw(hx,hy,hz)*(mq(1,2,hx,hy,hz)*x(1) + mq(2,2,hx,hy,hz)*x(2) + mq(3,2,hx,hy,hz)*x(3))
                   case default;  Ug3(i,j,k) = fvw(hx,hy,hz)*(mq(1,3,hx,hy,hz)*x(1) + mq(2,3,hx,hy,hz)*x(2) + mq(3,3,hx,hy,hz)*x(3))
                   end select
                end do
             end do
          end do
       end do
    end do
    call sync_face_jacobian(.true.)
    step_start = .false.
    n3 = 0;  vG = 0.0_dp
  contains
    ! one Jacobi sweep of the displacement of level l along the block's tangential directions
    ! (periodic across the seam of the around direction of a periodic sheet)
    subroutine diffuse_level(lv)
      integer, intent(in) :: lv
      integer :: a, bq, ii, jj, kk, nn_, id(3), jd(3), cnt
      real(dp) :: acc(3)
      dtmp = dmb
      do kk = 1, nk
         do jj = 1, nj
            do ii = 1, ni
               id = [ii, jj, kk]
               if (id(ax) /= lv) cycle
               do a = 1, 2                                ! per direction (a field uniform along one stays uniform)
                  acc = 0.0_dp;  cnt = 0
                  bq = merge(s%t1, s%t2, a == 1)
                  nn_ = nnode(b, bq)
                  do c1 = -1, 1, 2
                     jd = id;  jd(bq) = id(bq) + c1
                     if (s%periodic .and. s%ar == a) jd(bq) = modulo(jd(bq) - 1, nn_ - 1) + 1
                     if (jd(bq) < 1 .or. jd(bq) > nn_) cycle
                     acc = acc + dtmp(:, jd(1), jd(2), jd(3));  cnt = cnt + 1
                  end do
                  if (cnt > 0) dmb(:, ii, jj, kk) = dmb(:, ii, jj, kk) + 0.5_dp*(acc/cnt - dtmp(:, ii, jj, kk))
               end do
            end do
         end do
      end do
    end subroutine diffuse_level
    ! neighbour (a1, a2) of the wall-rotation diffusion along direction dir (periodic across the
    ! seam of the around direction of a periodic sheet, otherwise only inside the patch)
    subroutine nbr(a1, a2, dir)
      integer, intent(in) :: a1, a2, dir
      integer :: b1, b2
      b1 = a1;  b2 = a2
      if (s%periodic .and. s%ar == dir) then
         if (dir == 1) b1 = modulo(a1 - 1, s%n1 - 1) + 1
         if (dir == 2) b2 = modulo(a2 - 1, s%n2 - 1) + 1
      end if
      if (b1 < 1 .or. b1 > s%n1 .or. b2 < 1 .or. b2 > s%n2) return
      vn = vn + tt(:, b1, b2);  c1 = c1 + 1
    end subroutine nbr
  end subroutine fsi_morph

  ! level-n fields of every owned sub-domain (saved at the step start, restored before every
  ! sub-iteration after the first: the repeated fluid step then starts from the same state on
  ! every block, the morphed blocks' fluxes being re-projected from the level-n grid)
  subroutine save_level_n
    integer :: ib, d
    do ib = 1, nown
       d = own_list(ib)
       if (.not. allocated(sv(d)%u)) then
          allocate(sv(d)%u, mold=dom(d)%u);   allocate(sv(d)%v, mold=dom(d)%v);   allocate(sv(d)%w, mold=dom(d)%w)
          allocate(sv(d)%U1, mold=dom(d)%U1); allocate(sv(d)%U2, mold=dom(d)%U2); allocate(sv(d)%U3, mold=dom(d)%U3)
          allocate(sv(d)%p, mold=dom(d)%p)
       end if
       sv(d)%u = dom(d)%u;   sv(d)%v = dom(d)%v;   sv(d)%w = dom(d)%w
       sv(d)%U1 = dom(d)%U1; sv(d)%U2 = dom(d)%U2; sv(d)%U3 = dom(d)%U3
       sv(d)%p = dom(d)%p
    end do
  end subroutine save_level_n

  subroutine restore_level_n
    integer :: ib, d
    do ib = 1, nown
       d = own_list(ib)
       if (.not. allocated(sv(d)%u)) cycle
       ! the face fluxes always: the advecting fluxes of the repeated step are those of level n
       ! (the re-projection onto the next grid starts from the level-n fluxes, fsi_flux_reproject)
       dom(d)%U1 = sv(d)%U1; dom(d)%U2 = sv(d)%U2; dom(d)%U3 = sv(d)%U3
       if (fsi_restore_mode == 1) then
          dom(d)%u = sv(d)%u;   dom(d)%v = sv(d)%v;   dom(d)%w = sv(d)%w
          ! the pressure: only with the incremental form (the predictor uses grad p^n); with the full
          ! pressure (Prs_form = 1) it is just the initial guess of the Poisson solve, and the latest
          ! sub-iterate is a far better one than p^n (24c: about half the Krylov iterations of the
          ! repeated steps, the same converged pressure)
          if (Prs_form == 2) dom(d)%p = sv(d)%p
       end if
    end do
  end subroutine restore_level_n

  ! HOOK: discrete grid velocity (non-dim, grid frame) at half index h of morphing block b
  function fsi_hvel(b, h) result(ug)
    integer, intent(in) :: b, h(3)
    real(dp) :: ug(3)
    integer :: hh(3)
    if (frames_blk(b)) then                             ! prescribed deformation (H1)
       ug = frames_hvel(b, h);  return
    end if
    hh(1) = min(max(h(1), 1), 2*blk(b)%ni - 1)
    hh(2) = min(max(h(2), 1), 2*blk(b)%nj - 1)
    hh(3) = min(max(h(3), 1), 2*blk(b)%nk - 1)
    ug = mb(b)%vh(:, hh(1), hh(2), hh(3))
  end function fsi_hvel

  ! HOOK (sixdof.after_motion): the stored face fluxes of a morphed sub-domain with the new metrics,
  ! U_new = U_old + fvw_new (m_new . u) - fvw_old (m_old . u), the Cartesian face velocity carried
  subroutine fsi_flux_reproject(d)
    integer, intent(in) :: d
    integer :: i, j, k, lo(3), hi(3), h1, h2, h3
    real(dp) :: uc(3), fo, fn
    real(dp), pointer :: u(:,:,:), v(:,:,:), w(:,:,:)
    if (.not. allocated(om(d)%Jh)) return
    call set_block(d)
    ! the Cartesian velocity carried: that of level n, the one the restored fluxes belong to
    ! (saved by fsi_predict; the current field is the latest sub-iterate in restore mode 2)
    if (allocated(sv(d)%u)) then
       u => sv(d)%u;  v => sv(d)%v;  w => sv(d)%w
    else
       u => dom(d)%u;  v => dom(d)%v;  w => dom(d)%w
    end if
    lo = [lbound(U1,1), lbound(U1,2), lbound(U1,3)];  hi = [ubound(U1,1), ubound(U1,2), ubound(U1,3)]
    do k = lo(3)+1, hi(3)-1
     do j = lo(2)+1, hi(2)-1
      do i = lo(1), hi(1)-1
         h1 = 2*i+1;  h2 = 2*j;  h3 = 2*k
         if (h1 <= hlo(1) .or. h1 >= hhi(1)) cycle
         uc(1) = u(i,j,k)
         uc(2) = 0.25_dp*(v(i,j-1,k) + v(i+1,j-1,k) + v(i,j,k) + v(i+1,j,k))
         uc(3) = 0.25_dp*(w(i,j,k-1) + w(i+1,j,k-1) + w(i,j,k) + w(i+1,j,k))
         fo = fvw_of(om(d)%Jh(h1,h2,h3));  fn = fvw(h1,h2,h3)
         U1(i,j,k) = U1(i,j,k) + fn*sum(mh(1:3,h1,h2,h3)*uc) - fo*sum(om(d)%mh(1:3,h1,h2,h3)*uc)
      end do
     end do
    end do
    do k = lo(3)+1, hi(3)-1
     do j = lo(2), hi(2)-1
      do i = lo(1)+1, hi(1)-1
         h1 = 2*i;  h2 = 2*j+1;  h3 = 2*k
         if (h2 <= hlo(2) .or. h2 >= hhi(2)) cycle
         uc(2) = v(i,j,k)
         uc(1) = 0.25_dp*(u(i-1,j,k) + u(i-1,j+1,k) + u(i,j,k) + u(i,j+1,k))
         uc(3) = 0.25_dp*(w(i,j,k-1) + w(i,j+1,k-1) + w(i,j,k) + w(i,j+1,k))
         fo = fvw_of(om(d)%Jh(h1,h2,h3));  fn = fvw(h1,h2,h3)
         U2(i,j,k) = U2(i,j,k) + fn*sum(mh(4:6,h1,h2,h3)*uc) - fo*sum(om(d)%mh(4:6,h1,h2,h3)*uc)
      end do
     end do
    end do
    do k = lo(3), hi(3)-1
     do j = lo(2)+1, hi(2)-1
      do i = lo(1)+1, hi(1)-1
         h1 = 2*i;  h2 = 2*j;  h3 = 2*k+1
         if (h3 <= hlo(3) .or. h3 >= hhi(3)) cycle
         uc(3) = w(i,j,k)
         uc(1) = 0.25_dp*(u(i-1,j,k) + u(i-1,j,k+1) + u(i,j,k) + u(i,j,k+1))
         uc(2) = 0.25_dp*(v(i,j-1,k) + v(i,j-1,k+1) + v(i,j,k) + v(i,j,k+1))
         fo = fvw_of(om(d)%Jh(h1,h2,h3));  fn = fvw(h1,h2,h3)
         U3(i,j,k) = U3(i,j,k) + fn*sum(mh(7:9,h1,h2,h3)*uc) - fo*sum(om(d)%mh(7:9,h1,h2,h3)*uc)
      end do
     end do
    end do
  contains
    ! the face-variable weight of the OLD metrics (global_vars.fvw evaluates the current Jh)
    pure real(dp) function fvw_of(jv)
      real(dp), intent(in) :: jv
      if (curvi_form == 2) then
         fvw_of = jv
         if (fvw_of == 0.0_dp) fvw_of = 1.0_dp
      else
         fvw_of = 1.0_dp
      end if
    end function fvw_of
  end subroutine fsi_flux_reproject

  !=====================================================================
  ! Step start (before sixdof_advance): predictor of the interface
  ! displacement, activity, re-assembly request
  !=====================================================================
  subroutine fsi_predict(istep, t)
    integer,  intent(in) :: istep
    real(dp), intent(in) :: t
    integer :: ip, i
    real(dp) :: dmax
    type(part_t), pointer :: p
    if (.not. inited) return
    step_start = .true.
    if (any(t >= pt(1:npart)%def%start_t - 1.0e-12_dp)) call save_level_n
    do ip = 1, npart
       p => pt(ip)
       if (frames_part(p%k) > 0) then                   ! prescribed deformation (H1): no structure, no coupling
          p%active = .false.;  p%converged = .true.;  p%iter = 0
          cycle
       end if
       p%active = (t >= p%def%start_t - 1.0e-12_dp)
       p%iter = 0;  p%converged = .false.;  p%omega = p%def%relax;  p%rprev = 0.0_dp
       if (p%active .and. p%def%mode == 1) then
          ! the DISCRETE velocities of the last two steps (24c: the generalized-alpha velocity can be
          ! far from the displacement increment - at the release of the pluck it jumped to 0.9 m/s
          ! while the tip moved 0.08 mm, and the predicted grid motion blew the fluid step up)
          do i = 1, p%sh%nn
             p%drel(1:3, i) = p%sh%u(1:3, i) + dt_dim*(1.5_dp*p%vd(1:3, i) - 0.5_dp*p%vprev(1:3, i))
             p%drel(4:6, i) = dt_dim*(1.5_dp*p%vd(4:6, i) - 0.5_dp*p%vprev(4:6, i))
          end do
          if (.not. has_prev_vel) then
             do i = 1, p%sh%nn
                p%drel(1:3, i) = p%sh%u(1:3, i) + dt_dim*p%vd(1:3, i)
                p%drel(4:6, i) = dt_dim*p%vd(4:6, i)
             end do
          end if
       else
          p%drel(1:3, :) = p%sh%u(1:3, :);  p%drel(4:6, :) = 0.0_dp
       end if
       if (p%active .and. testing(p, t)) call test_predict(p, t)
       call wall_from_structure(p)
       dmax = maxval(abs(p%drel(1:3,:) - p%dasm(1:3,:)))
       if (dmax >= ovst_reasm_tol) then
          fsi_reasm_request = .true.;  p%dasm = p%drel
       end if
    end do
    i = istep
  end subroutine fsi_predict

  !---------------------------------------------------------------------
  ! FSI_Test_Motion = A, f [, mode [, shape]] (verification / start-up
  ! aid, update 24c): the structure is not solved; its nodes follow the
  ! imposed bending mode
  !   d = A r xi^2 sin(w tau) e_n,  theta = A r (2 xi / L) sin(w tau) e_s,
  ! tau = t - FSI_Start_T, r = min(1, f tau) (one-period ramp), xi = the
  ! distance from the clamped nodes along the chord / L (0 at the clamp,
  ! 1 at the free end), e_s = the FSI_Plane_Strain axis, e_n = e_s x chord.
  ! The fluid then sees a known wall motion (amplitude and speed of the
  ! benchmark) without any structural feedback: a test of the moving-grid
  ! treatment of the fluid alone (fluid loads are computed and written).
  ! 3rd value: 0 = the oscillation (default), 1 = ramp to the static
  ! deflection A xi^2 over one period and hold, 2 = PLUCK: the oscillation
  ! d = A r phi sin(w tau) in the shape phi (4th value: 0 = the static
  ! tip-load shape of the structure, m >= 1 = its m-th vibration mode in
  ! vacuum), applied inextensibly, up to its maximum at tau = 2.25 / f
  ! (zero velocity), then the structure is released and solved - a start
  ! of a flow-induced oscillation from a finite amplitude and a developed
  ! flapping wake instead of from the round-off asymmetry of the wake.
  !---------------------------------------------------------------------
  subroutine test_setup(p)
    type(part_t), intent(inout) :: p
    integer :: i, q, n
    real(dp) :: xc(3), cov(3,3), ec(3), sc, smax, sm, x(3)
    if (allocated(p%txi)) return
    if (p%def%plane < 1) then
       if (myrank == 0) write(*,'(A)') ' ERROR: FSI_Test_Motion needs FSI_Plane_Strain (the span axis of the bending mode)'
       call MPI_Abort(MPI_COMM_WORLD, 71, i)
    end if
    allocate(p%txi(p%sh%nn))
    xc = 0.0_dp;  n = 0
    do i = 1, p%sh%nn
       if (all(p%sh%fixed(:, i))) then
          xc = xc + p%sh%X0(:, i);  n = n + 1
       end if
    end do
    if (n > 0) xc = xc/n
    if (n == 0) xc = sum(p%sh%X0, 2)/p%sh%nn
    p%txc = xc
    cov = 0.0_dp
    x = sum(p%sh%X0, 2)/p%sh%nn                  ! the chord: principal axis about the centroid
    do i = 1, p%sh%nn
       ec = p%sh%X0(:, i) - x
       do q = 1, 3
          cov(:, q) = cov(:, q) + ec*ec(q)
       end do
    end do
    call principal_axis(cov, ec)
    p%tes = 0.0_dp;  p%tes(p%def%plane) = 1.0_dp
    ec = ec - sum(ec*p%tes)*p%tes;  ec = ec/sqrt(sum(ec*ec))
    sm = 0.0_dp
    do i = 1, p%sh%nn
       sm = sm + sum((p%sh%X0(:, i) - xc)*ec)
    end do
    if (sm < 0.0_dp) ec = -ec
    sc = -huge(1.0_dp);  smax = -huge(1.0_dp)
    do i = 1, p%sh%nn
       x(1) = sum((p%sh%X0(:, i) - xc)*ec)
       smax = max(smax, x(1))
       if (all(p%sh%fixed(:, i))) sc = max(sc, x(1))
    end do
    if (sc < -0.5_dp*huge(1.0_dp)) sc = -smax
    p%tL = max(smax - sc, 1.0e-12_dp)
    do i = 1, p%sh%nn
       p%txi(i) = max(0.0_dp, (sum((p%sh%X0(:, i) - xc)*ec) - sc)/p%tL)
    end do
    p%ten = cross3(p%tes, ec)
    if (myrank == 0) write(*,'(A,ES10.3,A,F8.3,A,ES10.3,A,3F7.3)') '    TEST (FSI_Test_Motion): imposed bending mode, tip amplitude ', &
       p%def%test_amp, ' m at ', p%def%test_freq, ' Hz (one-period ramp), free length ', p%tL, ' m, direction', p%ten
    if (myrank == 0 .and. p%def%test_mode == 1) write(*,'(A)') '    TEST: mode 1 - ramp to the static deflection over one period, then held'
  end subroutine test_setup

  ! imposed state at time t: translation u, rotation angle th (about e_s), and their rates
  subroutine test_state(p, t, i, u, th, v, thd)
    type(part_t), intent(in) :: p
    real(dp), intent(in) :: t
    integer, intent(in) :: i
    real(dp), intent(out) :: u(3), th, v(3), thd
    real(dp) :: tau, w, r, rd, a
    tau = max(t - p%def%start_t, 0.0_dp);  w = 2.0_dp*pi*p%def%test_freq
    r = min(1.0_dp, p%def%test_freq*tau);  rd = merge(p%def%test_freq, 0.0_dp, p%def%test_freq*tau < 1.0_dp)
    a = p%def%test_amp
    if (p%def%test_mode == 1) then              ! ramp to the static deflection A xi^2 over one period, then hold
       u = a*r*p%txi(i)**2*p%ten;  v = a*rd*p%txi(i)**2*p%ten
       th = a*r*(2.0_dp*p%txi(i)/p%tL);  thd = a*rd*(2.0_dp*p%txi(i)/p%tL)
       return
    end if
    if (p%def%test_mode == 2) then              ! the pluck: the shape phi, inextensible
       call pluck_state(p, i, a*r*sin(w*tau), a*(rd*sin(w*tau) + r*w*cos(w*tau)), u, th, v, thd)
       return
    end if
    u = a*r*p%txi(i)**2*sin(w*tau)*p%ten
    v = a*p%txi(i)**2*(rd*sin(w*tau) + r*w*cos(w*tau))*p%ten
    th = a*r*(2.0_dp*p%txi(i)/p%tL)*sin(w*tau)
    thd = a*(2.0_dp*p%txi(i)/p%tL)*(rd*sin(w*tau) + r*w*cos(w*tau))
  end subroutine test_state

  subroutine test_predict(p, t)
    type(part_t), intent(inout) :: p
    real(dp), intent(in) :: t
    integer :: i
    real(dp) :: u(3), th, v(3), thd, u0(3), th0, v0(3), thd0
    call test_setup(p)
    if (p%def%test_mode == 2) call test_shape(p)
    do i = 1, p%sh%nn
       call test_state(p, t, i, u, th, v, thd)
       call test_state(p, t - dt_dim, i, u0, th0, v0, thd0)
       p%drel(1:3, i) = u;  p%drel(4:6, i) = (th - th0)*p%tes
    end do
  end subroutine test_predict

  ! the structure follows the imposed motion in the step ending at t: FSI_Test_Motion set, and for
  ! mode 2 (the pluck) up to its release at 2.25 periods (a maximum: deflection A phi, zero velocity)
  logical function testing(p, t)
    type(part_t), intent(in) :: p
    real(dp), intent(in) :: t
    testing = (p%def%test_amp > 0.0_dp)
    if (testing .and. p%def%test_mode == 2) &
       testing = (t - p%def%start_t < 2.25_dp/p%def%test_freq - 0.5_dp*dt_dim)
  end function testing

  ! FSI_Test_Motion mode 2: the shape phi of the pluck, normalized to a unit tip deflection along
  ! e_n - the linear static deflection of the structure under a tip load (4th value 0; a copy of
  ! the shell solved statically) or its m-th vibration mode in vacuum (4th value m, shell_mode).
  ! The released structure then carries only the elastic force of that shape: the mode-0 shape
  ! A xi^2 also bends the steel rear mass, whose elastic force snapped the released structure (24c:
  ! the tip velocity jumped to 1.7 m/s in the first free step and the run broke down in the second).
  subroutine test_shape(p)
    type(part_t), intent(inout) :: p
    type(shell_t) :: tmp
    real(dp), allocatable :: f(:,:)
    integer, allocatable :: cnt(:)
    real(dp), allocatable :: phi(:,:)
    real(dp) :: ft, del, xm, sv_, fmode
    integer :: i, it, ntip, ns, j, k
    logical :: conv
    if (allocated(p%tphu)) return
    allocate(p%tphu(3, p%sh%nn), p%tpht(p%sh%nn), f(6, p%sh%nn))
    xm = maxval(p%txi)
    ntip = count(p%txi >= xm - 1.0e-6_dp)
    ft = 1.0e-6_dp;  del = 0.0_dp
    do it = 1, 4
       tmp = p%sh
       tmp%u1 = 0.0_dp;  tmp%dth = 0.0_dp           ! from the reference configuration (also after a restart)
       do i = 1, p%sh%nn
          tmp%R(:, :, i) = reshape([1,0,0, 0,1,0, 0,0,1], [3,3])*1.0_dp
       end do
       f = 0.0_dp
       do i = 1, p%sh%nn
          if (p%txi(i) >= xm - 1.0e-6_dp) f(1:3, i) = ft/ntip*p%ten
       end do
       call shell_static(tmp, f, conv)
       del = 0.0_dp
       do i = 1, p%sh%nn
          if (p%txi(i) >= xm - 1.0e-6_dp) del = max(del, sum(tmp%u1(1:3, i)*p%ten))
       end do
       if (del > 1.0e-7_dp*p%tL .and. del < 1.0e-3_dp*p%tL) exit
       ft = ft*1.0e-4_dp*p%tL/max(abs(del), 1.0e-30_dp)
    end do
    del = max(del, 1.0e-30_dp)
    do i = 1, p%sh%nn
       p%tphu(:, i) = tmp%u1(1:3, i)/del
       p%tpht(i) = sum(tmp%dth(:, i)*p%tes)/del
    end do
    if (myrank == 0) write(*,'(A,ES10.3,A,F8.4,A)') '    TEST: mode 2 (pluck) - static stiffness at the tip ', &
       ft/del, ' N/m; the structure is released at t = ', p%def%start_t + 2.25_dp/p%def%test_freq, &
       ' s (a maximum, zero velocity) and solved from there on'
    if (p%def%test_shp > 0) then                    ! the m-th vibration mode in vacuum instead
       allocate(phi(6, p%sh%nn))
       tmp = p%sh
       do i = 1, p%sh%nn
          tmp%R(:, :, i) = reshape([1,0,0, 0,1,0, 0,0,1], [3,3])*1.0_dp
       end do
       call shell_mode(tmp, p%def%test_shp, phi, fmode)
       ! normalized to a unit MAXIMUM deflection along e_n (a higher mode has its largest deflection
       ! inside the free length - normalized by the tip, the 2nd mode of the benchmark moved the rubber
       ! by 1.4 x the tip amplitude in the opposite direction and folded the morphed grid), the sign
       ! such that the tip moves along +e_n
       del = 0.0_dp;  ft = 0.0_dp
       do i = 1, p%sh%nn
          if (abs(sum(phi(1:3, i)*p%ten)) > abs(del)) del = abs(sum(phi(1:3, i)*p%ten))
          if (p%txi(i) >= xm - 1.0e-6_dp .and. abs(sum(phi(1:3, i)*p%ten)) > abs(ft)) ft = sum(phi(1:3, i)*p%ten)
       end do
       if (ft < 0.0_dp) del = -del
       if (abs(del) < 1.0e-30_dp) del = 1.0_dp
       do i = 1, p%sh%nn
          p%tphu(:, i) = phi(1:3, i)/del
          p%tpht(i) = sum(phi(4:6, i)*p%tes)/del
       end do
       if (myrank == 0) write(*,'(A,I0,A,F9.3,A)') '    TEST: pluck shape = vibration mode ', p%def%test_shp, ' in vacuum (', fmode, &
          ' Hz), the amplitude = its largest deflection'
       deallocate(phi)
    else
       if (myrank == 0) write(*,'(A)') '    TEST: pluck shape = the static deflection under a tip load'
    end if
    ! the chord stations of the nodes (sorted, unique) and the rotation of the shape at each
    ns = 0
    allocate(p%tsg(p%sh%nn), p%tthg(p%sh%nn), p%tix(p%sh%nn), p%th0(p%sh%nn), cnt(p%sh%nn))
    do i = 1, p%sh%nn
       sv_ = p%txi(i)*p%tL;  k = 0
       do j = 1, ns
          if (abs(p%tsg(j) - sv_) <= 1.0e-9_dp*p%tL) k = j
       end do
       if (k == 0) then
          ns = ns + 1;  k = ns;  p%tsg(k) = sv_;  p%tthg(k) = 0.0_dp;  cnt(k) = 0
       end if
       p%tix(i) = k;  p%tthg(k) = p%tthg(k) + p%tpht(i);  cnt(k) = cnt(k) + 1
    end do
    p%tthg(1:ns) = p%tthg(1:ns)/cnt(1:ns)
    do j = 1, ns                                     ! sort by station (exchange), carrying the node indices
       do k = j + 1, ns
          if (p%tsg(k) < p%tsg(j)) then
             sv_ = p%tsg(j);  p%tsg(j) = p%tsg(k);  p%tsg(k) = sv_
             sv_ = p%tthg(j);  p%tthg(j) = p%tthg(k);  p%tthg(k) = sv_
             do i = 1, p%sh%nn
                if (p%tix(i) == j) then
                   p%tix(i) = -k
                else if (p%tix(i) == k) then
                   p%tix(i) = j
                end if
             end do
             p%tix = abs(p%tix)
          end if
       end do
    end do
    p%tns = ns
    do i = 1, p%sh%nn                                ! offset of the node from the chord line along e_n
       p%th0(i) = sum((p%sh%X0(:, i) - p%txc)*p%ten)
    end do
    deallocate(cnt, f)
    ! the shape for the record: chord station, deflection along e_n and rotation per unit amplitude
    if (myrank == 0) then
       open(newunit=k, file=out_dir//trim(casename)//'_fsi_pluck_shape.dat', status='replace')
       write(k, '(A)') '# FSI_Test_Motion pluck shape (24c): node, xi (distance from the clamp / free length), phi_n, theta_phi [rad per m]'
       do i = 1, p%sh%nn
          write(k, '(I6, 3ES14.5)') i, p%txi(i), sum(p%tphu(:, i)*p%ten), p%tpht(i)
       end do
       close(k)
    end if
  end subroutine test_shape

  ! the pluck state of node i for the tip-deflection factor c (rate cd): the rotation of the
  ! shape, theta(s) = c theta_phi(s), is integrated along the chord - the INEXTENSIBLE deformed
  ! line P(s) = int_0^s (cos theta e_c + sin theta e_n) ds, the node at its offset h from the line
  ! along the rotated normal.  (The linear shape c phi itself stretches the line by ~ w'^2 / 2:
  ! 2 % at 11 mm, 4e9 Pa in the steel rear mass, which snapped the released structure - 24c.)
  subroutine pluck_state(p, i, c, cd, u, th, v, thd)
    type(part_t), intent(in) :: p
    integer, intent(in) :: i
    real(dp), intent(in) :: c, cd
    real(dp), intent(out) :: u(3), th, v(3), thd
    real(dp) :: pl(3), pld(3), ec(3), t0, t1, ds, sn, cs, thl
    integer :: j
    ec = cross3(p%ten, p%tes)                        ! the chord direction (e_n = e_s x e_c)
    pl = 0.0_dp;  pld = 0.0_dp
    do j = 2, p%tix(i)
       ds = p%tsg(j) - p%tsg(j-1)
       t0 = c*p%tthg(j-1);  t1 = c*p%tthg(j)
       pl = pl + 0.5_dp*ds*((cos(t0) + cos(t1))*ec + (sin(t0) + sin(t1))*p%ten)
       pld = pld + 0.5_dp*ds*cd*(p%tthg(j-1)*(-sin(t0)*ec + cos(t0)*p%ten) + p%tthg(j)*(-sin(t1)*ec + cos(t1)*p%ten))
    end do
    pl = pl + p%tsg(1)*ec                            ! the first station (the clamp) does not move
    thl = c*p%tthg(p%tix(i));  sn = sin(thl);  cs = cos(thl)
    u = pl - p%tsg(p%tix(i))*ec + p%th0(i)*((cs - 1.0_dp)*p%ten - sn*ec)
    v = pld + p%th0(i)*cd*p%tthg(p%tix(i))*(-sn*p%ten - cs*ec)
    th = thl;  thd = cd*p%tthg(p%tix(i))          ! the station's rotation: the same on every span node
  end subroutine pluck_state

  ! the imposed trial state of the step (in place of the structural solve)
  subroutine test_trial(p, t)
    type(part_t), intent(inout) :: p
    real(dp), intent(in) :: t
    integer :: i
    real(dp) :: u(3), th, v(3), thd, u0(3), th0, v0(3), thd0
    call test_setup(p)
    if (p%def%test_mode == 2) call test_shape(p)
    do i = 1, p%sh%nn
       call test_state(p, t, i, u, th, v, thd)
       call test_state(p, t - dt_dim, i, u0, th0, v0, thd0)
       p%sh%u1(1:3, i) = u
       p%sh%dth(:, i) = (th - th0)*p%tes
       p%sh%R1(:, :, i) = matmul(rot_exp(p%sh%dth(:, i)), p%sh%R(:, :, i))
       p%sh%v1(1:3, i) = v;  p%sh%v1(4:6, i) = thd*p%tes
       ! the acceleration of the imposed motion (backward difference): the structural state the
       ! pluck releases (mode 2) is then kinematically consistent (0 there disturbed the released
       ! generalized-alpha step)
       p%sh%a1(1:3, i) = (v - v0)/dt_dim;  p%sh%a1(4:6, i) = (thd - thd0)/dt_dim*p%tes
    end do
    p%sh%newton_it = 0;  p%sh%newton_res = 0.0_dp;  p%sh%substeps = 1
  end subroutine test_trial

  ! FSI_Plane_Strain (24c): the held components of a restored structural state are zeroed (the span
  ! displacement, the rotations about the two other axes; the triads keep the span-axis rotation)
  subroutine plane_project(p)
    type(part_t), intent(inout) :: p
    integer :: i, q
    real(dp) :: th(3)
    if (p%def%plane < 1) return
    do i = 1, p%sh%nn
       p%sh%u(p%def%plane, i) = 0.0_dp;  p%sh%v(p%def%plane, i) = 0.0_dp;  p%sh%a(p%def%plane, i) = 0.0_dp
       do q = 1, 3
          if (q /= p%def%plane) then
             p%sh%u(3 + q, i) = 0.0_dp;  p%sh%v(3 + q, i) = 0.0_dp;  p%sh%a(3 + q, i) = 0.0_dp
          end if
       end do
       th = rot_log(p%sh%R(:, :, i))
       do q = 1, 3
          if (q /= p%def%plane) th(q) = 0.0_dp
       end do
       p%sh%R(:, :, i) = rot_exp(th)
    end do
  end subroutine plane_project

  ! FSI_Plane_Strain (24c): the span groups - nodes whose coordinates differ only along the span axis
  subroutine plane_groups(p)
    type(part_t), intent(inout) :: p
    integer :: i, j, q
    real(dp) :: tolg, d
    allocate(p%pgrp(p%sh%nn));  p%pgrp = 0;  p%ngrp = 0
    tolg = 1.0e-6_dp*maxval(abs(p%sh%X0))
    do i = 1, p%sh%nn
       if (p%pgrp(i) > 0) cycle
       p%ngrp = p%ngrp + 1;  p%pgrp(i) = p%ngrp
       do j = i + 1, p%sh%nn
          if (p%pgrp(j) > 0) cycle
          d = 0.0_dp
          do q = 1, 3
             if (q /= p%def%plane) d = max(d, abs(p%sh%X0(q, i) - p%sh%X0(q, j)))
          end do
          if (d <= tolg) p%pgrp(j) = p%ngrp
       end do
    end do
    if (myrank == 0) write(*,'(A,I0,A,I0,A)') '    FSI_Plane_Strain: ', p%sh%nn, ' nodes in ', p%ngrp, &
       ' span groups (the trial states are averaged over every group: exact cylindrical bending)'
  end subroutine plane_groups

  ! the trial state averaged over every span group (24c: solved independently, the span nodes of the
  ! free benchmark run drifted apart by 0.1 mm in 0.3 s and the fluid grid was no longer extruded)
  subroutine plane_average(p)
    type(part_t), intent(inout) :: p
    integer :: g, i
    real(dp), allocatable :: su(:,:), sv_(:,:), sa(:,:), sd(:,:)
    integer, allocatable :: cnt(:)
    if (p%def%plane < 1 .or. .not. allocated(p%pgrp)) return
    allocate(su(6, p%ngrp), sv_(6, p%ngrp), sa(6, p%ngrp), sd(3, p%ngrp), cnt(p%ngrp))
    su = 0.0_dp;  sv_ = 0.0_dp;  sa = 0.0_dp;  sd = 0.0_dp;  cnt = 0
    do i = 1, p%sh%nn
       g = p%pgrp(i)
       su(:, g) = su(:, g) + p%sh%u1(:, i);  sv_(:, g) = sv_(:, g) + p%sh%v1(:, i)
       sa(:, g) = sa(:, g) + p%sh%a1(:, i);  sd(:, g) = sd(:, g) + p%sh%dth(:, i);  cnt(g) = cnt(g) + 1
    end do
    do i = 1, p%sh%nn
       g = p%pgrp(i)
       p%sh%u1(:, i) = su(:, g)/cnt(g);  p%sh%v1(:, i) = sv_(:, g)/cnt(g);  p%sh%a1(:, i) = sa(:, g)/cnt(g)
       p%sh%dth(:, i) = sd(:, g)/cnt(g)
       p%sh%R1(:, :, i) = matmul(rot_exp(p%sh%dth(:, i)), p%sh%R(:, :, i))
    end do
    deallocate(su, sv_, sa, sd, cnt)
  end subroutine plane_average

  ! the interface displacement averaged over every span group (restart)
  subroutine plane_average_drel(p)
    type(part_t), intent(inout) :: p
    integer :: g, i
    real(dp), allocatable :: sd(:,:)
    integer, allocatable :: cnt(:)
    if (p%def%plane < 1 .or. .not. allocated(p%pgrp)) return
    allocate(sd(6, p%ngrp), cnt(p%ngrp));  sd = 0.0_dp;  cnt = 0
    do i = 1, p%sh%nn
       g = p%pgrp(i);  sd(:, g) = sd(:, g) + p%drel(:, i);  cnt(g) = cnt(g) + 1
    end do
    do i = 1, p%sh%nn
       g = p%pgrp(i);  p%drel(:, i) = sd(:, g)/cnt(g)
    end do
    deallocate(sd, cnt)
  end subroutine plane_average_drel

  ! the same for a restored state (restart): u, v, a and the triads (their rotation about the span axis)
  subroutine plane_average_state(p)
    type(part_t), intent(inout) :: p
    integer :: g, i
    real(dp), allocatable :: su(:,:), sv_(:,:), sa(:,:), st(:,:)
    integer, allocatable :: cnt(:)
    if (p%def%plane < 1 .or. .not. allocated(p%pgrp)) return
    allocate(su(6, p%ngrp), sv_(6, p%ngrp), sa(6, p%ngrp), st(3, p%ngrp), cnt(p%ngrp))
    su = 0.0_dp;  sv_ = 0.0_dp;  sa = 0.0_dp;  st = 0.0_dp;  cnt = 0
    do i = 1, p%sh%nn
       g = p%pgrp(i)
       su(:, g) = su(:, g) + p%sh%u(:, i);  sv_(:, g) = sv_(:, g) + p%sh%v(:, i)
       sa(:, g) = sa(:, g) + p%sh%a(:, i);  st(:, g) = st(:, g) + rot_log(p%sh%R(:, :, i));  cnt(g) = cnt(g) + 1
    end do
    do i = 1, p%sh%nn
       g = p%pgrp(i)
       p%sh%u(:, i) = su(:, g)/cnt(g);  p%sh%v(:, i) = sv_(:, g)/cnt(g);  p%sh%a(:, i) = sa(:, g)/cnt(g)
       p%sh%R(:, :, i) = rot_exp(st(:, g)/cnt(g))
    end do
    deallocate(su, sv_, sa, st, cnt)
  end subroutine plane_average_state

  ! debugging aid (24c): PROBCFD_DUMP_GRID=<step> writes the grid of that step (all blocks,
  ! formatted PLOT3D, double precision) to <case>_grid_NNNNNNNN.grd
  subroutine dump_grid(istep)
    integer, intent(in) :: istep
    character(len=32) :: ev
    character(len=256) :: fn
    integer :: st, ln, ios, b, iu, i, j, k
    call get_environment_variable('PROBCFD_DUMP_GRID', ev, ln, st)
    if (st /= 0 .or. ln == 0) return
    read(ev, *, iostat=ios) i
    if (ios /= 0 .or. i /= istep .or. myrank /= 0) return
    write(fn, '(A,"_grid_",I8.8,".grd")') trim(casename), istep
    open(newunit=iu, file=trim(fn), form='formatted', status='replace')
    write(iu, *) nblk
    write(iu, *) (blk(b)%ni, blk(b)%nj, blk(b)%nk, b = 1, nblk)
    do b = 1, nblk
       write(iu, '(4ES24.16)') (((blk(b)%xn(i,j,k), i = 1, blk(b)%ni), j = 1, blk(b)%nj), k = 1, blk(b)%nk)
       write(iu, '(4ES24.16)') (((blk(b)%yn(i,j,k), i = 1, blk(b)%ni), j = 1, blk(b)%nj), k = 1, blk(b)%nk)
       write(iu, '(4ES24.16)') (((blk(b)%zn(i,j,k), i = 1, blk(b)%ni), j = 1, blk(b)%nj), k = 1, blk(b)%nk)
    end do
    close(iu)
    write(*,'(A)') '  wrote the grid of this step: '//trim(fn)//' (PROBCFD_DUMP_GRID)'
  end subroutine dump_grid

  !=====================================================================
  ! Inner-iteration hook: wall forces -> structure -> relaxation ->
  ! convergence; the grid is morphed again unless converged (or the last
  ! permitted iteration).  all_done = every part converged.
  !=====================================================================
  subroutine fsi_iterate(istep, t, m, all_done, last)
    integer,  intent(in)  :: istep
    real(dp), intent(in)  :: t
    integer,  intent(in)  :: m
    logical,  intent(out) :: all_done
    logical,  intent(in)  :: last
    integer :: ip, i, c, ierr
    real(dp) :: rinf, dinf, dmax
    real(dp), allocatable :: dtil(:,:), rcur(:,:)
    logical :: conv, remesh
    type(part_t), pointer :: p
    all_done = .true.
    if (.not. inited) return
    call wall_forces
    remesh = .false.
    do ip = 1, npart
       p => pt(ip)
       if (.not. p%active) cycle
       if (p%converged) cycle
       p%iter = p%iter + 1
       call build_loads(p)
       if (testing(p, t)) then
          call test_trial(p, t);  conv = .true.
       else if (p%def%mode == 1) then
          call shell_step(p%sh, dt_dim, p%fext, conv)
       else
          call shell_static(p%sh, p%fext, conv)
       end if
       call plane_average(p)
       p%nnewton = p%nnewton + abs(p%sh%newton_it)
       if (conv .and. p%sh%substeps > 1 .and. myrank == 0) write(*,'(A,I0,A,I0,A)') '  NOTE (FSI): '//sixdof_part_name(p%k)// &
          ' structural step ', istep, ' solved in ', p%sh%substeps, ' sub-steps (the Newton iterations of the whole step failed)'
       if (.not. conv .and. myrank == 0) write(*,'(A,I0,A,I0,A,ES9.2,A)') '  WARNING (FSI): '//sixdof_part_name(p%k)// &
          ' structural Newton not converged at step ', istep, ' (', p%sh%newton_it, ' its, residual ', p%sh%newton_res, ')'
       allocate(dtil(6, p%sh%nn), rcur(6, p%sh%nn))
       dtil(1:3, :) = p%sh%u1(1:3, :);  dtil(4:6, :) = p%sh%dth
       rcur = dtil - p%drel
       rinf = maxval(abs(rcur(1:3, :)))
       dinf = max(maxval(abs(dtil(1:3, :) - p%sh%u(1:3, :))), 1.0e-14_dp)
       p%rprev = rcur
       p%resid = rinf/dinf
       if (p%iter == 1) p%resid0 = p%resid
       ! converged: relative to the displacement change of the step (Fabbri et al. Eq. 23), with the
       ! absolute floor FSI_Tol x shell thickness (a smaller change of the wall position is below
       ! what the fluid solution resolves, and the pressure-solver tolerance leaves load noise of
       ! that order on a very compliant structure)
       p%converged = (p%resid <= p%def%tol .or. rinf <= p%def%tol*p%sh%tmin)
       if (.not. p%converged) then
          ! interface quasi-Newton (IQN-ILS, Degroote et al. 2009): the first sub-iteration relaxes
          ! with FSI_Relax, the later ones take the least-squares secant step on the history of
          ! (residual, structural output) pairs of this time step
          call qn_update(p, rcur, dtil)
       else
          p%drel = dtil
       end if
       p%dmax = maxval(abs(dtil(1:3, :)))
       if (myrank == 0 .and. log_interval == 1) write(*,'(A,I3,A,ES10.3,A,ES10.3,A,F6.3,A,I3)') '  [ FSI ] '//sixdof_part_name(p%k)// &
          ' sub-iteration ', p%iter, ': residual ', p%resid, ', max |d| ', p%dmax, ' m, omega ', p%omega, ', Newton its ', abs(p%sh%newton_it)
       if (.not. p%converged .and. .not. last) then
          remesh = .true.
          dmax = maxval(abs(p%drel(1:3,:) - p%dasm(1:3,:)))
          if (dmax >= ovst_reasm_tol) then
             fsi_reasm_request = .true.;  p%dasm = p%drel
          end if
       end if
       all_done = all_done .and. p%converged
       deallocate(dtil, rcur)
       call sync_part(p)
    end do
    if (remesh) then
       ! the next sub-iteration repeats the time step on the moved grid from the level-n STATE
       ! (24b): the face fluxes (saved at the step start, re-projected from the level-n grid by
       ! after_motion - the old metrics are kept through the step; rigid moving blocks are rotated
       ! again), the Cartesian field and, with Prs_form = 2, the pressure (mode 1; the full pressure
       ! of Prs_form = 1 is kept as the Poisson initial guess, restore_level_n).  The fluid step is
       ! then the same function of the interface displacement in every sub-iteration, which the
       ! quasi-Newton needs.  The first version reset only the Cartesian field and re-projected the
       ! latest iterate's fluxes incrementally: the fluid solution drifted between sub-iterations
       ! and the residual stalled at 1e-2 of the step increment.
       if (fsi_restore_mode >= 1) then
          call restore_level_n
          call refresh_uvw(1)
       end if
       call sixdof_remesh
    end if
    c = m;  ierr = 0;  i = nint(t)
  end subroutine fsi_iterate

  ! IQN-ILS update of the relaxed interface displacement
  subroutine qn_update(p, rcur, dtil)
    type(part_t), intent(inout), target :: p
    real(dp), intent(in) :: rcur(:,:), dtil(:,:)
    integer :: n, k, nc, i, j, c
    real(dp), allocatable :: V(:,:), W(:,:), Q(:,:), R(:,:), alpha(:), rk(:), qr(:), sc(:)
    real(dp) :: nrm, vmax, lrot
    n = 6*p%sh%nn
    ! the history is kept in a SCALED space: the rotation entries [rad] are multiplied by the wall
    ! offset (half the shell thickness) so that they are commensurate with the translations [m] -
    ! the wall position the fluid sees is u + theta x r0 (24b; without the scaling the least squares
    ! was dominated by the rotations and the secant model of the fluid response was poor)
    allocate(sc(n))
    do i = 1, p%sh%nn
       lrot = 0.5_dp*p%sh%tmax
       if (allocated(p%thick)) lrot = max(0.5_dp*p%thick(i), 1.0e-6_dp)
       sc(6*(i-1)+1:6*(i-1)+3) = 1.0_dp;  sc(6*(i-1)+4:6*i) = lrot
    end do
    ! shift the history (newest first), store the current pair
    k = min(p%iter, MAXQN)
    do i = min(k, MAXQN) - 1, 1, -1
       if (i + 1 <= MAXQN) then
          p%rhist(:, i+1) = p%rhist(:, i);  p%dhist(:, i+1) = p%dhist(:, i)
       end if
    end do
    p%rhist(:, 1) = reshape(rcur, [n])*sc;  p%dhist(:, 1) = reshape(dtil, [n])*sc
    if (p%iter == 1) then
       p%omega = p%def%relax
       p%drel = p%drel + p%omega*rcur
       deallocate(sc)
       return
    end if
    nc = k - 1
    allocate(V(n, nc), W(n, nc), Q(n, nc), R(nc, nc), alpha(nc), rk(n), qr(nc))
    rk = p%rhist(:, 1)
    do i = 1, nc
       V(:, i) = p%rhist(:, i+1) - rk
       W(:, i) = p%dhist(:, i+1) - p%dhist(:, 1)
    end do
    ! modified Gram-Schmidt QR with the filtering of dependent columns
    vmax = maxval(sqrt(sum(V*V, dim=1)))
    R = 0.0_dp;  Q = 0.0_dp;  c = 0
    do i = 1, nc
       Q(:, i) = V(:, i)
       do j = 1, i - 1
          R(j, i) = sum(Q(:, j)*Q(:, i))
          Q(:, i) = Q(:, i) - R(j, i)*Q(:, j)
       end do
       nrm = sqrt(sum(Q(:, i)**2))
       if (nrm <= 1.0e-10_dp*max(vmax, 1.0e-300_dp)) then
          Q(:, i) = 0.0_dp;  R(i, i) = 0.0_dp
       else
          R(i, i) = nrm;  Q(:, i) = Q(:, i)/nrm;  c = c + 1
       end if
    end do
    if (c == 0) then
       p%drel = p%drel + p%omega*rcur
    else
       ! least squares  min || V alpha + r_k ||  ->  R alpha = -Q^T r_k  (back substitution, zero columns skipped)
       do i = 1, nc
          qr(i) = -sum(Q(:, i)*rk)
       end do
       alpha = 0.0_dp
       do i = nc, 1, -1
          if (R(i, i) == 0.0_dp) cycle
          nrm = qr(i)
          do j = i + 1, nc
             nrm = nrm - R(i, j)*alpha(j)
          end do
          alpha(i) = nrm/R(i, i)
       end do
       rk = (p%dhist(:, 1) + matmul(W, alpha))/sc       ! d~_k + W alpha (unscaled)
       p%drel = reshape(rk, [6, p%sh%nn])
       p%omega = 1.0_dp
    end if
    deallocate(V, W, Q, R, alpha, rk, qr, sc)
  end subroutine qn_update

  ! identical structural state on every rank (the solve is replicated; rank 0's values are taken)
  subroutine sync_part(p)
    type(part_t), intent(inout), target :: p
    integer :: ierr, n, i
    real(dp) :: sc(3)
    n = 6*p%sh%nn
    call MPI_Bcast(p%drel, n, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(p%rprev, n, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(p%sh%u1, n, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(p%sh%v1, n, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(p%sh%a1, n, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(p%sh%dth, 3*p%sh%nn, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    sc = [p%omega, p%resid, p%dmax]
    call MPI_Bcast(sc, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    p%omega = sc(1);  p%resid = sc(2);  p%dmax = sc(3)
    call MPI_Bcast(p%converged, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    do i = 1, p%sh%nn
       p%sh%R1(:,:,i) = matmul(rot_exp(p%sh%dth(:,i)), p%sh%R(:,:,i))
    end do
    if (.not. p%converged) call wall_from_structure(p)
  end subroutine sync_part

  !=====================================================================
  ! End of the time step: commit the structural states, statistics,
  ! stresses, previous coordinates of the morphing blocks, exports
  !=====================================================================
  subroutine fsi_finish_step(istep, t)
    integer,  intent(in) :: istep
    real(dp), intent(in) :: t
    integer :: ip, b, i, j, k
    type(part_t), pointer :: p
    if (.not. inited) return
    do ip = 1, npart
       p => pt(ip)
       if (p%active) then
          if (.not. p%converged) p%nfail = p%nfail + 1
          p%vprev = p%vd
          p%vd(1:3, :) = (p%sh%u1(1:3, :) - p%sh%u(1:3, :))/max(dt_dim, 1.0e-300_dp)
          p%vd(4:6, :) = p%sh%dth/max(dt_dim, 1.0e-300_dp)
          call shell_accept(p%sh, p%fext)
          call shell_stress(p%sh, p%vm)
          p%nsteps = p%nsteps + 1;  p%nit_total = p%nit_total + p%iter
       end if
    end do
    has_prev_vel = .true.
    call frames_commit                                  ! prescribed deformation (H1): previous coordinates
    do b = 1, nblk
       if (.not. mb(b)%on) cycle
       do k = 1, blk(b)%nk
          do j = 1, blk(b)%nj
             do i = 1, blk(b)%ni
                mb(b)%xprev(:, i, j, k) = [blk(b)%xn(i,j,k), blk(b)%yn(i,j,k), blk(b)%zn(i,j,k)]
             end do
          end do
       end do
    end do
    call fsi_export(istep, t)
    call dump_grid(istep)
  end subroutine fsi_finish_step

  !---------------------------------------------------------------------
  ! wall node forces of every sheet (grid frame, [N]); all ranks fill
  ! their faces, the sum is replicated
  !---------------------------------------------------------------------
  subroutine wall_forces
    integer :: ip, is, ib, d, b, n, ax, sd, t1, t2, rr(2,3), c1s, c1e, c2s, c2e, m1, m2, i1, i2, ierr
    real(dp), allocatable :: wgt(:,:)
    real(dp) :: Av(3), amag, nrm(3), xf(3), pw, tauw(3), wt, f(3)
    logical :: ok
    type(sheet_t), pointer :: s
    do ip = 1, npart
       do is = 1, pt(ip)%nsheet
          pt(ip)%sheet(is)%fw = 0.0_dp
       end do
    end do
    do ib = 1, nown
       d = own_list(ib);  b = dom(d)%b
       if (.not. fsi_blk(b)) cycle
       call set_block(d)
       do ip = 1, npart
          do is = 1, pt(ip)%nsheet
             s => pt(ip)%sheet(is)
             if (s%b /= b) cycle
             n = s%n
             if (.not. wall_patch_range(d, b, n, ax, sd, t1, t2, rr, c1s, c1e, c2s, c2e)) cycle
             allocate(wgt(s%n1 - 1, s%n2 - 1))
             call ovs_wall_weights(b, n, wgt)
             do m2 = c2s, c2e
                do m1 = c1s, c1e
                   i1 = m1 - rr(1,t1) + 1;  i2 = m2 - rr(1,t2) + 1
                   wt = wgt(i1, i2)
                   if (wt <= 0.0_dp) cycle
                   call wall_face_stress(b, n, ax, sd, t1, t2, m1, m2, ok, Av, amag, nrm, xf, pw, tauw)
                   if (.not. ok) cycle
                   f = 0.25_dp*wt*(-pw*Av + tauw*amag)
                   s%fw(:, i1, i2) = s%fw(:, i1, i2) + f;      s%fw(:, i1+1, i2) = s%fw(:, i1+1, i2) + f
                   s%fw(:, i1+1, i2+1) = s%fw(:, i1+1, i2+1) + f;  s%fw(:, i1, i2+1) = s%fw(:, i1, i2+1) + f
                end do
             end do
             deallocate(wgt)
          end do
       end do
    end do
    do ip = 1, npart
       do is = 1, pt(ip)%nsheet
          s => pt(ip)%sheet(is)
          call MPI_Allreduce(MPI_IN_PLACE, s%fw, 3*s%n1*s%n2, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
       end do
    end do
  end subroutine wall_forces

  !---------------------------------------------------------------------
  ! structural load vector of a part (part frame): fluid wall forces,
  ! gravity, inertial loads of the moving frame
  !---------------------------------------------------------------------
  subroutine build_loads(p)
    type(part_t), intent(inout), target :: p
    integer :: is, i1, i2, a, nd, i
    real(dp) :: Rg(3,3), O0(3), dg(3), vG(3), aG(3), omG(3), alG(3), fs(3), r(3), Ra(3,3), gB(3), rpos(3), vrel(3), af(3), w
    type(sheet_t), pointer :: s
    call sixdof_frame_kin(p%k, Rg, O0, dg, vG, aG, omG, alG)
    p%fext = 0.0_dp
    do is = 1, p%nsheet
       s => p%sheet(is)
       do i2 = 1, s%n2
          do i1 = 1, s%n1
             fs = matmul(transpose(Rg), s%fw(:, i1, i2))
             do a = 1, 4
                nd = s%att(a, i1, i2)
                if (nd == 0) cycle
                w = s%wgt(a, i1, i2)
                Ra = matmul(rot_exp(p%sh%dth(:, nd)), p%sh%R(:,:,nd))
                r = matmul(Ra, s%off(:, a, i1, i2))
                p%fext(1:3, nd) = p%fext(1:3, nd) + w*fs
                p%fext(4:6, nd) = p%fext(4:6, nd) + w*cross3(r, fs)
             end do
          end do
       end do
    end do
    ! gravity and frame inertia (translational dofs)
    gB = matmul(transpose(Rg), grav)
    do i = 1, p%sh%nn
       p%fext(1:3, i) = p%fext(1:3, i) + p%sh%mass(1, i)*gB
       if (any(abs(aG) > 0.0_dp) .or. any(abs(omG) > 0.0_dp) .or. any(abs(alG) > 0.0_dp)) then
          rpos = matmul(Rg, p%sh%X0(:, i) + p%sh%u1(1:3, i) - O0)
          vrel = matmul(Rg, p%sh%v1(1:3, i))
          af = aG + cross3(alG, rpos) + cross3(omG, cross3(omG, rpos)) + 2.0_dp*cross3(omG, vrel)
          p%fext(1:3, i) = p%fext(1:3, i) - p%sh%mass(1, i)*matmul(transpose(Rg), af)
       end if
    end do
    dg = 0.0_dp;  vG = 0.0_dp
  end subroutine build_loads

  !=====================================================================
  ! displacement [m] (grid frame) and von Mises stress [Pa] of wall node
  ! (i, j, k) of block b (0 when the node is not on a flexible sheet)
  !=====================================================================
  logical function fsi_wall_node_value(b, i, j, k, d, vm) result(is)
    integer,  intent(in)  :: b, i, j, k
    real(dp), intent(out) :: d(3), vm
    integer :: ip, isheet, idx(3), i1, i2, a, nd
    real(dp) :: Rg(3,3), O0(3), dg(3), vG(3), aG(3), omG(3), alG(3)
    type(sheet_t), pointer :: s
    d = 0.0_dp;  vm = 0.0_dp;  is = .false.
    if (.not. inited) return
    if (.not. fsi_blk(b)) return
    do ip = 1, npart
       do isheet = 1, pt(ip)%nsheet
          s => pt(ip)%sheet(isheet)
          if (s%b /= b) cycle
          idx = [i, j, k]
          if (idx(s%ax) /= merge(1, ncell(b, s%ax) + 1, s%sd == 1)) cycle
          if (idx(s%t1) < s%rr(1,s%t1) .or. idx(s%t1) > s%rr(2,s%t1)) cycle
          if (idx(s%t2) < s%rr(1,s%t2) .or. idx(s%t2) > s%rr(2,s%t2)) cycle
          i1 = idx(s%t1) - s%rr(1,s%t1) + 1;  i2 = idx(s%t2) - s%rr(1,s%t2) + 1
          call sixdof_frame_kin(pt(ip)%k, Rg, O0, dg, vG, aG, omG, alG)
          d = matmul(Rg, s%dw(:, i1, i2))
          vm = 0.0_dp
          do a = 1, 4
             nd = s%att(a, i1, i2)
             if (nd > 0) vm = vm + s%wgt(a, i1, i2)*pt(ip)%vm(nd)
          end do
          is = .true.;  return
       end do
    end do
  end function fsi_wall_node_value

  !=====================================================================
  ! exports: monitor node files (Data_Intvl), log line
  !=====================================================================
  subroutine open_monitors(p)
    type(part_t), intent(inout), target :: p
    integer :: q, iu
    character(len=256) :: fname
    logical :: ex
    inquire(file=out_dir//'.', exist=ex)
    if (.not. ex) call execute_command_line('mkdir -p '//out_dir, wait=.true.)
    do q = 1, p%def%nmon
       if (p%mon_sheet(q) == 0) cycle
       fname = out_dir//trim(casename)//'_fsi_'//sixdof_part_stem(p%k)//'_node'//itoa(q)//'.dat'
       open(newunit=iu, file=trim(fname), status=merge('unknown', 'replace', run_mode == 1), position='append', action='write')
       p%iu_mon(q) = iu
       if (run_mode /= 1) then
          write(iu,'(A,4I5,A)') '# probCFD FSI monitor node (update 24): '//sixdof_part_name(p%k)//', wall node block i j k =', &
             p%def%mon(:,q), '  (every Data_Intvl steps)'
          write(iu,'(A)') '#   displacement d of the wall node from its rigid-body position, in the PART frame (grid axes at t = 0) [m],'
          write(iu,'(A)') '#   velocity of the deformation [m/s], von Mises stress of the shell at the node [Pa], FSI sub-iterations,'
          write(iu,'(A)') '#   relative interface residual at the end of the step'
          write(iu,'(A)') '#   step      time[s]            dx              dy              dz            |d|             vx'// &
             '              vy              vz        vonMises[Pa]  iters   residual'
       end if
    end do
  end subroutine open_monitors

  subroutine fsi_export(istep, t)
    integer,  intent(in) :: istep
    real(dp), intent(in) :: t
    integer :: ip, q, is, i1, i2, a, nd
    real(dp) :: d(3), v(3), vm, w
    type(part_t), pointer :: p
    type(sheet_t), pointer :: s
    if (.not. inited .or. myrank /= 0) return
    if (mod(istep, data_interval) /= 0 .and. istep /= 0) return
    do ip = 1, npart
       p => pt(ip)
       do q = 1, p%def%nmon
          if (p%iu_mon(q) == 0) cycle
          is = p%mon_sheet(q);  i1 = p%mon_i1(q);  i2 = p%mon_i2(q)
          s => p%sheet(is)
          d = s%dw(:, i1, i2)
          v = 0.0_dp;  vm = 0.0_dp
          do a = 1, 4
             nd = s%att(a, i1, i2)
             if (nd == 0) cycle
             w = s%wgt(a, i1, i2)
             v = v + w*(p%sh%v(1:3, nd) + cross3(p%sh%v(4:6, nd), matmul(p%sh%R(:,:,nd), s%off(:, a, i1, i2))))
             vm = vm + w*p%vm(nd)
          end do
          write(p%iu_mon(q),'(I7,ES13.5,8ES16.8,I6,ES11.3)') istep, t, d, sqrt(sum(d*d)), v, vm, p%iter, p%resid
          flush(p%iu_mon(q))
       end do
    end do
  end subroutine fsi_export

  ! screen monitor (rank 0, log interval)
  subroutine fsi_monitor
    integer :: ip
    real(dp) :: ek, es
    type(part_t), pointer :: p
    if (.not. inited .or. myrank /= 0) return
    do ip = 1, npart
       p => pt(ip)
       write(*,'(A)') ' '
       if (frames_part(p%k) > 0) then
          write(*,'(A,ES12.5,A)') ' [ FSI ] '//sixdof_part_name(p%k)//': prescribed deformation (frame file), grid at t = ', t_now, ' s'
          cycle
       end if
       if (.not. p%active) then
          write(*,'(A,ES10.3,A)') ' [ FSI ] '//sixdof_part_name(p%k)//': structure held until FSI_Start_T = ', p%def%start_t, ' s'
          cycle
       end if
       call shell_energy(p%sh, ek, es)
       write(*,'(A,I3,A,ES9.2,A,F5.2,A,ES10.3,A,ES10.3,A)') ' [ FSI ] '//sixdof_part_name(p%k)//': ', p%iter, ' sub-iterations, residual ', &
          p%resid, ' (omega ', p%omega, '), max |d| ', p%dmax, ' m, max von Mises ', maxval(p%vm), ' Pa'
       write(*,'(A,I3,A,ES10.3,A,ES10.3,A)') '         Newton its of the last solve ', abs(p%sh%newton_it), ', kinetic / strain energy ', &
          ek, ' / ', es, ' J'
       if (.not. p%converged) write(*,'(A,I0,A)') '         NOTE: not converged within FSI_Max_Iter = ', p%def%max_iter, ' sub-iterations'
    end do
  end subroutine fsi_monitor

  subroutine fsi_summary
    integer :: ip
    type(part_t), pointer :: p
    if (.not. inited .or. myrank /= 0) return
    write(*,'(A)') ' '
    do ip = 1, npart
       p => pt(ip)
       if (p%nsteps == 0) cycle
       write(*,'(A,I0,A,F6.2,A,F6.2,A,I0,A)') '  FSI '//sixdof_part_name(p%k)//': ', p%nsteps, ' coupled steps, ', &
          real(p%nit_total)/p%nsteps, ' sub-iterations and ', real(p%nnewton)/max(p%nit_total,1), ' Newton its per structural solve on average, ', &
          p%nfail, ' steps not converged within FSI_Max_Iter'
    end do
  end subroutine fsi_summary

  subroutine fsi_finalize
    integer :: ip, q
    if (.not. inited) return
    if (myrank == 0) then
       do ip = 1, npart
          do q = 1, pt(ip)%def%nmon
             if (pt(ip)%iu_mon(q) /= 0) close(pt(ip)%iu_mon(q))
          end do
       end do
    end if
  end subroutine fsi_finalize

  !---------------------------------------------------------------------
  ! the mid-surface mesh of a part (PLOT3D, one block per sheet: nup x nsp x 1)
  !---------------------------------------------------------------------
  subroutine write_shell_mesh(p)
    type(part_t), intent(in), target :: p
    integer :: iu, is, m, sidx, l, e, a, nact, act(64)
    real(dp), allocatable :: mat(:), clamp(:)
    character(len=256) :: fname
    logical :: ex
    inquire(file=out_dir//'.', exist=ex)
    if (.not. ex) call execute_command_line('mkdir -p '//out_dir, wait=.true.)
    allocate(mat(p%sh%nn), clamp(p%sh%nn))
    mat = 0.0_dp
    do e = 1, p%sh%ne
       do a = 1, 4
          mat(p%sh%conn(a, e)) = max(mat(p%sh%conn(a, e)), real(p%matid(e), dp))
       end do
    end do
    do m = 1, p%sh%nn
       clamp(m) = merge(1.0_dp, 0.0_dp, p%sh%fixed(1, m))
    end do
    nact = 0
    do is = 1, p%nsheet
       if (p%sheet(is)%passive) cycle
       nact = nact + 1;  act(nact) = is
    end do
    fname = out_dir//trim(casename)//'_fsi_'//sixdof_part_stem(p%k)//'_shell.grd'
    open(newunit=iu, file=trim(fname), form='unformatted', status='replace', action='write')
    write(iu) nact
    write(iu) (p%sheet(act(is))%nup, p%sheet(act(is))%nsp, 1, is=1,nact)
    do is = 1, p%nsheet
       if (p%sheet(is)%passive) cycle
       write(iu) (((real(p%sh%X0(l, p%sheet(is)%snode(m, sidx)), 4), m=1,p%sheet(is)%nup), sidx=1,p%sheet(is)%nsp), l=1,3)
    end do
    close(iu)
    fname = out_dir//trim(casename)//'_fsi_'//sixdof_part_stem(p%k)//'_shell.p3d'
    open(newunit=iu, file=trim(fname), form='unformatted', status='replace', action='write')
    write(iu) nact
    write(iu) (p%sheet(act(is))%nup, p%sheet(act(is))%nsp, 1, 3, is=1,nact)
    do is = 1, p%nsheet
       if (p%sheet(is)%passive) cycle
       write(iu) ((real(p%thick(p%sheet(is)%snode(m, sidx)), 4), m=1,p%sheet(is)%nup), sidx=1,p%sheet(is)%nsp), &
                 ((real(mat(p%sheet(is)%snode(m, sidx)), 4), m=1,p%sheet(is)%nup), sidx=1,p%sheet(is)%nsp), &
                 ((real(clamp(p%sheet(is)%snode(m, sidx)), 4), m=1,p%sheet(is)%nup), sidx=1,p%sheet(is)%nsp)
    end do
    close(iu)
    open(newunit=iu, file=out_dir//trim(casename)//'_fsi_'//sixdof_part_stem(p%k)//'_shell.nam', status='replace', action='write')
    write(iu,'(A)') 'Thickness [m]';  write(iu,'(A)') 'Submat_id [-]';  write(iu,'(A)') 'Clamped [-]'
    close(iu)
    write(*,'(A)') '  FSI: wrote the mid-surface mesh of '//sixdof_part_name(p%k)//' -> '//trim(fname)//' (+ .grd, .nam)'
    deallocate(mat, clamp)
  end subroutine write_shell_mesh

  !=====================================================================
  ! restart: structural states (rank 0 writes; all ranks receive)
  !=====================================================================
  subroutine fsi_restart_write(iu)
    integer, intent(in) :: iu
    integer :: ip, i
    if (myrank /= 0) return
    if (.not. inited) then
       write(iu) 0;  return
    end if
    write(iu) npart
    do ip = 1, npart
       if (frames_part(pt(ip)%k) > 0) then              ! prescribed deformation (H1): the grid is a function of time
          write(iu) 0, merge(1, 0, has_prev_vel);  cycle
       end if
       write(iu) pt(ip)%sh%nn, merge(1, 0, has_prev_vel)
       write(iu) pt(ip)%sh%u;  write(iu) pt(ip)%sh%v;  write(iu) pt(ip)%sh%a;  write(iu) pt(ip)%sh%fext_n
       write(iu) pt(ip)%vprev
       write(iu) (pt(ip)%sh%R(:,:,i), i=1,pt(ip)%sh%nn)
    end do
    ! 24b: the relaxed interface displacement the grid was morphed with at the end of the step (it
    ! differs from the structural output by the unconverged interface residual): the restored grid
    ! is morphed with it, so that a restart continues from the same grid
    write(iu) 2401
    do ip = 1, npart
       if (frames_part(pt(ip)%k) > 0) cycle
       write(iu) pt(ip)%drel
    end do
    ! 24c: the discrete velocity of the last step (predictor)
    write(iu) 2402
    do ip = 1, npart
       if (frames_part(pt(ip)%k) > 0) cycle
       write(iu) pt(ip)%vd
    end do
  end subroutine fsi_restart_write

  subroutine fsi_restart_read(iu, ok)
    integer, intent(in) :: iu
    logical, intent(out) :: ok
    integer :: np, ip, nn, ios, ierr, i, hp
    ok = .false.
    np = 0
    if (myrank == 0) then
       read(iu, iostat=ios) np
       if (ios /= 0) np = 0
       if (np /= npart) then
          if (np > 0 .or. npart > 0) write(*,'(A,I0,A,I0,A)') '  NOTE: the restart file holds ', np, &
             ' flexible parts, this run has ', npart, ' - the structures start undeformed'
          np = -1
       end if
    end if
    call MPI_Bcast(np, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if (np <= 0) return
    do ip = 1, npart
       if (myrank == 0) then
          read(iu) nn, hp
          if (nn /= pt(ip)%sh%nn) then
             write(*,'(A)') '  NOTE: the mid-surface of '//sixdof_part_name(pt(ip)%k)//' differs from the restart file - undeformed start'
             nn = -1
          end if
       end if
       call MPI_Bcast(nn, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
       call MPI_Bcast(hp, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
       if (nn < 0) return
       if (frames_part(pt(ip)%k) > 0) then              ! prescribed deformation (H1): nothing stored
          has_prev_vel = (hp == 1);  cycle
       end if
       if (myrank == 0) then
          read(iu) pt(ip)%sh%u;  read(iu) pt(ip)%sh%v;  read(iu) pt(ip)%sh%a;  read(iu) pt(ip)%sh%fext_n
          read(iu) pt(ip)%vprev
          read(iu) (pt(ip)%sh%R(:,:,i), i=1,nn)
       end if
       call MPI_Bcast(pt(ip)%sh%u, 6*nn, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
       call MPI_Bcast(pt(ip)%sh%v, 6*nn, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
       call MPI_Bcast(pt(ip)%sh%a, 6*nn, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
       call MPI_Bcast(pt(ip)%sh%fext_n, 6*nn, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
       call MPI_Bcast(pt(ip)%vprev, 6*nn, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
       call MPI_Bcast(pt(ip)%sh%R, 9*nn, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
       has_prev_vel = (hp == 1)
       call plane_project(pt(ip))
       call plane_average_state(pt(ip))
       call shell_set_state_from_trial(pt(ip)%sh)
       call shell_stress(pt(ip)%sh, pt(ip)%vm)
       pt(ip)%drel(1:3, :) = pt(ip)%sh%u(1:3, :);  pt(ip)%drel(4:6, :) = 0.0_dp
    end do
    ! 24b: the interface displacement of the grid (older files end before it: the structural output is used)
    hp = 0
    if (myrank == 0) then
       read(iu, iostat=ios) hp
       if (ios /= 0 .or. hp /= 2401) hp = 0
    end if
    call MPI_Bcast(hp, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    do ip = 1, npart
       if (frames_part(pt(ip)%k) > 0) cycle
       if (hp == 2401) then
          if (myrank == 0) read(iu) pt(ip)%drel
          call MPI_Bcast(pt(ip)%drel, 6*pt(ip)%sh%nn, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
       end if
       call plane_average_drel(pt(ip))
       pt(ip)%dasm = pt(ip)%drel
       call wall_from_structure(pt(ip))
    end do
    ! 24c: the discrete velocity of the last step (older files: the structural velocity)
    hp = 0
    if (myrank == 0) then
       read(iu, iostat=ios) hp
       if (ios /= 0 .or. hp /= 2402) hp = 0
    end if
    call MPI_Bcast(hp, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    do ip = 1, npart
       if (frames_part(pt(ip)%k) > 0) cycle
       if (hp == 2402) then
          if (myrank == 0) read(iu) pt(ip)%vd
          call MPI_Bcast(pt(ip)%vd, 6*pt(ip)%sh%nn, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
       else
          pt(ip)%vd = pt(ip)%sh%v
       end if
    end do
    ok = .true.
    if (myrank == 0) write(*,'(A,I0,A)') '  FSI: structural states of ', npart, ' flexible part(s) restored'
  end subroutine fsi_restart_read

end module fsi
