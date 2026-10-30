!=======================================================================
! sixdof.F90
!   6-DoF RIGID-BODY MOTION of the overset objects (update 12) with
!   CHILD and SUB-CHILD bodies (update 14: three-level hierarchy
!   body > child body > sub-child body, e.g. hull > stern drive >
!   propeller).
!
!   CONVENTIONS (SNAME 1950 / Fossen 1994)
!   --------------------------------------
!   * Grid frame G (the .grd coordinates): x downstream (= aft when the
!     hull faces the inflow), y starboard, z up, gravity (0,0,-g).
!   * SNAME earth-fixed frame E: x forward, y starboard, z DOWN, with the
!     same origin as the grid:   p_E = F p_G,  F = diag(-1, 1, -1)
!     (a proper rotation, det F = +1, so vectors and moments map alike).
!   * Body-fixed frame B_k of every rigid part k (body, child or sub-
!     child): aligned with E at t = 0 (the initial attitude of the grid
!     is the zero attitude), origin at its centre of rotation O_k
!     (Rot_Center of the namelist), so the rigid-body equations are
!     Fossen's general form (2.89)/(2.91) with the centre of gravity at
!     r_G = CG - O (parallel-axes theorem 2.108).  Inertia, Prod_Inertia
!     of a child / sub-child are given in ITS OWN axes (B_k).
!   * BODY (level 1): generalized position eta = [x y z phi theta psi]
!     = displacement of O from its initial position in E [m] and the
!     Euler angles (yaw-pitch-roll, Z-Y-X order, R = Rz(psi) Ry(theta)
!     Rx(phi)) [rad].  Body velocities nu = [u v w p q r].
!   * CHILD (level 2, mother = a body) and SUB-CHILD (level 3, mother =
!     a child): eta is the RELATIVE motion with respect to the mother:
!     eta(1:3) = displacement of O_k from its initial position measured
!     in the MOTHER frame B_m, eta(4:6) = Euler angles of B_k relative
!     to B_m (R_k = R_m R_rel); nu = relative body velocities (velocity
!     of O_k relative to the mother in B_k, angular velocity of k
!     relative to the mother in B_k).  The absolute pose follows by
!     composition (update_transform):
!        R_k = R_m R_rel,   p_k = p_m + R_m (r0 + eta(1:3)),
!        v_k = v_m + R_m (omega_m x (r0 + eta(1:3)) + etadot(1:3)),
!        Omega_k = Omega_m + R_k omega_rel,
!     r0 = initial offset O_k - O_m in B_m, omega_rel = J2^-1 etadot(4:6).
!     A child whose six DOFs are 0 (fixed) moves rigidly with its mother.
!   * Positive directions: surge forward, sway to starboard, heave DOWN,
!     roll starboard-down, pitch bow-up, yaw bow-to-starboard.
!
!   DEGREES OF FREEDOM (per part, Motion_DOF, SNAME order)
!     0 fixed   : the component of eta is held (child: rigid with the
!                 mother in that DOF)
!     1 imposed : from the file  <case>_body<k>_imposed.dat
!                 (child)        <case>_body<k>_child<c>_imposed.dat
!                 (sub-child)    <case>_body<k>_child<c>_sub<s>_imposed.dat
!                 columns t, x, y, z [m], phi, theta, psi [deg] (relative
!                 to the mother for a child / sub-child); cubic Hermite
!                 interpolation in time; the values outside the table
!                 are held at the end values.  Any rotation / translation
!                 / combined motion can be imposed through the Euler
!                 angles (e.g. a propeller: phi = n t about the shaft
!                 axis = roll axis of the sub-child frame)
!     2 free    : integrated from the fluid forces + weight (from
!                 Motion_Start_T on; before that the DOF is held)
!
!   DYNAMICS of a part k (Newton-Euler about O_k in B_k, written for the
!   RELATIVE acceleration, Fossen's M_RB with the COMPOSITE mass, CG and
!   inertia of k and all its descendants at their current relative
!   poses - the parts below k are carried along rigidly):
!     M_RB x = b,  x = [R_rel^T a_rel ; omegadot_rel]
!     b(1:3) = F - M (a_fr + alpha_fr x r_G + Omega x (Omega x r_G))
!     b(4:6) = tau - Omega x (I_O Omega) - I_O alpha_fr - M r_G x a_fr
!   F, tau: the fluid wrench of k's subtree about O_k + weight + the
!   artificial damping; Omega: absolute angular velocity in B_k; a_fr,
!   alpha_fr: acceleration of the mother frame at O_k and its angular
!   acceleration (centripetal, Coriolis and Euler terms, zero for a
!   body).  Written in the eta representation M_eta etadd = b_eta so
!   that fixed / imposed DOFs are simple rows; the free block is solved
!   directly.  Loose coupling: the forces of the converged previous time
!   step drive the motion of the current step (explicit); an added-mass
!   estimate Added_Mass (per DOF, body frame) is added to both sides,
!     (M_eta + M_A) a^{n+1} = b + M_A a^n,
!   which is consistent at convergence and stabilizes the explicit
!   fluid-body coupling of light bodies.  Time integration: velocity
!   Verlet on eta.  The parts are advanced mothers first (bodies, then
!   children, then sub-children) so that a child sees the mother's
!   motion of the current step.
!
!   MOTION COUPLING child <-> mother: the fluid wrench of a child (and of
!   its own sub-children) is transferred to the mother at O_m (force +
!   moment + lever arm) EXCEPT the components along the child's FREE
!   DOFs (a free joint transmits nothing along its axis: a free-wheeling
!   propeller passes its thrust but not its torque; fixed and imposed
!   DOFs transmit everything - an imposed propeller rotation is a
!   motor-driven shaft, its torque is reacted by the mother).  The
!   mother's inertia is the composite inertia of its subtree at the
!   current relative poses.  Neglected: the inertial reaction of the
!   RELATIVE motion of a child on its mother (I_child alpha_rel and the
!   gyroscopic terms), valid for appendages light compared with the
!   mother (propeller, rudder, stern drive on a hull).
!
!   GRID MOTION: every block of a part (rows 6-8 of the csv) is moved
!   rigidly from its ORIGINAL coordinates,  x = Rg (x0 - O0) + O0 + d,
!   Rg = F R_k F, d = F p_k - O0 (absolute pose of the part); the half-
!   index coordinates rotate with it and the metric vectors m(:,q)
!   rotate (J is invariant).  The grid velocity of every face is stored
!   in Ug1..3 (units of U1..3) from the absolute velocity of O_k and the
!   absolute angular velocity; the advecting fluxes become U - Ug (ALE);
!   the no-slip wall velocity equals the grid velocity (boundary.F90).
!   After a displacement the overset assembly is redone from the new
!   coordinates (overset.F90, ovst_reassemble); every part is its own
!   rigid motion group of the donor cache.
!
!   FORCES: F = sum over the wall faces of  -p n A + tau_w A, n = normal
!   into the fluid; p of the first cell (zero-gradient wall), tau_w =
!   mu_eff (u_t - u_wall,t)/d with the first-cell tangential velocity
!   relative to the wall, mu_eff including the wall-function factor;
!   each face weighted by the fraction of its triangles kept by the
!   surface-overset trimming and skipped when its cell is a hole.
!   Moments about the CURRENT centre of rotation O_k of the part owning
!   the wall block.
!
!   FILES (every Data_Intvl steps, rank 0) in the folder Output_data/:
!     body k (own walls only):
!       <case>_force_body<k>.dat     X Y Z [N] in B_k: friction, pressure, total
!       <case>_moment_body<k>.dat   K M N [N m] about O_k in B_k
!       <case>_force_earth<k>.dat   in E (SNAME earth: x fwd, z down)
!       <case>_moment_earth<k>.dat
!       <case>_motion<k>.dat        eta [m, deg], nu [m/s, deg/s],
!                                   eta_dot, nu_dot, O and CG in G [m]
!     child c of body k (own walls; force / moment / motion RELATIVE to
!     the mother body, i.e. in the mother frame B_k, moment about O_c):
!       <case>_body<k>_child<c>_force.dat, _moment.dat, _motion.dat
!     sub-child s of child c of body k (relative to the mother child):
!       <case>_body<k>_child<c>_sub<s>_force.dat, _moment.dat, _motion.dat
!=======================================================================
module sixdof
  use mpi
  use pid_control, only: pidc, npid, pid_nml, pid_check, pid_find, pid_step, pid_open, pid_export, pid_finalize, &
                         pid_axis_unit, pid_summary, pid_var_name
  use global_vars
  use mpi_vars
  use curvi, only: cart_to_contra, contra_to_cart
  use b2b_runtime, only: refresh_all, refresh_uvw
  use overset, only: ovst_reassemble, ovs_wall_weights, ovs_update_faces, ovs_update_cells, ovs_finish_status, ovs_counts
  use stop_criteria, only: stop_clean
  implicit none
  private
  public :: sixdof_nml, sixdof_precheck, sixdof_check_movement, sixdof_init, sixdof_advance, sixdof_forces, sixdof_export, &
            sixdof_wall_velocity, sixdof_monitor, sixdof_finalize, sixdof_summary, &
            sixdof_nbody, sixdof_get_state, sixdof_set_state, sixdof_restore_pose, sixdof_on_flag, sixdof_clear_ext, &
            sixdof_node_of, sixdof_node_pose, sixdof_set_ext, sixdof_part_name, &
            wall_patch_range, wall_face_stress, is_wall_patch, tang_axes, sixdof_part_stem, &
            sixdof_fsi_def, sixdof_mass_data, sixdof_frame_kin, sixdof_remesh, &
            Moving_BGR, Moving_BGR_DoF, BGR_Following_Body, Moving_BGR_Filter_T

  integer, parameter :: MAXB = 64
  real(dp), parameter :: pi = 3.14159265358979323846_dp
  real(dp), parameter :: Fm(3,3) = reshape([-1.0_dp, 0.0_dp, 0.0_dp,  0.0_dp, 1.0_dp, 0.0_dp,  0.0_dp, 0.0_dp, -1.0_dp], [3,3])
  character(len=*), parameter :: out_dir = 'Output_data/'

  ! flexible part (update 24, SW_FSI = 1 in the part's group; evaluated by fsi.F90)
  integer, parameter, public :: MAXFG = 16
  type, public :: fsidef_t
     integer  :: sw = 0                       ! SW_FSI
     real(dp) :: E = 0.0_dp, nu = 0.3_dp      ! Youngs_modulus [Pa], Poisson_ratio
     real(dp) :: rho = -1.0_dp                ! Density [kg/m3] (-1: from Mass and the enclosed volume)
     integer  :: nnof = 0                     ! Num_NoFSI_Group
     integer  :: nof_blk(MAXFG) = 0, nof_rng(2,3,MAXFG) = 0     ! NoFSI_block_num<n>, NoFSI_i/j/k_range<n> (0 = whole range)
     integer  :: nsub = 0                     ! Num_FSI_Submat_regions
     integer  :: sub_blk(MAXFG) = 0, sub_rng(2,3,MAXFG) = 0     ! Submat_block_num<n>, Submat_i/j/k_range<n>
     real(dp) :: sub_rho(MAXFG) = -1.0_dp, sub_E(MAXFG) = -1.0_dp, sub_nu(MAXFG) = -1.0_dp
     real(dp) :: start_t = 0.0_dp             ! FSI_Start_T [s]
     real(dp) :: zeta = 0.0_dp                ! Damping_ratio
     real(dp) :: tmin_frac = 0.01_dp          ! FSI_Min_Thickness (fraction of the maximum)
     real(dp) :: tol = 1.0e-5_dp              ! FSI_Tol
     real(dp) :: relax = 0.5_dp               ! FSI_Relax (initial Aitken factor)
     integer  :: max_iter = 20                ! FSI_Max_Iter
     integer  :: mode = 1                     ! FSI_Mode: 1 dynamic, 0 quasi-static
     integer  :: nmon = 0, mon(4,MAXFG) = 0   ! FSI_Monitor_node<n>: block i j k
     integer  :: plane = 0                    ! FSI_Plane_Strain (24c): span axis 1/2/3 of a quasi-2-D case, 0 = off
     real(dp) :: test_amp = 0.0_dp, test_freq = 0.0_dp   ! FSI_Test_Motion: imposed bending mode (verification / start-up aid, 24c)
     integer  :: test_mode = 0                ! FSI_Test_Motion 3rd value: 0 oscillation, 1 ramp to a static deflection and hold, 2 pluck (release at 2.25 periods)
     integer  :: test_shp = 0                 ! FSI_Test_Motion 4th value (mode 2): 0 the static tip-load shape, m >= 1 the m-th vibration mode in vacuum
  end type fsidef_t

  type :: body_t
     ! ---- flexible part (update 24) ---------------------------------------
     type(fsidef_t) :: fsi
     ! ---- identity (update 14) -------------------------------------------
     integer  :: lev = 1                    ! 1 body, 2 child body, 3 sub-child body
     integer  :: bid = 0, cid = 0, sid = 0  ! body id, child id (0: none), sub-child id (0: none)
     integer  :: par = 0                    ! mother node (0: a body)
     ! ---- input --------------------------------------------------------
     logical  :: defined = .false.
     integer  :: dof(6) = 0                 ! 0 fixed, 1 imposed, 2 free (surge sway heave roll pitch yaw)
     real(dp) :: mass = 0.0_dp
     real(dp) :: cg0(3) = 0.0_dp            ! centre of gravity, grid frame [m]
     real(dp) :: rot0(3) = 0.0_dp           ! centre of rotation O, grid frame [m]
     logical  :: rot_given = .false.
     real(dp) :: Icg(3) = 0.0_dp            ! Ixx Iyy Izz about the CG, own body axes [kg m2]
     real(dp) :: Ipr(3) = 0.0_dp            ! Ixy Ixz Iyz [kg m2]
     real(dp) :: madd(6) = 0.0_dp           ! added-mass estimate per DOF (stabilization) [kg, kg m2]
     real(dp) :: damp(6) = 0.0_dp           ! artificial damping coefficient per DOF [1/s]
     character(len=128) :: mfile = ' '      ! imposed-motion file
     integer  :: mother_b = 0, mother_c = 0 ! Mother_body / Mother_child_body of the namelist (checked)
     ! ---- derived --------------------------------------------------------
     logical  :: own_moving = .false.       ! any own DOF /= 0
     logical  :: moving = .false.           ! own_moving or a moving mother (the grids move)
     integer  :: nblk = 0                   ! blocks of the part
     real(dp) :: rcg(3) = 0.0_dp            ! own CG - O in B_k [m]
     real(dp) :: Ic(3,3) = 0.0_dp           ! own inertia tensor about the CG in B_k
     real(dp) :: r0(3) = 0.0_dp             ! initial offset O_k - O_mother in B_mother (E-aligned) [m]
     ! composite (k + descendants at the current relative poses)
     real(dp) :: Mc = 0.0_dp, rcgc(3) = 0.0_dp, IOc(3,3) = 0.0_dp, MRBc(6,6) = 0.0_dp
     ! ---- state (bodies: absolute; children: relative to the mother) ------
     real(dp) :: eta(6) = 0.0_dp, etad(6) = 0.0_dp, etadd(6) = 0.0_dp
     real(dp) :: nu(6) = 0.0_dp, nud(6) = 0.0_dp
     ! ---- absolute pose / velocities (composition) ---------------------
     real(dp) :: R(3,3)                     ! absolute rotation B_k -> E
     real(dp) :: pE(3) = 0.0_dp             ! position of O_k in E [m]
     real(dp) :: vE(3) = 0.0_dp             ! absolute velocity of O_k in E [m/s]
     real(dp) :: OmE(3) = 0.0_dp            ! absolute angular velocity in E [rad/s]
     real(dp) :: aE(3) = 0.0_dp             ! absolute acceleration of O_k in E [m/s2]
     real(dp) :: alB(3) = 0.0_dp            ! absolute angular acceleration in B_k [rad/s2]
     real(dp) :: Rg(3,3), dg(3) = 0.0_dp    ! grid frame: x = Rg (x0 - O0) + O0 + dg
     real(dp) :: Rg_prev(3,3)               ! grid rotation before the current step
     real(dp) :: dg_prev(3) = 0.0_dp         ! grid translation before the current step (update 20: discrete grid velocity)
     real(dp) :: vG(3) = 0.0_dp, omG(3) = 0.0_dp   ! velocity of O_k and angular velocity in G [m/s, rad/s]
     ! ---- forces (grid frame, about O_k) [N], [N m] -----------------------
     real(dp) :: fp(3) = 0.0_dp, ff(3) = 0.0_dp, mp(3) = 0.0_dp, mf(3) = 0.0_dp
     real(dp) :: area = 0.0_dp
     real(dp) :: Fs(3) = 0.0_dp, Ms(3) = 0.0_dp   ! subtree fluid wrench in E about O_k (free joint components removed)
     real(dp) :: fext(3) = 0.0_dp, mext(3) = 0.0_dp   ! external wrench of the force models (grid frame, about O_k; update 19)
     logical  :: vhalf(6) = .false.         ! the DOF was advanced as free in the previous step (velocity Verlet half-kick pending)
     ! ---- PID-controlled rotation (update 21): the part turns about the axis pid_axB
     !      (unit, B_mother / E-aligned frame) through O by pid_ang; replaces its own DOFs
     logical  :: pid_rot = .false.
     integer  :: npid = 0, ipid(4) = 0     ! its controllers (pid_control indices; the first rotates it)
     real(dp) :: pid_axB(3) = 0.0_dp, pid_ang = 0.0_dp, pid_angd = 0.0_dp
     ! ---- imposed motion table ------------------------------------------
     integer  :: nt = 0
     real(dp), allocatable :: tt(:), et(:,:)   ! (6, nt) [m], [rad]
     ! ---- output units --------------------------------------------------
     integer  :: iu(5) = -1
  end type body_t

  type(body_t), allocatable, target :: bd(:)     ! the nodes (bodies, children, sub-children)

  ! probe points assigned to a part (update 22, read here, evaluated by probes.F90): grid-frame
  ! position at t = 0, point number within the part, the part's ids
  type, public :: prdef_t
     integer  :: b = 0, c = 0, s = 0, n = 0
     real(dp) :: x(3) = 0.0_dp
  end type prdef_t
  integer, parameter, public :: MAXPRDEF = 256
  type(prdef_t), public :: prdef(MAXPRDEF)
  integer, public :: nprdef = 0
  integer :: nbd = 0
  type(body_t), allocatable :: ent(:)            ! namelist entries (before the csv node table is known)
  integer :: nent = 0
  integer :: cur_b = 0, cur_c = 0, cur_s = 0, cure = 0   ! namelist context: current body / child / sub-child / entry
  real(dp), public :: Motion_Start_T = 0.0_dp   ! free DOFs released from this time on [s]
  logical :: nml_seen = .false.
  logical :: restoring = .false.            ! restart: after_motion applies a saved pose (no flux rotation / fill)
  logical :: nodes_ready = .false.          ! sixdof_precheck done
  ! original grid of the moving blocks
  type :: coord_t
     real(dp), allocatable :: x(:,:,:), y(:,:,:), z(:,:,:), xg(:,:,:,:)
  end type coord_t
  type(coord_t), allocatable :: g0(:)
  real(dp) :: t_reassemble = 0.0_dp
  integer  :: n_reassemble = 0, n_skipped = 0
  ! ---- moving background (update 20) ------------------------------------
  integer  :: Moving_BGR = 0                   ! 0 fixed background, 1 the Object_ID = 0 blocks follow a Body
  integer  :: Moving_BGR_DoF(6) = 1            ! permitted DOFs of the following (surge sway heave roll pitch yaw)
  integer  :: BGR_Following_Body = 0           ! the Body (level 1) whose motion is followed
  real(dp) :: Moving_BGR_Filter_T = 0.0_dp     ! first-order low-pass of the followed motion [s] (0 = exact)
  integer  :: kbg = 0, kfol = 0                ! background node (nbd + 1) and the followed node
  real(dp) :: bg_prev(6) = 0.0_dp              ! followed (masked, filtered) pose of the previous step
  logical  :: bg_first = .true.
  integer  :: iu_bg = -1
  logical  :: bg_open = .false.
  ! re-assembly tolerance (update 20): pose of every node at the last overset assembly
  real(dp), allocatable :: asm_Rg(:,:,:), asm_dg(:,:)
  logical  :: asm_valid = .false.

contains

  !=====================================================================
  ! Namelist group $SixDOF (called from read_namelist for every line of
  ! the group).  Returns ok = .false. for an unknown name.
  !   Body = k                 starts the data of body k
  !   Child_body = c           starts child c of the current body
  !   Sub_child_body = s       starts sub-child s of the current child
  !   Mother_body = k          (child / sub-child) must equal the enclosing body
  !   Mother_child_body = c    (sub-child) must equal the enclosing child
  ! then Mass, CG, Rot_Center, Motion_DOF, Inertia, Prod_Inertia,
  ! Added_Mass, Motion_Damp, Motion_File of the current part.
  !=====================================================================
  subroutine sixdof_nml(nm, vals, nv, ok)
    character(len=*), intent(in) :: nm
    character(len=*), intent(in) :: vals(:)
    integer, intent(in) :: nv
    logical, intent(out) :: ok
    integer :: i, io, n
    real(dp) :: r(8)
    if (.not. allocated(ent)) allocate(ent(MAXB))
    nml_seen = .true.
    ok = .true.
    if (len(nm) >= 4) then
       if (nm(1:4) == 'pid_') then                        ! PID controller of the current part (update 21)
          if (cure == 0) call nml_abort('"'//trim(nm)//'" must follow a Body / Child_body / Sub_child_body line')
          call pid_nml(1, ent(cure)%bid, ent(cure)%cid, ent(cure)%sid, 0, nm, vals, nv, ok)
          return
       end if
    end if
    if (len(nm) > 11) then
       if (nm(1:11) == 'probe_point') then                ! probe point following the current part (update 22)
          if (cure == 0) call nml_abort('"'//trim(nm)//'" must follow a Body / Child_body / Sub_child_body line')
          read(nm(12:), *, iostat=io) n
          if (io /= 0 .or. n < 1) call nml_abort('"'//trim(nm)//'": the point number must follow Probe_point (Probe_point1, Probe_point2, ...)')
          if (nv < 3) call nml_abort('"'//trim(nm)//'" needs x y z [m] (grid frame at t = 0)')
          do i = 1, 3
             read(vals(i), *, iostat=io) r(i)
             if (io /= 0) call nml_abort('"'//trim(nm)//'": value '//trim(vals(i))//' is not a number')
          end do
          if (nprdef >= MAXPRDEF) call nml_abort('too many probe points on the parts (max 256)')
          nprdef = nprdef + 1
          prdef(nprdef)%b = ent(cure)%bid;  prdef(nprdef)%c = ent(cure)%cid;  prdef(nprdef)%s = ent(cure)%sid
          prdef(nprdef)%n = n;  prdef(nprdef)%x = r(1:3)
          return
       end if
    end if
    if (fsi_key(nm)) then                                   ! flexible part (update 24)
       if (cure == 0) call nml_abort('"'//trim(nm)//'" must follow a Body / Child_body / Sub_child_body line')
       call fsi_nml(ent(cure)%fsi, nm, vals, nv, ok)
       return
    end if
    r = 0.0_dp
    do i = 1, min(nv, 8)
       if (nm == 'motion_file') exit
       read(vals(i), *, iostat=io) r(i)
       if (io /= 0) then
          ok = .false.;  return
       end if
    end do
    n = nint(r(1))
    select case (nm)
    case ('body')
       if (n < 1 .or. n > MAXB) call nml_abort('Body must be 1..64')
       cur_b = n;  cur_c = 0;  cur_s = 0
       call select_entry(1, cur_b, 0, 0)
    case ('child_body')
       if (cur_b == 0) call nml_abort('Child_body must follow the Body it belongs to')
       if (n < 1 .or. n > MAXB) call nml_abort('Child_body must be 1..64')
       cur_c = n;  cur_s = 0
       call select_entry(2, cur_b, cur_c, 0)
    case ('sub_child_body')
       if (cur_c == 0) call nml_abort('Sub_child_body must follow the Child_body it belongs to')
       if (n < 1 .or. n > MAXB) call nml_abort('Sub_child_body must be 1..64')
       cur_s = n
       call select_entry(3, cur_b, cur_c, cur_s)
    case ('mother_body')
       if (cure == 0 .or. ent(max(cure,1))%lev == 1) call nml_abort('Mother_body is a property of a Child_body / Sub_child_body')
       if (n /= ent(cure)%bid) call nml_abort('Mother_body differs from the Body the child is written under')
       ent(cure)%mother_b = n
    case ('mother_child_body')
       if (cure == 0 .or. ent(max(cure,1))%lev /= 3) call nml_abort('Mother_child_body is a property of a Sub_child_body')
       if (n /= ent(cure)%cid) call nml_abort('Mother_child_body differs from the Child_body the sub-child is written under')
       ent(cure)%mother_c = n
    case ('motion_start_t');      Motion_Start_T = r(1)
    case ('body_end')                                    ! end marks (update 19): close the current part;
       cure = 0;  cur_c = 0;  cur_s = 0                   ! the Body stays the mother of the Child_body /
    case ('child_body_end')                              ! Force_model groups written after it (update 20:
       cure = 0;  cur_s = 0                               ! before, Child_body after Body_End was refused)
    case ('sub_child_body_end')
       cure = 0
    case default
       if (cure == 0) call nml_abort('"'//trim(nm)//'" must follow a Body / Child_body / Sub_child_body line')
       select case (nm)
       case ('mass');                ent(cure)%mass = r(1)
       case ('cg', 'center_of_gravity');   ent(cure)%cg0 = r(1:3)
       case ('rot_center', 'center_of_rotation');  ent(cure)%rot0 = r(1:3);  ent(cure)%rot_given = .true.
       case ('motion_dof', 'dof');   ent(cure)%dof = nint(r(1:6))
       case ('inertia', 'moment_of_inertia');  ent(cure)%Icg = r(1:3)
       case ('prod_inertia', 'product_of_inertia');  ent(cure)%Ipr = r(1:3)
       case ('added_mass');          ent(cure)%madd = r(1:6)
       case ('motion_damp', 'damping');  ent(cure)%damp = r(1:6)
       case ('motion_file');         ent(cure)%mfile = adjustl(vals(1))
          if (ent(cure)%mfile(1:1) == '"' .or. ent(cure)%mfile(1:1) == "'") &
             ent(cure)%mfile = ent(cure)%mfile(2:len_trim(ent(cure)%mfile)-1)
       case default
          ok = .false.
       end select
    end select
  end subroutine sixdof_nml

  ! the namelist entry of (level, body, child, sub): existing or new
  subroutine select_entry(lev, b, c, s)
    integer, intent(in) :: lev, b, c, s
    integer :: i
    do i = 1, nent
       if (ent(i)%lev == lev .and. ent(i)%bid == b .and. ent(i)%cid == c .and. ent(i)%sid == s) then
          cure = i;  return
       end if
    end do
    if (nent >= MAXB) call nml_abort('more than 64 bodies / children / sub-children')
    nent = nent + 1;  cure = nent
    ent(cure)%lev = lev;  ent(cure)%bid = b;  ent(cure)%cid = c;  ent(cure)%sid = s
    ent(cure)%defined = .true.
  end subroutine select_entry

  subroutine nml_abort(msg)
    character(len=*), intent(in) :: msg
    integer :: io
    if (myrank == 0) write(*,'(A)') ' ERROR (probCFD.nml, $SixDOF): '//msg
    flush(6)
    call MPI_Abort(MPI_COMM_WORLD, 70, io)
  end subroutine nml_abort

  ! is nm one of the flexible-part keys (update 24)?
  logical function fsi_key(nm)
    character(len=*), intent(in) :: nm
    fsi_key = .false.
    select case (nm)
    case ('sw_fsi', 'youngs_modulus', 'poisson_ratio', 'density', 'num_nofsi_group', 'num_fsi_submat_regions', &
          'fsi_start_t', 'damping_ratio', 'fsi_min_thickness', 'fsi_max_iter', 'fsi_tol', 'fsi_relax', 'fsi_mode', &
          'fsi_plane_strain', 'fsi_test_motion')
       fsi_key = .true.
    case default
       if (len(nm) > 6) then
          if (nm(1:6) == 'nofsi_') fsi_key = .true.
       end if
       if (len(nm) > 7) then
          if (nm(1:7) == 'submat_') fsi_key = .true.
       end if
       if (len(nm) > 16) then
          if (nm(1:16) == 'fsi_monitor_node') fsi_key = .true.
       end if
    end select
  end function fsi_key

  ! the flexible-part keys of the current part (update 24)
  subroutine fsi_nml(f, nm, vals, nv, ok)
    type(fsidef_t), intent(inout) :: f
    character(len=*), intent(in) :: nm, vals(:)
    integer, intent(in) :: nv
    logical, intent(out) :: ok
    real(dp) :: r(8)
    integer :: i, io, n, ax, p
    character(len=64) :: key
    ok = .true.
    r = 0.0_dp
    do i = 1, min(nv, 8)
       read(vals(i), *, iostat=io) r(i)
       if (io /= 0) call nml_abort('"'//trim(nm)//'": value '//trim(vals(i))//' is not a number')
    end do
    select case (nm)
    case ('sw_fsi');                 f%sw = nint(r(1))
       if (f%sw /= 0 .and. f%sw /= 1) call nml_abort('SW_FSI must be 0 or 1')
    case ('youngs_modulus');         f%E = r(1)
    case ('poisson_ratio');          f%nu = r(1)
    case ('density');                f%rho = r(1)
    case ('num_nofsi_group');        f%nnof = nint(r(1))
       if (f%nnof < 0 .or. f%nnof > MAXFG) call nml_abort('Num_NoFSI_Group must be 0..16')
    case ('num_fsi_submat_regions'); f%nsub = nint(r(1))
       if (f%nsub < 0 .or. f%nsub > MAXFG) call nml_abort('Num_FSI_Submat_regions must be 0..16')
    case ('fsi_start_t');            f%start_t = r(1)
    case ('damping_ratio');          f%zeta = r(1)
    case ('fsi_min_thickness');      f%tmin_frac = r(1)
    case ('fsi_max_iter');           f%max_iter = nint(r(1))
    case ('fsi_tol');                f%tol = r(1)
    case ('fsi_relax');              f%relax = r(1)
    case ('fsi_mode');               f%mode = nint(r(1))
       if (f%mode /= 0 .and. f%mode /= 1) call nml_abort('FSI_Mode must be 0 (quasi-static) or 1 (dynamic)')
    case ('fsi_plane_strain');       f%plane = nint(r(1))
       if (f%plane < 0 .or. f%plane > 3) call nml_abort('FSI_Plane_Strain must be 0 (off) or the span axis 1 / 2 / 3')
    case ('fsi_test_motion')
       if (nv < 2) call nml_abort('FSI_Test_Motion = A, f [, mode 0 / 1 / 2 [, shape 0 / m]]')
       f%test_amp = r(1);  f%test_freq = r(2)
       if (nv >= 3) f%test_mode = nint(r(3))
       if (nv >= 4) f%test_shp = nint(r(4))
       if (f%test_mode < 0 .or. f%test_mode > 2 .or. (f%test_amp > 0.0_dp .and. f%test_freq <= 0.0_dp)) &
          call nml_abort('FSI_Test_Motion = A, f [, mode 0 / 1 / 2 [, shape 0 / m]] with f > 0')
       if (f%test_shp < 0 .or. f%test_shp > 6) call nml_abort('FSI_Test_Motion: the 4th value (pluck shape) is 0 or a vibration mode 1..6')
    case default
       ! numbered keys: NoFSI_block_num<n>, NoFSI_i_range<n>, Submat_*<n>, FSI_Monitor_node<n>
       p = len_trim(nm)
       do while (p > 1)
          if (nm(p:p) < '0' .or. nm(p:p) > '9') exit
          p = p - 1
       end do
       key = nm(1:p)
       read(nm(p+1:), *, iostat=io) n
       if (io /= 0 .or. n < 1 .or. n > MAXFG) call nml_abort('"'//trim(nm)//'": the group / region number 1..16 must end the name')
       ax = 0
       select case (trim(key))
       case ('nofsi_block_num');  f%nof_blk(n) = nint(r(1))
       case ('nofsi_i_range');    ax = 1
       case ('nofsi_j_range');    ax = 2
       case ('nofsi_k_range');    ax = 3
       case ('submat_block_num'); f%sub_blk(n) = nint(r(1))
       case ('submat_i_range');   ax = -1
       case ('submat_j_range');   ax = -2
       case ('submat_k_range');   ax = -3
       case ('submat_density');   f%sub_rho(n) = r(1)
       case ('submat_youngs_modulus');  f%sub_E(n) = r(1)
       case ('submat_poisson_ratio');   f%sub_nu(n) = r(1)
       case ('fsi_monitor_node')
          if (nv < 4) call nml_abort('"'//trim(nm)//'" needs block i j k')
          f%mon(:,n) = nint(r(1:4));  f%nmon = max(f%nmon, n)
       case default
          ok = .false.;  return
       end select
       if (ax /= 0) then
          if (nv < 1) call nml_abort('"'//trim(nm)//'" needs one index or a range i1, i2')
          if (ax > 0) then
             f%nof_rng(1,ax,n) = nint(r(1));  f%nof_rng(2,ax,n) = nint(r(1))
             if (nv >= 2) f%nof_rng(2,ax,n) = nint(r(2))
          else
             f%sub_rng(1,-ax,n) = nint(r(1));  f%sub_rng(2,-ax,n) = nint(r(1))
             if (nv >= 2) f%sub_rng(2,-ax,n) = nint(r(2))
          end if
       end if
    end select
  end subroutine fsi_nml

  ! name of a part for messages / file names: "body 1", "body 1 child 2", "body 1 child 2 sub-child 1"
  function part_name(k) result(s)
    integer, intent(in) :: k
    character(len=:), allocatable :: s
    character(len=16) :: a, b, c
    write(a,'(I0)') bd(k)%bid;  write(b,'(I0)') bd(k)%cid;  write(c,'(I0)') bd(k)%sid
    select case (bd(k)%lev)
    case (1);  s = 'body '//trim(a)
    case (2);  s = 'body '//trim(a)//' child '//trim(b)
    case default;  s = 'body '//trim(a)//' child '//trim(b)//' sub-child '//trim(c)
    end select
  end function part_name

  ! file-name stem: body<k>, body<k>_child<c>, body<k>_child<c>_sub<s>
  function part_stem(k) result(s)
    integer, intent(in) :: k
    character(len=:), allocatable :: s
    character(len=16) :: a, b, c
    write(a,'(I0)') bd(k)%bid;  write(b,'(I0)') bd(k)%cid;  write(c,'(I0)') bd(k)%sid
    select case (bd(k)%lev)
    case (1);  s = 'body'//trim(a)
    case (2);  s = 'body'//trim(a)//'_child'//trim(b)
    case default;  s = 'body'//trim(a)//'_child'//trim(b)//'_sub'//trim(c)
    end select
  end function part_stem

  ! is node a an ancestor of node k (mother or grand-mother)?
  pure logical function is_ancestor(a, k)
    integer, intent(in) :: a, k
    integer :: m
    is_ancestor = .false.
    m = bd(k)%par
    do while (m > 0)
       if (m == a) then
          is_ancestor = .true.;  return
       end if
       m = bd(m)%par
    end do
  end function is_ancestor

  !=====================================================================
  ! Set-up: the nodes from the csv node table (read_split_csv) +
  ! namelist entries, original grid copies, imposed-motion tables,
  ! output files.  Called after the overset assembly and the initial
  ! condition.
  !=====================================================================
  subroutine sixdof_precheck
    ! the node / name-list consistency, mass properties and imposed-
    ! motion tables (called right after the csv, BEFORE the overset
    ! assembly, so that an input error stops the run at once)
    integer :: k, ierr, i, e, m
    real(dp) :: Ic(3,3)

    nodes_ready = .true.
    nbd = n6node
    if (run_mode == 2) then
       if (nbd == 0 .or. .not. any_imposed_entry()) call stop_clean(5, 'There is no impose motion: Mode = 2 needs a body / '// &
          'child body / sub-child body with an imposed DOF (Motion_DOF = 1) in $SixDOF and its impose motion file.')
    end if
    if (nbd == 0) then
       sixdof_on = .false.;  any_moving = .false.
       if (myrank == 0 .and. nml_seen) write(*,'(A)') &
          '  6-DoF            : $SixDOF group present but every block has object id 0 (row 6 of the csv): no bodies'
       return
    end if
    if (nbd > MAXB) then
       if (myrank == 0) write(*,'(A,I0,A)') ' ERROR (6-DoF): ', nbd, ' parts (bodies + children + sub-children) in the csv (max 64)'
       call MPI_Abort(MPI_COMM_WORLD, 71, ierr)
    end if
    allocate(bd(MAXB))
    if (.not. allocated(ent)) allocate(ent(MAXB))
    sixdof_on = .true.
    wall_velocity_fn => sixdof_wall_velocity
    ! ---- nodes <- namelist entries ---------------------------------------
    do k = 1, nbd
       e = 0
       do i = 1, nent
          if (ent(i)%lev == node_lev(k) .and. ent(i)%bid == node_bid(k) .and. ent(i)%cid == node_cid(k) .and. &
              ent(i)%sid == node_sid(k)) e = i
       end do
       if (e == 0) then
          bd(k)%lev = node_lev(k);  bd(k)%bid = node_bid(k);  bd(k)%cid = node_cid(k);  bd(k)%sid = node_sid(k)
          if (myrank == 0) write(*,'(A)') ' ERROR (6-DoF): no data of '//part_name(k)//' in the $SixDOF group of probCFD.nml'// &
             ' (Body = k / Child_body = c / Sub_child_body = s lines)'
          call MPI_Abort(MPI_COMM_WORLD, 73, ierr)
       end if
       bd(k) = ent(e)
       bd(k)%par = node_par(k)
       bd(k)%nblk = count(blk_obj == k)
    end do
    do i = 1, nent                                      ! namelist parts without blocks
       e = 0
       do k = 1, nbd
          if (ent(i)%lev == bd(k)%lev .and. ent(i)%bid == bd(k)%bid .and. ent(i)%cid == bd(k)%cid .and. ent(i)%sid == bd(k)%sid) e = k
       end do
       if (e == 0) then
          if (myrank == 0) write(*,'(A,I0,A,I0,A,I0,A)') ' ERROR (6-DoF): the $SixDOF part (Body ', ent(i)%bid, ', Child_body ', &
             ent(i)%cid, ', Sub_child_body ', ent(i)%sid, ') has no block in rows 6-8 of the csv'
          call MPI_Abort(MPI_COMM_WORLD, 72, ierr)
       end if
    end do
    do k = 1, nbd
       if (bd(k)%lev >= 2 .and. bd(k)%mother_b == 0) then
          if (myrank == 0) write(*,'(A)') ' ERROR (6-DoF): '//part_name(k)//' needs "Mother_body = '//itoa(bd(k)%bid)//'" in probCFD.nml'
          call MPI_Abort(MPI_COMM_WORLD, 73, ierr)
       end if
       if (bd(k)%lev == 3 .and. bd(k)%mother_c == 0) then
          if (myrank == 0) write(*,'(A)') ' ERROR (6-DoF): '//part_name(k)//' needs "Mother_child_body = '//itoa(bd(k)%cid)//'" in probCFD.nml'
          call MPI_Abort(MPI_COMM_WORLD, 73, ierr)
       end if
       if (any(bd(k)%dof < 0 .or. bd(k)%dof > 2)) then
          if (myrank == 0) write(*,'(A)') ' ERROR (6-DoF): Motion_DOF of '//part_name(k)//' must be 0, 1 or 2 per DOF'
          call MPI_Abort(MPI_COMM_WORLD, 74, ierr)
       end if
       if (any(bd(k)%dof == 2) .and. bd(k)%mass <= 0.0_dp) then
          if (myrank == 0) write(*,'(A)') ' ERROR (6-DoF): '//part_name(k)//' has free DOFs but Mass <= 0'
          call MPI_Abort(MPI_COMM_WORLD, 75, ierr)
       end if
       if (.not. bd(k)%rot_given) bd(k)%rot0 = bd(k)%cg0
       bd(k)%own_moving = any(bd(k)%dof /= 0)
       ! ---- own inertia about the CG in B_k, CG offset -----------------------
       bd(k)%rcg = matmul(Fm, bd(k)%cg0 - bd(k)%rot0)
       Ic = 0.0_dp
       Ic(1,1) = bd(k)%Icg(1);  Ic(2,2) = bd(k)%Icg(2);  Ic(3,3) = bd(k)%Icg(3)
       Ic(1,2) = -bd(k)%Ipr(1);  Ic(2,1) = Ic(1,2)
       Ic(1,3) = -bd(k)%Ipr(2);  Ic(3,1) = Ic(1,3)
       Ic(2,3) = -bd(k)%Ipr(3);  Ic(3,2) = Ic(2,3)
       bd(k)%Ic = Ic
       ! ---- initial offset from the mother (B_mother = E at t = 0) -----------
       m = bd(k)%par
       if (m > 0) bd(k)%r0 = matmul(Fm, bd(k)%rot0 - bd(m)%rot0)
       bd(k)%R = eye3();  bd(k)%Rg = eye3();  bd(k)%dg = 0.0_dp;  bd(k)%Rg_prev = eye3();  bd(k)%dg_prev = 0.0_dp
       bd(k)%pE = matmul(Fm, bd(k)%rot0)
       ! ---- imposed motion table ------------------------------------------
       if (any(bd(k)%dof == 1)) call read_imposed(k)
    end do
    ! PID controllers of the grid parts (update 21): the part rotates about PID_Axis through O
    call pid_setup_parts
    ! a part moves when it or one of its mothers moves (nodes are ordered mothers first)
    do k = 1, nbd
       bd(k)%moving = bd(k)%own_moving .or. (bd(k)%fsi%sw == 1)     ! a flexible part moves its grids (update 24)
       if (bd(k)%par > 0) bd(k)%moving = bd(k)%moving .or. bd(bd(k)%par)%moving
    end do
    any_moving = any(bd(1:nbd)%moving)
    fsi_on = any(bd(1:nbd)%fsi%sw == 1)
    call composite_properties
    call bgr_setup
  end subroutine sixdof_precheck

  !---------------------------------------------------------------------
  ! Moving background (update 20): the blocks with Object_ID = 0 become
  ! the node kbg = nbd + 1 (its own rigid motion group of the overset),
  ! which follows the Body BGR_Following_Body in the DOFs permitted by
  ! Moving_BGR_DoF: its eta = the body's absolute eta with the other
  ! components set to 0 (then the yaw-only rotation is about the
  ! vertical through the body's Rot_Center), optionally low-pass
  ! filtered (Moving_BGR_Filter_T).  Dirichlet velocities on the
  ! background blocks are EARTH-frame values (no grid velocity added):
  ! the far field sees still water / the incident wave.
  !---------------------------------------------------------------------
  subroutine bgr_setup
    integer :: k, b, ierr, nb0
    if (Moving_BGR == 0) return
    if (.not. sixdof_on) then
       if (myrank == 0) write(*,'(A)') ' ERROR (Moving_BGR = 1): no 6-DoF body (every block has Object_ID 0)'
       call MPI_Abort(MPI_COMM_WORLD, 77, ierr)
    end if
    if (any(Moving_BGR_DoF < 0 .or. Moving_BGR_DoF > 1)) then
       if (myrank == 0) write(*,'(A)') ' ERROR (Moving_BGR): Moving_BGR_DoF must be six 0/1 flags (surge sway heave roll pitch yaw)'
       call MPI_Abort(MPI_COMM_WORLD, 77, ierr)
    end if
    if (all(Moving_BGR_DoF == 0)) then
       if (myrank == 0) write(*,'(A)') ' ERROR (Moving_BGR): Moving_BGR_DoF = 0 0 0 0 0 0 - the background would not move; set Moving_BGR = 0'
       call MPI_Abort(MPI_COMM_WORLD, 77, ierr)
    end if
    if (Moving_BGR_Filter_T < 0.0_dp) then
       if (myrank == 0) write(*,'(A)') ' ERROR (Moving_BGR): Moving_BGR_Filter_T must be >= 0'
       call MPI_Abort(MPI_COMM_WORLD, 77, ierr)
    end if
    kfol = 0
    do k = 1, nbd
       if (bd(k)%lev == 1 .and. bd(k)%bid == BGR_Following_Body) kfol = k
    end do
    if (kfol == 0) then
       if (myrank == 0) write(*,'(A,I0,A)') ' ERROR (Moving_BGR): BGR_Following_Body = ', BGR_Following_Body, &
          ' is not a Body of the run (a child / sub-child cannot be followed)'
       call MPI_Abort(MPI_COMM_WORLD, 77, ierr)
    end if
    nb0 = count(blk_obj == 0)
    if (nb0 == 0) then
       if (myrank == 0) write(*,'(A)') ' ERROR (Moving_BGR = 1): no background block (Object_ID = 0) in the csv'
       call MPI_Abort(MPI_COMM_WORLD, 77, ierr)
    end if
    if (nbd + 1 > MAXB) then
       if (myrank == 0) write(*,'(A)') ' ERROR (Moving_BGR): no node left for the background (64 parts)'
       call MPI_Abort(MPI_COMM_WORLD, 77, ierr)
    end if
    kbg = nbd + 1
    bd(kbg)%lev = 0;  bd(kbg)%bid = 0;  bd(kbg)%cid = 0;  bd(kbg)%sid = 0;  bd(kbg)%par = 0
    bd(kbg)%defined = .true.;  bd(kbg)%dof = 0;  bd(kbg)%mass = 0.0_dp
    bd(kbg)%rot0 = bd(kfol)%rot0;  bd(kbg)%cg0 = bd(kfol)%rot0;  bd(kbg)%rot_given = .true.
    bd(kbg)%own_moving = .false.;  bd(kbg)%moving = .true.
    bd(kbg)%R = eye3();  bd(kbg)%Rg = eye3();  bd(kbg)%dg = 0.0_dp;  bd(kbg)%Rg_prev = eye3();  bd(kbg)%dg_prev = 0.0_dp
    bd(kbg)%pE = matmul(Fm, bd(kbg)%rot0)
    bd(kbg)%nblk = nb0
    do b = 1, nblk
       if (blk_obj(b) == 0) blk_obj(b) = kbg
    end do
    any_moving = .true.
    bgr_on = .true.;  bgr_O0 = bd(kbg)%rot0;  bgr_Rg = eye3();  bgr_dg = 0.0_dp
    if (myrank == 0) then
       write(*,'(A)') ' '
       write(*,'(A,I0,A,I0,A,6I2,A)') '  Moving background: ', nb0, ' block(s) with Object_ID 0 follow body ', BGR_Following_Body, &
          ' in the DOFs (surge sway heave roll pitch yaw) ', Moving_BGR_DoF, '  (update 20)'
       if (Moving_BGR_Filter_T > 0.0_dp) then
          write(*,'(A,ES10.3,A)') '    followed motion low-pass filtered, Moving_BGR_Filter_T = ', Moving_BGR_Filter_T, ' s'
       else
          write(*,'(A)') '    followed motion not filtered (exact body motion in the permitted DOFs)'
       end if
       write(*,'(A,3F10.4,A)') '    rotation centre of the background = Rot_Center of the body ', bd(kbg)%rot0, ' m (grid frame)'
       write(*,'(A)') '    Dirichlet velocities / wave values on the background blocks are earth-frame values (no grid velocity added)'
    end if
  end subroutine bgr_setup

  ! the followed pose of the background for the step ending at t (after the followed body moved)
  subroutine bgr_follow(dtt, init)
    real(dp), intent(in) :: dtt
    logical,  intent(in) :: init
    real(dp) :: raw(6), rawd(6), a, f(6)
    integer :: c
    raw = 0.0_dp;  rawd = 0.0_dp
    do c = 1, 6
       if (Moving_BGR_DoF(c) == 1) then
          raw(c) = bd(kfol)%eta(c);  rawd(c) = bd(kfol)%etad(c)
       end if
    end do
    if (init .or. bg_first .or. Moving_BGR_Filter_T <= 0.0_dp) then
       f = raw
       if (init .or. bg_first) then
          bd(kbg)%etad = rawd
       else
          bd(kbg)%etad = rawd
       end if
       bg_first = .false.
    else
       a = exp(-dtt/Moving_BGR_Filter_T)
       f = a*bd(kbg)%eta + (1.0_dp - a)*raw
       bd(kbg)%etad = (f - bd(kbg)%eta)/dtt
    end if
    bd(kbg)%eta = f;  bd(kbg)%etadd = 0.0_dp;  bd(kbg)%nud = 0.0_dp
  end subroutine bgr_follow

  ! has the background pose changed since the last call (grids must move)?
  logical function bgr_changed()
    bgr_changed = any(abs(bd(kbg)%eta - bg_prev) > 0.0_dp)
    bg_prev = bd(kbg)%eta
  end function bgr_changed

  subroutine bgr_publish
    bgr_Rg = bd(kbg)%Rg;  bgr_dg = bd(kbg)%dg;  bgr_O0 = bd(kbg)%rot0
  end subroutine bgr_publish

  subroutine sixdof_init
    integer :: b, k, ib, d, j, n
    logical :: ex
    if (.not. nodes_ready) call sixdof_precheck
    if (.not. sixdof_on) return

    ! ---- original grid of the moving blocks -----------------------------
    allocate(g0(nblk))
    do b = 1, nblk
       k = blk_obj(b)
       if (k == 0) cycle
       if (.not. bd(k)%moving) cycle
       allocate(g0(b)%x(blk(b)%ni, blk(b)%nj, blk(b)%nk), g0(b)%y(blk(b)%ni, blk(b)%nj, blk(b)%nk), &
                g0(b)%z(blk(b)%ni, blk(b)%nj, blk(b)%nk))
       g0(b)%x = blk(b)%xn;  g0(b)%y = blk(b)%yn;  g0(b)%z = blk(b)%zn
       allocate(g0(b)%xg(3, 1-ng:blk(b)%ni+ng, 1-ng:blk(b)%nj+ng, 1-ng:blk(b)%nk+ng))
       g0(b)%xg = blk(b)%xg
    end do
    do ib = 1, nown
       d = own_list(ib);  b = dom(d)%b;  k = blk_obj(b)
       if (k == 0) cycle
       if (.not. bd(k)%moving) cycle
       dom(d)%moving = .true.
       allocate(dom(d)%Xh0(3, dom(d)%hlo(1):dom(d)%hhi(1), dom(d)%hlo(2):dom(d)%hhi(2), dom(d)%hlo(3):dom(d)%hhi(3)))
       allocate(dom(d)%mh0(9, dom(d)%hlo(1):dom(d)%hhi(1), dom(d)%hlo(2):dom(d)%hhi(2), dom(d)%hlo(3):dom(d)%hhi(3)))
       dom(d)%Xh0 = dom(d)%Xh;  dom(d)%mh0 = dom(d)%mh
       allocate(dom(d)%Ug1(lbound(dom(d)%U1,1):ubound(dom(d)%U1,1), lbound(dom(d)%U1,2):ubound(dom(d)%U1,2), &
                           lbound(dom(d)%U1,3):ubound(dom(d)%U1,3)))
       allocate(dom(d)%Ug2, dom(d)%Ug3, mold=dom(d)%Ug1)
       dom(d)%Ug1 = 0.0_dp;  dom(d)%Ug2 = 0.0_dp;  dom(d)%Ug3 = 0.0_dp
    end do

    ! ---- output files (folder Output_data/) -------------------------------
    if (myrank == 0) then
       inquire(file=out_dir//'.', exist=ex)
       if (.not. ex) call execute_command_line('mkdir -p '//out_dir, wait=.true.)
       do k = 1, nbd                                        ! PID controllers of the parts (update 21)
          do j = 1, bd(k)%npid
             call pid_open(bd(k)%ipid(j), part_stem(k), part_name(k))
          end do
       end do
       do k = 1, nbd
          if (bd(k)%lev == 1) then
             call open_out(k, bd(k)%iu(1), out_dir//trim(casename)//'_force_'//part_stem(k)//'.dat',  'force',  'body')
             call open_out(k, bd(k)%iu(2), out_dir//trim(casename)//'_moment_'//part_stem(k)//'.dat', 'moment', 'body')
             call open_out(k, bd(k)%iu(3), out_dir//trim(casename)//'_force_earth'//itoa(bd(k)%bid)//'.dat', 'force',  'earth')
             call open_out(k, bd(k)%iu(4), out_dir//trim(casename)//'_moment_earth'//itoa(bd(k)%bid)//'.dat','moment', 'earth')
             call open_motion(k, out_dir//trim(casename)//'_motion'//itoa(bd(k)%bid)//'.dat')
          else
             call open_out(k, bd(k)%iu(1), out_dir//trim(casename)//'_'//part_stem(k)//'_force.dat',  'force',  'mother')
             call open_out(k, bd(k)%iu(2), out_dir//trim(casename)//'_'//part_stem(k)//'_moment.dat', 'moment', 'mother')
             call open_motion(k, out_dir//trim(casename)//'_'//part_stem(k)//'_motion.dat')
          end if
       end do
    end if

    if (myrank == 0) then
       write(*,'(A)') ' '
       write(*,'(A)') ' ---------------- 6-DoF rigid-body motion ----------------'
       write(*,'(A,I0,A,I0,A,I0,A)') '  parts            : ', count(bd(1:nbd)%lev == 1), ' body(ies), ', count(bd(1:nbd)%lev == 2), &
          ' child body(ies), ', count(bd(1:nbd)%lev == 3), ' sub-child body(ies) (rows 6-8 of Input_polyMesh_probCFD.csv)'
       do k = 1, nbd
          write(*,'(A,I0,A,6I2,A)') '  '//part_name(k)//': ', bd(k)%nblk, ' blocks, Motion_DOF (surge sway heave roll pitch yaw) = ', bd(k)%dof, &
             merge('  moving', '  static', bd(k)%moving)
          if (bd(k)%lev >= 2) write(*,'(A,3F10.4,A)') '    motion RELATIVE to '//part_name(bd(k)%par)//'; initial offset of O in the mother frame [m]', &
             bd(k)%r0, '  (SNAME: x fwd, y stbd, z down)'
          write(*,'(A,ES12.4,A,3F10.4,A,3F10.4)') '    mass [kg] ', bd(k)%mass, '   CG (grid) [m]', bd(k)%cg0, '   O (grid) [m]', bd(k)%rot0
          write(*,'(A,3ES11.3,A,3ES11.3)') '    I about CG, own axes [kg m2]: Ixx Iyy Izz', bd(k)%Icg, '  Ixy Ixz Iyz', bd(k)%Ipr
          if (bd(k)%Mc /= bd(k)%mass) write(*,'(A,ES12.4,A,3ES11.3)') '    composite with its children: mass ', bd(k)%Mc, &
             '  I about O (diagonal)', bd(k)%IOc(1,1), bd(k)%IOc(2,2), bd(k)%IOc(3,3)
          if (any(bd(k)%madd /= 0.0_dp)) write(*,'(A,6ES10.2)') '    added-mass stabilization        ', bd(k)%madd
          if (any(bd(k)%damp /= 0.0_dp)) write(*,'(A,6ES10.2)') '    artificial damping [1/s]        ', bd(k)%damp
          if (any(bd(k)%dof == 1)) write(*,'(A,A,A,I0,A)') '    imposed motion from ', trim(bd(k)%mfile), ' (', bd(k)%nt, ' rows)'
          if (bd(k)%pid_rot) write(*,'(A,3F9.4,A)') '    rotated by its PID controller about the axis', &
             matmul(Fm, bd(k)%pid_axB), ' (grid frame at t = 0) through O (update 21)'
       end do
       call pid_summary
       if (any(bd(1:nbd)%dof(1) == 2 .or. bd(1:nbd)%dof(2) == 2 .or. bd(1:nbd)%dof(3) == 2 .or. &
               bd(1:nbd)%dof(4) == 2 .or. bd(1:nbd)%dof(5) == 2 .or. bd(1:nbd)%dof(6) == 2)) &
          write(*,'(A,ES12.4,A)') '  free DOFs released at t >= ', Motion_Start_T, ' s (Motion_Start_T)'
       write(*,'(A,I0,A)') '  data files every ', data_interval, ' steps (Data_Intvl) in '//out_dir// &
          ': <case>_{force,moment}_{body,earth}<k>.dat, <case>_motion<k>.dat, <case>_body<k>_child<c>[_sub<s>]_{force,moment,motion}.dat'
    end if
    ! imposed DOFs: initial position of the table (t = 0); a non-zero
    ! initial value displaces the grid before the first step (a child of
    ! a displaced mother moves with it)
    n = 0
    do k = 1, nbd
       if (any(bd(k)%dof == 1)) then
          call imposed_at(k, 0.0_dp, bd(k)%eta, bd(k)%etad, bd(k)%etadd)
          if (any(abs(bd(k)%eta) > 0.0_dp)) n = n + 1
       end if
       call update_transform(k)
       bd(k)%Rg_prev = bd(k)%Rg;  bd(k)%dg_prev = bd(k)%dg
    end do
    if (kbg > 0) then                            ! moving background (update 20): initial followed pose
       call bgr_follow(dt_dim, .true.)
       call update_transform(kbg);  bd(kbg)%Rg_prev = bd(kbg)%Rg;  bd(kbg)%dg_prev = bd(kbg)%dg
       if (bgr_changed()) n = n + 1
       call bgr_publish
       if (myrank == 0) then
          open(newunit=iu_bg, file=out_dir//trim(casename)//'_bgr_motion.dat', status=merge('unknown', 'replace', run_mode == 1), &
               position='append', action='write')
          bg_open = .true.
          if (run_mode /= 1) then
             write(iu_bg,'(A,I0,A,6I2)') '# probCFD moving background (update 20): follows body ', BGR_Following_Body, &
                ', permitted DOFs (surge sway heave roll pitch yaw) ', Moving_BGR_DoF
             write(iu_bg,'(A)') '#   eta = displacement of the background rotation centre in the SNAME earth frame [m] and its'// &
                ' Euler angles [deg] (filtered when Moving_BGR_Filter_T > 0), eta_dot [m/s, deg/s], O = rotation centre in the grid frame [m]'
             write(iu_bg,'(A)') '#   step      time[s]            x               y               z            roll           pitch'// &
                '             yaw            xdot            ydot            zdot         rolldot        pitchdot          yawdot'// &
                '              Ox              Oy              Oz'
          end if
       end if
    end if
    if (n > 0) then
       call move_grids
       call after_motion
    end if
    j = 0
  end subroutine sixdof_init

  ! any name-list part with an imposed DOF?
  logical function any_imposed_entry()
    integer :: i
    any_imposed_entry = .false.
    do i = 1, nent
       if (any(ent(i)%dof == 1)) any_imposed_entry = .true.
    end do
  end function any_imposed_entry

  !=====================================================================
  ! Mode = 2: GRID-MOVEMENT CHECK of the imposed motions (update 14b).
  ! No flow solution, no overset assembly: every part with imposed DOFs
  ! follows its table (free DOFs are held, fixed children ride with
  ! their mothers), the node coordinates of the moving blocks are moved
  ! rigidly and the grid is written every Output_Intvl steps (and at the
  ! last step) as Check_grid_movement/<case>_NNNNNNNN.grd + a companion
  ! PLOT3D solution file <case>_NNNNNNNN.q carrying the time (solution-
  ! time bar of Tecplot), with the relative / absolute attitude of every
  ! imposed part in Check_grid_movement/<case>_movement.dat.  Called right after the csv
  ! (no field arrays: it runs for any grid size), then the run stops.
  !=====================================================================
  subroutine sixdof_check_movement
    character(len=*), parameter :: cdir = 'Check_grid_movement/'
    integer :: k, b, istep, nst, nout, iu, ierr, nimp
    real(dp) :: t, r2d, ea(3)
    logical :: ex
    character(len=32) :: sk
    if (.not. nodes_ready) call sixdof_precheck
    r2d = 180.0_dp/pi
    nimp = 0
    do k = 1, nbd
       if (any(bd(k)%dof == 1)) nimp = nimp + 1
    end do
    ! original node coordinates of the moving blocks
    allocate(g0(nblk))
    do b = 1, nblk
       k = blk_obj(b)
       if (k == 0) cycle
       if (.not. bd(k)%moving) cycle
       allocate(g0(b)%x(blk(b)%ni, blk(b)%nj, blk(b)%nk), g0(b)%y(blk(b)%ni, blk(b)%nj, blk(b)%nk), &
                g0(b)%z(blk(b)%ni, blk(b)%nj, blk(b)%nk))
       g0(b)%x = blk(b)%xn;  g0(b)%y = blk(b)%yn;  g0(b)%z = blk(b)%zn
    end do
    nst = max(1, nint(End_T/dt_dim))
    if (myrank == 0) then
       write(*,'(A)') ' '
       write(*,'(A)') ' ============== Mode = 2: grid-movement check of the imposed motions =============='
       write(*,'(A,I0,A,ES11.4,A,ES11.4,A)') '  steps            : ', nst, ' (Del_T = ', dt_dim, ' s, End_T = ', End_T, ' s)'
       write(*,'(A,I0,A)') '  grid export      : every ', out_interval, ' steps (Output_Intvl) + the last step -> '//cdir// &
          trim(casename)//'_NNNNNNNN.grd + .q (PLOT3D solution file with TIME = t: load the pairs together in Tecplot)'
       write(*,'(A,I0,A,I0,A)') '  imposed parts    : ', nimp, ' of ', nbd, ' (free DOFs are held at zero; fixed children move with their mothers)'
       do k = 1, nbd
          if (any(bd(k)%dof == 1)) write(*,'(A,6I2,A)') '    '//part_name(k)//': Motion_DOF ', bd(k)%dof, '  table '//trim(bd(k)%mfile)
       end do
       inquire(file=cdir//'.', exist=ex)
       if (.not. ex) call execute_command_line('mkdir -p '//cdir, wait=.true.)
       open(newunit=iu, file=cdir//trim(casename)//'_movement.dat', status='replace', action='write')
       write(iu,'(A)') '# probCFD Mode = 2 grid-movement check: attitude of every imposed part at the exported steps'
       write(iu,'(A)') '# per part: x y z [m] (displacement of O, in the mother frame; a body: earth frame), roll pitch yaw [deg]'
       write(iu,'(A)') '#   relative to the mother, then ABSOLUTE roll pitch yaw [deg], then O position in the grid frame [m]'
       write(iu,'(A)',advance='no') '# columns: 1 step  2 time[s]'
       nout = 3
       do k = 1, nbd
          if (.not. any(bd(k)%dof == 1)) cycle
          write(sk,'(I0,A,I0)') nout, '-', nout + 11
          write(iu,'(A)',advance='no') '  '//trim(sk)//' '//part_name(k)
          nout = nout + 12
       end do
       write(iu,'(A)') ' '
    end if
    nout = 0
    do istep = 0, nst
       t = istep*dt_dim
       do k = 1, nbd
          if (any(bd(k)%dof == 1)) call imposed_at(k, t, bd(k)%eta, bd(k)%etad, bd(k)%etadd)
          call update_transform(k)
       end do
       if (mod(istep, out_interval) /= 0 .and. istep /= nst) cycle
       call move_nodes_only
       call write_grid_to(istep, cdir)
       call write_check_q(istep, t, cdir)
       nout = nout + 1
       if (myrank == 0) then
          write(iu,'(I8,ES16.8)',advance='no') istep, t
          write(*,'(A,I8,A,ES12.5,A)',advance='no') '  step ', istep, '  t = ', t, ' s :'
          do k = 1, nbd
             if (.not. any(bd(k)%dof == 1)) cycle
             ea = euler_of(bd(k)%R)
             write(iu,'(12ES16.8)',advance='no') bd(k)%eta(1:3), bd(k)%eta(4:6)*r2d, ea*r2d, bd(k)%rot0 + bd(k)%dg
             write(*,'(A,3F9.3,A)',advance='no') '  '//part_name(k)//' rel. roll/pitch/yaw', bd(k)%eta(4:6)*r2d, ' deg'
          end do
          write(iu,'(A)') ' ';  write(*,'(A)') '  -> grid written'
          flush(iu)
       end if
    end do
    if (myrank == 0) close(iu)
    write(sk,'(I0)') nout
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    call stop_clean(6, 'Mode = 2: grid-movement check completed, '//trim(sk)//' grid files written to '//cdir// &
       ' ('//trim(casename)//'_NNNNNNNN.grd + .q pairs, attitudes in '//trim(casename)//'_movement.dat). Set Mode = 0 for the flow simulation.')
  end subroutine sixdof_check_movement

  ! Mode = 2 companion PLOT3D SOLUTION file (q format, 5 variables, with
  ! the time in the block headers) of every exported grid: Tecplot's
  ! PLOT3D loader takes the solution time of a grid / solution pair from
  ! the TIME header value, so the pairs <case>_NNNNNNNN.grd / .q loaded
  ! together animate with the solution-time bar.  Variables (named in
  ! <case>.nam): Part_ID (0 background, k = 6-DoF part / node of the
  ! block), dX dY dZ (displacement of every node from the initial grid
  ! [m]) and Disp (its magnitude).  Header: FSMACH = 0, ALPHA = 0,
  ! RE = step, TIME = t [s].
  subroutine write_check_q(istep, t, dir)
    integer,  intent(in) :: istep
    real(dp), intent(in) :: t
    character(len=*), intent(in) :: dir
    integer :: iu, b, i, j, k, n, ni, nj, nk
    real(4), allocatable :: q(:,:,:,:)
    character(len=256) :: fname
    logical :: ex
    if (myrank /= 0) return
    inquire(file=trim(dir)//trim(casename)//'.nam', exist=ex)
    if (.not. ex) then
       open(newunit=iu, file=trim(dir)//trim(casename)//'.nam', status='replace', action='write')
       write(iu,'(A)') 'Part_ID';  write(iu,'(A)') 'dX';  write(iu,'(A)') 'dY';  write(iu,'(A)') 'dZ';  write(iu,'(A)') 'Disp'
       close(iu)
    end if
    write(fname,'(A,A,A,I8.8,A)') trim(dir), trim(casename), '_', istep, '.q'
    open(newunit=iu, file=trim(fname), form='unformatted', status='replace', action='write')
    write(iu) nblk
    write(iu) (blk(b)%ni, blk(b)%nj, blk(b)%nk, b=1,nblk)
    do b = 1, nblk
       ni = blk(b)%ni;  nj = blk(b)%nj;  nk = blk(b)%nk
       allocate(q(ni,nj,nk,5));  q = 0.0
       q(:,:,:,1) = real(blk_obj(b), 4)
       if (allocated(g0(b)%x)) then
          q(:,:,:,2) = real(blk(b)%xn - g0(b)%x, 4)
          q(:,:,:,3) = real(blk(b)%yn - g0(b)%y, 4)
          q(:,:,:,4) = real(blk(b)%zn - g0(b)%z, 4)
          q(:,:,:,5) = sqrt(q(:,:,:,2)**2 + q(:,:,:,3)**2 + q(:,:,:,4)**2)
       end if
       write(iu) 0.0_4, 0.0_4, real(istep, 4), real(t, 4)        ! FSMACH, ALPHA, RE (= step), TIME [s]
       write(iu) ((((q(i,j,k,n), i=1,ni), j=1,nj), k=1,nk), n=1,5)
       deallocate(q)
    end do
    close(iu)
  end subroutine write_check_q

  ! node coordinates of the moving blocks from the original grid (Mode = 2:
  ! no half-index coordinates / metrics / grid velocities exist)
  subroutine move_nodes_only
    integer :: b, k, i, j, kk
    real(dp) :: Rg(3,3), c0(3), dg(3), x(3)
    do b = 1, nblk
       k = blk_obj(b)
       if (k == 0) cycle
       if (.not. bd(k)%moving) cycle
       Rg = bd(k)%Rg;  c0 = bd(k)%rot0;  dg = bd(k)%dg
       do kk = 1, blk(b)%nk
        do j = 1, blk(b)%nj
         do i = 1, blk(b)%ni
            x = matmul(Rg, [g0(b)%x(i,j,kk), g0(b)%y(i,j,kk), g0(b)%z(i,j,kk)] - c0) + c0 + dg
            blk(b)%xn(i,j,kk) = x(1);  blk(b)%yn(i,j,kk) = x(2);  blk(b)%zn(i,j,kk) = x(3)
         end do
        end do
       end do
    end do
  end subroutine move_nodes_only

  function itoa(i) result(s)
    integer, intent(in) :: i
    character(len=:), allocatable :: s
    character(len=16) :: buf
    write(buf,'(I0)') i
    s = trim(buf)
  end function itoa

  subroutine open_out(k, iu, fname, what, frame)
    integer, intent(in)  :: k
    integer, intent(out) :: iu
    character(len=*), intent(in) :: fname, what, frame
    if (run_mode == 1) then                   ! restart: continue the existing file
       open(newunit=iu, file=fname, status='unknown', position='append', action='write')
       return
    end if
    open(newunit=iu, file=fname, status='replace', action='write')
    if (what == 'force') then
       write(iu,'(A)') '# probCFD 6-DoF: fluid force on '//part_name(k)//' (its own walls) [N], '//frame//'-fixed frame'
    else
       write(iu,'(A)') '# probCFD 6-DoF: fluid moment on '//part_name(k)//' (its own walls) about its centre of rotation O [N m], '// &
          frame//'-fixed frame'
    end if
    select case (frame)
    case ('body')
       write(iu,'(A)') '# body frame (SNAME): x forward, y starboard, z down, aligned with the earth frame at t = 0'
    case ('mother')
       write(iu,'(A)') '# frame of the mother ('//part_name(bd(k)%par)//'): x forward, y starboard, z down, aligned with the earth frame at t = 0'
    case default
       write(iu,'(A)') '# earth frame (SNAME): x forward (= -x grid), y starboard (= y grid), z down (= -z grid)'
    end select
    if (what == 'force') then
       write(iu,'(A)') '# columns: 1 step  2 time[s]  3-5 friction X Y Z  6-8 pressure X Y Z  9-11 total X Y Z   (weight NOT included)'
    else
       write(iu,'(A)') '# columns: 1 step  2 time[s]  3-5 friction K M N  6-8 pressure K M N  9-11 total K M N   (weight NOT included)'
    end if
  end subroutine open_out

  subroutine open_motion(k, fname)
    integer, intent(in) :: k
    character(len=*), intent(in) :: fname
    integer :: iu
    if (run_mode == 1) then
       open(newunit=bd(k)%iu(5), file=fname, status='unknown', position='append', action='write')
       return
    end if
    open(newunit=bd(k)%iu(5), file=fname, status='replace', action='write')
    iu = bd(k)%iu(5)
    write(iu,'(A)') '# probCFD 6-DoF motion of '//part_name(k)//' (SNAME / Fossen convention)'
    write(iu,'(A)') '# earth frame E: x forward (= -x grid), y starboard (= y grid), z down (= -z grid); origin = grid origin'
    if (bd(k)%lev == 1) then
       write(iu,'(A)') '# body frame B: aligned with E at t = 0, origin at the centre of rotation O'
       write(iu,'(A)') '# Euler angles Z-Y-X (yaw, pitch, roll): R = Rz(psi) Ry(theta) Rx(phi); positive roll = starboard down,'
       write(iu,'(A)') '#   positive pitch = bow up, positive yaw = bow to starboard; heave positive DOWN'
       write(iu,'(A)') '# columns: 1 step  2 time[s]  3-5 surge sway heave = displacement of O in E [m]  6-8 roll pitch yaw [deg]'
       write(iu,'(A)') '#   9-11 u v w body velocity [m/s]  12-14 p q r body angular velocity [deg/s]'
       write(iu,'(A)') '#   15-17 xdot ydot zdot velocity of O in E [m/s]  18-20 Euler-angle rates [deg/s]'
       write(iu,'(A)') '#   21-23 udot vdot wdot body acceleration [m/s2]  24-26 pdot qdot rdot [deg/s2]'
       write(iu,'(A)') '#   27-29 position of O in the GRID frame [m]  30-32 position of the CG in the GRID frame [m]'
    else
       write(iu,'(A)') '# ALL motion quantities are RELATIVE to the mother ('//part_name(bd(k)%par)//'): the mother frame B_m is'
       write(iu,'(A)') '#   aligned with E at t = 0 (x forward, y starboard, z down); the own frame B_k is aligned with B_m at t = 0,'
       write(iu,'(A)') '#   origin at the own centre of rotation O; R_k = R_m Rz(psi) Ry(theta) Rx(phi) (relative Euler angles)'
       write(iu,'(A)') '# columns: 1 step  2 time[s]  3-5 displacement of O from its initial position, in the MOTHER frame [m]'
       write(iu,'(A)') '#   6-8 relative roll pitch yaw [deg]  9-11 relative velocity of O in the own frame [m/s]'
       write(iu,'(A)') '#   12-14 relative angular velocity p q r in the own frame [deg/s]'
       write(iu,'(A)') '#   15-17 velocity of O relative to the mother, in the mother frame [m/s]  18-20 relative Euler-angle rates [deg/s]'
       write(iu,'(A)') '#   21-23 relative acceleration of O in the own frame [m/s2]  24-26 pdot qdot rdot relative [deg/s2]'
       write(iu,'(A)') '#   27-29 position of O in the GRID frame [m]  30-32 position of the CG in the GRID frame [m]'
       write(iu,'(A)') '#   33-35 ABSOLUTE roll pitch yaw of the part [deg]  36-38 ABSOLUTE angular velocity in the own frame [deg/s]'
    end if
  end subroutine open_motion

  !=====================================================================
  ! Imposed-motion file of part k: "t x y z phi theta psi" per row.
  !=====================================================================
  subroutine read_imposed(k)
    integer, intent(in) :: k
    integer :: iu, ios, n, ierr, i
    character(len=512) :: line
    real(dp) :: v(7)
    if (len_trim(bd(k)%mfile) == 0) bd(k)%mfile = trim(casename)//'_'//part_stem(k)//'_imposed.dat'
    open(newunit=iu, file=trim(bd(k)%mfile), status='old', action='read', iostat=ios)
    if (ios /= 0 .and. run_mode == 2) call stop_clean(5, 'There is no impose motion file: '//trim(bd(k)%mfile)// &
       ' (imposed DOFs of '//part_name(k)//').')
    if (ios /= 0) then
       if (myrank == 0) write(*,'(A)') ' ERROR (6-DoF): '//part_name(k)//' has imposed DOFs but the motion file '// &
          trim(bd(k)%mfile)//' cannot be opened'
       call MPI_Abort(MPI_COMM_WORLD, 76, ierr)
    end if
    n = 0
    do
       read(iu,'(A)', iostat=ios) line
       if (ios /= 0) exit
       if (len_trim(line) == 0) cycle
       line = adjustl(line)
       if (line(1:1) == '#' .or. line(1:1) == '!') cycle
       read(line, *, iostat=ios) v
       if (ios /= 0) cycle
       n = n + 1
    end do
    if (n < 1) then
       if (myrank == 0) write(*,'(A)') ' ERROR (6-DoF): no data rows in '//trim(bd(k)%mfile)
       call MPI_Abort(MPI_COMM_WORLD, 77, ierr)
    end if
    allocate(bd(k)%tt(n), bd(k)%et(6,n));  bd(k)%nt = n
    rewind(iu);  n = 0
    do
       read(iu,'(A)', iostat=ios) line
       if (ios /= 0) exit
       if (len_trim(line) == 0) cycle
       line = adjustl(line)
       if (line(1:1) == '#' .or. line(1:1) == '!') cycle
       read(line, *, iostat=ios) v
       if (ios /= 0) cycle
       n = n + 1
       bd(k)%tt(n) = v(1);  bd(k)%et(1:3,n) = v(2:4);  bd(k)%et(4:6,n) = v(5:7)*pi/180.0_dp
    end do
    close(iu)
    do i = 2, bd(k)%nt
       if (bd(k)%tt(i) <= bd(k)%tt(i-1)) then
          if (myrank == 0) write(*,'(A)') ' ERROR (6-DoF): the times in '//trim(bd(k)%mfile)//' must increase'
          call MPI_Abort(MPI_COMM_WORLD, 78, ierr)
       end if
    end do
  end subroutine read_imposed

  !---------------------------------------------------------------------
  ! Imposed eta (and its first / second time derivatives) of body k at
  ! time t: cubic Hermite (Catmull-Rom) interpolation of the table,
  ! held constant outside its time range.  Only the imposed components
  ! of e, ed, edd are written (the others are untouched) unless all_ is set.
  !---------------------------------------------------------------------
  subroutine imposed_at(k, t, e, ed, edd)
    integer,  intent(in)    :: k
    real(dp), intent(in)    :: t
    real(dp), intent(inout) :: e(6), ed(6), edd(6)
    integer :: n, i, c
    real(dp) :: s, h, p0, p1, m0, m1, h00, h10, h01, h11, d00, d10, d01, d11, dd00, dd10, dd01, dd11
    n = bd(k)%nt
    do c = 1, 6
       if (bd(k)%dof(c) /= 1) cycle
       if (n == 1 .or. t <= bd(k)%tt(1)) then
          e(c) = bd(k)%et(c,1);  ed(c) = 0.0_dp;  edd(c) = 0.0_dp;  cycle
       end if
       if (t >= bd(k)%tt(n)) then
          e(c) = bd(k)%et(c,n);  ed(c) = 0.0_dp;  edd(c) = 0.0_dp;  cycle
       end if
       i = 1
       do while (i < n-1 .and. t > bd(k)%tt(i+1))
          i = i + 1
       end do
       h = bd(k)%tt(i+1) - bd(k)%tt(i);  s = (t - bd(k)%tt(i))/h
       p0 = bd(k)%et(c,i);  p1 = bd(k)%et(c,i+1)
       if (i > 1) then
          m0 = (p1 - bd(k)%et(c,i-1))/(bd(k)%tt(i+1) - bd(k)%tt(i-1))*h
       else
          m0 = p1 - p0
       end if
       if (i < n-1) then
          m1 = (bd(k)%et(c,i+2) - p0)/(bd(k)%tt(i+2) - bd(k)%tt(i))*h
       else
          m1 = p1 - p0
       end if
       h00 = 2*s**3 - 3*s**2 + 1;  h10 = s**3 - 2*s**2 + s;  h01 = -2*s**3 + 3*s**2;  h11 = s**3 - s**2
       d00 = 6*s**2 - 6*s;  d10 = 3*s**2 - 4*s + 1;  d01 = -6*s**2 + 6*s;  d11 = 3*s**2 - 2*s
       dd00 = 12*s - 6;  dd10 = 6*s - 4;  dd01 = -12*s + 6;  dd11 = 6*s - 2
       e(c)   = h00*p0 + h10*m0 + h01*p1 + h11*m1
       ed(c)  = (d00*p0 + d10*m0 + d01*p1 + d11*m1)/h
       edd(c) = (dd00*p0 + dd10*m0 + dd01*p1 + dd11*m1)/h**2
    end do
  end subroutine imposed_at

  !=====================================================================
  ! Small linear algebra
  !=====================================================================
  pure function eye3() result(e)
    real(dp) :: e(3,3)
    e = 0.0_dp;  e(1,1) = 1.0_dp;  e(2,2) = 1.0_dp;  e(3,3) = 1.0_dp
  end function eye3

  pure function skew(a) result(s)
    real(dp), intent(in) :: a(3)
    real(dp) :: s(3,3)
    s = 0.0_dp
    s(1,2) = -a(3);  s(1,3) =  a(2)
    s(2,1) =  a(3);  s(2,3) = -a(1)
    s(3,1) = -a(2);  s(3,2) =  a(1)
  end function skew

  pure function cross3(a, b) result(c)
    real(dp), intent(in) :: a(3), b(3)
    real(dp) :: c(3)
    c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
  end function cross3

  ! R = Rz(psi) Ry(theta) Rx(phi): body -> earth
  pure function rotmat(e2) result(R)
    real(dp), intent(in) :: e2(3)
    real(dp) :: R(3,3), cf, sf, ct, st, cp, sp
    cf = cos(e2(1));  sf = sin(e2(1));  ct = cos(e2(2));  st = sin(e2(2));  cp = cos(e2(3));  sp = sin(e2(3))
    R(1,1) = cp*ct;  R(1,2) = -sp*cf + cp*st*sf;  R(1,3) =  sp*sf + cp*cf*st
    R(2,1) = sp*ct;  R(2,2) =  cp*cf + sf*st*sp;  R(2,3) = -cp*sf + st*sp*cf
    R(3,1) = -st;    R(3,2) =  ct*sf;             R(3,3) =  ct*cf
  end function rotmat

  ! J2: Euler-angle rates = J2 * body angular velocity
  pure function j2mat(e2) result(J)
    real(dp), intent(in) :: e2(3)
    real(dp) :: J(3,3), cf, sf, ct, st, tt
    cf = cos(e2(1));  sf = sin(e2(1));  ct = cos(e2(2));  st = sin(e2(2))
    tt = st/max(ct, 1.0e-12_dp)
    J(1,1) = 1.0_dp;  J(1,2) = sf*tt;  J(1,3) = cf*tt
    J(2,1) = 0.0_dp;  J(2,2) = cf;     J(2,3) = -sf
    J(3,1) = 0.0_dp;  J(3,2) = sf/max(ct, 1.0e-12_dp);  J(3,3) = cf/max(ct, 1.0e-12_dp)
  end function j2mat

  pure function j2inv(e2) result(J)
    real(dp), intent(in) :: e2(3)
    real(dp) :: J(3,3), cf, sf, ct, st
    cf = cos(e2(1));  sf = sin(e2(1));  ct = cos(e2(2));  st = sin(e2(2))
    J(1,1) = 1.0_dp;  J(1,2) = 0.0_dp;  J(1,3) = -st
    J(2,1) = 0.0_dp;  J(2,2) = cf;      J(2,3) = ct*sf
    J(3,1) = 0.0_dp;  J(3,2) = -sf;     J(3,3) = ct*cf
  end function j2inv

  ! rotation about the unit axis a by the angle th (Rodrigues)
  pure function rodrigues(a, th) result(R)
    real(dp), intent(in) :: a(3), th
    real(dp) :: R(3,3), K(3,3), c, s
    c = cos(th);  s = sin(th)
    K = 0.0_dp
    K(1,2) = -a(3);  K(1,3) =  a(2);  K(2,1) =  a(3);  K(2,3) = -a(1);  K(3,1) = -a(2);  K(3,2) =  a(1)
    R = eye3() + s*K + (1.0_dp - c)*matmul(K, K)
  end function rodrigues


  !---------------------------------------------------------------------
  ! PID controllers of the grid parts (update 21): resolve the owner and
  ! the target of every controller written in a Body / Child_body /
  ! Sub_child_body group; the first controller of a part rotates it
  ! (throttle: rate, steering: angle), its Motion_DOF is ignored.
  !---------------------------------------------------------------------
  subroutine pid_setup_parts
    integer :: i, k, n, lst(4), ierr
    real(dp) :: a(3)
    do k = 1, nbd
       call pid_find(1, bd(k)%bid, bd(k)%cid, bd(k)%sid, 0, n, lst)
       if (n == 0) cycle
       if (n > 1) then
          if (myrank == 0) write(*,'(A)') ' ERROR (6-DoF): '//part_name(k)//' has more than one PID controller'// &
             ' (one rotation axis per grid part)'
          call MPI_Abort(MPI_COMM_WORLD, 83, ierr)
       end if
       bd(k)%npid = n;  bd(k)%ipid(1:n) = lst(1:n)
       i = lst(1)
       call pid_check(i)
       pidc(i)%node = k
       pidc(i)%tnode = sixdof_node_of(pidc(i)%tb, pidc(i)%tc, pidc(i)%ts)
       if (pidc(i)%tnode == 0) then
          if (myrank == 0) write(*,'(A,3I3,A)') ' ERROR (6-DoF): the PID target of '//part_name(k)//' (Body / Child_body / '// &
             'Sub_child_body', pidc(i)%tb, pidc(i)%tc, pidc(i)%ts, ') is not a part of the csv'
          call MPI_Abort(MPI_COMM_WORLD, 83, ierr)
       end if
       a = pid_axis_unit(i)
       bd(k)%pid_axB = matmul(Fm, a)                       ! grid frame -> E-aligned mother / body frame
       bd(k)%pid_rot = .true.
       if (any(bd(k)%dof /= 0) .and. myrank == 0) write(*,'(A)') '  NOTE (6-DoF): Motion_DOF of '//part_name(k)// &
          ' is ignored - the part is rotated by its PID controller (Added_Mass / Motion_Damp unused)'
       bd(k)%dof = 0
       bd(k)%own_moving = .true.
       bd(k)%pid_ang = pidc(i)%angle;  bd(k)%pid_angd = 0.0_dp
       if (pidc(i)%typ == 1) bd(k)%pid_angd = 2.0_dp*pi*pidc(i)%out
       if (pidc(i)%typ == 2) bd(k)%pid_ang = pidc(i)%out*pi/180.0_dp
    end do
    do i = 1, npid                                           ! controllers of force models: target node
       if (pidc(i)%kind /= 2 .or. pidc(i)%use /= 1) cycle
       pidc(i)%tnode = sixdof_node_of(pidc(i)%tb, pidc(i)%tc, pidc(i)%ts)
    end do
  end subroutine pid_setup_parts

  ! measured value of controller i: eta / nu of its target part [m, deg, m/s, deg/s]
  real(dp) function pid_measure(i)
    integer, intent(in) :: i
    integer :: k, c
    real(dp) :: r2d
    r2d = 180.0_dp/pi
    k = pidc(i)%tnode;  c = pidc(i)%tvar
    if (pidc(i)%ttype == 1) then
       pid_measure = bd(k)%eta(c)
    else
       pid_measure = bd(k)%nu(c)
    end if
    if (c >= 4) pid_measure = pid_measure*r2d
  end function pid_measure

  ! one controller step of the PID-rotated part k: new angle / rate, state for the output files
  subroutine pid_rotate(k, t)
    integer,  intent(in) :: k
    real(dp), intent(in) :: t
    integer :: i
    real(dp) :: out, ang_old, R(3,3)
    i = bd(k)%ipid(1)
    out = pid_step(i, t, pid_measure(i))
    ang_old = bd(k)%pid_ang
    if (pidc(i)%typ == 1) then
       bd(k)%pid_angd = 2.0_dp*pi*out                       ! rev/s -> rad/s
       bd(k)%pid_ang  = ang_old + bd(k)%pid_angd*dt_dim
    else
       bd(k)%pid_ang  = out*pi/180.0_dp
       bd(k)%pid_angd = (bd(k)%pid_ang - ang_old)/max(dt_dim, 1.0e-300_dp)
    end if
    pidc(i)%angle = bd(k)%pid_ang
    ! state for the motion files (relative pose: Euler angles of the rotation, rates)
    R = rodrigues(bd(k)%pid_axB, bd(k)%pid_ang)
    bd(k)%eta(4:6) = euler_of(R)
    bd(k)%nu(4:6) = bd(k)%pid_axB*bd(k)%pid_angd
    bd(k)%etad(4:6) = matmul(j2mat(bd(k)%eta(4:6)), bd(k)%nu(4:6))
    bd(k)%etadd = 0.0_dp;  bd(k)%nud = 0.0_dp
  end subroutine pid_rotate

  ! solve A x = b (n <= 6) by Gaussian elimination with partial pivoting
  subroutine solve_small(n, A, b, x, ok)
    integer,  intent(in)  :: n
    real(dp), intent(in)  :: A(n,n), b(n)
    real(dp), intent(out) :: x(n)
    logical,  intent(out) :: ok
    real(dp) :: M(n,n+1), f, tmp(n+1)
    integer :: i, j, p
    ok = .true.
    M(:,1:n) = A;  M(:,n+1) = b
    do i = 1, n
       p = i
       do j = i+1, n
          if (abs(M(j,i)) > abs(M(p,i))) p = j
       end do
       if (abs(M(p,i)) < 1.0e-300_dp) then
          ok = .false.;  x = 0.0_dp;  return
       end if
       if (p /= i) then
          tmp = M(i,:);  M(i,:) = M(p,:);  M(p,:) = tmp
       end if
       do j = i+1, n
          f = M(j,i)/M(i,i)
          M(j,i:n+1) = M(j,i:n+1) - f*M(i,i:n+1)
       end do
    end do
    do i = n, 1, -1
       x(i) = (M(i,n+1) - sum(M(i,i+1:n)*x(i+1:n)))/M(i,i)
    end do
  end subroutine solve_small

  !=====================================================================
  ! Absolute pose / velocities / accelerations of part k from its state
  ! (absolute for a body, relative to the mother for a child / sub-
  ! child) by composition with the mother's ABSOLUTE quantities, which
  ! must be up to date (the nodes are ordered mothers first).  Sets
  ! R, pE, vE, OmE, aE, alB, the grid transform Rg, dg and the grid
  ! velocities vG, omG.
  !=====================================================================
  subroutine update_transform(k)
    integer, intent(in) :: k
    integer :: m
    real(dp) :: Rrel(3,3), r(3), om_m(3), vrel(3), om_rel(3), arel(3), Rm(3,3), al_m(3)
    m = bd(k)%par
    if (bd(k)%pid_rot) then                   ! PID-controlled rotation about the axis (update 21)
       Rrel = rodrigues(bd(k)%pid_axB, bd(k)%pid_ang)
       om_rel = bd(k)%pid_axB*bd(k)%pid_angd
    else
       Rrel = rotmat(bd(k)%eta(4:6))
       om_rel = matmul(j2inv(bd(k)%eta(4:6)), bd(k)%etad(4:6))
    end if
    if (m == 0) then
       bd(k)%R   = Rrel
       bd(k)%pE  = matmul(Fm, bd(k)%rot0) + bd(k)%eta(1:3)
       bd(k)%vE  = bd(k)%etad(1:3)
       bd(k)%OmE = matmul(Rrel, om_rel)
       bd(k)%aE  = bd(k)%etadd(1:3)
       bd(k)%alB = bd(k)%nud(4:6)
       bd(k)%nu(1:3) = matmul(transpose(Rrel), bd(k)%etad(1:3));  bd(k)%nu(4:6) = om_rel
    else
       Rm = bd(m)%R
       r = bd(k)%r0 + bd(k)%eta(1:3);  vrel = bd(k)%etad(1:3);  arel = bd(k)%etadd(1:3)
       om_m = matmul(transpose(Rm), bd(m)%OmE)                     ! mother's angular velocity in B_m
       al_m = bd(m)%alB
       bd(k)%R   = matmul(Rm, Rrel)
       bd(k)%pE  = bd(m)%pE + matmul(Rm, r)
       bd(k)%vE  = bd(m)%vE + matmul(Rm, cross3(om_m, r) + vrel)
       bd(k)%OmE = bd(m)%OmE + matmul(bd(k)%R, om_rel)
       bd(k)%aE  = bd(m)%aE + matmul(Rm, cross3(al_m, r) + cross3(om_m, cross3(om_m, r)) + 2.0_dp*cross3(om_m, vrel) + arel)
       bd(k)%alB = matmul(transpose(Rrel), al_m) - cross3(om_rel, matmul(transpose(Rrel), om_m)) + bd(k)%nud(4:6)
       bd(k)%nu(1:3) = matmul(transpose(Rrel), vrel);  bd(k)%nu(4:6) = om_rel
    end if
    bd(k)%Rg  = matmul(Fm, matmul(bd(k)%R, Fm))
    bd(k)%dg  = matmul(Fm, bd(k)%pE) - bd(k)%rot0
    bd(k)%vG  = matmul(Fm, bd(k)%vE)
    bd(k)%omG = matmul(Fm, bd(k)%OmE)
  end subroutine update_transform

  !---------------------------------------------------------------------
  ! Composite mass, CG and inertia about O_k in B_k of every part with
  ! its descendants at the current relative poses (bottom-up), and the
  ! rigid-body mass matrix MRBc (Fossen 2.91) of the composite.
  !---------------------------------------------------------------------
  subroutine composite_properties
    integer :: k, j, i, l
    real(dp) :: S(3), rj(3), Rjk(3,3), pjk(3,3), rr(3)
    do k = nbd, 1, -1
       bd(k)%Mc = bd(k)%mass
       S = bd(k)%mass*bd(k)%rcg
       rr = bd(k)%rcg
       do i = 1, 3
          do l = 1, 3
             bd(k)%IOc(i,l) = bd(k)%Ic(i,l) + bd(k)%mass*(merge(sum(rr*rr), 0.0_dp, i == l) - rr(i)*rr(l))
          end do
       end do
       do j = k+1, nbd
          if (.not. is_ancestor(k, j)) cycle
          Rjk = matmul(transpose(bd(k)%R), bd(j)%R)                        ! B_j -> B_k
          rj  = matmul(transpose(bd(k)%R), bd(j)%pE - bd(k)%pE) + matmul(Rjk, bd(j)%rcg)   ! CG_j - O_k in B_k
          bd(k)%Mc = bd(k)%Mc + bd(j)%mass
          S = S + bd(j)%mass*rj
          pjk = matmul(Rjk, matmul(bd(j)%Ic, transpose(Rjk)))
          do i = 1, 3
             do l = 1, 3
                bd(k)%IOc(i,l) = bd(k)%IOc(i,l) + pjk(i,l) + bd(j)%mass*(merge(sum(rj*rj), 0.0_dp, i == l) - rj(i)*rj(l))
             end do
          end do
       end do
       if (bd(k)%Mc > 0.0_dp) then
          bd(k)%rcgc = S/bd(k)%Mc
       else
          bd(k)%rcgc = bd(k)%rcg
       end if
       bd(k)%MRBc = 0.0_dp
       do i = 1, 3
          bd(k)%MRBc(i,i) = bd(k)%Mc
       end do
       bd(k)%MRBc(1:3,4:6) = -bd(k)%Mc*skew(bd(k)%rcgc)
       bd(k)%MRBc(4:6,1:3) =  bd(k)%Mc*skew(bd(k)%rcgc)
       bd(k)%MRBc(4:6,4:6) =  bd(k)%IOc
    end do
  end subroutine composite_properties

  !---------------------------------------------------------------------
  ! Subtree fluid wrench of every part in E about its O (bottom-up): its
  ! own wall forces + those of its children transferred to O_k, without
  ! the components along the children's free (released) joint DOFs -
  ! translations: the mother-frame axes of eta(1:3); rotations: the
  ! Euler-angle axes (columns of J2^-1 in B_child).
  !---------------------------------------------------------------------
  subroutine transfer_wrenches(t)
    real(dp), intent(in) :: t
    integer :: k, j, c
    real(dp) :: Ftr(3), Mtr(3), e(3), J2i(3,3)
    logical :: released
    released = (t >= Motion_Start_T - 1.0e-12_dp)
    do k = nbd, 1, -1
       bd(k)%Fs = matmul(Fm, bd(k)%fp + bd(k)%ff + bd(k)%fext)      ! fluid + force models (update 19)
       bd(k)%Ms = matmul(Fm, bd(k)%mp + bd(k)%mf + bd(k)%mext)
       do j = k+1, nbd
          if (bd(j)%par /= k) cycle
          Ftr = bd(j)%Fs;  Mtr = bd(j)%Ms
          if (released) then
             J2i = j2inv(bd(j)%eta(4:6))
             do c = 1, 6
                if (bd(j)%dof(c) /= 2) cycle
                if (c <= 3) then
                   e = bd(k)%R(:,c)                                      ! mother axis in E
                   Ftr = Ftr - sum(Ftr*e)*e
                else
                   e = matmul(bd(j)%R, J2i(:,c-3))
                   e = e/max(sqrt(sum(e*e)), 1.0e-12_dp)
                   Mtr = Mtr - sum(Mtr*e)*e
                end if
             end do
          end if
          bd(k)%Fs = bd(k)%Fs + Ftr
          bd(k)%Ms = bd(k)%Ms + Mtr + cross3(bd(j)%pE - bd(k)%pE, Ftr)
       end do
    end do
  end subroutine transfer_wrenches

  !=====================================================================
  ! One time step of every part: t_old -> t (forces of the previous
  ! step), mothers first, then the grid update and the overset
  ! re-assembly.
  !=====================================================================
  subroutine sixdof_advance(istep, t)
    integer,  intent(in) :: istep
    real(dp), intent(in) :: t
    integer :: k, ierr
    logical :: moved
    real(dp) :: state(30)
    if (.not. sixdof_on .or. .not. any_moving) return
    moved = .false.
    call composite_properties                 ! at the poses of the previous step
    call transfer_wrenches(t)
    do k = 1, nbd
       if (.not. bd(k)%moving) cycle
       if (bd(k)%pid_rot) then
          call pid_rotate(k, t)                ! identical on every rank (replicated states)
       else if (bd(k)%own_moving) then
          if (myrank == 0) call body_step(k, t)
          ! identical state on every rank (the grids are replicated)
          state(1:6) = bd(k)%eta;  state(7:12) = bd(k)%etad;  state(13:18) = bd(k)%etadd
          state(19:24) = bd(k)%nu;  state(25:30) = bd(k)%nud
          call MPI_Bcast(state, 30, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
          bd(k)%eta = state(1:6);  bd(k)%etad = state(7:12);  bd(k)%etadd = state(13:18)
          bd(k)%nu = state(19:24);  bd(k)%nud = state(25:30)
       end if
       bd(k)%Rg_prev = bd(k)%Rg;  bd(k)%dg_prev = bd(k)%dg
       call update_transform(k)               ! composition with the mother's new pose
       moved = .true.
    end do
    if (kbg > 0) then                         ! moving background (update 20): follows the body's new pose
       call bgr_follow(dt_dim, .false.)
       bd(kbg)%Rg_prev = bd(kbg)%Rg;  bd(kbg)%dg_prev = bd(kbg)%dg
       call update_transform(kbg)
       if (bgr_changed()) moved = .true.
       call bgr_publish
    end if
    if (moved) then
       call move_grids
       call after_motion
    end if
    if (istep < 0) ierr = 0
  end subroutine sixdof_advance

  !---------------------------------------------------------------------
  ! Rigid-body equations of part k over one step dt_dim ending at t
  ! (Newton-Euler about O_k in B_k for the acceleration RELATIVE to the
  ! mother frame; composite mass properties; the mother's motion of the
  ! current step enters through the frame accelerations).
  !---------------------------------------------------------------------
  subroutine body_step(k, t)
    integer,  intent(in) :: k
    real(dp), intent(in) :: t
    real(dp) :: Rrel(3,3), Rabs(3,3), J2(3,3), J2i(3,3), Ji(6,6), Meta(6,6), MA(6,6), A(6,6)
    real(dp) :: tau(6), bvec(6), bb(6), g(6), eps, J2p(3,3), om_rel(3), om_k(3), afr(3), alfr(3)
    real(dp) :: gE(3), W(3), ei(6), edi(6), eddi(6), acc(6), etad_n(6), eta_n(6), nu_n(6)
    real(dp) :: Af(6,6), bf(6), xf(6), mi(6), r(3), vrel(3), om_m(3), afr_m(3)
    integer :: i, j, nf, idf(6), c, m
    logical :: free(6), ok, released
    real(dp) :: dtt

    dtt = dt_dim
    m = bd(k)%par
    released = (t >= Motion_Start_T - 1.0e-12_dp)
    free = (bd(k)%dof == 2) .and. released
    ! imposed components at the new time
    ei = bd(k)%eta;  edi = bd(k)%etad;  eddi = bd(k)%etadd
    if (any(bd(k)%dof == 1)) call imposed_at(k, t, ei, edi, eddi)

    acc = 0.0_dp
    if (any(free)) then
       Rrel = rotmat(bd(k)%eta(4:6));  J2 = j2mat(bd(k)%eta(4:6));  J2i = j2inv(bd(k)%eta(4:6))
       Ji = 0.0_dp
       Ji(1:3,1:3) = transpose(Rrel);  Ji(4:6,4:6) = J2i
       om_rel = matmul(J2i, bd(k)%etad(4:6))
       if (m == 0) then
          Rabs = Rrel;  om_k = om_rel;  afr = 0.0_dp;  alfr = 0.0_dp
       else
          ! frame terms from the mother's absolute motion (already advanced to t)
          Rabs = matmul(bd(m)%R, Rrel)
          om_m = matmul(transpose(bd(m)%R), bd(m)%OmE)
          r = bd(k)%r0 + bd(k)%eta(1:3);  vrel = bd(k)%etad(1:3)
          afr_m = matmul(transpose(bd(m)%R), bd(m)%aE) + cross3(bd(m)%alB, r) + cross3(om_m, cross3(om_m, r)) &
                  + 2.0_dp*cross3(om_m, vrel)
          afr  = matmul(transpose(Rrel), afr_m)
          alfr = matmul(transpose(Rrel), bd(m)%alB) - cross3(om_rel, matmul(transpose(Rrel), om_m))
          om_k = matmul(transpose(Rrel), om_m) + om_rel
       end if
       ! ---- generalized force in B_k about O_k ------------------------------
       tau(1:3) = matmul(transpose(Rabs), bd(k)%Fs);  tau(4:6) = matmul(transpose(Rabs), bd(k)%Ms)
       gE = matmul(Fm, grav)                                   ! (0, 0, +g) in E
       W  = matmul(transpose(Rabs), bd(k)%Mc*gE)               ! weight of the composite in B_k
       tau(1:3) = tau(1:3) + W
       tau(4:6) = tau(4:6) + cross3(bd(k)%rcgc, W)
       mi(1:3) = bd(k)%Mc;  mi(4) = bd(k)%IOc(1,1);  mi(5) = bd(k)%IOc(2,2);  mi(6) = bd(k)%IOc(3,3)
       do c = 1, 6
          tau(c) = tau(c) - bd(k)%damp(c)*mi(c)*bd(k)%nu(c)
       end do
       ! ---- Newton-Euler right-hand side for x = [Rrel^T a_rel ; omegadot_rel] ---
       bb(1:3) = tau(1:3) - bd(k)%Mc*(afr + cross3(alfr, bd(k)%rcgc) + cross3(om_k, cross3(om_k, bd(k)%rcgc)))
       bb(4:6) = tau(4:6) - cross3(om_k, matmul(bd(k)%IOc, om_k)) - matmul(bd(k)%IOc, alfr) - bd(k)%Mc*cross3(bd(k)%rcgc, afr)
       ! ---- eta representation: x = Ji etadd - g, g = [0 ; J2^-1 J2dot omega_rel] ---
       eps = 1.0e-6_dp
       J2p = j2mat(bd(k)%eta(4:6) + eps*bd(k)%etad(4:6))
       g = 0.0_dp
       g(4:6) = matmul(J2i, matmul((J2p - J2)/eps, om_rel))
       Meta = matmul(transpose(Ji), matmul(bd(k)%MRBc, Ji))
       bvec = matmul(transpose(Ji), bb + matmul(bd(k)%MRBc, g))
       MA = 0.0_dp
       do c = 1, 6
          MA(c,c) = bd(k)%madd(c)
       end do
       MA = matmul(transpose(Ji), matmul(MA, Ji))
       A = Meta + MA
       bvec = bvec + matmul(MA, bd(k)%etadd)
       ! known accelerations: fixed -> 0, imposed -> table, held free -> 0
       do c = 1, 6
          if (bd(k)%dof(c) == 1) acc(c) = eddi(c)
       end do
       nf = 0
       do c = 1, 6
          if (free(c)) then
             nf = nf + 1;  idf(nf) = c
          end if
       end do
       do i = 1, nf
          bf(i) = bvec(idf(i))
          do c = 1, 6
             if (.not. free(c)) bf(i) = bf(i) - A(idf(i),c)*acc(c)
          end do
          do j = 1, nf
             Af(i,j) = A(idf(i),idf(j))
          end do
       end do
       call solve_small(nf, Af(1:nf,1:nf), bf(1:nf), xf(1:nf), ok)
       if (.not. ok) then
          write(*,'(A)') ' WARNING (6-DoF): singular mass matrix of '//part_name(k)//' (pitch = 90 deg?), accelerations set to 0'
          xf = 0.0_dp
       end if
       do i = 1, nf
          acc(idf(i)) = xf(i)
       end do
    else
       do c = 1, 6
          if (bd(k)%dof(c) == 1) acc(c) = eddi(c)
       end do
    end if
    ! ---- integration (velocity Verlet on eta) ----------------------------
    ! The forces of this step were evaluated at the pose reached by the previous step, so acc
    ! is a^n at eta^n.  The velocity of the previous free step was advanced with a^(n-1) only;
    ! its second half-kick 1/2 dt (a^n - a^(n-1)) is applied now (update 19: without it the
    ! scheme was first order in energy - visible as a 5 % amplitude error per half period of a
    ! stiff point-force spring resolved with 50 steps per period; with the usual thousands of
    ! steps per motion period of a fluid-driven body the difference is negligible).
    do c = 1, 6
       if (free(c) .and. bd(k)%vhalf(c)) bd(k)%etad(c) = bd(k)%etad(c) + 0.5_dp*dtt*(acc(c) - bd(k)%etadd(c))
    end do
    eta_n = bd(k)%eta;  etad_n = bd(k)%etad
    do c = 1, 6
       bd(k)%vhalf(c) = .false.
       select case (bd(k)%dof(c))
       case (1)
          eta_n(c) = ei(c);  etad_n(c) = edi(c)
       case (2)
          if (free(c)) then
             eta_n(c)  = bd(k)%eta(c) + dtt*bd(k)%etad(c) + 0.5_dp*dtt*dtt*acc(c)
             etad_n(c) = bd(k)%etad(c) + dtt*acc(c)
             bd(k)%vhalf(c) = .true.
          else
             etad_n(c) = 0.0_dp
          end if
       case default
          etad_n(c) = 0.0_dp
       end select
    end do
    nu_n = bd(k)%nu
    bd(k)%eta = eta_n;  bd(k)%etad = etad_n;  bd(k)%etadd = acc
    Rrel = rotmat(bd(k)%eta(4:6))
    bd(k)%nu(1:3) = matmul(transpose(Rrel), bd(k)%etad(1:3))
    bd(k)%nu(4:6) = matmul(j2inv(bd(k)%eta(4:6)), bd(k)%etad(4:6))
    bd(k)%nud = (bd(k)%nu - nu_n)/dtt
  end subroutine body_step

  !=====================================================================
  ! Grid update of the moving blocks from the original coordinates.
  !=====================================================================
  subroutine move_grids
    integer :: b, k, ib, d, i, j, kk, c, h1, h2, h3, q
    real(dp) :: Rg(3,3), c0(3), dg(3), c0n(3), dgn(3), x(3), xf(3), ug(3), v(3)
    integer  :: hx, hy, hz
    do b = 1, nblk
       k = blk_obj(b)
       if (k == 0) cycle
       if (.not. bd(k)%moving) cycle
       Rg = bd(k)%Rg;  c0 = bd(k)%rot0;  dg = bd(k)%dg
       do kk = 1, blk(b)%nk
        do j = 1, blk(b)%nj
         do i = 1, blk(b)%ni
            x = matmul(Rg, [g0(b)%x(i,j,kk), g0(b)%y(i,j,kk), g0(b)%z(i,j,kk)] - c0) + c0 + dg
            blk(b)%xn(i,j,kk) = x(1);  blk(b)%yn(i,j,kk) = x(2);  blk(b)%zn(i,j,kk) = x(3)
         end do
        end do
       end do
       c0n = c0/Lch;  dgn = dg/Lch
       do kk = lbound(blk(b)%xg,4), ubound(blk(b)%xg,4)
        do j = lbound(blk(b)%xg,3), ubound(blk(b)%xg,3)
         do i = lbound(blk(b)%xg,2), ubound(blk(b)%xg,2)
            blk(b)%xg(:,i,j,kk) = matmul(Rg, g0(b)%xg(:,i,j,kk) - c0n) + c0n + dgn
         end do
        end do
       end do
    end do
    ! owned sub-domains: half-index coordinates, rotated metric vectors, grid fluxes
    do ib = 1, nown
       d = own_list(ib);  b = dom(d)%b;  k = blk_obj(b)
       if (.not. dom(d)%moving) cycle
       ! morphing blocks (update 24c): fsi_morph builds their coordinates, metrics and grid fluxes;
       ! the rigid rotation of the UNDEFORMED Xh0 / mh0 must not overwrite the level-n metrics here,
       ! which fsi_morph keeps as the old metrics of the flux re-projection (with the overwrite the
       ! re-projection added fvw (m_deformed - m_undeformed) . u to the fluxes in every step)
       if (fsi_on .and. allocated(fsi_blk)) then
          if (fsi_blk(b)) cycle
       end if
       call set_block(d)
       Rg = bd(k)%Rg;  c0n = bd(k)%rot0/Lch;  dgn = bd(k)%dg/Lch
       do h3 = hlo(3), hhi(3)
        do h2 = hlo(2), hhi(2)
         do h1 = hlo(1), hhi(1)
            Xh(:,h1,h2,h3) = matmul(Rg, dom(d)%Xh0(:,h1,h2,h3) - c0n) + c0n + dgn
            do q = 1, 3
               mh(3*q-2:3*q, h1,h2,h3) = matmul(Rg, dom(d)%mh0(3*q-2:3*q, h1,h2,h3))
            end do
         end do
        end do
       end do
       ! grid velocity (grid frame) = displacement of the face over the step / dt (update 20:
       ! the discrete grid velocity, consistent with the moved coordinates the ALE time
       ! derivative sees; it equals v_O + omega x (x - x_O) for a uniform motion and the mean
       ! velocity over the step for an accelerating one - with v^(n+1) a body accelerating at a
       ! convected its near field with the spurious velocity a dt / 2)
       do kk = lbound(U1,3), ubound(U1,3)
        do j = lbound(U1,2), ubound(U1,2)
         do i = lbound(U1,1), ubound(U1,1)
            do q = 1, 3
               hx = 2*i;  hy = 2*j;  hz = 2*kk
               select case (q)
               case (1);  hx = hx + 1
               case (2);  hy = hy + 1
               case default;  hz = hz + 1
               end select
               hx = min(max(hx, hlo(1)+1), hhi(1)-1)
               hy = min(max(hy, hlo(2)+1), hhi(2)-1)
               hz = min(max(hz, hlo(3)+1), hhi(3)-1)
               x  = grid_velocity(k, Xh(:,hx,hy,hz))                     ! non-dimensional grid velocity
               select case (q)
               case (1);  Ug1(i,j,kk) = fvw(hx,hy,hz)*(mq(1,1,hx,hy,hz)*x(1) + mq(2,1,hx,hy,hz)*x(2) + mq(3,1,hx,hy,hz)*x(3))
               case (2);  Ug2(i,j,kk) = fvw(hx,hy,hz)*(mq(1,2,hx,hy,hz)*x(1) + mq(2,2,hx,hy,hz)*x(2) + mq(3,2,hx,hy,hz)*x(3))
               case default;  Ug3(i,j,kk) = fvw(hx,hy,hz)*(mq(1,3,hx,hy,hz)*x(1) + mq(2,3,hx,hy,hz)*x(2) + mq(3,3,hx,hy,hz)*x(3))
               end select
            end do
         end do
        end do
       end do
    end do
    ! flexible parts (update 24): the morphing of their blocks, metrics, grid velocities
    if (fsi_on .and. associated(fsi_morph_fn)) call fsi_morph_fn(restoring)
    c = 0
  end subroutine move_grids

  !---------------------------------------------------------------------
  ! Cartesian grid velocity (non-dimensional, grid frame) of block b at
  ! the non-dimensional position xnd: the wall velocity of the no-slip
  ! patches of moving parts (boundary.F90) and the wall-function
  ! relative velocity (turbulence.F90).  Zero on static blocks.
  !---------------------------------------------------------------------
  function sixdof_wall_velocity(b, xnd, h) result(ug)
    integer,  intent(in) :: b
    real(dp), intent(in) :: xnd(3)
    integer,  intent(in), optional :: h(3)
    real(dp) :: ug(3)
    integer :: k
    ug = 0.0_dp
    if (.not. sixdof_on) return
    if (.not. allocated(blk_obj)) return
    k = blk_obj(b)
    if (k == 0) return
    if (k == kbg) return                      ! background (update 20): Dirichlet values are earth-frame velocities
    if (.not. bd(k)%moving) return
    if (fsi_on .and. present(h) .and. allocated(fsi_blk)) then
       if (fsi_blk(b) .and. associated(fsi_hvel_fn)) then   ! morphing block (update 24): discrete velocity of the position
          ug = fsi_hvel_fn(b, h);  return
       end if
    end if
    ug = grid_velocity(k, xnd)
  end function sixdof_wall_velocity

  !---------------------------------------------------------------------
  ! Discrete grid velocity (non-dimensional, grid frame) of part k at the
  ! non-dimensional CURRENT position xnd: the displacement of that grid
  ! point over the step divided by the time step,
  !   x^n = Rg_prev (x0 - O0) + O0 + dg_prev,  x0 = Rg^T (x - O0 - dg) + O0,
  ! so that the relative flux U - Ug of a face moving with the grid is
  ! exactly zero and the ALE convection is consistent with the moved
  ! coordinates (update 20).  Equal to v_O + omega x (x - x_O) for a
  ! uniform translation / rotation to O(dt).
  !---------------------------------------------------------------------
  function grid_velocity(k, xnd) result(ug)
    integer,  intent(in) :: k
    real(dp), intent(in) :: xnd(3)
    real(dp) :: ug(3), x(3), x0(3), xp(3)
    x  = xnd*Lch
    x0 = matmul(transpose(bd(k)%Rg), x - bd(k)%rot0 - bd(k)%dg) + bd(k)%rot0
    xp = matmul(bd(k)%Rg_prev, x0 - bd(k)%rot0) + bd(k)%rot0 + bd(k)%dg_prev
    ug = (x - xp)/max(dt_dim, 1.0e-300_dp)/Uch
  end function grid_velocity

  !=====================================================================
  ! After the grids moved: overset re-assembly, fill of the uncovered
  ! cells, transformed components from the (Lagrangian) Cartesian
  ! field with the new metrics, fringe update.
  !=====================================================================
  subroutine after_motion(sub)
    logical, intent(in), optional :: sub       ! FSI sub-iteration (update 24): the rigid poses are those of the step start
    integer :: ib, d, i, j, k, ierr
    logical :: subit
    type :: st_t
       integer, allocatable :: c(:,:,:), f1(:,:,:), f2(:,:,:), f3(:,:,:)
    end type st_t
    type(st_t), allocatable :: old(:)
    real(dp) :: tsec
    subit = .false.
    if (present(sub)) subit = sub
    if (ovst_on .and. reassembly_needed()) then
       allocate(old(nown))
       do ib = 1, nown
          d = own_list(ib)
          allocate(old(ib)%c, source=dom(d)%cs);  allocate(old(ib)%f1, source=dom(d)%fs1)
          allocate(old(ib)%f2, source=dom(d)%fs2);  allocate(old(ib)%f3, source=dom(d)%fs3)
       end do
       call ovst_reassemble(tsec)
       t_reassemble = t_reassemble + tsec;  n_reassemble = n_reassemble + 1
       do ib = 1, nown
          d = own_list(ib)
          call set_block(d)
          if (.not. restoring) call fill_uncovered(old(ib)%c, old(ib)%f1, old(ib)%f2, old(ib)%f3)
          deallocate(old(ib)%c, old(ib)%f1, old(ib)%f2, old(ib)%f3)
       end do
       deallocate(old)
    end if
    ! transformed components on the rotated faces: the Cartesian velocity
    ! of a face is carried with the grid (Lagrangian), so the face
    ! component changes by fvw (m_new - m_old) . u_vec; only that
    ! INCREMENT uses the interpolated foreign Cartesian components (a
    ! translation leaves the stored fluxes exactly as they are)
    do ib = 1, nown
       d = own_list(ib)
       if (.not. dom(d)%moving) cycle
       if (restoring) cycle                     ! restart: the fluxes were saved on the moved grid
       k = blk_obj(dom(d)%b)
       if (fsi_on .and. allocated(fsi_blk)) then
          if (fsi_blk(dom(d)%b)) then                ! morphed block (update 24): re-projection with the new metrics
             if (associated(fsi_flux_fn)) call fsi_flux_fn(d)
             cycle
          end if
       end if
       ! rigid blocks: rotated from the level-n fluxes - at the step start, and again in every FSI
       ! sub-iteration, which restores the level-n fluxes before the re-mesh (24b)
       call rotate_fluxes(d, bd(k)%Rg - bd(k)%Rg_prev)
    end do
    call refresh_uvw(1)
    call contra_to_cart(.false., .true.)
    if (ovst_on) then
       call ovs_update_faces(1)
       call refresh_uvw(1)
       call contra_to_cart(.false., .true.)
       call ovs_update_cells(8)
       call ovs_finish_status
    end if
    call refresh_all(4)
    if (Num_of_Phase == 2) call refresh_all(8)
    if (cav_on) call refresh_all(19)
    if (energy_on) call refresh_all(20)
    if (turb_on) then
       call refresh_all(14);  call refresh_all(16)
       if (TURB_model >= 2) call refresh_all(15)
    end if
    ! star field = current field (the predictor starts from it)
    do ib = 1, nown
       d = own_list(ib)
       dom(d)%us = dom(d)%u;  dom(d)%vs = dom(d)%v;  dom(d)%ws = dom(d)%w
       dom(d)%Us1 = dom(d)%U1;  dom(d)%Us2 = dom(d)%U2;  dom(d)%Us3 = dom(d)%U3
    end do
    i = 0;  j = 0;  k = 0;  ierr = 0
  end subroutine after_motion

  !---------------------------------------------------------------------
  ! Re-assembly tolerance (update 20): the overset is re-assembled when a
  ! moving part has moved more than Ovst_Reassemble_Tol [m] RELATIVE to
  ! the background (relative to the fixed grid when the background does
  ! not move) since the last assembly - measured at its rotation centre
  ! and at a point 1 m away along each axis (rotation).  0 = every step.
  !---------------------------------------------------------------------
  logical function reassembly_needed()
    integer :: k, c
    real(dp) :: dmax, Rb(3,3), db(3), Ob(3), p0(3), p1(3), q0(3), q1(3), e(3)
    reassembly_needed = .true.
    if (ovst_reasm_tol <= 0.0_dp .or. restoring .or. fsi_reasm_request) then
       call remember_asm_pose;  fsi_reasm_request = .false.;  return
    end if
    if (.not. asm_valid) then
       call remember_asm_pose;  return
    end if
    Rb = eye3();  db = 0.0_dp;  Ob = 0.0_dp
    if (kbg > 0) then
       Rb = bd(kbg)%Rg;  db = bd(kbg)%dg;  Ob = bd(kbg)%rot0
    end if
    dmax = 0.0_dp
    do k = 1, nbd
       if (.not. bd(k)%moving) cycle
       ! the part's rotation centre and unit offsets, now and at the last assembly, in the background frame
       do c = 0, 3
          e = 0.0_dp
          if (c > 0) e(c) = 1.0_dp
          p1 = matmul(bd(k)%Rg, e) + bd(k)%rot0 + bd(k)%dg                    ! current grid position
          p0 = matmul(asm_Rg(:,:,k), e) + bd(k)%rot0 + asm_dg(:,k)            ! at the last assembly
          q1 = matmul(transpose(Rb), p1 - Ob - db)                            ! in the background frame
          if (kbg > 0) then
             q0 = matmul(transpose(asm_Rg(:,:,kbg)), p0 - Ob - asm_dg(:,kbg))
          else
             q0 = p0
          end if
          dmax = max(dmax, sqrt(sum((q1 - q0)**2)))
       end do
    end do
    if (dmax < ovst_reasm_tol) then
       reassembly_needed = .false.;  n_skipped = n_skipped + 1
    else
       call remember_asm_pose
    end if
  end function reassembly_needed

  subroutine remember_asm_pose
    integer :: k, kmax
    kmax = max(nbd, kbg)
    if (.not. allocated(asm_Rg)) allocate(asm_Rg(3,3,MAXB), asm_dg(3,MAXB))
    do k = 1, kmax
       asm_Rg(:,:,k) = bd(k)%Rg;  asm_dg(:,k) = bd(k)%dg
    end do
    asm_valid = .true.
  end subroutine remember_asm_pose

  !---------------------------------------------------------------------
  ! U^q <- U^q + fvw * ((dRg m0_q) . u_vec) on every face of domain d
  ! (dRg = change of the grid rotation matrix over the step; the metric
  ! vectors are m = Rg m0).  The foreign Cartesian components at the
  ! face are the same 4-point averages as in cart_to_contra.
  !---------------------------------------------------------------------
  subroutine rotate_fluxes(d, dRg)
    integer,  intent(in) :: d
    real(dp), intent(in) :: dRg(3,3)
    integer :: i, j, k, lo(3), hi(3), h1, h2, h3
    real(dp) :: uc(3), dm(3)
    if (all(abs(dRg) < 1.0e-300_dp)) return
    call set_block(d)
    lo = [lbound(U1,1), lbound(U1,2), lbound(U1,3)];  hi = [ubound(U1,1), ubound(U1,2), ubound(U1,3)]
    ! family 1
    do k = lo(3)+1, hi(3)-1
     do j = lo(2)+1, hi(2)-1
      do i = lo(1), hi(1)-1
         h1 = 2*i+1;  h2 = 2*j;  h3 = 2*k
         if (h1 <= hlo(1) .or. h1 >= hhi(1)) cycle
         uc(1) = u(i,j,k)
         uc(2) = 0.25_dp*(v(i,j-1,k) + v(i+1,j-1,k) + v(i,j,k) + v(i+1,j,k))
         uc(3) = 0.25_dp*(w(i,j,k-1) + w(i+1,j,k-1) + w(i,j,k) + w(i+1,j,k))
         dm = matmul(dRg, dom(d)%mh0(1:3,h1,h2,h3))
         U1(i,j,k) = U1(i,j,k) + fvw(h1,h2,h3)*sum(dm*uc)
      end do
     end do
    end do
    ! family 2
    do k = lo(3)+1, hi(3)-1
     do j = lo(2), hi(2)-1
      do i = lo(1)+1, hi(1)-1
         h1 = 2*i;  h2 = 2*j+1;  h3 = 2*k
         if (h2 <= hlo(2) .or. h2 >= hhi(2)) cycle
         uc(2) = v(i,j,k)
         uc(1) = 0.25_dp*(u(i-1,j,k) + u(i-1,j+1,k) + u(i,j,k) + u(i,j+1,k))
         uc(3) = 0.25_dp*(w(i,j,k-1) + w(i,j+1,k-1) + w(i,j,k) + w(i,j+1,k))
         dm = matmul(dRg, dom(d)%mh0(4:6,h1,h2,h3))
         U2(i,j,k) = U2(i,j,k) + fvw(h1,h2,h3)*sum(dm*uc)
      end do
     end do
    end do
    ! family 3
    do k = lo(3), hi(3)-1
     do j = lo(2)+1, hi(2)-1
      do i = lo(1)+1, hi(1)-1
         h1 = 2*i;  h2 = 2*j;  h3 = 2*k+1
         if (h3 <= hlo(3) .or. h3 >= hhi(3)) cycle
         uc(3) = w(i,j,k)
         uc(1) = 0.25_dp*(u(i-1,j,k) + u(i-1,j,k+1) + u(i,j,k) + u(i,j,k+1))
         uc(2) = 0.25_dp*(v(i,j-1,k) + v(i,j-1,k+1) + v(i,j,k) + v(i,j,k+1))
         dm = matmul(dRg, dom(d)%mh0(7:9,h1,h2,h3))
         U3(i,j,k) = U3(i,j,k) + fvw(h1,h2,h3)*sum(dm*uc)
      end do
     end do
    end do
  end subroutine rotate_fluxes

  !---------------------------------------------------------------------
  ! Cells / faces of the current domain that were holes and are not any
  ! more take the value of the nearest previously valid position
  ! (Chebyshev index search up to 3 layers) - a first estimate that the
  ! fringe interpolation / the solver then overwrites.
  !---------------------------------------------------------------------
  subroutine fill_uncovered(oc, of1, of2, of3)
    integer, intent(in) :: oc(is-ng:,js-ng:,ks-ng:), of1(is-ng:,js-ng:,ks-ng:), of2(is-ng:,js-ng:,ks-ng:), of3(is-ng:,js-ng:,ks-ng:)
    integer :: i, j, k, r, ii, jj, kk, ib, jb, kb, lo(3), hi(3)
    real(dp) :: dbest, dd
    logical :: found
    lo = [is-ng, js-ng, ks-ng];  hi = [ie+ng, je+ng, ke+ng]
    ! cells
    do k = lo(3), hi(3)
     do j = lo(2), hi(2)
      do i = lo(1), hi(1)
         if (.not. (oc(i,j,k) < 0 .and. ocs(i,j,k) >= 0)) cycle
         found = .false.;  dbest = huge(1.0_dp)
         do r = 1, 3
            do kk = max(k-r, lo(3)), min(k+r, hi(3))
             do jj = max(j-r, lo(2)), min(j+r, hi(2))
              do ii = max(i-r, lo(1)), min(i+r, hi(1))
                 if (oc(ii,jj,kk) < 0) cycle
                 dd = real((ii-i)**2 + (jj-j)**2 + (kk-k)**2, dp)
                 if (dd < dbest) then
                    dbest = dd;  ib = ii;  jb = jj;  kb = kk;  found = .true.
                 end if
              end do
             end do
            end do
            if (found) exit
         end do
         if (.not. found) cycle
         p(i,j,k) = p(ib,jb,kb);  alpha(i,j,k) = alpha(ib,jb,kb);  alpn(i,j,k) = alpn(ib,jb,kb)
         if (cav_on) then
            vap(i,j,k) = vap(ib,jb,kb);  vapn(i,j,k) = vapn(ib,jb,kb)
         end if
         if (energy_on) then
            temp(i,j,k) = temp(ib,jb,kb);  tempn(i,j,k) = tempn(ib,jb,kb)
         end if
         if (allocated(dom(cur_dom)%tk)) then
            tk(i,j,k) = tk(ib,jb,kb);  tw(i,j,k) = tw(ib,jb,kb);  mut(i,j,k) = mut(ib,jb,kb)
         end if
      end do
     end do
    end do
    ! faces (Cartesian field; the transformed one is rebuilt afterwards)
    call fill_faces(1, of1, u, un)
    call fill_faces(2, of2, v, vn)
    call fill_faces(3, of3, w, wn)
  contains
    subroutine fill_faces(q, ofs, f, fn)
      integer, intent(in) :: q
      integer, intent(in) :: ofs(is-ng:,js-ng:,ks-ng:)
      real(dp), intent(inout) :: f(is-ng:,js-ng:,ks-ng:), fn(is-ng:,js-ng:,ks-ng:)
      integer :: nfs(3), e(3)
      integer, pointer :: nfs_ptr(:,:,:)
      e = 0;  e(q) = 1
      select case (q)
      case (1);  nfs_ptr => ofs1
      case (2);  nfs_ptr => ofs2
      case default;  nfs_ptr => ofs3
      end select
      nfs = hi - e
      do k = lo(3), nfs(3)
       do j = lo(2), nfs(2)
        do i = lo(1), nfs(1)
           if (.not. (ofs(i,j,k) < 0 .and. nfs_ptr(i,j,k) >= 0)) cycle
           found = .false.;  dbest = huge(1.0_dp)
           do r = 1, 3
              do kk = max(k-r, lo(3)), min(k+r, nfs(3))
               do jj = max(j-r, lo(2)), min(j+r, nfs(2))
                do ii = max(i-r, lo(1)), min(i+r, nfs(1))
                   if (ofs(ii,jj,kk) < 0) cycle
                   dd = real((ii-i)**2 + (jj-j)**2 + (kk-k)**2, dp)
                   if (dd < dbest) then
                      dbest = dd;  ib = ii;  jb = jj;  kb = kk;  found = .true.
                   end if
                end do
               end do
              end do
              if (found) exit
           end do
           if (.not. found) cycle
           f(i,j,k) = f(ib,jb,kb);  fn(i,j,k) = fn(ib,jb,kb)
        end do
       end do
      end do
    end subroutine fill_faces
  end subroutine fill_uncovered

  !=====================================================================
  ! Fluid forces and moments on every body from the current field.
  !=====================================================================
  subroutine sixdof_forces
    integer :: ib, d, b, n, k, ax, sd, t1, t2, rr(2,3), m1, m2, n1, n2, ierr
    integer :: c1s, c1e, c2s, c2e
    real(dp), allocatable :: wgt(:,:)
    real(dp) :: acc(13, MAXB), Av(3), xf(3), amag, nrm(3), pw, tauw(3), fpv(3), ffv(3), r(3), wt
    logical :: ok
    if (.not. sixdof_on) return
    acc = 0.0_dp
    do ib = 1, nown
       d = own_list(ib);  b = dom(d)%b;  k = blk_obj(b)
       if (k == 0 .or. k == kbg) cycle            ! (the background's tank walls carry no body)
       call set_block(d)
       do n = 1, nbc
          if (.not. is_wall_patch(b, n)) cycle
          if (.not. wall_patch_range(d, b, n, ax, sd, t1, t2, rr, c1s, c1e, c2s, c2e)) cycle
          n1 = rr(2,t1) - rr(1,t1);  n2 = rr(2,t2) - rr(1,t2)          ! face cells
          allocate(wgt(n1, n2))
          call ovs_wall_weights(b, n, wgt)
          do m2 = c2s, c2e
           do m1 = c1s, c1e
              wt = wgt(m1 - rr(1,t1) + 1, m2 - rr(1,t2) + 1)
              if (wt <= 0.0_dp) cycle
              call wall_face_stress(b, n, ax, sd, t1, t2, m1, m2, ok, Av, amag, nrm, xf, pw, tauw)
              if (.not. ok) cycle
              fpv = -pw*Av
              ffv = tauw*amag
              r = xf - (bd(k)%rot0 + bd(k)%dg)
              acc(1:3,k)   = acc(1:3,k)   + wt*fpv
              acc(4:6,k)   = acc(4:6,k)   + wt*ffv
              acc(7:9,k)   = acc(7:9,k)   + wt*cross3(r, fpv)
              acc(10:12,k) = acc(10:12,k) + wt*cross3(r, ffv)
              acc(13,k)    = acc(13,k)    + wt*amag
           end do
          end do
          deallocate(wgt)
       end do
    end do
    call MPI_Allreduce(MPI_IN_PLACE, acc, 13*MAXB, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    do k = 1, nbd
       bd(k)%fp = acc(1:3,k);  bd(k)%ff = acc(4:6,k);  bd(k)%mp = acc(7:9,k);  bd(k)%mf = acc(10:12,k)
       bd(k)%area = acc(13,k)
    end do
  end subroutine sixdof_forces

  !---------------------------------------------------------------------
  ! the part of wall patch n of block b that sub-domain d holds (update 22: shared by the
  ! force integration and the surface output): .false. when the domain does not touch the
  ! patch.  ax / sd = face direction / side, t1 t2 the tangential index directions, rr the node
  ! ranges of the patch, c1s..c2e the face-cell ranges of this domain
  !---------------------------------------------------------------------
  logical function wall_patch_range(d, b, n, ax, sd, t1, t2, rr, c1s, c1e, c2s, c2e) result(has)
    integer, intent(in)  :: d, b, n
    integer, intent(out) :: ax, sd, t1, t2, rr(2,3), c1s, c1e, c2s, c2e
    has = .false.
    ax = abs(bc_dir(n));  sd = merge(1, 2, bc_dir(n) > 0)
    rr(:,1) = bc_irng(:,n);  rr(:,2) = bc_jrng(:,n);  rr(:,3) = bc_krng(:,n)
    call tang_axes(ax, t1, t2)
    if (rr(2,t1) - rr(1,t1) < 1 .or. rr(2,t2) - rr(1,t2) < 1) return
    if (sd == 1) then
       if (dom(d)%r1(ax) /= 1) return
    else
       if (dom(d)%r2(ax) /= ncell(b, ax)) return
    end if
    c1s = max(rr(1,t1), dom(d)%r1(t1));  c1e = min(rr(2,t1)-1, dom(d)%r2(t1))
    c2s = max(rr(1,t2), dom(d)%r1(t2));  c2e = min(rr(2,t2)-1, dom(d)%r2(t2))
    has = (c1e >= c1s .and. c2e >= c2s)
  end function wall_patch_range

  !---------------------------------------------------------------------
  ! pressure and friction on ONE wall face (face cell m1, m2 of patch n of block b, the current
  ! block set by set_block; update 22: factored out of sixdof_forces, used unchanged by the force
  ! integration and by the solid-wall surface output).  Av = area vector into the fluid [m2],
  ! amag its magnitude, nrm the unit normal, xf the face centre [m], pw the wall pressure [Pa]
  ! (first cell, zero gradient), tauw the friction stress vector ON THE WALL [Pa] (the first-cell
  ! velocity relative to the wall, tangential part, times mu_eff / dn incl. the wall function).
  ! ok = .false. for a hole cell or a degenerate face.  The pressure force on the wall is -pw Av,
  ! the friction force tauw amag.
  !---------------------------------------------------------------------
  subroutine wall_face_stress(b, n, ax, sd, t1, t2, m1, m2, ok, Av, amag, nrm, xf, pw, tauw)
    integer,  intent(in)  :: b, n, ax, sd, t1, t2, m1, m2
    logical,  intent(out) :: ok
    real(dp), intent(out) :: Av(3), amag, nrm(3), xf(3), pw, tauw(3)
    integer  :: idx(3), i, j, kk, nn(4,3), hx, hy, hz
    real(dp) :: x1(3), x2(3), x3(3), x4(3), xc(3), uc(3), uw(3), urel(3), ut(3), dn, mu_eff, xnd(3)
    ok = .false.;  Av = 0.0_dp;  amag = 0.0_dp;  nrm = 0.0_dp;  xf = 0.0_dp;  pw = 0.0_dp;  tauw = 0.0_dp
    idx(t1) = m1;  idx(t2) = m2
    idx(ax) = merge(1, ncell(b, ax), sd == 1)                 ! adjacent cell
    i = idx(1);  j = idx(2);  kk = idx(3)
    if (cstat(i,j,kk) < 0) return
    ! wall face nodes (block node indices)
    nn(1,:) = idx;  nn(1,ax) = merge(1, ncell(b, ax) + 1, sd == 1)
    nn(2,:) = nn(1,:);  nn(2,t1) = nn(2,t1) + 1
    nn(3,:) = nn(1,:);  nn(3,t1) = nn(3,t1) + 1;  nn(3,t2) = nn(3,t2) + 1
    nn(4,:) = nn(1,:);  nn(4,t2) = nn(4,t2) + 1
    x1 = [blk(b)%xn(nn(1,1),nn(1,2),nn(1,3)), blk(b)%yn(nn(1,1),nn(1,2),nn(1,3)), blk(b)%zn(nn(1,1),nn(1,2),nn(1,3))]
    x2 = [blk(b)%xn(nn(2,1),nn(2,2),nn(2,3)), blk(b)%yn(nn(2,1),nn(2,2),nn(2,3)), blk(b)%zn(nn(2,1),nn(2,2),nn(2,3))]
    x3 = [blk(b)%xn(nn(3,1),nn(3,2),nn(3,3)), blk(b)%yn(nn(3,1),nn(3,2),nn(3,3)), blk(b)%zn(nn(3,1),nn(3,2),nn(3,3))]
    x4 = [blk(b)%xn(nn(4,1),nn(4,2),nn(4,3)), blk(b)%yn(nn(4,1),nn(4,2),nn(4,3)), blk(b)%zn(nn(4,1),nn(4,2),nn(4,3))]
    Av = 0.5_dp*cross3(x3 - x1, x4 - x2)                       ! area vector [m2]
    amag = sqrt(sum(Av*Av))
    if (amag <= 0.0_dp) return
    xf = 0.25_dp*(x1 + x2 + x3 + x4)
    xc = Xh(:, 2*i, 2*j, 2*kk)*Lch                             ! adjacent cell centre [m]
    if (sum(Av*(xc - xf)) < 0.0_dp) Av = -Av                   ! into the fluid
    nrm = Av/amag
    ! pressure (first cell, zero gradient) [Pa]
    pw = p(i,j,kk)*rho1*Uch*Uch
    ! friction: first-cell velocity relative to the wall, tangential part
    uc(1) = 0.5_dp*(u(i-1,j,kk) + u(i,j,kk))
    uc(2) = 0.5_dp*(v(i,j-1,kk) + v(i,j,kk))
    uc(3) = 0.5_dp*(w(i,j,kk-1) + w(i,j,kk))
    xnd = xf/Lch
    ! wall velocity: the ramped Dirichlet value of the patch (u_init -> bc_val over t_ramp,
    ! as boundary.F90 applies it: a towing-tank start moves the wall with the fluid during
    ! the ramp; update 20, before it the ramp term was missing here) + the body motion
    hx = 2*i;  hy = 2*j;  hz = 2*kk
    select case (ax)
    case (1);  hx = merge(2*i-1, 2*i+1, sd == 1)
    case (2);  hy = merge(2*j-1, 2*j+1, sd == 1)
    case default;  hz = merge(2*kk-1, 2*kk+1, sd == 1)
    end select
    uw = sixdof_wall_velocity(b, xnd, [hx, hy, hz]) + [u_init, v_init, w_init]/Uch &
         + bc_ramp*([bc_val(1,n), bc_val(2,n), bc_val(3,n)] - [u_init, v_init, w_init])/Uch
    urel = (uc - uw)*Uch                                       ! [m/s]
    ut = urel - sum(urel*nrm)*nrm
    dn = max(abs(sum((xc - xf)*nrm)), 1.0e-12_dp)              ! [m]
    mu_eff = mu_h(hx,hy,hz)*mu1                                ! [kg/(m s)], incl. wall function
    tauw = mu_eff*ut/dn
    ok = .true.
  end subroutine wall_face_stress

  pure logical function is_wall_patch(b, n)
    integer, intent(in) :: b, n
    is_wall_patch = (all(blk(b)%bc_flag(1:3,n) == 1) .and. all(blk(b)%bc_val(1:3,n) == 0.0_dp)) &
                    .or. all(blk(b)%bc_flag(1:3,n) == 11)
  end function is_wall_patch

  pure subroutine tang_axes(ax, t1, t2)
    integer, intent(in)  :: ax
    integer, intent(out) :: t1, t2
    select case (ax)
    case (1);  t1 = 2;  t2 = 3
    case (2);  t1 = 1;  t2 = 3
    case default;  t1 = 1;  t2 = 2
    end select
  end subroutine tang_axes

  !=====================================================================
  ! Export (rank 0) every Data_Intvl steps.  Bodies: own wrench in B and
  ! E, absolute motion.  Children / sub-children: own wrench in the
  ! MOTHER frame (moment about the own O), motion relative to the mother.
  !=====================================================================
  subroutine sixdof_export(istep, t)
    integer,  intent(in) :: istep
    real(dp), intent(in) :: t
    integer :: k, m
    real(dp) :: fE(3,3), mE(3,3), fB(3,3), mB(3,3), R(3,3), xO(3), xG(3), r2d, ea(3), oma(3)
    if (.not. sixdof_on .or. myrank /= 0) return
    if (mod(istep, data_interval) /= 0 .and. istep /= 0) return
    r2d = 180.0_dp/pi
    do k = 1, nbd                                           ! PID controllers of the parts (update 21)
       do m = 1, bd(k)%npid
          call pid_export(bd(k)%ipid(m), istep, t)
       end do
    end do
    do k = 1, nbd
       m = bd(k)%par
       fE(:,1) = matmul(Fm, bd(k)%ff);  fE(:,2) = matmul(Fm, bd(k)%fp);  fE(:,3) = fE(:,1) + fE(:,2)
       mE(:,1) = matmul(Fm, bd(k)%mf);  mE(:,2) = matmul(Fm, bd(k)%mp);  mE(:,3) = mE(:,1) + mE(:,2)
       if (m == 0) then
          R = bd(k)%R                                   ! body: own frame
       else
          R = bd(m)%R                                   ! child: the mother's frame
       end if
       fB(:,1) = matmul(transpose(R), fE(:,1));  fB(:,2) = matmul(transpose(R), fE(:,2));  fB(:,3) = fB(:,1) + fB(:,2)
       mB(:,1) = matmul(transpose(R), mE(:,1));  mB(:,2) = matmul(transpose(R), mE(:,2));  mB(:,3) = mB(:,1) + mB(:,2)
       write(bd(k)%iu(1),'(I8,ES16.8,9ES16.8)') istep, t, fB(:,1), fB(:,2), fB(:,3)
       write(bd(k)%iu(2),'(I8,ES16.8,9ES16.8)') istep, t, mB(:,1), mB(:,2), mB(:,3)
       xO = bd(k)%rot0 + bd(k)%dg
       xG = matmul(bd(k)%Rg, bd(k)%cg0 - bd(k)%rot0) + xO
       if (m == 0) then
          write(bd(k)%iu(3),'(I8,ES16.8,9ES16.8)') istep, t, fE(:,1), fE(:,2), fE(:,3)
          write(bd(k)%iu(4),'(I8,ES16.8,9ES16.8)') istep, t, mE(:,1), mE(:,2), mE(:,3)
          write(bd(k)%iu(5),'(I8,ES16.8,30ES16.8)') istep, t, bd(k)%eta(1:3), bd(k)%eta(4:6)*r2d, &
             bd(k)%nu(1:3), bd(k)%nu(4:6)*r2d, bd(k)%etad(1:3), bd(k)%etad(4:6)*r2d, &
             bd(k)%nud(1:3), bd(k)%nud(4:6)*r2d, xO, xG
          flush(bd(k)%iu(3));  flush(bd(k)%iu(4))
       else
          ea = euler_of(bd(k)%R);  oma = matmul(transpose(bd(k)%R), bd(k)%OmE)
          write(bd(k)%iu(5),'(I8,ES16.8,36ES16.8)') istep, t, bd(k)%eta(1:3), bd(k)%eta(4:6)*r2d, &
             bd(k)%nu(1:3), bd(k)%nu(4:6)*r2d, bd(k)%etad(1:3), bd(k)%etad(4:6)*r2d, &
             bd(k)%nud(1:3), bd(k)%nud(4:6)*r2d, xO, xG, ea*r2d, oma*r2d
       end if
       flush(bd(k)%iu(1));  flush(bd(k)%iu(2));  flush(bd(k)%iu(5))
    end do
    if (kbg > 0 .and. bg_open) then
       write(iu_bg,'(I8,ES16.8,15ES16.8)') istep, t, bd(kbg)%eta(1:3), bd(kbg)%eta(4:6)*r2d, bd(kbg)%etad(1:3), &
          bd(kbg)%etad(4:6)*r2d, bd(kbg)%rot0 + bd(kbg)%dg
       flush(iu_bg)
    end if
  end subroutine sixdof_export

  ! Euler angles (phi, theta, psi) of R = Rz(psi) Ry(theta) Rx(phi)
  pure function euler_of(R) result(e)
    real(dp), intent(in) :: R(3,3)
    real(dp) :: e(3)
    e(2) = asin(max(-1.0_dp, min(1.0_dp, -R(3,1))))
    e(1) = atan2(R(3,2), R(3,3))
    e(3) = atan2(R(2,1), R(1,1))
  end function euler_of

  ! ---- restart support (update 13b) ------------------------------------
  ! number of state records in the restart file: the parts + the background node (update 20)
  integer function sixdof_nbody()
    sixdof_nbody = merge(nbd, 0, sixdof_on)
    if (sixdof_on .and. kbg > 0) sixdof_nbody = nbd + 1       ! the background node kbg = nbd + 1
  end function sixdof_nbody
  logical function sixdof_on_flag()
    sixdof_on_flag = sixdof_on
  end function sixdof_on_flag
  ! state vector of part k (42 values): eta, etad, etadd, nu, nud, fp, ff, mp, mf
  subroutine sixdof_get_state(k, st)
    integer,  intent(in)  :: k
    real(dp), intent(out) :: st(42)
    st(1:6) = bd(k)%eta;  st(7:12) = bd(k)%etad;  st(13:18) = bd(k)%etadd
    st(19:24) = bd(k)%nu;  st(25:30) = bd(k)%nud
    st(31:33) = bd(k)%fp;  st(34:36) = bd(k)%ff;  st(37:39) = bd(k)%mp;  st(40:42) = bd(k)%mf
  end subroutine sixdof_get_state
  subroutine sixdof_set_state(k, st)
    integer,  intent(in) :: k
    real(dp), intent(in) :: st(42)
    bd(k)%eta = st(1:6);  bd(k)%etad = st(7:12);  bd(k)%etadd = st(13:18)
    bd(k)%nu = st(19:24);  bd(k)%nud = st(25:30)
    bd(k)%fp = st(31:33);  bd(k)%ff = st(34:36);  bd(k)%mp = st(37:39);  bd(k)%mf = st(40:42)
    if (k == kbg) bg_first = .false.            ! the filtered background pose continues from the file
  end subroutine sixdof_set_state
  ! restart: move the grids to the restored poses and re-assemble the
  ! overset (the restored fields already live on the moved grids: no
  ! flux rotation, no fill of uncovered cells)
  subroutine sixdof_restore_pose(t)
    real(dp), intent(in) :: t
    integer :: k, c
    logical :: moved
    if (.not. sixdof_on .or. .not. any_moving) return
    moved = .false.
    do k = 1, nbd                             ! mothers first (composition)
       ! velocity-Verlet half-kick flag (update 21 fix): a free DOF was advanced as free in the
       ! step before the restart, so its second half-kick with the next acceleration is pending
       do c = 1, 6
          bd(k)%vhalf(c) = (bd(k)%dof(c) == 2 .and. t >= Motion_Start_T - 1.0e-12_dp .and. t > 0.0_dp)
       end do
       if (bd(k)%pid_rot) then               ! PID-rotated part (update 21): angle and rate from the restored states
          bd(k)%pid_ang  = pidc(bd(k)%ipid(1))%angle
          bd(k)%pid_angd = dot_product(bd(k)%pid_axB, bd(k)%nu(4:6))
       end if
       if (.not. bd(k)%moving) cycle
       call update_transform(k)
       bd(k)%Rg_prev = bd(k)%Rg;  bd(k)%dg_prev = bd(k)%dg
       moved = .true.
    end do
    if (kbg > 0) then                         ! background (update 20): restored filtered pose
       call update_transform(kbg);  bd(kbg)%Rg_prev = bd(kbg)%Rg;  bd(kbg)%dg_prev = bd(kbg)%dg
       bg_prev = bd(kbg)%eta;  call bgr_publish
       moved = .true.
    end if
    if (.not. moved) return
    restoring = .true.
    call move_grids
    call after_motion
    restoring = .false.
    if (myrank == 0) write(*,'(A)') '  6-DoF: grids moved to the restored body poses, overset re-assembled'
  end subroutine sixdof_restore_pose

  !---------------------------------------------------------------------
  ! Screen monitor (rank 0, at the log interval): one quantity per
  ! line, labels right-aligned to the colon, an empty line between the
  ! parts (update 13 layout); children / sub-children indented, their
  ! motion relative to the mother.
  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  ! Force-model interface (update 19)
  !---------------------------------------------------------------------
  ! node index of (body b, child c, sub-child s); 0 when it does not exist
  integer function sixdof_node_of(b, c, s)
    integer, intent(in) :: b, c, s
    integer :: k, lev
    sixdof_node_of = 0
    if (.not. sixdof_on) return
    lev = 1
    if (c > 0) lev = 2
    if (s > 0) lev = 3
    do k = 1, nbd
       if (bd(k)%lev == lev .and. bd(k)%bid == b .and. bd(k)%cid == c .and. bd(k)%sid == s) then
          sixdof_node_of = k;  return
       end if
    end do
  end function sixdof_node_of

  ! current pose of node k: grid motion x = Rg (x0 - O0) + O0 + dg [m], velocity of a point
  ! vG + omG x (x - O) [m/s], rotation R (B_k -> E), state eta / nu (bodies absolute, children relative)
  subroutine sixdof_node_pose(k, Rg, O0, dg, vG, omG, R, eta, nu)
    integer,  intent(in)  :: k
    real(dp), intent(out) :: Rg(3,3), O0(3), dg(3), vG(3), omG(3), R(3,3), eta(6), nu(6)
    Rg = bd(k)%Rg;  O0 = bd(k)%rot0;  dg = bd(k)%dg;  vG = bd(k)%vG;  omG = bd(k)%omG
    R = bd(k)%R;  eta = bd(k)%eta;  nu = bd(k)%nu
  end subroutine sixdof_node_pose

  ! external wrench on node k (grid frame [N], [N m] about the current O_k), replaced every step
  ! external wrench of the force models (update 19): cleared by force_compute at the start of
  ! every evaluation, then ADDED model by model (update 21 fix: several models on one mother)
  subroutine sixdof_clear_ext
    if (allocated(bd)) then
       bd(:)%fext(1) = 0.0_dp;  bd(:)%fext(2) = 0.0_dp;  bd(:)%fext(3) = 0.0_dp
       bd(:)%mext(1) = 0.0_dp;  bd(:)%mext(2) = 0.0_dp;  bd(:)%mext(3) = 0.0_dp
    end if
  end subroutine sixdof_clear_ext
  subroutine sixdof_set_ext(k, F, M)
    integer,  intent(in) :: k
    real(dp), intent(in) :: F(3), M(3)
    bd(k)%fext = bd(k)%fext + F;  bd(k)%mext = bd(k)%mext + M
  end subroutine sixdof_set_ext

  function sixdof_part_name(k) result(s)
    integer, intent(in) :: k
    character(len=:), allocatable :: s
    s = part_name(k)
  end function sixdof_part_name

  ! FSI sub-iteration (update 24): the grids again from the original coordinates with the rigid poses
  ! of the step and the new wall deformation (the morph hook), then the overset / flux update
  subroutine sixdof_remesh
    call move_grids
    call after_motion(.true.)
  end subroutine sixdof_remesh

  ! flexible-part definition of node k (update 24)
  function sixdof_fsi_def(k) result(f)
    integer, intent(in) :: k
    type(fsidef_t) :: f
    f = bd(k)%fsi
  end function sixdof_fsi_def
  ! mass data of node k: Mass of the namelist (set when computed from Density), CG, O
  subroutine sixdof_mass_data(k, mass, cg, O, set_mass)
    integer,  intent(in) :: k
    real(dp), intent(inout) :: mass
    real(dp), intent(out) :: cg(3), O(3)
    logical,  intent(in) :: set_mass
    if (set_mass) bd(k)%mass = mass
    mass = bd(k)%mass;  cg = bd(k)%cg0;  O = bd(k)%rot0
  end subroutine sixdof_mass_data
  ! frame kinematics of node k in the GRID frame (update 24): rotation Rg, O0, dg (x = Rg (x0 - O0) + O0 + dg),
  ! velocity / acceleration of O [m/s, m/s2], angular velocity / acceleration [rad/s, rad/s2]
  subroutine sixdof_frame_kin(k, Rg, O0, dg, vG, aG, omG, alG)
    integer,  intent(in)  :: k
    real(dp), intent(out) :: Rg(3,3), O0(3), dg(3), vG(3), aG(3), omG(3), alG(3)
    Rg = bd(k)%Rg;  O0 = bd(k)%rot0;  dg = bd(k)%dg;  vG = bd(k)%vG;  omG = bd(k)%omG
    aG = matmul(Fm, bd(k)%aE)
    alG = matmul(Fm, matmul(bd(k)%R, bd(k)%alB))
  end subroutine sixdof_frame_kin
  function sixdof_part_stem(k) result(s)          ! body1 / body1_child1 / body1_child1_sub1 (file names)
    integer, intent(in) :: k
    character(len=:), allocatable :: s
    s = part_stem(k)
  end function sixdof_part_stem

  subroutine sixdof_monitor
    integer :: k, na, nf, nh, no, ne
    real(dp) :: r2d, fE(3), mE(3)
    character(len=6) :: ind
    if (.not. sixdof_on .or. myrank /= 0) return
    r2d = 180.0_dp/pi
    write(*,'(A)') ' '
    do k = 1, nbd
       fE = matmul(Fm, bd(k)%ff + bd(k)%fp);  mE = matmul(Fm, bd(k)%mf + bd(k)%mp)
       ind = repeat(' ', 2*(bd(k)%lev - 1))
       write(*,'(A,3ES13.4)')   ' [ 6-DoF ] '//trim(ind)//part_name(k)//'   force (E, SNAME) [N]: ', fE
       write(*,'(A,3ES13.4,A)') '                   moment about O [N m]: ', mE, merge('  (moving)', '  (static)', bd(k)%moving)
       if (bd(k)%lev == 1) then
          if (bd(k)%Mc /= bd(k)%mass .or. any(bd(k)%Fs /= fE)) &
             write(*,'(A,3ES13.4)') '          with the children, force [N]: ', bd(k)%Fs
          write(*,'(A,3ES13.4)')   '                   surge/sway/heave [m]: ', bd(k)%eta(1:3)
          write(*,'(A,3F13.4)')    '                   roll/pitch/yaw [deg]: ', bd(k)%eta(4:6)*r2d
          write(*,'(A,3ES13.4)')   '                            u v w [m/s]: ', bd(k)%nu(1:3)
          write(*,'(A,3ES13.4)')   '                          p q r [deg/s]: ', bd(k)%nu(4:6)*r2d
       else
          write(*,'(A,3ES13.4)')   '      relative to the mother, x y z [m]: ', bd(k)%eta(1:3)
          write(*,'(A,3F13.4)')    '          relative roll/pitch/yaw [deg]: ', bd(k)%eta(4:6)*r2d
          write(*,'(A,3ES13.4)')   '                   relative u v w [m/s]: ', bd(k)%nu(1:3)
          write(*,'(A,3ES13.4)')   '                 relative p q r [deg/s]: ', bd(k)%nu(4:6)*r2d
       end if
       write(*,'(A)') ' '
    end do
    if (n_reassemble > 0) then
       write(*,'(A,I0,A,F10.3,A)') ' overset re-assemblies : ', n_reassemble, &
          ' (average ', t_reassemble/n_reassemble, ' s each)'
       call ovs_counts(na, nf, nh, no, ne)
       write(*,'(A,I0,A,I0,A,I0)') ' cells active / fringe / hole : ', na, ' / ', nf, ' / ', nh
       write(*,'(A,I0,A,I0)') ' orphans: ', no, '   fringe entries: ', ne
    end if
  end subroutine sixdof_monitor

  subroutine sixdof_summary
    if (.not. sixdof_on .or. myrank /= 0) return
    if (n_reassemble > 0) write(*,'(A,I0,A,F10.2,A,I0,A)') '  6-DoF: ', n_reassemble, ' overset re-assemblies, ', &
       t_reassemble, ' s wall time in total (', n_skipped, ' moving steps re-used the previous assembly, Ovst_Reassemble_Tol)'
  end subroutine sixdof_summary

  subroutine sixdof_finalize
    integer :: k, i
    if (.not. sixdof_on .or. myrank /= 0) return
    do k = 1, nbd
       do i = 1, 5
          if (bd(k)%iu(i) > 0) close(bd(k)%iu(i))
       end do
    end do
    if (bg_open) close(iu_bg)
    call pid_finalize
  end subroutine sixdof_finalize

end module sixdof
