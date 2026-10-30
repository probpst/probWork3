!=======================================================================
! mod_global.F90
!   Global data module for probCFD (multi-block, curvilinear, two-phase
!   capable, intra-block MPI splitting).
!
!   TWO-LEVEL DOMAIN DECOMPOSITION
!   ------------------------------
!   Level 1: the grid consists of nblk structured BLOCKS (PLOT3D).
!   Level 2: every block b is split into sp(1) x sp(2) x sp(3)
!            SUB-DOMAINS in its own i/j/k directions, as prescribed
!            per block in 'Input_polyMesh_probCFD.csv'.  One CPU
!            works on one sub-domain when nprocs equals the csv
!            total; with fewer ranks the sub-domains are distributed
!            round-robin.
!
!   A "domain" (type mdom_t) carries the field storage of one
!   sub-domain, allocated over its BLOCK-GLOBAL cell range plus ng
!   ghost layers, so all numerical kernels keep indexing with the
!   block-global i, j, k.  set_block(d) points the familiar module
!   names (u, v, w, p, alpha, nx, is, ie, ...) at domain d.
!
!   Ghost updates (intra-block halos AND block-to-block interfaces,
!   including differing block orientations) are handled by one
!   unified exchange-record engine in b2b_exchange.F90.
!
!   FIELD IDS used by the exchange engine (field_ptr / comp_axis):
!      1..3  U1,U2,U3   transformed face velocity components
!      4     p          pressure (cell)
!      5..7  Us1..Us3   intermediate (star) face components
!      8     alpha      liquid volume fraction (cell)
!      9,10  kv1, kv2   Krylov work vectors of the pressure solver
!                       (cell, homogeneous pressure BCs)
!      11..13 G1..G3    transformed pressure-gradient face fluxes of the
!                       Poisson operator (stored in rhsu/rhsv/rhsw,
!                       face-type like U; exchanged, no BCs)
!      14..16 tk,tw,mut turbulence variables / eddy viscosity (cell)
!      17,18  alps,kap  surface tension (cell)
!      19     vap       vapour volume fraction (cavitation, cell)
!      20     temp      temperature [K] (energy equation, cell; update 16)
!=======================================================================
module global_vars
  implicit none
  save

  integer, parameter :: dp = selected_real_kind(15, 307)
  integer, parameter :: ng = 3          ! ghost layers (3: 2nd-order upwind
                                        !   evaluated ON b2b interface faces
                                        !   reaches 2 faces + 1 cell outside)

  !------------------ namelist inputs (probCFD.nml) --------------------
  character(len=64) :: casename            ! case name (grid = casename.grd, bc = casename.bcs)
  real(dp) :: Lch                          ! characteristic length  [m]
  real(dp) :: Uch                          ! characteristic velocity [m/s]
  integer  :: Cnvt_scheme                  ! 1 = 1st-order upwind, 2 = 2nd-order upwind,
                                           !   3 = QUICK, 4 = 2nd-order ENO (convection.F90)
  ! face-variable formulation on the curvilinear grid:
  !   1 = Park's unweighted transformed components U^q = u.grad(xi_q)
  !   2 = J-weighted volume fluxes F^q = J U^q = b_q.u (Rosenfeld /
  !       CFDShip form: the stored face dof is the volume flow through
  !       the face; discrete continuity telescopes structurally and the
  !       advecting mass fluxes are interpolated as fluxes)
  integer  :: curvi_form
  ! time integration of the momentum terms (per term, per direction):
  !   Cnvt_mode : convection  1 = explicit, 2 = implicit,
  !               3 = hybrid: per i/j/k direction by the local CFL
  !                   number (explicit where CFL_q < cfl_switch)
  !   Dffs_mode : diffusion   1 = explicit, 2 = implicit,
  !               3 = hybrid: per direction by the local Von Neumann
  !                   number (explicit where VN_q < vn_switch)
  integer  :: Cnvt_mode, Dffs_mode
  integer  :: quasi2d = 0                   ! SW_Quasi_2D (update 24c): 3 = strictly 2-D across the grid's k (span) direction
  real(dp) :: cfl_switch, vn_switch
  ! pressure form of the fractional step (update 16):
  !   1 = non-incremental (existing): the predictor carries no pressure
  !       gradient, the Poisson equation gives the full pressure
  !   2 = incremental (Van Kan / Goda): the predictor includes the
  !       gradient of the latest pressure, the projection corrects with
  !       the pressure INCREMENT.  The body force then never passes
  !       through the implicit viscous operator of the predictor: with
  !       form 1 a fluid at rest under gravity develops spurious
  !       velocities next to no-slip walls wherever the wall-normal
  !       diffusion is implicit (VN >= vn_switch), e.g. on wall-clustered
  !       grids - use 2 for gravity / buoyancy-driven flows
  integer  :: Prs_form = 1
  ! phases and fluid properties
  integer  :: Num_of_Phase                 ! 1 = single phase, 2 = two-phase (VoF)
  real(dp) :: rho1, mu1                    ! fluid 1 = liquid (alpha = 1) [kg/m3], [kg/(m s)]
  real(dp) :: rho2, mu2                    ! fluid 2 = gas    (alpha = 0)
  real(dp) :: rho, mu                      ! legacy aliases of rho1, mu1
  real(dp) :: grav(3)                      ! gravity vector [m/s2] (SI)
  real(dp) :: z_free_surface               ! initial interface height [m]:
                                           !   liquid fills z < z_free_surface
  real(dp) :: vof_cfl_max                  ! VoF sub-cycling target Courant number
  ! time stepping and solvers
  real(dp) :: dt_dim                       ! time step [s]
  real(dp) :: t_ramp                       ! Dirichlet velocity BCs move linearly from
                                           !   the initial velocity (u_init,v_init,w_init)
                                           !   to their values over t_ramp [s] (0 = no ramp)
  real(dp) :: bc_ramp = 1.0_dp             ! current ramp factor (set per step)
  real(dp) :: u_init, v_init, w_init       ! initial velocities [m/s]
  real(dp) :: p_init                       ! initial pressure [Pa]
  real(dp) :: End_T                        ! end time [s] (nstep = End_T/Del_T)
  integer  :: nstep                        ! number of time steps to run = nint(End_T/Del_T)
  real(dp) :: alpha_init                   ! a_cel_initial: initial alpha ABOVE FreeSurf_Z
  ! ---- inner (outer-Picard) iterations per time step -----------------
  integer  :: Iner_Itn                     ! momentum + Poisson + correction repeated
                                           !   Iner_Itn times per time step, the advecting
                                           !   fluxes and the explicit terms re-evaluated
                                           !   from the latest iterate (1 = classic
                                           !   non-iterative fractional step)
  ! ---- momentum predictor (Jacobi sweeps) ----------------------------
  integer  :: mom_maxiter, mom_miniter     ! Max_Itr_Vel / Min_Itr_Vel
  real(dp) :: mom_tol                      ! Rsd_Lmt_Vel: ||b - A u||_inf / ||b||_inf
  real(dp) :: mom_sor                      ! SOR_fct_Vel: weighted-Jacobi factor (<= 1)
  ! ---- pressure-Poisson solver ----------------------------------------
  integer  :: poi_solver                   ! 1 = Jacobi, 2 = FGMRES + 7-point PCG preconditioner
  integer  :: poi_precond                  ! PCG preconditioner (solver 2): 1 = Jacobi,
                                           !   2 = line symmetric Gauss-Seidel (tridiagonal
                                           !   solves along the strongest-coupled direction)
  integer, parameter :: poi_mkry = 20      ! FGMRES restart length
  integer, parameter :: ovs_mkry = 20      ! overset: restart length of the inner GMRES (coupled operator)
  real(dp), parameter :: ovs_prelax = 0.2_dp ! overset: weight of the same-block neighbour average in the
                                           !   fringe rows of the pressure system (regularization of
                                           !   interpolation cycles, see poisson.F90 fringe_row)
  integer  :: poi_maxiter, poi_miniter     ! Max_Itr_Prs / Min_Itr_Prs
  real(dp) :: poi_tol                      ! Trl_Lmt_Prs: update norm max|p_new - p_old|
                                           !   (non-dim.), criterion of solver 1
  real(dp) :: poi_rtol                     ! Rsd_Lmt_Prs: ||b - A p||_inf / ||b||_inf,
                                           !   criterion of solver 2
  real(dp) :: poi_sor                      ! SOR_fct_Prs: weighted-Jacobi factor (solver 1)
  ! ---- VoF solver -------------------------------------------------------
  integer  :: vof_scheme                   ! 1 = explicit CICSAM, 2 = implicit CICSAM, 3 = geometric PLIC
  !------------------ turbulence (turbulence.F90) ------------------------
  integer  :: TURB_model = 0               ! 0 none, 1 Spalart-Allmaras, 2 k-epsilon, 3 k-omega, 4 k-omega SST
  integer  :: IDDES_model = 0              ! 1 = IDDES hybrid RANS/LES on top of TURB_model (0: WMLES only)
  integer  :: Max_Itr_Turb = 20, Min_Itr_Turb = 1
  real(dp) :: Rsd_Lmt_Turb = 1.0e-6_dp, SOR_fct_Turb = 1.0_dp
  integer  :: Cnvt_scheme_Turb = 1         ! 1 = 1st-order upwind, 2 = 2nd-order upwind (turbulence equations)
  integer  :: Cnvt_mode_Turb = 2           ! 1 = explicit, 2 = implicit (deferred correction for scheme 2)
  integer  :: Dffs_mode_Turb = 2           ! 1 = explicit, 2 = implicit
  real(dp) :: cfl_switch_Turb = 0.5_dp     ! explicit convection: sub-cycles so that the cell CFL <= this
  real(dp) :: vn_switch_Turb  = 0.5_dp     ! explicit diffusion : sub-cycles so that the cell VN  <= this
  logical  :: turb_on = .false.            ! any eddy viscosity active (TURB_model > 0 or IDDES_model > 0)
  !------------------ surface tension (surface_tension.F90) --------------
  integer  :: Surf_tension = 0             ! 1 = CSF surface-tension force (two-phase runs only)
  real(dp) :: Sigma_ST = 0.0728_dp         ! surface-tension coefficient [N/m]
  real(dp) :: sigma_nd = 0.0_dp            ! Sigma_ST / (rho1 Uch^2 Lch)  (1 / Weber number)
  integer  :: Init_shape = 0               ! initial alpha: 0 flat free surface, 1 cylinder (y), 2 sphere, 3 square prism (y)
  real(dp) :: Drop_Center(3) = 0.0_dp, Drop_Radius = 0.01_dp   ! [m]
  logical  :: st_on = .false.              ! Surf_tension = 1 and Num_of_Phase = 2
  !------------------ calculation-stop criteria (stop_criteria.F90) ------
  real(dp) :: Stop_Rsd_Max = 1.0e9_dp       ! stop when a residual |b-Ap|/|b| or an update norm max|dp| exceeds this
  real(dp) :: Stop_CFL_Max = 1.0e4_dp       ! stop when the maximum local CFL number exceeds this
  ! rank-local maxima with their location (block, i, j, k) recorded by the solvers
  real(dp) :: poi_res_local = 0.0_dp, poi_upd_local = 0.0_dp, cfl_local = 0.0_dp
  integer  :: poi_res_loc(4) = 0, poi_upd_loc(4) = 0, cfl_loc(4) = 0
  real(dp) :: mom_res_local(3) = 0.0_dp
  integer  :: mom_res_loc(4,3) = 0
  logical  :: wallfn_on = .false.          ! any flag-11 (wall function) patch present
  !------------------ 6-DoF rigid-body motion (sixdof.F90) ---------------
  integer  :: data_interval = 1            ! Data_Intvl: export interval of the force / moment / motion files
  real(dp) :: grid_scale = 1.0_dp          ! Grid_scale: the .grd coordinates are multiplied by this factor (-> metres)
  integer, allocatable :: blk_obj(:)       ! 6-DoF NODE id of every block (0 = background), built by read_split_csv
                                           ! from rows 6-8 of the csv (Object_ID, Child_Object_ID, Sub_Child_Object_ID)
  ! node table (update 14): one node per rigid part - the bodies first
  ! (node = body id), then the child bodies (sorted by body, child),
  ! then the sub-child bodies (sorted by body, child, sub)
  integer :: n6node = 0
  integer, allocatable :: node_bid(:), node_cid(:), node_sid(:)   ! body / child / sub-child ids of a node (0 = none)
  integer, allocatable :: node_lev(:), node_par(:)                ! level 1 body, 2 child, 3 sub-child; parent node (0: body)
  logical  :: sixdof_on = .false.          ! any body defined (forces are integrated; bodies may move)
  logical  :: any_moving = .false.         ! any body with an imposed or free degree of freedom
  ! moving background (update 20): the blocks with Object_ID = 0 follow a Body; current rigid pose of the
  ! background x = bgr_Rg (x0 - bgr_O0) + bgr_O0 + bgr_dg [m] (grid frame), kept here for the modules
  ! compiled before sixdof.F90 (waves.F90: relaxation zones attached to the grid faces)
  logical  :: bgr_on = .false.
  real(dp) :: bgr_Rg(3,3) = reshape([1.0_dp,0.0_dp,0.0_dp, 0.0_dp,1.0_dp,0.0_dp, 0.0_dp,0.0_dp,1.0_dp],[3,3])
  real(dp) :: bgr_O0(3) = 0.0_dp, bgr_dg(3) = 0.0_dp
  ! wall velocity of the moving bodies at a grid position (set by
  ! sixdof_init; curvi.F90 uses it for the Cartesian Dirichlet faces
  ! without a module dependency on sixdof)
  abstract interface
     function wallvel_if(b, xnd, h) result(ug)
       import :: dp
       integer,  intent(in) :: b
       real(dp), intent(in) :: xnd(3)
       integer,  intent(in), optional :: h(3)     ! half index of the position (FSI blocks, update 24)
       real(dp) :: ug(3)
     end function wallvel_if
  end interface
  procedure(wallvel_if), pointer :: wall_velocity_fn => null()
  ! fluid-structure interaction (update 24, fsi.F90): hooks called by sixdof.F90 without a
  ! module dependency - the mesh morphing of the flexible parts after the rigid grid motion,
  ! the discrete grid velocity of a morphing block at a half-index position, the flux
  ! re-projection of a morphed sub-domain
  logical :: fsi_on = .false.
  logical, allocatable :: fsi_blk(:)             ! (nblk) the block morphs with a flexible part
  logical :: fsi_reasm_request = .false.         ! the flexible walls moved more than the re-assembly tolerance
  abstract interface
     subroutine fsi_morph_if(restoring)
       logical, intent(in) :: restoring
     end subroutine fsi_morph_if
     function fsi_hvel_if(b, h) result(ug)
       import :: dp
       integer, intent(in) :: b, h(3)
       real(dp) :: ug(3)
     end function fsi_hvel_if
     subroutine fsi_flux_if(d)
       integer, intent(in) :: d
     end subroutine fsi_flux_if
  end interface
  procedure(fsi_morph_if), pointer :: fsi_morph_fn => null()
  procedure(fsi_hvel_if),  pointer :: fsi_hvel_fn  => null()
  procedure(fsi_flux_if),  pointer :: fsi_flux_fn  => null()
  ! prescribed grid deformation from frame files (heart extension H1, morph_frames.F90): the node
  ! coordinates of the frame parts at time t (Mode = 2 grid-movement check, sixdof.F90)
  logical :: frames_on = .false.
  abstract interface
     subroutine frames_nodes_if(t)
       import :: dp
       real(dp), intent(in) :: t
     end subroutine frames_nodes_if
  end interface
  procedure(frames_nodes_if), pointer :: frames_nodes_fn => null()
  ! energy equation (update 16): the physical-boundary ghost values of T
  ! (Dirichlet / heat flux / convective) are set by energy.F90 after the
  ! generic apply_bc (zero gradient) inside refresh_all(20)
  abstract interface
     subroutine hook_if()
     end subroutine hook_if
  end interface
  procedure(hook_if), pointer :: temp_bc_hook => null()
  logical, allocatable :: obj_moving(:)    ! (object id) .true. for the objects that move
                                           !   (iterative, residual criterion)
  integer  :: vof_maxiter, vof_miniter     ! Max_Itr_VoF / Min_Itr_VoF (scheme 2)
  real(dp) :: vof_tol                      ! Rsd_Lmt_VoF: max |alpha^n - alpha - dt/V sum F alpha_f|
  real(dp) :: vof_sor                      ! SOR_fct_VoF: weighted-Jacobi factor (scheme 2)
  integer  :: out_interval                 ! Output_Intvl: export interval (time steps)
  integer  :: sw_prs_frc_surf = 0          ! SW_Prs_Frc_surf (update 22): 1 = write the solid-wall surface files
                                           !   <case>_surf_NNNNNNNN.grd / .p3d (pressure and friction vectors [Pa]) with the volume files
  integer  :: frame_prs_frc_surf = 0       ! Frame_Prs_Frc_surf: 0 earth-fixed frame, 1 background (moving-grid) frame
  ! ---- turbulence statistics (update 23, $Statistics; evaluated by statistics.F90) --------------
  integer  :: sw_statistics = 0            ! SW_statistics: 1 = on
  integer  :: stat_type(4) = 0             ! Statistics_type1..4: original grid, boxes, cylinders, planes
  integer  :: stat_group(2) = 0            ! Statistics_first_group / _second_group
  real(dp) :: stat_start_t = 0.0_dp        ! Statistics_Start_T [s]
  real(dp) :: stat_end_t = 0.0_dp          ! Statistics_End_T [s] (0 = End_T)
  integer  :: stat_sample_intvl = 1        ! Statistics_Sample_Intvl (steps)
  integer  :: stat_output_intvl = 0        ! Statistics_Output_Intvl (steps; 0 = end of the run only)
  integer  :: stat_type1_frame = 1         ! Statistics_Type1_Frame: 1 = each block's part frame, 0 = grid axes
  integer, parameter :: MAXLATT = 64
  type :: latt_def_t                       ! a user-defined lattice (box 2 / cylinder 3 / plane 4)
     integer  :: kind = 0, id = 0, mb = 0, mc = 0, ms = 0
     real(dp) :: org(3) = 0.0_dp, ax1(3) = [1.0_dp, 0.0_dp, 0.0_dp], ax2(3) = [0.0_dp, 1.0_dp, 0.0_dp]
     real(dp) :: len(3) = 0.0_dp, rad(2) = 0.0_dp, theta(2) = [0.0_dp, 360.0_dp]
     integer  :: np(3) = 1, polar = 0
     logical  :: open = .false.
  end type latt_def_t
  type(latt_def_t) :: latt_def(MAXLATT)
  integer  :: nlatt_def = 0
  ! ---- probe points of the background (update 22, $Probe_background; the parts' points are in sixdof.F90)
  integer, parameter :: MAXPRBG = 256
  integer  :: nprobe_bg = 0
  integer  :: probe_bg_n(MAXPRBG) = 0       ! point number of Probe_point<n>
  real(dp) :: probe_bg_x(3, MAXPRBG) = 0.0_dp   ! grid frame at t = 0 [m]
  integer  :: log_interval                 ! Log_Intvl: step-monitor interval (time steps)
  real(dp) :: b2b_tol                      ! block-to-block point matching tolerance
  ! ---- overset (Chimera) grids, see overset.F90 -----------------------
  integer  :: ovst_switch = -1             ! SW_Ovst: -1 auto (on when a block face is
                                           !   neither in .bcs nor in .bcs_b2b), 0 off, 1 on
  real(dp) :: surfovl_tol = 1.0e-4_dp      ! Surfovl_tol [m]: tolerance of the wall-surface
                                           !   overlap / domain-boundary contact tests
  integer  :: ovst_nf = 2                  ! Ovst_fringe: fringe (interpolated) cell layers
  integer  :: ovst_wall_layers = 8         ! Ovst_wall_layers: boundary-layer cells of a wall block that stay active
  integer  :: ovst_active_overlap = 1      ! Ovst_active_overlap: 1 = active regions of overlapping grids overlap by one cell
  real(dp) :: ovst_reasm_tol = 0.0_dp       ! Ovst_Reassemble_Tol [m] (update 20): the dynamic overset is re-assembled only
                                           !   when a moving part has moved more than this relative to the background
                                           !   since the last assembly (0 = every step)
  real(dp) :: ovst_vtol = 1.0e-6_dp        ! Ovst_vol_tol: relative cell-volume tolerance of
                                           !   the "finer grid wins" decision (ties: both active)
  character(len=*), parameter :: vol_dir = 'Volume_solutions_grids/'   ! sub-folder of the volume grid / solution files (update 13)
  integer  :: run_mode = 0                 ! Mode: 0 = standard (from step 0), 1 = restart from the latest <case>_NNNNNNNN.rsto (update 13b)
  integer  :: restart_interval = 0         ! Restart_Intvl: restart file every N steps (0 = none; the last step is always written when > 0)
  integer  :: restart_step = 0             ! step of the restart file that was read (Mode = 1)
  integer  :: Time_scheme = 1              ! Time_scheme: 1 = existing (first order in time), 2 = 2nd-order SSP (Heun) (update 13)
  real(dp) :: t_now = 0.0_dp               ! current physical time [s] (accumulated: Del_T may change at run time)
  integer  :: ovst_only = 0                ! Ovst_only: 1 = stop after the assembly (writes
                                           !   the report and the iblank solution file)
  logical  :: ovst_on = .false.            ! overset assembly active (set by ovst_setup)
  logical  :: poi_verbose = .false.        ! print the FGMRES restart history (overset)
  integer  :: ovst_inner_max = 150         ! Ovst_inner_max: BiCGSTAB iterations per FGMRES
                                           !   preconditioner application (overset pressure solve)
  real(dp) :: ovst_inner_tol = 1.0e-2_dp   ! Ovst_inner_tol: its relative residual reduction
  logical  :: p_ref_switch                 ! .true. when ALL pressure BCs are Neumann
  real(dp) :: p_ref_x, p_ref_y, p_ref_z    ! reference pressure location [m]
  real(dp) :: p_ref_val                    ! reference pressure value [Pa]
  integer  :: p_ref_blk                    ! block containing the reference location
  integer  :: p_ref_i, p_ref_j, p_ref_k    ! cell containing the reference location
  integer  :: p_ref_dom = 0                ! sub-domain containing that cell

  !------------------ derived non-dimensional parameters ---------------
  real(dp) :: tch                          ! characteristic time = Lch/Uch [s]
  real(dp) :: Re                           ! Reynolds number = rho1*Uch*Lch/mu1
  real(dp) :: dt                           ! non-dimensional time step = dt_dim/tch
  real(dp) :: rho_ratio, mu_ratio          ! rho2/rho1, mu2/mu1
  !------------------ cavitation (update 15, cavitation.F90) --------------
  logical  :: cav_on = .false.             ! SW_Cav = 1: vapour volume fraction + mass transfer
  integer  :: Cav_model = 1                ! 1 Schnerr-Sauer, 2 Zwart-Gerber-Belamri, 3 Kunz
  integer  :: Vapor_scheme = 2             ! vapour transport: 1 explicit CICSAM, 2 implicit CICSAM, 3 ENO (minmod) + Cnvt_mode
  real(dp) :: Cav_Start_T = 0.0_dp         ! mass transfer switched on from this time [s] (50-step ramp)
  real(dp) :: Prs_Vapor = 2340.0_dp        ! vapour pressure [Pa]
  real(dp) :: rho_v = 0.0173_dp, mu_v = 9.8e-6_dp   ! vapour density [kg/m3], viscosity [kg/(m s)]
  real(dp) :: SS_Nuclei = 1.0e13_dp, SS_Radius = 1.0e-6_dp               ! Schnerr-Sauer n0 [1/m3], R0 [m]
  real(dp) :: ZGB_Fvap = 50.0_dp, ZGB_Fcond = 0.01_dp, ZGB_Rnuc = 5.0e-4_dp, ZGB_Rb = 1.0e-6_dp   ! Zwart
  real(dp) :: Kunz_Cdest = 100.0_dp, Kunz_Cprod = 100.0_dp               ! Kunz (t_inf = Lch/Uch)
  integer  :: Reboud = 0                   ! 1 = Reboud eddy-viscosity correction (n = 10) in cavitating cells
  real(dp) :: rv_ratio = 1.0_dp, muv_ratio = 1.0_dp   ! rho_v/rho1, mu_v/mu1
  real(dp) :: pv_nd = 0.0_dp               ! vapour pressure / (rho1 Uch^2)
  real(dp) :: p_level = 0.0_dp             ! absolute pressure level of the case [Pa] (BC / RefPrs / p_init; sigma print, Dym_Prs)
  real(dp) :: cav_ramp = 1.0_dp            ! current source ramp factor (0..1)
  real(dp) :: Latent_Heat = 2.2567e6_dp    ! latent heat of vaporisation [J/kg] (water; energy equation sink -mdot L)
  integer  :: Pv_of_T = 0                  ! 1 = vapour pressure from the saturation law p_v(T) (energy on; else Prs_Vapor)
  ! ---- update 17: compressible mixture and thermal coupling ------------
  integer  :: Cav_compressible = 0         ! 1 = barotropic compressibility of the liquid-vapour mixture (Wood / Wallis
                                           !   sound speed: liquid Snd_liquid, vapour = isothermal ideal gas c_v^2 = p/rho_v)
                                           !   in the pressure equation: (1/(rho c^2)) dp/dt + div u = S
  real(dp) :: Snd_liquid = 1480.0_dp       ! speed of sound of the liquid [m/s] (water 20 C: 1482)
  integer  :: Cav_thermal = 0              ! 1 = thermal cavitation (needs SW_Energy = 1): Clausius-Clapeyron p_v(T) and
                                           !   rho_v(T) through (Prs_Vapor, Lou_v, Latent_Heat at T_init), thermodynamic
                                           !   delay of the evaporation (B-factor), implicit latent-heat sink
  real(dp) :: Cav_thermal_t = 0.0_dp       ! bubble-growth time [s] of the delay dT* = 1/2 (rho_v L/(rho_l c_p)) R_B/sqrt(alpha t) (0 = 1 ms)
  real(dp) :: vap_init = 0.0_dp            ! vap_cel_initial: initial vapour fraction of the liquid cells (tests)
  !------------------ heat / energy equation (update 16, energy.F90) -----
  !   T [K] is carried DIMENSIONAL at the cell centres (field id 20):
  !     dT/dt* + u* . grad* T = div*( Gam* grad* T ) + S*,
  !     Gam* = k_eff / ((rho c_p) Uch Lch),  k_eff = k + mu_t c_p / Pr_t
  !   (one-fluid mixture: rho c_p and k volume-weighted over liquid /
  !   air / vapour).  Buoyancy by the Boussinesq term -beta (T - T_init) g
  !   in the momentum equation (Beta_l / Beta_a > 0).  Cavitation: latent
  !   heat sink -mdot L / (rho c_p) and, with Pv_of_T = 1, p_v(T).
  logical  :: energy_on = .false.          ! SW_Energy = 1
  real(dp) :: T_init = 293.15_dp           ! initial temperature [K] (also the Boussinesq reference)
  real(dp) :: Cp_l = 4182.0_dp, k_l = 0.6_dp          ! liquid  specific heat [J/(kg K)], conductivity [W/(m K)]
  real(dp) :: Cp_a = 1005.0_dp, k_a = 0.0257_dp       ! air (two-phase)
  real(dp) :: Cp_v = 1996.0_dp, k_v = 0.0184_dp       ! vapour (cavitation)
  real(dp) :: Beta_l = 0.0_dp, Beta_a = 0.0_dp        ! thermal expansion [1/K] (0 = no buoyancy coupling)
  real(dp) :: Pr_t = 0.85_dp               ! turbulent Prandtl number
  integer  :: Max_Itr_Engy = 50, Min_Itr_Engy = 1     ! SOR sweeps of the temperature equation
  real(dp) :: Rsd_Lmt_Engy = 1.0e-6_dp, SOR_fct_Engy = 1.0_dp
  integer  :: Cnvt_scheme_Engy = 2         ! 1 1st-order upwind, 2 2nd-order upwind, 3 QUICK, 4 2nd-order ENO (minmod)
  integer  :: Cnvt_mode_Engy = 3           ! convection 1 explicit (sub-cycled), 2 implicit, 3 hybrid per direction (cfl_switch_Engy)
  integer  :: Dffs_mode_Engy = 3           ! diffusion  1 explicit (sub-cycled), 2 implicit, 3 hybrid per direction (vn_switch_Engy)
  real(dp) :: cfl_switch_Engy = 0.3_dp, vn_switch_Engy = 0.5_dp
  real(dp) :: rcp_a = 1.0_dp, rcp_v = 1.0_dp          ! (rho c_p)_a / (rho c_p)_l, (rho c_p)_v / (rho c_p)_l
  real(dp) :: kr_a = 1.0_dp, kr_v = 1.0_dp            ! k_a / k_l, k_v / k_l
  real(dp) :: Pe_l = 1.0_dp                ! Peclet number rho_l c_p,l Uch Lch / k_l
  real(dp) :: Pr_l = 1.0_dp, Pr_a = 1.0_dp, Pr_v = 1.0_dp   ! molecular Prandtl numbers (thermal wall function)
  logical  :: buoy_on = .false.            ! any Beta /= 0 (Boussinesq term active)
  real(dp) :: gnd(3)                       ! gravity * Lch / Uch^2 (= 1/Fr^2 per axis)
  real(dp) :: zfs                          ! z_free_surface / Lch

  !------------------ one grid block (bc + split info) -----------------
  type :: mblk_t
     integer :: ni = 0, nj = 0, nk = 0     ! node counts
     integer :: nx = 0, ny = 0, nz = 0     ! cell counts
     ! MPI splitting of this block (from the csv file)
     integer :: sp(3) = 1                  ! i/j/k split counts
     integer :: ncpu  = 1                  ! CPUs for this block = product(sp)
     integer :: dom0  = 0                  ! global id of its first sub-domain
     ! node coordinates in meters (replicated on every rank)
     real(dp), allocatable :: xn(:,:,:), yn(:,:,:), zn(:,:,:)
     ! ghost-extended NON-DIMENSIONAL node coordinates (1-ng:ni+ng, ...):
     ! filled by extrapolation at physical boundaries and from the donor
     ! block (plus the translation vector for periodic pairs) at
     ! interfaces
     real(dp), allocatable :: xg(:,:,:,:)          ! (3, i,j,k)
     ! physical boundary condition patches of this block
     integer :: nbc = 0
     integer,  allocatable :: bc_dir(:)
     integer,  allocatable :: bc_irng(:,:), bc_jrng(:,:), bc_krng(:,:)
     integer,  allocatable :: bc_flag(:,:)         ! (6, patch): u,v,w,p,alpha,T (Cartesian)
     character(len=48), allocatable :: bc_label(:) ! (patch): the comment of the direction line of the .bcs record (patch_flux.F90, H1)
     real(dp), allocatable :: bc_tinf(:)           ! (patch): T line flag 4: ambient temperature T_inf [K] (bc_val(6) = h)
     logical :: pole(2,3) = .false.                ! (side, axis): degenerate (pole)
                                                   !   face, from the .bcs_b2b file
     integer,  allocatable :: bc_eflag(:,:)        ! (3, patch): effective flag of the
                                                   !   transformed families U1,U2,U3 =
                                                   !   flag of the Cartesian axis each
                                                   !   grid direction points along
     real(dp), allocatable :: bc_val(:,:)
     ! b2b interface-face ownership masks: .true. where THIS block
     ! updates the shared normal-velocity dof on its min/max face
     logical, allocatable :: mu_lo(:,:), mu_hi(:,:)        ! (ny,nz) for u
     logical, allocatable :: mv_lo(:,:), mv_hi(:,:)        ! (nx,nz) for v
     logical, allocatable :: mw_lo(:,:), mw_hi(:,:)        ! (nx,ny) for w
     ! block-boundary face TYPE for the projection (same layouts):
     !   0 = block-to-block / periodic interface face (normal face)
     !   1 = physical boundary with a FIXED face flux (Dirichlet
     !       velocity, or zero-gradient velocity with Neumann p): the
     !       Poisson operator carries no pressure flux through it and
     !       the corrector leaves it at its star value
     !   2 = pressure outlet (Dirichlet p, zero-gradient velocity): the
     !       face is corrected by the projection like an interior face
     !   3 = mass-corrected outflow (velocity flag 2, Neumann p): fixed
     !       star flux like type 1, but the star fluxes of all type-3
     !       faces are scaled before the projection so that their total
     !       outflow equals the net inflow through the type-1 faces
     integer, allocatable :: tu_lo(:,:), tu_hi(:,:)
     integer, allocatable :: tv_lo(:,:), tv_hi(:,:)
     integer, allocatable :: tw_lo(:,:), tw_hi(:,:)
  end type mblk_t

  !------------------ one MPI sub-domain --------------------------------
  ! Field arrays use BLOCK-GLOBAL indices (r1-ng : r2+ng per axis).
  type :: mdom_t
     integer :: b  = 0                     ! parent block
     integer :: pc(3) = 0                  ! position in the block's split grid (0-based)
     integer :: r1(3) = 0, r2(3) = 0       ! owned cell range (block-global)
     integer :: owner = 0                  ! owning MPI rank
     ! ---- primary unknowns (curvilinear formulation) ------------------
     ! U1,U2,U3 : transformed (contravariant) velocity components
     !            U^q = u_vec . grad(xi_q) (times J for curvi_form 2),
     !            stored at the q-family faces; THESE are the variables
     !            coupled across sub-domain halos and block-to-block /
     !            periodic interfaces.
     ! p        : pressure at cell centers.
     ! alpha    : liquid volume fraction at cell centers (two-phase).
     real(dp), allocatable :: U1(:,:,:), U2(:,:,:), U3(:,:,:), p(:,:,:)
     real(dp), allocatable :: Us1(:,:,:), Us2(:,:,:), Us3(:,:,:)
     real(dp), allocatable :: alpha(:,:,:)
     ! ---- cavitation (cell centred, ghosted like alpha; update 15) ----
     !   vap   : vapour volume fraction of the liquid-vapour part
     !   vapn  : its time-level-n copy; cavsc / cavsp: explicit source and
     !           implicit sink coefficient of the vapour equation [1/tch];
     !   cavs  : volume source of the continuity equation (div u target)
     real(dp), allocatable :: vap(:,:,:), vapn(:,:,:), cavsc(:,:,:), cavsp(:,:,:), cavs(:,:,:), s_vap(:,:,:)
     real(dp), allocatable :: cavdp(:,:,:)                 ! dS/dp (<= 0): pressure-implicit part of the volume source (Poisson diagonal)
     real(dp), allocatable :: cavp0(:,:,:)                 ! pressure the sources were evaluated at (linearisation point)
     real(dp), allocatable :: s_cavs(:,:,:)                ! SSP copy of the continuity source
     real(dp), allocatable :: cavbc(:,:,:)                 ! compressibility coefficient beta_c/dt of the pressure equation (update 17)
     real(dp), allocatable :: presn(:,:,:)                 ! pressure at the start of the step (reference of the compressibility term)
     ! ---- energy (update 16): temperature [K] (cell centred, ghosted) --
     !   temp : T, tempn : its time-level-n copy, s_temp : SSP copy
     real(dp), allocatable :: temp(:,:,:), tempn(:,:,:), s_temp(:,:,:)
     ! ---- turbulence (cell centred, ghosted like alpha) ----------------
     real(dp), allocatable :: tk(:,:,:), tw(:,:,:)        ! k / nu_tilde ; omega / epsilon
     real(dp), allocatable :: mut(:,:,:)                  ! mu_t / mu1 (dimensionless eddy viscosity)
     real(dp), allocatable :: wdist(:,:,:)                ! wall distance [m]
     real(dp), allocatable :: wfac(:,:,:)                 ! wall-function factor mu_eff/mu at the wall face of the cell (1 = none)
     real(dp), allocatable :: alps(:,:,:), kap(:,:,:)     ! smoothed alpha and interface curvature (surface tension)
     ! ---- moving grid (rigid body motion, sixdof.F90) ------------------
     !   Ug1..3 : grid velocity in the units of U1..U3 (fvw * m . u_grid),
     !            subtracted from the fluid face components wherever a face
     !            flux advects something (ALE: relative flux); allocated for
     !            the sub-domains of moving blocks only
     !   Xh0, mh0 : half-index coordinates / metrics of the ORIGINAL grid
     logical :: moving = .false.
     real(dp), allocatable :: Ug1(:,:,:), Ug2(:,:,:), Ug3(:,:,:)
     real(dp), allocatable :: Xh0(:,:,:,:), mh0(:,:,:,:)
     integer,  allocatable :: wface(:,:,:)                ! wall face of the cell for the factor: 0 none, +-1..+-3
     ! ---- derived Cartesian face velocities (local, not exchanged) ----
     real(dp), allocatable :: u(:,:,:),  v(:,:,:),  w(:,:,:)
     real(dp), allocatable :: us(:,:,:), vs(:,:,:), ws(:,:,:)
     real(dp), allocatable :: rhsu(:,:,:), rhsv(:,:,:), rhsw(:,:,:), rhsp(:,:,:)
     real(dp), allocatable :: wrk(:,:,:)
     ! ---- time level n copies for the inner iterations ----------------
     real(dp), allocatable :: un(:,:,:), vn(:,:,:), wn(:,:,:)   ! Cartesian face vel.
     real(dp), allocatable :: alpn(:,:,:)                        ! volume fraction
     ! ---- 2nd-order SSP time scheme: time-level-n copies (Time_scheme = 2)
     real(dp), allocatable :: s_u(:,:,:), s_v(:,:,:), s_w(:,:,:), s_U1(:,:,:), s_U2(:,:,:), s_U3(:,:,:)
     real(dp), allocatable :: s_alpha(:,:,:), s_tk(:,:,:), s_tw(:,:,:), s_mut(:,:,:)
     real(dp), allocatable :: kpy(:,:,:), kpt(:,:,:)  ! line-SGS preconditioner work
     real(dp), allocatable :: mco(:,:,:,:,:)          ! momentum coefficients (26, i,j,k, c)
     ! ---- overset status (allocated when ovst_on) ---------------------
     !   cs  : cell status 1 active / 0 fringe (interpolated) / -1 hole
     !   fs1..3 : face status of the i/j/k families (same codes; 1 =
     !         solved face = at least one active neighbour cell)
     !   ovd : pressure-preconditioner diagonal compensation of the
     !         couplings dropped at fringe/hole cells (poisson.F90)
     integer,  allocatable :: cs(:,:,:), fs1(:,:,:), fs2(:,:,:), fs3(:,:,:)
     real(dp), allocatable :: ovd(:,:,:)
     integer :: qline = 2                  ! line direction of the SGS preconditioner
     ! ---- Krylov work vectors of the BiCGSTAB pressure solver ----------
     ! kv1/kv2 are exchanged (ghosts needed for the operator); the
     ! others are owned-cell only.
     real(dp), allocatable :: kv1(:,:,:), kv2(:,:,:)
     real(dp), allocatable :: kr(:,:,:), kr0(:,:,:), kpv(:,:,:)
     real(dp), allocatable :: kc(:,:,:,:)          ! 7-point face coefficients of the
                                                    !   pressure preconditioner (i,j,k, 1..3);
                                                    !   4..6: overset, unmasked (coupled operator)
     real(dp), allocatable :: krr(:,:,:)            ! overset: BiCGSTAB residual
     real(dp), allocatable :: kvv(:,:,:), ksv(:,:,:), kt(:,:,:)
     real(dp), allocatable :: kV(:,:,:,:), kZ(:,:,:,:)   ! FGMRES bases (owned cells)
     real(dp), allocatable :: kW(:,:,:,:)                 ! overset: inner GMRES basis (owned cells)
     ! ---- curvilinear metrics on the HALF-INDEX lattice ----------------
     ! Half index h: cell i center -> h=2i, face f -> h=2f+1 (node n ->
     ! h=2n-1).  Xh = coordinates at every half location; Jh = Jacobian
     ! det(dx/dxi); mh(c,q) = d(xi_q)/d(x_c) (Park's m matrix).
     integer :: hlo(3) = 0, hhi(3) = 0
     real(dp), allocatable :: Xh(:,:,:,:)          ! (3, h1,h2,h3)
     real(dp), allocatable :: Jh(:,:,:)
     real(dp), allocatable :: mh(:,:,:,:)          ! (9, h1,h2,h3): m(c,q)=mh(c+3(q-1))
  end type mdom_t

  integer :: nblk = 0                          ! number of grid blocks
  type(mblk_t), allocatable, target :: blk(:)  ! all blocks
  integer :: ndom = 0                          ! total number of sub-domains
  type(mdom_t), allocatable, target :: dom(:)  ! all sub-domains (metadata
                                               ! replicated; fields on owner)
  integer :: nown = 0                          ! # sub-domains owned by this rank
  integer, allocatable :: own_list(:)          ! their global domain ids

  !------------------ current domain view (set_block) ------------------
  integer :: cur_dom = 0, cur_blk = 0
  integer :: ni, nj, nk                    ! nodes of the current block
  integer :: nx, ny, nz                    ! cells of the current block
  integer :: is, ie, js, je, ks, ke        ! owned cell range of the domain
  real(dp), pointer :: u(:,:,:),  v(:,:,:),  w(:,:,:),  p(:,:,:)
  real(dp), pointer :: us(:,:,:), vs(:,:,:), ws(:,:,:)
  real(dp), pointer :: U1(:,:,:), U2(:,:,:), U3(:,:,:)
  real(dp), pointer :: Us1(:,:,:), Us2(:,:,:), Us3(:,:,:)
  real(dp), pointer :: alpha(:,:,:)
  real(dp), pointer :: vap(:,:,:), vapn(:,:,:), cavsc(:,:,:), cavsp(:,:,:), cavs(:,:,:), cavdp(:,:,:), cavp0(:,:,:)
  real(dp), pointer :: cavbc(:,:,:), presn(:,:,:)
  real(dp), pointer :: tk(:,:,:), tw(:,:,:), mut(:,:,:), wdist(:,:,:), wfac(:,:,:)
  real(dp), pointer :: temp(:,:,:), tempn(:,:,:)
  real(dp), pointer :: alps(:,:,:), kap(:,:,:)
  integer,  pointer :: wface(:,:,:)
  real(dp), pointer :: kv1(:,:,:), kv2(:,:,:), kr(:,:,:), kr0(:,:,:)
  real(dp), pointer :: kpv(:,:,:), kvv(:,:,:), ksv(:,:,:), kt(:,:,:)
  real(dp), pointer :: kV(:,:,:,:), kZ(:,:,:,:), kc(:,:,:,:), krr(:,:,:), kW(:,:,:,:)
  real(dp), pointer :: Xh(:,:,:,:), Jh(:,:,:), mh(:,:,:,:)
  real(dp), pointer :: Ug1(:,:,:), Ug2(:,:,:), Ug3(:,:,:)
  logical :: cur_moving = .false.          ! current domain belongs to a moving body
  logical :: cur_morph = .false.           ! current domain belongs to a morphing (flexible) block (update 24)
  integer :: hlo(3), hhi(3)
  real(dp), pointer :: rhsu(:,:,:), rhsv(:,:,:), rhsw(:,:,:), rhsp(:,:,:)
  real(dp), pointer :: wrk(:,:,:)
  real(dp), pointer :: un(:,:,:), vn(:,:,:), wn(:,:,:), alpn(:,:,:)
  real(dp), pointer :: kpy(:,:,:), kpt(:,:,:)
  integer :: nbc
  integer,  pointer :: bc_dir(:)
  integer,  pointer :: bc_irng(:,:), bc_jrng(:,:), bc_krng(:,:)
  integer,  pointer :: bc_flag(:,:), bc_eflag(:,:)
  real(dp), pointer :: bc_val(:,:), bc_tinf(:)
  logical,  pointer :: mu_lo(:,:), mu_hi(:,:)
  logical,  pointer :: mv_lo(:,:), mv_hi(:,:)
  logical,  pointer :: mw_lo(:,:), mw_hi(:,:)
  integer,  pointer :: tu_lo(:,:), tu_hi(:,:), tv_lo(:,:), tv_hi(:,:)
  integer,  pointer :: tw_lo(:,:), tw_hi(:,:)
  logical :: p_has_dirichlet               ! any Dirichlet pressure patch (all blocks)
  integer,  pointer :: ocs(:,:,:), ofs1(:,:,:), ofs2(:,:,:), ofs3(:,:,:)   ! overset status
  real(dp), pointer :: ovd(:,:,:)

  !------------------ b2b interface patches (block-level form) ---------
  ! Index maps are affine:  donor_index(q(a)) = sgn(a)*recv_index + off
  ! with off = coff(a) for cell-type indices, foff(a) for face-type.
  type :: xpatch_t
     integer :: rblk = 0, dblk = 0         ! receiver / donor block
     integer :: na = 0                     ! receiver face normal axis (1,2,3)
     integer :: rside = 0                  ! receiver side: 1 = min, 2 = max
     logical :: r_owner = .false.          ! receiver owns the shared face dof
     integer :: q(3) = 0                   ! donor axis for receiver axis a
     integer :: sgn(3) = 0                 ! direction sign
     integer :: coff(3) = 0                ! cell-index affine offsets
     integer :: foff(3) = 0                ! face-index affine offsets
     integer :: rn(2,3) = 0                ! receiver node window (from file)
  end type xpatch_t

  integer :: nxp = 0                       ! number of exchange patches
                                           !   (b2b + periodic)
  type(xpatch_t), allocatable :: xp(:)
  logical :: b2b_permuted = .false.        ! any non-identity axis map

  !------------------ periodic boundary patches (flag 3) ---------------
  type :: pbc_t
     integer :: rblk = 0, dirs = 0         ! local block, signed face dir
     integer :: rn(2,3) = 0                ! local surface node ranges
     integer :: dblk = 0, ddir = 0         ! pair block, signed face dir
     integer :: dn(2,3) = 0                ! pair surface node ranges
     integer :: q(3) = 0, sgn(3) = 0       ! resolved axis map
     real(dp) :: tv(3) = 0.0_dp            ! node-to-node translation [m]
  end type pbc_t
  integer :: nper = 0
  type(pbc_t), allocatable :: pbc(:)

  !------------------ unified exchange records (domain level) ----------
  type :: xrec_t
     integer :: kind = 0                   ! 1 = b2b patch, 2 = intra-block halo
     integer :: rd = 0, dd = 0             ! receiver / donor domain id
     integer :: ip = 0                     ! parent b2b patch (kind 1)
     integer :: na = 0                     ! receiver normal axis (pass group)
     integer :: nrm1 = 0, nrm2 = 0         ! kind 2: ghost index layers
     integer :: ds(2,3) = 0                ! donor supply (donor axis order)
  end type xrec_t

  integer :: nxr = 0
  type(xrec_t), allocatable :: xr(:)

contains

  !---------------------------------------------------------------------
  ! Point the module-level view at sub-domain d.
  !---------------------------------------------------------------------
  subroutine set_block(d)
    integer, intent(in) :: d
    integer :: b
    b = dom(d)%b
    cur_dom = d;  cur_blk = b
    ni = blk(b)%ni;  nj = blk(b)%nj;  nk = blk(b)%nk
    nx = blk(b)%nx;  ny = blk(b)%ny;  nz = blk(b)%nz
    is = dom(d)%r1(1);  ie = dom(d)%r2(1)
    js = dom(d)%r1(2);  je = dom(d)%r2(2)
    ks = dom(d)%r1(3);  ke = dom(d)%r2(3)
    u  => dom(d)%u;   v  => dom(d)%v;   w  => dom(d)%w;   p => dom(d)%p
    us => dom(d)%us;  vs => dom(d)%vs;  ws => dom(d)%ws
    U1 => dom(d)%U1;  U2 => dom(d)%U2;  U3 => dom(d)%U3
    Us1 => dom(d)%Us1;  Us2 => dom(d)%Us2;  Us3 => dom(d)%Us3
    alpha => dom(d)%alpha
    if (allocated(dom(d)%vap)) then
       vap => dom(d)%vap;  vapn => dom(d)%vapn;  cavsc => dom(d)%cavsc;  cavsp => dom(d)%cavsp;  cavs => dom(d)%cavs
       cavdp => dom(d)%cavdp;  cavp0 => dom(d)%cavp0;  cavbc => dom(d)%cavbc;  presn => dom(d)%presn
    end if
    if (allocated(dom(d)%tk)) then
       tk => dom(d)%tk;  tw => dom(d)%tw;  mut => dom(d)%mut;  wdist => dom(d)%wdist;  wfac => dom(d)%wfac;  wface => dom(d)%wface
    end if
    if (allocated(dom(d)%temp)) then
       temp => dom(d)%temp;  tempn => dom(d)%tempn
    end if
    if (allocated(dom(d)%alps)) then
       alps => dom(d)%alps;  kap => dom(d)%kap
    end if
    kv1 => dom(d)%kv1;  kv2 => dom(d)%kv2
    kr => dom(d)%kr;  kr0 => dom(d)%kr0;  kpv => dom(d)%kpv
    kvv => dom(d)%kvv;  ksv => dom(d)%ksv;  kt => dom(d)%kt
    kV => dom(d)%kV;  kZ => dom(d)%kZ
    if (allocated(dom(d)%kc)) kc => dom(d)%kc
    if (allocated(dom(d)%krr)) krr => dom(d)%krr
    if (allocated(dom(d)%kW))  kW  => dom(d)%kW
    Xh => dom(d)%Xh;  Jh => dom(d)%Jh;  mh => dom(d)%mh
    hlo = dom(d)%hlo;  hhi = dom(d)%hhi
    cur_moving = dom(d)%moving
    cur_morph = .false.
    if (fsi_on .and. allocated(fsi_blk)) cur_morph = fsi_blk(b)      ! morphing block (update 24)
    if (cur_moving) then
       Ug1 => dom(d)%Ug1;  Ug2 => dom(d)%Ug2;  Ug3 => dom(d)%Ug3
    else
       Ug1 => null();  Ug2 => null();  Ug3 => null()
    end if
    rhsu => dom(d)%rhsu;  rhsv => dom(d)%rhsv
    rhsw => dom(d)%rhsw;  rhsp => dom(d)%rhsp
    wrk  => dom(d)%wrk
    un => dom(d)%un;  vn => dom(d)%vn;  wn => dom(d)%wn;  alpn => dom(d)%alpn
    kpy => dom(d)%kpy;  kpt => dom(d)%kpt
    nbc = blk(b)%nbc
    bc_dir => blk(b)%bc_dir
    bc_irng => blk(b)%bc_irng;  bc_jrng => blk(b)%bc_jrng
    bc_krng => blk(b)%bc_krng
    bc_flag => blk(b)%bc_flag;  bc_val => blk(b)%bc_val;  bc_tinf => blk(b)%bc_tinf
    bc_eflag => blk(b)%bc_eflag
    mu_lo => blk(b)%mu_lo;  mu_hi => blk(b)%mu_hi
    mv_lo => blk(b)%mv_lo;  mv_hi => blk(b)%mv_hi
    mw_lo => blk(b)%mw_lo;  mw_hi => blk(b)%mw_hi
    tu_lo => blk(b)%tu_lo;  tu_hi => blk(b)%tu_hi
    tv_lo => blk(b)%tv_lo;  tv_hi => blk(b)%tv_hi
    tw_lo => blk(b)%tw_lo;  tw_hi => blk(b)%tw_hi
    if (ovst_on) then
       ocs => dom(d)%cs;  ofs1 => dom(d)%fs1;  ofs2 => dom(d)%fs2;  ofs3 => dom(d)%fs3
       ovd => dom(d)%ovd
    end if
  end subroutine set_block

  !---------------------------------------------------------------------
  ! Grid-velocity face component of family q at array position (i,j,k)
  ! of the current domain, in the units of U1..U3 (zero on static
  ! blocks).  The ADVECTING flux of any transported quantity is the
  ! relative flux  U^q - Ug^q  (ALE form on the rigidly moving grids).
  !---------------------------------------------------------------------
  pure real(dp) function ugq(q, i, j, k)
    integer, intent(in) :: q, i, j, k
    ugq = 0.0_dp
    if (.not. cur_moving) return
    select case (q)
    case (1);  ugq = Ug1(i,j,k)
    case (2);  ugq = Ug2(i,j,k)
    case default;  ugq = Ug3(i,j,k)
    end select
  end function ugq

  !---------------------------------------------------------------------
  ! Overset status helpers on the current domain (1 = solved/active
  ! everywhere when the overset assembly is off).
  !---------------------------------------------------------------------
  pure integer function fstat(q, i, j, k)
    integer, intent(in) :: q, i, j, k
    if (.not. ovst_on) then
       fstat = 1;  return
    end if
    select case (q)
    case (1);  fstat = ofs1(i,j,k)
    case (2);  fstat = ofs2(i,j,k)
    case default;  fstat = ofs3(i,j,k)
    end select
  end function fstat

  pure integer function cstat(i, j, k)
    integer, intent(in) :: i, j, k
    if (.not. ovst_on) then
       cstat = 1;  return
    end if
    cstat = ocs(i,j,k)
  end function cstat

  !---------------------------------------------------------------------
  ! Projection type of the family-q face with face index fq and the
  ! two other (cell) indices taken from pt: 0 interior/interface,
  ! 1 fixed-flux physical boundary, 2 pressure outlet (see mblk_t).
  !---------------------------------------------------------------------
  pure integer function face_type(q, fq, pt)
    integer, intent(in) :: q, fq, pt(3)
    face_type = 0
    select case (q)
    case (1)
       if (fq == 0) then
          face_type = tu_lo(pt(2),pt(3))
       else if (fq == nx) then
          face_type = tu_hi(pt(2),pt(3))
       end if
    case (2)
       if (fq == 0) then
          face_type = tv_lo(pt(1),pt(3))
       else if (fq == ny) then
          face_type = tv_hi(pt(1),pt(3))
       end if
    case default
       if (fq == 0) then
          face_type = tw_lo(pt(1),pt(2))
       else if (fq == nz) then
          face_type = tw_hi(pt(1),pt(2))
       end if
    end select
  end function face_type

  !---------------------------------------------------------------------
  ! Pointer to the field array 'iv' of sub-domain d (see the id table
  ! in the header).
  !---------------------------------------------------------------------
  function field_ptr(d, iv) result(f)
    integer, intent(in) :: d, iv
    real(dp), pointer   :: f(:,:,:)
    select case (iv)
    case (1);  f => dom(d)%U1             ! transformed components are the
    case (2);  f => dom(d)%U2             ! exchanged face unknowns (a face
    case (3);  f => dom(d)%U3             ! dof is geometric: 1:1 across
    case (4);  f => dom(d)%p              ! interfaces of any orientation)
    case (5);  f => dom(d)%Us1
    case (6);  f => dom(d)%Us2
    case (7);  f => dom(d)%Us3
    case (8);  f => dom(d)%alpha
    case (9);  f => dom(d)%kv1
    case (10); f => dom(d)%kv2
    case (11); f => dom(d)%rhsu
    case (12); f => dom(d)%rhsv
    case (13); f => dom(d)%rhsw
    case (14); f => dom(d)%tk
    case (15); f => dom(d)%tw
    case (16); f => dom(d)%mut
    case (17); f => dom(d)%alps
    case (18); f => dom(d)%kap
    case (19); f => dom(d)%vap             ! vapour volume fraction (cavitation)
    case (20); f => dom(d)%temp            ! temperature [K] (energy equation, update 16)
    case default;  f => null()
    end select
  end function field_ptr

  !---------------------------------------------------------------------
  ! Component axis of field id iv (0 for cell-centered fields).
  !---------------------------------------------------------------------
  pure integer function comp_axis(iv)
    integer, intent(in) :: iv
    if (iv == 4 .or. (iv >= 8 .and. iv <= 10) .or. iv >= 14) then
       comp_axis = 0                        ! cell-centred: p, alpha, Krylov, turbulence (14-16), alps / kappa (17-18), vapour (19), T (20)
    else if (iv <= 3) then
       comp_axis = iv
    else if (iv >= 11) then
       comp_axis = iv - 10
    else
       comp_axis = iv - 4
    end if
  end function comp_axis

  !---------------------------------------------------------------------
  ! Cell count of block b along axis a.
  !---------------------------------------------------------------------
  pure integer function ncell(b, a)
    integer, intent(in) :: b, a
    select case (a)
    case (1);  ncell = blk(b)%nx
    case (2);  ncell = blk(b)%ny
    case default;  ncell = blk(b)%nz
    end select
  end function ncell

  !---------------------------------------------------------------------
  ! Curvilinear metric accessors on the current domain's half lattice.
  ! Half index: cell i -> 2i, face f -> 2f+1.  All are O(1) lookups /
  ! small combinations of the stored mh, Jh, Xh.
  !---------------------------------------------------------------------
  ! face-variable weight: the stored face unknown is  fvw * (m . u)
  !   curvi_form = 1 -> 1     (Park's unweighted component U^q)
  !   curvi_form = 2 -> J     (volume flux F^q = J U^q)
  pure real(dp) function fvw(h1, h2, h3)
    integer, intent(in) :: h1, h2, h3
    if (curvi_form == 2) then
       fvw = Jh(h1, h2, h3)
       if (fvw == 0.0_dp) fvw = 1.0_dp       ! pole faces (zero Jacobian)
    else
       fvw = 1.0_dp
    end if
  end function fvw

  pure real(dp) function mq(c, q, h1, h2, h3)      ! d(xi_q)/d(x_c)
    integer, intent(in) :: c, q, h1, h2, h3
    mq = mh(c + 3*(q-1), h1, h2, h3)
  end function mq

  pure real(dp) function gql(q, l, h1, h2, h3)     ! grad(xi_q).grad(xi_l)
    integer, intent(in) :: q, l, h1, h2, h3
    gql = mh(1+3*(q-1),h1,h2,h3)*mh(1+3*(l-1),h1,h2,h3)                  &
        + mh(2+3*(q-1),h1,h2,h3)*mh(2+3*(l-1),h1,h2,h3)                  &
        + mh(3+3*(q-1),h1,h2,h3)*mh(3+3*(l-1),h1,h2,h3)
  end function gql

  pure real(dp) function Gcap(q, l, h1, h2, h3)    ! J * g^{ql}
    integer, intent(in) :: q, l, h1, h2, h3
    Gcap = Jh(h1,h2,h3) * gql(q, l, h1, h2, h3)
  end function Gcap

  ! forward tangent d(x_c)/d(xi_q) = inverse of the m matrix at the half
  ! position (T = M^-1, T(c,q) = J * cofactor of M(q,c)).  Where m was
  ! built from the coordinate differences this IS the central difference
  ! X(h+e_q) - X(h-e_q) (round-off apart); where a degenerate ghost /
  ! wedge-edge metric was replaced by the nearest interior one (curvi.F90)
  ! it is the inverse of THAT metric, so the round trip Cartesian ->
  ! transformed -> Cartesian of the ghost / boundary faces is the identity.
  ! (Until update 12 tq was the raw coordinate difference: at replaced
  ! positions m and tq were not inverses and the momentum iteration
  ! amplified the ghost velocities of the wall edges of the Box_6DoF
  ! grids by ~20 per sweep.)
  pure real(dp) function tq(c, q, h1, h2, h3)
    integer, intent(in) :: c, q, h1, h2, h3
    real(dp) :: a(3,3), cofq(3), det
    integer :: qq, cc, q1, q2, c1, c2
    ! a(q,c) = m(c,q) = d(xi_q)/d(x_c)
    do qq = 1, 3
       do cc = 1, 3
          a(qq,cc) = mh(cc + 3*(qq-1), h1, h2, h3)
       end do
    end do
    ! cofactors of column c: cof(q',c) for q' = 1..3
    c1 = mod(c, 3) + 1;  c2 = mod(c+1, 3) + 1
    do qq = 1, 3
       q1 = mod(qq, 3) + 1;  q2 = mod(qq+1, 3) + 1
       cofq(qq) = a(q1,c1)*a(q2,c2) - a(q1,c2)*a(q2,c1)
    end do
    det = a(1,c)*cofq(1) + a(2,c)*cofq(2) + a(3,c)*cofq(3)
    if (abs(det) > 1.0e-300_dp) then
       tq = cofq(q) / det                 ! T(c,q) = cof(M)(q,c) / det(M)
    else
       select case (q)                     ! pole / zero metric: raw difference
       case (1);  tq = Xh(c, h1+1, h2, h3) - Xh(c, h1-1, h2, h3)
       case (2);  tq = Xh(c, h1, h2+1, h3) - Xh(c, h1, h2-1, h3)
       case default;  tq = Xh(c, h1, h2, h3+1) - Xh(c, h1, h2, h3-1)
       end select
    end if
  end function tq

  !---------------------------------------------------------------------
  ! Fluid properties on the half lattice (two-phase, one-fluid model).
  ! The liquid volume fraction alpha lives at cell centers (even half
  ! indices); at any half location it is the arithmetic average of the
  ! 1/2/4/8 cells touching that location.  Non-dimensional density and
  ! viscosity are linear in alpha:
  !    rho* = alpha + rho_ratio (1 - alpha),  mu* = alpha + mu_ratio (1 - alpha)
  ! For Num_of_Phase = 1 all three functions return exactly 1.
  !---------------------------------------------------------------------
  pure real(dp) function alph_h(h1, h2, h3)
    integer, intent(in) :: h1, h2, h3
    integer :: c1a, c1b, c2a, c2b, c3a, c3b
    real(dp) :: s
    if (Num_of_Phase == 1) then
       alph_h = 1.0_dp
       return
    end if
    if (mod(h1,2) == 0) then
       c1a = h1/2;  c1b = c1a
    else
       c1a = (h1-1)/2;  c1b = c1a + 1
    end if
    if (mod(h2,2) == 0) then
       c2a = h2/2;  c2b = c2a
    else
       c2a = (h2-1)/2;  c2b = c2a + 1
    end if
    if (mod(h3,2) == 0) then
       c3a = h3/2;  c3b = c3a
    else
       c3a = (h3-1)/2;  c3b = c3a + 1
    end if
    s = alpha(c1a,c2a,c3a)
    if (c1b /= c1a) s = s + alpha(c1b,c2a,c3a)
    if (c2b /= c2a) then
       s = s + alpha(c1a,c2b,c3a)
       if (c1b /= c1a) s = s + alpha(c1b,c2b,c3a)
    end if
    if (c3b /= c3a) then
       s = s + alpha(c1a,c2a,c3b)
       if (c1b /= c1a) s = s + alpha(c1b,c2a,c3b)
       if (c2b /= c2a) then
          s = s + alpha(c1a,c2b,c3b)
          if (c1b /= c1a) s = s + alpha(c1b,c2b,c3b)
       end if
    end if
    alph_h = s / real((c1b-c1a+1)*(c2b-c2a+1)*(c3b-c3a+1), dp)
  end function alph_h

  ! vapour fraction at a half location (average of the touching cells)
  pure real(dp) function vap_h(h1, h2, h3)
    integer, intent(in) :: h1, h2, h3
    integer :: c1a, c1b, c2a, c2b, c3a, c3b, i, j, k
    real(dp) :: s
    c1a = h1/2;  c1b = (h1+1)/2;  c2a = h2/2;  c2b = (h2+1)/2;  c3a = h3/2;  c3b = (h3+1)/2
    c1a = max(c1a, lbound(vap,1));  c1b = min(c1b, ubound(vap,1))
    c2a = max(c2a, lbound(vap,2));  c2b = min(c2b, ubound(vap,2))
    c3a = max(c3a, lbound(vap,3));  c3b = min(c3b, ubound(vap,3))
    s = 0.0_dp
    do k = c3a, c3b
       do j = c2a, c2b
          do i = c1a, c1b
             s = s + vap(i,j,k)
          end do
       end do
    end do
    vap_h = s / real((c1b-c1a+1)*(c2b-c2a+1)*(c3b-c3a+1), dp)
  end function vap_h

  ! three-phase (update 15): rho* = alpha [(1 - av) + av rho_v/rho1] + (1 - alpha) rho2/rho1
  pure real(dp) function rho_h(h1, h2, h3)         ! rho / rho1
    integer, intent(in) :: h1, h2, h3
    real(dp) :: a, av
    a = 1.0_dp
    if (Num_of_Phase == 2) a = alph_h(h1, h2, h3)
    if (cav_on) then
       av = vap_h(h1, h2, h3)
       rho_h = a*((1.0_dp - av) + av*rv_ratio) + rho_ratio*(1.0_dp - a)
    else if (Num_of_Phase == 1) then
       rho_h = 1.0_dp
    else
       rho_h = a + rho_ratio*(1.0_dp - a)
    end if
  end function rho_h

  pure real(dp) function mu_h(h1, h2, h3)          ! (mu + mu_t) / mu1
    integer, intent(in) :: h1, h2, h3
    real(dp) :: a, av
    integer  :: c1a, c1b, c2a, c2b, c3a, c3b, i, j, k
    a = 1.0_dp
    if (Num_of_Phase == 2) a = alph_h(h1, h2, h3)
    if (cav_on) then
       av = vap_h(h1, h2, h3)
       mu_h = a*((1.0_dp - av) + av*muv_ratio) + mu_ratio*(1.0_dp - a)
    else if (Num_of_Phase == 1) then
       mu_h = 1.0_dp
    else
       mu_h = a + mu_ratio*(1.0_dp - a)
    end if
    if (turb_on) then
       ! eddy viscosity: average of the cells sharing the half-index position
       c1a = h1/2;  c1b = (h1+1)/2;  c2a = h2/2;  c2b = (h2+1)/2;  c3a = h3/2;  c3b = (h3+1)/2
       c1a = max(c1a, lbound(mut,1));  c1b = min(c1b, ubound(mut,1))
       c2a = max(c2a, lbound(mut,2));  c2b = min(c2b, ubound(mut,2))
       c3a = max(c3a, lbound(mut,3));  c3b = min(c3b, ubound(mut,3))
       mu_h = mu_h + sum(mut(c1a:c1b, c2a:c2b, c3a:c3b)) / real((c1b-c1a+1)*(c2b-c2a+1)*(c3b-c3a+1), dp)
    end if
    if (wallfn_on) then
       ! wall-function faces: the viscous flux through the wall face of a
       ! flag-11 cell is mu_eff = wfac * mu (laminar), so that
       ! mu_eff (u1 - 0)/y1 = tau_w of the wall function
       i = (h1+1)/2;  j = (h2+1)/2;  k = (h3+1)/2
       i = min(max(i, lbound(wface,1)), ubound(wface,1))
       j = min(max(j, lbound(wface,2)), ubound(wface,2))
       k = min(max(k, lbound(wface,3)), ubound(wface,3))
       if (wface(i,j,k) /= 0) then
          select case (abs(wface(i,j,k)))
          case (1);  if (h1 == 2*i - merge(1, -1, wface(i,j,k) > 0)) mu_h = wfac(i,j,k)
          case (2);  if (h2 == 2*j - merge(1, -1, wface(i,j,k) > 0)) mu_h = wfac(i,j,k)
          case (3);  if (h3 == 2*k - merge(1, -1, wface(i,j,k) > 0)) mu_h = wfac(i,j,k)
          end select
       end if
    end if
  end function mu_h

  !---------------------------------------------------------------------
  ! Donor supply range of domain d along block axis a: the cell/face
  ! indices this domain provides to others = its owned range, extended
  ! by ng past block boundaries it touches (block-ghost storage).
  !---------------------------------------------------------------------
  subroutine dom_supply(d, a, s1, s2)
    integer, intent(in)  :: d, a
    integer, intent(out) :: s1, s2
    s1 = dom(d)%r1(a);  s2 = dom(d)%r2(a)
    if (s1 == 1)                 s1 = 1 - ng
    if (s2 == ncell(dom(d)%b,a)) s2 = s2 + ng
  end subroutine dom_supply

end module global_vars
