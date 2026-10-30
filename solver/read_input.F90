!=======================================================================
! read_input.F90
!   Input routines for probCFD (multi-block version):
!     read_namelist   : probCFD.nml
!     read_split_csv  : Input_polyMesh_probCFD.csv (per-block MPI split)
!     read_grid       : casename.grd  (PLOT3D, ASCII, MULTI-block)
!     read_bcs        : casename.bcs  (boundary conditions, per block)
!=======================================================================

!-----------------------------------------------------------------------
! read_namelist : probCFD.nml in the GROUPED format
!
!   $Group_name
!    Variable = value(s)        ! comment
!   /end
!
!   Groups: $Case_title, $Control, $Fluid_Properties, $Reference_Pressure,
!   $Cell_Initial_Condition, $Solver_Iteration, $Export_Interval and the
!   probCFD-specific $Numerical_Scheme (discretization / solver switches).
!   Variable names are case-insensitive; everything after "!" is a
!   comment; values are separated by commas or blanks; strings may be
!   quoted.  Unknown variable names abort with a message (misspelled
!   names are never silently ignored).  Every variable has a default
!   (listed below); a group may be omitted entirely.
!-----------------------------------------------------------------------
subroutine read_namelist
  use global_vars
  use mpi_vars
  use sixdof, only: sixdof_nml, Moving_BGR, Moving_BGR_DoF, BGR_Following_Body, Moving_BGR_Filter_T
  use waves, only: wave_nml
  use force_model, only: force_nml, force_in_group
  use stop_criteria, only: stop_report_input, stop_clean
  implicit none
  logical :: ok6, ex
  integer :: iu, ios, ierr, lnum, pe, nv, io, i
  character(len=512) :: line, work, name, group
  character(len=256) :: vals(32)      ! up to 32 values per line (update 19: Prop_* station lines)
  integer  :: nml_nsteps

  ! ---- defaults --------------------------------------------------------
  casename    = 'Case1'
  Lch         = 1.0_dp
  Uch         = 1.0_dp
  dt_dim      = 1.0e-3_dp
  End_T       = 0.1_dp
  z_free_surface = 0.0_dp
  grav        = 0.0_dp
  Num_of_Phase= 1
  rho1        = 1000.0_dp
  mu1         = 1.0e-3_dp
  rho2        = 1.2_dp
  mu2         = 1.8e-5_dp
  rho         = -1.0_dp           ! legacy aliases (< 0: unused)
  mu          = -1.0_dp
  p_ref_switch= .false.
  p_ref_x     = 0.0_dp
  p_ref_y     = 0.0_dp
  p_ref_z     = 0.0_dp
  p_ref_val   = 0.0_dp
  u_init      = 0.0_dp
  v_init      = 0.0_dp
  w_init      = 0.0_dp
  alpha_init  = 0.0_dp
  p_init      = 0.0_dp
  Iner_Itn    = 1
  mom_maxiter = 1000
  mom_miniter = 1
  mom_tol     = 1.0e-9_dp
  mom_sor     = 1.0_dp
  poi_maxiter = 500
  poi_miniter = 1
  poi_tol     = 1.0e-7_dp
  poi_rtol    = 1.0e-6_dp
  poi_sor     = 1.0_dp
  vof_maxiter = 1000
  vof_miniter = 1
  vof_tol     = 1.0e-9_dp
  vof_sor     = 1.0_dp
  out_interval= 100
  log_interval= 1
  ! numerical scheme switches
  Cnvt_scheme = 1
  curvi_form  = 1
  Cnvt_mode   = 3
  Dffs_mode   = 3
  cfl_switch  = 1.0_dp
  vn_switch   = 0.5_dp
  poi_solver  = 2
  poi_precond = 2
  vof_scheme  = 2
  vof_cfl_max = 0.5_dp
  t_ramp      = 0.0_dp
  b2b_tol     = 1.0e-8_dp
  ovst_switch = -1
  surfovl_tol = 1.0e-4_dp
  ovst_nf     = 2
  ovst_wall_layers = 8
  ovst_active_overlap = 1
  ovst_vtol   = 1.0e-6_dp
  ovst_only   = 0

  inquire(file='probCFD.nml', exist=ex)
  if (.not. ex) call stop_clean(5, 'There is no name list file (probCFD.nml must be in the working folder).')
  open(newunit=iu, file='probCFD.nml', status='old', action='read', iostat=ios)
  if (ios /= 0) then
     if (myrank == 0) write(*,*) ' ERROR: cannot open probCFD.nml'
     call MPI_Abort(MPI_COMM_WORLD, 10, ierr)
  end if

  group = ' ';  lnum = 0
  do
     read(iu, '(A)', iostat=ios) line
     if (ios /= 0) exit
     lnum = lnum + 1
     work = line
     pe = index(work, char(13));  if (pe > 0) work(pe:) = ' '     ! CR (Windows)
     pe = index(work, '!');       if (pe > 0) work(pe:) = ' '     ! comment
     work = adjustl(work)
     if (len_trim(work) == 0) cycle
     if (work(1:1) == '$') then                                 ! group start
        group = lower(trim(work(2:)))
        cycle
     end if
     if (lower(work(1:4)) == '/end' .or. work(1:1) == '/') then ! group end
        group = ' '
        cycle
     end if
     pe = index(work, '=')
     if (pe == 0) call nml_err_line('expected "name = value"')
     name = lower(trim(adjustl(work(1:pe-1))))
     call split_values(work(pe+1:), vals, nv)
     if (nv == 0) call nml_err_line('no value given for '//trim(name))
     call assign(trim(name))
  end do
  close(iu)

  ! ---- legacy aliases and validation ------------------------------------
  if (rho > 0.0_dp) rho1 = rho
  if (mu  > 0.0_dp) mu1  = mu
  call chk_int(curvi_form,   1, 2, 'curvi_form')
  call chk_int(Cnvt_scheme,  1, 4, 'Cnvt_scheme')
  call chk_int(Cnvt_mode,    1, 3, 'Cnvt_mode')
  call chk_int(Dffs_mode,    1, 3, 'Dffs_mode')
  call chk_int(Prs_form,     1, 2, 'Prs_form')
  call chk_int(quasi2d,      0, 3, 'SW_Quasi_2D')
  if (quasi2d == 1 .or. quasi2d == 2) then
     call nml_err('SW_Quasi_2D: only 3 (the span along k) is supported')
  end if
  call chk_int(Num_of_Phase, 1, 2, 'Num_of_Phase')
  call chk_int(poi_solver,   1, 2, 'poi_solver')
  call chk_int(poi_precond,  1, 2, 'poi_precond')
  call chk_int(vof_scheme,   1, 3, 'vof_scheme')
  call chk_int(Surf_tension, 0, 1, 'Surf_tension')
  call chk_int(Init_shape, 0, 3, 'Init_shape')
  if (Sigma_ST < 0.0_dp) call nml_err('Sigma_ST must be >= 0')
  st_on = (Surf_tension == 1 .and. Num_of_Phase == 2)
  sigma_nd = 0.0_dp
  if (st_on) sigma_nd = Sigma_ST/(rho1*Uch*Uch*Lch)
  call chk_int(TURB_model,   0, 4, 'TURB_model')
  call chk_int(IDDES_model,  0, 1, 'IDDES_model')
  call chk_int(Cnvt_scheme_Turb, 1, 2, 'Cnvt_scheme_Turb')
  call chk_int(Cnvt_mode_Turb,   1, 2, 'Cnvt_mode_Turb')
  call chk_int(Dffs_mode_Turb,   1, 2, 'Dffs_mode_Turb')
  if (Max_Itr_Turb < 1 .or. Min_Itr_Turb < 0 .or. Rsd_Lmt_Turb <= 0.0_dp .or. SOR_fct_Turb <= 0.0_dp &
      .or. cfl_switch_Turb <= 0.0_dp .or. vn_switch_Turb <= 0.0_dp) call nml_err('turbulence solver controls must be positive')
  turb_on = (TURB_model > 0 .or. IDDES_model > 0)
  if (Stop_Rsd_Max <= 0.0_dp .or. Stop_CFL_Max <= 0.0_dp) call nml_err('Stop_Rsd_Max and Stop_CFL_Max must be positive')
  call chk_int(Iner_Itn,     1, 100, 'Iner_Itn')
  call chk_int(Time_scheme,  1, 2, 'Time_scheme')
  call chk_int(run_mode,     0, 2, 'Mode')
  inquire(file=trim(casename)//'.grd', exist=ex)
  if (.not. ex) call stop_clean(5, 'There is no grid file ('//trim(casename)//'.grd must be in the working folder).')
  if (restart_interval < 0) call nml_err('Restart_Intvl must be >= 0')
  if (Lch <= 0.0_dp .or. Uch <= 0.0_dp) call nml_err('Char_L and Char_U must be positive')
  if (dt_dim <= 0.0_dp) call nml_err('Del_T must be positive')
  if (End_T <= 0.0_dp) call nml_err('End_T must be positive')
  if (rho1 <= 0.0_dp .or. mu1 <= 0.0_dp) call nml_err('Lou_w/mu_w must be positive')
  if (Num_of_Phase == 2 .and. (rho2 <= 0.0_dp .or. mu2 <= 0.0_dp))       &
     call nml_err('Lou_a/mu_a must be positive for Num_of_Phase = 2')
  if (cfl_switch <= 0.0_dp .or. vn_switch <= 0.0_dp)                     &
     call nml_err('cfl_switch and vn_switch must be positive')
  ! explicit 2nd-order upwind is stable only up to CFL 1/2 (1st-order:
  ! 1), and the explicit directions add up in multi-D
  if (myrank == 0) then
     if (Cnvt_scheme >= 2 .and. Cnvt_mode == 3 .and. cfl_switch > 0.5_dp) then
        write(*,'(A,I0,A)') '  WARNING: Cnvt_scheme = ', Cnvt_scheme, ' with Cnvt_mode = 3 and cfl_switch > 0.5:'
        write(*,'(A)') '           the explicit higher-order schemes are stable only below CFL 1/2'
        write(*,'(A)') '           per direction (the explicit directions add up); use'
        write(*,'(A)') '           cfl_switch <= 0.3 or Cnvt_mode = 2.'
     end if
     if (Cnvt_scheme == 3 .and. Cnvt_mode /= 2) then
        write(*,'(A)') '  WARNING: QUICK (Cnvt_scheme = 3) with explicit forward-Euler directions'
        write(*,'(A)') '           is unstable for pure convection (only viscosity damps it):'
        write(*,'(A)') '           Cnvt_mode = 2 (implicit) with Iner_Itn = 2 is recommended.'
     end if
     if (Cnvt_mode == 1) write(*,'(A)')                                    &
        '  NOTE: Cnvt_mode = 1 (explicit convection): keep the total CFL below 1 (1st order) / 0.5 (higher order)'
  end if
  if (vof_cfl_max <= 0.0_dp .or. vof_cfl_max > 1.0_dp)                   &
     call nml_err('vof_cfl_max must lie in (0, 1]')
  if (mom_sor <= 0.0_dp .or. mom_sor > 1.0_dp) call nml_err('SOR_fct_Vel must lie in (0, 1]')
  if (poi_sor <= 0.0_dp .or. poi_sor > 1.0_dp) call nml_err('SOR_fct_Prs must lie in (0, 1]')
  if (vof_sor <= 0.0_dp .or. vof_sor > 1.0_dp) call nml_err('SOR_fct_VoF must lie in (0, 1]')
  if (mom_maxiter < 1 .or. poi_maxiter < 1 .or. vof_maxiter < 1)          &
     call nml_err('Max_Itr_Vel/Prs/VoF must be >= 1')
  mom_miniter = max(1, min(mom_miniter, mom_maxiter))
  poi_miniter = max(1, min(poi_miniter, poi_maxiter))
  vof_miniter = max(1, min(vof_miniter, vof_maxiter))
  if (out_interval < 1) call nml_err('Output_Intvl must be >= 1')
  if (data_interval < 1) call nml_err('Data_Intvl must be >= 1')
  if (grid_scale <= 0.0_dp) call nml_err('Grid_scale must be positive')
  if (log_interval < 1) log_interval = 1
  if (alpha_init < 0.0_dp .or. alpha_init > 1.0_dp) call nml_err('a_cel_initial must lie in [0, 1]')
  call chk_int(ovst_switch, -1, 1, 'SW_Ovst')
  if (ovst_nf < 2 .or. ovst_nf > ng) call nml_err('Ovst_fringe must be 2 .. ng (= 3)')
  if (surfovl_tol < 0.0_dp) call nml_err('Surfovl_tol must be >= 0')
  if (ovst_vtol < 0.0_dp .or. ovst_vtol >= 0.5_dp) call nml_err('Ovst_vol_tol must lie in [0, 0.5)')
  if (ovst_reasm_tol < 0.0_dp) call nml_err('Ovst_Reassemble_Tol must be >= 0')
  call chk_int(Moving_BGR, 0, 1, 'Moving_BGR')
  if (Moving_BGR == 1 .and. BGR_Following_Body < 1) call nml_err('Moving_BGR = 1 needs BGR_Following_Body = <Body id>')

  ! ---- derived quantities -----------------------------------------------
  nml_nsteps = nint(End_T / dt_dim)
  if (nml_nsteps < 1) nml_nsteps = 1
  nstep = nml_nsteps
  tch = Lch / Uch
  Re  = rho1 * Uch * Lch / mu1
  dt  = dt_dim / tch
  if (Num_of_Phase == 2) then
     rho_ratio = rho2 / rho1
     mu_ratio  = mu2  / mu1
  else
     rho_ratio = 1.0_dp
     mu_ratio  = 1.0_dp
  end if
  gnd = grav * Lch / (Uch*Uch)
  zfs = z_free_surface / Lch
  ! cavitation (update 15): vapour property ratios, non-dimensional vapour pressure
  if (cav_on) then
     call chk_int(Cav_model, 1, 3, 'Cav_model')
     call chk_int(Vapor_scheme, 1, 3, 'Vapor_scheme')
     if (rho_v <= 0.0_dp .or. mu_v <= 0.0_dp) call nml_err('Lou_v and mu_v must be > 0')
     if (SS_Nuclei <= 0.0_dp .or. SS_Radius <= 0.0_dp) call nml_err('SS_Nuclei and SS_Radius must be > 0')
     if (ZGB_Rb <= 0.0_dp .or. ZGB_Rnuc <= 0.0_dp) call nml_err('ZGB_Rb and ZGB_Rnuc must be > 0')
     rv_ratio  = rho_v / rho1
     muv_ratio = mu_v / mu1
     pv_nd     = Prs_Vapor / (rho1*Uch*Uch)
     call chk_int(Pv_of_T, 0, 1, 'Pv_of_T')
     if (Latent_Heat < 0.0_dp) call nml_err('Latent_Heat must be >= 0')
     call chk_int(Cav_compressible, 0, 1, 'Cav_compressible')
     call chk_int(Cav_thermal, 0, 1, 'Cav_thermal')
     if (Snd_liquid <= 0.0_dp) call nml_err('Snd_liquid must be > 0')
     if (Cav_thermal == 1 .and. .not. energy_on) call nml_err('Cav_thermal = 1 needs SW_Energy = 1')
     if (Cav_thermal == 1) Pv_of_T = 1
     if (Cav_thermal_t < 0.0_dp) call nml_err('Cav_thermal_t must be >= 0 (0 = Char_L/Char_U)')
     if (Cav_thermal_t == 0.0_dp) Cav_thermal_t = 1.0e-3_dp
     if (vap_init < 0.0_dp .or. vap_init > 1.0_dp) call nml_err('vap_cel_initial must lie in [0, 1]')
  end if
  ! energy equation (update 16): property ratios, Peclet / Prandtl numbers
  if (energy_on) then
     if (T_init <= 0.0_dp) call nml_err('T_init must be > 0 (temperatures in K)')
     if (Cp_l <= 0.0_dp .or. k_l <= 0.0_dp) call nml_err('Cp_l and k_l must be > 0')
     if (Num_of_Phase == 2 .and. (Cp_a <= 0.0_dp .or. k_a <= 0.0_dp)) call nml_err('Cp_a and k_a must be > 0 (two-phase)')
     if (cav_on .and. (Cp_v <= 0.0_dp .or. k_v <= 0.0_dp)) call nml_err('Cp_v and k_v must be > 0 (cavitation)')
     if (Pr_t <= 0.0_dp) call nml_err('Pr_t must be > 0')
     call chk_int(Cnvt_scheme_Engy, 1, 4, 'Cnvt_scheme_Engy')
     call chk_int(Cnvt_mode_Engy,   1, 3, 'Cnvt_mode_Engy')
     call chk_int(Dffs_mode_Engy,   1, 3, 'Dffs_mode_Engy')
     if (Max_Itr_Engy < 1 .or. Min_Itr_Engy < 0 .or. Rsd_Lmt_Engy <= 0.0_dp .or. SOR_fct_Engy <= 0.0_dp &
         .or. cfl_switch_Engy <= 0.0_dp .or. vn_switch_Engy <= 0.0_dp) call nml_err('energy solver controls (_Engy) must be positive')
     Min_Itr_Engy = min(Min_Itr_Engy, Max_Itr_Engy)
     Pe_l = rho1*Cp_l*Uch*Lch/k_l
     Pr_l = mu1*Cp_l/k_l
     rcp_a = 1.0_dp;  kr_a = 1.0_dp;  Pr_a = Pr_l
     if (Num_of_Phase == 2) then
        rcp_a = rho2*Cp_a/(rho1*Cp_l);  kr_a = k_a/k_l;  Pr_a = mu2*Cp_a/k_a
     end if
     rcp_v = 1.0_dp;  kr_v = 1.0_dp;  Pr_v = Pr_l
     if (cav_on) then
        rcp_v = rho_v*Cp_v/(rho1*Cp_l);  kr_v = k_v/k_l;  Pr_v = mu_v*Cp_v/k_v
     end if
     buoy_on = (Beta_l /= 0.0_dp .or. (Num_of_Phase == 2 .and. Beta_a /= 0.0_dp)) .and. any(grav /= 0.0_dp)
  end if

  if (myrank == 0) then
     write(*,'(A)')        ' ---------------- input summary (probCFD.nml) ----------------'
     write(*,'(A,A)')      '  Case_name        : ', trim(casename)
     write(*,'(A,ES12.4)') '  Char_L  [m]      : ', Lch
     write(*,'(A,ES12.4)') '  Char_U  [m/s]    : ', Uch
     write(*,'(A,ES12.4)') '  Re (fluid 1)     : ', Re
     write(*,'(A,ES12.4)') '  Del_T   [s]      : ', dt_dim
     write(*,'(A,ES12.4)') '  dt*     [-]      : ', dt
     write(*,'(A,ES12.4,A,I0,A)') '  End_T   [s]      : ', End_T, '  (', nml_nsteps, ' steps)'
     write(*,'(A,I0)')     '  time steps       : ', nstep
     if (t_ramp > 0.0_dp) write(*,'(A,ES12.4)') '  BC ramp [s]      : ', t_ramp
     write(*,'(A,I0)')     '  Num_of_Phase     : ', Num_of_Phase
     write(*,'(A,2ES12.4)')'  fluid 1 rho, mu  : ', rho1, mu1
     write(*,'(A,ES9.2,A,ES9.2)') '  stop criteria    : NaN anywhere; residual / update norm > ', Stop_Rsd_Max, '; max CFL > ', Stop_CFL_Max
     select case (TURB_model)
     case (0);  write(*,'(A)') '  turbulence       : none'//merge(' (IDDES WMLES branch only)', '                          ', IDDES_model == 1)
     case (1);  write(*,'(A)') '  turbulence       : Spalart-Allmaras'//merge(' + IDDES', '        ', IDDES_model == 1)
     case (2);  write(*,'(A)') '  turbulence       : k-epsilon'//merge(' + IDDES', '        ', IDDES_model == 1)
     case (3);  write(*,'(A)') '  turbulence       : k-omega (Wilcox)'//merge(' + IDDES', '        ', IDDES_model == 1)
     case (4);  write(*,'(A)') '  turbulence       : k-omega SST (Menter)'//merge(' + IDDES', '        ', IDDES_model == 1)
     end select
     if (turb_on) write(*,'(A,I0,A,I0,A,I0,A,I0,A,ES9.2,A,F5.2)') '  turbulence solver: scheme ', Cnvt_scheme_Turb, &
        ', conv mode ', Cnvt_mode_Turb, ', diff mode ', Dffs_mode_Turb, ', max itr ', Max_Itr_Turb, ', rsd ', Rsd_Lmt_Turb, ', SOR ', SOR_fct_Turb
     if (Num_of_Phase == 2) then
        write(*,'(A,2ES12.4)')'  fluid 2 rho, mu  : ', rho2, mu2
        write(*,'(A,ES12.4)') '  FreeSurf_Z [m]   : ', z_free_surface
        if (Init_shape > 0) write(*,'(A,I0,A,3ES11.3,A,ES11.3)') '  initial alpha    : drop shape ', Init_shape, ', centre [m]', Drop_Center, ', radius / half side [m]', Drop_Radius
        if (st_on) then
           write(*,'(A,ES10.3,A,ES10.3)') '  surface tension  : CSF on, sigma [N/m] = ', Sigma_ST, ', 1/We = sigma/(rho1 U^2 L) = ', sigma_nd
        else if (Surf_tension == 1) then
           write(*,'(A)') '  surface tension  : requested but inactive (single-phase run)'
        else
           write(*,'(A)') '  surface tension  : off'
        end if
        select case (vof_scheme)
        case (1);  write(*,'(A,I0,A)') '  VoF scheme       : ', vof_scheme, ' (explicit CICSAM, sub-cycled)'
        case (2);  write(*,'(A,I0,A)') '  VoF scheme       : ', vof_scheme, ' (implicit CICSAM, iterative)'
        case default; write(*,'(A,I0,A,F5.2,A)') '  VoF scheme       : ', vof_scheme, ' (geometric PLIC, direction-split, Weymouth-Yue; face CFL limit ', min(vof_cfl_max, 0.5_dp), ')'
        end select
        write(*,'(A,ES12.4)') '  VoF Courant max  : ', vof_cfl_max
     end if
     write(*,'(A,3ES12.4)')'  gravity [m/s2]   : ', grav
     write(*,'(A,I0,A)')   '  Cnvt_scheme      : ', Cnvt_scheme, trim(scheme_name(Cnvt_scheme))
     write(*,'(A,I0,A)')   '  curvi form       : ', curvi_form,          &
          merge(' (unweighted U^q)  ', ' (J-weighted flux) ', curvi_form == 1)
     write(*,'(A,I0,A)')   '  Cnvt_mode        : ', Cnvt_mode, trim(mode_name(Cnvt_mode, 'CFL'))
     write(*,'(A,I0,A)')   '  Dffs_mode        : ', Dffs_mode, trim(mode_name(Dffs_mode, 'VN '))
     if (Cnvt_mode == 3) write(*,'(A,ES12.4)') '  cfl_switch       : ', cfl_switch
     if (Dffs_mode == 3) write(*,'(A,ES12.4)') '  vn_switch        : ', vn_switch
     if (quasi2d > 0) write(*,'(A,I0,A)') '  SW_Quasi_2D      : ', quasi2d, ' (strictly 2-D: w = 0, the fields averaged along k after every correction)'
     if (Prs_form == 1) then
        write(*,'(A,I0,A)') '  Prs_form         : ', Prs_form, ' (non-incremental: full pressure from the Poisson equation)'
     else
        write(*,'(A,I0,A)') '  Prs_form         : ', Prs_form, ' (incremental: grad p^n in the predictor, the increment is projected)'
     end if
     write(*,'(A,I0)')     '  Iner_Itn         : ', Iner_Itn
     if (Time_scheme == 2) then
        write(*,'(A,I0,A)') '  Time_scheme      : ', Time_scheme, ' (2nd-order SSP, Heun: two stages per step)'
     else
        write(*,'(A,I0,A)') '  Time_scheme      : ', Time_scheme, ' (1st order in time)'
     end if
     write(*,'(A,I0,A,I0,A,ES9.2,A,F5.2)') '  momentum         : Jacobi, Max/Min_Itr_Vel ', &
          mom_maxiter, '/', mom_miniter, ', Rsd_Lmt_Vel ', mom_tol, ', SOR ', mom_sor
     write(*,'(A,I0,A)')   '  Poisson solver   : ', poi_solver,          &
          merge(' (Jacobi: criterion = update norm Trl_Lmt_Prs)           ', &
                ' (FGMRES + 7-point PCG: criterion = residual Rsd_Lmt_Prs)', poi_solver == 1)
     if (poi_solver == 2) write(*,'(A,I0,A)') '  PCG precond.     : ', poi_precond, &
          merge(' (Jacobi)  ', ' (line SGS)', poi_precond == 1)
     write(*,'(A,I0,A,I0,A,ES9.2,A,ES9.2,A,F5.2)')                         &
          '  Poisson controls : Max/Min_Itr_Prs ', poi_maxiter, '/', poi_miniter, &
          ', Trl_Lmt_Prs ', poi_tol, ', Rsd_Lmt_Prs ', poi_rtol, ', SOR ', poi_sor
     if (Num_of_Phase == 2) write(*,'(A,I0,A,I0,A,ES9.2,A,F5.2)')          &
          '  VoF controls     : Max/Min_Itr_VoF ', vof_maxiter, '/', vof_miniter, &
          ', Rsd_Lmt_VoF ', vof_tol, ', SOR ', vof_sor
     write(*,'(A,L1)')     '  SW_Ref_Prs       : ', p_ref_switch
     write(*,'(A,I0)')     '  Output_Intvl     : ', out_interval
     write(*,'(A,I0)')     '  Data_Intvl       : ', data_interval
     if (restart_interval > 0) then
        write(*,'(A,I0,A)') '  Restart_Intvl    : ', restart_interval, ' (restart files <case>_NNNNNNNN.rsto)'
     else
        write(*,'(A)')      '  Restart_Intvl    : 0 (no restart files)'
     end if
     if (run_mode == 2) then
        write(*,'(A,I0,A)') '  Mode             : ', run_mode, ' (grid-movement check of the imposed motions: no flow solution)'
     else if (run_mode == 1) then
        write(*,'(A,I0,A)') '  Mode             : ', run_mode, ' (restart from the latest .rsto file)'
     else
        write(*,'(A,I0,A)') '  Mode             : ', run_mode, ' (standard: from time step 0)'
     end if
     if (grid_scale /= 1.0_dp) write(*,'(A,ES12.4)') '  Grid_scale       : ', grid_scale
     if (cav_on) then
        select case (Cav_model)
        case (1);  write(*,'(A,ES11.3,A,ES10.2,A)') '  Cavitation       : ON, Schnerr-Sauer (n0 ', SS_Nuclei, ' 1/m3, R0 ', SS_Radius, ' m)'
        case (2);  write(*,'(A,2F8.3,A,2ES10.2,A)') '  Cavitation       : ON, Zwart-Gerber-Belamri (F_vap F_cond ', ZGB_Fvap, ZGB_Fcond, &
                      ', r_nuc R_b ', ZGB_Rnuc, ZGB_Rb, ')'
        case default;  write(*,'(A,2F8.1,A)') '  Cavitation       : ON, Kunz (C_dest C_prod ', Kunz_Cdest, Kunz_Cprod, ', t_inf = Char_L/Char_U)'
        end select
        write(*,'(A,F9.1,A,ES10.3,A,ES10.3,A)') '    vapour: p_v ', Prs_Vapor, ' Pa, rho_v ', rho_v, ' kg/m3, mu_v ', mu_v, ' kg/(m s)'
        select case (Vapor_scheme)
        case (1);  write(*,'(A)') '    vapour transport: explicit CICSAM (sub-cycled)'
        case (2);  write(*,'(A)') '    vapour transport: implicit CICSAM'
        case default;  write(*,'(A)') '    vapour transport: 2nd-order ENO (minmod) face values, time integration by Cnvt_mode'
        end select
        write(*,'(A,ES11.4,A)') '    mass transfer from t = ', Cav_Start_T, ' s (50-step ramp); Reboud correction '// &
           merge('on ', 'off', Reboud == 1)
     end if
     if (cav_on .and. energy_on) write(*,'(A,ES10.3,A)') '    latent heat ', Latent_Heat, ' J/kg; vapour pressure '// &
        merge('p_v(T) saturation law (Pv_of_T = 1)', 'constant Prs_Vapor                 ', Pv_of_T == 1)
     if (cav_on .and. Cav_compressible == 1) write(*,'(A,F8.1,A,F8.1,A,ES10.3)') &
        '    compressible mixture (Wood law): c_liquid ', Snd_liquid, ' m/s, c_vapour (isothermal gas at p_v) ', &
        sqrt(Prs_Vapor/rho_v), ' m/s; Ma^2 of the liquid = ', (Uch/Snd_liquid)**2
     if (cav_on .and. Cav_thermal == 1) write(*,'(A,ES10.3,A)') &
        '    thermal cavitation: Clausius-Clapeyron p_v(T), rho_v(T) through (Prs_Vapor, Lou_v, T_init), B-factor '// &
        'delay dT* = 1/2 (rho_v L/(rho_l c_p)) R_B/sqrt(alpha t) with t = ', Cav_thermal_t, ' s, implicit latent heat'
     if (cav_on .and. Cav_thermal == 1) write(*,'(A,ES10.3,A,ES10.3,A)') '    rho_v L/(rho_l c_p) = ', &
        rho_v*Latent_Heat/(rho1*Cp_l), ' K (B-factor scale); Brennen Sigma(T_init) = ', &
        (rho_v*Latent_Heat)**2/(rho1**2*Cp_l*T_init*sqrt(k_l/(rho1*Cp_l))), ' m/s^1.5 (thermal effects matter when Sigma sqrt(t) ~ U)'
     if (cav_on .and. vap_init > 0.0_dp) write(*,'(A,F8.4)') '    initial vapour fraction of the liquid cells: ', vap_init
     if (energy_on) then
        write(*,'(A,F9.3,A)') '  Energy equation  : ON, T_init ', T_init, ' K; temperature in K at the cells, Temperature in the .p3d files'
        write(*,'(A,F8.1,A,F8.4,A,ES10.3,A,F7.3)') '    liquid: c_p ', Cp_l, ' J/(kg K), k ', k_l, ' W/(m K), Pe = rho c_p U L/k = ', Pe_l, ', Pr = ', Pr_l
        if (Num_of_Phase == 2) write(*,'(A,F8.1,A,F8.4,A,F7.3)') '    air   : c_p ', Cp_a, ' J/(kg K), k ', k_a, ' W/(m K), Pr = ', Pr_a
        if (cav_on) write(*,'(A,F8.1,A,F8.4,A,F7.3)') '    vapour: c_p ', Cp_v, ' J/(kg K), k ', k_v, ' W/(m K), Pr = ', Pr_v
        if (buoy_on) then
           write(*,'(A,ES10.3,A,ES10.3,A)') '    buoyancy: Boussinesq -beta (T - T_init) g, Beta_l ', Beta_l, ' 1/K, Beta_a ', Beta_a, ' 1/K'
        else
           write(*,'(A)') '    buoyancy: off (Beta_l = Beta_a = 0 or no gravity): T is a passive scalar'
        end if
        write(*,'(A,F6.3,A,I0,A,I0,A,I0,A,I0,A,I0,A,ES9.2,A,F5.2)') '    Pr_t ', Pr_t, '; solver: scheme ', Cnvt_scheme_Engy, &
           ', conv mode ', Cnvt_mode_Engy, ', diff mode ', Dffs_mode_Engy, ', max/min itr ', Max_Itr_Engy, '/', Min_Itr_Engy, &
           ', rsd ', Rsd_Lmt_Engy, ', SOR ', SOR_fct_Engy
        if (Cnvt_mode_Engy /= 2) write(*,'(A,F6.3,A,F6.3)') '    explicit convection below CFL ', cfl_switch_Engy, &
           ', explicit diffusion below VN ', vn_switch_Engy
     end if
     write(*,'(A,ES12.4)') '  b2b tolerance    : ', b2b_tol
     write(*,'(A,I0,A,ES10.3,A,I0,A,ES9.2)') '  overset          : SW_Ovst ', ovst_switch,   &
          ' (-1 auto), Surfovl_tol ', surfovl_tol, ' m, fringe layers ', ovst_nf,          &
          ', vol tol ', ovst_vtol
     if (ovst_reasm_tol > 0.0_dp) write(*,'(A,ES10.3,A)') '                     re-assembly only when a part moved more than ', &
          ovst_reasm_tol, ' m relative to the background (Ovst_Reassemble_Tol)'
     if (Moving_BGR == 1) write(*,'(A,I0,A,6I2,A,ES10.3,A)') '  moving background: follows Body ', BGR_Following_Body, &
          ' in the DOFs ', Moving_BGR_DoF, ', filter ', Moving_BGR_Filter_T, ' s (Moving_BGR, update 20)'
  end if

contains
  !--------------------------------------------------------------------
  subroutine assign(nm)
    character(len=*), intent(in) :: nm
    select case (nm)
    ! $Case_title
    case ('case_name', 'casename');           casename = unquote(vals(1))
    case ('mode', 'run_mode');                run_mode = nint(r(1))
    ! moving background (update 20)
    case ('moving_bgr');                      Moving_BGR = nint(r(1))
    case ('moving_bgr_dof')
       if (nv < 6) call nml_err_line('Moving_BGR_DoF needs six 0/1 flags (surge sway heave roll pitch yaw)')
       Moving_BGR_DoF = nint([r(1), r(2), r(3), r(4), r(5), r(6)])
    case ('bgr_following_body');              BGR_Following_Body = nint(r(1))
    case ('moving_bgr_filter_t');             Moving_BGR_Filter_T = r(1)
    ! $Control
    case ('char_l', 'lch');                   Lch = r(1)
    case ('char_u', 'uch');                   Uch = r(1)
    case ('del_t', 'dt_dim');                 dt_dim = r(1)
    case ('end_t');                           End_T = r(1)
    case ('freesurf_z', 'z_free_surface');    z_free_surface = r(1)
    case ('grvty', 'grav')
       if (nv < 3) call nml_err_line('Grvty needs three components')
       grav = [r(1), r(2), r(3)]
    ! $Fluid_Properties
    case ('num_of_phase');                    Num_of_Phase = nint(r(1))
    case ('lou_w', 'rho1');                   rho1 = r(1)
    case ('mu_w', 'mu1');                     mu1  = r(1)
    case ('lou_a', 'rho2');                   rho2 = r(1)
    case ('mu_a', 'mu2');                     mu2  = r(1)
    case ('rho');                             rho  = r(1)
    case ('mu');                              mu   = r(1)
    ! $Reference_Pressure
    case ('sw_ref_prs');                      p_ref_switch = (nint(r(1)) /= 0)
    case ('p_ref_switch');                    p_ref_switch = lg(vals(1))
    case ('loc_x', 'p_ref_x');                p_ref_x = r(1)
    case ('loc_y', 'p_ref_y');                p_ref_y = r(1)
    case ('loc_z', 'p_ref_z');                p_ref_z = r(1)
    case ('refprs', 'p_ref_val');             p_ref_val = r(1)
    ! $Cell_Initial_Condition
    case ('u_cel_initial', 'u_init');         u_init = r(1)
    case ('v_cel_initial', 'v_init');         v_init = r(1)
    case ('w_cel_initial', 'w_init');         w_init = r(1)
    case ('a_cel_initial');                   alpha_init = r(1)
    case ('p_cel_initial', 'p_init');         p_init = r(1)
    ! $Solver_Iteration
    case ('iner_itn');                        Iner_Itn = nint(r(1))
    case ('max_itr_vel', 'mom_maxiter');      mom_maxiter = nint(r(1))
    case ('min_itr_vel');                     mom_miniter = nint(r(1))
    case ('rsd_lmt_vel', 'mom_tol');          mom_tol = r(1)
    case ('sor_fct_vel');                     mom_sor = r(1)
    case ('max_itr_prs', 'poi_maxiter');      poi_maxiter = nint(r(1))
    case ('min_itr_prs');                     poi_miniter = nint(r(1))
    case ('trl_lmt_prs', 'poi_tol');          poi_tol = r(1)
    case ('rsd_lmt_prs');                     poi_rtol = r(1)
    case ('sor_fct_prs');                     poi_sor = r(1)
    case ('max_itr_vof');                     vof_maxiter = nint(r(1))
    case ('min_itr_vof');                     vof_miniter = nint(r(1))
    case ('rsd_lmt_vof');                     vof_tol = r(1)
    case ('sor_fct_vof');                     vof_sor = r(1)
    ! $Export_Interval
    case ('output_intvl', 'out_interval');    out_interval = nint(r(1))
    case ('log_intvl');                       log_interval = nint(r(1))
    case ('data_intvl', 'data_interval');     data_interval = nint(r(1))
    case ('restart_intvl', 'restart_interval'); restart_interval = nint(r(1))
    case ('grid_scale');                      grid_scale = r(1)
    ! $Numerical_Scheme (probCFD switches)
    case ('cnvt_scheme', 'conv_scheme');      Cnvt_scheme = nint(r(1))
    case ('curvi_form');                      curvi_form = nint(r(1))
    case ('cnvt_mode');                       Cnvt_mode = nint(r(1))
    case ('prs_form');                        Prs_form = nint(r(1))
    case ('sw_quasi_2d');                     quasi2d = nint(r(1))
    case ('dffs_mode');                       Dffs_mode = nint(r(1))
    case ('cfl_switch');                      cfl_switch = r(1)
    case ('vn_switch');                       vn_switch = r(1)
    case ('poi_solver');                      poi_solver = nint(r(1))
    case ('poi_precond');                     poi_precond = nint(r(1))
    case ('vof_scheme');                      vof_scheme = nint(r(1))
    case ('surf_tension');                    Surf_tension = nint(r(1))
    case ('init_shape');                      Init_shape = nint(r(1))
    case ('drop_radius');                     Drop_Radius = r(1)
    case ('drop_center')
       if (nv < 3) call nml_err_line('Drop_Center needs three components')
       Drop_Center = [r(1), r(2), r(3)]
    case ('sigma_st');                        Sigma_ST = r(1)
    case ('stop_rsd_max');                    Stop_Rsd_Max = r(1)
    case ('stop_cfl_max');                    Stop_CFL_Max = r(1)
    case ('turb_model');                      TURB_model = nint(r(1))
    case ('iddes_model');                     IDDES_model = nint(r(1))
    case ('max_itr_turb');                    Max_Itr_Turb = nint(r(1))
    case ('min_itr_turb');                    Min_Itr_Turb = nint(r(1))
    case ('rsd_lmt_turb');                    Rsd_Lmt_Turb = r(1)
    case ('sor_fct_turb');                    SOR_fct_Turb = r(1)
    case ('cnvt_scheme_turb');                Cnvt_scheme_Turb = nint(r(1))
    case ('cnvt_mode_turb');                  Cnvt_mode_Turb = nint(r(1))
    case ('dffs_mode_turb');                  Dffs_mode_Turb = nint(r(1))
    case ('cfl_switch_turb');                 cfl_switch_Turb = r(1)
    case ('vn_switch_turb');                  vn_switch_Turb = r(1)
    case ('vof_cfl_max');                     vof_cfl_max = r(1)
    case ('t_ramp');                          t_ramp = r(1)
    case ('time_scheme');                     Time_scheme = nint(r(1))
    ! $Cavitation (update 15)
    case ('sw_cav');                          cav_on = (nint(r(1)) /= 0)
    case ('cav_model');                       Cav_model = nint(r(1))
    case ('vapor_scheme', 'vapour_scheme');   Vapor_scheme = nint(r(1))
    case ('cav_start_t');                     Cav_Start_T = r(1)
    case ('prs_vapor', 'prs_vapour', 'p_vapor');  Prs_Vapor = r(1)
    case ('lou_v', 'rho_v');                  rho_v = r(1)
    case ('mu_v');                            mu_v = r(1)
    case ('ss_nuclei');                       SS_Nuclei = r(1)
    case ('ss_radius');                       SS_Radius = r(1)
    case ('zgb_fvap');                        ZGB_Fvap = r(1)
    case ('zgb_fcond');                       ZGB_Fcond = r(1)
    case ('zgb_rnuc');                        ZGB_Rnuc = r(1)
    case ('zgb_rb');                          ZGB_Rb = r(1)
    case ('kunz_cdest');                      Kunz_Cdest = r(1)
    case ('kunz_cprod');                      Kunz_Cprod = r(1)
    case ('reboud');                          Reboud = nint(r(1))
    case ('latent_heat');                     Latent_Heat = r(1)
    case ('pv_of_t');                         Pv_of_T = nint(r(1))
    case ('cav_compressible');                Cav_compressible = nint(r(1))
    case ('snd_liquid');                      Snd_liquid = r(1)
    case ('cav_thermal');                     Cav_thermal = nint(r(1))
    case ('cav_thermal_t');                   Cav_thermal_t = r(1)
    case ('vap_cel_initial', 'vap_init');     vap_init = r(1)
    ! $Energy (update 16)
    case ('sw_energy');                       energy_on = (nint(r(1)) /= 0)
    case ('t_init', 't_cel_initial');         T_init = r(1)
    case ('cp_l', 'cp_w');                    Cp_l = r(1)
    case ('k_l', 'k_w');                      k_l = r(1)
    case ('cp_a');                            Cp_a = r(1)
    case ('k_a');                             k_a = r(1)
    case ('cp_v');                            Cp_v = r(1)
    case ('k_v');                             k_v = r(1)
    case ('beta_l', 'beta_w');                Beta_l = r(1)
    case ('beta_a');                          Beta_a = r(1)
    case ('pr_t');                            Pr_t = r(1)
    case ('max_itr_engy');                    Max_Itr_Engy = nint(r(1))
    case ('min_itr_engy');                    Min_Itr_Engy = nint(r(1))
    case ('rsd_lmt_engy');                    Rsd_Lmt_Engy = r(1)
    case ('sor_fct_engy');                    SOR_fct_Engy = r(1)
    case ('cnvt_scheme_engy');                Cnvt_scheme_Engy = nint(r(1))
    case ('cnvt_mode_engy');                  Cnvt_mode_Engy = nint(r(1))
    case ('dffs_mode_engy');                  Dffs_mode_Engy = nint(r(1))
    case ('cfl_switch_engy');                 cfl_switch_Engy = r(1)
    case ('vn_switch_engy');                  vn_switch_Engy = r(1)
    case ('b2b_tol');                         b2b_tol = r(1)
    ! $Overset
    case ('sw_ovst', 'ovst_switch');          ovst_switch = nint(r(1))
    case ('surfovl_tol');                     surfovl_tol = r(1)
    case ('ovst_fringe', 'ovst_nf');          ovst_nf = nint(r(1))
    case ('ovst_wall_layers');                ovst_wall_layers = nint(r(1))
    case ('ovst_active_overlap');             ovst_active_overlap = nint(r(1))
    case ('ovst_vol_tol', 'ovst_vtol');       ovst_vtol = r(1)
    case ('ovst_reassemble_tol');             ovst_reasm_tol = r(1)
    case ('ovst_only');                       ovst_only = nint(r(1))
    case ('poi_verbose');                     poi_verbose = (nint(r(1)) /= 0)
    case ('ovst_inner_max');                  ovst_inner_max = nint(r(1))
    case ('ovst_inner_tol');                  ovst_inner_tol = r(1)
    ! $Export_Interval: wall surface files (update 22)
    case ('sw_prs_frc_surf');                 sw_prs_frc_surf = nint(r(1))
    case ('frame_prs_frc_surf');              frame_prs_frc_surf = nint(r(1))
    ! $Statistics (update 23, statistics.F90)
    case ('sw_statistics');                   sw_statistics = nint(r(1))
    case ('statistics_type1');                stat_type(1) = nint(r(1))
    case ('statistics_type2');                stat_type(2) = nint(r(1))
    case ('statistics_type3');                stat_type(3) = nint(r(1))
    case ('statistics_type4');                stat_type(4) = nint(r(1))
    case ('statistics_first_group');          stat_group(1) = nint(r(1))
    case ('statistics_second_group');         stat_group(2) = nint(r(1))
    case ('statistics_start_t');              stat_start_t = r(1)
    case ('statistics_end_t');                stat_end_t = r(1)
    case ('statistics_sample_intvl');         stat_sample_intvl = max(nint(r(1)), 1)
    case ('statistics_output_intvl');         stat_output_intvl = nint(r(1))
    case ('statistics_type1_frame');          stat_type1_frame = nint(r(1))
    case ('statistics_box', 'statistics_cyl', 'statistics_plane')
       if (nlatt_def >= MAXLATT) call nml_err_line('too many statistics lattices (max 64)')
       nlatt_def = nlatt_def + 1
       latt_def(nlatt_def)%id = nint(r(1));  latt_def(nlatt_def)%open = .true.
       select case (nm)
       case ('statistics_box');  latt_def(nlatt_def)%kind = 2
       case ('statistics_cyl');  latt_def(nlatt_def)%kind = 3
       case default;             latt_def(nlatt_def)%kind = 4
       end select
    case ('statistics_box_end', 'statistics_cyl_end', 'statistics_plane_end')
       if (nlatt_def == 0) call nml_err_line('"'//nm//'" without an open Statistics_box / _cyl / _plane group')
       latt_def(nlatt_def)%open = .false.
    case ('statistics_box_origin', 'statistics_plane_origin', 'statistics_cyl_facecenter')
       call latt_chk(nm);  latt_def(nlatt_def)%org = [r(1), r(2), r(3)]
    case ('statistics_box_axis1', 'statistics_plane_axis1', 'statistics_cyl_axis')
       call latt_chk(nm);  latt_def(nlatt_def)%ax1 = [r(1), r(2), r(3)]
    case ('statistics_box_axis2', 'statistics_plane_axis2', 'statistics_cyl_ref')
       call latt_chk(nm);  latt_def(nlatt_def)%ax2 = [r(1), r(2), r(3)]
    case ('statistics_box_length');   call latt_chk(nm);  latt_def(nlatt_def)%len = [r(1), r(2), r(3)]
    case ('statistics_plane_length'); call latt_chk(nm);  latt_def(nlatt_def)%len = [r(1), r(2), 0.0_dp]
    case ('statistics_cyl_length');   call latt_chk(nm);  latt_def(nlatt_def)%len = [r(1), 0.0_dp, 0.0_dp]
    case ('statistics_box_npoints', 'statistics_cyl_npoints')
       call latt_chk(nm);  latt_def(nlatt_def)%np = [nint(r(1)), nint(r(2)), nint(r(3))]
    case ('statistics_plane_npoints'); call latt_chk(nm);  latt_def(nlatt_def)%np = [nint(r(1)), nint(r(2)), 1]
    case ('statistics_cyl_radius');   call latt_chk(nm);  latt_def(nlatt_def)%rad = [r(1), r(2)]
    case ('statistics_cyl_theta');    call latt_chk(nm);  latt_def(nlatt_def)%theta = [r(1), r(2)]
    case ('statistics_cyl_polar_vel'); call latt_chk(nm);  latt_def(nlatt_def)%polar = nint(r(1))
    case default
       ! the mother of a statistics lattice (inside its group; outside, the keys belong to the force models)
       if (trim(group) == 'statistics' .and. nlatt_def > 0) then
          if (latt_def(nlatt_def)%open) then
             select case (nm)
             case ('mother_body');            latt_def(nlatt_def)%mb = nint(r(1));  return
             case ('mother_child_body');      latt_def(nlatt_def)%mc = nint(r(1));  return
             case ('mother_sub_child_body');  latt_def(nlatt_def)%ms = nint(r(1));  return
             end select
          end if
       end if
       ! $Probe_background (update 22): Probe_point<n> = x y z of a point following the background
       if (trim(group) == 'probe_background' .and. len(nm) > 11) then
          if (nm(1:11) == 'probe_point') then
             read(nm(12:), *, iostat=io) i
             if (io /= 0 .or. i < 1) call nml_err_line('"'//nm//'": the point number must follow Probe_point (Probe_point1, ...)')
             if (nv < 3) call nml_err_line('"'//nm//'" needs x y z [m]')
             if (nprobe_bg >= MAXPRBG) call nml_err_line('too many background probe points (max 256)')
             nprobe_bg = nprobe_bg + 1
             probe_bg_n(nprobe_bg) = i;  probe_bg_x(:, nprobe_bg) = [r(1), r(2), r(3)]
             return
          end if
       end if
       ! $SixDOF (per-body rigid-body data, sixdof.F90)
       ok6 = .false.
       ! Force_model groups inside $SixDOF (update 19): every line between Force_model = n and Force_model_End = -1
       if (nm == 'force_model' .or. force_in_group()) call force_nml(nm, vals, nv, ok6)
       if (.not. ok6 .and. (trim(group) == 'sixdof' .or. nm == 'motion_start_t')) call sixdof_nml(nm, vals, nv, ok6)
       ! $Wave (update 18, waves.F90)
       if (.not. ok6) call wave_nml(nm, vals, nv, ok6)
       if (.not. ok6) call nml_err_line('unknown variable "'//nm//'" (group $'//trim(group)//')')
    end select
  end subroutine assign

  subroutine latt_chk(nm)
    character(len=*), intent(in) :: nm
    if (nlatt_def == 0) call nml_err_line('"'//nm//'" must follow a Statistics_box / Statistics_cyl / Statistics_plane line')
    if (.not. latt_def(nlatt_def)%open) call nml_err_line('"'//nm//'" after the _End mark of the last statistics lattice')
  end subroutine latt_chk

  ! numeric value number n of the current line
  real(dp) function r(n)
    integer, intent(in) :: n
    integer :: io
    r = 0.0_dp
    if (n > nv) call nml_err_line('too few values for '//trim(name))
    read(vals(n), *, iostat=io) r
    if (io /= 0) call nml_err_line('cannot read a number from "'//trim(vals(n))//'" for '//trim(name))
  end function r

  logical function lg(s)
    character(len=*), intent(in) :: s
    character(len=16) :: t
    t = lower(adjustl(s))
    lg = (t(1:1) == 't' .or. t(1:2) == '.t' .or. t(1:1) == '1')
  end function lg

  function unquote(s) result(t)
    character(len=*), intent(in) :: s
    character(len=len(s)) :: t
    integer :: n
    t = adjustl(s)
    n = len_trim(t)
    if (n >= 2) then
       if ((t(1:1) == '"' .and. t(n:n) == '"') .or. (t(1:1) == "'" .and. t(n:n) == "'")) &
          t = t(2:n-1)
    end if
  end function unquote

  function lower(s) result(t)
    character(len=*), intent(in) :: s
    character(len=len(s)) :: t
    integer :: i, c
    t = s
    do i = 1, len(s)
       c = iachar(s(i:i))
       if (c >= iachar('A') .and. c <= iachar('Z')) t(i:i) = achar(c + 32)
    end do
  end function lower

  ! split "v1, v2 v3" into tokens (quoted strings kept whole)
  subroutine split_values(s, out, n)
    character(len=*), intent(in)  :: s
    character(len=*), intent(out) :: out(:)
    integer,          intent(out) :: n
    integer :: i, L
    character(len=1) :: ch, qc
    logical :: inq
    n = 0;  out = ' ';  inq = .false.;  qc = ' '
    L = len_trim(s)
    i = 1
    do while (i <= L)
       ch = s(i:i)
       if (.not. inq .and. (ch == ' ' .or. ch == ',' .or. ch == char(9))) then
          i = i + 1;  cycle
       end if
       ! start of a token
       if (n >= size(out)) exit                 ! more values than slots: the rest is ignored
       n = n + 1
       if (ch == '"' .or. ch == "'") then
          qc = ch;  inq = .true.
          out(n) = ch
          i = i + 1
          do while (i <= L)
             out(n) = trim(out(n))//s(i:i)
             if (s(i:i) == qc) exit
             i = i + 1
          end do
          inq = .false.
          i = i + 1
       else
          do while (i <= L)
             if (s(i:i) == ' ' .or. s(i:i) == ',' .or. s(i:i) == char(9)) exit
             out(n) = trim(out(n))//s(i:i)
             i = i + 1
          end do
       end if
    end do
  end subroutine split_values

  subroutine chk_int(val, lo, hi, nm)
    integer, intent(in) :: val, lo, hi
    character(len=*), intent(in) :: nm
    if (val < lo .or. val > hi) then
       if (myrank == 0) write(*,'(A,A,A,I0,A,I0,A,I0,A)') ' ERROR: ',    &
          nm, ' = ', val, ' (must be ', lo, '..', hi, ')'
       call MPI_Abort(MPI_COMM_WORLD, 11, ierr)
    end if
  end subroutine chk_int
  subroutine nml_err(msg)
    character(len=*), intent(in) :: msg
    if (myrank == 0) write(*,'(A)') ' ERROR (probCFD.nml): '//msg
    call stop_report_input('probCFD.nml: '//msg)
  end subroutine nml_err
  subroutine nml_err_line(msg)
    character(len=*), intent(in) :: msg
    character(len=16) :: sl
    write(sl,'(I0)') lnum
    if (myrank == 0) write(*,'(A,I0,A)') ' ERROR (probCFD.nml, line ', lnum, '): '//msg
    call stop_report_input('probCFD.nml, line '//trim(sl)//': '//msg)
  end subroutine nml_err_line
  function scheme_name(m) result(s)
    integer, intent(in) :: m
    character(len=40) :: s
    select case (m)
    case (1);  s = ' (1st-order upwind)'
    case (2);  s = ' (2nd-order upwind)'
    case (3);  s = ' (QUICK)'
    case default;  s = ' (2nd-order ENO)'
    end select
  end function scheme_name
  function mode_name(m, crit) result(s)
    integer, intent(in) :: m
    character(len=*), intent(in) :: crit
    character(len=40) :: s
    select case (m)
    case (1);  s = ' (explicit)'
    case (2);  s = ' (implicit)'
    case default;  s = ' (hybrid per direction by local '//crit//')'
    end select
  end function mode_name
end subroutine read_namelist

!-----------------------------------------------------------------------
! Read the per-block MPI splitting from 'Input_polyMesh_probCFD.csv'.
! Must be called AFTER read_grid (the number of grid blocks N is taken
! from the PLOT3D grid file).  The csv contains one column per block:
!
!   probCFD polyMesh Input, block-1, block-2, ..., block-N
!   i_split , i1, i2, ..., iN      ! sub-domains in the block i-direction
!   j_split , j1, j2, ..., jN      ! sub-domains in the block j-direction
!   k_split , k1, k2, ..., kN      ! sub-domains in the block k-direction
!   n_of_BLK, c1, c2, ..., cN      ! CPUs per block = i*j*k (checked)
!   Object_ID, o1, o2, ..., oN     ! (row 6, optional) 6-DoF body of the
!                                  ! block: 0 background, 1 first body and
!                                  ! its refinement grids, 2 second body...
!   Child_Object_ID, c1, ..., cN   ! (row 7, optional) child body of that
!                                  ! body the block belongs to (0: the body)
!   Sub_Child_Object_ID, s1,.., sN ! (row 8, optional) sub-child body of
!                                  ! that child (0: the child itself)
!
! The total number of CPUs for the simulation is the SUM of row 5.
!-----------------------------------------------------------------------
subroutine read_split_csv
  use global_vars
  use mpi_vars
  use stop_criteria, only: stop_report_input
  implicit none
  integer :: iu, ios, ierr, b, a, ntot
  integer, allocatable :: vals(:), obid(:), chid(:), suid(:)
  character(len=8192)  :: line
  character(len=16)    :: labels(4)

  labels = [character(len=16) :: 'i_split', 'j_split', 'k_split', 'n_of_BLK']

  open(newunit=iu, file='Input_polyMesh_probCFD.csv', status='old',      &
       action='read', iostat=ios)
  if (ios /= 0) then
     if (myrank == 0) write(*,'(A)')                                     &
        ' ERROR: cannot open Input_polyMesh_probCFD.csv (it must have'// &
        ' one splitting column per grid block).'
     call MPI_Abort(MPI_COMM_WORLD, 32, ierr)
  end if

  allocate(vals(nblk))
  read(iu,'(A)', iostat=ios) line                     ! header row
  do a = 1, 4                                         ! i/j/k splits, cpus
     read(iu,'(A)', iostat=ios) line
     if (ios /= 0) call csv_err('missing row "'//trim(labels(a))//'"')
     call parse_row(line, trim(labels(a)), vals)
     if (a <= 3) then
        do b = 1, nblk
           blk(b)%sp(a) = vals(b)
        end do
     else
        do b = 1, nblk
           blk(b)%ncpu = vals(b)
           if (blk(b)%ncpu /= product(blk(b)%sp)) then
              if (myrank == 0) write(*,'(A,I0,A,I0,A,I0)')               &
                 ' ERROR: csv row 5 (CPUs) of block ', b, ' is ',        &
                 blk(b)%ncpu, ' but i*j*k of rows 2-4 gives ',           &
                 product(blk(b)%sp)
              call MPI_Abort(MPI_COMM_WORLD, 33, ierr)
           end if
        end do
     end if
  end do
  ! rows 6-8 (optional): object id (body), child id, sub-child id of
  ! every block (6-DoF, update 14): Object_ID = body the block belongs
  ! to (0 background); Child_Object_ID = child body of that body the
  ! block belongs to (0: the body itself); Sub_Child_Object_ID = sub-
  ! child of that child (0: the child itself)
  allocate(obid(nblk), chid(nblk), suid(nblk));  obid = 0;  chid = 0;  suid = 0
  do a = 1, 3
     read(iu,'(A)', iostat=ios) line
     if (ios /= 0) exit
     if (len_trim(line) == 0 .or. index(line, ',') == 0) exit
     call parse_row(line, ' ', vals)                   ! any label accepted
     select case (a)
     case (1);  obid = vals
     case (2);  chid = vals
     case default;  suid = vals
     end select
  end do
  close(iu)
  if (any(obid < 0)) call csv_err('row 6 (Object_ID) must be >= 0')
  if (any(chid < 0)) call csv_err('row 7 (Child_Object_ID) must be >= 0')
  if (any(suid < 0)) call csv_err('row 8 (Sub_Child_Object_ID) must be >= 0')
  do b = 1, nblk
     if (chid(b) > 0 .and. obid(b) == 0) call csv_err('a block with Child_Object_ID > 0 needs Object_ID > 0 (its mother body)')
     if (suid(b) > 0 .and. chid(b) == 0) call csv_err('a block with Sub_Child_Object_ID > 0 needs Child_Object_ID > 0 (its mother child body)')
  end do
  call build_node_table

  ntot = sum(blk(1:nblk)%ncpu)
  if (myrank == 0) then
     write(*,'(A,I0,A)') '  csv splitting   : read for ', nblk, ' blocks'
     write(*,'(A,I0)')   '  total CPUs (csv): ', ntot
     if (n6node > 0) then
        write(*,'(A,I0,A,I0,A,I0,A,I0,A)') '  csv object ids  : ', n6node, ' 6-DoF node(s) = ', count(node_lev == 1), ' body(ies), ', &
           count(node_lev == 2), ' child body(ies), ', count(node_lev == 3), ' sub-child body(ies)'
        write(*,'(A,100I3)') '    body id per block     : ', obid(1:min(nblk,100))
        if (maxval(chid) > 0) write(*,'(A,100I3)') '    child id per block    : ', chid(1:min(nblk,100))
        if (maxval(suid) > 0) write(*,'(A,100I3)') '    sub-child id per block: ', suid(1:min(nblk,100))
        if (maxval(chid) > 0) write(*,'(A,100I3)') '    node id per block     : ', blk_obj(1:min(nblk,100))
     end if
  end if
  deallocate(obid, chid, suid)

contains

  !-------------------------------------------------------------------
  ! node table from the (body, child, sub) ids of the blocks: the
  ! bodies first (node = body id, every id 1..nb must exist), then the
  ! children sorted by (body, child), then the sub-children sorted by
  ! (body, child, sub); the child / sub-child ids of one parent must be
  ! 1, 2, ... without gaps.  blk_obj(b) becomes the node id.
  !-------------------------------------------------------------------
  subroutine build_node_table
    integer :: nb, k, c, s, n, lev, ncmax, nsmax, nc, ns, par
    nb = maxval(obid);  ncmax = maxval(chid);  nsmax = maxval(suid)
    if (allocated(blk_obj)) deallocate(blk_obj)
    allocate(blk_obj(nblk));  blk_obj = 0
    if (allocated(node_bid)) deallocate(node_bid, node_cid, node_sid, node_lev, node_par)
    n6node = 0
    if (nb == 0) then
       allocate(node_bid(0), node_cid(0), node_sid(0), node_lev(0), node_par(0))
       return
    end if
    n = nblk + nb                                      ! upper bound of the node count
    allocate(node_bid(n), node_cid(n), node_sid(n), node_lev(n), node_par(n))
    node_bid = 0;  node_cid = 0;  node_sid = 0;  node_lev = 0;  node_par = 0
    do lev = 1, 3
       do k = 1, nb
          if (lev == 1) then
             if (count(obid == k) == 0) call csv_err('row 6: body id '//itoa(k)//' has no block (ids must be 1, 2, ... without gaps)')
             n6node = n6node + 1
             node_bid(n6node) = k;  node_lev(n6node) = 1;  node_par(n6node) = 0
             cycle
          end if
          nc = 0
          do c = 1, ncmax
             if (count(obid == k .and. chid == c) > 0) nc = c
          end do
          do c = 1, nc
             if (count(obid == k .and. chid == c) == 0) call csv_err('row 7: body '//itoa(k)//' has no block of child '// &
                itoa(c)//' (child ids of one body must be 1, 2, ... without gaps)')
             if (lev == 2) then
                n6node = n6node + 1
                node_bid(n6node) = k;  node_cid(n6node) = c;  node_lev(n6node) = 2;  node_par(n6node) = k
                cycle
             end if
             ns = 0
             do s = 1, nsmax
                if (count(obid == k .and. chid == c .and. suid == s) > 0) ns = s
             end do
             par = 0
             do n = 1, n6node
                if (node_lev(n) == 2 .and. node_bid(n) == k .and. node_cid(n) == c) par = n
             end do
             do s = 1, ns
                if (count(obid == k .and. chid == c .and. suid == s) == 0) call csv_err('row 8: child '//itoa(c)//' of body '// &
                   itoa(k)//' has no block of sub-child '//itoa(s)//' (sub-child ids must be 1, 2, ... without gaps)')
                n6node = n6node + 1
                node_bid(n6node) = k;  node_cid(n6node) = c;  node_sid(n6node) = s;  node_lev(n6node) = 3;  node_par(n6node) = par
             end do
          end do
       end do
    end do
    do b = 1, nblk
       if (obid(b) == 0) cycle
       do n = 1, n6node
          if (node_bid(n) == obid(b) .and. node_cid(n) == chid(b) .and. node_sid(n) == suid(b)) blk_obj(b) = n
       end do
    end do
  end subroutine build_node_table

  function itoa(i) result(s)
    integer, intent(in) :: i
    character(len=:), allocatable :: s
    character(len=12) :: buf
    write(buf,'(I0)') i
    s = trim(buf)
  end function itoa

  subroutine parse_row(row, label, out)
    character(len=*), intent(in)  :: row
    character(len=*), intent(in)  :: label
    integer,          intent(out) :: out(:)
    character(len=len(row)) :: work
    integer :: p, q, n, io_
    work = row
    ! strip a possible carriage return (files written on Windows)
    p = index(work, char(13));  if (p > 0) work(p:) = ' '
    p = index(work, ',')
    if (p <= 1) call csv_err('malformed row (expected "'//label//'")')
    if (len_trim(label) > 0) then
       if (index(work(1:p-1), label) == 0) call csv_err(                 &
          'expected row "'//label//'" but found "'//trim(work(1:p-1))//'"')
    end if
    work = work(p+1:)
    do n = 1, size(out)
       q = index(work, ',')
       if (q == 0) q = len_trim(work) + 1
       if (len_trim(work(1:q-1)) == 0) call csv_err(                     &
          'row "'//label//'" has fewer values than grid blocks')
       read(work(1:q-1), *, iostat=io_) out(n)
       if (io_ /= 0) call csv_err('cannot read a value in row "'//label//'"')
       work = work(q+1:)
    end do
  end subroutine parse_row

  subroutine csv_err(msg)
    character(len=*), intent(in) :: msg
    if (myrank == 0) write(*,'(A)') ' ERROR (csv): '//msg
    call stop_report_input('Input_polyMesh_probCFD.csv: '//msg)
  end subroutine csv_err
end subroutine read_split_csv

!-----------------------------------------------------------------------
! Read a MULTI-block ASCII PLOT3D grid 'casename.grd'.  Node
! coordinates (meters) are kept on every rank: the curvilinear metrics
! of each owned sub-domain are built from them (curvi.F90), the
! reference pressure cell is located on rank 0, and rank 0 writes the
! PLOT3D grid file.
!-----------------------------------------------------------------------
subroutine read_grid
  use global_vars
  use mpi_vars
  implicit none
  integer :: iu, ios, ierr, b, i, j, k
  character(len=128) :: fname
  integer :: refdat(4)

  fname = trim(casename)//'.grd'
  open(newunit=iu, file=trim(fname), status='old', action='read', iostat=ios)
  if (ios /= 0) then
     if (myrank == 0) write(*,*) ' ERROR: cannot open grid file ', trim(fname)
     call MPI_Abort(MPI_COMM_WORLD, 13, ierr)
  end if
  read(iu,*) nblk
  if (nblk < 1) then
     if (myrank == 0) write(*,*) ' ERROR: bad number of blocks:', nblk
     call MPI_Abort(MPI_COMM_WORLD, 14, ierr)
  end if
  allocate(blk(nblk))
  read(iu,*) (blk(b)%ni, blk(b)%nj, blk(b)%nk, b=1,nblk)

  if (p_ref_switch) then
     p_ref_blk = 1;  p_ref_i = 1;  p_ref_j = 1;  p_ref_k = 1
  end if

  do b = 1, nblk
     associate(ni_=>blk(b)%ni, nj_=>blk(b)%nj, nk_=>blk(b)%nk)
       if (min(ni_, nj_, nk_) < 3) then
          if (myrank == 0) write(*,'(A,I0,A)') ' ERROR: block ', b,      &
             ' has fewer than 3 nodes in a direction'
          call MPI_Abort(MPI_COMM_WORLD, 14, ierr)
       end if
       blk(b)%nx = ni_ - 1;  blk(b)%ny = nj_ - 1;  blk(b)%nz = nk_ - 1
       allocate(blk(b)%xn(ni_,nj_,nk_), blk(b)%yn(ni_,nj_,nk_),           &
                blk(b)%zn(ni_,nj_,nk_))
       read(iu,*) (((blk(b)%xn(i,j,k), i=1,ni_), j=1,nj_), k=1,nk_)
       read(iu,*) (((blk(b)%yn(i,j,k), i=1,ni_), j=1,nj_), k=1,nk_)
       read(iu,*) (((blk(b)%zn(i,j,k), i=1,ni_), j=1,nj_), k=1,nk_)
       if (grid_scale /= 1.0_dp) then              ! Grid_scale: .grd units -> metres
          blk(b)%xn = blk(b)%xn*grid_scale;  blk(b)%yn = blk(b)%yn*grid_scale;  blk(b)%zn = blk(b)%zn*grid_scale
       end if
     end associate
  end do
  close(iu)

  if (myrank == 0) then
     write(*,'(A,I0)') '  grid blocks     : ', nblk
     write(*,'(A,I0,A,I0,A,I0,A)') '  block 1 (nodes) : ',                &
          blk(1)%ni,' x ',blk(1)%nj,' x ',blk(1)%nk
  end if

  ! locate the reference pressure cell (search on rank 0, broadcast)
  if (p_ref_switch) then
     if (myrank == 0) call find_ref_cell
     refdat = [p_ref_blk, p_ref_i, p_ref_j, p_ref_k]
     call MPI_Bcast(refdat, 4, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     p_ref_blk = refdat(1);  p_ref_i = refdat(2)
     p_ref_j   = refdat(3);  p_ref_k = refdat(4)
  end if

contains

  !---------------------------------------------------------------------
  ! Cell whose center (in meters) is closest to the reference pressure
  ! location, searched over all blocks (rank 0; node coords available).
  !---------------------------------------------------------------------
  subroutine find_ref_cell
    integer  :: b_, ic, jc, kc
    real(dp) :: cx, cy, cz, d2, d2min
    real(dp) :: cb(3), tolr
    logical  :: take
    d2min = huge(1.0_dp);  cb = huge(1.0_dp);  tolr = 0.0_dp
    do b_ = 1, nblk
      associate(bb=>blk(b_))
        do kc = 1, bb%nz
         do jc = 1, bb%ny
          do ic = 1, bb%nx
             cx = 0.125_dp*( bb%xn(ic,jc,kc)     + bb%xn(ic+1,jc,kc)      &
                           + bb%xn(ic,jc+1,kc)   + bb%xn(ic+1,jc+1,kc)    &
                           + bb%xn(ic,jc,kc+1)   + bb%xn(ic+1,jc,kc+1)    &
                           + bb%xn(ic,jc+1,kc+1) + bb%xn(ic+1,jc+1,kc+1) )
             cy = 0.125_dp*( bb%yn(ic,jc,kc)     + bb%yn(ic+1,jc,kc)      &
                           + bb%yn(ic,jc+1,kc)   + bb%yn(ic+1,jc+1,kc)    &
                           + bb%yn(ic,jc,kc+1)   + bb%yn(ic+1,jc,kc+1)    &
                           + bb%yn(ic,jc+1,kc+1) + bb%yn(ic+1,jc+1,kc+1) )
             cz = 0.125_dp*( bb%zn(ic,jc,kc)     + bb%zn(ic+1,jc,kc)      &
                           + bb%zn(ic,jc+1,kc)   + bb%zn(ic+1,jc+1,kc)    &
                           + bb%zn(ic,jc,kc+1)   + bb%zn(ic+1,jc,kc+1)    &
                           + bb%zn(ic,jc+1,kc+1) + bb%zn(ic+1,jc+1,kc+1) )
             d2 = (cx-p_ref_x)**2 + (cy-p_ref_y)**2 + (cz-p_ref_z)**2
             ! deterministic, orientation-independent choice: nearest
             ! cell center; ties broken by the smallest physical center
             ! coordinates (lexicographic in x, y, z)
             tolr = 1.0e-10_dp * max(1.0_dp, abs(d2), abs(d2min))
             take = .false.
             if (d2 < d2min - tolr) then
                take = .true.
             else if (d2 < d2min + tolr) then
                if (cx < cb(1) - tolr) then
                   take = .true.
                else if (cx < cb(1) + tolr) then
                   if (cy < cb(2) - tolr) then
                      take = .true.
                   else if (cy < cb(2) + tolr) then
                      if (cz < cb(3) - tolr) take = .true.
                   end if
                end if
             end if
             if (take) then
                d2min = d2;  cb = [cx, cy, cz]
                p_ref_blk = b_;  p_ref_i = ic;  p_ref_j = jc;  p_ref_k = kc
             end if
          end do
         end do
        end do
      end associate
    end do
    write(*,'(A,3ES11.3,A)') '  reference p location [m] :',              &
         p_ref_x, p_ref_y, p_ref_z
    write(*,'(A,I0,A,I0,A,I0,A,I0,A,ES12.4,A)')                           &
         '  -> pinned at block ', p_ref_blk, ' cell (', p_ref_i, ',',     &
         p_ref_j, ',', p_ref_k, ') to ', p_ref_val, ' Pa'
  end subroutine find_ref_cell
end subroutine read_grid

!-----------------------------------------------------------------------
! Read the boundary condition file 'casename.bcs' (multi-block).
!
! Format (comment lines starting with '!' and blank lines are skipped;
! text after the numbers on a line is ignored):
!
!   nblk               ! number of grid blocks (must match the grid)
!   --- for each block ---
!   nbc                ! number of boundary surfaces of this block
!     --- for each boundary surface ---
!     d                ! direction normal to surface
!                      !   |d| = 1(i), 2(j), 3(k);  + = min face, - = max face
!     i1  i2           ! start, end NODE index in i-direction
!     j1  j2           ! start, end NODE index in j-direction
!     k1  k2           ! start, end NODE index in k-direction
!     flag  value      ! u : 1 = Dirichlet [m/s], 0 = Neumann (zero gradient),
!     flag  value      ! v      2 = outflow: zero gradient + global mass-flux
!     flag  value      ! w      correction (3 on the u line = periodic surface)
!     flag  value      ! p : value in [Pa] for Dirichlet
!     flag  value      ! alpha (OPTIONAL 5th line, detected automatically;
!                      !   used when Num_of_Phase = 2, ignored otherwise):
!                      !   0 = zero gradient, 1 = Dirichlet value (0..1),
!                      !   2 = still-water profile (liquid below
!                      !       z_free_surface, same rule as the initial
!                      !       condition; use it at inlets)
!     flag  value(s)   ! T (update 16; REQUIRED on every surface when
!                      !   SW_Energy = 1, after the alpha line - which is
!                      !   then required too in a two-phase run; ignored
!                      !   when SW_Energy = 0):
!                      !   0 = zero gradient (adiabatic wall, outlet, symmetry)
!                      !   1 = Dirichlet T [K] (inlet, isothermal wall)
!                      !   3 = prescribed heat flux q_w [W/m2] INTO the fluid
!                      !   4 = convective (Robin): "4  h [W/(m2 K)]  T_inf [K]"
!                      !       (two values on the line)
!
! Dirichlet velocity values are GLOBAL CARTESIAN components [m/s] on
! every block; the solver projects them onto the transformed face
! variables with the local metric matrix.  Block-to-block interfaces
! must NOT be listed here; they are detected automatically and handled
! through 'casename.bcs_b2b'.
!-----------------------------------------------------------------------
subroutine read_bcs
  use global_vars
  use mpi_vars
  use b2b_mod, only: b2b_read_file, n_pole, pole_blk, pole_dir, pole_rng
  implicit none
  integer :: iu, ios, ierr, nblk_f, b, n, m, nbtot, nvar, n_alpha_lines
  logical :: eof
  integer :: nbf, nb, d, fl1, rr(2,3), npb, ip, nbp, tolb, nxl, xfl(2), n_t_lines, n_wave_patches
  real(dp) :: tol_file, xv1(2), xv2(2)
  character(len=256) :: line
  character(len=128) :: fname

  ! pole (degenerate) faces found by the interface search: they get an
  ! automatic zero-gradient / zero-flux patch each (no .bcs record)
  call b2b_read_file(trim(casename)//'.bcs_b2b', tolb, tol_file, ierr)
  if (ierr /= 0) call MPI_Abort(MPI_COMM_WORLD, 20, ierr)
  do b = 1, nblk
     blk(b)%pole = .false.
  end do
  nbp = 0
  do ip = 1, n_pole
     b = pole_blk(ip)
     blk(b)%pole(merge(1, 2, pole_dir(ip) > 0), abs(pole_dir(ip))) = .true.
  end do

  fname = trim(casename)//'.bcs'
  open(newunit=iu, file=trim(fname), status='old', action='read', iostat=ios)
  if (ios /= 0) then
     if (myrank == 0) write(*,*) ' ERROR: cannot open bc file ', trim(fname)
     call MPI_Abort(MPI_COMM_WORLD, 15, ierr)
  end if

  call next_line(iu, line);  read(line,*) nblk_f
  if (nblk_f /= nblk) then
     if (myrank == 0) write(*,'(A,I0,A,I0)') ' ERROR: bc file has ',      &
        nblk_f, ' blocks but the grid has ', nblk
     call MPI_Abort(MPI_COMM_WORLD, 16, ierr)
  end if

  nvar = merge(5, 4, Num_of_Phase == 2)     ! u,v,w,p (+ alpha) (+ T, row 6, energy_on)
  n_alpha_lines = 0;  n_t_lines = 0;  n_wave_patches = 0
  nbtot = 0
  nper  = 0
  if (allocated(pbc)) deallocate(pbc)
  allocate(pbc(8))
  p_has_dirichlet = .false.
  do b = 1, nblk
     call next_line(iu, line);  read(line,*) nbf
     ! nbf surfaces follow; Dirichlet/Neumann ones go into the bc arrays,
     ! periodic ones (flag 3) into the pbc(:) list; pole faces of this
     ! block are appended as automatic patches (npb of them)
     npb = count(pole_blk(1:n_pole) == b)
     allocate(blk(b)%bc_dir(max(nbf+npb,1)), blk(b)%bc_irng(2,max(nbf+npb,1)),    &
              blk(b)%bc_jrng(2,max(nbf+npb,1)), blk(b)%bc_krng(2,max(nbf+npb,1)))
     allocate(blk(b)%bc_flag(6,max(nbf+npb,1)), blk(b)%bc_val(6,max(nbf+npb,1)), blk(b)%bc_tinf(max(nbf+npb,1)))
     allocate(blk(b)%bc_eflag(3,max(nbf+npb,1)))
     blk(b)%bc_flag = 0;  blk(b)%bc_val = 0.0_dp;  blk(b)%bc_eflag = 0;  blk(b)%bc_tinf = 0.0_dp
     nb = 0
     do n = 1, nbf
        call next_line(iu, line);  read(line,*) d
        call next_line(iu, line);  read(line,*) rr(:,1)
        call next_line(iu, line);  read(line,*) rr(:,2)
        call next_line(iu, line);  read(line,*) rr(:,3)
        call next_line(iu, line);  read(line,*) fl1
        if (fl1 == 3) then
           !--------------- periodic boundary surface (flag 3) ---------
           ! record layout after the flag line:
           !   pair block / pair face dir / pair i, j, k node ranges
           nper = nper + 1
           call grow_pbc(nper)
           pbc(nper)%rblk = b
           pbc(nper)%dirs = d
           pbc(nper)%rn   = rr
           call next_line(iu, line);  read(line,*) pbc(nper)%dblk
           call next_line(iu, line);  read(line,*) pbc(nper)%ddir
           call next_line(iu, line);  read(line,*) pbc(nper)%dn(:,1)
           call next_line(iu, line);  read(line,*) pbc(nper)%dn(:,2)
           call next_line(iu, line);  read(line,*) pbc(nper)%dn(:,3)
        else
           !--------------- Dirichlet / Neumann surface -----------------
           nb = nb + 1
           blk(b)%bc_dir(nb)    = d
           blk(b)%bc_irng(:,nb) = rr(:,1)
           blk(b)%bc_jrng(:,nb) = rr(:,2)
           blk(b)%bc_krng(:,nb) = rr(:,3)
           read(line,*) blk(b)%bc_flag(1,nb), blk(b)%bc_val(1,nb)
           do m = 2, 4
              call next_line(iu, line)
              read(line,*) blk(b)%bc_flag(m,nb), blk(b)%bc_val(m,nb)
           end do
           ! optional lines after p (alpha, T): a line that holds TWO or
           ! more numbers ("flag value [value]") is a variable line; a
           ! single number is the direction of the next surface or the
           ! surface count of the next block.  Up to two extra lines:
           !   two-phase run : alpha [, T]      (T required with SW_Energy = 1)
           !   single-phase  : T  or  alpha (ignored), T  (energy on);
           !                   alpha (ignored)         (energy off)
           nxl = 0
           do
              call next_line(iu, line, eof)
              if (eof) exit
              if (ntok(line) < 2 .or. nxl == 2) then
                 backspace(iu);  exit
              end if
              nxl = nxl + 1
              xfl(nxl) = 0;  xv1(nxl) = 0.0_dp;  xv2(nxl) = 0.0_dp
              if (ntok(line) >= 3) then
                 read(line,*) xfl(nxl), xv1(nxl), xv2(nxl)
              else
                 read(line,*) xfl(nxl), xv1(nxl)
              end if
           end do
           if (energy_on) then
              if (Num_of_Phase == 2) then
                 if (nxl < 2) call bcs_err(b, n, 'SW_Energy = 1 and Num_of_Phase = 2: every surface needs an alpha line AND a T line after p')
                 blk(b)%bc_flag(5,nb) = xfl(1);  blk(b)%bc_val(5,nb) = xv1(1);  n_alpha_lines = n_alpha_lines + 1
                 blk(b)%bc_flag(6,nb) = xfl(2);  blk(b)%bc_val(6,nb) = xv1(2);  blk(b)%bc_tinf(nb) = xv2(2)
              else
                 if (nxl < 1) call bcs_err(b, n, 'SW_Energy = 1: every surface needs a T line after p')
                 if (nxl == 2) n_alpha_lines = n_alpha_lines + 1
                 blk(b)%bc_flag(6,nb) = xfl(nxl);  blk(b)%bc_val(6,nb) = xv1(nxl);  blk(b)%bc_tinf(nb) = xv2(nxl)
              end if
              n_t_lines = n_t_lines + 1
           else
              if (nxl >= 1) then
                 blk(b)%bc_flag(5,nb) = xfl(1);  blk(b)%bc_val(5,nb) = xv1(1);  n_alpha_lines = n_alpha_lines + 1
              end if
              if (nxl == 2) n_t_lines = n_t_lines + 1
           end if
           call check_flags(b, nb, nvar)
        end if
     end do
     ! automatic pole patches: zero gradient for u, v, w, p, alpha; the
     ! face flux is zero through the vanishing face Jacobian
     do ip = 1, n_pole
        if (pole_blk(ip) /= b) cycle
        do n = 1, nb
           if (blk(b)%bc_dir(n) == pole_dir(ip)) then
              if (myrank == 0) write(*,'(A,I0,A,I0,A)') '  WARNING: block ', b, &
                 ' face dir ', pole_dir(ip), ' is a pole face but is also listed'// &
                 ' in the .bcs file; the .bcs patch is applied first, the pole'//  &
                 ' patch (zero gradient) last and overrides it'
           end if
        end do
        nb = nb + 1;  nbp = nbp + 1
        blk(b)%bc_dir(nb)    = pole_dir(ip)
        blk(b)%bc_irng(:,nb) = pole_rng(:,1,ip)
        blk(b)%bc_jrng(:,nb) = pole_rng(:,2,ip)
        blk(b)%bc_krng(:,nb) = pole_rng(:,3,ip)
        blk(b)%bc_flag(:,nb) = 0
        blk(b)%bc_val(:,nb)  = 0.0_dp;  blk(b)%bc_tinf(nb) = 0.0_dp
     end do
     blk(b)%nbc = nb
     nbtot = nbtot + nb
     ! non-dimensionalize: velocity by Uch, pressure by rho1*Uch^2
     do n = 1, nb
        blk(b)%bc_val(1:3,n) = blk(b)%bc_val(1:3,n) / Uch
        blk(b)%bc_val(4,n)   = blk(b)%bc_val(4,n)   / (rho1*Uch*Uch)
     end do
     if (nb > 0) then
        if (any(blk(b)%bc_flag(4,1:nb) == 1)) p_has_dirichlet = .true.
     end if
  end do
  close(iu)

  ! effective flags of the transformed face families on every patch
  call effective_flags

  ! resolve and check the periodic pairings (axis map, directions)
  if (nper > 0) call resolve_periodic

  ! consistency between the bc file and the reference pressure switch
  if (.not. p_has_dirichlet .and. .not. p_ref_switch) then
     if (myrank == 0) then
        write(*,'(A)') ' ERROR: all pressure BCs are Neumann (zero gradient), so the'
        write(*,'(A)') '        pressure level is undetermined. Set p_ref_switch = .true.'
        write(*,'(A)') '        in probCFD.nml and give p_ref_x/y/z [m] and p_ref_val [Pa].'
     end if
     call MPI_Abort(MPI_COMM_WORLD, 18, ierr)
  end if
  if (p_has_dirichlet .and. p_ref_switch) then
     if (myrank == 0) then
        write(*,'(A)') '  WARNING: a Dirichlet pressure BC exists; the reference'
        write(*,'(A)') '           pressure (p_ref_switch) is ignored.'
     end if
  end if

  if (myrank == 0) then
     if (nbp > 0) write(*,'(A,I0,A)') '  pole faces      : ', nbp,       &
          ' (automatic zero-gradient / zero-flux axis patches, from .bcs_b2b)'
     write(*,'(A,I0,A)') '  boundary faces  : ', nbtot,                   &
          ' (Dirichlet/Neumann, all blocks)'
     if (Num_of_Phase == 2 .and. n_alpha_lines < nbtot) write(*,'(A)')   &
          '  NOTE: two-phase run but the bcs file has no alpha line on'// &
          ' some surfaces: zero-gradient alpha assumed there'
     if (Num_of_Phase == 1 .and. n_alpha_lines > 0) write(*,'(A)')       &
          '  NOTE: single-phase run: the alpha lines of the bcs file are ignored'
     if (energy_on) write(*,'(A,I0,A)') '  temperature BCs : ', n_t_lines, &
          ' T lines (flags 0 adiabatic / 1 Dirichlet / 3 heat flux / 4 convective)'
     if (.not. energy_on .and. n_t_lines > 0) write(*,'(A)') '  NOTE: SW_Energy = 0: the T lines of the bcs file are ignored'
     if (nper > 0) write(*,'(A,I0,A)') '  periodic faces  : ', nper,      &
          ' (flag 3; paired through the exchange layer)'
     if (n_wave_patches > 0) write(*,'(A,I0,A)') '  wave patches    : ', n_wave_patches, &
          ' (flag 12: incident wave on u, v, w, p, alpha; mean velocity from the u, v, w values)'
  end if

contains

  !--------------------------------------------------------------------
  ! The .bcs flags are given per CARTESIAN component (u, v, w), while
  ! the solver's face unknowns are the transformed components U^q along
  ! the grid directions.  On a block whose axes are not aligned with
  ! x,y,z (e.g. an O-grid whose i-axis is vertical) the condition that
  ! applies to U^q is the one of the Cartesian axis the grid direction
  ! q points along.  For every patch the three grid directions are
  ! evaluated at the patch centre from the node coordinates (the face
  ! normal for the normal family, the grid-line tangents for the two
  ! tangential families) and the dominant Cartesian axis of each
  ! selects the flag: bc_eflag(q) = bc_flag(axis(q)).
  ! On axis-aligned planar patches this is exact; on curved patches it
  ! holds for uniform flags (all Dirichlet / all zero-gradient) and is
  ! a patch-averaged choice for mixed (slip-type) conditions.
  !--------------------------------------------------------------------
  subroutine effective_flags
    integer :: b_, n_, na, a, ic, jc, kc, q, i0, i1, j0, j1, k0, k1
    real(dp) :: t(3,3), nrm(3), v(3)
    integer :: dom_ax(3)
    do b_ = 1, nblk
       associate(bb => blk(b_))
         do n_ = 1, bb%nbc
            na = abs(bb%bc_dir(n_))
            ic = (bb%bc_irng(1,n_) + bb%bc_irng(2,n_))/2
            jc = (bb%bc_jrng(1,n_) + bb%bc_jrng(2,n_))/2
            kc = (bb%bc_krng(1,n_) + bb%bc_krng(2,n_))/2
            ! grid-line tangents at the patch centre (one-sided at the
            ! block boundary)
            i0 = max(ic-1,1);  i1 = min(ic+1,bb%ni)
            j0 = max(jc-1,1);  j1 = min(jc+1,bb%nj)
            k0 = max(kc-1,1);  k1 = min(kc+1,bb%nk)
            t(:,1) = [bb%xn(i1,jc,kc)-bb%xn(i0,jc,kc), bb%yn(i1,jc,kc)-bb%yn(i0,jc,kc), &
                      bb%zn(i1,jc,kc)-bb%zn(i0,jc,kc)]
            t(:,2) = [bb%xn(ic,j1,kc)-bb%xn(ic,j0,kc), bb%yn(ic,j1,kc)-bb%yn(ic,j0,kc), &
                      bb%zn(ic,j1,kc)-bb%zn(ic,j0,kc)]
            t(:,3) = [bb%xn(ic,jc,k1)-bb%xn(ic,jc,k0), bb%yn(ic,jc,k1)-bb%yn(ic,jc,k0), &
                      bb%zn(ic,jc,k1)-bb%zn(ic,jc,k0)]
            do q = 1, 3
               if (q == na) then
                  ! face normal = cross product of the two tangential tangents
                  select case (na)
                  case (1);  nrm = cross(t(:,2), t(:,3))
                  case (2);  nrm = cross(t(:,3), t(:,1))
                  case default;  nrm = cross(t(:,1), t(:,2))
                  end select
                  v = nrm
               else
                  v = t(:,q)
               end if
               dom_ax(q) = maxloc(abs(v), 1)
               bb%bc_eflag(q,n_) = bb%bc_flag(dom_ax(q), n_)
            end do
            do a = 1, 3
               if (bb%bc_flag(a,n_) == 3) bb%bc_eflag(:,n_) = 3   ! periodic marker (unused here)
               if (bb%bc_flag(a,n_) == 11) then
                  if (.not. all(bb%bc_flag(1:3,n_) == 11)) then
                     if (myrank == 0) write(*,'(A,I0,A,I0,A,3I4)') ' ERROR (bcs): block ', b_, ' surface ', n_, &
                        ' - flag 11 (no-slip wall + wall function) must be given for u, v AND w; flags found:', bb%bc_flag(1:3,n_)
                     call MPI_Abort(MPI_COMM_WORLD, 61, ierr)
                  end if
                  wallfn_on = .true.
               end if
            end do
         end do
       end associate
    end do
  end subroutine effective_flags

  pure function cross(a, b) result(c)
    real(dp), intent(in) :: a(3), b(3)
    real(dp) :: c(3)
    c = [a(2)*b(3)-a(3)*b(2), a(3)*b(1)-a(1)*b(3), a(1)*b(2)-a(2)*b(1)]
  end function cross

  !--------------------------------------------------------------------
  subroutine check_flags(b_, n_, nv)
    integer, intent(in) :: b_, n_, nv
    integer :: m_, ie
    ! flag 11 (no-slip wall + wall function) is valid for u, v, w only and
    ! must be given for all three components together
    if (any(blk(b_)%bc_flag(1:3,n_) == 11) .and. .not. all(blk(b_)%bc_flag(1:3,n_) == 11)) then
       if (myrank == 0) write(*,'(A,I0,A,I0,A,3I4)') ' ERROR (bcs): block ', b_, ' surface ', n_, &
          ' - flag 11 (no-slip wall + wall function) must be set for u, v AND w together; flags found (u v w):', &
          blk(b_)%bc_flag(1:3,n_)
       call MPI_Abort(MPI_COMM_WORLD, 61, ie)
    end if
    if (all(blk(b_)%bc_flag(1:3,n_) == 11)) wallfn_on = .true.
    ! flag 12 (incident wave, update 18): u, v, w, p AND alpha together (two-phase only)
    if (any(blk(b_)%bc_flag(1:5,n_) == 12)) then
       if (nv < 5) call bcs_err(b_, n_, 'flag 12 (incident wave) needs a two-phase run (Num_of_Phase = 2) with an alpha line')
       if (.not. all(blk(b_)%bc_flag(1:5,n_) == 12)) call bcs_err(b_, n_, &
          'flag 12 (incident wave) must be set on the u, v, w, p AND alpha lines together')
       n_wave_patches = n_wave_patches + 1
    end if
    do m_ = 1, nv
       if (m_ <= 3 .and. blk(b_)%bc_flag(m_,n_) == 11) cycle
       if (m_ <= 5 .and. blk(b_)%bc_flag(m_,n_) == 12) cycle
       if (blk(b_)%bc_flag(m_,n_) < 0 .or. blk(b_)%bc_flag(m_,n_) > 2 .or. &
           (m_ == 4 .and. blk(b_)%bc_flag(m_,n_) == 2)) then
          if (myrank == 0) write(*,'(A,I0,A,I0,A,I0,A,I0)')              &
             ' ERROR (bcs): block ', b_, ' surface ', n_, ' variable ',  &
             m_, ' has an invalid flag ', blk(b_)%bc_flag(m_,n_)
          call MPI_Abort(MPI_COMM_WORLD, 19, ie)
       end if
    end do
    if (nv == 5) then
       if (blk(b_)%bc_flag(5,n_) == 1 .and. (blk(b_)%bc_val(5,n_) < 0.0_dp &
           .or. blk(b_)%bc_val(5,n_) > 1.0_dp)) then
          if (myrank == 0) write(*,'(A,I0,A,I0,A)') ' ERROR (bcs): block ', &
             b_, ' surface ', n_, ': Dirichlet alpha must lie in [0,1]'
          call MPI_Abort(MPI_COMM_WORLD, 19, ie)
       end if
    end if
    if (energy_on) then                      ! T line (row 6): flags 0, 1, 3, 4
       select case (blk(b_)%bc_flag(6,n_))
       case (0, 3)
       case (1)
          if (blk(b_)%bc_val(6,n_) <= 0.0_dp) call bcs_err(b_, n_, 'Dirichlet T must be > 0 K')
       case (4)
          if (blk(b_)%bc_val(6,n_) < 0.0_dp) call bcs_err(b_, n_, 'convective T BC: h must be >= 0')
          if (blk(b_)%bc_tinf(n_) <= 0.0_dp) call bcs_err(b_, n_, 'convective T BC needs "4  h  T_inf" with T_inf > 0 K')
       case default
          call bcs_err(b_, n_, 'invalid T flag (0 zero gradient, 1 Dirichlet, 3 heat flux, 4 convective)')
       end select
    end if
  end subroutine check_flags

  subroutine bcs_err(b_, n_, msg)
    integer, intent(in) :: b_, n_
    character(len=*), intent(in) :: msg
    integer :: ie
    if (myrank == 0) write(*,'(A,I0,A,I0,A)') ' ERROR (bcs): block ', b_, ' surface ', n_, ': '//msg
    flush(6)
    call MPI_Abort(MPI_COMM_WORLD, 19, ie)
  end subroutine bcs_err

  !--------------------------------------------------------------------
  subroutine grow_pbc(need)
    integer, intent(in) :: need
    type(pbc_t), allocatable :: tmp(:)
    if (need <= size(pbc)) return
    allocate(tmp(2*size(pbc)))
    tmp(1:size(pbc)) = pbc
    call move_alloc(tmp, pbc)
  end subroutine grow_pbc

  !--------------------------------------------------------------------
  ! Resolve the axis map of every periodic patch and validate it.
  !
  ! The two tangential axes of the local face are matched to the two
  ! tangential axes of the pair face by their index-range EXTENTS (the
  ! running direction comes from the ordering of the pair range, e.g.
  ! "10 1" = opposite direction).  When both tangential extents are
  ! equal (square faces) the pairing is ambiguous from the ranges
  ! alone; in that case the grid coordinates decide: for a (rigid,
  ! translational) periodic pairing the vector from each local node to
  ! its pair node must be THE SAME across the whole surface, and only
  ! one of the two candidate pairings satisfies this.  The geometric
  ! test runs on rank 0 and the result is broadcast.
  !--------------------------------------------------------------------
  subroutine resolve_periodic
    integer :: ip, na, qd, t(2), s(2), eL(2), eP(2), a, k2
    integer :: opt, nok, okopt, qtry(3,2), styr(3,2)
    integer, allocatable :: qs(:,:)
    real(dp), allocatable :: tvs(:,:)
    logical :: fits(2), geo_ok
    integer :: ierr2

    do ip = 1, nper
       associate(pp => pbc(ip))
         na = abs(pp%dirs);  qd = abs(pp%ddir)
         !---- basic validation ---------------------------------------
         call chk(pp%rblk >= 1 .and. pp%rblk <= nblk, ip, 'local block out of range')
         call chk(pp%dblk >= 1 .and. pp%dblk <= nblk, ip, 'pair block out of range')
         call chk(na >= 1 .and. na <= 3 .and. qd >= 1 .and. qd <= 3, ip, &
                  'face direction must be +-1, +-2 or +-3')
         call chk(pp%rn(1,na) == pp%rn(2,na), ip, 'local normal range must be fixed')
         call chk(pp%rn(1,na) == merge(1, nnode_b(pp%rblk,na), pp%dirs > 0), ip, &
                  'local range does not lie on the given face')
         call chk(pp%dn(1,qd) == pp%dn(2,qd), ip, 'pair normal range must be fixed')
         call chk(pp%dn(1,qd) == merge(1, nnode_b(pp%dblk,qd), pp%ddir > 0), ip, &
                  'pair range does not lie on the given pair face')
         !---- tangential axes and extents -----------------------------
         k2 = 0
         do a = 1, 3
            if (a /= na) then
               k2 = k2 + 1;  t(k2) = a
               call chk(pp%rn(1,a) <= pp%rn(2,a), ip,                    &
                        'local tangential ranges must be ascending')
            end if
         end do
         k2 = 0
         do a = 1, 3
            if (a /= qd) then
               k2 = k2 + 1;  s(k2) = a
            end if
         end do
         eL = pp%rn(2,t) - pp%rn(1,t) + 1
         eP = abs(pp%dn(2,s) - pp%dn(1,s)) + 1
         !---- the two candidate pairings ------------------------------
         fits(1) = (eL(1) == eP(1) .and. eL(2) == eP(2))   ! t1-s1, t2-s2
         fits(2) = (eL(1) == eP(2) .and. eL(2) == eP(1))   ! t1-s2, t2-s1
         call chk(fits(1) .or. fits(2), ip,                              &
                  'pair surface extents do not match the local surface')
         do opt = 1, 2
            qtry(na,opt) = qd
            styr(na,opt) = merge(-1, 1, pp%dirs > 0) * merge(1, -1, pp%ddir > 0)
            do k2 = 1, 2
               a = t(k2)
               qtry(a,opt) = merge(s(k2), s(3-k2), opt == 1)
               styr(a,opt) = rng_sgn(pp%dn(:, qtry(a,opt)))
            end do
         end do
         !---- choose the pairing ---------------------------------------
         nok = count(fits)
         pp%tv = 0.0_dp
         if (nok == 1) then
            okopt = merge(1, 2, fits(1))
            if (myrank == 0) then
               call geo_check(ip, qtry(:,okopt), styr(:,okopt), geo_ok, pp%tv)
               if (.not. geo_ok) write(*,'(A,I0,A)')                     &
                  '  WARNING (periodic): patch ', ip, ' - the node-to-'//&
                  'node translation is not uniform; check the bcs setup'
            end if
         else
            ! square face: rank 0 decides geometrically
            okopt = 1
            if (myrank == 0) then
               call geo_check(ip, qtry(:,1), styr(:,1), geo_ok, pp%tv)
               if (.not. geo_ok) then
                  call geo_check(ip, qtry(:,2), styr(:,2), geo_ok, pp%tv)
                  if (geo_ok) then
                     okopt = 2
                  else
                     write(*,'(A,I0,A)') '  WARNING (periodic): patch ', &
                        ip, ' - no axis pairing gives a uniform node-'// &
                        'to-node translation; using the ascending-'//    &
                        'axis-order pairing'
                  end if
               end if
            end if
         end if
         pp%q   = qtry(:,okopt)
         pp%sgn = styr(:,okopt)
       end associate
    end do

    ! rank 0 made the (possibly geometric) choices: broadcast them
    allocate(qs(6, nper), tvs(3, nper))
    do ip = 1, nper
       qs(1:3,ip) = pbc(ip)%q
       qs(4:6,ip) = pbc(ip)%sgn
       tvs(:,ip)  = pbc(ip)%tv
    end do
    call MPI_Bcast(qs, 6*nper, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr2)
    call MPI_Bcast(tvs, 3*nper, MPI_DOUBLE_PRECISION, 0,                 &
                   MPI_COMM_WORLD, ierr2)
    do ip = 1, nper
       pbc(ip)%q   = qs(1:3,ip)
       pbc(ip)%sgn = qs(4:6,ip)
       pbc(ip)%tv  = tvs(:,ip)
    end do
    deallocate(qs, tvs)
  end subroutine resolve_periodic

  !--------------------------------------------------------------------
  pure integer function rng_sgn(r)
    integer, intent(in) :: r(2)
    rng_sgn = merge(-1, 1, r(2) < r(1))
  end function rng_sgn

  !--------------------------------------------------------------------
  pure integer function nnode_b(b_, a_)
    integer, intent(in) :: b_, a_
    select case (a_)
    case (1);  nnode_b = blk(b_)%ni
    case (2);  nnode_b = blk(b_)%nj
    case default;  nnode_b = blk(b_)%nk
    end select
  end function nnode_b

  !--------------------------------------------------------------------
  subroutine chk(cond, ip, msg)
    logical, intent(in) :: cond
    integer, intent(in) :: ip
    character(len=*), intent(in) :: msg
    integer :: ie
    if (cond) return
    if (myrank == 0) write(*,'(A,I0,A)') ' ERROR (periodic bc, patch ',  &
       ip, '): '//msg
    call MPI_Abort(MPI_COMM_WORLD, 40, ie)
  end subroutine chk

  !--------------------------------------------------------------------
  ! Rank-0 geometric test of one candidate axis map: the translation
  ! vector from every local surface node to its mapped pair node must
  ! be uniform over the surface (within a spacing-scaled tolerance).
  !--------------------------------------------------------------------
  subroutine geo_check(ip, qm, sm, ok, tvout)
    integer, intent(in)  :: ip, qm(3), sm(3)
    logical, intent(out) :: ok
    real(dp), intent(out) :: tvout(3)
    integer  :: na, a, i, j, k, li(3), di(3), noffs(3)
    real(dp) :: t0(3), tv(3), tolp, spanm
    logical  :: first

    associate(pp => pbc(ip), rb => blk(pbc(ip)%rblk), db => blk(pbc(ip)%dblk))
      na = abs(pp%dirs)
      do a = 1, 3
         noffs(a) = pp%dn(1,qm(a)) - sm(a)*pp%rn(1,a)
      end do
      spanm = max( maxval(rb%xn)-minval(rb%xn),                          &
                   maxval(rb%yn)-minval(rb%yn),                          &
                   maxval(rb%zn)-minval(rb%zn) )
      tolp = 1.0e-6_dp * max(spanm, 1.0e-30_dp)
      ok = .true.
      first = .true.
      tvout = 0.0_dp
      do k = pp%rn(1,3), pp%rn(2,3)
       do j = pp%rn(1,2), pp%rn(2,2)
        do i = pp%rn(1,1), pp%rn(2,1)
           li = [i, j, k]
           do a = 1, 3
              di(qm(a)) = sm(a)*li(a) + noffs(a)
           end do
           tv = [ db%xn(di(1),di(2),di(3)) - rb%xn(i,j,k),               &
                  db%yn(di(1),di(2),di(3)) - rb%yn(i,j,k),               &
                  db%zn(di(1),di(2),di(3)) - rb%zn(i,j,k) ]
           if (first) then
              t0 = tv;  tvout = tv;  first = .false.
           else if (any(abs(tv - t0) > tolp)) then
              ok = .false.
              return
           end if
        end do
       end do
      end do
    end associate
  end subroutine geo_check

  ! number of numeric tokens before any '!' comment
  integer function ntok(s)
    character(len=*), intent(in) :: s
    integer :: p, l
    logical :: intok
    p = index(s, '!');  l = len_trim(s)
    if (p > 0) l = p - 1
    ntok = 0;  intok = .false.
    do p = 1, l
       if (s(p:p) == ' ' .or. s(p:p) == ',' .or. s(p:p) == char(9)) then
          intok = .false.
       else if (.not. intok) then
          intok = .true.;  ntok = ntok + 1
       end if
    end do
  end function ntok

  subroutine next_line(unit, out, at_eof)
    integer, intent(in) :: unit
    character(len=*), intent(out) :: out
    logical, intent(out), optional :: at_eof
    integer :: io
    if (present(at_eof)) at_eof = .false.
    do
       read(unit,'(A)', iostat=io) out
       if (io /= 0) then
          if (present(at_eof)) then
             at_eof = .true.;  out = ' '
             return
          end if
          if (myrank == 0) write(*,*) ' ERROR: unexpected end of bc file'
          call MPI_Abort(MPI_COMM_WORLD, 17, io)
       end if
       if (len_trim(out) == 0) cycle
       out = adjustl(out)
       if (out(1:1) == '!') cycle
       exit
    end do
  end subroutine next_line
end subroutine read_bcs

!-----------------------------------------------------------------------
! Runtime name-list updates (update 13): every time step rank 0 reads
! probCFD.nml again; when the file changed (content hash) the lines are
! broadcast and every rank applies the RUNTIME-CHANGEABLE variables:
!   Del_T, End_T, Iner_Itn, Time_scheme, the momentum / pressure / VoF /
!   turbulence / energy iteration controls (incl. Pr_t), Output_Intvl, Log_Intvl, Data_Intvl,
!   Cnvt_scheme, Cnvt_mode, Dffs_mode, cfl_switch, vn_switch,
!   vof_cfl_max, t_ramp, Stop_Rsd_Max, Stop_CFL_Max, poi_verbose,
!   Motion_Start_T, Wave_Ramp_T, Relax_In_L / Relax_Out_L / Relax_Side_L.
! Everything else (grid, phases, fluids, references, overset set-up,
! turbulence model, VoF scheme, 6-DoF body data) needs a restart and is
! ignored here.  The values that changed are written to the log; the
! number of remaining steps follows a new Del_T / End_T from the current
! time on.
!-----------------------------------------------------------------------
subroutine reread_namelist(istep)
  use global_vars
  use mpi_vars
  use sixdof, only: Motion_Start_T
  use waves, only: wave_runtime
  implicit none
  integer, intent(in) :: istep
  integer, parameter :: MXL = 4000, LL = 256
  character(len=LL), allocatable, save :: buf(:)
  integer, save :: hash_old = -1
  integer :: iu, ios, nl, hash, ierr, flag, l, pe, nv, i, c
  character(len=LL) :: line, work, name
  character(len=128) :: vals(8)
  character(len=LL) :: msg(40)
  integer :: nmsg
  real(dp) :: rv(8)
  integer  :: nstep_new
  logical  :: wchg

  if (.not. allocated(buf)) allocate(buf(MXL))
  flag = 0
  if (myrank == 0) then
     open(newunit=iu, file='probCFD.nml', status='old', action='read', iostat=ios)
     if (ios == 0) then
        nl = 0;  hash = 7
        do
           read(iu, '(A)', iostat=ios) line
           if (ios /= 0) exit
           if (nl < MXL) then
              nl = nl + 1;  buf(nl) = line
           end if
           do i = 1, len_trim(line)
              hash = mod(hash*31 + iachar(line(i:i)), 1000000007)
           end do
        end do
        close(iu)
        if (hash_old < 0) then
           hash_old = hash                        ! first call: the file as it was read at start-up
        else if (hash /= hash_old) then
           hash_old = hash;  flag = nl
        end if
     end if
  end if
  call MPI_Bcast(flag, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  if (flag == 0) return
  nl = flag
  call MPI_Bcast(buf, MXL*LL, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

  nmsg = 0
  do l = 1, nl
     work = buf(l)
     pe = index(work, char(13));  if (pe > 0) work(pe:) = ' '
     pe = index(work, '!');       if (pe > 0) work(pe:) = ' '
     work = adjustl(work)
     if (len_trim(work) == 0 .or. work(1:1) == '$' .or. work(1:1) == '/') cycle
     pe = index(work, '=')
     if (pe == 0) cycle
     name = lower2(trim(adjustl(work(1:pe-1))))
     call split2(work(pe+1:), vals, nv)
     if (nv == 0) cycle
     rv = 0.0_dp
     do i = 1, nv
        read(vals(i), *, iostat=ios) rv(i)
        if (ios /= 0) rv(i) = 0.0_dp
     end do
     select case (trim(name))
     case ('del_t', 'dt_dim')
        if (rv(1) > 0.0_dp .and. rv(1) /= dt_dim) then
           call note('Del_T', dt_dim, rv(1));  dt_dim = rv(1);  dt = dt_dim/tch
        end if
     case ('end_t')
        if (rv(1) > 0.0_dp .and. rv(1) /= End_T) then
           call note('End_T', End_T, rv(1));  End_T = rv(1)
        end if
     case ('iner_itn');           call seti('Iner_Itn', Iner_Itn, nint(rv(1)), 1, 100)
     case ('time_scheme');        call seti('Time_scheme', Time_scheme, nint(rv(1)), 1, 2)
     case ('max_itr_vel', 'mom_maxiter');  call seti('Max_Itr_Vel', mom_maxiter, nint(rv(1)), 1, 1000000)
     case ('min_itr_vel');        call seti('Min_Itr_Vel', mom_miniter, nint(rv(1)), 1, 1000000)
     case ('rsd_lmt_vel', 'mom_tol');      call setr('Rsd_Lmt_Vel', mom_tol, rv(1), 0.0_dp, 1.0e30_dp)
     case ('sor_fct_vel');        call setr('SOR_fct_Vel', mom_sor, rv(1), 0.0_dp, 1.0_dp)
     case ('max_itr_prs', 'poi_maxiter');  call seti('Max_Itr_Prs', poi_maxiter, nint(rv(1)), 1, 1000000)
     case ('min_itr_prs');        call seti('Min_Itr_Prs', poi_miniter, nint(rv(1)), 1, 1000000)
     case ('trl_lmt_prs', 'poi_tol');      call setr('Trl_Lmt_Prs', poi_tol, rv(1), 0.0_dp, 1.0e30_dp)
     case ('rsd_lmt_prs');        call setr('Rsd_Lmt_Prs', poi_rtol, rv(1), 0.0_dp, 1.0e30_dp)
     case ('sor_fct_prs');        call setr('SOR_fct_Prs', poi_sor, rv(1), 0.0_dp, 1.0_dp)
     case ('max_itr_vof');        call seti('Max_Itr_VoF', vof_maxiter, nint(rv(1)), 1, 1000000)
     case ('min_itr_vof');        call seti('Min_Itr_VoF', vof_miniter, nint(rv(1)), 1, 1000000)
     case ('rsd_lmt_vof');        call setr('Rsd_Lmt_VoF', vof_tol, rv(1), 0.0_dp, 1.0e30_dp)
     case ('sor_fct_vof');        call setr('SOR_fct_VoF', vof_sor, rv(1), 0.0_dp, 1.0_dp)
     case ('output_intvl', 'out_interval');  call seti('Output_Intvl', out_interval, nint(rv(1)), 1, 100000000)
    case ('sw_prs_frc_surf');               call seti('SW_Prs_Frc_surf', sw_prs_frc_surf, nint(rv(1)), 0, 1)
    case ('statistics_output_intvl');       call seti('Statistics_Output_Intvl', stat_output_intvl, nint(rv(1)), 0, 100000000)
    case ('frame_prs_frc_surf');            call seti('Frame_Prs_Frc_surf', frame_prs_frc_surf, nint(rv(1)), 0, 1)
     case ('cav_start_t');                   call setr('Cav_Start_T', Cav_Start_T, rv(1), 0.0_dp, 1.0e30_dp)
     case ('cav_thermal_t');                 call setr('Cav_thermal_t', Cav_thermal_t, rv(1), 0.0_dp, 1.0e30_dp)
     case ('log_intvl');          call seti('Log_Intvl', log_interval, nint(rv(1)), 1, 100000000)
     case ('data_intvl', 'data_interval');   call seti('Data_Intvl', data_interval, nint(rv(1)), 1, 100000000)
     case ('restart_intvl', 'restart_interval'); call seti('Restart_Intvl', restart_interval, nint(rv(1)), 0, 100000000)
     case ('cnvt_scheme', 'conv_scheme');    call seti('Cnvt_scheme', Cnvt_scheme, nint(rv(1)), 1, 4)
     case ('cnvt_mode');          call seti('Cnvt_mode', Cnvt_mode, nint(rv(1)), 1, 3)
     case ('dffs_mode');          call seti('Dffs_mode', Dffs_mode, nint(rv(1)), 1, 3)
     case ('cfl_switch');         call setr('cfl_switch', cfl_switch, rv(1), 0.0_dp, 1.0e30_dp)
     case ('vn_switch');          call setr('vn_switch', vn_switch, rv(1), 0.0_dp, 1.0e30_dp)
     case ('vof_cfl_max');        call setr('vof_cfl_max', vof_cfl_max, rv(1), 0.0_dp, 1.0_dp)
     case ('t_ramp');             call setr('t_ramp', t_ramp, rv(1), -1.0_dp, 1.0e30_dp)
     case ('stop_rsd_max');       call setr('Stop_Rsd_Max', Stop_Rsd_Max, rv(1), 0.0_dp, 1.0e300_dp)
     case ('stop_cfl_max');       call setr('Stop_CFL_Max', Stop_CFL_Max, rv(1), 0.0_dp, 1.0e300_dp)
     case ('max_itr_turb');       call seti('Max_Itr_Turb', Max_Itr_Turb, nint(rv(1)), 1, 1000000)
     case ('min_itr_turb');       call seti('Min_Itr_Turb', Min_Itr_Turb, nint(rv(1)), 0, 1000000)
     case ('rsd_lmt_turb');       call setr('Rsd_Lmt_Turb', Rsd_Lmt_Turb, rv(1), 0.0_dp, 1.0e30_dp)
     case ('sor_fct_turb');       call setr('SOR_fct_Turb', SOR_fct_Turb, rv(1), 0.0_dp, 1.0e30_dp)
     case ('cnvt_scheme_turb');   call seti('Cnvt_scheme_Turb', Cnvt_scheme_Turb, nint(rv(1)), 1, 2)
     case ('cnvt_mode_turb');     call seti('Cnvt_mode_Turb', Cnvt_mode_Turb, nint(rv(1)), 1, 2)
     case ('dffs_mode_turb');     call seti('Dffs_mode_Turb', Dffs_mode_Turb, nint(rv(1)), 1, 2)
     case ('cfl_switch_turb');    call setr('cfl_switch_Turb', cfl_switch_Turb, rv(1), 0.0_dp, 1.0e30_dp)
     case ('vn_switch_turb');     call setr('vn_switch_Turb', vn_switch_Turb, rv(1), 0.0_dp, 1.0e30_dp)
     case ('max_itr_engy');       call seti('Max_Itr_Engy', Max_Itr_Engy, nint(rv(1)), 1, 1000000)
     case ('min_itr_engy');       call seti('Min_Itr_Engy', Min_Itr_Engy, nint(rv(1)), 0, 1000000)
     case ('rsd_lmt_engy');       call setr('Rsd_Lmt_Engy', Rsd_Lmt_Engy, rv(1), 0.0_dp, 1.0e30_dp)
     case ('sor_fct_engy');       call setr('SOR_fct_Engy', SOR_fct_Engy, rv(1), 0.0_dp, 1.0e30_dp)
     case ('cnvt_scheme_engy');   call seti('Cnvt_scheme_Engy', Cnvt_scheme_Engy, nint(rv(1)), 1, 4)
     case ('cnvt_mode_engy');     call seti('Cnvt_mode_Engy', Cnvt_mode_Engy, nint(rv(1)), 1, 3)
     case ('dffs_mode_engy');     call seti('Dffs_mode_Engy', Dffs_mode_Engy, nint(rv(1)), 1, 3)
     case ('cfl_switch_engy');    call setr('cfl_switch_Engy', cfl_switch_Engy, rv(1), 0.0_dp, 1.0e30_dp)
     case ('vn_switch_engy');     call setr('vn_switch_Engy', vn_switch_Engy, rv(1), 0.0_dp, 1.0e30_dp)
     case ('pr_t');               call setr('Pr_t', Pr_t, rv(1), 0.0_dp, 1.0e30_dp)
     case ('poi_verbose')
        if ((nint(rv(1)) /= 0) .neqv. poi_verbose) then
           poi_verbose = (nint(rv(1)) /= 0)
           if (nmsg < 40) then
              nmsg = nmsg + 1;  write(msg(nmsg),'(A,L1)') 'poi_verbose -> ', poi_verbose
           end if
        end if
     case ('motion_start_t');     call setr('Motion_Start_T', Motion_Start_T, rv(1), -1.0e30_dp, 1.0e30_dp)
     case ('wave_ramp_t', 'relax_in_l', 'relax_out_l', 'relax_side_l')   ! incident waves (update 18)
        call wave_runtime(trim(name), rv(1), wchg)
        if (wchg .and. nmsg < 40) then
           nmsg = nmsg + 1;  write(msg(nmsg),'(A,ES12.4)') trim(name)//' -> ', rv(1)
        end if
     case default
        ! not changeable at run time (or a static input): ignored
     end select
  end do
  mom_miniter = max(1, min(mom_miniter, mom_maxiter))
  poi_miniter = max(1, min(poi_miniter, poi_maxiter))
  vof_miniter = max(1, min(vof_miniter, vof_maxiter))
  ! remaining steps from the current time with the (new) Del_T / End_T
  nstep_new = (istep - 1) + max(1, nint((End_T - t_now)/dt_dim))   ! t_now = time before this step
  if (nstep_new /= nstep) then
     if (nmsg < 40) then
        nmsg = nmsg + 1;  write(msg(nmsg),'(A,I0,A,I0)') 'time steps (from End_T and Del_T) ', nstep, ' -> ', nstep_new
     end if
     nstep = nstep_new
  end if
  if (myrank == 0) then
     write(*,'(A)') ' '
     write(*,'(A,I0,A)') ' [ name list ] probCFD.nml changed - re-read at step ', istep, ':'
     if (nmsg == 0) then
        write(*,'(A)') '     (no runtime-changeable value differs)'
     else
        do i = 1, nmsg
           write(*,'(A)') '     '//trim(msg(i))
        end do
     end if
     flush(6)
  end if
contains
  subroutine note(nm, vold, vnew)
    character(len=*), intent(in) :: nm
    real(dp), intent(in) :: vold, vnew
    if (nmsg < 40) then
       nmsg = nmsg + 1;  write(msg(nmsg),'(A,ES12.4,A,ES12.4)') nm//' ', vold, ' -> ', vnew
    end if
  end subroutine note
  subroutine seti(nm, var, vnew, lo, hi)
    character(len=*), intent(in) :: nm
    integer, intent(inout) :: var
    integer, intent(in) :: vnew, lo, hi
    if (vnew == var) return
    if (vnew < lo .or. vnew > hi) then
       if (nmsg < 40) then
          nmsg = nmsg + 1;  write(msg(nmsg),'(A,I0,A)') nm//' = ', vnew, ' ignored (out of range)'
       end if
       return
    end if
    if (nmsg < 40) then
       nmsg = nmsg + 1;  write(msg(nmsg),'(A,I0,A,I0)') nm//' ', var, ' -> ', vnew
    end if
    var = vnew
  end subroutine seti
  subroutine setr(nm, var, vnew, lo, hi)
    character(len=*), intent(in) :: nm
    real(dp), intent(inout) :: var
    real(dp), intent(in) :: vnew, lo, hi
    if (vnew == var) return
    if (vnew < lo .or. vnew > hi) then
       if (nmsg < 40) then
          nmsg = nmsg + 1;  write(msg(nmsg),'(A,ES12.4,A)') nm//' = ', vnew, ' ignored (out of range)'
       end if
       return
    end if
    call note(nm, var, vnew)
    var = vnew
  end subroutine setr
  function lower2(s) result(t)
    character(len=*), intent(in) :: s
    character(len=len(s)) :: t
    integer :: i, c
    t = s
    do i = 1, len(s)
       c = iachar(s(i:i))
       if (c >= iachar('A') .and. c <= iachar('Z')) t(i:i) = achar(c + 32)
    end do
  end function lower2
  subroutine split2(s, out, n)
    character(len=*), intent(in)  :: s
    character(len=*), intent(out) :: out(:)
    integer,          intent(out) :: n
    integer :: i, L
    character(len=1) :: ch
    n = 0;  out = ' '
    L = len_trim(s)
    i = 1
    do while (i <= L)
       ch = s(i:i)
       if (ch == ' ' .or. ch == ',' .or. ch == char(9)) then
          i = i + 1;  cycle
       end if
       n = n + 1
       if (n > size(out)) exit
       do while (i <= L)
          if (s(i:i) == ' ' .or. s(i:i) == ',' .or. s(i:i) == char(9)) exit
          out(n) = trim(out(n))//s(i:i)
          i = i + 1
       end do
    end do
  end subroutine split2
end subroutine reread_namelist
