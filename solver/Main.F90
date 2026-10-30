!=======================================================================
! Main.F90
!   probCFD_부산 : 3-D incompressible Navier-Stokes solver (MULTI-BLOCK,
!                  CURVILINEAR, single- or two-phase)
!
!   - Fortran 90, finite difference, fully-staggered structured grid
!     with a general 3-D curvilinear coordinate transformation
!   - non-dimensionalized by Char_L, Char_U, tch = L/U, Lou_w (fluid 1)
!   - projection (fractional time step) method; per time step the
!     sequence below is repeated Iner_Itn times (inner iterations):
!       (0) two-phase: advect the liquid volume fraction alpha from
!           alpha^n with the latest divergence-free face fluxes
!           (CICSAM, explicit or implicit), rho and mu follow from it
!       (1) predictor : u* from  u*/dt + CONV_imp(u*) - DIFF_imp(u*)
!                           = u^n/dt - CONV_exp(u^m) + DIFF_exp(u^m) + g
!                           [ - grad(p^m)/rho with Prs_form = 2: incremental
!                             form, the projection then applies the pressure
!                             increment; update 16 ]
!           (Cnvt_mode / Dffs_mode choose explicit, implicit or
!            per-direction hybrid treatment; the advecting fluxes and
!            the explicit terms use the latest iterate u^m, m = the
!            inner iteration; u^m = u^n for m = 1)
!       (2) Poisson   : div( (1/rho) grad p ) = div(u*)/dt
!                       (Jacobi or FGMRES+PCG, poi_solver)
!       (3) corrector : u^(m+1) = u* - dt (1/rho) grad p
!     With Iner_Itn = 1 this is the classic non-iterative fractional
!     step; Iner_Itn > 1 is a Picard iteration towards the fully
!     implicit (backward Euler) coupled solution of the time step.
!   - laminar, multi-block structured grid, block-to-block (patched)
!     interfaces detected at start-up (b2b_tol) with full support for
!     differing block orientations and partial-face windows
!   - MPI (OpenMPI): per-block sub-domain splitting from the csv file
!
!   Input files:
!     probCFD.nml                  grouped name list (annotated example)
!     Input_polyMesh_probCFD.csv   per-block i/j/k splitting
!     <Case_name>.grd              PLOT3D grid (ASCII, multi-block)
!     <Case_name>.bcs              boundary conditions (per block)
!     <Case_name>.bcs_b2b          written + read back at start-up
!
!   Output files (binary PLOT3D, multi-block):
!     Volume_solutions_grids/<Case_name>_NNNNNNNN.grd   grid, node coordinates [m]
!     Volume_solutions_grids/<Case_name>_NNNNNNNN.p3d   u,v,w [m/s], p [Pa] (+ alpha, p_dyn, vapour, temperature) at nodes
!     <Case_name>_monitor.dat      one line per time step (history)
!     Reason_for_stop.md           why the run ended (stop criteria / input error / normal end)
!     <Case_name>_NNNNNNNN.rsto    restart file every Restart_Intvl steps (Mode = 1 restarts from the latest)
!     Output_data/<Case_name>_*.dat 6-DoF forces / moments / motions of the bodies, child and
!                                  sub-child bodies (sixdof.F90, update 14)
!     Output_data/<Case_name>_heat.dat  energy equation: T range, heat flow of every heat surface (update 16)
!     Output_data/<Case_name>_wave.dat  incident waves: probe elevations (computed / exact), _components.wav (update 18)
!     Output_data/<Case_name>_force_model<n>.dat  force models: wrench on the mother, disk n, J, K_T, K_Q (update 19)
!     Check_grid_movement/<Case_name>_NNNNNNNN.grd  Mode = 2: the moved grids of the imposed motions (update 14b)
!
!   Screen log per time step (Log_Intvl): current time, maximum local
!   CFL numbers (total and per i/j/k direction), min/max/avg cell
!   |divergence|, the u/v/w face monitors (sweeps, max/avg residual),
!   the pressure monitor per inner iteration (iterations, residual,
!   update norm), the VoF monitor and the wall-clock time of the step.
!=======================================================================
program probCFD_main
  use global_vars
  use mpi_vars
  use curvi, only: fill_ghost_coords, build_metrics, cart_to_contra,      &
                   contra_to_cart
  use b2b_runtime, only: b2b_setup, refresh_all, refresh_uvw
  use mom_ops, only: sync_face_jacobian
  use vof, only: init_alpha, advect_alpha, liquid_volume, alpha_range
  use turbulence, only: init_turbulence, solve_turbulence, turb_monitor, write_wallfn_data
  use stop_criteria
  use surface_tension, only: compute_curvature, st_monitor
  use overset, only: ovst_setup, ovs_update_faces, ovs_finish_status, ovs_selftest
  use sixdof, only: sixdof_precheck, sixdof_check_movement, sixdof_init, sixdof_advance, sixdof_forces, sixdof_export, sixdof_monitor, &
                    sixdof_finalize, sixdof_summary
  use restart, only: write_restart, read_restart, restart_available
  use cavitation, only: init_cavitation, cav_init_field, cav_sources, advect_vapor, cav_pressure_update, vapor_volume, &
                        cav_monitor, cav_div_target
  use energy, only: init_energy, energy_start, solve_energy, energy_monitor, energy_export, temp_range, energy_finalize, &
                    energy_iters, energy_resid
  use waves, only: wave_on, wave_setup, wave_start, wave_update, wave_export, wave_monitor, wave_finalize
  use force_model, only: fm_on, force_init, force_compute, force_monitor, force_finalize
  use pid_control, only: pid_hold
  use probes, only: probe_init, probe_export, probe_finalize
  use surface_output, only: surf_write
  use statistics, only: stat_init, stat_sample, stat_write, stat_finalize, stat_summary
  use morph_frames, only: frames_hooks
  use patch_flux, only: patch_flux_init, patch_flux_export, patch_flux_finalize
  use valves, only: valve_init, valve_update, valve_export, valve_finalize
  use fsi, only: fsi_init, fsi_predict, fsi_iterate, fsi_finish_step, fsi_monitor, fsi_summary, fsi_finalize, fsi_max_iterations
  implicit none

  integer  :: istep, ierr, ib, b, m, imon, bb, stage, nstage, istep0, nit
  real(dp) :: rvliq0, rclip
  integer  :: itm, itp, nsub, nvof
  real(dp) :: resm(3), resa(3), resp, updp, t, vliq, vliq0, amin, amax, clip, cliptot
  real(dp) :: qin, qout, qscale, resv, resva, vvap, vmaxv, tmn, tmx, tav
  real(dp) :: cfl_tot, cfl_dir(3), dmin, dmax, davg, twall0, twall1
  real(dp) :: tw_vof, tw_mom, tw_poi, tw_a
  logical  :: dolog, fsi_done

  !=================== initialization ==================================
  call mpi_setup
  call print_banner
  call read_namelist
  call make_bcs_b2b_file    ! step 1: detect block-to-block interfaces
                            ! and export <casename>.bcs_b2b
  call read_grid            ! multi-block grid (node coordinates)
  call read_split_csv       ! per-block i/j/k MPI splitting: N columns,
                            ! CPUs per block in row 5, total = sum(row 5)
  call sixdof_precheck      ! 6-DoF parts (rows 6-8 + $SixDOF): data, imposed files (before the assembly)
  call frames_hooks         ! prescribed grid deformation from frame files (heart extension H1): headers, initial grid
  if (run_mode == 2) call sixdof_check_movement   ! Mode = 2: grid-movement check of the imposed motions, then stop
  call mpi_assign_domains   ! build the sub-domains and allocate the CPUs
  call read_bcs             ! physical BCs, per block (+ periodic)
  call wave_setup           ! incident waves (update 18): components, wave patches, checks
  call b2b_setup            ! step 2: read <casename>.bcs_b2b, build the
                            ! exchange maps / masks
  call fill_ghost_coords    ! ghost node coordinates (extrapolation,
                            ! donor blocks, periodic translation)
  if (ovst_only /= 0) then  ! GRID CHECK MODE: overset assembly only,
     call ovst_setup        ! without the field arrays (large grids fit)
     call write_grid_plot3d(0)
     call write_iblank_plot3d
     if (myrank == 0) write(*,'(A)') ' Ovst_only = 1: overset assembly written ('//trim(casename)// &
        '.ovst report, Volume_solutions_grids/'//trim(casename)//'_00000000.grd grid, _00000000_iblank.p3d), stopping.'
     call MPI_Finalize(ierr)
     stop
  end if
  call allocate_fields
  call init_cavitation      ! cavitation (update 15): vapour arrays (the mixture properties need them
                            ! from here on), pressure level, sigma
  call init_energy          ! energy equation (update 16): temperature arrays, T = T_init
  call build_metrics        ! half-index curvilinear metrics (m, J)
  call sync_face_jacobian   ! shared interface faces: owner's J_f
  call ovst_setup           ! overset (Chimera) assembly: hole cutting,
                            ! active / fringe / hole cells, donor search
  bc_ramp = 1.0_dp
  if (t_ramp > 0.0_dp) bc_ramp = 0.0_dp   ! ramped BCs start from zero
  call initialize_fields
  if (ovst_on .and. poi_solver == 2) then   ! inner GMRES basis of the coupled overset pressure operator
     do ib = 1, nown
        bb = own_list(ib)
        allocate(dom(bb)%kW(dom(bb)%r1(1):dom(bb)%r2(1), dom(bb)%r1(2):dom(bb)%r2(2), &
                            dom(bb)%r1(3):dom(bb)%r2(3), ovs_mkry+1))
        dom(bb)%kW = 0.0_dp
     end do
  end if
  call init_turbulence          ! turbulence models / wall functions (update 9)
  call energy_start             ! energy equation (update 16): boundary-face records, heat-surface file, T ghosts
  call cav_init_field           ! cavitation (update 17): initial vapour fraction vap_cel_initial of the liquid cells
  call fsi_init                 ! flexible parts (update 24): mid-surfaces, shells, morph data (grids at t = 0)
  call sixdof_init              ! 6-DoF bodies (update 12): object ids, mass data, imposed motions
  call wave_start               ! incident waves (update 18): relaxation-zone box, probe / component files
  call wave_init_fields         !   Wave_Init = 1: alpha and velocity from the exact wave at t = 0
  call force_init               ! force models (update 19): Hough-Ordway / BEM disks, point forces
  call force_compute(0, 0.0_dp) !   sources and wrenches of the first step from the initial field
  call probe_init               ! volume-field probe points (update 22)
  call patch_flux_init          ! volume flux through the inlet / outlet patches (heart extension H1)
  call valve_init               ! heart valves as porous zones (H2)
  call stat_init                ! turbulence statistics (update 23)
  ! --- Mode = 1: restart from the latest <case>_NNNNNNNN.rsto (update 13b)
  istep0 = 0;  t_now = 0.0_dp;  rvliq0 = -1.0_dp;  rclip = 0.0_dp
  if (run_mode == 1) then
     istep0 = restart_available()
     call MPI_Bcast(istep0, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     if (istep0 <= 0) then
        if (myrank == 0) write(*,'(A)') '  Mode = 1 but no restart file is listed in '//trim(casename)// &
           '_restart.log - starting from step 0'
        run_mode = 0
     else
        call read_restart(istep0, t_now, rvliq0, rclip)
        ! force models of the restored state (update 21): momentum sources and wrenches of the
        ! first step from the restored field, poses and controller outputs (the controllers hold)
        pid_hold = .true.
        call force_compute(istep0, t_now)
        pid_hold = .false.
     end if
  end if
  if (run_mode == 0) then
     call sixdof_forces            ! forces of the initial field
     call sixdof_export(0, 0.0_dp)
     call energy_export(0, 0.0_dp)
     call probe_export(0, 0.0_dp)
     call patch_flux_export(0, 0.0_dp)
     call valve_export(0, 0.0_dp)
     call write_grid_plot3d(0)
     call write_solution(0)
     call surf_write(0)
  else
     if (any_moving) call write_grid_plot3d(istep0)   ! the restored grid pose
     call write_solution(istep0)   ! the restored field (also refreshes the .nam name file)
  end if
  if (poi_verbose) then
     call ovs_selftest                ! interpolation exchange check
     call initialize_fields           ! (the test overwrote u, v, w)
  end if

  cliptot = 0.0_dp
  vliq = 0.0_dp;  amin = 0.0_dp;  amax = 0.0_dp
  if (Num_of_Phase == 2) then
     vliq0 = liquid_volume() * Lch**3
     if (myrank == 0) write(*,'(A,ES14.6,A)') '  initial liquid volume : ', &
        vliq0, ' m3'
  end if
  if (run_mode == 1) then
     if (rvliq0 >= 0.0_dp) vliq0 = rvliq0       ! monitor sums of the restarted run
     cliptot = rclip
  end if
  if (myrank == 0 .and. run_mode == 1) then
     open(newunit=imon, file=trim(casename)//'_monitor.dat', status='unknown', position='append', action='write')
     write(imon,'(A,I0,A,ES13.5,A)') '# restart at step ', istep0, ' (t = ', t_now, ' s)'
     write(*,'(A)') ' '
     write(*,'(A,I0,A,I0,A)') ' ============================== time marching (restart: steps ', istep0+1, ' .. ', nstep, ') =========='
  else if (myrank == 0) then
     open(newunit=imon, file=trim(casename)//'_monitor.dat', status='replace', action='write')
     if (Num_of_Phase == 2) then
        write(imon,'(A)', advance='no') '#   step      time[s]   CFL_tot   CFL_i   CFL_j   CFL_k'// &
           '   |div|min   |div|max   |div|avg  itm  res(u*)  itp  res(p)  upd(p)'// &
           '  itVoF res(VoF) nsub    Vliq[m3]    amin    amax  clipped[m3]'
     else
        write(imon,'(A)', advance='no') '#   step      time[s]   CFL_tot   CFL_i   CFL_j   CFL_k'// &
           '   |div|min   |div|max   |div|avg  itm  res(u*)  itp  res(p)  upd(p)'
     end if
     if (cav_on) write(imon,'(A)', advance='no') '     Vvap[m3]  avmax'
     if (energy_on) write(imon,'(A)', advance='no') '  itT  res(T)     Tmin[K]     Tmax[K]'
     write(imon,'(A)') ' '
     write(*,'(A)') ' '
     write(*,'(A)') ' ============================== time marching =============================='
  end if

  !=================== time marching loop ==============================
  istep = istep0
  if (myrank == 0 .and. istep0 >= nstep) write(*,'(A,I0,A,I0,A)') '  NOTE: the restart step ', istep0, &
     ' is not before the last step ', nstep, ' - raise End_T to continue'
  do while (istep < nstep)                    ! nstep may change at run time (name list re-read)
     istep = istep + 1
     twall0 = MPI_Wtime()
     call reread_namelist(istep)              ! runtime name-list updates (update 13)
     t = t_now + dt_dim;  t_now = t           ! accumulated (Del_T may change at run time)
     dolog = (myrank == 0 .and. (mod(istep, log_interval) == 0 .or. istep == 1))
     ! --- 6-DoF: body motion over this step (forces of the previous step),
     !     grid update, overset re-assembly, grid fluxes (the boundary
     !     values of the time-level-n field keep the ramp of the previous
     !     step: a body at rest reproduces the static solver exactly) -----
     call fsi_predict(istep, t)              ! flexible parts (update 24): predicted interface displacement
     call sixdof_advance(istep, t)
     call valve_update(t)                    ! valve closure at the target time (H2)

     ! start-up ramp of the Dirichlet velocity BCs (inlets, moving walls)
     bc_ramp = 1.0_dp
     if (t_ramp > 0.0_dp) bc_ramp = min(t/t_ramp, 1.0_dp)
     if (wave_on) call wave_update(t)         ! incident wave: BC time t^(n+1), ramp, frame shift (update 18)

     ! --- time level n copies (inner iterations restart from them) ------
     do ib = 1, nown
        b = own_list(ib)
        dom(b)%un = dom(b)%u;  dom(b)%vn = dom(b)%v;  dom(b)%wn = dom(b)%w
        if (Num_of_Phase == 2) dom(b)%alpn = dom(b)%alpha
        if (cav_on) then
           dom(b)%vapn = dom(b)%vap;  dom(b)%presn = dom(b)%p
        end if
        if (energy_on) dom(b)%tempn = dom(b)%temp
     end do

     ! --- CFL / divergence of the field the step starts from (u^n) ------
     call cfl_numbers(cfl_tot, cfl_dir)
     call check_cfl(cfl_tot, istep)                 ! stop criterion 3
     call divergence_stats(dmin, dmax, davg)
     if (dolog) then
        do ib = 1, 21                         ! 21 empty lines between step blocks (3 x the former 7; update 13)
           write(*,'(A)') ' '
        end do
        write(*,'(A,I8,A,F24.10,A)') ' Step ', istep, '   Current Time(real) : ', t, ' s'
        write(*,'(A,F37.4)')   ' Max Local Total CFL#        : ', cfl_tot
        write(*,'(A,F37.4)')   ' Max Local <i> CFL #         : ', cfl_dir(1)
        write(*,'(A,F37.4)')   ' Max Local <j> CFL #         : ', cfl_dir(2)
        write(*,'(A,F37.4)')   ' Max Local <k> CFL #         : ', cfl_dir(3)
        write(*,'(A,3ES11.2)') ' Cell abs(DIV) {min/max/avg} : ', dmin, dmax, davg
     end if

     nsub = 0;  clip = 0.0_dp;  nvof = 0;  resv = 0.0_dp;  resva = 0.0_dp
     tw_vof = 0.0_dp;  tw_mom = 0.0_dp;  tw_poi = 0.0_dp
     ! --- Time_scheme = 2: 2nd-order SSP (Heun) - two forward stages of
     !     the fractional step, U^(1) = S(U^n), U^(2) = S(U^(1)), then
     !     U^(n+1) = (U^n + U^(2))/2 for every transported field (velocity
     !     faces and contravariant fluxes, volume fraction, turbulence); the
     !     average of two discretely divergence-free fields is divergence-
     !     free, the pressure is that of the second stage (update 13) -----
     nstage = merge(2, 1, Time_scheme == 2)
     if (nstage == 2) call ssp_save
     do stage = 1, nstage
     if (stage == 2) then
        do ib = 1, nown                          ! stage 2 starts from stage 1
           b = own_list(ib)
           dom(b)%un = dom(b)%u;  dom(b)%vn = dom(b)%v;  dom(b)%wn = dom(b)%w
           if (Num_of_Phase == 2) dom(b)%alpn = dom(b)%alpha
           if (cav_on) then
              dom(b)%vapn = dom(b)%vap;  dom(b)%presn = dom(b)%p
           end if
           if (energy_on) dom(b)%tempn = dom(b)%temp
        end do
     end if
     nit = Iner_Itn
     if (fsi_on) nit = max(Iner_Itn, fsi_max_iterations())    ! FSI sub-iterations (update 24)
     do m = 1, nit
        ! --- step 0: volume-fraction transport (two-phase) ---------------
        tw_a = MPI_Wtime()
        if (Num_of_Phase == 2) then
           if (m > 1) then
              do ib = 1, nown
                 b = own_list(ib)
                 dom(b)%alpha = dom(b)%alpn
              end do
           end if
           call advect_alpha(nsub, clip, nvof, resv, resva)
        end if
        ! --- cavitation: mass transfer from the latest pressure, vapour transport (update 15)
        if (cav_on) then
           if (m > 1) then
              do ib = 1, nown
                 b = own_list(ib)
                 dom(b)%vap = dom(b)%vapn
              end do
           end if
           call cav_sources(istep, t)
           call advect_vapor
        end if
        tw_vof = tw_vof + MPI_Wtime() - tw_a

        ! --- explicit RHS (u^n/dt, explicit terms of the latest iterate) --
        tw_a = MPI_Wtime()
        call refresh_uvw(1)                ! u, v, w combined (b2b safe)
        if (st_on) call compute_curvature          ! surface tension: CSF curvature of the interface
        call compute_rhs_momentum
        call add_pgrad(1)                  ! Prs_form = 2: - grad p^n / rho in the predictor (update 16)

        ! --- step 1: predictor (intermediate velocities) -----------------
        do ib = 1, nown
           b = own_list(ib)                    ! b = sub-domain id
           dom(b)%us = dom(b)%u;  dom(b)%vs = dom(b)%v;  dom(b)%ws = dom(b)%w
        end do
        call solve_momentum(itm, resm, resa)
        call add_pgrad(2)                  ! Prs_form = 2: u* + dt grad p^n / rho -> the projection corrects with the increment
        call check_residual('momentum', maxval(resm), 0.0_dp, istep, 1)   ! stop criterion 2
        call check_nan('after momentum', istep)                          ! stop criterion 1   ! u, v, w together (one combined
                                               ! Jacobi loop: required for b2b
                                               ! interfaces with permuted axes)
        call outflow_correction(qin, qout, qscale)   ! mass-corrected outflow faces
        tw_mom = tw_mom + MPI_Wtime() - tw_a

        ! --- step 2: pressure-Poisson equation ---------------------------
        tw_a = MPI_Wtime()
        call solve_poisson(itp, resp, updp)
        call check_residual('pressure', resp, updp, istep, 2)            ! stop criterion 2
        call check_nan('after the pressure solve', istep)                ! stop criterion 1

        ! --- step 3: velocity correction ----------------------------------
        call correct_velocity
        if (quasi2d > 0) call quasi2d_project(istep)      ! SW_Quasi_2D (update 24c): w = 0, fields averaged across the span
        if (cav_on) call cav_pressure_update    ! vapour fraction consistent with the pressure-implicit source
        call solve_turbulence           ! turbulence transport + eddy viscosity (update 9)
        ! --- energy equation (update 16): T from T^n with the corrected fluxes, the latest
        !     eddy viscosity and the cavitation mass transfer (its buoyancy enters the next
        !     momentum right-hand side)
        if (energy_on) then
           if (m > 1) then
              do ib = 1, nown
                 b = own_list(ib)
                 dom(b)%temp = dom(b)%tempn
              end do
           end if
           call solve_energy
        end if
        call check_nan('after the projection / VoF / turbulence / energy', istep)   ! stop criterion 1
        tw_poi = tw_poi + MPI_Wtime() - tw_a

        ! --- inner-iteration monitors --------------------------------------
        if (dolog) then
           write(*,'(A)') ' '
           if (Iner_Itn > 1) write(*,'(A,I4,A,I4,A)') ' ----- Inner Itr# ', m, ' of ', Iner_Itn, ' -----'
           write(*,'(A)')            ' [ I ] Face Monitor :          < u Vel >'
           write(*,'(A,I15)')        '  --> TTL # of Itr  : ', itm
           write(*,'(A,ES19.4)')     '  --> Max Residual  : ', resm(1)
           write(*,'(A,ES19.4)')     '  --> Avg Residual  : ', resa(1)
           write(*,'(A)') ' '
           write(*,'(A)')            ' [ J ] Face Monitor :          < v Vel >'
           write(*,'(A,I15)')        '  --> TTL # of Itr  : ', itm
           write(*,'(A,ES19.4)')     '  --> Max Residual  : ', resm(2)
           write(*,'(A,ES19.4)')     '  --> Avg Residual  : ', resa(2)
           write(*,'(A)') ' '
           write(*,'(A)')            ' [ K ] Face Monitor :          < w Vel >'
           write(*,'(A,I15)')        '  --> TTL # of Itr  : ', itm
           write(*,'(A,ES19.4)')     '  --> Max Residual  : ', resm(3)
           write(*,'(A,ES19.4)')     '  --> Avg Residual  : ', resa(3)
           write(*,'(A)') ' '
           write(*,'(A,I12)')        ' [ Cell ] Center Monitor : Pressure - Inner Itr#: ', m
           write(*,'(A,I26)')        '  ---> Total Number of Iteration      ', itp
           write(*,'(A,ES26.6)')     '  ---> Pressure Residual |b-Ap|/|b|   ', resp
           write(*,'(A,ES26.6)')     '  ---> Pressure Update norm max|dp|   ', updp
           if (Num_of_Phase == 2) then
              write(*,'(A)') ' '
              select case (vof_scheme)
              case (1)
                 write(*,'(A)')      ' [ Cell ] Center Monitor :    < CICSAM-VoF (explicit) >'
                 write(*,'(A,I15)')  '  --> sub-cycles    : ', nsub
              case (2)
                 write(*,'(A)')      ' [ Cell ] Center Monitor :    < CICSAM-VoF (implicit) >'
                 write(*,'(A,I15)')  '  --> TTL # of Itr  : ', nvof
                 write(*,'(A,ES19.4)')'  --> Max Residual  : ', resv
                 write(*,'(A,ES19.4)')'  --> Avg Residual  : ', resva
                 write(*,'(A,I15)')  '  --> sub-cycles    : ', nsub
              case default
                 write(*,'(A)')      ' [ Cell ] Center Monitor :    < Geometric PLIC-VoF (direction-split, Weymouth-Yue) >'
                 write(*,'(A,I15)')  '  --> sub-cycles    : ', nsub
                 write(*,'(A,F15.3)')'  --> face CFL limit: ', min(vof_cfl_max, 0.5_dp)
              end select
           end if
           call turb_monitor
           call st_monitor
           call cav_monitor
           call energy_monitor
        end if
        ! --- fluid-structure interaction (update 24): wall forces -> structure -> relaxation ->
        !     grid morphing; the next inner iteration solves the flow on the moved grid from the
        !     level-n field; the loop ends when every flexible part has converged (and the flow's
        !     own Iner_Itn iterations are done)
        if (fsi_on) then
           call fsi_iterate(istep, t, m, fsi_done, m == nit)
           if (fsi_done .and. m >= Iner_Itn) exit
        end if
     end do
     end do                                      ! SSP stages
     if (nstage == 2) call ssp_average
     if (wave_on) call wave_relax_fields(t)         ! incident waves: relaxation zones (update 18)

     ! --- step monitoring: the corrected field u^(n+1) --------------------
     call sixdof_forces                             ! fluid forces on the bodies
     call force_compute(istep, t)                   ! force models: sources of the next step, wrench on the mothers (update 19)
     call fsi_finish_step(istep, t)                 ! flexible parts (update 24): commit the structural states, monitor files
     call sixdof_export(istep, t)
     call energy_export(istep, t)                   ! heat flows of the heat surfaces (update 16)
     if (wave_on) call wave_export(istep, t)        ! wave probes (update 18)
     call probe_export(istep, t)                    ! volume-field probe points (update 22)
     call patch_flux_export(istep, t)               ! inlet / outlet volume fluxes (H1)
     call valve_export(istep, t)
     call stat_sample(istep, t)                     ! turbulence statistics (update 23)
     call stat_write(istep, .false.)
     call cfl_numbers(cfl_tot, cfl_dir)
     call divergence_stats(dmin, dmax, davg)
     if (Num_of_Phase == 2) then
        vliq = liquid_volume() * Lch**3
        call alpha_range(amin, amax)
        cliptot = cliptot + clip * Lch**3
     end if
     vvap = 0.0_dp;  vmaxv = 0.0_dp
     if (cav_on) then
        vvap = vapor_volume() * Lch**3
        vmaxv = vap_max()
     end if
     tmn = 0.0_dp;  tmx = 0.0_dp;  tav = 0.0_dp
     if (energy_on) call temp_range(tmn, tmx, tav)
     twall1 = MPI_Wtime()
     if (dolog) then
        write(*,'(A)') ' '
        write(*,'(A,3ES11.2,A)') ' Cell abs(DIV) {min/max/avg} : ', dmin, dmax, davg, '  (after correction)'
        if (Num_of_Phase == 2) write(*,'(A,ES13.5,A,2F9.5,A,ES11.3,A)')          &
           ' Liquid volume [m3] : ', vliq, '   alpha min/max :', amin, amax,      &
           '   clipped [m3] :', cliptot, ' (cumulative)'
        if (cav_on) write(*,'(A,ES13.5,A,F9.5)') ' Vapour volume [m3] : ', vvap, '   max vapour fraction :', vmaxv
        if (energy_on) write(*,'(A,3F11.4)') ' Temperature [K] min / max / avg : ', tmn, tmx, tav
        write(*,'(A)') ' '
        write(*,'(A,F16.3,A)') ' ClockTime[s] for ONE time step approach: ', twall1 - twall0, ' s'
        write(*,'(A,3F10.3,A)') '   (VoF / momentum / Poisson+correction : ', tw_vof, tw_mom, tw_poi, ' s)'
        call sixdof_monitor
        call fsi_monitor
        call wave_monitor
        call force_monitor
     end if
     if (myrank == 0) call record_step(istep, t, cfl_tot, itm, maxval(resm), itp, resp, updp, vliq, twall1 - twall0)
     if (myrank == 0) then
        if (Num_of_Phase == 2) then
           write(imon,'(I7, ES13.4, 4F8.4, 3ES11.3, I5, ES9.2, I6, 2ES9.2, I6, ES9.2, I5, ES13.5, 2F8.4, ES11.3)', advance='no') &
              istep, t, cfl_tot, cfl_dir, dmin, dmax, davg, itm, maxval(resm), itp, resp, updp, &
              nvof, resv, nsub, vliq, amin, amax, cliptot
        else
           write(imon,'(I7, ES13.4, 4F8.4, 3ES11.3, I5, ES9.2, I6, 2ES9.2)', advance='no') &
              istep, t, cfl_tot, cfl_dir, dmin, dmax, davg, itm, maxval(resm), itp, resp, updp
        end if
        if (cav_on) write(imon,'(ES13.5, F8.4)', advance='no') vvap, vmaxv
        if (energy_on) write(imon,'(I5, ES9.2, 2ES12.4)', advance='no') energy_iters(), energy_resid(), tmn, tmx
        write(imon,'(A)') ' '
        flush(imon)
     end if

     ! --- solution export -----------------------------------------------
     if (mod(istep, out_interval) == 0) then
        if (any_moving) call write_grid_plot3d(istep)   ! moving grids: <case>_NNNNNNNN.grd
        call write_solution(istep);  call write_wallfn_data(istep)
        call surf_write(istep)                      ! solid-wall surface files (update 22)
     end if
     ! --- restart file (Restart_Intvl) -------------------------------------
     if (restart_interval > 0) then
        if (mod(istep, restart_interval) == 0 .and. istep < nstep) call write_restart(istep, t, vliq0, cliptot)
     end if
  end do

  if (mod(nstep, out_interval) /= 0) then
     if (any_moving) call write_grid_plot3d(nstep)
     call write_solution(nstep);  call write_wallfn_data(nstep)
     call surf_write(nstep)
  end if
  call stat_write(nstep, .true.)                        ! final statistics (update 23)
  if (restart_interval > 0 .and. nstep > istep0) call write_restart(nstep, t_now, vliq0, cliptot)
  call sixdof_summary
  call fsi_summary
  call stat_summary
  call fsi_finalize
  call sixdof_finalize
  call energy_finalize
  call wave_finalize
  call force_finalize
  call probe_finalize
  call patch_flux_finalize
  call valve_finalize
  call stat_finalize
  call stop_report_normal(nstep, t_now)         ! Reason_for_stop.md: finished normally
  if (myrank == 0) then
     close(imon)
     write(*,'(A)') ' probCFD_부산 finished normally.'
  end if
  call MPI_Finalize(ierr)

contains

  !---------------------------------------------------------------------
  ! 2nd-order SSP (Time_scheme = 2): copies of the time-level-n fields
  subroutine ssp_save
    integer :: ib, b
    do ib = 1, nown
       b = own_list(ib)
       if (.not. allocated(dom(b)%s_u)) then
          allocate(dom(b)%s_u,  mold=dom(b)%u);   allocate(dom(b)%s_v,  mold=dom(b)%v);   allocate(dom(b)%s_w,  mold=dom(b)%w)
          allocate(dom(b)%s_U1, mold=dom(b)%U1);  allocate(dom(b)%s_U2, mold=dom(b)%U2);  allocate(dom(b)%s_U3, mold=dom(b)%U3)
          if (Num_of_Phase == 2) allocate(dom(b)%s_alpha, mold=dom(b)%alpha)
          if (cav_on) then
             allocate(dom(b)%s_vap, mold=dom(b)%vap);  allocate(dom(b)%s_cavs, mold=dom(b)%cavs)
          end if
          if (turb_on) then
             allocate(dom(b)%s_tk, mold=dom(b)%tk);  allocate(dom(b)%s_tw, mold=dom(b)%tw);  allocate(dom(b)%s_mut, mold=dom(b)%mut)
          end if
          if (energy_on) allocate(dom(b)%s_temp, mold=dom(b)%temp)
       end if
       dom(b)%s_u = dom(b)%u;   dom(b)%s_v = dom(b)%v;   dom(b)%s_w = dom(b)%w
       dom(b)%s_U1 = dom(b)%U1; dom(b)%s_U2 = dom(b)%U2; dom(b)%s_U3 = dom(b)%U3
       if (Num_of_Phase == 2) dom(b)%s_alpha = dom(b)%alpha
       if (cav_on) then
          dom(b)%s_vap = dom(b)%vap;  dom(b)%s_cavs = dom(b)%cavs
       end if
       if (turb_on) then
          dom(b)%s_tk = dom(b)%tk;  dom(b)%s_tw = dom(b)%tw;  dom(b)%s_mut = dom(b)%mut
       end if
       if (energy_on) dom(b)%s_temp = dom(b)%temp
    end do
  end subroutine ssp_save

  ! U^(n+1) = (U^n + U^(2))/2
  subroutine ssp_average
    integer :: ib, b
    do ib = 1, nown
       b = own_list(ib)
       dom(b)%u  = 0.5_dp*(dom(b)%s_u  + dom(b)%u);   dom(b)%v  = 0.5_dp*(dom(b)%s_v  + dom(b)%v)
       dom(b)%w  = 0.5_dp*(dom(b)%s_w  + dom(b)%w)
       dom(b)%U1 = 0.5_dp*(dom(b)%s_U1 + dom(b)%U1);  dom(b)%U2 = 0.5_dp*(dom(b)%s_U2 + dom(b)%U2)
       dom(b)%U3 = 0.5_dp*(dom(b)%s_U3 + dom(b)%U3)
       if (Num_of_Phase == 2) dom(b)%alpha = 0.5_dp*(dom(b)%s_alpha + dom(b)%alpha)
       if (cav_on) then
          dom(b)%vap  = 0.5_dp*(dom(b)%s_vap + dom(b)%vap)
          dom(b)%cavs = 0.5_dp*(dom(b)%s_cavs + dom(b)%cavs)   ! the averaged field's divergence target
       end if
       if (turb_on) then
          dom(b)%tk  = 0.5_dp*(dom(b)%s_tk  + dom(b)%tk);   dom(b)%tw  = 0.5_dp*(dom(b)%s_tw  + dom(b)%tw)
          dom(b)%mut = 0.5_dp*(dom(b)%s_mut + dom(b)%mut)
       end if
       if (energy_on) dom(b)%temp = 0.5_dp*(dom(b)%s_temp + dom(b)%temp)
    end do
    if (energy_on) call refresh_all(20)
  end subroutine ssp_average

  !---------------------------------------------------------------------
  subroutine print_banner
    if (myrank /= 0) return
    write(*,'(A)') ' ==============================================================================='
    write(*,'(A)') ' ||                                                                           ||'
    write(*,'(A)') ' ||                                                                           ||'
    write(*,'(A)') ' ||     pppp   r  r   ooo     b       CCCCC  FFFFFF  DDDDD                    ||'
    write(*,'(A)') ' ||     p   p  r r   o   o    b      C       F       D    D     probCFD_부산  ||'
    write(*,'(A)') ' ||     p  p   rr   o     o   b      C       FFFF    D     D                  ||'
    write(*,'(A)') ' ||     ppp    r    o     o   bbbb   C       F       D     D   multi-block    ||'
    write(*,'(A)') ' ||     p      r    o     o   b   b  C       F       D     D   curvilinear    ||'
    write(*,'(A)') ' ||     p      r     o   o    b   b   C      F       D    D    two-phase      ||'
    write(*,'(A)') ' ||     p      r      ooo     bbbb     CCCCC F       DDDDD     MPI            ||'
    write(*,'(A)') ' ||                                                                           ||'
    write(*,'(A)') ' ||                                                                           ||'
    write(*,'(A)') ' ==============================================================================='
    write(*,'(A)') ' '
    write(*,'(A)') ' probCFD_부산'
    write(*,'(A)') ' '
  end subroutine print_banner

  !---------------------------------------------------------------------
  ! Maximum local CFL numbers of the corrected face fluxes:
  !   CFL_q = dt |U^q| per direction (computational, d(xi) = 1);
  !   total = max over the cells of sum_q CFL_q (face averages).
  !---------------------------------------------------------------------
  subroutine cfl_numbers(ctot, cdir)
    real(dp), intent(out) :: ctot, cdir(3)
    integer  :: bb, i, j, k
    real(dp) :: c1, c2, c3
    ctot = 0.0_dp;  cdir = 0.0_dp
    do ib = 1, nown
       bb = own_list(ib)
       call set_block(bb)
       do k = ks, ke
        do j = js, je
         do i = is, ie
            if (cstat(i,j,k) /= 1) cycle                 ! overset fringe / hole
            ! (relative to the grid on moving blocks: the convective CFL of the ALE form)
            c1 = 0.5_dp*dt*( abs(U1(i,j,k)-ugq(1,i,j,k))/fvw(2*i+1,2*j,2*k) + abs(U1(i-1,j,k)-ugq(1,i-1,j,k))/fvw(2*i-1,2*j,2*k) )
            c2 = 0.5_dp*dt*( abs(U2(i,j,k)-ugq(2,i,j,k))/fvw(2*i,2*j+1,2*k) + abs(U2(i,j-1,k)-ugq(2,i,j-1,k))/fvw(2*i,2*j-1,2*k) )
            c3 = 0.5_dp*dt*( abs(U3(i,j,k)-ugq(3,i,j,k))/fvw(2*i,2*j,2*k+1) + abs(U3(i,j,k-1)-ugq(3,i,j,k-1))/fvw(2*i,2*j,2*k-1) )
            cdir(1) = max(cdir(1), c1);  cdir(2) = max(cdir(2), c2);  cdir(3) = max(cdir(3), c3)
            if (c1 + c2 + c3 > ctot) then
               ctot = c1 + c2 + c3;  cfl_loc = [dom(bb)%b, i, j, k]
            end if
         end do
        end do
       end do
    end do
    cfl_local = ctot
    ctot = global_max(ctot)
    cdir(1) = global_max(cdir(1));  cdir(2) = global_max(cdir(2));  cdir(3) = global_max(cdir(3))
  end subroutine cfl_numbers

  !---------------------------------------------------------------------
  ! min / max / average of |div u| over all cells [1/s]
  !---------------------------------------------------------------------
  subroutine divergence_stats(dmn, dmx, dav)
    real(dp), intent(out) :: dmn, dmx, dav
    integer  :: bb, i, j, k
    real(dp) :: div, s, n
    dmn = huge(1.0_dp);  dmx = 0.0_dp;  s = 0.0_dp;  n = 0.0_dp
    do ib = 1, nown
       bb = own_list(ib)
       call set_block(bb)
       do k = ks, ke
        do j = js, je
         do i = is, ie
            if (cstat(i,j,k) /= 1) cycle                 ! overset fringe / hole
            div = ( Jh(2*i+1,2*j,2*k)/fvw(2*i+1,2*j,2*k)*U1(i,j,k) - Jh(2*i-1,2*j,2*k)/fvw(2*i-1,2*j,2*k)*U1(i-1,j,k)  &
                  + Jh(2*i,2*j+1,2*k)/fvw(2*i,2*j+1,2*k)*U2(i,j,k) - Jh(2*i,2*j-1,2*k)/fvw(2*i,2*j-1,2*k)*U2(i,j-1,k)  &
                  + Jh(2*i,2*j,2*k+1)/fvw(2*i,2*j,2*k+1)*U3(i,j,k) - Jh(2*i,2*j,2*k-1)/fvw(2*i,2*j,2*k-1)*U3(i,j,k-1)  &
                  ) / Jh(2*i,2*j,2*k)
            if (cav_on) div = div - cav_div_target(i,j,k)   ! cavitation: the divergence target is the volume source
                                                            !   (minus the compression of the compressible mixture)
            div = abs(div) * Uch / Lch
            dmn = min(dmn, div);  dmx = max(dmx, div);  s = s + div;  n = n + 1.0_dp
         end do
        end do
       end do
    end do
    dmn = -global_max(-dmn);  dmx = global_max(dmx)
    s = global_sum(s);  n = global_sum(n)
    dav = s / max(n, 1.0_dp)
  end subroutine divergence_stats

  ! maximum vapour fraction over the active cells (cavitation monitor)
  real(dp) function vap_max()
    integer :: bb, i, j, k
    real(dp) :: vm
    vm = 0.0_dp
    do ib = 1, nown
       bb = own_list(ib)
       call set_block(bb)
       do k = ks, ke
        do j = js, je
         do i = is, ie
            if (cstat(i,j,k) /= 1) cycle
            vm = max(vm, vap(i,j,k))
         end do
        end do
       end do
    end do
    vap_max = global_max(vm)
  end function vap_max

  !---------------------------------------------------------------------
  subroutine init_contra(d)
    integer, intent(in) :: d
    integer :: i, j, k, hx, hy, hz
    real(dp) :: u0, v0, w0
    u0 = u_init/Uch;  v0 = v_init/Uch;  w0 = w_init/Uch
    call set_block(d)
    do k = lbound(U1,3), ubound(U1,3)
     do j = lbound(U1,2), ubound(U1,2)
      do i = lbound(U1,1), ubound(U1,1)
         hx = min(max(2*i+1, hlo(1)+1), hhi(1)-1)
         hy = min(max(2*j,   hlo(2)+1), hhi(2)-1)
         hz = min(max(2*k,   hlo(3)+1), hhi(3)-1)
         U1(i,j,k) = fvw(hx,hy,hz)*( mq(1,1,hx,hy,hz)*u0                   &
                   + mq(2,1,hx,hy,hz)*v0 + mq(3,1,hx,hy,hz)*w0 )
         hx = min(max(2*i,   hlo(1)+1), hhi(1)-1)
         hy = min(max(2*j+1, hlo(2)+1), hhi(2)-1)
         U2(i,j,k) = fvw(hx,hy,hz)*( mq(1,2,hx,hy,hz)*u0                   &
                   + mq(2,2,hx,hy,hz)*v0 + mq(3,2,hx,hy,hz)*w0 )
         hy = min(max(2*j,   hlo(2)+1), hhi(2)-1)
         hz = min(max(2*k+1, hlo(3)+1), hhi(3)-1)
         U3(i,j,k) = fvw(hx,hy,hz)*( mq(1,3,hx,hy,hz)*u0                   &
                   + mq(2,3,hx,hy,hz)*v0 + mq(3,3,hx,hy,hz)*w0 )
      end do
     end do
    end do
  end subroutine init_contra

  subroutine allocate_fields
    integer :: bb, i1, i2, j1, j2, k1, k2
    do ib = 1, nown
       bb = own_list(ib)                    ! bb = sub-domain id
       i1 = dom(bb)%r1(1)-ng;  i2 = dom(bb)%r2(1)+ng
       j1 = dom(bb)%r1(2)-ng;  j2 = dom(bb)%r2(2)+ng
       k1 = dom(bb)%r1(3)-ng;  k2 = dom(bb)%r2(3)+ng
       allocate(dom(bb)%U1  (i1:i2, j1:j2, k1:k2))
       allocate(dom(bb)%U2  (i1:i2, j1:j2, k1:k2))
       allocate(dom(bb)%U3  (i1:i2, j1:j2, k1:k2))
       allocate(dom(bb)%Us1 (i1:i2, j1:j2, k1:k2))
       allocate(dom(bb)%Us2 (i1:i2, j1:j2, k1:k2))
       allocate(dom(bb)%Us3 (i1:i2, j1:j2, k1:k2))
       allocate(dom(bb)%u   (i1:i2, j1:j2, k1:k2))
       allocate(dom(bb)%v   (i1:i2, j1:j2, k1:k2))
       allocate(dom(bb)%w   (i1:i2, j1:j2, k1:k2))
       allocate(dom(bb)%un  (i1:i2, j1:j2, k1:k2))
       allocate(dom(bb)%vn  (i1:i2, j1:j2, k1:k2))
       allocate(dom(bb)%wn  (i1:i2, j1:j2, k1:k2))
       allocate(dom(bb)%p   (i1:i2, j1:j2, k1:k2))
       allocate(dom(bb)%us  (i1:i2, j1:j2, k1:k2))
       allocate(dom(bb)%vs  (i1:i2, j1:j2, k1:k2))
       allocate(dom(bb)%ws  (i1:i2, j1:j2, k1:k2))
       allocate(dom(bb)%rhsu(i1:i2, j1:j2, k1:k2))
       allocate(dom(bb)%rhsv(i1:i2, j1:j2, k1:k2))
       allocate(dom(bb)%rhsw(i1:i2, j1:j2, k1:k2))
       allocate(dom(bb)%rhsp(i1:i2, j1:j2, k1:k2))
       allocate(dom(bb)%wrk (i1:i2, j1:j2, k1:k2))
       allocate(dom(bb)%mco(26, i1:i2, j1:j2, k1:k2, 3))
       dom(bb)%mco = 0.0_dp
       allocate(dom(bb)%alpha(i1:i2, j1:j2, k1:k2))
       allocate(dom(bb)%alpn (i1:i2, j1:j2, k1:k2))
       allocate(dom(bb)%tk(i1:i2, j1:j2, k1:k2), dom(bb)%tw(i1:i2, j1:j2, k1:k2), dom(bb)%mut(i1:i2, j1:j2, k1:k2), &
                dom(bb)%wdist(i1:i2, j1:j2, k1:k2), dom(bb)%wfac(i1:i2, j1:j2, k1:k2), dom(bb)%wface(i1:i2, j1:j2, k1:k2))
       dom(bb)%tk = 0.0_dp;  dom(bb)%tw = 0.0_dp;  dom(bb)%mut = 0.0_dp;  dom(bb)%wdist = 1.0e30_dp
       dom(bb)%wfac = 1.0_dp;  dom(bb)%wface = 0
       allocate(dom(bb)%alps(i1:i2, j1:j2, k1:k2), dom(bb)%kap(i1:i2, j1:j2, k1:k2))
       dom(bb)%alps = 0.0_dp;  dom(bb)%kap = 0.0_dp
       dom(bb)%alpha = 1.0_dp;  dom(bb)%alpn = 1.0_dp
       dom(bb)%un = 0.0_dp;  dom(bb)%vn = 0.0_dp;  dom(bb)%wn = 0.0_dp
       allocate(dom(bb)%kr(i1:i2, j1:j2, k1:k2))        ! residual (both solvers)
       dom(bb)%kr = 0.0_dp
       if (poi_solver == 2) then
          allocate(dom(bb)%kv1(i1:i2, j1:j2, k1:k2), dom(bb)%kv2(i1:i2, j1:j2, k1:k2))
          allocate(dom(bb)%kr0(i1:i2, j1:j2, k1:k2))
          allocate(dom(bb)%kpv(i1:i2, j1:j2, k1:k2), dom(bb)%kvv(i1:i2, j1:j2, k1:k2))
          allocate(dom(bb)%ksv(i1:i2, j1:j2, k1:k2), dom(bb)%kt (i1:i2, j1:j2, k1:k2))
          dom(bb)%kv1 = 0.0_dp;  dom(bb)%kv2 = 0.0_dp
          dom(bb)%kr0 = 0.0_dp;  dom(bb)%kpv = 0.0_dp
          dom(bb)%kvv = 0.0_dp;  dom(bb)%ksv = 0.0_dp;  dom(bb)%kt  = 0.0_dp
          ! FGMRES bases over the OWNED cells only
          allocate(dom(bb)%kV(dom(bb)%r1(1):dom(bb)%r2(1), dom(bb)%r1(2):dom(bb)%r2(2), &
                              dom(bb)%r1(3):dom(bb)%r2(3), poi_mkry+1))
          allocate(dom(bb)%kZ(dom(bb)%r1(1):dom(bb)%r2(1), dom(bb)%r1(2):dom(bb)%r2(2), &
                              dom(bb)%r1(3):dom(bb)%r2(3), poi_mkry))
          dom(bb)%kV = 0.0_dp;  dom(bb)%kZ = 0.0_dp
          allocate(dom(bb)%kc(i1:i2, j1:j2, k1:k2, 6));  dom(bb)%kc = 0.0_dp
          allocate(dom(bb)%krr(i1:i2, j1:j2, k1:k2));    dom(bb)%krr = 0.0_dp
          if (poi_precond == 2) then
             allocate(dom(bb)%kpy(i1:i2, j1:j2, k1:k2), dom(bb)%kpt(i1:i2, j1:j2, k1:k2))
             dom(bb)%kpy = 0.0_dp;  dom(bb)%kpt = 0.0_dp
          end if
       end if
    end do
  end subroutine allocate_fields

  !---------------------------------------------------------------------
  subroutine initialize_fields
    integer :: bb, k1, k2, k3, h1, h2, h3
    ! volume fraction first: the properties of everything else depend
    ! on it (two-phase: still water below FreeSurf_Z)
    if (Num_of_Phase == 2) call init_alpha
    do ib = 1, nown
       bb = own_list(ib)                    ! bb = sub-domain id
       ! non-dimensional initial conditions
       dom(bb)%u = u_init / Uch
       dom(bb)%v = v_init / Uch
       dom(bb)%w = w_init / Uch
       dom(bb)%p = p_init / (rho1*Uch*Uch)
       ! two-phase: hydrostatic initial pressure p = p_init + rho g (z_fs - z)
       ! (the first projection then starts from a balanced field instead of
       ! having to build the hydrostatic pressure from zero; update 12)
       ! (single phase too with the incremental pressure form, update 16: p^n of the first predictor
       ! is then already hydrostatic and the fluid at rest stays at rest)
       if ((Num_of_Phase == 2 .or. Prs_form == 2) .and. any(grav /= 0.0_dp)) then
          call set_block(bb)
          do k3 = lbound(p,3), ubound(p,3)
           do k2 = lbound(p,2), ubound(p,2)
            do k1 = lbound(p,1), ubound(p,1)
               h1 = min(max(2*k1, hlo(1)+1), hhi(1)-1)
               h2 = min(max(2*k2, hlo(2)+1), hhi(2)-1)
               h3 = min(max(2*k3, hlo(3)+1), hhi(3)-1)
               p(k1,k2,k3) = p(k1,k2,k3) - rho_h(h1,h2,h3)*gnd(3)*(zfs - Xh(3,h1,h2,h3))
            end do
           end do
          end do
       end if
       dom(bb)%rhsu = 0.0_dp;  dom(bb)%rhsv = 0.0_dp
       dom(bb)%rhsw = 0.0_dp;  dom(bb)%rhsp = 0.0_dp
       dom(bb)%wrk  = 0.0_dp
       dom(bb)%us = dom(bb)%u;  dom(bb)%vs = dom(bb)%v;  dom(bb)%ws = dom(bb)%w
       ! transformed components of the (uniform) initial velocity:
       ! U^q = sum_c m(c,q) * u0_c, evaluated with the local metrics
       call init_contra(bb)
       dom(bb)%Us1 = dom(bb)%U1;  dom(bb)%Us2 = dom(bb)%U2
       dom(bb)%Us3 = dom(bb)%U3
    end do
    call ovs_finish_status                 ! overset: blank the hole cells/faces
    call refresh_uvw(1)
    call contra_to_cart(.false., .true.)   ! Cartesian ghost/BC faces
                                           ! consistent with the BCs
    if (ovst_on) then
       call ovs_update_faces(1)            ! fringe faces from the donors
       call refresh_uvw(1)
       call contra_to_cart(.false., .true.)
       call ovs_update_faces(5)
    end if
    call refresh_all(4)
  end subroutine initialize_fields

end program probCFD_main
