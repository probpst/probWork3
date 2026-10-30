!=======================================================================
! stop_criteria.F90 -- calculation-stop criteria (update 10)
!
!   1. NaN: after the momentum solve, after the pressure solve and after
!      the projection / VoF / turbulence step every field (u, v, w, p,
!      alpha, k / nu_tilde, omega / eps) of every own cell is scanned.
!      The first NaN found is reported with the field name, grid block,
!      global cell (face) index and position, and the run stops.
!   2. Residual / update norm: whenever a solver residual |b - A p|/|b|
!      (momentum or pressure) or the pressure update norm max|dp|
!      exceeds Stop_Rsd_Max (name list, default 1e9), the run stops and
!      the block / index / position of the cell that carries the maximum
!      residual (or maximum update) is written to the log.
!   3. CFL: whenever the maximum local CFL number exceeds Stop_CFL_Max
!      (name list, default 1e4), the run stops and the block / index /
!      position of the cell with the maximum CFL number is written.
!
!   Locations are recorded by the solvers on the rank that owns the
!   cell (poi_res_loc, poi_upd_loc, mom_res_loc, cfl_loc in global_vars);
!   the rank whose local maximum equals the global maximum prints it.
!   The stop is collective (every rank reaches the check), so the message
!   is flushed before MPI_Abort.
!=======================================================================
module stop_criteria
  use global_vars
  use mpi_vars
  implicit none
  private
  public :: check_nan, check_residual, check_cfl, record_step, stop_report_input, stop_report_normal, stop_clean

  ! ---- Reason_for_stop.md (update 13) ----------------------------------
  ! every stop of the run - stop criteria, input errors, normal end -
  ! writes Reason_for_stop.md in the working folder: the reason, the
  ! location and values reported by the criterion, the history of the
  ! last steps (CFL, residuals, iterations, liquid volume) and, per kind
  ! of stop, the likely causes and the settings to change.
  integer, parameter :: MXDET = 60, MXHIST = 40
  character(len=200) :: det(MXDET)                 ! detail lines of this rank
  integer :: ndet = 0
  integer :: nhist = 0
  integer :: h_step(MXHIST), h_itm(MXHIST), h_itp(MXHIST)
  real(dp) :: h_t(MXHIST), h_cfl(MXHIST), h_resm(MXHIST), h_resp(MXHIST), h_updp(MXHIST), h_vliq(MXHIST), h_wall(MXHIST)

contains

  ! called by Main after every step (rank 0): ring buffer of the history
  subroutine record_step(istep, t, cfl, itm, resm, itp, resp, updp, vliq, wall)
    integer,  intent(in) :: istep, itm, itp
    real(dp), intent(in) :: t, cfl, resm, resp, updp, vliq, wall
    integer :: m
    if (nhist < MXHIST) then
       nhist = nhist + 1;  m = nhist
    else
       h_step(1:MXHIST-1) = h_step(2:MXHIST);  h_t(1:MXHIST-1) = h_t(2:MXHIST);  h_cfl(1:MXHIST-1) = h_cfl(2:MXHIST)
       h_itm(1:MXHIST-1) = h_itm(2:MXHIST);  h_resm(1:MXHIST-1) = h_resm(2:MXHIST);  h_itp(1:MXHIST-1) = h_itp(2:MXHIST)
       h_resp(1:MXHIST-1) = h_resp(2:MXHIST);  h_updp(1:MXHIST-1) = h_updp(2:MXHIST);  h_vliq(1:MXHIST-1) = h_vliq(2:MXHIST)
       h_wall(1:MXHIST-1) = h_wall(2:MXHIST)
       m = MXHIST
    end if
    h_step(m) = istep;  h_t(m) = t;  h_cfl(m) = cfl;  h_itm(m) = itm;  h_resm(m) = resm
    h_itp(m) = itp;  h_resp(m) = resp;  h_updp(m) = updp;  h_vliq(m) = vliq;  h_wall(m) = wall
  end subroutine record_step

  subroutine add_det(line)
    character(len=*), intent(in) :: line
    if (ndet < MXDET) then
       ndet = ndet + 1;  det(ndet) = line
    end if
  end subroutine add_det

  ! the markdown report (rank 0); kind: 1 NaN, 2 residual, 3 CFL, 4 input error, 0 normal end
  subroutine write_reason_file(kind, title, istep, t)
    integer, intent(in) :: kind, istep
    character(len=*), intent(in) :: title
    real(dp), intent(in) :: t
    integer :: iu, m, r, n, ierr, nrec
    integer :: stat(MPI_STATUS_SIZE)
    character(len=200) :: buf
    ! ---- collect the detail lines of all ranks on rank 0 ----------------
    if (kind >= 1 .and. kind <= 3) then
       if (myrank == 0) then
          do r = 1, nprocs-1
             call MPI_Recv(nrec, 1, MPI_INTEGER, r, 900, MPI_COMM_WORLD, stat, ierr)
             do n = 1, nrec
                call MPI_Recv(buf, 200, MPI_CHARACTER, r, 901, MPI_COMM_WORLD, stat, ierr)
                call add_det(buf)
             end do
          end do
       else
          call MPI_Send(ndet, 1, MPI_INTEGER, 0, 900, MPI_COMM_WORLD, ierr)
          do n = 1, ndet
             call MPI_Send(det(n), 200, MPI_CHARACTER, 0, 901, MPI_COMM_WORLD, ierr)
          end do
       end if
    end if
    if (myrank /= 0) return
    open(newunit=iu, file='Reason_for_stop.md', status='replace', action='write')
    write(iu,'(A)') '# Reason for stop - probCFD_부산 case '//trim(casename)
    write(iu,'(A)') ' '
    select case (kind)
    case (0);  write(iu,'(A)') '## Status: finished normally'
    case (4);  write(iu,'(A)') '## Status: STOPPED before the time marching - input error'
    case (5);  write(iu,'(A)') '## Status: STOPPED before the time marching - missing input'
    case (6);  write(iu,'(A)') '## Status: Mode = 2 (grid-movement check) finished'
    case default;  write(iu,'(A)') '## Status: STOPPED during the time marching - calculation diverged'
    end select
    write(iu,'(A)') ' '
    write(iu,'(A)') '**'//trim(title)//'**'
    write(iu,'(A)') ' '
    if (kind /= 4 .and. kind /= 5) then
       write(iu,'(A,I0,A,I0,A,ES11.4,A)') '* time step ', istep, ' of ', nstep, ' (t = ', t, ' s)'
       write(iu,'(A,ES11.4,A,ES11.4,A)') '* Del_T = ', dt_dim, ' s, End_T = ', End_T, ' s'
    end if
    if (ndet > 0) then
       write(iu,'(A)') ' '
       write(iu,'(A)') '## What the stop criterion found'
       write(iu,'(A)') ' '
       write(iu,'(A)') '```'
       do n = 1, ndet
          write(iu,'(A)') trim(det(n))
       end do
       write(iu,'(A)') '```'
    end if
    if (nhist > 0) then
       write(iu,'(A)') ' '
       write(iu,'(A,I0,A)') '## History of the last ', nhist, ' completed steps'
       write(iu,'(A)') ' '
       write(iu,'(A)') '(the step that stopped is not in this table: it did not complete)'
       write(iu,'(A)') ' '
       write(iu,'(A)') '| step | time [s] | max CFL | it(u*) | res(u*) | it(p) | res(p) | max dp | Vliq [m3] | wall [s] |'
       write(iu,'(A)') '|---|---|---|---|---|---|---|---|---|---|'
       do m = 1, nhist
          write(iu,'(A,I0,A,ES11.4,A,F9.3,A,I0,A,ES9.2,A,I0,A,ES9.2,A,ES9.2,A,ES11.4,A,F8.2,A)') '| ', h_step(m), ' | ', &
             h_t(m), ' | ', h_cfl(m), ' | ', h_itm(m), ' | ', h_resm(m), ' | ', h_itp(m), ' | ', h_resp(m), ' | ', &
             h_updp(m), ' | ', h_vliq(m), ' | ', h_wall(m), ' |'
       end do
    end if
    write(iu,'(A)') ' '
    select case (kind)
    case (1)
       write(iu,'(A)') '## Diagnosis'
       write(iu,'(A)') ' '
       write(iu,'(A)') 'A NaN (not a number) appeared in the field named above; the run cannot continue.'
       write(iu,'(A)') 'A NaN is the end of a divergence that usually shows in the history table first:'
       write(iu,'(A)') 'a CFL number that grows from step to step, a momentum or pressure residual that'
       write(iu,'(A)') 'stops converging, a pressure update norm that jumps.  Look at the reported cell:'
       write(iu,'(A)') ' '
       write(iu,'(A)') '* **Near a wall, a sharp edge or a grid singularity (pole, wake cut)** - the grid'
       write(iu,'(A)') '  cell there is probably too skewed / too small for the time step: reduce `Del_T`,'
       write(iu,'(A)') '  or use the implicit convection (`Cnvt_mode = 2`, `Iner_Itn = 2`).'
       write(iu,'(A)') '* **In the first steps of a run with a moving body or an inlet velocity** - the'
       write(iu,'(A)') '  impulsive start: use `t_ramp` (a few hundred steps) or `Motion_Start_T`.'
       write(iu,'(A)') '* **At an overset fringe (two-phase)** - an air fringe cell next to water donors:'
       write(iu,'(A)') '  check the blank status of the slice through the cell (`Ovst_only = 1`,'
       write(iu,'(A)') '  `tools/blank_slice.py`), a smaller `Del_T`, `Iner_Itn = 2`.'
       write(iu,'(A)') '* **In the free surface region** - a VoF over/undershoot: `vof_cfl_max` <= 0.3,'
       write(iu,'(A)') '  `vof_scheme = 3` (geometric), surface tension only with a resolved interface.'
       write(iu,'(A)') '* **Everywhere at once** - a boundary condition mismatch (check `<case>.bcs`, an'
       write(iu,'(A)') '  outflow face with inflow, a wall flag on an open face) or a wrong `Char_U`/`Char_L`.'
    case (2)
       write(iu,'(A)') '## Diagnosis'
       write(iu,'(A)') ' '
       write(iu,'(A,ES9.2,A)') 'A solver residual or update norm exceeded `Stop_Rsd_Max` = ', Stop_Rsd_Max, ': the linear'
       write(iu,'(A)') 'solver (Jacobi momentum sweeps, or the pressure FGMRES/PCG) diverged instead of converging.'
       write(iu,'(A)') ' '
       write(iu,'(A)') '* **Momentum residual** - the explicit part of the convection is unstable at this CFL'
       write(iu,'(A)') '  (see the history): reduce `Del_T` or `cfl_switch`, or switch to `Cnvt_mode = 2`'
       write(iu,'(A)') '  (implicit) with `Iner_Itn = 2`; `SOR_fct_Vel` < 1 damps the sweeps.'
       write(iu,'(A)') '* **Pressure residual / update** - a singular or ill-conditioned Poisson system:'
       write(iu,'(A)') '  an all-Neumann case without `SW_Ref_Prs`, an overset fringe cycle (check the'
       write(iu,'(A)') '  `.ovst` report for orphans), an extreme density ratio at the interface, or a'
       write(iu,'(A)') '  grid cell with a nearly zero volume at the reported position.  Try'
       write(iu,'(A)') '  `poi_precond = 2`, more `Max_Itr_Prs`, a smaller `Del_T`.'
       write(iu,'(A)') '* The reported cell / face position tells where the divergence starts: inspect'
       write(iu,'(A)') '  the grid quality and the boundary condition there.'
    case (3)
       write(iu,'(A)') '## Diagnosis'
       write(iu,'(A)') ' '
       write(iu,'(A,ES9.2,A)') 'The maximum local CFL number exceeded `Stop_CFL_Max` = ', Stop_CFL_Max, ': the velocity'
       write(iu,'(A)') 'at the reported cell is far larger than the time step allows - a velocity blow-up.'
       write(iu,'(A)') ' '
       write(iu,'(A)') '* If the CFL was already large at step 1: `Del_T` is too large for the smallest'
       write(iu,'(A)') '  cells (an explicit convection direction needs CFL < 0.5 - 1), or the initial /'
       write(iu,'(A)') '  boundary velocity is wrong (units: `Char_U`, `u_cel_initial`, the `.bcs` values).'
       write(iu,'(A)') '* If the CFL grew over the last steps: the same causes as a NaN (see the history):'
       write(iu,'(A)') '  an impulsive start (`t_ramp`), an air cell next to interpolated water values'
       write(iu,'(A)') '  (two-phase overset), a skewed cell near a wall edge.'
       write(iu,'(A)') '* `Stop_CFL_Max` itself can be raised if a large local CFL is expected and the'
       write(iu,'(A)') '  residuals are fine (implicit convection tolerates it).'
    case (4)
       write(iu,'(A)') '## What to do'
       write(iu,'(A)') ' '
       write(iu,'(A)') 'Fix the input named above (probCFD.nml, Input_polyMesh_probCFD.csv, the grid or'
       write(iu,'(A)') 'the .bcs file) and restart; nothing was computed.'
    case default
       write(iu,'(A)') '## Note'
       write(iu,'(A)') ' '
       write(iu,'(A)') 'The run reached End_T.  The last volume solution is in Volume_solutions_grids/.'
    end select
    if (kind /= 4) then
       write(iu,'(A)') ' '
       write(iu,'(A)') '## Settings of this run'
       write(iu,'(A)') ' '
       write(iu,'(A,I0,A,I0,A,I0,A,ES9.2,A,ES9.2)') '* Cnvt_scheme ', Cnvt_scheme, ', Cnvt_mode ', Cnvt_mode, ', Dffs_mode ', Dffs_mode, &
          ', cfl_switch ', cfl_switch, ', vn_switch ', vn_switch
       write(iu,'(A,I0,A,I0,A,I0,A,I0)') '* Iner_Itn ', Iner_Itn, ', Time_scheme ', Time_scheme, ', Num_of_Phase ', Num_of_Phase, &
          ', vof_scheme ', vof_scheme
       write(iu,'(A,I0,A,I0,A,ES9.2,A,F5.2)') '* momentum: Max_Itr_Vel ', mom_maxiter, ', Min_Itr_Vel ', mom_miniter, ', Rsd_Lmt_Vel ', mom_tol, &
          ', SOR ', mom_sor
       write(iu,'(A,I0,A,I0,A,I0,A,ES9.2,A,ES9.2)') '* pressure: poi_solver ', poi_solver, ', poi_precond ', poi_precond, ', Max_Itr_Prs ', &
          poi_maxiter, ', Trl_Lmt_Prs ', poi_tol, ', Rsd_Lmt_Prs ', poi_rtol
       write(iu,'(A,ES9.2,A,ES9.2,A,ES9.2)') '* stop limits: Stop_Rsd_Max ', Stop_Rsd_Max, ', Stop_CFL_Max ', Stop_CFL_Max, ', t_ramp ', t_ramp
       write(iu,'(A,I0,A,I0,A,L1)') '* turbulence model ', TURB_model, ', IDDES ', IDDES_model, ', overset ', ovst_on
    end if
    close(iu)
    write(*,'(A)') ' wrote Reason_for_stop.md'
    flush(6)
  end subroutine write_reason_file

  ! input error (before the time marching): rank 0 writes the report and aborts
  subroutine stop_report_input(msg)
    character(len=*), intent(in) :: msg
    integer :: ierr
    if (myrank == 0) then
       ndet = 0
       call add_det(msg)
       call write_reason_file(4, 'Input error: '//msg, 0, 0.0_dp)
    end if
    call MPI_Abort(MPI_COMM_WORLD, 11, ierr)
  end subroutine stop_report_input

  ! clean stop (all ranks): a missing input file (kind 5) or the end of
  ! the Mode = 2 grid-movement check (kind 6); the message goes to the
  ! log and to Reason_for_stop.md, then MPI is finalized
  subroutine stop_clean(kind, msg)
    integer, intent(in) :: kind
    character(len=*), intent(in) :: msg
    integer :: ierr
    if (myrank == 0) then
       write(*,'(A)') ' '//msg
       ndet = 0
       call add_det(msg)
       call write_reason_file(kind, msg, 0, 0.0_dp)
       flush(6)
    end if
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    call MPI_Finalize(ierr)
    stop
  end subroutine stop_clean

  ! normal end
  subroutine stop_report_normal(istep, t)
    integer,  intent(in) :: istep
    real(dp), intent(in) :: t
    if (myrank /= 0) return
    ndet = 0
    call write_reason_file(0, 'The run completed all time steps (End_T reached).', istep, t)
  end subroutine stop_report_normal

  !---------------------------------------------------------------------
  subroutine check_nan(stage, istep)
    character(len=*), intent(in) :: stage
    integer, intent(in) :: istep
    integer  :: ib, d, i, j, k, nf, found, gfound
    character(len=16) :: names(8)
    real(dp) :: p3(3)
    found = 0
    names = ['u (i-face)  ', 'v (j-face)  ', 'w (k-face)  ', 'p           ', 'alpha       ', 'k/nu_tilde  ', 'omega/eps   ', &
             'T           ']
    outer: do ib = 1, nown
       d = own_list(ib);  call set_block(d)
       nf = 4
       if (Num_of_Phase == 2) nf = 5
       if (turb_on) nf = 7
       do k = ks, ke
        do j = js, je
         do i = is, ie
            if (u(i,j,k) /= u(i,j,k)) then
               found = 1;  call report(1, d, i, j, k);  exit outer
            else if (v(i,j,k) /= v(i,j,k)) then
               found = 1;  call report(2, d, i, j, k);  exit outer
            else if (w(i,j,k) /= w(i,j,k)) then
               found = 1;  call report(3, d, i, j, k);  exit outer
            else if (p(i,j,k) /= p(i,j,k)) then
               found = 1;  call report(4, d, i, j, k);  exit outer
            end if
            if (nf >= 5) then
               if (alpha(i,j,k) /= alpha(i,j,k)) then
                  found = 1;  call report(5, d, i, j, k);  exit outer
               end if
            end if
            if (nf >= 7) then
               if (tk(i,j,k) /= tk(i,j,k)) then
                  found = 1;  call report(6, d, i, j, k);  exit outer
               end if
               if (tw(i,j,k) /= tw(i,j,k)) then
                  found = 1;  call report(7, d, i, j, k);  exit outer
               end if
            end if
            if (energy_on) then
               if (temp(i,j,k) /= temp(i,j,k)) then
                  found = 1;  call report(8, d, i, j, k);  exit outer
               end if
            end if
         end do
        end do
       end do
    end do outer
    gfound = nint(global_max(real(found, dp)))
    if (gfound > 0) call stop_run(1, 'STOP CRITERION 1 (NaN) at step '//itoa(istep)//', '//trim(stage), istep)
  contains
    subroutine report(iv, dd, ii, jj, kk)
      integer, intent(in) :: iv, dd, ii, jj, kk
      character(len=200) :: l1, l2
      p3 = cell_pos(dd, ii, jj, kk)
      write(l1,'(A,I0,A)') 'stop criterion 1: NaN in field '//trim(names(iv))//' - grid block ', dom(dd)%b, &
         ', index (i, j, k) = (' // itoa(ii) // ', ' // itoa(jj) // ', ' // itoa(kk) // ')'
      write(l2,'(A,3ES13.5,A,I0,A)') '   cell position [m] :', p3, '   (found on rank ', myrank, ')'
      write(*,'(A)') ' '
      write(*,'(A)') ' ERROR ('//trim(l1)//')'
      write(*,'(A)') '     '//trim(l2)
      call add_det(l1);  call add_det(l2)
      flush(6)
    end subroutine report
  end subroutine check_nan

  !---------------------------------------------------------------------
  ! residual (|b-Ap|/|b|) and, for the pressure solver, the update norm
  ! max|dp|; kind 1 = momentum, 2 = pressure
  !---------------------------------------------------------------------
  subroutine check_residual(name, res, upd, istep, kind)
    character(len=*), intent(in) :: name
    real(dp), intent(in) :: res, upd
    integer,  intent(in) :: istep, kind
    integer  :: loc(4), c, cmax
    real(dp) :: p3(3), lres
    logical  :: bad
    character(len=200) :: l1
    bad = (res > Stop_Rsd_Max) .or. (upd > Stop_Rsd_Max) .or. (res /= res) .or. (upd /= upd)
    if (.not. bad) return
    write(*,'(A)') ' '
    write(l1,'(A,ES12.4,A,ES12.4,A,ES9.2,A)') 'stop criterion 2: '//trim(name)//' residual |b-Ap|/|b| = ', res, &
       ', update norm max|dp| = ', upd, '  (limit ', Stop_Rsd_Max, ')'
    write(*,'(A)') ' ERROR ('//trim(l1)//')'
    if (myrank == 0) call add_det(l1)
    if (kind == 1) then
       ! momentum: component with the largest rank-local residual
       cmax = maxloc(mom_res_local, 1)
       lres = mom_res_local(cmax)
       if (lres >= res*(1.0_dp - 1.0e-12_dp) .or. res /= res) then
          loc = mom_res_loc(:,cmax)
          p3 = face_pos(cmax, loc)
          write(l1,'(A,I0,A)') '   maximum residual of the '//merge('u', merge('v', 'w', cmax == 2), cmax == 1)// &
             ' equation: grid block ', loc(1), ', face index (i, j, k) = (' // itoa(loc(2)) // ', ' // itoa(loc(3)) // &
             ', ' // itoa(loc(4)) // ')'
          write(*,'(A)') '     '//trim(l1);  call add_det(l1)
          write(l1,'(A,3ES13.5,A,I0,A)') '   position [m] :', p3, '   (rank ', myrank, ')'
          write(*,'(A)') '     '//trim(l1);  call add_det(l1)
       end if
    else
       if (poi_res_local >= res*(1.0_dp - 1.0e-12_dp) .and. poi_res_loc(1) > 0) then
          p3 = cell_pos_blk(poi_res_loc)
          write(l1,'(A,I0,A)') '   maximum residual : grid block ', poi_res_loc(1), ', cell index (i, j, k) = (' // &
             itoa(poi_res_loc(2)) // ', ' // itoa(poi_res_loc(3)) // ', ' // itoa(poi_res_loc(4)) // ')'
          write(*,'(A)') '     '//trim(l1);  call add_det(l1)
          write(l1,'(A,3ES13.5,A,I0,A)') '   position [m] :', p3, '   (rank ', myrank, ')'
          write(*,'(A)') '     '//trim(l1);  call add_det(l1)
       end if
       if (poi_upd_local >= upd*(1.0_dp - 1.0e-12_dp) .and. poi_upd_loc(1) > 0) then
          p3 = cell_pos_blk(poi_upd_loc)
          write(l1,'(A,I0,A)') '   maximum update   : grid block ', poi_upd_loc(1), ', cell index (i, j, k) = (' // &
             itoa(poi_upd_loc(2)) // ', ' // itoa(poi_upd_loc(3)) // ', ' // itoa(poi_upd_loc(4)) // ')'
          write(*,'(A)') '     '//trim(l1);  call add_det(l1)
          write(l1,'(A,3ES13.5,A,I0,A)') '   position [m] :', p3, '   (rank ', myrank, ')'
          write(*,'(A)') '     '//trim(l1);  call add_det(l1)
       end if
    end if
    flush(6)
    call stop_run(2, 'STOP CRITERION 2 (residual / update norm) at step '//itoa(istep), istep)
  end subroutine check_residual

  !---------------------------------------------------------------------
  subroutine check_cfl(cfl, istep)
    real(dp), intent(in) :: cfl
    integer,  intent(in) :: istep
    real(dp) :: p3(3), uc(3)
    integer  :: d
    character(len=200) :: l1
    if (.not. (cfl > Stop_CFL_Max .or. cfl /= cfl)) return
    write(*,'(A)') ' '
    write(l1,'(A,ES12.4,A,ES9.2,A)') 'stop criterion 3: maximum local CFL number = ', cfl, '  (limit ', Stop_CFL_Max, ')'
    write(*,'(A)') ' ERROR ('//trim(l1)//')'
    if (myrank == 0) call add_det(l1)
    if (cfl_local >= cfl*(1.0_dp - 1.0e-12_dp) .and. cfl_loc(1) > 0) then
       p3 = cell_pos_blk(cfl_loc)
       d  = dom_of_block_cell(cfl_loc)
       uc = 0.0_dp
       if (d > 0) then
          call set_block(d)
          uc = [0.5_dp*(u(cfl_loc(2)-1,cfl_loc(3),cfl_loc(4)) + u(cfl_loc(2),cfl_loc(3),cfl_loc(4))), &
                0.5_dp*(v(cfl_loc(2),cfl_loc(3)-1,cfl_loc(4)) + v(cfl_loc(2),cfl_loc(3),cfl_loc(4))), &
                0.5_dp*(w(cfl_loc(2),cfl_loc(3),cfl_loc(4)-1) + w(cfl_loc(2),cfl_loc(3),cfl_loc(4)))]
       end if
       write(l1,'(A,I0,A)') '   maximum CFL at grid block ', cfl_loc(1), ', cell index (i, j, k) = (' // &
          itoa(cfl_loc(2)) // ', ' // itoa(cfl_loc(3)) // ', ' // itoa(cfl_loc(4)) // ')'
       write(*,'(A)') '     '//trim(l1);  call add_det(l1)
       write(l1,'(A,3ES13.5,A,3ES12.4,A,I0,A)') '   position [m] :', p3, '   velocity [m/s] :', uc*Uch, '   (rank ', myrank, ')'
       write(*,'(A)') '     '//trim(l1);  call add_det(l1)
    end if
    flush(6)
    call stop_run(3, 'STOP CRITERION 3 (CFL) at step '//itoa(istep), istep)
  end subroutine check_cfl

  !---------------------------------------------------------------------
  ! helpers
  !---------------------------------------------------------------------
  subroutine stop_run(kind, msg, istep)
    integer, intent(in) :: kind, istep
    character(len=*), intent(in) :: msg
    integer :: ierr
    if (myrank == 0) then
       write(*,'(A)') ' '
       write(*,'(A)') ' ***** CALCULATION STOPPED: '//trim(msg)//' *****'
       flush(6)
    end if
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    call write_reason_file(kind, msg, istep, t_now)     ! Reason_for_stop.md (collective: gathers the details)
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    call MPI_Abort(MPI_COMM_WORLD, 77, ierr)
  end subroutine stop_run

  ! domain id (on this rank) that owns the cell (block, i, j, k); 0 if none
  integer function dom_of_block_cell(loc) result(d)
    integer, intent(in) :: loc(4)
    integer :: ib, dd
    d = 0
    do ib = 1, nown
       dd = own_list(ib)
       if (dom(dd)%b /= loc(1)) cycle
       if (all(loc(2:4) >= dom(dd)%r1) .and. all(loc(2:4) <= dom(dd)%r2)) then
          d = dd;  return
       end if
    end do
  end function dom_of_block_cell

  ! cell centre of the global cell (i,j,k) of block b (node average)
  function cell_pos_blk(loc) result(p3)
    integer, intent(in) :: loc(4)
    real(dp) :: p3(3)
    integer :: b, i, j, k
    b = loc(1);  i = min(max(loc(2),1),blk(b)%nx);  j = min(max(loc(3),1),blk(b)%ny);  k = min(max(loc(4),1),blk(b)%nz)
    p3 = [sum(blk(b)%xn(i:i+1,j:j+1,k:k+1)), sum(blk(b)%yn(i:i+1,j:j+1,k:k+1)), sum(blk(b)%zn(i:i+1,j:j+1,k:k+1))]/8.0_dp
  end function cell_pos_blk

  function cell_pos(d, i, j, k) result(p3)
    integer, intent(in) :: d, i, j, k
    real(dp) :: p3(3)
    p3 = cell_pos_blk([dom(d)%b, i, j, k])
  end function cell_pos

  ! centre of the face (c; i,j,k): the + face of cell (i,j,k) along c
  function face_pos(c, loc) result(p3)
    integer, intent(in) :: c, loc(4)
    real(dp) :: p3(3), pa(3), pb(3)
    integer :: e(3)
    e = 0;  e(c) = 1
    pa = cell_pos_blk(loc)
    pb = cell_pos_blk([loc(1), loc(2)+e(1), loc(3)+e(2), loc(4)+e(3)])
    p3 = 0.5_dp*(pa + pb)
  end function face_pos

  function itoa(n) result(s)
    integer, intent(in) :: n
    character(len=:), allocatable :: s
    character(len=16) :: buf
    write(buf,'(I0)') n
    s = trim(buf)
  end function itoa

end module stop_criteria
