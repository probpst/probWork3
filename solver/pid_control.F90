!=======================================================================
! pid_control.F90  (update 21)
!   PID controllers of the 6-DoF parts and of the force models.
!
!   A controller is a block of "PID_*" lines inside a Body / Child_body /
!   Sub_child_body group of $SixDOF or inside a Force_model group; a
!   "PID_use = 1" line opens a new controller of the group it is written
!   in (two or more controllers per group are allowed, e.g. throttle and
!   steering of the same disk).  The lines of a block:
!
!     PID_use            1 on (0: the block is ignored)
!     PID_type           1 throttle: rotation rate [rev/s] about the axis
!                        2 steering: angle [deg] about the axis
!                        3 direct force: a component of a point force
!                          (Force_model with Force_type = 3 only)
!     PID_Axis           direction of the rotation axis, grid frame at
!                        t = 0 (turns with the mother): grid parts (both
!                        types) and the steering of a disk; the throttle
!                        of a disk spins about Prop_P1 -> Prop_P2
!     PID_Axis_Center    centre of the steering axis of a disk [m]
!                        (default: the disk centre); grid parts: Rot_Center
!     PID_Component      type 3: Force_Fx .. Force_Mz
!     PID_P, PID_I, PID_D   gains: out = out0 + P e + I int(e) dt - D d(meas)/dt
!                        (derivative on the measurement; e = target - meas)
!     PID_Target_body, PID_Target_child_body, PID_Target_sub_child_body
!                        the part whose motion is the controlled variable
!     PID_Target_type    1 location (eta: x y z [m], roll pitch yaw [deg];
!                        a Body in the earth frame, a child relative to its
!                        mother), 2 velocity (nu: u v w [m/s], roll_rate
!                        pitch_rate yaw_rate [deg/s], body frame)
!     PID_Target_var     x y z roll pitch yaw | u v w roll_rate pitch_rate yaw_rate
!     PID_Target_value   up to 32 targets: the next one is taken when
!                        |e| < PID_Target_tol for PID_Target_hold_T [s];
!                        the last one is held to the end
!     PID_Out_init       output before the controller starts (throttle of
!                        a disk: Prop_RPS is used instead)
!     PID_Out_limit      min max of the output
!     PID_Rate_limit     max |d out/dt| (0 = none)
!                        (anti-windup: the error is not integrated while
!                        either limit cuts the output in its direction)
!     PID_Start_T        the controller is frozen at the initial output
!                        before this time [s]
!
!   The owners (sixdof.F90 for the grid parts, force_model.F90 for the
!   force models) look their controllers up with pid_find, feed the
!   measured value of the target part to pid_step every time step and
!   apply the output: a grid part rotates about its axis through
!   Rot_Center (its Motion_DOF is ignored), a disk gets its rps or its
!   steering angle, a point force gets the output added to the chosen
!   component.  Output: Output_data/<case>_<group>_pid<i>.dat every
!   Data_Intvl steps (target index and value, measurement, error, the
!   three terms, the output and the accumulated angle of a throttle).
!   The controller state (integral, previous measurement, output,
!   angle, target index, hold time) is part of the restart file.
!=======================================================================
module pid_control
  use global_vars
  use mpi_vars, only: myrank
  implicit none
  private
  public :: pid_t, pidc, npid, pid_hold, pid_nml, pid_check, pid_find, pid_step, pid_open, pid_export, pid_finalize, &
            pid_nstate, pid_get_state, pid_set_state, pid_summary, pid_axis_unit, pid_var_name, MAXPID

  integer, parameter :: MAXPID = 64, MAXTGT = 32
  real(dp), parameter :: pi_ = 3.14159265358979323846_dp

  type :: pid_t
     ! ---- owner ------------------------------------------------------------
     integer  :: kind = 0                  ! 1 grid part, 2 force model
     integer  :: b = 0, c = 0, s = 0       ! owner part (kind 1): Body / Child_body / Sub_child_body ids
     integer  :: fid = 0                   ! owner force model id (kind 2)
     integer  :: idx = 1                   ! index among the controllers of the same owner
     ! ---- input --------------------------------------------------------------
     integer  :: use = 0, typ = 0
     real(dp) :: axis(3) = 0.0_dp, axc(3) = 0.0_dp
     logical  :: axis_given = .false., axc_given = .false.
     integer  :: comp = 0                  ! type 3: 1..6 = Force_Fx Force_Fy Force_Fz Force_Mx Force_My Force_Mz
     real(dp) :: kp = 0.0_dp, ki = 0.0_dp, kd = 0.0_dp
     integer  :: tb = 0, tc = 0, ts = 0    ! target part ids
     integer  :: ttype = 0, tvar = 0       ! 1 location / 2 velocity; 1..6
     real(dp) :: tval(MAXTGT) = 0.0_dp
     integer  :: ntgt = 0
     real(dp) :: tol = 0.0_dp, hold = 0.0_dp, out0 = 0.0_dp, olim(2) = [-huge(1.0_dp), huge(1.0_dp)]
     real(dp) :: rlim = 0.0_dp, tstart = 0.0_dp
     logical  :: out0_given = .false.
     ! ---- resolved -----------------------------------------------------------
     integer  :: node = 0                  ! owner node (kind 1)
     integer  :: tnode = 0                 ! target node
     ! ---- state (replicated on every rank; in the restart file) -------------
     integer  :: itgt = 1
     real(dp) :: held = 0.0_dp, integ = 0.0_dp, prev = 0.0_dp, out = 0.0_dp, angle = 0.0_dp
     logical  :: started = .false.
     real(dp) :: meas = 0.0_dp, err = 0.0_dp, pterm = 0.0_dp, iterm = 0.0_dp, dterm = 0.0_dp
     integer  :: iused = 1                 ! target index the last error was computed against (export)
     integer  :: iu = -1
     logical  :: open = .false.
  end type pid_t

  type(pid_t), allocatable, target :: pidc(:)
  integer :: npid = 0
  logical :: pid_hold = .false.         ! .true. while the force models are re-evaluated after a restart
  integer :: cur = 0                       ! the controller the PID_* lines are filling

  character(len=10), parameter :: locname(6) = [character(len=10) :: 'x', 'y', 'z', 'roll', 'pitch', 'yaw']
  character(len=10), parameter :: velname(6) = [character(len=10) :: 'u', 'v', 'w', 'roll_rate', 'pitch_rate', 'yaw_rate']
  character(len=8),  parameter :: compname(6) = [character(len=8) :: 'Force_Fx', 'Force_Fy', 'Force_Fz', 'Force_Mx', &
                                                 'Force_My', 'Force_Mz']
  character(len=12), parameter :: typname(3) = [character(len=12) :: 'throttle', 'steering', 'direct force']

contains

  !=====================================================================
  ! Name-list hook: a "PID_*" line of the group (kind, b, c, s / fid).
  ! Called by sixdof_nml and force_nml; ok = .false. for an unknown key.
  !=====================================================================
  subroutine pid_nml(kind, b, c, s, fid, nm, vals, nv, ok)
    integer, intent(in) :: kind, b, c, s, fid, nv
    character(len=*), intent(in) :: nm, vals(:)
    logical, intent(out) :: ok
    integer :: i, io, n
    real(dp) :: r(MAXTGT)
    character(len=32) :: w
    if (.not. allocated(pidc)) allocate(pidc(MAXPID))
    ok = .true.
    r = 0.0_dp
    if (nm /= 'pid_target_var' .and. nm /= 'pid_component') then
       do i = 1, min(nv, MAXTGT)
          read(vals(i), *, iostat=io) r(i)
          if (io /= 0) call perr('cannot read a number from "'//trim(vals(i))//'" for '//trim(nm))
       end do
    end if
    n = nint(r(1))
    if (nm == 'pid_use') then
       if (npid >= MAXPID) call perr('more than 64 PID controllers')
       npid = npid + 1;  cur = npid
       pidc(cur)%kind = kind;  pidc(cur)%b = b;  pidc(cur)%c = c;  pidc(cur)%s = s;  pidc(cur)%fid = fid
       pidc(cur)%use = n
       pidc(cur)%idx = 1
       do i = 1, npid - 1
          if (pidc(i)%kind == kind .and. pidc(i)%b == b .and. pidc(i)%c == c .and. pidc(i)%s == s .and. pidc(i)%fid == fid) &
             pidc(cur)%idx = pidc(cur)%idx + 1
       end do
       return
    end if
    if (cur == 0) call perr('"'//trim(nm)//'" must follow a "PID_use = 1" line of its group')
    ! the line must belong to the group the open controller was written in
    if (pidc(cur)%kind /= kind .or. pidc(cur)%b /= b .or. pidc(cur)%c /= c .or. pidc(cur)%s /= s .or. pidc(cur)%fid /= fid) &
       call perr('"'//trim(nm)//'" is written in another group than its "PID_use" line')
    associate(p => pidc(cur))
      select case (nm)
      case ('pid_type');                   p%typ = n
      case ('pid_axis');                   p%axis = r(1:3);  p%axis_given = .true.
      case ('pid_axis_center');            p%axc = r(1:3);  p%axc_given = .true.
      case ('pid_component')
         w = lower(trim(adjustl(vals(1))))
         p%comp = 0
         do i = 1, 6
            if (w == lower(compname(i))) p%comp = i
         end do
         if (p%comp == 0) call perr('PID_Component must be one of Force_Fx .. Force_Mz')
      case ('pid_p');                      p%kp = r(1)
      case ('pid_i');                      p%ki = r(1)
      case ('pid_d');                      p%kd = r(1)
      case ('pid_target_body');            p%tb = n
      case ('pid_target_child_body');      p%tc = n
      case ('pid_target_sub_child_body');  p%ts = n
      case ('pid_target_type');            p%ttype = n
      case ('pid_target_var')
         w = lower(trim(adjustl(vals(1))))
         p%tvar = 0
         do i = 1, 6
            if (w == lower(locname(i))) then
               p%tvar = i;  if (p%ttype == 0) p%ttype = 1
            end if
            if (w == lower(velname(i))) then
               p%tvar = i;  if (p%ttype == 0) p%ttype = 2
            end if
         end do
         if (w == 'p') p%tvar = 4
         if (w == 'q') p%tvar = 5
         if (w == 'r') p%tvar = 6
         if (p%tvar == 0) call perr('PID_Target_var "'//trim(w)//'" unknown (x y z roll pitch yaw / u v w roll_rate pitch_rate yaw_rate)')
      case ('pid_target_value')
         p%ntgt = min(nv, MAXTGT);  p%tval(1:p%ntgt) = r(1:p%ntgt)
      case ('pid_target_tol');             p%tol = r(1)
      case ('pid_target_hold_t');          p%hold = r(1)
      case ('pid_out_init');               p%out0 = r(1);  p%out0_given = .true.
      case ('pid_out_limit');              p%olim = r(1:2)
      case ('pid_rate_limit');             p%rlim = r(1)
      case ('pid_start_t');                p%tstart = r(1)
      case default
         ok = .false.
      end select
    end associate
  end subroutine pid_nml

  subroutine perr(msg)
    use mpi
    character(len=*), intent(in) :: msg
    integer :: ierr
    if (myrank == 0) write(*,'(A)') ' ERROR (probCFD.nml, PID controller): '//msg
    call MPI_Abort(MPI_COMM_WORLD, 82, ierr)
  end subroutine perr

  pure function lower(s) result(t)
    character(len=*), intent(in) :: s
    character(len=len(s)) :: t
    integer :: i, ic
    t = s
    do i = 1, len(s)
       ic = iachar(s(i:i))
       if (ic >= 65 .and. ic <= 90) t(i:i) = achar(ic + 32)
    end do
  end function lower

  !=====================================================================
  ! Checks common to every controller (the owner-specific ones are made
  ! by the owner): called by the owners after their own set-up.
  !=====================================================================
  subroutine pid_check(i)
    integer, intent(in) :: i
    character(len=48) :: a
    a = pid_owner_name(i)
    associate(p => pidc(i))
      if (p%typ < 1 .or. p%typ > 3) call perr(trim(a)//': PID_type must be 1 (throttle), 2 (steering) or 3 (direct force)')
      if (p%tb < 1) call perr(trim(a)//': PID_Target_body must be given (>= 1)')
      if (p%ts > 0 .and. p%tc == 0) call perr(trim(a)//': PID_Target_sub_child_body needs PID_Target_child_body')
      if (p%ttype < 1 .or. p%ttype > 2) call perr(trim(a)//': PID_Target_type must be 1 (location) or 2 (velocity)')
      if (p%tvar < 1) call perr(trim(a)//': PID_Target_var must be given')
      if (p%ntgt < 1) call perr(trim(a)//': PID_Target_value must be given (at least one value)')
      if (p%olim(2) < p%olim(1)) call perr(trim(a)//': PID_Out_limit must be min max')
      if (p%rlim < 0.0_dp) call perr(trim(a)//': PID_Rate_limit must be >= 0')
      if (p%tol < 0.0_dp .or. p%hold < 0.0_dp) call perr(trim(a)//': PID_Target_tol / PID_Target_hold_T must be >= 0')
      if (p%typ /= 3 .and. p%comp /= 0) call perr(trim(a)//': PID_Component belongs to PID_type = 3')
      if (p%kind == 1 .and. p%typ == 3) call perr(trim(a)//': PID_type = 3 (direct force) belongs to a Force_model of Force_type = 3')
      if (p%kind == 1 .and. .not. p%axis_given) call perr(trim(a)//': PID_Axis must be given for a grid part')
      if (p%axis_given .and. sum(p%axis**2) <= 0.0_dp) call perr(trim(a)//': PID_Axis must be a non-zero vector')
      p%out = min(max(p%out0, p%olim(1)), p%olim(2))
      p%angle = merge(p%out*pi_/180.0_dp, 0.0_dp, p%typ == 2)
    end associate
  end subroutine pid_check

  function pid_owner_name(i) result(s)
    integer, intent(in) :: i
    character(len=48) :: s
    character(len=12) :: a, b, c
    write(a,'(I0)') pidc(i)%b;  write(b,'(I0)') pidc(i)%c;  write(c,'(I0)') pidc(i)%s
    if (pidc(i)%kind == 2) then
       write(a,'(I0)') pidc(i)%fid
       s = 'Force_model '//trim(a)
    else if (pidc(i)%s > 0) then
       s = 'Body '//trim(a)//' Child_body '//trim(b)//' Sub_child_body '//trim(c)
    else if (pidc(i)%c > 0) then
       s = 'Body '//trim(a)//' Child_body '//trim(b)
    else
       s = 'Body '//trim(a)
    end if
    write(a,'(I0)') pidc(i)%idx
    s = trim(s)//' PID '//trim(a)
  end function pid_owner_name

  ! the controllers of an owner: n of them, indices in list(1:n)
  subroutine pid_find(kind, b, c, s, fid, n, list)
    integer, intent(in)  :: kind, b, c, s, fid
    integer, intent(out) :: n, list(:)
    integer :: i
    n = 0;  list = 0
    if (.not. allocated(pidc)) return
    do i = 1, npid
       if (pidc(i)%use /= 1) cycle
       if (pidc(i)%kind /= kind) cycle
       if (kind == 1 .and. .not. (pidc(i)%b == b .and. pidc(i)%c == c .and. pidc(i)%s == s)) cycle
       if (kind == 2 .and. pidc(i)%fid /= fid) cycle
       n = n + 1
       if (n <= size(list)) list(n) = i
    end do
  end subroutine pid_find

  pure function pid_axis_unit(i) result(a)
    integer, intent(in) :: i
    real(dp) :: a(3)
    a = pidc(i)%axis/max(sqrt(sum(pidc(i)%axis**2)), 1.0e-300_dp)
  end function pid_axis_unit

  function pid_var_name(i) result(s)
    integer, intent(in) :: i
    character(len=10) :: s
    if (pidc(i)%ttype == 1) then
       s = locname(pidc(i)%tvar)
    else
       s = velname(pidc(i)%tvar)
    end if
  end function pid_var_name

  !=====================================================================
  ! One controller step at time t with the measured value meas (SI, deg):
  ! returns the output (rev/s, deg, N or N m).  Called every time step by
  ! the owner, on every rank (identical states).
  !=====================================================================
  function pid_step(i, t, meas) result(out)
    integer,  intent(in) :: i
    real(dp), intent(in) :: t, meas
    real(dp) :: out, e, dm, raw, lo, hi, dtt
    associate(p => pidc(i))
      if (pid_hold) then                        ! re-evaluation after a restart: the restored output, no state change
         out = p%out;  return
      end if
      dtt = max(dt_dim, 1.0e-300_dp)
      p%meas = meas;  p%iused = p%itgt
      e = p%tval(p%itgt) - meas
      if (p%ttype == 1 .and. p%tvar >= 4) e = modulo(e + 180.0_dp, 360.0_dp) - 180.0_dp      ! angles: shortest way
      p%err = e
      if (t < p%tstart - 1.0e-12_dp .or. .not. p%started) then
         ! frozen at the initial output; the derivative starts from the first measurement
         p%started = (t >= p%tstart - 1.0e-12_dp)
         p%prev = meas
         p%pterm = 0.0_dp;  p%iterm = 0.0_dp;  p%dterm = 0.0_dp
         out = p%out
      else
         dm = (meas - p%prev)/dtt
         p%prev = meas
         lo = p%olim(1);  hi = p%olim(2)
         if (p%rlim > 0.0_dp) then
            lo = max(lo, p%out - p%rlim*dtt);  hi = min(hi, p%out + p%rlim*dtt)
         end if
         p%pterm = p%kp*e;  p%iterm = p%ki*p%integ;  p%dterm = -p%kd*dm
         raw = p%out0 + p%pterm + p%iterm + p%dterm
         out = min(max(raw, lo), hi)
         ! anti-windup (conditional integration): the error of this step is added to the integral
         ! for the next step unless a limiter - output or rate limit - cut the output in the
         ! error's direction
         if (.not. ((out < raw .and. e > 0.0_dp) .or. (out > raw .and. e < 0.0_dp))) p%integ = p%integ + e*dtt
         p%out = out
      end if
      ! target sequence: the next value once the error stayed inside the tolerance for the hold time
      if (abs(e) <= p%tol) then
         p%held = p%held + dtt
      else
         p%held = 0.0_dp
      end if
      if (p%held >= p%hold - 1.0e-12_dp .and. abs(e) <= p%tol .and. p%itgt < p%ntgt) then
         p%itgt = p%itgt + 1;  p%held = 0.0_dp
      end if
    end associate
  end function pid_step

  !=====================================================================
  ! Output file of controller i: Output_data/<case>_<stem>_pid<idx>.dat
  !=====================================================================
  subroutine pid_open(i, stem, ownertxt)
    integer, intent(in) :: i
    character(len=*), intent(in) :: stem, ownertxt
    character(len=256) :: fname
    character(len=12) :: a
    logical :: ex
    if (myrank /= 0) return
    write(a,'(I0)') pidc(i)%idx
    fname = 'Output_data/'//trim(casename)//'_'//trim(stem)//'_pid'//trim(a)//'.dat'
    inquire(file='Output_data/.', exist=ex)
    if (.not. ex) call execute_command_line('mkdir -p Output_data', wait=.true.)
    associate(p => pidc(i))
      if (run_mode == 1) then
         open(newunit=p%iu, file=trim(fname), status='unknown', position='append', action='write')
         p%open = .true.;  return
      end if
      open(newunit=p%iu, file=trim(fname), status='replace', action='write')
      p%open = .true.
      write(p%iu,'(A)') '# probCFD PID controller (update 21): '//trim(ownertxt)//', type '//trim(typname(p%typ))
      if (p%typ == 1) then
         write(p%iu,'(A,3F9.4,A)') '#   output = rotation rate [rev/s] about the axis', pid_axis_unit(i), &
            ' (grid frame at t = 0); angle = accumulated rotation [deg]'
      else if (p%typ == 2) then
         write(p%iu,'(A,3F9.4,A)') '#   output = steering angle [deg] about the axis', pid_axis_unit(i), ' (grid frame at t = 0)'
      else
         write(p%iu,'(A)') '#   output = '//trim(compname(p%comp))//' [N or N m] added to the point force (mother body frame)'
      end if
      write(p%iu,'(A,I0,A,I0,A,I0,A)') '#   target: '//trim(pid_var_name(i))//' of Body ', p%tb, ' Child_body ', p%tc, &
         ' Sub_child_body ', p%ts, merge(' (location: eta [m, deg])   ', ' (velocity: nu [m/s, deg/s])', p%ttype == 1)
      write(p%iu,'(A,3ES12.4,A,ES12.4,A,ES12.4,A)') '#   gains P I D =', p%kp, p%ki, p%kd, ', tolerance ', p%tol, &
         ', hold ', p%hold, ' s; targets:'
      write(p%iu,'(A,32ES12.4)') '#  ', p%tval(1:p%ntgt)
      write(p%iu,'(A)') '# columns: 1 step  2 time[s]  3 target_index  4 target  5 measured  6 error  7 P_term  8 I_term'// &
         '  9 D_term  10 output  11 angle[deg]'
    end associate
  end subroutine pid_open

  subroutine pid_export(i, istep, t)
    integer,  intent(in) :: i, istep
    real(dp), intent(in) :: t
    if (myrank /= 0 .or. .not. pidc(i)%open) return
    if (mod(istep, data_interval) /= 0 .and. istep /= 0) return
    if (pid_hold .or. (run_mode == 1 .and. istep == 0)) return   ! restart: the step-0 / re-evaluation rows exist
    associate(p => pidc(i))
      write(p%iu,'(I8,ES16.8,I6,8ES16.8)') istep, t, p%iused, p%tval(p%iused), p%meas, p%err, p%pterm, p%iterm, p%dterm, &
         p%out, p%angle*180.0_dp/pi_
      flush(p%iu)
    end associate
  end subroutine pid_export

  subroutine pid_finalize
    integer :: i
    if (myrank /= 0 .or. .not. allocated(pidc)) return
    do i = 1, npid
       if (pidc(i)%open) close(pidc(i)%iu)
    end do
  end subroutine pid_finalize

  subroutine pid_summary
    integer :: i
    if (myrank /= 0 .or. npid == 0) return
    do i = 1, npid
       if (pidc(i)%use /= 1) cycle
       associate(p => pidc(i))
         write(*,'(A,ES10.3,A,ES10.3,A,ES10.3)') '  PID: '//trim(pid_owner_name(i))//', '//trim(typname(p%typ))// &
            ', target '//trim(pid_var_name(i))//' of Body '//trim(itoa(p%tb))//'/child '//trim(itoa(p%tc))// &
            '/sub '//trim(itoa(p%ts))//', P I D ', p%kp, ' ', p%ki, ' ', p%kd
         write(*,'(A,32F10.3)') '       targets', p%tval(1:p%ntgt)
         write(*,'(A,ES10.3,A,ES10.3,A,ES10.3,A,2ES10.3,A,ES10.3)') '       tol ', p%tol, ', hold ', p%hold, ' s, initial out ', &
            p%out, ', limits', p%olim, ', rate limit ', p%rlim
       end associate
    end do
  end subroutine pid_summary

  function itoa(n) result(s)
    integer, intent(in) :: n
    character(len=:), allocatable :: s
    character(len=16) :: a
    write(a,'(I0)') n
    s = trim(a)
  end function itoa

  ! ---- restart state -------------------------------------------------------
  integer function pid_nstate()
    pid_nstate = 0
    if (allocated(pidc)) pid_nstate = npid
  end function pid_nstate

  subroutine pid_get_state(i, st)
    integer,  intent(in)  :: i
    real(dp), intent(out) :: st(8)
    st = [real(pidc(i)%itgt, dp), pidc(i)%held, pidc(i)%integ, pidc(i)%prev, pidc(i)%out, pidc(i)%angle, &
          merge(1.0_dp, 0.0_dp, pidc(i)%started), 0.0_dp]
  end subroutine pid_get_state

  subroutine pid_set_state(i, st)
    integer,  intent(in) :: i
    real(dp), intent(in) :: st(8)
    pidc(i)%itgt = max(1, min(nint(st(1)), pidc(i)%ntgt));  pidc(i)%held = st(2);  pidc(i)%integ = st(3)
    pidc(i)%prev = st(4);  pidc(i)%out = st(5);  pidc(i)%angle = st(6);  pidc(i)%started = (st(7) > 0.5_dp)
  end subroutine pid_set_state
end module pid_control
