!=======================================================================
! force_model.F90  -  force models (update 19)
!
!   Three force models that act as "bodies without a grid" in the 6-DoF
!   hierarchy (name-list group $SixDOF, Force_model = n ... Force_model_End
!   = -1; a force model hangs on a Body, a Child_body or a Sub_child_body
!   through Mother_body / Mother_child_body / Mother_sub_child_body, or
!   on nothing, Mother_body = 0 = fixed in space; it cannot have children):
!
!   Force_type = 1  Hough-Ordway actuator disk: the open-water curves
!       K_T(J) = c0 + c1 J + c2 J^2, K_Q(J) (Prop_KT, Prop_KQ) at J =
!       U_a/(n D) of the liquid-weighted disk-averaged axial inflow
!       (relative to the moving disk) give T = K_T rho n^2 D^4 and Q =
!       K_Q rho n^2 D^5; the axial and tangential body forces follow the
!       Hough & Ordway (1965) radial distributions
!          f_x ~ r* sqrt(1 - r*),   f_t ~ r* sqrt(1 - r*) / (r* (1 - r_h) + r_h),
!          r* = (r/R - r_h)/(1 - r_h),
!       normalised on the DISCRETE disk cells so that sum f_x V = T and
!       sum f_t r V = Q exactly (CFDShip-Iowa's two-pass rescaling, done
!       in one pass here), then weighted by the liquid fraction of the
!       cell (an emerging disk loses thrust).
!   Force_type = 2  blade-element momentum (BEM) body force, pure (no
!       calibration): at every disk cell the local relative inflow of
!       the CFD (U_a axial, U_t tangential, which already contains the
!       induced velocities) and the blade speed omega r give the inflow
!       angle beta = atan2(U_a, omega r - U_t), the angle of attack
!       alpha = phi(r) - beta + 2 f/c (phi = atan((P/D)/(pi r/R)) the
!       geometric pitch angle, alpha_0 = -2 f/c the zero-lift angle of a
!       circular-arc camber line), C_L = C_La alpha (|C_L| <= 1.2), C_D =
!       C_D0 + C_Dk C_L^2, the element loads 1/2 rho W^2 c C_L, C_D (W^2
!       = U_a^2 + (omega r - U_t)^2) with the Prandtl tip / hub loss
!       factor and Z blades, smeared over the annulus of the disk
!       thickness: f_x = Z F q (C_L cos beta - C_D sin beta)/(2 pi r dx),
!       f_t = Z F q (C_L sin beta + C_D cos beta)/(2 pi r dx), weighted
!       by the liquid fraction.  Sections: Prop_Sections stations with
!       Prop_rR, Prop_cD, Prop_PD, Prop_fc, Prop_CD0 (+ optional Prop_CLa,
!       Prop_CDk), linear interpolation in r/R.  K_T, K_Q of the result
!       are reported.
!   Force_type = 3  point force: Fx, Fy, Fz [N] and the optional couples
!       Mx, My, Mz [N m] given as text formulas of the mother's state
!       (x y z roll pitch yaw = eta of the mother part, SNAME: bodies in
!       the earth frame [m, rad], children / sub-children relative to
!       their mother; u v w p q r = nu [m/s, rad/s]; t [s]; constants pi
!       g rho; functions sin cos tan asin acos atan atan2 sqrt exp log
!       abs min max sign step, operators + - * / ^, parentheses).  The
!       force is in the MOTHER's body frame, applied at Force_Point
!       (moving with the mother), and acts on the mother only.  Without
!       a mother (Mother_body = 0) the formulas see the fixed point
!       (x y z = Force_Point [m], the rest 0) and the force is put into
!       the fluid at the cell that contains the point (grid frame).
!
!   Coupling: the disk forces are momentum sources of the fluid (per
!   cell, Cartesian, grid frame; force_source adds (f/rho) to the
!   momentum right-hand side at the faces) and their reaction - the sum
!   over the disk cells of -f V and of -(x - O_m) x f V - is the wrench
!   on the mother part about its centre of rotation O_m (the torque is
!   the moment of the tangential forces), passed to sixdof through
!   sixdof_set_ext and transferred up the hierarchy with the fluid
!   forces (transfer_wrenches).  The forces are evaluated once per time
!   step from the corrected field (after sixdof_forces) and drive the
!   next step (explicit, like the fluid forces).
!   Geometry: Prop_P1, Prop_P2, Force_Point are grid coordinates [m] at
!   t = 0 and move with the mother part (x = Rg (x0 - O0) + O0 + dg).
!   The disk velocity (vG + omG x (x - O)) is subtracted from the fluid
!   velocity for the inflow.
!   Output (rank 0, every Data_Intvl steps): Output_data/<case>_force_model<n>.dat
!   with the wrench in the mother's body frame B_m about O_m, the moment
!   about the model's own centre ((P1+P2)/2 or Force_Point), and for the
!   disks n, J, K_T, K_Q, T, Q, the liquid fraction and the number of
!   disk cells (no mother: grid frame, moment about the own centre).
!=======================================================================
module force_model
  use global_vars
  use mpi_vars
  use sixdof, only: sixdof_node_of, sixdof_node_pose, sixdof_set_ext, sixdof_clear_ext, sixdof_on_flag, sixdof_part_name
  use pid_control, only: pidc, pid_nml, pid_check, pid_find, pid_step, pid_open, pid_export, pid_axis_unit, pid_hold
  implicit none
  private
  public :: fm_on, force_nml, force_in_group, force_init, force_compute, force_source, force_monitor, force_finalize

  integer,  parameter :: MAXF = 64, MAXS = 20, MAXTOK = 256
  real(dp), parameter :: pi = 3.14159265358979323846_dp
  real(dp), parameter :: Fm(3,3) = reshape([-1.0_dp, 0.0_dp, 0.0_dp,  0.0_dp, 1.0_dp, 0.0_dp,  0.0_dp, 0.0_dp, -1.0_dp], [3,3])
  real(dp), parameter :: CL_MAX = 1.2_dp

  ! ---- compiled formula (reverse Polish) -----------------------------------
  type :: prog_t
     integer :: n = 0
     integer :: kind(MAXTOK) = 0          ! 1 number, 2 variable, 3 operator, 4 function
     integer :: code(MAXTOK) = 0          ! variable id / operator id / function id
     real(dp) :: val(MAXTOK) = 0.0_dp
     character(len=256) :: text = ' '
  end type prog_t

  type :: fmod_t
     integer  :: id = 0, typ = 0
     integer  :: mb = 0, mc = 0, ms = 0    ! Mother_body / Mother_child_body / Mother_sub_child_body (0 = none)
     integer  :: node = 0                  ! sixdof node (0 = fixed in space)
     logical  :: defined = .false.
     ! propeller disks (types 1, 2)
     real(dp) :: p1(3) = 0.0_dp, p2(3) = 0.0_dp  ! [m], grid frame at t = 0
     real(dp) :: D = 0.0_dp, hub = 0.0_dp, rps = 0.0_dp
     real(dp) :: kt(3) = 0.0_dp, kq(3) = 0.0_dp
     logical  :: has_kt = .false., has_kq = .false.
     integer  :: Z = 0, ns = 0
     real(dp) :: rR(MAXS) = 0.0_dp, cD(MAXS) = 0.0_dp, PD(MAXS) = 0.0_dp, fc(MAXS) = 0.0_dp
     real(dp) :: cd0(MAXS) = 0.0_dp, cla(MAXS) = 2.0_dp*pi, cdk(MAXS) = 0.01_dp
     logical  :: has_fc = .false., has_cd0 = .false.
     ! point force (type 3)
     real(dp) :: fpt(3) = 0.0_dp
     type(prog_t) :: prg(6)
     logical  :: has_f(6) = .false.
     ! ---- results of the last evaluation (SI) ------------------------------
     real(dp) :: p1c(3) = 0.0_dp, p2c(3) = 0.0_dp, ptc(3) = 0.0_dp   ! current positions [m]
     real(dp) :: FG(3) = 0.0_dp, MG(3) = 0.0_dp   ! wrench on the mother, grid frame, about O_m
     real(dp) :: MGo(3) = 0.0_dp                  ! moment about the own centre, grid frame
     real(dp) :: T = 0.0_dp, Q = 0.0_dp, J = 0.0_dp, ktv = 0.0_dp, kqv = 0.0_dp, liq = 0.0_dp, ua = 0.0_dp
     integer  :: ncell = 0
     integer  :: iu = -1
     ! ---- PID controllers (update 21) --------------------------------------
     integer  :: npid = 0, ipid(4) = 0
     real(dp) :: steer = 0.0_dp                   ! steering angle [rad] about (steer_ax, steer_c), grid frame at t = 0
     real(dp) :: steer_ax(3) = 0.0_dp, steer_c(3) = 0.0_dp
     real(dp) :: pid_add(6) = 0.0_dp              ! direct-force outputs added to Force_Fx .. Force_Mz
  end type fmod_t

  type(fmod_t), allocatable :: fmd(:)
  integer :: nfm = 0, curf = 0
  logical :: in_group = .false.
  logical :: fm_on = .false.
  ! momentum source per owned domain: non-dimensional force per unit volume / (rho1 Uch^2/Lch), cells
  type :: src_t
     real(dp), allocatable :: f(:,:,:,:)   ! (3, is:ie, js:je, ks:ke)
  end type src_t
  type(src_t), allocatable :: src(:)
  real(dp) :: t_last = 0.0_dp

  ! formula variables: 1 x 2 y 3 z 4 roll 5 pitch 6 yaw 7 u 8 v 9 w 10 p 11 q 12 r 13 t 14 pi 15 g 16 rho
  character(len=5), parameter :: varname(16) = [character(len=5) :: 'x', 'y', 'z', 'roll', 'pitch', 'yaw', &
       'u', 'v', 'w', 'p', 'q', 'r', 't', 'pi', 'g', 'rho']
  ! functions: 1 sin 2 cos 3 tan 4 asin 5 acos 6 atan 7 atan2 8 sqrt 9 exp 10 log 11 abs 12 min 13 max 14 sign 15 step 16 neg
  character(len=5), parameter :: fnname(15) = [character(len=5) :: 'sin', 'cos', 'tan', 'asin', 'acos', 'atan', &
       'atan2', 'sqrt', 'exp', 'log', 'abs', 'min', 'max', 'sign', 'step']
  integer, parameter :: fnargs(15) = [1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 2, 2, 1, 1]
  character(len=2), parameter :: cmpname(6) = ['Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz']

contains

  logical function force_in_group()
    force_in_group = in_group
  end function force_in_group

  !=====================================================================
  ! Name-list hook: every line between "Force_model = n" and
  ! "Force_model_End = -1" (called by read_namelist before sixdof_nml)
  !=====================================================================
  subroutine force_nml(nm, vals, nv, ok)
    character(len=*), intent(in) :: nm, vals(:)
    integer, intent(in) :: nv
    logical, intent(out) :: ok
    integer :: i, io, n
    real(dp) :: r(MAXS)
    if (.not. allocated(fmd)) allocate(fmd(MAXF))
    ok = .true.
    if (len(nm) >= 4 .and. in_group) then
       if (nm(1:4) == 'pid_') then                        ! PID controller of this force model (update 21)
          call pid_nml(2, 0, 0, 0, fmd(curf)%id, nm, vals, nv, ok)
          return
       end if
    end if
    r = 0.0_dp
    if (.not. (nm == 'force_fx' .or. nm == 'force_fy' .or. nm == 'force_fz' .or. &
               nm == 'force_mx' .or. nm == 'force_my' .or. nm == 'force_mz')) then
       do i = 1, min(nv, MAXS)
          read(vals(i), *, iostat=io) r(i)
          if (io /= 0) call ferr('cannot read a number from "'//trim(vals(i))//'" for '//trim(nm))
       end do
    end if
    n = nint(r(1))
    if (nm == 'force_model') then
       if (in_group) call ferr('Force_model = n inside another Force_model group (missing Force_model_End = -1)')
       if (n < 1 .or. n > MAXF) call ferr('Force_model must be 1..64')
       do i = 1, nfm
          if (fmd(i)%id == n) call ferr('Force_model id given twice')
       end do
       nfm = nfm + 1;  curf = nfm;  in_group = .true.
       fmd(curf)%id = n;  fmd(curf)%defined = .true.
       return
    end if
    if (.not. in_group) then
       ok = .false.;  return
    end if
    select case (nm)
    case ('force_model_end');          in_group = .false.;  curf = 0
    case ('force_type');               fmd(curf)%typ = n
    case ('mother_body');              fmd(curf)%mb = n
    case ('mother_child_body');        fmd(curf)%mc = n
    case ('mother_sub_child_body');    fmd(curf)%ms = n
    case ('prop_p1');                  fmd(curf)%p1 = r(1:3)
    case ('prop_p2');                  fmd(curf)%p2 = r(1:3)
    case ('prop_d');                   fmd(curf)%D = r(1)
    case ('prop_hub');                 fmd(curf)%hub = r(1)
    case ('prop_rps');                 fmd(curf)%rps = r(1)
    case ('prop_kt');                  fmd(curf)%kt = r(1:3);  fmd(curf)%has_kt = .true.
    case ('prop_kq');                  fmd(curf)%kq = r(1:3);  fmd(curf)%has_kq = .true.
    case ('prop_blades');              fmd(curf)%Z = n
    case ('prop_sections')
       if (n < 2 .or. n > MAXS) call ferr('Prop_Sections must be 2..20')
       fmd(curf)%ns = n
    case ('prop_rr');                  fmd(curf)%rR(1:min(nv,MAXS)) = r(1:min(nv,MAXS))
    case ('prop_cd');                  fmd(curf)%cD(1:min(nv,MAXS)) = r(1:min(nv,MAXS))
    case ('prop_pd');                  fmd(curf)%PD(1:min(nv,MAXS)) = r(1:min(nv,MAXS))
    case ('prop_fc');                  fmd(curf)%fc(1:min(nv,MAXS)) = r(1:min(nv,MAXS));  fmd(curf)%has_fc = .true.
    case ('prop_cd0');                 fmd(curf)%cd0(1:min(nv,MAXS)) = r(1:min(nv,MAXS));  fmd(curf)%has_cd0 = .true.
    case ('prop_cla');                 fmd(curf)%cla(1:min(nv,MAXS)) = r(1:min(nv,MAXS))
    case ('prop_cdk');                 fmd(curf)%cdk(1:min(nv,MAXS)) = r(1:min(nv,MAXS))
    case ('force_point');              fmd(curf)%fpt = r(1:3)
    case ('force_fx');                 call compile_formula(fmd(curf)%prg(1), vals(1));  fmd(curf)%has_f(1) = .true.
    case ('force_fy');                 call compile_formula(fmd(curf)%prg(2), vals(1));  fmd(curf)%has_f(2) = .true.
    case ('force_fz');                 call compile_formula(fmd(curf)%prg(3), vals(1));  fmd(curf)%has_f(3) = .true.
    case ('force_mx');                 call compile_formula(fmd(curf)%prg(4), vals(1));  fmd(curf)%has_f(4) = .true.
    case ('force_my');                 call compile_formula(fmd(curf)%prg(5), vals(1));  fmd(curf)%has_f(5) = .true.
    case ('force_mz');                 call compile_formula(fmd(curf)%prg(6), vals(1));  fmd(curf)%has_f(6) = .true.
    case default
       ok = .false.
    end select
  end subroutine force_nml

  subroutine ferr(msg)
    character(len=*), intent(in) :: msg
    integer :: io
    if (myrank == 0) write(*,'(A)') ' ERROR (probCFD.nml, Force_model): '//msg
    flush(6)
    call MPI_Abort(MPI_COMM_WORLD, 76, io)
  end subroutine ferr

  !=====================================================================
  ! Formula compiler: infix text -> reverse Polish (shunting yard)
  !=====================================================================
  subroutine compile_formula(pg, txt)
    type(prog_t), intent(out) :: pg
    character(len=*), intent(in) :: txt
    character(len=256) :: s
    integer :: i, l, j, io, nst, ntok, k
    integer :: st_kind(MAXTOK), st_code(MAXTOK)      ! operator stack: kind 3 op, 4 function, 5 '(' , 6 ','
    integer :: tk_kind(MAXTOK), tk_code(MAXTOK)      ! token list
    real(dp) :: tk_val(MAXTOK), v
    logical :: expect_operand
    character(len=1) :: ch
    character(len=16) :: word
    s = adjustl(txt)
    l = len_trim(s)
    if (l >= 2) then
       if ((s(1:1) == '"' .and. s(l:l) == '"') .or. (s(1:1) == "'" .and. s(l:l) == "'")) then
          s = s(2:l-1);  l = len_trim(s)
       end if
    end if
    pg%text = s
    ! ---- tokenize ---------------------------------------------------------
    ntok = 0;  i = 1;  expect_operand = .true.
    do while (i <= l)
       ch = s(i:i)
       if (ch == ' ') then
          i = i + 1;  cycle
       end if
       if ((ch >= '0' .and. ch <= '9') .or. ch == '.') then
          j = i
          do while (j <= l)
             ch = s(j:j)
             if ((ch >= '0' .and. ch <= '9') .or. ch == '.' .or. ch == 'e' .or. ch == 'E' .or. ch == 'd' .or. ch == 'D') then
                j = j + 1
             else if ((ch == '+' .or. ch == '-') .and. j > i) then
                if (s(j-1:j-1) == 'e' .or. s(j-1:j-1) == 'E' .or. s(j-1:j-1) == 'd' .or. s(j-1:j-1) == 'D') then
                   j = j + 1
                else
                   exit
                end if
             else
                exit
             end if
          end do
          read(s(i:j-1), *, iostat=io) v
          if (io /= 0) call ferr('bad number "'//s(i:j-1)//'" in formula '//trim(s))
          call push_tok(1, 0, v)
          expect_operand = .false.;  i = j;  cycle
       end if
       if ((ch >= 'a' .and. ch <= 'z') .or. (ch >= 'A' .and. ch <= 'Z') .or. ch == '_') then
          j = i
          do while (j <= l)
             ch = s(j:j)
             if ((ch >= 'a' .and. ch <= 'z') .or. (ch >= 'A' .and. ch <= 'Z') .or. (ch >= '0' .and. ch <= '9') .or. ch == '_') then
                j = j + 1
             else
                exit
             end if
          end do
          word = lower(s(i:j-1))
          k = 0
          do io = 1, size(fnname)
             if (trim(word) == trim(fnname(io))) k = io
          end do
          if (k > 0) then
             call push_tok(4, k, 0.0_dp)              ! function (must be followed by '(')
             expect_operand = .true.
          else
             k = 0
             do io = 1, size(varname)
                if (trim(word) == trim(varname(io))) k = io
             end do
             if (k == 0) call ferr('unknown name "'//trim(word)//'" in formula '//trim(s))
             call push_tok(2, k, 0.0_dp)
             expect_operand = .false.
          end if
          i = j;  cycle
       end if
       select case (ch)
       case ('(');  call push_tok(5, 0, 0.0_dp);  expect_operand = .true.
       case (')');  call push_tok(7, 0, 0.0_dp);  expect_operand = .false.
       case (',');  call push_tok(6, 0, 0.0_dp);  expect_operand = .true.
       case ('+', '-', '*', '/', '^')
          if (expect_operand .and. ch == '-') then
             call push_tok(4, 16, 0.0_dp)            ! unary minus = function neg
          else if (expect_operand .and. ch == '+') then
             ! unary plus: ignore
          else
             select case (ch)
             case ('+');  call push_tok(3, 1, 0.0_dp)
             case ('-');  call push_tok(3, 2, 0.0_dp)
             case ('*');  call push_tok(3, 3, 0.0_dp)
             case ('/');  call push_tok(3, 4, 0.0_dp)
             case ('^');  call push_tok(3, 5, 0.0_dp)
             end select
          end if
          expect_operand = .true.
       case default
          call ferr('unexpected character "'//ch//'" in formula '//trim(s))
       end select
       i = i + 1
    end do
    ! ---- shunting yard ---------------------------------------------------
    pg%n = 0;  nst = 0
    do k = 1, ntok
       select case (tk_kind(k))
       case (1, 2)
          call out_tok(tk_kind(k), tk_code(k), tk_val(k))
       case (4)
          nst = nst + 1;  st_kind(nst) = 4;  st_code(nst) = tk_code(k)
       case (6)
          do while (nst > 0)
             if (st_kind(nst) == 5) exit
             call out_tok(st_kind(nst), st_code(nst), 0.0_dp);  nst = nst - 1
          end do
          if (nst == 0) call ferr('misplaced comma in formula '//trim(s))
       case (3)
          do while (nst > 0)
             if (st_kind(nst) == 5) exit
             if (st_kind(nst) == 4) then
                if (st_code(nst) == 16 .and. tk_code(k) == 5) exit     ! -x^2 = -(x^2)
                call out_tok(4, st_code(nst), 0.0_dp);  nst = nst - 1;  cycle
             end if
             if (prec(st_code(nst)) > prec(tk_code(k)) .or. (prec(st_code(nst)) == prec(tk_code(k)) .and. tk_code(k) /= 5)) then
                call out_tok(3, st_code(nst), 0.0_dp);  nst = nst - 1
             else
                exit
             end if
          end do
          nst = nst + 1;  st_kind(nst) = 3;  st_code(nst) = tk_code(k)
       case (5)
          nst = nst + 1;  st_kind(nst) = 5;  st_code(nst) = 0
       case (7)
          do while (nst > 0)
             if (st_kind(nst) == 5) exit
             call out_tok(st_kind(nst), st_code(nst), 0.0_dp);  nst = nst - 1
          end do
          if (nst == 0) call ferr('unbalanced parentheses in formula '//trim(s))
          nst = nst - 1                                                 ! pop '('
          if (nst > 0) then
             if (st_kind(nst) == 4) then
                call out_tok(4, st_code(nst), 0.0_dp);  nst = nst - 1
             end if
          end if
       end select
    end do
    do while (nst > 0)
       if (st_kind(nst) == 5) call ferr('unbalanced parentheses in formula '//trim(s))
       call out_tok(st_kind(nst), st_code(nst), 0.0_dp);  nst = nst - 1
    end do
    ! ---- check by a trial evaluation ---------------------------------------
    v = eval_formula(pg, [0.1_dp, 0.2_dp, 0.3_dp, 0.01_dp, 0.02_dp, 0.03_dp, 1.0_dp, 0.1_dp, 0.1_dp, 0.01_dp, 0.01_dp, 0.01_dp, &
                          1.0_dp, pi, 9.81_dp, 1000.0_dp], io)
    if (io /= 0) call ferr('formula does not evaluate (operands / arguments mismatch): '//trim(s))
  contains
    subroutine push_tok(kd, cd, vv)
      integer, intent(in) :: kd, cd
      real(dp), intent(in) :: vv
      if (ntok >= MAXTOK) call ferr('formula too long: '//trim(s))
      ntok = ntok + 1;  tk_kind(ntok) = kd;  tk_code(ntok) = cd;  tk_val(ntok) = vv
    end subroutine push_tok
    subroutine out_tok(kd, cd, vv)
      integer, intent(in) :: kd, cd
      real(dp), intent(in) :: vv
      if (pg%n >= MAXTOK) call ferr('formula too long: '//trim(s))
      pg%n = pg%n + 1;  pg%kind(pg%n) = kd;  pg%code(pg%n) = cd;  pg%val(pg%n) = vv
    end subroutine out_tok
    pure integer function prec(op)
      integer, intent(in) :: op
      select case (op)
      case (1, 2);  prec = 1
      case (3, 4);  prec = 2
      case default; prec = 3
      end select
    end function prec
  end subroutine compile_formula

  pure function lower(s) result(t)
    character(len=*), intent(in) :: s
    character(len=len(s)) :: t
    integer :: i, c
    t = s
    do i = 1, len(s)
       c = iachar(s(i:i))
       if (c >= 65 .and. c <= 90) t(i:i) = achar(c + 32)
    end do
  end function lower

  ! evaluate a compiled formula with the 16 variable values
  real(dp) function eval_formula(pg, var, ierr) result(res)
    type(prog_t), intent(in) :: pg
    real(dp), intent(in) :: var(16)
    integer, intent(out) :: ierr
    real(dp) :: st(MAXTOK), a, b
    integer :: k, ns
    ierr = 0;  ns = 0;  res = 0.0_dp
    do k = 1, pg%n
       select case (pg%kind(k))
       case (1);  ns = ns + 1;  st(ns) = pg%val(k)
       case (2);  ns = ns + 1;  st(ns) = var(pg%code(k))
       case (3)
          if (ns < 2) then; ierr = 1; return; end if
          a = st(ns-1);  b = st(ns);  ns = ns - 1
          select case (pg%code(k))
          case (1);  st(ns) = a + b
          case (2);  st(ns) = a - b
          case (3);  st(ns) = a*b
          case (4);  st(ns) = a/sign(max(abs(b), 1.0e-300_dp), b)
          case (5)
             if (a < 0.0_dp .and. abs(b - nint(b)) < 1.0e-12_dp) then
                st(ns) = a**nint(b)
             else
                st(ns) = sign(abs(a)**b, 1.0_dp)
             end if
          end select
       case (4)
          if (pg%code(k) == 16) then
             if (ns < 1) then; ierr = 1; return; end if
             st(ns) = -st(ns);  cycle
          end if
          if (fnargs(pg%code(k)) == 2) then
             if (ns < 2) then; ierr = 1; return; end if
             a = st(ns-1);  b = st(ns);  ns = ns - 1
             select case (pg%code(k))
             case (7);   st(ns) = atan2(a, b)
             case (12);  st(ns) = min(a, b)
             case (13);  st(ns) = max(a, b)
             end select
          else
             if (ns < 1) then; ierr = 1; return; end if
             a = st(ns)
             select case (pg%code(k))
             case (1);   st(ns) = sin(a)
             case (2);   st(ns) = cos(a)
             case (3);   st(ns) = tan(a)
             case (4);   st(ns) = asin(max(-1.0_dp, min(1.0_dp, a)))
             case (5);   st(ns) = acos(max(-1.0_dp, min(1.0_dp, a)))
             case (6);   st(ns) = atan(a)
             case (8);   st(ns) = sqrt(max(a, 0.0_dp))
             case (9);   st(ns) = exp(min(a, 700.0_dp))
             case (10);  st(ns) = log(max(a, 1.0e-300_dp))
             case (11);  st(ns) = abs(a)
             case (14);  st(ns) = sign(1.0_dp, a)
             case (15);  st(ns) = merge(1.0_dp, 0.0_dp, a >= 0.0_dp)
             end select
          end if
       end select
    end do
    if (ns /= 1) then
       ierr = 1;  return
    end if
    res = st(1)
  end function eval_formula

  !=====================================================================
  ! Set-up after sixdof_init / initialize_fields: checks, node mapping,
  ! source arrays, output files
  !=====================================================================
  subroutine force_init
    integer :: n, k, ib, d, ierr, i
    character(len=16) :: a
    character(len=:), allocatable :: mname
    fm_on = (nfm > 0)
    if (.not. fm_on) return
    if (in_group) call ferr('the last Force_model group has no Force_model_End = -1')
    do n = 1, nfm
       associate(f => fmd(n))
         write(a,'(I0)') f%id
         if (f%typ < 1 .or. f%typ > 3) call ferr('Force_model '//trim(a)//': Force_type must be 1, 2 or 3')
         if (f%mb < 0 .or. f%mc < 0 .or. f%ms < 0) call ferr('Force_model '//trim(a)//': negative mother id')
         if (f%mb == 0 .and. (f%mc /= 0 .or. f%ms /= 0)) call ferr('Force_model '//trim(a)//': Mother_child_body / '// &
            'Mother_sub_child_body need Mother_body')
         if (f%ms /= 0 .and. f%mc == 0) call ferr('Force_model '//trim(a)//': Mother_sub_child_body needs Mother_child_body')
         if (f%mb > 0) then
            if (.not. sixdof_on_flag()) call ferr('Force_model '//trim(a)//': Mother_body given but the run has no 6-DoF bodies')
            f%node = sixdof_node_of(f%mb, f%mc, f%ms)
            if (f%node == 0) call ferr('Force_model '//trim(a)//': the mother (Body '//trim(itoa(f%mb))//', Child_body '// &
               trim(itoa(f%mc))//', Sub_child_body '//trim(itoa(f%ms))//') does not exist')
         end if
         select case (f%typ)
         case (1, 2)
            if (f%D <= 0.0_dp) call ferr('Force_model '//trim(a)//': Prop_D must be > 0')
            if (f%hub < 0.0_dp .or. f%hub >= 1.0_dp) call ferr('Force_model '//trim(a)//': Prop_Hub must be in [0, 1)')
            if (sum((f%p2 - f%p1)**2) <= 0.0_dp) call ferr('Force_model '//trim(a)//': Prop_P1 and Prop_P2 must differ')
            if (f%typ == 1) then
               if (.not. (f%has_kt .and. f%has_kq)) call ferr('Force_model '//trim(a)//' (Hough-Ordway): Prop_KT and Prop_KQ needed')
               if (f%Z /= 0 .or. f%ns /= 0 .or. f%has_fc .or. f%has_cd0) call ferr('Force_model '//trim(a)// &
                  ' (Hough-Ordway): the blade data (Prop_Blades, Prop_Sections, ...) belong to Force_type = 2')
            else
               if (f%has_kt .or. f%has_kq) call ferr('Force_model '//trim(a)//' (BEM): Prop_KT / Prop_KQ belong to Force_type = 1')
               if (f%Z < 1) call ferr('Force_model '//trim(a)//' (BEM): Prop_Blades must be >= 1')
               if (f%ns < 2) call ferr('Force_model '//trim(a)//' (BEM): Prop_Sections (2..20) and the station lines needed')
               if (.not. (f%has_fc .and. f%has_cd0)) call ferr('Force_model '//trim(a)//' (BEM): Prop_fc and Prop_CD0 needed')
               do i = 2, f%ns
                  if (f%rR(i) <= f%rR(i-1)) call ferr('Force_model '//trim(a)//' (BEM): Prop_rR must increase')
               end do
               if (any(f%cD(1:f%ns) < 0.0_dp) .or. any(f%PD(1:f%ns) <= 0.0_dp)) call ferr('Force_model '//trim(a)// &
                  ' (BEM): Prop_cD >= 0 and Prop_PD > 0 needed')
            end if
            if (any(f%has_f)) call ferr('Force_model '//trim(a)//': Force_Fx ... belong to Force_type = 3')
         case (3)
            call pid_find(2, 0, 0, 0, f%id, f%npid, f%ipid)
            if (.not. any(f%has_f) .and. f%npid == 0) call ferr('Force_model '//trim(a)// &
               ' (point force): at least one of Force_Fx .. Force_Mz (or a PID controller) needed')
            if (f%D /= 0.0_dp .or. f%rps /= 0.0_dp .or. f%has_kt .or. f%ns /= 0) call ferr('Force_model '//trim(a)// &
               ' (point force): the Prop_* lines belong to Force_type = 1 / 2')
         end select
         ! ---- PID controllers of this model (update 21) ------------------------
         call pid_find(2, 0, 0, 0, f%id, f%npid, f%ipid)
         if (f%npid > 4) call ferr('Force_model '//trim(a)//': more than 4 PID controllers')
         do i = 1, f%npid
            k = f%ipid(i)
            call pid_check(k)
            if (.not. sixdof_on_flag() .or. pidc(k)%tnode == 0) call ferr('Force_model '//trim(a)// &
               ': the PID target part (PID_Target_body / _child_body / _sub_child_body) does not exist')
            select case (pidc(k)%typ)
            case (1)
               if (f%typ == 3) call ferr('Force_model '//trim(a)//': a throttle PID (type 1) needs a disk (Force_type 1 / 2)')
               pidc(k)%out0 = f%rps                                   ! Prop_RPS is the initial rps
               pidc(k)%out = min(max(f%rps, pidc(k)%olim(1)), pidc(k)%olim(2))
               if (pidc(k)%out0_given .and. myrank == 0) write(*,'(A)') '  NOTE (Force_model '//trim(a)// &
                  '): PID_Out_init is ignored for the throttle of a disk - Prop_RPS is the initial rps'
            case (2)
               if (f%typ == 3) call ferr('Force_model '//trim(a)//': a steering PID (type 2) needs a disk (Force_type 1 / 2)')
               if (.not. pidc(k)%axis_given) call ferr('Force_model '//trim(a)//': PID_Axis (steering axis) must be given')
               f%steer_ax = pid_axis_unit(k)
               f%steer_c = 0.5_dp*(f%p1 + f%p2)
               if (pidc(k)%axc_given) f%steer_c = pidc(k)%axc
               f%steer = pidc(k)%out*pi/180.0_dp
            case (3)
               if (f%typ /= 3) call ferr('Force_model '//trim(a)//': a direct-force PID (type 3) needs Force_type = 3')
               if (pidc(k)%comp == 0) call ferr('Force_model '//trim(a)//': PID_Component (Force_Fx .. Force_Mz) must be given')
               f%pid_add(pidc(k)%comp) = pidc(k)%out
            end select
         end do
         f%p1c = f%p1;  f%p2c = f%p2;  f%ptc = f%fpt
       end associate
    end do
    ! momentum-source arrays
    allocate(src(ndom))
    do ib = 1, nown
       d = own_list(ib)
       allocate(src(d)%f(3, dom(d)%r1(1):dom(d)%r2(1), dom(d)%r1(2):dom(d)%r2(2), dom(d)%r1(3):dom(d)%r2(3)))
       src(d)%f = 0.0_dp
    end do
    ! output files
    if (myrank == 0) then
       call make_dir('Output_data')
       do n = 1, nfm
          write(a,'(I0)') fmd(n)%id
          do i = 1, fmd(n)%npid                                   ! PID controllers (update 21)
             call pid_open(fmd(n)%ipid(i), 'force_model'//trim(a), 'Force_model '//trim(a))
          end do
          open(newunit=fmd(n)%iu, file='Output_data/'//trim(casename)//'_force_model'//trim(a)//'.dat', &
               status=merge('unknown', 'replace', run_mode == 1), position='append', action='write')
          if (run_mode /= 1) then
             if (fmd(n)%node > 0) then
                mname = sixdof_part_name(fmd(n)%node)
             else
                mname = 'none (fixed in space)'
             end if
             write(fmd(n)%iu,'(A,I0,A,I0,A)') '# probCFD force model ', fmd(n)%id, ' (update 19): Force_type ', fmd(n)%typ, &
                '  (1 Hough-Ordway disk, 2 BEM disk, 3 point force), mother: '//mname
             if (fmd(n)%node > 0) then
                write(fmd(n)%iu,'(A)') '#   Fx Fy Fz [N], Mx My Mz [N m]: wrench ON THE MOTHER in its body frame B_m (SNAME: x fwd,'// &
                   ' y stbd, z down), moment about the mother''s Rot_Center;'
                write(fmd(n)%iu,'(A)') '#   Mox Moy Moz: moment about the model''s own centre (disk centre / Force_Point) in B_m'
             else
                write(fmd(n)%iu,'(A)') '#   Fx Fy Fz [N], Mx My Mz [N m]: force on the fluid (point force) / on the disk mount'// &
                   ' (reaction, disks) in the GRID frame, moment about the model''s own centre; Mox.. = the same'
             end if
             if (fmd(n)%typ /= 3) write(fmd(n)%iu,'(A)') '#   n [rev/s], J = U_a/(n D) of the liquid-weighted disk-averaged relative'// &
                ' axial inflow, K_T = T/(rho n^2 D^4), K_Q = Q/(rho n^2 D^5), T [N], Q [N m] in the liquid, liquid fraction of'// &
                ' the disk, disk cells'
             if (fmd(n)%typ /= 3) then
                write(fmd(n)%iu,'(A)') '#   step      time[s]            Fx              Fy              Fz              Mx'// &
                   '              My              Mz             Mox             Moy             Moz'// &
                   '           n[1/s]               J             K_T             K_Q            T[N]          Q[N m]        liquid   cells'
             else
                write(fmd(n)%iu,'(A)') '#   step      time[s]            Fx              Fy              Fz              Mx'// &
                   '              My              Mz             Mox             Moy             Moz'
             end if
          end if
       end do
       write(*,'(A)') ' '
       write(*,'(A,I0,A)') '  Force models     : ', nfm, ' (update 19)'
       do n = 1, nfm
          if (fmd(n)%node > 0) then
             mname = sixdof_part_name(fmd(n)%node)
          else
             mname = 'none (fixed in space)'
          end if
          select case (fmd(n)%typ)
          case (1)
             write(*,'(A,I0,A,F8.4,A,F6.3,A,F8.3,A,3F8.4,A,3F8.4,A)') '    model ', fmd(n)%id, ': Hough-Ordway disk, D ', fmd(n)%D, &
                ' m, hub ', fmd(n)%hub, ', n ', fmd(n)%rps, ' 1/s, K_T(J) = ', fmd(n)%kt, ', K_Q(J) = ', fmd(n)%kq, '; mother '//mname
          case (2)
             write(*,'(A,I0,A,F8.4,A,F6.3,A,F8.3,A,I0,A,I0,A)') '    model ', fmd(n)%id, ': BEM disk, D ', fmd(n)%D, &
                ' m, hub ', fmd(n)%hub, ', n ', fmd(n)%rps, ' 1/s, ', fmd(n)%Z, ' blades, ', fmd(n)%ns, ' sections; mother '//mname
          case (3)
             write(*,'(A,I0,A,3F9.4,A)') '    model ', fmd(n)%id, ': point force at ', fmd(n)%fpt, ' m; mother '//mname
             do i = 1, 6
                if (fmd(n)%has_f(i)) write(*,'(A)') '          '//cmpname(i)//' = '//trim(fmd(n)%prg(i)%text)
             end do
          end select
       end do
    end if
    ierr = 0
  end subroutine force_init

  function itoa(i) result(s)
    integer, intent(in) :: i
    character(len=16) :: s
    write(s,'(I0)') i
  end function itoa

  !=====================================================================
  ! Evaluate every force model from the current field: momentum sources
  ! for the next step, wrench on the mother, output.
  !=====================================================================
  subroutine force_compute(istep, t)
    integer,  intent(in) :: istep
    real(dp), intent(in) :: t
    integer :: n, ib, d
    if (.not. fm_on) return
    t_last = t
    do ib = 1, nown
       d = own_list(ib)
       src(d)%f = 0.0_dp
    end do
    call sixdof_clear_ext                      ! the wrenches are added model by model (update 21 fix)
    do n = 1, nfm
       call pid_apply(n, t)                    ! PID controllers (update 21): rps / steering / component
       call move_points(n)
       select case (fmd(n)%typ)
       case (1, 2);  call disk_force(n)
       case (3);     call point_force(n, t)
       end select
       if (fmd(n)%node > 0) call sixdof_set_ext(fmd(n)%node, fmd(n)%FG, fmd(n)%MG)
    end do
    call force_export(istep, t)
  end subroutine force_compute

  ! current positions of the geometry points (mother's pose), [m]
  subroutine move_points(n)
    integer, intent(in) :: n
    real(dp) :: Rg(3,3), O0(3), dg(3), vG(3), omG(3), R(3,3), eta(6), nu(6)
    real(dp) :: q1(3), q2(3), Rs(3,3)
    associate(f => fmd(n))
      q1 = f%p1;  q2 = f%p2
      if (f%steer /= 0.0_dp) then                 ! steered disk (update 21): the axis points turn about the steering axis
         Rs = rodrigues3(f%steer_ax, f%steer)
         q1 = matmul(Rs, f%p1 - f%steer_c) + f%steer_c
         q2 = matmul(Rs, f%p2 - f%steer_c) + f%steer_c
      end if
      if (f%node == 0) then
         f%p1c = q1;  f%p2c = q2;  f%ptc = f%fpt
         return
      end if
      call sixdof_node_pose(f%node, Rg, O0, dg, vG, omG, R, eta, nu)
      f%p1c = matmul(Rg, q1 - O0) + O0 + dg
      f%p2c = matmul(Rg, q2 - O0) + O0 + dg
      f%ptc = matmul(Rg, f%fpt - O0) + O0 + dg
    end associate
  end subroutine move_points

  pure function rodrigues3(a, th) result(R)
    real(dp), intent(in) :: a(3), th
    real(dp) :: R(3,3), K(3,3), c, s
    integer :: i
    c = cos(th);  s = sin(th)
    K = 0.0_dp
    K(1,2) = -a(3);  K(1,3) =  a(2);  K(2,1) =  a(3);  K(2,3) = -a(1);  K(3,1) = -a(2);  K(3,2) =  a(1)
    R = s*K + (1.0_dp - c)*matmul(K, K)
    do i = 1, 3
       R(i,i) = R(i,i) + 1.0_dp
    end do
  end function rodrigues3

  !---------------------------------------------------------------------
  ! PID controllers of force model n (update 21): the measured value of
  ! the target part -> controller step -> rps (throttle), steering angle
  ! (steering) or the direct-force component.
  !---------------------------------------------------------------------
  subroutine pid_apply(n, t)
    integer,  intent(in) :: n
    real(dp), intent(in) :: t
    integer :: j, i, c
    real(dp) :: Rg(3,3), O0(3), dg(3), vG(3), omG(3), R(3,3), eta(6), nu(6), meas, out
    associate(f => fmd(n))
      do j = 1, f%npid
         i = f%ipid(j)
         call sixdof_node_pose(pidc(i)%tnode, Rg, O0, dg, vG, omG, R, eta, nu)
         c = pidc(i)%tvar
         if (pidc(i)%ttype == 1) then
            meas = eta(c)
         else
            meas = nu(c)
         end if
         if (c >= 4) meas = meas*180.0_dp/pi
         out = pid_step(i, t, meas)
         select case (pidc(i)%typ)
         case (1);  f%rps = out
         case (2);  f%steer = out*pi/180.0_dp;  pidc(i)%angle = f%steer
         case (3);  f%pid_add(pidc(i)%comp) = out
         end select
      end do
    end associate
  end subroutine pid_apply

  ! velocity of a point x [m] (grid frame) that moves with the mother of model n [m/s]; the centre O [m]
  subroutine mother_motion(n, Oc, vG, omG)
    integer,  intent(in)  :: n
    real(dp), intent(out) :: Oc(3), vG(3), omG(3)
    real(dp) :: Rg(3,3), O0(3), dg(3), R(3,3), eta(6), nu(6)
    Oc = 0.0_dp;  vG = 0.0_dp;  omG = 0.0_dp
    if (fmd(n)%node == 0) return
    call sixdof_node_pose(fmd(n)%node, Rg, O0, dg, vG, omG, R, eta, nu)
    Oc = O0 + dg
  end subroutine mother_motion

  ! cell-centre Cartesian velocity (non-dim) of cell (i,j,k)
  pure function cell_vel(i, j, k) result(uc)
    integer, intent(in) :: i, j, k
    real(dp) :: uc(3)
    uc(1) = 0.5_dp*(u(i-1,j,k) + u(i,j,k))
    uc(2) = 0.5_dp*(v(i,j-1,k) + v(i,j,k))
    uc(3) = 0.5_dp*(w(i,j,k-1) + w(i,j,k))
  end function cell_vel

  ! liquid fraction of cell (i,j,k)
  pure real(dp) function cell_liq(i, j, k)
    integer, intent(in) :: i, j, k
    cell_liq = 1.0_dp
    if (Num_of_Phase == 2) cell_liq = max(0.0_dp, min(1.0_dp, alpha(i,j,k)))
  end function cell_liq

  !---------------------------------------------------------------------
  ! Actuator disk (types 1 and 2)
  !---------------------------------------------------------------------
  subroutine disk_force(n)
    integer, intent(in) :: n
    integer :: ib, d, i, j, k, ierr, ncl, pass
    real(dp) :: e(3), dx, Rp, rh, xc(3), s, rv(3), r, rr, rs, a, b, V, al, Oc(3), vG(3), omG(3), vdisk(3), urel(3), uc(3)
    real(dp) :: et(3), ua_l, ut_l, loc(5), glb(5), uavg, Jadv, nn, T, Qt, Ax, At, fx, ft, fv(3), FGs(6), FGg(6), Mo(3), dirn
    real(dp) :: om, x, cD, PD, fcam, cd0, cla, cdk, phi, beta, alph, CL, CDs, W2, qdyn, Fl, c, dT, dF, ex(3)
    real(dp) :: sumT, sumQ, sumTg, sumQg
    associate(f => fmd(n))
      e = (f%p2c - f%p1c);  dx = sqrt(sum(e*e));  e = e/dx          ! axis (grid frame), disk thickness [m]
      Rp = 0.5_dp*f%D;  rh = f%hub
      nn = abs(f%rps);  dirn = sign(1.0_dp, f%rps)
      om = 2.0_dp*pi*nn
      call mother_motion(n, Oc, vG, omG)
      f%ncell = 0;  f%T = 0.0_dp;  f%Q = 0.0_dp;  f%J = 0.0_dp;  f%ktv = 0.0_dp;  f%kqv = 0.0_dp;  f%liq = 0.0_dp
      f%FG = 0.0_dp;  f%MG = 0.0_dp;  f%MGo = 0.0_dp
      if (nn <= 0.0_dp) return
      ! ---- pass 1: disk-averaged inflow (and the Hough-Ordway normalisation sums)
      loc = 0.0_dp;  ncl = 0
      do ib = 1, nown
         d = own_list(ib);  call set_block(d)
         do k = ks, ke
          do j = js, je
           do i = is, ie
              if (cstat(i,j,k) /= 1) cycle
              xc = Xh(:,2*i,2*j,2*k)*Lch
              if (.not. in_disk(xc, s, r, rv)) cycle
              ncl = ncl + 1
              V = Jh(2*i,2*j,2*k)*Lch**3
              al = cell_liq(i, j, k)
              vdisk = vG + cross3(omG, xc - Oc)
              uc = cell_vel(i, j, k)*Uch
              urel = uc - vdisk
              loc(1) = loc(1) + sum(urel*e)*V*al
              loc(2) = loc(2) + V*al
              rs = (r/Rp - rh)/(1.0_dp - rh)
              b = rs*sqrt(max(1.0_dp - rs, 0.0_dp))
              a = b/(rs*(1.0_dp - rh) + rh + 1.0e-12_dp)
              loc(3) = loc(3) + b*V
              loc(4) = loc(4) + a*r*V
              loc(5) = loc(5) + V
           end do
          end do
         end do
      end do
      call MPI_Allreduce(loc, glb, 5, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_Allreduce(ncl, f%ncell, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
      if (f%ncell == 0) then
         if (myrank == 0 .and. t_last <= 0.0_dp) write(*,'(A,I0,A)') '  WARNING: force model ', f%id, ': no active cell inside the disk'
         return
      end if
      uavg = 0.0_dp
      if (glb(2) > 0.0_dp) uavg = glb(1)/glb(2)
      f%ua = uavg
      Jadv = uavg/(nn*f%D)
      f%J = Jadv
      Ax = 0.0_dp;  At = 0.0_dp
      if (f%typ == 1) then
         T = (f%kt(1) + f%kt(2)*Jadv + f%kt(3)*Jadv**2)*rho1*nn**2*f%D**4
         Qt = (f%kq(1) + f%kq(2)*Jadv + f%kq(3)*Jadv**2)*rho1*nn**2*f%D**5
         if (glb(3) > 0.0_dp) Ax = T/glb(3)
         if (glb(4) > 0.0_dp) At = Qt/glb(4)
      end if
      ! ---- pass 2: the forces
      FGs = 0.0_dp;  sumT = 0.0_dp;  sumQ = 0.0_dp
      do ib = 1, nown
         d = own_list(ib);  call set_block(d)
         do k = ks, ke
          do j = js, je
           do i = is, ie
              if (cstat(i,j,k) /= 1) cycle
              xc = Xh(:,2*i,2*j,2*k)*Lch
              if (.not. in_disk(xc, s, r, rv)) cycle
              V = Jh(2*i,2*j,2*k)*Lch**3
              al = cell_liq(i, j, k)
              et = dirn*cross3(e, rv)                                   ! tangential unit vector in the rotation direction
              if (f%typ == 1) then
                 rs = (r/Rp - rh)/(1.0_dp - rh)
                 b = rs*sqrt(max(1.0_dp - rs, 0.0_dp))
                 a = b/(rs*(1.0_dp - rh) + rh + 1.0e-12_dp)
                 fx = Ax*b;  ft = At*a                                  ! At = Qt / sum(a r V): sum f_t r V = Qt
              else
                 vdisk = vG + cross3(omG, xc - Oc)
                 uc = cell_vel(i, j, k)*Uch
                 urel = uc - vdisk
                 ua_l = sum(urel*e);  ut_l = sum(urel*et)
                 x = r/Rp
                 call section(f, x, cD, PD, fcam, cd0, cla, cdk)
                 c = cD*f%D
                 phi = atan2(PD, pi*x)
                 beta = atan2(ua_l, om*r - ut_l)
                 alph = phi - beta + 2.0_dp*fcam
                 CL = max(-CL_MAX, min(CL_MAX, cla*alph))
                 CDs = cd0 + cdk*CL*CL
                 W2 = ua_l**2 + (om*r - ut_l)**2
                 qdyn = 0.5_dp*rho1*W2*c
                 Fl = tiploss(f%Z, x, rh, beta)
                 dT = f%Z*Fl*qdyn*(CL*cos(beta) - CDs*sin(beta))            ! per unit radius
                 dF = f%Z*Fl*qdyn*(CL*sin(beta) + CDs*cos(beta))
                 fx = dT/(2.0_dp*pi*max(r, 1.0e-6_dp)*dx)
                 ft = dF/(2.0_dp*pi*max(r, 1.0e-6_dp)*dx)
              end if
              fv = (fx*e + ft*et)*al                                    ! force on the fluid per volume [N/m3]
              src(d)%f(:,i,j,k) = src(d)%f(:,i,j,k) + fv/(rho1*Uch*Uch/Lch)
              FGs(1:3) = FGs(1:3) - fv*V                                ! reaction on the mount
              FGs(4:6) = FGs(4:6) - cross3(xc - Oc, fv*V)
              sumT = sumT + fx*al*V
              sumQ = sumQ + ft*al*V*r
           end do
          end do
         end do
      end do
      call MPI_Allreduce(FGs, FGg, 6, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_Allreduce(sumT, sumTg, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_Allreduce(sumQ, sumQg, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      f%FG = FGg(1:3);  f%MG = FGg(4:6)
      f%T = sumTg;  f%Q = sumQg
      f%ktv = sumTg/(rho1*nn**2*f%D**4);  f%kqv = sumQg/(rho1*nn**2*f%D**5)
      if (glb(5) > 0.0_dp) f%liq = glb(2)/glb(5)
      ! moment about the own centre: M_O = M_own + (c - O) x F  ->  M_own = M_O - (c - O) x F
      ex = 0.5_dp*(f%p1c + f%p2c)
      f%MGo = f%MG - cross3(ex - Oc, f%FG)
    end associate
  contains
    logical function in_disk(xc, s, r, rv)
      real(dp), intent(in)  :: xc(3)
      real(dp), intent(out) :: s, r, rv(3)
      real(dp) :: d3(3)
      in_disk = .false.
      d3 = xc - fmd(n)%p1c
      s = sum(d3*e)
      if (s < 0.0_dp .or. s > dx) return
      rv = d3 - s*e
      r = sqrt(sum(rv*rv))
      if (r > Rp .or. r < rh*Rp) return
      rv = rv/max(r, 1.0e-12_dp)
      in_disk = .true.
    end function in_disk
  end subroutine disk_force

  pure real(dp) function disk_volume(f)
    type(fmod_t), intent(in) :: f
    real(dp) :: dx
    dx = sqrt(sum((f%p2c - f%p1c)**2))
    disk_volume = pi*(0.5_dp*f%D)**2*(1.0_dp - f%hub**2)*dx
  end function disk_volume

  ! blade section data at x = r/R (linear interpolation, clamped)
  pure subroutine section(f, x, cD, PD, fcam, cd0, cla, cdk)
    type(fmod_t), intent(in) :: f
    real(dp), intent(in)  :: x
    real(dp), intent(out) :: cD, PD, fcam, cd0, cla, cdk
    integer :: i
    real(dp) :: w
    if (x <= f%rR(1)) then
       i = 1;  w = 0.0_dp
    else if (x >= f%rR(f%ns)) then
       i = f%ns - 1;  w = 1.0_dp
    else
       i = 1
       do while (i < f%ns - 1 .and. x > f%rR(i+1))
          i = i + 1
       end do
       w = (x - f%rR(i))/(f%rR(i+1) - f%rR(i))
    end if
    cD   = (1.0_dp - w)*f%cD(i)  + w*f%cD(i+1)
    PD   = (1.0_dp - w)*f%PD(i)  + w*f%PD(i+1)
    fcam = (1.0_dp - w)*f%fc(i)  + w*f%fc(i+1)
    cd0  = (1.0_dp - w)*f%cd0(i) + w*f%cd0(i+1)
    cla  = (1.0_dp - w)*f%cla(i) + w*f%cla(i+1)
    cdk  = (1.0_dp - w)*f%cdk(i) + w*f%cdk(i+1)
  end subroutine section

  ! Prandtl tip and hub loss factor
  pure real(dp) function tiploss(Z, x, xh, beta)
    integer,  intent(in) :: Z
    real(dp), intent(in) :: x, xh, beta
    real(dp) :: sb, ft, fh
    sb = max(abs(sin(beta)), 1.0e-3_dp)
    ft = 2.0_dp/pi*acos(max(-1.0_dp, min(1.0_dp, exp(-0.5_dp*Z*(1.0_dp - x)/(x*sb)))))
    fh = 1.0_dp
    if (xh > 0.0_dp) fh = 2.0_dp/pi*acos(max(-1.0_dp, min(1.0_dp, exp(-0.5_dp*Z*(x - xh)/(xh*sb)))))
    tiploss = max(ft*fh, 0.0_dp)
  end function tiploss

  !---------------------------------------------------------------------
  ! Point force (type 3)
  !---------------------------------------------------------------------
  subroutine point_force(n, t)
    integer,  intent(in) :: n
    real(dp), intent(in) :: t
    real(dp) :: var(16), FB(3), MB(3), Rg(3,3), O0(3), dg(3), vG(3), omG(3), R(3,3), eta(6), nu(6), FE(3), ME(3), Oc(3)
    integer :: i, ierr, ib, d, ii, jj, kk, a, b, c
    real(dp) :: xmn, xmx, ymn, ymx, zmn, zmx, xp(3), V
    associate(f => fmd(n))
      var = 0.0_dp
      var(14) = pi;  var(15) = sqrt(sum(grav**2));  var(16) = rho1;  var(13) = t
      if (f%node > 0) then
         call sixdof_node_pose(f%node, Rg, O0, dg, vG, omG, R, eta, nu)
         var(1:6) = eta;  var(7:12) = nu
      else
         var(1:3) = f%fpt
      end if
      FB = 0.0_dp;  MB = 0.0_dp
      do i = 1, 3
         if (f%has_f(i))   FB(i) = eval_formula(f%prg(i), var, ierr)
         if (f%has_f(i+3)) MB(i) = eval_formula(f%prg(i+3), var, ierr)
      end do
      FB = FB + f%pid_add(1:3);  MB = MB + f%pid_add(4:6)          ! PID direct-force outputs (update 21)
      if (f%node > 0) then
         ! mother body frame -> earth -> grid
         FE = matmul(R, FB);  ME = matmul(R, MB)
         f%FG = matmul(Fm, FE);  f%MGo = matmul(Fm, ME)
         Oc = O0 + dg
         f%MG = f%MGo + cross3(f%ptc - Oc, f%FG)
      else
         ! fixed point: force into the fluid (grid frame) at the cell containing the point
         f%FG = FB;  f%MGo = MB;  f%MG = MB
         xp = f%fpt/Lch
         do ib = 1, nown
            d = own_list(ib);  call set_block(d)
            do kk = ks, ke
             do jj = js, je
              do ii = is, ie
                 if (cstat(ii,jj,kk) /= 1) cycle
                 xmn = huge(1.0_dp);  xmx = -xmn;  ymn = xmn;  ymx = -xmn;  zmn = xmn;  zmx = -xmn
                 do c = 2*kk-1, 2*kk+1, 2
                  do b = 2*jj-1, 2*jj+1, 2
                   do a = 2*ii-1, 2*ii+1, 2
                      xmn = min(xmn, Xh(1,a,b,c));  xmx = max(xmx, Xh(1,a,b,c))
                      ymn = min(ymn, Xh(2,a,b,c));  ymx = max(ymx, Xh(2,a,b,c))
                      zmn = min(zmn, Xh(3,a,b,c));  zmx = max(zmx, Xh(3,a,b,c))
                   end do
                  end do
                 end do
                 if (xp(1) < xmn .or. xp(1) >= xmx .or. xp(2) < ymn .or. xp(2) >= ymx .or. xp(3) < zmn .or. xp(3) >= zmx) cycle
                 V = Jh(2*ii,2*jj,2*kk)*Lch**3
                 src(d)%f(:,ii,jj,kk) = src(d)%f(:,ii,jj,kk) + (FB/V)/(rho1*Uch*Uch/Lch)
              end do
             end do
            end do
         end do
      end if
    end associate
  end subroutine point_force

  !---------------------------------------------------------------------
  ! Momentum right-hand side: body force / rho at the c-face between the
  ! cells pt and pt + e_c (non-dimensional acceleration)
  !---------------------------------------------------------------------
  real(dp) function force_source(c, pt, H)
    integer, intent(in) :: c, pt(3), H(3)
    integer :: e(3), p2(3), d
    real(dp) :: s
    force_source = 0.0_dp
    d = cur_dom
    if (.not. allocated(src(d)%f)) return
    e = 0;  e(c) = 1;  p2 = pt + e
    s = 0.0_dp
    if (inside(pt)) s = s + 0.5_dp*src(d)%f(c, pt(1), pt(2), pt(3))
    if (inside(p2)) s = s + 0.5_dp*src(d)%f(c, p2(1), p2(2), p2(3))
    force_source = s/rho_h(H(1), H(2), H(3))
  contains
    pure logical function inside(q)
      integer, intent(in) :: q(3)
      inside = (q(1) >= is .and. q(1) <= ie .and. q(2) >= js .and. q(2) <= je .and. q(3) >= ks .and. q(3) <= ke)
    end function inside
  end function force_source

  !---------------------------------------------------------------------
  ! Output every Data_Intvl steps (rank 0)
  !---------------------------------------------------------------------
  subroutine force_export(istep, t)
    integer,  intent(in) :: istep
    real(dp), intent(in) :: t
    integer :: n, j
    real(dp) :: FB(3), MB(3), MoB(3), Rg(3,3), O0(3), dg(3), vG(3), omG(3), R(3,3), eta(6), nu(6)
    if (myrank /= 0) return
    if (mod(istep, data_interval) /= 0 .and. istep /= 0) return
    if (run_mode == 1 .and. (istep == 0 .or. pid_hold)) return   ! restart (update 21): no duplicate rows
    do n = 1, nfm
       do j = 1, fmd(n)%npid
          call pid_export(fmd(n)%ipid(j), istep, t)
       end do
    end do
    do n = 1, nfm
       associate(f => fmd(n))
         if (f%node > 0) then
            call sixdof_node_pose(f%node, Rg, O0, dg, vG, omG, R, eta, nu)
            FB  = matmul(transpose(R), matmul(Fm, f%FG))
            MB  = matmul(transpose(R), matmul(Fm, f%MG))
            MoB = matmul(transpose(R), matmul(Fm, f%MGo))
         else
            FB = f%FG;  MB = f%MG;  MoB = f%MGo
         end if
         if (f%typ /= 3) then
            write(f%iu,'(I8,ES16.8,9ES16.8,7ES16.8,I8)') istep, t, FB, MB, MoB, f%rps, f%J, f%ktv, f%kqv, f%T, f%Q, f%liq, f%ncell
         else
            write(f%iu,'(I8,ES16.8,9ES16.8)') istep, t, FB, MB, MoB
         end if
         flush(f%iu)
       end associate
    end do
  end subroutine force_export

  subroutine force_monitor
    integer :: n
    if (.not. fm_on .or. myrank /= 0) return
    write(*,'(A)') ' '
    do n = 1, nfm
       associate(f => fmd(n))
         if (f%typ /= 3) then
            write(*,'(A,I0,A,I0,A,F8.3,A,F7.3,A,F8.4,A,F8.4,A,ES11.3,A,ES11.3,A,F6.3,A,I0)') ' [ Force model ', f%id, &
               ' ] type ', f%typ, ': n ', f%rps, '  J ', f%J, '  K_T ', f%ktv, '  K_Q ', f%kqv, '  T ', f%T, ' N  Q ', f%Q, &
               ' N m  liquid ', f%liq, '  cells ', f%ncell
         else
            write(*,'(A,I0,A,3ES11.3,A,3ES11.3)') ' [ Force model ', f%id, ' ] point force F (grid) ', f%FG, '  M about O ', f%MG
         end if
       end associate
    end do
  end subroutine force_monitor

  subroutine force_finalize
    integer :: n
    if (.not. fm_on .or. myrank /= 0) return
    do n = 1, nfm
       if (fmd(n)%iu > 0) close(fmd(n)%iu)
    end do
  end subroutine force_finalize

  pure function cross3(a, b) result(c)
    real(dp), intent(in) :: a(3), b(3)
    real(dp) :: c(3)
    c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
  end function cross3

end module force_model
