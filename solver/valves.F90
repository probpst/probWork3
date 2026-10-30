!=======================================================================
! valves.F90  -  HEART VALVES AS TIME-DEPENDENT POROUS (DARCY) ZONES
!                (heart CFD extension, update H2)
!
!   A valve is a thin slab of cells around its annulus plane (a zone
!   file lists the cells: block i j k weight, written by
!   tools/make_valve_zones.py from the reference grid; the zone moves
!   with the grid).  A closed valve carries the resistance K [1/s]:
!     momentum   : ... + K u = 0  (implicit, on the face unknowns)
!     projection : u^(n+1) = u* - dt/(1 + K dt) grad p / rho
!   so the pressure fluxes G^q of the zone faces are scaled by the
!   mobility m_f = 1 / (1 + K_f dt) in the Poisson operator, the
!   preconditioner and the corrector (consistent), K_f = max of the two
!   cells of the face.  K(t) = K_max c(t), c = 1 closed, 0 open, with
!   smooth ramps of Valve_ramp seconds at the opening / closing times.
!   Mode 1: the open window is given as phase fractions of the period
!   (Valve_open<n> = s_open, s_close; Valve_Period).  Mode 2: pressure
!   driven (not yet implemented).
!
!   $Valves group of probCFD.nml:
!     Num_Valves = 4,  Valve_Period = 0.857142857 [s]
!     Valve_name1 = "mitral", Valve_file1 = "valve_mitral.dat",
!     Valve_K1 = 1.0e6 [1/s], Valve_open1 = 0.43, 1.0, Valve_ramp1 = 0.01 [s], Valve_mode1 = 1
!   Output: Output_data/<case>_valves.dat (step, t, phase, closure c of every valve).
!=======================================================================
module valves
  use mpi
  use global_vars
  use mpi_vars
  implicit none
  private
  public :: valve_nml, valve_init, valve_update, valve_export, valve_finalize, valve_kf, valve_mob, valve_mob_h, valves_on

  integer, parameter :: MAXV = 16
  type :: vdef_t
     character(len=64) :: name = ' ', file = ' '
     real(dp) :: K = 1.0e6_dp                    ! [1/s]
     real(dp) :: s_open(2) = [0.0_dp, 1.0_dp]    ! open window (phase fractions)
     real(dp) :: ramp = 0.01_dp                  ! [s]
     integer  :: mode = 1
     real(dp) :: c = 0.0_dp                      ! closure 0 (open) .. 1 (closed)
     integer  :: ncell = 0                       ! cells in the zone file (all blocks)
  end type vdef_t
  type :: karr_t
     real(dp), allocatable :: w(:,:,:,:)         ! (nv, is-ng:ie+ng, js-ng:je+ng, ks-ng:ke+ng) zone weights
  end type karr_t
  integer :: nv = 0, iu = -1
  type(vdef_t) :: vd(MAXV)
  real(dp) :: period = 0.0_dp
  real(dp), allocatable :: kcur(:)               ! (nv) current resistance, non-dimensional (K tch c)
  type(karr_t), allocatable :: kc(:)             ! per domain
  logical :: valves_on = .false.
  character(len=*), parameter :: out_dir = 'Output_data/'

contains

  !---------------------------------------------------------------------
  ! $Valves keys (read_input.F90 dispatches every key of the group here)
  !---------------------------------------------------------------------
  subroutine valve_nml(nm, vals, nvals, ok)
    character(len=*), intent(in) :: nm, vals(:)
    integer, intent(in) :: nvals
    logical, intent(out) :: ok
    integer :: n, p, io, i
    real(dp) :: r(4)
    character(len=64) :: key
    ok = .true.
    if (nm == 'num_valves') then
       read(vals(1), *, iostat=io) nv
       if (io /= 0 .or. nv < 0 .or. nv > MAXV) call vabort('Num_Valves must be 0..16')
       return
    end if
    if (nm == 'valve_period') then
       read(vals(1), *, iostat=io) period
       if (io /= 0) call vabort('Valve_Period needs a number [s]')
       return
    end if
    p = len_trim(nm)
    do while (p > 1)
       if (nm(p:p) < '0' .or. nm(p:p) > '9') exit
       p = p - 1
    end do
    key = nm(1:p)
    read(nm(p+1:), *, iostat=io) n
    if (io /= 0 .or. n < 1 .or. n > MAXV) then
       ok = .false.;  return
    end if
    select case (trim(key))
    case ('valve_name');  vd(n)%name = unq(vals(1))
    case ('valve_file');  vd(n)%file = unq(vals(1))
    case ('valve_k', 'valve_mode', 'valve_ramp', 'valve_open')
       r = 0.0_dp
       do i = 1, min(nvals, 4)
          read(vals(i), *, iostat=io) r(i)
          if (io /= 0) call vabort('"'//trim(nm)//'": value '//trim(vals(i))//' is not a number')
       end do
       select case (trim(key))
       case ('valve_k');     vd(n)%K = r(1)
       case ('valve_mode');  vd(n)%mode = nint(r(1))
          if (vd(n)%mode /= 1) call vabort('Valve_mode: only 1 (phase table) is implemented')
       case ('valve_ramp');  vd(n)%ramp = r(1)
       case ('valve_open')
          if (nvals < 2) call vabort('"'//trim(nm)//'" needs s_open, s_close (phase fractions 0..1)')
          vd(n)%s_open = r(1:2)
       end select
    case default
       ok = .false.
    end select
  contains
    function unq(s) result(t)
      character(len=*), intent(in) :: s
      character(len=64) :: t
      integer :: l
      t = adjustl(s);  l = len_trim(t)
      if (l >= 2) then
         if ((t(1:1) == '"' .and. t(l:l) == '"') .or. (t(1:1) == "'" .and. t(l:l) == "'")) t = t(2:l-1)
      end if
    end function unq
  end subroutine valve_nml

  subroutine vabort(msg)
    character(len=*), intent(in) :: msg
    integer :: ierr
    if (myrank == 0) write(*,'(A)') ' ERROR (valves): '//msg
    call MPI_Abort(MPI_COMM_WORLD, 82, ierr)
  end subroutine vabort

  !---------------------------------------------------------------------
  ! zone files -> per-domain weights; output file
  !---------------------------------------------------------------------
  subroutine valve_init
    integer :: v, ib, d, b, i, j, k, ios, iuz, nl, nmine, ierr
    real(dp) :: w
    character(len=256) :: line
    logical :: ex
    if (nv == 0) return
    if (period <= 0.0_dp) call vabort('Valve_Period [s] must be > 0')
    allocate(kc(ndom), kcur(nv));  kcur = 0.0_dp
    do ib = 1, nown
       d = own_list(ib)
       allocate(kc(d)%w(nv, dom(d)%r1(1)-ng:dom(d)%r2(1)+ng, dom(d)%r1(2)-ng:dom(d)%r2(2)+ng, dom(d)%r1(3)-ng:dom(d)%r2(3)+ng))
       kc(d)%w = 0.0_dp
    end do
    do v = 1, nv
       if (len_trim(vd(v)%file) == 0) call vabort('Valve_file'//itoa(v)//' is missing')
       if (vd(v)%K <= 0.0_dp) call vabort('Valve_K'//itoa(v)//' must be > 0 [1/s]')
       open(newunit=iuz, file=trim(vd(v)%file), status='old', action='read', iostat=ios)
       if (ios /= 0) call vabort('cannot open the zone file '//trim(vd(v)%file))
       nl = 0;  nmine = 0
       do
          read(iuz, '(A)', iostat=ios) line
          if (ios /= 0) exit
          line = adjustl(line)
          if (len_trim(line) == 0 .or. line(1:1) == '#' .or. line(1:1) == '!') cycle
          read(line, *, iostat=ios) b, i, j, k, w
          if (ios /= 0) call vabort('bad line in '//trim(vd(v)%file)//': '//trim(line))
          if (b < 1 .or. b > nblk) call vabort('block out of range in '//trim(vd(v)%file))
          nl = nl + 1
          do ib = 1, nown
             d = own_list(ib)
             if (dom(d)%b /= b) cycle
             if (i < dom(d)%r1(1) .or. i > dom(d)%r2(1)) cycle
             if (j < dom(d)%r1(2) .or. j > dom(d)%r2(2)) cycle
             if (k < dom(d)%r1(3) .or. k > dom(d)%r2(3)) cycle
             kc(d)%w(v, i, j, k) = max(0.0_dp, min(1.0_dp, w));  nmine = nmine + 1
          end do
       end do
       close(iuz)
       vd(v)%ncell = nl
       call MPI_Allreduce(MPI_IN_PLACE, nmine, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
       if (nmine /= nl .and. myrank == 0) write(*,'(A,I0,A,I0,A)') '  WARNING (valves): ', nl - nmine, ' of the ', nl, &
          ' zone cells of '//trim(vd(v)%name)//' lie outside the grid cells'
    end do
    valves_on = .true.
    call valve_update(t_now)
    if (myrank == 0) then
       inquire(file=out_dir//'.', exist=ex)
       if (.not. ex) call execute_command_line('mkdir -p '//out_dir, wait=.true.)
       open(newunit=iu, file=out_dir//trim(casename)//'_valves.dat', status=merge('unknown', 'replace', run_mode == 1), &
            position='append', action='write')
       if (run_mode /= 1) then
          write(iu,'(A)') '# probCFD valves (H2): closure c of every valve (1 closed, 0 open), resistance K_max c [1/s]'
          write(iu,'(A)',advance='no') '#   step      time[s]        phase'
          do v = 1, nv
             write(iu,'(A)',advance='no') '  '//adjustr(vd(v)%name(1:16))
          end do
          write(iu,'(A)') ' '
       end if
       write(*,'(A)') ' '
       write(*,'(A)') ' ---------------- heart valves: porous zones (H2) ----------------'
       write(*,'(A,I0,A,F10.6,A)') '  valves           : ', nv, ', period ', period, ' s'
       do v = 1, nv
          write(*,'(A,I0,A,ES10.3,A,F6.3,A,F6.3,A,F7.4,A)') '    '//trim(vd(v)%name)//': ', vd(v)%ncell, ' zone cells ('// &
             trim(vd(v)%file)//'), K_max ', vd(v)%K, ' 1/s, open for phase ', vd(v)%s_open(1), ' .. ', vd(v)%s_open(2), &
             ', ramp ', vd(v)%ramp, ' s'
       end do
       write(*,'(A)') '  closed valve: momentum sink K u (implicit), pressure flux x 1/(1 + K dt) (Poisson, preconditioner, corrector)'
    end if
  end subroutine valve_init

  function itoa(i) result(s)
    integer, intent(in) :: i
    character(len=12) :: s
    write(s,'(I0)') i
  end function itoa

  !---------------------------------------------------------------------
  ! closure of every valve at time t (the target time of the step)
  !---------------------------------------------------------------------
  subroutine valve_update(t)
    real(dp), intent(in) :: t
    integer :: v
    real(dp) :: s, u, L, tau, o
    if (.not. valves_on) return
    s = modulo(t/period, 1.0_dp)
    do v = 1, nv
       tau = max(vd(v)%ramp/period, 1.0e-9_dp)
       u = modulo(s - vd(v)%s_open(1), 1.0_dp)                 ! phase since the opening time
       L = modulo(vd(v)%s_open(2) - vd(v)%s_open(1), 1.0_dp)   ! open length
       if (L == 0.0_dp) L = 1.0_dp
       if (u < L) then
          o = sstep(u/tau)                                      ! opening ramp
       else
          o = 1.0_dp - sstep((u - L)/tau)                       ! closing ramp
       end if
       vd(v)%c = 1.0_dp - o
       kcur(v) = vd(v)%K*tch*vd(v)%c                            ! non-dimensional (K tch)
    end do
  contains
    pure real(dp) function sstep(x)
      real(dp), intent(in) :: x
      real(dp) :: y
      y = min(max(x, 0.0_dp), 1.0_dp)
      sstep = y*y*(3.0_dp - 2.0_dp*y)
    end function sstep
  end subroutine valve_update

  !---------------------------------------------------------------------
  ! resistance (non-dim) of the q-face between cell pt and pt + e_q of
  ! the current domain (set_block); mobility 1 / (1 + K dt)
  !---------------------------------------------------------------------
  pure real(dp) function valve_kf(q, pt)
    integer, intent(in) :: q, pt(3)
    integer :: v, e(3)
    real(dp) :: wa, wb
    valve_kf = 0.0_dp
    if (.not. valves_on) return
    e = 0;  e(q) = 1
    do v = 1, nv
       if (kcur(v) == 0.0_dp) cycle
       wa = kc(cur_dom)%w(v, pt(1), pt(2), pt(3))
       wb = kc(cur_dom)%w(v, pt(1)+e(1), pt(2)+e(2), pt(3)+e(3))
       valve_kf = valve_kf + kcur(v)*max(wa, wb)
    end do
  end function valve_kf

  pure real(dp) function valve_mob(q, pt)
    integer, intent(in) :: q, pt(3)
    valve_mob = 1.0_dp
    if (.not. valves_on) return
    valve_mob = 1.0_dp/(1.0_dp + valve_kf(q, pt)*dt)
  end function valve_mob

  ! mobility at half index H of the q-face family (H(q) odd: the face between the cells H(q)/2 and H(q)/2 + 1)
  pure real(dp) function valve_mob_h(q, H)
    integer, intent(in) :: q, H(3)
    integer :: pt(3)
    valve_mob_h = 1.0_dp
    if (.not. valves_on) return
    pt = H/2;  pt(q) = (H(q) - 1)/2
    valve_mob_h = 1.0_dp/(1.0_dp + valve_kf(q, pt)*dt)
  end function valve_mob_h

  subroutine valve_export(istep, t)
    integer,  intent(in) :: istep
    real(dp), intent(in) :: t
    integer :: v
    if (.not. valves_on .or. myrank /= 0) return
    if (mod(istep, data_interval) /= 0 .and. istep /= 0) return
    write(iu,'(I8,ES13.5,F13.6)',advance='no') istep, t, modulo(t/period, 1.0_dp)
    do v = 1, nv
       write(iu,'(F18.6)',advance='no') vd(v)%c
    end do
    write(iu,'(A)') ' '
    flush(iu)
  end subroutine valve_export

  subroutine valve_finalize
    if (valves_on .and. myrank == 0 .and. iu > 0) close(iu)
  end subroutine valve_finalize

end module valves
