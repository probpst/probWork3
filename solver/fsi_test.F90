!=======================================================================
! fsi_test.F90 - stand-alone checks of the shell solver (update 24)
!   fsi_test cant      linear cantilever plate, tip load (beam theory)
!   fsi_test large     large-deflection cantilever, P L^2/EI = 10 (Mattiasson)
!   fsi_test csm3 [dt] Turek-Hron CSM3 (plane strain, gravity 2 m/s2)
!   fsi_test free      free vibration frequency of a cantilever vs Euler-Bernoulli
! build: make fsi_test
!=======================================================================
program fsi_test
  use fsi_shell
  implicit none
  character(len=32) :: arg
  real(dp) :: dtv, dyp = 0.0_dp
  integer :: nargs
  nargs = command_argument_count()
  arg = 'cant'
  if (nargs >= 1) call get_command_argument(1, arg)
  select case (trim(arg))
  case ('cant');   call test_cantilever
  case ('large');  call test_large
  case ('free');   call test_free
  case ('csm3')
     dtv = 0.005_dp
     if (nargs >= 2) then
        call get_command_argument(2, arg);  read(arg, *) dtv
     end if
     call test_csm3(dtv)
  case default
     print *, 'unknown test ', trim(arg)
  end select

contains

  ! rectangular plate 0..L x -W/2..W/2 (nx x ny elements), clamped at x = 0
  subroutine make_plate(sh, L, W, t, Em, nu, rho, nx, ny, plane_strain)
    type(shell_t), intent(inout) :: sh
    real(dp), intent(in) :: L, W, t, Em, nu, rho
    integer, intent(in) :: nx, ny
    logical, intent(in) :: plane_strain
    integer :: i, j, n, e
    sh%nn = (nx+1)*(ny+1);  sh%ne = nx*ny
    allocate(sh%X0(3,sh%nn), sh%conn(4,sh%ne), sh%te(sh%ne), sh%Ee(sh%ne), sh%nue(sh%ne), sh%rhoe(sh%ne), sh%fixed(6,sh%nn))
    n = 0
    do j = 0, ny
       do i = 0, nx
          n = n + 1
          sh%X0(:,n) = [L*i/real(nx,dp), W*(j/real(ny,dp) - 0.5_dp), 0.0_dp]
       end do
    end do
    e = 0
    do j = 1, ny
       do i = 1, nx
          e = e + 1
          sh%conn(:,e) = [(j-1)*(nx+1)+i, (j-1)*(nx+1)+i+1, j*(nx+1)+i+1, j*(nx+1)+i]
       end do
    end do
    sh%te = t;  sh%Ee = Em;  sh%nue = nu;  sh%rhoe = rho
    sh%fixed = .false.
    do n = 1, sh%nn
       if (sh%X0(1,n) == 0.0_dp) sh%fixed(:,n) = .true.
       if (plane_strain) then       ! cylindrical bending: no lateral displacement, no rotation about x / z
          sh%fixed(2,n) = .true.;  sh%fixed(4,n) = .true.;  sh%fixed(6,n) = .true.
       end if
    end do
  end subroutine make_plate

  subroutine test_cantilever
    type(shell_t) :: sh
    real(dp), allocatable :: f(:,:), vm(:)
    real(dp) :: L, W, t, E, nu, P, I2, dref, wtip
    integer :: n, ntip
    logical :: conv
    L = 1.0_dp;  W = 0.1_dp;  t = 0.01_dp;  E = 2.0e11_dp;  nu = 0.3_dp;  P = 100.0_dp
    call make_plate(sh, L, W, t, E, nu, 7800.0_dp, 20, 2, .false.)
    call shell_setup(sh)
    allocate(f(6,sh%nn), vm(sh%nn));  f = 0.0_dp
    ntip = 0
    do n = 1, sh%nn
       if (sh%X0(1,n) == L) ntip = ntip + 1
    end do
    do n = 1, sh%nn
       if (sh%X0(1,n) == L) f(3,n) = P/ntip
    end do
    call shell_static(sh, f, conv)
    wtip = 0.0_dp
    do n = 1, sh%nn
       if (sh%X0(1,n) == L) wtip = wtip + sh%u1(3,n)/ntip
    end do
    I2 = W*t**3/12.0_dp;  dref = P*L**3/(3.0_dp*E*I2)
    call shell_stress(sh, vm)
    print '(A,I0,A,I0,A,ES10.3)', ' cantilever 20 x 2, Newton its ', sh%newton_it, ' (conv ', merge(1,0,conv), '), band ', real(sh%kd)
    print '(A,F10.6,A,F10.6,A,F7.3,A)', ' tip deflection ', wtip, ' m   beam theory P L^3 / 3 E I = ', dref, '   ratio ', wtip/dref
    print '(A,ES11.4,A,ES11.4)', ' max von Mises at the root ', maxval(vm), ' Pa   beam: 6 P L / (W t^2) = ', 6.0_dp*P*L/(W*t*t)
    print '(A,F10.4,A,F10.4)', ' first frequency [Hz] ', sh%omega1/(2.0_dp*pi_()), '   Euler-Bernoulli 0.5596 sqrt(EI/(rho A L^4)) = ', &
       0.5596_dp*sqrt(E*I2/(7800.0_dp*W*t*L**4))
  end subroutine test_cantilever

  subroutine test_free
    type(shell_t) :: sh
    real(dp) :: L, W, t, E, nu, I2
    L = 0.35_dp;  W = 0.02_dp;  t = 0.02_dp;  E = 1.4e6_dp;  nu = 0.4_dp
    call make_plate(sh, L, W, t, E, nu, 1000.0_dp, 40, 2, .true.)
    call shell_setup(sh)
    I2 = W*t**3/12.0_dp
    print '(A,F10.4,A,F10.4)', ' CSM3 rod, plane strain: f1 [Hz] ', sh%omega1/(2.0_dp*pi_()), '   beam (E/(1-nu^2)) ', &
       0.5596_dp*sqrt(E/(1.0_dp - nu*nu)*I2/(1000.0_dp*W*t*L**4))
  end subroutine test_free

  subroutine test_large
    type(shell_t) :: sh
    real(dp), allocatable :: f(:,:)
    real(dp) :: L, W, t, E, nu, P, I2, wtip, utip
    integer :: n, ntip, step, nstep
    logical :: conv
    L = 1.0_dp;  W = 0.1_dp;  t = 0.01_dp;  E = 2.0e11_dp;  nu = 0.0_dp
    I2 = W*t**3/12.0_dp
    P = 10.0_dp*E*I2/L**2
    call make_plate(sh, L, W, t, E, nu, 7800.0_dp, 40, 2, .false.)
    sh%max_newton = 50
    call shell_setup(sh)
    allocate(f(6,sh%nn))
    ntip = 0
    do n = 1, sh%nn
       if (sh%X0(1,n) == L) ntip = ntip + 1
    end do
    nstep = 20
    do step = 1, nstep
       f = 0.0_dp
       do n = 1, sh%nn
          if (sh%X0(1,n) == L) f(3,n) = P*step/real(nstep,dp)/ntip
       end do
       call shell_static(sh, f, conv)
       call shell_accept(sh, f)
       wtip = 0.0_dp;  utip = 0.0_dp
       do n = 1, sh%nn
          if (sh%X0(1,n) == L) then
             wtip = wtip + sh%u1(3,n)/ntip;  utip = utip - sh%u1(1,n)/ntip
          end if
       end do
       print '(A,F6.2,A,F8.4,A,F8.4,A,I3,A,I0)', ' P L^2/EI = ', 10.0_dp*step/nstep, '   w/L = ', wtip/L, '   u/L = ', utip/L, &
          '   Newton ', sh%newton_it, '  conv ', merge(1,0,conv)
    end do
    print '(A)', ' Mattiasson (1981), P L^2/EI = 10: w/L = 0.8109, u/L = 0.5552'
  end subroutine test_large

  subroutine test_csm3(dt)
    real(dp), intent(in) :: dt
    type(shell_t) :: sh
    real(dp), allocatable :: f(:,:)
    real(dp) :: L, W, t, E, nu, rho, g, tnow, tend, dy, dx, dymin, dymax, dxmin, dxmax, tlast, tprev
    real(dp) :: ekin, estr
    integer :: n, ntip, step, nstep, ncross, iu, its
    logical :: conv
    L = 0.35_dp;  W = 0.02_dp;  t = 0.02_dp;  E = 1.4e6_dp;  nu = 0.4_dp;  rho = 1000.0_dp;  g = 2.0_dp
    call make_plate(sh, L, W, t, E, nu, rho, 40, 2, .true.)
    sh%rho_inf = 0.8_dp;  sh%max_newton = 30
    call shell_setup(sh)
    allocate(f(6,sh%nn));  f = 0.0_dp
    do n = 1, sh%nn
       f(3,n) = -sh%mass(3,n)*g               ! gravity along -z (the plate's normal) = the benchmark's -y
    end do
    ntip = 0
    do n = 1, sh%nn
       if (sh%X0(1,n) == L) ntip = ntip + 1
    end do
    tend = 10.0_dp;  nstep = nint(tend/dt)
    open(newunit=iu, file='csm3_tip.dat', status='replace')
    write(iu,'(A)') '# t  dx  dy (tip, point A)'
    dymin = 0.0_dp;  dymax = -1.0_dp;  dxmin = 0.0_dp;  dxmax = -1.0_dp;  ncross = 0;  tlast = 0.0_dp;  tprev = -1.0_dp
    its = 0
    do step = 1, nstep
       tnow = step*dt
       call shell_step(sh, dt, f, conv)
       its = its + abs(sh%newton_it)
       if (.not. conv) print *, 'WARNING: Newton not converged at step ', step, sh%newton_it, sh%newton_res
       call shell_accept(sh, f)
       dx = 0.0_dp;  dy = 0.0_dp
       do n = 1, sh%nn
          if (sh%X0(1,n) == L) then
             dx = dx + sh%u(1,n)/ntip;  dy = dy + sh%u(3,n)/ntip
          end if
       end do
       write(iu,'(3ES16.8)') tnow, dx, dy
       if (tnow > tend - 4.0_dp) then          ! statistics over the last 4 s
          dymin = min(dymin, dy);  dymax = max(dymax, dy);  dxmin = min(dxmin, dx);  dxmax = max(dxmax, dx)
       end if
       ! frequency from the upward crossings of the reference mean (-63.6 mm) over the last 6 s
       if (tnow > tend - 6.0_dp) then
          if (dy > -0.0636_dp .and. dy_prev() <= -0.0636_dp) then
             ncross = ncross + 1
             if (ncross == 1) tprev = tnow
             tlast = tnow
          end if
       end if
       call store_prev(dy)
    end do
    close(iu)
    call shell_energy(sh, ekin, estr)
    print '(A,F7.4,A,I0,A,F6.2)', ' CSM3 (40 x 2 shell, plane strain, rho_inf 0.8), dt = ', dt, ', steps ', nstep, &
       ', mean Newton its/step ', real(its)/nstep
    print '(A,F9.3,A,F9.3,A)', '   dy: mean ', 500.0_dp*(dymax + dymin), ' mm, amplitude ', 500.0_dp*(dymax - dymin), ' mm'
    print '(A,F9.3,A,F9.3,A)', '   dx: mean ', 500.0_dp*(dxmax + dxmin), ' mm, amplitude ', 500.0_dp*(dxmax - dxmin), ' mm'
    if (ncross > 1) print '(A,F9.4,A)', '   frequency ', (ncross - 1)/(tlast - tprev), ' Hz'
    print '(A)', '   Turek-Hron reference: dy -63.607 +- 65.160 mm, dx -14.305 +- 14.305 mm, f = 1.0995 Hz'
    print '(A,2ES12.4)', '   final kinetic / strain energy ', ekin, estr
  end subroutine test_csm3

  real(dp) function dy_prev()
    dy_prev = dyp
  end function dy_prev
  subroutine store_prev(d)
    real(dp), intent(in) :: d
    dyp = d
  end subroutine store_prev
  real(dp) function pi_()
    pi_ = 3.14159265358979323846_dp
  end function pi_
end program fsi_test
