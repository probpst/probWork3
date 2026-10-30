!=======================================================================
! fsi_shell.F90
!   STRUCTURAL SOLVER of the fluid-structure interaction (update 24):
!   geometrically nonlinear (large rotation, small strain) shell on a
!   quadrilateral mid-surface mesh, independent of the flow code (no MPI,
!   no global_vars) so that it can be verified on its own (fsi_test.F90).
!
!   ELEMENT: 4-node flat shell in a CO-ROTATIONAL frame.  The local
!   (small-strain) element is the sum of
!     * the bilinear plane-stress membrane (2 x 2 Gauss),
!     * the MITC4 Mindlin-Reissner plate (Bathe-Dvorkin assumed
!       transverse shear strains, 2 x 2 Gauss, no shear locking),
!     * the Hughes-Brezzi drilling term  gamma G t A (theta_z - 1/2(v,x - u,y))^2
!       (one point) that gives the rotation about the normal a stiffness
!       tied to the in-plane rotation of the displacement field.
!   The local stiffness K_l (24 x 24, dofs u v w thx thy thz per node) is
!   constant; the nonlinearity is in the kinematics: an element frame Ef
!   follows the current nodal positions (diagonal construction), the
!   local displacements are the current local coordinates minus the
!   initial ones, the local rotations the logarithm of Ef^T R_a E0 with
!   R_a the nodal triad (updated multiplicatively, R <- exp(dtheta) R).
!   Internal force f = T f_l, tangent K = T K_l T^T with T = diag(Ef):
!   the simplified co-rotational tangent (the frame-spin terms are
!   dropped), symmetric, exact in the limit of small increments - the
!   Newton iterations of a time step converge linearly, which is
!   sufficient with the hundreds of steps per structural period of an
!   FSI run.
!   MATERIAL: isotropic, per element (Ef, nu, rho, thickness t) so that
!   a part can carry several materials (Submat regions) and the
!   thickness measured from the fluid surface grid.
!   MASS: lumped, m_a = sum rho t A_e / 4, rotational inertia
!   m_a (t^2/12 + A_e/48).  DAMPING: Rayleigh, stiffness proportional,
!   a_K = 2 zeta / omega_1 (omega_1 from inverse iteration at set-up).
!   TIME: generalized-alpha (Chung-Hulbert, rho_inf = 0.8 by default)
!   with Newton iterations on the nonlinear residual; the step can be
!   repeated from the state at t_n with new loads (FSI sub-iterations,
!   shell_step) and committed (shell_accept).  Quasi-static solve:
!   shell_static.  Linear algebra: symmetric banded Cholesky after a
!   reverse Cuthill-McKee numbering of the nodes.
!=======================================================================
module fsi_shell
  implicit none
  private
  integer, parameter, public :: dp = kind(1.0d0)
  real(dp), parameter :: pi = 3.14159265358979323846_dp

  type, public :: shell_t
     ! ---- mesh / material (filled by the user before shell_setup) -------
     integer :: nn = 0, ne = 0
     real(dp), allocatable :: X0(:,:)          ! (3,nn) reference node coordinates
     integer,  allocatable :: conn(:,:)        ! (4,ne) counter-clockwise node numbers
     real(dp), allocatable :: te(:), Ee(:), nue(:), rhoe(:)   ! per element: thickness, Ef, nu, rho
     logical,  allocatable :: fixed(:,:)       ! (6,nn) .true. = dof held at 0
     real(dp) :: zeta = 0.0_dp                 ! damping ratio of the first mode
     real(dp) :: rho_inf = 0.8_dp              ! generalized-alpha spectral radius
     real(dp) :: gdrill = 0.1_dp               ! drilling penalty factor (Hughes-Brezzi)
     integer  :: max_newton = 20
     real(dp) :: tol_newton = 1.0e-6_dp
     ! ---- derived (shell_setup) ------------------------------------------
     real(dp), allocatable :: E0(:,:,:)        ! (3,3,ne) initial element frames (columns e1 e2 e3)
     real(dp), allocatable :: xl0(:,:,:)       ! (2,4,ne) initial local node coordinates
     real(dp), allocatable :: Ae(:)            ! element areas
     real(dp), allocatable :: Kl(:,:,:)        ! (24,24,ne) local stiffness
     real(dp), allocatable :: mass(:,:)        ! (6,nn) lumped mass
     real(dp), allocatable :: anode(:)         ! nodal area (stress averaging)
     integer,  allocatable :: perm(:)          ! perm(node) = position in the numbering
     integer  :: ndof = 0, kd = 0              ! dofs, half bandwidth
     real(dp), allocatable :: Ab(:,:)          ! (0:kd, ndof) banded lower factor
     real(dp) :: omega1 = 0.0_dp, aK = 0.0_dp, aM = 0.0_dp
     real(dp) :: tmin = 0.0_dp, tmax = 0.0_dp, mtot = 0.0_dp
     ! ---- state at t_n ------------------------------------------------------
     real(dp), allocatable :: u(:,:), v(:,:), a(:,:)     ! (6,nn): translations / angular velocity, acceleration
     real(dp), allocatable :: R(:,:,:)                   ! (3,3,nn) nodal triads
     real(dp), allocatable :: fext_n(:,:)                ! (6,nn) external load at t_n
     ! ---- trial state at t_n+1 ----------------------------------------------
     real(dp), allocatable :: u1(:,:), v1(:,:), a1(:,:), R1(:,:,:), dth(:,:)
     integer  :: newton_it = 0
     integer  :: substeps = 1                  ! sub-steps the last shell_step needed (24b fallback)
     real(dp) :: newton_res = 0.0_dp, newton_res0 = 0.0_dp
     logical  :: ready = .false.
  end type shell_t

  public :: shell_setup, shell_step, shell_accept, shell_static, shell_stress, shell_free
  public :: rot_exp, rot_log, shell_energy, shell_set_state_from_trial, shell_mode

contains

  !=====================================================================
  ! rotation vector -> matrix (Rodrigues) and back
  !=====================================================================
  pure function rot_exp(th) result(Rm)
    real(dp), intent(in) :: th(3)
    real(dp) :: Rm(3,3), ang, K(3,3), s, c1
    integer :: i
    ang = sqrt(sum(th*th))
    K = 0.0_dp
    K(1,2) = -th(3);  K(1,3) =  th(2)
    K(2,1) =  th(3);  K(2,3) = -th(1)
    K(3,1) = -th(2);  K(3,2) =  th(1)
    if (ang < 1.0e-8_dp) then
       s = 1.0_dp - ang*ang/6.0_dp;  c1 = 0.5_dp - ang*ang/24.0_dp
    else
       s = sin(ang)/ang;  c1 = (1.0_dp - cos(ang))/(ang*ang)
    end if
    Rm = s*K + c1*matmul(K, K)
    do i = 1, 3
       Rm(i,i) = Rm(i,i) + 1.0_dp
    end do
  end function rot_exp

  pure function rot_log(Rm) result(th)
    real(dp), intent(in) :: Rm(3,3)
    real(dp) :: th(3), c, ang, s, w(3)
    c = 0.5_dp*(Rm(1,1) + Rm(2,2) + Rm(3,3) - 1.0_dp)
    c = max(-1.0_dp, min(1.0_dp, c))
    ang = acos(c)
    w = 0.5_dp*[Rm(3,2) - Rm(2,3), Rm(1,3) - Rm(3,1), Rm(2,1) - Rm(1,2)]
    if (ang < 1.0e-6_dp) then
       th = w*(1.0_dp + ang*ang/6.0_dp)
    else if (ang < pi - 1.0e-4_dp) then
       s = sin(ang)
       th = w*ang/s
    else
       ! near pi: axis from the symmetric part
       th(1) = sqrt(max(0.0_dp, 0.5_dp*(Rm(1,1) + 1.0_dp)))
       th(2) = sqrt(max(0.0_dp, 0.5_dp*(Rm(2,2) + 1.0_dp)))
       th(3) = sqrt(max(0.0_dp, 0.5_dp*(Rm(3,3) + 1.0_dp)))
       if (w(1) < 0.0_dp) th(1) = -th(1)
       if (w(2) < 0.0_dp) th(2) = -th(2)
       if (w(3) < 0.0_dp) th(3) = -th(3)
       th = th/max(sqrt(sum(th*th)), 1.0e-30_dp)*ang
    end if
  end function rot_log

  pure function cross3(a, b) result(c)
    real(dp), intent(in) :: a(3), b(3)
    real(dp) :: c(3)
    c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
  end function cross3

  !---------------------------------------------------------------------
  ! element frame from four node positions: e1 along the xi diagonal
  ! direction, e3 the normal, origin at the centroid
  !---------------------------------------------------------------------
  pure subroutine elem_frame(x, Ef, c)
    real(dp), intent(in)  :: x(3,4)
    real(dp), intent(out) :: Ef(3,3), c(3)
    real(dp) :: v1(3), v2(3), e1(3), e2(3), e3(3)
    v1 = x(:,2) + x(:,3) - x(:,1) - x(:,4)
    v2 = x(:,3) + x(:,4) - x(:,1) - x(:,2)
    e3 = cross3(v1, v2);  e3 = e3/max(sqrt(sum(e3*e3)), 1.0e-300_dp)
    e1 = v1 - sum(v1*e3)*e3;  e1 = e1/max(sqrt(sum(e1*e1)), 1.0e-300_dp)
    e2 = cross3(e3, e1)
    Ef(:,1) = e1;  Ef(:,2) = e2;  Ef(:,3) = e3
    c = 0.25_dp*(x(:,1) + x(:,2) + x(:,3) + x(:,4))
  end subroutine elem_frame

  !---------------------------------------------------------------------
  ! bilinear shape functions and derivatives at (r,s)
  !---------------------------------------------------------------------
  pure subroutine shape4(r, s, N, dNr, dNs)
    real(dp), intent(in)  :: r, s
    real(dp), intent(out) :: N(4), dNr(4), dNs(4)
    N   = 0.25_dp*[(1-r)*(1-s), (1+r)*(1-s), (1+r)*(1+s), (1-r)*(1+s)]
    dNr = 0.25_dp*[-(1-s), (1-s), (1+s), -(1+s)]
    dNs = 0.25_dp*[-(1-r), -(1+r), (1+r), (1-r)]
  end subroutine shape4

  !---------------------------------------------------------------------
  ! Local stiffness of element e (local coordinates xl(2,4), thickness t,
  ! Ef, nu): membrane + MITC4 plate + drilling.  Dof order per node
  ! (u, v, w, thx, thy, thz), node a at 6(a-1)+1 ... 6a.
  ! Plate kinematics: u_z-plane displacement (z thy, -z thx, w),
  !   curvatures kx = thy,x, ky = -thx,y, kxy = thy,y - thx,x
  !   shear gxz = w,x + thy, gyz = w,y - thx.
  !---------------------------------------------------------------------
  subroutine local_stiffness(xl, t, Em, nu, gdr, K, area)
    real(dp), intent(in)  :: xl(2,4), t, Em, nu, gdr
    real(dp), intent(out) :: K(24,24), area
    real(dp) :: gp(2), N(4), dNr(4), dNs(4), J(2,2), Ji(2,2), detJ, dNx(4), dNy(4)
    real(dp) :: Bm(3,24), Bb(3,24), Bs(2,24), Dm(3,3), Db(3,3), Ds(2,2), G, r, s, w
    real(dp) :: gtA(24), gtB(24), gtC(24), gtD(24), gr(24), gs(24), Jt(2,2), bdr(24), kdr
    integer :: ig, jg, a, c
    real(dp), parameter :: gpt(2) = [-0.577350269189626_dp, 0.577350269189626_dp]
    K = 0.0_dp;  area = 0.0_dp
    G = Em/(2.0_dp*(1.0_dp + nu))
    Dm = reshape([1.0_dp, nu, 0.0_dp,  nu, 1.0_dp, 0.0_dp,  0.0_dp, 0.0_dp, 0.5_dp*(1.0_dp - nu)], [3,3])
    Db = Dm*Em*t**3/(12.0_dp*(1.0_dp - nu*nu))
    Dm = Dm*Em*t/(1.0_dp - nu*nu)
    Ds = 0.0_dp;  Ds(1,1) = 5.0_dp/6.0_dp*G*t;  Ds(2,2) = Ds(1,1)
    ! MITC4 tying points: A = (0,1), C = (0,-1) for g_rz;  D = (1,0), B = (-1,0) for g_sz
    call shear_cov(xl, 0.0_dp,  1.0_dp, 1, gtA)
    call shear_cov(xl, 0.0_dp, -1.0_dp, 1, gtC)
    call shear_cov(xl,  1.0_dp, 0.0_dp, 2, gtD)
    call shear_cov(xl, -1.0_dp, 0.0_dp, 2, gtB)
    do jg = 1, 2
     do ig = 1, 2
        r = gpt(ig);  s = gpt(jg)
        call shape4(r, s, N, dNr, dNs)
        J(1,1) = sum(dNr*xl(1,:));  J(1,2) = sum(dNr*xl(2,:))
        J(2,1) = sum(dNs*xl(1,:));  J(2,2) = sum(dNs*xl(2,:))
        detJ = J(1,1)*J(2,2) - J(1,2)*J(2,1)
        Ji(1,1) =  J(2,2)/detJ;  Ji(1,2) = -J(1,2)/detJ
        Ji(2,1) = -J(2,1)/detJ;  Ji(2,2) =  J(1,1)/detJ
        dNx = Ji(1,1)*dNr + Ji(1,2)*dNs
        dNy = Ji(2,1)*dNr + Ji(2,2)*dNs
        w = detJ;  area = area + w
        Bm = 0.0_dp;  Bb = 0.0_dp
        do a = 1, 4
           c = 6*(a-1)
           Bm(1, c+1) = dNx(a);  Bm(2, c+2) = dNy(a)
           Bm(3, c+1) = dNy(a);  Bm(3, c+2) = dNx(a)
           Bb(1, c+5) =  dNx(a)
           Bb(2, c+4) = -dNy(a)
           Bb(3, c+5) =  dNy(a);  Bb(3, c+4) = -dNx(a)
        end do
        ! assumed covariant shear strains, transformed to Cartesian: (gx,gy) = J^-1 (gr,gs)
        gr = 0.5_dp*(1.0_dp + s)*gtA + 0.5_dp*(1.0_dp - s)*gtC
        gs = 0.5_dp*(1.0_dp + r)*gtD + 0.5_dp*(1.0_dp - r)*gtB
        Jt = Ji          ! (gx,gy)^T = Ji (gr,gs)^T  with Ji = J^-1 as above (J rows = covariant base vectors)
        Bs(1,:) = Jt(1,1)*gr + Jt(1,2)*gs
        Bs(2,:) = Jt(2,1)*gr + Jt(2,2)*gs
        K = K + w*(matmul(transpose(Bm), matmul(Dm, Bm)) + matmul(transpose(Bb), matmul(Db, Bb)) &
                 + matmul(transpose(Bs), matmul(Ds, Bs)))
     end do
    end do
    ! drilling: gdr G t (thz - 1/2 (v,x - u,y))^2 integrated with 2 x 2 points (a one-point
    ! evaluation leaves three zero-energy drilling modes per element)
    do jg = 1, 2
     do ig = 1, 2
        r = gpt(ig);  s = gpt(jg)
        call shape4(r, s, N, dNr, dNs)
        J(1,1) = sum(dNr*xl(1,:));  J(1,2) = sum(dNr*xl(2,:))
        J(2,1) = sum(dNs*xl(1,:));  J(2,2) = sum(dNs*xl(2,:))
        detJ = J(1,1)*J(2,2) - J(1,2)*J(2,1)
        Ji(1,1) =  J(2,2)/detJ;  Ji(1,2) = -J(1,2)/detJ
        Ji(2,1) = -J(2,1)/detJ;  Ji(2,2) =  J(1,1)/detJ
        dNx = Ji(1,1)*dNr + Ji(1,2)*dNs
        dNy = Ji(2,1)*dNr + Ji(2,2)*dNs
        bdr = 0.0_dp
        do a = 1, 4
           c = 6*(a-1)
           bdr(c+6) = N(a)
           bdr(c+2) = -0.5_dp*dNx(a)
           bdr(c+1) =  0.5_dp*dNy(a)
        end do
        kdr = gdr*G*t*detJ
        do a = 1, 24
           K(:,a) = K(:,a) + kdr*bdr*bdr(a)
        end do
     end do
    end do
    K = 0.5_dp*(K + transpose(K))
  contains
    ! covariant transverse shear strain g_rz (dir = 1) or g_sz (dir = 2) at (r,s) as a row over the 24 dofs:
    !   g_rz = w,r + x,r thy - y,r thx
    subroutine shear_cov(xl_, r_, s_, dir, g)
      real(dp), intent(in)  :: xl_(2,4), r_, s_
      integer,  intent(in)  :: dir
      real(dp), intent(out) :: g(24)
      real(dp) :: N_(4), dNr_(4), dNs_(4), xr, yr, dN(4)
      integer :: a_, c_
      call shape4(r_, s_, N_, dNr_, dNs_)
      if (dir == 1) then
         dN = dNr_
      else
         dN = dNs_
      end if
      xr = sum(dN*xl_(1,:));  yr = sum(dN*xl_(2,:))
      g = 0.0_dp
      do a_ = 1, 4
         c_ = 6*(a_-1)
         g(c_+3) = dN(a_)
         g(c_+5) =  xr*N_(a_)
         g(c_+4) = -yr*N_(a_)
      end do
    end subroutine shear_cov
  end subroutine local_stiffness

  !=====================================================================
  ! Set-up: frames, local stiffness, lumped mass, numbering, band,
  ! first eigenfrequency, damping coefficient, zero state.
  !=====================================================================
  subroutine shell_setup(sh, info)
    type(shell_t), intent(inout) :: sh
    integer, intent(out), optional :: info
    integer :: e, a, n, i, c, maxb
    real(dp) :: x(3,4), Ef(3,3), cc(3), K(24,24), area, xl(2,4), ma
    integer, allocatable :: deg(:)
    n = sh%nn
    if (allocated(sh%E0)) deallocate(sh%E0, sh%xl0, sh%Ae, sh%Kl, sh%mass, sh%anode)
    allocate(sh%E0(3,3,sh%ne), sh%xl0(2,4,sh%ne), sh%Ae(sh%ne), sh%Kl(24,24,sh%ne), sh%mass(6,n), sh%anode(n))
    sh%mass = 0.0_dp;  sh%anode = 0.0_dp;  sh%mtot = 0.0_dp
    sh%tmin = huge(1.0_dp);  sh%tmax = 0.0_dp
    do e = 1, sh%ne
       do a = 1, 4
          x(:,a) = sh%X0(:, sh%conn(a,e))
       end do
       call elem_frame(x, Ef, cc)
       sh%E0(:,:,e) = Ef
       do a = 1, 4
          xl(1,a) = sum(Ef(:,1)*(x(:,a) - cc));  xl(2,a) = sum(Ef(:,2)*(x(:,a) - cc))
       end do
       sh%xl0(:,:,e) = xl
       call local_stiffness(xl, sh%te(e), sh%Ee(e), sh%nue(e), sh%gdrill, K, area)
       sh%Kl(:,:,e) = K;  sh%Ae(e) = area
       sh%tmin = min(sh%tmin, sh%te(e));  sh%tmax = max(sh%tmax, sh%te(e))
       ma = sh%rhoe(e)*sh%te(e)*area/4.0_dp
       sh%mtot = sh%mtot + 4.0_dp*ma
       do a = 1, 4
          c = sh%conn(a,e)
          sh%mass(1:3,c) = sh%mass(1:3,c) + ma
          sh%mass(4:6,c) = sh%mass(4:6,c) + ma*(sh%te(e)**2/12.0_dp + area/48.0_dp)
          sh%anode(c) = sh%anode(c) + area/4.0_dp
       end do
    end do
    ! numbering (reverse Cuthill-McKee on the node graph) and bandwidth
    call rcm_numbering(sh)
    maxb = 0
    do e = 1, sh%ne
       do a = 1, 4
          do i = 1, 4
             maxb = max(maxb, abs(sh%perm(sh%conn(a,e)) - sh%perm(sh%conn(i,e))))
          end do
       end do
    end do
    sh%ndof = 6*n;  sh%kd = 6*maxb + 5
    if (allocated(sh%Ab)) deallocate(sh%Ab)
    allocate(sh%Ab(0:sh%kd, sh%ndof))
    ! state
    if (allocated(sh%u)) deallocate(sh%u, sh%v, sh%a, sh%R, sh%fext_n, sh%u1, sh%v1, sh%a1, sh%R1, sh%dth)
    allocate(sh%u(6,n), sh%v(6,n), sh%a(6,n), sh%R(3,3,n), sh%fext_n(6,n), sh%u1(6,n), sh%v1(6,n), sh%a1(6,n), &
             sh%R1(3,3,n), sh%dth(3,n))
    sh%u = 0.0_dp;  sh%v = 0.0_dp;  sh%a = 0.0_dp;  sh%fext_n = 0.0_dp
    sh%u1 = 0.0_dp;  sh%v1 = 0.0_dp;  sh%a1 = 0.0_dp;  sh%dth = 0.0_dp
    do i = 1, n
       sh%R(:,:,i) = reshape([1,0,0, 0,1,0, 0,0,1], [3,3])*1.0_dp
       sh%R1(:,:,i) = sh%R(:,:,i)
    end do
    sh%ready = .true.
    ! first natural frequency (inverse iteration) and the damping coefficient
    call first_frequency(sh)
    sh%aK = 0.0_dp;  sh%aM = 0.0_dp
    if (sh%zeta > 0.0_dp .and. sh%omega1 > 0.0_dp) sh%aK = 2.0_dp*sh%zeta/sh%omega1
    if (present(info)) info = 0
    allocate(deg(1));  deallocate(deg)
  end subroutine shell_setup

  !---------------------------------------------------------------------
  ! reverse Cuthill-McKee node numbering (perm(node) = new position)
  !---------------------------------------------------------------------
  subroutine rcm_numbering(sh)
    type(shell_t), intent(inout) :: sh
    integer :: n, e, a, b, i, j, k, cnt, start, head, tail, mind, nb, p, q, tmp
    integer, allocatable :: deg(:), ptr(:), adj(:), order(:), queue(:), nbr(:)
    logical, allocatable :: seen(:)
    n = sh%nn
    allocate(deg(n), ptr(n+1), seen(n), order(n), queue(n))
    deg = 0
    do e = 1, sh%ne
       do a = 1, 4
          do b = 1, 4
             if (a /= b) deg(sh%conn(a,e)) = deg(sh%conn(a,e)) + 1
          end do
       end do
    end do
    ptr(1) = 1
    do i = 1, n
       ptr(i+1) = ptr(i) + deg(i)
    end do
    allocate(adj(ptr(n+1)))
    deg = 0
    do e = 1, sh%ne
       do a = 1, 4
          do b = 1, 4
             if (a == b) cycle
             i = sh%conn(a,e)
             adj(ptr(i) + deg(i)) = sh%conn(b,e);  deg(i) = deg(i) + 1
          end do
       end do
    end do
    ! BFS from a node of minimum degree (per connected component), neighbours by increasing degree
    seen = .false.;  cnt = 0
    allocate(nbr(maxval(deg) + 1))
    do while (cnt < n)
       start = 0;  mind = huge(1)
       do i = 1, n
          if (.not. seen(i) .and. deg(i) < mind) then
             mind = deg(i);  start = i
          end if
       end do
       head = 1;  tail = 1;  queue(1) = start;  seen(start) = .true.
       do while (head <= tail)
          i = queue(head);  head = head + 1
          cnt = cnt + 1;  order(cnt) = i
          nb = 0
          do k = ptr(i), ptr(i+1) - 1
             j = adj(k)
             if (seen(j)) cycle
             seen(j) = .true.;  nb = nb + 1;  nbr(nb) = j
          end do
          ! sort the new neighbours by degree
          do p = 2, nb
             tmp = nbr(p);  q = p - 1
             do while (q >= 1)
                if (deg(nbr(q)) <= deg(tmp)) exit
                nbr(q+1) = nbr(q);  q = q - 1
             end do
             nbr(q+1) = tmp
          end do
          do p = 1, nb
             tail = tail + 1;  queue(tail) = nbr(p)
          end do
       end do
    end do
    if (allocated(sh%perm)) deallocate(sh%perm)
    allocate(sh%perm(n))
    do i = 1, n
       sh%perm(order(i)) = n + 1 - i          ! reverse
    end do
    deallocate(deg, ptr, adj, seen, order, queue, nbr)
  end subroutine rcm_numbering

  pure integer function gdof(sh, node, c)
    type(shell_t), intent(in) :: sh
    integer, intent(in) :: node, c
    gdof = 6*(sh%perm(node) - 1) + c
  end function gdof

  !---------------------------------------------------------------------
  ! Banded symmetric assembly helpers (lower band storage A(r,j) = K(j+r,j))
  !---------------------------------------------------------------------
  subroutine band_add(sh, i, j, val)
    type(shell_t), intent(inout) :: sh
    integer,  intent(in) :: i, j
    real(dp), intent(in) :: val
    if (i >= j) then
       sh%Ab(i-j, j) = sh%Ab(i-j, j) + val
    end if
  end subroutine band_add

  subroutine band_cholesky(sh, ok)
    type(shell_t), intent(inout) :: sh
    logical, intent(out) :: ok
    integer :: j, k, i, n, kd, m
    real(dp) :: ajj
    n = sh%ndof;  kd = sh%kd;  ok = .true.
    do j = 1, n
       ajj = sh%Ab(0,j)
       if (ajj <= 0.0_dp) then
          ok = .false.;  return
       end if
       ajj = sqrt(ajj);  sh%Ab(0,j) = ajj
       m = min(kd, n - j)
       if (m > 0) then
          sh%Ab(1:m, j) = sh%Ab(1:m, j)/ajj
          do k = 1, m
             do i = k, m
                sh%Ab(i-k, j+k) = sh%Ab(i-k, j+k) - sh%Ab(i,j)*sh%Ab(k,j)
             end do
          end do
       end if
    end do
  end subroutine band_cholesky

  subroutine band_solve(sh, b)
    type(shell_t), intent(in) :: sh
    real(dp), intent(inout) :: b(:)
    integer :: j, n, kd, m, i
    n = sh%ndof;  kd = sh%kd
    do j = 1, n                              ! L y = b
       b(j) = b(j)/sh%Ab(0,j)
       m = min(kd, n - j)
       do i = 1, m
          b(j+i) = b(j+i) - sh%Ab(i,j)*b(j)
       end do
    end do
    do j = n, 1, -1                          ! L^T x = y
       m = min(kd, n - j)
       do i = 1, m
          b(j) = b(j) - sh%Ab(i,j)*b(j+i)
       end do
       b(j) = b(j)/sh%Ab(0,j)
    end do
  end subroutine band_solve

  !---------------------------------------------------------------------
  ! Element kinematics at a configuration (x = X0 + ut, triads Rt):
  ! frame Ef, local displacement vector p_l(24)
  !---------------------------------------------------------------------
  subroutine elem_local(sh, e, ut, Rt, Ef, pl)
    type(shell_t), intent(in) :: sh
    integer,  intent(in)  :: e
    real(dp), intent(in)  :: ut(:,:), Rt(:,:,:)
    real(dp), intent(out) :: Ef(3,3), pl(24)
    real(dp) :: x(3,4), cc(3), y(3), Rl(3,3)
    integer :: a, c, nd
    do a = 1, 4
       nd = sh%conn(a,e)
       x(:,a) = sh%X0(:,nd) + ut(1:3,nd)
    end do
    call elem_frame(x, Ef, cc)
    do a = 1, 4
       nd = sh%conn(a,e);  c = 6*(a-1)
       y = matmul(transpose(Ef), x(:,a) - cc)
       pl(c+1) = y(1) - sh%xl0(1,a,e);  pl(c+2) = y(2) - sh%xl0(2,a,e);  pl(c+3) = y(3)
       Rl = matmul(transpose(Ef), matmul(Rt(:,:,nd), sh%E0(:,:,e)))
       pl(c+4:c+6) = rot_log(Rl)
    end do
  end subroutine elem_local

  !---------------------------------------------------------------------
  ! internal force (6,nn) at a configuration; optionally assemble the
  ! tangent (scaled by fac) into the band matrix
  !---------------------------------------------------------------------
  subroutine internal_force(sh, ut, Rt, fint, assemble, fac)
    type(shell_t), intent(inout) :: sh
    real(dp), intent(in)  :: ut(:,:), Rt(:,:,:)
    real(dp), intent(out) :: fint(:,:)
    logical,  intent(in)  :: assemble
    real(dp), intent(in)  :: fac
    real(dp) :: Ef(3,3), pl(24), fl(24), T(24,24), Kg(24,24)
    integer :: e, a, b, c, nd, i, j, gi, gj
    fint = 0.0_dp
    do e = 1, sh%ne
       call elem_local(sh, e, ut, Rt, Ef, pl)
       fl = matmul(sh%Kl(:,:,e), pl)
       do a = 1, 4
          nd = sh%conn(a,e);  c = 6*(a-1)
          fint(1:3,nd) = fint(1:3,nd) + matmul(Ef, fl(c+1:c+3))
          fint(4:6,nd) = fint(4:6,nd) + matmul(Ef, fl(c+4:c+6))
       end do
       if (assemble) then
          T = 0.0_dp
          do a = 1, 4
             c = 6*(a-1)
             T(c+1:c+3, c+1:c+3) = Ef;  T(c+4:c+6, c+4:c+6) = Ef
          end do
          Kg = matmul(T, matmul(sh%Kl(:,:,e), transpose(T)))
          do a = 1, 4
             do i = 1, 6
                gi = gdof(sh, sh%conn(a,e), i)
                do b = 1, 4
                   do j = 1, 6
                      gj = gdof(sh, sh%conn(b,e), j)
                      if (gi >= gj) call band_add(sh, gi, gj, fac*Kg(6*(a-1)+i, 6*(b-1)+j))
                   end do
                end do
             end do
          end do
       end if
    end do
  end subroutine internal_force

  !=====================================================================
  ! One generalized-alpha step from the state at t_n to t_n+1 with the
  ! external load fext1 (6,nn: forces and moments at t_n+1); repeatable:
  ! every call starts again from (u, v, a, R) at t_n.  The result is the
  ! trial state (u1, v1, a1, R1); shell_accept commits it.
  !=====================================================================
  subroutine shell_step(sh, dt, fext1, converged)
    type(shell_t), intent(inout) :: sh
    real(dp), intent(in) :: dt
    real(dp), intent(in) :: fext1(:,:)
    logical, intent(out) :: converged
    real(dp), allocatable :: u0(:,:), v0(:,:), a0(:,:), R0(:,:,:), f0(:,:), fj(:,:)
    integer :: ns, j, i, nit
    logical :: ok
    sh%substeps = 1
    call step_once(sh, dt, fext1, converged)
    if (converged) return
    ! fallback (24b): the step in 2, 4, 8 sub-steps with the load interpolated linearly from t_n -
    ! the Newton iterations of a whole step can fail when the load changes abruptly (a sub-iteration
    ! of the coupling far from the converged interface position); the result is again a trial state
    ! relative to the state at t_n (shell_accept commits it)
    allocate(u0, source=sh%u);  allocate(v0, source=sh%v);  allocate(a0, source=sh%a)
    allocate(R0, source=sh%R);  allocate(f0, source=sh%fext_n);  allocate(fj, mold=fext1)
    nit = 0
    do ns = 2, 8, 2
       if (ns == 6) cycle
       ok = .true.
       do j = 1, ns
          fj = f0 + (fext1 - f0)*(real(j, dp)/ns)
          call step_once(sh, dt/ns, fj, converged)
          nit = nit + abs(sh%newton_it)
          if (.not. converged) then
             ok = .false.;  exit
          end if
          call shell_accept(sh, fj)
       end do
       if (ok) then
          sh%u1 = sh%u;  sh%v1 = sh%v;  sh%a1 = sh%a;  sh%R1 = sh%R
          do i = 1, sh%nn
             sh%dth(:, i) = rot_log(matmul(sh%R1(:,:,i), transpose(R0(:,:,i))))
          end do
          sh%u = u0;  sh%v = v0;  sh%a = a0;  sh%R = R0;  sh%fext_n = f0
          sh%substeps = ns;  sh%newton_it = nit;  converged = .true.
          return
       end if
       sh%u = u0;  sh%v = v0;  sh%a = a0;  sh%R = R0;  sh%fext_n = f0
    end do
    ! no success: the trial state of the whole step (not converged) is returned
    sh%u = u0;  sh%v = v0;  sh%a = a0;  sh%R = R0;  sh%fext_n = f0
    call step_once(sh, dt, fext1, converged)
    sh%substeps = 1
  end subroutine shell_step

  subroutine step_once(sh, dt, fext1, converged)
    type(shell_t), intent(inout) :: sh
    real(dp), intent(in) :: dt
    real(dp), intent(in) :: fext1(:,:)
    logical, intent(out) :: converged
    real(dp) :: am, af, beta, gam, c0, c1, c2
    real(dp), allocatable :: ua(:,:), va(:,:), aa(:,:), Ra(:,:,:), fint(:,:), res(:), fa(:,:), dv(:,:), ddth(:,:)
    integer :: n, it, i, c, g
    real(dp) :: rnorm, dnorm, unorm, dtha(3)
    logical :: ok
    n = sh%nn
    allocate(ua(6,n), va(6,n), aa(6,n), Ra(3,3,n), fint(6,n), res(sh%ndof), fa(6,n), dv(6,n), ddth(3,n))
    am = (2.0_dp*sh%rho_inf - 1.0_dp)/(sh%rho_inf + 1.0_dp)
    af = sh%rho_inf/(sh%rho_inf + 1.0_dp)
    beta = 0.25_dp*(1.0_dp - am + af)**2
    gam = 0.5_dp - am + af
    c0 = (1.0_dp - am)/(beta*dt*dt);  c1 = (1.0_dp - af)*gam/(beta*dt);  c2 = 1.0_dp - af
    ! predictor: constant acceleration (a1 = a_n)
    sh%a1 = sh%a
    sh%u1(1:3,:) = sh%u(1:3,:) + dt*sh%v(1:3,:) + dt*dt*(0.5_dp*sh%a(1:3,:))
    sh%dth = dt*sh%v(4:6,:) + dt*dt*0.5_dp*sh%a(4:6,:)
    do i = 1, n
       do c = 1, 6
          if (sh%fixed(c,i)) then
             if (c <= 3) then
                sh%u1(c,i) = 0.0_dp
             else
                sh%dth(c-3,i) = 0.0_dp
             end if
          end if
       end do
    end do
    converged = .false.
    do it = 1, sh%max_newton
       ! accelerations / velocities of the trial state (Newmark)
       sh%a1(1:3,:) = (sh%u1(1:3,:) - sh%u(1:3,:) - dt*sh%v(1:3,:) - dt*dt*(0.5_dp - beta)*sh%a(1:3,:))/(beta*dt*dt)
       sh%a1(4:6,:) = (sh%dth - dt*sh%v(4:6,:) - dt*dt*(0.5_dp - beta)*sh%a(4:6,:))/(beta*dt*dt)
       sh%v1 = sh%v + dt*((1.0_dp - gam)*sh%a + gam*sh%a1)
       do i = 1, n
          sh%R1(:,:,i) = matmul(rot_exp(sh%dth(:,i)), sh%R(:,:,i))
       end do
       ! alpha-weighted configuration
       ua = sh%u + (1.0_dp - af)*(sh%u1 - sh%u)
       va = sh%v + (1.0_dp - af)*(sh%v1 - sh%v)
       aa = sh%a + (1.0_dp - am)*(sh%a1 - sh%a)
       do i = 1, n
          dtha = (1.0_dp - af)*sh%dth(:,i)
          Ra(:,:,i) = matmul(rot_exp(dtha), sh%R(:,:,i))
       end do
       fa = sh%fext_n + (1.0_dp - af)*(fext1 - sh%fext_n)
       ! tangent: c0 M + c1 C + c2 K  (C = aK K + aM M)
       sh%Ab = 0.0_dp
       call internal_force(sh, ua, Ra, fint, .true., c2 + c1*sh%aK)
       do i = 1, n
          do c = 1, 6
             g = gdof(sh, i, c)
             sh%Ab(0,g) = sh%Ab(0,g) + (c0 + c1*sh%aM)*sh%mass(c,i)
          end do
       end do
       ! residual r = fa - M aa - C va - fint;  C va = aK K va + aM M va: the stiffness
       ! damping force is approximated by aK * fint-increment rate, i.e. aK K_T va (tangent)
       call damping_force(sh, ua, Ra, va, dv)
       res = 0.0_dp;  rnorm = 0.0_dp
       do i = 1, n
          do c = 1, 6
             g = gdof(sh, i, c)
             if (sh%fixed(c,i)) then
                res(g) = 0.0_dp
                sh%Ab(0,g) = sh%Ab(0,g) + 1.0e30_dp
             else
                res(g) = fa(c,i) - sh%mass(c,i)*aa(c,i) - sh%aM*sh%mass(c,i)*va(c,i) - dv(c,i) - fint(c,i)
                rnorm = max(rnorm, abs(res(g)))
             end if
          end do
       end do
       if (it == 1) sh%newton_res0 = rnorm
       call band_cholesky(sh, ok)
       if (.not. ok) then
          sh%newton_it = -it;  sh%newton_res = rnorm
          return
       end if
       call band_solve(sh, res)
       dnorm = 0.0_dp;  unorm = 0.0_dp
       do i = 1, n
          do c = 1, 3
             g = gdof(sh, i, c)
             sh%u1(c,i) = sh%u1(c,i) + res(g)
             dnorm = max(dnorm, abs(res(g)));  unorm = max(unorm, abs(sh%u1(c,i) - sh%u(c,i)))
          end do
          do c = 4, 6
             g = gdof(sh, i, c)
             sh%dth(c-3,i) = sh%dth(c-3,i) + res(g)
             dnorm = max(dnorm, abs(res(g)));  unorm = max(unorm, abs(sh%dth(c-3,i)))
          end do
       end do
       sh%newton_it = it;  sh%newton_res = rnorm
       if (dnorm <= sh%tol_newton*max(unorm, 1.0e-12_dp) .or. dnorm < 1.0e-14_dp) then
          converged = .true.;  exit
       end if
    end do
    ! final trial velocities / accelerations / triads for the converged displacements
    sh%a1(1:3,:) = (sh%u1(1:3,:) - sh%u(1:3,:) - dt*sh%v(1:3,:) - dt*dt*(0.5_dp - beta)*sh%a(1:3,:))/(beta*dt*dt)
    sh%a1(4:6,:) = (sh%dth - dt*sh%v(4:6,:) - dt*dt*(0.5_dp - beta)*sh%a(4:6,:))/(beta*dt*dt)
    sh%v1 = sh%v + dt*((1.0_dp - gam)*sh%a + gam*sh%a1)
    do i = 1, n
       sh%R1(:,:,i) = matmul(rot_exp(sh%dth(:,i)), sh%R(:,:,i))
       sh%u1(4:6,i) = sh%u(4:6,i) + sh%dth(:,i)
    end do
    deallocate(ua, va, aa, Ra, fint, res, fa, dv, ddth)
  end subroutine step_once

  ! stiffness-proportional damping force aK K_T v at a configuration (tangent times velocity,
  ! element by element in the co-rotational frame)
  subroutine damping_force(sh, ut, Rt, vt, fd)
    type(shell_t), intent(inout) :: sh
    real(dp), intent(in)  :: ut(:,:), Rt(:,:,:), vt(:,:)
    real(dp), intent(out) :: fd(:,:)
    real(dp) :: Ef(3,3), pl(24), vl(24), fl(24)
    integer :: e, a, c, nd
    fd = 0.0_dp
    if (sh%aK == 0.0_dp) return
    do e = 1, sh%ne
       call elem_local(sh, e, ut, Rt, Ef, pl)
       do a = 1, 4
          nd = sh%conn(a,e);  c = 6*(a-1)
          vl(c+1:c+3) = matmul(transpose(Ef), vt(1:3,nd))
          vl(c+4:c+6) = matmul(transpose(Ef), vt(4:6,nd))
       end do
       fl = sh%aK*matmul(sh%Kl(:,:,e), vl)
       do a = 1, 4
          nd = sh%conn(a,e);  c = 6*(a-1)
          fd(1:3,nd) = fd(1:3,nd) + matmul(Ef, fl(c+1:c+3))
          fd(4:6,nd) = fd(4:6,nd) + matmul(Ef, fl(c+4:c+6))
       end do
    end do
  end subroutine damping_force

  ! commit the trial state
  subroutine shell_accept(sh, fext1)
    type(shell_t), intent(inout) :: sh
    real(dp), intent(in) :: fext1(:,:)
    sh%u = sh%u1;  sh%v = sh%v1;  sh%a = sh%a1;  sh%R = sh%R1
    sh%fext_n = fext1
    sh%dth = 0.0_dp
  end subroutine shell_accept

  ! trial state <- state (e.g. after a restart) so that the outputs of the trial are valid
  subroutine shell_set_state_from_trial(sh)
    type(shell_t), intent(inout) :: sh
    sh%u1 = sh%u;  sh%v1 = sh%v;  sh%a1 = sh%a;  sh%R1 = sh%R;  sh%dth = 0.0_dp
  end subroutine shell_set_state_from_trial

  !=====================================================================
  ! Quasi-static equilibrium with the load fext (Newton from the current
  ! trial state); the result is in the trial state (velocities 0).
  !=====================================================================
  subroutine shell_static(sh, fext, converged)
    type(shell_t), intent(inout) :: sh
    real(dp), intent(in) :: fext(:,:)
    logical, intent(out) :: converged
    real(dp), allocatable :: fint(:,:), res(:)
    integer :: n, it, i, c, g
    real(dp) :: rnorm, dnorm, unorm
    logical :: ok
    n = sh%nn
    allocate(fint(6,n), res(sh%ndof))
    converged = .false.
    do i = 1, n
       sh%R1(:,:,i) = matmul(rot_exp(sh%dth(:,i)), sh%R(:,:,i))
    end do
    do it = 1, sh%max_newton
       sh%Ab = 0.0_dp
       call internal_force(sh, sh%u1, sh%R1, fint, .true., 1.0_dp)
       res = 0.0_dp;  rnorm = 0.0_dp
       do i = 1, n
          do c = 1, 6
             g = gdof(sh, i, c)
             if (sh%fixed(c,i)) then
                res(g) = 0.0_dp;  sh%Ab(0,g) = sh%Ab(0,g) + 1.0e30_dp
             else
                res(g) = fext(c,i) - fint(c,i);  rnorm = max(rnorm, abs(res(g)))
             end if
          end do
       end do
       if (it == 1) sh%newton_res0 = rnorm
       call band_cholesky(sh, ok)
       if (.not. ok) then
          sh%newton_it = -it;  sh%newton_res = rnorm;  return
       end if
       call band_solve(sh, res)
       ! increment limiter (the simplified tangent is only exact for small increments): at most
       ! 0.15 rad of rotation per iteration
       dnorm = 0.0_dp
       do i = 1, n
          do c = 4, 6
             dnorm = max(dnorm, abs(res(gdof(sh, i, c))))
          end do
       end do
       if (dnorm > 0.15_dp) res = res*(0.15_dp/dnorm)
       dnorm = 0.0_dp;  unorm = 0.0_dp
       do i = 1, n
          do c = 1, 3
             g = gdof(sh, i, c);  sh%u1(c,i) = sh%u1(c,i) + res(g)
             dnorm = max(dnorm, abs(res(g)));  unorm = max(unorm, abs(sh%u1(c,i)))
          end do
          do c = 4, 6
             g = gdof(sh, i, c);  sh%dth(c-3,i) = sh%dth(c-3,i) + res(g)
             dnorm = max(dnorm, abs(res(g)));  unorm = max(unorm, abs(sh%dth(c-3,i)))
          end do
          sh%R1(:,:,i) = matmul(rot_exp(sh%dth(:,i)), sh%R(:,:,i))
       end do
       sh%newton_it = it;  sh%newton_res = rnorm
       if (dnorm <= sh%tol_newton*max(unorm, 1.0e-12_dp) .or. dnorm < 1.0e-14_dp) then
          converged = .true.;  exit
       end if
    end do
    sh%v1 = 0.0_dp;  sh%a1 = 0.0_dp
    do i = 1, n
       sh%u1(4:6,i) = sh%u(4:6,i) + sh%dth(:,i)
    end do
    deallocate(fint, res)
  end subroutine shell_static

  !=====================================================================
  ! von Mises stress at the nodes (max of the two surfaces z = +-t/2) of
  ! the trial configuration, area-averaged from the element centres
  !=====================================================================
  subroutine shell_stress(sh, vm)
    type(shell_t), intent(inout) :: sh
    real(dp), intent(out) :: vm(:)
    real(dp) :: Ef(3,3), pl(24), N(4), dNr(4), dNs(4), J(2,2), Ji(2,2), detJ, dNx(4), dNy(4)
    real(dp) :: eps(3), kap(3), sg(3), s1, s2, z, Em, nu, fac, t
    integer :: e, a, c, nd, side
    vm = 0.0_dp
    do e = 1, sh%ne
       call elem_local(sh, e, sh%u1, sh%R1, Ef, pl)
       call shape4(0.0_dp, 0.0_dp, N, dNr, dNs)
       J(1,1) = sum(dNr*sh%xl0(1,:,e));  J(1,2) = sum(dNr*sh%xl0(2,:,e))
       J(2,1) = sum(dNs*sh%xl0(1,:,e));  J(2,2) = sum(dNs*sh%xl0(2,:,e))
       detJ = J(1,1)*J(2,2) - J(1,2)*J(2,1)
       Ji(1,1) =  J(2,2)/detJ;  Ji(1,2) = -J(1,2)/detJ
       Ji(2,1) = -J(2,1)/detJ;  Ji(2,2) =  J(1,1)/detJ
       dNx = Ji(1,1)*dNr + Ji(1,2)*dNs;  dNy = Ji(2,1)*dNr + Ji(2,2)*dNs
       eps = 0.0_dp;  kap = 0.0_dp
       do a = 1, 4
          c = 6*(a-1)
          eps(1) = eps(1) + dNx(a)*pl(c+1);  eps(2) = eps(2) + dNy(a)*pl(c+2)
          eps(3) = eps(3) + dNy(a)*pl(c+1) + dNx(a)*pl(c+2)
          kap(1) = kap(1) + dNx(a)*pl(c+5);  kap(2) = kap(2) - dNy(a)*pl(c+4)
          kap(3) = kap(3) + dNy(a)*pl(c+5) - dNx(a)*pl(c+4)
       end do
       Em = sh%Ee(e);  nu = sh%nue(e);  t = sh%te(e);  fac = Em/(1.0_dp - nu*nu)
       s1 = 0.0_dp
       do side = 1, 2
          z = merge(0.5_dp*t, -0.5_dp*t, side == 1)
          sg(1) = fac*((eps(1) + z*kap(1)) + nu*(eps(2) + z*kap(2)))
          sg(2) = fac*(nu*(eps(1) + z*kap(1)) + (eps(2) + z*kap(2)))
          sg(3) = fac*0.5_dp*(1.0_dp - nu)*(eps(3) + z*kap(3))
          s2 = sqrt(max(sg(1)**2 - sg(1)*sg(2) + sg(2)**2 + 3.0_dp*sg(3)**2, 0.0_dp))
          s1 = max(s1, s2)
       end do
       do a = 1, 4
          nd = sh%conn(a,e)
          vm(nd) = vm(nd) + s1*sh%Ae(e)/4.0_dp
       end do
    end do
    do nd = 1, sh%nn
       if (sh%anode(nd) > 0.0_dp) vm(nd) = vm(nd)/sh%anode(nd)
    end do
  end subroutine shell_stress

  ! kinetic + strain energy of the trial state (checks)
  subroutine shell_energy(sh, ekin, estr)
    type(shell_t), intent(inout) :: sh
    real(dp), intent(out) :: ekin, estr
    real(dp) :: Ef(3,3), pl(24)
    integer :: e, i, c
    ekin = 0.0_dp;  estr = 0.0_dp
    do i = 1, sh%nn
       do c = 1, 6
          ekin = ekin + 0.5_dp*sh%mass(c,i)*sh%v1(c,i)**2
       end do
    end do
    do e = 1, sh%ne
       call elem_local(sh, e, sh%u1, sh%R1, Ef, pl)
       estr = estr + 0.5_dp*sum(pl*matmul(sh%Kl(:,:,e), pl))
    end do
  end subroutine shell_energy

  !---------------------------------------------------------------------
  ! first natural frequency by inverse iteration on the linear system
  ! K x = lambda M x (fixed dofs removed); 0 when the structure is free
  !---------------------------------------------------------------------
  subroutine first_frequency(sh)
    type(shell_t), intent(inout) :: sh
    real(dp), allocatable :: x(:), y(:), fint(:,:), ut(:,:), kx(:)
    integer :: n, i, c, g, it
    real(dp) :: xm, lam, xk
    logical :: ok
    n = sh%nn
    sh%omega1 = 0.0_dp
    if (.not. any(sh%fixed)) return
    allocate(x(sh%ndof), y(sh%ndof), fint(6,n), ut(6,n), kx(sh%ndof))
    sh%Ab = 0.0_dp
    ut = 0.0_dp
    call internal_force(sh, ut, sh%R, fint, .true., 1.0_dp)
    do i = 1, n
       do c = 1, 6
          g = gdof(sh, i, c)
          if (sh%fixed(c,i)) sh%Ab(0,g) = sh%Ab(0,g) + 1.0e30_dp
       end do
    end do
    call band_cholesky(sh, ok)
    if (.not. ok) return
    ! start vector: a smooth translation pattern
    do i = 1, n
       do c = 1, 6
          g = gdof(sh, i, c)
          x(g) = merge(0.0_dp, 1.0_dp + 0.1_dp*mod(i, 7), sh%fixed(c,i) .or. c > 3)
       end do
    end do
    lam = 0.0_dp
    do it = 1, 60
       do i = 1, n
          do c = 1, 6
             g = gdof(sh, i, c)
             y(g) = sh%mass(c,i)*x(g)
             if (sh%fixed(c,i)) y(g) = 0.0_dp
          end do
       end do
       kx = y                                       ! M x
       call band_solve(sh, y)                       ! y = K^-1 M x
       xm = 0.0_dp;  xk = 0.0_dp
       do i = 1, n
          do c = 1, 6
             g = gdof(sh, i, c)
             xm = xm + y(g)*sh%mass(c,i)*y(g)
             xk = xk + y(g)*kx(g)                   ! y . M x = y . K y
          end do
       end do
       if (xm <= 0.0_dp) exit
       lam = xk/xm
       x = y/sqrt(xm)
    end do
    if (lam > 0.0_dp) sh%omega1 = sqrt(lam)
    deallocate(x, y, fint, ut, kx)
  end subroutine first_frequency

  !=====================================================================
  ! Mode m of the linear free vibration about the reference configuration
  ! (inverse iteration with M-orthogonal deflation of the lower modes;
  ! lumped mass, the fixed dofs held): phi (6,nn) [translations, rotations],
  ! M-normalized, and its frequency [Hz] (24c: the pluck of FSI_Test_Motion).
  ! Overwrites the factor in sh%Ab (every structural solve re-assembles it).
  !=====================================================================
  subroutine shell_mode(sh, m, phi, freq)
    type(shell_t), intent(inout) :: sh
    integer, intent(in) :: m
    real(dp), intent(out) :: phi(:,:), freq
    real(dp), allocatable :: x(:,:), y(:), fint(:,:), ut(:,:), mv(:)
    integer :: n, i, c, g, it, k, l
    real(dp) :: xm, lam, pr, lam_old
    logical :: ok
    n = sh%nn;  phi = 0.0_dp;  freq = 0.0_dp
    allocate(x(sh%ndof, m), y(sh%ndof), fint(6,n), ut(6,n), mv(sh%ndof))
    sh%Ab = 0.0_dp;  ut = 0.0_dp
    call internal_force(sh, ut, sh%R, fint, .true., 1.0_dp)
    do i = 1, n
       do c = 1, 6
          g = gdof(sh, i, c);  mv(g) = merge(0.0_dp, sh%mass(c,i), sh%fixed(c,i))
          if (sh%fixed(c,i)) sh%Ab(0,g) = sh%Ab(0,g) + 1.0e30_dp
       end do
    end do
    call band_cholesky(sh, ok)
    if (.not. ok) return
    lam = 0.0_dp
    do k = 1, m
       do i = 1, n                                  ! start: translations with k - 1 sign changes along the numbering
          do c = 1, 6
             g = gdof(sh, i, c)
             x(g, k) = merge(0.0_dp, 1.0_dp + 0.3_dp*cos(k*3.14159_dp*real(i, dp)/n) + 0.05_dp*mod(i, 5), mv(g) == 0.0_dp)
          end do
       end do
       lam = 0.0_dp
       do it = 1, 400
          do l = 1, k - 1                           ! M-orthogonal to the lower modes
             pr = sum(x(:, l)*mv*x(:, k));  x(:, k) = x(:, k) - pr*x(:, l)
          end do
          y = mv*x(:, k)
          call band_solve(sh, y)                    ! y = K^-1 M x
          do l = 1, k - 1
             pr = sum(x(:, l)*mv*y);  y = y - pr*x(:, l)
          end do
          xm = sqrt(max(sum(y*mv*y), 1.0e-300_dp))
          lam_old = lam
          lam = 1.0_dp/max(sum(y*mv*x(:, k))/max(sum(x(:, k)*mv*x(:, k)), 1.0e-300_dp), 1.0e-300_dp)   ! omega^2 (Rayleigh)
          x(:, k) = y/xm
          if (it > 20 .and. abs(lam - lam_old) <= 1.0e-10_dp*lam) exit
       end do
    end do
    freq = sqrt(max(lam, 0.0_dp))/(2.0_dp*3.141592653589793_dp)
    do i = 1, n
       do c = 1, 6
          phi(c, i) = x(gdof(sh, i, c), m)
       end do
    end do
  end subroutine shell_mode

  subroutine shell_free(sh)
    type(shell_t), intent(inout) :: sh
    if (allocated(sh%X0)) deallocate(sh%X0)
    if (allocated(sh%conn)) deallocate(sh%conn, sh%te, sh%Ee, sh%nue, sh%rhoe, sh%fixed)
    if (allocated(sh%E0)) deallocate(sh%E0, sh%xl0, sh%Ae, sh%Kl, sh%mass, sh%anode)
    if (allocated(sh%perm)) deallocate(sh%perm)
    if (allocated(sh%Ab)) deallocate(sh%Ab)
    if (allocated(sh%u)) deallocate(sh%u, sh%v, sh%a, sh%R, sh%fext_n, sh%u1, sh%v1, sh%a1, sh%R1, sh%dth)
    sh%nn = 0;  sh%ne = 0;  sh%ready = .false.
  end subroutine shell_free

end module fsi_shell
