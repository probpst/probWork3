!=======================================================================
! turbulence.F90 -- RANS turbulence models, IDDES hybrid RANS/LES and
! wall functions for probCFD_부산 (update 9).
!
!   TURB_model : 0 none, 1 Spalart-Allmaras (1994, negative-free form),
!                2 k-epsilon (Launder-Spalding, standard), 3 k-omega
!                (Wilcox 1988), 4 k-omega SST (Menter 2003)
!   IDDES_model: 1 = IDDES (Shur, Spalart, Strelets, Travin 2008) on top
!                of the RANS model: the RANS length scale of the
!                destruction term is replaced by the hybrid length
!                l_IDDES = fd_tilde (1 + fe) l_RANS + (1 - fd_tilde) l_LES,
!                l_LES = C_DES Delta_IDDES.  With TURB_model = 0 the
!                WMLES branch alone gives an algebraic subgrid viscosity
!                nu_t = (C_s Delta_IDDES)^2 |S| with the IDDES grid scale
!                (no RANS model in the attached boundary layers).
!   Wall function (bcs flag 11 on u, v, w): tau_w from the law of the
!     wall at the first cell centre, imposed as an effective wall
!     viscosity mu_eff = tau_w y1 / |u1| on the wall face (mu_h); the
!     first-cell values of the turbulence variables follow the model:
!       SA        : Spalding's profile for u_tau; nu_tilde1 = kappa u_tau y1
!       k-epsilon : log law; k1 = u_tau^2/sqrt(C_mu), eps1 = u_tau^3/(kappa y1)
!       k-omega/SST: Menter's automatic wall treatment - u_tau blended
!                   between the viscous and the log solution,
!                   omega1 = sqrt(omega_vis^2 + omega_log^2), k1 blended
!       none      : k-omega SST wall function (Spalding-type u_tau only)
!
!   Transport of the model variables: cell-centred, contravariant face
!   fluxes of the projected velocity, 1st- or 2nd-order upwind
!   convection (Cnvt_scheme_Turb), diffusion with the orthogonal part of
!   the metric (cross terms dropped), implicit Euler with a linearised
!   destruction term, solved by SOR sweeps (Max/Min_Itr_Turb,
!   Rsd_Lmt_Turb, SOR_fct_Turb).  Cnvt_mode_Turb / Dffs_mode_Turb = 1
!   evaluate the term explicitly (sub-cycled to cfl_switch_Turb /
!   vn_switch_Turb), = 2 implicitly.  Ghost / block-interface / overset
!   values are refreshed through the common exchange (field ids 14, 15).
!
!   Quantities are non-dimensional as in the rest of the solver
!   (velocity / Uch, length / Lch, time / tch); nu = 1/Re.
!=======================================================================
module turbulence
  use global_vars
  use sixdof, only: sixdof_wall_velocity
  use mpi_vars
  use ovs_geom
  use b2b_runtime, only: refresh_all
  implicit none
  private
  public :: init_turbulence, solve_turbulence, turb_monitor, write_wallfn_data, mu_lam, rho_lam, wall_utau, &
            velocity_gradient, velocity_gradient_compact, strain_vorticity

  ! ---- model constants ---------------------------------------------
  real(dp), parameter :: kappa = 0.41_dp, Bwall = 5.0_dp
  ! SA
  real(dp), parameter :: cb1 = 0.1355_dp, cb2 = 0.622_dp, sig_sa = 2.0_dp/3.0_dp, cv1 = 7.1_dp, &
                         cw2 = 0.3_dp, cw3 = 2.0_dp, cw1 = cb1/kappa**2 + (1.0_dp + cb2)/sig_sa
  ! k-epsilon
  real(dp), parameter :: cmu = 0.09_dp, ce1 = 1.44_dp, ce2 = 1.92_dp, sig_k_e = 1.0_dp, sig_e = 1.3_dp
  ! k-omega (Wilcox 1988)
  real(dp), parameter :: bstar = 0.09_dp, beta_w = 0.075_dp, alpha_w = 5.0_dp/9.0_dp, sig_k_w = 0.5_dp, sig_w_w = 0.5_dp
  ! SST (Menter 2003)
  real(dp), parameter :: sk1 = 0.85_dp, sw1 = 0.5_dp, b1 = 0.075_dp, g1 = 5.0_dp/9.0_dp
  real(dp), parameter :: sk2 = 1.0_dp,  sw2 = 0.856_dp, b2 = 0.0828_dp, g2 = 0.44_dp, a1 = 0.31_dp
  ! IDDES
  real(dp), parameter :: cdes_sa = 0.65_dp, cdes_kw1 = 0.78_dp, cdes_kw2 = 0.61_dp, cw_idd = 0.15_dp, &
                         cl_sa = 3.55_dp, ct_sa = 1.63_dp, cl_kw = 5.0_dp, ct_kw = 1.87_dp, c_smag = 0.12_dp

  ! ---- wall surface (for the wall distance) --------------------------
  type(surf_t) :: wsurf_t
  logical :: have_wall = .false.
  ! ---- monitor -------------------------------------------------------
  integer  :: mon_it = 0
  real(dp) :: mon_res = 0.0_dp, mon_mutmax = 0.0_dp, mon_yp_min = 0.0_dp, mon_yp_max = 0.0_dp, mon_cf = 0.0_dp
  integer  :: nsub_turb = 1
  integer  :: ierr_t
  ! freestream values
  real(dp) :: k_inf, w_inf, e_inf, nut_inf
  ! per-domain work arrays of the scalar transport (update 16: the sweeps run over all domains)
  type :: twk_t
     real(dp), allocatable :: sc(:,:,:,:), sp(:,:,:,:), gam(:,:,:,:), fo(:,:,:), ac(:,:,:), an(:,:,:,:), bb(:,:,:)
  end type twk_t
  type(twk_t), allocatable :: twk(:)

contains

  !=====================================================================
  ! set-up: wall surface, wall distance, wall-function faces, initial
  ! turbulence field
  !=====================================================================
  subroutine init_turbulence
    integer :: ib, d, i, j, k, b, n, ax, sd
    real(dp) :: p(3), dist, q(3), wbc(3), nu, dmax, hmx
    integer :: tb, region, t
    real(dp), allocatable :: vd(:), hmax_t(:)
    if (.not. (turb_on .or. wallfn_on)) return
    nu = 1.0_dp/Re
    if (myrank == 0) write(*,'(A)') '  turbulence set-up: building the wall surface ...';  flush(6)
    call build_wall_surface
    if (myrank == 0) write(*,'(A,I0,A)') '  turbulence set-up: wall surface done (', wsurf_t%nt, ' triangles), wall distance ...';  flush(6)
    ! ---- wall distance of every cell (including ghosts) ---------------
    ! nearest wall VERTEX by brute force (a few thousand vertices), then
    ! the exact distance to the triangles that have a vertex within
    ! d_vertex + h_max of the point (all others are farther away)
    if (have_wall) then
       allocate(vd(wsurf_t%nv), hmax_t(wsurf_t%nt))
       hmax_t = 0.0_dp
       do t = 1, wsurf_t%nt
          hmax_t(t) = max(norm2(wsurf_t%v(:,wsurf_t%t(1,t)) - wsurf_t%v(:,wsurf_t%t(2,t))), &
                          norm2(wsurf_t%v(:,wsurf_t%t(2,t)) - wsurf_t%v(:,wsurf_t%t(3,t))), &
                          norm2(wsurf_t%v(:,wsurf_t%t(3,t)) - wsurf_t%v(:,wsurf_t%t(1,t))))
       end do
       hmx = maxval(hmax_t)
    end if
    do ib = 1, nown
       d = own_list(ib);  call set_block(d)
       do k = lbound(wdist,3), ubound(wdist,3)
        do j = lbound(wdist,2), ubound(wdist,2)
         do i = lbound(wdist,1), ubound(wdist,1)
            if (.not. have_wall) then
               wdist(i,j,k) = 1.0e30_dp;  cycle
            end if
            p = cell_center_local(i, j, k)
            do t = 1, wsurf_t%nv
               vd(t) = (wsurf_t%v(1,t)-p(1))**2 + (wsurf_t%v(2,t)-p(2))**2 + (wsurf_t%v(3,t)-p(3))**2
            end do
            dist = sqrt(minval(vd))
            dmax = dist                                    ! exact distance <= vertex distance
            do t = 1, wsurf_t%nt
               if (min(vd(wsurf_t%t(1,t)), vd(wsurf_t%t(2,t)), vd(wsurf_t%t(3,t))) > (dmax + hmax_t(t))**2) cycle
               call tri_closest(p, wsurf_t%v(:,wsurf_t%t(1,t)), wsurf_t%v(:,wsurf_t%t(2,t)), wsurf_t%v(:,wsurf_t%t(3,t)), q, wbc, region)
               dmax = min(dmax, norm2(p - q))
            end do
            wdist(i,j,k) = max(dmax, 1.0e-12_dp)/Lch          ! non-dimensional (the node coordinates are in metres; update 16 fix)
         end do
        end do
       end do
       ! ---- wall-function faces (flag 11): wall direction of the cell --
       wface = 0;  wfac = 1.0_dp
       b = dom(d)%b
       do n = 1, blk(b)%nbc
          if (.not. all(blk(b)%bc_flag(1:3,n) == 11)) cycle
          ax = abs(blk(b)%bc_dir(n));  sd = merge(1, -1, blk(b)%bc_dir(n) > 0)
          call mark_wall_cells(d, n, ax, sd)
       end do
    end do
    ! ---- freestream / initial values ---------------------------------
    ! turbulence intensity 1 %, viscosity ratio 10 (non-dimensional)
    k_inf   = 1.5_dp*(0.01_dp)**2
    nut_inf = 10.0_dp*nu
    w_inf   = k_inf/nut_inf
    e_inf   = max(cmu*k_inf**2/nut_inf, cmu**0.75_dp*k_inf**1.5_dp/0.07_dp)
    do ib = 1, nown
       d = own_list(ib);  call set_block(d)
       select case (TURB_model)
       case (1);  tk = 3.0_dp*nu;  tw = 0.0_dp
       case (2);  tk = k_inf;  tw = e_inf
       case (3, 4);  tk = k_inf;  tw = w_inf
       case default;  tk = 0.0_dp;  tw = 0.0_dp
       end select
       mut = 0.0_dp
    end do
    call eddy_viscosity
    if (allocated(vd)) deallocate(vd, hmax_t)
    if (myrank == 0) then
       write(*,'(A,L1,A,I0,A)') '  turbulence set-up: wall surface ', have_wall, ' (', wsurf_t%nt, ' triangles), wall distance computed'
       if (wallfn_on) write(*,'(A)') '  wall function    : flag-11 patches active (tau_w from the law of the wall at the first cell)'
    end if
  end subroutine init_turbulence

  ! cell centre from the node coordinates of the domain's block
  function cell_center_local(i, j, k) result(p)
    integer, intent(in) :: i, j, k
    real(dp) :: p(3)
    integer :: b, ii, jj, kk
    b = dom(cur_dom)%b
    ii = min(max(i, 1), blk(b)%nx);  jj = min(max(j, 1), blk(b)%ny);  kk = min(max(k, 1), blk(b)%nz)
    p = [sum(blk(b)%xn(ii:ii+1,jj:jj+1,kk:kk+1)), sum(blk(b)%yn(ii:ii+1,jj:jj+1,kk:kk+1)), &
         sum(blk(b)%zn(ii:ii+1,jj:jj+1,kk:kk+1))]/8.0_dp
  end function cell_center_local

  ! all no-slip patches (flag 1 with zero velocity, or flag 11) -> triangles
  subroutine build_wall_surface
    integer :: b, n, ax, sd, i1, i2, j1, j2, k1, k2, ntri, n1, n2, a1, a2, v0, ip, off(3), t1, t2
    real(dp), allocatable :: fx(:,:), fy(:,:), fz(:,:), fd(:,:,:)
    integer :: m(3), mi(3)
    ntri = 0
    do b = 1, nblk
       do n = 1, blk(b)%nbc
          if (.not. patch_is_wall(b, n)) cycle
          i1 = blk(b)%bc_irng(1,n); i2 = blk(b)%bc_irng(2,n)
          j1 = blk(b)%bc_jrng(1,n); j2 = blk(b)%bc_jrng(2,n)
          k1 = blk(b)%bc_krng(1,n); k2 = blk(b)%bc_krng(2,n)
          ntri = ntri + 2*max(i2-i1,1)*max(j2-j1,1)*max(k2-k1,1)
       end do
    end do
    have_wall = (ntri > 0)
    if (.not. have_wall) return
    call surf_init(wsurf_t, 2*ntri + 8, ntri)
    ip = 0
    do b = 1, nblk
       do n = 1, blk(b)%nbc
          if (.not. patch_is_wall(b, n)) cycle
          ax = abs(blk(b)%bc_dir(n));  sd = merge(1, -1, blk(b)%bc_dir(n) > 0)
          i1 = blk(b)%bc_irng(1,n); i2 = blk(b)%bc_irng(2,n)
          j1 = blk(b)%bc_jrng(1,n); j2 = blk(b)%bc_jrng(2,n)
          k1 = blk(b)%bc_krng(1,n); k2 = blk(b)%bc_krng(2,n)
          select case (ax)
          case (1);  t1 = 2;  t2 = 3;  n1 = j2-j1+1;  n2 = k2-k1+1
          case (2);  t1 = 1;  t2 = 3;  n1 = i2-i1+1;  n2 = k2-k1+1
          case default;  t1 = 1;  t2 = 2;  n1 = i2-i1+1;  n2 = j2-j1+1
          end select
          allocate(fx(n1,n2), fy(n1,n2), fz(n1,n2), fd(3,n1,n2))
          off = 0;  off(ax) = sd
          do a2 = 1, n2
             do a1 = 1, n1
                m = [i1, j1, k1];  m(t1) = m(t1) + a1 - 1;  m(t2) = m(t2) + a2 - 1
                mi = m + off
                fx(a1,a2) = blk(b)%xn(m(1),m(2),m(3));  fy(a1,a2) = blk(b)%yn(m(1),m(2),m(3));  fz(a1,a2) = blk(b)%zn(m(1),m(2),m(3))
                fd(:,a1,a2) = [blk(b)%xn(mi(1),mi(2),mi(3)) - fx(a1,a2), blk(b)%yn(mi(1),mi(2),mi(3)) - fy(a1,a2), &
                               blk(b)%zn(mi(1),mi(2),mi(3)) - fz(a1,a2)]
             end do
          end do
          ip = ip + 1
          call surf_add_patch(wsurf_t, ip, n1, n2, fx, fy, fz, fd, v0)
          deallocate(fx, fy, fz, fd)
       end do
    end do
    wsurf_t%tb(1:wsurf_t%nt) = 1
    call surf_finish(wsurf_t)
  end subroutine build_wall_surface

  logical function patch_is_wall(b, n)
    integer, intent(in) :: b, n
    patch_is_wall = (all(blk(b)%bc_flag(1:3,n) == 1) .and. all(blk(b)%bc_val(1:3,n) == 0.0_dp)) &
                    .or. all(blk(b)%bc_flag(1:3,n) == 11)
  end function patch_is_wall

  ! cells adjacent to the flag-11 patch n of the current domain's block
  subroutine mark_wall_cells(d, n, ax, sd)
    integer, intent(in) :: d, n, ax, sd
    integer :: b, lo(3), hi(3), i, j, k, rr(2,3)
    b = dom(d)%b
    rr(:,1) = blk(b)%bc_irng(:,n);  rr(:,2) = blk(b)%bc_jrng(:,n);  rr(:,3) = blk(b)%bc_krng(:,n)
    lo = [max(rr(1,1), dom(d)%r1(1)), max(rr(1,2), dom(d)%r1(2)), max(rr(1,3), dom(d)%r1(3))]
    hi = [min(rr(2,1)-1, dom(d)%r2(1)), min(rr(2,2)-1, dom(d)%r2(2)), min(rr(2,3)-1, dom(d)%r2(3))]
    if (sd > 0) then
       lo(ax) = 1;  hi(ax) = 1
    else
       lo(ax) = ncell(b, ax);  hi(ax) = lo(ax)
    end if
    if (lo(ax) < dom(d)%r1(ax) .or. lo(ax) > dom(d)%r2(ax)) return
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             wface(i,j,k) = sd*ax
          end do
       end do
    end do
  end subroutine mark_wall_cells

  !=====================================================================
  ! one time step: wall functions, sources, transport, eddy viscosity
  !=====================================================================
  subroutine solve_turbulence
    if (.not. (turb_on .or. wallfn_on)) return
    if (wallfn_on) call wall_functions
    if (TURB_model > 0) call transport_models
    call eddy_viscosity
    call refresh_all(16)
  end subroutine solve_turbulence

  !---------------------------------------------------------------------
  ! wall functions on the flag-11 cells: u_tau from the law of the wall
  ! at the first cell centre (tangential velocity, wall distance), the
  ! wall-face viscosity factor, and the first-cell turbulence values
  !---------------------------------------------------------------------
  subroutine wall_functions
    integer :: ib, d, i, j, k, ax
    real(dp) :: nu, uc(3), nrm(3), ut(3), umag, y, utau, yp, mulam, ypmin, ypmax, r1, r2, wv, wl, kk, om, cfsum, ncf
    nu = 1.0_dp/Re
    ypmin = huge(1.0_dp);  ypmax = 0.0_dp;  cfsum = 0.0_dp;  ncf = 0.0_dp
    do ib = 1, nown
       d = own_list(ib);  call set_block(d)
       do k = ks, ke
        do j = js, je
         do i = is, ie
            if (wface(i,j,k) == 0) cycle
            ax = abs(wface(i,j,k))
            uc = cell_velocity(i, j, k)
            uc = uc - (1.0_dp - bc_ramp)*[u_init, v_init, w_init]/Uch        ! relative to the ramped wall (towing-tank start; update 20)
            if (any_moving) uc = uc - sixdof_wall_velocity(cur_blk, Xh(:,2*i,2*j,2*k), [2*i, 2*j, 2*k])   ! relative to a moving wall
            nrm = face_normal(ax, i, j, k)
            ut = uc - dot_product(uc, nrm)*nrm
            umag = max(norm2(ut), 1.0e-12_dp)
            y = wdist(i,j,k)
            mulam = mu_lam(i, j, k)
            select case (TURB_model)
            case (2)
               ! k-epsilon: u_tau from k (log law), viscous limit for small y+
               utau = cmu**0.25_dp*sqrt(max(tk(i,j,k), 1.0e-20_dp))
               yp = utau*y/(nu*mulam)
               if (yp < 11.0_dp) utau = sqrt(umag*nu*mulam/y)
            case default
               utau = spalding_utau(umag, y, nu*mulam)
            end select
            yp = utau*y/(nu*mulam)
            ypmin = min(ypmin, yp);  ypmax = max(ypmax, yp)
            cfsum = cfsum + 2.0_dp*rho_lam(i,j,k)*utau**2;  ncf = ncf + 1.0_dp     ! Cf = tau_w / (1/2 rho U^2), non-dim
            ! effective wall viscosity: mu_eff u1/y1 = rho u_tau^2 (both / rho1 mu1 units -> mu_eff/mu1 = Re u_tau^2 y / u1 * rho)
            wfac(i,j,k) = max(Re*rho_lam(i,j,k)*utau**2*y/umag, mulam)
            ! first-cell values of the model variables
            select case (TURB_model)
            case (1)
               tk(i,j,k) = kappa*utau*y                          ! nu_tilde (log region)
            case (2)
               tk(i,j,k) = utau**2/sqrt(cmu)
               tw(i,j,k) = utau**3/(kappa*y)
            case (3, 4)
               ! Menter's automatic wall treatment
               wv = 6.0_dp*nu*mulam/(b1*y**2)
               wl = utau/(sqrt(bstar)*kappa*y)
               om = sqrt(wv**2 + wl**2)
               kk = utau**2/sqrt(bstar)
               r1 = min(yp/11.0_dp, 1.0_dp)
               tw(i,j,k) = om
               tk(i,j,k) = r1*kk + (1.0_dp - r1)*tk(i,j,k)
            end select
         end do
        end do
       end do
    end do
    mon_yp_min = -global_max(-ypmin);  mon_yp_max = global_max(ypmax)
    mon_cf = global_sum(cfsum)/max(global_sum(ncf), 1.0_dp)
    if (.not. turb_on) return
    call refresh_all(14);  if (TURB_model >= 2) call refresh_all(15)
  end subroutine wall_functions

  ! u_tau from Spalding's law of the wall (Newton iteration on u_tau)
  real(dp) function spalding_utau(u, y, nu) result(ut)
    real(dp), intent(in) :: u, y, nu
    real(dp) :: up, yp, f, df, ku, ekb
    integer :: it
    ekb = exp(-kappa*Bwall)
    ut = max(sqrt(u*nu/y), 1.0e-10_dp)                ! viscous guess
    do it = 1, 30
       up = u/ut;  yp = y*ut/nu;  ku = kappa*up
       f  = up + ekb*(exp(ku) - 1.0_dp - ku - ku**2/2.0_dp - ku**3/6.0_dp) - yp
       ! df/dut : d(up)/dut = -up/ut ; d(yp)/dut = y/nu
       df = (-up/ut)*(1.0_dp + ekb*kappa*(exp(ku) - 1.0_dp - ku - ku**2/2.0_dp)) - y/nu
       if (abs(df) < 1.0e-30_dp) exit
       ut = ut - f/df
       ut = max(ut, 1.0e-10_dp)
       if (abs(f) < 1.0e-8_dp*max(yp, 1.0_dp)) exit
    end do
  end function spalding_utau

  ! Cartesian velocity at the cell centre (average of the face values)
  function cell_velocity(i, j, k) result(uc)
    integer, intent(in) :: i, j, k
    real(dp) :: uc(3)
    uc = [0.5_dp*(u(i-1,j,k) + u(i,j,k)), 0.5_dp*(v(i,j-1,k) + v(i,j,k)), 0.5_dp*(w(i,j,k-1) + w(i,j,k))]
  end function cell_velocity

  ! unit normal of the wall face of the cell (from the metric of the face)
  function face_normal(ax, i, j, k) result(nrm)
    integer, intent(in) :: ax, i, j, k
    real(dp) :: nrm(3), s
    nrm = [mh(1+3*(ax-1),2*i,2*j,2*k), mh(2+3*(ax-1),2*i,2*j,2*k), mh(3+3*(ax-1),2*i,2*j,2*k)]
    s = norm2(nrm)
    if (s > 0.0_dp) then
       nrm = nrm/s
    else
       nrm = [0.0_dp, 0.0_dp, 1.0_dp]
    end if
  end function face_normal

  ! laminar viscosity / density of the cell (three-phase with the
  ! cavitation vapour fraction, update 15)
  real(dp) function mu_lam(i, j, k)
    integer, intent(in) :: i, j, k
    real(dp) :: a, av
    a = 1.0_dp
    if (Num_of_Phase == 2) a = alpha(i,j,k)
    mu_lam = a + mu_ratio*(1.0_dp - a)
    if (cav_on) then
       av = min(max(vap(i,j,k), 0.0_dp), 1.0_dp)
       mu_lam = a*((1.0_dp - av) + av*muv_ratio) + mu_ratio*(1.0_dp - a)
    end if
  end function mu_lam

  real(dp) function rho_lam(i, j, k)
    integer, intent(in) :: i, j, k
    real(dp) :: a, av
    a = 1.0_dp
    if (Num_of_Phase == 2) a = alpha(i,j,k)
    rho_lam = a + rho_ratio*(1.0_dp - a)
    if (cav_on) then
       av = min(max(vap(i,j,k), 0.0_dp), 1.0_dp)
       rho_lam = a*((1.0_dp - av) + av*rv_ratio) + rho_ratio*(1.0_dp - a)
    end if
  end function rho_lam

  ! Reboud correction (Reboud = 1): the density of the eddy viscosity in
  ! a cavitating cell is f(rho) = rho_v + (rho - rho_v)^n/(rho_l - rho_v)^(n-1),
  ! n = 10 (rho_l = 1 non-dimensional), which damps mu_t in the mixture
  ! and lets the re-entrant jet / cloud shedding develop
  real(dp) function rho_reboud(rho)
    real(dp), intent(in) :: rho
    real(dp), parameter :: nexp = 10.0_dp
    real(dp) :: rl, rvv
    rl = 1.0_dp;  rvv = rv_ratio
    if (Reboud /= 1 .or. .not. cav_on .or. rho >= rl .or. rl <= rvv) then
       rho_reboud = rho;  return
    end if
    rho_reboud = rvv + (max(rho, rvv) - rvv)**nexp/(rl - rvv)**(nexp - 1.0_dp)
  end function rho_reboud

  !---------------------------------------------------------------------
  ! velocity gradient tensor at the cell centre: g(c,l) = du_c/dx_l
  !---------------------------------------------------------------------
  subroutine velocity_gradient(i, j, k, g)
    integer, intent(in) :: i, j, k
    real(dp), intent(out) :: g(3,3)
    real(dp) :: dq(3,3), up(3), um(3)
    integer :: q, c, l
    ! d u_c / d xi_q by central differences of the cell-centred velocity
    up = cell_velocity(i+1,j,k);  um = cell_velocity(i-1,j,k);  dq(:,1) = 0.5_dp*(up - um)
    up = cell_velocity(i,j+1,k);  um = cell_velocity(i,j-1,k);  dq(:,2) = 0.5_dp*(up - um)
    up = cell_velocity(i,j,k+1);  um = cell_velocity(i,j,k-1);  dq(:,3) = 0.5_dp*(up - um)
    do l = 1, 3
       do c = 1, 3
          g(c,l) = 0.0_dp
          do q = 1, 3
             g(c,l) = g(c,l) + dq(c,q)*mh(l+3*(q-1),2*i,2*j,2*k)
          end do
       end do
    end do
  end subroutine velocity_gradient

  !---------------------------------------------------------------------
  ! velocity gradient tensor at the cell centre from the FACE DIFFERENCES of the staggered
  ! velocities (update 23, the statistics): d u_c / d xi_q is the compact difference of the two
  ! faces of the cell when the component lives on those faces (u in xi, v in eta, w in zeta -
  ! the difference the viscous term and the divergence of the solver use), and the central
  ! difference of the cell-centred values across the two neighbouring cells otherwise (the
  ! average of the compact differences at the four edges of the cell in that plane).  The
  ! compact own-direction derivatives resolve the 2-cell waves the central difference loses.
  !---------------------------------------------------------------------
  subroutine velocity_gradient_compact(i, j, k, g)
    integer, intent(in) :: i, j, k
    real(dp), intent(out) :: g(3,3)
    real(dp) :: dq(3,3), up(3), um(3)
    integer :: q, c, l
    up = cell_velocity(i+1,j,k);  um = cell_velocity(i-1,j,k);  dq(:,1) = 0.5_dp*(up - um)
    up = cell_velocity(i,j+1,k);  um = cell_velocity(i,j-1,k);  dq(:,2) = 0.5_dp*(up - um)
    up = cell_velocity(i,j,k+1);  um = cell_velocity(i,j,k-1);  dq(:,3) = 0.5_dp*(up - um)
    dq(1,1) = u(i,j,k) - u(i-1,j,k)                       ! d u / d xi   : the cell's two xi faces
    dq(2,2) = v(i,j,k) - v(i,j-1,k)                       ! d v / d eta  : the cell's two eta faces
    dq(3,3) = w(i,j,k) - w(i,j,k-1)                       ! d w / d zeta : the cell's two zeta faces
    do l = 1, 3
       do c = 1, 3
          g(c,l) = 0.0_dp
          do q = 1, 3
             g(c,l) = g(c,l) + dq(c,q)*mh(l+3*(q-1),2*i,2*j,2*k)
          end do
       end do
    end do
  end subroutine velocity_gradient_compact

  ! |S| = sqrt(2 S_ij S_ij), |Omega| = sqrt(2 W_ij W_ij), S_ij S_ij
  subroutine strain_vorticity(g, smag, wmag)
    real(dp), intent(in) :: g(3,3)
    real(dp), intent(out) :: smag, wmag
    real(dp) :: s2, w2
    integer :: c, l
    s2 = 0.0_dp;  w2 = 0.0_dp
    do l = 1, 3
       do c = 1, 3
          s2 = s2 + (0.5_dp*(g(c,l) + g(l,c)))**2
          w2 = w2 + (0.5_dp*(g(c,l) - g(l,c)))**2
       end do
    end do
    smag = sqrt(2.0_dp*s2);  wmag = sqrt(2.0_dp*w2)
  end subroutine strain_vorticity

  ! largest and wall-normal cell sizes of the cell (index-space metric)
  subroutine cell_sizes(i, j, k, hmax, hwn, ax)
    integer, intent(in) :: i, j, k, ax
    real(dp), intent(out) :: hmax, hwn
    real(dp) :: h(3)
    integer :: q
    do q = 1, 3
       h(q) = 1.0_dp/max(sqrt(gql(q, q, 2*i, 2*j, 2*k)), 1.0e-30_dp)
    end do
    hmax = maxval(h)
    hwn  = h(min(max(ax,1),3))
  end subroutine cell_sizes

  ! IDDES grid scale and blending: returns the hybrid length scale
  ! l_hyb given the RANS length scale l_rans, the wall distance dw, the
  ! eddy viscosity nut (kinematic, non-dim), the laminar nu, |grad u|
  real(dp) function iddes_length(l_rans, dw, nut, nu, gradu, hmax, hwn, cdes, cl, ct) result(lh)
    real(dp), intent(in) :: l_rans, dw, nut, nu, gradu, hmax, hwn, cdes, cl, ct
    real(dp) :: delta, rdt, rdl, fd, fdt, fb, alph, fe1, fe2, fe, ft, fl, lles
    delta = min(max(cw_idd*dw, cw_idd*hmax, hwn), hmax)
    lles  = cdes*delta
    alph  = 0.25_dp - dw/max(hmax, 1.0e-30_dp)
    fb    = min(2.0_dp*exp(-9.0_dp*alph**2), 1.0_dp)
    rdt   = nut/max(kappa**2*dw**2*max(gradu, 1.0e-10_dp), 1.0e-30_dp)
    rdl   = nu /max(kappa**2*dw**2*max(gradu, 1.0e-10_dp), 1.0e-30_dp)
    fdt   = 1.0_dp - tanh((8.0_dp*rdt)**3)
    fd    = max(1.0_dp - fdt, fb)
    ft    = tanh((ct**2*rdt)**3)
    fl    = tanh((cl**2*rdl)**10)
    fe2   = 1.0_dp - max(ft, fl)
    if (alph >= 0.0_dp) then
       fe1 = 2.0_dp*exp(-11.09_dp*alph**2)
    else
       fe1 = 2.0_dp*exp(-9.0_dp*alph**2)
    end if
    fe = max(fe1 - 1.0_dp, 0.0_dp)*fe2
    lh = fd*(1.0_dp + fe)*l_rans + (1.0_dp - fd)*lles
  end function iddes_length

  !---------------------------------------------------------------------
  ! eddy viscosity from the model variables (mu_t / mu1, i.e. Re nu_t
  ! in non-dimensional form: nu_t is in Uch Lch units)
  !---------------------------------------------------------------------
  subroutine eddy_viscosity
    integer :: ib, d, i, j, k
    real(dp) :: nu, g(3,3), smag, wmag, chi, fv1, nut, f2, arg, dw, hmax, hwn, delta, mutmax, rho
    mutmax = 0.0_dp
    nu = 1.0_dp/Re
    do ib = 1, nown
       d = own_list(ib);  call set_block(d)
       do k = ks, ke
        do j = js, je
         do i = is, ie
            if (cstat(i,j,k) < 0) then
               mut(i,j,k) = 0.0_dp;  cycle
            end if
            rho = rho_reboud(rho_lam(i, j, k))
            select case (TURB_model)
            case (1)
               chi = max(tk(i,j,k), 0.0_dp)/nu
               fv1 = chi**3/(chi**3 + cv1**3)
               nut = max(tk(i,j,k), 0.0_dp)*fv1
            case (2)
               nut = cmu*max(tk(i,j,k), 0.0_dp)**2/max(tw(i,j,k), 1.0e-20_dp)
            case (3)
               nut = max(tk(i,j,k), 0.0_dp)/max(tw(i,j,k), 1.0e-20_dp)
            case (4)
               call velocity_gradient(i, j, k, g);  call strain_vorticity(g, smag, wmag)
               dw = wdist(i,j,k)
               arg = max(2.0_dp*sqrt(max(tk(i,j,k),0.0_dp))/(bstar*max(tw(i,j,k),1.0e-20_dp)*dw), &
                         500.0_dp*nu/(dw**2*max(tw(i,j,k),1.0e-20_dp)))
               f2 = tanh(arg**2)
               nut = a1*max(tk(i,j,k), 0.0_dp)/max(a1*max(tw(i,j,k),1.0e-20_dp), smag*f2)
            case default
               nut = 0.0_dp
               if (IDDES_model == 1) then
                  ! WMLES branch only: Smagorinsky viscosity with the IDDES grid scale
                  call velocity_gradient(i, j, k, g);  call strain_vorticity(g, smag, wmag)
                  call cell_sizes(i, j, k, hmax, hwn, wall_axis(i, j, k))
                  dw = wdist(i,j,k)
                  delta = min(max(cw_idd*dw, cw_idd*hmax, hwn), hmax)
                  nut = (c_smag*delta)**2*smag*(1.0_dp - exp(-dw*sqrt(max(smag,0.0_dp)/nu)/26.0_dp))**2
               end if
            end select
            mut(i,j,k) = Re*rho*max(nut, 0.0_dp)
            mutmax = max(mutmax, mut(i,j,k))
         end do
        end do
       end do
    end do
    mon_mutmax = global_max(mutmax)
  end subroutine eddy_viscosity

  ! wall-normal index direction of a cell (the axis along which the wall
  ! distance changes fastest); used for the IDDES wall-normal spacing
  integer function wall_axis(i, j, k) result(ax)
    integer, intent(in) :: i, j, k
    real(dp) :: dd(3)
    dd = [abs(wdist(i+1,j,k) - wdist(i-1,j,k)), abs(wdist(i,j+1,k) - wdist(i,j-1,k)), abs(wdist(i,j,k+1) - wdist(i,j,k-1))]
    ax = maxloc(dd, 1)
  end function wall_axis

  !---------------------------------------------------------------------
  ! transport of the model variables: sources per model, then the
  ! common scalar solver for each variable
  !---------------------------------------------------------------------
  subroutine transport_models
    integer :: ib, d, i, j, k, nvar, iv
    real(dp) :: nu, g(3,3), smag, wmag, dw, hmax, hwn, nut, kk, om, ep, f1, arg, cdkw, gam, bet, sk, sw, lrans, lhyb, &
                chi, fv1, fv2, stil, r, gg, fw, dnt(3), gk(3), gw(3), rho, prodk, dk, dwv, pk, xd
    ! source arrays (per domain, twk): Sc (explicit source, per unit volume), Sp (implicit
    ! destruction coefficient: -Sp * phi), diffusivity per variable
    nu = 1.0_dp/Re
    nvar = merge(1, 2, TURB_model == 1)
    if (.not. allocated(twk)) allocate(twk(ndom))
    do ib = 1, nown
       d = own_list(ib);  call set_block(d)
       allocate(twk(d)%sc(is:ie,js:je,ks:ke,2), twk(d)%sp(is:ie,js:je,ks:ke,2), &
                twk(d)%gam(lbound(tk,1):ubound(tk,1),lbound(tk,2):ubound(tk,2),lbound(tk,3):ubound(tk,3),2))
       associate(sc => twk(d)%sc, sp => twk(d)%sp, gam_v => twk(d)%gam)
       sc = 0.0_dp;  sp = 0.0_dp;  gam_v = nu
       do k = ks, ke
        do j = js, je
         do i = is, ie
            if (cstat(i,j,k) /= 1) cycle
            rho = rho_lam(i, j, k)
            call velocity_gradient(i, j, k, g);  call strain_vorticity(g, smag, wmag)
            dw = wdist(i,j,k)
            call cell_sizes(i, j, k, hmax, hwn, wall_axis(i, j, k))
            nut = mut(i,j,k)/(Re*rho)
            select case (TURB_model)
            case (1)
               !---------------- Spalart-Allmaras ------------------------
               chi = max(tk(i,j,k), 0.0_dp)/nu
               fv1 = chi**3/(chi**3 + cv1**3);  fv2 = 1.0_dp - chi/(1.0_dp + chi*fv1)
               lrans = dw
               if (IDDES_model == 1) lrans = iddes_length(dw, dw, nut, nu, max(smag, wmag), hmax, hwn, cdes_sa, cl_sa, ct_sa)
               stil = max(wmag + tk(i,j,k)*fv2/(kappa**2*lrans**2), 0.3_dp*wmag)
               r  = min(tk(i,j,k)/max(stil*kappa**2*lrans**2, 1.0e-30_dp), 10.0_dp)
               gg = r + cw2*(r**6 - r)
               fw = gg*((1.0_dp + cw3**6)/(gg**6 + cw3**6))**(1.0_dp/6.0_dp)
               dnt = grad_scalar(tk, i, j, k)
               sc(i,j,k,1) = cb1*stil*tk(i,j,k) + cb2/sig_sa*dot_product(dnt, dnt)
               sp(i,j,k,1) = cw1*fw*max(tk(i,j,k), 0.0_dp)/lrans**2
               gam_v(i,j,k,1) = (nu + max(tk(i,j,k),0.0_dp))/sig_sa
            case (2)
               !---------------- k-epsilon -------------------------------
               kk = max(tk(i,j,k), 1.0e-20_dp);  ep = max(tw(i,j,k), 1.0e-20_dp)
               pk = min(nut*smag**2, 10.0_dp*ep)
               lrans = kk**1.5_dp/ep
               dk = ep
               if (IDDES_model == 1) then
                  lhyb = iddes_length(lrans, dw, nut, nu, max(smag, wmag), hmax, hwn, cdes_kw1, cl_kw, ct_kw)
                  dk = kk**1.5_dp/max(lhyb, 1.0e-30_dp)
               end if
               sc(i,j,k,1) = pk;                sp(i,j,k,1) = dk/kk
               sc(i,j,k,2) = ce1*pk*ep/kk;      sp(i,j,k,2) = ce2*ep/kk
               gam_v(i,j,k,1) = nu + nut/sig_k_e;  gam_v(i,j,k,2) = nu + nut/sig_e
            case (3)
               !---------------- k-omega (Wilcox) ------------------------
               kk = max(tk(i,j,k), 1.0e-20_dp);  om = max(tw(i,j,k), 1.0e-20_dp)
               pk = min(nut*smag**2, 10.0_dp*bstar*kk*om)
               lrans = sqrt(kk)/(bstar*om)
               dk = bstar*om
               if (IDDES_model == 1) then
                  lhyb = iddes_length(lrans, dw, nut, nu, max(smag, wmag), hmax, hwn, cdes_kw1, cl_kw, ct_kw)
                  dk = sqrt(kk)/max(lhyb, 1.0e-30_dp)
               end if
               sc(i,j,k,1) = pk;                     sp(i,j,k,1) = dk
               sc(i,j,k,2) = alpha_w*om/kk*pk;       sp(i,j,k,2) = beta_w*om
               gam_v(i,j,k,1) = nu + sig_k_w*nut;    gam_v(i,j,k,2) = nu + sig_w_w*nut
            case (4)
               !---------------- k-omega SST (Menter 2003) ---------------
               kk = max(tk(i,j,k), 1.0e-20_dp);  om = max(tw(i,j,k), 1.0e-20_dp)
               gk = grad_scalar(tk, i, j, k);  gw = grad_scalar(tw, i, j, k)
               cdkw = max(2.0_dp*sw2/om*dot_product(gk, gw), 1.0e-10_dp)
               arg = min(max(sqrt(kk)/(bstar*om*dw), 500.0_dp*nu/(dw**2*om)), 4.0_dp*sw2*kk/(cdkw*dw**2))
               f1 = tanh(arg**4)
               gam = f1*g1 + (1.0_dp - f1)*g2;  bet = f1*b1 + (1.0_dp - f1)*b2
               sk = f1*sk1 + (1.0_dp - f1)*sk2;  sw = f1*sw1 + (1.0_dp - f1)*sw2
               pk = min(nut*smag**2, 10.0_dp*bstar*kk*om)
               lrans = sqrt(kk)/(bstar*om)
               dk = bstar*om
               if (IDDES_model == 1) then
                  lhyb = iddes_length(lrans, dw, nut, nu, max(smag, wmag), hmax, hwn, f1*cdes_kw1 + (1.0_dp-f1)*cdes_kw2, cl_kw, ct_kw)
                  dk = sqrt(kk)/max(lhyb, 1.0e-30_dp)
               end if
               sc(i,j,k,1) = pk;                                      sp(i,j,k,1) = dk
               sc(i,j,k,2) = gam*pk/max(nut, 1.0e-20_dp) + (1.0_dp - f1)*cdkw*0.0_dp + 2.0_dp*(1.0_dp - f1)*sw2/om*max(dot_product(gk, gw), 0.0_dp)
               sp(i,j,k,2) = bet*om + 2.0_dp*(1.0_dp - f1)*sw2/om**2*max(-dot_product(gk, gw), 0.0_dp)
               gam_v(i,j,k,1) = nu + sk*nut;  gam_v(i,j,k,2) = nu + sw*nut
            end select
         end do
        end do
       end do
       ! production treated implicitly (linearised in phi, P = (P/phi_old) phi),
       ! with the growth per time step limited to a factor 10: the explicit
       ! form explodes when (P/phi) dt >> 1 in the start-up transient
       do iv = 1, nvar
          do k = ks, ke
           do j = js, je
            do i = is, ie
               if (cstat(i,j,k) /= 1) cycle
               if (iv == 1) then
                  xd = max(tk(i,j,k), 1.0e-20_dp)
               else
                  xd = max(tw(i,j,k), 1.0e-20_dp)
               end if
               prodk = sc(i,j,k,iv)/xd                     ! growth rate
               sc(i,j,k,iv) = 0.0_dp
               sp(i,j,k,iv) = max(sp(i,j,k,iv) - prodk, -0.9_dp/dt)
            end do
           end do
          end do
       end do
       end associate
    end do
    ! the transport of every variable over ALL owned domains (one collective per
    ! sweep: a per-domain solve would deadlock with uneven domain counts; update 16)
    call scalar_transport(14, 1)
    if (nvar == 2) call scalar_transport(15, 2)
    do ib = 1, nown
       d = own_list(ib)
       deallocate(twk(d)%sc, twk(d)%sp, twk(d)%gam)
    end do
    call refresh_all(14);  if (nvar == 2) call refresh_all(15)
    ! positivity and wall / inlet values
    do ib = 1, nown
       d = own_list(ib);  call set_block(d)
       call wall_and_inlet_values
    end do
    call refresh_all(14);  if (nvar == 2) call refresh_all(15)
  end subroutine transport_models

  ! gradient of a cell-centred scalar at the cell centre
  function grad_scalar(f, i, j, k) result(gr)
    real(dp), intent(in) :: f(lbound(alpha,1):, lbound(alpha,2):, lbound(alpha,3):)
    integer, intent(in) :: i, j, k
    real(dp) :: gr(3), dq(3)
    integer :: l, q
    dq(1) = 0.5_dp*(f(i+1,j,k) - f(i-1,j,k))
    dq(2) = 0.5_dp*(f(i,j+1,k) - f(i,j-1,k))
    dq(3) = 0.5_dp*(f(i,j,k+1) - f(i,j,k-1))
    do l = 1, 3
       gr(l) = 0.0_dp
       do q = 1, 3
          gr(l) = gr(l) + dq(q)*mh(l+3*(q-1),2*i,2*j,2*k)
       end do
    end do
  end function grad_scalar

  !---------------------------------------------------------------------
  ! wall (no-slip, no wall function) and inlet values of the variables
  !---------------------------------------------------------------------
  subroutine wall_and_inlet_values
    integer :: i, j, k, b, n, ax, sd, lo(3), hi(3), rr(2,3), d, ig(3), e(3)
    real(dp) :: nu, mulam, y
    logical :: wall, inlet
    nu = 1.0_dp/Re
    d = cur_dom;  b = dom(d)%b
    ! positivity
    tk(is:ie,js:je,ks:ke) = max(tk(is:ie,js:je,ks:ke), 0.0_dp)
    if (TURB_model >= 2) tw(is:ie,js:je,ks:ke) = max(tw(is:ie,js:je,ks:ke), 1.0e-20_dp)
    do n = 1, blk(b)%nbc
       wall  = (all(blk(b)%bc_flag(1:3,n) == 1) .and. all(blk(b)%bc_val(1:3,n) == 0.0_dp))
       inlet = (all(blk(b)%bc_flag(1:3,n) == 1) .and. any(blk(b)%bc_val(1:3,n) /= 0.0_dp)) &
               .or. all(blk(b)%bc_flag(1:3,n) == 12)          ! incident-wave patch: free-stream values (update 18)
       if (.not. (wall .or. inlet)) cycle
       ax = abs(blk(b)%bc_dir(n));  sd = merge(1, -1, blk(b)%bc_dir(n) > 0)
       rr(:,1) = blk(b)%bc_irng(:,n);  rr(:,2) = blk(b)%bc_jrng(:,n);  rr(:,3) = blk(b)%bc_krng(:,n)
       lo = [max(rr(1,1), dom(d)%r1(1)), max(rr(1,2), dom(d)%r1(2)), max(rr(1,3), dom(d)%r1(3))]
       hi = [min(rr(2,1)-1, dom(d)%r2(1)), min(rr(2,2)-1, dom(d)%r2(2)), min(rr(2,3)-1, dom(d)%r2(3))]
       if (sd > 0) then
          lo(ax) = 1;  hi(ax) = 1
       else
          lo(ax) = ncell(b, ax);  hi(ax) = lo(ax)
       end if
       if (lo(ax) < dom(d)%r1(ax) .or. lo(ax) > dom(d)%r2(ax)) cycle
       e = 0;  e(ax) = -sd                                   ! ghost direction
       do k = lo(3), hi(3)
        do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            ig = [i,j,k] + e
            if (wall) then
               mulam = mu_lam(i, j, k);  y = wdist(i,j,k)
               select case (TURB_model)
               case (1)
                  tk(ig(1),ig(2),ig(3)) = -tk(i,j,k)                     ! nu_tilde = 0 at the wall
               case (2)
                  ! standard k-epsilon has no wall-resolved form: log-law wall values
                  tk(ig(1),ig(2),ig(3)) = tk(i,j,k);  tw(ig(1),ig(2),ig(3)) = tw(i,j,k)
               case (3, 4)
                  tk(ig(1),ig(2),ig(3)) = -tk(i,j,k)                     ! k = 0 at the wall
                  tw(i,j,k) = max(tw(i,j,k), 6.0_dp*nu*mulam/(b1*y**2))  ! omega at the first cell (Menter)
                  tw(ig(1),ig(2),ig(3)) = tw(i,j,k)
               end select
            else
               select case (TURB_model)
               case (1);     tk(ig(1),ig(2),ig(3)) = 2.0_dp*3.0_dp*nu - tk(i,j,k)
               case (2);     tk(ig(1),ig(2),ig(3)) = 2.0_dp*k_inf - tk(i,j,k);  tw(ig(1),ig(2),ig(3)) = 2.0_dp*e_inf - tw(i,j,k)
               case (3, 4);  tk(ig(1),ig(2),ig(3)) = 2.0_dp*k_inf - tk(i,j,k);  tw(ig(1),ig(2),ig(3)) = 2.0_dp*w_inf - tw(i,j,k)
               end select
            end if
         end do
        end do
       end do
    end do
  end subroutine wall_and_inlet_values

  !---------------------------------------------------------------------
  ! implicit-Euler transport of one cell-centred variable (field id iv,
  ! sources / diffusivity number m of twk):
  !   (phi - phi_n)/dt V + conv(phi) - diff(phi) = Sc V - Sp phi V
  ! convection: contravariant face fluxes, 1st / 2nd-order upwind
  ! (2nd order as deferred correction); diffusion: orthogonal metric
  ! part; SOR sweeps on the linearised system.  Explicit modes evaluate
  ! the term with the old field and sub-cycle to the CFL / VN limits.
  ! The sweeps run over all owned domains between two exchanges with
  ! one global residual per sweep; the sub-cycle count is global.
  !---------------------------------------------------------------------
  subroutine scalar_transport(iv, m)
    integer,  intent(in) :: iv, m
    real(dp), pointer :: f(:,:,:)
    integer  :: ib, d, i, j, k, q, it, nsub, isub, e(3), sgn, ii(3)
    real(dp) :: vol, ff, gf, dts, comax, vnmax, res, res0, resn, dph, sumd, cfl, vn, phi_up, phi_face
    ! ---- sub-cycles for explicit terms (global) --------------------------
    comax = 0.0_dp;  vnmax = 0.0_dp
    do ib = 1, nown
       d = own_list(ib);  call set_block(d)
       allocate(twk(d)%fo(is-ng:ie+ng, js-ng:je+ng, ks-ng:ke+ng))
       allocate(twk(d)%ac(is:ie,js:je,ks:ke), twk(d)%an(6,is:ie,js:je,ks:ke), twk(d)%bb(is:ie,js:je,ks:ke))
       associate(gam => twk(d)%gam)
       do k = ks, ke
        do j = js, je
         do i = is, ie
            if (cstat(i,j,k) /= 1) cycle
            vol = Jh(2*i,2*j,2*k)
            cfl = 0.0_dp;  vn = 0.0_dp
            do q = 1, 3
               e = 0;  e(q) = 1
               cfl = cfl + max(face_flux(q,i,j,k), 0.0_dp) + max(-face_flux(q,i-e(1),j-e(2),k-e(3)), 0.0_dp)
               vn  = vn + 2.0_dp*gam(i,j,k,m)*gql(q,q,2*i,2*j,2*k)
            end do
            comax = max(comax, cfl*dt/vol);  vnmax = max(vnmax, vn*dt)
         end do
        end do
       end do
       end associate
    end do
    comax = global_max(comax);  vnmax = global_max(vnmax)
    nsub = 1
    if (Cnvt_mode_Turb == 1) nsub = max(nsub, ceiling(comax/cfl_switch_Turb))
    if (Dffs_mode_Turb == 1) nsub = max(nsub, ceiling(vnmax/vn_switch_Turb))
    nsub = min(nsub, 200);  nsub_turb = max(nsub_turb, nsub)
    dts = dt/real(nsub, dp)
    do isub = 1, nsub
       ! ---- coefficients (every owned domain) -------------------------
       do ib = 1, nown
          d = own_list(ib);  call set_block(d)
          f => field_ptr(d, iv)
          associate(fo => twk(d)%fo, ac => twk(d)%ac, an => twk(d)%an, bb => twk(d)%bb, &
                    sc => twk(d)%sc, sp => twk(d)%sp, gam => twk(d)%gam)
          fo = f
          do k = ks, ke
           do j = js, je
            do i = is, ie
               ac(i,j,k) = 1.0_dp;  an(:,i,j,k) = 0.0_dp;  bb(i,j,k) = fo(i,j,k)
               if (cstat(i,j,k) /= 1) cycle
               vol = Jh(2*i,2*j,2*k)
               ac(i,j,k) = vol/dts + sp(i,j,k,m)*vol
               bb(i,j,k) = vol/dts*fo(i,j,k) + sc(i,j,k,m)*vol
               do q = 1, 3
                  e = 0;  e(q) = 1
                  do sgn = -1, 1, 2
                     if (sgn > 0) then
                        ff = face_flux(q, i, j, k)                       ! out of the cell through the + face
                        ii = [i,j,k] + e
                        gf = 0.5_dp*(gam(i,j,k,m) + gam(ii(1),ii(2),ii(3),m))*gql(q,q,2*i+1,2*j,2*k)*Jh(2*i+e(1),2*j+e(2),2*k+e(3)) &
                             * merge(1.0_dp, 0.0_dp, q == 1)
                        if (q == 2) gf = 0.5_dp*(gam(i,j,k,m) + gam(ii(1),ii(2),ii(3),m))*gql(q,q,2*i,2*j+1,2*k)*Jh(2*i,2*j+1,2*k)
                        if (q == 3) gf = 0.5_dp*(gam(i,j,k,m) + gam(ii(1),ii(2),ii(3),m))*gql(q,q,2*i,2*j,2*k+1)*Jh(2*i,2*j,2*k+1)
                     else
                        ff = -face_flux(q, i-e(1), j-e(2), k-e(3))        ! out through the - face
                        ii = [i,j,k] - e
                        if (q == 1) gf = 0.5_dp*(gam(i,j,k,m) + gam(ii(1),ii(2),ii(3),m))*gql(q,q,2*i-1,2*j,2*k)*Jh(2*i-1,2*j,2*k)
                        if (q == 2) gf = 0.5_dp*(gam(i,j,k,m) + gam(ii(1),ii(2),ii(3),m))*gql(q,q,2*i,2*j-1,2*k)*Jh(2*i,2*j-1,2*k)
                        if (q == 3) gf = 0.5_dp*(gam(i,j,k,m) + gam(ii(1),ii(2),ii(3),m))*gql(q,q,2*i,2*j,2*k-1)*Jh(2*i,2*j,2*k-1)
                     end if
                     ! diffusion (orthogonal part): gf (phi_nb - phi_c)
                     if (Dffs_mode_Turb == 2) then
                        ac(i,j,k) = ac(i,j,k) + gf
                        an(2*q-1+(sgn+1)/2,i,j,k) = an(2*q-1+(sgn+1)/2,i,j,k) + gf
                     else
                        bb(i,j,k) = bb(i,j,k) + gf*(fo(ii(1),ii(2),ii(3)) - fo(i,j,k))
                     end if
                     ! convection: upwind
                     if (ff >= 0.0_dp) then                                 ! outflow: phi_face = phi_c (1st) + correction
                        phi_up = fo(i,j,k)
                        phi_face = phi_up
                        if (Cnvt_scheme_Turb == 2) phi_face = phi_up + 0.5_dp*minmod2(fo(i,j,k) - fo(ii(1)-2*e(1)*sgn, ii(2)-2*e(2)*sgn, ii(3)-2*e(3)*sgn), fo(ii(1),ii(2),ii(3)) - fo(i,j,k))
                        if (Cnvt_mode_Turb == 2) then
                           ac(i,j,k) = ac(i,j,k) + ff
                           bb(i,j,k) = bb(i,j,k) - ff*(phi_face - phi_up)
                        else
                           bb(i,j,k) = bb(i,j,k) - ff*phi_face
                        end if
                     else                                                   ! inflow from the neighbour
                        phi_up = fo(ii(1),ii(2),ii(3))
                        phi_face = phi_up
                        if (Cnvt_scheme_Turb == 2) phi_face = phi_up + 0.5_dp*minmod2(fo(ii(1),ii(2),ii(3)) - fo(ii(1)+e(1)*sgn, ii(2)+e(2)*sgn, ii(3)+e(3)*sgn), fo(i,j,k) - fo(ii(1),ii(2),ii(3)))
                        if (Cnvt_mode_Turb == 2) then
                           an(2*q-1+(sgn+1)/2,i,j,k) = an(2*q-1+(sgn+1)/2,i,j,k) - ff
                           bb(i,j,k) = bb(i,j,k) - ff*(phi_face - phi_up)
                        else
                           bb(i,j,k) = bb(i,j,k) - ff*phi_face
                        end if
                     end if
                  end do
               end do
            end do
           end do
          end do
          end associate
       end do
       ! ---- SOR sweeps (all owned domains per sweep) --------------------
       res0 = -1.0_dp
       do it = 1, Max_Itr_Turb
          res = 0.0_dp;  sumd = 0.0_dp
          do ib = 1, nown
             d = own_list(ib);  call set_block(d)
             f => field_ptr(d, iv)
             associate(ac => twk(d)%ac, an => twk(d)%an, bb => twk(d)%bb)
             do k = ks, ke
              do j = js, je
               do i = is, ie
                  if (cstat(i,j,k) /= 1) cycle
                  if (wface(i,j,k) /= 0) cycle              ! wall-function cell: Dirichlet value from the wall function
                  dph = bb(i,j,k) + an(1,i,j,k)*f(i-1,j,k) + an(2,i,j,k)*f(i+1,j,k) + an(3,i,j,k)*f(i,j-1,k) &
                        + an(4,i,j,k)*f(i,j+1,k) + an(5,i,j,k)*f(i,j,k-1) + an(6,i,j,k)*f(i,j,k+1)
                  dph = dph/max(ac(i,j,k), 1.0e-30_dp) - f(i,j,k)
                  res = res + dph**2;  sumd = sumd + 1.0_dp
                  f(i,j,k) = f(i,j,k) + SOR_fct_Turb*dph
               end do
              end do
             end do
             end associate
          end do
          resn = sqrt(global_sum(res)/max(global_sum(sumd), 1.0_dp))
          if (res0 < 0.0_dp) res0 = max(resn, 1.0e-30_dp)
          mon_it = it;  mon_res = resn/res0
          call refresh_all(iv)
          if (it >= Min_Itr_Turb .and. resn <= Rsd_Lmt_Turb*res0) exit
       end do
    end do
    do ib = 1, nown
       d = own_list(ib)
       deallocate(twk(d)%fo, twk(d)%ac, twk(d)%an, twk(d)%bb)
    end do
  contains
    pure real(dp) function minmod2(a, b)
      real(dp), intent(in) :: a, b
      if (a*b <= 0.0_dp) then
         minmod2 = 0.0_dp
      else
         minmod2 = sign(min(abs(a), abs(b)), a)
      end if
    end function minmod2
  end subroutine scalar_transport

  ! face volume flux J U^q / fvw of face (q; i,j,k)
  real(dp) function face_flux(q, i, j, k)
    integer, intent(in) :: q, i, j, k
    select case (q)
    case (1);  face_flux = Jh(2*i+1,2*j,2*k)/fvw(2*i+1,2*j,2*k)*(U1(i,j,k) - ugq(1,i,j,k))
    case (2);  face_flux = Jh(2*i,2*j+1,2*k)/fvw(2*i,2*j+1,2*k)*(U2(i,j,k) - ugq(2,i,j,k))
    case default; face_flux = Jh(2*i,2*j,2*k+1)/fvw(2*i,2*j,2*k+1)*(U3(i,j,k) - ugq(3,i,j,k))
    end select
  end function face_flux

  !---------------------------------------------------------------------
  ! wall-function data of every flag-11 cell: <case>_wallfn_NNNNNNNN.dat
  ! (block i j k  x y z  y+  Cf  u_tau/U); serial (rank 0 writes its own
  ! cells; with MPI each rank appends its part)
  !---------------------------------------------------------------------
  subroutine write_wallfn_data(istep)
    integer, intent(in) :: istep
    integer :: iu, ib, d, i, j, k, ax, r
    real(dp) :: uc(3), nrm(3), ut(3), umag, y, utau, yp, nu, mulam, p(3), cf
    character(len=256) :: fname
    if (.not. wallfn_on) return
    nu = 1.0_dp/Re
    write(fname,'(A,A,I8.8,A)') trim(casename), '_wallfn_', istep, '.dat'
    do r = 0, nprocs-1
       if (r == myrank) then
          if (r == 0) then
             open(newunit=iu, file=trim(fname), status='replace', action='write')
             write(iu,'(A)') '# block  i  j  k   x  y  z   y+   Cf   u_tau/U   (wall-function cells, flag 11)'
          else
             open(newunit=iu, file=trim(fname), status='old', position='append', action='write')
          end if
          do ib = 1, nown
             d = own_list(ib);  call set_block(d)
             do k = ks, ke
              do j = js, je
               do i = is, ie
                  if (wface(i,j,k) == 0) cycle
                  ax = abs(wface(i,j,k))
                  uc = cell_velocity(i, j, k);  nrm = face_normal(ax, i, j, k)
                  ut = uc - dot_product(uc, nrm)*nrm;  umag = max(norm2(ut), 1.0e-12_dp)
                  y = wdist(i,j,k);  mulam = mu_lam(i, j, k)
                  ! tau_w = mu_eff u1/y1 = (wfac mu1) u1 / y1  ->  u_tau^2 = wfac u1 / (Re y1 rho)
                  utau = sqrt(max(wfac(i,j,k)*umag/(Re*y*rho_lam(i,j,k)), 0.0_dp))
                  yp = utau*y/(nu*mulam);  cf = 2.0_dp*rho_lam(i,j,k)*utau**2
                  p = cell_center_local(i, j, k)                 ! metres
                  write(iu,'(4I6,3ES14.6,3ES13.5)') dom(d)%b, i, j, k, p, yp, cf, utau
               end do
              end do
             end do
          end do
          close(iu)
       end if
       call MPI_Barrier(MPI_COMM_WORLD, ierr_t)
    end do
  end subroutine write_wallfn_data

  !---------------------------------------------------------------------
  ! u_tau / Uch and y+ of a wall-function (flag-11) cell from its wall
  ! factor (energy.F90: thermal wall function); 0 for other cells
  !---------------------------------------------------------------------
  real(dp) function wall_utau(i, j, k, yp)
    integer,  intent(in)  :: i, j, k
    real(dp), intent(out) :: yp
    integer  :: ax
    real(dp) :: uc(3), nrm(3), ut(3), umag, y, mulam
    wall_utau = 0.0_dp;  yp = 0.0_dp
    if (.not. wallfn_on) return
    if (wface(i,j,k) == 0) return
    ax = abs(wface(i,j,k))
    uc = cell_velocity(i, j, k)
    uc = uc - (1.0_dp - bc_ramp)*[u_init, v_init, w_init]/Uch        ! relative to the ramped wall (update 20)
    if (any_moving) uc = uc - sixdof_wall_velocity(cur_blk, Xh(:,2*i,2*j,2*k), [2*i, 2*j, 2*k])
    nrm = face_normal(ax, i, j, k)
    ut = uc - dot_product(uc, nrm)*nrm;  umag = max(norm2(ut), 1.0e-12_dp)
    y = wdist(i,j,k);  mulam = mu_lam(i, j, k)
    wall_utau = sqrt(max(wfac(i,j,k)*umag/(Re*y*rho_lam(i,j,k)), 0.0_dp))
    yp = wall_utau*y*Re/mulam
  end function wall_utau

  !---------------------------------------------------------------------
  ! monitor line
  !---------------------------------------------------------------------
  subroutine turb_monitor
    if (myrank /= 0 .or. .not. (turb_on .or. wallfn_on)) return
    write(*,'(A)') ' '
    select case (TURB_model)
    case (0);  write(*,'(A)') ' [ Cell ] Center Monitor :    < Turbulence: none / IDDES WMLES >'
    case (1);  write(*,'(A)') ' [ Cell ] Center Monitor :    < Turbulence: Spalart-Allmaras'//merge(' + IDDES >', ' >        ', IDDES_model == 1)
    case (2);  write(*,'(A)') ' [ Cell ] Center Monitor :    < Turbulence: k-epsilon'//merge(' + IDDES >', ' >        ', IDDES_model == 1)
    case (3);  write(*,'(A)') ' [ Cell ] Center Monitor :    < Turbulence: k-omega'//merge(' + IDDES >', ' >        ', IDDES_model == 1)
    case (4);  write(*,'(A)') ' [ Cell ] Center Monitor :    < Turbulence: k-omega SST'//merge(' + IDDES >', ' >        ', IDDES_model == 1)
    end select
    if (TURB_model > 0) then
       write(*,'(A,I15)')    '  --> SOR iterations: ', mon_it
       write(*,'(A,ES19.4)') '  --> residual ratio: ', mon_res
       write(*,'(A,I15)')    '  --> sub-cycles    : ', nsub_turb
    end if
    write(*,'(A,ES19.4)') '  --> max mu_t/mu   : ', mon_mutmax
    if (wallfn_on) write(*,'(A,2F10.2,A,ES11.3)') '  --> wall-function y+ min/max: ', mon_yp_min, mon_yp_max, '   mean Cf: ', mon_cf
    nsub_turb = 1
  end subroutine turb_monitor

end module turbulence
