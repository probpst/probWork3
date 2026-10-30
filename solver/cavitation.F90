!=======================================================================
! cavitation.F90 -- transport-equation cavitation models (update 15)
!
!   A third phase, the VAPOUR, is carried as the volume fraction av of
!   the liquid-vapour part of a cell (the air / water fraction alpha of
!   the two-phase VoF is untouched):
!       rho = alpha [(1 - av) rho_l + av rho_v] + (1 - alpha) rho_a
!   (mu alike).  The vapour fraction obeys
!       d(av)/dt + div(av u) = mdot / rho_v                     (1)
!   and the mixture is no longer divergence-free,
!       div(u) = mdot (1/rho_v - 1/rho_l)                        (2)
!   mdot [kg/(m3 s)] > 0 = evaporation.  (1) is solved by the VoF
!   transport of vof.F90 (CICSAM explicit / implicit, or ENO face values
!   with the time integration of Cnvt_mode; Vapor_scheme) with the mass
!   transfer split into an explicit production and an IMPLICIT
!   destruction (sink proportional to av, so that av stays in [0,1]),
!   and (2) enters the pressure-Poisson equation: S = (mdot/rho_v)
!   (1 - rho_v/rho_l) on the right-hand side (div u* - S)/dt, with its
!   pressure derivative dS/dp (<= 0) taken IMPLICITLY in the pressure
!   (a diagonal term of the Poisson operator, the usual "compressibility"
!   of cavitating projection solvers: without it the explicit volume
!   source makes the pressure overshoot and the cavity flash on and off
!   from step to step).  The slope used is the SECANT towards p_v,
!   S0/(p0 - p_v) - the source vanishes when the pressure reaches p_v -
!   or the tangent of the rate law when that is steeper; after the
!   pressure solve the vapour fraction receives the increment the
!   implicit part produced, dt dS/dp (p - p0)/(1 - rho_v/rho_l), so that
!   the vapour created and the dilatation of the velocity stay
!   consistent (cav_pressure_update).  With Iner_Itn > 1 the sources are
!   re-evaluated from the new pressure in every inner iteration (Picard
!   iteration of the pressure / mass-transfer coupling; 2 recommended).
!   The outflow mass-flux correction includes the volume created.
!   Everything is non-dimensional: rates in 1/tch,
!   pressures by rho_l Uch^2; the SOLVER PRESSURE IS THE ABSOLUTE PRESSURE
!   (the level is set by the Dirichlet p of the .bcs file or RefPrs, and
!   the initial p_cel_initial), so the driving pressure is p - p_v.
!
!   MODELS (Cav_model), rates written as  mdot / rho_v  [1/s]:
!   1 Schnerr-Sauer (bubble dynamics, nuclei density n0, radius R0):
!       R_B = [ (av + a_nuc)/(1 - av) 3/(4 pi n0) ]^(1/3),  a_nuc = V/(1+V),
!       V = n0 4/3 pi R0^3;  C = 3 rho_l / (rho_m R_B) sqrt(2|p - p_v|/(3 rho_l))
!       evaporation (p < p_v):  C (av + a_nuc)(1 - av)
!       condensation (p > p_v): C  av (1 - av)
!   2 Zwart-Gerber-Belamri (F_vap, F_cond, r_nuc, R_b):
!       evaporation:  F_vap  3 r_nuc (1 - av)/R_b sqrt(2 (p_v - p)/(3 rho_l))
!       condensation: F_cond 3 av /R_b          sqrt(2 (p - p_v)/(3 rho_l))
!   3 Kunz (C_dest, C_prod, t_inf = Char_L/Char_U):
!       evaporation:  C_dest (1 - av) (p_v - p)/(1/2 rho_l U^2 t_inf)
!       condensation: C_prod av (1 - av)^2 / t_inf
!   In every model the factor (1 - av) of the evaporation and the factor
!   av of the condensation are taken implicitly; the explicit production
!   is limited to a fraction change of 0.2 per step (the implicit part
!   below supplies the rest).  Fixed numerics: the
!   sources are switched on at Cav_Start_T with a 50-step linear ramp,
!   zero where the air fraction exceeds 0.5 (1 - alpha > 0.5), and the
!   transported av is clipped to [0,1] (the clipped volume is reported).
!   Reboud = 1: the eddy viscosity of the turbulence model uses
!   f(rho) = rho_v + (rho - rho_v)^n/(rho_l - rho_v)^(n-1), n = 10,
!   instead of rho in cavitating cells (turbulence.F90).
!
!   OUTPUT: 'Vapor' after 'VoF' in the .p3d solution files; vapour
!   volume, max av and the evaporating / condensing cell counts in the
!   monitor; sigma and the pressure level in the start-up summary.
!=======================================================================
module cavitation
  use global_vars
  use mpi_vars
  use b2b_runtime, only: refresh_all
  use vof, only: advect_alpha
  implicit none
  private
  public :: init_cavitation, cav_init_field, cav_sources, advect_vapor, cav_divergence_source, cav_pressure_update, &
            vapor_volume, cav_monitor, cav_sigma, pv_sat, rv_cell, lat_cell, dpvdT_cell, cav_div_target, sigma_brennen

  real(dp), parameter :: pi = 3.14159265358979323846_dp
  real(dp), parameter :: AIR_CUT = 0.5_dp            ! no mass transfer where 1 - alpha > AIR_CUT
  integer,  parameter :: NRAMP = 50                  ! source ramp [steps] after Cav_Start_T
  integer :: istep_on = -1                           ! first step with mass transfer
  ! monitor
  real(dp) :: mon_vvap = 0.0_dp, mon_vmax = 0.0_dp, mon_clip = 0.0_dp, mon_src = 0.0_dp
  integer  :: mon_nev = 0, mon_nco = 0, mon_nsub = 0, mon_nit = 0
  real(dp) :: mon_res = 0.0_dp
  real(dp) :: mon_acfl = 0.0_dp, mon_cmin = 0.0_dp, mon_dti = 0.0_dp     ! update 17: acoustic CFL, min sound speed, max delay dT*

contains

  !---------------------------------------------------------------------
  ! allocation of the vapour arrays, pressure level, sigma print
  !---------------------------------------------------------------------
  subroutine init_cavitation
    integer :: ib, d, b, n, lo(3), hi(3)
    logical :: found
    real(dp) :: q, sig
    if (.not. cav_on) return
    do ib = 1, nown
       d = own_list(ib)
       lo = lbound(dom(d)%p);  hi = ubound(dom(d)%p)
       allocate(dom(d)%vap(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3)), dom(d)%vapn(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3)), &
                dom(d)%cavsc(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3)), dom(d)%cavsp(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3)), &
                dom(d)%cavs(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3)), dom(d)%cavdp(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3)), &
                dom(d)%cavp0(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3)), dom(d)%cavbc(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3)), &
                dom(d)%presn(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3)))
       dom(d)%vap = 0.0_dp;  dom(d)%vapn = 0.0_dp;  dom(d)%cavbc = 0.0_dp;  dom(d)%presn = 0.0_dp
       dom(d)%cavsc = 0.0_dp;  dom(d)%cavsp = 0.0_dp;  dom(d)%cavs = 0.0_dp;  dom(d)%cavdp = 0.0_dp;  dom(d)%cavp0 = 0.0_dp
    end do
    ! the absolute pressure level of the case: RefPrs when the level is
    ! pinned, else the (first) Dirichlet pressure of the .bcs file, else
    ! p_cel_initial
    found = .false.;  p_level = p_init
    if (p_ref_switch) then
       p_level = p_ref_val;  found = .true.
    else
       do b = 1, nblk
          do n = 1, size(blk(b)%bc_flag, 2)
             if (blk(b)%bc_flag(4,n) == 1) then
                p_level = blk(b)%bc_val(4,n)*(rho1*Uch*Uch);  found = .true.;  exit
             end if
          end do
          if (found) exit
       end do
    end if
    q = 0.5_dp*rho1*Uch*Uch
    sig = (p_level - Prs_Vapor)/q
    if (myrank == 0) then
       write(*,'(A)') ' '
       write(*,'(A)') ' ---------------- cavitation ----------------'
       if (p_ref_switch) then
          write(*,'(A,ES12.4,A)') '  pressure level   : ', p_level, ' Pa (RefPrs, absolute)'
       else if (found) then
          write(*,'(A,ES12.4,A)') '  pressure level   : ', p_level, ' Pa (Dirichlet p of the .bcs file, absolute)'
       else
          write(*,'(A,ES12.4,A)') '  pressure level   : ', p_level, ' Pa (p_cel_initial; no Dirichlet p and no RefPrs: the level is NOT fixed)'
       end if
       write(*,'(A,ES12.4,A,ES12.4,A)') '  vapour pressure  : ', Prs_Vapor, ' Pa   q = 1/2 rho U^2 = ', q, ' Pa'
       write(*,'(A,F10.4,A)') '  sigma            : ', sig, '  = (p_level - p_v)/q  (at the level; two-phase: + rho g depth/q below the free surface)'
       if (sig <= 0.0_dp) write(*,'(A)') '  WARNING: the pressure level is below the vapour pressure - the whole liquid cavitates'
       if (.not. (p_ref_switch .or. found)) write(*,'(A)') &
          '  WARNING: no Dirichlet pressure and SW_Ref_Prs = 0: the pressure is known up to a constant, p - p_v is meaningless'
    end if
  end subroutine init_cavitation

  real(dp) function cav_sigma()
    cav_sigma = (p_level - Prs_Vapor)/(0.5_dp*rho1*Uch*Uch)
  end function cav_sigma

  ! non-dimensional vapour pressure of the cell: the constant Prs_Vapor,
  ! or, with the energy equation and Pv_of_T = 1, the saturation
  ! pressure at the cell temperature (thermodynamic effect; update 16)
  real(dp) function pv_cell(i, j, k)
    integer, intent(in) :: i, j, k
    pv_cell = pv_nd
    if (energy_on .and. Pv_of_T == 1) pv_cell = pv_T(temp(i,j,k))/(rho1*Uch*Uch)
  end function pv_cell

  ! saturation pressure of water [Pa] at T [K] (Buck 1996, 0.1 Pa-level
  ! accuracy between 0 and 100 C: 2339 Pa at 20 C, 3170 Pa at 25 C)
  pure real(dp) function pv_sat(T)
    real(dp), intent(in) :: T
    real(dp) :: tc
    tc = min(max(T - 273.15_dp, -40.0_dp), 150.0_dp)
    pv_sat = 611.21_dp*exp((18.678_dp - tc/234.5_dp)*(tc/(257.14_dp + tc)))
  end function pv_sat

  ! ---- update 17: saturation properties, thermal delay -------------------
  ! Cav_thermal = 1: the vapour pressure follows the Clausius-Clapeyron law
  ! through the reference point (Prs_Vapor, Lou_v at T_init) with the
  ! constant Latent_Heat - fluid-independent (water, cryogens, ...): the
  ! molar mass M = Lou_v R T_init / Prs_Vapor (ideal gas) and
  !   p_v(T) = Prs_Vapor exp( L M / R (1/T_init - 1/T) ),  rho_v(T) = p_v M / (R T)
  ! Pv_of_T = 1 alone (water): the Buck law pv_sat of 20 C water.
  pure real(dp) function pv_T(T)
    real(dp), intent(in) :: T
    real(dp) :: mm
    if (Cav_thermal == 1) then
       mm = rho_v*8.314_dp*T_init/Prs_Vapor
       pv_T = Prs_Vapor*exp(Latent_Heat*mm/8.314_dp*(1.0_dp/T_init - 1.0_dp/max(T, 1.0_dp)))
    else
       pv_T = pv_sat(T)
    end if
  end function pv_T
  ! vapour density [kg/m3] at saturation (ideal gas at p_v(T))
  pure real(dp) function rhov_T(T)
    real(dp), intent(in) :: T
    real(dp) :: mm
    mm = 0.018015_dp
    if (Cav_thermal == 1) mm = rho_v*8.314_dp*T_init/Prs_Vapor
    rhov_T = pv_T(T)*mm/(8.314_dp*max(T, 1.0_dp))
  end function rhov_T
  ! Brennen's thermodynamic parameter Sigma(T) [m/s^1.5] = rho_v^2 L^2 / (rho_l^2 c_p,l T sqrt(alpha_l))
  pure real(dp) function sigma_brennen(T)
    real(dp), intent(in) :: T
    real(dp) :: al
    al = k_l/(rho1*Cp_l)
    sigma_brennen = (rhov_T(T)*Latent_Heat)**2/(rho1**2*Cp_l*max(T, 1.0_dp)*sqrt(al))
  end function sigma_brennen
  ! vapour density ratio rho_v/rho_l of the cell (saturation value at T with Cav_thermal = 1)
  real(dp) function rv_cell(i, j, k)
    integer, intent(in) :: i, j, k
    rv_cell = rv_ratio
    if (Cav_thermal == 1) rv_cell = rhov_T(temp(i,j,k))/rho1
  end function rv_cell
  ! latent heat of the cell [J/kg]
  real(dp) function lat_cell(i, j, k)
    integer, intent(in) :: i, j, k
    lat_cell = Latent_Heat
  end function lat_cell
  ! d(p_v)/dT of the cell [non-dimensional pressure per K] (0 without the T dependence)
  real(dp) function dpvdT_cell(i, j, k)
    integer, intent(in) :: i, j, k
    real(dp) :: T
    dpvdT_cell = 0.0_dp
    if (.not. (energy_on .and. Pv_of_T == 1)) return
    T = temp(i,j,k)
    dpvdT_cell = (pv_T(T + 0.05_dp) - pv_T(T - 0.05_dp))/0.1_dp/(rho1*Uch*Uch)
  end function dpvdT_cell
  ! compressibility coefficient beta_c = rho_l U^2 sum_k alpha_k/(rho_k c_k^2) of the liquid-vapour part of
  ! the cell (Wood / Wallis mixture law; the air of a two-phase cell stays incompressible), and the mixture
  ! speed of sound [m/s]; 1/(rho c^2) of the vapour = 1/max(p, p_v) (isothermal ideal gas)
  subroutine cav_compress(i, j, k, beta, cmix)
    integer,  intent(in)  :: i, j, k
    real(dp), intent(out) :: beta, cmix
    real(dp) :: av, aw, rv, rhom, inv
    beta = 0.0_dp;  cmix = Snd_liquid
    if (Cav_compressible /= 1) return
    aw = 1.0_dp
    if (Num_of_Phase == 2) aw = min(max(alpha(i,j,k), 0.0_dp), 1.0_dp)
    av = min(max(vap(i,j,k), 0.0_dp), 1.0_dp)
    rv = rv_cell(i, j, k)
    ! the vapour is an isothermal ideal gas: rho_v c_v^2 = p (at least p_v; c_v = sqrt(p_v/rho_v) ~ 370 m/s
    ! for water vapour at 20 C), the liquid has the constant Snd_liquid
    inv = av/max(p(i,j,k)*rho1*Uch*Uch, pv_cell(i,j,k)*rho1*Uch*Uch, 1.0_dp) + (1.0_dp - av)/(rho1*Snd_liquid**2)
    rhom = ((1.0_dp - av) + av*rv)*rho1
    cmix = 1.0_dp/sqrt(rhom*inv)
    beta = aw*rho1*Uch*Uch*inv
  end subroutine cav_compress
  ! divergence target of the cell: the mass-transfer source minus the compression (cavbc = beta_c/dt)
  real(dp) function cav_div_target(i, j, k)
    integer, intent(in) :: i, j, k
    cav_div_target = cavs(i,j,k) - cavbc(i,j,k)*(p(i,j,k) - presn(i,j,k))
  end function cav_div_target

  !---------------------------------------------------------------------
  ! initial vapour fraction of the liquid cells (vap_cel_initial; after the
  ! volume fraction is initialised)
  !---------------------------------------------------------------------
  subroutine cav_init_field
    integer :: ib, d, i, j, k
    if (.not. cav_on .or. vap_init <= 0.0_dp) return
    do ib = 1, nown
       d = own_list(ib);  call set_block(d)
       do k = lbound(vap,3), ubound(vap,3)
        do j = lbound(vap,2), ubound(vap,2)
         do i = lbound(vap,1), ubound(vap,1)
            vap(i,j,k) = vap_init
            if (Num_of_Phase == 2) then
               if (alpha(i,j,k) < 0.5_dp) vap(i,j,k) = 0.0_dp
            end if
         end do
        end do
       end do
       vapn = vap
    end do
    call refresh_all(19)
  end subroutine cav_init_field

  !---------------------------------------------------------------------
  ! Explicit source cavsc and implicit sink coefficient cavsp of the
  ! vapour equation (non-dimensional rates, 1/tch) from the current
  ! pressure and fractions:  d(av)/dt = cavsc - cavsp av
  !---------------------------------------------------------------------
  subroutine cav_sources(istep, t)
    integer,  intent(in) :: istep
    real(dp), intent(in) :: t
    integer :: ib, d, i, j, k, nev, nco
    real(dp) :: av, dp_, aw, rb, rhom, c, anuc, vnuc, ev, co, ramp, s3, dpd, dev, dco, rv, beta, cmix, acfl, cmn, dti, pvl, h2
    real(dp), parameter :: DPMIN = 1.0e-2_dp     ! |p - p_v| floor of the sqrt-rate linearisation [rho U^2]
    real(dp), parameter :: CAPDAV = 0.2_dp       ! maximum explicit vapour-fraction production per step
    if (.not. cav_on) return
    ! ramp: 0 before Cav_Start_T, then 1/NRAMP per step
    if (t < Cav_Start_T - 1.0e-12_dp) then
       ramp = 0.0_dp;  istep_on = -1
    else
       if (istep_on < 0) istep_on = istep
       ramp = min(1.0_dp, real(istep - istep_on + 1, dp)/real(NRAMP, dp))
    end if
    cav_ramp = ramp
    vnuc = SS_Nuclei*4.0_dp/3.0_dp*pi*SS_Radius**3
    anuc = vnuc/(1.0_dp + vnuc)
    nev = 0;  nco = 0;  acfl = 0.0_dp;  cmn = huge(1.0_dp);  dti = 0.0_dp
    do ib = 1, nown
       d = own_list(ib)
       call set_block(d)
       cavsc = 0.0_dp;  cavsp = 0.0_dp;  cavdp = 0.0_dp;  cavp0 = p;  cavbc = 0.0_dp
       ! compressibility of the mixture (update 17): independent of the mass transfer
       if (Cav_compressible == 1) then
          do k = ks, ke
           do j = js, je
            do i = is, ie
               if (cstat(i,j,k) /= 1) cycle
               call cav_compress(i, j, k, beta, cmix)
               ! the vapour volume responds to the pressure only above p_v (below it the phase change
               ! governs, cav_pressure_update clips the gas law at p_v): SECANT of the clipped response
               ! between p^n and the latest pressure iterate (p = p^n in the first inner iteration)
               pvl = pv_cell(i, j, k)
               if (abs(p(i,j,k) - presn(i,j,k)) > 1.0e-12_dp) then
                  beta = beta*(max(p(i,j,k), pvl) - max(presn(i,j,k), pvl))/(p(i,j,k) - presn(i,j,k))
               else if (presn(i,j,k) <= pvl) then
                  beta = 0.0_dp
               end if
               cavbc(i,j,k) = beta/dt
               cmn = min(cmn, cmix)
               ! acoustic CFL: (|u| + c) dt / h along the grid directions
               h2 = max(gql(1,1,2*i,2*j,2*k), gql(2,2,2*i,2*j,2*k), gql(3,3,2*i,2*j,2*k))
               acfl = max(acfl, (cmix/Uch + sqrt((u(i,j,k)**2 + v(i,j,k)**2 + w(i,j,k)**2)))*dt*sqrt(h2))
            end do
           end do
          end do
       end if
       if (ramp <= 0.0_dp) cycle
       do k = ks, ke
        do j = js, je
         do i = is, ie
            if (cstat(i,j,k) /= 1) cycle
            aw = 1.0_dp
            if (Num_of_Phase == 2) aw = alpha(i,j,k)
            if (1.0_dp - aw > AIR_CUT) cycle                      ! air: no mass transfer
            av = min(max(vap(i,j,k), 0.0_dp), 1.0_dp)
            rv = rv_cell(i, j, k)
            pvl = pv_cell(i, j, k)
            if (Cav_thermal == 1) then
               ! thermodynamic delay (B-factor / Plesset-Zwick): the interface of a bubble that has
               ! grown for the time t is cooler than the cell by
               !    dT* = 1/2 (rho_v L / (rho_l c_p,l)) R_B / sqrt(alpha_l t)        (<= 0.1 T)
               ! (the heat of evaporation is conducted through the thermal layer sqrt(alpha t) around
               ! the bubble of radius R_B, the Schnerr-Sauer radius of the cell's vapour fraction), so
               ! the evaporation is driven by p_v(T - dT*)
               if (p(i,j,k) < pvl) then
                  rb = max(((av + anuc)/max(1.0_dp - av, 1.0e-6_dp)*3.0_dp/(4.0_dp*pi*SS_Nuclei))**(1.0_dp/3.0_dp), SS_Radius)
                  dp_ = 0.5_dp*rv*rho1*Latent_Heat/(rho1*Cp_l)*rb/sqrt(k_l/(rho1*Cp_l)*Cav_thermal_t)
                  dp_ = min(dp_, 0.1_dp*temp(i,j,k))
                  dti = max(dti, dp_)
                  pvl = pv_T(temp(i,j,k) - dp_)/(rho1*Uch*Uch)
               end if
            end if
            dp_ = p(i,j,k) - pvl                                   ! (p - p_v)/(rho_l U^2), > 0 condensation
            ev = 0.0_dp;  co = 0.0_dp                              ! evaporation: d(av)/dt = ev (1 - av); condensation: co av
            dev = 0.0_dp;  dco = 0.0_dp                            ! their derivatives with respect to p (non-dimensional)
            select case (Cav_model)
            case (1)                                               ! Schnerr-Sauer
               rhom = (1.0_dp - av) + av*rv
               rb = ((av + anuc)/max(1.0_dp - av, 1.0e-6_dp)*3.0_dp/(4.0_dp*pi*SS_Nuclei))**(1.0_dp/3.0_dp)
               rb = max(rb, SS_Radius)/Lch                         ! non-dimensional bubble radius
               c = 3.0_dp/(rhom*rb)*sqrt(2.0_dp*abs(dp_)/3.0_dp)
               if (dp_ < 0.0_dp) then
                  ev = c*(av + anuc);  dev = -ev/(2.0_dp*max(abs(dp_), DPMIN))
               else
                  co = c*(1.0_dp - av);  dco = co/(2.0_dp*max(abs(dp_), DPMIN))
               end if
            case (2)                                               ! Zwart-Gerber-Belamri
               s3 = sqrt(2.0_dp*abs(dp_)/3.0_dp)
               if (dp_ < 0.0_dp) then
                  ev = ZGB_Fvap*3.0_dp*ZGB_Rnuc/(ZGB_Rb/Lch)*s3;  dev = -ev/(2.0_dp*max(abs(dp_), DPMIN))
               else
                  co = ZGB_Fcond*3.0_dp/(ZGB_Rb/Lch)*s3;  dco = co/(2.0_dp*max(abs(dp_), DPMIN))
               end if
            case default                                           ! Kunz (t_inf = tch)
               dpd = 2.0_dp*dp_                                    ! (p - p_v)/(1/2 rho_l U^2)
               if (dp_ < 0.0_dp) then
                  ev = Kunz_Cdest*(-dpd);  dev = -2.0_dp*Kunz_Cdest
               else
                  co = Kunz_Cprod*(1.0_dp - av)**2
               end if
            end select
            ! the explicit production is limited to a vapour-fraction change of
            ! CAPDAV per step (the bubble-growth rate of the models at a few
            ! hundred Pa of tension fills a cell in a fraction of a millisecond:
            ! with ordinary time steps the pressure-implicit part of the source
            ! below, not the rate law, decides how much vapour forms - it drives
            ! the pressure to p_v); the condensation sink is implicit and unbounded
            ev = min(ev, CAPDAV/(ramp*dt))
            ! explicit production ev (1 - av) -> ev - ev av: the av part implicit
            cavsc(i,j,k) = ramp*ev
            cavsp(i,j,k) = ramp*(ev + co)
            ! pressure-implicit part of the continuity source S = [ev (1 - av)
            ! - co av](1 - rv): the SECANT slope towards p_v, S'(p) = S0/(p0 - p_v)
            ! (<= 0; the source vanishes when the pressure reaches p_v), always
            ! at least the tangent of the rate law; goes to the Poisson diagonal
            cavdp(i,j,k) = ramp*(ev*(1.0_dp - av) - co*av)*(1.0_dp - rv)/(dp_ - sign(DPMIN, dp_))
            if (abs(dp_) > DPMIN) cavdp(i,j,k) = ramp*(ev*(1.0_dp - av) - co*av)*(1.0_dp - rv)/dp_
            cavdp(i,j,k) = min(cavdp(i,j,k), ramp*(dev*(1.0_dp - av) - dco*av)*(1.0_dp - rv), 0.0_dp)
            if (ev > 0.0_dp) nev = nev + 1
            if (co > 0.0_dp .and. av > 1.0e-6_dp) nco = nco + 1
         end do
        end do
       end do
    end do
    mon_nev = nint(global_sum(real(nev, dp)));  mon_nco = nint(global_sum(real(nco, dp)))
    mon_acfl = global_max(acfl);  mon_cmin = -global_max(-cmn);  mon_dti = global_max(dti)
  end subroutine cav_sources

  !---------------------------------------------------------------------
  ! transport of the vapour fraction over one step (from vapn) with the
  ! sources of cav_sources; afterwards the continuity source cavs
  !---------------------------------------------------------------------
  subroutine advect_vapor
    integer :: nsub, nit
    real(dp) :: clipped, resmax, resavg
    if (.not. cav_on) return
    call advect_alpha(nsub, clipped, nit, resmax, resavg, 19)
    mon_nsub = nsub;  mon_nit = nit;  mon_res = resmax;  mon_clip = clipped*Lch**3
    call cav_divergence_source
  end subroutine advect_vapor

  !---------------------------------------------------------------------
  ! S = (mdot/rho_v)(1 - rho_v/rho_l) at the cells from the final
  ! vapour fraction of the step (the same rate the vapour equation used)
  !---------------------------------------------------------------------
  subroutine cav_divergence_source
    integer :: ib, d, i, j, k
    real(dp) :: r, smax
    if (.not. cav_on) return
    smax = 0.0_dp
    do ib = 1, nown
       d = own_list(ib)
       call set_block(d)
       cavs = 0.0_dp
       do k = ks, ke
        do j = js, je
         do i = is, ie
            if (cstat(i,j,k) /= 1) cycle
            r = cavsc(i,j,k) - cavsp(i,j,k)*vap(i,j,k)          ! d(av)/dt by the mass transfer [1/tch]
            cavs(i,j,k) = r*(1.0_dp - rv_cell(i, j, k))
            smax = max(smax, abs(cavs(i,j,k)))
         end do
        end do
       end do
    end do
    mon_src = global_max(smax)
  end subroutine cav_divergence_source

  !---------------------------------------------------------------------
  ! After the pressure solve: the implicit part of the volume source,
  ! dS/dp (p - p0), acted in the projection; the vapour fraction gets the
  ! same increment (dt dS/dp (p - p0)/(1 - rv), clipped to [0,1]) and the
  ! continuity source is updated, so that the vapour created and the
  ! dilatation of the velocity field stay consistent
  !---------------------------------------------------------------------
  subroutine cav_pressure_update
    integer :: ib, d, i, j, k
    real(dp) :: ds, dav, av, pvl, r
    if (.not. cav_on) return
    do ib = 1, nown
       d = own_list(ib)
       call set_block(d)
       do k = ks, ke
        do j = js, je
         do i = is, ie
            if (cstat(i,j,k) /= 1) cycle
            if (Cav_compressible == 1 .and. cavbc(i,j,k) /= 0.0_dp) then
               ! compressible mixture (update 17): the vapour volume follows the isothermal gas law
               ! V_v ~ 1/p between the pressure at the start of the step and the new pressure, both
               ! taken at least p_v (below p_v the mass transfer, not the compression, governs):
               ! av_new = av r / (1 - av + av r), r = p_old/p_new  (bounded; the same compressibility
               ! 1/(rho_v c_v^2) = 1/max(p, p_v) that the Wood law of the pressure equation uses)
               av = min(max(vap(i,j,k), 0.0_dp), 1.0_dp)
               pvl = pv_cell(i, j, k)
               r = max(presn(i,j,k), pvl)/max(p(i,j,k), pvl)
               r = min(max(r, 0.5_dp), 2.0_dp)
               vap(i,j,k) = av*r/max(1.0_dp - av + av*r, 1.0e-30_dp)
            end if
            if (cavdp(i,j,k) == 0.0_dp) cycle
            ds = cavdp(i,j,k)*(p(i,j,k) - cavp0(i,j,k))
            dav = dt*ds/(1.0_dp - rv_cell(i, j, k))
            vap(i,j,k) = min(max(vap(i,j,k) + dav, 0.0_dp), 1.0_dp)
            cavs(i,j,k) = cavs(i,j,k) + ds
         end do
        end do
       end do
    end do
    call refresh_all(19)
  end subroutine cav_pressure_update

  !---------------------------------------------------------------------
  ! vapour volume [non-dimensional] over the active cells
  !---------------------------------------------------------------------
  real(dp) function vapor_volume()
    integer :: ib, d, i, j, k
    real(dp) :: s, vm
    s = 0.0_dp;  vm = 0.0_dp
    if (cav_on) then
       do ib = 1, nown
          d = own_list(ib)
          call set_block(d)
          do k = ks, ke
           do j = js, je
            do i = is, ie
               if (cstat(i,j,k) /= 1) cycle
               s = s + vap(i,j,k)*Jh(2*i,2*j,2*k)*merge(alpha(i,j,k), 1.0_dp, Num_of_Phase == 2)
               vm = max(vm, vap(i,j,k))
            end do
           end do
          end do
       end do
    end if
    vapor_volume = global_sum(s)
    mon_vvap = vapor_volume*Lch**3;  mon_vmax = global_max(vm)
  end function vapor_volume

  !---------------------------------------------------------------------
  ! log block
  !---------------------------------------------------------------------
  subroutine cav_monitor
    if (.not. cav_on .or. myrank /= 0) return
    write(*,'(A)') ' '
    select case (Cav_model)
    case (1);  write(*,'(A)') ' [ Cell ] Center Monitor :    < Cavitation: Schnerr-Sauer >'
    case (2);  write(*,'(A)') ' [ Cell ] Center Monitor :    < Cavitation: Zwart-Gerber-Belamri >'
    case default;  write(*,'(A)') ' [ Cell ] Center Monitor :    < Cavitation: Kunz >'
    end select
    write(*,'(A,ES19.4,A)') '  --> vapour volume : ', mon_vvap, ' m3'
    write(*,'(A,F19.4)')    '  --> max vapour fr.: ', mon_vmax
    write(*,'(A,I9,A,I9)')  '  --> evaporating / condensing cells: ', mon_nev, ' /', mon_nco
    write(*,'(A,ES19.4,A)') '  --> max |div u| src: ', mon_src/tch, ' 1/s'
    write(*,'(A,F19.3)')    '  --> source ramp   : ', cav_ramp
    if (Vapor_scheme == 2 .or. (Vapor_scheme == 3 .and. Cnvt_mode /= 1)) then
       write(*,'(A,I15)')   '  --> TTL # of Itr  : ', mon_nit
       write(*,'(A,ES19.4)')'  --> Max Residual  : ', mon_res
    end if
    write(*,'(A,I15)')      '  --> sub-cycles    : ', mon_nsub
    write(*,'(A,ES19.4,A)') '  --> clipped       : ', mon_clip, ' m3 (this step)'
    if (Cav_compressible == 1) write(*,'(A,F10.3,A,F10.3,A)') '  --> compressible: max acoustic CFL ', mon_acfl, &
       '   min mixture sound speed ', mon_cmin, ' m/s'
    if (Cav_thermal == 1) write(*,'(A,F19.4,A)') '  --> thermal delay dT*: ', mon_dti, ' K (max)'
  end subroutine cav_monitor

end module cavitation
