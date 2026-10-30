!=======================================================================
! energy.F90 -- heat / energy equation for probCFD_부산 (update 16)
!
!   The temperature T [K] is carried DIMENSIONAL at the cell centres
!   (field id 20, ghosted and exchanged like alpha).  One-fluid mixture
!   form, non-conservative, non-dimensional in space / time / velocity
!   (x / Lch, t / tch, u / Uch):
!
!       dT/dt + u . grad T  =  div( Gam grad T ) + S                (1)
!
!       Gam = k_eff / ( (rho c_p) Uch Lch ),  k_eff = k + mu_t c_p / Pr_t
!       (rho c_p) = alpha [ (1 - av) (rho c_p)_l + av (rho c_p)_v ]
!                 + (1 - alpha) (rho c_p)_a               (volume weighted)
!       k         = alpha [ (1 - av) k_l + av k_v ] + (1 - alpha) k_a
!
!   S : latent heat of the cavitation mass transfer, - mdot L / (rho c_p)
!       (mdot > 0 evaporation cools the liquid).  Viscous dissipation and
!       pressure work are neglected (incompressible liquid at low Mach).
!
!   Discretisation (cell centred, contravariant face fluxes of the
!   projected velocity, ALE on moving grids):
!     convection : Cnvt_scheme_Engy 1 1st-order upwind, 2 2nd-order
!                  upwind, 3 QUICK, 4 2nd-order ENO (minmod); the higher-
!                  order part is a deferred correction of the 1st-order
!                  upwind flux.  Written as sum_f F_f (T_f - T_c) - the
!                  NON-CONSERVATIVE form - so that the dilatation of a
!                  cavitating field (div u = S_cav) and the fringe cells
!                  of an overset grid cannot create or destroy heat: the
!                  update is a convex combination of the neighbours.
!     diffusion  : orthogonal part implicit / explicit per direction
!                  (face conductivity = harmonic mean, divided by the
!                  cell's rho c_p), cross-derivative part explicit.
!     time       : implicit Euler per inner iteration (Picard with
!                  Iner_Itn) and the SSP stages of Time_scheme = 2;
!                  Cnvt_mode_Engy / Dffs_mode_Engy 1 explicit (sub-
!                  cycled to cfl_switch_Engy / vn_switch_Engy), 2
!                  implicit, 3 hybrid per direction by the local cell
!                  CFL / Von Neumann number (like Cnvt_mode / Dffs_mode
!                  of the momentum equation).
!     solver     : SOR sweeps (Max_Itr_Engy, Min_Itr_Engy, Rsd_Lmt_Engy
!                  relative to the first sweep, SOR_fct_Engy).
!
!   Boundary conditions (6th line of every .bcs surface, see read_bcs):
!     0 zero gradient (adiabatic wall, outlet, symmetry): no heat flux,
!       convective flux by the cell value
!     1 Dirichlet T_w [K]: implicit 2 Gam_f (T_w - T_c) on the half-cell
!       distance; on a wall-function (flag 11) wall the flux is
!       q_w = rho c_p u_tau (T_w - T_c) / T+(y+, Pr), T+ by Kader's
!       law of the wall for temperature (Pr_t / kappa log slope);
!       an inflow face carries T_w with the mass flux
!     3 heat flux q_w [W/m2] into the fluid (explicit source)
!     4 convective (Robin) h [W/(m2 K)], T_inf [K]: implicit
!       U (T_inf - T_c), U = 1 / (1/h + delta_n / k_f)
!   The ghost cells carry the equivalent mirrored values (output,
!   gradients, upwind stencils); block interfaces / halos / overset
!   fringes are refreshed by the common exchange (refresh_all(20)).
!
!   Buoyancy (Beta_l / Beta_a /= 0 and gravity): Boussinesq acceleration
!   - beta_f (T_f - T_init) g added to the momentum right-hand side
!   (buoyancy_force), beta_f volume-weighted at the velocity face.
!
!   Output: 'Temperature' [K] in the .p3d files; T min / max and the
!   total wall heat in the monitor; Output_data/<case>_heat.dat every
!   Data_Intvl steps with the heat flow [W] into the fluid and the mean
!   wall temperature of every Dirichlet / flux / convective surface.
!=======================================================================
module energy
  use global_vars
  use mpi_vars
  use b2b_runtime, only: refresh_all
  use turbulence, only: mu_lam, rho_lam, wall_utau
  use cavitation, only: rv_cell, lat_cell, dpvdT_cell
  implicit none
  private
  public :: init_energy, energy_start, solve_energy, energy_monitor, energy_export, buoyancy_force, temp_range, &
            energy_finalize, energy_wall_heat, energy_iters, energy_resid

  real(dp), parameter :: kappa = 0.41_dp

  ! ---- physical-boundary faces of the owned domains ------------------
  type :: ebc_t
     integer :: n = 0                                  ! number of boundary faces
     integer, allocatable :: ijk(:,:)                  ! (3,n) cell
     integer, allocatable :: q(:), sd(:), pn(:)        ! face family, side (-1 min / +1 max), .bcs patch
     integer, allocatable :: map(:,:,:,:)              ! (6, is:ie, js:je, ks:ke): record id of the cell faces (0 = none)
     real(dp), allocatable :: qw(:)                    ! heat flow into the fluid through the face [W] (export)
  end type ebc_t
  type(ebc_t), allocatable :: ebc(:)
  ! ---- per-domain work arrays of the transport solver -----------------
  type :: ework_t
     real(dp), allocatable :: fo(:,:,:), gam(:,:,:), ac(:,:,:), an(:,:,:,:), bb(:,:,:)
     integer,  allocatable :: cm(:,:,:,:), dm(:,:,:,:)
  end type ework_t
  type(ework_t), allocatable :: wk(:)
  ! ---- heat surfaces (global list, every rank): block, patch ----------
  integer :: nhs = 0
  integer, allocatable :: hs_blk(:), hs_pat(:)
  integer :: iu_heat = -1
  logical :: heat_open = .false.
  ! ---- monitor ---------------------------------------------------------
  integer  :: mon_it = 0, nsub_e = 1
  real(dp) :: mon_res = 0.0_dp, mon_tmin = 0.0_dp, mon_tmax = 0.0_dp, mon_q = 0.0_dp
  real(dp) :: lat_nd = 0.0_dp                          ! rho_v L / (rho_l c_p,l) [K]

contains

  !=====================================================================
  ! allocation and initial field (after allocate_fields)
  !=====================================================================
  subroutine init_energy
    integer :: ib, d, lo(3), hi(3)
    if (.not. energy_on) return
    do ib = 1, nown
       d = own_list(ib)
       lo = lbound(dom(d)%p);  hi = ubound(dom(d)%p)
       allocate(dom(d)%temp(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3)), dom(d)%tempn(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3)))
       dom(d)%temp = T_init;  dom(d)%tempn = T_init
    end do
    lat_nd = 0.0_dp
    if (cav_on) lat_nd = rho_v*Latent_Heat/(rho1*Cp_l)
    temp_bc_hook => energy_bc
  end subroutine init_energy

  !=====================================================================
  ! boundary-face records, heat-surface list, output file, first ghost
  ! fill (after the metrics / overset assembly / initial field)
  !=====================================================================
  subroutine energy_start
    integer :: ib, d, b, n, ax, sd, lo(3), hi(3), rr(2,3), i, j, k, cnt, pass, idx
    if (.not. energy_on) return
    allocate(ebc(ndom))
    do ib = 1, nown
       d = own_list(ib);  call set_block(d);  b = dom(d)%b
       allocate(ebc(d)%map(6, is:ie, js:je, ks:ke));  ebc(d)%map = 0
       do pass = 1, 2
          cnt = 0
          do n = 1, nbc
             ax = abs(bc_dir(n));  sd = merge(1, -1, bc_dir(n) > 0)
             rr(:,1) = bc_irng(:,n);  rr(:,2) = bc_jrng(:,n);  rr(:,3) = bc_krng(:,n)
             lo = [max(rr(1,1), is), max(rr(1,2), js), max(rr(1,3), ks)]
             hi = [min(rr(2,1)-1, ie), min(rr(2,2)-1, je), min(rr(2,3)-1, ke)]
             if (sd > 0) then
                lo(ax) = 1;  hi(ax) = 1
             else
                lo(ax) = ncell(b, ax);  hi(ax) = lo(ax)
             end if
             if (lo(ax) < dom(d)%r1(ax) .or. lo(ax) > dom(d)%r2(ax)) cycle   ! the domain does not touch this block face
             if (any(hi < lo)) cycle
             do k = lo(3), hi(3)
              do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  cnt = cnt + 1
                  if (pass == 2) then
                     ebc(d)%ijk(:,cnt) = [i, j, k];  ebc(d)%q(cnt) = ax;  ebc(d)%sd(cnt) = -sd;  ebc(d)%pn(cnt) = n
                     ! the boundary face of the cell is on the side opposite to the interior: min face for a min patch
                     idx = 2*ax - 1 + merge(1, 0, sd < 0)
                     ebc(d)%map(idx, i, j, k) = cnt
                  end if
               end do
              end do
             end do
          end do
          if (pass == 1) then
             ebc(d)%n = cnt
             allocate(ebc(d)%ijk(3, max(cnt,1)), ebc(d)%q(max(cnt,1)), ebc(d)%sd(max(cnt,1)), ebc(d)%pn(max(cnt,1)), &
                      ebc(d)%qw(max(cnt,1)))
             ebc(d)%qw = 0.0_dp
          end if
       end do
    end do
    ! heat surfaces: every patch with a T flag 1, 3 or 4 (all blocks, same order on every rank)
    nhs = 0
    do b = 1, nblk
       do n = 1, blk(b)%nbc
          if (blk(b)%bc_flag(6,n) == 1 .or. blk(b)%bc_flag(6,n) == 3 .or. blk(b)%bc_flag(6,n) == 4) nhs = nhs + 1
       end do
    end do
    allocate(hs_blk(max(nhs,1)), hs_pat(max(nhs,1)))
    nhs = 0
    do b = 1, nblk
       do n = 1, blk(b)%nbc
          if (blk(b)%bc_flag(6,n) == 1 .or. blk(b)%bc_flag(6,n) == 3 .or. blk(b)%bc_flag(6,n) == 4) then
             nhs = nhs + 1;  hs_blk(nhs) = b;  hs_pat(nhs) = n
          end if
       end do
    end do
    if (myrank == 0) then
       call make_dir('Output_data')
       open(newunit=iu_heat, file='Output_data/'//trim(casename)//'_heat.dat', status=merge('unknown', 'replace', run_mode == 1), &
            position='append', action='write')
       heat_open = .true.
       if (run_mode /= 1) then
          write(iu_heat,'(A)') '# probCFD energy equation (update 16): temperature range and the heat flow INTO the fluid through'
          write(iu_heat,'(A)') '#   every surface with a T boundary condition of flag 1 (Dirichlet), 3 (heat flux) or 4 (convective)'
          write(iu_heat,'(A,I0)') '# heat surfaces: ', nhs
          do n = 1, nhs
             write(iu_heat,'(A,I0,A,I0,A,I0,A,I0,A,2ES12.4)') '#   surface ', n, ': block ', hs_blk(n), ' patch ', hs_pat(n), &
                '  T flag ', blk(hs_blk(n))%bc_flag(6,hs_pat(n)), '  values ', blk(hs_blk(n))%bc_val(6,hs_pat(n)), &
                blk(hs_blk(n))%bc_tinf(hs_pat(n))
          end do
          write(iu_heat,'(A)', advance='no') '#   step      time[s]      Tmin[K]      Tmax[K]      Tavg[K]'
          do n = 1, nhs
             write(iu_heat,'(A,I0,A,I0,A)', advance='no') '     Q', n, '[W]      Tw', n, '[K]'
          end do
          write(iu_heat,'(A)') ' '
       else
          write(iu_heat,'(A)') '# restart'
       end if
    end if
    call refresh_all(20)
    if (myrank == 0) then
       write(*,'(A,I0,A)') '  energy equation : boundary-face records built, ', nhs, ' heat surface(s) listed in Output_data/'// &
          trim(casename)//'_heat.dat'
    end if
  end subroutine energy_start

  !---------------------------------------------------------------------
  ! mixture properties of a cell (ratios to the liquid): rho c_p, k, and
  ! the effective conductivity k_eff/k_l = k/k_l + mu_t c_p/(k_l Pr_t)
  !---------------------------------------------------------------------
  pure real(dp) function rcp_cell(i, j, k)
    integer, intent(in) :: i, j, k
    real(dp) :: a, av
    a = 1.0_dp;  av = 0.0_dp
    if (Num_of_Phase == 2) a = min(max(alpha(i,j,k), 0.0_dp), 1.0_dp)
    if (cav_on) av = min(max(vap(i,j,k), 0.0_dp), 1.0_dp)
    rcp_cell = a*((1.0_dp - av) + av*rcp_v) + (1.0_dp - a)*rcp_a
  end function rcp_cell

  pure real(dp) function kr_cell(i, j, k)
    integer, intent(in) :: i, j, k
    real(dp) :: a, av
    a = 1.0_dp;  av = 0.0_dp
    if (Num_of_Phase == 2) a = min(max(alpha(i,j,k), 0.0_dp), 1.0_dp)
    if (cav_on) av = min(max(vap(i,j,k), 0.0_dp), 1.0_dp)
    kr_cell = a*((1.0_dp - av) + av*kr_v) + (1.0_dp - a)*kr_a
  end function kr_cell

  real(dp) function kef_cell(i, j, k)
    integer, intent(in) :: i, j, k
    real(dp) :: rho
    kef_cell = kr_cell(i, j, k)
    if (turb_on) then
       rho = rho_lam(i, j, k)
       ! k_t / k_l = (mu_t/mu_l) (mu_l c_p,mix / k_l) / Pr_t,  c_p,mix / c_p,l = rcp / rho
       kef_cell = kef_cell + max(mut(i,j,k), 0.0_dp)*Pr_l/Pr_t*rcp_cell(i,j,k)/max(rho, 1.0e-30_dp)
    end if
  end function kef_cell

  ! molecular Prandtl number of the cell mixture
  real(dp) function pr_cell(i, j, k)
    integer, intent(in) :: i, j, k
    pr_cell = mu_lam(i, j, k)*(rcp_cell(i,j,k)/max(rho_lam(i,j,k), 1.0e-30_dp))/max(kr_cell(i,j,k), 1.0e-30_dp)*Pr_l
  end function pr_cell

  ! face volume flux (relative to the grid) through face (q; i,j,k)
  real(dp) function fflux(q, i, j, k)
    integer, intent(in) :: q, i, j, k
    select case (q)
    case (1);  fflux = Jh(2*i+1,2*j,2*k)/fvw(2*i+1,2*j,2*k)*(U1(i,j,k) - ugq(1,i,j,k))
    case (2);  fflux = Jh(2*i,2*j+1,2*k)/fvw(2*i,2*j+1,2*k)*(U2(i,j,k) - ugq(2,i,j,k))
    case default; fflux = Jh(2*i,2*j,2*k+1)/fvw(2*i,2*j,2*k+1)*(U3(i,j,k) - ugq(3,i,j,k))
    end select
  end function fflux

  !---------------------------------------------------------------------
  ! Kader's law of the wall for the temperature: T+ (y+, Pr)
  !---------------------------------------------------------------------
  pure real(dp) function tplus_kader(yp, pr)
    real(dp), intent(in) :: yp, pr
    real(dp) :: g, bet
    g = 0.01_dp*(pr*yp)**4/(1.0_dp + 5.0_dp*pr**3*yp)
    bet = (3.85_dp*pr**(1.0_dp/3.0_dp) - 1.3_dp)**2 + 2.12_dp*log(pr)
    tplus_kader = pr*yp*exp(-g) + (Pr_t/kappa*log(1.0_dp + yp) + bet)*exp(-1.0_dp/max(g, 1.0e-30_dp))
    tplus_kader = max(tplus_kader, 1.0e-6_dp)
  end function tplus_kader

  !---------------------------------------------------------------------
  ! Boundary face r of the current domain: the diffusion coefficient
  ! cfb (non-dimensional, per unit rho c_p of the cell, so that the
  ! flux into the cell is cfb (T_b - T_c)), the equivalent neighbour
  ! value T_b, and the explicit heat source qsrc (flag 3), all in the
  ! units of the cell equation (K volume / time).  The face "area" is
  ! J |grad xi_q| (per unit computational face).
  !---------------------------------------------------------------------
  subroutine bface_coef(d, r, cfb, tb, qsrc, flag)
    integer,  intent(in)  :: d, r
    real(dp), intent(out) :: cfb, tb, qsrc
    integer,  intent(out) :: flag
    integer  :: i, j, k, q, sd, n, H(3)
    real(dp) :: gq, area, gam, rcp, dn, utau, yp, tp, hh, kface, ueff
    i = ebc(d)%ijk(1,r);  j = ebc(d)%ijk(2,r);  k = ebc(d)%ijk(3,r)
    q = ebc(d)%q(r);  sd = ebc(d)%sd(r);  n = ebc(d)%pn(r)
    flag = bc_flag(6, n)
    cfb = 0.0_dp;  tb = temp(i,j,k);  qsrc = 0.0_dp
    if (flag == 0) return
    H = 2*[i, j, k];  H(q) = H(q) + sd                 ! half index of the boundary face
    gq = gql(q, q, H(1), H(2), H(3))
    area = Jh(H(1),H(2),H(3))*sqrt(max(gq, 0.0_dp))   ! non-dimensional face area
    rcp = rcp_cell(i, j, k)
    kface = kr_cell(i, j, k)                           ! laminar conductivity at the wall / inlet face
    gam = kface/(Pe_l*rcp)                             ! laminar Gam of the cell
    dn = 0.5_dp/max(sqrt(gql(q, q, 2*i, 2*j, 2*k)), 1.0e-30_dp)   ! centre-to-face distance (non-dim)
    select case (flag)
    case (1)
       tb = bc_val(6, n)
       cfb = gam*area/dn                               ! = 2 Gam gql J (half-cell distance)
       if (wallfn_on) then
          if (wface(i,j,k) == -sd*q) then              ! wall-function wall (wface = +ax min patch / -ax max patch): Kader's T+
             utau = wall_utau(i, j, k, yp)
             if (utau > 0.0_dp .and. yp > 0.0_dp) then
                tp = tplus_kader(yp, pr_cell(i, j, k))
                cfb = utau*area/tp
             end if
          end if
       end if
    case (3)
       ! q_w A / ((rho c_p)_c Uch), A non-dimensional: K volume* / time*
       qsrc = bc_val(6, n)*area/(rcp*rho1*Cp_l*Uch)
    case (4)
       tb = bc_tinf(n)
       hh = bc_val(6, n)                               ! [W/(m2 K)]
       ! U = 1/(1/h + delta_n/k_f), delta_n = dn Lch, k_f = kface k_l  ->  non-dim coefficient U A / ((rho c_p)_c Uch)
       ueff = 0.0_dp
       if (hh > 0.0_dp) ueff = 1.0_dp/(1.0_dp/hh + dn*Lch/(kface*k_l))
       cfb = ueff*area/(rcp*rho1*Cp_l*Uch)
    end select
  end subroutine bface_coef

  !---------------------------------------------------------------------
  ! ghost values of T at the physical boundaries of the owned domains
  ! (called by refresh_all(20) after the zero-gradient apply_bc):
  ! Dirichlet -> 2 T_w - T_c (the mirror), flux / convective -> the
  ! value that reproduces the boundary flux with the interior stencil
  !---------------------------------------------------------------------
  subroutine energy_bc
    integer  :: ib, d, r, i, j, k, q, sd, g(3), g2(3), e(3), flag, H(3)
    real(dp) :: cfb, tb, qsrc, gf, tg
    if (.not. energy_on .or. .not. allocated(ebc)) return
    do ib = 1, nown
       d = own_list(ib);  call set_block(d)
       do r = 1, ebc(d)%n
          i = ebc(d)%ijk(1,r);  j = ebc(d)%ijk(2,r);  k = ebc(d)%ijk(3,r)
          q = ebc(d)%q(r);  sd = ebc(d)%sd(r)
          e = 0;  e(q) = sd
          g = [i,j,k] + e;  g2 = g + e
          call bface_coef(d, r, cfb, tb, qsrc, flag)
          select case (flag)
          case (1)
             tg = 2.0_dp*tb - temp(i,j,k)
          case (3, 4)
             ! interior stencil flux through the face: gf (T_g - T_c), gf = Gam gql J
             H = 2*[i,j,k];  H(q) = H(q) + sd
             gf = kr_cell(i,j,k)/(Pe_l*rcp_cell(i,j,k))*gql(q,q,H(1),H(2),H(3))*Jh(H(1),H(2),H(3))
             tg = temp(i,j,k) + (cfb*(tb - temp(i,j,k)) + qsrc)/max(gf, 1.0e-30_dp)
          case default
             tg = temp(i,j,k)
          end select
          temp(g(1),g(2),g(3)) = tg
          if (all(g2 >= lbound(temp)) .and. all(g2 <= ubound(temp))) temp(g2(1),g2(2),g2(3)) = tg
       end do
    end do
  end subroutine energy_bc

  !=====================================================================
  ! one solve of the temperature equation (per inner iteration).  The
  ! SOR sweeps run over ALL owned domains between two exchanges, with
  ! ONE global residual per sweep (a per-domain loop with collectives
  ! inside would deadlock when the ranks own different numbers of
  ! domains); the sub-cycle count is global for the same reason.
  !=====================================================================
  subroutine solve_energy
    integer :: ib, d, it, isub, nsub
    real(dp) :: comax, vnmax, res, sumd, resn, res0, dts
    if (.not. energy_on) return
    if (.not. allocated(wk)) allocate(wk(ndom))
    ! ---- masks and sub-cycles (global) ----------------------------------
    comax = 0.0_dp;  vnmax = 0.0_dp
    do ib = 1, nown
       d = own_list(ib);  call set_block(d)
       call setup_masks(d, comax, vnmax)
    end do
    comax = global_max(comax);  vnmax = global_max(vnmax)
    nsub = 1
    if (Cnvt_mode_Engy == 1) nsub = max(nsub, ceiling(comax/cfl_switch_Engy))
    if (Dffs_mode_Engy == 1) nsub = max(nsub, ceiling(vnmax/vn_switch_Engy))
    nsub = min(nsub, 200);  nsub_e = max(nsub_e, nsub)
    dts = dt/real(nsub, dp)
    do isub = 1, nsub
       do ib = 1, nown
          d = own_list(ib);  call set_block(d)
          call assemble(d, dts)
       end do
       ! ---- SOR sweeps -----------------------------------------------
       res0 = -1.0_dp
       do it = 1, Max_Itr_Engy
          res = 0.0_dp;  sumd = 0.0_dp
          do ib = 1, nown
             d = own_list(ib);  call set_block(d)
             call sweep(d, res, sumd)
          end do
          resn = sqrt(global_sum(res)/max(global_sum(sumd), 1.0_dp))
          if (res0 < 0.0_dp) res0 = max(resn, 1.0e-30_dp)
          mon_it = it;  mon_res = resn/res0
          call refresh_all(20)
          if (it >= Min_Itr_Engy .and. resn <= Rsd_Lmt_Engy*res0) exit
       end do
    end do
    do ib = 1, nown
       d = own_list(ib)
       deallocate(wk(d)%fo, wk(d)%gam, wk(d)%ac, wk(d)%an, wk(d)%bb, wk(d)%cm, wk(d)%dm)
    end do
  end subroutine solve_energy

  !---------------------------------------------------------------------
  ! effective conductivity of every cell (ghosts included), implicit /
  ! explicit masks per cell and direction, explicit CFL / VN maxima
  !---------------------------------------------------------------------
  subroutine setup_masks(d, comax, vnmax)
    integer,  intent(in)    :: d
    real(dp), intent(inout) :: comax, vnmax
    integer  :: i, j, k, q, e(3)
    real(dp) :: vol, rcp, cfl(3), vn(3), cfe, vne
    allocate(wk(d)%fo(is-ng:ie+ng, js-ng:je+ng, ks-ng:ke+ng), wk(d)%gam(is-ng:ie+ng, js-ng:je+ng, ks-ng:ke+ng))
    allocate(wk(d)%ac(is:ie,js:je,ks:ke), wk(d)%an(6,is:ie,js:je,ks:ke), wk(d)%bb(is:ie,js:je,ks:ke))
    allocate(wk(d)%cm(3,is:ie,js:je,ks:ke), wk(d)%dm(3,is:ie,js:je,ks:ke))
    associate(gam => wk(d)%gam, cmask => wk(d)%cm, dmask => wk(d)%dm)
    do k = lbound(gam,3), ubound(gam,3)
     do j = lbound(gam,2), ubound(gam,2)
      do i = lbound(gam,1), ubound(gam,1)
         gam(i,j,k) = kef_cell(i, j, k)
      end do
     end do
    end do
    do k = ks, ke
     do j = js, je
      do i = is, ie
         cmask(:,i,j,k) = 1;  dmask(:,i,j,k) = 1
         if (cstat(i,j,k) /= 1) cycle
         vol = Jh(2*i,2*j,2*k);  rcp = rcp_cell(i, j, k)
         do q = 1, 3
            e = 0;  e(q) = 1
            cfl(q) = dt*max(abs(fflux(q,i,j,k)), abs(fflux(q,i-e(1),j-e(2),k-e(3))))/max(vol, 1.0e-300_dp)
            vn(q)  = 2.0_dp*dt*gam(i,j,k)/(Pe_l*rcp)*gql(q,q,2*i,2*j,2*k)
            select case (Cnvt_mode_Engy)
            case (1);  cmask(q,i,j,k) = 0
            case (2);  cmask(q,i,j,k) = 1
            case default;  cmask(q,i,j,k) = merge(1, 0, cfl(q) >= cfl_switch_Engy)
            end select
            select case (Dffs_mode_Engy)
            case (1);  dmask(q,i,j,k) = 0
            case (2);  dmask(q,i,j,k) = 1
            case default;  dmask(q,i,j,k) = merge(1, 0, vn(q) >= vn_switch_Engy)
            end select
         end do
         cfe = 0.0_dp;  vne = 0.0_dp
         do q = 1, 3
            if (cmask(q,i,j,k) == 0) cfe = cfe + cfl(q)
            if (dmask(q,i,j,k) == 0) vne = vne + vn(q)
         end do
         comax = max(comax, cfe);  vnmax = max(vnmax, vne)
      end do
     end do
    end do
    end associate
  end subroutine setup_masks

  !---------------------------------------------------------------------
  ! coefficients of the linear system of one domain (see the header)
  !---------------------------------------------------------------------
  subroutine assemble(d, dts)
    integer,  intent(in) :: d
    real(dp), intent(in) :: dts
    integer  :: i, j, k, q, l, e(3), el(3), sgn, ii(3), uu(3), H(3), idx, r, flag
    real(dp) :: vol, ff, gf, kf, kc, kn, rcp, src, cfb, tb, qsrc, xs, gq, gl, tfl, dsrc, rvc
    real(dp) :: phi_up, phi_uu, phi_dn, phi_face
    associate(fo => wk(d)%fo, gam => wk(d)%gam, ac => wk(d)%ac, an => wk(d)%an, bb => wk(d)%bb, &
              cmask => wk(d)%cm, dmask => wk(d)%dm)
    fo = temp
    do k = ks, ke
     do j = js, je
      do i = is, ie
         ac(i,j,k) = 1.0_dp;  an(:,i,j,k) = 0.0_dp;  bb(i,j,k) = fo(i,j,k)
         if (cstat(i,j,k) /= 1) cycle
         vol = Jh(2*i,2*j,2*k);  rcp = rcp_cell(i, j, k)
         kc = gam(i,j,k)
         ! sources: latent heat of the cavitation mass transfer, - mdot L / (rho c_p), with the
         ! saturation properties of the cell (update 17) and, when p_v depends on T, its
         ! derivative dS/dT = -dS/dp dp_v/dT taken implicitly (a stiff sink for cryogens)
         src = 0.0_dp;  dsrc = 0.0_dp
         if (cav_on) then
            rvc = rv_cell(i, j, k)
            src = -cavs(i,j,k)/(1.0_dp - rvc)*rvc*lat_cell(i, j, k)/(Cp_l*rcp)
            dsrc = cavdp(i,j,k)*dpvdT_cell(i, j, k)*rvc*lat_cell(i, j, k)/((1.0_dp - rvc)*Cp_l*rcp)   ! <= 0
         end if
         ac(i,j,k) = vol/dts - dsrc*vol
         bb(i,j,k) = vol/dts*fo(i,j,k) + (src - dsrc*fo(i,j,k))*vol
         do q = 1, 3
            e = 0;  e(q) = 1
            do sgn = -1, 1, 2
               idx = 2*q - 1 + (sgn+1)/2                              ! 1- 1+ 2- 2+ 3- 3+
               ii = [i,j,k] + sgn*e                                   ! neighbour across the face
               uu = [i,j,k] - sgn*e                                   ! upstream-upstream of the cell
               H = 2*[i,j,k];  H(q) = H(q) + sgn                      ! face half index
               if (sgn > 0) then
                  ff = fflux(q, i, j, k)                              ! out of the cell through the + face
               else
                  ff = -fflux(q, i-e(1), j-e(2), k-e(3))              ! out through the - face
               end if
               r = ebc(d)%map(idx, i, j, k)
               if (r == 0) then
                  !------------- interior / interface / overset face -----------
                  kn = gam(ii(1),ii(2),ii(3))
                  kf = 2.0_dp*kc*kn/max(kc + kn, 1.0e-30_dp)         ! harmonic mean
                  gq = gql(q, q, H(1), H(2), H(3))
                  gf = kf/(Pe_l*rcp)*gq*Jh(H(1),H(2),H(3))
                  ! diffusion, orthogonal part
                  if (dmask(q,i,j,k) == 1) then
                     ac(i,j,k) = ac(i,j,k) + gf
                     an(idx,i,j,k) = an(idx,i,j,k) + gf
                  else
                     bb(i,j,k) = bb(i,j,k) + gf*(fo(ii(1),ii(2),ii(3)) - fo(i,j,k))
                  end if
                  ! diffusion, cross-derivative part (explicit): kf/(Pe rcp) J sum_l g^{ql} dT/dxi_l at the face
                  xs = 0.0_dp
                  do l = 1, 3
                     if (l == q) cycle
                     gl = gql(q, l, H(1), H(2), H(3))
                     if (gl == 0.0_dp) cycle
                     el = 0;  el(l) = 1
                     tfl = 0.25_dp*( fo(i+el(1),j+el(2),k+el(3)) - fo(i-el(1),j-el(2),k-el(3)) &
                                   + fo(ii(1)+el(1),ii(2)+el(2),ii(3)+el(3)) - fo(ii(1)-el(1),ii(2)-el(2),ii(3)-el(3)) )
                     xs = xs + gl*tfl
                  end do
                  bb(i,j,k) = bb(i,j,k) + real(sgn, dp)*kf/(Pe_l*rcp)*Jh(H(1),H(2),H(3))*xs
                  ! convection, non-conservative: ff (T_face - T_c)
                  if (ff >= 0.0_dp) then                              ! outflow: upwind = this cell
                     phi_up = fo(i,j,k);  phi_uu = fo(uu(1),uu(2),uu(3));  phi_dn = fo(ii(1),ii(2),ii(3))
                     phi_face = face_value(phi_uu, phi_up, phi_dn)
                     bb(i,j,k) = bb(i,j,k) - ff*(phi_face - phi_up)  ! deferred higher-order part (T_c part = 0)
                  else                                                ! inflow: upwind = the neighbour
                     phi_up = fo(ii(1),ii(2),ii(3));  phi_dn = fo(i,j,k)
                     phi_uu = fo(ii(1)+sgn*e(1), ii(2)+sgn*e(2), ii(3)+sgn*e(3))
                     phi_face = face_value(phi_uu, phi_up, phi_dn)
                     if (cmask(q,i,j,k) == 1) then
                        ac(i,j,k) = ac(i,j,k) - ff                    ! |ff| (T_nb - T_c) implicit
                        an(idx,i,j,k) = an(idx,i,j,k) - ff
                        bb(i,j,k) = bb(i,j,k) - ff*(phi_face - phi_up)
                     else
                        bb(i,j,k) = bb(i,j,k) - ff*(phi_face - fo(i,j,k))
                     end if
                  end if
               else
                  !------------- physical boundary face -------------------------
                  call bface_coef(d, r, cfb, tb, qsrc, flag)
                  ! diffusion / Robin: cfb (T_b - T_c), heat flux: qsrc
                  if (dmask(q,i,j,k) == 1) then
                     ac(i,j,k) = ac(i,j,k) + cfb
                     bb(i,j,k) = bb(i,j,k) + cfb*tb + qsrc
                  else
                     bb(i,j,k) = bb(i,j,k) + cfb*(tb - fo(i,j,k)) + qsrc
                  end if
                  ! convection: an inflow face carries T_w (flag 1) or the cell value (0, 3, 4);
                  ! outflow faces carry the cell value -> no contribution in the non-conservative form
                  if (ff < 0.0_dp .and. flag == 1) then
                     if (cmask(q,i,j,k) == 1) then
                        ac(i,j,k) = ac(i,j,k) - ff
                        bb(i,j,k) = bb(i,j,k) - ff*tb
                     else
                        bb(i,j,k) = bb(i,j,k) - ff*(tb - fo(i,j,k))
                     end if
                  end if
               end if
            end do
         end do
      end do
     end do
    end do
    end associate
  contains
    ! face value of the convected variable from the upstream-upstream,
    ! upwind and downwind cell values (Cnvt_scheme_Engy)
    pure real(dp) function face_value(puu, pup, pdn)
      real(dp), intent(in) :: puu, pup, pdn
      select case (Cnvt_scheme_Engy)
      case (2);  face_value = pup + 0.5_dp*(pup - puu)                       ! 2nd-order upwind
      case (3);  face_value = pup + (3.0_dp*pdn - 2.0_dp*pup - puu)/8.0_dp   ! QUICK
      case (4);  face_value = pup + 0.5_dp*minmod2(pup - puu, pdn - pup)     ! 2nd-order ENO (minmod)
      case default;  face_value = pup                                        ! 1st-order upwind
      end select
    end function face_value
    pure real(dp) function minmod2(a, b)
      real(dp), intent(in) :: a, b
      if (a*b <= 0.0_dp) then
         minmod2 = 0.0_dp
      else
         minmod2 = sign(min(abs(a), abs(b)), a)
      end if
    end function minmod2
  end subroutine assemble

  ! one SOR sweep over the active cells of domain d
  subroutine sweep(d, res, sumd)
    integer,  intent(in)    :: d
    real(dp), intent(inout) :: res, sumd
    integer  :: i, j, k
    real(dp) :: dph
    associate(ac => wk(d)%ac, an => wk(d)%an, bb => wk(d)%bb)
    do k = ks, ke
     do j = js, je
      do i = is, ie
         if (cstat(i,j,k) /= 1) cycle
         dph = bb(i,j,k) + an(1,i,j,k)*temp(i-1,j,k) + an(2,i,j,k)*temp(i+1,j,k) + an(3,i,j,k)*temp(i,j-1,k) &
               + an(4,i,j,k)*temp(i,j+1,k) + an(5,i,j,k)*temp(i,j,k-1) + an(6,i,j,k)*temp(i,j,k+1)
         dph = dph/max(ac(i,j,k), 1.0e-30_dp) - temp(i,j,k)
         res = res + dph**2;  sumd = sumd + 1.0_dp
         temp(i,j,k) = temp(i,j,k) + SOR_fct_Engy*dph
      end do
     end do
    end do
    end associate
  end subroutine sweep

  !---------------------------------------------------------------------
  ! Boussinesq buoyancy acceleration of the family-c face at pt (half
  ! index H): - beta_f (T_f - T_init) g_c, T_f and beta_f averaged over
  ! the two cells sharing the face
  !---------------------------------------------------------------------
  real(dp) function buoyancy_force(c, pt, H)
    integer, intent(in) :: c, pt(3), H(3)
    integer  :: e(3)
    real(dp) :: tf, af, bet
    buoyancy_force = 0.0_dp
    if (.not. (energy_on .and. buoy_on)) return
    e = 0;  e(c) = 1
    tf = 0.5_dp*(temp(pt(1),pt(2),pt(3)) + temp(pt(1)+e(1),pt(2)+e(2),pt(3)+e(3)))
    bet = Beta_l
    if (Num_of_Phase == 2) then
       af = alph_h(H(1), H(2), H(3))
       bet = af*Beta_l + (1.0_dp - af)*Beta_a
    end if
    buoyancy_force = -bet*(tf - T_init)*gnd(c)
  end function buoyancy_force

  !---------------------------------------------------------------------
  ! heat flow into the fluid through the heat surfaces [W] and their mean
  ! wall temperature (global sums); the per-face flows are kept for the
  ! monitor
  !---------------------------------------------------------------------
  subroutine energy_wall_heat(qs, tw, qtot)
    real(dp), intent(out) :: qs(:), tw(:), qtot
    integer  :: ib, d, r, i, j, k, q, sd, n, b, m, flag, H(3)
    real(dp) :: cfb, tb, qsrc, qf, rcp, area, twf, gq
    real(dp), allocatable :: asum(:)
    allocate(asum(max(nhs,1)))
    qs = 0.0_dp;  tw = 0.0_dp;  asum = 0.0_dp;  qtot = 0.0_dp
    do ib = 1, nown
       d = own_list(ib);  call set_block(d);  b = dom(d)%b
       do r = 1, ebc(d)%n
          i = ebc(d)%ijk(1,r);  j = ebc(d)%ijk(2,r);  k = ebc(d)%ijk(3,r)
          q = ebc(d)%q(r);  sd = ebc(d)%sd(r);  n = ebc(d)%pn(r)
          ebc(d)%qw(r) = 0.0_dp
          flag = bc_flag(6, n)
          if (flag /= 1 .and. flag /= 3 .and. flag /= 4) cycle
          if (cstat(i,j,k) /= 1) cycle
          call bface_coef(d, r, cfb, tb, qsrc, flag)
          rcp = rcp_cell(i, j, k)
          ! non-dimensional K volume/time -> W: x (rho c_p)_c Uch Lch^2
          qf = (cfb*(tb - temp(i,j,k)) + qsrc)*rcp*rho1*Cp_l*Uch*Lch**2
          ebc(d)%qw(r) = qf
          H = 2*[i,j,k];  H(q) = H(q) + sd
          gq = gql(q, q, H(1), H(2), H(3))
          area = Jh(H(1),H(2),H(3))*sqrt(max(gq, 0.0_dp))*Lch**2
          ! wall temperature: Dirichlet value; flux / convective: from the ghost mirror
          twf = 0.5_dp*(temp(i,j,k) + temp(i+merge(sd,0,q==1), j+merge(sd,0,q==2), k+merge(sd,0,q==3)))
          if (flag == 1) twf = tb
          do m = 1, nhs
             if (hs_blk(m) == b .and. hs_pat(m) == n) then
                qs(m) = qs(m) + qf;  tw(m) = tw(m) + twf*area;  asum(m) = asum(m) + area
                exit
             end if
          end do
       end do
    end do
    do m = 1, nhs
       qs(m) = global_sum(qs(m));  tw(m) = global_sum(tw(m));  asum(m) = global_sum(asum(m))
       if (asum(m) > 0.0_dp) tw(m) = tw(m)/asum(m)
       qtot = qtot + qs(m)
    end do
    deallocate(asum)
  end subroutine energy_wall_heat

  ! min / max / volume average of T over the active cells
  subroutine temp_range(tmin, tmax, tavg)
    real(dp), intent(out) :: tmin, tmax, tavg
    integer  :: ib, d, i, j, k
    real(dp) :: s, v, vol
    tmin = huge(1.0_dp);  tmax = -huge(1.0_dp);  s = 0.0_dp;  v = 0.0_dp
    if (.not. energy_on) then
       tmin = 0.0_dp;  tmax = 0.0_dp;  tavg = 0.0_dp;  return
    end if
    do ib = 1, nown
       d = own_list(ib);  call set_block(d)
       do k = ks, ke
        do j = js, je
         do i = is, ie
            if (cstat(i,j,k) /= 1) cycle
            vol = Jh(2*i,2*j,2*k)
            tmin = min(tmin, temp(i,j,k));  tmax = max(tmax, temp(i,j,k))
            s = s + temp(i,j,k)*vol;  v = v + vol
         end do
        end do
       end do
    end do
    tmin = -global_max(-tmin);  tmax = global_max(tmax)
    s = global_sum(s);  v = global_sum(v)
    tavg = s/max(v, 1.0e-300_dp)
  end subroutine temp_range

  !---------------------------------------------------------------------
  ! export (every Data_Intvl steps) and monitor
  !---------------------------------------------------------------------
  subroutine energy_export(istep, t)
    integer,  intent(in) :: istep
    real(dp), intent(in) :: t
    real(dp), allocatable :: qs(:), tw(:)
    real(dp) :: tmin, tmax, tavg, qtot
    integer :: m
    if (.not. energy_on) return
    if (mod(istep, data_interval) /= 0 .and. istep /= 0) return
    allocate(qs(max(nhs,1)), tw(max(nhs,1)))
    call energy_wall_heat(qs, tw, qtot)
    call temp_range(tmin, tmax, tavg)
    mon_tmin = tmin;  mon_tmax = tmax;  mon_q = qtot
    if (myrank == 0 .and. heat_open) then
       write(iu_heat,'(I8,ES16.8,3ES13.5)', advance='no') istep, t, tmin, tmax, tavg
       do m = 1, nhs
          write(iu_heat,'(2ES13.5)', advance='no') qs(m), tw(m)
       end do
       write(iu_heat,'(A)') ' '
       flush(iu_heat)
    end if
    deallocate(qs, tw)
  end subroutine energy_export

  subroutine energy_monitor
    if (myrank /= 0 .or. .not. energy_on) return
    write(*,'(A)') ' '
    write(*,'(A)') ' [ Cell ] Center Monitor :    < Energy: temperature [K] >'
    write(*,'(A,I15)')    '  --> SOR iterations: ', mon_it
    write(*,'(A,ES19.4)') '  --> residual ratio: ', mon_res
    write(*,'(A,I15)')    '  --> sub-cycles    : ', nsub_e
    write(*,'(A,2F10.3)') '  --> T min / max   : ', mon_tmin, mon_tmax
    write(*,'(A,ES19.4,A)') '  --> wall heat flow: ', mon_q, ' W (into the fluid, all heat surfaces)'
    nsub_e = 1
  end subroutine energy_monitor

  integer function energy_iters()
    energy_iters = mon_it
  end function energy_iters
  real(dp) function energy_resid()
    energy_resid = mon_res
  end function energy_resid

  subroutine energy_finalize
    if (myrank == 0 .and. heat_open) close(iu_heat)
  end subroutine energy_finalize

end module energy
