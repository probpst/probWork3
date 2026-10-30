!=======================================================================
! waves.F90  -  incident waves (update 18)
!
!   Regular and irregular incident waves for the two-phase (VoF) solver,
!   following the CFDShip-Iowa approach (exact potential-flow wave field
!   imposed on the far-field boundaries) with three changes:
!
!   1. The wave is defined ONCE in the earth-fixed frame as a function of
!      the absolute position and time, eta(X,Y,t), u_w(X,Y,Z,t), p_w.
!      The solver's grid coordinates x are converted with one shift,
!        X = x + X_shift(t),  X_shift = -integral of the mean inflow
!      velocity of the wave patches (ramped like every Dirichlet inlet).
!      A body held fixed in a stream U therefore sees the encounter
!      frequency omega + k cos(theta) U automatically, a time-varying
!      stream is handled exactly, and the later moving-background mode
!      only has to move the grid (update 20: the moved grid coordinates
!      ARE the earth positions; the relaxation zones stay attached to
!      the grid faces through the background pose of mod_global).
!   2. The boundary condition is a far-field INCIDENT-WAVE condition
!      (.bcs flag 12 on u, v, w, p and alpha together) that can be put on
!      any face: velocity = mean + wave orbital velocity (Dirichlet,
!      frozen at the free-surface value above the surface), alpha = fill
!      fraction of the ghost cell below eta, pressure ghost = interior +
!      exact wave-pressure difference (zero-gradient-type row of the
!      Poisson operator: the face flux is the prescribed velocity).
!   3. Relaxation zones (waves2Foam style) blend the computed velocity
!      and volume fraction towards the INCIDENT wave over Relax_L(1:6)
!      (x-min x-max y-min y-max z-min z-max; Relax_In_L / Relax_Out_L /
!      Relax_Side_L are aliases of the x-min / x-max / both y faces) of the
!      bounding box of the wave blocks: they absorb the body-generated
!      waves at the inlet (no reflection from the Dirichlet face) and let
!      the incident wave pass out of the outlet.  phi <- (1-w) phi + w
!      phi_wave, w = (exp(chi^3.5) - 1)/(e - 1) (Jacobsen et al. 2012),
!      chi = 1 - distance from the boundary / zone length (w = 1 at the
!      boundary, 0 at the inner edge, small over the inner half of the
!      zone).  The velocity weight is multiplied by the liquid fraction
!      of the face: the air above the surface is left to the solver (the
!      frozen orbital velocity in the air is not divergence-free).
!
!   Wave models (SI inside this module, non-dimensional values returned):
!     SW_Wave = 1 regular wave: Wave_Amp [m], Wave_Period [s], Wave_Phase
!               [deg], Wave_Heading [deg from +x], Water_Depth [m]
!               (0 = deep water), Wave_Order 1 Airy / 2 Stokes 2nd order
!     SW_Wave = 2 irregular long- or short-crested sea from a spectrum:
!               Spectrum 1 Bretschneider (ITTC two-parameter), 2 JONSWAP
!               (Gamma_JONSWAP), 3 Pierson-Moskowitz (Hs only); Hs [m],
!               Tp [s]; N_Wave_Comp components with equal d(omega) over
!               [0.5, 3.5] omega_p, a_l = sqrt(2 S(omega_l) d(omega));
!               N_Wave_Dir > 1: cos^2 spreading over +-90 deg about the
!               heading; random phases from Wave_Seed (own generator,
!               reproducible on every machine); or a component list read
!               from Wave_File ("amp[m] omega[rad/s] heading[deg]
!               phase[deg]" per line).  The components actually used are
!               written to Output_data/<case>_components.wav.
!   Wave_Ramp_T [s]: half-cosine amplitude ramp from the start.
!   Wave_Init = 1: the interior field (alpha, velocity) starts from the
!               exact wave at t = 0 (no ramp).
!   Wave probes: Wave_Probe_N points (Wave_Probe_X/Y [m]); the elevation
!               (column integral of alpha) and the exact incident
!               elevation are written to Output_data/<case>_wave.dat.
!
!   Non-dimensionalisation: positions / Lch, velocities / Uch, pressure
!   / (rho1 Uch^2), elevation / Lch; z is measured from FreeSurf_Z.
!=======================================================================
module waves
  use global_vars
  use mpi_vars
  implicit none
  private
  real(dp), parameter :: pi = 3.14159265358979323846_dp
  public :: wave_on, SW_Wave, wave_nml, wave_setup, wave_start, wave_update, wave_t,                 &
            wave_vel, wave_eta, wave_pdyn, wave_fill, wave_weight, wave_export, wave_monitor,        &
            wave_finalize, wave_umean, Wave_Init, Wave_Ramp_T, wave_ramp,          &
            Relax_L, wave_zones_on, wave_mean_vel, wave_runtime

  ! ---- name list ($Wave) ---------------------------------------------
  integer  :: SW_Wave = 0                 ! 0 off, 1 regular, 2 irregular
  real(dp) :: Wave_Amp = 0.0_dp           ! regular: amplitude [m]
  real(dp) :: Wave_Period = 0.0_dp        ! regular: period [s]
  real(dp) :: Wave_Phase = 0.0_dp         ! regular: phase at X = 0, t = 0 [deg]
  real(dp) :: Wave_Heading = 0.0_dp       ! propagation direction from +x [deg]
  real(dp) :: Water_Depth = 0.0_dp        ! [m] below FreeSurf_Z, 0 = deep water
  integer  :: Wave_Order = 1              ! regular: 1 Airy, 2 Stokes 2nd order
  integer  :: Spectrum = 2                ! irregular: 1 Bretschneider, 2 JONSWAP, 3 PM
  real(dp) :: Hs = 0.0_dp, Tp = 0.0_dp    ! irregular: significant height [m], peak period [s]
  real(dp) :: Gamma_JONSWAP = 3.3_dp
  integer  :: N_Wave_Comp = 100
  integer  :: N_Wave_Dir = 1
  integer  :: Wave_Seed = 1
  character(len=128) :: Wave_File = ' '
  real(dp) :: Wave_Ramp_T = 0.0_dp        ! [s]
  integer  :: Wave_Init = 0
  real(dp) :: Relax_L(6) = 0.0_dp         ! zone lengths [m] from the x-min x-max y-min y-max z-min z-max faces
  integer  :: Wave_Probe_N = 0
  real(dp) :: Wave_Probe_X(8) = 0.0_dp, Wave_Probe_Y(8) = 0.0_dp

  logical  :: wave_on = .false.
  real(dp) :: wave_t = 0.0_dp             ! time [s] the boundary values are evaluated at
  real(dp) :: wave_umean(3) = 0.0_dp      ! mean velocity of the wave patches [m/s]
  logical  :: wave_zones_on = .false.

  ! ---- components (SI) ------------------------------------------------
  integer :: ncmp = 0
  real(dp), allocatable :: wa(:), wk(:), wom(:), wph(:), wcs(:), wsn(:)
  real(dp) :: grv = 9.81_dp               ! |g| [m/s2]
  real(dp) :: depth = 0.0_dp              ! [m], 0 = deep
  real(dp) :: shift_n(3) = 0.0_dp         ! current X_shift / Lch
  real(dp) :: rampf = 1.0_dp              ! current ramp factor
  ! relaxation-zone box (non-dimensional) of the wave blocks
  real(dp) :: box_lo(3), box_hi(3)
  ! probes
  integer  :: iu_wave = -1
  logical  :: wave_open = .false.
  real(dp), allocatable :: probe_eta(:), probe_exact(:)
  integer  :: lcg_state = 12345

contains

  !---------------------------------------------------------------------
  ! name-list hook (called by read_input for unknown keys of group $Wave
  ! and for keys starting with wave_ / relax_)
  !---------------------------------------------------------------------
  subroutine wave_nml(nm, vals, nv, ok)
    character(len=*), intent(in) :: nm, vals(:)
    integer, intent(in) :: nv
    logical, intent(out) :: ok
    integer :: i
    ok = .true.
    select case (nm)
    case ('sw_wave');            SW_Wave = nint(rv(1))
    case ('wave_amp');           Wave_Amp = rv(1)
    case ('wave_height');        Wave_Amp = 0.5_dp*rv(1)
    case ('wave_period');        Wave_Period = rv(1)
    case ('wave_phase');         Wave_Phase = rv(1)
    case ('wave_heading');       Wave_Heading = rv(1)
    case ('water_depth');        Water_Depth = rv(1)
    case ('wave_order');         Wave_Order = nint(rv(1))
    case ('spectrum');           Spectrum = nint(rv(1))
    case ('hs');                 Hs = rv(1)
    case ('tp');                 Tp = rv(1)
    case ('gamma_jonswap');      Gamma_JONSWAP = rv(1)
    case ('n_wave_comp');        N_Wave_Comp = nint(rv(1))
    case ('n_wave_dir');         N_Wave_Dir = nint(rv(1))
    case ('wave_seed');          Wave_Seed = nint(rv(1))
    case ('wave_file');          Wave_File = unq(vals(1))
    case ('wave_ramp_t');        Wave_Ramp_T = rv(1)
    case ('wave_init');          Wave_Init = nint(rv(1))
    case ('relax_l')
       do i = 1, min(nv, 6);  Relax_L(i) = rv(i);  end do
    case ('relax_in_l');         Relax_L(1) = rv(1)
    case ('relax_out_l');        Relax_L(2) = rv(1)
    case ('relax_side_l');       Relax_L(3) = rv(1);  Relax_L(4) = rv(1)
    case ('wave_probe_n');       Wave_Probe_N = min(nint(rv(1)), 8)
    case ('wave_probe_x')
       do i = 1, min(nv, 8);  Wave_Probe_X(i) = rv(i);  end do
       if (Wave_Probe_N == 0) Wave_Probe_N = min(nv, 8)
    case ('wave_probe_y')
       do i = 1, min(nv, 8);  Wave_Probe_Y(i) = rv(i);  end do
    case default
       ok = .false.
    end select
  contains
    real(dp) function rv(n)
      integer, intent(in) :: n
      integer :: io
      rv = 0.0_dp
      if (n > nv) return
      read(vals(n), *, iostat=io) rv
      if (io /= 0) rv = 0.0_dp
    end function rv
    function unq(s) result(tt)
      character(len=*), intent(in) :: s
      character(len=len(s)) :: tt
      integer :: l
      tt = adjustl(s)
      l = len_trim(tt)
      if (l >= 2) then
         if ((tt(1:1) == '''' .and. tt(l:l) == '''') .or. (tt(1:1) == '"' .and. tt(l:l) == '"')) tt = tt(2:l-1)
      end if
    end function unq
  end subroutine wave_nml

  ! runtime-changeable keys (reread_namelist): Wave_Ramp_T and the zone lengths
  subroutine wave_runtime(nm, val, changed)
    character(len=*), intent(in) :: nm
    real(dp), intent(in) :: val
    logical, intent(out) :: changed
    changed = .false.
    select case (nm)
    case ('wave_ramp_t');  if (val /= Wave_Ramp_T .and. val >= 0.0_dp) then; Wave_Ramp_T = val; changed = .true.; end if
    case ('relax_in_l');   if (val /= Relax_L(1) .and. val >= 0.0_dp) then; Relax_L(1) = val; changed = .true.; end if
    case ('relax_out_l');  if (val /= Relax_L(2) .and. val >= 0.0_dp) then; Relax_L(2) = val; changed = .true.; end if
    case ('relax_side_l'); if (val /= Relax_L(3) .and. val >= 0.0_dp) then; Relax_L(3) = val; Relax_L(4) = val; changed = .true.; end if
    end select
    if (changed) wave_zones_on = any(Relax_L > 0.0_dp)
  end subroutine wave_runtime

  !---------------------------------------------------------------------
  ! Set-up after read_input and read_bcs: checks, wave components, mean
  ! velocity of the wave patches, summary.
  !---------------------------------------------------------------------
  subroutine wave_setup
    integer :: b, n, nwp, ierr, l, m, i
    real(dp) :: um(3), om, kk, cs, sn, omp, dom_, wmin, wmax, s, a2, th, dth, dsp, m0, sumk
    wave_on = (SW_Wave /= 0)
    if (.not. wave_on) return
    if (Num_of_Phase /= 2) call werr('SW_Wave needs Num_of_Phase = 2 (free surface)')
    if (grav(1) /= 0.0_dp .or. grav(2) /= 0.0_dp .or. grav(3) >= 0.0_dp) &
       call werr('SW_Wave needs gravity along -z (Grvty = 0, 0, -g)')
    grv = -grav(3)
    depth = max(Water_Depth, 0.0_dp)
    if (Wave_Ramp_T < 0.0_dp) call werr('Wave_Ramp_T must be >= 0')
    if (Wave_Init /= 0) Wave_Ramp_T = 0.0_dp
    ! --- wave patches: count, mean velocity
    nwp = 0;  um = 0.0_dp
    do b = 1, nblk
       do n = 1, blk(b)%nbc
          if (blk(b)%bc_flag(1,n) == 12) then
             nwp = nwp + 1
             if (nwp == 1) then
                um = blk(b)%bc_val(1:3,n)*Uch
             else if (any(abs(blk(b)%bc_val(1:3,n)*Uch - um) > 1.0e-12_dp)) then
                call werr('all flag-12 (wave) patches must carry the same mean velocity on the u, v, w lines')
             end if
          end if
       end do
    end do
    if (nwp == 0) call werr('SW_Wave /= 0 but no surface of the .bcs file has flag 12')
    wave_umean = um
    ! --- components
    select case (SW_Wave)
    case (1)
       if (Wave_Amp <= 0.0_dp .or. Wave_Period <= 0.0_dp) call werr('regular wave: Wave_Amp and Wave_Period must be > 0')
       if (Wave_Order < 1 .or. Wave_Order > 2) call werr('Wave_Order must be 1 (Airy) or 2 (Stokes 2nd order)')
       ncmp = 1
       allocate(wa(1), wk(1), wom(1), wph(1), wcs(1), wsn(1))
       wa(1) = Wave_Amp;  wom(1) = 2.0_dp*pi/Wave_Period
       wk(1) = wavenumber(wom(1))
       wph(1) = Wave_Phase*pi/180.0_dp
       wcs(1) = cos(Wave_Heading*pi/180.0_dp);  wsn(1) = sin(Wave_Heading*pi/180.0_dp)
    case (2)
       Wave_Order = 1
       if (len_trim(Wave_File) > 0) then
          call read_components
       else
          if (Hs <= 0.0_dp) call werr('irregular wave: Hs must be > 0')
          if (Spectrum /= 3 .and. Tp <= 0.0_dp) call werr('irregular wave: Tp must be > 0')
          if (N_Wave_Comp < 1) call werr('N_Wave_Comp must be >= 1')
          if (N_Wave_Dir < 1) call werr('N_Wave_Dir must be >= 1')
          if (Spectrum == 3) then                 ! PM: peak from Hs
             omp = (4.0_dp*0.0324_dp*grv*grv/(Hs*Hs)/5.0_dp)**0.25_dp
             Tp = 2.0_dp*pi/omp
          end if
          omp = 2.0_dp*pi/Tp
          wmin = 0.5_dp*omp;  wmax = 3.5_dp*omp
          dom_ = (wmax - wmin)/real(N_Wave_Comp, dp)
          ncmp = N_Wave_Comp*N_Wave_Dir
          allocate(wa(ncmp), wk(ncmp), wom(ncmp), wph(ncmp), wcs(ncmp), wsn(ncmp))
          lcg_state = max(Wave_Seed, 1)
          dth = pi/real(N_Wave_Dir, dp)           ! spreading band -90 .. +90 deg
          i = 0;  m0 = 0.0_dp
          do l = 1, N_Wave_Comp
             om = wmin + (real(l, dp) - 0.5_dp)*dom_
             s = spectrum_S(om, omp)
             m0 = m0 + s*dom_
             do m = 1, N_Wave_Dir
                i = i + 1
                if (N_Wave_Dir == 1) then
                   th = 0.0_dp;  dsp = 1.0_dp
                else
                   th = -0.5_dp*pi + (real(m, dp) - 0.5_dp)*dth
                   dsp = 2.0_dp/pi*cos(th)**2*dth    ! cos^2 spreading, integral 1
                end if
                a2 = 2.0_dp*s*dom_*dsp
                wa(i) = sqrt(max(a2, 0.0_dp))
                wom(i) = om
                wk(i) = wavenumber(om)
                wph(i) = 2.0_dp*pi*lcg_rand()
                wcs(i) = cos(Wave_Heading*pi/180.0_dp + th)
                wsn(i) = sin(Wave_Heading*pi/180.0_dp + th)
             end do
          end do
          if (myrank == 0) write(*,'(A,F8.4,A,F8.4,A)') '  wave spectrum    : discretised Hs = 4 sqrt(m0) = ', &
             4.0_dp*sqrt(m0), ' m (target ', Hs, ' m)'
       end if
    case default
       call werr('SW_Wave must be 0, 1 or 2')
    end select
    if (any(Relax_L < 0.0_dp)) call werr('Relax_L / Relax_In_L / Relax_Out_L / Relax_Side_L must be >= 0')
    wave_zones_on = any(Relax_L > 0.0_dp)
    if (Wave_Probe_N > 0) allocate(probe_eta(Wave_Probe_N), probe_exact(Wave_Probe_N))
    call wave_update(0.0_dp)                 ! ramp / shift of the initial boundary values
    ! --- summary
    if (myrank == 0) then
       write(*,'(A)') ' '
       select case (SW_Wave)
       case (1)
          write(*,'(A,ES10.3,A,F8.4,A,F8.4,A,F8.4,A)') '  Incident wave    : regular, amplitude ', Wave_Amp, ' m, period ', &
             Wave_Period, ' s, length ', 2.0_dp*pi/wk(1), ' m, celerity ', wom(1)/wk(1), ' m/s'
          write(*,'(A,F7.2,A,F7.2,A,I0,A)') '                     heading ', Wave_Heading, ' deg, phase ', Wave_Phase, &
             ' deg, order ', Wave_Order, merge(' (Stokes 2nd)', ' (Airy)      ', Wave_Order == 2)
          write(*,'(A,F8.4,A,F8.4)') '                     steepness k a = ', wk(1)*Wave_Amp, ', k d = ', &
             merge(wk(1)*depth, 1.0e30_dp, depth > 0.0_dp)
       case (2)
          if (len_trim(Wave_File) > 0) then
             write(*,'(A,I0,A)') '  Incident wave    : irregular, ', ncmp, ' components from '//trim(Wave_File)
          else
             write(*,'(A,I0,A,I0,A,F7.3,A,F7.3,A,I0,A,I0)') '  Incident wave    : irregular, spectrum ', Spectrum, &
                ' (1 Bretschneider, 2 JONSWAP, 3 PM), Hs ', nint(Hs*1000.0_dp), ' mm, Tp ', Tp, ' s, gamma ', Gamma_JONSWAP, &
                ', components ', N_Wave_Comp, ' x directions ', N_Wave_Dir
             write(*,'(A,F7.2,A,I0)') '                     heading ', Wave_Heading, ' deg, seed ', Wave_Seed
          end if
          sumk = 0.0_dp
          do i = 1, ncmp
             sumk = sumk + wa(i)**2
          end do
          write(*,'(A,F8.4,A)') '                     Hs of the component list 4 sqrt(sum a^2/2) = ', 4.0_dp*sqrt(0.5_dp*sumk), ' m'
       end select
       if (depth > 0.0_dp) then
          write(*,'(A,F8.3,A)') '                     water depth ', depth, ' m (finite-depth dispersion)'
       else
          write(*,'(A)') '                     deep water'
       end if
       write(*,'(A,3F8.3,A,I0,A)') '                     mean velocity of the wave patches ', wave_umean, ' m/s (', nwp, ' patches)'
       write(*,'(A,F7.3,A,I0)') '                     ramp ', Wave_Ramp_T, ' s (half cosine), Wave_Init ', Wave_Init
       if (wave_zones_on) then
          write(*,'(A,6F8.3,A)') '                     relaxation zones x-min x-max y-min y-max z-min z-max ', Relax_L, ' m'
       else
          write(*,'(A)') '                     no relaxation zones'
       end if
       if (Wave_Probe_N > 0) write(*,'(A,I0,A)') '                     wave probes: ', Wave_Probe_N, ' (Output_data/'// &
          trim(casename)//'_wave.dat)'
    end if
  end subroutine wave_setup

  ! finite-depth / deep-water wave number of a circular frequency
  real(dp) function wavenumber(om) result(k)
    real(dp), intent(in) :: om
    real(dp) :: f, df
    integer :: it
    k = om*om/grv
    if (depth <= 0.0_dp) return
    do it = 1, 50                              ! Newton on om^2 = g k tanh(k d)
       f  = grv*k*tanh(k*depth) - om*om
       df = grv*(tanh(k*depth) + k*depth/cosh(k*depth)**2)
       if (abs(f) < 1.0e-12_dp*om*om) exit
       k = k - f/df
    end do
  end function wavenumber

  ! spectral density S(omega) [m2 s]
  real(dp) function spectrum_S(om, omp) result(s)
    real(dp), intent(in) :: om, omp
    real(dp) :: sg, r
    select case (Spectrum)
    case (1)                                     ! Bretschneider / ITTC two-parameter
       s = 5.0_dp/16.0_dp*Hs*Hs*omp**4/om**5*exp(-1.25_dp*(omp/om)**4)
    case (2)                                     ! JONSWAP (ITTC form)
       sg = merge(0.07_dp, 0.09_dp, om <= omp)
       r = exp(-(om - omp)**2/(2.0_dp*sg*sg*omp*omp))
       s = 5.0_dp/16.0_dp*(1.0_dp - 0.287_dp*log(Gamma_JONSWAP))*Hs*Hs*omp**4/om**5*exp(-1.25_dp*(omp/om)**4)*Gamma_JONSWAP**r
    case default                                 ! Pierson-Moskowitz (A/om^5 exp(-B/om^4), B from Hs)
       s = 8.1e-3_dp*grv*grv/om**5*exp(-0.0324_dp*grv*grv/(Hs*Hs)/om**4)
    end select
  end function spectrum_S

  ! portable linear congruential generator in (0,1)
  real(dp) function lcg_rand()
    lcg_state = mod(1103515245*lcg_state + 12345, 2147483647)
    if (lcg_state < 0) lcg_state = lcg_state + 2147483647
    lcg_rand = real(lcg_state, dp)/2147483647.0_dp
  end function lcg_rand

  subroutine read_components
    integer :: iu, ios, n, i
    real(dp) :: a, om, th, ph
    character(len=256) :: line
    open(newunit=iu, file=trim(Wave_File), status='old', action='read', iostat=ios)
    if (ios /= 0) call werr('cannot open Wave_File '//trim(Wave_File))
    n = 0
    do
       read(iu, '(A)', iostat=ios) line
       if (ios /= 0) exit
       if (len_trim(line) == 0 .or. line(1:1) == '#' .or. line(1:1) == '!') cycle
       n = n + 1
    end do
    if (n == 0) call werr('Wave_File has no components')
    ncmp = n
    allocate(wa(n), wk(n), wom(n), wph(n), wcs(n), wsn(n))
    rewind(iu);  i = 0
    do
       read(iu, '(A)', iostat=ios) line
       if (ios /= 0) exit
       if (len_trim(line) == 0 .or. line(1:1) == '#' .or. line(1:1) == '!') cycle
       read(line, *, iostat=ios) a, om, th, ph
       if (ios /= 0) call werr('Wave_File: expected "amp omega heading phase" on every line')
       i = i + 1
       wa(i) = a;  wom(i) = om;  wk(i) = wavenumber(om)
       wcs(i) = cos(th*pi/180.0_dp);  wsn(i) = sin(th*pi/180.0_dp);  wph(i) = ph*pi/180.0_dp
    end do
    close(iu)
  end subroutine read_components

  subroutine werr(msg)
    character(len=*), intent(in) :: msg
    integer :: ierr
    if (myrank == 0) write(*,'(A)') ' ERROR ($Wave): '//msg
    flush(6)
    call MPI_Abort(MPI_COMM_WORLD, 71, ierr)
  end subroutine werr

  !---------------------------------------------------------------------
  ! After initialize_fields: relaxation-zone box, probe file, component
  ! file.  (The Wave_Init field itself is set by wave_init_fields in
  ! boundary.F90, which has the metric / exchange routines.)
  !---------------------------------------------------------------------
  subroutine wave_start
    integer :: b, n, i, ierr
    logical :: isw
    real(dp) :: lo(3), hi(3)
    if (.not. wave_on) return
    lo = huge(1.0_dp);  hi = -huge(1.0_dp)
    do b = 1, nblk
       isw = .false.
       do n = 1, blk(b)%nbc
          if (blk(b)%bc_flag(1,n) == 12) isw = .true.
       end do
       if (.not. isw) cycle
       do i = 1, 3
          lo(i) = min(lo(i), minval(blk(b)%xg(i,1:blk(b)%ni,1:blk(b)%nj,1:blk(b)%nk)))
          hi(i) = max(hi(i), maxval(blk(b)%xg(i,1:blk(b)%ni,1:blk(b)%nj,1:blk(b)%nk)))
       end do
    end do
    box_lo = lo;  box_hi = hi
    call wave_update(t_now)
    if (myrank == 0) then
       call make_dir('Output_data')
       if (wave_zones_on) write(*,'(A,3F9.3,A,3F9.3,A)') '  wave blocks box  : ', lo*Lch, '  to ', hi*Lch, ' m (relaxation zones measured from it)'
       open(newunit=iu_wave, file='Output_data/'//trim(casename)//'_components.wav', status='replace', action='write')
       write(iu_wave,'(A)') '# incident-wave components: amplitude [m]  omega [rad/s]  heading [deg from +x]  phase [deg]  (k [1/m])'
       do i = 1, ncmp
          write(iu_wave,'(ES14.6, ES14.6, F10.3, F10.3, ES14.6)') wa(i), wom(i), atan2(wsn(i), wcs(i))*180.0_dp/pi, &
             wph(i)*180.0_dp/pi, wk(i)
       end do
       close(iu_wave)
       if (Wave_Probe_N > 0) then
          open(newunit=iu_wave, file='Output_data/'//trim(casename)//'_wave.dat', status=merge('unknown', 'replace', run_mode == 1), &
               position='append', action='write')
          wave_open = .true.
          if (run_mode /= 1) then
             write(iu_wave,'(A)') '# probCFD incident waves (update 18): free-surface elevation [m] above FreeSurf_Z at the wave probes'
             write(iu_wave,'(A)') '#   eta_i = column integral of alpha over the owned cells whose x-y footprint contains the probe;'
             write(iu_wave,'(A)') '#   exact_i = incident-wave elevation at the probe (ramp included, in the grid frame)'
             write(iu_wave,'(A)', advance='no') '# probes (x, y) [m]:'
             do i = 1, Wave_Probe_N
                write(iu_wave,'(A,F9.4,A,F9.4,A)', advance='no') ' (', Wave_Probe_X(i), ',', Wave_Probe_Y(i), ')'
             end do
             write(iu_wave,'(A)') ' '
             write(iu_wave,'(A)', advance='no') '#   step      time[s]   ramp   X_shift[m]'
             do i = 1, Wave_Probe_N
                write(iu_wave,'(A,I0,A,I0)', advance='no') '       eta_', i, '     exact_', i
             end do
             write(iu_wave,'(A)') ' '
          end if
       end if
    end if
    deallocate(probe_eta, stat=ierr);  deallocate(probe_exact, stat=ierr)
    if (Wave_Probe_N > 0) allocate(probe_eta(Wave_Probe_N), probe_exact(Wave_Probe_N))
  end subroutine wave_start

  !---------------------------------------------------------------------
  ! Per step: time of the boundary values, ramp factor and frame shift.
  !   X_shift(t) = -[u_init t + (U_mean - u_init) R(t)] - extra
  !   R(t) = integral of bc_ramp: t^2/(2 t_ramp) (t < t_ramp), t - t_ramp/2
  !---------------------------------------------------------------------
  subroutine wave_update(tt)
    real(dp), intent(in) :: tt
    real(dp) :: rint, ui(3)
    wave_t = tt
    rampf = wave_ramp(tt)
    if (t_ramp > 0.0_dp) then
       if (tt < t_ramp) then
          rint = tt*tt/(2.0_dp*t_ramp)
       else
          rint = tt - 0.5_dp*t_ramp
       end if
    else
       rint = tt
    end if
    ui = [u_init, v_init, w_init]
    shift_n = -(ui*tt + (wave_umean - ui)*rint)/Lch
  end subroutine wave_update

  real(dp) function wave_ramp(tt)
    real(dp), intent(in) :: tt
    if (Wave_Ramp_T <= 0.0_dp) then
       wave_ramp = 1.0_dp
    else if (tt >= Wave_Ramp_T) then
       wave_ramp = 1.0_dp
    else
       wave_ramp = 0.5_dp*(1.0_dp - cos(pi*max(tt, 0.0_dp)/Wave_Ramp_T))
    end if
  end function wave_ramp

  ! ramped mean velocity (non-dim) of the wave patches, as the Dirichlet inlets
  function wave_mean_vel() result(ub)
    real(dp) :: ub(3)
    ub(1) = u_init/Uch + bc_ramp*(wave_umean(1)/Uch - u_init/Uch)
    ub(2) = v_init/Uch + bc_ramp*(wave_umean(2)/Uch - v_init/Uch)
    ub(3) = w_init/Uch + bc_ramp*(wave_umean(3)/Uch - w_init/Uch)
  end function wave_mean_vel

  !---------------------------------------------------------------------
  ! Wave field at the grid position xn (non-dim) and time tt [s]:
  !   eta  (non-dim, above FreeSurf_Z), uw(3) (non-dim orbital velocity,
  !   frozen at the surface value above the surface), pw (non-dim dynamic
  !   pressure / rho* - multiply by the local density ratio).
  !---------------------------------------------------------------------
  subroutine wave_field(xn, tt, eta, uw, pw, need_uw)
    real(dp), intent(in)  :: xn(3), tt
    real(dp), intent(out) :: eta, uw(3), pw
    logical, intent(in)   :: need_uw
    real(dp) :: X, Y, Z, zz, ps, cs, sn, ek, ch, sh, chp, uh, ww, pp, a, k, om, e2, c2, s2, s4, thd
    integer :: l
    X = (xn(1) + shift_n(1))*Lch;  Y = (xn(2) + shift_n(2))*Lch
    Z = (xn(3) - zfs)*Lch
    eta = 0.0_dp;  uw = 0.0_dp;  pw = 0.0_dp
    if (ncmp == 0) return
    ! elevation first (the velocity above the surface is frozen at eta)
    do l = 1, ncmp
       a = wa(l)*rampf;  k = wk(l)
       ps = k*(X*wcs(l) + Y*wsn(l)) - wom(l)*tt + wph(l)
       eta = eta + a*cos(ps)
       if (Wave_Order == 2 .and. depth > 0.0_dp) then
          eta = eta + a*a*k/4.0_dp*cosh(k*depth)*(2.0_dp + cosh(2.0_dp*k*depth))/sinh(k*depth)**3*cos(2.0_dp*ps)
       else if (Wave_Order == 2) then
          eta = eta + 0.5_dp*a*a*k*cos(2.0_dp*ps)
       end if
    end do
    if (.not. need_uw) then
       eta = eta/Lch
       return
    end if
    zz = min(Z, eta)
    if (depth > 0.0_dp) zz = max(zz, -depth)
    do l = 1, ncmp
       a = wa(l)*rampf;  k = wk(l);  om = wom(l)
       ps = k*(X*wcs(l) + Y*wsn(l)) - om*tt + wph(l)
       cs = cos(ps);  sn = sin(ps)
       if (depth > 0.0_dp) then
          thd = k*depth
          sh = sinh(thd)
          ch  = cosh(k*(zz + depth))/sh
          chp = sinh(k*(zz + depth))/sh
          uh = a*om*ch*cs;  ww = a*om*chp*sn
          pp = grv*a*cosh(k*(zz + depth))/cosh(thd)*cs
          if (Wave_Order == 2) then
             s2 = sh*sh;  s4 = s2*s2
             c2 = cosh(2.0_dp*k*(zz + depth))
             e2 = cos(2.0_dp*ps)
             uh = uh + 0.75_dp*a*a*om*k*c2/s4*e2
             ww = ww + 0.75_dp*a*a*om*k*sinh(2.0_dp*k*(zz + depth))/s4*sin(2.0_dp*ps)
             pp = pp + grv*(0.75_dp*a*a*k*tanh(thd)/s2*(c2/s2 - 1.0_dp/3.0_dp)*e2 - 0.25_dp*a*a*k*tanh(thd)/s2*(c2 - 1.0_dp))
          end if
       else
          ek = exp(k*zz)
          uh = a*om*ek*cs;  ww = a*om*ek*sn
          pp = grv*a*ek*(cs - 0.5_dp*a*ek*k)          ! linear + second-order mean term (as CFDShip-Iowa; the
                                                      ! deep-water limit of the Stokes pressure, no 2nd-order velocity)
       end if
       uw(1) = uw(1) + uh*wcs(l)
       uw(2) = uw(2) + uh*wsn(l)
       uw(3) = uw(3) + ww
       pw = pw + pp
    end do
    eta = eta/Lch
    uw = uw/Uch
    pw = pw/(Uch*Uch)                               ! dynamic pressure / rho1 / Uch^2 (times rho* at the point)
  end subroutine wave_field

  function wave_vel(xn, tt) result(uw)
    real(dp), intent(in) :: xn(3), tt
    real(dp) :: uw(3), eta, pw
    call wave_field(xn, tt, eta, uw, pw, .true.)
  end function wave_vel

  real(dp) function wave_eta(xn, tt)
    real(dp), intent(in) :: xn(3), tt
    real(dp) :: uw(3), pw
    call wave_field(xn, tt, wave_eta, uw, pw, .false.)
  end function wave_eta

  ! non-dim dynamic wave pressure at the half-index position H (density of the point included)
  real(dp) function wave_pdyn(H, tt)
    integer,  intent(in) :: H(3)
    real(dp), intent(in) :: tt
    real(dp) :: eta, uw(3), pw
    integer :: h1, h2, h3
    h1 = min(max(H(1), hlo(1)+1), hhi(1)-1)
    h2 = min(max(H(2), hlo(2)+1), hhi(2)-1)
    h3 = min(max(H(3), hlo(3)+1), hhi(3)-1)
    call wave_field(Xh(:,H(1),H(2),H(3)), tt, eta, uw, pw, .true.)
    wave_pdyn = rho_h(h1, h2, h3)*pw
  end function wave_pdyn

  ! liquid fill fraction of cell (i,j,k) below the incident surface at
  ! the cell centre (x, y) - the ghost / zone target of alpha
  real(dp) function wave_fill(i, j, k, tt)
    integer,  intent(in) :: i, j, k
    real(dp), intent(in) :: tt
    real(dp) :: zmin, zmax, zs
    integer :: a, b, c
    zmin = huge(1.0_dp);  zmax = -huge(1.0_dp)
    do c = 2*k-1, 2*k+1, 2
     do b = 2*j-1, 2*j+1, 2
      do a = 2*i-1, 2*i+1, 2
         zmin = min(zmin, Xh(3,a,b,c))
         zmax = max(zmax, Xh(3,a,b,c))
      end do
     end do
    end do
    zs = zfs + wave_eta(Xh(:,2*i,2*j,2*k), tt)
    if (zmax <= zs) then
       wave_fill = 1.0_dp
    else if (zmin >= zs) then
       wave_fill = 0.0_dp
    else
       wave_fill = (zs - zmin)/max(zmax - zmin, tiny(1.0_dp))
    end if
  end function wave_fill

  !---------------------------------------------------------------------
  ! Relaxation weight (1 at the boundary, 0 inside) at the position xn
  !---------------------------------------------------------------------
  real(dp) function wave_weight(xn)
    real(dp), intent(in) :: xn(3)
    real(dp) :: s, x0(3)
    ! the zones are attached to the grid faces: with a moving background (update 20) the position is
    ! taken back to the original grid coordinates, x0 = Rg^T (x - O0 - dg) + O0
    x0 = xn
    if (bgr_on) x0 = (matmul(transpose(bgr_Rg), xn*Lch - bgr_O0 - bgr_dg) + bgr_O0)/Lch
    s = huge(1.0_dp)
    if (Relax_L(1) > 0.0_dp) s = min(s, (x0(1) - box_lo(1))/(Relax_L(1)/Lch))
    if (Relax_L(2) > 0.0_dp) s = min(s, (box_hi(1) - x0(1))/(Relax_L(2)/Lch))
    if (Relax_L(3) > 0.0_dp) s = min(s, (x0(2) - box_lo(2))/(Relax_L(3)/Lch))
    if (Relax_L(4) > 0.0_dp) s = min(s, (box_hi(2) - x0(2))/(Relax_L(4)/Lch))
    if (Relax_L(5) > 0.0_dp) s = min(s, (x0(3) - box_lo(3))/(Relax_L(5)/Lch))
    if (Relax_L(6) > 0.0_dp) s = min(s, (box_hi(3) - x0(3))/(Relax_L(6)/Lch))
    s = max(s, 0.0_dp)
    if (s >= 1.0_dp) then
       wave_weight = 0.0_dp
    else                                   ! Jacobsen et al. (2012): target weight (exp(chi^3.5) - 1)/(e - 1),
       wave_weight = (exp((1.0_dp - s)**3.5_dp) - 1.0_dp)/(exp(1.0_dp) - 1.0_dp)   ! chi = 1 - s = 1 at the boundary
    end if
  end function wave_weight

  !---------------------------------------------------------------------
  ! Wave probes: elevation from the column integral of alpha
  !---------------------------------------------------------------------
  subroutine wave_export(istep, tt)
    integer,  intent(in) :: istep
    real(dp), intent(in) :: tt
    integer :: ib, d, i, j, k, n, a, b, c
    real(dp) :: zmn, zmx, xp, yp, q(2,4)
    real(dp), allocatable :: s_al(:), z_lo(:), s_alg(:), z_log(:)
    integer :: ierr
    if (.not. wave_on .or. Wave_Probe_N == 0) return
    allocate(s_al(Wave_Probe_N), z_lo(Wave_Probe_N), s_alg(Wave_Probe_N), z_log(Wave_Probe_N))
    s_al = 0.0_dp;  z_lo = huge(1.0_dp)
    do ib = 1, nown
       d = own_list(ib)
       call set_block(d)
       do k = ks, ke
        do j = js, je
         do i = is, ie
            if (cstat(i,j,k) /= 1) cycle
            zmn = huge(1.0_dp);  zmx = -zmn
            do c = 2*k-1, 2*k+1, 2
             do b = 2*j-1, 2*j+1, 2
              do a = 2*i-1, 2*i+1, 2
                 zmn = min(zmn, Xh(3,a,b,c));  zmx = max(zmx, Xh(3,a,b,c))
              end do
             end do
            end do
            ! x-y footprint of the cell = the quadrilateral of its four bottom nodes
            ! (point-in-cell by the inverse bilinear map, half-open in both directions so
            ! that a probe on a shared edge is counted once; update 20: on a yawed grid the
            ! former bounding-box test counted the probe in two neighbouring columns)
            q(:,1) = Xh(1:2, 2*i-1, 2*j-1, 2*k-1);  q(:,2) = Xh(1:2, 2*i+1, 2*j-1, 2*k-1)
            q(:,3) = Xh(1:2, 2*i+1, 2*j+1, 2*k-1);  q(:,4) = Xh(1:2, 2*i-1, 2*j+1, 2*k-1)
            do n = 1, Wave_Probe_N
               xp = Wave_Probe_X(n)/Lch;  yp = Wave_Probe_Y(n)/Lch
               if (xp < minval(q(1,:)) .or. xp > maxval(q(1,:)) .or. yp < minval(q(2,:)) .or. yp > maxval(q(2,:))) cycle
               if (.not. in_quad(q, xp, yp)) cycle
               s_al(n) = s_al(n) + alpha(i,j,k)*(zmx - zmn)
               z_lo(n) = min(z_lo(n), zmn)
            end do
         end do
        end do
       end do
    end do
    call MPI_Allreduce(s_al, s_alg, Wave_Probe_N, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_Allreduce(z_lo, z_log, Wave_Probe_N, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
    do n = 1, Wave_Probe_N
       if (z_log(n) < 0.5_dp*huge(1.0_dp)) then
          probe_eta(n) = (z_log(n) + s_alg(n) - zfs)*Lch
       else
          probe_eta(n) = 0.0_dp
       end if
       probe_exact(n) = wave_eta([Wave_Probe_X(n)/Lch, Wave_Probe_Y(n)/Lch, zfs], tt)*Lch
    end do
    if (myrank == 0 .and. wave_open) then
       write(iu_wave,'(I7, ES13.5, F7.4, ES13.5)', advance='no') istep, tt, rampf, -shift_n(1)*Lch
       do n = 1, Wave_Probe_N
          write(iu_wave,'(2ES13.5)', advance='no') probe_eta(n), probe_exact(n)
       end do
       write(iu_wave,'(A)') ' '
       flush(iu_wave)
    end if
  contains
    ! point (xp, yp) inside the quadrilateral q(:,1..4) (counter-clockwise or clockwise):
    ! inverse bilinear map by Newton iteration, inside = 0 <= xi < 1 and 0 <= et < 1
    logical function in_quad(qq, px, py)
      real(dp), intent(in) :: qq(2,4), px, py
      real(dp) :: xi, et, f(2), jm(2,2), det, dxi, det_, r(2)
      integer :: it
      xi = 0.5_dp;  et = 0.5_dp
      do it = 1, 8
         r = (1-xi)*(1-et)*qq(:,1) + xi*(1-et)*qq(:,2) + xi*et*qq(:,3) + (1-xi)*et*qq(:,4)
         f = r - [px, py]
         jm(:,1) = -(1-et)*qq(:,1) + (1-et)*qq(:,2) + et*qq(:,3) - et*qq(:,4)      ! d r / d xi
         jm(:,2) = -(1-xi)*qq(:,1) - xi*qq(:,2) + xi*qq(:,3) + (1-xi)*qq(:,4)      ! d r / d et
         det = jm(1,1)*jm(2,2) - jm(1,2)*jm(2,1)
         if (abs(det) < 1.0e-300_dp) exit
         dxi  = ( jm(2,2)*f(1) - jm(1,2)*f(2))/det
         det_ = (-jm(2,1)*f(1) + jm(1,1)*f(2))/det
         xi = xi - dxi;  et = et - det_
         if (abs(dxi) + abs(det_) < 1.0e-12_dp) exit
      end do
      ! half-open with a round-off margin: a probe exactly on a node or an edge belongs to the
      ! cell on its + side only (xi, et in [-1e-9, 1 - 1e-9))
      in_quad = (xi >= -1.0e-9_dp .and. xi < 1.0_dp - 1.0e-9_dp .and. et >= -1.0e-9_dp .and. et < 1.0_dp - 1.0e-9_dp)
    end function in_quad
  end subroutine wave_export

  subroutine wave_monitor
    integer :: n
    if (.not. wave_on .or. myrank /= 0) return
    write(*,'(A)') ' '
    write(*,'(A,F7.4,A,3ES11.3,A)') ' [ Wave ] ramp ', rampf, '   frame shift X_shift [m] ', -shift_n*Lch, '   (earth X = x - X_shift)'
    if (Wave_Probe_N > 0) then
       write(*,'(A)', advance='no') '  --> probe elevation / exact [m]:'
       do n = 1, Wave_Probe_N
          write(*,'(A,ES10.3,A,ES10.3)', advance='no') '  ', probe_eta(n), ' /', probe_exact(n)
       end do
       write(*,'(A)') ' '
    end if
  end subroutine wave_monitor

  subroutine wave_finalize
    if (wave_open) close(iu_wave)
  end subroutine wave_finalize

end module waves
