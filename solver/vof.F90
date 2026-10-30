!=======================================================================
! vof.F90
!   Algebraic Volume-of-Fluid transport of the liquid volume fraction
!   alpha (1 = liquid / fluid 1, 0 = gas / fluid 2) on the fully
!   staggered curvilinear grid, one-fluid formulation.
!
!   Transport equation (divergence-free velocity):
!        d(alpha)/dt + div(alpha u) = 0
!   Finite-volume update per cell with the SAME face volume fluxes
!   F^q = J U^q that make the discrete continuity equation telescope
!   (poisson.F90), so the scheme is conservative and, with a converged
!   pressure solve, mass errors stay at the Poisson residual level.
!
!   Face values: CICSAM (Ubbink & Issa, JCP 153, 1999) - the compressive
!   HYPER-C bound and ULTIMATE-QUICKEST are blended by the angle between
!   the interface normal (grad alpha at the donor cell, from the half-
!   lattice metrics) and the face normal grad(xi_q).  Time integration
!   (vof_scheme): 1 = explicit, 2 = implicit with the CICSAM weights
!   lagged and weighted-Jacobi sweeps to a true-residual criterion
!   (Rsd_Lmt_VoF); both with automatic sub-cycling so that the donor-
!   cell Courant number stays below vof_cfl_max (namelist).  A final
!   clipping to [0,1] guards against multi-dimensional overshoots; the
!   clipped volume is reported in the step monitor.
!
!   vof_scheme = 3: GEOMETRIC VoF (PLIC).  The interface is reconstructed
!   in every mixed cell as a plane in the cell's computational unit cube
!   (Youngs normal from the 3x3x3 index stencil, plane constant from the
!   analytical volume relation of Scardovelli & Zaleski, JCP 164, 2000)
!   and advected by direction-split sweeps (order rotated every sub-
!   step) with the exact liquid volume of the flux slab through every
!   face and the divergence correction of Weymouth & Yue (JCP 229, 2010),
!   so that the three sweeps together are conservative for the discretely
!   divergence-free face fluxes and bounded for face Courant numbers
!   <= 0.5 (enforced by the sub-cycling).  Ghost / interface / overset
!   fringe values are refreshed after every sweep.
!
!   Initial condition / inlet profile: the liquid fills z < zfs
!   (z_free_surface / Lch); a cell cut by the plane gets the fraction of
!   its node z-extent that lies below the plane.
!=======================================================================
module vof
  use global_vars
  use mpi_vars
  use b2b_runtime, only: refresh_all
  use overset, only: ovs_overlap_weight
  use ovs_geom, only: plic_cube_frac, plic_plane_const, plic_box_frac, plic_normal27
  implicit none
  private
  real(dp), pointer :: avf(:,:,:) => null()   ! the fraction being advected: alpha (id 8) or the vapour fraction (id 19)
  integer :: cur_iv = 8
  public :: fill_fraction, init_alpha, advect_alpha, liquid_volume,     &
            alpha_range
  logical :: mirrored(3) = .false.        ! axes mirrored by the last plic_plane call

contains

  !---------------------------------------------------------------------
  ! Liquid fraction of cell (i,j,k) of the CURRENT domain for the still-
  ! water initial state: fraction of the cell's node z-range below zfs.
  ! Valid for every cell whose 8 nodes lie inside the half lattice
  ! (all cells of the ghost-extended storage).
  !---------------------------------------------------------------------
  ! set_block + the advected-fraction pointer (advect_alpha works on
  ! alpha or, for the cavitation model, on the vapour fraction)
  subroutine sbv(d)
    integer, intent(in) :: d
    call set_block(d)
    if (cur_iv == 19) then
       avf => vap
    else
       avf => alpha
    end if
  end subroutine sbv

  pure real(dp) function fill_fraction(i, j, k)
    integer, intent(in) :: i, j, k
    real(dp) :: zmin, zmax
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
    if (zmax <= zfs) then
       fill_fraction = 1.0_dp
    else if (zmin >= zfs) then
       fill_fraction = 0.0_dp
    else
       fill_fraction = (zfs - zmin) / max(zmax - zmin, tiny(1.0_dp))
    end if
  end function fill_fraction

  !---------------------------------------------------------------------
  ! Liquid fraction of cell (i,j,k) for a liquid DROP in gas (Init_shape
  ! 1 = circular cylinder along y, 2 = sphere, 3 = square prism along y
  ! with half side Drop_Radius), centre Drop_Center; 5x5x5 sub-samples of
  ! the trilinear cell map.
  !---------------------------------------------------------------------
  pure real(dp) function drop_fraction(i, j, k) result(fr)
    integer, intent(in) :: i, j, k
    integer  :: a, b, c, t, ins
    real(dp) :: s(3), wn(0:7), x(3), dxc(3), rr, cen(3)
    integer, parameter :: ns = 5
    cen = Drop_Center/Lch;  rr = Drop_Radius/Lch
    ins = 0
    do c = 1, ns
     do b = 1, ns
      do a = 1, ns
         s = ([a, b, c] - 0.5_dp)/real(ns, dp)
         x = 0.0_dp
         do t = 0, 7
            wn(t) = merge(s(1), 1.0_dp-s(1), iand(t,1) == 1)*merge(s(2), 1.0_dp-s(2), iand(t/2,1) == 1)* &
                    merge(s(3), 1.0_dp-s(3), iand(t/4,1) == 1)
            x = x + wn(t)*Xh(:, 2*i-1+2*iand(t,1), 2*j-1+2*iand(t/2,1), 2*k-1+2*iand(t/4,1))
         end do
         dxc = x - cen
         select case (Init_shape)
         case (1);  if (dxc(1)**2 + dxc(3)**2 <= rr**2) ins = ins + 1
         case (2);  if (sum(dxc**2) <= rr**2) ins = ins + 1
         case default;  if (abs(dxc(1)) <= rr .and. abs(dxc(3)) <= rr) ins = ins + 1
         end select
      end do
     end do
    end do
    fr = real(ins, dp)/real(ns**3, dp)
  end function drop_fraction

  !---------------------------------------------------------------------
  ! Initial volume fraction on every owned domain (ghost cells included,
  ! then refreshed through exchange + BCs).
  !---------------------------------------------------------------------
  subroutine init_alpha
    integer :: ib, d, i, j, k
    do ib = 1, nown
       d = own_list(ib)
       call set_block(d)
       do k = lbound(alpha,3), ubound(alpha,3)
        do j = lbound(alpha,2), ubound(alpha,2)
         do i = lbound(alpha,1), ubound(alpha,1)
            if (Init_shape == 0) then
               alpha(i,j,k) = fill_fraction(i, j, k)
               ! a_cel_initial: volume fraction of the cells above FreeSurf_Z
               alpha(i,j,k) = alpha(i,j,k) + (1.0_dp - alpha(i,j,k))*alpha_init
            else
               alpha(i,j,k) = drop_fraction(i, j, k)
            end if
         end do
        end do
       end do
    end do
    call refresh_all(8)
  end subroutine init_alpha

  !---------------------------------------------------------------------
  ! Advect alpha over one time step dt with the current (divergence-
  ! free) face velocities U1..U3, starting from the cell values held in
  ! alpha on entry.  Both schemes sub-cycle so that the donor-cell
  ! Courant number stays below vof_cfl_max (nsub sub-steps):
  !   vof_scheme = 1 : explicit CICSAM (face values from the start-of-
  !                    sub-step field), nit = 0, residuals = 0;
  !   vof_scheme = 2 : implicit CICSAM.  The face value is written as
  !                    alpha_f = (1 - beta_f) alpha_D + beta_f alpha_A
  !                    with the CICSAM weight beta_f (0 = upwind, 1 =
  !                    downwind) lagged from the previous iterate, and
  !                    the linear system
  !                      alpha_P + (dt/V_P) sum_f F_f alpha_f = alpha_P^n
  !                    is solved by weighted Jacobi sweeps (SOR_fct_VoF)
  !                    until the TRUE ALGEBRAIC RESIDUAL of that scaled
  !                    equation, max_P |alpha^n - alpha - dt/V sum F alpha_f|,
  !                    drops below Rsd_Lmt_VoF (Min/Max_Itr_VoF honoured).
  ! Returns the sub-cycle count, the volume removed by the final [0,1]
  ! clipping (non-dimensional, all ranks), the number of sweeps (summed
  ! over the sub-steps) and the max / average residual of the last sweep.
  ! Work storage: rhsu/rhsv/rhsw hold the face values (scheme 1) or the
  ! weights beta (scheme 2); rhsp holds alpha at the start of the sub-
  ! step; wrk the new iterate (all rebuilt by the momentum/Poisson
  ! stages afterwards).
  !---------------------------------------------------------------------
  subroutine advect_alpha(nsub, clipped, nit, resmax, resavg, iv)
    integer,  intent(out) :: nsub, nit
    real(dp), intent(out) :: clipped, resmax, resavg
    integer,  intent(in), optional :: iv      ! field id: 8 alpha (default), 19 vapour fraction (cavitation, with sources)
    integer  :: sch
    logical  :: eno, with_src
    integer  :: ib, d, i, j, k, isub, it, ncell_tot
    integer, parameter :: n_beta = 8
    real(dp) :: comax, co, dts, fe, fw, fn, fs, ft, fb, anew, cl, vol
    real(dp) :: diag, off, rr, rsum, rmax, ncl

    cur_iv = 8
    if (present(iv)) cur_iv = iv
    with_src = (cur_iv == 19);  eno = .false.
    if (cur_iv == 8) then
       sch = vof_scheme
    else
       ! vapour (update 15): 1 explicit CICSAM, 2 implicit CICSAM, 3 ENO
       ! (minmod) face values with the time integration of Cnvt_mode
       ! (1 explicit sub-cycled, 2 implicit, 3 hybrid: implicit when the
       ! Courant number exceeds cfl_switch)
       sch = min(Vapor_scheme, 2)
       if (Vapor_scheme == 3) then
          eno = .true.
          sch = merge(1, 2, Cnvt_mode == 1)
       end if
    end if
    if (sch == 3) then
       call advect_alpha_plic(nsub, clipped)
       nit = 0;  resmax = 0.0_dp;  resavg = 0.0_dp
       return
    end if
    !---- sub-cycle count from the maximum donor-cell outflow Courant ---
    comax = 0.0_dp
    do ib = 1, nown
       d = own_list(ib)
       call sbv(d)
       do k = ks, ke
        do j = js, je
         do i = is, ie
            if (cstat(i,j,k) /= 1) cycle
            co = ( max(flx(1,i,j,k),0.0_dp)   + max(-flx(1,i-1,j,k),0.0_dp)  &
                 + max(flx(2,i,j,k),0.0_dp)   + max(-flx(2,i,j-1,k),0.0_dp)  &
                 + max(flx(3,i,j,k),0.0_dp)   + max(-flx(3,i,j,k-1),0.0_dp) )&
                 * dt / Jh(2*i,2*j,2*k)
            comax = max(comax, co)
         end do
        end do
       end do
    end do
    comax = global_max(comax)
    if (eno .and. Cnvt_mode == 3) sch = merge(1, 2, comax <= cfl_switch)    ! hybrid: explicit below the switch
    nsub  = max(1, ceiling(comax / vof_cfl_max))
    if (sch == 2 .and. eno) nsub = 1                                       ! implicit ENO: no sub-cycling
    dts   = dt / real(nsub, dp)
    cl    = 0.0_dp
    nit   = 0;  resmax = 0.0_dp;  resavg = 0.0_dp

    do isub = 1, nsub
       if (sch == 1) then
          !================ explicit CICSAM ==============================
          call face_values(dts, .false.)
          do ib = 1, nown
             d = own_list(ib)
             call sbv(d)
             do k = ks, ke
              do j = js, je
               do i = is, ie
                  if (cstat(i,j,k) /= 1) then
                     wrk(i,j,k) = avf(i,j,k);  cycle
                  end if
                  fe = flx(1,i,  j,k)*rhsu(i,  j,k)
                  fw = flx(1,i-1,j,k)*rhsu(i-1,j,k)
                  fn = flx(2,i,j,  k)*rhsv(i,j,  k)
                  fs = flx(2,i,j-1,k)*rhsv(i,j-1,k)
                  ft = flx(3,i,j,k  )*rhsw(i,j,k  )
                  fb = flx(3,i,j,k-1)*rhsw(i,j,k-1)
                  wrk(i,j,k) = avf(i,j,k) - dts*(fe - fw + fn - fs + ft - fb) &
                                            / Jh(2*i,2*j,2*k)
                  if (with_src) wrk(i,j,k) = (wrk(i,j,k) + dts*cavsc(i,j,k))/(1.0_dp + dts*cavsp(i,j,k))
               end do
              end do
             end do
             avf(is:ie, js:je, ks:ke) = wrk(is:ie, js:je, ks:ke)
          end do
       else
          !================ implicit CICSAM (iterative) ==================
          do ib = 1, nown
             d = own_list(ib)
             call sbv(d)
             rhsp(is:ie,js:je,ks:ke) = avf(is:ie,js:je,ks:ke)     ! alpha^n of the sub-step
             if (with_src) rhsp(is:ie,js:je,ks:ke) = rhsp(is:ie,js:je,ks:ke) + dts*cavsc(is:ie,js:je,ks:ke)
          end do
          do it = 1, vof_maxiter
             ! CICSAM weights from the current iterate: re-evaluated during
             ! the first n_beta sweeps (Picard on the non-linear weights),
             ! then frozen so the remaining sweeps converge the linear
             ! system to the residual limit
             if (it <= n_beta) call face_values(dts, .true.)
             rmax = 0.0_dp;  rsum = 0.0_dp;  ncell_tot = 0
             do ib = 1, nown
                d = own_list(ib)
                call sbv(d)
                do k = ks, ke
                 do j = js, je
                  do i = is, ie
                     if (cstat(i,j,k) /= 1) then
                        wrk(i,j,k) = avf(i,j,k);  cycle
                     end if
                     vol = Jh(2*i,2*j,2*k)
                     diag = 1.0_dp;  off = 0.0_dp
                     if (with_src) diag = diag + dts*cavsp(i,j,k)
                     call face_terms( flx(1,i,  j,k), rhsu(i,  j,k), avf(i+1,j,k))
                     call face_terms(-flx(1,i-1,j,k), rhsu(i-1,j,k), avf(i-1,j,k))
                     call face_terms( flx(2,i,j,  k), rhsv(i,j,  k), avf(i,j+1,k))
                     call face_terms(-flx(2,i,j-1,k), rhsv(i,j-1,k), avf(i,j-1,k))
                     call face_terms( flx(3,i,j,k  ), rhsw(i,j,k  ), avf(i,j,k+1))
                     call face_terms(-flx(3,i,j,k-1), rhsw(i,j,k-1), avf(i,j,k-1))
                     rr = rhsp(i,j,k) - diag*avf(i,j,k) - off       ! true residual
                     rmax = max(rmax, abs(rr));  rsum = rsum + abs(rr)
                     ncell_tot = ncell_tot + 1
                     anew = (rhsp(i,j,k) - off) / diag
                     wrk(i,j,k) = avf(i,j,k) + vof_sor*(anew - avf(i,j,k))
                  end do
                 end do
                end do
                avf(is:ie, js:je, ks:ke) = wrk(is:ie, js:je, ks:ke)
             end do
             call refresh_all(cur_iv)
             rmax = global_max(rmax)
             ncl  = global_sum(real(ncell_tot, dp))
             rsum = global_sum(rsum) / max(ncl, 1.0_dp)
             nit = nit + 1
             resmax = rmax;  resavg = rsum
             if (it >= vof_miniter .and. rmax < vof_tol) exit
          end do
       end if
       !---- boundedness: clip to [0,1] (volume accounted) -----------------
       do ib = 1, nown
          d = own_list(ib)
          call sbv(d)
          do k = ks, ke
           do j = js, je
            do i = is, ie
               anew = avf(i,j,k)
               if (anew < 0.0_dp) then
                  cl = cl + (0.0_dp - anew)*Jh(2*i,2*j,2*k);  anew = 0.0_dp
               else if (anew > 1.0_dp) then
                  cl = cl + (anew - 1.0_dp)*Jh(2*i,2*j,2*k);  anew = 1.0_dp
               end if
               avf(i,j,k) = anew
            end do
           end do
          end do
       end do
       call refresh_all(cur_iv)
    end do
    clipped = global_sum(cl)

  contains

    !-------- diagonal / off-diagonal contribution of one face ------------
    ! fout: outward volume flux of the cell through the face, bt: CICSAM
    ! weight (0 upwind, 1 downwind), an: neighbour value across the face
    subroutine face_terms(fout, bt, an)
      real(dp), intent(in) :: fout, bt, an
      real(dp) :: c
      c = dts * fout / vol
      if (fout >= 0.0_dp) then          ! cell is the donor
         diag = diag + c*(1.0_dp - bt)
         off  = off  + c*bt*an
      else                              ! cell is the acceptor
         diag = diag + c*bt
         off  = off  + c*(1.0_dp - bt)*an
      end if
    end subroutine face_terms

    !-------- face values (want_beta = .false.) or weights (.true.) ------
    subroutine face_values(dtsub, want_beta)
      real(dp), intent(in) :: dtsub
      logical,  intent(in) :: want_beta
      do ib = 1, nown
         d = own_list(ib)
         call sbv(d)
         do k = ks, ke
          do j = js, je
           do i = is-1, ie
              rhsu(i,j,k) = face_alpha(1, i, j, k, dtsub, want_beta)
           end do
          end do
         end do
         do k = ks, ke
          do j = js-1, je
           do i = is, ie
              rhsv(i,j,k) = face_alpha(2, i, j, k, dtsub, want_beta)
           end do
          end do
         end do
         do k = ks-1, ke
          do j = js, je
           do i = is, ie
              rhsw(i,j,k) = face_alpha(3, i, j, k, dtsub, want_beta)
           end do
          end do
         end do
      end do
    end subroutine face_values

    !-------- volume flux J U^q through face (i,j,k) of family q --------
    pure real(dp) function flx(q, i_, j_, k_)
      integer, intent(in) :: q, i_, j_, k_
      select case (q)
      case (1);  flx = Jh(2*i_+1,2*j_,2*k_)/fvw(2*i_+1,2*j_,2*k_)*(U1(i_,j_,k_) - ugq(1,i_,j_,k_))
      case (2);  flx = Jh(2*i_,2*j_+1,2*k_)/fvw(2*i_,2*j_+1,2*k_)*(U2(i_,j_,k_) - ugq(2,i_,j_,k_))
      case default; flx = Jh(2*i_,2*j_,2*k_+1)/fvw(2*i_,2*j_,2*k_+1)*(U3(i_,j_,k_) - ugq(3,i_,j_,k_))
      end select
    end function flx

    !-------- CICSAM face value (or weight beta) at face (i,j,k), family q
    pure real(dp) function face_alpha(q, i_, j_, k_, dtsub, want_beta)
      integer,  intent(in) :: q, i_, j_, k_
      real(dp), intent(in) :: dtsub
      logical,  intent(in) :: want_beta
      integer  :: e(3), cD(3), cA(3), cU(3), hf(3), s
      real(dp) :: F, aD, aA, aU, at, cf, cbc, uq, gam, aft, den, bt
      real(dp) :: gxi(3), gx(3), nf(3), gn, fn, cth

      F = flx(q, i_, j_, k_)
      e = 0;  e(q) = 1
      hf = [2*i_, 2*j_, 2*k_];  hf(q) = hf(q) + 1       ! face half index
      if (F >= 0.0_dp) then
         cD = [i_, j_, k_];  s = 1
      else
         cD = [i_, j_, k_] + e;  s = -1
      end if
      cA = cD + s*e
      cU = cD - s*e
      aD = avf(cD(1),cD(2),cD(3))
      aA = avf(cA(1),cA(2),cA(3))
      aU = avf(cU(1),cU(2),cU(3))

      face_alpha = merge(0.0_dp, aD, want_beta)          ! upwind default
      if (eno) then
         ! 2nd-order ENO (minmod-limited upwind) face value for the vapour
         ! fraction: between the donor value and the donor/acceptor mean
         aft = aD + 0.5_dp*minmod_(aD - aU, aA - aD)
         if (want_beta) then
            if (abs(aA - aD) > 1.0e-12_dp) face_alpha = min(max((aft - aD)/(aA - aD), 0.0_dp), 1.0_dp)
         else
            face_alpha = aft
         end if
         return
      end if
      den = aA - aU
      if (abs(den) < 1.0e-12_dp) return
      at = (aD - aU) / den
      if (at <= 0.0_dp .or. at >= 1.0_dp) return         ! outside the NVD band

      ! donor-cell Courant number
      cf = abs(F)*dtsub / Jh(2*cD(1),2*cD(2),2*cD(3))
      cf = max(cf, 1.0e-12_dp)
      ! HYPER-C (compressive bound) and ULTIMATE-QUICKEST
      cbc = min(1.0_dp, at/cf)
      uq  = min( (8.0_dp*cf*at + (1.0_dp-cf)*(6.0_dp*at + 3.0_dp))/8.0_dp, cbc )
      ! blending by the interface / face orientation at the donor cell
      gxi(1) = 0.5_dp*(avf(cD(1)+1,cD(2),cD(3)) - avf(cD(1)-1,cD(2),cD(3)))
      gxi(2) = 0.5_dp*(avf(cD(1),cD(2)+1,cD(3)) - avf(cD(1),cD(2)-1,cD(3)))
      gxi(3) = 0.5_dp*(avf(cD(1),cD(2),cD(3)+1) - avf(cD(1),cD(2),cD(3)-1))
      gx(1) = mq(1,1,2*cD(1),2*cD(2),2*cD(3))*gxi(1) + mq(1,2,2*cD(1),2*cD(2),2*cD(3))*gxi(2) &
            + mq(1,3,2*cD(1),2*cD(2),2*cD(3))*gxi(3)
      gx(2) = mq(2,1,2*cD(1),2*cD(2),2*cD(3))*gxi(1) + mq(2,2,2*cD(1),2*cD(2),2*cD(3))*gxi(2) &
            + mq(2,3,2*cD(1),2*cD(2),2*cD(3))*gxi(3)
      gx(3) = mq(3,1,2*cD(1),2*cD(2),2*cD(3))*gxi(1) + mq(3,2,2*cD(1),2*cD(2),2*cD(3))*gxi(2) &
            + mq(3,3,2*cD(1),2*cD(2),2*cD(3))*gxi(3)
      nf(1) = mq(1,q,hf(1),hf(2),hf(3))
      nf(2) = mq(2,q,hf(1),hf(2),hf(3))
      nf(3) = mq(3,q,hf(1),hf(2),hf(3))
      gn = sqrt(gx(1)**2 + gx(2)**2 + gx(3)**2)
      fn = sqrt(nf(1)**2 + nf(2)**2 + nf(3)**2)
      if (gn*fn > 1.0e-30_dp) then
         cth = abs(gx(1)*nf(1) + gx(2)*nf(2) + gx(3)*nf(3)) / (gn*fn)
         cth = min(cth, 1.0_dp)
      else
         cth = 1.0_dp
      end if
      gam = cth*cth                                  ! = min(k (cos2th+1)/2, 1), k = 1
      aft = gam*cbc + (1.0_dp - gam)*uq
      aft = min(max(aft, at), 1.0_dp)                ! donor .. acceptor band
      if (want_beta) then
         bt = (aft - at) / (1.0_dp - at)             ! weight of the acceptor
         face_alpha = min(max(bt, 0.0_dp), 1.0_dp)
      else
         face_alpha = aft*den + aU
         face_alpha = min(max(face_alpha, min(aD,aA)), max(aD,aA))
      end if
    end function face_alpha

    pure real(dp) function minmod_(a, b)
      real(dp), intent(in) :: a, b
      if (a*b <= 0.0_dp) then
         minmod_ = 0.0_dp
      else
         minmod_ = sign(min(abs(a), abs(b)), a)
      end if
    end function minmod_
  end subroutine advect_alpha

  !---------------------------------------------------------------------
  ! Total liquid volume (non-dimensional, all ranks).
  !---------------------------------------------------------------------
  !=====================================================================
  ! GEOMETRIC VoF (vof_scheme = 3): PLIC reconstruction + direction-
  ! split advection with the Weymouth-Yue divergence correction.
  !=====================================================================
  subroutine advect_alpha_plic(nsub, clipped)
    integer,  intent(out) :: nsub
    real(dp), intent(out) :: clipped
    integer  :: ib, d, i, j, k, isub, q, sw, order(3,3), anew_cnt
    real(dp) :: comax, dts, cl, anew, vol, cfl_lim
    real(dp) :: fp, fm, lp, lm, cc
    order = reshape([1,2,3, 2,3,1, 3,1,2], [3,3])
    !---- sub-cycles: face Courant number <= min(vof_cfl_max, 0.5) -----
    comax = 0.0_dp
    do ib = 1, nown
       d = own_list(ib);  call set_block(d)
       do k = ks, ke
        do j = js, je
         do i = is, ie
            if (cstat(i,j,k) /= 1) cycle
            vol = Jh(2*i,2*j,2*k)
            comax = max(comax, abs(flx(1,i,j,k))*dt/vol, abs(flx(1,i-1,j,k))*dt/vol,       &
                        abs(flx(2,i,j,k))*dt/vol, abs(flx(2,i,j-1,k))*dt/vol,               &
                        abs(flx(3,i,j,k))*dt/vol, abs(flx(3,i,j,k-1))*dt/vol)
         end do
        end do
       end do
    end do
    comax   = global_max(comax)
    cfl_lim = min(vof_cfl_max, 0.5_dp)
    nsub    = max(1, ceiling(comax / cfl_lim))
    dts     = dt / real(nsub, dp)
    cl      = 0.0_dp
    do isub = 1, nsub
       sw = mod(isub-1, 3) + 1
       do q = 1, 3
          call plic_sweep(order(q, sw), dts, cl)
          call refresh_all(8)
       end do
    end do
    clipped = global_sum(cl)
  end subroutine advect_alpha_plic

  !---------------------------------------------------------------------
  ! one direction-split sweep along axis q:
  !   alpha_i <- alpha_i - (L+ - L-)/V_i + c_i (F+ - F-) dts / V_i
  ! L+/L- = liquid volume crossing the +/- face (from the PLIC plane of
  ! the donor cell), F = face volume flux, c_i = 1 if alpha_i > 1/2 at
  ! the start of the sweep (Weymouth & Yue).  Faces are evaluated once
  ! per cell (its + face) plus the - face of the first cell of the range.
  !---------------------------------------------------------------------
  subroutine plic_sweep(q, dts, cl)
    integer,  intent(in)    :: q
    real(dp), intent(in)    :: dts
    real(dp), intent(inout) :: cl
    integer  :: ib, d, i, j, k, e(3), lo(3), hi(3), idx(3)
    real(dp) :: fq, vol, lq, anew, div
    real(dp), allocatable :: lf(:,:,:)      ! liquid volume through the + face of cell (i,j,k) along q (physical volume, signed with the flux)
    e = 0;  e(q) = 1
    do ib = 1, nown
       d = own_list(ib);  call set_block(d)
       lo = [is, js, ks];  hi = [ie, je, ke]
       allocate(lf(lo(1)-e(1):hi(1), lo(2)-e(2):hi(2), lo(3)-e(3):hi(3)))
       ! (a) liquid volume through every q-face of the range (+ the - face of the first cell)
       do k = lo(3)-e(3), hi(3)
        do j = lo(2)-e(2), hi(2)
         do i = lo(1)-e(1), hi(1)
            fq = flx(q, i, j, k)*dts                ! physical volume through the face, + along q
            if (fq >= 0.0_dp) then
               lf(i,j,k) = liquid_through(q, [i,j,k], fq, .true.)      ! donor = cell (i,j,k), slab at its + side
            else
               lf(i,j,k) = -liquid_through(q, [i,j,k]+e, -fq, .false.) ! donor = cell + e, slab at its - side
            end if
         end do
        end do
       end do
       ! (b) update of the active cells (wrk = new alpha; c_i from the current alpha)
       do k = ks, ke
        do j = js, je
         do i = is, ie
            if (cstat(i,j,k) /= 1) then
               wrk(i,j,k) = alpha(i,j,k);  cycle
            end if
            vol = Jh(2*i,2*j,2*k)
            idx = [i,j,k]
            div = (flx(q,i,j,k) - flx(q,i-e(1),j-e(2),k-e(3)))*dts        ! net volume out along q
            lq  = lf(i,j,k) - lf(i-e(1),j-e(2),k-e(3))                    ! net liquid out along q
            anew = alpha(i,j,k) - lq/vol + merge(1.0_dp, 0.0_dp, alpha(i,j,k) > 0.5_dp)*div/vol
            wrk(i,j,k) = anew
         end do
        end do
       end do
       do k = ks, ke
        do j = js, je
         do i = is, ie
            anew = wrk(i,j,k)
            if (anew < 0.0_dp) then
               cl = cl - anew*Jh(2*i,2*j,2*k);  anew = 0.0_dp
            else if (anew > 1.0_dp) then
               cl = cl + (anew - 1.0_dp)*Jh(2*i,2*j,2*k);  anew = 1.0_dp
            end if
            if (anew < 1.0e-8_dp) anew = 0.0_dp
            if (anew > 1.0_dp - 1.0e-8_dp) anew = 1.0_dp
            alpha(i,j,k) = anew
         end do
        end do
       end do
       deallocate(lf)
    end do
  end subroutine plic_sweep

  !---------------------------------------------------------------------
  ! liquid volume (physical) carried through a q-face of donor cell c by
  ! the volume fq (>= 0): the liquid inside the slab of the cell's unit
  ! cube adjacent to that face, of width fq / V_c.  plus = .true.: the
  ! + face (slab x_q in [1-w, 1]); .false.: the - face (slab [0, w]).
  ! Pure cells need no reconstruction.
  !---------------------------------------------------------------------
  real(dp) function liquid_through(q, c, fq, plus) result(lv)
    integer,  intent(in) :: q, c(3)
    real(dp), intent(in) :: fq
    logical,  intent(in) :: plus
    real(dp) :: a, vol, w, m(3), acst, x0
    a = alpha(c(1),c(2),c(3))
    if (fq <= 0.0_dp .or. a <= 0.0_dp) then
       lv = 0.0_dp;  return
    end if
    if (a >= 1.0_dp) then
       lv = fq;  return
    end if
    vol = Jh(2*c(1),2*c(2),2*c(3))
    w   = min(fq/vol, 1.0_dp)
    call plic_plane(c, m, acst)                        ! m >= 0, sum m = 1, liquid: m.x <= acst
    ! the slab in the (possibly mirrored) cube: the + face of the cell is
    ! at x_q = 1 unless the axis was mirrored (flag in m's sign handling
    ! below through plic_plane's mirror array)
    if (plus .neqv. mirrored(q)) then
       x0 = 1.0_dp - w
    else
       x0 = 0.0_dp
    end if
    lv = vol*slab_volume(m, acst, q, x0, w)
  end function liquid_through

  !---------------------------------------------------------------------
  ! PLIC plane of the mixed cell c in its computational unit cube:
  ! Youngs normal from the 3x3x3 stencil of alpha (index space), axes
  ! mirrored so that all components are >= 0 (mirrored(:) records it),
  ! normalized to sum 1; acst from the inverse of the volume relation.
  !---------------------------------------------------------------------
  subroutine plic_plane(c, m, acst)
    integer,  intent(in)  :: c(3)
    real(dp), intent(out) :: m(3), acst
    real(dp) :: st(-1:1,-1:1,-1:1)
    integer  :: d1, d2, d3
    do d3 = -1, 1
       do d2 = -1, 1
          do d1 = -1, 1
             st(d1,d2,d3) = alpha(c(1)+d1, c(2)+d2, c(3)+d3)
          end do
       end do
    end do
    call plic_normal27(st, m, mirrored)
    acst = plic_plane_const(m, alpha(c(1),c(2),c(3)))
  end subroutine plic_plane

  ! liquid fraction (of the unit cube) inside the slab x_q in [x0, x0+w]
  pure real(dp) function slab_volume(m, acst, q, x0, w) result(v)
    real(dp), intent(in) :: m(3), acst, x0, w
    integer,  intent(in) :: q
    real(dp) :: x0v(3), wv(3)
    x0v = 0.0_dp;  wv = 1.0_dp;  x0v(q) = x0;  wv(q) = w
    v = w*plic_box_frac(m, acst, x0v, wv)
  end function slab_volume

  ! face volume flux J_f U_f / fvw_f of face (q; i,j,k) (module-level copy)
  pure real(dp) function flx(q, i_, j_, k_)
    integer, intent(in) :: q, i_, j_, k_
    select case (q)
    case (1);  flx = Jh(2*i_+1,2*j_,2*k_)/fvw(2*i_+1,2*j_,2*k_)*(U1(i_,j_,k_) - ugq(1,i_,j_,k_))
    case (2);  flx = Jh(2*i_,2*j_+1,2*k_)/fvw(2*i_,2*j_+1,2*k_)*(U2(i_,j_,k_) - ugq(2,i_,j_,k_))
    case default; flx = Jh(2*i_,2*j_,2*k_+1)/fvw(2*i_,2*j_,2*k_+1)*(U3(i_,j_,k_) - ugq(3,i_,j_,k_))
    end select
  end function flx

  real(dp) function liquid_volume()
    integer :: ib, d, i, j, k
    real(dp) :: s
    s = 0.0_dp
    do ib = 1, nown
       d = own_list(ib)
       call set_block(d)
       do k = ks, ke
        do j = js, je
         do i = is, ie
            if (cstat(i,j,k) /= 1) cycle
            s = s + alpha(i,j,k)*Jh(2*i,2*j,2*k)*ovs_overlap_weight(dom(d)%b, i, j, k)
         end do
        end do
       end do
    end do
    liquid_volume = global_sum(s)
  end function liquid_volume

  !---------------------------------------------------------------------
  subroutine alpha_range(amin, amax)
    real(dp), intent(out) :: amin, amax
    integer :: ib, d
    real(dp) :: lo, hi
    lo = huge(1.0_dp);  hi = -huge(1.0_dp)
    do ib = 1, nown
       d = own_list(ib)
       call set_block(d)
       lo = min(lo, minval(alpha(is:ie,js:je,ks:ke)))
       hi = max(hi, maxval(alpha(is:ie,js:je,ks:ke)))
    end do
    amin = -global_max(-lo)
    amax =  global_max(hi)
  end subroutine alpha_range

end module vof
