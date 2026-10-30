!=======================================================================
! convection.F90
!   Shared momentum-operator pieces (module mom_ops) and the explicit
!   part of the momentum right-hand side.
!
!   Momentum equation solved for the Cartesian face velocities
!   (one-fluid formulation, divided by the local face density):
!
!     du/dt + CONV(u) = -(1/rho) grad p + (1/(rho Re)) div(mu grad u) + g
!
!   CONV in conservative flux form in the transformed space:
!     CONV = (1/J_f) sum_q [ (J U^q phi)_{+1/2 along q} - (J U^q phi)_{-1/2} ]
!   with the advecting fluxes J U^q interpolated to the half positions
!   around the face and the advected value upwinded along xi_q
!   (Cnvt_scheme = 1: 1st-order upwind, 2: 2nd-order upwind, 3: QUICK,
!   4: 2nd-order ENO; see upw).
!   DIFF = (1/(J_f rho_f Re)) sum_q D_q [ mu_h sum_l G_ql d(phi)/d(xi_l) ],
!   G_ql = J g^{ql}; the q = l parts are the compact 7-point terms, the
!   q /= l parts are corner-averaged cross derivatives (Park Eq. 11).
!
!   TIME INTEGRATION SWITCHES (per term and per direction q = i,j,k):
!     Cnvt_mode / Dffs_mode = 1 explicit, 2 implicit, 3 hybrid.
!   face_setup returns the masks cm(q), dm(q) = 1 (implicit) or 0
!   (explicit) for the unknown at hand.  In hybrid mode
!     cm(q) = 1  iff  CFL_q = dt |U^q| >= cfl_switch   (local, per direction)
!     dm(q) = 1  iff  VN_q = dt nu_eff g^{qq} >= vn_switch
!   Explicit parts are evaluated from u^n here (compute_rhs_momentum);
!   implicit parts enter the Jacobi matrix in u_momentum.F90.  Every
!   term is multiplied by its mask so the default path (explicit
!   convection, implicit diffusion, single phase) reproduces the
!   previous code bit for bit.
!
!   compute_rhs_momentum fills
!     rhs = u^n/dt - CONV_explicit(u^m) + DIFF_explicit(u^m) + g*
!   where u^n (un, vn, wn) is the time level n field and u^m the latest
!   iterate of the time step (= u^n in the first inner iteration; the
!   corrected velocity of the previous inner iteration afterwards, see
!   Iner_Itn in Main.F90).  The advecting fluxes ap/am of face_setup
!   come from U1..U3, which likewise hold the latest divergence-free
!   field.
!=======================================================================
module mom_ops
  use global_vars
  use mpi_vars, only: myrank
  implicit none
contains
  !--------------------------------------------------------------------
  ! upwind interpolation of the advected value (Cnvt_scheme)
  !--------------------------------------------------------------------
  ! Face value between the cells/faces fm (index i) and fp (i+1) for the
  ! advecting flux vel; fmm = i-1, fpp = i+2 (unit spacing in xi):
  !   Cnvt_scheme = 1  1st-order upwind      phi_f = phi_D
  !                2  2nd-order upwind      phi_f = 3/2 phi_D - 1/2 phi_U
  !                3  QUICK (3rd order)     phi_f = 3/8 phi_A + 3/4 phi_D - 1/8 phi_U
  !                4  2nd-order ENO         phi_f = phi_D + 1/2 m(phi_A - phi_D, phi_D - phi_U)
  !                   (D = donor, A = acceptor, U = upstream of D; m = minmod:
  !                    the smoother of the central and the upwind slope, 0
  !                    at an extremum -> essentially non-oscillatory,
  !                    conservative flux form of the 2nd-order ENO scheme)
  pure real(dp) function upw(vel, fmm, fm, fp, fpp)
    real(dp), intent(in) :: vel, fmm, fm, fp, fpp
    select case (Cnvt_scheme)
    case (2)                                   ! 2nd-order upwind
       if (vel >= 0.0_dp) then
          upw = 1.5_dp*fm - 0.5_dp*fmm
       else
          upw = 1.5_dp*fp - 0.5_dp*fpp
       end if
    case (3)                                   ! QUICK
       if (vel >= 0.0_dp) then
          upw = 0.375_dp*fp + 0.75_dp*fm - 0.125_dp*fmm
       else
          upw = 0.375_dp*fm + 0.75_dp*fp - 0.125_dp*fpp
       end if
    case (4)                                   ! 2nd-order ENO (minmod)
       if (vel >= 0.0_dp) then
          upw = fm + 0.5_dp*minmod(fp - fm, fm - fmm)
       else
          upw = fp - 0.5_dp*minmod(fp - fm, fpp - fp)
       end if
    case default                               ! 1st-order upwind
       if (vel >= 0.0_dp) then
          upw = fm
       else
          upw = fp
       end if
    end select
  end function upw

  ! m(x,y) = x if |x| <= |y|, y if |x| > |y|, 0 if x y < 0
  pure real(dp) function minmod(x, y)
    real(dp), intent(in) :: x, y
    if (x*y <= 0.0_dp) then
       minmod = 0.0_dp
    else if (abs(x) <= abs(y)) then
       minmod = x
    else
       minmod = y
    end if
  end function minmod

  pure real(dp) function upw1(vel, fm, fp)      ! 1st-order upwind
    real(dp), intent(in) :: vel, fm, fp
    if (vel >= 0.0_dp) then
       upw1 = fm
    else
       upw1 = fp
    end if
  end function upw1

  !--------------------------------------------------------------------
  ! Loop range of the face unknowns of family c handled by this
  ! sub-domain (owned interior faces + block boundary faces).
  !--------------------------------------------------------------------
  pure subroutine face_range(c, lo, hi)
    integer, intent(in)  :: c
    integer, intent(out) :: lo(3), hi(3)
    integer :: n(3)
    lo = [is, js, ks];  hi = [ie, je, ke];  n = [nx, ny, nz]
    if (lo(c) == 1)    lo(c) = 0
    if (hi(c) == n(c)) hi(c) = n(c)
  end subroutine face_range

  !--------------------------------------------------------------------
  ! .true. when the face unknown of family c at (i,j,k) is computed by
  ! this block (interior face, or an interface face this side owns).
  !--------------------------------------------------------------------
  pure logical function face_active(c, i, j, k)
    integer, intent(in) :: c, i, j, k
    face_active = face_owned(c, i, j, k)
    if (fstat(c, i, j, k) /= 1) face_active = .false.   ! overset: fringe / hole face
  end function face_active

  ! ownership only (interior face, or an interface face this side owns),
  ! regardless of the overset status
  pure logical function face_owned(c, i, j, k)
    integer, intent(in) :: c, i, j, k
    face_owned = .true.
    select case (c)
    case (1)
       if (i == 0) then
          face_owned = mu_lo(j,k)
       else if (i == nx) then
          face_owned = mu_hi(j,k)
       end if
    case (2)
       if (j == 0) then
          face_owned = mv_lo(i,k)
       else if (j == ny) then
          face_owned = mv_hi(i,k)
       end if
    case default
       if (k == 0) then
          face_owned = mw_lo(i,j)
       else if (k == nz) then
          face_owned = mw_hi(i,j)
       end if
    end select
  end function face_owned

  !--------------------------------------------------------------------
  ! Transformed component U^q at array position (a1,a2,a3)
  !--------------------------------------------------------------------
  pure real(dp) function Uq(q, a1, a2, a3)
    integer, intent(in) :: q, a1, a2, a3
    select case (q)
    case (1);  Uq = U1(a1,a2,a3)
    case (2);  Uq = U2(a1,a2,a3)
    case default;  Uq = U3(a1,a2,a3)
    end select
  end function Uq

  !--------------------------------------------------------------------
  ! Per-face setup for the unknown of family c at array position pt:
  !   H        : half index of the face
  !   ap(q)/am(q): advecting volume fluxes J U^q at the half positions
  !              H + e_q and H - e_q (built from U^n, 2-point averages)
  !   cm(q), dm(q): implicit masks (1.0 / 0.0) of convection/diffusion
  !   nuf      : 1/(Re rho_f) at the face
  !--------------------------------------------------------------------
  subroutine face_setup(c, pt, H, ap, am, cm, dm, nuf)
    integer,  intent(in)  :: c, pt(3)
    integer,  intent(out) :: H(3)
    real(dp), intent(out) :: ap(3), am(3), cm(3), dm(3), nuf
    integer  :: q, e(3), ec(3), Hp(3), Hm(3), Pm(3)
    real(dp) :: cfl, vn

    H = 2*pt;  H(c) = H(c) + 1
    ec = 0;  ec(c) = 1
    do q = 1, 3
       e = 0;  e(q) = 1
       Hp = H + e;  Hm = H - e;  Pm = pt - e
       ! advecting VOLUME fluxes at the half positions H +- e_q: the mean
       ! of the two q-face fluxes J_f U^q_f / fvw_f that bound the
       ! position - the same fluxes discrete continuity is built from,
       ! so a divergence-free field convects nothing.  (Until update 4
       ! the cell Jacobian times the mean UNWEIGHTED U^q was used; on
       ! stretched grids that is not the mean face flux, 6 % off in the
       ! first cell of a stretched inlet where the ghost spacing is a
       ! copy of the first cell: a uniform inflow then got a spurious
       ! convective source there and the projection answered with a
       ! pressure wiggle in the inlet cells.)
       ap(q) = 0.5_dp * ( fflux(q, pt) + fflux(q, pt + ec) )
       am(q) = 0.5_dp * ( fflux(q, Pm) + fflux(q, Pm + ec) )
    end do
    nuf = (1.0_dp/Re) / rho_h(H(1),H(2),H(3))

    select case (Cnvt_mode)
    case (1);  cm = 0.0_dp
    case (2);  cm = 1.0_dp
    case default
       do q = 1, 3
          e = 0;  e(q) = 1
          Hp = H + e;  Hm = H - e
          cfl = dt * max( abs(ap(q))/max(abs(Jh(Hp(1),Hp(2),Hp(3))), tiny(1.0_dp)), &
                          abs(am(q))/max(abs(Jh(Hm(1),Hm(2),Hm(3))), tiny(1.0_dp)) )
          cm(q) = merge(1.0_dp, 0.0_dp, cfl >= cfl_switch)
       end do
    end select
    select case (Dffs_mode)
    case (1);  dm = 0.0_dp
    case (2);  dm = 1.0_dp
    case default
       do q = 1, 3
          e = 0;  e(q) = 1
          Hp = H + e;  Hm = H - e
          vn = dt * nuf * max( mu_h(Hp(1),Hp(2),Hp(3))*gql(q,q,Hp(1),Hp(2),Hp(3)), &
                               mu_h(Hm(1),Hm(2),Hm(3))*gql(q,q,Hm(1),Hm(2),Hm(3)) )
          dm(q) = merge(1.0_dp, 0.0_dp, vn >= vn_switch)
       end do
    end select
  end subroutine face_setup

  ! net grid-velocity flux out of the control volume of the face unknown at pt (half index H),
  ! divided by J_f: the discrete divergence of the grid velocity (deforming grids, update 24)
  pure real(dp) function grid_div(pt, H)
    integer, intent(in) :: pt(3), H(3)
    integer :: q, e(3), ec(3), Pm(3)
    real(dp) :: s
    s = 0.0_dp
    ec = 0
    do q = 1, 3
       if (H(q) /= 2*pt(q)) ec(q) = 1              ! the family of the unknown (H = 2 pt + e_c)
    end do
    do q = 1, 3
       e = 0;  e(q) = 1;  Pm = pt - e
       s = s + 0.5_dp*(gflux(q, pt) + gflux(q, pt + ec)) - 0.5_dp*(gflux(q, Pm) + gflux(q, Pm + ec))
    end do
    grid_div = s/Jh(H(1),H(2),H(3))
  end function grid_div

  pure real(dp) function gflux(q, A)
    integer, intent(in) :: q, A(3)
    integer :: Hf(3)
    Hf = 2*A;  Hf(q) = Hf(q) + 1
    gflux = Jh(Hf(1),Hf(2),Hf(3)) / fvw(Hf(1),Hf(2),Hf(3)) * ugq(q, A(1), A(2), A(3))
  end function gflux

  ! volume flux through the q-face at index A: J_f U^q_f / fvw_f
  ! (RELATIVE to the grid on a moving block: U^q - Ug^q, ALE form)
  pure real(dp) function fflux(q, A)
    integer, intent(in) :: q, A(3)
    integer :: Hf(3)
    Hf = 2*A;  Hf(q) = Hf(q) + 1
    fflux = Jh(Hf(1),Hf(2),Hf(3)) / fvw(Hf(1),Hf(2),Hf(3)) * (Uq(q, A(1), A(2), A(3)) - ugq(q, A(1), A(2), A(3)))
  end function fflux

  !--------------------------------------------------------------------
  ! Viscosity-weighted diffusion coefficient mu_h * J g^{ql} at half
  ! position H.
  !--------------------------------------------------------------------
  pure real(dp) function Gm(q, l, H)
    integer, intent(in) :: q, l, H(3)
    Gm = Gcap(q, l, H(1),H(2),H(3)) * mu_h(H(1),H(2),H(3))
  end function Gm

  !--------------------------------------------------------------------
  ! Cross derivative d(f)/d(xi_l) at the half position between the
  ! array positions A and B (corner average, Park Eq. 11 style).
  !--------------------------------------------------------------------
  pure real(dp) function xder(f, A, B, l)
    real(dp), intent(in) :: f(is-ng:, js-ng:, ks-ng:)
    integer,  intent(in) :: A(3), B(3), l
    integer :: e(3)
    e = 0;  e(l) = 1
    xder = 0.25_dp*( f(A(1)+e(1),A(2)+e(2),A(3)+e(3))                     &
                   + f(B(1)+e(1),B(2)+e(2),B(3)+e(3))                     &
                   - f(A(1)-e(1),A(2)-e(2),A(3)-e(3))                     &
                   - f(B(1)-e(1),B(2)-e(2),B(3)-e(3)) )
  end function xder

  !--------------------------------------------------------------------
  ! Sum of the masked cross-derivative fluxes around the unknown at pt
  ! (face half index H):
  !   sum_q w(q) sum_{l/=q} [ Gm(q,l,H+e_q) xder(f,pt,pt+e_q,l)
  !                         - Gm(q,l,H-e_q) xder(f,pt-e_q,pt,l) ]
  ! accumulated in the fixed order (1,2),(1,3),(2,1),(2,3),(3,1),(3,2)
  ! so the all-implicit path matches the previous code exactly.
  !--------------------------------------------------------------------
  pure real(dp) function cross_sum(f, pt, H, w)
    real(dp), intent(in) :: f(is-ng:, js-ng:, ks-ng:)
    integer,  intent(in) :: pt(3), H(3)
    real(dp), intent(in) :: w(3)
    integer :: q, l, e(3)
    real(dp) :: s
    s = 0.0_dp
    do q = 1, 3
       e = 0;  e(q) = 1
       do l = 1, 3
          if (l == q) cycle
          s = s + w(q)*Gm(q,l,H+e)*xder(f, pt, pt+e, l)
          s = s - w(q)*Gm(q,l,H-e)*xder(f, pt-e, pt, l)
       end do
    end do
    cross_sum = s
  end function cross_sum

  !--------------------------------------------------------------------
  ! Masked compact (q = l) diffusion pieces around the unknown at pt:
  !   dsum  = sum_q w(q) [Gm(q,q,H+e_q) + Gm(q,q,H-e_q)]        (diagonal)
  !   osum  = sum_q w(q) [Gm(q,q,H+e_q) f(pt+e_q) + Gm(q,q,H-e_q) f(pt-e_q)]
  ! accumulated in the order 1+,1-,2+,2-,3+,3- (previous code order).
  !--------------------------------------------------------------------
  pure subroutine diag_diff(f, pt, H, w, dsum, osum)
    real(dp), intent(in)  :: f(is-ng:, js-ng:, ks-ng:)
    integer,  intent(in)  :: pt(3), H(3)
    real(dp), intent(in)  :: w(3)
    real(dp), intent(out) :: dsum, osum
    integer :: q, e(3)
    real(dp) :: gp, gmn
    dsum = 0.0_dp;  osum = 0.0_dp
    do q = 1, 3
       e = 0;  e(q) = 1
       gp  = w(q)*Gm(q,q,H+e)
       gmn = w(q)*Gm(q,q,H-e)
       dsum = dsum + gp + gmn
       osum = osum + gp*f(pt(1)+e(1),pt(2)+e(2),pt(3)+e(3))                  &
                   + gmn*f(pt(1)-e(1),pt(2)-e(2),pt(3)-e(3))
    end do
  end subroutine diag_diff

  !--------------------------------------------------------------------
  ! Masked explicit convection of f around the unknown at pt, divided
  ! by J_f:  sum_q w(q) [ ap(q) f_{+}(upw) - am(q) f_{-}(upw) ] / J_f
  ! in the order 1+,1-,2+,2-,3+,3- (previous code order).
  !--------------------------------------------------------------------
  pure real(dp) function conv_sum(f, pt, H, ap, am, w)
    real(dp), intent(in) :: f(is-ng:, js-ng:, ks-ng:)
    integer,  intent(in) :: pt(3), H(3)
    real(dp), intent(in) :: ap(3), am(3), w(3)
    integer :: q, e(3)
    real(dp) :: s, fp, fm
    s = 0.0_dp
    do q = 1, 3
       e = 0;  e(q) = 1
       fp = upw(ap(q), fv(f,pt-e), fv(f,pt), fv(f,pt+e), fv(f,pt+2*e))
       fm = upw(am(q), fv(f,pt-2*e), fv(f,pt-e), fv(f,pt), fv(f,pt+e))
       s = s + w(q)*ap(q)*fp - w(q)*am(q)*fm
    end do
    conv_sum = s / Jh(H(1),H(2),H(3))
  end function conv_sum

  pure real(dp) function fv(f, A)
    real(dp), intent(in) :: f(is-ng:, js-ng:, ks-ng:)
    integer,  intent(in) :: A(3)
    fv = f(A(1),A(2),A(3))
  end function fv

  !--------------------------------------------------------------------
  ! Make the face Jacobian of every shared interface face single-
  ! valued: non-owned interface faces take the OWNER's Jh at the face.
  ! Where two windows of a partial-face interface meet (a face that
  ! bends from one neighbour onto another) the receiver's ghost-node
  ! column at the junction is ambiguous and its ghost cells are
  ! distorted, so the two blocks would otherwise multiply the same
  ! exchanged face component U^q by different J_f: the volume flux
  ! through the face - and with it discrete continuity and VoF mass
  ! conservation - would differ between the two sides.  Called once
  ! after build_metrics; uses rhsu/rhsv/rhsw (ids 11..13) as carriers.
  !--------------------------------------------------------------------
  subroutine sync_face_jacobian(quiet)
    use b2b_runtime, only: exchange_faces
    logical, intent(in), optional :: quiet
    integer :: ib, d, q, i, j, k, lo(3), hi(3), pt(3), H(3)
    real(dp), pointer :: G(:,:,:)
    do ib = 1, nown
       d = own_list(ib)
       call set_block(d)
       do q = 1, 3
          G => field_ptr(d, 10 + q)
          G = 0.0_dp
          do k = lbound(G,3), ubound(G,3)
           do j = lbound(G,2), ubound(G,2)
            do i = lbound(G,1), ubound(G,1)
               H = 2*[i, j, k];  H(q) = H(q) + 1
               if (any(H <= hlo) .or. any(H >= hhi)) cycle
               G(i,j,k) = Jh(H(1),H(2),H(3))
            end do
           end do
          end do
       end do
    end do
    call exchange_faces(11)
    do ib = 1, nown
       d = own_list(ib)
       call set_block(d)
       do q = 1, 3
          G => field_ptr(d, 10 + q)
          call face_range(q, lo, hi)
          do k = lo(3), hi(3)
           do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               pt = [i, j, k]
               if (face_type(q, pt(q), pt) /= 0) cycle       ! physical boundary
               if (face_active(q, i, j, k)) cycle           ! owned / interior
               H = 2*pt;  H(q) = H(q) + 1
               Jh(H(1),H(2),H(3)) = abs(G(i,j,k))            ! owner's J (sign
            end do                                           !  factor removed)
           end do
          end do
       end do
    end do
    if (present(quiet)) then
       if (quiet) return
    end if
    if (myrank == 0) write(*,'(A)')                                       &
       '  curvilinear     : interface face Jacobians synchronized (owner values)'
  end subroutine sync_face_jacobian
end module mom_ops

!-----------------------------------------------------------------------
subroutine compute_rhs_momentum
  use global_vars
  use mpi_vars
  use mom_ops
  use surface_tension, only: st_force
  use energy, only: buoyancy_force
  use force_model, only: fm_on, force_source
  implicit none
  integer  :: i, j, k, ib, b, c, lo(3), hi(3), pt(3), H(3)
  real(dp) :: ap(3), am(3), cm(3), dm(3), ce(3), de(3), nuf
  real(dp) :: conv, dsum, osum, dexp
  real(dp), pointer :: f(:,:,:), r(:,:,:), fn(:,:,:)

  do ib = 1, nown
   b = own_list(ib)
   call set_block(b)
   do c = 1, 3
      select case (c)
      case (1);  f => u;  r => rhsu;  fn => un
      case (2);  f => v;  r => rhsv;  fn => vn
      case default;  f => w;  r => rhsw;  fn => wn
      end select
      call face_range(c, lo, hi)
      do k = lo(3), hi(3)
       do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           if (.not. face_active(c, i, j, k)) cycle
           pt = [i, j, k]
           call face_setup(c, pt, H, ap, am, cm, dm, nuf)
           ce = 1.0_dp - cm;  de = 1.0_dp - dm
           ! explicit convection (masked per direction)
           conv = conv_sum(f, pt, H, ap, am, ce)
           ! explicit diffusion (masked per direction): compact + cross
           call diag_diff(f, pt, H, de, dsum, osum)
           dexp = (osum - dsum*f(i,j,k)) + cross_sum(f, pt, H, de)
           r(i,j,k) = fn(i,j,k)/dt - conv + nuf*dexp/Jh(H(1),H(2),H(3))    &
                    + gnd(c)
           ! deforming grid (FSI mesh morphing, update 24): the conservative ALE convection of a
           ! cell whose volume changes carries the term -phi div(u_grid); it is removed here so
           ! that a uniform flow stays uniform on the morphing blocks (geometric conservation)
           if (cur_morph) r(i,j,k) = r(i,j,k) + f(i,j,k)*grid_div(pt, H)
           if (st_on) r(i,j,k) = r(i,j,k) + st_force(c, i, j, k, H)   ! surface tension (CSF)
           if (buoy_on) r(i,j,k) = r(i,j,k) + buoyancy_force(c, pt, H)   ! Boussinesq buoyancy (energy equation, update 16)
           if (fm_on) r(i,j,k) = r(i,j,k) + force_source(c, pt, H)     ! actuator disks / fixed point force (update 19)
        end do
       end do
      end do
   end do
  end do
end subroutine compute_rhs_momentum
