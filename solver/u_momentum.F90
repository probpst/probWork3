!=======================================================================
! u_momentum.F90
!   Predictor of the projection method on the fully staggered
!   CURVILINEAR grid, Cartesian face velocities as unknowns (no
!   Christoffel symbols).  See convection.F90 for the operators and
!   the time-integration switches (Cnvt_mode, Dffs_mode).
!
!   Implicit pieces (masks cm/dm = 1 per direction), Jacobi iteration:
!     * diffusion, compact (q = l) parts: 7-point implicit diagonal /
!       neighbour coefficients mu_h G_qq / (J_f rho_f Re);
!     * diffusion, cross (q /= l) parts: lagged (deferred) inside the
!       iteration (Park Eq. 11 corner averages);
!     * convection: FIRST-ORDER UPWIND with the advecting fluxes of
!       time level n frozen (linearized) enters the matrix
!       (diagonal part max(a+,0) - min(a-,0), diagonally dominant);
!       for Cnvt_scheme = 2, 3, 4 (2nd-order upwind, QUICK, 2nd-order
!       ENO) the difference between the higher-order and the 1st-order
!       upwind fluxes is a deferred correction evaluated from the
!       latest velocity iterate and frozen during the sweeps (a source
!       term); with inner iterations (Iner_Itn > 1) the converged
!       solution is the fully implicit higher-order scheme, with
!       Iner_Itn = 1 the classic explicit deferred correction.
!
!   Convergence: the TRUE ALGEBRAIC RESIDUAL r = b - A u* of every face
!   equation, ||r||_inf / ||b||_inf < Rsd_Lmt_Vel (b = the momentum
!   right-hand side of all three components), after at least
!   Min_Itr_Vel and at most Max_Itr_Vel weighted-Jacobi sweeps
!   (SOR_fct_Vel).  Max and average residuals are reported per
!   component.
!
!   All three components iterate in ONE Jacobi loop; every sweep the
!   current iterate is transformed to the face components
!   (cart_to_contra), exchanged through the unified b2b/halo engine
!   with the boundary conditions applied in the transformed space
!   (refresh_uvw(5)), and the Cartesian ghost values are reconstructed
!   (contra_to_cart) - so arbitrary block orientations and curved
!   interfaces are handled by the same machinery.
!=======================================================================
subroutine solve_momentum(iter_out, resmax, resavg)
  use global_vars
  use mpi_vars
  use b2b_runtime, only: refresh_uvw
  use curvi, only: cart_to_contra, contra_to_cart
  use overset, only: ovs_update_faces
  use mom_ops
  implicit none
  integer,  intent(out) :: iter_out
  real(dp), intent(out) :: resmax(3), resavg(3)   ! per component u, v, w

  integer  :: i, j, k, it, ib, b, c, q, l, np, lo(3), hi(3), pt(3), H(3), e(3)
  real(dp) :: ap(3), am(3), cm(3), dm(3), nuf, apc(3), amc(3)
  real(dp) :: diag, off, qnew, dsum, cdiag, coff, dc, Jf, rr, nuJ
  real(dp) :: f2p, f1p, f2m, f1m, bnorm, rmx(3), rsm(3), ncnt(3)
  real(dp), pointer :: fs(:,:,:), r(:,:,:), mc(:,:,:,:)

  ! ---- normalization: ||b||_inf over the three right-hand sides ---------
  bnorm = 0.0_dp
  do ib = 1, nown
     b = own_list(ib)
     call set_block(b)
     do c = 1, 3
        select case (c)
        case (1);  r => rhsu
        case (2);  r => rhsv
        case default;  r => rhsw
        end select
        call face_range(c, lo, hi)
        do k = lo(3), hi(3)
         do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             if (face_active(c, i, j, k)) bnorm = max(bnorm, abs(r(i,j,k)))
          end do
         end do
        end do
     end do
  end do
  bnorm = global_max(bnorm)
  if (bnorm <= 0.0_dp) bnorm = 1.0_dp

  ! ---- coefficient assembly (once per solve: the implicit operator is
  !      linear in u* with the advecting fluxes and masks frozen) --------
  !  mc( 1)      : 1/dt + implicit compact diffusion diagonal
  !  mc( 2: 7)   : compact neighbour coefficients f(pt +/- e_q), q = 1..3
  !  mc( 8:19)   : cross-derivative coefficients of the pairs (q,l) =
  !                (1,2),(1,3),(2,1),(2,3),(3,1),(3,2): plus face (times
  !                xder(f,pt,pt+e_q,l)) and minus face (times
  !                xder(f,pt-e_q,pt,l), sign folded in)
  !  mc(20:22)   : cm(q) ap(q)/J_f,  mc(23:25) : cm(q) am(q)/J_f
  !  mc(26)      : Cnvt_scheme >= 2: deferred correction (high - 1st order
  !                upwind fluxes) of the implicit directions, evaluated
  !                ONCE from the latest iterate u^m (frozen source term
  !                during the sweeps; converges to the fully implicit
  !                2nd-order scheme through the inner iterations).  A
  !                correction updated inside the Jacobi loop is not a
  !                contraction for CFL > 1/2 (the Jacobi sweeps diverge).
  do ib = 1, nown
     b = own_list(ib)
     call set_block(b)
     do c = 1, 3
        mc(1:, lbound(dom(b)%mco,2):, lbound(dom(b)%mco,3):, lbound(dom(b)%mco,4):) => dom(b)%mco(:,:,:,:,c)
        select case (c)
        case (1);  fs => u
        case (2);  fs => v
        case default;  fs => w
        end select
        call face_range(c, lo, hi)
        do k = lo(3), hi(3)
         do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             if (.not. face_active(c, i, j, k)) cycle
             pt = [i, j, k]
             call face_setup(c, pt, H, ap, am, cm, dm, nuf)
             Jf = Jh(H(1),H(2),H(3));  nuJ = nuf / Jf
             dsum = 0.0_dp
             do q = 1, 3
                e = 0;  e(q) = 1
                mc(2*q,   i,j,k) = nuJ*dm(q)*Gm(q,q,H+e)
                mc(2*q+1, i,j,k) = nuJ*dm(q)*Gm(q,q,H-e)
                dsum = dsum + mc(2*q,i,j,k) + mc(2*q+1,i,j,k)
                mc(19+q, i,j,k) = cm(q)*ap(q)/Jf
                mc(22+q, i,j,k) = cm(q)*am(q)/Jf
             end do
             mc(1,i,j,k) = 1.0_dp/dt + dsum
             dc = 0.0_dp
             if (Cnvt_scheme >= 2) then
                do q = 1, 3
                   e = 0;  e(q) = 1
                   f2p = upw (ap(q), fv(fs,pt-e), fv(fs,pt), fv(fs,pt+e), fv(fs,pt+2*e))
                   f1p = upw1(ap(q), fv(fs,pt), fv(fs,pt+e))
                   f2m = upw (am(q), fv(fs,pt-2*e), fv(fs,pt-e), fv(fs,pt), fv(fs,pt+e))
                   f1m = upw1(am(q), fv(fs,pt-e), fv(fs,pt))
                   dc = dc + cm(q)*( ap(q)*(f2p - f1p) - am(q)*(f2m - f1m) )/Jf
                end do
             end if
             mc(26,i,j,k) = dc
             np = 0
             do q = 1, 3
                e = 0;  e(q) = 1
                do l = 1, 3
                   if (l == q) cycle
                   mc(8+2*np, i,j,k) =  nuJ*dm(q)*Gm(q,l,H+e)
                   mc(9+2*np, i,j,k) = -nuJ*dm(q)*Gm(q,l,H-e)
                   np = np + 1
                end do
             end do
          end do
         end do
        end do
     end do
  end do

  resmax = 0.0_dp;  resavg = 0.0_dp
  iter_out = 0
  do it = 1, mom_maxiter
     call cart_to_contra(.true.)               ! us,vs,ws -> Us1,Us2,Us3
     call refresh_uvw(5)                       ! exchange Us + BCs
     call contra_to_cart(.true., .true.)       ! ghost us,vs,ws back
     if (ovst_on) call ovs_update_faces(5)     ! overset fringe faces <- donors
     rmx = 0.0_dp;  rsm = 0.0_dp;  ncnt = 0.0_dp

     do ib = 1, nown
        b = own_list(ib)
        call set_block(b)
        do c = 1, 3
           select case (c)
           case (1);  fs => us;  r => rhsu
           case (2);  fs => vs;  r => rhsv
           case default;  fs => ws;  r => rhsw
           end select
           mc(1:, lbound(dom(b)%mco,2):, lbound(dom(b)%mco,3):, lbound(dom(b)%mco,4):) => dom(b)%mco(:,:,:,:,c)
           call face_range(c, lo, hi)
           !----- Jacobi sweep of family c ---------------------------
           do k = lo(3), hi(3)
            do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                if (.not. face_active(c, i, j, k)) cycle
                pt = [i, j, k]
                ! implicit diffusion, compact parts + lagged cross fluxes
                off = mc(2,i,j,k)*fs(i+1,j,k) + mc(3,i,j,k)*fs(i-1,j,k)          &
                    + mc(4,i,j,k)*fs(i,j+1,k) + mc(5,i,j,k)*fs(i,j-1,k)          &
                    + mc(6,i,j,k)*fs(i,j,k+1) + mc(7,i,j,k)*fs(i,j,k-1)
                np = 0
                do q = 1, 3
                   e = 0;  e(q) = 1
                   do l = 1, 3
                      if (l == q) cycle
                      off = off + mc(8+2*np,i,j,k)*xder(fs, pt, pt+e, l)          &
                                + mc(9+2*np,i,j,k)*xder(fs, pt-e, pt, l)
                      np = np + 1
                   end do
                end do
                ! implicit convection: linearized 1st-order upwind (+ the
                ! frozen 2nd-order deferred correction mc(26))
                cdiag = 0.0_dp;  coff = 0.0_dp
                do q = 1, 3
                   e = 0;  e(q) = 1
                   apc(q) = mc(19+q,i,j,k);  amc(q) = mc(22+q,i,j,k)
                   cdiag = cdiag + max(apc(q),0.0_dp) - min(amc(q),0.0_dp)
                   coff  = coff  - min(apc(q),0.0_dp)*fv(fs,pt+e)                &
                                 + max(amc(q),0.0_dp)*fv(fs,pt-e)
                end do
                diag = mc(1,i,j,k) + cdiag
                off  = off + coff - mc(26,i,j,k)
                ! true algebraic residual of this face equation
                rr = r(i,j,k) + off - diag*fs(i,j,k)
                if (abs(rr) > rmx(c)) then
                   rmx(c) = abs(rr);  mom_res_loc(:,c) = [dom(b)%b, i, j, k]
                end if
                rsm(c) = rsm(c) + abs(rr)
                ncnt(c) = ncnt(c) + 1.0_dp
                qnew = (r(i,j,k) + off) / diag
                wrk(i,j,k) = fs(i,j,k) + mom_sor*(qnew - fs(i,j,k))
             end do
            end do
           end do
           do k = lo(3), hi(3)
            do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                if (.not. face_active(c, i, j, k)) cycle
                fs(i,j,k) = wrk(i,j,k)
             end do
            end do
           end do
        end do
     end do

     mom_res_local = rmx / bnorm
     do c = 1, 3
        rmx(c)  = global_max(rmx(c)) / bnorm
        rsm(c)  = global_sum(rsm(c))
        ncnt(c) = global_sum(ncnt(c))
        rsm(c)  = rsm(c) / max(ncnt(c), 1.0_dp) / bnorm
     end do
     iter_out = it
     resmax = rmx;  resavg = rsm
     if (it >= mom_miniter .and. maxval(rmx) < mom_tol) exit
  end do

  call cart_to_contra(.true.)
  call refresh_uvw(5)
  call contra_to_cart(.true., .true.)
  if (ovst_on) call ovs_update_faces(5)
end subroutine solve_momentum
