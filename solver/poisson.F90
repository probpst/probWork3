!=======================================================================
! poisson.F90
!   Pressure-Poisson equation, velocity correction and divergence
!   monitor on the fully staggered CURVILINEAR grid, variable density
!   (one-fluid two-phase formulation).
!
!   Divergence (exact, via the metric identity):
!     div(u) = (1/J_c) sum_q [ (J U^q)_face+ - (J U^q)_face- ]
!   Poisson operator (projection with the face density rho_f):
!     A p = (1/J_c) sum_q D_q [ (1/rho_f) sum_l G_ql d(p)/d(xi_l) ] = div(U*)/dt
!     G_ql = J * grad(xi_q).grad(xi_l) at the faces.
!   The l = q part is the 7-point core; the l /= q cross terms use
!   corner-averaged tangential p-derivatives (Park Eq. 11 style) and
!   vanish on rectilinear grids.
!   Corrector (transformed components AND the Cartesian face
!   velocities, both directly):
!     U^q_new = U^q* - dt (1/rho_f) sum_l g^{ql} d(p)/d(xi_l)   (x J for form 2)
!     u_new   = u*   - dt (1/rho_f) sum_l m(c,l)  d(p)/d(xi_l)
!   The SAME rho_f is used in the operator and in the corrector, so a
!   hydrostatic pressure field balances gravity exactly on the faces
!   (balanced-force property; no spurious currents from the density
!   jump at rest).
!
!   Solvers (poi_solver); both report the update norm max|p_new - p_old|
!   AND the true residual ||b - A p||_inf / ||b||_inf every time step:
!     1 = weighted point Jacobi (SOR_fct_Prs; cross terms lagged through
!         the previous iterate); criterion = update norm < Trl_Lmt_Prs
!     2 = flexible GMRES (restart poi_mkry = 20) on the full 19-point
!         operator (exchanged fluxes, true residual at every restart),
!         preconditioned by an inexact Jacobi-PCG solve of the SYMMETRIC
!         7-point part (face coefficients exchanged so that matrix is
!         exactly symmetric across interfaces; 1e-2 reduction, <= 200
!         CG iterations per application; CG preconditioner = Jacobi or
!         line symmetric Gauss-Seidel, poi_precond).  Robust for the
!         1000:1 density jump of the two-phase Poisson equation.  itp =
!         total inner CG iterations; criterion = residual < Rsd_Lmt_Prs
!         (checked at every restart, Min/Max_Itr_Prs honoured).  All-
!         Neumann problems: b is made compatible (J-weighted mean
!         removed) and the level is pinned.
!   The previous pressure is the initial guess of every solve.
!=======================================================================

!-----------------------------------------------------------------------
! Operator pieces shared by both solvers (module poisson_ops).
!
! The Poisson operator is assembled in FLUX FORM:
!   (1) face_grad(iv): the transformed, density-weighted pressure-
!       gradient flux of the cell field iv at every face this domain
!       computes (owned interior faces, owned interface faces, pressure
!       outlets; zero on fixed-flux physical boundaries)
!         G^q = G_qq Q_q + sum_{l/=q} G_ql <Q_l>,   Q_l = (dp/dxi_l)/rho_(l-face)
!       (balanced-force form, see Ql) stored in rhsu/rhsv/rhsw (ids 11..13);
!   (2) exchange_faces(11): non-owner interface faces and sub-domain
!       halo faces receive the value computed by the OWNING side;
!   (3) cell_div_G: (A x)_P = (1/J_c) sum_q [G^q(+) - G^q(-)].
! The corrector uses the SAME exchanged G^q, so the discrete continuity
! of the corrected velocity equals the Poisson residual cell by cell -
! also where two blocks discretize a shared face differently (window
! edges of partial-face interfaces, kinked grid lines).
!-----------------------------------------------------------------------
module poisson_ops
  use global_vars
  use mpi_vars
  use valves, only: valve_mob, valve_mob_h
  use mom_ops, only: xder, face_range, face_active, face_owned
  use b2b_runtime, only: exchange_faces
  use overset, only: dcel, e_ridx, ovs_val, ovs_dump_cell, ovs_ph
  implicit none
contains
  ! face types with a prescribed (star) flux
  pure logical function fixed_face(ft)
    integer, intent(in) :: ft
    fixed_face = (ft == 1 .or. ft == 3)
  end function fixed_face

  ! 1/rho_f-weighted coefficient at half position H
  pure real(dp) function Gr(q, l, H)
    integer, intent(in) :: q, l, H(3)
    Gr = Gcap(q, l, H(1),H(2),H(3)) / rho_h(H(1),H(2),H(3)) * valve_mob_h(q, H)   ! valve zone mobility (H2)
  end function Gr

  ! BALANCED directional pressure difference at the l-face of cell pt:
  !   Q_l = (x(pt+e_l) - x(pt)) / rho_(l-face)
  ! Every pressure difference is divided by the density of ITS OWN face
  ! before any averaging, so a hydrostatic pressure field (Q_l = g.t_l
  ! on every face) is reproduced exactly by the corner averages used in
  ! the cross terms and in the Cartesian correction - no spurious
  ! currents from the density jump (balanced-force formulation).
  pure real(dp) function Ql(x, pt, l)
    real(dp), intent(in) :: x(is-ng:, js-ng:, ks-ng:)
    integer,  intent(in) :: pt(3), l
    integer :: e(3), H(3)
    e = 0;  e(l) = 1
    H = 2*pt + e
    Ql = ( x(pt(1)+e(1),pt(2)+e(2),pt(3)+e(3)) - x(pt(1),pt(2),pt(3)) )     &
       / rho_h(H(1),H(2),H(3))
  end function Ql

  ! corner average of Q_l around the q-face between pt and pt+e_q
  pure real(dp) function Qavg(x, pt, q, l)
    real(dp), intent(in) :: x(is-ng:, js-ng:, ks-ng:)
    integer,  intent(in) :: pt(3), q, l
    integer :: eq(3), el(3)
    eq = 0;  eq(q) = 1;  el = 0;  el(l) = 1
    Qavg = 0.25_dp*( Ql(x, pt, l) + Ql(x, pt-el, l)                        &
                   + Ql(x, pt+eq, l) + Ql(x, pt+eq-el, l) )
  end function Qavg

  ! Jacobi diagonal of cell pt (compact parts, fixed-flux faces excluded)
  pure real(dp) function poi_diag(pt)
    integer, intent(in) :: pt(3)
    integer :: q, e(3), H(3)
    real(dp) :: s
    H = 2*pt;  s = 0.0_dp
    do q = 1, 3
       e = 0;  e(q) = 1
       if (.not. fixed_face(face_type(q, pt(q),   pt))) s = s + Gr(q,q,H+e)
       if (.not. fixed_face(face_type(q, pt(q)-1, pt))) s = s + Gr(q,q,H-e)
    end do
    poi_diag = s / Jh(H(1),H(2),H(3))
  end function poi_diag

  ! G^q of the cell field iv at all faces of every owned domain, then
  ! exchanged (ids 11..13 = rhsu, rhsv, rhsw)
  subroutine face_grad(iv)
    integer, intent(in) :: iv
    integer :: ib, d, q, i, j, k, lo(3), hi(3), pt(3), H(3), e(3), l, ft
    real(dp), pointer :: x(:,:,:), G(:,:,:)
    real(dp) :: s
    do ib = 1, nown
       d = own_list(ib)
       call set_block(d)
       x => field_ptr(d, iv)
       do q = 1, 3
          G => field_ptr(d, 10 + q)
          e = 0;  e(q) = 1
          call face_range(q, lo, hi)
          do k = lo(3), hi(3)
           do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               pt = [i, j, k]
               ft = face_type(q, pt(q), pt)
               if (fixed_face(ft)) then
                  G(i,j,k) = 0.0_dp                    ! fixed-flux boundary
                  cycle
               end if
               if (ft == 0 .and. .not. face_active(q, i, j, k)) then
                  G(i,j,k) = 0.0_dp                    ! received from the owner
                  cycle
               end if
               if (fstat(q, i, j, k) /= 1) then
                  G(i,j,k) = 0.0_dp                    ! overset fringe / hole face
                  cycle
               end if
               H = 2*pt;  H(q) = H(q) + 1
               s = Gcap(q,q,H(1),H(2),H(3))*Ql(x, pt, q)
               do l = 1, 3
                  if (l == q) cycle
                  s = s + Gcap(q,l,H(1),H(2),H(3))*Qavg(x, pt, q, l)
               end do
               G(i,j,k) = s*valve_mob(q, pt)                 ! valve zone: pressure flux x 1/(1 + K dt) (H2)
            end do
           end do
          end do
       end do
    end do
    call exchange_faces(11)
  end subroutine face_grad

  ! Fringe row of the pressure system (update 12):
  !   (1 + eps) x_P - I x - eps <x>_nb = c_P
  ! <x>_nb = average of the non-hole 6-neighbours of P in its own block
  ! (ghosts included), eps = ovs_prelax, c_P = ovs_ph (two-phase
  ! hydrostatic offset, iv = 4 only).  With eps = 0 the rows are the
  ! pure interpolation identities x_P = I x; chains of fringe cells
  ! interpolating each other across three or more grids (Box_6DoF:
  ! tube - cap - refinement - background) can close into cycles on
  ! which the identities are singular - the Krylov solve then returns
  ! arbitrary values there (|dp| ~ 1e19 with a converged residual).
  ! The neighbour average couples every fringe cell to its own block,
  ! which reaches solved cells, and removes the singularity; for a
  ! smooth field <x>_nb and I x agree to O(h^2), so the blend keeps
  ! the interpolation accuracy.
  pure real(dp) function fringe_row(x, i, j, k, ival)
    real(dp), intent(in) :: x(is-ng:, js-ng:, ks-ng:)
    integer,  intent(in) :: i, j, k
    real(dp), intent(in) :: ival           ! I x (interpolated value)
    real(dp) :: s, w
    s = 0.0_dp;  w = 0.0_dp
    if (ocs(i+1,j,k) >= 0) then;  s = s + x(i+1,j,k);  w = w + 1.0_dp;  end if
    if (ocs(i-1,j,k) >= 0) then;  s = s + x(i-1,j,k);  w = w + 1.0_dp;  end if
    if (ocs(i,j+1,k) >= 0) then;  s = s + x(i,j+1,k);  w = w + 1.0_dp;  end if
    if (ocs(i,j-1,k) >= 0) then;  s = s + x(i,j-1,k);  w = w + 1.0_dp;  end if
    if (ocs(i,j,k+1) >= 0) then;  s = s + x(i,j,k+1);  w = w + 1.0_dp;  end if
    if (ocs(i,j,k-1) >= 0) then;  s = s + x(i,j,k-1);  w = w + 1.0_dp;  end if
    if (w > 0.0_dp) then
       fringe_row = (1.0_dp + ovs_prelax)*x(i,j,k) - ival - ovs_prelax*s/w
    else
       fringe_row = x(i,j,k) - ival
    end if
  end function fringe_row

  ! ACTIVE cell without a single active face (all six faces fringe /
  ! hole / fixed-flux): its Poisson row is empty (no pressure flux
  ! through any face) and the system is singular there.  Such islands
  ! appear where the overset status leaves one active cell between
  ! fringe bands.  The row is replaced by the neighbour-average
  ! identity  x_P - <x>_nb = 0  (Laplace fill), see fringe_row.
  pure logical function isolated_active(i, j, k)
    integer, intent(in) :: i, j, k
    isolated_active = .false.
    if (ocs(i,j,k) /= 1) return
    if (fstat(1,i,j,k) == 1 .and. face_type(1,i,[i,j,k]) /= 1 .and. face_type(1,i,[i,j,k]) /= 3) return
    if (fstat(1,i-1,j,k) == 1 .and. face_type(1,i-1,[i,j,k]) /= 1 .and. face_type(1,i-1,[i,j,k]) /= 3) return
    if (fstat(2,i,j,k) == 1 .and. face_type(2,j,[i,j,k]) /= 1 .and. face_type(2,j,[i,j,k]) /= 3) return
    if (fstat(2,i,j-1,k) == 1 .and. face_type(2,j-1,[i,j,k]) /= 1 .and. face_type(2,j-1,[i,j,k]) /= 3) return
    if (fstat(3,i,j,k) == 1 .and. face_type(3,k,[i,j,k]) /= 1 .and. face_type(3,k,[i,j,k]) /= 3) return
    if (fstat(3,i,j,k-1) == 1 .and. face_type(3,k-1,[i,j,k]) /= 1 .and. face_type(3,k-1,[i,j,k]) /= 3) return
    isolated_active = .true.
  end function isolated_active

  ! (A x)_P = divergence of the exchanged G^q, owned cells of the
  ! current domain, into y
  ! Overset: the rows of the fringe cells are the interpolation
  ! identities  -diag_P (x_P - I x) = 0  (diag_P = the Jacobi diagonal,
  ! so the rows are scaled like the Poisson rows; I x = donor
  ! interpolation of the same vector, delivered by ovs_update_cells into
  ! ovs_val); hole rows are  x_P = 0.  iv = id of the cell field x.
  subroutine cell_div_G(y, iv)
    real(dp), intent(inout) :: y(is-ng:, js-ng:, ks-ng:)
    integer,  intent(in)    :: iv
    real(dp), pointer :: x(:,:,:)
    integer :: i, j, k, m, e
    do k = ks, ke
     do j = js, je
      do i = is, ie
         y(i,j,k) = ( rhsu(i,j,k) - rhsu(i-1,j,k) + rhsv(i,j,k) - rhsv(i,j-1,k) &
                    + rhsw(i,j,k) - rhsw(i,j,k-1) ) / Jh(2*i,2*j,2*k)
      end do
     end do
    end do
    x => field_ptr(cur_dom, iv)
    if (cav_on) then                        ! cavitation: pressure-implicit volume source (update 15)
       do k = ks, ke
        do j = js, je
         do i = is, ie
            if (cstat(i,j,k) == 1) y(i,j,k) = y(i,j,k) + (cavdp(i,j,k) - cavbc(i,j,k))/dt*x(i,j,k)
         end do
        end do
       end do
    end if
    if (.not. ovst_on) return
    do k = ks, ke
     do j = js, je
      do i = is, ie
         select case (ocs(i,j,k))
         case (0);  y(i,j,k) = -poi_diag([i,j,k]) * x(i,j,k)      ! fringe (orphan default)
         case (:-1); y(i,j,k) = -x(i,j,k)                        ! hole (same sign as the Poisson diagonal)
         case default
            if (isolated_active(i,j,k)) y(i,j,k) = -poi_diag([i,j,k]) * fringe_row(x, i, j, k, x(i,j,k))
         end select
      end do
     end do
    end do
    if (iv == 4) then                       ! the pressure: hydrostatic offset (two-phase)
       do m = 1, dcel(cur_dom)%n
          e = dcel(cur_dom)%e(m)
          i = e_ridx(1,e);  j = e_ridx(2,e);  k = e_ridx(3,e)
          y(i,j,k) = -poi_diag([i,j,k]) * (fringe_row(x, i, j, k, ovs_val(1,e)) - ovs_ph(e))
       end do
    else                                    ! Krylov vectors: homogeneous rows
       do m = 1, dcel(cur_dom)%n
          e = dcel(cur_dom)%e(m)
          i = e_ridx(1,e);  j = e_ridx(2,e);  k = e_ridx(3,e)
          y(i,j,k) = -poi_diag([i,j,k]) * fringe_row(x, i, j, k, ovs_val(1,e))
       end do
    end if
  end subroutine cell_div_G
end module poisson_ops

subroutine compute_rhs_poisson
  use global_vars
  use mpi_vars
  use poisson_ops, only: isolated_active
  implicit none
  integer  :: i, j, k, ib, b
  real(dp) :: div

  do ib = 1, nown
     b = own_list(ib)
     call set_block(b)
     do k = ks, ke
      do j = js, je
       do i = is, ie
          if (cstat(i,j,k) /= 1) then           ! overset fringe / hole cell
             rhsp(i,j,k) = 0.0_dp;  cycle
          end if
          if (ovst_on) then
             if (isolated_active(i,j,k)) then    ! active cell without an active face
                rhsp(i,j,k) = 0.0_dp;  cycle
             end if
          end if
          div = ( Jh(2*i+1,2*j,2*k)/fvw(2*i+1,2*j,2*k)*Us1(i,j,k) - Jh(2*i-1,2*j,2*k)/fvw(2*i-1,2*j,2*k)*Us1(i-1,j,k)  &
                + Jh(2*i,2*j+1,2*k)/fvw(2*i,2*j+1,2*k)*Us2(i,j,k) - Jh(2*i,2*j-1,2*k)/fvw(2*i,2*j-1,2*k)*Us2(i,j-1,k)  &
                + Jh(2*i,2*j,2*k+1)/fvw(2*i,2*j,2*k+1)*Us3(i,j,k) - Jh(2*i,2*j,2*k-1)/fvw(2*i,2*j,2*k-1)*Us3(i,j,k-1)  &
                ) / Jh(2*i,2*j,2*k)
          if (cav_on) div = div - cavs(i,j,k) + cavdp(i,j,k)*p(i,j,k)   & ! cavitation (update 15): div u = S(p) linearised,
                              - cavbc(i,j,k)*presn(i,j,k)                 !   S = S0 + dS/dp (p - p0); the dS/dp p part is in the
          rhsp(i,j,k) = div / dt                                        !   operator; update 17: - beta_c/dt (p - p^n) (compressible mixture)
       end do
      end do
     end do
  end do
end subroutine compute_rhs_poisson

!-----------------------------------------------------------------------
subroutine solve_poisson(iter_out, res_out, upd_out)
  use global_vars
  use mpi_vars
  use b2b_runtime, only: refresh_all, exchange_faces
  use mom_ops, only: face_range, face_active, face_owned
  use poisson_ops
  implicit none
  integer,  intent(out) :: iter_out
  real(dp), intent(out) :: res_out     ! residual ||b - A p||_inf / ||b||_inf
  real(dp), intent(out) :: upd_out     ! update norm max |p_new - p_old| (non-dim.)

  integer  :: i, j, k, it, ib, b, iout, ik, jk, nkry, ncg, ncell_tot
  real(dp) :: res, pshift, upd, dpk
  real(dp) :: bnorm, rnorm, sJ, sb, tmp, hsub
  real(dp) :: hg(poi_mkry+1, poi_mkry), gv(poi_mkry+1), yk(poi_mkry)
  real(dp) :: cs(poi_mkry), sn(poi_mkry)
  logical  :: pin
  integer  :: nan_cycles
  ! line-SGS preconditioner state (host-associated by line_solve)
  integer  :: lq, la, lb, nl, lo3(3), hi3(3)
  real(dp), allocatable :: laa(:), lbd(:), lcc(:), lrr(:), lxx(:)

  res_out = 0.0_dp;  upd_out = 0.0_dp
  iter_out = 0;  nan_cycles = 0
  pin = (p_ref_switch .and. .not. p_has_dirichlet)
  ! compressible mixture (update 17): the compressibility term fixes the level of an all-Neumann
  ! problem (the pressure rises with the compression) - no pinning, no compatibility correction
  if (cav_on .and. Cav_compressible == 1) pin = .false.
  ncell_tot = 0
  do ib = 1, nown
     b = own_list(ib)
     ncell_tot = ncell_tot + (dom(b)%r2(1)-dom(b)%r1(1)+1)*(dom(b)%r2(2)-dom(b)%r1(2)+1)*(dom(b)%r2(3)-dom(b)%r1(3)+1)
  end do
  ncell_tot = nint(global_sum(real(ncell_tot, dp)))
  call compute_rhs_poisson                     ! div(U*)/dt at the cells
  ! Jacobi diagonal (both solvers): wrk = 1/diag
  do ib = 1, nown
     b = own_list(ib)
     call set_block(b)
     do k = ks, ke
      do j = js, je
       do i = is, ie
          tmp = poi_diag([i, j, k])
          if (cav_on) tmp = tmp - (cavdp(i,j,k) - cavbc(i,j,k))/dt
          if (cstat(i,j,k) < 0) then
             wrk(i,j,k) = 1.0_dp                            ! overset hole row: x = 0
          else if (abs(tmp) > 1.0e-300_dp) then
             wrk(i,j,k) = 1.0_dp / tmp
          else
             wrk(i,j,k) = 0.0_dp                            ! isolated cell (all faces blocked): no update
             if (poi_verbose) write(*,'(A,4I6)') '  WARNING (poisson): zero diagonal at block/cell', b, i, j, k
          end if
       end do
      end do
     end do
  end do

  if (poi_solver == 1) then
     !=================== point Jacobi ===================================
     !  p_new = p + SOR_fct_Prs (rhs - A p)/diag ; the cross terms are
     !  inside A p (lagged through the previous iterate).  Both the
     !  update norm (criterion, Trl_Lmt_Prs) and the true residual
     !  (reported) come out of the same sweep.
     bnorm = vnorm_inf(0)
     if (bnorm <= 0.0_dp) bnorm = 1.0_dp
     do it = 1, poi_maxiter
        call refresh_all(4)                    ! p: ghosts + BCs
        call face_grad(4)
        res = 0.0_dp;  upd = 0.0_dp
        do ib = 1, nown
           b = own_list(ib)
           call set_block(b)
           call cell_div_G(kr, 4)              ! kr = A p
           do k = ks, ke
            do j = js, je
             do i = is, ie
                kr(i,j,k) = rhsp(i,j,k) - kr(i,j,k)          ! residual b - A p
                if (abs(kr(i,j,k)) > res) then
                   res = abs(kr(i,j,k));  poi_res_loc = [dom(b)%b, i, j, k]
                end if
                ! weighted Jacobi update (the Laplacian diagonal is -1/wrk)
                kr(i,j,k) = -poi_sor * kr(i,j,k) * wrk(i,j,k)
                if (abs(kr(i,j,k)) > upd) then
                   upd = abs(kr(i,j,k));  poi_upd_loc = [dom(b)%b, i, j, k]
                end if
             end do
            end do
           end do
           p(is:ie,js:je,ks:ke) = p(is:ie,js:je,ks:ke) + kr(is:ie,js:je,ks:ke)
        end do
        poi_res_local = res / bnorm;  poi_upd_local = upd
        res      = global_max(res) / bnorm
        upd      = global_max(upd)
        iter_out = it
        res_out  = res;  upd_out = upd
        if (pin) call pin_level
        if (it >= poi_miniter .and. upd < poi_tol) exit
     end do
     call refresh_all(4, 3)

  else
     !=================== FGMRES + 7-point PCG preconditioner ===========
     ! all-Neumann: make the right-hand side compatible (not with
     ! overset grids: the fringe rows are not conservative, the level
     ! is only pinned)
     if (pin .and. .not. ovst_on) then
        sJ = 0.0_dp;  sb = 0.0_dp
        do ib = 1, nown
           b = own_list(ib)
           call set_block(b)
           do k = ks, ke
            do j = js, je
             do i = is, ie
                sJ = sJ + Jh(2*i,2*j,2*k)
                sb = sb + Jh(2*i,2*j,2*k)*rhsp(i,j,k)
             end do
            end do
           end do
        end do
        sJ = global_sum(sJ);  sb = global_sum(sb)
        do ib = 1, nown
           b = own_list(ib)
           call set_block(b)
           rhsp(is:ie,js:je,ks:ke) = rhsp(is:ie,js:je,ks:ke) - sb/sJ
        end do
     end if
     call face_coefs                 ! symmetric face coefficients c_q, diag
     bnorm = vnorm_inf(0)
     if (bnorm <= 0.0_dp) bnorm = 1.0_dp
     iter_out = 0
     res_out = 1.0_dp
     restart: do iout = 1, 1000
        ! ---- true residual r = b - A p (full operator) -----------------
        call refresh_all(4, 3)
        call face_grad(4)
        do ib = 1, nown
           b = own_list(ib)
           call set_block(b)
           call cell_div_G(kr, 4)
           kr(is:ie,js:je,ks:ke) = rhsp(is:ie,js:je,ks:ke) - kr(is:ie,js:je,ks:ke)
        end do
        res_out = vnorm_inf(1) / bnorm
        if (ovst_on .and. myrank == 0 .and. poi_verbose) write(*,'(A,I5,A,I8,A,ES12.4)') &
           '    FGMRES restart ', iout, '  CG its ', iter_out, '  |b-Ap|/|b| ', res_out
        if ((iout > 1 .and. iter_out >= poi_miniter .and. res_out < poi_rtol)   &
            .or. iter_out >= poi_maxiter) exit restart
        if (iout == 1 .and. res_out < poi_rtol .and. poi_miniter <= 1) exit restart
        ! ---- Arnoldi with flexible preconditioning ---------------------
        rnorm = sqrt(vdot(1, 1))
        if (rnorm <= 0.0_dp) exit restart
        do ib = 1, nown
           b = own_list(ib)
           call set_block(b)
           kV(:,:,:,1) = kr(is:ie,js:je,ks:ke) / rnorm
        end do
        hg = 0.0_dp;  gv = 0.0_dp;  gv(1) = rnorm
        nkry = 0
        do jk = 1, poi_mkry
           ! z_j = M^-1 v_j  (inexact PCG on the 7-point part) -> kZ(:,j), kv1
           do ib = 1, nown
              b = own_list(ib)
              call set_block(b)
              do k = ks, ke
               do j = js, je
                do i = is, ie
                   kr(i,j,k) = -Jh(2*i,2*j,2*k)*kV(i,j,k,jk)
                end do
               end do
              end do
           end do
           if (ovst_on) then
              call gmres_c(ncg)            ! coupled inner operator (overset)
           else
              call pcg7(ncg)
           end if
           iter_out = iter_out + ncg
           ! w = A z_j
           call refresh_all(9, 3)
           call face_grad(9)
           do ib = 1, nown
              b = own_list(ib)
              call set_block(b)
              kZ(:,:,:,jk) = kv1(is:ie,js:je,ks:ke)
              call cell_div_G(kr, 9)
           end do
           ! modified Gram-Schmidt
           do ik = 1, jk
              hg(ik,jk) = vdotV(ik)
              do ib = 1, nown
                 b = own_list(ib)
                 call set_block(b)
                 kr(is:ie,js:je,ks:ke) = kr(is:ie,js:je,ks:ke) - hg(ik,jk)*kV(:,:,:,ik)
              end do
           end do
           hg(jk+1,jk) = sqrt(vdot(1, 1))
           hsub = hg(jk+1,jk)
           nkry = jk
           if (hsub /= hsub) then                      ! NaN from the inner solve: drop this direction
              nkry = jk - 1;  nan_cycles = nan_cycles + 1;  exit
           end if
           if (hsub > 1.0e-300_dp) then
              do ib = 1, nown
                 b = own_list(ib)
                 call set_block(b)
                 kV(:,:,:,jk+1) = kr(is:ie,js:je,ks:ke) / hg(jk+1,jk)
              end do
           end if
           ! Givens rotations
           do ik = 1, jk-1
              tmp          =  cs(ik)*hg(ik,jk) + sn(ik)*hg(ik+1,jk)
              hg(ik+1,jk)  = -sn(ik)*hg(ik,jk) + cs(ik)*hg(ik+1,jk)
              hg(ik,jk)    =  tmp
           end do
           tmp = sqrt(hg(jk,jk)**2 + hg(jk+1,jk)**2)
           if (tmp <= 1.0e-300_dp .or. tmp /= tmp) then
              nkry = jk - 1;  exit                     ! degenerate / NaN direction: discard it
           end if
           cs(jk) = hg(jk,jk)/tmp;  sn(jk) = hg(jk+1,jk)/tmp
           hg(jk,jk) = tmp;  hg(jk+1,jk) = 0.0_dp
           gv(jk+1) = -sn(jk)*gv(jk)
           gv(jk)   =  cs(jk)*gv(jk)
           ! |gv(jk+1)| = 2-norm of the current FGMRES residual (>= its
           ! inf-norm): leave the Arnoldi loop as soon as the criterion
           ! can be met, the true residual is checked at the restart
           if (abs(gv(jk+1)) < 0.5_dp*poi_rtol*bnorm .or. hsub <= 1.0e-300_dp     &
               .or. iter_out >= poi_maxiter) exit
        end do
        ! ---- update p = p + Z y,  H y = g ------------------------------
        do ik = nkry, 1, -1
           yk(ik) = gv(ik)
           do jk = ik+1, nkry
              yk(ik) = yk(ik) - hg(ik,jk)*yk(jk)
           end do
           if (abs(hg(ik,ik)) > 1.0e-300_dp) then
              yk(ik) = yk(ik) / hg(ik,ik)
           else
              yk(ik) = 0.0_dp
           end if
        end do
        if (any(yk(1:nkry) /= yk(1:nkry))) then       ! NaN in the Krylov solution: skip this cycle's update
           yk(1:nkry) = 0.0_dp
           nan_cycles = nan_cycles + 1
        end if
        upd = 0.0_dp
        do ib = 1, nown
           b = own_list(ib)
           call set_block(b)
           do k = ks, ke
            do j = js, je
             do i = is, ie
                dpk = 0.0_dp
                do ik = 1, nkry
                   dpk = dpk + yk(ik)*kZ(i,j,k,ik)
                end do
                p(i,j,k) = p(i,j,k) + dpk
                if (abs(dpk) > upd) then
                   upd = abs(dpk);  poi_upd_loc = [dom(b)%b, i, j, k]
                end if
             end do
            end do
           end do
        end do
        poi_upd_local = upd
        upd_out = global_max(upd)               ! update norm of this cycle
        if (pin) call pin_level
        if (nan_cycles > 0 .and. myrank == 0 .and. iout == 1) write(*,'(A)') &
           '  WARNING (poisson): a NaN Krylov direction was discarded (inner solver breakdown)'
        if (ovst_on .and. myrank == 0 .and. poi_verbose) write(*,'(A,I3,A,ES12.4,A,ES12.4,A,ES12.4)') &
           '      Arnoldi steps ', nkry, '  est. |r|_2 ', abs(gv(nkry+1)), '  target ', 0.5_dp*poi_rtol*bnorm, &
           '  bnorm ', bnorm
     end do restart
     call refresh_all(4, 3)
     if (ovst_on .and. poi_verbose) call ovs_diag
  end if

contains

  ! overset debugging: residual split by row type, fringe mismatch,
  ! linearity check of the operator paths (p vs kv1)
  subroutine ovs_diag
    real(dp) :: ra, rf, rh, mm, lin, sc
    integer :: la(4), lf(4), m, e
    ra = 0.0_dp;  rf = 0.0_dp;  rh = 0.0_dp;  mm = 0.0_dp;  la = 0;  lf = 0
    call face_grad(4)
    do ib = 1, nown
       b = own_list(ib)
       call set_block(b)
       call cell_div_G(kr, 4)
       do k = ks, ke
        do j = js, je
         do i = is, ie
            kr(i,j,k) = rhsp(i,j,k) - kr(i,j,k)
            select case (ocs(i,j,k))
            case (1)
               if (abs(kr(i,j,k)) > ra) then
                  ra = abs(kr(i,j,k));  la = [cur_blk, i, j, k]
               end if
            case (0)
               if (abs(kr(i,j,k)) > rf) then
                  rf = abs(kr(i,j,k));  lf = [cur_blk, i, j, k]
               end if
            case default;  rh = max(rh, abs(kr(i,j,k)))
            end select
         end do
        end do
       end do
       do m = 1, dcel(b)%n
          e = dcel(b)%e(m)
          mm = max(mm, abs(fringe_row(p, e_ridx(1,e),e_ridx(2,e),e_ridx(3,e), ovs_val(1,e)) - ovs_ph(e)))
       end do
    end do
    ! linearity: A(kv1 := p) vs A(p)
    do ib = 1, nown
       b = own_list(ib)
       call set_block(b)
       kv1 = p;  kt(is:ie,js:je,ks:ke) = kr(is:ie,js:je,ks:ke)
    end do
    call refresh_all(9, 3)
    call face_grad(9)
    lin = 0.0_dp;  sc = 0.0_dp
    do ib = 1, nown
       b = own_list(ib)
       call set_block(b)
       call cell_div_G(kr, 9)
       kr(is:ie,js:je,ks:ke) = rhsp(is:ie,js:je,ks:ke) - kr(is:ie,js:je,ks:ke)
       lin = max(lin, maxval(abs(kr(is:ie,js:je,ks:ke) - kt(is:ie,js:je,ks:ke))))
       sc  = max(sc, maxval(abs(kt(is:ie,js:je,ks:ke))))
    end do
    block
      integer :: cnt(3,-1:1), bb
      cnt = 0
      do ib = 1, nown
         b = own_list(ib)
         call set_block(b)
         do k = ks, ke
          do j = js, je
           do i = is, ie
              if (abs(kt(i,j,k)) > 1.0e-3_dp*bnorm) cnt(cur_blk, ocs(i,j,k)) = cnt(cur_blk, ocs(i,j,k)) + 1
           end do
          end do
         end do
      end do
      do bb = 1, 3
         write(*,'(A,I2,A,3I8)') '    diag: rows with |r|/|b| > 1e-3 in block ', bb, ' (hole/fringe/active):', cnt(bb,-1), cnt(bb,0), cnt(bb,1)
      end do
    end block
    write(*,'(A,ES11.3,A,4I5)') '    diag: max res active ', ra/bnorm, ' at blk,i,j,k', la
    write(*,'(A,ES11.3,A,4I5)') '          max res fringe ', rf/bnorm, ' at blk,i,j,k', lf
    write(*,'(A,ES11.3,A,ES11.3)') '          max res hole   ', rh/bnorm, '   fringe mismatch max|p_f - Ip| ', mm
    write(*,'(A,ES11.3,A,ES11.3)') '          |A(p) via kv1 - A(p) via p|max ', lin, '  (res max ', sc, ')'
    if (lf(1) > 0) call ovs_dump_cell(lf(1), lf(2), lf(3), lf(4))
    if (la(1) > 0) call ovs_dump_cell(la(1), la(2), la(3), la(4))
  end subroutine ovs_diag

  ! Symmetric face coefficients c_q = Gr(q,q,H) of the 7-point part,
  ! exchanged so both sides of an interface use the OWNER's value
  ! (carriers rhsu/rhsv/rhsw -> stored in kr0/kpv/kvv), zero on fixed-
  ! flux boundaries; and the Jacobi diagonal of Mpos in wrk.
  subroutine face_coefs
    integer :: q, lo(3), hi(3), pt(3), H(3), ft
    real(dp), pointer :: G(:,:,:), C(:,:,:)
    real(dp) :: sq(3)
    do ib = 1, nown
       b = own_list(ib)
       call set_block(b)
       do q = 1, 3
          G => field_ptr(b, 10 + q)
          G = 0.0_dp
          call face_range(q, lo, hi)
          do k = lo(3), hi(3)
           do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               pt = [i, j, k]
               ft = face_type(q, pt(q), pt)
               if (fixed_face(ft)) cycle
               if (ft == 0 .and. .not. face_owned(q, i, j, k)) cycle
               ! (overset: fringe / hole faces included - the overlapping
               !  Schwarz preconditioner needs them, see mask_coefs)
               H = 2*pt;  H(q) = H(q) + 1
               G(i,j,k) = Gr(q, q, H)
            end do
           end do
          end do
       end do
    end do
    call exchange_faces(11)
    do ib = 1, nown
       b = own_list(ib)
       call set_block(b)
       do q = 1, 3
          G => field_ptr(b, 10 + q)
          C => cptr(q)
          C = abs(G)                                   ! sign factor removed
          if (ovst_on) then
             C => coptr(q);  C = abs(G)                ! unmasked copy (coupled operator)
          end if
       end do
       do k = ks, ke
        do j = js, je
         do i = is, ie
            pt = [i, j, k]
            wrk(i,j,k) = 0.0_dp
            do q = 1, 3
               wrk(i,j,k) = wrk(i,j,k) + cw(q, pt, 1) + cw(q, pt, -1)
            end do
            if (cav_on) wrk(i,j,k) = wrk(i,j,k) + dcav(i, j, k)
            if (ovst_on) then
               ! fringe / hole rows of Mpos = -J A: J diag_P (operator row
               ! -diag_P (x_P - I x)) resp. J (identity row -x_P)
               if (ocs(i,j,k) == 0) wrk(i,j,k) = Jh(2*i,2*j,2*k)*poi_diag(pt)
               if (ocs(i,j,k) < 0)  wrk(i,j,k) = Jh(2*i,2*j,2*k)
            end if
            wrk(i,j,k) = 1.0_dp / wrk(i,j,k)
         end do
        end do
       end do

       ! line direction of the SGS preconditioner: strongest mean coupling
       if (poi_precond == 2) then
          sq = 0.0_dp
          do q = 1, 3
             C => cptr(q)
             sq(q) = sum(abs(C(is:ie,js:je,ks:ke)))
          end do
          dom(b)%qline = maxloc(sq, 1)
       end if
    end do
    ! overset: overlapping-Schwarz preconditioner, see mask_coefs
    if (ovst_on) call mask_coefs
  end subroutine face_coefs

  ! overset: the preconditioner couples ACTIVE cells only (block
  ! Jacobi across the grids with a Dirichlet closure at the fringes):
  ! the coefficients of every face without two active neighbour cells
  ! are dropped and their diagonal contribution is kept through ovd
  ! (apply_M7) / 1/wrk (line SGS); fringe and hole rows are diagonal.
  ! (Variants tried on the NACA66mod overset case: overlapping Schwarz
  ! with the dominated cells as overlap, restricted or not - both
  ! converged more slowly than this block-Jacobi form.)
  subroutine mask_coefs
    integer :: q, pt(3), e(3), nb(3)
    real(dp), pointer :: C(:,:,:)
    real(dp) :: dsum
    do ib = 1, nown
       b = own_list(ib)
       call set_block(b)
       ovd(:,:,:) = 0.0_dp
       do k = ks, ke
        do j = js, je
         do i = is, ie
            pt = [i, j, k]
            if (ocs(i,j,k) /= 1) then
               ovd(i,j,k) = 1.0_dp / wrk(i,j,k)
               cycle
            end if
            dsum = 0.0_dp
            do q = 1, 3
               e = 0;  e(q) = 1
               nb = pt + e
               if (ocs(nb(1),nb(2),nb(3)) /= 1) dsum = dsum + cw(q, pt, 1)
               nb = pt - e
               if (ocs(nb(1),nb(2),nb(3)) /= 1) dsum = dsum + cw(q, pt, -1)
            end do
            ovd(i,j,k) = dsum
         end do
        end do
       end do
       do q = 1, 3
          C => cptr(q)
          e = 0;  e(q) = 1
          do k = lbound(C,3), ubound(C,3)
           do j = lbound(C,2), ubound(C,2)
            do i = lbound(C,1), ubound(C,1)
               if (i+e(1) > ubound(ocs,1) .or. j+e(2) > ubound(ocs,2) .or. k+e(3) > ubound(ocs,3)) cycle
               if (ocs(i,j,k) /= 1 .or. ocs(i+e(1),j+e(2),k+e(3)) /= 1) C(i,j,k) = 0.0_dp
            end do
           end do
          end do
       end do
    end do
  end subroutine mask_coefs

  ! diagonal weight of face (+1: pt(q), -1: pt(q)-1) of cell pt:
  ! c (interior/interface), 2c at a pressure outlet (ghost = -x)
  real(dp) function cw(q, pt, side)
    integer, intent(in) :: q, pt(3), side
    integer :: fq
    real(dp), pointer :: C(:,:,:)
    real(dp) :: c0
    C => cptr(q)
    fq = merge(pt(q), pt(q)-1, side > 0)
    select case (q)
    case (1);  c0 = C(fq, pt(2), pt(3))
    case (2);  c0 = C(pt(1), fq, pt(3))
    case default;  c0 = C(pt(1), pt(2), fq)
    end select
    cw = c0
    if (face_type(q, fq, pt) == 2) cw = 2.0_dp*c0
  end function cw

  ! face coefficients of the preconditioner (masked when overset)
  function cptr(q) result(g)
    integer, intent(in) :: q
    real(dp), pointer :: g(:,:,:)
    g(lbound(kc,1):, lbound(kc,2):, lbound(kc,3):) => kc(:,:,:,q)
  end function cptr
  ! unmasked coefficients (overset: the coupled inner operator)
  function coptr(q) result(g)
    integer, intent(in) :: q
    real(dp), pointer :: g(:,:,:)
    g(lbound(kc,1):, lbound(kc,2):, lbound(kc,3):) => kc(:,:,:,3+q)
  end function coptr

  ! kt = Mpos x  for the search direction x = kv2 (ghosts refreshed),
  ! Mpos = -sum_q [ c_q(+)(x(+) - x) - c_q(-)(x - x(-)) ]  (owned cells)
  subroutine apply_M7
    real(dp), pointer :: C1(:,:,:), C2(:,:,:), C3(:,:,:)
    do ib = 1, nown
       b = own_list(ib)
       call set_block(b)
       C1 => cptr(1);  C2 => cptr(2);  C3 => cptr(3)
       do k = ks, ke
        do j = js, je
         do i = is, ie
            kt(i,j,k) = -( C1(i,j,k)*(kv2(i+1,j,k) - kv2(i,j,k))            &
                         - C1(i-1,j,k)*(kv2(i,j,k) - kv2(i-1,j,k))          &
                         + C2(i,j,k)*(kv2(i,j+1,k) - kv2(i,j,k))            &
                         - C2(i,j-1,k)*(kv2(i,j,k) - kv2(i,j-1,k))          &
                         + C3(i,j,k)*(kv2(i,j,k+1) - kv2(i,j,k))            &
                         - C3(i,j,k-1)*(kv2(i,j,k) - kv2(i,j,k-1)) )
            if (ovst_on) kt(i,j,k) = kt(i,j,k) + ovd(i,j,k)*kv2(i,j,k)
            if (cav_on) kt(i,j,k) = kt(i,j,k) + dcav(i, j, k)*kv2(i,j,k)
         end do
        end do
       end do
    end do
  end subroutine apply_M7

  ! cavitation (update 15): diagonal of Mpos = -J A from the pressure-
  ! implicit volume source, -J dS/dp / dt >= 0 (active cells only)
  real(dp) function dcav(i, j, k)
    integer, intent(in) :: i, j, k
    dcav = 0.0_dp
    if (cstat(i,j,k) /= 1) return
    dcav = -Jh(2*i,2*j,2*k)*(cavdp(i,j,k) - cavbc(i,j,k))/dt
  end function dcav

  ! overset: kt = Mfull x for x = kv2 (ghosts refreshed, fringe
  ! interpolation in ovs_val), Mfull = -J A restricted to the 7-point
  ! part: active rows with the UNMASKED coefficients (fringe neighbours
  ! included), fringe rows J diag_P (x_P - I x), hole rows J x_P.
  subroutine apply_Mfull
    real(dp), pointer :: C1(:,:,:), C2(:,:,:), C3(:,:,:)
    integer :: m, e
    do ib = 1, nown
       b = own_list(ib)
       call set_block(b)
       C1 => coptr(1);  C2 => coptr(2);  C3 => coptr(3)
       do k = ks, ke
        do j = js, je
         do i = is, ie
            if (ocs(i,j,k) == 1) then
               kt(i,j,k) = -( C1(i,j,k)*(kv2(i+1,j,k) - kv2(i,j,k))            &
                            - C1(i-1,j,k)*(kv2(i,j,k) - kv2(i-1,j,k))          &
                            + C2(i,j,k)*(kv2(i,j+1,k) - kv2(i,j,k))            &
                            - C2(i,j-1,k)*(kv2(i,j,k) - kv2(i,j-1,k))          &
                            + C3(i,j,k)*(kv2(i,j,k+1) - kv2(i,j,k))            &
                            - C3(i,j,k-1)*(kv2(i,j,k) - kv2(i,j,k-1)) )
               if (cav_on) kt(i,j,k) = kt(i,j,k) + dcav(i, j, k)*kv2(i,j,k)
               if (isolated_active(i,j,k)) kt(i,j,k) = Jh(2*i,2*j,2*k)*poi_diag([i,j,k])*fringe_row(kv2, i, j, k, kv2(i,j,k))
            else if (ocs(i,j,k) == 0) then
               kt(i,j,k) = Jh(2*i,2*j,2*k)*poi_diag([i,j,k])*kv2(i,j,k)   ! fringe (orphan default)
            else
               kt(i,j,k) = Jh(2*i,2*j,2*k)*kv2(i,j,k)                     ! hole
            end if
         end do
        end do
       end do
       do m = 1, dcel(b)%n
          e = dcel(b)%e(m)
          i = e_ridx(1,e);  j = e_ridx(2,e);  k = e_ridx(3,e)
          kt(i,j,k) = Jh(2*i,2*j,2*k)*poi_diag([i,j,k])*fringe_row(kv2, i, j, k, ovs_val(1,e))
       end do
    end do
  end subroutine apply_Mfull

  ! overset: inexact BiCGSTAB solve of  Mfull d = kr  (rhs in kr, owned
  ! cells) -> d in kv1, right-preconditioned by the masked line-SGS /
  ! Jacobi (precond: kr -> ksv).  The interpolation coupling between the
  ! grids is part of the operator, so the inner solve resolves it (the
  ! block-Jacobi PCG alone converges like a minimal-overlap Schwarz
  ! iteration: far too slowly).  Work: krr (r/s), kr0, kpv, kvv, ksv, kt.
  subroutine bicgstab_c(ncg)
    integer, intent(out) :: ncg
    real(dp) :: r0n, rn, alpha, beta, omega, rho, rho_new, r0v, tt, ts
    integer :: itc
    do ib = 1, nown
       b = own_list(ib)
       call set_block(b)
       kv1(is:ie,js:je,ks:ke) = 0.0_dp
       krr(is:ie,js:je,ks:ke) = kr(is:ie,js:je,ks:ke)
       kr0(is:ie,js:je,ks:ke) = kr(is:ie,js:je,ks:ke)
       kpv(is:ie,js:je,ks:ke) = kr(is:ie,js:je,ks:ke)
    end do
    r0n = sqrt(vdot(8, 8))
    ncg = 0
    if (r0n <= 0.0_dp) return
    rho = r0n*r0n
    do itc = 1, ovst_inner_max
       ! y = K^-1 p ;  v = Mfull y
       do ib = 1, nown
          b = own_list(ib);  call set_block(b)
          kr(is:ie,js:je,ks:ke) = kpv(is:ie,js:je,ks:ke)
       end do
       call precond
       do ib = 1, nown
          b = own_list(ib);  call set_block(b)
          kv2(is:ie,js:je,ks:ke) = ksv(is:ie,js:je,ks:ke)
       end do
       call refresh_all(10)
       call apply_Mfull
       do ib = 1, nown
          b = own_list(ib);  call set_block(b)
          kvv(is:ie,js:je,ks:ke) = kt(is:ie,js:je,ks:ke)
       end do
       r0v = vdot(2, 6)
       if (abs(r0v) < 1.0e-30_dp*max(rho, tiny(1.0_dp)) .or. r0v /= r0v) exit   ! breakdown / NaN: keep the last iterate
       alpha = rho / r0v
       do ib = 1, nown
          b = own_list(ib);  call set_block(b)
          krr(is:ie,js:je,ks:ke) = krr(is:ie,js:je,ks:ke) - alpha*kvv(is:ie,js:je,ks:ke)   ! s
          kv1(is:ie,js:je,ks:ke) = kv1(is:ie,js:je,ks:ke) + alpha*ksv(is:ie,js:je,ks:ke)
       end do
       ncg = itc
       rn = sqrt(vdot(8, 8))
       if (poi_verbose .and. myrank == 0 .and. iter_out == 0 .and.                    &
           (itc <= 5 .or. mod(itc,25) == 0)) write(*,'(A,I5,A,ES12.4)') '        bicgstab it ', itc, '  |s|/|r0| ', rn/r0n
       if (rn < ovst_inner_tol*r0n) exit
       ! z = K^-1 s ;  t = Mfull z
       do ib = 1, nown
          b = own_list(ib);  call set_block(b)
          kr(is:ie,js:je,ks:ke) = krr(is:ie,js:je,ks:ke)
       end do
       call precond
       do ib = 1, nown
          b = own_list(ib);  call set_block(b)
          kv2(is:ie,js:je,ks:ke) = ksv(is:ie,js:je,ks:ke)
       end do
       call refresh_all(10)
       call apply_Mfull
       tt = vdot(4, 4);  ts = vdot(4, 8)
       if (tt <= 0.0_dp .or. tt /= tt .or. ts /= ts) exit
       omega = ts / tt
       if (abs(omega) < 1.0e-30_dp) exit
       do ib = 1, nown
          b = own_list(ib);  call set_block(b)
          kv1(is:ie,js:je,ks:ke) = kv1(is:ie,js:je,ks:ke) + omega*ksv(is:ie,js:je,ks:ke)
          krr(is:ie,js:je,ks:ke) = krr(is:ie,js:je,ks:ke) - omega*kt(is:ie,js:je,ks:ke)
       end do
       rn = sqrt(vdot(8, 8))
       if (rn < ovst_inner_tol*r0n) exit
       rho_new = vdot(2, 8)
       if (abs(rho_new) < 1.0e-30_dp*r0n*r0n .or. rho_new /= rho_new) exit    ! near-breakdown: stop, FGMRES continues
       beta = (rho_new/rho)*(alpha/omega)
       if (abs(beta) > 1.0e30_dp) exit
       rho = rho_new
       do ib = 1, nown
          b = own_list(ib);  call set_block(b)
          kpv(is:ie,js:je,ks:ke) = krr(is:ie,js:je,ks:ke) + beta*(kpv(is:ie,js:je,ks:ke) - omega*kvv(is:ie,js:je,ks:ke))
       end do
    end do
  end subroutine bicgstab_c

  ! overset: inexact right-preconditioned GMRES(ovs_mkry) solve of
  ! Mfull d = kr (rhs in kr, owned cells) -> d in kv1, preconditioner
  ! K = masked line-SGS / Jacobi (precond: kr -> ksv).  Replaces the
  ! BiCGSTAB inner solve, which diverged on two-phase overset systems
  ! (density ratio 1000 + interpolation rows): the GMRES residual is
  ! monotone and there is no breakdown.  Restarted until ovst_inner_tol
  ! or ovst_inner_max operator applications.  Work: kW (basis), kr0
  ! (rhs), krr (residual), kv2/kt (operator in/out), ksv (precond out).
  subroutine gmres_c(ncg)
    integer, intent(out) :: ncg
    real(dp) :: hh(ovs_mkry+1, ovs_mkry), gg(ovs_mkry+1), yy(ovs_mkry)
    real(dp) :: cc(ovs_mkry), ss(ovs_mkry)
    real(dp) :: r0n, rn, hij, hsub_, t_
    integer :: itc, jj, ii, nk, irs
    do ib = 1, nown
       b = own_list(ib);  call set_block(b)
       kv1(is:ie,js:je,ks:ke) = 0.0_dp
       kr0(is:ie,js:je,ks:ke) = kr(is:ie,js:je,ks:ke)
       krr(is:ie,js:je,ks:ke) = kr(is:ie,js:je,ks:ke)
    end do
    r0n = sqrt(vdot(8, 8))
    ncg = 0
    if (r0n <= 0.0_dp) return
    restarts: do irs = 1, 1000
       rn = sqrt(vdot(8, 8))
       if (rn <= ovst_inner_tol*r0n .or. ncg >= ovst_inner_max .or. rn /= rn) exit restarts
       do ib = 1, nown
          b = own_list(ib);  call set_block(b)
          kW(:,:,:,1) = krr(is:ie,js:je,ks:ke) / rn
       end do
       hh = 0.0_dp;  gg = 0.0_dp;  gg(1) = rn;  nk = 0
       arnoldi: do jj = 1, ovs_mkry
          ! w = Mfull K^-1 v_j
          do ib = 1, nown
             b = own_list(ib);  call set_block(b)
             kr(is:ie,js:je,ks:ke) = kW(:,:,:,jj)
          end do
          call precond
          do ib = 1, nown
             b = own_list(ib);  call set_block(b)
             kv2(is:ie,js:je,ks:ke) = ksv(is:ie,js:je,ks:ke)
          end do
          call refresh_all(10)
          call apply_Mfull
          ncg = ncg + 1
          do ib = 1, nown
             b = own_list(ib);  call set_block(b)
             kr(is:ie,js:je,ks:ke) = kt(is:ie,js:je,ks:ke)
          end do
          ! modified Gram-Schmidt against W_1..W_j
          do ii = 1, jj
             hij = vdotW(ii)
             hh(ii,jj) = hij
             do ib = 1, nown
                b = own_list(ib);  call set_block(b)
                kr(is:ie,js:je,ks:ke) = kr(is:ie,js:je,ks:ke) - hij*kW(:,:,:,ii)
             end do
          end do
          hsub_ = sqrt(vdot(1, 1))
          hh(jj+1,jj) = hsub_
          if (hsub_ /= hsub_) exit arnoldi                    ! NaN: keep the previous directions
          nk = jj
          if (hsub_ > 1.0e-300_dp) then
             do ib = 1, nown
                b = own_list(ib);  call set_block(b)
                kW(:,:,:,jj+1) = kr(is:ie,js:je,ks:ke) / hsub_
             end do
          end if
          do ii = 1, jj-1
             t_          =  cc(ii)*hh(ii,jj) + ss(ii)*hh(ii+1,jj)
             hh(ii+1,jj) = -ss(ii)*hh(ii,jj) + cc(ii)*hh(ii+1,jj)
             hh(ii,jj)   =  t_
          end do
          t_ = sqrt(hh(jj,jj)**2 + hh(jj+1,jj)**2)
          if (t_ <= 1.0e-300_dp .or. t_ /= t_) then
             nk = jj - 1;  exit arnoldi
          end if
          cc(jj) = hh(jj,jj)/t_;  ss(jj) = hh(jj+1,jj)/t_
          hh(jj,jj) = t_;  hh(jj+1,jj) = 0.0_dp
          gg(jj+1) = -ss(jj)*gg(jj)
          gg(jj)   =  cc(jj)*gg(jj)
          if (poi_verbose .and. myrank == 0 .and. iter_out == 0 .and.                    &
              (ncg <= 5 .or. mod(ncg,25) == 0)) write(*,'(A,I5,A,ES12.4)') '        inner gmres it ', ncg, '  |r|/|r0| ', abs(gg(jj+1))/r0n
          if (abs(gg(jj+1)) <= ovst_inner_tol*r0n .or. hsub_ <= 1.0e-300_dp .or. ncg >= ovst_inner_max) exit arnoldi
       end do arnoldi
       if (nk == 0) exit restarts
       ! y = H^-1 g ;  u = W y ;  d = d + K^-1 u
       do ii = nk, 1, -1
          yy(ii) = gg(ii)
          do jj = ii+1, nk
             yy(ii) = yy(ii) - hh(ii,jj)*yy(jj)
          end do
          if (abs(hh(ii,ii)) > 1.0e-300_dp) then
             yy(ii) = yy(ii) / hh(ii,ii)
          else
             yy(ii) = 0.0_dp
          end if
       end do
       if (any(yy(1:nk) /= yy(1:nk))) exit restarts
       do ib = 1, nown
          b = own_list(ib);  call set_block(b)
          kr(is:ie,js:je,ks:ke) = 0.0_dp
          do ii = 1, nk
             kr(is:ie,js:je,ks:ke) = kr(is:ie,js:je,ks:ke) + yy(ii)*kW(:,:,:,ii)
          end do
       end do
       call precond
       do ib = 1, nown
          b = own_list(ib);  call set_block(b)
          kv1(is:ie,js:je,ks:ke) = kv1(is:ie,js:je,ks:ke) + ksv(is:ie,js:je,ks:ke)
       end do
       ! true residual r = b - Mfull d
       do ib = 1, nown
          b = own_list(ib);  call set_block(b)
          kv2(is:ie,js:je,ks:ke) = kv1(is:ie,js:je,ks:ke)
       end do
       call refresh_all(10)
       call apply_Mfull
       do ib = 1, nown
          b = own_list(ib);  call set_block(b)
          krr(is:ie,js:je,ks:ke) = kr0(is:ie,js:je,ks:ke) - kt(is:ie,js:je,ks:ke)
       end do
    end do restarts
    ! restore the rhs in kr for the caller (FGMRES continues with kv1)
    do ib = 1, nown
       b = own_list(ib);  call set_block(b)
       kr(is:ie,js:je,ks:ke) = kr0(is:ie,js:je,ks:ke)
    end do
  end subroutine gmres_c

  real(dp) function vdotW(iw)
    integer, intent(in) :: iw
    real(dp) :: s
    s = 0.0_dp
    do ib = 1, nown
       b = own_list(ib)
       call set_block(b)
       s = s + sum( kr(is:ie,js:je,ks:ke) * kW(:,:,:,iw) )
    end do
    vdotW = global_sum(s)
  end function vdotW

  ! inexact PCG solve of  Mpos d = kr  (rhs in kr, owned cells)
  ! -> d in kv1; uses kv2 (search dir), ksv (z), kt (q).  Stops at a
  ! relative 2-norm reduction of 1e-2 or 200 iterations: the result is
  ! only a preconditioner for the outer FGMRES.  The CG preconditioner
  ! is Jacobi (poi_precond = 1) or the line symmetric Gauss-Seidel
  ! sweep of line_sgs (poi_precond = 2).
  subroutine pcg7(ncg)
    integer, intent(out) :: ncg
    real(dp) :: r0n, rn, a_k, b_k, rz, rzn, pq
    integer :: itc
    do ib = 1, nown
       b = own_list(ib)
       call set_block(b)
       kv1(is:ie,js:je,ks:ke) = 0.0_dp
    end do
    call precond
    do ib = 1, nown
       b = own_list(ib)
       call set_block(b)
       kv2(is:ie,js:je,ks:ke) = ksv(is:ie,js:je,ks:ke)
    end do
    rz  = vdot(1, 5)
    r0n = sqrt(vdot(1, 1))
    ncg = 0
    if (r0n <= 0.0_dp) return
    do itc = 1, 200
       call refresh_all(10)          ! 7-point operator: one hop suffices
       call apply_M7
       pq = vdot(4, 6)
       if (abs(pq) < 1.0e-300_dp) exit
       a_k = rz / pq
       do ib = 1, nown
          b = own_list(ib)
          call set_block(b)
          kv1(is:ie,js:je,ks:ke) = kv1(is:ie,js:je,ks:ke) + a_k*kv2(is:ie,js:je,ks:ke)
          kr(is:ie,js:je,ks:ke)  = kr(is:ie,js:je,ks:ke)  - a_k*kt(is:ie,js:je,ks:ke)
       end do
       ncg = itc
       rn = sqrt(vdot(1, 1))
       if (rn < 1.0e-2_dp*r0n) exit
       call precond
       rzn = vdot(1, 5)
       b_k = rzn / rz;  rz = rzn
       do ib = 1, nown
          b = own_list(ib)
          call set_block(b)
          kv2(is:ie,js:je,ks:ke) = ksv(is:ie,js:je,ks:ke) + b_k*kv2(is:ie,js:je,ks:ke)
       end do
    end do
  end subroutine pcg7

  ! ksv = M^-1 kr  (owned cells of every domain)
  subroutine precond
    do ib = 1, nown
       b = own_list(ib)
       call set_block(b)
       if (poi_precond == 2) then
          call line_sgs(dom(b)%qline)
       else
          ksv(is:ie,js:je,ks:ke) = wrk(is:ie,js:je,ks:ke)*kr(is:ie,js:je,ks:ke)
       end if
    end do
  end subroutine precond

  ! Line symmetric Gauss-Seidel preconditioner of the 7-point operator
  ! Mpos = D_L + L + U on the current domain (block-Jacobi across sub-
  ! domains and interfaces: neighbours outside the owned range are
  ! dropped, which keeps M symmetric positive definite for the CG):
  !   M = (D_L + L) D_L^-1 (D_L + U),  D_L = tridiagonal part along ql
  ! Forward sweep: y = (D_L + L)^-1 r, keeping t = D_L y = r - L y;
  ! backward sweep: z = (D_L + U)^-1 t.  Lines are ordered by the two
  ! remaining indices; each line is a Thomas solve.  Work: kpy, kpt.
  subroutine line_sgs(ql)
    integer, intent(in) :: ql
    integer :: ib_, ic
    lq = ql
    lo3 = [is, js, ks];  hi3 = [ie, je, ke]
    nl = hi3(lq) - lo3(lq) + 1
    if (allocated(laa)) then
       if (size(laa) < nl) deallocate(laa, lbd, lcc, lrr, lxx)
    end if
    if (.not. allocated(laa)) allocate(laa(nl), lbd(nl), lcc(nl), lrr(nl), lxx(nl))
    ! the two line-ordering axes la < lb (both /= lq)
    select case (lq)
    case (1);  la = 2;  lb = 3
    case (2);  la = 1;  lb = 3
    case default;  la = 1;  lb = 2
    end select
    kpy(is:ie,js:je,ks:ke) = 0.0_dp
    ! ---- forward: lines in increasing (lb, la) order --------------------
    do ic = lo3(lb), hi3(lb)
     do ib_ = lo3(la), hi3(la)
        call line_solve(ic, ib_, .true.)
     end do
    end do
    ! ---- backward: lines in decreasing order ----------------------------
    ksv(is:ie,js:je,ks:ke) = 0.0_dp
    do ic = hi3(lb), lo3(lb), -1
     do ib_ = hi3(la), lo3(la), -1
        call line_solve(ic, ib_, .false.)
     end do
    end do
  end subroutine line_sgs

  ! one line of the SGS sweep (see line_sgs); Thomas algorithm
  subroutine line_solve(ic_, ib2, forward)
    integer, intent(in) :: ic_, ib2
    logical, intent(in) :: forward
    integer :: m_, q, nb(3), idx(3)
    real(dp) :: cf
    do m_ = 1, nl
       idx(lq) = lo3(lq) + m_ - 1;  idx(la) = ib2;  idx(lb) = ic_
       lbd(m_) = 1.0_dp / wrk(idx(1),idx(2),idx(3))       ! diagonal (incl. outlet 2c)
       ! tridiagonal couplings along lq (inside the owned range only)
       laa(m_) = 0.0_dp;  lcc(m_) = 0.0_dp
       if (m_ > 1)  laa(m_) = -cface(lq, idx, -1)
       if (m_ < nl) lcc(m_) = -cface(lq, idx, +1)
       if (forward) then
          lrr(m_) = kr(idx(1),idx(2),idx(3))
          ! L y: neighbours of lower line order (already solved)
          do q = 1, 3
             if (q == lq) cycle
             nb = idx;  nb(q) = nb(q) - 1
             if (nb(q) >= lo3(q)) lrr(m_) = lrr(m_) + cface(q, idx, -1)*kpy(nb(1),nb(2),nb(3))
          end do
          kpt(idx(1),idx(2),idx(3)) = lrr(m_)              ! t = r - L y (stored)
       else
          lrr(m_) = kpt(idx(1),idx(2),idx(3))
          ! U z: neighbours of higher line order (already solved)
          do q = 1, 3
             if (q == lq) cycle
             nb = idx;  nb(q) = nb(q) + 1
             if (nb(q) <= hi3(q)) lrr(m_) = lrr(m_) + cface(q, idx, +1)*ksv(nb(1),nb(2),nb(3))
          end do
       end if
    end do
    do m_ = 2, nl
       cf = laa(m_) / lbd(m_-1)
       lbd(m_) = lbd(m_) - cf*lcc(m_-1)
       lrr(m_) = lrr(m_) - cf*lrr(m_-1)
    end do
    lxx(nl) = lrr(nl) / lbd(nl)
    do m_ = nl-1, 1, -1
       lxx(m_) = (lrr(m_) - lcc(m_)*lxx(m_+1)) / lbd(m_)
    end do
    do m_ = 1, nl
       idx(lq) = lo3(lq) + m_ - 1;  idx(la) = ib2;  idx(lb) = ic_
       if (forward) then
          kpy(idx(1),idx(2),idx(3)) = lxx(m_)
       else
          ksv(idx(1),idx(2),idx(3)) = lxx(m_)
       end if
    end do
  end subroutine line_solve

  ! coefficient of the face on side (+1 / -1) of cell id along q
  real(dp) function cface(q, id, side)
    integer, intent(in) :: q, id(3), side
    real(dp), pointer :: C(:,:,:)
    integer :: fq
    C => cptr(q)
    fq = merge(id(q), id(q)-1, side > 0)
    select case (q)
    case (1);  cface = C(fq, id(2), id(3))
    case (2);  cface = C(id(1), fq, id(3))
    case default;  cface = C(id(1), id(2), fq)
    end select
  end function cface

  ! <kr, V_i> over owned cells (all ranks)
  real(dp) function vdotV(iv_)
    integer, intent(in) :: iv_
    real(dp) :: s
    s = 0.0_dp
    do ib = 1, nown
       b = own_list(ib)
       call set_block(b)
       s = s + sum( kr(is:ie,js:je,ks:ke) * kV(:,:,:,iv_) )
    end do
    vdotV = global_sum(s)
  end function vdotV

  ! global dot product over owned cells of two solver fields
  ! (codes: 1 = kr, 4 = kt, 5 = ksv, 6 = kv2, 7 = kv1)
  real(dp) function vdot(ca, cb)
    integer, intent(in) :: ca, cb
    real(dp), pointer :: fa(:,:,:), fb(:,:,:)
    real(dp) :: s
    s = 0.0_dp
    do ib = 1, nown
       b = own_list(ib)
       call set_block(b)
       fa => vptr(ca);  fb => vptr(cb)
       s = s + sum( fa(is:ie,js:je,ks:ke) * fb(is:ie,js:je,ks:ke) )
    end do
    vdot = global_sum(s)
  end function vdot

  function vptr(code) result(g)
    integer, intent(in) :: code
    real(dp), pointer :: g(:,:,:)
    select case (code)
    case (1);  g => kr
    case (2);  g => kr0
    case (4);  g => kt
    case (5);  g => ksv
    case (6);  g => kv2
    case (8);  g => krr
    case default;  g => kv1
    end select
  end function vptr

  ! infinity norm over owned cells: 0 = rhsp, 1 = kr
  real(dp) function vnorm_inf(which)
    integer, intent(in) :: which
    real(dp) :: m
    m = 0.0_dp
    do ib = 1, nown
       b = own_list(ib)
       call set_block(b)
       select case (which)
       case (0);  m = max(m, maxval(abs(rhsp(is:ie,js:je,ks:ke))))
       case default
          if (maxval(abs(kr(is:ie,js:je,ks:ke))) > m) then
             m = maxval(abs(kr(is:ie,js:je,ks:ke)))
             poi_res_loc(1) = dom(b)%b
             poi_res_loc(2:4) = maxloc(abs(kr(is:ie,js:je,ks:ke))) + [is, js, ks] - 1
          end if
       end select
    end do
    if (which /= 0) poi_res_local = m / max(bnorm, tiny(1.0_dp))
    vnorm_inf = global_max(m)
  end function vnorm_inf

  ! all-Neumann problem: pin the pressure level at the reference cell
  subroutine pin_level
    pshift = 0.0_dp
    if (dom(p_ref_dom)%owner == myrank)                                  &
       pshift = dom(p_ref_dom)%p(p_ref_i, p_ref_j, p_ref_k)               &
              - p_ref_val/(rho1*Uch*Uch)
    pshift = global_sum(pshift)               ! only the owner contributes
    do ib = 1, nown
       b = own_list(ib)                       ! b = sub-domain id
       dom(b)%p(dom(b)%r1(1):dom(b)%r2(1),                                &
                dom(b)%r1(2):dom(b)%r2(2),                                &
                dom(b)%r1(3):dom(b)%r2(3)) =                              &
       dom(b)%p(dom(b)%r1(1):dom(b)%r2(1),                                &
                dom(b)%r1(2):dom(b)%r2(2),                                &
                dom(b)%r1(3):dom(b)%r2(3)) - pshift
    end do
  end subroutine pin_level
end subroutine solve_poisson

!-----------------------------------------------------------------------
! Velocity correction of the transformed face components AND of the
! Cartesian face velocities (both directly with the same density-
! weighted pressure gradient; the solved interior Cartesian values are
! never rebuilt by interpolation).  Physical boundary faces keep the
! star (BC) values.
!-----------------------------------------------------------------------
subroutine correct_velocity
  use valves, only: valve_mob
  use global_vars
  use mpi_vars
  use b2b_runtime, only: refresh_uvw
  use curvi, only: contra_to_cart
  use mom_ops, only: face_range, face_active
  use poisson_ops, only: face_grad, Ql, Qavg, fixed_face
  use overset, only: ovs_update_faces
  implicit none
  integer :: i, j, k, ib, b, q, lo(3), hi(3), pt(3), H(3), e(3), ft
  real(dp) :: dp1, dp2, dp3
  real(dp), pointer :: Uq(:,:,:), Usq(:,:,:), uc(:,:,:), usc(:,:,:), G(:,:,:)

  ! transformed pressure-gradient fluxes of the converged p (exchanged)
  call face_grad(4)

  do ib = 1, nown
     b = own_list(ib)
     call set_block(b)
     do q = 1, 3
        select case (q)
        case (1);  Uq => U1;  Usq => Us1;  uc => u;  usc => us
        case (2);  Uq => U2;  Usq => Us2;  uc => v;  usc => vs
        case default;  Uq => U3;  Usq => Us3;  uc => w;  usc => ws
        end select
        G => field_ptr(b, 10 + q)
        e = 0;  e(q) = 1
        call face_range(q, lo, hi)
        do k = lo(3), hi(3)
         do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             pt = [i, j, k]
             ft = face_type(q, pt(q), pt)
             if (fixed_face(ft)) then
                ! fixed-flux physical boundary: keep the star (BC) value
                Uq(i,j,k) = Usq(i,j,k);  uc(i,j,k) = usc(i,j,k)
                cycle
             end if
             if (ft == 0 .and. .not. face_active(q, i, j, k)) cycle   ! owner corrects
             H = 2*pt;  H(q) = H(q) + 1
             ! transformed component: the same flux the Poisson operator used
             Uq(i,j,k) = Usq(i,j,k) - dt*fvw(H(1),H(2),H(3))*G(i,j,k)/Jh(H(1),H(2),H(3))
             ! Cartesian component: chain-rule gradient from the balanced
             ! directional differences Q_l (own-face density first, then
             ! corner average)
             dp1 = pbal(1);  dp2 = pbal(2);  dp3 = pbal(3)
             uc(i,j,k) = usc(i,j,k) - dt*valve_mob(q, pt)*( mq(q,1,H(1),H(2),H(3))*dp1      &
                       + mq(q,2,H(1),H(2),H(3))*dp2 + mq(q,3,H(1),H(2),H(3))*dp3 )     ! valve mobility (H2)
          end do
         end do
        end do
     end do
  end do

  if (poi_verbose) call outlet_balance
  call refresh_uvw(1)                    ! exchange U1,U2,U3 + BCs
  call contra_to_cart(.false., .true.)   ! ghost/BC Cartesian values only:
                                         ! the solved interior components
                                         ! were corrected directly above
  if (ovst_on) then
     call ovs_update_faces(1)            ! fringe faces <- corrected donors
     call refresh_uvw(1)                 ! ... and their halo / b2b copies
     call contra_to_cart(.false., .true.)
  end if
contains
  ! diagnostic (poi_verbose): volume fluxes through the fixed-flux
  ! (inlet/wall) faces and the pressure-outlet faces, star and corrected
  subroutine outlet_balance
    real(dp) :: fin, fout_s, fout, sgn, tmp(3)
    integer :: ierr
    fin = 0.0_dp;  fout_s = 0.0_dp;  fout = 0.0_dp
    do ib = 1, nown
       b = own_list(ib);  call set_block(b)
       do q = 1, 3
          select case (q)
          case (1);  Uq => U1;  Usq => Us1
          case (2);  Uq => U2;  Usq => Us2
          case default;  Uq => U3;  Usq => Us3
          end select
          call face_range(q, lo, hi)
          do k = lo(3), hi(3)
           do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               pt = [i, j, k];  ft = face_type(q, pt(q), pt)
               if (ft == 0) cycle
               H = 2*pt;  H(q) = H(q) + 1
               sgn = merge(-1.0_dp, 1.0_dp, pt(q) == 0)      ! outward normal
               if (fixed_face(ft)) then
                  fin = fin - sgn*Jh(H(1),H(2),H(3))/fvw(H(1),H(2),H(3))*Uq(i,j,k)
               else if (ft == 2) then
                  fout_s = fout_s + sgn*Jh(H(1),H(2),H(3))/fvw(H(1),H(2),H(3))*Usq(i,j,k)
                  fout   = fout   + sgn*Jh(H(1),H(2),H(3))/fvw(H(1),H(2),H(3))*Uq(i,j,k)
               end if
            end do
           end do
          end do
       end do
    end do
    tmp = [fin, fout_s, fout]
    call MPI_Allreduce(MPI_IN_PLACE, tmp, 3, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    if (myrank == 0) write(*,'(A,ES14.6,A,ES14.6,A,ES14.6,A,F8.4)') '    mass balance: inflow(fixed faces) ', tmp(1), &
       '  outlet star ', tmp(2), '  outlet corrected ', tmp(3), '  star/inflow-1 = ', tmp(2)/tmp(1)-1.0_dp
  end subroutine outlet_balance

  ! (1/rho) d(p)/d(xi_l) at the face between pt and pt+e (l = q: the
  ! face's own balanced difference, l /= q: corner average of Q_l)
  real(dp) function pbal(l)
    integer, intent(in) :: l
    if (l == q) then
       pbal = Ql(p, pt, q)
    else
       pbal = Qavg(p, pt, q, l)
    end if
  end function pbal
end subroutine correct_velocity

!-----------------------------------------------------------------------
! Incremental pressure form (Prs_form = 2, update 16).  The SAME
! gradient operators as the corrector (G^q fluxes for the transformed
! components, balanced chain-rule differences for the Cartesian ones):
!   mode 1 : add - (1/rho) grad p^n to the momentum right-hand sides
!            (rhsu/rhsv/rhsw, after compute_rhs_momentum) - the
!            predictor then produces the increment-form u*
!   mode 2 : add + dt (1/rho) grad p^n to the star fields (us, vs, ws and
!            Us1..3, after solve_momentum), so that the unchanged Poisson
!            equation / corrector (full pressure p^(n+1)) yield
!            u = u* - dt (1/rho) grad (p^(n+1) - p^n): the projection of
!            the increment.
! Fixed-flux faces (Dirichlet velocity) are left alone, like the
! corrector does.
!-----------------------------------------------------------------------
subroutine add_pgrad(mode)
  use global_vars
  use mpi_vars
  use b2b_runtime, only: refresh_uvw
  use curvi, only: contra_to_cart
  use mom_ops, only: face_range, face_active
  use poisson_ops, only: face_grad, Ql, Qavg, fixed_face
  use overset, only: ovs_update_faces
  implicit none
  integer, intent(in) :: mode
  integer :: i, j, k, ib, b, q, lo(3), hi(3), pt(3), H(3), e(3), ft
  real(dp) :: dp1, dp2, dp3, fac
  real(dp), pointer :: Usq(:,:,:), usc(:,:,:), G(:,:,:), r(:,:,:)
  if (Prs_form /= 2) return
  if (mode == 2) call face_grad(4)          ! transformed gradient fluxes of p^n (uses rhsu..rhsw as carriers)
  do ib = 1, nown
     b = own_list(ib)
     call set_block(b)
     do q = 1, 3
        select case (q)
        case (1);  Usq => Us1;  usc => us;  r => rhsu
        case (2);  Usq => Us2;  usc => vs;  r => rhsv
        case default;  Usq => Us3;  usc => ws;  r => rhsw
        end select
        if (mode == 2) G => field_ptr(b, 10 + q)
        e = 0;  e(q) = 1
        call face_range(q, lo, hi)
        do k = lo(3), hi(3)
         do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             pt = [i, j, k]
             ft = face_type(q, pt(q), pt)
             if (fixed_face(ft)) cycle
             if (ft == 0 .and. .not. face_active(q, i, j, k)) cycle
             if (fstat(q, i, j, k) /= 1) cycle
             H = 2*pt;  H(q) = H(q) + 1
             dp1 = pb(1);  dp2 = pb(2);  dp3 = pb(3)
             if (mode == 1) then
                r(i,j,k) = r(i,j,k) - ( mq(q,1,H(1),H(2),H(3))*dp1 + mq(q,2,H(1),H(2),H(3))*dp2 + mq(q,3,H(1),H(2),H(3))*dp3 )
             else
                fac = dt
                Usq(i,j,k) = Usq(i,j,k) + fac*fvw(H(1),H(2),H(3))*G(i,j,k)/Jh(H(1),H(2),H(3))
                usc(i,j,k) = usc(i,j,k) + fac*( mq(q,1,H(1),H(2),H(3))*dp1 + mq(q,2,H(1),H(2),H(3))*dp2 &
                           + mq(q,3,H(1),H(2),H(3))*dp3 )
             end if
          end do
         end do
        end do
     end do
  end do
  if (mode == 2) then
     call refresh_uvw(5)
     call contra_to_cart(.true., .true.)
     if (ovst_on) call ovs_update_faces(5)
  end if
contains
  real(dp) function pb(l)
    integer, intent(in) :: l
    if (l == q) then
       pb = Ql(p, pt, q)
    else
       pb = Qavg(p, pt, q, l)
    end if
  end function pb
end subroutine add_pgrad

!-----------------------------------------------------------------------
real(dp) function max_divergence()
  use global_vars
  use mpi_vars
  implicit none
  integer  :: i, j, k, ib, b
  real(dp) :: div, dmax

  dmax = 0.0_dp
  do ib = 1, nown
     b = own_list(ib)
     call set_block(b)
     do k = ks, ke
      do j = js, je
       do i = is, ie
          div = ( Jh(2*i+1,2*j,2*k)/fvw(2*i+1,2*j,2*k)*U1(i,j,k) - Jh(2*i-1,2*j,2*k)/fvw(2*i-1,2*j,2*k)*U1(i-1,j,k)  &
                + Jh(2*i,2*j+1,2*k)/fvw(2*i,2*j+1,2*k)*U2(i,j,k) - Jh(2*i,2*j-1,2*k)/fvw(2*i,2*j-1,2*k)*U2(i,j-1,k)  &
                + Jh(2*i,2*j,2*k+1)/fvw(2*i,2*j,2*k+1)*U3(i,j,k) - Jh(2*i,2*j,2*k-1)/fvw(2*i,2*j,2*k-1)*U3(i,j,k-1)  &
                ) / Jh(2*i,2*j,2*k)
          dmax = max(dmax, abs(div))
       end do
      end do
     end do
  end do
  max_divergence = global_max(dmax)
end function max_divergence

!=======================================================================
! Global mass-flux correction of the outflow faces (velocity flag 2 =
! face type 3): before the projection the star fluxes of all type-3
! faces are scaled so that their total outflow equals the net inflow
! through the fixed-flux (type-1) faces, i.e. inlets, walls and
! zero-gradient sides.  Pressure outlets (type 2) are left to the
! projection.  Both the transformed star component and the Cartesian
! star components at those faces are scaled (same face velocity
! vector).  Fluxes are the volume fluxes J U^q / fvw.
!=======================================================================
subroutine outflow_correction(qin, qout, scale)
  use global_vars
  use mpi_vars
  use mom_ops, only: face_range
  implicit none
  real(dp), intent(out) :: qin, qout, scale
  integer  :: ib, d, q, i, j, k, lo(3), hi(3), pt(3), H(3), ft
  real(dp) :: sgn, fl, s_in, s_out
  real(dp), pointer :: Usq(:,:,:), usc(:,:,:)

  s_in = 0.0_dp;  s_out = 0.0_dp
  do ib = 1, nown
     d = own_list(ib)
     call set_block(d)
     do q = 1, 3
        select case (q)
        case (1);  Usq => Us1
        case (2);  Usq => Us2
        case default;  Usq => Us3
        end select
        call face_range(q, lo, hi)
        do k = lo(3), hi(3)
         do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             pt = [i, j, k]
             if (pt(q) /= 0 .and. pt(q) /= ncell(cur_blk, q)) cycle   ! block boundary faces only
             ft = face_type(q, pt(q), pt)
             if (ft /= 1 .and. ft /= 3) cycle
             H = 2*pt;  H(q) = H(q) + 1
             sgn = merge(-1.0_dp, 1.0_dp, pt(q) == 0)       ! outward sense
             fl  = sgn * Jh(H(1),H(2),H(3))/fvw(H(1),H(2),H(3)) * Usq(i,j,k)
             if (ft == 1) then
                s_in = s_in - fl                             ! net inflow
             else
                s_out = s_out + fl                           ! outflow
             end if
          end do
         end do
        end do
     end do
  end do
  ! cavitation (update 15): the volume created inside the domain by the
  ! mass transfer leaves through the outflow faces, qout = qin + sum S V
  if (cav_on) then
     do ib = 1, nown
        d = own_list(ib)
        call set_block(d)
        do k = ks, ke
         do j = js, je
          do i = is, ie
             if (cstat(i,j,k) /= 1) cycle
             s_in = s_in + cavs(i,j,k)*Jh(2*i,2*j,2*k)
          end do
         end do
        end do
     end do
  end if
  qin  = global_sum(s_in)
  qout = global_sum(s_out)
  scale = 1.0_dp
  if (qout > 1.0e-12_dp*max(abs(qin), 1.0_dp) .and. qin > 0.0_dp) scale = qin / qout
  if (scale == 1.0_dp) return

  do ib = 1, nown
     d = own_list(ib)
     call set_block(d)
     do q = 1, 3
        select case (q)
        case (1);  Usq => Us1;  usc => us
        case (2);  Usq => Us2;  usc => vs
        case default;  Usq => Us3;  usc => ws
        end select
        call face_range(q, lo, hi)
        do k = lo(3), hi(3)
         do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             pt = [i, j, k]
             if (pt(q) /= 0 .and. pt(q) /= ncell(cur_blk, q)) cycle
             if (face_type(q, pt(q), pt) /= 3) cycle
             Usq(i,j,k) = scale*Usq(i,j,k)
             usc(i,j,k) = scale*usc(i,j,k)
          end do
         end do
        end do
     end do
  end do
end subroutine outflow_correction

!-----------------------------------------------------------------------
! SW_Quasi_2D = 3 (update 24c): a strictly two-dimensional solution on
! an extruded quasi-2-D grid (the span along the grid index k, every
! block a single sub-domain in k: k_split = 1).  After the velocity
! correction the spanwise velocity is set to zero (w and U3 on every
! k-face) and the other fields are averaged across the span: u / U1 and
! v / U2 in the face columns whose faces are all solved (fstat = 1), the
! pressure and the scalars in the cell columns whose cells are all solved.
! On an extruded grid with slip span walls (U3 = 0 there) the averaged
! fluxes keep the discrete continuity of the column: the sum of the
! U3-differences over a column vanishes.  (24c: a thin-span laminar run
! of FSI-PfS-2a still grew a spanwise mode from the cylinder's separation
! line; the line-direction choice of the preconditioner was tried
! against the thin span and reverted - see face_coefs.)
!-----------------------------------------------------------------------
subroutine quasi2d_project(istep_q2)
  use global_vars
  use mpi_vars
  use b2b_runtime, only: refresh_uvw, refresh_all
  use curvi, only: contra_to_cart
  use overset, only: ovs_update_faces
  implicit none
  integer :: ib, b, i, j, k, q, lo(3), hi(3), nk_, ierr
  real(dp) :: q2s1, q2s2, q2tol, q2xr
  integer, intent(in) :: istep_q2
  logical :: all_on
  logical, save :: checked = .false.
  real(dp), pointer :: Fq(:,:,:), Gq(:,:,:)
  if (.not. checked) then
     ierr = 0
     do ib = 1, nown
        b = own_list(ib)
        call set_block(b)
        if (ks /= 1 .or. ke /= nz) ierr = 1                       ! k_split = 1
        ! extruded: the in-plane coordinates of the half lattice do not change along k
        ! (1e-3 of the extent: a morphed grid restored from a restart file can differ slightly between the
        ! k planes - the structure's span nodes are solved independently)
        q2tol = 1.0e-3_dp*max(maxval(abs(Xh(1:2,:,:,:))), 1.0e-30_dp)
        do k = lbound(Xh,4)+1, ubound(Xh,4)
           q2xr = maxval(abs(Xh(1:2,:,:,k) - Xh(1:2,:,:,lbound(Xh,4))))
           if (q2xr > q2tol .and. k > 2*ks - 1 .and. k < 2*ke + 1) then
              ierr = max(ierr, 2)
              write(*,'(A,I0,A,ES10.3,A,ES10.3)') '  SW_Quasi_2D check: domain ', b, ': in-plane deviation between k planes ', q2xr, ' > ', q2tol
           end if
        end do
     end do
     call MPI_Allreduce(MPI_IN_PLACE, ierr, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, i)
     if (ierr == 2 .and. istep_q2 > 1) then
        ! a run restarted on a morphed grid (the first call of this run is not the first step): the
        ! deviation is reported, not fatal
        if (myrank == 0) write(*,'(A)') '  NOTE: SW_Quasi_2D - the restored (morphed) grid differs between its k planes (see above)'
        ierr = 0
     end if
     if (ierr > 0) then
        if (myrank == 0) then
           if (ierr == 1) write(*,'(A)') ' ERROR: SW_Quasi_2D needs every block undivided in k (k_split = 1)'
           if (ierr == 2) write(*,'(A)') ' ERROR: SW_Quasi_2D needs an extruded grid (the span along k, the same x-y in every k plane)'
        end if
        call MPI_Abort(MPI_COMM_WORLD, 73, i)
     end if
     checked = .true.
  end if
  do ib = 1, nown
     b = own_list(ib)
     call set_block(b)
     nk_ = ke - ks + 1
     ! the in-plane face families: u / U1 (q = 1), v / U2 (q = 2)
     do q = 1, 2
        if (q == 1) then
           Fq => u;  Gq => U1
        else
           Fq => v;  Gq => U2
        end if
        call face_range(q, lo, hi)
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              all_on = .true.
              do k = ks, ke
                 if (fstat(q, i, j, k) /= 1) all_on = .false.
              end do
              if (.not. all_on) cycle
              q2s1 = sum(Fq(i, j, ks:ke))/nk_;  q2s2 = sum(Gq(i, j, ks:ke))/nk_
              Fq(i, j, ks:ke) = q2s1;  Gq(i, j, ks:ke) = q2s2
           end do
        end do
     end do
     w(:, :, ks-1:ke) = 0.0_dp;  U3(:, :, ks-1:ke) = 0.0_dp
     ! cell fields
     do j = js, je
        do i = is, ie
           all_on = .true.
           do k = ks, ke
              if (cstat(i, j, k) /= 1) all_on = .false.
           end do
           if (.not. all_on) cycle
           p(i, j, ks:ke) = sum(p(i, j, ks:ke))/nk_
           if (allocated(dom(b)%tk))   tk(i, j, ks:ke) = sum(tk(i, j, ks:ke))/nk_
           if (allocated(dom(b)%tw))   tw(i, j, ks:ke) = sum(tw(i, j, ks:ke))/nk_
           if (allocated(dom(b)%mut))  mut(i, j, ks:ke) = sum(mut(i, j, ks:ke))/nk_
           if (allocated(dom(b)%temp)) temp(i, j, ks:ke) = sum(temp(i, j, ks:ke))/nk_
           if (Num_of_Phase == 2)      alpha(i, j, ks:ke) = sum(alpha(i, j, ks:ke))/nk_
        end do
     end do
  end do
  call refresh_uvw(1)
  call contra_to_cart(.false., .true.)
  if (ovst_on) then
     call ovs_update_faces(1)
     call refresh_uvw(1)
     call contra_to_cart(.false., .true.)
  end if
  call refresh_all(4)
  if (Num_of_Phase == 2) call refresh_all(8)
  if (energy_on) call refresh_all(20)
  if (turb_on) then
     call refresh_all(14);  call refresh_all(16)
     if (TURB_model >= 2) call refresh_all(15)
  end if
contains
  pure subroutine face_range(c, lo_, hi_)
    integer, intent(in)  :: c
    integer, intent(out) :: lo_(3), hi_(3)
    integer :: n(3)
    lo_ = [is, js, ks];  hi_ = [ie, je, ke];  n = [nx, ny, nz]
    if (lo_(c) == 1)    lo_(c) = 0
    if (hi_(c) == n(c)) hi_(c) = n(c)
  end subroutine face_range
end subroutine quasi2d_project
