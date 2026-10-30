!=======================================================================
! boundary.F90
!   Boundary condition application on the fully-staggered grid.
!
!   apply_bc(f, iv), iv = field id (mod_global table):
!     iv = 1..3 : U1..U3 corrected face components   flag rows 1..3
!     iv = 5..7 : Us1..Us3 star (predictor) components, same rows
!     iv = 4    : p  (cell centered)                  flag row 4
!     iv = 8    : alpha (cell centered)               flag row 5
!     iv = 14..20: turbulence, surface tension, vapour, T (cell): zero
!                 gradient (the T ghosts are completed by energy.F90)
!     iv = 9,10 : Krylov vectors of the pressure solver (cell
!                 centered, pressure flags with HOMOGENEOUS values)
!
!   Zero-gradient / outflow (flag 0 / 2) on the face-normal family: the
!   boundary face value is NEVER overwritten - it is the momentum
!   solution of that face (the face is in the sweep set), corrected by
!   the projection on a pressure outlet or scaled by the global mass-
!   flux correction on an outflow face; only the two outer ghost faces
!   are filled by copying it (2nd-order stencils of the boundary face).
!   (Until update 4 the STAR field had the boundary face copied from the
!   first interior face; the projection then had to remove the
!   difference between that extrapolated flux and the momentum-
!   consistent one every step, which costs rho du dx/(2 dt) in the
!   outlet cell and left a uniform pressure offset of kPa size on
!   grids with large outlet cells - the "pressure level" and "outlet
!   wiggle" problems of the NACA66mod tutorial cases.)
!
!   Velocity families: the flag that applies to the transformed family
!   U^q is the EFFECTIVE flag bc_eflag(q) (the .bcs flag of the
!   Cartesian axis the grid direction q points along on that patch,
!   see read_bcs/effective_flags).  Dirichlet (flag 1): the .bcs values
!   are the global Cartesian components; the transformed face value is
!   fvw * (m . u_bc(t)) at the face with the start-up ramp
!   u_bc(t) = u_init + bc_ramp*(u_bc - u_init), bc_ramp = min(t/t_ramp,1)
!   (t_ramp in the namelist; bc_ramp = 1 when no ramp).  With u_init =
!   the free-stream speed this emulates a towed body (steady far field,
!   wall velocity decelerating from the free stream to rest); with
!   u_init = 0 it is a ramped inflow.
!     - field staggered NORMAL to the patch: the boundary face value is
!       set directly, the outer ghost face is linearly extrapolated;
!     - otherwise (cell-like across the patch): velocity-space mirror
!       f_ghost = 2*val - bw*f_int (bw = J_ghost/J_int for curvi_form 2).
!   Pressure (flag 1): mirror 2*val - f_int.
!   Volume fraction: flag 1 -> ghost = value (no mirror, stays in [0,1]);
!                    flag 2 -> ghost = still-water fill fraction of the
!                              ghost cell (inlet profile).
!   Neumann / zero-gradient (flag 0): ghost = mirrored interior value.
!   Incident wave (flag 12, update 18, waves.F90): velocity = ramped mean
!   + wave orbital velocity at the face (Dirichlet handling), alpha ghost
!   = fill fraction below the incident surface, pressure ghost = interior
!   + exact wave-pressure difference (homogeneous: zero gradient).
!
!   The patch tangential extent is widened by ng cells before clipping
!   to the local array bounds so that domain-edge/corner ghost values
!   are also filled.
!=======================================================================
subroutine apply_bc(f, iv)
  use global_vars
  use mpi_vars
  use vof, only: fill_fraction
  use sixdof, only: sixdof_wall_velocity
  use waves, only: wave_vel, wave_fill, wave_pdyn, wave_mean_vel, wave_t
  implicit none
  real(dp), intent(inout) :: f(is-ng:ie+ng, js-ng:je+ng, ks-ng:ke+ng)
  integer,  intent(in)    :: iv

  integer :: n, d, flag, row, ivar
  integer :: i1, i2, j1, j2, k1, k2      ! loop ranges (global indices)
  integer :: i, j, k
  logical :: homog, star, wave

  select case (iv)
  case (1:4);  ivar = iv;      row = iv
  case (5:7);  ivar = iv - 4;  row = ivar
  case (8);    ivar = 8;       row = 5
  case (14:20); ivar = iv;     row = 0   ! turbulence / surface-tension / vapour / temperature fields: zero gradient on
                                         ! every physical face (T: the Dirichlet / flux / convective ghosts are then set by
                                         ! energy.F90 through temp_bc_hook; update 16 - before, 19 fell through to the
                                         ! pressure flags and read the metric out of bounds at Dirichlet-p faces)
  case default; ivar = iv;     row = 4   ! Krylov vectors: pressure flags
  end select
  homog = (iv >= 9 .and. iv <= 13)
  star  = (iv >= 5 .and. iv <= 7)

  do n = 1, nbc
     if (ivar <= 3) then
        flag = bc_eflag(ivar, n)          ! transformed family: effective flag
        if (flag == 2) flag = 0           ! outflow = zero gradient here
     else if (row == 0) then
        flag = 0                          ! turbulence: zero gradient (walls / inlets are set by turbulence.F90)
     else
        flag = bc_flag(row, n)
     end if
     if (flag == 11) flag = 1             ! wall function patch: no-slip Dirichlet here
     wave = (flag == 12)                  ! incident-wave patch (update 18)
     if (wave .and. ivar <= 3) flag = 1   !   velocity: Dirichlet mean + wave (bv / mirror)
     d    = bc_dir(n)

     ! ---- tangential loop ranges: node range -> index range of ivar ----
     call tan_range(1, ivar, bc_irng(1,n), bc_irng(2,n), ni, is, ie, i1, i2)
     call tan_range(2, ivar, bc_jrng(1,n), bc_jrng(2,n), nj, js, je, j1, j2)
     call tan_range(3, ivar, bc_krng(1,n), bc_krng(2,n), nk, ks, ke, k1, k2)

     select case (d)

     case (1)                              ! ---- i-min face ----
        if (is /= 1) cycle
        if (ivar == 1) then                ! U1 : normal staggered
           do k = k1, k2;  do j = j1, j2
              if (flag == 1) then
                 f(0,j,k)  = bv(0,j,k)
                 f(-1,j,k) = mirror(-1,j,k, 1,j,k)
              else
                 ! zero gradient / outflow: the boundary face value is the
                 ! momentum solution of that face (the projection adds the
                 ! pressure gradient of a pressure outlet, outflow faces are
                 ! mass-flux corrected); only the ghost faces are filled
                 f(-1,j,k) = f(0,j,k)
                 f(-2,j,k) = f(0,j,k)
              end if
           end do;  end do
        else                               ! cell-like in i
           do k = k1, k2;  do j = j1, j2
              f(0,j,k)  = cellval(flag, 0,j,k, 1,j,k)
              f(-1,j,k) = cellval(flag, -1,j,k, 2,j,k)
           end do;  end do
        end if

     case (-1)                             ! ---- i-max face ----
        if (ie /= nx) cycle
        if (ivar == 1) then
           do k = k1, k2;  do j = j1, j2
              if (flag == 1) then
                 f(nx,j,k)   = bv(nx,j,k)
                 f(nx+1,j,k) = mirror(nx+1,j,k, nx-1,j,k)
              else
                 f(nx+1,j,k) = f(nx,j,k)
                 f(nx+2,j,k) = f(nx,j,k)
              end if
           end do;  end do
        else
           do k = k1, k2;  do j = j1, j2
              f(nx+1,j,k) = cellval(flag, nx+1,j,k, nx,j,k)
              f(nx+2,j,k) = cellval(flag, nx+2,j,k, nx-1,j,k)
           end do;  end do
        end if

     case (2)                              ! ---- j-min face ----
        if (js /= 1) cycle
        if (ivar == 2) then                ! U2 : normal staggered
           do k = k1, k2;  do i = i1, i2
              if (flag == 1) then
                 f(i,0,k)  = bv(i,0,k)
                 f(i,-1,k) = mirror(i,-1,k, i,1,k)
              else
                 f(i,-1,k) = f(i,0,k)
                 f(i,-2,k) = f(i,0,k)
              end if
           end do;  end do
        else
           do k = k1, k2;  do i = i1, i2
              f(i,0,k)  = cellval(flag, i,0,k, i,1,k)
              f(i,-1,k) = cellval(flag, i,-1,k, i,2,k)
           end do;  end do
        end if

     case (-2)                             ! ---- j-max face ----
        if (je /= ny) cycle
        if (ivar == 2) then
           do k = k1, k2;  do i = i1, i2
              if (flag == 1) then
                 f(i,ny,k)   = bv(i,ny,k)
                 f(i,ny+1,k) = mirror(i,ny+1,k, i,ny-1,k)
              else
                 f(i,ny+1,k) = f(i,ny,k)
                 f(i,ny+2,k) = f(i,ny,k)
              end if
           end do;  end do
        else
           do k = k1, k2;  do i = i1, i2
              f(i,ny+1,k) = cellval(flag, i,ny+1,k, i,ny,k)
              f(i,ny+2,k) = cellval(flag, i,ny+2,k, i,ny-1,k)
           end do;  end do
        end if

     case (3)                              ! ---- k-min face ----
        if (ks /= 1) cycle
        if (ivar == 3) then                ! U3 : normal staggered
           do j = j1, j2;  do i = i1, i2
              if (flag == 1) then
                 f(i,j,0)  = bv(i,j,0)
                 f(i,j,-1) = mirror(i,j,-1, i,j,1)
              else
                 f(i,j,-1) = f(i,j,0)
                 f(i,j,-2) = f(i,j,0)
              end if
           end do;  end do
        else
           do j = j1, j2;  do i = i1, i2
              f(i,j,0)  = cellval(flag, i,j,0, i,j,1)
              f(i,j,-1) = cellval(flag, i,j,-1, i,j,2)
           end do;  end do
        end if

     case (-3)                             ! ---- k-max face ----
        if (ke /= nz) cycle
        if (ivar == 3) then
           do j = j1, j2;  do i = i1, i2
              if (flag == 1) then
                 f(i,j,nz)   = bv(i,j,nz)
                 f(i,j,nz+1) = mirror(i,j,nz+1, i,j,nz-1)
              else
                 f(i,j,nz+1) = f(i,j,nz)
                 f(i,j,nz+2) = f(i,j,nz)
              end if
           end do;  end do
        else
           do j = j1, j2;  do i = i1, i2
              f(i,j,nz+1) = cellval(flag, i,j,nz+1, i,j,nz)
              f(i,j,nz+2) = cellval(flag, i,j,nz+2, i,j,nz-1)
           end do;  end do
        end if

     end select
  end do

contains

  !--------------------------------------------------------------------
  ! Ghost value at (g1,g2,g3) of a cell-like variable, mirrored about
  ! the boundary from the interior position (b1,b2,b3).
  !--------------------------------------------------------------------
  real(dp) function cellval(fl, g1, g2, g3, b1, b2, b3)
    integer, intent(in) :: fl, g1, g2, g3, b1, b2, b3
    select case (fl)
    case (1)
       if (ivar == 8) then
          cellval = bc_val(5, n)
       else
          cellval = mirror(g1,g2,g3, b1,b2,b3)
       end if
    case (2)
       cellval = fill_fraction(g1, g2, g3)          ! alpha only (checked)
    case (12)                                       ! incident wave (update 18)
       if (ivar == 8) then
          cellval = wave_fill(g1, g2, g3, wave_t)   ! alpha: fill below the incident surface
       else if (ivar == 4 .and. .not. homog) then   ! p: interior + exact wave-pressure difference
          cellval = f(b1,b2,b3) + wave_pdyn([2*g1,2*g2,2*g3], wave_t) - wave_pdyn([2*b1,2*b2,2*b3], wave_t)
       else
          cellval = f(b1,b2,b3)                     ! Krylov vectors / others: zero gradient
       end if
    case default
       cellval = f(b1,b2,b3)
    end select
  end function cellval

  !--------------------------------------------------------------------
  ! Boundary value of variable ivar at the staggered array position
  ! (a1,a2,a3).  For the velocity families the .bcs Dirichlet values
  ! are the CARTESIAN components (u,v,w); the face variable is the
  ! transformed component, so the value is projected pointwise with
  ! the local m matrix:  val_q = fvw * sum_c m(c,q) * u_bc(c).
  ! Pressure keeps its scalar value; Krylov vectors get 0.
  !--------------------------------------------------------------------
  real(dp) function bv(a1, a2, a3)
    integer, intent(in) :: a1, a2, a3
    integer :: hx, hy, hz
    real(dp) :: ub(3)
    if (homog) then
       bv = 0.0_dp
       return
    end if
    if (ivar == 4) then
       bv = bc_val(4, n)
       return
    end if
    hx = 2*a1;  hy = 2*a2;  hz = 2*a3
    select case (ivar)
    case (1);  hx = hx + 1
    case (2);  hy = hy + 1
    case (3);  hz = hz + 1
    end select
    ! start-up ramp: the Dirichlet value moves from the initial velocity
    ! to the .bcs value,  u(t) = u_init + bc_ramp*(u_bc - u_init);
    ! on the blocks of a moving body the grid velocity is added (the
    ! no-slip wall moves with the body, sixdof.F90)
    ub(1) = u_init/Uch + bc_ramp*(bc_val(1,n) - u_init/Uch)
    ub(2) = v_init/Uch + bc_ramp*(bc_val(2,n) - v_init/Uch)
    ub(3) = w_init/Uch + bc_ramp*(bc_val(3,n) - w_init/Uch)
    if (wave) ub = wave_mean_vel() + wave_vel(Xh(:,hx,hy,hz), wave_t)   ! incident wave at the face (update 18)
    if (any_moving) ub = ub + sixdof_wall_velocity(cur_blk, Xh(:,hx,hy,hz), [hx, hy, hz])
    bv = fvw(hx,hy,hz) * ( mq(1,ivar,hx,hy,hz)*ub(1) + mq(2,ivar,hx,hy,hz)*ub(2) + mq(3,ivar,hx,hy,hz)*ub(3) )
  end function bv

  !--------------------------------------------------------------------
  ! Dirichlet ghost value at the staggered position (g1,g2,g3) mirrored
  ! about the boundary from the interior position (b1,b2,b3).
  !   Velocity families (ivar 1..3): the CARTESIAN velocity vector is
  !   mirrored, u_g = 2 u_bc - u(b), with u(b) assembled at the interior
  !   face exactly as cart_to_contra does (own component from the
  !   family's Cartesian array, the two foreign components averaged from
  !   the neighbouring faces of their families), and projected with the
  !   metrics of the ghost position: U_g = fvw_g m_g . u_g.  (Until
  !   update 4 the mirror was done in transformed space with the
  !   extrapolated ghost metrics, 2 bv_g - (J_g/J_b) U(b); on curved or
  !   stretched boundaries - the O-grid outer circle, the C-grid arc -
  !   the ghost metrics differ from the interior ones by percents, the
  !   tangential Cartesian ghost velocity was off by the same amount,
  !   and the momentum convected through the inlet with it: a pressure
  !   wiggle in the first cells behind every Dirichlet inlet.)
  !   Cell fields: the old formula (exact).
  !--------------------------------------------------------------------
  real(dp) function mirror(g1, g2, g3, b1, b2, b3)
    integer, intent(in) :: g1, g2, g3, b1, b2, b3
    real(dp) :: uc(3), ub(3)
    integer :: hx, hy, hz
    if (ivar > 3) then
       mirror = 2.0_dp*bv(g1,g2,g3) - bw(g1,g2,g3, b1,b2,b3)*f(b1,b2,b3)
       return
    end if
    hx = 2*g1;  hy = 2*g2;  hz = 2*g3
    select case (ivar)
    case (1);  hx = hx + 1
    case (2);  hy = hy + 1
    case (3);  hz = hz + 1
    end select
    if (homog) then
       ub = 0.0_dp
    else
       ub(1) = u_init/Uch + bc_ramp*(bc_val(1,n) - u_init/Uch)
       ub(2) = v_init/Uch + bc_ramp*(bc_val(2,n) - v_init/Uch)
       ub(3) = w_init/Uch + bc_ramp*(bc_val(3,n) - w_init/Uch)
       ! incident wave: mean + orbital velocity at the boundary position
       ! (midway between the ghost and the interior position; update 18)
       if (wave) ub = wave_mean_vel() + wave_vel(0.5_dp*(Xh(:,hx,hy,hz) + Xh(:,2*b1+merge(1,0,ivar==1), &
                                                  2*b2+merge(1,0,ivar==2), 2*b3+merge(1,0,ivar==3))), wave_t)
       ! moving body: wall velocity at the boundary position (midway
       ! between the ghost and the interior position)
       if (any_moving) ub = ub + sixdof_wall_velocity(cur_blk, 0.5_dp*(Xh(:,hx,hy,hz) + Xh(:,2*b1+merge(1,0,ivar==1), &
                                                        2*b2+merge(1,0,ivar==2), 2*b3+merge(1,0,ivar==3))), &
                                                        ([hx, hy, hz] + [2*b1+merge(1,0,ivar==1), 2*b2+merge(1,0,ivar==2), 2*b3+merge(1,0,ivar==3)])/2)
    end if
    uc = cart_at_face(ivar, b1, b2, b3)
    mirror = fvw(hx,hy,hz) * ( mq(1,ivar,hx,hy,hz)*(2.0_dp*ub(1) - uc(1))    &
                             + mq(2,ivar,hx,hy,hz)*(2.0_dp*ub(2) - uc(2))    &
                             + mq(3,ivar,hx,hy,hz)*(2.0_dp*ub(3) - uc(3)) )
  end function mirror

  ! Cartesian velocity vector at the staggered position (a1,a2,a3) of
  ! family q (same averaging as cart_to_contra), star or final field
  function cart_at_face(q, a1, a2, a3) result(uc)
    integer, intent(in) :: q, a1, a2, a3
    real(dp) :: uc(3)
    select case (q)
    case (1)
       uc(1) = cf(1, a1, a2, a3)
       uc(2) = 0.25_dp*( cf(2,a1,a2-1,a3) + cf(2,a1+1,a2-1,a3) + cf(2,a1,a2,a3) + cf(2,a1+1,a2,a3) )
       uc(3) = 0.25_dp*( cf(3,a1,a2,a3-1) + cf(3,a1+1,a2,a3-1) + cf(3,a1,a2,a3) + cf(3,a1+1,a2,a3) )
    case (2)
       uc(2) = cf(2, a1, a2, a3)
       uc(1) = 0.25_dp*( cf(1,a1-1,a2,a3) + cf(1,a1-1,a2+1,a3) + cf(1,a1,a2,a3) + cf(1,a1,a2+1,a3) )
       uc(3) = 0.25_dp*( cf(3,a1,a2,a3-1) + cf(3,a1,a2+1,a3-1) + cf(3,a1,a2,a3) + cf(3,a1,a2+1,a3) )
    case default
       uc(3) = cf(3, a1, a2, a3)
       uc(1) = 0.25_dp*( cf(1,a1-1,a2,a3) + cf(1,a1-1,a2,a3+1) + cf(1,a1,a2,a3) + cf(1,a1,a2,a3+1) )
       uc(2) = 0.25_dp*( cf(2,a1,a2-1,a3) + cf(2,a1,a2-1,a3+1) + cf(2,a1,a2,a3) + cf(2,a1,a2,a3+1) )
    end select
  end function cart_at_face

  real(dp) function cf(c, a1, a2, a3)
    integer, intent(in) :: c, a1, a2, a3
    integer :: b1, b2, b3
    ! (the transverse average at the outermost ghost layer reaches one
    !  index beyond the storage: clamp to the array bounds)
    if (star) then
       select case (c)
       case (1);  b1 = clampi(a1, lbound(us,1), ubound(us,1)); b2 = clampi(a2, lbound(us,2), ubound(us,2)); b3 = clampi(a3, lbound(us,3), ubound(us,3)); cf = us(b1,b2,b3)
       case (2);  b1 = clampi(a1, lbound(vs,1), ubound(vs,1)); b2 = clampi(a2, lbound(vs,2), ubound(vs,2)); b3 = clampi(a3, lbound(vs,3), ubound(vs,3)); cf = vs(b1,b2,b3)
       case default;  b1 = clampi(a1, lbound(ws,1), ubound(ws,1)); b2 = clampi(a2, lbound(ws,2), ubound(ws,2)); b3 = clampi(a3, lbound(ws,3), ubound(ws,3)); cf = ws(b1,b2,b3)
       end select
    else
       select case (c)
       case (1);  b1 = clampi(a1, lbound(u,1), ubound(u,1)); b2 = clampi(a2, lbound(u,2), ubound(u,2)); b3 = clampi(a3, lbound(u,3), ubound(u,3)); cf = u(b1,b2,b3)
       case (2);  b1 = clampi(a1, lbound(v,1), ubound(v,1)); b2 = clampi(a2, lbound(v,2), ubound(v,2)); b3 = clampi(a3, lbound(v,3), ubound(v,3)); cf = v(b1,b2,b3)
       case default;  b1 = clampi(a1, lbound(w,1), ubound(w,1)); b2 = clampi(a2, lbound(w,2), ubound(w,2)); b3 = clampi(a3, lbound(w,3), ubound(w,3)); cf = w(b1,b2,b3)
       end select
    end if
  end function cf

  pure integer function clampi(a, lo, hi)
    integer, intent(in) :: a, lo, hi
    clampi = min(max(a, lo), hi)
  end function clampi

  !--------------------------------------------------------------------
  ! Mirror weight between the ghost position (a1..a3) and the interior
  ! position (b1..b3): fvw(ghost)/fvw(interior).  With it the Dirichlet
  ! mirror reproduces the velocity-space mirror for BOTH formulations
  ! (for the unweighted form and for the cell fields it is exactly 1).
  !--------------------------------------------------------------------
  real(dp) function bw(a1, a2, a3, b1, b2, b3)
    integer, intent(in) :: a1, a2, a3, b1, b2, b3
    integer :: hx, hy, hz, gx, gy, gz
    if (ivar > 3 .or. curvi_form == 1) then
       bw = 1.0_dp
       return
    end if
    gx = 2*a1;  gy = 2*a2;  gz = 2*a3
    hx = 2*b1;  hy = 2*b2;  hz = 2*b3
    select case (ivar)
    case (1);  gx = gx + 1;  hx = hx + 1
    case (2);  gy = gy + 1;  hy = hy + 1
    case (3);  gz = gz + 1;  hz = hz + 1
    end select
    bw = Jh(gx,gy,gz) / Jh(hx,hy,hz)
  end function bw

  !--------------------------------------------------------------------
  ! Tangential index range of variable ivar covered by the patch node
  ! range [n1,n2] in direction dir, widened by ng and clipped to the
  ! local sub-domain extent (including ghost layers).
  !   cell-like index range : n1   .. n2-1
  !   staggered (face) range: n1-1 .. n2-1
  !--------------------------------------------------------------------
  subroutine tan_range(dir, ivar_, n1, n2, nnode, ls, le, r1, r2)
    integer, intent(in)  :: dir, ivar_, n1, n2, nnode, ls, le
    integer, intent(out) :: r1, r2
    if (ivar_ == dir) then         ! field staggered in this direction
       r1 = n1 - 1
       r2 = n2 - 1
    else                           ! cell-like in this direction
       r1 = n1
       r2 = n2 - 1
    end if
    ! extend into ghost region only where the patch touches the global
    ! domain edge (fills corner/edge ghost values without bleeding into
    ! neighboring patches that share the same boundary face)
    if (n1 <= 1)     r1 = r1 - ng
    if (n2 >= nnode) r2 = r2 + ng
    ! clip to the local sub-domain extent (incl. ghost layers)
    r1 = max(r1, ls - ng)
    r2 = min(r2, le + ng)
  end subroutine tan_range
end subroutine apply_bc


