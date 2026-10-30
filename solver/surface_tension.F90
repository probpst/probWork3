!=======================================================================
! surface_tension.F90 -- continuum surface force (CSF, Brackbill, Kothe
! & Zemach 1992) for the two-phase solver (update 11).
!
!   Name list ($VoF): Surf_tension = 0 / 1, Sigma_ST [N/m] (default
!   0.0728).  Active only when Num_of_Phase = 2; independent of the VoF
!   scheme (the force depends on alpha only, so it works with the
!   explicit CICSAM, the implicit CICSAM and the geometric VoF).
!
!   Force per unit volume  f = sigma kappa grad(alpha), added to the
!   momentum equations as an acceleration f / rho on every velocity face
!   (compute_rhs_momentum), with the face gradient of alpha discretised
!   exactly like the face pressure gradient (balanced-force form: the
!   pressure jump across a static interface is reproduced with small
!   spurious currents).
!
!   Curvature  kappa = -div( grad(alpha_s) / |grad(alpha_s)| )  at the
!   cell centres, from alpha smoothed twice with a 27-point kernel
!   (weights 1-2-1 per axis); the unit normal is evaluated at the cell
!   centres, averaged to the faces, and its divergence is formed from the
!   contravariant fluxes  J m(:,q).n  through the six faces.  kappa is
!   set to zero where |grad(alpha_s)| is negligible (pure cells).
!
!   Non-dimensional: sigma_nd = Sigma_ST / (rho1 Uch^2 Lch) = 1 / We.
!   Capillary time-step limit  dt < sqrt((rho1 + rho2)/2 h^3 / (2 pi
!   sigma))  is evaluated on the smallest interface cell and reported
!   in the monitor (a warning when dt exceeds it).
!=======================================================================
module surface_tension
  use global_vars
  use mpi_vars
  use b2b_runtime, only: refresh_all
  implicit none
  private
  public :: compute_curvature, st_force, st_monitor

  real(dp) :: mon_kmax = 0.0_dp, mon_dtcap = huge(1.0_dp)
  integer  :: mon_nint = 0

contains

  !---------------------------------------------------------------------
  ! smoothed alpha, unit normal, curvature (once per time step)
  !---------------------------------------------------------------------
  subroutine compute_curvature
    integer  :: ib, d, i, j, k, pass, q, l, sgn, e(3), ii(3)
    real(dp) :: g(3), gn, nf(3), nq(3), flux, kmax, hmin, dtc, hcell, rhoavg
    real(dp), allocatable :: nx(:,:,:), ny(:,:,:), nz(:,:,:), tmp(:,:,:)
    real(dp), parameter :: wgt(-1:1) = [0.25_dp, 0.5_dp, 0.25_dp]
    real(dp), parameter :: eps_g = 1.0e-8_dp
    ! ---- 1. smoothing: alps = K * alpha (two passes, ghosts refreshed) --
    do ib = 1, nown
       d = own_list(ib);  call set_block(d)
       alps = alpha
    end do
    call refresh_all(17)
    do pass = 1, 2
       do ib = 1, nown
          d = own_list(ib);  call set_block(d)
          allocate(tmp(is:ie, js:je, ks:ke))
          do k = ks, ke
           do j = js, je
            do i = is, ie
               tmp(i,j,k) = 0.0_dp
               do q = -1, 1
                  do l = -1, 1
                     do sgn = -1, 1
                        tmp(i,j,k) = tmp(i,j,k) + wgt(q)*wgt(l)*wgt(sgn)*alps(i+q, j+l, k+sgn)
                     end do
                  end do
               end do
            end do
           end do
          end do
          alps(is:ie,js:je,ks:ke) = tmp
          deallocate(tmp)
       end do
       call refresh_all(17)
    end do
    ! ---- 2. unit normal at the cell centres (own cells + one ghost layer) --
    kmax = 0.0_dp;  hmin = huge(1.0_dp);  mon_nint = 0
    do ib = 1, nown
       d = own_list(ib);  call set_block(d)
       allocate(nx(is-1:ie+1, js-1:je+1, ks-1:ke+1), ny(is-1:ie+1, js-1:je+1, ks-1:ke+1), nz(is-1:ie+1, js-1:je+1, ks-1:ke+1))
       do k = ks-1, ke+1
        do j = js-1, je+1
         do i = is-1, ie+1
            g = grad_cell(alps, i, j, k)
            gn = norm2(g)
            if (gn > eps_g) then
               nx(i,j,k) = g(1)/gn;  ny(i,j,k) = g(2)/gn;  nz(i,j,k) = g(3)/gn
            else
               nx(i,j,k) = 0.0_dp;  ny(i,j,k) = 0.0_dp;  nz(i,j,k) = 0.0_dp
            end if
         end do
        end do
       end do
       ! ---- 3. curvature = - div(n): contravariant flux J m(:,q).n through the faces
       do k = ks, ke
        do j = js, je
         do i = is, ie
            kap(i,j,k) = 0.0_dp
            if (cstat(i,j,k) < 0) cycle
            g = grad_cell(alps, i, j, k)
            if (norm2(g) <= eps_g) cycle
            flux = 0.0_dp
            do q = 1, 3
               e = 0;  e(q) = 1
               do sgn = -1, 1, 2
                  ii = [i,j,k] + sgn*e
                  nf = 0.5_dp*([nx(i,j,k), ny(i,j,k), nz(i,j,k)] + [nx(ii(1),ii(2),ii(3)), ny(ii(1),ii(2),ii(3)), nz(ii(1),ii(2),ii(3))])
                  nq = [mh(1+3*(q-1), 2*i+sgn*e(1), 2*j+sgn*e(2), 2*k+sgn*e(3)), &
                        mh(2+3*(q-1), 2*i+sgn*e(1), 2*j+sgn*e(2), 2*k+sgn*e(3)), &
                        mh(3+3*(q-1), 2*i+sgn*e(1), 2*j+sgn*e(2), 2*k+sgn*e(3))]
                  flux = flux + real(sgn, dp)*Jh(2*i+sgn*e(1), 2*j+sgn*e(2), 2*k+sgn*e(3))*dot_product(nq, nf)
               end do
            end do
            kap(i,j,k) = -flux/Jh(2*i,2*j,2*k)
            if (alpha(i,j,k) > 0.01_dp .and. alpha(i,j,k) < 0.99_dp) then
               mon_nint = mon_nint + 1
               kmax = max(kmax, abs(kap(i,j,k)))
               hcell = Jh(2*i,2*j,2*k)**(1.0_dp/3.0_dp)
               hmin = min(hmin, hcell)
            end if
         end do
        end do
       end do
       deallocate(nx, ny, nz)
    end do
    call refresh_all(18)
    ! ---- monitor: max |kappa| (1/Lch), capillary time-step limit -------
    mon_kmax = global_max(kmax)
    hmin = -global_max(-hmin)
    mon_dtcap = huge(1.0_dp)
    if (hmin < huge(1.0_dp) .and. Sigma_ST > 0.0_dp) then
       rhoavg = 0.5_dp*(rho1 + rho2)
       dtc = sqrt(rhoavg*(hmin*Lch)**3/(2.0_dp*acos(-1.0_dp)*Sigma_ST))    ! [s]
       mon_dtcap = dtc
    end if
  end subroutine compute_curvature

  ! Cartesian gradient of a cell-centred field at the cell centre
  function grad_cell(f, i, j, k) result(g)
    real(dp), intent(in) :: f(lbound(alpha,1):, lbound(alpha,2):, lbound(alpha,3):)
    integer,  intent(in) :: i, j, k
    real(dp) :: g(3), dq(3)
    integer  :: l, q
    dq(1) = 0.5_dp*(f(i+1,j,k) - f(i-1,j,k))
    dq(2) = 0.5_dp*(f(i,j+1,k) - f(i,j-1,k))
    dq(3) = 0.5_dp*(f(i,j,k+1) - f(i,j,k-1))
    do l = 1, 3
       g(l) = 0.0_dp
       do q = 1, 3
          g(l) = g(l) + dq(q)*mh(l+3*(q-1), 2*i, 2*j, 2*k)
       end do
    end do
  end function grad_cell

  !---------------------------------------------------------------------
  ! surface-tension acceleration on the velocity face (c; i,j,k) at the
  ! half-index position H:  sigma_nd kappa_f (d alpha / d x_c)_f / rho_f
  ! The face gradient uses the direct difference across the face along
  ! the face axis and the average of the two cells' central differences
  ! in the tangential index directions (same stencil as the pressure
  ! gradient on the face).
  !---------------------------------------------------------------------
  real(dp) function st_force(c, i, j, k, H) result(acc)
    integer, intent(in) :: c, i, j, k, H(3)
    integer  :: e(3), q, ia(3), ib(3)
    real(dp) :: dq(3), gc, kf, rf
    acc = 0.0_dp
    e = 0;  e(c) = 1
    ia = [i, j, k];  ib = ia + e                       ! the two cells sharing the face
    if (ia(1) < lbound(alpha,1)+1 .or. ib(1) > ubound(alpha,1)-1 .or. ia(2) < lbound(alpha,2)+1 .or. &
        ib(2) > ubound(alpha,2)-1 .or. ia(3) < lbound(alpha,3)+1 .or. ib(3) > ubound(alpha,3)-1) return
    kf = 0.5_dp*(kap(ia(1),ia(2),ia(3)) + kap(ib(1),ib(2),ib(3)))
    if (kf == 0.0_dp) return
    do q = 1, 3
       if (q == c) then
          dq(q) = alpha(ib(1),ib(2),ib(3)) - alpha(ia(1),ia(2),ia(3))
       else
          dq(q) = 0.25_dp*(alpha(ia(1)+merge(1,0,q==1), ia(2)+merge(1,0,q==2), ia(3)+merge(1,0,q==3)) - &
                           alpha(ia(1)-merge(1,0,q==1), ia(2)-merge(1,0,q==2), ia(3)-merge(1,0,q==3)) + &
                           alpha(ib(1)+merge(1,0,q==1), ib(2)+merge(1,0,q==2), ib(3)+merge(1,0,q==3)) - &
                           alpha(ib(1)-merge(1,0,q==1), ib(2)-merge(1,0,q==2), ib(3)-merge(1,0,q==3)))
       end if
    end do
    gc = 0.0_dp
    do q = 1, 3
       gc = gc + dq(q)*mh(c+3*(q-1), H(1), H(2), H(3))     ! d alpha / d x_c at the face
    end do
    if (gc == 0.0_dp) return
    rf  = rho_h(H(1), H(2), H(3))
    acc = sigma_nd*kf*gc/rf
  end function st_force

  !---------------------------------------------------------------------
  subroutine st_monitor
    if (myrank /= 0 .or. .not. st_on) return
    write(*,'(A)') ' '
    write(*,'(A)') ' [ Cell ] Center Monitor :    < Surface tension (CSF) >'
    write(*,'(A,I15)')    '  --> interface cells: ', mon_nint
    write(*,'(A,ES19.4)') '  --> max |kappa| [1/m]: ', mon_kmax/Lch
    if (mon_dtcap < huge(1.0_dp)) then
       write(*,'(A,ES19.4)') '  --> capillary dt limit [s]: ', mon_dtcap
       if (dt*Lch/Uch > mon_dtcap) write(*,'(A)') '  WARNING: the time step exceeds the capillary limit (Brackbill); ' // &
          'expect spurious currents / instability at the interface'
    end if
  end subroutine st_monitor

end module surface_tension
