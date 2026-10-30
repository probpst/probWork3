!=======================================================================
! wave_fields.F90  -  incident waves (update 18): field operations that
!   need the momentum / exchange / metric modules (waves.F90 itself is
!   compiled before them so that the boundary routines can use it).
!   wave_relax_fields(tt): after the time step, blend the volume
!   fraction (cells) and the velocity (faces: transformed AND Cartesian
!   component, same weight) towards the incident wave inside the zones,
!   phi <- (1 - w) phi + w phi_wave, then refresh the ghost / interface
!   values.  Boundary faces are left to the boundary conditions, non-
!   owned interface faces and overset fringe / hole faces are skipped
!   (they receive the owner's values).
!   wave_init_fields: Wave_Init = 1 - alpha and the velocity of the whole
!   field from the exact wave at t = 0 (after initialize_fields).
!=======================================================================
subroutine wave_relax_fields(tt)
  use global_vars
  use mpi_vars
  use waves, only: wave_on, wave_zones_on, wave_weight, wave_vel, wave_fill, wave_mean_vel
  use mom_ops, only: face_range, face_active
  use b2b_runtime, only: refresh_uvw, refresh_all
  use curvi, only: contra_to_cart
  use overset, only: ovs_update_faces
  implicit none
  real(dp), intent(in) :: tt
  integer :: ib, d, q, i, j, k, lo(3), hi(3), pt(3), H(3)
  real(dp) :: wgt, um(3), uw(3), ut
  real(dp), pointer :: Uq(:,:,:), uc(:,:,:)
  if (.not. wave_on .or. .not. wave_zones_on) return
  um = wave_mean_vel()
  do ib = 1, nown
     d = own_list(ib)
     call set_block(d)
     do k = ks, ke
      do j = js, je
       do i = is, ie
          if (cstat(i,j,k) /= 1) cycle
          wgt = wave_weight(Xh(:,2*i,2*j,2*k))
          if (wgt <= 0.0_dp) cycle
          alpha(i,j,k) = (1.0_dp - wgt)*alpha(i,j,k) + wgt*wave_fill(i, j, k, tt)
       end do
      end do
     end do
     do q = 1, 3
        select case (q)
        case (1);  Uq => U1;  uc => u
        case (2);  Uq => U2;  uc => v
        case default;  Uq => U3;  uc => w
        end select
        call face_range(q, lo, hi)
        do k = lo(3), hi(3)
         do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             pt = [i, j, k]
             if (pt(q) == 0 .or. pt(q) == ncell(cur_blk, q)) cycle     ! block boundary faces: the BCs
             if (.not. face_active(q, i, j, k)) cycle
             H = 2*pt;  H(q) = H(q) + 1
             wgt = wave_weight(Xh(:,H(1),H(2),H(3)))*alph_h(H(1),H(2),H(3))   ! water only (liquid fraction of the face)
             if (wgt <= 0.0_dp) cycle
             uw = um + wave_vel(Xh(:,H(1),H(2),H(3)), tt)
             ut = fvw(H(1),H(2),H(3))*( mq(1,q,H(1),H(2),H(3))*uw(1) + mq(2,q,H(1),H(2),H(3))*uw(2) &
                                      + mq(3,q,H(1),H(2),H(3))*uw(3) )
             Uq(i,j,k) = (1.0_dp - wgt)*Uq(i,j,k) + wgt*ut
             uc(i,j,k) = (1.0_dp - wgt)*uc(i,j,k) + wgt*uw(q)
          end do
         end do
        end do
     end do
  end do
  call refresh_all(8)
  call refresh_uvw(1)
  call contra_to_cart(.false., .true.)
  if (ovst_on) then
     call ovs_update_faces(1)
     call refresh_uvw(1)
     call contra_to_cart(.false., .true.)
  end if
end subroutine wave_relax_fields

subroutine wave_init_fields
  use global_vars
  use mpi_vars
  use waves, only: wave_on, Wave_Init, wave_vel, wave_fill, wave_mean_vel, wave_t
  use b2b_runtime, only: refresh_uvw, refresh_all
  use curvi, only: contra_to_cart, cart_to_contra
  use overset, only: ovs_update_faces
  implicit none
  integer :: ib, d, q, i, j, k, H(3)
  real(dp) :: um(3), uw(3)
  real(dp), pointer :: uc(:,:,:)
  if (.not. wave_on .or. Wave_Init == 0) return
  um = wave_mean_vel()
  do ib = 1, nown
     d = own_list(ib)
     call set_block(d)
     do k = ks, ke
      do j = js, je
       do i = is, ie
          alpha(i,j,k) = wave_fill(i, j, k, wave_t)
       end do
      end do
     end do
     ! Cartesian face velocities everywhere (ghost layers included: the
     ! transformed components are then built from complete stencils)
     do q = 1, 3
        select case (q)
        case (1);  uc => u
        case (2);  uc => v
        case default;  uc => w
        end select
        do k = lbound(uc,3), ubound(uc,3)
         do j = lbound(uc,2), ubound(uc,2)
          do i = lbound(uc,1), ubound(uc,1)
             H = [2*i, 2*j, 2*k];  H(q) = H(q) + 1
             if (any(H < hlo) .or. any(H > hhi)) cycle
             uw = um + wave_vel(Xh(:,H(1),H(2),H(3)), wave_t)
             uc(i,j,k) = uw(q)
          end do
         end do
        end do
     end do
  end do
  call refresh_all(8)
  call cart_to_contra(.false.)            ! U^q from the Cartesian faces
  do ib = 1, nown
     d = own_list(ib)
     dom(d)%us = dom(d)%u;  dom(d)%vs = dom(d)%v;  dom(d)%ws = dom(d)%w
     dom(d)%Us1 = dom(d)%U1;  dom(d)%Us2 = dom(d)%U2;  dom(d)%Us3 = dom(d)%U3
  end do
  call refresh_uvw(1)
  call contra_to_cart(.false., .true.)
  if (ovst_on) then
     call ovs_update_faces(1)
     call refresh_uvw(1)
     call contra_to_cart(.false., .true.)
     call ovs_update_faces(5)
  end if
  if (myrank == 0) write(*,'(A)') '  Wave_Init = 1: alpha and the velocity initialised from the incident wave at t = 0'
end subroutine wave_init_fields
