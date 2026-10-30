!=======================================================================
! curvi.F90
!   3-D coordinate transformation for general curvilinear grids.
!
!   Formulation (per Park, Phys. Fluids 36, 095107 (2024), "fully
!   staggered curvilinear grid", with the unweighted transformed
!   components chosen as the primary face variables):
!
!     * U^q = u_vec . grad(xi_q)  (q = 1,2,3) stored at the q-family
!       faces are the PRIMARY unknowns coupled across halos and
!       block-to-block / periodic interfaces.  A face dof is geometric,
!       so the existing exchange machinery (index maps, sign factors,
!       ownership masks) applies unchanged for any block orientation.
!     * Cartesian face velocities u, v, w (x/y/z components at the
!       i/j/k faces) are LOCAL derived fields reconstructed from U^q
!       with the forward tangents; the momentum equations are solved
!       for them (no Christoffel symbols; Ge & Sotiropoulos style).
!     * All metrics live on the HALF-INDEX lattice (CFDShip-Iowa
!       technique): cell i -> h=2i, face f -> h=2f+1, node n -> h=2n-1.
!       Metrics are exact central differences of the shared half
!       coordinates, so the discrete metric identity holds (freestream
!       preservation) and every staggered location uses one code path.
!     * Ghost node coordinates: extrapolated at physical boundaries and
!       copied from the donor block (plus the translation vector for
!       periodic pairs) at interfaces, so interface stencils see the
!       true neighbor geometry (generalizes the old 1-D ghost metrics).
!=======================================================================
module curvi
  use global_vars
  use mpi_vars
  implicit none

contains

  !=====================================================================
  ! Ghost-extended node coordinates blk(b)%xg for every block
  ! (replicated computation; deterministic on all ranks).
  !=====================================================================
  subroutine fill_ghost_coords
    integer :: b, c, i, j, k, g, m, pass, sub, wid
    integer :: na, noffn(3), a
    integer :: r(3), dnode(3), w1(3), w2(3)
    real(dp) :: tv(3)

    !---------------- interior + extrapolated ghosts --------------------
    do b = 1, nblk
       associate(bb => blk(b))
         if (allocated(bb%xg)) deallocate(bb%xg)
         allocate(bb%xg(3, 1-ng:bb%ni+ng, 1-ng:bb%nj+ng, 1-ng:bb%nk+ng))
         bb%xg = 0.0_dp
         do k = 1, bb%nk
          do j = 1, bb%nj
           do i = 1, bb%ni
              ! non-dimensional coordinates (the solver's equations
              ! are scaled by Lch, like the old 1-D metrics)
              bb%xg(1,i,j,k) = bb%xn(i,j,k) / Lch
              bb%xg(2,i,j,k) = bb%yn(i,j,k) / Lch
              bb%xg(3,i,j,k) = bb%zn(i,j,k) / Lch
           end do
          end do
         end do
         ! linear extrapolation, axis by axis (covers edges/corners too)
         do g = 1, ng
            bb%xg(:,1-g,:,:)      = 2.0_dp*bb%xg(:,2-g,:,:)      - bb%xg(:,3-g,:,:)
            bb%xg(:,bb%ni+g,:,:)  = 2.0_dp*bb%xg(:,bb%ni+g-1,:,:)- bb%xg(:,bb%ni+g-2,:,:)
         end do
         do g = 1, ng
            bb%xg(:,:,1-g,:)      = 2.0_dp*bb%xg(:,:,2-g,:)      - bb%xg(:,:,3-g,:)
            bb%xg(:,:,bb%nj+g,:)  = 2.0_dp*bb%xg(:,:,bb%nj+g-1,:)- bb%xg(:,:,bb%nj+g-2,:)
         end do
         do g = 1, ng
            bb%xg(:,:,:,1-g)      = 2.0_dp*bb%xg(:,:,:,2-g)      - bb%xg(:,:,:,3-g)
            bb%xg(:,:,:,bb%nk+g)  = 2.0_dp*bb%xg(:,:,:,bb%nk+g-1)- bb%xg(:,:,:,bb%nk+g-2)
         end do
       end associate
    end do

    !---------------- overwrite from interface donors -------------------
    ! Two passes so ghost corners can relay through the donor's own
    ! (previously filled) ghost coordinates.  Within a pass, first all
    ! patches with their tangential windows WIDENED by ng (reaches the
    ! edge/corner ghosts), then all patches with their UN-WIDENED core
    ! windows, overwriting: where several patches of one face meet the
    ! core data of the covering patch is authoritative (same rule as
    ! the field exchange in b2b_exchange.F90).
    do pass = 1, 2
      do sub = 1, 2
       wid = merge(ng, 0, sub == 1)
       do m = 1, nxp
          na = xp(m)%na
          tv = 0.0_dp
          if (m > nxp - nper) tv = pbc(m - (nxp - nper))%tv / Lch
          do a = 1, 3
             ! node-index affine map: d = sgn*r + noffn,
             ! noffn = coff - (sgn-1)/2
             noffn(a) = xp(m)%coff(a) - (xp(m)%sgn(a)-1)/2
          end do
          ! receiver node windows
          do a = 1, 3
             if (a == na) cycle
             w1(a) = max(xp(m)%rn(1,a)-wid, 1-ng)
             w2(a) = min(xp(m)%rn(2,a)+wid, nnode_of(xp(m)%rblk,a)+ng)
          end do
          if (xp(m)%rside == 1) then
             w1(na) = 1-ng;  w2(na) = 0
          else
             w1(na) = nnode_of(xp(m)%rblk,na)+1
             w2(na) = nnode_of(xp(m)%rblk,na)+ng
          end if
          do k = w1(3), w2(3)
           do j = w1(2), w2(2)
            do i = w1(1), w2(1)
               r = [i, j, k]
               do a = 1, 3
                  dnode(xp(m)%q(a)) = xp(m)%sgn(a)*r(a) + noffn(a)
               end do
               if (any(dnode < 1-ng)) cycle
               if (dnode(1) > nnode_of(xp(m)%dblk,1)+ng) cycle
               if (dnode(2) > nnode_of(xp(m)%dblk,2)+ng) cycle
               if (dnode(3) > nnode_of(xp(m)%dblk,3)+ng) cycle
               do c = 1, 3
                  blk(xp(m)%rblk)%xg(c,i,j,k) =                          &
                     blk(xp(m)%dblk)%xg(c,dnode(1),dnode(2),dnode(3)) + tv(c)
               end do
            end do
           end do
          end do
       end do
      end do
    end do

  end subroutine fill_ghost_coords

  !--------------------------------------------------------------------
  pure integer function nnode_of(b, a)
    integer, intent(in) :: b, a
    select case (a)
    case (1);  nnode_of = blk(b)%ni
    case (2);  nnode_of = blk(b)%nj
    case default;  nnode_of = blk(b)%nk
    end select
  end function nnode_of

  !=====================================================================
  ! Half-index coordinates and metrics of every OWNED sub-domain.
  !=====================================================================
  subroutine build_metrics
    integer :: ib, d, ierr, ndeg
    ndeg = 0
    do ib = 1, nown
       d = own_list(ib)
       call metrics_domain(d, ndeg)
    end do
    call MPI_Allreduce(MPI_IN_PLACE, ndeg, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    if (myrank == 0) write(*,'(A)')                                       &
       '  curvilinear     : half-index metrics built (Park m-matrix + J)'
    if (myrank == 0 .and. ndeg > 0) write(*,'(A,I0,A)') '  curvilinear     : ', ndeg, &
       ' degenerate ghost / wedge-edge metrics replaced by the nearest interior cell metric'
  end subroutine build_metrics

  ! metrics of one owned sub-domain (also called by the mesh morphing of the FSI, update 24)
  subroutine metrics_domain(d, ndeg)
    integer, intent(in) :: d
    integer, intent(inout) :: ndeg
    integer :: b, a, h1, h2, h3, n1, n2, n3, c
    real(dp) :: t(3,3), bb1(3), bb2(3), bb3(3), det
    integer  :: ierr, c1, c2, c3
       b = dom(d)%b
       do a = 1, 3
          dom(d)%hlo(a) = 2*(dom(d)%r1(a) - ng) - 1
          dom(d)%hhi(a) = 2*(dom(d)%r2(a) + ng) + 1
       end do
       associate(dd => dom(d), lo => dom(d)%hlo, hi => dom(d)%hhi)
         if (allocated(dd%Xh)) deallocate(dd%Xh, dd%Jh, dd%mh)
         allocate(dd%Xh(3, lo(1):hi(1), lo(2):hi(2), lo(3):hi(3)))
         allocate(dd%Jh(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3)))
         allocate(dd%mh(9, lo(1):hi(1), lo(2):hi(2), lo(3):hi(3)))
         dd%Jh = 1.0_dp;  dd%mh = 0.0_dp

         ! (1) odd-odd-odd = nodes (node n -> h = 2n-1)
         do h3 = lo(3), hi(3), 2
          n3 = (h3+1)/2
          do h2 = lo(2), hi(2), 2
           n2 = (h2+1)/2
           do h1 = lo(1), hi(1), 2
              n1 = (h1+1)/2
              dd%Xh(:,h1,h2,h3) = blk(b)%xg(:,n1,n2,n3)
           end do
          end do
         end do
         ! (2) fill even entries by averaging, axis by axis
         do h3 = lo(3), hi(3), 2
          do h2 = lo(2), hi(2), 2
           do h1 = lo(1)+1, hi(1)-1, 2
              dd%Xh(:,h1,h2,h3) = 0.5_dp*(dd%Xh(:,h1-1,h2,h3)+dd%Xh(:,h1+1,h2,h3))
           end do
          end do
         end do
         do h3 = lo(3), hi(3), 2
          do h2 = lo(2)+1, hi(2)-1, 2
           do h1 = lo(1), hi(1)
              dd%Xh(:,h1,h2,h3) = 0.5_dp*(dd%Xh(:,h1,h2-1,h3)+dd%Xh(:,h1,h2+1,h3))
           end do
          end do
         end do
         do h3 = lo(3)+1, hi(3)-1, 2
          do h2 = lo(2), hi(2)
           do h1 = lo(1), hi(1)
              dd%Xh(:,h1,h2,h3) = 0.5_dp*(dd%Xh(:,h1,h2,h3-1)+dd%Xh(:,h1,h2,h3+1))
           end do
          end do
         end do

         ! (3) metrics: forward tangents t_q = X(h+e_q) - X(h-e_q)
         !     (one full index unit), Jacobian, cofactors, m = b/J
         do h3 = lo(3)+1, hi(3)-1
          do h2 = lo(2)+1, hi(2)-1
           do h1 = lo(1)+1, hi(1)-1
              do c = 1, 3
                 t(c,1) = dd%Xh(c,h1+1,h2,h3) - dd%Xh(c,h1-1,h2,h3)
                 t(c,2) = dd%Xh(c,h1,h2+1,h3) - dd%Xh(c,h1,h2-1,h3)
                 t(c,3) = dd%Xh(c,h1,h2,h3+1) - dd%Xh(c,h1,h2,h3-1)
              end do
              bb1 = [ t(2,2)*t(3,3)-t(2,3)*t(3,2),                       &
                      t(3,2)*t(1,3)-t(3,3)*t(1,2),                       &
                      t(1,2)*t(2,3)-t(1,3)*t(2,2) ]
              bb2 = [ t(2,3)*t(3,1)-t(2,1)*t(3,3),                       &
                      t(3,3)*t(1,1)-t(3,1)*t(1,3),                       &
                      t(1,3)*t(2,1)-t(1,1)*t(2,3) ]
              bb3 = [ t(2,1)*t(3,2)-t(2,2)*t(3,1),                       &
                      t(3,1)*t(1,2)-t(3,2)*t(1,1),                       &
                      t(1,1)*t(2,2)-t(1,2)*t(2,1) ]
              det = t(1,1)*bb1(1) + t(1,2)*bb2(1) + t(1,3)*bb3(1)
              ! degenerate (pole) face plane: zero face Jacobian and
              ! metrics -> zero face area, no flux, no coefficient
              if (any(blk(b)%pole)) then
                 if (on_pole(b, h1, h2, h3)) then
                    dd%Jh(h1,h2,h3) = 0.0_dp
                    dd%mh(:,h1,h2,h3) = 0.0_dp
                    cycle
                 end if
              end if
              if (abs(det) < 1.0e-30_dp) then
                 ! degenerate metric.  In a GHOST position (outside the
                 ! interior half range of the domain) this happens where
                 ! the extrapolated corner-ghost coordinates collapse, e.g.
                 ! the ghost layers below a wall at the edge of a block
                 ! whose other faces are overset boundaries; those cells
                 ! are never solved - flag it (Jh = 0) and fill it from the
                 ! nearest interior position after the loop.  In the
                 ! interior it is a real grid defect: abort.
                 ! (a degenerate NODE / EDGE / FACE position on the block
                 ! boundary is legitimate as well: a wedge edge where two
                 ! block faces lie in one plane, e.g. two faces on a
                 ! symmetry plane - only a degenerate interior CELL centre
                 ! is a grid defect)
                 if (mod(h1,2) == 0 .and. mod(h2,2) == 0 .and. mod(h3,2) == 0 .and. &
                     h1 >= 2*dd%r1(1) .and. h1 <= 2*dd%r2(1) .and.            &
                     h2 >= 2*dd%r1(2) .and. h2 <= 2*dd%r2(2) .and.            &
                     h3 >= 2*dd%r1(3) .and. h3 <= 2*dd%r2(3)) then
                    write(*,'(A,I0,A,3I5)') ' ERROR (curvi): zero/negative'// &
                       ' Jacobian in block ', b, ' at half index ', h1,h2,h3
                    call MPI_Abort(MPI_COMM_WORLD, 50, ierr)
                 end if
                 dd%Jh(h1,h2,h3) = 0.0_dp;  dd%mh(:,h1,h2,h3) = 0.0_dp
                 ndeg = ndeg + 1
                 cycle
              end if
              dd%Jh(h1,h2,h3)   = det
              dd%mh(1:3,h1,h2,h3) = bb1 / det       ! m(c,1), c=1..3
              dd%mh(4:6,h1,h2,h3) = bb2 / det       ! m(c,2)
              dd%mh(7:9,h1,h2,h3) = bb3 / det       ! m(c,3)
           end do
          end do
         end do
         ! (4) degenerate ghost positions: copy the metric of the nearest
         !     interior cell position (clamped half index)
         if (ndeg > 0) then
            do h3 = lo(3)+1, hi(3)-1
             do h2 = lo(2)+1, hi(2)-1
              do h1 = lo(1)+1, hi(1)-1
                 if (dd%Jh(h1,h2,h3) /= 0.0_dp) cycle
                 if (any(blk(b)%pole)) then
                    if (on_pole(b, h1, h2, h3)) cycle
                 end if
                 c1 = min(max(h1, 2*dd%r1(1)), 2*dd%r2(1))
                 c2 = min(max(h2, 2*dd%r1(2)), 2*dd%r2(2))
                 c3 = min(max(h3, 2*dd%r1(3)), 2*dd%r2(3))
                 if (dd%Jh(c1,c2,c3) == 0.0_dp) then
                    write(*,'(A,I0,A,3I5)') ' ERROR (curvi): degenerate ghost'// &
                       ' metric without interior neighbour in block ', b, ' at ', h1,h2,h3
                    call MPI_Abort(MPI_COMM_WORLD, 50, ierr)
                 end if
                 dd%Jh(h1,h2,h3)   = dd%Jh(c1,c2,c3)
                 dd%mh(:,h1,h2,h3) = dd%mh(:,c1,c2,c3)
              end do
             end do
            end do
         end if
       end associate
  contains
    ! .true. when the half position lies on a pole plane of block bb_
    ! (min face: h = 1, max face: h = 2 n - 1 of that axis)
    pure logical function on_pole(bb_, k1, k2, k3)
      integer, intent(in) :: bb_, k1, k2, k3
      integer :: h(3), nn(3), a
      h = [k1, k2, k3]
      nn = [blk(bb_)%ni, blk(bb_)%nj, blk(bb_)%nk]
      on_pole = .false.
      do a = 1, 3
         if (blk(bb_)%pole(1,a) .and. h(a) == 1)         on_pole = .true.
         if (blk(bb_)%pole(2,a) .and. h(a) == 2*nn(a)-1) on_pole = .true.
      end do
    end function on_pole
  end subroutine metrics_domain

  !=====================================================================
  ! Cartesian -> transformed components (Park Eq. 6) at the faces this
  ! sub-domain is responsible for (owned faces + owned interface /
  ! boundary faces); star = .true. acts on us,vs,ws -> Us1..3, else on
  ! u,v,w -> U1..3.
  !=====================================================================
  subroutine cart_to_contra(star)
    logical, intent(in) :: star
    integer :: ib, d, i, j, k, f
    real(dp) :: cu, cv, cw

    do ib = 1, nown
       d = own_list(ib)
       call set_block(d)
       ! --- family 1: i-faces ------------------------------------------
       do k = ks, ke
        do j = js, je
         do f = merge(0, is, is == 1), merge(nx, ie, ie == nx)
            if (fstat(1, f, j, k) /= 1) cycle          ! overset fringe/hole face
            cu = pick(star, 1, f, j, k)
            cv = 0.25_dp*( pick(star,2,f,j-1,k) + pick(star,2,f+1,j-1,k)  &
                         + pick(star,2,f,j  ,k) + pick(star,2,f+1,j  ,k) )
            cw = 0.25_dp*( pick(star,3,f,j,k-1) + pick(star,3,f+1,j,k-1)  &
                         + pick(star,3,f,j,k  ) + pick(star,3,f+1,j,k  ) )
            call put(star, 1, f, j, k, fvw(2*f+1,2*j,2*k) *               &
               ( mq(1,1,2*f+1,2*j,2*k)*cu + mq(2,1,2*f+1,2*j,2*k)*cv      &
                                          + mq(3,1,2*f+1,2*j,2*k)*cw ))
         end do
        end do
       end do
       ! --- family 2: j-faces ------------------------------------------
       do k = ks, ke
        do f = merge(0, js, js == 1), merge(ny, je, je == ny)
         do i = is, ie
            if (fstat(2, i, f, k) /= 1) cycle
            cv = pick(star, 2, i, f, k)
            cu = 0.25_dp*( pick(star,1,i-1,f,k) + pick(star,1,i-1,f+1,k)  &
                         + pick(star,1,i  ,f,k) + pick(star,1,i  ,f+1,k) )
            cw = 0.25_dp*( pick(star,3,i,f,k-1) + pick(star,3,i,f+1,k-1)  &
                         + pick(star,3,i,f,k  ) + pick(star,3,i,f+1,k  ) )
            call put(star, 2, i, f, k, fvw(2*i,2*f+1,2*k) *               &
               ( mq(1,2,2*i,2*f+1,2*k)*cu + mq(2,2,2*i,2*f+1,2*k)*cv      &
                                          + mq(3,2,2*i,2*f+1,2*k)*cw ))
         end do
        end do
       end do
       ! --- family 3: k-faces ------------------------------------------
       do f = merge(0, ks, ks == 1), merge(nz, ke, ke == nz)
        do j = js, je
         do i = is, ie
            if (fstat(3, i, j, f) /= 1) cycle
            cw = pick(star, 3, i, j, f)
            cu = 0.25_dp*( pick(star,1,i-1,j,f) + pick(star,1,i-1,j,f+1)  &
                         + pick(star,1,i  ,j,f) + pick(star,1,i  ,j,f+1) )
            cv = 0.25_dp*( pick(star,2,i,j-1,f) + pick(star,2,i,j-1,f+1)  &
                         + pick(star,2,i,j  ,f) + pick(star,2,i,j  ,f+1) )
            call put(star, 3, i, j, f, fvw(2*i,2*j,2*f+1) *               &
               ( mq(1,3,2*i,2*j,2*f+1)*cu + mq(2,3,2*i,2*j,2*f+1)*cv      &
                                          + mq(3,3,2*i,2*j,2*f+1)*cw ))
         end do
        end do
       end do
    end do

  contains
    pure real(dp) function pick(st, c, a1, a2, a3)
      logical, intent(in) :: st
      integer, intent(in) :: c, a1, a2, a3
      if (st) then
         select case (c)
         case (1);  pick = us(a1,a2,a3)
         case (2);  pick = vs(a1,a2,a3)
         case default;  pick = ws(a1,a2,a3)
         end select
      else
         select case (c)
         case (1);  pick = u(a1,a2,a3)
         case (2);  pick = v(a1,a2,a3)
         case default;  pick = w(a1,a2,a3)
         end select
      end if
    end function pick
    subroutine put(st, q, a1, a2, a3, val)
      logical, intent(in) :: st
      integer, intent(in) :: q, a1, a2, a3
      real(dp), intent(in) :: val
      if (st) then
         select case (q)
         case (1);  Us1(a1,a2,a3) = val
         case (2);  Us2(a1,a2,a3) = val
         case default;  Us3(a1,a2,a3) = val
         end select
      else
         select case (q)
         case (1);  U1(a1,a2,a3) = val
         case (2);  U2(a1,a2,a3) = val
         case default;  U3(a1,a2,a3) = val
         end select
      end if
    end subroutine put
  end subroutine cart_to_contra

  !=====================================================================
  ! Transformed components -> Cartesian face velocities over the full
  ! ghost-extended range (u_c = sum_q t_q(c) U^q with the forward
  ! tangents; interpolation of the two foreign families to the face).
  !=====================================================================
  subroutine contra_to_cart(star, ghosts_only)
    logical, intent(in) :: star
    logical, intent(in) :: ghosts_only   ! .true.: do not touch the
                                         ! positions the momentum sweep
                                         ! iterates (interior solution)
    integer :: ib, d, i, j, k, f
    integer :: f1, f2, t1a, t1b, t2a, t2b, t3a, t3b
    real(dp) :: q1, q2, q3

    do ib = 1, nown
       d = own_list(ib)
       call set_block(d)
       f1 = dom(d)%r1(1)-ng;  f2 = dom(d)%r2(1)+ng-1
       t2a = dom(d)%r1(2)-ng+1;  t2b = dom(d)%r2(2)+ng-1
       t3a = dom(d)%r1(3)-ng+1;  t3b = dom(d)%r2(3)+ng-1
       do k = t3a, t3b
        do j = t2a, t2b
         do f = f1, f2
            ! keep the momentum iterate: skip exactly the sweep set
            ! (interior owned faces + OWNED interface faces; BC-set
            ! boundary faces are reconstructed from the transformed BC)
            if (ghosts_only .and. j >= js .and. j <= je .and.             &
                k >= ks .and. k <= ke .and. f >= merge(0,is,is==1) .and.  &
                f <= merge(nx,ie,ie==nx)) then
               if (f == 0) then
                  if (mu_lo(j,k)) cycle
               else if (f == nx) then
                  if (mu_hi(j,k)) cycle
               else
                  cycle
               end if
            end if
            q1 = gv(star,1,f,j,k)
            q2 = 0.25_dp*( gv(star,2,f,j-1,k)+gv(star,2,f+1,j-1,k)        &
                         + gv(star,2,f,j  ,k)+gv(star,2,f+1,j  ,k) )
            q3 = 0.25_dp*( gv(star,3,f,j,k-1)+gv(star,3,f+1,j,k-1)        &
                         + gv(star,3,f,j,k  )+gv(star,3,f+1,j,k  ) )
            call pv(star, 1, f, j, k,                                     &
                 tq(1,1,2*f+1,2*j,2*k)*q1 + tq(1,2,2*f+1,2*j,2*k)*q2      &
                                          + tq(1,3,2*f+1,2*j,2*k)*q3)
         end do
        end do
       end do
       f1 = dom(d)%r1(2)-ng;  f2 = dom(d)%r2(2)+ng-1
       t1a = dom(d)%r1(1)-ng+1;  t1b = dom(d)%r2(1)+ng-1
       do k = t3a, t3b
        do f = f1, f2
         do i = t1a, t1b
            if (ghosts_only .and. i >= is .and. i <= ie .and.             &
                k >= ks .and. k <= ke .and. f >= merge(0,js,js==1) .and.  &
                f <= merge(ny,je,je==ny)) then
               if (f == 0) then
                  if (mv_lo(i,k)) cycle
               else if (f == ny) then
                  if (mv_hi(i,k)) cycle
               else
                  cycle
               end if
            end if
            q2 = gv(star,2,i,f,k)
            q1 = 0.25_dp*( gv(star,1,i-1,f,k)+gv(star,1,i-1,f+1,k)        &
                         + gv(star,1,i  ,f,k)+gv(star,1,i  ,f+1,k) )
            q3 = 0.25_dp*( gv(star,3,i,f,k-1)+gv(star,3,i,f+1,k-1)        &
                         + gv(star,3,i,f,k  )+gv(star,3,i,f+1,k  ) )
            call pv(star, 2, i, f, k,                                     &
                 tq(2,1,2*i,2*f+1,2*k)*q1 + tq(2,2,2*i,2*f+1,2*k)*q2      &
                                          + tq(2,3,2*i,2*f+1,2*k)*q3)
         end do
        end do
       end do
       f1 = dom(d)%r1(3)-ng;  f2 = dom(d)%r2(3)+ng-1
       do f = f1, f2
        do j = t2a, t2b
         do i = t1a, t1b
            if (ghosts_only .and. i >= is .and. i <= ie .and.             &
                j >= js .and. j <= je .and. f >= merge(0,ks,ks==1) .and.  &
                f <= merge(nz,ke,ke==nz)) then
               if (f == 0) then
                  if (mw_lo(i,j)) cycle
               else if (f == nz) then
                  if (mw_hi(i,j)) cycle
               else
                  cycle
               end if
            end if
            q3 = gv(star,3,i,j,f)
            q1 = 0.25_dp*( gv(star,1,i-1,j,f)+gv(star,1,i-1,j,f+1)        &
                         + gv(star,1,i  ,j,f)+gv(star,1,i  ,j,f+1) )
            q2 = 0.25_dp*( gv(star,2,i,j-1,f)+gv(star,2,i,j-1,f+1)        &
                         + gv(star,2,i,j  ,f)+gv(star,2,i,j  ,f+1) )
            call pv(star, 3, i, j, f,                                     &
                 tq(3,1,2*i,2*j,2*f+1)*q1 + tq(3,2,2*i,2*j,2*f+1)*q2      &
                                          + tq(3,3,2*i,2*j,2*f+1)*q3)
         end do
        end do
       end do
    end do
    call dirichlet_cart(star)

  contains
    pure real(dp) function gv(st, q, a1, a2, a3)
      logical, intent(in) :: st
      integer, intent(in) :: q, a1, a2, a3
      integer :: hx, hy, hz
      if (st) then
         select case (q)
         case (1);  gv = Us1(a1,a2,a3)
         case (2);  gv = Us2(a1,a2,a3)
         case default;  gv = Us3(a1,a2,a3)
         end select
      else
         select case (q)
         case (1);  gv = U1(a1,a2,a3)
         case (2);  gv = U2(a1,a2,a3)
         case default;  gv = U3(a1,a2,a3)
         end select
      end if
      ! pointwise conversion to the unweighted component at the value's
      ! own face position (no-op for curvi_form = 1)
      if (curvi_form == 2) then
         hx = 2*a1;  hy = 2*a2;  hz = 2*a3
         select case (q)
         case (1);  hx = hx + 1
         case (2);  hy = hy + 1
         case default;  hz = hz + 1
         end select
         if (Jh(hx, hy, hz) /= 0.0_dp) then
            gv = gv / Jh(hx, hy, hz)
         else
            gv = 0.0_dp                     ! pole face (zero area)
         end if
      end if
    end function gv
    subroutine pv(st, c, a1, a2, a3, val)
      logical, intent(in) :: st
      integer, intent(in) :: c, a1, a2, a3
      real(dp), intent(in) :: val
      if (st) then
         select case (c)
         case (1);  us(a1,a2,a3) = val
         case (2);  vs(a1,a2,a3) = val
         case default;  ws(a1,a2,a3) = val
         end select
      else
         select case (c)
         case (1);  u(a1,a2,a3) = val
         case (2);  v(a1,a2,a3) = val
         case default;  w(a1,a2,a3) = val
         end select
      end if
    end subroutine pv
  end subroutine contra_to_cart

  !=====================================================================
  ! Dirichlet velocity patches (.bcs flag 1 / 11): the CARTESIAN face
  ! values on the boundary face (normal family) and in the ghost layers
  ! (mirror 2 u_bc - u_interior, own component of every family) are the
  ! prescribed wall / inlet velocity - not the reconstruction t_q U^q of
  ! contra_to_cart.  At a wall edge where two grid families are strongly
  ! non-orthogonal (the tube grids of Box_6DoF wrapping a box edge, or a
  ! degenerate ghost corner) the round trip transformed -> Cartesian ->
  ! mirror -> transformed of the ghost/boundary faces has a gain > 1 and
  ! the momentum sweeps amplified the boundary velocities by ~20 per
  ! iteration.  With the known Cartesian values there is nothing to
  ! reconstruct and the loop is broken.  (update 12)
  !=====================================================================
  subroutine dirichlet_cart(star)
    use waves, only: wave_vel, wave_mean_vel, wave_t
    logical, intent(in) :: star
    integer :: ib, d, n, c, dir, ax, sd, flag, i, j, k
    integer :: r1(3), r2(3), nn(3), idx(3), g(3), mm, hh(3), gh(3), lo3(3), hi3(3)
    real(dp) :: ub(3), ubc(3), xw(3)
    real(dp), pointer :: fc(:,:,:)
    logical :: wave
    do ib = 1, nown
       d = own_list(ib)
       call set_block(d)
       nn = [nx, ny, nz]
       do n = 1, nbc
          dir = bc_dir(n);  ax = abs(dir);  sd = merge(1, 2, dir > 0)
          lo3 = [is, js, ks];  hi3 = [ie, je, ke]
          if (sd == 1 .and. lo3(ax) /= 1) cycle
          if (sd == 2 .and. hi3(ax) /= nn(ax)) cycle
          do c = 1, 3
             flag = bc_flag(c, n)
             if (flag == 11) flag = 1
             wave = (flag == 12)                        ! incident-wave patch (update 18): mean + wave at the position
             if (wave) flag = 1
             if (flag /= 1) cycle
             ubc(1) = u_init/Uch + bc_ramp*(bc_val(1,n) - u_init/Uch)
             ubc(2) = v_init/Uch + bc_ramp*(bc_val(2,n) - v_init/Uch)
             ubc(3) = w_init/Uch + bc_ramp*(bc_val(3,n) - w_init/Uch)
             if (wave) ubc = wave_mean_vel()
             if (star) then
                select case (c)
                case (1);  fc => us
                case (2);  fc => vs
                case default;  fc => ws
                end select
             else
                select case (c)
                case (1);  fc => u
                case (2);  fc => v
                case default;  fc => w
                end select
             end if
             ! tangential index ranges of family c on this patch (as apply_bc)
             call trange(1, c, bc_irng(1,n), bc_irng(2,n), ni, is, ie, r1(1), r2(1))
             call trange(2, c, bc_jrng(1,n), bc_jrng(2,n), nj, js, je, r1(2), r2(2))
             call trange(3, c, bc_krng(1,n), bc_krng(2,n), nk, ks, ke, r1(3), r2(3))
             do k = r1(3), r2(3)
              do j = r1(2), r2(2)
               do i = r1(1), r2(1)
                  idx = [i, j, k]
                  if (idx(ax) /= merge(0, nn(ax), sd == 1) .and. c == ax) cycle   ! loop hits every layer: act once
                  if (idx(ax) /= merge(1, nn(ax), sd == 1) .and. c /= ax) cycle
                  if (c == ax) then
                     ! normal family: boundary face + two mirrored ghost faces
                     g = idx                                            ! boundary face position
                     hh = 2*g;  hh(c) = hh(c) + 1
                     xw = Xh(:, hh(1), hh(2), hh(3))
                     ub = ubc
                     if (wave) ub = ub + wave_vel(xw, wave_t)
                     if (any_moving .and. associated(wall_velocity_fn)) ub = ub + wall_velocity_fn(cur_blk, xw, hh)
                     fc(g(1),g(2),g(3)) = ub(c)
                     do mm = 1, 2
                        g = idx;  g(ax) = g(ax) + merge(-mm, mm, sd == 1)      ! ghost face
                        gh = idx;  gh(ax) = gh(ax) + merge(mm, -mm, sd == 1)   ! mirror partner
                        if (g(ax) < lbound(fc,ax) .or. g(ax) > ubound(fc,ax)) cycle
                        fc(g(1),g(2),g(3)) = 2.0_dp*ub(c) - fc(gh(1),gh(2),gh(3))
                     end do
                  else
                     ! tangential family: ghost cell levels mirrored about the wall
                     do mm = 1, 2
                        g = idx;  g(ax) = merge(1 - mm, nn(ax) + mm, sd == 1)          ! ghost level
                        gh = idx; gh(ax) = merge(mm, nn(ax) + 1 - mm, sd == 1)         ! interior partner
                        if (g(ax) < lbound(fc,ax) .or. g(ax) > ubound(fc,ax)) cycle
                        hh = 2*g;  hh(c) = hh(c) + 1
                        xw = 0.5_dp*(Xh(:, hh(1), hh(2), hh(3)) + Xh(:, 2*gh(1)+merge(1,0,c==1), 2*gh(2)+merge(1,0,c==2), 2*gh(3)+merge(1,0,c==3)))
                        ub = ubc
                        if (wave) ub = ub + wave_vel(xw, wave_t)
                        if (any_moving .and. associated(wall_velocity_fn)) ub = ub + wall_velocity_fn(cur_blk, xw, &
                           (hh + [2*gh(1)+merge(1,0,c==1), 2*gh(2)+merge(1,0,c==2), 2*gh(3)+merge(1,0,c==3)])/2)
                        fc(g(1),g(2),g(3)) = 2.0_dp*ub(c) - fc(gh(1),gh(2),gh(3))
                     end do
                  end if
               end do
              end do
             end do
          end do
       end do
    end do
  contains
    subroutine trange(dr, ivar_, n1, n2, nnode, ls, le, ra, rb)
      integer, intent(in)  :: dr, ivar_, n1, n2, nnode, ls, le
      integer, intent(out) :: ra, rb
      if (ivar_ == dr) then
         ra = n1 - 1;  rb = n2 - 1
      else
         ra = n1;  rb = n2 - 1
      end if
      if (n1 <= 1)     ra = ra - ng
      if (n2 >= nnode) rb = rb + ng
      ra = max(ra, ls - ng);  rb = min(rb, le + ng)
    end subroutine trange
  end subroutine dirichlet_cart

end module curvi
