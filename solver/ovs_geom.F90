!=======================================================================
! ovs_geom.F90
!   Geometry kernels of the overset (Chimera) grid assembly:
!
!   (1) cell-search buckets of one grid block + inverse-trilinear point
!       location (which hexahedral cell contains a point, and where in
!       its local (r,s,t) coordinates);
!   (2) triangle surfaces (the no-slip wall patches and the physical
!       domain-boundary patches) with a bucket search, closest-point
!       queries and the angle-weighted pseudo-normal inside/outside test
!       (Baerentzen & Aanaes 2005) that is robust for CLOSED surfaces
!       made of OVERLAPPING patches (hydrofoil / KCS type grids) and for
!       sharp edges (foil trailing edge);
!   (3) cell volumes in the solver's own definition (determinant of the
!       three face-centre difference vectors = J at the cell centre).
!=======================================================================
module ovs_geom
  implicit none
  private
  public :: bkt_t, build_buckets, locate_point, cell_volumes, inv_trilinear
  public :: surf_wind_prep, surf_wind
  integer(8), public :: wind_nq = 0, wind_ntri = 0, wind_ncl = 0   ! statistics: queries, exact triangles, clusters visited
  public :: plic_cube_frac, plic_plane_const, plic_box_frac, plic_normal27
  public :: surf_t, surf_init, surf_add_patch, surf_finish, surf_closest,   &
            surf_inside, surf_free, surf_free_all, tri_closest

  integer, parameter :: dp = selected_real_kind(15, 307)

  !------------------------------------------------ cell search buckets --
  type :: bkt_t
     integer  :: ni = 0, nj = 0, nk = 0      ! node counts of the block
     integer  :: nx = 0, ny = 0, nz = 0      ! cell counts
     integer  :: n(3) = 1                    ! buckets per axis
     real(dp) :: o(3) = 0.0_dp, h(3) = 1.0_dp, bmin(3), bmax(3)
     integer, allocatable :: start(:)        ! CSR bucket -> cell list
     integer, allocatable :: list(:)
  end type bkt_t

  !------------------------------------------------ triangle surface -----
  ! one level of the cluster hierarchy of the winding number: the fine
  ! buckets coarsened 2^l times per axis; per cluster and body the
  ! area-weighted dipole, the area centroid, the area, and the radius
  ! (bounding sphere of the vertices about the centroid)
  type :: wlev_t
     integer :: n(3) = 1
     real(dp), allocatable :: d(:,:), c(:,:), a(:), r(:)   ! dipole (3,ncl), centroid, area, radius
     real(dp), allocatable :: q(:,:,:)           ! (3,3,ncl) second moment  sum_t a_t n_t (c_t - c)^T
  end type wlev_t
  ! cluster tree of one body: base grid over its bounding box + levels
  type :: wtree_t
     integer  :: nt = 0, n(3) = 1, nlev = 0
     real(dp) :: o(3) = 0.0_dp, h(3) = 1.0_dp
     real(dp), allocatable :: fd(:,:), fc(:,:), fa(:), fr(:), fq(:,:,:)   ! base cells
     integer,  allocatable :: cstart(:), clist(:)   ! CSR: triangles of each base cell (each once)
     type(wlev_t), allocatable :: lev(:)            ! (nlev) levels, 1 = coarsest
  end type wtree_t

  type :: surf_t
     integer :: nv = 0, nt = 0
     real(dp), allocatable :: v(:,:)         ! (3, nv) vertices
     real(dp), allocatable :: vn(:,:)        ! (3, nv) pseudo normals (angle weighted)
     integer,  allocatable :: t(:,:)         ! (3, nt) vertex ids
     real(dp), allocatable :: tn(:,:)        ! (3, nt) unit normals (into the fluid)
     integer,  allocatable :: tp(:)          ! (nt) patch id of the triangle
     integer,  allocatable :: tb(:)          ! (nt) body id (set by the caller)
     integer  :: n(3) = 1                    ! buckets
     real(dp) :: o(3), h(3), bmin(3), bmax(3)
     integer, allocatable :: start(:), list(:)
     logical :: ready = .false.
     !--- generalized winding number (surf_wind_prep / surf_wind) -------
     logical,  allocatable :: keep(:)        ! (nt) .false. = trimmed (covered by an overlapping patch)
     real(dp), allocatable :: ta(:), tcen(:,:)  ! (nt) area, (3,nt) centroid
     integer  :: nbody = 0
     type(wtree_t), allocatable :: wt(:)         ! (nbody) cluster tree of every body
     integer,  allocatable :: stamp(:), tbk(:)   ! query stamp, base cell (of its body's tree) of each triangle
     integer  :: qid = 0
     logical  :: wind_ready = .false.
  end type surf_t

contains

  !=====================================================================
  ! (1) BUCKETS + POINT LOCATION
  !=====================================================================
  subroutine build_buckets(x, y, z, ni, nj, nk, bk)
    integer,  intent(in)  :: ni, nj, nk
    real(dp), intent(in)  :: x(ni,nj,nk), y(ni,nj,nk), z(ni,nj,nk)
    type(bkt_t), intent(out) :: bk
    integer  :: nc, a, i, j, k, c, b1, b2, b3, l1(3), l2(3), nb, m
    real(dp) :: L(3), Lg, cmin(3), cmax(3)
    integer, allocatable :: cnt(:)

    bk%ni = ni;  bk%nj = nj;  bk%nk = nk
    bk%nx = ni-1;  bk%ny = nj-1;  bk%nz = nk-1
    nc = bk%nx*bk%ny*bk%nz
    bk%bmin = [minval(x), minval(y), minval(z)]
    bk%bmax = [maxval(x), maxval(y), maxval(z)]
    L = bk%bmax - bk%bmin
    do a = 1, 3
       if (L(a) <= 0.0_dp) L(a) = 1.0e-12_dp
    end do
    Lg = (L(1)*L(2)*L(3))**(1.0_dp/3.0_dp)
    do a = 1, 3
       bk%n(a) = max(1, min(nint(real(nc,dp)**(1.0_dp/3.0_dp)*L(a)/Lg), 4096))
    end do
    ! keep the bucket count near the cell count
    do while (product(bk%n) > 2*nc + 8)
       a = maxloc(bk%n, 1)
       bk%n(a) = max(1, bk%n(a)/2)
    end do
    bk%o = bk%bmin - 1.0e-9_dp*L
    bk%h = (L + 2.0e-9_dp*L) / real(bk%n, dp)
    nb = product(bk%n)
    allocate(cnt(nb), bk%start(nb+1))
    cnt = 0
    ! pass 1: count
    do k = 1, bk%nz
     do j = 1, bk%ny
      do i = 1, bk%nx
         call cell_box(x, y, z, ni, nj, nk, i, j, k, cmin, cmax)
         call box_range(bk, cmin, cmax, l1, l2)
         do b3 = l1(3), l2(3)
          do b2 = l1(2), l2(2)
           do b1 = l1(1), l2(1)
              m = bid(bk, b1, b2, b3)
              cnt(m) = cnt(m) + 1
           end do
          end do
         end do
      end do
     end do
    end do
    bk%start(1) = 1
    do m = 1, nb
       bk%start(m+1) = bk%start(m) + cnt(m)
    end do
    allocate(bk%list(bk%start(nb+1)-1))
    cnt = 0
    ! pass 2: fill
    do k = 1, bk%nz
     do j = 1, bk%ny
      do i = 1, bk%nx
         c = i + bk%nx*((j-1) + bk%ny*(k-1))
         call cell_box(x, y, z, ni, nj, nk, i, j, k, cmin, cmax)
         call box_range(bk, cmin, cmax, l1, l2)
         do b3 = l1(3), l2(3)
          do b2 = l1(2), l2(2)
           do b1 = l1(1), l2(1)
              m = bid(bk, b1, b2, b3)
              bk%list(bk%start(m) + cnt(m)) = c
              cnt(m) = cnt(m) + 1
           end do
          end do
         end do
      end do
     end do
    end do
    deallocate(cnt)
  end subroutine build_buckets

  pure integer function bid(bk, b1, b2, b3)
    type(bkt_t), intent(in) :: bk
    integer, intent(in) :: b1, b2, b3
    bid = b1 + bk%n(1)*((b2-1) + bk%n(2)*(b3-1))
  end function bid

  pure subroutine box_range(bk, cmin, cmax, l1, l2)
    type(bkt_t), intent(in) :: bk
    real(dp), intent(in) :: cmin(3), cmax(3)
    integer, intent(out) :: l1(3), l2(3)
    integer :: a
    do a = 1, 3
       l1(a) = max(1, min(bk%n(a), int(floor((cmin(a) - bk%o(a))/bk%h(a))) + 1))
       l2(a) = max(1, min(bk%n(a), int(floor((cmax(a) - bk%o(a))/bk%h(a))) + 1))
    end do
  end subroutine box_range

  pure subroutine cell_box(x, y, z, ni, nj, nk, i, j, k, cmin, cmax)
    integer,  intent(in) :: ni, nj, nk, i, j, k
    real(dp), intent(in) :: x(ni,nj,nk), y(ni,nj,nk), z(ni,nj,nk)
    real(dp), intent(out) :: cmin(3), cmax(3)
    cmin(1) = minval(x(i:i+1,j:j+1,k:k+1));  cmax(1) = maxval(x(i:i+1,j:j+1,k:k+1))
    cmin(2) = minval(y(i:i+1,j:j+1,k:k+1));  cmax(2) = maxval(y(i:i+1,j:j+1,k:k+1))
    cmin(3) = minval(z(i:i+1,j:j+1,k:k+1));  cmax(3) = maxval(z(i:i+1,j:j+1,k:k+1))
  end subroutine cell_box

  !---------------------------------------------------------------------
  ! Locate point p in the block: returns the containing cell (ic,jc,kc)
  ! and local coordinates rst in [0,1]^3 (trilinear), found = .false.
  ! when no cell contains it.  If several candidates contain the point
  ! (shared faces), the one where the point is most interior is taken -
  ! symmetric under mirror images of the grid.
  !---------------------------------------------------------------------
  ! xi / xd (optional): exclude the cells within index distance xd of
  ! cell xi - the search of a block in ITSELF (self-overlap: the seam of
  ! a C- or O-type grid whose two ends overlap without a point match)
  subroutine locate_point(bk, x, y, z, p, ic, jc, kc, rst, found, xi, xd)
    type(bkt_t), intent(in) :: bk
    real(dp), intent(in) :: x(bk%ni,bk%nj,bk%nk), y(bk%ni,bk%nj,bk%nk), z(bk%ni,bk%nj,bk%nk)
    real(dp), intent(in) :: p(3)
    integer,  intent(out) :: ic, jc, kc
    real(dp), intent(out) :: rst(3)
    logical,  intent(out) :: found
    integer,  intent(in), optional :: xi(3), xd
    integer  :: l1(3), l2(3), m, n, c, i, j, k
    real(dp) :: cmin(3), cmax(3), r(3), q, qbest, eps, tolr, cc(3), cb(3), key(3), keyb(3), scale
    logical  :: ok

    found = .false.;  ic = 0;  jc = 0;  kc = 0;  rst = 0.0_dp
    qbest = huge(1.0_dp)
    if (any(p < bk%bmin - 1.0e-9_dp*(bk%bmax-bk%bmin)) .or.                &
        any(p > bk%bmax + 1.0e-9_dp*(bk%bmax-bk%bmin))) return
    call box_range(bk, p, p, l1, l2)
    m = bid(bk, l1(1), l1(2), l1(3))
    tolr = 1.0e-7_dp
    do n = bk%start(m), bk%start(m+1)-1
       c = bk%list(n)
       k = (c-1)/(bk%nx*bk%ny) + 1
       j = (c-1 - (k-1)*bk%nx*bk%ny)/bk%nx + 1
       i = c - bk%nx*((j-1) + bk%ny*(k-1))
       if (present(xi)) then
          if (maxval(abs([i, j, k] - xi)) <= xd) cycle
       end if
       call cell_box(x, y, z, bk%ni, bk%nj, bk%nk, i, j, k, cmin, cmax)
       eps = 1.0e-9_dp*maxval(cmax - cmin)
       if (any(p < cmin - eps) .or. any(p > cmax + eps)) cycle
       call inv_trilinear(x, y, z, bk%ni, bk%nj, bk%nk, i, j, k, p, r, ok)
       if (.not. ok) cycle
       if (any(r < -tolr) .or. any(r > 1.0_dp + tolr)) cycle
       q = maxval(abs(r - 0.5_dp))
       if (q < qbest - 1.0e-9_dp) then
          qbest = q;  found = .true.
          ic = i;  jc = j;  kc = k
          rst = min(max(r, 0.0_dp), 1.0_dp)
       else if (found .and. abs(q - qbest) <= 1.0e-9_dp) then
          ! tie (the point lies on a face / edge shared by two cells): a
          ! MIRROR-INVARIANT canonical choice - the cell whose centre has
          ! the larger x, then the larger |y|, then the larger z - so that
          ! mirrored grids get mirrored donors independently of the
          ! bucket / index order
          cc = cell_ctr(i, j, k);  cb = cell_ctr(ic, jc, kc)
          key = [cc(1), abs(cc(2)), cc(3)];  keyb = [cb(1), abs(cb(2)), cb(3)]
          scale = max(maxval(abs(key)), maxval(abs(keyb)), 1.0e-30_dp)
          if (key(1) > keyb(1) + 1.0e-12_dp*scale .or.                                       &
              (abs(key(1)-keyb(1)) <= 1.0e-12_dp*scale .and. (key(2) > keyb(2) + 1.0e-12_dp*scale .or. &
              (abs(key(2)-keyb(2)) <= 1.0e-12_dp*scale .and. key(3) > keyb(3) + 1.0e-12_dp*scale)))) then
             qbest = q;  ic = i;  jc = j;  kc = k
             rst = min(max(r, 0.0_dp), 1.0_dp)
          end if
       end if
    end do
  contains
    function cell_ctr(a, b, c) result(cc)
      integer, intent(in) :: a, b, c
      real(dp) :: cc(3)
      cc(1) = 0.125_dp*sum(x(a:a+1,b:b+1,c:c+1))
      cc(2) = 0.125_dp*sum(y(a:a+1,b:b+1,c:c+1))
      cc(3) = 0.125_dp*sum(z(a:a+1,b:b+1,c:c+1))
    end function cell_ctr
  end subroutine locate_point

  !---------------------------------------------------------------------
  ! Newton inversion of the trilinear map of cell (i,j,k):
  !   x(r,s,t) = sum_n N_n(r,s,t) X_n,  N: standard trilinear shapes.
  !---------------------------------------------------------------------
  pure subroutine inv_trilinear(x, y, z, ni, nj, nk, i, j, k, p, r, ok)
    integer,  intent(in)  :: ni, nj, nk, i, j, k
    real(dp), intent(in)  :: x(ni,nj,nk), y(ni,nj,nk), z(ni,nj,nk), p(3)
    real(dp), intent(out) :: r(3)
    logical,  intent(out) :: ok
    real(dp) :: X8(3,8), f(3), Jm(3,3), d(3), det, Nsh(8), dN(8,3)
    integer  :: it, n, a, b
    real(dp) :: rr, ss, tt, sc

    X8(:,1) = [x(i,j,k),       y(i,j,k),       z(i,j,k)]
    X8(:,2) = [x(i+1,j,k),     y(i+1,j,k),     z(i+1,j,k)]
    X8(:,3) = [x(i,j+1,k),     y(i,j+1,k),     z(i,j+1,k)]
    X8(:,4) = [x(i+1,j+1,k),   y(i+1,j+1,k),   z(i+1,j+1,k)]
    X8(:,5) = [x(i,j,k+1),     y(i,j,k+1),     z(i,j,k+1)]
    X8(:,6) = [x(i+1,j,k+1),   y(i+1,j,k+1),   z(i+1,j,k+1)]
    X8(:,7) = [x(i,j+1,k+1),   y(i,j+1,k+1),   z(i,j+1,k+1)]
    X8(:,8) = [x(i+1,j+1,k+1), y(i+1,j+1,k+1), z(i+1,j+1,k+1)]
    sc = 0.0_dp
    do n = 2, 8
       sc = max(sc, maxval(abs(X8(:,n) - X8(:,1))))
    end do
    if (sc <= 0.0_dp) sc = 1.0_dp
    r = 0.5_dp;  ok = .false.
    do it = 1, 30
       rr = r(1);  ss = r(2);  tt = r(3)
       Nsh(1) = (1-rr)*(1-ss)*(1-tt);  Nsh(2) = rr*(1-ss)*(1-tt)
       Nsh(3) = (1-rr)*ss*(1-tt);      Nsh(4) = rr*ss*(1-tt)
       Nsh(5) = (1-rr)*(1-ss)*tt;      Nsh(6) = rr*(1-ss)*tt
       Nsh(7) = (1-rr)*ss*tt;          Nsh(8) = rr*ss*tt
       dN(1,:) = [-(1-ss)*(1-tt), -(1-rr)*(1-tt), -(1-rr)*(1-ss)]
       dN(2,:) = [ (1-ss)*(1-tt), -rr*(1-tt),     -rr*(1-ss)]
       dN(3,:) = [-ss*(1-tt),      (1-rr)*(1-tt), -(1-rr)*ss]
       dN(4,:) = [ ss*(1-tt),      rr*(1-tt),     -rr*ss]
       dN(5,:) = [-(1-ss)*tt,     -(1-rr)*tt,      (1-rr)*(1-ss)]
       dN(6,:) = [ (1-ss)*tt,     -rr*tt,          rr*(1-ss)]
       dN(7,:) = [-ss*tt,          (1-rr)*tt,      (1-rr)*ss]
       dN(8,:) = [ ss*tt,          rr*tt,          rr*ss]
       f = 0.0_dp;  Jm = 0.0_dp
       do n = 1, 8
          f = f + Nsh(n)*X8(:,n)
          do b = 1, 3
             do a = 1, 3
                Jm(a,b) = Jm(a,b) + dN(n,b)*X8(a,n)
             end do
          end do
       end do
       f = p - f
       det = Jm(1,1)*(Jm(2,2)*Jm(3,3)-Jm(2,3)*Jm(3,2))                    &
           - Jm(1,2)*(Jm(2,1)*Jm(3,3)-Jm(2,3)*Jm(3,1))                    &
           + Jm(1,3)*(Jm(2,1)*Jm(3,2)-Jm(2,2)*Jm(3,1))
       if (abs(det) < 1.0e-300_dp) return
       d(1) = ( f(1)*(Jm(2,2)*Jm(3,3)-Jm(2,3)*Jm(3,2))                    &
              - Jm(1,2)*(f(2)*Jm(3,3)-Jm(2,3)*f(3))                       &
              + Jm(1,3)*(f(2)*Jm(3,2)-Jm(2,2)*f(3)) ) / det
       d(2) = ( Jm(1,1)*(f(2)*Jm(3,3)-Jm(2,3)*f(3))                       &
              - f(1)*(Jm(2,1)*Jm(3,3)-Jm(2,3)*Jm(3,1))                    &
              + Jm(1,3)*(Jm(2,1)*f(3)-f(2)*Jm(3,1)) ) / det
       d(3) = ( Jm(1,1)*(Jm(2,2)*f(3)-f(2)*Jm(3,2))                       &
              - Jm(1,2)*(Jm(2,1)*f(3)-f(2)*Jm(3,1))                       &
              + f(1)*(Jm(2,1)*Jm(3,2)-Jm(2,2)*Jm(3,1)) ) / det
       r = r + d
       if (any(abs(r) > 50.0_dp)) return
       if (maxval(abs(d)) < 1.0e-11_dp .or. maxval(abs(f)) < 1.0e-13_dp*sc) then
          ok = .true.;  return
       end if
    end do
    ok = (maxval(abs(d)) < 1.0e-8_dp)
  end subroutine inv_trilinear

  !---------------------------------------------------------------------
  ! Cell volumes: det[X(i+) - X(i-), X(j+) - X(j-), X(k+) - X(k-)] with
  ! the face centres X(+/-) (the solver's J at the cell centre).
  !---------------------------------------------------------------------
  subroutine cell_volumes(x, y, z, ni, nj, nk, vol)
    integer,  intent(in)  :: ni, nj, nk
    real(dp), intent(in)  :: x(ni,nj,nk), y(ni,nj,nk), z(ni,nj,nk)
    real(dp), intent(out) :: vol(ni-1,nj-1,nk-1)
    integer  :: i, j, k
    real(dp) :: t1(3), t2(3), t3(3)
    do k = 1, nk-1
     do j = 1, nj-1
      do i = 1, ni-1
         t1 = fc(i+1,j,k,2,3) - fc(i,j,k,2,3)
         t2 = fc(i,j+1,k,1,3) - fc(i,j,k,1,3)
         t3 = fc(i,j,k+1,1,2) - fc(i,j,k,1,2)
         vol(i,j,k) = abs( t1(1)*(t2(2)*t3(3)-t2(3)*t3(2))                 &
                         - t1(2)*(t2(1)*t3(3)-t2(3)*t3(1))                 &
                         + t1(3)*(t2(1)*t3(2)-t2(2)*t3(1)) )
      end do
     end do
    end do
  contains
    ! centre of the face with corner node (i0,j0,k0) spanning axes a,b
    pure function fc(i0, j0, k0, a, b) result(c)
      integer, intent(in) :: i0, j0, k0, a, b
      real(dp) :: c(3)
      integer :: e(3), g(3), n1(3), n2(3), n3(3), n0(3)
      e = 0;  e(a) = 1;  g = 0;  g(b) = 1
      n0 = [i0,j0,k0];  n1 = n0 + e;  n2 = n0 + g;  n3 = n0 + e + g
      c(1) = 0.25_dp*(x(n0(1),n0(2),n0(3)) + x(n1(1),n1(2),n1(3)) + x(n2(1),n2(2),n2(3)) + x(n3(1),n3(2),n3(3)))
      c(2) = 0.25_dp*(y(n0(1),n0(2),n0(3)) + y(n1(1),n1(2),n1(3)) + y(n2(1),n2(2),n2(3)) + y(n3(1),n3(2),n3(3)))
      c(3) = 0.25_dp*(z(n0(1),n0(2),n0(3)) + z(n1(1),n1(2),n1(3)) + z(n2(1),n2(2),n2(3)) + z(n3(1),n3(2),n3(3)))
    end function fc
  end subroutine cell_volumes

  !=====================================================================
  ! (2) TRIANGLE SURFACES
  !=====================================================================
  subroutine surf_init(s, nvmax, ntmax)
    type(surf_t), intent(out) :: s
    integer, intent(in) :: nvmax, ntmax
    allocate(s%v(3,max(nvmax,1)), s%vn(3,max(nvmax,1)), s%t(3,max(ntmax,1)),  &
             s%tn(3,max(ntmax,1)), s%tp(max(ntmax,1)), s%tb(max(ntmax,1)))
    s%nv = 0;  s%nt = 0;  s%tb = 0;  s%ready = .false.
  end subroutine surf_init

  subroutine surf_free(s)
    type(surf_t), intent(inout) :: s
    if (allocated(s%v)) deallocate(s%v, s%vn, s%t, s%tn, s%tp, s%tb)
    if (allocated(s%start)) deallocate(s%start, s%list)
    s%nv = 0;  s%nt = 0;  s%ready = .false.
  end subroutine surf_free

  ! release everything, including the trimming flags and the winding-
  ! number data (re-entrant overset assembly)
  subroutine surf_free_all(s)
    type(surf_t), intent(inout) :: s
    call surf_free(s)
    if (allocated(s%keep)) deallocate(s%keep)
    if (allocated(s%ta)) deallocate(s%ta, s%tcen, s%stamp)
    if (allocated(s%tbk)) deallocate(s%tbk)
    if (allocated(s%wt)) deallocate(s%wt)
    s%wind_ready = .false.;  s%nbody = 0;  s%qid = 0
  end subroutine surf_free_all

  !---------------------------------------------------------------------
  ! Add a structured n1 x n2 node sheet (fx,fy,fz) as one patch 'ip'.
  ! fd = "fluid direction" per node (vector pointing from the surface
  ! into the fluid; used to orient the triangle normals).  Returns the
  ! vertex offset v0 (vertex id of node (m1,m2) = v0 + m1 + (m2-1)*n1).
  !---------------------------------------------------------------------
  subroutine surf_add_patch(s, ip, n1, n2, fx, fy, fz, fd, v0)
    type(surf_t), intent(inout) :: s
    integer,  intent(in) :: ip, n1, n2
    real(dp), intent(in) :: fx(n1,n2), fy(n1,n2), fz(n1,n2), fd(3,n1,n2)
    integer,  intent(out) :: v0
    integer  :: m1, m2, a, b, c, d, t
    real(dp) :: nr(3), fdir(3)
    v0 = s%nv
    do m2 = 1, n2
       do m1 = 1, n1
          s%nv = s%nv + 1
          s%v(:,s%nv) = [fx(m1,m2), fy(m1,m2), fz(m1,m2)]
       end do
    end do
    do m2 = 1, n2-1
       do m1 = 1, n1-1
          a = v0 + m1 + (m2-1)*n1;  b = a + 1;  c = a + n1 + 1;  d = a + n1
          fdir = 0.25_dp*(fd(:,m1,m2) + fd(:,m1+1,m2) + fd(:,m1+1,m2+1) + fd(:,m1,m2+1))
          call add_tri(a, b, c)
          call add_tri(a, c, d)
       end do
    end do
  contains
    subroutine add_tri(p1, p2, p3)
      integer, intent(in) :: p1, p2, p3
      real(dp) :: nl
      nr = cross(s%v(:,p2)-s%v(:,p1), s%v(:,p3)-s%v(:,p1))
      nl = sqrt(sum(nr*nr))
      if (nl <= 0.0_dp) return                    ! degenerate (collapsed edge)
      nr = nr / nl
      s%nt = s%nt + 1;  t = s%nt
      if (sum(nr*fdir) >= 0.0_dp) then
         s%t(:,t) = [p1, p2, p3];  s%tn(:,t) = nr
      else
         s%t(:,t) = [p1, p3, p2];  s%tn(:,t) = -nr
      end if
      s%tp(t) = ip
    end subroutine add_tri
  end subroutine surf_add_patch

  !---------------------------------------------------------------------
  ! Pseudo normals (angle weighted) and the bucket search structure.
  !---------------------------------------------------------------------
  subroutine surf_finish(s)
    type(surf_t), intent(inout) :: s
    integer  :: t, a, n, k, l1(3), l2(3), b1, b2, b3, m, nb
    real(dp) :: e1(3), e2(3), ang, L(3), Lg, cmin(3), cmax(3), nl
    integer, allocatable :: cnt(:)
    integer :: vv(3), prv, nxt
    s%vn = 0.0_dp
    do t = 1, s%nt
       vv = s%t(:,t)
       do k = 1, 3
          prv = vv(mod(k+1,3)+1);  nxt = vv(mod(k,3)+1)
          e1 = s%v(:,nxt) - s%v(:,vv(k));  e2 = s%v(:,prv) - s%v(:,vv(k))
          nl = sqrt(sum(e1*e1))*sqrt(sum(e2*e2))
          if (nl <= 0.0_dp) cycle
          ang = acos(max(-1.0_dp, min(1.0_dp, sum(e1*e2)/nl)))
          s%vn(:,vv(k)) = s%vn(:,vv(k)) + ang*s%tn(:,t)
       end do
    end do
    do n = 1, s%nv
       nl = sqrt(sum(s%vn(:,n)**2))
       if (nl > 0.0_dp) s%vn(:,n) = s%vn(:,n)/nl
    end do
    ! buckets over the triangles
    s%bmin = huge(1.0_dp);  s%bmax = -huge(1.0_dp)
    do n = 1, s%nv
       s%bmin = min(s%bmin, s%v(:,n));  s%bmax = max(s%bmax, s%v(:,n))
    end do
    if (s%nv == 0) then                      ! empty surface (no wall patch): unit box
       s%bmin = 0.0_dp;  s%bmax = 1.0_dp
    end if
    L = s%bmax - s%bmin
    do a = 1, 3
       if (L(a) <= 0.0_dp) L(a) = 1.0e-12_dp*max(1.0_dp, maxval(L))
    end do
    Lg = (L(1)*L(2)*L(3))**(1.0_dp/3.0_dp)
    do a = 1, 3
       s%n(a) = max(1, min(nint(real(max(s%nt,1),dp)**(1.0_dp/3.0_dp)*L(a)/Lg), 2048))
    end do
    do while (product(s%n) > 2*s%nt + 8)
       a = maxloc(s%n, 1)
       s%n(a) = max(1, s%n(a)/2)
    end do
    s%o = s%bmin - 1.0e-6_dp*L
    s%h = (L + 2.0e-6_dp*L) / real(s%n, dp)
    nb = product(s%n)
    if (allocated(s%start)) deallocate(s%start, s%list)
    allocate(cnt(nb), s%start(nb+1));  cnt = 0
    do t = 1, s%nt
       call tri_box(s, t, cmin, cmax);  call sbox_range(s, cmin, cmax, l1, l2)
       do b3 = l1(3), l2(3);  do b2 = l1(2), l2(2);  do b1 = l1(1), l2(1)
          m = b1 + s%n(1)*((b2-1) + s%n(2)*(b3-1));  cnt(m) = cnt(m) + 1
       end do;  end do;  end do
    end do
    s%start(1) = 1
    do m = 1, nb
       s%start(m+1) = s%start(m) + cnt(m)
    end do
    allocate(s%list(max(1, s%start(nb+1)-1)));  cnt = 0
    do t = 1, s%nt
       call tri_box(s, t, cmin, cmax);  call sbox_range(s, cmin, cmax, l1, l2)
       do b3 = l1(3), l2(3);  do b2 = l1(2), l2(2);  do b1 = l1(1), l2(1)
          m = b1 + s%n(1)*((b2-1) + s%n(2)*(b3-1))
          s%list(s%start(m) + cnt(m)) = t;  cnt(m) = cnt(m) + 1
       end do;  end do;  end do
    end do
    deallocate(cnt)
    s%ready = .true.
  end subroutine surf_finish

  pure subroutine tri_box(s, t, cmin, cmax)
    type(surf_t), intent(in) :: s
    integer, intent(in) :: t
    real(dp), intent(out) :: cmin(3), cmax(3)
    cmin = min(s%v(:,s%t(1,t)), s%v(:,s%t(2,t)), s%v(:,s%t(3,t)))
    cmax = max(s%v(:,s%t(1,t)), s%v(:,s%t(2,t)), s%v(:,s%t(3,t)))
  end subroutine tri_box

  pure subroutine sbox_range(s, cmin, cmax, l1, l2)
    type(surf_t), intent(in) :: s
    real(dp), intent(in) :: cmin(3), cmax(3)
    integer, intent(out) :: l1(3), l2(3)
    integer :: a
    do a = 1, 3
       l1(a) = max(1, min(s%n(a), int(floor((cmin(a) - s%o(a))/s%h(a))) + 1))
       l2(a) = max(1, min(s%n(a), int(floor((cmax(a) - s%o(a))/s%h(a))) + 1))
    end do
  end subroutine sbox_range

  !---------------------------------------------------------------------
  ! Closest point on triangle (a,b,c) to p (Ericson, RTCD 5.1.5).
  ! region: 1 = interior, 2 = edge, 3 = vertex; the barycentric weights
  ! w(3) tell which vertices/edges are involved.
  !---------------------------------------------------------------------
  pure subroutine tri_closest(p, a, b, c, q, w, region)
    real(dp), intent(in)  :: p(3), a(3), b(3), c(3)
    real(dp), intent(out) :: q(3), w(3)
    integer,  intent(out) :: region
    real(dp) :: ab(3), ac(3), ap(3), bp(3), cp(3), d1, d2, d3, d4, d5, d6
    real(dp) :: vc, vb, va, v, ww, denom
    ab = b - a;  ac = c - a;  ap = p - a
    d1 = sum(ab*ap);  d2 = sum(ac*ap)
    if (d1 <= 0.0_dp .and. d2 <= 0.0_dp) then
       q = a;  w = [1.0_dp, 0.0_dp, 0.0_dp];  region = 3;  return
    end if
    bp = p - b;  d3 = sum(ab*bp);  d4 = sum(ac*bp)
    if (d3 >= 0.0_dp .and. d4 <= d3) then
       q = b;  w = [0.0_dp, 1.0_dp, 0.0_dp];  region = 3;  return
    end if
    vc = d1*d4 - d3*d2
    if (vc <= 0.0_dp .and. d1 >= 0.0_dp .and. d3 <= 0.0_dp) then
       v = d1 / (d1 - d3);  q = a + v*ab;  w = [1.0_dp-v, v, 0.0_dp];  region = 2;  return
    end if
    cp = p - c;  d5 = sum(ab*cp);  d6 = sum(ac*cp)
    if (d6 >= 0.0_dp .and. d5 <= d6) then
       q = c;  w = [0.0_dp, 0.0_dp, 1.0_dp];  region = 3;  return
    end if
    vb = d5*d2 - d1*d6
    if (vb <= 0.0_dp .and. d2 >= 0.0_dp .and. d6 <= 0.0_dp) then
       ww = d2 / (d2 - d6);  q = a + ww*ac;  w = [1.0_dp-ww, 0.0_dp, ww];  region = 2;  return
    end if
    va = d3*d6 - d5*d4
    if (va <= 0.0_dp .and. (d4 - d3) >= 0.0_dp .and. (d5 - d6) >= 0.0_dp) then
       ww = (d4 - d3) / ((d4 - d3) + (d5 - d6))
       q = b + ww*(c - b);  w = [0.0_dp, 1.0_dp-ww, ww];  region = 2;  return
    end if
    denom = 1.0_dp / (va + vb + vc)
    v = vb*denom;  ww = vc*denom
    q = a + ab*v + ac*ww;  w = [1.0_dp-v-ww, v, ww];  region = 1
  end subroutine tri_closest

  !---------------------------------------------------------------------
  ! Closest surface point to p over the triangles selected by the
  ! filter (mode 0: all; 1: exclude patch 'sel'; 2: only body 'sel').
  ! Returns the distance, the closest point q, the triangle tbest and
  ! its region / barycentric weights.  dmax bounds the search (only
  ! triangles within dmax are considered; dist = huge if none).
  !---------------------------------------------------------------------
  subroutine surf_closest(s, p, mode, sel, dmax, dist, q, tbest, wb, region)
    type(surf_t), intent(in)  :: s
    real(dp),     intent(in)  :: p(3), dmax
    integer,      intent(in)  :: mode, sel
    real(dp),     intent(out) :: dist, q(3), wb(3)
    integer,      intent(out) :: tbest, region
    integer  :: c(3), ring, b1, b2, b3, m, n, t, rg, maxring, a
    real(dp) :: qq(3), w(3), d, dring, lim
    logical  :: any_in

    dist = huge(1.0_dp);  q = p;  tbest = 0;  wb = 0.0_dp;  region = 0
    if (.not. s%ready .or. s%nt == 0) return
    do a = 1, 3
       c(a) = int(floor((p(a) - s%o(a))/s%h(a))) + 1
    end do
    maxring = maxval(s%n) + maxval(abs(c)) + 2
    lim = dmax
    do ring = 0, maxring
       ! distance from p to the nearest point outside the (ring-1) shell
       dring = 0.0_dp
       if (ring > 0) dring = ring_dist(ring)
       if (dring > min(dist, lim)) exit
       any_in = .false.
       do b3 = c(3)-ring, c(3)+ring
        do b2 = c(2)-ring, c(2)+ring
         do b1 = c(1)-ring, c(1)+ring
            if (max(abs(b1-c(1)), abs(b2-c(2)), abs(b3-c(3))) /= ring) cycle
            if (b1 < 1 .or. b1 > s%n(1) .or. b2 < 1 .or. b2 > s%n(2) .or. b3 < 1 .or. b3 > s%n(3)) cycle
            any_in = .true.
            m = b1 + s%n(1)*((b2-1) + s%n(2)*(b3-1))
            do n = s%start(m), s%start(m+1)-1
               t = s%list(n)
               if (mode == 1 .and. s%tp(t) == sel) cycle
               if (mode == 2 .and. s%tb(t) /= sel) cycle
               if (mode == 3 .and. s%tp(t) /= sel) cycle       ! only patch sel
               call tri_closest(p, s%v(:,s%t(1,t)), s%v(:,s%t(2,t)), s%v(:,s%t(3,t)), qq, w, rg)
               d = sqrt(sum((p-qq)**2))
               if (d < dist) then
                  dist = d;  q = qq;  tbest = t;  wb = w;  region = rg
               end if
            end do
         end do
        end do
       end do
       if (.not. any_in .and. ring > 0) then
          ! the whole ring is outside the bucket grid: stop when the ring
          ! distance exceeds the best found, or the grid is exhausted
          if (ring > maxval(s%n) + 1) exit
       end if
    end do
  contains
    ! lower bound of the distance from p to any point in ring 'rg'
    real(dp) function ring_dist(rg)
      integer, intent(in) :: rg
      real(dp) :: dd, lo, hi
      integer :: aa
      ring_dist = huge(1.0_dp)
      do aa = 1, 3
         lo = s%o(aa) + real(c(aa)-rg, dp)*s%h(aa)
         hi = s%o(aa) + real(c(aa)+rg-1, dp)*s%h(aa)
         dd = min(p(aa) - lo, hi - p(aa))
         ring_dist = min(ring_dist, dd)
      end do
      ring_dist = max(ring_dist, 0.0_dp)
    end function ring_dist
  end subroutine surf_closest

  !---------------------------------------------------------------------
  ! Inside test of p with respect to the closed surface of body 'ib'
  ! (pseudo-normal sign test at the closest point):  sd < 0 inside.
  ! Returns the signed distance.
  !---------------------------------------------------------------------
  subroutine surf_inside(s, p, ib, sd, dist)
    type(surf_t), intent(in) :: s
    real(dp), intent(in)  :: p(3)
    integer,  intent(in)  :: ib
    real(dp), intent(out) :: sd, dist
    real(dp) :: q(3), wb(3), nr(3), nl
    integer  :: t, region, k
    call surf_closest(s, p, 2, ib, huge(1.0_dp), dist, q, t, wb, region)
    if (t == 0) then
       sd = huge(1.0_dp);  return
    end if
    select case (region)
    case (1)
       nr = s%tn(:,t)
    case default
       nr = 0.0_dp
       do k = 1, 3
          if (wb(k) > 0.0_dp) nr = nr + s%vn(:,s%t(k,t))
       end do
       nl = sqrt(sum(nr*nr))
       if (nl > 0.0_dp) then
          nr = nr / nl
       else
          nr = s%tn(:,t)
       end if
    end select
    sd = sum((p - q)*nr)
    if (abs(sd) < dist) sd = sign(dist, sd)
  end subroutine surf_inside

  pure function cross(a, b) result(c)
    real(dp), intent(in) :: a(3), b(3)
    real(dp) :: c(3)
    c = [a(2)*b(3)-a(3)*b(2), a(3)*b(1)-a(1)*b(3), a(1)*b(2)-a(2)*b(1)]
  end function cross


  !=====================================================================
  ! Generalized winding number of the wall surface (Jacobson et al.):
  !   w(p) = sum_t Omega_t(p) / 4 pi  with the solid angle Omega_t signed
  !   by the triangle's fluid-side normal tn (positive when tn points
  !   away from p).  Closed surface: w = 1 inside, 0 outside, robust to
  !   overlapping sheets and small gaps.  Only triangles with keep =
  !   .true. and body id ib contribute.
  !=====================================================================
  ! (4) GENERALIZED WINDING NUMBER (per-body cluster trees)
  !=====================================================================
  ! surf_wind_prep builds, for every body, a cluster tree over the kept
  ! triangles of that body: a base grid over the body's own bounding box
  ! (about nt^(1/3) cells per axis, cells at least as large as the mean
  ! triangle), each triangle in the cell of its centroid, and levels of
  ! 2x coarsening above it.  Every cluster carries the area-weighted
  ! dipole sum(a n), the area centroid c, the area, the bounding radius
  ! of its vertices about c, and the second moment sum(a n (c_t - c)^T).
  subroutine surf_wind_prep(s, nbody)
    type(surf_t), intent(inout) :: s
    integer, intent(in) :: nbody
    integer :: t, m, nb, ib, b(3), mc, k, l, ncl, nl(3), nc(3), k1, k2, k3, nt_, a
    real(dp) :: v1(3), v2(3), v3(3), ng(3), ar, cen(3), nrm(3), bmin(3), bmax(3), bl(3), blg, emean
    integer, allocatable :: fill(:)
    if (.not. s%ready) return
    s%nbody = nbody
    if (.not. allocated(s%keep)) then
       allocate(s%keep(s%nt));  s%keep = .true.
    end if
    if (allocated(s%ta)) deallocate(s%ta, s%tcen, s%stamp)
    if (allocated(s%tbk)) deallocate(s%tbk)
    allocate(s%ta(s%nt), s%tcen(3,s%nt), s%stamp(s%nt), s%tbk(s%nt));  s%stamp = 0;  s%qid = 0;  s%tbk = 0
    do t = 1, s%nt
       v1 = s%v(:,s%t(1,t));  v2 = s%v(:,s%t(2,t));  v3 = s%v(:,s%t(3,t))
       ng = cross3(v2 - v1, v3 - v1)
       s%ta(t) = 0.5_dp*sqrt(sum(ng*ng))
       s%tcen(:,t) = (v1 + v2 + v3)/3.0_dp
    end do
    if (allocated(s%wt)) deallocate(s%wt)
    allocate(s%wt(max(nbody,1)))
    do ib = 1, nbody
       associate (w => s%wt(ib))
       ! ---- base grid over the body's bounding box ------------------
       bmin = huge(1.0_dp);  bmax = -huge(1.0_dp);  nt_ = 0;  emean = 0.0_dp
       do t = 1, s%nt
          if (.not. s%keep(t) .or. s%tb(t) /= ib) cycle
          v1 = s%v(:,s%t(1,t));  v2 = s%v(:,s%t(2,t));  v3 = s%v(:,s%t(3,t))
          bmin = min(bmin, v1, v2, v3);  bmax = max(bmax, v1, v2, v3)
          emean = emean + max(sqrt(sum((v2-v1)**2)), sqrt(sum((v3-v2)**2)), sqrt(sum((v1-v3)**2)))
          nt_ = nt_ + 1
       end do
       w%nt = nt_
       if (nt_ == 0) cycle
       emean = emean/nt_
       bl = bmax - bmin
       do a = 1, 3
          if (bl(a) <= 0.0_dp) bl(a) = 1.0e-12_dp*max(1.0_dp, maxval(bl))
       end do
       blg = (bl(1)*bl(2)*bl(3))**(1.0_dp/3.0_dp)
       do a = 1, 3
          w%n(a) = max(1, min(nint(real(nt_,dp)**(1.0_dp/3.0_dp)*bl(a)/blg), 1024))
          ! cells at least as large as the mean triangle (longest edge):
          ! with tiny cells under large triangles every cell's radius spans
          ! many cells and the exact near-field sums explode
          do while (w%n(a) > 1 .and. bl(a)/w%n(a) < emean)
             w%n(a) = max(1, w%n(a)/2)
          end do
       end do
       w%o = bmin - 1.0e-6_dp*bl;  w%h = (bl + 2.0e-6_dp*bl)/real(w%n, dp)
       nb = product(w%n)
       allocate(w%fd(3,nb), w%fc(3,nb), w%fa(nb), w%fr(nb), w%fq(3,3,nb), w%cstart(nb+1), w%clist(nt_))
       w%fd = 0.0_dp;  w%fc = 0.0_dp;  w%fa = 0.0_dp;  w%fr = 0.0_dp;  w%fq = 0.0_dp;  w%cstart = 0
       ! ---- triangles -> base cells (centroid) ------------------------
       do t = 1, s%nt
          if (.not. s%keep(t) .or. s%tb(t) /= ib) cycle
          cen = s%tcen(:,t);  ar = s%ta(t);  nrm = s%tn(:,t)
          do a = 1, 3
             b(a) = min(max(int((cen(a)-w%o(a))/w%h(a)) + 1, 1), w%n(a))
          end do
          m = b(1) + w%n(1)*((b(2)-1) + w%n(2)*(b(3)-1))
          s%tbk(t) = m
          w%fd(:,m) = w%fd(:,m) + ar*nrm
          w%fc(:,m) = w%fc(:,m) + ar*cen
          w%fa(m)   = w%fa(m) + ar
          w%cstart(m+1) = w%cstart(m+1) + 1
       end do
       do m = 1, nb
          if (w%fa(m) > 0.0_dp) w%fc(:,m) = w%fc(:,m)/w%fa(m)
       end do
       w%cstart(1) = 1
       do m = 1, nb
          w%cstart(m+1) = w%cstart(m+1) + w%cstart(m)
       end do
       allocate(fill(nb));  fill = 0
       do t = 1, s%nt
          if (.not. s%keep(t) .or. s%tb(t) /= ib) cycle
          m = s%tbk(t)
          w%clist(w%cstart(m) + fill(m)) = t;  fill(m) = fill(m) + 1
          ! radius (farthest vertex) and second moment about the cell centroid
          do k = 1, 3
             v1 = s%v(:,s%t(k,t))
             w%fr(m) = max(w%fr(m), sqrt(sum((v1 - w%fc(:,m))**2)))
          end do
          cen = s%tcen(:,t) - w%fc(:,m)
          do k = 1, 3
             w%fq(:,k,m) = w%fq(:,k,m) + s%ta(t)*s%tn(:,t)*cen(k)
          end do
       end do
       deallocate(fill)
       ! ---- levels: 2x coarsening up to 1-2 clusters per axis ---------
       w%nlev = 0
       nl = w%n
       do while (maxval(nl) > 2)
          nl = (nl + 1)/2;  w%nlev = w%nlev + 1
       end do
       allocate(w%lev(max(w%nlev,1)))
       nl = w%n
       do l = w%nlev, 1, -1
          nl = (nl + 1)/2
          w%lev(l)%n = nl
          ncl = product(nl)
          allocate(w%lev(l)%d(3,ncl), w%lev(l)%c(3,ncl), w%lev(l)%a(ncl), w%lev(l)%r(ncl), w%lev(l)%q(3,3,ncl))
          w%lev(l)%d = 0.0_dp;  w%lev(l)%c = 0.0_dp;  w%lev(l)%a = 0.0_dp;  w%lev(l)%r = 0.0_dp;  w%lev(l)%q = 0.0_dp
       end do
       do l = w%nlev, 1, -1
          if (l == w%nlev) then
             nc = w%n
          else
             nc = w%lev(l+1)%n
          end if
          do k3 = 1, nc(3);  do k2 = 1, nc(2);  do k1 = 1, nc(1)
             m  = k1 + nc(1)*((k2-1) + nc(2)*(k3-1))
             mc = (k1+1)/2 + w%lev(l)%n(1)*(((k2+1)/2-1) + w%lev(l)%n(2)*((k3+1)/2-1))
             if (l == w%nlev) then
                ar = w%fa(m);  if (ar <= 0.0_dp) cycle
                w%lev(l)%d(:,mc) = w%lev(l)%d(:,mc) + w%fd(:,m)
                w%lev(l)%c(:,mc) = w%lev(l)%c(:,mc) + ar*w%fc(:,m)
             else
                ar = w%lev(l+1)%a(m);  if (ar <= 0.0_dp) cycle
                w%lev(l)%d(:,mc) = w%lev(l)%d(:,mc) + w%lev(l+1)%d(:,m)
                w%lev(l)%c(:,mc) = w%lev(l)%c(:,mc) + ar*w%lev(l+1)%c(:,m)
             end if
             w%lev(l)%a(mc) = w%lev(l)%a(mc) + ar
          end do;  end do;  end do
          do mc = 1, product(w%lev(l)%n)
             if (w%lev(l)%a(mc) > 0.0_dp) w%lev(l)%c(:,mc) = w%lev(l)%c(:,mc)/w%lev(l)%a(mc)
          end do
          ! radius = max over the children of (centroid distance + child
          ! radius); second moment about the parent centroid:
          !   sum_t a_t n_t (c_t - cp)^T = Q_child + d_child (c_child - cp)^T
          do k3 = 1, nc(3);  do k2 = 1, nc(2);  do k1 = 1, nc(1)
             m  = k1 + nc(1)*((k2-1) + nc(2)*(k3-1))
             mc = (k1+1)/2 + w%lev(l)%n(1)*(((k2+1)/2-1) + w%lev(l)%n(2)*((k3+1)/2-1))
             if (l == w%nlev) then
                if (w%fa(m) <= 0.0_dp) cycle
                cen = w%fc(:,m) - w%lev(l)%c(:,mc)
                w%lev(l)%r(mc) = max(w%lev(l)%r(mc), sqrt(sum(cen*cen)) + w%fr(m))
                do k = 1, 3
                   w%lev(l)%q(:,k,mc) = w%lev(l)%q(:,k,mc) + w%fq(:,k,m) + w%fd(:,m)*cen(k)
                end do
             else
                if (w%lev(l+1)%a(m) <= 0.0_dp) cycle
                cen = w%lev(l+1)%c(:,m) - w%lev(l)%c(:,mc)
                w%lev(l)%r(mc) = max(w%lev(l)%r(mc), sqrt(sum(cen*cen)) + w%lev(l+1)%r(m))
                do k = 1, 3
                   w%lev(l)%q(:,k,mc) = w%lev(l)%q(:,k,mc) + w%lev(l+1)%q(:,k,m) + w%lev(l+1)%d(:,m)*cen(k)
                end do
             end if
          end do;  end do;  end do
       end do
       end associate
    end do
    s%wind_ready = .true.
  end subroutine surf_wind_prep

  ! generalized winding number of point p with respect to the (kept)
  ! triangles of body ib: exact solid angles (van Oosterom - Strackee)
  ! for the triangles near p, multipole sums (dipole + second moment)
  ! for the clusters far from p, descending the body's cluster tree from
  ! the coarsest level.  A cluster is "far" when the distance of p from
  ! its centroid exceeds BETA times its radius (bounding sphere of its
  ! vertices about the centroid) - the far-field criterion of Barill et
  ! al. (2018, "Fast winding numbers"); otherwise its 8 children are
  ! visited, and at the base level the triangles of a near cell are
  ! summed exactly.  (A fixed bucket-index distance, used before update
  ! 12c, failed for thin bodies with triangles larger than the buckets:
  ! DelftTwist foil, w = 0.4 - 0.7 next to the surface, 1865 O-grid cells
  ! cut as "inside".)
  real(dp) function surf_wind(s, p, ib) result(w)
    type(surf_t), intent(inout) :: s
    real(dp), intent(in) :: p(3)
    integer, intent(in)  :: ib
    real(dp), parameter :: BETA = 2.0_dp
    integer  :: m, t, l, lv, k1, k2, k3, c(3), nsp, nc(3), i1, i2, i3
    integer  :: stk(4, 4096)                     ! (level, k1, k2, k3): level 0 = base cells
    real(dp) :: om, d(3), r2, a(3), b(3), cv(3), la, lb, lc, det, den, ng(3), v1(3), v2(3), v3(3)
    real(dp), parameter :: fourpi = 12.566370614359172_dp
    w = 0.0_dp
    if (.not. s%wind_ready) return
    if (ib < 1 .or. ib > size(s%wt)) return
    if (s%wt(ib)%nt == 0) return
    wind_nq = wind_nq + 1
    associate (wt => s%wt(ib))
    om = 0.0_dp
    nsp = 0
    if (wt%nlev == 0) then                       ! tiny surface: all base cells
       do k3 = 1, wt%n(3);  do k2 = 1, wt%n(2);  do k1 = 1, wt%n(1)
          nsp = nsp + 1;  stk(:,nsp) = [0, k1, k2, k3]
       end do;  end do;  end do
    else
       do k3 = 1, wt%lev(1)%n(3);  do k2 = 1, wt%lev(1)%n(2);  do k1 = 1, wt%lev(1)%n(1)
          nsp = nsp + 1;  stk(:,nsp) = [1, k1, k2, k3]
       end do;  end do;  end do
    end if
    do while (nsp > 0)
       lv = stk(1,nsp);  c = stk(2:4,nsp);  nsp = nsp - 1
       wind_ncl = wind_ncl + 1
       if (lv == 0) then
          m = c(1) + wt%n(1)*((c(2)-1) + wt%n(2)*(c(3)-1))
          if (wt%fa(m) <= 0.0_dp) cycle
          d = wt%fc(:,m) - p;  r2 = sum(d*d)
          if (r2 > (BETA*wt%fr(m))**2) then               ! far base cell: multipole
             om = om + far_term(wt%fd(:,m), wt%fq(:,:,m), d, r2)
             cycle
          end if
          ! near: every triangle whose centroid lies in this cell, exactly
          do l = wt%cstart(m), wt%cstart(m+1)-1
             t = wt%clist(l)
             wind_ntri = wind_ntri + 1
             v1 = s%v(:,s%t(1,t));  v2 = s%v(:,s%t(2,t));  v3 = s%v(:,s%t(3,t))
             a = v1 - p;  b = v2 - p;  cv = v3 - p
             la = sqrt(sum(a*a));  lb = sqrt(sum(b*b));  lc = sqrt(sum(cv*cv))
             det = sum(a*cross3(b, cv))
             den = la*lb*lc + sum(a*b)*lc + sum(a*cv)*lb + sum(b*cv)*la
             ng  = cross3(v2 - v1, v3 - v1)
             om  = om + 2.0_dp*atan2(det, den)*sign(1.0_dp, sum(ng*s%tn(:,t)))
          end do
       else
          m = c(1) + wt%lev(lv)%n(1)*((c(2)-1) + wt%lev(lv)%n(2)*(c(3)-1))
          if (wt%lev(lv)%a(m) <= 0.0_dp) cycle
          d = wt%lev(lv)%c(:,m) - p;  r2 = sum(d*d)
          if (r2 > (BETA*wt%lev(lv)%r(m))**2) then        ! far cluster: multipole
             om = om + far_term(wt%lev(lv)%d(:,m), wt%lev(lv)%q(:,:,m), d, r2)
             cycle
          end if
          ! near: push the (up to 8) children
          if (lv == wt%nlev) then
             nc = wt%n
          else
             nc = wt%lev(lv+1)%n
          end if
          do i3 = 2*c(3)-1, min(2*c(3), nc(3))
           do i2 = 2*c(2)-1, min(2*c(2), nc(2))
            do i1 = 2*c(1)-1, min(2*c(1), nc(1))
               if (nsp >= size(stk,2)) cycle               ! cannot happen for a depth-first descent
               nsp = nsp + 1;  stk(:,nsp) = [merge(0, lv+1, lv == wt%nlev), i1, i2, i3]
            end do
           end do
          end do
       end if
    end do
    end associate
    w = om/fourpi
  contains
    ! solid angle of a far cluster: Taylor expansion of f(x) = (x - p)/|x - p|^3
    ! about the cluster centroid c (d = c - p):  n.f(c) + sum_ij G_ij Q_ij,
    ! G = grad f(c) = I/|d|^3 - 3 d d^T/|d|^5  (the first two orders of
    ! Barill et al. 2018)
    pure real(dp) function far_term(dn, qm, d, r2)
      real(dp), intent(in) :: dn(3), qm(3,3), d(3), r2
      real(dp) :: r3, r5, g(3,3)
      integer :: i_, j_
      r3 = r2*sqrt(r2);  r5 = r3*r2
      do j_ = 1, 3
         do i_ = 1, 3
            g(i_,j_) = -3.0_dp*d(i_)*d(j_)/r5
         end do
         g(j_,j_) = g(j_,j_) + 1.0_dp/r3
      end do
      far_term = sum(dn*d)/r3 + sum(g*qm)
    end function far_term
  end function surf_wind

  pure function cross3(u, v) result(r)
    real(dp), intent(in) :: u(3), v(3)
    real(dp) :: r(3)
    r = [u(2)*v(3)-u(3)*v(2), u(3)*v(1)-u(1)*v(3), u(1)*v(2)-u(2)*v(1)]
  end function cross3


  !=====================================================================
  ! PLIC kernels (unit cube, normal m >= 0 with sum m = 1, liquid side
  ! m.x <= a), shared by vof.F90 and by the overset fringe update.
  !=====================================================================
  ! liquid fraction of the unit cube cut by m.x <= a
  pure real(dp) function plic_cube_frac(m, a) result(v)
    real(dp), intent(in) :: m(3), a
    real(dp) :: t
    integer :: l
    if (a <= 0.0_dp) then
       v = 0.0_dp;  return
    else if (a >= 1.0_dp) then
       v = 1.0_dp;  return
    end if
    v = a**3
    do l = 1, 3
       t = a - m(l);          if (t > 0.0_dp) v = v - t**3
       t = a - 1.0_dp + m(l); if (t > 0.0_dp) v = v + t**3
    end do
    v = v/(6.0_dp*m(1)*m(2)*m(3))
    v = min(max(v, 0.0_dp), 1.0_dp)
  end function plic_cube_frac

  ! plane constant of the liquid fraction a (bisection)
  pure real(dp) function plic_plane_const(m, a) result(acst)
    real(dp), intent(in) :: m(3), a
    real(dp) :: lo, hi, mid
    integer :: it
    lo = 0.0_dp;  hi = 1.0_dp
    do it = 1, 48
       mid = 0.5_dp*(lo + hi)
       if (plic_cube_frac(m, mid) < a) then
          lo = mid
       else
          hi = mid
       end if
    end do
    acst = 0.5_dp*(lo + hi)
  end function plic_plane_const

  ! liquid fraction of the box [x0, x0+w] (per axis, in the cube's
  ! coordinates; the plane m.x <= a is extended beyond the cube)
  pure real(dp) function plic_box_frac(m, acst, x0, w) result(v)
    real(dp), intent(in) :: m(3), acst, x0(3), w(3)
    real(dp) :: mp(3), ap, s
    if (any(w <= 0.0_dp)) then
       v = 0.0_dp;  return
    end if
    mp = m*w
    s  = sum(mp)
    ap = (acst - sum(m*x0))/s
    mp = mp/s
    mp = max(mp, 1.0e-9_dp);  mp = mp/sum(mp)
    v  = plic_cube_frac(mp, ap)
  end function plic_box_frac

  ! Youngs normal from the 3x3x3 stencil c(-1:1,-1:1,-1:1) of alpha:
  ! m = normal from liquid to gas, mirrored so that m >= 0 (mir), sum 1
  pure subroutine plic_normal27(c, m, mir)
    real(dp), intent(in)  :: c(-1:1,-1:1,-1:1)
    real(dp), intent(out) :: m(3)
    logical,  intent(out) :: mir(3)
    real(dp), parameter :: wgt(-1:1) = [0.25_dp, 0.5_dp, 0.25_dp]
    real(dp) :: g(3), s
    integer :: d1, d2
    g = 0.0_dp
    do d2 = -1, 1
       do d1 = -1, 1
          g(1) = g(1) + wgt(d1)*wgt(d2)*(c( 1,d1,d2) - c(-1,d1,d2))
          g(2) = g(2) + wgt(d1)*wgt(d2)*(c(d1, 1,d2) - c(d1,-1,d2))
          g(3) = g(3) + wgt(d1)*wgt(d2)*(c(d1,d2, 1) - c(d1,d2,-1))
       end do
    end do
    m = -0.5_dp*g
    s = abs(m(1)) + abs(m(2)) + abs(m(3))
    if (s < 1.0e-30_dp) then
       m = [0.0_dp, 0.0_dp, 1.0_dp];  s = 1.0_dp
    end if
    mir = (m < 0.0_dp)
    m = abs(m)/s
    m = max(m, 1.0e-9_dp);  m = m/sum(m)
  end subroutine plic_normal27
end module ovs_geom
