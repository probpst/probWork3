!=======================================================================
! bcs_b2b.F90
!   Step 1 of the multi-block extension of probCFD:
!
!   Detect block-to-block (patched, type #92) interfaces in a
!   multi-block PLOT3D grid 'casename.grd' and export them to
!   'casename.bcs_b2b'.
!
!   Method
!   ------
!   Two block faces are "connected" when every node of one face patch
!   coincides with a node of the other face patch within a tolerance
!   b2b_tol, applied to x, y and z separately:
!        |x1-x2| <= tol  .and.  |y1-y2| <= tol  .and.  |z1-z2| <= tol
!
!   For every block, all six faces (i-min, i-max, j-min, j-max,
!   k-min, k-max) are tested against all faces of all other blocks
!   (and against the other faces of the same block, so O-grid /
!   periodic-in-index topologies are found as well):
!
!     1. fast bounding-box overlap test (boxes inflated by tol),
!     2. every node of the receiver face is searched in a spatial
!        hash table of the donor-face nodes (tolerance aware),
!     3. the matched nodes are split into 4-connected components; each
!        must fill a full rectangle in the receiver face index space
!        (>= 2 x 2 nodes, so mere edge- or corner-touching blocks are
!        NOT reported as interfaces).  Several components = several
!        interface windows with the same neighbour (e.g. a face that
!        meets its neighbour on both sides of a body),
!     4. the receiver->donor index correspondence must be affine,
!        i.e. a signed permutation of the two tangential indices;
!        this determines the relative orientation of the two faces.
!
!   Partial-face contact is supported: the exported index ranges are
!   the actual overlap window, which may be a sub-range of the face.
!
!   POLE (degenerate) faces: before the pair search every face is
!   tested for degeneracy - all its nodes along one tangential index
!   coincide (within the same tolerance b2b_tol), i.e. the face has
!   collapsed to a line (the axis of an O-grid pipe / polar grid) or
!   to a point.  Such a face is a "pole" boundary: it is excluded from
!   the interface matching and exported as a type-93 record (same
!   12-line layout, donor = the face itself, map 1 2 3).  The solver
!   treats a pole face as a zero-area boundary: zero-gradient ghost
!   values for all variables, zero face flux (the face Jacobian and
!   metrics are zero there), no pressure flux; the pole face needs NO
!   record in the .bcs file.
!
!   File format of 'casename.bcs_b2b'
!   ---------------------------------
!   Modeled on the type-92 (patched) records of the example .bcs
!   files, plus one trailing "index map" line that makes the record
!   unambiguous even when the tangential axes are permuted between
!   the two blocks:
!
!     nblk                 ! number of grid blocks
!     tol                  ! tolerance used for the search
!     --- for each block b = 1..nblk ---
!     n_b2b                ! number of b2b + pole surfaces of block b
!       --- for each surface ---
!       92                 ! bc type: 92 patched block-to-block, 93 pole
!       dir                ! normal direction: |dir|=1,2,3 (i,j,k);
!                          !   + = min face, - = max face
!       i1  i2             ! receiver node range, i
!       j1  j2             ! receiver node range, j
!       k1  k2             ! receiver node range, k
!       0                  ! flag for order of first derivative bc
!       dblk               ! donor block number
!       ddir               ! donor face normal direction (same coding)
!       di1 di2            ! donor node range, i   (start > end when the
!       dj1 dj2            ! donor node range, j    running direction is
!       dk1 dk2            ! donor node range, k    reversed, as in the
!                          !                        example .bcs files)
!       t1  t2  t3         ! index map: receiver axis a (=1,2,3) runs
!                          !   along signed donor axis t(a); e.g.
!                          !   "1 2 3" = same orientation,
!                          !   "2 -1 3" = receiver i runs along donor j,
!                          !              receiver j along donor -i, ...
!
!   The donor range of a donor axis goes from the donor index that
!   corresponds to the START corner of the receiver window to the
!   donor index that corresponds to the END of the receiver range that
!   maps onto that donor axis.  Records are written twice, once from
!   each side of the interface, exactly like in the example files.
!=======================================================================
module b2b_mod
  implicit none
  private
  public :: build_bcs_b2b, b2b_read_file
  public :: n_b2b_tot, b2b_blk, b2b_dir, b2b_irng, b2b_jrng, b2b_krng,  &
            b2b_dblk, b2b_ddir, b2b_dirng, b2b_djrng, b2b_dkrng, b2b_map
  public :: n_pole, pole_blk, pole_dir, pole_rng

  integer, parameter :: dp = selected_real_kind(15, 307)

  !--- one block of the PLOT3D grid ------------------------------------
  type :: blk_t
     integer :: ni, nj, nk
     real(dp), allocatable :: x(:,:,:), y(:,:,:), z(:,:,:)
  end type blk_t

  !--- one face of a block, unpacked to a 2-D node sheet ---------------
  ! tangential axes t1 < t2; node sheet f*(1:n1,1:n2);
  ! spatial hash table over the sheet nodes for tolerance searches
  type :: face_t
     integer :: blk                    ! parent block
     integer :: axis                   ! normal axis 1,2,3
     integer :: side                   ! 1 = min face, 2 = max face
     integer :: dir                    ! signed code: +axis (min), -axis (max)
     integer :: fidx                   ! node index on the normal axis (1 or n)
     integer :: t1, t2                 ! tangential axes (ascending)
     integer :: n1, n2                 ! node counts along t1, t2
     real(dp), allocatable :: fx(:,:), fy(:,:), fz(:,:)
     real(dp) :: bmin(3), bmax(3)      ! bounding box
     integer :: hsize                  ! hash table size (power of two)
     integer, allocatable :: hhead(:)  ! hash heads (0 = empty)
     integer, allocatable :: hnext(:)  ! chain links over linear node index
     logical :: pole = .false.         ! degenerate face (pole line / point)
     integer :: pole_axis = 0          ! tangential axis along which the
                                       !   nodes coincide (0 = both: point)
  end type face_t

  !--- one detected block-to-block patch --------------------------------
  type :: patch_t
     integer :: typ = 92               ! 92 = patched b2b, 93 = pole
     integer :: blk, dir
     integer :: irng(2), jrng(2), krng(2)
     integer :: dblk, ddir
     integer :: dirng(2), djrng(2), dkrng(2)
     integer :: map(3)
  end type patch_t

  !--- connectivity as read back by b2b_read_file (for the solver) -----
  integer :: n_b2b_tot = 0
  integer, allocatable :: b2b_blk(:), b2b_dir(:)
  integer, allocatable :: b2b_irng(:,:), b2b_jrng(:,:), b2b_krng(:,:)
  integer, allocatable :: b2b_dblk(:), b2b_ddir(:)
  integer, allocatable :: b2b_dirng(:,:), b2b_djrng(:,:), b2b_dkrng(:,:)
  integer, allocatable :: b2b_map(:,:)
  !--- pole (degenerate) faces as read back (type 93 records) -----------
  integer :: n_pole = 0
  integer, allocatable :: pole_blk(:), pole_dir(:), pole_rng(:,:,:)   ! (2,3,n)

contains

  !=====================================================================
  ! Read 'gridfile' (multi-block ASCII PLOT3D), find all block-to-block
  ! interfaces within tolerance 'tol' and write them to 'outfile'.
  ! verbose  : print a summary to stdout
  ! npatches : number of one-sided patch records written
  ! ierr     : 0 = ok, /=0 = error (message already printed)
  !=====================================================================
  subroutine build_bcs_b2b(gridfile, outfile, tol, verbose, npatches, ierr)
    character(len=*), intent(in)  :: gridfile, outfile
    real(dp),         intent(in)  :: tol
    logical,          intent(in)  :: verbose
    integer,          intent(out) :: npatches, ierr

    type(blk_t),   allocatable :: blk(:)
    type(face_t),  allocatable :: fac(:)
    type(patch_t), allocatable :: pat(:), tmp(:)
    integer :: nblk, nfac, npat
    integer :: b, f, fa, fb, ax, sd
    integer, allocatable :: nb2b(:)

    npatches = 0
    ierr     = 0

    !------------------------------------------------- read the grid --
    call read_plot3d(gridfile, blk, nblk, ierr)
    if (ierr /= 0) return
    if (verbose) then
       write(*,'(A)')      '  ---- block-to-block interface search ----'
       write(*,'(A,A)')    '  grid file       : ', trim(gridfile)
       write(*,'(A,I0)')   '  grid blocks     : ', nblk
       write(*,'(A,ES12.4)') '  tolerance       : ', tol
    end if

    !----------------------------------- unpack the 6 faces per block --
    nfac = 6*nblk
    allocate(fac(nfac))
    f = 0
    do b = 1, nblk
       do ax = 1, 3
          do sd = 1, 2
             f = f + 1
             call make_face(blk(b), b, ax, sd, tol, fac(f))
             call detect_pole(fac(f), tol)
          end do
       end do
    end do

    !------------------------------------- search all ordered pairs ---
    allocate(pat(64))
    npat = 0
    do fa = 1, nfac
       if (fac(fa)%pole) cycle                 ! degenerate faces are poles
       do fb = 1, nfac
          if (fa == fb) cycle
          if (fac(fb)%pole) cycle
          if (fac(fa)%blk == fac(fb)%blk .and.                          &
              fac(fa)%axis == fac(fb)%axis .and.                        &
              fac(fa)%side == fac(fb)%side) cycle
          call match_faces(fac(fa), fac(fb), tol, verbose, pat, npat)
       end do
    end do

    !------------------------------------- pole records (type 93) -----
    do fa = 1, nfac
       if (.not. fac(fa)%pole) cycle
       call pole_patch(fac(fa), blk(fac(fa)%blk), pat, npat)
       if (verbose) write(*,'(A,I0,A,A,A,A)') '  block ', fac(fa)%blk,  &
            ' face ', trim(dirname(fac(fa)%dir)), ' is a POLE (degenerate) face', &
            merge(' (point)', ' (line) ', fac(fa)%pole_axis == 0)
    end do

    !----------------------------------------------- write the file ---
    allocate(nb2b(nblk));  nb2b = 0
    do f = 1, npat
       nb2b(pat(f)%blk) = nb2b(pat(f)%blk) + 1
    end do
    call write_b2b_file(outfile, nblk, tol, pat, npat, nb2b, ierr)
    if (ierr /= 0) return

    if (verbose) then
       do b = 1, nblk
          if (nb2b(b) > 0)                                              &
             write(*,'(A,I4,A,I0)') '  block', b,                       &
                  '  b2b surfaces : ', nb2b(b)
       end do
       write(*,'(A,I0,A,I0,A)') '  total           : ', npat,          &
            ' surface records (of which pole: ', count(pat(1:npat)%typ == 93), ')'
       write(*,'(A,A)')    '  written to      : ', trim(outfile)
       write(*,'(A)')      '  ------------------------------------------'
    end if
    npatches = npat

    deallocate(blk, fac, pat, nb2b)
    if (allocated(tmp)) deallocate(tmp)
  end subroutine build_bcs_b2b

  !=====================================================================
  ! Read a multi-block ASCII PLOT3D grid (whole format:
  !   nblk / (ni nj nk per block) / per block: all x, all y, all z).
  !=====================================================================
  subroutine read_plot3d(fname, blk, nblk, ierr)
    character(len=*), intent(in)               :: fname
    type(blk_t), allocatable, intent(out)      :: blk(:)
    integer, intent(out)                       :: nblk, ierr
    integer :: iu, ios, b, i, j, k

    ierr = 0
    open(newunit=iu, file=trim(fname), status='old', action='read', iostat=ios)
    if (ios /= 0) then
       write(*,*) ' ERROR (b2b): cannot open grid file ', trim(fname)
       ierr = 1;  return
    end if
    read(iu,*,iostat=ios) nblk
    if (ios /= 0 .or. nblk < 1) then
       write(*,*) ' ERROR (b2b): cannot read number of blocks from ', trim(fname)
       ierr = 2;  close(iu);  return
    end if
    allocate(blk(nblk))
    read(iu,*,iostat=ios) (blk(b)%ni, blk(b)%nj, blk(b)%nk, b=1,nblk)
    if (ios /= 0) then
       write(*,*) ' ERROR (b2b): cannot read block dimensions from ', trim(fname)
       ierr = 3;  close(iu);  return
    end if
    do b = 1, nblk
       if (min(blk(b)%ni, blk(b)%nj, blk(b)%nk) < 2) then
          write(*,*) ' ERROR (b2b): block ', b, ' has fewer than 2 nodes in a direction'
          ierr = 4;  close(iu);  return
       end if
       associate(ni=>blk(b)%ni, nj=>blk(b)%nj, nk=>blk(b)%nk)
         allocate(blk(b)%x(ni,nj,nk), blk(b)%y(ni,nj,nk), blk(b)%z(ni,nj,nk))
         read(iu,*,iostat=ios) (((blk(b)%x(i,j,k), i=1,ni), j=1,nj), k=1,nk)
         if (ios == 0) read(iu,*,iostat=ios) (((blk(b)%y(i,j,k), i=1,ni), j=1,nj), k=1,nk)
         if (ios == 0) read(iu,*,iostat=ios) (((blk(b)%z(i,j,k), i=1,ni), j=1,nj), k=1,nk)
         if (ios /= 0) then
            write(*,*) ' ERROR (b2b): failure reading coordinates of block ', b
            ierr = 5;  close(iu);  return
         end if
       end associate
    end do
    close(iu)
  end subroutine read_plot3d

  !=====================================================================
  ! Unpack one face of a block into a 2-D node sheet and build its
  ! spatial hash table (cell size tied to the tolerance).
  !=====================================================================
  subroutine make_face(bk, b, ax, sd, tol, fc)
    type(blk_t), intent(in)    :: bk
    integer,     intent(in)    :: b, ax, sd
    real(dp),    intent(in)    :: tol
    type(face_t), intent(out)  :: fc
    integer :: n(3), m1, m2, i, j, k, idx(3), lin, h

    n = [bk%ni, bk%nj, bk%nk]
    fc%blk  = b
    fc%axis = ax
    fc%side = sd
    fc%fidx = merge(1, n(ax), sd == 1)
    fc%dir  = merge(ax, -ax, sd == 1)     ! + = min face, - = max face
    select case (ax)
    case (1);  fc%t1 = 2;  fc%t2 = 3
    case (2);  fc%t1 = 1;  fc%t2 = 3
    case (3);  fc%t1 = 1;  fc%t2 = 2
    end select
    fc%n1 = n(fc%t1);  fc%n2 = n(fc%t2)

    allocate(fc%fx(fc%n1,fc%n2), fc%fy(fc%n1,fc%n2), fc%fz(fc%n1,fc%n2))
    do m2 = 1, fc%n2
       do m1 = 1, fc%n1
          idx(ax)     = fc%fidx
          idx(fc%t1)  = m1
          idx(fc%t2)  = m2
          i = idx(1);  j = idx(2);  k = idx(3)
          fc%fx(m1,m2) = bk%x(i,j,k)
          fc%fy(m1,m2) = bk%y(i,j,k)
          fc%fz(m1,m2) = bk%z(i,j,k)
       end do
    end do
    fc%bmin(1) = minval(fc%fx);  fc%bmax(1) = maxval(fc%fx)
    fc%bmin(2) = minval(fc%fy);  fc%bmax(2) = maxval(fc%fy)
    fc%bmin(3) = minval(fc%fz);  fc%bmax(3) = maxval(fc%fz)

    ! hash table: size = next power of two >= 2 * (# nodes)
    fc%hsize = 16
    do while (fc%hsize < 2*fc%n1*fc%n2)
       fc%hsize = 2*fc%hsize
    end do
    allocate(fc%hhead(fc%hsize), fc%hnext(fc%n1*fc%n2))
    fc%hhead = 0;  fc%hnext = 0
    do m2 = 1, fc%n2
       do m1 = 1, fc%n1
          lin = (m2-1)*fc%n1 + m1
          h   = hash_cell(cell_of(fc%fx(m1,m2), tol),                   &
                          cell_of(fc%fy(m1,m2), tol),                   &
                          cell_of(fc%fz(m1,m2), tol), fc%hsize)
          fc%hnext(lin) = fc%hhead(h)
          fc%hhead(h)   = lin
       end do
    end do
  end subroutine make_face

  !=====================================================================
  ! Degenerate-face test: all nodes along tangential axis t1 (for every
  ! t2) or along t2 (for every t1) coincide within tol -> pole.
  !=====================================================================
  subroutine detect_pole(fc, tol)
    type(face_t), intent(inout) :: fc
    real(dp),     intent(in)    :: tol
    logical :: deg1, deg2
    integer :: m1, m2
    deg1 = .true.                          ! nodes coincide along t1
    do m2 = 1, fc%n2
       do m1 = 2, fc%n1
          if (abs(fc%fx(m1,m2)-fc%fx(1,m2)) > tol .or.                  &
              abs(fc%fy(m1,m2)-fc%fy(1,m2)) > tol .or.                  &
              abs(fc%fz(m1,m2)-fc%fz(1,m2)) > tol) then
             deg1 = .false.;  exit
          end if
       end do
       if (.not. deg1) exit
    end do
    deg2 = .true.                          ! nodes coincide along t2
    do m1 = 1, fc%n1
       do m2 = 2, fc%n2
          if (abs(fc%fx(m1,m2)-fc%fx(m1,1)) > tol .or.                  &
              abs(fc%fy(m1,m2)-fc%fy(m1,1)) > tol .or.                  &
              abs(fc%fz(m1,m2)-fc%fz(m1,1)) > tol) then
             deg2 = .false.;  exit
          end if
       end do
       if (.not. deg2) exit
    end do
    fc%pole = (deg1 .or. deg2)
    if (deg1 .and. deg2) then
       fc%pole_axis = 0
    else if (deg1) then
       fc%pole_axis = fc%t1
    else if (deg2) then
       fc%pole_axis = fc%t2
    end if
  end subroutine detect_pole

  ! type-93 record of a pole face (full face window, donor = itself)
  subroutine pole_patch(fc, bk, pat, npat)
    type(face_t), intent(in) :: fc
    type(blk_t),  intent(in) :: bk
    type(patch_t), allocatable, intent(inout) :: pat(:)
    integer, intent(inout) :: npat
    type(patch_t) :: p
    integer :: n(3), r(2,3), a
    n = [bk%ni, bk%nj, bk%nk]
    do a = 1, 3
       r(:,a) = [1, n(a)]
    end do
    r(:,fc%axis) = fc%fidx
    p%typ = 93
    p%blk = fc%blk;  p%dir = fc%dir
    p%irng = r(:,1);  p%jrng = r(:,2);  p%krng = r(:,3)
    p%dblk = fc%blk;  p%ddir = fc%dir
    p%dirng = r(:,1);  p%djrng = r(:,2);  p%dkrng = r(:,3)
    p%map = [1, 2, 3]
    call push_patch(pat, npat, p)
  end subroutine pole_patch

  function dirname(d) result(s)
    integer, intent(in) :: d
    character(len=5) :: s
    select case (d)
    case ( 1);  s = 'i-min'
    case (-1);  s = 'i-max'
    case ( 2);  s = 'j-min'
    case (-2);  s = 'j-max'
    case ( 3);  s = 'k-min'
    case default;  s = 'k-max'
    end select
  end function dirname

  !--- hash helpers ----------------------------------------------------
  pure function cell_of(c, tol) result(ic)
    real(dp), intent(in) :: c, tol
    integer(kind=8)      :: ic
    ic = floor(c / max(4.0_dp*tol, tiny(1.0_dp)*1.0e10_dp), kind=8)
  end function cell_of

  pure function hash_cell(c1, c2, c3, hsize) result(h)
    integer(kind=8), intent(in) :: c1, c2, c3
    integer, intent(in)         :: hsize
    integer                     :: h
    integer(kind=8)             :: k
    k = ieor(c1*73856093_8, ieor(c2*19349663_8, c3*83492791_8))
    h = int(iand(k, int(hsize-1, kind=8))) + 1
  end function hash_cell

  !=====================================================================
  ! Find the node of face fc that coincides with (px,py,pz) within tol
  ! (per coordinate).  Returns tangential indices (m1,m2), or (0,0).
  ! All hash cells that a tol-ball around the point can touch are
  ! scanned, so no match is lost at cell boundaries.  If several nodes
  ! qualify (degenerate face), the closest one is returned.
  !=====================================================================
  subroutine find_node(fc, tol, px, py, pz, m1, m2)
    type(face_t), intent(in)  :: fc
    real(dp),     intent(in)  :: tol, px, py, pz
    integer,      intent(out) :: m1, m2
    integer(kind=8) :: c1a, c1b, c2a, c2b, c3a, c3b, c1, c2, c3
    integer  :: h, lin, i1, i2
    real(dp) :: d, dbest

    m1 = 0;  m2 = 0;  dbest = huge(1.0_dp)
    c1a = cell_of(px-tol, tol);  c1b = cell_of(px+tol, tol)
    c2a = cell_of(py-tol, tol);  c2b = cell_of(py+tol, tol)
    c3a = cell_of(pz-tol, tol);  c3b = cell_of(pz+tol, tol)
    do c3 = c3a, c3b
     do c2 = c2a, c2b
      do c1 = c1a, c1b
         h   = hash_cell(c1, c2, c3, fc%hsize)
         lin = fc%hhead(h)
         do while (lin /= 0)
            i2 = (lin-1)/fc%n1 + 1
            i1 = lin - (i2-1)*fc%n1
            if (abs(fc%fx(i1,i2)-px) <= tol .and.                       &
                abs(fc%fy(i1,i2)-py) <= tol .and.                       &
                abs(fc%fz(i1,i2)-pz) <= tol) then
               d = max(abs(fc%fx(i1,i2)-px), abs(fc%fy(i1,i2)-py),      &
                       abs(fc%fz(i1,i2)-pz))
               if (d < dbest) then
                  dbest = d;  m1 = i1;  m2 = i2
               end if
            end if
            lin = fc%hnext(lin)
         end do
      end do
     end do
    end do
  end subroutine find_node

  !=====================================================================
  ! Test receiver face fa against donor face fb.  If a >=2x2 node
  ! rectangle of fa coincides with nodes of fb (within tol) through a
  ! consistent signed-permutation index map, append a patch record.
  !=====================================================================
  subroutine match_faces(fa, fb, tol, verbose, pat, npat)
    type(face_t), intent(in)      :: fa, fb
    real(dp),     intent(in)      :: tol
    logical,      intent(in)      :: verbose
    type(patch_t), allocatable, intent(inout) :: pat(:)
    integer,       intent(inout)  :: npat

    integer, allocatable :: c1(:,:), c2(:,:), lab(:,:), st(:,:)
    integer :: a1, a2, nm, r1s, r1e, r2s, r2e
    integer :: ncomp, ic, nst, q1, q2, dd, n1, n2, nm2, a1c, a2c, nw, nsplit, e1, e2
    integer, allocatable :: win(:,:)
    logical, allocatable :: cm(:,:)
    integer :: o(2), d1(2), d2(2), pr(2)
    integer :: rrng(2,3), drng(2,3), map(3)
    integer :: in_d, out_r
    type(patch_t) :: p
    logical :: ok

    !--- 1. bounding boxes must overlap (inflated by tol) --------------
    if (any(fa%bmin > fb%bmax + tol) .or. any(fb%bmin > fa%bmax + tol)) return

    !--- 2. node correspondence receiver -> donor ----------------------
    allocate(c1(fa%n1,fa%n2), c2(fa%n1,fa%n2))
    c1 = 0;  c2 = 0
    nm = 0
    do a2 = 1, fa%n2
       do a1 = 1, fa%n1
          call find_node(fb, tol, fa%fx(a1,a2), fa%fy(a1,a2), fa%fz(a1,a2), &
                         c1(a1,a2), c2(a1,a2))
          if (c1(a1,a2) > 0) nm = nm + 1
       end do
    end do
    if (nm < 4) then
       ! near miss?  Faces whose nodes coincide only within 1000 x tol are
       ! NOT connected (a resampled interface, a slightly different grid
       ! generation on the two sides): say so, because the faces then fall
       ! back to overset boundaries with fringe cells along the interface
       nm2 = 0
       do a2 = 1, fa%n2
          do a1 = 1, fa%n1
             call find_node(fb, 1000.0_dp*tol, fa%fx(a1,a2), fa%fy(a1,a2), fa%fz(a1,a2), q1, q2)
             if (q1 > 0) nm2 = nm2 + 1
          end do
       end do
       if (nm2 >= 4) write(*,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,ES9.2,A)') '  WARNING (b2b): block ', fa%blk, &
          ' face dir ', fa%dir, ' nearly coincides with block ', fb%blk, ' face dir ', fb%dir, ' (', nm2, &
          ' nodes within 1000 x b2b_tol = ', 1000.0_dp*tol, ' m, none within b2b_tol): NOT connected ->'// &
          ' overset boundary; regenerate the interface or raise b2b_tol'
       deallocate(c1, c2);  return
    end if

    !--- 3. split the matched nodes into 4-connected components; each
    !       component must fill a full index rectangle >= 2 x 2 (a face
    !       may touch the same neighbour in several disjoint windows,
    !       e.g. on both sides of a body) ------------------------------
    allocate(lab(fa%n1,fa%n2), st(2, fa%n1*fa%n2))
    lab = 0
    ncomp = 0
    do a2 = 1, fa%n2
       do a1 = 1, fa%n1
          if (c1(a1,a2) <= 0 .or. lab(a1,a2) /= 0) cycle
          ncomp = ncomp + 1
          nst = 1;  st(:,1) = [a1, a2];  lab(a1,a2) = ncomp
          do while (nst > 0)
             q1 = st(1,nst);  q2 = st(2,nst);  nst = nst - 1
             do dd = 1, 4
                n1 = q1 + merge(1, merge(-1, 0, dd == 2), dd == 1)
                n2 = q2 + merge(1, merge(-1, 0, dd == 4), dd == 3)
                if (n1 < 1 .or. n1 > fa%n1 .or. n2 < 1 .or. n2 > fa%n2) cycle
                if (c1(n1,n2) <= 0 .or. lab(n1,n2) /= 0) cycle
                lab(n1,n2) = ncomp
                nst = nst + 1;  st(:,nst) = [n1, n2]
             end do
          end do
       end do
    end do

    !--- 3b. windows: a rectangular component is one window; a non-
    !        rectangular one (L-shaped overlap around a wall patch that
    !        lies on the same face, as in the source grid's split 92
    !        windows) is decomposed into rectangles of matched FACE CELLS
    !        (a cell counts when its four corner nodes are matched), each
    !        becoming its own window --------------------------------------
    allocate(win(4, max(fa%n1*fa%n2, 1)));  nw = 0
    do ic = 1, ncomp
       r1s = fa%n1;  r1e = 1;  r2s = fa%n2;  r2e = 1;  nm = 0
       do a2 = 1, fa%n2
          do a1 = 1, fa%n1
             if (lab(a1,a2) == ic) then
                nm = nm + 1
                r1s = min(r1s,a1);  r1e = max(r1e,a1)
                r2s = min(r2s,a2);  r2e = max(r2e,a2)
             end if
          end do
       end do
       if (r1e-r1s < 1 .or. r2e-r2s < 1) cycle       ! only an edge or a corner
       if (nm == (r1e-r1s+1)*(r2e-r2s+1)) then
          nw = nw + 1;  win(:,nw) = [r1s, r1e, r2s, r2e]
          cycle
       end if
       allocate(cm(r1s:r1e-1, r2s:r2e-1))
       do a2 = r2s, r2e-1
          do a1 = r1s, r1e-1
             cm(a1,a2) = lab(a1,a2) == ic .and. lab(a1+1,a2) == ic .and. lab(a1,a2+1) == ic .and. lab(a1+1,a2+1) == ic
          end do
       end do
       nsplit = 0
       do a2 = r2s, r2e-1
          do a1 = r1s, r1e-1
             if (.not. cm(a1,a2)) cycle
             e1 = a1
             do while (e1 + 1 <= r1e-1)
                if (.not. cm(e1+1,a2)) exit
                e1 = e1 + 1
             end do
             e2 = a2
             do while (e2 + 1 <= r2e-1)
                if (.not. all(cm(a1:e1,e2+1))) exit
                e2 = e2 + 1
             end do
             cm(a1:e1, a2:e2) = .false.
             nw = nw + 1;  win(:,nw) = [a1, e1+1, a2, e2+1];  nsplit = nsplit + 1
          end do
       end do
       deallocate(cm)
       if (verbose) write(*,'(A,I0,A,I0,A,I0,A,I0,A,I0,A)')                  &
          '  b2b: block ', fa%blk, ' face dir ', fa%dir, ' / block ', fb%blk, &
          ' face dir ', fb%dir, ': non-rectangular overlap split into ', nsplit, ' windows'
    end do

    do ic = 1, nw
       r1s = win(1,ic);  r1e = win(2,ic);  r2s = win(3,ic);  r2e = win(4,ic)

       !--- 4. index map must be a signed permutation -------------------
       ! derived from an interior node of the window (a corner node can
       ! coincide with two donor nodes where the donor face folds on
       ! itself, e.g. the wake cut of a C-grid); a node whose matched
       ! index differs from the map is still accepted when the donor node
       ! predicted by the map lies at the same position (ambiguous node)
       a1c = min(r1s+1, r1e-1);  a2c = min(r2s+1, r2e-1)
       o  = [c1(a1c,a2c), c2(a1c,a2c)]
       d1 = [c1(a1c+1,a2c), c2(a1c+1,a2c)] - o
       d2 = [c1(a1c,a2c+1), c2(a1c,a2c+1)] - o
       o  = o - (a1c-r1s)*d1 - (a2c-r2s)*d2
       ok = (abs(d1(1))+abs(d1(2)) == 1) .and. (abs(d2(1))+abs(d2(2)) == 1) &
            .and. (abs(d1(1)*d2(2) - d1(2)*d2(1)) == 1)
       if (ok) then
          chk: do a2 = r2s, r2e
             do a1 = r1s, r1e
                pr = o + (a1-r1s)*d1 + (a2-r2s)*d2
                if (c1(a1,a2) /= pr(1) .or. c2(a1,a2) /= pr(2)) then
                   if (pr(1) >= 1 .and. pr(1) <= fb%n1 .and. pr(2) >= 1 .and. pr(2) <= fb%n2) then
                      if ((fa%fx(a1,a2)-fb%fx(pr(1),pr(2)))**2 + (fa%fy(a1,a2)-fb%fy(pr(1),pr(2)))**2 + &
                          (fa%fz(a1,a2)-fb%fz(pr(1),pr(2)))**2 <= tol*tol) then
                         c1(a1,a2) = pr(1);  c2(a1,a2) = pr(2);  cycle
                      end if
                   end if
                   ok = .false.;  exit chk
                end if
             end do
          end do chk
       end if
       if (.not. ok) then
          if (verbose) write(*,'(A,I0,A,I0,A,I0,A,I0,A)')                  &
             '  WARNING (b2b): block ', fa%blk, ' face dir ', fa%dir,      &
             ' overlaps block ', fb%blk, ' face dir ', fb%dir,             &
             ' with an inconsistent index mapping -> skipped'
          cycle
       end if

       !--- 5. assemble the patch record --------------------------------
       ! receiver node ranges (i,j,k order)
       rrng(:, fa%axis) = [fa%fidx, fa%fidx]
       rrng(:, fa%t1)   = [r1s, r1e]
       rrng(:, fa%t2)   = [r2s, r2e]

       ! index map, tangential axes
       map = 0
       if (d1(1) /= 0) then          ! receiver t1 runs along donor t1
          map(fa%t1) = sign(fb%t1, d1(1))
          map(fa%t2) = sign(fb%t2, d2(2))
       else                          ! receiver t1 runs along donor t2
          map(fa%t1) = sign(fb%t2, d1(2))
          map(fa%t2) = sign(fb%t1, d2(1))
       end if
       ! index map, normal axis: stepping into the receiver interior
       ! corresponds to stepping into the donor interior from its face
       out_r = merge(-1, 1, fa%side == 1)      ! receiver outward index sense
       in_d  = merge( 1,-1, fb%side == 1)      ! donor inward index sense
       map(fa%axis) = sign(fb%axis, out_r*in_d)

       ! donor node ranges (donor i,j,k order)
       drng(:, fb%axis) = [fb%fidx, fb%fidx]
       if (d1(1) /= 0) then
          drng(:, fb%t1) = [o(1), o(1) + (r1e-r1s)*d1(1)]
          drng(:, fb%t2) = [o(2), o(2) + (r2e-r2s)*d2(2)]
       else
          drng(:, fb%t2) = [o(2), o(2) + (r1e-r1s)*d1(2)]
          drng(:, fb%t1) = [o(1), o(1) + (r2e-r2s)*d2(1)]
       end if

       p%blk   = fa%blk;   p%dir  = fa%dir
       p%irng  = rrng(:,1);  p%jrng = rrng(:,2);  p%krng = rrng(:,3)
       p%dblk  = fb%blk;   p%ddir = fb%dir
       p%dirng = drng(:,1);  p%djrng = drng(:,2);  p%dkrng = drng(:,3)
       p%map   = map

       call push_patch(pat, npat, p)
    end do
    deallocate(c1, c2, lab, st, win)
  end subroutine match_faces

  !--- grow-on-demand storage for patches -------------------------------
  subroutine push_patch(pat, npat, p)
    type(patch_t), allocatable, intent(inout) :: pat(:)
    integer,       intent(inout)              :: npat
    type(patch_t), intent(in)                 :: p
    type(patch_t), allocatable :: tmp(:)
    if (npat == size(pat)) then
       allocate(tmp(2*size(pat)))
       tmp(1:npat) = pat(1:npat)
       call move_alloc(tmp, pat)
    end if
    npat = npat + 1
    pat(npat) = p
  end subroutine push_patch

  !=====================================================================
  ! Write 'casename.bcs_b2b' (records grouped by receiver block).
  !=====================================================================
  subroutine write_b2b_file(outfile, nblk, tol, pat, npat, nb2b, ierr)
    character(len=*), intent(in) :: outfile
    integer,          intent(in) :: nblk, npat, nb2b(:)
    real(dp),         intent(in) :: tol
    type(patch_t),    intent(in) :: pat(:)
    integer,          intent(out):: ierr
    integer :: iu, ios, b, n
    character(len=16) :: axname(-3:3)

    axname(-3) = 'k-max';  axname(-2) = 'j-max';  axname(-1) = 'i-max'
    axname( 1) = 'i-min';  axname( 2) = 'j-min';  axname( 3) = 'k-min'
    axname( 0) = '?'

    ierr = 0
    open(newunit=iu, file=trim(outfile), status='replace', action='write', iostat=ios)
    if (ios /= 0) then
       write(*,*) ' ERROR (b2b): cannot open output file ', trim(outfile)
       ierr = 6;  return
    end if

    write(iu,'(I8,16X,A)')      nblk, '! number of grid blocks'
    write(iu,'(ES16.8,8X,A)')   tol,  '! block-to-block matching tolerance'
    do b = 1, nblk
       write(iu,'(I8,16X,A,I0)') nb2b(b),                               &
            '! number of block-to-block surfaces, block ', b
       do n = 1, npat
          if (pat(n)%blk /= b) cycle
          if (pat(n)%typ == 93) then
             write(iu,'(I8,16X,A,I0,A,A,A)') 93,                        &
                  '!*** pole (degenerate face): block ', pat(n)%blk,    &
                  ' ', trim(axname(pat(n)%dir)), ' (zero-area boundary)'
          else
             write(iu,'(I8,16X,A,I0,A,A,A,I0,A,A)') 92,                 &
                  '!*** patched: block ', pat(n)%blk, ' ',              &
                  trim(axname(pat(n)%dir)), ' <-> block ', pat(n)%dblk, &
                  ' ', trim(axname(pat(n)%ddir))
          end if
          write(iu,'(I8,16X,A)')   pat(n)%dir,                          &
               '! coordinate direction normal to surface (+:min, -:max)'
          write(iu,'(2I8,8X,A)')   pat(n)%irng, '! start, end in i-direction'
          write(iu,'(2I8,8X,A)')   pat(n)%jrng, '! start, end in j-direction'
          write(iu,'(2I8,8X,A)')   pat(n)%krng, '! start, end in k-direction'
          write(iu,'(I8,16X,A)')   0, '! flag for order of first derivative bc'
          write(iu,'(I8,16X,A)')   pat(n)%dblk, '! donor block'
          write(iu,'(I8,16X,A)')   pat(n)%ddir,                         &
               '! donor coordinate direction normal to surface'
          write(iu,'(2I8,8X,A)')   pat(n)%dirng, '! donor start, end in i-direction'
          write(iu,'(2I8,8X,A)')   pat(n)%djrng, '! donor start, end in j-direction'
          write(iu,'(2I8,8X,A)')   pat(n)%dkrng, '! donor start, end in k-direction'
          write(iu,'(3I8,A)')      pat(n)%map,                          &
               '        ! index map: signed donor axis for receiver i, j, k'
       end do
    end do
    close(iu)
  end subroutine write_b2b_file

  !=====================================================================
  ! Read 'casename.bcs_b2b' back into the module arrays (this is what
  ! the multi-block solver will use at start-up in the next steps).
  ! Comment lines (starting with '!') and blank lines are skipped;
  ! trailing comments after the numbers are ignored by list reads.
  !=====================================================================
  subroutine b2b_read_file(fname, nblk_out, tol_out, ierr)
    character(len=*), intent(in)  :: fname
    integer,          intent(out) :: nblk_out, ierr
    real(dp),         intent(out) :: tol_out
    integer :: iu, ios, b, n, m, nb, ntot, ityp, mp, npl
    integer :: rd(3), rr(2,3), dd(3), dr(2,3), mp3(3)
    character(len=256) :: line

    ierr = 0
    open(newunit=iu, file=trim(fname), status='old', action='read', iostat=ios)
    if (ios /= 0) then
       write(*,*) ' ERROR (b2b): cannot open ', trim(fname)
       ierr = 7;  return
    end if
    call nxt(iu, line, ierr);  if (ierr /= 0) return
    read(line,*) nblk_out
    call nxt(iu, line, ierr);  if (ierr /= 0) return
    read(line,*) tol_out

    ! first pass: count records per type (92 = b2b, 93 = pole)
    ntot = 0;  npl = 0
    do b = 1, nblk_out
       call nxt(iu, line, ierr);  if (ierr /= 0) return
       read(line,*) nb
       do n = 1, nb
          call nxt(iu, line, ierr);  if (ierr /= 0) return
          read(line,*) ityp
          if (ityp == 93) then
             npl = npl + 1
          else
             ntot = ntot + 1
          end if
          do m = 2, 12
             call nxt(iu, line, ierr);  if (ierr /= 0) return
          end do
       end do
    end do
    rewind(iu)

    n_b2b_tot = ntot;  n_pole = npl
    if (allocated(b2b_blk)) deallocate(b2b_blk, b2b_dir, b2b_irng,       &
         b2b_jrng, b2b_krng, b2b_dblk, b2b_ddir, b2b_dirng, b2b_djrng,   &
         b2b_dkrng, b2b_map)
    if (allocated(pole_blk)) deallocate(pole_blk, pole_dir, pole_rng)
    allocate(b2b_blk(ntot), b2b_dir(ntot), b2b_irng(2,ntot),             &
             b2b_jrng(2,ntot), b2b_krng(2,ntot), b2b_dblk(ntot),         &
             b2b_ddir(ntot), b2b_dirng(2,ntot), b2b_djrng(2,ntot),       &
             b2b_dkrng(2,ntot), b2b_map(3,ntot))
    allocate(pole_blk(npl), pole_dir(npl), pole_rng(2,3,npl))

    ! second pass: read records
    call nxt(iu, line, ierr);  call nxt(iu, line, ierr)   ! nblk, tol
    m = 0;  mp = 0
    do b = 1, nblk_out
       call nxt(iu, line, ierr);  if (ierr /= 0) return
       read(line,*) nb
       do n = 1, nb
          call nxt(iu, line, ierr);  read(line,*) ityp          ! 92 / 93
          call nxt(iu, line, ierr);  read(line,*) rd(1)
          call nxt(iu, line, ierr);  read(line,*) rr(:,1)
          call nxt(iu, line, ierr);  read(line,*) rr(:,2)
          call nxt(iu, line, ierr);  read(line,*) rr(:,3)
          call nxt(iu, line, ierr)                              ! deriv flag
          call nxt(iu, line, ierr);  read(line,*) dd(1)
          call nxt(iu, line, ierr);  read(line,*) dd(2)
          call nxt(iu, line, ierr);  read(line,*) dr(:,1)
          call nxt(iu, line, ierr);  read(line,*) dr(:,2)
          call nxt(iu, line, ierr);  read(line,*) dr(:,3)
          call nxt(iu, line, ierr);  read(line,*) mp3
          if (ierr /= 0) return
          if (ityp == 93) then
             mp = mp + 1
             pole_blk(mp) = b;  pole_dir(mp) = rd(1);  pole_rng(:,:,mp) = rr
          else
             m = m + 1
             b2b_blk(m) = b;  b2b_dir(m) = rd(1)
             b2b_irng(:,m) = rr(:,1);  b2b_jrng(:,m) = rr(:,2);  b2b_krng(:,m) = rr(:,3)
             b2b_dblk(m) = dd(1);  b2b_ddir(m) = dd(2)
             b2b_dirng(:,m) = dr(:,1);  b2b_djrng(:,m) = dr(:,2);  b2b_dkrng(:,m) = dr(:,3)
             b2b_map(:,m) = mp3
          end if
       end do
    end do
    close(iu)

  contains
    subroutine nxt(unit, out, jerr)
      integer, intent(in)           :: unit
      character(len=*), intent(out) :: out
      integer, intent(out)          :: jerr
      integer :: io
      jerr = 0
      do
         read(unit,'(A)', iostat=io) out
         if (io /= 0) then
            write(*,*) ' ERROR (b2b): unexpected end of file ', trim(fname)
            jerr = 8;  return
         end if
         out = adjustl(out)
         if (len_trim(out) == 0) cycle
         if (out(1:1) == '!') cycle
         exit
      end do
    end subroutine nxt
  end subroutine b2b_read_file

end module b2b_mod

!=======================================================================
! Wrapper called from Main: rank 0 performs the search and writes
! 'casename.bcs_b2b'; the other ranks wait at the barrier.  This runs
! every start-up, before the grid is used by the solver, as the first
! step of the multi-block extension.
!=======================================================================
subroutine make_bcs_b2b_file
  use global_vars, only: casename, b2b_tol, grid_scale
  use mpi_vars
  use b2b_mod,     only: build_bcs_b2b
  implicit none
  integer :: npatches, ierr, ierr_mpi

  if (myrank == 0) then
     ! (the search reads the raw .grd coordinates: the metre tolerance is
     !  converted to the file units when Grid_scale /= 1)
     call build_bcs_b2b(trim(casename)//'.grd',                          &
                        trim(casename)//'.bcs_b2b',                      &
                        b2b_tol/grid_scale, .true., npatches, ierr)
     if (ierr /= 0) call MPI_Abort(MPI_COMM_WORLD, 19, ierr_mpi)
  end if
  call MPI_Barrier(MPI_COMM_WORLD, ierr_mpi)
end subroutine make_bcs_b2b_file
