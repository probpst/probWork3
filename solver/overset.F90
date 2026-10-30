!=======================================================================
! overset.F90
!   OVERSET (Chimera) GRID ASSEMBLY and runtime fringe exchange for
!   probCFD.
!
!   Block faces that are neither listed in <case>.bcs nor detected as
!   block-to-block / pole faces in <case>.bcs_b2b are OVERSET faces:
!   the block ends inside another block's field there.  When at least
!   one overset face exists (or SW_Ovst = 1) the assembly runs once at
!   start-up (replicated on every rank - deterministic):
!
!   1. Coverage: per block face a cell mask of the covered part (.bcs
!      patches, .bcs_b2b interfaces, periodic pairs, poles); uncovered
!      face cells are overset boundary cells.  The nf (= Ovst_fringe)
!      cell layers behind them can never be solved ("forced").
!   2. Solid bodies: every .bcs patch with a Dirichlet velocity of zero
!      magnitude (u = v = w = 0) is a NO-SLIP WALL.  Wall patches are
!      grouped into bodies through their boundary edges: an edge node
!      is "connected" when it lies within Surfovl_tol of another wall
!      patch (surface overlap, hydrofoil/KCS type), or on a block-to-
!      block interface of its block (the wall continues in the
!      neighbour block / O-grid cut), or within Surfovl_tol of a
!      physical domain-boundary patch (symmetry plane, free surface,
!      background grid end: NACA66mod / NSWC15E type).  A body whose
!      patches have no unconnected edge node is CLOSED.
!   3. Hole cutting: for every closed body all cells of all blocks whose
!      centre or any node lies inside the body (closest-point pseudo-
!      normal test on the triangulated walls, tolerant to the overlap
!      of the patches) are HOLE cells.
!   4. Cell classification ("implicit hole cutting", Lee & Baeder):
!      every cell centre is located in every other block; the finest
!      (smallest-volume) containing non-hole cell is its donor.  A cell
!      is DOMINATED when the donor is finer than itself (relative
!      tolerance Ovst_vol_tol; ties keep both active, hence symmetric
!      grids give symmetric results).  ACTIVE = neither dominated nor
!      forced nor hole.  FRINGE = non-active, non-hole cells within a
!      Chebyshev index distance nf of an active cell (they carry the
!      interpolated boundary data of the active region; nf = 2 covers
!      the 2nd-order upwind / QUICK / ENO stencils).  Rules iterated to
!      a fixed point: a dominated fringe whose donor stencil contains a
!      forced cell of the other block is made active (no fringe-from-
!      boundary-fringe), donor-stencil members of fringe cells are never
!      holes (they become fringes), active cells with a hole within nf
!      become fringes (or are reported as orphans).  Everything else
!      is a HOLE (blanked).
!   5. Faces: a face with an active neighbour cell is SOLVED (momentum
!      + projection), a face without an active but with a fringe
!      neighbour is a FRINGE face: the Cartesian velocity vector at its
!      centre is interpolated from the donor block's staggered face
!      lattices (u from the i-faces, v from the j-faces, w from the
!      k-faces, trilinear in the donor's computational space) and
!      projected onto the transformed component with the local metrics.
!      Fringe cells receive p (and alpha) from the donor cell-centre
!      lattice.  In the pressure solver the fringe rows are part of the
!      operator (p_f - I p = 0, implicit coupling), see poisson.F90.
!
!   Output: <case>.ovst (report), iblank (1 active, 0 fringe, -1 hole
!   at the nodes, CFDShip-Iowa convention) as the last variable of every
!   .q file.
!=======================================================================
module overset
  use mpi
  use global_vars
  use mpi_vars
  use ovs_geom
  implicit none
  private
  public :: ovst_setup, ovs_update_cells, ovs_update_faces, ovs_node_iblank, ovs_overlap_weight, &
            ovs_finish_status, nent, e_kind, e_rd, e_ridx, ovs_val, dcel, ovs_dump_cell, ovs_selftest, &
            ovst_reassemble, ovs_wall_weights, ovs_n_assembly, ovs_counts, ovs_ph, ovs_build_phyd

  !--- block-level status arrays (replicated) ---------------------------
  type :: iarr_t
     integer, allocatable :: a(:,:,:)
  end type iarr_t
  type :: larr_t
     logical, allocatable :: a(:,:)
  end type larr_t
  type :: rarr_t
     real(dp), allocatable :: a(:,:,:)
  end type rarr_t
  type :: ilist_t
     integer, allocatable :: e(:)
     integer :: n = 0
  end type ilist_t

  type(iarr_t), allocatable :: cs(:)          ! final cell status (nx,ny,nz)
  type(iarr_t), allocatable :: gcs(:)         ! ghost-extended cell status
  type(iarr_t), allocatable :: gf(:,:)        ! (3,nblk) ghost-extended face status
  type(larr_t), allocatable :: cov(:,:,:)     ! (2,3,nblk) face coverage masks
  type(rarr_t), allocatable :: vol(:)
  type(bkt_t),  allocatable :: bk(:)

  !--- work arrays of the classification -------------------------------
  integer, parameter :: T_FREE = 0, T_FORCED = 1, T_DOM = 2, T_SOLID = 3, T_BLH = 4   ! T_BLH: hole inside another block's boundary layer
  type(iarr_t), allocatable :: typ(:), act(:), frg(:), need(:)
  type(iarr_t), allocatable :: dnb(:)          ! donor block (0 none)
  type(iarr_t), allocatable :: dfb(:)          ! 1 = fallback (nearest-neighbour) donor
  type(iarr_t), allocatable :: ovl(:)          ! 1 = active cell of the one-cell overlap (also active in the other grid)
  type(iarr_t), allocatable :: dnc(:)          ! donor cell (3, ncell, 1)
  type(rarr_t), allocatable :: dnr(:)          ! donor rst  (3, ncell, 1)

  !--- wall / body data --------------------------------------------------
  integer :: nwp = 0                           ! wall patches
  integer, allocatable :: wp_blk(:), wp_bc(:), wp_body(:), wp_v0(:), wp_n1(:), wp_n2(:)
  integer :: nbody = 0
  logical, allocatable :: body_closed(:)
  integer, allocatable :: body_np(:), body_unc(:), body_ovl(:), body_b2b(:), body_dom(:)
  integer, allocatable :: body_nsol(:)
  type(surf_t) :: wsurf, dsurf
  type(surf_t) :: hsurf          ! hole-cutting surface: kept wall triangles + mirror images across the planar
                                 ! domain-boundary patches that close half bodies (symmetry planes)

  !--- fringe entries (replicated, sorted by (sender rank, receiver rank))
  integer :: nent = 0
  integer,  allocatable :: e_kind(:), e_rd(:), e_dd(:), e_ridx(:,:), e_base(:,:,:)
  integer,  allocatable :: e_dcel(:,:)        ! (3,nent) donor cell containing the receiver centre (0 = none)
  real(dp), allocatable :: e_box(:,:)         ! (6,nent) receiver cell box in the donor cell's natural coordinates
  real(dp), allocatable :: e_w(:,:,:)
  real(dp), allocatable :: ovs_val(:,:)        ! (3, nent) received values
  real(dp), allocatable :: ovs_ph(:)           ! (nent) two-phase: hydrostatic offset of the
                                               !   pressure interpolation, p_f = I p + ovs_ph
                                               !   (= p_hyd(z_f) - I p_hyd, see ovs_build_phyd)
  integer :: npair = 0
  integer, allocatable :: p_s(:), p_r(:), p_e1(:), p_e2(:), p_nc(:), p_nfc(:)
  type(ilist_t), allocatable :: dcel(:)        ! per domain: receiver cell entries
  type(ilist_t), allocatable :: dfac(:)        ! per domain: receiver face entries

  integer :: n_orphan_cells = 0, n_orphan_faces = 0, n_ff = 0, n_ffo = 0
  integer :: n_fill_cells = 0, n_fill_faces = 0     ! orphans closed by a same-block neighbour
  integer :: n_gapfill = 0                           ! solid-bounded active cells whose only cover is forced fringes
  integer :: n_cov_self = 0, n_cov_other = 0, n_cov_none = 0, n_cov_thin = 0   ! coverage pass counts (report)
  integer :: n_perm_cen = 0, n_perm_pt = 0                    ! permeation check: active cells reaching into a body
  integer :: n_fb_cells = 0, n_fb_faces = 0, n_fb_used = 0   ! fallback (nearest-neighbour) donors
  integer :: n_resv = 0                                        ! cells in the reserved region of another wall block
  integer :: n_blh = 0                                         ! holes inside another block's boundary layer
  integer :: n_wallb = 0                                       ! active cells bounded by solid holes without another grid
  integer :: nwf = 0                                           ! fringe faces next to solid holes -> wall faces
  integer :: n_bcf = 0                                         ! fringe faces on covered block faces (set by BC / exchange)
  integer :: n_ovl = 0                                         ! dominated cells kept active for the one-cell overlap
  integer, allocatable :: wf_list(:,:)
  integer, parameter :: FB_R = 4                              ! search radius (index layers) of the fallback
  integer, allocatable :: orph_blk(:)
  integer :: n_uncov_faces = 0
  !--- dynamic overset (moving bodies, sixdof.F90) -----------------------
  logical :: rank0v = .false.                   ! rank 0 and not quiet: assembly messages
  logical :: ovs_quiet = .false.                ! suppress the assembly messages (re-assemblies)
  integer :: ovs_n_assembly = 0                 ! number of assemblies done
  ! wall-face force weights: fraction of the two triangles of every wall
  ! face cell kept by the surface-overset trimming (0, 0.5, 1); with them
  ! the wall faces of overlapping patches of one surface integrate every
  ! surface point exactly once (forces and moments, sixdof.F90)
  type :: rarr2_t
     real(dp), allocatable :: a(:,:)
  end type rarr2_t
  integer :: nwpw = 0
  integer, allocatable :: wpw_blk(:), wpw_bc(:)
  type(rarr2_t), allocatable :: wpw(:)
  ! caches of the motion-invariant parts of the assembly (rigid bodies):
  !   selfsol(b)%a : body index whose OWN surface makes the cell solid (the
  !                  blocks of the body's object move with the surface)
  !   keep0        : trimming flags of the wall triangles (patches of one
  !                  body move together)
  !   donor cache  : per cell, the finest donor among the blocks of the SAME
  !                  motion group (static blocks, or the blocks of one moving
  !                  object) - see find_cell_donors
  logical :: cache_ok = .false.
  type(iarr_t), allocatable :: selfsol(:)
  logical, allocatable :: keep0(:)
  type :: dcache_t
     integer,  allocatable :: dnb(:,:,:), dnc(:,:,:,:), fbb(:,:,:), fbc(:,:,:,:), vdc(:,:,:,:), flag(:,:,:)
     real(dp), allocatable :: dnr(:,:,:,:), vbest(:,:,:), vdom(:,:,:), vfb(:,:,:)
  end type dcache_t
  type(dcache_t), allocatable :: dch(:)

contains

  !=====================================================================
  ! Driver (all ranks; rank 0 reports).  Called after build_metrics.
  !=====================================================================
  subroutine ovst_setup
    integer :: nunc

    rank0v = (myrank == 0)
    call face_coverage(nunc)
    n_uncov_faces = nunc
    if (ovst_switch == 0) then
       ovst_on = .false.
       if (nunc > 0 .and. myrank == 0) then
          write(*,'(A,I0,A)') '  WARNING: ', nunc, ' block faces are neither in .bcs nor in .bcs_b2b'
          write(*,'(A)')      '           but SW_Ovst = 0: their ghost values are undefined!'
       end if
    else if (ovst_switch == 1) then
       ovst_on = .true.
    else
       ovst_on = (nunc > 0)
    end if
    if (.not. ovst_on) then
       if (myrank == 0) write(*,'(A)') '  overset         : off (every block face is covered by .bcs / .bcs_b2b)'
       return
    end if
    if (myrank == 0) then
       write(*,'(A)') ' '
       write(*,'(A)') ' ---------------- overset grid assembly ----------------'
       write(*,'(A,I0,A,I0,A)') '  overset faces   : ', nunc, ' uncovered block-face windows (fringe layers: ', ovst_nf, ')'
       write(*,'(A,ES12.4)') '  Surfovl_tol [m] : ', surfovl_tol
    end if
    call add_auto_patches
    call ovst_assemble(.true.)
  end subroutine ovst_setup

  !=====================================================================
  ! Re-assembly after the blocks of moving bodies were displaced
  ! (sixdof.F90): everything from the wall surfaces on is rebuilt from
  ! the current node coordinates blk(b)%xn/yn/zn; the face coverage and
  ! the automatic overset patches are topological and stay.  Quiet
  ! (one summary line); the report file is rewritten.
  !=====================================================================
  subroutine ovst_reassemble(tsec)
    real(dp), intent(out) :: tsec
    real(dp) :: t0
    if (.not. ovst_on) then
       tsec = 0.0_dp;  return
    end if
    t0 = MPI_Wtime()
    ovs_quiet = .not. poi_verbose                 ! poi_verbose = 1: full messages of every re-assembly
    rank0v = (myrank == 0 .and. .not. ovs_quiet)
    call ovst_clear
    call ovst_assemble(.false.)
    ovs_quiet = .false.
    rank0v = (myrank == 0)
    tsec = MPI_Wtime() - t0
  end subroutine ovst_reassemble

  !---------------------------------------------------------------------
  ! Release every result of a previous assembly (re-entrant assembly).
  !---------------------------------------------------------------------
  subroutine ovst_clear
    integer :: b, d, q, ib
    if (allocated(cs)) then
       do b = 1, nblk
          if (allocated(cs(b)%a)) deallocate(cs(b)%a)
       end do
       deallocate(cs)
    end if
    if (allocated(gcs)) then
       do b = 1, nblk
          if (allocated(gcs(b)%a)) deallocate(gcs(b)%a)
          do q = 1, 3
             if (allocated(gf(q,b)%a)) deallocate(gf(q,b)%a)
          end do
       end do
       deallocate(gcs, gf)
    end if
    if (allocated(ovl)) then
       do b = 1, nblk
          if (allocated(ovl(b)%a)) deallocate(ovl(b)%a)
       end do
       deallocate(ovl)
    end if
    if (allocated(wp_blk)) deallocate(wp_blk, wp_bc, wp_body, wp_v0, wp_n1, wp_n2)
    if (allocated(body_np)) deallocate(body_np, body_unc, body_ovl, body_b2b, body_dom, body_closed, body_nsol)
    call surf_free_all(wsurf);  call surf_free_all(dsurf);  call surf_free_all(hsurf)
    if (allocated(e_kind)) deallocate(e_kind, e_rd, e_dd, e_ridx, e_base, e_w, ovs_val)
    if (allocated(ovs_ph)) deallocate(ovs_ph)
    if (allocated(e_dcel)) deallocate(e_dcel, e_box)
    if (allocated(p_s)) deallocate(p_s, p_r, p_e1, p_e2, p_nc, p_nfc)
    if (allocated(dcel)) then
       do d = 1, ndom
          if (allocated(dcel(d)%e)) deallocate(dcel(d)%e)
          if (allocated(dfac(d)%e)) deallocate(dfac(d)%e)
       end do
       deallocate(dcel, dfac)
    end if
    if (allocated(wf_list)) deallocate(wf_list)
    if (allocated(orph_blk)) deallocate(orph_blk)
    do ib = 1, nown
       d = own_list(ib)
       if (allocated(dom(d)%cs)) deallocate(dom(d)%cs, dom(d)%fs1, dom(d)%fs2, dom(d)%fs3, dom(d)%ovd)
    end do
    nent = 0;  npair = 0;  nwp = 0;  nbody = 0;  nwf = 0
  end subroutine ovst_clear

  !---------------------------------------------------------------------
  ! Assembly steps 2..5 from the current node coordinates.
  !---------------------------------------------------------------------
  subroutine ovst_assemble(first)
    logical, intent(in) :: first
    integer :: b, ierr
    real(dp) :: t0, t1, tstage(7)

    t0 = MPI_Wtime()
    allocate(vol(nblk), bk(nblk), cs(nblk), typ(nblk), act(nblk), frg(nblk),   &
             need(nblk), dnb(nblk), dfb(nblk), dnc(nblk), dnr(nblk))
    do b = 1, nblk
       allocate(vol(b)%a(blk(b)%nx, blk(b)%ny, blk(b)%nz))
       call cell_volumes(blk(b)%xn, blk(b)%yn, blk(b)%zn, blk(b)%ni, blk(b)%nj, blk(b)%nk, vol(b)%a)
       call build_buckets(blk(b)%xn, blk(b)%yn, blk(b)%zn, blk(b)%ni, blk(b)%nj, blk(b)%nk, bk(b))
       allocate(cs(b)%a(blk(b)%nx,blk(b)%ny,blk(b)%nz), typ(b)%a(blk(b)%nx,blk(b)%ny,blk(b)%nz))
       allocate(act(b)%a(blk(b)%nx,blk(b)%ny,blk(b)%nz), frg(b)%a(blk(b)%nx,blk(b)%ny,blk(b)%nz))
       allocate(need(b)%a(blk(b)%nx,blk(b)%ny,blk(b)%nz))
       allocate(dnb(b)%a(blk(b)%nx,blk(b)%ny,blk(b)%nz), dnc(b)%a(3, blk(b)%nx*blk(b)%ny*blk(b)%nz, 1))
       allocate(dfb(b)%a(blk(b)%nx,blk(b)%ny,blk(b)%nz));  dfb(b)%a = 0
       if (.not. allocated(ovl)) allocate(ovl(nblk))
       if (allocated(ovl(b)%a)) deallocate(ovl(b)%a)
       allocate(ovl(b)%a(blk(b)%nx,blk(b)%ny,blk(b)%nz));  ovl(b)%a = 0
       allocate(dnr(b)%a(3, blk(b)%nx*blk(b)%ny*blk(b)%nz, 1))
       typ(b)%a = T_FREE;  dnb(b)%a = 0;  cs(b)%a = 1
    end do
    allocate(orph_blk(nblk));  orph_blk = 0
    tstage(1) = MPI_Wtime()
    call wall_bodies                       ! 2. bodies + closedness
    tstage(2) = MPI_Wtime()
    call cut_holes                         ! 3. solid holes
    tstage(3) = MPI_Wtime()
    call blank_in_boundary_layers          ! 3b. holes inside another wall block's boundary layer
    tstage(4) = MPI_Wtime()
    call mark_forced                       ! 1b. cells behind overset faces
    call find_cell_donors                  ! 4a. finest donor of every cell
    tstage(5) = MPI_Wtime()
    call classify                          ! 4b. active / fringe / hole
    if (first .and. ovst_only /= 0) call check_permeation   ! grid check: no solved cell may reach into a body
    tstage(6) = MPI_Wtime()
    call ghost_status                      ! ghost-extended status, faces
    call build_entries                     ! 5. fringe cells + faces -> donors
    tstage(7) = MPI_Wtime()
    call domain_status                     ! per-domain slices
    call ovs_build_phyd                    ! 6. two-phase: hydrostatic offset of the p interpolation
    ovs_n_assembly = ovs_n_assembly + 1
    cache_ok = .true.                      ! the motion-invariant caches are filled
    t1 = MPI_Wtime()
    if (myrank == 0 .and. (first .or. poi_verbose)) write(*,'(A,7F8.2,A)') '  assembly timing : geometry / bodies / holes / BL-holes / donors / classify / entries / domains =', &
       tstage(1)-t0, tstage(2)-tstage(1), tstage(3)-tstage(2), tstage(4)-tstage(3), tstage(5)-tstage(4), tstage(6)-tstage(5), tstage(7)-tstage(6), ' s'
    if (myrank == 0 .and. (first .or. n_orphan_cells > 0)) call write_report(t1 - t0)
    if (first) call MPI_Barrier(MPI_COMM_WORLD, ierr)
    do b = 1, nblk
       deallocate(typ(b)%a, act(b)%a, frg(b)%a, need(b)%a, dnb(b)%a, dfb(b)%a, dnc(b)%a, dnr(b)%a, vol(b)%a)
       if (allocated(bk(b)%start)) deallocate(bk(b)%start, bk(b)%list)
    end do
    deallocate(typ, act, frg, need, dnb, dfb, dnc, dnr, vol, bk)
    call surf_free(wsurf);  call surf_free(dsurf)      ! hsurf is kept for the node-iblank output
  end subroutine ovst_assemble

  !---------------------------------------------------------------------
  ! Global cell counts of the current assembly (monitor of the dynamic
  ! overset): active, fringe, hole, orphan cells, fringe entries.
  !---------------------------------------------------------------------
  subroutine ovs_counts(nact, nfr, nhole, norph, nentr)
    integer, intent(out) :: nact, nfr, nhole, norph, nentr
    integer :: b
    nact = 0;  nfr = 0;  nhole = 0;  norph = n_orphan_cells;  nentr = nent
    if (.not. allocated(cs)) return
    do b = 1, nblk
       nact = nact + count(cs(b)%a == 1);  nfr = nfr + count(cs(b)%a == 0);  nhole = nhole + count(cs(b)%a < 0)
    end do
  end subroutine ovs_counts

  !---------------------------------------------------------------------
  ! Force-integration weights of the wall face cells of patch n of
  ! block b: w(m1,m2) = fraction of the face's two triangles kept by the
  ! surface-overset trimming (1.0 when the assembly is off or the patch
  ! is not a wall patch).
  !---------------------------------------------------------------------
  subroutine ovs_wall_weights(b, n, w)
    integer,  intent(in)  :: b, n
    real(dp), intent(out) :: w(:,:)
    integer :: ip
    w = 1.0_dp
    if (.not. ovst_on .or. nwpw == 0) return
    do ip = 1, nwpw
       if (wpw_blk(ip) == b .and. wpw_bc(ip) == n) then
          if (size(w,1) == size(wpw(ip)%a,1) .and. size(w,2) == size(wpw(ip)%a,2)) w = wpw(ip)%a
          return
       end if
    end do
  end subroutine ovs_wall_weights

  !=====================================================================
  ! 1. Coverage masks of the 6 faces of every block.
  !=====================================================================
  subroutine face_coverage(nunc)
    integer, intent(out) :: nunc
    integer :: b, ax, sd, n, m, t1, t2, n1, n2, c1s, c1e, c2s, c2e, rr(2,3)
    allocate(cov(2,3,nblk))
    nunc = 0
    do b = 1, nblk
       do ax = 1, 3
          call tang(ax, t1, t2)
          n1 = ncell(b,t1);  n2 = ncell(b,t2)
          do sd = 1, 2
             allocate(cov(sd,ax,b)%a(n1,n2));  cov(sd,ax,b)%a = .false.
             if (blk(b)%pole(sd,ax)) cov(sd,ax,b)%a = .true.
             do n = 1, blk(b)%nbc
                if (blk(b)%bc_dir(n) /= merge(ax, -ax, sd == 1)) cycle
                rr(:,1) = blk(b)%bc_irng(:,n);  rr(:,2) = blk(b)%bc_jrng(:,n);  rr(:,3) = blk(b)%bc_krng(:,n)
                c1s = max(1, rr(1,t1));  c1e = min(n1, rr(2,t1)-1)
                c2s = max(1, rr(1,t2));  c2e = min(n2, rr(2,t2)-1)
                if (c1s <= c1e .and. c2s <= c2e) cov(sd,ax,b)%a(c1s:c1e, c2s:c2e) = .true.
             end do
             do m = 1, nxp
                if (xp(m)%rblk /= b .or. xp(m)%na /= ax .or. xp(m)%rside /= sd) cycle
                c1s = max(1, xp(m)%rn(1,t1));  c1e = min(n1, xp(m)%rn(2,t1)-1)
                c2s = max(1, xp(m)%rn(1,t2));  c2e = min(n2, xp(m)%rn(2,t2)-1)
                if (c1s <= c1e .and. c2s <= c2e) cov(sd,ax,b)%a(c1s:c1e, c2s:c2e) = .true.
             end do
             if (.not. all(cov(sd,ax,b)%a)) then
                nunc = nunc + 1
                ! a face that is PARTLY block-to-block: the nodes of the two
                ! sides match only on part of the face (e.g. different
                ! wall-normal distributions after a grid modification) -
                ! the uncovered rest becomes an overset boundary
                do m = 1, nxp
                   if (xp(m)%rblk /= b .or. xp(m)%na /= ax .or. xp(m)%rside /= sd) cycle
                   if (myrank == 0) write(*,'(A,I0,A,I0,A,I0,A,I0,A)') '  WARNING: block ', b, ' face axis ', ax, &
                      ' side ', sd, ' is only PARTLY block-to-block (', count(.not. cov(sd,ax,b)%a), &
                      ' uncovered face cells become an overset boundary) - check the node match of the interface'
                   exit
                end do
             end if
          end do
       end do
    end do
  end subroutine face_coverage

  pure subroutine tang(ax, t1, t2)
    integer, intent(in) :: ax
    integer, intent(out) :: t1, t2
    select case (ax)
    case (1);  t1 = 2;  t2 = 3
    case (2);  t1 = 1;  t2 = 3
    case default;  t1 = 1;  t2 = 2
    end select
  end subroutine tang

  !=====================================================================
  ! Automatic zero-gradient patches (all variables) on every block face
  ! with an uncovered part, PREPENDED to the bc list so that the .bcs
  ! patches of the covered part override them.  They only keep the
  ! ghost values finite (fringe cells behind an overset face are
  ! interpolated, the solved stencils never reach those ghosts).
  !=====================================================================
  subroutine add_auto_patches
    integer :: b, ax, sd, nadd, n, nb, k
    integer, allocatable :: dir(:), ir(:,:), jr(:,:), kr(:,:), fl(:,:), ef(:,:)
    real(dp), allocatable :: vl(:,:), ti(:)
    do b = 1, nblk
       nadd = 0
       do ax = 1, 3
          do sd = 1, 2
             if (.not. all(cov(sd,ax,b)%a)) nadd = nadd + 1
          end do
       end do
       if (nadd == 0) cycle
       nb = blk(b)%nbc
       allocate(dir(nb+nadd), ir(2,nb+nadd), jr(2,nb+nadd), kr(2,nb+nadd),  &
                fl(6,nb+nadd), vl(6,nb+nadd), ef(3,nb+nadd), ti(nb+nadd))     ! 6 rows: u v w p alpha T (update 16)
       k = 0
       do ax = 1, 3
          do sd = 1, 2
             if (all(cov(sd,ax,b)%a)) cycle
             k = k + 1
             dir(k) = merge(ax, -ax, sd == 1)
             ir(:,k) = [1, blk(b)%ni];  jr(:,k) = [1, blk(b)%nj];  kr(:,k) = [1, blk(b)%nk]
             select case (ax)
             case (1);  ir(:,k) = merge(1, blk(b)%ni, sd == 1)
             case (2);  jr(:,k) = merge(1, blk(b)%nj, sd == 1)
             case (3);  kr(:,k) = merge(1, blk(b)%nk, sd == 1)
             end select
             fl(:,k) = 0;  vl(:,k) = 0.0_dp;  ef(:,k) = 0;  ti(k) = 0.0_dp
          end do
       end do
       do n = 1, nb
          dir(nadd+n) = blk(b)%bc_dir(n)
          ir(:,nadd+n) = blk(b)%bc_irng(:,n);  jr(:,nadd+n) = blk(b)%bc_jrng(:,n)
          kr(:,nadd+n) = blk(b)%bc_krng(:,n)
          fl(:,nadd+n) = blk(b)%bc_flag(:,n);  vl(:,nadd+n) = blk(b)%bc_val(:,n)
          ef(:,nadd+n) = blk(b)%bc_eflag(:,n);  ti(nadd+n) = blk(b)%bc_tinf(n)
       end do
       deallocate(blk(b)%bc_dir, blk(b)%bc_irng, blk(b)%bc_jrng, blk(b)%bc_krng, &
                  blk(b)%bc_flag, blk(b)%bc_val, blk(b)%bc_eflag, blk(b)%bc_tinf)
       call move_alloc(dir, blk(b)%bc_dir);  call move_alloc(ir, blk(b)%bc_irng)
       call move_alloc(jr, blk(b)%bc_jrng);  call move_alloc(kr, blk(b)%bc_krng)
       call move_alloc(fl, blk(b)%bc_flag);  call move_alloc(vl, blk(b)%bc_val)
       call move_alloc(ef, blk(b)%bc_eflag);  call move_alloc(ti, blk(b)%bc_tinf)
       blk(b)%nbc = nb + nadd
    end do
  end subroutine add_auto_patches

  !=====================================================================
  ! 2. Wall patches -> triangle surface, bodies, closedness.
  !=====================================================================
  subroutine wall_bodies
    integer :: b, n, ax, sd, t1, t2, n1, n2, m1, m2, idx(3), inw(3), i, j, k
    integer :: nvmax, ntmax, ip, v0, iq, ib
    real(dp), allocatable :: fx(:,:), fy(:,:), fz(:,:), fd(:,:,:)
    integer, allocatable :: root(:), bid(:), p_unc(:), p_ovl(:), p_b2b(:), p_dom(:)
    real(dp) :: p(3), dist, q(3), wb(3)
    integer  :: tb, region, cnt_ovl, cnt_b2b, cnt_dom, cnt_unc
    logical  :: conn

    ! ---- count and size ---------------------------------------------------
    nwp = 0;  nvmax = 0;  ntmax = 0
    do b = 1, nblk
       do n = 1, blk(b)%nbc
          if (is_wall(b, n)) nwp = nwp + 1
       end do
    end do
    allocate(wp_blk(max(nwp,1)), wp_bc(max(nwp,1)), wp_body(max(nwp,1)), wp_v0(max(nwp,1)), &
             wp_n1(max(nwp,1)), wp_n2(max(nwp,1)))
    do b = 1, nblk
       do n = 1, blk(b)%nbc
          call patch_dims(b, n, ax, sd, t1, t2, n1, n2)
          if (is_wall(b, n)) then
             nvmax = nvmax + n1*n2;  ntmax = ntmax + 2*(n1-1)*(n2-1)
          end if
       end do
    end do
    call surf_init(wsurf, nvmax, ntmax)
    nvmax = 0;  ntmax = 0
    do b = 1, nblk
       do n = 1, blk(b)%nbc
          call patch_dims(b, n, ax, sd, t1, t2, n1, n2)
          if (is_phys(b, n) .and. .not. is_wall(b, n)) then
             nvmax = nvmax + n1*n2;  ntmax = ntmax + 2*(n1-1)*(n2-1)
          end if
       end do
    end do
    call surf_init(dsurf, nvmax, ntmax)

    ! ---- build the sheets ---------------------------------------------
    ip = 0;  iq = 0
    do b = 1, nblk
       do n = 1, blk(b)%nbc
          if (.not. is_phys(b, n)) cycle
          call patch_dims(b, n, ax, sd, t1, t2, n1, n2)
          allocate(fx(n1,n2), fy(n1,n2), fz(n1,n2), fd(3,n1,n2))
          do m2 = 1, n2
             do m1 = 1, n1
                idx(t1) = rng1(b,n,t1) + m1 - 1
                idx(t2) = rng1(b,n,t2) + m2 - 1
                idx(ax) = merge(1, nnode(b,ax), sd == 1)
                i = idx(1);  j = idx(2);  k = idx(3)
                fx(m1,m2) = blk(b)%xn(i,j,k);  fy(m1,m2) = blk(b)%yn(i,j,k);  fz(m1,m2) = blk(b)%zn(i,j,k)
                inw = idx;  inw(ax) = inw(ax) + merge(1, -1, sd == 1)
                fd(:,m1,m2) = [blk(b)%xn(inw(1),inw(2),inw(3)) - fx(m1,m2),  &
                               blk(b)%yn(inw(1),inw(2),inw(3)) - fy(m1,m2),  &
                               blk(b)%zn(inw(1),inw(2),inw(3)) - fz(m1,m2)]
             end do
          end do
          if (is_wall(b, n)) then
             ip = ip + 1
             wp_blk(ip) = b;  wp_bc(ip) = n;  wp_n1(ip) = n1;  wp_n2(ip) = n2
             call surf_add_patch(wsurf, ip, n1, n2, fx, fy, fz, fd, v0)
             wp_v0(ip) = v0
          else
             iq = iq + 1
             call surf_add_patch(dsurf, iq, n1, n2, fx, fy, fz, fd, v0)
          end if
          deallocate(fx, fy, fz, fd)
       end do
    end do
    call surf_finish(wsurf)
    call surf_finish(dsurf)

    ! ---- connections through the boundary edge nodes -------------------
    allocate(root(max(nwp,1)), p_unc(max(nwp,1)), p_ovl(max(nwp,1)), p_b2b(max(nwp,1)), p_dom(max(nwp,1)))
    do ip = 1, nwp
       root(ip) = ip
    end do
    p_unc = 0;  p_ovl = 0;  p_b2b = 0;  p_dom = 0
    do ip = 1, nwp
       b = wp_blk(ip);  n = wp_bc(ip)
       call patch_dims(b, n, ax, sd, t1, t2, n1, n2)
       cnt_ovl = 0;  cnt_b2b = 0;  cnt_dom = 0;  cnt_unc = 0
       do m2 = 1, n2
          do m1 = 1, n1
             if (.not. (m1 == 1 .or. m1 == n1 .or. m2 == 1 .or. m2 == n2)) cycle
             p = wsurf%v(:, wp_v0(ip) + m1 + (m2-1)*n1)
             conn = .false.
             ! (a) overlap with another wall patch
             if (nwp > 1) then
                call surf_closest(wsurf, p, 1, ip, surfovl_tol*1.0001_dp, dist, q, tb, wb, region)
                if (tb > 0) then
                   if (dist <= surfovl_tol) then
                      conn = .true.;  cnt_ovl = cnt_ovl + 1
                      call unite(ip, wsurf%tp(tb))
                   end if
                end if
             end if
             ! (c) on a b2b / periodic interface of the block
             if (.not. conn) then
                idx(t1) = rng1(b,n,t1) + m1 - 1;  idx(t2) = rng1(b,n,t2) + m2 - 1
                idx(ax) = merge(1, nnode(b,ax), sd == 1)
                if (on_b2b(b, idx)) then
                   conn = .true.;  cnt_b2b = cnt_b2b + 1
                end if
             end if
             ! (d) self-closure of the patch: the edge node coincides (within
             !     Surfovl_tol) with an edge node of another edge of the SAME
             !     patch - a degenerate pole line (all nodes of the edge are one
             !     point) or an O-cut seam (the m1 = 1 and m1 = n1 lines are the
             !     same curve) - the surface closes on itself there
             if (.not. conn) then
                if (m1 == 1 .or. m1 == n1) then
                   q = wsurf%v(:, wp_v0(ip) + merge(n1, 1, m1 == 1) + (m2-1)*n1)
                   if (sqrt(sum((q - p)**2)) <= surfovl_tol) conn = .true.
                   if (.not. conn .and. m2 > 1) then
                      q = wsurf%v(:, wp_v0(ip) + m1 + (m2-2)*n1)
                      if (sqrt(sum((q - p)**2)) <= surfovl_tol) conn = .true.
                   end if
                   if (.not. conn .and. m2 < n2) then
                      q = wsurf%v(:, wp_v0(ip) + m1 + m2*n1)
                      if (sqrt(sum((q - p)**2)) <= surfovl_tol) conn = .true.
                   end if
                end if
                if (.not. conn .and. (m2 == 1 .or. m2 == n2)) then
                   q = wsurf%v(:, wp_v0(ip) + m1 + (merge(n2, 1, m2 == 1)-1)*n1)
                   if (sqrt(sum((q - p)**2)) <= surfovl_tol) conn = .true.
                   if (.not. conn .and. m1 > 1) then
                      q = wsurf%v(:, wp_v0(ip) + m1 - 1 + (m2-1)*n1)
                      if (sqrt(sum((q - p)**2)) <= surfovl_tol) conn = .true.
                   end if
                   if (.not. conn .and. m1 < n1) then
                      q = wsurf%v(:, wp_v0(ip) + m1 + 1 + (m2-1)*n1)
                      if (sqrt(sum((q - p)**2)) <= surfovl_tol) conn = .true.
                   end if
                end if
                if (conn) cnt_b2b = cnt_b2b + 1
             end if
             ! (b) on a physical domain boundary
             if (.not. conn .and. dsurf%nt > 0) then
                call surf_closest(dsurf, p, 0, 0, surfovl_tol*1.0001_dp, dist, q, tb, wb, region)
                if (tb > 0) then
                   if (dist <= surfovl_tol) then
                      conn = .true.;  cnt_dom = cnt_dom + 1
                   end if
                end if
             end if
             if (.not. conn) then
                cnt_unc = cnt_unc + 1
                if (rank0v .and. cnt_unc <= 20) then
                   call surf_closest(wsurf, p, 1, ip, huge(1.0_dp), dist, q, tb, wb, region)
                   write(*,'(A,I0,A,2I5,A,3F11.6,A,ES10.3,A,I0)') '    unconnected edge node of wall ', ip, &
                      ' (m1,m2) ', m1, m2, '  xyz ', p, '  nearest other wall at ', dist, ' m, wall ', &
                      merge(wsurf%tp(max(tb,1)), 0, tb > 0)
                end if
             end if
          end do
       end do
       p_unc(ip) = cnt_unc;  p_ovl(ip) = cnt_ovl;  p_b2b(ip) = cnt_b2b;  p_dom(ip) = cnt_dom
    end do
    ! bodies = union-find components, renumbered
    allocate(bid(max(nwp,1)));  bid = 0
    nbody = 0
    do ip = 1, nwp
       if (find(ip) == ip) then
          nbody = nbody + 1;  bid(ip) = nbody
       end if
    end do
    do ip = 1, nwp
       wp_body(ip) = bid(find(ip))
    end do
    allocate(body_np(max(nbody,1)), body_unc(max(nbody,1)), body_ovl(max(nbody,1)), &
             body_b2b(max(nbody,1)), body_dom(max(nbody,1)), body_closed(max(nbody,1)), body_nsol(max(nbody,1)))
    body_np = 0;  body_unc = 0;  body_ovl = 0;  body_b2b = 0;  body_dom = 0;  body_nsol = 0
    do ip = 1, nwp
       ib = wp_body(ip)
       body_np(ib) = body_np(ib) + 1;  body_unc(ib) = body_unc(ib) + p_unc(ip)
       body_ovl(ib) = body_ovl(ib) + p_ovl(ip);  body_b2b(ib) = body_b2b(ib) + p_b2b(ip)
       body_dom(ib) = body_dom(ib) + p_dom(ip)
    end do
    do ib = 1, nbody
       body_closed(ib) = (body_unc(ib) == 0)
    end do
    do n = 1, wsurf%nt
       wsurf%tb(n) = wp_body(wsurf%tp(n))
    end do
    deallocate(root, bid, p_unc, p_ovl, p_b2b, p_dom)
    if (rank0v) then
       write(*,'(A,I0,A,I0,A)') '  wall patches    : ', nwp, ' no-slip surfaces -> ', nbody, ' solid bodies'
       do ib = 1, nbody
          write(*,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A)') '    body ', ib, ': ', body_np(ib),   &
             ' patches; edge nodes connected by overlap / b2b / domain boundary: ',           &
             body_ovl(ib), ' / ', body_b2b(ib), ' / ', body_dom(ib), ', unconnected: ', body_unc(ib), &
             merge('  -> CLOSED', '  -> OPEN  ', body_closed(ib))
       end do
    end if

  contains
    recursive integer function find(x) result(r)
      integer, intent(in) :: x
      r = x
      do while (root(r) /= r)
         root(r) = root(root(r));  r = root(r)
      end do
    end function find
    subroutine unite(a, c)
      integer, intent(in) :: a, c
      integer :: ra, rc
      ra = find(a);  rc = find(c)
      if (ra /= rc) root(max(ra,rc)) = min(ra,rc)
    end subroutine unite
  end subroutine wall_bodies

  ! motion group of block b = its object id (0 background): the blocks
  ! of one group keep their relative position (rigid motion), so every
  ! geometric relation inside a group is invariant in time
  pure integer function mgroup(b)
    integer, intent(in) :: b
    mgroup = 0
    if (.not. allocated(blk_obj)) return
    mgroup = blk_obj(b)
  end function mgroup

  ! object id carrying overset body ib (0 background; -1 when its patches
  ! lie on blocks of different objects)
  integer function body_object(ib)
    integer, intent(in) :: ib
    integer :: ip, o
    body_object = -2
    do ip = 1, nwp
       if (wp_body(ip) /= ib) cycle
       o = 0
       if (allocated(blk_obj)) o = blk_obj(wp_blk(ip))
       if (body_object == -2) then
          body_object = o
       else if (body_object /= o) then
          body_object = -1
       end if
    end do
    if (body_object == -2) body_object = 0
  end function body_object

  ! node range start of patch n of block b along axis a
  pure integer function rng1(b, n, a)
    integer, intent(in) :: b, n, a
    select case (a)
    case (1);  rng1 = blk(b)%bc_irng(1,n)
    case (2);  rng1 = blk(b)%bc_jrng(1,n)
    case default;  rng1 = blk(b)%bc_krng(1,n)
    end select
  end function rng1

  pure integer function nnode(b, a)
    integer, intent(in) :: b, a
    select case (a)
    case (1);  nnode = blk(b)%ni
    case (2);  nnode = blk(b)%nj
    case default;  nnode = blk(b)%nk
    end select
  end function nnode

  pure subroutine patch_dims(b, n, ax, sd, t1, t2, n1, n2)
    integer, intent(in)  :: b, n
    integer, intent(out) :: ax, sd, t1, t2, n1, n2
    integer :: rr(2,3)
    ax = abs(blk(b)%bc_dir(n));  sd = merge(1, 2, blk(b)%bc_dir(n) > 0)
    call tang(ax, t1, t2)
    rr(:,1) = blk(b)%bc_irng(:,n);  rr(:,2) = blk(b)%bc_jrng(:,n);  rr(:,3) = blk(b)%bc_krng(:,n)
    n1 = rr(2,t1) - rr(1,t1) + 1;  n2 = rr(2,t2) - rr(1,t2) + 1
  end subroutine patch_dims

  ! no-slip wall: Dirichlet velocity of zero magnitude
  pure logical function is_wall(b, n)
    integer, intent(in) :: b, n
    is_wall = (all(blk(b)%bc_flag(1:3,n) == 1) .and. all(blk(b)%bc_val(1:3,n) == 0.0_dp)) .or. all(blk(b)%bc_flag(1:3,n) == 11)
  end function is_wall

  ! genuine .bcs patch: not a pole face and not an automatic overset
  ! patch (all-zero flags on a face with an uncovered part)
  logical function is_phys(b, n)
    integer, intent(in) :: b, n
    integer :: ax, sd
    ax = abs(blk(b)%bc_dir(n));  sd = merge(1, 2, blk(b)%bc_dir(n) > 0)
    is_phys = .true.
    if (blk(b)%pole(sd,ax)) is_phys = .false.
    if (all(blk(b)%bc_flag(:,n) == 0) .and. .not. all(cov(sd,ax,b)%a)) is_phys = .false.
  end function is_phys

  ! node idx of block b lies on a face window covered by a b2b/periodic patch
  logical function on_b2b(b, idx)
    integer, intent(in) :: b, idx(3)
    integer :: a, sd, t1, t2, c1, c2, m
    on_b2b = .false.
    do a = 1, 3
       do sd = 1, 2
          if (idx(a) /= merge(1, nnode(b,a), sd == 1)) cycle
          call tang(a, t1, t2)
          do m = 1, nxp
             if (xp(m)%rblk /= b .or. xp(m)%na /= a .or. xp(m)%rside /= sd) cycle
             c1 = idx(t1);  c2 = idx(t2)
             if (c1 >= xp(m)%rn(1,t1) .and. c1 <= xp(m)%rn(2,t1) .and.       &
                 c2 >= xp(m)%rn(1,t2) .and. c2 <= xp(m)%rn(2,t2)) then
                on_b2b = .true.;  return
             end if
          end do
       end do
    end do
  end function on_b2b

  !=====================================================================
  ! 3. Hole cutting: cells with the centre or a node inside a closed
  !    body.
  !=====================================================================
  subroutine cut_holes
    ! Solid hole cutting, update 6:
    !  (a) overlap trimming: a triangle of patch P whose centroid lies
    !      within Surfovl_tol of a FINER triangle of another patch Q of the
    !      same body, at a point interior to Q (farther than Surfovl_tol
    !      from Q's boundary nodes), is dropped from the hole-cutting
    !      surface (surface overset: one sheet per surface point);
    !  (b) inside test by the generalized winding number w of the trimmed
    !      body surface (cell centre or any node with w > 0.5 -> hole);
    !      no distance margin: Surfovl_tol only enters the closure and
    !      trimming tests;
    !  (c) flood fill of the remaining cells of every block from its non-
    !      wall boundary faces: a region that is cut off from them by
    !      hole cells and the block's own wall (an appendage root, a duct
    !      penetration) is inside the body as well.
    integer :: ib, b, i, j, k, t, nn, nsol, nflood, ntrim, ip, iq, lin, tb, region
    real(dp) :: bmin(3), bmax(3), p(3), marg, w, wmax, wmin, dist, q(3), wbc(3), cen(3), diag, twind, rq(3)
    integer  :: ncross
    logical  :: okq
    logical  :: inside
    integer(1), allocatable :: lab(:,:,:)
    integer, allocatable :: queue(:,:)
    integer :: nq, iq0, ii, jj, kk, ni_, nj_, nk_, nb_dir, dd(3), ax, sd, t1, t2, n1, n2, v0, m1, m2
    integer :: kobj
    logical :: selfc
    real(dp) :: dq
    ! ---- (a) trim covered triangles --------------------------------------
    allocate(wsurf%keep(wsurf%nt));  wsurf%keep = .true.;  ntrim = 0
    if (cache_ok .and. allocated(keep0)) then
       if (size(keep0) == wsurf%nt) then
          wsurf%keep = keep0;  ntrim = count(.not. keep0)
       end if
    end if
    do t = 1, wsurf%nt
       if (cache_ok .and. allocated(keep0)) exit                     ! rigid patches: trimming is invariant
       ip = wsurf%tp(t);  ib = wsurf%tb(t)
       if (.not. body_closed(ib)) cycle
       cen = (wsurf%v(:,wsurf%t(1,t)) + wsurf%v(:,wsurf%t(2,t)) + wsurf%v(:,wsurf%t(3,t)))/3.0_dp
       call surf_closest(wsurf, cen, 1, ip, surfovl_tol*1.0001_dp, dist, q, tb, wbc, region)
       if (tb <= 0) cycle
       if (dist > surfovl_tol) cycle
       if (wsurf%tb(tb) /= ib) cycle
       iq = wsurf%tp(tb)
       if (abs(tri_area(tb) - tri_area(t)) <= 1.0e-6_dp*tri_area(t)) then
          ! equal size: an exact DUPLICATE (a closing strip copied from a
          ! block, two patches of one surface) - keep the triangle of the
          ! patch with more triangles, then the lower patch index
          if (.not. (np_tri(iq) > np_tri(ip) .or. (np_tri(iq) == np_tri(ip) .and. iq < ip))) cycle
          if (dist > 1.0e-6_dp*surfovl_tol .and. dist > 1.0e-9_dp) then
             if (dist > surfovl_tol) cycle
          end if
       else if (.not. (tri_area(tb) < tri_area(t))) then
          cycle                                                          ! only a finer sheet covers
       end if
       ! interior of Q: farther than Surfovl_tol from Q's boundary nodes
       iq = wsurf%tp(tb)
       call patch_dims(wp_blk(iq), wp_bc(iq), ax, sd, t1, t2, n1, n2)
       v0 = wp_v0(iq);  dq = huge(1.0_dp)
       do m2 = 1, n2
          do m1 = 1, n1
             if (m1 /= 1 .and. m1 /= n1 .and. m2 /= 1 .and. m2 /= n2) cycle
             dq = min(dq, sum((wsurf%v(:, v0 + m1 + n1*(m2-1)) - q)**2))         ! (index fixed in update 24: was v0 + m1 - 1)
          end do
       end do
       if (sqrt(dq) <= surfovl_tol) cycle
       wsurf%keep(t) = .false.;  ntrim = ntrim + 1
    end do
    if (.not. cache_ok) then
       if (allocated(keep0)) deallocate(keep0)
       allocate(keep0(wsurf%nt));  keep0 = wsurf%keep
    end if
    ! wall-face force weights from the kept triangles (sixdof.F90)
    if (allocated(wpw)) then
       do ip = 1, nwpw
          if (allocated(wpw(ip)%a)) deallocate(wpw(ip)%a)
       end do
       deallocate(wpw, wpw_blk, wpw_bc)
    end if
    nwpw = nwp
    allocate(wpw(max(nwp,1)), wpw_blk(max(nwp,1)), wpw_bc(max(nwp,1)))
    do ip = 1, nwp
       wpw_blk(ip) = wp_blk(ip);  wpw_bc(ip) = wp_bc(ip)
       allocate(wpw(ip)%a(max(wp_n1(ip)-1,1), max(wp_n2(ip)-1,1)));  wpw(ip)%a = 0.0_dp
    end do
    do t = 1, wsurf%nt
       ip = wsurf%tp(t)
       lin = wsurf%t(1,t) - wp_v0(ip) - 1          ! first vertex = node (m1,m2) of the quad
       m1 = mod(lin, wp_n1(ip)) + 1;  m2 = lin/wp_n1(ip) + 1
       if (m1 < 1 .or. m1 > wp_n1(ip)-1 .or. m2 < 1 .or. m2 > wp_n2(ip)-1) cycle
       if (wsurf%keep(t)) wpw(ip)%a(m1,m2) = wpw(ip)%a(m1,m2) + 0.5_dp
    end do
    call build_hsurf
    call surf_wind_prep(hsurf, nbody)
    if (rank0v) write(*,'(A,I0,A,I0,A)') '  surface overset : ', ntrim, ' of ', wsurf%nt, &
       ' wall triangles covered by a finer overlapping patch are excluded from the hole cutting'
    ! ---- (b) winding-number test ------------------------------------------
    if (.not. allocated(selfsol)) then
       allocate(selfsol(nblk))
       do b = 1, nblk
          allocate(selfsol(b)%a(blk(b)%nx, blk(b)%ny, blk(b)%nz));  selfsol(b)%a = 0
       end do
    end if
    twind = MPI_Wtime();  wind_nq = 0;  wind_ntri = 0;  wind_ncl = 0
    do ib = 1, nbody
       if (.not. body_closed(ib)) cycle
       kobj = body_object(ib)
       bmin = huge(1.0_dp);  bmax = -huge(1.0_dp)
       do t = 1, wsurf%nt
          if (wsurf%tb(t) /= ib) cycle
          do nn = 1, 3
             bmin = min(bmin, wsurf%v(:,wsurf%t(nn,t)));  bmax = max(bmax, wsurf%v(:,wsurf%t(nn,t)))
          end do
       end do
       marg = 1.0e-6_dp*maxval(bmax - bmin)
       nsol = 0;  ncross = 0
       do b = 1, nblk
          ! blocks of the body's own (moving) object: the cut is rigid ->
          ! restore the cached marks of the first assembly
          selfc = (kobj > 0 .and. mgroup(b) == kobj)
          if (selfc .and. cache_ok) then
             do k = 1, blk(b)%nz
              do j = 1, blk(b)%ny
               do i = 1, blk(b)%nx
                  if (selfsol(b)%a(i,j,k) == ib) then
                     if (typ(b)%a(i,j,k) /= T_SOLID) nsol = nsol + 1
                     typ(b)%a(i,j,k) = T_SOLID
                  end if
               end do
              end do
             end do
             cycle
          end if
          do k = 1, blk(b)%nz
           do j = 1, blk(b)%ny
            do i = 1, blk(b)%nx
               if (typ(b)%a(i,j,k) == T_SOLID) cycle
               if (maxval(blk(b)%xn(i:i+1,j:j+1,k:k+1)) < bmin(1)-marg .or. minval(blk(b)%xn(i:i+1,j:j+1,k:k+1)) > bmax(1)+marg) cycle
               if (maxval(blk(b)%yn(i:i+1,j:j+1,k:k+1)) < bmin(2)-marg .or. minval(blk(b)%yn(i:i+1,j:j+1,k:k+1)) > bmax(2)+marg) cycle
               if (maxval(blk(b)%zn(i:i+1,j:j+1,k:k+1)) < bmin(3)-marg .or. minval(blk(b)%zn(i:i+1,j:j+1,k:k+1)) > bmax(3)+marg) cycle
               p = cell_center(b, i, j, k)
               ! only the cut layer (centre within 0.6 cell diagonals of the
               ! surface) is decided by the winding number; the cells deeper
               ! inside the body are found by the flood fill in step (c)
               diag = sqrt((maxval(blk(b)%xn(i:i+1,j:j+1,k:k+1)) - minval(blk(b)%xn(i:i+1,j:j+1,k:k+1)))**2 + &
                           (maxval(blk(b)%yn(i:i+1,j:j+1,k:k+1)) - minval(blk(b)%yn(i:i+1,j:j+1,k:k+1)))**2 + &
                           (maxval(blk(b)%zn(i:i+1,j:j+1,k:k+1)) - minval(blk(b)%zn(i:i+1,j:j+1,k:k+1)))**2)
               ! a cell that the surface passes through has its centre
               ! within half a diagonal of it: that is the cut layer; the
               ! flood fill takes care of everything behind it
               call surf_closest(hsurf, p, 2, ib, 0.6_dp*diag, dist, q, tb, wbc, region)
               if (tb <= 0) cycle
               if (dist > 0.6_dp*diag) cycle
               ! a cell next to its OWN wall whose nearest surface point lies
               ! on that wall (or on another patch of the same surface within
               ! Surfovl_tol - a coincident overlap) is outside its body by
               ! construction: never cut by a micron-level sign decision
               ! (mirrored grids would otherwise get different answers)
               if (near_own_wall(b, i, j, k)) then
                  if (own_surface_at(b, q)) cycle
               end if
               ! a body THINNER than the cell (a blade tip inside a
               ! background cell): the centre and the nodes are outside, yet
               ! the surface passes through the cell - the closest surface
               ! point q lies strictly inside it.  Such a cell is a hole too:
               ! a solved cell must never reach into a body (update 12c; the
               ! body's own grid solves around the body there)
               call inv_trilinear(blk(b)%xn, blk(b)%yn, blk(b)%zn, blk(b)%ni, blk(b)%nj, blk(b)%nk, i, j, k, q, rq, okq)
               if (okq) then
                  if (all(rq > 1.0e-3_dp) .and. all(rq < 1.0_dp - 1.0e-3_dp)) then
                     typ(b)%a(i,j,k) = T_SOLID;  nsol = nsol + 1;  ncross = ncross + 1
                     if (selfc) selfsol(b)%a(i,j,k) = ib
                     cycle
                  end if
               end if
               w = surf_wind(hsurf, p, ib)
               inside = (w > 0.5_dp)
               if (.not. inside) then
                  nodes: do t = 0, 7
                     p = [blk(b)%xn(i+iand(t,1), j+iand(t/2,1), k+iand(t/4,1)),   &
                          blk(b)%yn(i+iand(t,1), j+iand(t/2,1), k+iand(t/4,1)),   &
                          blk(b)%zn(i+iand(t,1), j+iand(t/2,1), k+iand(t/4,1))]
                     w = surf_wind(hsurf, p, ib)
                     if (w > 0.5_dp + 1.0e-2_dp) then     ! a node inside the body (nodes ON the surface give 0.5; the multipole
                                                          ! sums are accurate to ~1e-3, so a margin of 1e-2 keeps them out)
                        inside = .true.;  exit nodes
                     end if
                  end do nodes
               end if
               if (inside) then
                  typ(b)%a(i,j,k) = T_SOLID;  nsol = nsol + 1
                  if (selfc) selfsol(b)%a(i,j,k) = ib
               end if
            end do
           end do
          end do
       end do
       body_nsol(ib) = nsol
       if (rank0v) write(*,'(A,I0,A,I0,A,I0,A)') '    body ', ib, ': ', nsol, &
          ' cut-layer cells inside the solid -> holes (winding number; ', ncross, ' of them crossed by the surface, thinner than the cell)'
       if (rank0v) write(*,'(A,I0,A,I0,A,I0,A,F8.2,A)') '      winding number: ', wind_nq, ' queries, ', wind_ntri, &
          ' exact triangles, ', wind_ncl, ' clusters visited, ', MPI_Wtime() - twind, ' s'
       wind_nq = 0;  wind_ntri = 0;  wind_ncl = 0;  twind = MPI_Wtime()
    end do
    ! ---- (c) flood fill from the non-wall boundary faces of every block ----
    nflood = 0
    do b = 1, nblk
       ni_ = blk(b)%nx;  nj_ = blk(b)%ny;  nk_ = blk(b)%nz
       allocate(lab(ni_,nj_,nk_));  lab = 0                  ! 0 unvisited, 1 reached, 2 hole
       do k = 1, nk_;  do j = 1, nj_;  do i = 1, ni_
          if (typ(b)%a(i,j,k) == T_SOLID) lab(i,j,k) = 2
       end do;  end do;  end do
       allocate(queue(3, ni_*nj_*nk_));  nq = 0
       ! seeds: cells on the six block faces, except cells whose face lies
       ! on a no-slip wall patch of this block
       do k = 1, nk_;  do j = 1, nj_;  do i = 1, ni_
          if (i /= 1 .and. i /= ni_ .and. j /= 1 .and. j /= nj_ .and. k /= 1 .and. k /= nk_) cycle
          if (lab(i,j,k) /= 0) cycle
          if (seed_ok(b, i, j, k)) then
             ! a boundary cell deep inside a body (e.g. on the symmetry
             ! plane that closes a half body, or an outer face of an
             ! appendage block inside the hub) must not seed the fill
             inside = .false.
             p = cell_center(b, i, j, k)
             do ib = 1, nbody
                if (.not. body_closed(ib)) cycle
                if (surf_wind(hsurf, p, ib) > 0.5_dp) then
                   inside = .true.;  exit
                end if
             end do
             if (inside) then
                lab(i,j,k) = 2;  typ(b)%a(i,j,k) = T_SOLID;  nflood = nflood + 1
             else
                nq = nq + 1;  queue(:,nq) = [i, j, k];  lab(i,j,k) = 1
             end if
          end if
       end do;  end do;  end do
       iq0 = 0
       do while (iq0 < nq)
          iq0 = iq0 + 1;  dd = queue(:,iq0)
          do nb_dir = 1, 6
             ii = dd(1);  jj = dd(2);  kk = dd(3)
             select case (nb_dir)
             case (1);  ii = ii - 1
             case (2);  ii = ii + 1
             case (3);  jj = jj - 1
             case (4);  jj = jj + 1
             case (5);  kk = kk - 1
             case default;  kk = kk + 1
             end select
             if (ii < 1 .or. ii > ni_ .or. jj < 1 .or. jj > nj_ .or. kk < 1 .or. kk > nk_) cycle
             if (lab(ii,jj,kk) /= 0) cycle
             lab(ii,jj,kk) = 1;  nq = nq + 1;  queue(:,nq) = [ii, jj, kk]
          end do
       end do
       do k = 1, nk_;  do j = 1, nj_;  do i = 1, ni_
          if (lab(i,j,k) == 0) then
             typ(b)%a(i,j,k) = T_SOLID;  nflood = nflood + 1
          end if
       end do;  end do;  end do
       deallocate(lab, queue)
    end do
    if (rank0v) write(*,'(A,I0,A)') '    flood fill      : ', nflood, &
       ' cells enclosed by the cut layer and the block walls -> holes (body interiors, appendage roots, duct penetrations)'
  contains
    ! hole-cutting surface = kept wall triangles + their mirror images
    ! across every PLANAR domain-boundary patch (symmetry plane, ...)
    ! touched by the body: a half body closed by a symmetry plane then
    ! has a closed surface for the winding number
    subroutine build_hsurf
      integer :: t, iq, nq, nmir, nv, k, nplane, ipl, ib2, nbt, ncap, nslab, c
      real(dp) :: pn(3,64), pc(3,64), e1(3), e2(3), nl, L(3), v1(3), v2(3), v3(3), n_(3), d
      logical :: planar(64), touch(64), slab(64)
      logical, allocatable :: btouch(:,:)
      real(dp), allocatable :: cv(:,:), ctn(:,:)
      integer, allocatable :: ctp(:), ctb(:)
      integer :: pl_of(64), vv(3)
      ! planar domain patches
      nplane = 0;  nq = 0
      do t = 1, dsurf%nt
         nq = max(nq, dsurf%tp(t))
      end do
      nq = min(nq, 64)
      do iq = 1, nq
         planar(iq) = .false.
         do t = 1, dsurf%nt
            if (dsurf%tp(t) == iq) then
               v1 = dsurf%v(:,dsurf%t(1,t));  v2 = dsurf%v(:,dsurf%t(2,t));  v3 = dsurf%v(:,dsurf%t(3,t))
               e1 = v2 - v1;  e2 = v3 - v1
               n_ = [e1(2)*e2(3)-e1(3)*e2(2), e1(3)*e2(1)-e1(1)*e2(3), e1(1)*e2(2)-e1(2)*e2(1)]
               nl = sqrt(sum(n_*n_))
               if (nl > 0.0_dp) then
                  pn(:,iq) = n_/nl;  pc(:,iq) = v1;  planar(iq) = .true.;  exit
               end if
            end if
         end do
         if (.not. planar(iq)) cycle
         L = dsurf%bmax - dsurf%bmin
         do t = 1, dsurf%nt
            if (dsurf%tp(t) /= iq) cycle
            do k = 1, 3
               d = abs(sum((dsurf%v(:,dsurf%t(k,t)) - pc(:,iq))*pn(:,iq)))
               if (d > 1.0e-6_dp*maxval(L)) planar(iq) = .false.
            end do
         end do
      end do
      ! distinct planes only (several block faces lie on the same
      ! symmetry plane): keep the first patch of every plane
      do iq = 1, nq
         if (.not. planar(iq)) cycle
         do ipl = 1, iq-1
            if (.not. planar(ipl)) cycle
            if (abs(abs(sum(pn(:,iq)*pn(:,ipl))) - 1.0_dp) < 1.0e-6_dp .and.                 &
                abs(sum((pc(:,iq) - pc(:,ipl))*pn(:,ipl))) < 1.0e-6_dp*maxval(dsurf%bmax - dsurf%bmin)) then
               planar(iq) = .false.;  exit
            end if
         end do
      end do
      ! which planes does each body touch (a wall vertex within Surfovl_tol)
      nbt = max(1, maxval(wsurf%tb(1:max(wsurf%nt,1))))
      allocate(btouch(nbt, 64));  btouch = .false.
      do iq = 1, nq
         touch(iq) = .false.
         if (.not. planar(iq)) cycle
         do t = 1, wsurf%nt
            if (.not. wsurf%keep(t)) cycle
            if (wsurf%tb(t) < 1) cycle
            do k = 1, 3
               if (abs(sum((wsurf%v(:,wsurf%t(k,t)) - pc(:,iq))*pn(:,iq))) <= surfovl_tol) then
                  touch(iq) = .true.;  btouch(wsurf%tb(t), iq) = .true.;  exit
               end if
            end do
         end do
      end do
      ! SLAB planes (24c): a body spanning two PARALLEL touched planes with the same touching bodies -
      ! the two span walls of an extruded quasi-2-D grid - is closed by end caps on both planes
      ! (slab_caps) instead of being mirrored: the mirror images across one span wall overlap the
      ! domain across the other when the span is thin, and the body was left open (24c: on the
      ! 0.4-mm-span FSI-PfS-2a grid 148 active cells remained inside the plate)
      slab = .false.
      do iq = 1, nq
         if (.not. touch(iq)) cycle
         do ipl = 1, nq
            if (ipl == iq .or. .not. touch(ipl)) cycle
            if (abs(abs(sum(pn(:,iq)*pn(:,ipl))) - 1.0_dp) > 1.0e-6_dp) cycle
            if (abs(sum((pc(:,iq) - pc(:,ipl))*pn(:,ipl))) <= 1.0e-6_dp*maxval(dsurf%bmax - dsurf%bmin)) cycle
            if (all(btouch(:, iq) .eqv. btouch(:, ipl))) slab(iq) = .true.
         end do
      end do
      call slab_caps(nq, slab, pn, pc, btouch, nbt, ncap, cv, ctn, ctp, ctb)
      nmir = 0;  nslab = 0
      do iq = 1, nq
         if (slab(iq)) nslab = nslab + 1
         if (slab(iq)) touch(iq) = .false.                    ! closed by the caps, not mirrored
         if (touch(iq)) nmir = nmir + 1
      end do
      nv = count(wsurf%keep)
      call surf_init(hsurf, (1+nmir)*wsurf%nv + 3*ncap, (1+nmir)*nv + ncap)
      ! original kept triangles (vertex numbering preserved: copy all vertices)
      hsurf%nv = wsurf%nv;  hsurf%v(:,1:wsurf%nv) = wsurf%v(:,1:wsurf%nv)
      hsurf%nt = 0
      do t = 1, wsurf%nt
         if (.not. wsurf%keep(t)) cycle
         hsurf%nt = hsurf%nt + 1
         hsurf%t(:,hsurf%nt) = wsurf%t(:,t);  hsurf%tn(:,hsurf%nt) = wsurf%tn(:,t)
         hsurf%tp(hsurf%nt) = wsurf%tp(t);  hsurf%tb(hsurf%nt) = wsurf%tb(t)
      end do
      do iq = 1, nq
         if (.not. touch(iq)) cycle
         ! mirrored vertices
         do k = 1, wsurf%nv
            d = sum((wsurf%v(:,k) - pc(:,iq))*pn(:,iq))
            hsurf%v(:,hsurf%nv+k) = wsurf%v(:,k) - 2.0_dp*d*pn(:,iq)
         end do
         do t = 1, wsurf%nt
            if (.not. wsurf%keep(t)) cycle
            hsurf%nt = hsurf%nt + 1
            hsurf%t(:,hsurf%nt) = wsurf%t(:,t) + hsurf%nv
            n_ = wsurf%tn(:,t);  d = sum(n_*pn(:,iq))
            hsurf%tn(:,hsurf%nt) = n_ - 2.0_dp*d*pn(:,iq)
            hsurf%tp(hsurf%nt) = wsurf%tp(t);  hsurf%tb(hsurf%nt) = wsurf%tb(t)
         end do
         hsurf%nv = hsurf%nv + wsurf%nv
      end do
      ! end caps of the slab bodies
      do c = 1, ncap
         hsurf%v(:, hsurf%nv+1:hsurf%nv+3) = cv(:, 3*c-2:3*c)
         hsurf%nt = hsurf%nt + 1
         hsurf%t(:, hsurf%nt) = [hsurf%nv+1, hsurf%nv+2, hsurf%nv+3]
         hsurf%tn(:, hsurf%nt) = ctn(:, c);  hsurf%tp(hsurf%nt) = ctp(c);  hsurf%tb(hsurf%nt) = ctb(c)
         hsurf%nv = hsurf%nv + 3
      end do
      if (rank0v .and. ncap > 0) write(*,'(A,I0,A,I0,A)') '  hole-cutting surface: bodies closed by ', ncap, &
         ' end-cap triangles on ', nslab, ' slab planes (quasi-2-D span walls, update 24c)'
      deallocate(btouch, cv, ctn, ctp, ctb)
      call surf_finish(hsurf)
      if (allocated(hsurf%keep)) deallocate(hsurf%keep)
      allocate(hsurf%keep(hsurf%nt));  hsurf%keep = .true.
      if (rank0v .and. nmir > 0) write(*,'(A,I0,A)') '  hole-cutting surface: body mirrored across ', nmir, &
         ' planar domain-boundary patch(es) (symmetry plane closure)'
    end subroutine build_hsurf

    integer function np_tri(ip)
      integer, intent(in) :: ip
      np_tri = count(wsurf%tp(1:wsurf%nt) == ip)
    end function np_tri

    real(dp) function tri_area(t)
      integer, intent(in) :: t
      real(dp) :: e1(3), e2(3), c(3)
      e1 = wsurf%v(:,wsurf%t(2,t)) - wsurf%v(:,wsurf%t(1,t));  e2 = wsurf%v(:,wsurf%t(3,t)) - wsurf%v(:,wsurf%t(1,t))
      c = [e1(2)*e2(3)-e1(3)*e2(2), e1(3)*e2(1)-e1(1)*e2(3), e1(1)*e2(2)-e1(2)*e2(1)]
      tri_area = 0.5_dp*sqrt(sum(c*c))
    end function tri_area
    ! a boundary cell may seed the flood fill unless every block face it
    ! touches is covered by a no-slip wall patch at that position
    logical function seed_ok(b, i, j, k)
      integer, intent(in) :: b, i, j, k
      integer :: n, ax, sd, t1, t2, n1, n2, idx(3), rr(2,3)
      logical :: onface(6), walled(6)
      integer :: f
      idx = [i, j, k]
      onface = [i == 1, i == blk(b)%nx, j == 1, j == blk(b)%ny, k == 1, k == blk(b)%nz]
      walled = .false.
      do n = 1, blk(b)%nbc
         if (.not. is_wall(b, n)) cycle
         ax = abs(blk(b)%bc_dir(n));  sd = merge(1, 2, blk(b)%bc_dir(n) > 0)
         f = 2*(ax-1) + sd
         if (.not. onface(f)) cycle
         rr(:,1) = blk(b)%bc_irng(:,n);  rr(:,2) = blk(b)%bc_jrng(:,n);  rr(:,3) = blk(b)%bc_krng(:,n)
         if (all(idx >= rr(1,:)) .and. all(idx <= max(rr(2,:)-1, rr(1,:)))) walled(f) = .true.
      end do
      seed_ok = any(onface .and. .not. walled)
    end function seed_ok
  end subroutine cut_holes

  pure function cell_center(b, i, j, k) result(c)
    integer, intent(in) :: b, i, j, k
    real(dp) :: c(3)
    c(1) = 0.125_dp*sum(blk(b)%xn(i:i+1,j:j+1,k:k+1))
    c(2) = 0.125_dp*sum(blk(b)%yn(i:i+1,j:j+1,k:k+1))
    c(3) = 0.125_dp*sum(blk(b)%zn(i:i+1,j:j+1,k:k+1))
  end function cell_center

  !=====================================================================
  ! 1b. Cells behind uncovered face cells are forced non-active.
  !=====================================================================
  !---------------------------------------------------------------------
  ! End caps of the bodies on the slab planes (24c, build_hsurf): on every
  ! slab plane iq and for every body touching it, the edges of the body's
  ! kept wall triangles lying in the plane (two vertices within
  ! Surfovl_tol) are chained by coincidence into closed loops - the
  ! outlines of the body's cross-section - and every loop is triangulated
  ! by ear clipping in plane coordinates.  The normal of a cap triangle
  ! ("into the fluid", the winding-number orientation) points out of the
  ! body: along the plane normal away from the body's vertices off the
  ! plane.  A loop that does not close is reported and skipped.
  !---------------------------------------------------------------------
  subroutine slab_caps(nq, slab, pn, pc, btouch, nbt, ncap, cv, ctn, ctp, ctb)
    integer, intent(in) :: nq, nbt
    logical, intent(in) :: slab(:), btouch(:,:)
    real(dp), intent(in) :: pn(:,:), pc(:,:)
    integer, intent(out) :: ncap
    real(dp), allocatable, intent(out) :: cv(:,:), ctn(:,:)
    integer, allocatable, intent(out) :: ctp(:), ctb(:)
    integer :: iq, bb, t, k, ne, nin, a(3), e, cur, nl, np_, m, i0, i1, i2, guard, nfail
    integer, allocatable :: ea(:,:), lp(:)
    logical, allocatable :: used(:)
    real(dp) :: tolc, e1(3), e2(3), dsum, sgn, xe(3), tn_(3), dmax
    real(dp), allocatable :: P2(:,:), P3(:,:), cvt(:,:), ctnt(:,:)
    integer, allocatable :: ctpt(:), ctbt(:)
    logical, allocatable :: alive(:)
    integer :: cap_max
    cap_max = max(16, 4*wsurf%nt)
    allocate(cvt(3, 3*cap_max), ctnt(3, cap_max), ctpt(cap_max), ctbt(cap_max))
    ncap = 0;  nfail = 0
    tolc = max(surfovl_tol, 1.0e-9_dp*maxval(wsurf%bmax - wsurf%bmin))
    allocate(ea(2, max(1, wsurf%nt)), used(max(1, wsurf%nt)), lp(max(1, wsurf%nt) + 1))
    do iq = 1, nq
       if (.not. slab(iq)) cycle
       ! plane basis
       e1 = [1.0_dp, 0.0_dp, 0.0_dp]
       if (abs(pn(1,iq)) > 0.9_dp) e1 = [0.0_dp, 1.0_dp, 0.0_dp]
       e1 = e1 - sum(e1*pn(:,iq))*pn(:,iq);  e1 = e1/sqrt(sum(e1*e1))
       e2 = [pn(2,iq)*e1(3) - pn(3,iq)*e1(2), pn(3,iq)*e1(1) - pn(1,iq)*e1(3), pn(1,iq)*e1(2) - pn(2,iq)*e1(1)]
       do bb = 1, nbt
          if (.not. btouch(bb, iq)) cycle
          ! in-plane edges of the body's kept triangles; the side of the body (vertices off the plane)
          ne = 0;  dsum = 0.0_dp;  dmax = 0.0_dp
          do t = 1, wsurf%nt
             if (.not. wsurf%keep(t) .or. wsurf%tb(t) /= bb) cycle
             nin = 0
             do k = 1, 3
                xe(1) = sum((wsurf%v(:,wsurf%t(k,t)) - pc(:,iq))*pn(:,iq))
                if (abs(xe(1)) <= tolc) then
                   nin = nin + 1;  a(nin) = wsurf%t(k,t)
                else if (abs(xe(1)) > dmax) then
                   dmax = abs(xe(1));  dsum = xe(1)
                end if
             end do
             if (nin == 2) then
                ne = ne + 1;  ea(:, ne) = a(1:2)
             end if
          end do
          if (ne < 3) cycle
          sgn = merge(-1.0_dp, 1.0_dp, dsum > 0.0_dp)        ! the cap normal points away from the body
          tn_ = sgn*pn(:,iq)
          used = .false.
          do e = 1, ne
             if (used(e)) cycle
             ! chain a loop starting with edge e
             nl = 2;  lp(1) = ea(1,e);  lp(2) = ea(2,e);  used(e) = .true.
             guard = 0
             do
                guard = guard + 1
                if (guard > ne + 2) exit
                if (coinc(lp(nl), lp(1)) .and. nl > 3) exit
                cur = 0
                do m = 1, ne
                   if (used(m)) cycle
                   if (coinc(ea(1,m), lp(nl))) then
                      cur = m;  i0 = ea(2,m);  exit
                   else if (coinc(ea(2,m), lp(nl))) then
                      cur = m;  i0 = ea(1,m);  exit
                   end if
                end do
                if (cur == 0) exit
                used(cur) = .true.;  nl = nl + 1;  lp(nl) = i0
             end do
             if (.not. coinc(lp(nl), lp(1)) .or. nl < 4) then
                nfail = nfail + 1;  cycle
             end if
             np_ = nl - 1                                   ! the closing vertex repeats the first
             allocate(P2(2, np_), P3(3, np_), alive(np_))
             do m = 1, np_
                P3(:, m) = wsurf%v(:, lp(m))
                P2(1, m) = sum((P3(:, m) - pc(:,iq))*e1);  P2(2, m) = sum((P3(:, m) - pc(:,iq))*e2)
             end do
             call ear_clip(np_, P2, P3, tn_, iq, bb)
             deallocate(P2, P3, alive)
          end do
       end do
    end do
    if (nfail > 0 .and. rank0v) write(*,'(A,I0,A)') '  NOTE (hole cutting): ', nfail, &
       ' outline loop(s) of a slab body on a span plane did not close - no end cap for them'
    allocate(cv(3, max(1, 3*ncap)), ctn(3, max(1, ncap)), ctp(max(1, ncap)), ctb(max(1, ncap)))
    if (ncap > 0) then
       cv(:, 1:3*ncap) = cvt(:, 1:3*ncap);  ctn(:, 1:ncap) = ctnt(:, 1:ncap)
       ctp(1:ncap) = ctpt(1:ncap);  ctb(1:ncap) = ctbt(1:ncap)
    end if
    deallocate(cvt, ctnt, ctpt, ctbt, ea, used, lp)
  contains
    logical function coinc(ia, ib_)
      integer, intent(in) :: ia, ib_
      coinc = (ia == ib_)
      if (.not. coinc) coinc = (maxval(abs(wsurf%v(:, ia) - wsurf%v(:, ib_))) <= tolc)
    end function coinc
    ! ear clipping of a simple polygon (np vertices, plane coordinates P2)
    subroutine ear_clip(np, P2, P3, tnc, iqq, bbb)
      integer, intent(in) :: np, iqq, bbb
      real(dp), intent(in) :: P2(:,:), P3(:,:), tnc(3)
      integer :: nleft, ia, ib0, ic, j, tries
      integer, allocatable :: nx(:), pv(:)
      real(dp) :: area2, orient, eps2
      logical :: ear
      allocate(nx(np), pv(np))
      do j = 1, np
         nx(j) = mod(j, np) + 1;  pv(j) = mod(j - 2 + np, np) + 1
      end do
      area2 = 0.0_dp
      do j = 1, np
         area2 = area2 + P2(1,j)*P2(2,nx(j)) - P2(1,nx(j))*P2(2,j)
      end do
      orient = sign(1.0_dp, area2)
      nleft = np;  ib0 = 1;  tries = 0
      eps2 = 1.0e-12_dp*abs(area2)
      do while (nleft > 3 .and. tries < 2*nleft)
         ia = pv(ib0);  ic = nx(ib0)
         if (abs(cr(P2(:,ia), P2(:,ib0), P2(:,ic))) <= eps2) then     ! collinear vertex: dropped, no triangle
            nx(ia) = ic;  pv(ic) = ia;  nleft = nleft - 1;  ib0 = ic;  tries = 0
            cycle
         end if
         ear = (orient*cr(P2(:,ia), P2(:,ib0), P2(:,ic)) > 0.0_dp)
         if (ear) then
            j = nx(ic)
            do while (j /= ia)
               if (in_tri(P2(:,j), P2(:,ia), P2(:,ib0), P2(:,ic))) then
                  ear = .false.;  exit
               end if
               j = nx(j)
            end do
         end if
         if (ear) then
            call add_cap(P3(:,ia), P3(:,ib0), P3(:,ic), tnc, iqq, bbb)
            nx(ia) = ic;  pv(ic) = ia;  nleft = nleft - 1;  ib0 = ic;  tries = 0
         else
            ib0 = nx(ib0);  tries = tries + 1
         end if
      end do
      if (nleft == 3) then
         ia = pv(ib0);  ic = nx(ib0)
         call add_cap(P3(:,ia), P3(:,ib0), P3(:,ic), tnc, iqq, bbb)
      else
         nfail = nfail + 1
      end if
      deallocate(nx, pv)
    end subroutine ear_clip
    pure real(dp) function cr(a_, b_, c_)
      real(dp), intent(in) :: a_(2), b_(2), c_(2)
      cr = (b_(1) - a_(1))*(c_(2) - a_(2)) - (b_(2) - a_(2))*(c_(1) - a_(1))
    end function cr
    pure logical function in_tri(q_, a_, b_, c_)
      real(dp), intent(in) :: q_(2), a_(2), b_(2), c_(2)
      real(dp) :: d1, d2, d3
      d1 = cr(a_, b_, q_);  d2 = cr(b_, c_, q_);  d3 = cr(c_, a_, q_)
      in_tri = ((d1 > 0.0_dp .and. d2 > 0.0_dp .and. d3 > 0.0_dp) .or. (d1 < 0.0_dp .and. d2 < 0.0_dp .and. d3 < 0.0_dp))
    end function in_tri
    subroutine add_cap(xa, xb, xc, tnc, iqq, bbb)
      real(dp), intent(in) :: xa(3), xb(3), xc(3), tnc(3)
      integer, intent(in) :: iqq, bbb
      if (ncap >= cap_max) return
      ncap = ncap + 1
      cvt(:, 3*ncap-2) = xa;  cvt(:, 3*ncap-1) = xb;  cvt(:, 3*ncap) = xc
      ctnt(:, ncap) = tnc;  ctpt(ncap) = iqq;  ctbt(ncap) = bbb
    end subroutine add_cap
  end subroutine slab_caps

  subroutine mark_forced
    integer :: b, ax, sd, t1, t2, c1, c2, l, idx(3), a1, a2, rad
    rad = ovst_nf - 1
    do b = 1, nblk
       do ax = 1, 3
          call tang(ax, t1, t2)
          do sd = 1, 2
             if (all(cov(sd,ax,b)%a)) cycle
             do c2 = 1, ncell(b,t2)
                do c1 = 1, ncell(b,t1)
                   if (cov(sd,ax,b)%a(c1,c2)) cycle
                   do l = 1, ovst_nf
                      idx(ax) = merge(l, ncell(b,ax)-l+1, sd == 1)
                      do a2 = max(1,c2-rad), min(ncell(b,t2),c2+rad)
                         do a1 = max(1,c1-rad), min(ncell(b,t1),c1+rad)
                            idx(t1) = a1;  idx(t2) = a2
                            if (typ(b)%a(idx(1),idx(2),idx(3)) /= T_SOLID .and. typ(b)%a(idx(1),idx(2),idx(3)) /= T_BLH) &
                               typ(b)%a(idx(1),idx(2),idx(3)) = T_FORCED
                         end do
                      end do
                   end do
                end do
             end do
          end do
       end do
    end do
  end subroutine mark_forced

  !=====================================================================
  ! 4a. Finest containing non-solid cell of every cell centre in the
  !     other blocks -> donor; dominated when the donor is finer.
  !=====================================================================

  ! locate p in block ob for a cell (i,j,k) of block b: the search of the
  ! block in itself excludes the cells within 2 Ovst_fringe + 1 index
  ! layers of (i,j,k) - self-overlap (the two ends of a C / O grid that
  ! overlap without a point match: the KP505 tip-cap grids) is then an
  ! ordinary overset overlap, the far part of the block acting as donor
  subroutine locate_from(b, i, j, k, ob, p, ic, jc, kc, rst, found)
    integer,  intent(in)  :: b, i, j, k, ob
    real(dp), intent(in)  :: p(3)
    integer,  intent(out) :: ic, jc, kc
    real(dp), intent(out) :: rst(3)
    logical,  intent(out) :: found
    if (ob == b) then
       call locate_point(bk(ob), blk(ob)%xn, blk(ob)%yn, blk(ob)%zn, p, ic, jc, kc, rst, found, [i, j, k], 2*ovst_nf + 1)
    else
       call locate_point(bk(ob), blk(ob)%xn, blk(ob)%yn, blk(ob)%zn, p, ic, jc, kc, rst, found)
    end if
  end subroutine locate_from
  !---------------------------------------------------------------------
  ! Cells of any block that lie inside the BOUNDARY-LAYER region (the
  ! Ovst_wall_layers cells next to a no-slip wall) of ANOTHER wall block
  ! are holes, never fringe cells: that region is computed by the grid
  ! that carries the wall, and blanking instead of fringing there keeps
  ! the fringe boundary of coarse Cartesian grids one layer away from the
  ! body surface (no background fringe cell touches the near-wall band).
  ! Cells next to their own wall are exempt (coincident-wall overlaps).
  !---------------------------------------------------------------------
  subroutine blank_in_boundary_layers
    integer :: b, i, j, k, ob, ic, jc, kc, nbl
    real(dp) :: p(3), rst(3)
    logical :: found
    nbl = 0
    do b = 1, nblk
       do k = 1, blk(b)%nz
        do j = 1, blk(b)%ny
         do i = 1, blk(b)%nx
            if (typ(b)%a(i,j,k) /= T_FREE) cycle
            if (near_own_wall(b, i, j, k)) cycle
            p = cell_center(b, i, j, k)
            do ob = 1, nblk
               if (ob == b) cycle
               if (nwall_blk(ob) == 0) cycle
               call locate_point(bk(ob), blk(ob)%xn, blk(ob)%yn, blk(ob)%zn, p, ic, jc, kc, rst, found)
               if (.not. found) cycle
               if (typ(ob)%a(ic,jc,kc) == T_SOLID) cycle
               ! the wall block's own cell there must be solvable: behind
               ! its overset face (FORCED, e.g. the last layers of a blade
               ! O-grid under the tip-cap grid) the wall block cannot
               ! compute the boundary layer, so this cell is not blanked
               ! (update 12c: the near-wall sliver at such an interface
               ! was solved by nobody)
               if (typ(ob)%a(ic,jc,kc) == T_FORCED) cycle
               if (near_own_wall(ob, ic, jc, kc)) then
                  typ(b)%a(i,j,k) = T_BLH;  nbl = nbl + 1;  exit
               end if
            end do
         end do
        end do
       end do
    end do
    n_blh = nbl
    if (rank0v) write(*,'(A,I0,A)') '  boundary layers  : ', nbl, &
       ' cells inside the near-wall band of another wall block -> holes (never fringe)'
  end subroutine blank_in_boundary_layers

  subroutine find_cell_donors
    integer :: b, i, j, k, ob, ic, jc, kc, lin, fbb, fbc(3), vdc(4), gb, ngm
    real(dp) :: p(3), rst(3), vbest, vd, vfb, vdom
    logical :: found, use_cache, build_cache, stale
    n_fb_cells = 0;  n_resv = 0
    ! ---- donor cache: the blocks of the same motion group as b -----------
    if (.not. allocated(dch)) then
       allocate(dch(nblk))
       do b = 1, nblk
          ngm = 0
          do ob = 1, nblk
             if (mgroup(ob) == mgroup(b)) ngm = ngm + 1
          end do
          if (ngm < 2) cycle                                          ! nothing to cache
          allocate(dch(b)%dnb(blk(b)%nx,blk(b)%ny,blk(b)%nz), dch(b)%dnc(3,blk(b)%nx,blk(b)%ny,blk(b)%nz), &
                   dch(b)%fbb(blk(b)%nx,blk(b)%ny,blk(b)%nz), dch(b)%fbc(3,blk(b)%nx,blk(b)%ny,blk(b)%nz), &
                   dch(b)%vdc(4,blk(b)%nx,blk(b)%ny,blk(b)%nz), dch(b)%flag(blk(b)%nx,blk(b)%ny,blk(b)%nz), &
                   dch(b)%dnr(3,blk(b)%nx,blk(b)%ny,blk(b)%nz), dch(b)%vbest(blk(b)%nx,blk(b)%ny,blk(b)%nz), &
                   dch(b)%vdom(blk(b)%nx,blk(b)%ny,blk(b)%nz), dch(b)%vfb(blk(b)%nx,blk(b)%ny,blk(b)%nz))
          dch(b)%flag = 0
       end do
    end if
    do b = 1, nblk
       gb = mgroup(b)
       build_cache = allocated(dch(b)%flag) .and. .not. cache_ok
       use_cache   = allocated(dch(b)%flag) .and. cache_ok
       do k = 1, blk(b)%nz
        do j = 1, blk(b)%ny
         do i = 1, blk(b)%nx
            lin = i + blk(b)%nx*((j-1) + blk(b)%ny*(k-1))
            dnb(b)%a(i,j,k) = 0
            if (typ(b)%a(i,j,k) == T_SOLID .or. typ(b)%a(i,j,k) == T_BLH) cycle
            dfb(b)%a(i,j,k) = 0
            p = cell_center(b, i, j, k)
            vbest = huge(1.0_dp);  vfb = huge(1.0_dp);  fbb = 0;  vdom = huge(1.0_dp);  vdc = 0
            stale = .true.
            if (use_cache) then
               ! flag 1 = exact partial result over the same-group blocks
               ! (no containing cell was solid / in a boundary-layer band
               ! then); still valid when neither the cached donor nor the
               ! cached domination cell became solid / BLH since
               if (dch(b)%flag(i,j,k) == 1) then
                  stale = .false.
                  ob = dch(b)%dnb(i,j,k)
                  if (ob > 0) then
                     if (typ(ob)%a(dch(b)%dnc(1,i,j,k), dch(b)%dnc(2,i,j,k), dch(b)%dnc(3,i,j,k)) == T_SOLID .or. &
                         typ(ob)%a(dch(b)%dnc(1,i,j,k), dch(b)%dnc(2,i,j,k), dch(b)%dnc(3,i,j,k)) == T_BLH) stale = .true.
                  end if
                  ob = dch(b)%vdc(1,i,j,k)
                  if (ob > 0 .and. .not. stale) then
                     if (typ(ob)%a(dch(b)%vdc(2,i,j,k), dch(b)%vdc(3,i,j,k), dch(b)%vdc(4,i,j,k)) == T_SOLID .or. &
                         typ(ob)%a(dch(b)%vdc(2,i,j,k), dch(b)%vdc(3,i,j,k), dch(b)%vdc(4,i,j,k)) == T_BLH) stale = .true.
                  end if
               end if
               if (.not. stale) then
                  vbest = dch(b)%vbest(i,j,k);  vdom = dch(b)%vdom(i,j,k);  vfb = dch(b)%vfb(i,j,k)
                  fbb = dch(b)%fbb(i,j,k);  fbc = dch(b)%fbc(:,i,j,k);  vdc = dch(b)%vdc(:,i,j,k)
                  dnb(b)%a(i,j,k) = dch(b)%dnb(i,j,k)
                  if (dnb(b)%a(i,j,k) > 0) then
                     dnc(b)%a(:,lin,1) = dch(b)%dnc(:,i,j,k);  dnr(b)%a(:,lin,1) = dch(b)%dnr(:,i,j,k)
                  end if
               end if
            end if
            if (build_cache) dch(b)%flag(i,j,k) = 1
            do ob = 1, nblk
               if (.not. stale .and. mgroup(ob) == gb) cycle           ! taken from the cache
               call locate_from(b, i, j, k, ob, p, ic, jc, kc, rst, found)
               if (.not. found) cycle
               if (typ(ob)%a(ic,jc,kc) == T_SOLID .or. typ(ob)%a(ic,jc,kc) == T_BLH) then
                  ! the containing cell lies inside a body: FALLBACK donor =
                  ! nearest non-solid cell of that block within FB_R index
                  ! layers, used as a nearest-neighbour (0th order) donor.
                  ! Never used for the domination test.
                  if (build_cache .and. mgroup(ob) == gb) dch(b)%flag(i,j,k) = 2   ! not cacheable
                  call nearest_nonsolid(ob, ic, jc, kc, p, found)
                  if (.not. found) cycle
                  vd = vol(ob)%a(ic,jc,kc)
                  if (vd < vfb) then
                     vfb = vd;  fbb = ob;  fbc = [ic, jc, kc]
                  end if
                  cycle
               end if
               vd = vol(ob)%a(ic,jc,kc)
               if (vd < vbest*(1.0_dp - ovst_vtol)) then
                  vbest = vd
                  dnb(b)%a(i,j,k) = ob
                  dnc(b)%a(:,lin,1) = [ic, jc, kc];  dnr(b)%a(:,lin,1) = rst
               end if
               ! domination is decided by the finest containing cell that can
               ! be ACTIVE: a FORCED fringe cell of the other block (the layers
               ! behind its overset face) never dominates - the cells under
               ! those layers stay active here, so the two grids overlap
               ! instead of leaving a blanked gap
               if (typ(ob)%a(ic,jc,kc) /= T_FORCED .and. vd < vdom*(1.0_dp - ovst_vtol)) then
                  vdom = vd;  vdc = [ob, ic, jc, kc]
               end if
            end do
            if (build_cache) then
               ! the partial result over the same-group blocks: rebuild it
               ! from those blocks only (the loop above mixed all groups)
               call group_partial(b, i, j, k, p, gb)
            end if
            if (dnb(b)%a(i,j,k) > 0 .and. typ(b)%a(i,j,k) == T_FREE .and. .not. near_own_wall(b, i, j, k)) then
               if (vdom < vol(b)%a(i,j,k)*(1.0_dp - ovst_vtol)) typ(b)%a(i,j,k) = T_DOM
            end if
            ! RESERVED REGION of another wall block: a cell of wall block b
            ! that lies inside the grid of another wall block ob, at a place
            ! where ob's wall is NOT the same surface as b's wall (the
            ! closest point of b's own wall is farther than Surfovl_tol from
            ! ob's wall: the passage of a duct mouth, the region above an
            ! appendage), is never solved by b - regardless of the cell
            ! sizes: it is dominated (fringe next to active cells of b,
            ! hole elsewhere) and ob computes the flow there.  In a
            ! coincident-wall overlap (two patches of one hull surface) the
            ! normal rules apply.
            if (typ(b)%a(i,j,k) == T_FREE .and. nwall_blk(b) > 0 .and. dnb(b)%a(i,j,k) > 0) then
               if (reserved_by_other(b, dnb(b)%a(i,j,k), p)) then
                  typ(b)%a(i,j,k) = T_DOM;  n_resv = n_resv + 1
               end if
            end if
            if (dnb(b)%a(i,j,k) == 0 .and. fbb > 0) then     ! no regular donor: fallback
               dnb(b)%a(i,j,k) = fbb;  dfb(b)%a(i,j,k) = 1
               dnc(b)%a(:,lin,1) = fbc;  dnr(b)%a(:,lin,1) = 0.5_dp    ! cell centre -> weight 1
               n_fb_cells = n_fb_cells + 1
            end if
         end do
        end do
       end do
    end do
  contains
    ! cache entry of cell (i,j,k) of block b: the search restricted to the
    ! blocks of its own motion group gb
    subroutine group_partial(b, i, j, k, p, gb)
      integer,  intent(in) :: b, i, j, k, gb
      real(dp), intent(in) :: p(3)
      integer :: ob, ic, jc, kc, fbb, fbc(3), vdc(4)
      real(dp) :: rst(3), vbest, vd, vfb, vdom
      logical :: found
      vbest = huge(1.0_dp);  vfb = huge(1.0_dp);  fbb = 0;  vdom = huge(1.0_dp);  vdc = 0
      dch(b)%dnb(i,j,k) = 0;  dch(b)%dnc(:,i,j,k) = 0;  dch(b)%dnr(:,i,j,k) = 0.0_dp;  fbc = 0
      do ob = 1, nblk
         if (mgroup(ob) /= gb) cycle
         call locate_from(b, i, j, k, ob, p, ic, jc, kc, rst, found)
         if (.not. found) cycle
         if (typ(ob)%a(ic,jc,kc) == T_SOLID .or. typ(ob)%a(ic,jc,kc) == T_BLH) then
            dch(b)%flag(i,j,k) = 2
            call nearest_nonsolid(ob, ic, jc, kc, p, found)
            if (.not. found) cycle
            vd = vol(ob)%a(ic,jc,kc)
            if (vd < vfb) then
               vfb = vd;  fbb = ob;  fbc = [ic, jc, kc]
            end if
            cycle
         end if
         vd = vol(ob)%a(ic,jc,kc)
         if (vd < vbest*(1.0_dp - ovst_vtol)) then
            vbest = vd;  dch(b)%dnb(i,j,k) = ob;  dch(b)%dnc(:,i,j,k) = [ic, jc, kc];  dch(b)%dnr(:,i,j,k) = rst
         end if
         if (typ(ob)%a(ic,jc,kc) /= T_FORCED .and. vd < vdom*(1.0_dp - ovst_vtol)) then
            vdom = vd;  vdc = [ob, ic, jc, kc]
         end if
      end do
      dch(b)%vbest(i,j,k) = vbest;  dch(b)%vdom(i,j,k) = vdom;  dch(b)%vfb(i,j,k) = vfb
      dch(b)%fbb(i,j,k) = fbb;  dch(b)%fbc(:,i,j,k) = fbc;  dch(b)%vdc(:,i,j,k) = vdc
    end subroutine group_partial
  end subroutine find_cell_donors

  ! weight of an active cell in global sums: 0.5 for the cells of the
  ! one-cell overlap (counted by both grids), 1 otherwise
  real(dp) function ovs_overlap_weight(b, i, j, k)
    integer, intent(in) :: b, i, j, k
    ovs_overlap_weight = 1.0_dp
    if (.not. allocated(ovl)) return
    if (.not. allocated(ovl(b)%a)) return
    if (ovl(b)%a(i,j,k) == 1) ovs_overlap_weight = 0.5_dp
  end function ovs_overlap_weight

  ! number of no-slip wall patches of block b
  integer function nwall_blk(b)
    integer, intent(in) :: b
    integer :: n
    nwall_blk = 0
    do n = 1, blk(b)%nbc
       if (is_wall(b, n)) nwall_blk = nwall_blk + 1
    end do
  end function nwall_blk

  ! cell centre p of wall block b lies in wall block ob: the cell is
  ! RESERVED for ob when ob's wall is closer to p than b's own wall by
  ! more than Surfovl_tol (the region above an appendage seen from the
  ! hull grid, the mouth of a duct seen from the hull grid); the block
  ! whose wall is nearest owns the cell.  Coincident walls (two patches of
  ! one hull surface) give equal distances -> normal rules.
  logical function reserved_by_other(b, ob, p)
    integer, intent(in) :: b, ob
    real(dp), intent(in) :: p(3)
    integer :: ip, tb, region
    real(dp) :: dist, down, doth, q(3), wbc(3)
    reserved_by_other = .false.
    if (nwall_blk(ob) == 0) return
    down = huge(1.0_dp);  doth = huge(1.0_dp)
    do ip = 1, nwp                                   ! own wall first
       if (wp_blk(ip) /= b) cycle
       call surf_closest(wsurf, p, 3, ip, down, dist, q, tb, wbc, region)
       if (tb > 0) down = min(down, dist)
    end do
    if (down == huge(1.0_dp)) return
    do ip = 1, nwp                                   ! other wall: only within the own distance
       if (wp_blk(ip) /= ob) cycle
       call surf_closest(wsurf, p, 3, ip, min(doth, down), dist, q, tb, wbc, region)
       if (tb > 0) doth = min(doth, dist)
    end do
    if (doth == huge(1.0_dp)) return
    reserved_by_other = (doth < down - surfovl_tol)
  end function reserved_by_other

  ! is the surface point q within Surfovl_tol of a wall patch of block b?
  logical function own_surface_at(b, q)
    integer, intent(in) :: b
    real(dp), intent(in) :: q(3)
    integer :: ip, tb, region
    real(dp) :: dist, qq(3), wbc(3)
    own_surface_at = .false.
    do ip = 1, nwp
       if (wp_blk(ip) /= b) cycle
       call surf_closest(wsurf, q, 3, ip, surfovl_tol*1.0001_dp, dist, qq, tb, wbc, region)
       if (tb > 0) then
          if (dist <= surfovl_tol) then
             own_surface_at = .true.;  return
          end if
       end if
    end do
  end function own_surface_at

  ! cell within Ovst_wall_layers of a no-slip wall patch of its own block:
  ! never dominated by another grid (the wall boundary layer stays with
  ! the grid that carries the wall)
  logical function near_own_wall(b, i, j, k)
    integer, intent(in) :: b, i, j, k
    integer :: n, ax, sd, rr(2,3), idx(3), dist, n3(3)
    near_own_wall = .false.
    idx = [i, j, k]
    do n = 1, blk(b)%nbc
       if (.not. is_wall(b, n)) cycle
       ax = abs(blk(b)%bc_dir(n));  sd = merge(1, 2, blk(b)%bc_dir(n) > 0)
       rr(:,1) = blk(b)%bc_irng(:,n);  rr(:,2) = blk(b)%bc_jrng(:,n);  rr(:,3) = blk(b)%bc_krng(:,n)
       n3 = [blk(b)%nx, blk(b)%ny, blk(b)%nz]
       if (sd == 1) then
          dist = idx(ax) - 1
       else
          dist = n3(ax) - idx(ax)
       end if
       if (dist >= ovst_wall_layers) cycle
       if (any(idx < rr(1,:)) .or. any(idx > max(rr(2,:)-1, rr(1,:)))) cycle
       near_own_wall = .true.;  return
    end do
  end function near_own_wall

  !---------------------------------------------------------------------
  ! Nearest non-solid cell of block ob within +-2 index layers of the
  ! (solid) cell (ic,jc,kc), by distance of the cell centre to p.  On
  ! return (ic,jc,kc) is that cell and found = .true., or found = .false.
  !---------------------------------------------------------------------
  subroutine nearest_nonsolid(ob, ic, jc, kc, p, found)
    integer, intent(in)    :: ob
    integer, intent(inout) :: ic, jc, kc
    real(dp), intent(in)   :: p(3)
    logical, intent(out)   :: found
    integer :: di, dj, dk, ii, jj, kk, best(3)
    real(dp) :: d2, dbest
    found = .false.;  dbest = huge(1.0_dp);  best = 0
    do dk = -FB_R, FB_R
       kk = kc + dk;  if (kk < 1 .or. kk > blk(ob)%nz) cycle
       do dj = -FB_R, FB_R
          jj = jc + dj;  if (jj < 1 .or. jj > blk(ob)%ny) cycle
          do di = -FB_R, FB_R
             ii = ic + di;  if (ii < 1 .or. ii > blk(ob)%nx) cycle
             if (typ(ob)%a(ii,jj,kk) == T_SOLID .or. typ(ob)%a(ii,jj,kk) == T_BLH) cycle
             d2 = sum((cell_center(ob, ii, jj, kk) - p)**2)
             if (d2 < dbest) then
                dbest = d2;  best = [ii, jj, kk];  found = .true.
             end if
          end do
       end do
    end do
    if (found) then
       ic = best(1);  jc = best(2);  kc = best(3)
    end if
  end subroutine nearest_nonsolid

  !---------------------------------------------------------------------
  ! Nearest non-hole cell (cs >= 0) of block ob within +-2 index layers
  ! of cell (ic,jc,kc), by distance of the cell centre to p.
  !---------------------------------------------------------------------
  subroutine nearest_nonhole(ob, ic, jc, kc, p, found)
    integer, intent(in)    :: ob
    integer, intent(inout) :: ic, jc, kc
    real(dp), intent(in)   :: p(3)
    logical, intent(out)   :: found
    integer :: di, dj, dk, ii, jj, kk, best(3)
    real(dp) :: d2, dbest
    found = .false.;  dbest = huge(1.0_dp);  best = 0
    do dk = -FB_R, FB_R
       kk = kc + dk;  if (kk < 1 .or. kk > blk(ob)%nz) cycle
       do dj = -FB_R, FB_R
          jj = jc + dj;  if (jj < 1 .or. jj > blk(ob)%ny) cycle
          do di = -FB_R, FB_R
             ii = ic + di;  if (ii < 1 .or. ii > blk(ob)%nx) cycle
             if (cs(ob)%a(ii,jj,kk) < 0) cycle
             d2 = sum((cell_center(ob, ii, jj, kk) - p)**2)
             if (d2 < dbest) then
                dbest = d2;  best = [ii, jj, kk];  found = .true.
             end if
          end do
       end do
    end do
    if (found) then
       ic = best(1);  jc = best(2);  kc = best(3)
    end if
  end subroutine nearest_nonhole

  !=====================================================================
  ! Donor stencil (8 cells of the donor cell-centre lattice around the
  ! centre of cell (i,j,k) of block b): base index and validity mask.
  !=====================================================================
  subroutine cell_stencil(b, i, j, k, ob, base, okc)
    integer, intent(in)  :: b, i, j, k
    integer, intent(out) :: ob, base(3)
    logical, intent(out) :: okc(8)
    integer :: lin, t, id(3)
    real(dp) :: xi(3), w(8)
    ob = dnb(b)%a(i,j,k)
    okc = .false.;  base = 0
    if (ob == 0) return
    lin = i + blk(b)%nx*((j-1) + blk(b)%ny*(k-1))
    xi = real(dnc(b)%a(:,lin,1), dp) + dnr(b)%a(:,lin,1)
    call lattice_weights(xi, [0.5_dp, 0.5_dp, 0.5_dp], base, w)
    do t = 0, 7
       id = base + [iand(t,1), iand(t/2,1), iand(t/4,1)]
       okc(t+1) = cell_valid(ob, id) .and. w(t+1) > 0.0_dp
    end do
  end subroutine cell_stencil

  ! interior cell: not solid / not hole; ghost cell (one layer): valid
  ! when the block face there is covered (BC or b2b ghost data)
  logical function cell_valid(b, id)
    integer, intent(in) :: b, id(3)
    integer :: a, sd, t1, t2, n(3), nin, ax
    n = [blk(b)%nx, blk(b)%ny, blk(b)%nz]
    cell_valid = .false.
    nin = 0;  ax = 0
    do a = 1, 3
       if (id(a) >= 1 .and. id(a) <= n(a)) then
          nin = nin + 1
       else if (id(a) == 0 .or. id(a) == n(a)+1) then
          ax = a
       else
          return
       end if
    end do
    if (nin == 3) then
       cell_valid = (typ(b)%a(id(1),id(2),id(3)) /= T_SOLID) .and. (typ(b)%a(id(1),id(2),id(3)) /= T_BLH) .and. (cs(b)%a(id(1),id(2),id(3)) >= 0)
    else if (nin == 2) then
       sd = merge(1, 2, id(ax) == 0)
       call tang(ax, t1, t2)
       cell_valid = cov(sd,ax,b)%a(id(t1), id(t2))
    end if
  end function cell_valid

  !=====================================================================
  ! 4b. Classification to a fixed point.
  !=====================================================================
  subroutine classify
    integer :: b, i, j, k, it, nch, ob, base(3), t, id(3), nfr, nact, nhol, nff, nfo
    logical :: okc(8), changed
    integer :: n(3), n_prom, nfd, nfc, n_unc, ic, jc, kc, i2, j2, k2, n_cover, n_cover_f, n_nocov, n_cov_body, lin
    real(dp) :: p(3), rst(3)
    logical :: found, covered

    do b = 1, nblk
       act(b)%a = merge(1, 0, typ(b)%a == T_FREE)
       frg(b)%a = 0;  need(b)%a = 0
    end do
    ! ---- fixed point of the two rules ----------------------------------
    !  A: a dominated (not forced) fringe whose donor stencil contains a
    !     FORCED cell of the donor block becomes active - the boundary
    !     fringe of the other block then interpolates from active cells
    !     (no mutual fringe-from-fringe cycles between two boundaries);
    !  B: every non-active, non-solid donor-stencil member of a fringe
    !     is needed and becomes a fringe itself (never a hole).
    ! Both rules only add active cells / needed cells -> terminates.
    do it = 1, 200
       call mark_fringes
       changed = .false.
       do b = 1, nblk
          do k = 1, blk(b)%nz
           do j = 1, blk(b)%ny
            do i = 1, blk(b)%nx
               if (frg(b)%a(i,j,k) /= 1) cycle
               call cell_stencil(b, i, j, k, ob, base, okc)
               if (ob == 0) cycle
               do t = 0, 7
                  if (.not. okc(t+1)) cycle
                  id = base + [iand(t,1), iand(t/2,1), iand(t/4,1)]
                  if (any(id < 1) .or. id(1) > blk(ob)%nx .or. id(2) > blk(ob)%ny .or. id(3) > blk(ob)%nz) cycle
                  if (act(ob)%a(id(1),id(2),id(3)) == 1) cycle
                  if (typ(ob)%a(id(1),id(2),id(3)) == T_FORCED .and. typ(b)%a(i,j,k) == T_DOM) then
                     act(b)%a(i,j,k) = 1;  changed = .true.;  exit          ! rule A
                  end if
                  if (need(ob)%a(id(1),id(2),id(3)) == 0) then
                     need(ob)%a(id(1),id(2),id(3)) = 1;  changed = .true.    ! rule B
                  end if
               end do
            end do
           end do
          end do
       end do
       if (.not. changed) exit
    end do
    ! ---- rule C (once): a dominated fringe cell next to an active cell of
    !      its own block becomes ACTIVE as well, so that the active regions
    !      of the two grids overlap by one cell (Ovst_active_overlap = 1);
    !      the interpolated band then starts one layer further out and the
    !      rules A / B are iterated again for the new fringes
    if (ovst_active_overlap /= 0) then
       call mark_fringes
       nact = 0
       do b = 1, nblk
          n = [blk(b)%nx, blk(b)%ny, blk(b)%nz]
          do k = 1, blk(b)%nz
           do j = 1, blk(b)%ny
            do i = 1, blk(b)%nx
               if (frg(b)%a(i,j,k) /= 1 .or. typ(b)%a(i,j,k) /= T_DOM) cycle
               if (dnb(b)%a(i,j,k) == 0) cycle
               if (any(act(b)%a(max(1,i-1):min(n(1),i+1), j, k) == 1) .or.          &
                   any(act(b)%a(i, max(1,j-1):min(n(2),j+1), k) == 1) .or.          &
                   any(act(b)%a(i, j, max(1,k-1):min(n(3),k+1)) == 1)) then
                  need(b)%a(i,j,k) = 2                         ! mark, activate below
               end if
            end do
           end do
          end do
       end do
       do b = 1, nblk
          where (need(b)%a == 2)
             act(b)%a = 1;  need(b)%a = 0;  ovl(b)%a = 1
          end where
          nact = nact + count(ovl(b)%a == 1)
       end do
       n_ovl = nact
    end if
    ! ---- rules B + A' (fixed point): every donor-stencil member of a
    !      fringe is needed (B); a DOMINATED (not forced) stencil member
    !      becomes ACTIVE (A', update 12) - a fringe cell then interpolates
    !      from solved cells (or from forced boundary fringes, which
    !      cannot be solved), never from a dominated fringe.  Chains of
    !      dominated fringes interpolating each other across three or more
    !      grids (tube - cap - background of Box_6DoF) formed closed cycles
    !      of interpolation identities; the pressure system is singular on
    !      such a cycle and the solve produced arbitrary values there.  The
    !      6 neighbours of a promoted cell are made needed so that its
    !      stencil never reaches a hole.
    n_prom = 0
    do it = 1, 200
       call mark_fringes
       changed = .false.
       do b = 1, nblk
          do k = 1, blk(b)%nz
           do j = 1, blk(b)%ny
            do i = 1, blk(b)%nx
               if (frg(b)%a(i,j,k) /= 1) cycle
               call cell_stencil(b, i, j, k, ob, base, okc)
               if (ob == 0) cycle
               do t = 0, 7
                  if (.not. okc(t+1)) cycle
                  id = base + [iand(t,1), iand(t/2,1), iand(t/4,1)]
                  if (any(id < 1) .or. id(1) > blk(ob)%nx .or. id(2) > blk(ob)%ny .or. id(3) > blk(ob)%nz) cycle
                  if (act(ob)%a(id(1),id(2),id(3)) == 1) cycle
                  if (typ(ob)%a(id(1),id(2),id(3)) == T_DOM .and. dnb(ob)%a(id(1),id(2),id(3)) /= 0) then
                     act(ob)%a(id(1),id(2),id(3)) = 1;  ovl(ob)%a(id(1),id(2),id(3)) = 1     ! rule A'
                     n_prom = n_prom + 1;  changed = .true.
                     call need_neighbours(ob, id)
                  else if (need(ob)%a(id(1),id(2),id(3)) == 0) then
                     need(ob)%a(id(1),id(2),id(3)) = 1;  changed = .true.                    ! rule B
                  end if
               end do
            end do
           end do
          end do
       end do
       if (.not. changed) exit
    end do
    call mark_fringes
    ! ---- final status -------------------------------------------------
    n_orphan_cells = 0
    do b = 1, nblk
       do k = 1, blk(b)%nz
        do j = 1, blk(b)%ny
         do i = 1, blk(b)%nx
            if (act(b)%a(i,j,k) == 1) then
               cs(b)%a(i,j,k) = 1
            else if (frg(b)%a(i,j,k) == 1) then
               cs(b)%a(i,j,k) = 0
            else if (typ(b)%a(i,j,k) == T_SOLID) then
               cs(b)%a(i,j,k) = -2                      ! hole inside a solid body
            else
               cs(b)%a(i,j,k) = -1                      ! hole (dominated by a finer grid)
            end if
         end do
        end do
       end do
    end do
    n_wallb = 0;  n_gapfill = 0
    ! ---- safety: active cells with a (solid) hole within nf ------------
    do it = 1, ovst_nf + 1
       nch = 0
       do b = 1, nblk
          n = [blk(b)%nx, blk(b)%ny, blk(b)%nz]
          do k = 1, blk(b)%nz
           do j = 1, blk(b)%ny
            do i = 1, blk(b)%nx
               if (cs(b)%a(i,j,k) /= 1) cycle
               if (any(cs(b)%a(max(1,i-ovst_nf):min(n(1),i+ovst_nf),          &
                               max(1,j-ovst_nf):min(n(2),j+ovst_nf),          &
                               max(1,k-ovst_nf):min(n(3),k+ovst_nf)) < 0)) then
                  ! a boundary-layer cell of a wall block (within
                  ! Ovst_wall_layers of its own no-slip wall) next to SOLID
                  ! holes only (cells cut by another body's surface) stays
                  ! ACTIVE: the solid holes act as the wall of that body
                  if (near_own_wall(b, i, j, k)) then
                     if (.not. any(cs(b)%a(max(1,i-ovst_nf):min(n(1),i+ovst_nf),          &
                                           max(1,j-ovst_nf):min(n(2),j+ovst_nf),          &
                                           max(1,k-ovst_nf):min(n(3),k+ovst_nf)) == -1)) cycle
                  end if
                  ! covered = another grid has an ACTIVE (solved) cell in the
                  ! donor stencil of this cell.  A donor stencil made only of
                  ! FORCED boundary fringes of the finer grid (the corner cells
                  ! of a body grid wrapping a box edge, the outer layers of the
                  ! blade grids) does not count: those cells cannot be solved
                  ! and would interpolate back from this very cell (update 12c)
                  if (dnb(b)%a(i,j,k) > 0 .and. stencil_active(b, i, j, k)) then
                     cs(b)%a(i,j,k) = 0;  nch = nch + 1
                     cycle
                  end if
                  ! not covered: the DOMINATED holes (-1) within nf that have a
                  ! donor become fringes (interpolated from the finer grid), so
                  ! that the stencil of this active cell never reads a hole
                  if (dnb(b)%a(i,j,k) > 0) then
                     do k2 = max(1,k-ovst_nf), min(n(3),k+ovst_nf)
                      do j2 = max(1,j-ovst_nf), min(n(2),j+ovst_nf)
                       do i2 = max(1,i-ovst_nf), min(n(1),i+ovst_nf)
                          if (cs(b)%a(i2,j2,k2) == -1 .and. dnb(b)%a(i2,j2,k2) > 0) cs(b)%a(i2,j2,k2) = 0
                       end do
                      end do
                     end do
                  end if
                  if (.not. any(cs(b)%a(max(1,i-ovst_nf):min(n(1),i+ovst_nf),          &
                                        max(1,j-ovst_nf):min(n(2),j+ovst_nf),          &
                                        max(1,k-ovst_nf):min(n(3),k+ovst_nf)) == -1)) then
                     ! not covered by a solved cell of another grid, and the
                     ! holes around are SOLID holes only: it stays active with
                     ! the solid holes acting as the wall of that body (the
                     ! background above a deck, the wedge between a body grid's
                     ! corner and its refinement box); not an orphan.  Solid
                     ! holes are never reopened: no cell inside a body is ever
                     ! solved, the wall is the staircase of hole faces there.
                     n_wallb = n_wallb + 1
                     if (dnb(b)%a(i,j,k) > 0) n_gapfill = n_gapfill + 1
                  else
                     n_orphan_cells = n_orphan_cells + 1;  orph_blk(b) = orph_blk(b) + 1
                  end if
               end if
            end do
           end do
          end do
       end do
       if (nch == 0) exit
    end do
    ! ---- uncovered holes (update 12c): a DOMINATED hole (-1) must lie in
    !      an ACTIVE cell of another grid - domination was decided by the
    !      finest non-forced containing cell, whose own status was not
    !      known at that time (it can itself have become a hole or a fringe
    !      of a third grid: propeller blade grid - refinement - background).
    !      Such a hole is nobody's cell: it becomes ACTIVE here (the finest
    !      covering non-hole cell, if any, is only a fringe), and the safety
    !      rule (active cells next to holes -> fringes with a donor) is
    !      applied to it again; iterated until no uncovered hole is left.
    n_unc = 0
    do it = 1, 10
       nch = 0
       do b = 1, nblk
          n = [blk(b)%nx, blk(b)%ny, blk(b)%nz]
          do k = 1, n(3)
           do j = 1, n(2)
            do i = 1, n(1)
               if (cs(b)%a(i,j,k) /= -1) cycle
               if (typ(b)%a(i,j,k) == T_BLH) cycle
               if (typ(b)%a(i,j,k) == T_FORCED) cycle          ! never solvable: see the coverage pass
               p = cell_center(b, i, j, k)
               covered = .false.
               do ob = 1, nblk
                  call locate_from(b, i, j, k, ob, p, ic, jc, kc, rst, found)
                  if (.not. found) cycle
                  if (cs(ob)%a(ic,jc,kc) == 1) then
                     covered = .true.;  exit
                  end if
                  if (cs(ob)%a(ic,jc,kc) == -2 .or. typ(ob)%a(ic,jc,kc) == T_SOLID .or. typ(ob)%a(ic,jc,kc) == T_BLH) then
                     if (point_in_body(p)) then                 ! a body thinner than this cell: stays a hole
                        covered = .true.;  exit
                     end if
                  end if
               end do
               if (covered) cycle
               cs(b)%a(i,j,k) = 1;  nch = nch + 1;  n_unc = n_unc + 1
            end do
           end do
          end do
       end do
       if (nch == 0) exit
       ! the re-activated cells next to remaining holes: fringes (safety rule)
       do b = 1, nblk
          n = [blk(b)%nx, blk(b)%ny, blk(b)%nz]
          do k = 1, n(3)
           do j = 1, n(2)
            do i = 1, n(1)
               if (cs(b)%a(i,j,k) /= 1) cycle
               if (.not. any(cs(b)%a(max(1,i-ovst_nf):min(n(1),i+ovst_nf),          &
                                     max(1,j-ovst_nf):min(n(2),j+ovst_nf),          &
                                     max(1,k-ovst_nf):min(n(3),k+ovst_nf)) < 0)) cycle
               if (near_own_wall(b, i, j, k)) then
                  if (.not. any(cs(b)%a(max(1,i-ovst_nf):min(n(1),i+ovst_nf),          &
                                        max(1,j-ovst_nf):min(n(2),j+ovst_nf),          &
                                        max(1,k-ovst_nf):min(n(3),k+ovst_nf)) == -1)) cycle
               end if
               if (dnb(b)%a(i,j,k) > 0) cs(b)%a(i,j,k) = 0
            end do
           end do
          end do
       end do
    end do
    ! ---- coverage pass (update 12c): every fluid cell must lie inside at
    !      least one SOLVED cell.  A non-active, non-solid cell whose finest
    !      containing cell (dnb / dnc) is not active and whose donor stencil
    !      has no active corner is covered by nobody (the cut layers of two
    !      body grids meeting at a hub, the corner of a body grid under its
    !      refinement box): it becomes active itself, the dominated holes
    !      within nf around it become fringes (interpolated from the finer
    !      grid) so that its stencil never reads a hole.  Solid holes and
    !      boundary-layer holes are never touched: no cell inside a body is
    !      ever solved, solid holes next to such a cell act as its wall.
    !      A FORCED cell (the layers behind an overset face of its own block)
    !      can never be solved - its stencil would read the block's outer
    !      ghost layer, which carries no boundary condition: when such a cell
    !      is not covered, the finest non-solid, non-forced cell of another
    !      grid containing each of its nodes is made active instead (the
    !      coarse grid is solved under the outer layers of the fine one).
    !      A forced cell whose nodes lie only in solid cells of the other
    !      grids (a blade root inside the hub) is left as it is: the wall of
    !      that body is there.
    n_cover = 0;  n_cover_f = 0;  n_nocov = 0;  n_cov_body = 0
    do b = 1, nblk
       need(b)%a = 0            ! reused here: 1 = examined by the coverage pass (active cells only
    end do                      ! ever increase, so a covered cell stays covered: examine each cell once)
    do it = 1, 50
       nch = 0
       do b = 1, nblk
          n = [blk(b)%nx, blk(b)%ny, blk(b)%nz]
          do k = 1, n(3)
           do j = 1, n(2)
            do i = 1, n(1)
               if (cs(b)%a(i,j,k) == 1 .or. cs(b)%a(i,j,k) == -2) cycle
               if (typ(b)%a(i,j,k) == T_SOLID .or. typ(b)%a(i,j,k) == T_BLH) cycle
               if (need(b)%a(i,j,k) == 1) cycle
               ob = dnb(b)%a(i,j,k)
               if (ob == 0) cycle                                   ! no other grid here at all
               need(b)%a(i,j,k) = 1
               lin = i + n(1)*((j-1) + n(2)*(k-1))
               id = dnc(b)%a(:,lin,1)
               if (cs(ob)%a(id(1),id(2),id(3)) == 1) then          ! centre covered by a solved cell:
                  if (nodes_covered(b, i, j, k, ob, nfd)) cycle    ! ... and every test point of the cell too
               else
                  covered = nodes_covered(b, i, j, k, ob, nfd)
               end if
               ! a body thinner than this cell passes through it (a test
               ! point lies inside a body: the own hole cut missed it because
               ! the centre and the nodes are outside): the finer grid of that
               ! body solves here, the cell stays a hole / fringe - a solved
               ! cell must never reach into a body (update 12c)
               if (nfd > 0) then
                  n_cov_body = n_cov_body + 1;  cycle
               end if
               if (typ(b)%a(i,j,k) == T_FORCED) then
                  n_prom = cover_forced(b, i, j, k, nfd)
                  if (n_prom > 0) nch = nch + n_prom
                  n_cover_f = n_cover_f + n_prom
                  if (it == 1 .and. nfd > 0) n_nocov = n_nocov + 1
                  cycle
               end if
               cs(b)%a(i,j,k) = 1;  ovl(b)%a(i,j,k) = 1;  nch = nch + 1;  n_cover = n_cover + 1
               do k2 = max(1,k-ovst_nf), min(n(3),k+ovst_nf)
                do j2 = max(1,j-ovst_nf), min(n(2),j+ovst_nf)
                 do i2 = max(1,i-ovst_nf), min(n(1),i+ovst_nf)
                    if (cs(b)%a(i2,j2,k2) == -1 .and. dnb(b)%a(i2,j2,k2) > 0) cs(b)%a(i2,j2,k2) = 0
                 end do
                end do
               end do
            end do
           end do
          end do
       end do
       if (nch == 0) exit
    end do
    n_cov_self = n_cover;  n_cov_other = n_cover_f;  n_cov_none = n_nocov;  n_cov_thin = n_cov_body
    if (rank0v .and. n_cov_body > 0) write(*,'(A,I0,A)') '  coverage         : ', n_cov_body, &
       ' uncovered cells crossed by a body thinner than the cell left to the finer grid (never solved)'
    if (rank0v .and. n_cover > 0) write(*,'(A,I0,A)') '  coverage         : ', n_cover, &
       ' fluid cells covered by no solved cell of any grid made active (their dominated holes within nf -> fringes)'
    if (rank0v .and. n_cover_f > 0) write(*,'(A,I0,A)') '  coverage         : ', n_cover_f, &
       ' cells of coarser grids made active under uncovered forced (overset-face) layers of finer grids'
    if (rank0v .and. n_nocov > 0) write(*,'(A,I0,A)') '  coverage         : ', n_nocov, &
       ' forced cells with a node covered by nothing but forced cells of other grids (or outside all grids)'
    ! ---- rule A' on the final status (update 12): the cut-layer fringes
    !      created by the safety step above were not seen by the fixed
    !      point; a dominated fringe cell that serves as a donor of any
    !      fringe cell and has no hole within ovst_nf becomes active
    !      (the same promotion as in the fixed point, on cs) -------------
    do it = 1, 20
       nch = 0
       do b = 1, nblk
          do k = 1, blk(b)%nz
           do j = 1, blk(b)%ny
            do i = 1, blk(b)%nx
               if (cs(b)%a(i,j,k) /= 0) cycle
               call cell_stencil(b, i, j, k, ob, base, okc)
               if (ob == 0) cycle
               n = [blk(ob)%nx, blk(ob)%ny, blk(ob)%nz]
               do t = 0, 7
                  if (.not. okc(t+1)) cycle
                  id = base + [iand(t,1), iand(t/2,1), iand(t/4,1)]
                  if (any(id < 1) .or. any(id > n)) cycle
                  if (cs(ob)%a(id(1),id(2),id(3)) /= 0) cycle
                  if (typ(ob)%a(id(1),id(2),id(3)) /= T_DOM .or. dnb(ob)%a(id(1),id(2),id(3)) == 0) cycle
                  if (any(cs(ob)%a(max(1,id(1)-ovst_nf):min(n(1),id(1)+ovst_nf),          &
                                   max(1,id(2)-ovst_nf):min(n(2),id(2)+ovst_nf),          &
                                   max(1,id(3)-ovst_nf):min(n(3),id(3)+ovst_nf)) < 0)) cycle
                  cs(ob)%a(id(1),id(2),id(3)) = 1;  ovl(ob)%a(id(1),id(2),id(3)) = 1
                  n_prom = n_prom + 1;  nch = nch + 1
               end do
            end do
           end do
          end do
       end do
       if (nch == 0) exit
    end do
    ! ---- statistics: fringes interpolated from fringe / forced cells ----
    nff = 0;  nfo = 0;  nfd = 0;  nfc = 0
    do b = 1, nblk
       do k = 1, blk(b)%nz
        do j = 1, blk(b)%ny
         do i = 1, blk(b)%nx
            if (cs(b)%a(i,j,k) /= 0) cycle
            call cell_stencil(b, i, j, k, ob, base, okc)
            if (ob == 0) cycle
            do t = 0, 7
               if (.not. okc(t+1)) cycle
               id = base + [iand(t,1), iand(t/2,1), iand(t/4,1)]
               if (any(id < 1) .or. id(1) > blk(ob)%nx .or. id(2) > blk(ob)%ny .or. id(3) > blk(ob)%nz) cycle
               if (cs(ob)%a(id(1),id(2),id(3)) == 0) then
                  nff = nff + 1
                  if (typ(ob)%a(id(1),id(2),id(3)) == T_FORCED) nfo = nfo + 1
                  if (typ(ob)%a(id(1),id(2),id(3)) == T_DOM) nfd = nfd + 1
                  if (typ(ob)%a(id(1),id(2),id(3)) == T_FREE) nfc = nfc + 1
                  exit
               end if
            end do
         end do
        end do
       end do
    end do
    n_ff = nff;  n_ffo = nfo
    if (rank0v) then
       do b = 1, nblk
          nact = count(cs(b)%a == 1);  nfr = count(cs(b)%a == 0);  nhol = count(cs(b)%a < 0)
          write(*,'(A,I3,A,I0,A,I0,A,I0,A,I0,A)') '  block ', b, ': active ', nact, ', fringe ', nfr, &
             ', hole ', nhol, ' cells  (of ', blk(b)%nx*blk(b)%ny*blk(b)%nz, ')'
       end do
       write(*,'(A,I0,A)') '  reserved regions : ', n_resv, ' cells of wall blocks inside the passage / region of another wall block -> not solved there'
       write(*,'(A,I0,A)') '  active overlap   : ', n_ovl, ' dominated cells kept active next to the active region (one-cell overlap of the grids)'
       write(*,'(A,I0,A)') '  donor promotion  : ', n_prom, ' dominated donor cells of fringe cells made active (no dominated-fringe donors)'
       write(*,'(A,I0,A)') '  uncovered holes  : ', n_unc, ' dominated holes lying in no active cell of another grid -> active / fringe again'
       write(*,'(A,I0,A,I0,A,I0,A,I0,A)') '  fringe cells with fringe donors: ', nff, ' (donor: boundary fringe ', nfo, &
          ', dominated fringe ', nfd, ', cut-layer fringe ', nfc, ')'
       if (n_wallb > 0) write(*,'(A,I0,A,I0,A)') '  solid-bounded    : ', n_wallb, &
          ' active cells next to solid holes with no solved cell of another grid there (the holes act as the wall; ', &
          n_gapfill, ' of them covered only by forced fringes of a finer grid)'
       if (n_orphan_cells > 0) write(*,'(A,I0,A,7I8)') '  WARNING: ', n_orphan_cells, &
          ' active cells next to holes have no donor (orphans); per block:', orph_blk(1:min(nblk,7))
    end if

  contains
    ! fringe = non-active non-solid cell with an active cell within nf,
    ! or a needed donor-stencil member
    subroutine mark_fringes
      integer :: bb, ii, jj, kk, nn(3)
      do bb = 1, nblk
         nn = [blk(bb)%nx, blk(bb)%ny, blk(bb)%nz]
         frg(bb)%a = 0
         do kk = 1, nn(3)
          do jj = 1, nn(2)
           do ii = 1, nn(1)
              if (act(bb)%a(ii,jj,kk) == 1 .or. typ(bb)%a(ii,jj,kk) == T_SOLID) cycle
              if (need(bb)%a(ii,jj,kk) == 1) then
                 frg(bb)%a(ii,jj,kk) = 1;  cycle
              end if
              if (any(act(bb)%a(max(1,ii-ovst_nf):min(nn(1),ii+ovst_nf),        &
                                max(1,jj-ovst_nf):min(nn(2),jj+ovst_nf),        &
                                max(1,kk-ovst_nf):min(nn(3),kk+ovst_nf)) == 1)) &
                 frg(bb)%a(ii,jj,kk) = 1
           end do
          end do
         end do
      end do
    end subroutine mark_fringes
    ! .true. when all 8 nodes of cell (i,j,k) of block b lie inside ACTIVE
    ! cells of other blocks (the finest block ob first, then any other):
    ! a coarse cell straddling the forced boundary layers of a finer block
    ! is only partly covered by solved cells even when its centre is
    ! (NSWC15E bow: background cells 3x larger than the refinement's
    ! two forced layers) - such a cell stays solved itself
    ! (the 8 nodes, the 12 edge midpoints and the 6 face centres are
    ! tested - 26 points of the trilinear cell: a coarse cell whose nodes
    ! all lie in solved fine cells can still hold a sliver of the fine
    ! grid's fringe layer between two of its nodes)
    ! nin_ = number of test points that lie inside a body (in a solid cell
    ! of another block and with winding number > 1/2): such a cell is
    ! crossed by a body thinner than itself and must not be solved
    logical function nodes_covered(bb_, i_, j_, k_, ob_, nin_)
      integer, intent(in)  :: bb_, i_, j_, k_, ob_
      integer, intent(out) :: nin_
      integer :: t_, ic_, jc_, kc_, m_, iu_, iv_, iw_, q_
      real(dp) :: pn_(3), r_(3), xn_(3,8), u_, v_, w_, wt_
      logical :: fnd_, sol_
      nodes_covered = .true.;  nin_ = 0
      do t_ = 0, 7
         xn_(:,t_+1) = [blk(bb_)%xn(i_+iand(t_,1), j_+iand(t_/2,1), k_+iand(t_/4,1)), &
                        blk(bb_)%yn(i_+iand(t_,1), j_+iand(t_/2,1), k_+iand(t_/4,1)), &
                        blk(bb_)%zn(i_+iand(t_,1), j_+iand(t_/2,1), k_+iand(t_/4,1))]
      end do
      do t_ = 0, 26
         iu_ = mod(t_,3);  iv_ = mod(t_/3,3);  iw_ = t_/9
         if (iu_ == 1 .and. iv_ == 1 .and. iw_ == 1) cycle          ! the centre is tested by the caller
         u_ = 0.5_dp*iu_;  v_ = 0.5_dp*iv_;  w_ = 0.5_dp*iw_
         pn_ = 0.0_dp
         do q_ = 0, 7
            wt_ = merge(u_, 1.0_dp-u_, iand(q_,1) == 1)*merge(v_, 1.0_dp-v_, iand(q_/2,1) == 1)* &
                  merge(w_, 1.0_dp-w_, iand(q_/4,1) == 1)
            pn_ = pn_ + wt_*xn_(:,q_+1)
         end do
         sol_ = .false.
         call locate_from(bb_, i_, j_, k_, ob_, pn_, ic_, jc_, kc_, r_, fnd_)
         if (fnd_) then
            if (cs(ob_)%a(ic_,jc_,kc_) == 1) cycle
            if (cs(ob_)%a(ic_,jc_,kc_) == -2 .or. typ(ob_)%a(ic_,jc_,kc_) == T_SOLID .or. &
                typ(ob_)%a(ic_,jc_,kc_) == T_BLH) sol_ = .true.
         end if
         fnd_ = .false.
         do m_ = 1, nblk
            if (m_ == ob_) cycle
            call locate_from(bb_, i_, j_, k_, m_, pn_, ic_, jc_, kc_, r_, fnd_)
            if (.not. fnd_) cycle
            if (cs(m_)%a(ic_,jc_,kc_) == 1) exit
            if (cs(m_)%a(ic_,jc_,kc_) == -2 .or. typ(m_)%a(ic_,jc_,kc_) == T_SOLID .or. &
                typ(m_)%a(ic_,jc_,kc_) == T_BLH) sol_ = .true.
            fnd_ = .false.
         end do
         if (.not. fnd_) then
            nodes_covered = .false.
            if (sol_) then
               if (point_in_body(pn_)) nin_ = nin_ + 1
            end if
         end if
      end do
    end function nodes_covered
    ! is point p inside one of the closed bodies (winding number > 1/2)?
    logical function point_in_body(p_)
      real(dp), intent(in) :: p_(3)
      integer :: ib_
      point_in_body = .false.
      if (.not. hsurf%wind_ready) return
      do ib_ = 1, nbody
         if (.not. body_closed(ib_)) cycle
         if (surf_wind(hsurf, p_, ib_) > 0.5_dp + 1.0e-2_dp) then
            point_in_body = .true.;  return
         end if
      end do
    end function point_in_body
    ! cover a FORCED cell (i,j,k) of block b by the other grids: for each of
    ! its 27 test points (nodes, edge midpoints, face centres, centre) that does not lie in an ACTIVE cell of
    ! another block, the finest containing cell of another block that can
    ! be solved (not solid, not a boundary-layer hole, not forced itself)
    ! is made active (overlap cell), and the dominated holes within nf
    ! around that cell become fringes.  Returns the number of cells made
    ! active (0: every node already covered, or only solid / forced cells
    ! of the other grids around).  nbad_ = number of nodes that lie neither
    ! in an active nor in a solid cell of another grid and for which no
    ! cell could be promoted (only forced cells of other grids there, or
    ! outside every other grid): a place solved by nobody.
    integer function cover_forced(bb_, i_, j_, k_, nbad_)
      integer, intent(in)  :: bb_, i_, j_, k_
      integer, intent(out) :: nbad_
      integer :: t_, m_, ic_, jc_, kc_, pb_, pc_(3), nn_(3), i2_, j2_, k2_, iu_, iv_, iw_, q_
      real(dp) :: pn_(3), r_(3), vb_, xn_(3,8), u_, v_, w_, wt_
      logical :: fnd_, act_, sol_
      cover_forced = 0;  nbad_ = 0
      do t_ = 0, 7
         xn_(:,t_+1) = [blk(bb_)%xn(i_+iand(t_,1), j_+iand(t_/2,1), k_+iand(t_/4,1)), &
                        blk(bb_)%yn(i_+iand(t_,1), j_+iand(t_/2,1), k_+iand(t_/4,1)), &
                        blk(bb_)%zn(i_+iand(t_,1), j_+iand(t_/2,1), k_+iand(t_/4,1))]
      end do
      do t_ = 0, 26                                   ! nodes, edge midpoints, face centres, centre
         iu_ = mod(t_,3);  iv_ = mod(t_/3,3);  iw_ = t_/9
         u_ = 0.5_dp*iu_;  v_ = 0.5_dp*iv_;  w_ = 0.5_dp*iw_
         pn_ = 0.0_dp
         do q_ = 0, 7
            wt_ = merge(u_, 1.0_dp-u_, iand(q_,1) == 1)*merge(v_, 1.0_dp-v_, iand(q_/2,1) == 1)* &
                  merge(w_, 1.0_dp-w_, iand(q_/4,1) == 1)
            pn_ = pn_ + wt_*xn_(:,q_+1)
         end do
         act_ = .false.;  sol_ = .false.;  pb_ = 0;  vb_ = huge(1.0_dp);  pc_ = 0
         do m_ = 1, nblk
            call locate_from(bb_, i_, j_, k_, m_, pn_, ic_, jc_, kc_, r_, fnd_)
            if (.not. fnd_) cycle
            if (cs(m_)%a(ic_,jc_,kc_) == 1) then
               act_ = .true.;  exit
            end if
            if (cs(m_)%a(ic_,jc_,kc_) == -2 .or. typ(m_)%a(ic_,jc_,kc_) == T_SOLID .or. &
                typ(m_)%a(ic_,jc_,kc_) == T_BLH) then
               sol_ = .true.;  cycle
            end if
            if (typ(m_)%a(ic_,jc_,kc_) == T_FORCED) cycle
            if (vol(m_)%a(ic_,jc_,kc_) < vb_) then
               vb_ = vol(m_)%a(ic_,jc_,kc_);  pb_ = m_;  pc_ = [ic_, jc_, kc_]
            end if
         end do
         if (act_) cycle
         if (pb_ == 0) then
            if (.not. sol_) nbad_ = nbad_ + 1
            cycle
         end if
         cs(pb_)%a(pc_(1),pc_(2),pc_(3)) = 1;  ovl(pb_)%a(pc_(1),pc_(2),pc_(3)) = 1
         cover_forced = cover_forced + 1
         nn_ = [blk(pb_)%nx, blk(pb_)%ny, blk(pb_)%nz]
         do k2_ = max(1,pc_(3)-ovst_nf), min(nn_(3),pc_(3)+ovst_nf)
          do j2_ = max(1,pc_(2)-ovst_nf), min(nn_(2),pc_(2)+ovst_nf)
           do i2_ = max(1,pc_(1)-ovst_nf), min(nn_(1),pc_(1)+ovst_nf)
              if (cs(pb_)%a(i2_,j2_,k2_) == -1 .and. dnb(pb_)%a(i2_,j2_,k2_) > 0) cs(pb_)%a(i2_,j2_,k2_) = 0
           end do
          end do
         end do
      end do
    end function cover_forced
    ! .true. when the donor stencil of cell (i,j,k) of block b contains at
    ! least one ACTIVE cell of the donor block (final status cs)
    logical function stencil_active(bb_, i_, j_, k_)
      integer, intent(in) :: bb_, i_, j_, k_
      integer :: ob_, base_(3), t_, id_(3)
      logical :: okc_(8)
      stencil_active = .false.
      call cell_stencil(bb_, i_, j_, k_, ob_, base_, okc_)
      if (ob_ == 0) return
      do t_ = 0, 7
         if (.not. okc_(t_+1)) cycle
         id_ = base_ + [iand(t_,1), iand(t_/2,1), iand(t_/4,1)]
         if (any(id_ < 1) .or. id_(1) > blk(ob_)%nx .or. id_(2) > blk(ob_)%ny .or. id_(3) > blk(ob_)%nz) cycle
         if (cs(ob_)%a(id_(1),id_(2),id_(3)) == 1) then
            stencil_active = .true.;  return
         end if
      end do
    end function stencil_active
    ! the 6 neighbours of a promoted (active) cell must not be holes:
    ! mark the non-active, non-solid ones as needed (they become fringes)
    subroutine need_neighbours(bb_, c0)
      integer, intent(in) :: bb_, c0(3)
      integer :: dd_, c1(3), nn_(3)
      nn_ = [blk(bb_)%nx, blk(bb_)%ny, blk(bb_)%nz]
      do dd_ = 1, 6
         c1 = c0;  c1((dd_+1)/2) = c1((dd_+1)/2) + merge(1, -1, mod(dd_,2) == 1)
         if (any(c1 < 1) .or. any(c1 > nn_)) cycle
         if (act(bb_)%a(c1(1),c1(2),c1(3)) == 1) cycle
         if (typ(bb_)%a(c1(1),c1(2),c1(3)) == T_SOLID .or. typ(bb_)%a(c1(1),c1(2),c1(3)) == T_BLH) cycle
         if (need(bb_)%a(c1(1),c1(2),c1(3)) == 0) need(bb_)%a(c1(1),c1(2),c1(3)) = 1
      end do
    end subroutine need_neighbours
  end subroutine classify

  !=====================================================================
  ! Ghost-extended cell status per block (halo / b2b through the patch
  ! maps, physical and overset boundaries mirrored) and face statuses.
  !=====================================================================
  subroutine ghost_status
    integer :: b, i, j, k, m, a, na, g, w1(3), w2(3), r(3), dcell(3), n(3), q, e(3), sa, sb
    allocate(gcs(nblk), gf(3,nblk))
    do b = 1, nblk
       n = [blk(b)%nx, blk(b)%ny, blk(b)%nz]
       allocate(gcs(b)%a(1-ng:n(1)+ng, 1-ng:n(2)+ng, 1-ng:n(3)+ng))
       gcs(b)%a = -1
       gcs(b)%a(1:n(1),1:n(2),1:n(3)) = cs(b)%a
       do g = 1, ng
          gcs(b)%a(1-g,:,:)      = gcs(b)%a(2-g,:,:)
          gcs(b)%a(n(1)+g,:,:)   = gcs(b)%a(n(1)+g-1,:,:)
       end do
       do g = 1, ng
          gcs(b)%a(:,1-g,:)      = gcs(b)%a(:,2-g,:)
          gcs(b)%a(:,n(2)+g,:)   = gcs(b)%a(:,n(2)+g-1,:)
       end do
       do g = 1, ng
          gcs(b)%a(:,:,1-g)      = gcs(b)%a(:,:,2-g)
          gcs(b)%a(:,:,n(3)+g)   = gcs(b)%a(:,:,n(3)+g-1)
       end do
    end do
    ! b2b / periodic ghosts from the donor block (cell index map)
    do m = 1, nxp
       na = xp(m)%na
       do a = 1, 3
          if (a == na) cycle
          w1(a) = max(xp(m)%rn(1,a) - ng, 1-ng)
          w2(a) = min(xp(m)%rn(2,a) - 1 + ng, ncell(xp(m)%rblk,a)+ng)
       end do
       if (xp(m)%rside == 1) then
          w1(na) = 1-ng;  w2(na) = 0
       else
          w1(na) = ncell(xp(m)%rblk,na)+1;  w2(na) = ncell(xp(m)%rblk,na)+ng
       end if
       do k = w1(3), w2(3)
        do j = w1(2), w2(2)
         do i = w1(1), w2(1)
            r = [i, j, k]
            do a = 1, 3
               dcell(xp(m)%q(a)) = xp(m)%sgn(a)*r(a) + xp(m)%coff(a)
            end do
            if (any(dcell < 1)) cycle
            if (dcell(1) > blk(xp(m)%dblk)%nx .or. dcell(2) > blk(xp(m)%dblk)%ny .or. &
                dcell(3) > blk(xp(m)%dblk)%nz) cycle
            gcs(xp(m)%rblk)%a(i,j,k) = cs(xp(m)%dblk)%a(dcell(1),dcell(2),dcell(3))
         end do
        end do
       end do
    end do
    ! face statuses: face f of family q between cells f and f+1
    do b = 1, nblk
       n = [blk(b)%nx, blk(b)%ny, blk(b)%nz]
       do q = 1, 3
          e = 0;  e(q) = 1
          allocate(gf(q,b)%a(1-ng:n(1)+ng-e(1), 1-ng:n(2)+ng-e(2), 1-ng:n(3)+ng-e(3)))
          do k = 1-ng, n(3)+ng-e(3)
           do j = 1-ng, n(2)+ng-e(2)
            do i = 1-ng, n(1)+ng-e(1)
               sa = gcs(b)%a(i,j,k);  sb = gcs(b)%a(i+e(1),j+e(2),k+e(3))
               if ((sa == 1 .or. sb == 1) .and. sa >= 0 .and. sb >= 0) then
                  gf(q,b)%a(i,j,k) = 1
               else if (sa == 0 .or. sb == 0) then
                  gf(q,b)%a(i,j,k) = 0
               else
                  gf(q,b)%a(i,j,k) = -1
               end if
            end do
           end do
          end do
       end do
    end do
  end subroutine ghost_status

  !=====================================================================
  ! 5. Fringe entries: cells and faces with their donor stencils.
  !=====================================================================
  subroutine build_entries
    integer :: b, i, j, k, q, e(3), cap, ob, base(3), t, id(3), lin, ic, jc, kc, dd, l, idx3(3), t1, t2
    integer :: n(3), ntot, m, bs(3,3)
    real(dp) :: xi(3), p(3), rst(3), wsum, wt(8), r(3), ws(8,3), off(3)
    logical :: okc(8), found, ok
    integer, allocatable :: key(:), ord(:), cnt(:)
    integer,  allocatable :: t_kind(:), t_rd(:), t_dd(:), t_ridx(:,:), t_base(:,:,:), t_dcel(:,:)
    real(dp), allocatable :: t_box(:,:)
    integer  :: cn, gi, gj, gk, dc(3)
    real(dp) :: pn(3), rn(3), bmn(3), bmx(3)
    logical  :: fnd
    real(dp), allocatable :: t_w(:,:,:)
    ! orphan fill (update 12): fringe cells / faces without any donor take
    ! the value of a neighbouring active cell / face (or of a neighbour
    ! that has a regular donor) of their OWN block - a zero-gradient
    ! closure instead of an undefined value (p = 0, stale velocity),
    ! which blew up the momentum sweeps at the box edges of Box_6DoF
    integer, allocatable :: oc_list(:,:), of_list(:,:)
    type(iarr_t), allocatable :: hasc(:)
    type(iarr_t), allocatable :: hasf(:,:)
    integer :: noc, nof, npass, nfill, nb3(3), c, d2, cand(3), tq1, tq2
    logical :: done

    cap = 1024;  nent = 0
    noc = 0;  nof = 0
    allocate(oc_list(4, 1024), of_list(5, 1024))
    allocate(hasc(nblk), hasf(3, nblk))
    do b = 1, nblk
       n = [blk(b)%nx, blk(b)%ny, blk(b)%nz]
       allocate(hasc(b)%a(n(1), n(2), n(3)));  hasc(b)%a = 0
       do q = 1, 3
          e = 0;  e(q) = 1
          allocate(hasf(q,b)%a(1-e(1):n(1), 1-e(2):n(2), 1-e(3):n(3)));  hasf(q,b)%a = 0
       end do
    end do
    allocate(t_kind(cap), t_rd(cap), t_dd(cap), t_ridx(3,cap), t_base(3,3,cap), t_w(8,3,cap), t_dcel(3,cap), t_box(6,cap))
    t_dcel = 0;  t_box = 0.0_dp
    n_orphan_faces = 0;  n_fb_faces = 0;  n_fb_used = 0;  nwf = 0;  n_bcf = 0
    if (allocated(wf_list)) deallocate(wf_list);  allocate(wf_list(5, 1024))
    do b = 1, nblk
       n = [blk(b)%nx, blk(b)%ny, blk(b)%nz]
       ! ---- fringe cells --------------------------------------------------
       do k = 1, n(3)
        do j = 1, n(2)
         do i = 1, n(1)
            if (cs(b)%a(i,j,k) /= 0) cycle
            call cell_stencil(b, i, j, k, ob, base, okc)
            if (ob == 0) then
               call add_oc(b, i, j, k);  cycle
            end if
            if (dfb(b)%a(i,j,k) == 1) n_fb_used = n_fb_used + 1
            lin = i + n(1)*((j-1) + n(2)*(k-1))
            xi = real(dnc(b)%a(:,lin,1), dp) + dnr(b)%a(:,lin,1)
            call lattice_weights(xi, [0.5_dp, 0.5_dp, 0.5_dp], base, wt)
            do t = 1, 8
               if (.not. okc(t)) wt(t) = 0.0_dp
            end do
            wsum = sum(wt)
            if (wsum <= 1.0e-12_dp) then
               call add_oc(b, i, j, k);  cycle
            end if
            wt = wt / wsum
            dd = dom_of_cell(ob, dnc(b)%a(:,lin,1))
            do l = 1, 3
               bs(:,l) = base;  ws(:,l) = wt
            end do
            call push(0, dom_of_cell(b, [i,j,k]), [i,j,k], dd, bs, ws)
            hasc(b)%a(i,j,k) = 2                           ! regular donor
            ! geometric-VoF payload: the receiver cell's node box in the
            ! natural coordinates of the donor cell that contains its centre
            ! (only for regular, non-fallback donors)
            if (dfb(b)%a(i,j,k) == 0) then
               dc = dnc(b)%a(:,lin,1);  bmn = huge(1.0_dp);  bmx = -huge(1.0_dp);  fnd = .true.
               do cn = 0, 7
                  gi = i + iand(cn,1);  gj = j + iand(cn/2,1);  gk = k + iand(cn/4,1)
                  pn = [blk(b)%xn(gi,gj,gk), blk(b)%yn(gi,gj,gk), blk(b)%zn(gi,gj,gk)]
                  call locate_from(b, i, j, k, ob, pn, ic, jc, kc, rn, found)
                  if (.not. found) then
                     fnd = .false.;  exit
                  end if
                  rn = rn + real([ic, jc, kc] - dc, dp)
                  bmn = min(bmn, rn);  bmx = max(bmx, rn)
               end do
               if (fnd) then
                  t_dcel(:,nent) = dc;  t_box(1:3,nent) = bmn;  t_box(4:6,nent) = bmx
               end if
            end if
         end do
        end do
       end do
       ! ---- fringe faces --------------------------------------------------
       do q = 1, 3
          e = 0;  e(q) = 1
          do k = 1-e(3), n(3)
           do j = 1-e(2), n(2)
            do i = 1-e(1), n(1)
               if (gf(q,b)%a(i,j,k) /= 0) cycle
               ! a fringe face lying ON a covered block face (physical boundary
               ! patch, block-to-block window) gets its value from that boundary
               ! condition / exchange: no donor needed
               idx3 = [i, j, k]
               if (idx3(q) == 0 .or. idx3(q) == n(q)) then
                  call tang(q, t1, t2)
                  if (idx3(t1) >= 1 .and. idx3(t1) <= n(t1) .and. idx3(t2) >= 1 .and. idx3(t2) <= n(t2)) then
                     if (cov(merge(1, 2, idx3(q) == 0), q, b)%a(idx3(t1), idx3(t2))) then
                        n_bcf = n_bcf + 1;  cycle
                     end if
                  end if
               end if
               p = face_center(b, q, i, j, k)
               ob = 0;  rst = 0.0_dp;  wsum = huge(1.0_dp);  id = 0
               do m = 1, nblk
                  call locate_from(b, i, j, k, m, p, ic, jc, kc, r, found)
                  if (.not. found) cycle
                  if (cs(m)%a(ic,jc,kc) < 0) cycle
                  if (vol(m)%a(ic,jc,kc) < wsum*(1.0_dp - ovst_vtol)) then
                     wsum = vol(m)%a(ic,jc,kc);  ob = m;  id = [ic, jc, kc];  rst = r
                  end if
               end do
               ok = .false.
               if (ob > 0) then
                  xi = real(id, dp) + rst
                  call face_weights(ob, xi, bs, ws, ok)
               end if
               if (.not. ok) then
                  ! FALLBACK: the face centre lies in no non-hole cell of any
                  ! other block, or the staggered donor lattices there are
                  ! invalid (hole faces) -> take the nearest non-hole cell of
                  ! the finest block within 2 index layers and interpolate at
                  ! its centre (the two faces of that cell, 0th order)
                  ob = 0;  wsum = huge(1.0_dp)
                  do m = 1, nblk
                     call locate_from(b, i, j, k, m, p, ic, jc, kc, r, found)
                     if (.not. found) cycle
                     call nearest_nonhole(m, ic, jc, kc, p, found)
                     if (.not. found) cycle
                     if (vol(m)%a(ic,jc,kc) < wsum) then
                        wsum = vol(m)%a(ic,jc,kc);  ob = m;  id = [ic, jc, kc]
                     end if
                  end do
                  if (ob > 0) then
                     xi = real(id, dp) + 0.5_dp
                     call face_weights(ob, xi, bs, ws, ok)
                     if (ok) n_fb_faces = n_fb_faces + 1
                  end if
               end if
               if (.not. ok) then
                  ! no donor at all: a fringe face next to a SOLID hole is a
                  ! wall face of that body (zero velocity) - not an orphan;
                  ! it is turned into a hole face after the donor search
                  if (gcs(b)%a(i,j,k) == -2 .or. gcs(b)%a(i+e(1),j+e(2),k+e(3)) == -2) then
                     nwf = nwf + 1
                     if (nwf > size(wf_list, 2)) call grow_wf()
                     wf_list(:,nwf) = [q, b, i, j, k]
                  else
                     nof = nof + 1
                     if (nof > size(of_list, 2)) call grow_lists()
                     of_list(:,nof) = [q, b, i, j, k]
                     if (poi_verbose .and. nof <= 40) write(*,'(A,5I6,A,3F11.6,A,2I3)') &
                        '    orphan face (q,blk,i,j,k)', q, b, i, j, k, ' xyz', p, ' cells', &
                        gcs(b)%a(i,j,k), gcs(b)%a(i+e(1),j+e(2),k+e(3))
                  end if
               else
                  dd = dom_of_cell(ob, id)
                  call push(q, dom_of_face(b, q, [i,j,k]), [i,j,k], dd, bs, ws)
                  hasf(q,b)%a(i,j,k) = 2                   ! regular donor
               end if
            end do
           end do
          end do
       end do
    end do

    ! ---- orphan fill: zero-gradient closure from a neighbour of the same
    !      block (active first, then a neighbour with a regular donor, then
    !      a neighbour filled in an earlier pass); iterated so that chains
    !      end in a solved value ---------------------------------------------
    n_orphan_faces = 0;  n_fill_cells = 0;  n_fill_faces = 0
    do npass = 1, 8
       nfill = 0
       do m = 1, noc
          if (oc_list(1,m) == 0) cycle                     ! already filled
          b = oc_list(1,m);  i = oc_list(2,m);  j = oc_list(3,m);  k = oc_list(4,m)
          n = [blk(b)%nx, blk(b)%ny, blk(b)%nz]
          done = .false.
          do c = 1, 3                                      ! 1: active, 2: regular donor, 3: filled earlier
             do d2 = 1, 6
                cand = [i, j, k];  cand((d2+1)/2) = cand((d2+1)/2) + merge(1, -1, mod(d2,2) == 1)
                if (any(cand < 1) .or. any(cand > n)) cycle
                if (c == 1) then
                   if (cs(b)%a(cand(1),cand(2),cand(3)) /= 1) cycle
                else
                   if (hasc(b)%a(cand(1),cand(2),cand(3)) /= c) cycle
                end if
                done = .true.;  exit
             end do
             if (done) exit
          end do
          if (.not. done) cycle
          ws = 0.0_dp;  ws(1,:) = 1.0_dp
          do l = 1, 3
             bs(:,l) = cand
          end do
          call push(0, dom_of_cell(b, [i,j,k]), [i,j,k], dom_of_cell(b, cand), bs, ws)
          hasc(b)%a(i,j,k) = 3;  oc_list(1,m) = 0;  nfill = nfill + 1;  n_fill_cells = n_fill_cells + 1
       end do
       do m = 1, nof
          if (of_list(1,m) == 0) cycle
          q = of_list(1,m);  b = of_list(2,m);  i = of_list(3,m);  j = of_list(4,m);  k = of_list(5,m)
          n = [blk(b)%nx, blk(b)%ny, blk(b)%nz]
          e = 0;  e(q) = 1
          call tang(q, tq1, tq2)
          done = .false.
          do c = 1, 3
             do d2 = 1, 6
                cand = [i, j, k]
                select case ((d2+1)/2)                     ! along q first, then the tangential axes
                case (1);  cand(q)   = cand(q)   + merge(1, -1, mod(d2,2) == 1)
                case (2);  cand(tq1) = cand(tq1) + merge(1, -1, mod(d2,2) == 1)
                case (3);  cand(tq2) = cand(tq2) + merge(1, -1, mod(d2,2) == 1)
                end select
                if (cand(q) < 0 .or. cand(q) > n(q)) cycle
                if (cand(tq1) < 1 .or. cand(tq1) > n(tq1) .or. cand(tq2) < 1 .or. cand(tq2) > n(tq2)) cycle
                if (c == 1) then
                   if (gf(q,b)%a(cand(1),cand(2),cand(3)) /= 1) cycle
                else
                   if (hasf(q,b)%a(cand(1),cand(2),cand(3)) /= c) cycle
                end if
                done = .true.;  exit
             end do
             if (done) exit
          end do
          if (.not. done) cycle
          ! the face value is taken from the two cells of the neighbour face
          ! through the staggered lattice of the same block: weight 1 on the
          ! face itself = 0.5 / 0.5 on its two cells for the q-component and
          ! the cell-centred values for the others is what face_weights does
          ! for a point ON the face; copy the face entry instead (0th order)
          xi = real(cand, dp) + 0.5_dp;  xi(q) = real(cand(q), dp) + 1.0_dp
          call face_weights(b, xi, bs, ws, ok)
          if (.not. ok) cycle
          nb3 = cand;  if (nb3(q) < 1) nb3 = nb3 + e
          call push(q, dom_of_face(b, q, [i,j,k]), [i,j,k], dom_of_cell(b, nb3), bs, ws)
          hasf(q,b)%a(i,j,k) = 3;  of_list(1,m) = 0;  nfill = nfill + 1;  n_fill_faces = n_fill_faces + 1
       end do
       if (nfill == 0) exit
    end do
    do m = 1, noc
       if (oc_list(1,m) /= 0) n_orphan_cells = n_orphan_cells + 1
    end do
    do m = 1, nof
       if (of_list(1,m) /= 0) n_orphan_faces = n_orphan_faces + 1
    end do
    deallocate(oc_list, of_list)
    do b = 1, nblk
       deallocate(hasc(b)%a)
       do q = 1, 3
          deallocate(hasf(q,b)%a)
       end do
    end do
    deallocate(hasc, hasf)

    ! ---- sort by (sender rank, receiver rank) --------------------------
    allocate(key(max(nent,1)), ord(max(nent,1)), cnt(nprocs*nprocs+1))
    cnt = 0
    do m = 1, nent
       key(m) = dom(t_dd(m))%owner*nprocs + dom(t_rd(m))%owner + 1
       cnt(key(m)+1) = cnt(key(m)+1) + 1
    end do
    cnt(1) = 1
    do m = 2, nprocs*nprocs+1
       cnt(m) = cnt(m) + cnt(m-1)
    end do
    do m = 1, nent
       ord(cnt(key(m))) = m;  cnt(key(m)) = cnt(key(m)) + 1
    end do
    allocate(e_kind(max(nent,1)), e_rd(max(nent,1)), e_dd(max(nent,1)), e_ridx(3,max(nent,1)), &
             e_base(3,3,max(nent,1)), e_w(8,3,max(nent,1)), ovs_val(3,max(nent,1)))
    if (allocated(e_dcel)) deallocate(e_dcel, e_box)
    allocate(e_dcel(3,max(nent,1)), e_box(6,max(nent,1)));  e_dcel = 0;  e_box = 0.0_dp
    do m = 1, nent
       e_dcel(:,m) = t_dcel(:,ord(m));  e_box(:,m) = t_box(:,ord(m))
    end do
    do m = 1, nent
       e_kind(m) = t_kind(ord(m));  e_rd(m) = t_rd(ord(m));  e_dd(m) = t_dd(ord(m))
       e_ridx(:,m) = t_ridx(:,ord(m));  e_base(:,:,m) = t_base(:,:,ord(m));  e_w(:,:,m) = t_w(:,:,ord(m))
    end do
    ovs_val = 0.0_dp
    deallocate(t_kind, t_rd, t_dd, t_ridx, t_base, t_w, t_dcel, t_box, key, ord, cnt)
    ! ---- rank pairs -----------------------------------------------------
    npair = 0
    do m = 1, nent
       if (m == 1) then
          npair = 1
       else if (dom(e_dd(m))%owner /= dom(e_dd(m-1))%owner .or. dom(e_rd(m))%owner /= dom(e_rd(m-1))%owner) then
          npair = npair + 1
       end if
    end do
    allocate(p_s(max(npair,1)), p_r(max(npair,1)), p_e1(max(npair,1)), p_e2(max(npair,1)), &
             p_nc(max(npair,1)), p_nfc(max(npair,1)))
    p_nc = 0;  p_nfc = 0;  p_e1 = 1;  p_e2 = 0
    npair = 0
    do m = 1, nent
       if (m == 1) then
          npair = 1;  p_e1(1) = 1
       else if (dom(e_dd(m))%owner /= dom(e_dd(m-1))%owner .or. dom(e_rd(m))%owner /= dom(e_rd(m-1))%owner) then
          p_e2(npair) = m - 1
          npair = npair + 1;  p_e1(npair) = m
       end if
       p_s(npair) = dom(e_dd(m))%owner;  p_r(npair) = dom(e_rd(m))%owner
       if (e_kind(m) == 0) then
          p_nc(npair) = p_nc(npair) + 1
       else
          p_nfc(npair) = p_nfc(npair) + 1
       end if
    end do
    if (nent > 0) p_e2(npair) = nent
    ! ---- per-domain receiver lists --------------------------------------
    allocate(dcel(ndom), dfac(ndom))
    do m = 1, nent
       if (e_kind(m) == 0) then
          dcel(e_rd(m))%n = dcel(e_rd(m))%n + 1
       else
          dfac(e_rd(m))%n = dfac(e_rd(m))%n + 1
       end if
    end do
    do dd = 1, ndom
       allocate(dcel(dd)%e(max(dcel(dd)%n,1)), dfac(dd)%e(max(dfac(dd)%n,1)))
       dcel(dd)%n = 0;  dfac(dd)%n = 0
    end do
    do m = 1, nent
       if (e_kind(m) == 0) then
          dcel(e_rd(m))%n = dcel(e_rd(m))%n + 1;  dcel(e_rd(m))%e(dcel(e_rd(m))%n) = m
       else
          dfac(e_rd(m))%n = dfac(e_rd(m))%n + 1;  dfac(e_rd(m))%e(dfac(e_rd(m))%n) = m
       end if
    end do
    ntot = count(e_kind(1:nent) == 0)
    do t = 1, nwf
       gf(wf_list(1,t), wf_list(2,t))%a(wf_list(3,t), wf_list(4,t), wf_list(5,t)) = -1
    end do
    if (rank0v) then
       write(*,'(A,I0,A,I0,A)') '  fringe entries  : ', ntot, ' cells + ', nent - ntot, ' faces interpolated'
       write(*,'(A,I0,A,I0,A)') '  fallback (nearest-neighbour) donors: ', n_fb_used, ' fringe cells, ', n_fb_faces, ' fringe faces'
       if (nwf > 0) write(*,'(A,I0,A)') '  solid-bounded    : ', nwf, &
          ' fringe faces next to solid holes without donor -> wall faces (zero velocity)'
       if (n_bcf > 0) write(*,'(A,I0,A)') '  boundary faces   : ', n_bcf, &
          ' fringe faces on covered block faces are set by the boundary condition / exchange (no donor)'
       if (n_orphan_faces > 0) write(*,'(A,I0,A)') '  WARNING: ', n_orphan_faces, ' fringe faces without donor (orphans)'
    end if

  contains
    subroutine add_oc(b_, i_, j_, k_)
      integer, intent(in) :: b_, i_, j_, k_
      noc = noc + 1
      if (noc > size(oc_list, 2)) call grow_lists()
      oc_list(:,noc) = [b_, i_, j_, k_]
    end subroutine add_oc
    subroutine grow_lists()
      integer, allocatable :: a1(:,:), a2(:,:)
      allocate(a1(4, 2*size(oc_list,2)), a2(5, 2*size(of_list,2)));  a1 = 0;  a2 = 0
      a1(:,1:size(oc_list,2)) = oc_list;  a2(:,1:size(of_list,2)) = of_list
      call move_alloc(a1, oc_list);  call move_alloc(a2, of_list)
    end subroutine grow_lists
    subroutine push(kd, rd, ridx, dd_, bs_, ws_)
      integer,  intent(in) :: kd, rd, ridx(3), dd_, bs_(3,3)
      real(dp), intent(in) :: ws_(8,3)
      integer,  allocatable :: a1(:), a2(:), a3(:), a4(:,:), a5(:,:,:), a7(:,:)
      real(dp), allocatable :: a6(:,:,:), a8(:,:)
      if (nent == cap) then
         allocate(a1(2*cap), a2(2*cap), a3(2*cap), a4(3,2*cap), a5(3,3,2*cap), a6(8,3,2*cap))
         a1(1:cap) = t_kind;  a2(1:cap) = t_rd;  a3(1:cap) = t_dd
         a4(:,1:cap) = t_ridx;  a5(:,:,1:cap) = t_base;  a6(:,:,1:cap) = t_w
         call move_alloc(a1, t_kind);  call move_alloc(a2, t_rd);  call move_alloc(a3, t_dd)
         call move_alloc(a4, t_ridx);  call move_alloc(a5, t_base);  call move_alloc(a6, t_w)
         allocate(a7(3,2*cap), a8(6,2*cap));  a7 = 0;  a8 = 0.0_dp
         a7(:,1:cap) = t_dcel;  a8(:,1:cap) = t_box
         call move_alloc(a7, t_dcel);  call move_alloc(a8, t_box)
         cap = 2*cap
      end if
      nent = nent + 1
      t_dcel(:,nent) = 0;  t_box(:,nent) = 0.0_dp
      t_kind(nent) = kd;  t_rd(nent) = rd;  t_dd(nent) = dd_
      t_ridx(:,nent) = ridx;  t_base(:,:,nent) = bs_;  t_w(:,:,nent) = ws_
    end subroutine push
  end subroutine build_entries

  ! trilinear weights on a lattice with offsets off (0.5 = cell centres,
  ! 1.0 = faces of that family) in the donor's node-index coordinates
  pure subroutine lattice_weights(xi, off, base, w)
    real(dp), intent(in)  :: xi(3), off(3)
    integer,  intent(out) :: base(3)
    real(dp), intent(out) :: w(8)
    real(dp) :: f(3)
    integer :: t
    base = floor(xi - off)
    f = xi - off - real(base, dp)
    do t = 0, 7
       w(t+1) = merge(f(1), 1.0_dp-f(1), iand(t,1) == 1)      &
              * merge(f(2), 1.0_dp-f(2), iand(t/2,1) == 1)    &
              * merge(f(3), 1.0_dp-f(3), iand(t/4,1) == 1)
    end do
  end subroutine lattice_weights

  ! face (family q, index id) of block b usable as a donor: inside the
  ! face range (ghost rows only across covered faces) and not a hole
  subroutine grow_wf()
    integer, allocatable :: tmp(:,:)
    allocate(tmp(5, 2*size(wf_list,2)));  tmp(:,1:size(wf_list,2)) = wf_list
    call move_alloc(tmp, wf_list)
  end subroutine grow_wf

  ! weights of the three staggered donor lattices of block ob at the
  ! natural position xi (index + local coordinate); ok = every component
  ! has at least one valid face
  subroutine face_weights(ob, xi, bs, ws, ok)
    integer, intent(in)   :: ob
    real(dp), intent(in)  :: xi(3)
    integer, intent(out)  :: bs(3,3)
    real(dp), intent(out) :: ws(8,3)
    logical, intent(out)  :: ok
    integer :: l, t
    real(dp) :: off(3)
    ok = .true.
    do l = 1, 3
       off = 0.5_dp;  off(l) = 1.0_dp
       call lattice_weights(xi, off, bs(:,l), ws(:,l))
       do t = 0, 7
          if (ws(t+1,l) <= 0.0_dp) cycle
          if (.not. face_valid(ob, l, bs(:,l) + [iand(t,1), iand(t/2,1), iand(t/4,1)])) ws(t+1,l) = 0.0_dp
       end do
       if (sum(ws(:,l)) <= 1.0e-12_dp) then
          ok = .false.
       else
          ws(:,l) = ws(:,l) / sum(ws(:,l))
       end if
    end do
  end subroutine face_weights

  logical function face_valid(b, q, id)
    integer, intent(in) :: b, q, id(3)
    integer :: n(3), a, sd, t1, t2, cid(3), ng_dir
    n = [blk(b)%nx, blk(b)%ny, blk(b)%nz]
    face_valid = .false.
    do a = 1, 3
       if (a == q) then
          if (id(a) < 0 .or. id(a) > n(a)) return
       else
          if (id(a) < 0 .or. id(a) > n(a)+1) return
       end if
    end do
    if (gf(q,b)%a(id(1),id(2),id(3)) == -1) return
    ! ghost rows (cell-type directions): at most ONE ghost direction
    ! (corner ghosts are never used), and only across a covered face
    ng_dir = 0
    do a = 1, 3
       if (a == q) cycle
       if (id(a) == 0 .or. id(a) == n(a)+1) ng_dir = ng_dir + 1
    end do
    if (ng_dir > 1) return
    do a = 1, 3
       if (a == q) cycle
       if (id(a) == 0 .or. id(a) == n(a)+1) then
          sd = merge(1, 2, id(a) == 0)
          call tang(a, t1, t2)
          cid = id;  cid(a) = merge(1, n(a), sd == 1)
          cid(q) = min(max(id(q), 1), n(q))
          if (.not. cov(sd,a,b)%a(cid(t1), cid(t2))) return
       end if
    end do
    face_valid = .true.
  end function face_valid

  pure function face_center(b, q, i, j, k) result(c)
    integer, intent(in) :: b, q, i, j, k
    real(dp) :: c(3)
    integer :: n0(3), e(3), g(3), t1, t2, n1(3), n2(3), n3(3)
    call tang(q, t1, t2)
    e = 0;  e(t1) = 1;  g = 0;  g(t2) = 1
    n0 = [i, j, k];  n0(q) = n0(q) + 1              ! face f -> node f+1
    n1 = n0 + e;  n2 = n0 + g;  n3 = n0 + e + g
    c(1) = 0.25_dp*(blk(b)%xn(n0(1),n0(2),n0(3)) + blk(b)%xn(n1(1),n1(2),n1(3)) + blk(b)%xn(n2(1),n2(2),n2(3)) + blk(b)%xn(n3(1),n3(2),n3(3)))
    c(2) = 0.25_dp*(blk(b)%yn(n0(1),n0(2),n0(3)) + blk(b)%yn(n1(1),n1(2),n1(3)) + blk(b)%yn(n2(1),n2(2),n2(3)) + blk(b)%yn(n3(1),n3(2),n3(3)))
    c(3) = 0.25_dp*(blk(b)%zn(n0(1),n0(2),n0(3)) + blk(b)%zn(n1(1),n1(2),n1(3)) + blk(b)%zn(n2(1),n2(2),n2(3)) + blk(b)%zn(n3(1),n3(2),n3(3)))
  end function face_center

  ! sub-domain of block b owning cell id (clipped into the block)
  pure integer function dom_of_cell(b, id)
    integer, intent(in) :: b, id(3)
    integer :: d, c(3)
    c = min(max(id, 1), [blk(b)%nx, blk(b)%ny, blk(b)%nz])
    dom_of_cell = blk(b)%dom0
    do d = blk(b)%dom0, blk(b)%dom0 + blk(b)%ncpu - 1
       if (all(c >= dom(d)%r1) .and. all(c <= dom(d)%r2)) then
          dom_of_cell = d;  return
       end if
    end do
  end function dom_of_cell

  ! sub-domain responsible for face (q, id): the one owning cell id with
  ! the face index mapped to a cell (face 0 -> cell 1)
  pure integer function dom_of_face(b, q, id)
    integer, intent(in) :: b, q, id(3)
    integer :: c(3)
    c = id;  c(q) = max(id(q), 1)
    dom_of_face = dom_of_cell(b, c)
  end function dom_of_face

  !=====================================================================
  ! Per-domain status slices (ghost-extended).
  !=====================================================================
  subroutine domain_status
    integer :: ib, d, b, i1, i2, j1, j2, k1, k2
    do ib = 1, nown
       d = own_list(ib);  b = dom(d)%b
       i1 = dom(d)%r1(1)-ng;  i2 = dom(d)%r2(1)+ng
       j1 = dom(d)%r1(2)-ng;  j2 = dom(d)%r2(2)+ng
       k1 = dom(d)%r1(3)-ng;  k2 = dom(d)%r2(3)+ng
       allocate(dom(d)%cs(i1:i2,j1:j2,k1:k2))
       dom(d)%cs = gcs(b)%a(i1:i2,j1:j2,k1:k2)
       allocate(dom(d)%fs1(i1:i2-1,j1:j2,k1:k2), dom(d)%fs2(i1:i2,j1:j2-1,k1:k2), dom(d)%fs3(i1:i2,j1:j2,k1:k2-1))
       dom(d)%fs1 = gf(1,b)%a(i1:i2-1,j1:j2,k1:k2)
       dom(d)%fs2 = gf(2,b)%a(i1:i2,j1:j2-1,k1:k2)
       dom(d)%fs3 = gf(3,b)%a(i1:i2,j1:j2,k1:k2-1)
       allocate(dom(d)%ovd(i1:i2,j1:j2,k1:k2));  dom(d)%ovd = 0.0_dp
    end do
  end subroutine domain_status

  !=====================================================================
  ! Zero the hole cells / faces of every owned domain (after the initial
  ! condition), so blanked positions never carry data.
  !=====================================================================
  subroutine ovs_finish_status
    integer :: ib, d, i, j, k
    if (.not. ovst_on) return
    do ib = 1, nown
       d = own_list(ib)
       call set_block(d)
       do k = lbound(p,3), ubound(p,3)
        do j = lbound(p,2), ubound(p,2)
         do i = lbound(p,1), ubound(p,1)
            if (ocs(i,j,k) < 0) then
               p(i,j,k) = 0.0_dp;  alpha(i,j,k) = 0.0_dp
            end if
            if (i <= ubound(ofs1,1)) then
               if (ofs1(i,j,k) == -1) then
                  u(i,j,k) = 0.0_dp;  us(i,j,k) = 0.0_dp;  U1(i,j,k) = 0.0_dp;  Us1(i,j,k) = 0.0_dp
               end if
            end if
            if (j <= ubound(ofs2,2)) then
               if (ofs2(i,j,k) == -1) then
                  v(i,j,k) = 0.0_dp;  vs(i,j,k) = 0.0_dp;  U2(i,j,k) = 0.0_dp;  Us2(i,j,k) = 0.0_dp
               end if
            end if
            if (k <= ubound(ofs3,3)) then
               if (ofs3(i,j,k) == -1) then
                  w(i,j,k) = 0.0_dp;  ws(i,j,k) = 0.0_dp;  U3(i,j,k) = 0.0_dp;  Us3(i,j,k) = 0.0_dp
               end if
            end if
         end do
        end do
       end do
    end do
  end subroutine ovs_finish_status

  !=====================================================================
  ! RUNTIME EXCHANGE
  !   ovs_update_cells(iv): interpolate the cell field iv of the donors
  !     into ovs_val(1,:) of every cell entry; alpha (8) is written into
  !     the fringe cells, p (4) and the Krylov vector kv1 (9) are kept in
  !     ovs_val for the implicit fringe rows of the pressure solver.
  !   ovs_update_faces(ivbase): Cartesian velocity vector at every fringe
  !     face from the donor's staggered lattices (ivbase 1: u,v,w;
  !     5: us,vs,ws) -> Cartesian component + transformed component.
  !=====================================================================
  subroutine ovs_update_cells(iv)
    integer, intent(in) :: iv
    integer :: ie_, m, e
    real(dp), pointer :: f(:,:,:)
    if (.not. ovst_on .or. nent == 0) return
    call exchange(0, iv)
    ! written into the fringe cells: alpha (8), the vapour fraction (19,
    ! update 15), the temperature (20, update 16) and the turbulence fields (14-16, update 15: before,
    ! the fringe cells of k / omega / mu_t kept their initial values)
    if (iv == 8 .or. iv == 19 .or. iv == 20 .or. (iv >= 14 .and. iv <= 16)) then
       do m = 1, ndom
          if (dom(m)%owner /= myrank) cycle
          f => field_ptr(m, iv)
          do ie_ = 1, dcel(m)%n
             e = dcel(m)%e(ie_)
             f(e_ridx(1,e), e_ridx(2,e), e_ridx(3,e)) = ovs_val(1, e)
          end do
       end do
    end if
  end subroutine ovs_update_cells

  !---------------------------------------------------------------------
  ! Two-phase: the pressure is interpolated to a fringe cell as
  !   p_f = I(p - p_hyd) + p_hyd(z_f),   p_hyd(z) = rho_ref(z) |g| (z_fs - z)
  ! (rho_ref = 1 below, rho2/rho1 above the still free surface), i.e.
  ! the still-water hydrostatic part is subtracted before the donor
  ! interpolation and added back at the receiver's own height.  With
  ! the raw interpolation a fringe cell in air next to the free surface
  ! takes a water-like pressure from a donor cell that straddles the
  ! interface (or from a nearest-neighbour fallback donor), and the
  ! neighbouring active air cells then see a pressure gradient of the
  ! order rho_w g dz / rho_a -> spurious velocities of tens of m/s in
  ! the first step (Box_6DoF).  A still hydrostatic field is now
  ! reproduced exactly by the fringe rows whatever the donors are.
  ! ovs_ph(e) = p_hyd(z_f) - I p_hyd is built once per assembly by
  ! interpolating the cell field p_hyd(z) (carried in kv1, id 9).
  !---------------------------------------------------------------------
  subroutine ovs_build_phyd
    integer :: m, ie_, e, i, j, k
    real(dp), pointer :: f(:,:,:)
    if (allocated(ovs_ph)) deallocate(ovs_ph)
    allocate(ovs_ph(max(nent,1)));  ovs_ph = 0.0_dp
    if (.not. ovst_on .or. nent == 0) return
    if (Num_of_Phase /= 2 .or. gnd(3) == 0.0_dp) return
    do m = 1, ndom
       if (dom(m)%owner /= myrank) cycle
       if (.not. allocated(dom(m)%kv1)) return          ! grid-check mode: no fields
    end do
    do m = 1, ndom
       if (dom(m)%owner /= myrank) cycle
       call set_block(m)
       f => field_ptr(m, 9)
       do k = lbound(f,3), ubound(f,3)
        do j = lbound(f,2), ubound(f,2)
         do i = lbound(f,1), ubound(f,1)
            f(i,j,k) = phyd(Xh(3, min(max(2*i,hlo(1)+1),hhi(1)-1), min(max(2*j,hlo(2)+1),hhi(2)-1), &
                                    min(max(2*k,hlo(3)+1),hhi(3)-1)))
         end do
        end do
       end do
    end do
    call ovs_update_cells(9)                 ! ovs_val(1,e) = I p_hyd
    do m = 1, ndom
       if (dom(m)%owner /= myrank) cycle
       call set_block(m)
       do ie_ = 1, dcel(m)%n
          e = dcel(m)%e(ie_)
          i = e_ridx(1,e);  j = e_ridx(2,e);  k = e_ridx(3,e)
          ovs_ph(e) = phyd(Xh(3,2*i,2*j,2*k)) - ovs_val(1,e)
       end do
       dom(m)%kv1 = 0.0_dp
    end do
  contains
    pure real(dp) function phyd(z)
      real(dp), intent(in) :: z
      phyd = -gnd(3)*(zfs - z)*merge(1.0_dp, rho_ratio, z <= zfs)
    end function phyd
  end subroutine ovs_build_phyd

  subroutine ovs_update_faces(ivbase)
    integer, intent(in) :: ivbase
    integer :: ie_, m, e, q, H(3), id(3)
    real(dp) :: uc(3), val
    real(dp), pointer :: Uqp(:,:,:), cq(:,:,:)
    if (.not. ovst_on .or. nent == 0) return
    call exchange(1, ivbase)
    do m = 1, ndom
       if (dom(m)%owner /= myrank) cycle
       call set_block(m)
       do ie_ = 1, dfac(m)%n
          e = dfac(m)%e(ie_)
          q = e_kind(e);  id = e_ridx(:,e)
          uc = ovs_val(:,e)
          H = 2*id;  H(q) = H(q) + 1
          val = fvw(H(1),H(2),H(3))*( mq(1,q,H(1),H(2),H(3))*uc(1)              &
              + mq(2,q,H(1),H(2),H(3))*uc(2) + mq(3,q,H(1),H(2),H(3))*uc(3) )
          Uqp => field_ptr(m, ivbase - 1 + q)
          Uqp(id(1),id(2),id(3)) = val
          cq => cart_ptr(m, ivbase, q)
          cq(id(1),id(2),id(3)) = uc(q)
       end do
    end do
  end subroutine ovs_update_faces

  !---------------------------------------------------------------------
  ! Pack (donor side) / exchange / unpack into ovs_val.
  !---------------------------------------------------------------------
  subroutine exchange(kd, iv)
    integer, intent(in) :: kd, iv
    type :: buf_t
       real(dp), allocatable :: a(:)
    end type buf_t
    type(buf_t), allocatable :: sb(:), rb(:)
    integer, allocatable :: reqs(:)
    integer :: ip, nv, nsel, sz, nreq, ierr, e, n, l
    real(dp), pointer :: f(:,:,:)

    nv = merge(1, 3, kd == 0)
    allocate(sb(max(npair,1)), rb(max(npair,1)), reqs(2*max(npair,1)))
    nreq = 0
    do ip = 1, npair
       nsel = merge(p_nc(ip), p_nfc(ip), kd == 0)
       if (nsel == 0) cycle
       sz = nsel*nv
       if (p_s(ip) == myrank) then
          allocate(sb(ip)%a(sz))
          n = 0
          do e = p_e1(ip), p_e2(ip)
             if ((e_kind(e) == 0) .neqv. (kd == 0)) cycle
             if (kd == 0) then
                f => field_ptr(e_dd(e), iv)
                n = n + 1;  sb(ip)%a(n) = interp(f, dom(e_dd(e))%r1 - ng, e_base(:,1,e), e_w(:,1,e))
                if (iv == 8 .and. vof_scheme == 3) sb(ip)%a(n) = alpha_geometric(f, e_dd(e), e, sb(ip)%a(n))
             else
                do l = 1, 3
                   f => cart_ptr(e_dd(e), iv, l)
                   n = n + 1;  sb(ip)%a(n) = interp(f, dom(e_dd(e))%r1 - ng, e_base(:,l,e), e_w(:,l,e))
                end do
             end if
          end do
          if (p_r(ip) /= myrank) then
             nreq = nreq + 1
             call MPI_Isend(sb(ip)%a, sz, MPI_DOUBLE_PRECISION, p_r(ip), 100000+ip, MPI_COMM_WORLD, reqs(nreq), ierr)
          end if
       end if
       if (p_r(ip) == myrank .and. p_s(ip) /= myrank) then
          allocate(rb(ip)%a(sz))
          nreq = nreq + 1
          call MPI_Irecv(rb(ip)%a, sz, MPI_DOUBLE_PRECISION, p_s(ip), 100000+ip, MPI_COMM_WORLD, reqs(nreq), ierr)
       end if
    end do
    if (nreq > 0) call MPI_Waitall(nreq, reqs, MPI_STATUSES_IGNORE, ierr)
    do ip = 1, npair
       if (p_r(ip) /= myrank) cycle
       nsel = merge(p_nc(ip), p_nfc(ip), kd == 0)
       if (nsel == 0) cycle
       n = 0
       do e = p_e1(ip), p_e2(ip)
          if ((e_kind(e) == 0) .neqv. (kd == 0)) cycle
          do l = 1, nv
             n = n + 1
             if (p_s(ip) == myrank) then
                ovs_val(l,e) = sb(ip)%a(n)
             else
                ovs_val(l,e) = rb(ip)%a(n)
             end if
          end do
       end do
    end do
    deallocate(sb, rb, reqs)
  end subroutine exchange

  !---------------------------------------------------------------------
  ! GEOMETRIC fringe value of alpha (vof_scheme = 3): the liquid fraction
  ! of the receiver cell's box (e_box, natural coordinates of the donor
  ! cell e_dcel that contains the receiver centre) on the liquid side of
  ! the donor cell's PLIC plane.  Used when the donor cell is mixed and
  ! the box spans at most two donor cells per axis; otherwise the
  ! trilinear value atri is kept.
  !---------------------------------------------------------------------
  real(dp) function alpha_geometric(f, d, e, atri) result(av)
    real(dp), intent(in) :: f(:,:,:)
    integer,  intent(in) :: d, e
    real(dp), intent(in) :: atri
    integer  :: o(3), c(3), cc(3), d1, d2, d3, ax, k1, k2, k3, klo(3), khi(3)
    real(dp) :: st(-1:1,-1:1,-1:1), m(3), acst, x0(3), w(3), ad, blo(3), bhi(3), slo(3), shi(3), vtot, vliq, vsub, fr
    logical  :: mir(3)
    av = atri
    if (e_dcel(1,e) == 0) return
    c = e_dcel(:,e)
    o = dom(d)%r1 - ng - 1                        ! f(i - o) holds global cell i
    blo = e_box(1:3,e);  bhi = e_box(4:6,e);  w = bhi - blo
    if (any(w <= 0.0_dp) .or. any(w > 2.0_dp)) return
    ! donor cells covered by the receiver box (relative index offsets)
    do ax = 1, 3
       klo(ax) = floor(blo(ax) + 1.0e-9_dp);  khi(ax) = ceiling(bhi(ax) - 1.0e-9_dp) - 1
       if (khi(ax) < klo(ax)) khi(ax) = klo(ax)
    end do
    ! all covered donor cells (plus the reconstruction stencil) must be stored on this rank
    if (any(c + klo - 2 - o < 1) .or. c(1)+khi(1)+2-o(1) > size(f,1) .or. c(2)+khi(2)+2-o(2) > size(f,2) &
        .or. c(3)+khi(3)+2-o(3) > size(f,3)) return
    ! the receiver box is split at the donor cell boundaries: every donor
    ! cell contributes its own liquid fraction of the overlap - the PLIC
    ! plane of a mixed cell restricted to that cell, 0 / 1 for pure cells.
    ! (No plane is extended beyond its own cell, so an interface cell can
    !  never put liquid into the neighbouring air cells or vice versa.)
    vtot = 0.0_dp;  vliq = 0.0_dp
    do k3 = klo(3), khi(3)
     do k2 = klo(2), khi(2)
      do k1 = klo(1), khi(1)
         cc = c + [k1, k2, k3]
         slo = max(blo - real([k1,k2,k3],dp), 0.0_dp);  shi = min(bhi - real([k1,k2,k3],dp), 1.0_dp)
         if (any(shi <= slo)) cycle
         vsub = product(shi - slo)
         ad = f(cc(1)-o(1), cc(2)-o(2), cc(3)-o(3))
         if (ad <= 0.0_dp) then
            fr = 0.0_dp
         else if (ad >= 1.0_dp) then
            fr = 1.0_dp
         else
            do d3 = -1, 1
               do d2 = -1, 1
                  do d1 = -1, 1
                     st(d1,d2,d3) = f(cc(1)+d1-o(1), cc(2)+d2-o(2), cc(3)+d3-o(3))
                  end do
               end do
            end do
            call plic_normal27(st, m, mir)
            acst = plic_plane_const(m, ad)
            x0 = slo
            do ax = 1, 3
               if (mir(ax)) x0(ax) = 1.0_dp - shi(ax)
            end do
            fr = plic_box_frac(m, acst, x0, shi - slo)
         end if
         vtot = vtot + vsub;  vliq = vliq + fr*vsub
      end do
     end do
    end do
    if (vtot <= 0.0_dp) return
    av = min(max(vliq/vtot, 0.0_dp), 1.0_dp)
  end function alpha_geometric

  function cart_ptr(d, ivbase, l) result(f)
    integer, intent(in) :: d, ivbase, l
    real(dp), pointer :: f(:,:,:)
    if (ivbase == 1) then
       select case (l)
       case (1);  f => dom(d)%u
       case (2);  f => dom(d)%v
       case default;  f => dom(d)%w
       end select
    else
       select case (l)
       case (1);  f => dom(d)%us
       case (2);  f => dom(d)%vs
       case default;  f => dom(d)%ws
       end select
    end if
  end function cart_ptr

  ! trilinear sum over the donor stencil; lo = storage lower bounds of
  ! the donor domain (assumed-shape dummies lose them)
  real(dp) function interp(f, lo, base, w)
    real(dp), intent(in) :: f(:,:,:)
    integer,  intent(in) :: lo(3), base(3)
    real(dp), intent(in) :: w(8)
    integer :: t, id(3)
    interp = 0.0_dp
    do t = 0, 7
       if (w(t+1) == 0.0_dp) cycle
       id = base + [iand(t,1), iand(t/2,1), iand(t/4,1)] - lo + 1
       interp = interp + w(t+1)*f(id(1),id(2),id(3))
    end do
  end function interp

  !=====================================================================
  ! Node Blank of block b (rank 0, output; update 16b): the EXACT cell
  ! status without any combination of neighbouring cells - node (i,j,k)
  ! carries the status of cell (i,j,k) (the cell whose low corner is the
  ! node; the nodes on the max faces i = ni / j = nj / k = nk repeat the
  ! last cell): 1 active (solved), 0 fringe (interpolated), -1 hole
  ! (dominated or solid).  Replaces the max / solid rule of update 12c
  ! and the separate <case>_NNNNNNNN_cblank.p3d file.
  !=====================================================================
  subroutine ovs_node_iblank(b, ibl)
    integer, intent(in) :: b
    integer, intent(out) :: ibl(blk(b)%ni, blk(b)%nj, blk(b)%nk)
    integer :: i, j, k
    do k = 1, blk(b)%nk
     do j = 1, blk(b)%nj
      do i = 1, blk(b)%ni
         ibl(i,j,k) = max(-1, cs(b)%a(min(i, blk(b)%nx), min(j, blk(b)%ny), min(k, blk(b)%nz)))
      end do
     end do
    end do
  end subroutine ovs_node_iblank

  !=====================================================================
  ! Self-test of the interpolation exchange (poi_verbose): the cell
  ! field kv1 and the face fields u,v,w are set to analytic linear
  ! functions of the coordinates, exchanged, and checksums / errors are
  ! printed - they must agree between different rank counts (the
  ! overset tables are replicated; only the message path differs).
  !=====================================================================
  subroutine ovs_selftest
    integer :: ib, d, i, j, k, e, m, ierr, H(3)
    real(dp) :: sc, sf, ec, ef, x(3), xr(3), tmp(4)
    if (.not. ovst_on) return
    do ib = 1, nown
       d = own_list(ib);  call set_block(d)
       do k = lbound(kv1,3), ubound(kv1,3)
        do j = lbound(kv1,2), ubound(kv1,2)
         do i = lbound(kv1,1), ubound(kv1,1)
            H = [2*i, 2*j, 2*k]
            if (any(H < hlo) .or. any(H > hhi)) cycle
            kv1(i,j,k) = lin(Xh(:,2*i,2*j,2*k))
            if (2*i+1 <= hhi(1)) u(i,j,k) = lin(Xh(:,2*i+1,2*j,2*k))
            if (2*j+1 <= hhi(2)) v(i,j,k) = 2.0_dp*lin(Xh(:,2*i,2*j+1,2*k))
            if (2*k+1 <= hhi(3)) w(i,j,k) = 3.0_dp*lin(Xh(:,2*i,2*j,2*k+1))
         end do
        end do
       end do
    end do
    call exchange(0, 9)
    sc = 0.0_dp;  ec = 0.0_dp
    do m = 1, ndom
       if (dom(m)%owner /= myrank) cycle
       call set_block(m)
       do e = 1, dcel(m)%n
          i = e_ridx(1,dcel(m)%e(e));  j = e_ridx(2,dcel(m)%e(e));  k = e_ridx(3,dcel(m)%e(e))
          sc = sc + ovs_val(1,dcel(m)%e(e))
          ec = max(ec, abs(ovs_val(1,dcel(m)%e(e)) - lin(Xh(:,2*i,2*j,2*k))))
       end do
    end do
    call exchange(1, 1)
    sf = 0.0_dp;  ef = 0.0_dp
    do m = 1, ndom
       if (dom(m)%owner /= myrank) cycle
       call set_block(m)
       do e = 1, dfac(m)%n
          i = e_ridx(1,dfac(m)%e(e));  j = e_ridx(2,dfac(m)%e(e));  k = e_ridx(3,dfac(m)%e(e))
          H = [2*i, 2*j, 2*k];  H(e_kind(dfac(m)%e(e))) = H(e_kind(dfac(m)%e(e))) + 1
          xr = Xh(:,H(1),H(2),H(3))
          sf = sf + sum(ovs_val(:,dfac(m)%e(e)))
          ef = max(ef, abs(ovs_val(1,dfac(m)%e(e)) - lin(xr)), abs(ovs_val(2,dfac(m)%e(e)) - 2.0_dp*lin(xr)), &
                       abs(ovs_val(3,dfac(m)%e(e)) - 3.0_dp*lin(xr)))
       end do
    end do
    tmp = [sc, sf, ec, ef]
    call MPI_Allreduce(MPI_IN_PLACE, tmp(1:2), 2, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_Allreduce(MPI_IN_PLACE, tmp(3:4), 2, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    if (myrank == 0) then
       write(*,'(A,2ES22.14)') '  overset self-test: checksum cells / faces ', tmp(1), tmp(2)
       write(*,'(A,2ES12.4)')  '                     max interpolation error of a linear field, cells / faces ', tmp(3), tmp(4)
    end if
    ! restore the initial condition of the test fields
    do ib = 1, nown
       d = own_list(ib);  call set_block(d)
       kv1 = 0.0_dp
    end do
  contains
    pure real(dp) function lin(xx)
      real(dp), intent(in) :: xx(3)
      lin = 0.7_dp + xx(1) + 2.0_dp*xx(2) + 3.0_dp*xx(3)
    end function lin
  end subroutine ovs_selftest

  !=====================================================================
  ! Debug: print the entry (donor stencil) of cell (i,j,k) of block b and
  ! the status of the donor cells (rank 0 only, replicated data).
  !=====================================================================
  subroutine ovs_dump_cell(b, i, j, k)
    integer, intent(in) :: b, i, j, k
    integer :: e, t, id(3), ob
    write(*,'(A,4I5,A,I3)') '    dump cell blk,i,j,k', b, i, j, k, ' status ', cs(b)%a(i,j,k)
    do t = max(1,k-2), min(blk(b)%nz,k+2)
       write(*,'(A,I4,A,27I3)') '      status box i-2..i+2 x j-1..j+1 at k=', t, ':', &
          cs(b)%a(max(1,i-2):min(blk(b)%nx,i+2), max(1,j-1):min(blk(b)%ny,j+1), t)
    end do
    do e = 1, nent
       if (e_kind(e) /= 0) cycle
       if (dom(e_rd(e))%b /= b .or. any(e_ridx(:,e) /= [i,j,k])) cycle
       ob = dom(e_dd(e))%b
       write(*,'(A,I3,A,3I5)') '      donor block ', ob, ' base ', e_base(:,1,e)
       do t = 0, 7
          id = e_base(:,1,e) + [iand(t,1), iand(t/2,1), iand(t/4,1)]
          if (any(id < 1) .or. id(1) > blk(ob)%nx .or. id(2) > blk(ob)%ny .or. id(3) > blk(ob)%nz) then
             write(*,'(A,3I5,A,F8.4,A)') '        cell ', id, ' w ', e_w(t+1,1,e), '  (ghost)'
          else
             write(*,'(A,3I5,A,F8.4,A,I3)') '        cell ', id, ' w ', e_w(t+1,1,e), '  status ', cs(ob)%a(id(1),id(2),id(3))
          end if
       end do
    end do
  end subroutine ovs_dump_cell

  !=====================================================================
  ! Report file <case>.ovst (rank 0).
  !=====================================================================
  ! ---------------------------------------------------------------------
  ! Verification of the user's condition "no grid point inside a solid
  ! body": every ACTIVE cell within the cut layer of a closed body (its
  ! centre within 0.6 cell diagonals of the surface) is tested with the
  ! winding number at its centre and at its 8 nodes, 12 edge midpoints
  ! and 6 face centres.  Counts are reported in the log and in the .ovst
  ! file and the first offending cells are listed.  Grid-check mode
  ! (Ovst_only = 1) only: the winding-number queries cost as much as the
  ! hole cut itself (ONRT: 25 min).
  ! ---------------------------------------------------------------------
  subroutine check_permeation
    integer :: b, i, j, k, ib, t, q, iu, iv, iw, tb, region, nlist
    real(dp) :: p(3), q_(3), dist, wbc(3), diag, bmin(3), bmax(3), marg, xn(3,8), u, v, w, wt, wn
    logical :: cen_in, any_in
    n_perm_cen = 0;  n_perm_pt = 0;  nlist = 0
    if (.not. hsurf%wind_ready) return
    do ib = 1, nbody
       if (.not. body_closed(ib)) cycle
       bmin = huge(1.0_dp);  bmax = -huge(1.0_dp)
       do t = 1, wsurf%nt
          if (wsurf%tb(t) /= ib) cycle
          do q = 1, 3
             bmin = min(bmin, wsurf%v(:,wsurf%t(q,t)));  bmax = max(bmax, wsurf%v(:,wsurf%t(q,t)))
          end do
       end do
       marg = 1.0e-6_dp*maxval(bmax - bmin)
       do b = 1, nblk
          do k = 1, blk(b)%nz
           do j = 1, blk(b)%ny
            do i = 1, blk(b)%nx
               if (cs(b)%a(i,j,k) /= 1) cycle
               if (maxval(blk(b)%xn(i:i+1,j:j+1,k:k+1)) < bmin(1)-marg .or. minval(blk(b)%xn(i:i+1,j:j+1,k:k+1)) > bmax(1)+marg) cycle
               if (maxval(blk(b)%yn(i:i+1,j:j+1,k:k+1)) < bmin(2)-marg .or. minval(blk(b)%yn(i:i+1,j:j+1,k:k+1)) > bmax(2)+marg) cycle
               if (maxval(blk(b)%zn(i:i+1,j:j+1,k:k+1)) < bmin(3)-marg .or. minval(blk(b)%zn(i:i+1,j:j+1,k:k+1)) > bmax(3)+marg) cycle
               p = cell_center(b, i, j, k)
               diag = sqrt((maxval(blk(b)%xn(i:i+1,j:j+1,k:k+1)) - minval(blk(b)%xn(i:i+1,j:j+1,k:k+1)))**2 + &
                           (maxval(blk(b)%yn(i:i+1,j:j+1,k:k+1)) - minval(blk(b)%yn(i:i+1,j:j+1,k:k+1)))**2 + &
                           (maxval(blk(b)%zn(i:i+1,j:j+1,k:k+1)) - minval(blk(b)%zn(i:i+1,j:j+1,k:k+1)))**2)
               call surf_closest(hsurf, p, 2, ib, 0.6_dp*diag, dist, q_, tb, wbc, region)
               if (tb <= 0) cycle
               if (dist > 0.6_dp*diag) cycle
               ! the body's own wall cells (nodes ON the surface) are outside by construction
               if (near_own_wall(b, i, j, k)) then
                  if (own_surface_at(b, q_)) cycle
               end if
               do t = 0, 7
                  xn(:,t+1) = [blk(b)%xn(i+iand(t,1), j+iand(t/2,1), k+iand(t/4,1)),   &
                               blk(b)%yn(i+iand(t,1), j+iand(t/2,1), k+iand(t/4,1)),   &
                               blk(b)%zn(i+iand(t,1), j+iand(t/2,1), k+iand(t/4,1))]
               end do
               cen_in = (surf_wind(hsurf, p, ib) > 0.5_dp + 1.0e-2_dp)
               any_in = cen_in
               if (.not. any_in) then
                  do t = 0, 26
                     iu = mod(t,3);  iv = mod(t/3,3);  iw = t/9
                     if (iu == 1 .and. iv == 1 .and. iw == 1) cycle
                     u = 0.5_dp*iu;  v = 0.5_dp*iv;  w = 0.5_dp*iw
                     p = 0.0_dp
                     do q = 0, 7
                        wt = merge(u, 1.0_dp-u, iand(q,1) == 1)*merge(v, 1.0_dp-v, iand(q/2,1) == 1)* &
                             merge(w, 1.0_dp-w, iand(q/4,1) == 1)
                        p = p + wt*xn(:,q+1)
                     end do
                     wn = surf_wind(hsurf, p, ib)
                     if (wn > 0.5_dp + 1.0e-2_dp) then
                        any_in = .true.;  exit
                     end if
                  end do
               end if
               if (cen_in) n_perm_cen = n_perm_cen + 1
               if (any_in) then
                  n_perm_pt = n_perm_pt + 1
                  if (rank0v .and. ovst_only /= 0 .and. nlist < 40 .and. (cen_in .or. nlist < 10)) then
                     nlist = nlist + 1
                     write(*,'(A,I0,A,3I5,A,I0,A,3ES12.4,A,L1,A,I0,A,L1)') '    PERMEATION: block ', b, ' cell ', i, j, k, &
                        ' body ', ib, ' centre ', cell_center(b, i, j, k), ' centre inside ', cen_in, &
                        ' typ ', typ(b)%a(i,j,k), ' near own wall ', near_own_wall(b, i, j, k)
                  end if
               end if
            end do
           end do
          end do
       end do
    end do
    if (rank0v) write(*,'(A,I0,A,I0,A)') '  permeation check : ', n_perm_pt, &
       ' active cells with a point (centre / node / edge / face) inside a body, ', n_perm_cen, ' with the centre inside'
  end subroutine check_permeation

  subroutine write_report(tsec)
    real(dp), intent(in) :: tsec
    integer :: iu, b, ax, sd, ib, ip, nun
    character(len=5) :: fn(-3:3)
    fn(-3) = 'k-max';  fn(-2) = 'j-max';  fn(-1) = 'i-max';  fn(0) = '?'
    fn(1) = 'i-min';  fn(2) = 'j-min';  fn(3) = 'k-min'
    open(newunit=iu, file=trim(casename)//'.ovst', status='replace', action='write')
    write(iu,'(A)') '! probCFD overset grid assembly report, case '//trim(casename)
    write(iu,'(A,I0)')     '! blocks              : ', nblk
    write(iu,'(A,I0)')     '! fringe layers       : ', ovst_nf
    write(iu,'(A,ES12.4)') '! Surfovl_tol [m]     : ', surfovl_tol
    write(iu,'(A,ES12.4)') '! Ovst_vol_tol        : ', ovst_vtol
    write(iu,'(A,F10.2,A)')'! assembly wall time  : ', tsec, ' s'
    write(iu,'(A)') '!'
    write(iu,'(A)') '! ---- overset boundary faces (not in .bcs / .bcs_b2b) ----'
    do b = 1, nblk
       do ax = 1, 3
          do sd = 1, 2
             if (all(cov(sd,ax,b)%a)) cycle
             nun = count(.not. cov(sd,ax,b)%a)
             write(iu,'(A,I0,A,A,A,I0,A,I0,A)') '  block ', b, ' face ', fn(merge(ax,-ax,sd==1)), &
                ' : ', nun, ' of ', size(cov(sd,ax,b)%a), ' face cells are overset boundary'
          end do
       end do
    end do
    write(iu,'(A)') '!'
    write(iu,'(A)') '! ---- no-slip wall patches and solid bodies ----'
    do ip = 1, nwp
       write(iu,'(A,I0,A,I0,A,A,A,I0,A,I0,A,I0)') '  wall ', ip, ': block ', wp_blk(ip), ' face ', &
          fn(blk(wp_blk(ip))%bc_dir(wp_bc(ip))), ' (', wp_n1(ip), ' x ', wp_n2(ip), ' nodes), body ', wp_body(ip)
    end do
    do ib = 1, nbody
       write(iu,'(A,I0,A,I0,A)') '  body ', ib, ': ', body_np(ib), ' wall patches'
       write(iu,'(A,I0)')  '     edge nodes connected by surface overlap (Surfovl_tol) : ', body_ovl(ib)
       write(iu,'(A,I0)')  '     edge nodes on block-to-block interfaces              : ', body_b2b(ib)
       write(iu,'(A,I0)')  '     edge nodes on physical domain boundaries (Surfovl_tol): ', body_dom(ib)
       write(iu,'(A,I0)')  '     unconnected edge nodes                               : ', body_unc(ib)
       if (body_closed(ib)) then
          write(iu,'(A,I0,A)') '     -> CLOSED body: ', body_nsol(ib), ' cells inside are blanked (holes)'
       else
          write(iu,'(A)')      '     -> OPEN surface: no hole cutting for this body'
       end if
    end do
    write(iu,'(A)') '!'
    write(iu,'(A)') '! ---- cell classification (1 active, 0 fringe, -1 hole) ----'
    do b = 1, nblk
       write(iu,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0)') '  block ', b, ': active ', count(cs(b)%a == 1), &
          '  fringe ', count(cs(b)%a == 0), '  hole ', count(cs(b)%a < 0), ' (inside solid ', count(cs(b)%a == -2), &
          ')  total ', size(cs(b)%a), '  orphan active cells ', orph_blk(b)
    end do
    write(iu,'(A,I0,A,I0)') '  interpolated fringe cells / faces : ', count(e_kind(1:nent) == 0), ' / ', &
         nent - count(e_kind(1:nent) == 0)
    write(iu,'(A,I0,A,I0)') '  orphans (no donor) cells / faces  : ', n_orphan_cells, ' / ', n_orphan_faces
    write(iu,'(A,I0,A,I0)') '  orphans closed by a same-block neighbour (zero gradient) cells / faces : ', n_fill_cells, ' / ', n_fill_faces
    write(iu,'(A,I0,A,I0,A)') '  fringe cells interpolated from fringe cells : ', n_ff, ' (of which from boundary fringes: ', n_ffo, ')'
    write(iu,'(A,I0,A,I0,A)') '  solid-bounded active cells (solid holes act as the wall) : ', n_wallb, &
         ' (covered only by forced fringes of a finer grid: ', n_gapfill, ')'
    write(iu,'(A,I0,A,I0,A,I0,A,I0)') '  coverage pass: uncovered fluid cells made active ', n_cov_self, &
         ', coarse cells made active under uncovered forced layers ', n_cov_other, &
         ', forced cells with a node covered by nothing ', n_cov_none, &
         ', uncovered cells crossed by a thinner body (left to the finer grid) ', n_cov_thin
    if (ovst_only /= 0) write(iu,'(A,I0,A,I0,A)') '  permeation check (no solved cell inside a body): ', n_perm_pt, &
         ' active cells with a test point inside a body (', n_perm_cen, ' with the centre inside)'
    close(iu)
    write(*,'(A)') '  wrote overset report: '//trim(casename)//'.ovst'
    write(*,'(A,F8.2,A)') '  assembly time   : ', tsec, ' s'
    write(*,'(A)') ' --------------------------------------------------------'
  end subroutine write_report

end module overset
