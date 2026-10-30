!=======================================================================
! b2b_exchange.F90
!   Runtime block-to-block (b2b) coupling for the fully-staggered
!   multi-block solver.
!
!   b2b_setup    : reads 'casename.bcs_b2b' (written at start-up by
!                  make_bcs_b2b_file), converts every record into the
!                  affine index maps used by the exchange, decides the
!                  ownership of the shared interface-face velocity
!                  dofs, builds the interface masks, and overwrites
!                  ghost metric entries with the donor block values.
!
!   build_xrecs  : expands the block-level patches into DOMAIN-level
!                  exchange records and adds the intra-block halo
!                  records between the MPI sub-domains of each block
!                  (identity index maps).  Every ghost update -
!                  intra-block halo and block-to-block interface -
!                  runs through the same record engine.
!
!   refresh_all  : ghost update of one field over all sub-domains =
!                  record exchange (3 axis passes, repeated three times
!                  when any interface permutes the axes) followed by the
!                  physical boundary conditions.  Works for every field
!                  id of the table in mod_global (velocities, pressure,
!                  volume fraction, Krylov vectors).
!
!   Staggering & orientation
!   ------------------------
!   Receiver index r along axis a maps to donor index along donor
!   axis q(a):
!       node index :  d = sgn(a)*r + noff(a)
!       cell index :  d = sgn(a)*r + coff(a),  coff = noff + (sgn-1)/2
!       face index :  d = sgn(a)*r + foff(a),  foff = noff + sgn - 1
!   A velocity component along receiver axis c corresponds to the
!   donor velocity component along axis q(c) multiplied by sgn(c)
!   (grid-direction reversal flips the component sign).  Face-type
!   index maps are used along the component axis, cell-type maps
!   along the other two - this is the "interpolation" that makes the
!   staggered data land on the correct physical locations even when
!   an i-face meets a k-face.
!
!   Shared interface-face dofs (e.g. u on a shared i-face): exactly
!   one side (the "owner", see own_key) computes them in the momentum
!   and corrector steps; the other side receives the value through
!   the exchange, so the dof stays single-valued like in the
!   equivalent single-block grid.
!=======================================================================
module b2b_runtime
  use mpi
  use global_vars
  use mpi_vars
  use b2b_mod, only: b2b_read_file, n_b2b_tot, b2b_blk, b2b_dir,        &
                     b2b_irng, b2b_jrng, b2b_krng, b2b_dblk, b2b_ddir,  &
                     b2b_dirng, b2b_djrng, b2b_dkrng, b2b_map
  use overset, only: ovs_update_cells
  implicit none
  private
  public :: b2b_setup, refresh_all, refresh_uvw, exchange_faces

  type :: buf_t
     real(dp), allocatable :: a(:)
  end type buf_t

contains

  !=====================================================================
  ! Build the runtime patch table from 'casename.bcs_b2b'.
  !=====================================================================
  subroutine b2b_setup
    integer  :: nb_file, ierr, ierr_mpi
    real(dp) :: tol_file
    integer  :: m, a, qa, b, na
    integer  :: R0(3), D0(3), noff
    integer  :: t1, t2, c1s, c1e, c2s, c2e, ftyp
    integer  :: nwarn

    call b2b_read_file(trim(casename)//'.bcs_b2b', nb_file, tol_file, ierr)
    if (ierr /= 0) call MPI_Abort(MPI_COMM_WORLD, 20, ierr_mpi)
    if (nb_file /= nblk) then
       if (myrank == 0) write(*,'(A,I0,A,I0)')                           &
          ' ERROR: bcs_b2b file has ', nb_file, ' blocks but grid has ', nblk
       call MPI_Abort(MPI_COMM_WORLD, 21, ierr_mpi)
    end if

    nxp = n_b2b_tot + nper        ! auto-detected b2b + periodic (flag 3)
    allocate(xp(nxp))
    b2b_permuted = .false.

    do m = 1, n_b2b_tot           ! auto-detected b2b patches (from file)
       xp(m)%rblk  = b2b_blk(m)
       xp(m)%dblk  = b2b_dblk(m)
       xp(m)%na    = abs(b2b_dir(m))
       xp(m)%rside = merge(1, 2, b2b_dir(m) > 0)      ! + = min face
       xp(m)%rn(:,1) = b2b_irng(:,m)
       xp(m)%rn(:,2) = b2b_jrng(:,m)
       xp(m)%rn(:,3) = b2b_krng(:,m)
       R0 = xp(m)%rn(1,:)
       D0 = [ b2b_dirng(1,m), b2b_djrng(1,m), b2b_dkrng(1,m) ]
       do a = 1, 3
          qa = abs(b2b_map(a,m))
          xp(m)%q(a)   = qa
          xp(m)%sgn(a) = sign(1, b2b_map(a,m))
          if (qa /= a) b2b_permuted = .true.
          noff = D0(qa) - xp(m)%sgn(a)*R0(a)
          xp(m)%coff(a) = noff + (xp(m)%sgn(a) - 1)/2
          xp(m)%foff(a) = noff +  xp(m)%sgn(a) - 1
       end do
       ! ownership of the shared face dofs: the side with the smaller
       ! key owns; preference max face > min face, then block number,
       ! then face axis (deterministic and identical for both mirror
       ! records of an interface)
       xp(m)%r_owner = own_key(b2b_blk(m),  b2b_dir(m)) <               &
                       own_key(b2b_dblk(m), b2b_ddir(m))
    end do

    !------------------- periodic patches (flag 3 in the bcs file) ------
    ! A periodic pairing is the same ghost-exchange coupling as a b2b
    ! interface (the two boundary faces are IDENTIFIED, sharing the
    ! face-normal velocity dofs), except the faces are not coincident
    ! in space.  The user-given pair ranges define the affine index
    ! map; the axis pairing/signs were resolved in read_bcs.
    do a = 1, nper
       m = n_b2b_tot + a
       xp(m)%rblk  = pbc(a)%rblk
       xp(m)%dblk  = pbc(a)%dblk
       xp(m)%na    = abs(pbc(a)%dirs)
       xp(m)%rside = merge(1, 2, pbc(a)%dirs > 0)
       xp(m)%rn    = pbc(a)%rn
       do qa = 1, 3
          xp(m)%q(qa)   = pbc(a)%q(qa)
          xp(m)%sgn(qa) = pbc(a)%sgn(qa)
          if (xp(m)%q(qa) /= qa) b2b_permuted = .true.
          noff = pbc(a)%dn(1, pbc(a)%q(qa)) - pbc(a)%sgn(qa)*pbc(a)%rn(1,qa)
          xp(m)%coff(qa) = noff + (xp(m)%sgn(qa) - 1)/2
          xp(m)%foff(qa) = noff +  xp(m)%sgn(qa) - 1
       end do
       xp(m)%r_owner = own_key(pbc(a)%rblk, pbc(a)%dirs) <               &
                       own_key(pbc(a)%dblk, pbc(a)%ddir)
    end do
    ! every periodic patch needs its mirror (the pair side must also be
    ! configured in the bcs file, or its ghosts would never be filled)
    do a = 1, nper
       nwarn = 0
       do m = 1, nper
          if (pbc(m)%rblk == pbc(a)%dblk .and. pbc(m)%dirs == pbc(a)%ddir &
              .and. pbc(m)%dblk == pbc(a)%rblk                            &
              .and. pbc(m)%ddir == pbc(a)%dirs) nwarn = 1
       end do
       if (nwarn == 0 .and. myrank == 0) write(*,'(A,I0,A,I0,A)')         &
          '  WARNING (periodic): block ', pbc(a)%rblk, ' face dir ',      &
          pbc(a)%dirs, ' has no mirror periodic surface on the pair'//    &
          ' block - its ghosts will not be updated'
    end do

    !--------------------------- interface-face masks -------------------
    do b = 1, nblk
       allocate(blk(b)%mu_lo(blk(b)%ny, blk(b)%nz));  blk(b)%mu_lo = .false.
       allocate(blk(b)%mu_hi(blk(b)%ny, blk(b)%nz));  blk(b)%mu_hi = .false.
       allocate(blk(b)%mv_lo(blk(b)%nx, blk(b)%nz));  blk(b)%mv_lo = .false.
       allocate(blk(b)%mv_hi(blk(b)%nx, blk(b)%nz));  blk(b)%mv_hi = .false.
       allocate(blk(b)%mw_lo(blk(b)%nx, blk(b)%ny));  blk(b)%mw_lo = .false.
       allocate(blk(b)%mw_hi(blk(b)%nx, blk(b)%ny));  blk(b)%mw_hi = .false.
    end do
    do b = 1, nblk
       allocate(blk(b)%tu_lo(blk(b)%ny, blk(b)%nz));  blk(b)%tu_lo = 1
       allocate(blk(b)%tu_hi(blk(b)%ny, blk(b)%nz));  blk(b)%tu_hi = 1
       allocate(blk(b)%tv_lo(blk(b)%nx, blk(b)%nz));  blk(b)%tv_lo = 1
       allocate(blk(b)%tv_hi(blk(b)%nx, blk(b)%nz));  blk(b)%tv_hi = 1
       allocate(blk(b)%tw_lo(blk(b)%nx, blk(b)%ny));  blk(b)%tw_lo = 1
       allocate(blk(b)%tw_hi(blk(b)%nx, blk(b)%ny));  blk(b)%tw_hi = 1
    end do
    ! interface faces (both sides) are type 0
    do m = 1, nxp
       b  = xp(m)%rblk
       na = xp(m)%na
       select case (na)
       case (1);  t1 = 2;  t2 = 3
       case (2);  t1 = 1;  t2 = 3
       case (3);  t1 = 1;  t2 = 2
       end select
       c1s = xp(m)%rn(1,t1);  c1e = xp(m)%rn(2,t1) - 1
       c2s = xp(m)%rn(1,t2);  c2e = xp(m)%rn(2,t2) - 1
       select case (na)
       case (1)
          if (xp(m)%rside == 1) then
             blk(b)%tu_lo(c1s:c1e, c2s:c2e) = 0
          else
             blk(b)%tu_hi(c1s:c1e, c2s:c2e) = 0
          end if
       case (2)
          if (xp(m)%rside == 1) then
             blk(b)%tv_lo(c1s:c1e, c2s:c2e) = 0
          else
             blk(b)%tv_hi(c1s:c1e, c2s:c2e) = 0
          end if
       case (3)
          if (xp(m)%rside == 1) then
             blk(b)%tw_lo(c1s:c1e, c2s:c2e) = 0
          else
             blk(b)%tw_hi(c1s:c1e, c2s:c2e) = 0
          end if
       end select
    end do
    ! pressure-outlet patches (Dirichlet p + zero-gradient normal
    ! velocity) are type 2; mass-corrected outflow patches (normal
    ! velocity flag 2, Neumann p) are type 3
    do b = 1, nblk
       do m = 1, blk(b)%nbc
          na = abs(blk(b)%bc_dir(m))
          if (blk(b)%bc_flag(4,m) == 1 .and. blk(b)%bc_eflag(na,m) /= 1) then
             ftyp = 2
          else if (blk(b)%bc_flag(4,m) == 0 .and. blk(b)%bc_eflag(na,m) == 2) then
             ftyp = 3
          else
             cycle
          end if
          select case (na)
          case (1);  t1 = 2;  t2 = 3
          case (2);  t1 = 1;  t2 = 3
          case (3);  t1 = 1;  t2 = 2
          end select
          if (na == 1) then
             c1s = blk(b)%bc_jrng(1,m);  c1e = blk(b)%bc_jrng(2,m) - 1
             c2s = blk(b)%bc_krng(1,m);  c2e = blk(b)%bc_krng(2,m) - 1
             if (blk(b)%bc_dir(m) > 0) then
                blk(b)%tu_lo(c1s:c1e, c2s:c2e) = ftyp
             else
                blk(b)%tu_hi(c1s:c1e, c2s:c2e) = ftyp
             end if
          else if (na == 2) then
             c1s = blk(b)%bc_irng(1,m);  c1e = blk(b)%bc_irng(2,m) - 1
             c2s = blk(b)%bc_krng(1,m);  c2e = blk(b)%bc_krng(2,m) - 1
             if (blk(b)%bc_dir(m) > 0) then
                blk(b)%tv_lo(c1s:c1e, c2s:c2e) = ftyp
             else
                blk(b)%tv_hi(c1s:c1e, c2s:c2e) = ftyp
             end if
          else
             c1s = blk(b)%bc_irng(1,m);  c1e = blk(b)%bc_irng(2,m) - 1
             c2s = blk(b)%bc_jrng(1,m);  c2e = blk(b)%bc_jrng(2,m) - 1
             if (blk(b)%bc_dir(m) > 0) then
                blk(b)%tw_lo(c1s:c1e, c2s:c2e) = ftyp
             else
                blk(b)%tw_hi(c1s:c1e, c2s:c2e) = ftyp
             end if
          end if
       end do
    end do
    !------------------ shared-face ownership masks ---------------------
    do m = 1, nxp
       if (.not. xp(m)%r_owner) cycle
       b  = xp(m)%rblk
       na = xp(m)%na
       select case (na)
       case (1);  t1 = 2;  t2 = 3
       case (2);  t1 = 1;  t2 = 3
       case (3);  t1 = 1;  t2 = 2
       end select
       c1s = xp(m)%rn(1,t1);  c1e = xp(m)%rn(2,t1) - 1    ! cell windows
       c2s = xp(m)%rn(1,t2);  c2e = xp(m)%rn(2,t2) - 1
       select case (na)
       case (1)
          if (xp(m)%rside == 1) then
             blk(b)%mu_lo(c1s:c1e, c2s:c2e) = .true.
          else
             blk(b)%mu_hi(c1s:c1e, c2s:c2e) = .true.
          end if
       case (2)
          if (xp(m)%rside == 1) then
             blk(b)%mv_lo(c1s:c1e, c2s:c2e) = .true.
          else
             blk(b)%mv_hi(c1s:c1e, c2s:c2e) = .true.
          end if
       case (3)
          if (xp(m)%rside == 1) then
             blk(b)%mw_lo(c1s:c1e, c2s:c2e) = .true.
          else
             blk(b)%mw_hi(c1s:c1e, c2s:c2e) = .true.
          end if
       end select
    end do

    !------------- expand to domain-level exchange records --------------
    call build_xrecs

    if (myrank == 0) then
       write(*,'(A,I0,A,I0,A)') '  b2b interfaces  : ', nxp,             &
            ' patch records active (of which periodic: ', nper, ')'
       write(*,'(A,I0,A)') '  exchange records: ', nxr,                  &
            ' (intra-block halos + b2b, domain level)'
       if (b2b_permuted) write(*,'(A)')                                  &
          '  b2b orientation : permuted axis maps present (triple exchange sweep)'
    end if
  end subroutine b2b_setup

  !=====================================================================
  ! Expand the block-level b2b patches into per-(receiver domain,
  ! donor domain) records and generate the intra-block halo records.
  ! All ranks build the identical table (deterministic, from the
  ! replicated domain metadata) so message tags and sizes agree.
  !=====================================================================
  subroutine build_xrecs
    integer :: cap, m, rd, dd, b, a, dir, iv
    integer :: g1, g2, c1, c2, dn, step, nc
    integer :: tr(2,3)
    logical :: ok, any_ok
    type(xrec_t) :: r

    cap = 64;  allocate(xr(cap));  nxr = 0

    !---------------- kind 1: block-to-block sub-records ---------------
    do m = 1, nxp
       do rd = blk(xp(m)%rblk)%dom0, blk(xp(m)%rblk)%dom0 + blk(xp(m)%rblk)%ncpu - 1
          ! receiver sub-domain must touch the interface face
          if (xp(m)%rside == 1) then
             if (dom(rd)%r1(xp(m)%na) /= 1) cycle
          else
             if (dom(rd)%r2(xp(m)%na) /= ncell(xp(m)%rblk, xp(m)%na)) cycle
          end if
          do dd = blk(xp(m)%dblk)%dom0, blk(xp(m)%dblk)%dom0 + blk(xp(m)%dblk)%ncpu - 1
             r = xrec_t(kind=1, rd=rd, dd=dd, ip=m, na=xp(m)%na)
             do a = 1, 3
                call dom_supply(dd, xp(m)%q(a), r%ds(1,xp(m)%q(a)),      &
                                                r%ds(2,xp(m)%q(a)))
             end do
             ! keep the record only if any field class has work
             any_ok = .false.
             call push_rec(r, cap)               ! provisional
             do iv = 1, 4
                call xr_window(nxr, iv, tr, ok)
                if (ok) any_ok = .true.
             end do
             if (.not. any_ok) nxr = nxr - 1     ! drop it again
          end do
       end do
    end do

    !---------------- kind 2: intra-block halo records ------------------
    do rd = 1, ndom
       b = dom(rd)%b
       do a = 1, 3
          do dir = -1, 1, 2
             nc = ncell(b, a)
             if (dir == -1 .and. dom(rd)%r1(a) == 1)  cycle   ! block face
             if (dir ==  1 .and. dom(rd)%r2(a) == nc) cycle   ! block face
             if (dir == 1) then
                g1 = dom(rd)%r2(a) + 1;  g2 = dom(rd)%r2(a) + ng
             else
                g1 = dom(rd)%r1(a) - ng;  g2 = dom(rd)%r1(a) - 1
             end if
             ! walk neighboring sub-domains until the layers are covered
             dn   = dom(rd)%pc(a) + dir
             step = merge(g1, g2, dir == 1)      ! next uncovered layer
             do while ( (dir ==  1 .and. step <= g2) .or.                &
                        (dir == -1 .and. step >= g1) )
                if (dn < 0 .or. dn > blk(b)%sp(a)-1) exit
                dd = dom_at(b, dom(rd)%pc, a, dn)
                r  = xrec_t(kind=2, rd=rd, dd=dd, na=a)
                call dom_supply(dd, a, c1, c2)
                if (dir == 1) then
                   r%nrm1 = max(g1, dom(dd)%r1(a));  r%nrm2 = min(g2, c2)
                else
                   r%nrm1 = max(g1, c1);  r%nrm2 = min(g2, dom(dd)%r2(a))
                end if
                if (r%nrm1 <= r%nrm2) then
                   call push_rec(r, cap)
                   if (dir == 1) then
                      step = r%nrm2 + 1
                   else
                      step = r%nrm1 - 1
                   end if
                else
                   exit
                end if
                dn = dn + dir
             end do
          end do
       end do
    end do
  end subroutine build_xrecs

  !--------------------------------------------------------------------
  ! Domain id of the sub-domain of block b at the split-grid position
  ! obtained from 'pc' with coordinate 'a' replaced by 'val'.
  !--------------------------------------------------------------------
  pure integer function dom_at(b, pc, a, val)
    integer, intent(in) :: b, pc(3), a, val
    integer :: q(3)
    q = pc;  q(a) = val
    dom_at = blk(b)%dom0 + q(1) + blk(b)%sp(1)*(q(2) + blk(b)%sp(2)*q(3))
  end function dom_at

  !--------------------------------------------------------------------
  subroutine push_rec(r, cap)
    type(xrec_t), intent(in)    :: r
    integer,      intent(inout) :: cap
    type(xrec_t), allocatable :: tmp(:)
    if (nxr == cap) then
       allocate(tmp(2*cap))
       tmp(1:nxr) = xr(1:nxr)
       call move_alloc(tmp, xr)
       cap = 2*cap
    end if
    nxr = nxr + 1
    xr(nxr) = r
  end subroutine push_rec

  !--------------------------------------------------------------------
  ! Total ordering key used to pick the owner side of an interface.
  !--------------------------------------------------------------------
  pure integer function own_key(b, dir)
    integer, intent(in) :: b, dir
    integer :: side
    side = merge(0, 1, dir < 0)              ! prefer max faces
    own_key = side*1000000 + b*10 + abs(dir)
  end function own_key


  !=====================================================================
  ! Ghost refresh of field iv (1..3 = U, 4 = p, 5..7 = Us, 8 = alpha,
  ! 9/10 = Krylov vectors) on all blocks: b2b exchange in three axis
  ! passes (three sweeps when permuted interfaces exist, so edge/corner
  ! ghosts settle), then the physical boundary conditions of every
  ! owned block (apply_bc receives the field id itself).
  !=====================================================================
  ! Optional reps > 1 repeats the whole (exchange + BC) cycle: edge and
  ! corner ghost cells that are relayed through a third block or through
  ! a donor's BC-set ghosts are one cycle stale after a single cycle
  ! (harmless inside fixed-point iterations), and the repetition makes
  ! the ghost state an exact, iteration-independent function of the
  ! owned values - required for the Krylov pressure solver, whose
  ! operator must be strictly linear.
  subroutine refresh_all(iv, reps)
    integer, intent(in) :: iv
    integer, intent(in), optional :: reps
    integer :: sw, ax, ib, b, base, rep, nrep
    real(dp), pointer :: f(:,:,:)

    nrep = 1
    if (present(reps)) nrep = reps
    base = iv                              ! apply_bc maps the id itself
    do rep = 1, nrep
       do sw = 1, merge(3, 1, b2b_permuted)
          do ax = 1, 3
             call b2b_pass(iv, ax)
          end do
       end do
       do ib = 1, nown
          b = own_list(ib)
          call set_block(b)
          f => field_ptr(b, iv)
          call apply_bc(f, base)
       end do
       if (iv == 20 .and. associated(temp_bc_hook)) call temp_bc_hook()   ! T: Dirichlet / flux / convective ghosts (update 16)
    end do
    ! overset: donor interpolation of the cell field at the fringe cells
    ! (alpha is written into the field; p / kv1 are kept for the
    ! implicit fringe rows of the pressure solver)
    if (ovst_on) call ovs_update_cells(iv)
  end subroutine refresh_all

  !=====================================================================
  ! Combined ghost refresh of a velocity TRIPLET: ivbase = 1 refreshes
  ! (u, v, w), ivbase = 5 refreshes (us, vs, ws).  All three components
  ! are exchanged inside the same sweep/pass structure.  This is
  ! REQUIRED for b2b interfaces with permuted axis maps: there, the
  ! ghost values of one component are sourced from a DIFFERENT
  ! component of the donor block (e.g. u ghosts from the donor's w),
  ! so per-component refreshes would relay data that is one refresh
  ! cycle stale through intermediate blocks.  For identity-oriented
  ! interfaces this is equivalent to three separate refreshes.
  !=====================================================================
  subroutine refresh_uvw(ivbase)
    integer, intent(in) :: ivbase       ! 1 -> u,v,w ; 5 -> us,vs,ws
    integer :: sw, ax, ib, b, c, base
    real(dp), pointer :: f(:,:,:)

    do sw = 1, merge(3, 1, b2b_permuted)
       do ax = 1, 3
          do c = 0, 2
             call b2b_pass(ivbase + c, ax)
          end do
       end do
    end do

    base = ivbase                          ! apply_bc maps the id itself
    do ib = 1, nown
       b = own_list(ib)
       call set_block(b)
       do c = 0, 2
          f => field_ptr(b, ivbase + c)
          call apply_bc(f, base + c)
       end do
    end do
  end subroutine refresh_uvw

  !=====================================================================
  ! Exchange-only refresh of the face-type triplet ivbase..ivbase+2
  ! (no boundary conditions): used for the pressure-gradient fluxes
  ! G^q (ids 11..13) so that non-owner interface faces and sub-domain
  ! halo faces receive the value computed by the owning side.
  !=====================================================================
  subroutine exchange_faces(ivbase)
    integer, intent(in) :: ivbase
    integer :: sw, ax, c
    do sw = 1, merge(3, 1, b2b_permuted)
       do ax = 1, 3
          do c = 0, 2
             call b2b_pass(ivbase + c, ax)
          end do
       end do
    end do
  end subroutine exchange_faces

  !=====================================================================
  ! One exchange pass: all records whose receiver normal axis is
  ! 'axis' (intra-block halos and b2b interfaces alike).  Same-rank
  ! pairs copy directly; cross-rank pairs use non-blocking messages
  ! (tag = record index).
  !=====================================================================
  ! Two sub-passes per call: (1) every record with its tangential
  ! window WIDENED by ng - reaches the edge/corner ghosts that no patch
  ! covers directly; (2) the block-to-block records again with their
  ! UN-WIDENED core windows, overwriting.  Where several patches of one
  ! face meet (a face that bends from one neighbour onto another, e.g.
  ! an H-block wrapping an O-grid), the widened data of a patch is only
  ! the geometric continuation of ITS donor and may be wrong beyond
  ! the junction; the core data of the patch that actually covers the
  ! cell is authoritative and wins.  The same rule builds the ghost
  ! coordinates (curvi.F90), so metrics and fields stay consistent.
  subroutine b2b_pass(iv, axis)
    integer, intent(in) :: iv, axis
    integer :: sub
    do sub = 1, 2
       call b2b_pass_sub(iv, axis, sub == 2)
    end do
  end subroutine b2b_pass

  subroutine b2b_pass_sub(iv, axis, core)
    integer, intent(in) :: iv, axis
    logical, intent(in) :: core
    type(buf_t), allocatable :: sb(:), rb(:)
    integer, allocatable :: reqs(:), upk(:)
    integer, allocatable :: trs(:,:,:)
    integer :: m, ro, do_, tr(2,3), sz, nreq, nupk, ierr, i
    logical :: ok

    allocate(sb(nxr), rb(nxr), reqs(2*nxr), upk(nxr), trs(2,3,nxr))
    nreq = 0;  nupk = 0

    do m = 1, nxr
       if (xr(m)%na /= axis) cycle
       if (core .and. xr(m)%kind /= 1) cycle
       ro  = dom(xr(m)%rd)%owner
       do_ = dom(xr(m)%dd)%owner
       if (myrank /= ro .and. myrank /= do_) cycle
       call xr_window(m, iv, tr, ok, core)
       if (.not. ok) cycle
       trs(:,:,m) = tr
       sz = (tr(2,1)-tr(1,1)+1)*(tr(2,2)-tr(1,2)+1)*(tr(2,3)-tr(1,3)+1)

       ! ALL packing happens before ANY unpacking, and the unpacking is
       ! performed in ascending record order: every rank distribution
       ! (including np = 1) then produces bit-identical ghost values,
       ! because packs read the pass-start state and overlapping windows
       ! are resolved by the same (record-id) write order everywhere.
       if (myrank == do_) then
          allocate(sb(m)%a(sz))
          call xr_pack(m, iv, tr, sb(m)%a)
          if (ro /= myrank) then
             nreq = nreq + 1
             call MPI_Isend(sb(m)%a, sz, MPI_DOUBLE_PRECISION, ro, m,    &
                            MPI_COMM_WORLD, reqs(nreq), ierr)
          end if
       end if
       if (myrank == ro) then
          if (do_ /= myrank) then
             allocate(rb(m)%a(sz))
             nreq = nreq + 1
             call MPI_Irecv(rb(m)%a, sz, MPI_DOUBLE_PRECISION, do_, m,   &
                            MPI_COMM_WORLD, reqs(nreq), ierr)
          end if
          nupk = nupk + 1
          upk(nupk) = m
       end if
    end do

    if (nreq > 0) call MPI_Waitall(nreq, reqs, MPI_STATUSES_IGNORE, ierr)
    do i = 1, nupk
       m = upk(i)
       if (dom(xr(m)%dd)%owner == myrank) then
          call xr_unpack(m, iv, trs(:,:,m), sb(m)%a)
       else
          call xr_unpack(m, iv, trs(:,:,m), rb(m)%a)
       end if
    end do
    deallocate(sb, rb, reqs, upk, trs)
  end subroutine b2b_pass_sub

  !=====================================================================
  ! Receiver-side target index window of exchange record m for field iv.
  !
  ! kind 1 (block-to-block): the block-level ghost-window logic of the
  !   parent patch (normal direction: ghost cells, or the shared face
  !   plus ghost faces for the face-type component, honoring the dof
  !   ownership; tangential: patch node window widened by ng with the
  !   face/cell index distinction), clipped to the receiver DOMAIN
  !   storage and to receiver indices whose mapped donor indices stay
  !   inside the donor DOMAIN's supply range.
  !
  ! kind 2 (intra-block halo, identity map): normal direction = the
  !   ghost layers nrm1..nrm2 assigned to this donor at setup (same
  !   index range for every variable class, replicating the original
  !   single-block halo exchange); tangential = the common sub-domain
  !   extent widened by ng, clipped to both storages.
  !=====================================================================
  subroutine xr_window(m, iv, tr, ok, core)
    integer, intent(in)  :: m, iv
    integer, intent(out) :: tr(2,3)
    logical, intent(out) :: ok
    logical, intent(in), optional :: core   ! .true.: tangential window
                                            ! NOT widened (authoritative
                                            ! data of the patch itself)
    integer :: comp, a, nr, off, tlo, thi, s1, s2, qa
    integer :: rd, dd, ip, wid

    rd = xr(m)%rd;  dd = xr(m)%dd;  ip = xr(m)%ip
    comp = comp_axis(iv)
    ok = .true.
    wid = ng
    if (present(core)) then
       if (core) wid = 0
    end if

    do a = 1, 3
       if (xr(m)%kind == 2) then
          !--------------- intra-block halo (identity) ------------------
          if (a == xr(m)%na) then
             tlo = xr(m)%nrm1;  thi = xr(m)%nrm2
          else
             tlo = max(dom(rd)%r1(a), dom(dd)%r1(a)) - ng
             thi = min(dom(rd)%r2(a), dom(dd)%r2(a)) + ng
             ! clip to both storages
             tlo = max(tlo, dom(rd)%r1(a)-ng, dom(dd)%r1(a)-ng)
             thi = min(thi, dom(rd)%r2(a)+ng, dom(dd)%r2(a)+ng)
          end if
       else
          !--------------- block-to-block sub-record --------------------
          nr = ncell(xp(ip)%rblk, a)
          if (a == xp(ip)%na) then
             if (comp == a) then                       ! face-type normal
                ! the shared face plane itself (index 0 / nr) is written
                ! only by the CORE window of the patch covering it, and
                ! only when this side does not own the dof; the widened
                ! sub-pass fills the outer ghost faces only
                if (xp(ip)%rside == 1) then
                   tlo = -(ng-1)
                   thi = merge(-1, 0, xp(ip)%r_owner .or. wid > 0)
                else
                   tlo = merge(nr+1, nr, xp(ip)%r_owner .or. wid > 0)
                   thi = nr + (ng-1)
                end if
             else                                      ! ghost cells
                if (xp(ip)%rside == 1) then
                   tlo = 1 - ng;  thi = 0
                else
                   tlo = nr + 1;  thi = nr + ng
                end if
             end if
          else                                         ! tangential
             if (comp == a) then                       ! face-type index
                tlo = xp(ip)%rn(1,a) - 1 - wid
                thi = xp(ip)%rn(2,a) - 1 + wid
             else                                      ! cell-type index
                tlo = xp(ip)%rn(1,a) - wid
                thi = xp(ip)%rn(2,a) - 1 + wid
             end if
          end if

          ! clip to the receiver domain storage
          tlo = max(tlo, dom(rd)%r1(a) - ng)
          thi = min(thi, dom(rd)%r2(a) + ng)

          ! clip so the mapped donor indices stay inside the donor
          ! domain's supply range
          qa  = xp(ip)%q(a)
          off = merge(xp(ip)%foff(a), xp(ip)%coff(a), comp == a)
          s1  = xr(m)%ds(1,qa);  s2 = xr(m)%ds(2,qa)
          if (xp(ip)%sgn(a) > 0) then
             tlo = max(tlo, s1 - off)
             thi = min(thi, s2 - off)
          else
             tlo = max(tlo, off - s2)
             thi = min(thi, off - s1)
          end if
       end if

       if (tlo > thi) ok = .false.
       tr(1,a) = tlo;  tr(2,a) = thi
    end do
  end subroutine xr_window

  !--------------------------------------------------------------------
  ! Donor field id and component sign of b2b patch ip for receiver
  ! field iv (kind-1 records; kind-2 uses div = iv, factor 1).  Cell-
  ! centered fields (p, alpha, Krylov vectors) map onto themselves.
  !--------------------------------------------------------------------
  subroutine xp_dvar(ip, iv, div, fac)
    integer, intent(in)  :: ip, iv
    integer, intent(out) :: div
    real(dp), intent(out):: fac
    integer :: comp
    comp = comp_axis(iv)
    if (comp == 0) then
       div = iv                                 ! cell field -> same field
       fac = 1.0_dp
    else
       div = xp(ip)%q(comp)
       if (iv > 10) then
          div = div + 10                        ! pressure-gradient fluxes
       else if (iv > 4) then
          div = div + 4                         ! star velocities family
       end if
       fac = real(xp(ip)%sgn(comp), dp)
    end if
  end subroutine xp_dvar

  !=====================================================================
  ! Pack donor data of record m for receiver window tr into buf
  ! (canonical r1-fastest order).
  !=====================================================================
  subroutine xr_pack(m, iv, tr, buf)
    integer,  intent(in)  :: m, iv, tr(2,3)
    real(dp), intent(out) :: buf(:)
    real(dp), pointer :: fd(:,:,:)
    integer  :: comp, div, off(3), ddx(3), r1, r2, r3, n, a, ip
    real(dp) :: fac

    if (xr(m)%kind == 2) then
       !------------- intra-block: identity copy -------------------------
       fd => field_ptr(xr(m)%dd, iv)
       n = 0
       do r3 = tr(1,3), tr(2,3)
          do r2 = tr(1,2), tr(2,2)
             do r1 = tr(1,1), tr(2,1)
                n = n + 1
                buf(n) = fd(r1, r2, r3)
             end do
          end do
       end do
    else
       !------------- block-to-block: affine index + component map -------
       ip = xr(m)%ip
       comp = comp_axis(iv)
       call xp_dvar(ip, iv, div, fac)
       fd => field_ptr(xr(m)%dd, div)
       do a = 1, 3
          off(a) = merge(xp(ip)%foff(a), xp(ip)%coff(a), comp == a)
       end do
       n = 0
       do r3 = tr(1,3), tr(2,3)
          ddx(xp(ip)%q(3)) = xp(ip)%sgn(3)*r3 + off(3)
          do r2 = tr(1,2), tr(2,2)
             ddx(xp(ip)%q(2)) = xp(ip)%sgn(2)*r2 + off(2)
             do r1 = tr(1,1), tr(2,1)
                ddx(xp(ip)%q(1)) = xp(ip)%sgn(1)*r1 + off(1)
                n = n + 1
                buf(n) = fac * fd(ddx(1), ddx(2), ddx(3))
             end do
          end do
       end do
    end if
  end subroutine xr_pack

  !=====================================================================
  ! Unpack buf into the receiver domain field over window tr.
  !=====================================================================
  subroutine xr_unpack(m, iv, tr, buf)
    integer,  intent(in) :: m, iv, tr(2,3)
    real(dp), intent(in) :: buf(:)
    real(dp), pointer :: fr(:,:,:)
    integer :: r1, r2, r3, n

    fr => field_ptr(xr(m)%rd, iv)
    n = 0
    do r3 = tr(1,3), tr(2,3)
       do r2 = tr(1,2), tr(2,2)
          do r1 = tr(1,1), tr(2,1)
             n = n + 1
             fr(r1, r2, r3) = buf(n)
          end do
       end do
    end do
  end subroutine xr_unpack

end module b2b_runtime
