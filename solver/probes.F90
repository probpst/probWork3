!=======================================================================
! probes.F90 - volume-field probe points (update 22)
!
!   A probe point exports the time history of every volume variable at
!   a point: u v w [m/s], TTL_Prs and Dym_Prs = p - p_hyd [Pa], VoF
!   (two-phase), Vapor (cavitation), Temperature [K] (energy), the
!   shear stress magnitude [Pa] (surface_output.cell_shear) and the
!   cell that holds the point (block i j k, overset status 1 active /
!   0 fringe).  Points:
!     $Probe_background      Probe_point<n> = x y z   (grid frame [m] at t = 0)
!       follow the BACKGROUND: fixed in the grid frame, or moving with
!       the background blocks (Moving_BGR = 1, the published pose
!       bgr_Rg / bgr_O0 / bgr_dg of mod_global)
!     inside a Body / Child_body / Sub_child_body group of $SixDOF
!       Probe_point<n> = x y z   follow the PART's 6-DoF motion
!       (x = Rg (x0 - O0) + O0 + dg of the part, as the force-model
!       points): a point on a propeller rotates with it.
!   Values: trilinear interpolation in the cell that contains the
!   point (located every export step in the current, moved grid: the
!   last cell's neighbourhood first, then the bucket search of every
!   block) of the node values of the volume files (face / cell averages
!   with the hole cells excluded).  With overlapping grids the ACTIVE
!   cell wins over a fringe cell; a point in a hole, outside the grid or
!   in a cell that no rank holds writes "NoData" for the step.
!   Files (rank 0, every Data_Intvl steps, Output_data/):
!     <case>_probe<n>.dat                       background point n
!     <case>_body1_probe<n>.dat                 point n of Body 1
!     <case>_body1_child1_probe<n>.dat          point n of Child_body 1 of Body 1
!     <case>_body1_child1_sub1_probe<n>.dat     point n of Sub_child_body 1 of it
!   columns: step time x y z (current position, grid frame = earth
!   frame unless the background moves) u v w TTL_Prs Dym_Prs [VoF]
!   [Vapor] [Temperature] Shear block i j k status
!=======================================================================
module probes
  use global_vars
  use mpi_vars
  use mpi
  use ovs_geom, only: bkt_t, build_buckets, locate_point, inv_trilinear
  use sixdof, only: sixdof_on_flag, sixdof_node_of, sixdof_node_pose, sixdof_part_stem, sixdof_part_name, prdef, nprdef
  use surface_output, only: cell_shear
  use turbulence, only: velocity_gradient_compact, strain_vorticity, mu_lam
  implicit none
  private
  public :: probe_init, probe_export, probe_finalize, probe_count
  ! point sampler shared with the statistics (update 23): pt_locate finds the candidate cells
  ! of a point (one per block), pt_owner_key the best candidate this rank owns, pt_sample the
  ! node-interpolated values of the fixed layout below (all dimensional)
  public :: pt_locate, pt_owner_key, pt_sample, node_sample, cell_grad_sample, NSAMP, IS_U, IS_P, IS_PD, IS_TAU, IS_ALPHA, IS_VAP, &
            IS_T, IS_NUT, IS_TK, IS_TW, IS_G, IS_NU, IS_DX
  integer, parameter :: NSAMP = 24
  integer, parameter :: IS_U = 1, IS_P = 4, IS_PD = 5, IS_TAU = 6, IS_ALPHA = 7, IS_VAP = 8, IS_T = 9, IS_NUT = 10, &
                        IS_TK = 11, IS_TW = 12, IS_G = 13, IS_NU = 22, IS_DX = 23
  ! IS_U..IS_U+2 u v w [m/s]; IS_P total pressure [Pa]; IS_PD p - p_hyd [Pa]; IS_TAU shear stress magnitude [Pa];
  ! IS_ALPHA liquid fraction; IS_VAP vapour fraction; IS_T temperature [K]; IS_NUT eddy viscosity [m2/s];
  ! IS_TK model k [m2/s2] (SA: nu_tilde [m2/s]); IS_TW model omega [1/s] (k-eps: eps [m2/s3]);
  ! IS_G..IS_G+8 velocity gradient g(c,l) = du_c/dx_l [1/s], index IS_G + (c-1) + 3(l-1), from the face
  !   differences of the staggered velocities (velocity_gradient_compact, update 23);
  ! IS_NU laminar kinematic viscosity [m2/s]; IS_DX local cell size (volume^(1/3)) [m]

  type :: probe_t
     integer  :: node = 0, n = 0                     ! owner node (0 = background), point number
     real(dp) :: x0(3) = 0.0_dp                       ! position at t = 0, grid frame [m]
     real(dp) :: x(3) = 0.0_dp                        ! current position [m]
     character(len=96) :: stem = ' '
     integer  :: iu = -1
     integer  :: lb = 0, li = 0, lj = 0, lk = 0       ! last cell (replicated search: identical on every rank)
  end type probe_t
  type(probe_t), allocatable :: pr(:)
  integer :: npr = 0, nval = 0
  type(bkt_t), allocatable :: bk(:)
  logical, allocatable :: bk_ok(:)
  integer, allocatable :: bk_step(:)
  character(len=*), parameter :: out_dir = 'Output_data/'

contains

  integer function probe_count()
    probe_count = npr
  end function probe_count

  !---------------------------------------------------------------------
  subroutine probe_init
    integer :: i, k, iu, ierr
    character(len=256) :: fname
    character(len=16) :: a
    npr = nprobe_bg + nprdef
    if (npr == 0) return
    allocate(pr(npr))
    do i = 1, nprobe_bg
       pr(i)%node = 0;  pr(i)%n = probe_bg_n(i);  pr(i)%x0 = probe_bg_x(:,i)
       write(a,'(I0)') pr(i)%n
       pr(i)%stem = 'probe'//trim(a)
    end do
    do i = 1, nprdef
       k = sixdof_node_of(prdef(i)%b, prdef(i)%c, prdef(i)%s)
       if (k == 0) then
          if (myrank == 0) write(*,'(A,I0,A,I0,A,I0,A)') ' ERROR (probes): Probe_point of Body ', prdef(i)%b, ' / Child_body ', &
             prdef(i)%c, ' / Sub_child_body ', prdef(i)%s, ': the part does not exist (csv / 6-DoF off)'
          call MPI_Abort(MPI_COMM_WORLD, 22, ierr)
       end if
       pr(nprobe_bg + i)%node = k;  pr(nprobe_bg + i)%n = prdef(i)%n;  pr(nprobe_bg + i)%x0 = prdef(i)%x
       write(a,'(I0)') prdef(i)%n
       pr(nprobe_bg + i)%stem = sixdof_part_stem(k)//'_probe'//trim(a)
    end do
    ! duplicate names (the same point number twice in one owner)
    do i = 2, npr
       do k = 1, i - 1
          if (pr(i)%stem == pr(k)%stem) then
             if (myrank == 0) write(*,'(A)') ' ERROR (probes): the probe point '//trim(pr(i)%stem)//' is defined twice'
             call MPI_Abort(MPI_COMM_WORLD, 22, ierr)
          end if
       end do
    end do
    nval = 6                                         ! u v w TTL_Prs Dym_Prs + [VoF] [Vapor] [Temperature] + Shear
    if (Num_of_Phase == 2) nval = nval + 1
    if (cav_on) nval = nval + 1
    if (energy_on) nval = nval + 1
    if (myrank == 0) then
       call make_dir(out_dir)
       do i = 1, npr
          fname = out_dir//trim(casename)//'_'//trim(pr(i)%stem)//'.dat'
          if (run_mode == 1) then
             open(newunit=iu, file=trim(fname), status='unknown', position='append', action='write')
          else
             open(newunit=iu, file=trim(fname), status='replace', action='write')
             write(iu,'(A)') '# probCFD volume-field probe point (update 22): '//trim(pr(i)%stem)
             if (pr(i)%node == 0) then
                write(iu,'(A,3ES14.6,A)') '#   follows the background grid; position at t = 0:', pr(i)%x0, ' m (grid frame)'
             else
                write(iu,'(A,3ES14.6,A)') '#   follows the 6-DoF motion of '//sixdof_part_name(pr(i)%node)// &
                   '; position at t = 0:', pr(i)%x0, ' m (grid frame)'
             end if
             write(iu,'(A)') '#   values: trilinear interpolation of the node values of the volume files in the cell that holds the'
             write(iu,'(A)') '#   point (active cell preferred over a fringe cell); "NoData" when the point lies in a hole or outside the grid'
             write(iu,'(A)') '# columns: 1 step  2 time[s]  3-5 x y z [m] (current position, grid frame)  6-8 u v w [m/s]'// &
                '  9 TTL_Prs [Pa]  10 Dym_Prs = p - p_hyd [Pa]'//merge('  11 VoF', '        ', Num_of_Phase == 2)// &
                merge('  Vapor', '       ', cav_on)//merge('  Temperature [K]', '                 ', energy_on)// &
                '  Shear_Stress [Pa]  block i j k  status (1 active, 0 fringe)'
          end if
          pr(i)%iu = iu
       end do
       write(*,'(A,I0,A)') '  probe points (update 22): ', npr, ' (files Output_data/<case>_<owner>_probe<n>.dat)'
       do i = 1, npr
          if (pr(i)%node == 0) then
             write(*,'(A,3F10.4,A)') '    '//trim(pr(i)%stem)//' at', pr(i)%x0, ' m, following the background'
          else
             write(*,'(A,3F10.4,A)') '    '//trim(pr(i)%stem)//' at', pr(i)%x0, ' m, following '//sixdof_part_name(pr(i)%node)
          end if
       end do
    end if
  end subroutine probe_init

  ! create a sub-folder (rank 0) when it does not exist
  subroutine make_dir(path)
    character(len=*), intent(in) :: path
    logical :: ex
    inquire(file=trim(path)//'.', exist=ex)
    if (.not. ex) call execute_command_line('mkdir -p '//trim(path), wait=.true.)
  end subroutine make_dir

  !---------------------------------------------------------------------
  ! export of step istep (every Data_Intvl steps and step 0)
  !---------------------------------------------------------------------
  subroutine probe_export(istep, t)
    integer,  intent(in) :: istep
    real(dp), intent(in) :: t
    integer  :: i, ncand, cand_b(8), cand_i(8), cand_j(8), cand_k(8), key, keymax, ierr, dwin, qwin, st
    real(dp) :: cand_r(3,8), vals(NSAMP + 5), vloc(NSAMP + 5), smp(NSAMP), Rg(3,3), O0(3), dg(3), vG(3), omG(3), R(3,3), eta(6), nu(6)
    real(dp) :: out(16)
    integer  :: nv
    if (npr == 0) return
    if (mod(istep, data_interval) /= 0 .and. istep /= 0) return
    do i = 1, npr
       ! current position
       if (pr(i)%node > 0) then
          call sixdof_node_pose(pr(i)%node, Rg, O0, dg, vG, omG, R, eta, nu)
          pr(i)%x = matmul(Rg, pr(i)%x0 - O0) + O0 + dg
       else if (bgr_on) then
          pr(i)%x = matmul(bgr_Rg, pr(i)%x0 - bgr_O0) + bgr_O0 + bgr_dg
       else
          pr(i)%x = pr(i)%x0
       end if
       call pt_locate(pr(i)%x, istep, pr(i)%lb, pr(i)%li, pr(i)%lj, pr(i)%lk, ncand, cand_b, cand_i, cand_j, cand_k, cand_r)
       call pt_owner_key(ncand, cand_b, cand_i, cand_j, cand_k, key, dwin, qwin, st)
       vloc = 0.0_dp
       if (key > 0) then
          call pt_sample(dwin, cand_i(qwin), cand_j(qwin), cand_k(qwin), cand_r(:,qwin), smp)
          vloc(1:NSAMP) = smp
          vloc(NSAMP+1) = real(cand_b(qwin), dp);  vloc(NSAMP+2) = real(cand_i(qwin), dp)
          vloc(NSAMP+3) = real(cand_j(qwin), dp);  vloc(NSAMP+4) = real(cand_k(qwin), dp);  vloc(NSAMP+5) = real(st, dp)
       end if
       call MPI_Allreduce(key, keymax, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
       if (key /= keymax .or. keymax == 0) vloc = 0.0_dp
       call MPI_Allreduce(vloc, vals, NSAMP + 5, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
       if (keymax > 0) then
          pr(i)%lb = nint(vals(NSAMP+1));  pr(i)%li = nint(vals(NSAMP+2));  pr(i)%lj = nint(vals(NSAMP+3));  pr(i)%lk = nint(vals(NSAMP+4))
       else
          pr(i)%lb = 0
       end if
       if (myrank == 0) then
          if (keymax > 0) then
             ! file columns: u v w, TTL_Prs, Dym_Prs, [VoF], [Vapor], [Temperature], Shear
             out(1:5) = vals(IS_U:IS_PD);  nv = 5
             if (Num_of_Phase == 2) then
                nv = nv + 1;  out(nv) = vals(IS_ALPHA)
             end if
             if (cav_on) then
                nv = nv + 1;  out(nv) = vals(IS_VAP)
             end if
             if (energy_on) then
                nv = nv + 1;  out(nv) = vals(IS_T)
             end if
             nv = nv + 1;  out(nv) = vals(IS_TAU)
             write(pr(i)%iu,'(I8,ES16.8,3ES15.7,20ES15.7)', advance='no') istep, t, pr(i)%x, out(1:nv)
             write(pr(i)%iu,'(I5,3I6,I4)') nint(vals(NSAMP+1)), nint(vals(NSAMP+2)), nint(vals(NSAMP+3)), nint(vals(NSAMP+4)), &
                nint(vals(NSAMP+5))
          else
             write(pr(i)%iu,'(I8,ES16.8,3ES15.7,A)') istep, t, pr(i)%x, '   NoData'
          end if
          flush(pr(i)%iu)
       end if
    end do
  end subroutine probe_export

  !---------------------------------------------------------------------
  ! candidate cells of point x: one per block that contains it (the grids are replicated, so
  ! every rank finds the same cells); the last cell (lb, li, lj, lk) is tried first
  !---------------------------------------------------------------------
  subroutine pt_locate(x, istep, lb, li, lj, lk, ncand, cb, ci, cj, ck, cr)
    real(dp), intent(in)  :: x(3)
    integer,  intent(in)  :: istep, lb, li, lj, lk
    integer,  intent(out) :: ncand, cb(8), ci(8), cj(8), ck(8)
    real(dp), intent(out) :: cr(3,8)
    integer  :: b, ic, jc, kc
    real(dp) :: rst(3)
    logical  :: found
    if (.not. allocated(bk)) then
       allocate(bk(nblk), bk_ok(nblk), bk_step(nblk));  bk_ok = .false.;  bk_step = -1
    end if
    ncand = 0
    do b = 1, nblk
       found = .false.
       if (lb == b) call near_search(b, li, lj, lk, x, ic, jc, kc, rst, found)
       if (.not. found) call bucket_search(b, istep, x, ic, jc, kc, rst, found)
       if (found .and. ncand < 8) then
          ncand = ncand + 1
          cb(ncand) = b;  ci(ncand) = ic;  cj(ncand) = jc;  ck(ncand) = kc;  cr(:,ncand) = rst
       end if
    end do
  end subroutine pt_locate

  !---------------------------------------------------------------------
  ! the best candidate this rank owns: key = 0 none; otherwise the owner domain dwin, the
  ! candidate qwin and its status st (active 1 / fringe 0); the key is comparable across ranks
  ! (active > fringe > lower block > lower rank)
  !---------------------------------------------------------------------
  subroutine pt_owner_key(ncand, cb, ci, cj, ck, key, dwin, qwin, st)
    integer, intent(in)  :: ncand, cb(8), ci(8), cj(8), ck(8)
    integer, intent(out) :: key, dwin, qwin, st
    integer :: q, ib, d, s2, k2
    key = 0;  dwin = 0;  qwin = 0;  st = -1
    do q = 1, ncand
       do ib = 1, nown
          d = own_list(ib)
          if (dom(d)%b /= cb(q)) cycle
          if (ci(q) < dom(d)%r1(1) .or. ci(q) > dom(d)%r2(1)) cycle
          if (cj(q) < dom(d)%r1(2) .or. cj(q) > dom(d)%r2(2)) cycle
          if (ck(q) < dom(d)%r1(3) .or. ck(q) > dom(d)%r2(3)) cycle
          call set_block(d)
          s2 = cstat(ci(q), cj(q), ck(q))
          if (s2 < 0) cycle
          k2 = merge(2, 1, s2 > 0)*1000000 + (nblk - cb(q))*1000 + (nprocs - myrank)
          if (k2 > key) then
             key = k2;  dwin = d;  qwin = q;  st = s2
          end if
       end do
    end do
  end subroutine pt_owner_key

  ! the cell holding p near the last cell (i0, j0, k0) of block b: +-2 cells
  subroutine near_search(b, i0, j0, k0, p, ic, jc, kc, rst, found)
    integer,  intent(in)  :: b, i0, j0, k0
    real(dp), intent(in)  :: p(3)
    integer,  intent(out) :: ic, jc, kc
    real(dp), intent(out) :: rst(3)
    logical,  intent(out) :: found
    integer  :: i, j, k, r
    real(dp) :: rr(3), q, qbest
    logical  :: ok
    found = .false.;  qbest = huge(1.0_dp);  ic = 0;  jc = 0;  kc = 0;  rst = 0.0_dp
    do r = 0, 2
       do k = max(k0-r, 1), min(k0+r, blk(b)%nk - 1)
        do j = max(j0-r, 1), min(j0+r, blk(b)%nj - 1)
         do i = max(i0-r, 1), min(i0+r, blk(b)%ni - 1)
            if (max(abs(i-i0), abs(j-j0), abs(k-k0)) /= r) cycle
            call inv_trilinear(blk(b)%xn, blk(b)%yn, blk(b)%zn, blk(b)%ni, blk(b)%nj, blk(b)%nk, i, j, k, p, rr, ok)
            if (.not. ok) cycle
            if (any(rr < -1.0e-7_dp) .or. any(rr > 1.0_dp + 1.0e-7_dp)) cycle
            q = maxval(abs(rr - 0.5_dp))
            if (q < qbest) then
               qbest = q;  found = .true.;  ic = i;  jc = j;  kc = k;  rst = min(max(rr, 0.0_dp), 1.0_dp)
            end if
         end do
        end do
       end do
       if (found) return
    end do
  end subroutine near_search

  ! bucket search of block b (the bucket is rebuilt when the block has moved since its build)
  subroutine bucket_search(b, istep, p, ic, jc, kc, rst, found)
    integer,  intent(in)  :: b, istep
    real(dp), intent(in)  :: p(3)
    integer,  intent(out) :: ic, jc, kc
    real(dp), intent(out) :: rst(3)
    logical,  intent(out) :: found
    logical :: moving
    moving = (blk_obj(b) /= 0 .and. sixdof_on_flag()) .or. (blk_obj(b) == 0 .and. bgr_on)
    if (.not. bk_ok(b) .or. (moving .and. bk_step(b) /= istep)) then
       call build_buckets(blk(b)%xn, blk(b)%yn, blk(b)%zn, blk(b)%ni, blk(b)%nj, blk(b)%nk, bk(b))
       bk_ok(b) = .true.;  bk_step(b) = istep
    end if
    call locate_point(bk(b), blk(b)%xn, blk(b)%yn, blk(b)%zn, p, ic, jc, kc, rst, found)
  end subroutine bucket_search

  !---------------------------------------------------------------------
  ! values at a point: trilinear interpolation (r s t of the inverse map) of the node values of
  ! cell (i, j, k) of domain d (set_block(d) is done here); the layout of the module header
  !---------------------------------------------------------------------
  subroutine pt_sample(d, i, j, k, r, v)
    integer,  intent(in)  :: d, i, j, k
    real(dp), intent(in)  :: r(3)
    real(dp), intent(out) :: v(NSAMP)
    integer  :: a, c, e
    real(dp) :: wn, vn(NSAMP), zn, phyd, z
    call set_block(d)
    v = 0.0_dp;  z = 0.0_dp
    do e = 0, 1
     do c = 0, 1
      do a = 0, 1
         wn = merge(r(1), 1.0_dp - r(1), a == 1)*merge(r(2), 1.0_dp - r(2), c == 1)*merge(r(3), 1.0_dp - r(3), e == 1)
         call node_sample(i+a-1, j+c-1, k+e-1, vn)
         v = v + wn*vn
         z = z + wn*blk(cur_blk)%zn(i+a, j+c, k+e)
      end do
     end do
    end do
    ! dynamic pressure: p - p_hyd(z) at the point's height (the Dym_Prs rule of the volume files)
    zn = z/Lch
    phyd = -gnd(3)*(zfs - zn)*merge(1.0_dp, rho_ratio, zn <= zfs)*(rho1*Uch*Uch)
    v(IS_PD) = v(IS_P) - phyd - merge(p_level, 0.0_dp, cav_on)
  end subroutine pt_sample

  !---------------------------------------------------------------------
  ! the group-2 quantities of ONE cell (update 23, type-1 statistics at the cells): the compact
  ! gradient tensor, nu_t, k, omega, nu, Delta of cell (i, j, k) of the current block, dimensional
  !---------------------------------------------------------------------
  subroutine cell_grad_sample(i, j, k, vn)
    integer,  intent(in)  :: i, j, k
    real(dp), intent(out) :: vn(NSAMP)
    real(dp) :: g(3,3)
    integer :: cc, l
    vn = 0.0_dp
    vn(IS_U)   = 0.5_dp*(u(i-1,j,k) + u(i,j,k))*Uch                 ! cell-centred velocity (for k_res)
    vn(IS_U+1) = 0.5_dp*(v(i,j-1,k) + v(i,j,k))*Uch
    vn(IS_U+2) = 0.5_dp*(w(i,j,k-1) + w(i,j,k))*Uch
    call velocity_gradient_compact(i, j, k, g)
    do l = 1, 3
       do cc = 1, 3
          vn(IS_G + (cc-1) + 3*(l-1)) = g(cc,l)*Uch/Lch
       end do
    end do
    vn(IS_NU) = mu_lam(i, j, k)*mu1/rho1
    vn(IS_DX) = Jh(2*i, 2*j, 2*k)**(1.0_dp/3.0_dp)*Lch
    if (turb_on) vn(IS_NUT) = mut(i,j,k)*mu1/rho1
    if (TURB_model == 1) then
       vn(IS_TK) = tk(i,j,k)*mu1/rho1
    else if (TURB_model >= 2) then
       vn(IS_TK) = tk(i,j,k)*Uch*Uch
       if (TURB_model == 2) then
          vn(IS_TW) = tw(i,j,k)*Uch**3/Lch
       else
          vn(IS_TW) = tw(i,j,k)*Uch/Lch
       end if
    end if
  end subroutine cell_grad_sample

  !---------------------------------------------------------------------
  ! node values at node (im+1, jm+1, km+1) of the current block (set_block done): the averages
  ! of the volume files (faces / cells around the node, holes excluded), dimensional
  !---------------------------------------------------------------------
  subroutine node_sample(im, jm, km, vn)
    integer,  intent(in)  :: im, jm, km
    real(dp), intent(out) :: vn(NSAMP)
    integer  :: a, c, e, n, ng2, cc, l
    real(dp) :: s, g(3,3), sg(3,3), smag, wmag, snu, sdx, smut, stk, stw
    vn = 0.0_dp
    vn(IS_U)   = favg(1, im, jm, km)*Uch
    vn(IS_U+1) = favg(2, im, jm, km)*Uch
    vn(IS_U+2) = favg(3, im, jm, km)*Uch
    vn(IS_P)   = cavg(p, im, jm, km)*(rho1*Uch*Uch)
    vn(IS_ALPHA) = 1.0_dp
    if (Num_of_Phase == 2) vn(IS_ALPHA) = cavg(alpha, im, jm, km)
    if (cav_on)    vn(IS_VAP) = cavg(vap, im, jm, km)
    if (energy_on) vn(IS_T)   = cavg(temp, im, jm, km)
    ! cell-based quantities whose stencil must lie inside the ghost range: shear, gradient, nu
    s = 0.0_dp;  sg = 0.0_dp;  snu = 0.0_dp;  sdx = 0.0_dp;  smut = 0.0_dp;  stk = 0.0_dp;  stw = 0.0_dp;  n = 0;  ng2 = 0
    do e = 0, 1
     do c = 0, 1
      do a = 0, 1
         if (cstat(im+a, jm+c, km+e) < 0) cycle
         snu = snu + mu_lam(im+a, jm+c, km+e);  sdx = sdx + Jh(2*(im+a), 2*(jm+c), 2*(km+e))
         if (turb_on) then
            smut = smut + mut(im+a, jm+c, km+e)
            if (TURB_model > 0) then
               stk = stk + tk(im+a, jm+c, km+e);  stw = stw + tw(im+a, jm+c, km+e)
            end if
         end if
         n = n + 1
         if (im+a-1 < is-ng .or. im+a+1 > ie+ng .or. jm+c-1 < js-ng .or. jm+c+1 > je+ng .or. &
             km+e-1 < ks-ng .or. km+e+1 > ke+ng) cycle
         call velocity_gradient_compact(im+a, jm+c, km+e, g)                ! face differences (update 23)
         call strain_vorticity(g, smag, wmag)
         s = s + cell_shear(im+a, jm+c, km+e)
         sg = sg + g;  ng2 = ng2 + 1
      end do
     end do
    end do
    if (n > 0) then
       vn(IS_NU)  = snu/real(n, dp)*mu1/rho1
       vn(IS_DX)  = (sdx/real(n, dp))**(1.0_dp/3.0_dp)*Lch
       if (turb_on) vn(IS_NUT) = smut/real(n, dp)*mu1/rho1
       if (TURB_model == 1) then
          vn(IS_TK) = stk/real(n, dp)*mu1/rho1                      ! nu_tilde [m2/s]
       else if (TURB_model >= 2) then
          vn(IS_TK) = stk/real(n, dp)*Uch*Uch                       ! k [m2/s2]
          if (TURB_model == 2) then
             vn(IS_TW) = stw/real(n, dp)*Uch**3/Lch                 ! epsilon [m2/s3]
          else
             vn(IS_TW) = stw/real(n, dp)*Uch/Lch                    ! omega [1/s]
          end if
       end if
    end if
    if (ng2 > 0) then
       vn(IS_TAU) = s/real(ng2, dp)
       do l = 1, 3
          do cc = 1, 3
             vn(IS_G + (cc-1) + 3*(l-1)) = sg(cc,l)/real(ng2, dp)*Uch/Lch
          end do
       end do
    end if
  end subroutine node_sample

  ! Cartesian component q of the 4 faces of family q around the node, hole faces excluded
  real(dp) function favg(q, im, jm, km)
    integer, intent(in) :: q, im, jm, km
    integer :: a, c, id(3), e1(3), e2(3), t1, t2, n
    real(dp) :: s
    real(dp), pointer :: f(:,:,:)
    select case (q)
    case (1);  f => u;  t1 = 2;  t2 = 3
    case (2);  f => v;  t1 = 1;  t2 = 3
    case default;  f => w;  t1 = 1;  t2 = 2
    end select
    e1 = 0;  e1(t1) = 1;  e2 = 0;  e2(t2) = 1
    s = 0.0_dp;  n = 0
    do a = 0, 1
       do c = 0, 1
          id = [im, jm, km] + a*e1 + c*e2
          if (fstat(q, id(1), id(2), id(3)) /= -1) then
             s = s + f(id(1),id(2),id(3));  n = n + 1
          end if
       end do
    end do
    favg = 0.0_dp
    if (n > 0) favg = s/real(n, dp)
  end function favg

  ! cell field over the 8 cells around the node, holes excluded
  real(dp) function cavg(f, im, jm, km)
    real(dp), intent(in) :: f(is-ng:, js-ng:, ks-ng:)
    integer,  intent(in) :: im, jm, km
    integer :: a, c, e, n
    real(dp) :: s
    s = 0.0_dp;  n = 0
    do e = 0, 1
     do c = 0, 1
      do a = 0, 1
         if (cstat(im+a, jm+c, km+e) >= 0) then
            s = s + f(im+a, jm+c, km+e);  n = n + 1
         end if
      end do
     end do
    end do
    cavg = 0.0_dp
    if (n > 0) cavg = s/real(n, dp)
  end function cavg

  subroutine probe_finalize
    integer :: i
    if (myrank /= 0 .or. npr == 0) return
    do i = 1, npr
       if (pr(i)%iu > 0 .or. pr(i)%iu < -1) close(pr(i)%iu)
    end do
  end subroutine probe_finalize
end module probes
