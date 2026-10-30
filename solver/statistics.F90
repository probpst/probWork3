!=======================================================================
! statistics.F90 - on-the-fly turbulence statistics (update 23)
!
!   $Statistics (SW_statistics = 1): time statistics of the volume
!   fields over the window Statistics_Start_T .. Statistics_End_T (0 =
!   End_T), one sample every Statistics_Sample_Intvl steps, written
!   every Statistics_Output_Intvl steps (0 = at the end only) and at the
!   end of the run into Statistics_volume_solutions/.
!   Method (raw moments, the fluctuation statistics formed at output):
!     mean  = S/n,  rms = sqrt(max(Q/n - mean^2, 0)),  <x'y'> = C/n - mean_x mean_y
!   with per-point counters n that count only the steps at which the
!   point was valid (not a hole), a second counter for the gradient
!   terms (stencil fully valid); points never sampled are written as 0
!   with Coverage = 0.  The accumulators are in the restart file.
!   Where:
!     Statistics_type1 = 1  the original grid blocks (node based, written
!       on the t = 0 grid <case>_stat_grid.grd); with Statistics_Type1_Frame
!       = 1 the vectors / tensors of a moving block are turned into its
!       part's frame (R^T u, R^T G R: a propeller block gives blade-fixed
!       statistics), with 0 the grid axes
!     Statistics_type2/3/4 = 1  user-defined boxes / cylinders / planes:
!       structured lattices attached to a frame (Mother_body = 0: the
!       background; > 0: that part, with Mother_child_body /
!       Mother_sub_child_body), sampled with the probe-point machinery of
!       update 22 (trilinear interpolation in the active cell of the
!       moved grid); vectors in the frame's axes, a cylinder with
!       Statistics_cyl_Polar_Vel = 1 in axial / radial / tangential
!       components
!   Group 1 (Statistics_first_group): Coverage, the mean and RMS of
!     u v w TTL_Prs Dym_Prs Shear_Stress [VoF] [Vapor] [Temperature]
!     [nu_t] [k_model] [omega/eps_model / nu_tilde], the Reynolds
!     stresses uu vv ww uv vw wu and k_res
!   Group 2 (Statistics_second_group, velocity only): Coverage_grad, the 9
!     mean velocity gradients, the 9 RMS fluctuating gradients, the 3
!     cross terms <du'/dy dv'/dx> ..., eps_res = nu <(du'_i/dx_j +
!     du'_j/dx_i)^2>/2, with a model eps_model / eps_total / k_model /
!     k_total, the Kolmogorov scale eta, the Taylor scale lambda_g, Re_lambda,
!     Delta/eta and the mean vorticity
!=======================================================================
module statistics
  use global_vars
  use mpi_vars
  use mpi
  use sixdof, only: sixdof_on_flag, sixdof_node_of, sixdof_node_pose, sixdof_part_name
  use probes, only: pt_locate, pt_owner_key, pt_sample, node_sample, cell_grad_sample, NSAMP, IS_U, IS_P, IS_PD, IS_TAU, IS_ALPHA, IS_VAP, &
                    IS_T, IS_NUT, IS_TK, IS_TW, IS_G, IS_NU, IS_DX
  use node_groups
  implicit none
  private
  public :: stat_on, stat_init, stat_sample, stat_write, stat_finalize, stat_restart_write, stat_restart_read, stat_summary

  logical :: stat_on = .false.
  character(len=*), parameter :: sdir = 'Statistics_volume_solutions/'
  real(dp), parameter :: bstar = 0.09_dp
  ! ---- field list of group 1 (sample indices, names, units) ----
  integer :: nf1 = 0
  integer :: f1idx(12)
  character(len=16) :: f1name(12), f1unit(12)
  ! ---- accumulator layout (per point): group 1 then group 2 ----
  integer :: na1 = 0, na2 = 0, nacc = 0, na2c = 0       ! sizes (na2c = the cell accumulator of type 1: 6 + na2)
  integer :: iS1, iQ1, iC1                            ! sums, squares (nf1 each), cross uv vw wu (3)
  integer :: iS2, iQ2, iX2, iEM, iKM, iNU, iDX        ! gradient sums (9), squares (9), cross (3), eps_mod, k_mod, nu, dx
  ! ---- type 1: node accumulators per owned domain ----
  type :: dacc_t
     real(dp), allocatable :: a(:,:,:,:)              ! (nacc, i, j, k) node range r1..r2+1 (group 1 at the nodes)
     real(dp), allocatable :: n1(:,:,:), n2(:,:,:)    ! counts (n2 unused at the nodes since update 23)
     real(dp), allocatable :: c(:,:,:,:)              ! (na2, i, j, k) cell range r1..r2: group 2 at the CELLS
     real(dp), allocatable :: nc(:,:,:)               ! cell counts of group 2
  end type dacc_t
  type(dacc_t), allocatable :: da(:)                  ! (ndom), allocated for owned domains
  ! ---- lattices ----
  type :: latt_t
     integer  :: kind = 0, id = 0, node = 0, npt = 0, n(3) = 1, polar = 0
     character(len=32) :: stem = ' '
     real(dp), allocatable :: x0(:,:)                  ! (3, npt) grid frame at t = 0
     real(dp), allocatable :: pb(:,:,:)                ! (3, 3, npt) local basis (columns) for the polar components
     real(dp), allocatable :: a(:,:), n1(:), n2(:)    ! accumulators (this rank's share)
     integer,  allocatable :: lb(:), li(:), lj(:), lk(:)
     integer,  allocatable :: key(:), keymax(:)
  end type latt_t
  type(latt_t), allocatable :: lt(:)
  integer :: nlt = 0
  integer :: nsmp = 0                                ! samples taken (all points share the schedule)
  real(dp) :: t_first = -1.0_dp, t_last = -1.0_dp
  logical :: grid_written = .false., sampled_any = .false.
  type :: g0_t
     real(dp), allocatable :: x(:,:,:,:)              ! (3, ni, nj, nk) the t = 0 grid (rank 0)
  end type g0_t
  type(g0_t), allocatable :: g0(:)
  type :: nbf_t
     real(dp), allocatable :: f(:,:,:,:)
  end type nbf_t
  logical, allocatable :: lt_grid_written(:)
  integer :: last_written = -1
  real(dp), allocatable :: stash_v(:,:)

contains

  !=====================================================================
  subroutine stat_init
    integer :: ib, d, i, k, ierr, n1n, n2n, n3n
    character(len=16) :: a
    stat_on = (sw_statistics == 1)
    if (.not. stat_on) return
    if (stat_group(1) == 0 .and. stat_group(2) == 0) then
       if (myrank == 0) write(*,'(A)') '  NOTE (statistics): both groups are off (Statistics_first_group = Statistics_second_group = 0)'
    end if
    ! ---- group-1 fields ----
    nf1 = 0
    call addf(IS_U,   'u',            'm/s');  call addf(IS_U+1, 'v', 'm/s');  call addf(IS_U+2, 'w', 'm/s')
    call addf(IS_P,   'TTL_Prs',      'Pa');   call addf(IS_PD,  'Dym_Prs', 'Pa');  call addf(IS_TAU, 'Shear_Stress', 'Pa')
    if (Num_of_Phase == 2) call addf(IS_ALPHA, 'VoF', '-')
    if (cav_on)    call addf(IS_VAP, 'Vapor', '-')
    if (energy_on) call addf(IS_T, 'Temperature', 'K')
    if (turb_on)   call addf(IS_NUT, 'nu_t', 'm2/s')
    select case (TURB_model)
    case (1);     call addf(IS_TK, 'nu_tilde', 'm2/s')
    case (2);     call addf(IS_TK, 'k_model', 'm2/s2');  call addf(IS_TW, 'eps_model', 'm2/s3')
    case (3, 4);  call addf(IS_TK, 'k_model', 'm2/s2');  call addf(IS_TW, 'omega_model', '1/s')
    end select
    ! ---- accumulator layout ----
    iS1 = 1;  iQ1 = iS1 + nf1;  iC1 = iQ1 + nf1;  na1 = 2*nf1 + 3
    iS2 = na1 + 1;  iQ2 = iS2 + 9;  iX2 = iQ2 + 9;  iEM = iX2 + 3;  iKM = iEM + 1;  iNU = iKM + 1;  iDX = iNU + 1
    na2 = 9 + 9 + 3 + 4;  nacc = na1 + na2;  na2c = 6 + na2
    ! ---- type 1 ----
    if (stat_type(1) == 1) then
       allocate(da(ndom))
       do ib = 1, nown
          d = own_list(ib)
          n1n = dom(d)%r2(1) - dom(d)%r1(1) + 2;  n2n = dom(d)%r2(2) - dom(d)%r1(2) + 2;  n3n = dom(d)%r2(3) - dom(d)%r1(3) + 2
          allocate(da(d)%a(nacc, n1n, n2n, n3n), da(d)%n1(n1n, n2n, n3n), da(d)%n2(n1n, n2n, n3n))
          da(d)%a = 0.0_dp;  da(d)%n1 = 0.0_dp;  da(d)%n2 = 0.0_dp
          allocate(da(d)%c(na2c, n1n-1, n2n-1, n3n-1), da(d)%nc(n1n-1, n2n-1, n3n-1))
          da(d)%c = 0.0_dp;  da(d)%nc = 0.0_dp
       end do
    end if
    ! ---- lattices ----
    nlt = 0
    do i = 1, nlatt_def
       if (stat_type(latt_def(i)%kind) == 1) nlt = nlt + 1
    end do
    allocate(lt(max(nlt, 1)), lt_grid_written(max(nlt, 1)))
    lt_grid_written = .false.
    nlt = 0
    do i = 1, nlatt_def
       if (stat_type(latt_def(i)%kind) /= 1) cycle
       nlt = nlt + 1
       call build_lattice(latt_def(i), lt(nlt))
    end do
    if (myrank == 0 .and. stat_type(1) == 1) then          ! the t = 0 grid for the type-1 files
       allocate(g0(nblk))
       do k = 1, nblk
          allocate(g0(k)%x(3, blk(k)%ni, blk(k)%nj, blk(k)%nk))
          g0(k)%x(1,:,:,:) = blk(k)%xn;  g0(k)%x(2,:,:,:) = blk(k)%yn;  g0(k)%x(3,:,:,:) = blk(k)%zn
       end do
    end if
    if (myrank == 0) then
       call make_dir(sdir)
       write(*,'(A)') '  statistics (update 23): ON'
       write(*,'(A,F10.4,A,F10.4,A,I0,A,I0,A)') '    window ', stat_start_t, ' .. ', merge(stat_end_t, End_T, stat_end_t > 0.0_dp), &
          ' s, one sample every ', stat_sample_intvl, ' steps, output every ', stat_output_intvl, ' steps (0 = end of the run)'
       write(*,'(A,I0,A,I0,A)') '    groups: first ', stat_group(1), ', second ', stat_group(2)
       if (stat_type(1) == 1) write(*,'(A)') '    type 1: the original grid blocks ('// &
          merge('part frames of the moving blocks', 'grid axes                       ', stat_type1_frame == 1)//')'
       do k = 1, nlt
          write(a,'(I0)') lt(k)%npt
          write(*,'(A)') '    '//trim(lt(k)%stem)//': '//trim(a)//' points, following '//trim(owner_name(lt(k)%node))
       end do
       write(*,'(A,I0,A)') '    group-1 fields: ', nf1, ' ('//trim(fields_line())//')'
    end if
    ierr = 0
  contains
    subroutine addf(idx, nm, un)
      integer, intent(in) :: idx
      character(len=*), intent(in) :: nm, un
      nf1 = nf1 + 1;  f1idx(nf1) = idx;  f1name(nf1) = nm;  f1unit(nf1) = un
    end subroutine addf
    function fields_line() result(s)
      character(len=256) :: s
      integer :: q
      s = ' '
      do q = 1, nf1
         s = trim(s)//' '//trim(f1name(q))
      end do
    end function fields_line
  end subroutine stat_init

  function owner_name(node) result(s)
    integer, intent(in) :: node
    character(len=64) :: s
    if (node > 0) then
       s = sixdof_part_name(node)
    else if (bgr_on) then
       s = 'the moving background'
    else
       s = 'the earth-fixed frame'
    end if
  end function owner_name

  subroutine make_dir(path)
    character(len=*), intent(in) :: path
    logical :: ex
    inquire(file=trim(path)//'.', exist=ex)
    if (.not. ex) call execute_command_line('mkdir -p '//trim(path), wait=.true.)
  end subroutine make_dir

  !---------------------------------------------------------------------
  ! lattice points of a box / cylinder / plane definition
  !---------------------------------------------------------------------
  subroutine build_lattice(df, L)
    type(latt_def_t), intent(in) :: df
    type(latt_t), intent(out) :: L
    integer  :: i, j, k, q, ierr
    real(dp) :: e1(3), e2(3), e3(3), s1, s2, s3, rr, th, er(3), et(3)
    character(len=16) :: a
    L%kind = df%kind;  L%id = df%id;  L%polar = df%polar
    write(a,'(I0)') df%id
    select case (df%kind)
    case (2);  L%stem = 'box'//trim(a)
    case (3);  L%stem = 'cyl'//trim(a)
    case default;  L%stem = 'plane'//trim(a)
    end select
    L%node = 0
    if (df%mb > 0) then
       L%node = sixdof_node_of(df%mb, df%mc, df%ms)
       if (L%node == 0) then
          if (myrank == 0) write(*,'(A,I0,A)') ' ERROR (statistics): Statistics_'//trim(L%stem)//': the mother part (Body ', df%mb, &
             ' ...) does not exist'
          call MPI_Abort(MPI_COMM_WORLD, 23, ierr)
       end if
    end if
    L%n = df%np
    if (df%kind == 4) L%n(3) = 1
    if (any(L%n < 1)) call lerr('nPoints must be >= 1')
    L%npt = L%n(1)*L%n(2)*L%n(3)
    allocate(L%x0(3, L%npt), L%pb(3, 3, L%npt))
    e1 = unit(df%ax1);  e2 = df%ax2 - dot_product(df%ax2, e1)*e1
    if (norm2(e2) < 1.0e-12_dp) call lerr('Axis1 and Axis2 must not be parallel')
    if (abs(dot_product(unit(df%ax2), e1)) > 1.0e-6_dp .and. df%kind /= 3) call lerr('Axis1 and Axis2 must be perpendicular')
    e2 = unit(e2);  e3 = cross(e1, e2)
    q = 0
    do k = 1, L%n(3)
     do j = 1, L%n(2)
      do i = 1, L%n(1)
         q = q + 1
         L%pb(:,:,q) = reshape([1.0_dp,0.0_dp,0.0_dp, 0.0_dp,1.0_dp,0.0_dp, 0.0_dp,0.0_dp,1.0_dp], [3,3])
         select case (df%kind)
         case (2, 4)
            s1 = frac(i, L%n(1))*df%len(1);  s2 = frac(j, L%n(2))*df%len(2)
            s3 = 0.0_dp
            if (df%kind == 2) s3 = frac(k, L%n(3))*df%len(3)
            L%x0(:,q) = df%org + s1*e1 + s2*e2 + s3*e3
         case (3)
            ! e1 = axis, e2 = reference direction (theta = 0), e3 = e1 x e2; i radial, j angular, k axial
            rr = df%rad(1) + frac(i, L%n(1))*(df%rad(2) - df%rad(1))
            th = (df%theta(1) + frac(j, L%n(2))*(df%theta(2) - df%theta(1)))*acos(-1.0_dp)/180.0_dp
            s3 = frac(k, L%n(3))*df%len(1)
            er = cos(th)*e2 + sin(th)*e3;  et = -sin(th)*e2 + cos(th)*e3
            L%x0(:,q) = df%org + s3*e1 + rr*er
            if (df%polar == 1) then
               L%pb(:,1,q) = e1;  L%pb(:,2,q) = er;  L%pb(:,3,q) = et       ! axial, radial, tangential
            end if
         end select
      end do
     end do
    end do
    allocate(L%a(nacc, L%npt), L%n1(L%npt), L%n2(L%npt), L%lb(L%npt), L%li(L%npt), L%lj(L%npt), L%lk(L%npt), &
             L%key(L%npt), L%keymax(L%npt))
    L%a = 0.0_dp;  L%n1 = 0.0_dp;  L%n2 = 0.0_dp;  L%lb = 0;  L%li = 0;  L%lj = 0;  L%lk = 0
  contains
    real(dp) function frac(i, n)
      integer, intent(in) :: i, n
      frac = 0.0_dp
      if (n > 1) frac = real(i-1, dp)/real(n-1, dp)
    end function frac
    subroutine lerr(msg)
      character(len=*), intent(in) :: msg
      integer :: ierr2
      if (myrank == 0) write(*,'(A)') ' ERROR (statistics): Statistics_'//trim(L%stem)//': '//msg
      call MPI_Abort(MPI_COMM_WORLD, 23, ierr2)
    end subroutine lerr
  end subroutine build_lattice

  pure function unit(v) result(u)
    real(dp), intent(in) :: v(3)
    real(dp) :: u(3), s
    s = norm2(v)
    u = v
    if (s > 0.0_dp) u = v/s
  end function unit
  pure function cross(a, b) result(c)
    real(dp), intent(in) :: a(3), b(3)
    real(dp) :: c(3)
    c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
  end function cross

  !---------------------------------------------------------------------
  ! rotation of the sample into a frame: u -> R^T u, G -> R^T G R (R = the frame's rotation
  ! in the grid axes), then into the local basis P (columns): P^T u, P^T G P
  !---------------------------------------------------------------------
  subroutine rotate_sample(v, R, P)
    real(dp), intent(inout) :: v(NSAMP)
    real(dp), intent(in) :: R(3,3), P(3,3)
    real(dp) :: M(3,3), G(3,3), u(3)
    integer :: c, l
    M = matmul(R, P)                                    ! combined basis: columns = frame axes in the grid axes
    u = matmul(transpose(M), v(IS_U:IS_U+2));  v(IS_U:IS_U+2) = u
    do l = 1, 3
       do c = 1, 3
          G(c,l) = v(IS_G + (c-1) + 3*(l-1))
       end do
    end do
    G = matmul(transpose(M), matmul(G, M))
    do l = 1, 3
       do c = 1, 3
          v(IS_G + (c-1) + 3*(l-1)) = G(c,l)
       end do
    end do
  end subroutine rotate_sample

  ! accumulate one sample into a(:)
  subroutine accumulate(v, a, n1, n2, do2)
    real(dp), intent(in) :: v(NSAMP)
    real(dp), intent(inout) :: a(nacc), n1, n2
    logical, intent(in) :: do2
    integer :: q, c, l
    real(dp) :: G(3,3), smag, s2, em, km, x
    n1 = n1 + 1.0_dp
    do q = 1, nf1
       x = v(f1idx(q))
       a(iS1 + q - 1) = a(iS1 + q - 1) + x
       a(iQ1 + q - 1) = a(iQ1 + q - 1) + x*x
    end do
    a(iC1)   = a(iC1)   + v(IS_U)*v(IS_U+1)
    a(iC1+1) = a(iC1+1) + v(IS_U+1)*v(IS_U+2)
    a(iC1+2) = a(iC1+2) + v(IS_U+2)*v(IS_U)
    if (stat_group(2) == 0 .or. .not. do2) return
    n2 = n2 + 1.0_dp
    do l = 1, 3
       do c = 1, 3
          G(c,l) = v(IS_G + (c-1) + 3*(l-1))
       end do
    end do
    do q = 1, 9
       a(iS2 + q - 1) = a(iS2 + q - 1) + v(IS_G + q - 1)
       a(iQ2 + q - 1) = a(iQ2 + q - 1) + v(IS_G + q - 1)**2
    end do
    a(iX2)   = a(iX2)   + G(1,2)*G(2,1)
    a(iX2+1) = a(iX2+1) + G(2,3)*G(3,2)
    a(iX2+2) = a(iX2+2) + G(3,1)*G(1,3)
    ! modelled dissipation / kinetic energy of the sample
    s2 = 0.0_dp
    do l = 1, 3
       do c = 1, 3
          s2 = s2 + (0.5_dp*(G(c,l) + G(l,c)))**2
       end do
    end do
    smag = sqrt(2.0_dp*s2)
    em = 0.0_dp;  km = 0.0_dp
    select case (TURB_model)
    case (2);     km = v(IS_TK);  em = v(IS_TW)
    case (3, 4);  km = v(IS_TK);  em = bstar*v(IS_TK)*v(IS_TW)
    case default
       if (turb_on) then                              ! SA or the WMLES branch: production = dissipation estimate
          em = v(IS_NUT)*smag**2;  km = v(IS_NUT)*smag/sqrt(bstar)
       end if
    end select
    a(iEM) = a(iEM) + em;  a(iKM) = a(iKM) + km
    a(iNU) = a(iNU) + v(IS_NU);  a(iDX) = a(iDX) + v(IS_DX)
  end subroutine accumulate

  ! the cell accumulator of type 1 (update 23): c(1:3) = sum u v w, c(4:6) = sum of squares (for k_res),
  ! c(7:na2c) = the group-2 layout a(iS2:nacc)
  subroutine accumulate2(v, c, n2)
    real(dp), intent(in) :: v(NSAMP)
    real(dp), intent(inout) :: c(na2c), n2
    real(dp) :: tmp(nacc), n1d
    c(1:3) = c(1:3) + v(IS_U:IS_U+2);  c(4:6) = c(4:6) + v(IS_U:IS_U+2)**2
    tmp = 0.0_dp;  tmp(iS2:nacc) = c(7:na2c);  n1d = 0.0_dp
    call accumulate(v, tmp, n1d, n2, .true.)
    c(7:na2c) = tmp(iS2:nacc)
  end subroutine accumulate2

  ! the group-2 statistics of one cell from its accumulator
  subroutine cell_stats(c, n, out, nv)
    real(dp), intent(in) :: c(na2c), n
    real(dp), intent(out) :: out(80)
    integer, intent(out) :: nv
    real(dp) :: tmp(nacc)
    tmp = 0.0_dp
    tmp(iS1:iS1+2) = c(1:3);  tmp(iQ1:iQ1+2) = c(4:6);  tmp(iS2:nacc) = c(7:na2c)
    call point_stats(tmp, n, n, 2, out, nv)
  end subroutine cell_stats

  !=====================================================================
  ! one time step: sample when inside the window
  !=====================================================================
  subroutine stat_sample(istep, t)
    integer,  intent(in) :: istep
    real(dp), intent(in) :: t
    integer  :: ib, d, b, i, j, k, in1, in2, in3, q, ierr, ncand, cb(8), ci(8), cj(8), ck(8), key, dwin, qwin, st, L
    real(dp) :: v(NSAMP), R(3,3), Id(3,3), cr(3,8), Rg(3,3), O0(3), dg(3), vG(3), omG(3), Rr(3,3), eta(6), nu(6), x(3)
    logical  :: ok1, ok2
    if (.not. stat_on) return
    if (istep == 0) return
    if (t < stat_start_t - 1.0e-12_dp) return
    if (stat_end_t > 0.0_dp .and. t > stat_end_t + 1.0e-12_dp) return
    if (mod(istep, max(stat_sample_intvl, 1)) /= 0) return
    nsmp = nsmp + 1;  sampled_any = .true.
    if (t_first < 0.0_dp) t_first = t
    t_last = t
    Id = reshape([1.0_dp,0.0_dp,0.0_dp, 0.0_dp,1.0_dp,0.0_dp, 0.0_dp,0.0_dp,1.0_dp], [3,3])
    ! ---- type 1: nodes of the owned domains ----
    if (stat_type(1) == 1) then
       do ib = 1, nown
          d = own_list(ib);  b = dom(d)%b
          call set_block(d)
          R = Id
          if (stat_type1_frame == 1 .and. blk_obj(b) > 0 .and. sixdof_on_flag()) then
             call sixdof_node_pose(blk_obj(b), Rg, O0, dg, vG, omG, Rr, eta, nu);  R = Rg
          end if
          do k = dom(d)%r1(3), dom(d)%r2(3) + 1
           do j = dom(d)%r1(2), dom(d)%r2(2) + 1
            do i = dom(d)%r1(1), dom(d)%r2(1) + 1
               in1 = i - dom(d)%r1(1) + 1;  in2 = j - dom(d)%r1(2) + 1;  in3 = k - dom(d)%r1(3) + 1
               call node_valid(i-1, j-1, k-1, ok1, ok2)
               if (.not. ok1) cycle
               call node_sample(i-1, j-1, k-1, v)
               if (stat_type1_frame == 1) call rotate_sample(v, R, Id)
               call accumulate(v, da(d)%a(:, in1, in2, in3), da(d)%n1(in1, in2, in3), da(d)%n2(in1, in2, in3), .false.)
            end do
           end do
          end do
          ! group 2 at the cells (update 23): the compact gradient of the cell itself, no averaging
          if (stat_group(2) == 1) then
             do k = dom(d)%r1(3), dom(d)%r2(3)
              do j = dom(d)%r1(2), dom(d)%r2(2)
               do i = dom(d)%r1(1), dom(d)%r2(1)
                  in1 = i - dom(d)%r1(1) + 1;  in2 = j - dom(d)%r1(2) + 1;  in3 = k - dom(d)%r1(3) + 1
                  if (cstat(i, j, k) < 0) cycle
                  call cell_grad_sample(i, j, k, v)
                  if (stat_type1_frame == 1) call rotate_sample(v, R, Id)
                  call accumulate2(v, da(d)%c(:, in1, in2, in3), da(d)%nc(in1, in2, in3))
               end do
              end do
             end do
          end if
       end do
    end if
    ! ---- lattices ----
    do L = 1, nlt
       R = Id
       if (lt(L)%node > 0) then
          call sixdof_node_pose(lt(L)%node, Rg, O0, dg, vG, omG, Rr, eta, nu);  R = Rg
       else if (bgr_on) then
          R = bgr_Rg;  O0 = bgr_O0;  dg = bgr_dg
       end if
       do q = 1, lt(L)%npt
          if (lt(L)%node > 0 .or. bgr_on) then
             x = matmul(R, lt(L)%x0(:,q) - O0) + O0 + dg
          else
             x = lt(L)%x0(:,q)
          end if
          call pt_locate(x, istep, lt(L)%lb(q), lt(L)%li(q), lt(L)%lj(q), lt(L)%lk(q), ncand, cb, ci, cj, ck, cr)
          call pt_owner_key(ncand, cb, ci, cj, ck, key, dwin, qwin, st)
          lt(L)%key(q) = key
          if (key > 0) then
             lt(L)%lb(q) = cb(qwin);  lt(L)%li(q) = ci(qwin);  lt(L)%lj(q) = cj(qwin);  lt(L)%lk(q) = ck(qwin)
             call pt_sample(dwin, ci(qwin), cj(qwin), ck(qwin), cr(:,qwin), v)
             call rotate_sample(v, R, lt(L)%pb(:,:,q))
             call stash(L, q, v)
          end if
       end do
       call MPI_Allreduce(lt(L)%key, lt(L)%keymax, lt(L)%npt, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
       do q = 1, lt(L)%npt
          if (lt(L)%keymax(q) > 0 .and. lt(L)%key(q) == lt(L)%keymax(q)) then
             call accumulate(stash_v(:, q), lt(L)%a(:, q), lt(L)%n1(q), lt(L)%n2(q), .true.)
          end if
       end do
       ! the last cell of a point found by another rank is unknown here: keep the local guess
    end do
    ierr = 0
  end subroutine stat_sample

  ! the sampled values of the lattice points this rank evaluated (kept until the owner is known)
  subroutine stash(L, q, v)
    integer, intent(in) :: L, q
    real(dp), intent(in) :: v(NSAMP)
    if (.not. allocated(stash_v)) allocate(stash_v(NSAMP, maxval(lt(1:nlt)%npt)))
    if (size(stash_v, 2) < lt(L)%npt) then
       deallocate(stash_v);  allocate(stash_v(NSAMP, lt(L)%npt))
    end if
    stash_v(:, q) = v
  end subroutine stash

  ! validity of node (im+1, jm+1, km+1): ok1 = a non-hole cell around it, ok2 = a cell with a full gradient stencil
  subroutine node_valid(im, jm, km, ok1, ok2)
    integer, intent(in) :: im, jm, km
    logical, intent(out) :: ok1, ok2
    integer :: a, c, e
    ok1 = .false.;  ok2 = .false.
    do e = 0, 1
     do c = 0, 1
      do a = 0, 1
         if (cstat(im+a, jm+c, km+e) < 0) cycle
         ok1 = .true.
         if (im+a-1 < is-ng .or. im+a+1 > ie+ng .or. jm+c-1 < js-ng .or. jm+c+1 > je+ng .or. &
             km+e-1 < ks-ng .or. km+e+1 > ke+ng) cycle
         ok2 = .true.
      end do
     end do
    end do
  end subroutine node_valid

  !=====================================================================
  ! output of the running statistics (rank 0 writes)
  !=====================================================================
  subroutine stat_write(istep, final)
    integer, intent(in) :: istep
    logical, intent(in) :: final
    if (.not. stat_on .or. .not. sampled_any) return
    if (.not. final) then
       if (stat_output_intvl <= 0) return
       if (mod(istep, stat_output_intvl) /= 0) return
    end if
    if (istep == last_written) return
    last_written = istep
    if (myrank == 0) call make_dir(sdir)
    if (stat_type(1) == 1) call write_type1(istep)
    call write_lattices(istep)
  end subroutine stat_write

  ! ---- the variable lists of the two groups (names with units) ----
  subroutine group_names(g, polar, names, nv)
    integer, intent(in) :: g, polar
    character(len=40), intent(out) :: names(80)
    integer, intent(out) :: nv
    integer :: q
    character(len=8) :: cn(3), gx(3)
    cn = ['u', 'v', 'w'];  gx = ['x', 'y', 'z']
    if (polar == 1) then
       cn = ['u_axial', 'u_rad  ', 'u_tan  '];  gx = ['a', 'r', 't']
    end if
    nv = 0
    if (g == 1) then
       call addn('Coverage [-]')
       do q = 1, nf1
          call addn(trim(fname(q))//'_mean ['//trim(f1unit(q))//']')
       end do
       do q = 1, nf1
          call addn(trim(fname(q))//'_rms ['//trim(f1unit(q))//']')
       end do
       call addn('Rss_'//trim(cn(1))//trim(cn(1))//' [m2/s2]');  call addn('Rss_'//trim(cn(2))//trim(cn(2))//' [m2/s2]')
       call addn('Rss_'//trim(cn(3))//trim(cn(3))//' [m2/s2]');  call addn('Rss_'//trim(cn(1))//trim(cn(2))//' [m2/s2]')
       call addn('Rss_'//trim(cn(2))//trim(cn(3))//' [m2/s2]');  call addn('Rss_'//trim(cn(3))//trim(cn(1))//' [m2/s2]')
       call addn('k_res [m2/s2]')
    else
       call addn('Coverage_grad [-]')
       do q = 1, 9
          call addn('d'//trim(cn(mod(q-1,3)+1))//'_d'//trim(gx((q-1)/3+1))//'_mean [1/s]')
       end do
       do q = 1, 9
          call addn('d'//trim(cn(mod(q-1,3)+1))//'_d'//trim(gx((q-1)/3+1))//'_rms [1/s]')
       end do
       call addn('<d'//trim(cn(1))//'_d'//trim(gx(2))//'*d'//trim(cn(2))//'_d'//trim(gx(1))//'> [1/s2]')
       call addn('<d'//trim(cn(2))//'_d'//trim(gx(3))//'*d'//trim(cn(3))//'_d'//trim(gx(2))//'> [1/s2]')
       call addn('<d'//trim(cn(3))//'_d'//trim(gx(1))//'*d'//trim(cn(1))//'_d'//trim(gx(3))//'> [1/s2]')
       call addn('eps_res [m2/s3]')
       if (turb_on) then
          call addn('eps_model [m2/s3]');  call addn('eps_total [m2/s3]');  call addn('k_model [m2/s2]');  call addn('k_total [m2/s2]')
       end if
       call addn('eta_Kolmogorov [m]');  call addn('lambda_Taylor [m]');  call addn('Re_lambda [-]');  call addn('Delta_over_eta [-]')
       call addn('vort_'//trim(gx(1))//'_mean [1/s]');  call addn('vort_'//trim(gx(2))//'_mean [1/s]');  call addn('vort_'//trim(gx(3))//'_mean [1/s]')
    end if
  contains
    subroutine addn(s)
      character(len=*), intent(in) :: s
      nv = nv + 1;  names(nv) = s
    end subroutine addn
    function fname(q) result(s)
      integer, intent(in) :: q
      character(len=16) :: s
      s = f1name(q)
      if (polar == 1 .and. q <= 3) s = cn(q)
    end function fname
  end subroutine group_names

  ! ---- statistics of one point from its accumulators ----
  subroutine point_stats(a, n1, n2, g, out, nv)
    real(dp), intent(in) :: a(nacc), n1, n2
    integer, intent(in) :: g
    real(dp), intent(out) :: out(80)
    integer, intent(out) :: nv
    integer :: q, c, l
    real(dp) :: m(12), r(12), uu, vv, ww, uv, vw, wu, Gm(3,3), Gq(3,3), X(3), eres, emod, etot, kres, kmod, ktot, nu, dx, &
                eta, lam, rel, s
    out = 0.0_dp;  nv = 0
    if (g == 1) then
       nv = 1 + 2*nf1 + 7
       if (n1 <= 0.0_dp) return
       out(1) = n1/real(max(nsmp, 1), dp)
       do q = 1, nf1
          m(q) = a(iS1 + q - 1)/n1
          r(q) = sqrt(max(a(iQ1 + q - 1)/n1 - m(q)**2, 0.0_dp))
          out(1 + q) = m(q);  out(1 + nf1 + q) = r(q)
       end do
       uu = r(1)**2;  vv = r(2)**2;  ww = r(3)**2
       uv = a(iC1)/n1 - m(1)*m(2);  vw = a(iC1+1)/n1 - m(2)*m(3);  wu = a(iC1+2)/n1 - m(3)*m(1)
       out(2 + 2*nf1 : 7 + 2*nf1) = [uu, vv, ww, uv, vw, wu]
       out(8 + 2*nf1) = 0.5_dp*(uu + vv + ww)
    else
       nv = 1 + 9 + 9 + 3 + 1 + merge(4, 0, turb_on) + 4 + 3
       if (n2 <= 0.0_dp .or. n1 <= 0.0_dp) return
       out(1) = n2/real(max(nsmp, 1), dp)
       do l = 1, 3
          do c = 1, 3
             q = (c-1) + 3*(l-1)
             Gm(c,l) = a(iS2 + q)/n2
             Gq(c,l) = max(a(iQ2 + q)/n2 - Gm(c,l)**2, 0.0_dp)        ! <g'^2>
             out(1 + 1 + q) = Gm(c,l);  out(1 + 9 + 1 + q) = sqrt(Gq(c,l))
          end do
       end do
       X(1) = a(iX2)/n2   - Gm(1,2)*Gm(2,1)
       X(2) = a(iX2+1)/n2 - Gm(2,3)*Gm(3,2)
       X(3) = a(iX2+2)/n2 - Gm(3,1)*Gm(1,3)
       out(20:22) = X
       nu = a(iNU)/n2;  dx = a(iDX)/n2
       eres = nu*(2.0_dp*(Gq(1,1) + Gq(2,2) + Gq(3,3)) + Gq(1,2) + Gq(1,3) + Gq(2,1) + Gq(2,3) + Gq(3,1) + Gq(3,2) &
                  + 2.0_dp*(X(1) + X(2) + X(3)))
       eres = max(eres, 0.0_dp)
       kres = 0.5_dp*(max(a(iQ1)/n1 - (a(iS1)/n1)**2, 0.0_dp) + max(a(iQ1+1)/n1 - (a(iS1+1)/n1)**2, 0.0_dp) &
                      + max(a(iQ1+2)/n1 - (a(iS1+2)/n1)**2, 0.0_dp))
       q = 23;  out(q) = eres
       emod = 0.0_dp;  kmod = 0.0_dp
       if (turb_on) then
          emod = max(a(iEM)/n2, 0.0_dp);  kmod = max(a(iKM)/n2, 0.0_dp)
          out(q+1) = emod;  out(q+2) = eres + emod;  out(q+3) = kmod;  out(q+4) = kres + kmod
          q = q + 4
       end if
       etot = eres + emod;  ktot = kres + kmod
       eta = 0.0_dp;  lam = 0.0_dp;  rel = 0.0_dp;  s = 0.0_dp
       if (etot > 0.0_dp .and. nu > 0.0_dp) then
          eta = (nu**3/etot)**0.25_dp
          lam = sqrt(10.0_dp*nu*ktot/etot)
          rel = sqrt(2.0_dp*ktot/3.0_dp)*lam/nu
          s = dx/eta
       end if
       out(q+1) = eta;  out(q+2) = lam;  out(q+3) = rel;  out(q+4) = s
       out(q+5) = Gm(3,2) - Gm(2,3);  out(q+6) = Gm(1,3) - Gm(3,1);  out(q+7) = Gm(2,1) - Gm(1,2)
    end if
  end subroutine point_stats

  ! ---- type 1 output ----
  subroutine write_type1(istep)
    integer, intent(in) :: istep
    type :: nb_t
       real(dp), allocatable :: a(:,:,:,:), n1(:,:,:), n2(:,:,:), c(:,:,:,:), nc(:,:,:)
    end type nb_t
    type(nb_t), allocatable :: nb(:)
    type(nbf_t), allocatable :: fbs(:)
    real(dp), allocatable :: piece(:,:,:,:), pn(:,:,:,:), fb(:,:,:,:), co(:,:,:,:)
    integer :: b, d, ib, i, j, k, ierr, n1(3), n2(3), ncnt, g, nv, c, gq, iu, l, a2, c2, e2, nn
    integer :: stat(MPI_STATUS_SIZE)
    real(dp) :: out(80), avg(80)
    logical :: isq(80)
    character(len=40) :: names(80)
    character(len=256) :: fname
    if (myrank == 0) then
       if (.not. groups_built) call build_node_groups
       allocate(nb(nblk))
       do b = 1, nblk
          allocate(nb(b)%a(nacc, blk(b)%ni, blk(b)%nj, blk(b)%nk), nb(b)%n1(blk(b)%ni, blk(b)%nj, blk(b)%nk), &
                   nb(b)%n2(blk(b)%ni, blk(b)%nj, blk(b)%nk))
          allocate(nb(b)%c(na2c, blk(b)%nx, blk(b)%ny, blk(b)%nz), nb(b)%nc(blk(b)%nx, blk(b)%ny, blk(b)%nz))
          nb(b)%a = 0.0_dp;  nb(b)%n1 = 0.0_dp;  nb(b)%n2 = 0.0_dp;  nb(b)%c = 0.0_dp;  nb(b)%nc = 0.0_dp
          do d = blk(b)%dom0, blk(b)%dom0 + blk(b)%ncpu - 1
             n1 = dom(d)%r1;  n2 = dom(d)%r2 + 1
             allocate(piece(nacc, n1(1):n2(1), n1(2):n2(2), n1(3):n2(3)), pn(2, n1(1):n2(1), n1(2):n2(2), n1(3):n2(3)))
             allocate(co(na2c + 1, n1(1):n2(1)-1, n1(2):n2(2)-1, n1(3):n2(3)-1))
             if (dom(d)%owner == 0) then
                piece = da(d)%a;  pn(1,:,:,:) = da(d)%n1;  pn(2,:,:,:) = da(d)%n2
                co(1:na2c,:,:,:) = da(d)%c;  co(na2c+1,:,:,:) = da(d)%nc
             else
                ncnt = (n2(1)-n1(1)+1)*(n2(2)-n1(2)+1)*(n2(3)-n1(3)+1)
                call MPI_Recv(piece, ncnt*nacc, MPI_DOUBLE_PRECISION, dom(d)%owner, 2300+d, MPI_COMM_WORLD, stat, ierr)
                call MPI_Recv(pn, ncnt*2, MPI_DOUBLE_PRECISION, dom(d)%owner, 2300+ndom+d, MPI_COMM_WORLD, stat, ierr)
                ncnt = (n2(1)-n1(1))*(n2(2)-n1(2))*(n2(3)-n1(3))
                call MPI_Recv(co, ncnt*(na2c+1), MPI_DOUBLE_PRECISION, dom(d)%owner, 2300+2*ndom+d, MPI_COMM_WORLD, stat, ierr)
             end if
             ! domain boundary nodes are shared by neighbouring domains (identical sums): take any
             nb(b)%a(:, n1(1):n2(1), n1(2):n2(2), n1(3):n2(3)) = piece
             nb(b)%n1(n1(1):n2(1), n1(2):n2(2), n1(3):n2(3)) = pn(1,:,:,:)
             nb(b)%n2(n1(1):n2(1), n1(2):n2(2), n1(3):n2(3)) = pn(2,:,:,:)
             nb(b)%c(:, n1(1):n2(1)-1, n1(2):n2(2)-1, n1(3):n2(3)-1) = co(1:na2c,:,:,:)
             nb(b)%nc(n1(1):n2(1)-1, n1(2):n2(2)-1, n1(3):n2(3)-1) = co(na2c+1,:,:,:)
             deallocate(piece, pn, co)
          end do
       end do
       ! grid (once): the t = 0 coordinates are in the initial grid file; here the current blk arrays may have
       ! moved, so the statistics grid is written at the first output from the stored copy
       call write_stat_grid
       do g = 1, 2
          if (stat_group(g) == 0) cycle
          call group_names(g, 0, names, nv)
          do b = 1, nblk
             allocate(fb(nv, blk(b)%ni, blk(b)%nj, blk(b)%nk))
             if (g == 1) then
                do k = 1, blk(b)%nk
                 do j = 1, blk(b)%nj
                  do i = 1, blk(b)%ni
                     call point_stats(nb(b)%a(:,i,j,k), nb(b)%n1(i,j,k), nb(b)%n2(i,j,k), g, out, nv)
                     fb(:,i,j,k) = out(1:nv)
                  end do
                 end do
                end do
             else
                ! group 2 (update 23): the statistics of every cell, then the node = the average of its
                ! sampled cells (none: 0)
                allocate(co(nv, blk(b)%nx, blk(b)%ny, blk(b)%nz))
                do k = 1, blk(b)%nz
                 do j = 1, blk(b)%ny
                  do i = 1, blk(b)%nx
                     call cell_stats(nb(b)%c(:,i,j,k), nb(b)%nc(i,j,k), out, nv)
                     co(:,i,j,k) = out(1:nv)
                  end do
                 end do
                end do
                ! the rms columns are averaged as variances (sqrt afterwards), so that eps_res of a node
                ! recomputed from its exported components equals the averaged eps_res of its cells
                do l = 1, nv
                   isq(l) = (index(names(l), '_rms [1/s]') > 0)
                   if (isq(l)) co(l,:,:,:) = co(l,:,:,:)**2
                end do
                do k = 1, blk(b)%nk
                 do j = 1, blk(b)%nj
                  do i = 1, blk(b)%ni
                     avg = 0.0_dp;  nn = 0
                     do e2 = max(k-1, 1), min(k, blk(b)%nz)
                      do c2 = max(j-1, 1), min(j, blk(b)%ny)
                       do a2 = max(i-1, 1), min(i, blk(b)%nx)
                          if (nb(b)%nc(a2, c2, e2) <= 0.0_dp) cycle
                          avg(1:nv) = avg(1:nv) + co(:, a2, c2, e2);  nn = nn + 1
                       end do
                      end do
                     end do
                     if (nn > 0) avg(1:nv) = avg(1:nv)/real(nn, dp)
                     do l = 1, nv
                        if (isq(l)) avg(l) = sqrt(max(avg(l), 0.0_dp))
                     end do
                     fb(:,i,j,k) = avg(1:nv)
                  end do
                 end do
                end do
                deallocate(co)
             end if
             if (b == 1) allocate(fbs(nblk))
             call move_alloc(fb, fbs(b)%f)
          end do
          ! coincident nodes of different blocks: the average of their statistics (the rms columns as
          ! variances, so that the exported components stay consistent with k_res and eps_res)
          do l = 1, nv
             isq(l) = (index(names(l), '_rms') > 0)
          end do
          do gq = 1, ngrp
             avg = 0.0_dp
             do c = gstart(gq), gstart(gq+1)-1
                out(1:nv) = fbs(gblk(c))%f(:, gi(c), gj(c), gk(c))
                do l = 1, nv
                   if (isq(l)) out(l) = out(l)**2
                end do
                avg(1:nv) = avg(1:nv) + out(1:nv)
             end do
             avg = avg/real(gstart(gq+1) - gstart(gq), dp)
             do l = 1, nv
                if (isq(l)) avg(l) = sqrt(max(avg(l), 0.0_dp))
             end do
             do c = gstart(gq), gstart(gq+1)-1
                fbs(gblk(c))%f(:, gi(c), gj(c), gk(c)) = avg(1:nv)
             end do
          end do
          write(fname,'(A,A,A,I1,A,I8.8,A)') sdir, trim(casename), '_stat_grid_g', g, '_', istep, '.p3d'
          open(newunit=iu, file=trim(fname), form='unformatted', status='replace', action='write')
          write(iu) nblk
          write(iu) (blk(b)%ni, blk(b)%nj, blk(b)%nk, nv, b=1,nblk)
          do b = 1, nblk
             write(iu) ((((real(fbs(b)%f(l,i,j,k),4), i=1,blk(b)%ni), j=1,blk(b)%nj), k=1,blk(b)%nk), l=1,nv)
             deallocate(fbs(b)%f)
          end do
          close(iu)
          deallocate(fbs)
          write(fname,'(A,A,A,I1,A)') sdir, trim(casename), '_stat_grid_g', g, '.nam'
          open(newunit=iu, file=trim(fname), status='replace', action='write')
          do l = 1, nv
             write(iu,'(A)') trim(names(l))
          end do
          close(iu)
       end do
       write(*,'(A,I0,A,I0,A)') '  wrote statistics (type 1, grid blocks) of step ', istep, ': ', nsmp, ' samples -> '// &
          sdir//trim(casename)//'_stat_grid_g*_NNNNNNNN.p3d'
       deallocate(nb)
    else
       do ib = 1, nown
          d = own_list(ib)
          n1 = dom(d)%r1;  n2 = dom(d)%r2 + 1
          ncnt = (n2(1)-n1(1)+1)*(n2(2)-n1(2)+1)*(n2(3)-n1(3)+1)
          allocate(pn(2, n1(1):n2(1), n1(2):n2(2), n1(3):n2(3)))
          pn(1,:,:,:) = da(d)%n1;  pn(2,:,:,:) = da(d)%n2
          call MPI_Send(da(d)%a, ncnt*nacc, MPI_DOUBLE_PRECISION, 0, 2300+d, MPI_COMM_WORLD, ierr)
          call MPI_Send(pn, ncnt*2, MPI_DOUBLE_PRECISION, 0, 2300+ndom+d, MPI_COMM_WORLD, ierr)
          deallocate(pn)
          allocate(co(na2c + 1, n1(1):n2(1)-1, n1(2):n2(2)-1, n1(3):n2(3)-1))
          co(1:na2c,:,:,:) = da(d)%c;  co(na2c+1,:,:,:) = da(d)%nc
          ncnt = (n2(1)-n1(1))*(n2(2)-n1(2))*(n2(3)-n1(3))
          call MPI_Send(co, ncnt*(na2c+1), MPI_DOUBLE_PRECISION, 0, 2300+2*ndom+d, MPI_COMM_WORLD, ierr)
          deallocate(co)
       end do
    end if
  end subroutine write_type1

  ! the t = 0 grid of the statistics (stored at init, written once)
  subroutine write_stat_grid
    integer :: iu, b, i, j, k
    if (grid_written) return
    open(newunit=iu, file=sdir//trim(casename)//'_stat_grid.grd', form='unformatted', status='replace', action='write')
    write(iu) nblk
    write(iu) (blk(b)%ni, blk(b)%nj, blk(b)%nk, b=1,nblk)
    do b = 1, nblk
       write(iu) (((real(g0(b)%x(1,i,j,k),4), i=1,blk(b)%ni), j=1,blk(b)%nj), k=1,blk(b)%nk), &
                 (((real(g0(b)%x(2,i,j,k),4), i=1,blk(b)%ni), j=1,blk(b)%nj), k=1,blk(b)%nk), &
                 (((real(g0(b)%x(3,i,j,k),4), i=1,blk(b)%ni), j=1,blk(b)%nj), k=1,blk(b)%nk)
    end do
    close(iu)
    grid_written = .true.
    write(*,'(A)') '  wrote the statistics grid (t = 0): '//sdir//trim(casename)//'_stat_grid.grd'
  end subroutine write_stat_grid

  ! ---- lattice output ----
  subroutine write_lattices(istep)
    integer, intent(in) :: istep
    integer :: L, q, g, nv, iu, ierr, l2, i, j, k
    real(dp), allocatable :: asum(:,:), n1s(:), n2s(:), fb(:,:)
    real(dp) :: out(80)
    character(len=40) :: names(80)
    character(len=256) :: fname
    do L = 1, nlt
       allocate(asum(nacc, lt(L)%npt), n1s(lt(L)%npt), n2s(lt(L)%npt))
       call MPI_Reduce(lt(L)%a, asum, nacc*lt(L)%npt, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
       call MPI_Reduce(lt(L)%n1, n1s, lt(L)%npt, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
       call MPI_Reduce(lt(L)%n2, n2s, lt(L)%npt, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
       if (myrank == 0) then
          ! grid of the lattice (t = 0 positions in the grid frame), once
          write(fname,'(A,A,A,A,A)') sdir, trim(casename), '_stat_', trim(lt(L)%stem), '.grd'
          if (.not. lt_grid_written(L)) then
             open(newunit=iu, file=trim(fname), form='unformatted', status='replace', action='write')
             write(iu) 1
             write(iu) lt(L)%n(1), lt(L)%n(2), lt(L)%n(3)
             write(iu) (real(lt(L)%x0(1,q),4), q=1,lt(L)%npt), (real(lt(L)%x0(2,q),4), q=1,lt(L)%npt), &
                       (real(lt(L)%x0(3,q),4), q=1,lt(L)%npt)
             close(iu)
             lt_grid_written(L) = .true.
          end if
          do g = 1, 2
             if (stat_group(g) == 0) cycle
             call group_names(g, lt(L)%polar, names, nv)
             allocate(fb(nv, lt(L)%npt))
             do q = 1, lt(L)%npt
                call point_stats(asum(:,q), n1s(q), n2s(q), g, out, nv)
                fb(:,q) = out(1:nv)
             end do
             write(fname,'(A,A,A,A,A,I1,A,I8.8,A)') sdir, trim(casename), '_stat_', trim(lt(L)%stem), '_g', g, '_', istep, '.p3d'
             open(newunit=iu, file=trim(fname), form='unformatted', status='replace', action='write')
             write(iu) 1
             write(iu) lt(L)%n(1), lt(L)%n(2), lt(L)%n(3), nv
             write(iu) ((real(fb(l2,q),4), q=1,lt(L)%npt), l2=1,nv)
             close(iu)
             deallocate(fb)
             write(fname,'(A,A,A,A,A,I1,A)') sdir, trim(casename), '_stat_', trim(lt(L)%stem), '_g', g, '.nam'
             open(newunit=iu, file=trim(fname), status='replace', action='write')
             do l2 = 1, nv
                write(iu,'(A)') trim(names(l2))
             end do
             close(iu)
          end do
          write(*,'(A,I0,A,I0,A,F6.3)') '  wrote statistics of '//trim(lt(L)%stem)//' of step ', istep, ': ', nsmp, &
             ' samples, mean coverage ', sum(n1s)/real(max(nsmp, 1), dp)/real(lt(L)%npt, dp)
       end if
       deallocate(asum, n1s, n2s)
    end do
    i = 0;  j = 0;  k = 0
  end subroutine write_lattices

  !=====================================================================
  ! restart: the accumulators after the PID states (rank 0 writes; the lattice sums reduced)
  !=====================================================================
  subroutine stat_restart_write(iu)
    integer, intent(in) :: iu
    integer :: ib, d, ierr, L, n1(3), n2(3), ncnt
    real(dp), allocatable :: piece(:,:,:,:), pn(:,:,:,:), asum(:,:), n1s(:), n2s(:), co(:,:,:,:)
    integer :: stat(MPI_STATUS_SIZE)
    if (.not. stat_on) then
       if (myrank == 0) write(iu) 0
       return
    end if
    if (myrank == 0) write(iu) 1, nacc + 1000*na2c, nsmp, t_first, t_last
    if (stat_type(1) == 1) then
       do d = 1, ndom
          n1 = dom(d)%r1;  n2 = dom(d)%r2 + 1
          ncnt = (n2(1)-n1(1)+1)*(n2(2)-n1(2)+1)*(n2(3)-n1(3)+1)
          if (myrank == 0) then
             allocate(piece(nacc, n1(1):n2(1), n1(2):n2(2), n1(3):n2(3)), pn(2, n1(1):n2(1), n1(2):n2(2), n1(3):n2(3)))
             if (dom(d)%owner == 0) then
                piece = da(d)%a;  pn(1,:,:,:) = da(d)%n1;  pn(2,:,:,:) = da(d)%n2
             else
                call MPI_Recv(piece, ncnt*nacc, MPI_DOUBLE_PRECISION, dom(d)%owner, 2400+d, MPI_COMM_WORLD, stat, ierr)
                call MPI_Recv(pn, ncnt*2, MPI_DOUBLE_PRECISION, dom(d)%owner, 2400+ndom+d, MPI_COMM_WORLD, stat, ierr)
             end if
             write(iu) piece;  write(iu) pn
             deallocate(piece, pn)
             ! the cell accumulators of group 2 (update 23)
             ncnt = (n2(1)-n1(1))*(n2(2)-n1(2))*(n2(3)-n1(3))
             allocate(co(na2c+1, n1(1):n2(1)-1, n1(2):n2(2)-1, n1(3):n2(3)-1))
             if (dom(d)%owner == 0) then
                co(1:na2c,:,:,:) = da(d)%c;  co(na2c+1,:,:,:) = da(d)%nc
             else
                call MPI_Recv(co, ncnt*(na2c+1), MPI_DOUBLE_PRECISION, dom(d)%owner, 2400+2*ndom+d, MPI_COMM_WORLD, stat, ierr)
             end if
             write(iu) co
             deallocate(co)
          else if (dom(d)%owner == myrank) then
             allocate(pn(2, n1(1):n2(1), n1(2):n2(2), n1(3):n2(3)))
             pn(1,:,:,:) = da(d)%n1;  pn(2,:,:,:) = da(d)%n2
             call MPI_Send(da(d)%a, ncnt*nacc, MPI_DOUBLE_PRECISION, 0, 2400+d, MPI_COMM_WORLD, ierr)
             call MPI_Send(pn, ncnt*2, MPI_DOUBLE_PRECISION, 0, 2400+ndom+d, MPI_COMM_WORLD, ierr)
             deallocate(pn)
             ncnt = (n2(1)-n1(1))*(n2(2)-n1(2))*(n2(3)-n1(3))
             allocate(co(na2c+1, n1(1):n2(1)-1, n1(2):n2(2)-1, n1(3):n2(3)-1))
             co(1:na2c,:,:,:) = da(d)%c;  co(na2c+1,:,:,:) = da(d)%nc
             call MPI_Send(co, ncnt*(na2c+1), MPI_DOUBLE_PRECISION, 0, 2400+2*ndom+d, MPI_COMM_WORLD, ierr)
             deallocate(co)
          end if
       end do
    end if
    do L = 1, nlt
       allocate(asum(nacc, lt(L)%npt), n1s(lt(L)%npt), n2s(lt(L)%npt))
       call MPI_Reduce(lt(L)%a, asum, nacc*lt(L)%npt, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
       call MPI_Reduce(lt(L)%n1, n1s, lt(L)%npt, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
       call MPI_Reduce(lt(L)%n2, n2s, lt(L)%npt, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
       if (myrank == 0) then
          write(iu) lt(L)%npt;  write(iu) asum;  write(iu) n1s;  write(iu) n2s
       end if
       deallocate(asum, n1s, n2s)
    end do
    ib = 0
  end subroutine stat_restart_write

  subroutine stat_restart_read(iu, ok)
    integer, intent(in) :: iu
    logical, intent(out) :: ok
    integer :: flag, nacc_f, ns_f, d, L, ierr, npt_f, n1(3), n2(3), ncnt, ios
    real(dp) :: tf, tl
    real(dp), allocatable :: piece(:,:,:,:), pn(:,:,:,:), asum(:,:), n1s(:), n2s(:), co(:,:,:,:)
    integer :: stat(MPI_STATUS_SIZE)
    ok = .false.
    flag = 0
    if (myrank == 0) then
       read(iu, iostat=ios) flag
       if (ios /= 0) flag = 0
       if (flag == 1) then
          backspace(iu)
          read(iu) flag, nacc_f, ns_f, tf, tl
          if (nacc_f /= nacc + 1000*na2c .or. .not. stat_on) then
             if (stat_on) write(*,'(A)') '  NOTE (statistics): the restart file holds statistics of another set-up - starting fresh'
             flag = -1
          end if
       end if
    end if
    call MPI_Bcast(flag, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if (flag /= 1) return
    call MPI_Bcast(ns_f, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(tf, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(tl, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    nsmp = ns_f;  t_first = tf;  t_last = tl;  sampled_any = (nsmp > 0)
    if (stat_type(1) == 1) then
       do d = 1, ndom
          n1 = dom(d)%r1;  n2 = dom(d)%r2 + 1
          ncnt = (n2(1)-n1(1)+1)*(n2(2)-n1(2)+1)*(n2(3)-n1(3)+1)
          if (myrank == 0) then
             allocate(piece(nacc, n1(1):n2(1), n1(2):n2(2), n1(3):n2(3)), pn(2, n1(1):n2(1), n1(2):n2(2), n1(3):n2(3)))
             read(iu) piece;  read(iu) pn
             if (dom(d)%owner == 0) then
                da(d)%a = piece;  da(d)%n1 = pn(1,:,:,:);  da(d)%n2 = pn(2,:,:,:)
             else
                call MPI_Send(piece, ncnt*nacc, MPI_DOUBLE_PRECISION, dom(d)%owner, 2500+d, MPI_COMM_WORLD, ierr)
                call MPI_Send(pn, ncnt*2, MPI_DOUBLE_PRECISION, dom(d)%owner, 2500+ndom+d, MPI_COMM_WORLD, ierr)
             end if
             deallocate(piece, pn)
             ncnt = (n2(1)-n1(1))*(n2(2)-n1(2))*(n2(3)-n1(3))
             allocate(co(na2c+1, n1(1):n2(1)-1, n1(2):n2(2)-1, n1(3):n2(3)-1))
             read(iu) co
             if (dom(d)%owner == 0) then
                da(d)%c = co(1:na2c,:,:,:);  da(d)%nc = co(na2c+1,:,:,:)
             else
                call MPI_Send(co, ncnt*(na2c+1), MPI_DOUBLE_PRECISION, dom(d)%owner, 2500+2*ndom+d, MPI_COMM_WORLD, ierr)
             end if
             deallocate(co)
          else if (dom(d)%owner == myrank) then
             allocate(pn(2, n1(1):n2(1), n1(2):n2(2), n1(3):n2(3)))
             call MPI_Recv(da(d)%a, ncnt*nacc, MPI_DOUBLE_PRECISION, 0, 2500+d, MPI_COMM_WORLD, stat, ierr)
             call MPI_Recv(pn, ncnt*2, MPI_DOUBLE_PRECISION, 0, 2500+ndom+d, MPI_COMM_WORLD, stat, ierr)
             da(d)%n1 = pn(1,:,:,:);  da(d)%n2 = pn(2,:,:,:)
             deallocate(pn)
             ncnt = (n2(1)-n1(1))*(n2(2)-n1(2))*(n2(3)-n1(3))
             allocate(co(na2c+1, n1(1):n2(1)-1, n1(2):n2(2)-1, n1(3):n2(3)-1))
             call MPI_Recv(co, ncnt*(na2c+1), MPI_DOUBLE_PRECISION, 0, 2500+2*ndom+d, MPI_COMM_WORLD, stat, ierr)
             da(d)%c = co(1:na2c,:,:,:);  da(d)%nc = co(na2c+1,:,:,:)
             deallocate(co)
          end if
       end do
    end if
    do L = 1, nlt
       if (myrank == 0) then
          read(iu) npt_f
          if (npt_f /= lt(L)%npt) then
             write(*,'(A)') '  NOTE (statistics): the restart file holds a lattice of another size - its sums are dropped'
             allocate(asum(nacc, npt_f), n1s(npt_f), n2s(npt_f))
             read(iu) asum;  read(iu) n1s;  read(iu) n2s
             deallocate(asum, n1s, n2s)
          else
             allocate(asum(nacc, lt(L)%npt), n1s(lt(L)%npt), n2s(lt(L)%npt))
             read(iu) asum;  read(iu) n1s;  read(iu) n2s
             lt(L)%a = asum;  lt(L)%n1 = n1s;  lt(L)%n2 = n2s          ! rank 0 keeps the restored sums
             deallocate(asum, n1s, n2s)
          end if
       end if
    end do
    ok = .true.
    if (myrank == 0) write(*,'(A,I0,A)') '  statistics restored from the restart file: ', nsmp, ' samples'
  end subroutine stat_restart_read

  subroutine stat_summary
    if (.not. stat_on .or. myrank /= 0) return
    write(*,'(A,I0,A,F10.4,A,F10.4,A)') '  statistics: ', nsmp, ' samples between t = ', t_first, ' and ', t_last, ' s'
  end subroutine stat_summary

  subroutine stat_finalize
    integer :: d
    if (.not. stat_on) return
    if (allocated(da)) then
       do d = 1, size(da)
          if (allocated(da(d)%a)) deallocate(da(d)%a, da(d)%n1, da(d)%n2, da(d)%c, da(d)%nc)
       end do
       deallocate(da)
    end if
  end subroutine stat_finalize
end module statistics
