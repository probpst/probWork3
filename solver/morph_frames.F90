!=======================================================================
! morph_frames.F90  -  PRESCRIBED GRID DEFORMATION FROM FRAME FILES
!                      (heart CFD extension, update H1)
!
!   A flexible part (SW_FSI = 1) with  FSI_Frames = "<file>"  has no
!   structure: the node coordinates of ALL its blocks follow a sequence
!   of pre-computed deformation frames,  x(t) = X0 + scale * d(t),  d(t)
!   interpolated in time between the frames (periodic cubic Hermite by
!   default).  The rest of the moving-grid machinery of update 24 is
!   reused unchanged: metrics rebuilt every step, discrete grid velocity
!   (x^(n+1) - x^n) / dt at the half lattice -> ALE grid fluxes, wall
!   velocity of the no-slip patches, flux re-projection, GCL term.
!
!   Frame file (little-endian stream, written by tools/write_frames_bin.py):
!     int32 magic 0x4652414D ('FRAM'), int32 version (1), int32 nblk,
!     int32 nframes, int32 periodic (0/1), float64 period [s],
!     float64 scale (file length unit -> m), int32 dims(3, nblk),
!     float64 times(nframes) [s], then for every frame, for every block:
!     float32 dx(ni,nj,nk), dy(ni,nj,nk), dz(ni,nj,nk) (i fastest), the
!     displacement of every node from the INITIAL grid (<case>.grd).
!
!   Keys of the part group ($SixDOF):
!     FSI_Frames        = "file"   (turns the mode on; needs SW_FSI = 1, Motion_DOF = 0 0 0 0 0 0)
!     FSI_Frames_Interp = 3        (3 cubic Hermite, 1 linear)
!     FSI_Frames_Period = 0        (> 0 overrides the file's period; 0 = file)
!     FSI_Frames_Scale  = 0        (> 0 overrides the file's length scale)
!     FSI_Frames_Time0  = 0        (time offset: the grid at solver time t is the frame state at t + Time0)
!   Mode = 2 writes the deformed grids of the frame motion (grid-movement check).
!=======================================================================
module morph_frames
  use mpi
  use global_vars
  use mpi_vars
  use sixdof, only: sixdof_nbody, sixdof_on_flag, sixdof_fsi_def, fsidef_t, sixdof_part_name
  implicit none
  private
  public :: frames_hooks, frames_present, frames_part, frames_blk, frames_morph, frames_hvel, frames_commit, &
            frames_check_nodes, frames_info, frames_nparts, frames_time

  integer, parameter :: NSLOT = 4

  type :: fblk_t
     logical :: on = .false.
     integer :: part = 0
     real(dp), allocatable :: X0(:,:,:,:)          ! (3, ni, nj, nk) initial coordinates [m]
     real(dp), allocatable :: xprev(:,:,:,:)       ! (3, ni, nj, nk) coordinates at the end of the previous step [m]
     real(dp), allocatable :: vh(:,:,:,:)          ! (3, 2ni-1, 2nj-1, 2nk-1) discrete grid velocity, non-dim (owned blocks)
     real(4),  allocatable :: d(:,:,:,:,:)         ! (3, ni, nj, nk, NSLOT) frame displacements held in memory (file units)
  end type fblk_t

  type :: fpart_t
     integer :: k = 0                              ! 6-DoF node (part)
     character(len=256) :: file = ' '
     integer :: iu = -1
     integer :: nf = 0, interp = 3
     logical :: periodic = .true.
     real(dp) :: period = 0.0_dp, scale = 1.0_dp, t0 = 0.0_dp
     real(dp), allocatable :: tf(:)                ! (nf) frame times [s]
     integer :: slot_f(NSLOT) = -1                 ! frame index (1-based) held in every slot, -1 = empty
     integer :: next_slot = 1
     integer(8) :: hdr_bytes = 0, frame_bytes = 0
     integer(8), allocatable :: boff(:)            ! (nblk) byte offset of a block inside a frame
     integer :: nblk_on = 0
     real(dp) :: t_last = -1.0e300_dp              ! time of the last evaluation (log)
     integer  :: nread = 0                         ! frames read from the file
  end type fpart_t

  type(fblk_t), allocatable :: fb(:)
  type(fpart_t), allocatable :: fp(:)
  integer :: np = 0
  logical :: inited = .false., own_ready = .false.
  logical, allocatable :: own_blk(:)

contains

  integer function frames_nparts()
    frames_nparts = np
  end function frames_nparts

  logical function frames_present()
    frames_present = (np > 0)
  end function frames_present

  ! the part index of 6-DoF node k (0: not a frame part)
  integer function frames_part(k)
    integer, intent(in) :: k
    integer :: ip
    frames_part = 0
    do ip = 1, np
       if (fp(ip)%k == k) frames_part = ip
    end do
  end function frames_part

  logical function frames_blk(b)
    integer, intent(in) :: b
    frames_blk = .false.
    if (.not. allocated(fb)) return
    frames_blk = fb(b)%on
  end function frames_blk

  ! the time the grid is evaluated at (the target time of the step: Main sets t_now = t^(n+1) before the motion)
  real(dp) function frames_time()
    frames_time = t_now
  end function frames_time

  subroutine frames_abort(msg)
    character(len=*), intent(in) :: msg
    integer :: ierr
    if (myrank == 0) write(*,'(A)') ' ERROR (frames): '//msg
    call MPI_Abort(MPI_COMM_WORLD, 81, ierr)
  end subroutine frames_abort

  !=====================================================================
  ! Called by Main right after sixdof_precheck (the grid is the initial
  ! one, the 6-DoF nodes are known): find the frame parts, read the file
  ! headers, keep the initial coordinates, set the Mode-2 hook.
  !=====================================================================
  subroutine frames_hooks
    integer :: nb, k, ip, b, ni, nj, nk, ierr
    type(fsidef_t) :: def
    if (inited) return
    inited = .true.
    np = 0
    if (.not. sixdof_on_flag()) return
    nb = sixdof_nbody()
    do k = 1, nb
       def = sixdof_fsi_def(k)
       if (len_trim(def%frames) > 0) np = np + 1
    end do
    if (np == 0) return
    allocate(fp(np), fb(nblk))
    ip = 0
    do k = 1, nb
       def = sixdof_fsi_def(k)
       if (len_trim(def%frames) == 0) cycle
       ip = ip + 1
       fp(ip)%k = k;  fp(ip)%file = def%frames;  fp(ip)%interp = def%fr_interp;  fp(ip)%t0 = def%fr_t0
       if (def%sw /= 1) call frames_abort(sixdof_part_name(k)//': FSI_Frames needs SW_FSI = 1')
       if (fp(ip)%interp /= 1 .and. fp(ip)%interp /= 3) call frames_abort('FSI_Frames_Interp must be 1 (linear) or 3 (cubic)')
       call open_frames(ip, def)
       do b = 1, nblk
          if (blk_obj(b) /= k) cycle
          fb(b)%on = .true.;  fb(b)%part = ip
          ni = blk(b)%ni;  nj = blk(b)%nj;  nk = blk(b)%nk
          allocate(fb(b)%X0(3, ni, nj, nk), fb(b)%xprev(3, ni, nj, nk), fb(b)%d(3, ni, nj, nk, NSLOT))
          fb(b)%X0(1,:,:,:) = blk(b)%xn;  fb(b)%X0(2,:,:,:) = blk(b)%yn;  fb(b)%X0(3,:,:,:) = blk(b)%zn
          fb(b)%xprev = fb(b)%X0
          fb(b)%d = 0.0_4
          fp(ip)%nblk_on = fp(ip)%nblk_on + 1
       end do
       if (fp(ip)%nblk_on == 0) call frames_abort(sixdof_part_name(k)//': FSI_Frames but no block of the part (rows 6-8 of the csv)')
    end do
    frames_nodes_fn => frames_check_nodes
    frames_on = .true.
    ierr = 0
  end subroutine frames_hooks

  ! the file header of part ip (every rank reads the file; the blocks are replicated)
  subroutine open_frames(ip, def)
    integer, intent(in) :: ip
    type(fsidef_t), intent(in) :: def
    integer :: ios, magic, ver, nbf, nf, per, b
    integer, allocatable :: dims(:,:)
    real(dp) :: period, scale
    integer(8) :: pos
    character(len=16) :: s1, s2
    open(newunit=fp(ip)%iu, file=trim(fp(ip)%file), access='stream', form='unformatted', status='old', action='read', iostat=ios)
    if (ios /= 0) call frames_abort('cannot open the frame file '//trim(fp(ip)%file))
    read(fp(ip)%iu, iostat=ios) magic, ver, nbf, nf, per, period, scale
    if (ios /= 0 .or. magic /= 1179210061) call frames_abort(trim(fp(ip)%file)//' is not a probCFD frame file (magic FRAM)')
    if (ver /= 1) call frames_abort(trim(fp(ip)%file)//': unknown version')
    if (nbf /= nblk) then
       write(s1,'(I0)') nbf;  write(s2,'(I0)') nblk
       call frames_abort(trim(fp(ip)%file)//' holds '//trim(s1)//' blocks, the grid has '//trim(s2))
    end if
    if (nf < 2) call frames_abort(trim(fp(ip)%file)//': at least 2 frames are needed')
    allocate(dims(3, nblk), fp(ip)%tf(nf), fp(ip)%boff(nblk))
    read(fp(ip)%iu) dims
    read(fp(ip)%iu) fp(ip)%tf
    do b = 1, nblk
       if (dims(1,b) /= blk(b)%ni .or. dims(2,b) /= blk(b)%nj .or. dims(3,b) /= blk(b)%nk) then
          write(s1,'(I0)') b
          call frames_abort(trim(fp(ip)%file)//': the dimensions of block '//trim(s1)//' differ from the grid')
       end if
    end do
    fp(ip)%nf = nf;  fp(ip)%periodic = (per == 1);  fp(ip)%period = period;  fp(ip)%scale = scale
    if (def%fr_period > 0.0_dp) fp(ip)%period = def%fr_period
    if (def%fr_scale > 0.0_dp) fp(ip)%scale = def%fr_scale
    do b = 2, nf
       if (fp(ip)%tf(b) <= fp(ip)%tf(b-1)) call frames_abort(trim(fp(ip)%file)//': the frame times must increase')
    end do
    if (fp(ip)%periodic) then
       if (fp(ip)%period <= fp(ip)%tf(nf) - fp(ip)%tf(1)) call frames_abort(trim(fp(ip)%file)//': the period must exceed the frame span')
    end if
    fp(ip)%hdr_bytes = 4*5 + 8*2 + 4*3*nblk + 8*nf
    pos = 0
    do b = 1, nblk
       fp(ip)%boff(b) = pos
       pos = pos + 3_8*4_8*int(blk(b)%ni,8)*int(blk(b)%nj,8)*int(blk(b)%nk,8)
    end do
    fp(ip)%frame_bytes = pos
    ! the file length
    inquire(unit=fp(ip)%iu, size=pos)
    if (pos < fp(ip)%hdr_bytes + fp(ip)%frame_bytes*int(nf,8)) call frames_abort(trim(fp(ip)%file)//' is shorter than its header announces')
    deallocate(dims)
  end subroutine open_frames

  subroutine frames_info(k)
    integer, intent(in) :: k
    integer :: ip
    ip = frames_part(k)
    if (ip == 0 .or. myrank /= 0) return
    write(*,'(A)') '  '//sixdof_part_name(k)//': PRESCRIBED DEFORMATION from frame file '//trim(fp(ip)%file)
    write(*,'(A,I0,A,I0,A,F10.6,A,F10.6,A)') '    frames ', fp(ip)%nf, ' on ', fp(ip)%nblk_on, ' blocks, times ', fp(ip)%tf(1), ' .. ', &
       fp(ip)%tf(fp(ip)%nf), ' s'
    if (fp(ip)%periodic) then
       write(*,'(A,F10.6,A)') '    periodic, period ', fp(ip)%period, ' s'
    else
       write(*,'(A)') '    not periodic (the grid holds the first / last frame outside the frame span)'
    end if
    write(*,'(A,ES10.3,A,A,A,F10.6,A)') '    length scale ', fp(ip)%scale, ' m per file unit, interpolation ', &
       merge('cubic Hermite', 'linear       ', fp(ip)%interp == 3), ', time offset ', fp(ip)%t0, ' s'
    write(*,'(A)') '    the node coordinates x(t) = X0 + scale d(t); grid velocity = (x(n+1) - x(n)) / dt (ALE, GCL term on)'
  end subroutine frames_info

  !---------------------------------------------------------------------
  ! slot holding frame f (1-based; loaded from the file when absent)
  !---------------------------------------------------------------------
  integer function slot_of(ip, f)
    integer, intent(in) :: ip, f
    integer :: s, b, ios
    integer(8) :: pos
    do s = 1, NSLOT
       if (fp(ip)%slot_f(s) == f) then
          slot_of = s;  return
       end if
    end do
    s = fp(ip)%next_slot
    fp(ip)%next_slot = mod(fp(ip)%next_slot, NSLOT) + 1
    do b = 1, nblk
       if (.not. fb(b)%on) cycle
       if (fb(b)%part /= ip) cycle
       pos = fp(ip)%hdr_bytes + fp(ip)%frame_bytes*int(f - 1, 8) + fp(ip)%boff(b) + 1_8
       read(fp(ip)%iu, pos=pos, iostat=ios) fb(b)%d(1,:,:,:,s), fb(b)%d(2,:,:,:,s), fb(b)%d(3,:,:,:,s)
       if (ios /= 0) call frames_abort('read error in '//trim(fp(ip)%file))
    end do
    fp(ip)%slot_f(s) = f;  fp(ip)%nread = fp(ip)%nread + 1
    slot_of = s
  end function slot_of

  !---------------------------------------------------------------------
  ! interpolation weights at time t: frames j0..j3 (slots) and weights
  ! w(4) such that d(t) = sum w(m) d_frame(m)
  !---------------------------------------------------------------------
  subroutine weights_at(ip, t, sl, w)
    integer,  intent(in)  :: ip
    real(dp), intent(in)  :: t
    integer,  intent(out) :: sl(4)
    real(dp), intent(out) :: w(4)
    integer :: nf, j, jm, jp, jpp
    real(dp) :: s, tj, tj1, h, hm, hp, u, h00, h10, h01, h11, am(4), ap(4), tjm, tjpp
    nf = fp(ip)%nf
    s = t + fp(ip)%t0
    if (fp(ip)%periodic) then
       s = fp(ip)%tf(1) + modulo(s - fp(ip)%tf(1), fp(ip)%period)
    else
       s = min(max(s, fp(ip)%tf(1)), fp(ip)%tf(nf))
    end if
    ! interval j: tf(j) <= s < tf(j+1)  (j = nf: the wrap interval tf(nf) .. tf(1) + period)
    j = nf
    do jp = 1, nf - 1
       if (s < fp(ip)%tf(jp+1)) then
          j = jp;  exit
       end if
    end do
    if (j == nf .and. .not. fp(ip)%periodic) then
       sl = slot_of(ip, nf);  w = [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
       return
    end if
    tj = fp(ip)%tf(j)
    if (j < nf) then
       jp = j + 1;  tj1 = fp(ip)%tf(jp)
    else
       jp = 1;  tj1 = fp(ip)%tf(1) + fp(ip)%period
    end if
    h = tj1 - tj;  u = (s - tj)/h
    if (fp(ip)%interp == 1) then
       sl(1) = slot_of(ip, j);  sl(2) = slot_of(ip, jp);  sl(3) = sl(1);  sl(4) = sl(2)
       w = [1.0_dp - u, u, 0.0_dp, 0.0_dp]
       return
    end if
    ! cubic Hermite with three-point (non-uniform) finite-difference tangents; one-sided at the ends
    ! of a non-periodic set; wrapped around the period of a periodic one
    if (j > 1) then
       jm = j - 1;  tjm = fp(ip)%tf(jm)
    else if (fp(ip)%periodic) then
       jm = nf;  tjm = fp(ip)%tf(nf) - fp(ip)%period
    else
       jm = j;  tjm = tj
    end if
    if (jp < nf .and. j < nf) then
       jpp = jp + 1;  tjpp = fp(ip)%tf(jpp)
    else if (fp(ip)%periodic) then
       if (j == nf) then
          jpp = 2;  tjpp = fp(ip)%tf(2) + fp(ip)%period
       else                                    ! jp == nf
          jpp = 1;  tjpp = fp(ip)%tf(1) + fp(ip)%period
       end if
    else
       jpp = jp;  tjpp = tj1
    end if
    hm = tj - tjm;  hp = tjpp - tj1
    ! tangent at j (as weights on d_jm, d_j, d_jp) and at jp (on d_j, d_jp, d_jpp)
    am = 0.0_dp;  ap = 0.0_dp
    if (hm > 0.0_dp) then
       am(1) = -h/(hm*(hm + h));  am(2) = (h - hm)/(hm*h);  am(3) = hm/(h*(hm + h))
    else
       am(2) = -1.0_dp/h;  am(3) = 1.0_dp/h
    end if
    if (hp > 0.0_dp) then
       ap(2) = -hp/(h*(h + hp));  ap(3) = (hp - h)/(h*hp);  ap(4) = h/(hp*(h + hp))
    else
       ap(2) = -1.0_dp/h;  ap(3) = 1.0_dp/h
    end if
    h00 = (1.0_dp + 2.0_dp*u)*(1.0_dp - u)**2;  h10 = u*(1.0_dp - u)**2
    h01 = u*u*(3.0_dp - 2.0_dp*u);               h11 = u*u*(u - 1.0_dp)
    w = 0.0_dp
    w(2) = w(2) + h00;  w(3) = w(3) + h01
    w(1:3) = w(1:3) + h10*h*am(1:3)
    w(2:4) = w(2:4) + h11*h*ap(2:4)
    sl(1) = slot_of(ip, jm);  sl(2) = slot_of(ip, j);  sl(3) = slot_of(ip, jp);  sl(4) = slot_of(ip, jpp)
  end subroutine weights_at

  !---------------------------------------------------------------------
  ! node coordinates of the part's blocks at time t
  !---------------------------------------------------------------------
  subroutine nodes_at(ip, t)
    integer,  intent(in) :: ip
    real(dp), intent(in) :: t
    integer :: sl(4), b, i, j, k, m
    real(dp) :: w(4), x(3)
    call weights_at(ip, t, sl, w)
    do b = 1, nblk
       if (.not. fb(b)%on) cycle
       if (fb(b)%part /= ip) cycle
       do k = 1, blk(b)%nk
          do j = 1, blk(b)%nj
             do i = 1, blk(b)%ni
                x = 0.0_dp
                do m = 1, 4
                   if (w(m) /= 0.0_dp) x = x + w(m)*real(fb(b)%d(:, i, j, k, sl(m)), dp)
                end do
                x = fb(b)%X0(:, i, j, k) + fp(ip)%scale*x
                blk(b)%xn(i,j,k) = x(1);  blk(b)%yn(i,j,k) = x(2);  blk(b)%zn(i,j,k) = x(3)
             end do
          end do
       end do
    end do
    fp(ip)%t_last = t
  end subroutine nodes_at

  !=====================================================================
  ! Mode = 2 hook (sixdof_check_movement): nodes of every frame part at t
  !=====================================================================
  subroutine frames_check_nodes(t)
    real(dp), intent(in) :: t
    integer :: ip
    do ip = 1, np
       call nodes_at(ip, t)
    end do
  end subroutine frames_check_nodes

  !=====================================================================
  ! HOOK (fsi_morph): coordinates of the part's blocks at the target time
  ! t_now, discrete grid velocity at the half lattice of the owned blocks
  !=====================================================================
  subroutine frames_morph(k, restoring)
    integer, intent(in) :: k
    logical, intent(in) :: restoring
    integer :: ip, b, ib, ni, nj, nk, h1, h2, h3, hx, hy, hz, q
    real(dp) :: dtt, vn(3)
    ip = frames_part(k)
    if (ip == 0) return
    if (.not. own_ready) then
       allocate(own_blk(nblk));  own_blk = .false.
       do ib = 1, nown
          own_blk(dom(own_list(ib))%b) = .true.
       end do
       do b = 1, nblk
          if (fb(b)%on .and. own_blk(b)) then
             allocate(fb(b)%vh(3, 2*blk(b)%ni-1, 2*blk(b)%nj-1, 2*blk(b)%nk-1))
             fb(b)%vh = 0.0_dp
          end if
       end do
       own_ready = .true.
    end if
    call nodes_at(ip, t_now)
    dtt = max(dt_dim, 1.0e-300_dp)
    do b = 1, nblk
       if (.not. fb(b)%on) cycle
       if (fb(b)%part /= ip) cycle
       ni = blk(b)%ni;  nj = blk(b)%nj;  nk = blk(b)%nk
       if (restoring) then
          ! restart: the restored grid is the level-n grid of the first step; no velocity yet
          fb(b)%xprev(1,:,:,:) = blk(b)%xn;  fb(b)%xprev(2,:,:,:) = blk(b)%yn;  fb(b)%xprev(3,:,:,:) = blk(b)%zn
          if (allocated(fb(b)%vh)) fb(b)%vh = 0.0_dp
          cycle
       end if
       if (.not. allocated(fb(b)%vh)) cycle
       do h3 = 1, 2*nk - 1
          do h2 = 1, 2*nj - 1
             do h1 = 1, 2*ni - 1
                vn = 0.0_dp;  q = 0
                do hz = (h3+1)/2, min((h3+2)/2, nk)
                   do hy = (h2+1)/2, min((h2+2)/2, nj)
                      do hx = (h1+1)/2, min((h1+2)/2, ni)
                         vn = vn + ([blk(b)%xn(hx,hy,hz), blk(b)%yn(hx,hy,hz), blk(b)%zn(hx,hy,hz)] - fb(b)%xprev(:, hx, hy, hz))
                         q = q + 1
                      end do
                   end do
                end do
                fb(b)%vh(:, h1, h2, h3) = vn/q/dtt/Uch
             end do
          end do
       end do
    end do
  end subroutine frames_morph

  ! HOOK (fsi_hvel): discrete grid velocity (non-dim) at half index h of block b
  function frames_hvel(b, h) result(ug)
    integer, intent(in) :: b, h(3)
    real(dp) :: ug(3)
    integer :: hh(3)
    ug = 0.0_dp
    if (.not. allocated(fb(b)%vh)) return
    hh(1) = min(max(h(1), 1), 2*blk(b)%ni - 1)
    hh(2) = min(max(h(2), 1), 2*blk(b)%nj - 1)
    hh(3) = min(max(h(3), 1), 2*blk(b)%nk - 1)
    ug = fb(b)%vh(:, hh(1), hh(2), hh(3))
  end function frames_hvel

  ! end of the step (fsi_finish_step): the coordinates become the previous ones
  subroutine frames_commit
    integer :: b
    if (np == 0) return
    do b = 1, nblk
       if (.not. fb(b)%on) cycle
       fb(b)%xprev(1,:,:,:) = blk(b)%xn;  fb(b)%xprev(2,:,:,:) = blk(b)%yn;  fb(b)%xprev(3,:,:,:) = blk(b)%zn
    end do
  end subroutine frames_commit

end module morph_frames
