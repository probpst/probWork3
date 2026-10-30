!=======================================================================
! restart.F90 -- restart files (update 13b)
!
!   <case>_NNNNNNNN.rsto  (working folder, binary, one file written by
!   rank 0) holds everything the time marching needs to continue from
!   step NNNNNNNN exactly:
!     header : format tag, case name, step, time, nblk, ndom, phases,
!              turbulence flag, number of 6-DoF bodies, monitor sums
!     bodies : the 42-value state of every body (position, velocities,
!              accelerations, last forces / moments)
!     domains: for every sub-domain the bounds and the arrays
!              U1 U2 U3 (contravariant fluxes), u v w (Cartesian faces),
!              p, [alpha], [vapour fraction], [tk tw mut]  - with their ghost layers
!   The sub-domain decomposition is that of the csv file, so the file is
!   independent of the number of MPI ranks (a domain is sent to whoever
!   owns it on reading).
!
!   Mode = 1 (name list): the latest file is found through
!   <case>_restart.log (one line per restart file written) and read
!   after the initial assembly; the 6-DoF grids are then moved to the
!   restored poses and the overset is re-assembled; the monitor and
!   data files are continued (append).  Restart_Intvl (name list,
!   $Export_Interval) sets the writing interval; the last step of a run
!   is always written when Restart_Intvl > 0.
!=======================================================================
module restart
  use global_vars
  use mpi_vars
  use sixdof, only: sixdof_nbody, sixdof_get_state, sixdof_set_state, sixdof_restore_pose
  use waves, only: wave_update, wave_on
  use pid_control, only: pid_nstate, pid_get_state, pid_set_state
  use statistics, only: stat_restart_write, stat_restart_read
  use fsi, only: fsi_restart_write, fsi_restart_read
  implicit none
  private
  public :: write_restart, read_restart, restart_available

  character(len=*), parameter :: TAG = 'probCFD_rsto_v3'   ! v3 (update 16): + energy flag / temperature (last array)
  character(len=*), parameter :: TAG2 = 'probCFD_rsto_v2'  ! v2 (update 15): + cavitation flag / vapour fraction
  character(len=*), parameter :: TAG1 = 'probCFD_rsto_v1'
  integer :: file_energy = 0                               ! energy flag of the file being read

contains

  !---------------------------------------------------------------------
  subroutine restart_name(istep, fname)
    integer, intent(in) :: istep
    character(len=*), intent(out) :: fname
    write(fname,'(A,A,I8.8,A)') trim(casename), '_', istep, '.rsto'
  end subroutine restart_name

  !---------------------------------------------------------------------
  ! the step of the latest restart file (0: none) - last line of
  ! <case>_restart.log whose file still exists
  !---------------------------------------------------------------------
  integer function restart_available() result(istep)
    integer :: iu, ios, n
    logical :: ex
    character(len=256) :: fname
    istep = 0
    if (myrank /= 0) return
    open(newunit=iu, file=trim(casename)//'_restart.log', status='old', action='read', iostat=ios)
    if (ios /= 0) return
    do
       read(iu, *, iostat=ios) n
       if (ios /= 0) exit
       call restart_name(n, fname)
       inquire(file=trim(fname), exist=ex)
       if (ex) istep = n
    end do
    close(iu)
  end function restart_available

  !---------------------------------------------------------------------
  subroutine write_restart(istep, t, vliq0, cliptot)
    integer,  intent(in) :: istep
    real(dp), intent(in) :: t, vliq0, cliptot
    integer :: iu, d, k, nb, ierr, lo(3), hi(3), n, ia, np
    integer :: stat(MPI_STATUS_SIZE)
    real(dp) :: st(42), pst(8)
    real(dp), allocatable :: buf(:,:,:)
    real(dp), pointer :: a(:,:,:)
    character(len=256) :: fname
    nb = sixdof_nbody()
    call restart_name(istep, fname)
    if (myrank == 0) then
       open(newunit=iu, file=trim(fname), form='unformatted', status='replace', action='write')
       write(iu) TAG
       write(iu) casename
       write(iu) istep, t, nblk, ndom, Num_of_Phase, merge(1, 0, turb_on), nb, dt_dim, merge(1, 0, cav_on), merge(1, 0, energy_on)
       write(iu) vliq0, cliptot
       do k = 1, nb
          call sixdof_get_state(k, st)
          write(iu) st
       end do
    end if
    file_energy = merge(1, 0, energy_on)
    do d = 1, ndom
       do ia = 1, narrays()
          if (dom(d)%owner == myrank) then
             a => buf_target(d, ia);  lo = lbound(a);  hi = ubound(a)
             if (myrank /= 0) then
                call MPI_Send(lo, 3, MPI_INTEGER, 0, 700, MPI_COMM_WORLD, ierr)
                call MPI_Send(hi, 3, MPI_INTEGER, 0, 701, MPI_COMM_WORLD, ierr)
                call MPI_Send(a, product(hi-lo+1), MPI_DOUBLE_PRECISION, 0, 702, MPI_COMM_WORLD, ierr)
             end if
          end if
          if (myrank == 0) then
             if (dom(d)%owner == 0) then
                a => buf_target(d, ia);  lo = lbound(a);  hi = ubound(a)
                write(iu) d, ia, lo, hi
                write(iu) a
             else
                call MPI_Recv(lo, 3, MPI_INTEGER, dom(d)%owner, 700, MPI_COMM_WORLD, stat, ierr)
                call MPI_Recv(hi, 3, MPI_INTEGER, dom(d)%owner, 701, MPI_COMM_WORLD, stat, ierr)
                allocate(buf(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3)))
                n = product(hi-lo+1)
                call MPI_Recv(buf, n, MPI_DOUBLE_PRECISION, dom(d)%owner, 702, MPI_COMM_WORLD, stat, ierr)
                write(iu) d, ia, lo, hi
                write(iu) buf
                deallocate(buf)
             end if
          end if
       end do
    end do
    if (myrank == 0) then
       ! PID controller states (update 21): appended after the fields; absent in older files
       np = pid_nstate()
       write(iu) np
       do k = 1, np
          call pid_get_state(k, pst)
          write(iu) pst
       end do
    end if
    call stat_restart_write(iu)                          ! statistics accumulators (update 23; all ranks take part)
    call fsi_restart_write(iu)                           ! structural states of the flexible parts (update 24)
    if (myrank == 0) then
       close(iu)
       open(newunit=iu, file=trim(casename)//'_restart.log', status='unknown', position='append', action='write')
       write(iu,'(I0)') istep
       close(iu)
       write(*,'(A)') '  wrote restart file: '//trim(fname)
    end if
  end subroutine write_restart

  !---------------------------------------------------------------------
  ! read the restart file of step istep0 (found by restart_available on
  ! rank 0, broadcast by the caller); returns t, the monitor sums; the
  ! bodies and the fields are restored, the grids moved
  !---------------------------------------------------------------------
  subroutine read_restart(istep0, t, vliq0, cliptot)
    integer,  intent(in)  :: istep0
    real(dp), intent(out) :: t, vliq0, cliptot
    integer :: iu, d, k, nb, ierr, lo(3), hi(3), n, ia, ios, dd, iaa, lo2(3), hi2(3)
    integer :: hdr(8)
    integer :: stat(MPI_STATUS_SIZE)
    real(dp) :: st(42), rr(4), dt_file, pst(8)
    integer :: np
    logical :: sok, fok
    real(dp), allocatable :: buf(:,:,:)
    real(dp), pointer :: a(:,:,:)
    character(len=len(TAG)) :: tag_in
    character(len=len(casename)) :: cn
    character(len=256) :: fname
    nb = sixdof_nbody()
    call restart_name(istep0, fname)
    if (myrank == 0) then
       open(newunit=iu, file=trim(fname), form='unformatted', status='old', action='read', iostat=ios)
       if (ios /= 0) call fail('cannot open '//trim(fname))
       read(iu, iostat=ios) tag_in
       if (ios /= 0 .or. (tag_in /= TAG .and. tag_in /= TAG2 .and. tag_in /= TAG1)) &
          call fail(trim(fname)//' is not a probCFD restart file (tag)')
       read(iu) cn
       hdr(7) = 0;  hdr(8) = 0
       if (tag_in == TAG) then
          read(iu) hdr(1), rr(1), hdr(2), hdr(3), hdr(4), hdr(5), hdr(6), dt_file, hdr(7), hdr(8)
       else if (tag_in == TAG2) then
          read(iu) hdr(1), rr(1), hdr(2), hdr(3), hdr(4), hdr(5), hdr(6), dt_file, hdr(7)
       else
          read(iu) hdr(1), rr(1), hdr(2), hdr(3), hdr(4), hdr(5), hdr(6), dt_file
       end if
       if (hdr(7) /= merge(1, 0, cav_on)) call fail('cavitation on/off differs from the restart file')
       ! energy (update 16): a file without T may be restarted with SW_Energy = 1 (T starts at T_init);
       ! a file with T restarted with SW_Energy = 0 ignores the temperature records (they come last)
       if (hdr(8) == 0 .and. energy_on .and. myrank == 0) write(*,'(A)') &
          '  NOTE: the restart file has no temperature field - the energy equation starts from T = T_init'
       if (hdr(8) == 1 .and. .not. energy_on .and. myrank == 0) write(*,'(A)') &
          '  NOTE: the restart file carries a temperature field that is not used (SW_Energy = 0)'
       read(iu) rr(2), rr(3)
       if (hdr(1) /= istep0) call fail('step in the file differs from the file name')
       if (hdr(2) /= nblk .or. hdr(3) /= ndom) call fail('grid / csv splitting differs from the restart file')
       if (hdr(4) /= Num_of_Phase) call fail('Num_of_Phase differs from the restart file')
       if (hdr(5) /= merge(1, 0, turb_on)) call fail('turbulence on/off differs from the restart file')
       if (hdr(6) /= nb) call fail('number of 6-DoF bodies differs from the restart file')
       rr(4) = 0.0_dp
    end if
    call MPI_Bcast(rr, 4, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(hdr, 8, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    file_energy = merge(1, 0, hdr(8) == 1 .and. energy_on)      ! temperature records to read
    t = rr(1);  vliq0 = rr(2);  cliptot = rr(3)
    do k = 1, nb
       if (myrank == 0) read(iu) st
       call MPI_Bcast(st, 42, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
       call sixdof_set_state(k, st)
    end do
    do d = 1, ndom
       do ia = 1, narrays()
          if (myrank == 0) then
             read(iu, iostat=ios) dd, iaa, lo, hi
             if (ios /= 0 .or. dd /= d .or. iaa /= ia) call fail('unexpected record order in '//trim(fname))
             if (dom(d)%owner == 0) then
                a => buf_target(d, ia);  lo2 = lbound(a);  hi2 = ubound(a)
                if (any(lo /= lo2) .or. any(hi /= hi2)) call fail('array bounds differ from the restart file')
                read(iu) a
             else
                allocate(buf(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3)))
                read(iu) buf
                call MPI_Send(lo, 3, MPI_INTEGER, dom(d)%owner, 710, MPI_COMM_WORLD, ierr)
                call MPI_Send(hi, 3, MPI_INTEGER, dom(d)%owner, 711, MPI_COMM_WORLD, ierr)
                call MPI_Send(buf, product(hi-lo+1), MPI_DOUBLE_PRECISION, dom(d)%owner, 712, MPI_COMM_WORLD, ierr)
                deallocate(buf)
             end if
          else if (dom(d)%owner == myrank) then
             call MPI_Recv(lo, 3, MPI_INTEGER, 0, 710, MPI_COMM_WORLD, stat, ierr)
             call MPI_Recv(hi, 3, MPI_INTEGER, 0, 711, MPI_COMM_WORLD, stat, ierr)
             a => buf_target(d, ia);  lo2 = lbound(a);  hi2 = ubound(a)
             if (any(lo /= lo2) .or. any(hi /= hi2)) call fail('array bounds differ from the restart file')
             n = product(hi-lo+1)
             call MPI_Recv(a, n, MPI_DOUBLE_PRECISION, 0, 712, MPI_COMM_WORLD, stat, ierr)
          end if
       end do
    end do
    ! PID controller states (update 21), appended after the fields; an older file has none
    np = 0
    if (myrank == 0) then
       read(iu, iostat=ios) np
       if (ios /= 0) np = 0
       if (np /= pid_nstate()) then
          if (np > 0 .or. pid_nstate() > 0) write(*,'(A,I0,A,I0,A)') '  NOTE: the restart file holds ', np, &
             ' PID controller states, this run has ', pid_nstate(), ' - the controllers start from their initial state'
          np = 0
       end if
    end if
    call MPI_Bcast(np, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    do k = 1, np
       if (myrank == 0) read(iu) pst
       call MPI_Bcast(pst, 8, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
       call pid_set_state(k, pst)
    end do
    call stat_restart_read(iu, sok)                       ! statistics accumulators (update 23)
    call fsi_restart_read(iu, fok)                        ! structural states of the flexible parts (update 24)
    if (myrank == 0) then
       close(iu)
       write(*,'(A,I0,A,ES14.6,A)') '  restart file read: '//trim(fname)//' (step ', istep0, ', t = ', t, ' s)'
       if (abs(dt_file - dt_dim) > 1.0e-12_dp*dt_dim) write(*,'(A,ES11.4,A,ES11.4,A)') &
          '  NOTE: Del_T of the restart file ', dt_file, ' s, of this run ', dt_dim, ' s'
       if (np > 0) write(*,'(A,I0,A)') '  PID controller states restored: ', np
    end if
    ! boundary values of the restored field at the restored time (update 20):
    ! the start-up ramp and the incident wave are those of t, not of t = 0 -
    ! the restore applies the Dirichlet values (inlets, moving walls u_init
    ! (1 - ramp) + body velocity) to the level-n field of the first step
    bc_ramp = 1.0_dp
    if (t_ramp > 0.0_dp) bc_ramp = min(t/t_ramp, 1.0_dp)
    if (wave_on) call wave_update(t)
    ! the grids of moving bodies: to the restored poses (+ overset re-assembly)
    call sixdof_restore_pose(t)
  end subroutine read_restart

  !---------------------------------------------------------------------
  ! the arrays saved per domain, in order
  !---------------------------------------------------------------------
  integer function narrays()
    narrays = 7
    if (Num_of_Phase == 2) narrays = narrays + 1
    if (cav_on) narrays = narrays + 1
    if (turb_on) narrays = narrays + 3
    if (file_energy == 1) narrays = narrays + 1       ! temperature: the LAST array (update 16)
  end function narrays

  ! pointer-like access: which array of domain d is number ia
  function buf_target(d, ia) result(a)
    integer, intent(in) :: d, ia
    real(dp), pointer :: a(:,:,:)
    integer :: m
    m = ia
    select case (m)
    case (1);  a => dom(d)%U1
    case (2);  a => dom(d)%U2
    case (3);  a => dom(d)%U3
    case (4);  a => dom(d)%u
    case (5);  a => dom(d)%v
    case (6);  a => dom(d)%w
    case (7);  a => dom(d)%p
    case default
       if (Num_of_Phase == 2) then
          if (m == 8) then
             a => dom(d)%alpha;  return
          end if
          m = m - 1
       end if
       if (cav_on) then
          if (m == 8) then
             a => dom(d)%vap;  return
          end if
          m = m - 1
       end if
       if (turb_on) then
          select case (m)
          case (8);  a => dom(d)%tk;  return
          case (9);  a => dom(d)%tw;  return
          case (10); a => dom(d)%mut;  return
          end select
          m = m - 3
       end if
       a => dom(d)%temp                                ! m == 8: temperature (update 16)
    end select
  end function buf_target

  subroutine fail(msg)
    character(len=*), intent(in) :: msg
    integer :: ierr
    write(*,'(A)') ' ERROR (restart): '//msg
    flush(6)
    call MPI_Abort(MPI_COMM_WORLD, 41, ierr)
  end subroutine fail

end module restart
