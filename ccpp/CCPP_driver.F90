module CCPP_driver

  use ccpp_types,         only: ccpp_t

  use ccpp_static_api,    only: ccpp_physics_init,                   &
                                ccpp_physics_timestep_init,          &
                                ccpp_physics_run,                    &
                                ccpp_physics_timestep_finalize,      &
                                ccpp_physics_finalize

  use CCPP_data,          only: cdata_tile,                          &
                                cdata_domain,                        &
                                cdata_block,                         &
                                ccpp_suite,                          &
                                GFS_control,                         &
                                GFS_Intdiag,                         &
                                GFS_Interstitial

  implicit none

!--------------------------------------------------------!
!  Pointer to CCPP containers defined in CCPP_data       !
!--------------------------------------------------------!
  type(ccpp_t), pointer :: cdata => null()

!--------------------------------------------------------!
!  Flag for non-uniform block sizes (last block smaller) !
!  and number of OpenMP threads (with special thread     !
!  number nthrdsX in case of non-uniform block sizes)    !
!--------------------------------------------------------!
  logical :: non_uniform_blocks
  integer :: nthrds, nthrdsX

!----------------
! Public Entities
!----------------
! functions
  public CCPP_step
! module variables
  public non_uniform_blocks

  CONTAINS
!*******************************************************************************************

!-------------------------------
!  CCPP step
!-------------------------------
  subroutine CCPP_step (step, nblks, ierr, dycore)

#ifdef _OPENMP
    use omp_lib
#endif

    implicit none

    character(len=*),         intent(in)  :: step
    integer,                  intent(in)  :: nblks
    integer,                  intent(out) :: ierr
    character(len=*),         intent(in)  :: dycore
    ! Local variables
    integer :: nb, nt
    integer :: ierr2
    integer :: kdt_iau
    logical :: iauwindow_center
    ! DH* 20210104 - remove kdt_rad when code to clear diagnostic buckets is removed
    integer :: kdt_rad

    ierr = 0

    ! CCPP Framework init (same for all dynamical cores)
    if (trim(step)=="init") then

      ! Get and set number of OpenMP threads (module
      ! variable) that are available to run physics
#ifdef _OPENMP
      nthrds = omp_get_max_threads()
#else
      nthrds = 1
#endif

      ! For non-uniform blocksizes, we use index nthrds+1
      ! for the interstitial data type with different length
      if (non_uniform_blocks) then
        nthrdsX = nthrds+1
      else
        nthrdsX = nthrds
      end if

      ! For physics running over the entire domain, block, chunk and thread
      ! numbers are not used; set to safe values
      cdata_domain%blk_no = 1
      cdata_domain%chunk_no = 1
      cdata_domain%thrd_no = 1
      cdata_domain%thrd_cnt = 1

      ! Allocate cdata structures for blocks and threads
      if (.not.allocated(cdata_block)) allocate(cdata_block(1:nblks,1:nthrdsX))

      ! Loop over all blocks and threads
      do nt=1,nthrdsX
        do nb=1,nblks
          ! Assign the correct block, chunk and thread numbers
          ! Note that we can use block number as chunk number
          cdata_block(nb,nt)%blk_no = nb
          cdata_block(nb,nt)%chunk_no = nb
          cdata_block(nb,nt)%thrd_no = nt
          cdata_block(nb,nt)%thrd_cnt = nthrdsX
        end do
      end do
    ! Physics init (same for all dynamical cores)
    else if (trim(step)=="physics_init") then

      ! Since the physics init step is independent of the blocking structure,
      ! we can use cdata_domain. And since we don't use threading on the host
      ! model side, we can allow threading inside the physics init routines.
      GFS_control%nthreads = nthrds

      call ccpp_physics_init(cdata_domain, suite_name=trim(ccpp_suite), ierr=ierr)
      if (ierr/=0) then
        write(0,'(a)') "An error occurred in ccpp_physics_init"
        write(0,'(a)') trim(cdata_domain%errmsg)
        return
      end if

    ! Timestep init = time_vary (dycore specific)
    else if (trim(step)=="timestep_init") then

      ! Since the physics timestep init step is independent of the blocking structure,
      ! we can use cdata_domain. And since we don't use threading on the host
      ! model side, we can allow threading inside the timestep init (time_vary) routines.
      GFS_control%nthreads = nthrds

      call ccpp_physics_timestep_init(cdata_domain, suite_name=trim(ccpp_suite), group_name="time_vary", ierr=ierr)
      if (ierr/=0) then
        write(0,'(a)') "An error occurred in ccpp_physics_timestep_init for group time_vary"
        write(0,'(a)') trim(cdata_domain%errmsg)
        return
      end if

      if (trim(dycore)=='fv3') then
         ! call timestep_init for "phys_ps"---required for Land IAU
         call ccpp_physics_timestep_init(cdata_domain, suite_name=trim(ccpp_suite),group_name="phys_ps", ierr=ierr)
         if (ierr/=0) then
            write(0,'(a)') "An error occurred in ccpp_physics_timestep_init for group phys_ps"
            write(0,'(a)') trim(cdata_domain%errmsg)
            return
         end if

         ! call timestep_init for "phys_ts"---required for Land IAU
         call ccpp_physics_timestep_init(cdata_domain, suite_name=trim(ccpp_suite),group_name="phys_ts", ierr=ierr)
         if (ierr/=0) then
            write(0,'(a)') "An error occurred in ccpp_physics_timestep_init for group phys_ts"
            write(0,'(a)') trim(cdata_domain%errmsg)
            return
         end if
      endif

      if (trim(dycore)=='mpas') then
         ! Physics group
         call ccpp_physics_timestep_init(cdata_domain, suite_name=trim(ccpp_suite),group_name="physics", ierr=ierr)
         if (ierr/=0) then
            write(0,'(a)') "An error occurred in ccpp_physics_timestep_init for group physics"
            write(0,'(a)') trim(cdata_domain%errmsg)
            return
         end if

         call ccpp_physics_timestep_init(cdata_domain, suite_name=trim(ccpp_suite),group_name="microphysics", ierr=ierr)
         if (ierr/=0) then
            write(0,'(a)') "An error occurred in ccpp_physics_timestep_init for group microphysics"
            write(0,'(a)') trim(cdata_domain%errmsg)
            return
         end if
      end if

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! DH* 20210104 - this block of code will be removed once the CCPP framework    !
      ! fully supports handling diagnostics through its metadata, work in progress   !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !--- determine if radiation diagnostics buckets need to be cleared
      if (nint(GFS_control%fhzero*3600) >= nint(max(GFS_control%fhswr,GFS_control%fhlwr))) then
        if (mod(GFS_control%kdt,GFS_control%nszero) == 1) then
          call GFS_Intdiag%rad_zero(GFS_control)
        endif
      else
        kdt_rad = nint(min(GFS_control%fhswr,GFS_control%fhlwr)/GFS_control%dtp)
        if (mod(GFS_control%kdt,kdt_rad) == 1) then
          call GFS_Intdiag%rad_zero(GFS_control)
        endif
      endif

      !--- determine if physics diagnostics buckets need to be cleared
      iauwindow_center = .false.
      if (GFS_control%iau_offset > 0) then
        kdt_iau = nint(GFS_control%iau_offset*3600./GFS_control%dtp)
        if (GFS_control%kdt-1 == kdt_iau) then
          iauwindow_center = .true.
          if( GFS_control%me == 0)print *,'in ccpp step vary, iauwindow_center=',iauwindow_center,&
            'kdt=',GFS_control%kdt,'dtp=',GFS_control%dtp,'iau_offset=',GFS_control%iau_offset
        endif
      endif
      if ((mod(GFS_control%kdt-1,GFS_control%nszero)) == 0) then
        call GFS_Intdiag%phys_zero(GFS_control, iauwindow_center=iauwindow_center)
      endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! *DH 20210104                                                                 !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Radiation, physics and and stochastic physics - threaded regions using blocked data structures
    else if (trim(step)=="radiation" .or. trim(step)=="physics" .or. trim(step)=="stochastics" .or. trim(step)=="microphysics") then

      ! Set number of threads available to physics schemes to one,
      ! because threads are used on the host model side for blocking
      GFS_control%nthreads = 1

!$OMP parallel num_threads (nthrds)                        &
!$OMP          default (none)                              &
!$OMP          shared (nblks, cdata_block, ccpp_suite,     &
!$OMP                  step, GFS_Control, GFS_Interstitial,&
!$OMP                  dycore)                             &
!$OMP          private (nb, nt, ierr2)                     &
!$OMP          reduction (+:ierr)
#ifdef _OPENMP
      nt = omp_get_thread_num()+1
#else
      nt = 1
#endif
!$OMP do schedule (dynamic,1)
      do nb = 1,nblks
         ! Allocate physics interstitals for current thread
         call GFS_Interstitial(nt)%create(ixs=GFS_control%chunk_begin(nb), ixe=GFS_control%chunk_end(nb), model=GFS_control)
        !--- Call CCPP radiation/physics/stochastics group
        if (trim(step)=="physics") then
           if (trim(dycore)=="fv3") then
              ! Reset GFS_Interstitial DDT fields for this thread
              call GFS_Interstitial(nt)%reset(GFS_control)
              ! Process-split physics
              call ccpp_physics_run(cdata_block(nb,nt), suite_name=trim(ccpp_suite), group_name="phys_ps", ierr=ierr2)
              if (ierr2/=0) then
                 write(0,'(2a,3(a,i4),a)') "An error occurred in ccpp_physics_run for group ", "phys_ps", &
                                           ", block/chunk ", nb, " and thread ", nt, " (nt=", nt, "):"
                 write(0,'(a)') trim(cdata_block(nb,nt)%errmsg)
                 ierr = ierr + ierr2
              endif
              ! Time-split physics
              call ccpp_physics_run(cdata_block(nb,nt), suite_name=trim(ccpp_suite), group_name="phys_ts", ierr=ierr2)
              if (ierr2/=0) then
                 write(0,'(2a,3(a,i4),a)') "An error occurred in ccpp_physics_run for group ", "phys_ts", &
                                           ", block/chunk ", nb, " and thread ", nt, " (nt=", nt, "):"
                 write(0,'(a)') trim(cdata_block(nb,nt)%errmsg)
                 ierr = ierr + ierr2
              endif
           endif
           if (trim(dycore)=="mpas") then
              ! Physics
              call ccpp_physics_run(cdata_block(nb,nt), suite_name=trim(ccpp_suite), group_name="physics", ierr=ierr2)
              if (ierr2/=0) then
                 write(0,'(2a,3(a,i4),a)') "An error occurred in ccpp_physics_run for group ", "physics", &
                                           ", block/chunk ", nb, " and thread ", nt, " (nt=", nt, "):"
                 write(0,'(a)') trim(cdata_block(nb,nt)%errmsg)
                 ierr = ierr + ierr2
              endif
           endif
        else
           if (trim(step)=="radiation") then
              ! Reset GFS_Interstitial DDT fields for this thread
              call GFS_Interstitial(nt)%reset(GFS_control)
           endif
           ! Radiation
           call ccpp_physics_run(cdata_block(nb,nt), suite_name=trim(ccpp_suite), group_name=trim(step), ierr=ierr2)
           if (ierr2/=0) then
              write(0,'(2a,3(a,i4),a)') "An error occurred in ccpp_physics_run for group ", trim(step), &
                   ", block/chunk ", nb, " and thread ", nt, " (nt=", nt, "):"
              write(0,'(a)') trim(cdata_block(nb,nt)%errmsg)
              ierr = ierr + ierr2
           endif
           ! Microphysics (MPAS only)
           if (trim(step)=="microphysics") then
              if (trim(dycore)=="mpas") then
                 call ccpp_physics_run(cdata_block(nb,nt), suite_name=trim(ccpp_suite), group_name="microphysics", ierr=ierr2)
                 if (ierr2/=0) then
                    write(0,'(2a,3(a,i4),a)') "An error occurred in ccpp_physics_run for group ", "microphysics", &
                                              ", block/chunk ", nb, " and thread ", nt, " (nt=", nt, "):"
                    write(0,'(a)') trim(cdata_block(nb,nt)%errmsg)
                    ierr = ierr + ierr2
                 endif
              else
                 write(0,'(a)') "An error occurred in ccpp_physics_run for group microphysics. Group microphysics only valid with MPAS dycore."
                 ierr = ierr + 1
              endif
           endif
        endif
        call GFS_Interstitial(nt)%destroy(GFS_control)
     end do
!$OMP end do

!$OMP end parallel
      if (ierr/=0) return

    ! Timestep finalize = time_vary
    else if (trim(step)=="timestep_finalize") then

      ! Since the physics timestep finalize step is independent of the blocking structure,
      ! we can use cdata_domain. And since we don't use threading on the host model side,
      ! we can allow threading inside the timestep finalize (time_vary) routines.
      GFS_control%nthreads = nthrds

      call ccpp_physics_timestep_finalize(cdata_domain, suite_name=trim(ccpp_suite), group_name="time_vary", ierr=ierr)
      if (ierr/=0) then
        write(0,'(a)') "An error occurred in ccpp_physics_timestep_finalize for group time_vary"
        write(0,'(a)') trim(cdata_domain%errmsg)
        return
      end if

      if (trim(dycore)=='fv3') then
         ! call timestep_finalize for "phys_ps"---required for Land IAU
         call ccpp_physics_timestep_finalize(cdata_domain, suite_name=trim(ccpp_suite), group_name="phys_ps", ierr=ierr)
         if (ierr/=0) then
            write(0,'(a)') "An error occurred in ccpp_physics_timestep_finalize for group phys_ps"
            write(0,'(a)') trim(cdata_domain%errmsg)
            return
         end if

         ! call timestep_finalize for "phys_ts"---required for Land IAU
         call ccpp_physics_timestep_finalize(cdata_domain, suite_name=trim(ccpp_suite), group_name="phys_ts", ierr=ierr)
         if (ierr/=0) then
            write(0,'(a)') "An error occurred in ccpp_physics_timestep_finalize for group phys_ts"
            write(0,'(a)') trim(cdata_domain%errmsg)
            return
         end if
      endif
      if (trim(dycore)=='mpas') then
         call ccpp_physics_timestep_finalize(cdata_domain, suite_name=trim(ccpp_suite), group_name="physics", ierr=ierr)
         if (ierr/=0) then
            write(0,'(a)') "An error occurred in ccpp_physics_timestep_finalize for group physics"
            write(0,'(a)') trim(cdata_domain%errmsg)
            return
         end if

         call ccpp_physics_timestep_finalize(cdata_domain, suite_name=trim(ccpp_suite), group_name="microphysics", ierr=ierr)
         if (ierr/=0) then
            write(0,'(a)') "An error occurred in ccpp_physics_timestep_finalize for group microphysics"
            write(0,'(a)') trim(cdata_domain%errmsg)
            return
         end if
      endif

    ! Physics finalize (same for all dynamical cores)
    else if (trim(step)=="physics_finalize") then

      ! Since the physics finalize step is independent of the blocking structure,
      ! we can use cdata_domain. And since we don't use threading on the host
      ! model side, we can allow threading inside the physics finalize routines.
      GFS_control%nthreads = nthrds

      call ccpp_physics_finalize(cdata_domain, suite_name=trim(ccpp_suite), ierr=ierr)
      if (ierr/=0) then
        write(0,'(a)') "An error occurred in ccpp_physics_finalize"
        write(0,'(a)') trim(cdata_domain%errmsg)
        return
      end if

    ! Finalize (same for all dynamical cores)
    else if (trim(step)=="finalize") then
      ! Deallocate cdata structure for blocks and threads
      if (allocated(cdata_block)) deallocate(cdata_block)

    else

      write(0,'(2a)') 'Error, undefined CCPP step ', trim(step)
      ierr = 1
      return

    end if

  end subroutine CCPP_step

end module CCPP_driver
