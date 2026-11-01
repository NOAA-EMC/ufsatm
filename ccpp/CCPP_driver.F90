module CCPP_driver

  use ufs_ccpp_cap,       only: ccpp_register,               &
                                ccpp_init,                   &
                                ccpp_physics_init,           &
                                ccpp_physics_timestep_init,  &
                                ccpp_physics_run,            &
                                ccpp_physics_timestep_final, &
                                ccpp_physics_final,          &
                                ccpp_final
  use CCPP_data,          only: GFS_control,                 &
                                GFS_Intdiag,                 &
                                GFS_Interstitial
  use iso_fortran_env,    only: error_unit

  implicit none

  !--------------------------------------------------------!
  ! CCPP control (mandatory) data.
  !--------------------------------------------------------!
  integer :: mythread
  integer :: nthreads
  integer :: nphys_threads
  integer :: lb
  integer :: ub
  integer :: errflg
  character(len=512) :: errmsg
  character(len=256) :: ccpp_suite='undefined'
  character(len=256) :: group_name='undefined'

!--------------------------------------------------------!
!  Number of OpenMP threads                             !
!--------------------------------------------------------!
  integer :: nthrds

  !--------------------------------------------------------!
  ! Public Entities
  !--------------------------------------------------------!
  public CCPP_step

CONTAINS
  !--------------------------------------------------------!
  !  CCPP step
  !--------------------------------------------------------!
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
    integer :: kdt_iau
    logical :: iauwindow_center
    ! DH* 20210104 - remove kdt_rad when code to clear diagnostic buckets is removed
    integer :: kdt_rad

    ierr = 0

    ! CCPP register (same for all dynamical cores)
    if (trim(step)=="register") then
       call ccpp_register(ccpp_suite=trim(ccpp_suite), errmsg=errmsg, errflg=errflg)
       if (errflg/=0) then
          write(error_unit,'(a,i0,a)') 'An error occurred in ccpp_register: ' // trim(errmsg) // '. Exiting...'
          ierr=errflg
          return
       end if
       
    ! CCPP Framework init (same for all dynamical cores)
    else if (trim(step)=="init") then

       ! Get and set number of OpenMP threads (module
       ! variable) that are available to run physics
#ifdef _OPENMP
       nthrds = omp_get_max_threads()
#else
       nthrds = 1
#endif
      
       call ccpp_init(ccpp_suite=trim(ccpp_suite), errmsg=errmsg, errflg=errflg)
       if (errflg/=0) then
          write(error_unit,'(a,i0,a)') 'An error occurred in ccpp_init: ' // trim(errmsg) // '. Exiting...'
          ierr=errflg
          return
       end if

    ! Physics init (same for all dynamical cores)
    else if (trim(step)=="physics_init") then

       call ccpp_physics_init( ccpp_suite=trim(ccpp_suite), group_name='all', &
            errmsg=errmsg, errflg=errflg, lb=1, ub=GFS_control%ncols,         &
            mythread=1, nthreads=1, nphys_threads=nthrds)
       if (errflg/=0) then
          write(error_unit,'(a,i0,a)') 'An error occurred in ccpp_physics_init: ' // trim(errmsg) // '. Exiting...'
          ierr=errflg
          return
       end if

    ! Timestep init = time_vary (dycore specific)
    else if (trim(step)=="timestep_init") then

       call ccpp_physics_timestep_init(ccpp_suite=trim(ccpp_suite), group_name='all', &
            errmsg=errmsg, errflg=errflg, lb=1, ub=GFS_control%ncols,         &
            mythread=1, nthreads=1, nphys_threads=nthrds)
       if (errflg/=0) then
          write(error_unit,'(a,i0,a)') 'An error occurred in ccpp_physics_timestep_init for group time_vary: ' // trim(errmsg) // '. Exiting...'
          ierr=errflg
          return
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

!$OMP parallel num_threads (nthrds)                        &
!$OMP          default (none)                              &
!$OMP          shared (nblks, ccpp_suite,                  &
!$OMP                  step, GFS_Control, GFS_Interstitial,&
!$OMP                  dycore, nthrds)                     &
!$OMP          private (nb, nt, errmsg, errflg)            &
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
                call ccpp_physics_run(ccpp_suite=trim(ccpp_suite), group_name="phys_ps", &
                     errmsg=errmsg, errflg=errflg, lb=GFS_control%chunk_begin(nb), ub=GFS_control%chunk_end(nb), &
                     mythread=nt, nthreads=nthrds, nphys_threads=1)
                if (errflg/=0) then
                   write(error_unit,'(a,i0,a)') 'An error occurred in ccpp_physics_run for group phys_ps: ' // trim(errmsg) // '. Exiting...'
                   ierr = ierr + errflg
                end if
                ! Time-split physics
                call ccpp_physics_run(ccpp_suite=trim(ccpp_suite), group_name="phys_ts", &
                     errmsg=errmsg, errflg=errflg, lb=GFS_control%chunk_begin(nb), ub=GFS_control%chunk_end(nb), &
                     mythread=nt, nthreads=nthrds, nphys_threads=1)
                if (errflg/=0) then
                   write(error_unit,'(a,i0,a)') 'An error occurred in ccpp_physics_run for group phys_ts: ' // trim(errmsg) // '. Exiting...'
                   ierr = ierr + errflg
                end if
             endif
             if (trim(dycore)=="mpas") then
                ! Physics
                call ccpp_physics_run(ccpp_suite=trim(ccpp_suite), group_name="physics", &
                     errmsg=errmsg, errflg=errflg, lb=GFS_control%chunk_begin(nb), ub=GFS_control%chunk_end(nb), &
                     mythread=nt, nthreads=nthrds, nphys_threads=1)
                if (errflg/=0) then
                   write(error_unit,'(a,i0,a)') 'An error occurred in ccpp_physics_run for group physics: ' // trim(errmsg) // '. Exiting...'
                   ierr = ierr + errflg
                endif
             endif
          else
             ! DH* WHY WAS THIS FOR RADIATION ONLY??? PERFORMANCE?
             ! Reset GFS_Interstitial DDT fields for this thread
             call GFS_Interstitial(nt)%reset(GFS_control)
             ! *DH
             ! Radiation
             if (trim(step)=="radiation") then
                call ccpp_physics_run(ccpp_suite=trim(ccpp_suite), group_name="radiation", &
                     errmsg=errmsg, errflg=errflg, lb=GFS_control%chunk_begin(nb), ub=GFS_control%chunk_end(nb), &
                     mythread=nt, nthreads=nthrds, nphys_threads=1)
                if (errflg/=0) then
                   write(error_unit,'(a,i0,a)') 'An error occurred in ccpp_physics_run for group radiation: ' // trim(errmsg) // '. Exiting...'
                   ierr = ierr + errflg
                end if
             ! Microphysics (MPAS only)
             else if (trim(step)=="microphysics") then
                if (trim(dycore)=="mpas") then
                   call ccpp_physics_run(ccpp_suite=trim(ccpp_suite), group_name="microphysics", &
                        errmsg=errmsg, errflg=errflg, lb=GFS_control%chunk_begin(nb), ub=GFS_control%chunk_end(nb), &
                        mythread=nt, nthreads=nthrds, nphys_threads=1)
                   if (errflg/=0) then
                      write(error_unit,'(a,i0,a)') 'An error occurred in ccpp_physics_run for group microphysics: ' // trim(errmsg) // '. Exiting...'
                      ierr = ierr + errflg
                   end if
                else
                   write(error_unit,'(a)') "An error occurred in ccpp_physics_run for group microphysics. Group microphysics only valid with MPAS dycore."
                   ierr = ierr + errflg
                endif
             ! Stochastic physics
             else if (trim(step)=="stochastics") then
                call ccpp_physics_run(ccpp_suite=trim(ccpp_suite), group_name="stochastics", &
                     errmsg=errmsg, errflg=errflg, lb=GFS_control%chunk_begin(nb), ub=GFS_control%chunk_end(nb), &
                     mythread=nt, nthreads=nthrds, nphys_threads=1)
                if (errflg/=0) then
                   write(error_unit,'(a,i0,a)') 'An error occurred in ccpp_physics_run for group stochastics: ' // trim(errmsg) // '. Exiting...'
                   ierr = ierr + errflg
                end if
             ! Catchall
             else
                write(error_unit,'(a,i0,a)') 'An error occurred in ccpp_physics_run: unknown group ' // trim(step) // '. Exiting...'
                ierr = ierr + 1
             end if
          endif
          call GFS_Interstitial(nt)%destroy(GFS_control)
       end do
!$OMP end do

!$OMP end parallel
       if (ierr/=0) return

    ! Timestep final = time_vary (same for all dynamical cores)
    else if (trim(step)=="timestep_final") then

       call ccpp_physics_timestep_final(ccpp_suite=trim(ccpp_suite), group_name="all", &
            errmsg=errmsg, errflg=errflg, lb=1, ub=GFS_control%ncols,         &
            mythread=1, nthreads=1, nphys_threads=nthrds)
       if (errflg/=0) then
          write(error_unit,'(a,i0,a)') 'An error occurred in ccpp_physics_timestep_final group time_vary: ' // trim(errmsg) // '. Exiting...'
          ierr=errflg
          return
       end if

    ! Physics final (same for all dynamical cores)
    else if (trim(step)=="physics_final") then
       call ccpp_physics_final(ccpp_suite=trim(ccpp_suite), group_name='all', &
            errmsg=errmsg, errflg=errflg, lb=1, ub=GFS_control%ncols,         &
            mythread=1, nthreads=1, nphys_threads=nthrds)
       if (errflg/=0) then
          write(error_unit,'(a,i0,a)') 'An error occurred in ccpp_physics_final: ' // trim(errmsg) // '. Exiting...'
          ierr=errflg
          return
       end if
       
    ! Frameowrk final (same for all dynamical cores)
    else if (trim(step)=="final") then
       call ccpp_final(ccpp_suite=trim(ccpp_suite), errmsg=errmsg, errflg=errflg)
       if (errflg/=0) then
          write(error_unit,'(a,i0,a)') 'An error occurred in ccpp_final: ' // trim(errmsg) // '. Exiting...'
          ierr=errflg
          return
       end if
    ! Undefined CCPP step
    else
       write(error_unit,'(2a)') 'Error, undefined CCPP step ', trim(step)
       ierr=errflg
       return
    end if

  end subroutine CCPP_step

end module CCPP_driver
