module ufs_mpas_subdriver
  use mpas_derived_types, only : core_type, domain_type, mpas_Clock_type, mpas_TimeInterval_type
  use module_mpas_config, only : pio_subsystem, pio_stride, pio_numiotasks, pio_iodesc
  use module_mpas_config, only : ic_filename, lbc_filename
  use module_mpas_config, only : pio_iotype, fcst_mpi_comm, pioid
  use module_mpas_config, only : zref, zref_edge, sphere_radius, pref, pref_edge
  use module_mpas_config, only : maxNCells, maxEdges, nVertLevels
  use module_mpas_config, only : nCellsGlobal, nEdgesGlobal, nVerticesGlobal
  use module_mpas_config, only : nCellsSolve, nEdgesSolve, nVerticesSolve, nVertLevelsSolve
  use module_mpas_config, only : dt_atmos, n_atmos
  implicit none
  
  private

  public :: ufs_mpas_init_phase1
  public :: ufs_mpas_init_phase2
  public :: ufs_mpas_run
  public :: ufs_mpas_dyn_set
  public :: ufs_mpas_open_init
  public :: ufs_mpas_read_init
  public :: corelist, domain_ptr

  type(core_type),       pointer :: corelist   => null()
  type(domain_type),     pointer :: domain_ptr => null()
  type(mpas_Clock_type), pointer :: clock      => null()
  type (mpas_TimeInterval_type)  :: timeStep
  
contains
  !> #########################################################################################
  !> Procedure to initialize UWM with MPAS dynamical core.
  !>
  !> #########################################################################################
  subroutine ufs_mpas_init_phase1(Cfg, time_start, time_end, total_time, calendar, logUnits)
    ! MPAS
    use mpas_pool_routines,         only : mpas_pool_add_config, mpas_pool_get_subpool
    use mpas_pool_routines,         only : mpas_pool_add_dimension, mpas_pool_get_field
    use mpas_pool_routines,         only : mpas_pool_get_array, mpas_pool_get_config
    use mpas_framework,             only : mpas_framework_init_phase1, mpas_framework_init_phase2
    use mpas_domain_routines,       only : mpas_allocate_domain, mpas_pool_get_dimension
    use mpas_bootstrapping,         only : mpas_bootstrap_framework_phase1
    use mpas_bootstrapping,         only : mpas_bootstrap_framework_phase2
    use mpas_stream_inquiry,        only : mpas_stream_inquiry_new_streaminfo
    use mpas_derived_types,         only : mpas_pool_type, mpas_IO_NETCDF, field3dReal
    use mpas_kind_types,            only : StrKIND, RKIND
    use mpas_log,                   only : mpas_log_write
    use atm_core_interface,         only : atm_setup_core, atm_setup_domain
    use mpas_constants,             only : mpas_constants_compute_derived
    use mpas_attlist,               only : mpas_add_att
    ! FMS
    use field_manager_mod,          only : MODEL_ATMOS
    use fms2_io_mod,                only : file_exists
    use mpp_mod,                    only : FATAL, mpp_error
    ! UFSATM
    use mpas_typedefs,              only : mpas_control_type
    ! PIO
    use pio,                        only : pio_global, pio_get_att
    ! Arguments
    type(mpas_control_type), intent(inout) :: Cfg
    integer,                 intent(in   ) :: time_start(6), time_end(6), logUnits(2)
    integer,                 intent(in   ) :: total_time
    character(17),           intent(in   ) :: calendar
    ! Locals
    character(len=*), parameter :: subname = 'ufs_mpas_subdriver::ufs_mpas_init_phase1'
    integer :: i, ndate1, ndate2, tod, ierr, ik, kk
    type (mpas_pool_type), pointer :: state, mesh, tend
    type (field3dReal), pointer :: scalarsField
    character (len=StrKIND), pointer :: initial_time, config_start_time

    ! Setup MPAS infrastructure
    allocate(corelist, stat=ierr)
    if ( ierr /= 0 ) call mpp_error(FATAL,subname//": failed to allocate corelist array")
    nullify(corelist % next)

    allocate(corelist % domainlist, stat=ierr)
    if ( ierr /= 0 ) call mpp_error(FATAL,subname//": failed to allocate corelist%domainlist%next")
    nullify(corelist % domainlist % next)

    domain_ptr => corelist % domainlist
    domain_ptr % core => corelist

    call mpas_allocate_domain(domain_ptr)
    domain_ptr % domainID = 0

    ! Initialize MPAS infrastructure
    call mpas_framework_init_phase1(domain_ptr % dminfo, external_comm=fcst_mpi_comm)

    call atm_setup_core(corelist)
    call atm_setup_domain(domain_ptr)

    ! Set up the log manager as early as possible so we can use it for any errors/messages
    ! during subsequent init steps.  We need:
    ! 1) domain_ptr to be allocated,
    ! 2) dmpar_init complete to access dminfo,
    ! 3) *_setup_core to assign the setup_log function pointer
    domain_ptr % core % git_version = 'unknown'
    domain_ptr % core % build_target = 'N/A'
    ierr = domain_ptr % core % setup_log(domain_ptr % logInfo, domain_ptr, unitNumbers=logUnits)
    if ( ierr /= 0 ) then
       call mpp_error(FATAL,subname//": Log setup failed for MPAS-A dycore")
    end if

    ! Read MPAS namelist.
    if (file_exists('input.nml')) then
       call read_mpas_namelist('input.nml', domain_ptr % configs, Cfg % mpi_comm, Cfg % master, Cfg % me)
    else
       call mpp_error(FATAL,subname//": Cannot find MPAS namelist file, input.nml")
    end if

    ! Set forecast start time (config_start_time)
    ndate1 = time_start(1)*10000 + time_start(2)*100 + time_start(3)
    tod    = time_start(4)*3600  + time_start(5)*60  + time_start(6)
    call mpas_pool_add_config(domain_ptr % configs, 'config_start_time', date2yyyymmdd(ndate1)//'_'//sec2hms(tod))

    ! Set forecast end time (config_stop_time)
    ndate2 = time_end(1)*10000   + time_end(2)*100   + time_end(3)
    tod	   = time_end(4)*3600    + time_end(5)*60    + time_end(6)
    call mpas_pool_add_config(domain_ptr % configs, 'config_stop_time', date2yyyymmdd(ndate2)//'_'//sec2hms(tod))

    ! Set forecaste run time (config_run_duration) #DJS2025 this is not correct. need to fix, but works for current test.
    tod = ndate2 - ndate1 -1
    call mpas_pool_add_config(domain_ptr % configs, 'config_run_duration', trim(int2str(tod))//'_'//sec2hms(total_time))

    ! Set other MPAS required configuration information.
    call mpas_pool_add_config(domain_ptr % configs, 'config_restart_timestamp_name', 'restart_timestamp')
    call mpas_pool_add_config(domain_ptr % configs, 'config_IAU_option',             'off')
    call mpas_pool_add_config(domain_ptr % configs, 'config_do_DAcycling',           .false.)
    call mpas_pool_add_config(domain_ptr % configs, 'config_halo_exch_method',       'mpas_halo')

    ! Initialize MPAS infrastructure (phase 2)
    call mpas_framework_init_phase2(domain_ptr, io_system=pio_subsystem, calendar = trim(calendar))

    ! Before defining packages, initialize the stream inquiry instance for the domain
    domain_ptr % streamInfo => mpas_stream_inquiry_new_streaminfo()
    if (.not. associated(domain_ptr % streamInfo)) then
       call mpp_error(FATAL,subname//": Failed to instantiate streamInfo object for "//trim(domain_ptr % core % coreName))
    end if

    ierr = domain_ptr % core % define_packages(domain_ptr % packages)
    if (ierr /= 0) then
       call mpp_error(FATAL,subname//": Package definition failed for "//trim(domain_ptr % core % coreName))
    end if

    ierr = domain_ptr % core % setup_packages(domain_ptr % configs,  domain_ptr % streamInfo,       &
                                              domain_ptr % packages, domain_ptr % iocontext)
    if (ierr /= 0) then
       call mpp_error(FATAL,subname//": Package setup failed for "//trim(domain_ptr % core % coreName))
    end if

    ierr = domain_ptr % core % setup_decompositions(domain_ptr % decompositions)
    if (ierr /= 0) then
       call mpp_error(FATAL,subname//": Decomposition setup failed for "//trim(domain_ptr % core % coreName))
    end if

    ierr = domain_ptr % core % setup_clock(domain_ptr % clock, domain_ptr % configs)
    if (ierr /= 0) then
       call mpp_error(FATAL,subname//": Clock setup failed for "//trim(domain_ptr % core % coreName))
    end if

    ! Adding a config named 'cam_pcnst' with the number of constituents will indicate to
    ! MPAS-A setup code that it is operating as a UFS dycore, and that it is necessary to
    ! allocate scalars separately from other Registry-defined fields
    call mpas_pool_add_config(domain_ptr % configs, 'cam_pcnst', Cfg % nwat)

    ! Call MPAS framework bootstrap phase 1
    call mpas_bootstrap_framework_phase1(domain_ptr, "external mesh file", mpas_IO_NETCDF, pio_file_desc=pioid)

    ! Finalize the setup of blocks and fields
    call mpas_bootstrap_framework_phase2(domain_ptr, pio_file_desc=pioid)

    ! Set up scalars
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'state', state)
    call mpas_pool_get_field(state, 'scalars', scalarsField, timeLevel=1)
    !
    call mpas_pool_add_dimension(state, 'index_qv', 1)
    scalarsField % constituentNames(1) = 'qv'
    !
    call mpas_pool_add_dimension(state, 'index_qc', 2)
    scalarsField % constituentNames(2) = 'qc'
    !
    call mpas_pool_add_dimension(state, 'index_qr', 3)
    scalarsField % constituentNames(3) = 'qr'
    !
    call mpas_pool_add_dimension(state, 'index_qs', 4)
    scalarsField % constituentNames(4) = 'qs'
    !
    call mpas_pool_add_dimension(state, 'index_qi', 5)
    scalarsField % constituentNames(5) = 'qi'
    !
    call mpas_pool_add_dimension(state, 'index_qh', 6)
    scalarsField % constituentNames(6) = 'qh'
    !
    call mpas_pool_add_dimension(state, 'moist_start', 1)
    call mpas_pool_add_dimension(state, 'moist_end', Cfg % nwat)

    ! Set inital_time
    call mpas_pool_get_config(domain_ptr % blocklist % configs, 'config_start_time', config_start_time)
    call mpas_pool_get_array(state, 'initial_time', initial_time,1)
    initial_time = config_start_time
    
    ! Define scalars_tend
    nullify(tend)
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'tend', tend)
    
    if (.not. associated(tend)) then
       call mpp_error(FATAL,subname//': ERROR: The ''tend'' pool was not found.')
       ierr = 1
       return
    end if

    nullify(scalarsField)
    call mpas_pool_get_field(tend, 'scalars_tend', scalarsField, timeLevel=1)

    if (.not. associated(scalarsField)) then
       call mpp_error(FATAL,subname//': ERROR: The ''scalars_tend'' field was not found in the ''tend'' pool')
       ierr = 1
       return
    end if
    call mpas_pool_add_dimension(tend, 'index_qv', 1)
    scalarsField % constituentNames(1) = 'tend_qv'
    call mpas_add_att(scalarsField % attLists(1) % attList, 'units', 'kg m^{-3} s^{-1}')
    call mpas_add_att(scalarsField % attLists(1) % attList, 'long_name', 'Tendency of water vapor mixing ratio')
    scalarsField % constituentNames(2) = 'tend_qc'
    scalarsField % constituentNames(3) = 'tend_qh'
    scalarsField % constituentNames(4) = 'tend_qr'
    scalarsField % constituentNames(5) = 'tend_qi'
    scalarsField % constituentNames(6) = 'tend_qs'
    
    
    ! Read in static (invariant) data
    call ufs_mpas_read_invariant()

    ! Compute unit vectors giving the local north and east directions as well as
    ! the unit normal vector for edges
    call ufs_mpas_compute_unit_vectors()
    
    ! Access dimensions that are made public via this module
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh', mesh)
    call mpas_pool_get_dimension(mesh, 'nCellsSolve',    nCellsSolve)
    call mpas_pool_get_dimension(mesh, 'nEdgesSolve',    nEdgesSolve)
    call mpas_pool_get_dimension(mesh, 'nVerticesSolve', nVerticesSolve)
    call mpas_pool_get_dimension(mesh, 'nVertLevels',    nVertLevelsSolve) ! MPAS always solves over the full column

    ! Read the global sphere_radius attribute.  This is needed to normalize the cell areas.
    ierr = pio_get_att(pioid, pio_global, 'sphere_radius', sphere_radius)
    if( ierr /= 0 ) then
       call mpp_error(FATAL,subname//": Could not find sphere_radius PIO attribute")
    endif

    ! Query global grid dimensions from MPAS
    call ufs_mpas_get_global_dims(nCellsGlobal, nEdgesGlobal, nVerticesGlobal, maxEdges, nVertLevels, maxNCells)

    ! Setup constants
    call mpas_constants_compute_derived()

  end subroutine ufs_mpas_init_phase1

  !> ########################################################################################
  !> Procedure to initialize UWM with MPAS dynamical core.
  !> 
  !> ########################################################################################
  subroutine ufs_mpas_init_phase2(Cfg, Statein)
    use mpas_typedefs,              only : mpas_control_type
    use mpas_typedefs,              only : mpas_statein_type
    use mpas_kind_types,            only : StrKIND, RKIND
    use mpas_derived_types,         only : mpas_pool_type, mpas_Time_Type
    use mpas_domain_routines,       only : mpas_pool_get_dimension
    use mpas_pool_routines,         only : mpas_pool_get_subpool
    use mpas_pool_routines,         only : mpas_pool_initialize_time_levels, mpas_pool_get_config
    use mpas_pool_routines,         only : mpas_pool_get_array
    use mpas_atm_dimensions,        only : mpas_atm_set_dims
    use mpas_atm_threading,         only : mpas_atm_threading_init
    use mpp_mod,                    only : FATAL, mpp_error
    use mpas_atm_halos,             only : atm_build_halo_groups, exchange_halo_group
    use atm_core,                   only : atm_mpas_init_block, core_clock => clock
    use atm_time_integration,       only : mpas_atm_dynamics_init
    use mpas_timekeeping,           only : mpas_get_clock_time, mpas_get_time, mpas_START_TIME
    use mpas_timekeeping,           only : mpas_set_timeInterval
    ! Arguments
    type(mpas_control_type), intent(inout) :: Cfg
    type(mpas_statein_type), intent(inout) :: Statein
    type(mpas_pool_type), pointer :: tend_physics_pool
    ! Locals
    character(len=*), parameter :: subname = 'ufs_mpas_subdriver::ufs_mpas_init_phase2'
    type (mpas_pool_type), pointer :: state, mesh
    integer :: ierr
    integer, pointer :: nVertLevels1, maxEdges1, maxEdges2, num_scalars
    real (kind=RKIND), pointer :: dt
    logical, pointer :: config_do_restart
    type (mpas_Time_Type) :: startTime
    character(len=StrKIND) :: startTimeStamp
    character (len=StrKIND), pointer :: xtime
    character (len=StrKIND), pointer :: initial_time1, initial_time2

    !
    ! Setup threading
    !
    call mpas_atm_threading_init(domain_ptr%blocklist, ierr)
    if ( ierr /= 0 ) then
       call mpp_error(FATAL,subname//": Threading setup failed for core "//trim(domain_ptr % core % coreName))
    end if

    !
    ! Set up inner dimensions used by arrays in optimized dynamics routines
    !
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'state', state)
    call mpas_pool_get_dimension(state, 'nVertLevels', nVertLevels1)
    call mpas_pool_get_dimension(state, 'maxEdges', maxEdges1)
    call mpas_pool_get_dimension(state, 'maxEdges2', maxEdges2)
    call mpas_pool_get_dimension(state, 'num_scalars', num_scalars)
    call mpas_atm_set_dims(nVertLevels1, maxEdges1, maxEdges2, num_scalars)
    Cfg % levs = nVertLevels1

    !
    ! Set "local" clock to point to the clock contained in the domain type
    !
    clock => domain_ptr % clock
    core_clock => domain_ptr % clock

    !
    ! Build halo exchange groups and set method for exchanging halos in a group
    !
    call atm_build_halo_groups(domain_ptr, ierr)
    if (ierr /= 0) then
       call mpp_error(FATAL,subname//": failed to build MPAS-A halo exchange groups.")
    end if

    call mpas_pool_get_config(domain_ptr % blocklist % configs, 'config_do_restart', config_do_restart)
    call mpas_pool_get_config(domain_ptr % blocklist % configs, 'config_dt', dt)

    ! How many calls to MPAS dycore for each ATMosphere time step?
    Cfg%dt_dycore = dt
    n_atmos = dt_atmos/dt

    if (.not. config_do_restart) then
       call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'state', state)
       call mpas_pool_initialize_time_levels(state)
    end if

    !
    ! Set startTimeStamp based on the start time of the simulation clock
    !
    startTime = mpas_get_clock_time(clock, mpas_START_TIME, ierr)
    if ( ierr /= 0 ) then
       call mpp_error(FATAL,subname//": failed to get mpas_START_TIME")
    end if
    call mpas_get_time(startTime, dateTimeString=startTimeStamp)

    call exchange_halo_group(domain_ptr, 'initialization:u')

    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh', mesh)
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'state', state)

    call atm_mpas_init_block(domain_ptr % dminfo, domain_ptr % streamManager, domain_ptr % blocklist, mesh, dt)

    call mpas_pool_get_array(state, 'xtime', xtime, 1)
    xtime = startTimeStamp

    ! Initialize initial_time in second time level. We need to do this because initial state
    ! is read into time level 1, and if we write output from the set of state arrays that
    ! represent the original time level 2, the initial_time field will be invalid.

    call mpas_pool_get_array(state, 'initial_time', initial_time1, 1)
    call mpas_pool_get_array(state, 'initial_time', initial_time2, 2)
    initial_time2 = initial_time1
    
    call exchange_halo_group(domain_ptr, 'initialization:pv_edge,ru,rw')

    !
    ! Prepare the dynamics for integration
    !
    call mpas_atm_dynamics_init(domain_ptr)

    !
    ! Initialize physics tendencies needed by MPAS.
    !
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'tend_physics', tend_physics_pool)
    call mpas_pool_get_array(tend_physics_pool, 'tend_ru_physics',     Statein % ru_tend)
    call mpas_pool_get_array(tend_physics_pool, 'tend_rtheta_physics', Statein % rtheta_tend)
    call mpas_pool_get_array(tend_physics_pool, 'tend_rho_physics',    Statein % rho_tend)

    ! Set dycore time interval.
    call mpas_set_timeInterval(timeStep, dt=dt, ierr=ierr)

  end subroutine ufs_mpas_init_phase2

  !> #########################################################################################
  !> Routine to call MPAS dynamical core
  !> Loop over dynamical time-step(s) and increment MPAS state (timelevel 1->2)
  !>
  !> #########################################################################################
  subroutine ufs_mpas_run(statein, stateout)
    ! UFSATM
    use mpas_typedefs,        only : mpas_statein_type, mpas_stateout_type
    ! MPAS
    use atm_core,             only : atm_do_timestep
    use mpas_domain_routines, only : mpas_pool_get_dimension
    use mpas_derived_types,   only : mpas_Time_type, mpas_pool_type
    use mpas_timekeeping,     only : mpas_set_timeInterval
    use mpas_kind_types,      only : StrKIND, RKIND, R8KIND
    use mpas_constants,       only : rvord
    use mpas_pool_routines,   only : mpas_pool_get_config, mpas_pool_get_subpool
    use mpas_pool_routines,   only : mpas_pool_shift_time_levels, mpas_pool_get_array
    use mpas_timekeeping,     only : mpas_advance_clock, mpas_get_clock_time, mpas_get_time
    use mpas_timekeeping,     only : mpas_NOW, mpas_is_clock_stop_time, mpas_dmpar_get_time
    use mpas_log,             only : mpas_log_write
    use mpas_timer,           only : mpas_timer_start, mpas_timer_stop
    ! Arguments
    type(mpas_statein_type),  intent(inout) :: statein
    type(mpas_stateout_type), intent(inout) :: stateout
    ! Locals
    real (kind=RKIND), pointer :: dt
    type (mpas_pool_type), pointer :: state, diag, mesh
    type (mpas_Time_type) :: timeNow, timeStop
    character(len=StrKIND) :: timeStamp
    integer :: ierr, itime
    integer, pointer :: index_qv
    integer, pointer :: nCellsSolve
    real(kind=RKIND), dimension(:,:), pointer :: theta_m, rho_zz, zz, theta, rho
    real(kind=RKIND), dimension(:,:,:), pointer :: scalars
    integer :: itimestep
    real (kind=R8KIND) :: integ_start_time, integ_stop_time 
    integer, parameter :: id = 40
    logical, pointer :: config_apply_lbcs
    
    clock => domain_ptr % clock

    ! Eventually, dt should be domain specific
    call mpas_pool_get_config( domain_ptr % blocklist % configs, 'config_dt', dt)

    !
    ! Read initial boundary state
    ! NOT YET IMPLEMENTED (Follow src/core_atmosphere/mpas_atm_core.F:atm_core_run())
    call mpas_pool_get_config( domain_ptr % blocklist % configs, 'config_apply_lbcs', config_apply_lbcs)
    if (config_apply_lbcs) then
       
    endif
    
    ! During integration, time level 1 stores the model state at the beginning of the
    !   time step, and time level 2 stores the state advanced dt in time by timestep(...)
    ! ONLY RUNNING SINGLE TIMESTEP.
    itimestep=1
    do itime = 1, 1!n_atmos
       ! Get current time.
       timeNow  = mpas_get_clock_time(clock, mpas_NOW, ierr)
       call mpas_get_time(curr_time=timeNow, dateTimeString=timeStamp, ierr=ierr)
       call mpas_log_write('') 
       call mpas_log_write(' MPAS dynamics start timestep '//trim(timeStamp))

       ! Integrate forward one dycore time step
       call mpas_timer_start('time integration')
       call mpas_dmpar_get_time(integ_start_time)
       call atm_do_timestep(domain_ptr, dt, itimestep)
       call mpas_dmpar_get_time(integ_stop_time)
       call mpas_timer_stop('time integration')
       call mpas_log_write(' Timing for integration step: $r s', realArgs=(/real(integ_stop_time - integ_start_time, kind=RKIND)/))
       
       ! Move time level 2 fields back into time level 1 for next time step
       call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'state', state)
       call mpas_pool_shift_time_levels(state)

       ! Advance clock.
       itimestep = itimestep + 1
       call mpas_advance_clock(clock)
       timeNow = mpas_get_clock_time(clock, mpas_NOW, ierr)

       ! Print IN/OUT state (DEBUGGING)
       call mpas_pool_get_array(state, 'u',       stateout % uperp,   timeLevel=1)
       call mpas_pool_get_array(state, 'w',       stateout % w,       timeLevel=1)
       call mpas_pool_get_array(state, 'theta_m', stateout % theta_m, timeLevel=1)
       call mpas_pool_get_array(state, 'rho_zz',  stateout % rho_zz,  timeLevel=1)
       call mpas_pool_get_array(state, 'scalars', stateout % tracers, timeLevel=1)
       print*,'u(IN,OUT):       ',statein % uperp(1,id),     stateout % uperp(1,id)
       print*,'w(IN,OUT):       ',statein % w(1,id),         stateout % w(1,id)
       print*,'theta_m(IN,OUT): ',statein % theta_m(1,id),   stateout % theta_m(1,id)
       print*,'rho_zz(IN,OUT):  ',statein % rho_zz(1,id),    stateout % rho_zz(1,id)
       print*,'scalars(IN,OUT): ',statein % tracers(1,1,id), stateout % tracers(1,1,id)
    end do

    !
    ! Update final prognostic state
    !
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'state', state)
    call mpas_pool_get_array(state, 'u',       stateout % uperp,   timeLevel=1)
    call mpas_pool_get_array(state, 'w',       stateout % w,       timeLevel=1)
    call mpas_pool_get_array(state, 'theta_m', stateout % theta_m, timeLevel=1)
    call mpas_pool_get_array(state, 'rho_zz',  stateout % rho_zz,  timeLevel=1)
    call mpas_pool_get_array(state, 'scalars', stateout % tracers, timeLevel=1)
    print*,'#####################################################################'
    print*,'#u(IN,OUT):       ',statein % uperp(1,id),     stateout % uperp(1,id)
    print*,'#w(IN,OUT):       ',statein % w(1,id),         stateout % w(1,id)
    print*,'#theta_m(IN,OUT): ',statein % theta_m(1,id),   stateout % theta_m(1,id)
    print*,'#rho_zz(IN,OUT):  ',statein % rho_zz(1,id),    stateout % rho_zz(1,id)
    print*,'#scalars(IN,OUT): ',statein % tracers(1,1,id), stateout % tracers(1,1,id)

    !
    ! Compute diagnostic fields from the final prognostic state
    ! From mpas_atm_core.F:_atm_compute_output_diagnostics()
    !
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'diag', diag)
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh', mesh)
    call mpas_pool_get_dimension(state, 'index_qv', index_qv)
    call mpas_pool_get_dimension(state, 'nCellsSolve', nCellsSolve)
    call mpas_pool_get_array(state, 'theta_m', theta_m, timeLevel=1)
    call mpas_pool_get_array(state, 'rho_zz',  rho_zz,  timeLevel=1)
    call mpas_pool_get_array(state, 'scalars', scalars, timeLevel=1)
    call mpas_pool_get_array(mesh,  'zz',      zz)
    call mpas_pool_get_array(diag,  'theta',   theta)
    call mpas_pool_get_array(diag,  'rho',     rho)

    rho(:,1:nCellsSolve) = rho_zz(:,1:nCellsSolve) * zz(:,1:nCellsSolve)
    theta(:,1:nCellsSolve) = theta_m(:,1:nCellsSolve) / (1.0_RKIND + rvord * scalars(index_qv,:,1:nCellsSolve))
    
  end subroutine ufs_mpas_run
  
  !> #########################################################################################
  !> Procedure to populate input/output state for MPAS dycore.
  !>
  !> #########################################################################################
  subroutine ufs_mpas_dyn_set(statein, stateout)
    ! UFSATM
    use mpas_typedefs,        only : mpas_statein_type, mpas_stateout_type
    ! MPAS
    use mpas_derived_types,   only : mpas_pool_type
    use mpas_pool_routines,   only : mpas_pool_get_subpool
    use mpas_pool_routines,   only : mpas_pool_get_array
    use mpas_domain_routines, only : mpas_pool_get_dimension
    use mpas_kind_types,      only : RKIND
    ! FMS
    use mpp_mod,              only : FATAL, mpp_error
    ! Arguments
    type(mpas_statein_type), intent(inout) :: statein
    type(mpas_stateout_type), intent(inout) :: stateout
    ! Locals
    character(len=*), parameter :: subname = 'ufs_mpas_subdriver::ufs_mpas_dyn_set'
    type(mpas_pool_type), pointer :: mesh_pool
    type(mpas_pool_type), pointer :: state_pool
    type(mpas_pool_type), pointer :: diag_pool
    integer, pointer :: nCells, nEdges, nVertices, nVertLevels, nCellsSolve, nEdgesSolve, &
         nVerticesSolve, index_qv, index_qc, index_qh, index_qr, index_qi, index_qs
    integer :: i1, i2, ierr
    real(kind=RKIND),dimension(:,:),pointer  :: rtheta_p

    !
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh',         mesh_pool)
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'state',        state_pool)
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'diag',         diag_pool)

    ! Get dimensions
    call mpas_pool_get_dimension(mesh_pool,  'nCells',         nCells)
    call mpas_pool_get_dimension(mesh_pool,  'nEdges',         nEdges)
    call mpas_pool_get_dimension(mesh_pool,  'nVertices',      nVertices)
    call mpas_pool_get_dimension(mesh_pool,  'nVertLevels',    nVertLevels)
    call mpas_pool_get_dimension(mesh_pool,  'nCellsSolve',    nCellsSolve)
    call mpas_pool_get_dimension(mesh_pool,  'nEdgesSolve',    nEdgesSolve)
    call mpas_pool_get_dimension(mesh_pool,  'nVerticesSolve', nVerticesSolve)
    call mpas_pool_get_dimension(state_pool, 'index_qv',       index_qv)
    call mpas_pool_get_dimension(state_pool, 'index_qc',       index_qc)
    call mpas_pool_get_dimension(state_pool, 'index_qh',       index_qh)
    call mpas_pool_get_dimension(state_pool, 'index_qr',       index_qr)
    call mpas_pool_get_dimension(state_pool, 'index_qi',       index_qi)
    call mpas_pool_get_dimension(state_pool, 'index_qs',       index_qs)

    ! Set dimensions
    statein % nCells         = nCells
    statein % nEdges         = nEdges
    statein % nVertices      = nVertices
    statein % nVertLevels    = nVertLevels
    statein % nCellsSolve    = nCellsSolve
    statein % nEdgesSolve    = nEdgesSolve
    statein % nVerticesSolve = nVerticesSolve
    statein % index_qv       = index_qv
    statein % index_qc       = index_qc
    statein % index_qh       = index_qh
    statein % index_qr       = index_qr
    statein % index_qi       = index_qi
    statein % index_qs       = index_qs

    ! In MPAS timeLevel=1 is the current state.  So the fields input to the dycore should
    ! be in timeLevel=1.
    call mpas_pool_get_array(state_pool, 'u',                      statein % uperp,   timeLevel=1)
    call mpas_pool_get_array(state_pool, 'w',                      statein % w,       timeLevel=1)
    call mpas_pool_get_array(state_pool, 'theta_m',                statein % theta_m, timeLevel=1)
    call mpas_pool_get_array(state_pool, 'rho_zz',                 statein % rho_zz,  timeLevel=1)
    call mpas_pool_get_array(state_pool, 'scalars',                statein % tracers, timeLevel=1)
    !
    call mpas_pool_get_array(diag_pool, 'rho_base',                statein % rho_base)
    call mpas_pool_get_array(diag_pool, 'theta_base',              statein % theta_base)
    !
    call mpas_pool_get_array(mesh_pool,  'zgrid',                  statein % zint)
    call mpas_pool_get_array(mesh_pool,  'zz',                     statein % zz)
    call mpas_pool_get_array(mesh_pool,  'fzm',                    statein % fzm)
    call mpas_pool_get_array(mesh_pool,  'fzp',                    statein % fzp)
    call mpas_pool_get_array(mesh_pool,  'areaCell',               statein % areaCell)
    !
    call mpas_pool_get_array(mesh_pool,  'east',                   statein % east)
    call mpas_pool_get_array(mesh_pool,  'north',                  statein % north)
    call mpas_pool_get_array(mesh_pool,  'edgeNormalVectors',      statein % normal)
    call mpas_pool_get_array(mesh_pool,  'cellsOnEdge',            statein % cellsOnEdge)
    !
    call mpas_pool_get_array(diag_pool,  'theta',                  statein % theta)
    call mpas_pool_get_array(diag_pool,  'exner',                  statein % exner)
    call mpas_pool_get_array(diag_pool,  'rho',                    statein % rho)
    call mpas_pool_get_array(diag_pool,  'uReconstructZonal',      statein % ux)
    call mpas_pool_get_array(diag_pool,  'uReconstructMeridional', statein % uy)

    ! Compute variables needed in the MPAS dynamical core.
    
    ! density-weighted perturbation potential temperature:
    !call mpas_pool_get_array(diag_pool,  'rtheta_p',               rtheta_p)
    !rtheta_p = statein % rho_zz * statein % theta_m - (statein % rho_base * statein % theta_base)
    
    ! Let dynamics export state point to memory managed by MPAS-Atmosphere
    ! Exception: pmiddry and pintdry are not managed by the MPAS infrastructure
    stateout % nCells         = statein % nCells
    stateout % nEdges         = statein % nEdges
    stateout % nVertices      = statein % nVertices
    stateout % nVertLevels    = statein % nVertLevels
    stateout % nCellsSolve    = statein % nCellsSolve
    stateout % nEdgesSolve    = statein % nEdgesSolve
    stateout % nVerticesSolve = statein % nVerticesSolve
    stateout % index_qv       = statein % index_qv
    stateout % index_qc       = statein % index_qc
    stateout % index_qh       = statein % index_qh
    stateout % index_qr       = statein % index_qr
    stateout % index_qi       = statein % index_qi
    stateout % index_qs       = statein % index_qs

    ! MPAS swaps pointers internally so that after a dycore timestep, the updated state is
    ! in timeLevel=1.  Thus we want stateout to also point to timeLevel=1.  Can just copy
    ! the pointers from statein.
    stateout % uperp   => statein % uperp
    stateout % w       => statein % w
    stateout % theta_m => statein % theta_m
    stateout % rho_zz  => statein % rho_zz
    stateout % tracers => statein % tracers

    ! These components don't have a time level index.
    stateout % zint  => statein % zint
    stateout % zz    => statein % zz
    stateout % fzm   => statein % fzm
    stateout % fzp   => statein % fzp
    
    stateout % theta => statein % theta
    stateout % exner => statein % exner
    stateout % rho   => statein % rho
    stateout % ux    => statein % ux
    stateout % uy    => statein % uy

    ! Hydrostatic pressure
    allocate(stateout % pmiddry(stateout % nVertLevels,   stateout % nCells), stat=ierr)
    if( ierr /= 0 ) call mpp_error(FATAL,subname//': failed to allocate stateout%pmiddry array')

    allocate(stateout % pintdry(stateout % nVertLevels+1, stateout % nCells), stat=ierr)
    if( ierr /= 0 ) call mpp_error(FATAL,subname//': failed to allocate stateout%pintdry array')

    call mpas_pool_get_array(diag_pool, 'vorticity',  stateout % vorticity)
    call mpas_pool_get_array(diag_pool, 'divergence', stateout % divergence)

  end subroutine ufs_mpas_dyn_set

  !> #########################################################################################
  !> Procedure to open MPAS IC file.
  !>
  !> ######################################################################################### 
  subroutine ufs_mpas_open_init()
    ! PIO
    use pio,         only : pio_openfile, pio_nowrite
    ! FMS
    use fms2_io_mod, only : file_exists
    use mpp_mod,     only : FATAL, mpp_error
    ! Arguments
    ! Locals
    integer :: ierr
    character(len=*), parameter :: subname = 'ufs_mpas_subdriver::ufs_mpas_open_init'

    ! Open MPAS Initial Condition file.
    if (file_exists(ic_filename)) then
       ierr = pio_openfile(pio_subsystem, pioid, pio_iotype, ic_filename, pio_nowrite)
       if (ierr /= 0) then
          call mpp_error(FATAL,subname//": Failed opening MPAS IC File, "//trim(ic_filename))
       end if
    else
       call mpp_error(FATAL,subname//": Cannot find MPAS IC file: "//trim(ic_filename))
    end if
  end subroutine ufs_mpas_open_init
  
  !> #########################################################################################
  !> Procedure to read MPAS namelist(s).
  !>
  !> The namelist for MPAS are described in MPAS-Model/src/core_atmosphere/Registry.xml, this
  !> is also where the default values defined below originate.
  !>
  !> #########################################################################################
  subroutine read_mpas_namelist(nml_file, configPool, mpicomm, master, me)
    use mpi_f08,            only: MPI_Comm, MPI_CHARACTER, MPI_INTEGER, MPI_REAL8,  MPI_LOGICAL
    use mpi_f08,            only: mpi_bcast, mpi_barrier
    use mpas_derived_types, only: mpas_pool_type
    use mpas_kind_types,    only: StrKIND, RKIND
    use mpas_pool_routines, only: mpas_pool_add_config
    use mpas_typedefs,      only: r8 => kind_dbl_prec
    use fms_mod,            only: check_nml_error
    use mpp_mod,            only: input_nml_file
    ! Inputs
    type(MPI_Comm),       intent(in   ) :: mpicomm
    integer,              intent(in   ) :: master, me
    character(len=*),     intent(in   ) :: nml_file
    type(mpas_pool_type), intent(inout) :: configPool

    ! Namelist nhyd_model
    character (len=StrKIND) :: mpas_time_integration               = 'SRK3'
    integer                 :: mpas_time_integration_order         = 2
    real(r8)                :: mpas_dt                             = 720.0_r8
    logical                 :: mpas_split_dynamics_transport       = .true.
    integer                 :: mpas_number_of_sub_steps            = 2
    integer                 :: mpas_dynamics_split_steps           = 3
    real(r8)                :: mpas_h_mom_eddy_visc2               = 0.0_r8
    real(r8)                :: mpas_h_mom_eddy_visc4               = 0.0_r8
    real(r8)                :: mpas_v_mom_eddy_visc2               = 0.0_r8
    real(r8)                :: mpas_h_theta_eddy_visc2             = 0.0_r8
    real(r8)                :: mpas_h_theta_eddy_visc4             = 0.0_r8
    real(r8)                :: mpas_v_theta_eddy_visc2             = 0.0_r8
    character (len=StrKIND) :: mpas_horiz_mixing                   = '2d_smagorinsky'
    real(r8)                :: mpas_len_disp                       = 120000.0_r8
    real(r8)                :: mpas_visc4_2dsmag                   = 0.05_r8
    real(r8)                :: mpas_del4u_div_factor               = 10.0_r8
    integer                 :: mpas_w_adv_order                    = 3
    integer                 :: mpas_theta_adv_order                = 3
    integer                 :: mpas_scalar_adv_order               = 3
    integer                 :: mpas_u_vadv_order                   = 3
    integer                 :: mpas_w_vadv_order                   = 3
    integer                 :: mpas_theta_vadv_order               = 3
    integer                 :: mpas_scalar_vadv_order              = 3
    logical                 :: mpas_scalar_advection               = .true.
    logical                 :: mpas_positive_definite              = .false.
    logical                 :: mpas_monotonic                      = .true.
    real(r8)                :: mpas_coef_3rd_order                 = 0.25_r8
    real(r8)                :: mpas_smagorinsky_coef               = 0.125_r8
    logical                 :: mpas_mix_full                       = .true.
    real(r8)                :: mpas_epssm                          = 0.1_r8
    real(r8)                :: mpas_smdiv                          = 0.1_r8
    real(r8)                :: mpas_apvm_upwinding                 = 0.5_r8
    logical                 :: mpas_h_ScaleWithMesh                = .true.
    ! Namelist damping
    real(r8)                :: mpas_zd                             = 22000.0_r8
    real(r8)                :: mpas_xnutr                          = 0.2_r8
    real(r8)                :: mpas_cam_coef                       = 0.0_r8
    integer                 :: mpas_cam_damping_levels             = 0
    logical                 :: mpas_rayleigh_damp_u                = .true.
    real(r8)                :: mpas_rayleigh_damp_u_timescale_days = 5.0_r8
    integer                 :: mpas_number_rayleigh_damp_u_levels  = 3
    ! Namelist limited_area
    logical                 :: mpas_apply_lbcs                     = .false.
    ! Namelist PIO
    integer                 :: mpas_pio_num_iotasks                = 1
    integer                 :: mpas_pio_stride                     = 1
    ! Namelist assimilation
    logical                 :: mpas_jedi_da                        = .false.
    ! Namelist decomposition
    character (len=StrKIND) :: mpas_block_decomp_file_prefix       = 'x1.40962.graph.info.part.'
    ! Namelist restart
    logical                 :: mpas_do_restart                     = .false.
    ! Namelist printout
    logical                 :: mpas_print_global_minmax_vel        = .true.
    logical                 :: mpas_print_detailed_minmax_vel      = .true.
    logical                 :: mpas_print_global_minmax_sca        = .true.

    namelist /mpas_nhyd_model/ mpas_time_integration, mpas_time_integration_order, mpas_dt,   &
         mpas_split_dynamics_transport, mpas_number_of_sub_steps, mpas_dynamics_split_steps,  &
         mpas_h_mom_eddy_visc2, mpas_h_mom_eddy_visc4, mpas_v_mom_eddy_visc2,                 &
         mpas_h_theta_eddy_visc2, mpas_h_theta_eddy_visc4, mpas_v_theta_eddy_visc2,           &
         mpas_horiz_mixing, mpas_len_disp, mpas_visc4_2dsmag, mpas_del4u_div_factor,          &
         mpas_w_adv_order, mpas_theta_adv_order, mpas_scalar_adv_order, mpas_u_vadv_order,    & 
         mpas_w_vadv_order, mpas_theta_vadv_order, mpas_scalar_vadv_order,                    &
         mpas_scalar_advection, mpas_positive_definite, mpas_monotonic, mpas_coef_3rd_order,  &
         mpas_smagorinsky_coef, mpas_mix_full, mpas_epssm, mpas_smdiv, mpas_apvm_upwinding,   &
         mpas_h_ScaleWithMesh
    !
    namelist /mpas_damping/ mpas_zd, mpas_xnutr, mpas_cam_coef, mpas_cam_damping_levels,      &
         mpas_rayleigh_damp_u, mpas_rayleigh_damp_u_timescale_days,                           &
         mpas_number_rayleigh_damp_u_levels
    !
    namelist /mpas_limited_area/  mpas_apply_lbcs
    !
    namelist /mpas_io/ mpas_pio_num_iotasks, mpas_pio_stride
    !
    namelist /mpas_assimilation/ mpas_jedi_da
    !
    namelist /mpas_decomposition/ mpas_block_decomp_file_prefix
    !
    namelist /mpas_restart/ mpas_do_restart
    !
    namelist /mpas_printout/ mpas_print_global_minmax_vel, mpas_print_detailed_minmax_vel,    &
         mpas_print_global_minmax_sca

    ! These configuration parameters must be set in the MPAS configPool, but can't be changed
    ! in UFS. *From CAM src/dynamics/mpas/dyn_comp.F90*
    integer                :: config_num_halos = 2
    integer                :: config_number_of_blocks = 0
    logical                :: config_explicit_proc_decomp = .false.
    character(len=StrKIND) :: config_proc_decomp_file_prefix = 'graph.info.part'

    ! Locals
    integer :: ierr, io, mpierr

    ! Read in namelists...
    if (me == master) then
       print*,'Reading MPAS-A dynamical core namelist'
       ! nhyd_model
       read(input_nml_file, nml=mpas_nhyd_model, iostat=io)
       ierr = check_nml_error(io, 'mpas_nhyd_model')
       ! damping
       read(input_nml_file, nml=mpas_damping, iostat=io)
       ierr = check_nml_error(io, 'mpas_damping')
       ! limited_area
       read(input_nml_file, nml=mpas_limited_area, iostat=io)
       ierr = check_nml_error(io, 'mpas_limited_area')
       ! PIO
       read(input_nml_file, nml=mpas_io, iostat=io)
       ierr = check_nml_error(io, 'mpas_io')
       ! assimilation
       read(input_nml_file, nml=mpas_assimilation, iostat=io)
       ierr = check_nml_error(io, 'mpas_assimilation')
       ! decomposition
       read(input_nml_file, nml=mpas_decomposition, iostat=io)
       ierr = check_nml_error(io, 'mpas_decomposition')
       ! restart
       read(input_nml_file, nml=mpas_restart, iostat=io)
       ierr = check_nml_error(io, 'mpas_restart')
       ! printout
       read(input_nml_file, nml=mpas_printout, iostat=io)
       ierr = check_nml_error(io, 'mpas_printout')
    endif

    ! Other processors waiting...
    call mpi_barrier(mpicomm, mpierr)

    !
    ! MPI Broadcast to all
    !    
    call mpi_bcast(mpas_time_integration,         StrKIND, mpi_character, master, mpicomm, mpierr)
    call mpi_bcast(mpas_time_integration_order,         1, mpi_integer,   master, mpicomm, mpierr)
    call mpi_bcast(mpas_dt,                             1, mpi_real8,     master, mpicomm, mpierr)
    call mpi_bcast(mpas_split_dynamics_transport,       1, mpi_logical,   master, mpicomm, mpierr)
    call mpi_bcast(mpas_number_of_sub_steps,            1, mpi_integer,   master, mpicomm, mpierr)
    call mpi_bcast(mpas_dynamics_split_steps,           1, mpi_integer,   master, mpicomm, mpierr)
    call mpi_bcast(mpas_h_mom_eddy_visc2,               1, mpi_real8,     master, mpicomm, mpierr)
    call mpi_bcast(mpas_h_mom_eddy_visc4,               1, mpi_real8,     master, mpicomm, mpierr)
    call mpi_bcast(mpas_v_mom_eddy_visc2,               1, mpi_real8,     master, mpicomm, mpierr)
    call mpi_bcast(mpas_h_theta_eddy_visc2,             1, mpi_real8,     master, mpicomm, mpierr)
    call mpi_bcast(mpas_h_theta_eddy_visc4,             1, mpi_real8,     master, mpicomm, mpierr)
    call mpi_bcast(mpas_v_theta_eddy_visc2,             1, mpi_real8,     master, mpicomm, mpierr)
    call mpi_bcast(mpas_horiz_mixing,             StrKIND, mpi_character, master, mpicomm, mpierr)
    call mpi_bcast(mpas_len_disp,                       1, mpi_real8,     master, mpicomm, mpierr)
    call mpi_bcast(mpas_visc4_2dsmag,                   1, mpi_real8,     master, mpicomm, mpierr)
    call mpi_bcast(mpas_del4u_div_factor,               1, mpi_real8,     master, mpicomm, mpierr)
    call mpi_bcast(mpas_w_adv_order,                    1, mpi_integer,   master, mpicomm, mpierr)
    call mpi_bcast(mpas_theta_adv_order,                1, mpi_integer,   master, mpicomm, mpierr)
    call mpi_bcast(mpas_scalar_adv_order,               1, mpi_integer,   master, mpicomm, mpierr)
    call mpi_bcast(mpas_u_vadv_order,                   1, mpi_integer,   master, mpicomm, mpierr)
    call mpi_bcast(mpas_w_vadv_order,                   1, mpi_integer,   master, mpicomm, mpierr)
    call mpi_bcast(mpas_theta_vadv_order,               1, mpi_integer,   master, mpicomm, mpierr)
    call mpi_bcast(mpas_scalar_vadv_order,              1, mpi_integer,   master, mpicomm, mpierr)
    call mpi_bcast(mpas_scalar_advection,               1, mpi_logical,   master, mpicomm, mpierr)
    call mpi_bcast(mpas_positive_definite,              1, mpi_logical,   master, mpicomm, mpierr)
    call mpi_bcast(mpas_monotonic,                      1, mpi_logical,   master, mpicomm, mpierr)
    call mpi_bcast(mpas_coef_3rd_order,                 1, mpi_real8,     master, mpicomm, mpierr)
    call mpi_bcast(mpas_smagorinsky_coef,               1, mpi_real8,     master, mpicomm, mpierr)
    call mpi_bcast(mpas_mix_full,                       1, mpi_logical,   master, mpicomm, mpierr)
    call mpi_bcast(mpas_epssm,                          1, mpi_real8,     master, mpicomm, mpierr)
    call mpi_bcast(mpas_smdiv,                          1, mpi_real8,     master, mpicomm, mpierr)
    call mpi_bcast(mpas_apvm_upwinding,                 1, mpi_real8,     master, mpicomm, mpierr)
    call mpi_bcast(mpas_h_ScaleWithMesh,                1, mpi_logical,   master, mpicomm, mpierr)
    !
    call mpi_bcast(mpas_zd,                             1, mpi_real8,     master, mpicomm, mpierr)
    call mpi_bcast(mpas_xnutr,                          1, mpi_real8,     master, mpicomm, mpierr)
    call mpi_bcast(mpas_cam_coef,                       1, mpi_real8,     master, mpicomm, mpierr)
    call mpi_bcast(mpas_cam_damping_levels,             1, mpi_integer,   master, mpicomm, mpierr)
    call mpi_bcast(mpas_rayleigh_damp_u,                1, mpi_logical,   master, mpicomm, mpierr)
    call mpi_bcast(mpas_rayleigh_damp_u_timescale_days, 1, mpi_real8,     master, mpicomm, mpierr)
    call mpi_bcast(mpas_number_rayleigh_damp_u_levels,  1, mpi_integer,   master, mpicomm, mpierr)
    !
    call mpi_bcast(mpas_apply_lbcs,                     1, mpi_logical,   master, mpicomm, mpierr)
    !
    call mpi_bcast(mpas_pio_num_iotasks,                1, mpi_integer,   master, mpicomm, mpierr)
    call mpi_bcast(mpas_pio_stride,                     1, mpi_integer,   master, mpicomm, mpierr)
    !
    call mpi_bcast(mpas_jedi_da,                        1, mpi_logical,   master, mpicomm, mpierr)
    !
    call mpi_bcast(mpas_block_decomp_file_prefix, StrKIND, mpi_character, master, mpicomm, mpierr)
    !
    call mpi_bcast(mpas_do_restart,                     1, mpi_logical,   master, mpicomm, mpierr)
    !
    call mpi_bcast(mpas_print_global_minmax_vel,        1, mpi_logical,   master, mpicomm, mpierr)
    call mpi_bcast(mpas_print_detailed_minmax_vel,      1, mpi_logical,   master, mpicomm, mpierr)
    call mpi_bcast(mpas_print_global_minmax_sca,        1, mpi_logical,   master, mpicomm, mpierr)
    
    !
    ! Set MPAS configuration information pool variables
    !
    call mpas_pool_add_config(configPool, 'config_time_integration',               mpas_time_integration)
    call mpas_pool_add_config(configPool, 'config_time_integration_order',         mpas_time_integration_order)
    call mpas_pool_add_config(configPool, 'config_dt',                             real(mpas_dt,kind=RKIND))
    call mpas_pool_add_config(configPool, 'config_split_dynamics_transport',       mpas_split_dynamics_transport)
    call mpas_pool_add_config(configPool, 'config_number_of_sub_steps',            mpas_number_of_sub_steps)
    call mpas_pool_add_config(configPool, 'config_dynamics_split_steps',           mpas_dynamics_split_steps)
    call mpas_pool_add_config(configPool, 'config_h_mom_eddy_visc2',               real(mpas_h_mom_eddy_visc2,kind=RKIND))
    call mpas_pool_add_config(configPool, 'config_h_mom_eddy_visc4',               real(mpas_h_mom_eddy_visc4,kind=RKIND))
    call mpas_pool_add_config(configPool, 'config_v_mom_eddy_visc2',               real(mpas_v_mom_eddy_visc2,kind=RKIND))
    call mpas_pool_add_config(configPool, 'config_h_theta_eddy_visc2',             real(mpas_h_theta_eddy_visc2,kind=RKIND))
    call mpas_pool_add_config(configPool, 'config_h_theta_eddy_visc4',             real(mpas_h_theta_eddy_visc4,kind=RKIND))
    call mpas_pool_add_config(configPool, 'config_v_theta_eddy_visc2',             real(mpas_v_theta_eddy_visc2,kind=RKIND))
    call mpas_pool_add_config(configPool, 'config_horiz_mixing',                   mpas_horiz_mixing)
    call mpas_pool_add_config(configPool, 'config_len_disp',                       real(mpas_len_disp,kind=RKIND))
    call mpas_pool_add_config(configPool, 'config_visc4_2dsmag',                   real(mpas_visc4_2dsmag,kind=RKIND))
    call mpas_pool_add_config(configPool, 'config_del4u_div_factor',               real(mpas_del4u_div_factor,kind=RKIND))
    call mpas_pool_add_config(configPool, 'config_w_adv_order',                    mpas_w_adv_order)
    call mpas_pool_add_config(configPool, 'config_theta_adv_order',                mpas_theta_adv_order)
    call mpas_pool_add_config(configPool, 'config_scalar_adv_order',               mpas_scalar_adv_order)
    call mpas_pool_add_config(configPool, 'config_u_vadv_order',                   mpas_u_vadv_order)
    call mpas_pool_add_config(configPool, 'config_w_vadv_order',                   mpas_w_vadv_order)
    call mpas_pool_add_config(configPool, 'config_theta_vadv_order',               mpas_theta_vadv_order)
    call mpas_pool_add_config(configPool, 'config_scalar_vadv_order',              mpas_scalar_vadv_order)
    call mpas_pool_add_config(configPool, 'config_scalar_advection',               mpas_scalar_advection)
    call mpas_pool_add_config(configPool, 'config_positive_definite',              mpas_positive_definite)
    call mpas_pool_add_config(configPool, 'config_monotonic',                      mpas_monotonic)
    call mpas_pool_add_config(configPool, 'config_coef_3rd_order',                 real(mpas_coef_3rd_order,kind=RKIND))
    call mpas_pool_add_config(configPool, 'config_smagorinsky_coef',               real(mpas_smagorinsky_coef,kind=RKIND))
    call mpas_pool_add_config(configPool, 'config_mix_full',                       mpas_mix_full)
    call mpas_pool_add_config(configPool, 'config_epssm',                          real(mpas_epssm,kind=RKIND))
    call mpas_pool_add_config(configPool, 'config_smdiv',                          real(mpas_smdiv,kind=RKIND))
    call mpas_pool_add_config(configPool, 'config_apvm_upwinding',                 real(mpas_apvm_upwinding,kind=RKIND))
    call mpas_pool_add_config(configPool, 'config_h_ScaleWithMesh',                mpas_h_ScaleWithMesh)
    !
    call mpas_pool_add_config(configPool, 'config_zd',                             real(mpas_zd,kind=RKIND))
    call mpas_pool_add_config(configPool, 'config_xnutr',                          real(mpas_xnutr,kind=RKIND))
    call mpas_pool_add_config(configPool, 'config_mpas_cam_coef',                  real(mpas_cam_coef,kind=RKIND))
    call mpas_pool_add_config(configPool, 'config_number_cam_damping_levels',      mpas_cam_damping_levels)
    call mpas_pool_add_config(configPool, 'config_rayleigh_damp_u',                mpas_rayleigh_damp_u)
    call mpas_pool_add_config(configPool, 'config_rayleigh_damp_u_timescale_days', real(mpas_rayleigh_damp_u_timescale_days,kind=RKIND))
    call mpas_pool_add_config(configPool, 'config_number_rayleigh_damp_u_levels',  mpas_number_rayleigh_damp_u_levels)
    !
    call mpas_pool_add_config(configPool, 'config_apply_lbcs',                     mpas_apply_lbcs)
    !
    call mpas_pool_add_config(configPool, 'config_pio_num_iotasks',                mpas_pio_num_iotasks)
    call mpas_pool_add_config(configPool, 'config_pio_stride',                     mpas_pio_stride)
    !
    call mpas_pool_add_config(configPool, 'config_jedi_da',                        mpas_jedi_da)
    !
    call mpas_pool_add_config(configPool, 'config_block_decomp_file_prefix',       mpas_block_decomp_file_prefix)
    !
    call mpas_pool_add_config(configPool, 'config_do_restart',                     mpas_do_restart)
    !
    call mpas_pool_add_config(configPool, 'config_print_global_minmax_vel',        mpas_print_global_minmax_vel)
    call mpas_pool_add_config(configPool, 'config_print_detailed_minmax_vel',      mpas_print_detailed_minmax_vel)
    call mpas_pool_add_config(configPool, 'config_print_global_minmax_sca',        mpas_print_global_minmax_sca)

    ! Set some configuration parameters that cannot be changed by UFSATM. *From CAM src/dynamics/mpas/dyn_comp.F90*
    call mpas_pool_add_config(configPool, 'config_num_halos',                      config_num_halos)
    call mpas_pool_add_config(configPool, 'config_number_of_blocks',               config_number_of_blocks)
    call mpas_pool_add_config(configPool, 'config_explicit_proc_decomp',           config_explicit_proc_decomp)
    call mpas_pool_add_config(configPool, 'config_proc_decomp_file_prefix',        config_proc_decomp_file_prefix)

    ! Display namelist information (master processor only)
    if (me == master) then
      write(*,*) 'MPAS-A dycore configuration:'
      write(*,*) '   mpas_time_integration               = ', trim(mpas_time_integration)
      write(*,*) '   mpas_time_integration_order         = ', mpas_time_integration_order
      write(*,*) '   mpas_dt                             = ', mpas_dt
      write(*,*) '   mpas_split_dynamics_transport       = ', mpas_split_dynamics_transport
      write(*,*) '   mpas_number_of_sub_steps            = ', mpas_number_of_sub_steps
      write(*,*) '   mpas_dynamics_split_steps           = ', mpas_dynamics_split_steps
      write(*,*) '   mpas_h_mom_eddy_visc2               = ', mpas_h_mom_eddy_visc2
      write(*,*) '   mpas_h_mom_eddy_visc4               = ', mpas_h_mom_eddy_visc4
      write(*,*) '   mpas_v_mom_eddy_visc2               = ', mpas_v_mom_eddy_visc2
      write(*,*) '   mpas_h_theta_eddy_visc2             = ', mpas_h_theta_eddy_visc2
      write(*,*) '   mpas_h_theta_eddy_visc4             = ', mpas_h_theta_eddy_visc4
      write(*,*) '   mpas_v_theta_eddy_visc2             = ', mpas_v_theta_eddy_visc2
      write(*,*) '   mpas_horiz_mixing                   = ', trim(mpas_horiz_mixing)
      write(*,*) '   mpas_len_disp                       = ', mpas_len_disp
      write(*,*) '   mpas_visc4_2dsmag                   = ', mpas_visc4_2dsmag
      write(*,*) '   mpas_del4u_div_factor               = ', mpas_del4u_div_factor
      write(*,*) '   mpas_w_adv_order                    = ', mpas_w_adv_order
      write(*,*) '   mpas_theta_adv_order                = ', mpas_theta_adv_order
      write(*,*) '   mpas_scalar_adv_order               = ', mpas_scalar_adv_order
      write(*,*) '   mpas_u_vadv_order                   = ', mpas_u_vadv_order
      write(*,*) '   mpas_w_vadv_order                   = ', mpas_w_vadv_order
      write(*,*) '   mpas_theta_vadv_order               = ', mpas_theta_vadv_order
      write(*,*) '   mpas_scalar_vadv_order              = ', mpas_scalar_vadv_order
      write(*,*) '   mpas_scalar_advection               = ', mpas_scalar_advection
      write(*,*) '   mpas_positive_definite              = ', mpas_positive_definite
      write(*,*) '   mpas_monotonic                      = ', mpas_monotonic
      write(*,*) '   mpas_coef_3rd_order                 = ', mpas_coef_3rd_order
      write(*,*) '   mpas_smagorinsky_coef               = ', mpas_smagorinsky_coef
      write(*,*) '   mpas_mix_full                       = ', mpas_mix_full
      write(*,*) '   mpas_epssm                          = ', mpas_epssm
      write(*,*) '   mpas_smdiv                          = ', mpas_smdiv
      write(*,*) '   mpas_apvm_upwinding                 = ', mpas_apvm_upwinding
      write(*,*) '   mpas_h_ScaleWithMesh                = ', mpas_h_ScaleWithMesh
      write(*,*) '   mpas_zd                             = ', mpas_zd
      write(*,*) '   mpas_xnutr                          = ', mpas_xnutr
      write(*,*) '   mpas_cam_coef                       = ', mpas_cam_coef
      write(*,*) '   mpas_cam_damping_levels             = ', mpas_cam_damping_levels
      write(*,*) '   mpas_rayleigh_damp_u                = ', mpas_rayleigh_damp_u
      write(*,*) '   mpas_rayleigh_damp_u_timescale_days = ', mpas_rayleigh_damp_u_timescale_days
      write(*,*) '   mpas_number_rayleigh_damp_u_levels  = ', mpas_number_rayleigh_damp_u_levels
      write(*,*) '   mpas_apply_lbcs                     = ', mpas_apply_lbcs
      write(*,*) '   mpas_pio_num_iotasks                = ', mpas_pio_num_iotasks
      write(*,*) '   mpas_pio_stride                     = ', mpas_pio_stride
      write(*,*) '   mpas_jedi_da                        = ', mpas_jedi_da
      write(*,*) '   mpas_block_decomp_file_prefix       = ', trim(mpas_block_decomp_file_prefix)
      write(*,*) '   mpas_do_restart                     = ', mpas_do_restart
      write(*,*) '   mpas_print_global_minmax_vel        = ', mpas_print_global_minmax_vel
      write(*,*) '   mpas_print_detailed_minmax_vel      = ', mpas_print_detailed_minmax_vel
      write(*,*) '   mpas_print_global_minmax_sca        = ', mpas_print_global_minmax_sca
   end if

 end subroutine read_mpas_namelist

 !> ########################################################################################
 !>
 !> \brief  Reads time-invariant ("static") fields from an MPAS-A mesh file
 !> \author Michael Duda
 !> \date   6 January 2020
 !> \details
 !>  This routine takes as input an opened PIO file descriptor and a routine
 !>  to call if catastrophic errors are encountered. An MPAS stream is constructed
 !>  from this file descriptor, and most of the fields that exist in MPAS's
 !>  "mesh" pool are read from this stream.
 !>  Upon successful completion, valid mesh fields may be accessed from the mesh
 !>  pool.
 !>
 !> \update: Dustin Swales April 2025 - Modified for use in UWM
 !>
 !> ########################################################################################
 subroutine ufs_mpas_read_invariant()
   ! MPAS
   use mpas_kind_types,     only : StrKIND
   use mpas_io_streams,     only : mpas_createStream, mpas_closeStream, mpas_streamAddField
   use mpas_io_streams,     only : mpas_readStream
   use mpas_derived_types,  only : mpas_IO_READ, mpas_IO_NETCDF, mpas_Stream_type, mpas_pool_type
   use mpas_derived_types,  only : field0DReal, field1DReal, field2DReal, field3DReal
   use mpas_derived_types,  only : field1DInteger, field2DInteger, mpas_STREAM_NOERR
   use mpas_pool_routines,  only : mpas_pool_get_subpool, mpas_pool_get_field, mpas_pool_create_pool
   use mpas_pool_routines,  only : mpas_pool_destroy_pool, mpas_pool_add_config
   use mpas_dmpar,          only : mpas_dmpar_exch_halo_field
   use mpas_stream_manager, only : postread_reindex
   ! FMS
   use mpp_mod,             only : FATAL, mpp_error
   ! Arguments
   ! Locals
   character(len=*), parameter :: subname = 'ufs_mpas_subdriver::ufs_mpas_read_invariant'
   character(len=StrKIND) :: errString
   integer :: ierr
   integer :: ierr_total
   type (mpas_pool_type), pointer :: meshPool
   type (mpas_pool_type), pointer :: reindexPool
   type (mpas_pool_type), pointer :: allPackages, reindexPkgs
   type (field1DReal), pointer :: latCell, lonCell, xCell, yCell, zCell
   type (field1DReal), pointer :: latEdge, lonEdge, xEdge, yEdge, zEdge
   type (field1DReal), pointer :: latVertex, lonVertex, xVertex, yVertex, zVertex
   type (field1DInteger), pointer :: indexToCellID, indexToEdgeID, indexToVertexID
   type (field1DReal), pointer :: fEdge, fVertex
   type (field1DReal), pointer :: areaCell, areaTriangle, dcEdge, dvEdge, angleEdge
   type (field2DReal), pointer :: kiteAreasOnVertex, weightsOnEdge
   type (field1DReal), pointer :: meshDensity
   type (field1DInteger), pointer :: nEdgesOnCell, nEdgesOnEdge
   type (field2DInteger), pointer :: cellsOnEdge, edgesOnCell, edgesOnEdge, cellsOnCell, verticesOnCell, &
                                     verticesOnEdge, edgesOnVertex, cellsOnVertex
   type (field0DReal), pointer :: cf1, cf2, cf3
   type (field1DReal), pointer :: rdzw, dzu, rdzu, fzm, fzp
   type (field2DReal), pointer :: zgrid, zxu, zz
   type (field3DReal), pointer :: zb, zb3, deriv_two, cellTangentPlane, coeffs_reconstruct

   type (field2DReal), pointer :: edgeNormalVectors, localVerticalUnitVectors, defc_a, defc_b
   type (field2DReal), pointer :: cell_gradient_coef_x, cell_gradient_coef_y

   type (mpas_Stream_type) :: mesh_stream

   nullify(cell_gradient_coef_x)
   nullify(cell_gradient_coef_y)

   call mpas_createStream(mesh_stream, domain_ptr % ioContext, 'not_used', mpas_IO_NETCDF, mpas_IO_READ, &
                           pio_file_desc=pioid, ierr=ierr)
   if (ierr /= mpas_STREAM_NOERR) then
      call mpp_error(FATAL,subname//': FATAL: Failed to create static input stream.')
   end if

   call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh', meshPool)

   call mpas_pool_get_field(meshPool, 'latCell', latCell)
   call mpas_pool_get_field(meshPool, 'lonCell', lonCell)
   call mpas_pool_get_field(meshPool, 'xCell', xCell)
   call mpas_pool_get_field(meshPool, 'yCell', yCell)
   call mpas_pool_get_field(meshPool, 'zCell', zCell)

   call mpas_pool_get_field(meshPool, 'latEdge', latEdge)
   call mpas_pool_get_field(meshPool, 'lonEdge', lonEdge)
   call mpas_pool_get_field(meshPool, 'xEdge', xEdge)
   call mpas_pool_get_field(meshPool, 'yEdge', yEdge)
   call mpas_pool_get_field(meshPool, 'zEdge', zEdge)

   call mpas_pool_get_field(meshPool, 'latVertex', latVertex)
   call mpas_pool_get_field(meshPool, 'lonVertex', lonVertex)
   call mpas_pool_get_field(meshPool, 'xVertex', xVertex)
   call mpas_pool_get_field(meshPool, 'yVertex', yVertex)
   call mpas_pool_get_field(meshPool, 'zVertex', zVertex)

   call mpas_pool_get_field(meshPool, 'indexToCellID', indexToCellID)
   call mpas_pool_get_field(meshPool, 'indexToEdgeID', indexToEdgeID)
   call mpas_pool_get_field(meshPool, 'indexToVertexID', indexToVertexID)

   call mpas_pool_get_field(meshPool, 'fEdge', fEdge)
   call mpas_pool_get_field(meshPool, 'fVertex', fVertex)

   call mpas_pool_get_field(meshPool, 'areaCell', areaCell)
   call mpas_pool_get_field(meshPool, 'areaTriangle', areaTriangle)
   call mpas_pool_get_field(meshPool, 'dcEdge', dcEdge)
   call mpas_pool_get_field(meshPool, 'dvEdge', dvEdge)
   call mpas_pool_get_field(meshPool, 'angleEdge', angleEdge)
   call mpas_pool_get_field(meshPool, 'kiteAreasOnVertex', kiteAreasOnVertex)
   call mpas_pool_get_field(meshPool, 'weightsOnEdge', weightsOnEdge)

   call mpas_pool_get_field(meshPool, 'meshDensity', meshDensity)

   call mpas_pool_get_field(meshPool, 'nEdgesOnCell', nEdgesOnCell)
   call mpas_pool_get_field(meshPool, 'nEdgesOnEdge', nEdgesOnEdge)

   call mpas_pool_get_field(meshPool, 'cellsOnEdge', cellsOnEdge)
   call mpas_pool_get_field(meshPool, 'edgesOnCell', edgesOnCell)
   call mpas_pool_get_field(meshPool, 'edgesOnEdge', edgesOnEdge)
   call mpas_pool_get_field(meshPool, 'cellsOnCell', cellsOnCell)
   call mpas_pool_get_field(meshPool, 'verticesOnCell', verticesOnCell)
   call mpas_pool_get_field(meshPool, 'verticesOnEdge', verticesOnEdge)
   call mpas_pool_get_field(meshPool, 'edgesOnVertex', edgesOnVertex)
   call mpas_pool_get_field(meshPool, 'cellsOnVertex', cellsOnVertex)

   call mpas_pool_get_field(meshPool, 'cf1', cf1)
   call mpas_pool_get_field(meshPool, 'cf2', cf2)
   call mpas_pool_get_field(meshPool, 'cf3', cf3)

   call mpas_pool_get_field(meshPool, 'rdzw', rdzw)
   call mpas_pool_get_field(meshPool, 'dzu', dzu)
   call mpas_pool_get_field(meshPool, 'rdzu', rdzu)
   call mpas_pool_get_field(meshPool, 'fzm', fzm)
   call mpas_pool_get_field(meshPool, 'fzp', fzp)

   call mpas_pool_get_field(meshPool, 'zgrid', zgrid)
   call mpas_pool_get_field(meshPool, 'zxu', zxu)
   call mpas_pool_get_field(meshPool, 'zz', zz)
   call mpas_pool_get_field(meshPool, 'zb', zb)
   call mpas_pool_get_field(meshPool, 'zb3', zb3)

   call mpas_pool_get_field(meshPool, 'deriv_two', deriv_two)
   call mpas_pool_get_field(meshPool, 'cellTangentPlane', cellTangentPlane)
   call mpas_pool_get_field(meshPool, 'coeffs_reconstruct', coeffs_reconstruct)

   call mpas_pool_get_field(meshPool, 'edgeNormalVectors', edgeNormalVectors)
   call mpas_pool_get_field(meshPool, 'localVerticalUnitVectors', localVerticalUnitVectors)

   call mpas_pool_get_field(meshPool, 'defc_a', defc_a)
   call mpas_pool_get_field(meshPool, 'defc_b', defc_b)

   ierr_total = 0

   call mpas_streamAddField(mesh_stream, latCell, ierr=ierr)
   if (ierr /= mpas_STREAM_NOERR) ierr_total = ierr_total + 1
   call mpas_streamAddField(mesh_stream, lonCell, ierr=ierr)
   if (ierr /= mpas_STREAM_NOERR) ierr_total = ierr_total + 1
   call mpas_streamAddField(mesh_stream, xCell, ierr=ierr)
   if (ierr /= mpas_STREAM_NOERR) ierr_total = ierr_total + 1
   call mpas_streamAddField(mesh_stream, yCell, ierr=ierr)
   if (ierr /= mpas_STREAM_NOERR) ierr_total = ierr_total + 1
   call mpas_streamAddField(mesh_stream, zCell, ierr=ierr)
   if (ierr /= mpas_STREAM_NOERR) ierr_total = ierr_total + 1

   call mpas_streamAddField(mesh_stream, latEdge, ierr=ierr)
   if (ierr /= mpas_STREAM_NOERR) ierr_total = ierr_total + 1
   call mpas_streamAddField(mesh_stream, lonEdge, ierr=ierr)
   if (ierr /= mpas_STREAM_NOERR) ierr_total = ierr_total + 1
   call mpas_streamAddField(mesh_stream, xEdge, ierr=ierr)
   if (ierr /= mpas_STREAM_NOERR) ierr_total = ierr_total + 1
   call mpas_streamAddField(mesh_stream, yEdge, ierr=ierr)
   if (ierr /= mpas_STREAM_NOERR) ierr_total = ierr_total + 1
   call mpas_streamAddField(mesh_stream, zEdge, ierr=ierr)
   if (ierr /= mpas_STREAM_NOERR) ierr_total = ierr_total + 1

   call mpas_streamAddField(mesh_stream, latVertex, ierr=ierr)
   if (ierr /= mpas_STREAM_NOERR) ierr_total = ierr_total + 1
   call mpas_streamAddField(mesh_stream, lonVertex, ierr=ierr)
   if (ierr /= mpas_STREAM_NOERR) ierr_total = ierr_total + 1
   call mpas_streamAddField(mesh_stream, xVertex, ierr=ierr)
   if (ierr /= mpas_STREAM_NOERR) ierr_total = ierr_total + 1
   call mpas_streamAddField(mesh_stream, yVertex, ierr=ierr)
   if (ierr /= mpas_STREAM_NOERR) ierr_total = ierr_total + 1
   call mpas_streamAddField(mesh_stream, zVertex, ierr=ierr)
   if (ierr /= mpas_STREAM_NOERR) ierr_total = ierr_total + 1

   call mpas_streamAddField(mesh_stream, indexToCellID, ierr=ierr)
   if (ierr /= mpas_STREAM_NOERR) ierr_total = ierr_total + 1
   call mpas_streamAddField(mesh_stream, indexToEdgeID, ierr=ierr)
   if (ierr /= mpas_STREAM_NOERR) ierr_total = ierr_total + 1
   call mpas_streamAddField(mesh_stream, indexToVertexID, ierr=ierr)
   if (ierr /= mpas_STREAM_NOERR) ierr_total = ierr_total + 1

   call mpas_streamAddField(mesh_stream, fEdge, ierr=ierr)
   if (ierr /= mpas_STREAM_NOERR) ierr_total = ierr_total + 1
   call mpas_streamAddField(mesh_stream, fVertex, ierr=ierr)
   if (ierr /= mpas_STREAM_NOERR) ierr_total = ierr_total + 1

   call mpas_streamAddField(mesh_stream, areaCell, ierr=ierr)
   if (ierr /= mpas_STREAM_NOERR) ierr_total = ierr_total + 1
   call mpas_streamAddField(mesh_stream, areaTriangle, ierr=ierr)
   if (ierr /= mpas_STREAM_NOERR) ierr_total = ierr_total + 1
   call mpas_streamAddField(mesh_stream, dcEdge, ierr=ierr)
   if (ierr /= mpas_STREAM_NOERR) ierr_total = ierr_total + 1
   call mpas_streamAddField(mesh_stream, dvEdge, ierr=ierr)
   if (ierr /= mpas_STREAM_NOERR) ierr_total = ierr_total + 1
   call mpas_streamAddField(mesh_stream, angleEdge, ierr=ierr)
   if (ierr /= mpas_STREAM_NOERR) ierr_total = ierr_total + 1
   call mpas_streamAddField(mesh_stream, kiteAreasOnVertex, ierr=ierr)
   if (ierr /= mpas_STREAM_NOERR) ierr_total = ierr_total + 1
   call mpas_streamAddField(mesh_stream, weightsOnEdge, ierr=ierr)
   if (ierr /= mpas_STREAM_NOERR) ierr_total = ierr_total + 1

   call mpas_streamAddField(mesh_stream, meshDensity, ierr=ierr)
   if (ierr /= mpas_STREAM_NOERR) ierr_total = ierr_total + 1

   call mpas_streamAddField(mesh_stream, nEdgesOnCell, ierr=ierr)
   if (ierr /= mpas_STREAM_NOERR) ierr_total = ierr_total + 1
   call mpas_streamAddField(mesh_stream, nEdgesOnEdge, ierr=ierr)
   if (ierr /= mpas_STREAM_NOERR) ierr_total = ierr_total + 1

   call mpas_streamAddField(mesh_stream, cellsOnEdge, ierr=ierr)
   if (ierr /= mpas_STREAM_NOERR) ierr_total = ierr_total + 1
   call mpas_streamAddField(mesh_stream, edgesOnCell, ierr=ierr)
   if (ierr /= mpas_STREAM_NOERR) ierr_total = ierr_total + 1
   call mpas_streamAddField(mesh_stream, edgesOnEdge, ierr=ierr)
   if (ierr /= mpas_STREAM_NOERR) ierr_total = ierr_total + 1
   call mpas_streamAddField(mesh_stream, cellsOnCell, ierr=ierr)
   if (ierr /= mpas_STREAM_NOERR) ierr_total = ierr_total + 1
   call mpas_streamAddField(mesh_stream, verticesOnCell, ierr=ierr)
   if (ierr /= mpas_STREAM_NOERR) ierr_total = ierr_total + 1
   call mpas_streamAddField(mesh_stream, verticesOnEdge, ierr=ierr)
   if (ierr /= mpas_STREAM_NOERR) ierr_total = ierr_total + 1
   call mpas_streamAddField(mesh_stream, edgesOnVertex, ierr=ierr)
   if (ierr /= mpas_STREAM_NOERR) ierr_total = ierr_total + 1
   call mpas_streamAddField(mesh_stream, cellsOnVertex, ierr=ierr)
   if (ierr /= mpas_STREAM_NOERR) ierr_total = ierr_total + 1

   call mpas_streamAddField(mesh_stream, cf1, ierr=ierr)
   if (ierr /= mpas_STREAM_NOERR) ierr_total = ierr_total + 1
   call mpas_streamAddField(mesh_stream, cf2, ierr=ierr)
   if (ierr /= mpas_STREAM_NOERR) ierr_total = ierr_total + 1
   call mpas_streamAddField(mesh_stream, cf3, ierr=ierr)
   if (ierr /= mpas_STREAM_NOERR) ierr_total = ierr_total + 1

   call mpas_streamAddField(mesh_stream, rdzw, ierr=ierr)
   if (ierr /= mpas_STREAM_NOERR) ierr_total = ierr_total + 1
   call mpas_streamAddField(mesh_stream, dzu, ierr=ierr)
   if (ierr /= mpas_STREAM_NOERR) ierr_total = ierr_total + 1
   call mpas_streamAddField(mesh_stream, rdzu, ierr=ierr)
   if (ierr /= mpas_STREAM_NOERR) ierr_total = ierr_total + 1
   call mpas_streamAddField(mesh_stream, fzm, ierr=ierr)
   if (ierr /= mpas_STREAM_NOERR) ierr_total = ierr_total + 1
   call mpas_streamAddField(mesh_stream, fzp, ierr=ierr)
   if (ierr /= mpas_STREAM_NOERR) ierr_total = ierr_total + 1

   call mpas_streamAddField(mesh_stream, zgrid, ierr=ierr)
   if (ierr /= mpas_STREAM_NOERR) ierr_total = ierr_total + 1
   call mpas_streamAddField(mesh_stream, zxu, ierr=ierr)
   if (ierr /= mpas_STREAM_NOERR) ierr_total = ierr_total + 1
   call mpas_streamAddField(mesh_stream, zz, ierr=ierr)
   if (ierr /= mpas_STREAM_NOERR) ierr_total = ierr_total + 1
   call mpas_streamAddField(mesh_stream, zb, ierr=ierr)
   if (ierr /= mpas_STREAM_NOERR) ierr_total = ierr_total + 1
   call mpas_streamAddField(mesh_stream, zb3, ierr=ierr)
   if (ierr /= mpas_STREAM_NOERR) ierr_total = ierr_total + 1

   call mpas_streamAddField(mesh_stream, deriv_two, ierr=ierr)
   if (ierr /= mpas_STREAM_NOERR) ierr_total = ierr_total + 1
   call mpas_streamAddField(mesh_stream, cellTangentPlane, ierr=ierr)
   if (ierr /= mpas_STREAM_NOERR) ierr_total = ierr_total + 1
   call mpas_streamAddField(mesh_stream, coeffs_reconstruct, ierr=ierr)
   if (ierr /= mpas_STREAM_NOERR) ierr_total = ierr_total + 1

   call mpas_streamAddField(mesh_stream, edgeNormalVectors, ierr=ierr)
   if (ierr /= mpas_STREAM_NOERR) ierr_total = ierr_total + 1
   call mpas_streamAddField(mesh_stream, localVerticalUnitVectors, ierr=ierr)
   if (ierr /= mpas_STREAM_NOERR) ierr_total = ierr_total + 1
   call mpas_streamAddField(mesh_stream, defc_a, ierr=ierr)
   if (ierr /= mpas_STREAM_NOERR) ierr_total = ierr_total + 1
   call mpas_streamAddField(mesh_stream, defc_b, ierr=ierr)
   if (ierr /= mpas_STREAM_NOERR) ierr_total = ierr_total + 1

   if (ierr_total > 0) then
      write(errString, '(a,i0,a)') subname//': FATAL: Failed to add ', ierr_total, ' fields to static input stream.'
      call mpp_error(FATAL,trim(errString))
   end if

   call mpas_readStream(mesh_stream, 1, ierr=ierr)
   if (ierr /= mpas_STREAM_NOERR) then
      call mpp_error(FATAL,subname//': FATAL: Failed to read static input stream.')
   end if

   call mpas_closeStream(mesh_stream, ierr=ierr)
   if (ierr /= mpas_STREAM_NOERR) then
      call mpp_error(FATAL,subname//': FATAL: Failed to close static input stream.')
   end if

   !
   ! Perform halo updates for all decomposed fields (i.e., fields with
   ! an outermost dimension of nCells, nVertices, or nEdges)
   !
   call mpas_dmpar_exch_halo_field(latCell)
   call mpas_dmpar_exch_halo_field(lonCell)
   call mpas_dmpar_exch_halo_field(xCell)
   call mpas_dmpar_exch_halo_field(yCell)
   call mpas_dmpar_exch_halo_field(zCell)

   call mpas_dmpar_exch_halo_field(latEdge)
   call mpas_dmpar_exch_halo_field(lonEdge)
   call mpas_dmpar_exch_halo_field(xEdge)
   call mpas_dmpar_exch_halo_field(yEdge)
   call mpas_dmpar_exch_halo_field(zEdge)

   call mpas_dmpar_exch_halo_field(latVertex)
   call mpas_dmpar_exch_halo_field(lonVertex)
   call mpas_dmpar_exch_halo_field(xVertex)
   call mpas_dmpar_exch_halo_field(yVertex)
   call mpas_dmpar_exch_halo_field(zVertex)

   call mpas_dmpar_exch_halo_field(indexToCellID)
   call mpas_dmpar_exch_halo_field(indexToEdgeID)
   call mpas_dmpar_exch_halo_field(indexToVertexID)

   call mpas_dmpar_exch_halo_field(fEdge)
   call mpas_dmpar_exch_halo_field(fVertex)

   call mpas_dmpar_exch_halo_field(areaCell)
   call mpas_dmpar_exch_halo_field(areaTriangle)
   call mpas_dmpar_exch_halo_field(dcEdge)
   call mpas_dmpar_exch_halo_field(dvEdge)
   call mpas_dmpar_exch_halo_field(angleEdge)
   call mpas_dmpar_exch_halo_field(kiteAreasOnVertex)
   call mpas_dmpar_exch_halo_field(weightsOnEdge)

   call mpas_dmpar_exch_halo_field(meshDensity)

   call mpas_dmpar_exch_halo_field(nEdgesOnCell)
   call mpas_dmpar_exch_halo_field(nEdgesOnEdge)

   call mpas_dmpar_exch_halo_field(cellsOnEdge)
   call mpas_dmpar_exch_halo_field(edgesOnCell)
   call mpas_dmpar_exch_halo_field(edgesOnEdge)
   call mpas_dmpar_exch_halo_field(cellsOnCell)
   call mpas_dmpar_exch_halo_field(verticesOnCell)
   call mpas_dmpar_exch_halo_field(verticesOnEdge)
   call mpas_dmpar_exch_halo_field(edgesOnVertex)
   call mpas_dmpar_exch_halo_field(cellsOnVertex)

   call mpas_dmpar_exch_halo_field(zgrid)
   call mpas_dmpar_exch_halo_field(zxu)
   call mpas_dmpar_exch_halo_field(zz)
   call mpas_dmpar_exch_halo_field(zb)
   call mpas_dmpar_exch_halo_field(zb3)

   call mpas_dmpar_exch_halo_field(deriv_two)
   call mpas_dmpar_exch_halo_field(cellTangentPlane)
   call mpas_dmpar_exch_halo_field(coeffs_reconstruct)

   call mpas_dmpar_exch_halo_field(edgeNormalVectors)
   call mpas_dmpar_exch_halo_field(localVerticalUnitVectors)
   call mpas_dmpar_exch_halo_field(defc_a)
   call mpas_dmpar_exch_halo_field(defc_b)
   !
   ! Re-index from global index space to local index space
   !
   call mpas_pool_create_pool(reindexPool)

   call mpas_pool_add_config(reindexPool, 'cellsOnEdge', 1)
   call mpas_pool_add_config(reindexPool, 'edgesOnCell', 1)
   call mpas_pool_add_config(reindexPool, 'edgesOnEdge', 1)
   call mpas_pool_add_config(reindexPool, 'cellsOnCell', 1)
   call mpas_pool_add_config(reindexPool, 'verticesOnCell', 1)
   call mpas_pool_add_config(reindexPool, 'verticesOnEdge', 1)
   call mpas_pool_add_config(reindexPool, 'edgesOnVertex', 1)
   call mpas_pool_add_config(reindexPool, 'cellsOnVertex', 1)

   ! Use an empty package list for reindexPool
   call mpas_pool_create_pool(reindexPkgs)

   call postread_reindex(meshPool, domain_ptr % streamManager % allPackages, &
                         reindexPool, reindexPkgs)

   call mpas_pool_destroy_pool(reindexPool)
   call mpas_pool_destroy_pool(reindexPkgs)

 end subroutine ufs_mpas_read_invariant

 !> ########################################################################################
 !> Procedure to read MPAS IC file and populate UWM data containers.
 !>
 !> ########################################################################################
 subroutine ufs_mpas_read_init(statein)
   ! UFSATM
   use mpas_typedefs,      only : mpas_statein_type
   ! MPAS
   use mpas_kind_types,    only : RKIND
   use mpas_constants,     only : rvord
   use mpas_derived_types, only : field2DReal, field3DReal, mpas_pool_type
   use mpas_vector_reconstruction, only : mpas_reconstruct
   use mpas_pool_routines, only : mpas_pool_get_field, mpas_pool_get_subpool, mpas_pool_get_array
   use mpas_dmpar,         only : mpas_dmpar_exch_halo_field

   ! Arguments
   type (mpas_statein_type), intent(inout), target :: statein
   ! Locals
   character(len=*), parameter :: subname = 'ufs_mpas_subdriver::ufs_mpas_read_init'
   integer :: index_qv, index_qc, index_qh, index_qr, index_qi, index_qs, ierr
   ! Local MPAS pointers
   real(RKIND), pointer :: uperp(:,:)     ! Normal velocity at edges [m/s]  (nlev,nedge)
   real(RKIND), pointer :: w(:,:)         ! Vertical velocity [m/s]        (nlev+1,ncol)
   real(RKIND), pointer :: theta_m(:,:)   ! Moist potential temperature [K]  (nlev,ncol)
   real(RKIND), pointer :: rho_zz(:,:)    ! Dry density [kg/m^3]
                                          ! divided by d(zeta)/dz            (nlev,ncol)
   real(RKIND), pointer :: tracers(:,:,:) ! Tracers [kg/kg dry air]       (nq,nlev,ncol)
   real(RKIND), pointer :: zint(:,:)      ! Geometric height [m]
                                          ! at layer interfaces            (nlev+1,ncol)
   real(RKIND), pointer :: zz(:,:)        ! Vertical coordinate metric [1]
                                          ! at layer midpoints               (nlev,ncol)
   real(RKIND), pointer :: theta(:,:)     ! Potential temperature [K]        (nlev,ncol)
   real(RKIND), pointer :: rho(:,:)       ! Dry density [kg/m^3]             (nlev,ncol)
   real(RKIND), pointer :: ux(:,:)        ! Zonal veloc at center [m/s]      (nlev,ncol)
   real(RKIND), pointer :: uy(:,:)        ! Meridional veloc at center [m/s] (nlev,ncol)
   real(RKIND), pointer :: theta_base(:,:)
   real(RKIND), pointer :: rho_base(:,:)
   type (field2DReal), pointer :: field_real2d
   type (field3DReal), pointer :: field_real3d
   real(RKIND), allocatable :: mpas3d(:,:,:)
   integer :: ij
   type(mpas_pool_type), pointer :: mesh_pool
   type(mpas_pool_type), pointer :: diag_pool
   real(RKIND), pointer :: uReconstructX(:,:)
   real(RKIND), pointer :: uReconstructY(:,:)
   real(RKIND), pointer :: uReconstructZ(:,:)

   ! Local pointers
   uperp      => statein % uperp
   w          => statein % w
   theta_m    => statein % theta_m
   rho_zz     => statein % rho_zz
   tracers    => statein % tracers
   zz         => statein % zz
   theta      => statein % theta
   rho        => statein % rho
   ux         => statein % ux
   uy         => statein % uy
   rho_base   => statein % rho_base
   theta_base => statein % theta_base

   ! Tracer indices
   index_qv = statein % index_qv
   index_qc = statein % index_qc
   index_qh = statein % index_qh
   index_qr = statein % index_qr
   index_qi = statein % index_qi
   index_qs = statein % index_qs

   ! Read fields
   call ufs_mpas_read_init_field('u',          (/statein % nVertLevels,   statein % nEdgesSolve, 1/), uperp)
   call ufs_mpas_read_init_field('w',          (/statein % nVertLevels+1, statein % nCellsSolve, 1/), w)
   call ufs_mpas_read_init_field('theta',      (/statein % nVertLevels,   statein % nCellsSolve, 1/), theta)
   call ufs_mpas_read_init_field('rho',        (/statein % nVertLevels,   statein % nCellsSolve, 1/), rho)
   call ufs_mpas_read_init_field('theta_base', (/statein % nVertLevels,   statein % nCellsSolve, 1/), theta_base)
   call ufs_mpas_read_init_field('rho_base',   (/statein % nVertLevels,   statein % nCellsSolve, 1/), rho_base)
   call ufs_mpas_read_init_field('scalars',    (/statein % nVertLevels,   statein % nCellsSolve, 1/), tracers(index_qv,:,:), tracer_name='qv')
   call ufs_mpas_read_init_field('scalars',    (/statein % nVertLevels,   statein % nCellsSolve, 1/), tracers(index_qc,:,:), tracer_name='qc')
   call ufs_mpas_read_init_field('scalars',    (/statein % nVertLevels,   statein % nCellsSolve, 1/), tracers(index_qh,:,:), tracer_name='qh')
   call ufs_mpas_read_init_field('scalars',    (/statein % nVertLevels,   statein % nCellsSolve, 1/), tracers(index_qr,:,:), tracer_name='qr')
   call ufs_mpas_read_init_field('scalars',    (/statein % nVertLevels,   statein % nCellsSolve, 1/), tracers(index_qi,:,:), tracer_name='qi')
   call ufs_mpas_read_init_field('scalars',    (/statein % nVertLevels,   statein % nCellsSolve, 1/), tracers(index_qs,:,:), tracer_name='qs')   

   ! Compute derived quantities.
   theta_m(:,1:statein % nCellsSolve) = theta(:,1:statein % nCellsSolve) * (1.0_RKIND + rvord * tracers(index_qv,:,1:statein % nCellsSolve))
   rho_zz(:,1:statein % nCellsSolve)  = rho(:,1:statein % nCellsSolve) / zz(:,1:statein % nCellsSolve)
   
   ! Update halos for initial state fields
   nullify(field_real2d)
   call mpas_pool_get_field(domain_ptr % blocklist % allFields, 'u', field_real2d, timeLevel=1)
   uperp(:,statein % nEdges + 1) = uperp(:,statein % nEdges)
   call mpas_dmpar_exch_halo_field(field_real2d)

   ! Reconstruct ux and uy from uperp.
   ! This is only needed because during CAM's initialization the physics package
   ! is called before the dycore advances a step.
   nullify(mesh_pool)
   nullify(diag_pool)
   call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh', mesh_pool)
   call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'diag', diag_pool)

   ! The uReconstruct{X,Y,Z} arguments to mpas_reconstruct are required, but these
   ! field already exist in the diag pool
   nullify(uReconstructX)
   nullify(uReconstructY)
   nullify(uReconstructZ)
   call mpas_pool_get_array(diag_pool, 'uReconstructX', uReconstructX)
   call mpas_pool_get_array(diag_pool, 'uReconstructY', uReconstructY)
   call mpas_pool_get_array(diag_pool, 'uReconstructZ', uReconstructZ)

      call mpas_reconstruct(mesh_pool, uperp, &
         uReconstructX, uReconstructY, uReconstructZ, &
         ux, uy)

   

   nullify(field_real2d)
   call mpas_pool_get_field(domain_ptr % blocklist % allFields, 'w', field_real2d, timeLevel=1)
   w(:,statein % nCells + 1) = w(:,statein % nCells)
   call mpas_dmpar_exch_halo_field(field_real2d)

   nullify(field_real3d)
   call mpas_pool_get_field(domain_ptr % blocklist % allFields, 'scalars', field_real3d, timeLevel=1)
   tracers(:,:,statein % nCells + 1) =	tracers(:,:,statein % nCells)
   call mpas_dmpar_exch_halo_field(field_real3d)

   nullify(field_real2d)
   call mpas_pool_get_field(domain_ptr % blocklist % allFields, 'theta_m', field_real2d, timeLevel=1)
   theta_m(:,statein % nCells + 1) =	theta_m(:,statein % nCells)
   call mpas_dmpar_exch_halo_field(field_real2d)

   nullify(field_real2d)
   call mpas_pool_get_field(domain_ptr % blocklist % allFields, 'rho_zz', field_real2d, timeLevel=1)
   rho_zz(:,statein % nCells + 1) =	rho_zz(:,statein % nCells)
   call mpas_dmpar_exch_halo_field(field_real2d)
   
   nullify(field_real2d)
   call mpas_pool_get_field(domain_ptr % blocklist % allFields, 'rho', field_real2d)
   rho(:,statein % nCells + 1) =  rho(:,statein % nCells)
   call mpas_dmpar_exch_halo_field(field_real2d)

   theta(:,statein % nCells + 1) = theta(:,statein % nCells)
   nullify(field_real2d)
   call mpas_pool_get_field(domain_ptr % blocklist % allFields, 'theta', field_real2d)
   call mpas_dmpar_exch_halo_field(field_real2d)

   nullify(field_real2d)
   call mpas_pool_get_field(domain_ptr % blocklist % allFields, 'theta_base', field_real2d)
   theta_base(:,statein % nCells + 1) = theta_base(:,statein % nCells)
   call mpas_dmpar_exch_halo_field(field_real2d)

   nullify(field_real2d)
   call mpas_pool_get_field(domain_ptr % blocklist % allFields, 'rho_base', field_real2d)
   rho_base(:,statein % nCells + 1) =  rho_base(:,statein % nCells)
   call mpas_dmpar_exch_halo_field(field_real2d)

 end subroutine ufs_mpas_read_init

 !> ########################################################################################
 !> Procedure to read MPAS initial-condition data from opened PIO file.
 !>
 !> ########################################################################################
 subroutine ufs_mpas_read_init_field(varname, dims, varOUT, tracer_name)
   ! PIO
   use pio,                  only : var_desc_t, PIO_NOERR, PIO_inq_varid, pio_get_var
   use pio,                  only : PIO_inq_varndims, PIO_inq_vardimid, PIO_inq_dimlen
   use pio,                  only : io_desc_t, pio_initdecomp, pio_real, pio_read_darray
   ! FMS
   use mpp_mod,              only : FATAL, mpp_error
   ! MPAS
   use mpas_kind_types,      only : StrKIND, RKIND
   use mpas_pool_routines,   only : mpas_pool_get_field_info, mpas_pool_get_field
   use mpas_derived_types,   only : mpas_pool_field_info_type, field3DReal
   ! Arguments
   character(len=*), intent(in   ) :: varname
   integer,          intent(in   ) :: dims(3)
   real(RKIND),      intent(inout) :: varOUT(:,:)
   character(len=*), intent(in   ), optional :: tracer_name
   ! Locals
   character(len=*), parameter :: subname = 'ufs_mpas_subdriver::ufs_mpas_read_init_field'
   integer :: ierr, i1, i2, i3, indx, pd
   type(var_desc_t) :: varid
   type(io_desc_t) :: iodesc
   real(RKIND), allocatable :: field(:,:)
   integer :: i, ndims
   integer, dimension(:), allocatable :: dimlist, dimids
   integer, dimension(:), pointer :: indices
   integer, dimension(:), pointer :: dof
   character(len=64) :: varname_local
 
   ! Tracers are stored in 3D MPAS variable "scalars". Here we read in the tracers as
   ! 2D fields from the MPAS IC file to populate the 3D array.
   if (trim(varname) == 'scalars') then
      if (present(tracer_name)) then
         varname_local = tracer_name
      else
         varname_local = 'qv'
      endif
   else
      varname_local = varname
   endif

   ! Check that variable exists in file.
   ierr = PIO_inq_varid(pioid, trim(varname_local), varid)
   if (ierr /= PIO_NOERR) then
      call mpp_error(FATAL,subname//": variable "//trim(varname_local)//" is not on file")
   else
      ! Get dimensions
      ndims = 0
      ierr = PIO_inq_varndims(pioid, varid, ndims)
      if (ierr /= 0) call mpp_error(FATAL,subname//": Error with PIO_inq_varndims")
      allocate(dimids(ndims))
      ierr = PIO_inq_vardimid(pioid, varid, dimids(1:ndims))
      allocate(dimlist(ndims))
      do i = 1, ndims
         ierr = PIO_inq_dimlen(pioid, dimids(i), dimlist(i))
         if (ierr /= 0) call mpp_error(FATAL,subname//": Error with PIO_inq_dimlen")
      end do

      ! Get MPAS domain decomposition.
      call get_mpas_pio_decomp(varname, indices)

      ! Initialize domain decomp.
      allocate(dof(dimlist(1)*size(indices)))
      indx=1
      do i2=1,size(indices)
         do i1=1,dimlist(1)
            dof(indx) = i1 + int(indices(i2)-1)*int(dimlist(1))
            indx = indx + 1
         end do
      end do
      call pio_initdecomp(pio_subsystem, pio_real, dims, dof, iodesc)
      deallocate(dof)

      ! Read in distributed array data.
      allocate(field(dims(1), dims(2)))
      call pio_read_darray(pioid, varid, iodesc, field, ierr)
      if (ierr /= 0) call mpp_error(FATAL,subname//": Error with PIO_read_darray for "//trim(varname_local))
      varOUT(:,1:dims(2)) = field(:,:dims(2))
      deallocate(field)

   endif
   
 end subroutine ufs_mpas_read_init_field
 
 !> ########################################################################################
 !>
 !> \brief  Computes local unit north, east, and edge-normal vectors
 !> \author Michael Duda
 !> \date   15 January 2020
 !> \details
 !>  This routine computes the local unit north and east vectors at all cell
 !>  centers, storing the resulting fields in the mesh pool as 'north' and
 !>  'east'. It also computes the edge-normal unit vectors by calling
 !>  the mpas_initialize_vectors routine. Before this routine is called,
 !>  the mesh pool must contain 'latCell' and 'lonCell' fields that are valid
 !>  for all cells (not just solve cells), plus any fields that are required
 !>  by the mpas_initialize_vectors routine.
 !>
 !> \update: Dustin Swales April 2025 - Modified for use in UWM
 !>
 !> ########################################################################################
 subroutine ufs_mpas_compute_unit_vectors()
   use mpas_pool_routines,     only : mpas_pool_get_subpool, mpas_pool_get_dimension, mpas_pool_get_array
   use mpas_derived_types,     only : mpas_pool_type
   use mpas_kind_types,        only : RKIND
   use mpas_vector_operations, only : mpas_initialize_vectors

   type (mpas_pool_type), pointer :: meshPool
   real(kind=RKIND), dimension(:), pointer :: latCell, lonCell
   real(kind=RKIND), dimension(:,:), pointer :: east, north
   integer, pointer :: nCells
   integer :: iCell

   call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh', meshPool)
   call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
   call mpas_pool_get_array(meshPool, 'latCell', latCell)
   call mpas_pool_get_array(meshPool, 'lonCell', lonCell)
   call mpas_pool_get_array(meshPool, 'east', east)
   call mpas_pool_get_array(meshPool, 'north', north)

   do iCell = 1, nCells
      east(1,iCell) = -sin(lonCell(iCell))
      east(2,iCell) =  cos(lonCell(iCell))
      east(3,iCell) =  0.0_RKIND

      ! Normalize
      east(1:3,iCell) = east(1:3,iCell) / sqrt(sum(east(1:3,iCell) * east(1:3,iCell)))

      north(1,iCell) = -cos(lonCell(iCell))*sin(latCell(iCell))
      north(2,iCell) = -sin(lonCell(iCell))*sin(latCell(iCell))
      north(3,iCell) =  cos(latCell(iCell))

      ! Normalize
      north(1:3,iCell) = north(1:3,iCell) / sqrt(sum(north(1:3,iCell) * north(1:3,iCell)))

   end do

   call mpas_initialize_vectors(meshPool)

 end subroutine ufs_mpas_compute_unit_vectors
 !> ########################################################################################
 !>
 !> \brief  Returns global mesh dimensions
 !> \author Michael Duda
 !> \date   22 August 2019
 !> \details
 !>  This routine returns on all tasks the number of global cells, edges,
 !>  vertices, maxEdges, vertical layers, and the maximum number of cells owned by any task.
 !>
 !> \update: Dustin Swales April 2025 - Modified for use in UWM
 !>
 !> ########################################################################################
 subroutine ufs_mpas_get_global_dims(nCellsGlobal, nEdgesGlobal, nVerticesGlobal, maxEdges,&
      nVertLevels, maxNCells)
   use mpas_pool_routines, only : mpas_pool_get_subpool, mpas_pool_get_dimension
   use mpas_derived_types, only : mpas_pool_type
   use mpas_dmpar,         only : mpas_dmpar_sum_int, mpas_dmpar_max_int

   integer, intent(out) :: nCellsGlobal
   integer, intent(out) :: nEdgesGlobal
   integer, intent(out) :: nVerticesGlobal
   integer, intent(out) :: maxEdges
   integer, intent(out) :: nVertLevels
   integer, intent(out) :: maxNCells

   integer, pointer :: nCellsSolve
   integer, pointer :: nEdgesSolve
   integer, pointer :: nVerticesSolve
   integer, pointer :: maxEdgesLocal
   integer, pointer :: nVertLevelsLocal

   type (mpas_pool_type), pointer :: meshPool

   call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh', meshPool)
   call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
   call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve)
   call mpas_pool_get_dimension(meshPool, 'nVerticesSolve', nVerticesSolve)
   call mpas_pool_get_dimension(meshPool, 'maxEdges', maxEdgesLocal)
   call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevelsLocal)

   call mpas_dmpar_sum_int(domain_ptr % dminfo, nCellsSolve, nCellsGlobal)
   call mpas_dmpar_sum_int(domain_ptr % dminfo, nEdgesSolve, nEdgesGlobal)
   call mpas_dmpar_sum_int(domain_ptr % dminfo, nVerticesSolve, nVerticesGlobal)

   maxEdges = maxEdgesLocal
   nVertLevels = nVertLevelsLocal

   call mpas_dmpar_max_int(domain_ptr % dminfo, nCellsSolve, maxNCells)

 end subroutine ufs_mpas_get_global_dims

 ! ##########################################################################################
 !
 ! ##########################################################################################
 character(len=10) function date2yyyymmdd (date)
   ! Input arguments
   integer, intent(in) :: date

   ! Local workspace
   integer :: year    ! year of yyyy-mm-dd
   integer :: month   ! month of yyyy-mm-dd
   integer :: day     ! day of yyyy-mm-dd

   year  = date / 10000
   month = (date - year*10000) / 100
   day   = date - year*10000 - month*100

   write(date2yyyymmdd,80) year, month, day
80 format(i4.4,'-',i2.2,'-',i2.2)

 end function date2yyyymmdd
 ! #########################################################################################
 !
 ! #########################################################################################
 character(len=8) function sec2hms (seconds)
   ! Input arguments
   integer, intent(in) :: seconds

   ! Local workspace
   integer :: hours     ! hours of hh:mm:ss
   integer :: minutes   ! minutes of hh:mm:ss
   integer :: secs      ! seconds of hh:mm:ss

   hours   = seconds / 3600
   minutes = (seconds - hours*3600) / 60
   secs    = (seconds - hours*3600 - minutes*60)

   write(sec2hms,80) hours, minutes, secs
80 format(i2.2,':',i2.2,':',i2.2)

 end function sec2hms

 ! #########################################################################################
 !
 ! #########################################################################################
 character(len=10) function int2str(n)
   ! return default integer as a left justified string
   ! arguments
   integer, intent(in) :: n

   write(int2str,'(i0)') n
     
 end function int2str

 !> #########################################################################################
 !> Procedure to retreieve MPAS domain decomposition <indices>, for <varname>.
 !>
 !> ######################################################################################### 
 subroutine get_mpas_pio_decomp(varname, indices)
   use mpas_kind_types,      only : StrKIND, RKIND
   use mpas_pool_routines,   only : mpas_pool_get_field_info, mpas_pool_get_field
   use mpas_pool_routines,   only : mpas_pool_get_subpool, mpas_pool_get_array
   use mpas_pool_routines,   only : mpas_pool_get_dimension
   use mpas_derived_types,   only : mpas_pool_field_info_type, field2DReal, field3DReal
   use mpas_derived_types,   only : mpas_pool_type
   ! Arguments
   character(len=*), intent(in)  :: varname
   integer, dimension(:), pointer, intent(inout) :: indices
   ! Locals
   character(len=*), parameter :: subname = 'ufs_mpas_subdriver::get_mpas_pio_decomp'
   integer, dimension(:), pointer :: indexArray
   integer, pointer :: indexDimension
   type (field2DReal), pointer :: field2d
   type (field3DReal), pointer :: field3d
   type (mpas_pool_field_info_type) :: fieldInfo
   character (len=StrKIND) :: elementName, elementNamePlural
   logical :: meshFieldDim
   integer :: i
   
   !
   call mpas_pool_get_field_info(domain_ptr % blocklist % allFields, trim(varname), fieldInfo)
   if (trim(varname) == 'scalars') then
      nullify(field3d)
      if (fieldInfo % nTimeLevels > 1) then
         call mpas_pool_get_field(domain_ptr % blocklist % allFields, trim(varname), field3d, &
                                  timeLevel=fieldInfo % nTimeLevels )
      else
         call mpas_pool_get_field(domain_ptr % blocklist % allFields, trim(varname), field3d)
      endif
      if ( field3d % isDecomposed ) then
         meshFieldDim = .false.
         if (trim(field3d % dimNames(fieldInfo % nDims)) == 'nCells') then
            elementName = 'Cell'
            elementNamePlural = 'Cells'
            meshFieldDim = .true.
         else if (trim(field3d % dimNames(fieldInfo % nDims)) == 'nEdges') then
            elementName = 'Edge'
            elementNamePlural = 'Edges'
            meshFieldDim = .true.
         else if (trim(field3d % dimNames(fieldInfo % nDims)) == 'nVertices') then
            elementName = 'Vertex'
            elementNamePlural = 'Vertices'
            meshFieldDim = .true.
         end if
      endif
      nullify(field3d)
   else
      nullify(field2d)
      if (fieldInfo % nTimeLevels > 1) then
         call mpas_pool_get_field(domain_ptr % blocklist % allFields, trim(varname), field2d, &
                                  timeLevel=fieldInfo % nTimeLevels )
      else
         call mpas_pool_get_field(domain_ptr % blocklist % allFields, trim(varname), field2d)
      endif
      !
      if ( field2d % isDecomposed ) then
         meshFieldDim = .false.
         if (trim(field2d % dimNames(fieldInfo % nDims)) == 'nCells') then
            elementName = 'Cell'
            elementNamePlural = 'Cells'
            meshFieldDim = .true.
         else if (trim(field2d % dimNames(fieldInfo % nDims)) == 'nEdges') then
            elementName = 'Edge'
            elementNamePlural = 'Edges'
            meshFieldDim = .true.
         else if (trim(field2d % dimNames(fieldInfo % nDims)) == 'nVertices') then
            elementName = 'Vertex'
            elementNamePlural = 'Vertices'
            meshFieldDim = .true.
         end if
      endif
      nullify(field2d)
   endif
   !
   if ( meshFieldDim ) then
      allocate(indices(0))
      call mpas_pool_get_array(domain_ptr % blocklist % allFields, 'indexTo' // &
                               trim(elementName) // 'ID', indexArray)
      call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions, 'n' //  &
                                   trim(elementNamePlural) // 'Solve', indexDimension)
      call mergeArrays(indices, indexArray(1:indexDimension))

!      call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions, 'n' // &
!                                    trim(elementNamePlural) // 'Solve', indexDimension)
!      allocate(indexArray(indexDimension))
!      do i = 1, indexDimension
!         indexArray(i) = i
!      end do
!      call mergeArrays(indices, indexArray(1:indexDimension))
!      deallocate(indexArray)
   endif
   
 end subroutine get_mpas_pio_decomp
 
 subroutine mergeArrays(array1, array2)
   implicit none
   integer, dimension(:), pointer :: array1
   integer, dimension(:), intent(in) :: array2
   integer :: n1, n2
   integer, dimension(:), pointer :: newArray

   n1 = size(array1)
   n2 = size(array2)

   allocate(newArray(n1+n2))

   newArray(1:n1) = array1(:)
   newArray(n1+1:n1+n2) = array2(:)

   deallocate(array1)
   array1 => newArray

 end subroutine mergeArrays
 
end module ufs_mpas_subdriver
