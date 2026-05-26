!> ###########################################################################################
!> \file ufs_mpas_subdriver.F90
!> UFSATM subdriver for MPAS dynamical core.
!>
!> Overview:
!> Initialization is broken down into two phases, with ufs_mpas_define_scalars() called in
!> between:
!> ufs_mpas_init          :  Initialize MPAS framework, Read in namelist, Read static data.
!> ufs_mpas_atm_core_init :  Complete MPAS initialization
!>
!> Forward integration of the dycore is handled in ufs_mpas_run. The current forecast time,
!> forecast interval, and MPAS dycore time step are used to integrate the model forward in
!> time. Afterwards, atm_compute_output_diagnostics() is called to compute fields needed by
!> the Physics.
!>
!> ###########################################################################################
module ufs_mpas_subdriver
  use mpi_f08
  use mpas_kind_types,    only : StrKIND, rkind
  use mpas_derived_types, only : MPAS_LOG_ERR, MPAS_LOG_CRIT
  use mpas_log,           only : mpas_log_write
  use module_mpas_config, only : pioid_ic, pioid_restart
  use module_mpas_config, only : fcst_mpi_comm
  use module_mpas_config, only : zref, zref_edge, sphere_radius, pref, pref_edge
  use module_mpas_config, only : maxNCells, maxEdges, nVertLevels
  use module_mpas_config, only : nCellsGlobal, nEdgesGlobal, nVerticesGlobal
  use module_mpas_config, only : nEdgesSolve, nVerticesSolve, nVertLevelsSolve
  use module_mpas_config, only : dt_atmos, n_atmos
  use module_mpas_config, only : latCellGlobal, lonCellGlobal, areaCellGlobal
  use module_mpas_config, only : nml_filename, nml_funit
  use module_mpas_config, only : mpas_output_times, mpas_restart_times
  use module_mpas_config, only : out_file_index, restart_file_index
  use ufs_mpas_tools
  use ufs_mpas_io
  use ufs_mpas_boundaries
  use ufs_mpas_constituents

  implicit none

  private

  public :: MPAS_control_type
  public :: ufs_mpas_init
  public :: ufs_mpas_run

  logical :: init_lbc    = .true.
  integer :: nRecord_lbc = 1
  !> #########################################################################################
  !>
  !> #########################################################################################
  type MPAS_control_type

     ! Namelist filename
     character(len=64) :: fn_nml

     ! Full namelist for use with internal file reads
     ! This is not needed, but maintains the same interface with GFS_typedefs.F90:control_initialize()
     character(len=:), pointer, dimension(:) :: input_nml_file => null()

     ! MPI Bookkeeping
     integer          :: me           !< current MPI-rank
     integer          :: master       !< master MPI-rank
     type(MPI_Comm)   :: mpi_comm     !< forecast tasks mpi communicator

     ! ESMF
     integer          :: fcst_ntasks  !< total number of forecast tasks

     ! Log file identifier
     integer          :: nlunit       !< fortran unit number for file opens
     integer          :: logunit      !< fortran unit number for writing logfile

     ! UFS date(s) for model time.
     integer          :: bdat(8)      !< model begin date in GFS format   (same as idat)
     integer          :: cdat(8)      !< model current date in GFS format (same as jdat)

     ! Spatial/Temporal parameters for physics/dynamics coupling.
     real(rkind)      :: dt_dycore    !< dynamics time step in seconds
     real(rkind)      :: dt_phys      !< physics  time step in seconds
     integer          :: nblks        !< Number of data (physics) blocks.
     integer, pointer :: blksz(:)     !< Block size for  data blocking (default blksz(1)=[nCells])
     integer          :: levs         !< number of vertical levels

     !
     integer          :: iau_offset   !< iau running window length
     logical          :: restart      !< flag whether this is a coldstart (.false.) or a warmstart/restart (.true.)

     ! Tracers
     integer                    :: nConstituents   !< Number of constituents (tracers).
     integer                    :: nwat            !< number of hydrometeors in dcyore (including water vapor)
     character(len=32), pointer :: tracer_names(:) !< tracers names to dereference tracer id
     integer,           pointer :: tracer_types(:) !< tracers types: 0=generic, 1=chem,prog, 2=chem,diag

  end type MPAS_control_type

contains

  !> #########################################################################################
  !> Procedure to initialize UWM with MPAS dynamical core.
  !>
  !> Follows mpas_init() in MPAS-Model/src/driver/mpas_subdriver.F
  !>
  !> #########################################################################################
  subroutine ufs_mpas_init(Cfg, time_start, time_end, total_time, calendar, logUnits,        &
                           mpas_from_ufs_cnst, ufs_from_mpas_cnst, debug)
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
    use mpas_derived_types,         only : MPAS_STREAM_MGR_NOERR
    use mpas_kind_types,            only : StrKIND, RKIND
    use atm_core_interface,         only : atm_setup_core, atm_setup_domain
    use mpas_constants,             only : mpas_constants_compute_derived, pi => pii
    use mpas_attlist,               only : mpas_add_att
    use mpas_rbf_interpolation,     only : mpas_rbf_interp_initialize
    use mpas_vector_reconstruction, only : mpas_init_reconstruct
    use mpas_timekeeping,           only : mpas_NOW
    ! PIO
    use pio,                        only : pio_global, pio_get_att
    ! Arguments
    type(mpas_control_type), intent(inout) :: Cfg
    integer,                 intent(in   ) :: time_start(6), time_end(6), logUnits(2)
    integer,                 intent(in   ) :: total_time
    character(17),           intent(in   ) :: calendar
    integer, pointer,        intent(in   ) :: mpas_from_ufs_cnst(:), ufs_from_mpas_cnst(:)
    logical,                 intent(in   ) :: debug
    ! Locals
    character(len=*), parameter :: subname = 'ufs_mpas_subdriver::ufs_mpas_init'
    integer :: i, ndate1, ndate2, tod, ierr, ik, kk
    type (mpas_pool_type), pointer :: state, mesh, tend, lbc
    type (field3dReal), pointer :: scalarsField
    character (len=StrKIND), pointer :: initial_time, config_start_time
    integer, pointer :: num_scalars, mpas_from_ufs_cnst2(:), ufs_from_mpas_cnst2(:)
    logical, pointer :: config_apply_lbcs
    logical :: file_exists

    ! Setup MPAS infrastructure
    allocate(corelist, stat=ierr)
    !if ( ierr /= 0 ) call mpas_log_write(subname // " failed to allocate corelist array", messageType=MPAS_LOG_CRIT)
    nullify(corelist % next)

    allocate(corelist % domainlist, stat=ierr)
    !if ( ierr /= 0 ) call mpas_log_write(subname // " failed to allocate corelist%domainlist%next", messageType=MPAS_LOG_CRIT)
    nullify(corelist % domainlist % next)

    domain_ptr => corelist % domainlist
    domain_ptr % core => corelist

    call mpas_allocate_domain(domain_ptr)
    domain_ptr % domainID = 0

    !
    ! Initialize MPAS infrastructure (phase 1)
    !
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
       call mpas_log_write(subname // " Log setup failed for MPAS-A dycore", messageType=MPAS_LOG_CRIT)
    end if

    !
    ! Read MPAS namelist.
    !
    INQUIRE(FILE='input.nml', EXIST=file_exists)
    if (file_exists) then
       call read_mpas_namelist('input.nml', domain_ptr % configs, Cfg % mpi_comm, Cfg % master, Cfg % me)
    else
       call mpas_log_write(subname // " Cannot find MPAS namelist file, input.nml", messageType=MPAS_LOG_CRIT)
    end if

    ! Set forecast start time (config_start_time)
    ndate1 = time_start(1)*10000 + time_start(2)*100 + time_start(3)
    tod    = time_start(4)*3600  + time_start(5)*60  + time_start(6)
    call mpas_pool_add_config(domain_ptr % configs, 'config_start_time', date2yyyymmdd(ndate1)//'_'//sec2hms(tod))
    call mpas_log_write('config_start_time = '//date2yyyymmdd(ndate1)//'_'//sec2hms(tod))

    ! Set forecast end time (config_stop_time)
    ndate2 = time_end(1)*10000   + time_end(2)*100   + time_end(3)
    tod	   = time_end(4)*3600    + time_end(5)*60    + time_end(6)
    call mpas_pool_add_config(domain_ptr % configs, 'config_stop_time', date2yyyymmdd(ndate2)//'_'//sec2hms(tod))
    call mpas_log_write('config_stop_time  = '//date2yyyymmdd(ndate2)//'_'//sec2hms(tod))

    ! Set forecaste run time (config_run_duration)
    tod = max(ndate2 - ndate1 - 1,0)
    call mpas_pool_add_config(domain_ptr % configs, 'config_run_duration', trim(int2str(tod))//'_'//sec2hms(total_time))
    call mpas_log_write('config_run_duration = '//trim(int2str(tod))//'_'//sec2hms(total_time))

    ! Set other MPAS required configuration information.
    call mpas_pool_add_config(domain_ptr % configs, 'config_restart_timestamp_name', 'restart_timestamp')
    call mpas_pool_add_config(domain_ptr % configs, 'config_IAU_option',             'off')
    call mpas_pool_add_config(domain_ptr % configs, 'config_do_DAcycling',           .false.)
    call mpas_pool_add_config(domain_ptr % configs, 'config_halo_exch_method',       'mpas_halo')

    !
    ! Initialize MPAS infrastructure (phase 2)
    !
    call mpas_framework_init_phase2(domain_ptr, io_system=pio_subsystem_ic, calendar = trim(calendar))

    !
    ! Before defining packages, initialize the stream inquiry instance for the domain
    !
    domain_ptr % streamInfo => mpas_stream_inquiry_new_streaminfo()
    if (.not. associated(domain_ptr % streamInfo)) then
       call mpas_log_write(subname // " Failed to instantiate streamInfo object for "// &
                           trim(domain_ptr % core % coreName), messageType=MPAS_LOG_CRIT)
    end if

    ierr = domain_ptr % core % define_packages(domain_ptr % packages)
    if (ierr /= 0) then
       call mpas_log_write(subname // " Package definition failed for "// &
                           trim(domain_ptr % core % coreName), messageType=MPAS_LOG_CRIT)
    end if

    ierr = domain_ptr % core % setup_packages(domain_ptr % configs,  domain_ptr % streamInfo,       &
                                              domain_ptr % packages, domain_ptr % iocontext)
    if (ierr /= 0) then
       call mpas_log_write(subname // " Package setup failed for "// &
            trim(domain_ptr % core % coreName), messageType=MPAS_LOG_CRIT)
    end if

    ierr = domain_ptr % core % setup_decompositions(domain_ptr % decompositions)
    if (ierr /= 0) then
       call mpas_log_write(subname // " Decomposition setup failed for "// &
            trim(domain_ptr % core % coreName), messageType=MPAS_LOG_CRIT)
    end if

    ierr = domain_ptr % core % setup_clock(domain_ptr % clock, domain_ptr % configs)
    if (ierr /= 0) then
       call mpas_log_write(subname // " Clock setup failed for "// &
            trim(domain_ptr % core % coreName), messageType=MPAS_LOG_CRIT)
    end if

    ! Adding a config named 'cam_pcnst' with the number of constituents will indicate to
    ! MPAS-A setup code that it is operating as a UFS dycore, and that it is necessary to
    ! allocate scalars separately from other Registry-defined fields
    call mpas_pool_add_config(domain_ptr % configs, 'cam_pcnst', Cfg % nConstituents)

    ! Call MPAS framework bootstrap (phase 1)
    call mpas_bootstrap_framework_phase1(domain_ptr, "external mesh file", mpas_IO_NETCDF, pio_file_desc=pioid_ic)

    !
    ! Finalize the setup of blocks and fields
    !
    call mpas_bootstrap_framework_phase2(domain_ptr, pio_file_desc=pioid_ic)

    !
    ! END OF MPAS-Model/src/driver/mpas_subdriver.F:mpas_init()
    !

    ! Add num_scalars from "state" pool to "dimensions".
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'state', state)
    call mpas_pool_get_dimension(state, 'num_scalars', num_scalars)
    call mpas_pool_add_dimension(domain_ptr % blocklist % dimensions, 'num_scalars', num_scalars)

    !
    ! Setup scalars for State pool and scalars for Tend pool
    !
    call mpas_pool_add_dimension(state, 'index_qv', 1)
    call mpas_pool_add_dimension(state, 'moist_start', 1)
    call mpas_pool_add_dimension(state, 'moist_end', Cfg % nwat)
    nullify (state)
    call ufs_mpas_define_scalars(mpas_from_ufs_cnst, ufs_from_mpas_cnst, ierr)
    if (ierr /= 0) then
       call mpas_log_write(subname // " Set-up of constituents for MPAS-A dycore failed.", messageType=MPAS_LOG_CRIT)
    end if

    !
    ! Setup scalars for LBC pool and scalars_tend for LBC pool.
    !
    call mpas_pool_get_config( domain_ptr % blocklist % configs, 'config_apply_lbcs', config_apply_lbcs)
    if (config_apply_lbcs) then
       call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'lbc', lbc)
       call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions, 'num_scalars', num_scalars)
       call mpas_pool_add_dimension(lbc, 'num_scalars', num_scalars)
       call mpas_pool_add_dimension(lbc, 'moist_start', 1)
       call mpas_pool_add_dimension(lbc, 'moist_end', Cfg % nwat)
       call mpas_pool_add_dimension(lbc, 'index_qv', 1)
       nullify (lbc)
       call ufs_mpas_define_lbc_scalars(mpas_from_ufs_cnst, ufs_from_mpas_cnst, ierr)
       if (ierr /= 0) then
          call mpas_log_write(subname // " Set-up of LBC constituents for MPAS-A dycore failed.", messageType=MPAS_LOG_CRIT)
       end if
    end if

    !
    ! Read in static (invariant) data
    !
    call dyn_mpas_read_write_stream(domain_ptr % clock,  'r', 'invariant', pio_file_desc=pioid_ic, &
                                    ierr=ierr, timeLevel=1, whence=mpas_NOW, nRecord=1, debug=debug)
    if (ierr /= MPAS_STREAM_MGR_NOERR) then
       call mpas_log_write(subname // " Could not read from ''invariant'' stream ",messageType=MPAS_LOG_CRIT)
    end if

    ! FROM CAM/driver/cam_mpas_subdriver.F90
    ! Compute unit vectors giving the local north and east directions as well as
    ! the unit normal vector for edges
    call ufs_mpas_compute_unit_vectors()

    ! FROM CAM/dyn_grid.F90:setup_time_invariant()
    ! Initialize fields needed for reconstruction of cell-centered winds from edge-normal winds
    ! Note: This same pair of calls happens a second time later in the initialization of
    !       the MPAS-A dycore (in atm_mpas_init_block), but the redundant calls do no harm
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh', mesh)
    call mpas_rbf_interp_initialize(mesh)
    call mpas_init_reconstruct(mesh)
    nullify (mesh)

    ! Read the global sphere_radius attribute.  This is needed to normalize the cell areas.
    ierr = pio_get_att(pioid_ic, pio_global, 'sphere_radius', domain_ptr % sphere_radius)
    if( ierr /= 0 ) then
       call mpas_log_write(subname // " Could not find sphere_radius PIO attribute",messageType=MPAS_LOG_CRIT)
    endif

    ! FROM CAM/dyn_grid.F90:dyn_grid_init()
    ! Query global grid dimensions from MPAS
    call ufs_mpas_get_global_dims(nCellsGlobal, nEdgesGlobal, nVerticesGlobal, maxEdges, nVertLevels, maxNCells)

    ! Setup constants
    call mpas_constants_compute_derived()

    ! Set MPAS mesh lon/lat/area.
    allocate(latCellGlobal(nCellsGlobal), lonCellGlobal(nCellsGlobal), areaCellGlobal(nCellsGlobal))
    call ufs_mpas_get_global_coords(latCellGlobal, lonCellGlobal, areaCellGlobal)

    !
    ! Initialize core
    !
    call ufs_mpas_atm_core_init(Cfg, debug)

  end subroutine ufs_mpas_init

  !> ########################################################################################
  !> Procedure to initialize UWM with MPAS dynamical core.
  !>
  !> Follows atm_core_init() in MPAS-Model/src/core_atmosphere/mpas_atm_core.F.
  !>
  !> ########################################################################################
  subroutine ufs_mpas_atm_core_init(Cfg, debug)
    use mpas_kind_types,            only : StrKIND, RKIND
    use mpas_derived_types,         only : mpas_pool_type, mpas_Time_Type, field0DReal, field2dreal
    use mpas_derived_types,         only : block_type, field3dreal, MPAS_STREAM_MGR_NOERR
    use mpas_domain_routines,       only : mpas_pool_get_dimension
    use mpas_pool_routines,         only : mpas_pool_get_subpool
    use mpas_pool_routines,         only : mpas_pool_initialize_time_levels, mpas_pool_get_config
    use mpas_pool_routines,         only : mpas_pool_get_array, mpas_pool_get_field
    use mpas_atm_dimensions,        only : mpas_atm_set_dims
    use mpas_atm_threading,         only : mpas_atm_threading_init
    use mpas_atm_halos,             only : atm_build_halo_groups, exchange_halo_group
    use atm_core,                   only : atm_mpas_init_block
    use atm_time_integration,       only : mpas_atm_dynamics_checks
    use atm_time_integration,       only : mpas_atm_dynamics_init
    use mpas_timekeeping,           only : mpas_get_clock_time, mpas_get_time, mpas_START_TIME
    use mpas_timekeeping,           only : mpas_NOW, mpas_set_timeInterval, operator(+)
    use mpas_attlist,               only : mpas_modify_att
    use mpas_string_utils,          only : mpas_string_replace
    use mpas_field_routines,        only : mpas_allocate_scratch_field
    ! Arguments
    type(mpas_control_type), intent(inout) :: Cfg
    logical,                 intent(in   ) :: debug
    type(mpas_pool_type), pointer :: tend_physics_pool
    ! Locals
    character(len=*), parameter :: subname = 'ufs_mpas_subdriver::ufs_mpas_atm_core_init'
    type (mpas_pool_type), pointer :: state, mesh, diag
    integer :: ierr
    integer, pointer :: nVertLevels1, maxEdges1, maxEdges2, num_scalars
    real (kind=RKIND), pointer :: dt
    type (block_type), pointer :: block
    logical, pointer :: config_do_restart
    type (mpas_Time_Type) :: startTime
    character(len=StrKIND) :: startTimeStamp
    character (len=StrKIND), pointer :: xtime
    character (len=StrKIND), pointer :: initial_time1, initial_time2
    real(RKIND), dimension(:,:,:), pointer :: field_3d_real
    logical, pointer :: config_apply_lbcs
    real(RKIND), dimension(:,:), pointer :: theta1

    !
    ! Setup threading
    !
    call mpas_log_write('Setting up OpenMP threading')
    call mpas_atm_threading_init(domain_ptr%blocklist, ierr)
    if ( ierr /= 0 ) then
       call mpas_log_write(subname // " Threading setup failed for core "//trim(domain_ptr % core % coreName))
    end if

    !
    ! Set up inner dimensions used by arrays in optimized dynamics routines
    !
    call mpas_log_write('Setting up dimensions')
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'state', state)
    call mpas_pool_get_dimension(state, 'nVertLevels', nVertLevels1)
    call mpas_pool_get_dimension(state, 'maxEdges', maxEdges1)
    call mpas_pool_get_dimension(state, 'maxEdges2', maxEdges2)
    call mpas_pool_get_dimension(state, 'num_scalars', num_scalars)
    call mpas_atm_set_dims(nVertLevels1, maxEdges1, maxEdges2, num_scalars)
    Cfg % levs = nVertLevels1 !DJS: Do we need this?
    nullify (state)

    !
    ! Set "local" clock to point to the clock contained in the domain type
    !
    clock => domain_ptr % clock

    !
    ! Build halo exchange groups and set method for exchanging halos in a group
    !
    call mpas_log_write('Building halo exchange groups.')
    nullify(exchange_halo_group)
    call atm_build_halo_groups(domain_ptr, ierr)
    if (ierr /= 0) then
       call mpas_log_write(subname // " failed to build MPAS-A halo exchange groups.",messageType=MPAS_LOG_CRIT)
    end if
    if (.not. associated(exchange_halo_group)) then
       call mpas_log_write(subname // " failed to build MPAS-A halo exchange groups.",messageType=MPAS_LOG_CRIT)
    endif

    !
    call mpas_pool_get_config(domain_ptr % blocklist % configs, 'config_do_restart', config_do_restart)
    call mpas_pool_get_config(domain_ptr % blocklist % configs, 'config_dt', dt)

    !
    ! Read in initial-conditions
    !
    call mpas_log_write('Reading in MPAS initial condition stream.')
    call dyn_mpas_read_write_stream(clock, 'r', 'input', pio_file_desc=pioid_ic, ierr=ierr, &
                                    timeLevel=1, whence=mpas_NOW, nRecord=1, debug=debug)
    if (ierr /= MPAS_STREAM_MGR_NOERR) then
       call mpas_log_write(subname // " Could not read from ''input'' stream ",messageType=MPAS_LOG_CRIT)
    end if
    call dyn_mpas_read_write_stream(clock, 'r', 'sfc_input', pio_file_desc=pioid_ic, ierr=ierr, &
                                    timeLevel=1, whence=mpas_NOW, nRecord=1, debug=debug)
    if (ierr /= MPAS_STREAM_MGR_NOERR) then
       call mpas_log_write(subname // " Could not read from ''sfc_input'' stream ",messageType=MPAS_LOG_CRIT)
    end if

    !
    ! Read in restart data.
    !
    if (config_do_restart) then
       call mpas_log_write('Reading in MPAS restart stream.')
       call dyn_mpas_read_write_stream(clock, 'r', 'restart', pio_file_desc=pioid_restart, ierr=ierr, &
                                       timeLevel=1, whence=mpas_NOW, nRecord=1, debug=debug)
    end if

    if (.not. config_do_restart) then
       call mpas_log_write('Initializing time levels')
       call mpas_pool_get_subpool(domain_ptr % blocklist  % structs, 'state', state)
       call mpas_pool_initialize_time_levels(state)
       nullify (state)
    end if

    call mpas_log_write('Initializing atmospheric variables')

    ! How many calls to MPAS dycore for each ATMosphere time step?
    Cfg%dt_dycore = dt    ! DJS: Does this need to be here?
    n_atmos = dt_atmos/dt ! DJS: Does this need to be here?

    !
    ! Set startTimeStamp based on the start time of the simulation clock
    !
    startTime = mpas_get_clock_time(clock, mpas_START_TIME, ierr)
    if ( ierr /= 0 ) then
       call mpas_log_write(subname // " Failed to get clock_time mpas_START_TIME",messageType=MPAS_LOG_CRIT)
    end if
    call mpas_get_time(startTime, dateTimeString=startTimeStamp, ierr=ierr)
    if ( ierr /= 0 ) then
       call mpas_log_write(subname // " Failed to get time mpas_START_TIME",messageType=MPAS_LOG_CRIT)
    end if
    call mpas_log_write('Setting simulation start time :'//startTimeStamp)

    !
    call exchange_halo_group(domain_ptr, 'initialization:u',ierr=ierr)
    if ( ierr /= 0 ) then
       call mpas_log_write(subname // ' Failed to exchange halo layers for group "initialization:u"',messageType=MPAS_LOG_CRIT)
    end if

    !
    ! Perform basic compatibility checks among the fields that were read and the run-time options that were selected
    !
    call mpas_atm_dynamics_checks(domain_ptr % dminfo, domain_ptr % blocklist, domain_ptr % streamManager, ierr)
    if (ierr /= 0) then
       call mpas_log_write(subname // " Failed dynamics compatibility test.",messageType=MPAS_LOG_CRIT)
    end if
    call mpas_pool_get_config( domain_ptr % blocklist % configs, 'config_apply_lbcs', config_apply_lbcs)
    if (config_apply_lbcs) then
       call ufs_mpas_atm_bdy_checks(domain_ptr % dminfo, domain_ptr % blocklist, ierr)
       if (ierr /= 0) then
          call mpas_log_write(subname // " Failed regional compatibility test.",messageType=MPAS_LOG_CRIT)
       end if
    end if

    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh', mesh)
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'state', state)

    call atm_mpas_init_block(domain_ptr % dminfo, domain_ptr % streamManager, domain_ptr % blocklist, mesh, dt)
    nullify (mesh)

    call mpas_pool_get_array(state, 'xtime', xtime, timelevel=1)
    xtime = startTimeStamp

    ! Initialize initial_time in second time level. We need to do this because initial state
    ! is read into time level 1, and if we write output from the set of state arrays that
    ! represent the original time level 2, the initial_time field will be invalid.
    call mpas_pool_get_array(state, 'initial_time', initial_time1, timelevel=1)
    call mpas_pool_get_array(state, 'initial_time', initial_time2, timelevel=2)
    initial_time2 = initial_time1
    nullify (state)

    call exchange_halo_group(domain_ptr, 'initialization:pv_edge,ru,rw',ierr=ierr)
    if ( ierr /= 0 ) then
       call mpas_log_write(subname // ' Failed to exchange halo layers for group "initialization:ru,rw"',messageType=MPAS_LOG_CRIT)
    end if

    !
    ! Prepare the dynamics for integration
    !
    call mpas_log_write('Initializing the dynamics')
    call mpas_atm_dynamics_init(domain_ptr)

    call mpas_log_write('Successful initialization of MPAS dynamical core')

  end subroutine ufs_mpas_atm_core_init

  !> #########################################################################################
  !> Routine to call MPAS dynamical core
  !> Loop over dynamical time-step(s) and increment MPAS state (timelevel 1->2)
  !>
  !> #########################################################################################
  subroutine ufs_mpas_run(mpasClock, outClock, debug)
    ! MPAS
    use atm_core,             only : atm_do_timestep, atm_compute_output_diagnostics
    use mpas_domain_routines, only : mpas_pool_get_dimension
    use mpas_derived_types,   only : mpas_Time_type, mpas_pool_type, MPAS_TimeInterval_type, field2DReal
    use mpas_kind_types,      only : StrKIND, RKIND, R8KIND
    use mpas_constants,       only : rvord
    use mpas_pool_routines,   only : mpas_pool_get_config, mpas_pool_get_subpool
    use mpas_pool_routines,   only : mpas_pool_shift_time_levels, mpas_pool_get_array
    use mpas_timer,           only : mpas_timer_start, mpas_timer_stop
    use mpas_timekeeping,     only : mpas_advance_clock, mpas_get_clock_time, mpas_get_time
    use mpas_timekeeping,     only : mpas_NOW, mpas_is_clock_stop_time, mpas_dmpar_get_time
    use mpas_timekeeping,     only : mpas_set_timeInterval, operator(+), operator(.LT.), operator(.GT.), operator(.LE.), operator(.EQ.)
    ! Arguments
    real(kind=R8KIND), intent(inout) :: mpasClock,outClock
    logical, intent(in   ) :: debug
    ! Locals
    character(len=*), parameter :: subname = 'ufs_mpas_run::ufs_mpas_run'
    real (kind=RKIND), pointer :: config_dt
    type (mpas_pool_type), pointer :: state, diag, mesh
    type (mpas_Time_type) :: timeNow, timeStop,timeLBCnew
    character(len=StrKIND) :: timeStamp, timeStampOutFile
    integer :: ierr, itime, itimestep, iout
    real (kind=R8KIND) :: integ_start_time, integ_stop_time
    logical, pointer :: config_apply_lbcs
    type(mpas_timeinterval_type) :: mpas_time_interval, mpas_output_interval, mpas_restart_interval
    real (kind=RKIND), dimension(:,:,:), pointer :: scalars
    real (kind=RKIND) :: start_time, stop_time
    
    ! Start dynamics timer
    start_time = MPI_Wtime()

    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'state', state)
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'diag',  diag)
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh',  mesh)
    call atm_compute_output_diagnostics(state, 1, diag, mesh)

    ! Grab runtime configuration
    call mpas_pool_get_config( domain_ptr % blocklist % configs, 'config_dt',         config_dt)
    call mpas_pool_get_config( domain_ptr % blocklist % configs, 'config_apply_lbcs', config_apply_lbcs)
    
    ! Set up clock
    timeNow  = mpas_get_clock_time(clock, mpas_NOW, ierr=ierr)
    if (ierr /= 0) then
       call mpas_log_write(subname // ' Failed to get clock_time for "mpas_NOW"',messageType=MPAS_LOG_CRIT)
       
    endif

    call mpas_get_time(curr_time=timeNow, dateTimeString=timeStamp, ierr=ierr)
    if (ierr /= 0) then
       call mpas_log_write(subname // ' Failed to get clock_time for "mpas_NOW"',messageType=MPAS_LOG_CRIT)
    endif

    ! Set dycore interval
    call MPAS_set_timeInterval(mpas_time_interval, S=dt_atmos, ierr=ierr)
    if (ierr /= 0) then
       call mpas_log_write(subname // ' Failed to set dynamics time step',messageType=MPAS_LOG_CRIT)
    endif
    
    !
    ! Read initial boundary state
    ! During integration, time level 1 stores the boundary tendencies (next-current) file records,
    ! and time level 2 stores the state at the next file record.
    if (config_apply_lbcs .and. init_lbc) then
       call mpas_log_write('--------------------------------------------------')
       call mpas_log_write('Compute initial lateral boundary conditions for timestep '//trim(timeStamp))
       call ufs_mpas_atm_update_bdy_tend(clock, domain_ptr % blocklist, .true., nRecord_lbc, ierr, debug)
       if (ierr /= 0) then
          call mpas_log_write('Failed to process LBC data at next time after '//trim(timeStamp), messageType=MPAS_LOG_ERR)
          return
       end if
       init_lbc = .false.
    end if

    ! During integration, time level 1 stores the model state at the beginning of the
    !   time step, and time level 2 stores the state advanced config_dt in time by timestep
    timeStop = timeNow + mpas_time_interval
    itimestep =	0
    call mpas_log_write('--------------------------------------------------')
    call mpas_log_write('MPAS dynamics start timestep')
    do while (timeNow .LT. timeStop)
       itimestep = itimestep + 1

       call mpas_get_time(curr_time=timeNow, dateTimeString=timeStamp, ierr=ierr)
       if ( ierr /= 0 ) then
          call mpas_log_write(subname // ' Failed to get time "mpas_NOW"',messageType=MPAS_LOG_CRIT)
       end if
       call mpas_log_write(' Start timestep at '//trim(timeStamp))

       !
       ! Read future boundary state and compute boundary tendencies
       !
       if (config_apply_lbcs) then
          if (LBC_intv_end .LE. timeNow) then
             nRecord_lbc = nRecord_lbc + 1
             call mpas_log_write('--------------------------------------------------')
             call mpas_log_write('Update lateral boundary conditions for timestep '//trim(timeStamp))
             call ufs_mpas_atm_update_bdy_tend(clock, domain_ptr % blocklist, .false., nRecord_lbc, ierr, debug)
             if (ierr /= 0) then
                call mpas_log_write('Failed to process LBC data at next time after '//trim(timeStamp), messageType=MPAS_LOG_ERR)
                return
             end if
          end if
       end if

       ! Integrate forward one dycore time step
       call mpas_timer_start('time integration')
       call mpas_dmpar_get_time(integ_start_time)
       call atm_do_timestep(domain_ptr, config_dt, itimestep)
       call mpas_dmpar_get_time(integ_stop_time)
       call mpas_timer_stop('time integration')
       call mpas_log_write(' Timing for integration step: $r s', realArgs=(/real(integ_stop_time - integ_start_time, kind=RKIND)/))

       ! Move time level 2 fields back into time level 1 for next time step
       call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'state', state)
       call mpas_pool_shift_time_levels(state)

       ! Advance clock.
       call mpas_advance_clock(clock, ierr=ierr)
       if (ierr /= 0) then
          call mpas_log_write(subname // ' Failed to advance clock',messageType=MPAS_LOG_CRIT)
       endif
       timeNow = mpas_get_clock_time(clock, mpas_NOW, ierr=ierr)
       if (ierr /= 0) then
          call mpas_log_write(subname // ' Failed to get clock_time for "mpas_NOW"',messageType=MPAS_LOG_CRIT)
       endif
    end do
    call mpas_log_write('MPAS dynamics stop timestep')
    stop_time = MPI_Wtime()
    mpasClock = mpasCLock + + (stop_time - start_time)

    !
    ! Compute diagnostic fields  (theta, rho, pres) from
    ! the final prognostic state (theta_m, rho_zz, zz)
    !
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'state', state)
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'diag',  diag)
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh',  mesh)
    call atm_compute_output_diagnostics(state, 1, diag, mesh)

    !
    ! Write any output streams
    !
    start_time = MPI_Wtime()
    call mpas_get_time(curr_time=timeStop, dateTimeString=timeStamp, ierr=ierr)
    if ( ierr /= 0 ) then
       call mpas_log_write(subname // ' Failed to get time timeStop"',messageType=MPAS_LOG_CRIT)
    end if
    call create_file_timeStamp(timeStop,timeStampOutFile)

    ! Output stream
    if (timeStop .EQ. mpas_output_times(out_file_index)) then
       call ufs_mpas_write("output", timeStampOutFile, debug)
       out_file_index = out_file_index + 1
    end if

    ! Restart stream
    if (allocated(mpas_restart_times)) then
       if (timeStop .EQ. mpas_restart_times(restart_file_index)) then
          call ufs_mpas_write("restart", timeStampOutFile, debug)
          restart_file_index = restart_file_index + 1
       end if
    end if

    ! Stop timer outClock
    stop_time = MPI_Wtime()
    outClock = outClock + (stop_time - start_time)

  end subroutine ufs_mpas_run

  !> #########################################################################################
  !> Procedure to read MPAS namelist(s).
  !>
  !> The namelist for MPAS are described in MPAS-Model/src/core_atmosphere/Registry.xml, this
  !> is also where the default values defined below originate.
  !>
  !> #########################################################################################
  subroutine read_mpas_namelist(nml_file,configPool, mpicomm, master, me)
    use mpi_f08,            only: MPI_Comm, MPI_CHARACTER, MPI_INTEGER, MPI_REAL8,  MPI_LOGICAL
    use mpi_f08,            only: mpi_bcast, mpi_barrier
    use mpas_derived_types, only: mpas_pool_type
    use mpas_kind_types,    only: StrKIND, RKIND
    use mpas_pool_routines, only: mpas_pool_add_config
    use mpas_typedefs,      only: r8 => kind_dbl_prec

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
    real(r8)                :: mpas_len_disp                       = 0.0_r8
    real(r8)                :: mpas_visc4_2dsmag                   = 0.05_r8
    real(r8)                :: mpas_del4u_div_factor               = 10.0_r8
    integer                 :: mpas_w_adv_order                    = 3
    integer                 :: mpas_theta_adv_order                = 3
    integer                 :: mpas_scalar_adv_order               = 3
    real(r8)                :: mpas_h_scalar_filter4               = 0.0_r8
    logical                 :: mpas_scalar_eddy_mix                = .false.
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
    integer                 :: mpas_cam_damping_levels             = 4
    logical                 :: mpas_rayleigh_damp_u                = .false.
    real(r8)                :: mpas_rayleigh_damp_u_timescale_days = 5.0_r8
    integer                 :: mpas_number_rayleigh_damp_u_levels  = 6
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
         mpas_w_adv_order, mpas_theta_adv_order, mpas_scalar_adv_order, mpas_h_scalar_filter4,&
         mpas_scalar_eddy_mix, mpas_u_vadv_order,                                             &
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
    real(RKIND)            :: config_relax_zone_divdamp_coef = 6

    ! Locals
    integer :: ierr, io, mpierr
    character(len=*), parameter :: subname = 'ufs_mpas_subdriver::read_mpas_namelist'
    logical :: file_exists

    ! Read in namelists...
    if (me == master) then
       inquire(file = trim(nml_filename), exist=file_exists)
       if (file_exists) then
          call mpas_log_write('Reading MPAS-A dynamical core namelist')
          open(newunit=nml_funit,file=trim(nml_filename),status='unknown')
          ! nhyd_model
          read(nml_funit, nml=mpas_nhyd_model, iostat=io)
          if (io .ne. 0) call mpas_log_write(subname // ' Reading in MPAS namelist mpas_nhyd_model',messageType=MPAS_LOG_CRIT)
          ! damping
          read(nml_funit, nml=mpas_damping, iostat=io)
          if (io .ne. 0) call mpas_log_write(subname // ' Reading in MPAS namelist mpas_damping',messageType=MPAS_LOG_CRIT)
          ! limited_area
          read(nml_funit, nml=mpas_limited_area, iostat=io)
          if (io .ne. 0) call mpas_log_write(subname // ' Reading in MPAS namelist mpas_limited_area',messageType=MPAS_LOG_CRIT)
          ! PIO
          read(nml_funit, nml=mpas_io, iostat=io)
          if (io .ne. 0) call mpas_log_write(subname // ' Reading in MPAS namelist mpas_io',messageType=MPAS_LOG_CRIT)
          ! assimilation
          read(nml_funit, nml=mpas_assimilation, iostat=io)
          if (io .ne. 0) call mpas_log_write(subname // ' Reading in MPAS namelist mpas_assimilation',messageType=MPAS_LOG_CRIT)
          ! decomposition
          read(nml_funit, nml=mpas_decomposition, iostat=io)
          if (io .ne. 0) call mpas_log_write(subname // ' Reading in MPAS namelist mpas_decomposition',messageType=MPAS_LOG_CRIT)
          ! restart
          read(nml_funit, nml=mpas_restart, iostat=io)
          if (io .ne. 0) call mpas_log_write(subname // ' Reading in MPAS namelist mpas_restart',messageType=MPAS_LOG_CRIT)
          ! printout
          read(nml_funit, nml=mpas_printout, iostat=io)
          if (io .ne. 0) call mpas_log_write(subname // ' Reading in MPAS namelist mpas_printout',messageType=MPAS_LOG_CRIT)
       endif
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
    call mpi_bcast(mpas_h_scalar_filter4,               1, mpi_real8,     master, mpicomm, mpierr)
    call mpi_bcast(mpas_scalar_eddy_mix,                1, mpi_logical,   master, mpicomm, mpierr)
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
    call mpas_pool_add_config(configPool, 'config_h_scalar_filter4',               real(mpas_h_scalar_filter4))
    call mpas_pool_add_config(configPool, 'config_scalar_eddy_mix',                mpas_scalar_eddy_mix)
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
    call mpas_pool_add_config(configPool, 'config_relax_zone_divdamp_coef',        config_relax_zone_divdamp_coef)

    ! Display namelist information (master processor only)
    if (me == master) then
       call mpas_log_write('-------------------------------- MPAS-A dycore namelist ---------------------------------')
       call mpas_log_write('')
       call mpas_log_write('   mpas_time_integration               = '//trim(mpas_time_integration))
       call mpas_log_write('   mpas_time_integration_order         = '//int2str(mpas_time_integration_order))
       call mpas_log_write('   mpas_dt                             = '//int2str(int(mpas_dt)))
       call mpas_log_write('   mpas_split_dynamics_transport       = '//log2str(mpas_split_dynamics_transport))
       call mpas_log_write('   mpas_number_of_sub_steps            = '//int2str(mpas_number_of_sub_steps))
       call mpas_log_write('   mpas_dynamics_split_steps           = '//int2str(mpas_dynamics_split_steps))
       call mpas_log_write('   mpas_h_mom_eddy_visc2               = '//int2str(int(mpas_h_mom_eddy_visc2)))
       call mpas_log_write('   mpas_h_mom_eddy_visc4               = '//int2str(int(mpas_h_mom_eddy_visc4)))
       call mpas_log_write('   mpas_v_mom_eddy_visc2               = '//int2str(int(mpas_v_mom_eddy_visc2)))
       call mpas_log_write('   mpas_h_theta_eddy_visc2             = '//int2str(int(mpas_h_theta_eddy_visc2)))
       call mpas_log_write('   mpas_h_theta_eddy_visc4             = '//int2str(int(mpas_h_theta_eddy_visc4)))
       call mpas_log_write('   mpas_v_theta_eddy_visc2             = '//int2str(int(mpas_v_theta_eddy_visc2)))
       call mpas_log_write('   mpas_horiz_mixing                   = '//trim(mpas_horiz_mixing))
       call mpas_log_write('   mpas_len_disp                       = '//int2str(int(mpas_len_disp)))
       call mpas_log_write('   mpas_visc4_2dsmag                   = '//int2str(int(mpas_visc4_2dsmag)))
       call mpas_log_write('   mpas_del4u_div_factor               = '//int2str(int(mpas_del4u_div_factor)))
       call mpas_log_write('   mpas_w_adv_order                    = '//int2str(mpas_w_adv_order))
       call mpas_log_write('   mpas_theta_adv_order                = '//int2str(mpas_theta_adv_order))
       call mpas_log_write('   mpas_scalar_adv_order               = '//int2str(mpas_scalar_adv_order))
       call mpas_log_write('   mpas_u_vadv_order                   = '//int2str(mpas_u_vadv_order))
       call mpas_log_write('   mpas_w_vadv_order                   = '//int2str(mpas_w_vadv_order))
       call mpas_log_write('   mpas_theta_vadv_order               = '//int2str(mpas_theta_vadv_order))
       call mpas_log_write('   mpas_scalar_vadv_order              = '//int2str(mpas_scalar_vadv_order))
       call mpas_log_write('   mpas_scalar_advection               = '//log2str(mpas_scalar_advection))
       call mpas_log_write('   mpas_positive_definite              = '//log2str(mpas_positive_definite))
       call mpas_log_write('   mpas_monotonic                      = '//log2str(mpas_monotonic))
       call mpas_log_write('   mpas_coef_3rd_order                 = '//int2str(int(mpas_coef_3rd_order)))
       call mpas_log_write('   mpas_smagorinsky_coef               = '//int2str(int(mpas_smagorinsky_coef)))
       call mpas_log_write('   mpas_mix_full                       = '//log2str(mpas_mix_full))
       call mpas_log_write('   mpas_epssm                          = '//int2str(int(mpas_epssm)))
       call mpas_log_write('   mpas_smdiv                          = '//int2str(int(mpas_smdiv)))
       call mpas_log_write('   mpas_apvm_upwinding                 = '//int2str(int(mpas_apvm_upwinding)))
       call mpas_log_write('   mpas_h_ScaleWithMesh                = '//log2str(mpas_h_ScaleWithMesh))
       call mpas_log_write('   mpas_zd                             = '//int2str(int(mpas_zd)))
       call mpas_log_write('   mpas_xnutr                          = '//int2str(int(mpas_xnutr)))
       call mpas_log_write('   mpas_cam_coef                       = '//int2str(int(mpas_cam_coef)))
       call mpas_log_write('   mpas_cam_damping_levels             = '//int2str(mpas_cam_damping_levels))
       call mpas_log_write('   mpas_rayleigh_damp_u                = '//log2str(mpas_rayleigh_damp_u))
       call mpas_log_write('   mpas_rayleigh_damp_u_timescale_days = '//int2str(int(mpas_rayleigh_damp_u_timescale_days)))
       call mpas_log_write('   mpas_number_rayleigh_damp_u_levels  = '//int2str(mpas_number_rayleigh_damp_u_levels))
       call mpas_log_write('   mpas_apply_lbcs                     = '//log2str(mpas_apply_lbcs))
       call mpas_log_write('   mpas_pio_num_iotasks                = '//int2str(mpas_pio_num_iotasks))
       call mpas_log_write('   mpas_pio_stride                     = '//int2str(mpas_pio_stride))
       call mpas_log_write('   mpas_jedi_da                        = '//log2str(mpas_jedi_da))
       call mpas_log_write('   mpas_block_decomp_file_prefix       = '//trim(mpas_block_decomp_file_prefix))
       call mpas_log_write('   mpas_do_restart                     = '//log2str(mpas_do_restart))
       call mpas_log_write('   mpas_print_global_minmax_vel        = '//log2str(mpas_print_global_minmax_vel))
       call mpas_log_write('   mpas_print_detailed_minmax_vel      = '//log2str(mpas_print_detailed_minmax_vel))
       call mpas_log_write('   mpas_print_global_minmax_sca        = '//log2str(mpas_print_global_minmax_sca))
    end if
 end subroutine read_mpas_namelist

end module ufs_mpas_subdriver
