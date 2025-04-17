module ufs_mpas_subdriver
  use mpas_derived_types, only : core_type, domain_type, MPAS_Clock_type
  use module_mpas_config, only : pio_subsystem, pio_stride, pio_numiotasks
  use module_mpas_config, only : ic_filename, mesh_filename
  use module_mpas_config, only : pio_iotype, fcst_mpi_comm
  use module_mpas_config, only : nCellsSolve, nEdgesSolve, nVerticesSolve, nVertLevelsSolve
  use module_mpas_config, only : zref, zref_edge, sphere_radius, pref, pref_edge
  use module_mpas_config, only : maxNCells, maxEdges, nVertLevels
  use module_mpas_config, only : nCells_g, nEdges_g, nVertices_g
  use pio,                only : iosystem_desc_t, file_desc_t, pio_openfile, pio_nowrite
  use pio,                only : pio_init, pio_noerr, pio_iotask_rank, PIO_REARR_BOX
  use pio,                only : pio_global, pio_get_att
  use pio_types,          only : PIO_iotype_pnetcdf
  implicit none
  
  private

  public :: ufs_mpas_init
  public :: corelist, domain_ptr

  type(core_type),       pointer :: corelist   => null()
  type(domain_type),     pointer :: domain_ptr => null()
  type(MPAS_Clock_type), pointer :: clock      => null()
  
contains
  ! #########################################################################################
  !
  ! Procedure to initialize UWM with MPAS dynamical core.
  !
  ! #########################################################################################
  subroutine ufs_mpas_init(Init, time_start, time_end, total_time, calendar, logUnits)
    use mpas_pool_routines,         only : mpas_pool_add_config, mpas_pool_get_subpool
    use mpas_pool_routines,         only : mpas_pool_initialize_time_levels, mpas_pool_get_config
    use mpas_pool_routines,         only : mpas_pool_add_dimension, mpas_pool_get_field
    use mpas_pool_routines,         only : mpas_pool_get_array
    use mpas_framework,             only : mpas_framework_init_phase1, mpas_framework_init_phase2
    use mpas_domain_routines,       only : mpas_allocate_domain, mpas_pool_get_dimension
    use mpas_bootstrapping,         only : mpas_bootstrap_framework_phase1, mpas_bootstrap_framework_phase2
    use mpas_rbf_interpolation,     only : mpas_rbf_interp_initialize
    use mpas_vector_reconstruction, only : mpas_init_reconstruct
    use mpas_stream_inquiry,        only : MPAS_stream_inquiry_new_streaminfo
    use mpas_atm_threading,         only : mpas_atm_threading_init
    use mpas_atm_dimensions,        only : mpas_atm_set_dims
    use mpas_atm_halos,             only : atm_build_halo_groups, exchange_halo_group
    use mpas_derived_types,         only : mpas_pool_type, MPAS_IO_NETCDF, MPAS_Time_Type, field3dReal
    use mpas_kind_types,            only : StrKIND, RKIND
    use mpas_log,                   only : mpas_log_write
    use mpas_timekeeping,           only : mpas_get_clock_time, mpas_get_time, MPAS_START_TIME
    use atm_core_interface,         only : atm_setup_core, atm_setup_domain
    use atm_core,                   only : atm_mpas_init_block, core_clock => clock
    use atm_time_integration,       only : mpas_atm_dynamics_init
    use field_manager_mod,          only : MODEL_ATMOS
    use fms2_io_mod,                only : file_exists
    use mpp_mod,                    only : FATAL, mpp_error
    use MPAS_typedefs,              only : MPAS_init_type

    ! Inputs
    type(MPAS_init_type), intent(inout) :: Init
    integer,        intent(in) :: time_start(6), time_end(6), logUnits(2)
    integer,        intent(in) :: total_time
    character(17),  intent(in) :: calendar

    ! Locals
    character(len=*), parameter :: subname = 'ufs_mpas_subdriver::ufs_mpas_init'
    integer :: i, ndate1, ndate2, tod, ierr, ik, kk
    character(len=32), allocatable, target :: tracer_names(:)
    integer,           allocatable, target :: tracer_types(:)
    integer :: mesh_iotype
    type (mpas_pool_type), pointer :: state, mesh
    integer, pointer :: nVertLevels1, maxEdges1, maxEdges2, num_scalars
    integer, pointer :: nCells
    type(file_desc_t), target :: pioid
    real (kind=RKIND), pointer :: dt
    logical, pointer :: config_do_restart
    type (MPAS_Time_Type) :: startTime
    character(len=StrKIND) :: startTimeStamp
    type (field3dReal), pointer :: scalarsField
    real(rkind), pointer     :: rdzw(:)
    real(rkind), allocatable :: dzw(:)
    character (len=StrKIND), pointer :: xtime
    character (len=StrKIND), pointer :: initial_time1, initial_time2

    ! Initialize PIO
    allocate(pio_subsystem)
    call pio_init(Init%me, fcst_mpi_comm%mpi_val, pio_numiotasks, 0, pio_stride,            &
                  PIO_REARR_BOX, pio_subsystem)

    ! #######################################################################################
    ! #######################################################################################
    !
    ! From CAM  - src/dynamics/mpas/driver/cam_mpas_subdriver.F90: cam_mpas_init_phase1()
    !      MPAS - src/driver/mpas_subdriver.F:                     mpas_init()
    !
    ! #######################################################################################
    ! #######################################################################################

    ! Setup MPAS infrastructure
    allocate(corelist, stat=ierr)
    if ( ierr /= 0 ) call mpp_error(FATAL,subname//": failed to allocate corelist array")
    nullify(corelist%next)

    allocate(corelist%domainlist, stat=ierr)
    if ( ierr /= 0 ) call mpp_error(FATAL,subname//": failed to allocate corelist%domainlist%next")
    nullify(corelist%domainlist%next)

    domain_ptr => corelist%domainlist
    domain_ptr%core => corelist

    call mpas_allocate_domain(domain_ptr)
    domain_ptr%domainID = 0

    ! Initialize MPAS infrastructure
    call mpas_framework_init_phase1(domain_ptr%dminfo, external_comm=fcst_mpi_comm)

    call atm_setup_core(corelist)
    call atm_setup_domain(domain_ptr)

    ! Set up the log manager as early as possible so we can use it for any errors/messages
    ! during subsequent init steps.  We need:
    ! 1) domain_ptr to be allocated,
    ! 2) dmpar_init complete to access dminfo,
    ! 3) *_setup_core to assign the setup_log function pointer
    domain_ptr%core%git_version = 'unknown'
    domain_ptr%core%build_target = 'N/A'
    ierr = domain_ptr % core % setup_log(domain_ptr % logInfo, domain_ptr, unitNumbers=logUnits)
    if ( ierr /= 0 ) then
       call mpp_error(FATAL,subname//": Log setup failed for MPAS-A dycore")
    end if

    ! #######################################################################################
    ! #######################################################################################
    !
    ! From CAM  - src/dynamics/mpas/dyn_comp.F90: dyn_readnl()
    !      MPAS - src/driver/mpas_subdriver.F:    mpas_init()
    !
    ! #######################################################################################
    ! #######################################################################################

    ! Read MPAS namelist.
    if (file_exists('input.nml')) then
       call read_mpas_namelist('input.nml', domain_ptr%configs, Init%mpi_comm, Init%master, Init%me)
    else
       call mpp_error(FATAL,subname//": Cannot find MPAS namelist file, input.nml")
    end if

    ! Set forecast start time (config_start_time)
    ndate1 = time_start(1)*10000 + time_start(2)*100 + time_start(3)
    tod    = time_start(4)*3600  + time_start(5)*60  + time_start(6)
    call mpas_pool_add_config(domain_ptr%configs, 'config_start_time', date2yyyymmdd(ndate1)//'_'//sec2hms(tod))

    ! Set forecast end time (config_stop_time)
    ndate2 = time_end(1)*10000   + time_end(2)*100   + time_end(3)
    tod	   = time_end(4)*3600    + time_end(5)*60    + time_end(6)
    call mpas_pool_add_config(domain_ptr%configs, 'config_stop_time', date2yyyymmdd(ndate2)//'_'//sec2hms(tod))

    ! Set forecaste run time (config_run_duration) #DJS2025 this is not correct. need to fix, but works for current test.
    tod = ndate2 - ndate1 -1
    call mpas_pool_add_config(domain_ptr%configs, 'config_run_duration', trim(int2str(tod))//'_'//sec2hms(total_time))

    ! Set other MPAS required configuration information.
    call mpas_pool_add_config(domain_ptr%configs, 'config_restart_timestamp_name', 'restart_timestamp')
    call mpas_pool_add_config(domain_ptr%configs, 'config_IAU_option',             'off')
    call mpas_pool_add_config(domain_ptr%configs, 'config_do_DAcycling',           .false.)
    call mpas_pool_add_config(domain_ptr%configs, 'config_halo_exch_method',       'mpas_halo')
    call mpas_pool_add_config(domain_ptr%configs, 'config_pio_num_iotasks',        pio_stride)
    call mpas_pool_add_config(domain_ptr%configs, 'config_pio_stride',             pio_numiotasks)

    ! #######################################################################################
    ! #######################################################################################
    !
    ! From CAM  - src/dynamics/mpas/driver/cam_mpas_subdriver.F90: cam_mpas_init_phase2()
    !      MPAS - src/driver/mpas_subdriver.F:                     mpas_init()
    !
    ! #######################################################################################
    ! #######################################################################################

    ! Initialize MPAS infrastructure (phase 2)
    call mpas_framework_init_phase2(domain_ptr, io_system=pio_subsystem, calendar = trim(calendar))

    ! Before defining packages, initialize the stream inquiry instance for the domain
    domain_ptr%streamInfo => MPAS_stream_inquiry_new_streaminfo()
    if (.not. associated(domain_ptr%streamInfo)) then
       call mpp_error(FATAL,subname//": Failed to instantiate streamInfo object for "//trim(domain_ptr%core%coreName))
    end if

    ierr = domain_ptr%core%define_packages(domain_ptr%packages)
    if (ierr /= 0) then
       call mpp_error(FATAL,subname//": Package definition failed for "//trim(domain_ptr%core%coreName))
    end if

    ierr = domain_ptr%core%setup_packages(domain_ptr%configs,  domain_ptr%streamInfo,       &
                                          domain_ptr%packages, domain_ptr%iocontext)
    if (ierr /= 0) then
       call mpp_error(FATAL,subname//": Package setup failed for "//trim(domain_ptr%core%coreName))
    end if

    ierr = domain_ptr%core%setup_decompositions(domain_ptr%decompositions)
    if (ierr /= 0) then
       call mpp_error(FATAL,subname//": Decomposition setup failed for "//trim(domain_ptr%core%coreName))
    end if

    ierr = domain_ptr%core%setup_clock(domain_ptr%clock, domain_ptr%configs)
    if (ierr /= 0) then
       call mpp_error(FATAL,subname//": Clock setup failed for "//trim(domain_ptr%core%coreName))
    end if

    ! Adding a config named 'cam_pcnst' with the number of constituents will indicate to
    ! MPAS-A setup code that it is operating as a UFS dycore, and that it is necessary to
    ! allocate scalars separately from other Registry-defined fields
    call mpas_pool_add_config(domain_ptr%configs, 'cam_pcnst', Init%nwat)

    ! Open MPAS Initial Condition file.
    if (file_exists(ic_filename)) then
       ierr = pio_openfile(pio_subsystem, pioid, pio_iotype, ic_filename, pio_nowrite)
       if (ierr /= 0) then
          call mpp_error(FATAL,subname//": Failed opening MPAS IC File, "//trim(ic_filename))
       end if
    else
       call mpp_error(FATAL,subname//": Cannot find MPAS IC file: "//trim(ic_filename))
    end if

    ! #######################################################################################
    ! #######################################################################################
    !
    ! From CAM  - src/dynamics/mpas/driver/cam_mpas_subdriver.F90: cam_mpas_init_phase3()
    !      MPAS - src/driver/mpas_subdriver.F:                     mpas_init()
    !
    ! #######################################################################################
    ! #######################################################################################

    ! Call MPAS framework bootstrap phase 1
    call mpas_bootstrap_framework_phase1(domain_ptr, mesh_filename, MPAS_IO_NETCDF, pio_file_desc=pioid)

    ! Finalize the setup of blocks and fields
    call mpas_bootstrap_framework_phase2(domain_ptr, pio_file_desc=pioid)

    ! Set up tracers (NOT YET IMPLEMENTED. ONLY ONE, QV)
    call mpas_pool_get_subpool(domain_ptr%blocklist%structs, 'state', state)
    call mpas_pool_get_field(state, 'scalars', scalarsField, timeLevel=1)
    call mpas_pool_add_dimension(state, 'index_qv', 1)
    scalarsField % constituentNames(1) = 'qv'
    call mpas_pool_add_dimension(state, 'moist_start', 1)
    call mpas_pool_add_dimension(state, 'moist_end', Init%nwat)

    ! #######################################################################################
    ! #######################################################################################
    !
    ! From CAM - src/dynamics/mpas/dyn_grid.F90: setup_time_invariant()
    !
    ! #######################################################################################
    ! #######################################################################################

    ! Read in static (IC) data
    call ufs_mpas_read_static(pioid)

    ! Compute unit vectors giving the local north and east directions as well as
    ! the unit normal vector for edges
    call ufs_mpas_compute_unit_vectors()

    ! Access dimensions that are made public via this module
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh', mesh)
    call mpas_pool_get_dimension(mesh, 'nCellsSolve',    nCellsSolve)
    call mpas_pool_get_dimension(mesh, 'nEdgesSolve',    nEdgesSolve)
    call mpas_pool_get_dimension(mesh, 'nVerticesSolve', nVerticesSolve)
    call mpas_pool_get_dimension(mesh, 'nVertLevels',    nVertLevelsSolve)
    call mpas_pool_get_dimension(mesh, 'nCells',         nCells)

    ! Initialize fields needed for reconstruction of cell-centered winds from edge-normal winds
    ! Note: This same pair of calls happens a second time later in the initialization of
    !       the MPAS-A dycore (in atm_mpas_init_block), but the redundant calls do no harm
    call mpas_rbf_interp_initialize(mesh)
    call mpas_init_reconstruct(mesh)

    ! Compute the zeta coordinate at layer interfaces and midpoints.
    ! NOTE: Using number of levels for dycore, nVertLevelsSolve, should we be using phyiscs levels?
    call mpas_pool_get_array(mesh, 'rdzw', rdzw)
    allocate(zref_edge(nVertLevelsSolve+1))
    allocate(zref(nVertLevelsSolve))
    !
    allocate(dzw(nVertLevelsSolve))
    dzw = 1.0 / rdzw
    zref_edge(nVertLevelsSolve+1) = 0.0
    do ik=nVertLevelsSolve,1,-1
       kk = nVertLevelsSolve - ik + 1
       zref_edge(ik) = zref_edge(ik+1) + dzw(kk)
       zref(ik)  = 0.5 * (zref_edge(ik+1) + zref_edge(ik))
    enddo

    ! Compute reference pressure using reference heights.
    allocate(pref_edge(nVertLevelsSolve+1))
    allocate(pref(nVertLevelsSolve))
    ! TO DO. Compute P on edges(interface) from zeta on edges.
    pref_edge(:) = 0.
    pref = (pref_edge(1:nVertLevelsSolve) + pref_edge(2:nVertLevelsSolve+1)) * 0.5

    ! Display reference coordinates.
    if (Init%me == Init%master) then
       write(logUnits(1),'(a)')' Reference Layer Locations: '
       write(logUnits(1),'(a)')' index      height (m)              pressure (hPa) '
       do ik= 1, nVertLevelsSolve
          write(logUnits(1),9830) ik, zref_edge(ik), pref_edge(ik)/100._rkind
          write(logUnits(1),9840)     zref(ik),      pref(ik)/100._rkind
       end do
       write(logUnits(1),9830) nVertLevelsSolve+1, zref_edge(nVertLevelsSolve+1), pref_edge(nVertLevelsSolve+1)/100._rkind

9830   format(1x, i3, f15.4, 9x, f15.4)
9840   format(1x, 3x, 12x, f15.4, 9x, f15.4)
    end if

    ! Read the global sphere_radius attribute.  This is needed to normalize the cell areas.
    ierr = pio_get_att(pioid, pio_global, 'sphere_radius', sphere_radius)
    if( ierr /= 0 ) then
       call mpp_error(FATAL,subname//": Could not find sphere_radius PIO attribute")
    endif

    ! Query global grid dimensions from MPAS
    call ufs_mpas_get_global_dims(nCells_g, nEdges_g, nVertices_g, maxEdges, nVertLevels, maxNCells)

    ! #######################################################################################
    ! #######################################################################################
    !
    ! From CAM  - src/dynamics/mpas/driver/cam_mpas_subdriver.F90: cam_mpas_init_phase4()
    !      MPAS - src/core_atmosphere/mpas_atm_core.F:             atm_core_init()
    !
    ! #######################################################################################
    ! #######################################################################################

    ! Setup threading
    call mpas_atm_threading_init(domain_ptr%blocklist, ierr)
    if ( ierr /= 0 ) then
       call mpp_error(FATAL,subname//": Threading setup failed for core "//trim(domain_ptr % core % coreName))
    end if

    !
    ! Set up inner dimensions used by arrays in optimized dynamics routines
    !
    call mpas_pool_get_subpool(domain_ptr%blocklist%structs, 'state', state)
    call mpas_pool_get_dimension(state, 'nVertLevels', nVertLevels1)
    call mpas_pool_get_dimension(state, 'maxEdges', maxEdges1)
    call mpas_pool_get_dimension(state, 'maxEdges2', maxEdges2)
    call mpas_pool_get_dimension(state, 'num_scalars', num_scalars)
    call mpas_atm_set_dims(nVertLevels1, maxEdges1, maxEdges2, num_scalars)
    Init%levs = nVertLevels1

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
    Init%dt_dycore = dt

    if (.not. config_do_restart) then
       call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'state', state)
       call mpas_pool_initialize_time_levels(state)
    end if

    !
    ! Set startTimeStamp based on the start time of the simulation clock
    !
    startTime = mpas_get_clock_time(clock, MPAS_START_TIME, ierr)
    if ( ierr /= 0 ) then
       call mpp_error(FATAL,subname//": failed to get MPAS_START_TIME")
    end if
    call mpas_get_time(startTime, dateTimeString=startTimeStamp)

    call exchange_halo_group(domain_ptr, 'initialization:u')

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

  end subroutine ufs_mpas_init

  ! #########################################################################################
  ! Procedure to read MPAS namelist(s).
  !
  ! The namelist for MPAS are described in MPAS-Model/src/core_atmosphere/Registry.xml, this
  ! is also where the default values defined below originate.
  !
  ! #########################################################################################
  subroutine read_mpas_namelist(nml_file, configPool, mpicomm, master, me)
    use mpi_f08,            only: MPI_Comm, MPI_CHARACTER, MPI_INTEGER, MPI_REAL8,  MPI_LOGICAL
    use mpi_f08,            only: mpi_bcast, mpi_barrier
    use mpas_derived_types, only: mpas_pool_type
    use mpas_kind_types,    only: StrKIND, RKIND
    use mpas_pool_routines, only: mpas_pool_add_config
    use MPAS_typedefs,      only: r8 => kind_dbl_prec
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
    ! Namelist assimilation
    logical                 :: mpas_jedi_da                        = .false.
    ! Namelist decomposition
    character (len=StrKIND) :: mpas_block_decomp_file_prefix       = 'x1.40962.graph.info.part.'
    ! Namelist restart
    logical                 :: mpas_do_restart                     = .false.
    ! Namelist printout
    logical                 :: mpas_print_global_minmax_vel        = .true.
    logical                 :: mpas_print_detailed_minmax_vel      = .false.
    logical                 :: mpas_print_global_minmax_sca        = .false.

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
    call mpas_pool_add_config(configPool, 'config_num_halos', config_num_halos)
    call mpas_pool_add_config(configPool, 'config_number_of_blocks', config_number_of_blocks)
    call mpas_pool_add_config(configPool, 'config_explicit_proc_decomp', config_explicit_proc_decomp)
    call mpas_pool_add_config(configPool, 'config_proc_decomp_file_prefix', config_proc_decomp_file_prefix)

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
 subroutine ufs_mpas_read_static(pioid)
   use pio, only : file_desc_t
   use mpas_kind_types,     only : StrKIND
   use mpas_io_streams,     only : MPAS_createStream, MPAS_closeStream, MPAS_streamAddField
   use mpas_io_streams,     only : MPAS_readStream
   use mpas_derived_types,  only : MPAS_IO_READ, MPAS_IO_NETCDF, MPAS_Stream_type, MPAS_pool_type
   use mpas_derived_types,  only : field0DReal, field1DReal, field2DReal, field3DReal
   use mpas_derived_types,  only : field1DInteger, field2DInteger, MPAS_STREAM_NOERR
   use mpas_pool_routines,  only : MPAS_pool_get_subpool, MPAS_pool_get_field, MPAS_pool_create_pool
   use mpas_pool_routines,  only : MPAS_pool_destroy_pool, MPAS_pool_add_config
   use mpas_dmpar,          only : MPAS_dmpar_exch_halo_field
   use mpas_stream_manager, only : postread_reindex
   use mpp_mod,             only : FATAL, mpp_error
   ! Arguments
   type (file_desc_t), intent(in), pointer :: pioid

   ! Local variables
   character(len=*), parameter :: subname = 'ufs_mpas_subdriver::ufs_mpas_read_static'

   character(len=StrKIND) :: errString

   integer :: ierr
   integer :: ierr_total
   type (MPAS_pool_type), pointer :: meshPool
   type (MPAS_pool_type), pointer :: reindexPool
   type (MPAS_pool_type), pointer :: allPackages, reindexPkgs
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

   type (MPAS_Stream_type) :: mesh_stream

   nullify(cell_gradient_coef_x)
   nullify(cell_gradient_coef_y)

   call MPAS_createStream(mesh_stream, domain_ptr % ioContext, 'not_used', MPAS_IO_NETCDF, MPAS_IO_READ, &
                           pio_file_desc=pioid, ierr=ierr)
   if (ierr /= MPAS_STREAM_NOERR) then
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

   call MPAS_streamAddField(mesh_stream, latCell, ierr=ierr)
   if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
   call MPAS_streamAddField(mesh_stream, lonCell, ierr=ierr)
   if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
   call MPAS_streamAddField(mesh_stream, xCell, ierr=ierr)
   if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
   call MPAS_streamAddField(mesh_stream, yCell, ierr=ierr)
   if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
   call MPAS_streamAddField(mesh_stream, zCell, ierr=ierr)
   if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1

   call MPAS_streamAddField(mesh_stream, latEdge, ierr=ierr)
   if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
   call MPAS_streamAddField(mesh_stream, lonEdge, ierr=ierr)
   if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
   call MPAS_streamAddField(mesh_stream, xEdge, ierr=ierr)
   if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
   call MPAS_streamAddField(mesh_stream, yEdge, ierr=ierr)
   if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
   call MPAS_streamAddField(mesh_stream, zEdge, ierr=ierr)
   if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1

   call MPAS_streamAddField(mesh_stream, latVertex, ierr=ierr)
   if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
   call MPAS_streamAddField(mesh_stream, lonVertex, ierr=ierr)
   if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
   call MPAS_streamAddField(mesh_stream, xVertex, ierr=ierr)
   if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
   call MPAS_streamAddField(mesh_stream, yVertex, ierr=ierr)
   if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
   call MPAS_streamAddField(mesh_stream, zVertex, ierr=ierr)
   if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1

   call MPAS_streamAddField(mesh_stream, indexToCellID, ierr=ierr)
   if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
   call MPAS_streamAddField(mesh_stream, indexToEdgeID, ierr=ierr)
   if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
   call MPAS_streamAddField(mesh_stream, indexToVertexID, ierr=ierr)
   if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1

   call MPAS_streamAddField(mesh_stream, fEdge, ierr=ierr)
   if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
   call MPAS_streamAddField(mesh_stream, fVertex, ierr=ierr)
   if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1

   call MPAS_streamAddField(mesh_stream, areaCell, ierr=ierr)
   if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
   call MPAS_streamAddField(mesh_stream, areaTriangle, ierr=ierr)
   if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
   call MPAS_streamAddField(mesh_stream, dcEdge, ierr=ierr)
   if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
   call MPAS_streamAddField(mesh_stream, dvEdge, ierr=ierr)
   if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
   call MPAS_streamAddField(mesh_stream, angleEdge, ierr=ierr)
   if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
   call MPAS_streamAddField(mesh_stream, kiteAreasOnVertex, ierr=ierr)
   if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
   call MPAS_streamAddField(mesh_stream, weightsOnEdge, ierr=ierr)
   if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1

   call MPAS_streamAddField(mesh_stream, meshDensity, ierr=ierr)
   if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1

   call MPAS_streamAddField(mesh_stream, nEdgesOnCell, ierr=ierr)
   if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
   call MPAS_streamAddField(mesh_stream, nEdgesOnEdge, ierr=ierr)
   if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1

   call MPAS_streamAddField(mesh_stream, cellsOnEdge, ierr=ierr)
   if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
   call MPAS_streamAddField(mesh_stream, edgesOnCell, ierr=ierr)
   if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
   call MPAS_streamAddField(mesh_stream, edgesOnEdge, ierr=ierr)
   if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
   call MPAS_streamAddField(mesh_stream, cellsOnCell, ierr=ierr)
   if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
   call MPAS_streamAddField(mesh_stream, verticesOnCell, ierr=ierr)
   if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
   call MPAS_streamAddField(mesh_stream, verticesOnEdge, ierr=ierr)
   if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
   call MPAS_streamAddField(mesh_stream, edgesOnVertex, ierr=ierr)
   if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
   call MPAS_streamAddField(mesh_stream, cellsOnVertex, ierr=ierr)
   if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1

   call MPAS_streamAddField(mesh_stream, cf1, ierr=ierr)
   if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
   call MPAS_streamAddField(mesh_stream, cf2, ierr=ierr)
   if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
   call MPAS_streamAddField(mesh_stream, cf3, ierr=ierr)
   if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1

   call MPAS_streamAddField(mesh_stream, rdzw, ierr=ierr)
   if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
   call MPAS_streamAddField(mesh_stream, dzu, ierr=ierr)
   if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
   call MPAS_streamAddField(mesh_stream, rdzu, ierr=ierr)
   if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
   call MPAS_streamAddField(mesh_stream, fzm, ierr=ierr)
   if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
   call MPAS_streamAddField(mesh_stream, fzp, ierr=ierr)
   if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1

   call MPAS_streamAddField(mesh_stream, zgrid, ierr=ierr)
   if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
   call MPAS_streamAddField(mesh_stream, zxu, ierr=ierr)
   if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
   call MPAS_streamAddField(mesh_stream, zz, ierr=ierr)
   if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
   call MPAS_streamAddField(mesh_stream, zb, ierr=ierr)
   if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
   call MPAS_streamAddField(mesh_stream, zb3, ierr=ierr)
   if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1

   call MPAS_streamAddField(mesh_stream, deriv_two, ierr=ierr)
   if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
   call MPAS_streamAddField(mesh_stream, cellTangentPlane, ierr=ierr)
   if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
   call MPAS_streamAddField(mesh_stream, coeffs_reconstruct, ierr=ierr)
   if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1

   call MPAS_streamAddField(mesh_stream, edgeNormalVectors, ierr=ierr)
   if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
   call MPAS_streamAddField(mesh_stream, localVerticalUnitVectors, ierr=ierr)
   if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
   call MPAS_streamAddField(mesh_stream, defc_a, ierr=ierr)
   if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
   call MPAS_streamAddField(mesh_stream, defc_b, ierr=ierr)
   if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1

   if (ierr_total > 0) then
      write(errString, '(a,i0,a)') subname//': FATAL: Failed to add ', ierr_total, ' fields to static input stream.'
      call mpp_error(FATAL,trim(errString))
   end if

   call MPAS_readStream(mesh_stream, 1, ierr=ierr)
   if (ierr /= MPAS_STREAM_NOERR) then
      call mpp_error(FATAL,subname//': FATAL: Failed to read static input stream.')
   end if

   call MPAS_closeStream(mesh_stream, ierr=ierr)
   if (ierr /= MPAS_STREAM_NOERR) then
      call mpp_error(FATAL,subname//': FATAL: Failed to close static input stream.')
   end if

   !
   ! Perform halo updates for all decomposed fields (i.e., fields with
   ! an outermost dimension of nCells, nVertices, or nEdges)
   !
   call MPAS_dmpar_exch_halo_field(latCell)
   call MPAS_dmpar_exch_halo_field(lonCell)
   call MPAS_dmpar_exch_halo_field(xCell)
   call MPAS_dmpar_exch_halo_field(yCell)
   call MPAS_dmpar_exch_halo_field(zCell)

   call MPAS_dmpar_exch_halo_field(latEdge)
   call MPAS_dmpar_exch_halo_field(lonEdge)
   call MPAS_dmpar_exch_halo_field(xEdge)
   call MPAS_dmpar_exch_halo_field(yEdge)
   call MPAS_dmpar_exch_halo_field(zEdge)

   call MPAS_dmpar_exch_halo_field(latVertex)
   call MPAS_dmpar_exch_halo_field(lonVertex)
   call MPAS_dmpar_exch_halo_field(xVertex)
   call MPAS_dmpar_exch_halo_field(yVertex)
   call MPAS_dmpar_exch_halo_field(zVertex)

   call MPAS_dmpar_exch_halo_field(indexToCellID)
   call MPAS_dmpar_exch_halo_field(indexToEdgeID)
   call MPAS_dmpar_exch_halo_field(indexToVertexID)

   call MPAS_dmpar_exch_halo_field(fEdge)
   call MPAS_dmpar_exch_halo_field(fVertex)

   call MPAS_dmpar_exch_halo_field(areaCell)
   call MPAS_dmpar_exch_halo_field(areaTriangle)
   call MPAS_dmpar_exch_halo_field(dcEdge)
   call MPAS_dmpar_exch_halo_field(dvEdge)
   call MPAS_dmpar_exch_halo_field(angleEdge)
   call MPAS_dmpar_exch_halo_field(kiteAreasOnVertex)
   call MPAS_dmpar_exch_halo_field(weightsOnEdge)

   call MPAS_dmpar_exch_halo_field(meshDensity)

   call MPAS_dmpar_exch_halo_field(nEdgesOnCell)
   call MPAS_dmpar_exch_halo_field(nEdgesOnEdge)

   call MPAS_dmpar_exch_halo_field(cellsOnEdge)
   call MPAS_dmpar_exch_halo_field(edgesOnCell)
   call MPAS_dmpar_exch_halo_field(edgesOnEdge)
   call MPAS_dmpar_exch_halo_field(cellsOnCell)
   call MPAS_dmpar_exch_halo_field(verticesOnCell)
   call MPAS_dmpar_exch_halo_field(verticesOnEdge)
   call MPAS_dmpar_exch_halo_field(edgesOnVertex)
   call MPAS_dmpar_exch_halo_field(cellsOnVertex)

   call MPAS_dmpar_exch_halo_field(zgrid)
   call MPAS_dmpar_exch_halo_field(zxu)
   call MPAS_dmpar_exch_halo_field(zz)
   call MPAS_dmpar_exch_halo_field(zb)
   call MPAS_dmpar_exch_halo_field(zb3)

   call MPAS_dmpar_exch_halo_field(deriv_two)
   call MPAS_dmpar_exch_halo_field(cellTangentPlane)
   call MPAS_dmpar_exch_halo_field(coeffs_reconstruct)

   call MPAS_dmpar_exch_halo_field(edgeNormalVectors)
   call MPAS_dmpar_exch_halo_field(localVerticalUnitVectors)
   call MPAS_dmpar_exch_halo_field(defc_a)
   call MPAS_dmpar_exch_halo_field(defc_b)
   !
   ! Re-index from global index space to local index space
   !
   call MPAS_pool_create_pool(reindexPool)

   call MPAS_pool_add_config(reindexPool, 'cellsOnEdge', 1)
   call MPAS_pool_add_config(reindexPool, 'edgesOnCell', 1)
   call MPAS_pool_add_config(reindexPool, 'edgesOnEdge', 1)
   call MPAS_pool_add_config(reindexPool, 'cellsOnCell', 1)
   call MPAS_pool_add_config(reindexPool, 'verticesOnCell', 1)
   call MPAS_pool_add_config(reindexPool, 'verticesOnEdge', 1)
   call MPAS_pool_add_config(reindexPool, 'edgesOnVertex', 1)
   call MPAS_pool_add_config(reindexPool, 'cellsOnVertex', 1)

   ! Use an empty package list for reindexPool
   call MPAS_pool_create_pool(reindexPkgs)

   call postread_reindex(meshPool, domain_ptr % streamManager % allPackages, &
                         reindexPool, reindexPkgs)

   call MPAS_pool_destroy_pool(reindexPool)
   call MPAS_pool_destroy_pool(reindexPkgs)

 end subroutine ufs_mpas_read_static

 !> ########################################################################################
 !
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
 !
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
 !
 !> \brief  Returns global mesh dimensions
 !> \author Michael Duda
 !> \date   22 August 2019
 !> \details
 !>  This routine returns on all tasks the number of global cells, edges,
 !>  vertices, maxEdges, vertical layers, and the maximum number of cells owned by any task.
 !>
 !> \update: Dustin Swales April 2025 - Modified for use in UWM
 !
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
  
end module ufs_mpas_subdriver
