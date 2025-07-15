!> ###########################################################################################
!> \file ufs_mpas_subdriver.F90
!> UFSATM subdriver for MPAS dynamical core.
!>
!> Routines from the subdrivers for MPAS-A and CAM-SIMA have been adopted/modified here for use
!> within the UFS Weather Model.
!> MPAS-A Subdriver:    MPAS-Model/src/driver/mpas_subdriver.F
!> CAM-SIMA (external): src/dynamics/mpas/driver/dyn_mpas_subdriver.F90
!>                      (https://github.com/ESCOMP/CAM-SIMA/blob/development/)
!>
!> Overview:
!> Initialization is broken down into two phases, with ufs_mpas_define_scalars() called in
!> between:
!> ufs_mpas_init_phase1:    Initialize MPAS framework, Read in namelist, Read static data.
!> ufs_mpas_define_scalars: Set up scalars/tracers/constituents/... 
!> ufs_mpas_init_phase2:    Complete MPAS initialization
!>
!> Forward integration of the dycore is handled in ufs_mpas_run. The current forecast time,
!> forecast interval, and MPAS dycore time step are used to integrate the model forward in
!> time. Afterwards, atm_compute_output_diagnostics() is called to compute fields needed by
!> the Physics.
!>
!> Other public routines used the UFSATM driver
!> ufs_mpas_open_init:      Open MPAS Initial Condition file, return PIO file handle.
!>
!> ###########################################################################################
module ufs_mpas_subdriver
  use mpi_f08
  use mpas_derived_types, only : core_type, domain_type, mpas_Clock_type
  use mpas_kind_types,    only : StrKIND, rkind
  use module_mpas_config, only : pio_subsystem, pio_stride, pio_numiotasks, pio_iodesc
  use module_mpas_config, only : ic_filename, lbc_filename
  use module_mpas_config, only : pio_iotype, fcst_mpi_comm, pioid
  use module_mpas_config, only : zref, zref_edge, sphere_radius, pref, pref_edge
  use module_mpas_config, only : maxNCells, maxEdges, nVertLevels
  use module_mpas_config, only : nCellsGlobal, nEdgesGlobal, nVerticesGlobal
  use module_mpas_config, only : nCellsSolve, nEdgesSolve, nVerticesSolve, nVertLevelsSolve
  use module_mpas_config, only : dt_atmos, n_atmos
  use module_mpas_config, only : latCellGlobal, lonCellGlobal, areaCellGlobal
  implicit none
  
  private

  public :: MPAS_control_type
  public :: ufs_mpas_init_phase1
  public :: ufs_mpas_define_scalars
  public :: ufs_mpas_init_phase2
  public :: ufs_mpas_run
  public :: ufs_mpas_open_init
  public :: corelist, domain_ptr
  public :: constituent_name
  public :: is_water_species
  public :: dyn_mpas_read_write_stream

  !> #########################################################################################
  !>
  !> #########################################################################################
  type MPAS_control_type

     ! Namelist filename
     character(len=64) :: fn_nml

     ! Full namelist for use with internal file reads
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

  !> #########################################################################################
  !
  !> #########################################################################################
  type :: var_info_type
     private
     character(64) :: name = ''
     character(10) :: type = ''
     integer :: rank = 0
  end type var_info_type

  !> #########################################################################################
  !> This list corresponds to the "invariant" stream in MPAS registry.
  !> It consists of variables that are members of the "mesh" struct.
  !> #########################################################################################
  type(var_info_type), parameter :: invariant_var_info_list(*) = [ &
       var_info_type('angleEdge'                       , 'real'      , 1), &
       var_info_type('areaCell'                        , 'real'      , 1), &
       var_info_type('areaTriangle'                    , 'real'      , 1), &
       var_info_type('bdyMaskCell'                     , 'integer'   , 1), &
       var_info_type('bdyMaskEdge'                     , 'integer'   , 1), &
       var_info_type('bdyMaskVertex'                   , 'integer'   , 1), &
       var_info_type('cellTangentPlane'                , 'real'      , 3), &
       var_info_type('cell_gradient_coef_x'            , 'real'      , 2), &
       var_info_type('cell_gradient_coef_y'            , 'real'      , 2), &
       var_info_type('cellsOnCell'                     , 'integer'   , 2), &
       var_info_type('cellsOnEdge'                     , 'integer'   , 2), &
       var_info_type('cellsOnVertex'                   , 'integer'   , 2), &
       var_info_type('cf1'                             , 'real'      , 0), &
       var_info_type('cf2'                             , 'real'      , 0), &
       var_info_type('cf3'                             , 'real'      , 0), &
       var_info_type('coeffs_reconstruct'              , 'real'      , 3), &
       var_info_type('dcEdge'                          , 'real'      , 1), &
       var_info_type('defc_a'                          , 'real'      , 2), &
       var_info_type('defc_b'                          , 'real'      , 2), &
       var_info_type('deriv_two'                       , 'real'      , 3), &
       var_info_type('dss'                             , 'real'      , 2), &
       var_info_type('dvEdge'                          , 'real'      , 1), &
       var_info_type('dzu'                             , 'real'      , 1), &
       var_info_type('edgeNormalVectors'               , 'real'      , 2), &
       var_info_type('edgesOnCell'                     , 'integer'   , 2), &
       var_info_type('edgesOnEdge'                     , 'integer'   , 2), &
       var_info_type('edgesOnVertex'                   , 'integer'   , 2), &
       var_info_type('fEdge'                           , 'real'      , 1), &
       var_info_type('fVertex'                         , 'real'      , 1), &
       var_info_type('fzm'                             , 'real'      , 1), &
       var_info_type('fzp'                             , 'real'      , 1), &
       var_info_type('indexToCellID'                   , 'integer'   , 1), &
       var_info_type('indexToEdgeID'                   , 'integer'   , 1), &
       var_info_type('indexToVertexID'                 , 'integer'   , 1), &
       var_info_type('kiteAreasOnVertex'               , 'real'      , 2), &
       var_info_type('latCell'                         , 'real'      , 1), &
       var_info_type('latEdge'                         , 'real'      , 1), &
       var_info_type('latVertex'                       , 'real'      , 1), &
       var_info_type('localVerticalUnitVectors'        , 'real'      , 2), &
       var_info_type('lonCell'                         , 'real'      , 1), &
       var_info_type('lonEdge'                         , 'real'      , 1), &
       var_info_type('lonVertex'                       , 'real'      , 1), &
       var_info_type('meshDensity'                     , 'real'      , 1), &
       var_info_type('nEdgesOnCell'                    , 'integer'   , 1), &
       var_info_type('nEdgesOnEdge'                    , 'integer'   , 1), &
       var_info_type('nominalMinDc'                    , 'real'      , 0), &
       var_info_type('qv_init'                         , 'real'      , 1), &
       var_info_type('rdzu'                            , 'real'      , 1), &
       var_info_type('rdzw'                            , 'real'      , 1), &
       var_info_type('t_init'                          , 'real'      , 2), &
       var_info_type('u_init'                          , 'real'      , 1), &
       var_info_type('v_init'                          , 'real'      , 1), &
       var_info_type('verticesOnCell'                  , 'integer'   , 2), &
       var_info_type('verticesOnEdge'                  , 'integer'   , 2), &
       var_info_type('weightsOnEdge'                   , 'real'      , 2), &
       var_info_type('xCell'                           , 'real'      , 1), &
       var_info_type('xEdge'                           , 'real'      , 1), &
       var_info_type('xVertex'                         , 'real'      , 1), &
       var_info_type('yCell'                           , 'real'      , 1), &
       var_info_type('yEdge'                           , 'real'      , 1), &
       var_info_type('yVertex'                         , 'real'      , 1), &
       var_info_type('zCell'                           , 'real'      , 1), &
       var_info_type('zEdge'                           , 'real'      , 1), &
       var_info_type('zVertex'                         , 'real'      , 1), &
       var_info_type('zb'                              , 'real'      , 3), &
       var_info_type('zb3'                             , 'real'      , 3), &
       var_info_type('zgrid'                           , 'real'      , 2), &
       var_info_type('zxu'                             , 'real'      , 2), &
       var_info_type('zz'                              , 'real'      , 2)  &
    ]

  ! Whether a variable should be in input or restart can be determined by looking at
  ! the `atm_init_coupled_diagnostics` subroutine in MPAS.
  ! If a variable first appears on the LHS of an equation, it should be in restart.
  ! If a variable first appears on the RHS of an equation, it should be in input.
  ! The remaining ones of interest should be in output.

  !> #########################################################################################
  !> This list corresponds to the "input" stream in MPAS registry.
  !> It consists of variables that are members of the "diag" and "state" struct.
  !> Only variables that are specific to the "input" stream are included.
  !> #########################################################################################
  type(var_info_type), parameter :: input_var_info_list(*) = [ &
       var_info_type('Time'                            , 'real'      , 0), &
       var_info_type('initial_time'                    , 'character' , 0), &
       var_info_type('rho'                             , 'real'      , 2), &
       var_info_type('rho_base'                        , 'real'      , 2), &
       var_info_type('scalars'                         , 'real'      , 3), &
       var_info_type('theta'                           , 'real'      , 2), &
       var_info_type('theta_base'                      , 'real'      , 2), &
       var_info_type('u'                               , 'real'      , 2), &
       var_info_type('w'                               , 'real'      , 2), &
       var_info_type('xtime'                           , 'character' , 0)  &
    ]

  !> #########################################################################################
  !> This list corresponds to the "restart" stream in MPAS registry.
  !> It consists of variables that are members of the "diag" and "state" struct.
  !> Only variables that are specific to the "restart" stream are included.
  !> #########################################################################################
  type(var_info_type), parameter :: restart_var_info_list(*) = [ &
       var_info_type('exner'                           , 'real'      , 2), &
       var_info_type('exner_base'                      , 'real'      , 2), &
       var_info_type('pressure_base'                   , 'real'      , 2), &
       var_info_type('pressure_p'                      , 'real'      , 2), &
       var_info_type('rho_p'                           , 'real'      , 2), &
       var_info_type('rho_zz'                          , 'real'      , 2), &
       var_info_type('rtheta_base'                     , 'real'      , 2), &
       var_info_type('rtheta_p'                        , 'real'      , 2), &
       var_info_type('ru'                              , 'real'      , 2), &
       var_info_type('ru_p'                            , 'real'      , 2), &
       var_info_type('rw'                              , 'real'      , 2), &
       var_info_type('rw_p'                            , 'real'      , 2), &
       var_info_type('theta_m'                         , 'real'      , 2)  &
    ]

  !> #########################################################################################
  !> This list corresponds to the "output" stream in MPAS registry.
  !> It consists of variables that are members of the "diag" struct.
  !> Only variables that are specific to the "output" stream are included.
  !> #########################################################################################
  type(var_info_type), parameter :: output_var_info_list(*) = [ &
       var_info_type('divergence'                      , 'real'      , 2), &
       var_info_type('pressure'                        , 'real'      , 2), &
       var_info_type('relhum'                          , 'real'      , 2), &
       var_info_type('surface_pressure'                , 'real'      , 1), &
       var_info_type('uReconstructMeridional'          , 'real'      , 2), &
       var_info_type('uReconstructZonal'               , 'real'      , 2), &
       var_info_type('vorticity'                       , 'real'      , 2)  &
    ]

  !> #########################################################################################
  !>
  !> #########################################################################################
  type(core_type),       pointer :: corelist   => null()
  type(domain_type),     pointer :: domain_ptr => null()
  type(mpas_Clock_type), pointer :: clock      => null()
  
  character(StrKIND), allocatable :: constituent_name(:)
  integer, allocatable :: index_constituent_to_mpas_scalar(:)
  integer, allocatable :: index_mpas_scalar_to_constituent(:)
  logical, allocatable :: is_water_species(:)  
contains
  !> #########################################################################################
  !> Convert one or more values of any intrinsic data types to a character string for pretty
  !> printing.
  !> If `value` contains more than one element, the elements will be stringified, delimited by `separator`, then concatenated.
  !> If `value` contains exactly one element, the element will be stringified without using `separator`.
  !> If `value` contains zero element or is of unsupported data types, an empty character string is produced.
  !> If `separator` is not supplied, it defaults to ", " (i.e., a comma and a space).
  !> (KCW, 2024-02-04)
  !> Ported for UWM (DJS: 2025)
  !> #########################################################################################
  pure function stringify(value, separator)
    use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64

    class(*), intent(in) :: value(:)
    character(*), optional, intent(in) :: separator
    character(:), allocatable :: stringify

    integer, parameter :: sizelimit = 1024

    character(:), allocatable :: buffer, delimiter, format
    character(:), allocatable :: value_c(:)
    integer :: i, n, offset

    if (present(separator)) then
       delimiter = separator
    else
       delimiter = ', '
    end if

    n = min(size(value), sizelimit)

    if (n == 0) then
       stringify = ''

       return
    end if

    select type (value)
    type is (character(*))
       allocate(character(len(value) * n + len(delimiter) * (n - 1)) :: buffer)

       buffer(:) = ''
       offset = 0

       ! Workaround for a bug in GNU Fortran >= 12. This is perhaps the manifestation of GCC Bugzilla Bug 100819.
       ! When a character string array is passed as the actual argument to an unlimited polymorphic dummy argument,
       ! its array index and length parameter are mishandled.
       allocate(character(len(value)) :: value_c(size(value)))

       value_c(:) = value(:)

       do i = 1, n
          if (len(delimiter) > 0 .and. i > 1) then
             buffer(offset + 1:offset + len(delimiter)) = delimiter
             offset = offset + len(delimiter)
          end if

          if (len_trim(adjustl(value_c(i))) > 0) then
             buffer(offset + 1:offset + len_trim(adjustl(value_c(i)))) = trim(adjustl(value_c(i)))
             offset = offset + len_trim(adjustl(value_c(i)))
          end if
       end do

       deallocate(value_c)
    type is (integer(int32))
       allocate(character(11 * n + len(delimiter) * (n - 1)) :: buffer)
       allocate(character(17 + len(delimiter) + floor(log10(real(n))) + 1) :: format)

       write(format, '(a, i0, 3a)') '(ss, ', n, '(i0, :, "', delimiter, '"))'
       write(buffer, format) value
    type is (integer(int64))
       allocate(character(20 * n + len(delimiter) * (n - 1)) :: buffer)
       allocate(character(17 + len(delimiter) + floor(log10(real(n))) + 1) :: format)

       write(format, '(a, i0, 3a)') '(ss, ', n, '(i0, :, "', delimiter, '"))'
       write(buffer, format) value
    type is (logical)
       allocate(character(1 * n + len(delimiter) * (n - 1)) :: buffer)
       allocate(character(13 + len(delimiter) + floor(log10(real(n))) + 1) :: format)

       write(format, '(a, i0, 3a)') '(', n, '(l1, :, "', delimiter, '"))'
       write(buffer, format) value
    type is (real(real32))
       allocate(character(13 * n + len(delimiter) * (n - 1)) :: buffer)

       if (maxval(abs(value)) < 1.0e5_real32) then
          allocate(character(20 + len(delimiter) + floor(log10(real(n))) + 1) :: format)
          write(format, '(a, i0, 3a)') '(ss, ', n, '(f13.6, :, "', delimiter, '"))'
       else
          allocate(character(23 + len(delimiter) + floor(log10(real(n))) + 1) :: format)
          write(format, '(a, i0, 3a)') '(ss, ', n, '(es13.6e2, :, "', delimiter, '"))'
       end if

       write(buffer, format) value
    type is (real(real64))
       allocate(character(13 * n + len(delimiter) * (n - 1)) :: buffer)

       if (maxval(abs(value)) < 1.0e5_real64) then
          allocate(character(20 + len(delimiter) + floor(log10(real(n))) + 1) :: format)
          write(format, '(a, i0, 3a)') '(ss, ', n, '(f13.6, :, "', delimiter, '"))'
       else
          allocate(character(23 + len(delimiter) + floor(log10(real(n))) + 1) :: format)
          write(format, '(a, i0, 3a)') '(ss, ', n, '(es13.6e2, :, "', delimiter, '"))'
       end if

       write(buffer, format) value
    class default
       stringify = ''

       return
    end select

    stringify = trim(buffer)
  end function stringify
  
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
    use mpas_constants,             only : mpas_constants_compute_derived, pi => pii
    use mpas_attlist,               only : mpas_add_att
    ! FMS
    use field_manager_mod,          only : MODEL_ATMOS
    use fms2_io_mod,                only : file_exists
    use mpp_mod,                    only : FATAL, mpp_error
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
    integer, pointer :: num_scalars

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
    call mpas_log_write('config_start_time = '//date2yyyymmdd(ndate1)//'_'//sec2hms(tod))

    ! Set forecast end time (config_stop_time)
    ndate2 = time_end(1)*10000   + time_end(2)*100   + time_end(3)
    tod	   = time_end(4)*3600    + time_end(5)*60    + time_end(6)
    call mpas_pool_add_config(domain_ptr % configs, 'config_stop_time', date2yyyymmdd(ndate2)//'_'//sec2hms(tod))
    call mpas_log_write('config_stop_time  = '//date2yyyymmdd(ndate2)//'_'//sec2hms(tod))

    ! Set forecaste run time (config_run_duration) #DJS2025 this is not correct. need to fix, but works for current test.
    tod = max(ndate2 - ndate1 - 1,0)
    call mpas_pool_add_config(domain_ptr % configs, 'config_run_duration', trim(int2str(tod))//'_'//sec2hms(total_time))
    call mpas_log_write('config_run_duration = '//trim(int2str(tod))//'_'//sec2hms(total_time))
    
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
    call mpas_pool_add_config(domain_ptr % configs, 'cam_pcnst', Cfg % nConstituents)

    ! Call MPAS framework bootstrap phase 1
    call mpas_bootstrap_framework_phase1(domain_ptr, "external mesh file", mpas_IO_NETCDF, pio_file_desc=pioid)

    ! Finalize the setup of blocks and fields
    call mpas_bootstrap_framework_phase2(domain_ptr, pio_file_desc=pioid)
    
    ! Add num_scalars from "state" pool to "dimensions".
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'state', state)
    call mpas_pool_get_dimension(state, 'num_scalars', num_scalars)
    call mpas_pool_add_dimension(domain_ptr % blocklist % dimensions, 'num_scalars', num_scalars)
    nullify(num_scalars)
    call mpas_pool_add_dimension(state, 'index_qv', 1)
    call mpas_pool_add_dimension(state, 'moist_start', 1)
    call mpas_pool_add_dimension(state, 'moist_end', Cfg % nwat)

    ! Read in static (invariant) data
    call dyn_mpas_read_write_stream( 'r', 'invariant')

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

    ! Set MPAS mesh lon/lat/area.
    allocate(latCellGlobal(nCellsGlobal), lonCellGlobal(nCellsGlobal), areaCellGlobal(nCellsGlobal))
    call ufs_mpas_get_global_coords(latCellGlobal, lonCellGlobal, areaCellGlobal)

  end subroutine ufs_mpas_init_phase1

  !> ########################################################################################
  !> Procedure to initialize UWM with MPAS dynamical core.
  !> 
  !> ########################################################################################
  subroutine ufs_mpas_init_phase2(Cfg)
    use mpas_kind_types,            only : StrKIND, RKIND
    use mpas_derived_types,         only : mpas_pool_type, mpas_Time_Type, field0DReal, field2dreal
    use mpas_domain_routines,       only : mpas_pool_get_dimension
    use mpas_pool_routines,         only : mpas_pool_get_subpool
    use mpas_pool_routines,         only : mpas_pool_initialize_time_levels, mpas_pool_get_config
    use mpas_pool_routines,         only : mpas_pool_get_array, mpas_pool_get_field
    use mpas_atm_dimensions,        only : mpas_atm_set_dims
    use mpas_atm_threading,         only : mpas_atm_threading_init
    use mpp_mod,                    only : FATAL, mpp_error
    use mpas_atm_halos,             only : atm_build_halo_groups, exchange_halo_group
    use atm_core,                   only : atm_mpas_init_block, core_clock => clock
    use atm_time_integration,       only : mpas_atm_dynamics_init
    use mpas_timekeeping,           only : mpas_get_clock_time, mpas_get_time, mpas_START_TIME
    use mpas_log,                   only : mpas_log_write
    use mpas_attlist,               only : mpas_modify_att
    use mpas_string_utils,          only : mpas_string_replace
    use mpas_field_routines,        only : mpas_allocate_scratch_field
    ! Arguments
    type(mpas_control_type), intent(inout) :: Cfg
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
    type(field0dreal), pointer :: field_0d_real
    type(field2dreal), pointer :: field_2d_real

    !
    ! Setup threading
    !
    call mpas_log_write('Setting up OpenMP threading')
    call mpas_atm_threading_init(domain_ptr%blocklist, ierr)
    if ( ierr /= 0 ) then
       call mpp_error(FATAL,subname//": Threading setup failed for core "//trim(domain_ptr % core % coreName))
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
    Cfg % levs = nVertLevels1

    !
    ! Set "local" clock to point to the clock contained in the domain type
    !
    clock => domain_ptr % clock
    core_clock => domain_ptr % clock

    !
    ! Build halo exchange groups and set method for exchanging halos in a group
    !
    call mpas_log_write('Building halo exchange groups.')

    nullify(exchange_halo_group)
    call atm_build_halo_groups(domain_ptr, ierr)

    if (ierr /= 0) then
       call mpp_error(FATAL,subname//": failed to build MPAS-A halo exchange groups.")
    end if
    if (.not. associated(exchange_halo_group)) then
       call mpp_error(FATAL,subname//": failed to build MPAS-A halo exchange groups.")
    end if

    ! Variables in MPAS "state" pool have more than one time level. Copy the values from the first time level of
    ! such variables into all subsequent time levels to initialize them.
    call mpas_pool_get_config(domain_ptr % blocklist % configs, 'config_do_restart', config_do_restart)
    call mpas_pool_get_config(domain_ptr % blocklist % configs, 'config_dt', dt)

    if (.not. config_do_restart) then
       call mpas_log_write('Initializing time levels')
       call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'state', state)
       call mpas_pool_initialize_time_levels(state)
       nullify(state)
    end if
    nullify (config_do_restart)

    call exchange_halo_group(domain_ptr, 'initialization:u',ierr=ierr)
    if ( ierr /= 0 ) then
       call mpp_error(FATAL,subname//'Failed to exchange halo layers for group "initialization:u"')
    end if

    call mpas_log_write('Initializing atmospheric variables')
    
    ! How many calls to MPAS dycore for each ATMosphere time step?
    Cfg%dt_dycore = dt
    n_atmos = dt_atmos/dt
    
    !
    ! Set startTimeStamp based on the start time of the simulation clock
    !
    startTime = mpas_get_clock_time(clock, mpas_START_TIME, ierr)
    if ( ierr /= 0 ) then
       call mpp_error(FATAL,subname//': Failed to get clock_time "mpas_START_TIME"')
    end if
    call mpas_get_time(startTime, dateTimeString=startTimeStamp, ierr=ierr)
    if ( ierr /= 0 ) then
       call mpp_error(FATAL,subname//': Failed to get time mpas_START_TIME"')
    end if

    
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh', mesh)
    !call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'state', state)

    call atm_mpas_init_block(domain_ptr % dminfo, domain_ptr % streamManager, domain_ptr % blocklist, mesh, dt)

    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'state', state)
    call mpas_pool_get_array(state, 'xtime', xtime, timelevel=1)
    xtime = startTimeStamp

    ! Initialize initial_time in second time level. We need to do this because initial state
    ! is read into time level 1, and if we write output from the set of state arrays that
    ! represent the original time level 2, the initial_time field will be invalid.
    call mpas_pool_get_array(state, 'initial_time', initial_time1, timelevel=1)
    call mpas_pool_get_array(state, 'initial_time', initial_time2, timelevel=2)
    initial_time2 = initial_time1

    !
    ! Set time units to CF-compliant "seconds since <date and time>".
    !
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'state', state)
    call mpas_pool_get_field(state, 'Time', field_0d_real, timelevel=1)
    
    if (.not. associated(field_0d_real)) then
       call mpp_error(FATAL,subname//'Failed to find variable "Time"')
    end if

    call mpas_modify_att(field_0d_real % attlists(1) % attlist, 'units', &
         'seconds since ' // mpas_string_replace(initial_time1, '_', ' '), ierr=ierr)

    if (ierr /= 0) then
       call mpp_error(FATAL,subname//'Failed to set time units')
    end if

    call exchange_halo_group(domain_ptr, 'initialization:pv_edge,ru,rw',ierr=ierr)
    if ( ierr /= 0 ) then
       call mpp_error(FATAL,subname//'Failed to exchange halo layers for group "initialization:ru,rw"')
    end if

    !
    ! Prepare the dynamics for integration
    !
    call mpas_log_write('Initializing the dynamics')
    call mpas_atm_dynamics_init(domain_ptr)

    !
    ! Some additional "scratch" fields are needed for interoperability with CAM-SIMA, but they are not initialized by
    ! `mpas_atm_dynamics_init`. Initialize them below.
    !
!    call mpas_pool_get_field(domain_ptr % blocklist % allfields, 'tend_uzonal', field_2d_real, timelevel=1)
!    call mpas_allocate_scratch_field(field_2d_real)
!    nullify(field_2d_real)
    
!    call mpas_pool_get_field(domain_ptr % blocklist % allfields, 'tend_umerid', field_2d_real, timelevel=1)
!    call mpas_allocate_scratch_field(field_2d_real)
!    nullify(field_2d_real)

    call mpas_log_write('Successful initialization of MPAS dynamical core')

  end subroutine ufs_mpas_init_phase2

  !> #########################################################################################
  !> Routine to call MPAS dynamical core
  !> Loop over dynamical time-step(s) and increment MPAS state (timelevel 1->2)
  !>
  !> #########################################################################################
  subroutine ufs_mpas_run()
    ! MPAS
    use atm_core,             only : atm_do_timestep, atm_compute_output_diagnostics
    use mpas_domain_routines, only : mpas_pool_get_dimension
    use mpas_derived_types,   only : mpas_Time_type, mpas_pool_type, MPAS_TimeInterval_type
    use mpas_kind_types,      only : StrKIND, RKIND, R8KIND
    use mpas_constants,       only : rvord
    use mpas_pool_routines,   only : mpas_pool_get_config, mpas_pool_get_subpool
    use mpas_pool_routines,   only : mpas_pool_shift_time_levels, mpas_pool_get_array
    use mpas_log,             only : mpas_log_write
    use mpas_timer,           only : mpas_timer_start, mpas_timer_stop
    use mpas_timekeeping,     only : mpas_advance_clock, mpas_get_clock_time, mpas_get_time
    use mpas_timekeeping,     only : mpas_NOW, mpas_is_clock_stop_time, mpas_dmpar_get_time
    use mpas_timekeeping,     only : mpas_set_timeInterval, operator(+), operator(<)
    ! FMS
    use mpp_mod,              only : FATAL, mpp_error
    ! Locals
    character(len=*), parameter :: subname = 'ufs_mpas_run::ufs_mpas_run'
    real (kind=RKIND), pointer :: config_dt
    type (mpas_pool_type), pointer :: state, diag, mesh
    type (mpas_Time_type) :: timeNow, timeStop
    character(len=StrKIND) :: timeStamp
    integer :: ierr, itime, itimestep
    integer, pointer :: index_qv
    integer, pointer :: nCellsSolve
    real(kind=RKIND), dimension(:,:), pointer :: theta_m, rho_zz, zz, theta, rho
    real(kind=RKIND), dimension(:,:,:), pointer :: scalars
    real (kind=R8KIND) :: integ_start_time, integ_stop_time 
    logical, pointer :: config_apply_lbcs
    type(mpas_timeinterval_type) :: mpas_time_interval

    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'state', state)
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'diag',  diag)
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh',  mesh)
    
    ! Eventually, dt should be domain specific
    call mpas_pool_get_config( domain_ptr % blocklist % configs, 'config_dt', config_dt)
    call MPAS_set_timeInterval(mpas_time_interval, S=dt_atmos, ierr=ierr)
    if (ierr /= 0) then
       call mpp_error(FATAL,subname//'Failed to set dynamics time step')
    endif

    !
    ! Read initial boundary state
    ! NOT YET IMPLEMENTED (Follow src/core_atmosphere/mpas_atm_core.F:atm_core_run())
    !
    call mpas_pool_get_config( domain_ptr % blocklist % configs, 'config_apply_lbcs', config_apply_lbcs)
    if (config_apply_lbcs) then
       
    endif
    
    ! During integration, time level 1 stores the model state at the beginning of the
    !   time step, and time level 2 stores the state advanced config_dt in time by timestep(...)
    timeNow  = mpas_get_clock_time(clock, mpas_NOW, ierr)
    if (ierr /= 0) then
        call mpp_error(FATAL,subname//': Failed to get clock_time for "mpas_NOW"')
    endif

    timeStop = timeNow + mpas_time_interval
    itimestep =	0
    do while (timeNow < timeStop)
       itimestep = itimestep + 1
       !
       call mpas_get_time(curr_time=timeNow, dateTimeString=timeStamp, ierr=ierr)
       if ( ierr /= 0 ) then
          call mpp_error(FATAL,subname//': Failed to get time mpas_NOW"')
       end if
       call mpas_log_write('') 
       call mpas_log_write(' MPAS dynamics start timestep '//trim(timeStamp))

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
       call mpas_advance_clock(clock)
       timeNow = mpas_get_clock_time(clock, mpas_NOW, ierr)
       if (ierr /= 0) then
          call mpp_error(FATAL,subname//': Failed to get clock_time for "mpas_NOW"')
       endif

    end do

    !
    ! Compute diagnostic fields from the final prognostic state
    !
    call atm_compute_output_diagnostics(state, 1, diag, mesh)

  end subroutine ufs_mpas_run

  
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
    use mpas_log,           only : mpas_log_write
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
    real(RKIND)            :: config_relax_zone_divdamp_coef = 6

    ! Locals
    integer :: ierr, io, mpierr

    ! Read in namelists...
    if (me == master) then
       !print*,'Reading MPAS-A dynamical core namelist'
       call mpas_log_write('Reading MPAS-A dynamical core namelist')
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
    call mpas_pool_add_config(configPool, 'config_relax_zone_divdamp_coef',        config_relax_zone_divdamp_coef)

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
 ! subroutine dyn_mpas_read_write_stream
 !
 !> summary: Read or write an MPAS stream.
 !> author: Kuan-Chih Wang
 !> date: 2024-03-15
 !>
 !> In the context of MPAS, the concept of a "pool" resembles a group of
 !> (related) variables, while the concept of a "stream" resembles a file.
 !> This subroutine reads or writes an MPAS stream. It provides the mechanism
 !> for CAM-SIMA to input/output data to/from MPAS dynamical core.
 !> Analogous to the `{read,write}_stream` subroutines in MPAS stream manager.
 !
 !> ########################################################################################
 subroutine dyn_mpas_read_write_stream(stream_mode, stream_name)
   ! Module(s) from external libraries.
   use pio, only: file_desc_t
   use mpp_mod,             only : FATAL, mpp_error
   ! Module(s) from MPAS.
   use mpas_derived_types,  only : mpas_pool_type, mpas_stream_noerr, mpas_stream_type
   use mpas_io_streams,     only : mpas_closestream, mpas_readstream, mpas_writestream
   use mpas_pool_routines,  only : mpas_pool_destroy_pool
   use mpas_stream_manager, only : postread_reindex, prewrite_reindex, postwrite_reindex
   use mpas_log,            only : mpas_log_write
   use mpas_atm_halos,      only : exchange_halo_group

   character(*), intent(in) :: stream_mode
   character(*), intent(in) :: stream_name

   character(*), parameter :: subname = 'dyn_mpas_subdriver::dyn_mpas_read_write_stream'
   integer :: i, ierr
   type(mpas_pool_type), pointer :: mpas_pool
   type(mpas_stream_type), pointer :: mpas_stream
   type(var_info_type), allocatable :: var_info_list(:)

   call mpas_log_write('')

   nullify(mpas_pool)
   nullify(mpas_stream)

   call	mpas_log_write( 'Initializing stream "' // trim(adjustl(stream_name)) // '"')

   call dyn_mpas_init_stream_with_pool(mpas_pool, mpas_stream, pioid, stream_mode, stream_name)

   if (.not. associated(mpas_pool)) then
      call mpp_error(FATAL,subname//'Failed to initialize stream "' // trim(adjustl(stream_name)) // '"')
   end if

   if (.not. associated(mpas_stream)) then
      call mpp_error(FATAL,subname//'Failed to initialize stream "' // trim(adjustl(stream_name)) // '"')
   end if

   select case (trim(adjustl(stream_mode)))
   case ('r', 'read')
      call mpas_log_write('Reading stream "' // trim(adjustl(stream_name)) // '"')

      call mpas_readstream(mpas_stream, 1, ierr=ierr)

      if (ierr /= mpas_stream_noerr) then
         call mpp_error(FATAL,subname//'Failed to read stream "' // trim(adjustl(stream_name)) // '"')
      end if

      ! Exchange halo layers because new data have just been read.
      var_info_list = parse_stream_name(stream_name)

      do i = 1, size(var_info_list)
         call dyn_mpas_exchange_halo(var_info_list(i) % name)
         if ( ierr /= 0 ) then
            call mpp_error(FATAL,subname//'Failed to exchange halo layers for group '//var_info_list(i) % name)
         end if
      end do

      ! For any connectivity arrays in this stream, convert global indexes to local indexes.
      call postread_reindex(domain_ptr % blocklist % allfields, domain_ptr % packages, &
           mpas_pool, mpas_pool)
   case ('w', 'write')
      call mpas_log_write('Writing stream "' // trim(adjustl(stream_name)) // '"')

      ! WARNING:
      ! The `{pre,post}write_reindex` subroutines are STATEFUL because they store information inside their module
      ! (i.e., module variables). They MUST be called in pairs, like below, to prevent undefined behaviors.

      ! For any connectivity arrays in this stream, temporarily convert local indexes to global indexes.
      call prewrite_reindex(domain_ptr % blocklist % allfields, domain_ptr % packages, &
           mpas_pool, mpas_pool)

      call mpas_writestream(mpas_stream, 1, ierr=ierr)

      if (ierr /= mpas_stream_noerr) then
         call mpp_error(FATAL,subname//'Failed to write stream "' // trim(adjustl(stream_name)) // '"')
      end if

      ! For any connectivity arrays in this stream, reset global indexes back to local indexes.
      call postwrite_reindex(domain_ptr % blocklist % allfields, mpas_pool)
   case default
      call mpp_error(FATAL,subname//'Unsupported stream mode "' // trim(adjustl(stream_mode)) // '"')
   end select

   call mpas_log_write('Closing stream "' // trim(adjustl(stream_name)) // '"')

   call mpas_closestream(mpas_stream, ierr=ierr)

   if (ierr /= mpas_stream_noerr) then
      call mpp_error(FATAL,subname//'Failed to close stream "' // trim(adjustl(stream_name)) // '"')
   end if

   ! Deallocate temporary pointers to avoid memory leaks.
   call mpas_pool_destroy_pool(mpas_pool)
   nullify(mpas_pool)
   
   deallocate(mpas_stream)
   nullify(mpas_stream)

   call mpas_log_write(subname // ' completed')
 end subroutine dyn_mpas_read_write_stream

 !> ########################################################################################
 ! subroutine dyn_mpas_exchange_halo
 !
 !> summary: Update the halo layers of the named field.
 !> author: Michael Duda
 !> date: 16 January 2020
 !>
 !> Given a field name that is defined in MPAS registry, this subroutine updates
 !> the halo layers for that field.
 !> Ported and refactored for CAM-SIMA. (KCW, 2024-03-18)
 !> Ported and refactored for UWM (DJS: 2025)
 !
 !> ########################################################################################
 subroutine dyn_mpas_exchange_halo(field_name)
   ! Module(s) from MPAS.
   use mpas_derived_types, only : field1dinteger, field2dinteger, field3dinteger,           &
                                  field1dreal, field2dreal, field3dreal, field4dreal,       &
                                  field5dreal, mpas_pool_field_info_type, mpas_pool_integer,&
                                  mpas_pool_real
   use mpas_dmpar,         only : mpas_dmpar_exch_halo_field
   use mpas_pool_routines, only : mpas_pool_get_field, mpas_pool_get_field_info
   use mpp_mod,            only : FATAL, mpp_error
   use mpas_log,           only : mpas_log_write
   character(*), intent(in) :: field_name

   character(*), parameter :: subname = 'dyn_mpas_subdriver::dyn_mpas_exchange_halo'
   type(field1dinteger), pointer :: field_1d_integer
   type(field2dinteger), pointer :: field_2d_integer
   type(field3dinteger), pointer :: field_3d_integer
   type(field1dreal), pointer :: field_1d_real
   type(field2dreal), pointer :: field_2d_real
   type(field3dreal), pointer :: field_3d_real
   type(field4dreal), pointer :: field_4d_real
   type(field5dreal), pointer :: field_5d_real
   type(mpas_pool_field_info_type) :: mpas_pool_field_info

   call mpas_log_write(subname // ' entered')

   nullify(field_1d_integer)
   nullify(field_2d_integer)
   nullify(field_3d_integer)
   nullify(field_1d_real)
   nullify(field_2d_real)
   nullify(field_3d_real)
   nullify(field_4d_real)
   nullify(field_5d_real)

   call mpas_log_write('Inquiring field information for "' // trim(adjustl(field_name)) // '"')

   call mpas_pool_get_field_info(domain_ptr % blocklist % allfields, &
        trim(adjustl(field_name)), mpas_pool_field_info)

   if (mpas_pool_field_info % fieldtype == -1 .or. &
        mpas_pool_field_info % ndims == -1 .or. &
        mpas_pool_field_info % nhalolayers == -1) then
      call mpp_error(FATAL,subname//'Invalid field information for "' // trim(adjustl(field_name)) // '"')
   end if
   
   ! No halo layers to exchange. This field is not decomposed.
   if (mpas_pool_field_info % nhalolayers == 0) then
      call mpas_log_write('Skipping field "' // trim(adjustl(field_name)) // '" due to not decomposed')
      
      return
   end if
   
   call mpas_log_write('Exchanging halo layers for "' // trim(adjustl(field_name)) // '"')
   
   select case (mpas_pool_field_info % fieldtype)
   case (mpas_pool_integer)
      select case (mpas_pool_field_info % ndims)
      case (1)
         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(field_name)), field_1d_integer, timelevel=1)
         
         if (.not. associated(field_1d_integer)) then
            call mpp_error(FATAL,subname//'Failed to find field "' // trim(adjustl(field_name)) // '"')
         end if
         
         call mpas_dmpar_exch_halo_field(field_1d_integer)
         
         nullify(field_1d_integer)
      case (2)
         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(field_name)), field_2d_integer, timelevel=1)
                        
         if (.not. associated(field_2d_integer)) then
            call mpp_error(FATAL,subname//'Failed to find field "' // trim(adjustl(field_name)) // '"')
         end if

         call mpas_dmpar_exch_halo_field(field_2d_integer)
         
         nullify(field_2d_integer)
      case (3)
         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(field_name)), field_3d_integer, timelevel=1)

         if (.not. associated(field_3d_integer)) then
            call mpp_error(FATAL,subname//'Failed to find field "' // trim(adjustl(field_name)) // '"')
         end if

         call mpas_dmpar_exch_halo_field(field_3d_integer)

         nullify(field_3d_integer)
      case default
         call mpp_error(FATAL,subname//'Unsupported field rank ' // stringify([mpas_pool_field_info % ndims]))
      end select
   case (mpas_pool_real)
      select case (mpas_pool_field_info % ndims)
      case (1)
         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(field_name)), field_1d_real, timelevel=1)

         if (.not. associated(field_1d_real)) then
            call mpp_error(FATAL,subname//'Failed to find field "' // trim(adjustl(field_name)) // '"')
         end if

         call mpas_dmpar_exch_halo_field(field_1d_real)

         nullify(field_1d_real)
      case (2)
         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(field_name)), field_2d_real, timelevel=1)

         if (.not. associated(field_2d_real)) then
            call mpp_error(FATAL,subname//'Failed to find field "' // trim(adjustl(field_name)) // '"')
         end if

         call mpas_dmpar_exch_halo_field(field_2d_real)

         nullify(field_2d_real)
      case (3)
         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(field_name)), field_3d_real, timelevel=1)

         if (.not. associated(field_3d_real)) then
            call mpp_error(FATAL,subname//'Failed to find field "' // trim(adjustl(field_name)) // '"')
         end if

         call mpas_dmpar_exch_halo_field(field_3d_real)

         nullify(field_3d_real)
      case (4)
         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(field_name)), field_4d_real, timelevel=1)

         if (.not. associated(field_4d_real)) then
            call mpp_error(FATAL,subname//'Failed to find field "' // trim(adjustl(field_name)) // '"')
         end if

         call mpas_dmpar_exch_halo_field(field_4d_real)

         nullify(field_4d_real)
      case (5)
         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(field_name)), field_5d_real, timelevel=1)

         if (.not. associated(field_5d_real)) then
            call mpp_error(FATAL,subname//'Failed to find field "' // trim(adjustl(field_name)) // '"')
         end if

         call mpas_dmpar_exch_halo_field(field_5d_real)

         nullify(field_5d_real)
      case default
         call mpp_error(FATAL,subname//'Unsupported field rank ' // stringify([mpas_pool_field_info % ndims]))
      end select
   case default
      call mpp_error(FATAL,subname//'Unsupported field type (Must be one of: integer, real)')
   end select

   call mpas_log_write(subname // ' completed')
 end subroutine dyn_mpas_exchange_halo

 !> ######################################################################################## 
 ! subroutine dyn_mpas_init_stream_with_pool
 !
 !> summary: Initialize an MPAS stream with an accompanying MPAS pool.
 !> author: Kuan-Chih Wang
 !> date: 2024-03-14
 !>
 !> In the context of MPAS, the concept of a "pool" resembles a group of
 !> (related) variables, while the concept of a "stream" resembles a file.
 !> This subroutine initializes an MPAS stream with an accompanying MPAS pool by
 !> adding variable and attribute information to them. After that, MPAS is ready
 !> to perform IO on them.
 !> Analogous to the `build_stream` and `mpas_stream_mgr_add_field`
 !> subroutines in MPAS stream manager.
 !>
 !> Ported and refactored for UWM (DJS: 2025)
 !
 !> ######################################################################################## 
 subroutine dyn_mpas_init_stream_with_pool(mpas_pool, mpas_stream, pio_file, stream_mode,  &
                                           stream_name)
   ! Module(s) from external libraries.
   use pio, only: file_desc_t, pio_file_is_open
   ! Module(s) from MPAS.
   use mpas_derived_types, only : field0dchar, field1dchar, field0dinteger, field1dinteger,&
                                  field2dinteger, field3dinteger, field0dreal, field1dreal,&
                                  field2dreal, field3dreal, field4dreal, field5dreal,      &
                                  mpas_io_native_precision, mpas_io_pnetcdf, mpas_io_read, &
                                  mpas_io_write, mpas_pool_type, mpas_stream_noerr,        &
                                  mpas_stream_type
   use mpas_io_streams,    only : mpas_createstream, mpas_streamaddfield
   use mpas_pool_routines, only : mpas_pool_add_config, mpas_pool_create_pool, mpas_pool_get_field
   use mpas_kind_types,    only : StrKIND, RKIND
   use mpp_mod,            only : FATAL, mpp_error
   use mpas_log,           only : mpas_log_write

   type(mpas_pool_type), pointer, intent(out) :: mpas_pool
   type(mpas_stream_type), pointer, intent(out) :: mpas_stream
   type(file_desc_t), pointer, intent(in) :: pio_file
   character(*), intent(in) :: stream_mode
   character(*), intent(in) :: stream_name

   interface add_stream_attribute
      procedure :: add_stream_attribute_0d
      procedure :: add_stream_attribute_1d
   end interface add_stream_attribute

   character(*), parameter :: subname = 'dyn_mpas_subdriver::dyn_mpas_init_stream_with_pool'
   character(strkind) :: stream_filename
   integer :: i, ierr, stream_format
   !> Whether a variable is present on the file (i.e., `pio_file`).
   logical, allocatable :: var_is_present(:)
   !> Whether a variable is type, kind, and rank compatible with what MPAS expects on the file (i.e., `pio_file`).
   logical, allocatable :: var_is_tkr_compatible(:)
   type(field0dchar), pointer :: field_0d_char
   type(field1dchar), pointer :: field_1d_char
   type(field0dinteger), pointer :: field_0d_integer
   type(field1dinteger), pointer :: field_1d_integer
   type(field2dinteger), pointer :: field_2d_integer
   type(field3dinteger), pointer :: field_3d_integer
   type(field0dreal), pointer :: field_0d_real
   type(field1dreal), pointer :: field_1d_real
   type(field2dreal), pointer :: field_2d_real
   type(field3dreal), pointer :: field_3d_real
   type(field4dreal), pointer :: field_4d_real
   type(field5dreal), pointer :: field_5d_real
   type(var_info_type), allocatable :: var_info_list(:)
   
   call mpas_log_write(subname // ' entered')
   
   nullify(field_0d_char)
   nullify(field_1d_char)
   nullify(field_0d_integer)
   nullify(field_1d_integer)
   nullify(field_2d_integer)
   nullify(field_3d_integer)
   nullify(field_0d_real)
   nullify(field_1d_real)
   nullify(field_2d_real)
   nullify(field_3d_real)
   nullify(field_4d_real)
   nullify(field_5d_real)

   call mpas_pool_create_pool(mpas_pool)

   allocate(mpas_stream, stat=ierr)

   if (ierr /= 0) then
      call mpp_error(FATAL,subname//'Failed to allocate stream "' // trim(adjustl(stream_name)) // '"')
   end if

   ! Not actually used because a PIO file descriptor is directly supplied.
   stream_filename = 'external stream'
   stream_format = mpas_io_pnetcdf

   call mpas_log_write('Checking PIO file descriptor')

   if (.not. associated(pio_file)) then
      call mpp_error(FATAL,subname//'Invalid PIO file descriptor')
   end if

   if (.not. pio_file_is_open(pio_file)) then
      call mpp_error(FATAL,subname//'Invalid PIO file descriptor')
   end if

   select case (trim(adjustl(stream_mode)))
   case ('r', 'read')
      call mpas_log_write('Creating stream "' // trim(adjustl(stream_name)) // '" for reading')

      call mpas_createstream( &
           mpas_stream, domain_ptr % iocontext, stream_filename, stream_format, mpas_io_read,  &
           clobberrecords=.false., clobberfiles=.false., truncatefiles=.false., &
           precision=mpas_io_native_precision, pio_file_desc=pio_file, ierr=ierr)
   case ('w', 'write')
      call mpas_log_write('Creating stream "' // trim(adjustl(stream_name)) // '" for writing')

      call mpas_createstream( &
           mpas_stream, domain_ptr % iocontext, stream_filename, stream_format, mpas_io_write, &
           clobberrecords=.false., clobberfiles=.false., truncatefiles=.false., &
           precision=mpas_io_native_precision, pio_file_desc=pio_file, ierr=ierr)
   case default
      call mpp_error(FATAL,subname//'Unsupported stream mode "' // trim(adjustl(stream_mode)) // '"')
   end select

   if (ierr /= mpas_stream_noerr) then
      call mpp_error(FATAL,subname//'Failed to create stream "' // trim(adjustl(stream_name)) // '"')
   end if

   var_info_list = parse_stream_name(stream_name)

   ! Add variables contained in `var_info_list` to stream.
   do i = 1, size(var_info_list)
      call mpas_log_write('var_info_list(' // stringify([i]) // ') % name = ' // stringify([var_info_list(i) % name]))
      call mpas_log_write('var_info_list(' // stringify([i]) // ') % type = ' // stringify([var_info_list(i) % type]))
      call mpas_log_write('var_info_list(' // stringify([i]) // ') % rank = ' // stringify([var_info_list(i) % rank]))

      if (trim(adjustl(stream_mode)) == 'r' .or. trim(adjustl(stream_mode)) == 'read') then
         call dyn_mpas_check_variable_status(var_is_present, var_is_tkr_compatible, pio_file, var_info_list(i))

         ! Do not hard crash the model if a variable is missing and cannot be read.
         ! This can happen if users attempt to initialize/restart the model with data generated by
         ! older versions of MPAS. Print a debug message to let users decide if this is acceptable.
         if (.not. any(var_is_present)) then
            call mpas_log_write('Skipping variable "' // trim(adjustl(var_info_list(i) % name)) // '" due to not present')

            cycle
         end if

         if (any(var_is_present .and. .not. var_is_tkr_compatible)) then
            call mpas_log_write('Skipping variable "' // trim(adjustl(var_info_list(i) % name)) // '" due to not TKR compatible')

            !cycle
         end if
      end if

      ! Add "<variable name>" to pool with the value of `1`.
      ! The existence of "<variable name>" in pool causes it to be considered for IO in MPAS.
      call mpas_pool_add_config(mpas_pool, trim(adjustl(var_info_list(i) % name)), 1)
      ! Add "<variable name>:packages" to pool with the value of an empty character string.
      ! This causes "<variable name>" to be always considered active for IO in MPAS.
      !call mpas_pool_add_config(mpas_pool, trim(adjustl(var_info_list(i) % name) // ':packages'), '')

      ! Add "<variable name>" to stream.
      call mpas_log_write('Adding variable "' // trim(adjustl(var_info_list(i) % name)) // &
           '" to stream "' // trim(adjustl(stream_name)) // '"')

      select case (trim(adjustl(var_info_list(i) % type)))
      case ('character')
         select case (var_info_list(i) % rank)
         case (0)
            call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
                 trim(adjustl(var_info_list(i) % name)), field_0d_char, timelevel=1)

            if (.not. associated(field_0d_char)) then
               call mpp_error(FATAL,subname//'Failed to find variable "' // trim(adjustl(var_info_list(i) % name)) // '"')
            end if

            call mpas_streamaddfield(mpas_stream, field_0d_char, ierr=ierr)

            nullify(field_0d_char)
         case (1)
            call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
                 trim(adjustl(var_info_list(i) % name)), field_1d_char, timelevel=1)

            if (.not. associated(field_1d_char)) then
               call mpp_error(FATAL,subname//'Failed to find variable "' // trim(adjustl(var_info_list(i) % name)) // '"')
            end if

            call mpas_streamaddfield(mpas_stream, field_1d_char, ierr=ierr)

            nullify(field_1d_char)
         case default
            call mpp_error(FATAL,subname//'Unsupported variable rank ' // stringify([var_info_list(i) % rank]) // &
                 ' for "' // trim(adjustl(var_info_list(i) % name)) // '"')
         end select
      case ('integer')
         select case (var_info_list(i) % rank)
         case (0)
            call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
                 trim(adjustl(var_info_list(i) % name)), field_0d_integer, timelevel=1)

            if (.not. associated(field_0d_integer)) then
               call mpp_error(FATAL,subname//'Failed to find variable "' // trim(adjustl(var_info_list(i) % name)) // '"')
            end if

            call mpas_streamaddfield(mpas_stream, field_0d_integer, ierr=ierr)

            nullify(field_0d_integer)
         case (1)
            call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
                 trim(adjustl(var_info_list(i) % name)), field_1d_integer, timelevel=1)

            if (.not. associated(field_1d_integer)) then
               call mpp_error(FATAL,subname//'Failed to find variable "' // trim(adjustl(var_info_list(i) % name)) // '"')
            end if

            call mpas_streamaddfield(mpas_stream, field_1d_integer, ierr=ierr)

            nullify(field_1d_integer)
         case (2)
            call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
                 trim(adjustl(var_info_list(i) % name)), field_2d_integer, timelevel=1)
            
            if (.not. associated(field_2d_integer)) then
               call mpp_error(FATAL,subname//'Failed to find variable "' // trim(adjustl(var_info_list(i) % name)) // '"')
            end if

            call mpas_streamaddfield(mpas_stream, field_2d_integer, ierr=ierr)

            nullify(field_2d_integer)
         case (3)
            call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
                 trim(adjustl(var_info_list(i) % name)), field_3d_integer, timelevel=1)

            if (.not. associated(field_3d_integer)) then
               call mpp_error(FATAL,subname//'Failed to find variable "' // trim(adjustl(var_info_list(i) % name)) // '"')
            end if

            call mpas_streamaddfield(mpas_stream, field_3d_integer, ierr=ierr)

            nullify(field_3d_integer)
         case default
            call mpp_error(FATAL,subname//'Unsupported variable rank ' // stringify([var_info_list(i) % rank]) // &
                 ' for "' // trim(adjustl(var_info_list(i) % name)) // '"')
         end select
      case ('real')
         select case (var_info_list(i) % rank)
         case (0)
            call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
                 trim(adjustl(var_info_list(i) % name)), field_0d_real, timelevel=1)

            if (.not. associated(field_0d_real)) then
               call mpp_error(FATAL,subname//'Failed to find variable "' // trim(adjustl(var_info_list(i) % name)) // '"')
            end if

            call mpas_streamaddfield(mpas_stream, field_0d_real, ierr=ierr)
            
            nullify(field_0d_real)
         case (1)
            call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
                 trim(adjustl(var_info_list(i) % name)), field_1d_real, timelevel=1)

            if (.not. associated(field_1d_real)) then
               call mpp_error(FATAL,subname//'Failed to find variable "' // trim(adjustl(var_info_list(i) % name)) // '"')
            end if

            call mpas_streamaddfield(mpas_stream, field_1d_real, ierr=ierr)

            nullify(field_1d_real)
         case (2)
            call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
                 trim(adjustl(var_info_list(i) % name)), field_2d_real, timelevel=1)

            if (.not. associated(field_2d_real)) then
               call mpp_error(FATAL,subname//'Failed to find variable "' // trim(adjustl(var_info_list(i) % name)) // '"')
            end if
            
            call mpas_streamaddfield(mpas_stream, field_2d_real, ierr=ierr)

            nullify(field_2d_real)
         case (3)
            call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
                 trim(adjustl(var_info_list(i) % name)), field_3d_real, timelevel=1)

            if (.not. associated(field_3d_real)) then
               call mpp_error(FATAL,subname//'Failed to find variable "' // trim(adjustl(var_info_list(i) % name)) // '"')
            end if

            call mpas_streamaddfield(mpas_stream, field_3d_real, ierr=ierr)

            nullify(field_3d_real)
         case (4)
            call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
                 trim(adjustl(var_info_list(i) % name)), field_4d_real, timelevel=1)

            if (.not. associated(field_4d_real)) then
               call mpp_error(FATAL,subname//'Failed to find variable "' // trim(adjustl(var_info_list(i) % name)) // '"')
            end if

            call mpas_streamaddfield(mpas_stream, field_4d_real, ierr=ierr)

            nullify(field_4d_real)
         case (5)
            call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
                 trim(adjustl(var_info_list(i) % name)), field_5d_real, timelevel=1)

            if (.not. associated(field_5d_real)) then
               call mpp_error(FATAL,subname//'Failed to find variable "' // trim(adjustl(var_info_list(i) % name)) // '"')
            end if

            call mpas_streamaddfield(mpas_stream, field_5d_real, ierr=ierr)
            
            nullify(field_5d_real)
         case default
            call mpp_error(FATAL,subname//'Unsupported variable rank ' // stringify([var_info_list(i) % rank]) // &
                 ' for "' // trim(adjustl(var_info_list(i) % name)) // '"')
         end select
      case default
         call mpp_error(FATAL,subname//'Unsupported variable type "' // trim(adjustl(var_info_list(i) % type)) // &
              '" for "' // trim(adjustl(var_info_list(i) % name)) // '"')
      end select

      if (ierr /= mpas_stream_noerr) then
         call mpp_error(FATAL,subname//'Failed to add variable "' // trim(adjustl(var_info_list(i) % name)) // &
              '" to stream "' // trim(adjustl(stream_name)) // '"')
      end if
   end do

   if (trim(adjustl(stream_mode)) == 'w' .or. trim(adjustl(stream_mode)) == 'write') then
      ! Add MPAS-specific attributes to stream.

      ! Attributes related to MPAS core (i.e., `core_type`).
      call add_stream_attribute('conventions', domain_ptr % core % conventions)
      call add_stream_attribute('core_name', domain_ptr % core % corename)
      call add_stream_attribute('git_version', domain_ptr % core % git_version)
      call add_stream_attribute('model_name', domain_ptr % core % modelname)
      call add_stream_attribute('source', domain_ptr % core % source)

      ! Attributes related to MPAS domain (i.e., `domain_type`).
      call add_stream_attribute('is_periodic', domain_ptr % is_periodic)
      call add_stream_attribute('mesh_spec', domain_ptr % mesh_spec)
      call add_stream_attribute('on_a_sphere', domain_ptr % on_a_sphere)
      call add_stream_attribute('parent_id',  domain_ptr % parent_id)
      call add_stream_attribute('sphere_radius', domain_ptr % sphere_radius)
      call add_stream_attribute('x_period',  domain_ptr % x_period)
      call add_stream_attribute('y_period',  domain_ptr % y_period)
   end if

   call mpas_log_write(subname // ' completed')
 contains
   !> Helper subroutine for adding a 0-d stream attribute by calling `mpas_writestreamatt` with error checking.
   !> (KCW, 2024-03-14)
   subroutine add_stream_attribute_0d(attribute_name, attribute_value)
     ! Module(s) from MPAS.
     use mpas_io_streams, only : mpas_writestreamatt
     use mpas_log,        only : mpas_log_write
     character(*), intent(in) :: attribute_name
     class(*), intent(in) :: attribute_value

     call mpas_log_write('Adding attribute "' // trim(adjustl(attribute_name)) // &
          '" to stream "' // trim(adjustl(stream_name)) // '"')

     select type (attribute_value)
     type is (character(*))
        call mpas_writestreamatt(mpas_stream, &
             trim(adjustl(attribute_name)), trim(adjustl(attribute_value)), syncval=.false., ierr=ierr)
     type is (integer)
        call mpas_writestreamatt(mpas_stream, &
             trim(adjustl(attribute_name)), attribute_value, syncval=.false., ierr=ierr)
     type is (logical)
        if (attribute_value) then
           ! Logical `.true.` becomes character string "YES".
           call mpas_writestreamatt(mpas_stream, &
                trim(adjustl(attribute_name)), 'YES', syncval=.false., ierr=ierr)
        else
           ! Logical `.false.` becomes character string "NO".
           call mpas_writestreamatt(mpas_stream, &
                trim(adjustl(attribute_name)), 'NO', syncval=.false., ierr=ierr)
        end if
     type is (real(rkind))
        call mpas_writestreamatt(mpas_stream, &
             trim(adjustl(attribute_name)), attribute_value, syncval=.false., ierr=ierr)
     class default
        call mpp_error(FATAL,subname//'Unsupported attribute type (Must be one of: character, integer, logical, real)')
     end select

     if (ierr /= mpas_stream_noerr) then
        call mpp_error(FATAL,subname//'Failed to add attribute "' // trim(adjustl(attribute_name)) // &
             '" to stream "' // trim(adjustl(stream_name)) // '"')
     end if
   end subroutine add_stream_attribute_0d

   !> Helper subroutine for adding a 1-d stream attribute by calling `mpas_writestreamatt` with error checking.
   !> (KCW, 2024-03-14)
   subroutine add_stream_attribute_1d(attribute_name, attribute_value)
     ! Module(s) from MPAS.
     use mpas_io_streams, only : mpas_writestreamatt
     use mpas_log,        only : mpas_log_write
     character(*), intent(in) :: attribute_name
     class(*), intent(in) :: attribute_value(:)

     call mpas_log_write('Adding attribute "' // trim(adjustl(attribute_name)) // &
          '" to stream "' // trim(adjustl(stream_name)) // '"')
     
     select type (attribute_value)
     type is (integer)
        call mpas_writestreamatt(mpas_stream, &
             trim(adjustl(attribute_name)), attribute_value, syncval=.false., ierr=ierr)
     type is (real(rkind))
        call mpas_writestreamatt(mpas_stream, &
             trim(adjustl(attribute_name)), attribute_value, syncval=.false., ierr=ierr)
     class default
        call mpp_error(FATAL,subname//'Unsupported attribute type (Must be one of: integer, real)')
     end select

     if (ierr /= mpas_stream_noerr) then
        call mpp_error(FATAL,subname//'Failed to add attribute "' // trim(adjustl(attribute_name)) // &
             '" to stream "' // trim(adjustl(stream_name)) // '"')
     end if
   end subroutine add_stream_attribute_1d
 end subroutine dyn_mpas_init_stream_with_pool
 
 !> Parse a stream name, which consists of one or more stream name fragments, and return the corresponding variable information
 !> as a list of `var_info_type`. Multiple stream name fragments should be separated by "+" (i.e., a plus, meaning "addition"
 !> operation) or "-" (i.e., a minus, meaning "subtraction" operation).
 !> A stream name fragment can be a predefined stream name (e.g., "invariant", "input", etc.) or a single variable name.
 !> For example, a stream name of "invariant+input+restart" means the union of variables in the "invariant", "input", and
 !> "restart" streams.
 !> Duplicate variable information in the resulting list is discarded.
 !> (KCW, 2024-06-01)
 pure function parse_stream_name(stream_name) result(var_info_list)
   character(*), intent(in) :: stream_name
   type(var_info_type), allocatable :: var_info_list(:)
        
   character(*), parameter :: supported_stream_name_operator = '+-'
   character(1) :: stream_name_operator
   character(:), allocatable :: stream_name_fragment
   character(len(invariant_var_info_list % name)), allocatable :: var_name_list(:)
   integer :: i, j, n, offset
   type(var_info_type), allocatable :: var_info_list_buffer(:)

   n = len_trim(stream_name)

   if (n == 0) then
      ! Empty character string means empty list.
      var_info_list = parse_stream_name_fragment('')

      return
   end if

   i = scan(stream_name, supported_stream_name_operator)

   if (i == 0) then
      ! No operators are present in the stream name. It is just a single stream name fragment.
      stream_name_fragment = stream_name
      var_info_list = parse_stream_name_fragment(stream_name_fragment)

      return
   end if

   offset = 0
   var_info_list = parse_stream_name_fragment('')

   do while (.true.)
      ! Extract operator from the stream name.
      if (offset > 0) then
         stream_name_operator = stream_name(offset:offset)
      else
         stream_name_operator = '+'
      end if

      ! Extract stream name fragment from the stream name.
      if (i > 1) then
         stream_name_fragment = stream_name(offset + 1:offset + i - 1)
      else
         stream_name_fragment = ''
      end if

      ! Process the stream name fragment according to the operator.
      if (len_trim(stream_name_fragment) > 0) then
         var_info_list_buffer = parse_stream_name_fragment(stream_name_fragment)
         
         select case (stream_name_operator)
         case ('+')
            var_info_list = [var_info_list, var_info_list_buffer]
         case ('-')
            do j = 1, size(var_info_list_buffer)
               var_name_list = var_info_list % name
               var_info_list = pack(var_info_list, var_name_list /= var_info_list_buffer(j) % name)
            end do
         case default
            ! Do nothing for unknown operators. Should not happen at all.
         end select
      end if

      offset = offset + i

      ! Terminate loop when everything in the stream name has been processed.
      if (offset + 1 > n) then
         exit
      end if

      i = scan(stream_name(offset + 1:), supported_stream_name_operator)
      
      ! Run the loop one last time for the remaining stream name fragment.
      if (i == 0) then
         i = n - offset + 1
      end if
   end do

   ! Discard duplicate variable information by names.
   var_name_list = var_info_list % name
   var_info_list = var_info_list(index_unique(var_name_list))
 end function parse_stream_name

 !> Parse a stream name fragment and return the corresponding variable information as a list of `var_info_type`.
 !> A stream name fragment can be a predefined stream name (e.g., "invariant", "input", etc.) or a single variable name.
 !> (KCW, 2024-06-01)
 pure function parse_stream_name_fragment(stream_name_fragment) result(var_info_list)
   character(*), intent(in) :: stream_name_fragment
   type(var_info_type), allocatable :: var_info_list(:)

   character(len(invariant_var_info_list % name)), allocatable :: var_name_list(:)
   type(var_info_type), allocatable :: var_info_list_buffer(:)

   select case (trim(adjustl(stream_name_fragment)))
   case ('')
      allocate(var_info_list(0))
   case ('invariant')
      allocate(var_info_list, source=invariant_var_info_list)
   case ('input')
      allocate(var_info_list, source=input_var_info_list)
   case ('restart')
      allocate(var_info_list, source=restart_var_info_list)
   case ('output')
      allocate(var_info_list, source=output_var_info_list)
   case default
      allocate(var_info_list(0))
      
      var_name_list = invariant_var_info_list % name
      
      if (any(var_name_list == trim(adjustl(stream_name_fragment)))) then
         var_info_list_buffer = pack(invariant_var_info_list, var_name_list == trim(adjustl(stream_name_fragment)))
         var_info_list = [var_info_list, var_info_list_buffer]
      end if
      
      var_name_list = input_var_info_list % name
      
      if (any(var_name_list == trim(adjustl(stream_name_fragment)))) then
         var_info_list_buffer = pack(input_var_info_list, var_name_list == trim(adjustl(stream_name_fragment)))
         var_info_list = [var_info_list, var_info_list_buffer]
      end if
      
      var_name_list = restart_var_info_list % name
      
      if (any(var_name_list == trim(adjustl(stream_name_fragment)))) then
         var_info_list_buffer = pack(restart_var_info_list, var_name_list == trim(adjustl(stream_name_fragment)))
         var_info_list = [var_info_list, var_info_list_buffer]
      end if
      
      var_name_list = output_var_info_list % name
      
      if (any(var_name_list == trim(adjustl(stream_name_fragment)))) then
         var_info_list_buffer = pack(output_var_info_list, var_name_list == trim(adjustl(stream_name_fragment)))
         var_info_list = [var_info_list, var_info_list_buffer]
      end if
   end select
 end function parse_stream_name_fragment

 !> Return the index of unique elements in `array`, which can be any intrinsic data types, as an integer array.
 !> If `array` contains zero element or is of unsupported data types, an empty integer array is produced.
 !> For example, `index_unique([1, 2, 3, 1, 2, 3, 4, 5])` returns `[1, 2, 3, 7, 8]`.
 !> (KCW, 2024-03-22)
 pure function index_unique(array)
   use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64

   class(*), intent(in) :: array(:)
   integer, allocatable :: index_unique(:)

   character(:), allocatable :: array_c(:)
   integer :: i, n
   logical :: mask_unique(size(array))

   n = size(array)
   
   if (n == 0) then
      allocate(index_unique(0))

      return
   end if

   mask_unique = .false.

   select type (array)
   type is (character(*))
      ! Workaround for a bug in GNU Fortran >= 12. This is perhaps the manifestation of GCC Bugzilla Bug 100819.
      ! When a character string array is passed as the actual argument to an unlimited polymorphic dummy argument,
      ! its array index and length parameter are mishandled.
      allocate(character(len(array)) :: array_c(size(array)))
         
      array_c(:) = array(:)
         
      do i = 1, n
         if (.not. any(array_c(i) == array_c .and. mask_unique)) then
            mask_unique(i) = .true.
         end if
      end do
         
      deallocate(array_c)
   type is (integer(int32))
      do i = 1, n
         if (.not. any(array(i) == array .and. mask_unique)) then
            mask_unique(i) = .true.
         end if
      end do
   type is (integer(int64))
      do i = 1, n
         if (.not. any(array(i) == array .and. mask_unique)) then
            mask_unique(i) = .true.
         end if
      end do
   type is (logical)
      do i = 1, n
         if (.not. any((array(i) .eqv. array) .and. mask_unique)) then
            mask_unique(i) = .true.
         end if
      end do
   type is (real(real32))
      do i = 1, n
         if (.not. any(array(i) == array .and. mask_unique)) then
            mask_unique(i) = .true.
         end if
      end do
   type is (real(real64))
      do i = 1, n
         if (.not. any(array(i) == array .and. mask_unique)) then
            mask_unique(i) = .true.
         end if
      end do
   class default
      allocate(index_unique(0))

      return
   end select
      
   index_unique = pack([(i, i = 1, n)], mask_unique)
 end function index_unique

 !> ######################################################################################## 
 ! subroutine dyn_mpas_check_variable_status
 !
 !> summary: Check and return variable status on the given file.
 !> author: Kuan-Chih Wang
 !> date: 2024-06-04
 !>
 !> On the given file (i.e., `pio_file`), this subroutine checks whether the
 !> given variable (i.e., `var_info`) is present, and whether it is "TKR"
 !> compatible with what MPAS expects. "TKR" means type, kind, and rank.
 !> This subroutine can handle both ordinary variables and variable arrays.
 !> They are indicated by the `var` and `var_array` elements, respectively,
 !> in MPAS registry. For an ordinary variable, the checks are performed on
 !> itself. Otherwise, for a variable array, the checks are performed on its
 !> constituent parts instead.
 !
 !> ######################################################################################## 
 subroutine dyn_mpas_check_variable_status(var_is_present, var_is_tkr_compatible, pio_file,&
                                           var_info)
   ! Module(s) from external libraries.
   use pio, only: file_desc_t, pio_file_is_open, pio_char, pio_int, pio_real, pio_double,  &
                  pio_inq_varid, pio_inq_varndims, pio_inq_vartype, pio_noerr
   ! Module(s) from MPAS.
   use mpas_derived_types, only : field0dchar, field1dchar, field0dinteger, field1dinteger,&
                                  field2dinteger, field3dinteger, field0dreal, field1dreal,&
                                  field2dreal, field3dreal, field4dreal, field5dreal
   use mpas_kind_types,    only : r4kind, r8kind
   use mpas_pool_routines, only : mpas_pool_get_field
   use mpas_log,           only : mpas_log_write
   use mpas_kind_types,    only : StrKIND, RKIND
   use mpp_mod,            only : FATAL, mpp_error
   
   logical, allocatable, intent(out) :: var_is_present(:)
   logical, allocatable, intent(out) :: var_is_tkr_compatible(:)
   type(file_desc_t), pointer, intent(in) :: pio_file
   type(var_info_type), intent(in) :: var_info

   character(*), parameter :: subname = 'dyn_mpas_subdriver::dyn_mpas_check_variable_status'
   character(strkind), allocatable :: var_name_list(:)
   integer :: i, ierr, varid, varndims, vartype
   type(field0dchar), pointer :: field_0d_char
   type(field1dchar), pointer :: field_1d_char
   type(field0dinteger), pointer :: field_0d_integer
   type(field1dinteger), pointer :: field_1d_integer
   type(field2dinteger), pointer :: field_2d_integer
   type(field3dinteger), pointer :: field_3d_integer
   type(field0dreal), pointer :: field_0d_real
   type(field1dreal), pointer :: field_1d_real
   type(field2dreal), pointer :: field_2d_real
   type(field3dreal), pointer :: field_3d_real
   type(field4dreal), pointer :: field_4d_real
   type(field5dreal), pointer :: field_5d_real

   call mpas_log_write(subname // ' entered')

   nullify(field_0d_char)
   nullify(field_1d_char)
   nullify(field_0d_integer)
   nullify(field_1d_integer)
   nullify(field_2d_integer)
   nullify(field_3d_integer)
   nullify(field_0d_real)
   nullify(field_1d_real)
   nullify(field_2d_real)
   nullify(field_3d_real)
   nullify(field_4d_real)
   nullify(field_5d_real)

   ! Extract a list of variable names to check on the file.
   ! For an ordinary variable, this list just contains its name.
   ! For a variable array, this list contains the names of its constituent parts.
   select case (trim(adjustl(var_info % type)))
   case ('character')
      select case (var_info % rank)
      case (0)
         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(var_info % name)), field_0d_char, timelevel=1)

         if (.not. associated(field_0d_char)) then
            call mpp_error(FATAL,subname//'Failed to find variable "' // trim(adjustl(var_info % name)))
         end if

         if (field_0d_char % isvararray .and. associated(field_0d_char % constituentnames)) then
            allocate(var_name_list(size(field_0d_char % constituentnames)), stat=ierr)

            if (ierr /= 0) then
               call mpp_error(FATAL,subname//'Failed to allocate var_name_list')
            end if
            
            var_name_list(:) = field_0d_char % constituentnames(:)
         end if

         nullify(field_0d_char)
      case (1)
         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(var_info % name)), field_1d_char, timelevel=1)

         if (.not. associated(field_1d_char)) then
            call mpp_error(FATAL,subname//'Failed to find variable "' // trim(adjustl(var_info % name)))
         end if

         if (field_1d_char % isvararray .and. associated(field_1d_char % constituentnames)) then
            allocate(var_name_list(size(field_1d_char % constituentnames)), stat=ierr)

            if (ierr /= 0) then
               call mpp_error(FATAL,subname//'Failed to allocate var_name_list')
            end if

            var_name_list(:) = field_1d_char % constituentnames(:)
         end if

         nullify(field_1d_char)
      case default
         call mpp_error(FATAL,subname//'Unsupported variable rank ' // stringify([var_info % rank]) // &
              ' for "' // trim(adjustl(var_info % name)) // '"')
      end select
   case ('integer')
      select case (var_info % rank)
      case (0)
         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(var_info % name)), field_0d_integer, timelevel=1)

         if (.not. associated(field_0d_integer)) then
            call mpp_error(FATAL,subname//'Failed to find variable "' // trim(adjustl(var_info % name)) // '"')
         end if

         if (field_0d_integer % isvararray .and. associated(field_0d_integer % constituentnames)) then
            allocate(var_name_list(size(field_0d_integer % constituentnames)), stat=ierr)

            if (ierr /= 0) then
               call mpp_error(FATAL,subname//'Failed to allocate var_name_list')
            end if

            var_name_list(:) = field_0d_integer % constituentnames(:)
         end if

         nullify(field_0d_integer)
      case (1)
         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(var_info % name)), field_1d_integer, timelevel=1)

         if (.not. associated(field_1d_integer)) then
            call mpp_error(FATAL,subname//'Failed to find variable "' // trim(adjustl(var_info % name)) // '"')
         end if

         if (field_1d_integer % isvararray .and. associated(field_1d_integer % constituentnames)) then
            allocate(var_name_list(size(field_1d_integer % constituentnames)), stat=ierr)

            if (ierr /= 0) then
               call mpp_error(FATAL,subname//'Failed to allocate var_name_list')
            end if

            var_name_list(:) = field_1d_integer % constituentnames(:)
         end if

         nullify(field_1d_integer)
      case (2)
         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(var_info % name)), field_2d_integer, timelevel=1)

         if (.not. associated(field_2d_integer)) then
            call mpp_error(FATAL,subname//'Failed to find variable "' // trim(adjustl(var_info % name)) // '"')
         end if

         if (field_2d_integer % isvararray .and. associated(field_2d_integer % constituentnames)) then
            allocate(var_name_list(size(field_2d_integer % constituentnames)), stat=ierr)

            if (ierr /= 0) then
               call mpp_error(FATAL,subname//'Failed to allocate var_name_list')
            end if

            var_name_list(:) = field_2d_integer % constituentnames(:)
         end if

         nullify(field_2d_integer)
      case (3)
         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(var_info % name)), field_3d_integer, timelevel=1)

         if (.not. associated(field_3d_integer)) then
            call mpp_error(FATAL,subname//'Failed to find variable "' // trim(adjustl(var_info % name)) // '"')
         end if

         if (field_3d_integer % isvararray .and. associated(field_3d_integer % constituentnames)) then
            allocate(var_name_list(size(field_3d_integer % constituentnames)), stat=ierr)

            if (ierr /= 0) then
               call mpp_error(FATAL,subname//'Failed to allocate var_name_list')
            end if

            var_name_list(:) = field_3d_integer % constituentnames(:)
         end if

         nullify(field_3d_integer)
      case default
         call mpp_error(FATAL,subname//'Unsupported variable rank ' // stringify([var_info % rank]) // &
              ' for "' // trim(adjustl(var_info % name)) // '"')
      end select
   case ('real')
      select case (var_info % rank)
      case (0)
         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(var_info % name)), field_0d_real, timelevel=1)

         if (.not. associated(field_0d_real)) then
            call mpp_error(FATAL,subname//'Failed to find variable "' // trim(adjustl(var_info % name)) // '"')
         end if

         if (field_0d_real % isvararray .and. associated(field_0d_real % constituentnames)) then
            allocate(var_name_list(size(field_0d_real % constituentnames)), stat=ierr)

            if (ierr /= 0) then
               call mpp_error(FATAL,subname//'Failed to allocate var_name_list')
            end if

            var_name_list(:) = field_0d_real % constituentnames(:)
         end if

         nullify(field_0d_real)
      case (1)
         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(var_info % name)), field_1d_real, timelevel=1)

         if (.not. associated(field_1d_real)) then
            call mpp_error(FATAL,subname//'Failed to find variable "' // trim(adjustl(var_info % name)) // '"')
         end if

         if (field_1d_real % isvararray .and. associated(field_1d_real % constituentnames)) then
            allocate(var_name_list(size(field_1d_real % constituentnames)), stat=ierr)

            if (ierr /= 0) then
               call mpp_error(FATAL,subname//'Failed to allocate var_name_list')
            end if

            var_name_list(:) = field_1d_real % constituentnames(:)
         end if

         nullify(field_1d_real)
      case (2)
         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(var_info % name)), field_2d_real, timelevel=1)

         if (.not. associated(field_2d_real)) then
            call mpp_error(FATAL,subname//'Failed to find variable "' // trim(adjustl(var_info % name)) // '"')
         end if

         if (field_2d_real % isvararray .and. associated(field_2d_real % constituentnames)) then
            allocate(var_name_list(size(field_2d_real % constituentnames)), stat=ierr)

            if (ierr /= 0) then
               call mpp_error(FATAL,subname//'Failed to allocate var_name_list')
            end if

            var_name_list(:) = field_2d_real % constituentnames(:)
         end if

         nullify(field_2d_real)
      case (3)
         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(var_info % name)), field_3d_real, timelevel=1)

         if (.not. associated(field_3d_real)) then
            call mpp_error(FATAL,subname//'Failed to find variable "' // trim(adjustl(var_info % name)) // '"')
         end if

         if (field_3d_real % isvararray .and. associated(field_3d_real % constituentnames)) then
            allocate(var_name_list(size(field_3d_real % constituentnames)), stat=ierr)

            if (ierr /= 0) then
               call mpp_error(FATAL,subname//'Failed to allocate var_name_list')
            end if

            var_name_list(:) = field_3d_real % constituentnames(:)
         end if

         nullify(field_3d_real)
      case (4)
         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(var_info % name)), field_4d_real, timelevel=1)

         if (.not. associated(field_4d_real)) then
            call mpp_error(FATAL,subname//'Failed to find variable "' // trim(adjustl(var_info % name)) // '"')
         end if

         if (field_4d_real % isvararray .and. associated(field_4d_real % constituentnames)) then
            allocate(var_name_list(size(field_4d_real % constituentnames)), stat=ierr)

            if (ierr /= 0) then
               call mpp_error(FATAL,subname//'Failed to allocate var_name_list')
            end if

            var_name_list(:) = field_4d_real % constituentnames(:)
         end if

         nullify(field_4d_real)
      case (5)
         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(var_info % name)), field_5d_real, timelevel=1)

         if (.not. associated(field_5d_real)) then
            call mpp_error(FATAL,subname//'Failed to find variable "' // trim(adjustl(var_info % name)) // '"')
         end if

         if (field_5d_real % isvararray .and. associated(field_5d_real % constituentnames)) then
            allocate(var_name_list(size(field_5d_real % constituentnames)), stat=ierr)

            if (ierr /= 0) then
               call mpp_error(FATAL,subname//'Failed to allocate var_name_list')
            end if

            var_name_list(:) = field_5d_real % constituentnames(:)
         end if

         nullify(field_5d_real)
      case default
         call mpp_error(FATAL,subname//'Unsupported variable rank ' // stringify([var_info % rank]) // &
              ' for "' // trim(adjustl(var_info % name)) // '"')
      end select
   case default
      call mpp_error(FATAL,subname//'Unsupported variable type "' // trim(adjustl(var_info % type)) // &
           '" for "' // trim(adjustl(var_info % name)) // '"')
   end select

   if (.not. allocated(var_name_list)) then
      allocate(var_name_list(1), stat=ierr)

      if (ierr /= 0) then
         call mpp_error(FATAL,subname//'Failed to allocate var_name_list')
      end if

      var_name_list(1) = var_info % name
   end if

   allocate(var_is_present(size(var_name_list)), stat=ierr)

   if (ierr /= 0) then
      call mpp_error(FATAL,subname//'Failed to allocate var_is_present')
   end if

   var_is_present(:) = .false.

   allocate(var_is_tkr_compatible(size(var_name_list)), stat=ierr)

   if (ierr /= 0) then
      call mpp_error(FATAL,subname//'Failed to allocate var_is_tkr_compatible')
   end if

   var_is_tkr_compatible(:) = .false.

   if (.not. associated(pio_file)) then
      return
   end if

   if (.not. pio_file_is_open(pio_file)) then
      return
   end if

   call mpas_log_write('Checking variable "' // trim(adjustl(var_info % name)) // &
        '" for presence and TKR compatibility')

   do i = 1, size(var_name_list)
      ! Check if the variable is present on the file.
      ierr = pio_inq_varid(pio_file, trim(adjustl(var_name_list(i))), varid)

      if (ierr /= pio_noerr) then
         cycle
      end if

      var_is_present(i) = .true.

      ! Check if the variable is "TK"R compatible between MPAS and the file.
      ierr = pio_inq_vartype(pio_file, varid, vartype)

      if (ierr /= pio_noerr) then
         cycle
      end if

      select case (trim(adjustl(var_info % type)))
      case ('character')
         if (vartype /= pio_char) then
            cycle
         end if
      case ('integer')
         if (vartype /= pio_int) then
            cycle
         end if
      case ('real')
         ! When MPAS dynamical core is compiled at single precision, pairing it with double precision input data
         ! is not allowed to prevent loss of precision.
         if (rkind == r4kind .and. vartype /= pio_real) then
            cycle
         end if

         ! When MPAS dynamical core is compiled at double precision, pairing it with single and double precision
         ! input data is allowed.
         if (rkind == r8kind .and. vartype /= pio_real .and. vartype /= pio_double) then
            cycle
         end if
      case default
         cycle
      end select

      ! Check if the variable is TK"R" compatible between MPAS and the file.
      ierr = pio_inq_varndims(pio_file, varid, varndims)

      if (ierr /= pio_noerr) then
         cycle
      end if

      if (varndims /= var_info % rank) then
         cycle
      end if

      var_is_tkr_compatible(i) = .true.
   end do

   call mpas_log_write('var_name_list = ' // stringify(var_name_list))
   call mpas_log_write('var_is_present = ' // stringify(var_is_present))
   call mpas_log_write('var_is_tkr_compatible = ' // stringify(var_is_tkr_compatible))

   call mpas_log_write(subname // ' completed')
 end subroutine dyn_mpas_check_variable_status
    
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
 !> \brief  Define the names of constituents at run-time
 !> \author Michael Duda
 !> \date   21 May 2020
 !> \details
 !>  Given an array of constituent names, which must have size equal to the number
 !>  of scalars that were set in the call to ufs_mpas_init_phase1, and given
 !>  a function to identify which scalars are moisture species, this routine defines
 !>  scalar constituents for the MPAS-A dycore.
 !>  Because the MPAS-A dycore expects all moisture constituents to appear in
 !>  a contiguous range of constituent indices, this routine may in general need
 !>  to reorder the constituents; to allow for mapping of indices between UFS
 !>  physics and the MPAS-A dycore, this routine returns index mapping arrays
 !>  mpas_from_ufs_cnst and ufs_from_mpas_cnst.
 !>
 !> \update: Dustin Swales April 2025 - Modified for use in UWM  
 !>
 !> ########################################################################################
 subroutine ufs_mpas_define_scalars(mpas_from_ufs_cnst, ufs_from_mpas_cnst, ierr)
   use mpas_derived_types, only : mpas_pool_type, field3dReal
   use mpas_pool_routines, only : mpas_pool_get_subpool, mpas_pool_get_field, &
                                  mpas_pool_get_dimension, mpas_pool_add_dimension
   use mpas_attlist,       only : mpas_add_att
   use mpas_log,           only : mpas_log_write
   use mpas_derived_types, only : MPAS_LOG_ERR
   ! FMS
   use mpp_mod,              only : FATAL, mpp_error
   
   ! Arguments
   integer, dimension(:), pointer :: mpas_from_ufs_cnst, ufs_from_mpas_cnst
   integer, intent(out) :: ierr

   ! Local variables
   character(len=*), parameter :: subname = 'ufs_mpas_subdriver::ufs_mpas_define_scalars'
   integer :: i, j, timeLevs
   integer, pointer :: num_scalars
   integer :: num_moist
   integer :: idx_passive
   type (mpas_pool_type), pointer :: statePool
   type (mpas_pool_type), pointer :: tendPool
   type (field3dReal), pointer :: scalarsField
   character(len=128) :: tempstr
   character :: moisture_char

   ierr = 0

   !
   ! Define scalars
   !
   nullify(statePool)
   call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'state', statePool)

   if (.not. associated(statePool)) then
      call mpas_log_write(trim(subname)//': ERROR: The ''state'' pool was not found.', &
                          messageType=MPAS_LOG_ERR)
      ierr = 1
      return
   end if

   nullify(num_scalars)
   call mpas_pool_get_dimension(statePool, 'num_scalars', num_scalars)

   !
   ! The num_scalars dimension should have been defined by atm_core_interface::atm_allocate_scalars, and
   ! if this dimension does not exist, something has gone wrong
   !
   if (.not. associated(num_scalars)) then
      call mpas_log_write(trim(subname)//': ERROR: The ''num_scalars'' dimension does not exist in the ''state'' pool.', &
                          messageType=MPAS_LOG_ERR)
      ierr = 1
      return
   end if

   !
   ! If at runtime there are not num_scalars names in the array of constituent names provided by UFS,
   ! something has gone wrong
   !
   if (size(constituent_name) /= num_scalars) then
      call mpas_log_write(trim(subname)//': ERROR: The number of constituent names is not equal to the num_scalars dimension', &
                          messageType=MPAS_LOG_ERR)
      call mpas_log_write('size(constituent_name) = $i, num_scalars = $i', intArgs=[size(constituent_name), num_scalars], &
                          messageType=MPAS_LOG_ERR)
      ierr = 1
      return
   end if

   !
   ! In UFS, the first scalar (if there are any) is always sphum (specific humidity); if this is not
   ! the case, something has gone wrong
   !
   if (size(constituent_name) > 0) then
      if (trim(constituent_name(1)) /= 'sphum') then
         call mpas_log_write(trim(subname)//': ERROR: The first constituent is not sphum', messageType=MPAS_LOG_ERR)
         ierr = 1
         return
      end if
   end if

   !
   ! Determine which of the constituents are moisture species
   !
   allocate(mpas_from_ufs_cnst(num_scalars), stat=ierr)
   if( ierr /= 0 ) call mpp_error(FATAL,subname//':failed to allocate mpas_from_ufs_cnst array')
   mpas_from_ufs_cnst(:) = 0
   num_moist = 0
   do i = 1, size(constituent_name)
      if (is_water_species(i)) then
         num_moist = num_moist + 1
         mpas_from_ufs_cnst(num_moist) = i
      end if
   end do

   !
   ! If UFS has no scalars, let the only scalar in MPAS be 'qv' (a moisture species)
   !
   if (num_scalars == 1 .and. size(constituent_name) == 0) then
      num_moist = 1
   end if

   !
   ! Assign non-moisture constituents to mpas_from_ufs_cnst(num_moist+1:size(constituent_name))
   !
   idx_passive = num_moist + 1
   do i = 1, size(constituent_name)

      ! If UFS constituent i is not already mapped as a moist constituent
      if (.not. is_water_species(i)) then
         mpas_from_ufs_cnst(idx_passive) = i
         idx_passive = idx_passive + 1
      end if
   end do

   !
   ! Create inverse map, ufs_from_mpas_cnst
   !
   allocate(ufs_from_mpas_cnst(num_scalars), stat=ierr)
   if( ierr /= 0 ) call mpp_error(FATAL,subname//':failed to allocate ufs_from_mpas_cnst array')
   ufs_from_mpas_cnst(:) = 0

   do i = 1, size(constituent_name)
      ufs_from_mpas_cnst(mpas_from_ufs_cnst(i)) = i
   end do

   timeLevs = 2

   do i = 1, timeLevs
      nullify(scalarsField)
      call mpas_pool_get_field(statePool, 'scalars', scalarsField, timeLevel=i)

      if (.not. associated(scalarsField)) then
         call mpas_log_write(trim(subname)//': ERROR: The ''scalars'' field was not found in the ''state'' pool', &
                             messageType=MPAS_LOG_ERR)
         ierr = 1
         return
      end if

      if (i == 1) call mpas_pool_add_dimension(statePool, 'index_qv', 1)
      scalarsField % constituentNames(1) = 'qv'
      call mpas_add_att(scalarsField % attLists(1) % attList, 'units', 'kg kg^{-1}')
      call mpas_add_att(scalarsField % attLists(1) % attList, 'long_name', 'Water vapor mixing ratio')

      do j = 2, size(constituent_name)
         scalarsField % constituentNames(j) = trim(constituent_name(mpas_from_ufs_cnst(j)))
      end do

   end do

   call mpas_pool_add_dimension(statePool, 'moist_start', 1)
   call mpas_pool_add_dimension(statePool, 'moist_end', num_moist)

   !
   ! Print a tabular summary of the mapping between constituent indices
   !
   call mpas_log_write('')
   call mpas_log_write('  i MPAS constituent mpas_from_ufs_cnst(i)       i UFS constituent  ufs_from_mpas_cnst(i)')
   call mpas_log_write('------------------------------------------     ------------------------------------------')
   do i = 1, min(num_scalars, size(constituent_name))
      if (i <= num_moist) then
         moisture_char = '*'
      else
         moisture_char = ' '
      end if
      write(tempstr, '(i3,1x,a16,1x,i18,8x,i3,1x,a16,1x,i18)') i, trim(scalarsField % constituentNames(i))//moisture_char, &
                                                               mpas_from_ufs_cnst(i), &
                                                               i, trim(constituent_name(i)), &
                                                               ufs_from_mpas_cnst(i)
      call mpas_log_write(trim(tempstr))
   end do
   call mpas_log_write('------------------------------------------     ------------------------------------------')
   call mpas_log_write('* = constituent used as a moisture species in MPAS-A dycore')
   call mpas_log_write('')


   !
   ! Define scalars_tend
   !
   nullify(tendPool)
   call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'tend', tendPool)

   if (.not. associated(tendPool)) then
      call mpas_log_write(trim(subname)//': ERROR: The ''tend'' pool was not found.', &
                          messageType=MPAS_LOG_ERR)
      ierr = 1
      return
   end if

   timeLevs = 1

   do i = 1, timeLevs
      nullify(scalarsField)
      call mpas_pool_get_field(tendPool, 'scalars_tend', scalarsField, timeLevel=i)

      if (.not. associated(scalarsField)) then
         call mpas_log_write(trim(subname)//': ERROR: The ''scalars_tend'' field was not found in the ''tend'' pool', &
                             messageType=MPAS_LOG_ERR)
         ierr = 1
         return
      end if

      if (i == 1) call mpas_pool_add_dimension(tendPool, 'index_qv', 1)
      scalarsField % constituentNames(1) = 'tend_qv'
      call mpas_add_att(scalarsField % attLists(1) % attList, 'units', 'kg m^{-3} s^{-1}')
      call mpas_add_att(scalarsField % attLists(1) % attList, 'long_name', 'Tendency of water vapor mixing ratio')

      do j = 2, size(constituent_name)
         scalarsField % constituentNames(j) = 'tend_'//trim(constituent_name(mpas_from_ufs_cnst(j)))
      end do
   end do

   call mpas_pool_add_dimension(tendPool, 'moist_start', 1)
   call mpas_pool_add_dimension(tendPool, 'moist_end', num_moist)

 end subroutine ufs_mpas_define_scalars
 
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

 !> ########################################################################################
 !>
 !> \brief  Returns global coordinate arrays
 !> \author Michael Duda
 !> \date   22 August 2019
 !> \details
 !>  This routine returns on all tasks arrays of latitude, longitude, and cell
 !>  area for all (global) cells.
 !>
 !>  It is assumed that latCellGlobal, lonCellGlobal, and areaCellGlobal have
 !>  been allocated by the caller with a size equal to the global number of
 !>  cells in the mesh.
 !>
 !> \update: Dustin Swales April 2025 - Modified for use in UWM
 !>
 !> ########################################################################################
 subroutine ufs_mpas_get_global_coords(latCellGlobal, lonCellGlobal, areaCellGlobal)
   use mpas_pool_routines, only : mpas_pool_get_subpool, mpas_pool_get_dimension, mpas_pool_get_array
   use mpas_derived_types, only : mpas_pool_type
   use mpas_kind_types,    only : RKIND
   use mpas_dmpar,         only : mpas_dmpar_sum_int, mpas_dmpar_max_real_array
   use mpp_mod,            only : FATAL, mpp_error
   real (kind=RKIND), dimension(:), intent(out) :: latCellGlobal
   real (kind=RKIND), dimension(:), intent(out) :: lonCellGlobal
   real (kind=RKIND), dimension(:), intent(out) :: areaCellGlobal

   integer :: iCell

   integer, pointer :: nCellsSolve
   integer, dimension(:), pointer :: indexToCellID

   type (mpas_pool_type), pointer :: meshPool
   integer :: nCellsGlobal,ierr

   real (kind=RKIND), dimension(:), pointer :: latCell
   real (kind=RKIND), dimension(:), pointer :: lonCell
   real (kind=RKIND), dimension(:), pointer :: areaCell
   real (kind=RKIND), dimension(:), pointer :: temp

   character(len=*), parameter :: subname = 'ufs_mpas_subdriver::ufs_mpas_get_global_coords'


   call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh', meshPool)
   call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
   call mpas_pool_get_array(meshPool, 'indexToCellID', indexToCellID)
   call mpas_pool_get_array(meshPool, 'latCell', latCell)
   call mpas_pool_get_array(meshPool, 'lonCell', lonCell)
   call mpas_pool_get_array(meshPool, 'areaCell', areaCell)

   call mpas_dmpar_sum_int(domain_ptr % dminfo, nCellsSolve, nCellsGlobal)

   ! check: size(latCellGlobal) ?= nCellsGlobal
   allocate(temp(nCellsGlobal), stat=ierr)
   if( ierr /= 0 ) call mpp_error(FATAL,subname//':failed to allocate temp array')

   !
   ! latCellGlobal
   !
   temp(:) = -huge(temp(0))
   do iCell=1,nCellsSolve
      temp(indexToCellID(iCell)) = latCell(iCell)
   end do

   call mpas_dmpar_max_real_array(domain_ptr % dminfo, nCellsGlobal, temp, latCellGlobal)

   !
   ! lonCellGlobal
   !
   temp(:) = -huge(temp(0))
   do iCell=1,nCellsSolve
      temp(indexToCellID(iCell)) = lonCell(iCell)
   end do

   call mpas_dmpar_max_real_array(domain_ptr % dminfo, nCellsGlobal, temp, lonCellGlobal)

   !
   ! areaCellGlobal
   !
   temp(:) = -huge(temp(0))
   do iCell=1,nCellsSolve
      temp(indexToCellID(iCell)) = areaCell(iCell)
   end do

   call mpas_dmpar_max_real_array(domain_ptr % dminfo, nCellsGlobal, temp, areaCellGlobal)

   deallocate(temp)

 end subroutine ufs_mpas_get_global_coords
 
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
