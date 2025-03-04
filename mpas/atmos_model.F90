! ###########################################################################################
!> \file atmos_model.F90
!>  Driver for the UFS atmospheric model with MPAS dynamical core.
!>  Contains routines to advance the atmospheric model state by one time step.
!>
! ###########################################################################################
module atmos_model_mod
  use fms2_io_mod,        only: file_exists
  use fms_mod,            only: check_nml_error
  use time_manager_mod,   only: time_type, get_time, get_date, operator(+), operator(-)
  use field_manager_mod,  only: MODEL_ATMOS
  use tracer_manager_mod, only: get_number_tracers, get_tracer_names, get_tracer_index, NO_TRACER
  use mpp_mod,            only: input_nml_file, mpp_error, FATAL
  use mpi_f08,            only: MPI_Comm, MPI_CHARACTER, MPI_INTEGER, MPI_REAL8, MPI_LOGICAL
  use MPAS_typedefs,      only: kind_phys, r8 => kind_dbl_prec
  use mpas_derived_types, only: core_type, domain_type, MPAS_Clock_type
  use atm_core_interface
  use module_mpas_config, only: pio_subsystems, fh_init
  implicit none

  private

  public :: atmos_model_init, atmos_model_end, atmos_data_type
  public :: corelist, domain_ptr

  ! #########################################################################################
  !
  ! #########################################################################################
  type atmos_data_type
     integer          :: iau_offset         ! iau running window length
     type(time_type)  :: Time               ! current time
     type(time_type)  :: Time_init          ! reference time.
  end type atmos_data_type

  type(core_type),       pointer :: corelist   => null()
  type(domain_type),     pointer :: domain_ptr => null()
  type(MPAS_Clock_type), pointer :: clock      => null()
  
contains
  ! #########################################################################################
  !
  ! Procedure to initialize UWM with MPAS dynamical core.
  !
  ! #########################################################################################
  subroutine atmos_model_init(mpicomm, master, me, time_start, time_end, total_time, calendar)
    use mpas_pool_routines,   only : mpas_pool_add_config
    use mpas_framework,       only : mpas_framework_init_phase1, mpas_framework_init_phase2
    use mpas_domain_routines, only : mpas_allocate_domain, mpas_pool_get_dimension
    use mpas_bootstrapping,   only : mpas_bootstrap_framework_phase1, mpas_bootstrap_framework_phase2
    use atm_core_interface,   only : atm_setup_core, atm_setup_domain
    use mpas_stream_inquiry,  only : MPAS_stream_inquiry_new_streaminfo
    use mpas_atm_threading,   only : mpas_atm_threading_init
    use mpas_derived_types,   only : mpas_pool_type
    use pio_types,            only : iosystem_desc_t, file_desc_t
    
    ! Inputs
    integer,        intent(in) :: time_start(6), time_end(6)
    integer,        intent(in) :: total_time, master, me
    type(MPI_Comm), intent(in) :: mpicomm
    character(17),  intent(in) :: calendar 

    ! Locals
    integer, dimension(2) :: logUnits
    integer :: i, ndate1, ndate2, tod, ierr, ntracers
    type (iosystem_desc_t), pointer :: pio_subsystem
    type (file_desc_t),     pointer :: fh_ini
    character(len=32), allocatable, target :: tracer_names(:)
    integer,           allocatable, target :: tracer_types(:)
    character(len=StrKIND) :: mesh_filename
    integer :: mesh_iotype
    type (mpas_pool_type), pointer :: state
    integer, pointer :: nVertLevels, maxEdges, maxEdges2, num_scalars
    

    allocate(corelist, stat=ierr)
    if( ierr /= 0 ) stop
    nullify(corelist%next)

    allocate(corelist%domainlist, stat=ierr)
    if( ierr /= 0 ) stop
    nullify(corelist%domainlist%next)

    domain_ptr => corelist%domainlist
    domain_ptr%core => corelist

    call mpas_allocate_domain(domain_ptr)
    domain_ptr%domainID = 0

    !
    ! Initialize MPAS infrastructure (phase 1)
    !
    call mpas_framework_init_phase1(domain_ptr%dminfo, external_comm=mpicomm)
    call atm_setup_core(corelist)
    call atm_setup_domain(domain_ptr)

    !
    ! Read MPAS namelist.
    !
    if (file_exists('input.nml')) then
       call read_mpas_namelist('input.nml', domain_ptr%configs, mpicomm, master, me)
    else
       call mpp_error(FATAL,"Cannot find namelist file: input.nml")
    end if

    ! Set forecast start time (config_start_time)
    ndate1 = time_start(1)*10000 + time_start(2)*100 + time_start(3)
    tod    = time_start(4)*3600  + time_start(5)*60  + time_start(6)
    call mpas_pool_add_config(domain_ptr%configs, 'config_start_time', date2yyyymmdd(ndate1)//'_'//sec2hms(tod))
    if (me == master) print*,'SWALES start_time ',ndate1,tod,date2yyyymmdd(ndate1)//'_'//sec2hms(tod)

    ! Set forecast end time (config_stop_time)
    ndate2 = time_end(1)*10000   + time_end(2)*100   + time_end(3)
    tod	   = time_end(4)*3600	+ time_end(5)*60    + time_end(6)
    call mpas_pool_add_config(domain_ptr%configs, 'config_stop_time', date2yyyymmdd(ndate2)//'_'//sec2hms(tod))
    if (me == master) print*,'SWALES stop_time ',ndate2,tod,date2yyyymmdd(ndate2)//'_'//sec2hms(tod)

    ! Set forecaste run time (config_run_duration) #DJS2025 this is not correct. need to fix, but works for current test.
    tod = ndate2 - ndate1 -1
    call mpas_pool_add_config(domain_ptr%configs, 'config_run_duration', trim(int2str(tod))//'_'//sec2hms(total_time))
    if (me == master) print*,'SWALES run_time ',trim(int2str(1))//'_'//sec2hms(total_time)

    ! Set other MPAS required configuration information.
    call mpas_pool_add_config(domain_ptr%configs, 'config_restart_timestamp_name', 'restart_timestamp')
    call mpas_pool_add_config(domain_ptr%configs, 'config_IAU_option',             'off')
    call mpas_pool_add_config(domain_ptr%configs, 'config_do_DAcycling',           .false.)
    call mpas_pool_add_config(domain_ptr%configs, 'config_halo_exch_method',       'mpas_halo')

    !
    ! Initialize MPAS infrastructure (phase 2)
    !
    pio_subsystem => pio_subsystems(1)
    call mpas_framework_init_phase2(domain_ptr, io_system=pio_subsystem, calendar = trim(calendar))
    print*,'SWALES POST mpas_framework_init_phase2:'

    !
    ! Before defining packages, initialize the stream inquiry instance for the domain
    !
    domain_ptr%streamInfo => MPAS_stream_inquiry_new_streaminfo()
    if (.not. associated(domain_ptr%streamInfo)) then
       call mpp_error(FATAL, 'Failed to instantiate streamInfo object for '//trim(domain_ptr%core%coreName))
    end if

    ierr = domain_ptr%core%define_packages(domain_ptr%packages)
    if (ierr /= 0) then
       call mpp_error(FATAL, 'Package definition failed for '//trim(domain_ptr%core%coreName))
    end if

    ierr = domain_ptr%core%setup_packages(domain_ptr%configs,  domain_ptr%streamInfo,       &
                                          domain_ptr%packages, domain_ptr%iocontext)
    if (ierr /= 0) then
       call mpp_error(FATAL, 'Package setup failed for '//trim(domain_ptr%core%coreName))
    end if

    ierr = domain_ptr%core%setup_decompositions(domain_ptr%decompositions)
    if (ierr /= 0) then
       call mpp_error(FATAL, 'Decomposition setup failed for '//trim(domain_ptr%core%coreName))
    end if

    ierr = domain_ptr%core%setup_clock(domain_ptr%clock, domain_ptr%configs)
    if (ierr /= 0) then
       call mpp_error(FATAL, 'Clock setup failed for '//trim(domain_ptr%core%coreName))
    end if

    !
    ! Get the number of tracers.
    !
    call get_number_tracers(MODEL_ATMOS, num_tracers=ntracers)
    allocate (tracer_names(ntracers), tracer_types(ntracers))
    do i = 1, ntracers
       call get_tracer_names(MODEL_ATMOS, i, tracer_names(i))
    enddo
    print*,'SWALES ntracers = ',ntracers
    
    !
    ! Adding a config named 'cam_pcnst' with the number of constituents will indicate to
    ! MPAS-A setup code that it is operating as a CAM dycore, and that it is necessary to
    ! allocate scalars separately from other Registry-defined fields
    !
    call mpas_pool_add_config(domain_ptr%configs, 'cam_pcnst', ntracers)

    mesh_iotype   = MPAS_IO_NETCDF  ! Not actually used
    mesh_filename = 'external mesh file'
    fh_ini => fh_init(1)
!    call mpas_bootstrap_framework_phase1(domain_ptr, mesh_filename, mesh_iotype, pio_file_desc=fh_ini)

    !
    ! Finalize the setup of blocks and fields
    !
!    call mpas_bootstrap_framework_phase2(domain_ptr, pio_file_desc=fh_ini)

    !
    ! Setup threading
    !
!    call mpas_atm_threading_init(domain_ptr%blocklist, ierr)
    if ( ierr /= 0 ) then
       call mpp_error(FATAL, 'Threading setup failed for core '//trim(domain_ptr % core % coreName))
    end if

!    !
!    ! Set up inner dimensions used by arrays in optimized dynamics routines
!    !
!    call mpas_pool_get_subpool(domain%blocklist%structs, 'state', state)
!    call mpas_pool_get_dimension(state, 'nVertLevels', nVertLevels)
!    call mpas_pool_get_dimension(state, 'maxEdges', maxEdges)
!    call mpas_pool_get_dimension(state, 'maxEdges2', maxEdges2)
!    call mpas_pool_get_dimension(state, 'num_scalars', num_scalars)
!    call mpas_atm_set_dims(nVertLevels, maxEdges, maxEdges2, num_scalars)
!
!    !
!    ! Set "local" clock to point to the clock contained in the domain type
!    !
!    clock => domain % clock
!      mpas_log_info => domain % logInfo
    
  end subroutine atmos_model_init
  
  ! #########################################################################################
  ! Procedure to finalize model.
  ! #########################################################################################
  subroutine atmos_model_end(Atmos)
    type (atmos_data_type), intent(inout) :: Atmos
  end subroutine atmos_model_end
  
  ! #########################################################################################
  ! Procedure to read MPAS namelist(s).
  !
  ! The namelist for MPAS are described in MPAS-Model/src/core_atmosphere/Registry.xml, this
  ! is also where the default values defined below originate.
  !
  ! #########################################################################################
  subroutine read_mpas_namelist(nml_file, configPool, mpicomm, master, me)
    use mpas_derived_types, only: mpas_pool_type
    use mpas_kind_types,    only: StrKIND, RKIND
    use mpas_pool_routines, only: mpas_pool_add_config

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
    ! in UFS.
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
  
end module atmos_model_mod
