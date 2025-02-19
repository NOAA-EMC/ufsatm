! ###########################################################################################
!> \file mpas_model.F90
!>  Driver for the UFS atmospheric model with MPAS dynamical core.
!>  Contains routines to advance the atmospheric model state by one time step.
!>
! ###########################################################################################
module mpas_model_mod
  use fms2_io_mod,        only: file_exists
  use fms_mod,            only: check_nml_error
  use time_manager_mod,   only: time_type, get_time, get_date, operator(+), operator(-)
  use mpp_mod,            only: input_nml_file
  use mpi_f08,            only: MPI_Comm
  use MPAS_typedefs,      only: kind_phys, r8 => kind_dbl_prec
  use mpas_derived_types, only: core_type, domain_type, MPAS_Clock_type
  use atm_core_interface

  implicit none

  private

  public :: mpas_model_init, mpas_model_end, atmos_data_type
  public :: corelist, domain_ptr
  
  type atmos_data_type
     integer                       :: axes(4)            ! axis indices (returned by diag_manager) for the atmospheric grid
                                                         ! (they correspond to the x, y, pfull, phalf axes)
     integer, pointer              :: pelist(:) =>null() ! pelist where atmosphere is running.
     integer                       :: layout(2)          ! computer task laytout
     logical                       :: regional           ! true if domain is regional
     logical                       :: nested             ! true if there is a nest
     logical                       :: moving_nest_parent ! true if this grid has a moving nest child
     logical                       :: is_moving_nest     ! true if this is a moving nest grid
     logical                       :: isAtCapTime        ! true if currTime is at the cap driverClock's currTime
     integer                       :: ngrids             !
     integer                       :: mygrid             !
     integer                       :: mlon, mlat
     integer                       :: iau_offset         ! iau running window length
     logical                       :: pe                 ! current pe.
     real(kind=kind_phys), pointer, dimension(:)     :: ak, bk
     real(kind=kind_phys), pointer, dimension(:,:)   :: lon_bnd  => null() ! local longitude axis grid box corners in radians.
     real(kind=kind_phys), pointer, dimension(:,:)   :: lat_bnd  => null() ! local latitude axis grid box corners in radians.
     real(kind=kind_phys), pointer, dimension(:,:)   :: lon      => null() ! local longitude axis grid box centers in radians.
     real(kind=kind_phys), pointer, dimension(:,:)   :: lat      => null() ! local latitude axis grid box centers in radians.
     real(kind=kind_phys), pointer, dimension(:,:)   :: dx, dy
     real(kind=kind_phys), pointer, dimension(:,:)   :: area
     real(kind=kind_phys), pointer, dimension(:,:,:) :: layer_hgt, level_hgt
     !type(domain2d)                :: domain             ! domain decomposition
     !type(domain2d)                :: domain_for_read    ! domain decomposition
     type(time_type)               :: Time               ! current time
     type(time_type)               :: Time_step          ! atmospheric time step.
     type(time_type)               :: Time_init          ! reference time.
     !type(grid_box_type)           :: grid               ! hold grid information needed for 2nd order conservative flux exchange
     !type(GFS_externaldiag_type), pointer, dimension(:) :: Diag
  end type atmos_data_type

  type(core_type),       pointer :: corelist   => null()
  type(domain_type),     pointer :: domain_ptr => null()
  type(MPAS_Clock_type), pointer :: clock      => null()
  
contains
  ! #########################################################################################
  ! Procedure to initialize UWM with MPAS dynamical core.
  ! #########################################################################################
  subroutine mpas_model_init(mpicomm, time_start, time_end, total_time)
    use mpas_pool_routines, only : mpas_pool_add_config

    ! Inputs
    integer, intent(in) :: time_start(6), time_end(6)
    integer, intent(in) :: total_time
    type(MPI_Comm), intent(in) :: mpicomm

    ! Locals
    integer, dimension(2) :: logUnits
    integer :: ndate, tod

    print*,'SWALES mpas_init_phase1 A'
    ! Set up MPAS framework/infrastructure.
    call mpas_init_phase1(mpicomm, logUnits)
    print*,'SWALES mpas_init_phase1 B'

    ! Read MPAS namelist.
    if (file_exists('input.nml')) then
       call read_mpas_namelist('input.nml', domain_ptr%configs)
    end if
    print*,'SWALES mpas_init_phase1 C'

    ! Set config_start_time
    ndate = time_start(1)*10000 + time_start(2)*100 + time_start(3)
    tod   = time_start(4)*3600  + time_start(5)*60  + time_start(6)
    call mpas_pool_add_config(domain_ptr%configs, 'config_start_time', date2yyyymmdd(ndate)//'_'//sec2hms(tod))

    ! Set config_stop_time
    ndate = time_end(1)*10000   + time_end(2)*100   + time_end(3)
    tod	  = time_end(4)*3600	+ time_end(5)*60    + time_end(6)
    call mpas_pool_add_config(domain_ptr%configs, 'config_stop_time', date2yyyymmdd(ndate)//'_'//sec2hms(tod))

    ! Set config_run_duration
    call mpas_pool_add_config(domain_ptr%configs, 'config_run_duration', trim(int2str(1))//'_'//sec2hms(total_time))
    print*,'SWALES mpas_init_phase1 D'

  end subroutine mpas_model_init
  
  ! #########################################################################################
  ! Procedure to finalize model.
  ! #########################################################################################
  subroutine mpas_model_end(Atmos)
    type (atmos_data_type), intent(inout) :: Atmos
  end subroutine mpas_model_end
  
  ! #########################################################################################
  ! Procedure to read MPAS namelist(s).
  !
  ! The namelist for MPAS are described in MPAS-Model/src/core_atmosphere/Registry.xml, this
  ! is also where the default values below come from.
  !
  ! #########################################################################################
  subroutine read_mpas_namelist(nml_file, configPool)
    use mpas_derived_types, only: mpas_pool_type
    use mpas_kind_types,    only: StrKIND
    use mpas_pool_routines, only: mpas_pool_add_config

    character(len=*), intent(in) :: nml_file
    type(mpas_pool_type), intent(inout) :: configPool

    ! Namelist
    character (len=StrKIND) :: mpas_time_integration = 'SRK3'
    integer                 :: mpas_time_integration_order = 2
    real(r8)                :: mpas_dt = 720.0_r8
    logical                 :: mpas_split_dynamics_transport = .true.
    integer                 :: mpas_number_of_sub_steps = 2
    integer                 :: mpas_dynamics_split_steps = 3
    real(r8)                :: mpas_h_mom_eddy_visc2 = 0.0_r8
    real(r8)                :: mpas_h_mom_eddy_visc4 = 0.0_r8
    real(r8)                :: mpas_v_mom_eddy_visc2 = 0.0_r8
    real(r8)                :: mpas_h_theta_eddy_visc2 = 0.0_r8
    real(r8)                :: mpas_h_theta_eddy_visc4 = 0.0_r8
    real(r8)                :: mpas_v_theta_eddy_visc2 = 0.0_r8
    character (len=StrKIND) :: mpas_horiz_mixing = '2d_smagorinsky'
    real(r8)                :: mpas_len_disp = 120000.0_r8
    real(r8)                :: mpas_visc4_2dsmag = 0.05_r8
    real(r8)                :: mpas_del4u_div_factor = 10.0_r8
    integer                 :: mpas_w_adv_order = 3
    integer                 :: mpas_theta_adv_order = 3
    integer                 :: mpas_scalar_adv_order = 3
    integer                 :: mpas_u_vadv_order = 3
    integer                 :: mpas_w_vadv_order = 3
    integer                 :: mpas_theta_vadv_order = 3
    integer                 :: mpas_scalar_vadv_order = 3
    logical                 :: mpas_scalar_advection = .true.
    logical                 :: mpas_positive_definite = .false.
    logical                 :: mpas_monotonic = .true.
    real(r8)                :: mpas_coef_3rd_order = 0.25_r8
    real(r8)                :: mpas_smagorinsky_coef = 0.125_r8
    logical                 :: mpas_mix_full = .true.
    real(r8)                :: mpas_epssm = 0.1_r8
    real(r8)                :: mpas_smdiv = 0.1_r8
    real(r8)                :: mpas_apvm_upwinding = 0.5_r8
    logical                 :: mpas_h_ScaleWithMesh = .true.
    real(r8)                :: mpas_zd = 22000.0_r8
    real(r8)                :: mpas_xnutr = 0.2_r8
    real(r8)                :: mpas_cam_coef = 0.0_r8
    integer                 :: mpas_cam_damping_levels = 0
    logical                 :: mpas_rayleigh_damp_u = .true.
    real(r8)                :: mpas_rayleigh_damp_u_timescale_days = 5.0_r8
    integer                 :: mpas_number_rayleigh_damp_u_levels = 3
    logical                 :: mpas_apply_lbcs = .false.
    logical                 :: mpas_jedi_da = .false.
    character (len=StrKIND) :: mpas_block_decomp_file_prefix = 'x1.40962.graph.info.part.'
    logical                 :: mpas_do_restart = .false.
    logical                 :: mpas_print_global_minmax_vel = .true.
    logical                 :: mpas_print_detailed_minmax_vel = .false.
    logical                 :: mpas_print_global_minmax_sca = .false.
    
    namelist /nhyd_model/ &
             mpas_time_integration, &
             mpas_time_integration_order, &
             mpas_dt, &
             mpas_split_dynamics_transport, &
             mpas_number_of_sub_steps, &
             mpas_dynamics_split_steps, &
             mpas_h_mom_eddy_visc2, &
             mpas_h_mom_eddy_visc4, &
             mpas_v_mom_eddy_visc2, &
             mpas_h_theta_eddy_visc2, &
             mpas_h_theta_eddy_visc4, &
             mpas_v_theta_eddy_visc2, &
             mpas_horiz_mixing, &
             mpas_len_disp, &
             mpas_visc4_2dsmag, &
             mpas_del4u_div_factor, &
             mpas_w_adv_order, &
             mpas_theta_adv_order, &
             mpas_scalar_adv_order, &
             mpas_u_vadv_order, &
             mpas_w_vadv_order, &
             mpas_theta_vadv_order, &
             mpas_scalar_vadv_order, &
             mpas_scalar_advection, &
             mpas_positive_definite, &
             mpas_monotonic, &
             mpas_coef_3rd_order, &
             mpas_smagorinsky_coef, &
             mpas_mix_full, &
             mpas_epssm, &
             mpas_smdiv, &
             mpas_apvm_upwinding, &
             mpas_h_ScaleWithMesh

    namelist /damping/ &
             mpas_zd, &
             mpas_xnutr, &
             mpas_cam_coef, &
             mpas_cam_damping_levels, &
             mpas_rayleigh_damp_u, &
             mpas_rayleigh_damp_u_timescale_days, &
             mpas_number_rayleigh_damp_u_levels

    namelist /limited_area/ &
             mpas_apply_lbcs

    namelist /assimilation/ &
             mpas_jedi_da

    namelist /decomposition/ &
             mpas_block_decomp_file_prefix

    namelist /restart/ &
             mpas_do_restart

    namelist /printout/ &
             mpas_print_global_minmax_vel, &
             mpas_print_detailed_minmax_vel, &
             mpas_print_global_minmax_sca

    ! These configuration parameters must be set in the MPAS configPool, but can't
    ! be changed in UFS.
    integer                :: config_num_halos = 2
    integer                :: config_number_of_blocks = 0
    logical                :: config_explicit_proc_decomp = .false.
    character(len=StrKIND) :: config_proc_decomp_file_prefix = 'graph.info.part'

    ! Locals
    integer :: ierr, io

    ! Read in namelists...
    if (file_exists(nml_file)) then
       read(input_nml_file, nml=nhyd_model, iostat=io)
       ierr = check_nml_error(io, 'nhyd_model')
    endif

  end subroutine read_mpas_namelist

  ! #########################################################################################
  !
  ! #########################################################################################
  subroutine mpas_init_phase1(mpicomm, logUnits)
    use mpas_domain_routines, only : mpas_allocate_domain
    use mpas_framework,       only : mpas_framework_init_phase1
    use atm_core_interface,   only : atm_setup_core, atm_setup_domain

    ! Inputs
    type(MPI_Comm), intent(in) :: mpicomm
    integer, dimension(2), intent(in) :: logUnits
    ! Locals
    integer :: ierr

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

    ! Initialize MPAS infrastructure
    call mpas_framework_init_phase1(domain_ptr%dminfo, external_comm=mpicomm)
    call atm_setup_core(corelist)
    call atm_setup_domain(domain_ptr)

  end subroutine mpas_init_phase1

  
  ! #########################################################################################
  !
  ! #########################################################################################
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
80  format(i4.4,'-',i2.2,'-',i2.2)

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
80  format(i2.2,':',i2.2,':',i2.2)

  end function sec2hms
  
  ! #########################################################################################
  !
  ! #########################################################################################
  character(len=10) function int2str(n)

    ! return default integer as a left justified string
    ! arguments
    integer, intent(in) :: n
    !----------------------------------------------------------------------------

    write(int2str,'(i0)') n
     
  end function int2str
  
end module mpas_model_mod
