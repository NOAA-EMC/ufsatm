! ###########################################################################################
!> \file atmos_model.F90
!>  Driver for the UFS ATMospheric model with MPAS dynamical core and CCPP Physics.
!>  Contains routines to advance the atmospheric model state by one forecast time step.
!>
! ###########################################################################################
module atmos_model_mod
  ! Fortran
  use mpi_f08,               only : MPI_Comm, MPI_CHARACTER, MPI_INTEGER, MPI_REAL8, MPI_LOGICAL
  ! MPAS
  use MPAS_typedefs,         only : MPAS_kind_phys => kind_phys
  use atmos_coupling_mod,    only : MPAS_statein_type, MPAS_stateout_type
  ! CCPP
  use CCPP_data,             only : UFSATM_control      => GFS_control
  use CCPP_data,             only : UFSATM_intdiag      => GFS_intdiag
  use CCPP_data,             only : UFSATM_interstitial => GFS_interstitial
  use CCPP_data,             only : UFSATM_grid         => GFS_grid
  use CCPP_data,             only : UFSATM_tbd          => GFS_tbd
  use CCPP_data,             only : UFSATM_sfcprop      => GFS_sfcprop
  use CCPP_data,             only : UFSATM_statein      => GFS_statein
  use CCPP_data,             only : UFSATM_stateout     => GFS_stateout
  use CCPP_data,             only : UFSATM_cldprop      => GFS_cldprop
  use CCPP_data,             only : UFSATM_radtend      => GFS_radtend
  use CCPP_data,             only : UFSATM_coupling     => GFS_coupling
  use CCPP_data,             only : ccpp_suite
  use CCPP_driver,           only : CCPP_step
  ! FMS
  use time_manager_mod,      only : time_type, get_time, get_date, operator(+), operator(-)
  use field_manager_mod,     only : MODEL_ATMOS
  use tracer_manager_mod,    only : get_number_tracers, get_tracer_names, get_tracer_index
  use fms_mod,               only : check_nml_error
  use fms2_io_mod,           only : file_exists
  use mpp_mod,               only : input_nml_file, mpp_error, FATAL
  use mpp_mod,               only : mpp_pe, mpp_root_pe, mpp_clock_id, mpp_clock_begin
  use mpp_mod,               only : mpp_clock_end, CLOCK_COMPONENT, MPP_CLOCK_SYNC
  use fms_mod,               only : clock_flag_default
  use fms_mod,               only : stdlog
  use mpp_mod,               only : stdout
  ! UFSATM
  use module_mpas_config,    only : nCellsGlobal, ic_filename, lbc_filename, nCellsSolve
  use module_mpas_config,    only : lonCell, latCell, areaCellGlobal
  use module_mpas_config,    only : pi
  use mod_ufsatm_util,       only : get_atmos_tracer_types
#ifdef _OPENMP
  use omp_lib
#endif
  implicit none

  private

  public :: atmos_control_type
  public :: atmos_model_init
  public :: atmos_model_end
  public :: atmos_model_radiation_physics
  public :: atmos_model_microphysics
  public :: atmos_model_dynamics
  public :: update_atmos_model_state

  !> #########################################################################################
  !> Type containing information on MPAS enabled UFSATM forecast.
  !>
  !> #########################################################################################
  type atmos_control_type
     type(time_type)  :: Time       ! current time
     type(time_type)  :: Time_step  ! atmospheric time step.
     type(time_type)  :: Time_init  ! reference time.
     logical          :: isAtCapTime ! true if currTime is at the cap driverClock's currTime 
     integer          :: nblks      ! Number of physics blocks.
  end type atmos_control_type
  
  ! Index map between MPAS tracers and UFS constituents
  integer, dimension(:), pointer :: mpas_from_ufs_cnst => null() ! indices into UFS constituent array
  ! Index map between UFS tracers and MPAS constituents
  integer, dimension(:), pointer :: ufs_from_mpas_cnst => null() ! indices into MPAS tracers array  
  
  ! Namelist
  integer :: blocksize    = 1
  logical :: dycore_only  = .false.
  logical :: debug        = .false.
  logical :: regional     = .false.

  namelist /atmos_model_nml/ blocksize, dycore_only, debug, ccpp_suite, ic_filename, lbc_filename, &
       regional

  ! Component Timers
  integer :: setupClock, radClock, physClock, mpasClock, mpClock, atmiClock

  ! DJS2025: For UFS WM RTs unitl output is setup for MPAS.
  integer, parameter :: mpas_logfile_handle = 42323

  type(MPAS_statein_type)  :: MPAS_statein
  type(MPAS_stateout_type) :: MPAS_stateout
contains
  !> #########################################################################################
  !> Procedure to initialize UWM ATMosphere with MPAS dynamical core.
  !>
  !> - Read in ATMosphere namelist
  !> - Initialize MPAS framework
  !> - Read in MPAS namelist
  !> - Initialize MPAS dynamical core
  !>   - Read in MPAS initial conditions
  !> - Read in physics namelist
  !> - Initialize CCPP framework
  !> - Initialize CCPP Physics
  !>
  !> #########################################################################################
  subroutine atmos_model_init(Atmos, Time_init, Time, Time_end, Time_step, mpicomm, calendar)
    use ufs_mpas_subdriver, only : MPAS_control_type
    use ufs_mpas_subdriver, only : ufs_mpas_init
    use ufs_mpas_subdriver, only : ufs_mpas_open_init, ufs_mpas_open_lbc
    use ufs_mpas_module,    only : constituent_name, is_water_species
    use atmos_coupling_mod, only : ufs_mpas_to_physics, ufs_mpas_grid_to_physics
    use MPAS_init,          only : MPAS_initialize

    ! Arguments
    type(atmos_control_type), intent(inout) :: Atmos
    type(time_type),          intent(in   ) :: Time_init, Time, Time_step, Time_end
    type(MPI_Comm),           intent(in   ) :: mpicomm
    character(17),            intent(in   ) :: calendar 

    ! Locals
    integer :: i, io, ierr, nConstituents, sec, iCol
    type(MPAS_control_type) :: Cfg
    integer :: times(6), timee(6), ttime, logUnits(2), nthrds
    
    ! Set up timers
    setupClock = mpp_clock_id( 'Time-Step Setup       ', flags=clock_flag_default, grain=CLOCK_COMPONENT )
    atmiClock  = mpp_clock_id( 'ATMosphere Setup      ', flags=clock_flag_default, grain=CLOCK_COMPONENT )
    radClock   = mpp_clock_id( 'Radiation             ', flags=clock_flag_default, grain=CLOCK_COMPONENT )
    physClock  = mpp_clock_id( 'Physics               ', flags=clock_flag_default, grain=CLOCK_COMPONENT )
    mpasClock  = mpp_clock_id( 'MPAS Dycore           ', flags=clock_flag_default, grain=CLOCK_COMPONENT )
    mpClock    = mpp_clock_id( 'Microphysics          ', flags=clock_flag_default, grain=CLOCK_COMPONENT )

    ! Start timer for this procedure (init).
    call mpp_clock_begin(atmiClock)

    ! Set atmospheric model time.
    Atmos % isAtCapTime = .false.
    Atmos % Time_init = Time_init
    Atmos % Time      = Time
    Atmos % Time_step = Time_step
    
    call get_time (Atmos % Time_step, sec)
    Cfg%dt_phys   = real(sec)
    
    ! Get forecast start/stop times (year/month/day/hour/minute/second)
    call get_date(Time_init,times(1),times(2),times(3),times(4),times(5),times(6))
    call get_date(Time_end, timee(1),timee(2),timee(3),timee(4),timee(5),timee(6))
    call get_time(Time_end - Time_init, ttime)
    
    ! Set MPI bookeeping parameters.
    Cfg%me        = mpp_pe()
    Cfg%master    = mpp_root_pe()
    Cfg%mpi_comm  = mpicomm
    
    ! Read in ATMosphere namelist.
    if (file_exists('input.nml')) then
       read(input_nml_file, nml=atmos_model_nml, iostat=io)
       ierr = check_nml_error(io, 'atmos_model_nml')
    endif

    !
    ! Handle constituents (scalars/tracers)
    !

    ! Get constituent name(s) and type(s).
    ! Active constituents are defined in the FMS "field_table".
    call get_number_tracers(MODEL_ATMOS, num_tracers=Cfg % nConstituents)
    allocate (Cfg % tracer_names(Cfg % nConstituents), Cfg % tracer_types(Cfg % nConstituents))
    do i = 1, Cfg % nConstituents
       call get_tracer_names(MODEL_ATMOS, i, Cfg % tracer_names(i))
    enddo
    call get_atmos_tracer_types(Cfg % tracer_types)

    ! Get number of water species.
    ! DJS Asks? With FV3, this is set during dycore initialization. How do we get this information
    ! here? Does MPAS have a routine for this?
    !
    ! It would be simple, albeit not the most elegant thing, but we could create a simple routine
    ! that has a list of "known MPAS water species" and compare each "tracer_name" to that.
    ! A more robust solution IMO would be to quiery the field table entries for a "water-species"
    ! attribute, or something along those lines. Actually, I think this is straightforward if we
    ! extend ../ufsatm_util.F90.
    
    !
    ! From field_tables:
    ! For RRFS   MPAS we have: 11 water tracers (ql,qc,qi,qr,qs,qg,nc,nc,ni,nr,ng)
    !                           2 prog. tracers (o3,sgs-tke)
    ! For GFSv17 MPAS we have:  6 water species (ql,qc,qi,qr,qs,qg)
    !                           4 prog. tracers (o3,sgs-tke,cld_amt,sigma_b)
    Cfg % nwat = 6

    call get_number_tracers(MODEL_ATMOS, num_tracers=Cfg % nConstituents)
    allocate (constituent_name(Cfg % nConstituents), is_water_species(Cfg % nConstituents))
    do i = 1, Cfg % nConstituents
       call get_tracer_names(MODEL_ATMOS, i, constituent_name(i))
    enddo
    is_water_species(:) = .false.
    is_water_species(1:Cfg % nwat) = .true.

    ! Open (PIO) MPAS Initial Condition (IC) file.
    call ufs_mpas_open_init()

    ! Open (PIO) MPAS Lateral Boundary Condition (LBC) file.
    if (regional) then
       call ufs_mpas_open_lbc()
    endif

    ! Call MPAS initialization.
    ! - Set up MPAS framework
    ! - Read in MPAS namelists
    ! - Set up MPAS logging
    ! - Read in static data, setup MPAS invariant stream
    ! - Setup physical constants used by MPAS dycore
    logUnits(1) = stdout()
    logUnits(2) = stdlog()

    ! DJS2025: This is for UWM RT logging only. Can be removed when MPAS output is added.
    if (Cfg % master == Cfg % me) then
       open(unit=mpas_logfile_handle, file='mpas_log.txt', action='write', status='unknown')
       logunits(1) = mpas_logfile_handle
       logunits(2) = mpas_logfile_handle
    endif

    call ufs_mpas_init(Cfg, times, timee, ttime, calendar, logUnits, mpas_from_ufs_cnst, ufs_from_mpas_cnst)

    !> #########################################################################################
    !> #########################################################################################
    !> END MPAS DYCORE INITIALIZATION
    !> #########################################################################################
    !> #########################################################################################

    !> #########################################################################################
    !> #########################################################################################
    !> BEGIN CCPP PHYSICS INITIALIZATION
    !> #########################################################################################
    !> #########################################################################################
#ifdef _OPENMP
    nthrds = omp_get_max_threads()
#else
    nthrds = 1
#endif
    ! Set file ID for log file
    Cfg%nlunit = stdlog()
    
    ! Number of physics blocks
    Atmos % nblks = nCellsSolve / blocksize
    if (mod(nCellsSolve, blocksize) .gt. 0) Atmos % nblks = Atmos % nblks + 1

    ! Physics block sizes.
    Cfg % nblks = Atmos % nblks
    allocate(Cfg % blksz(Atmos % nblks))
    Cfg % blksz(:) = blocksize
    Cfg % blksz(Atmos % nblks) = nCellsSolve - (Atmos % nblks - 1)*blocksize

    allocate(UFSATM_interstitial(nthrds+1))
    
    ! Update time (UFS specific time formatting array)
    Cfg%bdat(:) = 0
    call get_date (Time_init, Cfg%bdat(1), Cfg%bdat(2), Cfg%bdat(3), Cfg%bdat(5), Cfg%bdat(6), Cfg%bdat(7))
    Cfg%cdat(:) = 0
    call get_date (Time,      Cfg%cdat(1), Cfg%cdat(2), Cfg%cdat(3), Cfg%cdat(5), Cfg%cdat(6), Cfg%cdat(7))

    ! Allocate required to work around GNU compiler bug 100886 https://gcc.gnu.org/bugzilla/show_bug.cgi?id=100886
    allocate(Cfg%input_nml_file, mold=input_nml_file)
    Cfg%input_nml_file  => input_nml_file
    Cfg%fn_nml='using internal file'

    ! Read in physics namelist and allocate data containers.
    call MPAS_initialize(UFSATM_control, UFSATM_intdiag, UFSATM_grid, UFSATM_tbd, UFSATM_sfcprop, &
         UFSATM_statein, UFSATM_cldprop, UFSATM_radtend, UFSATM_coupling, Cfg)
    
    call ufs_mpas_grid_to_physics(UFSATM_grid)

    ! Populate UFSATM data containers with MPAS "input" stream. We need to do this becuase
    ! we are calling the physics before the MPAS dynamical core.
    !
    ! DJS to GJF: See fcst_run_phase_1 in module_fcst_grid_comp.F90. That is where we call the
    ! "pieces" of the Atmospheric timestep defined below.
    ! Since we are calling the radiation/physics first, we need to take the MPAS Initial state
    ! and map it to the physics data containers (e.g. Typdefs). We will use a similar routine
    ! in a different "piece" later, but copying the Updated state from the dycore before calling
    ! the microphsyics.
    !
    call ufs_mpas_to_physics(UFSATM_statein)

    ! Initialize the CCPP framework
    call CCPP_step (step="init", nblks=Atmos % nblks, ierr=ierr, dycore='mpas')
    if (ierr/=0)  call mpp_error(FATAL, 'Call to CCPP init step failed')

    ! Initialize the CCPP physics
    call CCPP_step (step="physics_init", nblks=Atmos % nblks, ierr=ierr, dycore='mpas')
    if (ierr/=0)  call mpp_error(FATAL, 'Call to CCPP physics_init step failed')

    ! Initialize stochastic physics pattern generation / cellular automata
    ! NOT YET IMPLEMENTED

    ! Initialize three-dimensional physics.
    ! NOT YET IMPLEMENTED
    
    call mpp_clock_end(atmiClock)
    !
  end subroutine atmos_model_init

  !> #########################################################################################
  !> Procedure to finalize atmospheric forecast.
  !>
  !> #########################################################################################
  subroutine atmos_model_end(Atmos)
    type (atmos_control_type), intent(inout) :: Atmos
    ! Locals
    integer :: ierr

    close(unit=mpas_logfile_handle)

    ! Finalize the CCPP physics.
    call CCPP_step (step="finalize", nblks=Atmos % nblks, ierr=ierr, dycore='mpas')
    if (ierr/=0)  call mpp_error(FATAL, 'Call to CCPP finalize step failed')

  end subroutine atmos_model_end

  !> #########################################################################################
  !> Procedure to call atmospheric radiation and physics groups (CCPP).
  !>
  !> #########################################################################################
  subroutine atmos_model_radiation_physics(Atmos)
    type (atmos_control_type), intent(inout) :: Atmos
    ! Locals
    integer :: ierr

    ! Call CCPP Timestep_initialize Group
    call mpp_clock_begin(setupClock)
    call CCPP_step (step="timestep_init", nblks=Atmos % nblks, ierr=ierr, dycore='mpas')
    if (ierr/=0)  call mpp_error(FATAL, 'Call to CCPP timestep_init step failed')
    call mpp_clock_end(setupClock)

    ! Call CCPP Radiation Group
    call mpp_clock_begin(radClock)
    if (UFSATM_control%lsswr .or. UFSATM_control%lslwr) then
       ! DJS to GJF: If you un comment this line, you will get an error in the RRTMG radiation.
       ! Needless to say, I didn't see why, but I assume it is due to one of the many instances
       ! that we will need to identify as being FV3/MPAS specifc. Mostly in the Typedefs I suspect,
       ! but there may be interstitial schemes (NOTE that I added an new MPAS specific interstital file
       ! already, GFS_rad_time_vary.mpas.F90. I don't think it is complete.
       ! 
       !call CCPP_step (step="radiation", nblks=Atmos % nblks, ierr=ierr, dycore='mpas')
       if (ierr/=0)  call mpp_error(FATAL, 'Call to CCPP radiation step failed')
    endif
    call mpp_clock_end(radClock)

    ! Call CCPP Physics Group
    ! NOT YET IMPLEMENTED in SDF
    call mpp_clock_begin(physClock)
    call CCPP_step (step="physics", nblks=Atmos % nblks, ierr=ierr, dycore='mpas')
    if (ierr/=0)  call mpp_error(FATAL, 'Call to CCPP physics step failed')
    call mpp_clock_end(physClock)

  end subroutine atmos_model_radiation_physics

  !> #########################################################################################
  !> Procedure to call atmospheric dynamics (MPAS).
  !>
  !> #########################################################################################
  subroutine atmos_model_dynamics(Atmos)
    use ufs_mpas_subdriver, only : ufs_mpas_run
    use atmos_coupling_mod, only : ufs_microphysics_to_mpas
    use MPAS_init,          only : MPAS_initialize
    
    type (atmos_control_type), intent(inout) :: Atmos

    ! Prepare MPAS dycore inputs with CCPP physics outputs.
    ! NOT YET IMPLEMENTED
    call ufs_microphysics_to_mpas(UFSATM_stateout)
    
    ! Call MPAS dycore
    call mpp_clock_begin(mpasClock)
    call ufs_mpas_run()
    call mpp_clock_end(mpasClock)
    
  end subroutine atmos_model_dynamics

  !> #########################################################################################
  !> Procedure to call microphysics group (CCPP).
  !>
  !> #########################################################################################
  subroutine atmos_model_microphysics(Atmos)
    use atmos_coupling_mod, only : ufs_mpas_to_microphysics
    type (atmos_control_type), intent(inout) :: Atmos
    ! Locals
    integer :: ierr

    ! Prepare CCPP physics inputs with MPAS dycore outputs.
    ! NOT YET IMPLEMENTED
    call ufs_mpas_to_microphysics(UFSATM_statein)

    ! Call CCPP Microphysics Group
    ! NOT YET IMPLEMENTED in SDF
    call mpp_clock_begin(mpClock)
    call CCPP_step (step="microphysics", nblks=Atmos % nblks, ierr=ierr, dycore='mpas')
    if (ierr/=0)  call mpp_error(FATAL, 'Call to CCPP microphysics step failed')
    call mpp_clock_end(mpClock)

    ! Call CCPP Timestep_finalize Group
    call mpp_clock_begin(setupClock)
    call CCPP_step (step="timestep_finalize", nblks=Atmos % nblks, ierr=ierr, dycore='mpas')
    if (ierr/=0)  call mpp_error(FATAL, 'Call to CCPP timestep_finalize step failed')
    call mpp_clock_end(setupClock)

  end subroutine atmos_model_microphysics

  !> #########################################################################################
  !> Procedure to advance the model forecast time
  !>
  !> #########################################################################################
  subroutine update_atmos_model_state(Atmos)
    type (atmos_control_type), intent(inout) :: Atmos

    ! Advance time
    Atmos % Time = Atmos % Time + Atmos % Time_step
  end subroutine update_atmos_model_state
  
end module atmos_model_mod
