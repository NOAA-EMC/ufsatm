! ###########################################################################################
!> \file atmos_model.F90
!>  Driver for the UFS atmospheric model with MPAS dynamical core.
!>  Contains routines to advance the atmospheric model state by one time step.
!>
! ###########################################################################################
module atmos_model_mod
  ! Fortran
  use mpi_f08,               only : MPI_Comm, MPI_CHARACTER, MPI_INTEGER, MPI_REAL8, MPI_LOGICAL
  ! MPAS
  use MPAS_typedefs,         only : MPAS_control_type, MPAS_kind_phys => kind_phys
  ! CCPP
  use MPAS_data,             only : MPAS_statein, MPAS_stateint, MPAS_stateout
  use CCPP_data,             only : GFS_control, GFS_intdiag, ccpp_suite
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
  use module_mpas_config,    only : pio_numiotasks, nCellsGlobal
  implicit none

  private

  public :: atmos_model_init, atmos_model_end, atmos_model_radiation_physics, atmos_data_type,&
       atmos_model_microphysics, atmos_model_dynamics

  !> #########################################################################################
  !> Type containing information on MPAS enabled UFSATM forecast.
  !>
  !> #########################################################################################
  type atmos_data_type
     type(time_type)  :: Time       ! current time
     type(time_type)  :: Time_step  ! atmospheric time step.
     type(time_type)  :: Time_init  ! reference time.
     integer          :: nblks      ! Number of physics blocks.
  end type atmos_data_type

  ! Namelist
  integer :: blocksize    = 1
  logical :: dycore_only  = .false.
  logical :: debug        = .false.

  namelist /atmos_model_nml/ blocksize, dycore_only, debug, ccpp_suite

  ! Component Timers
  integer :: setupClock, radClock, physClock, mpasClock, mpClock, atmiClock

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
    use ufs_mpas_subdriver, only : ufs_mpas_init_phase1, ufs_mpas_init_phase2, ufs_mpas_dyn_set
    use ufs_mpas_subdriver, only : ufs_mpas_open_init, ufs_mpas_read_init
    use MPAS_init,          only : MPAS_initialize

    ! Arguments
    type(atmos_data_type), intent(inout) :: Atmos
    type(time_type),       intent(in   ) :: Time_init, Time, Time_step, Time_end
    type(MPI_Comm),        intent(in   ) :: mpicomm
    character(17),         intent(in   ) :: calendar 

    ! Locals
    integer :: i, io, ierr, nConstituents, sec
    type(MPAS_control_type) :: Cfg
    integer :: times(6), timee(6), ttime, logUnits(2)
    
    ! Set up timers
    setupClock = mpp_clock_id( 'Time-Step Setup       ', flags=clock_flag_default, grain=CLOCK_COMPONENT )
    atmiClock  = mpp_clock_id( 'ATMosphere Setup      ', flags=clock_flag_default, grain=CLOCK_COMPONENT )
    radClock   = mpp_clock_id( 'Radiation             ', flags=clock_flag_default, grain=CLOCK_COMPONENT )
    physClock  = mpp_clock_id( 'Physics               ', flags=clock_flag_default, grain=CLOCK_COMPONENT )
    mpasClock  = mpp_clock_id( 'MPAS Dycore           ', flags=clock_flag_default, grain=CLOCK_COMPONENT )
    mpClock    = mpp_clock_id( 'Microphysics          ', flags=clock_flag_default, grain=CLOCK_COMPONENT )

    ! Start timer for this procedure (init).
    call mpp_clock_begin(atmiClock)

    ! Set model time
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

    ! Get the number of tracers.
    call get_number_tracers(MODEL_ATMOS, num_tracers=Cfg % nConstituents)
    allocate (Cfg % tracer_names(Cfg % nConstituents), Cfg % tracer_types(Cfg % nConstituents))
    do i = 1, Cfg % nConstituents
       call get_tracer_names(MODEL_ATMOS, i, Cfg % tracer_names(i))
    enddo
    ! DJS2025: There are 9 tracers, but only 6 are water. How do we get to 6?
    ! With FV3, this is set during dycore initialization.
    Cfg % nwat = 6

    ! Open (PIO) MPAS IC data file.
    call ufs_mpas_open_init()

    ! Call MPAS initialization phase 1.
    ! - Set up MPAS framework
    ! - Read in MPAS namelists
    ! - Set up MPAS logging
    ! - Read in static data, setup MPAS invariant stream
    ! - Setup physical constants used by MPAS dycore
    logUnits(1) = stdout()
    logUnits(2) = stdlog()
    call ufs_mpas_init_phase1(Cfg, times, timee, ttime, calendar, logUnits)

    ! Create MPAS data containers
    ! (Associate UWM data containers with MPAS pool variables)
    call ufs_mpas_dyn_set(MPAS_Statein, MPAS_Stateout)

    ! Read in MPAS IC data. Populate UWM data containers and MPAS "input" stream.
    call ufs_mpas_read_init(MPAS_Statein)

    ! Complete the MPAS dycore initialization.
    ! - Set up threading.
    ! - Call MPAS core_atmosphere init.
    call ufs_mpas_init_phase2(Cfg, MPAS_Statein)

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

    ! Set file ID for log file
    Cfg%nlunit = stdlog()
    
    ! Number of physics blocks
    Atmos % nblks = nCellsGlobal / blocksize
    if (mod(nCellsGlobal, blocksize) .gt. 0) Atmos % nblks = Atmos % nblks + 1
    
    ! Physics block sizes.
    Cfg % nblks = Atmos % nblks
    allocate(Cfg % blksz(Atmos % nblks))
    Cfg % blksz(:) = blocksize
    Cfg % blksz(Atmos % nblks) = nCellsGlobal - (Atmos % nblks - 1)*blocksize
    
    ! ### Remove when implementing GFS physics and partitioning GFS data containers ###
    ! The hybrid-sigma coordinates are included in a GFS data container that is set up
    ! during the Physics init/namelist-read, below when calling MPAS_initialize(). These
    ! fields should be in FV3 data container, allowing for both MPAS and FV3 to share
    ! the GFS data container.
    !allocate(Cfg%ak(Cfg%levs+1))
    !allocate(Cfg%bk(Cfg%levs+1))
    !Cfg%ak(:) = 0.0
    !Cfg%bk(:) = 0.0
    
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
    call MPAS_initialize(GFS_control, GFS_intdiag, MPAS_Statein, MPAS_Stateint, MPAS_Stateout, Cfg)

    ! Initialize the CCPP framework
    call CCPP_step (step="init", nblks=Atmos % nblks, ierr=ierr, dynamics='mpas')
    if (ierr/=0)  call mpp_error(FATAL, 'Call to CCPP init step failed')

    ! Initialize the CCPP physics
    call CCPP_step (step="physics_init", nblks=Atmos % nblks, ierr=ierr, dynamics='mpas')
    if (ierr/=0)  call mpp_error(FATAL, 'Call to CCPP physics_init step failed')

    ! Initialize stochastic physics pattern generation / cellular automata
    ! NOT YET IMPLEMENTED

    call mpp_clock_end(atmiClock)
    !
  end subroutine atmos_model_init

  !> #########################################################################################
  !> Procedure to finalize model.
  !>
  !> #########################################################################################
  subroutine atmos_model_end(Atmos)
    type (atmos_data_type), intent(inout) :: Atmos
    ! Locals
    integer :: ierr

    ! Finalize the CCPP physics.
    call CCPP_step (step="finalize", nblks=Atmos % nblks, ierr=ierr, dynamics='mpas')
    if (ierr/=0)  call mpp_error(FATAL, 'Call to CCPP finalize step failed')

  end subroutine atmos_model_end

  !> #########################################################################################
  !> Procedure to call atmospheric radiation and physics (CCPP).
  !>
  !> #########################################################################################
  subroutine atmos_model_radiation_physics(Atmos)
    type (atmos_data_type), intent(inout) :: Atmos
    ! Locals
    integer :: ierr

    ! Call CCPP Timestep_initialize Group
    call mpp_clock_begin(setupClock)
    call CCPP_step (step="timestep_init", nblks=Atmos % nblks, ierr=ierr, dynamics='mpas')
    if (ierr/=0)  call mpp_error(FATAL, 'Call to CCPP timestep_init step failed')
    call mpp_clock_end(setupClock)
    
    ! Call CCPP Group Radiation
    call mpp_clock_begin(radClock)
    if (GFS_control%lsswr .or. GFS_control%lslwr) then
       !call CCPP_step (step="radiation", nblks=Atmos % nblks, ierr=ierr, dynamics='mpas')
       if (ierr/=0)  call mpp_error(FATAL, 'Call to CCPP radiation step failed')
    endif
    call mpp_clock_end(radClock)

    ! Call CCPP Group Physics
    call mpp_clock_begin(physClock)
    call CCPP_step (step="physics", nblks=Atmos % nblks, ierr=ierr, dynamics='mpas')
    if (ierr/=0)  call mpp_error(FATAL, 'Call to CCPP physics step failed')
    call mpp_clock_end(physClock)

    ! Call CCPP Timestep_finalize Group
    call mpp_clock_begin(setupClock)
    call CCPP_step (step="timestep_finalize", nblks=Atmos % nblks, ierr=ierr, dynamics='mpas')
    if (ierr/=0)  call mpp_error(FATAL, 'Call to CCPP timestep_finalize step failed')
    call mpp_clock_end(setupClock)
    
  end subroutine atmos_model_radiation_physics

  !> #########################################################################################
  !> Procedure to call atmospheric dynamics (MPAS).
  !>
  !> #########################################################################################
  subroutine atmos_model_dynamics(Atmos)
    use ufs_mpas_subdriver, only: ufs_mpas_run
    use MPAS_init,          only: MPAS_initialize
    
    type (atmos_data_type), intent(inout) :: Atmos
    
    ! Call MPAS dycore
    call mpp_clock_begin(mpasClock)
    call ufs_mpas_run(MPAS_Statein, MPAS_Stateout)
    call mpp_clock_end(mpasClock)
    
  end subroutine atmos_model_dynamics

  !> #########################################################################################
  !> Procedure to call microphysics (CCPP)
  !>
  !> #########################################################################################
  subroutine atmos_model_microphysics(Atmos)
    type (atmos_data_type), intent(inout) :: Atmos
    ! Locals
    integer :: ierr
    
    ! Call CCPP Group Microphysics
    call mpp_clock_begin(mpClock)
    call CCPP_step (step="microphysics", nblks=Atmos % nblks, ierr=ierr, dynamics='mpas')
    if (ierr/=0)  call mpp_error(FATAL, 'Call to CCPP microphysics step failed')
    call mpp_clock_end(mpClock)

   end subroutine atmos_model_microphysics
  !
end module atmos_model_mod
