! ###########################################################################################
!> \file atmos_model.F90
!>  Driver for the UFS ATMospheric model with MPAS dynamical core and CCPP Physics.
!>  Contains routines to advance the atmospheric model state by one forecast time step.
!>
! ###########################################################################################
module atmos_model_mod
  ! Fortran
  use mpi_f08
  ! MPAS
  use MPAS_typedefs,         only : MPAS_kind_phys => kind_phys
  use atmos_coupling_mod,    only : MPAS_statein_type, MPAS_stateout_type
  use ufs_mpas_constituents, only : constituent_name, is_water_species
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
  ! MPAS
  use mpas_log,              only : mpas_log_write
  use mpas_derived_types,    only : MPAS_LOG_CRIT
  ! FMS
  use time_manager_mod,      only : time_type, get_time, get_date, operator(+), operator(-)

  ! UFSATM
  use module_mpas_config,    only : nCellsGlobal, ic_filename, lbc_filename, nCellsSolve
  use module_mpas_config,    only : lonCell, latCell, areaCellGlobal
  use module_mpas_config,    only : pi, mpas_errfile_handle, mpas_logfile_handle
  use module_mpas_config,    only : nml_filename, nml_funit
  use module_mpas_config,    only : tracer_funit, tracer_filename
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
  real(MPAS_kind_phys) :: setupClock, atmiClock, radClock, physClock,mpasClock, mpClock, outClock

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
    use ufs_mpas_subdriver,     only : MPAS_control_type
    use ufs_mpas_subdriver,     only : ufs_mpas_init
    use ufs_mpas_io,            only : ufs_mpas_open_init, ufs_mpas_open_lbc
    use atmos_coupling_mod,     only : ufs_mpas_to_physics, ufs_mpas_grid_to_physics
    use MPAS_init,              only : MPAS_initialize

    ! Arguments
    type(atmos_control_type), intent(inout) :: Atmos
    type(time_type),          intent(in   ) :: Time_init, Time, Time_step, Time_end
    type(MPI_Comm),           intent(in   ) :: mpicomm
    character(17),            intent(in   ) :: calendar 

    ! Locals
    integer :: i, io, ierr, nConstituents, sec, iCol, mpi_size, mpi_rank
    type(MPAS_control_type) :: Cfg
    integer :: times(6), timee(6), ttime, logUnits(2), nthrds
    logical :: file_exists
    real(MPAS_kind_phys) :: start_time, stop_time
    character(len=*), parameter :: subname = 'atmos_model::atmos_model_init'
    
    ! Start timer for this procedure (init).
    start_time = MPI_Wtime()

    ! Open log files.
    logunits(1) = mpas_logfile_handle
    logunits(2) = mpas_errfile_handle
    if (Cfg % master == Cfg % me) then
       open(unit=mpas_logfile_handle, file='mpas_log.txt', action='write', status='unknown')
       open(unit=mpas_errfile_handle, file='mpas_err.txt', action='write', status='unknown')
    endif
    
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
    Cfg%master    = 0
    Cfg%mpi_comm  = mpicomm

    call MPI_Comm_rank(MPI_COMM_WORLD, Cfg%me, ierr)
    print*,'SWALES MPI size = ',Cfg%master,Cfg%me
    !
    ! Read in ATMosphere namelist (master processor only)
    !
    if ( Cfg%me == Cfg%master) then
       inquire(file = trim(nml_filename), exist=file_exists)
       if (file_exists) then
          close(nml_funit)
          open(nml_funit,file=trim(nml_filename),status='unknown')
          read(nml_funit, nml=atmos_model_nml, iostat=ierr)
          if (ierr/=0) call mpas_log_write(subname // " ERROR: When Reading in ATM Namelist",messageType=MPAS_LOG_CRIT)
       endif
    end if
    ! Broadcast ATMosphere namelist to all processors.
    call mpi_barrier(Cfg%mpi_comm, ierr)
    call mpi_bcast(regional,     1,                 MPI_LOGICAL,   Cfg%master, Cfg%mpi_comm, ierr)
    call mpi_bcast(dycore_only,  1,                 MPI_LOGICAL,   Cfg%master, Cfg%mpi_comm, ierr)
    call mpi_bcast(debug,        1,                 MPI_LOGICAL,   Cfg%master, Cfg%mpi_comm, ierr)
    call mpi_bcast(ccpp_suite,   len(ccpp_suite),   MPI_CHARACTER, Cfg%master, Cfg%mpi_comm, ierr)
    call mpi_bcast(blocksize,    1,                 MPI_INTEGER,   Cfg%master, Cfg%mpi_comm, ierr)
    call mpi_bcast(ic_filename,  len(ic_filename),  MPI_CHARACTER, Cfg%master, Cfg%mpi_comm, ierr)
    call mpi_bcast(lbc_filename, len(lbc_filename), MPI_CHARACTER, Cfg%master, Cfg%mpi_comm, ierr)
    
    !
    ! Handle constituents (scalars/tracers)
    !
    Cfg % nwat = 6
    call get_number_tracers(tracer_funit, tracer_filename, Cfg % nConstituents)
    allocate (constituent_name(Cfg % nConstituents), is_water_species(Cfg % nConstituents))
    allocate (Cfg % tracer_names(Cfg % nConstituents), Cfg % tracer_types(Cfg % nConstituents))
    call get_tracer_names(tracer_funit, tracer_filename, Cfg % nConstituents, Cfg % nwat)
    do i = 1, Cfg % nConstituents
       Cfg % tracer_names(i) = trim(constituent_name(i))
    enddo

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
    call ufs_mpas_init(Cfg, times, timee, ttime, calendar, logUnits, mpas_from_ufs_cnst, ufs_from_mpas_cnst, debug)

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
    ! Set file ID for namelist file
    Cfg%nlunit = nml_funit
    
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

    ! Read in physics namelist and allocate data containers.
    Cfg%fn_nml = nml_filename
    call MPAS_initialize(UFSATM_control, UFSATM_intdiag, UFSATM_grid, UFSATM_tbd, UFSATM_sfcprop, &
         UFSATM_statein, UFSATM_stateout, UFSATM_cldprop, UFSATM_radtend, UFSATM_coupling, Cfg)
    
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
    call ufs_mpas_to_physics(UFSATM_statein, UFSATM_sfcprop)

    ! Initialize the CCPP framework
    call CCPP_step (step="init", nblks=Atmos % nblks, ierr=ierr, dycore='mpas')
    if (ierr/=0) call mpas_log_write(subname // " ERROR: Call to CCPP init step failed",messageType=MPAS_LOG_CRIT)

    ! Initialize the CCPP physics
    call CCPP_step (step="physics_init", nblks=Atmos % nblks, ierr=ierr, dycore='mpas')
    if (ierr/=0) call mpas_log_write(subname // " ERROR: Call to CCPP physics_init step failed",messageType=MPAS_LOG_CRIT)

    ! Initialize stochastic physics pattern generation / cellular automata
    ! NOT YET IMPLEMENTED

    ! Initialize three-dimensional physics.
    ! NOT YET IMPLEMENTED
    
    stop_time = MPI_Wtime()
    atmiClock = atmiClock + (stop_time - start_time)
    !
  end subroutine atmos_model_init

  !> #########################################################################################
  !> Procedure to finalize atmospheric forecast.
  !>
  !> #########################################################################################
  subroutine atmos_model_end(Atmos)
    use ufs_mpas_tools,      only : stringify
    type (atmos_control_type), intent(inout) :: Atmos
    ! Locals
    integer :: ierr
    character(len=*), parameter :: subname = 'atmos_model::atmos_model_end'

    ! Finalize the CCPP physics.
    call CCPP_step (step="finalize", nblks=Atmos % nblks, ierr=ierr, dycore='mpas')
    if (ierr/=0) call mpas_log_write(subname // " ERROR: Call to CCPP finalize step failed",messageType=MPAS_LOG_CRIT)

    call mpas_log_write('------------------------------------------------------------------')
    call mpas_log_write('UFSATM-MPAS Timing Information (seconds):')
    call mpas_log_write('Total runtime:             '// stringify([setupClock+atmiClock+radClock+physClock+mpasClock+mpClock+outClock]))
    call mpas_log_write('Time-Step Setup:           '// stringify([setupClock]))
    call mpas_log_write('ATMosphere Initialization: '// stringify([atmiClock]))
    call mpas_log_write('CCPP Radiation:            '// stringify([radClock]))
    call mpas_log_write('CCPP Physics:              '// stringify([physClock]))
    call mpas_log_write('MPAS Dynamics:             '// stringify([mpasClock]))
    call mpas_log_write('CCPP Microphysics:         '// stringify([mpClock]))
    call mpas_log_write('MPAS Output                '// stringify([outClock]))
    call mpas_log_write('------------------------------------------------------------------')
    close(unit=mpas_logfile_handle)
    close(unit=mpas_errfile_handle)
  end subroutine atmos_model_end

  !> #########################################################################################
  !> Procedure to call atmospheric radiation and physics groups (CCPP).
  !>
  !> #########################################################################################
  subroutine atmos_model_radiation_physics(Atmos)
    use atmos_coupling_mod,     only : ufs_mpas_to_physics
    type (atmos_control_type), intent(inout) :: Atmos
    ! Locals
    integer :: ierr
    real(MPAS_kind_phys) :: start_time, stop_time
    character(len=*), parameter :: subname = 'atmos_model::atmos_model_radiation_physics'

    ! Populate physics inputs with MPAS data.
    call ufs_mpas_to_physics(UFSATM_statein, UFSATM_sfcprop)

    ! Call CCPP Timestep_initialize Group
    start_time = MPI_Wtime()
    call CCPP_step (step="timestep_init", nblks=Atmos % nblks, ierr=ierr, dycore='mpas')
    if (ierr/=0) call mpas_log_write(subname // " ERROR: Call to CCPP timestep_init step failed",messageType=MPAS_LOG_CRIT)
    stop_time = MPI_Wtime()
    setupClock = setupClock + (stop_time - start_time)

    ! Call CCPP Radiation Group
    start_time = MPI_Wtime()
    if (UFSATM_control%lsswr .or. UFSATM_control%lslwr) then
       ! DJS to GJF: If you un comment this line, you will get an error in the RRTMG radiation.
       ! Needless to say, I didn't see why, but I assume it is due to one of the many instances
       ! that we will need to identify as being FV3/MPAS specifc. Mostly in the Typedefs I suspect,
       ! but there may be interstitial schemes (NOTE that I added an new MPAS specific interstital file
       ! already, GFS_rad_time_vary.mpas.F90. I don't think it is complete.
       ! 
       !call CCPP_step (step="radiation", nblks=Atmos % nblks, ierr=ierr, dycore='mpas')
       if (ierr/=0) call mpas_log_write(subname // " ERROR: Call to CCPP radiation step failed",messageType=MPAS_LOG_CRIT)
    endif
    stop_time = MPI_Wtime()
    radClock = radClock + (stop_time - start_time)

    ! Call CCPP Physics Group
    ! NOT YET IMPLEMENTED in SDF
    start_time = MPI_Wtime()
    call CCPP_step (step="physics", nblks=Atmos % nblks, ierr=ierr, dycore='mpas')
    if (ierr/=0) call mpas_log_write(subname // " ERROR: Call to CCPP physics step failed",messageType=MPAS_LOG_CRIT)
    stop_time = MPI_Wtime()
    physClock = physClock + (stop_time - start_time)

  end subroutine atmos_model_radiation_physics

  !> #########################################################################################
  !> Procedure to call atmospheric dynamics (MPAS).
  !>
  !> #########################################################################################
  subroutine atmos_model_dynamics(Atmos)
    use ufs_mpas_subdriver, only : ufs_mpas_run
    use atmos_coupling_mod, only : ufs_physics_to_mpas
    use MPAS_init,          only : MPAS_initialize
    
    type (atmos_control_type), intent(inout) :: Atmos
    real(MPAS_kind_phys) :: start_time, stop_time
    
    ! Prepare MPAS dycore inputs with CCPP physics outputs.
    ! NOT YET IMPLEMENTED
    call ufs_physics_to_mpas()
    
    ! Call MPAS dycore
    call ufs_mpas_run(mpasClock, outClock, debug)
    
  end subroutine atmos_model_dynamics

  !> #########################################################################################
  !> Procedure to call microphysics group (CCPP).
  !>
  !> #########################################################################################
  subroutine atmos_model_microphysics(Atmos)
    use atmos_coupling_mod, only : ufs_mpas_to_microphysics, ufs_microphysics_to_mpas
    type (atmos_control_type), intent(inout) :: Atmos
    ! Locals
    integer :: ierr
    character(len=*), parameter :: subname = 'atmos_model::atmos_model_microphysics'
    real(MPAS_kind_phys) :: start_time, stop_time
 
    ! Prepare CCPP physics inputs with MPAS dycore outputs.
    ! NOT YET IMPLEMENTED
    call ufs_mpas_to_microphysics(UFSATM_statein)

    ! Call CCPP Microphysics Group
    ! NOT YET IMPLEMENTED in SDF
    start_time = MPI_Wtime()
    call CCPP_step (step="microphysics", nblks=Atmos % nblks, ierr=ierr, dycore='mpas')
    if (ierr/=0) call mpas_log_write(subname // " ERROR: Call to CCPP microphysics step failed",messageType=MPAS_LOG_CRIT)
    stop_time = MPI_Wtime()
    mpClock = mpClock + (stop_time - start_time)

    ! Call CCPP Timestep_finalize Group
    start_time = MPI_Wtime()
    call CCPP_step (step="timestep_finalize", nblks=Atmos % nblks, ierr=ierr, dycore='mpas')
    if (ierr/=0) call mpas_log_write(subname // " ERROR: Call to CCPP timestep_finalize step failed",messageType=MPAS_LOG_CRIT)
    stop_time = MPI_Wtime()
    setupClock = setupClock + (stop_time - start_time)
  
    ! Prepare MPAS dycore inputs with CCPP physics outputs.
    call ufs_microphysics_to_mpas(UFSATM_stateout)

  end subroutine atmos_model_microphysics

  !> #########################################################################################
  !> Procedure to advance the model forecast time
  !>
  !> #########################################################################################
  subroutine update_atmos_model_state(Atmos)
    type (atmos_control_type), intent(inout) :: Atmos
    character(len=*), parameter :: subname = 'atmos_model::update_atmos_model_state'

    ! Advance time
    Atmos % Time = Atmos % Time + Atmos % Time_step
  end subroutine update_atmos_model_state

  !> #########################################################################################
  !> Internal procedure to get the number of tracers (lines) in the tracer table file.
  !>
  !> #########################################################################################
  subroutine get_number_tracers(funit, fname, flines)
    integer,          intent(in)  :: funit
    character(len=*), intent(in)  :: fname
    integer,          intent(out) :: flines
    character(len=1) :: dummy
    integer :: status

    ! Get number of lines (tracers) in file
    flines = 0
    open(funit,file=trim(fname),status='unknown')
    do 
       read(funit, "(a)",iostat=status) dummy
       if (status /= 0) exit
       flines = flines + 1
    enddo
    close(funit)
  end subroutine get_number_tracers
  !> #########################################################################################
  !> Internal procedure to get tracer names from the tracer table file.
  !> ach line of the tracer table is of this format: (a10,a,a40,a,a10,a,i1)
  !>
  !> #########################################################################################
  subroutine get_tracer_names(funit, fname, ntracers, nwat)
    integer,          intent(in)  :: funit
    character(len=*), intent(in)  :: fname
    integer,          intent(in)  :: ntracers
    integer,          intent(out) :: nwat

    integer :: itracer, status
    character(len=10) :: tracer_name
    character(len=1) :: c1,c2,c3
    character(len=40) :: tracer_long_name
    character(len=10) :: tracer_unit
    integer :: tracer_type

    nwat = 0
    is_water_species(:) = .false.
    open(funit,file=trim(fname),status='unknown')
    do itracer=1,ntracers
       read(funit, "(a10,a,a40,a,a10,a,i1)",iostat=status) tracer_name,c1,tracer_long_name,c2,tracer_unit,c3,tracer_type
       constituent_name(itracer) = tracer_name
       if (tracer_type == 0) then
          is_water_species(itracer) = .true.
          nwat = nwat+1
       endif
    enddo
    close(funit)

    !Thompson/TEMPO are hardcoded to only accept 6 tracers.???
    nwat = 6
  end subroutine get_tracer_names

end module atmos_model_mod
