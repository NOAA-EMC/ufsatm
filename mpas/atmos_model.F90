! ###########################################################################################
!> \file atmos_model.F90
!>  Driver for the UFS atmospheric model with MPAS dynamical core.
!>  Contains routines to advance the atmospheric model state by one time step.
!>
! ###########################################################################################
module atmos_model_mod
  ! Fortran
  use mpi_f08,               only: MPI_Comm, MPI_CHARACTER, MPI_INTEGER, MPI_REAL8, MPI_LOGICAL
  ! MPAS
  use MPAS_typedefs,         only: MPAS_init_type, MPAS_kind_phys => kind_phys
  ! CCPP
  use CCPP_data,             only: GFS_control, MPAS_statein, MPAS_stateint, MPAS_stateout
  use CCPP_data,             only: GFS_intdiag, ccpp_suite
  use CCPP_driver,           only: CCPP_step
  ! FMS
  use time_manager_mod,      only: time_type, get_time, get_date, operator(+), operator(-)
  use field_manager_mod,     only: MODEL_ATMOS
  use tracer_manager_mod,    only: get_number_tracers, get_tracer_names, get_tracer_index
  use fms_mod,               only: check_nml_error
  use fms2_io_mod,           only: file_exists
  use block_control_mod,     only: block_control_type, define_blocks_packed
  use mpp_mod,               only: input_nml_file, mpp_error, FATAL
  use mpp_mod,               only: mpp_pe, mpp_root_pe, mpp_clock_id, mpp_clock_begin
  use mpp_mod,               only: mpp_clock_end, CLOCK_COMPONENT, MPP_CLOCK_SYNC
  use fms_mod,               only: clock_flag_default
  use fms_mod,               only: stdlog
  use mpp_mod,               only: stdout
  implicit none

  private

  public :: atmos_model_init, atmos_model_end, atmos_model_radiation_physics, atmos_data_type,&
       atmos_model_microphysics, atmos_model_dynamics

  ! #########################################################################################
  !
  ! #########################################################################################
  type atmos_data_type
     integer          :: iau_offset         ! iau running window length
     type(time_type)  :: Time               ! current time
     type(time_type)  :: Time_step          ! atmospheric time step.
     type(time_type)  :: Time_init          ! reference time.
  end type atmos_data_type

  ! Namelist
  integer :: blocksize    = 1
  logical :: dycore_only  = .false.
  logical :: debug        = .false.

  namelist /atmos_model_nml/ blocksize, dycore_only, debug, ccpp_suite

  type (block_control_type), target   :: Atm_block

  ! Component Timers
  integer :: setupClock, radClock, physClock, mpasClock, mpClock, atmiClock
contains
  ! #########################################################################################
  !
  ! Procedure to initialize UWM ATMosphere with MPAS dynamical core.
  !
  ! - Read in ATMosphere namelist
  ! - Initialize MPAS framework
  ! - Read in MPAS namelist
  ! - Initialize MPAS dynamical core
  !   - Read in MPAS initial conditions
  !   - Read in MPAS mesh
  ! - Read in physics namelist
  ! - Initialize CCPP framework
  ! - Initialize CCPP Physics
  !
  ! #########################################################################################
  subroutine atmos_model_init(Atmos, Time_init, Time, Time_end, Time_step, mpicomm, calendar)
    use ufs_mpas_subdriver, only: ufs_mpas_init
    use MPAS_init,          only: MPAS_initialize

    ! Inputs
    type(atmos_data_type), intent(inout) :: Atmos
    type(time_type), intent(in) :: Time_init, Time, Time_step, Time_end
    type(MPI_Comm), intent(in) :: mpicomm
    character(17),  intent(in) :: calendar 

    ! Locals
    integer :: i, io, ierr, nConstituents, sec
    type(MPAS_init_type)  :: Init
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

    ! Get forecast start/stop times (year/month/day/hour/minute/second)
    call get_date(Time_init,times(1),times(2),times(3),times(4),times(5),times(6))
    call get_date(Time_end, timee(1),timee(2),timee(3),timee(4),timee(5),timee(6))
    call get_time(Time_end - Time_init, ttime)
    
    ! Parameters needed for initialization
    Init%me        = mpp_pe()
    Init%master    = mpp_root_pe()
    Init%mpi_comm  = mpicomm
    Init%dt_phys   = real(sec)
    
    ! Read in ATMosphere namelist.
    if (file_exists('input.nml')) then
       read(input_nml_file, nml=atmos_model_nml, iostat=io)
       ierr = check_nml_error(io, 'atmos_model_nml')
    endif

    ! Get the number of tracers.
    call get_number_tracers(MODEL_ATMOS, num_tracers=Init%nConstituents)
    allocate (Init%tracer_names(Init%nConstituents), Init%tracer_types(Init%nConstituents))
    do i = 1, Init%nConstituents
       call get_tracer_names(MODEL_ATMOS, i, Init%tracer_names(i))
    enddo
    ! DJS2025: There are 9 tracers, but only 6 are water. How do we get to 6?
    ! With FV3, this is set during dycore initialization.
    Init%nwat = 6

    ! Initialize the MPAS dynamical core. Read in MPAS dycore namelist.
    logUnits(1) = stdout()
    logUnits(2) = stdlog()
    call ufs_mpas_init(Init, times, timee, ttime, calendar, logUnits)

    ! Set domain description.
    ! ### Work in Progress. Set to 1 for MPAS initialization testing. ###
    Atm_block%nblks = 1
    allocate(Init%blksz(Atm_block%nblks))
    allocate(Init%ak(Init%levs+1))
    allocate(Init%bk(Init%levs+1))
    Init%ak(:)    = 0.0
    Init%bk(:)    = 0.0
    Init%blksz(:) = blocksize
    Init%nlunit = stdlog()

    ! Update time (UFS specific time formatting array)
    Init%bdat(:) = 0
    call get_date (Time_init, Init%bdat(1), Init%bdat(2), Init%bdat(3), Init%bdat(5),       &
         Init%bdat(6), Init%bdat(7))
    Init%cdat(:) = 0
    call get_date (Time,      Init%cdat(1), Init%cdat(2), Init%cdat(3), Init%cdat(5),       &
         Init%cdat(6), Init%cdat(7))

    ! Allocate required to work around GNU compiler bug 100886 https://gcc.gnu.org/bugzilla/show_bug.cgi?id=100886
    allocate(Init%input_nml_file, mold=input_nml_file)
    Init%input_nml_file  => input_nml_file
    Init%fn_nml='using internal file'

    ! Read in physics namelist and allocate data containers.
    call MPAS_initialize(GFS_control, GFS_intdiag, MPAS_Statein, MPAS_Stateint, MPAS_Stateout, Init)

    ! Initialize the CCPP framework
    call CCPP_step (step="init", nblks=Atm_block%nblks, ierr=ierr)
    if (ierr/=0)  call mpp_error(FATAL, 'Call to CCPP init step failed')

    ! Initialize the CCPP physics
    call CCPP_step (step="physics_init", nblks=Atm_block%nblks, ierr=ierr)
    if (ierr/=0)  call mpp_error(FATAL, 'Call to CCPP physics_init step failed')

    ! Initialize stochastic physics pattern generation / cellular automata
    ! NOT YET IMPLEMENTED

    call mpp_clock_end(atmiClock)
    !
  end subroutine atmos_model_init

  ! #########################################################################################
  !
  ! Procedure to finalize model.
  !
  ! #########################################################################################
  subroutine atmos_model_end(Atmos)
    type (atmos_data_type), intent(inout) :: Atmos
    ! Locals
    integer :: ierr

    ! Finalize the CCPP physics.
    call CCPP_step (step="finalize", nblks=Atm_block%nblks, ierr=ierr)
    if (ierr/=0)  call mpp_error(FATAL, 'Call to CCPP finalize step failed')

  end subroutine atmos_model_end

  ! #########################################################################################
  !
  ! Procedure to call atmospheric radiation and physics.
  !
  ! #########################################################################################
  subroutine atmos_model_radiation_physics(Atmos)
    type (atmos_data_type), intent(inout) :: Atmos
    ! Locals
    integer :: ierr

    ! Call CCPP Timestep_initialize Group
    call mpp_clock_begin(setupClock)
    call CCPP_step (step="timestep_init", nblks=Atm_block%nblks, ierr=ierr)
    if (ierr/=0)  call mpp_error(FATAL, 'Call to CCPP timestep_init step failed')
    call mpp_clock_end(setupClock)
    
    ! Call CCPP Group Radiation
    call mpp_clock_begin(radClock)
    if (GFS_control%lsswr .or. GFS_control%lslwr) then
       !call CCPP_step (step="radiation", nblks=Atm_block%nblks, ierr=ierr)
       if (ierr/=0)  call mpp_error(FATAL, 'Call to CCPP radiation step failed')
    endif
    call mpp_clock_end(radClock)

    ! Call CCPP Group Physics
    call mpp_clock_begin(physClock)
    !call CCPP_step (step="phys_ps", nblks=Atm_block%nblks, ierr=ierr)
    if (ierr/=0)  call mpp_error(FATAL, 'Call to CCPP physics step failed')
    call mpp_clock_end(physClock)

    ! Call CCPP Timestep_finalize Group
    call mpp_clock_begin(setupClock)
    !call CCPP_step (step="timestep_finalize", nblks=Atm_block%nblks, ierr=ierr)
    if (ierr/=0)  call mpp_error(FATAL, 'Call to CCPP timestep_finalize step failed')
    call mpp_clock_end(setupClock)
    
  end subroutine atmos_model_radiation_physics

  ! #########################################################################################
  !
  ! Procedure to call atmospheric dynamics
  !
  ! #########################################################################################
  subroutine atmos_model_dynamics(Atmos)
    type (atmos_data_type), intent(inout) :: Atmos
    
    ! Call MPAS dycore
    call mpp_clock_begin(mpasClock)
    !!! NOT YET IMPLEMENTED!!!
    call mpp_clock_end(mpasClock)
    
  end subroutine atmos_model_dynamics

  ! #########################################################################################
  !
  ! Procedure to call microphysics
  !
  ! #########################################################################################
  subroutine atmos_model_microphysics(Atmos)
    type (atmos_data_type), intent(inout) :: Atmos
    ! Locals
    integer :: ierr
    
    ! Call CCPP Group Microphysics
    call mpp_clock_begin(mpClock)
    call CCPP_step (step="microphysics", nblks=Atm_block%nblks, ierr=ierr)
    if (ierr/=0)  call mpp_error(FATAL, 'Call to CCPP microphysics step failed')
    call mpp_clock_end(mpClock)

   end subroutine atmos_model_microphysics
  !
end module atmos_model_mod
