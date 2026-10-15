! ###########################################################################################
!> \file atmos_model.F90
!>  Driver for the UFS ATMospheric model with MPAS dynamical core and CCPP Physics.
!>  Contains routines to advance the atmospheric model state by one forecast time step.
!>
! ###########################################################################################
module atmos_model_mod
  use esmf
  use mpi_f08
  ! MPAS 
  use MPAS_typedefs,         only : MPAS_kind_phys => kind_phys
  use mpas_kind_types,       only : RKIND
  use ufs_mpas_constituents, only : constituent_name, is_water_species, constituent_type
  ! CCPP
  use CCPP_data,             only : control      => GFS_control
  use CCPP_data,             only : diag         => GFS_intdiag
  use CCPP_data,             only : interstitial => GFS_interstitial
  use CCPP_data,             only : grid         => GFS_grid
  use CCPP_data,             only : tbd          => GFS_tbd
  use CCPP_data,             only : surface      => GFS_sfcprop
  use CCPP_data,             only : statein      => GFS_statein
  use CCPP_data,             only : stateout     => GFS_stateout
  use CCPP_data,             only : cldprop      => GFS_cldprop
  use CCPP_data,             only : radiation    => GFS_radtend
  use CCPP_data,             only : coupling     => GFS_coupling
  use CCPP_driver,           only : ccpp_suite
  use CCPP_driver,           only : CCPP_step
  ! MPAS
  use mpas_log,              only : mpas_log_write
  use mpas_derived_types,    only : MPAS_LOG_CRIT
  ! UFSATM
  use module_mpas_config,    only : ic_filename, lbc_filename, oro_filename, nCellsSolve, nCellsGlobal
  use module_mpas_config,    only : stream_list_history, stream_list_restart, stream_list_diag
  use module_mpas_config,    only : lonCell, latCell, areaCellGlobal
  use module_mpas_config,    only : mpas_errfile_funit, mpas_errfilename
  use module_mpas_config,    only : mpas_logfile_funit, mpas_logfilename
  use module_mpas_config,    only : nml_filename, nml_funit
  use module_mpas_config,    only : tracer_funit, tracer_filename, constituents_file
  use module_mpas_config,    only : pi, dt_atmos, fcst_ntasks
#ifdef _OPENMP
  use omp_lib
#endif
  implicit none

  ! Day of year. Surface properties are updated daily using this time index.
  integer :: doyc

  private

  public :: dycore_only
  public :: phys_diag
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
     logical          :: isAtCapTime ! true if currTime is at the cap driverClock's currTime
     integer          :: nblks      ! Number of physics blocks.
     type(ESMF_Time)  :: CurrTime, StartTime, StopTime
     type(ESMF_TimeInterval) :: timeStep
  end type atmos_control_type

  ! Index map between MPAS tracers and UFS constituents
  integer, dimension(:), pointer :: mpas_from_ufs_cnst => null() ! indices into UFS constituent array
  ! Index map between UFS tracers and MPAS constituents
  integer, dimension(:), pointer :: ufs_from_mpas_cnst => null() ! indices into MPAS tracers array

  ! Namelist
  integer :: blocksize        = 1
  logical :: dycore_only      = .false.
  logical :: debug            = .false.
  logical :: regional         = .false.
  logical :: phys_diag        = .false.

  namelist /atmos_model_nml/ blocksize, dycore_only, phys_diag, debug, ccpp_suite, ic_filename,&
       lbc_filename, oro_filename, regional, stream_list_history, stream_list_restart,       &
       stream_list_diag, constituents_file

  ! Component Timers
  real(MPAS_kind_phys) :: setupClock = 0., atmiClock = 0., radClock = 0., physClock = 0., &
                          mpasClock = 0., mpClock = 0., outClock = 0.

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
  subroutine atmos_model_init(Atmos, mpicomm, calendar, CurrTime, StartTime, StopTime)
    use ufs_mpas_subdriver,     only : ufs_mpas_init
    use ufs_mpas_io,            only : ufs_mpas_open_init, ufs_mpas_open_lbc, ufs_mpas_open_oro
    use ufs_mpas_io,            only : ufs_mpas_read_stream_lists, ufs_mpas_landuse_read
    use ufs_mpas_io,            only : use_mpas_slopedata_read
    use atmos_coupling_mod,     only : ufs_mpas_to_physics, ufs_mpas_grid_to_physics, ufs_mpas_sfc_to_physics
    use atmos_coupling_mod,     only : ufs_mpas_landuse_update, ufs_mpas_gwd_to_physics, ufs_mpas_surface_update
    use atmos_coupling_mod,     only : ufs_mpas_reference_pressure
    use MPAS_init,              only : MPAS_initialize

    ! Arguments
    type(atmos_control_type), intent(inout) :: Atmos
    type(MPI_Comm),           intent(in   ) :: mpicomm
    character(17),            intent(in   ) :: calendar
    type(ESMF_Time),          intent(in   ) :: CurrTime, StartTime, StopTime

    ! Locals
    integer :: ierr,  rc, dt_dyn, dt_phys
    integer :: times(6), timee(6), ttime,  nthrds, me, master, nlevs
    integer :: logUnits(2)=-1
    ! Horizontal-resolution proxy (lonr/latr) and dycore-neutral reference pressure
    ! profile (p_ref) for UGWPv1 under MPAS; see the block below and MPAS_init.F90.
    ! p_ref is built in MPAS's own RKIND (matching ufs_mpas_reference_pressure) and
    ! converted to CCPP's kind_phys only at the MPAS_initialize() call.
    integer :: lonr, latr
    real(RKIND), allocatable :: p_ref(:)
    logical :: file_exists
    real(MPAS_kind_phys) :: start_time, stop_time
    integer              :: nConstituents   !< Number of constituents (tracers).
    integer              :: nwat            !< number of hydrometeors in dcyore (including water vapor)
    integer              :: bdat(8)         !< model begin date in GFS format   (same as idat)
    integer              :: cdat(8)         !< model current date in GFS format (same as jdat)
    integer              :: nblks           !< Number of data (physics) blocks
    integer, pointer     :: blksz(:)        !< Block size for  data blocking (default blksz(1)=[nCells])
    logical, parameter   :: restart=.false. !< flag whether this is a coldstart (.false.) or a warmstart/restart (.true.)
    character(len=:), pointer, dimension(:) :: input_nml_file => null()
    character(len=*), parameter :: subname = 'atmos_model::atmos_model_init'

    ! Start timer for this procedure (init).
    start_time = MPI_Wtime()

    ! Set MPI bookeeping parameters.
    master = 0
    call MPI_Comm_rank(MPI_COMM_WORLD, me, ierr)

    ! Open log files.
    if ( me == master) then
       open(newunit=mpas_logfile_funit, file=trim(mpas_logfilename), action='write', status='unknown')
       open(newunit=mpas_errfile_funit, file=trim(mpas_errfilename), action='write', status='unknown')
       logunits(1) = mpas_logfile_funit
       logunits(2) = mpas_errfile_funit
    endif

    ! Set atmospheric model time.
    Atmos % isAtCapTime = .false.
    Atmos % StartTime = StartTime
    Atmos % CurrTime  = CurrTime
    Atmos % StopTime  = StopTime

    dt_phys = real(dt_atmos)

    ! Get forecast start/stop times (year/month/day/hour/minute/second)
    call ESMF_TimeIntervalGet(StopTime-StartTime, s=ttime, rc=rc)
    call ESMF_TimeGet (StartTime, YY=times(1),MM=times(2),DD=times(3),H=times(4),M=times(5),S=times(6),rc=rc)
    call ESMF_TimeGet (StopTime,  YY=timee(1),MM=timee(2),DD=timee(3),H=timee(4),M=timee(5),S=timee(6),rc=rc)

    ! Set forecast time interval
    call ESMF_TimeIntervalSet(Atmos % timeStep, s=dt_atmos, rc=rc)

    !
    ! Read in ATMosphere namelist (master processor only)
    !
    if ( me == master) then
       inquire(file = trim(nml_filename), exist=file_exists)
       if (file_exists) then
          open(newunit=nml_funit,file=trim(nml_filename),status='unknown')
          read(nml_funit, nml=atmos_model_nml, iostat=ierr)
          if (ierr/=0) then
             print*,'ERROR: When Reading in ATM Namelist'
             stop
          endif
       endif
    end if
    ! Broadcast ATMosphere namelist to all processors.
    call mpi_barrier(mpicomm, ierr)
    call mpi_bcast(regional,            1,                        MPI_LOGICAL,   master, mpicomm, ierr)
    call mpi_bcast(dycore_only,         1,                        MPI_LOGICAL,   master, mpicomm, ierr)
    call mpi_bcast(debug,               1,                        MPI_LOGICAL,   master, mpicomm, ierr)
    call mpi_bcast(phys_diag,           1,                        MPI_LOGICAL,   master, mpicomm, ierr)
    call mpi_bcast(ccpp_suite,          len(ccpp_suite),          MPI_CHARACTER, master, mpicomm, ierr)
    call mpi_bcast(blocksize,           1,                        MPI_INTEGER,   master, mpicomm, ierr)
    call mpi_bcast(ic_filename,         len(ic_filename),         MPI_CHARACTER, master, mpicomm, ierr)
    call mpi_bcast(lbc_filename,        len(lbc_filename),        MPI_CHARACTER, master, mpicomm, ierr)
    call mpi_bcast(oro_filename,        len(oro_filename),        MPI_CHARACTER, master, mpicomm, ierr)
    call mpi_bcast(stream_list_history, len(stream_list_history), MPI_CHARACTER, master, mpicomm, ierr)
    call mpi_bcast(stream_list_restart, len(stream_list_restart), MPI_CHARACTER, master, mpicomm, ierr)
    call mpi_bcast(stream_list_diag,    len(stream_list_diag),    MPI_CHARACTER, master, mpicomm, ierr)
    call mpi_bcast(constituents_file,   len(constituents_file),   MPI_CHARACTER, master, mpicomm, ierr)

    !
    ! Handle constituents (scalars/tracers) XML
    !
    call get_tracers(constituents_file, nConstituents, nwat, debug, ierr)
    if (ierr/=0) then
       print*,'ERROR: Could not parse xml file: ',constituents_file
       stop
    end if

    ! Open (PIO) MPAS Initial Condition (IC) file.
    call ufs_mpas_open_init(ierr)
    if (ierr/=0) then
       print*,'ERROR: Could not open MPAS IC file'
       stop
    end if

    ! Open (PIO) MPAS Lateral Boundary Condition (LBC) file.
    if (regional) then
       call ufs_mpas_open_lbc(ierr)
       if (ierr/=0) then
          print*,'ERROR: Could not open MPAS LBC file'
          stop
       endif
    endif

    ! Open (PIO) MPAS orography file for GWD parameterization(s).
    if (trim(oro_filename) .ne. "none") then
       call ufs_mpas_open_oro(ierr)
       if (ierr/=0) then
          print*,'ERROR: Could not open MPAS Orography file'
          stop
       end if
    end if

    ! Call MPAS initialization.
    ! - Set up MPAS framework
    ! - Read in MPAS namelists
    ! - Set up MPAS logging
    ! - Read in static data, setup MPAS invariant stream
    ! - Setup physical constants used by MPAS dycore
    call ufs_mpas_init(me, master, mpicomm, nConstituents, nwat, times, timee, ttime, calendar,&
                       logUnits, mpas_from_ufs_cnst, ufs_from_mpas_cnst, debug, nlevs, dt_dyn)

    !
    ! Read in MPAS Stream_list file(s) (master processor only in ufs_mpas_read_stream_lists)
    !
    call ufs_mpas_read_stream_lists(me, master, mpicomm)

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

    if (dycore_only) then
       stop_time = MPI_Wtime()
       atmiClock = atmiClock + (stop_time - start_time)
       return
    end if

    ! Number of physics blocks
    Atmos % nblks = nCellsSolve / blocksize
    if (mod(nCellsSolve, blocksize) .gt. 0) Atmos % nblks = Atmos % nblks + 1

    ! Physics block sizes.
    nblks = Atmos % nblks
    allocate(blksz(Atmos % nblks))
    blksz(:) = blocksize
    blksz(Atmos % nblks) = nCellsSolve - (Atmos % nblks - 1)*blocksize

    allocate(interstitial(nthrds+1))

    ! Update time (UFS specific time formatting array)
    bdat(:) = 0
    call ESMF_TimeGet (StartTime, YY=bdat(1),MM=bdat(2),DD=bdat(3),H=bdat(5),M=bdat(6),S=bdat(7),rc=rc)
    cdat(:) = 0
    call ESMF_TimeGet (CurrTime,  YY=cdat(1),MM=cdat(2),DD=cdat(3),H=cdat(5),M=cdat(6),S=cdat(7),rc=rc)

    ! Reference vertical pressure profile for UGWPv1. Must be built before
    ! MPAS_initialize(), which is where control_initialize() consumes p_ref.
    ! Safe here: the MPAS dycore init above (ufs_mpas_init) has already populated the base
    ! state and set nlevs.
    allocate(p_ref(nlevs + 1))
    call ufs_mpas_reference_pressure(nlevs, p_ref)

    ! Horizontal-resolution proxy for the GFS physics (lonr/latr = Gaussian-grid points
    ! around the equator / pole to pole). For a quasi-uniform mesh the number of cells
    ! around the equator is sqrt(pi*N): 359 for x1.40962 (120 km). Undefined for
    ! variable-resolution meshes. Consumed by control_initialize via gnx/gny.
    lonr = nint(sqrt(pi * real(nCellsGlobal, kind=MPAS_kind_phys)))
    latr = lonr / 2
    if (me == master) call mpas_log_write('UGWP resolution proxy: lonr = $i, latr = $i from $i global cells', &
                                          intArgs=(/lonr, latr, nCellsGlobal/))

    ! Read in physics namelist and allocate data containers.
    call MPAS_initialize(control, diag, grid, tbd, surface, statein, stateout, cldprop, radiation,&
         coupling, me, master, mpicomm, nlevs, dt_dyn, dt_phys, nml_funit, nml_filename, bdat,    &
         cdat, nwat, fcst_ntasks, blksz, input_nml_file, constituent_name, constituent_type,      &
         restart, gnx=lonr, gny=latr, p_ref=real(p_ref, MPAS_kind_phys))

    ! Radiation intervals must be integer multiples of the physics time step.
    call check_radiation_interval('fhswr', control%fhswr, control%dtp)
    call check_radiation_interval('fhlwr', control%fhlwr, control%dtp)
    
    !> Read and initialize landuse fields needed by surface physics.
    call ufs_mpas_landuse_read(mpicomm, me, master)
    call ESMF_TimeGet(CurrTime, dayOfYear=doyc, rc=rc)
    call ufs_mpas_landuse_update(doyc)
    doyc = -999

    !> Read RUC LSM slope data.
    call use_mpas_slopedata_read(mpicomm, me, master)

    ! Populate UFSATM data containers with MPAS "input" stream. We need to do this becuase
    ! we are calling the physics before the MPAS dynamical core.
    call ufs_mpas_grid_to_physics(grid)
    call ufs_mpas_sfc_to_physics(surface, control)
    call ufs_mpas_gwd_to_physics(control, surface)
    call ufs_mpas_to_physics(statein, stateout, surface, radiation, mpas_from_ufs_cnst)

    ! Register CCPP
    call CCPP_step (step="register", nblks=Atmos % nblks, ierr=ierr, dycore='mpas')
    if (ierr/=0) call mpas_log_write(subname // " ERROR: Call to CCPP register step failed",messageType=MPAS_LOG_CRIT)

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
    call CCPP_step (step="final", nblks=Atmos % nblks, ierr=ierr, dycore='mpas')
    if (ierr/=0) call mpas_log_write(subname // " ERROR: Call to CCPP final step failed",messageType=MPAS_LOG_CRIT)

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
    close(unit=mpas_logfile_funit)
    close(unit=mpas_errfile_funit)
  end subroutine atmos_model_end

  !> #########################################################################################
  !> Procedure to call atmospheric radiation and physics groups (CCPP).
  !>
  !> #########################################################################################
  subroutine atmos_model_radiation_physics(Atmos)
    use atmos_coupling_mod,     only : ufs_mpas_to_physics, ufs_physics_to_mpas
    use atmos_coupling_mod,     only : ufs_mpas_landuse_update, ufs_mpas_surface_update
    use atmos_coupling_mod,     only : ufs_mpas_sfc_rad_update, ufs_mpas_ozone_init
    use ufs_mpas_io,            only : ufs_mpas_phys_diag
    type (atmos_control_type), intent(inout) :: Atmos
    ! Locals
    integer :: ierr
    real(MPAS_kind_phys) :: start_time, stop_time
    character(len=*), parameter :: subname = 'atmos_model::atmos_model_radiation_physics'
    integer :: jdat(8), rc, doy, moy, dom
    integer :: i, k
    real(MPAS_kind_phys) :: rho1, E, Qp
    real(MPAS_kind_phys), parameter :: g = 9.80665, rd = 287.05, fv = 0.6078
 
    ! Update physics time
    jdat(:) = 0
    call ESMF_TimeGet (Atmos%CurrTime, YY=jdat(1),MM=jdat(2),DD=jdat(3),H=jdat(5),M=jdat(6),S=jdat(7),rc=rc)
    control%jdat(:) = jdat(:)

    ! Update surface properties for this day? (doyc is set to -999 in init to trigger here)
    call ESMF_TimeGet(Atmos%CurrTime, dayOfYear=doy, rc=rc)
    if ((doy .gt. doyc) .or. doyc==-999) then
       call ESMF_TimeGet(Atmos%CurrTime, MM=moy, rc=rc)
       call ESMF_TimeGet(Atmos%CurrTime, DD=dom, rc=rc)
       call ufs_mpas_surface_update(dom, moy, control, surface)
       doyc = doy
    endif

    ! Populate physics inputs with MPAS data.
    call ufs_mpas_to_physics(statein, stateout, surface, radiation, mpas_from_ufs_cnst)

    ! Call CCPP Timestep_initialize Group
    start_time = MPI_Wtime()
    call CCPP_step (step="timestep_init", nblks=Atmos % nblks, ierr=ierr, dycore='mpas')
    if (ierr/=0) call mpas_log_write(subname // " ERROR: Call to CCPP timestep_init step failed",messageType=MPAS_LOG_CRIT)
    stop_time = MPI_Wtime()
    setupClock = setupClock + (stop_time - start_time)

    ! Initialize prognostic ozone tracer. Not needed for restart, since o3mmr comes from restart files.
    ! DJS2026: UFS-MPAS ICs should include ozone?
    if (control%first_time_step .and. .not. control%restart) then
       call ufs_mpas_ozone_init(control, tbd, statein, stateout, mpas_from_ufs_cnst)
    end if
    
    ! Call CCPP Radiation Group
    start_time = MPI_Wtime()
    if (control%lsswr .or. control%lslwr) then
       call CCPP_step (step="radiation", nblks=Atmos % nblks, ierr=ierr, dycore='mpas')
       if (ierr/=0) call mpas_log_write(subname // " ERROR: Call to CCPP radiation step failed",messageType=MPAS_LOG_CRIT)
    endif
    stop_time = MPI_Wtime()
    radClock = radClock + (stop_time - start_time)
    
    ! Call CCPP Physics Group
    start_time = MPI_Wtime()
    call CCPP_step (step="physics", nblks=Atmos % nblks, ierr=ierr, dycore='mpas')
    if (ierr/=0) call mpas_log_write(subname // " ERROR: Call to CCPP physics step failed",messageType=MPAS_LOG_CRIT)
    stop_time = MPI_Wtime()
    physClock = physClock + (stop_time - start_time)

    ! Update radiative surface properties from LSM
    call ufs_mpas_sfc_rad_update(surface, control)

    ! Populate MPAS pools with physics data (for diagnostics).
    call ufs_mpas_phys_diag(statein, control, radiation, diag, tbd, surface)

    ! Prepare MPAS dycore inputs with CCPP physics outputs.
    call ufs_physics_to_mpas(stateout, mpas_from_ufs_cnst)

  end subroutine atmos_model_radiation_physics

  !> #########################################################################################
  !> Procedure to call atmospheric dynamics (MPAS).
  !>
  !> #########################################################################################
  subroutine atmos_model_dynamics(Atmos)
    use ufs_mpas_subdriver, only : ufs_mpas_run

    type (atmos_control_type), intent(inout) :: Atmos

    ! Call MPAS dycore
    call ufs_mpas_run(mpasClock, outClock, debug, phys_diag)

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

    ! Prepare CCPP microphysics inputs with MPAS dycore outputs.
    call ufs_mpas_to_microphysics(stateout, statein, mpas_from_ufs_cnst)

    ! Call CCPP Microphysics Group
    ! NOT YET IMPLEMENTED in SDF
    start_time = MPI_Wtime()
    call CCPP_step (step="microphysics", nblks=Atmos % nblks, ierr=ierr, dycore='mpas')
    if (ierr/=0) call mpas_log_write(subname // " ERROR: Call to CCPP microphysics step failed",messageType=MPAS_LOG_CRIT)
    stop_time = MPI_Wtime()
    mpClock = mpClock + (stop_time - start_time)

    ! Call CCPP Timestep_final Group
    start_time = MPI_Wtime()
    call CCPP_step (step="timestep_final", nblks=Atmos % nblks, ierr=ierr, dycore='mpas')
    if (ierr/=0) call mpas_log_write(subname // " ERROR: Call to CCPP timestep_final step failed",messageType=MPAS_LOG_CRIT)
    stop_time = MPI_Wtime()
    setupClock = setupClock + (stop_time - start_time)

    ! Prepare MPAS dycore inputs with CCPP physics outputs.
    call ufs_microphysics_to_mpas(stateout, mpas_from_ufs_cnst)

    control % first_time_step = .false.

  end subroutine atmos_model_microphysics

  !> #########################################################################################
  !> Procedure to advance the model forecast time
  !>
  !> #########################################################################################
  subroutine update_atmos_model_state(Atmos)
    type (atmos_control_type), intent(inout) :: Atmos
    character(len=*), parameter :: subname = 'atmos_model::update_atmos_model_state'

    ! Advance time
    Atmos % CurrTime = Atmos % CurrTime + Atmos % TimeStep

  end subroutine update_atmos_model_state

  !> ########################################################################################
  !> Procedure to parse constituents.XML file.
  !>
  !> ########################################################################################
  subroutine get_tracers(constituents_file, nvars, nvarsw, debug, ierr)
    use iso_c_binding
    use mpas_log,  only : mpas_log_write
    use mpas_derived_types,  only : MPAS_LOG_CRIT
    use ezxml_mod
    implicit none
    character(len=*), intent(in) :: constituents_file
    logical, intent(in) :: debug
    integer, intent(out) :: ierr, nvars, nvarsw
    type(xml_node) :: root, variable
    integer :: ivar, ivarw
    character(len=200) :: name, standard_name, units, type, kind, allocatble, dimensions
    character(len=200) :: water_species
    character(len=*), parameter :: subname = 'atmos_model:get_tracers'

    ! Initialize
    ierr   = 0
    nvars  = 0
    nvarsw = 0

    ! Open and parse the XML file.
    root = xml_parse_file(trim(constituents_file))
    if (.not. xml_is_valid(root)) then
       call mpas_log_write(subname // " could not find xml file: "//trim(constituents_file), messageType=MPAS_LOG_CRIT)
       ierr = -1
       return
    end if

    !
    ! First parse to get number of constituents and number of water-species.
    !
    variable = xml_child(root, "var")
    do while (xml_is_valid(variable))
       water_species = xml_attr(variable, "water_species")
       variable = xml_next(variable)

       ! Increment counters
       nvars = nvars + 1
       if (trim(water_species) == "yes") then
          nvarsw = nvarsw + 1
       end if
    end do
    if (debug) then
       print*, " Number of constituents  : ",nvars
       print*, " Number of water_species : ",nvarsw
    end if

    ! Allocate space for tracer names/attributes
    allocate(constituent_name(nvars))
    allocate(constituent_type(nvars))
    allocate(is_water_species(nvars))

    !
    ! Second parse, save fields and attributes
    !
    ivar  = 0
    ivarw = 0
    is_water_species(:) = .false.
    constituent_type(:) = 0
    ! Move to the <var> child element. Process each <var> in <constituents>.
    variable = xml_child(root, "var")
    do while (xml_is_valid(variable))

       name          = xml_attr(variable, "name")
       standard_name = xml_attr(variable, "standard_name")
       units         = xml_attr(variable, "units")
       type          = xml_attr(variable, "type")
       kind          = xml_attr(variable, "kind")
       allocatble    = xml_attr(variable, "allocatable")
       dimensions    = xml_attr(variable, "dimensions")
       water_species = xml_attr(variable, "water_species")
       if (debug) then
          print '(A,A)', "  name          : ", trim(name)
          print '(A,A)', "  units         : ", trim(units)
          print '(A,A)', "  std_name      : ", trim(standard_name)
          print '(A,A)', "  type          : ", trim(type)
          print '(A,A)', "  kind          : ", trim(kind)
          print '(A,A)', "  alloc         : ", trim(allocatble)
          print '(A,A)', "  dimensions    : ", trim(dimensions)
          print '(A,A)', "  water_species : ", trim(water_species)
       end if

       ! Move to next variable in file
       variable = xml_next(variable)

       ! Save
       ivar = ivar + 1
       if (trim(water_species) == "yes") then
          is_water_species(ivar) = .true.
          constituent_type(ivar) = 1
       end if
       constituent_name(ivar) = trim(name)

    end do

    ! Clean up memory
    call xml_free(root)

  end subroutine get_tracers

  !> ########################################################################################
  !> Abort if a radiation interval (s) is not an exact multiple of the physics step (s).
  !>
  !> ########################################################################################
  subroutine check_radiation_interval(name, interval, dtp)
    character(len=*),     intent(in) :: name
    real(MPAS_kind_phys), intent(in) :: interval, dtp
    real(MPAS_kind_phys) :: nsteps
    character(len=*), parameter :: subname = 'atmos_model::check_radiation_interval'
 
    if (interval <= 0.0_MPAS_kind_phys .or. dtp <= 0.0_MPAS_kind_phys) return  ! not set / not used
 
    nsteps = interval / dtp
    if (abs(nsteps - real(nint(nsteps), MPAS_kind_phys)) > 1.0e-6_MPAS_kind_phys * max(1.0_MPAS_kind_phys, nsteps)) then
       call mpas_log_write(subname // ': '// trim(name) // ' = $r s is not a multiple of the physics '// &
            'time step dt_atmos = $r s ($r steps). Radiation would be called every $i steps ($r s) '// &
            'while the zenith angle is averaged over '// trim(name) // '. Choose '// trim(name) // &
            ' as an integer multiple of dt_atmos.', &
            realArgs=(/real(interval,RKIND), real(dtp,RKIND), real(nsteps,RKIND), &
                       real(max(1,nint(nsteps))*dtp,RKIND)/), &
            intArgs=(/max(1,nint(nsteps))/), messageType=MPAS_LOG_CRIT)
    end if
  end subroutine check_radiation_interval
end module atmos_model_mod
