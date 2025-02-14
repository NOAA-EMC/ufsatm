#define ESMF_ERR_ABORT(rc) \
if (rc /= ESMF_SUCCESS) write(0,*) 'rc=',rc,__FILE__,__LINE__; if(ESMF_LogFoundError(rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
! ###########################################################################################
!> \file module_fcst_grid_comp.F90
!>
!> ESMF forecast gridded component for MPAS Atmosphere.
!>
! ###########################################################################################
module module_fcst_grid_comp

  use mpi_f08
  use esmf
  use nuopc

  use time_manager_mod,   only: time_type, set_calendar_type, set_time,    &
                                set_date, month_name,                      &
                                operator(+), operator(-), operator (<),    &
                                operator (>), operator (/=), operator (/), &
                                operator (==), operator (*),               &
                                THIRTY_DAY_MONTHS, JULIAN, GREGORIAN,      &
                                NOLEAP, NO_CALENDAR,                       &
                                date_to_string, get_date, get_time
  use mpas_model_mod,     only: mpas_model_init, atmos_data_type
  use constants_mod,      only: constants_init
  use fms_mod,            only: error_mesg, fms_init, fms_end,             &
                                write_version_number, uppercase
  use mpp_mod,            only: mpp_init, mpp_pe, mpp_npes, mpp_root_pe, mpp_set_current_pelist,  &
                                mpp_error, FATAL, WARNING, NOTE
  use mpp_mod,            only: mpp_clock_id, mpp_clock_begin
  use sat_vapor_pres_mod, only: sat_vapor_pres_init
  use diag_manager_mod,   only: diag_manager_init, diag_manager_end,       &
                                diag_manager_set_time_end
  use fms2_io_mod,        only: FmsNetcdfFile_t, open_file, close_file, variable_exists, read_data
  use module_mpas_config, only: dt_atmos, fcst_mpi_comm, fcst_ntasks,      &
                                quilting, quilting_restart,                &
                                calendar, cpl_grid_id,                     &
                                cplprint_flag

  implicit none
  private

  !---- model defined-types ----
  type(atmos_data_type), save :: Atmos
  type(ESMF_GridComp),dimension(:),allocatable    :: fcstGridComp
  integer                                         :: ngrids, mygrid

  integer                     :: n_atmsteps

  !----- coupled model data -----
  integer :: calendar_type = -99
  integer :: date_init(6)
  integer :: numLevels     = 0
  integer :: numSoilLayers = 0
  integer :: numTracers    = 0

  integer :: frestart(999)

  integer :: mype
  integer, parameter :: iau_offset = 0

  public SetServices

contains

  ! #########################################################################################
  ! 
  ! #########################################################################################
  subroutine SetServices(fcst_comp, rc)
    type(ESMF_GridComp)  :: fcst_comp
    integer, intent(out) :: rc

    rc = ESMF_SUCCESS

    call ESMF_GridCompSetEntryPoint(fcst_comp, ESMF_METHOD_INITIALIZE, &
                                    userRoutine=fcst_initialize, phase=1, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

  end subroutine SetServices
  
  ! #########################################################################################
  ! Initialize the ESMF forecast grid component.
  ! #########################################################################################
  subroutine fcst_initialize(fcst_comp, importState, exportState, clock, rc)
    type(esmf_GridComp)                    :: fcst_comp
    type(ESMF_State)                       :: importState, exportState
    type(esmf_Clock)                       :: clock
    integer,intent(out)                    :: rc

    ! Locals
    integer                                :: i, j

    type(ESMF_VM)                          :: VM
    type(ESMF_Time)                        :: CurrTime, StartTime, StopTime
    type(ESMF_Config)                      :: cf

    real(kind=8) :: mpi_wtime, timeis
    integer :: n, k
    logical :: fexist

    integer :: initClock, total_inttime, io_unit, calendar_type_res, date_res(6), date_init_res(6)
    integer,dimension(6)                   :: date, date_end
    type(time_type)                        :: Time_init, Time, Time_step, Time_end, &
                                              Time_restart, Time_step_restart
    
    ! #######################################################################################
    ! #######################################################################################
    timeis = mpi_wtime()
    rc     = ESMF_SUCCESS

    call ESMF_VMGetCurrent(vm=vm,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_VMGet(vm=vm, localPet=mype, mpiCommunicator=fcst_mpi_comm%mpi_val, &
                    petCount=fcst_ntasks, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    if (mype == 0) write(*,*)'in fcst comp init, fcst_ntasks=',fcst_ntasks

    CF = ESMF_ConfigCreate(rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_ConfigLoadFile(config=CF ,filename='model_configure' ,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call fms_init(fcst_mpi_comm%mpi_val)
    call mpp_init()
    initClock = mpp_clock_id( 'Initialization' )
    call mpp_clock_begin (initClock) !nesting problem

    call constants_init
    call sat_vapor_pres_init

    select case( uppercase(trim(calendar)) )
    case( 'JULIAN' )
        calendar_type = JULIAN
    case( 'GREGORIAN' )
        calendar_type = GREGORIAN
    case( 'NOLEAP' )
        calendar_type = NOLEAP
    case( 'THIRTY_DAY' )
        calendar_type = THIRTY_DAY_MONTHS
    case( 'NO_CALENDAR' )
        calendar_type = NO_CALENDAR
    case default
        call mpp_error ( FATAL, 'fcst_initialize: calendar must be one of '// &
                                'JULIAN|GREGORIAN|NOLEAP|THIRTY_DAY|NO_CALENDAR.' )
    end select

    call set_calendar_type (calendar_type)
!
!-----------------------------------------------------------------------
!***  set atmos time
!-----------------------------------------------------------------------
!
    call ESMF_ClockGet(clock, CurrTime=CurrTime, StartTime=StartTime, &
                       StopTime=StopTime, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    date_init = 0
    call ESMF_TimeGet (StartTime,                      &
                       YY=date_init(1), MM=date_init(2), DD=date_init(3), &
                       H=date_init(4),  M =date_init(5), S =date_init(6), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    Time_init  = set_date (date_init(1), date_init(2), date_init(3), &
                           date_init(4), date_init(5), date_init(6))
    if (mype == 0) write(*,'(A,6I5)') 'StartTime=',date_init

    date=0
    call ESMF_TimeGet (CurrTime,                           &
                       YY=date(1), MM=date(2), DD=date(3), &
                       H=date(4),  M =date(5), S =date(6), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    Time = set_date (date(1), date(2), date(3),  &
                     date(4), date(5), date(6))
    if (mype == 0) write(*,'(A,6I5)') 'CurrTime =',date

    date_end=0
    call ESMF_TimeGet (StopTime,                                       &
                       YY=date_end(1), MM=date_end(2), DD=date_end(3), &
                       H=date_end(4),  M =date_end(5), S =date_end(6), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    Time_end   = set_date (date_end(1), date_end(2), date_end(3),  &
                           date_end(4), date_end(5), date_end(6))
    if (mype == 0) write(*,'(A,6I5)') 'StopTime =',date_end

!------------------------------------------------------------------------
!   If this is a restarted run ('INPUT/coupler.res' file exists),
!   compare date and date_init to the values in 'coupler.res'

    if (mype == 0) then
      inquire(FILE='INPUT/coupler.res', EXIST=fexist)
      if (fexist) then  ! file exists, this is a restart run

        open(newunit=io_unit, file='INPUT/coupler.res', status='old', action='read', err=998)
        read (io_unit,*,err=999) calendar_type_res
        read (io_unit,*) date_init_res
        read (io_unit,*) date_res
        close(io_unit)

        if(date_res(1) == 0 .and. date_init_res(1) /= 0) date_res = date_init_res

        if(mype == 0) write(*,'(A,6(I4))') 'INPUT/coupler.res: date_init=',date_init_res
        if(mype == 0) write(*,'(A,6(I4))') 'INPUT/coupler.res: date     =',date_res

        if (calendar_type /= calendar_type_res) then
          write(0,'(A)')      'fcst_initialize ERROR: calendar_type /= calendar_type_res'
          write(0,'(A,6(I4))')'                       calendar_type     = ', calendar_type
          write(0,'(A,6(I4))')'                       calendar_type_res = ', calendar_type_res
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
        endif

        if (.not. ALL(date_init.EQ.date_init_res)) then
          write(0,'(A)')      'fcst_initialize ERROR: date_init /= date_init_res'
          write(0,'(A,6(I4))')'                       date_init     = ', date_init
          write(0,'(A,6(I4))')'                       date_init_res = ', date_init_res
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
        endif

        if (.not. ALL(date.EQ.date_res)) then
          write(0,'(A)')      'fcst_initialize ERROR: date /= date_res'
          write(0,'(A,6(I4))')'                       date     = ', date
          write(0,'(A,6(I4))')'                       date_res = ', date_res
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
        endif

  999 continue
  998 continue

      endif ! fexist
    endif ! mype == 0

    call diag_manager_init (TIME_INIT=date)
    call diag_manager_set_time_end(Time_end)
!
    Time_step = set_time (dt_atmos,0)
    if (mype == 0) write(*,*)'time_init=', date_init,'time=',date,'time_end=',date_end,'dt_atmos=',dt_atmos

! set up forecast time array that controls when to write out restart files
    frestart = 0
    call get_time(Time_end - Time_init, total_inttime)
! set iau offset time
    Atmos%iau_offset    = iau_offset

!------ initialize component models ------

     call mpas_model_init(fcst_mpi_comm)

   end subroutine fcst_initialize

 end module  module_fcst_grid_comp
