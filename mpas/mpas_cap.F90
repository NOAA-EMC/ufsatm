! ###########################################################################################
!> \file mpas_cap.F90
!>
!> This file contains the NUOPC Cap for the UWM Atmosphere with the MPAS dynamical core.
!>
! ###########################################################################################
module mpasatm_cap_mod

  use ESMF
  use NUOPC
  use NUOPC_Model,            only: model_routine_SS => SetServices,                        &
                                    SetVM,                                                  &
                                    routine_Run,                                            &
                                    label_Advertise,                                        &
                                    label_RealizeProvided,                                  &
                                    label_Advance,                                          &
                                    label_CheckImport,                                      &
                                    label_SetRunClock,                                      &
                                    label_TimestampExport,                                  &
                                    label_Finalize,                                         &
                                    NUOPC_ModelGet
  use module_mpas_config,     only: output_fh, dt_atmos, calendar, fcst_mpi_comm
  use module_mpas_config,     only: pio_ioformat, pio_iotype, pio_subsystem, pio_stride, pio_numiotasks, pio_iodesc
  use module_fcst_grid_comp,  only: fcstSS => SetServices

  implicit none
  private

  integer :: iau_offset = 0

  public SetServices

  type(ESMF_GridComp)               :: fcstComp
  type(ESMF_State)                  :: fcstState
  integer,dimension(:), allocatable :: fcstPetList
  integer, save                     :: FBCount
  logical                           :: profile_memory = .true.
  logical                           :: write_runtimelog = .false.
  logical                           :: lprint = .false.
  integer                           :: mype = 0
  integer                           :: dbug = 0
  real(kind=8)                      :: timere, timep2re
contains

  ! #########################################################################################
  ! ESMF entrypoints.
  ! #########################################################################################
  subroutine SetServices(gcomp, rc)

    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    character(len=*),parameter  :: subname='(mpasatm_cap:SetServices)'

    ! Initialize ESMF error message.
    rc = ESMF_SUCCESS

    ! The NUOPC model component will register the generic methods
    call NUOPC_CompDerive(gcomp, model_routine_SS, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! ESMF Initialize and Advertise
    call NUOPC_CompSpecialize(gcomp, specLabel=label_Advertise, specRoutine=InitializeAdvertise, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! ESMF Realize
    call NUOPC_CompSpecialize(gcomp, specLabel=label_RealizeProvided, specRoutine=InitializeRealize, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! ESMF Advance
    call NUOPC_CompSpecialize(gcomp, specLabel=label_Advance, specRoutine=ModelAdvance, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! Setup ESMF Run/Advance phase: phase1
    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_RUN, phaseLabelList=(/"phase1"/), userRoutine=routine_Run, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! ESMF Advance Phase 1.
    call NUOPC_CompSpecialize(gcomp, specLabel=label_Advance, specPhaseLabel="phase1", specRoutine=ModelAdvance_phase1, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! ESMF Finalize
    call NUOPC_CompSpecialize(gcomp, specLabel=label_Finalize, specRoutine=ModelFinalize, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

  end subroutine SetServices

  ! #########################################################################################
  !
  ! #########################################################################################
  subroutine InitializeAdvertise(gcomp, rc)
    use pio,                    only: PIO_64BIT_OFFSET, PIO_64BIT_DATA
    use pio,                    only: PIO_IOTYPE_NETCDF, PIO_IOTYPE_PNETCDF
    use pio,                    only: PIO_IOTYPE_NETCDF4C, PIO_IOTYPE_NETCDF4P
    use pio,                    only: PIO_REARR_BOX, PIO_REARR_SUBSET
    use pio,                    only: PIO_REARR_COMM_P2P, PIO_REARR_COMM_COLL
    use pio,                    only: PIO_REARR_COMM_FC_2D_ENABLE, PIO_REARR_COMM_FC_2D_DISABLE
    use pio,                    only: PIO_REARR_COMM_FC_1D_COMP2IO, PIO_REARR_COMM_FC_1D_IO2COMP
    use pio,                    only: PIO_DOUBLE, PIO_REAL
    use pio,                    only: pio_init, pio_setdebuglevel, pio_set_rearr_opts, pio_openfile, pio_nowrite, pio_iotask_rank
    type(ESMF_GridComp)                    :: gcomp
    integer, intent(out)                   :: rc

    ! local variables
    type(ESMF_State)                       :: importState, exportState
    type(ESMF_Clock)                       :: clock
    character(len=10)                      :: value
    character(240)                         :: msgString
    logical                                :: isPresent, isSet
    type(ESMF_VM)                          :: vm
    type(ESMF_TimeInterval)                :: timeStep
    type(ESMF_Config)                      :: cf
    integer                                :: i, j, k, urc, petcount
    real                                   :: nfhmax
    character(ESMF_MAXSTR)                 :: gc_name
    type(ESMF_Info)                        :: parentInfo, childInfo, info
    character(len=*),parameter             :: subname='(mpas_nuopc_cap:InitializeAdvertise)'
    real(kind=8)                           :: MPI_Wtime, timeis, timerhs
    integer                                :: num_threads
    integer                                :: num_pes_fcst
    ! PIO
    character(len=20)       :: cvalue
    integer                 :: pio_root
    integer                 :: pio_rearranger
    integer                 :: pio_debug_level
    integer                 :: pio_rearr_comm_type
    integer                 :: pio_rearr_comm_fcd
    logical                 :: pio_rearr_comm_enable_hs_comp2io
    logical                 :: pio_rearr_comm_enable_isend_comp2io
    integer                 :: pio_rearr_comm_max_pend_req_comp2io
    logical                 :: pio_rearr_comm_enable_hs_io2comp
    logical                 :: pio_rearr_comm_enable_isend_io2comp
    integer                 :: pio_rearr_comm_max_pend_req_io2comp
    integer                 :: ret

    ! Initialize ESMF error message.
    rc = ESMF_SUCCESS
    
    ! Timing info (debug mode)
    timeis = MPI_Wtime()

    ! ########################################################################################
    !
    ! Setup: Get information about the Component.
    !
    ! ########################################################################################
    ! Get this Component's VM 
    call ESMF_GridCompGet(gcomp, name=gc_name, vm=vm,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! Get internal information for this Component's VM.
    !
    ! - petcount is the number of PETs
    ! - mype is the master PET (000).
    !
    !call ESMF_VMGet(vm, petCount=petcount, localpet=mype, mpiCommunicator=fcst_mpi_comm%mpi_val, rc=rc)
    call ESMF_VMGet(vm=vm, localPet=mype, mpiCommunicator=fcst_mpi_comm%mpi_val, &
                    petCount=petcount, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! Get this Component's ESMF_infor object handle (info) and query the object handle to get
    ! the number of thread (num_threads).
    !
    ! This is needed for ESMF write-grid component, where num_threads is needed to compute
    ! actual wrttasks_per_group_from_parent
    call ESMF_InfoGetFromHost(gcomp, info=info, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc,  msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_InfoGet(info, key="/NUOPC/Hint/PePerPet/MaxCount", value=num_threads, default=1, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc,  msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! Query for importState and exportState
    call NUOPC_ModelGet(gcomp, driverClock=clock, importState=importState, exportState=exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! Setup memory profiler.
    call ESMF_AttributeGet(gcomp, name="ProfileMemory", value=value, defaultValue="false", &
                           convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    profile_memory = (trim(value)/="false")

    ! Setup runtime logger.
    call ESMF_AttributeGet(gcomp, name="RunTimeLog", value=value, defaultValue="false", &
                           convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    write_runtimelog = (trim(value)=="true")

    ! Read in cap debug flag.
    call NUOPC_CompAttributeGet(gcomp, name='dbug_flag', value=value, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    if (isPresent .and. isSet) then
     read(value,*) dbug
    end if
    write(msgString,'(A,i6)') trim(subname)//' dbug = ',dbug
    call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO)

    ! #######################################################################################
    !
    ! PIO
    !
    ! #######################################################################################
    ! pio_netcdf_format
    call NUOPC_CompAttributeGet(gcomp, name='pio_netcdf_format', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    if (isPresent .and. isSet) then
       cvalue = ESMF_UtilStringUpperCase(cvalue)
       if (trim(cvalue) .eq. 'CLASSIC') then
          pio_ioformat = 0
       else if (trim(cvalue) .eq. '64BIT_OFFSET') then
          pio_ioformat = PIO_64BIT_OFFSET
       else if (trim(cvalue) .eq. '64BIT_DATA') then
          pio_ioformat = PIO_64BIT_DATA
       else
          call ESMF_LogWrite(trim("need to provide valid option for pio_ioformat (CLASSIC|64BIT_OFFSET|64BIT_DATA)"), ESMF_LOGMSG_INFO)
          return
       end if
    else
       cvalue = '64BIT_OFFSET'
       pio_ioformat = PIO_64BIT_OFFSET
    end if

    ! pio_typename
    call NUOPC_CompAttributeGet(gcomp, name='pio_typename', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    if (isPresent .and. isSet) then
       cvalue = ESMF_UtilStringUpperCase(cvalue)
       if (trim(cvalue) .eq. 'NETCDF') then
          pio_iotype = PIO_IOTYPE_NETCDF
       else if (trim(cvalue) .eq. 'PNETCDF') then
          pio_iotype = PIO_IOTYPE_PNETCDF
       else if (trim(cvalue) .eq. 'NETCDF4C') then
          pio_iotype = PIO_IOTYPE_NETCDF4C
       else if (trim(cvalue) .eq. 'NETCDF4P') then
          pio_iotype = PIO_IOTYPE_NETCDF4P
       else
          call ESMF_LogWrite(trim("need to provide valid option for pio_typename (NETCDF|PNETCDF|NETCDF4C|NETCDF4P)"), ESMF_LOGMSG_INFO)
          return
       end if
    else
       cvalue = 'NETCDF'
       pio_iotype = PIO_IOTYPE_NETCDF
    end if

    ! pio_root
    call NUOPC_CompAttributeGet(gcomp, name='pio_root', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    if (isPresent .and. isSet) then
       read(cvalue,*) pio_root
       if (pio_root < 0) then
          pio_root = 1
       endif
       pio_root = min(pio_root, petCount-1)
    else
       pio_root = 1
    end if

    ! pio_stride
    call NUOPC_CompAttributeGet(gcomp, name='pio_stride', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    if (isPresent .and. isSet) then
       read(cvalue,*) pio_stride
    else
       pio_stride = -99
    end if

    ! pio_numiotasks
    call NUOPC_CompAttributeGet(gcomp, name='pio_numiotasks', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    if (isPresent .and. isSet) then
       read(cvalue,*) pio_numiotasks
    else
       pio_numiotasks = -99
    end if

    ! check for parallel IO, it requires at least two io pes
    if (petCount > 1 .and. pio_numiotasks == 1 .and. &
       (pio_iotype .eq. PIO_IOTYPE_PNETCDF .or. pio_iotype .eq. PIO_IOTYPE_NETCDF4P)) then
       pio_numiotasks = 2
       pio_stride = min(pio_stride, petCount/2)
    endif

    if (pio_root + (pio_stride)*(pio_numiotasks-1) >= petCount .or. &
        pio_stride <= 0 .or. pio_numiotasks <= 0 .or. pio_root < 0 .or. pio_root > petCount-1) then
       if (petCount < 100) then
          pio_stride = max(1, petCount/4)
       else if(petCount < 1000) then
          pio_stride = max(1, petCount/8)
       else
          pio_stride = max(1, petCount/16)
       end if
       if(pio_stride > 1) then
          pio_numiotasks = petCount/pio_stride
          pio_root = min(1, petCount-1)
       else
          pio_numiotasks = petCount
          pio_root = 0
       end if
    end if

    ! pio_rearranger
    call NUOPC_CompAttributeGet(gcomp, name='pio_rearranger', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    if (isPresent .and. isSet) then
       cvalue = ESMF_UtilStringUpperCase(cvalue)
       if (trim(cvalue) .eq. 'BOX') then
          pio_rearranger = PIO_REARR_BOX
       else if (trim(cvalue) .eq. 'SUBSET') then
          pio_rearranger = PIO_REARR_SUBSET
       else
          call ESMF_LogWrite(trim("need to provide valid option for pio_rearranger (BOX|SUBSET)"), ESMF_LOGMSG_INFO)
          return
       end if
    else
       cvalue = 'SUBSET'
       pio_rearranger = PIO_REARR_SUBSET
    end if

    ! Initialize PIO
    allocate(pio_subsystem)
    call pio_init(mype, fcst_mpi_comm%mpi_val, pio_numiotasks, 0, pio_stride, pio_rearranger, pio_subsystem, base=pio_root)
    
    ! PIO debug related options
    ! pio_debug_level
    call NUOPC_CompAttributeGet(gcomp, name='pio_debug_level', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    if (isPresent .and. isSet) then
       read(cvalue,*) pio_debug_level
       if (pio_debug_level < 0 .or. pio_debug_level > 6) then
          call ESMF_LogWrite(trim("MPAS_NUOPC_CAP: need to provide valid option for pio_debug_level (0-6)"), ESMF_LOGMSG_INFO)
          return
       end if
    else
       pio_debug_level = 0
    end if

    ! set PIO debug level
    call pio_setdebuglevel(pio_debug_level)

    ! query shared PIO rearranger attributes
    ! pio_rearr_comm_type
    call NUOPC_CompAttributeGet(gcomp, name='pio_rearr_comm_type', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    if (isPresent .and. isSet) then
       cvalue = ESMF_UtilStringUpperCase(cvalue)
       if (trim(cvalue) .eq. 'P2P') then
          pio_rearr_comm_type = PIO_REARR_COMM_P2P
       else if (trim(cvalue) .eq. 'COLL') then
          pio_rearr_comm_type = PIO_REARR_COMM_COLL
       else
          call ESMF_LogWrite("MPAS_NUOPC_CAP: need to provide valid option for pio_rearr_comm_type (P2P|COLL)",ESMF_LOGMSG_INFO)
          return
       end if
    else
       cvalue = 'P2P'
       pio_rearr_comm_type = PIO_REARR_COMM_P2P
    end if

    ! pio_rearr_comm_fcd
    call NUOPC_CompAttributeGet(gcomp, name='pio_rearr_comm_fcd', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    if (isPresent .and. isSet) then
       cvalue = ESMF_UtilStringUpperCase(cvalue)
       if (trim(cvalue) .eq. '2DENABLE') then
          pio_rearr_comm_fcd = PIO_REARR_COMM_FC_2D_ENABLE
       else if (trim(cvalue) .eq. 'IO2COMP') then
          pio_rearr_comm_fcd = PIO_REARR_COMM_FC_1D_IO2COMP
       else if (trim(cvalue) .eq. 'COMP2IO') then
          pio_rearr_comm_fcd = PIO_REARR_COMM_FC_1D_COMP2IO
       else if (trim(cvalue) .eq. '2DDISABLE') then
          pio_rearr_comm_fcd = PIO_REARR_COMM_FC_2D_DISABLE
       else
          call ESMF_LogWrite("MPAS_NUOPC_CAP: need to provide valid option for pio_rearr_comm_fcd (2DENABLE|IO2COMP|COMP2IO|2DDISABLE)",ESMF_LOGMSG_INFO)
          return
       end if
    else
       cvalue = '2DENABLE'
       pio_rearr_comm_fcd = PIO_REARR_COMM_FC_2D_ENABLE
    end if

    ! pio_rearr_comm_enable_hs_comp2io
    call NUOPC_CompAttributeGet(gcomp, name='pio_rearr_comm_enable_hs_comp2io', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    if (isPresent .and. isSet) then
       read(cvalue,*) pio_rearr_comm_enable_hs_comp2io
    else
       pio_rearr_comm_enable_hs_comp2io = .true.
    end if

    ! pio_rearr_comm_enable_isend_comp2io
    call NUOPC_CompAttributeGet(gcomp, name='pio_rearr_comm_enable_isend_comp2io', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    if (isPresent .and. isSet) then
       read(cvalue,*) pio_rearr_comm_enable_isend_comp2io
    else
       pio_rearr_comm_enable_isend_comp2io = .false.
    end if

    ! pio_rearr_comm_max_pend_req_comp2io
    call NUOPC_CompAttributeGet(gcomp, name='pio_rearr_comm_max_pend_req_comp2io', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    if (isPresent .and. isSet) then
       read(cvalue,*) pio_rearr_comm_max_pend_req_comp2io
    else
       pio_rearr_comm_max_pend_req_comp2io = 0
    end if

    ! pio_rearr_comm_enable_hs_io2comp
    call NUOPC_CompAttributeGet(gcomp, name='pio_rearr_comm_enable_hs_io2comp', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    if (isPresent .and. isSet) then
       read(cvalue,*) pio_rearr_comm_enable_hs_io2comp
    else
       pio_rearr_comm_enable_hs_io2comp = .false.
    end if

    ! pio_rearr_comm_enable_isend_io2comp
    call NUOPC_CompAttributeGet(gcomp, name='pio_rearr_comm_enable_isend_io2comp', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    if (isPresent .and. isSet) then
       read(cvalue,*) pio_rearr_comm_enable_isend_io2comp
    else
       pio_rearr_comm_enable_isend_io2comp = .true.
    end if

    ! pio_rearr_comm_max_pend_req_io2comp
    call NUOPC_CompAttributeGet(gcomp, name='pio_rearr_comm_max_pend_req_io2comp', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    if (isPresent .and. isSet) then
       read(cvalue,*) pio_rearr_comm_max_pend_req_io2comp
    else
       pio_rearr_comm_max_pend_req_io2comp = 64
    end if
    
    ! #######################################################################################
    !
    ! Get configuration variables.
    !
    ! #######################################################################################
    CF = ESMF_ConfigCreate(rc=rc)
    call ESMF_ConfigLoadFile(config=CF ,filename='model_configure' ,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! Calendar
    call ESMF_ConfigGetAttribute(config=CF,value=calendar, label ='calendar:', default='gregorian',rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! IAU offset
    call ESMF_ConfigGetAttribute(config=CF,value=iau_offset,default=0,label ='iau_offset:',rc=rc)
    if (iau_offset < 0) iau_offset=0

    ! Timestep (dt_atmos) and forecast length (nfhmax)
    call ESMF_ConfigGetAttribute(config=CF, value=dt_atmos, label ='dt_atmos:',   rc=rc)
    call ESMF_ConfigGetAttribute(config=CF, value=nfhmax,   label ='nhours_fcst:',rc=rc)
    if(mype == 0) print *,'MPAS_NUOPC_CAP: dt_atmos = ',dt_atmos,', nfhmax = ',nfhmax

    ! Set ESMF time interval.
    call ESMF_TimeIntervalSet(timeStep, s=dt_atmos, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! Print messages on master processor?
    if( mype == 0) lprint = .true.

    ! #######################################################################################
    !
    ! Initialize fcst grid component
    !
    ! #######################################################################################
    ! Create ESMF grid component for each MPI process.
    num_pes_fcst = petcount
    allocate(fcstPetList(num_pes_fcst))
    do j=1, num_pes_fcst
      fcstPetList(j) = j - 1
    enddo
    fcstComp = ESMF_GridCompCreate(petList=fcstPetList, name='mpas_fcst', rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! Attributes from forecast Component
    call ESMF_InfoGetFromHost(gcomp, info=parentInfo, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc,  msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! Attributes from MPAS component
    call ESMF_InfoGetFromHost(fcstComp, info=childInfo, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc,  msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! Copy attributes from MPAS component to forecast Component
    call ESMF_InfoUpdate(lhs=childInfo, rhs=parentInfo, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc,  msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! Use the generic SetVM method to do resource and threading control
    call ESMF_GridCompSetVM(fcstComp, SetVM, userRc=urc, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc,  msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

    ! Use MPAS forecast Component routines (fcstSS) for Initialize(), Run(), and Finalize() services.
    call ESMF_GridCompSetServices(fcstComp, fcstSS, userRc=urc, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc,  msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

    ! Create fcst state
    fcstState = ESMF_StateCreate(rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! ########################################################################################
    !
    ! Call fcst_initialize (including creating fcstgrid and fcst fieldbundle)
    !
    ! ########################################################################################
    call ESMF_GridCompInitialize(fcstComp, exportState=fcstState,  clock=clock, phase=1, userRc=urc, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc,  msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

    ! Reconcile the fcstComp's export state
    call ESMF_StateReconcile(fcstState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! Determine number elements in fcstState
    call ESMF_StateGet(fcstState, itemCount=FBCount, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    if(mype == 0) print *,'MPAS_NUOPC_CAP: field bundles in fcstComp export state, FBCount= ',FBcount

    ! ########################################################################################
    !
    ! Call fcst_advertise
    !
    ! ########################################################################################
    call ESMF_GridCompInitialize(fcstComp, importState=importState, exportState=exportState, &
                                 clock=clock, phase=2, userRc=urc, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc,  msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

    ! Timing info (debug mode)
    if (mype == 0) print *,'MPAS_NUOPC_CAP, initAdvertise time=',MPI_Wtime()-timeis,mype

  end subroutine InitializeAdvertise
  
  ! ########################################################################################
  !
  ! ########################################################################################
  subroutine InitializeRealize(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! Locals
    integer :: urc
    real(8) :: mpi_wtime, timeirs
    type(ESMF_Clock) :: clock
    type(ESMF_State) :: importState, exportState

    ! Initialize ESMF error message.
    rc = ESMF_SUCCESS

    ! Timing info (debug mode)
    timeirs = MPI_Wtime()

    ! Query for importState and exportState
    call NUOPC_ModelGet(gcomp, driverClock=clock, importState=importState, exportState=exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! ########################################################################################
    !
    ! Call fcst_Realize
    !
    ! ########################################################################################
    call ESMF_GridCompInitialize(fcstComp, importState=importState, exportState=exportState, &
                                 clock=clock, phase=3, userRc=urc, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc,  msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

    ! Timing info (debug mode)
    if (mype == 0) print *,'MPAS_NUOPC_CAP, InitializeRealize time=', MPI_Wtime()-timeirs, mype

  end subroutine InitializeRealize

  ! ########################################################################################
  !
  ! This routine calls ALL run phase(s) of the ESMF forecast grid component.
  !
  ! ########################################################################################
  subroutine ModelAdvance(gcomp, rc)
    type(ESMF_GridComp)   :: gcomp
    integer, intent(out)  :: rc

    ! Locals
    character(len=*),parameter  :: subname='(mpas_nuopc_cap:ModelAdvance)'
    real(kind=8)                :: MPI_Wtime, timers

    ! Initialize ESMF error message.
    rc = ESMF_SUCCESS

    ! Timing info (debug mode)
    timers = MPI_Wtime()
    if(write_runtimelog .and. timere>0. .and. lprint) print *,'MPAS_NUOPC_CAP, time between ModelAdvance phase=', timers-timere, mype

    ! Begin memory profiling.
    if (profile_memory) call ESMF_VMLogMemInfo("Entering MPAS ModelAdvance: ")

    ! Call Run phases...
    call ModelAdvance_phase1(gcomp, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! Finish memory profiling.
    if (profile_memory) call ESMF_VMLogMemInfo("Leaving MPAS ModelAdvance: ")

    ! Timing info (debug mode)
    timere = MPI_Wtime()
    if (mype == 0) print *,'MPAS_NUOPC_CAP, ModelAdvance time=', timere-timers, mype

  end subroutine ModelAdvance
  
  ! ########################################################################################
  !
  ! This routine calls the first run phase of the ESMF forecast grid component.
  !
  ! ########################################################################################
  subroutine ModelAdvance_phase1(gcomp, rc)
    type(ESMF_GridComp)         :: gcomp
    integer, intent(out)        :: rc

    ! Locals
    character(len=*),parameter  :: subname='(mpas_nuopc_cap:ModelAdvance_phase1)'
    real(kind=8)                :: MPI_Wtime, timep1rs, timep1re
    type(ESMF_Clock)            :: clock
    integer                     :: urc
    character(240)              :: msgString

    ! Initialize ESMF error message.
    rc = ESMF_SUCCESS

    ! Timing info (debug mode)
    timep1rs = MPI_Wtime()
    if(write_runtimelog .and. timep2re>0. .and. lprint) print *,'MPAS_NUOPC_CAP, time between mpas run phase2 and phase1 ', timep1rs-timep2re,mype

    ! Begin memory profiling. 
    if(profile_memory) call ESMF_VMLogMemInfo("Entering MPAS ModelAdvance_phase1: ")

    ! Get information on grid component.
    call ESMF_GridCompGet(gcomp, clock=clock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! Display clock information.
    call ESMF_ClockPrint(clock, options="currTime", &
                         preString="entering MPAS_ADVANCE phase1 with clock current: ", &
                         unit=msgString)
    call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO)
    call ESMF_ClockPrint(clock, options="startTime", &
                         preString="entering MPAS_ADVANCE phase1 with clock start:   ", &
                         unit=msgString)
    call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO)
    call ESMF_ClockPrint(clock, options="stopTime", &
                         preString="entering MPAS_ADVANCE phase1 with clock stop:    ", &
                         unit=msgString)
    call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO)

    ! Call Run phase 1...
    call ESMF_GridCompRun(fcstComp, exportState=fcstState, clock=clock, phase=1, userRc=urc, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

    
    ! Finish memory profiling.
    if (profile_memory) call ESMF_VMLogMemInfo("Leaving MPAS ModelAdvance_phase1: ")

    ! Timing info (debug mode)
    timep1re = MPI_Wtime()
    if(write_runtimelog .and. lprint) print *,'MPAS_NUOPC_CAP, ModelAdvance phase1 time=', timep1re-timep1rs,mype

  end subroutine ModelAdvance_phase1

  ! ########################################################################################
  !
  ! ########################################################################################
  subroutine ModelFinalize(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! Locals
    character(len=*),parameter :: subname='(mpasatm_cap:ModelFinalize)'
    type(ESMF_VM)              :: vm
    real(kind=8)               :: MPI_Wtime, timeffs
    integer                    :: urc

    ! Initialize ESMF error message.
    rc = ESMF_SUCCESS

    ! Timing info (debug mode) 
    timeffs = MPI_Wtime()

    ! Get information on grid component. 
    call ESMF_GridCompGet(gcomp,vm=vm,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! Call fcst_finalize
    call ESMF_GridCompFinalize(fcstComp, exportState=fcststate,userRc=urc, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc,  msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

    ! Destroy forecast grid components
    call ESMF_StateDestroy(fcstState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_GridCompDestroy(fcstComp, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! Timing info (debug mode) 
    if(write_runtimelog .and. lprint) print *,'MPAS_NUOPC_CAP, ModelFinalize time=',MPI_Wtime()-timeffs,mype

  end subroutine ModelFinalize

end module mpasatm_cap_mod
