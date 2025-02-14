! ###########################################################################################
!> \file mpas_cap.F90
!>
!> This file contains the NUOPC Cap for the UWM Atmosphere with the MPAS dynamical core.
!>
! ###########################################################################################
module mpasatm_cap_mod

  use ESMF
  use NUOPC
  use NUOPC_Model,            only: model_routine_SS => SetServices,         &
                                    SetVM,                                   &
                                    routine_Run,                             &
                                    label_Advertise,                         &
                                    label_RealizeProvided,                   &
                                    label_Advance,                           &
                                    label_CheckImport,                       &
                                    label_SetRunClock,                       &
                                    label_TimestampExport,                   &
                                    label_Finalize,                          &
                                    NUOPC_ModelGet

  use module_mpas_config,     only: quilting, quilting_restart, output_fh,   &
                                    dt_atmos,                                &
                                    calendar, cpl_grid_id,                   &
                                    cplprint_flag, first_kdt

  use module_fcst_grid_comp,  only: fcstSS => SetServices

  use module_cplscalars,      only: flds_scalar_name, flds_scalar_num,          &
                                    flds_scalar_index_nx, flds_scalar_index_ny, &
                                    flds_scalar_index_ntile

  implicit none
  private

  integer :: iau_offset = 0

  public SetServices

  type(ESMF_GridComp)                         :: fcstComp
  type(ESMF_State)                            :: fcstState
  integer,dimension(:), allocatable           :: fcstPetList
  integer, save                               :: FBCount
  logical                                     :: profile_memory = .true.
  logical                                     :: write_runtimelog = .false.
  logical                                     :: lprint = .false.
  integer                                     :: mype = -1
  integer                                     :: dbug = 0

contains

  ! #########################################################################################
  !
  ! #########################################################################################
  subroutine SetServices(gcomp, rc)

    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    character(len=*),parameter  :: subname='(mpasatm_cap:SetServices)'

    rc = ESMF_SUCCESS

    ! the NUOPC model component will register the generic methods
    call NUOPC_CompDerive(gcomp, model_routine_SS, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! set entry point for methods that require specific implementation
    call NUOPC_CompSpecialize(gcomp, specLabel=label_Advertise, specRoutine=InitializeAdvertise, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

  end subroutine SetServices

  ! #########################################################################################
  !
  ! #########################################################################################
  subroutine InitializeAdvertise(gcomp, rc)
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
    character(len=*),parameter             :: subname='(mpas_cap:InitializeAdvertise)'
    real(kind=8)                           :: MPI_Wtime, timeis, timerhs
    integer                                :: num_threads
    character(len=20)                      :: cvalue
    integer                                :: num_pes_fcst

    rc = ESMF_SUCCESS
    timeis = MPI_Wtime()

    call ESMF_GridCompGet(gcomp, name=gc_name, vm=vm,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_VMGet(vm, petCount=petcount, localpet=mype, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! num_threads is needed to compute actual wrttasks_per_group_from_parent
    call ESMF_InfoGetFromHost(gcomp, info=info, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc,  msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_InfoGet(info, key="/NUOPC/Hint/PePerPet/MaxCount", value=num_threads, default=1, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc,  msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! query for importState and exportState
    call NUOPC_ModelGet(gcomp, driverClock=clock, importState=importState, exportState=exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_AttributeGet(gcomp, name="cpl_grid_id", value=value, defaultValue="1", &
                           convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    cpl_grid_id = ESMF_UtilString2Int(value, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_AttributeGet(gcomp, name="ProfileMemory", value=value, defaultValue="false", &
                           convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    profile_memory = (trim(value)/="false")

    call ESMF_AttributeGet(gcomp, name="RunTimeLog", value=value, defaultValue="false", &
                           convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    write_runtimelog = (trim(value)=="true")

    call ESMF_AttributeGet(gcomp, name="DumpFields", value=value, defaultValue="false", &
                           convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    cplprint_flag = (trim(value)=="true")
    write(msgString,'(A,l6)') trim(subname)//' cplprint_flag = ',cplprint_flag
    call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO)

    ! Read in cap debug flag
    call NUOPC_CompAttributeGet(gcomp, name='dbug_flag', value=value, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    if (isPresent .and. isSet) then
     read(value,*) dbug
    end if
    write(msgString,'(A,i6)') trim(subname)//' dbug = ',dbug
    call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO)

    ! set cpl_scalars from config. Default to null values for standalone
    flds_scalar_name = ''
    flds_scalar_num = 0
    flds_scalar_index_nx = 0
    flds_scalar_index_ny = 0
    flds_scalar_index_ntile = 0
    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldName", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    if (isPresent .and. isSet) then
       flds_scalar_name = trim(cvalue)
       call ESMF_LogWrite(trim(subname)//' flds_scalar_name = '//trim(flds_scalar_name), ESMF_LOGMSG_INFO)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    endif
    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldCount", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    if (isPresent .and. isSet) then
       read(cvalue, *) flds_scalar_num
       write(msgString,*) flds_scalar_num
       call ESMF_LogWrite(trim(subname)//' flds_scalar_num = '//trim(msgString), ESMF_LOGMSG_INFO)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    endif
    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldIdxGridNX", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) flds_scalar_index_nx
       write(msgString,*) flds_scalar_index_nx
       call ESMF_LogWrite(trim(subname)//' : flds_scalar_index_nx = '//trim(msgString), ESMF_LOGMSG_INFO)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    endif
    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldIdxGridNY", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) flds_scalar_index_ny
       write(msgString,*) flds_scalar_index_ny
       call ESMF_LogWrite(trim(subname)//' : flds_scalar_index_ny = '//trim(msgString), ESMF_LOGMSG_INFO)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    endif
    ! tile index must be present if indices for nx and ny are non-zero
    if (flds_scalar_index_nx /= 0 .and. flds_scalar_index_ny /=0 ) then
       call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldIdxGridNTile", isPresent=isPresent, isSet=isSet, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
       if (.not. isPresent .and. .not. isSet) then
          if (mype == 0)write(*,*)'ERROR : ScalarFieldIdxGridNTile must be set'
          call ESMF_LogWrite('ERROR : ScalarFieldIdxGridNTile must be set', ESMF_LOGMSG_ERROR)
          rc = ESMF_FAILURE
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
       else
          call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldIdxGridNTile", value=cvalue, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
          read(cvalue,*) flds_scalar_index_ntile
          write(msgString,*) flds_scalar_index_ntile
          call ESMF_LogWrite(trim(subname)//' : flds_scalar_index_ntile = '//trim(msgString), ESMF_LOGMSG_INFO)
       endif
    end if

    ! #######################################################################################
    ! Get configuration variables.
    ! #######################################################################################
    CF = ESMF_ConfigCreate(rc=rc)
    call ESMF_ConfigLoadFile(config=CF ,filename='model_configure' ,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_ConfigGetAttribute(config=CF,value=calendar, &
                                 label ='calendar:', &
                                 default='gregorian',rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return


    call ESMF_ConfigGetAttribute(config=CF,value=iau_offset,default=0,label ='iau_offset:',rc=rc)
    if (iau_offset < 0) iau_offset=0

    call ESMF_ConfigGetAttribute(config=CF, value=dt_atmos, label ='dt_atmos:',   rc=rc)
    call ESMF_ConfigGetAttribute(config=CF, value=nfhmax,   label ='nhours_fcst:',rc=rc)
    if(mype == 0) print *,'af ufs config,dt_atmos=',dt_atmos,'nfhmax=',nfhmax

    call ESMF_TimeIntervalSet(timeStep, s=dt_atmos, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    first_kdt = 1
    if( mype == 0) lprint = .true.

    ! #######################################################################################
    ! Initialize fcst grid component
    ! #######################################################################################

    num_pes_fcst = petcount
    allocate(fcstPetList(num_pes_fcst))
    do j=1, num_pes_fcst
      fcstPetList(j) = j - 1
    enddo
    fcstComp = ESMF_GridCompCreate(petList=fcstPetList, name='mpas_fcst', rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! Copy attributes from mpascap component to fcstComp
    call ESMF_InfoGetFromHost(gcomp, info=parentInfo, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc,  msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_InfoGetFromHost(fcstComp, info=childInfo, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc,  msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_InfoUpdate(lhs=childInfo, rhs=parentInfo, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc,  msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! Use the generic SetVM method to do resource and threading control
    call ESMF_GridCompSetVM(fcstComp, SetVM, userRc=urc, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc,  msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
    call ESMF_GridCompSetServices(fcstComp, fcstSS, userRc=urc, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc,  msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

    ! Create fcst state
    fcstState = ESMF_StateCreate(rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! Call fcst Initialize (including creating fcstgrid and fcst fieldbundle)
    call ESMF_GridCompInitialize(fcstComp, exportState=fcstState,    &
                                 clock=clock, phase=1, userRc=urc, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

    ! reconcile the fcstComp's export state
    call ESMF_StateReconcile(fcstState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! determine number elements in fcstState
    call ESMF_StateGet(fcstState, itemCount=FBCount, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    if(mype == 0) print *,'mpas_cap: field bundles in fcstComp export state, FBCount= ',FBcount

    ! #######################################################################################
    ! Advertise Fields in importState and exportState.
    ! #######################################################################################
!    call ESMF_GridCompInitialize(fcstComp, importState=importState, exportState=exportState, &
!                                 clock=clock, phase=2, userRc=urc, rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
!
!    if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
!
!    if(write_runtimelog .and. lprint) print *,'in mpas_cap, init time=',MPI_Wtime()-timeis,mype

  end subroutine InitializeAdvertise

end module mpasatm_cap_mod
