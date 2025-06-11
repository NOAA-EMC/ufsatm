module module_cap_cpl
!
!*** this module contains the debug subroutines for fv3 coupled run
!
! revision history
!  12 Mar 2018: J. Wang       Pull coupled subroutines from fv3_cap.F90 to this module
!
  use ESMF

  implicit none

  private
  public diagnose_cplFields
  contains

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------

    subroutine diagnose_cplFields(gcomp, clock_fv3, fcstpe, &
                                  statewrite_flag, stdiagnose_flag, state_tag, rc)

      type(ESMF_GridComp), intent(in)       :: gcomp
      type(ESMF_Clock),intent(in)           :: clock_fv3
      logical, intent(in)                   :: fcstpe
      logical, intent(in)                   :: statewrite_flag
      integer, intent(in)                   :: stdiagnose_flag
      character(len=*), intent(in)          :: state_tag                        !< "import" or "export".
      integer, intent(out)                  :: rc

      character(len=*),parameter :: subname='(module_cap_cpl:diagnose_cplFields)'
      type(ESMF_Time) :: currTime
      type(ESMF_State) :: state
      type(ESMF_TimeInterval) :: timeStep
      character(len=240) :: import_timestr, export_timestr
      character(len=160) :: nuopcMsg
      character(len=160) :: filename
!
      call ESMF_ClockGet(clock_fv3, currTime=currTime, timeStep=timestep, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      call ESMF_TimeGet(currTime, timestring=import_timestr, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      call ESMF_TimeGet(currTime+timestep, timestring=export_timestr, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call ESMF_ClockPrint(clock_fv3, options="currTime", preString="current time: ", unit=nuopcMsg)
      call ESMF_LogWrite(trim(subname)//' '//trim(state_tag)//' '//trim(nuopcMsg), ESMF_LOGMSG_INFO)

      if(trim(state_tag) .eq. 'import')then
        call ESMF_GridCompGet(gcomp, importState=state, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        if(stdiagnose_flag > 0 .and. fcstpe)then
         call state_diagnose(state, ':IS', rc=rc)
        end if

        ! Dump Fields out
        if (statewrite_flag) then
          write(filename,'(A)') 'fv3_cap_import_'//trim(import_timestr)//'.tile*.nc'
          call State_RWFields_tiles(state,trim(filename), rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
        end if
      end if

      if(trim(state_tag) .eq. 'export')then
        call ESMF_GridCompGet(gcomp, exportState=state, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        if(stdiagnose_flag > 0 .and. fcstpe)then
         call state_diagnose(state, ':ES', rc=rc)
        end if

        ! Dump Fields out
        if (statewrite_flag) then
          write(filename,'(A)') 'fv3_cap_export_'//trim(export_timestr)//'.tile*.nc'
          call State_RWFields_tiles(state,trim(filename), rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
        end if
      end if

    end subroutine diagnose_cplFields

  !-----------------------------------------------------------------------------

  ! This subroutine requires ESMFv8 - for coupled FV3
    subroutine State_RWFields_tiles(state,filename,rc)

      type(ESMF_State), intent(in)          :: state
      character(len=*), intent(in)          :: fileName
      integer, intent(out)                  :: rc

      ! local variables
      type(ESMF_Array)                       :: array
      type(ESMF_Grid)                        :: grid
      type(ESMF_FieldBundle)                 :: fieldbundle
      type(ESMF_Field), allocatable          :: flds(:)
      type(ESMF_DistGrid)                    :: distgrid
      integer                                :: i, icount, ifld, id
      integer                                :: fieldcount, firstfld
      integer                                :: fieldDimCount, gridDimCount, dimCount, tileCount, ungriddedDimCount
      character(64), allocatable             :: itemNameList(:), fldNameList(:)
      type(ESMF_StateItem_Flag), allocatable :: typeList(:)
      integer, allocatable                   :: minIndexPTile(:,:), maxIndexPTile(:,:)
      integer, allocatable                   :: ungriddedLBound(:), ungriddedUBound(:)
      integer, allocatable                   :: fieldDimLen(:)
      character(len=32), allocatable         :: gridded_dim_labels(:), ungridded_dim_labels(:)

      character(16), parameter :: convention = 'NetCDF'
      character(16), parameter :: purpose = 'FV3'

      integer, parameter :: max_n_axes = 4
      integer, parameter :: max_n_dim = 16
      integer, dimension(max_n_axes, max_n_dim) :: axes_dimcount = 0

      character(len=*),parameter :: subname='(module_cap_cpl:State_RWFields_tiles)'

      rc = ESMF_SUCCESS
      !call ESMF_LogWrite(trim(subname)//trim(filename)//": called", ESMF_LOGMSG_INFO, rc=rc)

      call ESMF_StateGet(state, itemCount=icount, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      allocate(typeList(icount), itemNameList(icount))
      call ESMF_StateGet(state, itemTypeList=typeList, itemNameList=itemNameList, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      ! find first stateitem that is a field and the count of fields
      firstfld = 0; fieldcount = 0
      do i = icount,1,-1
        if(typeList(i) == ESMF_STATEITEM_FIELD) firstfld = i
        if(typeList(i) == ESMF_STATEITEM_FIELD) fieldcount = fieldcount + 1
      enddo

      allocate(flds(fieldCount),fldNameList(fieldCount))
      ifld = 1
      do i = 1, icount
        if(typeList(i) == ESMF_STATEITEM_FIELD) then
          fldNameList(ifld) = itemNameList(i)
          ifld = ifld + 1
        endif
      enddo

      fieldbundle = ESMF_FieldBundleCreate(rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call ESMF_LogWrite(trim(subname)//": write "//trim(filename), ESMF_LOGMSG_INFO, rc=rc)

      do ifld=1, fieldCount
        call ESMF_StateGet(state, itemName=fldNameList(ifld), field=flds(ifld), rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        call ESMF_FieldGet(flds(ifld), grid=grid, dimCount=fieldDimCount, array=array, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        if (fieldDimCount > 4) then
          call ESMF_LogWrite(trim(subname)//": fieldDimCount > 4 unsupported", ESMF_LOGMSG_ERROR, rc=rc)
        end if

        call ESMF_GridGet(grid, dimCount=gridDimCount, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        if (gridDimCount > 2) then
          call ESMF_LogWrite(trim(subname)//": gridDimCount > 2 unsupported", ESMF_LOGMSG_ERROR, rc=rc)
        end if

        call ESMF_ArrayGet(array, distgrid=distgrid, dimCount=dimCount, tileCount=tileCount, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        ! skip 'cpl_scalars' field because it has tileCount == 1, while all other fields have 6.
        ! This causes the following error:
        ! 20240705 134459.788 ERROR            PET000 ESMCI_IO.C:1614 ESMCI::IO::checkNtiles() Wrong data value  - New number of tiles (6) does not match previously-set number of tiles (1) for this IO object. All arrays handled by a given IO object must have the same number of tiles.
        if (trim(fldNameList(ifld)) == 'cpl_scalars') then
          cycle
        endif

        allocate(fieldDimLen(fieldDimCount))

        allocate(minIndexPTile(dimCount, tileCount))
        allocate(maxIndexPTile(dimCount, tileCount))
        call ESMF_DistGridGet(distgrid, minIndexPTile=minIndexPTile, maxIndexPTile=maxIndexPTile, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        allocate(gridded_dim_labels(gridDimCount))
        do i = 1, gridDimCount
          fieldDimLen(i) = maxIndexPTile(i,1) - minIndexPTile(i,1) + 1
          id = find_axis_id_for_axis_count(i,fieldDimLen(i))
          if (id < 1) then
            call ESMF_LogWrite(trim(subname)//": id < 1", ESMF_LOGMSG_ERROR, rc=rc)
          endif
          if (i == 1) write(gridded_dim_labels(i),'(A,I0)') 'xaxis_',id
          if (i == 2) write(gridded_dim_labels(i),'(A,I0)') 'yaxis_',id
        end do

        deallocate(minIndexPTile)
        deallocate(maxIndexPTile)

        ungriddedDimCount = fieldDimCount - gridDimCount
        allocate(ungridded_dim_labels(ungriddedDimCount))
        if (fieldDimCount > gridDimCount) then
          allocate(ungriddedLBound(ungriddedDimCount))
          allocate(ungriddedUBound(ungriddedDimCount))
          call ESMF_FieldGet(flds(ifld), ungriddedLBound=ungriddedLBound, ungriddedUBound=ungriddedUBound, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          do i=1,ungriddedDimCount
            fieldDimLen(i+gridDimCount) = ungriddedUBound(i) - ungriddedLBound(i) + 1
            id = find_axis_id_for_axis_count(i+gridDimCount, fieldDimLen(i+gridDimCount))
            if (id < 1) then
              write(0,*)'stop error', id, i, fieldDimLen(i+gridDimCount)
            endif
            if (i==1) write(ungridded_dim_labels(i),'(A,I0)') 'zaxis_',id
            if (i==2) write(ungridded_dim_labels(i),'(A,I0)') 'taxis_',id
          end do
          deallocate(ungriddedLBound)
          deallocate(ungriddedUBound)
        end if

        call ESMF_AttributeAdd(grid, convention=convention, purpose=purpose, attrList=(/ ESMF_ATT_GRIDDED_DIM_LABELS /), rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        call ESMF_AttributeSet(grid, convention=convention, purpose=purpose, &
                               name=ESMF_ATT_GRIDDED_DIM_LABELS, valueList=gridded_dim_labels, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        if (ungriddedDimCount > 0) then
          call ESMF_AttributeAdd(flds(ifld), convention=convention, purpose=purpose, &
                                 attrList=(/ ESMF_ATT_UNGRIDDED_DIM_LABELS /), rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          call ESMF_AttributeSet(flds(ifld), convention=convention, purpose=purpose, &
                                 name=ESMF_ATT_UNGRIDDED_DIM_LABELS, valueList=ungridded_dim_labels, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
        end if

        deallocate(fieldDimLen)
        deallocate(gridded_dim_labels)
        deallocate(ungridded_dim_labels)

        call ESMF_FieldBundleAdd(fieldbundle, (/flds(ifld)/), rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      enddo

      call ESMF_FieldBundleWrite(fieldbundle, fileName=trim(filename), convention=convention, purpose=purpose, &
           timeslice=1, overwrite=.true., rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

! -- Finalize
      deallocate(flds)
      deallocate(fldNameList)

      call ESMF_FieldBundleDestroy(fieldbundle, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      !call ESMF_LogWrite(trim(subname)//trim(filename)//": finished", ESMF_LOGMSG_INFO, rc=rc)

    contains

      function find_axis_id_for_axis_count(axis, count) result(id)
        integer, intent(in) :: axis, count

        integer :: id
        integer :: i

        id = -1 ! not found

        if (axis > max_n_axes) then
          call ESMF_LogWrite('axis > max_n_axes. Increase max_n_axes in '//trim(subname), ESMF_LOGMSG_ERROR)
          return
        end if

        do i =1, max_n_dim
           if (axes_dimcount(axis, i) == 0) then
             axes_dimcount(axis, i) = count
             id = i
             return
           else
             if (axes_dimcount(axis, i) == count) then
                 id = i
                 return
             end if
           end if
        end do

        call ESMF_LogWrite('Increase max_n_dim in '//trim(subname), ESMF_LOGMSG_ERROR)

      end function find_axis_id_for_axis_count

    end subroutine State_RWFields_tiles

  !-----------------------------------------------------------------------------

  subroutine state_diagnose(State,string, rc)
    ! ----------------------------------------------
    ! Diagnose status of state
    ! ----------------------------------------------
    type(ESMF_State), intent(inout) :: State
    character(len=*), intent(in), optional :: string
    integer, intent(out), optional  :: rc

    ! local variables
    integer                       :: i,j,n
    integer                       :: itemCount
    character(len=64) ,pointer    :: itemNameList(:)
    character(len=64)             :: lstring
    character(len=256)            :: tmpstr

    type(ESMF_Field)            :: lfield
    type(ESMF_StateItem_Flag)   :: itemType
    real(ESMF_KIND_R8), pointer :: dataPtr2d(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr3d(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr4d(:,:,:,:)
    integer                     :: lrc, localDeCount, dimCount
    character(len=*),parameter  :: subname='(FV3: state_diagnose)'

    lstring = ''
    if (present(string)) then
       lstring = trim(string)
    endif

    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    call ESMF_StateGet(State, itemCount=itemCount, rc=lrc)
    if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    allocate(itemNameList(itemCount))

    call ESMF_StateGet(State, itemNameList=itemNameList, rc=lrc)
    if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    do n = 1, itemCount
       call ESMF_StateGet(State, itemName=trim(itemNameList(n)), itemType=itemType, rc=lrc)
       if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

     if(itemType == ESMF_STATEITEM_FIELD)then
       call ESMF_StateGet(State, itemName=trim(itemNameList(n)), field=lfield, rc=lrc)
       if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

       call ESMF_FieldGet(lfield, localDeCount=localDeCount, dimCount=dimcount, rc=lrc)
       if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

       if(localDeCount.gt.0) then
         if(dimcount == 2)then
           call ESMF_FieldGet(lfield, farrayPtr=dataPtr2d, rc=lrc)
           if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

           write(tmpstr,'(A,3g14.7)') trim(subname)//' '//trim(lstring)//':'//trim(itemNameList(n))//'  ', &
             minval(dataPtr2d),maxval(dataPtr2d),sum(dataPtr2d)
           call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=lrc)
         else if(dimcount == 3) then
           call ESMF_FieldGet(lfield, farrayPtr=dataPtr3d, rc=lrc)
           if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

           write(tmpstr,'(A,3g14.7)') trim(subname)//' '//trim(lstring)//':'//trim(itemNameList(n))//'  ', &
             minval(dataPtr3d),maxval(dataPtr3d),sum(dataPtr3d)
           call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=lrc)
         else if(dimcount == 4) then
           call ESMF_FieldGet(lfield, farrayPtr=dataPtr4d, rc=lrc)
           if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

           write(tmpstr,'(A,3g14.7)') trim(subname)//' '//trim(lstring)//':'//trim(itemNameList(n))//'  ', &
             minval(dataPtr4d),maxval(dataPtr4d),sum(dataPtr4d)
           call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=lrc)
         else
           call ESMF_LogWrite('dimcount of >4 currently unsupported', ESMF_LOGMSG_INFO, rc=lrc)
         end if
       end if
     end if
    enddo
     deallocate(itemNameList)

    if (present(rc)) rc = lrc
    call ESMF_LogWrite(subname//' exit', ESMF_LOGMSG_INFO)

  end subroutine state_diagnose

  !-----------------------------------------------------------------------------

end module module_cap_cpl
