module esmf_utils

  use ESMF

  implicit none

  contains

  subroutine init_field_to_missing_value(field, rc)

    type(ESMF_Field), intent(inout)           :: field
    integer,          intent(out),   optional :: rc

    integer dimCount, rank
    type(ESMF_TypeKind_Flag) typekind
    character(len=ESMF_MAXSTR) fieldName
    real(ESMF_KIND_R4), dimension(:,:),     pointer  :: ptr2d_r4
    real(ESMF_KIND_R4), dimension(:,:,:),   pointer  :: ptr3d_r4
    real(ESMF_KIND_R4), dimension(:,:,:,:), pointer  :: ptr4d_r4
    real(ESMF_KIND_R8), dimension(:,:),     pointer  :: ptr2d_r8
    real(ESMF_KIND_R8), dimension(:,:,:),   pointer  :: ptr3d_r8
    real(ESMF_KIND_R8), dimension(:,:,:,:), pointer  :: ptr4d_r8

    real(ESMF_KIND_R4) :: missing_value_r4=9.99e20
    real(ESMF_KIND_R8) :: missing_value_r8=9.99e20

    call ESMF_FieldGet(field, name=fieldName, typekind=typekind, dimCount=dimCount, rank=rank, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

    ! Only work on ESMF_TYPEKIND_R4 fields for now
    if (typekind == ESMF_TYPEKIND_R4) then
      if (dimCount == 2) then
        call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr2d_r4, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
        ptr2d_r4 = missing_value_r4
      else if (dimCount == 3) then
        call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr3d_r4, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
        ptr3d_r4 = missing_value_r4
      else if (dimCount == 4) then
        call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr4d_r4, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
        ptr4d_r4 = missing_value_r4
      else
        write(0,*)' Unsupported dimCount = ', dimCount
        stop
      endif
    else if (typekind == ESMF_TYPEKIND_R8) then
      if (dimCount == 2) then
        call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr2d_r8, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
        ptr2d_r8 = missing_value_r8
      else if (dimCount == 3) then
        call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr3d_r8, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
        ptr3d_r8 = missing_value_r8
      else if (dimCount == 4) then
        call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr4d_r8, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
        ptr4d_r8 = missing_value_r8
      else
        write(0,*)' Unsupported dimCount = ', dimCount
        stop
      endif
    else
      write(0,*)' Unsupported typekind = ', typekind
      stop
    endif

    rc = 0

  end subroutine init_field_to_missing_value

  subroutine generate_dst_field_mask(src_grid, dst_grid, dst_field, rc)

    type(ESMF_Grid), intent(in)    :: src_grid, dst_grid
    type(ESMF_Field), intent(inout):: dst_field
    integer, intent(out)           :: rc

    type(ESMF_Field)               :: src_field
    type(ESMF_Field)               :: dst_status_field
    real(ESMF_KIND_R4), pointer    :: ptr(:,:)
    integer(ESMF_KIND_I4), pointer :: maskPtr(:,:)
    integer(ESMF_KIND_I4), pointer :: ptr_dst_status(:,:)
    type(ESMF_RouteHandle)         :: routehandle_mask
    character(ESMF_MAXSTR)         :: itemName
    integer                        :: localDeCount
    integer                        :: ig,jg, istart,iend, jstart,jend
    integer                        :: srcTermProcessing


    src_field = ESMF_FieldCreate(src_grid, &
                                 typekind=ESMF_TYPEKIND_R4, &
                                 staggerloc=ESMF_STAGGERLOC_CENTER, &
                                 rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_FieldGet(src_field, localDeCount=localDeCount, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    if (localDeCount > 0) then
      call ESMF_FieldGet(src_field, farrayPtr=ptr, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      ptr = 1.0
    end if


    dst_field = ESMF_FieldCreate(dst_grid, &
                                 typekind=ESMF_TYPEKIND_R4, &
                                 staggerloc=ESMF_STAGGERLOC_CENTER, &
                                 rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_FieldGet(dst_field, localDeCount=localDeCount, rc=rc);
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    if (localDeCount > 0) then
      call ESMF_FieldGet(dst_field, farrayPtr=ptr, rc=rc);
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      ptr = 0.0
    end if

    dst_status_field = ESMF_FieldCreate(dst_grid, ESMF_TYPEKIND_I4, name='dst_status_field', rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    srcTermProcessing = 0

    call ESMF_FieldRegridStore(srcField=src_field, &
                               dstField=dst_field, &
                               regridmethod=ESMF_REGRIDMETHOD_BILINEAR, &
                               routehandle=routehandle_mask, &
                               unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, &
                               srcTermProcessing=srcTermProcessing, &
                               dstStatusField=dst_status_field, &
                               rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_FieldRegrid(src_field, dst_field, &
                          routehandle=routehandle_mask, &
                          termorderflag=ESMF_TERMORDER_SRCSEQ, &
                          zeroregion=ESMF_REGION_SELECT, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_FieldGet(dst_field, localDeCount=localDeCount, rc=rc);
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    if (localDeCount > 0) then
      call ESMF_FieldGet(dst_field, farrayPtr=ptr, rc=rc);
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call ESMF_FieldGet(dst_status_field, farrayPtr=ptr_dst_status, rc=rc);
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      istart = lbound(ptr,1)
      iend   = ubound(ptr,1)
      jstart = lbound(ptr,2)
      jend   = ubound(ptr,2)
      do jg=jstart, jend
        do ig=istart, iend
          if (ptr(ig,jg) == 0.0 .and. ptr_dst_status(ig,jg) /= ESMF_REGRIDSTATUS_OUTSIDE) then
            write(0,*)'remap warning: ', ig,jg,ptr(ig,jg), ptr_dst_status(ig,jg)
          else if (ptr(ig,jg) /= 0.0 .and. ptr_dst_status(ig,jg) /= ESMF_REGRIDSTATUS_MAPPED) then
            write(0,*)'remap warning: ', ig,jg,ptr(ig,jg), ptr_dst_status(ig,jg)
          endif
        enddo
      enddo
    end if

    call ESMF_RouteHandleDestroy(routehandle_mask, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_FieldDestroy(src_field, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_FieldDestroy(dst_status_field, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    rc = 0
  end subroutine generate_dst_field_mask

  subroutine add_dst_mask(dst_grid, dst_field, dstOutsideMaskValue, rc)

    type(ESMF_Grid), intent(inout) :: dst_grid
    type(ESMF_Field), intent(in)   :: dst_field
    integer, intent(in)            :: dstOutsideMaskValue
    integer, intent(out)           :: rc

    real(ESMF_KIND_R4), pointer    :: ptr(:,:)
    integer(ESMF_KIND_I4), pointer :: maskPtr(:,:)
    integer                        :: localDeCount


    call ESMF_FieldGet(dst_field, localDeCount=localDeCount, rc=rc);
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! Set destination grid mask to dstOutsideMaskValue where dst_field is 0 which will be for all destination
    ! points outside of source (forecast, computational) grid
    if (localDeCount > 0) then
      call ESMF_FieldGet(dst_field, farrayPtr=ptr, rc=rc);
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call ESMF_GridGetItem(dst_grid, itemflag=ESMF_GRIDITEM_MASK,   &
                            staggerloc=ESMF_STAGGERLOC_CENTER, farrayPtr=maskPtr, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      maskPtr = 0
      where (ptr == 0.0)
        maskPtr = dstOutsideMaskValue
      endwhere
    end if

    rc = 0

  end subroutine add_dst_mask

end module esmf_utils
