!> Return error to ESMF and finalize it.
#define ESMF_ERR_RETURN(rc) \
    if (ESMF_LogFoundError(rc, msg="Breaking out of subroutine", line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

module module_gather_test

  use mpi_f08
  use esmf

  implicit none
  private
  public gather_test

contains

  subroutine gather_test(wrtfb, rc)
!
    type(ESMF_FieldBundle), intent(in) :: wrtfb
    integer, optional,intent(out)      :: rc

!** local vars
    type(ESMF_VM)                      :: vm
    type(MPI_Comm)                     :: comm
    integer                            :: mype
    integer                            :: i,j,t, istart,iend,jstart,jend
    integer                            :: im, jm, lm, lsoil
    integer                            :: nproc, nproc_per_tile

    integer, dimension(:), allocatable              :: fldlev

    real(ESMF_KIND_R4), dimension(:,:), pointer     :: array_r4
    real(ESMF_KIND_R4), dimension(:,:,:), pointer   :: array_r4_3d

    real(ESMF_KIND_R8), dimension(:,:), pointer     :: array_r8
    real(ESMF_KIND_R8), dimension(:,:,:), pointer   :: array_r8_3d

    integer :: fieldCount, fieldDimCount, gridDimCount
    integer, dimension(:), allocatable   :: ungriddedLBound, ungriddedUBound

    type(ESMF_Field), allocatable        :: fcstField(:)
    type(ESMF_TypeKind_Flag)             :: typekind
    type(ESMF_TypeKind_Flag)             :: attTypeKind
    type(ESMF_Grid)                      :: wrtgrid
    type(ESMF_Array)                     :: array
    type(ESMF_DistGrid)                  :: distgrid

    integer :: attCount
    character(len=ESMF_MAXSTR) :: attName, fldName

    logical :: is_cubed_sphere
    logical :: is_cubed_sphere_tiled
    integer :: rank, deCount, localDeCount, dimCount, tileCount, tile_number, rootPet
    integer :: my_tile
    integer, dimension(:,:), allocatable :: minIndexPDe, maxIndexPDe
    integer, dimension(:,:), allocatable :: minIndexPTile, maxIndexPTile
    integer, dimension(:), allocatable :: deToTileMap, localDeToDeMap

    call ESMF_LogWrite("Enter gather_test", ESMF_LOGMSG_INFO, rc=rc)

    call ESMF_VMGetCurrent(vm=vm, rc=rc); ; ESMF_ERR_RETURN(rc)
    call ESMF_VMGet(vm=VM, localPet=mype, petCount=nproc, mpiCommunicator=comm%mpi_val, rc=rc); ESMF_ERR_RETURN(rc)

    is_cubed_sphere = .false.
    is_cubed_sphere_tiled = .false.
    tileCount = 0
    my_tile = 0
    rootPet = 0

    call ESMF_FieldBundleGet(wrtfb, fieldCount=fieldCount, rc=rc); ESMF_ERR_RETURN(rc)

    allocate(fldlev(fieldCount)) ; fldlev = 0
    allocate(fcstField(fieldCount))

    call ESMF_FieldBundleGet(wrtfb, fieldList=fcstField, grid=wrtgrid, &
                             rc=rc); ESMF_ERR_RETURN(rc)

    call ESMF_GridGet(wrtgrid, dimCount=gridDimCount, rc=rc); ESMF_ERR_RETURN(rc)

    do i=1,fieldCount
       call ESMF_FieldGet(fcstField(i), dimCount=fieldDimCount, array=array, rc=rc); ESMF_ERR_RETURN(rc)

       if (fieldDimCount > 3) then
          if (mype==0) write(0,*)"write_netcdf: Only 2D and 3D fields are supported!"
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
       end if

       ! use first field to determine tile number, grid size, start index etc.
       if (i == 1) then
          call ESMF_ArrayGet(array, &
                             distgrid=distgrid, &
                             dimCount=dimCount, &
                             deCount=deCount, &
                             localDeCount=localDeCount, &
                             tileCount=tileCount, &
                             rc=rc); ESMF_ERR_RETURN(rc)

          allocate(minIndexPDe(dimCount,deCount))
          allocate(maxIndexPDe(dimCount,deCount))
          allocate(minIndexPTile(dimCount, tileCount))
          allocate(maxIndexPTile(dimCount, tileCount))
          call ESMF_DistGridGet(distgrid, &
                                minIndexPDe=minIndexPDe, maxIndexPDe=maxIndexPDe, &
                                minIndexPTile=minIndexPTile, maxIndexPTile=maxIndexPTile, &
                                rc=rc); ESMF_ERR_RETURN(rc)

          allocate(deToTileMap(deCount))
          allocate(localDeToDeMap(localDeCount))
          call ESMF_ArrayGet(array, &
                             deToTileMap=deToTileMap, &
                             localDeToDeMap=localDeToDeMap, &
                             rc=rc); ESMF_ERR_RETURN(rc)

          is_cubed_sphere = (tileCount == 6)
          my_tile = deToTileMap(localDeToDeMap(1)+1)
          im = maxIndexPTile(1,1)
          jm = maxIndexPTile(2,1)
       end if

       if (fieldDimCount > gridDimCount) then
         allocate(ungriddedLBound(fieldDimCount-gridDimCount))
         allocate(ungriddedUBound(fieldDimCount-gridDimCount))
         call ESMF_FieldGet(fcstField(i), &
                            ungriddedLBound=ungriddedLBound, &
                            ungriddedUBound=ungriddedUBound, rc=rc); ESMF_ERR_RETURN(rc)
         fldlev(i) = ungriddedUBound(fieldDimCount-gridDimCount) - &
                     ungriddedLBound(fieldDimCount-gridDimCount) + 1
         deallocate(ungriddedLBound)
         deallocate(ungriddedUBound)
       else if (fieldDimCount == 2) then
         fldlev(i) = 1
       end if
    end do

    if (is_cubed_sphere) then
       is_cubed_sphere = .false.
       is_cubed_sphere_tiled = .true.
    end if

    tile_number = 1
    if (is_cubed_sphere_tiled) then
       if (mod(nproc,6) /=0) then
         if (mype==0) write(0,*)'For cubed sphere restarts, nproc must be divisible by 6. nproc = ', nproc
         call ESMF_Finalize(endflag=ESMF_END_ABORT)
       end if
       nproc_per_tile = nproc / 6
       tile_number = mype / nproc_per_tile + 1
       rootPet = (tile_number - 1) * nproc_per_tile
       if (tile_number /= my_tile) then
         if (mype==0) write(0,*)'Internal error: tile_number /= my_tile ', tile_number, my_tile
         call ESMF_Finalize(endflag=ESMF_END_ABORT)
       end if
    end if

    ! write variables (fields)
    do i=1, fieldCount

       call ESMF_FieldGet(fcstField(i),name=fldName,rank=rank,typekind=typekind, rc=rc); ESMF_ERR_RETURN(rc)

       call ESMF_LogWrite("gather "//trim(fldName), ESMF_LOGMSG_INFO, rc=rc)
       if (rank == 2) then

         if (typekind == ESMF_TYPEKIND_R4) then
                  allocate(array_r4(im,jm))
                  call ESMF_LogWrite(" Before FieldGather: "//trim(fldName), ESMF_LOGMSG_INFO, rc=rc)
                  call ESMF_VMBarrier(VM, rc=rc)
                  call ESMF_FieldGather(fcstField(i), array_r4, rootPet=rootPet, tile=my_tile, rc=rc); ESMF_ERR_RETURN(rc)
                  call ESMF_LogWrite(" After  FieldGather: "//trim(fldName), ESMF_LOGMSG_INFO, rc=rc)
                  deallocate(array_r4)
         else if (typekind == ESMF_TYPEKIND_R8) then
                  allocate(array_r8(im,jm))
                  call ESMF_LogWrite(" Before FieldGather: "//trim(fldName), ESMF_LOGMSG_INFO, rc=rc)
                  call ESMF_VMBarrier(VM, rc=rc)
                  call ESMF_FieldGather(fcstField(i), array_r8, rootPet=rootPet, tile=my_tile, rc=rc); ESMF_ERR_RETURN(rc)
                  call ESMF_LogWrite(" After  FieldGather: "//trim(fldName), ESMF_LOGMSG_INFO, rc=rc)
                  deallocate(array_r8)
         end if

      else if (rank == 3) then

         if (typekind == ESMF_TYPEKIND_R4) then
                  allocate(array_r4_3d(im,jm,fldlev(i)))
                  call ESMF_LogWrite(" Before FieldGather: "//trim(fldName), ESMF_LOGMSG_INFO, rc=rc)
                  call ESMF_VMBarrier(VM, rc=rc)
                  call ESMF_FieldGather(fcstField(i), array_r4_3d, rootPet=rootPet, tile=my_tile, rc=rc); ESMF_ERR_RETURN(rc)
                  call ESMF_LogWrite(" After  FieldGather: "//trim(fldName), ESMF_LOGMSG_INFO, rc=rc)
                  deallocate(array_r4_3d)
         else if (typekind == ESMF_TYPEKIND_R8) then
                  allocate(array_r8_3d(im,jm,fldlev(i)))
                  call ESMF_LogWrite(" Before FieldGather: "//trim(fldName), ESMF_LOGMSG_INFO, rc=rc)
                  call ESMF_VMBarrier(VM, rc=rc)
                  call ESMF_FieldGather(fcstField(i), array_r8_3d, rootPet=rootPet, tile=my_tile, rc=rc); ESMF_ERR_RETURN(rc)
                  call ESMF_LogWrite(" After  FieldGather: "//trim(fldName), ESMF_LOGMSG_INFO, rc=rc)
                  deallocate(array_r8_3d)
         end if ! end typekind

      else

         if (mype==0) write(0,*)'Unsupported rank ', rank
         call ESMF_Finalize(endflag=ESMF_END_ABORT)

      end if ! end rank

    end do ! end fieldCount
    deallocate(fcstField)

  end subroutine gather_test

!----------------------------------------------------------------------------------------
end module module_gather_test
