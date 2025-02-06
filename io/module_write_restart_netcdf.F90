#define ESMF_ERR_RETURN(rc) \
    if (ESMF_LogFoundError(rc, msg="Breaking out of subroutine", line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

#define NC_ERR_STOP(status) \
    if (status /= nf90_noerr) write(0,*) "file: ", __FILE__, " line: ", __LINE__, trim(nf90_strerror(status)); \
    if (status /= nf90_noerr) call ESMF_Finalize(endflag=ESMF_END_ABORT)

module module_write_restart_netcdf

  use mpi_f08
  use esmf
  use fms
  use mpp_mod, only : mpp_chksum   ! needed for fms 2023.02
  use netcdf
  use module_fv3_io_def,only : zstandard_level

  implicit none
  private
  public write_restart_netcdf

  logical :: par

  contains

!----------------------------------------------------------------------------------------
  subroutine write_restart_netcdf(wrtfb, filename, &
                                  use_parallel_netcdf, comm, mype, &
                                  rc)
!
    type(ESMF_FieldBundle), intent(in) :: wrtfb
    character(*), intent(in)           :: filename
    logical, intent(in)                :: use_parallel_netcdf
    type(MPI_Comm), intent(in)         :: comm
    integer, intent(in)                :: mype
    integer, optional,intent(out)      :: rc
!
!** local vars
    integer, parameter :: i8_kind=8
    integer :: i,j,k,t, istart,iend,jstart,jend
    integer :: im, jm, lm

    integer, dimension(:), allocatable              :: fldlev

    real(ESMF_KIND_R4), dimension(:,:), pointer     :: array_r4
    real(ESMF_KIND_R4), dimension(:,:,:), pointer   :: array_r4_cube
    real(ESMF_KIND_R4), dimension(:,:,:), pointer   :: array_r4_3d
    real(ESMF_KIND_R4), dimension(:,:,:,:), pointer :: array_r4_3d_cube

    real(ESMF_KIND_R8), dimension(:,:), pointer     :: array_r8
    real(ESMF_KIND_R8), dimension(:,:,:), pointer   :: array_r8_cube
    real(ESMF_KIND_R8), dimension(:,:,:), pointer   :: array_r8_3d
    real(ESMF_KIND_R8), dimension(:,:,:,:), pointer :: array_r8_3d_cube

    real(8), dimension(:), allocatable :: x,y
    integer :: fieldCount, fieldDimCount, gridDimCount
    integer, dimension(:), allocatable   :: ungriddedLBound, ungriddedUBound
    integer, dimension(:), allocatable   :: start_idx

    type(ESMF_Field), allocatable        :: fcstField(:)
    type(ESMF_TypeKind_Flag)             :: typekind
    type(ESMF_StaggerLoc)                :: staggerloc
    type(ESMF_TypeKind_Flag)             :: attTypeKind
    type(ESMF_Grid)                      :: wrtgrid
    type(ESMF_Array)                     :: array
    type(ESMF_DistGrid)                  :: distgrid

    integer :: attCount
    character(len=ESMF_MAXSTR) :: attName, fldName

    integer :: varival
    real(4) :: varr4val, dataMin, dataMax
    real(4), allocatable, dimension(:) :: compress_err
    real(8) :: varr8val
    character(len=ESMF_MAXSTR) :: varcval

    integer :: ncerr,ierr
    integer :: ncid
    integer :: oldMode
    integer :: dimid, dimtype
    integer :: im_dimid, im_p1_dimid, jm_dimid, jm_p1_dimid, time_dimid
    integer :: im_varid, im_p1_varid, jm_varid, jm_p1_varid, time_varid
    integer, dimension(:), allocatable :: dimids_2d, dimids_3d, chunksizes
    integer, dimension(:), allocatable :: varids, zaxis_dimids
    logical shuffle

    logical :: get_im_jm, found_im_jm
    integer :: rank, deCount, localDeCount, dimCount, tileCount
    integer :: my_tile, start_i, start_j
    integer, dimension(:,:), allocatable :: minIndexPDe, maxIndexPDe
    integer, dimension(:,:), allocatable :: minIndexPTile, maxIndexPTile
    integer, dimension(:), allocatable :: deToTileMap, localDeToDeMap
    logical :: do_io
    integer :: par_access

    real(ESMF_KIND_R8), allocatable  :: valueListr8(:)
    logical :: isPresent, thereAreVerticals, is_restart_core, dynamics_restart_file
    integer :: udimCount
    character(80), allocatable :: udimList(:)
    ! character(32), allocatable :: field_checksums(:)
    character(32) :: axis_attr_name
    character(32) :: field_checksum
!
    is_restart_core = .false.
    if ( index(trim(filename),'fv_core.res') > 0 ) is_restart_core = .true.
    tileCount = 0
    my_tile = 0
    start_i = -10000000
    start_j = -10000000

    par = use_parallel_netcdf
    do_io = par .or. (mype==0)

    call ESMF_FieldBundleGet(wrtfb, fieldCount=fieldCount, rc=rc); ESMF_ERR_RETURN(rc)

    allocate(compress_err(fieldCount)); compress_err=-999.
    allocate(fldlev(fieldCount)) ; fldlev = 0
    allocate(fcstField(fieldCount))
    allocate(varids(fieldCount))
    allocate(zaxis_dimids(fieldCount))
    ! allocate(field_checksums(fieldCount))

    call ESMF_FieldBundleGet(wrtfb, fieldList=fcstField, grid=wrtGrid, &
                             itemorderflag=ESMF_ITEMORDER_ADDORDER, &
                             rc=rc); ESMF_ERR_RETURN(rc)

    call ESMF_GridGet(wrtgrid, dimCount=gridDimCount, rc=rc); ESMF_ERR_RETURN(rc)

    found_im_jm = .false.
    do i=1,fieldCount
       call ESMF_FieldGet(fcstField(i), dimCount=fieldDimCount, array=array, rc=rc); ESMF_ERR_RETURN(rc)
       call ESMF_FieldGet(fcstField(i), name=fldName, rank=rank, typekind=typekind, rc=rc); ESMF_ERR_RETURN(rc)

       if (fieldDimCount > 3) then
          if (mype==0) write(0,*)"write_netcdf: Only 2D and 3D fields are supported!"
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
       end if

       ! use first field to determine tile number, grid size, start index etc.
       if (is_restart_core) then
          get_im_jm = trim(fldName) /= 'u' .and. trim(fldName) /= 'v' ! skip staggered fields
       else
          get_im_jm = (i==1)
       end if
       if ( .not.found_im_jm .and. get_im_jm) then
          call ESMF_ArrayGet(array, &
                             distgrid=distgrid, &
                             dimCount=dimCount, &
                             deCount=deCount, &
                             localDeCount=localDeCount, &
                             tileCount=tileCount, &
                             rc=rc); ESMF_ERR_RETURN(rc)

          if (tileCount /= 1) then
             if (mype==0) write(0,*)"write_restart_netcdf: Only tileCount==1 fields are supported!"
             call ESMF_Finalize(endflag=ESMF_END_ABORT)
          endif

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

          my_tile = deToTileMap(localDeToDeMap(1)+1)
          im = maxIndexPTile(1,1) - minIndexPTile(1,1) + 1
          jm = maxIndexPTile(2,1) - minIndexPTile(2,1) + 1
          start_i = minIndexPDe(1,localDeToDeMap(1)+1)
          start_j = minIndexPDe(2,localDeToDeMap(1)+1)
          if (.not. par) then
             start_i = 1
             start_j = 1
          end if
          found_im_jm = .true.
          start_i = mod(start_i, im) ! for cubed sphere multi-tile grid
          start_j = mod(start_j, jm) ! for cubed sphere multi-tile grid

          deallocate(minIndexPDe)
          deallocate(maxIndexPDe)
          deallocate(minIndexPTile)
          deallocate(maxIndexPTile)
          deallocate(deToTileMap)
          deallocate(localDeToDeMap)

          if (typekind == ESMF_TYPEKIND_R4) then
             dimtype = NF90_FLOAT
          else if (typekind == ESMF_TYPEKIND_R8) then
             dimtype = NF90_DOUBLE
          else
             if (mype==0) write(0,*)'Unsupported typekind ', typekind
             call ESMF_Finalize(endflag=ESMF_END_ABORT)
          end if
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

    ! create netcdf file and enter define mode
    if (do_io) then

       if (par) then
          ncerr = nf90_create(trim(filename),&
                  cmode=IOR(NF90_CLOBBER,NF90_NETCDF4),&
                  comm=comm%mpi_val, info = MPI_INFO_NULL%mpi_val, ncid=ncid); NC_ERR_STOP(ncerr)
       else
          ncerr = nf90_create(trim(filename),&
                  ! cmode=IOR(NF90_CLOBBER,NF90_64BIT_OFFSET),&
                  cmode=IOR(NF90_CLOBBER,NF90_NETCDF4),&
                  ncid=ncid); NC_ERR_STOP(ncerr)
       end if

       ! disable auto filling.
       ncerr = nf90_set_fill(ncid, NF90_NOFILL, oldMode); NC_ERR_STOP(ncerr)

       dynamics_restart_file = index(trim(filename),"fv_") > 0

       ! Unnecessary naming inconsistency
       if (dynamics_restart_file) then
          axis_attr_name = "axis"
       else
          axis_attr_name = "cartesian_axis"
       end if

       ! define dimensions [xaxis_1, yaxis_1 ,(zaxis_1,...), Time]
       if ( .not.is_restart_core ) then

          ncerr = nf90_def_dim(ncid, "xaxis_1", im, im_dimid); NC_ERR_STOP(ncerr)
          ncerr = nf90_def_var(ncid, "xaxis_1", dimtype, im_dimid, im_varid); NC_ERR_STOP(ncerr)
          ncerr = nf90_put_att(ncid, im_varid, trim(axis_attr_name), "X"); NC_ERR_STOP(ncerr)

          ncerr = nf90_def_dim(ncid, "yaxis_1", jm, jm_dimid); NC_ERR_STOP(ncerr)
          ncerr = nf90_def_var(ncid, "yaxis_1", dimtype, jm_dimid, jm_varid); NC_ERR_STOP(ncerr)
          ncerr = nf90_put_att(ncid, jm_varid, trim(axis_attr_name), "Y"); NC_ERR_STOP(ncerr)

       else

          ncerr = nf90_def_dim(ncid, "xaxis_1", im, im_dimid); NC_ERR_STOP(ncerr)
          ncerr = nf90_def_var(ncid, "xaxis_1", dimtype, im_dimid, im_varid); NC_ERR_STOP(ncerr)
          ncerr = nf90_put_att(ncid, im_varid, trim(axis_attr_name), "X"); NC_ERR_STOP(ncerr)

          ncerr = nf90_def_dim(ncid, "xaxis_2", im+1, im_p1_dimid); NC_ERR_STOP(ncerr)
          ncerr = nf90_def_var(ncid, "xaxis_2", dimtype, im_p1_dimid, im_p1_varid); NC_ERR_STOP(ncerr)
          ncerr = nf90_put_att(ncid, im_p1_varid, trim(axis_attr_name), "X"); NC_ERR_STOP(ncerr)

          ncerr = nf90_def_dim(ncid, "yaxis_1", jm+1, jm_p1_dimid); NC_ERR_STOP(ncerr)
          ncerr = nf90_def_var(ncid, "yaxis_1", dimtype, jm_p1_dimid, jm_p1_varid); NC_ERR_STOP(ncerr)
          ncerr = nf90_put_att(ncid, jm_p1_varid, trim(axis_attr_name), "Y"); NC_ERR_STOP(ncerr)

          ncerr = nf90_def_dim(ncid, "yaxis_2", jm, jm_dimid); NC_ERR_STOP(ncerr)
          ncerr = nf90_def_var(ncid, "yaxis_2", dimtype, jm_dimid, jm_varid); NC_ERR_STOP(ncerr)
          ncerr = nf90_put_att(ncid, jm_varid, trim(axis_attr_name), "Y"); NC_ERR_STOP(ncerr)

       end if

       ! define ungridded (vertical) coordinate variables
       do i=1,fieldCount
          call ESMF_FieldGet(fcstField(i), name=fldName, rank=rank, typekind=typekind, rc=rc); ESMF_ERR_RETURN(rc)
          call ESMF_AttributeGetAttPack(fcstField(i), convention="NetCDF", purpose="FV3", isPresent=isPresent, rc=rc); ESMF_ERR_RETURN(rc)
          if (isPresent) then
            call ESMF_AttributeGet(fcstField(i), convention="NetCDF", purpose="FV3", &
                                   name="ESMF:ungridded_dim_labels", isPresent=isPresent, &
                                   itemCount=udimCount, rc=rc); ESMF_ERR_RETURN(rc)
            if (udimCount>1) then
              if (mype==0) write(0,*)'udimCount>1 for ', trim(fldName)
              ESMF_ERR_RETURN(-1)
            end if
            if (udimCount>0 .and. isPresent) then
              thereAreVerticals = .true.
              allocate(udimList(udimCount))
              call ESMF_AttributeGet(fcstField(i), convention="NetCDF", purpose="FV3", name="ESMF:ungridded_dim_labels", valueList=udimList, rc=rc); ESMF_ERR_RETURN(rc)
              ! loop over all ungridded dimension labels
              do k=1, udimCount
                call write_out_ungridded_dim_atts_from_field(fcstField(i), trim(udimList(k)), dimid, rc=rc); ESMF_ERR_RETURN(rc)
                zaxis_dimids(i) = dimid
              enddo
              deallocate(udimList)
            end if
          end if
       end do

       ncerr = nf90_def_dim(ncid, "Time", NF90_UNLIMITED, time_dimid); NC_ERR_STOP(ncerr)
       ! ncerr = nf90_def_dim(ncid, "Time", 1, time_dimid); NC_ERR_STOP(ncerr)
       ncerr = nf90_def_var(ncid, "Time", dimtype, time_dimid, time_varid); NC_ERR_STOP(ncerr)
       if (par) then
          ncerr = nf90_var_par_access(ncid, time_varid, NF90_COLLECTIVE); NC_ERR_STOP(ncerr)
       end if

       ! Again, unnecessary naming inconsistency
       if (dynamics_restart_file) then
          ncerr = nf90_put_att(ncid, time_varid, "cartesian_axis", "T"); NC_ERR_STOP(ncerr)
          ncerr = nf90_put_att(ncid, time_varid, "units", "time level"); NC_ERR_STOP(ncerr)
          ncerr = nf90_put_att(ncid, time_varid, "long_name", "Time"); NC_ERR_STOP(ncerr)
       else
          ncerr = nf90_put_att(ncid, time_varid, trim(axis_attr_name), "T"); NC_ERR_STOP(ncerr)
       end if

       ncerr = nf90_redef(ncid=ncid)

       ! define variables (fields)
       allocate(dimids_2d(3))
       allocate(dimids_3d(4))

       do i=1, fieldCount
         call ESMF_FieldGet(fcstField(i), name=fldName, rank=rank, typekind=typekind, staggerloc=staggerloc, rc=rc); ESMF_ERR_RETURN(rc)

         ! FIXME remove this once ESMF_FieldGet above returns correct staggerloc
         staggerloc = ESMF_STAGGERLOC_CENTER
         if (is_restart_core .and. trim(fldName) == 'u' ) staggerloc = ESMF_STAGGERLOC_EDGE2
         if (is_restart_core .and. trim(fldName) == 'v' ) staggerloc = ESMF_STAGGERLOC_EDGE1

         ! par_access = NF90_INDEPENDENT
         par_access = NF90_COLLECTIVE   ! because of time unlimited
         ! define variables
         if (rank == 2) then
           dimids_2d =             [im_dimid,jm_dimid,                       time_dimid]
           chunksizes =            [im, jm, 1]
           if (typekind == ESMF_TYPEKIND_R4) then
             ncerr = nf90_def_var(ncid, trim(fldName), NF90_FLOAT, dimids_2d, varids(i)); NC_ERR_STOP(ncerr)
           else if (typekind == ESMF_TYPEKIND_R8) then
             ncerr = nf90_def_var(ncid, trim(fldName), NF90_DOUBLE, dimids_2d, varids(i)); NC_ERR_STOP(ncerr)
           else
             if (mype==0) write(0,*)'Unsupported typekind ', typekind
             call ESMF_Finalize(endflag=ESMF_END_ABORT)
           end if
         else if (rank == 3) then
           if ( .not.is_restart_core ) then
             dimids_3d = [im_dimid,jm_dimid,zaxis_dimids(i),time_dimid]
             chunksizes = [im, jm, 1, 1]
           else
             if (staggerloc == ESMF_STAGGERLOC_CENTER) then
                dimids_3d = [im_dimid,jm_dimid,zaxis_dimids(i),time_dimid]
                chunksizes = [im, jm, 1, 1]
             else if (staggerloc == ESMF_STAGGERLOC_EDGE1) then  ! east
                dimids_3d = [im_p1_dimid,jm_dimid,zaxis_dimids(i),time_dimid]
                chunksizes = [im+1, jm, 1, 1]
             else if (staggerloc == ESMF_STAGGERLOC_EDGE2) then  ! south
                dimids_3d = [im_dimid,jm_p1_dimid,zaxis_dimids(i),time_dimid]
                chunksizes = [im, jm+1, 1, 1]
             else
               if (mype==0) write(0,*)'Unsupported staggerloc ', staggerloc
               call ESMF_Finalize(endflag=ESMF_END_ABORT)
             end if
           end if
           if (typekind == ESMF_TYPEKIND_R4) then
             ncerr = nf90_def_var(ncid, trim(fldName), NF90_FLOAT, dimids_3d, varids(i)); NC_ERR_STOP(ncerr)
           else if (typekind == ESMF_TYPEKIND_R8) then
             ncerr = nf90_def_var(ncid, trim(fldName), NF90_DOUBLE, dimids_3d, varids(i)); NC_ERR_STOP(ncerr)
           else
             if (mype==0) write(0,*)'Unsupported typekind ', typekind
             call ESMF_Finalize(endflag=ESMF_END_ABORT)
           end if
         else
           if (mype==0) write(0,*)'Unsupported rank ', rank
           call ESMF_Finalize(endflag=ESMF_END_ABORT)
         end if
         if (par) then
             ncerr = nf90_var_par_access(ncid, varids(i), par_access); NC_ERR_STOP(ncerr)
         end if

         ncerr = nf90_def_var_chunking(ncid, varids(i), NF90_CHUNKED, chunksizes) ; NC_ERR_STOP(ncerr)

         if (zstandard_level(1) > 0) then
            ncerr = nf90_def_var_zstandard(ncid, varids(i), zstandard_level(1))
            if (ncerr /= nf90_noerr) then
               if (ncerr == nf90_enofilter) then
                  if (mype==0) write(0,*) 'Zstandard filter not found.'
               end if
               NC_ERR_STOP(ncerr)
            end if
         end if

       end do   ! i=1,fieldCount

       ncerr = nf90_put_att(ncid, NF90_GLOBAL, "NumFilesInSet", 1); NC_ERR_STOP(ncerr)

! end of define mode
       ncerr = nf90_enddef(ncid); NC_ERR_STOP(ncerr)

       if (allocated(valueListr8)) deallocate(valueListr8)
       allocate (valueListr8(im))
       valueListr8 = (/(i, i=1,im)/)
       ncerr = nf90_put_var(ncid, im_varid, values=valueListr8); NC_ERR_STOP(ncerr)

       if (allocated(valueListr8)) deallocate(valueListr8)
       allocate (valueListr8(jm))
       valueListr8 = (/(i, i=1,jm)/)
       ncerr = nf90_put_var(ncid, jm_varid, values=valueListr8); NC_ERR_STOP(ncerr)

       if ( is_restart_core ) then
          if (allocated(valueListr8)) deallocate(valueListr8)
          allocate (valueListr8(im+1))
          valueListr8 = (/(i, i=1,im+1)/)
          ncerr = nf90_put_var(ncid, im_p1_varid, values=valueListr8); NC_ERR_STOP(ncerr)

          if (allocated(valueListr8)) deallocate(valueListr8)
          allocate (valueListr8(jm+1))
          valueListr8 = (/(i, i=1,jm+1)/)
          ncerr = nf90_put_var(ncid, jm_p1_varid, values=valueListr8); NC_ERR_STOP(ncerr)
       end if

       ncerr = nf90_put_var(ncid, time_varid, values=[1]); NC_ERR_STOP(ncerr)
    end if

    ! write variables (fields)
    do i=1, fieldCount

       call ESMF_FieldGet(fcstField(i),name=fldName,rank=rank,typekind=typekind,staggerloc=staggerloc, rc=rc); ESMF_ERR_RETURN(rc)

       ! FIXME remove this once ESMF_FieldGet above returns correct staggerloc
       staggerloc = ESMF_STAGGERLOC_CENTER
       if (is_restart_core .and. trim(fldName) == 'u' ) staggerloc = ESMF_STAGGERLOC_EDGE2
       if (is_restart_core .and. trim(fldName) == 'v' ) staggerloc = ESMF_STAGGERLOC_EDGE1

       if (rank == 2) then

         if (allocated(start_idx)) deallocate(start_idx)
         allocate(start_idx(3))
         start_idx = [start_i,start_j,        1]

         if (typekind == ESMF_TYPEKIND_R4) then
            call ESMF_FieldGet(fcstField(i), localDe=0, farrayPtr=array_r4, rc=rc); ESMF_ERR_RETURN(rc)
            write(field_checksum,'(Z16)') mpp_chksum(array_r4)
            if (par) then
               ncerr = nf90_put_var(ncid, varids(i), values=array_r4, start=start_idx); NC_ERR_STOP(ncerr)
            else
               allocate(array_r4(im,jm))
               call ESMF_FieldGather(fcstField(i), array_r4, rootPet=0, rc=rc); ESMF_ERR_RETURN(rc)
               if (do_io) then
                  ncerr = nf90_put_var(ncid, varids(i), values=array_r4, start=start_idx); NC_ERR_STOP(ncerr)
               end if
               deallocate(array_r4)
            end if
         else if (typekind == ESMF_TYPEKIND_R8) then
            call ESMF_FieldGet(fcstField(i), localDe=0, farrayPtr=array_r8, rc=rc); ESMF_ERR_RETURN(rc)
            write(field_checksum,'(Z16)') mpp_chksum(array_r8)
            if (par) then
               ncerr = nf90_put_var(ncid, varids(i), values=array_r8, start=start_idx); NC_ERR_STOP(ncerr)
            else
               allocate(array_r8(im,jm))
               call ESMF_FieldGather(fcstField(i), array_r8, rootPet=0, rc=rc); ESMF_ERR_RETURN(rc)
               if (do_io) then
                  ncerr = nf90_put_var(ncid, varids(i), values=array_r8, start=start_idx); NC_ERR_STOP(ncerr)
               end if
               deallocate(array_r8)
            end if
         end if

      else if (rank == 3) then

         if (allocated(start_idx)) deallocate(start_idx)
         allocate(start_idx(4))
         start_idx = [start_i,start_j,1,        1]

         if (typekind == ESMF_TYPEKIND_R4) then
            call ESMF_FieldGet(fcstField(i), localDe=0, farrayPtr=array_r4_3d, rc=rc); ESMF_ERR_RETURN(rc)
            write(field_checksum,'(Z16)') mpp_chksum(array_r4_3d)
            if (par) then
               ncerr = nf90_put_var(ncid, varids(i), values=array_r4_3d, start=start_idx); NC_ERR_STOP(ncerr)
            else
               if (staggerloc == ESMF_STAGGERLOC_CENTER) then
                 allocate(array_r4_3d(im,jm,fldlev(i)))
               else if (staggerloc == ESMF_STAGGERLOC_EDGE1) then  ! east
                 allocate(array_r4_3d(im+1,jm,fldlev(i)))
               else if (staggerloc == ESMF_STAGGERLOC_EDGE2) then  ! south
                 allocate(array_r4_3d(im,jm+1,fldlev(i)))
               else
                 if (mype==0) write(0,*)'Unsupported staggerloc ', staggerloc
                 call ESMF_Finalize(endflag=ESMF_END_ABORT)
               end if
               call ESMF_FieldGather(fcstField(i), array_r4_3d, rootPet=0, rc=rc); ESMF_ERR_RETURN(rc)
               if (mype==0) then
                  ncerr = nf90_put_var(ncid, varids(i), values=array_r4_3d, start=start_idx); NC_ERR_STOP(ncerr)
               end if
               deallocate(array_r4_3d)
            end if
         else if (typekind == ESMF_TYPEKIND_R8) then
            call ESMF_FieldGet(fcstField(i), localDe=0, farrayPtr=array_r8_3d, rc=rc); ESMF_ERR_RETURN(rc)
            write(field_checksum,'(Z16)') mpp_chksum(array_r8_3d)
            if (par) then
               ncerr = nf90_put_var(ncid, varids(i), values=array_r8_3d, start=start_idx); NC_ERR_STOP(ncerr)
            else
               if (staggerloc == ESMF_STAGGERLOC_CENTER) then
                 allocate(array_r8_3d(im,jm,fldlev(i)))
               else if (staggerloc == ESMF_STAGGERLOC_EDGE1) then  ! east
                 allocate(array_r8_3d(im+1,jm,fldlev(i)))
               else if (staggerloc == ESMF_STAGGERLOC_EDGE2) then  ! south
                 allocate(array_r8_3d(im,jm+1,fldlev(i)))
               else
                 if (mype==0) write(0,*)'Unsupported staggerloc ', staggerloc
                 call ESMF_Finalize(endflag=ESMF_END_ABORT)
               end if
               call ESMF_FieldGather(fcstField(i), array_r8_3d, rootPet=0, rc=rc); ESMF_ERR_RETURN(rc)
               if (mype==0) then
                  ncerr = nf90_put_var(ncid, varids(i), values=array_r8_3d, start=start_idx); NC_ERR_STOP(ncerr)
               end if
               deallocate(array_r8_3d)
            end if
         end if ! end typekind

      else

         if (mype==0) write(0,*)'Unsupported rank ', rank
         call ESMF_Finalize(endflag=ESMF_END_ABORT)

      end if ! end rank
      if (do_io) then
         ncerr = nf90_redef(ncid=ncid); NC_ERR_STOP(ncerr)
         ncerr = nf90_put_att(ncid, varids(i), 'checksum', field_checksum); NC_ERR_STOP(ncerr)
         ncerr = nf90_enddef(ncid=ncid); NC_ERR_STOP(ncerr)
      end if

    end do ! end fieldCount

    if (do_io) then
       ncerr = nf90_close(ncid=ncid); NC_ERR_STOP(ncerr)
    end if

contains

    subroutine write_out_ungridded_dim_atts_from_field(field, dimLabel, dimid, rc)

      type(ESMF_Field),intent(in) :: field
      character(len=*),intent(in) :: dimLabel
      integer, intent(out) :: dimid
      integer, intent(out)  :: rc

      real(ESMF_KIND_R4), allocatable  :: valueListr4(:)
      real(ESMF_KIND_R8), allocatable  :: valueListr8(:)
      integer                          :: valueCount, udimCount
      integer                          :: ncerr, varid, ind
      integer                          :: itemCount, attCount
      character(len=80),  allocatable  :: attNameList(:)
      character(len=80)                :: attName
      type(ESMF_TypeKind_Flag)         :: typekind
      character(len=80)                :: valueS
      integer                          :: valueI4
      real(ESMF_KIND_R4)               :: valueR4
      real(ESMF_KIND_R8)               :: valueR8

      ! inquire if NetCDF file already contains this ungridded dimension variable
      ncerr = nf90_inq_varid(ncid, trim(dimLabel), varid=varid)
      if (ncerr == NF90_NOERR) then
         ! if it does inquire dimid, it must be defined already
         ncerr = nf90_inq_dimid(ncid, trim(dimLabel), dimid=dimid)
         if (ncerr == NF90_NOERR) then
           return
         else
           if (mype==0) write(0,*) 'in write_out_ungridded_dim_atts: ERROR missing dimid for already defined ungridded dimension variable '
           ESMF_ERR_RETURN(-1)
        endif
      endif

      ! the variable does not exist in the NetCDF file yet -> add it
      ! access the undistributed dimension attribute on the grid
      call ESMF_AttributeGet(field, convention="NetCDF", purpose="FV3-dim", name=trim(dimLabel), itemCount=valueCount, typekind=typekind, rc=rc); ESMF_ERR_RETURN(rc)
      if( typekind == ESMF_TYPEKIND_R4 ) then
        allocate(valueListr4(valueCount))
        call ESMF_AttributeGet(field, convention="NetCDF", purpose="FV3-dim", name=trim(dimLabel), valueList=valueListr4, rc=rc); ESMF_ERR_RETURN(rc)
      else if ( typekind == ESMF_TYPEKIND_R8) then
        allocate(valueListr8(valueCount))
        call ESMF_AttributeGet(field, convention="NetCDF", purpose="FV3-dim", name=trim(dimLabel), valueList=valueListr8, rc=rc); ESMF_ERR_RETURN(rc)
      else
        if (mype==0) write(0,*) 'in write_out_ungridded_dim_atts: ERROR unknown typekind'
        ESMF_ERR_RETURN(-1)
      endif

      ! now add it to the NetCDF file
      ncerr = nf90_inq_dimid(ncid, trim(dimLabel), dimid=dimid)
      if (ncerr /= NF90_NOERR) then
        ! dimension does not yet exist, and must be defined
        ncerr = nf90_def_dim(ncid, trim(dimLabel), valueCount, dimid=dimid); NC_ERR_STOP(ncerr); NC_ERR_STOP(ncerr)
      endif
      if( typekind == ESMF_TYPEKIND_R4 ) then
        ncerr = nf90_def_var(ncid, trim(dimLabel), NF90_FLOAT, dimids=(/dimid/), varid=varid); NC_ERR_STOP(ncerr)
        ncerr = nf90_put_att(ncid, varid, trim(axis_attr_name), "Z"); NC_ERR_STOP(ncerr)
        ncerr = nf90_enddef(ncid=ncid); NC_ERR_STOP(ncerr)
        ncerr = nf90_put_var(ncid, varid, values=valueListr4); NC_ERR_STOP(ncerr)
        ncerr = nf90_redef(ncid=ncid); NC_ERR_STOP(ncerr)
        deallocate(valueListr4)
      else if(typekind == ESMF_TYPEKIND_R8) then
        ncerr = nf90_def_var(ncid, trim(dimLabel), NF90_DOUBLE,  dimids=(/dimid/), varid=varid); NC_ERR_STOP(ncerr)
        ncerr = nf90_put_att(ncid, varid, trim(axis_attr_name), "Z"); NC_ERR_STOP(ncerr)
        ncerr = nf90_enddef(ncid=ncid); NC_ERR_STOP(ncerr)
        ncerr = nf90_put_var(ncid, varid, values=valueListr8); NC_ERR_STOP(ncerr)
        ncerr = nf90_redef(ncid=ncid); NC_ERR_STOP(ncerr)
        deallocate(valueListr8)
      endif
    end subroutine write_out_ungridded_dim_atts_from_field

  end subroutine write_restart_netcdf

!----------------------------------------------------------------------------------------
end module module_write_restart_netcdf
