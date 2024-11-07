!> @file
!> @brief Module containing history files output routines.
!> @author Dusan Jovic @date Nov 1, 2017

!> Return error to ESMF and finalize it.
#define ESMF_ERR_RETURN(rc) \
    if (ESMF_LogFoundError(rc, msg="Breaking out of subroutine", line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

!> Return error to ESMF and finalize it.
#define NC_ERR_STOP(status) \
    if (status /= nf90_noerr) write(0,*) "file: ", __FILE__, " line: ", __LINE__, trim(nf90_strerror(status)); \
    if (status /= nf90_noerr) call ESMF_Finalize(endflag=ESMF_END_ABORT)

!> @brief Output routines for writing history files.
!>
!> @author Dusan Jovic @date Nov 1, 2017
module module_write_netcdf

  use mpi_f08
  use esmf
  use netcdf
  use module_fv3_io_def,only : ideflate, quantize_mode, quantize_nsd, zstandard_level, &
                               ichunk2d,jchunk2d,ichunk3d,jchunk3d,kchunk3d, &
                               dx,dy,lon1,lat1,lon2,lat2, &
                               time_unlimited

  implicit none
  private
  public write_netcdf

  logical :: par !< True if parallel I/O should be used.

  integer, parameter :: netcdf_file_type = NF90_NETCDF4 !< NetCDF file type HDF5
  ! integer, parameter :: netcdf_file_type = NF90_64BIT_DATA !< NetCDF file type CDF5
  ! integer, parameter :: netcdf_file_type = NF90_64BIT_OFFSET !< NetCDF file type CDF2

contains

  !> Write netCDF file.
  !>
  !> @param[in] wrtfb ESMF write field bundle.
  !> @param[in] filename NetCDF filename.
  !> @param[in] use_parallel_netcdf True if parallel I/O should be used.
  !> @param[in] comm MPI communicator for parallel I/O.
  !> @param[in] mype MPI rank.
  !> @param[in] grid_id Output grid identifier.
  !> @param[out] rc Return code - 0 for success, ESMF error code otherwise.
  !>
  !> @author Dusan Jovic @date Nov 1, 2017
  subroutine write_netcdf(wrtfb, filename, &
                          use_parallel_netcdf, comm, mype, &
                          grid_id, rc)
!
    use mpi_f08

    type(ESMF_FieldBundle), intent(in) :: wrtfb
    character(*), intent(in)           :: filename
    logical, intent(in)                :: use_parallel_netcdf
    type(MPI_Comm), intent(in)         :: comm
    integer, intent(in)                :: mype
    integer, intent(in)                :: grid_id
    integer, optional,intent(out)      :: rc
!
!** local vars
    integer, parameter :: NF90_NODIMSCALE_ATTACH = int(Z'40000')
    integer :: i,j,t, istart,iend,jstart,jend
    integer :: im, jm, lm, lsoil

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
    type(ESMF_TypeKind_Flag)             :: attTypeKind
    type(ESMF_Grid)                      :: wrtgrid
    type(ESMF_Array)                     :: array
    type(ESMF_DistGrid)                  :: distgrid

    integer :: attCount
    character(len=ESMF_MAXSTR) :: attName, fldName

    integer :: varival
    real(4) :: varr4val
    real(8) :: varr8val
    character(len=ESMF_MAXSTR) :: varcval

    integer :: ncerr,ierr
    integer :: ncid
    integer :: oldMode
    integer :: dim_len
    integer :: im_dimid, jm_dimid, tile_dimid, pfull_dimid, phalf_dimid, time_dimid, lsoil_dimid, ch_dimid
    integer :: im_varid, jm_varid, tile_varid, pfull_varid, phalf_varid, time_varid, lsoil_varid
    integer :: lon_varid, lat_varid, timeiso_varid
    integer, dimension(:), allocatable :: dimids_2d, dimids_3d, dimids_soil, dimids, chunksizes
    integer, dimension(:), allocatable :: varids
    integer :: xtype
    integer :: quant_mode
    integer :: ishuffle
    logical shuffle

    logical :: is_cubed_sphere
    integer :: rank, deCount, localDeCount, dimCount, tileCount
    integer :: my_tile, start_i, start_j
    integer, dimension(:,:), allocatable :: minIndexPDe, maxIndexPDe
    integer, dimension(:,:), allocatable :: minIndexPTile, maxIndexPTile
    integer, dimension(:), allocatable :: deToTileMap, localDeToDeMap
    logical :: do_io
    integer :: par_access
    character(len=ESMF_MAXSTR) :: output_grid_name
!
    interface
      function nf_set_log_level(new_level) result(status)
        integer, intent(in) :: new_level
        integer             :: status
      end function nf_set_log_level
    end interface

    ! ncerr = nf_set_log_level(3); NC_ERR_STOP(ncerr)

    is_cubed_sphere = .false.
    tileCount = 0
    my_tile = 0
    start_i = -10000000
    start_j = -10000000

    par = use_parallel_netcdf

    if (netcdf_file_type /= NF90_NETCDF4) then
       par = .false.
       if (ideflate(grid_id) > 0 .or. zstandard_level(grid_id) > 0) then
          write(0,*)'Compression is unsupporeted in classic netcdf'
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
       end if
    end if

    do_io = par .or. (mype==0)

    call ESMF_FieldBundleGet(wrtfb, fieldCount=fieldCount, rc=rc); ESMF_ERR_RETURN(rc)
    call ESMF_AttributeGet(wrtfb, convention="NetCDF", purpose="FV3", &
                           name='grid', value=output_grid_name, rc=rc); ESMF_ERR_RETURN(rc)

    allocate(fldlev(fieldCount)) ; fldlev = 0
    allocate(fcstField(fieldCount))
    allocate(varids(fieldCount))

    call ESMF_FieldBundleGet(wrtfb, fieldList=fcstField, grid=wrtgrid, &
!                             itemorderflag=ESMF_ITEMORDER_ADDORDER, &
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
          start_i = minIndexPDe(1,localDeToDeMap(1)+1)
          start_j = minIndexPDe(2,localDeToDeMap(1)+1)
          if (.not. par) then
             start_i = 1
             start_j = 1
          end if
          if (is_cubed_sphere) then
             start_i = mod(start_i, im)
             start_j = mod(start_j, jm)
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

    lm = maxval(fldlev(:))
    call get_dimlen_if_exists(ncid, "zsoil", wrtgrid, lsoil, rc)
    if (lsoil > 0 .and. (.not. any(fldlev(:) == lsoil))) then
       lsoil = 0
    end if
    if (lsoil > 0 .and. (.not. any(fldlev(:) > 1 .and. fldlev(:) /= lsoil))) then
       lm = 1
    end if

    ! for serial output allocate 'global' arrays
    if (.not. par) then
       ! allocate(array_r4(im,jm))
       allocate(array_r8(im,jm))
       ! allocate(array_r4_3d(im,jm,lm))
       allocate(array_r8_3d(im,jm,lm))
       if (is_cubed_sphere) then
          ! allocate(array_r4_cube(im,jm,tileCount))
          allocate(array_r8_cube(im,jm,tileCount))
          ! allocate(array_r4_3d_cube(im,jm,lm,tileCount))
          allocate(array_r8_3d_cube(im,jm,lm,tileCount))
       end if
    end if

    ! create netcdf file and enter define mode
    if (do_io) then

       if (par) then
          ncerr = nf90_create(trim(filename),&
                  cmode=IOR(IOR(NF90_CLOBBER,netcdf_file_type),NF90_NODIMSCALE_ATTACH),&
                  comm=comm%mpi_val, info = MPI_INFO_NULL%mpi_val, ncid=ncid); NC_ERR_STOP(ncerr)
       else
          ncerr = nf90_create(trim(filename),&
                  cmode=IOR(IOR(NF90_CLOBBER,netcdf_file_type),NF90_NODIMSCALE_ATTACH),&
                  ncid=ncid); NC_ERR_STOP(ncerr)
       end if

       ! disable auto filling.
       ncerr = nf90_set_fill(ncid, NF90_NOFILL, oldMode); NC_ERR_STOP(ncerr)

       ! define dimensions [grid_xt, grid_yt, nchars, (pfull/phalf), (tile), time]
       ncerr = nf90_def_dim(ncid, "grid_xt", im, im_dimid); NC_ERR_STOP(ncerr)
       ncerr = nf90_def_dim(ncid, "grid_yt", jm, jm_dimid); NC_ERR_STOP(ncerr)
       ncerr = nf90_def_dim(ncid, "nchars", 20, ch_dimid); NC_ERR_STOP(ncerr)
       if (lm > 1) then
         call add_dim(ncid, "pfull", pfull_dimid, pfull_varid, wrtgrid, mype, rc)
         call add_dim(ncid, "phalf", phalf_dimid, phalf_varid, wrtgrid, mype, rc)
       end if
       if (lsoil > 1) then
         call add_dim(ncid, "zsoil", lsoil_dimid, lsoil_varid, wrtgrid, mype, rc)
       end if
       if (is_cubed_sphere) then
          ncerr = nf90_def_dim(ncid, "tile", tileCount, tile_dimid); NC_ERR_STOP(ncerr)
       end if
       call add_dim(ncid, "time", time_dimid, time_varid, wrtgrid, mype, rc)

       ! define coordinate variables
       ncerr = nf90_def_var(ncid, "grid_xt", NF90_DOUBLE, im_dimid, im_varid); NC_ERR_STOP(ncerr)
       ncerr = nf90_put_att(ncid, im_varid, "cartesian_axis", "X"); NC_ERR_STOP(ncerr)
       ncerr = nf90_def_var(ncid, "grid_yt", NF90_DOUBLE, jm_dimid, jm_varid); NC_ERR_STOP(ncerr)
       ncerr = nf90_put_att(ncid, jm_varid, "cartesian_axis", "Y"); NC_ERR_STOP(ncerr)
       if (is_cubed_sphere) then
          ncerr = nf90_def_var(ncid, "tile", NF90_INT, tile_dimid, tile_varid); NC_ERR_STOP(ncerr)
          ncerr = nf90_put_att(ncid, tile_varid, "long_name", "cubed-sphere face"); NC_ERR_STOP(ncerr)
       end if

       ncerr = nf90_def_var(ncid, "time_iso", NF90_CHAR, [ch_dimid,time_dimid], timeiso_varid); NC_ERR_STOP(ncerr)
       ncerr = nf90_put_att(ncid, timeiso_varid, "long_name", "valid time"); NC_ERR_STOP(ncerr)
       ncerr = nf90_put_att(ncid, timeiso_varid, "description", "ISO 8601 datetime string"); NC_ERR_STOP(ncerr)
       ncerr = nf90_put_att(ncid, timeiso_varid, "_Encoding", "UTF-8"); NC_ERR_STOP(ncerr)

       ! coordinate variable attributes based on output_grid type
       if (trim(output_grid_name) == 'gaussian' .or. &
           trim(output_grid_name) == 'latlon') then
          ncerr = nf90_put_att(ncid, im_varid, "long_name", "T-cell longitude"); NC_ERR_STOP(ncerr)
          ncerr = nf90_put_att(ncid, im_varid, "units", "degrees_E"); NC_ERR_STOP(ncerr)
          ncerr = nf90_put_att(ncid, jm_varid, "long_name", "T-cell latiitude"); NC_ERR_STOP(ncerr)
          ncerr = nf90_put_att(ncid, jm_varid, "units", "degrees_N"); NC_ERR_STOP(ncerr)
       else if (trim(output_grid_name) == 'rotated_latlon') then
          ncerr = nf90_put_att(ncid, im_varid, "long_name", "rotated T-cell longiitude"); NC_ERR_STOP(ncerr)
          ncerr = nf90_put_att(ncid, im_varid, "units", "degrees"); NC_ERR_STOP(ncerr)
          ncerr = nf90_put_att(ncid, jm_varid, "long_name", "rotated T-cell latiitude"); NC_ERR_STOP(ncerr)
          ncerr = nf90_put_att(ncid, jm_varid, "units", "degrees"); NC_ERR_STOP(ncerr)
       else if (trim(output_grid_name) == 'lambert_conformal') then
          ncerr = nf90_put_att(ncid, im_varid, "long_name", "x-coordinate of projection"); NC_ERR_STOP(ncerr)
          ncerr = nf90_put_att(ncid, im_varid, "units", "meters"); NC_ERR_STOP(ncerr)
          ncerr = nf90_put_att(ncid, jm_varid, "long_name", "y-coordinate of projection"); NC_ERR_STOP(ncerr)
          ncerr = nf90_put_att(ncid, jm_varid, "units", "meters"); NC_ERR_STOP(ncerr)
       end if

       ! define longitude variable
       if (is_cubed_sphere) then
          ncerr = nf90_def_var(ncid, "lon", NF90_DOUBLE, [im_dimid,jm_dimid,tile_dimid], lon_varid); NC_ERR_STOP(ncerr)
       else
          ncerr = nf90_def_var(ncid, "lon", NF90_DOUBLE, [im_dimid,jm_dimid           ], lon_varid); NC_ERR_STOP(ncerr)
       end if
       ncerr = nf90_put_att(ncid, lon_varid, "long_name", "T-cell longitude"); NC_ERR_STOP(ncerr)
       ncerr = nf90_put_att(ncid, lon_varid, "units", "degrees_E"); NC_ERR_STOP(ncerr)

       ! define latitude variable
       if (is_cubed_sphere) then
          ncerr = nf90_def_var(ncid, "lat", NF90_DOUBLE, [im_dimid,jm_dimid,tile_dimid], lat_varid); NC_ERR_STOP(ncerr)
       else
          ncerr = nf90_def_var(ncid, "lat", NF90_DOUBLE, [im_dimid,jm_dimid           ], lat_varid); NC_ERR_STOP(ncerr)
       end if
       ncerr = nf90_put_att(ncid, lat_varid, "long_name", "T-cell latitude"); NC_ERR_STOP(ncerr)
       ncerr = nf90_put_att(ncid, lat_varid, "units", "degrees_N"); NC_ERR_STOP(ncerr)

       if (par) then
          ncerr = nf90_var_par_access(ncid, im_varid, NF90_COLLECTIVE); NC_ERR_STOP(ncerr)
          ncerr = nf90_var_par_access(ncid, lon_varid, NF90_COLLECTIVE); NC_ERR_STOP(ncerr)
          ncerr = nf90_var_par_access(ncid, jm_varid, NF90_COLLECTIVE); NC_ERR_STOP(ncerr)
          ncerr = nf90_var_par_access(ncid, lat_varid, NF90_COLLECTIVE); NC_ERR_STOP(ncerr)
          ncerr = nf90_var_par_access(ncid, timeiso_varid, NF90_COLLECTIVE); NC_ERR_STOP(ncerr)
          if (is_cubed_sphere) then
             ncerr = nf90_var_par_access(ncid, tile_varid, NF90_COLLECTIVE); NC_ERR_STOP(ncerr)
          end if
       end if

       call get_global_attr(wrtfb, ncid, mype, rc)

       ! define variables (fields)
       if (is_cubed_sphere) then
          allocate(dimids_2d(4))
          allocate(dimids_3d(5), dimids_soil(5))
          dimids_2d =                  [im_dimid,jm_dimid,            tile_dimid,time_dimid]
          if (lm > 1) dimids_3d =      [im_dimid,jm_dimid,pfull_dimid,tile_dimid,time_dimid]
          if (lsoil > 1) dimids_soil = [im_dimid,jm_dimid,lsoil_dimid,tile_dimid,time_dimid]
       else
          allocate(dimids_2d(3))
          allocate(dimids_3d(4), dimids_soil(4))
          dimids_2d =                  [im_dimid,jm_dimid,                       time_dimid]
          if (lm > 1) dimids_3d =      [im_dimid,jm_dimid,pfull_dimid,           time_dimid]
          if (lsoil > 1) dimids_soil = [im_dimid,jm_dimid,lsoil_dimid,           time_dimid]
       end if

       do i=1, fieldCount
         call ESMF_FieldGet(fcstField(i), name=fldName, rank=rank, typekind=typekind, rc=rc); ESMF_ERR_RETURN(rc)

         par_access = NF90_COLLECTIVE

         if (rank == 2) then
           dimids = dimids_2d
         else if (rank == 3) then
           if (fldlev(i) == lsoil) then
             dimids = dimids_soil
           else if (fldlev(i) == lm) then
             dimids = dimids_3d
           endif
         else
           if (mype==0) write(0,*)'Unsupported rank ', rank
           call ESMF_Finalize(endflag=ESMF_END_ABORT)
         end if

         if (typekind == ESMF_TYPEKIND_R4) then
           xtype = NF90_FLOAT
         else if (typekind == ESMF_TYPEKIND_R8) then
           xtype = NF90_DOUBLE
         else
           if (mype==0) write(0,*)'Unsupported typekind ', typekind
           call ESMF_Finalize(endflag=ESMF_END_ABORT)
         end if

         ! define variable
         ncerr = nf90_def_var(ncid, trim(fldName), xtype, dimids, varids(i)) ; NC_ERR_STOP(ncerr)

         ! compression, shuffling  and chunking
         if (ideflate(grid_id) > 0 .or. zstandard_level(grid_id) > 0) then
            par_access = NF90_COLLECTIVE
            if (rank == 2 .and. ichunk2d(grid_id) > 0 .and. jchunk2d(grid_id) > 0) then
               if (is_cubed_sphere) then
                  chunksizes = [im, jm, 1, 1]
               else
                  chunksizes = [ichunk2d(grid_id), jchunk2d(grid_id),            1]
               end if
               ncerr = nf90_def_var_chunking(ncid, varids(i), NF90_CHUNKED, chunksizes) ; NC_ERR_STOP(ncerr)
            else if (rank == 3 .and. ichunk3d(grid_id) > 0 .and. jchunk3d(grid_id) > 0 .and. kchunk3d(grid_id) > 0) then
               if (is_cubed_sphere) then
                  chunksizes = [im, jm, 1, 1, 1]
               else
                  chunksizes = [ichunk3d(grid_id), jchunk3d(grid_id), min(kchunk3d(grid_id),fldlev(i)), 1]
               end if
               ncerr = nf90_def_var_chunking(ncid, varids(i), NF90_CHUNKED, chunksizes) ; NC_ERR_STOP(ncerr)
            end if

            ishuffle = NF90_NOSHUFFLE
            ! shuffle filter on when using lossy compression
            if (quantize_nsd(grid_id) > 0) then
                ishuffle = NF90_SHUFFLE
            end if
            if (ideflate(grid_id) > 0) then
              ncerr = nf90_def_var_deflate(ncid, varids(i), ishuffle, 1, ideflate(grid_id)) ; NC_ERR_STOP(ncerr)
            else if (zstandard_level(grid_id) > 0) then
              ncerr = nf90_def_var_deflate(ncid, varids(i), ishuffle, 0, 0) ; NC_ERR_STOP(ncerr)
              ncerr = nf90_def_var_zstandard(ncid, varids(i), zstandard_level(grid_id)) ; NC_ERR_STOP(ncerr)
            end if

            ! turn on quantize only for 3d variables and if requested
            if (rank == 3 .and. quantize_nsd(grid_id) > 0) then
              ! nf90_quantize_bitgroom = 1
              ! nf90_quantize_granularbr = 2
              ! nf90_quantize_bitround = 3  (nsd is number of bits)
              if (trim(quantize_mode(grid_id)) == 'quantize_bitgroom') then
                quant_mode = 1
              else if (trim(quantize_mode(grid_id)) == 'quantize_granularbr') then
                quant_mode = 2
              else if (trim(quantize_mode(grid_id)) == 'quantize_bitround') then
                quant_mode = 3
              else
                if (mype==0) write(0,*)'Unknown quantize_mode ', trim(quantize_mode(grid_id))
                call ESMF_Finalize(endflag=ESMF_END_ABORT)
              endif

              ncerr = nf90_def_var_quantize(ncid, varids(i), quant_mode, quantize_nsd(grid_id)) ; NC_ERR_STOP(ncerr)
            end if
         end if

         if (par) then
             ncerr = nf90_var_par_access(ncid, varids(i), par_access); NC_ERR_STOP(ncerr)
         end if

         ! define variable attributes
         call ESMF_AttributeGet(fcstField(i), convention="NetCDF", purpose="FV3", &
                                attnestflag=ESMF_ATTNEST_OFF, count=attCount, &
                                rc=rc); ESMF_ERR_RETURN(rc)

         do j=1,attCount
           call ESMF_AttributeGet(fcstField(i), convention="NetCDF", purpose="FV3", &
                                  attnestflag=ESMF_ATTNEST_OFF, attributeIndex=j, &
                                  name=attName, typekind=attTypeKind, &
                                  rc=rc); ESMF_ERR_RETURN(rc)

           if (index(trim(attName),"ESMF") /= 0) then
              cycle
           end if

           if (attTypeKind==ESMF_TYPEKIND_I4) then
              call ESMF_AttributeGet(fcstField(i), convention="NetCDF", purpose="FV3", &
                                     name=trim(attName), value=varival, &
                                     rc=rc); ESMF_ERR_RETURN(rc)
              ncerr = nf90_put_att(ncid, varids(i), trim(attName), varival); NC_ERR_STOP(ncerr)

           else if (attTypeKind==ESMF_TYPEKIND_R4) then
              call ESMF_AttributeGet(fcstField(i), convention="NetCDF", purpose="FV3", &
                                     name=trim(attName), value=varr4val, &
                                     rc=rc); ESMF_ERR_RETURN(rc)
              ncerr = nf90_put_att(ncid, varids(i), trim(attName), varr4val); NC_ERR_STOP(ncerr)

           else if (attTypeKind==ESMF_TYPEKIND_R8) then
              call ESMF_AttributeGet(fcstField(i), convention="NetCDF", purpose="FV3", &
                                     name=trim(attName), value=varr8val, &
                                     rc=rc); ESMF_ERR_RETURN(rc)
              if (trim(attName) /= '_FillValue') then
                 ! FIXME:  _FillValue must be cast to var type when using NF90_NETCDF4
                 ncerr = nf90_put_att(ncid, varids(i), trim(attName), varr8val); NC_ERR_STOP(ncerr)
              end if

           else if (attTypeKind==ESMF_TYPEKIND_CHARACTER) then
              call ESMF_AttributeGet(fcstField(i), convention="NetCDF", purpose="FV3", &
                                     name=trim(attName), value=varcval, &
                                     rc=rc); ESMF_ERR_RETURN(rc)
              ncerr = nf90_put_att(ncid, varids(i), trim(attName), trim(varcval)); NC_ERR_STOP(ncerr)

           end if

         end do ! j=1,attCount

         if (is_cubed_sphere) then
            ncerr = nf90_put_att(ncid, varids(i), 'coordinates', 'lon lat'); NC_ERR_STOP(ncerr)
            ncerr = nf90_put_att(ncid, varids(i), 'grid_mapping', 'cubed_sphere'); NC_ERR_STOP(ncerr)
         end if

       end do   ! i=1,fieldCount

       ncerr = nf90_enddef(ncid); NC_ERR_STOP(ncerr)
       ! end of define mode

       ! write dimension variables, except grid_xt, grid_yt
       ! those will be written later with lon,lat variables
       if (lm > 1) then
         call write_dim(ncid, "pfull", pfull_dimid, pfull_varid, wrtgrid, mype, rc)
         call write_dim(ncid, "phalf", phalf_dimid, phalf_varid, wrtgrid, mype, rc)
       end if
       if (lsoil > 1) then
         call write_dim(ncid, "zsoil", lsoil_dimid, lsoil_varid, wrtgrid, mype, rc)
       end if
       call write_dim(ncid, "time", time_dimid, time_varid, wrtgrid, mype, rc)

    end if

    !
    ! write lon,lat variables
    !
    if (allocated(start_idx)) deallocate(start_idx)
    if (is_cubed_sphere) then
       allocate(start_idx(3))
       start_idx = [start_i, start_j, my_tile]
    else
       allocate(start_idx(2))
       start_idx = [start_i, start_j]
    end if

    ! write lon (lon_varid)
    if (par) then
       call ESMF_GridGetCoord(wrtgrid, coordDim=1, farrayPtr=array_r8, rc=rc); ESMF_ERR_RETURN(rc)
       ncerr = nf90_put_var(ncid, lon_varid, values=array_r8, start=start_idx); NC_ERR_STOP(ncerr)
    else
       call ESMF_GridGetCoord(wrtgrid, coordDim=1, array=array, rc=rc); ESMF_ERR_RETURN(rc)
       if (is_cubed_sphere) then
          do t=1,tileCount
             call ESMF_ArrayGather(array, array_r8_cube(:,:,t), rootPet=0, tile=t, rc=rc); ESMF_ERR_RETURN(rc)
          end do
          if (do_io) then
             ncerr = nf90_put_var(ncid, lon_varid, values=array_r8_cube, start=start_idx); NC_ERR_STOP(ncerr)
          end if
       else
          call ESMF_ArrayGather(array, array_r8, rootPet=0, rc=rc); ESMF_ERR_RETURN(rc)
          if (do_io) then
             ncerr = nf90_put_var(ncid, lon_varid, values=array_r8, start=start_idx); NC_ERR_STOP(ncerr)
          end if
       endif
    end if

    istart = lbound(array_r8,1); iend   = ubound(array_r8,1)
    jstart = lbound(array_r8,2); jend   = ubound(array_r8,2)

    ! write grid_xt (im_varid)
    if (do_io) then
       allocate (x(im))
       if (trim(output_grid_name) == 'gaussian' .or. trim(output_grid_name) == 'latlon') then
          ncerr = nf90_put_var(ncid, im_varid, values=array_r8(:,jstart), start=[istart], count=[iend-istart+1]); NC_ERR_STOP(ncerr)
       else if (trim(output_grid_name) == 'rotated_latlon') then
          do i=1,im
             x(i) = lon1(grid_id) + (lon2(grid_id)-lon1(grid_id))/(im-1) * (i-1)
          end do
          ncerr = nf90_put_var(ncid, im_varid, values=x); NC_ERR_STOP(ncerr)
       else if (trim(output_grid_name) == 'lambert_conformal') then
          do i=1,im
             x(i) = dx(grid_id) * (i-1)
          end do
          ncerr = nf90_put_var(ncid, im_varid, values=x); NC_ERR_STOP(ncerr)
       else if (trim(output_grid_name) == 'cubed_sphere') then
          do i=1,im
             x(i) = i
          end do
          ncerr = nf90_put_var(ncid, im_varid, values=x); NC_ERR_STOP(ncerr)
       else
          if (mype==0) write(0,*)'unknown output_grid ', trim(output_grid_name)
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
       end if
    end if

    ! write lat (lat_varid)
    if (par) then
       call ESMF_GridGetCoord(wrtgrid, coordDim=2, farrayPtr=array_r8, rc=rc); ESMF_ERR_RETURN(rc)
       ncerr = nf90_put_var(ncid, lat_varid, values=array_r8, start=start_idx); NC_ERR_STOP(ncerr)
    else
       call ESMF_GridGetCoord(wrtgrid, coordDim=2, array=array, rc=rc); ESMF_ERR_RETURN(rc)
       if (is_cubed_sphere) then
          do t=1,tileCount
             call ESMF_ArrayGather(array, array_r8_cube(:,:,t), rootPet=0, tile=t, rc=rc); ESMF_ERR_RETURN(rc)
          end do
          if (do_io) then
             ncerr = nf90_put_var(ncid, lat_varid, values=array_r8_cube, start=start_idx); NC_ERR_STOP(ncerr)
          end if
       else
          call ESMF_ArrayGather(array, array_r8, rootPet=0, rc=rc); ESMF_ERR_RETURN(rc)
          if (do_io) then
             ncerr = nf90_put_var(ncid, lat_varid, values=array_r8, start=start_idx); NC_ERR_STOP(ncerr)
          end if
       endif
    end if

    ! write grid_yt (jm_varid)
    if (do_io) then
       allocate (y(jm))
       if (trim(output_grid_name) == 'gaussian' .or. trim(output_grid_name) == 'latlon') then
          ncerr = nf90_put_var(ncid, jm_varid, values=array_r8(istart,:), start=[jstart], count=[jend-jstart+1]); NC_ERR_STOP(ncerr)
       else if (trim(output_grid_name) == 'rotated_latlon') then
          do j=1,jm
             y(j) = lat1(grid_id) + (lat2(grid_id)-lat1(grid_id))/(jm-1) * (j-1)
          end do
          ncerr = nf90_put_var(ncid, jm_varid, values=y); NC_ERR_STOP(ncerr)
       else if (trim(output_grid_name) == 'lambert_conformal') then
          do j=1,jm
             y(j) = dy(grid_id) * (j-1)
          end do
          ncerr = nf90_put_var(ncid, jm_varid, values=y); NC_ERR_STOP(ncerr)
       else if (trim(output_grid_name) == 'cubed_sphere') then
          do j=1,jm
             y(j) = j
          end do
          ncerr = nf90_put_var(ncid, jm_varid, values=y); NC_ERR_STOP(ncerr)
       else
          if (mype==0) write(0,*)'unknown output_grid ', trim(output_grid_name)
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
       end if
    end if

    ! write tile (tile_varid)
    if (do_io .and. is_cubed_sphere) then
       ncerr = nf90_put_var(ncid, tile_varid, values=[1,2,3,4,5,6]); NC_ERR_STOP(ncerr)
    end if

    ! write time_iso (timeiso_varid)
    if (do_io) then
       call ESMF_AttributeGet(wrtgrid, convention="NetCDF", purpose="FV3", &
                              name="time_iso", value=varcval, rc=rc); ESMF_ERR_RETURN(rc)
       ncerr = nf90_put_var(ncid, timeiso_varid, values=[trim(varcval)]); NC_ERR_STOP(ncerr)
    end if

    ! write variables (fields)
    do i=1, fieldCount

       call ESMF_FieldGet(fcstField(i),name=fldName,rank=rank,typekind=typekind, rc=rc); ESMF_ERR_RETURN(rc)

       if (rank == 2) then

         if (allocated(start_idx)) deallocate(start_idx)
         if (is_cubed_sphere) then
            allocate(start_idx(4))
            start_idx = [start_i,start_j,my_tile,1]
         else
            allocate(start_idx(3))
            start_idx = [start_i,start_j,        1]
         end if

         if (typekind == ESMF_TYPEKIND_R4) then
            if (par) then
               call ESMF_FieldGet(fcstField(i), localDe=0, farrayPtr=array_r4, rc=rc); ESMF_ERR_RETURN(rc)
               ncerr = nf90_put_var(ncid, varids(i), values=array_r4, start=start_idx); NC_ERR_STOP(ncerr)
            else
               if (is_cubed_sphere) then
                  call ESMF_FieldGet(fcstField(i), array=array, rc=rc); ESMF_ERR_RETURN(rc)
                  allocate(array_r4_cube(im,jm,tileCount))
                  do t=1,tileCount
                     call ESMF_ArrayGather(array, array_r4_cube(:,:,t), rootPet=0, tile=t, rc=rc); ESMF_ERR_RETURN(rc)
                  end do
                  if (do_io) then
                     ncerr = nf90_put_var(ncid, varids(i), values=array_r4_cube, start=start_idx); NC_ERR_STOP(ncerr)
                  end if
                  deallocate(array_r4_cube)
               else
                  allocate(array_r4(im,jm))
                  call ESMF_FieldGather(fcstField(i), array_r4, rootPet=0, rc=rc); ESMF_ERR_RETURN(rc)
                  if (do_io) then
                     ncerr = nf90_put_var(ncid, varids(i), values=array_r4, start=start_idx); NC_ERR_STOP(ncerr)
                  end if
                  deallocate(array_r4)
               end if
            end if
         else if (typekind == ESMF_TYPEKIND_R8) then
            if (par) then
               call ESMF_FieldGet(fcstField(i), localDe=0, farrayPtr=array_r8, rc=rc); ESMF_ERR_RETURN(rc)
               ncerr = nf90_put_var(ncid, varids(i), values=array_r8, start=start_idx); NC_ERR_STOP(ncerr)
            else
               if (is_cubed_sphere) then
                  call ESMF_FieldGet(fcstField(i), array=array, rc=rc); ESMF_ERR_RETURN(rc)
                  do t=1,tileCount
                     call ESMF_ArrayGather(array, array_r8_cube(:,:,t), rootPet=0, tile=t, rc=rc); ESMF_ERR_RETURN(rc)
                  end do
                  if (do_io) then
                     ncerr = nf90_put_var(ncid, varids(i), values=array_r8_cube, start=start_idx); NC_ERR_STOP(ncerr)
                  end if
               else
                  call ESMF_FieldGather(fcstField(i), array_r8, rootPet=0, rc=rc); ESMF_ERR_RETURN(rc)
                  if (do_io) then
                     ncerr = nf90_put_var(ncid, varids(i), values=array_r8, start=start_idx); NC_ERR_STOP(ncerr)
                  end if
               end if
            end if
         end if

      else if (rank == 3) then

         if (allocated(start_idx)) deallocate(start_idx)
         if (is_cubed_sphere) then
            allocate(start_idx(5))
            start_idx = [start_i,start_j,1,my_tile,1]
         else
            allocate(start_idx(4))
            start_idx = [start_i,start_j,1,        1]
         end if

         if (typekind == ESMF_TYPEKIND_R4) then
            if (par) then
               call ESMF_FieldGet(fcstField(i), localDe=0, farrayPtr=array_r4_3d, rc=rc); ESMF_ERR_RETURN(rc)
               ncerr = nf90_put_var(ncid, varids(i), values=array_r4_3d, start=start_idx); NC_ERR_STOP(ncerr)
            else
               if (is_cubed_sphere) then
                  call ESMF_FieldGet(fcstField(i), array=array, rc=rc); ESMF_ERR_RETURN(rc)
                  allocate(array_r4_3d_cube(im,jm,fldlev(i),tileCount))
                  do t=1,tileCount
                     call ESMF_ArrayGather(array, array_r4_3d_cube(:,:,:,t), rootPet=0, tile=t, rc=rc); ESMF_ERR_RETURN(rc)
                  end do
                  if (mype==0) then
                     ncerr = nf90_put_var(ncid, varids(i), values=array_r4_3d_cube, start=start_idx); NC_ERR_STOP(ncerr)
                  end if
                  deallocate(array_r4_3d_cube)
               else
                  allocate(array_r4_3d(im,jm,fldlev(i)))
                  call ESMF_FieldGather(fcstField(i), array_r4_3d, rootPet=0, rc=rc); ESMF_ERR_RETURN(rc)
                  if (mype==0) then
                     ncerr = nf90_put_var(ncid, varids(i), values=array_r4_3d, start=start_idx); NC_ERR_STOP(ncerr)
                  end if
                  deallocate(array_r4_3d)
               end if
            end if
         else if (typekind == ESMF_TYPEKIND_R8) then
            if (par) then
               call ESMF_FieldGet(fcstField(i), localDe=0, farrayPtr=array_r8_3d, rc=rc); ESMF_ERR_RETURN(rc)
               ncerr = nf90_put_var(ncid, varids(i), values=array_r8_3d, start=start_idx); NC_ERR_STOP(ncerr)
            else
               if (is_cubed_sphere) then
                  call ESMF_FieldGet(fcstField(i), array=array, rc=rc); ESMF_ERR_RETURN(rc)
                  do t=1,tileCount
                     call ESMF_ArrayGather(array, array_r8_3d_cube(:,:,:,t), rootPet=0, tile=t, rc=rc); ESMF_ERR_RETURN(rc)
                  end do
                  if (mype==0) then
                     ncerr = nf90_put_var(ncid, varids(i), values=array_r8_3d_cube, start=start_idx); NC_ERR_STOP(ncerr)
                  end if
               else
                  call ESMF_FieldGather(fcstField(i), array_r8_3d, rootPet=0, rc=rc); ESMF_ERR_RETURN(rc)
                  if (mype==0) then
                     ncerr = nf90_put_var(ncid, varids(i), values=array_r8_3d, start=start_idx); NC_ERR_STOP(ncerr)
                  end if
               end if
            end if
         end if ! end typekind

      else

         if (mype==0) write(0,*)'Unsupported rank ', rank
         call ESMF_Finalize(endflag=ESMF_END_ABORT)

      end if ! end rank

    end do ! end fieldCount

    if (.not. par) then
       ! deallocate(array_r4)
       deallocate(array_r8)
       ! deallocate(array_r4_3d)
       deallocate(array_r8_3d)
       if (is_cubed_sphere) then
          ! deallocate(array_r4_cube)
          deallocate(array_r8_cube)
          ! deallocate(array_r4_3d_cube)
          deallocate(array_r8_3d_cube)
       end if
    end if

    if (do_io) then
       deallocate(dimids_2d)
       deallocate(dimids_3d)
       deallocate(dimids_soil)
    end if

    deallocate(fcstField)
    deallocate(varids)

    if (do_io) then
       ncerr = nf90_close(ncid=ncid); NC_ERR_STOP(ncerr)
    end if

  end subroutine write_netcdf

  !> Get global attribute.
  !>
  !> @param[in] fldbundle ESMF field bundle.
  !> @param[in] ncid NetCDF file ID.
  !> @param[in] mype MPI rank.
  !> @param[out] rc Return code - 0 for success, ESMF error code otherwise.
  !>
  !> @author Dusan Jovic @date Nov 1, 2017
  subroutine get_global_attr(fldbundle, ncid, mype, rc)
    type(ESMF_FieldBundle), intent(in) :: fldbundle
    integer, intent(in)                :: ncid
    integer, intent(in)                :: mype
    integer, intent(out)               :: rc

! local variable
    integer :: i, attCount
    integer :: ncerr
    character(len=ESMF_MAXSTR) :: attName
    type(ESMF_TypeKind_Flag)   :: typekind

    integer(ESMF_KIND_I4) :: varival_i4
    integer(ESMF_KIND_I8) :: varival_i8
    real(ESMF_KIND_R4), dimension(:), allocatable :: varr4list
    real(ESMF_KIND_R8), dimension(:), allocatable :: varr8list
    integer :: itemCount
    character(len=ESMF_MAXSTR) :: varcval
!
    call ESMF_AttributeGet(fldbundle, convention="NetCDF", purpose="FV3", &
                           attnestflag=ESMF_ATTNEST_OFF, count=attCount, &
                           rc=rc); ESMF_ERR_RETURN(rc)

    do i=1,attCount

      call ESMF_AttributeGet(fldbundle, convention="NetCDF", purpose="FV3", &
                             attnestflag=ESMF_ATTNEST_OFF, attributeIndex=i, name=attName, &
                             typekind=typekind, itemCount=itemCount, rc=rc); ESMF_ERR_RETURN(rc)

      if (typekind==ESMF_TYPEKIND_I4) then
         call ESMF_AttributeGet(fldbundle, convention="NetCDF", purpose="FV3", &
                                name=trim(attname), value=varival_i4, rc=rc); ESMF_ERR_RETURN(rc)
         ncerr = nf90_put_att(ncid, nf90_global, trim(attname), varival_i4); NC_ERR_STOP(ncerr)

      else if (typekind==ESMF_TYPEKIND_I8) then
         call ESMF_AttributeGet(fldbundle, convention="NetCDF", purpose="FV3", &
                                name=trim(attname), value=varival_i8, rc=rc); ESMF_ERR_RETURN(rc)
         if (netcdf_file_type == NF90_64BIT_OFFSET) then
            ! NetCDF NF90_64BIT_OFFSET (CDF2) does not support int64 attributes
            ! Currently only one global attribute is int64 (:grid_id = 1LL)
            varival_i4 = varival_i8
            ncerr = nf90_put_att(ncid, nf90_global, trim(attname), varival_i4); NC_ERR_STOP(ncerr)
         else
            ncerr = nf90_put_att(ncid, nf90_global, trim(attname), varival_i8); NC_ERR_STOP(ncerr)
         end if

      else if (typekind==ESMF_TYPEKIND_R4) then
         allocate (varr4list(itemCount))
         call ESMF_AttributeGet(fldbundle, convention="NetCDF", purpose="FV3", &
                                name=trim(attName), valueList=varr4list, rc=rc); ESMF_ERR_RETURN(rc)
         ncerr = nf90_put_att(ncid, NF90_GLOBAL, trim(attName), varr4list); NC_ERR_STOP(ncerr)
         deallocate(varr4list)

      else if (typekind==ESMF_TYPEKIND_R8) then
         allocate (varr8list(itemCount))
         call ESMF_AttributeGet(fldbundle, convention="NetCDF", purpose="FV3", &
                                name=trim(attName), valueList=varr8list, rc=rc); ESMF_ERR_RETURN(rc)
         ncerr = nf90_put_att(ncid, NF90_GLOBAL, trim(attName), varr8list); NC_ERR_STOP(ncerr)
         deallocate(varr8list)

      else if (typekind==ESMF_TYPEKIND_CHARACTER) then
         call ESMF_AttributeGet(fldbundle, convention="NetCDF", purpose="FV3", &
                                name=trim(attName), value=varcval, rc=rc); ESMF_ERR_RETURN(rc)
         ncerr = nf90_put_att(ncid, NF90_GLOBAL, trim(attName), trim(varcval)); NC_ERR_STOP(ncerr)

      else

         if (mype==0) write(0,*)'Unsupported typekind ', typekind
         call ESMF_Finalize(endflag=ESMF_END_ABORT)
      end if

    end do

  end subroutine get_global_attr

  !> Get grid attribute.
  !>
  !> @param[in] grid ESMF output grid.
  !> @param[in] prefix grid attribute prefix.
  !> @param[in] ncid NetCDF file ID.
  !> @param[in] varid NetCDF variable ID.
  !> @param[out] rc Return code - 0 for success, ESMF error code otherwise.
  !>
  !> @author Dusan Jovic @date Nov 1, 2017
  subroutine get_grid_attr(grid, prefix, ncid, varid, rc)
    type(ESMF_Grid), intent(in)  :: grid
    character(len=*), intent(in) :: prefix
    integer, intent(in)          :: ncid
    integer, intent(in)          :: varid
    integer, intent(out)         :: rc

! local variable
    integer :: i, attCount, n, ind
    integer :: ncerr
    character(len=ESMF_MAXSTR) :: attName
    type(ESMF_TypeKind_Flag)   :: typekind

    integer :: varival
    real(ESMF_KIND_R4) :: varr4val
    real(ESMF_KIND_R8) :: varr8val
    character(len=ESMF_MAXSTR) :: varcval
!
    call ESMF_AttributeGet(grid, convention="NetCDF", purpose="FV3", &
                           attnestflag=ESMF_ATTNEST_OFF, count=attCount, &
                           rc=rc); ESMF_ERR_RETURN(rc)

    do i=1,attCount

      call ESMF_AttributeGet(grid, convention="NetCDF", purpose="FV3", &
                             attnestflag=ESMF_ATTNEST_OFF, attributeIndex=i, name=attName, &
                             typekind=typekind, itemCount=n, rc=rc); ESMF_ERR_RETURN(rc)

      if (index(trim(attName), trim(prefix)//":")==1) then
         ind = len(trim(prefix)//":")

         if (typekind==ESMF_TYPEKIND_I4) then
            call ESMF_AttributeGet(grid, convention="NetCDF", purpose="FV3", &
                                   name=trim(attName), value=varival, rc=rc); ESMF_ERR_RETURN(rc)
            ncerr = nf90_put_att(ncid, varid, trim(attName(ind+1:len(attName))), varival); NC_ERR_STOP(ncerr)

         else if (typekind==ESMF_TYPEKIND_R4) then
            call ESMF_AttributeGet(grid, convention="NetCDF", purpose="FV3", &
                                   name=trim(attName), value=varr4val, rc=rc); ESMF_ERR_RETURN(rc)
            ncerr = nf90_put_att(ncid, varid, trim(attName(ind+1:len(attName))), varr4val); NC_ERR_STOP(ncerr)

         else if (typekind==ESMF_TYPEKIND_R8) then
            call ESMF_AttributeGet(grid, convention="NetCDF", purpose="FV3", &
                                   name=trim(attName), value=varr8val, rc=rc); ESMF_ERR_RETURN(rc)
            if (trim(attName) /= '_FillValue') then
              ! FIXME:  _FillValue must be cast to var type when using
              ! NF90_NETCDF4. Until this is fixed, using netCDF default _FillValue.
              ncerr = nf90_put_att(ncid, varid, trim(attName(ind+1:len(attName))), varr8val); NC_ERR_STOP(ncerr)
            end if

         else if (typekind==ESMF_TYPEKIND_CHARACTER) then
            call ESMF_AttributeGet(grid, convention="NetCDF", purpose="FV3", &
                                   name=trim(attName), value=varcval, rc=rc); ESMF_ERR_RETURN(rc)
            ncerr = nf90_put_att(ncid, varid, trim(attName(ind+1:len(attName))), trim(varcval)); NC_ERR_STOP(ncerr)

         end if

      end if

    end do

  end subroutine get_grid_attr

!----------------------------------------------------------------------------------------
  subroutine get_dimlen_if_exists(ncid, dim_name, grid, dim_len, rc)

    integer, intent(in)             :: ncid
    character(len=*), intent(in)    :: dim_name
    type(ESMF_Grid), intent(in)     :: grid
    integer, intent(out)            :: dim_len
    integer, intent(out)            :: rc

    dim_len = 0

    call ESMF_AttributeGet(grid, convention="NetCDF", purpose="FV3", &
                           name=dim_name, itemCount=dim_len, rc=rc); ESMF_ERR_RETURN(rc)

  end subroutine get_dimlen_if_exists

  !> Add a dimension.
  !>
  !> @param[in] ncid NetCDF file ID.
  !> @param[in] dim_name Dimension name.
  !> @param[in] dimid Dimension ID.
  !> @param[in] grpid Group ID.
  !> @param[in] grid ESMF output grid.
  !> @param[in] mype MPI rank.
  !> @param[out] rc Return code - 0 for success, ESMF error code otherwise.
  !>
  !> @author Dusan Jovic @date Nov 1, 2017
  subroutine add_dim(ncid, dim_name, dimid, dim_varid, grid, mype, rc)
    integer, intent(in)             :: ncid
    character(len=*), intent(in)    :: dim_name
    integer, intent(inout)          :: dimid
    integer, intent(inout)          :: dim_varid
    type(ESMF_Grid), intent(in)     :: grid
    integer, intent(in)             :: mype
    integer, intent(out)            :: rc

! local variable
    integer :: n
    integer :: ncerr
    type(ESMF_TypeKind_Flag)   :: typekind
!
    call ESMF_AttributeGet(grid, convention="NetCDF", purpose="FV3", &
                           attnestflag=ESMF_ATTNEST_OFF, name=dim_name, &
                           typekind=typekind, itemCount=n, rc=rc); ESMF_ERR_RETURN(rc)

    if (trim(dim_name) == "time") then
      ! using an unlimited dim requires collective mode (NF90_COLLECTIVE)
      ! for parallel writes, which seems to slow things down on hera.
      if (time_unlimited) then
        ncerr = nf90_def_dim(ncid, trim(dim_name), NF90_UNLIMITED, dimid); NC_ERR_STOP(ncerr)
      else
        ncerr = nf90_def_dim(ncid, trim(dim_name), 1, dimid); NC_ERR_STOP(ncerr)
      end if
    else
      ncerr = nf90_def_dim(ncid, trim(dim_name), n, dimid); NC_ERR_STOP(ncerr)
    end if

    if (typekind==ESMF_TYPEKIND_R8) then
      ncerr = nf90_def_var(ncid, dim_name, NF90_REAL8, dimids=[dimid], varid=dim_varid); NC_ERR_STOP(ncerr)
    else if (typekind==ESMF_TYPEKIND_R4) then
      ncerr = nf90_def_var(ncid, dim_name, NF90_REAL4, dimids=[dimid], varid=dim_varid); NC_ERR_STOP(ncerr)
    else if (typekind==ESMF_TYPEKIND_I4) then
      ncerr = nf90_def_var(ncid, dim_name, NF90_INT4, dimids=[dimid], varid=dim_varid); NC_ERR_STOP(ncerr)
    else
      if (mype==0) write(0,*)'Error in module_write_netcdf.F90(add_dim) unknown typekind for ',trim(dim_name)
      call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if
    if (par) then
      ncerr = nf90_var_par_access(ncid, dim_varid, NF90_COLLECTIVE); NC_ERR_STOP(ncerr)
    end if

    call get_grid_attr(grid, dim_name, ncid, dim_varid, rc)

  end subroutine add_dim

  !> Write a dimension variable.
  !>
  !> @param[in] ncid NetCDF file ID.
  !> @param[in] dim_name Dimension name.
  !> @param[in] dimid Dimension ID.
  !> @param[in] dim_varid Dimension variable ID.
  !> @param[in] grid ESMF output grid.
  !> @param[in] mype MPI rank.
  !> @param[out] rc Return code - 0 for success, ESMF error code otherwise.
  !>
  !> @author Dusan Jovic @date Nov 1, 2017
  subroutine write_dim(ncid, dim_name, dimid, dim_varid, grid, mype, rc)
    integer, intent(in)             :: ncid
    character(len=*), intent(in)    :: dim_name
    integer, intent(in)             :: dimid
    integer, intent(in)             :: dim_varid
    type(ESMF_Grid), intent(in)     :: grid
    integer, intent(in)             :: mype
    integer, intent(out)            :: rc

! local variable
    integer :: n
    integer :: ncerr
    type(ESMF_TypeKind_Flag)   :: typekind

    real(ESMF_KIND_I4), allocatable  :: valueListI4(:)
    real(ESMF_KIND_R4), allocatable  :: valueListR4(:)
    real(ESMF_KIND_R8), allocatable  :: valueListR8(:)
!
    call ESMF_AttributeGet(grid, convention="NetCDF", purpose="FV3", &
                           attnestflag=ESMF_ATTNEST_OFF, name=dim_name, &
                           typekind=typekind, itemCount=n, rc=rc); ESMF_ERR_RETURN(rc)

    if (typekind==ESMF_TYPEKIND_R8) then
       allocate(valueListR8(n))
       call ESMF_AttributeGet(grid, convention="NetCDF", purpose="FV3", &
                              name=trim(dim_name), valueList=valueListR8, rc=rc); ESMF_ERR_RETURN(rc)
       ncerr = nf90_put_var(ncid, dim_varid, values=valueListR8); NC_ERR_STOP(ncerr)
       deallocate(valueListR8)
    else if (typekind==ESMF_TYPEKIND_R4) then
       allocate(valueListR4(n))
       call ESMF_AttributeGet(grid, convention="NetCDF", purpose="FV3", &
                              name=trim(dim_name), valueList=valueListR4, rc=rc); ESMF_ERR_RETURN(rc)
       ncerr = nf90_put_var(ncid, dim_varid, values=valueListR4); NC_ERR_STOP(ncerr)
       deallocate(valueListR4)
    else if (typekind==ESMF_TYPEKIND_I4) then
       allocate(valueListI4(n))
       call ESMF_AttributeGet(grid, convention="NetCDF", purpose="FV3", &
                              name=trim(dim_name), valueList=valueListI4, rc=rc); ESMF_ERR_RETURN(rc)
       ncerr = nf90_put_var(ncid, dim_varid, values=valueListI4); NC_ERR_STOP(ncerr)
       deallocate(valueListI4)
    else
       if (mype==0) write(0,*)'Error in module_write_netcdf.F90(write_dim) unknown typekind for ',trim(dim_name)
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if

  end subroutine write_dim

!----------------------------------------------------------------------------------------
end module module_write_netcdf
