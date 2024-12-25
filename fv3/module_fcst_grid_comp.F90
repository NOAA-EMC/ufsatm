#define ESMF_ERR_ABORT(rc) \
if (rc /= ESMF_SUCCESS) write(0,*) 'rc=',rc,__FILE__,__LINE__; if(ESMF_LogFoundError(rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

!-----------------------------------------------------------------------
!
  module module_fcst_grid_comp
!
!-----------------------------------------------------------------------
!***  Forecast gridded component.
!-----------------------------------------------------------------------
!***
!***  HISTORY
!***
!       Apr 2017:  J. Wang  - initial code for forecast grid component
!
!---------------------------------------------------------------------------------
!
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

  use  atmos_model_mod,   only: atmos_model_init, atmos_model_end,         &
                                get_atmos_model_ungridded_dim,             &
                                update_atmos_model_dynamics,               &
                                update_atmos_radiation_physics,            &
                                update_atmos_model_state,                  &
                                atmos_data_type, atmos_model_restart,      &
                                atmos_model_exchange_phase_1,              &
                                atmos_model_exchange_phase_2,              &
                                addLsmask2grid, atmos_model_get_nth_domain_info

  use GFS_typedefs,       only: kind_phys, kind_sngl_prec

  use constants_mod,      only: constants_init
  use fms_mod,            only: error_mesg, fms_init, fms_end,             &
                                write_version_number, uppercase

  use mpp_mod,            only: mpp_init, mpp_pe, mpp_npes, mpp_root_pe, mpp_set_current_pelist,  &
                                mpp_error, FATAL, WARNING, NOTE
  use mpp_mod,            only: mpp_clock_id, mpp_clock_begin
  use mpp_domains_mod,    only: mpp_get_compute_domains, domain2D
  use sat_vapor_pres_mod, only: sat_vapor_pres_init

  use diag_manager_mod,   only: diag_manager_init, diag_manager_end,       &
                                diag_manager_set_time_end

  use data_override_mod,  only: data_override_init
  use fv_nggps_diags_mod, only: fv_dyn_bundle_setup
  use fv3atm_history_io_mod,  only: fv_phys_bundle_setup
  use fv3atm_restart_io_mod,  only: fv_phy_restart_bundle_setup, fv_sfc_restart_bundle_setup
  use fv_ufs_restart_io_mod,  only: fv_core_restart_bundle_setup, &
                                    fv_srf_wnd_restart_bundle_setup, &
                                    fv_tracer_restart_bundle_setup

  use fms2_io_mod,        only: FmsNetcdfFile_t, open_file, close_file, variable_exists, read_data

  use atmosphere_mod,     only: atmosphere_control_data

  use module_fv3_io_def,  only: num_pes_fcst, num_files, filename_base,    &
                                nbdlphys, iau_offset
  use module_fv3_config,  only: dt_atmos, fcst_mpi_comm, fcst_ntasks,      &
                                quilting, quilting_restart,                &
                                calendar, cpl_grid_id,                     &
                                cplprint_flag

  use get_stochy_pattern_mod, only: write_stoch_restart_atm
  use module_cplfields,       only: nExportFields, exportFields, exportFieldsInfo, &
                                    nImportFields, importFields, importFieldsInfo
  use module_cplfields,       only: realizeConnectedCplFields

  use atmos_model_mod,        only: setup_exportdata
  use CCPP_data,              only: GFS_control
!
!-----------------------------------------------------------------------
!
  implicit none
!
!-----------------------------------------------------------------------
!
  private
!
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
!
!-----------------------------------------------------------------------
!
  public SetServices
!
  contains
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
  subroutine SetServices(fcst_comp, rc)
!
    type(ESMF_GridComp)  :: fcst_comp
    integer, intent(out) :: rc

    rc = ESMF_SUCCESS

    call ESMF_GridCompSetEntryPoint(fcst_comp, ESMF_METHOD_INITIALIZE, &
                                    userRoutine=fcst_initialize, phase=1, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_GridCompSetEntryPoint(fcst_comp, ESMF_METHOD_INITIALIZE, &
                                    userRoutine=fcst_advertise, phase=2, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_GridCompSetEntryPoint(fcst_comp, ESMF_METHOD_INITIALIZE, &
                                    userRoutine=fcst_realize, phase=3, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
!
    call ESMF_GridCompSetEntryPoint(fcst_comp, ESMF_METHOD_RUN, &
                                    userRoutine=fcst_run_phase_1, phase=1, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_GridCompSetEntryPoint(fcst_comp, ESMF_METHOD_RUN, &
                                    userRoutine=fcst_run_phase_2, phase=2, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
!
    call ESMF_GridCompSetEntryPoint(fcst_comp, ESMF_METHOD_FINALIZE, &
                                    userRoutine=fcst_finalize, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

  end subroutine SetServices
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
  subroutine SetServicesNest(nest, rc)
!
    type(ESMF_GridComp)   :: nest
    integer, intent(out)  :: rc

    character(len=80)     :: name
    type(ESMF_Grid)       :: grid
    type(ESMF_Info)       :: info
    integer               :: layout(2), tilesize
    integer               :: tl, nx, ny
    integer,dimension(2,6):: decomptile                  !define delayout for the 6 cubed-sphere tiles
    integer,dimension(2)  :: regdecomp                   !define delayout for the nest grid
    type(ESMF_Decomp_Flag):: decompflagPTile(2,6)
    type(ESMF_TypeKind_Flag) :: grid_typekind
    character(3)          :: myGridStr
    type(ESMF_DistGrid)   :: distgrid
    type(ESMF_Array)      :: array

    rc = ESMF_SUCCESS

    call ESMF_GridCompSetEntryPoint(nest, ESMF_METHOD_INITIALIZE, userRoutine=init_dyn_fb, phase=1, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_GridCompSetEntryPoint(nest, ESMF_METHOD_INITIALIZE, userRoutine=init_phys_fb, phase=2, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_GridCompSetEntryPoint(nest, ESMF_METHOD_INITIALIZE, userRoutine=init_advertise, phase=3, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_GridCompSetEntryPoint(nest, ESMF_METHOD_INITIALIZE, userRoutine=init_realize, phase=4, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_GridCompGet(nest, name=name, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_InfoGetFromHost(nest, info=info, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_InfoGet(info, key="layout", values=layout, rc=rc); ESMF_ERR_ABORT(rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    if (kind_phys == kind_sngl_prec) then
      grid_typekind = ESMF_TYPEKIND_R4
    else
      grid_typekind = ESMF_TYPEKIND_R8
    endif

    if (trim(name)=="global") then
      ! global domain
      call ESMF_InfoGet(info, key="tilesize", value=tilesize, rc=rc); ESMF_ERR_ABORT(rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      do tl=1,6
        decomptile(1,tl) = layout(1)
        decomptile(2,tl) = layout(2)
        decompflagPTile(:,tl) = (/ESMF_DECOMP_SYMMEDGEMAX,ESMF_DECOMP_SYMMEDGEMAX/)
      enddo
      grid = ESMF_GridCreateCubedSphere(tileSize=tilesize, &
                                        coordSys=ESMF_COORDSYS_SPH_RAD, &
                                        coordTypeKind=grid_typekind, &
                                        regDecompPTile=decomptile, &
                                        decompflagPTile=decompflagPTile, &
                                        name="fcst_grid", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    else
      ! nest domain
      call ESMF_InfoGet(info, key="nx", value=nx, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      call ESMF_InfoGet(info, key="ny", value=ny, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      grid = ESMF_GridCreateNoPeriDim(regDecomp=(/layout(1),layout(2)/), &
                                      minIndex=(/1,1/), &
                                      maxIndex=(/nx,ny/), &
                                      gridAlign=(/-1,-1/), &
                                      coordSys=ESMF_COORDSYS_SPH_RAD, &
                                      coordTypeKind=grid_typekind, &
                                      decompflag=(/ESMF_DECOMP_SYMMEDGEMAX,ESMF_DECOMP_SYMMEDGEMAX/), &
                                      name="fcst_grid", &
                                      indexflag=ESMF_INDEX_DELOCAL, &
                                      rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    endif

    ! - Create coordinate arrays around allocations held within Atmos data structure and set in Grid

    call ESMF_GridGet(grid, staggerloc=ESMF_STAGGERLOC_CENTER, distgrid=distgrid, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    array = ESMF_ArrayCreate(distgrid, farray=Atmos%lon, indexflag=ESMF_INDEX_DELOCAL, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_GridSetCoord(grid, coordDim=1, staggerLoc=ESMF_STAGGERLOC_CENTER, array=array, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    array = ESMF_ArrayCreate(distgrid, farray=Atmos%lat, indexflag=ESMF_INDEX_DELOCAL, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_GridSetCoord(grid, coordDim=2, staggerLoc=ESMF_STAGGERLOC_CENTER, array=array, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_GridGet(grid, staggerloc=ESMF_STAGGERLOC_CORNER, distgrid=distgrid, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    array = ESMF_ArrayCreate(distgrid, farray=Atmos%lon_bnd, indexflag=ESMF_INDEX_DELOCAL, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_GridSetCoord(grid, coordDim=1, staggerLoc=ESMF_STAGGERLOC_CORNER, array=array, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    array = ESMF_ArrayCreate(distgrid, farray=Atmos%lat_bnd, indexflag=ESMF_INDEX_DELOCAL, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_GridSetCoord(grid, coordDim=2, staggerLoc=ESMF_STAGGERLOC_CORNER, array=array, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    !TODO: Consider aligning mask treatment with coordinates... especially if it requires updates for moving
    call addLsmask2grid(grid, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! - Add Attributes used by output

    call ESMF_AttributeAdd(grid, convention="NetCDF", purpose="FV3", &
                          attrList=(/"ESMF:gridded_dim_labels"/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_AttributeSet(grid, convention="NetCDF", purpose="FV3", &
                         name="ESMF:gridded_dim_labels", valueList=(/"grid_xt", "grid_yt"/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

!test to write out vtk file:
!    if( cplprint_flag ) then
!      call ESMF_GridWriteVTK(grid, staggerloc=ESMF_STAGGERLOC_CENTER,  &
!                             filename='fv3cap_fv3Grid', rc=rc)
!      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
!    endif
!
! Write grid to netcdf file
    if( cplprint_flag ) then
      write (myGridStr,"(I0)") mygrid
      if (trim(name)=="global") then
        call wrt_fcst_grid(grid, "diagnostic_FV3_fcstGrid"//trim(mygridStr)//".tile*.nc", &
                           regridArea=.TRUE., rc=rc)
      else
        call wrt_fcst_grid(grid, "diagnostic_FV3_fcstGrid"//trim(mygridStr)//".nc", &
                           regridArea=.TRUE., rc=rc)
      end if
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    endif

    ! - Hold on to the grid by GridComp

    call ESMF_GridCompSet(nest, grid=grid, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

  end subroutine SetServicesNest
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
  subroutine init_dyn_fb(nest, importState, exportState, clock, rc)
!
    type(ESMF_GridComp)                    :: nest
    type(ESMF_State)                       :: importState, exportState
    type(ESMF_Clock)                       :: clock
    integer,intent(out)                    :: rc

    type(ESMF_Grid)                        :: grid
    integer                                :: itemCount, i
    character(len=ESMF_MAXSTR), allocatable :: itemNameList(:)
    character(len=ESMF_MAXSTR)              :: fb_name
    type(ESMF_FieldBundle), allocatable     :: fbList(:)
    type(ESMF_FieldBundle)                  :: fcstFB

    call ESMF_GridCompGet(nest, grid=grid, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_StateGet(importState, itemCount=itemCount, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    allocate(itemNameList(itemCount), fbList(itemCount))

    call ESMF_StateGet(importState, itemNameList=itemNameList, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    do i=1, itemCount
      call ESMF_StateGet(importState, itemName=itemNameList(i), fieldbundle=fcstFB, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      fbList(i) = ESMF_FieldBundleCreate(name=itemNameList(i), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call ESMF_AttributeCopy(fcstFB, fbList(i), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call ESMF_StateAdd(exportState, (/fbList(i)/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    enddo

    ! get the name of the first field bundle and based on that determine if it's a history or restart bundles
    call ESMF_FieldBundleGet(fbList(1), name=fb_name, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    if (fb_name(1:19) == 'restart_fv_core.res') then
      call fv_core_restart_bundle_setup(fbList(1), grid, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    else if (fb_name(1:22) == 'restart_fv_srf_wnd.res') then
      call fv_srf_wnd_restart_bundle_setup(fbList(1), grid, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    else if (fb_name(1:21) == 'restart_fv_tracer.res') then
      call fv_tracer_restart_bundle_setup(fbList(1), grid, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    else
      do i=1, itemCount
        call fv_dyn_bundle_setup(Atmos%axes, fbList(i), grid, quilting=.true., rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      end do
    endif

  end subroutine init_dyn_fb
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
  subroutine init_phys_fb(nest, importState, exportState, clock, rc)
!
    type(ESMF_GridComp)                    :: nest
    type(ESMF_State)                       :: importState, exportState
    type(ESMF_Clock)                       :: clock
    integer,intent(out)                    :: rc

    type(ESMF_Grid)                        :: grid
    integer                                :: itemCount, i
    character(len=ESMF_MAXSTR), allocatable :: itemNameList(:)
    character(len=ESMF_MAXSTR)              :: fb_name
    type(ESMF_FieldBundle), allocatable     :: fbList(:)
    type(ESMF_FieldBundle)                  :: fcstFB

    call ESMF_GridCompGet(nest, grid=grid, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_StateGet(importState, itemCount=itemCount, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    allocate(itemNameList(itemCount), fbList(itemCount))

    call ESMF_StateGet(importState, itemNameList=itemNameList, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    do i=1, itemCount
      call ESMF_StateGet(importState, itemName=itemNameList(i), fieldbundle=fcstFB, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      fbList(i) = ESMF_FieldBundleCreate(name=itemNameList(i), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call ESMF_AttributeCopy(fcstFB, fbList(i), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call ESMF_StateAdd(exportState, (/fbList(i)/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    enddo

    ! get the name of the first field bundle and based on that determine if it's a history or restart bundles
    call ESMF_FieldBundleGet(fbList(1), name=fb_name, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    if (fb_name(1:16) == 'restart_phy_data') then
      call fv_phy_restart_bundle_setup(fbList(1), grid, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    elseif (fb_name(1:16) == 'restart_sfc_data') then
      call fv_sfc_restart_bundle_setup(fbList(1), grid, GFS_control, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    else
      call fv_phys_bundle_setup(Atmos%diag, Atmos%axes, fbList, grid, quilting=.true., nbdlphys=itemCount, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    endif

  end subroutine init_phys_fb
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
  subroutine init_advertise(nest, importState, exportState, clock, rc)
!
    type(ESMF_GridComp)                    :: nest
    type(ESMF_State)                       :: importState, exportState
    type(ESMF_Clock)                       :: clock
    integer,intent(out)                    :: rc
!
!***  local variables
!
    integer       :: i

    rc     = ESMF_SUCCESS
!
    ! importable fields:
    do i = 1, size(importFieldsInfo)
      call NUOPC_Advertise(importState, &
                           StandardName=trim(importFieldsInfo(i)%name), &
                           SharePolicyField='share', rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    end do

    ! exportable fields:
    do i = 1, size(exportFieldsInfo)
      call NUOPC_Advertise(exportState, &
                           StandardName=trim(exportFieldsInfo(i)%name), &
                           SharePolicyField='share', rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    end do

!
!-----------------------------------------------------------------------
!
   end subroutine init_advertise
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
  subroutine init_realize(nest, importState, exportState, clock, rc)
!

    use module_cplscalars, only : flds_scalar_name, flds_scalar_num,          &
                                  flds_scalar_index_nx, flds_scalar_index_ny, &
                                  flds_scalar_index_ntile
    use module_cplscalars, only : State_SetScalar

    type(ESMF_GridComp)                    :: nest
    type(ESMF_State)                       :: importState, exportState
    type(ESMF_Clock)                       :: clock
    integer,intent(out)                    :: rc
!
!***  local variables
!
    real(ESMF_KIND_R8)  :: scalardim(3)
    type(ESMF_Grid)     :: grid

    scalardim = 0.0
    ! cpl_scalars for export state
    scalardim(1) = real(Atmos%mlon,8)
    scalardim(2) = real(Atmos%mlat,8)
    scalardim(3) = 1.0
    if (.not. Atmos%regional)scalardim(3) = 6.0

    rc     = ESMF_SUCCESS

    ! access this domain grid
    call ESMF_GridCompGet(nest, grid=grid, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,  file=__FILE__)) return

    ! -- realize connected fields in exportState
    call realizeConnectedCplFields(exportState, grid, &
                                   numLevels, numSoilLayers, numTracers, &
                                   exportFieldsInfo, 'FV3 Export', exportFields, 0.0_ESMF_KIND_R8, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,  file=__FILE__)) return

    if (flds_scalar_num > 0) then
      ! Set the scalar data into the exportstate
      call State_SetScalar(scalardim(1), flds_scalar_index_nx, exportState, flds_scalar_name, flds_scalar_num, rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,  file=__FILE__)) return
      call State_SetScalar(scalardim(2), flds_scalar_index_ny, exportState, flds_scalar_name, flds_scalar_num, rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,  file=__FILE__)) return
      call State_SetScalar(scalardim(3), flds_scalar_index_ntile, exportState, flds_scalar_name, flds_scalar_num, rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,  file=__FILE__)) return
    end if

    ! -- initialize export fields if applicable
    call setup_exportdata(rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,  file=__FILE__)) return

    ! -- realize connected fields in importState
    call realizeConnectedCplFields(importState, grid, &
                                   numLevels, numSoilLayers, numTracers, &
                                   importFieldsInfo, 'FV3 Import', importFields, 9.99e20_ESMF_KIND_R8, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,  file=__FILE__)) return
!
!-----------------------------------------------------------------------
!
   end subroutine init_realize
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
  subroutine fcst_initialize(fcst_comp, importState, exportState, clock, rc)
!
!-----------------------------------------------------------------------
!***  INITIALIZE THE FORECAST GRIDDED COMPONENT.
!-----------------------------------------------------------------------
!
    type(esmf_GridComp)                    :: fcst_comp
    type(ESMF_State)                       :: importState, exportState
    type(esmf_Clock)                       :: clock
    integer,intent(out)                    :: rc
!
!***  local variables
!
    integer                                :: i, j
!
    type(ESMF_VM)                          :: VM
    type(ESMF_Time)                        :: CurrTime, StartTime, StopTime
    type(ESMF_Config)                      :: cf

    integer,dimension(6)                   :: date, date_end
!
    integer :: initClock, unit, total_inttime
    integer :: stat
    character(4) dateSY
    character(2) dateSM,dateSD,dateSH,dateSN,dateSS
    character(len=esmf_maxstr) name_FB, name_FB1
    character(len=80) :: dateS

    character(256)                         :: gridfile

    character(8) :: bundle_grid

    real(kind=8) :: mpi_wtime, timeis

    type(ESMF_DELayout) :: delayout
    type(ESMF_DistGrid) :: distgrid
    integer :: jsc, jec, isc, iec, nlev
    type(domain2D)  :: domain
    integer :: n, fcstNpes, tmpvar, k
    logical :: freq_restart, fexist
    integer, allocatable, dimension(:) :: isl, iel, jsl, jel
    integer, allocatable, dimension(:,:,:) :: deBlockList
    integer, allocatable, dimension(:) :: petListNest

    integer               :: globalTileLayout(2)
    integer               :: nestRootPet, peListSize(1)
    integer, allocatable  :: petMap(:)
    integer               :: layout(2), nx, ny
    integer, pointer      :: pelist(:) => null()
    logical               :: top_parent_is_global
    logical               :: history_file_on_native_grid

    integer                       :: num_restart_interval, restart_starttime
    real,dimension(:),allocatable :: restart_interval

    integer           :: urc
    type(ESMF_State)  :: tempState
    type(ESMF_Info)   :: info

    type(time_type)               :: Time_init, Time, Time_step, Time_end, &
                                     Time_restart, Time_step_restart
    type(time_type)               :: iautime
    integer                       :: io_unit, calendar_type_res, date_res(6), date_init_res(6)

    integer,allocatable           :: grid_number_on_all_pets(:)
    logical,allocatable           :: is_moving_on_all_pets(:), is_moving(:)
    character(len=7)              :: nest_suffix

    type(FmsNetcdfFile_t)         :: fileobj
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
    timeis = mpi_wtime()
    rc     = ESMF_SUCCESS
!
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

    num_restart_interval = ESMF_ConfigGetLen(config=CF, label ='restart_interval:',rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    if (mype == 0) print *,'af ufs config,num_restart_interval=',num_restart_interval
    if (num_restart_interval<=0) num_restart_interval = 1
    allocate(restart_interval(num_restart_interval))
    restart_interval = 0
    call ESMF_ConfigGetAttribute(CF,valueList=restart_interval,label='restart_interval:', &
                                 count=num_restart_interval, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    if (mype == 0) print *,'af ufs config,restart_interval=',restart_interval
!
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
    if(iau_offset > 0 ) then
      iautime =  set_time(iau_offset * 3600, 0)
    endif
! if the second item is -1, the first number is frequency
    freq_restart = .false.
    if(num_restart_interval == 2) then
      if(restart_interval(2)== -1) freq_restart = .true.
    endif
    if(freq_restart) then
      if(restart_interval(1) >= 0) then
        tmpvar = restart_interval(1) * 3600
        Time_step_restart = set_time (tmpvar, 0)
        if(iau_offset > 0 ) then
          Time_restart = Time_init + iautime + Time_step_restart
          frestart(1) = tmpvar + iau_offset *3600
        else
          Time_restart = Time_init + Time_step_restart
          frestart(1) = tmpvar
        endif
        if(restart_interval(1) > 0) then
          i = 2
          do while ( Time_restart < Time_end )
            frestart(i) = frestart(i-1) + tmpvar
            Time_restart = Time_restart + Time_step_restart
             i = i + 1
          enddo
        endif
      endif
! otherwise it is an array with forecast time at which the restart files will be written out
    else if(num_restart_interval >= 1) then
      if(num_restart_interval == 1 .and. restart_interval(1) == 0 ) then
        frestart(1) = total_inttime
      else
        if(iau_offset > 0 ) then
          restart_starttime = iau_offset *3600
        else
          restart_starttime = 0
        endif
        do i=1,num_restart_interval
          frestart(i) = restart_interval(i) * 3600. + restart_starttime
        enddo
      endif
    endif
! if to write out restart at the end of forecast
    if (mype == 0) print *,'frestart=',frestart(1:10)/3600, 'total_inttime=',total_inttime

!------ initialize component models ------

     call  atmos_model_init (Atmos, Time_init, Time, Time_step)
!
     inquire(FILE='data_table', EXIST=fexist)
     if (fexist) then
       call data_override_init()
     endif
!-----------------------------------------------------------------------
!---- open and close dummy file in restart dir to check if dir exists --

      if (mpp_pe() == 0 ) then
         open( newunit=unit, file='RESTART/file', iostat=stat )
         if (stat == 0) then
            close(unit, status='delete')
         else
            call mpp_error ( FATAL, 'fcst_initialize: RESTART subdirectory does not exist in the run directory' )
         endif
      endif
!
!-----------------------------------------------------------------------
!*** create grid for output fields, using FV3 parameters
!-----------------------------------------------------------------------
!
      call mpp_error(NOTE, 'before create fcst grid')

      gridfile = "grid_spec.nc" ! default

      if (open_file(fileobj, "INPUT/grid_spec.nc", "read")) then
        if (variable_exists(fileobj, "atm_mosaic_file")) then
          call read_data(fileobj, "atm_mosaic_file", gridfile)
        endif
        call close_file(fileobj)
      endif

      ngrids = Atmos%ngrids
      mygrid = Atmos%mygrid
      allocate(grid_number_on_all_pets(fcst_ntasks), is_moving_on_all_pets(fcst_ntasks))
      call mpi_allgather(mygrid, 1, MPI_INTEGER, &
                         grid_number_on_all_pets, 1, MPI_INTEGER, &
                         fcst_mpi_comm, rc)
      call mpi_allgather(Atmos%is_moving_nest, 1, MPI_LOGICAL, &
                         is_moving_on_all_pets, 1, MPI_LOGICAL, &
                         fcst_mpi_comm, rc)
      allocate(is_moving(ngrids))
      do n=1, fcst_ntasks
        is_moving(grid_number_on_all_pets(n)) = is_moving_on_all_pets(n)
      enddo
      deallocate(grid_number_on_all_pets, is_moving_on_all_pets)

      call ESMF_InfoGetFromHost(exportState, info=info, rc=rc); ESMF_ERR_ABORT(rc)
      call ESMF_InfoSet(info, key="is_moving", values=is_moving, rc=rc); ESMF_ERR_ABORT(rc)
      deallocate(is_moving)

      allocate (fcstGridComp(ngrids))
      do n=1,ngrids

        pelist => null()
        call atmos_model_get_nth_domain_info(n, layout, nx, ny, pelist)
        call ESMF_VMBroadcast(vm, bcstData=layout, count=2, rootPet=pelist(1), rc=rc); ESMF_ERR_ABORT(rc)

        if (n==1) then
           ! on grid==1 (top level parent) determine if the domain is global or regional
           top_parent_is_global = .true.
           if(mygrid==1) then
              if (Atmos%regional) top_parent_is_global = .false.
           endif
           call mpi_bcast(top_parent_is_global, 1, MPI_LOGICAL, 0, fcst_mpi_comm, rc)
        endif

        if (n==1 .and. top_parent_is_global) then

          fcstGridComp(n) = ESMF_GridCompCreate(name="global", petList=pelist, rc=rc); ESMF_ERR_ABORT(rc)

          call ESMF_InfoGetFromHost(fcstGridComp(n), info=info, rc=rc); ESMF_ERR_ABORT(rc)
          call ESMF_InfoSet(info, key="layout", values=layout, rc=rc); ESMF_ERR_ABORT(rc)
          call ESMF_InfoSet(info, key="tilesize", value=Atmos%mlon, rc=rc); ESMF_ERR_ABORT(rc)

          call ESMF_GridCompSetServices(fcstGridComp(n), SetServicesNest, userrc=urc, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
          if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

        else

          allocate(petListNest(layout(1)*layout(2)))
          k=pelist(1)
          do j=1,layout(2)
          do i=1,layout(1)
             petListNest(k-pelist(1)+1) = k
             k = k + 1
          end do
          end do

          fcstGridComp(n) = ESMF_GridCompCreate(name="nest", petList=petListNest, rc=rc); ESMF_ERR_ABORT(rc)

          call ESMF_InfoGetFromHost(fcstGridComp(n), info=info, rc=rc); ESMF_ERR_ABORT(rc)
          call ESMF_InfoSet(info, key="layout", values=layout, rc=rc); ESMF_ERR_ABORT(rc)
          call ESMF_InfoSet(info, key="nx", value=nx, rc=rc); ESMF_ERR_ABORT(rc)
          call ESMF_InfoSet(info, key="ny", value=ny, rc=rc); ESMF_ERR_ABORT(rc)

          call ESMF_GridCompSetServices(fcstGridComp(n), SetServicesNest, userrc=urc, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
          if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

          deallocate(petListNest)

        end if
      end do

! Add gridfile Attribute to the exportState
      call ESMF_AttributeAdd(exportState, convention="NetCDF", purpose="FV3", &
                               attrList=(/"gridfile"/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call ESMF_AttributeSet(exportState, convention="NetCDF", purpose="FV3", &
                               name="gridfile", value=trim(gridfile), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

! Add total number of domains(grids) Attribute to the exportState
      call ESMF_AttributeAdd(exportState, convention="NetCDF", purpose="FV3", &
                               attrList=(/"ngrids"/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call ESMF_AttributeSet(exportState, convention="NetCDF", purpose="FV3", &
                               name="ngrids", value=ngrids, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

! Add top_parent_is_global Attribute to the exportState
      call ESMF_AttributeAdd(exportState, convention="NetCDF", purpose="FV3", &
                               attrList=(/"top_parent_is_global"/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call ESMF_AttributeSet(exportState, convention="NetCDF", purpose="FV3", &
                               name="top_parent_is_global", value=top_parent_is_global, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

! Add time Attribute to the exportState
      call ESMF_AttributeAdd(exportState, convention="NetCDF", purpose="FV3", &
        attrList=(/ "time               ", &
                    "time:long_name     ", &
                    "time:units         ", &
                    "time:cartesian_axis", &
                    "time:calendar_type ", &
                    "time:calendar      " /), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call ESMF_AttributeSet(exportState, convention="NetCDF", purpose="FV3", &
                               name="time", value=real(0,ESMF_KIND_R8), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      write(dateSY,'(I4.4)')date_init(1)
      write(dateSM,'(I2.2)')date_init(2)
      write(dateSD,'(I2.2)')date_init(3)
      write(dateSH,'(I2.2)')date_init(4)
      write(dateSN,'(I2.2)')date_init(5)
      write(dateSS,'(I2.2)')date_init(6)

      dateS="hours since "//dateSY//'-'//dateSM//'-'//dateSD//' '//dateSH//':'//    &
            dateSN//":"//dateSS
      if (mype == 0) write(*,*)'dateS=',trim(dateS),'date_init=',date_init

      call ESMF_AttributeSet(exportState, convention="NetCDF", purpose="FV3", &
                               name="time:units", value=trim(dateS), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call ESMF_AttributeSet(exportState, convention="NetCDF", purpose="FV3", &
                               name="time:long_name", value="time", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call ESMF_AttributeSet(exportState, convention="NetCDF", purpose="FV3", &
                               name="time:cartesian_axis", value="T", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call ESMF_AttributeSet(exportState, convention="NetCDF", purpose="FV3", &
                               name="time:calendar_type", value=uppercase(trim(calendar)), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call ESMF_AttributeSet(exportState, convention="NetCDF", purpose="FV3", &
                               name="time:calendar", value=uppercase(trim(calendar)), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

! Add time_iso Attribute to the exportState
      call ESMF_AttributeAdd(exportState, convention="NetCDF", purpose="FV3", &
        attrList=(/ "time_iso               ", &
                    "time_iso:long_name     ", &
                    "time_iso:description   " /), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call ESMF_AttributeSet(exportState, convention="NetCDF", purpose="FV3", &
                               name="time_iso", value="yyyy-mm-ddThh:mm:ssZ", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call ESMF_AttributeSet(exportState, convention="NetCDF", purpose="FV3", &
                               name="time_iso:description", value="ISO 8601 Date String", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call ESMF_AttributeSet(exportState, convention="NetCDF", purpose="FV3", &
                               name="time_iso:long_name", value="valid time", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

! Create FieldBundle for Fields that need to be regridded bilinear
      if( quilting ) then

        call ESMF_ConfigGetAttribute(config=CF, value=history_file_on_native_grid, default=.false., label='history_file_on_native_grid:', rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        nbdlphys = 2

        do n=1,ngrids
        bundle_grid=''
        if (ngrids > 1 .and. n >= 2) then
          write(bundle_grid,'(A5,I2.2,A1)') '.nest', n, '.'
        endif

        do i=1,num_files
!
         tempState = ESMF_StateCreate(rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

         name_FB = trim(filename_base(i)) // trim(bundle_grid)
!
         if (i == 1) then ! for dyn
           name_FB1 = trim(name_FB)//'_bilinear'
           call create_bundle_and_add_it_to_state(trim(name_FB1), tempState, rc=rc)
           if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

           if (n == 1 .AND. top_parent_is_global .AND. history_file_on_native_grid) then
             call create_bundle_and_add_it_to_state('cubed_sphere_grid_'//trim(name_FB1), tempState, rc=rc)
             if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
           end if

           call ESMF_GridCompInitialize(fcstGridComp(n), importState=tempState,&
                                        exportState=exportState, phase=1, userrc=urc, rc=rc)
           if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
           if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

         else if (i == 2) then ! for phys

           do j=1, nbdlphys
             if (j == 1) then
               name_FB1 = trim(name_FB)//'_nearest_stod'
             else
               name_FB1 = trim(name_FB)//'_bilinear'
             endif
             call create_bundle_and_add_it_to_state(trim(name_FB1), tempState, rc=rc)
             if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

             if (n == 1 .AND. top_parent_is_global .AND. history_file_on_native_grid) then
               call create_bundle_and_add_it_to_state('cubed_sphere_grid_'//trim(name_FB1), tempState, rc=rc)
               if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
             endif

           enddo

           call ESMF_GridCompInitialize(fcstGridComp(n), importState=tempState,&
                                        exportState=exportState, phase=2, userrc=urc, rc=rc)
           if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
           if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

         else

           write(0,*)' unknown name_FB ', trim(name_FB)
           ESMF_ERR_ABORT(101)

         endif
!
         call ESMF_StateDestroy(tempState, noGarbage=.true., rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        enddo ! num_files history

        if ( quilting_restart ) then

          do i=1,3 ! 3 dynamics restart bundles

            tempState = ESMF_StateCreate(rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

            if (i == 1) then
              name_FB = 'restart_fv_core.res'
            elseif (i == 2) then
              name_FB = 'restart_fv_srf_wnd.res'
            elseif (i == 3) then
              name_FB = 'restart_fv_tracer.res'
            else
              write(0,*)' unknown name_dynamics restart bundle ', i
              ESMF_ERR_ABORT(101)
            endif

            if (n > 1) then
              write(nest_suffix,'(A5,I2.2)') '.nest', n
              name_FB = trim(name_FB)//nest_suffix
            endif

            call create_bundle_and_add_it_to_state(trim(name_FB), tempState, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

            call ESMF_GridCompInitialize(fcstGridComp(n), importState=tempState, &
                                         exportState=exportState, phase=1, userrc=urc, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
            if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

            call ESMF_StateDestroy(tempState, noGarbage=.true., rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          enddo ! 3 dynamics restart bundles

          do i=1,2 ! 2 physics restart bundles

            tempState = ESMF_StateCreate(rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

            if (i == 1) then
              name_FB = 'restart_phy_data'
            elseif (i == 2) then
              name_FB = 'restart_sfc_data'
            else
              write(0,*)' unknown name_physics restart bundle ', i
              ESMF_ERR_ABORT(101)
            endif

            if (n > 1) then
              write(nest_suffix,'(A5,I2.2)') '.nest', n
              name_FB = trim(name_FB)//nest_suffix
            endif

            call create_bundle_and_add_it_to_state(trim(name_FB), tempState, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

            call ESMF_GridCompInitialize(fcstGridComp(n), importState=tempState, &
                                         exportState=exportState, phase=2, userrc=urc, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
            if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

            call ESMF_StateDestroy(tempState, noGarbage=.true., rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          enddo ! 2 physics restart bundles
        endif ! quilting_restart

        enddo ! ngrids

      endif ! quilting

      call get_atmos_model_ungridded_dim(nlev=numLevels,         &
                                         nsoillev=numSoilLayers, &
                                         ntracers=numTracers)

      if (mype == 0) write(*,*)'fcst_initialize total time: ', mpi_wtime() - timeis
!
!-----------------------------------------------------------------------
!
   contains

     subroutine create_bundle_and_add_it_to_state(name_fb, state, rc)

       character(len=*), intent(in)    :: name_fb
       type(ESMF_State), intent(inout) :: state
       integer, intent(out)            :: rc

       type(ESMF_FieldBundle) :: fieldbundle

       fieldbundle = ESMF_FieldBundleCreate(name=trim(name_fb), rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

       call ESMF_AttributeAdd(fieldbundle, convention="NetCDF", purpose="FV3", attrList=(/"grid_id"/), rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

       call ESMF_AttributeSet(fieldbundle, convention="NetCDF", purpose="FV3", name="grid_id", value=n, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

       call ESMF_AttributeAdd(fieldbundle, convention="NetCDF", purpose="FV3-nooutput", attrList=(/"frestart"/), rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

       call ESMF_AttributeSet(fieldbundle, convention="NetCDF", purpose="FV3-nooutput", name="frestart", valueList=frestart, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

       call ESMF_StateAdd(state, (/fieldbundle/), rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

     end subroutine create_bundle_and_add_it_to_state

   end subroutine fcst_initialize
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
  subroutine fcst_advertise(fcst_comp, importState, exportState, clock, rc)
!
!-----------------------------------------------------------------------
!***  INITIALIZE THE FORECAST GRIDDED COMPONENT.
!-----------------------------------------------------------------------
!
    type(esmf_GridComp)                    :: fcst_comp
    type(ESMF_State)                       :: importState, exportState
    type(esmf_Clock)                       :: clock
    integer,intent(out)                    :: rc
!
!***  local variables
    type(ESMF_VM)     :: vm
    integer           :: n
    integer           :: urc

!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
    call ESMF_VMGetCurrent(vm=vm,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_VMGet(vm=vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    if (mype == 0) write(*,*)'fcst_advertise, cpl_grid_id=',cpl_grid_id

    call ESMF_GridCompInitialize(fcstGridComp(cpl_grid_id), importState=importState, &
                                 exportState=exportState, phase=3, userrc=urc, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
!
!-----------------------------------------------------------------------
!
   end subroutine fcst_advertise
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
  subroutine fcst_realize(fcst_comp, importState, exportState, clock, rc)
!
!-----------------------------------------------------------------------
!***  INITIALIZE THE FORECAST GRIDDED COMPONENT.
!-----------------------------------------------------------------------
!
    type(esmf_GridComp)                    :: fcst_comp
    type(ESMF_State)                       :: importState, exportState
    type(esmf_Clock)                       :: clock
    integer,intent(out)                    :: rc
!
!***  local variables
    type(ESMF_VM)     :: vm
    integer           :: n
    integer           :: urc

!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
    call ESMF_VMGetCurrent(vm=vm,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_VMGet(vm=vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    if (mype == 0) write(*,*)'fcst_realize, cpl_grid_id=',cpl_grid_id

    call ESMF_GridCompInitialize(fcstGridComp(cpl_grid_id), importState=importState, &
                                 exportState=exportState, phase=4, userrc=urc, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
!
!
!-----------------------------------------------------------------------
!
   end subroutine fcst_realize
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
   subroutine fcst_run_phase_1(fcst_comp, importState, exportState,clock,rc)
!
!-----------------------------------------------------------------------
!***  the run step for the fcst gridded component.
!-----------------------------------------------------------------------
!
      type(ESMF_GridComp)        :: fcst_comp
      type(ESMF_State)           :: importState, exportState
      type(ESMF_Clock)           :: clock
      integer,intent(out)        :: rc
!
!***  local variables
!
      logical,save               :: first=.true.
      integer,save               :: dt_cap=0
      type(ESMF_Time)            :: currTime,stopTime
      integer                    :: seconds
      real(kind=8)               :: mpi_wtime, tbeg1
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      tbeg1 = mpi_wtime()
      rc    = ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!
      call get_time(Atmos%Time - Atmos%Time_init, seconds)
      n_atmsteps = seconds/dt_atmos

      if (first) then
        call ESMF_ClockGet(clock, currTime=currTime, stopTime=stopTime, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        call ESMF_TimeIntervalGet(stopTime-currTime, s=dt_cap, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        first=.false.
      endif

      if ( dt_cap > 0 .and. mod(seconds, dt_cap) == 0 ) then
        Atmos%isAtCapTime = .true.
      else
        Atmos%isAtCapTime = .false.
      endif
!
!-----------------------------------------------------------------------
! *** call fcst integration subroutines

      call update_atmos_model_dynamics (Atmos)

      call update_atmos_radiation_physics (Atmos)

      call atmos_model_exchange_phase_1 (Atmos, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      if (mype == 0) write(*,'(A,I16,A,F16.6)')'PASS: fcstRUN phase 1, n_atmsteps = ', &
                                               n_atmsteps,' time is ',mpi_wtime()-tbeg1
!
!-----------------------------------------------------------------------
!
   end subroutine fcst_run_phase_1
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
   subroutine fcst_run_phase_2(fcst_comp, importState, exportState,clock,rc)
!
!-----------------------------------------------------------------------
!***  the run step for the fcst gridded component.
!-----------------------------------------------------------------------
!
      type(ESMF_GridComp)        :: fcst_comp
      type(ESMF_State)           :: importState, exportState
      type(ESMF_Clock)           :: clock
      integer,intent(out)        :: rc
!
!***  local variables
!
      integer                    :: date(6), seconds
      character(len=64)          :: timestamp
      integer                    :: unit
      real(kind=8)               :: mpi_wtime, tbeg1
!
      integer                                :: FBCount, i
      logical                                :: isPresent
      character(len=esmf_maxstr),allocatable :: itemNameList(:)
      type(ESMF_StateItem_Flag), allocatable :: itemTypeList(:)
      type(ESMF_FieldBundle)                 :: fcstExportFB
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      tbeg1 = mpi_wtime()
      rc    = ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!
! *** call fcst integration subroutines

      call atmos_model_exchange_phase_2 (Atmos, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call update_atmos_model_state (Atmos, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      !--- intermediate restart
      call get_time(Atmos%Time - Atmos%Time_init, seconds)
      if (ANY(frestart(:) == seconds)) then
          if (mype == 0) write(*,*)'write out restart at n_atmsteps=',n_atmsteps,' seconds=',seconds,  &
                                   'integration length=',n_atmsteps*dt_atmos/3600.

          timestamp = date_to_string (Atmos%Time)
          call atmos_model_restart(Atmos, timestamp)
          call write_stoch_restart_atm('RESTART/'//trim(timestamp)//'.atm_stoch.res.nc')

          !----- write coupler.res file ------
          if (.not. quilting_restart .and. mpp_pe() == mpp_root_pe()) then
              call get_date (Atmos%Time, date(1), date(2), date(3), date(4), date(5), date(6))
              open( newunit=unit, file='RESTART/'//trim(timestamp)//'.coupler.res' )
              write( unit, '(i6,8x,a)' )calendar_type, &
                   '(Calendar: no_calendar=0, thirty_day_months=1, julian=2, gregorian=3, noleap=4)'

              write( unit, '(6i6,8x,a)' )date_init, &
                   'Model start time:   year, month, day, hour, minute, second'
              write( unit, '(6i6,8x,a)' )date, &
                   'Current model time: year, month, day, hour, minute, second'
              close( unit )
          endif
      endif

      ! update fhzero
      call ESMF_StateGet(exportState, itemCount=FBCount, rc=rc)

      allocate (itemNameList(FBCount))
      allocate (itemTypeList(FBCount))
      call ESMF_StateGet(exportState, &
                         itemNameList=itemNameList, &
                         itemTypeList=itemTypeList, &
                         rc=rc)
      do i=1, FBcount
        if (itemTypeList(i) == ESMF_STATEITEM_FIELDBUNDLE) then
          call ESMF_StateGet(exportState, itemName=itemNameList(i), &
                             fieldbundle=fcstExportFB, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          call ESMF_AttributeGet(fcstExportFB, convention="NetCDF", purpose="FV3", &
                                 name="fhzero", isPresent=isPresent, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          if (isPresent) then
            call ESMF_AttributeSet(fcstExportFB, convention="NetCDF", purpose="FV3", name="fhzero", value=GFS_control%fhzero, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
          endif
        else
          !***### anything but a FieldBundle in the state is unexpected here
          call ESMF_LogSetError(ESMF_RC_ARG_BAD,                                 &
                                msg="Only FieldBundles supported in fcstState.", &
                                line=__LINE__, file=__FILE__, rcToReturn=rc)
          return
        endif

      enddo

      if (mype == 0) write(*,'(A,I16,A,F16.6)')'PASS: fcstRUN phase 2, n_atmsteps = ', &
                                               n_atmsteps,' time is ',mpi_wtime()-tbeg1
!
!-----------------------------------------------------------------------
!
   end subroutine fcst_run_phase_2
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
   subroutine fcst_finalize(fcst_comp, importState, exportState,clock,rc)
!
!-----------------------------------------------------------------------
!***  finalize the forecast grid component.
!-----------------------------------------------------------------------
!
      type(ESMF_GridComp)        :: fcst_comp
      type(ESMF_State)           :: importState, exportState
      type(ESMF_Clock)           :: clock
      integer,intent(out)        :: rc
!
!***  local variables
!
      integer                    :: unit
      integer,dimension(6)       :: date
      real(kind=8)               :: mpi_wtime, tbeg1
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      tbeg1 = mpi_wtime()
      rc    = ESMF_SUCCESS

      call atmos_model_end (Atmos)

      call diag_manager_end (Atmos%Time)

      call fms_end

      if (mype == 0) write(*,*)'fcst_finalize total time: ', mpi_wtime() - tbeg1
!
!-----------------------------------------------------------------------
!
  end subroutine fcst_finalize
!
!#######################################################################
!-- write forecast grid to NetCDF file for diagnostics
!
  subroutine wrt_fcst_grid(grid, fileName, relaxedflag, regridArea, rc)
    type(ESMF_Grid), intent(in)                      :: grid
    character(len=*), intent(in), optional           :: fileName
    logical, intent(in), optional                    :: relaxedflag
    logical, intent(in), optional                    :: regridArea
    integer, intent(out)                             :: rc
!
!***  local variables
!
    logical                     :: ioCapable
    logical                     :: doItFlag
    character(len=64)           :: lfileName
    character(len=64)           :: gridName
    type(ESMF_Array)            :: array
    type(ESMF_ArrayBundle)      :: arraybundle
    logical                     :: isPresent
    logical                     :: hasCorners
    logical                     :: lRegridArea
    type(ESMF_Field)            :: areaField
    type(ESMF_FieldStatus_Flag) :: areaFieldStatus

    ioCapable = (ESMF_IO_PIO_PRESENT .and. &
      (ESMF_IO_NETCDF_PRESENT .or. ESMF_IO_PNETCDF_PRESENT))
    doItFlag = .true.
    if (present(relaxedFlag)) then
      doItFlag = .not.relaxedflag .or. (relaxedflag.and.ioCapable)
    endif

    if (doItFlag) then
      ! Process optional arguments
      if (present(fileName)) then
        lfileName = trim(fileName)
      else
        call ESMF_GridGet(grid, name=gridName, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        lfileName = trim(gridName)//".nc"
      endif
      if (present(regridArea)) then
        lRegridArea = regridArea
      else
        lRegridArea = .FALSE.
      endif

      ! Create bundle for storing output
      arraybundle = ESMF_ArrayBundleCreate(rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return

      ! -- Centers --
      call ESMF_GridGetCoord(grid, staggerLoc=ESMF_STAGGERLOC_CENTER, &
        isPresent=isPresent, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      if (isPresent) then
        call ESMF_GridGetCoord(grid, coordDim=1, &
          staggerLoc=ESMF_STAGGERLOC_CENTER, array=array, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        call ESMF_ArraySet(array, name="lon_center", rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        call ESMF_ArrayBundleAdd(arraybundle,(/array/), rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        call ESMF_GridGetCoord(grid, coordDim=2, &
          staggerLoc=ESMF_STAGGERLOC_CENTER, array=array, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        call ESMF_ArraySet(array, name="lat_center", rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        call ESMF_ArrayBundleAdd(arraybundle,(/array/), rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
      endif

      ! -- Corners --
      call ESMF_GridGetCoord(grid, staggerLoc=ESMF_STAGGERLOC_CORNER, &
        isPresent=hasCorners, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      if (hasCorners) then
        call ESMF_GridGetCoord(grid, coordDim=1, &
          staggerLoc=ESMF_STAGGERLOC_CORNER, array=array, rc=rc)
        if (.not. ESMF_LogFoundError(rc, line=__LINE__, file=__FILE__)) then
          call ESMF_ArraySet(array, name="lon_corner", rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return
          call ESMF_ArrayBundleAdd(arraybundle,(/array/), rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return
        endif
        call ESMF_GridGetCoord(grid, coordDim=2, &
          staggerLoc=ESMF_STAGGERLOC_CORNER, array=array, rc=rc)
        if (.not. ESMF_LogFoundError(rc, line=__LINE__, file=__FILE__)) then
          call ESMF_ArraySet(array, name="lat_corner", rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return
          call ESMF_ArrayBundleAdd(arraybundle,(/array/), rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return
        endif
        if (lRegridArea) then
          areaField = ESMF_FieldCreate(grid=grid, &
            typekind=ESMF_TYPEKIND_R8, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
          call ESMF_FieldRegridGetArea(areaField, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return
          call ESMF_FieldGet(areaField, array=array, rc=rc)
          if (.not. ESMF_LogFoundError(rc, line=__LINE__, file=__FILE__)) then
            call ESMF_ArraySet(array, name="regrid_area", rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, file=__FILE__)) return
            call ESMF_ArrayBundleAdd(arraybundle,(/array/), rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, file=__FILE__)) return
          endif
        endif
      endif

      ! -- Mask --
      call ESMF_GridGetItem(grid, itemflag=ESMF_GRIDITEM_MASK, &
        staggerLoc=ESMF_STAGGERLOC_CENTER, isPresent=isPresent, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      if (isPresent) then
        call ESMF_GridGetItem(grid, staggerLoc=ESMF_STAGGERLOC_CENTER, &
          itemflag=ESMF_GRIDITEM_MASK, array=array, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        call ESMF_ArraySet(array, name="mask", rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        call ESMF_ArrayBundleAdd(arraybundle,(/array/), rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
      endif

      ! -- Area --
      call ESMF_GridGetItem(grid, itemflag=ESMF_GRIDITEM_AREA, &
        staggerLoc=ESMF_STAGGERLOC_CENTER, isPresent=isPresent, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      if (isPresent) then
        call ESMF_GridGetItem(grid, staggerLoc=ESMF_STAGGERLOC_CENTER, &
          itemflag=ESMF_GRIDITEM_AREA, array=array, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        call ESMF_ArraySet(array, name="area", rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        call ESMF_ArrayBundleAdd(arraybundle,(/array/), rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
      endif

      ! Write array bundle to grid file
      call ESMF_ArrayBundleWrite(arraybundle, fileName=trim(lfileName), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return

      ! Clean-up
      if (lRegridArea) then
        call ESMF_FieldGet(areaField, status=areaFieldStatus, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        if (areaFieldStatus.eq.ESMF_FIELDSTATUS_COMPLETE) then
          call ESMF_FieldDestroy(areaField, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        endif
      endif
      call ESMF_ArrayBundleDestroy(arraybundle,rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    endif
  end subroutine wrt_fcst_grid
!
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!
end module  module_fcst_grid_comp
!
!----------------------------------------------------------------------------
