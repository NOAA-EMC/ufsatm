module module_cplfields

  !-----------------------------------------------------------------------------
  ! This module contains the fv3 Coupling Fields: export and import
  !
  !-----------------------------------------------------------------------------

  use ESMF
  use NUOPC

  implicit none

  private

  type, public :: FieldInfo
    character(len=41) :: name
    character(len=1) :: type
  end type

! Export Fields ----------------------------------------

  ! Please specify fields as: FieldInfo("standard_name", "type")
  ! Field types should be provided according to the table below:
  !  g : soil levels (3D)
  !  i : interface (3D)
  !  l : model levels (3D)
  !  s : surface (2D)
  !  t : tracers (4D)
  integer,          public, parameter :: NexportFields = 121
  type(ESMF_Field), target, public    :: exportFields(NexportFields)

  type(FieldInfo), dimension(NexportFields), public, parameter :: exportFieldsInfo = [ &
    FieldInfo("inst_pres_interface                      ", "i"), &
    FieldInfo("inst_pres_levels                         ", "l"), &
    FieldInfo("inst_geop_interface                      ", "i"), &
    FieldInfo("inst_geop_levels                         ", "l"), &
    FieldInfo("inst_temp_levels                         ", "l"), &
    FieldInfo("inst_zonal_wind_levels                   ", "l"), &
    FieldInfo("inst_merid_wind_levels                   ", "l"), &
    FieldInfo("inst_omega_levels                        ", "l"), &
    FieldInfo("inst_tracer_mass_frac                    ", "t"), &
    FieldInfo("soil_type                                ", "s"), &
    FieldInfo("inst_pbl_height                          ", "s"), &
    FieldInfo("surface_cell_area                        ", "s"), &
    FieldInfo("inst_convective_rainfall_amount          ", "s"), &
    FieldInfo("inst_exchange_coefficient_heat_levels    ", "l"), &
    FieldInfo("inst_spec_humid_conv_tendency_levels     ", "l"), &
    FieldInfo("inst_ice_nonconv_tendency_levels         ", "l"), &
    FieldInfo("inst_liq_nonconv_tendency_levels         ", "l"), &
    FieldInfo("inst_cloud_frac_levels                   ", "l"), &
    FieldInfo("inst_friction_velocity                   ", "s"), &
    FieldInfo("inst_rainfall_amount                     ", "s"), &
    FieldInfo("inst_soil_moisture_content               ", "g"), &
    FieldInfo("inst_surface_soil_wetness                ", "s"), &
    FieldInfo("inst_up_sensi_heat_flx                   ", "s"), &
    FieldInfo("inst_lwe_snow_thickness                  ", "s"), &
    FieldInfo("vegetation_type                          ", "s"), &
    FieldInfo("inst_vegetation_area_frac                ", "s"), &
    FieldInfo("inst_surface_roughness                   ", "s"), &
    FieldInfo("mean_zonal_moment_flx_atm                ", "s"), &
    FieldInfo("mean_merid_moment_flx_atm                ", "s"), &
    FieldInfo("mean_sensi_heat_flx                      ", "s"), &
    FieldInfo("mean_laten_heat_flx                      ", "s"), &
    FieldInfo("mean_evap_rate                           ", "s"), &
    FieldInfo("mean_down_lw_flx                         ", "s"), &
    FieldInfo("mean_down_sw_flx                         ", "s"), &
    FieldInfo("mean_prec_rate                           ", "s"), &
    FieldInfo("inst_prec_rate                           ", "s"), &
    FieldInfo("inst_zonal_moment_flx                    ", "s"), &
    FieldInfo("inst_merid_moment_flx                    ", "s"), &
    FieldInfo("inst_sensi_heat_flx                      ", "s"), &
    FieldInfo("inst_laten_heat_flx                      ", "s"), &
    FieldInfo("inst_evap_rate                           ", "s"), &
    FieldInfo("inst_down_lw_flx                         ", "s"), &
    FieldInfo("inst_down_sw_flx                         ", "s"), &
    FieldInfo("inst_temp_height2m                       ", "s"), &
    FieldInfo("inst_spec_humid_height2m                 ", "s"), &
    FieldInfo("inst_zonal_wind_height10m                ", "s"), &
    FieldInfo("inst_merid_wind_height10m                ", "s"), &
    FieldInfo("inst_temp_height_surface                 ", "s"), &
    FieldInfo("inst_pres_height_surface                 ", "s"), &
    FieldInfo("inst_surface_height                      ", "s"), &
    FieldInfo("mean_net_lw_flx                          ", "s"), &
    FieldInfo("mean_net_sw_flx                          ", "s"), &
    FieldInfo("inst_net_lw_flx                          ", "s"), &
    FieldInfo("inst_net_sw_flx                          ", "s"), &
    FieldInfo("mean_down_sw_ir_dir_flx                  ", "s"), &
    FieldInfo("mean_down_sw_ir_dif_flx                  ", "s"), &
    FieldInfo("mean_down_sw_vis_dir_flx                 ", "s"), &
    FieldInfo("mean_down_sw_vis_dif_flx                 ", "s"), &
    FieldInfo("inst_down_sw_ir_dir_flx                  ", "s"), &
    FieldInfo("inst_down_sw_ir_dif_flx                  ", "s"), &
    FieldInfo("inst_down_sw_vis_dir_flx                 ", "s"), &
    FieldInfo("inst_down_sw_vis_dif_flx                 ", "s"), &
    FieldInfo("mean_net_sw_ir_dir_flx                   ", "s"), &
    FieldInfo("mean_net_sw_ir_dif_flx                   ", "s"), &
    FieldInfo("mean_net_sw_vis_dir_flx                  ", "s"), &
    FieldInfo("mean_net_sw_vis_dif_flx                  ", "s"), &
    FieldInfo("inst_net_sw_ir_dir_flx                   ", "s"), &
    FieldInfo("inst_net_sw_ir_dif_flx                   ", "s"), &
    FieldInfo("inst_net_sw_vis_dir_flx                  ", "s"), &
    FieldInfo("inst_net_sw_vis_dif_flx                  ", "s"), &
    FieldInfo("inst_land_sea_mask                       ", "s"), &
    FieldInfo("inst_temp_height_lowest                  ", "s"), &
    FieldInfo("inst_spec_humid_height_lowest            ", "s"), &
    FieldInfo("inst_zonal_wind_height_lowest            ", "s"), &
    FieldInfo("inst_merid_wind_height_lowest            ", "s"), &
    FieldInfo("inst_pres_height_lowest                  ", "s"), &
    FieldInfo("inst_height_lowest                       ", "s"), &
    FieldInfo("inst_fprec_rate                          ", "s"), &
    FieldInfo("openwater_frac_in_atm                    ", "s"), &
    FieldInfo("ice_fraction_in_atm                      ", "s"), &
    FieldInfo("lake_fraction                            ", "s"), &
    FieldInfo("ocean_fraction                           ", "s"), &
    FieldInfo("surface_snow_area_fraction               ", "s"), &
    FieldInfo("canopy_moisture_storage                  ", "s"), &
    FieldInfo("inst_aerodynamic_conductance             ", "s"), &
    FieldInfo("inst_canopy_resistance                   ", "s"), &
    FieldInfo("leaf_area_index                          ", "s"), &
    FieldInfo("temperature_of_soil_layer                ", "g"), &
    FieldInfo("height                                   ", "s"), &
    FieldInfo("inst_pres_height_lowest_from_phys        ", "s"), &
    FieldInfo("inst_spec_humid_height_lowest_from_phys  ", "s"), &
    FieldInfo("inst_prec_rate_conv                      ", "s"), &
    FieldInfo("inst_temp_height_lowest_from_phys        ", "s"), &
    FieldInfo("inst_exner_function_height_lowest        ", "s"), &
    FieldInfo("surface_friction_velocity                ", "s"), &
    ! FieldInfo("fraction_of_vegetation_category          ", "s"), &
    ! FieldInfo("number_of_vegetation_categories          ", "s"), &


    !  For JEDI
    ! dynamics
    FieldInfo("u                                        ", "l"), &
    FieldInfo("v                                        ", "l"), &
    FieldInfo("ua                                       ", "l"), &
    FieldInfo("va                                       ", "l"), &
    FieldInfo("t                                        ", "l"), &
    FieldInfo("delp                                     ", "l"), &
    FieldInfo("sphum                                    ", "l"), &
    FieldInfo("ice_wat                                  ", "l"), &
    FieldInfo("liq_wat                                  ", "l"), &
    FieldInfo("o3mr                                     ", "l"), &
    FieldInfo("phis                                     ", "s"), &
    FieldInfo("u_srf                                    ", "s"), &
    FieldInfo("v_srf                                    ", "s"), &
    ! physics
    FieldInfo("slmsk                                    ", "s"), &
    FieldInfo("weasd                                    ", "s"), &
    FieldInfo("tsea                                     ", "s"), &
    FieldInfo("vtype                                    ", "s"), &
    FieldInfo("stype                                    ", "s"), &
    FieldInfo("vfrac                                    ", "s"), &
    FieldInfo("stc                                      ", "g"), &
    FieldInfo("smc                                      ", "g"), &
    FieldInfo("snwdph                                   ", "s"), &
    FieldInfo("f10m                                     ", "s"), &
    FieldInfo("zorl                                     ", "s"), &
    FieldInfo("t2m                                      ", "s"), &
    FieldInfo("cpl_scalars                              ", "s")]

! Import Fields ----------------------------------------
  integer,          public, parameter :: NimportFields = 67
  logical,          public            :: importFieldsValid(NimportFields)
  type(ESMF_Field), target, public    :: importFields(NimportFields)

  type(FieldInfo), dimension(NimportFields), public, parameter :: importFieldsInfo = [ &
    FieldInfo("inst_tracer_mass_frac                    ", "t"), &
    FieldInfo("land_mask                                ", "s"), &
    FieldInfo("sea_ice_surface_temperature              ", "s"), &
    FieldInfo("sea_surface_temperature                  ", "s"), &
    FieldInfo("ice_fraction                             ", "s"), &
    FieldInfo("lwup_flx_ice                             ", "s"), &
    FieldInfo("laten_heat_flx_atm_into_ice              ", "s"), &
    FieldInfo("sensi_heat_flx_atm_into_ice              ", "s"), &
    FieldInfo("stress_on_air_ice_zonal                  ", "s"), &
    FieldInfo("stress_on_air_ice_merid                  ", "s"), &
    FieldInfo("sea_ice_volume                           ", "s"), &
    FieldInfo("snow_volume_on_sea_ice                   ", "s"), &
    FieldInfo("inst_ice_ir_dif_albedo                   ", "s"), &
    FieldInfo("inst_ice_ir_dir_albedo                   ", "s"), &
    FieldInfo("inst_ice_vis_dif_albedo                  ", "s"), &
    FieldInfo("inst_ice_vis_dir_albedo                  ", "s"), &
    FieldInfo("wave_z0_roughness_length                 ", "s"), &
    FieldInfo("inst_tracer_diag_aod                     ", "s"), &
    FieldInfo("ocn_current_zonal                        ", "s"), &
    FieldInfo("ocn_current_merid                        ", "s"), &

    ! For receiving fluxes from mediator
    FieldInfo("stress_on_air_ocn_zonal                  ", "s"), &
    FieldInfo("stress_on_air_ocn_merid                  ", "s"), &
    FieldInfo("laten_heat_flx_atm_into_ocn              ", "s"), &
    FieldInfo("sensi_heat_flx_atm_into_ocn              ", "s"), &
    FieldInfo("lwup_flx_ocn                             ", "s"), &

    ! For receiving fluxes from external land component
    FieldInfo("land_fraction                            ", "s"), &
    FieldInfo("inst_snow_area_fraction_lnd              ", "s"), &
    FieldInfo("inst_spec_humid_lnd                      ", "s"), &
    FieldInfo("inst_laten_heat_flx_lnd                  ", "s"), &
    FieldInfo("inst_sensi_heat_flx_lnd                  ", "s"), &
    FieldInfo("inst_potential_laten_heat_flx_lnd        ", "s"), &
    FieldInfo("inst_temp_height2m_lnd                   ", "s"), &
    FieldInfo("inst_spec_humid_height2m_lnd             ", "s"), &
    FieldInfo("inst_upward_heat_flux_lnd                ", "s"), &
    FieldInfo("inst_runoff_rate_lnd                     ", "s"), &
    FieldInfo("inst_subsurface_runoff_rate_lnd          ", "s"), &
    FieldInfo("inst_drag_wind_speed_for_momentum        ", "s"), &
    FieldInfo("inst_drag_mass_flux_for_heat_and_moisture", "s"), &
    FieldInfo("inst_func_of_roughness_length_and_vfrac  ", "s"), &

    !  For JEDI
    ! dynamics
    FieldInfo("u                                        ", "l"), &
    FieldInfo("v                                        ", "l"), &
    FieldInfo("ua                                       ", "l"), &
    FieldInfo("va                                       ", "l"), &
    FieldInfo("t                                        ", "l"), &
    FieldInfo("delp                                     ", "l"), &
    FieldInfo("sphum                                    ", "l"), &
    FieldInfo("ice_wat                                  ", "l"), &
    FieldInfo("liq_wat                                  ", "l"), &
    FieldInfo("o3mr                                     ", "l"), &
    FieldInfo("phis                                     ", "s"), &
    FieldInfo("u_srf                                    ", "s"), &
    FieldInfo("v_srf                                    ", "s"), &
    ! physics
    FieldInfo("slmsk                                    ", "s"), &
    FieldInfo("weasd                                    ", "s"), &
    FieldInfo("tsea                                     ", "s"), &
    FieldInfo("vtype                                    ", "s"), &
    FieldInfo("stype                                    ", "s"), &
    FieldInfo("vfrac                                    ", "s"), &
    FieldInfo("stc                                      ", "g"), &
    FieldInfo("smc                                      ", "g"), &
    FieldInfo("snwdph                                   ", "s"), &
    FieldInfo("f10m                                     ", "s"), &
    FieldInfo("zorl                                     ", "s"), &
    FieldInfo("t2m                                      ", "s"), &

    ! For FIRE
    FieldInfo("hflx_fire                                ", "s"), &
    FieldInfo("evap_fire                                ", "s"), &
    FieldInfo("smoke_fire                               ", "s") ]

! Fields exported exclusively for coupling with chemistry
  character(*), public, parameter :: chemistryFieldNames(*) = [ &
    "inst_pres_interface             ", &
    "inst_pres_levels                ", &
    "inst_geop_interface             ", &
    "inst_geop_levels                ", &
    "inst_temp_levels                ", &
    "inst_zonal_wind_levels          ", &
    "inst_merid_wind_levels          ", &
    "inst_tracer_mass_frac           ", &
    "inst_pbl_height                 ", &
    "surface_cell_area               ", &
    "inst_convective_rainfall_amount ", &
    "inst_friction_velocity          ", &
    "inst_rainfall_amount            ", &
    "inst_up_sensi_heat_flx          ", &
    "inst_surface_roughness          ", &
    "inst_soil_moisture_content      ", &
    "inst_liq_nonconv_tendency_levels", &
    "inst_ice_nonconv_tendency_levels", &
    "inst_cloud_frac_levels          ", &
    "inst_surface_soil_wetness       ", &
    "ice_fraction_in_atm             ", &
    "lake_fraction                   ", &
    "ocean_fraction                  ", &
    "surface_snow_area_fraction      ", &
    "inst_vegetation_area_frac       ", &
    "canopy_moisture_storage         ", &
    "inst_aerodynamic_conductance    ", &
    "inst_canopy_resistance          ", &
    "leaf_area_index                 ", &
    "soil_type                       ", &
    "temperature_of_soil_layer       ", &
    "height                          ", &
    "vegetation_type                 " &
    ! "number_of_vegetation_categories ", &
    ! "fraction_of_vegetation_category "  &
    ]

  ! Methods
  public queryImportFields, queryExportFields
  public cplFieldGet
  public realizeConnectedCplFields

!-----------------------------------------------------------------------------
  contains
!-----------------------------------------------------------------------------
  integer function queryExportFields(fieldname, abortflag)

    character(len=*),intent(in) :: fieldname
    logical, optional           :: abortflag

    queryExportFields = queryFieldList(exportFieldsInfo, fieldname, abortflag)

  end function queryExportFields

!-----------------------------------------------------------------------------

  integer function queryImportFields(fieldname, abortflag)

    character(len=*),intent(in) :: fieldname
    logical, optional           :: abortflag

    queryImportFields = queryFieldList(importFieldsInfo, fieldname, abortflag)

  end function queryImportFields

!-----------------------------------------------------------------------------

  integer function queryFieldList(fieldsInfo, fieldname, abortflag)
    ! returns integer index of first found fieldname in fieldlist
    ! by default, will abort if field not found, set abortflag to false
    ! to turn off the abort.
    ! return value of < 1 means the field was not found

    type(FieldInfo) ,intent(in) :: fieldsInfo(:)
    character(len=*),intent(in) :: fieldname
    logical, optional           :: abortflag

    integer :: n
    logical :: labort
    integer :: rc

    labort = .true.
    if (present(abortflag)) then
      labort = abortflag
    endif

    queryFieldList = 0
    n = 1
    do while (queryFieldList < 1 .and. n <= size(fieldsInfo))
      if (trim(fieldsInfo(n)%name) == trim(fieldname)) then
        queryFieldList = n
      else
        n = n + 1
      endif
    enddo

    if (labort .and. queryFieldList < 1) then
      call ESMF_LogWrite('queryFieldList ABORT on fieldname '//trim(fieldname), &
                          ESMF_LOGMSG_INFO, line=__LINE__, file=__FILE__, rc=rc)
      CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    endif
  end function queryFieldList

!-----------------------------------------------------------------------------

  subroutine cplStateGet(state, fieldList, fieldCount, rc)

    character(len=*), intent(in)            :: state
    type(ESMF_Field), pointer,     optional :: fieldList(:)
    integer,          intent(out), optional :: fieldCount
    integer,          intent(out), optional :: rc

    !--- begin
    if (present(rc)) rc = ESMF_SUCCESS

    select case (trim(state))
      case ('import','i')
        if (present(fieldList )) fieldList  => importFields
        if (present(fieldCount)) fieldCount =  size(importFields)
      case ('export','o')
        if (present(fieldList )) fieldList  => exportFields
        if (present(fieldCount)) fieldCount =  size(exportFields)
      case default
        call ESMF_LogSetError(ESMF_RC_ARG_OUTOFRANGE, &
          msg="state argument can only be import(i)/export(o).", &
          line=__LINE__, file=__FILE__, rcToReturn=rc)
        return
    end select

  end subroutine cplStateGet

!-----------------------------------------------------------------------------

  subroutine cplFieldGet(state, name, localDe, &
                         farrayPtr2d, farrayPtr3d, farrayPtr4d, rc)

    character(len=*),   intent(in)            :: state
    character(len=*),   intent(in)            :: name
    integer,            intent(in),  optional :: localDe
    real(ESMF_KIND_R8), pointer,     optional :: farrayPtr2d(:,:)
    real(ESMF_KIND_R8), pointer,     optional :: farrayPtr3d(:,:,:)
    real(ESMF_KIND_R8), pointer,     optional :: farrayPtr4d(:,:,:,:)
    integer,            intent(out), optional :: rc

    !--- local variables
    integer                    :: localrc
    integer                    :: de, item, fieldCount, rank
    logical                    :: isCreated
    type(ESMF_Field), pointer  :: fieldList(:)
    character(len=ESMF_MAXSTR) :: fieldName

    !--- begin
    if (present(rc)) rc = ESMF_SUCCESS

    if (present(farrayPtr2d)) nullify(farrayPtr2d)
    if (present(farrayPtr3d)) nullify(farrayPtr3d)
    if (present(farrayPtr4d)) nullify(farrayPtr4d)

    de = 0
    if (present(localDe)) de = localDe

    call cplStateGet(state, fieldList=fieldList, fieldCount=fieldCount, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

    do item = 1, fieldCount
      isCreated = ESMF_FieldIsCreated(fieldList(item), rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
      if (isCreated) then
        call ESMF_FieldGet(fieldList(item), name=fieldName, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
        if (trim(fieldName) == trim(name)) then
          call ESMF_FieldGet(fieldList(item), rank=rank, rc=localrc)
          if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
          select case (rank)
            case (2)
              if (present(farrayPtr2d)) then
                call ESMF_FieldGet(fieldList(item), localDe=de, farrayPtr=farrayPtr2d, rc=localrc)
                if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
              end if
            case (3)
              if (present(farrayPtr3d)) then
                call ESMF_FieldGet(fieldList(item), localDe=de, farrayPtr=farrayPtr3d, rc=localrc)
                if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
              end if
            case (4)
              if (present(farrayPtr4d)) then
                call ESMF_FieldGet(fieldList(item), localDe=de, farrayPtr=farrayPtr4d, rc=localrc)
                if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
              end if
            case default
              call ESMF_LogSetError(ESMF_RC_NOT_IMPL, msg="field rank should be 2, 3, or 4.", &
                                    line=__LINE__, file=__FILE__, rcToReturn=rc)
              return
          end select
          exit
        end if
      end if
    end do

  end subroutine cplFieldGet


  subroutine realizeConnectedCplFields(state, grid, &
                                       numLevels, numSoilLayers, numTracers, &
                                       fields_info, state_tag, fieldList, fill_value, rc)

    use field_manager_mod,  only: MODEL_ATMOS
    use tracer_manager_mod, only: get_number_tracers, get_tracer_names
    use module_cplscalars,  only: flds_scalar_name, flds_scalar_num, SetScalarField

    type(ESMF_State),            intent(inout)  :: state
    type(ESMF_Grid),                intent(in)  :: grid
    integer,                        intent(in)  :: numLevels
    integer,                        intent(in)  :: numSoilLayers
    integer,                        intent(in)  :: numTracers
    type(FieldInfo), dimension(:),  intent(in)  :: fields_info
    character(len=*),               intent(in)  :: state_tag                              !< Import or export.
    type(ESMF_Field), dimension(:), intent(out) :: fieldList
    real(ESMF_KIND_R8), optional  , intent(in)  :: fill_value
    integer,                        intent(out) :: rc

    ! local variables

    integer          :: item, pos, tracerCount
    logical          :: isConnected
    type(ESMF_Field) :: field
    real(ESMF_KIND_R8) :: l_fill_value
    real(ESMF_KIND_R8), parameter :: d_fill_value = 0._ESMF_KIND_R8
    type(ESMF_StateIntent_Flag) :: stateintent
    character(len=32), allocatable, dimension(:) :: tracerNames, tracerUnits

    ! begin
    rc = ESMF_SUCCESS

    if (present(fill_value)) then
      l_fill_value = fill_value
    else
      l_fill_value = d_fill_value
    end if

    ! attach list of tracer names to exported tracer field as metadata
    call ESMF_StateGet(state, stateintent=stateintent, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    if (stateintent == ESMF_STATEINTENT_EXPORT) then
      call get_number_tracers(MODEL_ATMOS, num_tracers=tracerCount)
      allocate(tracerNames(tracerCount), tracerUnits(tracerCount))
      do item = 1, tracerCount
        call get_tracer_names(MODEL_ATMOS, item, tracerNames(item), units=tracerUnits(item))
      end do
    end if

    do item = 1, size(fields_info)
      isConnected = NUOPC_IsConnected(state, fieldName=trim(fields_info(item)%name), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      if (isConnected) then
        if (trim(fields_info(item)%name) == trim(flds_scalar_name)) then
          ! Create the scalar field
          call SetScalarField(field, flds_scalar_name, flds_scalar_num, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
        else
          call ESMF_StateGet(state, field=field, itemName=trim(fields_info(item)%name), rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
          call ESMF_FieldEmptySet(field, grid=grid, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
          select case (fields_info(item)%type)
          case ('l','layer')
            call ESMF_FieldEmptyComplete(field, typekind=ESMF_TYPEKIND_R8, &
                 ungriddedLBound=(/1/), ungriddedUBound=(/numLevels/), rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
          case ('i','interface')
            call ESMF_FieldEmptyComplete(field, typekind=ESMF_TYPEKIND_R8, &
                 ungriddedLBound=(/1/), ungriddedUBound=(/numLevels+1/), rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
          case ('t','tracer')
            call ESMF_FieldEmptyComplete(field, typekind=ESMF_TYPEKIND_R8, &
                 ungriddedLBound=(/1, 1/), ungriddedUBound=(/numLevels, numTracers/), rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
            if (allocated(tracerNames)) then
              call addFieldMetadata(field, 'tracerNames', tracerNames, rc=rc)
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
            end if
            if (allocated(tracerUnits)) then
              call addFieldMetadata(field, 'tracerUnits', tracerUnits, rc=rc)
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
            end if
          case ('s','surface')
            call ESMF_FieldEmptyComplete(field, typekind=ESMF_TYPEKIND_R8, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
          case ('g','soil')
            call ESMF_FieldEmptyComplete(field, typekind=ESMF_TYPEKIND_R8, &
                 ungriddedLBound=(/1/), ungriddedUBound=(/numSoilLayers/), rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
          case default
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
                 msg="exportFieldType = '"//trim(fields_info(item)%type)//"' not recognized", &
                 line=__LINE__, file=__FILE__, rcToReturn=rc)
            return
          end select
        end if
        call NUOPC_Realize(state, field=field, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        ! -- initialize field value
        call ESMF_FieldFill(field, dataFillScheme="const", const1=l_fill_value, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        ! -- save field
        fieldList(item) = field
        call ESMF_LogWrite('realizeConnectedCplFields '//trim(state_tag)//' Field '//trim(fields_info(item)%name)  &
             // ' is connected ', ESMF_LOGMSG_INFO, line=__LINE__, file=__FILE__, rc=rc)
      else
        ! remove a not connected Field from State
        call ESMF_StateRemove(state, (/trim(fields_info(item)%name)/), rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
        call ESMF_LogWrite('realizeConnectedCplFields '//trim(state_tag)//' Field '//trim(fields_info(item)%name)  &
             // ' is not connected ', ESMF_LOGMSG_INFO, line=__LINE__, file=__FILE__, rc=rc)
      end if
    end do

    if (allocated(tracerNames)) deallocate(tracerNames)
    if (allocated(tracerUnits)) deallocate(tracerUnits)

  end subroutine realizeConnectedCplFields

!-----------------------------------------------------------------------------

  subroutine addFieldMetadata(field, key, values, rc)

    ! This subroutine implements a preliminary method to provide metadata to
    ! a coupled model that is accessing the field via reference sharing
    ! (NUOPC SharedStatusField=.true.). The method sets a (key, values) pair
    ! in the field's array ESMF_Info object to retrieve an array of strings
    ! encoding metadata.
    !
    ! Such a capability should be implemented in the standard NUOPC connector
    ! for more general applications, possibly providing access to the field's
    ! ESMF_Info object.

    type(ESMF_Field)               :: field
    character(len=*),  intent(in)  :: key
    character(len=*),  intent(in)  :: values(:)
    integer, optional, intent(out) :: rc

    ! local variable
    integer          :: localrc
    type(ESMF_Array) :: array
    type(ESMF_Info)  :: info

    ! begin
    if (present(rc)) rc = ESMF_SUCCESS

    call ESMF_FieldGet(field, array=array, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
    call ESMF_InfoGetFromHost(array, info, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
    call ESMF_InfoSet(info, key, values, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

  end subroutine addFieldMetadata
!
!------------------------------------------------------------------------------
!
end module module_cplfields
