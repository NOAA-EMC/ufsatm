module mod_atmos_util
  implicit none
  public :: get_atmos_tracer_types
contains
  !#######################################################################
  ! <SUBROUTINE NAME="get_atmos_tracer_types">
  ! <DESCRIPTION>
  !  Identify and return usage and type id of atmospheric tracers.
  !  Ids are defined as:
  !    0 = generic tracer
  !    1 = chemistry - prognostic
  !    2 = chemistry - diagnostic
  !
  !  Tracers are identified via the additional 'tracer_usage' keyword and
  !  their optional 'type' qualifier. A tracer is assumed prognostic if
  !  'type' is not provided. See examples from the field_table file below:
  !
  !  Prognostic tracer:
  !  ------------------
  !  "TRACER", "atmos_mod",    "so2"
  !            "longname",     "so2 mixing ratio"
  !            "units",        "ppm"
  !            "tracer_usage", "chemistry"
  !            "profile_type", "fixed", "surface_value=5.e-6" /
  !
  !  Diagnostic tracer:
  !  ------------------
  !  "TRACER", "atmos_mod",    "pm25"
  !            "longname",     "PM2.5"
  !            "units",        "ug/m3"
  !            "tracer_usage", "chemistry", "type=diagnostic"
  !            "profile_type", "fixed", "surface_value=5.e-6" /
  !
  !  For atmospheric chemistry, the order of both prognostic and diagnostic
  !  tracers is validated against the model's internal assumptions.
  !
  ! </DESCRIPTION>
  subroutine get_atmos_tracer_types(tracer_types)

    use field_manager_mod,  only: parse
    use tracer_manager_mod, only: query_method
    use field_manager_mod,  only: MODEL_ATMOS
    use mpp_mod,            only: mpp_error, FATAL
    use tracer_manager_mod, only: get_number_tracers
    
    integer, intent(out) :: tracer_types(:)

    !--- local variables
    logical :: found
    integer :: n, num_tracers, num_types
    integer :: id_max, id_min, id_num, ip_max, ip_min, ip_num
    character(len=32)  :: tracer_usage
    character(len=128) :: control, tracer_type

    !--- begin

    !--- validate array size
    call get_number_tracers(MODEL_ATMOS, num_tracers=num_tracers)

    if (size(tracer_types) < num_tracers) &
         call mpp_error(FATAL, 'insufficient size of tracer type array')

    !--- initialize tracer indices
    id_min = num_tracers + 1
    id_max = -id_min
    ip_min = id_min
    ip_max = id_max
    id_num = 0
    ip_num = 0

    do n = 1, num_tracers
       tracer_types(n) = 0
       found = query_method('tracer_usage',MODEL_ATMOS,n,tracer_usage,control)
       if (found) then
          if (trim(tracer_usage) == 'chemistry') then
             !--- set default to prognostic
             tracer_type = 'prognostic'
             num_types = parse(control, 'type', tracer_type)
             select case (trim(tracer_type))
             case ('diagnostic')
                tracer_types(n) = 2
                id_num = id_num + 1
                id_max = n
                if (id_num == 1) id_min = n
             case ('prognostic')
                tracer_types(n) = 1
                ip_num = ip_num + 1
                ip_max = n
                if (ip_num == 1) ip_min = n
             end select
          end if
       end if
    end do

    if (ip_num > 0) then
       !--- check if prognostic tracers are contiguous
       if (ip_num > ip_max - ip_min + 1) &
            call mpp_error(FATAL, 'prognostic chemistry tracers must be contiguous')
    end if

    if (id_num > 0) then
       !--- check if diagnostic tracers are contiguous
       if (id_num > id_max - id_min + 1) &
            call mpp_error(FATAL, 'diagnostic chemistry tracers must be contiguous')
    end if

    !--- prognostic tracers must precede diagnostic ones
    if (ip_max > id_min) &
         call mpp_error(FATAL, 'diagnostic chemistry tracers must follow prognostic ones')

  end subroutine get_atmos_tracer_types
  ! </SUBROUTINE>

end module mod_atmos_util
