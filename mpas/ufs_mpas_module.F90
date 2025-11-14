!> ###########################################################################################
!> \file ufs_mpas_module.F90
!>
!> Routines from the subdrivers for MPAS-A and CAM-SIMA have been adopted/modified here for
!> use within the UFS Weather Model.
!> MPAS-A Subdriver:    MPAS-Model/src/driver/mpas_subdriver.F
!> CAM-SIMA (external): src/dynamics/mpas/driver/dyn_mpas_subdriver.F90
!>                      (https://github.com/ESCOMP/CAM-SIMA/blob/development/)
!>
!> ###########################################################################################
module ufs_mpas_module
  use mpas_derived_types, only : core_type, domain_type, mpas_Clock_type
  use mpas_derived_types, only : mpas_time_type
  use mpas_kind_types,    only : StrKIND
  implicit none

  public

  !
  type(core_type),       pointer :: corelist   => null()
  type(domain_type),     pointer :: domain_ptr => null()
  type(mpas_Clock_type), pointer :: clock      => null()

  !
  type (MPAS_Time_Type) :: LBC_intv_end

  !
  character(StrKIND), allocatable :: constituent_name(:)
  integer, allocatable :: index_constituent_to_mpas_scalar(:)
  integer, allocatable :: index_mpas_scalar_to_constituent(:)
  logical, allocatable :: is_water_species(:)
  
  !> #########################################################################################
  !>
  !> #########################################################################################
  type :: var_info_type
     character(64) :: name = ''
     character(10) :: type = ''
     integer :: rank = 0
  end type var_info_type

  !> #########################################################################################
  !> This list corresponds to the "lbc_in" stream in core_atmosphere/Registry.xml
  !> It consists of variables that are members of the "lbc" structure.
  !> #########################################################################################
  type(var_info_type), parameter :: lbc_in_var_info_list(*) = [ &
       var_info_type('lbc_u'                           , 'real'      , 2), &
       var_info_type('lbc_w'                           , 'real'      , 2), &
       var_info_type('lbc_rho'                         , 'real'      , 2), &
       var_info_type('lbc_theta'                       , 'real'      , 2), &
       var_info_type('lbc_scalars'                     , 'real'      , 3)  &
       ]

  !> #########################################################################################
  !> This list corresponds to the "invariant" stream in MPAS registry.
  !> It consists of variables that are members of the "mesh" structure.
  !> #########################################################################################
  type(var_info_type), parameter :: invariant_var_info_list(*) = [ &
       var_info_type('angleEdge'                       , 'real'      , 1), &
       var_info_type('areaCell'                        , 'real'      , 1), &
       var_info_type('areaTriangle'                    , 'real'      , 1), &
       var_info_type('bdyMaskCell'                     , 'integer'   , 1), &
       var_info_type('bdyMaskEdge'                     , 'integer'   , 1), &
       var_info_type('bdyMaskVertex'                   , 'integer'   , 1), &
       var_info_type('cellTangentPlane'                , 'real'      , 3), &
       var_info_type('cell_gradient_coef_x'            , 'real'      , 2), &
       var_info_type('cell_gradient_coef_y'            , 'real'      , 2), &
       var_info_type('cellsOnCell'                     , 'integer'   , 2), &
       var_info_type('cellsOnEdge'                     , 'integer'   , 2), &
       var_info_type('cellsOnVertex'                   , 'integer'   , 2), &
       var_info_type('cf1'                             , 'real'      , 0), &
       var_info_type('cf2'                             , 'real'      , 0), &
       var_info_type('cf3'                             , 'real'      , 0), &
       var_info_type('coeffs_reconstruct'              , 'real'      , 3), &
       var_info_type('dcEdge'                          , 'real'      , 1), &
       var_info_type('defc_a'                          , 'real'      , 2), &
       var_info_type('defc_b'                          , 'real'      , 2), &
       var_info_type('deriv_two'                       , 'real'      , 3), &
       var_info_type('dss'                             , 'real'      , 2), &
       var_info_type('dvEdge'                          , 'real'      , 1), &
       var_info_type('dzu'                             , 'real'      , 1), &
       var_info_type('edgeNormalVectors'               , 'real'      , 2), &
       var_info_type('edgesOnCell'                     , 'integer'   , 2), &
       var_info_type('edgesOnEdge'                     , 'integer'   , 2), &
       var_info_type('edgesOnVertex'                   , 'integer'   , 2), &
       var_info_type('fEdge'                           , 'real'      , 1), &
       var_info_type('fVertex'                         , 'real'      , 1), &
       var_info_type('fzm'                             , 'real'      , 1), &
       var_info_type('fzp'                             , 'real'      , 1), &
       var_info_type('indexToCellID'                   , 'integer'   , 1), &
       var_info_type('indexToEdgeID'                   , 'integer'   , 1), &
       var_info_type('indexToVertexID'                 , 'integer'   , 1), &
       var_info_type('kiteAreasOnVertex'               , 'real'      , 2), &
       var_info_type('latCell'                         , 'real'      , 1), &
       var_info_type('latEdge'                         , 'real'      , 1), &
       var_info_type('latVertex'                       , 'real'      , 1), &
       var_info_type('localVerticalUnitVectors'        , 'real'      , 2), &
       var_info_type('lonCell'                         , 'real'      , 1), &
       var_info_type('lonEdge'                         , 'real'      , 1), &
       var_info_type('lonVertex'                       , 'real'      , 1), &
       var_info_type('meshDensity'                     , 'real'      , 1), &
       var_info_type('nEdgesOnCell'                    , 'integer'   , 1), &
       var_info_type('nEdgesOnEdge'                    , 'integer'   , 1), &
       var_info_type('nominalMinDc'                    , 'real'      , 0), &
       var_info_type('qv_init'                         , 'real'      , 1), &
       var_info_type('rdzu'                            , 'real'      , 1), &
       var_info_type('rdzw'                            , 'real'      , 1), &
       var_info_type('t_init'                          , 'real'      , 2), &
       var_info_type('u_init'                          , 'real'      , 1), &
       var_info_type('v_init'                          , 'real'      , 1), &
       var_info_type('verticesOnCell'                  , 'integer'   , 2), &
       var_info_type('verticesOnEdge'                  , 'integer'   , 2), &
       var_info_type('weightsOnEdge'                   , 'real'      , 2), &
       var_info_type('xCell'                           , 'real'      , 1), &
       var_info_type('xEdge'                           , 'real'      , 1), &
       var_info_type('xVertex'                         , 'real'      , 1), &
       var_info_type('yCell'                           , 'real'      , 1), &
       var_info_type('yEdge'                           , 'real'      , 1), &
       var_info_type('yVertex'                         , 'real'      , 1), &
       var_info_type('zCell'                           , 'real'      , 1), &
       var_info_type('zEdge'                           , 'real'      , 1), &
       var_info_type('zVertex'                         , 'real'      , 1), &
       var_info_type('zb'                              , 'real'      , 3), &
       var_info_type('zb3'                             , 'real'      , 3), &
       var_info_type('zgrid'                           , 'real'      , 2), &
       var_info_type('zxu'                             , 'real'      , 2), &
       var_info_type('zz'                              , 'real'      , 2)  &
    ]

  ! Whether a variable should be in input or restart can be determined by looking at
  ! the `atm_init_coupled_diagnostics` subroutine in MPAS.
  ! If a variable first appears on the LHS of an equation, it should be in restart.
  ! If a variable first appears on the RHS of an equation, it should be in input.
  ! The remaining ones of interest should be in output.

  !> #########################################################################################
  !> This list corresponds to the "input" stream in MPAS registry.
  !> It consists of variables that are members of the "diag" and "state" structure.
  !> Only variables that are specific to the "input" stream are included.
  !> #########################################################################################
  type(var_info_type), parameter :: input_var_info_list(*) = [ &
       var_info_type('Time'                            , 'real'      , 0), &
       var_info_type('initial_time'                    , 'character' , 0), &
       var_info_type('rho'                             , 'real'      , 2), &
       var_info_type('rho_base'                        , 'real'      , 2), &
       var_info_type('scalars'                         , 'real'      , 3), &
       var_info_type('theta'                           , 'real'      , 2), &
       var_info_type('theta_base'                      , 'real'      , 2), &
       var_info_type('u'                               , 'real'      , 2), &
       var_info_type('w'                               , 'real'      , 2), &
       var_info_type('xtime'                           , 'character' , 0)  &
       ]

  !> #########################################################################################
  !> This list corresponds to the "restart" stream in MPAS registry.
  !> It consists of variables that are members of the "diag" and "state" structure.
  !> Only variables that are specific to the "restart" stream are included.
  !> #########################################################################################
  type(var_info_type), parameter :: restart_var_info_list(*) = [ &
       var_info_type('exner'                           , 'real'      , 2), &
       var_info_type('exner_base'                      , 'real'      , 2), &
       var_info_type('pressure_base'                   , 'real'      , 2), &
       var_info_type('pressure_p'                      , 'real'      , 2), &
       var_info_type('rho_p'                           , 'real'      , 2), &
       var_info_type('rho_zz'                          , 'real'      , 2), &
       var_info_type('rtheta_base'                     , 'real'      , 2), &
       var_info_type('rtheta_p'                        , 'real'      , 2), &
       var_info_type('ru'                              , 'real'      , 2), &
       var_info_type('ru_p'                            , 'real'      , 2), &
       var_info_type('rw'                              , 'real'      , 2), &
       var_info_type('rw_p'                            , 'real'      , 2), &
       var_info_type('theta_m'                         , 'real'      , 2)  &
    ]

  !> #########################################################################################
  !> This list corresponds to the "output" stream in MPAS registry.
  !> It consists of variables that are members of the "diag" structure.
  !> Only variables that are specific to the "output" stream are included.
  !> #########################################################################################
  type(var_info_type), parameter :: output_var_info_list(*) = [ &
       var_info_type('divergence'                      , 'real'      , 2), &
       var_info_type('pressure'                        , 'real'      , 2), &
       var_info_type('relhum'                          , 'real'      , 2), &
       var_info_type('surface_pressure'                , 'real'      , 1), &
       var_info_type('uReconstructMeridional'          , 'real'      , 2), &
       var_info_type('uReconstructZonal'               , 'real'      , 2), &
       var_info_type('vorticity'                       , 'real'      , 2)  &
    ]
  
contains
  !> #########################################################################################
  !> Convert one or more values of any intrinsic data types to a character string for pretty
  !> printing.
  !> If `value` contains more than one element, the elements will be stringified, delimited by `separator`, then concatenated.
  !> If `value` contains exactly one element, the element will be stringified without using `separator`.
  !> If `value` contains zero element or is of unsupported data types, an empty character string is produced.
  !> If `separator` is not supplied, it defaults to ", " (i.e., a comma and a space).
  !> (KCW, 2024-02-04)
  !>
  !> \update: Dustin Swales April 2025 - Modified for use in UWM
  !>
  !> #########################################################################################
  pure function stringify(value, separator)
    use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64

    class(*), intent(in) :: value(:)
    character(*), optional, intent(in) :: separator
    character(:), allocatable :: stringify

    integer, parameter :: sizelimit = 1024

    character(:), allocatable :: buffer, delimiter, format
    character(:), allocatable :: value_c(:)
    integer :: i, n, offset

    if (present(separator)) then
       delimiter = separator
    else
       delimiter = ', '
    end if

    n = min(size(value), sizelimit)

    if (n == 0) then
       stringify = ''

       return
    end if

    select type (value)
    type is (character(*))
       allocate(character(len(value) * n + len(delimiter) * (n - 1)) :: buffer)

       buffer(:) = ''
       offset = 0

       ! Workaround for a bug in GNU Fortran >= 12. This is perhaps the manifestation of GCC Bugzilla Bug 100819.
       ! When a character string array is passed as the actual argument to an unlimited polymorphic dummy argument,
       ! its array index and length parameter are mishandled.
       allocate(character(len(value)) :: value_c(size(value)))

       value_c(:) = value(:)

       do i = 1, n
          if (len(delimiter) > 0 .and. i > 1) then
             buffer(offset + 1:offset + len(delimiter)) = delimiter
             offset = offset + len(delimiter)
          end if

          if (len_trim(adjustl(value_c(i))) > 0) then
             buffer(offset + 1:offset + len_trim(adjustl(value_c(i)))) = trim(adjustl(value_c(i)))
             offset = offset + len_trim(adjustl(value_c(i)))
          end if
       end do

       deallocate(value_c)
    type is (integer(int32))
       allocate(character(11 * n + len(delimiter) * (n - 1)) :: buffer)
       allocate(character(17 + len(delimiter) + floor(log10(real(n))) + 1) :: format)

       write(format, '(a, i0, 3a)') '(ss, ', n, '(i0, :, "', delimiter, '"))'
       write(buffer, format) value
    type is (integer(int64))
       allocate(character(20 * n + len(delimiter) * (n - 1)) :: buffer)
       allocate(character(17 + len(delimiter) + floor(log10(real(n))) + 1) :: format)

       write(format, '(a, i0, 3a)') '(ss, ', n, '(i0, :, "', delimiter, '"))'
       write(buffer, format) value
    type is (logical)
       allocate(character(1 * n + len(delimiter) * (n - 1)) :: buffer)
       allocate(character(13 + len(delimiter) + floor(log10(real(n))) + 1) :: format)

       write(format, '(a, i0, 3a)') '(', n, '(l1, :, "', delimiter, '"))'
       write(buffer, format) value
    type is (real(real32))
       allocate(character(13 * n + len(delimiter) * (n - 1)) :: buffer)

       if (maxval(abs(value)) < 1.0e5_real32) then
          allocate(character(20 + len(delimiter) + floor(log10(real(n))) + 1) :: format)
          write(format, '(a, i0, 3a)') '(ss, ', n, '(f13.6, :, "', delimiter, '"))'
       else
          allocate(character(23 + len(delimiter) + floor(log10(real(n))) + 1) :: format)
          write(format, '(a, i0, 3a)') '(ss, ', n, '(es13.6e2, :, "', delimiter, '"))'
       end if

       write(buffer, format) value
    type is (real(real64))
       allocate(character(13 * n + len(delimiter) * (n - 1)) :: buffer)

       if (maxval(abs(value)) < 1.0e5_real64) then
          allocate(character(20 + len(delimiter) + floor(log10(real(n))) + 1) :: format)
          write(format, '(a, i0, 3a)') '(ss, ', n, '(f13.6, :, "', delimiter, '"))'
       else
          allocate(character(23 + len(delimiter) + floor(log10(real(n))) + 1) :: format)
          write(format, '(a, i0, 3a)') '(ss, ', n, '(es13.6e2, :, "', delimiter, '"))'
       end if

       write(buffer, format) value
    class default
       stringify = ''

       return
    end select

    stringify = trim(buffer)
  end function stringify

  !> #########################################################################################
  !>
  !>  routine ufs_mpas_atm_update_bdy_tend
  !>
  !> \brief   Reads new boundary data and updates the LBC tendencies
  !> \author  Michael Duda
  !> \date    27 September 2016
  !> \details 
  !>  This routine reads from the 'lbc_in' stream all variables in the 'lbc'
  !>  pool. When called with firstCall=.true., the latest time before the
  !>  present is read into time level 2 of the lbc pool; otherwise, the
  !>  contents of time level 2 are shifted to time level 1, the earliest
  !>  time strictly later than the present is read into time level 2, and
  !>  the tendencies for all fields in the lbc pool are computed and stored
  !>  in time level 1.
  !>
  !> \update: Dustin Swales September 2025 - Modified for use in UWM
  !>
  !> #########################################################################################
  subroutine ufs_mpas_atm_update_bdy_tend(clock, block, firstCall, ierr)
    use mpas_constants,      only : rvord
    use mpas_stream_manager, only : mpas_stream_mgr_read
    use mpas_log,            only : mpas_log_write
    use mpas_derived_types,  only : MPAS_STREAM_MGR_NOERR, MPAS_LOG_ERR
    use mpas_derived_types,  only : mpas_pool_type, mpas_Clock_type, block_type
    use mpas_derived_types,  only : mpas_Time_type, MPAS_TimeInterval_type
    use mpas_timekeeping,    only : mpas_set_time
    use mpas_kind_types,     only : StrKIND, RKIND
    use mpas_derived_types,  only : MPAS_STREAM_LATEST_BEFORE
    use mpas_derived_types,  only : MPAS_STREAM_EARLIEST_STRICTLY_AFTER
    use mpas_timekeeping,    only : mpas_get_timeInterval, mpas_get_time, operator(-)
    use mpas_timekeeping,    only : mpas_get_clock_time, MPAS_NOW
    use mpas_pool_routines,  only : mpas_pool_get_config, mpas_pool_get_subpool
    use mpas_pool_routines,  only : mpas_pool_shift_time_levels, mpas_pool_get_array
    use mpas_pool_routines,  only : mpas_pool_get_dimension
    use module_mpas_config,  only : lbc_filename, pioid_lbc, pio_subsystem_lbc

    implicit none

    type (mpas_clock_type), intent(in) :: clock
    type (block_type), intent(inout) :: block
    logical, intent(in) :: firstCall
    integer, intent(out) :: ierr

    character(len=StrKIND) :: lbc_intv_start_string
    character(len=StrKIND) :: lbc_intv_end_string

    type (mpas_pool_type), pointer :: mesh
    type (mpas_pool_type), pointer :: state
    type (mpas_pool_type), pointer :: lbc
    real (kind=RKIND) :: dt

    integer, pointer :: nCells_ptr
    integer, pointer :: nEdges_ptr
    integer, pointer :: nVertLevels_ptr
    integer, pointer :: index_qv_ptr
    integer, pointer :: nScalars_ptr
    integer :: nCells, nEdges, nVertLevels, index_qv, nScalars

    real (kind=RKIND), dimension(:,:), pointer :: u
    real (kind=RKIND), dimension(:,:), pointer :: ru
    real (kind=RKIND), dimension(:,:), pointer :: rho_edge
    real (kind=RKIND), dimension(:,:), pointer :: w
    real (kind=RKIND), dimension(:,:), pointer :: theta
    real (kind=RKIND), dimension(:,:), pointer :: rtheta_m
    real (kind=RKIND), dimension(:,:), pointer :: rho_zz
    real (kind=RKIND), dimension(:,:), pointer :: rho
    real (kind=RKIND), dimension(:,:,:), pointer :: scalars
    real (kind=RKIND), dimension(:,:), pointer :: lbc_tend_u
    real (kind=RKIND), dimension(:,:), pointer :: lbc_tend_ru
    real (kind=RKIND), dimension(:,:), pointer :: lbc_tend_rho_edge
    real (kind=RKIND), dimension(:,:), pointer :: lbc_tend_w
    real (kind=RKIND), dimension(:,:), pointer :: lbc_tend_theta
    real (kind=RKIND), dimension(:,:), pointer :: lbc_tend_rtheta_m
    real (kind=RKIND), dimension(:,:), pointer :: lbc_tend_rho_zz
    real (kind=RKIND), dimension(:,:), pointer :: lbc_tend_rho
    real (kind=RKIND), dimension(:,:,:), pointer :: lbc_tend_scalars

    integer, dimension(:,:), pointer :: cellsOnEdge
    real (kind=RKIND), dimension(:,:), pointer :: zz

    integer :: dd_intv, s_intv, sn_intv, sd_intv
    type (MPAS_Time_Type) :: currTime
    type (MPAS_TimeInterval_Type) :: lbc_interval
    character(len=StrKIND) :: read_time
    integer :: iEdge, iCell, k, j
    integer :: cell1, cell2


    ierr = 0

    call mpas_pool_get_subpool(block % structs, 'mesh', mesh)
    call mpas_pool_get_subpool(block % structs, 'state', state)
    call mpas_pool_get_subpool(block % structs, 'lbc', lbc)

    if (firstCall) then
       call dyn_mpas_read_write_stream('r', 'lbc_in', pio_file_desc=pioid_lbc, ierr=ierr, timeLevel=2)
       if (ierr /= MPAS_STREAM_MGR_NOERR) then
          call mpas_log_write('Could not read from ''lbc_in'' stream on or before the current date '// &
                              'to update lateral boundary tendencies', messageType=MPAS_LOG_ERR)
          ierr = 1
       end if
    else
       call mpas_pool_shift_time_levels(lbc)
       call dyn_mpas_read_write_stream('r', 'lbc_in', pio_file_desc=pioid_lbc, ierr=ierr, timeLevel=2)
       if (ierr /= MPAS_STREAM_MGR_NOERR) then
          call mpas_log_write('Could not read from ''lbc_in'' stream after the current date '// &
                              'to update lateral boundary tendencies', messageType=MPAS_LOG_ERR)
          ierr = 1
       end if
    end if
    if (ierr /= 0) then
       return
    end if

    !read_time = '2023-03-10_00:00:00'
    !call mpas_set_time(currTime, dateTimeString=trim(read_time))
    currTime = mpas_get_clock_time(clock, MPAS_NOW, ierr)
    call mpas_get_time(currTime, dateTimeString=read_time, ierr=ierr)
    call mpas_set_time(currTime,dateTimeString=trim(read_time))
    !
    ! Compute any derived fields from those that were read from the lbc_in stream
    !
    call mpas_pool_get_array(lbc, 'lbc_u', u, 2)
    call mpas_pool_get_array(lbc, 'lbc_ru', ru, 2)
    call mpas_pool_get_array(lbc, 'lbc_rho_edge', rho_edge, 2)
    call mpas_pool_get_array(lbc, 'lbc_w', w, 2)
    call mpas_pool_get_array(lbc, 'lbc_theta', theta, 2)
    call mpas_pool_get_array(lbc, 'lbc_rtheta_m', rtheta_m, 2)
    call mpas_pool_get_array(lbc, 'lbc_rho_zz', rho_zz, 2)
    call mpas_pool_get_array(lbc, 'lbc_rho', rho, 2)
    call mpas_pool_get_array(lbc, 'lbc_scalars', scalars, 2)

    call mpas_pool_get_array(mesh, 'cellsOnEdge', cellsOnEdge)
    call mpas_pool_get_dimension(mesh, 'nCells', nCells_ptr)
    call mpas_pool_get_dimension(mesh, 'nEdges', nEdges_ptr)
    call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels_ptr)
    call mpas_pool_get_dimension(state, 'num_scalars', nScalars_ptr)
    call mpas_pool_get_dimension(lbc, 'index_qv', index_qv_ptr)
    call mpas_pool_get_array(mesh, 'zz', zz)

    if (.not. firstCall) then
       call mpas_pool_get_array(lbc, 'lbc_u', lbc_tend_u, 1)
       call mpas_pool_get_array(lbc, 'lbc_ru', lbc_tend_ru, 1)
       call mpas_pool_get_array(lbc, 'lbc_rho_edge', lbc_tend_rho_edge, 1)
       call mpas_pool_get_array(lbc, 'lbc_w', lbc_tend_w, 1)
       call mpas_pool_get_array(lbc, 'lbc_theta', lbc_tend_theta, 1)
       call mpas_pool_get_array(lbc, 'lbc_rtheta_m', lbc_tend_rtheta_m, 1)
       call mpas_pool_get_array(lbc, 'lbc_rho_zz', lbc_tend_rho_zz, 1)
       call mpas_pool_get_array(lbc, 'lbc_rho', lbc_tend_rho, 1)
       call mpas_pool_get_array(lbc, 'lbc_scalars', lbc_tend_scalars, 1)
    endif
    ! Dereference the pointers to avoid non-array pointer for OpenACC
    nCells = nCells_ptr
    nEdges = nEdges_ptr
    nVertLevels = nVertLevels_ptr
    nScalars = nScalars_ptr
    index_qv = index_qv_ptr

    ! Compute lbc_rho_zz
    do k=1,nVertLevels
       zz(k,nCells+1) = 1.0_RKIND          ! Avoid potential division by zero in the following line
    end do

    do iCell=1,nCells+1
       do k=1,nVertLevels
          rho_zz(k,iCell) = rho(k,iCell) / zz(k,iCell)
       end do
    end do

    ! Average lbc_rho_zz to edges
    do iEdge=1,nEdges
       cell1 = cellsOnEdge(1,iEdge)
       cell2 = cellsOnEdge(2,iEdge)
       if (cell1 > 0 .and. cell2 > 0) then
          do k = 1, nVertLevels
             rho_edge(k,iEdge) = 0.5_RKIND * (rho_zz(k,cell1) + rho_zz(k,cell2))
          end do
       end if
    end do

    do iEdge=1,nEdges+1
       do k=1,nVertLevels
          ru(k,iEdge) = u(k,iEdge) * rho_edge(k,iEdge)
       end do
    end do
    
    do iCell=1,nCells+1
       do k=1,nVertLevels
          rtheta_m(k,iCell) = theta(k,iCell) * rho_zz(k,iCell) * (1.0_RKIND + rvord * scalars(index_qv,k,iCell))
       end do
    end do

    if (.not. firstCall) then
       lbc_interval = currTime - LBC_intv_end
       call mpas_get_timeInterval(interval=lbc_interval, DD=dd_intv, S=s_intv, S_n=sn_intv, S_d=sd_intv, ierr=ierr)
       dt = 86400.0_RKIND * real(dd_intv, kind=RKIND) + real(s_intv, kind=RKIND) &
            + (real(sn_intv, kind=RKIND) / real(sd_intv, kind=RKIND))


       dt = 1.0_RKIND / dt

       do iEdge=1,nEdges+1
          do k=1,nVertLevels
             lbc_tend_u(k,iEdge) = (u(k,iEdge) - lbc_tend_u(k,iEdge)) * dt
             lbc_tend_ru(k,iEdge) = (ru(k,iEdge) - lbc_tend_ru(k,iEdge)) * dt
             lbc_tend_rho_edge(k,iEdge) = (rho_edge(k,iEdge) - lbc_tend_rho_edge(k,iEdge)) * dt
          end do
       end do

       do iCell=1,nCells+1
          do k=1,nVertLevels+1
             lbc_tend_w(k,iCell) = (w(k,iCell) - lbc_tend_w(k,iCell)) * dt
          end do
       end do

       do iCell=1,nCells+1
          do k=1,nVertLevels
             lbc_tend_theta(k,iCell) = (theta(k,iCell) - lbc_tend_theta(k,iCell)) * dt
             lbc_tend_rtheta_m(k,iCell) = (rtheta_m(k,iCell) - lbc_tend_rtheta_m(k,iCell)) * dt
             lbc_tend_rho_zz(k,iCell) = (rho_zz(k,iCell) - lbc_tend_rho_zz(k,iCell)) * dt
             lbc_tend_rho(k,iCell) = (rho(k,iCell) - lbc_tend_rho(k,iCell)) * dt
          end do
       end do
 
       do iCell=1,nCells+1
          do k=1,nVertLevels
             do j = 1,nScalars
                lbc_tend_scalars(j,k,iCell) = (scalars(j,k,iCell) - lbc_tend_scalars(j,k,iCell)) * dt
             end do
          end do
       end do

       !
       ! Logging the lbc start and end times appears to be backwards, but
       ! until the end of this function, LBC_intv_end == the last interval
       ! time and currTime == the next interval time.
       !
       call mpas_get_time(LBC_intv_end, dateTimeString=lbc_intv_start_string)
       call mpas_get_time(currTime, dateTimeString=lbc_intv_end_string)
       call mpas_log_write('----------------------------------------------------------------------')
       call mpas_log_write('Updated lateral boundary conditions. LBCs are now valid')
       call mpas_log_write('from '//trim(lbc_intv_start_string)//' to '//trim(lbc_intv_end_string))
       call mpas_log_write('----------------------------------------------------------------------')

    end if

    LBC_intv_end = currTime
    
  end subroutine ufs_mpas_atm_update_bdy_tend

!> ########################################################################################
 !>
 !> \brief  Computes local unit north, east, and edge-normal vectors
 !> \author Michael Duda
 !> \date   15 January 2020
 !> \details
 !>  This routine computes the local unit north and east vectors at all cell
 !>  centers, storing the resulting fields in the mesh pool as 'north' and
 !>  'east'. It also computes the edge-normal unit vectors by calling
 !>  the mpas_initialize_vectors routine. Before this routine is called,
 !>  the mesh pool must contain 'latCell' and 'lonCell' fields that are valid
 !>  for all cells (not just solve cells), plus any fields that are required
 !>  by the mpas_initialize_vectors routine.
 !>
 !> \update: Dustin Swales April 2025 - Modified for use in UWM
 !>
 !> ########################################################################################
 subroutine ufs_mpas_compute_unit_vectors()
   use mpas_pool_routines,     only : mpas_pool_get_subpool, mpas_pool_get_dimension, mpas_pool_get_array
   use mpas_derived_types,     only : mpas_pool_type
   use mpas_kind_types,        only : RKIND
   use mpas_vector_operations, only : mpas_initialize_vectors
   use module_mpas_config, only : nCellsSolve, latCell, lonCell
   
   type (mpas_pool_type), pointer :: meshPool
   real(kind=RKIND), dimension(:,:), pointer :: east, north
   integer, pointer :: nCells
   integer :: iCell

   call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh', meshPool)
   call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
   call mpas_pool_get_array(meshPool, 'latCell', latCell)
   call mpas_pool_get_array(meshPool, 'lonCell', lonCell)
   call mpas_pool_get_array(meshPool, 'east', east)
   call mpas_pool_get_array(meshPool, 'north', north)

   do iCell = 1, nCells
      east(1,iCell) = -sin(lonCell(iCell))
      east(2,iCell) =  cos(lonCell(iCell))
      east(3,iCell) =  0.0_RKIND

      ! Normalize
      east(1:3,iCell) = east(1:3,iCell) / sqrt(sum(east(1:3,iCell) * east(1:3,iCell)))

      north(1,iCell) = -cos(lonCell(iCell))*sin(latCell(iCell))
      north(2,iCell) = -sin(lonCell(iCell))*sin(latCell(iCell))
      north(3,iCell) =  cos(latCell(iCell))

      ! Normalize
      north(1:3,iCell) = north(1:3,iCell) / sqrt(sum(north(1:3,iCell) * north(1:3,iCell)))

   end do

   call mpas_initialize_vectors(meshPool)

 end subroutine ufs_mpas_compute_unit_vectors

 !> ########################################################################################
 !>
 !> \brief  Define the names of constituents at run-time
 !> \author Michael Duda
 !> \date   21 May 2020
 !> \details
 !>  Given an array of constituent names, which must have size equal to the number
 !>  of scalars that were set in the call to ufs_mpas_init_phase1, and given
 !>  a function to identify which scalars are moisture species, this routine defines
 !>  scalar constituents for the MPAS-A dycore.
 !>  Because the MPAS-A dycore expects all moisture constituents to appear in
 !>  a contiguous range of constituent indices, this routine may in general need
 !>  to reorder the constituents; to allow for mapping of indices between UFS
 !>  physics and the MPAS-A dycore, this routine returns index mapping arrays
 !>  mpas_from_ufs_cnst and ufs_from_mpas_cnst.
 !>
 !> \update: Dustin Swales April 2025 - Modified for use in UWM  
 !>
 !> ########################################################################################
 subroutine ufs_mpas_define_scalars(mpas_from_ufs_cnst, ufs_from_mpas_cnst, ierr)
   use mpas_derived_types, only : mpas_pool_type, field3dReal
   use mpas_pool_routines, only : mpas_pool_get_subpool, mpas_pool_get_field, &
                                  mpas_pool_get_dimension, mpas_pool_add_dimension
   use mpas_attlist,       only : mpas_add_att
   use mpas_log,           only : mpas_log_write
   use mpas_derived_types, only : MPAS_LOG_ERR
   ! FMS
   use mpp_mod,              only : FATAL, mpp_error
   
   ! Arguments
   integer, dimension(:), pointer :: mpas_from_ufs_cnst, ufs_from_mpas_cnst
   integer, intent(out) :: ierr

   ! Local variables
   character(len=*), parameter :: subname = 'ufs_mpas_subdriver::ufs_mpas_define_scalars'
   integer :: i, j, timeLevs
   integer, pointer :: num_scalars
   integer :: num_moist
   integer :: idx_passive
   type (mpas_pool_type), pointer :: statePool
   type (mpas_pool_type), pointer :: tendPool
   type (field3dReal), pointer :: scalarsField
   character(len=128) :: tempstr
   character :: moisture_char

   ierr = 0

   !
   ! Define scalars
   !
   nullify(statePool)
   call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'state', statePool)

   if (.not. associated(statePool)) then
      call mpas_log_write(trim(subname)//': ERROR: The ''state'' pool was not found.', &
                          messageType=MPAS_LOG_ERR)
      ierr = 1
      return
   end if

   nullify(num_scalars)
   call mpas_pool_get_dimension(statePool, 'num_scalars', num_scalars)

   !
   ! The num_scalars dimension should have been defined by atm_core_interface::atm_allocate_scalars, and
   ! if this dimension does not exist, something has gone wrong
   !
   if (.not. associated(num_scalars)) then
      call mpas_log_write(trim(subname)//': ERROR: The ''num_scalars'' dimension does not exist in the ''state'' pool.', &
                          messageType=MPAS_LOG_ERR)
      ierr = 1
      return
   end if

   !
   ! If at runtime there are not num_scalars names in the array of constituent names provided by UFS,
   ! something has gone wrong
   !
   if (size(constituent_name) /= num_scalars) then
      call mpas_log_write(trim(subname)//': ERROR: The number of constituent names is not equal to the num_scalars dimension', &
                          messageType=MPAS_LOG_ERR)
      call mpas_log_write('size(constituent_name) = $i, num_scalars = $i', intArgs=[size(constituent_name), num_scalars], &
                          messageType=MPAS_LOG_ERR)
      ierr = 1
      return
   end if

   !
   ! In UFS, the first scalar (if there are any) is always sphum (specific humidity); if this is not
   ! the case, something has gone wrong
   !
   if (size(constituent_name) > 0) then
      if (trim(constituent_name(1)) /= 'sphum') then
         call mpas_log_write(trim(subname)//': ERROR: The first constituent is not sphum', messageType=MPAS_LOG_ERR)
         ierr = 1
         return
      end if
   end if

   !
   ! Determine which of the constituents are moisture species
   !
   allocate(mpas_from_ufs_cnst(num_scalars), stat=ierr)
   if( ierr /= 0 ) call mpp_error(FATAL,subname//':failed to allocate mpas_from_ufs_cnst array')
   mpas_from_ufs_cnst(:) = 0
   num_moist = 0
   do i = 1, size(constituent_name)
      if (is_water_species(i)) then
         num_moist = num_moist + 1
         mpas_from_ufs_cnst(num_moist) = i
      end if
   end do

   !
   ! If UFS has no scalars, let the only scalar in MPAS be 'qv' (a moisture species)
   !
   if (num_scalars == 1 .and. size(constituent_name) == 0) then
      num_moist = 1
   end if

   !
   ! Assign non-moisture constituents to mpas_from_ufs_cnst(num_moist+1:size(constituent_name))
   !
   idx_passive = num_moist + 1
   do i = 1, size(constituent_name)

      ! If UFS constituent i is not already mapped as a moist constituent
      if (.not. is_water_species(i)) then
         mpas_from_ufs_cnst(idx_passive) = i
         idx_passive = idx_passive + 1
      end if
   end do

   !
   ! Create inverse map, ufs_from_mpas_cnst
   !
   allocate(ufs_from_mpas_cnst(num_scalars), stat=ierr)
   if( ierr /= 0 ) call mpp_error(FATAL,subname//':failed to allocate ufs_from_mpas_cnst array')
   ufs_from_mpas_cnst(:) = 0

   do i = 1, size(constituent_name)
      ufs_from_mpas_cnst(mpas_from_ufs_cnst(i)) = i
   end do

   timeLevs = 2

   do i = 1, timeLevs
      nullify(scalarsField)
      call mpas_pool_get_field(statePool, 'scalars', scalarsField, timeLevel=i)

      if (.not. associated(scalarsField)) then
         call mpas_log_write(trim(subname)//': ERROR: The ''scalars'' field was not found in the ''state'' pool', &
                             messageType=MPAS_LOG_ERR)
         ierr = 1
         return
      end if

      if (i == 1) call mpas_pool_add_dimension(statePool, 'index_qv', 1)
      scalarsField % constituentNames(1) = 'qv'
      call mpas_add_att(scalarsField % attLists(1) % attList, 'units', 'kg kg^{-1}')
      call mpas_add_att(scalarsField % attLists(1) % attList, 'long_name', 'Water vapor mixing ratio')

      do j = 2, size(constituent_name)
         scalarsField % constituentNames(j) = trim(constituent_name(mpas_from_ufs_cnst(j)))
      end do

   end do

   call mpas_pool_add_dimension(statePool, 'moist_start', 1)
   call mpas_pool_add_dimension(statePool, 'moist_end', num_moist)

   !
   ! Print a tabular summary of the mapping between constituent indices
   !
   call mpas_log_write('')
   call mpas_log_write('  i MPAS constituent mpas_from_ufs_cnst(i)       i UFS constituent  ufs_from_mpas_cnst(i)')
   call mpas_log_write('------------------------------------------     ------------------------------------------')
   do i = 1, min(num_scalars, size(constituent_name))
      if (i <= num_moist) then
         moisture_char = '*'
      else
         moisture_char = ' '
      end if
      write(tempstr, '(i3,1x,a16,1x,i18,8x,i3,1x,a16,1x,i18)') i, trim(scalarsField % constituentNames(i))//moisture_char, &
                                                               mpas_from_ufs_cnst(i), &
                                                               i, trim(constituent_name(i)), &
                                                               ufs_from_mpas_cnst(i)
      call mpas_log_write(trim(tempstr))
   end do
   call mpas_log_write('------------------------------------------     ------------------------------------------')
   call mpas_log_write('* = constituent used as a moisture species in MPAS-A dycore')
   call mpas_log_write('')


   !
   ! Define scalars_tend
   !
   nullify(tendPool)
   call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'tend', tendPool)

   if (.not. associated(tendPool)) then
      call mpas_log_write(trim(subname)//': ERROR: The ''tend'' pool was not found.', &
                          messageType=MPAS_LOG_ERR)
      ierr = 1
      return
   end if

   timeLevs = 1

   do i = 1, timeLevs
      nullify(scalarsField)
      call mpas_pool_get_field(tendPool, 'scalars_tend', scalarsField, timeLevel=i)

      if (.not. associated(scalarsField)) then
         call mpas_log_write(trim(subname)//': ERROR: The ''scalars_tend'' field was not found in the ''tend'' pool', &
                             messageType=MPAS_LOG_ERR)
         ierr = 1
         return
      end if

      if (i == 1) call mpas_pool_add_dimension(tendPool, 'index_qv', 1)
      scalarsField % constituentNames(1) = 'tend_qv'
      call mpas_add_att(scalarsField % attLists(1) % attList, 'units', 'kg m^{-3} s^{-1}')
      call mpas_add_att(scalarsField % attLists(1) % attList, 'long_name', 'Tendency of water vapor mixing ratio')

      do j = 2, size(constituent_name)
         scalarsField % constituentNames(j) = 'tend_'//trim(constituent_name(mpas_from_ufs_cnst(j)))
      end do
   end do

   call mpas_pool_add_dimension(tendPool, 'moist_start', 1)
   call mpas_pool_add_dimension(tendPool, 'moist_end', num_moist)

 end subroutine ufs_mpas_define_scalars
 
 !> ########################################################################################
 !>
 !> \brief  Returns global mesh dimensions
 !> \author Michael Duda
 !> \date   22 August 2019
 !> \details
 !>  This routine returns on all tasks the number of global cells, edges,
 !>  vertices, maxEdges, vertical layers, and the maximum number of cells owned by any task.
 !>
 !> \update: Dustin Swales April 2025 - Modified for use in UWM
 !>
 !> ########################################################################################
 subroutine ufs_mpas_get_global_dims(nCellsGlobal, nEdgesGlobal, nVerticesGlobal, maxEdges,&
      nVertLevels, maxNCells)
   use mpas_pool_routines, only : mpas_pool_get_subpool, mpas_pool_get_dimension
   use mpas_derived_types, only : mpas_pool_type
   use mpas_dmpar,         only : mpas_dmpar_sum_int, mpas_dmpar_max_int
   use module_mpas_config, only : nCellsSolve, nEdgesSolve, nVerticesSolve

   integer, intent(out) :: nCellsGlobal
   integer, intent(out) :: nEdgesGlobal
   integer, intent(out) :: nVerticesGlobal
   integer, intent(out) :: maxEdges
   integer, intent(out) :: nVertLevels
   integer, intent(out) :: maxNCells

   integer, pointer :: maxEdgesLocal
   integer, pointer :: nVertLevelsLocal

   type (mpas_pool_type), pointer :: meshPool

   call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh', meshPool)
   call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
   call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve)
   call mpas_pool_get_dimension(meshPool, 'nVerticesSolve', nVerticesSolve)
   call mpas_pool_get_dimension(meshPool, 'maxEdges', maxEdgesLocal)
   call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevelsLocal)

   call mpas_dmpar_sum_int(domain_ptr % dminfo, nCellsSolve, nCellsGlobal)
   call mpas_dmpar_sum_int(domain_ptr % dminfo, nEdgesSolve, nEdgesGlobal)
   call mpas_dmpar_sum_int(domain_ptr % dminfo, nVerticesSolve, nVerticesGlobal)

   maxEdges = maxEdgesLocal
   nVertLevels = nVertLevelsLocal

   call mpas_dmpar_max_int(domain_ptr % dminfo, nCellsSolve, maxNCells)

 end subroutine ufs_mpas_get_global_dims

 !> ########################################################################################
 !>
 !> \brief  Returns global coordinate arrays
 !> \author Michael Duda
 !> \date   22 August 2019
 !> \details
 !>  This routine returns on all tasks arrays of latitude, longitude, and cell
 !>  area for all (global) cells.
 !>
 !>  It is assumed that latCellGlobal, lonCellGlobal, and areaCellGlobal have
 !>  been allocated by the caller with a size equal to the global number of
 !>  cells in the mesh.
 !>
 !> \update: Dustin Swales April 2025 - Modified for use in UWM
 !>
 !> ########################################################################################
 subroutine ufs_mpas_get_global_coords(latCellGlobal, lonCellGlobal, areaCellGlobal)
   use mpas_pool_routines, only : mpas_pool_get_subpool, mpas_pool_get_dimension, mpas_pool_get_array
   use mpas_derived_types, only : mpas_pool_type
   use mpas_kind_types,    only : RKIND
   use mpas_dmpar,         only : mpas_dmpar_sum_int, mpas_dmpar_max_real_array
   use mpp_mod,            only : FATAL, mpp_error
   use module_mpas_config, only : nCellsSolve, latCell, lonCell
   real (kind=RKIND), dimension(:), intent(out) :: latCellGlobal
   real (kind=RKIND), dimension(:), intent(out) :: lonCellGlobal
   real (kind=RKIND), dimension(:), intent(out) :: areaCellGlobal

   integer :: iCell

   integer, dimension(:), pointer :: indexToCellID

   type (mpas_pool_type), pointer :: meshPool
   integer :: nCellsGlobal,ierr

   real (kind=RKIND), dimension(:), pointer :: areaCell
   real (kind=RKIND), dimension(:), pointer :: temp

   character(len=*), parameter :: subname = 'ufs_mpas_subdriver::ufs_mpas_get_global_coords'


   call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh', meshPool)
   call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
   call mpas_pool_get_array(meshPool, 'indexToCellID', indexToCellID)
   call mpas_pool_get_array(meshPool, 'latCell', latCell)
   call mpas_pool_get_array(meshPool, 'lonCell', lonCell)
   call mpas_pool_get_array(meshPool, 'areaCell', areaCell)

   call mpas_dmpar_sum_int(domain_ptr % dminfo, nCellsSolve, nCellsGlobal)

   ! check: size(latCellGlobal) ?= nCellsGlobal
   allocate(temp(nCellsGlobal), stat=ierr)
   if( ierr /= 0 ) call mpp_error(FATAL,subname//':failed to allocate temp array')

   !
   ! latCellGlobal
   !
   temp(:) = -huge(temp(0))
   do iCell=1,nCellsSolve
      temp(indexToCellID(iCell)) = latCell(iCell)
   end do

   call mpas_dmpar_max_real_array(domain_ptr % dminfo, nCellsGlobal, temp, latCellGlobal)

   !
   ! lonCellGlobal
   !
   temp(:) = -huge(temp(0))
   do iCell=1,nCellsSolve
      temp(indexToCellID(iCell)) = lonCell(iCell)
   end do

   call mpas_dmpar_max_real_array(domain_ptr % dminfo, nCellsGlobal, temp, lonCellGlobal)

   !
   ! areaCellGlobal
   !
   temp(:) = -huge(temp(0))
   do iCell=1,nCellsSolve
      temp(indexToCellID(iCell)) = areaCell(iCell)
   end do

   call mpas_dmpar_max_real_array(domain_ptr % dminfo, nCellsGlobal, temp, areaCellGlobal)

   deallocate(temp)

 end subroutine ufs_mpas_get_global_coords
 
 ! ##########################################################################################
 ! \update: Dustin Swales April 2025 - Modified for use in UWM
 ! ##########################################################################################
 character(len=10) function date2yyyymmdd (date)
   ! Input arguments
   integer, intent(in) :: date

   ! Local workspace
   integer :: year    ! year of yyyy-mm-dd
   integer :: month   ! month of yyyy-mm-dd
   integer :: day     ! day of yyyy-mm-dd

   year  = date / 10000
   month = (date - year*10000) / 100
   day   = date - year*10000 - month*100

   write(date2yyyymmdd,80) year, month, day
80 format(i4.4,'-',i2.2,'-',i2.2)

 end function date2yyyymmdd
 ! #########################################################################################
 ! \update: Dustin Swales April 2025 - Modified for use in UWM
 ! #########################################################################################
 character(len=8) function sec2hms (seconds)
   ! Input arguments
   integer, intent(in) :: seconds

   ! Local workspace
   integer :: hours     ! hours of hh:mm:ss
   integer :: minutes   ! minutes of hh:mm:ss
   integer :: secs      ! seconds of hh:mm:ss

   hours   = seconds / 3600
   minutes = (seconds - hours*3600) / 60
   secs    = (seconds - hours*3600 - minutes*60)

   write(sec2hms,80) hours, minutes, secs
80 format(i2.2,':',i2.2,':',i2.2)

 end function sec2hms

 ! #########################################################################################
 ! \update: Dustin Swales April 2025 - Modified for use in UWM
 ! #########################################################################################
 character(len=10) function int2str(n)
   ! return default integer as a left justified string
   ! arguments
   integer, intent(in) :: n

   write(int2str,'(i0)') n
     
 end function int2str

  character(len=10) function log2str(n)
   ! return default integer as a left justified string
   ! arguments
   logical, intent(in) :: n

   if (n) then
      write(log2str,'(a4)') 'TRUE'
   else
      write(log2str,'(a4)') 'FALSE'
   endif

 end function log2str
 !> ########################################################################################
 !>
 !> subroutine dyn_mpas_exchange_halo
 !>
 !> summary: Update the halo layers of the named field.
 !> author: Michael Duda
 !> date: 16 January 2020
 !>
 !> Given a field name that is defined in MPAS registry, this subroutine updates
 !> the halo layers for that field.
 !> Ported and refactored for CAM-SIMA. (KCW, 2024-03-18)
 !>
 !> \update: Dustin Swales April 2025 - Modified for use in UWM
 !>
 !> ########################################################################################
 subroutine dyn_mpas_exchange_halo(field_name)
   ! Module(s) from MPAS.
   use mpas_derived_types, only : field1dinteger, field2dinteger, field3dinteger,           &
                                  field1dreal, field2dreal, field3dreal, field4dreal,       &
                                  field5dreal, mpas_pool_field_info_type, mpas_pool_integer,&
                                  mpas_pool_real
   use mpas_dmpar,         only : mpas_dmpar_exch_halo_field
   use mpas_pool_routines, only : mpas_pool_get_field, mpas_pool_get_field_info
   use mpp_mod,            only : FATAL, mpp_error
   use mpas_log,           only : mpas_log_write
   character(*), intent(in) :: field_name

   character(*), parameter :: subname = 'dyn_mpas_subdriver::dyn_mpas_exchange_halo'
   type(field1dinteger), pointer :: field_1d_integer
   type(field2dinteger), pointer :: field_2d_integer
   type(field3dinteger), pointer :: field_3d_integer
   type(field1dreal), pointer :: field_1d_real
   type(field2dreal), pointer :: field_2d_real
   type(field3dreal), pointer :: field_3d_real
   type(field4dreal), pointer :: field_4d_real
   type(field5dreal), pointer :: field_5d_real
   type(mpas_pool_field_info_type) :: mpas_pool_field_info

   call mpas_log_write(subname // ' entered')

   nullify(field_1d_integer)
   nullify(field_2d_integer)
   nullify(field_3d_integer)
   nullify(field_1d_real)
   nullify(field_2d_real)
   nullify(field_3d_real)
   nullify(field_4d_real)
   nullify(field_5d_real)

   call mpas_log_write('Inquiring field information for "' // trim(adjustl(field_name)) // '"')

   call mpas_pool_get_field_info(domain_ptr % blocklist % allfields, &
        trim(adjustl(field_name)), mpas_pool_field_info)

   if (mpas_pool_field_info % fieldtype == -1 .or. &
        mpas_pool_field_info % ndims == -1 .or. &
        mpas_pool_field_info % nhalolayers == -1) then
      call mpp_error(FATAL,subname//'Invalid field information for "' // trim(adjustl(field_name)) // '"')
   end if

   ! No halo layers to exchange. This field is not decomposed.
   if (mpas_pool_field_info % nhalolayers == 0) then
      call mpas_log_write('Skipping field "' // trim(adjustl(field_name)) // '" due to not decomposed')
      return
   end if

   call mpas_log_write('Exchanging halo layers for "' // trim(adjustl(field_name)) // '"')

   select case (mpas_pool_field_info % fieldtype)
   case (mpas_pool_integer)
      select case (mpas_pool_field_info % ndims)
      case (1)
         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(field_name)), field_1d_integer, timelevel=1)

         if (.not. associated(field_1d_integer)) then
            call mpp_error(FATAL,subname//'Failed to find field "' // trim(adjustl(field_name)) // '"')
         end if

         call mpas_dmpar_exch_halo_field(field_1d_integer)

         nullify(field_1d_integer)
      case (2)
         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(field_name)), field_2d_integer, timelevel=1)

         if (.not. associated(field_2d_integer)) then
            call mpp_error(FATAL,subname//'Failed to find field "' // trim(adjustl(field_name)) // '"')
         end if

         call mpas_dmpar_exch_halo_field(field_2d_integer)

         nullify(field_2d_integer)
      case (3)
         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(field_name)), field_3d_integer, timelevel=1)

         if (.not. associated(field_3d_integer)) then
            call mpp_error(FATAL,subname//'Failed to find field "' // trim(adjustl(field_name)) // '"')
         end if

         call mpas_dmpar_exch_halo_field(field_3d_integer)

         nullify(field_3d_integer)
      case default
         call mpp_error(FATAL,subname//'Unsupported field rank ' // stringify([mpas_pool_field_info % ndims]))
      end select
   case (mpas_pool_real)
      select case (mpas_pool_field_info % ndims)
      case (1)
         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(field_name)), field_1d_real, timelevel=1)

         if (.not. associated(field_1d_real)) then
            call mpp_error(FATAL,subname//'Failed to find field "' // trim(adjustl(field_name)) // '"')
         end if

         call mpas_dmpar_exch_halo_field(field_1d_real)

         nullify(field_1d_real)
      case (2)
         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(field_name)), field_2d_real, timelevel=1)

         if (.not. associated(field_2d_real)) then
            call mpp_error(FATAL,subname//'Failed to find field "' // trim(adjustl(field_name)) // '"')
         end if

         call mpas_dmpar_exch_halo_field(field_2d_real)

         nullify(field_2d_real)
      case (3)
         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(field_name)), field_3d_real, timelevel=1)

         if (.not. associated(field_3d_real)) then
            call mpp_error(FATAL,subname//'Failed to find field "' // trim(adjustl(field_name)) // '"')
         end if

         call mpas_dmpar_exch_halo_field(field_3d_real)

         nullify(field_3d_real)
      case (4)
         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(field_name)), field_4d_real, timelevel=1)

         if (.not. associated(field_4d_real)) then
            call mpp_error(FATAL,subname//'Failed to find field "' // trim(adjustl(field_name)) // '"')
         end if

         call mpas_dmpar_exch_halo_field(field_4d_real)

         nullify(field_4d_real)
      case (5)
         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(field_name)), field_5d_real, timelevel=1)

         if (.not. associated(field_5d_real)) then
            call mpp_error(FATAL,subname//'Failed to find field "' // trim(adjustl(field_name)) // '"')
         end if

         call mpas_dmpar_exch_halo_field(field_5d_real)

         nullify(field_5d_real)
      case default
         call mpp_error(FATAL,subname//'Unsupported field rank ' // stringify([mpas_pool_field_info % ndims]))
      end select
   case default
      call mpp_error(FATAL,subname//'Unsupported field type (Must be one of: integer, real)')
   end select

   call mpas_log_write(subname // ' completed')
 end subroutine dyn_mpas_exchange_halo
 
 !> ########################################################################################
 !> subroutine dyn_mpas_read_write_stream
 !>
 !> summary: Read or write an MPAS stream.
 !> author: Kuan-Chih Wang
 !> date: 2024-03-15
 !>
 !> In the context of MPAS, the concept of a "pool" resembles a group of
 !> (related) variables, while the concept of a "stream" resembles a file.
 !> This subroutine reads or writes an MPAS stream. It provides the mechanism
 !> for CAM-SIMA to input/output data to/from MPAS dynamical core.
 !> Analogous to the `{read,write}_stream` subroutines in MPAS stream manager.
 !>
 !> \update: Dustin Swales April 2025 - Modified for use in UWM
 !>
 !> ########################################################################################
 subroutine dyn_mpas_read_write_stream(stream_mode, stream_name, pio_file_desc, timeLevel, ierr)
   ! Module(s) from external libraries.
   use pio, only: file_desc_t
   use mpp_mod,             only : FATAL, mpp_error
   ! Module(s) from MPAS.
   use mpas_derived_types,  only : mpas_pool_type, mpas_stream_noerr, mpas_stream_type
   use mpas_io_streams,     only : mpas_closestream, mpas_readstream, mpas_writestream
   use mpas_pool_routines,  only : mpas_pool_destroy_pool
   use mpas_stream_manager, only : postread_reindex, prewrite_reindex, postwrite_reindex
   use mpas_log,            only : mpas_log_write
   use mpas_atm_halos,      only : exchange_halo_group
   use mpas_io_streams,     only : MPAS_STREAM_EXACT_TIME

   character(*), intent(in) :: stream_mode
   character(*), intent(in) :: stream_name
   type(file_desc_t), pointer, intent(in) :: pio_file_desc
   integer, intent(in) :: timeLevel
   integer, intent(out) :: ierr
   
   character(*), parameter :: subname = 'dyn_mpas_subdriver::dyn_mpas_read_write_stream'
   integer :: i
   type(mpas_pool_type), pointer :: mpas_pool
   type(mpas_stream_type), pointer :: mpas_stream
   type(var_info_type), allocatable :: var_info_list(:)

   ierr = 0
   call mpas_log_write('')

   
   nullify(mpas_pool)
   nullify(mpas_stream)
   call mpas_log_write( '---------------------------------------------------------------------')
   call mpas_log_write( 'Initializing stream "' // trim(adjustl(stream_name)) // '"')

   call dyn_mpas_init_stream_with_pool(mpas_pool, mpas_stream, pio_file_desc, stream_mode, stream_name, timeLevel)

   if (.not. associated(mpas_pool)) then
      call mpp_error(FATAL,subname//'Failed to initialize stream "' // trim(adjustl(stream_name)) // '"')
   end if

   if (.not. associated(mpas_stream)) then
      call mpp_error(FATAL,subname//'Failed to initialize stream "' // trim(adjustl(stream_name)) // '"')
   end if

   select case (trim(adjustl(stream_mode)))
   case ('r', 'read')
      call mpas_log_write('Reading stream "' // trim(adjustl(stream_name)) // '"')

      call mpas_readstream(mpas_stream, timeLevel, ierr=ierr)

      if (ierr /= mpas_stream_noerr) then
         call mpp_error(FATAL,subname//'Failed to read stream "' // trim(adjustl(stream_name)) // '"')
      end if

      ! Exchange halo layers because new data have just been read.
      var_info_list = parse_stream_name(stream_name)

      do i = 1, size(var_info_list)
         call dyn_mpas_exchange_halo(var_info_list(i) % name)
         if ( ierr /= 0 ) then
            call mpp_error(FATAL,subname//'Failed to exchange halo layers for group '//var_info_list(i) % name)
         end if
      end do

      ! For any connectivity arrays in this stream, convert global indexes to local indexes.
      call postread_reindex(domain_ptr % blocklist % allfields, domain_ptr % packages, &
           mpas_pool, mpas_pool)
   case ('w', 'write')
      call mpas_log_write('Writing stream "' // trim(adjustl(stream_name)) // '"')

      ! WARNING:
      ! The `{pre,post}write_reindex` subroutines are STATEFUL because they store information inside their module
      ! (i.e., module variables). They MUST be called in pairs, like below, to prevent undefined behaviors.
      ! For any connectivity arrays in this stream, temporarily convert local indexes to global indexes.
      call prewrite_reindex(domain_ptr % blocklist % allfields, domain_ptr % packages, &
           mpas_pool, mpas_pool)

      call mpas_writestream(mpas_stream, timeLevel, ierr=ierr)

      if (ierr /= mpas_stream_noerr) then
         call mpp_error(FATAL,subname//'Failed to write stream "' // trim(adjustl(stream_name)) // '"')
      end if

      ! For any connectivity arrays in this stream, reset global indexes back to local indexes.
      call postwrite_reindex(domain_ptr % blocklist % allfields, mpas_pool)
   case default
      call mpp_error(FATAL,subname//'Unsupported stream mode "' // trim(adjustl(stream_mode)) // '"')
   end select

   call mpas_log_write('Closing stream "' // trim(adjustl(stream_name)) // '"')
   call mpas_log_write( '---------------------------------------------------------------------')

   call mpas_closestream(mpas_stream, ierr=ierr)

   if (ierr /= mpas_stream_noerr) then
      call mpp_error(FATAL,subname//'Failed to close stream "' // trim(adjustl(stream_name)) // '"')
   end if

   ! Deallocate temporary pointers to avoid memory leaks.
   call mpas_pool_destroy_pool(mpas_pool)
   nullify(mpas_pool)

   deallocate(mpas_stream)
   nullify(mpas_stream)
   call mpas_log_write(subname // ' completed')
   
 end subroutine dyn_mpas_read_write_stream
 
 !> ########################################################################################
 !> subroutine dyn_mpas_init_stream_with_pool
 !>
 !> summary: Initialize an MPAS stream with an accompanying MPAS pool.
 !> author: Kuan-Chih Wang
 !> date: 2024-03-14
 !>
 !> In the context of MPAS, the concept of a "pool" resembles a group of
 !> (related) variables, while the concept of a "stream" resembles a file.
 !> This subroutine initializes an MPAS stream with an accompanying MPAS pool by
 !> adding variable and attribute information to them. After that, MPAS is ready
 !> to perform IO on them.
 !> Analogous to the `build_stream` and `mpas_stream_mgr_add_field`
 !> subroutines in MPAS stream manager.
 !>
 !> \update: Dustin Swales April 2025 - Modified for use in UWM
 !>
 !> ########################################################################################
 subroutine dyn_mpas_init_stream_with_pool(mpas_pool, mpas_stream, pio_file, stream_mode,  &
                                           stream_name, timeLevel)
   ! Module(s) from external libraries.
   use pio, only: file_desc_t, pio_file_is_open
   ! Module(s) from MPAS.
   use mpas_derived_types, only : field0dchar, field1dchar, field0dinteger, field1dinteger,&
                                  field2dinteger, field3dinteger, field0dreal, field1dreal,&
                                  field2dreal, field3dreal, field4dreal, field5dreal,      &
                                  mpas_io_native_precision, mpas_io_pnetcdf, mpas_io_read, &
                                  mpas_io_write, mpas_pool_type, mpas_stream_noerr,        &
                                  mpas_stream_type
   use mpas_io_streams,    only : mpas_createstream, mpas_streamaddfield
   use mpas_pool_routines, only : mpas_pool_add_config, mpas_pool_create_pool, mpas_pool_get_field
   use mpas_kind_types,    only : StrKIND, RKIND
   use mpp_mod,            only : FATAL, mpp_error
   use mpas_log,           only : mpas_log_write

   type(mpas_pool_type), pointer, intent(out) :: mpas_pool
   type(mpas_stream_type), pointer, intent(out) :: mpas_stream
   type(file_desc_t), pointer, intent(in) :: pio_file
   character(*), intent(in) :: stream_mode
   character(*), intent(in) :: stream_name
   integer, intent(in) :: timeLevel

   interface add_stream_attribute
      procedure :: add_stream_attribute_0d
      procedure :: add_stream_attribute_1d
   end interface add_stream_attribute

   character(*), parameter :: subname = 'dyn_mpas_subdriver::dyn_mpas_init_stream_with_pool'
   character(strkind) :: stream_filename
   integer :: i, ierr, stream_format
   !> Whether a variable is present on the file (i.e., `pio_file`).
   logical, allocatable :: var_is_present(:)
   !> Whether a variable is type, kind, and rank compatible with what MPAS expects on the file (i.e., `pio_file`).
   logical, allocatable :: var_is_tkr_compatible(:)
   type(field0dchar), pointer :: field_0d_char
   type(field1dchar), pointer :: field_1d_char
   type(field0dinteger), pointer :: field_0d_integer
   type(field1dinteger), pointer :: field_1d_integer
   type(field2dinteger), pointer :: field_2d_integer
   type(field3dinteger), pointer :: field_3d_integer
   type(field0dreal), pointer :: field_0d_real
   type(field1dreal), pointer :: field_1d_real
   type(field2dreal), pointer :: field_2d_real
   type(field3dreal), pointer :: field_3d_real
   type(field4dreal), pointer :: field_4d_real
   type(field5dreal), pointer :: field_5d_real
   type(var_info_type), allocatable :: var_info_list(:)

   call mpas_log_write(subname // ' entered')

   nullify(field_0d_char)
   nullify(field_1d_char)
   nullify(field_0d_integer)
   nullify(field_1d_integer)
   nullify(field_2d_integer)
   nullify(field_3d_integer)
   nullify(field_0d_real)
   nullify(field_1d_real)
   nullify(field_2d_real)
   nullify(field_3d_real)
   nullify(field_4d_real)
   nullify(field_5d_real)

   call mpas_pool_create_pool(mpas_pool)

   allocate(mpas_stream, stat=ierr)

   if (ierr /= 0) then
      call mpp_error(FATAL,subname//'Failed to allocate stream "' // trim(adjustl(stream_name)) // '"')
   end if

   ! Not actually used because a PIO file descriptor is directly supplied.
   stream_filename = 'external stream'
   stream_format = mpas_io_pnetcdf

   call mpas_log_write('Checking PIO file descriptor')

   if (.not. associated(pio_file)) then
      call mpp_error(FATAL,subname//'Invalid PIO file descriptor')
   end if

   if (.not. pio_file_is_open(pio_file)) then
      call mpp_error(FATAL,subname//'Invalid PIO file descriptor')
   end if

   select case (trim(adjustl(stream_mode)))
   case ('r', 'read')
      call mpas_log_write('Creating stream "' // trim(adjustl(stream_name)) // '" for reading')

      call mpas_createstream( &
           mpas_stream, domain_ptr % iocontext, stream_filename, stream_format, mpas_io_read,  &
           clobberrecords=.false., clobberfiles=.false., truncatefiles=.false., &
           precision=mpas_io_native_precision, pio_file_desc=pio_file, ierr=ierr)
   case ('w', 'write')
      call mpas_log_write('Creating stream "' // trim(adjustl(stream_name)) // '" for writing')

      call mpas_createstream( &
           mpas_stream, domain_ptr % iocontext, stream_filename, stream_format, mpas_io_write, &
           clobberrecords=.false., clobberfiles=.false., truncatefiles=.false., &
           precision=mpas_io_native_precision, pio_file_desc=pio_file, ierr=ierr)
   case default
      call mpp_error(FATAL,subname//'Unsupported stream mode "' // trim(adjustl(stream_mode)) // '"')
   end select

   if (ierr /= mpas_stream_noerr) then
      call mpp_error(FATAL,subname//'Failed to create stream "' // trim(adjustl(stream_name)) // '"')
   end if

   var_info_list = parse_stream_name(stream_name)

   ! Add variables contained in `var_info_list` to stream.
   do i = 1, size(var_info_list)
      call mpas_log_write('var_info_list(' // stringify([i]) // ') % name = ' // stringify([var_info_list(i) % name]))
      call mpas_log_write('var_info_list(' // stringify([i]) // ') % type = ' // stringify([var_info_list(i) % type]))
      call mpas_log_write('var_info_list(' // stringify([i]) // ') % rank = ' // stringify([var_info_list(i) % rank]))

      if (trim(adjustl(stream_mode)) == 'r' .or. trim(adjustl(stream_mode)) == 'read') then
         call dyn_mpas_check_variable_status(var_is_present, var_is_tkr_compatible, pio_file, var_info_list(i))

         ! Do not hard crash the model if a variable is missing and cannot be read.
         ! This can happen if users attempt to initialize/restart the model with data generated by
         ! older versions of MPAS. Print a debug message to let users decide if this is acceptable.
         if (.not. any(var_is_present)) then
            call mpas_log_write('Skipping variable "' // trim(adjustl(var_info_list(i) % name)) // '" due to not present')

            !cycle
         end if

         if (any(var_is_present .and. .not. var_is_tkr_compatible)) then
            call mpas_log_write('Skipping variable "' // trim(adjustl(var_info_list(i) % name)) // '" due to not TKR compatible')

            !cycle
         end if
      end if

      ! Add "<variable name>" to pool with the value of `1`.
      ! The existence of "<variable name>" in pool causes it to be considered for IO in MPAS.
      call mpas_pool_add_config(mpas_pool, trim(adjustl(var_info_list(i) % name)), 1)
      ! Add "<variable name>:packages" to pool with the value of an empty character string.
      ! This causes "<variable name>" to be always considered active for IO in MPAS.
      !call mpas_pool_add_config(mpas_pool, trim(adjustl(var_info_list(i) % name) // ':packages'), '')

      ! Add "<variable name>" to stream.
      call mpas_log_write('Adding variable "' // trim(adjustl(var_info_list(i) % name)) // &
           '" to stream "' // trim(adjustl(stream_name)) // '"')

      select case (trim(adjustl(var_info_list(i) % type)))
      case ('character')
         select case (var_info_list(i) % rank)
         case (0)
            call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
                 trim(adjustl(var_info_list(i) % name)), field_0d_char, timelevel=timeLevel)

            if (.not. associated(field_0d_char)) then
               call mpp_error(FATAL,subname//'Failed to find variable "' // trim(adjustl(var_info_list(i) % name)) // '"')
            end if

            call mpas_streamaddfield(mpas_stream, field_0d_char, ierr=ierr)

            nullify(field_0d_char)
         case (1)
            call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
                 trim(adjustl(var_info_list(i) % name)), field_1d_char, timelevel=timeLevel)

            if (.not. associated(field_1d_char)) then
               call mpp_error(FATAL,subname//'Failed to find variable "' // trim(adjustl(var_info_list(i) % name)) // '"')
            end if

            call mpas_streamaddfield(mpas_stream, field_1d_char, ierr=ierr)

            nullify(field_1d_char)
         case default
            call mpp_error(FATAL,subname//'Unsupported variable rank ' // stringify([var_info_list(i) % rank]) // &
                 ' for "' // trim(adjustl(var_info_list(i) % name)) // '"')
         end select
      case ('integer')
         select case (var_info_list(i) % rank)
         case (0)
            call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
                 trim(adjustl(var_info_list(i) % name)), field_0d_integer, timelevel=timeLevel)

            if (.not. associated(field_0d_integer)) then
               call mpp_error(FATAL,subname//'Failed to find variable "' // trim(adjustl(var_info_list(i) % name)) // '"')
            end if

            call mpas_streamaddfield(mpas_stream, field_0d_integer, ierr=ierr)

            nullify(field_0d_integer)
         case (1)
            call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
                 trim(adjustl(var_info_list(i) % name)), field_1d_integer, timelevel=timeLevel)

            if (.not. associated(field_1d_integer)) then
               call mpp_error(FATAL,subname//'Failed to find variable "' // trim(adjustl(var_info_list(i) % name)) // '"')
            end if

            call mpas_streamaddfield(mpas_stream, field_1d_integer, ierr=ierr)

            nullify(field_1d_integer)
         case (2)
            call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
                 trim(adjustl(var_info_list(i) % name)), field_2d_integer, timelevel=timeLevel)

            if (.not. associated(field_2d_integer)) then
               call mpp_error(FATAL,subname//'Failed to find variable "' // trim(adjustl(var_info_list(i) % name)) // '"')
            end if

            call mpas_streamaddfield(mpas_stream, field_2d_integer, ierr=ierr)

            nullify(field_2d_integer)
         case (3)
            call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
                 trim(adjustl(var_info_list(i) % name)), field_3d_integer, timelevel=timeLevel)

            if (.not. associated(field_3d_integer)) then
               call mpp_error(FATAL,subname//'Failed to find variable "' // trim(adjustl(var_info_list(i) % name)) // '"')
            end if

            call mpas_streamaddfield(mpas_stream, field_3d_integer, ierr=ierr)

            nullify(field_3d_integer)
         case default
            call mpp_error(FATAL,subname//'Unsupported variable rank ' // stringify([var_info_list(i) % rank]) // &
                 ' for "' // trim(adjustl(var_info_list(i) % name)) // '"')
         end select
      case ('real')
         select case (var_info_list(i) % rank)
         case (0)
            call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
                 trim(adjustl(var_info_list(i) % name)), field_0d_real, timelevel=timeLevel)

            if (.not. associated(field_0d_real)) then
               call mpp_error(FATAL,subname//'Failed to find variable "' // trim(adjustl(var_info_list(i) % name)) // '"')
            end if

            call mpas_streamaddfield(mpas_stream, field_0d_real, ierr=ierr)

            nullify(field_0d_real)
         case (1)
            call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
                 trim(adjustl(var_info_list(i) % name)), field_1d_real, timelevel=timeLevel)

            if (.not. associated(field_1d_real)) then
               call mpp_error(FATAL,subname//'Failed to find variable "' // trim(adjustl(var_info_list(i) % name)) // '"')
            end if

            call mpas_streamaddfield(mpas_stream, field_1d_real, ierr=ierr)

            nullify(field_1d_real)
         case (2)
            call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
                 trim(adjustl(var_info_list(i) % name)), field_2d_real, timelevel=timeLevel)

            if (.not. associated(field_2d_real)) then
               call mpp_error(FATAL,subname//'Failed to find variable "' // trim(adjustl(var_info_list(i) % name)) // '"')
            end if

            call mpas_streamaddfield(mpas_stream, field_2d_real, ierr=ierr)

            nullify(field_2d_real)
         case (3)
            call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
                 trim(adjustl(var_info_list(i) % name)), field_3d_real, timelevel=timeLevel)

            if (.not. associated(field_3d_real)) then
               call mpp_error(FATAL,subname//'Failed to find variable "' // trim(adjustl(var_info_list(i) % name)) // '"')
            end if

            call mpas_streamaddfield(mpas_stream, field_3d_real, ierr=ierr)

            nullify(field_3d_real)
         case (4)
            call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
                 trim(adjustl(var_info_list(i) % name)), field_4d_real, timelevel=timeLevel)

            if (.not. associated(field_4d_real)) then
               call mpp_error(FATAL,subname//'Failed to find variable "' // trim(adjustl(var_info_list(i) % name)) // '"')
            end if

            call mpas_streamaddfield(mpas_stream, field_4d_real, ierr=ierr)

            nullify(field_4d_real)
         case (5)
            call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
                 trim(adjustl(var_info_list(i) % name)), field_5d_real, timelevel=timeLevel)

            if (.not. associated(field_5d_real)) then
               call mpp_error(FATAL,subname//'Failed to find variable "' // trim(adjustl(var_info_list(i) % name)) // '"')
            end if

            call mpas_streamaddfield(mpas_stream, field_5d_real, ierr=ierr)

            nullify(field_5d_real)
         case default
            call mpp_error(FATAL,subname//'Unsupported variable rank ' // stringify([var_info_list(i) % rank]) // &
                 ' for "' // trim(adjustl(var_info_list(i) % name)) // '"')
         end select
      case default
         call mpp_error(FATAL,subname//'Unsupported variable type "' // trim(adjustl(var_info_list(i) % type)) // &
              '" for "' // trim(adjustl(var_info_list(i) % name)) // '"')
      end select

      if (ierr /= mpas_stream_noerr) then
         call mpp_error(FATAL,subname//'Failed to add variable "' // trim(adjustl(var_info_list(i) % name)) // &
              '" to stream "' // trim(adjustl(stream_name)) // '"')
      end if
   end do

   if (trim(adjustl(stream_mode)) == 'w' .or. trim(adjustl(stream_mode)) == 'write') then
      ! Add MPAS-specific attributes to stream.

      ! Attributes related to MPAS core (i.e., `core_type`).
      call add_stream_attribute('conventions', domain_ptr % core % conventions)
      call add_stream_attribute('core_name', domain_ptr % core % corename)
      call add_stream_attribute('git_version', domain_ptr % core % git_version)
      call add_stream_attribute('model_name', domain_ptr % core % modelname)
      call add_stream_attribute('source', domain_ptr % core % source)

      ! Attributes related to MPAS domain (i.e., `domain_type`).
      call add_stream_attribute('is_periodic', domain_ptr % is_periodic)
      call add_stream_attribute('mesh_spec', domain_ptr % mesh_spec)
      call add_stream_attribute('on_a_sphere', domain_ptr % on_a_sphere)
      call add_stream_attribute('parent_id',  domain_ptr % parent_id)
      call add_stream_attribute('sphere_radius', domain_ptr % sphere_radius)
      call add_stream_attribute('x_period',  domain_ptr % x_period)
      call add_stream_attribute('y_period',  domain_ptr % y_period)
   end if

   call mpas_log_write(subname // ' completed')
 contains
   !> Helper subroutine for adding a 0-d stream attribute by calling `mpas_writestreamatt` with error checking.
   !> (KCW, 2024-03-14)
   subroutine add_stream_attribute_0d(attribute_name, attribute_value)
     ! Module(s) from MPAS.
     use mpas_io_streams, only : mpas_writestreamatt
     use mpas_log,        only : mpas_log_write
     character(*), intent(in) :: attribute_name
     class(*), intent(in) :: attribute_value

     call mpas_log_write('Adding attribute "' // trim(adjustl(attribute_name)) // &
          '" to stream "' // trim(adjustl(stream_name)) // '"')

     select type (attribute_value)
     type is (character(*))
        call mpas_writestreamatt(mpas_stream, &
             trim(adjustl(attribute_name)), trim(adjustl(attribute_value)), syncval=.false., ierr=ierr)
     type is (integer)
        call mpas_writestreamatt(mpas_stream, &
             trim(adjustl(attribute_name)), attribute_value, syncval=.false., ierr=ierr)
     type is (logical)
        if (attribute_value) then
           ! Logical `.true.` becomes character string "YES".
           call mpas_writestreamatt(mpas_stream, &
                trim(adjustl(attribute_name)), 'YES', syncval=.false., ierr=ierr)
        else
           ! Logical `.false.` becomes character string "NO".
           call mpas_writestreamatt(mpas_stream, &
                trim(adjustl(attribute_name)), 'NO', syncval=.false., ierr=ierr)
        end if
     type is (real(rkind))
        call mpas_writestreamatt(mpas_stream, &
             trim(adjustl(attribute_name)), attribute_value, syncval=.false., ierr=ierr)
     class default
        call mpp_error(FATAL,subname//'Unsupported attribute type (Must be one of: character, integer, logical, real)')
     end select

     if (ierr /= mpas_stream_noerr) then
        call mpp_error(FATAL,subname//'Failed to add attribute "' // trim(adjustl(attribute_name)) // &
             '" to stream "' // trim(adjustl(stream_name)) // '"')
     end if
   end subroutine add_stream_attribute_0d

   !> Helper subroutine for adding a 1-d stream attribute by calling `mpas_writestreamatt` with error checking.
   !> (KCW, 2024-03-14)
   subroutine add_stream_attribute_1d(attribute_name, attribute_value)
     ! Module(s) from MPAS.
     use mpas_io_streams, only : mpas_writestreamatt
     use mpas_log,        only : mpas_log_write
     character(*), intent(in) :: attribute_name
     class(*), intent(in) :: attribute_value(:)

     call mpas_log_write('Adding attribute "' // trim(adjustl(attribute_name)) // &
          '" to stream "' // trim(adjustl(stream_name)) // '"')

     select type (attribute_value)
     type is (integer)
        call mpas_writestreamatt(mpas_stream, &
             trim(adjustl(attribute_name)), attribute_value, syncval=.false., ierr=ierr)
     type is (real(rkind))
        call mpas_writestreamatt(mpas_stream, &
             trim(adjustl(attribute_name)), attribute_value, syncval=.false., ierr=ierr)
     class default
        call mpp_error(FATAL,subname//'Unsupported attribute type (Must be one of: integer, real)')
     end select

     if (ierr /= mpas_stream_noerr) then
        call mpp_error(FATAL,subname//'Failed to add attribute "' // trim(adjustl(attribute_name)) // &
             '" to stream "' // trim(adjustl(stream_name)) // '"')
     end if
   end subroutine add_stream_attribute_1d
 end subroutine dyn_mpas_init_stream_with_pool
 
 !> ########################################################################################
 !>
 !> Parse a stream name, which consists of one or more stream name fragments, and return the
 !> corresponding variable information as a list of `var_info_type`. Multiple stream name
 !> fragments should be separated by "+" (i.e., a plus, meaning "addition"
 !> operation) or "-" (i.e., a minus, meaning "subtraction" operation).
 !> A stream name fragment can be a predefined stream name (e.g., "invariant", "input", etc.)
 !> or a single variable name. For example, a stream name of "invariant+input+restart" means
 !> the union of variables in the "invariant", "input", and "restart" streams.
 !> Duplicate variable information in the resulting list is discarded.
 !>
 !> (KCW, 2024-06-01)
 !>
 !> \update: Dustin Swales April 2025 - Modified for use in UWM
 !>
 !> ######################################################################################## 
 pure function parse_stream_name(stream_name) result(var_info_list)
   character(*), intent(in) :: stream_name
   type(var_info_type), allocatable :: var_info_list(:)

   character(*), parameter :: supported_stream_name_operator = '+-'
   character(1) :: stream_name_operator
   character(:), allocatable :: stream_name_fragment
   character(len(invariant_var_info_list % name)), allocatable :: var_name_list(:)
   integer :: i, j, n, offset
   type(var_info_type), allocatable :: var_info_list_buffer(:)

   n = len_trim(stream_name)

   if (n == 0) then
      ! Empty character string means empty list.
      var_info_list = parse_stream_name_fragment('')

      return
   end if

   i = scan(stream_name, supported_stream_name_operator)

   if (i == 0) then
      ! No operators are present in the stream name. It is just a single stream name fragment.
      stream_name_fragment = stream_name
      var_info_list = parse_stream_name_fragment(stream_name_fragment)

      return
   end if

   offset = 0
   var_info_list = parse_stream_name_fragment('')

   do while (.true.)
      ! Extract operator from the stream name.
      if (offset > 0) then
         stream_name_operator = stream_name(offset:offset)
      else
         stream_name_operator = '+'
      end if

      ! Extract stream name fragment from the stream name.
      if (i > 1) then
         stream_name_fragment = stream_name(offset + 1:offset + i - 1)
      else
         stream_name_fragment = ''
      end if
      
      ! Process the stream name fragment according to the operator.
      if (len_trim(stream_name_fragment) > 0) then
         var_info_list_buffer = parse_stream_name_fragment(stream_name_fragment)

         select case (stream_name_operator)
         case ('+')
            var_info_list = [var_info_list, var_info_list_buffer]
         case ('-')
            do j = 1, size(var_info_list_buffer)
               var_name_list = var_info_list % name
               var_info_list = pack(var_info_list, var_name_list /= var_info_list_buffer(j) % name)
            end do
         case default
            ! Do nothing for unknown operators. Should not happen at all.
         end select
      end if

      offset = offset + i

      ! Terminate loop when everything in the stream name has been processed.
      if (offset + 1 > n) then
         exit
      end if

      i = scan(stream_name(offset + 1:), supported_stream_name_operator)

      ! Run the loop one last time for the remaining stream name fragment.
      if (i == 0) then
         i = n - offset + 1
      end if
   end do

   ! Discard duplicate variable information by names.
   var_name_list = var_info_list % name
   var_info_list = var_info_list(index_unique(var_name_list))
 end function parse_stream_name

 !> ########################################################################################
 !>
 !> Parse a stream name fragment and return the corresponding variable information as a list
 !> of `var_info_type`.
 !> A stream name fragment can be a predefined stream name (e.g., "invariant", "input", etc.)
 !> or a single variable name.
 !>
 !> (KCW, 2024-06-01)
 !>
 !> \update: Dustin Swales April 2025 - Modified for use in UWM
 !>
 !> ########################################################################################
 pure function parse_stream_name_fragment(stream_name_fragment) result(var_info_list)
   character(*), intent(in) :: stream_name_fragment
   type(var_info_type), allocatable :: var_info_list(:)

   character(len(invariant_var_info_list % name)), allocatable :: var_name_list(:)
   type(var_info_type), allocatable :: var_info_list_buffer(:)

   select case (trim(adjustl(stream_name_fragment)))
   case ('')
      allocate(var_info_list(0))
   case ('invariant')
      allocate(var_info_list, source=invariant_var_info_list)
   case ('input')
      allocate(var_info_list, source=input_var_info_list)
   case ('restart')
      allocate(var_info_list, source=restart_var_info_list)
   case ('output')
      allocate(var_info_list, source=output_var_info_list)
   case ('lbc_in')
      allocate(var_info_list, source=lbc_in_var_info_list)
   case default
      allocate(var_info_list(0))

      var_name_list = invariant_var_info_list % name

      if (any(var_name_list == trim(adjustl(stream_name_fragment)))) then
         var_info_list_buffer = pack(invariant_var_info_list, var_name_list == trim(adjustl(stream_name_fragment)))
         var_info_list = [var_info_list, var_info_list_buffer]
      end if

      var_name_list = input_var_info_list % name

      if (any(var_name_list == trim(adjustl(stream_name_fragment)))) then
         var_info_list_buffer = pack(input_var_info_list, var_name_list == trim(adjustl(stream_name_fragment)))
         var_info_list = [var_info_list, var_info_list_buffer]
      end if

      var_name_list = restart_var_info_list % name

      if (any(var_name_list == trim(adjustl(stream_name_fragment)))) then
         var_info_list_buffer = pack(restart_var_info_list, var_name_list == trim(adjustl(stream_name_fragment)))
         var_info_list = [var_info_list, var_info_list_buffer]
      end if

      var_name_list = output_var_info_list % name

      if (any(var_name_list == trim(adjustl(stream_name_fragment)))) then
         var_info_list_buffer = pack(output_var_info_list, var_name_list == trim(adjustl(stream_name_fragment)))
         var_info_list = [var_info_list, var_info_list_buffer]
      end if

      var_name_list = lbc_in_var_info_list % name

      if (any(var_name_list == trim(adjustl(stream_name_fragment)))) then
         var_info_list_buffer = pack(lbc_in_var_info_list, var_name_list == trim(adjustl(stream_name_fragment)))
         var_info_list = [var_info_list, var_info_list_buffer]
      end if
   end select
 end function parse_stream_name_fragment

 !> ########################################################################################
 !>
 !> Return the index of unique elements in `array`, which can be any intrinsic data types,
 !> as an integer array.
 !> If `array` contains zero element or is of unsupported data types, an empty integer array
 !> is produced. For example, `index_unique([1, 2, 3, 1, 2, 3, 4, 5])` returns `[1, 2, 3, 7, 8]`.
 !>
 !> (KCW, 2024-03-22)
 !>
 !> \update: Dustin Swales April 2025 - Modified for use in UWM
 !>
 !> ########################################################################################
 pure function index_unique(array)
   use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64

   class(*), intent(in) :: array(:)
   integer, allocatable :: index_unique(:)

   character(:), allocatable :: array_c(:)
   integer :: i, n
   logical :: mask_unique(size(array))

   n = size(array)

   if (n == 0) then
      allocate(index_unique(0))

      return
   end if

   mask_unique = .false.

   select type (array)
   type is (character(*))
      ! Workaround for a bug in GNU Fortran >= 12. This is perhaps the manifestation of GCC Bugzilla Bug 100819.
      ! When a character string array is passed as the actual argument to an unlimited polymorphic dummy argument,
      ! its array index and length parameter are mishandled.
      allocate(character(len(array)) :: array_c(size(array)))

      array_c(:) = array(:)

      do i = 1, n
         if (.not. any(array_c(i) == array_c .and. mask_unique)) then
            mask_unique(i) = .true.
         end if
      end do
      deallocate(array_c)
   type is (integer(int32))
      do i = 1, n
         if (.not. any(array(i) == array .and. mask_unique)) then
            mask_unique(i) = .true.
         end if
      end do
   type is (integer(int64))
      do i = 1, n
         if (.not. any(array(i) == array .and. mask_unique)) then
            mask_unique(i) = .true.
         end if
      end do
   type is (logical)
      do i = 1, n
         if (.not. any((array(i) .eqv. array) .and. mask_unique)) then
            mask_unique(i) = .true.
         end if
      end do
   type is (real(real32))
      do i = 1, n
         if (.not. any(array(i) == array .and. mask_unique)) then
            mask_unique(i) = .true.
         end if
      end do
   type is (real(real64))
      do i = 1, n
         if (.not. any(array(i) == array .and. mask_unique)) then
            mask_unique(i) = .true.
         end if
      end do
   class default
      allocate(index_unique(0))

      return
   end select
 
   index_unique = pack([(i, i = 1, n)], mask_unique)
 end function index_unique
 !> ########################################################################################
 !> subroutine dyn_mpas_check_variable_status
 !>
 !> summary: Check and return variable status on the given file.
 !> author: Kuan-Chih Wang
 !> date: 2024-06-04
 !>
 !> On the given file (i.e., `pio_file`), this subroutine checks whether the
 !> given variable (i.e., `var_info`) is present, and whether it is "TKR"
 !> compatible with what MPAS expects. "TKR" means type, kind, and rank.
 !> This subroutine can handle both ordinary variables and variable arrays.
 !> They are indicated by the `var` and `var_array` elements, respectively,
 !> in MPAS registry. For an ordinary variable, the checks are performed on
 !> itself. Otherwise, for a variable array, the checks are performed on its
 !> constituent parts instead.
 !>
 !> \update: Dustin Swales April 2025 - Modified for use in UWM
 !>
 !> ########################################################################################
 subroutine dyn_mpas_check_variable_status(var_is_present, var_is_tkr_compatible, pio_file,&
                                           var_info)
   ! Module(s) from external libraries.
   use pio, only: file_desc_t, pio_file_is_open, pio_char, pio_int, pio_real, pio_double,  &
                  pio_inq_varid, pio_inq_varndims, pio_inq_vartype, pio_noerr
   ! Module(s) from MPAS.
   use mpas_derived_types, only : field0dchar, field1dchar, field0dinteger, field1dinteger,&
                                  field2dinteger, field3dinteger, field0dreal, field1dreal,&
                                  field2dreal, field3dreal, field4dreal, field5dreal
   use mpas_kind_types,    only : r4kind, r8kind
   use mpas_pool_routines, only : mpas_pool_get_field
   use mpas_log,           only : mpas_log_write
   use mpas_kind_types,    only : StrKIND, RKIND
   use mpp_mod,            only : FATAL, mpp_error

   logical, allocatable, intent(out) :: var_is_present(:)
   logical, allocatable, intent(out) :: var_is_tkr_compatible(:)
   type(file_desc_t), pointer, intent(in) :: pio_file
   type(var_info_type), intent(in) :: var_info

   character(*), parameter :: subname = 'dyn_mpas_subdriver::dyn_mpas_check_variable_status'
   character(strkind), allocatable :: var_name_list(:)
   integer :: i, ierr, varid, varndims, vartype
   type(field0dchar), pointer :: field_0d_char
   type(field1dchar), pointer :: field_1d_char
   type(field0dinteger), pointer :: field_0d_integer
   type(field1dinteger), pointer :: field_1d_integer
      type(field2dinteger), pointer :: field_2d_integer
   type(field3dinteger), pointer :: field_3d_integer
   type(field0dreal), pointer :: field_0d_real
   type(field1dreal), pointer :: field_1d_real
   type(field2dreal), pointer :: field_2d_real
   type(field3dreal), pointer :: field_3d_real
   type(field4dreal), pointer :: field_4d_real
   type(field5dreal), pointer :: field_5d_real

   call mpas_log_write(subname // ' entered')

   nullify(field_0d_char)
   nullify(field_1d_char)
   nullify(field_0d_integer)
   nullify(field_1d_integer)
   nullify(field_2d_integer)
   nullify(field_3d_integer)
   nullify(field_0d_real)
   nullify(field_1d_real)
   nullify(field_2d_real)
   nullify(field_3d_real)
   nullify(field_4d_real)
   nullify(field_5d_real)

   ! Extract a list of variable names to check on the file.
   ! For an ordinary variable, this list just contains its name.
   ! For a variable array, this list contains the names of its constituent parts.
   select case (trim(adjustl(var_info % type)))
   case ('character')
      select case (var_info % rank)
      case (0)
         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(var_info % name)), field_0d_char, timelevel=1)

         if (.not. associated(field_0d_char)) then
            call mpp_error(FATAL,subname//'Failed to find variable "' // trim(adjustl(var_info % name)))
         end if

         if (field_0d_char % isvararray .and. associated(field_0d_char % constituentnames)) then
            allocate(var_name_list(size(field_0d_char % constituentnames)), stat=ierr)

            if (ierr /= 0) then
               call mpp_error(FATAL,subname//'Failed to allocate var_name_list')
            end if

            var_name_list(:) = field_0d_char % constituentnames(:)
         end if

         nullify(field_0d_char)
      case (1)
         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(var_info % name)), field_1d_char, timelevel=1)

         if (.not. associated(field_1d_char)) then
            call mpp_error(FATAL,subname//'Failed to find variable "' // trim(adjustl(var_info % name)))
         end if

         if (field_1d_char % isvararray .and. associated(field_1d_char % constituentnames)) then
            allocate(var_name_list(size(field_1d_char % constituentnames)), stat=ierr)

            if (ierr /= 0) then
               call mpp_error(FATAL,subname//'Failed to allocate var_name_list')
            end if

            var_name_list(:) = field_1d_char % constituentnames(:)
         end if

         nullify(field_1d_char)
      case default
         call mpp_error(FATAL,subname//'Unsupported variable rank ' // stringify([var_info % rank]) // &
              ' for "' // trim(adjustl(var_info % name)) // '"')
      end select
   case ('integer')
      select case (var_info % rank)
      case (0)
         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(var_info % name)), field_0d_integer, timelevel=1)

         if (.not. associated(field_0d_integer)) then
            call mpp_error(FATAL,subname//'Failed to find variable "' // trim(adjustl(var_info % name)) // '"')
         end if

         if (field_0d_integer % isvararray .and. associated(field_0d_integer % constituentnames)) then
            allocate(var_name_list(size(field_0d_integer % constituentnames)), stat=ierr)

            if (ierr /= 0) then
               call mpp_error(FATAL,subname//'Failed to allocate var_name_list')
            end if

            var_name_list(:) = field_0d_integer % constituentnames(:)
         end if

         nullify(field_0d_integer)
      case (1)
         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(var_info % name)), field_1d_integer, timelevel=1)

         if (.not. associated(field_1d_integer)) then
            call mpp_error(FATAL,subname//'Failed to find variable "' // trim(adjustl(var_info % name)) // '"')
         end if

         if (field_1d_integer % isvararray .and. associated(field_1d_integer % constituentnames)) then
            allocate(var_name_list(size(field_1d_integer % constituentnames)), stat=ierr)

            if (ierr /= 0) then
               call mpp_error(FATAL,subname//'Failed to allocate var_name_list')
            end if

            var_name_list(:) = field_1d_integer % constituentnames(:)
         end if

         nullify(field_1d_integer)
      case (2)
         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(var_info % name)), field_2d_integer, timelevel=1)

         if (.not. associated(field_2d_integer)) then
            call mpp_error(FATAL,subname//'Failed to find variable "' // trim(adjustl(var_info % name)) // '"')
         end if

         if (field_2d_integer % isvararray .and. associated(field_2d_integer % constituentnames)) then
            allocate(var_name_list(size(field_2d_integer % constituentnames)), stat=ierr)

            if (ierr /= 0) then
               call mpp_error(FATAL,subname//'Failed to allocate var_name_list')
            end if

            var_name_list(:) = field_2d_integer % constituentnames(:)
         end if

         nullify(field_2d_integer)
      case (3)
         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(var_info % name)), field_3d_integer, timelevel=1)

         if (.not. associated(field_3d_integer)) then
            call mpp_error(FATAL,subname//'Failed to find variable "' // trim(adjustl(var_info % name)) // '"')
         end if

         if (field_3d_integer % isvararray .and. associated(field_3d_integer % constituentnames)) then
            allocate(var_name_list(size(field_3d_integer % constituentnames)), stat=ierr)

            if (ierr /= 0) then
               call mpp_error(FATAL,subname//'Failed to allocate var_name_list')
            end if

            var_name_list(:) = field_3d_integer % constituentnames(:)
         end if

         nullify(field_3d_integer)
      case default
         call mpp_error(FATAL,subname//'Unsupported variable rank ' // stringify([var_info % rank]) // &
              ' for "' // trim(adjustl(var_info % name)) // '"')
      end select
   case ('real')
      select case (var_info % rank)
      case (0)

         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(var_info % name)), field_0d_real, timelevel=1)

         if (.not. associated(field_0d_real)) then
            call mpp_error(FATAL,subname//'Failed to find variable "' // trim(adjustl(var_info % name)) // '"')
         end if

         if (field_0d_real % isvararray .and. associated(field_0d_real % constituentnames)) then
            allocate(var_name_list(size(field_0d_real % constituentnames)), stat=ierr)

            if (ierr /= 0) then
               call mpp_error(FATAL,subname//'Failed to allocate var_name_list')
            end if

            var_name_list(:) = field_0d_real % constituentnames(:)
         end if

         nullify(field_0d_real)
      case (1)

         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(var_info % name)), field_1d_real, timelevel=1)

         if (.not. associated(field_1d_real)) then
            call mpas_log_write(subname//'Failed to find variable "' // trim(adjustl(var_info % name)) // '"')
            call mpp_error(FATAL,subname//'Failed to find variable "' // trim(adjustl(var_info % name)) // '"')
         end if

         if (field_1d_real % isvararray .and. associated(field_1d_real % constituentnames)) then
            allocate(var_name_list(size(field_1d_real % constituentnames)), stat=ierr)

            if (ierr /= 0) then
               call mpp_error(FATAL,subname//'Failed to allocate var_name_list')
            end if

            var_name_list(:) = field_1d_real % constituentnames(:)
         end if

         nullify(field_1d_real)
      case (2)
         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(var_info % name)), field_2d_real, timelevel=1)

         if (.not. associated(field_2d_real)) then
            call mpp_error(FATAL,subname//'Failed to find variable "' // trim(adjustl(var_info % name)) // '"')
         end if

         if (field_2d_real % isvararray .and. associated(field_2d_real % constituentnames)) then
            allocate(var_name_list(size(field_2d_real % constituentnames)), stat=ierr)

            if (ierr /= 0) then
               call mpp_error(FATAL,subname//'Failed to allocate var_name_list')
            end if

            var_name_list(:) = field_2d_real % constituentnames(:)
         end if

         nullify(field_2d_real)
      case (3)
         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(var_info % name)), field_3d_real, timelevel=1)

         if (.not. associated(field_3d_real)) then
            call mpp_error(FATAL,subname//'Failed to find variable "' // trim(adjustl(var_info % name)) // '"')
         end if

         if (field_3d_real % isvararray .and. associated(field_3d_real % constituentnames)) then
            allocate(var_name_list(size(field_3d_real % constituentnames)), stat=ierr)

            if (ierr /= 0) then
               call mpp_error(FATAL,subname//'Failed to allocate var_name_list')
            end if

            var_name_list(:) = field_3d_real % constituentnames(:)
         end if

         nullify(field_3d_real)
      case (4)
         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(var_info % name)), field_4d_real, timelevel=1)

         if (.not. associated(field_4d_real)) then
            call mpp_error(FATAL,subname//'Failed to find variable "' // trim(adjustl(var_info % name)) // '"')
         end if

         if (field_4d_real % isvararray .and. associated(field_4d_real % constituentnames)) then
            allocate(var_name_list(size(field_4d_real % constituentnames)), stat=ierr)

            if (ierr /= 0) then
               call mpp_error(FATAL,subname//'Failed to allocate var_name_list')
            end if

            var_name_list(:) = field_4d_real % constituentnames(:)
         end if

         nullify(field_4d_real)
      case (5)
         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(var_info % name)), field_5d_real, timelevel=1)

         if (.not. associated(field_5d_real)) then
            call mpp_error(FATAL,subname//'Failed to find variable "' // trim(adjustl(var_info % name)) // '"')
         end if

         if (field_5d_real % isvararray .and. associated(field_5d_real % constituentnames)) then
            allocate(var_name_list(size(field_5d_real % constituentnames)), stat=ierr)

            if (ierr /= 0) then
               call mpp_error(FATAL,subname//'Failed to allocate var_name_list')
            end if

            var_name_list(:) = field_5d_real % constituentnames(:)
         end if

         nullify(field_5d_real)
      case default
         call mpp_error(FATAL,subname//'Unsupported variable rank ' // stringify([var_info % rank]) // &
              ' for "' // trim(adjustl(var_info % name)) // '"')
      end select
   case default
      call mpp_error(FATAL,subname//'Unsupported variable type "' // trim(adjustl(var_info % type)) // &
           '" for "' // trim(adjustl(var_info % name)) // '"')
   end select

   if (.not. allocated(var_name_list)) then
      allocate(var_name_list(1), stat=ierr)

      if (ierr /= 0) then
         call mpp_error(FATAL,subname//'Failed to allocate var_name_list')
      end if

      var_name_list(1) = var_info % name
   end if

   allocate(var_is_present(size(var_name_list)), stat=ierr)

   if (ierr /= 0) then
      call mpp_error(FATAL,subname//'Failed to allocate var_is_present')
   end if

   var_is_present(:) = .false.

   allocate(var_is_tkr_compatible(size(var_name_list)), stat=ierr)
   if (ierr /= 0) then
      call mpp_error(FATAL,subname//'Failed to allocate var_is_tkr_compatible')
   end if

   var_is_tkr_compatible(:) = .false.

   if (.not. associated(pio_file)) then
      return
   end if

   if (.not. pio_file_is_open(pio_file)) then
      return
   end if

   call mpas_log_write('Checking variable "' // trim(adjustl(var_info % name)) // &
        '" for presence and TKR compatibility')

   do i = 1, size(var_name_list)
      ! Check if the variable is present on the file.
      ierr = pio_inq_varid(pio_file, trim(adjustl(var_name_list(i))), varid)

      if (ierr /= pio_noerr) then
         cycle
      end if

      var_is_present(i) = .true.

      ! Check if the variable is "TK"R compatible between MPAS and the file.
      ierr = pio_inq_vartype(pio_file, varid, vartype)

      if (ierr /= pio_noerr) then
         cycle
      end if

      select case (trim(adjustl(var_info % type)))
      case ('character')
         if (vartype /= pio_char) then
            cycle
         end if
      case ('integer')
         if (vartype /= pio_int) then
            cycle
        end if
      case ('real')
         ! When MPAS dynamical core is compiled at single precision, pairing it with double precision input data
         ! is not allowed to prevent loss of precision.
         if (rkind == r4kind .and. vartype /= pio_real) then

            cycle
         end if

         ! When MPAS dynamical core is compiled at double precision, pairing it with single and double precision
         ! input data is allowed.
         if (rkind == r8kind .and. vartype /= pio_real .and. vartype /= pio_double) then

            cycle
         end if
      case default
         cycle
      end select

      ! Check if the variable is TK"R" compatible between MPAS and the file.
      ierr = pio_inq_varndims(pio_file, varid, varndims)

      if (ierr /= pio_noerr) then
         cycle
      end if

      if (varndims /= var_info % rank) then
         cycle
      end if

      var_is_tkr_compatible(i) = .true.
   end do

   call mpas_log_write('var_name_list = ' // stringify(var_name_list))
   call mpas_log_write('var_is_present = ' // stringify(var_is_present))
   call mpas_log_write('var_is_tkr_compatible = ' // stringify(var_is_tkr_compatible))

   call mpas_log_write(subname // ' completed')
 end subroutine dyn_mpas_check_variable_status


end module ufs_mpas_module
