!> ###########################################################################################
!> \file ufs_mpas_io.F90
!>
!> Routines from the subdrivers for MPAS-A and CAM-SIMA have been adopted/modified here for
!> use within the UFS Weather Model for input/output.
!> 
!> MPAS-A Subdriver:    MPAS-Model/src/driver/mpas_subdriver.F
!> CAM-CESM (external): src/dynamics/mpas/driver/cam_mpas_subdriver.F90
!>                      (https://github.com/ESCOMP/CAM/blob/cam_development/)
!> CAM-SIMA (external): src/dynamics/mpas/driver/dyn_mpas_subdriver.F90
!>                      (https://github.com/ESCOMP/CAM-SIMA/blob/development/)
!>
!>
!> ###########################################################################################
module ufs_mpas_io
  use mpas_derived_types,  only : core_type, domain_type, mpas_Clock_type
  use mpas_derived_types,  only : MPAS_Time_Type
  use mpas_kind_types,     only : StrKIND
  use mpas_log,            only : mpas_log_write
  use mpas_derived_types,  only : MPAS_LOG_CRIT
  use module_mpas_config,  only : pio_iotype, pio_stride, pio_numiotasks, pio_iodesc
  use module_mpas_config,  only : lbc_filename, pioid_lbc, pio_subsystem_lbc
  use module_mpas_config,  only : ic_filename,  pioid_ic,  pio_subsystem_ic
  use module_mpas_config,  only :               pioid_output, pio_subsystem_output
  use module_mpas_config,  only : TIMELEVEL_NOW
  use ufs_mpas_tools,      only : stringify
  implicit none

  !
  type(core_type),       pointer :: corelist   => null()
  type(domain_type),     pointer :: domain_ptr => null()
  type(mpas_Clock_type), pointer :: clock      => null()

  !
  integer :: out_file_index
  type (MPAS_Time_Type), allocatable :: mpas_output_times(:)

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
  !> This list corresponds to the "ugwp_oro_data_in" stream in MPAS registry.
  !> It consists of variables that are members of the "sfc_input" structure.
  !> #########################################################################################
  type(var_info_type), parameter :: ugwp_oro_data_var_info_list(*) = [ &
       var_info_type('var2dls'                         , 'real'      , 1), &
       var_info_type('conls'                           , 'real'      , 1), &
       var_info_type('oa1ls'                           , 'real'      , 1), &
       var_info_type('oa2ls'                           , 'real'      , 1), &
       var_info_type('oa3ls'                           , 'real'      , 1), &
       var_info_type('oa4ls'                           , 'real'      , 1), &
       var_info_type('ol1ls'                           , 'real'      , 1), &
       var_info_type('ol2ls'                           , 'real'      , 1), &
       var_info_type('ol3ls'                           , 'real'      , 1), &
       var_info_type('ol4ls'                           , 'real'      , 1), &
       var_info_type('var2dss'                         , 'real'      , 1), &
       var_info_type('conss'                           , 'real'      , 1), &
       var_info_type('oa1ss'                           , 'real'      , 1), &
       var_info_type('oa2ss'                           , 'real'      , 1), &
       var_info_type('oa3ss'                           , 'real'      , 1), &
       var_info_type('oa4ss'                           , 'real'      , 1), &
       var_info_type('ol1ss'                           , 'real'      , 1), &
       var_info_type('ol2ss'                           , 'real'      , 1), &
       var_info_type('ol3ss'                           , 'real'      , 1), &
       var_info_type('ol4ss'                           , 'real'      , 1)  &
       ]

  !> #########################################################################################
  !> This list corresponds to the "sfc_input" stream in MPAS registry.
  !> It consists of variables that are members of the "sfc_input" structure.
  !> Only variables needed to initialize the CCPP physics surface schemes are included.
  !> #########################################################################################
  type(var_info_type), parameter :: sfc_input_var_info_list(*) = [ &
       var_info_type('isltyp'                          , 'integer'   , 1), &
       var_info_type('ivgtyp'                          , 'integer'   , 1), &
       var_info_type('sfc_albbck'                      , 'real'      , 1), &
       var_info_type('skintemp'                        , 'real'      , 1), &
       var_info_type('snow'                            , 'real'      , 1), &
       var_info_type('snowc'                           , 'real'      , 1), &
       var_info_type('snowh'                           , 'real'      , 1), &
       var_info_type('sst'                             , 'real'      , 1), &
       var_info_type('tmn'                             , 'real'      , 1), &
       var_info_type('vegfra'                          , 'real'      , 1), &
       var_info_type('seaice'                          , 'real'      , 1), &
       var_info_type('xice'                            , 'real'      , 1), &
       var_info_type('xland'                           , 'real'      , 1), &
       var_info_type('dzs'                             , 'real'      , 2), &
       var_info_type('sh2o'                            , 'real'      , 2), &
       var_info_type('smois'                           , 'real'      , 2), &
       var_info_type('tslb'                            , 'real'      , 2), &
       var_info_type('ter'                             , 'real'      , 1), &
       var_info_type('landmask'                        , 'integer'   , 1), &
       var_info_type('mminlu'                          , 'character' , 0), &
       var_info_type('isice_lu'                        , 'integer'   , 0), &
       var_info_type('iswater_lu'                      , 'integer'   , 0), &
       var_info_type('shdmin'                          , 'real'      , 1), &
       var_info_type('shdmax'                          , 'real'      , 1), &
       var_info_type('snoalb'                          , 'real'      , 1), &
       var_info_type('greenfrac'                       , 'real'      , 2), &
       var_info_type('albedo12m'                       , 'real'      , 2), &
       var_info_type('soilcomp'                        , 'real'      , 2), &
       var_info_type('soilcl1'                         , 'real'      , 1), &
       var_info_type('soilcl2'                         , 'real'      , 1), &
       var_info_type('soilcl3'                         , 'real'      , 1), &
       var_info_type('soilcl4'                         , 'real'      , 1), &
       var_info_type('var2d'                           , 'real'      , 1), &
       var_info_type('con'                             , 'real'      , 1), &
       var_info_type('oa1'                             , 'real'      , 1), &
       var_info_type('oa2'                             , 'real'      , 1), &
       var_info_type('oa3'                             , 'real'      , 1), &
       var_info_type('oa4'                             , 'real'      , 1), &
       var_info_type('ol1'                             , 'real'      , 1), &
       var_info_type('ol2'                             , 'real'      , 1), &
       var_info_type('ol3'                             , 'real'      , 1), &
       var_info_type('ol4'                             , 'real'      , 1)  &
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
       var_info_type('Time'                            , 'real'      , 0), &
       var_info_type('initial_time'                    , 'character' , 0), &
       !var_info_type('divergence'                      , 'real'      , 2), &
       !var_info_type('pressure'                        , 'real'      , 2), &
       !var_info_type('relhum'                          , 'real'      , 2), &
       !var_info_type('rho'                             , 'real'      , 2), &
       var_info_type('scalars'                         , 'real'      , 3), &
       var_info_type('surface_pressure'                , 'real'      , 1), &
       var_info_type('theta'                           , 'real'      , 2), &
       !var_info_type('u'                               , 'real'      , 2), &
       var_info_type('uReconstructMeridional'          , 'real'      , 2), &
       var_info_type('uReconstructZonal'               , 'real'      , 2), &
       var_info_type('vorticity'                       , 'real'      , 2), &
       var_info_type('w'                               , 'real'      , 2) &
       !var_info_type('zz'                              , 'real'      , 2)  &
    ]
  
contains

  !> #########################################################################################
  !> Procedure to open MPAS IC file.
  !>
  !> #########################################################################################
  subroutine ufs_mpas_open_init()
    use pio, only : pio_openfile, pio_nowrite
    ! Arguments
    ! Locals
    integer :: ierr
    logical :: file_exists
    character(len=*), parameter :: subname = 'ufs_mpas_subdriver::ufs_mpas_open_init'

    ! Open MPAS Initial Condition file.
    INQUIRE(FILE=ic_filename, EXIST=file_exists)
    if (file_exists) then
       ierr = pio_openfile(pio_subsystem_ic, pioid_ic, pio_iotype, ic_filename, pio_nowrite)
       if (ierr /= 0) then
          call mpas_log_write(subname // " ERROR: Opening initial condition file", &
                              messageType=MPAS_LOG_CRIT)
       end if
    else
       call mpas_log_write(subname // " ERROR: Initial condition file could not be found", &
                           messageType=MPAS_LOG_CRIT)
    end if
    
  end subroutine ufs_mpas_open_init

  !> #########################################################################################
  !> Procedure to open MPAS Lateral Boundary Condition file.
  !>
  !> #########################################################################################
  subroutine ufs_mpas_open_lbc()
    use pio, only : pio_openfile, pio_nowrite
    ! Arguments
    ! Locals
    integer :: ierr
    logical :: file_exists
    character(len=*), parameter :: subname = 'ufs_mpas_subdriver::ufs_mpas_open_lbc'

    ! Open MPAS Initial Condition file.
    INQUIRE(FILE=lbc_filename, EXIST=file_exists)
    if (file_exists) then
       ierr = pio_openfile(pio_subsystem_lbc, pioid_lbc, pio_iotype, lbc_filename, pio_nowrite)
       if (ierr /= 0) then
          call mpas_log_write(subname // " ERROR: Opening lateral boundary condition file", &
                              messageType=MPAS_LOG_CRIT)
       end if
    else
       call mpas_log_write(subname // " ERROR: Lateral boundary condition file could not be found", &
                           messageType=MPAS_LOG_CRIT)
    end if

  end subroutine ufs_mpas_open_lbc

  !> #########################################################################################
  !> Procedure to create and write to MPAS stream
  !>
  !> #########################################################################################
  subroutine ufs_mpas_write(stream_name, timestamp, debug)
    ! PIO
    use pio, only : pio_openfile, pio_createfile, PIO_WRITE, PIO_CLOBBER
    ! MPAS
    use mpas_timekeeping, only : MPAS_NOW, MPAS_STREAM_EARLIEST_STRICTLY_AFTER
    ! Arguments
    character(len=*), intent(in) :: stream_name
    character(len=*), intent(in) :: timestamp
    logical,          intent(in) :: debug
    ! Locals
    character(len=*), parameter :: subname = 'ufs_mpas_subdriver::ufs_mpas_write'
    character(len=:), allocatable :: filename
    integer :: ierr
    type(var_info_type), allocatable :: output_var_info_list(:)
    integer :: timelevel, whence

    if (trim(stream_name) == "output") then
       filename = 'history.'//trim(timestamp)//'.nc'
    else if (trim(stream_name) == "restart") then
       filename = 'restart.'//trim(timestamp)//'.nc'
    else if (trim(stream_name) == "input") then
       filename = 'input.'//trim(timestamp)//'.nc'
    else
       call mpas_log_write(subname // " ERROR: Invalid stream_name to ufs_mpas_write: stream_name ="// &
                           trim(stream_name), messageType=MPAS_LOG_CRIT)
    end if

    if (debug) call mpas_log_write(subname // "entering ufs_mpas_write")
    if (debug) call mpas_log_write(subname // "creating "//trim(stream_name)//" stream file: "//trim(filename))
    ierr = pio_createfile(pio_subsystem_output, pioid_output, pio_iotype, trim(filename))
    if ( ierr /= 0 ) then
       call mpas_log_write(subname // " ERROR: pio_createfile failed", messageType=MPAS_LOG_CRIT)
    endif

    output_var_info_list = parse_stream_name_fragment('output')
    timelevel = TIMELEVEL_NOW
    whence = MPAS_NOW

    call dyn_mpas_read_write_stream(clock, "write", stream_name, pioid_output, &
                                    timeLevel=timelevel, whence=whence,        &
                                    nRecord=1, ierr=ierr, debug=debug)
    if ( ierr /= 0 ) then
       call mpas_log_write(subname // " ERROR: dyn_mpas_read_write_stream failed ", &
                           messageType=MPAS_LOG_CRIT)
    endif
    
    if (debug) call mpas_log_write(subname // "exiting ufs_mpas_write")
  end subroutine ufs_mpas_write

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
   if ( ierr /= 0 ) then
      call mpas_log_write(subname // " ERROR: failed to allocate temp array", messageType=MPAS_LOG_CRIT)
   endif

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
 subroutine dyn_mpas_exchange_halo(field_name, debug)
   ! Module(s) from MPAS.
   use mpas_derived_types, only : field1dinteger, field2dinteger, field3dinteger,           &
                                  field1dreal, field2dreal, field3dreal, field4dreal,       &
                                  field5dreal, mpas_pool_field_info_type, mpas_pool_integer,&
                                  mpas_pool_real
   use mpas_dmpar,         only : mpas_dmpar_exch_halo_field
   use mpas_pool_routines, only : mpas_pool_get_field, mpas_pool_get_field_info
   character(*), intent(in) :: field_name
   logical, intent(in)      :: debug

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

   if (debug) call mpas_log_write(subname // ' entered')

   nullify(field_1d_integer)
   nullify(field_2d_integer)
   nullify(field_3d_integer)
   nullify(field_1d_real)
   nullify(field_2d_real)
   nullify(field_3d_real)
   nullify(field_4d_real)
   nullify(field_5d_real)

   if (debug) call mpas_log_write(subname // 'Inquiring field information for "' // trim(adjustl(field_name)) // '"')

   call mpas_pool_get_field_info(domain_ptr % blocklist % allfields, &
        trim(adjustl(field_name)), mpas_pool_field_info)

   if (mpas_pool_field_info % fieldtype == -1 .or. &
        mpas_pool_field_info % ndims == -1 .or. &
        mpas_pool_field_info % nhalolayers == -1) then
      call mpas_log_write(subname // ' ERROR: Invalid field information for "' // trim(adjustl(field_name)) // '"', messageType=MPAS_LOG_CRIT)
   end if

   ! No halo layers to exchange. This field is not decomposed.
   if (mpas_pool_field_info % nhalolayers == 0) then
      call mpas_log_write(subname // ' WARNING: Skipping field "' // trim(adjustl(field_name)) // '" due to not decomposed')
      return
   end if

   if (debug) call mpas_log_write(subname // 'Exchanging halo layers for "' // trim(adjustl(field_name)) // '"')

   select case (mpas_pool_field_info % fieldtype)
   case (mpas_pool_integer)
      select case (mpas_pool_field_info % ndims)
      case (1)
         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(field_name)), field_1d_integer, timelevel=1)

         if (.not. associated(field_1d_integer)) then
            call mpas_log_write(subname // ' ERROR: Failed to find field "' // trim(adjustl(field_name)) // '"', messageType=MPAS_LOG_CRIT)
         end if

         call mpas_dmpar_exch_halo_field(field_1d_integer)

         nullify(field_1d_integer)
      case (2)
         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(field_name)), field_2d_integer, timelevel=1)

         if (.not. associated(field_2d_integer)) then
            call mpas_log_write(subname // ' ERROR: Failed to find field "' // trim(adjustl(field_name)) // '"', messageType=MPAS_LOG_CRIT)
         end if

         call mpas_dmpar_exch_halo_field(field_2d_integer)

         nullify(field_2d_integer)
      case (3)
         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(field_name)), field_3d_integer, timelevel=1)

         if (.not. associated(field_3d_integer)) then
            call mpas_log_write(subname // ' ERROR: Failed to find field "' // trim(adjustl(field_name)) // '"', messageType=MPAS_LOG_CRIT)
         end if

         call mpas_dmpar_exch_halo_field(field_3d_integer)

         nullify(field_3d_integer)
      case default
         call mpas_log_write(subname // ' ERROR: Unsupported field rank "' // trim(adjustl(field_name)) // '"', messageType=MPAS_LOG_CRIT)
      end select
   case (mpas_pool_real)
      select case (mpas_pool_field_info % ndims)
      case (1)
         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(field_name)), field_1d_real, timelevel=1)

         if (.not. associated(field_1d_real)) then
            call mpas_log_write(subname // ' ERROR: Failed to find field "' // trim(adjustl(field_name)) // '"', messageType=MPAS_LOG_CRIT)
         end if

         call mpas_dmpar_exch_halo_field(field_1d_real)

         nullify(field_1d_real)
      case (2)
         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(field_name)), field_2d_real, timelevel=1)

         if (.not. associated(field_2d_real)) then
            call mpas_log_write(subname // ' ERROR: Failed to find field "' // trim(adjustl(field_name)) // '"', messageType=MPAS_LOG_CRIT)
         end if
         call mpas_dmpar_exch_halo_field(field_2d_real)
         nullify(field_2d_real)
      case (3)
         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(field_name)), field_3d_real, timelevel=1)

         if (.not. associated(field_3d_real)) then
            call mpas_log_write(subname // ' ERROR: Failed to find field "' // trim(adjustl(field_name)) // '"', messageType=MPAS_LOG_CRIT)
         end if

         call mpas_dmpar_exch_halo_field(field_3d_real)

         nullify(field_3d_real)
      case (4)
         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(field_name)), field_4d_real, timelevel=1)

         if (.not. associated(field_4d_real)) then
            call mpas_log_write(subname // ' ERROR: Failed to find field "' // trim(adjustl(field_name)) // '"', messageType=MPAS_LOG_CRIT)
         end if

         call mpas_dmpar_exch_halo_field(field_4d_real)

         nullify(field_4d_real)
      case (5)
         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(field_name)), field_5d_real, timelevel=1)

         if (.not. associated(field_5d_real)) then
            call mpas_log_write(subname // ' ERROR: Failed to find field "' // trim(adjustl(field_name)) // '"', messageType=MPAS_LOG_CRIT)
         end if

         call mpas_dmpar_exch_halo_field(field_5d_real)

         nullify(field_5d_real)
      case default
         call mpas_log_write(subname // ' ERROR: Unsupported field rank "' // trim(adjustl(field_name)) // '"', messageType=MPAS_LOG_CRIT)
      end select
   case default
      call mpas_log_write(subname // ' ERROR: Unsupported field type (Must be one of: integer, real)', messageType=MPAS_LOG_CRIT)
   end select

   if (debug) call mpas_log_write(subname // ' completed')
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
 subroutine dyn_mpas_read_write_stream(clock, stream_mode, stream_name, pio_file_desc,     &
      timeLevel, when, whence, actualWhen, nRecord, ierr, debug)
   ! Module(s) from external libraries.
   use pio,                 only : file_desc_t
   ! Module(s) from MPAS.
   use mpas_derived_types,  only : mpas_pool_type, mpas_stream_noerr, mpas_stream_type
   use mpas_io_streams,     only : mpas_closestream, mpas_writestream
   use mpas_pool_routines,  only : mpas_pool_destroy_pool
   use mpas_stream_manager, only : postread_reindex, prewrite_reindex, postwrite_reindex
   use mpas_io_streams,     only : MPAS_STREAM_EXACT_TIME
   use mpas_timekeeping,    only : mpas_get_clock_time, MPAS_NOW
   ! Arguments
   type (mpas_clock_type), intent(in) :: clock
   character(*), intent(in) :: stream_mode
   character(*), intent(in) :: stream_name
   type(file_desc_t), pointer, intent(in) :: pio_file_desc
   integer, intent(in) :: timeLevel
   character (len=*), intent(in), optional :: when
   integer, intent(in), optional :: whence
   character (len=*), intent(out), optional :: actualWhen
   integer, intent(in) :: nRecord
   integer, intent(out) :: ierr
   logical, intent(in) :: debug
   ! Local variables
   character(*), parameter :: subname = 'dyn_mpas_subdriver::dyn_mpas_read_write_stream'
   integer :: i
   type(mpas_pool_type), pointer :: mpas_pool
   type(mpas_stream_type), pointer :: mpas_stream
   type(var_info_type), allocatable :: var_info_list(:)
   
   ierr = 0

   nullify(mpas_pool)
   nullify(mpas_stream)
   if (debug) call mpas_log_write(subname // 'Initializing stream "' // trim(adjustl(stream_name)) // '"')

   call dyn_mpas_init_stream_with_pool(mpas_pool, mpas_stream, pio_file_desc, stream_mode, stream_name, timeLevel, debug)

   if (.not. associated(mpas_pool)) then
      call mpas_log_write(subname // ' ERROR: Failed to initialize stream "' // trim(adjustl(stream_name)) // '"', messageType=MPAS_LOG_CRIT)
   end if

   if (.not. associated(mpas_stream)) then
      call mpas_log_write(subname // ' ERROR: Failed to initialize stream "' // trim(adjustl(stream_name)) // '"', messageType=MPAS_LOG_CRIT)
   end if

   select case (trim(adjustl(stream_mode)))
   case ('r', 'read')
      if (debug) call mpas_log_write(subname // 'Reading stream "' // trim(adjustl(stream_name)) // '"')

      call read_stream(mpas_stream, actualWhen, nRecord, ierr)

      if (ierr /= mpas_stream_noerr) then
         call mpas_log_write(subname // ' ERROR: Failed to initialize stream "' // trim(adjustl(stream_name)) // '"', messageType=MPAS_LOG_CRIT)
      end if

      ! Exchange halo layers because new data have just been read.
      var_info_list = parse_stream_name(stream_name)

      do i = 1, size(var_info_list)
         call dyn_mpas_exchange_halo(var_info_list(i) % name, debug)
         if ( ierr /= 0 ) then
            call mpas_log_write(subname // ' ERROR: Failed to exchange halo layers for group '//var_info_list(i) % name, messageType=MPAS_LOG_CRIT)
         end if
      end do

      ! For any connectivity arrays in this stream, convert global indexes to local indexes.
      call postread_reindex(domain_ptr % blocklist % allfields, domain_ptr % packages, &
           mpas_pool, mpas_pool)
   case ('w', 'write')
      if (debug) call mpas_log_write(subname // 'Writing stream "' // trim(adjustl(stream_name)) // '"')

      ! WARNING:
      ! The `{pre,post}write_reindex` subroutines are STATEFUL because they store information inside their module
      ! (i.e., module variables). They MUST be called in pairs, like below, to prevent undefined behaviors.
      ! For any connectivity arrays in this stream, temporarily convert local indexes to global indexes.
      call prewrite_reindex(domain_ptr % blocklist % allfields, domain_ptr % packages, &
           mpas_pool, mpas_pool)

      call mpas_writestream(mpas_stream, 1, ierr=ierr)

      if (ierr /= mpas_stream_noerr) then
         call mpas_log_write(subname // ' ERROR: Failed to write stream "' // trim(adjustl(stream_name)) // '"', messageType=MPAS_LOG_CRIT)
      end if

      ! For any connectivity arrays in this stream, reset global indexes back to local indexes.
      call postwrite_reindex(domain_ptr % blocklist % allfields, mpas_pool)
   case default
      call mpas_log_write(subname // ' ERROR: Unsupported stream mode "' // trim(adjustl(stream_mode)) // '"', messageType=MPAS_LOG_CRIT)
   end select

   if (debug) call mpas_log_write('Closing stream "' // trim(adjustl(stream_name)) // '"')
   call mpas_log_write( '---------------------------------------------------------------------')

   call mpas_closestream(mpas_stream, ierr=ierr)

   if (ierr /= mpas_stream_noerr) then
      call mpas_log_write(subname // ' ERROR: Failed to close stream "' // trim(adjustl(stream_name)) // '"', messageType=MPAS_LOG_CRIT)
   end if

   ! Deallocate temporary pointers to avoid memory leaks.
   call mpas_pool_destroy_pool(mpas_pool)
   nullify(mpas_pool)

   deallocate(mpas_stream)
   nullify(mpas_stream)
   if (debug) call mpas_log_write(subname // ' completed')
   
 end subroutine dyn_mpas_read_write_stream

 !> ########################################################################################
 !> subroutine read_stream
 !>
 !>
 !> ########################################################################################
 subroutine read_stream(stream, actualWhen, nRecord, ierr)
   use mpas_io_streams,     only : MPAS_readStream, MPAS_streamTime
   use mpas_derived_types,  only : mpas_pool_type, mpas_stream_noerr, mpas_stream_type

   type(mpas_stream_type), pointer, intent(inout) :: stream
   integer, intent(in) :: nRecord
   character (len=*), intent(out), optional :: actualWhen
   integer, intent(out) :: ierr

   call MPAS_readStream(stream, nRecord, ierr=ierr)
   if (present(actualWhen)) then
      call MPAS_streamTime(stream, nRecord, actualWhen, ierr=ierr)
   endif
   
 end subroutine read_stream
     
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
                                           stream_name, timeLevel, debug)
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

   type(mpas_pool_type), pointer, intent(out) :: mpas_pool
   type(mpas_stream_type), pointer, intent(out) :: mpas_stream
   type(file_desc_t), pointer, intent(in) :: pio_file
   character(*), intent(in) :: stream_mode
   character(*), intent(in) :: stream_name
   integer, intent(in) :: timeLevel
   logical, intent(in) :: debug

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

   if (debug) call mpas_log_write(subname // ' entered')

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
      call mpas_log_write(subname // ' ERROR: Failed to allocate stream "' // trim(adjustl(stream_name)) // '"', messageType=MPAS_LOG_CRIT)
   end if

   ! Not actually used because a PIO file descriptor is directly supplied.
   stream_filename = 'external stream'
   stream_format = mpas_io_pnetcdf

   if (debug) call mpas_log_write('Checking PIO file descriptor')

   if (.not. associated(pio_file)) then
      call mpas_log_write(subname // ' ERROR: Invalid PIO file descriptor', messageType=MPAS_LOG_CRIT)
   end if

   if (.not. pio_file_is_open(pio_file)) then
      call mpas_log_write(subname // ' ERROR: Invalid PIO file descriptor', messageType=MPAS_LOG_CRIT)
   end if

   select case (trim(adjustl(stream_mode)))
   case ('r', 'read')
      if (debug) call mpas_log_write('Creating stream "' // trim(adjustl(stream_name)) // '" for reading')

      call mpas_createstream( &
           mpas_stream, domain_ptr % iocontext, stream_filename, stream_format, mpas_io_read,  &
           clobberrecords=.false., clobberfiles=.false., truncatefiles=.false., &
           precision=mpas_io_native_precision, pio_file_desc=pio_file, ierr=ierr)
   case ('w', 'write')
      if (debug) call mpas_log_write('Creating stream "' // trim(adjustl(stream_name)) // '" for writing')

      call mpas_createstream( &
           mpas_stream, domain_ptr % iocontext, stream_filename, stream_format, mpas_io_write, &
           clobberrecords=.false., clobberfiles=.false., truncatefiles=.false., &
           precision=mpas_io_native_precision, pio_file_desc=pio_file, ierr=ierr)
   case default
      call mpas_log_write(subname // ' ERROR: Unsupported stream mode "' // trim(adjustl(stream_mode)) // '"', messageType=MPAS_LOG_CRIT)
   end select

   if (ierr /= mpas_stream_noerr) then
      call mpas_log_write(subname // ' ERROR: Failed to create stream "' // trim(adjustl(stream_name)) // '"', messageType=MPAS_LOG_CRIT)
   end if

   var_info_list = parse_stream_name(stream_name)

   ! Add variables contained in `var_info_list` to stream.
   do i = 1, size(var_info_list)
      if (debug) then
         call mpas_log_write('var_info_list(' // stringify([i]) // ') % name = ' // stringify([var_info_list(i) % name]))
         call mpas_log_write('var_info_list(' // stringify([i]) // ') % type = ' // stringify([var_info_list(i) % type]))
         call mpas_log_write('var_info_list(' // stringify([i]) // ') % rank = ' // stringify([var_info_list(i) % rank]))
      endif

      if (trim(adjustl(stream_mode)) == 'r' .or. trim(adjustl(stream_mode)) == 'read') then
         call dyn_mpas_check_variable_status(var_is_present, var_is_tkr_compatible, pio_file, var_info_list(i), debug)

         ! Do not hard crash the model if a variable is missing and cannot be read.
         ! This can happen if users attempt to initialize/restart the model with data generated by
         ! older versions of MPAS. Print a debug message to let users decide if this is acceptable.
         if (.not. any(var_is_present)) then
            if (debug) call mpas_log_write('Skipping variable "' // trim(adjustl(var_info_list(i) % name)) // '" due to not present')
            cycle
         end if

         if (any(var_is_present .and. .not. var_is_tkr_compatible)) then
            if (debug) call mpas_log_write('Skipping variable "' // trim(adjustl(var_info_list(i) % name)) // '" due to not TKR compatible')
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
      if (debug) call mpas_log_write('Adding variable "' // trim(adjustl(var_info_list(i) % name)) // &
                                     '" to stream "' // trim(adjustl(stream_name)) // '"')

      select case (trim(adjustl(var_info_list(i) % type)))
      case ('character')
         select case (var_info_list(i) % rank)
         case (0)
            call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
                 trim(adjustl(var_info_list(i) % name)), field_0d_char, timelevel=timeLevel)

            if (.not. associated(field_0d_char)) then
               call mpas_log_write(subname // ' ERROR: Failed to find variable "' // &
                                  trim(adjustl(var_info_list(i) % name)) // '"',     &
                                  messageType=MPAS_LOG_CRIT)
            end if

            call mpas_streamaddfield(mpas_stream, field_0d_char, ierr=ierr)

            nullify(field_0d_char)
         case (1)
            call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
                 trim(adjustl(var_info_list(i) % name)), field_1d_char, timelevel=timeLevel)

            if (.not. associated(field_1d_char)) then
               call mpas_log_write(subname // ' ERROR: Failed to find variable "' // &
                                   trim(adjustl(var_info_list(i) % name)) // '"',    &
                                   messageType=MPAS_LOG_CRIT)
            end if

            call mpas_streamaddfield(mpas_stream, field_1d_char, ierr=ierr)

            nullify(field_1d_char)
         case default
            call mpas_log_write(subname // ' ERROR: Unsupported variable rank ' //         &
                                stringify([var_info_list(i) % rank]) //                    &
                                ' for "' // trim(adjustl(var_info_list(i) % name)) // '"', &
                                messageType=MPAS_LOG_CRIT)
         end select
      case ('integer')
         select case (var_info_list(i) % rank)
         case (0)
            call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
                 trim(adjustl(var_info_list(i) % name)), field_0d_integer, timelevel=timeLevel)

            if (.not. associated(field_0d_integer)) then
               call mpas_log_write(subname // ' ERROR: Failed to find variable "' //  &
                                   trim(adjustl(var_info_list(i) % name)) // '"',     &
                                   messageType=MPAS_LOG_CRIT)
            end if

            call mpas_streamaddfield(mpas_stream, field_0d_integer, ierr=ierr)

            nullify(field_0d_integer)
         case (1)
            call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
                 trim(adjustl(var_info_list(i) % name)), field_1d_integer, timelevel=timeLevel)

            if (.not. associated(field_1d_integer)) then
               call mpas_log_write(subname // ' ERROR: Failed to find variable "' // &
                                   trim(adjustl(var_info_list(i) % name)) // '"',    &
                                   messageType=MPAS_LOG_CRIT)
            end if

            call mpas_streamaddfield(mpas_stream, field_1d_integer, ierr=ierr)

            nullify(field_1d_integer)
         case (2)
            call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
                 trim(adjustl(var_info_list(i) % name)), field_2d_integer, timelevel=timeLevel)

            if (.not. associated(field_2d_integer)) then
               call mpas_log_write(subname // ' ERROR: Failed to find variable "' // &
                                   trim(adjustl(var_info_list(i) % name)) // '"',    &
                                   messageType=MPAS_LOG_CRIT)
            end if

            call mpas_streamaddfield(mpas_stream, field_2d_integer, ierr=ierr)

            nullify(field_2d_integer)
         case (3)
            call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
                 trim(adjustl(var_info_list(i) % name)), field_3d_integer, timelevel=timeLevel)

            if (.not. associated(field_3d_integer)) then
               call mpas_log_write(subname // ' ERROR: Failed to find variable "' // &
                                   trim(adjustl(var_info_list(i) % name)) // '"',    &
                                   messageType=MPAS_LOG_CRIT)
            end if

            call mpas_streamaddfield(mpas_stream, field_3d_integer, ierr=ierr)

            nullify(field_3d_integer)
         case default
            call mpas_log_write(subname // ' ERROR: Unsupported variable rank ' //         &
                                stringify([var_info_list(i) % rank]) //                    &
                                ' for "' // trim(adjustl(var_info_list(i) % name)) // '"', &
                                messageType=MPAS_LOG_CRIT)
         end select
      case ('real')
         select case (var_info_list(i) % rank)
         case (0)
            call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
                 trim(adjustl(var_info_list(i) % name)), field_0d_real, timelevel=timeLevel)

            if (.not. associated(field_0d_real)) then
               call mpas_log_write(subname // ' ERROR: Failed to find variable "' //  &
                                   trim(adjustl(var_info_list(i) % name)) // '"',     &
                                   messageType=MPAS_LOG_CRIT)
            end if

            call mpas_streamaddfield(mpas_stream, field_0d_real, ierr=ierr)

            nullify(field_0d_real)
         case (1)
            call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
                 trim(adjustl(var_info_list(i) % name)), field_1d_real, timelevel=timeLevel)

            if (.not. associated(field_1d_real)) then
               call mpas_log_write(subname // ' ERROR: Failed to find variable "' // &
                                   trim(adjustl(var_info_list(i) % name)) // '"',    &
                                   messageType=MPAS_LOG_CRIT)
            end if

            call mpas_streamaddfield(mpas_stream, field_1d_real, ierr=ierr)

            nullify(field_1d_real)
         case (2)
            call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
                 trim(adjustl(var_info_list(i) % name)), field_2d_real, timelevel=timeLevel)
            if (.not. associated(field_2d_real)) then
               call mpas_log_write(subname // ' ERROR: Failed to find variable "' // &
                                   trim(adjustl(var_info_list(i) % name)) // '"',    &
                                   messageType=MPAS_LOG_CRIT)
            end if
            call mpas_streamaddfield(mpas_stream, field_2d_real, ierr=ierr)

            nullify(field_2d_real)
         case (3)
            call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
                 trim(adjustl(var_info_list(i) % name)), field_3d_real, timelevel=timeLevel)

            if (.not. associated(field_3d_real)) then
               call mpas_log_write(subname // ' ERROR: Failed to find variable "' // &
                                   trim(adjustl(var_info_list(i) % name)) // '"',    &
                                   messageType=MPAS_LOG_CRIT)
            end if
            call mpas_streamaddfield(mpas_stream, field_3d_real, ierr=ierr)

            nullify(field_3d_real)
         case (4)
            call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
                 trim(adjustl(var_info_list(i) % name)), field_4d_real, timelevel=timeLevel)

            if (.not. associated(field_4d_real)) then
               call mpas_log_write(subname // ' ERROR: Failed to find variable "' // &
                                   trim(adjustl(var_info_list(i) % name)) // '"',    &
                                   messageType=MPAS_LOG_CRIT)
            end if

            call mpas_streamaddfield(mpas_stream, field_4d_real, ierr=ierr)

            nullify(field_4d_real)
         case (5)
            call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
                 trim(adjustl(var_info_list(i) % name)), field_5d_real, timelevel=timeLevel)

            if (.not. associated(field_5d_real)) then
               call mpas_log_write(subname // ' ERROR: Failed to find variable "' // &
                                   trim(adjustl(var_info_list(i) % name)) // '"',    &
                                   messageType=MPAS_LOG_CRIT)
            end if

            call mpas_streamaddfield(mpas_stream, field_5d_real, ierr=ierr)

            nullify(field_5d_real)
         case default
            call mpas_log_write(subname // ' ERROR: Unsupported variable rank ' //         &
                                stringify([var_info_list(i) % rank]) //                    &
                                ' for "' // trim(adjustl(var_info_list(i) % name)) // '"', &
                                messageType=MPAS_LOG_CRIT)
         end select
      case default
         call mpas_log_write(subname // ' ERROR: Unsupported variable type "' //         &
                             trim(adjustl(var_info_list(i) % type)) //                   &
                             '" for "' // trim(adjustl(var_info_list(i) % name)) // '"', &
                             messageType=MPAS_LOG_CRIT)
      end select

      if (ierr /= mpas_stream_noerr) then
         call mpas_log_write(subname // ' ERROR: Failed to add variable "' //      &
                             trim(adjustl(var_info_list(i) % name)) //             &
                             '" to stream "' // trim(adjustl(stream_name)) // '"', &
                             messageType=MPAS_LOG_CRIT)
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

   if (debug) call mpas_log_write(subname // ' completed')
 contains
   !> Helper subroutine for adding a 0-d stream attribute by calling `mpas_writestreamatt` with error checking.
   !> (KCW, 2024-03-14)
   subroutine add_stream_attribute_0d(attribute_name, attribute_value)
     ! Module(s) from MPAS.
     use mpas_io_streams, only : mpas_writestreamatt
     character(*), intent(in) :: attribute_name
     class(*), intent(in) :: attribute_value

     if (debug) call mpas_log_write('Adding attribute "' // trim(adjustl(attribute_name)) // &
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
        call mpas_log_write(subname // ' ERROR: Unsupported attribute type (Must be one of: character, integer, logical, real)', &
                            messageType=MPAS_LOG_CRIT)
     end select

     if (ierr /= mpas_stream_noerr) then
        call mpas_log_write(subname // ' ERROR: Failed to add attribute "' //     &
                            trim(adjustl(attribute_name)) //                      &
                            '" to stream "' // trim(adjustl(stream_name)) // '"', &
                            messageType=MPAS_LOG_CRIT)
     end if
   end subroutine add_stream_attribute_0d

   !> Helper subroutine for adding a 1-d stream attribute by calling `mpas_writestreamatt` with error checking.
   !> (KCW, 2024-03-14)
   subroutine add_stream_attribute_1d(attribute_name, attribute_value)
     ! Module(s) from MPAS.
     use mpas_io_streams, only : mpas_writestreamatt
     character(*), intent(in) :: attribute_name
     class(*), intent(in) :: attribute_value(:)

     if (debug) call mpas_log_write(subname // 'Adding attribute "' // trim(adjustl(attribute_name)) // &
                                    '" to stream "' // trim(adjustl(stream_name)) // '"')

     select type (attribute_value)
     type is (integer)
        call mpas_writestreamatt(mpas_stream, &
             trim(adjustl(attribute_name)), attribute_value, syncval=.false., ierr=ierr)
     type is (real(rkind))
        call mpas_writestreamatt(mpas_stream, &
             trim(adjustl(attribute_name)), attribute_value, syncval=.false., ierr=ierr)
     class default
        call mpas_log_write(subname // ' ERROR: Unsupported attribute type (Must be one of: integer, real)',&
                            messageType=MPAS_LOG_CRIT)
     end select

     if (ierr /= mpas_stream_noerr) then
        call mpas_log_write(subname // ' ERROR: Failed to add attribute "' //     &
                            trim(adjustl(attribute_name)) //                      &
                            '" to stream "' // trim(adjustl(stream_name)) // '"', &
                            messageType=MPAS_LOG_CRIT)
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
   case ('sfc_input')
      allocate(var_info_list, source=sfc_input_var_info_list)
   case ('ugwp_oro_data')
      allocate(var_info_list, source=ugwp_oro_data_var_info_list)
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
                                           var_info, debug)
   ! Module(s) from external libraries.
   use pio, only: file_desc_t, pio_file_is_open, pio_char, pio_int, pio_real, pio_double,  &
                  pio_inq_varid, pio_inq_varndims, pio_inq_vartype, pio_noerr
   ! Module(s) from MPAS.
   use mpas_derived_types, only : field0dchar, field1dchar, field0dinteger, field1dinteger,&
                                  field2dinteger, field3dinteger, field0dreal, field1dreal,&
                                  field2dreal, field3dreal, field4dreal, field5dreal
   use mpas_kind_types,    only : r4kind, r8kind
   use mpas_pool_routines, only : mpas_pool_get_field
   use mpas_kind_types,    only : StrKIND, RKIND

   logical, allocatable, intent(out) :: var_is_present(:)
   logical, allocatable, intent(out) :: var_is_tkr_compatible(:)
   type(file_desc_t), pointer, intent(in) :: pio_file
   type(var_info_type), intent(in) :: var_info
   logical, intent(in) :: debug

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

   if (debug) call mpas_log_write(subname // ' entered')

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
            call mpas_log_write(subname // ' ERROR: Failed to find variable "' // &
                                trim(adjustl(var_info % name)) // '"',            &
                                messageType=MPAS_LOG_CRIT)
         end if

         if (field_0d_char % isvararray .and. associated(field_0d_char % constituentnames)) then
            allocate(var_name_list(size(field_0d_char % constituentnames)), stat=ierr)

            if (ierr /= 0) then
               call mpas_log_write(subname // ' ERROR: Failed to allocate var_name_list', &
                                   messageType=MPAS_LOG_CRIT)
            end if

            var_name_list(:) = field_0d_char % constituentnames(:)
         end if

         nullify(field_0d_char)
      case (1)
         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(var_info % name)), field_1d_char, timelevel=1)

         if (.not. associated(field_1d_char)) then
            call mpas_log_write(subname // ' ERROR: Failed to find variable "' // &
                                trim(adjustl(var_info % name)) // '"',            &
                                messageType=MPAS_LOG_CRIT)
         end if

         if (field_1d_char % isvararray .and. associated(field_1d_char % constituentnames)) then
            allocate(var_name_list(size(field_1d_char % constituentnames)), stat=ierr)

            if (ierr /= 0) then
               call mpas_log_write(subname // ' ERROR: Failed to allocate var_name_list', &
                                   messageType=MPAS_LOG_CRIT)
            end if

            var_name_list(:) = field_1d_char % constituentnames(:)
         end if

         nullify(field_1d_char)
      case default
         call mpas_log_write(subname // ' ERROR: Unsupported variable rank ' // &
                             stringify([var_info % rank]) //                    &
                             ' for "' // trim(adjustl(var_info % name)) // '"', &
                             messageType=MPAS_LOG_CRIT)
      end select
   case ('integer')
      select case (var_info % rank)
      case (0)
         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(var_info % name)), field_0d_integer, timelevel=1)

         if (.not. associated(field_0d_integer)) then
            call mpas_log_write(subname // ' ERROR: Failed to find variable "' // &
                                trim(adjustl(var_info % name)) // '"',            &
                                messageType=MPAS_LOG_CRIT)
         end if

         if (field_0d_integer % isvararray .and. associated(field_0d_integer % constituentnames)) then
            allocate(var_name_list(size(field_0d_integer % constituentnames)), stat=ierr)

            if (ierr /= 0) then
               call mpas_log_write(subname // ' ERROR: Failed to allocate var_name_list', &
                                   messageType=MPAS_LOG_CRIT)
            end if

            var_name_list(:) = field_0d_integer % constituentnames(:)
         end if

         nullify(field_0d_integer)
      case (1)
         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(var_info % name)), field_1d_integer, timelevel=1)

         if (.not. associated(field_1d_integer)) then
            call mpas_log_write(subname // ' ERROR: Failed to find variable "' // &
                                trim(adjustl(var_info % name)) // '"',            &
                                messageType=MPAS_LOG_CRIT)
         end if

         if (field_1d_integer % isvararray .and. associated(field_1d_integer % constituentnames)) then
            allocate(var_name_list(size(field_1d_integer % constituentnames)), stat=ierr)

            if (ierr /= 0) then
               call mpas_log_write(subname // ' ERROR: Failed to allocate var_name_list', &
                                   messageType=MPAS_LOG_CRIT)
            end if

            var_name_list(:) = field_1d_integer % constituentnames(:)
         end if

         nullify(field_1d_integer)
      case (2)
         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(var_info % name)), field_2d_integer, timelevel=1)

         if (.not. associated(field_2d_integer)) then
            call mpas_log_write(subname // ' ERROR: Failed to find variable "' // &
                                trim(adjustl(var_info % name)) // '"',            &
                                messageType=MPAS_LOG_CRIT)
         end if

         if (field_2d_integer % isvararray .and. associated(field_2d_integer % constituentnames)) then
            allocate(var_name_list(size(field_2d_integer % constituentnames)), stat=ierr)

            if (ierr /= 0) then
               call mpas_log_write(subname // ' ERROR: Failed to allocate var_name_list', &
                                   messageType=MPAS_LOG_CRIT)
            end if

            var_name_list(:) = field_2d_integer % constituentnames(:)
         end if

         nullify(field_2d_integer)
      case (3)
         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(var_info % name)), field_3d_integer, timelevel=1)

         if (.not. associated(field_3d_integer)) then
            call mpas_log_write(subname // ' ERROR: Failed to find variable "' // &
                                trim(adjustl(var_info % name)) // '"',            &
                                messageType=MPAS_LOG_CRIT)
         end if

         if (field_3d_integer % isvararray .and. associated(field_3d_integer % constituentnames)) then
            allocate(var_name_list(size(field_3d_integer % constituentnames)), stat=ierr)

            if (ierr /= 0) then
               call mpas_log_write(subname // ' ERROR: Failed to allocate var_name_list', &
                                   messageType=MPAS_LOG_CRIT)
            end if

            var_name_list(:) = field_3d_integer % constituentnames(:)
         end if

         nullify(field_3d_integer)
      case default
         call mpas_log_write(subname // ' ERROR: Unsupported variable rank ' // &
                             stringify([var_info % rank]) //                    &
                             ' for "' // trim(adjustl(var_info % name)) // '"', &
                             messageType=MPAS_LOG_CRIT)
      end select
   case ('real')
      select case (var_info % rank)
      case (0)

         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(var_info % name)), field_0d_real, timelevel=1)

         if (.not. associated(field_0d_real)) then
            call mpas_log_write(subname // ' ERROR: Failed to find variable "' // &
                                trim(adjustl(var_info % name)) // '"',            &
                                messageType=MPAS_LOG_CRIT)
         end if

         if (field_0d_real % isvararray .and. associated(field_0d_real % constituentnames)) then
            allocate(var_name_list(size(field_0d_real % constituentnames)), stat=ierr)

            if (ierr /= 0) then
               call mpas_log_write(subname // ' ERROR: Failed to allocate var_name_list', &
                                   messageType=MPAS_LOG_CRIT)
            end if

            var_name_list(:) = field_0d_real % constituentnames(:)
         end if

         nullify(field_0d_real)
      case (1)

         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(var_info % name)), field_1d_real, timelevel=1)

         if (.not. associated(field_1d_real)) then
            call mpas_log_write(subname // ' ERROR: Failed to find variable "' // &
                                trim(adjustl(var_info % name)) // '"',            &
                                messageType=MPAS_LOG_CRIT)
         end if

         if (field_1d_real % isvararray .and. associated(field_1d_real % constituentnames)) then
            allocate(var_name_list(size(field_1d_real % constituentnames)), stat=ierr)

            if (ierr /= 0) then
               call mpas_log_write(subname // ' ERROR: Failed to allocate var_name_list', &
                                   messageType=MPAS_LOG_CRIT)
            end if

            var_name_list(:) = field_1d_real % constituentnames(:)
         end if

         nullify(field_1d_real)
      case (2)
         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(var_info % name)), field_2d_real, timelevel=1)

         if (.not. associated(field_2d_real)) then
            call mpas_log_write(subname // ' ERROR: Failed to find variable "' // &
                                trim(adjustl(var_info % name)) // '"',            &
                                messageType=MPAS_LOG_CRIT)
         end if

         if (field_2d_real % isvararray .and. associated(field_2d_real % constituentnames)) then
            allocate(var_name_list(size(field_2d_real % constituentnames)), stat=ierr)

            if (ierr /= 0) then
               call mpas_log_write(subname // ' ERROR: Failed to allocate var_name_list', &
                                   messageType=MPAS_LOG_CRIT)
            end if

            var_name_list(:) = field_2d_real % constituentnames(:)
         end if

         nullify(field_2d_real)
      case (3)
         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(var_info % name)), field_3d_real, timelevel=1)
         if (.not. associated(field_3d_real)) then
            call mpas_log_write(subname // ' ERROR: Failed to find variable "' // &
                                trim(adjustl(var_info % name)) // '"',            &
                                messageType=MPAS_LOG_CRIT)
         end if
         if (field_3d_real % isvararray .and. associated(field_3d_real % constituentnames)) then
            allocate(var_name_list(size(field_3d_real % constituentnames)), stat=ierr)

            if (ierr /= 0) then
               call mpas_log_write(subname // ' ERROR: Failed to allocate var_name_list', &
                                   messageType=MPAS_LOG_CRIT)
            end if

            var_name_list(:) = field_3d_real % constituentnames(:)
         end if
         nullify(field_3d_real)
      case (4)
         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(var_info % name)), field_4d_real, timelevel=1)

         if (.not. associated(field_4d_real)) then
            call mpas_log_write(subname // ' ERROR: Failed to find variable "' // &
                                trim(adjustl(var_info % name)) // '"',            &
                                messageType=MPAS_LOG_CRIT)
         end if

         if (field_4d_real % isvararray .and. associated(field_4d_real % constituentnames)) then
            allocate(var_name_list(size(field_4d_real % constituentnames)), stat=ierr)

            if (ierr /= 0) then
               call mpas_log_write(subname // ' ERROR: Failed to allocate var_name_list', &
                                   messageType=MPAS_LOG_CRIT)
            end if

            var_name_list(:) = field_4d_real % constituentnames(:)
         end if

         nullify(field_4d_real)
      case (5)
         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(var_info % name)), field_5d_real, timelevel=1)

         if (.not. associated(field_5d_real)) then
            call mpas_log_write(subname // ' ERROR: Failed to find variable "' // &
                                trim(adjustl(var_info % name)) // '"',            &
                                messageType=MPAS_LOG_CRIT)
         end if

         if (field_5d_real % isvararray .and. associated(field_5d_real % constituentnames)) then
            allocate(var_name_list(size(field_5d_real % constituentnames)), stat=ierr)

            if (ierr /= 0) then
               call mpas_log_write(subname // ' ERROR: Failed to allocate var_name_list', &
                                   messageType=MPAS_LOG_CRIT)
            end if

            var_name_list(:) = field_5d_real % constituentnames(:)
         end if

         nullify(field_5d_real)
      case default
         call mpas_log_write(subname // ' ERROR: Unsupported variable rank ' // &
                             stringify([var_info % rank]) //                    &
                             ' for "' // trim(adjustl(var_info % name)) // '"', &
                             messageType=MPAS_LOG_CRIT)
      end select
   case default
      call mpas_log_write(subname // ' ERROR: Unsupported variable type ' // &
                          stringify([var_info % type]) //                    &
                          ' for "' // trim(adjustl(var_info % name)) // '"', &
                          messageType=MPAS_LOG_CRIT)
   end select

   if (.not. allocated(var_name_list)) then
      allocate(var_name_list(1), stat=ierr)

      if (ierr /= 0) then
         call mpas_log_write(subname // ' ERROR: Failed to allocate var_name_list', &
                             messageType=MPAS_LOG_CRIT)
      end if

      var_name_list(1) = var_info % name
   end if

   allocate(var_is_present(size(var_name_list)), stat=ierr)

   if (ierr /= 0) then
      call mpas_log_write(subname // ' ERROR: Failed to allocate var_is_present', &
                          messageType=MPAS_LOG_CRIT)
   end if

   var_is_present(:) = .false.
   allocate(var_is_tkr_compatible(size(var_name_list)), stat=ierr)
   if (ierr /= 0) then
      call mpas_log_write(subname // ' ERROR: Failed to allocate var_is_tkr_compatible', &
                          messageType=MPAS_LOG_CRIT)
   end if

   var_is_tkr_compatible(:) = .false.
   if (.not. associated(pio_file)) then
      return
   end if

   if (.not. pio_file_is_open(pio_file)) then
      return
   end if

   if (debug) call mpas_log_write('Checking variable "' // trim(adjustl(var_info % name)) // &
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

   if (debug) then
      call mpas_log_write('var_name_list = ' // stringify(var_name_list))
      call mpas_log_write('var_is_present = ' // stringify(var_is_present))
      call mpas_log_write('var_is_tkr_compatible = ' // stringify(var_is_tkr_compatible))
      call mpas_log_write(subname // ' completed')
   end if
   
 end subroutine dyn_mpas_check_variable_status

 subroutine create_file_timeStamp(curr_time,timeStampOutFile)
   use esmf
   implicit none
   type (MPAS_Time_type),  intent(in ) :: curr_time
   character(len=StrKIND), intent(out) :: timeStampOutFile
   integer :: YY,MM,DD,H,M,S,S_n,S_d,ierr
   character(len=4) :: yy_str
   character(len=2) :: mm_str, dd_str, h_str, m_str, s_str

   call ESMF_TimeGet(curr_time % t, YY=YY, MM=MM, DD=DD, H=H, M=M, S=S, Sn=S_n, Sd=S_d, rc=ierr)
   write(yy_str, '(I4)') YY
   write(mm_str, '(I2.2)') MM
   write(dd_str, '(I2.2)') DD
   write(h_str,  '(I2.2)') H
   write(m_str,  '(I2.2)') M
   write(s_str,  '(I2.2)') S
   timeStampOutFile = yy_str//'-'//mm_str//'-'//dd_str//'_'//h_str//'.'//m_str//'.'//s_str

 end subroutine create_file_timeStamp
end module ufs_mpas_io
