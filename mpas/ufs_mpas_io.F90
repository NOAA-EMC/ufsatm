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
  use mpas_kind_types,     only : StrKIND,RKIND
  use mpas_log,            only : mpas_log_write
  use mpas_derived_types,  only : MPAS_LOG_CRIT, MPAS_LOG_WARN
  use module_mpas_config,  only : pio_iotype, pio_stride, pio_numiotasks, pio_iodesc
  use module_mpas_config,  only : lbc_filename,        pioid_lbc,      pio_subsystem_lbc
  use module_mpas_config,  only : ic_filename,         pioid_ic,       pio_subsystem_ic
  use module_mpas_config,  only : oro_filename,        pioid_oro,      pio_subsystem_oro
  use module_mpas_config,  only :                      pioid_restart,  pio_subsystem_restart
  use module_mpas_config,  only :                      pioid_output,   pio_subsystem_output
  use module_mpas_config,  only : stream_list_history, stream_list_history_funit
  use module_mpas_config,  only : stream_list_diag,    stream_list_diag_funit
  use module_mpas_config,  only : stream_list_restart, stream_list_restart_funit
  use module_mpas_config,  only : TIMELEVEL_NOW
  use ufs_mpas_tools,      only : stringify
  use mpi_f08
  implicit none

  !
  type(core_type),       pointer :: corelist   => null()
  type(domain_type),     pointer :: domain_ptr => null()
  type(mpas_Clock_type), pointer :: clock      => null()

  !> #########################################################################################
  !>
  !> #########################################################################################
  type :: var_info_type
     character(64) :: name = ''
     character(10) :: type = ''
     integer :: rank = 0
  end type var_info_type

  !> #########################################################################################
  !> These variable lists are set at runtime via the stream_list.atmosphere.STREAM files.
  !>
  !> #########################################################################################
  integer, allocatable :: stream_list_history_indices(:)
  integer, allocatable :: stream_list_restart_indices(:)
  integer, allocatable :: stream_list_diag_indices(:)

  !> #########################################################################################
  !> This list corresponds to the "lbc_in" stream in core_atmosphere/Registry.xml
  !> It consists of variables that are members of the "lbc" structure.
  !> #########################################################################################
  type(var_info_type), parameter :: lbc_in_var_info_list(*) = [ &
       var_info_type('xtime'                           , 'character' , 0), &
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
       var_info_type('cellsOnCell'                     , 'integer'   , 2), &
       var_info_type('cellsOnEdge'                     , 'integer'   , 2), &
       var_info_type('cellsOnVertex'                   , 'integer'   , 2), &
       var_info_type('cf1'                             , 'real'      , 0), &
       var_info_type('cf2'                             , 'real'      , 0), &
       var_info_type('cf3'                             , 'real'      , 0), &
       var_info_type('coeffs_reconstruct'              , 'real'      , 3), &
       var_info_type('dcEdge'                          , 'real'      , 1), &
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
       var_info_type('soiltemp'                        , 'real'      , 1), &
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
       var_info_type('relhum'                          , 'real'      , 2), &
       var_info_type('rho'                             , 'real'      , 2), &
       var_info_type('rho_base'                        , 'real'      , 2), &
       var_info_type('scalars'                         , 'real'      , 3), &
       var_info_type('surface_pressure'                , 'real'      , 1), &
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
       var_info_type('zs'                              , 'real'      , 2), &
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
       var_info_type('soiltemp'                        , 'real'      , 1), &
       var_info_type('greenfrac'                       , 'real'      , 2), &
       var_info_type('albedo12m'                       , 'real'      , 2), &
       var_info_type('landusef'                        , 'real'      , 2), &
       var_info_type('soilf'                           , 'real'      , 2), &
       var_info_type('lai12m'                          , 'real'      , 2), &
       var_info_type('canwat'                          , 'real'      , 1), &
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
       var_info_type('scalars'                         , 'real'      , 3), &
       var_info_type('initial_time'                    , 'character' , 0), &
       var_info_type('Time'                            , 'real'      , 0), &
       var_info_type('u'                               , 'real'      , 2), &
       var_info_type('w'                               , 'real'      , 2), &
       var_info_type('rho_zz'                          , 'real'      , 2), &
       var_info_type('theta_m'                         , 'real'      , 2), &
       var_info_type('pressure_p'                      , 'real'      , 2), &
       var_info_type('rho'                             , 'real'      , 2), &
       var_info_type('theta'                           , 'real'      , 2), &
       var_info_type('relhum'                          , 'real'      , 2), &
       var_info_type('circulation'                     , 'real'      , 2), &
       var_info_type('exner'                           , 'real'      , 2), &
       var_info_type('exner_base'                      , 'real'      , 2), &
       var_info_type('rtheta_base'                     , 'real'      , 2), &
       var_info_type('pressure_base'                   , 'real'      , 2), &
       var_info_type('rtheta_p'                        , 'real'      , 2), &
       var_info_type('ru'                              , 'real'      , 2), &
       var_info_type('ru_p'                            , 'real'      , 2), &
       var_info_type('rw'                              , 'real'      , 2), &
       var_info_type('rw_p'                            , 'real'      , 2), &
       var_info_type('rho_p'                           , 'real'      , 2), &
       var_info_type('surface_pressure'                , 'real'      , 1)  &
    ]

  !> #########################################################################################
  !> This list corresponds to the "output" stream in MPAS registry.
  !> It consists of variables that are members of the "diag" structure.
  !> Only variables that are specific to the "output" stream are included.
  !> #########################################################################################
  type(var_info_type), parameter :: history_var_info_list(*) = [ &
  var_info_type('stoch_pattern_sppt'              , 'real',       2), &
       var_info_type('Time'                            , 'real'      , 0), &
       var_info_type('initial_time'                    , 'character' , 0), &
       var_info_type('divergence'                      , 'real'      , 2), &
       var_info_type('pressure'                        , 'real'      , 2), &
       var_info_type('relhum'                          , 'real'      , 2), &
       var_info_type('rho'                             , 'real'      , 2), &
       var_info_type('scalars'                         , 'real'      , 3), &
       var_info_type('surface_pressure'                , 'real'      , 1), &
       var_info_type('theta'                           , 'real'      , 2), &
       var_info_type('u'                               , 'real'      , 2), &
       var_info_type('uReconstructMeridional'          , 'real'      , 2), &
       var_info_type('uReconstructZonal'               , 'real'      , 2), &
       var_info_type('vorticity'                       , 'real'      , 2), &
       var_info_type('w'                               , 'real'      , 2), &
       var_info_type('zz'                              , 'real'      , 2)  &
    ]

  !> #########################################################################################
  !> This list contains physics diagnostics.
  !> It consists of variables that are members of the MPAS "diag_phys" structure.
  !> For the UWM, we include these diagnostics in the MPAS "output" stream.
  !> #########################################################################################
  type(var_info_type), parameter :: diag_phys_var_info_list(*) = [ &
       var_info_type('Time'                            , 'real'      , 0), &
       var_info_type('initial_time'                    , 'character' , 0), &
       var_info_type('swdnb'                           , 'real'      , 1), & !RAD
       var_info_type('swdnbc'                          , 'real'      , 1), & !RAD
       var_info_type('swdnt'                           , 'real'      , 1), & !RAD
       var_info_type('swdntc'                          , 'real'      , 1), & !RAD
       var_info_type('swupb'                           , 'real'      , 1), & !RAD
       var_info_type('swupbc'                          , 'real'      , 1), & !RAD
       var_info_type('swupt'                           , 'real'      , 1), & !RAD
       var_info_type('swuptc'                          , 'real'      , 1), & !RAD
       var_info_type('lwdnb'                           , 'real'      , 1), & !RAD
       var_info_type('lwdnbc'                          , 'real'      , 1), & !RAD
       var_info_type('lwdnt'                           , 'real'      , 1), & !RAD
       var_info_type('lwdntc'                          , 'real'      , 1), & !RAD
       var_info_type('lwupb'                           , 'real'      , 1), & !RAD
       var_info_type('lwupbc'                          , 'real'      , 1), & !RAD
       var_info_type('lwupt'                           , 'real'      , 1), & !RAD
       var_info_type('lwuptc'                          , 'real'      , 1), & !RAD
       !
       var_info_type('refl10cm'                        , 'real'      , 2), & !MP
       var_info_type('rainncv'                         , 'real'      , 1), & !MP
       var_info_type('raincv'                          , 'real'      , 1), & !MP
       var_info_type('snowncv'                         , 'real'      , 1), & !MP
       var_info_type('graupelncv'                      , 'real'      , 1), & !MP
       var_info_type('rainnc'                          , 'real'      , 1), & !MP
       var_info_type('rainc'                           , 'real'      , 1), & !MP
       var_info_type('frainnc'                         , 'real'      , 1), & !MP
       var_info_type('sr'                              , 'real'      , 1), & !MP
       var_info_type('snownc'                          , 'real'      , 1), & !MP
       var_info_type('graupelnc'                       , 'real'      , 1), & !MP
       var_info_type('re_cloud'                        , 'real'      , 2), & !MP
       var_info_type('re_ice'                          , 'real'      , 2), & !MP
       var_info_type('re_snow'                         , 'real'      , 2), & !MP
       !
       var_info_type('sfc_emibck'                      , 'real'      , 1), & !LSM
       var_info_type('mavail'                          , 'real'      , 1), & !LSM
       var_info_type('sfc_albedo'                      , 'real'      , 1), & !LSM/RAD
       var_info_type('sfc_emiss'                       , 'real'      , 1), & !LSM/RAD
       var_info_type('thc'                             , 'real'      , 1), & !LSM
       var_info_type('xicem'                           , 'real'      , 1), & !LSM
       var_info_type('z0'                              , 'real'      , 1), & !LSM
       var_info_type('znt'                             , 'real'      , 1), & !LSM
       !
       var_info_type('ust'                             , 'real'      , 1), & !SFC
       var_info_type('ustm'                            , 'real'      , 1), & !SFC
       var_info_type('qsfc'                            , 'real'      , 1), & !SFC
       var_info_type('hfx'                             , 'real'      , 1), & !SFC
       var_info_type('qfx'                             , 'real'      , 1), & !SFC
       var_info_type('fm'                              , 'real'      , 1), & !SFC
       var_info_type('fh'     	      	      	       , 'real'      , 1), & !SFC
       var_info_type('chs'                             , 'real'      , 1), & !SFC
       var_info_type('cqs'                             , 'real'      , 1), & !SFC
       var_info_type('chs2'                            , 'real'      , 1), & !SFC
       var_info_type('cqs2'                            , 'real'      , 1), & !SFC
       var_info_type('zol'                             , 'real'      , 1), & !SFC
       var_info_type('mol'                             , 'real'      , 1), & !SFC
       var_info_type('rmol'                            , 'real'      , 1), & !SFC
       var_info_type('lh'                              , 'real'      , 1), & !SFC
       var_info_type('t2m'                             , 'real'      , 1), & !SFC
       var_info_type('q2'                              , 'real'      , 1), & !SFC
       var_info_type('skintemp'                        , 'real'      , 1), & !SFC
       var_info_type('sst'                             , 'real'      , 1), & !SFC
       var_info_type('precipw'                         , 'real'      , 1), & !DIAG
       !
       var_info_type('hpbl'                            , 'real'      , 1), & !PBL
       var_info_type('kpbl'                            , 'integer'   , 1), & !PBL
       var_info_type('kzm'                             , 'real'      , 2), & !PBL
       var_info_type('kzh'                             , 'real'      , 2), & !PBL
       var_info_type('kzq'                             , 'real'      , 2), & !PBL
       var_info_type('cldfrac_bl'                      , 'real'      , 2), & !PBL (MYNN)
       var_info_type('qc_bl'                           , 'real'      , 2), & !PBL (MYNN)
       var_info_type('qi_bl'                           , 'real'      , 2), & !PBL (MYNN)
       var_info_type('el_pbl'                          , 'real'      , 2), & !PBL (MYNN)
       var_info_type('sh3d'                            , 'real'      , 2), & !PBL (MYNN)
       var_info_type('sm3d'                            , 'real'      , 2), & !PBL (MYNN)
       var_info_type('qke'                             , 'real'      , 2), & !PBL (MYNN)
       var_info_type('tsq'                             , 'real'      , 2), & !PBL (MYNN)
       var_info_type('qsq'                             , 'real'      , 2), & !PBL (MYNN)
       var_info_type('cov'                             , 'real'      , 2), & !PBL (MYNN)
       !
       var_info_type('dusfcg'                          , 'real'      , 1), & !GWD
       var_info_type('dvsfcg'                          , 'real'      , 1), & !GWD
       var_info_type('dusfc_ls'                        , 'real'      , 1), & !GWD
       var_info_type('dvsfc_ls'                        , 'real'      , 1), & !GWD
       var_info_type('dusfc_bl'                        , 'real'      , 1), & !GWD
       var_info_type('dvsfc_bl'                        , 'real'      , 1), & !GWD
       var_info_type('dusfc_ss'                        , 'real'      , 1), & !GWD
       var_info_type('dvsfc_ss'                        , 'real'      , 1), & !GWD
       var_info_type('dusfc_fd'                        , 'real'      , 1), & !GWD
       var_info_type('dvsfc_fd'                        , 'real'      , 1), & !GWD
       var_info_type('dtaux3d'                         , 'real'      , 2), & !GWD
       var_info_type('dtauy3d'                         , 'real'      , 2), & !GWD
       var_info_type('dtaux3d_ls'                      , 'real'      , 2), & !GWD
       var_info_type('dtauy3d_ls'                      , 'real'      , 2), & !GWD
       var_info_type('dtauy3d_bl'                      , 'real'      , 2), & !GWD
       var_info_type('dtaux3d_bl'                      , 'real'      , 2), & !GWD
       var_info_type('dtaux3d_ss'                      , 'real'      , 2), & !GWD
       var_info_type('dtauy3d_ss'                      , 'real'      , 2), & !GWD
       var_info_type('dtaux3d_fd'                      , 'real'      , 2), & !GWD
       var_info_type('dtauy3d_fd'                      , 'real'      , 2)  & !GWD
       ]

  !> #########################################################################################
  !> Data from MPAS LANDUSE.TBL (read in during init)
  !> #########################################################################################
  integer :: lucats,lumatch,luseas,io,ic,is
  character(len=StrKIND) :: lutype
  real(kind=RKIND) :: li
  real(kind=RKIND),dimension(:,:),allocatable :: albd,slmo,sfem,sfz0,therin,scfx,sfhc

  !> #########################################################################################
  !> Data from MPAS GENPARM.TBL (read in during init)
  !> #########################################################################################
  integer :: num_slope
  real(kind=RKIND),dimension(:),allocatable :: slope_data
  real(kind=RKIND) :: sbeta_data, fxexp_data, csoil_data, salp_data, refdk_data, refkdt_data
  real(kind=RKIND) :: frzk_data, zbot_data, czil_data, smlow_data, smhigh_data

contains

  !> DJS: ToDo. Combine the PIO open routines below.
  
  !> #########################################################################################
  !> Procedure to open MPAS IC file.
  !>
  !> #########################################################################################
  subroutine ufs_mpas_open_init(ierr)
    use pio, only : pio_openfile, pio_nowrite
    integer, intent(out) :: ierr
    logical :: file_exists

    ! Open MPAS Initial Condition file.
    ierr = 0
    INQUIRE(FILE=ic_filename, EXIST=file_exists)
    if (file_exists) then
       ierr = pio_openfile(pio_subsystem_ic, pioid_ic, pio_iotype, ic_filename, pio_nowrite)
    else
       ierr = -1
    end if
    
  end subroutine ufs_mpas_open_init

  !> #########################################################################################
  !> Procedure to open MPAS Lateral Boundary Condition file.
  !>
  !> #########################################################################################
  subroutine ufs_mpas_open_lbc(ierr)
    use pio, only : pio_openfile, pio_nowrite, pio_seterrorhandling, PIO_BCAST_ERROR
    integer, intent(out) :: ierr
    logical :: file_exists

    ! Open MPAS Initial Condition file.
    ierr = 0
    INQUIRE(FILE=lbc_filename, EXIST=file_exists)
    if (file_exists) then
       ierr = pio_openfile(pio_subsystem_lbc, pioid_lbc, pio_iotype, lbc_filename, pio_nowrite)
       if (ierr == 0) call pio_seterrorhandling(pioid_lbc, PIO_BCAST_ERROR)
    else
       ierr = -1
    end if

  end subroutine ufs_mpas_open_lbc

  !> #########################################################################################
  !> Procedure to open MPAS GWD data file.
  !>
  !> #########################################################################################
  subroutine ufs_mpas_open_oro(ierr)
    use pio, only : pio_openfile, pio_nowrite
    integer, intent(out) :: ierr
    logical :: file_exists

    ! Open MPAS orography GWD file.
    ierr = 0
    INQUIRE(FILE=oro_filename, EXIST=file_exists)
    if (file_exists) then
       ierr = pio_openfile(pio_subsystem_oro, pioid_oro, pio_iotype, oro_filename, pio_nowrite)
    else
       ierr = -1
    end if

  end subroutine ufs_mpas_open_oro

  !> #########################################################################################
  !>  Procedure to read in file containing surface properties, LANDUSE.TBL.
  !>
  !> #########################################################################################
  subroutine ufs_mpas_landuse_read(mpicomm, me, master)
    use module_mpas_config,  only : mpas_landuse_file, mpas_land_funit
    use mpas_derived_types,  only : mpas_pool_type
    use mpas_pool_routines,  only : mpas_pool_get_subpool, mpas_pool_get_array
    use mpas_log,            only : mpas_log_write
    use mpas_derived_types,  only : MPAS_LOG_ERR, MPAS_LOG_WARN, MPAS_LOG_CRIT
    type(MPI_Comm), intent(in) :: mpicomm
    integer,        intent(in) :: me, master
    integer :: io,ic,is
    character(len=StrKIND),pointer :: mminlu
    type(mpas_pool_type), pointer :: sfc_input
    character(len=*), parameter :: subname = 'atmos_coupling::ufs_mpas_landuse_read'

    call mpas_log_write(subname //'   Reading/broadcasting MPAS landuse data ')

    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'sfc_input',    sfc_input)
    call mpas_pool_get_array(sfc_input,'mminlu'    ,mminlu  )

    ! Read in data dimensions (master processor)
    if (me == master) then
       open(newunit = mpas_land_funit,file=trim(mpas_landuse_file),action='READ',status='OLD',iostat=io)
       if (io /= 0) then
          call mpas_log_write(subname // " failed to open LANDUSE.TBL", messageType=MPAS_LOG_CRIT)
       endif

       lumatch=0
       find_lutype : do while (lumatch == 0)
          read(unit=mpas_land_funit,fmt='(a35)') lutype
          read(unit=mpas_land_funit,fmt=*) lucats,luseas

          if(lutype .eq. mminlu)then
             call mpas_log_write(subname //'   landuse type = ' // trim (lutype) // ' found '// stringify([lucats]) &
                  // ' categories '// stringify([luseas])//' seasons')
             lumatch=1
          else
             call mpas_log_write(subname //'   skipping over lutype = ' // trim (lutype))
             do is = 1,luseas
                read(unit=mpas_land_funit,fmt=*,iostat=io)
                do ic = 1,lucats
                   read(unit=mpas_land_funit,fmt=*,iostat=io)
                enddo
             enddo
          endif
       enddo find_lutype
    endif

    ! Broadcast data dimensions (all processors)
    call mpi_barrier(mpicomm, io)
    if (io /= 0) call mpas_log_write(subname // " failure waiting for other MPI processes", messageType=MPAS_LOG_CRIT)
    call mpi_bcast(lutype,  strKIND,      MPI_CHARACTER, master, mpicomm, io)
    if (io /= 0) call mpas_log_write(subname // " failure during MPI Broadcast of lutype",  messageType=MPAS_LOG_CRIT)
    call mpi_bcast(luseas,  1,            MPI_INTEGER,   master, mpicomm, io)
    if (io /= 0) call mpas_log_write(subname // " failure during MPI Broadcast of luseas",  messageType=MPAS_LOG_CRIT)
    call mpi_bcast(lucats,  1,            MPI_INTEGER,   master, mpicomm, io)
    if (io /= 0) call mpas_log_write(subname // " failure during MPI Broadcast of lucats",  messageType=MPAS_LOG_CRIT)

    ! Allocate space (all processors)
    allocate(albd(lucats,luseas))
    allocate(slmo(lucats,luseas))
    allocate(sfem(lucats,luseas))
    allocate(sfz0(lucats,luseas))
    allocate(therin(lucats,luseas))
    allocate(scfx(lucats,luseas))
    allocate(sfhc(lucats,luseas))

    ! Read data (master processor)
    if (me == master) then
       do is = 1, luseas
          read(unit=mpas_land_funit,fmt=*,iostat=io)
          if (io /= 0) then
             call mpas_log_write(subname // " failure reading LANDUSE.TBL", messageType=MPAS_LOG_CRIT)
          end if
          do ic = 1, lucats
             read(unit=mpas_land_funit,fmt=*) li,albd(ic,is),slmo(ic,is),sfem(ic,is),sfz0(ic,is), &
                  therin(ic,is),scfx(ic,is),sfhc(ic,is)
          end do
       end do
       close(mpas_land_funit)
    end if

    ! Broadcast data (all processors)
    call mpi_barrier(mpicomm, io)
    call mpi_bcast(albd,    size(albd),   mpi_real,     master, mpicomm, io)
    if (io /= 0) call mpas_log_write(subname // " failure during MPI Broadcast of albd",    messageType=MPAS_LOG_CRIT)
    call mpi_bcast(slmo,    size(slmo),   mpi_real,     master, mpicomm, io)
    if (io /= 0) call mpas_log_write(subname // " failure during MPI Broadcast of slmo",    messageType=MPAS_LOG_CRIT)
    call mpi_bcast(sfem,    size(sfem),   mpi_real,     master, mpicomm, io)
    if (io /= 0) call mpas_log_write(subname // " failure during MPI Broadcast of sfem",    messageType=MPAS_LOG_CRIT)
    call mpi_bcast(sfz0,    size(sfz0),   mpi_real,     master, mpicomm, io)
    if (io /= 0) call mpas_log_write(subname // " failure during MPI Broadcast of sfz0",    messageType=MPAS_LOG_CRIT)
    call mpi_bcast(therin,  size(therin), mpi_real,     master, mpicomm, io)
    if (io /= 0) call mpas_log_write(subname // " failure during MPI Broadcast of therin",  messageType=MPAS_LOG_CRIT)
    call mpi_bcast(sfhc,    size(sfhc),   mpi_real,     master, mpicomm, io)
    if (io /= 0) call mpas_log_write(subname // " failure during MPI Broadcast of sfhc",    messageType=MPAS_LOG_CRIT)
    call mpi_bcast(scfx,    size(scfx),   mpi_real,     master, mpicomm, io)
    if (io /= 0) call mpas_log_write(subname // " failure during MPI Broadcast of scfx",    messageType=MPAS_LOG_CRIT)

    call mpas_log_write(subname //'   Finished reading/broadcasting MPAS landuse data ')

  end subroutine ufs_mpas_landuse_read

  !> #########################################################################################
  !> Procedure to read in slope data needed by RUC LSM.
  !>
  !> #########################################################################################
  subroutine use_mpas_slopedata_read(mpicomm, me, master)
    use module_mpas_config,  only : mpas_genparm_file, mpas_gprm_funit
    use mpas_log,            only : mpas_log_write
    use mpas_derived_types,  only : MPAS_LOG_ERR, MPAS_LOG_WARN, MPAS_LOG_CRIT
    !
    type(MPI_Comm), intent(in) :: mpicomm
    integer,        intent(in) :: me, master
    ! Locals
    integer :: io, lc, slpcats
    character(len=*), parameter :: subname = 'atmos_coupling::ufs_mpas_genparm_read'

    call mpas_log_write(subname //'   Reading/broadcasting MPAS RUC LSM slope data ')

    ! Read in data dimensions (master processor)
    if (me == master) then
       open(newunit=mpas_gprm_funit,file=trim(mpas_genparm_file),form='FORMATTED',status='OLD',iostat=io)
       if (io /= 0) then
          call mpas_log_write(subname // " failed to open GENPARM.TBL", messageType=MPAS_LOG_CRIT)
       endif
       read(mpas_gprm_funit,*)
       read(mpas_gprm_funit,*)
       read(mpas_gprm_funit,*) num_slope

       slpcats=num_slope
    end if

    ! Broadcast data dimensions (all processors)
    call mpi_barrier(mpicomm, io)
    if (io /= 0) call mpas_log_write(subname // " failure waiting for other MPI processes",   messageType=MPAS_LOG_CRIT)
    call mpi_bcast(num_slope,  1, mpi_real, master, mpicomm, io)
    if (io /= 0) call mpas_log_write(subname // " failure during MPI Broadcast of num_slope", messageType=MPAS_LOG_CRIT)

    ! Allocate space (all processors)
    allocate(slope_data(num_slope))

    ! Read data (master processor)
    if (me == master) then
       do lc = 1, slpcats
          read(mpas_gprm_funit,*) slope_data(lc)
       enddo
       read(mpas_gprm_funit,*)
       read(mpas_gprm_funit,*) sbeta_data
       read(mpas_gprm_funit,*)
       read(mpas_gprm_funit,*) fxexp_data
       read(mpas_gprm_funit,*)
       read(mpas_gprm_funit,*) csoil_data
       read(mpas_gprm_funit,*)
       read(mpas_gprm_funit,*) salp_data
       read(mpas_gprm_funit,*)
       read(mpas_gprm_funit,*) refdk_data
       read(mpas_gprm_funit,*)
       read(mpas_gprm_funit,*) refkdt_data
       read(mpas_gprm_funit,*)
       read(mpas_gprm_funit,*) frzk_data
       read(mpas_gprm_funit,*)
       read(mpas_gprm_funit,*) zbot_data
       read(mpas_gprm_funit,*)
       read(mpas_gprm_funit,*) czil_data
       read(mpas_gprm_funit,*)
       read(mpas_gprm_funit,*) smlow_data
       read(mpas_gprm_funit,*)
       read(mpas_gprm_funit,*) smhigh_data
       close(mpas_gprm_funit)
    end if

    ! Broadcast data (all processors)
    call mpi_barrier(mpicomm, io)
    if (io /= 0) call mpas_log_write(subname // " failure waiting for other MPI processes",      messageType=MPAS_LOG_CRIT)
    call mpi_bcast(slope_data,   size(slope_data),  mpi_real,     master, mpicomm, io)
    if (io /= 0) call mpas_log_write(subname // " failure during MPI Broadcast of slope_data",   messageType=MPAS_LOG_CRIT)
    call mpi_bcast(sbeta_data,   1,  mpi_real,     master, mpicomm, io)
    if (io /= 0) call mpas_log_write(subname // " failure during MPI Broadcast of sbeta_data",   messageType=MPAS_LOG_CRIT)
    call mpi_bcast(fxexp_data,   1,  mpi_real,     master, mpicomm, io)
    if (io /= 0) call mpas_log_write(subname // " failure during MPI Broadcast of fxexp_data",   messageType=MPAS_LOG_CRIT)
    call mpi_bcast(csoil_data,   1,  mpi_real,     master, mpicomm, io)
    if (io /= 0) call mpas_log_write(subname // " failure during MPI Broadcast of csoil_data",   messageType=MPAS_LOG_CRIT)
    call mpi_bcast(salp_data,    1,  mpi_real,     master, mpicomm, io)
    if (io /= 0) call mpas_log_write(subname // " failure during MPI Broadcast of salp_data",    messageType=MPAS_LOG_CRIT)
    call mpi_bcast(refdk_data,   1,  mpi_real,     master, mpicomm, io)
    if (io /= 0) call mpas_log_write(subname // " failure during MPI Broadcast of refdk_data",   messageType=MPAS_LOG_CRIT)
    call mpi_bcast(refkdt_data,  1,  mpi_real,     master, mpicomm, io)
    if (io /= 0) call mpas_log_write(subname // " failure during MPI Broadcast of refkdt_data",  messageType=MPAS_LOG_CRIT)
    call mpi_bcast(frzk_data,    1,  mpi_real,     master, mpicomm, io)
    if (io /= 0) call mpas_log_write(subname // " failure during MPI Broadcast of frzk_data",    messageType=MPAS_LOG_CRIT)
    call mpi_bcast(zbot_data,    1,  mpi_real,     master, mpicomm, io)
    if (io /= 0) call mpas_log_write(subname // " failure during MPI Broadcast of zbot_data",    messageType=MPAS_LOG_CRIT)
    call mpi_bcast(czil_data,    1,  mpi_real,     master, mpicomm, io)
    if (io /= 0) call mpas_log_write(subname // " failure during MPI Broadcast of czil_data",    messageType=MPAS_LOG_CRIT)
    call mpi_bcast(smlow_data,   1,  mpi_real,     master, mpicomm, io)
    if (io /= 0) call mpas_log_write(subname // " failure during MPI Broadcast of smlow_data",   messageType=MPAS_LOG_CRIT)
    call mpi_bcast(smhigh_data,  1,  mpi_real,     master, mpicomm, io)
    if (io /= 0) call mpas_log_write(subname // " failure during MPI Broadcast of smhigh_data",  messageType=MPAS_LOG_CRIT)

    call mpas_log_write(subname //'   Finished reading/broadcasting MPAS RUC LSM slope data ')
  end subroutine use_mpas_slopedata_read

  !> #########################################################################################
  !> Procedure to populate MPAS diag_phys pool with CCPP data.
  !>
  !> The fields from diag_phys are allocated by MPAS, populated with CCPP Physics data, and
  !> written to output. So essentially we copy physics arrays back into MPAS memory to use
  !> MPAS's native output functionality.
  !>
  !> CAN BE REMOVED/NOT-CALLED WHEN ESMF WRITE GRID COMPONENT IMPLEMENTED.
  !> #########################################################################################
  subroutine ufs_mpas_phys_diag(statein, control, radiation, diagnostics, tbd, surface)
    use GFS_typedefs,         only : GFS_statein_type
    use GFS_typedefs,         only : GFS_control_type
    use GFS_typedefs,         only : GFS_radtend_type
    use GFS_typedefs,         only : GFS_diag_type
    use GFS_typedefs,         only : GFS_tbd_type
    use GFS_typedefs,         only : GFS_sfcprop_type
    use mpas_derived_types,   only : mpas_pool_type
    use mpas_pool_routines,   only : mpas_pool_get_subpool, mpas_pool_get_dimension
    use mpas_pool_routines,   only : mpas_pool_get_array, mpas_pool_get_config

    ! Arguments
    type(GFS_statein_type), intent(in) :: statein
    type(GFS_control_type), intent(in) :: control
    type(GFS_radtend_type), intent(in) :: radiation
    type(GFS_diag_type),    intent(in) :: diagnostics
    type(GFS_tbd_type),     intent(in) :: tbd
    type(GFS_sfcprop_type), intent(in) :: surface
    ! Locals
    type(mpas_pool_type), pointer :: diag_phys, sfc_input
    real(RKIND), pointer :: t2m(:), q2(:), skintemp(:), sst(:), precipw(:)
    real(RKIND), pointer :: swdnb(:),swdnbc(:),swupb(:),swupbc(:)
    real(RKIND), pointer :: lwdnb(:),lwdnbc(:),lwupb(:),lwupbc(:)
    real(RKIND), pointer :: swdnt(:),swdntc(:),swupt(:),swuptc(:)
    real(RKIND), pointer :: lwdnt(:),lwdntc(:),lwupt(:),lwuptc(:)
    real(RKIND), pointer :: re_cloud(:,:),re_ice(:,:),re_snow(:,:)
    real(RKIND), pointer :: sr(:),hfx(:),qfx(:),ust(:),znt(:),qsfc(:), ustm(:)
    real(RKIND), pointer :: fm(:),fh(:),chs(:),cqs(:),chs2(:),cqs2(:),lh(:)
    real(RKIND), pointer :: mol(:), rmol(:),zol(:),hpbl(:),kzm(:,:),kzh(:,:),kzq(:,:)
    real(RKIND), pointer :: sh3d(:,:),sm3d(:,:),cldfrac_bl(:,:),qc_bl(:,:),qi_bl(:,:),el_pbl(:,:),qke(:,:),tsq(:,:),qsq(:,:),cov(:,:)
    real(RKIND), pointer :: refl10cm(:,:),mavail(:)
    real(RKIND), pointer :: rainc(:),rainnc(:),frainnc(:),snownc(:),graupelnc(:)
    real(RKIND), pointer :: raincv(:),rainncv(:),snowncv(:),graupelncv(:)
    real(RKIND), pointer :: dusfcg(:),dvsfcg(:),dusfc_ls(:),dvsfc_ls(:)
    real(RKIND), pointer :: dusfc_bl(:),dvsfc_bl(:),dusfc_ss(:),dvsfc_ss(:)
    real(RKIND), pointer :: dusfc_fd(:),dvsfc_fd(:)
    real(RKIND), pointer :: dtaux3d(:,:), dtauy3d(:,:)
    real(RKIND), pointer :: dtaux3d_ls(:,:), dtauy3d_ls(:,:), dtaux3d_ss(:,:), dtauy3d_ss(:,:)
    real(RKIND), pointer :: dtaux3d_fd(:,:), dtauy3d_fd(:,:), dtaux3d_bl(:,:), dtauy3d_bl(:,:)
    integer,     pointer :: nThreads, cellSolveThreadStart(:), cellSolveThreadEnd(:), kpbl(:)
    integer :: iCol, ithread
    real(RKIND), parameter :: m2mm = 1000._RKIND, um2m = 1.e-6_RKIND
    real(RKIND), parameter :: cp = 1004.5_RKIND, rd = 287.05_RKIND, rv = 461.5_RKIND, hvap = 2.5e6_RKIND
    real(RKIND) :: rho_sfc
    character(len=*), parameter :: subname = 'ufs_mpas_io::ufs_mpas_phys_diag'

    ! Get openMP information
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'nThreads',             nThreads)
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'cellSolveThreadStart', cellSolveThreadStart)
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'cellSolveThreadEnd',   cellSolveThreadEnd)

    ! Access MPAS data pools.
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'diag_physics',  diag_phys)
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'sfc_input', sfc_input)

    ! Grab fields from MPAS pools
    call mpas_pool_get_array(diag_phys,'t2m'       , t2m       )
    call mpas_pool_get_array(diag_phys,'q2'        , q2        )
    call mpas_pool_get_array(diag_phys,'precipw'   , precipw   )
    call mpas_pool_get_array(sfc_input,'skintemp'  , skintemp  )
    call mpas_pool_get_array(sfc_input,'sst'       , sst       )
    call mpas_pool_get_array(diag_phys,'swdnb'     , swdnb     )
    call mpas_pool_get_array(diag_phys,'swdnbc'    , swdnbc    )
    call mpas_pool_get_array(diag_phys,'swupb'     , swupb     )
    call mpas_pool_get_array(diag_phys,'swupbc'    , swupbc    )
    call mpas_pool_get_array(diag_phys,'lwdnb'     , lwdnb     )
    call mpas_pool_get_array(diag_phys,'lwdnbc'    , lwdnbc    )
    call mpas_pool_get_array(diag_phys,'lwupb'     , lwupb     )
    call mpas_pool_get_array(diag_phys,'lwupbc'    , lwupbc    )
    call mpas_pool_get_array(diag_phys,'swdnt'     , swdnt     )
    call mpas_pool_get_array(diag_phys,'swdntc'    , swdntc    )
    call mpas_pool_get_array(diag_phys,'swupt'     , swupt     )
    call mpas_pool_get_array(diag_phys,'swuptc'    , swuptc    )
    call mpas_pool_get_array(diag_phys,'lwdnt'     , lwdnt     )
    call mpas_pool_get_array(diag_phys,'lwdntc'    , lwdntc    )
    call mpas_pool_get_array(diag_phys,'lwupt'     , lwupt     )
    call mpas_pool_get_array(diag_phys,'lwuptc'    , lwuptc    )
    
    call mpas_pool_get_array(diag_phys,'refl10cm'  , refl10cm  )
    call mpas_pool_get_array(diag_phys,'rainc'     , rainc     )
    call mpas_pool_get_array(diag_phys,'rainnc'    , rainnc    )
    call mpas_pool_get_array(diag_phys,'frainnc'   , frainnc   )
    call mpas_pool_get_array(diag_phys,'snownc'    , snownc    )
    call mpas_pool_get_array(diag_phys,'graupelnc' , graupelnc )
    call mpas_pool_get_array(diag_phys,'sr'        , sr        )
    call mpas_pool_get_array(diag_phys,'raincv'    , raincv    )
    call mpas_pool_get_array(diag_phys,'rainncv'   , rainncv   )
    call mpas_pool_get_array(diag_phys,'snowncv'   , snowncv   )
    call mpas_pool_get_array(diag_phys,'graupelncv', graupelncv)
    call mpas_pool_get_array(diag_phys,'re_cloud'  , re_cloud  )
    call mpas_pool_get_array(diag_phys,'re_ice'    , re_ice    )
    call mpas_pool_get_array(diag_phys,'re_snow'   , re_snow   )
    call mpas_pool_get_array(diag_phys,'hfx'       , hfx       )
    call mpas_pool_get_array(diag_phys,'qfx'       , qfx       )
    call mpas_pool_get_array(diag_phys,'ust'       , ust       )
    call mpas_pool_get_array(diag_phys,'ustm'      , ustm      )
    call mpas_pool_get_array(diag_phys,'znt'       , znt       )
    call mpas_pool_get_array(diag_phys,'qsfc'      , qsfc      )
    call mpas_pool_get_array(diag_phys,'mavail'    , mavail    )
    call mpas_pool_get_array(diag_phys,'fm'        , fm        )
    call mpas_pool_get_array(diag_phys,'fh'        , fh        )
    if (control % do_mynnedmf) then
       call mpas_pool_get_array(diag_phys,'sh3d'      , sh3d      )
       call mpas_pool_get_array(diag_phys,'sm3d'      , sm3d      )
       call mpas_pool_get_array(diag_phys,'cldfrac_bl', cldfrac_bl)
       call mpas_pool_get_array(diag_phys,'qc_bl'     , qc_bl     )
       call mpas_pool_get_array(diag_phys,'qi_bl'     , qi_bl     )
       call mpas_pool_get_array(diag_phys,'el_pbl'    , el_pbl    )
       call mpas_pool_get_array(diag_phys,'qke'       , qke       )
       call mpas_pool_get_array(diag_phys,'tsq'       , tsq       )
       call mpas_pool_get_array(diag_phys,'qsq'       , qsq       )
       call mpas_pool_get_array(diag_phys,'cov'       , cov       )
    endif
    call mpas_pool_get_array(diag_phys,'chs'       , chs       )
    call mpas_pool_get_array(diag_phys,'cqs'       , cqs       )
    call mpas_pool_get_array(diag_phys,'chs2'      , chs2      )
    call mpas_pool_get_array(diag_phys,'cqs2'      , cqs2      )
    call mpas_pool_get_array(diag_phys,'zol'       , zol       )
    call mpas_pool_get_array(diag_phys,'mol'       , mol       )
    call mpas_pool_get_array(diag_phys,'rmol'      , rmol      )
    call mpas_pool_get_array(diag_phys,'lh'        , lh        )
    call mpas_pool_get_array(diag_phys,'hpbl'      , hpbl      )
    call mpas_pool_get_array(diag_phys,'kpbl'      , kpbl      )
    call mpas_pool_get_array(diag_phys,'kzm'       , kzm       )
    call mpas_pool_get_array(diag_phys,'kzh'       , kzh       )
    call mpas_pool_get_array(diag_phys,'kzq'       , kzq       )
    ! UFS GWD diagnostics are conditionally allocated.
    if (control % ldiag_ugwp .or. control % do_ugwp_v1) then
       call mpas_pool_get_array(diag_phys,'dusfcg'    , dusfcg    )
       call mpas_pool_get_array(diag_phys,'dvsfcg'    , dvsfcg    )
       call mpas_pool_get_array(diag_phys,'dusfc_ls'  , dusfc_ls  )
       call mpas_pool_get_array(diag_phys,'dvsfc_ls'  , dvsfc_ls  )
       call mpas_pool_get_array(diag_phys,'dusfc_bl'  , dusfc_bl  )
       call mpas_pool_get_array(diag_phys,'dvsfc_bl'  , dvsfc_bl  )
       call mpas_pool_get_array(diag_phys,'dusfc_ss'  , dusfc_ss  )
       call mpas_pool_get_array(diag_phys,'dvsfc_ss'  , dvsfc_ss  )
       call mpas_pool_get_array(diag_phys,'dusfc_fd'  , dusfc_fd  )
       call mpas_pool_get_array(diag_phys,'dvsfc_fd'  , dvsfc_fd  )
       call mpas_pool_get_array(diag_phys,'dtaux3d'   , dtaux3d   )
       call mpas_pool_get_array(diag_phys,'dtauy3d'   , dtauy3d   )
       call mpas_pool_get_array(diag_phys,'dtaux3d_ls', dtaux3d_ls)
       call mpas_pool_get_array(diag_phys,'dtauy3d_ls', dtauy3d_ls)
       call mpas_pool_get_array(diag_phys,'dtaux3d_ss', dtaux3d_ss)
       call mpas_pool_get_array(diag_phys,'dtauy3d_ss', dtauy3d_ss)
       call mpas_pool_get_array(diag_phys,'dtaux3d_bl', dtaux3d_bl)
       call mpas_pool_get_array(diag_phys,'dtauy3d_bl', dtauy3d_bl)
       call mpas_pool_get_array(diag_phys,'dtaux3d_fd', dtaux3d_fd)
       call mpas_pool_get_array(diag_phys,'dtauy3d_fd', dtauy3d_fd)
    end if

    do ithread = 1,nThreads
       do iCol = cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread)
          ! --- Radiation fluxes at surface
          swdnb(iCol)  = radiation%sfcfsw(iCol)%dnfxc
          swdnbc(iCol) = radiation%sfcfsw(iCol)%dnfx0
          swupb(iCol)  = radiation%sfcfsw(iCol)%upfxc
          swupbc(iCol) = radiation%sfcfsw(iCol)%upfx0
          lwdnb(iCol)  = radiation%sfcflw(iCol)%dnfxc
          lwdnbc(iCol) = radiation%sfcflw(iCol)%dnfx0
          lwupb(iCol)  = radiation%sfcflw(iCol)%upfxc
          lwupbc(iCol) = radiation%sfcflw(iCol)%upfx0

          ! --- TOA fluxes (RRTMG: no clear-sky SW-down at TOA separate from all-sky; no LW-down at TOA)
          swdnt(iCol)  = diagnostics%topfsw(iCol)%dnfxc
          swdntc(iCol) = diagnostics%topfsw(iCol)%dnfxc
          swupt(iCol)  = diagnostics%topfsw(iCol)%upfxc
          swuptc(iCol) = diagnostics%topfsw(iCol)%upfx0
          lwdnt(iCol)  = 0._RKIND
          lwdntc(iCol) = 0._RKIND
          lwupt(iCol)  = diagnostics%topflw(iCol)%upfxc
          lwuptc(iCol) = diagnostics%topflw(iCol)%upfx0

          ! --- Other
          refl10cm(:,iCol) = diagnostics%refl_10cm(iCol,:)
          sr(iCol)         = diagnostics%sr(iCol)

          ! --- Precipitation: UFS m (lwe) -> MPAS mm
          !     (totsnw/totgrp are declared kg m-2 in GFS_typedefs.meta, but GFS_MP_generic_post
          !      accumulates Diag%snow/graupel, which are m -> same factor)
          raincv(iCol)     = m2mm * diagnostics%rainc(iCol)
          rainncv(iCol)    = m2mm * max(0._RKIND, diagnostics%rain(iCol) - diagnostics%rainc(iCol))
          snowncv(iCol)    = m2mm * diagnostics%snow(iCol)
          graupelncv(iCol) = m2mm * diagnostics%graupel(iCol)
          rainc(iCol)      = m2mm * diagnostics%cnvprcp(iCol)
          rainnc(iCol)     = m2mm * max(0._RKIND, diagnostics%totprcp(iCol) - diagnostics%cnvprcp(iCol))
          frainnc(iCol)    = m2mm * diagnostics%frzr(iCol)
          snownc(iCol)     = m2mm * diagnostics%totsnw(iCol)
          graupelnc(iCol)  = m2mm * diagnostics%totgrp(iCol)

          ! --- Effective radii: UFS micron -> MPAS m
          re_cloud(:,iCol) = um2m * tbd%phy_f3d(iCol,:,control%nleffr)
          re_ice(:,iCol)   = um2m * tbd%phy_f3d(iCol,:,control%nieffr)
          re_snow(:,iCol)  = um2m * tbd%phy_f3d(iCol,:,control%nseffr)

          ! Gravity-wave physics diagnostics (conditionally allocated)
          if (control % ldiag_ugwp .or. control % do_ugwp_v1) then
             dusfcg(iCol)       = diagnostics%dusfcg(iCol)
             dvsfcg(iCol)       = diagnostics%dvsfcg(iCol)
             dusfc_ls(iCol)     = diagnostics%du_ogwcol(iCol)
             dvsfc_ls(iCol)     = diagnostics%dv_ogwcol(iCol)
             dusfc_bl(iCol)     = diagnostics%du_oblcol(iCol)
             dvsfc_bl(iCol)     = diagnostics%dv_oblcol(iCol)
             dusfc_ss(iCol)     = diagnostics%du_osscol(iCol)
             dvsfc_ss(iCol)     = diagnostics%dv_osscol(iCol)
             dusfc_fd(iCol)     = diagnostics%du_ofdcol(iCol)
             dvsfc_fd(iCol)     = diagnostics%dv_ofdcol(iCol)
             dtaux3d(:,iCol)    = diagnostics%dudt_gw(iCol,:)
             dtauy3d(:,iCol)    = diagnostics%dvdt_gw(iCol,:)
             dtaux3d_ls(:,iCol) = diagnostics%dudt_ogw(iCol,:)
             dtauy3d_ls(:,iCol) = diagnostics%dvdt_ogw(iCol,:)
             dtaux3d_ss(:,iCol) = diagnostics%dudt_oss(iCol,:)
             dtauy3d_ss(:,iCol) = diagnostics%dvdt_oss(iCol,:)
             dtaux3d_bl(:,iCol) = diagnostics%dudt_obl(iCol,:)
             dtauy3d_bl(:,iCol) = diagnostics%dvdt_obl(iCol,:)
             dtaux3d_fd(:,iCol) = diagnostics%dudt_ofd(iCol,:)
             dtauy3d_fd(:,iCol) = diagnostics%dvdt_ofd(iCol,:)
          end if
          ! --- Surface layer: convert kinematic / rho-weighted UFS fields to MPAS units
          rho_sfc = statein%prsl(iCol,1) / (rd * statein%tgrs(iCol,1) * &
                    (1._RKIND + (rv/rd - 1._RKIND) * statein%qgrs(iCol,1,control%ntqv)))
          hfx(iCol)      = rho_sfc * cp * surface%hflx(iCol)       ! K m s-1        -> W m-2
          qfx(iCol)      = rho_sfc * surface%evap(iCol)            ! kg kg-1 m s-1  -> kg m-2 s-1
          znt(iCol)      = 0.01_RKIND * surface%zorl(iCol)         ! cm             -> m
          chs(iCol)      = surface%flhc(iCol) / (rho_sfc * cp)     ! W m-2 K-1      -> m s-1
          cqs(iCol)      = surface%flqc(iCol) / rho_sfc            ! kg m-2 s-1     -> m s-1
          rmol(iCol)     = surface%rmol(iCol)
          ust(iCol)      = surface%uustar(iCol)
          qsfc(iCol)     = surface%qss(iCol)
          t2m(iCol)      = surface%t2m(iCol)
          q2(iCol)       = surface%q2m(iCol)
          skintemp(iCol) = surface%tsfc(iCol)
          sst(iCol)      = surface%tsfco(iCol)
          precipw(iCol)  = diagnostics%pwat(iCol)
          ! LSM (RUC)
          if (control % lsm == control % lsm_ruc) then
             if (nint(surface%slmsk(iCol)) == 1) mavail(iCol) = surface%wetness(iCol)
          end if
          ! MONIN-OBUKHOV
          fm(iCol)     = surface%ffmm(iCol)
          fh(iCol)     = surface%ffhh(iCol)
          ! MYNN (surface-layer model)
          if (control % do_mynnsfclay) then
             ustm(iCol)   = surface%ustm(iCol)
             chs2(iCol)   = surface%chs2(iCol)
             cqs2(iCol)   = surface%cqs2(iCol)
             zol(iCol)    = surface%zol(iCol)
             mol(iCol)    = surface%mol(iCol)
             lh(iCol)     = 2.5e6_RKIND * qfx(iCol)
          endif
          ! MYNN (PBL)
          if (control % do_mynnedmf) then
             sh3d(:,iCol)       = tbd%sh3d(iCol,:)
             sm3d(:,iCol)       = tbd%sm3d(iCol,:)
             cldfrac_bl(:,iCol) = tbd%cldfra_bl(iCol,:)
             qc_bl(:,iCol)      = tbd%qc_bl(iCol,:)
             qi_bl(:,iCol)      = tbd%qi_bl(iCol,:)
             el_pbl(:,iCol)     = tbd%el_pbl(iCol,:)
             qke(:,iCol)        = tbd%qke(iCol,:)
             tsq(:,iCol)        = tbd%tsq(iCol,:)
             qsq(:,iCol)        = tbd%qsq(iCol,:)
             cov(:,iCol)        = tbd%cov(iCol,:)
          endif
          hpbl(iCol) = tbd%hpbl(iCol)
          kpbl(iCol) = tbd%kpbl(iCol)
          !kzm(:,iCol)  =
          !kzh(:,iCol)  =
          !kzq(:,iCol)  =
       end do
    end do
  end subroutine ufs_mpas_phys_diag

  !> #########################################################################################
  !> Procedure to read in stream_list (a.k.a File with fields to include in output stream)
  !> 
  !> #########################################################################################
  subroutine ufs_mpas_read_stream_lists(me, master, mpicomm)
    ! Arguments
    integer,        intent(in) :: me, master
    type(MPI_Comm), intent(in) :: mpicomm
    ! Locals
    character(len=*), parameter :: subname = 'ufs_mpas_io::ufs_mpas_read_stream_list'
    
    call read_stream_list(me, master, mpicomm, stream_list_history,  stream_list_history_funit, 'history')
    call read_stream_list(me, master, mpicomm, stream_list_diag,     stream_list_diag_funit,    'diag_phys')
    call read_stream_list(me, master, mpicomm, stream_list_restart,  stream_list_restart_funit, 'restart')

  end subroutine ufs_mpas_read_stream_lists

  !> #########################################################################################
  !> The procedure reads in a MPAS stream_list and compares the requested variables in the
  !> <stream_list_file> to the available varaibles in  <stream_name>_var_info_list.
  !>
  !> If no MPAS stream_list is provided, all fields included in <stream_name>_var_info_list
  !> will be included in the output file.
  !>
  !> #########################################################################################
  subroutine read_stream_list(me, master, mpicomm, stream_list_file, funit, stream_name)
    integer, intent(in) :: me, master
    type(MPI_Comm), intent(in) :: mpicomm
    character(len=*), intent(in) :: stream_list_file
    integer, intent(inout) :: funit
    character(len=*), intent(in) :: stream_name
    integer :: nvars, ivar, io, i, nvar_av, count, mpierr
    integer :: nvars_hist, nvars_rstr,nvars_diag
    logical :: file_exists, found
    character(len=128) :: line_buffer
    character(len=128), allocatable :: var_list(:)
    integer, allocatable :: indices_temp(:), indices(:)
    character(len=*), parameter :: subname = 'ufs_mpas_io::read_stream_list'
    type(var_info_type), allocatable :: var_info_list(:)

    ! Check if file exists before trying to read.
    file_exists = .false.
    INQUIRE(FILE=stream_list_file, EXIST=file_exists)
    if (.not. file_exists) then
       call mpas_log_write(subname // " No stream_list file provided for "//trim(stream_name)// &
                           ". All available fields will be output", messageType=MPAS_LOG_WARN)
       return
    end if
    
    ! Set var_info_list for given <stream_name>.
    select case (trim(adjustl(stream_name)))
    case ('')
       allocate(var_info_list(0))
    case ('restart')
       allocate(var_info_list, source=restart_var_info_list)
       nvars_rstr = size(restart_var_info_list)
    case ('history')
       allocate(var_info_list, source=history_var_info_list)
       nvars_hist = size(history_var_info_list)
    case ('diag_phys')
       allocate(var_info_list, source=diag_phys_var_info_list)
       nvars_diag = size(diag_phys_var_info_list)
    end select

    ! On master process...
    if (me == master) then
       ! Get number of lines (variables) in stream_list file.
       open(newunit=funit,file=trim(stream_list_file))
       nvars = 0
       do
          read(funit, *, iostat=io)
          if (io /= 0) exit  ! Exit loop on end-of-file or error
          nvars = nvars + 1
       end do
       close(funit)

       ! Read in stream_list from file.
       allocate(indices_temp(nvars))
       indices_temp(:) = -999
       allocate(var_list(nvars))
       open(newunit=funit,file=trim(stream_list_file))
       do iVar = 1,nvars
          read(funit, '(A)', iostat=io) line_buffer
          var_list(ivar) = line_buffer
       enddo
       close(funit)

       ! Are requested stream_list variables available in <var_info_list>?
       ! Loop over requested variables in <stream_list>, <nvars>, and check
       ! for existence in <var_info_list>.
       nvar_av = 0
       do iVar = 1,nvars
          do i = 1, size(var_info_list)
             found = .false.
             if (trim(var_list(ivar)) == trim(var_info_list(i)%name)) then
                found = .true.
                nvar_av = nvar_av + 1
                indices_temp(nvar_av) = i
             endif
          enddo
          ! If not found, requested variables is not supported. Print warning message.
          if (.not. found) then
             call mpas_log_write(subname // " Variable not supported, "//trim(var_list(ivar))// &
                                 ", skipping", messageType=MPAS_LOG_WARN)
          end if
       end do

       ! Handle case when fields requested in stream_list are not available.
       if (nvar_av .ne. nvars) then
          allocate(indices(nvar_av))
          count = 0
          do iVar = 1,nvars
             if (indices_temp(ivar) .ne. -999) then
                count = count + 1
                indices(count) = indices_temp(ivar)
             end if
          end do
          nvars = count
       ! Otherwise, use full requested variable list.
       else
          allocate(indices(nvars))
          indices = indices_temp
       end if
    end if

    ! Other processors waiting...
    call mpi_barrier(mpicomm, mpierr)

    ! Broadcast dimension
    select case (trim(adjustl(stream_name)))
    case ('restart')
       nvars_rstr = nvars
       call mpi_bcast(nvars_rstr, 1, MPI_INTEGER, master, mpicomm, mpierr)
    case ('history')
       nvars_hist = nvars
       call mpi_bcast(nvars_hist, 1, MPI_INTEGER, master, mpicomm, mpierr)
    case ('diag_phys')
       nvars_diag = nvars
       call mpi_bcast(nvars_diag, 1, MPI_INTEGER, master, mpicomm, mpierr)
    end select
    
    ! Allocate
    select case (trim(adjustl(stream_name)))
    case ('restart')
       allocate(stream_list_restart_indices(nvars_rstr))
    case ('history')
       allocate(stream_list_history_indices(nvars_hist))
    case ('diag_phys')
       allocate(stream_list_diag_indices(nvars_diag))
    end select
    
    ! Set
    if (me == master) then
       select case (trim(adjustl(stream_name)))
       case ('restart')
          stream_list_restart_indices = indices
       case ('history')
          stream_list_history_indices = indices
       case ('diag_phys')
          stream_list_diag_indices = indices
       end select
    end if
    
    ! Broadcast data
    select case (trim(adjustl(stream_name)))
    case ('restart')
       call mpi_bcast(stream_list_restart_indices, nvars_rstr, MPI_INTEGER, master, mpicomm, mpierr)
    case ('history')
       call mpi_bcast(stream_list_history_indices, nvars_hist, MPI_INTEGER, master, mpicomm, mpierr)
    case ('diag_phys')
       call mpi_bcast(stream_list_diag_indices,    nvars_diag, MPI_INTEGER, master, mpicomm, mpierr)
    end select

  end subroutine read_stream_list
  
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
    character(len=*), parameter :: subname = 'ufs_mpas_io::ufs_mpas_write'
    character(len=:), allocatable :: filename
    integer :: ierr
    !type(var_info_type), allocatable :: history_var_info_list(:)
    integer :: timelevel, whence

    if (trim(stream_name) == "output") then
       filename = 'history.'//trim(timestamp)//'.nc'
    else if (trim(stream_name) == "output+diag_phys") then
       filename = 'history.'//trim(timestamp)//'.nc'
    else if (trim(stream_name) == "restart") then
       filename = 'restart.'//trim(timestamp)//'.nc'
    else if (trim(stream_name) == "input") then
       filename = 'input.'//trim(timestamp)//'.nc'
    else if (trim(stream_name) == "diag_phys") then
       filename = 'diag_phys.'//trim(timestamp)//'.nc'
    else
       stop "Invalid stream_name to ufs_mpas_write: stream_name =" &
            //trim(stream_name)
    end if

    if (debug) call mpas_log_write(subname // " entering ufs_mpas_write")
    if (debug) call mpas_log_write(subname // " creating "//trim(stream_name)//" stream file: "//trim(filename))
    ierr = pio_createfile(pio_subsystem_output, pioid_output, pio_iotype, trim(filename))
    if ( ierr /= 0 ) then
       call mpas_log_write(subname // " pio_createfile failed", messageType=MPAS_LOG_CRIT)
    endif

    !history_var_info_list = parse_stream_name_fragment('output')
    timelevel = TIMELEVEL_NOW
    whence = MPAS_NOW

    call dyn_mpas_read_write_stream(clock, "write", stream_name, pioid_output, &
                                    timeLevel=timelevel, whence=whence,        &
                                    nRecord=1, ierr=ierr, debug=debug)
    if ( ierr /= 0 ) then
       call mpas_log_write(subname // " dyn_mpas_read_write_stream failed ", &
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
   use module_mpas_config,     only : latCell, lonCell
   
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

   character(len=*), parameter :: subname = 'ufs_mpas_io::ufs_mpas_get_global_coords'


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
      call mpas_log_write(subname // " failed to allocate temp array", messageType=MPAS_LOG_CRIT)
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
 subroutine dyn_mpas_exchange_halo(field_name, debug, timeLevel)
   ! Module(s) from MPAS.
   use mpas_derived_types, only : field1dinteger, field2dinteger, field3dinteger,           &
                                  field1dreal, field2dreal, field3dreal, field4dreal,       &
                                  field5dreal, mpas_pool_field_info_type, mpas_pool_integer,&
                                  mpas_pool_real
   use mpas_dmpar,         only : mpas_dmpar_exch_halo_field
   use mpas_pool_routines, only : mpas_pool_get_field, mpas_pool_get_field_info
   ! Arguments
   character(*), intent(in)           :: field_name
   logical,      intent(in)           :: debug
   integer,      intent(in), optional :: timeLevel
   ! Locals
   integer :: tl
   type(field1dinteger), pointer :: field_1d_integer
   type(field2dinteger), pointer :: field_2d_integer
   type(field3dinteger), pointer :: field_3d_integer
   type(field1dreal), pointer :: field_1d_real
   type(field2dreal), pointer :: field_2d_real
   type(field3dreal), pointer :: field_3d_real
   type(field4dreal), pointer :: field_4d_real
   type(field5dreal), pointer :: field_5d_real
   type(mpas_pool_field_info_type) :: mpas_pool_field_info
   character(*), parameter :: subname = 'dyn_mpas_io::dyn_mpas_exchange_halo'

   if (debug) call mpas_log_write(subname // ' entered')

   tl = 1
   if (present(timeLevel)) tl = timeLevel

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
      call mpas_log_write(subname // ' Invalid field information for "' // trim(adjustl(field_name)) // '"', messageType=MPAS_LOG_CRIT)
   end if

   ! No halo layers to exchange. This field is not decomposed.
   if (mpas_pool_field_info % nhalolayers == 0) then
      call mpas_log_write(subname // ' Skipping field "' // trim(adjustl(field_name)) // '" due to not decomposed')
      return
   end if

   if (debug) call mpas_log_write(subname // 'Exchanging halo layers for "' // trim(adjustl(field_name)) // '"')

   select case (mpas_pool_field_info % fieldtype)
   case (mpas_pool_integer)
      select case (mpas_pool_field_info % ndims)
      case (1)
         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(field_name)), field_1d_integer, timelevel=tl)

         if (.not. associated(field_1d_integer)) then
            call mpas_log_write(subname // ' Failed to find field "' // trim(adjustl(field_name)) // '"', messageType=MPAS_LOG_CRIT)
         end if

         call mpas_dmpar_exch_halo_field(field_1d_integer)

         nullify(field_1d_integer)
      case (2)
         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(field_name)), field_2d_integer, timelevel=tl)

         if (.not. associated(field_2d_integer)) then
            call mpas_log_write(subname // ' Failed to find field "' // trim(adjustl(field_name)) // '"', messageType=MPAS_LOG_CRIT)
         end if

         call mpas_dmpar_exch_halo_field(field_2d_integer)

         nullify(field_2d_integer)
      case (3)
         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(field_name)), field_3d_integer, timelevel=tl)

         if (.not. associated(field_3d_integer)) then
            call mpas_log_write(subname // ' Failed to find field "' // trim(adjustl(field_name)) // '"', messageType=MPAS_LOG_CRIT)
         end if

         call mpas_dmpar_exch_halo_field(field_3d_integer)

         nullify(field_3d_integer)
      case default
         call mpas_log_write(subname // ' Unsupported field rank "' // trim(adjustl(field_name)) // '"', messageType=MPAS_LOG_CRIT)
      end select
   case (mpas_pool_real)
      select case (mpas_pool_field_info % ndims)
      case (1)
         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(field_name)), field_1d_real, timelevel=tl)

         if (.not. associated(field_1d_real)) then
            call mpas_log_write(subname // ' Failed to find field "' // trim(adjustl(field_name)) // '"', messageType=MPAS_LOG_CRIT)
         end if

         call mpas_dmpar_exch_halo_field(field_1d_real)

         nullify(field_1d_real)
      case (2)
         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(field_name)), field_2d_real, timelevel=tl)

         if (.not. associated(field_2d_real)) then
            call mpas_log_write(subname // ' Failed to find field "' // trim(adjustl(field_name)) // '"', messageType=MPAS_LOG_CRIT)
         end if
         call mpas_dmpar_exch_halo_field(field_2d_real)
         nullify(field_2d_real)
      case (3)
         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(field_name)), field_3d_real, timelevel=tl)

         if (.not. associated(field_3d_real)) then
            call mpas_log_write(subname // ' Failed to find field "' // trim(adjustl(field_name)) // '"', messageType=MPAS_LOG_CRIT)
         end if

         call mpas_dmpar_exch_halo_field(field_3d_real)

         nullify(field_3d_real)
      case (4)
         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(field_name)), field_4d_real, timelevel=tl)

         if (.not. associated(field_4d_real)) then
            call mpas_log_write(subname // ' Failed to find field "' // trim(adjustl(field_name)) // '"', messageType=MPAS_LOG_CRIT)
         end if

         call mpas_dmpar_exch_halo_field(field_4d_real)

         nullify(field_4d_real)
      case (5)
         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(field_name)), field_5d_real, timelevel=tl)

         if (.not. associated(field_5d_real)) then
            call mpas_log_write(subname // ' Failed to find field "' // trim(adjustl(field_name)) // '"', messageType=MPAS_LOG_CRIT)
         end if

         call mpas_dmpar_exch_halo_field(field_5d_real)

         nullify(field_5d_real)
      case default
         call mpas_log_write(subname // ' Unsupported field rank "' // trim(adjustl(field_name)) // '"', messageType=MPAS_LOG_CRIT)
      end select
   case default
      call mpas_log_write(subname // ' Unsupported field type (Must be one of: integer, real)', messageType=MPAS_LOG_CRIT)
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
   character(*), parameter :: subname = 'dyn_mpas_io::dyn_mpas_read_write_stream'
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
      call mpas_log_write(subname // ' Failed to initialize stream "' // trim(adjustl(stream_name)) // '"', messageType=MPAS_LOG_CRIT)
   end if

   if (.not. associated(mpas_stream)) then
      call mpas_log_write(subname // ' Failed to initialize stream "' // trim(adjustl(stream_name)) // '"', messageType=MPAS_LOG_CRIT)
   end if

   select case (trim(adjustl(stream_mode)))
   case ('r', 'read')
      if (debug) call mpas_log_write(subname // 'Reading stream "' // trim(adjustl(stream_name)) // '"')

      call read_stream(mpas_stream, actualWhen, nRecord, ierr)

      if (ierr /= mpas_stream_noerr) then
         call mpas_log_write(subname // ' Failed to initialize stream "' // trim(adjustl(stream_name)) // '"', messageType=MPAS_LOG_CRIT)
      end if

      ! Exchange halo layers because new data have just been read.
      var_info_list = parse_stream_name(stream_name)

      do i = 1, size(var_info_list)
         call dyn_mpas_exchange_halo(var_info_list(i) % name, debug, timeLevel=timeLevel)
         if ( ierr /= 0 ) then
            call mpas_log_write(subname // ' Failed to exchange halo layers for group '//var_info_list(i) % name, messageType=MPAS_LOG_CRIT)
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
         call mpas_log_write(subname // ' Failed to write stream "' // trim(adjustl(stream_name)) // '"', messageType=MPAS_LOG_CRIT)
      end if

      ! For any connectivity arrays in this stream, reset global indexes back to local indexes.
      call postwrite_reindex(domain_ptr % blocklist % allfields, mpas_pool)
   case default
      call mpas_log_write(subname // ' Unsupported stream mode "' // trim(adjustl(stream_mode)) // '"', messageType=MPAS_LOG_CRIT)
   end select

   if (debug) call mpas_log_write('Closing stream "' // trim(adjustl(stream_name)) // '"')
   call mpas_log_write( '---------------------------------------------------------------------')

   call mpas_closestream(mpas_stream, ierr=ierr)

   if (ierr /= mpas_stream_noerr) then
      call mpas_log_write(subname // ' Failed to close stream "' // trim(adjustl(stream_name)) // '"', messageType=MPAS_LOG_CRIT)
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
      call mpas_log_write('MPAS_streamTime ierr = $i, nRecord = $i', intArgs=[ierr, nRecord])
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

   character(*), parameter :: subname = 'dyn_mpas_io::dyn_mpas_init_stream_with_pool'
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
      call mpas_log_write(subname // ' Failed to allocate stream "' // trim(adjustl(stream_name)) // '"', messageType=MPAS_LOG_CRIT)
   end if

   ! Not actually used because a PIO file descriptor is directly supplied.
   stream_filename = 'external stream'
   stream_format = mpas_io_pnetcdf

   if (debug) call mpas_log_write('Checking PIO file descriptor')

   if (.not. associated(pio_file)) then
      call mpas_log_write(subname // ' Invalid PIO file descriptor', messageType=MPAS_LOG_CRIT)
   end if

   if (.not. pio_file_is_open(pio_file)) then
      call mpas_log_write(subname // ' Invalid PIO file descriptor', messageType=MPAS_LOG_CRIT)
   end if

   select case (trim(adjustl(stream_mode)))
   case ('r', 'read')
      if (debug) call mpas_log_write(' Creating stream "' // trim(adjustl(stream_name)) // '" for reading')

      call mpas_createstream( &
           mpas_stream, domain_ptr % iocontext, stream_filename, stream_format, mpas_io_read,  &
           clobberrecords=.false., clobberfiles=.false., truncatefiles=.false., &
           precision=mpas_io_native_precision, pio_file_desc=pio_file, ierr=ierr)
   case ('w', 'write')
      if (debug) call mpas_log_write(' Creating stream "' // trim(adjustl(stream_name)) // '" for writing')

      call mpas_createstream( &
           mpas_stream, domain_ptr % iocontext, stream_filename, stream_format, mpas_io_write, &
           clobberrecords=.false., clobberfiles=.false., truncatefiles=.false., &
           precision=mpas_io_native_precision, pio_file_desc=pio_file, ierr=ierr)
   case default
      call mpas_log_write(subname // ' Unsupported stream mode "' // trim(adjustl(stream_mode)) // '"', messageType=MPAS_LOG_CRIT)
   end select

   if (ierr /= mpas_stream_noerr) then
      call mpas_log_write(subname // ' Failed to create stream "' // trim(adjustl(stream_name)) // '"', messageType=MPAS_LOG_CRIT)
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
               call mpas_log_write(subname // ' Failed to find variable "' //        &
                                  trim(adjustl(var_info_list(i) % name)) // '"',     &
                                  messageType=MPAS_LOG_CRIT)
            end if

            call mpas_streamaddfield(mpas_stream, field_0d_char, ierr=ierr)

            nullify(field_0d_char)
         case (1)
            call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
                 trim(adjustl(var_info_list(i) % name)), field_1d_char, timelevel=timeLevel)

            if (.not. associated(field_1d_char)) then
               call mpas_log_write(subname // ' Failed to find variable "' //        &
                                   trim(adjustl(var_info_list(i) % name)) // '"',    &
                                   messageType=MPAS_LOG_CRIT)
            end if

            call mpas_streamaddfield(mpas_stream, field_1d_char, ierr=ierr)

            nullify(field_1d_char)
         case default
            call mpas_log_write(subname // ' Unsupported variable rank ' //                &
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
               call mpas_log_write(subname // ' Failed to find variable "' //         &
                                   trim(adjustl(var_info_list(i) % name)) // '"',     &
                                   messageType=MPAS_LOG_CRIT)
            end if

            call mpas_streamaddfield(mpas_stream, field_0d_integer, ierr=ierr)

            nullify(field_0d_integer)
         case (1)
            call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
                 trim(adjustl(var_info_list(i) % name)), field_1d_integer, timelevel=timeLevel)

            if (.not. associated(field_1d_integer)) then
               call mpas_log_write(subname // ' Failed to find variable "' //        &
                                   trim(adjustl(var_info_list(i) % name)) // '"',    &
                                   messageType=MPAS_LOG_CRIT)
            end if

            call mpas_streamaddfield(mpas_stream, field_1d_integer, ierr=ierr)

            nullify(field_1d_integer)
         case (2)
            call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
                 trim(adjustl(var_info_list(i) % name)), field_2d_integer, timelevel=timeLevel)

            if (.not. associated(field_2d_integer)) then
               call mpas_log_write(subname // ' Failed to find variable "' //        &
                                   trim(adjustl(var_info_list(i) % name)) // '"',    &
                                   messageType=MPAS_LOG_CRIT)
            end if

            call mpas_streamaddfield(mpas_stream, field_2d_integer, ierr=ierr)

            nullify(field_2d_integer)
         case (3)
            call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
                 trim(adjustl(var_info_list(i) % name)), field_3d_integer, timelevel=timeLevel)

            if (.not. associated(field_3d_integer)) then
               call mpas_log_write(subname // ' Failed to find variable "' //        &
                                   trim(adjustl(var_info_list(i) % name)) // '"',    &
                                   messageType=MPAS_LOG_CRIT)
            end if

            call mpas_streamaddfield(mpas_stream, field_3d_integer, ierr=ierr)

            nullify(field_3d_integer)
         case default
            call mpas_log_write(subname // ' Unsupported variable rank ' //                &
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
               call mpas_log_write(subname // ' Failed to find variable "' //         &
                                   trim(adjustl(var_info_list(i) % name)) // '"',     &
                                   messageType=MPAS_LOG_CRIT)
            end if

            call mpas_streamaddfield(mpas_stream, field_0d_real, ierr=ierr)

            nullify(field_0d_real)
         case (1)
            call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
                 trim(adjustl(var_info_list(i) % name)), field_1d_real, timelevel=timeLevel)

            if (.not. associated(field_1d_real)) then
               call mpas_log_write(subname // ' Failed to find variable "' //        &
                                   trim(adjustl(var_info_list(i) % name)) // '"',    &
                                   messageType=MPAS_LOG_CRIT)
            end if

            call mpas_streamaddfield(mpas_stream, field_1d_real, ierr=ierr)

            nullify(field_1d_real)
         case (2)
            call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
                 trim(adjustl(var_info_list(i) % name)), field_2d_real, timelevel=timeLevel)
            if (.not. associated(field_2d_real)) then
               call mpas_log_write(subname // ' Failed to find variable "' //        &
                                   trim(adjustl(var_info_list(i) % name)) // '"',    &
                                   messageType=MPAS_LOG_CRIT)
            end if
            call mpas_streamaddfield(mpas_stream, field_2d_real, ierr=ierr)

            nullify(field_2d_real)
         case (3)
            call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
                 trim(adjustl(var_info_list(i) % name)), field_3d_real, timelevel=timeLevel)

            if (.not. associated(field_3d_real)) then
               call mpas_log_write(subname // ' Failed to find variable "' //        &
                                   trim(adjustl(var_info_list(i) % name)) // '"',    &
                                   messageType=MPAS_LOG_CRIT)
            end if
            call mpas_streamaddfield(mpas_stream, field_3d_real, ierr=ierr)

            nullify(field_3d_real)
         case (4)
            call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
                 trim(adjustl(var_info_list(i) % name)), field_4d_real, timelevel=timeLevel)

            if (.not. associated(field_4d_real)) then
               call mpas_log_write(subname // ' Failed to find variable "' //        &
                                   trim(adjustl(var_info_list(i) % name)) // '"',    &
                                   messageType=MPAS_LOG_CRIT)
            end if

            call mpas_streamaddfield(mpas_stream, field_4d_real, ierr=ierr)

            nullify(field_4d_real)
         case (5)
            call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
                 trim(adjustl(var_info_list(i) % name)), field_5d_real, timelevel=timeLevel)

            if (.not. associated(field_5d_real)) then
               call mpas_log_write(subname // ' Failed to find variable "' //        &
                                   trim(adjustl(var_info_list(i) % name)) // '"',    &
                                   messageType=MPAS_LOG_CRIT)
            end if

            call mpas_streamaddfield(mpas_stream, field_5d_real, ierr=ierr)

            nullify(field_5d_real)
         case default
            call mpas_log_write(subname // ' Unsupported variable rank ' //                &
                                stringify([var_info_list(i) % rank]) //                    &
                                ' for "' // trim(adjustl(var_info_list(i) % name)) // '"', &
                                messageType=MPAS_LOG_CRIT)
         end select
      case default
         call mpas_log_write(subname // ' Unsupported variable type "' //                &
                             trim(adjustl(var_info_list(i) % type)) //                   &
                             '" for "' // trim(adjustl(var_info_list(i) % name)) // '"', &
                             messageType=MPAS_LOG_CRIT)
      end select

      if (ierr /= mpas_stream_noerr) then
         call mpas_log_write(subname // ' Failed to add variable "' //             &
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
        call mpas_log_write(subname // ' Unsupported attribute type (Must be one of: character, integer, logical, real)', &
                            messageType=MPAS_LOG_CRIT)
     end select

     if (ierr /= mpas_stream_noerr) then
        call mpas_log_write(subname // ' Failed to add attribute "' //            &
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
        call mpas_log_write(subname // ' Unsupported attribute type (Must be one of: integer, real)',&
                            messageType=MPAS_LOG_CRIT)
     end select

     if (ierr /= mpas_stream_noerr) then
        call mpas_log_write(subname // ' Failed to add attribute "' //            &
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
      ! If stream_list provided at runtime, only include requested fields.
      if (allocated(stream_list_restart_indices)) then
         allocate(var_info_list, source=restart_var_info_list(stream_list_restart_indices))
      ! Otherwise, include all available fields from stream (default).
      else
         allocate(var_info_list, source=restart_var_info_list)
      end if
   case ('output')
      ! If stream_list provided at runtime, only include requested fields.
      if (allocated(stream_list_history_indices)) then
         allocate(var_info_list, source=history_var_info_list(stream_list_history_indices))
      ! Otherwise, include all available fields from stream (default).
      else
         allocate(var_info_list, source=history_var_info_list)
      end if
   case ('lbc_in')
      allocate(var_info_list, source=lbc_in_var_info_list)
   case ('sfc_input')
      allocate(var_info_list, source=sfc_input_var_info_list)
   case ('ugwp_oro_data')
      allocate(var_info_list, source=ugwp_oro_data_var_info_list)
   case ('diag_phys')
      ! If stream_list provided at runtime, only include requested fields.
      if (allocated(stream_list_diag_indices)) then
         allocate(var_info_list, source=diag_phys_var_info_list(stream_list_diag_indices))
         ! Otherwise, include all available fields from stream (default).
      else
         allocate(var_info_list, source=diag_phys_var_info_list)
      end if
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

      var_name_list = history_var_info_list % name

      if (any(var_name_list == trim(adjustl(stream_name_fragment)))) then
         var_info_list_buffer = pack(history_var_info_list, var_name_list == trim(adjustl(stream_name_fragment)))
         var_info_list = [var_info_list, var_info_list_buffer]
      end if

      var_name_list = lbc_in_var_info_list % name

      if (any(var_name_list == trim(adjustl(stream_name_fragment)))) then
         var_info_list_buffer = pack(lbc_in_var_info_list, var_name_list == trim(adjustl(stream_name_fragment)))
         var_info_list = [var_info_list, var_info_list_buffer]
      end if

      var_name_list = ugwp_oro_data_var_info_list % name

      if (any(var_name_list == trim(adjustl(stream_name_fragment)))) then
         var_info_list_buffer = pack(ugwp_oro_data_var_info_list, var_name_list == trim(adjustl(stream_name_fragment)))
         var_info_list = [var_info_list, var_info_list_buffer]
      end if

      var_name_list = diag_phys_var_info_list % name

      if (any(var_name_list == trim(adjustl(stream_name_fragment)))) then
         var_info_list_buffer = pack(diag_phys_var_info_list, var_name_list == trim(adjustl(stream_name_fragment)))
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

   logical, allocatable, intent(out) :: var_is_present(:)
   logical, allocatable, intent(out) :: var_is_tkr_compatible(:)
   type(file_desc_t), pointer, intent(in) :: pio_file
   type(var_info_type), intent(in) :: var_info
   logical, intent(in) :: debug

   character(*), parameter :: subname = 'dyn_mpas_io::dyn_mpas_check_variable_status'
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
            call mpas_log_write(subname // ' Failed to find variable "' //        &
                                trim(adjustl(var_info % name)) // '"',            &
                                messageType=MPAS_LOG_CRIT)
         end if

         if (field_0d_char % isvararray .and. associated(field_0d_char % constituentnames)) then
            allocate(var_name_list(size(field_0d_char % constituentnames)), stat=ierr)

            if (ierr /= 0) then
               call mpas_log_write(subname // ' Failed to allocate var_name_list', &
                                   messageType=MPAS_LOG_CRIT)
            end if

            var_name_list(:) = field_0d_char % constituentnames(:)
         end if

         nullify(field_0d_char)
      case (1)
         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(var_info % name)), field_1d_char, timelevel=1)

         if (.not. associated(field_1d_char)) then
            call mpas_log_write(subname // ' Failed to find variable "' //        &
                                trim(adjustl(var_info % name)) // '"',            &
                                messageType=MPAS_LOG_CRIT)
         end if

         if (field_1d_char % isvararray .and. associated(field_1d_char % constituentnames)) then
            allocate(var_name_list(size(field_1d_char % constituentnames)), stat=ierr)

            if (ierr /= 0) then
               call mpas_log_write(subname // ' Failed to allocate var_name_list', &
                                   messageType=MPAS_LOG_CRIT)
            end if

            var_name_list(:) = field_1d_char % constituentnames(:)
         end if

         nullify(field_1d_char)
      case default
         call mpas_log_write(subname // ' Unsupported variable rank ' //        &
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
            call mpas_log_write(subname // ' Failed to find variable "' //        &
                                trim(adjustl(var_info % name)) // '"',            &
                                messageType=MPAS_LOG_CRIT)
         end if

         if (field_0d_integer % isvararray .and. associated(field_0d_integer % constituentnames)) then
            allocate(var_name_list(size(field_0d_integer % constituentnames)), stat=ierr)

            if (ierr /= 0) then
               call mpas_log_write(subname // ' Failed to allocate var_name_list', &
                                   messageType=MPAS_LOG_CRIT)
            end if

            var_name_list(:) = field_0d_integer % constituentnames(:)
         end if

         nullify(field_0d_integer)
      case (1)
         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(var_info % name)), field_1d_integer, timelevel=1)

         if (.not. associated(field_1d_integer)) then
            call mpas_log_write(subname // ' Failed to find variable "' //        &
                                trim(adjustl(var_info % name)) // '"',            &
                                messageType=MPAS_LOG_CRIT)
         end if

         if (field_1d_integer % isvararray .and. associated(field_1d_integer % constituentnames)) then
            allocate(var_name_list(size(field_1d_integer % constituentnames)), stat=ierr)

            if (ierr /= 0) then
               call mpas_log_write(subname // ' Failed to allocate var_name_list', &
                                   messageType=MPAS_LOG_CRIT)
            end if

            var_name_list(:) = field_1d_integer % constituentnames(:)
         end if

         nullify(field_1d_integer)
      case (2)
         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(var_info % name)), field_2d_integer, timelevel=1)

         if (.not. associated(field_2d_integer)) then
            call mpas_log_write(subname // ' Failed to find variable "' //        &
                                trim(adjustl(var_info % name)) // '"',            &
                                messageType=MPAS_LOG_CRIT)
         end if

         if (field_2d_integer % isvararray .and. associated(field_2d_integer % constituentnames)) then
            allocate(var_name_list(size(field_2d_integer % constituentnames)), stat=ierr)

            if (ierr /= 0) then
               call mpas_log_write(subname // ' Failed to allocate var_name_list', &
                                   messageType=MPAS_LOG_CRIT)
            end if

            var_name_list(:) = field_2d_integer % constituentnames(:)
         end if

         nullify(field_2d_integer)
      case (3)
         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(var_info % name)), field_3d_integer, timelevel=1)

         if (.not. associated(field_3d_integer)) then
            call mpas_log_write(subname // ' Failed to find variable "' // &
                                trim(adjustl(var_info % name)) // '"',            &
                                messageType=MPAS_LOG_CRIT)
         end if

         if (field_3d_integer % isvararray .and. associated(field_3d_integer % constituentnames)) then
            allocate(var_name_list(size(field_3d_integer % constituentnames)), stat=ierr)

            if (ierr /= 0) then
               call mpas_log_write(subname // ' Failed to allocate var_name_list', &
                                   messageType=MPAS_LOG_CRIT)
            end if

            var_name_list(:) = field_3d_integer % constituentnames(:)
         end if

         nullify(field_3d_integer)
      case default
         call mpas_log_write(subname // ' Unsupported variable rank ' // &
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
            call mpas_log_write(subname // ' Failed to find variable "' //        &
                                trim(adjustl(var_info % name)) // '"',            &
                                messageType=MPAS_LOG_CRIT)
         end if

         if (field_0d_real % isvararray .and. associated(field_0d_real % constituentnames)) then
            allocate(var_name_list(size(field_0d_real % constituentnames)), stat=ierr)

            if (ierr /= 0) then
               call mpas_log_write(subname // ' Failed to allocate var_name_list', &
                                   messageType=MPAS_LOG_CRIT)
            end if

            var_name_list(:) = field_0d_real % constituentnames(:)
         end if

         nullify(field_0d_real)
      case (1)

         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(var_info % name)), field_1d_real, timelevel=1)

         if (.not. associated(field_1d_real)) then
            call mpas_log_write(subname // ' Failed to find variable "' //        &
                                trim(adjustl(var_info % name)) // '"',            &
                                messageType=MPAS_LOG_CRIT)
         end if

         if (field_1d_real % isvararray .and. associated(field_1d_real % constituentnames)) then
            allocate(var_name_list(size(field_1d_real % constituentnames)), stat=ierr)

            if (ierr /= 0) then
               call mpas_log_write(subname // ' Failed to allocate var_name_list', &
                                   messageType=MPAS_LOG_CRIT)
            end if

            var_name_list(:) = field_1d_real % constituentnames(:)
         end if

         nullify(field_1d_real)
      case (2)
         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(var_info % name)), field_2d_real, timelevel=1)

         if (.not. associated(field_2d_real)) then
            call mpas_log_write(subname // ' Failed to find variable "' //        &
                                trim(adjustl(var_info % name)) // '"',            &
                                messageType=MPAS_LOG_CRIT)
         end if

         if (field_2d_real % isvararray .and. associated(field_2d_real % constituentnames)) then
            allocate(var_name_list(size(field_2d_real % constituentnames)), stat=ierr)

            if (ierr /= 0) then
               call mpas_log_write(subname // ' Failed to allocate var_name_list', &
                                   messageType=MPAS_LOG_CRIT)
            end if

            var_name_list(:) = field_2d_real % constituentnames(:)
         end if

         nullify(field_2d_real)
      case (3)
         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(var_info % name)), field_3d_real, timelevel=1)
         if (.not. associated(field_3d_real)) then
            call mpas_log_write(subname // ' Failed to find variable "' //        &
                                trim(adjustl(var_info % name)) // '"',            &
                                messageType=MPAS_LOG_CRIT)
         end if
         if (field_3d_real % isvararray .and. associated(field_3d_real % constituentnames)) then
            allocate(var_name_list(size(field_3d_real % constituentnames)), stat=ierr)

            if (ierr /= 0) then
               call mpas_log_write(subname // ' Failed to allocate var_name_list', &
                                   messageType=MPAS_LOG_CRIT)
            end if

            var_name_list(:) = field_3d_real % constituentnames(:)
         end if
         nullify(field_3d_real)
      case (4)
         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(var_info % name)), field_4d_real, timelevel=1)

         if (.not. associated(field_4d_real)) then
            call mpas_log_write(subname // ' Failed to find variable "' //        &
                                trim(adjustl(var_info % name)) // '"',            &
                                messageType=MPAS_LOG_CRIT)
         end if

         if (field_4d_real % isvararray .and. associated(field_4d_real % constituentnames)) then
            allocate(var_name_list(size(field_4d_real % constituentnames)), stat=ierr)

            if (ierr /= 0) then
               call mpas_log_write(subname // ' Failed to allocate var_name_list', &
                                   messageType=MPAS_LOG_CRIT)
            end if

            var_name_list(:) = field_4d_real % constituentnames(:)
         end if

         nullify(field_4d_real)
      case (5)
         call mpas_pool_get_field(domain_ptr % blocklist % allfields, &
              trim(adjustl(var_info % name)), field_5d_real, timelevel=1)

         if (.not. associated(field_5d_real)) then
            call mpas_log_write(subname // ' Failed to find variable "' //        &
                                trim(adjustl(var_info % name)) // '"',            &
                                messageType=MPAS_LOG_CRIT)
         end if

         if (field_5d_real % isvararray .and. associated(field_5d_real % constituentnames)) then
            allocate(var_name_list(size(field_5d_real % constituentnames)), stat=ierr)

            if (ierr /= 0) then
               call mpas_log_write(subname // ' Failed to allocate var_name_list', &
                                   messageType=MPAS_LOG_CRIT)
            end if

            var_name_list(:) = field_5d_real % constituentnames(:)
         end if

         nullify(field_5d_real)
      case default
         call mpas_log_write(subname // ' Unsupported variable rank ' //        &
                             stringify([var_info % rank]) //                    &
                             ' for "' // trim(adjustl(var_info % name)) // '"', &
                             messageType=MPAS_LOG_CRIT)
      end select
   case default
      call mpas_log_write(subname // ' Unsupported variable type ' //        &
                          stringify([var_info % type]) //                    &
                          ' for "' // trim(adjustl(var_info % name)) // '"', &
                          messageType=MPAS_LOG_CRIT)
   end select

   if (.not. allocated(var_name_list)) then
      allocate(var_name_list(1), stat=ierr)

      if (ierr /= 0) then
         call mpas_log_write(subname // ' Failed to allocate var_name_list', &
                             messageType=MPAS_LOG_CRIT)
      end if

      var_name_list(1) = var_info % name
   end if

   allocate(var_is_present(size(var_name_list)), stat=ierr)

   if (ierr /= 0) then
      call mpas_log_write(subname // ' Failed to allocate var_is_present', &
                          messageType=MPAS_LOG_CRIT)
   end if

   var_is_present(:) = .false.
   allocate(var_is_tkr_compatible(size(var_name_list)), stat=ierr)
   if (ierr /= 0) then
      call mpas_log_write(subname // ' Failed to allocate var_is_tkr_compatible', &
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
