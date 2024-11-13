module GFS_typedefs

   use mpi_f08
   use machine,                  only: kind_phys, kind_dbl_prec, kind_sngl_prec
   use physcons,                 only: con_cp, con_fvirt, con_g, rholakeice,           &
                                       con_hvap, con_hfus, con_pi, con_rd, con_rv,     &
                                       con_t0c, con_cvap, con_cliq, con_eps, con_epsq, &
                                       con_epsm1, con_ttp, rlapse, con_jcal, con_rhw0, &
                                       con_sbc, con_tice, cimin, con_p0, rhowater,     &
                                       con_csol, con_epsqs, con_rocp, con_rog,         &
                                       con_omega, con_rerth, con_psat, karman, rainmin,&
                                       con_c, con_plnk, con_boltz, con_solr_2008,      &
                                       con_solr_2002, con_thgni, con_1ovg

   use module_radsw_parameters,  only: topfsw_type, sfcfsw_type
   use module_radlw_parameters,  only: topflw_type, sfcflw_type
   use module_ozphys,            only: ty_ozphys
   use module_h2ophys,           only: ty_h2ophys

   implicit none

   ! To ensure that these values match what's in the physics, array
   ! sizes are compared in the auto-generated physics caps in debug mode
   ! from aerclm_def
   integer, parameter, private :: ntrcaerm = 15

   ! This will be set later in GFS_Control%initialize, since
   ! it depends on the runtime config (Model%aero_in)
   integer, private  :: ntrcaer

   ! If these are changed to >99, need to adjust formatting string in GFS_diagnostics.F90 (and names in diag_tables)
   integer, parameter :: naux2dmax = 20 !< maximum number of auxiliary 2d arrays in output (for debugging)
   integer, parameter :: naux3dmax = 20 !< maximum number of auxiliary 3d arrays in output (for debugging)

   integer, parameter :: dfi_radar_max_intervals = 4 !< Number of radar-derived temperature tendency and/or convection suppression intervals. Do not change.

   real(kind=kind_phys), parameter :: limit_unspecified = 1e12 !< special constant for "namelist value was not provided" in radar-derived temperature tendency limit range


!> \section arg_table_GFS_typedefs
!! \htmlinclude GFS_typedefs.html
!!

  !--- version of physics
  character(len=64) :: phys_version = 'v2021 UFS PHYSICS'

  !--- parameter constants used for default initializations
  real(kind=kind_phys), parameter :: zero      = 0.0_kind_phys
  !real(kind=kind_phys), parameter :: huge      = 9.9692099683868690E36 ! NetCDF float FillValue
  real(kind=kind_phys), parameter :: clear_val = zero
  !real(kind=kind_phys), parameter :: clear_val = -9.9999e80
  real(kind=kind_phys), parameter :: rann_init = 0.6_kind_phys
  real(kind=kind_phys), parameter :: cn_one    = 1._kind_phys
  real(kind=kind_phys), parameter :: cn_100    = 100._kind_phys
  real(kind=kind_phys), parameter :: cn_th     = 1000._kind_phys
  real(kind=kind_phys), parameter :: cn_hr     = 3600._kind_phys

  ! optional extra top layer on top of low ceiling models
  ! this parameter was originally defined in the radiation driver
  ! (and is still for standard non-CCPP builds), but is required
  ! here for CCPP to allocate arrays used for the interstitial
  ! calculations previously in GFS_{physics,radiation}_driver.F90
  ! LTP=0: no extra top layer
  integer, parameter :: LTP = 0   ! no extra top layer
  !integer, parameter :: LTP = 1   ! add an extra top layer

!----------------
! Data Containers
!----------------
!    !--- GFS external initialization type
!    GFS_init_type
!    !--- GFS Derived Data Types (DDTs)
!    GFS_statein_type        !< prognostic state data in from dycore
!    GFS_stateout_type       !< prognostic state or tendencies return to dycore
!    GFS_sfcprop_type        !< surface fields
!    GFS_coupling_type       !< fields to/from coupling with other components (e.g. land/ice/ocean/etc.)
!    !---GFS specific containers
!    GFS_control_type        !< model control parameters
!    GFS_grid_type           !< grid and interpolation related data
!    GFS_tbd_type            !< to be determined data that doesn't fit in any one container
!    GFS_cldprop_type        !< cloud fields needed by radiation from physics
!    GFS_radtend_type        !< radiation tendencies needed in physics
!    GFS_diag_type           !< fields targetted for diagnostic output

!--------------------------------------------------------------------------------
! GFS_init_type
!--------------------------------------------------------------------------------
!   This container is the minimum set of data required from the dycore/atmosphere
!   component to allow proper initialization of the GFS physics
!--------------------------------------------------------------------------------
!! \section arg_table_GFS_init_type
!! \htmlinclude GFS_init_type.html
!!
  type GFS_init_type
    integer :: me                                !< my MPI-rank
    integer :: master                            !< master MPI-rank
    type(MPI_Comm) :: fcst_mpi_comm              !< forecast tasks mpi communicator
    integer :: fcst_ntasks                       !< total number of forecast tasks
    integer :: tile_num                          !< tile number for this MPI rank
    integer :: isc                               !< starting i-index for this MPI-domain
    integer :: jsc                               !< starting j-index for this MPI-domain
    integer :: nx                                !< number of points in i-dir for this MPI rank
    integer :: ny                                !< number of points in j-dir for this MPI rank
    integer :: levs                              !< number of vertical levels
    integer :: cnx                               !< number of points in i-dir for this cubed-sphere face
                                                 !< equal to gnx for lat-lon grids
    integer :: cny                               !< number of points in j-dir for this cubed-sphere face
                                                 !< equal to gny for lat-lon grids
    integer :: gnx                               !< number of global points in x-dir (i) along the equator
    integer :: gny                               !< number of global points in y-dir (j) along any meridian
    integer :: nlunit                            !< fortran unit number for file opens
    integer :: logunit                           !< fortran unit number for writing logfile
    integer :: bdat(8)                           !< model begin date in GFS format   (same as idat)
    integer :: cdat(8)                           !< model current date in GFS format (same as jdat)
    integer :: iau_offset                        !< iau running window length
    real(kind=kind_phys) :: dt_dycore            !< dynamics time step in seconds
    real(kind=kind_phys) :: dt_phys              !< physics  time step in seconds
!--- restart information
    logical :: restart                           !< flag whether this is a coldstart (.false.) or a warmstart/restart (.true.)
!--- hydrostatic/non-hydrostatic flag
    logical :: hydrostatic                       !< flag whether this is a hydrostatic or non-hydrostatic run
!--- blocking data
    integer, pointer :: blksz(:)                 !< for explicit data blocking
                                                 !< default blksz(1)=[nx*ny]
!--- ak/bk for pressure level calculations
    real(kind=kind_phys), pointer :: ak(:)       !< from surface (k=1) to TOA (k=levs)
    real(kind=kind_phys), pointer :: bk(:)       !< from surface (k=1) to TOA (k=levs)
!--- grid metrics
    real(kind=kind_phys), pointer :: xlon(:,:)   !< column longitude for MPI rank
    real(kind=kind_phys), pointer :: xlat(:,:)   !< column latitude  for MPI rank
    real(kind=kind_phys), pointer :: area(:,:)   !< column area for length scale calculations

    integer                    :: nwat            !< number of hydrometeors in dcyore (including water vapor)
    character(len=32), pointer :: tracer_names(:) !< tracers names to dereference tracer id
    integer,           pointer :: tracer_types(:) !< tracers types: 0=generic, 1=chem,prog, 2=chem,diag
    character(len=64) :: fn_nml                   !< namelist filename
    character(len=:), pointer, dimension(:) :: input_nml_file => null() !< character string containing full namelist
                                                                        !< for use with internal file reads
 end type GFS_init_type


!----------------------------------------------------------------
! GFS_statein_type
!   prognostic state variables with layer and level specific data
!----------------------------------------------------------------
!! \section arg_table_GFS_statein_type
!! \htmlinclude GFS_statein_type.html
!!
  type GFS_statein_type

!--- level geopotential and pressures
    real (kind=kind_phys), pointer :: phii  (:,:) => null()   !< interface geopotential height
    real (kind=kind_phys), pointer :: prsi  (:,:) => null()   !< model level pressure in Pa
    real (kind=kind_phys), pointer :: prsik (:,:) => null()   !< Exner function at interface

!--- layer geopotential and pressures
    real (kind=kind_phys), pointer :: phil  (:,:) => null()   !< layer geopotential height
    real (kind=kind_phys), pointer :: prsl  (:,:) => null()   !< model layer mean pressure Pa
    real (kind=kind_phys), pointer :: prslk (:,:) => null()   !< exner function = (p/p0)**rocp

!--- prognostic variables
    real (kind=kind_phys), pointer :: pgr  (:)     => null()  !< surface pressure (Pa) real
    real (kind=kind_phys), pointer :: ugrs (:,:)   => null()  !< u component of layer wind
    real (kind=kind_phys), pointer :: vgrs (:,:)   => null()  !< v component of layer wind
    real (kind=kind_phys), pointer :: wgrs (:,:)   => null()  !< w component of layer wind
    real (kind=kind_phys), pointer :: vvl  (:,:)   => null()  !< layer mean vertical velocity in pa/sec
    real (kind=kind_phys), pointer :: tgrs (:,:)   => null()  !< model layer mean temperature in k
    real (kind=kind_phys), pointer :: qgrs (:,:,:) => null()  !< layer mean tracer concentration
! dissipation estimate
    real (kind=kind_phys), pointer :: diss_est(:,:)   => null()  !< model layer mean temperature in k
    ! soil state variables - for soil SPPT - sfc-perts, mgehne
    real (kind=kind_phys), pointer :: smc (:,:)   => null()  !< soil moisture content
    real (kind=kind_phys), pointer :: stc (:,:)   => null()  !< soil temperature content
    real (kind=kind_phys), pointer :: slc (:,:)   => null()  !< soil liquid water content

    contains
      procedure :: create  => statein_create  !<   allocate array data
  end type GFS_statein_type


!------------------------------------------------------------------
! GFS_stateout_type
!   prognostic state or tendencies after physical parameterizations
!------------------------------------------------------------------
!! \section arg_table_GFS_stateout_type
!! \htmlinclude GFS_stateout_type.html
!!
  type GFS_stateout_type

    !-- Out (physics only)
    real (kind=kind_phys), pointer :: gu0 (:,:)   => null()  !< updated zonal wind
    real (kind=kind_phys), pointer :: gv0 (:,:)   => null()  !< updated meridional wind
    real (kind=kind_phys), pointer :: gt0 (:,:)   => null()  !< updated temperature
    real (kind=kind_phys), pointer :: gq0 (:,:,:) => null()  !< updated tracers

    contains
      procedure :: create  => stateout_create  !<   allocate array data
  end type GFS_stateout_type


!---------------------------------------------------------------------------------------
! GFS_sfcprop_type
!   surface properties that may be read in and/or updated by climatology or observations
!---------------------------------------------------------------------------------------
!! \section arg_table_GFS_sfcprop_type
!! \htmlinclude GFS_sfcprop_type.html
!!
  type GFS_sfcprop_type

!--- In (radiation and physics)
    real (kind=kind_phys), pointer :: slmsk  (:)   => null()  !< sea/land mask array (sea:0,land:1,sea-ice:2)
    real (kind=kind_phys), pointer :: oceanfrac(:) => null()  !< ocean fraction [0:1]
    real (kind=kind_phys), pointer :: landfrac(:)  => null()  !< land  fraction [0:1]

!--- In (lakes)
    real (kind=kind_phys), pointer :: lakefrac(:)  => null()  !< lake  fraction [0:1]
    real (kind=kind_phys), pointer :: lakedepth(:) => null()  !< lake  depth [ m ]
    real (kind=kind_phys), pointer :: clm_lakedepth(:) => null()  !< clm internal lake depth [ m ]
    integer,               pointer :: use_lake_model(:) => null()!1=run lake, 2=run lake&nsst, 0=no lake
    real (kind=kind_phys), pointer :: lake_t2m (:)   => null()  !< 2 meter temperature from CLM Lake model
    real (kind=kind_phys), pointer :: lake_q2m (:)   => null()  !< 2 meter humidity from CLM Lake model

    real (kind=kind_phys), pointer :: h_ML(:)      => null()  !Mixed Layer depth of lakes [m]
    real (kind=kind_phys), pointer :: t_ML(:)      => null()  !Mixing layer temperature in K
    real (kind=kind_phys), pointer :: t_mnw(:)     => null()  !Mean temperature of the water column [K]
    real (kind=kind_phys), pointer :: h_talb(:)    => null()  !the thermally active layer depth of the bottom sediments [m]
    real (kind=kind_phys), pointer :: t_talb(:)    => null()  !Temperature at the bottom of the sediment upper layer [K]
    real (kind=kind_phys), pointer :: t_bot1(:)    => null()  !Temperature at the water-bottom sediment interface [K]
    real (kind=kind_phys), pointer :: t_bot2(:)    => null()  !Temperature for bottom layer of water [K]
    real (kind=kind_phys), pointer :: c_t(:)       => null()  !Shape factor of water temperature vertical profile
    real (kind=kind_phys), pointer :: T_snow(:)    => null()  !temperature of snow on a lake [K]
    real (kind=kind_phys), pointer :: T_ice(:)     => null()  !temperature of ice on a lake [K]

    real (kind=kind_phys), pointer :: tsfc   (:)   => null()  !< surface air temperature in K
    real (kind=kind_phys), pointer :: vegtype_frac (:,:) => null()  !< fractions [0:1] of veg. categories
    real (kind=kind_phys), pointer :: soiltype_frac(:,:) => null()  !< fractions [0:1] of soil categories
                                                              !< [tsea in gbphys.f]
    real (kind=kind_phys), pointer :: tsfco  (:)   => null()  !< sst in K
    real (kind=kind_phys), pointer :: usfco  (:)   => null()  !< surface zonal current in m s-1
    real (kind=kind_phys), pointer :: vsfco  (:)   => null()  !< surface meridional current in m s-1
    real (kind=kind_phys), pointer :: tsfcl  (:)   => null()  !< surface land temperature in K
    real (kind=kind_phys), pointer :: tisfc  (:)   => null()  !< surface temperature over ice fraction
    real (kind=kind_phys), pointer :: tiice(:,:)   => null()  !< internal ice temperature
    real (kind=kind_phys), pointer :: snowd  (:)   => null()  !< snow depth water equivalent in mm ; same as snwdph
    real (kind=kind_phys), pointer :: zorl   (:)   => null()  !< composite surface roughness in cm
    real (kind=kind_phys), pointer :: zorlw  (:)   => null()  !< water surface roughness in cm
    real (kind=kind_phys), pointer :: zorll  (:)   => null()  !< land surface roughness in cm
    real (kind=kind_phys), pointer :: zorli  (:)   => null()  !< ice  surface roughness in cm
    real (kind=kind_phys), pointer :: zorlwav(:)   => null()  !< wave surface roughness in cm derived from wave model
    real (kind=kind_phys), pointer :: fice   (:)   => null()  !< ice fraction over open water grid
    real (kind=kind_phys), pointer :: snodl  (:)   => null()  !< snow depth over land
    real (kind=kind_phys), pointer :: weasdl (:)   => null()  !< weasd over land
    real (kind=kind_phys), pointer :: snodi  (:)   => null()  !< snow depth over ice
    real (kind=kind_phys), pointer :: weasdi (:)   => null()  !< weasd over ice
    real (kind=kind_phys), pointer :: hprime (:,:) => null()  !< orographic metrics
    real (kind=kind_phys), pointer :: dust12m_in  (:,:,:) => null()  !< fengsha dust input
    real (kind=kind_phys), pointer :: emi_in (:,:) => null()  !< anthropogenic background input
    real (kind=kind_phys), pointer :: smoke_RRFS(:,:,:) => null()  !< RRFS fire input hourly
    real (kind=kind_phys), pointer :: smoke2d_RRFS(:,:) => null()  !< RRFS fire input daily
    real (kind=kind_phys), pointer :: z0base (:)   => null()  !< background or baseline surface roughness length in m
    real (kind=kind_phys), pointer :: semisbase(:) => null()  !< background surface emissivity
    real (kind=kind_phys), pointer :: sfalb_lnd (:) => null() !< surface albedo over land for LSM
    real (kind=kind_phys), pointer :: sfalb_ice (:) => null() !< surface albedo over ice for LSM
    real (kind=kind_phys), pointer :: emis_lnd (:)  => null() !< surface emissivity over land for LSM
    real (kind=kind_phys), pointer :: emis_ice (:)  => null() !< surface emissivity over ice for LSM
    real (kind=kind_phys), pointer :: emis_wat (:)  => null() !< surface emissivity over water
    real (kind=kind_phys), pointer :: sfalb_lnd_bck (:) => null() !< snow-free albedo over land

!--- In (radiation only)
    real (kind=kind_phys), pointer :: sncovr (:)   => null()  !< snow cover in fraction over land
    real (kind=kind_phys), pointer :: sncovr_ice (:)  => null()  !< snow cover in fraction over ice (RUC LSM only)
    real (kind=kind_phys), pointer :: snoalb (:)   => null()  !< maximum snow albedo in fraction
    real (kind=kind_phys), pointer :: alvsf  (:)   => null()  !< mean vis albedo with strong cosz dependency
    real (kind=kind_phys), pointer :: alnsf  (:)   => null()  !< mean nir albedo with strong cosz dependency
    real (kind=kind_phys), pointer :: alvwf  (:)   => null()  !< mean vis albedo with weak cosz dependency
    real (kind=kind_phys), pointer :: alnwf  (:)   => null()  !< mean nir albedo with weak cosz dependency
    real (kind=kind_phys), pointer :: facsf  (:)   => null()  !< fractional coverage with strong cosz dependency
    real (kind=kind_phys), pointer :: facwf  (:)   => null()  !< fractional coverage with   weak cosz dependency

!--- In (physics only)
    integer,               pointer :: slope  (:)   => null()  !< sfc slope type for lsm
    integer,               pointer :: slope_save (:) => null()!< sfc slope type save
    real (kind=kind_phys), pointer :: shdmin (:)   => null()  !< min fractional coverage of green veg
    real (kind=kind_phys), pointer :: shdmax (:)   => null()  !< max fractnl cover of green veg (not used)
    real (kind=kind_phys), pointer :: tg3    (:)   => null()  !< deep soil temperature
    real (kind=kind_phys), pointer :: vfrac  (:)   => null()  !< vegetation fraction
    integer,               pointer :: vtype  (:)   => null()  !< vegetation type
    integer,               pointer :: stype  (:)   => null()  !< soil type
    integer,               pointer :: scolor  (:)   => null()  !< soil color
    integer,               pointer :: vtype_save (:) => null()!< vegetation type save
    integer,               pointer :: stype_save (:) => null()!< soil type save
    integer,               pointer :: scolor_save (:) => null()!< soil color save
    real (kind=kind_phys), pointer :: uustar (:)   => null()  !< boundary layer parameter
    real (kind=kind_phys), pointer :: oro    (:)   => null()  !< orography
    real (kind=kind_phys), pointer :: oro_uf (:)   => null()  !< unfiltered orography
    real (kind=kind_phys), pointer :: evap   (:)   => null()  !<
    real (kind=kind_phys), pointer :: hflx   (:)   => null()  !<
    real (kind=kind_phys), pointer :: qss    (:)   => null()  !<

!--- fire_behavior
    real (kind=kind_phys), pointer :: hflx_fire (:)   => null()  !< kinematic surface upward sensible heat flux of fire
    real (kind=kind_phys), pointer :: evap_fire (:)   => null()  !< kinematic surface upward latent heat flux of fire
    real (kind=kind_phys), pointer :: smoke_fire (:)   => null()  !< smoke emission of fire

!-- In/Out
    real (kind=kind_phys), pointer :: maxupmf(:)   => null()  !< maximum up draft mass flux for Grell-Freitas
    real (kind=kind_phys), pointer :: conv_act(:)  => null()  !< convective activity counter for Grell-Freitas
    real (kind=kind_phys), pointer :: conv_act_m(:)=> null()  !< midlevel convective activity counter for Grell-Freitas
    real (kind=kind_phys), pointer :: hice   (:)   => null()  !< sea ice thickness
    real (kind=kind_phys), pointer :: weasd  (:)   => null()  !< water equiv of accumulated snow depth (kg/m**2)
                                                              !< over land and sea ice
    real (kind=kind_phys), pointer :: canopy (:)   => null()  !< canopy water
    real (kind=kind_phys), pointer :: ffmm   (:)   => null()  !< fm parameter from PBL scheme
    real (kind=kind_phys), pointer :: ffhh   (:)   => null()  !< fh parameter from PBL scheme
    real (kind=kind_phys), pointer :: f10m   (:)   => null()  !< fm at 10m - Ratio of sigma level 1 wind and 10m wind
    real (kind=kind_phys), pointer :: rca    (:)   => null()  !< canopy resistance
    real (kind=kind_phys), pointer :: tprcp  (:)   => null()  !< sfc_fld%tprcp - total precipitation
    real (kind=kind_phys), pointer :: srflag (:)   => null()  !< sfc_fld%srflag - snow/rain flag for precipitation
    real (kind=kind_phys), pointer :: slc    (:,:) => null()  !< liquid soil moisture
    real (kind=kind_phys), pointer :: smc    (:,:) => null()  !< total soil moisture
    real (kind=kind_phys), pointer :: stc    (:,:) => null()  !< soil temperature

!--- Out
    real (kind=kind_phys), pointer :: t2m    (:)   => null()  !< 2 meter temperature
    real (kind=kind_phys), pointer :: th2m   (:)   => null()  !< 2 meter potential temperature
    real (kind=kind_phys), pointer :: q2m    (:)   => null()  !< 2 meter humidity

! -- In/Out for Noah MP
    real (kind=kind_phys), pointer :: snowxy  (:)  => null()  !<
    real (kind=kind_phys), pointer :: tvxy    (:)  => null()  !< veg temp
    real (kind=kind_phys), pointer :: tgxy    (:)  => null()  !< ground temp
    real (kind=kind_phys), pointer :: canicexy(:)  => null()  !<
    real (kind=kind_phys), pointer :: canliqxy(:)  => null()  !<
    real (kind=kind_phys), pointer :: eahxy   (:)  => null()  !<
    real (kind=kind_phys), pointer :: tahxy   (:)  => null()  !<
    real (kind=kind_phys), pointer :: cmxy    (:)  => null()  !<
    real (kind=kind_phys), pointer :: chxy    (:)  => null()  !<
    real (kind=kind_phys), pointer :: fwetxy  (:)  => null()  !<
    real (kind=kind_phys), pointer :: sneqvoxy(:)  => null()  !<
    real (kind=kind_phys), pointer :: alboldxy(:)  => null()  !<
    real (kind=kind_phys), pointer :: qsnowxy (:)  => null()  !<
    real (kind=kind_phys), pointer :: wslakexy(:)  => null()  !<
    real (kind=kind_phys), pointer :: zwtxy   (:)  => null()  !<
    real (kind=kind_phys), pointer :: waxy    (:)  => null()  !<
    real (kind=kind_phys), pointer :: wtxy    (:)  => null()  !<
    real (kind=kind_phys), pointer :: lfmassxy(:)  => null()  !<
    real (kind=kind_phys), pointer :: rtmassxy(:)  => null()  !<
    real (kind=kind_phys), pointer :: stmassxy(:)  => null()  !<
    real (kind=kind_phys), pointer :: woodxy  (:)  => null()  !<
    real (kind=kind_phys), pointer :: stblcpxy(:)  => null()  !<
    real (kind=kind_phys), pointer :: fastcpxy(:)  => null()  !<
    real (kind=kind_phys), pointer :: xsaixy  (:)  => null()  !<
    real (kind=kind_phys), pointer :: xlaixy  (:)  => null()  !<
    real (kind=kind_phys), pointer :: taussxy (:)  => null()  !<
    real (kind=kind_phys), pointer :: smcwtdxy(:)  => null()  !<
    real (kind=kind_phys), pointer :: deeprechxy(:)=> null()  !<
    real (kind=kind_phys), pointer :: rechxy  (:)  => null()  !<
    real (kind=kind_phys), pointer :: albdirvis_lnd (:) => null()  !<
    real (kind=kind_phys), pointer :: albdirnir_lnd (:) => null()  !<
    real (kind=kind_phys), pointer :: albdifvis_lnd (:) => null()  !<
    real (kind=kind_phys), pointer :: albdifnir_lnd (:) => null()  !<

    real (kind=kind_phys), pointer :: albdirvis_ice (:) => null()  !<
    real (kind=kind_phys), pointer :: albdifvis_ice (:) => null()  !<
    real (kind=kind_phys), pointer :: albdirnir_ice (:) => null()  !<
    real (kind=kind_phys), pointer :: albdifnir_ice (:) => null()  !<
!   real (kind=kind_phys), pointer :: sfalb_ice     (:) => null()  !<

    real (kind=kind_phys), pointer :: snicexy   (:,:) => null()  !<
    real (kind=kind_phys), pointer :: snliqxy   (:,:) => null()  !<
    real (kind=kind_phys), pointer :: tsnoxy    (:,:) => null()  !<
    real (kind=kind_phys), pointer :: smoiseq   (:,:) => null()  !<
    real (kind=kind_phys), pointer :: zsnsoxy   (:,:) => null()  !<

!--- NSSTM variables  (only allocated when [Model%nstf_name(1) > 0])
    real (kind=kind_phys), pointer :: tref   (:)   => null()  !< nst_fld%Tref - Reference Temperature
    real (kind=kind_phys), pointer :: z_c    (:)   => null()  !< nst_fld%z_c - Sub layer cooling thickness
    real (kind=kind_phys), pointer :: c_0    (:)   => null()  !< nst_fld%c_0 - coefficient1 to calculate d(Tz)/d(Ts)
    real (kind=kind_phys), pointer :: c_d    (:)   => null()  !< nst_fld%c_d - coefficient2 to calculate d(Tz)/d(Ts)
    real (kind=kind_phys), pointer :: w_0    (:)   => null()  !< nst_fld%w_0 - coefficient3 to calculate d(Tz)/d(Ts)
    real (kind=kind_phys), pointer :: w_d    (:)   => null()  !< nst_fld%w_d - coefficient4 to calculate d(Tz)/d(Ts)
    real (kind=kind_phys), pointer :: xt     (:)   => null()  !< nst_fld%xt      heat content in DTL
    real (kind=kind_phys), pointer :: xs     (:)   => null()  !< nst_fld%xs      salinity  content in DTL
    real (kind=kind_phys), pointer :: xu     (:)   => null()  !< nst_fld%xu      u current content in DTL
    real (kind=kind_phys), pointer :: xv     (:)   => null()  !< nst_fld%xv      v current content in DTL
    real (kind=kind_phys), pointer :: xz     (:)   => null()  !< nst_fld%xz      DTL thickness
    real (kind=kind_phys), pointer :: zm     (:)   => null()  !< nst_fld%zm      MXL thickness
    real (kind=kind_phys), pointer :: xtts   (:)   => null()  !< nst_fld%xtts    d(xt)/d(ts)
    real (kind=kind_phys), pointer :: xzts   (:)   => null()  !< nst_fld%xzts    d(xz)/d(ts)
    real (kind=kind_phys), pointer :: d_conv (:)   => null()  !< nst_fld%d_conv  thickness of Free Convection Layer (FCL)
    real (kind=kind_phys), pointer :: ifd    (:)   => null()  !< nst_fld%ifd     index to start DTM run or not
    real (kind=kind_phys), pointer :: dt_cool(:)   => null()  !< nst_fld%dt_cool Sub layer cooling amount
    real (kind=kind_phys), pointer :: qrain  (:)   => null()  !< nst_fld%qrain   sensible heat flux due to rainfall (watts)

    ! Soil properties for RUC LSM (number of levels different from NOAH 4-layer model)
    real (kind=kind_phys), pointer :: wetness(:)         => null()  !< normalized soil wetness for lsm
    real (kind=kind_phys), pointer :: sh2o(:,:)          => null()  !< volume fraction of unfrozen soil moisture for lsm
    real (kind=kind_phys), pointer :: keepsmfr(:,:)      => null()  !< RUC LSM: frozen moisture in soil
    real (kind=kind_phys), pointer :: smois(:,:)         => null()  !< volumetric fraction of soil moisture for lsm
    real (kind=kind_phys), pointer :: tslb(:,:)          => null()  !< soil temperature for land surface model
    real (kind=kind_phys), pointer :: flag_frsoil(:,:)   => null()  !< RUC LSM: flag for frozen soil physics
    !
    real (kind=kind_phys), pointer :: clw_surf_land(:)   => null()  !< RUC LSM: moist cloud water mixing ratio at surface over land
    real (kind=kind_phys), pointer :: clw_surf_ice(:)    => null()  !< RUC LSM: moist cloud water mixing ratio at surface over ice
    real (kind=kind_phys), pointer :: qwv_surf_land(:)   => null()  !< RUC LSM: water vapor mixing ratio at surface over land
    real (kind=kind_phys), pointer :: qwv_surf_ice(:)    => null()  !< RUC LSM: water vapor mixing ratio at surface over ice
    real (kind=kind_phys), pointer :: rhofr(:)           => null()  !< RUC LSM: internal density of frozen precipitation
    real (kind=kind_phys), pointer :: tsnow_land(:)      => null()  !< RUC LSM: snow temperature at the bottom of the first snow layer over land
    real (kind=kind_phys), pointer :: tsnow_ice(:)       => null()  !< RUC LSM: snow temperature at the bottom of the first snow layer over ice
    real (kind=kind_phys), pointer :: snowfallac_land(:) => null()  !< ruc lsm diagnostics over land
    real (kind=kind_phys), pointer :: snowfallac_ice(:)  => null()  !< ruc lsm diagnostics over ice
    real (kind=kind_phys), pointer :: acsnow_land(:)     => null()  !< ruc lsm diagnostics over land
    real (kind=kind_phys), pointer :: acsnow_ice(:)      => null()  !< ruc lsm diagnostics over ice

    !  MYNN surface layer
    real (kind=kind_phys), pointer :: ustm (:)         => null()  !u* including drag
    real (kind=kind_phys), pointer :: zol(:)           => null()  !surface stability parameter
    real (kind=kind_phys), pointer :: mol(:)           => null()  !theta star
    real (kind=kind_phys), pointer :: rmol(:)          => null()  !reciprocal of obukhov length
    real (kind=kind_phys), pointer :: flhc(:)          => null()  !drag coeff for heat
    real (kind=kind_phys), pointer :: flqc(:)          => null()  !drag coeff for moisture
    real (kind=kind_phys), pointer :: chs2(:)          => null()  !exch coeff for heat at 2m
    real (kind=kind_phys), pointer :: cqs2(:)          => null()  !exch coeff for moisture at 2m
    real (kind=kind_phys), pointer :: lh(:)            => null()  !latent heating at the surface

    !---- precipitation amounts from previous time step for RUC LSM/NoahMP LSM
    real (kind=kind_phys), pointer :: raincprv  (:)    => null()  !< explicit rainfall from previous timestep
    real (kind=kind_phys), pointer :: rainncprv (:)    => null()  !< convective_precipitation_amount from previous timestep
    real (kind=kind_phys), pointer :: iceprv    (:)    => null()  !< ice amount from previous timestep
    real (kind=kind_phys), pointer :: snowprv   (:)    => null()  !< snow amount from previous timestep
    real (kind=kind_phys), pointer :: graupelprv(:)    => null()  !< graupel amount from previous timestep

    !---- precipitation rates from previous time step for NoahMP LSM
    real (kind=kind_phys), pointer :: draincprv  (:)    => null()  !< convective precipitation rate from previous timestep
    real (kind=kind_phys), pointer :: drainncprv (:)    => null()  !< explicit rainfall rate from previous timestep
    real (kind=kind_phys), pointer :: diceprv    (:)    => null()  !< ice precipitation rate from previous timestep
    real (kind=kind_phys), pointer :: dsnowprv   (:)    => null()  !< snow precipitation rate from previous timestep
    real (kind=kind_phys), pointer :: dgraupelprv(:)    => null()  !< graupel precipitation rate from previous timestep

    ! CLM Lake model internal variables:
    real (kind=kind_phys), pointer :: lake_albedo(:)     => null()  !
    real (kind=kind_phys), pointer :: input_lakedepth(:) => null()  !
    real (kind=kind_phys), pointer :: lake_h2osno2d(:)   => null()  !
    real (kind=kind_phys), pointer :: lake_sndpth2d(:)   => null()  !
    real (kind=kind_phys), pointer :: lake_snl2d(:)      => null()  !
    real (kind=kind_phys), pointer :: lake_snow_z3d(:,:)      => null()  !
    real (kind=kind_phys), pointer :: lake_snow_dz3d(:,:)     => null()  !
    real (kind=kind_phys), pointer :: lake_snow_zi3d(:,:)     => null()  !
    real (kind=kind_phys), pointer :: lake_h2osoi_vol3d(:,:)  => null()  !
    real (kind=kind_phys), pointer :: lake_h2osoi_liq3d(:,:)   => null()  !
    real (kind=kind_phys), pointer :: lake_h2osoi_ice3d(:,:)   => null()  !
    real (kind=kind_phys), pointer :: lake_tsfc(:)   => null()  !
    real (kind=kind_phys), pointer :: lake_t_soisno3d(:,:) => null()  !
    real (kind=kind_phys), pointer :: lake_t_lake3d(:,:) => null()  !
    real (kind=kind_phys), pointer :: lake_savedtke12d(:)=> null()  !
    real (kind=kind_phys), pointer :: lake_icefrac3d(:,:)=> null()
    real (kind=kind_phys), pointer :: lake_rho0(:)=> null()
    real (kind=kind_phys), pointer :: lake_ht(:)=> null()
    integer, pointer :: lake_is_salty(:) => null()
    integer, pointer :: lake_cannot_freeze(:) => null()
    real (kind=kind_phys), pointer :: clm_lake_initialized(:) => null() !< lakeini was called
    !--- aerosol surface emissions for Thompson microphysics & smoke dust
    real (kind=kind_phys), pointer :: emdust  (:)     => null()  !< instantaneous dust emission
    real (kind=kind_phys), pointer :: emseas  (:)     => null()  !< instantaneous sea salt emission
    real (kind=kind_phys), pointer :: emanoc  (:)     => null()  !< instantaneous anthro. oc emission

    !--- Smoke. These 2 arrays are input smoke emission and frp
    real (kind=kind_phys), pointer :: ebb_smoke_in(:)   => null()  !< input smoke emission
    real (kind=kind_phys), pointer :: frp_output  (:)   => null()  !< output FRP
    !--- For fire diurnal cycle
    real (kind=kind_phys), pointer :: fhist       (:)   => null()  !< instantaneous fire coef_bb
    real (kind=kind_phys), pointer :: coef_bb_dc  (:)   => null()  !< instantaneous fire coef_bb
    integer, pointer               :: fire_type   (:)   => null()  !< fire type
    real (kind=kind_phys), pointer :: peak_hr     (:)   => null()  !< peak hour of fire emissions
    real (kind=kind_phys), pointer :: lu_nofire   (:)   => null()  !<lu_nofire pixels
    real (kind=kind_phys), pointer :: lu_qfire    (:)   => null()  !<lu_qfire pixels
    !--- wildfire heat flux
    real (kind=kind_phys), pointer :: fire_heat_flux   (:) => null() !< heat flux from wildfire
    real (kind=kind_phys), pointer :: frac_grid_burned (:) => null() !< fraction of grid cell burning

    !--- For smoke and dust auxiliary inputs
    real (kind=kind_phys), pointer :: fire_in   (:,:)   => null()  !< fire auxiliary inputs

    contains
      procedure :: create  => sfcprop_create  !<   allocate array data
  end type GFS_sfcprop_type


!---------------------------------------------------------------------
! GFS_coupling_type
!   fields to/from other coupled components (e.g. land/ice/ocean/etc.)
!---------------------------------------------------------------------
!! \section arg_table_GFS_coupling_type
!! \htmlinclude GFS_coupling_type.html
!!
  type GFS_coupling_type

!--- Out (radiation only)
    real (kind=kind_phys), pointer :: nirbmdi(:)     => null()   !< sfc nir beam sw downward flux (w/m2)
    real (kind=kind_phys), pointer :: nirdfdi(:)     => null()   !< sfc nir diff sw downward flux (w/m2)
    real (kind=kind_phys), pointer :: visbmdi(:)     => null()   !< sfc uv+vis beam sw downward flux (w/m2)
    real (kind=kind_phys), pointer :: visdfdi(:)     => null()   !< sfc uv+vis diff sw downward flux (w/m2)
    real (kind=kind_phys), pointer :: nirbmui(:)     => null()   !< sfc nir beam sw upward flux (w/m2)
    real (kind=kind_phys), pointer :: nirdfui(:)     => null()   !< sfc nir diff sw upward flux (w/m2)
    real (kind=kind_phys), pointer :: visbmui(:)     => null()   !< sfc uv+vis beam sw upward flux (w/m2)
    real (kind=kind_phys), pointer :: visdfui(:)     => null()   !< sfc uv+vis diff sw upward flux (w/m2)

    ! RRTMGP
    real (kind=kind_phys), pointer :: fluxlwUP_jac(:,:)       => null()  !< RRTMGP Jacobian of upward longwave all-sky flux
    real (kind=kind_phys), pointer :: htrlw(:,:)              => null()  !< RRTMGP updated LW heating rate
    real (kind=kind_phys), pointer :: tsfc_radtime(:)         => null()  !< RRTMGP surface temperature on radiation timestep
    real (kind=kind_phys), pointer :: fluxlwUP_radtime(:,:)   => null()  !< RRTMGP upward   longwave  all-sky flux profile
    real (kind=kind_phys), pointer :: fluxlwDOWN_radtime(:,:) => null()  !< RRTMGP downward  longwave  all-sky flux profile

    !--- In (physics only)
    real (kind=kind_phys), pointer :: sfcdsw(:)      => null()   !< total sky sfc downward sw flux ( w/m**2 )
                                                                 !< GFS_radtend_type%sfcfsw%dnfxc
    real (kind=kind_phys), pointer :: sfcnsw(:)      => null()   !< total sky sfc netsw flx into ground(w/m**2)
                                                                 !< difference of dnfxc & upfxc from GFS_radtend_type%sfcfsw
    real (kind=kind_phys), pointer :: sfcdlw(:)      => null()   !< total sky sfc downward lw flux ( w/m**2 )
                                                                 !< GFS_radtend_type%sfclsw%dnfxc
    real (kind=kind_phys), pointer :: sfculw(:)      => null()   !< total sky sfc upward lw flux ( w/m**2 )

!--- incoming quantities
    real (kind=kind_phys), pointer :: dusfcin_cpl(:)          => null()   !< aoi_fld%dusfcin(item,lan)
    real (kind=kind_phys), pointer :: dvsfcin_cpl(:)          => null()   !< aoi_fld%dvsfcin(item,lan)
    real (kind=kind_phys), pointer :: dtsfcin_cpl(:)          => null()   !< aoi_fld%dtsfcin(item,lan)
    real (kind=kind_phys), pointer :: dqsfcin_cpl(:)          => null()   !< aoi_fld%dqsfcin(item,lan)
    real (kind=kind_phys), pointer :: ulwsfcin_cpl(:)         => null()   !< aoi_fld%ulwsfcin(item,lan)
!   real (kind=kind_phys), pointer :: tseain_cpl(:)           => null()   !< aoi_fld%tseain(item,lan)
!   real (kind=kind_phys), pointer :: tisfcin_cpl(:)          => null()   !< aoi_fld%tisfcin(item,lan)
!   real (kind=kind_phys), pointer :: ficein_cpl(:)           => null()   !< aoi_fld%ficein(item,lan)
!   real (kind=kind_phys), pointer :: hicein_cpl(:)           => null()   !< aoi_fld%hicein(item,lan)
    real (kind=kind_phys), pointer :: hsnoin_cpl(:)           => null()   !< aoi_fld%hsnoin(item,lan)
!   real (kind=kind_phys), pointer :: sfc_alb_nir_dir_cpl(:)  => null()   !< sfc nir albedo for direct rad
!   real (kind=kind_phys), pointer :: sfc_alb_nir_dif_cpl(:)  => null()   !< sfc nir albedo for diffuse rad
!   real (kind=kind_phys), pointer :: sfc_alb_vis_dir_cpl(:)  => null()   !< sfc vis albedo for direct rad
!   real (kind=kind_phys), pointer :: sfc_alb_vis_dif_cpl(:)  => null()   !< sfc vis albedo for diffuse rad
    !--- only variable needed for cplwav2atm=.TRUE.
!   real (kind=kind_phys), pointer :: zorlwav_cpl(:)          => null()   !< roughness length from wave model
    !--- also needed for ice/ocn coupling
    real (kind=kind_phys), pointer :: slimskin_cpl(:)=> null()   !< aoi_fld%slimskin(item,lan)
    !--- variables needed for use_med_flux =.TRUE.
    real (kind=kind_phys), pointer :: dusfcin_med(:)         => null()   !< sfc u momentum flux over ocean
    real (kind=kind_phys), pointer :: dvsfcin_med(:)         => null()   !< sfc v momentum flux over ocean
    real (kind=kind_phys), pointer :: dtsfcin_med(:)         => null()   !< sfc latent heat flux over ocean
    real (kind=kind_phys), pointer :: dqsfcin_med(:)         => null()   !< sfc sensible heat flux over ocean
    real (kind=kind_phys), pointer :: ulwsfcin_med(:)        => null()   !< sfc upward lw flux over ocean
    !--- variables needed for cpllnd = .TRUE. and cpllnd2atm=.TRUE.
    real (kind=kind_phys), pointer :: sncovr1_lnd(:)         => null()   !< sfc snow area fraction over land
    real (kind=kind_phys), pointer :: qsurf_lnd(:)           => null()   !< sfc specific humidity
    real (kind=kind_phys), pointer :: evap_lnd(:)            => null()   !< sfc latent heat flux over land, converted to evaporative flux
    real (kind=kind_phys), pointer :: hflx_lnd(:)            => null()   !< sfc sensible heat flux over land
    real (kind=kind_phys), pointer :: ep_lnd(:)              => null()   !< sfc up pot latent heat flux over land
    real (kind=kind_phys), pointer :: t2mmp_lnd(:)           => null()   !< 2 meter temperature over land 
    real (kind=kind_phys), pointer :: q2mp_lnd(:)            => null()   !< 2 meter spec humidity over land
    real (kind=kind_phys), pointer :: gflux_lnd(:)           => null()   !< soil heat flux over land
    real (kind=kind_phys), pointer :: runoff_lnd(:)          => null()   !< surface runoff over land
    real (kind=kind_phys), pointer :: drain_lnd(:)           => null()   !< subsurface runoff over land
    real (kind=kind_phys), pointer :: cmm_lnd(:)             => null()   !< surface drag wind speed for momentum
    real (kind=kind_phys), pointer :: chh_lnd(:)             => null()   !< surface drag mass flux for heat and moisture 
    real (kind=kind_phys), pointer :: zvfun_lnd(:)           => null()   !< function of surface roughness length and green vegetation fraction 

!--- outgoing accumulated quantities
    real (kind=kind_phys), pointer :: rain_cpl  (:)  => null()   !< total rain precipitation
    real (kind=kind_phys), pointer :: rainc_cpl (:)  => null()   !< convective rain precipitation
    real (kind=kind_phys), pointer :: snow_cpl  (:)  => null()   !< total snow precipitation
    real (kind=kind_phys), pointer :: dusfc_cpl (:)  => null()   !< sfc u momentum flux
    real (kind=kind_phys), pointer :: dvsfc_cpl (:)  => null()   !< sfc v momentum flux
    real (kind=kind_phys), pointer :: dtsfc_cpl (:)  => null()   !< sfc sensible heat flux
    real (kind=kind_phys), pointer :: dqsfc_cpl (:)  => null()   !< sfc   latent heat flux
    real (kind=kind_phys), pointer :: dlwsfc_cpl(:)  => null()   !< sfc downward lw flux (w/m**2)
    real (kind=kind_phys), pointer :: dswsfc_cpl(:)  => null()   !< sfc downward sw flux (w/m**2)
    real (kind=kind_phys), pointer :: dnirbm_cpl(:)  => null()   !< sfc nir beam downward sw flux (w/m**2)
    real (kind=kind_phys), pointer :: dnirdf_cpl(:)  => null()   !< sfc nir diff downward sw flux (w/m**2)
    real (kind=kind_phys), pointer :: dvisbm_cpl(:)  => null()   !< sfc uv+vis beam dnwd sw flux (w/m**2)
    real (kind=kind_phys), pointer :: dvisdf_cpl(:)  => null()   !< sfc uv+vis diff dnwd sw flux (w/m**2)
    real (kind=kind_phys), pointer :: nlwsfc_cpl(:)  => null()   !< net downward lw flux (w/m**2)
    real (kind=kind_phys), pointer :: nswsfc_cpl(:)  => null()   !< net downward sw flux (w/m**2)
    real (kind=kind_phys), pointer :: nnirbm_cpl(:)  => null()   !< net nir beam downward sw flux (w/m**2)
    real (kind=kind_phys), pointer :: nnirdf_cpl(:)  => null()   !< net nir diff downward sw flux (w/m**2)
    real (kind=kind_phys), pointer :: nvisbm_cpl(:)  => null()   !< net uv+vis beam downward sw rad flux (w/m**2)
    real (kind=kind_phys), pointer :: nvisdf_cpl(:)  => null()   !< net uv+vis diff downward sw rad flux (w/m**2)

!--- outgoing instantaneous quantities
    real (kind=kind_phys), pointer :: dusfci_cpl (:) => null()   !< instantaneous sfc u momentum flux
    real (kind=kind_phys), pointer :: dvsfci_cpl (:) => null()   !< instantaneous sfc v momentum flux
    real (kind=kind_phys), pointer :: dtsfci_cpl (:) => null()   !< instantaneous sfc sensible heat flux
    real (kind=kind_phys), pointer :: dqsfci_cpl (:) => null()   !< instantaneous sfc   latent heat flux
    real (kind=kind_phys), pointer :: dlwsfci_cpl(:) => null()   !< instantaneous sfc downward lw flux
    real (kind=kind_phys), pointer :: dswsfci_cpl(:) => null()   !< instantaneous sfc downward sw flux
    real (kind=kind_phys), pointer :: dnirbmi_cpl(:) => null()   !< instantaneous sfc nir beam downward sw flux
    real (kind=kind_phys), pointer :: dnirdfi_cpl(:) => null()   !< instantaneous sfc nir diff downward sw flux
    real (kind=kind_phys), pointer :: dvisbmi_cpl(:) => null()   !< instantaneous sfc uv+vis beam downward sw flux
    real (kind=kind_phys), pointer :: dvisdfi_cpl(:) => null()   !< instantaneous sfc uv+vis diff downward sw flux
    real (kind=kind_phys), pointer :: nlwsfci_cpl(:) => null()   !< instantaneous net sfc downward lw flux
    real (kind=kind_phys), pointer :: nswsfci_cpl(:) => null()   !< instantaneous net sfc downward sw flux
    real (kind=kind_phys), pointer :: nnirbmi_cpl(:) => null()   !< instantaneous net nir beam sfc downward sw flux
    real (kind=kind_phys), pointer :: nnirdfi_cpl(:) => null()   !< instantaneous net nir diff sfc downward sw flux
    real (kind=kind_phys), pointer :: nvisbmi_cpl(:) => null()   !< instantaneous net uv+vis beam downward sw flux
    real (kind=kind_phys), pointer :: nvisdfi_cpl(:) => null()   !< instantaneous net uv+vis diff downward sw flux
    real (kind=kind_phys), pointer :: t2mi_cpl   (:) => null()   !< instantaneous T2m
    real (kind=kind_phys), pointer :: q2mi_cpl   (:) => null()   !< instantaneous Q2m
    real (kind=kind_phys), pointer :: u10mi_cpl  (:) => null()   !< instantaneous U10m
    real (kind=kind_phys), pointer :: v10mi_cpl  (:) => null()   !< instantaneous V10m
    real (kind=kind_phys), pointer :: tsfci_cpl  (:) => null()   !< instantaneous sfc temperature
    real (kind=kind_phys), pointer :: psurfi_cpl (:) => null()   !< instantaneous sfc pressure

    !--- topography-based information for the coupling system
    real (kind=kind_phys), pointer :: oro_cpl    (:) => null()   !< orography          (  oro from GFS_sfcprop_type)
    real (kind=kind_phys), pointer :: slmsk_cpl  (:) => null()   !< Land/Sea/Ice mask  (slmsk from GFS_sfcprop_type)

    !--- cellular automata
    real (kind=kind_phys), pointer :: ca1      (:)   => null() !
    real (kind=kind_phys), pointer :: ca2      (:)   => null() !
    real (kind=kind_phys), pointer :: ca3      (:)   => null() !
    real (kind=kind_phys), pointer :: ca_deep  (:)   => null() !
    real (kind=kind_phys), pointer :: ca_turb  (:)   => null() !
    real (kind=kind_phys), pointer :: ca_shal  (:)   => null() !
    real (kind=kind_phys), pointer :: ca_rad   (:)   => null() !
    real (kind=kind_phys), pointer :: ca_micro (:)   => null() !
    real (kind=kind_phys), pointer :: condition(:)   => null() !
    !--- stochastic physics
    real (kind=kind_phys), pointer :: shum_wts  (:,:) => null()  !
    real (kind=kind_phys), pointer :: sppt_wts  (:,:) => null()  !
    real (kind=kind_phys), pointer :: skebu_wts (:,:) => null()  !
    real (kind=kind_phys), pointer :: skebv_wts (:,:) => null()  !
    real (kind=kind_phys), pointer :: sfc_wts   (:,:) => null()  ! mg, sfc-perts
    real (kind=kind_phys), pointer :: spp_wts_pbl   (:,:) => null()  ! spp-pbl-perts
    real (kind=kind_phys), pointer :: spp_wts_sfc   (:,:) => null()  ! spp-sfc-perts
    real (kind=kind_phys), pointer :: spp_wts_mp    (:,:) => null()  ! spp-mp-perts
    real (kind=kind_phys), pointer :: spp_wts_gwd   (:,:) => null()  ! spp-gwd-perts
    real (kind=kind_phys), pointer :: spp_wts_rad   (:,:) => null()  ! spp-rad-perts
    real (kind=kind_phys), pointer :: spp_wts_cu_deep (:,:) => null()  ! spp-cu-deep-perts

    !--- aerosol surface emissions for Thompson microphysics
    real (kind=kind_phys), pointer :: nwfa2d  (:)     => null()  !< instantaneous water-friendly sfc aerosol source
    real (kind=kind_phys), pointer :: nifa2d  (:)     => null()  !< instantaneous ice-friendly sfc aerosol source

    !--- For fire diurnal cycle
    real (kind=kind_phys), pointer :: ebu_smoke (:,:) => null()  !< 3D ebu array

    !--- For MYNN PBL transport of  smoke and dust
    real (kind=kind_phys), pointer :: chem3d  (:,:,:)   => null()  !< 3D aod array
    real (kind=kind_phys), pointer :: ddvel   (:,:  )   => null()  !< 2D dry deposition velocity
    !--- For convective wet removal of  smoke and dust
    real (kind=kind_phys), pointer :: wetdpc_flux (:,:) => null()  !< 2D wet deposition array
    !--- For large-scale wet removal of  smoke and dust
    real (kind=kind_phys), pointer :: wetdpr_flux (:,:) => null()  !< 2D wet deposition array
    !--- For dry deposition of smoke and dust
    real (kind=kind_phys), pointer :: drydep_flux (:,:) => null()  !< 2D dry deposition flux of smoke

    !--- Fire plume rise diagnostics
    real (kind=kind_phys), pointer :: min_fplume (:)  => null()  !< minimum plume rise level
    real (kind=kind_phys), pointer :: max_fplume (:)  => null()  !< maximum plume rise level
    real (kind=kind_phys), pointer :: uspdavg (:)     => null()  !< BL average wind speed
    real (kind=kind_phys), pointer :: hpbl_thetav (:) => null()  !< BL depth parcel method
    !--- hourly fire potential index
    real (kind=kind_phys), pointer :: rrfs_hwp   (:)  => null()  !< hourly fire potential index
    real (kind=kind_phys), pointer :: rrfs_hwp_ave   (:)   => null()  !< *Average* hourly fire potential index

    !--- instantaneous quantities for chemistry coupling
    real (kind=kind_phys), pointer :: ushfsfci(:)     => null()  !< instantaneous upward sensible heat flux (w/m**2)
    real (kind=kind_phys), pointer :: qci_conv(:,:)   => null()  !< convective cloud condesate after rainout
    real (kind=kind_phys), pointer :: pfi_lsan(:,:)   => null()  !< instantaneous 3D flux of ice    nonconvective precipitation (kg m-2 s-1)
    real (kind=kind_phys), pointer :: pfl_lsan(:,:)   => null()  !< instantaneous 3D flux of liquid nonconvective precipitation (kg m-2 s-1)

    !-- prognostic updraft area fraction coupling in convection
    real (kind=kind_phys), pointer :: dqdt_qmicro(:,:) => null()  !< instantanious microphysics tendency to be passed from MP to convection

    contains
      procedure :: create  => coupling_create  !<   allocate array data
  end type GFS_coupling_type

!----------------------------------------------------------------
! dtend_var_label
!  Information about first dimension of dtidx
!----------------------------------------------------------------
  type dtend_var_label
    character(len=20) :: name
    character(len=44) :: desc
    character(len=32) :: unit
  end type dtend_var_label

!----------------------------------------------------------------
! dtend_process_label
!  Information about second dimension of dtidx
!----------------------------------------------------------------
  type dtend_process_label
    character(len=20) :: name
    character(len=44) :: desc
    logical :: time_avg
    character(len=20) :: mod_name
  end type dtend_process_label

!----------------------------------------------------------------------------------
! GFS_control_type
!   model control parameters input from a namelist and/or derived from others
!   list of those that can be modified during the run are at the bottom of the list
!----------------------------------------------------------------------------------
!! \section arg_table_GFS_control_type
!! \htmlinclude GFS_control_type.html
!!
  type GFS_control_type

    integer              :: me              !< MPI rank designator
    integer              :: master          !< MPI rank of master atmosphere processor
    type(MPI_Comm)       :: communicator    !< MPI communicator
    integer              :: ntasks          !< MPI size in communicator
    integer              :: nthreads        !< OpenMP threads available for physics
    integer              :: nlunit          !< unit for namelist
    character(len=64)    :: fn_nml          !< namelist filename for surface data cycling
    character(len=:), pointer, dimension(:) :: input_nml_file => null() !< character string containing full namelist
                                                                        !< for use with internal file reads
    integer              :: input_nml_file_length    !< length (number of lines) in namelist for internal reads
    integer              :: logunit
    real(kind=kind_phys) :: fhzero          !< hours between clearing of diagnostic buckets (current bucket)
    real(kind=kind_phys) :: fhzero_array(2) !< array to hold the the hours between clearing of diagnostic buckets
    real(kind=kind_phys) :: fhzero_fhour(2) !< the maximum forecast length for the hours between clearing of diagnostic buckets
    logical              :: ldiag3d         !< flag for 3d diagnostic fields
    logical              :: qdiag3d         !< flag for 3d tracer diagnostic fields
    logical              :: flag_for_gwd_generic_tend  !< true if GFS_GWD_generic should calculate tendencies
    logical              :: flag_for_pbl_generic_tend  !< true if GFS_PBL_generic should calculate tendencies
    logical              :: flag_for_scnv_generic_tend !< true if GFS_DCNV_generic should calculate tendencies
    logical              :: flag_for_dcnv_generic_tend !< true if GFS_DCNV_generic should calculate tendencies
    logical              :: lssav           !< logical flag for storing diagnostics
    integer              :: naux2d          !< number of auxiliary 2d arrays to output (for debugging)
    integer              :: naux3d          !< number of auxiliary 3d arrays to output (for debugging)
    logical, pointer     :: aux2d_time_avg(:) !< flags for time averaging of auxiliary 2d arrays
    logical, pointer     :: aux3d_time_avg(:) !< flags for time averaging of auxiliary 3d arrays

    real(kind=kind_phys) :: fhcyc           !< frequency for surface data cycling (hours)
    integer              :: thermodyn_id    !< valid for GFS only for get_prs/phi
    integer              :: sfcpress_id     !< valid for GFS only for get_prs/phi
    logical              :: gen_coord_hybrid!< for Henry's gen coord

!--- set some grid extent parameters
    integer              :: isc             !< starting i-index for this MPI-domain
    integer              :: jsc             !< starting j-index for this MPI-domain
    integer              :: nx              !< number of points in the i-dir for this MPI-domain
    integer              :: ny              !< number of points in the j-dir for this MPI-domain
    integer              :: levs            !< number of vertical levels
    !--- ak/bk for pressure level calculations
    real(kind=kind_phys), pointer :: ak(:)  !< from surface (k=1) to TOA (k=levs)
    real(kind=kind_phys), pointer :: bk(:)  !< from surface (k=1) to TOA (k=levs)
    integer              :: levsp1          !< number of vertical levels plus one
    integer              :: levsm1          !< number of vertical levels minus one
    integer              :: cnx             !< number of points in the i-dir for this cubed-sphere face
    integer              :: cny             !< number of points in the j-dir for this cubed-sphere face
    integer              :: lonr            !< number of global points in x-dir (i) along the equator
    integer              :: latr            !< number of global points in y-dir (j) along any meridian
    integer              :: tile_num
    integer              :: nblks           !< for explicit data blocking: number of blocks
    integer,     pointer :: blksz(:)        !< for explicit data blocking: block sizes of all blocks
    integer              :: ncols           !< total number of columns for all blocks
    !
    integer              :: nchunks         !< number of chunks of an array that are used in the CCPP run phase
    integer,     pointer :: chunk_begin(:)  !< first indices of chunks of an array for the CCPP run phase
    integer,     pointer :: chunk_end(:)    !< last indices of chunks of an array for the CCPP run phase
    !
    integer              :: fire_aux_data_levels !< vertical levels of fire auxiliary data

!--- coupling parameters
    logical              :: cplflx          !< default no cplflx collection
    logical              :: cplice          !< default no cplice collection (used together with cplflx)
    logical              :: cplocn2atm      !< default yes ocn->atm coupling
    logical              :: cplwav          !< default no cplwav collection
    logical              :: cplwav2atm      !< default no wav->atm coupling
    logical              :: cplaqm          !< default no cplaqm collection
    logical              :: cplchm          !< default no cplchm collection
    logical              :: cpllnd          !< default no cpllnd collection
    logical              :: cpllnd2atm      !< default no lnd->atm coupling 
    logical              :: rrfs_sd         !< default no rrfs_sd collection
    logical              :: cpl_fire        !< default no fire_behavior collection
    logical              :: use_cice_alb    !< default .false. - i.e. don't use albedo imported from the ice model
    logical              :: cpl_imp_mrg     !< default no merge import with internal forcings
    logical              :: cpl_imp_dbg     !< default no write import data to file post merge
    logical              :: use_med_flux    !< default .false. - i.e. don't use atmosphere-ocean fluxes imported from mediator

!--- integrated dynamics through earth's atmosphere
    logical              :: lsidea

!vay 2018  GW physics switches

    logical              :: ldiag_ugwp
    logical              :: ugwp_seq_update ! flag to update winds between UGWP steps
    logical              :: do_ugwp         ! do mesoscale UGWP + TOFD + RF
    logical              :: do_tofd         ! tofd flag in gwdps.f
    logical              :: do_gwd          ! logical for gravity wave drag (gwd)
    logical              :: do_cnvgwd       ! logical for convective gwd

!--- calendars and time parameters and activation triggers
    real(kind=kind_phys) :: dtp             !< physics timestep in seconds
    real(kind=kind_phys) :: dtf             !< dynamics timestep in seconds
    integer              :: nscyc           !< trigger for surface data cycling
    integer              :: nszero          !< trigger for zeroing diagnostic buckets
    integer              :: idat(1:8)       !< initialization date and time
                                            !< (yr, mon, day, t-zone, hr, min, sec, mil-sec)
    integer              :: idate(4)        !< initial date with different size and ordering
                                            !< (hr, mon, day, yr)
!--- radiation control parameters
    real(kind=kind_phys) :: fhswr           !< frequency for shortwave radiation (secs)
    real(kind=kind_phys) :: fhlwr           !< frequency for longwave radiation (secs)
    integer              :: nsswr           !< integer trigger for shortwave radiation
    integer              :: nslwr           !< integer trigger for longwave  radiation
    integer              :: nhfrad          !< number of timesteps for which to call radiation on physics timestep (coldstarts)
    integer              :: levr            !< number of vertical levels for radiation calculations
    integer              :: levrp1          !< number of vertical levels for radiation calculations plus one
    integer              :: nfxr            !< second dimension for fluxr diagnostic variable (radiation)
    logical              :: iaerclm         !< flag for initializing aerosol data
    integer              :: ntrcaer         !< number of aerosol tracers for Morrison-Gettelman microphysics
    logical              :: lmfshal         !< parameter for radiation
    logical              :: lmfdeep2        !< parameter for radiation
    integer              :: nrcm            !< second dimension of random number stream for RAS
    integer              :: iflip           !< iflip - is not the same as flipv
    integer              :: isol            !< use prescribed solar constant
                                            !< 0  => fixed value=1366.0\f$W/m^2\f$(old standard)
                                            !< 10 => fixed value=1360.8\f$W/m^2\f$(new standard)
                                            !< 1  => NOAA ABS-scale TSI table (yearly) w 11-yr cycle approx
                                            !< 2  => NOAA TIM-scale TSI table (yearly) w 11-yr cycle approx
                                            !< 3  => CMIP5 TIM-scale TSI table (yearly) w 11-yr cycle approx
                                            !< 4  => CMIP5 TIM-scale TSI table (monthly) w 11-yr cycle approx
    integer              :: ico2            !< prescribed global mean value (old opernl)
    integer              :: ialb            !< use climatology alb, based on sfc type
                                            !< 1 => use modis based alb
                                            !< 2 => use LSM alb
    integer              :: iems            !< 1 => use fixed value of 1.0
                                            !< 2 => use LSM emiss
    integer              :: iaer            !< default aerosol effect in sw only
    integer              :: iaermdl         !< tropospheric aerosol model scheme flag
    integer              :: iaerflg         !< aerosol effect control flag
    character(len=26)    :: aeros_file      !< external file: aerosol data file
    character(len=26)    :: solar_file      !< external file: solar constant data table
    character(len=26)    :: semis_file      !< external file: surface emissivity data for radiation
    character(len=26)    :: co2dat_file     !< external file: co2 monthly observation data table
    character(len=26)    :: co2gbl_file     !< external file: co2 global annual mean data table
    character(len=26)    :: co2usr_file     !< external file: co2 user defined data table
    character(len=26)    :: co2cyc_file     !< external file: co2 climotological monthly cycle data
    logical              :: lalw1bd         !< selects 1 band or multi bands for LW aerosol properties
    integer              :: icliq_sw        !< sw optical property for liquid clouds
    integer              :: icice_sw        !< sw optical property for ice clouds
    integer              :: icliq_lw        !< lw optical property for liquid clouds
    integer              :: icice_lw        !< lw optical property for ice clouds
    integer              :: iovr            !< cloud-overlap used in cloud-sampling by radiation scheme(s)
    integer              :: ictm            !< ictm=0 => use data at initial cond time, if not
                                            !<           available; use latest; no extrapolation.
                                            !< ictm=1 => use data at the forecast time, if not
                                            !<           available; use latest; do extrapolation.
                                            !< ictm=yyyy0 => use yyyy data for the forecast time;
                                            !<           no extrapolation.
                                            !< ictm=yyyy1 = > use yyyy data for the fcst. If needed,
                                            !<           do extrapolation to match the fcst time.
                                            !< ictm=-1 => use user provided external data for
                                            !<           the fcst time; no extrapolation.
                                            !< ictm=-2 => same as ictm=0, but add seasonal cycle
                                            !<           from climatology; no extrapolation.
    integer              :: isubc_sw        !< sw clouds without sub-grid approximation
    integer              :: isubc_lw        !< lw clouds without sub-grid approximation
                                            !< =1 => sub-grid cloud with prescribed seeds
                                            !< =2 => sub-grid cloud with randomly generated
                                            !< seeds
    integer              :: iswmode         !< SW control flag for scattering process approximation
                                            !< =1 => two-stream delta-eddington (Joseph et al. 1976)
                                            !< =2 => two-stream PIFM            (Zdunkowski et al. 1980)
                                            !< =3 => discrete ordinates         (Liou, 1973)
    integer              :: idcor           !< Decorrelation length type for overlap assumption
                                            !< =0 => Use constant decorrelation length, decorr_con
                                            !< =1 => Use spatially varying decorrelation length (Hogan et al. 2010)
                                            !< =2 => Use spatially and temporally varyint decorrelation length (Oreopoulos et al. 2012)
    real(kind_phys)      :: dcorr_con       !< Decorrelation length constant (km) (if idcor = 0)
    logical              :: lcrick          !< CRICK-Proof cloud water
    logical              :: lcnorm          !< Cloud condensate normalized by cloud cover
    logical              :: lnoprec         !< radiation precip flag for Ferrier/Moorthi
    logical              :: lwhtr           !< flag to output lw heating rate (Radtend%lwhc)
    logical              :: swhtr           !< flag to output sw heating rate (Radtend%swhc)
    integer              :: rad_hr_units    !< flag to control units of lw/sw heating rate
                                            !< 1: K day-1 - 2: K s-1
    logical              :: inc_minor_gas   !< Include minor trace gases in RRTMG radiation calculation?
    integer              :: ipsd0           !< initial permutaion seed for mcica radiation
    integer              :: ipsdlim         !< limit initial permutaion seed for mcica radiation
    logical              :: lrseeds         !< flag to use host-provided random seeds
    integer              :: nrstreams       !< number of random number streams in host-provided random seed array
    logical              :: lextop          !< flag for using an extra top layer for radiation

    ! RRTMGP
    logical              :: do_RRTMGP               !< Use RRTMGP
    character(len=128)   :: active_gases            !< Character list of active gases used in RRTMGP
    integer              :: nGases                  !< Number of active gases
    character(len=128)   :: rrtmgp_root             !< Directory of rte+rrtmgp source code
    character(len=128)   :: lw_file_gas             !< RRTMGP K-distribution file, coefficients to compute optics for gaseous atmosphere
    character(len=128)   :: lw_file_clouds          !< RRTMGP file containing coefficients used to compute clouds optical properties
    integer              :: rrtmgp_nBandsLW         !< Number of RRTMGP LW bands.
    integer              :: rrtmgp_nGptsLW          !< Number of RRTMGP LW spectral points.
    character(len=128)   :: sw_file_gas             !< RRTMGP K-distribution file, coefficients to compute optics for gaseous atmosphere
    character(len=128)   :: sw_file_clouds          !< RRTMGP file containing coefficients used to compute clouds optical properties
    integer              :: rrtmgp_nBandsSW         !< Number of RRTMGP SW bands.
    integer              :: rrtmgp_nGptsSW          !< Number of RRTMGP SW spectral points.
    logical              :: doG_cldoptics           !< Use legacy RRTMG cloud-optics?
    logical              :: doGP_cldoptics_PADE     !< Use RRTMGP cloud-optics: PADE approximation?
    logical              :: doGP_cldoptics_LUT      !< Use RRTMGP cloud-optics: LUTs?
    integer              :: iovr_convcld            !< Cloud-overlap assumption for convective-cloud
    integer              :: rrtmgp_nrghice          !< Number of ice-roughness categories
    integer              :: rrtmgp_nGauss_ang       !< Number of angles used in Gaussian quadrature
    logical              :: do_GPsw_Glw             !< If set to true use rrtmgp for SW calculation, rrtmg for LW.
    character(len=128), pointer :: active_gases_array(:) => null() !< character array for each trace gas name
    logical              :: use_LW_jacobian         !< If true, use Jacobian of LW to update radiation tendency.
    logical              :: damp_LW_fluxadj         !< If true, damp the LW flux adjustment using the Jacobian w/ height with logistic function
    real(kind_phys)      :: lfnc_k                  !<          Logistic function transition depth (Pa)
    real(kind_phys)      :: lfnc_p0                 !<          Logistic function transition level (Pa)
    logical              :: doGP_lwscat             !< If true, include scattering in longwave cloud-optics, only compatible w/ GP cloud-optics
    logical              :: doGP_sgs_cnv            !< If true, include SubGridScale convective cloud in RRTMGP
    logical              :: doGP_sgs_mynn           !< If true, include SubGridScale MYNN-EDMF cloud in RRTMGP
    integer              :: rrtmgp_lw_phys_blksz    !< Number of columns to pass to RRTMGP LW per block.
    integer              :: rrtmgp_sw_phys_blksz    !< Number of columns to pass to RRTMGP SW per block.
    logical              :: doGP_smearclds          !< If true, include implicit SubGridScale clouds in RRTMGP
    real(kind_phys)      :: minGPpres               !< Minimum pressure allowed in RRTMGP.
    real(kind_phys)      :: maxGPpres               !< Maximum pressure allowed in RRTMGP.
    real(kind_phys)      :: minGPtemp               !< Minimum temperature allowed in RRTMGP.
    real(kind_phys)      :: maxGPtemp               !< Maximum temperature allowed in RRTMGP.
    logical              :: top_at_1                !< Vertical ordering flag.
    integer              :: iSFC                    !< Vertical index for surface
    integer              :: iTOA                    !< Vertical index for TOA

!--- microphysical switch
    logical              :: convert_dry_rho = .true.       !< flag for converting mass/number concentrations from moist to dry
                                                           !< for physics options that expect dry mass/number concentrations;
                                                           !< this flag will no longer be needed once the CCPP standard
                                                           !< names and the CCPP framework logic have been augmented to
                                                           !< automatically determine whether such conversions are necessary
                                                           !< and if yes, perform them; hardcoded to .true. for now
    !--- new microphysical switch
    integer              :: imp_physics                    !< choice of microphysics scheme
    integer              :: imp_physics_gfdl          = 11 !< choice of GFDL     microphysics scheme
    integer              :: imp_physics_thompson      = 8  !< choice of Thompson microphysics scheme
    integer              :: imp_physics_wsm6          = 6  !< choice of WSMG     microphysics scheme
    integer              :: imp_physics_zhao_carr     = 99 !< choice of Zhao-Carr microphysics scheme
    integer              :: imp_physics_zhao_carr_pdf = 98 !< choice of Zhao-Carr microphysics scheme with PDF clouds
    integer              :: imp_physics_mg            = 10 !< choice of Morrison-Gettelman microphysics scheme
    integer              :: imp_physics_fer_hires     = 15 !< choice of Ferrier-Aligo microphysics scheme
    integer              :: imp_physics_nssl          = 17 !< choice of NSSL microphysics scheme with background CCN
    integer              :: imp_physics_nssl2mccn     = 18 !< choice of NSSL microphysics scheme with predicted CCN (compatibility)
    integer              :: iovr_rand                 = 0  !< choice of cloud-overlap: random
    integer              :: iovr_maxrand              = 1  !< choice of cloud-overlap: maximum random
    integer              :: iovr_max                  = 2  !< choice of cloud-overlap: maximum
    integer              :: iovr_dcorr                = 3  !< choice of cloud-overlap: decorrelation length
    integer              :: iovr_exp                  = 4  !< choice of cloud-overlap: exponential
    integer              :: iovr_exprand              = 5  !< choice of cloud-overlap: exponential random
    integer              :: idcor_con                 = 0  !< choice for decorrelation-length: Use constant value
    integer              :: idcor_hogan               = 1  !< choice for decorrelation-length: (https://rmets.onlinelibrary.wiley.com/doi/full/10.1002/qj.647)
    integer              :: idcor_oreopoulos          = 2  !< choice for decorrelation-length: (10.5194/acp-12-9097-2012)
    !--- Z-C microphysical parameters
    real(kind=kind_phys) :: psautco(2)         !< [in] auto conversion coeff from ice to snow
    real(kind=kind_phys) :: prautco(2)         !< [in] auto conversion coeff from cloud to rain
    real(kind=kind_phys) :: evpco              !< [in] coeff for evaporation of largescale rain
    real(kind=kind_phys) :: wminco(2)          !< [in] water and ice minimum threshold for Zhao
    real(kind=kind_phys) :: avg_max_length     !< reset time in seconds for max hourly fields
    !--- M-G microphysical parameters
    integer              :: fprcp              !< no prognostic rain and snow (MG)
    integer              :: pdfflag            !< pdf flag for MG macrophysics
    real(kind=kind_phys) :: mg_dcs             !< Morrison-Gettelman microphysics parameters
    real(kind=kind_phys) :: mg_qcvar
    real(kind=kind_phys) :: mg_ts_auto_ice(2)  !< ice auto conversion time scale
    real(kind=kind_phys) :: mg_rhmini          !< relative humidity threshold parameter for nucleating ice

    real(kind=kind_phys) :: mg_ncnst           !< constant droplet num concentration (m-3)
    real(kind=kind_phys) :: mg_ninst           !< constant ice num concentration (m-3)
    real(kind=kind_phys) :: mg_ngnst           !< constant graupel/hail num concentration (m-3)
    real(kind=kind_phys) :: mg_berg_eff_factor !< berg efficiency factor
    real(kind=kind_phys) :: mg_alf             !< tuning factor for alphs in MG macrophysics
    real(kind=kind_phys) :: mg_qcmin(2)        !< min liquid and ice mixing ratio in Mg macro clouds
    character(len=16)    :: mg_precip_frac_method ! type of precipitation fraction method
    real(kind=kind_phys) :: tf
    real(kind=kind_phys) :: tcr
    real(kind=kind_phys) :: tcrf
!
    integer              :: num_dfi_radar      !< number of timespans with radar-prescribed temperature tendencies
    real (kind=kind_phys) :: fh_dfi_radar(1+dfi_radar_max_intervals)   !< begin+end of timespans to receive radar-prescribed temperature tendencies
    logical              :: do_cap_suppress    !< enable convection suppression in GF scheme if fh_dfi_radar is specified
    real (kind=kind_phys) :: radar_tten_limits(2) !< radar_tten values outside this range (min,max) are discarded
    integer              :: ix_dfi_radar(dfi_radar_max_intervals) = -1 !< Index within dfi_radar_tten of each timespan (-1 means "none")
    integer              :: dfi_radar_max_intervals
    integer              :: dfi_radar_max_intervals_plus_one

    !
    logical              :: effr_in            !< eg to turn on ffective radii for MG
    logical              :: microp_uniform
    logical              :: do_cldliq
    logical              :: do_cldice
    logical              :: hetfrz_classnuc

    logical              :: mg_nccons
    logical              :: mg_nicons
    logical              :: mg_ngcons
    logical              :: sed_supersat
    logical              :: do_sb_physics
    logical              :: mg_do_graupel
    logical              :: mg_do_hail
    logical              :: mg_do_ice_gmao
    logical              :: mg_do_liq_liu

    real(kind=kind_phys) :: shoc_parm(5)    !< critical pressure in Pa for tke dissipation in shoc
    integer              :: ncnd            !< number of cloud condensate types

    !--- NSSL microphysics params
    real(kind=kind_phys) :: nssl_cccn      !<  CCN concentration (m-3)
    real(kind=kind_phys) :: nssl_alphah    !<  graupel shape parameter
    real(kind=kind_phys) :: nssl_alphahl   !<  hail shape parameter
    real(kind=kind_phys) :: nssl_alphar    ! shape parameter for rain (imurain=1 only)
    real(kind=kind_phys) :: nssl_ehw0      ! constant or max assumed graupel-droplet collection efficiency
    real(kind=kind_phys) :: nssl_ehlw0     ! constant or max assumed hail-droplet collection efficiency
    logical              :: nssl_hail_on   !<  NSSL flag to activate the hail category
    logical              :: nssl_ccn_on    !<  NSSL flag to activate the CCN category
    logical              :: nssl_invertccn !<  NSSL flag to treat CCN as activated (true) or unactivated (false)
    logical              :: nssl_3moment   !<  NSSL flag to turn on 3-moment for rain/graupel/hail

    !--- Thompson's microphysical parameters
    logical              :: ltaerosol       !< flag for aerosol version
    logical              :: mraerosol       !< flag for merra2_aerosol_aware
    logical              :: lradar          !< flag for radar reflectivity
    real(kind=kind_phys) :: nsfullradar_diag!< seconds between resetting radar reflectivity calculation
    real(kind=kind_phys) :: ttendlim        !< temperature tendency limiter per time step in K/s
    logical              :: ext_diag_thompson !< flag for extended diagnostic output from Thompson
    integer              :: thompson_ext_ndiag3d=37 !< number of 3d arrays for extended diagnostic output from Thompson
    real(kind=kind_phys) :: dt_inner        !< time step for the inner loop in s
    logical              :: sedi_semi       !< flag for semi Lagrangian sedi of rain
    integer              :: decfl           !< deformed CFL factor
    logical              :: thpsnmp_is_init !< Local scheme initialization flag

    !--- GFDL microphysical paramters
    logical              :: lgfdlmprad      !< flag for GFDL mp scheme and radiation consistency

    !--- Thompson,GFDL mp parameter
    logical              :: lrefres          !< flag for radar reflectivity in restart file

    !--- land/surface model parameters
    integer              :: lsm             !< flag for land surface model lsm=1 for noah lsm
    integer              :: lsm_noah=1      !< flag for NOAH land surface model
    integer              :: lsm_noahmp=2    !< flag for NOAH land surface model
    integer              :: lsm_ruc=3       !< flag for RUC land surface model
    integer              :: lsoil           !< number of soil layers
    integer              :: ivegsrc         !< ivegsrc = 0   => USGS,
                                            !< ivegsrc = 1   => IGBP (20 category)
                                            !< ivegsrc = 2   => UMD  (13 category)
                                            !< ivegsrc = 3   => NLCD40 (40 category, NOAH WRFv4 only)
                                            !< ivegsrc = 4   => USGS-RUC (28 category, NOAH WRFv4 only)
                                            !< ivegsrc = 5   => MODI-RUC (21 category, NOAH WRFv4 only)
    integer              :: nvegcat         !< nvegcat = 20 if ivegsrc = 1
    integer              :: isot            !< isot = 0   => Zobler soil type  ( 9 category)
                                            !< isot = 1   => STATSGO soil type (19 category, AKA 'STAS'(?))
                                            !< isot = 2   => STAS-RUC soil type (19 category, NOAH WRFv4 only)
    integer              :: nsoilcat        !< nsoilcat = 19 if isot = 1
    integer              :: kice            !< number of layers in sice
    integer              :: lsoil_lsm       !< number of soil layers internal to land surface model
    integer              :: lsnow_lsm       !< maximum number of snow layers internal to land surface model
    integer              :: lsnow_lsm_lbound!< lower bound for snow arrays, depending on lsnow_lsm
    integer              :: lsnow_lsm_ubound!< upper bound for snow arrays, depending on lsnow_lsm
    logical              :: exticeden       !< flag for calculating frozen precip ice density outside of the LSM
    real(kind=kind_phys), pointer :: zs(:)    => null() !< depth of soil levels for land surface model
    real(kind=kind_phys), pointer :: dzs(:)   => null() !< thickness of soil levels for land surface model
    real(kind=kind_phys), pointer :: pores(:) => null() !< max soil moisture for a given soil type for land surface model
    real(kind=kind_phys), pointer :: resid(:) => null() !< min soil moisture for a given soil type for land surface model
    logical              :: rdlai           !< read LAI from input file (for RUC LSM or NOAH LSM WRFv4)
    logical              :: ua_phys         !< flag for using University of Arizona? extension to NOAH LSM WRFv4
    logical              :: usemonalb       !< flag to read surface diffused shortwave albedo from input file for NOAH LSM WRFv4
    real(kind=kind_phys) :: aoasis          !< potential evaporation multiplication factor for NOAH LSM WRFv4
    integer              :: fasdas          !< flag to use "flux-adjusting surface data assimilation system"; 0 = OFF, 1 = ON
    integer              :: iopt_thcnd      !< option to treat thermal conductivity in Noah LSM (new in 3.8)
                                            !< = 1, original (default)
                                            !< = 2, McCumber and Pielke for silt loam and sandy loam

    ! -- the Noah MP options
    integer              :: iopt_dveg ! 1-> off table lai 2-> on 3-> off;4->off;5 -> on
    integer              :: iopt_crs  !canopy stomatal resistance (1-> ball-berry; 2->jarvis)
    integer              :: iopt_btr  !soil moisture factor for stomatal resistance (1-> noah; 2-> clm; 3-> ssib)
    integer              :: iopt_run  !runoff and groundwater (1->simgm; 2->simtop; 3->schaake96; 4->bats)
    integer              :: iopt_sfc  !surface layer drag coeff (ch & cm) (1->m-o; 2->chen97)
    integer              :: iopt_frz  !supercooled liquid water (1-> ny06; 2->koren99)
    integer              :: iopt_inf  !frozen soil permeability (1-> ny06; 2->koren99)
    integer              :: iopt_rad  !radiation transfer (1->gap=f(3d,cosz); 2->gap=0; 3->gap=1-fveg)
    integer              :: iopt_alb  !snow surface albedo (1->bats; 2->class)
    integer              :: iopt_snf  !rainfall & snowfall (1-jordan91; 2->bats; 3->noah)
    integer              :: iopt_tbot !lower boundary of soil temperature (1->zero-flux; 2->noah)
    integer              :: iopt_stc  !snow/soil temperature time scheme (only layer 1)
    integer              :: iopt_trs  !thermal roughness scheme (1-z0h=z0m; 2-czil; 3-ec;4-kb inversed)
    integer              :: iopt_diag !2m t/q diagnostic approach (1->external GFS sfc_diag 2->original NoahMP 2-title 3->NoahMP
                                      !2-title + internal GFS sfc_diag  )

    ! -- RUC LSM options
    integer              :: mosaic_lu=0     !< control for use of fractional landuse in RUC land surface model
    integer              :: mosaic_soil=0   !< control for use of fractional soil in RUC land surface model
    integer              :: isncond_opt=1   !< control for soil thermal conductivity option in RUC land surface model
    integer              :: isncovr_opt=1   !< control for snow cover fraction option in RUC land surface model

    ! -- Fire heat flux
    logical              :: add_fire_heat_flux=.false. !<control to add fireheat flux to RUC LSM

    logical              :: use_ufo         !< flag for gcycle surface option

    ! GFDL Surface Layer options
    logical              :: lcurr_sf        !< flag for taking ocean currents into account in GFDL surface layer
    logical              :: pert_cd         !< flag for perturbing the surface drag coefficient for momentum in surface layer scheme (1 = True)
    integer              :: ntsflg          !< flag for updating skin temperature in the GFDL surface layer scheme
    real(kind=kind_phys) :: sfenth          !< enthalpy flux factor 0 zot via charnock ..>0 zot enhanced>15m/s

!--- lake model parameters
    integer              :: lkm             !< =0 no lake, =1 lake, =2 lake&nsst
    integer              :: iopt_lake       !< =1 flake, =2 clm lake
    integer              :: iopt_lake_flake = 1
    integer              :: iopt_lake_clm = 2
    real(kind_phys)      :: lakedepth_threshold !< lakedepth must be GREATER than this value to enable a lake model
    real(kind_phys)      :: lakefrac_threshold  !< lakefrac must be GREATER than this value to enable a lake model
    logical              :: use_lake2m      !< use 2m T & Q calculated by the lake model

!--- clm lake model parameters
    integer              :: nlevlake_clm_lake !< Number of lake levels for clm lake model
    integer              :: nlevsoil_clm_lake !< Number of soil levels for clm lake model
    integer              :: nlevsnow_clm_lake !< Number of snow levels for clm lake model
    integer              :: nlevsnowsoil_clm_lake !< -nlevsnow:nlevsoil dimensioned variables
    integer              :: nlevsnowsoil1_clm_lake !< -nlevsnow+1:nlevsoil dimensioned variables
    real(kind_phys)      :: clm_lake_depth_default !< minimum lake elevation in clm lake model
    logical              :: clm_lake_use_lakedepth !< initialize lake from lakedepth
    logical              :: clm_lake_debug !< verbose debugging in clm_lake
    logical              :: clm_debug_print !< enables prints in clm_lakedebugging in clm_laki

!--- tuning parameters for physical parameterizations
    logical              :: ras             !< flag for ras convection scheme
    logical              :: flipv           !< flag for vertical direction flip (ras)
                                            !< .true. implies surface at k=1
    logical              :: trans_trac      !< flag for convective transport of tracers (RAS, CS, or SAMF)
    logical              :: old_monin       !< flag for diff monin schemes
    logical              :: cnvgwd          !< flag for conv gravity wave drag
    integer              :: gwd_opt         !< gwd_opt = 1  => original GFS gwd (gwdps.f)
                                            !< gwd_opt = 2  => unified ugwp GWD
                                            !< gwd_opt = 22 => unified ugwp GWD with extra output
                                            !< gwd_opt = 3  => GSL drag suite
                                            !< gwd_opt = 33 => GSL drag suite with extra output
    logical              :: do_ugwp_v0           !< flag for version 0 ugwp GWD
    logical              :: do_ugwp_v0_orog_only !< flag for version 0 ugwp GWD (orographic drag only)
    logical              :: do_ugwp_v0_nst_only  !< flag for version 0 ugwp GWD (non-stationary GWD only)
    logical              :: do_gsl_drag_ls_bl    !< flag for GSL drag (mesoscale GWD and blocking only)
    logical              :: do_gsl_drag_ss       !< flag for GSL drag (small-scale GWD only)
    logical              :: do_gsl_drag_tofd     !< flag for GSL drag (turbulent orog form drag only)
    logical              :: do_gwd_opt_psl       !< flag for PSL drag (mesoscale GWD and blocking only)
    logical              :: do_ugwp_v1           !< flag for version 1 ugwp GWD
    logical              :: do_ugwp_v1_orog_only !< flag for version 1 ugwp GWD (orographic drag only)
    logical              :: do_ugwp_v1_w_gsldrag !< flag for version 1 ugwp with OGWD of GSL
    logical              :: mstrat          !< flag for moorthi approach for stratus
    logical              :: moist_adj       !< flag for moist convective adjustment
    logical              :: cscnv           !< flag for Chikira-Sugiyama convection
    logical              :: cal_pre         !< flag controls precip type algorithm
    real(kind=kind_phys) :: rhgrd           !< fer_hires microphysics only
    logical              :: spec_adv        !< flag for individual cloud species advected
    integer              :: icloud          !< cloud effect to the optical depth in radiation; this also controls the cloud fraction options
                                            !<  3: with cloud effect, and use cloud fraction option 3, based on Sundqvist et al. (1989)
    logical              :: do_aw           !< AW scale-aware option in cs convection
    logical              :: do_awdd         !< AW scale-aware option in cs convection
    logical              :: flx_form        !< AW scale-aware option in cs convection
    logical              :: do_shoc         !< flag for SHOC
    logical              :: shocaftcnv      !< flag for SHOC
    logical              :: shoc_cld        !< flag for clouds
    logical              :: uni_cld         !< flag for clouds in grrad
    logical              :: h2o_phys        !< flag for stratosphere h2o
    logical              :: pdfcld          !< flag for pdfcld
    logical              :: shcnvcw         !< flag for shallow convective cloud
    logical              :: redrag          !< flag for reduced drag coeff. over sea
    logical              :: hybedmf         !< flag for hybrid edmf pbl scheme
    logical              :: satmedmf        !< flag for scale-aware TKE-based moist edmf
                                            !< vertical turbulent mixing scheme
    logical              :: shinhong        !< flag for scale-aware Shinhong vertical turbulent mixing scheme
    logical              :: do_ysu          !< flag for YSU turbulent mixing scheme
    logical              :: dspheat         !< flag for tke dissipative heating
    logical              :: hurr_pbl        !< flag for hurricane-specific options in PBL scheme
    logical              :: lheatstrg       !< flag for canopy heat storage parameterization
    logical              :: lseaspray       !< flag for sea spray parameterization
    logical              :: cnvcld
    logical              :: xr_cnvcld       !< flag for adding suspended convective clouds to Xu-Randall cloud fraction
    logical              :: random_clds     !< flag controls whether clouds are random
    logical              :: shal_cnv        !< flag for calling shallow convection
    logical              :: do_deep         !< whether to do deep convection
    integer              :: imfshalcnv      !< flag for mass-flux shallow convection scheme
                                            !<     1: July 2010 version of mass-flux shallow conv scheme
                                            !<         current operational version as of 2016
                                            !<     2: scale- & aerosol-aware mass-flux shallow conv scheme (2017)
                                            !<     3: scale- & aerosol-aware Grell-Freitas scheme (GSD)
                                            !<     4: New Tiedtke scheme (CAPS)
                                            !<     0: modified Tiedtke's eddy-diffusion shallow conv scheme
                                            !<    -1: no shallow convection used
    integer              :: imfshalcnv_sas      = 1 !< flag for SAS mass-flux shallow convection scheme
    integer              :: imfshalcnv_samf     = 2 !< flag for SAMF scale- & aerosol-aware mass-flux shallow convection scheme
    integer              :: imfshalcnv_gf       = 3 !< flag for scale- & aerosol-aware Grell-Freitas scheme (GSD)
    integer              :: imfshalcnv_ntiedtke = 4 !< flag for new Tiedtke scheme (CAPS)
    integer              :: imfshalcnv_c3       = 5 !< flag for the Community Convective Cloud (C3) scheme
    logical              :: hwrf_samfdeep           !< flag for HWRF SAMF deepcnv scheme (HWRF)
    logical              :: progsigma               !< flag for prognostic area fraction in samf ddepcnv scheme (GFS)
    integer              :: imfdeepcnv      !< flag for mass-flux deep convection scheme
                                            !<     1: July 2010 version of SAS conv scheme
                                            !<           current operational version as of 2016
                                            !<     2: scale- & aerosol-aware mass-flux deep conv scheme (2017)
                                            !<     3: scale- & aerosol-aware Grell-Freitas scheme (GSD)
                                            !<     4: New Tiedtke scheme (CAPS)
                                            !<     0: old SAS Convection scheme before July 2010
    integer              :: imfdeepcnv_sas      = 1 !< flag for SAS mass-flux deep convection scheme
    integer              :: imfdeepcnv_samf     = 2 !< flag for SAMF scale- & aerosol-aware mass-flux deep convection scheme
    integer              :: imfdeepcnv_gf       = 3 !< flag for scale- & aerosol-aware Grell-Freitas scheme (GSD)
    integer              :: imfdeepcnv_ntiedtke = 4 !< flag for new Tiedtke scheme (CAPS)
    integer              :: imfdeepcnv_c3       = 5 !< flag for the Community Convective Cloud (C3) scheme
    logical              :: hwrf_samfshal           !< flag for HWRF SAMF shalcnv scheme (HWRF)
    integer              :: isatmedmf       !< flag for scale-aware TKE-based moist edmf scheme
                                            !<     0: initial version of satmedmf (Nov. 2018)
                                            !<     1: updated version of satmedmf (as of May 2019)
    integer              :: isatmedmf_vdif  = 0 !< flag for initial version of satmedmf (Nov. 2018)
    integer              :: isatmedmf_vdifq = 1 !< flag for updated version of satmedmf (as of May 2019)
    integer              :: ichoice         = 0 !< flag for closure of C3/GF deep convection
    integer              :: ichoicem        = 13!< flag for closure of C3/GF mid convection
    integer              :: ichoice_s       = 3 !< flag for closure of C3/GF shallow convection

    integer              :: nmtvr           !< number of topographic variables such as variance etc
                                            !< used in the GWD parameterization - 10 more added if
                                            !< GSL orographic drag scheme is used
    integer              :: jcap            !< number of spectral wave trancation used only by sascnv shalcnv
    real(kind=kind_phys) :: cs_parm(10)     !< tunable parameters for Chikira-Sugiyama convection
    real(kind=kind_phys) :: flgmin(2)       !< [in] ice fraction bounds
    real(kind=kind_phys) :: cgwf(2)         !< multiplication factor for convective GWD
    real(kind=kind_phys) :: ccwf(2)         !< multiplication factor for critical cloud
                                            !< workfunction for RAS
    real(kind=kind_phys) :: cdmbgwd(4)      !< multiplication factors for cdmb, gwd and NS gwd, tke based enhancement
    real(kind=kind_phys) :: alpha_fd        !< alpha coefficient for turbulent orographic form drag
    real(kind=kind_phys) :: psl_gwd_dx_factor  !< multiplication factors for grid spacing
    real(kind=kind_phys) :: sup             !< supersaturation in pdf cloud when t is very low
    real(kind=kind_phys) :: ctei_rm(2)      !< critical cloud top entrainment instability criteria
                                            !< (used if mstrat=.true.)
    real(kind=kind_phys) :: crtrh(3)        !< critical relative humidity at the surface
                                            !< PBL top and at the top of the atmosphere
    real(kind=kind_phys) :: dlqf(2)         !< factor for cloud condensate detrainment
                                            !< from cloud edges for RAS
    real(kind=kind_phys) :: psauras(2)      !< [in] auto conversion coeff from ice to snow in ras
    real(kind=kind_phys) :: prauras(2)      !< [in] auto conversion coeff from cloud to rain in ras
    real(kind=kind_phys) :: wminras(2)      !< [in] water and ice minimum threshold for ras

    integer              :: seed0           !< random seed for radiation

    real(kind=kind_phys) :: rbcr            !< Critical Richardson Number in the PBL scheme
    real(kind=kind_phys) :: betascu         !< Tuning parameter for prog. closure shallow clouds
    real(kind=kind_phys) :: betamcu         !< Tuning parameter for prog. closure midlevel clouds 
    real(kind=kind_phys) :: betadcu         !< Tuning parameter for prog. closure deep clouds 

    !--- MYNN parameters/switches
    logical              :: do_mynnedmf
    logical              :: do_mynnsfclay
    ! DH* TODO - move this to MYNN namelist section
    integer              :: tke_budget         !< flag for activating TKE budget
    logical              :: bl_mynn_tkeadvect  !< activate computation of TKE advection (not yet in use for FV3)
    integer              :: bl_mynn_cloudpdf   !< flag to determine which cloud PDF to use
    integer              :: bl_mynn_mixlength  !< flag for different version of mixing length formulation
    integer              :: bl_mynn_edmf       !< flag to activate the mass-flux scheme
    integer              :: bl_mynn_edmf_mom   !< flag to activate the transport of momentum
    integer              :: bl_mynn_edmf_tke   !< flag to activate the transport of TKE
    integer              :: bl_mynn_cloudmix   !< flag to activate mixing of cloud species
    integer              :: bl_mynn_mixqt      !< flag to mix total water or individual species
    integer              :: bl_mynn_output     !< flag to initialize and write out extra 3D arrays
    integer              :: icloud_bl          !< flag for coupling sgs clouds to radiation
    real(kind=kind_phys) :: bl_mynn_closure    !< flag to determine closure level of MYNN
    logical              :: sfclay_compute_flux!< flag for thermal roughness lengths over water in mynnsfclay
    logical              :: sfclay_compute_diag!< flag for computing surface diagnostics in mynnsfclay
    integer              :: isftcflx           !< flag for thermal roughness lengths over water in mynnsfclay
    integer              :: iz0tlnd            !< flag for thermal roughness lengths over land in mynnsfclay
    real(kind=kind_phys) :: var_ric
    real(kind=kind_phys) :: coef_ric_l
    real(kind=kind_phys) :: coef_ric_s
    ! *DH
    ! MYJ switches
    logical              :: do_myjsfc          !< flag for MYJ surface layer scheme
    logical              :: do_myjpbl          !< flag for MYJ PBL scheme

!--- Rayleigh friction
    real(kind=kind_phys) :: prslrd0         !< pressure level from which Rayleigh Damping is applied
    real(kind=kind_phys) :: ral_ts          !< time scale for Rayleigh damping in days

!--- mass flux deep convection
    real(kind=kind_phys) :: clam_deep       !< c_e for deep convection (Han and Pan, 2011, eq(6))
    real(kind=kind_phys) :: c0s_deep        !< convective rain conversion parameter
    real(kind=kind_phys) :: c1_deep         !< conversion parameter of detrainment from liquid water into grid-scale cloud water
    real(kind=kind_phys) :: betal_deep      !< fraction factor of downdraft air mass reaching ground surface over land
    real(kind=kind_phys) :: betas_deep      !< fraction factor of downdraft air mass reaching ground surface over sea
    real(kind=kind_phys) :: evef            !< evaporation factor from convective rain
    real(kind=kind_phys) :: evfact_deep     !< evaporation factor from convective rain
    real(kind=kind_phys) :: evfactl_deep    !< evaporation factor from convective rain over land
    real(kind=kind_phys) :: pgcon_deep      !< reduction factor in momentum transport due to convection induced pressure gradient force
                                            !< 0.7 : Gregory et al. (1997, QJRMS)
                                            !< 0.55: Zhang & Wu (2003, JAS)
    real(kind=kind_phys) :: asolfac_deep    !< aerosol-aware parameter based on Lim (2011)
                                            !< asolfac= cx / c0s(=.002)
                                            !< cx = min([-0.7 ln(Nccn) + 24]*1.e-4, c0s)
                                            !< Nccn: CCN number concentration in cm^(-3)
                                            !< Until a realistic Nccn is provided, Nccns are assumed
                                            !< as Nccn=100 for sea and Nccn=1000 for land

!--- mass flux shallow convection
    real(kind=kind_phys) :: clam_shal       !< c_e for shallow convection (Han and Pan, 2011, eq(6))
    real(kind=kind_phys) :: c0s_shal        !< convective rain conversion parameter
    real(kind=kind_phys) :: c1_shal         !< conversion parameter of detrainment from liquid water into grid-scale cloud water
    real(kind=kind_phys) :: pgcon_shal      !< reduction factor in momentum transport due to convection induced pressure gradient force
                                            !< 0.7 : Gregory et al. (1997, QJRMS)
                                            !< 0.55: Zhang & Wu (2003, JAS)
    real(kind=kind_phys) :: asolfac_shal    !< aerosol-aware parameter based on Lim (2011)
                                            !< asolfac= cx / c0s(=.002)
                                            !< cx = min([-0.7 ln(Nccn) + 24]*1.e-4, c0s)
                                            !< Nccn: CCN number concentration in cm^(-3)
                                            !< Until a realistic Nccn is provided, Nccns are assumed
                                            !< as Nccn=100 for sea and Nccn=1000 for land

!--- near surface temperature model
    logical              :: nst_anl         !< flag for NSSTM analysis in gcycle/sfcsub
    integer              :: lsea
    integer              :: nstf_name(5)    !< flag 0 for no nst  1 for uncoupled nst  and 2 for coupled NST
                                            !< nstf_name contains the NSST related parameters
                                            !< nstf_name(1) : 0 = NSSTM off, 1 = NSSTM on but uncoupled
                                            !<                2 = NSSTM on and coupled
                                            !< nstf_name(2) : 1 = NSSTM spin up on, 0 = NSSTM spin up off
                                            !< nstf_name(3) : 1 = NSST analysis on, 0 = NSSTM analysis off
                                            !< nstf_name(4) : zsea1 in mm
                                            !< nstf_name(5) : zsea2 in mm
!--- fractional grid
    logical              :: frac_grid       !< flag for fractional grid
    logical              :: frac_ice        !< flag for fractional ice when fractional grid is not in use
    logical              :: ignore_lake     !< flag for ignoring lakes
    real(kind=kind_phys) :: min_lakeice     !< minimum lake ice value
    real(kind=kind_phys) :: min_seaice      !< minimum sea  ice value
    real(kind=kind_phys) :: min_lake_height !< minimum lake height value
    real(kind=kind_phys) :: rho_h2o         !< density of fresh water

!--- surface layer z0 scheme
    integer              :: sfc_z0_type     !< surface roughness options over ocean:
                                            !< 0=no change
                                            !< 6=areodynamical roughness over water with input 10-m wind
                                            !< 7=slightly decrease Cd for higher wind speed compare to 6
!--- air_sea_flux scheme
    integer              :: icplocn2atm     !< air_sea flux options over ocean:
                                            !< 0=no change
                                            !< l=including ocean current in the computation of air_sea fluxes

!--- potential temperature definition in surface layer physics
    logical              :: thsfc_loc       !< flag for local vs. standard potential temperature
!--- flux method in 2-m diagnostics
    logical              :: diag_flux       !< flag for flux method in 2-m diagnostics
!--- log method in 2-m diagnostics (for stable conditions)
    logical              :: diag_log        !< flag for log method in 2-m diagnostics (for stable conditions)

!--- vertical diffusion
    real(kind=kind_phys) :: xkzm_m          !< [in] bkgd_vdif_m  background vertical diffusion for momentum
    real(kind=kind_phys) :: xkzm_h          !< [in] bkgd_vdif_h  background vertical diffusion for heat q
    real(kind=kind_phys) :: xkzm_s          !< [in] bkgd_vdif_s  sigma threshold for background mom. diffusion
    real(kind=kind_phys) :: xkzminv         !< diffusivity in inversion layers
    real(kind=kind_phys) :: moninq_fac      !< turbulence diffusion coefficient factor
    real(kind=kind_phys) :: dspfac          !< tke dissipative heating factor
    real(kind=kind_phys) :: bl_upfr         !< updraft fraction in boundary layer mass flux scheme
    real(kind=kind_phys) :: bl_dnfr         !< downdraft fraction in boundary layer mass flux scheme
    real(kind=kind_phys) :: rlmx            !< maximum allowed mixing length in boundary layer mass flux scheme
    real(kind=kind_phys) :: elmx            !< maximum allowed dissipation mixing length in boundary layer mass flux scheme
    integer              :: sfc_rlm         !< choice of near surface mixing length in boundary layer mass flux scheme
    integer              :: tc_pbl          !< control for TC applications in the PBL scheme

!--- parameters for canopy heat storage (CHS) parameterization
    real(kind=kind_phys) :: h0facu          !< CHS factor for sensible heat flux in unstable surface layer
    real(kind=kind_phys) :: h0facs          !< CHS factor for sensible heat flux in stable surface layer

!---cellular automata control parameters
    integer              :: nca             !< number of independent cellular automata
    integer              :: nlives          !< cellular automata lifetime
    integer              :: ncells          !< cellular automata finer grid
    integer              :: nca_g           !< number of independent cellular automata
    integer              :: nlives_g        !< cellular automata lifetime
    integer              :: ncells_g        !< cellular automata finer grid
    real(kind=kind_phys) :: nfracseed       !< cellular automata seed probability
    integer              :: nseed           !< cellular automata seed frequency
    integer              :: nseed_g         !< cellular automata seed frequency
    logical              :: do_ca           !< cellular automata main switch
    logical              :: ca_advect       !< Advection of cellular automata
    logical              :: ca_sgs          !< switch for sgs ca
    logical              :: ca_global       !< switch for global ca
    logical              :: ca_smooth       !< switch for gaussian spatial filter
    integer(kind=kind_dbl_prec) :: iseed_ca        !< seed for random number generation in ca scheme
    integer              :: nspinup         !< number of iterations to spin up the ca
    real(kind=kind_phys) :: nthresh         !< threshold used for convection coupling
    real                 :: ca_amplitude    !< amplitude of ca trigger perturbation
    integer              :: nsmooth         !< number of passes through smoother
    logical              :: ca_closure      !< logical switch for ca on closure
    logical              :: ca_entr         !< logical switch for ca on entrainment
    logical              :: ca_trigger      !< logical switch for ca on trigger
    real (kind=kind_phys), allocatable :: vfact_ca(:) !< vertical tapering for ca_global

!--- stochastic physics control parameters
    logical              :: do_sppt
    logical              :: pert_clds
    logical              :: pert_radtend
    logical              :: pert_mp
    logical              :: use_zmtnblck
    logical              :: do_shum
    logical              :: do_skeb
    integer              :: skeb_npass
    integer              :: lndp_type         ! integer indicating land perturbation scheme type:
                                              ! 0 - none
                                              ! 1 - scheme from Gehne et al, MWR, 2019.  (Noah only, not maintained?)
                                              ! 2 - scheme from Draper, JHM, 2021.
    real(kind=kind_phys) :: sppt_amp          ! pjp cloud perturbations
    integer              :: n_var_lndp
    logical              :: lndp_each_step    ! flag to indicate that land perturbations are applied at every time step,
                                              ! otherwise they are applied only
                                              ! after gcycle is run

    ! next two are duplicated here to support lndp_type=1. If delete that scheme, could remove from GFS defs?
    character(len=3)    , pointer :: lndp_var_list(:)
    real(kind=kind_phys), pointer :: lndp_prt_list(:)
    logical              :: do_spp            ! Overall flag to turn on SPP or not
    integer              :: spp_pbl
    integer              :: spp_sfc
    integer              :: spp_mp
    integer              :: spp_rad
    integer              :: spp_gwd
    integer              :: spp_cu_deep
    integer              :: n_var_spp
    character(len=10)    , pointer :: spp_var_list(:)
    real(kind=kind_phys), pointer :: spp_prt_list(:)
    real(kind=kind_phys), pointer :: spp_stddev_cutoff(:)

!--- tracer handling
    character(len=32), pointer :: tracer_names(:) !< array of initialized tracers from dynamic core
    integer              :: ntrac                 !< number of tracers
    integer              :: ntracp1               !< number of tracers plus one
    integer              :: ntracp100             !< number of tracers plus one hundred
    integer              :: nqrimef               !< tracer index for mass weighted rime factor

    integer, pointer :: dtidx(:,:) => null()                                !< index in outermost dimension of dtend
    integer :: ndtend                                                       !< size of outermost dimension of dtend
    type(dtend_var_label), pointer :: dtend_var_labels(:) => null()         !< information about first dim of dtidx
    type(dtend_process_label), pointer :: dtend_process_labels(:) => null() !< information about second dim of dtidx

    ! Indices within inner dimension of dtidx for things that are not tracers:
    integer :: index_of_temperature  !< temperature in dtidx
    integer :: index_of_x_wind       !< x wind in dtidx
    integer :: index_of_y_wind       !< y wind in dtidx

    ! Indices within outer dimension of dtidx:
    integer :: nprocess                         !< maximum value of the below index_for_process_ variables
    integer :: nprocess_summed                  !< number of causes in dtend(:,:,dtidx(...)) to sum to make the physics tendency
    integer :: index_of_process_pbl              !< tracer changes caused by PBL scheme
    integer :: index_of_process_dcnv             !< tracer changes caused by deep convection scheme
    integer :: index_of_process_scnv             !< tracer changes caused by shallow convection scheme
    integer :: index_of_process_mp               !< tracer changes caused by microphysics scheme
    integer :: index_of_process_prod_loss        !< tracer changes caused by ozone production and loss
    integer :: index_of_process_ozmix            !< tracer changes caused by ozone mixing ratio
    integer :: index_of_process_temp             !< tracer changes caused by temperature
    integer :: index_of_process_longwave         !< tracer changes caused by long wave radiation
    integer :: index_of_process_shortwave        !< tracer changes caused by short wave radiation
    integer :: index_of_process_orographic_gwd   !< tracer changes caused by orographic gravity wave drag
    integer :: index_of_process_rayleigh_damping !< tracer changes caused by Rayleigh damping
    integer :: index_of_process_nonorographic_gwd   !< tracer changes caused by convective gravity wave drag
    integer :: index_of_process_overhead_ozone   !< tracer changes caused by overhead ozone column
    integer :: index_of_process_conv_trans       !< tracer changes caused by convective transport
    integer :: index_of_process_physics          !< tracer changes caused by physics schemes
    integer :: index_of_process_non_physics      !< tracer changes caused by everything except physics schemes
    integer :: index_of_process_dfi_radar        !< tracer changes caused by radar mp temperature tendency forcing
    integer :: index_of_process_photochem        !< all changes to ozone
    logical, pointer :: is_photochem(:) => null()!< flags for which processes should be summed as photochemical

    integer              :: ntqv            !< tracer index for water vapor (specific humidity)
    integer              :: ntoz            !< tracer index for ozone mixing ratio
    integer              :: ntcw            !< tracer index for cloud condensate (or liquid water)
    integer              :: ntiw            !< tracer index for ice water
    integer              :: ntrw            !< tracer index for rain water
    integer              :: ntsw            !< tracer index for snow water
    integer              :: ntgl            !< tracer index for graupel
    integer              :: nthl            !< tracer index for hail
    integer              :: ntclamt         !< tracer index for cloud amount
    integer              :: ntlnc           !< tracer index for liquid number concentration
    integer              :: ntinc           !< tracer index for ice    number concentration
    integer              :: ntrnc           !< tracer index for rain   number concentration
    integer              :: ntsnc           !< tracer index for snow   number concentration
    integer              :: ntgnc           !< tracer index for graupel number concentration
    integer              :: nthnc           !< tracer index for hail number concentration
    integer              :: ntccn           !< tracer index for CCN
    integer              :: ntccna          !< tracer index for activated CCN
    integer              :: ntgv            !< tracer index for graupel particle volume
    integer              :: nthv            !< tracer index for hail particle volume
    integer              :: ntrz            !< tracer index for rain reflectivity
    integer              :: ntgz            !< tracer index for graupel reflectivity
    integer              :: nthz            !< tracer index for hail reflectivity
    integer              :: ntke            !< tracer index for kinetic energy
    integer              :: ntsigma         !< tracer index for updraft area fraction
    integer              :: nto             !< tracer index for oxygen ion
    integer              :: nto2            !< tracer index for oxygen
    integer              :: ntwa            !< tracer index for water friendly aerosol
    integer              :: ntia            !< tracer index for ice friendly aerosol
    integer              :: ntfsmoke        !< tracer index for fire smoke
    integer              :: ntsmoke         !< tracer index for smoke
    integer              :: ntdust          !< tracer index for dust
    integer              :: ntcoarsepm      !< tracer index for coarse PM
    integer              :: nchem           !< number of prognostic chemical species (vertically mixied)
    integer              :: ndvel           !< number of prognostic chemical species (which are deposited, usually =nchem)
    integer              :: ntchm           !< number of prognostic chemical tracers (advected)
    integer              :: ntchs           !< tracer index for first prognostic chemical tracer
    integer              :: ntche           !< tracer index for last prognostic chemical tracer
    integer              :: ntdu1           !< tracer index for dust bin1
    integer              :: ntdu2           !< tracer index for dust bin2
    integer              :: ntdu3           !< tracer index for dust bin3
    integer              :: ntdu4           !< tracer index for dust bin4
    integer              :: ntdu5           !< tracer index for dust bin5
    integer              :: ntss1           !< tracer index for sea salt bin1
    integer              :: ntss2           !< tracer index for sea salt bin2
    integer              :: ntss3           !< tracer index for sea salt bin3
    integer              :: ntss4           !< tracer index for sea salt bin4
    integer              :: ntss5           !< tracer index for sea salt bin5
    integer              :: ntsu            !< tracer index for sulfate
    integer              :: ntbcl           !< tracer index for BCPHILIC
    integer              :: ntbcb           !< tracer index for BCPHOBIC
    integer              :: ntocl           !< tracer index for OCPHILIC
    integer              :: ntocb           !< tracer index for OCPHOBIC
    integer              :: ndchm           !< number of diagnostic chemical tracers (not advected)
    integer              :: ndchs           !< tracer index for first diagnostic chemical tracer
    integer              :: ndche           !< tracer index for last diagnostic chemical tracer
    logical, pointer     :: ntdiag(:) => null() !< array to control diagnostics for chemical tracers
    real(kind=kind_phys), pointer :: fscav(:)  => null() !< array of aerosol scavenging coefficients

    !--- derived totals for phy_f*d
    integer              :: ntot2d          !< total number of variables for phyf2d
    integer              :: ntot3d          !< total number of variables for phyf3d
    integer              :: indcld          !< location of cloud fraction in phyf3d (used only for SHOC or MG)
    integer              :: num_p2d         !< number of 2D arrays needed for microphysics
    integer              :: num_p3d         !< number of 3D arrays needed for microphysics
    integer              :: nshoc_2d        !< number of 2d fields for SHOC
    integer              :: nshoc_3d        !< number of 3d fields for SHOC
    integer              :: ncnvcld3d       !< number of convective 3d clouds fields
    integer              :: npdf3d          !< number of 3d arrays associated with pdf based clouds/microphysics
    integer              :: nctp            !< number of cloud types in Chikira-Sugiyama scheme
    integer              :: ncnvw           !< the index of cnvw in phy_f3d
    integer              :: ncnvc           !< the index of cnvc in phy_f3d
    integer              :: nleffr          !< the index of cloud liquid water effective radius in phy_f3d
    integer              :: nieffr          !< the index of ice effective radius in phy_f3d
    integer              :: nreffr          !< the index of rain effective radius in phy_f3d
    integer              :: nseffr          !< the index of snow effective radius in phy_f3d
    integer              :: ngeffr          !< the index of graupel effective radius in phy_f3d
    integer              :: nkbfshoc        !< the index of upward kinematic buoyancy flux from SHOC in phy_f3d
    integer              :: nahdshoc        !< the index of diffusivity for heat from from SHOC in phy_f3d
    integer              :: nscfshoc        !< the index of subgrid-scale cloud fraction from from SHOC in phy_f3d
    integer              :: nT2delt         !< the index of air temperature 2 timesteps back for Z-C MP in phy_f3d
    integer              :: nTdelt          !< the index of air temperature at the previous timestep for Z-C MP in phy_f3d
    integer              :: nqv2delt        !< the index of specific humidity 2 timesteps back for Z-C MP in phy_f3d
    integer              :: nqvdelt         !< the index of specific humidity at the previous timestep for Z-C MP in phy_f3d
    integer              :: nps2delt        !< the index of surface air pressure 2 timesteps back for Z-C MP in phy_f2d
    integer              :: npsdelt         !< the index of surface air pressure at the previous timestep for Z-C MP in phy_f2d
    integer              :: ncnvwind        !< the index of surface wind enhancement due to convection for MYNN SFC and RAS CNV in phy f2d

!-- nml variables for RRFS-SD
    real(kind=kind_phys) :: dust_drylimit_factor  !< factor for drylimit parameterization in fengsha
    real(kind=kind_phys) :: dust_moist_correction !< factor to tune volumetric soil moisture
    integer              :: dust_moist_opt        !< dust moisture option 1:fecan 2:shao
    real(kind=kind_phys) :: dust_alpha        !< alpha parameter for fengsha dust scheme
    real(kind=kind_phys) :: dust_gamma        !< gamma parameter for fengsha dust scheme
    real(kind=kind_phys) :: wetdep_ls_alpha   !< alpha parameter for wet deposition
    integer              :: ebb_dcycle        !< 1:retro; 2:forecast of fire emission
    integer              :: seas_opt
    integer              :: dust_opt
    integer              :: drydep_opt
    integer              :: coarsepm_settling
    integer              :: plume_wind_eff
    logical              :: extended_sd_diags
    integer              :: wetdep_ls_opt
    logical              :: do_plumerise
    integer              :: addsmoke_flag
    integer              :: plumerisefire_frq
    integer              :: n_dbg_lines
    integer              :: smoke_forecast
    logical              :: aero_ind_fdb    ! WFA/IFA indirect
    logical              :: aero_dir_fdb    ! smoke/dust direct
    logical              :: rrfs_smoke_debug
    logical              :: do_smoke_transport
    logical              :: mix_chem
    logical              :: enh_mix
    real(kind=kind_phys) :: smoke_dir_fdb_coef(7) !< smoke & dust direct feedbck coefficents
    real(kind=kind_phys) :: smoke_conv_wet_coef(3) !< smoke & dust convective wet removal coefficents

!--- debug flags
    logical              :: debug
    logical              :: pre_rad         !< flag for testing purpose
    logical              :: print_diff_pgr  !< print average change in pgr every timestep (does not need debug flag)

!--- variables modified at each time step
    integer              :: ipt             !< index for diagnostic printout point
    logical              :: lprnt           !< control flag for diagnostic print out
    logical              :: lsswr           !< logical flags for sw radiation calls
    logical              :: lslwr           !< logical flags for lw radiation calls
    real(kind=kind_phys) :: solhr           !< hour time after 00z at the t-step
    real(kind=kind_phys) :: solcon          !< solar constant (sun-earth distant adjusted)  [set via radupdate]
    real(kind=kind_phys) :: slag            !< equation of time ( radian )                  [set via radupdate]
    real(kind=kind_phys) :: sdec            !< sin of the solar declination angle           [set via radupdate]
    real(kind=kind_phys) :: cdec            !< cos of the solar declination angle           [set via radupdate]
    real(kind=kind_phys) :: clstp           !< index used by cnvc90 (for convective clouds)
                                            !< legacy stuff - does not affect forecast
    real(kind=kind_phys) :: phour           !< previous forecast hour
    real(kind=kind_phys) :: fhour           !< current forecast hour
    real(kind=kind_phys) :: zhour           !< previous hour diagnostic buckets emptied
    integer              :: kdt             !< current forecast iteration
    logical              :: first_time_step !< flag signaling first time step for time integration routine
    logical              :: restart         !< flag whether this is a coldstart (.false.) or a warmstart/restart (.true.)
    logical              :: lsm_cold_start
    logical              :: hydrostatic     !< flag whether this is a hydrostatic or non-hydrostatic run
    integer              :: jdat(1:8)       !< current forecast date and time
                                            !< (yr, mon, day, t-zone, hr, min, sec, mil-sec)
    integer              :: imn             !< initial forecast month
    real(kind=kind_phys) :: julian          !< julian day using midnight of January 1 of forecast year as initial epoch
    integer              :: yearlen         !< length of the current forecast year in days
!
    integer              :: iccn            !< using IN CCN forcing for MG2/3
    real(kind=kind_phys), pointer :: si(:)  !< vertical sigma coordinate for model initialization
    real(kind=kind_phys)          :: sec    !< seconds since model initialization

!--- Increment grid    
    logical              :: increment_file_on_native_grid ! increment on native grid else Gaussian grid
    
!--- IAU
    integer              :: iau_offset
    real(kind=kind_phys) :: iau_delthrs     ! iau time interval (to scale increments) in hours
    character(len=240)   :: iau_inc_files(7)! list of increment files
    real(kind=kind_phys) :: iaufhrs(7)      ! forecast hours associated with increment files
    logical :: iau_filter_increments, iau_drymassfixer

    ! From physcons.F90, updated/set in control_initialize
    real(kind=kind_phys) :: dxinv           ! inverse scaling factor for critical relative humidity, replaces dxinv in physcons.F90
    real(kind=kind_phys) :: dxmax           ! maximum scaling factor for critical relative humidity, replaces dxmax in physcons.F90
    real(kind=kind_phys) :: dxmin           ! minimum scaling factor for critical relative humidity, replaces dxmin in physcons.F90
    real(kind=kind_phys) :: rhcmax          ! maximum critical relative humidity, replaces rhc_max in physcons.F90
    real(kind=kind_phys) :: huge            !< huge fill value

!--- lightning threat and diagsnostics
    logical              :: lightning_threat !< report lightning threat indices

!--- NRL Ozone physics
    logical         :: oz_phys         !< Flag for old (2006) ozone physics
    logical         :: oz_phys_2015    !< Flag for new (2015) ozone physics
    type(ty_ozphys) :: ozphys          !< DDT with data needed by ozone physics
    integer         :: levozp          !< Number of vertical layers in ozone forcing data
    integer         :: oz_coeff        !< Number of coefficients in ozone forcing data
!--- NRL h2o photchemistry physics
    type(ty_h2ophys) :: h2ophys        !< DDT with data needed by h2o photchemistry physics.
    integer          :: levh2o         !< Number of vertical layers in stratospheric h2o data.
    integer          :: h2o_coeff      !< Number of coefficients in stratospheric h2o data.
    
    contains
      procedure :: init            => control_initialize
      procedure :: init_chemistry  => control_chemistry_initialize
      procedure :: init_scavenging => control_scavenging_initialize
      procedure :: print           => control_print
  end type GFS_control_type


!--------------------------------------------------------------------
! GFS_grid_type
!   grid data needed for interpolations and length-scale calculations
!--------------------------------------------------------------------
!! \section arg_table_GFS_grid_type
!! \htmlinclude GFS_grid_type.html
!!
  type GFS_grid_type

    real (kind=kind_phys), pointer :: xlon   (:)        !< grid longitude in radians, ok for both 0->2pi
                                                        !! or -pi -> +pi ranges
    real (kind=kind_phys), pointer :: xlat   (:)        !< grid latitude in radians, default to pi/2 ->
                                                        !! -pi/2 range, otherwise adj in subr called
    real (kind=kind_phys), pointer :: xlat_d (:)        !< grid latitude in degrees, default to 90 ->
                                                        !! -90 range, otherwise adj in subr called
    real (kind=kind_phys), pointer :: xlon_d (:)        !< grid longitude in degrees, default to 0 ->
                                                        !! 360 range, otherwise adj in subr called
    real (kind=kind_phys), pointer :: sinlat (:)        !< sine of the grids corresponding latitudes
    real (kind=kind_phys), pointer :: coslat (:)        !< cosine of the grids corresponding latitudes
    real (kind=kind_phys), pointer :: area   (:)        !< area of the grid cell
    real (kind=kind_phys), pointer :: dx     (:)        !< relative dx for the grid cell

!--- grid-related interpolation data for prognostic ozone
    real (kind=kind_phys), pointer :: ddy_o3    (:)     !< interpolation     weight for ozone
    integer,               pointer :: jindx1_o3 (:)     !< interpolation  low index for ozone
    integer,               pointer :: jindx2_o3 (:)     !< interpolation high index for ozone

!--- grid-related interpolation data for stratosphere water
    real (kind=kind_phys), pointer :: ddy_h     (:)     !< interpolation     weight for h2o
    integer,               pointer :: jindx1_h  (:)     !< interpolation  low index for h2o
    integer,               pointer :: jindx2_h  (:)     !< interpolation high index for h2o

!--- grid-related interpolation data for prognostic iccn
    real (kind=kind_phys), pointer :: ddy_ci    (:)     !< interpolation     weight for iccn
    integer,               pointer :: jindx1_ci (:)     !< interpolation  low index for iccn
    integer,               pointer :: jindx2_ci (:)     !< interpolation high index for iccn
    real (kind=kind_phys), pointer :: ddx_ci    (:)     !< interpolation     weight for iccn
    integer,               pointer :: iindx1_ci (:)     !< interpolation  low index for iccn
    integer,               pointer :: iindx2_ci (:)     !< interpolation high index for iccn

!--- grid-related interpolation data for prescribed aerosols
    real (kind=kind_phys), pointer :: ddy_aer    (:)    !< interpolation     weight for iaerclm
    integer,               pointer :: jindx1_aer (:)    !< interpolation  low index for iaerclm
    integer,               pointer :: jindx2_aer (:)    !< interpolation high index for iaerclm
    real (kind=kind_phys), pointer :: ddx_aer    (:)    !< interpolation     weight for iaerclm
    integer,               pointer :: iindx1_aer (:)    !< interpolation  low index for iaerclm
    integer,               pointer :: iindx2_aer (:)    !< interpolation high index for iaerclm

!--- grid-related interpolation data for cires_ugwp_v1
    real (kind=kind_phys), pointer :: ddy_j1tau  (:)    !< interpolation     weight for  tau_ugwp
    real (kind=kind_phys), pointer :: ddy_j2tau  (:)    !< interpolation     weight for  tau_ugwp
    integer,               pointer :: jindx1_tau (:)    !< interpolation  low index for tau_ugwp
    integer,               pointer :: jindx2_tau (:)    !< interpolation high index for tau_ugwp

    contains
      procedure :: create   => grid_create   !<   allocate array data
  end type GFS_grid_type


!-----------------------------------------------
! GFS_tbd_type
!   data not yet assigned to a defined container
!-----------------------------------------------
!! \section arg_table_GFS_tbd_type
!! \htmlinclude GFS_tbd_type.html
!!
  type GFS_tbd_type

!--- radiation random seeds
    integer,               pointer :: icsdsw   (:)     => null()  !< (rad. only) auxiliary cloud control arrays passed to main
    integer,               pointer :: icsdlw   (:)     => null()  !< (rad. only) radiations. if isubcsw/isubclw (input to init)
                                                                  !< (rad. only) are set to 2, the arrays contains provided
                                                                  !< (rad. only) random seeds for sub-column clouds generators
    integer,               pointer :: rseeds   (:,:)   => null()  !< (rad. only) random seeds provided by host

!--- In
    real (kind=kind_phys), pointer :: ozpl     (:,:,:) => null()  !< ozone forcing data
    real (kind=kind_phys), pointer :: h2opl    (:,:,:) => null()  !< water forcing data
    real (kind=kind_phys), pointer :: in_nm    (:,:)   => null()  !< IN number concentration
    real (kind=kind_phys), pointer :: ccn_nm   (:,:)   => null()  !< CCN number concentration
    real (kind=kind_phys), pointer :: aer_nm   (:,:,:) => null()  !< GOCART aerosol climo
    real (kind=kind_phys), pointer :: tau_amf  (:    ) => null()  !< nonsta-gw monthly data

    integer,               pointer :: imap     (:)     => null()  !< map of local index ix to global index i for this block
    integer,               pointer :: jmap     (:)     => null()  !< map of local index ix to global index j for this block

    !--- active when ((.not. newsas .or. cal_pre) .and. random_clds)
    real (kind=kind_phys), pointer :: rann     (:,:)   => null()  !< random number array (0-1)

!--- In/Out
    real (kind=kind_phys), pointer :: acv      (:)     => null()  !< array containing accumulated convective clouds
    real (kind=kind_phys), pointer :: acvb     (:)     => null()  !< arrays used by cnvc90 bottom
    real (kind=kind_phys), pointer :: acvt     (:)     => null()  !< arrays used by cnvc90 top (cnvc90.f)

!--- Stochastic physics properties calculated in physics_driver
    real (kind=kind_phys), pointer :: dtdtnp    (:,:)  => null()  !< temperature change from physics that should not be perturbed with SPPT (k)
    real (kind=kind_phys), pointer :: dtotprcp  (:)    => null()  !< change in totprcp  (diag_type)
    real (kind=kind_phys), pointer :: dcnvprcp  (:)    => null()  !< change in cnvprcp  (diag_type)
    real (kind=kind_phys), pointer :: drain_cpl (:)    => null()  !< change in rain_cpl (coupling_type)
    real (kind=kind_phys), pointer :: dsnow_cpl (:)    => null()  !< change in show_cpl (coupling_type)

!--- phy_f*d variables needed for seamless restarts and moving data between grrad and gbphys
    real (kind=kind_phys), pointer :: phy_fctd (:,:)   => null()  !< cloud base mass flux for CS convection
    real (kind=kind_phys), pointer :: phy_f2d  (:,:)   => null()  !< 2d arrays saved for restart
    real (kind=kind_phys), pointer :: phy_f3d  (:,:,:) => null()  !< 3d arrays saved for restart

!--- Diagnostic that needs to be carried over to the next time step (removed from diag_type)
    real (kind=kind_phys), pointer :: hpbl     (:)     => null()  !< Planetary boundary layer height
    real (kind=kind_phys), pointer :: ud_mf  (:,:)     => null()  !< updraft mass flux

    !--- dynamical forcing variables for Grell-Freitas convection
    real (kind=kind_phys), pointer :: forcet (:,:)     => null()  !<
    real (kind=kind_phys), pointer :: forceq (:,:)     => null()  !<
    real (kind=kind_phys), pointer :: prevst (:,:)     => null()  !<
    real (kind=kind_phys), pointer :: prevsq (:,:)     => null()  !<
    integer,               pointer :: cactiv   (:)     => null()  !< convective activity memory contour
    integer,               pointer :: cactiv_m (:)     => null()  !< mid-level convective activity memory contour
    real (kind=kind_phys), pointer :: aod_gf   (:)     => null()

    !--- MYNN prognostic variables that can't be in the Intdiag or Interstitial DDTs
    real (kind=kind_phys), pointer :: CLDFRA_BL  (:,:)   => null()  !
    real (kind=kind_phys), pointer :: QC_BL      (:,:)   => null()  !
    real (kind=kind_phys), pointer :: QI_BL      (:,:)   => null()  !
    real (kind=kind_phys), pointer :: el_pbl     (:,:)   => null()  !
    real (kind=kind_phys), pointer :: Sh3D       (:,:)   => null()  !
    real (kind=kind_phys), pointer :: Sm3D       (:,:)   => null()  !
    real (kind=kind_phys), pointer :: qke        (:,:)   => null()  !
    real (kind=kind_phys), pointer :: tsq        (:,:)   => null()  !
    real (kind=kind_phys), pointer :: qsq        (:,:)   => null()  !
    real (kind=kind_phys), pointer :: cov        (:,:)   => null()  !

    !--- MYJ schemes saved variables (from previous time step)
    real (kind=kind_phys), pointer :: phy_myj_qsfc(:)    => null()  !
    real (kind=kind_phys), pointer :: phy_myj_thz0(:)    => null()  !
    real (kind=kind_phys), pointer :: phy_myj_qz0(:)     => null()  !
    real (kind=kind_phys), pointer :: phy_myj_uz0(:)     => null()  !
    real (kind=kind_phys), pointer :: phy_myj_vz0(:)     => null()  !
    real (kind=kind_phys), pointer :: phy_myj_akhs(:)    => null()  !
    real (kind=kind_phys), pointer :: phy_myj_akms(:)    => null()  !
    real (kind=kind_phys), pointer :: phy_myj_chkqlm(:)  => null()  !
    real (kind=kind_phys), pointer :: phy_myj_elflx(:)   => null()  !
    real (kind=kind_phys), pointer :: phy_myj_a1u(:)     => null()  !
    real (kind=kind_phys), pointer :: phy_myj_a1t(:)     => null()  !
    real (kind=kind_phys), pointer :: phy_myj_a1q(:)     => null()  !

    !--- DFI Radar
    real (kind=kind_phys), pointer :: dfi_radar_tten(:,:,:) => null() !
    real (kind=kind_phys), pointer :: cap_suppress(:,:) => null() !

    contains
      procedure :: create  => tbd_create  !<   allocate array data
  end type GFS_tbd_type


!------------------------------------------------------------------
! GFS_cldprop_type
!  cloud properties and tendencies needed by radiation from physics
!------------------------------------------------------------------
!! \section arg_table_GFS_cldprop_type
!! \htmlinclude GFS_cldprop_type.html
!!
  type GFS_cldprop_type

!--- In     (radiation)
!--- In/Out (physics)
    real (kind=kind_phys), pointer :: cv  (:)     => null()  !< fraction of convective cloud ; phys
    real (kind=kind_phys), pointer :: cvt (:)     => null()  !< convective cloud top pressure in pa ; phys
    real (kind=kind_phys), pointer :: cvb (:)     => null()  !< convective cloud bottom pressure in pa ; phys, cnvc90

    contains
      procedure :: create  => cldprop_create  !<   allocate array data
  end type GFS_cldprop_type


!-----------------------------------------
! GFS_radtend_type
!   radiation tendencies needed by physics
!-----------------------------------------
!! \section arg_table_GFS_radtend_type
!! \htmlinclude GFS_radtend_type.html
!!
  type GFS_radtend_type

    type (sfcfsw_type),    pointer :: sfcfsw(:)   => null()   !< sw radiation fluxes at sfc
                                                              !< [dim(im): created in grrad.f], components:
                                                              !!     (check module_radsw_parameters for definition)
                                                              !!\n   %upfxc - total sky upward sw flux at sfc (w/m**2)
                                                              !!\n   %upfx0 - clear sky upward sw flux at sfc (w/m**2)
                                                              !!\n   %dnfxc - total sky downward sw flux at sfc (w/m**2)
                                                              !!\n   %dnfx0 - clear sky downward sw flux at sfc (w/m**2)

    type (sfcflw_type),    pointer :: sfcflw(:)    => null()  !< lw radiation fluxes at sfc
                                                              !< [dim(im): created in grrad.f], components:
                                                              !!     (check module_radlw_paramters for definition)
                                                              !!\n   %upfxc - total sky upward lw flux at sfc (w/m**2)
                                                              !!\n   %upfx0 - clear sky upward lw flux at sfc (w/m**2)
                                                              !!\n   %dnfxc - total sky downward lw flux at sfc (w/m**2)
                                                              !!\n   %dnfx0 - clear sky downward lw flux at sfc (w/m**2)

!--- Out (radiation only)
    real (kind=kind_phys), pointer :: htrsw (:,:)  => null()  !< swh  total sky sw heating rate in k/sec
    real (kind=kind_phys), pointer :: htrlw (:,:)  => null()  !< hlw  total sky lw heating rate in k/sec
    real (kind=kind_phys), pointer :: sfalb (:)    => null()  !< mean surface diffused sw albedo

    real (kind=kind_phys), pointer :: coszen(:)    => null()  !< mean cos of zenith angle over rad call period
    real (kind=kind_phys), pointer :: tsflw (:)    => null()  !< surface air temp during lw calculation in k
    real (kind=kind_phys), pointer :: semis (:)    => null()  !< surface lw emissivity in fraction
    real (kind=kind_phys), pointer :: ext550 (:,:) => null()  !< aerosol optical extinction from radiation

!--- In/Out (???) (radiaition only)
    real (kind=kind_phys), pointer :: coszdg(:)    => null()  !< daytime mean cosz over rad call period

!--- In/Out (???) (physics only)
    real (kind=kind_phys), pointer :: swhc (:,:)   => null()  !< clear sky sw heating rates ( k/s )
    real (kind=kind_phys), pointer :: lwhc (:,:)   => null()  !< clear sky lw heating rates ( k/s )
    real (kind=kind_phys), pointer :: lwhd (:,:,:) => null()  !< idea sky lw heating rates ( k/s ) !DJS2023 THIS IS NOT USED. IT IS REFERENCED, BUT NEVER SET?

    contains
      procedure :: create  => radtend_create   !<   allocate array data
  end type GFS_radtend_type


!----------------------------------------------------------------
! GFS_diag_type
!  internal diagnostic type used as arguments to gbphys and grrad
!----------------------------------------------------------------
!! \section arg_table_GFS_diag_type
!! \htmlinclude GFS_diag_type.html
!!
  type GFS_diag_type

!! Input/Output only in radiation
    real (kind=kind_phys), pointer :: fluxr(:,:)     => null()   !< to save time accumulated 2-d fields defined as:!
                                                                 !< hardcoded field indices, opt. includes aerosols!
    type (topfsw_type),    pointer :: topfsw(:)      => null()   !< sw radiation fluxes at toa, components:
                                               !       %upfxc    - total sky upward sw flux at toa (w/m**2)
                                               !       %dnfxc    - total sky downward sw flux at toa (w/m**2)
                                               !       %upfx0    - clear sky upward sw flux at toa (w/m**2)
    type (topflw_type),    pointer :: topflw(:)      => null()   !< lw radiation fluxes at top, component:
                                               !        %upfxc    - total sky upward lw flux at toa (w/m**2)
                                               !        %upfx0    - clear sky upward lw flux at toa (w/m**2)

! Input/output - used by physics
    real (kind=kind_phys), pointer :: srunoff(:)     => null()   !< accumulated surface storm runoff (from lsm)
    real (kind=kind_phys), pointer :: evbsa  (:)     => null()   !< accumulated direct evaporation
    real (kind=kind_phys), pointer :: evcwa  (:)     => null()   !< accumulated canopy evaporation
    real (kind=kind_phys), pointer :: snohfa (:)     => null()   !< heat flux for phase change of snow (melting)
    real (kind=kind_phys), pointer :: transa (:)     => null()   !< accumulated transpiration
    real (kind=kind_phys), pointer :: sbsnoa (:)     => null()   !< accumulated snow sublimation
    real (kind=kind_phys), pointer :: snowca (:)     => null()   !< snow cover
    real (kind=kind_phys), pointer :: sbsno  (:)     => null()   !< instantaneous snow sublimation
    real (kind=kind_phys), pointer :: evbs(:)        => null()   !< instantaneous direct evaporation
    real (kind=kind_phys), pointer :: trans  (:)     => null()   !< instantaneous transpiration
    real (kind=kind_phys), pointer :: evcw(:)        => null()   !< instantaneous canopy evaporation
    real (kind=kind_phys), pointer :: snowmt_land(:) => null()   !< ruc lsm diagnostics over land
    real (kind=kind_phys), pointer :: snowmt_ice(:)  => null()   !< ruc lsm diagnostics over ice
    real (kind=kind_phys), pointer :: soilm  (:)     => null()   !< integrated soil moisture
    real (kind=kind_phys), pointer :: paha   (:)     => null()   !< noah lsm diagnostics
    real (kind=kind_phys), pointer :: tmpmin (:)     => null()   !< min temperature at 2m height (k)
    real (kind=kind_phys), pointer :: tmpmax (:)     => null()   !< max temperature at 2m height (k)
    real (kind=kind_phys), pointer :: dusfc  (:)     => null()   !< u component of surface stress
    real (kind=kind_phys), pointer :: dvsfc  (:)     => null()   !< v component of surface stress
    real (kind=kind_phys), pointer :: dtsfc  (:)     => null()   !< sensible heat flux (w/m2)
    real (kind=kind_phys), pointer :: dqsfc  (:)     => null()   !< latent heat flux (w/m2)
    real (kind=kind_phys), pointer :: totprcp(:)     => null()   !< accumulated total precipitation (kg/m2)
    real (kind=kind_phys), pointer :: totprcpb(:)    => null()   !< accumulated total precipitation in bucket(kg/m2)
    real (kind=kind_phys), pointer :: gflux  (:)     => null()   !< groud conductive heat flux
    real (kind=kind_phys), pointer :: dlwsfc (:)     => null()   !< time accumulated sfc dn lw flux ( w/m**2 )
    real (kind=kind_phys), pointer :: ulwsfc (:)     => null()   !< time accumulated sfc up lw flux ( w/m**2 )
    real (kind=kind_phys), pointer :: suntim (:)     => null()   !< sunshine duration time (s)
    real (kind=kind_phys), pointer :: runoff (:)     => null()   !< total water runoff
    real (kind=kind_phys), pointer :: ep     (:)     => null()   !< potential evaporation
    real (kind=kind_phys), pointer :: tecan  (:)     => null()   !< total evaporation of intercepted water
    real (kind=kind_phys), pointer :: tetran (:)     => null()   !< total transpiration rate
    real (kind=kind_phys), pointer :: tedir  (:)     => null()   !< total soil surface evaporation rate
    real (kind=kind_phys), pointer :: twa    (:)     => null()   !< total water storage in aquifer
    real (kind=kind_phys), pointer :: cldwrk (:)     => null()   !< cloud workfunction (valid only with sas)
    real (kind=kind_phys), pointer :: dugwd  (:)     => null()   !< vertically integrated u change by OGWD
    real (kind=kind_phys), pointer :: dvgwd  (:)     => null()   !< vertically integrated v change by OGWD
    real (kind=kind_phys), pointer :: psmean (:)     => null()   !< surface pressure (kPa)
    real (kind=kind_phys), pointer :: cnvprcp(:)     => null()   !< accumulated convective precipitation (kg/m2)
    real (kind=kind_phys), pointer :: cnvprcpb(:)    => null()   !< accumulated convective precipitation in bucket (kg/m2)
    real (kind=kind_phys), pointer :: spfhmin(:)     => null()   !< minimum specific humidity
    real (kind=kind_phys), pointer :: spfhmax(:)     => null()   !< maximum specific humidity
    real (kind=kind_phys), pointer :: u10mmax(:)     => null()   !< maximum u-wind
    real (kind=kind_phys), pointer :: v10mmax(:)     => null()   !< maximum v-wind
    real (kind=kind_phys), pointer :: wind10mmax(:)  => null()   !< maximum wind speed
    real (kind=kind_phys), pointer :: u10max(:)      => null()   !< maximum u-wind used with avg_max_length
    real (kind=kind_phys), pointer :: v10max(:)      => null()   !< maximum v-wind used with avg_max_length
    real (kind=kind_phys), pointer :: spd10max(:)    => null()   !< maximum wind speed used with avg_max_length
    real (kind=kind_phys), pointer :: rain   (:)     => null()   !< total rain at this time step
    real (kind=kind_phys), pointer :: rainc  (:)     => null()   !< convective rain at this time step
    real (kind=kind_phys), pointer :: ice    (:)     => null()   !< ice fall at this time step
    real (kind=kind_phys), pointer :: snow   (:)     => null()   !< snow fall at this time step
    real (kind=kind_phys), pointer :: graupel(:)     => null()   !< graupel fall at this time step
    real (kind=kind_phys), pointer :: totice (:)     => null()   !< accumulated ice precipitation (kg/m2)
    real (kind=kind_phys), pointer :: totsnw (:)     => null()   !< accumulated snow precipitation (kg/m2)
    real (kind=kind_phys), pointer :: totgrp (:)     => null()   !< accumulated graupel precipitation (kg/m2)
    real (kind=kind_phys), pointer :: toticeb(:)     => null()   !< accumulated ice precipitation in bucket (kg/m2)
    real (kind=kind_phys), pointer :: totsnwb(:)     => null()   !< accumulated snow precipitation in bucket (kg/m2)
    real (kind=kind_phys), pointer :: totgrpb(:)     => null()   !< accumulated graupel precipitation in bucket (kg/m2)
    real (kind=kind_phys), pointer :: frzr   (:)     => null()   !< accumulated surface freezing rain (m)
    real (kind=kind_phys), pointer :: frzrb  (:)     => null()   !< accumulated surface freezing rain in bucket (m)
    real (kind=kind_phys), pointer :: frozr  (:)     => null()   !< accumulated surface graupel (m)
    real (kind=kind_phys), pointer :: frozrb (:)     => null()   !< accumulated surface graupel in bucket (m)
    real (kind=kind_phys), pointer :: tsnowp (:)     => null()   !< accumulated surface snowfall (m)
    real (kind=kind_phys), pointer :: tsnowpb(:)     => null()   !< accumulated surface snowfall in bucket (m)
    real (kind=kind_phys), pointer :: rhonewsn1(:)   => null()   !< precipitation ice density outside RUC LSM (kg/m3)

    !--- MYNN variables
    real (kind=kind_phys), pointer :: edmf_a     (:,:)   => null()  !
    real (kind=kind_phys), pointer :: edmf_w     (:,:)   => null()  !
    real (kind=kind_phys), pointer :: edmf_qt    (:,:)   => null()  !
    real (kind=kind_phys), pointer :: edmf_thl   (:,:)   => null()  !
    real (kind=kind_phys), pointer :: edmf_ent   (:,:)   => null()  !
    real (kind=kind_phys), pointer :: edmf_qc    (:,:)   => null()  !
    real (kind=kind_phys), pointer :: sub_thl    (:,:)   => null()  !
    real (kind=kind_phys), pointer :: sub_sqv    (:,:)   => null()  !
    real (kind=kind_phys), pointer :: det_thl    (:,:)   => null()  !
    real (kind=kind_phys), pointer :: det_sqv    (:,:)   => null()  !
    real (kind=kind_phys), pointer :: maxMF       (:)    => null()  !
    real (kind=kind_phys), pointer :: maxwidth    (:)    => null()  !
    real (kind=kind_phys), pointer :: ztop_plume  (:)    => null()  !
    integer, pointer               :: ktop_plume  (:)    => null()  !
    real (kind=kind_phys), pointer :: exch_h     (:,:)   => null()  !
    real (kind=kind_phys), pointer :: exch_m     (:,:)   => null()  !
    real (kind=kind_phys), pointer :: dqke       (:,:)   => null()  !< timestep change of tke
    real (kind=kind_phys), pointer :: qwt        (:,:)   => null()  !< vertical transport of tke
    real (kind=kind_phys), pointer :: qshear     (:,:)   => null()  !< shear production of tke
    real (kind=kind_phys), pointer :: qbuoy      (:,:)   => null()  !< buoyancy production of tke
    real (kind=kind_phys), pointer :: qdiss      (:,:)   => null()  !< dissipation of tke

! Output - only in physics
    real (kind=kind_phys), pointer :: u10m   (:)     => null()   !< 10 meter u/v wind speed
    real (kind=kind_phys), pointer :: v10m   (:)     => null()   !< 10 meter u/v wind speed
    real (kind=kind_phys), pointer :: dpt2m  (:)     => null()   !< 2 meter dew point temperature
    real (kind=kind_phys), pointer :: zlvl   (:)     => null()   !< layer 1 height (m)
    real (kind=kind_phys), pointer :: psurf  (:)     => null()   !< surface pressure (Pa)
    real (kind=kind_phys), pointer :: pwat   (:)     => null()   !< precipitable water
    real (kind=kind_phys), pointer :: t1     (:)     => null()   !< layer 1 temperature (K)
    real (kind=kind_phys), pointer :: q1     (:)     => null()   !< layer 1 specific humidity (kg/kg)
    real (kind=kind_phys), pointer :: u1     (:)     => null()   !< layer 1 zonal wind (m/s)
    real (kind=kind_phys), pointer :: v1     (:)     => null()   !< layer 1 merdional wind (m/s)
    real (kind=kind_phys), pointer :: chh    (:)     => null()   !< thermal exchange coefficient
    real (kind=kind_phys), pointer :: cmm    (:)     => null()   !< momentum exchange coefficient
    real (kind=kind_phys), pointer :: dlwsfci(:)     => null()   !< instantaneous sfc dnwd lw flux ( w/m**2 )
    real (kind=kind_phys), pointer :: ulwsfci(:)     => null()   !< instantaneous sfc upwd lw flux ( w/m**2 )
    real (kind=kind_phys), pointer :: dswsfci(:)     => null()   !< instantaneous sfc dnwd sw flux ( w/m**2 )
    real (kind=kind_phys), pointer :: nswsfci(:)     => null()   !< instantaneous sfc net dnwd sw flux ( w/m**2 )
    real (kind=kind_phys), pointer :: uswsfci(:)     => null()   !< instantaneous sfc upwd sw flux ( w/m**2 )
    real (kind=kind_phys), pointer :: dusfci (:)     => null()   !< instantaneous u component of surface stress
    real (kind=kind_phys), pointer :: dvsfci (:)     => null()   !< instantaneous v component of surface stress
    real (kind=kind_phys), pointer :: dtsfci (:)     => null()   !< instantaneous sfc sensible heat flux
    real (kind=kind_phys), pointer :: dqsfci (:)     => null()   !< instantaneous sfc latent heat flux
    real (kind=kind_phys), pointer :: gfluxi (:)     => null()   !< instantaneous sfc ground heat flux
    real (kind=kind_phys), pointer :: pahi   (:)     => null()   !< instantaneous precipitation advected heat flux
    real (kind=kind_phys), pointer :: epi    (:)     => null()   !< instantaneous sfc potential evaporation
    real (kind=kind_phys), pointer :: smcwlt2(:)     => null()   !< wilting point (volumetric)
    real (kind=kind_phys), pointer :: smcref2(:)     => null()   !< soil moisture threshold (volumetric)
    real (kind=kind_phys), pointer :: wet1   (:)     => null()   !< normalized soil wetness
    real (kind=kind_phys), pointer :: sr     (:)     => null()   !< snow ratio : ratio of snow to total precipitation
    real (kind=kind_phys), pointer :: tdomr  (:)     => null()   !< dominant accumulated rain type
    real (kind=kind_phys), pointer :: tdomzr (:)     => null()   !< dominant accumulated freezing rain type
    real (kind=kind_phys), pointer :: tdomip (:)     => null()   !< dominant accumulated sleet type
    real (kind=kind_phys), pointer :: tdoms  (:)     => null()   !< dominant accumulated snow type
    real (kind=kind_phys), pointer :: zmtnblck(:)    => null()   !<mountain blocking level

    ! dtend/dtidxt: Multitudinous 3d tendencies in a 4D array: (i,k,1:100+ntrac,nprocess)
    ! Sparse in outermost two dimensions. dtidx(1:100+ntrac,nprocess) maps to dtend
    ! outer dimension index.
    real (kind=kind_phys), pointer :: dtend (:,:,:)  => null()   !< tracer changes due to physics

    real (kind=kind_phys), pointer :: refdmax (:)    => null()   !< max hourly 1-km agl reflectivity
    real (kind=kind_phys), pointer :: refdmax263k(:) => null()   !< max hourly -10C reflectivity
    real (kind=kind_phys), pointer :: t02max  (:)    => null()   !< max hourly 2m T
    real (kind=kind_phys), pointer :: t02min  (:)    => null()   !< min hourly 2m T
    real (kind=kind_phys), pointer :: rh02max (:)    => null()   !< max hourly 2m RH
    real (kind=kind_phys), pointer :: rh02min (:)    => null()   !< min hourly 2m RH
    real (kind=kind_phys), pointer :: pratemax(:)    => null()   !< max hourly precipitation rate
!--- accumulated quantities for 3D diagnostics
    real (kind=kind_phys), pointer :: upd_mf (:,:)   => null()  !< instantaneous convective updraft mass flux
    real (kind=kind_phys), pointer :: dwn_mf (:,:)   => null()  !< instantaneous convective downdraft mass flux
    real (kind=kind_phys), pointer :: det_mf (:,:)   => null()  !< instantaneous convective detrainment mass flux
!--- F-A MP scheme
    real (kind=kind_phys), pointer :: train  (:,:)   => null()  !< accumulated stratiform T tendency (K s-1)
    real (kind=kind_phys), pointer :: cldfra (:,:)   => null()  !< instantaneous 3D cloud fraction
    !--- MP quantities for 3D diagnositics
    real (kind=kind_phys), pointer :: refl_10cm(:,:) => null()  !< instantaneous refl_10cm
    real (kind=kind_phys), pointer :: max_hail_diam_sfc(:) => null()  !< instantaneous max hail diameter sfc
    real (kind=kind_phys), pointer :: cldfra2d (:)   => null()  !< instantaneous 2D cloud fraction
    real (kind=kind_phys), pointer :: total_albedo (:)   => null()  !< total sky (with cloud) albedo at toa
    real (kind=kind_phys), pointer :: lwp_ex (:)     => null()  !< liquid water path from microphysics
    real (kind=kind_phys), pointer :: iwp_ex (:)     => null()  !< ice water path from microphysics
    real (kind=kind_phys), pointer :: lwp_fc (:)     => null()  !< liquid water path from cloud fraction scheme
    real (kind=kind_phys), pointer :: iwp_fc (:)     => null()  !< ice water path from cloud fraction scheme

    !--- Extra PBL diagnostics
    real (kind=kind_phys), pointer :: dkt(:,:)       => null()  !< Eddy diffusitivity for heat
    real (kind=kind_phys), pointer :: dku(:,:)       => null()  !< Eddy diffusitivity for momentum

!
!---vay-2018 UGWP-diagnostics instantaneous
!
! OGWs +NGWs
    real (kind=kind_phys), pointer :: dudt_gw(:,:)   => null()  !<
    real (kind=kind_phys), pointer :: dvdt_gw(:,:)   => null()  !<
    real (kind=kind_phys), pointer :: dtdt_gw(:,:)   => null()  !<
    real (kind=kind_phys), pointer :: kdis_gw(:,:)   => null()  !<
!oro-GWs
    real (kind=kind_phys), pointer :: dudt_ogw(:,:)  => null()  !<
    real (kind=kind_phys), pointer :: dvdt_ogw(:,:)  => null()  !<
    real (kind=kind_phys), pointer :: dudt_obl(:,:)  => null()  !<
    real (kind=kind_phys), pointer :: dvdt_obl(:,:)  => null()  !<
    real (kind=kind_phys), pointer :: dudt_oss(:,:)  => null()  !<
    real (kind=kind_phys), pointer :: dvdt_oss(:,:)  => null()  !<
    real (kind=kind_phys), pointer :: dudt_ofd(:,:)  => null()  !<
    real (kind=kind_phys), pointer :: dvdt_ofd(:,:)  => null()  !<

    real (kind=kind_phys), pointer :: du_ogwcol(:)   => null()  !< instantaneous sfc u-momentum flux from OGW
    real (kind=kind_phys), pointer :: dv_ogwcol(:)   => null()  !< instantaneous sfc v-momentum flux from OGW
    real (kind=kind_phys), pointer :: du_oblcol(:)   => null()  !< instantaneous sfc u-momentum flux from blocking
    real (kind=kind_phys), pointer :: dv_oblcol(:)   => null()  !< instantaneous sfc v-momentum flux from blocking
    real (kind=kind_phys), pointer :: du_osscol(:)   => null()  !< instantaneous sfc u-momentum flux from SSGWD
    real (kind=kind_phys), pointer :: dv_osscol(:)   => null()  !< instantaneous sfc v-momentum flux from SSGWD
    real (kind=kind_phys), pointer :: du_ofdcol(:)   => null()  !< instantaneous sfc u-momentum flux from TOFD
    real (kind=kind_phys), pointer :: dv_ofdcol(:)   => null()  !< instantaneous sfc v-momentum flux from TOFD
    real (kind=kind_phys), pointer :: du3_ogwcol(:)  => null()  !< time-averaged sfc u-momentum flux from OGW
    real (kind=kind_phys), pointer :: dv3_ogwcol(:)  => null()  !< time-averaged sfc v-momentum flux from OGW
    real (kind=kind_phys), pointer :: du3_oblcol(:)  => null()  !< time-averaged sfc u-momentum flux from blocking
    real (kind=kind_phys), pointer :: dv3_oblcol(:)  => null()  !< time-averaged sfc v-momentum flux from blocking
    real (kind=kind_phys), pointer :: du3_osscol(:)  => null()  !< time-averaged sfc u-momentum flux from SSGWD
    real (kind=kind_phys), pointer :: dv3_osscol(:)  => null()  !< time-averaged sfc v-momentum flux from SSGWD
    real (kind=kind_phys), pointer :: du3_ofdcol(:)  => null()  !< time-averaged sfc u-momentum flux from TOFD
    real (kind=kind_phys), pointer :: dv3_ofdcol(:)  => null()  !< time-averaged sfc v-momentum flux from TOFD
!
!---vay-2018 UGWP-diagnostics daily mean
!
    real (kind=kind_phys), pointer :: du3dt_ogw(:,:) => null()  !< daily aver GFS_phys tend for WE-U OGW
!
    real (kind=kind_phys), pointer :: ldu3dt_ogw(:,:) => null()  !< time aver GFS_phys tend for WE-U OGW
    real (kind=kind_phys), pointer :: ldu3dt_obl(:,:) => null()  !< time aver GFS_phys tend for WE-U OBL
    real (kind=kind_phys), pointer :: ldu3dt_oss(:,:) => null()  !< time aver GFS_phys tend for WE-U OSS
    real (kind=kind_phys), pointer :: ldu3dt_ofd(:,:) => null()  !< time aver GFS_phys tend for WE-U OFD
!
    real (kind=kind_phys), pointer :: du3dt_mtb(:,:) => null()  !< daily aver GFS_phys tend for WE-U MTB
!
    real (kind=kind_phys), pointer :: du3dt_tms(:,:) => null()  !< daily aver GFS_phys tend for WE-U TMS
!
    real (kind=kind_phys), pointer :: du3dt_ngw(:,:) => null()  !< daily aver GFS_phys tend for WE-U NGW
    real (kind=kind_phys), pointer :: dv3dt_ngw(:,:) => null()  !< daily aver GFS_phys tend for SN-V NGW
!
    real (kind=kind_phys), pointer :: dws3dt_ogw(:,:) => null()  !< time aver GFS_phys tend for windspeed OGW
    real (kind=kind_phys), pointer :: dws3dt_obl(:,:) => null()  !< time aver GFS_phys tend for windspeed OBL
    real (kind=kind_phys), pointer :: dws3dt_oss(:,:) => null()  !< time aver GFS_phys tend for windspeed OSS
    real (kind=kind_phys), pointer :: dws3dt_ofd(:,:) => null()  !< time aver GFS_phys tend for windspeed OFD
!
    real (kind=kind_phys), pointer :: ldu3dt_ngw(:,:) => null()  !< time aver GFS_phys tend for u wind NGW
    real (kind=kind_phys), pointer :: ldv3dt_ngw(:,:) => null()  !< time aver GFS_phys tend for v wind NGW
    real (kind=kind_phys), pointer :: ldt3dt_ngw(:,:) => null()  !< time aver GFS_phys tend for temperature NGW
!
!--- Instantaneous UGWP-diagnostics  16-variables
!       Diag%gwp_ax, Diag%gwp_axo, Diag%gwp_axc, Diag%gwp_axf,       &
!       Diag%gwp_ay, Diag%gwp_ayo, Diag%gwp_ayc, Diag%gwp_ayf,       &
!       Diag%gwp_dtdt,   Diag%gwp_kdis, Diag%gwp_okw, Diag%gwp_fgf,  &
!       Diag%gwp_dcheat, Diag%gwp_precip, Diag%gwp_klevs,            &
!       Diag%gwp_scheat

    real (kind=kind_phys), pointer :: gwp_scheat(:,:) => null()  ! instant shal-conv heat tendency
    real (kind=kind_phys), pointer :: gwp_dcheat(:,:) => null()  ! instant deep-conv heat tendency
    real (kind=kind_phys), pointer :: gwp_precip(:)   => null()  ! total precip rates
    integer , pointer              :: gwp_klevs(:,:)  => null()  ! instant levels for GW-launches
    real (kind=kind_phys), pointer :: gwp_fgf(:)      => null()  ! fgf triggers
    real (kind=kind_phys), pointer :: gwp_okw(:)      => null()  ! okw triggers

    real (kind=kind_phys), pointer :: gwp_ax(:,:)    => null()   ! instant total UGWP tend m/s/s EW
    real (kind=kind_phys), pointer :: gwp_ay(:,:)    => null()   ! instant total UGWP tend m/s/s NS
    real (kind=kind_phys), pointer :: gwp_dtdt(:,:)  => null()   ! instant total heat tend   K/s
    real (kind=kind_phys), pointer :: gwp_kdis(:,:)  => null()   ! instant total eddy mixing m2/s
    real (kind=kind_phys), pointer :: gwp_axc(:,:)   => null()   ! instant con-UGWP tend m/s/s EW
    real (kind=kind_phys), pointer :: gwp_ayc(:,:)   => null()   ! instant con-UGWP tend m/s/s NS
    real (kind=kind_phys), pointer :: gwp_axo(:,:)   => null()   ! instant oro-UGWP tend m/s/s EW
    real (kind=kind_phys), pointer :: gwp_ayo(:,:)   => null()   ! instant oro-UGWP tend m/s/s NS
    real (kind=kind_phys), pointer :: gwp_axf(:,:)   => null()   ! instant jet-UGWP tend m/s/s EW
    real (kind=kind_phys), pointer :: gwp_ayf(:,:)   => null()   ! instant jet-UGWP tend m/s/s NS

!--- COODRE ORO diagnostics
    real (kind=kind_phys), pointer :: zmtb(:)        => null()   !
    real (kind=kind_phys), pointer :: zogw(:)        => null()   !
    real (kind=kind_phys), pointer :: zlwb(:)        => null()   !
    real (kind=kind_phys), pointer :: tau_ogw(:)     => null()   !
    real (kind=kind_phys), pointer :: tau_ngw(:)     => null()   !
    real (kind=kind_phys), pointer :: tau_mtb(:)     => null()   !
    real (kind=kind_phys), pointer :: tau_tofd(:)    => null()   !
!---vay-2018 UGWP-diagnostics

    ! Diagnostic arrays for per-timestep diagnostics
    real (kind=kind_phys), pointer :: old_pgr(:) => null()     !< pgr at last timestep

    ! Extended output diagnostics for Thompson MP
    real (kind=kind_phys), pointer :: thompson_ext_diag3d (:,:,:) => null() ! extended diagnostic 3d output arrays from Thompson MP

    ! Diagnostics for coupled air quality model
    real (kind=kind_phys), pointer :: aod   (:)   => null()    !< instantaneous aerosol optical depth ( n/a )

    ! Auxiliary output arrays for debugging
    real (kind=kind_phys), pointer :: aux2d(:,:)  => null()    !< auxiliary 2d arrays in output (for debugging)
    real (kind=kind_phys), pointer :: aux3d(:,:,:)=> null()    !< auxiliary 2d arrays in output (for debugging)

    !--- Lightning threat indices
    real (kind=kind_phys), pointer :: ltg1_max(:)        => null()  !
    real (kind=kind_phys), pointer :: ltg2_max(:)        => null()  !
    real (kind=kind_phys), pointer :: ltg3_max(:)        => null()  !

    !--- NRL Ozone physics diagnostics
    real (kind=kind_phys), pointer :: do3_dt_prd(:,:)  => null()
    real (kind=kind_phys), pointer :: do3_dt_ozmx(:,:) => null()
    real (kind=kind_phys), pointer :: do3_dt_temp(:,:) => null()
    real (kind=kind_phys), pointer :: do3_dt_ohoz(:,:) => null()

    contains
      procedure :: create    => diag_create
      procedure :: rad_zero  => diag_rad_zero
      procedure :: phys_zero => diag_phys_zero
  end type GFS_diag_type

!----------------
! PUBLIC ENTITIES
!----------------
  public GFS_init_type
  public GFS_statein_type,  GFS_stateout_type, GFS_sfcprop_type, &
         GFS_coupling_type
  public GFS_control_type,  GFS_grid_type,     GFS_tbd_type, &
         GFS_cldprop_type,  GFS_radtend_type,  GFS_diag_type

!*******************************************************************************************
  CONTAINS

!------------------------
! GFS_statein_type%create
!------------------------
  subroutine statein_create (Statein, Model)
    implicit none

    class(GFS_statein_type)             :: Statein
    type(GFS_control_type),  intent(in) :: Model
    integer :: IM

    IM = Model%ncols

    !--- level geopotential and pressures
    allocate (Statein%phii  (IM,Model%levs+1))
    allocate (Statein%prsi  (IM,Model%levs+1))
    allocate (Statein%prsik (IM,Model%levs+1))

    Statein%phii  = clear_val
    Statein%prsi  = clear_val
    Statein%prsik = clear_val

    !--- layer geopotential and pressures
    allocate (Statein%phil  (IM,Model%levs))
    allocate (Statein%prsl  (IM,Model%levs))
    allocate (Statein%prslk (IM,Model%levs))

    Statein%phil  = clear_val
    Statein%prsl  = clear_val
    Statein%prslk = clear_val

    !--- shared radiation and physics variables
    allocate (Statein%vvl  (IM,Model%levs))
    allocate (Statein%tgrs (IM,Model%levs))

    Statein%vvl  = clear_val
    Statein%tgrs = clear_val
! stochastic physics SKEB variable
    allocate (Statein%diss_est(IM,Model%levs))
    Statein%diss_est= clear_val
    !--- physics only variables
    allocate (Statein%pgr    (IM))
    allocate (Statein%ugrs   (IM,Model%levs))
    allocate (Statein%vgrs   (IM,Model%levs))
    if(Model%lightning_threat) then
      allocate (Statein%wgrs   (IM,Model%levs))
    endif
    allocate (Statein%qgrs   (IM,Model%levs,Model%ntrac))

    Statein%qgrs   = clear_val
    Statein%pgr    = clear_val
    Statein%ugrs   = clear_val
    Statein%vgrs   = clear_val

    if(Model%lightning_threat) then
      Statein%wgrs = clear_val
    endif

    !--- soil state variables - for soil SPPT - sfc-perts, mgehne
    allocate (Statein%smc  (IM,Model%lsoil))
    allocate (Statein%stc  (IM,Model%lsoil))
    allocate (Statein%slc  (IM,Model%lsoil))

    Statein%smc   = clear_val
    Statein%stc   = clear_val
    Statein%slc   = clear_val

  end subroutine statein_create


!-------------------------
! GFS_stateout_type%create
!-------------------------
  subroutine stateout_create (Stateout, Model)

    implicit none

    class(GFS_stateout_type)           :: Stateout
    type(GFS_control_type), intent(in) :: Model
    integer :: IM

    IM = Model%ncols

    allocate (Stateout%gu0 (IM,Model%levs))
    allocate (Stateout%gv0 (IM,Model%levs))
    allocate (Stateout%gt0 (IM,Model%levs))
    allocate (Stateout%gq0 (IM,Model%levs,Model%ntrac))

    Stateout%gu0 = clear_val
    Stateout%gv0 = clear_val
    Stateout%gt0 = clear_val
    Stateout%gq0 = clear_val

 end subroutine stateout_create


!------------------------
! GFS_sfcprop_type%create
!------------------------
  subroutine sfcprop_create (Sfcprop, Model)

    implicit none

    class(GFS_sfcprop_type)            :: Sfcprop
    type(GFS_control_type), intent(in) :: Model
    integer :: IM

    IM = Model%ncols

    !--- physics and radiation
    allocate (Sfcprop%slmsk    (IM))
    allocate (Sfcprop%oceanfrac(IM))
    allocate (Sfcprop%landfrac (IM))
    allocate (Sfcprop%vegtype_frac (IM,Model%nvegcat))
    allocate (Sfcprop%soiltype_frac(IM,Model%nsoilcat))
    allocate (Sfcprop%lakefrac (IM))
    allocate (Sfcprop%lakedepth(IM))

    allocate (Sfcprop%use_lake_model(IM))

    if(Model%lkm > 0) then
      if(Model%iopt_lake==Model%iopt_lake_clm) then
        allocate (Sfcprop%clm_lakedepth(IM))
      else if(Model%iopt_lake==Model%iopt_lake_flake) then
        allocate (Sfcprop%h_ML     (IM))
        allocate (Sfcprop%t_ML     (IM))
        allocate (Sfcprop%t_mnw    (IM))
        allocate (Sfcprop%h_talb   (IM))
        allocate (Sfcprop%t_talb   (IM))
        allocate (Sfcprop%t_bot1   (IM))
        allocate (Sfcprop%t_bot2   (IM))
        allocate (Sfcprop%c_t      (IM))
      endif
      allocate (Sfcprop%T_snow   (IM))
      allocate (Sfcprop%T_ice    (IM))
    endif

    allocate (Sfcprop%tsfc     (IM))
    allocate (Sfcprop%usfco    (IM))
    allocate (Sfcprop%vsfco    (IM))
    allocate (Sfcprop%tsfco    (IM))
    allocate (Sfcprop%tsfcl    (IM))
    allocate (Sfcprop%tisfc    (IM))
    allocate (Sfcprop%tiice    (IM,Model%kice))
    allocate (Sfcprop%snowd    (IM))
    allocate (Sfcprop%zorl     (IM))
    allocate (Sfcprop%zorlw    (IM))
    allocate (Sfcprop%zorll    (IM))
    allocate (Sfcprop%zorli    (IM))
    allocate (Sfcprop%zorlwav  (IM))
    allocate (Sfcprop%fice     (IM))
    allocate (Sfcprop%snodl    (IM))
    allocate (Sfcprop%weasdl   (IM))
    allocate (Sfcprop%snodi    (IM))
    allocate (Sfcprop%weasdi   (IM))
    allocate (Sfcprop%hprime   (IM,Model%nmtvr))
    allocate (Sfcprop%dust12m_in  (IM,12,5))
    allocate (Sfcprop%smoke_RRFS(IM,24,2))
    allocate (Sfcprop%smoke2d_RRFS(IM,4))
    allocate (Sfcprop%emi_in   (IM,1))
    allocate(Sfcprop%albdirvis_lnd (IM))
    allocate(Sfcprop%albdirnir_lnd (IM))
    allocate(Sfcprop%albdifvis_lnd (IM))
    allocate(Sfcprop%albdifnir_lnd (IM))
    allocate (Sfcprop%emis_lnd (IM))
    allocate (Sfcprop%emis_ice (IM))
    allocate (Sfcprop%emis_wat (IM))

    Sfcprop%slmsk     = clear_val
    Sfcprop%oceanfrac = clear_val
    Sfcprop%landfrac  = clear_val
    Sfcprop%vegtype_frac  = clear_val
    Sfcprop%soiltype_frac = clear_val
    Sfcprop%lakefrac  = clear_val
    Sfcprop%lakedepth = clear_val

    Sfcprop%use_lake_model = zero
    if(Model%lkm > 0) then
      if(Model%iopt_lake==Model%iopt_lake_clm) then
        Sfcprop%clm_lakedepth = clear_val
      else if(Model%iopt_lake==Model%iopt_lake_flake) then
        Sfcprop%h_ML      = clear_val
        Sfcprop%t_ML      = clear_val
        Sfcprop%t_mnw     = clear_val
        Sfcprop%h_talb    = clear_val
        Sfcprop%t_talb    = clear_val
        Sfcprop%t_bot1    = clear_val
        Sfcprop%t_bot2    = clear_val
        Sfcprop%c_t       = clear_val
      endif
      Sfcprop%T_snow    = clear_val
      Sfcprop%T_ice     = clear_val
    endif

    Sfcprop%tsfc      = clear_val
    Sfcprop%usfco     = clear_val
    Sfcprop%vsfco     = clear_val
    Sfcprop%tsfco     = clear_val
    Sfcprop%tsfcl     = clear_val
    Sfcprop%tisfc     = clear_val
    Sfcprop%tiice     = clear_val
    Sfcprop%snowd     = clear_val
    Sfcprop%zorl      = clear_val
    Sfcprop%zorlw     = clear_val
    Sfcprop%zorll     = clear_val
    Sfcprop%zorli     = clear_val
    Sfcprop%zorlwav   = clear_val
    Sfcprop%fice      = clear_val
    Sfcprop%snodl     = clear_val
    Sfcprop%weasdl    = clear_val
    Sfcprop%snodi     = clear_val
    Sfcprop%weasdi    = clear_val
    Sfcprop%hprime    = clear_val
    Sfcprop%dust12m_in= clear_val
    Sfcprop%emi_in    = clear_val
    Sfcprop%smoke_RRFS= clear_val
    Sfcprop%smoke2d_RRFS= clear_val
    Sfcprop%albdirvis_lnd = clear_val
    Sfcprop%albdirnir_lnd = clear_val
    Sfcprop%albdifvis_lnd = clear_val
    Sfcprop%albdifnir_lnd = clear_val
    Sfcprop%emis_lnd  = clear_val
    Sfcprop%emis_ice  = clear_val
    Sfcprop%emis_wat  = clear_val

!--- In (radiation only)
    allocate (Sfcprop%snoalb (IM))
    allocate (Sfcprop%alvsf  (IM))
    allocate (Sfcprop%alnsf  (IM))
    allocate (Sfcprop%alvwf  (IM))
    allocate (Sfcprop%alnwf  (IM))
    allocate (Sfcprop%facsf  (IM))
    allocate (Sfcprop%facwf  (IM))

    Sfcprop%snoalb = clear_val
    Sfcprop%alvsf  = clear_val
    Sfcprop%alnsf  = clear_val
    Sfcprop%alvwf  = clear_val
    Sfcprop%alnwf  = clear_val
    Sfcprop%facsf  = clear_val
    Sfcprop%facwf  = clear_val

!--- physics surface props
!--- In
    allocate (Sfcprop%slope      (IM))
    allocate (Sfcprop%slope_save (IM))
    allocate (Sfcprop%shdmin     (IM))
    allocate (Sfcprop%shdmax     (IM))
    allocate (Sfcprop%snoalb     (IM))
    allocate (Sfcprop%tg3        (IM))
    allocate (Sfcprop%vfrac      (IM))
    allocate (Sfcprop%vtype      (IM))
    allocate (Sfcprop%vtype_save (IM))
    allocate (Sfcprop%stype      (IM))
    allocate (Sfcprop%stype_save (IM))
    allocate (Sfcprop%scolor     (IM))
    allocate (Sfcprop%scolor_save(IM))
    allocate (Sfcprop%uustar     (IM))
    allocate (Sfcprop%oro        (IM))
    allocate (Sfcprop%oro_uf     (IM))
    allocate (Sfcprop%evap       (IM))
    allocate (Sfcprop%hflx       (IM))
    allocate (Sfcprop%qss        (IM))

    Sfcprop%slope      = zero
    Sfcprop%slope_save = zero
    Sfcprop%shdmin     = clear_val
    Sfcprop%shdmax     = clear_val
    Sfcprop%snoalb     = clear_val
    Sfcprop%tg3        = clear_val
    Sfcprop%vfrac      = clear_val
    Sfcprop%vtype      = zero
    Sfcprop%vtype_save = zero
    Sfcprop%stype      = zero
    Sfcprop%stype_save = zero
    Sfcprop%scolor      = zero
    Sfcprop%scolor_save = zero
    Sfcprop%uustar     = clear_val
    Sfcprop%oro        = clear_val
    Sfcprop%oro_uf     = clear_val
    Sfcprop%evap       = clear_val
    Sfcprop%hflx       = clear_val
    Sfcprop%qss        = clear_val

!--- In/Out
    allocate (Sfcprop%hice   (IM))
    allocate (Sfcprop%weasd  (IM))
    allocate (Sfcprop%sncovr (IM))
    allocate (Sfcprop%sncovr_ice (IM))
    if (Model%use_cice_alb .or. Model%lsm == Model%lsm_ruc) then
      allocate (Sfcprop%albdirvis_ice (IM))
      allocate (Sfcprop%albdifvis_ice (IM))
      allocate (Sfcprop%albdirnir_ice (IM))
      allocate (Sfcprop%albdifnir_ice (IM))
    endif
    if (Model%lsm == Model%lsm_ruc) then
      allocate (Sfcprop%sfalb_lnd (IM))
      allocate (Sfcprop%sfalb_ice (IM))
      allocate (Sfcprop%sfalb_lnd_bck (IM))
    endif
    allocate (Sfcprop%canopy (IM))
    allocate (Sfcprop%ffmm   (IM))
    allocate (Sfcprop%ffhh   (IM))
    allocate (Sfcprop%f10m   (IM))
    allocate (Sfcprop%tprcp  (IM))
    allocate (Sfcprop%srflag (IM))
    allocate (Sfcprop%slc    (IM,Model%lsoil))
    allocate (Sfcprop%smc    (IM,Model%lsoil))
    allocate (Sfcprop%stc    (IM,Model%lsoil))

    Sfcprop%hice   = clear_val
    Sfcprop%weasd  = clear_val
    Sfcprop%sncovr = clear_val
    Sfcprop%sncovr_ice = clear_val
    if (Model%use_cice_alb .or. Model%lsm == Model%lsm_ruc) then
      Sfcprop%albdirvis_ice = clear_val
      Sfcprop%albdifvis_ice = clear_val
      Sfcprop%albdirnir_ice = clear_val
      Sfcprop%albdifnir_ice = clear_val
    endif
    if (Model%lsm == Model%lsm_ruc) then
      Sfcprop%sfalb_lnd     = clear_val
      Sfcprop%sfalb_ice     = clear_val
      Sfcprop%sfalb_lnd_bck = clear_val
    endif
    Sfcprop%canopy = clear_val
    Sfcprop%ffmm   = clear_val
    Sfcprop%ffhh   = clear_val
    Sfcprop%f10m   = clear_val
    Sfcprop%tprcp  = clear_val
    Sfcprop%srflag = clear_val
    Sfcprop%slc    = clear_val
    Sfcprop%smc    = clear_val
    Sfcprop%stc    = clear_val

!--- Out
    allocate (Sfcprop%t2m (IM))
    allocate (Sfcprop%th2m(IM))
    allocate (Sfcprop%q2m (IM))

    Sfcprop%t2m  = clear_val
    Sfcprop%th2m = clear_val
    Sfcprop%q2m  = clear_val

    if (Model%nstf_name(1) > 0) then
      allocate (Sfcprop%tref   (IM))
      allocate (Sfcprop%z_c    (IM))
      allocate (Sfcprop%c_0    (IM))
      allocate (Sfcprop%c_d    (IM))
      allocate (Sfcprop%w_0    (IM))
      allocate (Sfcprop%w_d    (IM))
      allocate (Sfcprop%xt     (IM))
      allocate (Sfcprop%xs     (IM))
      allocate (Sfcprop%xu     (IM))
      allocate (Sfcprop%xv     (IM))
      allocate (Sfcprop%xz     (IM))
      allocate (Sfcprop%zm     (IM))
      allocate (Sfcprop%xtts   (IM))
      allocate (Sfcprop%xzts   (IM))
      allocate (Sfcprop%d_conv (IM))
      allocate (Sfcprop%ifd    (IM))
      allocate (Sfcprop%dt_cool(IM))
      allocate (Sfcprop%qrain  (IM))

      Sfcprop%tref    = zero
      Sfcprop%z_c     = zero
      Sfcprop%c_0     = zero
      Sfcprop%c_d     = zero
      Sfcprop%w_0     = zero
      Sfcprop%w_d     = zero
      Sfcprop%xt      = zero
      Sfcprop%xs      = zero
      Sfcprop%xu      = zero
      Sfcprop%xv      = zero
      Sfcprop%xz      = zero
      Sfcprop%zm      = zero
      Sfcprop%xtts    = zero
      Sfcprop%xzts    = zero
      Sfcprop%d_conv  = zero
      Sfcprop%ifd     = zero
      Sfcprop%dt_cool = zero
      Sfcprop%qrain   = zero
    endif
    if (Model%lsm == Model%lsm_noah .or. Model%lsm == Model%lsm_noahmp .or. Model%lsm == Model%lsm_ruc) then
      allocate (Sfcprop%xlaixy   (IM))
      Sfcprop%xlaixy     = clear_val
    else
      allocate (Sfcprop%xlaixy   (0))
    end if
    if (Model%lsm == Model%lsm_noah .or. Model%lsm == Model%lsm_noahmp) then
      allocate (Sfcprop%rca      (IM))
      Sfcprop%rca        = clear_val
    else
      allocate (Sfcprop%rca      (0))
    end if
    if (Model%lsm == Model%lsm_ruc .or. Model%lsm == Model%lsm_noahmp .or. &
         (Model%lkm>0 .and. Model%iopt_lake==Model%iopt_lake_clm)) then
     allocate(Sfcprop%raincprv  (IM))
     allocate(Sfcprop%rainncprv (IM))
     Sfcprop%raincprv   = clear_val
     Sfcprop%rainncprv  = clear_val
     if (Model%lsm == Model%lsm_ruc .or. Model%lsm == Model%lsm_noahmp) then
      allocate(Sfcprop%iceprv    (IM))
      allocate(Sfcprop%snowprv   (IM))
      allocate(Sfcprop%graupelprv(IM))
      Sfcprop%iceprv     = clear_val
      Sfcprop%snowprv    = clear_val
      Sfcprop%graupelprv = clear_val
     end if
    end if
! Noah MP allocate and init when used
!
    if (Model%lsm == Model%lsm_noahmp ) then

      allocate (Sfcprop%snowxy   (IM))
      allocate (Sfcprop%tvxy     (IM))
      allocate (Sfcprop%tgxy     (IM))
      allocate (Sfcprop%canicexy (IM))
      allocate (Sfcprop%canliqxy (IM))
      allocate (Sfcprop%eahxy    (IM))
      allocate (Sfcprop%tahxy    (IM))
      allocate (Sfcprop%cmxy     (IM))
      allocate (Sfcprop%chxy     (IM))
      allocate (Sfcprop%fwetxy   (IM))
      allocate (Sfcprop%sneqvoxy (IM))
      allocate (Sfcprop%alboldxy (IM))
      allocate (Sfcprop%qsnowxy  (IM))
      allocate (Sfcprop%wslakexy (IM))
      allocate (Sfcprop%zwtxy    (IM))
      allocate (Sfcprop%waxy     (IM))
      allocate (Sfcprop%wtxy     (IM))
      allocate (Sfcprop%lfmassxy (IM))
      allocate (Sfcprop%rtmassxy (IM))
      allocate (Sfcprop%stmassxy (IM))
      allocate (Sfcprop%woodxy   (IM))
      allocate (Sfcprop%stblcpxy (IM))
      allocate (Sfcprop%fastcpxy (IM))
      allocate (Sfcprop%xsaixy   (IM))
      allocate (Sfcprop%taussxy  (IM))
      allocate (Sfcprop%smcwtdxy (IM))
      allocate (Sfcprop%deeprechxy (IM))
      allocate (Sfcprop%rechxy     (IM))
      allocate (Sfcprop%snicexy    (IM, Model%lsnow_lsm_lbound:Model%lsnow_lsm_ubound))
      allocate (Sfcprop%snliqxy    (IM, Model%lsnow_lsm_lbound:Model%lsnow_lsm_ubound))
      allocate (Sfcprop%tsnoxy     (IM, Model%lsnow_lsm_lbound:Model%lsnow_lsm_ubound))
      allocate (Sfcprop%smoiseq    (IM, Model%lsoil_lsm))
      allocate (Sfcprop%zsnsoxy    (IM, Model%lsnow_lsm_lbound:Model%lsoil_lsm))

      Sfcprop%snowxy     = clear_val
      Sfcprop%tvxy       = clear_val
      Sfcprop%tgxy       = clear_val
      Sfcprop%canicexy   = clear_val
      Sfcprop%canliqxy   = clear_val
      Sfcprop%eahxy      = clear_val
      Sfcprop%tahxy      = clear_val
      Sfcprop%cmxy       = clear_val
      Sfcprop%chxy       = clear_val
      Sfcprop%fwetxy     = clear_val
      Sfcprop%sneqvoxy   = clear_val
      Sfcprop%alboldxy   = clear_val
      Sfcprop%qsnowxy    = clear_val
      Sfcprop%wslakexy   = clear_val
      Sfcprop%zwtxy      = clear_val
      Sfcprop%waxy       = clear_val
      Sfcprop%wtxy       = clear_val
      Sfcprop%lfmassxy   = clear_val
      Sfcprop%rtmassxy   = clear_val
      Sfcprop%stmassxy   = clear_val
      Sfcprop%woodxy     = clear_val
      Sfcprop%stblcpxy   = clear_val
      Sfcprop%fastcpxy   = clear_val
      Sfcprop%xsaixy     = clear_val
      Sfcprop%taussxy    = clear_val
      Sfcprop%smcwtdxy   = clear_val
      Sfcprop%deeprechxy = clear_val
      Sfcprop%rechxy     = clear_val

      Sfcprop%snicexy    = clear_val
      Sfcprop%snliqxy    = clear_val
      Sfcprop%tsnoxy     = clear_val
      Sfcprop%smoiseq    = clear_val
      Sfcprop%zsnsoxy    = clear_val

      allocate(Sfcprop%draincprv  (IM))
      allocate(Sfcprop%drainncprv (IM))
      allocate(Sfcprop%diceprv    (IM))
      allocate(Sfcprop%dsnowprv   (IM))
      allocate(Sfcprop%dgraupelprv(IM))

      Sfcprop%draincprv   = clear_val
      Sfcprop%drainncprv  = clear_val
      Sfcprop%diceprv     = clear_val
      Sfcprop%dsnowprv    = clear_val
      Sfcprop%dgraupelprv = clear_val

    endif

    if (Model%do_myjsfc .or. Model%do_myjpbl) then
      allocate(Sfcprop%z0base(IM))
      Sfcprop%z0base = clear_val
    end if

    allocate(Sfcprop%semisbase(IM))
    Sfcprop%semisbase = clear_val

    if (Model%lsm == Model%lsm_ruc) then
       ! For land surface models with different numbers of levels than the four NOAH levels
       allocate (Sfcprop%wetness         (IM))
       allocate (Sfcprop%sh2o            (IM,Model%lsoil_lsm))
       allocate (Sfcprop%keepsmfr        (IM,Model%lsoil_lsm))
       allocate (Sfcprop%smois           (IM,Model%lsoil_lsm))
       allocate (Sfcprop%tslb            (IM,Model%lsoil_lsm))
       allocate (Sfcprop%flag_frsoil     (IM,Model%lsoil_lsm))
       allocate (Sfcprop%clw_surf_land   (IM))
       allocate (Sfcprop%clw_surf_ice    (IM))
       allocate (Sfcprop%qwv_surf_land   (IM))
       allocate (Sfcprop%qwv_surf_ice    (IM))
       allocate (Sfcprop%rhofr           (IM))
       allocate (Sfcprop%tsnow_land      (IM))
       allocate (Sfcprop%tsnow_ice       (IM))
       allocate (Sfcprop%snowfallac_land (IM))
       allocate (Sfcprop%snowfallac_ice  (IM))
       allocate (Sfcprop%acsnow_land     (IM))
       allocate (Sfcprop%acsnow_ice      (IM))
       allocate (Sfcprop%fire_heat_flux  (IM))
       allocate (Sfcprop%frac_grid_burned(IM))

       !
       Sfcprop%wetness         = clear_val
       Sfcprop%sh2o            = clear_val
       Sfcprop%keepsmfr        = clear_val
       Sfcprop%smois           = clear_val
       Sfcprop%tslb            = clear_val
       Sfcprop%clw_surf_land   = clear_val
       Sfcprop%clw_surf_ice    = clear_val
       Sfcprop%qwv_surf_land   = clear_val
       Sfcprop%qwv_surf_ice    = clear_val
       Sfcprop%flag_frsoil     = clear_val
       Sfcprop%rhofr           = -1.e3
       Sfcprop%tsnow_land      = clear_val
       Sfcprop%tsnow_ice       = clear_val
       Sfcprop%snowfallac_land = clear_val
       Sfcprop%snowfallac_ice  = clear_val
       Sfcprop%acsnow_land     = clear_val
       Sfcprop%acsnow_ice      = clear_val
       Sfcprop%fire_heat_flux  = clear_val
       Sfcprop%frac_grid_burned= clear_val
       !
    end if

    allocate (Sfcprop%rmol   (IM ))
    allocate (Sfcprop%flhc   (IM ))
    allocate (Sfcprop%flqc   (IM ))
    Sfcprop%rmol        = clear_val
    Sfcprop%flhc        = clear_val
    Sfcprop%flqc        = clear_val

    if (Model%do_mynnsfclay) then
    ! For MYNN surface layer scheme
       !print*,"Allocating all MYNN-sfclay variables"
       allocate (Sfcprop%ustm   (IM ))
       allocate (Sfcprop%zol    (IM ))
       allocate (Sfcprop%mol    (IM ))
       allocate (Sfcprop%chs2   (IM ))
       allocate (Sfcprop%cqs2   (IM ))
       allocate (Sfcprop%lh     (IM ))
       !
       !print*,"Initializing all MYNN-SfcLay variables with ",clear_val
       Sfcprop%ustm        = clear_val
       Sfcprop%zol         = clear_val
       Sfcprop%mol         = clear_val
       Sfcprop%chs2        = clear_val
       Sfcprop%cqs2        = clear_val
       Sfcprop%lh          = clear_val
    end if
    if (Model%imfdeepcnv == Model%imfdeepcnv_gf .or. Model%imfdeepcnv == Model%imfdeepcnv_c3) then
        allocate (Sfcprop%maxupmf(IM))
        allocate (Sfcprop%conv_act(IM))
        allocate (Sfcprop%conv_act_m(IM))
        Sfcprop%maxupmf = zero
        Sfcprop%conv_act = zero
        Sfcprop%conv_act_m = zero
    end if

    ! CLM Lake Model variables
    if (Model%lkm/=0 .and. Model%iopt_lake==Model%iopt_lake_clm) then
       allocate(Sfcprop%lake_t2m(IM))
       allocate(Sfcprop%lake_q2m(IM))
       allocate(Sfcprop%lake_albedo(IM))
       allocate(Sfcprop%input_lakedepth(IM))
       allocate(Sfcprop%lake_h2osno2d(IM))
       allocate(Sfcprop%lake_sndpth2d(IM))
       allocate(Sfcprop%lake_snl2d(IM))
       allocate(Sfcprop%lake_snow_z3d(IM,Model%nlevsnowsoil1_clm_lake))
       allocate(Sfcprop%lake_snow_dz3d(IM,Model%nlevsnowsoil1_clm_lake))
       allocate(Sfcprop%lake_snow_zi3d(IM,Model%nlevsnowsoil_clm_lake))
       allocate(Sfcprop%lake_h2osoi_vol3d(IM,Model%nlevsnowsoil1_clm_lake))
       allocate(Sfcprop%lake_h2osoi_liq3d(IM,Model%nlevsnowsoil1_clm_lake))
       allocate(Sfcprop%lake_h2osoi_ice3d(IM,Model%nlevsnowsoil1_clm_lake))
       allocate(Sfcprop%lake_tsfc(IM))
       allocate(Sfcprop%lake_t_soisno3d(IM,Model%nlevsnowsoil1_clm_lake))
       allocate(Sfcprop%lake_t_lake3d(IM,Model%nlevlake_clm_lake))
       allocate(Sfcprop%lake_savedtke12d(IM))
       allocate(Sfcprop%lake_icefrac3d(IM,Model%nlevlake_clm_lake))
       allocate(Sfcprop%lake_rho0(IM))
       allocate(Sfcprop%lake_ht(IM))
       allocate(Sfcprop%lake_is_salty(IM))
       allocate(Sfcprop%lake_cannot_freeze(IM))
       allocate(Sfcprop%clm_lake_initialized(IM))

       Sfcprop%lake_t2m = clear_val
       Sfcprop%lake_q2m = clear_val
       Sfcprop%lake_albedo = clear_val
       Sfcprop%input_lakedepth = clear_val
       Sfcprop%lake_h2osno2d = clear_val
       Sfcprop%lake_sndpth2d = clear_val
       Sfcprop%lake_snl2d = clear_val
       Sfcprop%lake_snow_z3d = clear_val
       Sfcprop%lake_snow_dz3d = clear_val
       Sfcprop%lake_snow_zi3d = clear_val
       Sfcprop%lake_h2osoi_vol3d = clear_val
       Sfcprop%lake_h2osoi_liq3d = clear_val
       Sfcprop%lake_h2osoi_ice3d = clear_val
       Sfcprop%lake_tsfc = clear_val
       Sfcprop%lake_t_soisno3d = clear_val
       Sfcprop%lake_t_lake3d = clear_val
       Sfcprop%lake_savedtke12d = clear_val
       Sfcprop%lake_icefrac3d = clear_val
       Sfcprop%lake_rho0 = -111
       Sfcprop%lake_ht = -111
       Sfcprop%lake_is_salty = zero
       Sfcprop%lake_cannot_freeze = zero
       Sfcprop%clm_lake_initialized = zero
    endif

    if(Model%rrfs_sd) then
    !--- needed for smoke aerosol option
      allocate (Sfcprop%emdust    (IM))
      allocate (Sfcprop%emseas    (IM))
      allocate (Sfcprop%emanoc    (IM))
      allocate (Sfcprop%ebb_smoke_in (IM))
      allocate (Sfcprop%frp_output (IM))
      allocate (Sfcprop%fhist     (IM))
      allocate (Sfcprop%coef_bb_dc(IM))
      allocate (Sfcprop%fire_type (IM))
      allocate (Sfcprop%peak_hr   (IM))
      allocate (Sfcprop%lu_nofire (IM))
      allocate (Sfcprop%lu_qfire  (IM))
      allocate (Sfcprop%fire_in   (IM,Model%fire_aux_data_levels))

      ! IMPORTANT: This initialization must match rrfs_sd_fill_data
      Sfcprop%emdust     = clear_val
      Sfcprop%emseas     = clear_val
      Sfcprop%emanoc     = clear_val
      Sfcprop%ebb_smoke_in = clear_val
      Sfcprop%frp_output  = clear_val
      Sfcprop%fhist      = 1.
      Sfcprop%coef_bb_dc = clear_val
      Sfcprop%fire_type  = 0
      Sfcprop%fire_in    = clear_val
      Sfcprop%peak_hr    = clear_val
      Sfcprop%lu_nofire  = clear_val
      Sfcprop%lu_qfire   = clear_val
    endif

    !--- if fire_behavior is on
    if(Model%cpl_fire) then
      allocate (Sfcprop%hflx_fire   (IM))
      allocate (Sfcprop%evap_fire   (IM))
      allocate (Sfcprop%smoke_fire  (IM))
      Sfcprop%hflx_fire = zero
      Sfcprop%evap_fire = zero
      Sfcprop%smoke_fire = zero
    endif

  end subroutine sfcprop_create


!-------------------------
! GFS_coupling_type%create
!-------------------------
  subroutine coupling_create (Coupling, Model)

    implicit none

    class(GFS_coupling_type)           :: Coupling
    type(GFS_control_type), intent(in) :: Model
    integer :: IM

    IM = Model%ncols

    !--- radiation out
    !--- physics in
    allocate (Coupling%nirbmdi (IM))
    allocate (Coupling%nirdfdi (IM))
    allocate (Coupling%visbmdi (IM))
    allocate (Coupling%visdfdi (IM))
    allocate (Coupling%nirbmui (IM))
    allocate (Coupling%nirdfui (IM))
    allocate (Coupling%visbmui (IM))
    allocate (Coupling%visdfui (IM))

    Coupling%nirbmdi = clear_val
    Coupling%nirdfdi = clear_val
    Coupling%visbmdi = clear_val
    Coupling%visdfdi = clear_val
    Coupling%nirbmui = clear_val
    Coupling%nirdfui = clear_val
    Coupling%visbmui = clear_val
    Coupling%visdfui = clear_val

    allocate (Coupling%sfcdsw (IM))
    allocate (Coupling%sfcnsw (IM))
    allocate (Coupling%sfcdlw (IM))
    allocate (Coupling%sfculw (IM))

    Coupling%sfcdsw = clear_val
    Coupling%sfcnsw = clear_val
    Coupling%sfcdlw = clear_val
    Coupling%sfculw = clear_val

    if (Model%do_RRTMGP) then
       allocate (Coupling%fluxlwUP_radtime   (IM, Model%levs+1))
       allocate (Coupling%fluxlwDOWN_radtime (IM, Model%levs+1))
       allocate (Coupling%fluxlwUP_jac       (IM, Model%levs+1))
       allocate (Coupling%htrlw              (IM, Model%levs))
       allocate (Coupling%tsfc_radtime       (IM))
       Coupling%fluxlwUP_radtime   = clear_val
       Coupling%fluxlwDOWN_radtime = clear_val
       Coupling%fluxlwUP_jac       = clear_val
       Coupling%htrlw              = clear_val
       Coupling%tsfc_radtime       = clear_val
    endif

    if (Model%cplflx .or. Model%do_sppt .or. Model%cplchm .or. Model%ca_global .or. Model%cpllnd .or. Model%cpl_fire) then
      allocate (Coupling%rain_cpl (IM))
      allocate (Coupling%snow_cpl (IM))
      Coupling%rain_cpl = clear_val
      Coupling%snow_cpl = clear_val
    endif

    if (Model%cplflx .or. Model%cplchm .or. Model%cplwav) then
      !--- instantaneous quantities
      allocate (Coupling%u10mi_cpl (IM))
      allocate (Coupling%v10mi_cpl (IM))

      Coupling%u10mi_cpl = clear_val
      Coupling%v10mi_cpl = clear_val
    endif

    if (Model%cplflx .or. Model%cplchm .or. Model%cpllnd) then
      !--- instantaneous quantities
      allocate (Coupling%tsfci_cpl (IM))
      Coupling%tsfci_cpl = clear_val
    endif

!   if (Model%cplwav2atm) then
      !--- incoming quantities
!     allocate (Coupling%zorlwav_cpl (IM))

!     Coupling%zorlwav_cpl  = clear_val
!   endif

    ! -- additional coupling options for air quality
    if (Model%cplflx .or. Model%cpllnd .or. Model%cpl_fire .or. (Model%cplaqm .and. .not.Model%cplflx)) then
      allocate (Coupling%psurfi_cpl  (IM))
      allocate (Coupling%nswsfci_cpl (IM))
      Coupling%psurfi_cpl  = clear_val
      Coupling%nswsfci_cpl = clear_val
    endif
     
    if (Model%cplflx .or. Model%cpl_fire .or. (Model%cplaqm .and. .not.Model%cplflx)) then
      allocate (Coupling%dtsfci_cpl  (IM))
      allocate (Coupling%dqsfci_cpl  (IM))
      allocate (Coupling%t2mi_cpl    (IM))
      allocate (Coupling%q2mi_cpl    (IM))
      Coupling%dtsfci_cpl  = clear_val
      Coupling%dqsfci_cpl  = clear_val
      Coupling%t2mi_cpl    = clear_val
      Coupling%q2mi_cpl    = clear_val
    endif
    
    if (Model%cplflx .or. Model%cpllnd) then
      allocate (Coupling%dlwsfci_cpl (IM))
      allocate (Coupling%dswsfci_cpl (IM))
      allocate (Coupling%dlwsfc_cpl  (IM))
      allocate (Coupling%dswsfc_cpl  (IM))
      allocate (Coupling%nswsfc_cpl  (IM))
      allocate (Coupling%nnirbmi_cpl (IM))
      allocate (Coupling%nnirdfi_cpl (IM))
      allocate (Coupling%nvisbmi_cpl (IM))
      allocate (Coupling%nvisdfi_cpl (IM))
      allocate (Coupling%nnirbm_cpl  (IM))
      allocate (Coupling%nnirdf_cpl  (IM))
      allocate (Coupling%nvisbm_cpl  (IM))
      allocate (Coupling%nvisdf_cpl  (IM))

      Coupling%dlwsfci_cpl  = clear_val
      Coupling%dswsfci_cpl  = clear_val
      Coupling%dlwsfc_cpl  = clear_val
      Coupling%dswsfc_cpl  = clear_val
      Coupling%nswsfc_cpl  = clear_val
      Coupling%nnirbmi_cpl = clear_val
      Coupling%nnirdfi_cpl = clear_val
      Coupling%nvisbmi_cpl = clear_val
      Coupling%nvisdfi_cpl = clear_val
      Coupling%nnirbm_cpl  = clear_val
      Coupling%nnirdf_cpl  = clear_val
      Coupling%nvisbm_cpl  = clear_val
      Coupling%nvisdf_cpl  = clear_val
    end if

    if (Model%cplflx .or. Model%cpl_fire) then
      !--- incoming quantities
      allocate (Coupling%slimskin_cpl        (IM))
      allocate (Coupling%dusfcin_cpl         (IM))
      allocate (Coupling%dvsfcin_cpl         (IM))
      allocate (Coupling%dtsfcin_cpl         (IM))
      allocate (Coupling%dqsfcin_cpl         (IM))
      allocate (Coupling%ulwsfcin_cpl        (IM))
!     allocate (Coupling%tseain_cpl          (IM))
!     allocate (Coupling%tisfcin_cpl         (IM))
!     allocate (Coupling%ficein_cpl          (IM))
!     allocate (Coupling%hicein_cpl          (IM))
      allocate (Coupling%hsnoin_cpl          (IM))
!     allocate (Coupling%sfc_alb_nir_dir_cpl (IM))
!     allocate (Coupling%sfc_alb_nir_dif_cpl (IM))
!     allocate (Coupling%sfc_alb_vis_dir_cpl (IM))
!     allocate (Coupling%sfc_alb_vis_dif_cpl (IM))

      Coupling%slimskin_cpl          = clear_val
      Coupling%dusfcin_cpl           = clear_val
      Coupling%dvsfcin_cpl           = clear_val
      Coupling%dtsfcin_cpl           = clear_val
      Coupling%dqsfcin_cpl           = clear_val
      Coupling%ulwsfcin_cpl          = clear_val
!     Coupling%tseain_cpl            = clear_val
!     Coupling%tisfcin_cpl           = clear_val
!     Coupling%ficein_cpl            = clear_val
!     Coupling%hicein_cpl            = clear_val
      Coupling%hsnoin_cpl            = clear_val
!     Coupling%sfc_alb_nir_dir_cpl   = clear_val
!     Coupling%sfc_alb_nir_dif_cpl   = clear_val
!     Coupling%sfc_alb_vis_dir_cpl   = clear_val
!     Coupling%sfc_alb_vis_dif_cpl   = clear_val

      ! -- Coupling options to retrive atmosphere-ocean fluxes from mediator
      if (Model%use_med_flux) then
        allocate (Coupling%dusfcin_med (IM))
        allocate (Coupling%dvsfcin_med (IM))
        allocate (Coupling%dtsfcin_med (IM))
        allocate (Coupling%dqsfcin_med (IM))
        allocate (Coupling%ulwsfcin_med(IM))

        Coupling%dusfcin_med  = clear_val
        Coupling%dvsfcin_med  = clear_val
        Coupling%dtsfcin_med  = clear_val
        Coupling%dqsfcin_med  = clear_val
        Coupling%ulwsfcin_med = clear_val
      end if

      !--- accumulated quantities
      allocate (Coupling%dusfc_cpl  (IM))
      allocate (Coupling%dvsfc_cpl  (IM))
      allocate (Coupling%dtsfc_cpl  (IM))
      allocate (Coupling%dqsfc_cpl  (IM))
      allocate (Coupling%dnirbm_cpl (IM))
      allocate (Coupling%dnirdf_cpl (IM))
      allocate (Coupling%dvisbm_cpl (IM))
      allocate (Coupling%dvisdf_cpl (IM))
      allocate (Coupling%nlwsfc_cpl (IM))

      Coupling%dusfc_cpl  = clear_val
      Coupling%dvsfc_cpl  = clear_val
      Coupling%dtsfc_cpl  = clear_val
      Coupling%dqsfc_cpl  = clear_val
      Coupling%dnirbm_cpl = clear_val
      Coupling%dnirdf_cpl = clear_val
      Coupling%dvisbm_cpl = clear_val
      Coupling%dvisdf_cpl = clear_val
      Coupling%nlwsfc_cpl = clear_val

      !--- instantaneous quantities
      allocate (Coupling%dusfci_cpl  (IM))
      allocate (Coupling%dvsfci_cpl  (IM))
      allocate (Coupling%dnirbmi_cpl (IM))
      allocate (Coupling%dnirdfi_cpl (IM))
      allocate (Coupling%dvisbmi_cpl (IM))
      allocate (Coupling%dvisdfi_cpl (IM))
      allocate (Coupling%nlwsfci_cpl (IM))
      allocate (Coupling%oro_cpl     (IM))
      allocate (Coupling%slmsk_cpl   (IM))

      Coupling%dusfci_cpl  = clear_val
      Coupling%dvsfci_cpl  = clear_val
      Coupling%dnirbmi_cpl = clear_val
      Coupling%dnirdfi_cpl = clear_val
      Coupling%dvisbmi_cpl = clear_val
      Coupling%dvisdfi_cpl = clear_val
      Coupling%nlwsfci_cpl = clear_val
      Coupling%oro_cpl     = clear_val  !< pointer to sfcprop%oro
      Coupling%slmsk_cpl   = clear_val  !< pointer to sfcprop%slmsk
    endif

    ! -- Coupling options to retrive land fluxes from external land component 
    if (Model%cpllnd .and. Model%cpllnd2atm) then
      allocate (Coupling%sncovr1_lnd (IM))
      allocate (Coupling%qsurf_lnd   (IM))
      allocate (Coupling%evap_lnd    (IM))
      allocate (Coupling%hflx_lnd    (IM))
      allocate (Coupling%ep_lnd      (IM))
      allocate (Coupling%t2mmp_lnd   (IM))
      allocate (Coupling%q2mp_lnd    (IM))
      allocate (Coupling%gflux_lnd   (IM))
      allocate (Coupling%runoff_lnd  (IM))
      allocate (Coupling%drain_lnd   (IM))
      allocate (Coupling%cmm_lnd     (IM))
      allocate (Coupling%chh_lnd     (IM))
      allocate (Coupling%zvfun_lnd   (IM))

      Coupling%sncovr1_lnd = clear_val
      Coupling%qsurf_lnd   = clear_val
      Coupling%evap_lnd    = clear_val
      Coupling%hflx_lnd    = clear_val
      Coupling%ep_lnd      = clear_val
      Coupling%t2mmp_lnd   = clear_val
      Coupling%q2mp_lnd    = clear_val
      Coupling%gflux_lnd   = clear_val
      Coupling%runoff_lnd  = clear_val
      Coupling%drain_lnd   = clear_val
      Coupling%cmm_lnd     = clear_val
      Coupling%chh_lnd     = clear_val
      Coupling%zvfun_lnd   = clear_val
    end if

    !-- cellular automata
    allocate (Coupling%condition(IM))
    Coupling%condition = clear_val
    if (Model%do_ca) then
      allocate (Coupling%ca1      (IM))
      allocate (Coupling%ca2      (IM))
      allocate (Coupling%ca3      (IM))
      allocate (Coupling%ca_deep  (IM))
      allocate (Coupling%ca_turb  (IM))
      allocate (Coupling%ca_shal  (IM))
      allocate (Coupling%ca_rad   (IM))
      allocate (Coupling%ca_micro (IM))
      Coupling%ca1       = clear_val
      Coupling%ca2       = clear_val
      Coupling%ca3       = clear_val
      Coupling%ca_deep   = clear_val
      Coupling%ca_turb   = clear_val
      Coupling%ca_shal   = clear_val
      Coupling%ca_rad    = clear_val
      Coupling%ca_micro  = clear_val
    endif

    ! -- Aerosols coupling options
    if (Model%cplchm) then
      !--- outgoing instantaneous quantities
      allocate (Coupling%ushfsfci  (IM))
      ! -- instantaneous 3d fluxes of nonconvective ice and liquid precipitations
      allocate (Coupling%pfi_lsan  (IM,Model%levs))
      allocate (Coupling%pfl_lsan  (IM,Model%levs))
      Coupling%ushfsfci  = clear_val
      Coupling%pfi_lsan  = clear_val
      Coupling%pfl_lsan  = clear_val
    endif

    if (Model%cplchm .or. Model%cplflx .or. Model%cpllnd .or. Model%cpl_fire) then
      !--- accumulated convective rainfall
      allocate (Coupling%rainc_cpl (IM))
      Coupling%rainc_cpl = clear_val
    end if

    if(Model%progsigma)then
       allocate(Coupling%dqdt_qmicro (IM,Model%levs))
       Coupling%dqdt_qmicro = clear_val
    endif

    !--- stochastic physics option
    if (Model%do_sppt .or. Model%ca_global) then
      allocate (Coupling%sppt_wts  (IM,Model%levs))
      Coupling%sppt_wts = clear_val
    endif

    !--- stochastic shum option
    if (Model%do_shum) then
      allocate (Coupling%shum_wts  (IM,Model%levs))
      Coupling%shum_wts = clear_val
    endif

    !--- stochastic skeb option
    if (Model%do_skeb) then
      allocate (Coupling%skebu_wts (IM,Model%levs))
      allocate (Coupling%skebv_wts (IM,Model%levs))

      Coupling%skebu_wts = clear_val
      Coupling%skebv_wts = clear_val
    endif

    !--- stochastic land perturbation option
    if (Model%lndp_type /= 0) then
      allocate (Coupling%sfc_wts  (IM,Model%n_var_lndp))
      Coupling%sfc_wts = clear_val
    endif

    !--- stochastic spp perturbation option
    if (Model%do_spp) then
      allocate (Coupling%spp_wts_pbl  (IM,Model%levs))
      Coupling%spp_wts_pbl = clear_val
      allocate (Coupling%spp_wts_sfc  (IM,Model%levs))
      Coupling%spp_wts_sfc = clear_val
      allocate (Coupling%spp_wts_mp   (IM,Model%levs))
      Coupling%spp_wts_mp = clear_val
      allocate (Coupling%spp_wts_gwd   (IM,Model%levs))
      Coupling%spp_wts_gwd = clear_val
      allocate (Coupling%spp_wts_rad   (IM,Model%levs))
      Coupling%spp_wts_rad = clear_val
      allocate (Coupling%spp_wts_cu_deep   (IM,Model%levs))
      Coupling%spp_wts_cu_deep = clear_val
    endif

    !--- needed for Thompson's aerosol option
    if(Model%imp_physics == Model%imp_physics_thompson .and. (Model%ltaerosol .or. Model%mraerosol)) then
      allocate (Coupling%nwfa2d (IM))
      allocate (Coupling%nifa2d (IM))
      Coupling%nwfa2d   = clear_val
      Coupling%nifa2d   = clear_val
    endif

    if(Model%rrfs_sd) then
    !--- needed for smoke aerosol option
      allocate (Coupling%ebu_smoke (IM,Model%levs))
      allocate (Coupling%chem3d    (IM,Model%levs,Model%nchem))
      allocate (Coupling%ddvel     (IM,Model%ndvel))
      allocate (Coupling%wetdpc_flux(IM,Model%nchem))
      allocate (Coupling%wetdpr_flux(IM,Model%nchem))
      allocate (Coupling%drydep_flux(IM,Model%ndvel))
      allocate (Coupling%min_fplume(IM))
      allocate (Coupling%max_fplume(IM))
      allocate (Coupling%uspdavg(IM))
      allocate (Coupling%hpbl_thetav(IM))
      allocate (Coupling%rrfs_hwp  (IM))
      allocate (Coupling%rrfs_hwp_ave  (IM))
      Coupling%ebu_smoke  = clear_val
      Coupling%chem3d     = clear_val
      Coupling%ddvel      = clear_val
      Coupling%wetdpc_flux = clear_val
      Coupling%wetdpr_flux = clear_val
      Coupling%drydep_flux = clear_val
      Coupling%min_fplume = clear_val
      Coupling%max_fplume = clear_val
      Coupling%uspdavg = clear_val
      Coupling%hpbl_thetav = clear_val
      Coupling%rrfs_hwp   = clear_val
      Coupling%rrfs_hwp_ave = clear_val
    endif

    if (Model%imfdeepcnv == Model%imfdeepcnv_gf .or. Model%imfdeepcnv == Model%imfdeepcnv_c3) then
      allocate (Coupling%qci_conv (IM,Model%levs))
      Coupling%qci_conv   = clear_val
    endif

  end subroutine coupling_create


!----------------------
! GFS_control_type%init
!----------------------
  subroutine control_initialize (Model, nlunit, fn_nml, me, master, &
                                 logunit, isc, jsc, nx, ny, levs,   &
                                 cnx, cny, gnx, gny, dt_dycore,     &
                                 dt_phys, iau_offset, idat, jdat,   &
                                 nwat, tracer_names, tracer_types,  &
                                 input_nml_file, tile_num, blksz,   &
                                 ak, bk, restart, hydrostatic,      &
                                 communicator, ntasks, nthreads)

!--- modules
    use physcons,         only: con_rerth, con_pi
    use mersenne_twister, only: random_setseed, random_number
    use parse_tracers,    only: get_tracer_index
!
    implicit none

!--- interface variables
    class(GFS_control_type)            :: Model
    integer,                intent(in) :: nlunit
    character(len=64),      intent(in) :: fn_nml
    integer,                intent(in) :: me
    integer,                intent(in) :: master
    integer,                intent(in) :: logunit
    integer,                intent(in) :: tile_num
    integer,                intent(in) :: isc
    integer,                intent(in) :: jsc
    integer,                intent(in) :: nx
    integer,                intent(in) :: ny
    integer,                intent(in) :: levs
    integer,                intent(in) :: cnx
    integer,                intent(in) :: cny
    integer,                intent(in) :: gnx
    integer,                intent(in) :: gny
    real(kind=kind_phys),   intent(in) :: dt_dycore
    real(kind=kind_phys),   intent(in) :: dt_phys
    integer,                intent(in) :: iau_offset
    integer,                intent(in) :: idat(8)
    integer,                intent(in) :: jdat(8)
    integer,                intent(in) :: nwat
    character(len=32),      intent(in) :: tracer_names(:)
    integer,                intent(in) :: tracer_types(:)
    character(len=:),       intent(in), dimension(:), pointer :: input_nml_file
    integer,                intent(in) :: blksz(:)
    real(kind=kind_phys), dimension(:), intent(in) :: ak
    real(kind=kind_phys), dimension(:), intent(in) :: bk
    logical,                intent(in) :: restart
    logical,                intent(in) :: hydrostatic
    type(MPI_Comm),         intent(in) :: communicator
    integer,                intent(in) :: ntasks
    integer,                intent(in) :: nthreads

    !--- local variables
    integer :: i, j, n
    integer :: ios
    integer :: seed0
    logical :: exists
    real(kind=kind_phys) :: tem
    real(kind=kind_dbl_prec) :: rinc(5)
    real(kind=kind_phys) :: wrk(1)
    real(kind=kind_phys), parameter :: con_hr = 3600.

!--- BEGIN NAMELIST VARIABLES
    real(kind=kind_phys) :: fhzero         = 0.0             !< hours between clearing of diagnostic buckets
    real(kind=kind_phys) :: fhzero_array(1:2)  = 0.0         !< array with hours between clearing of diagnostic buckets
    real(kind=kind_phys) :: fhzero_fhour(1:2)  = 0.0         !< the maximum forecast length for the hours between clearing of diagnostic buckets
    logical              :: ldiag3d        = .false.         !< flag for 3d diagnostic fields
    logical              :: qdiag3d        = .false.         !< flag for 3d tracer diagnostic fields
    logical              :: lssav          = .false.         !< logical flag for storing diagnostics
    integer, parameter :: pat_len = 60, pat_count=100        !< dimensions of dtend_select
    character(len=pat_len) :: dtend_select(pat_count)         !< fglob_list() patterns to decide which 3d diagnostic fields to enable
    integer              :: naux2d         = 0               !< number of auxiliary 2d arrays to output (for debugging)
    integer              :: naux3d         = 0               !< number of auxiliary 3d arrays to output (for debugging)
    logical              :: aux2d_time_avg(1:naux2dmax) = .false. !< flags for time averaging of auxiliary 2d arrays
    logical              :: aux3d_time_avg(1:naux3dmax) = .false. !< flags for time averaging of auxiliary 3d arrays

    real(kind=kind_phys) :: fhcyc          = 0.              !< frequency for surface data cycling (hours)
    integer              :: thermodyn_id   =  1              !< valid for GFS only for get_prs/phi
    integer              :: sfcpress_id    =  1              !< valid for GFS only for get_prs/phi

    !--- coupling parameters
    logical              :: cplflx         = .false.         !< default no cplflx collection
    logical              :: cplice         = .false.         !< default no cplice collection (used together with cplflx)
    logical              :: cplocn2atm     = .true.          !< default yes cplocn2atm coupling (turn on the feedback from ocn to atm)
    logical              :: cplwav         = .false.         !< default no cplwav collection
    logical              :: cplwav2atm     = .false.         !< default no cplwav2atm coupling
    logical              :: cplaqm         = .false.         !< default no cplaqm collection
    logical              :: cplchm         = .false.         !< default no cplchm collection
    logical              :: cpllnd         = .false.         !< default no cpllnd collection
    logical              :: cpllnd2atm     = .false.         !< default no cpllnd2atm coupling
    logical              :: rrfs_sd        = .false.         !< default no rrfs_sd collection
    logical              :: cpl_fire       = .false.         !< default no fire behavior colleciton
    logical              :: use_cice_alb   = .false.         !< default no cice albedo
    logical              :: cpl_imp_mrg    = .false.         !< default no merge import with internal forcings
    logical              :: cpl_imp_dbg    = .false.         !< default no write import data to file post merge
    logical              :: use_med_flux   = .false.         !< default no atmosphere-ocean fluxes from mediator

!--- integrated dynamics through earth's atmosphere
    logical              :: lsidea         = .false.

!--- radiation parameters
    real(kind=kind_phys) :: fhswr          = 3600.           !< frequency for shortwave radiation (secs)
    real(kind=kind_phys) :: fhlwr          = 3600.           !< frequency for longwave radiation (secs)
    integer              :: nhfrad         = 0               !< number of timesteps for which to call radiation on physics timestep (coldstarts)
    integer              :: levr           = -99             !< number of vertical levels for radiation calculations
    integer              :: nfxr           = 39+6            !< second dimension of input/output array fluxr
    logical              :: iaerclm        = .false.         !< flag for initializing aero data
    integer              :: iccn           =  0              !< logical to use IN CCN forcing for MG2/3
    integer              :: iflip          =  1              !< iflip - is not the same as flipv
    integer              :: isol           =  0              !< use prescribed solar constant
                                                             !< 0  => fixed value=1366.0\f$W/m^2\f$(old standard)
                                                             !< 10 => fixed value=1360.8\f$W/m^2\f$(new standard)
                                                             !< 1  => NOAA ABS-scale TSI table (yearly) w 11-yr cycle approx
                                                             !< 2  => NOAA TIM-scale TSI table (yearly) w 11-yr cycle approx
                                                             !< 3  => CMIP5 TIM-scale TSI table (yearly) w 11-yr cycle approx
                                                             !< 4  => CMIP5 TIM-scale TSI table (monthly) w 11-yr cycle approx
    integer              :: ico2           =  0              !< prescribed global mean value (old opernl)
    integer              :: ialb           =  0              !< use climatology alb, based on sfc type
                                                             !< 1 => use modis based alb (RUC lsm)
                                                             !< 2 => use LSM albedo (Noah MP lsm)
    integer              :: iems           =  0              !< 1.0 => Noah lsm
                                                             !< 2.0 => Noah MP and RUC lsms
    integer              :: iaer           =  1              !< default aerosol effect in sw only
    integer              :: iaermdl        =  0              !< default tropospheric aerosol model scheme flag
                                                             !< 0: seasonal global distributed OPAC aerosol climatology
                                                             !< 1: monthly global distributed GOCART aerosol climatology
                                                             !< 2: GOCART prognostic aerosol model
                                                             !< 5: OPAC climatoloy with new band mapping
    integer              :: iaerflg        =  0              !< aerosol effect control flag
                                                             !< 3-digit flag 'abc':
                                                             !< a-stratospheric volcanic aerols
                                                             !< b-tropospheric aerosols for LW
                                                             !< c-tropospheric aerosols for SW
                                                             !<  =0:aerosol effect is not included; =1:aerosol effect is included
    logical              :: lalw1bd        = .false.         !< selects 1 band or multi bands for LW aerosol properties
                                                             !< true.:  aerosol properties calculated in 1 broad LW band
                                                             !< false.: aerosol properties calculated for each LW bands
    character(len=26)    :: aeros_file     = 'aerosol.dat               '
    character(len=26)    :: solar_file     = 'solarconstant_noaa_a0.txt '
    character(len=26)    :: semis_file     = 'sfc_emissivity_idx.txt    '
    character(len=26)    :: co2dat_file    = 'co2historicaldata_2004.txt'
    character(len=26)    :: co2gbl_file    = 'co2historicaldata_glob.txt'
    character(len=26)    :: co2usr_file    = 'co2userdata.txt           '
    character(len=26)    :: co2cyc_file    = 'co2monthlycyc.txt         '
    integer              :: icliq_sw       =  1              !< sw optical property for liquid clouds
    integer              :: icice_sw       =  3              !< sw optical property for ice clouds
    integer              :: icliq_lw       =  1              !< lw optical property for liquid clouds
    integer              :: icice_lw       =  3              !< lw optical property for ice clouds
    integer              :: iovr           =  1              !< cloud-overlap used in cloud-sampling by radiation scheme(s)
    integer              :: ictm           =  1              !< ictm=0 => use data at initial cond time, if not
                                                             !<           available; use latest; no extrapolation.
                                                             !< ictm=1 => use data at the forecast time, if not
                                                             !<           available; use latest; do extrapolation.
                                                             !< ictm=yyyy0 => use yyyy data for the forecast time;
                                                             !<           no extrapolation.
                                                             !< ictm=yyyy1 = > use yyyy data for the fcst. If needed,
                                                             !<           do extrapolation to match the fcst time.
                                                             !< ictm=-1 => use user provided external data for
                                                             !<           the fcst time; no extrapolation.
                                                             !< ictm=-2 => same as ictm=0, but add seasonal cycle
                                                             !<           from climatology; no extrapolation.
    integer              :: isubc_sw          =  0           !< sw clouds without sub-grid approximation
    integer              :: isubc_lw          =  0           !< lw clouds without sub-grid approximation
                                                             !< =1 => sub-grid cloud with prescribed seeds
                                                             !< =2 => sub-grid cloud with randomly generated
                                                             !< seeds
    integer              :: iswmode           =  2           !< SW control flag for scattering process approximation
                                                             !< =1 => two-stream delta-eddington (Joseph et al. 1976)
                                                             !< =2 => two-stream PIFM            (Zdunkowski et al. 1980)
                                                             !< =3 => discrete ordinates         (Liou, 1973)
    integer              :: idcor = 1                        !< Decorrelation length type for overlap assumption
                                                             !< =0 => Use constant decorrelation length, decorr_con
                                                             !< =1 => Use spatially varying decorrelation length (Hogan et al. 2010)
                                                             !< =2 => Use spatially and temporally varyint decorrelation length (Oreopoulos et al. 2012)
    real(kind_phys)      :: dcorr_con         = 2.5          !< Decorrelation length constant (km) (if idcor = 0)
    logical              :: lcrick            = .false.      !< CRICK-Proof cloud water
    logical              :: lcnorm            = .false.      !< Cloud condensate normalized by cloud cover
    logical              :: lnoprec           = .false.      !< radiation precip flag for Ferrier/Moorthi
    logical              :: lwhtr             = .true.       !< flag to output lw heating rate (Radtend%lwhc)
    logical              :: swhtr             = .true.       !< flag to output sw heating rate (Radtend%swhc)
    integer              :: rad_hr_units      = 2            !< heating rate units are K s-1
    logical              :: inc_minor_gas     = .true.       !< Include minor trace gases in RRTMG radiation calculation
    integer              :: ipsd0             = 0            !< initial permutaion seed for mcica radiation
    integer              :: ipsdlim           = 1e8          !< limit initial permutaion seed for mcica radiation
    logical              :: lrseeds           = .false.      !< flag to use host-provided random seeds
    integer              :: nrstreams         = 2            !< number of random number streams in host-provided random seed array
    logical              :: lextop            = .false.      !< flag for using an extra top layer for radiation
    ! RRTMGP
    logical              :: do_RRTMGP           = .false.    !< Use RRTMGP?
    character(len=128)   :: active_gases        = ''         !< Character list of active gases used in RRTMGP
    integer              :: nGases              = 0          !< Number of active gases
    character(len=128)   :: rrtmgp_root         = ''         !< Directory of rte+rrtmgp source code
    character(len=128)   :: lw_file_gas         = ''         !< RRTMGP K-distribution file, coefficients to compute optics for gaseous atmosphere
    character(len=128)   :: lw_file_clouds      = ''         !< RRTMGP file containing coefficients used to compute clouds optical properties
    character(len=128)   :: sw_file_gas         = ''         !< RRTMGP K-distribution file, coefficients to compute optics for gaseous atmosphere
    character(len=128)   :: sw_file_clouds      = ''         !< RRTMGP file containing coefficients used to compute clouds optical properties
    integer              :: rrtmgp_nBandsSW     = -999       !< Number of RRTMGP SW bands.             # *NOTE*
    integer              :: rrtmgp_nGptsSW      = -999       !< Number of RRTMGP SW spectral points.   # The RRTMGP spectral dimensions in the files
    integer              :: rrtmgp_nBandsLW     = -999       !< Number of RRTMGP LW bands.             # need to be provided via namelsit.
    integer              :: rrtmgp_nGptsLW      = -999       !< Number of RRTMGP LW spectral points.   #
    logical              :: doG_cldoptics       = .false.    !< Use legacy RRTMG cloud-optics?
    logical              :: doGP_cldoptics_PADE = .false.    !< Use RRTMGP cloud-optics: PADE approximation?
    logical              :: doGP_cldoptics_LUT  = .false.    !< Use RRTMGP cloud-optics: LUTs?
    integer              :: iovr_convcld        = 1          !< Cloud-overlap assumption for convective-cloud (defaults to iovr if not set)
    integer              :: rrtmgp_nrghice      = 3          !< Number of ice-roughness categories
    integer              :: rrtmgp_nGauss_ang   = 1          !< Number of angles used in Gaussian quadrature
    logical              :: do_GPsw_Glw         = .false.
    logical              :: use_LW_jacobian     = .false.    !< Use Jacobian of LW to update LW radiation tendencies.
    logical              :: damp_LW_fluxadj     = .false.    !< Damp LW Jacobian flux adjustment with height.
    real(kind=kind_phys) :: lfnc_k              = -999       !<
    real(kind=kind_phys) :: lfnc_p0             = -999       !<
    logical              :: doGP_lwscat         = .false.    !< If true, include scattering in longwave cloud-optics, only compatible w/ GP cloud-optics
    logical              :: doGP_sgs_cnv        = .false.    !< If true, include SubGridScale convective cloud in RRTMGP
    logical              :: doGP_sgs_mynn       = .false.    !< If true, include SubGridScale MYNN-EDMF cloud in RRTMGP
    integer              :: rrtmgp_lw_phys_blksz= 1          !< Number of columns for RRTMGP LW scheme to process at each instance.
    integer              :: rrtmgp_sw_phys_blksz= 1          !< Number of columns for RRTMGP SW scheme to process at each instance.
    logical              :: doGP_smearclds      = .true.     !< If true, include implicit SubGridScale clouds in RRTMGP
!--- Z-C microphysical parameters
    integer              :: imp_physics       =  99                !< choice of cloud scheme
    real(kind=kind_phys) :: psautco(2)        = (/6.0d-4,3.0d-4/)  !< [in] auto conversion coeff from ice to snow
    real(kind=kind_phys) :: prautco(2)        = (/1.0d-4,1.0d-4/)  !< [in] auto conversion coeff from cloud to rain
    real(kind=kind_phys) :: evpco             = 2.0d-5             !< [in] coeff for evaporation of largescale rain
    real(kind=kind_phys) :: wminco(2)         = (/1.0d-5,1.0d-5/)  !< [in] water and ice minimum threshold for Zhao
!---Max hourly
    real(kind=kind_phys) :: avg_max_length = 3600.                 !< reset value in seconds for max hourly
!--- Ferrier-Aligo microphysical parameters
    real(kind=kind_phys) :: rhgrd             = 1.0                !< fer_hires microphysics only; for 3-km domain
    logical              :: spec_adv          = .true.             !< Individual cloud species advected
    integer              :: icloud            = 0                  !< cloud effect to the optical depth in radiation; this also controls the cloud fraction options
                                                                   !<  3: with cloud effect from FA, and use cloud fraction option 3, based on Sundqvist et al. (1989)
!--- M-G microphysical parameters
    integer              :: fprcp             =  2                 !< when "0" no prognostic rain and snow (MG)
                                                                   !< "1" for MG2 and "2" for MG3
    integer              :: pdfflag           =  4                 !< pdf flag for MG macro physics
    real(kind=kind_phys) :: mg_dcs            = 200.0              !< Morrison-Gettelman microphysics parameters
    real(kind=kind_phys) :: mg_qcvar          = 1.0
    real(kind=kind_phys) :: mg_ts_auto_ice(2) = (/180.0,180.0/)    !< ice auto conversion time scale
    real(kind=kind_phys) :: mg_rhmini         = 1.01               !< relative humidity threshold parameter for nucleating ice
    real(kind=kind_phys) :: mg_ncnst          = 100.e6             !< constant droplet num concentration (m-3)
    real(kind=kind_phys) :: mg_ninst          = 0.15e6             !< constant ice num concentration (m-3)
    real(kind=kind_phys) :: mg_ngnst          = 0.10e6             !< constant graupel/hail num concentration (m-3) = 0.1e6_r8
    real(kind=kind_phys) :: mg_alf            = 1.0                !< tuning factor for alphs in MG macrophysics
    real(kind=kind_phys) :: mg_qcmin(2)       = (/1.0d-9,1.0d-9/)  !< min liquid and ice mixing ratio in Mg macro clouds
    real(kind=kind_phys) :: mg_berg_eff_factor = 2.0               !< berg efficiency factor
    character(len=16)    :: mg_precip_frac_method = 'max_overlap'  !< type of precipitation fraction method
    real(kind=kind_phys) :: tf              = 258.16d0
    real(kind=kind_phys) :: tcr             = 273.16d0
!
    logical              :: effr_in         = .false.              !< flag to use effective radii of cloud species in radiation
    logical              :: microp_uniform  = .true.
    logical              :: do_cldliq       = .true.
    logical              :: do_cldice       = .true.
    logical              :: hetfrz_classnuc = .false.
    logical              :: mg_nccons       = .false.           !< set .true. to specify constant cloud droplet number
    logical              :: mg_nicons       = .false.           !< set .true. to specify constant cloud ice number
    logical              :: mg_ngcons       = .false.           !< set .true. to specify constant graupel/hail number
    logical              :: sed_supersat    = .true.
    logical              :: do_sb_physics   = .true.
    logical              :: mg_do_graupel   = .true.            !< set .true. to turn on prognostic grapuel (with fprcp=2)
    logical              :: mg_do_hail      = .false.           !< set .true. to turn on prognostic hail (with fprcp=2)
    logical              :: mg_do_ice_gmao  = .false.           !< set .true. to turn on gmao ice formulation
    logical              :: mg_do_liq_liu   = .true.            !< set .true. to turn on liu liquid treatment
    real(kind=kind_phys) :: fh_dfi_radar(1+dfi_radar_max_intervals) = -2e10             !< begin&end of four timespans over which radar_tten is applied
    logical              :: do_cap_suppress = .true.            !< set .true. to turn on convection suppression in GF scheme during limited intervals when fh_dfi_radar is enabled

    !--- NSSL microphysics params
    real(kind=kind_phys) :: nssl_cccn       = 0.6e9             !<  CCN concentration (m-3)
    real(kind=kind_phys) :: nssl_alphah     = 0.0               !<  graupel shape parameter
    real(kind=kind_phys) :: nssl_alphahl    = 1.0               !<  hail shape parameter
    real(kind=kind_phys) :: nssl_alphar     = 0.0               ! shape parameter for rain (imurain=1 only)
    real(kind=kind_phys) :: nssl_ehw0       = 0.9               ! constant or max assumed graupel-droplet collection efficiency
    real(kind=kind_phys) :: nssl_ehlw0      = 0.9               ! constant or max assumed hail-droplet collection efficiency
    logical              :: nssl_hail_on    = .false.           !<  NSSL flag to activate the hail category
    logical              :: nssl_ccn_on     = .true.            !<  NSSL flag to activate the CCN category
    logical              :: nssl_invertccn  = .true.            !<  NSSL flag to treat CCN as activated (true) or unactivated (false)
    logical              :: nssl_3moment    = .false.           !<  NSSL flag to turn on 3-moment for rain/graupel/hail

    !--- Thompson microphysical parameters
    logical              :: ltaerosol      = .false.            !< flag for aerosol version
    logical              :: mraerosol      = .false.            !< flag for merra2_aerosol_aware
    logical              :: lradar         = .false.            !< flag for radar reflectivity
    real(kind=kind_phys) :: nsfullradar_diag  = -999.0          !< seconds between resetting radar reflectivity calculation, set to <0 for every time step
    real(kind=kind_phys) :: ttendlim       = -999.0             !< temperature tendency limiter, set to <0 to deactivate
    logical              :: ext_diag_thompson = .false.         !< flag for extended diagnostic output from Thompson
    real(kind=kind_phys) :: dt_inner       = -999.0             !< time step for the inner loop
    logical              :: sedi_semi      = .false.            !< flag for semi Lagrangian sedi of rain
    integer              :: decfl          = 8                  !< deformed CFL factor
    logical              :: thpsnmp_is_init = .false.           !< Local scheme initialization flag 

    !--- GFDL microphysical parameters
    logical              :: lgfdlmprad     = .false.            !< flag for GFDLMP radiation interaction

    !--- Thompson,GFDL microphysical parameter
    logical              :: lrefres        = .false.            !< flag for radar reflectivity in restart file

    !--- CLM Lake Model parameters (MUST match clm_lake.F90)
    integer, parameter   :: nlevlake_clm_lake = 10           !< number of lake levels
    integer, parameter   :: nlevsoil_clm_lake = 10           !< number of soil levels
    integer, parameter   :: nlevsnow_clm_lake = 5            !< number of snow levels
    integer, parameter   :: nlevsnowsoil_clm_lake = nlevsnow_clm_lake+nlevsoil_clm_lake+1 !< -nlevsno:nlevsoil dimensioned variables
    integer, parameter   :: nlevsnowsoil1_clm_lake = nlevsnow_clm_lake+nlevsoil_clm_lake !< -nlevsno+1:nlevsoil dimensioned variables

    !--- CLM Lake configurables
    real(kind_phys)      :: clm_lake_depth_default = 50         !< default lake depth in clm lake model
    logical              :: clm_lake_use_lakedepth = .true.     !< initialize depth from lakedepth
    logical              :: clm_lake_debug = .false.            !< verbose debugging in clm_lake
    logical              :: clm_debug_print = .false.           !< enables prints in clm_lake

    !--- land/surface model parameters
    integer              :: lsm            =  1              !< flag for land surface model to use =0  for osu lsm; =1  for noah lsm; =2  for noah mp lsm; =3  for RUC lsm
    integer              :: lsoil          =  4              !< number of soil layers
    integer              :: lsoil_lsm      =  -1             !< number of soil layers internal to land surface model; -1 use lsoil
    integer              :: lsnow_lsm      =  3              !< maximum number of snow layers internal to land surface model
    logical              :: exticeden      = .false.         !< Use variable precip ice density for NOAH LSM if true or original formulation
    logical              :: rdlai          = .false.         !< read LAI from input file at cold start (for RUC LSM or NOAH LSM WRFv4)
    logical              :: ua_phys        = .false.         !< flag for using University of Arizona? extension to NOAH LSM WRFv4
    logical              :: usemonalb      = .true.          !< flag to read surface diffused shortwave albedo from input file for NOAH LSM WRFv4
    real(kind=kind_phys) :: aoasis         = 1.0             !< potential evaporation multiplication factor for NOAH LSM WRFv4
    integer              :: fasdas         = 0               !< flag to use "flux-adjusting surface data assimilation system"; 0 = OFF, 1 = ON
    integer              :: iopt_thcnd     = 1               !< option to treat thermal conductivity in Noah LSM (new in 3.8)
                                                             !< = 1, original (default)
                                                             !< = 2, McCumber and Pielke for silt loam and sandy loam
    integer              :: kice           =  2              !< number of layers in ice; default is 2 (GFS sice)
    integer              :: ivegsrc        =  2              !< ivegsrc = 0   => USGS,
                                                             !< ivegsrc = 1   => IGBP (20 category)
                                                             !< ivegsrc = 2   => UMD  (13 category)
    integer              :: nvegcat        =  20             !< number of veg.  categories depending on ivegsrc
    integer              :: isot           =  0              !< isot = 0   => Zobler soil type  ( 9 category)
                                                             !< isot = 1   => STATSGO soil type (19 category)
    integer              :: nsoilcat       =  16             !< number of soil categories depending on isot
    ! -- to use Noah MP, lsm needs to be set to 2 and both ivegsrc and isot are set
    ! to 1 - MODIS IGBP and STATSGO - the defaults are the same as in the
    ! scripts;change from namelist

    integer              :: iopt_dveg      =  4  ! 4 -> off (use table lai; use maximum vegetation fraction)
    integer              :: iopt_crs       =  1  !canopy stomatal resistance (1-> ball-berry; 2->jarvis)
    integer              :: iopt_btr       =  1  !soil moisture factor for stomatal resistance (1-> noah; 2-> clm; 3-> ssib)
    integer              :: iopt_run       =  3  !runoff and groundwater (1->simgm; 2->simtop; 3->schaake96; 4->bats)
    integer              :: iopt_sfc       =  1  !surface layer drag coeff (ch & cm) (1->m-o; 2->chen97)
    integer              :: iopt_frz       =  1  !supercooled liquid water (1-> ny06; 2->koren99)
    integer              :: iopt_inf       =  1  !frozen soil permeability (1-> ny06; 2->koren99)
    integer              :: iopt_rad       =  3  !radiation transfer (1->gap=f(3d,cosz); 2->gap=0; 3->gap=1-fveg)
    integer              :: iopt_alb       =  2  !snow surface albedo (1->bats; 2->class)
    integer              :: iopt_snf       =  1  !rainfall & snowfall (1-jordan91; 2->bats; 3->noah)
    integer              :: iopt_tbot      =  2  !lower boundary of soil temperature (1->zero-flux; 2->noah)
    integer              :: iopt_stc       =  1  !snow/soil temperature time scheme (only layer 1)
    integer              :: iopt_trs       =  2  !thermal roughness scheme (1-z0h=z0m; 2-czil; 3-ec;4-kb reversed)
    integer              :: iopt_diag      =  2  !2m t/q diagnostic approach (1->external GFS sfc_diag 2->original NoahMP 2-title
                                                 !3->NoahMP 2-title + internal GFS sfc_diag  )

    integer              :: mosaic_lu      =  0  ! 1 - used of fractional landuse in RUC lsm
    integer              :: mosaic_soil    =  0  ! 1 - used of fractional soil in RUC lsm
    integer              :: isncond_opt    =  1  ! 2 - Sturm (1997)
    integer              :: isncovr_opt    =  1  ! 2 - Niu-Yang (2007), 3-updated Niu-Yang similar to Noah MP

    logical              :: use_ufo        = .false.                  !< flag for gcycle surface option

    logical              :: add_fire_heat_flux = .false.              !< Flag for fire heat flux

    logical              :: lcurr_sf       = .false.                  !< flag for taking ocean currents into account in GFDL surface layer
    logical              :: pert_cd        = .false.                  !< flag for perturbing the surface drag coefficient for momentum in surface layer scheme
    integer              :: ntsflg         = 0                        !< flag for updating skin temperature in the GFDL surface layer scheme
    real(kind=kind_phys) :: sfenth         = 0.0                      !< enthalpy flux factor 0 zot via charnock ..>0 zot enhanced>15m/s

!--- flake model parameters
    integer              :: lkm            =  0                       !< =1 run lake, =2 run lake&nsst =0 no lake
    integer              :: iopt_lake      =  2                       !< =1 flake, =2 clm lake (default)
    real(kind_phys)      :: lakedepth_threshold = 1.0                 !< lakedepth must be GREATER than this value to enable a lake model
    real(kind_phys)      :: lakefrac_threshold  = 0.0                 !< lakefrac must be GREATER than this value to enable a lake model
    logical              :: use_lake2m     = .false.                  !< use 2m T & Q from clm lake model

!--- tuning parameters for physical parameterizations
    logical              :: ras            = .false.                  !< flag for ras convection scheme
    logical              :: flipv          = .true.                   !< flag for vertical direction flip (ras)
                                                                      !< .true. implies surface at k=1
    logical              :: trans_trac     = .false.                  !< flag for convective transport of tracers (RAS, CS, or SAMF)
    logical              :: old_monin      = .false.                  !< flag for diff monin schemes
    logical              :: cnvgwd         = .false.                  !< flag for conv gravity wave drag
    integer              :: gwd_opt        =  1                       !< flag for configuring gwd scheme
                                                                      !< gwd_opt = 2  => unified ugwp GWD
                                                                      !< gwd_opt = 22 => unified ugwp GWD with extra output
                                                                      !< gwd_opt = 3 : GSL drag suite
                                                                      !< gwd_opt = 33: GSL drag suite with extra output
    logical              :: do_ugwp_v0           = .true.       !< flag for version 0 ugwp GWD
    logical              :: do_ugwp_v0_orog_only = .false.      !< flag for version 0 ugwp GWD (orographic drag only)
    logical              :: do_ugwp_v0_nst_only  = .false.      !< flag for version 0 ugwp GWD (non-stationary GWD only)
    logical              :: do_gsl_drag_ls_bl    = .false.      !< flag for GSL drag (mesoscale GWD and blocking only)
    logical              :: do_gsl_drag_ss       = .false.      !< flag for GSL drag (small-scale GWD only)
    logical              :: do_gsl_drag_tofd     = .false.      !< flag for GSL drag (turbulent orog form drag only)
    logical              :: do_gwd_opt_psl       = .false.      !< flag for PSL drag (mesoscale GWD and blocking only)
    logical              :: do_ugwp_v1           = .false.      !< flag for version 1 ugwp GWD
    logical              :: do_ugwp_v1_orog_only = .false.      !< flag for version 1 ugwp GWD (orographic drag only)
    logical              :: do_ugwp_v1_w_gsldrag = .false.      !< flag for version 1 ugwp GWD (orographic drag only)
!--- vay-2018
    logical              :: ldiag_ugwp      = .false.                 !< flag for UGWP diag fields
    logical              :: ugwp_seq_update = .false.                 !< flag for updating winds between UGWP steps
    logical              :: do_ugwp         = .false.                 !< flag do UGWP+RF
    logical              :: do_tofd         = .false.                 !< flag do Turb oro Form Drag

    logical              :: do_gwd         = .false.                  !< flag for running gravity wave drag
    logical              :: do_cnvgwd      = .false.                  !< flag for running conv gravity wave drag

    logical              :: mstrat         = .false.                  !< flag for moorthi approach for stratus
    logical              :: moist_adj      = .false.                  !< flag for moist convective adjustment
    logical              :: cscnv          = .false.                  !< flag for Chikira-Sugiyama convection
    logical              :: cal_pre        = .false.                  !< flag controls precip type algorithm
    logical              :: do_aw          = .false.                  !< AW scale-aware option in cs convection
    logical              :: do_awdd        = .false.                  !< AW scale-aware option in cs convection
    logical              :: flx_form       = .false.                  !< AW scale-aware option in cs convection
    logical              :: do_shoc        = .false.                  !< flag for SHOC
    logical              :: shocaftcnv     = .false.                  !< flag for SHOC
    logical              :: shoc_cld       = .false.                  !< flag for SHOC in grrad
    logical              :: h2o_phys       = .false.                  !< flag for stratosphere h2o
    logical              :: pdfcld         = .false.                  !< flag for pdfcld
    logical              :: shcnvcw        = .false.                  !< flag for shallow convective cloud
    logical              :: redrag         = .false.                  !< flag for reduced drag coeff. over sea
    logical              :: hybedmf        = .false.                  !< flag for hybrid edmf pbl scheme
    logical              :: satmedmf       = .false.                  !< flag for scale-aware TKE-based moist edmf
                                                                      !< vertical turbulent mixing scheme
    logical              :: shinhong       = .false.                  !< flag for scale-aware Shinhong vertical turbulent mixing scheme
    logical              :: do_ysu         = .false.                  !< flag for YSU vertical turbulent mixing scheme
    logical              :: dspheat        = .false.                  !< flag for tke dissipative heating
    logical              :: hurr_pbl       = .false.                  !< flag for hurricane-specific options in PBL scheme
    logical              :: lheatstrg      = .false.                  !< flag for canopy heat storage parameterization
    logical              :: lseaspray      = .false.                  !< flag for sea spray parameterization
    logical              :: cnvcld         = .false.
    logical              :: xr_cnvcld      = .true.                   !< flag for including suspended convective clouds in Xu-Randall cloud fraction
    logical              :: random_clds    = .false.                  !< flag controls whether clouds are random
    logical              :: shal_cnv       = .false.                  !< flag for calling shallow convection
    integer              :: imfshalcnv     =  1                       !< flag for mass-flux shallow convection scheme
                                                                      !<     1: July 2010 version of mass-flux shallow conv scheme
                                                                      !<         current operational version as of 2016
                                                                      !<     2: scale- & aerosol-aware mass-flux shallow conv scheme (2017)
                                                                      !<     3: scale- & aerosol-aware Grell-Freitas scheme (GSD)
                                                                      !<     4: New Tiedtke scheme (CAPS)
                                                                      !<     0: modified Tiedtke's eddy-diffusion shallow conv scheme
                                                                      !<    -1: no shallow convection used
    integer              :: imfdeepcnv     =  1                       !< flag for mass-flux deep convection scheme
                                                                      !<     1: July 2010 version of SAS conv scheme
                                                                      !<           current operational version as of 2016
                                                                      !<     2: scale- & aerosol-aware mass-flux deep conv scheme (2017)
                                                                      !<     3: scale- & aerosol-aware Grell-Freitas scheme (GSD)
                                                                      !<     4: New Tiedtke scheme (CAPS)
    integer              :: isatmedmf      =  0                       !< flag for scale-aware TKE-based moist edmf scheme
                                                                      !<     0: initial version of satmedmf (Nov. 2018)
                                                                      !<     1: updated version of satmedmf (as of May 2019)
    logical              :: do_deep        = .true.                   !< whether to do deep convection

    logical              :: hwrf_samfdeep     = .false.               !< flag for HWRF SAMF deepcnv scheme
    logical              :: hwrf_samfshal     = .false.               !< flag for HWRF SAMF shalcnv scheme
    logical              :: progsigma         = .false.               !< flag for prognostic updraft area fraction closure in saSAS or Unified conv.
    logical              :: do_mynnedmf       = .false.               !< flag for MYNN-EDMF
    logical              :: do_mynnsfclay     = .false.               !< flag for MYNN Surface Layer Scheme
    ! DH* TODO - move to MYNN namelist section
    integer              :: tke_budget        = 0
    logical              :: bl_mynn_tkeadvect = .false.
    integer              :: bl_mynn_cloudpdf  = 2
    integer              :: bl_mynn_mixlength = 1
    integer              :: bl_mynn_edmf      = 1
    integer              :: bl_mynn_edmf_mom  = 1
    integer              :: bl_mynn_edmf_tke  = 0
    integer              :: bl_mynn_cloudmix  = 1
    integer              :: bl_mynn_mixqt     = 0
    integer              :: bl_mynn_output    = 0
    integer              :: icloud_bl         = 1
    real(kind=kind_phys) :: bl_mynn_closure   = 2.6                   !<   <= 2.5  only prognose tke
                                                                      !<   2.5 < and < 3.0, prognose tke and q'2
                                                                      !<   >= 3.0, prognose tke, q'2, T'2, and T'q'
    logical              :: sfclay_compute_diag = .false.
    logical              :: sfclay_compute_flux = .false.
    integer              :: isftcflx          = 0
    integer              :: iz0tlnd           = 0
    real(kind=kind_phys) :: var_ric           = 1.0
    real(kind=kind_phys) :: coef_ric_l        = 0.16
    real(kind=kind_phys) :: coef_ric_s        = 0.25
    !Prognostic convective closure
    real(kind=kind_phys) :: betascu           = 8.0 !< Tuning parameter for prog. closure shallow clouds
    real(kind=kind_phys) :: betamcu           = 1.0 !< Tuning parameter for prog. closure midlevel clouds
    real(kind=kind_phys) :: betadcu           = 2.0 !< Tuning parameter for prog. closure deep clouds
    ! *DH
    logical              :: do_myjsfc         = .false.               !< flag for MYJ surface layer scheme
    logical              :: do_myjpbl         = .false.               !< flag for MYJ PBL scheme

    integer              :: nmtvr          = 14                       !< number of topographic variables such as variance etc
                                                                      !< used in the GWD parameterization
    integer              :: jcap           =  1              !< number of spectral wave trancation used only by sascnv shalcnv
!   real(kind=kind_phys) :: cs_parm(10) = (/5.0,2.5,1.0e3,3.0e3,20.0,-999.,-999.,0.,0.,0./)
    real(kind=kind_phys) :: cs_parm(10) = (/8.0,4.0,1.0e3,3.5e3,20.0,1.0,-999.,1.,0.6,0./)
    real(kind=kind_phys) :: flgmin(2)      = (/0.180,0.220/)          !< [in] ice fraction bounds
    real(kind=kind_phys) :: cgwf(2)        = (/0.5d0,0.05d0/)         !< multiplication factor for convective GWD
    real(kind=kind_phys) :: ccwf(2)        = (/1.0d0,1.0d0/)          !< multiplication factor for critical cloud
                                                                      !< workfunction for RAS
    real(kind=kind_phys) :: cdmbgwd(4)     = (/2.0d0,0.25d0,1.0d0,1.0d0/)   !< multiplication factors for cdmb, gwd, and NS gwd, tke based enhancement
    real(kind=kind_phys) :: alpha_fd       = 12.0                     !< alpha coefficient for turbulent orographic form drag
    real(kind=kind_phys) :: psl_gwd_dx_factor  = 6.0                  !< multiplication factors for grid spacing
    real(kind=kind_phys) :: sup            = 1.0                      !< supersaturation in pdf cloud (IMP_physics=98) when t is very low
                                                                      !< or ice super saturation in SHOC (when do_shoc=.true.)
    real(kind=kind_phys) :: ctei_rm(2)     = (/10.0d0,10.0d0/)        !< critical cloud top entrainment instability criteria
                                                                      !< (used if mstrat=.true.)
    real(kind=kind_phys) :: crtrh(3)       = (/0.90d0,0.90d0,0.90d0/) !< critical relative humidity at the surface
                                                                      !< PBL top and at the top of the atmosphere
    real(kind=kind_phys) :: dlqf(2)        = (/0.15,0.15/)            !< factor for cloud condensate detrainment
                                                                      !< from cloud edges for RAS
    real(kind=kind_phys) :: psauras(2)     = (/1.0d-3,1.0d-3/)        !< [in] auto conversion coeff from ice to snow in ras
    real(kind=kind_phys) :: prauras(2)     = (/2.0d-3,2.0d-3/)        !< [in] auto conversion coeff from cloud to rain in ras
    real(kind=kind_phys) :: wminras(2)     = (/1.0d-6,1.0d-6/)        !< [in] water and ice minimum threshold for ras
    integer              :: nrcmax         = 32                       !< number of random numbers used in RAS

    real(kind=kind_phys) :: rbcr           = 0.25                     !< Critical Richardson Number in PBL scheme
    real(kind=kind_phys) :: shoc_parm(5)   = (/7000.0,1.0,4.2857143,0.7,-999.0/)  !< some tunable parameters for shoc

!--- Rayleigh friction
    real(kind=kind_phys) :: prslrd0        = 0.0d0           !< pressure level from which Rayleigh Damping is applied
    real(kind=kind_phys) :: ral_ts         = 0.0d0           !< time scale for Rayleigh damping in days

!--- mass flux deep convection
    real(kind=kind_phys) :: clam_deep      = 0.1             !< c_e for deep convection (Han and Pan, 2011, eq(6))
    real(kind=kind_phys) :: c0s_deep       = 0.002           !< convective rain conversion parameter
    real(kind=kind_phys) :: c1_deep        = 0.002           !< conversion parameter of detrainment from liquid water into grid-scale cloud water
    real(kind=kind_phys) :: betal_deep     = 0.05            !< fraction factor of downdraft air mass reaching ground surface over land
    real(kind=kind_phys) :: betas_deep     = 0.05            !< fraction factor of downdraft air mass reaching ground surface over sea
    real(kind=kind_phys) :: evef           = 0.09            !< evaporation factor from convective rain
    real(kind=kind_phys) :: evfact_deep    = 0.3             !< evaporation factor from convective rain
    real(kind=kind_phys) :: evfactl_deep   = 0.3             !< evaporation factor from convective rain over land
    real(kind=kind_phys) :: pgcon_deep     = 0.55            !< reduction factor in momentum transport due to convection induced pressure gradient force
                                                             !< 0.7 : Gregory et al. (1997, QJRMS)
                                                             !< 0.55: Zhang & Wu (2003, JAS)
    real(kind=kind_phys) :: asolfac_deep   = 0.958           !< aerosol-aware parameter based on Lim (2011)
                                                             !< asolfac= cx / c0s(=.002)
                                                             !< cx = min([-0.7 ln(Nccn) + 24]*1.e-4, c0s)
                                                             !< Nccn: CCN number concentration in cm^(-3)
                                                             !< Until a realistic Nccn is provided, Nccns are assumed
                                                             !< as Nccn=100 for sea and Nccn=1000 for land

!--- mass flux shallow convection
    real(kind=kind_phys) :: clam_shal      = 0.3             !< c_e for shallow convection (Han and Pan, 2011, eq(6))
    real(kind=kind_phys) :: c0s_shal       = 0.002           !< conversion parameter of detrainment from liquid water into convetive precipitaiton
    real(kind=kind_phys) :: c1_shal        = 5.e-4           !< conversion parameter of detrainment from liquid water into grid-scale cloud water
    real(kind=kind_phys) :: pgcon_shal     = 0.55            !< reduction factor in momentum transport due to convection induced pressure gradient force
                                                             !< 0.7 : Gregory et al. (1997, QJRMS)
                                                             !< 0.55: Zhang & Wu (2003, JAS)
    real(kind=kind_phys) :: asolfac_shal   = 0.958           !< aerosol-aware parameter based on Lim (2011)
                                                             !< asolfac= cx / c0s(=.002)
                                                             !< cx = min([-0.7 ln(Nccn) + 24]*1.e-4, c0s)
                                                             !< Nccn: CCN number concentration in cm^(-3)
                                                             !< Until a realistic Nccn is provided, Nccns are assumed
                                                             !< as Nccn=100 for sea and Nccn=1000 for land

!--- near surface sea temperature model
    logical              :: nst_anl        = .false.         !< flag for NSSTM analysis in gcycle/sfcsub
    integer              :: lsea           = 0
    integer              :: nstf_name(5)   = (/0,0,1,0,5/)   !< flag 0 for no nst  1 for uncoupled nst  and 2 for coupled NST
                                                             !< nstf_name(1) : 0 = NSSTM off, 1 = NSSTM on but uncoupled
                                                             !<                2 = NSSTM on and coupled
                                                             !< nstf_name(2) : 1 = NSSTM spin up on, 0 = NSSTM spin up off
                                                             !< nstf_name(3) : 1 = NSSTM analysis on, 0 = NSSTM analysis off
                                                             !< nstf_name(4) : zsea1 in mm
                                                             !< nstf_name(5) : zsea2 in mm
!--- fractional grid
    logical              :: frac_grid       = .false.         !< flag for fractional grid
    logical              :: frac_ice        = .true.          !< flag for lake fractional ice when fractional grid is not in use
    logical              :: ignore_lake     = .true.          !< flag for ignoring lakes
    real(kind=kind_phys) :: min_lakeice     = 0.15d0          !< minimum lake ice value
    real(kind=kind_phys) :: min_seaice      = 1.0d-11         !< minimum sea  ice value
    real(kind=kind_phys) :: min_lake_height = 250.0           !< minimum lake height value
    real(kind=kind_phys) :: rho_h2o         = rhowater        !< fresh water density

!--- surface layer z0 scheme
    integer              :: sfc_z0_type    = 0               !< surface roughness options over ocean
                                                             !< 0=no change
                                                             !< 6=areodynamical roughness over water with input 10-m wind
                                                             !< 7=slightly decrease Cd for higher wind speed compare to 6
                                                             !< negative when cplwav2atm=.true. - i.e. two way wave coupling
    integer              :: icplocn2atm    = 0               !< air_sea_flux options over ocean
                                                             !< 0=ocean current is not used in the computation of air_sea fluxes
                                                             !< 1=including ocean current in the computation of air_sea fluxes

!--- potential temperature definition in surface layer physics
    logical              :: thsfc_loc      = .true.          !< flag for local vs. standard potential temperature
!--- flux method in 2-m diagnostics
    logical              :: diag_flux      = .false.         !< flag for flux method in 2-m diagnostics
!--- flux method in 2-m diagnostics (for stable conditions)
    logical              :: diag_log       = .false.         !< flag for log method in 2-m diagnostics (for stable conditions)
                                                             !<.true. means use local (gridpoint) surface pressure to define potential temperature
                                                             !<       this is the current GFS physics approach
                                                             !<.false. means use reference pressure of 1000 hPa to define potential temperature
                                                             !<       this is the alternative method proposed by GSL

!--- vertical diffusion
    real(kind=kind_phys) :: xkzm_m         = 1.0d0           !< [in] bkgd_vdif_m  background vertical diffusion for momentum
    real(kind=kind_phys) :: xkzm_h         = 1.0d0           !< [in] bkgd_vdif_h  background vertical diffusion for heat q
    real(kind=kind_phys) :: xkzm_s         = 1.0d0           !< [in] bkgd_vdif_s  sigma threshold for background mom. diffusion
    real(kind=kind_phys) :: xkzminv        = 0.3             !< diffusivity in inversion layers
    real(kind=kind_phys) :: moninq_fac     = 1.0             !< turbulence diffusion coefficient factor
    real(kind=kind_phys) :: dspfac         = 1.0             !< tke dissipative heating factor
    real(kind=kind_phys) :: bl_upfr        = 0.13            !< updraft fraction in boundary layer mass flux scheme
    real(kind=kind_phys) :: bl_dnfr        = 0.1             !< downdraft fraction in boundary layer mass flux scheme
    real(kind=kind_phys) :: rlmx           = 300.            !< maximum allowed mixing length in boundary layer mass flux scheme
    real(kind=kind_phys) :: elmx           = 300.            !< maximum allowed dissipation mixing length in boundary layer mass flux scheme
    integer              :: sfc_rlm        = 0               !< choice of near surface mixing length in boundary layer mass flux scheme
    integer              :: tc_pbl         = 0               !< control for TC applications in the PBL scheme

!--- parameters for canopy heat storage (CHS) parameterization
    real(kind=kind_phys) :: h0facu         = 0.25
    real(kind=kind_phys) :: h0facs         = 1.0

!---Cellular automaton options
    integer              :: nca            = 1
    integer              :: ncells         = 5
    integer              :: nlives         = 12

    integer              :: nca_g          = 1
    integer              :: ncells_g       = 1
    integer              :: nlives_g       = 100
    real(kind=kind_phys) :: nfracseed      = 0.5
    integer              :: nseed          = 1
    integer              :: nseed_g        = 100
    integer              :: iseed_ca       = 1
    integer              :: nspinup        = 1
    logical              :: do_ca          = .false.
    logical              :: ca_advect      = .false.
    logical              :: ca_sgs         = .false.
    logical              :: ca_global      = .false.
    logical              :: ca_smooth      = .false.
    real(kind=kind_phys) :: nthresh        = 18
    real                 :: ca_amplitude   = 0.35
    integer              :: nsmooth        = 100
    logical              :: ca_closure     = .false.
    logical              :: ca_entr        = .false.
    logical              :: ca_trigger     = .false.

!--- Increment grid
    logical               :: increment_file_on_native_grid = .false. ! increment on native grid else Gaussian grid
    
!--- IAU options
    real(kind=kind_phys)  :: iau_delthrs      = 0           !< iau time interval (to scale increments)
    character(len=240)    :: iau_inc_files(7) = ''          !< list of increment files
    real(kind=kind_phys)  :: iaufhrs(7)       = -1          !< forecast hours associated with increment files
    logical  :: iau_filter_increments         = .false.     !< filter IAU increments
    logical  :: iau_drymassfixer              = .false.     !< IAU dry mass fixer

!--- debug flags
    logical              :: debug          = .false.
    logical              :: pre_rad        = .false.         !< flag for testing purpose
    logical              :: print_diff_pgr = .false.         !< print average change in pgr every timestep

!  max and min lon and lat for critical relative humidity
    integer :: max_lon=5000, max_lat=2000, min_lon=192, min_lat=94
    real(kind=kind_phys) :: rhcmax = 0.9999999               !< max critical rel. hum.
#ifdef SINGLE_PREC
    real(kind=kind_phys) :: huge   = 9.9692099683868690E30  !  NetCDF float FillValue
#else
    real(kind=kind_phys) :: huge   = 9.9692099683868690E36  !  NetCDF float FillValue
#endif


!--- stochastic physics control parameters
    logical :: do_sppt      = .false.
    logical :: pert_mp      = .false.
    logical :: pert_clds    = .false.
    logical :: pert_radtend = .true.
    logical :: use_zmtnblck = .false.
    logical :: do_shum      = .false.
    logical :: do_skeb      = .false.
    integer :: skeb_npass   = 11
    integer :: lndp_type      = 0
    integer :: n_var_lndp     = 0
    logical :: lndp_each_step = .false.
    integer :: n_var_spp    =  0
    integer :: spp_pbl      =  0
    integer :: spp_sfc      =  0
    integer :: spp_mp       =  0
    integer :: spp_rad      =  0
    integer :: spp_gwd      =  0
    integer :: spp_cu_deep  =  0
    logical :: do_spp       = .false.

    integer              :: ichoice         = 0 !< flag for closure of C3/GF deep convection
    integer              :: ichoicem        = 13!< flag for closure of C3/GF mid convection
    integer              :: ichoice_s       = 3 !< flag for closure of C3/GF shallow convection

!-- chem nml variables for RRFS-SD
    real(kind=kind_phys) :: dust_drylimit_factor  = 1.0
    real(kind=kind_phys) :: dust_moist_correction = 1.0
    real(kind=kind_phys) :: dust_alpha = 0.
    real(kind=kind_phys) :: dust_gamma = 0.
    real(kind=kind_phys) :: wetdep_ls_alpha = 0.5
    integer :: dust_moist_opt = 1         ! fecan :1  else shao
    integer :: ebb_dcycle = 1             ! 1:retro; 2:forecast
    integer :: seas_opt = 2
    integer :: dust_opt = 1
    integer :: drydep_opt  = 1
    integer :: coarsepm_settling  = 1
    integer :: plume_wind_eff = 1
    logical :: extended_sd_diags = .false.
    integer :: wetdep_ls_opt  = 1
    logical :: do_plumerise   = .false.
    integer :: addsmoke_flag  = 1
    integer :: plumerisefire_frq = 60
    integer :: n_dbg_lines = 3
    integer :: smoke_forecast = 0         ! RRFS-sd read in ebb_smoke
    logical :: aero_ind_fdb = .false.     ! RRFS-sd wfa/ifa emission
    logical :: aero_dir_fdb = .false.     ! RRFS-sd smoke/dust radiation feedback
    logical :: rrfs_smoke_debug = .false. ! RRFS-sd plumerise debug
    logical :: do_smoke_transport = .true.! RRFS-sd convective transport of smoke/dust
    logical :: mix_chem = .false.         ! tracer mixing option by MYNN PBL
    logical :: enh_mix  = .false.         ! enhance vertmix option by MYNN PBL
    real(kind=kind_phys) :: smoke_dir_fdb_coef(7) =(/ 0.33, 0.67, 0.02, 0.13, 0.85, 0.05, 0.95 /) !< smoke & dust direct feedbck coefficents
    real(kind=kind_phys) :: smoke_conv_wet_coef(3) =(/ 0.50, 0.50, 0.50 /) !< smoke & dust convective wet removal coefficents

!-- Lightning threat index
    logical :: lightning_threat = .false.

!--- NRL Ozone physics
    logical              :: oz_phys      = .false.  !< Flag for old (2006) ozone physics
    logical              :: oz_phys_2015 = .true.   !< Flag for new (2015) ozone physics
    integer              :: kozpl        = 28       !< File identifier for ozone forcing data
    integer              :: kozc         = 48       !< File identifier for ozone climotology data

!--- NRL h2o photochemistry physics
    integer              :: kh2oc        = 29       !< File identifier for h2o photochemistry data.
    
!--- aerosol scavenging factors
    integer, parameter :: max_scav_factors = 183
    character(len=40)  :: fscav_aero(max_scav_factors) = ''

    real(kind=kind_phys) :: radar_tten_limits(2) = (/ limit_unspecified, limit_unspecified /)
    integer :: itime

!--- END NAMELIST VARIABLES

    NAMELIST /gfs_physics_nml/                                                              &
                          !--- general parameters
                               fhzero, fhzero_array, fhzero_fhour, ldiag3d, qdiag3d, lssav, &
                               naux2d, dtend_select, naux3d, aux2d_time_avg,                &
                               aux3d_time_avg, fhcyc, thermodyn_id, sfcpress_id,            &
                          !--- coupling parameters
                               cplflx, cplice, cplocn2atm, cplwav, cplwav2atm, cplaqm,      &
                               cplchm, cpllnd, cpllnd2atm, cpl_imp_mrg, cpl_imp_dbg,        &
                               cpl_fire, rrfs_sd, use_cice_alb,                             & 
#ifdef IDEA_PHYS
                               lsidea, weimer_model, f107_kp_size, f107_kp_interval,        &
                               f107_kp_skip_size, f107_kp_data_size, f107_kp_read_in_start, &
                               ipe_to_wam_coupling,                                         &
#else
                               lsidea, use_med_flux,                                        &
#endif
                          !--- radiation parameters
                               fhswr, fhlwr, levr, nfxr, iaerclm, iflip, isol, ico2, ialb,  &
                               isot, iems, iaer, icliq_sw, iovr, ictm, isubc_sw,            &
                               isubc_lw, lcrick, lcnorm, lwhtr, swhtr,                      &
                               nhfrad, idcor, dcorr_con,                                    &
                          ! --- RRTMGP
                               do_RRTMGP, active_gases, nGases, rrtmgp_root,                &
                               lw_file_gas, lw_file_clouds, rrtmgp_nBandsLW, rrtmgp_nGptsLW,&
                               sw_file_gas, sw_file_clouds, rrtmgp_nBandsSW, rrtmgp_nGptsSW,&
                               doG_cldoptics, doGP_cldoptics_PADE, doGP_cldoptics_LUT,      &
                               rrtmgp_nrghice, rrtmgp_nGauss_ang, do_GPsw_Glw,              &
                               use_LW_jacobian, doGP_lwscat, damp_LW_fluxadj, lfnc_k,       &
                               lfnc_p0, iovr_convcld, doGP_sgs_cnv, doGP_sgs_mynn,          &
                               rrtmgp_lw_phys_blksz, rrtmgp_sw_phys_blksz,                  &
                          ! IN CCN forcing
                               iccn, mraerosol,                                             &
                          !--- microphysical parameterizations
                               imp_physics, psautco, prautco, evpco, wminco,                &
                               fprcp, pdfflag, mg_dcs, mg_qcvar, mg_ts_auto_ice, mg_rhmini, &
                               effr_in, tf, tcr,                                            &
                               microp_uniform, do_cldice, hetfrz_classnuc,                  &
                               mg_do_graupel, mg_do_hail, mg_nccons, mg_nicons, mg_ngcons,  &
                               mg_ncnst, mg_ninst, mg_ngnst, sed_supersat, do_sb_physics,   &
                               mg_alf,   mg_qcmin, mg_do_ice_gmao, mg_do_liq_liu,           &
                               ltaerosol, lradar, nsfullradar_diag, lrefres, ttendlim,      &
                               ext_diag_thompson, dt_inner, lgfdlmprad,                     &
                               sedi_semi, decfl,                                            &
                               nssl_cccn, nssl_alphah, nssl_alphahl,                        &
                               nssl_alphar, nssl_ehw0, nssl_ehlw0,                          &
                               nssl_invertccn, nssl_hail_on, nssl_ccn_on, nssl_3moment,     &
                          !--- max hourly
                               avg_max_length,                                              &
                          !--- land/surface model control
                               lsm, lsoil, lsoil_lsm, lsnow_lsm, kice, rdlai,               &
                               nmtvr, ivegsrc, use_ufo, iopt_thcnd, ua_phys, usemonalb,     &
                               aoasis, fasdas, exticeden, nvegcat, nsoilcat,                &
                          !    Noah MP options
                               iopt_dveg,iopt_crs,iopt_btr,iopt_run,iopt_sfc, iopt_frz,     &
                               iopt_inf, iopt_rad,iopt_alb,iopt_snf,iopt_tbot,iopt_stc,     &
                               iopt_trs, iopt_diag,                                         &
                          !    RUC lsm options
                               add_fire_heat_flux,                                          &
                               mosaic_lu, mosaic_soil, isncond_opt, isncovr_opt,            &
                          !    GFDL surface layer options
                               lcurr_sf, pert_cd, ntsflg, sfenth,                           &
                          !--- lake model control
                               lkm, iopt_lake, lakedepth_threshold, lakefrac_threshold,     &
                               clm_lake_depth_default, clm_lake_use_lakedepth,              &
                               clm_lake_debug, clm_debug_print, use_lake2m,                 &
                          !--- physical parameterizations
                               ras, trans_trac, old_monin, cnvgwd, mstrat, moist_adj,       &
                               cscnv, cal_pre, do_aw, do_shoc, shocaftcnv, shoc_cld,        &
                               oz_phys, oz_phys_2015,                                       &
                               do_mynnedmf, do_mynnsfclay,                                  &
                               ! DH* TODO - move to MYNN namelist section
                               bl_mynn_cloudpdf, bl_mynn_edmf, bl_mynn_edmf_mom,            &
                               bl_mynn_edmf_tke, bl_mynn_mixlength, bl_mynn_cloudmix,       &
                               bl_mynn_mixqt, bl_mynn_output, icloud_bl, bl_mynn_tkeadvect, &
                               bl_mynn_closure, tke_budget,                                 &
                               isftcflx, iz0tlnd, sfclay_compute_flux, sfclay_compute_diag, &
                               ! *DH
                               gwd_opt, do_ugwp_v0, do_ugwp_v0_orog_only,                   &
                               do_ugwp_v0_nst_only,                                         &
                               do_gsl_drag_ls_bl, do_gsl_drag_ss, do_gsl_drag_tofd,         &
                               do_gwd_opt_psl,                                              &
                               do_ugwp_v1, do_ugwp_v1_orog_only,  do_ugwp_v1_w_gsldrag,     &
                               ugwp_seq_update, var_ric, coef_ric_l, coef_ric_s, hurr_pbl,  &
                               do_myjsfc, do_myjpbl,                                        &
                               hwrf_samfdeep, hwrf_samfshal,progsigma,betascu,betamcu,      &
                               betadcu,h2o_phys, pdfcld, shcnvcw, redrag, hybedmf, satmedmf,&
                               shinhong, do_ysu, dspheat, lheatstrg, lseaspray, cnvcld,     &
                               xr_cnvcld, random_clds, shal_cnv, imfshalcnv, imfdeepcnv,    &
                               isatmedmf, do_deep, jcap,                                    &
                               cs_parm, flgmin, cgwf, ccwf, cdmbgwd, alpha_fd,              &
                               psl_gwd_dx_factor,                                           &
                               sup, ctei_rm, crtrh,                                         &
                               dlqf, rbcr, shoc_parm, psauras, prauras, wminras,            &
                               do_sppt, do_shum, do_skeb,                                   &
                               do_spp, n_var_spp,                                           &
                               lndp_type,  n_var_lndp, lndp_each_step,                      &
                               pert_mp,pert_clds,pert_radtend,                              &
                          !--- Rayleigh friction
                               prslrd0, ral_ts,  ldiag_ugwp, do_ugwp, do_tofd,              &
                          ! --- Ferrier-Aligo
                               spec_adv, rhgrd, icloud,                                     &
                          !--- mass flux deep convection
                               clam_deep, c0s_deep, c1_deep, betal_deep,                    &
                               betas_deep, evef, evfact_deep, evfactl_deep, pgcon_deep,     &
                               asolfac_deep,                                                &
                          !--- mass flux shallow convection
                               clam_shal, c0s_shal, c1_shal, pgcon_shal, asolfac_shal,      &
                          !--- near surface sea temperature model
                               nst_anl, lsea, nstf_name,                                    &
                               frac_grid, min_lakeice, min_seaice, min_lake_height,         &
                               ignore_lake, frac_ice,                                       &
                          !--- surface layer
                               sfc_z0_type, icplocn2atm,                                    &
                          !--- switch beteeen local and standard potential temperature
                               thsfc_loc,                                                   &
                          !--- switches in 2-m diagnostics
                               diag_flux, diag_log,                                         &
                          !    vertical diffusion
                               xkzm_m, xkzm_h, xkzm_s, xkzminv, moninq_fac, dspfac,         &
                               bl_upfr, bl_dnfr, rlmx, elmx, sfc_rlm, tc_pbl,               &
                          !--- canopy heat storage parameterization
                               h0facu, h0facs,                                              &
                          !--- cellular automata
                               nca, ncells, nlives, nca_g, ncells_g, nlives_g, nfracseed,   &
                               nseed,  nseed_g,  nthresh, do_ca, ca_advect,                 &
                               ca_sgs, ca_global,iseed_ca,ca_smooth,                        &
                               nspinup,ca_amplitude,nsmooth,ca_closure,ca_entr,ca_trigger,  &
                          !--- Increment grid
                               increment_file_on_native_grid,                               &
                          !--- IAU
                               iau_delthrs,iaufhrs,iau_inc_files,iau_filter_increments,     &
                               iau_drymassfixer,                                            &
                          !--- debug options
                               debug, pre_rad, print_diff_pgr,                              &
                          !--- parameter range for critical relative humidity
                               max_lon, max_lat, min_lon, min_lat, rhcmax, huge,            &
                               phys_version,                                                &
                          !--- aerosol scavenging factors ('name:value' string array)
                               fscav_aero,                                                  &
                          !--- RRFS-SD namelist
                               dust_drylimit_factor, dust_moist_correction, dust_moist_opt, &
                               dust_alpha, dust_gamma, wetdep_ls_alpha,                     &
                               seas_opt, dust_opt, drydep_opt, coarsepm_settling,           &
                               plume_wind_eff,ebb_dcycle, extended_sd_diags,                &
                               wetdep_ls_opt, smoke_forecast, aero_ind_fdb, aero_dir_fdb,   &
                               rrfs_smoke_debug, do_plumerise, plumerisefire_frq,           &
                               addsmoke_flag, enh_mix, mix_chem, smoke_dir_fdb_coef,        &
                               do_smoke_transport,smoke_conv_wet_coef,n_dbg_lines,          &
                          !--- C3/GF closures
                               ichoice,ichoicem,ichoice_s,                                  &
                          !--- (DFI) time ranges with radar-prescribed microphysics tendencies
                          !          and (maybe) convection suppression
                               fh_dfi_radar, radar_tten_limits, do_cap_suppress,            &
                          !--- GSL lightning threat indices
                               lightning_threat

!--- other parameters
    integer :: nctp    =  0                !< number of cloud types in CS scheme
    logical :: gen_coord_hybrid = .false.  !< for Henry's gen coord

!--- SHOC parameters
    integer :: nshoc_2d  = 0  !< number of 2d fields for SHOC
    integer :: nshoc_3d  = 0  !< number of 3d fields for SHOC

!--- convective clouds
    integer :: ncnvcld3d = 0       !< number of convective 3d clouds fields

    integer :: itrac, ipat, ichem
    logical :: have_pbl, have_dcnv, have_scnv, have_mp, have_oz_phys, have_samf, have_pbl_edmf, have_cnvtrans, have_rdamp
    character(len=20) :: namestr
    character(len=44) :: descstr

!--- NRL ozone physics
    character(len=128) :: err_message

    ! dtend selection: default is to match all variables:
    dtend_select(1)='*'
    do ipat=2,pat_count
       dtend_select(ipat)=' '
    enddo

!--- read in the namelist
#ifdef INTERNAL_FILE_NML
    ! allocate required to work around GNU compiler bug 100886 https://gcc.gnu.org/bugzilla/show_bug.cgi?id=100886
    allocate(Model%input_nml_file, mold=input_nml_file)
    Model%input_nml_file => input_nml_file
    read(Model%input_nml_file, nml=gfs_physics_nml)
    ! Set length (number of lines) in namelist for internal reads
    Model%input_nml_file_length = size(Model%input_nml_file)
#else
    inquire (file=trim(fn_nml), exist=exists)
    if (.not. exists) then
      write(6,*) 'GFS_namelist_read:: namelist file: ',trim(fn_nml),' does not exist'
      stop
    else
      open (unit=nlunit, file=fn_nml, action='READ', status='OLD', iostat=ios)
    endif
    rewind(nlunit)
    read (nlunit, nml=gfs_physics_nml)
    close (nlunit)
    ! Set length (number of lines) in namelist for internal reads
    Model%input_nml_file_length = 0
#endif
!--- write version number and namelist to log file ---
    if (me == master) then
      write(logunit, '(a80)') '================================================================================'
      write(logunit, '(a64)') phys_version
      write(logunit, nml=gfs_physics_nml)
    endif

!--- MPI parameters
    Model%me               = me
    Model%master           = master
    Model%communicator     = communicator
    Model%ntasks           = ntasks
    Model%nthreads         = nthreads
    Model%nlunit           = nlunit
    Model%fn_nml           = fn_nml
    Model%logunit          = logunit
    Model%fhzero           = fhzero
    Model%fhzero_array     = fhzero_array
    Model%fhzero_fhour     = fhzero_fhour
    if( Model%fhzero_array(1) > 0. ) then
      Model%fhzero = Model%fhzero_array(1)
    endif
    Model%ldiag3d          = ldiag3d
    Model%qdiag3d          = qdiag3d
    if (qdiag3d .and. .not. ldiag3d) then
      write(0,*) 'Logic error in GFS_typedefs.F90: qdiag3d requires ldiag3d'
      stop
    endif
    Model%flag_for_gwd_generic_tend = .true.
    Model%flag_for_pbl_generic_tend = .true.
    Model%flag_for_scnv_generic_tend = .true.
    Model%flag_for_dcnv_generic_tend = .true.

    Model%lightning_threat = lightning_threat

    Model%fh_dfi_radar     = fh_dfi_radar
    Model%num_dfi_radar    = 0
    Model%dfi_radar_max_intervals = dfi_radar_max_intervals ! module-level parameter, top of file
    Model%dfi_radar_max_intervals_plus_one = dfi_radar_max_intervals + 1
    Model%do_cap_suppress = do_cap_suppress

    call control_initialize_radar_tten(Model, radar_tten_limits)

    if(gwd_opt==1) then
      if(me==master) &
           write(*,*) 'FLAG: gwd_opt==1 so gwd not generic'
      Model%flag_for_gwd_generic_tend=.false.
    elseif(me==master) then
      write(*,*) 'NO FLAG: gwd is generic'
    endif

    if(satmedmf .and. isatmedmf==0) then
      if(me==master) &
           write(*,*) 'FLAG: satmedmf and isatedmf=0 so pbl not generic'
      Model%flag_for_pbl_generic_tend=.false.
    elseif(satmedmf .and. isatmedmf==1) then
      if(me==master) &
           write(*,*) 'FLAG: satmedmf and isatedmf=1 so pbl not generic'
      Model%flag_for_pbl_generic_tend=.false.
    else if(hybedmf) then
      if(me==master) &
           write(*,*) 'FLAG: hybedmf so pbl not generic'
      Model%flag_for_pbl_generic_tend=.false.
    else if(do_mynnedmf) then
      if(me==master) &
           write(*,*) 'FLAG: do_mynnedmf so pbl not generic'
      Model%flag_for_pbl_generic_tend=.false.
    elseif(me==master) then
      write(*,*) 'NO FLAG: pbl is generic'
    endif

    if(imfshalcnv == Model%imfshalcnv_gf .or. imfshalcnv == Model%imfshalcnv_c3) then
      if(me==master) &
           write(*,*) 'FLAG: imfshalcnv_gf or imfshalcnv_c3 so scnv not generic'
      Model%flag_for_scnv_generic_tend=.false.
    elseif(me==master) then
      write(*,*) 'NO FLAG: scnv is generic'
    endif

    if(imfdeepcnv == Model%imfdeepcnv_gf .or. imfdeepcnv == Model%imfdeepcnv_c3) then
      if(me==master) &
           write(*,*) 'FLAG: imfdeepcnv_gf or imfdeepcnv_c3 so dcnv not generic'
      Model%flag_for_dcnv_generic_tend=.false.
    elseif(me==master) then
      write(*,*) 'NO FLAG: dcnv is generic'
    endif

!
!VAY-ugwp  --- set some GW-related switches
!
    Model%ldiag_ugwp       = ldiag_ugwp
    Model%ugwp_seq_update  = ugwp_seq_update
    Model%do_ugwp          = do_ugwp
    Model%do_tofd          = do_tofd

    Model%lssav            = lssav
    !
    if (naux2d>naux2dmax) then
      write(0,*) "Error, number of requested auxiliary 2d arrays exceeds the maximum defined in GFS_typedefs.F90"
      stop
    endif
    if (naux3d>naux3dmax) then
      write(0,*) "Error, number of requested auxiliary 3d arrays exceeds the maximum defined in GFS_typedefs.F90"
      stop
    endif
    Model%naux2d           = naux2d
    Model%naux3d           = naux3d
    if (Model%naux2d>0) then
        allocate(Model%aux2d_time_avg(1:naux2d))
        Model%aux2d_time_avg(1:naux2d) = aux2d_time_avg(1:naux2d)
    end if
    if (Model%naux3d>0) then
        allocate(Model%aux3d_time_avg(1:naux3d))
        Model%aux3d_time_avg(1:naux3d) = aux3d_time_avg(1:naux3d)
    end if
    !
    Model%fhcyc            = fhcyc
    Model%thermodyn_id     = thermodyn_id
    Model%sfcpress_id      = sfcpress_id
    Model%gen_coord_hybrid = gen_coord_hybrid

    !--- set some grid extent parameters
    Model%tile_num         = tile_num
    Model%isc              = isc
    Model%jsc              = jsc
    Model%nx               = nx
    Model%ny               = ny
    Model%levs             = levs
    allocate(Model%ak(1:size(ak)))
    allocate(Model%bk(1:size(bk)))
    Model%ak               = ak
    Model%bk               = bk
    Model%levsp1           = Model%levs + 1
    Model%levsm1           = Model%levs - 1
    Model%cnx              = cnx
    Model%cny              = cny
    Model%lonr             = gnx         ! number longitudinal points
    Model%latr             = gny         ! number of latitudinal points from pole to pole
    Model%nblks            = size(blksz)
    allocate(Model%blksz(1:Model%nblks))
    Model%blksz            = blksz
    Model%ncols            = sum(Model%blksz)
    ! DH*
    Model%nchunks          = size(blksz)
    allocate(Model%chunk_begin(Model%nchunks))
    allocate(Model%chunk_end(Model%nchunks))
    Model%chunk_begin(1) = 1
    Model%chunk_end(1) = Model%chunk_begin(1) + blksz(1) - 1
    do i=2,Model%nchunks
        Model%chunk_begin(i) = Model%chunk_end(i-1) + 1
        Model%chunk_end(i) = Model%chunk_begin(i) + blksz(i) - 1
    end do

!--- coupling parameters
    Model%cplflx           = cplflx
    Model%cplice           = cplice
    ! Consistency check, currently allowed combinations are
    ! Model%cplflx == .false. and Model%cplice == .false. (uncoupled runs)
    ! Model%cplflx == .true.  and Model%cplice == .true.  (coupled S2S runs)
    ! Model%cplflx == .true.  and Model%cplice == .false. (HAFS FV3ATM-HYCOM)
    if (Model%cplice .and. .not. Model%cplflx) then
      print *,' Logic error: Model%cplflx==.false. and Model%cplice==.true. is currently not supported - shutting down'
      stop
    endif
    Model%cplocn2atm       = cplocn2atm
    Model%cplwav           = cplwav
    Model%cplwav2atm       = cplwav2atm
    Model%cplaqm           = cplaqm
    Model%cplchm           = cplchm .or. cplaqm
    Model%cpllnd           = cpllnd
    Model%cpllnd2atm       = cpllnd2atm
    Model%use_cice_alb     = use_cice_alb
    Model%cpl_imp_mrg      = cpl_imp_mrg
    Model%cpl_imp_dbg      = cpl_imp_dbg
    Model%use_med_flux     = use_med_flux

!--- RRFS-SD
    Model%rrfs_sd           = rrfs_sd
    Model%cpl_fire          = cpl_fire
    Model%dust_drylimit_factor = dust_drylimit_factor
    Model%dust_moist_correction = dust_moist_correction
    Model%dust_moist_opt    = dust_moist_opt
    Model%dust_alpha        = dust_alpha
    Model%dust_gamma        = dust_gamma
    Model%wetdep_ls_alpha   = wetdep_ls_alpha
    Model%ebb_dcycle        = ebb_dcycle
    Model%seas_opt          = seas_opt
    Model%dust_opt          = dust_opt
    Model%drydep_opt        = drydep_opt
    Model%coarsepm_settling = coarsepm_settling
    Model%plume_wind_eff    = plume_wind_eff
    Model%extended_sd_diags = extended_sd_diags
    Model%wetdep_ls_opt     = wetdep_ls_opt
    Model%do_plumerise      = do_plumerise
    Model%n_dbg_lines       = n_dbg_lines
    Model%plumerisefire_frq = plumerisefire_frq
    Model%addsmoke_flag     = addsmoke_flag
    Model%smoke_forecast    = smoke_forecast
    Model%aero_ind_fdb      = aero_ind_fdb
    Model%aero_dir_fdb      = aero_dir_fdb
    Model%rrfs_smoke_debug  = rrfs_smoke_debug
    Model%do_smoke_transport= do_smoke_transport
    Model%mix_chem          = mix_chem
    Model%enh_mix           = enh_mix
    Model%smoke_dir_fdb_coef  = smoke_dir_fdb_coef
    Model%smoke_conv_wet_coef = smoke_conv_wet_coef

    Model%fire_aux_data_levels = 1

    Model%ichoice_s = ichoice_s
    Model%ichoicem  = ichoicem
    Model%ichoice   = ichoice

!--- integrated dynamics through earth's atmosphere
    Model%lsidea           = lsidea
    if (Model%lsidea) then
      print *,' LSIDEA is active but needs to be reworked for FV3 - shutting down'
      stop
    endif
#ifdef IDEA_PHYS
!--- integrated dynamics through earth's atmosphere
    Model%weimer_model     = weimer_model
#endif

!--- calendars and time parameters and activation triggers
    Model%dtp              = dt_phys
    Model%dtf              = dt_dycore
    Model%nscyc            = nint(Model%fhcyc*con_hr/Model%dtp)
    Model%nszero           = nint(Model%fhzero*con_hr/Model%dtp)
    Model%idat(1:8)        = idat(1:8)
    Model%idate            = 0
    Model%idate(1)         = Model%idat(5)
    Model%idate(2)         = Model%idat(2)
    Model%idate(3)         = Model%idat(3)
    Model%idate(4)         = Model%idat(1)
    Model%iau_offset       = iau_offset

!--- radiation control parameters
    Model%fhswr            = fhswr
    Model%fhlwr            = fhlwr
    Model%nsswr            = nint(fhswr/Model%dtp)
    Model%nslwr            = nint(fhlwr/Model%dtp)
    if (restart) then
      Model%nhfrad         = 0
      if (Model%me == Model%master .and. nhfrad>0) &
        write(*,'(a)') 'Disable high-frequency radiation calls for restart run'
    else
      Model%nhfrad         = nhfrad
      if (Model%me == Model%master .and. nhfrad>0) &
        write(*,'(a,i0)') 'Number of high-frequency radiation calls for coldstart run: ', nhfrad
    endif

    if (levr < 0) then
      Model%levr           = levs
    else if (levr > levs) then
      write(0,*) "Logic error, number of radiation levels (levr) cannot exceed number of model levels (levs)"
      stop
    else
      Model%levr           = levr
    endif
    Model%levrp1           = Model%levr + 1

    if (isubc_sw < 0 .or. isubc_sw > 2) then
       write(0,'(a,i0)') 'ERROR: shortwave cloud-sampling (isubc_sw) scheme selected not valid: ',isubc_sw
       stop
    endif
    if (isubc_lw < 0 .or. isubc_lw > 2) then
       write(0,'(a,i0)') 'ERROR: longwave cloud-sampling (isubc_lw) scheme selected not valid: ',isubc_lw
       stop
    endif


    if ((iovr .ne. Model%iovr_rand) .and. (iovr .ne. Model%iovr_maxrand) .and.       &
        (iovr .ne. Model%iovr_max)  .and. (iovr .ne. Model%iovr_dcorr)   .and.       &
        (iovr .ne. Model%iovr_exp)  .and. (iovr .ne. Model%iovr_exprand)) then
       write(0,'(a,i0)') 'ERROR: cloud-overlap (iovr) scheme selected not valid: ',iovr
       stop
    endif

    if ((isubc_sw == 0 .or. isubc_lw == 0) .and. iovr > 2 ) then
        if (me == 0) then
           print *,'  *** IOVR=',iovr,' is not available for ISUBC_SW(LW)=0 setting!!'
           print *,'      The program will use maximum/random overlap instead.'
        endif
        iovr = 1
     endif

    Model%nfxr             = nfxr
    Model%iccn             = iccn
    ! further down: set Model%iccn to .false.
    ! for all microphysics schemes except
    ! MG2/3 (these are the only ones using ICCN)
    Model%iflip            = iflip
    Model%isol             = isol
    Model%ico2             = ico2
    Model%ialb             = ialb
    Model%iems             = iems
    Model%iaer             = iaer
    Model%iaerclm          = iaerclm
    if (iaer/1000 == 1 .or. Model%iccn == 2) then
      Model%iaerclm = .true.
      ntrcaer = ntrcaerm
    else if (iaer/1000 == 2) then
      ntrcaer = ntrcaerm
    else
      ntrcaer = 1
    endif
    Model%lalw1bd          = lalw1bd
    Model%iaerflg          = iaerflg
    Model%iaermdl          = iaermdl
    Model%aeros_file       = aeros_file
    Model%solar_file       = solar_file
    Model%semis_file       = semis_file
    Model%co2dat_file      = co2dat_file
    Model%co2gbl_file      = co2gbl_file
    Model%co2usr_file      = co2usr_file
    Model%co2cyc_file      = co2cyc_file
    Model%ntrcaer          = ntrcaer
    Model%idcor            = idcor
    Model%dcorr_con        = dcorr_con
    Model%icliq_sw         = icliq_sw
    Model%icice_sw         = icice_sw
    Model%icliq_lw         = icliq_lw
    Model%icice_lw         = icice_lw
    Model%iovr             = iovr
    Model%ictm             = ictm
    Model%isubc_sw         = isubc_sw
    Model%isubc_lw         = isubc_lw
    Model%iswmode          = iswmode
    Model%lcrick           = lcrick
    Model%lcnorm           = lcnorm
    Model%lwhtr            = lwhtr
    Model%swhtr            = swhtr
    Model%rad_hr_units     = rad_hr_units
    Model%inc_minor_gas    = inc_minor_gas
    Model%ipsd0            = ipsd0
    Model%ipsdlim          = ipsdlim
    Model%lrseeds          = lrseeds
    Model%nrstreams        = nrstreams
    Model%lextop           = (ltp > 0)

    ! RRTMGP
    Model%do_RRTMGP           = do_RRTMGP
    Model%rrtmgp_nrghice      = rrtmgp_nrghice
    Model%rrtmgp_nGauss_ang   = rrtmgp_nGauss_ang
    Model%do_GPsw_Glw         = do_GPsw_Glw
    Model%active_gases        = active_gases
    Model%ngases              = nGases
    if (Model%do_RRTMGP) then
      allocate (Model%active_gases_array(Model%nGases))
      ! Reset, will be populated by RRTMGP
      do ipat=1,Model%nGases
        Model%active_gases_array(ipat) = ''
      enddo
    endif
    Model%rrtmgp_root         = rrtmgp_root
    Model%lw_file_gas         = lw_file_gas
    Model%lw_file_clouds      = lw_file_clouds
    Model%rrtmgp_nBandsLW     = rrtmgp_nBandsLW
    Model%rrtmgp_nGptsLW      = rrtmgp_nGptsLW
    Model%sw_file_gas         = sw_file_gas
    Model%sw_file_clouds      = sw_file_clouds
    Model%rrtmgp_nBandsSW     = rrtmgp_nBandsSW
    Model%rrtmgp_nGptsSW      = rrtmgp_nGptsSW
    Model%doG_cldoptics       = doG_cldoptics
    Model%doGP_cldoptics_PADE = doGP_cldoptics_PADE
    Model%doGP_cldoptics_LUT  = doGP_cldoptics_LUT
    Model%iovr_convcld        = iovr_convcld
    Model%use_LW_jacobian     = use_LW_jacobian
    Model%damp_LW_fluxadj     = damp_LW_fluxadj
    Model%lfnc_k              = lfnc_k
    Model%lfnc_p0             = lfnc_p0
    Model%doGP_lwscat         = doGP_lwscat
    Model%doGP_sgs_cnv        = doGP_sgs_cnv
    Model%doGP_sgs_mynn       = doGP_sgs_mynn
    Model%rrtmgp_lw_phys_blksz   = rrtmgp_lw_phys_blksz
    Model%rrtmgp_sw_phys_blksz   = rrtmgp_sw_phys_blksz
    if (Model%do_RRTMGP) then
       ! RRTMGP incompatible with levr /= levs
       if (Model%levr /= Model%levs) then
          write(0,*) "Logic error, RRTMGP only works with levr = levs"
          stop
       end if
       ! RRTMGP LW scattering calculation not supported w/ RRTMG cloud-optics
       if (Model%doGP_lwscat .and. Model%doG_cldoptics) then
          write(0,*) "Logic error, RRTMGP Longwave cloud-scattering not supported with RRTMG cloud-optics."
          stop
       end if
       if (Model%doGP_sgs_mynn .and. .not. do_mynnedmf) then
          write(0,*) "Logic error, RRTMGP flag doGP_sgs_mynn only works with do_mynnedmf=.true."
          stop
       endif
       if (Model%doGP_sgs_cnv .or. Model%doGP_sgs_mynn) then
          write(0,*) "RRTMGP explicit cloud scheme being used."
          Model%doGP_smearclds = .false.
       else
           write(0,*) "RRTMGP implicit cloud scheme being used."
       endif

       if (Model%doGP_cldoptics_PADE .and. Model%doGP_cldoptics_LUT) then
          write(0,*) "Logic error, Both RRTMGP cloud-optics options cannot be selected. "
          stop
       end if
       if (.not. Model%doGP_cldoptics_PADE .and. .not. Model%doGP_cldoptics_LUT .and. .not. Model%doG_cldoptics) then
          write(0,*) "Logic error, No option for cloud-optics scheme provided. Using RRTMG cloud-optics"
          Model%doG_cldoptics = .true.
       end if
       if (Model%rrtmgp_nGptsSW  .lt. 0 .or. Model%rrtmgp_nGptsLW  .lt. 0 .or. &
           Model%rrtmgp_nBandsSW .lt. 0 .or. Model%rrtmgp_nBandsLW .lt. 0) then
          write(0,*) "Logic error, RRTMGP spectral dimensions (bands/gpts) need to be provided."
          stop
       endif
       else
          if (Model%use_LW_jacobian) then
             write(0,*) "Logic error, RRTMGP LW Jacobian adjustment cannot be used with RRTMG radiation."
             Model%use_LW_jacobian = .false.
             Model%damp_LW_fluxadj = .false.
          endif
    endif

    ! The CCPP versions of the RRTMG lw/sw schemes are configured
    ! such that lw and sw heating rate are output, i.e. they rely
    ! on the corresponding arrays to be allocated.
    if (.not.lwhtr .or. .not.swhtr) then
      write(0,*) "Logic error, the CCPP version of RRTMG lwrad/swrad require the output" // &
             " of the lw/sw heating rates to be turned on (namelist options lwhtr and swhtr)"
      stop
    end if

!--- microphysical switch
    Model%imp_physics      = imp_physics
!--- use effective radii in radiation, used by several microphysics options
    Model%effr_in          = effr_in
    ! turn off ICCN interpolation when MG2/3 are not used
    if (.not. Model%imp_physics==Model%imp_physics_mg) Model%iccn = 0
!--- Zhao-Carr MP parameters
    Model%psautco          = psautco
    Model%prautco          = prautco
    Model%evpco            = evpco
    Model%wminco           = wminco
!--- Max hourly
    Model%avg_max_length   = avg_max_length
!--- Morrison-Gettelman MP parameters
    Model%fprcp            = fprcp
    Model%pdfflag          = pdfflag
    Model%mg_dcs           = mg_dcs
    Model%mg_qcvar         = mg_qcvar
    Model%mg_ts_auto_ice   = mg_ts_auto_ice
    Model%mg_rhmini        = mg_rhmini
    Model%mg_alf           = mg_alf
    Model%mg_qcmin         = mg_qcmin
    Model%microp_uniform   = microp_uniform
    Model%do_cldice        = do_cldice
    Model%hetfrz_classnuc  = hetfrz_classnuc
    Model%mg_do_graupel    = mg_do_graupel
    Model%mg_do_hail       = mg_do_hail
    Model%mg_do_ice_gmao   = mg_do_ice_gmao
    Model%mg_do_liq_liu    = mg_do_liq_liu
    Model%mg_nccons        = mg_nccons
    Model%mg_nicons        = mg_nicons
    Model%mg_ngcons        = mg_ngcons
    Model%mg_ncnst         = mg_ncnst
    Model%mg_ninst         = mg_ninst
    Model%mg_ngnst         = mg_ngnst
    Model%sed_supersat     = sed_supersat
    Model%do_sb_physics    = do_sb_physics
    Model%mg_precip_frac_method  = mg_precip_frac_method
    Model%mg_berg_eff_factor     = mg_berg_eff_factor
    Model%tf               = tf
    Model%tcr              = tcr
    Model%tcrf             = 1.0/(tcr-tf)

!-- NSSL microphysics params
    Model%nssl_cccn        = nssl_cccn
    Model%nssl_alphah      = nssl_alphah
    Model%nssl_alphahl     = nssl_alphahl
    Model%nssl_alphar      = nssl_alphar
    Model%nssl_ehw0        = nssl_ehw0
    Model%nssl_ehlw0       = nssl_ehlw0
    Model%nssl_hail_on     = nssl_hail_on
    Model%nssl_ccn_on      = nssl_ccn_on
    Model%nssl_invertccn   = nssl_invertccn
    Model%nssl_3moment     = nssl_3moment

!--- Thompson MP parameters
    Model%ltaerosol        = ltaerosol
    Model%mraerosol        = mraerosol
    if (Model%ltaerosol .and. Model%mraerosol) then
      write(0,*) 'Logic error: Only one Thompson aerosol option can be true, either ltaerosol or mraerosol)'
      stop
    end if
    Model%lradar           = lradar
    Model%nsfullradar_diag = nsfullradar_diag
    Model%ttendlim         = ttendlim
    Model%ext_diag_thompson= ext_diag_thompson
    if (dt_inner>0) then
      Model%dt_inner       = dt_inner
    else
      Model%dt_inner       = Model%dtp
    endif
    Model%sedi_semi        = sedi_semi
    Model%decfl            = decfl
!--- F-A MP parameters
    Model%rhgrd            = rhgrd
    Model%spec_adv         = spec_adv
    Model%icloud           = icloud

!--- GFDL MP parameters
    Model%lgfdlmprad       = lgfdlmprad
!--- Thompson,GFDL,NSSL MP parameter
    Model%lrefres          = lrefres

!--- land/surface model parameters
    Model%lsm              = lsm
    Model%lsoil            = lsoil

    ! Flag to read leaf area index from input files (initial conditions)
    Model%rdlai = rdlai
    if (Model%rdlai .and. .not. Model%lsm == Model%lsm_ruc) then
      write(0,*) 'Logic error: rdlai = .true. only works with RUC LSM'
      stop
    end if

    ! Set surface layers for CCPP physics
    if (lsoil_lsm==-1) then
      Model%lsoil_lsm      = lsoil
    else
      Model%lsoil_lsm      = lsoil_lsm
    end if
    ! DH* TODO - need to clean up how different land surface models handle initializing zs and dzs
    ! For Noah and NoahMP, hardcode here for the moment; for RUC, these variables get initialized
    ! in the RUC LSM init calls.
    ! Allocate variables to store depth/thickness of soil layers
    allocate (Model%zs (Model%lsoil_lsm))
    allocate (Model%dzs(Model%lsoil_lsm))
    if (Model%lsm==Model%lsm_noah .or. Model%lsm==Model%lsm_noahmp) then
      if (Model%lsoil_lsm/=4) then
        write(0,*) 'Error in GFS_typedefs.F90, number of soil layers must be 4 for Noah/NoahMP'
        stop
      end if
      Model%zs  = (/-0.1_kind_phys, -0.4_kind_phys, -1.0_kind_phys, -2.0_kind_phys/)
      Model%dzs = (/ 0.1_kind_phys,  0.3_kind_phys,  0.6_kind_phys,  1.0_kind_phys/)
    elseif (Model%lsm==Model%lsm_ruc) then
      Model%zs  = clear_val
      Model%dzs = clear_val
    end if
    ! *DH

    if (Model%lsm==Model%lsm_ruc) then
      if (Model%lsoil_lsm/=9) then
        write(0,*) 'Error in GFS_typedefs.F90, number of soil layers must be 9 for RUC'
        stop
      end if
    end if

    ! Set number of ice model layers
    Model%kice      = kice

    if (Model%lsm==Model%lsm_noah .or. Model%lsm==Model%lsm_noahmp) then
      if (kice/=2) then
        write(0,*) 'Error in GFS_typedefs.F90, number of ice model layers must be 2 for Noah/NoahMP/Noah_WRFv4'
        stop
      end if
    elseif (Model%lsm==Model%lsm_ruc) then
      if (kice/=9) then
        write(0,*) 'Error in GFS_typedefs.F90, number of ice model layers must be 9 for RUC'
        stop
      end if
    end if

    ! Allocate variable for min/max soil moisture for a given soil type
    allocate (Model%pores(30))
    allocate (Model%resid(30))
    Model%pores    = clear_val
    Model%resid    = clear_val
    !
    if (Model%lsm==Model%lsm_noahmp) then
      if (lsnow_lsm/=3) then
        write(0,*) 'Logic error: NoahMP expects the maximum number of snow layers to be exactly 3 (see sfc_noahmp_drv.f)'
        stop
      else
        Model%lsnow_lsm        = lsnow_lsm
        ! Set lower bound for LSM model, runs from negative (above surface) to surface (zero)
        Model%lsnow_lsm_lbound = -Model%lsnow_lsm+1
        Model%lsnow_lsm_ubound = 0
      end if
    else
      ! Not used by any of the other LSM choices
      Model%lsnow_lsm        = 0
      Model%lsnow_lsm_lbound = 0
      Model%lsnow_lsm_ubound = 0
    end if
    Model%iopt_thcnd       = iopt_thcnd
    Model%ua_phys          = ua_phys
    Model%usemonalb        = usemonalb
    Model%aoasis           = aoasis
    Model%fasdas           = fasdas
    Model%ivegsrc          = ivegsrc
    Model%nvegcat          = nvegcat
    Model%isot             = isot
    Model%nsoilcat         = nsoilcat
    Model%use_ufo          = use_ufo
    Model%exticeden        = exticeden
    if (Model%exticeden .and. &
      (Model%imp_physics /= Model%imp_physics_gfdl .and. Model%imp_physics /= Model%imp_physics_thompson .and. &
       Model%imp_physics /= Model%imp_physics_nssl )) then
      !see GFS_MP_generic_post.F90; exticeden is only compatible with GFDL,
      !Thompson, or NSSL MP
      print *,' Using exticeden = T is only valid when using GFDL, Thompson, or NSSL microphysics.'
      stop
    end if
! GFDL surface layer options
    Model%lcurr_sf         = lcurr_sf
    Model%pert_cd          = pert_cd
    Model%ntsflg           = ntsflg
    Model%sfenth           = sfenth

!--- lake  model parameters
    Model%lkm              = lkm
    Model%iopt_lake        = iopt_lake
    Model%use_lake2m       = use_lake2m
    Model%lakedepth_threshold = lakedepth_threshold
    Model%lakefrac_threshold = lakefrac_threshold

!--- clm lake model parameters
    Model%nlevlake_clm_lake = nlevlake_clm_lake
    Model%nlevsoil_clm_lake = nlevsoil_clm_lake
    Model%nlevsnow_clm_lake = nlevsnow_clm_lake
    Model%nlevsnowsoil_clm_lake = nlevsnowsoil_clm_lake
    Model%nlevsnowsoil1_clm_lake = nlevsnowsoil1_clm_lake
    Model%clm_lake_depth_default = clm_lake_depth_default
    Model%clm_lake_use_lakedepth = clm_lake_use_lakedepth
    Model%clm_lake_debug = clm_lake_debug
    Model%clm_debug_print = clm_debug_print

! Noah MP options from namelist
!
    Model%iopt_dveg        = iopt_dveg
    Model%iopt_crs         = iopt_crs
    Model%iopt_btr         = iopt_btr
    Model%iopt_run         = iopt_run
    Model%iopt_sfc         = iopt_sfc
    Model%iopt_frz         = iopt_frz
    Model%iopt_inf         = iopt_inf
    Model%iopt_rad         = iopt_rad
    Model%iopt_alb         = iopt_alb
    if (Model%lsm==Model%lsm_noahmp .and. Model%exticeden .and. iopt_snf == 4) then
      Model%iopt_snf         = 5
    else
      Model%iopt_snf         = iopt_snf
    end if
    Model%iopt_tbot        = iopt_tbot
    Model%iopt_stc         = iopt_stc
    Model%iopt_trs         = iopt_trs
    Model%iopt_diag        = iopt_diag

! RUC lsm options
    Model%mosaic_lu        = mosaic_lu
    Model%mosaic_soil      = mosaic_soil
    Model%isncond_opt      = isncond_opt
    Model%isncovr_opt      = isncovr_opt
    Model%add_fire_heat_flux = add_fire_heat_flux ! JLS

!--- tuning parameters for physical parameterizations
    Model%ras              = ras
    Model%flipv            = flipv
    Model%trans_trac       = trans_trac
    Model%old_monin        = old_monin
    Model%cnvgwd           = cnvgwd
    Model%mstrat           = mstrat
    Model%moist_adj        = moist_adj
    Model%cscnv            = cscnv
    Model%cal_pre          = cal_pre
    Model%do_aw            = do_aw
    Model%cs_parm          = cs_parm
    Model%do_shoc          = do_shoc
    Model%shoc_parm        = shoc_parm
    Model%shocaftcnv       = shocaftcnv
    Model%shoc_cld         = shoc_cld

!HWRF physics suite
    if (hwrf_samfdeep .and. imfdeepcnv/=2) then
       write(*,*) 'Logic error: hwrf_samfdeep requires imfdeepcnv=2'
       stop
    end if
    if (hwrf_samfshal .and. imfshalcnv/=2) then
       write(*,*) 'Logic error: hwrf_samfshal requires imfshalcnv=2'
       stop
    end if
    Model%hwrf_samfdeep = hwrf_samfdeep
    Model%hwrf_samfshal = hwrf_samfshal

    !--prognostic closure - moisture coupling                                                                                 
    if ((progsigma .and. imfdeepcnv/=2) .and. (progsigma .and. imfdeepcnv/=5)) then
       write(*,*) 'Logic error: progsigma requires imfdeepcnv=2 or 5'
       stop
    end if
    Model%progsigma = progsigma
    Model%betascu = betascu
    Model%betamcu = betamcu
    Model%betadcu = betadcu

    if (oz_phys .and. oz_phys_2015) then
       write(*,*) 'Logic error: can only use one ozone physics option (oz_phys or oz_phys_2015), not both. Exiting.'
       stop
    end if
    Model%oz_phys          = oz_phys
    Model%oz_phys_2015     = oz_phys_2015
    Model%h2o_phys         = h2o_phys

    Model%pdfcld            = pdfcld
    Model%shcnvcw           = shcnvcw
    Model%redrag            = redrag
    Model%hybedmf           = hybedmf
    Model%satmedmf          = satmedmf
    Model%shinhong          = shinhong
    Model%do_ysu            = do_ysu
    Model%dspheat           = dspheat
    Model%hurr_pbl          = hurr_pbl
    Model%lheatstrg         = lheatstrg
    Model%lseaspray         = lseaspray
    Model%cnvcld            = cnvcld
    Model%xr_cnvcld         = xr_cnvcld
    Model%random_clds       = random_clds
    Model%shal_cnv          = shal_cnv
    Model%imfshalcnv        = imfshalcnv
    Model%imfdeepcnv        = imfdeepcnv
    Model%isatmedmf         = isatmedmf
    Model%do_deep           = do_deep
    Model%nmtvr             = nmtvr
    Model%jcap              = jcap
    Model%flgmin            = flgmin
    Model%cgwf              = cgwf
    Model%ccwf              = ccwf
    Model%cdmbgwd           = cdmbgwd
    Model%alpha_fd          = alpha_fd
    Model%psl_gwd_dx_factor = psl_gwd_dx_factor
    Model%sup               = sup
    Model%ctei_rm           = ctei_rm
    Model%crtrh             = crtrh
    Model%dlqf              = dlqf
    Model%psauras           = psauras
    Model%prauras           = prauras
    Model%wminras           = wminras
    Model%rbcr              = rbcr
    Model%do_gwd            = maxval(Model%cdmbgwd) > 0.0 ! flag to restore OGWs of GFS-v15
! OLD GFS-v12-15 conv scheme
    Model%do_cnvgwd         = Model%cnvgwd .and. maxval(Model%cdmbgwd(3:4)) == 0.0
    Model%do_mynnedmf       = do_mynnedmf
    Model%do_mynnsfclay     = do_mynnsfclay
    ! DH* TODO - move to MYNN namelist section
    Model%bl_mynn_cloudpdf  = bl_mynn_cloudpdf
    Model%bl_mynn_mixlength = bl_mynn_mixlength
    Model%bl_mynn_edmf      = bl_mynn_edmf
    Model%bl_mynn_edmf_mom  = bl_mynn_edmf_mom
    Model%bl_mynn_edmf_tke  = bl_mynn_edmf_tke
    Model%bl_mynn_cloudmix  = bl_mynn_cloudmix
    Model%bl_mynn_mixqt     = bl_mynn_mixqt
    Model%bl_mynn_output    = bl_mynn_output
    Model%bl_mynn_tkeadvect = bl_mynn_tkeadvect
    Model%bl_mynn_closure   = bl_mynn_closure
    Model%tke_budget        = tke_budget
    Model%icloud_bl         = icloud_bl
    Model%isftcflx          = isftcflx
    Model%iz0tlnd           = iz0tlnd
    Model%sfclay_compute_flux = sfclay_compute_flux
    Model%sfclay_compute_diag = sfclay_compute_diag
    Model%var_ric           = var_ric
    Model%coef_ric_l        = coef_ric_l
    Model%coef_ric_s        = coef_ric_s
    ! *DH

    Model%gwd_opt           = gwd_opt
    if (Model%gwd_opt==3 .or. Model%gwd_opt==33 .or. &
        Model%gwd_opt==2 .or. Model%gwd_opt==22) then
      ! Add 10 more orographic static fields for GSL drag scheme
      Model%nmtvr = 24
    end if
    Model%do_ugwp_v0           = do_ugwp_v0
    Model%do_ugwp_v0_orog_only = do_ugwp_v0_orog_only
    Model%do_ugwp_v0_nst_only  = do_ugwp_v0_nst_only
    Model%do_gsl_drag_ls_bl    = do_gsl_drag_ls_bl
    Model%do_gsl_drag_ss       = do_gsl_drag_ss
    Model%do_gsl_drag_tofd     = do_gsl_drag_tofd
    Model%do_gwd_opt_psl       = do_gwd_opt_psl
    Model%do_ugwp_v1           = do_ugwp_v1
    Model%do_ugwp_v1_orog_only = do_ugwp_v1_orog_only
    Model%do_ugwp_v1_w_gsldrag = do_ugwp_v1_w_gsldrag
!
! consistency in application of the combined ugwp-v1 and gsldrag
!
    if ( Model%do_ugwp_v1_w_gsldrag) then
       if(Model%gwd_opt == 1 )then
          Model%gwd_opt =2
	  Model%nmtvr = 24
       endif
       Model%do_gsl_drag_ls_bl    = .true.
       Model%do_gsl_drag_tofd     = .true.
       Model%do_gsl_drag_ss       = .true.
       Model%do_ugwp_v1_orog_only = .false.
       Model%do_gwd_opt_psl       = .true.
    endif

    Model%do_myjsfc            = do_myjsfc
    Model%do_myjpbl            = do_myjpbl

!--- Rayleigh friction
    Model%prslrd0          = prslrd0
    Model%ral_ts           = ral_ts

!--- mass flux deep convection
    Model%clam_deep        = clam_deep
    Model%c0s_deep         = c0s_deep
    Model%c1_deep          = c1_deep
    Model%betal_deep       = betal_deep
    Model%betas_deep       = betas_deep
    Model%evef             = evef
    Model%evfact_deep      = evfact_deep
    Model%evfactl_deep     = evfactl_deep
    Model%pgcon_deep       = pgcon_deep
    Model%asolfac_deep     = asolfac_deep

!--- mass flux shallow convection
    Model%clam_shal        = clam_shal
    Model%c0s_shal         = c0s_shal
    Model%c1_shal          = c1_shal
    Model%pgcon_shal       = pgcon_shal
    Model%asolfac_shal     = asolfac_shal

!--- near surface sea temperature model
    Model%nst_anl          = nst_anl
    Model%lsea             = lsea
    Model%nstf_name        = nstf_name

!--- fractional grid
    Model%frac_grid        = frac_grid
    Model%frac_ice         = frac_ice
    Model%ignore_lake      = ignore_lake
    Model%min_lakeice      = min_lakeice
    Model%min_seaice       = min_seaice
    Model%min_lake_height  = min_lake_height
    Model%rho_h2o          = rho_h2o

!--- surface layer
    Model%sfc_z0_type      = sfc_z0_type
    if (Model%cplwav2atm) Model%sfc_z0_type = -1
    Model%icplocn2atm      = icplocn2atm

!--- potential temperature reference in sfc layer
    Model%thsfc_loc        = thsfc_loc
!--- flux method in 2-m diagnostics
    Model%diag_flux        = diag_flux
!--- flux method in 2-m diagnostics (for stable conditions)
    Model%diag_log         = diag_log

!--- vertical diffusion
    Model%xkzm_m           = xkzm_m
    Model%xkzm_h           = xkzm_h
    Model%xkzm_s           = xkzm_s
    Model%xkzminv          = xkzminv
    Model%moninq_fac       = moninq_fac
    Model%dspfac           = dspfac
    Model%bl_upfr          = bl_upfr
    Model%bl_dnfr          = bl_dnfr
    Model%rlmx             = rlmx
    Model%elmx             = elmx
    Model%sfc_rlm          = sfc_rlm
    Model%tc_pbl           = tc_pbl

!--- canopy heat storage parametrization
    Model%h0facu           = h0facu
    Model%h0facs           = h0facs

!--- stochastic physics options
    ! do_sppt, do_shum, do_skeb and lndp_type are namelist variables in group
    ! physics that are parsed here and then compared in init_stochastic_physics
    ! to the stochastic physics namelist parametersto ensure consistency.
    Model%do_sppt          = do_sppt
    Model%pert_mp          = pert_mp
    Model%pert_clds        = pert_clds
    Model%pert_radtend     = pert_radtend
    Model%use_zmtnblck     = use_zmtnblck
    Model%do_shum          = do_shum
    Model%do_skeb          = do_skeb
    !--- stochastic surface perturbation options
    Model%lndp_type        = lndp_type
    Model%n_var_lndp       = n_var_lndp
    Model%lndp_each_step   = lndp_each_step
    Model%do_spp           = do_spp
    Model%n_var_spp        = n_var_spp

    if (Model%lndp_type/=0) then
      allocate(Model%lndp_var_list(Model%n_var_lndp))
      allocate(Model%lndp_prt_list(Model%n_var_lndp))
      Model%lndp_var_list(:) = ''
      Model%lndp_prt_list(:) = clear_val
    end if

    if (Model%do_spp) then
      allocate(Model%spp_var_list(Model%n_var_spp))
      allocate(Model%spp_prt_list(Model%n_var_spp))
      allocate(Model%spp_stddev_cutoff(Model%n_var_spp))
      Model%spp_var_list(:) = ''
      Model%spp_prt_list(:) = clear_val
      Model%spp_stddev_cutoff(:) = clear_val
    end if

    !--- cellular automata options
    ! force namelist constsitency
    allocate(Model%vfact_ca(levs))
    if ( .not. ca_global ) nca_g=0
    if ( .not. ca_sgs ) nca=0

    Model%nca              = nca
    Model%ncells           = ncells
    Model%nlives           = nlives
    Model%nca_g            = nca_g
    Model%ncells_g         = ncells_g
    Model%nlives_g         = nlives_g
    Model%nfracseed        = nfracseed
    Model%nseed            = nseed
    Model%nseed_g          = nseed_g
    Model%ca_global        = ca_global
    Model%do_ca            = do_ca
    Model%ca_advect        = ca_advect
    Model%ca_sgs           = ca_sgs
    Model%iseed_ca         = iseed_ca
    Model%ca_smooth        = ca_smooth
    Model%nspinup          = nspinup
    Model%nthresh          = nthresh
    Model%ca_amplitude     = ca_amplitude
    Model%nsmooth          = nsmooth
    Model%ca_closure       = ca_closure
    Model%ca_entr          = ca_entr
    Model%ca_trigger       = ca_trigger

    ! IAU flags
    !--- iau parameters
    Model%iaufhrs         = iaufhrs
    Model%iau_inc_files   = iau_inc_files
    Model%iau_delthrs     = iau_delthrs
    Model%iau_filter_increments = iau_filter_increments
    Model%iau_drymassfixer = iau_drymassfixer
    if(Model%me==0) print *,' model init,iaufhrs=',Model%iaufhrs

!--- debug flags
    Model%debug            = debug
    Model%pre_rad          = pre_rad
    Model%print_diff_pgr   = print_diff_pgr

!--- tracer handling
    Model%ntrac            = size(tracer_names)
    Model%ntracp1          = Model%ntrac + 1
    Model%ntracp100        = Model%ntrac + 100
    allocate (Model%tracer_names(Model%ntrac))
    Model%tracer_names(:)  = tracer_names(:)
    Model%ntqv             = 1
#ifdef MULTI_GASES
    Model%nto              = get_tracer_index(Model%tracer_names, 'spo',        Model%me, Model%master, Model%debug)
    Model%nto2             = get_tracer_index(Model%tracer_names, 'spo2',       Model%me, Model%master, Model%debug)
    Model%ntoz             = get_tracer_index(Model%tracer_names, 'spo3',       Model%me, Model%master, Model%debug)
#else
    Model%ntoz             = get_tracer_index(Model%tracer_names, 'o3mr',       Model%me, Model%master, Model%debug)
    if( Model%ntoz <= 0 )  &
    Model%ntoz             =  get_tracer_index(Model%tracer_names, 'spo3',       Model%me, Model%master, Model%debug)
#endif
    Model%ntcw             = get_tracer_index(Model%tracer_names, 'liq_wat',    Model%me, Model%master, Model%debug)
    Model%ntiw             = get_tracer_index(Model%tracer_names, 'ice_wat',    Model%me, Model%master, Model%debug)
    Model%ntrw             = get_tracer_index(Model%tracer_names, 'rainwat',    Model%me, Model%master, Model%debug)
    Model%ntsw             = get_tracer_index(Model%tracer_names, 'snowwat',    Model%me, Model%master, Model%debug)
    Model%ntgl             = get_tracer_index(Model%tracer_names, 'graupel',    Model%me, Model%master, Model%debug)
    Model%nthl             = get_tracer_index(Model%tracer_names, 'hailwat',    Model%me, Model%master, Model%debug)
    Model%ntclamt          = get_tracer_index(Model%tracer_names, 'cld_amt',    Model%me, Model%master, Model%debug)
    Model%ntlnc            = get_tracer_index(Model%tracer_names, 'water_nc',   Model%me, Model%master, Model%debug)
    Model%ntinc            = get_tracer_index(Model%tracer_names, 'ice_nc',     Model%me, Model%master, Model%debug)
    Model%ntrnc            = get_tracer_index(Model%tracer_names, 'rain_nc',    Model%me, Model%master, Model%debug)
    Model%ntsnc            = get_tracer_index(Model%tracer_names, 'snow_nc',    Model%me, Model%master, Model%debug)
    Model%ntgnc            = get_tracer_index(Model%tracer_names, 'graupel_nc', Model%me, Model%master, Model%debug)
    Model%nthnc            = get_tracer_index(Model%tracer_names, 'hail_nc',    Model%me, Model%master, Model%debug)
    Model%ntccn            = get_tracer_index(Model%tracer_names, 'ccn_nc',     Model%me, Model%master, Model%debug)
    Model%ntccna           = get_tracer_index(Model%tracer_names, 'ccna_nc',    Model%me, Model%master, Model%debug)
    Model%ntgv             = get_tracer_index(Model%tracer_names, 'graupel_vol',Model%me, Model%master, Model%debug)
    Model%nthv             = get_tracer_index(Model%tracer_names, 'hail_vol',   Model%me, Model%master, Model%debug)
    Model%ntrz             = get_tracer_index(Model%tracer_names, 'rain_ref',   Model%me, Model%master, Model%debug)
    Model%ntgz             = get_tracer_index(Model%tracer_names, 'graupel_ref',Model%me, Model%master, Model%debug)
    Model%nthz             = get_tracer_index(Model%tracer_names, 'hail_ref',   Model%me, Model%master, Model%debug)
    Model%ntke             = get_tracer_index(Model%tracer_names, 'sgs_tke',    Model%me, Model%master, Model%debug)
    Model%ntsigma          = get_tracer_index(Model%tracer_names, 'sigmab',     Model%me, Model%master, Model%debug)
    Model%nqrimef          = get_tracer_index(Model%tracer_names, 'q_rimef',    Model%me, Model%master, Model%debug)
    Model%ntwa             = get_tracer_index(Model%tracer_names, 'liq_aero',   Model%me, Model%master, Model%debug)
    Model%ntia             = get_tracer_index(Model%tracer_names, 'ice_aero',   Model%me, Model%master, Model%debug)
    if (Model%cpl_fire) then
    Model%ntfsmoke         = get_tracer_index(Model%tracer_names, 'fsmoke',   Model%me, Model%master, Model%debug)
    endif
    if (Model%rrfs_sd) then
    Model%ntsmoke          = get_tracer_index(Model%tracer_names, 'smoke',      Model%me, Model%master, Model%debug)
    Model%ntdust           = get_tracer_index(Model%tracer_names, 'dust',       Model%me, Model%master, Model%debug)
    Model%ntcoarsepm       = get_tracer_index(Model%tracer_names, 'coarsepm',   Model%me, Model%master, Model%debug)
    endif
!--- initialize parameters for atmospheric chemistry tracers
    call Model%init_chemistry(tracer_types)

!--- setup aerosol scavenging factors
    call Model%init_scavenging(fscav_aero)

    ! Tracer diagnostics indices and dimension size, which must be in
    ! Model to be forwarded to the right places.

    ! Individual processes:
    Model%index_of_process_pbl = 1
    Model%index_of_process_dcnv = 2
    Model%index_of_process_scnv = 3
    Model%index_of_process_mp = 4
    Model%index_of_process_prod_loss = 5
    Model%index_of_process_ozmix = 6
    Model%index_of_process_temp = 7
    Model%index_of_process_overhead_ozone = 8
    Model%index_of_process_longwave = 9
    Model%index_of_process_shortwave = 10
    Model%index_of_process_orographic_gwd = 11
    Model%index_of_process_rayleigh_damping = 12
    Model%index_of_process_nonorographic_gwd = 13
    Model%index_of_process_conv_trans = 14
    Model%index_of_process_dfi_radar = 15

    ! Number of processes to sum (last index of prior set)
    Model%nprocess_summed = Model%index_of_process_dfi_radar

    ! Sums of other processes, which must be after nprocess_summed:
    Model%index_of_process_physics = Model%nprocess_summed+1
    Model%index_of_process_non_physics = Model%nprocess_summed+2
    Model%index_of_process_photochem = Model%nprocess_summed+3

    ! Total number of processes (last index of prior set)
    Model%nprocess = Model%index_of_process_photochem

    ! List which processes should be summed as photochemical:
    allocate(Model%is_photochem(Model%nprocess))
    Model%is_photochem = .false.
    Model%is_photochem(Model%index_of_process_prod_loss) = .true.
    Model%is_photochem(Model%index_of_process_ozmix) = .true.
    Model%is_photochem(Model%index_of_process_temp) = .true.
    Model%is_photochem(Model%index_of_process_overhead_ozone) = .true.

    ! Non-tracers that appear in first dimension of dtidx:
    Model%index_of_temperature = 10
    Model%index_of_x_wind = 11
    Model%index_of_y_wind = 12

    ! Last index of outermost dimension of dtend
    Model%ndtend = 0
    allocate(Model%dtidx(Model%ntracp100,Model%nprocess))
    Model%dtidx = -99

    if(Model%ntchm>0) then
      Model%ntdu1 = get_tracer_index(Model%tracer_names, 'dust1', Model%me, Model%master, Model%debug)
      Model%ntdu2 = get_tracer_index(Model%tracer_names, 'dust2', Model%me, Model%master, Model%debug)
      Model%ntdu3 = get_tracer_index(Model%tracer_names, 'dust3', Model%me, Model%master, Model%debug)
      Model%ntdu4 = get_tracer_index(Model%tracer_names, 'dust4', Model%me, Model%master, Model%debug)
      Model%ntdu5 = get_tracer_index(Model%tracer_names, 'dust5', Model%me, Model%master, Model%debug)
      Model%ntss1 = get_tracer_index(Model%tracer_names, 'seas1', Model%me, Model%master, Model%debug)
      Model%ntss2 = get_tracer_index(Model%tracer_names, 'seas2', Model%me, Model%master, Model%debug)
      Model%ntss3 = get_tracer_index(Model%tracer_names, 'seas3', Model%me, Model%master, Model%debug)
      Model%ntss4 = get_tracer_index(Model%tracer_names, 'seas4', Model%me, Model%master, Model%debug)
      Model%ntss5 = get_tracer_index(Model%tracer_names, 'seas5', Model%me, Model%master, Model%debug)
      Model%ntsu  = get_tracer_index(Model%tracer_names, 'so4',   Model%me, Model%master, Model%debug)
      Model%ntbcb = get_tracer_index(Model%tracer_names, 'bc1',   Model%me, Model%master, Model%debug)
      Model%ntbcl = get_tracer_index(Model%tracer_names, 'bc2',   Model%me, Model%master, Model%debug)
      Model%ntocb = get_tracer_index(Model%tracer_names, 'oc1',   Model%me, Model%master, Model%debug)
      Model%ntocl = get_tracer_index(Model%tracer_names, 'oc2',   Model%me, Model%master, Model%debug)
    end if

    ! Lake & fractional grid safety checks
    if(Model%me==Model%master) then
      if(Model%lkm>0 .and. Model%frac_grid) then
        write(0,*) 'WARNING: Lake fractional grid support is experimental. Use at your own risk!'
      else if(Model%lkm>0 .and. Model%iopt_lake==Model%iopt_lake_clm .and. .not. Model%frac_ice) then
        write(0,*) 'WARNING: CLM Lake Model will not work without frac_ice=.true.'
      endif
      if(Model%lkm==2) then
        write(0,*) 'WARNING: Running both lake and nsst on lake points is experimental. Use at your own risk!'
      endif
    endif

    if(ldiag3d) then
       ! Flags used to turn on or off tracer "causes"
       have_pbl_edmf = Model%hybedmf .or. Model%satmedmf .or. Model%do_mynnedmf
       have_samf =  Model%satmedmf .or. Model%trans_trac .or. Model%ras .or. Model%do_shoc
       have_pbl = .true.
       have_dcnv = Model%imfdeepcnv>0 !Model%ras .or. Model%cscnv .or. Model%do_deep .or. Model%hwrf_samfdeep
       have_scnv = Model%imfshalcnv>0 !Model%shal_cnv
       have_mp = Model%imp_physics>0
       have_oz_phys = Model%oz_phys .or. Model%oz_phys_2015

       ! Rayleigh damping flag must match logic in rayleigh_damp.f
       have_rdamp = .not. (Model%lsidea .or. Model%ral_ts <= 0.0 .or. Model%prslrd0 == 0.0)

       ! have_cnvtrans flag must match logic elsewhere in GFS_typedefs and suite interstitials.
       have_cnvtrans = (have_dcnv .or. have_scnv) .and. &
            (cscnv .or. satmedmf .or. trans_trac .or. ras) &
            .and. Model%flag_for_scnv_generic_tend &
            .and. Model%flag_for_dcnv_generic_tend

       ! Increment idtend and fill dtidx:

        call allocate_dtend_labels_and_causes(Model)

        ! Default names of tracers just in case later code does not initialize them:
        do itrac=1,Model%ntrac
           write(namestr,'("tracer",I0)') itrac
           write(descstr,'("tracer ",I0," of ",I0)') itrac, Model%ntrac
           call label_dtend_tracer(Model,100+itrac,trim(namestr),trim(descstr),'kg kg-1 s-1')
        enddo

        if(Model%ntchs>0) then
           if(Model%ntchm>0) then
              ! Chemical tracers are first so more specific tracer names
              ! replace them. There is no straightforward way of getting
              ! chemical tracer short names or descriptions, so we use
              ! indices instead.
              do ichem=Model%ntchs,Model%ntchs+Model%ntchm-1
                 write(namestr,'("chem",I0)') ichem
                 write(descstr,'("chemical tracer ",I0," of ",I0)') ichem, Model%ntchm
                 call label_dtend_tracer(Model,100+ichem,trim(namestr),trim(descstr),'kg kg-1 s-1')
              enddo
           endif

           ! More specific chemical tracer names:
           call label_dtend_tracer(Model,100+Model%ntchs,'so2','sulfur dioxide concentration','kg kg-1 s-1')
           if(Model%ntchm>0) then
              ! Need better descriptions of these.
              call label_dtend_tracer(Model,100+Model%ntchm+Model%ntchs-1,'pp10','pp10 concentration','kg kg-1 s-1')

              itrac=get_tracer_index(Model%tracer_names, 'DMS', Model%me, Model%master, Model%debug)
              if(itrac>0) then
                 call label_dtend_tracer(Model,100+itrac,'DMS','DMS concentration','kg kg-1 s-1')
              endif
              itrac=get_tracer_index(Model%tracer_names, 'msa', Model%me, Model%master, Model%debug)
              if(itrac>0) then
                 call label_dtend_tracer(Model,100+itrac,'msa','msa concentration','kg kg-1 s-1')
              endif
           endif
        endif

        call label_dtend_tracer(Model,Model%index_of_temperature,'temp','temperature','K s-1')
        call label_dtend_tracer(Model,Model%index_of_x_wind,'u','x wind','m s-2')
        call label_dtend_tracer(Model,Model%index_of_y_wind,'v','y wind','m s-2')

        ! Other tracer names. These were taken from GFS_typedefs.F90 with descriptions from GFS_typedefs.meta
        call label_dtend_tracer(Model,100+Model%ntqv,'qv','water vapor specific humidity','kg kg-1 s-1')
        call label_dtend_tracer(Model,100+Model%ntoz,'o3','ozone concentration','kg kg-1 s-1')
        call label_dtend_tracer(Model,100+Model%ntcw,'liq_wat','cloud condensate (or liquid water)','kg kg-1 s-1')
        call label_dtend_tracer(Model,100+Model%ntiw,'ice_wat','ice water','kg kg-1 s-1')
        call label_dtend_tracer(Model,100+Model%ntrw,'rainwat','rain water','kg kg-1 s-1')
        call label_dtend_tracer(Model,100+Model%ntsw,'snowwat','snow water','kg kg-1 s-1')
        call label_dtend_tracer(Model,100+Model%ntgl,'graupel','graupel','kg kg-1 s-1')
        call label_dtend_tracer(Model,100+Model%nthl,'hailwat','hail','kg kg-1 s-1')
        call label_dtend_tracer(Model,100+Model%ntclamt,'cld_amt','cloud amount integer','kg kg-1 s-1')
        call label_dtend_tracer(Model,100+Model%ntlnc,'water_nc','liquid number concentration','kg-1 s-1')
        call label_dtend_tracer(Model,100+Model%ntinc,'ice_nc','ice number concentration','kg-1 s-1')
        call label_dtend_tracer(Model,100+Model%ntrnc,'rain_nc','rain number concentration','kg-1 s-1')
        call label_dtend_tracer(Model,100+Model%ntsnc,'snow_nc','snow number concentration','kg-1 s-1')
        call label_dtend_tracer(Model,100+Model%ntgnc,'graupel_nc','graupel number concentration','kg-1 s-1')
        call label_dtend_tracer(Model,100+Model%nthnc,'hail_nc','hail number concentration','kg-1 s-1')
        call label_dtend_tracer(Model,100+Model%ntccn,'ccn_nc','CCN number concentration','kg-1 s-1')
        call label_dtend_tracer(Model,100+Model%ntgv,'graupel_vol','graupel volume','m3 kg-1 s-1')
        call label_dtend_tracer(Model,100+Model%nthv,'hail_vol','hail volume','m3 kg-1 s-1')
        call label_dtend_tracer(Model,100+Model%ntrz,'rain_ref','rain reflectivity','m3 kg-1 s-1')
        call label_dtend_tracer(Model,100+Model%ntgz,'graupel_ref','graupel reflectivity','m3 kg-1 s-1')
        call label_dtend_tracer(Model,100+Model%nthz,'hail_ref','hail reflectivity','m3 kg-1 s-1')
        call label_dtend_tracer(Model,100+Model%ntke,'sgs_tke','turbulent kinetic energy','J s-1')
        call label_dtend_tracer(Model,100+Model%nqrimef,'q_rimef','mass weighted rime factor','kg-1 s-1')
        call label_dtend_tracer(Model,100+Model%ntwa,'liq_aero','number concentration of water-friendly aerosols','kg-1 s-1')
        call label_dtend_tracer(Model,100+Model%ntia,'ice_aero','number concentration of ice-friendly aerosols','kg-1 s-1')
        call label_dtend_tracer(Model,100+Model%nto,'o_ion','oxygen ion concentration','kg kg-1 s-1')
        call label_dtend_tracer(Model,100+Model%nto2,'o2','oxygen concentration','kg kg-1 s-1')

        call label_dtend_cause(Model,Model%index_of_process_pbl,'pbl','tendency due to PBL')
        call label_dtend_cause(Model,Model%index_of_process_dcnv,'deepcnv','tendency due to deep convection')
        call label_dtend_cause(Model,Model%index_of_process_scnv,'shalcnv','tendency due to shallow convection')
        call label_dtend_cause(Model,Model%index_of_process_mp,'mp','tendency due to microphysics')
        call label_dtend_cause(Model,Model%index_of_process_prod_loss,'prodloss','tendency due to production and loss rate')
        call label_dtend_cause(Model,Model%index_of_process_ozmix,'o3mix','tendency due to ozone mixing ratio')
        call label_dtend_cause(Model,Model%index_of_process_temp,'temp','tendency due to temperature')
        call label_dtend_cause(Model,Model%index_of_process_overhead_ozone,'o3column','tendency due to overhead ozone column')
        call label_dtend_cause(Model,Model%index_of_process_dfi_radar,'dfi_radar','tendency due to dfi radar mp temperature forcing')
        call label_dtend_cause(Model,Model%index_of_process_photochem,'photochem','tendency due to photochemical processes')
        call label_dtend_cause(Model,Model%index_of_process_physics,'phys','tendency due to physics')
        call label_dtend_cause(Model,Model%index_of_process_non_physics,'nophys','tendency due to non-physics processes', &
                               mod_name='gfs_dyn')
        call label_dtend_cause(Model,Model%index_of_process_conv_trans,'cnvtrans','tendency due to convective transport')
        call label_dtend_cause(Model,Model%index_of_process_longwave,'lw','tendency due to long wave radiation')
        call label_dtend_cause(Model,Model%index_of_process_shortwave,'sw','tendency due to short wave radiation')
        call label_dtend_cause(Model,Model%index_of_process_orographic_gwd,'orogwd','tendency due to orographic gravity wave drag')
        call label_dtend_cause(Model,Model%index_of_process_rayleigh_damping,'rdamp','tendency due to Rayleigh damping')
        call label_dtend_cause(Model,Model%index_of_process_nonorographic_gwd,'cnvgwd','tendency due to convective gravity wave drag')

       call fill_dtidx(Model,dtend_select,Model%index_of_temperature,Model%index_of_process_longwave)
       call fill_dtidx(Model,dtend_select,Model%index_of_temperature,Model%index_of_process_shortwave)
       call fill_dtidx(Model,dtend_select,Model%index_of_temperature,Model%index_of_process_pbl,have_pbl)
       call fill_dtidx(Model,dtend_select,Model%index_of_temperature,Model%index_of_process_dcnv,have_dcnv)
       call fill_dtidx(Model,dtend_select,Model%index_of_temperature,Model%index_of_process_scnv,have_scnv)
       call fill_dtidx(Model,dtend_select,Model%index_of_temperature,Model%index_of_process_mp,have_mp)
       call fill_dtidx(Model,dtend_select,Model%index_of_temperature,Model%index_of_process_dfi_radar,have_mp .and. Model%num_dfi_radar>0)
       call fill_dtidx(Model,dtend_select,Model%index_of_temperature,Model%index_of_process_orographic_gwd)
       call fill_dtidx(Model,dtend_select,Model%index_of_temperature,Model%index_of_process_rayleigh_damping,have_rdamp)
       call fill_dtidx(Model,dtend_select,Model%index_of_temperature,Model%index_of_process_nonorographic_gwd)
       call fill_dtidx(Model,dtend_select,Model%index_of_temperature,Model%index_of_process_physics)
       call fill_dtidx(Model,dtend_select,Model%index_of_temperature,Model%index_of_process_non_physics)

       call fill_dtidx(Model,dtend_select,Model%index_of_x_wind,Model%index_of_process_pbl,have_pbl)
       call fill_dtidx(Model,dtend_select,Model%index_of_y_wind,Model%index_of_process_pbl,have_pbl)
       call fill_dtidx(Model,dtend_select,Model%index_of_x_wind,Model%index_of_process_orographic_gwd)
       call fill_dtidx(Model,dtend_select,Model%index_of_y_wind,Model%index_of_process_orographic_gwd)
       call fill_dtidx(Model,dtend_select,Model%index_of_x_wind,Model%index_of_process_dcnv,have_dcnv)
       call fill_dtidx(Model,dtend_select,Model%index_of_y_wind,Model%index_of_process_dcnv,have_dcnv)
       call fill_dtidx(Model,dtend_select,Model%index_of_x_wind,Model%index_of_process_nonorographic_gwd)
       call fill_dtidx(Model,dtend_select,Model%index_of_y_wind,Model%index_of_process_nonorographic_gwd)
       call fill_dtidx(Model,dtend_select,Model%index_of_x_wind,Model%index_of_process_rayleigh_damping,have_rdamp)
       call fill_dtidx(Model,dtend_select,Model%index_of_y_wind,Model%index_of_process_rayleigh_damping,have_rdamp)
       call fill_dtidx(Model,dtend_select,Model%index_of_x_wind,Model%index_of_process_scnv,have_scnv)
       call fill_dtidx(Model,dtend_select,Model%index_of_y_wind,Model%index_of_process_scnv,have_scnv)
       call fill_dtidx(Model,dtend_select,Model%index_of_x_wind,Model%index_of_process_physics)
       call fill_dtidx(Model,dtend_select,Model%index_of_y_wind,Model%index_of_process_physics)
       call fill_dtidx(Model,dtend_select,Model%index_of_x_wind,Model%index_of_process_non_physics)
       call fill_dtidx(Model,dtend_select,Model%index_of_y_wind,Model%index_of_process_non_physics)

       if(qdiag3d) then
          call fill_dtidx(Model,dtend_select,100+Model%ntqv,Model%index_of_process_scnv,have_scnv)
          call fill_dtidx(Model,dtend_select,100+Model%ntqv,Model%index_of_process_dcnv,have_dcnv)

          if(have_cnvtrans) then
             do itrac=2,Model%ntrac
                if(itrac==Model%ntchs) exit ! remaining tracers are chemical
                if ( itrac /= Model%ntcw  .and. itrac /= Model%ntiw  .and. itrac /= Model%ntclamt .and. &
                     itrac /= Model%ntrw  .and. itrac /= Model%ntsw  .and. itrac /= Model%ntrnc   .and. &
                     itrac /= Model%ntsnc .and. itrac /= Model%ntgl  .and. itrac /= Model%ntgnc   .and. &
                     itrac /= Model%nthl  .and. itrac /= Model%nthnc .and. itrac /= Model%nthv    .and. &
                     itrac /= Model%ntgv  .and. itrac /= Model%ntrz  .and. itrac /= Model%ntgz    .and. &
                     itrac /= Model%nthz ) then
                   call fill_dtidx(Model,dtend_select,100+itrac,Model%index_of_process_scnv,have_scnv)
                   call fill_dtidx(Model,dtend_select,100+itrac,Model%index_of_process_dcnv,have_dcnv)
                else if(Model%ntchs<=0 .or. itrac<Model%ntchs) then
                   call fill_dtidx(Model,dtend_select,100+itrac,Model%index_of_process_conv_trans)
                endif
                if ( itrac == Model%ntlnc .or. itrac == Model%ntinc .or. itrac == Model%ntrnc .or. &
                     itrac == Model%ntsnc .or. itrac == Model%ntgnc .or. itrac == Model%nthnc .or. &
                     itrac == Model%nqrimef) then
                   call fill_dtidx(Model,dtend_select,100+itrac,Model%index_of_process_conv_trans)
                endif
             enddo
          else if(have_scnv .or. have_dcnv) then
             ! Scheme does its own tendency reporting, or does not use convective transport.
             do itrac=2,Model%ntrac
                if(itrac==Model%ntchs) exit ! remaining tracers are chemical
                call fill_dtidx(Model,dtend_select,100+itrac,Model%index_of_process_scnv,have_scnv)
                call fill_dtidx(Model,dtend_select,100+itrac,Model%index_of_process_dcnv,have_dcnv)
             enddo
          endif

          call fill_dtidx(Model,dtend_select,100+Model%ntoz,Model%index_of_process_pbl,have_pbl)
          call fill_dtidx(Model,dtend_select,100+Model%ntoz,Model%index_of_process_prod_loss,have_oz_phys)
          call fill_dtidx(Model,dtend_select,100+Model%ntoz,Model%index_of_process_ozmix,have_oz_phys)
          call fill_dtidx(Model,dtend_select,100+Model%ntoz,Model%index_of_process_temp,have_oz_phys)
          call fill_dtidx(Model,dtend_select,100+Model%ntoz,Model%index_of_process_overhead_ozone,have_oz_phys)
          call fill_dtidx(Model,dtend_select,100+Model%ntoz,Model%index_of_process_photochem,have_oz_phys)
          call fill_dtidx(Model,dtend_select,100+Model%ntoz,Model%index_of_process_physics,.true.)
          call fill_dtidx(Model,dtend_select,100+Model%ntoz,Model%index_of_process_non_physics,.true.)

          if(.not.Model%do_mynnedmf .and. .not. Model%satmedmf) then
            call fill_dtidx(Model,dtend_select,100+Model%ntqv,Model%index_of_process_pbl,have_pbl)
            call fill_dtidx(Model,dtend_select,100+Model%ntcw,Model%index_of_process_pbl,have_pbl)
            call fill_dtidx(Model,dtend_select,100+Model%ntiw,Model%index_of_process_pbl,have_pbl)
            call fill_dtidx(Model,dtend_select,100+Model%ntke,Model%index_of_process_pbl,have_pbl)
          endif

          do itrac=1,Model%ntrac
             if(itrac==Model%ntchs) exit ! remaining tracers are chemical
             if(itrac==Model%ntoz) cycle ! already took care of ozone
             if(Model%do_mynnedmf .or. Model%satmedmf) then
               call fill_dtidx(Model,dtend_select,100+itrac,Model%index_of_process_pbl,have_pbl)
             endif
             call fill_dtidx(Model,dtend_select,100+itrac,Model%index_of_process_physics,.true.)
             call fill_dtidx(Model,dtend_select,100+itrac,Model%index_of_process_non_physics,.true.)
             if(itrac/=Model%ntke) then
               call fill_dtidx(Model,dtend_select,100+itrac,Model%index_of_process_mp,have_mp)
             endif
          enddo
       endif
    end if

    IF ( Model%imp_physics == Model%imp_physics_nssl2mccn ) THEN ! recognize this option for compatibility
       Model%imp_physics = Model%imp_physics_nssl
       nssl_ccn_on = .true.
       Model%nssl_ccn_on = .true.
    ENDIF
    IF ( Model%imp_physics == Model%imp_physics_nssl ) THEN !{
        ! check if field_table and nssl_ccn_on flag are consistent
        ! Alternatively could simply rely on the field_table to set values for ccn_on and hail_on
        IF ( .not. Model%nssl_ccn_on ) THEN
          if (Model%me == Model%master) write(*,*) 'NSSL micro: CCN is OFF'
            IF ( Model%ntccn > 1 ) THEN
              IF (Model%me == Model%master) then
               write(*,*) 'NSSL micro: error! CCN is OFF (nssl_ccn_on = F) but ntccn > 1.'
               write(*,*) 'Should  either remove ccn_nc from field_table or set nssl_ccn_on = .true.'
               write(0,*) 'NSSL micro: error! CCN is OFF (nssl_ccn_on = F) but ntccn > 1.'
               write(0,*) 'Should  either remove ccn_nc from field_table or set nssl_ccn_on = .true.'
              ENDIF
              stop
            ENDIF
          Model%ntccn = -99
          Model%ntccna = -99
        ELSEIF ( Model%ntccn < 1 ) THEN
          if (Model%me == Model%master) then
            write(*,*) 'NSSL micro: error! CCN is ON but ntccn < 1. Must have ccn_nc in field_table if nssl_ccn_on=T'
            write(0,*) 'NSSL micro: error! CCN is ON but ntccn < 1. Must have ccn_nc in field_table if nssl_ccn_on=T'
          ENDIF
          stop
        ELSE
          if (Model%me == Model%master) then
            write(*,*) 'NSSL micro: CCN is ON'
          ENDIF
          IF ( Model%ntccna > 1 .and. Model%me == Model%master ) THEN
            write(*,*) 'NSSL micro: CCNA is ON'
          ENDIF
        ENDIF

      if (Model%me == Model%master) then
        write(*,*) 'Model%nthl = ',Model%nthl
      ENDIF
      IF ( ( Model%nthl < 1 ) ) THEN ! check if hail is in the field_table. If not, set flag so the microphysics knows.
        if (Model%me == Model%master) then
          write(*,*) 'NSSL micro: hail is OFF'
          IF ( nssl_hail_on ) write(*,*) 'Namelist had nssl_hail_on=true, but tracer config does not have hailwat'
        ENDIF
        nssl_hail_on = .false.
        Model%nssl_hail_on = .false.
        ! pretend that hail exists so that bad arrays are not passed to microphysics
!         Model%nthl =  Max( 1, Model%ntgl )
!         Model%nthv =  Max( 1, Model%ntgv )
!         Model%nthnc = Max( 1, Model%ntgnc )
      ELSE
        nssl_hail_on = .true.
        Model%nssl_hail_on = .true.
        if (Model%me == Model%master) then
          write(*,*) 'NSSL micro: hail is ON'
          IF ( .not. nssl_hail_on ) write(*,*) 'Namelist had nssl_hail_on=false, but tracer config has hailwat'
        ENDIF
        IF ( Model%nthv < 1 .or. Model%nthnc < 1 ) THEN
           if (Model%me == Model%master) THEN
             write(0,*) 'missing needed tracers for NSSL hail! nthl > 1 but either volume or number is not in field_table'
             write(0,*) 'nthv, nthnc = ', Model%nthv, Model%nthnc
           ENDIF
           stop
        ENDIF
      ENDIF

      Model%nssl_hail_on  = nssl_hail_on

      IF ( ( Model%ntccn < 1 ) ) THEN ! check if ccn is in the field_table. If not, set flag so the microphysics knows.
        if (Model%me == Model%master) then
          write(*,*) 'NSSL micro: CCN is OFF'
        ENDIF
        nssl_ccn_on = .false.
        Model%nssl_ccn_on = .false.
      ELSE
        nssl_ccn_on = .true.
        Model%nssl_ccn_on = .true.
        if (Model%me == Model%master) then
          write(*,*) 'NSSL micro: CCN is ON'
        ENDIF
      ENDIF
      
      ! add checks for nssl_3moment
      IF ( ( Model%nssl_3moment ) ) THEN 
        IF ( Model%ntrz < 1 ) THEN
          write(*,*) 'NSSL micro: 3-moment is ON, but rain_ref tracer is missing'
          stop
        ENDIF
        IF ( Model%ntgz < 1 ) THEN
          write(*,*) 'NSSL micro: 3-moment is ON, but graupel_ref tracer is missing'
          stop
        ENDIF
        IF ( nssl_hail_on ) THEN
        IF ( Model%nthz < 1 ) THEN
          write(*,*) 'NSSL micro: 3-moment is ON, but hail_ref tracer is missing'
          stop
        ENDIF
        ENDIF
      ENDIF

        IF ( Model%ntgl < 1 .or. Model%ntgv < 1 .or. Model%ntgnc < 1 .or. &
             Model%ntsw < 1 .or. Model%ntsnc < 1 .or. &
             Model%ntrw < 1 .or. Model%ntrnc < 1 .or. &
             Model%ntiw < 1 .or. Model%ntinc < 1 .or. &
             Model%ntcw < 1 .or. Model%ntlnc < 1      &
             ) THEN
          if (Model%me == Model%master)  write(0,*) 'missing needed tracers for NSSL!'
           stop
        ENDIF


    ENDIF !}

    !--- NRL ozone physics
    if (Model%ntoz > 0) then
       ! Load data for ozone physics into DDT ozphys
       err_message    = Model%ozphys%load_o3prog('global_o3prdlos.f77',kozpl)
       Model%levozp   = Model%ozphys%nlev
       Model%oz_coeff = Model%ozphys%ncf

       if (Model%me == Model%master) then
          write(*,*) 'Reading in o3data from global_o3prdlos.f77 '
          write(*,*) '      oz_coeff = ', Model%ozphys%ncf
          write(*,*) '       latsozp = ', Model%ozphys%nlat
          write(*,*) '        levozp = ', Model%ozphys%nlev
          write(*,*) '        timeoz = ', Model%ozphys%ntime
       endif
    else
       !--- Climatological ozone
       err_message = Model%ozphys%load_o3clim('global_o3prdlos.f77',kozc)
    end if

    !--- NRL h2o photochemistry physics.
    if (Model%h2o_phys) then
       ! Load data for h2o photochemistry physics.
       err_message     = Model%h2ophys%load('global_h2oprdlos.f77',kh2oc)
       Model%levh2o    = Model%h2ophys%nlev
       Model%h2o_coeff = Model%h2ophys%ncf
    else
       Model%levh2o    = 1
       Model%h2o_coeff = 1
    end if
    
!--- quantities to be used to derive phy_f*d totals
    Model%nshoc_2d         = nshoc_2d
    Model%nshoc_3d         = nshoc_3d
    Model%ncnvcld3d        = ncnvcld3d
    Model%nctp             = nctp

!--- set initial values for time varying properties
    Model%ipt              = 1
    Model%lprnt            = .false.
    Model%lsswr            = .false.
    Model%lslwr            = .false.
    Model%solhr            = -9999.
    Model%solcon           = -9999.
    Model%slag             = -9999.
    Model%sdec             = -9999.
    Model%cdec             = -9999.
    Model%clstp            = -9999
    rinc(1:5)              = 0
    call w3difdat(jdat,idat,4,rinc)
    Model%phour            = rinc(4)/con_hr
    Model%fhour            = (rinc(4) + Model%dtp)/con_hr
    Model%zhour            = mod(Model%phour,Model%fhzero)
    Model%kdt              = nint(Model%fhour*con_hr/Model%dtp)
    Model%first_time_step  = .true.
    Model%restart          = restart
    Model%lsm_cold_start   = .not. restart
    Model%hydrostatic      = hydrostatic
    if (Model%me == Model%master) then
      print *,'in atm phys init, phour=',Model%phour,'fhour=',Model%fhour,'zhour=',Model%zhour,'kdt=',Model%kdt
    endif
    

    if(Model%hydrostatic .and. Model%lightning_threat) then
      write(0,*) 'Turning off lightning threat index for hydrostatic run.'
      Model%lightning_threat = .false.
      lightning_threat = .false.
    endif

    Model%jdat(1:8)        = jdat(1:8)
    allocate(Model%si(Model%levs+1))
    !--- Define sigma level for radiation initialization
    !--- The formula converting hybrid sigma pressure coefficients to sigma coefficients follows Eckermann (2009, MWR)
    !--- ps is replaced with p0. The value of p0 uses that in http://www.emc.ncep.noaa.gov/officenotes/newernotes/on461.pdf
    !--- ak/bk have been flipped from their original FV3 orientation and are defined sfc -> toa
    Model%si(1:Model%levs+1) = (ak(1:Model%levs+1) + bk(1:Model%levs+1) * con_p0 - ak(Model%levs+1)) / (con_p0 - ak(Model%levs+1))
    Model%sec              = 0
    Model%yearlen          = 365
    Model%julian           = -9999.
    !--- Set vertical flag used by radiation schemes
    Model%top_at_1         = .false.
    if (Model%do_RRTMGP) then
       if (Model%top_at_1) then
          Model%iSFC = Model%levs
          Model%iTOA = 1
       else
          Model%iSFC = 1
          Model%iTOA = Model%levs
       endif
    endif

!--- BEGIN CODE FROM GFS_PHYSICS_INITIALIZE
!--- define physcons module variables
    tem          = con_rerth*con_rerth*(con_pi+con_pi)*con_pi
    Model%dxmax  = log(tem/(max_lon*max_lat))
    Model%dxmin  = log(tem/(min_lon*min_lat))
    Model%dxinv  = 1.0d0 / (Model%dxmax-Model%dxmin)
    Model%rhcmax = rhcmax
    Model%huge   = huge
    if (Model%me == Model%master) write(*,*)' dxmax=',Model%dxmax,' dxmin=',Model%dxmin,' dxinv=',Model%dxinv, &
       'max_lon=',max_lon,' max_lat=',max_lat,' min_lon=',min_lon,' min_lat=',min_lat,       &
       ' rhc_max=',Model%rhcmax,' huge=',huge

!--- set nrcm
    if (Model%ras) then
      Model%nrcm = min(nrcmax, Model%levs-1) * (Model%dtp/1200.d0) + 0.10001d0
    else
      Model%nrcm = 2
    endif

!--- cal_pre
    if (Model%cal_pre) then
      Model%random_clds = .true.
    endif
!--- END CODE FROM GFS_PHYSICS_INITIALIZE


!--- BEGIN CODE FROM COMPNS_PHYSICS
!--- shoc scheme
    if (do_shoc) then
      if (Model%imp_physics == Model%imp_physics_thompson) then
        print *,'SHOC is not currently compatible with Thompson MP -- shutting down'
        stop
      endif
      Model%nshoc_3d   = 3
      Model%nshoc_2d   = 0
      Model%shal_cnv   = .false.
      Model%imfshalcnv = -1
      Model%hybedmf    = .false.
      Model%satmedmf   = .false.
      if (Model%me == Model%master) print *,' Simplified Higher Order Closure Model used for', &
                                            ' Boundary layer and Shallow Convection',          &
                                            ' nshoc_3d=',Model%nshoc_3d,                       &
                                            ' nshoc_2d=',Model%nshoc_2d,                       &
                                            ' ntke=',Model%ntke,' shoc_parm=',shoc_parm
    endif

    !--- mynn-edmf scheme
    if (Model%do_mynnedmf) then
      if (Model%do_shoc .or. Model%hybedmf .or. Model%satmedmf) then
          print *,' Logic error: MYNN EDMF cannot be run with SHOC, HEDMF or SATMEDMF'
          stop
      end if
!      Model%shal_cnv   = .false.
!      Model%imfshalcnv = -1
      ! DH* substitute for MYNN namelist section
      Model%icloud_bl         = 1
      !Model%bl_mynn_tkeadvect = .true.
      Model%bl_mynn_edmf      = 1
      !Model%bl_mynn_edmf_mom  = 1
      ! *DH
      if (Model%me == Model%master) print *,' MYNN-EDMF scheme is used for both',                &
                                            ' boundary layer turbulence and shallow convection', &
                                            ' bl_mynn_cloudpdf=',Model%bl_mynn_cloudpdf,         &
                                            ' bl_mynn_mixlength=',Model%bl_mynn_mixlength,       &
                                            ' bl_mynn_edmf=',Model%bl_mynn_edmf,                 &
                                            ' bl_mynn_output=',Model%bl_mynn_output,             &
                                            ' bl_mynn_closure=',Model%bl_mynn_closure
    endif

    !--- mynn surface layer scheme
    if (Model%do_mynnsfclay) then
      if (Model%me == Model%master) print *,' MYNN surface layer scheme is used:',               &
                                            ' isftcflx=',Model%isftcflx,                         &
                                            ' iz0tlnd=',Model%iz0tlnd,                           &
                                            ' sfclay_compute_diag=',Model%sfclay_compute_diag,   &
                                            ' sfclay_compute_flux=',Model%sfclay_compute_flux
    end if

!--- set number of cloud types
    if (Model%cscnv) then
      Model%nctp = nint(Model%cs_parm(5))
      Model%nctp = max(Model%nctp,10)
      if (Model%cs_parm(7) < 0.0) Model%cs_parm(7) = Model%dtp
      Model%do_awdd  = Model%do_aw .and. Model%cs_parm(6) > 0.0
!     Model%flx_form = Model%do_aw .and. Model%cs_parm(8) > 0.0
      Model%flx_form = Model%cs_parm(8) > 0.0
    endif
!   Model%nctp = max(Model%nctp,1)

!--- output information about the run
    if (Model%me == Model%master) then
      if (Model%lsm == 1) then
        print *,' NOAH Land Surface Model used'

      elseif (Model%lsm == 0) then
        print *,' OSU no longer supported - job aborted'
        stop
      elseif (Model%lsm == Model%lsm_noahmp) then
        if (Model%ivegsrc /= 1) then
          print *,'Vegetation type must be IGBP if Noah MP is used'
          stop
        elseif (Model%isot /= 1) then
          print *,'Soil type must be STATSGO if Noah MP is used'
          stop
        endif
        print *, 'New Noah MP Land Surface Model will be used'
        print *, 'The Physics options are'

        print *,'iopt_dveg  =  ', Model%iopt_dveg
        print *,'iopt_crs   =  ', Model%iopt_crs
        print *,'iopt_btr   =  ', Model%iopt_btr
        print *,'iopt_run   =  ', Model%iopt_run
        print *,'iopt_sfc   =  ', Model%iopt_sfc
        print *,'iopt_frz   =  ', Model%iopt_frz
        print *,'iopt_inf   =  ', Model%iopt_inf
        print *,'iopt_rad   =  ', Model%iopt_rad
        print *,'iopt_alb   =  ', Model%iopt_alb
        print *,'iopt_snf   =  ', Model%iopt_snf
        print *,'iopt_tbot   =  ',Model%iopt_tbot
        print *,'iopt_stc   =  ', Model%iopt_stc
        print *,'iopt_trs   =  ', Model%iopt_trs
        print *,'iopt_diag  =  ', Model%iopt_diag
      elseif (Model%lsm == Model%lsm_ruc) then
        print *,' RUC Land Surface Model used'
        print *, 'The Physics options are'
        print *,' mosaic_lu   =  ',mosaic_lu
        print *,' mosaic_soil =  ',mosaic_soil
        print *,' isncond_opt =  ',isncond_opt
        print *,' isncovr_opt =  ',isncovr_opt
        print *,' add_fire_heat_flux = ',add_fire_heat_flux
      else
        print *,' Unsupported LSM type - job aborted - lsm=',Model%lsm
        stop
      endif

!      if (Model%lsm == Model%lsm_noahmp .and. Model%iopt_snf == 4) then
!        if (Model%imp_physics /= Model%imp_physics_gfdl) stop 'iopt_snf == 4 must use GFDL MP'
!      endif

      print *,' nst_anl=',Model%nst_anl,' use_ufo=',Model%use_ufo,' frac_grid=',Model%frac_grid,&
              ' ignore_lake=',ignore_lake,' frac_ice=',Model%frac_ice
      print *,' min_lakeice=',Model%min_lakeice,' min_seaice=',Model%min_seaice,                &
              'min_lake_height=',Model%min_lake_height

      print *, 'lake model parameters'
      print *, ' lake master flag lkm   : ', Model%lkm
      if(Model%lkm>0) then
        print *, ' lake model selection   : ', Model%iopt_lake
        if(Model%iopt_lake==Model%iopt_lake_clm) then
          print *,'  CLM Lake model configuration'
          print *,'   use_lake2m             = ',Model%use_lake2m
          print *,'   clm_lake_use_lakedepth = ',Model%clm_lake_use_lakedepth
          print *,'   clm_lake_depth_default = ',Model%clm_lake_depth_default
          print *,'   clm_lake_debug         = ',Model%clm_lake_debug
          print *,'   clm_debug_print        = ',Model%clm_debug_print
          print *,'   nlevlake_clm_lake      = ',Model%nlevlake_clm_lake
          print *,'   nlevsoil_clm_lake      = ',Model%nlevsoil_clm_lake
          print *,'   nlevsnow_clm_lake      = ',Model%nlevsnow_clm_lake
          print *,'   nlevsnowsoil_clm_lake  = ',Model%nlevsnowsoil_clm_lake
          print *,'   nlevsnowsoil1_clm_lake = ',Model%nlevsnowsoil1_clm_lake
        endif
      endif

      if (Model%nstf_name(1) > 0 ) then
        print *,' NSSTM is active '
        print *,' nstf_name(1)=',Model%nstf_name(1)
        print *,' nstf_name(2)=',Model%nstf_name(2)
        print *,' nstf_name(3)=',Model%nstf_name(3)
        print *,' nstf_name(4)=',Model%nstf_name(4)
        print *,' nstf_name(5)=',Model%nstf_name(5)
      endif
      if (Model%do_deep) then
        ! Consistency check for NTDK convection: deep and shallow convection are bundled
        ! and cannot be combined with any other deep or shallow convection scheme
        if ( (Model%imfdeepcnv == Model%imfdeepcnv_ntiedtke .or. Model%imfshalcnv == Model%imfshalcnv_ntiedtke) .and. &
            .not. (Model%imfdeepcnv == Model%imfdeepcnv_ntiedtke .and. Model%imfshalcnv == Model%imfshalcnv_ntiedtke) ) then
            write(0,*) "Logic error: if NTDK deep convection is used, must also use NTDK shallow convection (and vice versa)"
            stop
        end if

        if (.not. Model%cscnv) then
          if (Model%ras) then
            print *,' RAS Convection scheme used with ccwf=',Model%ccwf
            Model%imfdeepcnv = -1
          else
            if (Model%imfdeepcnv == 0) then
               print *,' old SAS Convection scheme before July 2010 used'
            elseif(Model%imfdeepcnv == Model%imfdeepcnv_sas) then
               print *,' July 2010 version of SAS conv scheme used'
            elseif(Model%imfdeepcnv == Model%imfdeepcnv_samf) then
               print *,' scale & aerosol-aware mass-flux deep conv scheme'
            elseif(Model%imfdeepcnv == Model%imfdeepcnv_gf) then
               print *,' Grell-Freitas scale & aerosol-aware mass-flux deep conv scheme'
            elseif(Model%imfdeepcnv == Model%imfdeepcnv_ntiedtke) then
               print *,' New Tiedtke cumulus scheme'
            elseif(Model%imfdeepcnv == Model%imfdeepcnv_c3) then
               print *,' New unified cumulus convection scheme'
            endif
          endif
        else
          if (Model%do_aw) then
            print *,'Chikira-Sugiyama convection scheme with Arakawa-Wu'&
     &,                ' unified parameterization used'
          else
              print *,'Chikira-Sugiyama convection scheme used'
          endif
          print *,' cs_parm=',Model%cs_parm,' nctp=',Model%nctp
        endif
      else
        print*, ' Deep convection scheme disabled'
      endif
      if (Model%satmedmf) then
        if (Model%isatmedmf == Model%isatmedmf_vdif) then
          print *,' initial version (Nov 2018) of sale-aware TKE-based moist EDMF scheme used'
        elseif(Model%isatmedmf == Model%isatmedmf_vdifq) then
          print *,' update version (May 2019) of sale-aware TKE-based moist EDMF scheme used'
        endif
      elseif (Model%hybedmf) then
        print *,' scale-aware hybrid edmf PBL scheme used'
      elseif (Model%old_monin) then
        print *,' old (old_monin) PBL scheme used'
      elseif (Model%do_mynnedmf) then
        print *,' MYNN PBL scheme used'
      elseif (Model%do_myjpbl)then
        print *,' MYJ PBL scheme used'
      endif
      if (.not. Model%shal_cnv) then
        Model%imfshalcnv = -1
        print *,' No shallow convection used'
      else
        if (Model%imfshalcnv == 0) then
          print *,' modified Tiedtke eddy-diffusion shallow conv scheme used'
        elseif (Model%imfshalcnv == Model%imfshalcnv_sas) then
          print *,' July 2010 version of mass-flux shallow conv scheme used'
        elseif (Model%imfshalcnv == Model%imfshalcnv_samf) then
          print *,' scale- & aerosol-aware mass-flux shallow conv scheme (2017)'
        elseif (Model%imfshalcnv == Model%imfshalcnv_gf) then
          print *,' Grell-Freitas scale- & aerosol-aware mass-flux shallow conv scheme (2013)'
        elseif (Model%imfshalcnv == Model%imfshalcnv_ntiedtke) then
          print *,' New Tiedtke cumulus scheme'
        elseif (Model%imfshalcnv == Model%imfshalcnv_c3) then
          print *,' New unified cumulus scheme'
        else
          print *,' unknown mass-flux scheme in use - defaulting to no shallow convection'
          Model%imfshalcnv = -1
        endif
      endif
      if (Model%do_gwd) then
        if (Model%do_ugwp) then
          print *,' Unified gravity wave drag parameterization used'
        elseif (Model%gwd_opt == 2) then
          print *,'GSL unified oragraphic gravity wave drag parameterization used'
        else
          print *,' Original mountain blocking and oragraphic gravity wave drag parameterization used'
          if (cdmbgwd(3) > 0.0) print *,' non-statioary gravity wave drag parameterization used'
        endif
          print *,' do_gwd=',Model%do_gwd
      endif
      if (Model%do_cnvgwd) then
        print *,' Convective GWD parameterization used, do_cnvgwd=',Model%do_cnvgwd
      endif
      if (Model%lcrick) print *,' CRICK-Proof cloud water used in radiation '
      if (Model%lcnorm) print *,' Cloud condensate normalized by cloud cover for radiation'
      if (Model%iovr == Model%iovr_rand) then
         print *,' random cloud overlap for Radiation IOVR=',            Model%iovr
      elseif (Model%iovr == Model%iovr_dcorr) then
         print *,' exponential-decorr cloud overlap for Radiation IOVR=',Model%iovr
      elseif (Model%iovr == Model%iovr_exp) then
         print *,' exponential cloud overlap for Radiation IOVR=',       Model%iovr
      elseif (Model%iovr == Model%iovr_exprand) then
         print *,' exponential-random cloud overlap for Radiation IOVR=',Model%iovr
      else
         print *,' max-random cloud overlap for Radiation IOVR=',        Model%iovr
      endif
      if (Model%isubc_sw == 0) then
        print *,' no sub-grid cloud for Shortwave ISUBC_SW=',Model%isubc_sw
      else
        print *,' sub-grid cloud for Shortwave ISUBC_SW=',Model%isubc_sw
      endif
      if (Model%isubc_lw == 0) then
        print *,' no sub-grid cloud for Longwave ISUBC_LW=',Model%isubc_lw
      else
        print *,' sub-grid cloud for Longwave ISUBC_LW=',Model%isubc_lw
      endif
    endif

    ! get_alpha routines for exponential and exponential-random overlap need this(!?!)
    if (Model%iovr == Model%iovr_exprand .or. Model%iovr == Model%iovr_exp) then
       Model%yearlen = 365
    endif

!--- set up cloud schemes and tracer elements
    Model%nleffr   = -999
    Model%nieffr   = -999
    Model%nreffr   = -999
    Model%nseffr   = -999
    Model%ngeffr   = -999
    Model%nT2delt  = -999
    Model%nTdelt   = -999
    Model%nqv2delt = -999
    Model%nqvdelt  = -999
    Model%nps2delt = -999
    Model%npsdelt  = -999
    Model%ncnd     = nwat - 1                   ! ncnd is the number of cloud condensate types
    if (Model%imp_physics == Model%imp_physics_zhao_carr) then
      Model%npdf3d   = 0
      Model%num_p3d  = 4
      Model%num_p2d  = 3
      Model%shcnvcw  = .false.
      Model%nT2delt  = 1
      Model%nqv2delt = 2
      Model%nTdelt   = 3
      Model%nqvdelt  = 4
      Model%nps2delt = 1
      Model%npsdelt  = 2
      if (nwat /= 2) then
        print *,' Zhao-Carr MP requires nwat to be set to 2 - job aborted'
        stop
      end if
      if (Model%me == Model%master) print *,' Using Zhao/Carr/Sundqvist Microphysics'

    elseif (Model%imp_physics == Model%imp_physics_zhao_carr_pdf) then !Zhao Microphysics with PDF cloud
      Model%npdf3d  = 3
      Model%num_p3d = 4
      Model%num_p2d = 3
      if (Model%me == Model%master) print *,'Using Zhao/Carr/Sundqvist Microphysics with PDF Cloud'

    else if (Model%imp_physics == Model%imp_physics_fer_hires) then     ! Ferrier-Aligo scheme
      Model%npdf3d  = 0
      Model%num_p3d = 3
      Model%num_p2d = 1
      Model%pdfcld  = .false.
      Model%shcnvcw = .false.
      Model%nleffr  = 1
      Model%nieffr  = 2
      Model%nseffr  = 3
      if (nwat /= 4) then
        print *,' Ferrier-Aligo MP requires nwat to be set to 4 - job aborted'
        stop
      end if
      if (Model%me == Model%master) print *,' Using Ferrier-Aligo MP scheme', &
                                          ' microphysics', &
                                          ' lradar =',Model%lradar

    elseif (Model%imp_physics == Model%imp_physics_wsm6) then !WSM6 microphysics
      print *,' Error, WSM6 no longer supported - job aborted'
      stop
      !Model%npdf3d  = 0
      !Model%num_p3d = 3
      !Model%num_p2d = 1
      !Model%pdfcld  = .false.
      !Model%shcnvcw = .false.
      !Model%nleffr  = 1
      !Model%nieffr  = 2
      !Model%nseffr  = 3
      !if (Model%me == Model%master) print *,' Using wsm6 microphysics'

    elseif (Model%imp_physics == Model%imp_physics_nssl) then !NSSL microphysics
      Model%npdf3d  = 0
      Model%num_p3d = 4 ! for size of phy3d
      Model%num_p2d = 1
      Model%pdfcld  = .false.
      Model%shcnvcw = .false.
      IF ( Model%nssl_hail_on ) THEN
        i = 1
      ELSE
        i = 0
      ENDIF
      if ( nwat /= 6+i ) then
        print *,' NSSL MP requires nwat to be set to ', 6+i,' - job aborted, nssl_hail_on = ',nssl_hail_on
        stop
      end if
      Model%nleffr = 1
      Model%nieffr = 2
      Model%nseffr = 3
      Model%nreffr = 4
      Model%lradar = .true.
      if (.not. Model%effr_in) then
        print *,' NSSL MP requires effr_in to be set to .true., changing value from false to true'
        Model%effr_in = .true.
        effr_in = .true.
      ENDIF
      if (Model%me == Model%master) print *,' Using NSSL double moment microphysics', &
                                          ' nssl_ccn_on =',Model%nssl_ccn_on,   &
                                          ' nssl_invertccn =',Model%nssl_invertccn,   &
                                          ' lradar =',Model%lradar, &
                                          ' num_p3d =',Model%num_p3d, &
                                          ' num_p2d =',Model%num_p2d


    elseif (Model%imp_physics == Model%imp_physics_thompson) then !Thompson microphysics
      Model%npdf3d  = 0
      Model%num_p3d = 3
      Model%num_p2d = 1
      Model%pdfcld  = .false.
      Model%shcnvcw = .false.
      Model%nleffr  = 1
      Model%nieffr  = 2
      Model%nseffr  = 3
      if (nwat /= 6) then
        print *,' Thompson MP requires nwat to be set to 6 - job aborted'
        stop
      end if
      if (.not. Model%effr_in) then
        print *,' Thompson MP requires effr_in to be set to .true. - job aborted'
        stop
      end if
      if (Model%me == Model%master) print *,' Using Thompson double moment microphysics', &
                                          ' ltaerosol = ',Model%ltaerosol, &
                                          ' mraerosol = ',Model%mraerosol, &
                                          ' ttendlim =',Model%ttendlim, &
                                          ' ext_diag_thompson =',Model%ext_diag_thompson, &
                                          ' dt_inner =',Model%dt_inner, &
                                          ' sedi_semi=',Model%sedi_semi, &
                                          ' decfl=',decfl, &
                                          ' effr_in =',Model%effr_in, &
                                          ' lradar =',Model%lradar, &
                                          ' nsfullradar_diag =',Model%nsfullradar_diag, &
                                          ' num_p3d =',Model%num_p3d, &
                                          ' num_p2d =',Model%num_p2d

    else if (Model%imp_physics == Model%imp_physics_mg) then        ! Morrison-Gettelman Microphysics
      Model%npdf3d  = 0
      Model%num_p3d = 5
      Model%num_p2d = 1
      Model%pdfcld  = .false.
      Model%shcnvcw = .false.
      Model%nleffr  = 2
      Model%nieffr  = 3
      Model%nreffr  = 4
      Model%nseffr  = 5
      if (Model%mg_do_graupel .or. Model%mg_do_hail) then
        Model%num_p3d = 6
        Model%ngeffr  = 6
      endif
      if (nwat /= 6 .and. Model%fprcp >= 2) then
        print *,' Morrison-Gettelman MP requires nwat to be set to 6 - job aborted'
        stop
      end if
      if (Model%me == Model%master)                                                                 &
         print *,' Using Morrison-Gettelman double moment microphysics',                            &
                 ' iaerclm=',         Model%iaerclm,         ' iccn=',          Model%iccn,         &
                 ' mg_dcs=',          Model%mg_dcs,          ' mg_qcvar=',      Model%mg_qcvar,     &
                 ' mg_ts_auto_ice=',  Model%mg_ts_auto_ice,  ' pdfflag=',       Model%pdfflag,      &
                 ' mg_do_graupel=',   Model%mg_do_graupel,   ' mg_do_hail=',    Model%mg_do_hail,   &
                 ' mg_nccons=',       Model%mg_nccons,       ' mg_nicon=',      Model%mg_nicons,    &
                 ' mg_ngcons=',       Model%mg_ngcons ,      ' mg_ncnst=',      Model%mg_ncnst,     &
                 ' mg_ninst=',        Model%mg_ninst ,       ' mg_ngnst=',      Model%mg_ngnst,     &
                 ' sed_supersat=',    Model%sed_supersat ,   ' do_sb_physics=', Model%do_sb_physics,&
                 ' microp_uniform=',  Model%microp_uniform,  ' do_cldice=',     Model%do_cldice,    &
                 ' hetfrz_classnuc=', Model%hetfrz_classnuc, ' ncnd=',          Model%ncnd,         &
                 ' mg_alf=',          Model%mg_alf,          ' mg_qcmin=',      Model%mg_qcmin,     &
                 ' mg_do_ice_gmao=',  Model%mg_do_ice_gmao,  ' mg_do_liq_liu=', Model%mg_do_liq_liu

    elseif (Model%imp_physics == Model%imp_physics_gfdl) then !GFDL microphysics
      Model%npdf3d  = 0
      if(Model%effr_in) then
        Model%num_p3d = 5
        Model%nleffr  = 1
        Model%nieffr  = 2
        Model%nreffr  = 3
        Model%nseffr  = 4
        Model%ngeffr  = 5
      else
        Model%num_p3d = 1
        ! Effective radii not used, point to valid index in dummy phy_f3d array
        Model%nleffr = 1
        Model%nieffr = 1
        Model%nreffr = 1
        Model%nseffr = 1
        Model%ngeffr = 1
      end if
      Model%num_p2d = 1
      Model%pdfcld  = .false.
      Model%shcnvcw = .false.
      if (nwat /= 6) then
        print *,' GFDL MP requires nwat to be set to 6 - job aborted'
        stop
      end if
      if (Model%me == Model%master) print *,' avg_max_length=',Model%avg_max_length
      if (Model%me == Model%master) print *,' Using GFDL Cloud Microphysics'

    else
      if (Model%me == Model%master) print *,'Wrong imp_physics value. Job abort.'
      stop
    endif

    if(Model%ras     .or. Model%cscnv)  Model%cnvcld = .false.
    if(Model%do_shoc .or. Model%pdfcld .or. Model%do_mynnedmf .or. Model%imfdeepcnv == Model%imfdeepcnv_gf .or. Model%imfdeepcnv == Model%imfdeepcnv_c3) Model%cnvcld = .false.
    if(Model%cnvcld) Model%ncnvcld3d = 1

!--- get cnvwind index in phy_f2d; last entry in phy_f2d array
    Model%ncnvwind = Model%num_p2d

!--- get cnvw and cnvc indices in phy_f3d
    Model%ncnvw = -999
    Model%ncnvc = -999
    if ((Model%npdf3d == 3) .and. (Model%num_p3d == 4)) then
      Model%ncnvw = Model%num_p3d + 2
      Model%ncnvc = Model%ncnvw + 1
    elseif ((Model%npdf3d == 0) .and. (Model%ncnvcld3d == 1)) then
      Model%ncnvw = Model%num_p3d + 1
    endif

!--- derived totals for phy_f*d
    Model%ntot2d = Model%num_p2d + Model%nshoc_2d
    Model%ntot3d = Model%num_p3d + Model%nshoc_3d + Model%npdf3d + Model%ncnvcld3d
!
!   Unified cloud for SHOC and/or MG3
    Model%uni_cld = .false.
    Model%indcld  = -1
    if (Model%imp_physics == Model%imp_physics_mg) then
      Model%uni_cld = .true.
      Model%indcld  = 1
    elseif (Model%shoc_cld) then
      Model%uni_cld = .true.
      Model%indcld  = Model%ntot3d - 2
    endif

    if (Model%do_shoc) then
      Model%nkbfshoc = Model%ntot3d   !< the index of upward kinematic buoyancy flux from SHOC in phy_f3d
      Model%nahdshoc = Model%ntot3d-1 !< the index of diffusivity for heat from from SHOC in phy_f3d
      Model%nscfshoc = Model%ntot3d-2 !< the index of subgrid-scale cloud fraction from from SHOC in phy_f3d
    else
      Model%nkbfshoc = -999
      Model%nahdshoc = -999
      Model%nscfshoc = -999
    endif

    if (me == Model%master)                                                     &
      write(*,*) ' num_p3d=',   Model%num_p3d,   ' num_p2d=',  Model%num_p2d,   &
                 ' crtrh=',     Model%crtrh,     ' npdf3d=',   Model%npdf3d,    &
                 ' pdfcld=',    Model%pdfcld,    ' shcnvcw=',  Model%shcnvcw,   &
                 ' cnvcld=',    Model%cnvcld,    ' ncnvcld3d=',Model%ncnvcld3d, &
                 ' do_shoc=',   Model%do_shoc,   ' nshoc3d=',  Model%nshoc_3d,  &
                 ' nshoc_2d=',  Model%nshoc_2d,  ' shoc_cld=', Model%shoc_cld,  &
                 ' nkbfshoc=',  Model%nkbfshoc,  ' nahdshoc=', Model%nahdshoc,  &
                 ' nscfshoc=',  Model%nscfshoc,  ' xr_cnvcld=',Model%xr_cnvcld, &
                 ' uni_cld=',   Model%uni_cld,                                  &
                 ' ntot3d=',    Model%ntot3d,    ' ntot2d=',   Model%ntot2d,    &
                 ' shocaftcnv=',Model%shocaftcnv,' indcld=',   Model%indcld,    &
                 ' shoc_parm=', Model%shoc_parm,                                &
                 ' ncnvw=',     Model%ncnvw,     ' ncnvc=',     Model%ncnvc

!--- END CODE FROM COMPNS_PHYSICS


!--- BEGIN CODE FROM GLOOPR
!--- set up parameters for Xu & Randall's cloudiness computation (Radiation)

    Model%lmfshal  = (Model%shal_cnv .and. Model%imfshalcnv > 0)
    Model%lmfdeep2 = (Model%imfdeepcnv == Model%imfdeepcnv_samf         &
                      .or. Model%imfdeepcnv == Model%imfdeepcnv_gf      &
                      .or. Model%imfdeepcnv == Model%imfdeepcnv_ntiedtke &
                      .or. Model%imfdeepcnv == Model%imfdeepcnv_c3)
!--- END CODE FROM GLOOPR

!--- BEGIN CODE FROM GLOOPB
!--- set up random number seed needed for RAS and old SAS and when cal_pre=.true.
!    Model%imfdeepcnv < 0 when Model%ras = .true.

    if (Model%imfdeepcnv <= 0 .or. Model%cal_pre ) then
      if (Model%random_clds) then
        seed0 = Model%idate(1) + Model%idate(2) + Model%idate(3) + Model%idate(4)
        call random_setseed(seed0)
        call random_number(wrk)
        Model%seed0 = seed0 + nint(wrk(1)*1000.0d0)
      endif
    endif
!--- END CODE FROM GLOOPB

    call Model%print ()

  end subroutine control_initialize

  subroutine control_initialize_radar_tten(Model, radar_tten_limits)
    implicit none

    ! Helper subroutine for initializing variables for radar-derived
    ! temperature tendency or convection suppression.

    class(GFS_control_type) :: Model
    real(kind_phys) :: radar_tten_limits(2)
    integer :: i

    Model%num_dfi_radar    = 0
    do i=1,dfi_radar_max_intervals
       if(Model%fh_dfi_radar(i)>-1e10 .and. Model%fh_dfi_radar(i+1)>-1e10) then
          Model%num_dfi_radar = Model%num_dfi_radar+1
          Model%ix_dfi_radar(i) = Model%num_dfi_radar
       else
          Model%ix_dfi_radar(i) = -1
       endif
    enddo

    if(Model%num_dfi_radar>0) then
       if(radar_tten_limits(1)==limit_unspecified) then
          if(radar_tten_limits(2)==limit_unspecified) then
             radar_tten_limits(1) = -19
             radar_tten_limits(2) = 19
             if(Model%me==Model%master) then
                write(0,*) 'Warning: using internal defaults for radar_tten_limits. If the oceans boil, try different values.'
                write(0,'(A,F12.4,A)') 'radar_tten_limits(1) = ',radar_tten_limits(1),' <-- lower limit'
                write(0,'(A,F12.4,A)') 'radar_tten_limits(2) = ',radar_tten_limits(2),' <-- upper limit'
             endif
          else
             radar_tten_limits(1) = -abs(radar_tten_limits(2))
             radar_tten_limits(2) = abs(radar_tten_limits(2))
          endif
       else if(radar_tten_limits(2)==limit_unspecified) then
           radar_tten_limits(1) = -abs(radar_tten_limits(1))
           radar_tten_limits(2) = abs(radar_tten_limits(1))
       else if(radar_tten_limits(1)>radar_tten_limits(2)) then
          if(Model%me==Model%master) then
             write(0,*) 'Error: radar_tten_limits lower limit is higher than upper!'
             write(0,'(A,F12.4,A)') 'radar_tten_limits(1) = ',radar_tten_limits(1),' <-- lower limit'
             write(0,'(A,F12.4,A)') 'radar_tten_limits(2) = ',radar_tten_limits(2),' <-- upper limit'
             write(0,*) "If you do not want me to apply the prescribed tendencies, just say so! Remove fh_dfi_radar from your namelist."
             stop
          endif
       else
          !o! Rejoice !o! Radar_tten_limits had lower and upper bounds.
       endif
       Model%radar_tten_limits = radar_tten_limits

       if(Model%do_cap_suppress) then
         if(Model%me==Model%master .and. Model%imfdeepcnv>=0) then
           if(Model%imfdeepcnv/=3) then
             write(0,*) 'Warning: untested configuration in use! Radar-derived convection suppression is only supported for the GF deep scheme. That feature will be inactive, but microphysics tendencies will still be enabled. This combination is untested. Beware!'
           else
             write(0,*) 'Warning: experimental configuration in use! Radar-derived convection suppression is experimental (GF deep scheme with fh_dfi_radar).'
           endif
         endif
       endif
    endif

  end subroutine control_initialize_radar_tten

!---------------------------
! GFS_control%init_chemistry
!---------------------------
  subroutine control_chemistry_initialize(Model, tracer_types)

    !--- Identify number and starting/ending indices of both
    !--- prognostic and diagnostic chemistry tracers.
    !--- Each tracer set is assumed to be contiguous.

    use parse_tracers, only: NO_TRACER

    !--- interface variables
    class(GFS_control_type) :: Model
    integer,     intent(in) :: tracer_types(:)

    !--- local variables
    integer :: n

    !--- begin
    Model%nchem = 0
    Model%ndvel = 0
    Model%ntchm = 0
    Model%ntchs = NO_TRACER
    Model%ntche = NO_TRACER
    Model%ndchm = 0
    Model%ndchs = NO_TRACER
    Model%ndche = NO_TRACER

    if (Model%rrfs_sd) then
      Model%nchem = 3
      Model%ndvel = 3
    endif

    do n = 1, size(tracer_types)
      select case (tracer_types(n))
        case (1)
          ! -- prognostic chemistry tracers
          Model%ntchm = Model%ntchm + 1
          if (Model%ntchm == 1) Model%ntchs = n
        case (2)
          ! -- diagnostic chemistry tracers
          Model%ndchm = Model%ndchm + 1
          if (Model%ndchm == 1) Model%ndchs = n
        case default
          ! -- generic tracers
      end select
    end do

    if (Model%ntchm > 0) Model%ntche = Model%ntchs + Model%ntchm - 1
    if (Model%ndchm > 0) Model%ndche = Model%ndchs + Model%ndchm - 1

  end subroutine control_chemistry_initialize


!----------------------------
! GFS_control%init_scavenging
!----------------------------
  subroutine control_scavenging_initialize(Model, fscav)

    use parse_tracers, only: get_tracer_index

    !--- interface variables
    class(GFS_control_type)      :: Model
    character(len=*), intent(in) :: fscav(:)

    !--- local variables
    integer              :: i, ios, j, n
    real(kind=kind_phys) :: tem

    !--- begin
    allocate(Model%fscav(Model%ntchm))

    if (Model%ntchm > 0) then
      !--- set default as no scavenging
      Model%fscav = zero
      ! -- read factors from namelist
      ! -- set default first, if available
      do i = 1, size(fscav)
        j = index(fscav(i),":")
        if (j > 1) then
          read(fscav(i)(j+1:), *, iostat=ios) tem
          if (ios /= 0) cycle
          if (adjustl(fscav(i)(:j-1)) == "*") then
            Model%fscav = tem
            exit
          endif
        endif
      enddo
      ! -- then read factors for each tracer
      do i = 1, size(fscav)
        j = index(fscav(i),":")
        if (j > 1) then
          read(fscav(i)(j+1:), *, iostat=ios) tem
          if (ios /= 0) cycle
          n = get_tracer_index(Model%tracer_names, adjustl(fscav(i)(:j-1)), Model%me, Model%master, Model%debug) &
              - Model%ntchs + 1
          if (n > 0) Model%fscav(n) = tem
        endif
      enddo
    endif

  end subroutine control_scavenging_initialize


!------------------
! GFS_control%print
!------------------
  subroutine control_print(Model)

    implicit none

!--- interface variables
    class(GFS_control_type) :: Model

!--- local variables
    integer :: i

    if (Model%me == Model%master) then
      print *, ' '
      print *, 'basic control parameters'
      print *, ' me                : ', Model%me
      print *, ' master            : ', Model%master
      print *, ' communicator      : ', Model%communicator%mpi_val
      print *, ' nlunit            : ', Model%nlunit
      print *, ' fn_nml            : ', trim(Model%fn_nml)
      print *, ' fhzero            : ', Model%fhzero
      print *, ' fhzero_array      : ', Model%fhzero_array
      print *, ' fhzero_fhour      : ', Model%fhzero_fhour
      print *, ' ldiag3d           : ', Model%ldiag3d
      print *, ' qdiag3d           : ', Model%qdiag3d
      print *, ' lssav             : ', Model%lssav
      print *, ' naux2d            : ', Model%naux2d
      print *, ' naux3d            : ', Model%naux3d
      if (Model%naux2d>0) then
        print *, ' aux2d_time_avg    : ', Model%aux2d_time_avg
      endif
      if (Model%naux3d>0) then
        print *, ' aux3d_time_avg    : ', Model%aux3d_time_avg
      endif
      print *, ' fhcyc             : ', Model%fhcyc
      print *, ' thermodyn_id      : ', Model%thermodyn_id
      print *, ' sfcpress_id       : ', Model%sfcpress_id
      print *, ' gen_coord_hybrid  : ', Model%gen_coord_hybrid
      print *, ' hydrostatic       : ', Model%hydrostatic
      print *, ' '
      print *, 'grid extent parameters'
      print *, ' isc               : ', Model%isc
      print *, ' jsc               : ', Model%jsc
      print *, ' nx                : ', Model%nx
      print *, ' ny                : ', Model%ny
      print *, ' levs              : ', Model%levs
      print *, ' cnx               : ', Model%cnx
      print *, ' cny               : ', Model%cny
      print *, ' lonr              : ', Model%lonr
      print *, ' latr              : ', Model%latr
      print *, ' blksz(1)          : ', Model%blksz(1)
      print *, ' blksz(nblks)      : ', Model%blksz(Model%nblks)
      print *, ' Model%ncols       : ', Model%ncols
      print *, ' '
      print *, 'coupling parameters'
      print *, ' cplflx            : ', Model%cplflx
      print *, ' cplice            : ', Model%cplice
      print *, ' cplocn2atm        : ', Model%cplocn2atm
      print *, ' cplwav            : ', Model%cplwav
      print *, ' cplwav2atm        : ', Model%cplwav2atm
      print *, ' cplaqm            : ', Model%cplaqm
      print *, ' cplchm            : ', Model%cplchm
      print *, ' cpllnd            : ', Model%cpllnd
      print *, ' cpllnd2atm        : ', Model%cpllnd2atm
      print *, ' rrfs_sd           : ', Model%rrfs_sd
      print *, ' cpl_fire          : ', Model%cpl_fire
      print *, ' use_cice_alb      : ', Model%use_cice_alb
      print *, ' cpl_imp_mrg       : ', Model%cpl_imp_mrg
      print *, ' cpl_imp_dbg       : ', Model%cpl_imp_dbg
      print *, ' use_med_flux      : ', Model%use_med_flux
      if(Model%imfdeepcnv == Model%imfdeepcnv_gf .or.Model%imfdeepcnv == Model%imfdeepcnv_c3) then
        print*,'ichoice_s          : ', Model%ichoice_s
        print*,'ichoicem           : ', Model%ichoicem
        print*,'ichoice            : ', Model%ichoice
      endif
      if(model%rrfs_sd) then
        print *, ' '
        print *, 'smoke parameters'
        print *, 'dust_drylimit_factor: ',Model%dust_drylimit_factor
        print *, 'dust_moist_correction: ',Model%dust_moist_correction
        print *, 'dust_moist_opt   : ',Model%dust_moist_opt
        print *, 'dust_alpha       : ',Model%dust_alpha
        print *, 'dust_gamma       : ',Model%dust_gamma
        print *, 'wetdep_ls_alpha  : ',Model%wetdep_ls_alpha
        print *, 'ebb_dcycle       : ',Model%ebb_dcycle
        print *, 'seas_opt         : ',Model%seas_opt
        print *, 'dust_opt         : ',Model%dust_opt
        print *, 'drydep_opt       : ',Model%drydep_opt
        print *, 'coarsepm_settling: ',Model%coarsepm_settling
        print *, 'plume_wind_eff   : ',Model%plume_wind_eff
        print *, 'extended_sd_diags: ',Model%extended_sd_diags
        print *, 'wetdep_ls_opt    : ',Model%wetdep_ls_opt
        print *, 'do_plumerise     : ',Model%do_plumerise
        print *, 'plumerisefire_frq: ',Model%plumerisefire_frq
        print *, 'addsmoke_flag    : ',Model%addsmoke_flag
        print *, 'smoke_forecast   : ',Model%smoke_forecast
        print *, 'aero_ind_fdb     : ',Model%aero_ind_fdb
        print *, 'aero_dir_fdb     : ',Model%aero_dir_fdb
        print *, 'rrfs_smoke_debug : ',Model%rrfs_smoke_debug
        print *, 'do_smoke_transport : ',Model%do_smoke_transport
        print *, 'mix_chem         : ',Model%mix_chem
        print *, 'enh_mix          : ',Model%enh_mix
        print *, 'smoke_dir_fdb_coef : ',Model%smoke_dir_fdb_coef
        print *, 'smoke_conv_wet_coef: ',Model%smoke_conv_wet_coef
      endif
      print *, ' '
      print *, ' lsidea            : ', Model%lsidea
      print *, ' '
      print *, 'calendars and time parameters and activation triggers'
      print *, ' dtp               : ', Model%dtp
      print *, ' dtf               : ', Model%dtf
      print *, ' nscyc             : ', Model%nscyc
      print *, ' nszero            : ', Model%nszero
      print *, ' idat              : ', Model%idat
      print *, ' idate             : ', Model%idate
      print *, ' '
      print *, 'radiation control parameters'
      print *, ' fhswr             : ', Model%fhswr
      print *, ' fhlwr             : ', Model%fhlwr
      print *, ' nsswr             : ', Model%nsswr
      print *, ' nslwr             : ', Model%nslwr
      print *, ' nhfrad            : ', Model%nhfrad
      print *, ' levr              : ', Model%levr
      print *, ' nfxr              : ', Model%nfxr
      print *, ' ntrcaer           : ', Model%ntrcaer
      print *, ' lmfshal           : ', Model%lmfshal
      print *, ' lmfdeep2          : ', Model%lmfdeep2
      print *, ' nrcm              : ', Model%nrcm
      print *, ' iflip             : ', Model%iflip
      print *, ' isol              : ', Model%isol
      print *, ' ico2              : ', Model%ico2
      print *, ' ialb              : ', Model%ialb
      print *, ' iems              : ', Model%iems
      print *, ' iaer              : ', Model%iaer
      print *, ' iaermdl           : ', Model%iaermdl
      print *, ' iaerflg           : ', Model%iaerflg
      print *, ' lalw1bd           : ', Model%lalw1bd
      print *, ' aeros_file        : ', Model%aeros_file
      print *, ' solar_file        : ', Model%solar_file
      print *, ' semis_file        : ', Model%semis_file
      print *, ' icliq_sw          : ', Model%icliq_sw
      print *, ' icice_sw          : ', Model%icice_sw
      print *, ' icliq_lw          : ', Model%icliq_lw
      print *, ' icice_lw          : ', Model%icice_lw
      print *, ' iovr              : ', Model%iovr
      print *, ' idcor             : ', Model%idcor
      print *, ' dcorr_con         : ', Model%dcorr_con
      print *, ' ictm              : ', Model%ictm
      print *, ' isubc_sw          : ', Model%isubc_sw
      print *, ' isubc_lw          : ', Model%isubc_lw
      print *, ' iswmode           : ', Model%iswmode
      print *, ' lcrick            : ', Model%lcrick
      print *, ' lcnorm            : ', Model%lcnorm
      print *, ' lnoprec           : ', Model%lnoprec
      print *, ' lwhtr             : ', Model%lwhtr
      print *, ' swhtr             : ', Model%swhtr
      print *, ' rad_hr_units      : ', Model%rad_hr_units
      print *, ' inc_minor_gas     : ', Model%inc_minor_gas
      print *, ' ipsd0             : ', Model%ipsd0
      print *, ' ipsdlim           : ', Model%ipsdlim
      print *, ' lrseeds           : ', Model%lrseeds
      print *, ' nrstreams         : ', Model%nrstreams
      print *, ' lextop            : ', Model%lextop
      if (Model%do_RRTMGP) then
        print *, ' rrtmgp_nrghice     : ', Model%rrtmgp_nrghice
        print *, ' do_GPsw_Glw        : ', Model%do_GPsw_Glw
        print *, ' active_gases       : ', Model%active_gases
        print *, ' nGases             : ', Model%ngases
        print *, ' rrtmgp_root        : ', Model%rrtmgp_root
        print *, ' lw_file_gas        : ', Model%lw_file_gas
        print *, ' lw_file_clouds     : ', Model%lw_file_clouds
        print *, ' rrtmgp_nBandsLW    : ', Model%rrtmgp_nBandsLW
        print *, ' rrtmgp_nGptsLW     : ', Model%rrtmgp_nGptsLW
        print *, ' sw_file_gas        : ', Model%sw_file_gas
        print *, ' sw_file_clouds     : ', Model%sw_file_clouds
        print *, ' rrtmgp_nBandsSW    : ', Model%rrtmgp_nBandsSW
        print *, ' rrtmgp_nGptsSW     : ', Model%rrtmgp_nGptsSW
        print *, ' doG_cldoptics      : ', Model%doG_cldoptics
        print *, ' doGP_cldoptics_PADE: ', Model%doGP_cldoptics_PADE
        print *, ' doGP_cldoptics_LUT : ', Model%doGP_cldoptics_LUT
        print *, ' use_LW_jacobian    : ', Model%use_LW_jacobian
        print *, ' damp_LW_fluxadj    : ', Model%damp_LW_fluxadj
        print *, ' lfnc_k             : ', Model%lfnc_k
        print *, ' lfnc_p0            : ', Model%lfnc_p0
        print *, ' doGP_lwscat        : ', Model%doGP_lwscat
        print *, ' doGP_sgs_cnv       : ', Model%doGP_sgs_cnv
        print *, ' doGP_sgs_mynn      : ', Model%doGP_sgs_cnv
        print *, ' doGP_smearclds     : ', Model%doGP_smearclds
        print *, ' iovr_convcld       : ', Model%iovr_convcld
        print *, ' rrtmgp_sw_phys_blksz  : ', Model%rrtmgp_sw_phys_blksz
        print *, ' rrtmgp_lw_phys_blksz  : ', Model%rrtmgp_lw_phys_blksz
      endif
      print *, ' '
      print *, 'microphysical switch'
      print *, ' imp_physics       : ', Model%imp_physics
      print *, ' '

      if (Model%imp_physics == Model%imp_physics_zhao_carr .or. Model%imp_physics == Model%imp_physics_zhao_carr_pdf) then
        print *, ' Z-C microphysical parameters'
        print *, ' psautco           : ', Model%psautco
        print *, ' prautco           : ', Model%prautco
        print *, ' evpco             : ', Model%evpco
        print *, ' wminco            : ', Model%wminco
        print *, ' '
      endif
      if (Model%imp_physics == Model%imp_physics_wsm6 .or. Model%imp_physics == Model%imp_physics_thompson) then
        print *, ' Thompson microphysical parameters'
        print *, ' ltaerosol         : ', Model%ltaerosol
        print *, ' mraerosol         : ', Model%mraerosol
        print *, ' lradar            : ', Model%lradar
        print *, ' nsfullradar_diag  : ', Model%nsfullradar_diag
        print *, ' lrefres           : ', Model%lrefres
        print *, ' ttendlim          : ', Model%ttendlim
        print *, ' ext_diag_thompson : ', Model%ext_diag_thompson
        print *, ' dt_inner          : ', Model%dt_inner
        print *, ' sedi_semi         : ', Model%sedi_semi
        print *, ' decfl             : ', Model%decfl
        print *, ' '
      endif
      if (Model%imp_physics == Model%imp_physics_nssl) then
        print *, ' NSSL microphysical parameters'
        print *, ' nssl_cccn - CCCN background CCN conc. : ', Model%nssl_cccn
        print *, ' nssl_alphah - graupel shape parameter : ', Model%nssl_alphah
        print *, ' nssl_alphahl - hail shape parameter   : ', Model%nssl_alphahl
        print *, ' nssl_alphar - rain shape parameter : ', Model%nssl_alphar
        print *, ' nssl_ehw0 - graupel-droplet collection effiency : ', Model%nssl_ehw0
        print *, ' nssl_ehlw0 - hail-droplet collection effiency : ', Model%nssl_ehlw0
        print *, ' nssl_hail_on - hail activation flag   : ', Model%nssl_hail_on
        print *, ' lradar - radar refl. flag             : ', Model%lradar
        print *, ' lrefres                : ', Model%lrefres
      endif
      if (Model%imp_physics == Model%imp_physics_mg) then
        print *, ' M-G microphysical parameters'
        print *, ' fprcp             : ', Model%fprcp
        print *, ' mg_dcs            : ', Model%mg_dcs
        print *, ' mg_qcvar          : ', Model%mg_qcvar
        print *, ' mg_ts_auto_ice    : ', Model%mg_ts_auto_ice
        print *, ' mg_alf            : ', Model%mg_alf
        print *, ' mg_qcmin          : ', Model%mg_qcmin
        print *, ' mg_rhmini         : ', Model%mg_rhmini
        print *, ' pdfflag           : ', Model%pdfflag
        print *, ' '
      endif
      if (Model%imp_physics == Model%imp_physics_gfdl) then
        print *, ' GFDL microphysical parameters'
        print *, ' GFDL MP radiation inter: ', Model%lgfdlmprad
        print *, ' lrefres                : ', Model%lrefres
        print *, ' '
      endif
      if (Model%imp_physics == Model%imp_physics_fer_hires) then
        print *, ' Ferrier-Aligo microphysical parameters'
        print *, ' spec_adv          : ', Model%spec_adv
        print *, ' rhgrd             : ', Model%rhgrd
        print *, ' icloud            : ', Model%icloud
        print *, ' '
      endif
      if (Model%num_dfi_radar>0) then
        print *, ' num_dfi_radar     : ', Model%num_dfi_radar
        print *, ' do_cap_suppress   : ', Model%do_cap_suppress
        do i = 1, dfi_radar_max_intervals+1
8888       format('  fh_dfi_radar(',I0,')   :',F12.4)
           if(Model%fh_dfi_radar(i)>-1e10) then
              print 8888,i,Model%fh_dfi_radar(i)
           endif
        enddo
9999    format('  radar_tten_limits: ', F12.4, ' ... ',F12.4)
        print 9999,Model%radar_tten_limits(1),Model%radar_tten_limits(2)
      endif
      print *, 'land/surface model parameters'
      print *, ' lsm               : ', Model%lsm
      print *, ' lsoil             : ', Model%lsoil
      print *, ' rdlai             : ', Model%rdlai
      print *, ' lsoil_lsm         : ', Model%lsoil_lsm
      if (Model%lsm==Model%lsm_noahmp) then
        print *, ' lsnow_lsm         : ', Model%lsnow_lsm
        print *, ' lsnow_lsm_lbound  : ', Model%lsnow_lsm_lbound
        print *, ' lsnow_lsm_ubound  : ', Model%lsnow_lsm_ubound
      end if
      print *, ' zs  (may be unset): ', Model%zs
      print *, ' dzs (may be unset): ', Model%dzs
      !
      print *, ' iopt_thcnd        : ', Model%iopt_thcnd
      print *, ' ua_phys           : ', Model%ua_phys
      print *, ' usemonalb         : ', Model%usemonalb
      print *, ' aoasis            : ', Model%aoasis
      print *, ' fasdas            : ', Model%fasdas
      print *, ' kice              : ', Model%kice
      print *, ' shape(pores)      : ', shape(Model%pores)
      print *, ' shape(resid)      : ', shape(Model%resid)
      print *, ' ivegsrc           : ', Model%ivegsrc
      print *, ' nvegcat           : ', Model%nvegcat
      print *, ' isot              : ', Model%isot
      print *, ' nsoilcat          : ', Model%nsoilcat

      if (Model%lsm == Model%lsm_noahmp) then
        print *, ' Noah MP LSM is used, the options are'
        print *, ' iopt_dveg         : ', Model%iopt_dveg
        print *, ' iopt_crs          : ', Model%iopt_crs
        print *, ' iopt_btr          : ', Model%iopt_btr
        print *, ' iopt_run          : ', Model%iopt_run
        print *, ' iopt_sfc          : ', Model%iopt_sfc
        print *, ' iopt_frz          : ', Model%iopt_frz
        print *, ' iopt_inf          : ', Model%iopt_inf
        print *, ' iopt_rad          : ', Model%iopt_rad
        print *, ' iopt_alb          : ', Model%iopt_alb
        print *, ' iopt_snf          : ', Model%iopt_snf
        print *, ' iopt_tbot         : ', Model%iopt_tbot
        print *, ' iopt_stc          : ', Model%iopt_stc
        print *, ' iopt_trs          : ', Model%iopt_trs
        print *, ' iopt_diag         : ', Model%iopt_diag
      elseif (Model%lsm == Model%lsm_ruc) then
        print *,' RUC Land Surface Model used'
        print *, 'The Physics options are'
        print *,' mosaic_lu   =  ',Model%mosaic_lu
        print *,' mosaic_soil =  ',Model%mosaic_soil
        print *,' isncond_opt =  ',Model%isncond_opt
        print *,' isncovr_opt =  ',Model%isncovr_opt
      endif
      print *, ' use_ufo           : ', Model%use_ufo
      print *, ' lcurr_sf          : ', Model%lcurr_sf
      print *, ' pert_cd           : ', Model%pert_cd
      print *, ' ntsflg            : ', Model%ntsflg
      print *, ' sfenth            : ', Model%sfenth
      print *, ' '
      print *, 'flake model parameters'
      print *, 'lkm                : ', Model%lkm
      print *, ' '
      print *, 'tuning parameters for physical parameterizations'
      print *, ' ras               : ', Model%ras
      if (Model%ras) then
        print *, ' psauras           : ', Model%psauras
        print *, ' prauras           : ', Model%prauras
        print *, ' wminras           : ', Model%wminras
      endif
      print *, ' flipv             : ', Model%flipv
      print *, ' trans_trac        : ', Model%trans_trac
      print *, ' old_monin         : ', Model%old_monin
      print *, ' do_gwd            : ', Model%do_gwd
      print *, ' cnvgwd            : ', Model%cnvgwd
      print *, ' do_cnvgwd         : ', Model%do_cnvgwd
      print *, ' mstrat            : ', Model%mstrat
      print *, ' moist_adj         : ', Model%moist_adj
      print *, ' cscnv             : ', Model%cscnv
      print *, ' cal_pre           : ', Model%cal_pre
      print *, ' do_aw             : ', Model%do_aw
      print *, ' flx_form          : ', Model%flx_form
      print *, ' do_shoc           : ', Model%do_shoc
      print *, ' shoc_parm         : ', Model%shoc_parm
      print *, ' shocaftcnv        : ', Model%shocaftcnv
      print *, ' shoc_cld          : ', Model%shoc_cld
      print *, ' uni_cld           : ', Model%uni_cld
      print *, ' oz_phys           : ', Model%oz_phys
      print *, ' oz_phys_2015      : ', Model%oz_phys_2015
      print *, ' h2o_phys          : ', Model%h2o_phys
      print *, ' pdfcld            : ', Model%pdfcld
      print *, ' shcnvcw           : ', Model%shcnvcw
      print *, ' redrag            : ', Model%redrag
      print *, ' hybedmf           : ', Model%hybedmf
      print *, ' satmedmf          : ', Model%satmedmf
      print *, ' isatmedmf         : ', Model%isatmedmf
      print *, ' shinhong          : ', Model%shinhong
      print *, ' do_ysu            : ', Model%do_ysu
      print *, ' dspheat           : ', Model%dspheat
      print *, ' lheatstrg         : ', Model%lheatstrg
      print *, ' lseaspray         : ', Model%lseaspray
      print *, ' cnvcld            : ', Model%cnvcld
      print *, ' random_clds       : ', Model%random_clds
      print *, ' shal_cnv          : ', Model%shal_cnv
      print *, ' imfshalcnv        : ', Model%imfshalcnv
      print *, ' imfdeepcnv        : ', Model%imfdeepcnv
      print *, ' do_deep           : ', Model%do_deep
      print *, ' nmtvr             : ', Model%nmtvr
      print *, ' jcap              : ', Model%jcap
      print *, ' cs_parm           : ', Model%cs_parm
      print *, ' flgmin            : ', Model%flgmin
      print *, ' cgwf              : ', Model%cgwf
      print *, ' ccwf              : ', Model%ccwf
      print *, ' cdmbgwd           : ', Model%cdmbgwd
      print *, ' alpha_fd          : ', Model%alpha_fd
      print *, ' psl_gwd_dx_factor : ', Model%psl_gwd_dx_factor
      print *, ' sup               : ', Model%sup
      print *, ' ctei_rm           : ', Model%ctei_rm
      print *, ' crtrh             : ', Model%crtrh
      print *, ' dlqf              : ', Model%dlqf
      print *, ' seed0             : ', Model%seed0
      print *, ' rbcr              : ', Model%rbcr
      print *, ' do_mynnedmf       : ', Model%do_mynnedmf
      print *, ' do_mynnsfclay     : ', Model%do_mynnsfclay
      print *, ' diag_flux         : ', Model%diag_flux
      print *, ' diag_log          : ', Model%diag_log
      print *, ' do_myjsfc         : ', Model%do_myjsfc
      print *, ' do_myjpbl         : ', Model%do_myjpbl
      print *, ' do_ugwp           : ', Model%do_ugwp
      print *, ' gwd_opt           : ', Model%gwd_opt
      print *, ' do_ugwp_v0           : ', Model%do_ugwp_v0
      print *, ' do_ugwp_v0_orog_only : ', Model%do_ugwp_v0_orog_only
      print *, ' do_ugwp_v0_nst_only  : ', Model%do_ugwp_v0_nst_only
      print *, ' do_gsl_drag_ls_bl    : ', Model%do_gsl_drag_ls_bl
      print *, ' do_gsl_drag_ss       : ', Model%do_gsl_drag_ss
      print *, ' do_gsl_drag_tofd     : ', Model%do_gsl_drag_tofd
      print *, ' do_gwd_opt_psl       : ', Model%do_gwd_opt_psl
      print *, ' do_ugwp_v1           : ', Model%do_ugwp_v1
      print *, ' do_ugwp_v1_orog_only : ', Model%do_ugwp_v1_orog_only
      print *, ' do_ugwp_v1_w_gsldrag : ', Model%do_ugwp_v1_w_gsldrag
      print *, ' hurr_pbl          : ', Model%hurr_pbl
      print *, ' var_ric           : ', Model%var_ric
      print *, ' coef_ric_l        : ', Model%coef_ric_l
      print *, ' coef_ric_s        : ', Model%coef_ric_s
      print *, ' '
      print *, 'Rayleigh friction'
      print *, ' prslrd0           : ', Model%prslrd0
      print *, ' ral_ts            : ', Model%ral_ts
      print *, ' '
      if (Model%imfdeepcnv >= 0) then
        print *, 'mass flux deep convection'
        print *, ' clam_deep         : ', Model%clam_deep
        print *, ' c0s_deep          : ', Model%c0s_deep
        print *, ' c1_deep           : ', Model%c1_deep
        print *, ' betal_deep        : ', Model%betal_deep
        print *, ' betas_deep        : ', Model%betas_deep
        print *, ' evef              : ', Model%evef
        print *, ' evfact_deep       : ', Model%evfact_deep
        print *, ' evfactl_deep      : ', Model%evfactl_deep
        print *, ' pgcon_deep        : ', Model%pgcon_deep
        print *, ' asolfac_deep      : ', Model%asolfac_deep
        print *, ' '
      endif
      if (Model%imfshalcnv >= 0) then
        print *, 'mass flux shallow convection'
        print *, ' clam_shal         : ', Model%clam_shal
        print *, ' c0s_shal          : ', Model%c0s_shal
        print *, ' c1_shal           : ', Model%c1_shal
        print *, ' pgcon_shal        : ', Model%pgcon_shal
        print *, ' asolfac_shal      : ', Model%asolfac_shal
      endif
      print *, ' '
      print *, 'near surface sea temperature model'
      print *, ' nst_anl           : ', Model%nst_anl
      print *, ' nstf_name         : ', Model%nstf_name
      print *, ' lsea              : ', Model%lsea
      print *, ' '
      print *, 'surface layer options'
      print *, ' sfc_z0_type       : ', Model%sfc_z0_type
      print *, ' icplocn2atm       : ', Model%icplocn2atm
      print *, ' '
      print *, 'vertical diffusion coefficients'
      print *, ' xkzm_m            : ', Model%xkzm_m
      print *, ' xkzm_h            : ', Model%xkzm_h
      print *, ' xkzm_s            : ', Model%xkzm_s
      print *, ' xkzminv           : ', Model%xkzminv
      print *, ' moninq_fac        : ', Model%moninq_fac
      print *, ' dspfac            : ', Model%dspfac
      print *, ' bl_upfr           : ', Model%bl_upfr
      print *, ' bl_dnfr           : ', Model%bl_dnfr
      print *, ' rlmx              : ', Model%rlmx
      print *, ' elmx              : ', Model%elmx
      print *, ' sfc_rlm           : ', Model%sfc_rlm
      print *, ' tc_pbl            : ', Model%tc_pbl
      print *, ' '
      print *, 'parameters for canopy heat storage parametrization'
      print *, ' h0facu            : ', Model%h0facu
      print *, ' h0facs            : ', Model%h0facs
      print *, ' '
      print *, 'stochastic physics'
      print *, ' do_sppt           : ', Model%do_sppt
      print *, ' pert_mp         : ', Model%pert_mp
      print *, ' pert_clds       : ', Model%pert_clds
      print *, ' pert_radtend    : ', Model%pert_radtend
      print *, ' do_shum           : ', Model%do_shum
      print *, ' do_skeb           : ', Model%do_skeb
      print *, ' lndp_type         : ', Model%lndp_type
      print *, ' n_var_lndp        : ', Model%n_var_lndp
      print *, ' lndp_each_step    : ', Model%lndp_each_step
      print *, ' do_spp            : ', Model%do_spp
      print *, ' n_var_spp         : ', Model%n_var_spp
      print *, ' '
      print *, 'convection'
      print *, 'betascu            : ', Model%betascu
      print *, 'betamcu            : ', Model%betamcu
      print *, 'betadcu            : ', Model%betadcu
      print *, ' '
      print *, 'cellular automata'
      print *, ' nca               : ', Model%nca
      print *, ' ncells            : ', Model%ncells
      print *, ' nlives            : ', Model%nlives
      print *, ' nca_g             : ', Model%nca_g
      print *, ' ncells_g          : ', Model%ncells_g
      print *, ' nlives_g          : ', Model%nlives_g
      print *, ' nfracseed         : ', Model%nfracseed
      print *, ' nseed_g           : ', Model%nseed_g
      print *, ' nseed             : ', Model%nseed
      print *, ' ca_global         : ', Model%ca_global
      print *, ' ca_sgs            : ', Model%ca_sgs
      print *, ' do_ca             : ', Model%do_ca
      print *, ' ca_advect         : ', Model%ca_advect
      print *, ' iseed_ca          : ', Model%iseed_ca
      print *, ' ca_smooth         : ', Model%ca_smooth
      print *, ' nspinup           : ', Model%nspinup
      print *, ' nthresh           : ', Model%nthresh
      print *, ' ca_amplitude      : ', Model%ca_amplitude
      print *, ' nsmooth           : ', Model%nsmooth
      print *, ' ca_closure        : ', Model%ca_closure
      print *, ' ca_entr           : ', Model%ca_entr
      print *, ' ca_trigger        : ', Model%ca_trigger
      print *, ' '
      print *, 'tracers'
      print *, ' tracer_names      : ', Model%tracer_names
      print *, ' ntrac             : ', Model%ntrac
      print *, ' nqrimef           : ', Model%nqrimef
      print *, ' ntqv              : ', Model%ntqv
      print *, ' ntoz              : ', Model%ntoz
      print *, ' ntcw              : ', Model%ntcw
      print *, ' ntiw              : ', Model%ntiw
      print *, ' ntrw              : ', Model%ntrw
      print *, ' ntsw              : ', Model%ntsw
      print *, ' ntgl              : ', Model%ntgl
      print *, ' nthl              : ', Model%nthl
      print *, ' ntclamt           : ', Model%ntclamt
      print *, ' ntlnc             : ', Model%ntlnc
      print *, ' ntinc             : ', Model%ntinc
      print *, ' ntrnc             : ', Model%ntrnc
      print *, ' ntsnc             : ', Model%ntsnc
      print *, ' ntgnc             : ', Model%ntgnc
      print *, ' nthnc             : ', Model%nthnc
      print *, ' ntccn             : ', Model%ntccn
      print *, ' ntccna            : ', Model%ntccna
      print *, ' ntgv              : ', Model%ntgv
      print *, ' nthv              : ', Model%nthv
      print *, ' ntrz              : ', Model%ntrz
      print *, ' ntgz              : ', Model%ntgz
      print *, ' nthz              : ', Model%nthz
      print *, ' ntke              : ', Model%ntke
      print *, ' ntsigma           : ', Model%ntsigma
      print *, ' nto               : ', Model%nto
      print *, ' nto2              : ', Model%nto2
      print *, ' ntwa              : ', Model%ntwa
      print *, ' ntia              : ', Model%ntia
      print *, ' ntfsmoke          : ', Model%ntfsmoke
      print *, ' ntsmoke           : ', Model%ntsmoke
      print *, ' ntdust            : ', Model%ntdust
      print *, ' ntcoarsepm        : ', Model%ntcoarsepm
      print *, ' nchem             : ', Model%nchem
      print *, ' ndvel             : ', Model%ndvel
      print *, ' ntchm             : ', Model%ntchm
      print *, ' ntchs             : ', Model%ntchs
      print *, ' ntche             : ', Model%ntche
      print *, ' ndchm             : ', Model%ndchm
      print *, ' ndchs             : ', Model%ndchs
      print *, ' ndche             : ', Model%ndche
      print *, ' fscav             : ', Model%fscav
      print *, ' '
      print *, 'derived totals for phy_f*d'
      print *, ' ntot2d            : ', Model%ntot2d
      print *, ' ntot3d            : ', Model%ntot3d
      print *, ' num_p2d           : ', Model%num_p2d
      print *, ' num_p3d           : ', Model%num_p3d
      print *, ' nshoc_2d          : ', Model%nshoc_2d
      print *, ' nshoc_3d          : ', Model%nshoc_3d
      print *, ' ncnvcld3d         : ', Model%ncnvcld3d
      print *, ' npdf3d            : ', Model%npdf3d
      print *, ' nctp              : ', Model%nctp
      print *, ' nkbfshoc          : ', Model%nkbfshoc
      print *, ' nahdshoc          : ', Model%nahdshoc
      print *, ' nscfshoc          : ', Model%nscfshoc
      print *, ' '
      print *, 'debug flags'
      print *, ' debug             : ', Model%debug
      print *, ' pre_rad           : ', Model%pre_rad
      print *, ' '
      print *, 'variables modified at each time step'
      print *, ' ipt               : ', Model%ipt
      print *, ' lprnt             : ', Model%lprnt
      print *, ' lsswr             : ', Model%lsswr
      print *, ' lslwr             : ', Model%lslwr
      print *, ' solhr             : ', Model%solhr
      print *, ' solcon            : ', Model%solcon
      print *, ' slag              : ', Model%slag
      print *, ' sdec              : ', Model%sdec
      print *, ' cdec              : ', Model%cdec
      print *, ' clstp             : ', Model%clstp
      print *, ' phour             : ', Model%phour
      print *, ' fhour             : ', Model%fhour
      print *, ' zhour             : ', Model%zhour
      print *, ' kdt               : ', Model%kdt
      print *, ' jdat              : ', Model%jdat
      print *, ' si                : ', Model%si
      print *, ' sec               : ', Model%sec
      print *, ' first_time_step   : ', Model%first_time_step
      print *, ' restart           : ', Model%restart
      print *, ' lsm_cold_start    : ', Model%lsm_cold_start
      print *, ' '
      print *, 'lightning threat indexes'
      print *, ' lightning_threat  : ', Model%lightning_threat
    endif

  end subroutine control_print


!----------------
! GFS_grid%create
!----------------
  subroutine grid_create (Grid, Model)

    implicit none

    class(GFS_grid_type)               :: Grid
    type(GFS_control_type), intent(in) :: Model
    integer :: IM

    IM = Model%ncols
    allocate (Grid%xlon   (IM))
    allocate (Grid%xlat   (IM))
    allocate (Grid%xlat_d (IM))
    allocate (Grid%xlon_d (IM))
    allocate (Grid%sinlat (IM))
    allocate (Grid%coslat (IM))
    allocate (Grid%area   (IM))
    allocate (Grid%dx     (IM))

    Grid%xlon   = clear_val
    Grid%xlat   = clear_val
    Grid%xlat_d = clear_val
    Grid%xlon_d = clear_val
    Grid%sinlat = clear_val
    Grid%coslat = clear_val
    Grid%area   = clear_val
    Grid%dx     = clear_val

!--- ozone active
    if ( Model%ntoz > 0 ) then
      allocate (Grid%ddy_o3    (IM))
      allocate (Grid%jindx1_o3 (IM))
      allocate (Grid%jindx2_o3 (IM))

      Grid%ddy_o3      = clear_val
      Grid%jindx1_o3   = clear_val
      Grid%jindx2_o3   = clear_val
    endif

!--- stratosphere h2o active
    if ( Model%h2o_phys ) then
      allocate (Grid%ddy_h    (IM))
      allocate (Grid%jindx1_h (IM))
      allocate (Grid%jindx2_h (IM))

      Grid%ddy_h       = clear_val
      Grid%jindx1_h    = clear_val
      Grid%jindx2_h    = clear_val
    endif

!--- iccn active
    if ( Model%iccn == 1) then
      allocate (Grid%ddy_ci    (IM))
      allocate (Grid%jindx1_ci (IM))
      allocate (Grid%jindx2_ci (IM))
      allocate (Grid%ddx_ci    (IM))
      allocate (Grid%iindx1_ci (IM))
      allocate (Grid%iindx2_ci (IM))

      Grid%ddy_ci      = clear_val
      Grid%jindx1_ci   = clear_val
      Grid%jindx2_ci   = clear_val
      Grid%ddx_ci      = clear_val
      Grid%iindx1_ci   = clear_val
      Grid%iindx2_ci   = clear_val
    endif

!--- iaerclm active
    if ( Model%iaerclm ) then
      allocate (Grid%ddy_aer   (IM))
      allocate (Grid%jindx1_aer(IM))
      allocate (Grid%jindx2_aer(IM))
      allocate (Grid%ddx_aer   (IM))
      allocate (Grid%iindx1_aer(IM))
      allocate (Grid%iindx2_aer(IM))

      Grid%ddy_aer     = clear_val
      Grid%jindx1_aer  = clear_val
      Grid%jindx2_aer  = clear_val
      Grid%ddx_aer     = clear_val
      Grid%iindx1_aer  = clear_val
      Grid%iindx2_aer  = clear_val
    endif

!---  Model%do_ugwpv1
   if ( Model%do_ugwp_v1 ) then
      allocate (Grid%ddy_j1tau  (IM))
      allocate (Grid%ddy_j2tau  (IM))
      allocate (Grid%jindx1_tau (IM))
      allocate (Grid%jindx2_tau (IM))

      Grid%ddy_j1tau   = clear_val
      Grid%ddy_j2tau   = clear_val
      Grid%jindx1_tau  = clear_val
      Grid%jindx2_tau  = clear_val
   endif

 end subroutine grid_create


!--------------------
! GFS_tbd_type%create
!--------------------
  subroutine tbd_create (Tbd, Model)

    implicit none

    class(GFS_tbd_type)                :: Tbd
    type(GFS_control_type), intent(in) :: Model
    integer :: IM

    IM = Model%ncols
!--- In
!--- sub-grid cloud radiation
    if ( Model%isubc_lw == 2 .or. Model%isubc_sw == 2 ) then
      allocate (Tbd%icsdsw (IM))
      allocate (Tbd%icsdlw (IM))
      Tbd%icsdsw = zero
      Tbd%icsdlw = zero
      if (Model%lrseeds) then
        allocate (Tbd%rseeds(IM,Model%nrstreams))
        Tbd%rseeds = zero
      endif
    endif

!--- DFI radar forcing
    nullify(Tbd%dfi_radar_tten)
    nullify(Tbd%cap_suppress)
    if(Model%num_dfi_radar>0) then
       allocate(Tbd%dfi_radar_tten(IM,Model%levs,Model%num_dfi_radar))
       Tbd%dfi_radar_tten = -20.0
       Tbd%dfi_radar_tten(:,1,:) = zero
       if(Model%do_cap_suppress) then
         allocate(Tbd%cap_suppress(IM,Model%num_dfi_radar))
         Tbd%cap_suppress(:,:) = zero
       endif
    endif

!--- ozone and stratosphere h2o needs
    allocate (Tbd%ozpl  (IM,Model%levozp,Model%oz_coeff))
    allocate (Tbd%h2opl (IM,Model%levh2o,Model%h2o_coeff))
    Tbd%h2opl = clear_val
    Tbd%ozpl  = clear_val

!--- ccn and in needs
    allocate (Tbd%in_nm  (IM,Model%levs))
    allocate (Tbd%ccn_nm (IM,Model%levs))
    Tbd%in_nm  = clear_val
    Tbd%ccn_nm = clear_val

!--- aerosol fields
    allocate (Tbd%aer_nm  (IM,Model%levs,ntrcaer))
    Tbd%aer_nm = clear_val

!--- tau_amf for  NGWs
    allocate (Tbd%tau_amf(im) )
    Tbd%tau_amf = clear_val

!--- maps of local index ix to global indices i and j for this block
    allocate (Tbd%imap (IM))
    allocate (Tbd%jmap (IM))
    Tbd%imap = 0
    Tbd%jmap = 0

    allocate (Tbd%rann (IM,Model%nrcm))
    Tbd%rann = rann_init

!--- In/Out
    allocate (Tbd%acv  (IM))
    allocate (Tbd%acvb (IM))
    allocate (Tbd%acvt (IM))

    Tbd%acv  = clear_val
    Tbd%acvb = clear_val
    Tbd%acvt = clear_val

    if (Model%cplflx .or. Model%cplchm .or. Model%cpllnd) then
      allocate (Tbd%drain_cpl (IM))
      allocate (Tbd%dsnow_cpl (IM))
      Tbd%drain_cpl = clear_val
      Tbd%dsnow_cpl = clear_val
    endif

    if (Model%do_sppt .or. Model%ca_global) then
      allocate (Tbd%dtdtnp    (IM,Model%levs))
      allocate (Tbd%dtotprcp  (IM))
      allocate (Tbd%dcnvprcp  (IM))
      Tbd%dtdtnp    = clear_val
      Tbd%dtotprcp  = clear_val
      Tbd%dcnvprcp  = clear_val
    endif

    allocate (Tbd%phy_f2d  (IM,Model%ntot2d))
    allocate (Tbd%phy_f3d  (IM,Model%levs,Model%ntot3d))
    Tbd%phy_f2d  = clear_val
    Tbd%phy_f3d  = clear_val

    if (Model%nctp > 0 .and. Model%cscnv) then
      allocate (Tbd%phy_fctd (IM,Model%nctp))
      Tbd%phy_fctd = clear_val
    endif

!   if (Model%do_shoc) Tbd%phy_f3d(:,1,Model%ntot3d-1) = 3.0
!   if (Model%do_shoc) Tbd%phy_f3d(:,:,Model%ntot3d-1) = 1.0

    allocate (Tbd%hpbl (IM))
    Tbd%hpbl     = clear_val

    if (Model%imfdeepcnv == Model%imfdeepcnv_gf .or. Model%imfdeepcnv == Model%imfdeepcnv_ntiedtke .or. Model%imfdeepcnv == Model%imfdeepcnv_samf .or. Model%imfshalcnv == Model%imfshalcnv_samf .or. Model%imfdeepcnv == Model%imfdeepcnv_c3 .or. Model%imfshalcnv == Model%imfshalcnv_c3) then
       allocate(Tbd%prevsq(IM, Model%levs))
       Tbd%prevsq = clear_val
    endif

    if (Model%imfdeepcnv .ge. 0 .or. Model%imfshalcnv .ge. 0) then
       allocate(Tbd%ud_mf(IM, Model%levs))
       Tbd%ud_mf = zero
    endif

    if (Model%imfdeepcnv == Model%imfdeepcnv_gf .or. Model%imfdeepcnv == Model%imfdeepcnv_ntiedtke .or.  Model%imfdeepcnv == Model%imfdeepcnv_c3) then
       allocate(Tbd%forcet(IM, Model%levs))
       allocate(Tbd%forceq(IM, Model%levs))
       allocate(Tbd%prevst(IM, Model%levs))
       Tbd%forcet = clear_val
       Tbd%forceq = clear_val
       Tbd%prevst = clear_val
    end if

    if (Model%imfdeepcnv == Model%imfdeepcnv_gf .or.  Model%imfdeepcnv == Model%imfdeepcnv_c3) then
       allocate(Tbd%cactiv(IM))
       allocate(Tbd%cactiv_m(IM))
       allocate(Tbd%aod_gf(IM))
       Tbd%cactiv = zero
       Tbd%cactiv_m = zero
       Tbd%aod_gf = zero
    end if

    !--- MYNN variables:
    if (Model%do_mynnedmf) then
       !print*,"Allocating all MYNN-EDMF variables:"
       allocate (Tbd%cldfra_bl (IM,Model%levs))
       allocate (Tbd%qc_bl     (IM,Model%levs))
       allocate (Tbd%qi_bl     (IM,Model%levs))
       allocate (Tbd%el_pbl    (IM,Model%levs))
       allocate (Tbd%sh3d      (IM,Model%levs))
       allocate (Tbd%sm3d      (IM,Model%levs))
       allocate (Tbd%qke       (IM,Model%levs))
       allocate (Tbd%tsq       (IM,Model%levs))
       allocate (Tbd%qsq       (IM,Model%levs))
       allocate (Tbd%cov       (IM,Model%levs))
       !print*,"Allocating all MYNN-EDMF variables:"
       Tbd%cldfra_bl     = clear_val
       Tbd%qc_bl         = clear_val
       Tbd%qi_bl         = clear_val
       Tbd%el_pbl        = clear_val
       Tbd%sh3d          = clear_val
       Tbd%sm3d          = clear_val
       Tbd%qke           = zero
       Tbd%tsq           = clear_val
       Tbd%qsq           = clear_val
       Tbd%cov           = clear_val
    end if

    ! MYJ variables
    if (Model%do_myjsfc.or.Model%do_myjpbl) then
       !print*,"Allocating all MYJ surface variables:"
       allocate (Tbd%phy_myj_qsfc   (IM))
       allocate (Tbd%phy_myj_thz0   (IM))
       allocate (Tbd%phy_myj_qz0    (IM))
       allocate (Tbd%phy_myj_uz0    (IM))
       allocate (Tbd%phy_myj_vz0    (IM))
       allocate (Tbd%phy_myj_akhs   (IM))
       allocate (Tbd%phy_myj_akms   (IM))
       allocate (Tbd%phy_myj_chkqlm (IM))
       allocate (Tbd%phy_myj_elflx  (IM))
       allocate (Tbd%phy_myj_a1u    (IM))
       allocate (Tbd%phy_myj_a1t    (IM))
       allocate (Tbd%phy_myj_a1q    (IM))
       !print*,"Allocating all MYJ schemes variables:"
       Tbd%phy_myj_qsfc   = clear_val
       Tbd%phy_myj_thz0   = clear_val
       Tbd%phy_myj_qz0    = clear_val
       Tbd%phy_myj_uz0    = clear_val
       Tbd%phy_myj_vz0    = clear_val
       Tbd%phy_myj_akhs   = clear_val
       Tbd%phy_myj_akms   = clear_val
       Tbd%phy_myj_chkqlm = clear_val
       Tbd%phy_myj_elflx  = clear_val
       Tbd%phy_myj_a1u    = clear_val
       Tbd%phy_myj_a1t    = clear_val
       Tbd%phy_myj_a1q    = clear_val
    end if

  end subroutine tbd_create


!------------------------
! GFS_cldprop_type%create
!------------------------
  subroutine cldprop_create (Cldprop, Model)

    implicit none

    class(GFS_cldprop_type)            :: Cldprop
    type(GFS_control_type), intent(in) :: Model
    integer :: IM

    IM = Model%ncols

    allocate (Cldprop%cv  (IM))
    allocate (Cldprop%cvt (IM))
    allocate (Cldprop%cvb (IM))

    Cldprop%cv  = clear_val
    Cldprop%cvt = clear_val
    Cldprop%cvb = clear_val

  end subroutine cldprop_create


!******************************************
! GFS_radtend_type%create
!******************************************
  subroutine radtend_create (Radtend, Model)

    implicit none

    class(GFS_radtend_type)            :: Radtend
    type(GFS_control_type), intent(in) :: Model
    integer :: IM

    IM = Model%ncols

    !--- Out (radiation only)
    allocate (Radtend%sfcfsw (IM))
    allocate (Radtend%sfcflw (IM))

    Radtend%sfcfsw%upfxc = clear_val
    Radtend%sfcfsw%upfx0 = clear_val
    Radtend%sfcfsw%dnfxc = clear_val
    Radtend%sfcfsw%dnfx0 = clear_val
    Radtend%sfcflw%upfxc = clear_val
    Radtend%sfcflw%upfx0 = clear_val
    Radtend%sfcflw%dnfxc = clear_val
    Radtend%sfcflw%dnfx0 = clear_val

    allocate (Radtend%htrsw  (IM,Model%levs))
    allocate (Radtend%htrlw  (IM,Model%levs))
    allocate (Radtend%sfalb  (IM))
    allocate (Radtend%coszen (IM))
    allocate (Radtend%tsflw  (IM))
    allocate (Radtend%semis  (IM))
    allocate (Radtend%ext550 (IM,Model%levs))

    Radtend%htrsw  = clear_val
    Radtend%htrlw  = clear_val
    Radtend%sfalb  = clear_val
    Radtend%coszen = clear_val
    Radtend%tsflw  = clear_val
    Radtend%semis  = clear_val
    Radtend%ext550 = clear_val

!--- In/Out (???) (radiation only)
    allocate (Radtend%coszdg (IM))

    Radtend%coszdg = clear_val

!--- In/Out (???) (physics only)
    allocate (Radtend%swhc  (IM,Model%levs))
    allocate (Radtend%lwhc  (IM,Model%levs))
    allocate (Radtend%lwhd  (IM,Model%levs,6))

    Radtend%lwhd  = clear_val
    Radtend%lwhc  = clear_val
    Radtend%swhc  = clear_val

  end subroutine radtend_create

  subroutine fill_dtidx(Model,dtend_select,itrac,icause,flag)
    implicit none
    class(GFS_control_type), intent(inout) :: Model
    character(len=*), intent(in) :: dtend_select(:)
    integer, intent(in) :: itrac
    integer, intent(in) :: icause
    logical, intent(in), optional :: flag

    character(len=100) :: name

    if(present(flag)) then
       if(.not. flag) return
    endif

    if(icause>0 .and. itrac>0) then
       if(Model%dtidx(itrac,icause)>0) then
          return ! This tendency is already allocated.
       endif

       name = 'dtend_'//trim(Model%dtend_var_labels(itrac)%name)//'_'//trim(Model%dtend_process_labels(icause)%name)

       if(fglob_list(dtend_select,trim(name))) then
          Model%ndtend = Model%ndtend+1
          Model%dtidx(itrac,icause) = Model%ndtend
          if(Model%me==Model%master) then
             print 308,'selected',trim(Model%dtend_process_labels(icause)%mod_name), trim(name), &
               trim(Model%dtend_var_labels(itrac)%desc), trim(Model%dtend_process_labels(icause)%desc), &
               trim(Model%dtend_var_labels(itrac)%unit)
          endif
       elseif(Model%me==Model%master) then
             print 308,'disabled',trim(Model%dtend_process_labels(icause)%mod_name), trim(name), &
               trim(Model%dtend_var_labels(itrac)%desc), trim(Model%dtend_process_labels(icause)%desc), &
               trim(Model%dtend_var_labels(itrac)%unit)
       endif
    endif
308 format('dtend ',A,': ',A,' ',A,' = ',A,' ',A,' (',A,')')
  end subroutine fill_dtidx

  recursive function fglob(pattern,string) result(match)
    ! Matches UNIX-style globs. A '*' matches 0 or more characters,
    ! and a '?' matches one character. Other characters must match
    ! exactly. The entire string must match, so if you want to match
    ! a substring in the middle, put '*' at the ends.
    !
    ! Spaces ARE significant, so make sure you trim() the inputs.
    !
    ! Examples:
    !
    !   fglob('dtend*_mp','dtend_temp_mp') => .true.
    !   fglob('dtend*_mp','dtend_cow_mp_dog') => .false. ! entire string must match
    !   fglob('c?w','cow') => .true.
    !   fglob('c?w','coow') => .false. ! "?" matches one char, not two
    !   fglob('c?w   ','cow  ') => .false. ! You forgot to trim() the inputs.
    implicit none
    logical :: match
    character(len=*), intent(in) :: pattern,string
    integer :: npat, nstr, ipat, istr, min_match, num_match
    logical :: match_infinity

    npat=len(pattern)
    nstr=len(string)
    ipat=1 ! Next pattern character to process
    istr=1 ! First string character not yet matched
    outer: do while(ipat<=npat)
       if_glob: if(pattern(ipat:ipat)=='*' .or. pattern(ipat:ipat)=='?') then
          ! Collect sequences of * and ? to avoid pathological cases.
          min_match=0 ! Number of "?" which is minimum number of chars to match
          match_infinity=.false. ! Do we see a "*"?
          glob_collect: do while(ipat<=npat)
             if(pattern(ipat:ipat)=='*') then
                match_infinity=.true.
             else if(pattern(ipat:ipat)=='?') then
                min_match=min_match+1
             else
                exit
             endif
             ipat=ipat+1
          end do glob_collect

          num_match=0
          glob_match: do while(istr<=len(string))
             if(num_match>=min_match) then
                if(match_infinity) then
                   if(fglob(pattern(ipat:npat),string(istr:nstr))) then
                      ! Remaining pattern matches remaining string.
                      match=.true.
                      return
                   else
                      ! Remaining pattern does NOT match, so we have
                      ! to consume another char.
                   endif
                else
                   ! This is a sequence of "?" and we matched them all.
                   cycle outer
                endif
             else
                ! Haven't consumed enough chars for all the "?" yet.
             endif
             istr=istr+1
             num_match=num_match+1
          enddo glob_match
          ! We get here if we hit the end of the string.
          if(num_match<min_match) then
             ! Number of "?" was greater than number of chars left, so match failed.
             match=.false.
             return
          elseif(ipat<=npat) then
             ! Not enough pattern to match the string.
             match=.false.
             return
          else
             ! Exact match. We're done.
             match=.true.
             return
          endif
       elseif(istr>nstr) then
          ! Not enough string left to match the pattern
          match=.false.
          return
       elseif(string(istr:istr)/=pattern(ipat:ipat)) then
          ! Exact character mismatch
          match=.false.
          return
       endif if_glob
       ! Exact character match
       istr=istr+1
       ipat=ipat+1
    end do outer
    ! We get here if we ran out of pattern. We must also hit the end of the string.
    match = istr>nstr
  end function fglob

  logical function fglob_list(patterns,string)
    ! Wrapper around fglob that returns .true. if ANY pattern
    ! matches. Unlike fglob(), patterns and strings ARE automatically
    ! trim()ed. Patterns are processed in order until one matches, one
    ! is empty, or one is '*'.
    implicit none
    character(len=*), intent(in) :: patterns(:)
    character(len=*), intent(in) :: string
    integer :: i,n,s
    fglob_list=.false.
    s=len_trim(string)
    do i=1,len(patterns)
       n=len_trim(patterns(i))
       if(n<1) then
          return ! end of pattern list
       elseif(n==1 .and. patterns(i)(1:1)=='*') then
          fglob_list=.true. ! A single "*" matches anything
          return
       else if(fglob(patterns(i)(1:n),string(1:s))) then
          fglob_list=.true.
          return
       else
       endif
    enddo
  end function fglob_list

  subroutine allocate_dtend_labels_and_causes(Model)
    implicit none
    type(GFS_control_type), intent(inout) :: Model
    integer :: i

    allocate(Model%dtend_var_labels(Model%ntracp100))
    allocate(Model%dtend_process_labels(Model%nprocess))

    Model%dtend_var_labels(1)%name = 'unallocated'
    Model%dtend_var_labels(1)%desc = 'unallocated tracer'
    Model%dtend_var_labels(1)%unit = 'kg kg-1 s-1'

    do i=2,Model%ntracp100
       Model%dtend_var_labels(i)%name = 'unknown'
       Model%dtend_var_labels(i)%desc = 'unspecified tracer'
       Model%dtend_var_labels(i)%unit = 'kg kg-1 s-1'
    enddo
    do i=1,Model%nprocess
       Model%dtend_process_labels(i)%name = 'unknown'
       Model%dtend_process_labels(i)%desc = 'unspecified tendency'
       Model%dtend_process_labels(i)%time_avg = .true.
       Model%dtend_process_labels(i)%mod_name = 'gfs_phys'
    enddo
  end subroutine allocate_dtend_labels_and_causes

  subroutine label_dtend_tracer(Model,itrac,name,desc,unit)
    implicit none
    type(GFS_control_type), intent(inout) :: Model
    integer, intent(in) :: itrac
    character(len=*), intent(in) :: name, desc
    character(len=*), intent(in) :: unit

    if(itrac<2) then
       ! Special index 1 is for unallocated tracers
       return
    endif

    Model%dtend_var_labels(itrac)%name = name
    Model%dtend_var_labels(itrac)%desc = desc
    Model%dtend_var_labels(itrac)%unit = unit
  end subroutine label_dtend_tracer

  subroutine label_dtend_cause(Model,icause,name,desc,mod_name,time_avg)
    implicit none
    type(GFS_control_type), intent(inout) :: Model
    integer, intent(in) :: icause
    character(len=*), intent(in) :: name, desc
    character(len=*), optional, intent(in) :: mod_name
    logical, optional, intent(in) :: time_avg

    Model%dtend_process_labels(icause)%name=name
    Model%dtend_process_labels(icause)%desc=desc
    if(present(mod_name)) then
       Model%dtend_process_labels(icause)%mod_name = mod_name
    else
       Model%dtend_process_labels(icause)%mod_name = "gfs_phys"
    endif
    if(present(time_avg)) then
       Model%dtend_process_labels(icause)%time_avg = time_avg
    else
       Model%dtend_process_labels(icause)%time_avg = .true.
    endif
  end subroutine label_dtend_cause

!----------------
! GFS_diag%create
!----------------
  subroutine diag_create (Diag, Model)
    use parse_tracers,    only: get_tracer_index
    class(GFS_diag_type)               :: Diag
    type(GFS_control_type), intent(in) :: Model
    integer :: IM
    logical, save :: linit
    logical :: have_pbl, have_dcnv, have_scnv, have_mp, have_oz_phys

    IM = Model%ncols

    if(Model%print_diff_pgr) then
      allocate(Diag%old_pgr(IM))
      Diag%old_pgr = clear_val
    endif

    if(Model%lightning_threat) then
       allocate (Diag%ltg1_max(IM))
       allocate (Diag%ltg2_max(IM))
       allocate (Diag%ltg3_max(IM))
       Diag%ltg1_max = zero
       Diag%ltg2_max = zero
       Diag%ltg3_max = zero
    endif

    !--- Radiation
    allocate (Diag%fluxr   (IM,Model%nfxr))
    allocate (Diag%topfsw  (IM))
    allocate (Diag%topflw  (IM))
!--- Physics
!--- In/Out
    allocate (Diag%srunoff (IM))
    allocate (Diag%evbsa   (IM))
    allocate (Diag%evcwa   (IM))
    allocate (Diag%snohfa  (IM))
    allocate (Diag%transa  (IM))
    allocate (Diag%sbsnoa  (IM))
    allocate (Diag%snowca  (IM))
    allocate (Diag%evbs    (IM))
    allocate (Diag%evcw    (IM))
    allocate (Diag%sbsno   (IM))
    allocate (Diag%trans   (IM))
    allocate (Diag%snowmt_land (IM))
    allocate (Diag%snowmt_ice  (IM))
    allocate (Diag%soilm   (IM))
    allocate (Diag%tmpmin  (IM))
    allocate (Diag%tmpmax  (IM))
    allocate (Diag%dusfc   (IM))
    allocate (Diag%dvsfc   (IM))
    allocate (Diag%dtsfc   (IM))
    allocate (Diag%dqsfc   (IM))
    allocate (Diag%totprcp (IM))
    allocate (Diag%totprcpb(IM))
    allocate (Diag%gflux   (IM))
    allocate (Diag%dlwsfc  (IM))
    allocate (Diag%ulwsfc  (IM))
    allocate (Diag%suntim  (IM))
    allocate (Diag%runoff  (IM))
    allocate (Diag%tecan   (IM))
    allocate (Diag%tetran  (IM))
    allocate (Diag%tedir   (IM))
    allocate (Diag%ep      (IM))
    allocate (Diag%cldwrk  (IM))
    allocate (Diag%dugwd   (IM))
    allocate (Diag%dvgwd   (IM))
    allocate (Diag%psmean  (IM))
    allocate (Diag%cnvprcp (IM))
    allocate (Diag%cnvprcpb(IM))
    allocate (Diag%spfhmin (IM))
    allocate (Diag%spfhmax (IM))
    allocate (Diag%u10mmax (IM))
    allocate (Diag%v10mmax (IM))
    allocate (Diag%wind10mmax (IM))
    allocate (Diag%u10max (IM))
    allocate (Diag%v10max (IM))
    allocate (Diag%spd10max (IM))
    allocate (Diag%rain    (IM))
    allocate (Diag%rainc   (IM))
    allocate (Diag%ice     (IM))
    allocate (Diag%snow    (IM))
    allocate (Diag%graupel (IM))
    allocate (Diag%totice  (IM))
    allocate (Diag%totsnw  (IM))
    allocate (Diag%totgrp  (IM))
    allocate (Diag%toticeb (IM))
    allocate (Diag%totsnwb (IM))
    allocate (Diag%totgrpb (IM))
    allocate (Diag%u10m    (IM))
    allocate (Diag%v10m    (IM))
    allocate (Diag%dpt2m   (IM))
    allocate (Diag%zlvl    (IM))
    allocate (Diag%psurf   (IM))
    allocate (Diag%pwat    (IM))
    allocate (Diag%t1      (IM))
    allocate (Diag%q1      (IM))
    allocate (Diag%u1      (IM))
    allocate (Diag%v1      (IM))
    allocate (Diag%chh     (IM))
    allocate (Diag%cmm     (IM))
    allocate (Diag%dlwsfci (IM))
    allocate (Diag%ulwsfci (IM))
    allocate (Diag%dswsfci (IM))
    allocate (Diag%nswsfci (IM))
    allocate (Diag%uswsfci (IM))
    allocate (Diag%dusfci  (IM))
    allocate (Diag%dvsfci  (IM))
    allocate (Diag%dtsfci  (IM))
    allocate (Diag%dqsfci  (IM))
    allocate (Diag%gfluxi  (IM))
    allocate (Diag%epi     (IM))
    allocate (Diag%smcwlt2 (IM))
    allocate (Diag%smcref2 (IM))
    allocate (Diag%rhonewsn1 (IM))
    allocate (Diag%frzr    (IM))
    allocate (Diag%frzrb   (IM))
    allocate (Diag%frozr   (IM))
    allocate (Diag%frozrb  (IM))
    allocate (Diag%tsnowp  (IM))
    allocate (Diag%tsnowpb (IM))
    if (.not. Model%lsm == Model%lsm_ruc) then
      allocate (Diag%wet1    (IM))
    end if
    allocate (Diag%sr       (IM))
    allocate (Diag%tdomr    (IM))
    allocate (Diag%tdomzr   (IM))
    allocate (Diag%tdomip   (IM))
    allocate (Diag%tdoms    (IM))
    allocate (Diag%zmtnblck (IM))

    if(Model%lsm == Model%lsm_noahmp) then
      allocate (Diag%paha    (IM))
      allocate (Diag%twa     (IM))
      allocate (Diag%pahi    (IM))
    endif

    ! F-A MP scheme
    if (Model%imp_physics == Model%imp_physics_fer_hires) then
     allocate (Diag%train     (IM,Model%levs))
    end if
    allocate (Diag%cldfra     (IM,Model%levr+LTP))
    allocate (Diag%cldfra2d   (IM))
    allocate (Diag%total_albedo (IM))
    allocate (Diag%lwp_ex (IM))
    allocate (Diag%iwp_ex (IM))
    allocate (Diag%lwp_fc (IM))
    allocate (Diag%iwp_fc (IM))

    !--- 3D diagnostics
    if (Model%ldiag3d) then
      allocate(Diag%dtend(IM,Model%levs,Model%ndtend))
      Diag%dtend = clear_val
      if (Model%qdiag3d) then
        allocate (Diag%upd_mf (IM,Model%levs))
        allocate (Diag%dwn_mf (IM,Model%levs))
        allocate (Diag%det_mf (IM,Model%levs))
      endif
      if (Model%oz_phys_2015) then
         allocate(Diag%do3_dt_prd( IM, Model%levs))
         allocate(Diag%do3_dt_ozmx(IM, Model%levs))
         allocate(Diag%do3_dt_temp(IM, Model%levs))
         allocate(Diag%do3_dt_ohoz(IM, Model%levs))
      endif
    endif

! UGWP
    allocate (Diag%zmtb      (IM)           )
    allocate (Diag%zogw      (IM)           )
    allocate (Diag%zlwb      (IM)           )
    allocate (Diag%tau_ogw   (IM)           )
    allocate (Diag%tau_ngw   (IM)           )
    allocate (Diag%tau_mtb   (IM)           )
    allocate (Diag%tau_tofd  (IM)           )
    allocate (Diag%dudt_gw   (IM,Model%levs))
    allocate (Diag%dvdt_gw   (IM,Model%levs))
    allocate (Diag%dtdt_gw   (IM,Model%levs))
    allocate (Diag%kdis_gw   (IM,Model%levs))

    if (Model%ldiag_ugwp) then
      allocate (Diag%du3dt_ogw  (IM,Model%levs) )
      allocate (Diag%du3dt_mtb  (IM,Model%levs) )
      allocate (Diag%du3dt_tms  (IM,Model%levs) )
      allocate (Diag%du3dt_ngw  (IM,Model%levs) )
      allocate (Diag%dv3dt_ngw  (IM,Model%levs) )
      allocate (Diag%dws3dt_ogw (IM,Model%levs) )
      allocate (Diag%dws3dt_obl (IM,Model%levs) )
      allocate (Diag%dws3dt_oss (IM,Model%levs) )
      allocate (Diag%dws3dt_ofd (IM,Model%levs) )
      allocate (Diag%ldu3dt_ogw  (IM,Model%levs) )
      allocate (Diag%ldu3dt_obl  (IM,Model%levs) )
      allocate (Diag%ldu3dt_oss  (IM,Model%levs) )
      allocate (Diag%ldu3dt_ofd  (IM,Model%levs) )
      allocate (Diag%ldu3dt_ngw (IM,Model%levs) )
      allocate (Diag%ldv3dt_ngw (IM,Model%levs) )
      allocate (Diag%ldt3dt_ngw (IM,Model%levs) )
    endif
    
    allocate (Diag%dudt_ogw  (IM,Model%levs))
    if (Model%do_ugwp_v1 .or. Model%ldiag_ugwp) then
      allocate (Diag%dvdt_ogw  (IM,Model%levs))
      allocate (Diag%dudt_obl  (IM,Model%levs))
      allocate (Diag%dvdt_obl  (IM,Model%levs))
      allocate (Diag%dudt_oss  (IM,Model%levs))
      allocate (Diag%dvdt_oss  (IM,Model%levs))
      allocate (Diag%dudt_ofd  (IM,Model%levs))
      allocate (Diag%dvdt_ofd  (IM,Model%levs))
      allocate (Diag%du_ogwcol (IM)           )
      allocate (Diag%dv_ogwcol (IM)           )
      allocate (Diag%du_oblcol (IM)           )
      allocate (Diag%dv_oblcol (IM)           )
      allocate (Diag%du_osscol (IM)           )
      allocate (Diag%dv_osscol (IM)           )
      allocate (Diag%du_ofdcol (IM)           )
      allocate (Diag%dv_ofdcol (IM)           )
      allocate (Diag%du3_ogwcol (IM)          )
      allocate (Diag%dv3_ogwcol (IM)          )
      allocate (Diag%du3_oblcol (IM)          )
      allocate (Diag%dv3_oblcol (IM)          )
      allocate (Diag%du3_osscol (IM)          )
      allocate (Diag%dv3_osscol (IM)          )
      allocate (Diag%du3_ofdcol (IM)          )
      allocate (Diag%dv3_ofdcol (IM)          )
    endif

    !--- 3D diagnostics for Thompson MP / GFDL MP
    allocate (Diag%refl_10cm(IM,Model%levs))
    allocate (Diag%max_hail_diam_sfc(IM))

    !--- New PBL Diagnostics
    allocate (Diag%dkt(IM,Model%levs))
    allocate (Diag%dku(IM,Model%levs))

    !--  New max hourly diag.
    allocate (Diag%refdmax(IM))
    allocate (Diag%refdmax263k(IM))
    allocate (Diag%t02max(IM))
    allocate (Diag%t02min(IM))
    allocate (Diag%rh02max(IM))
    allocate (Diag%rh02min(IM))
    allocate (Diag%pratemax(IM))

    !--- MYNN variables:
    if (Model%do_mynnedmf) then
      if (Model%bl_mynn_output .ne. 0) then
        allocate (Diag%edmf_a    (IM,Model%levs))
        allocate (Diag%edmf_w    (IM,Model%levs))
        allocate (Diag%edmf_qt   (IM,Model%levs))
        allocate (Diag%edmf_thl  (IM,Model%levs))
        allocate (Diag%edmf_ent  (IM,Model%levs))
        allocate (Diag%edmf_qc   (IM,Model%levs))
        allocate (Diag%sub_thl   (IM,Model%levs))
        allocate (Diag%sub_sqv   (IM,Model%levs))
        allocate (Diag%det_thl   (IM,Model%levs))
        allocate (Diag%det_sqv   (IM,Model%levs))
      endif
      if (Model%tke_budget .gt. 0) then
        allocate (Diag%dqke      (IM,Model%levs))
        allocate (Diag%qwt       (IM,Model%levs))
        allocate (Diag%qshear    (IM,Model%levs))
        allocate (Diag%qbuoy     (IM,Model%levs))
        allocate (Diag%qdiss     (IM,Model%levs))
      endif
      allocate (Diag%maxwidth  (IM))
      allocate (Diag%maxmf     (IM))
      allocate (Diag%ztop_plume(IM))
      allocate (Diag%ktop_plume(IM))
      allocate (Diag%exch_h    (IM,Model%levs))
      allocate (Diag%exch_m    (IM,Model%levs))
      if (Model%bl_mynn_output .ne. 0) then
        Diag%edmf_a        = clear_val
        Diag%edmf_w        = clear_val
        Diag%edmf_qt       = clear_val
        Diag%edmf_thl      = clear_val
        Diag%edmf_ent      = clear_val
        Diag%edmf_qc       = clear_val
        Diag%sub_thl       = clear_val
        Diag%sub_sqv       = clear_val
        Diag%det_thl       = clear_val
        Diag%det_sqv       = clear_val
      endif
      if (Model%tke_budget .gt. 0) then
        Diag%dqke          = clear_val
        Diag%qwt           = clear_val
        Diag%qshear        = clear_val
        Diag%qbuoy         = clear_val
        Diag%qdiss         = clear_val
      endif
      Diag%maxwidth      = clear_val
      Diag%maxmf         = clear_val
      Diag%ztop_plume    = clear_val
      Diag%ktop_plume    = 0
      Diag%exch_h        = clear_val
      Diag%exch_m        = clear_val
    endif

    ! Extended diagnostics for Thompson MP
    if (Model%ext_diag_thompson) then
      allocate (Diag%thompson_ext_diag3d(IM,Model%levs,Model%thompson_ext_ndiag3d))
      Diag%thompson_ext_diag3d = clear_val
    endif

    ! Air quality diagnostics
    ! -- initialize diagnostic variables
    if (Model%cplaqm) then
      allocate (Diag%aod(IM))
      Diag%aod = zero
    end if

    ! Auxiliary arrays in output for debugging
    if (Model%naux2d>0) then
      allocate (Diag%aux2d(IM,Model%naux2d))
      Diag%aux2d = clear_val
    endif
    if (Model%naux3d>0) then
      allocate (Diag%aux3d(IM,Model%levs,Model%naux3d))
      Diag%aux3d = clear_val
    endif

    call Diag%rad_zero  (Model)
!    if(Model%me==0) print *,'in diag_create, call rad_zero'
    linit = .true.
    call Diag%phys_zero (Model, linit=linit)
!    if(Model%me==0) print *,'in diag_create, call phys_zero'
    linit = .false.

  end subroutine diag_create

!-----------------------
! GFS_diag%rad_zero
!-----------------------
  subroutine diag_rad_zero(Diag, Model)
    class(GFS_diag_type)               :: Diag
    type(GFS_control_type), intent(in) :: Model

    Diag%fluxr        = zero
    Diag%topfsw%upfxc = zero
    Diag%topfsw%dnfxc = zero
    Diag%topfsw%upfx0 = zero
    Diag%topflw%upfxc = zero
    Diag%topflw%upfx0 = zero

  end subroutine diag_rad_zero

!------------------------
! GFS_diag%phys_zero
!------------------------
  subroutine diag_phys_zero (Diag, Model, linit, iauwindow_center)
    class(GFS_diag_type)               :: Diag
    type(GFS_control_type), intent(in) :: Model
    logical,optional, intent(in)       :: linit, iauwindow_center

    logical set_totprcp

    !--- In/Out
    Diag%srunoff    = zero
    Diag%evbsa      = zero
    Diag%evcwa      = zero
    Diag%snohfa     = zero
    Diag%transa     = zero
    Diag%snowca     = zero
    Diag%sbsnoa     = zero
    Diag%sbsno      = zero
    Diag%evbs       = zero
    Diag%evcw       = zero
    Diag%trans      = zero
    Diag%snowmt_land= zero
    Diag%snowmt_ice = zero
    Diag%soilm      = zero
    Diag%tmpmin     = Model%huge
    Diag%tmpmax     = zero
    Diag%dusfc      = zero
    Diag%dvsfc      = zero
    Diag%dtsfc      = zero
    Diag%dqsfc      = zero
    Diag%gflux      = zero
    Diag%dlwsfc     = zero
    Diag%ulwsfc     = zero
    Diag%suntim     = zero
    Diag%runoff     = zero
    Diag%tecan      = zero
    Diag%tetran     = zero
    Diag%tedir      = zero
    Diag%ep         = zero
    Diag%cldwrk     = zero
    Diag%dugwd      = zero
    Diag%dvgwd      = zero
    Diag%psmean     = zero
    Diag%spfhmin    = Model%huge
    Diag%spfhmax    = zero
    Diag%u10mmax    = zero
    Diag%v10mmax    = zero
    Diag%wind10mmax = zero
    Diag%u10max     = zero
    Diag%v10max     = zero
    Diag%spd10max   = zero
    Diag%rain       = zero
    Diag%rainc      = zero
    Diag%ice        = zero
    Diag%snow       = zero
    Diag%graupel    = zero

    !--- Out
    Diag%u10m       = zero
    Diag%v10m       = zero
    Diag%dpt2m      = zero
    Diag%zlvl       = zero
    Diag%psurf      = zero
    Diag%pwat       = zero
    Diag%t1         = zero
    Diag%q1         = zero
    Diag%u1         = zero
    Diag%v1         = zero
    Diag%chh        = zero
    Diag%cmm        = zero
    Diag%dlwsfci    = zero
    Diag%ulwsfci    = zero
    Diag%dswsfci    = zero
    Diag%nswsfci    = zero
    Diag%uswsfci    = zero
    Diag%dusfci     = zero
    Diag%dvsfci     = zero
    Diag%dtsfci     = zero
    Diag%dqsfci     = zero
    Diag%gfluxi     = zero
    Diag%epi        = zero
    Diag%smcwlt2    = zero
    Diag%smcref2    = zero
    if (.not. Model%lsm == Model%lsm_ruc) then
      Diag%wet1       = zero
    end if
    Diag%sr         = zero
    Diag%tdomr      = zero
    Diag%tdomzr     = zero
    Diag%tdomip     = zero
    Diag%tdoms      = zero
    Diag%zmtnblck   = zero

    if(Model%lsm == Model%lsm_noahmp)then
      Diag%paha       = zero
      Diag%twa        = zero
      Diag%pahi       = zero
    endif

    if (Model%imp_physics == Model%imp_physics_fer_hires) then
       Diag%train      = zero
    end if
    Diag%cldfra      = zero
    Diag%cldfra2d    = zero
    Diag%total_albedo = zero
    Diag%lwp_ex     = zero
    Diag%iwp_ex     = zero
    Diag%lwp_fc     = zero
    Diag%iwp_fc     = zero

    Diag%totprcpb   = zero
    Diag%cnvprcpb   = zero
    Diag%toticeb    = zero
    Diag%totsnwb    = zero
    Diag%totgrpb    = zero
    Diag%frzrb      = zero
    Diag%frozrb     = zero
    Diag%tsnowpb    = zero

    !--- MYNN variables:
    if (Model%do_mynnedmf) then
      if (Model%bl_mynn_output .ne. 0) then
        Diag%edmf_a        = clear_val
        Diag%edmf_w        = clear_val
        Diag%edmf_qt       = clear_val
        Diag%edmf_thl      = clear_val
        Diag%edmf_ent      = clear_val
        Diag%edmf_qc       = clear_val
        Diag%sub_thl       = clear_val
        Diag%sub_sqv       = clear_val
        Diag%det_thl       = clear_val
        Diag%det_sqv       = clear_val
      endif
      Diag%maxwidth      = clear_val
      Diag%maxmf         = clear_val
      Diag%ztop_plume    = clear_val
      Diag%ktop_plume    = 0
      Diag%exch_h        = clear_val
      Diag%exch_m        = clear_val
    endif

!    if(Model%me == Model%master) print *,'in diag_phys_zero, totprcpb set to 0,kdt=',Model%kdt

    if (Model%ldiag3d) then
       Diag%dtend    = zero
      if (Model%qdiag3d) then
        Diag%upd_mf   = zero
        Diag%dwn_mf   = zero
        Diag%det_mf   = zero
      endif
      if (Model%oz_phys_2015) then
         Diag%do3_dt_prd  = zero
         Diag%do3_dt_ozmx = zero
         Diag%do3_dt_temp = zero
         Diag%do3_dt_ohoz = zero
      endif
    endif

!
! UGWP
    Diag%zmtb        = zero
    Diag%zogw        = zero
    Diag%zlwb        = zero
    Diag%tau_mtb     = zero
    Diag%tau_ogw     = zero
    Diag%tau_ngw     = zero
    Diag%tau_tofd    = zero
    Diag%dudt_gw     = zero
    Diag%dvdt_gw     = zero
    Diag%dtdt_gw     = zero
    Diag%kdis_gw     = zero

    Diag%dudt_ogw    = zero
    if (Model%do_ugwp_v1 .or. Model%ldiag_ugwp) then
      Diag%dvdt_ogw    = zero
      Diag%dudt_obl    = zero
      Diag%dvdt_obl    = zero
      Diag%dudt_oss    = zero
      Diag%dvdt_oss    = zero
      Diag%dudt_ofd    = zero
      Diag%dvdt_ofd    = zero
      Diag%du_ogwcol   = zero
      Diag%dv_ogwcol   = zero
      Diag%du_oblcol   = zero
      Diag%dv_oblcol   = zero
      Diag%du_osscol   = zero
      Diag%dv_osscol   = zero
      Diag%du_ofdcol   = zero
      Diag%dv_ofdcol   = zero
      Diag%du3_ogwcol  = zero
      Diag%dv3_ogwcol  = zero
      Diag%du3_oblcol  = zero
      Diag%dv3_oblcol  = zero
      Diag%du3_osscol  = zero
      Diag%dv3_osscol  = zero
      Diag%du3_ofdcol  = zero
      Diag%dv3_ofdcol  = zero
    end if

    if (Model%ldiag_ugwp) then
      Diag%du3dt_ogw   = zero
      Diag%du3dt_mtb   = zero
      Diag%du3dt_tms   = zero
      Diag%du3dt_ngw   = zero
      Diag%dv3dt_ngw   = zero
      Diag%dws3dt_ogw  = zero
      Diag%dws3dt_obl  = zero
      Diag%dws3dt_oss  = zero
      Diag%dws3dt_ofd  = zero
      Diag%ldu3dt_ogw  = zero
      Diag%ldu3dt_obl  = zero
      Diag%ldu3dt_oss  = zero
      Diag%ldu3dt_ofd  = zero
      Diag%ldu3dt_ngw  = zero
      Diag%ldv3dt_ngw  = zero
      Diag%ldt3dt_ngw  = zero
    endif

!
!-----------------------------

! Extra PBL diagnostics
    Diag%dkt = zero
    Diag%dku = zero

! max hourly diagnostics
    Diag%refl_10cm   = -35.
    Diag%max_hail_diam_sfc = -999.
    Diag%refdmax     = -35.
    Diag%refdmax263k = -35.
    Diag%t02max      = -999.
    Diag%t02min      = 999.
    Diag%rh02max     = -999.
    Diag%rh02min     = 999.
    Diag%pratemax    = 0.
    Diag%rhonewsn1   = 200.
    set_totprcp      = .false.
    if (present(linit) ) set_totprcp = linit
    if (present(iauwindow_center) ) set_totprcp = iauwindow_center
    if (set_totprcp) then
      Diag%totprcp = zero
      Diag%cnvprcp = zero
      Diag%totice  = zero
      Diag%totsnw  = zero
      Diag%totgrp  = zero
      Diag%frzr    = zero
      Diag%frozr   = zero
      Diag%tsnowp  = zero
    endif

! GSL lightning threat indexes
    if(Model%lightning_threat) then
       Diag%ltg1_max = zero
       Diag%ltg2_max = zero
       Diag%ltg3_max = zero
    endif

  end subroutine diag_phys_zero

end module GFS_typedefs
