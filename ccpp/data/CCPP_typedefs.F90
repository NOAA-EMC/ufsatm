module CCPP_typedefs

!> \section arg_table_CCPP_typedefs Argument Table
!! \htmlinclude CCPP_typedefs.html
!!

    ! Physics kind defininitions needed for interstitial DDTs
    use machine,  only: kind_grid, kind_dyn, kind_phys

    ! Constants/dimensions needed for interstitial DDTs
    use GFS_typedefs,             only: clear_val, LTP

    ! Physics type defininitions needed for interstitial DDTs
    use module_radsw_parameters,  only: profsw_type, cmpfsw_type, NBDSW
    use module_radlw_parameters,  only: proflw_type, NBDLW
    use GFS_typedefs,             only: GFS_control_type

    implicit none

    private

    ! To ensure that these values match what's in the physics, array
    ! sizes are compared in the auto-generated physics caps in debug mode
    ! from module_radiation_aerosols
    integer, parameter :: NF_AESW = 3
    integer, parameter :: NF_AELW = 3
    integer, parameter :: NSPC    = 5
    integer, parameter :: NSPC1   = NSPC + 1
    ! from module_radiation_clouds
    integer, parameter :: NF_CLDS = 9
    ! from module_radiation_gases
    integer, parameter :: NF_VGAS = 10
    ! from module_radiation_surface
    integer, parameter :: NF_ALBD = 4

    ! GFS_interstitial_type         !< fields required to replace interstitial code in GFS_{physics,radiation}_driver.F90 in CCPP
    public GFS_interstitial_type

    ! GFDL_interstitial_type        !< fields required to replace interstitial code in FV3 dycore (fv_mapz.F90) in CCPP
    public GFDL_interstitial_type

!! \section arg_table_GFS_interstitial_type
!! \htmlinclude GFS_interstitial_type.html
!!
  type GFS_interstitial_type

    real (kind=kind_phys), pointer      :: adjsfculw_land(:)  => null()  !<
    real (kind=kind_phys), pointer      :: adjsfculw_ice(:)   => null()  !<
    real (kind=kind_phys), pointer      :: adjsfculw_water(:) => null()  !<
    real (kind=kind_phys), pointer      :: adjnirbmd(:)       => null()  !<
    real (kind=kind_phys), pointer      :: adjnirbmu(:)       => null()  !<
    real (kind=kind_phys), pointer      :: adjnirdfd(:)       => null()  !<
    real (kind=kind_phys), pointer      :: adjnirdfu(:)       => null()  !<
    real (kind=kind_phys), pointer      :: adjvisbmd(:)       => null()  !<
    real (kind=kind_phys), pointer      :: adjvisbmu(:)       => null()  !<
    real (kind=kind_phys), pointer      :: adjvisdfu(:)       => null()  !<
    real (kind=kind_phys), pointer      :: adjvisdfd(:)       => null()  !<
    real (kind=kind_phys), pointer      :: aerodp(:,:)        => null()  !<
    real (kind=kind_phys), pointer      :: alb1d(:)           => null()  !<
    real (kind=kind_phys), pointer      :: alpha(:,:)         => null()  !<
    real (kind=kind_phys), pointer      :: bexp1d(:)          => null()  !<
    real (kind=kind_phys), pointer      :: cd(:)              => null()  !<
    real (kind=kind_phys), pointer      :: cd_ice(:)          => null()  !<
    real (kind=kind_phys), pointer      :: cd_land(:)         => null()  !<
    real (kind=kind_phys), pointer      :: cd_water(:)        => null()  !<
    real (kind=kind_phys), pointer      :: cdq(:)             => null()  !<
    real (kind=kind_phys), pointer      :: cdq_ice(:)         => null()  !<
    real (kind=kind_phys), pointer      :: cdq_land(:)        => null()  !<
    real (kind=kind_phys), pointer      :: cdq_water(:)       => null()  !<
    real (kind=kind_phys), pointer      :: cf_upi(:,:)        => null()  !<
    real (kind=kind_phys), pointer      :: chh_ice(:)         => null()  !<
    real (kind=kind_phys), pointer      :: chh_land(:)        => null()  !<
    real (kind=kind_phys), pointer      :: chh_water(:)       => null()  !<
    real (kind=kind_phys), pointer      :: clcn(:,:)          => null()  !<
    real (kind=kind_phys), pointer      :: cldf(:)            => null()  !<
    real (kind=kind_phys), pointer      :: cldsa(:,:)         => null()  !<
    real (kind=kind_phys), pointer      :: cldtaulw(:,:)      => null()  !<
    real (kind=kind_phys), pointer      :: cldtausw(:,:)      => null()  !<
    real (kind=kind_phys), pointer      :: cld1d(:)           => null()  !<
    real (kind=kind_phys), pointer      :: clouds(:,:,:)      => null()  !<
    real (kind=kind_phys), pointer      :: clw(:,:,:)         => null()  !<
    real (kind=kind_phys), pointer      :: dclw(:,:,:)        => null()  !<
    real (kind=kind_phys), pointer      :: clx(:,:)           => null()  !<
    real (kind=kind_phys), pointer      :: cmm_ice(:)         => null()  !<
    real (kind=kind_phys), pointer      :: cmm_land(:)        => null()  !<
    real (kind=kind_phys), pointer      :: cmm_water(:)       => null()  !<
    real (kind=kind_phys), pointer      :: cnv_dqldt(:,:)     => null()  !<
    real (kind=kind_phys), pointer      :: cnv_fice(:,:)      => null()  !<
    real (kind=kind_phys), pointer      :: cnv_mfd(:,:)       => null()  !<
    real (kind=kind_phys), pointer      :: cnv_ndrop(:,:)     => null()  !<
    real (kind=kind_phys), pointer      :: cnv_nice(:,:)      => null()  !<
    real (kind=kind_phys), pointer      :: cnvc(:,:)          => null()  !<
    real (kind=kind_phys), pointer      :: ctei_r(:)          => null()  !<
    real (kind=kind_phys), pointer      :: ctei_rml(:)        => null()  !<
    real (kind=kind_phys), pointer      :: cumabs(:)          => null()  !<
    real (kind=kind_phys), pointer      :: dd_mf(:,:)         => null()  !<
    real (kind=kind_phys), pointer      :: de_lgth(:)         => null()  !<
    real (kind=kind_phys), pointer      :: del(:,:)           => null()  !<
    real (kind=kind_phys), pointer      :: del_gz(:,:)        => null()  !<
    real (kind=kind_phys), pointer      :: delr(:,:)          => null()  !<
    real (kind=kind_phys), pointer      :: dlength(:)         => null()  !<
    real (kind=kind_phys), pointer      :: dqdt(:,:,:)        => null()  !<
    real (kind=kind_phys), pointer      :: dqsfc1(:)          => null()  !<
    real (kind=kind_phys), pointer      :: drain(:)           => null()  !<
    real (kind=kind_phys), pointer      :: dtdt(:,:)          => null()  !<
    real (kind=kind_phys), pointer      :: dtsfc1(:)          => null()  !<
    real (kind=kind_phys), pointer      :: dtzm(:)            => null()  !<
    real (kind=kind_phys), pointer      :: dt_mf(:,:)         => null()  !<
    real (kind=kind_phys), pointer      :: dudt(:,:)          => null()  !<
    real (kind=kind_phys), pointer      :: dusfcg(:)          => null()  !<
    real (kind=kind_phys), pointer      :: dusfc1(:)          => null()  !<
    real (kind=kind_phys), pointer      :: dvdftra(:,:,:)     => null()  !<
    real (kind=kind_phys), pointer      :: dvdt(:,:)          => null()  !<
    real (kind=kind_phys), pointer      :: dvsfcg(:)          => null()  !<
    real (kind=kind_phys), pointer      :: dvsfc1(:)          => null()  !<
    real (kind=kind_phys), pointer      :: dzlyr(:,:)         => null()  !<
    real (kind=kind_phys), pointer      :: elvmax(:)          => null()  !<
    real (kind=kind_phys), pointer      :: ep1d(:)            => null()  !<
    real (kind=kind_phys), pointer      :: ep1d_ice(:)        => null()  !<
    real (kind=kind_phys), pointer      :: ep1d_land(:)       => null()  !<
    real (kind=kind_phys), pointer      :: ep1d_water(:)      => null()  !<
    real (kind=kind_phys), pointer      :: evap_ice(:)        => null()  !<
    real (kind=kind_phys), pointer      :: evap_land(:)       => null()  !<
    real (kind=kind_phys), pointer      :: evap_water(:)      => null()  !<
    real (kind=kind_phys), pointer      :: evbs(:)            => null()  !<
    real (kind=kind_phys), pointer      :: evcw(:)            => null()  !<
    real (kind=kind_phys), pointer      :: pah(:)             => null()  !<
    real (kind=kind_phys), pointer      :: ecan(:)            => null()  !<
    real (kind=kind_phys), pointer      :: etran(:)           => null()  !<
    real (kind=kind_phys), pointer      :: edir(:)            => null()  !<
    real (kind=kind_phys), pointer      :: faerlw(:,:,:,:)    => null()  !<
    real (kind=kind_phys), pointer      :: faersw(:,:,:,:)    => null()  !<
    real (kind=kind_phys), pointer      :: ffhh_ice(:)        => null()  !<
    real (kind=kind_phys), pointer      :: ffhh_land(:)       => null()  !<
    real (kind=kind_phys), pointer      :: ffhh_water(:)      => null()  !<
    real (kind=kind_phys), pointer      :: fh2(:)             => null()  !<
    real (kind=kind_phys), pointer      :: fh2_ice(:)         => null()  !<
    real (kind=kind_phys), pointer      :: fh2_land(:)        => null()  !<
    real (kind=kind_phys), pointer      :: fh2_water(:)       => null()  !<
    logical,               pointer      :: flag_cice(:)       => null()  !<
    logical,               pointer      :: flag_guess(:)      => null()  !<
    logical,               pointer      :: flag_iter(:)       => null()  !<
    logical,               pointer      :: flag_lakefreeze(:) => null()  !<
    real (kind=kind_phys), pointer      :: ffmm_ice(:)        => null()  !<
    real (kind=kind_phys), pointer      :: ffmm_land(:)       => null()  !<
    real (kind=kind_phys), pointer      :: ffmm_water(:)      => null()  !<
    real (kind=kind_phys), pointer      :: fm10(:)            => null()  !<
    real (kind=kind_phys), pointer      :: fm10_ice(:)        => null()  !<
    real (kind=kind_phys), pointer      :: fm10_land(:)       => null()  !<
    real (kind=kind_phys), pointer      :: fm10_water(:)      => null()  !<
    real (kind=kind_phys)               :: frain                         !<
    real (kind=kind_phys), pointer      :: frland(:)          => null()  !<
    real (kind=kind_phys), pointer      :: fscav(:)           => null()  !<
    real (kind=kind_phys), pointer      :: fswtr(:)           => null()  !<
    real (kind=kind_phys), pointer      :: gabsbdlw(:)        => null()  !<
    real (kind=kind_phys), pointer      :: gabsbdlw_ice(:)    => null()  !<
    real (kind=kind_phys), pointer      :: gabsbdlw_land(:)   => null()  !<
    real (kind=kind_phys), pointer      :: gabsbdlw_water(:)  => null()  !<
    real (kind=kind_phys), pointer      :: gamma(:)           => null()  !<
    real (kind=kind_phys), pointer      :: gamq(:)            => null()  !<
    real (kind=kind_phys), pointer      :: gamt(:)            => null()  !<
    real (kind=kind_phys), pointer      :: gasvmr(:,:,:)      => null()  !<
    real (kind=kind_phys), pointer      :: gflx(:)            => null()  !<
    real (kind=kind_phys), pointer      :: gflx_ice(:)        => null()  !<
    real (kind=kind_phys), pointer      :: gflx_land(:)       => null()  !<
    real (kind=kind_phys), pointer      :: gflx_water(:)      => null()  !<
    real (kind=kind_phys), pointer      :: graupelmp(:)       => null()  !<
    real (kind=kind_phys), pointer      :: gwdcu(:,:)         => null()  !<
    real (kind=kind_phys), pointer      :: gwdcv(:,:)         => null()  !<
    real (kind=kind_phys), pointer      :: zvfun(:)           => null()  !<
    real (kind=kind_phys), pointer      :: hffac(:)           => null()  !<
    real (kind=kind_phys), pointer      :: hflxq(:)           => null()  !<
    real (kind=kind_phys), pointer      :: hflx_ice(:)        => null()  !<
    real (kind=kind_phys), pointer      :: hflx_land(:)       => null()  !<
    real (kind=kind_phys), pointer      :: hflx_water(:)      => null()  !<
    !--- radiation variables that need to be carried over from radiation to physics
    real (kind=kind_phys), pointer      :: htlwc(:,:)         => null()  !<
    real (kind=kind_phys), pointer      :: htlw0(:,:)         => null()  !<
    real (kind=kind_phys), pointer      :: htswc(:,:)         => null()  !<
    real (kind=kind_phys), pointer      :: htsw0(:,:)         => null()  !<
    !
    real (kind=kind_phys), pointer      :: icemp(:)           => null()  !<
    logical,               pointer      :: dry(:)             => null()  !<
    integer,               pointer      :: idxday(:)          => null()  !<
    logical,               pointer      :: icy(:)             => null()  !<
    logical,               pointer      :: lake(:)            => null()  !<
    logical,               pointer      :: ocean(:)           => null()  !<
    integer                             :: ipr                           !<
    integer,               pointer      :: islmsk(:)          => null()  !<
    integer,               pointer      :: islmsk_cice(:)     => null()  !<
    integer                             :: itc                           !<
    logical,               pointer      :: wet(:)             => null()  !<
    integer                             :: kb                            !<
    integer,               pointer      :: kbot(:)            => null()  !<
    integer,               pointer      :: kcnv(:)            => null()  !<
    integer                             :: kd                            !<
    integer,               pointer      :: kinver(:)          => null()  !<
    integer,               pointer      :: kpbl(:)            => null()  !<
    integer                             :: kt                            !<
    integer,               pointer      :: ktop(:)            => null()  !<
    integer                             :: latidxprnt                    !<
    integer                             :: levi                          !<
    integer                             :: lmk                           !<
    integer                             :: lmp                           !<
    integer,               pointer      :: mbota(:,:)         => null()  !<
    logical                             :: mg3_as_mg2                    !<
    integer,               pointer      :: mtopa(:,:)         => null()  !<
    integer                             :: nbdlw                         !<
    integer                             :: nbdsw                         !<
    real (kind=kind_phys), pointer      :: ncgl(:,:)          => null()  !<
    real (kind=kind_phys), pointer      :: ncpr(:,:)          => null()  !<
    real (kind=kind_phys), pointer      :: ncps(:,:)          => null()  !<
    integer                             :: ncstrac                       !<
    integer                             :: nday                          !<
    integer                             :: nf_aelw                       !<
    integer                             :: nf_aesw                       !<
    integer                             :: nf_albd                       !<
    integer                             :: nn                            !<
    integer                             :: nsamftrac                     !<
    integer                             :: nscav                         !<
    integer                             :: nspc1                         !<
    integer                             :: ntcwx                         !<
    integer                             :: ntiwx                         !<
    integer                             :: ntrwx                         !<
    integer                             :: ntk                           !<
    integer                             :: ntkev                         !<
    integer                             :: nvdiff                        !<
    real (kind=kind_phys), pointer      :: oa4(:,:)           => null()  !<
    real (kind=kind_phys), pointer      :: oc(:)              => null()  !<
    real (kind=kind_phys), pointer      :: olyr(:,:)          => null()  !<
    logical              , pointer      :: otspt(:,:)         => null()  !<
    logical              , pointer      :: otsptflag(:)       => null()  !<
    integer                             :: oz_coeffp5                    !<
    logical                             :: phys_hydrostatic              !<
    real (kind=kind_phys), pointer      :: plvl(:,:)          => null()  !<
    real (kind=kind_phys), pointer      :: plyr(:,:)          => null()  !<
    real (kind=kind_phys), pointer      :: prcpmp(:)          => null()  !<
    real (kind=kind_phys), pointer      :: prnum(:,:)         => null()  !<
    real (kind=kind_phys), pointer      :: q2mp(:)            => null()  !<
    real (kind=kind_phys), pointer      :: qgl(:,:)           => null()  !<
    real (kind=kind_phys), pointer      :: qicn(:,:)          => null()  !<
    real (kind=kind_phys), pointer      :: qlcn(:,:)          => null()  !<
    real (kind=kind_phys), pointer      :: qlyr(:,:)          => null()  !<
    real (kind=kind_phys), pointer      :: qrn(:,:)           => null()  !<
    real (kind=kind_phys), pointer      :: qsnw(:,:)          => null()  !<
    real (kind=kind_phys), pointer      :: qss_ice(:)         => null()  !<
    real (kind=kind_phys), pointer      :: qss_land(:)        => null()  !<
    real (kind=kind_phys), pointer      :: qss_water(:)       => null()  !<
    logical                             :: fullradar_diag                !<
    real (kind=kind_phys)               :: raddt                         !<
    real (kind=kind_phys), pointer      :: rainmp(:)          => null()  !<
    real (kind=kind_phys), pointer      :: raincd(:)          => null()  !<
    real (kind=kind_phys), pointer      :: raincs(:)          => null()  !<
    real (kind=kind_phys), pointer      :: rainmcadj(:)       => null()  !<
    real (kind=kind_phys), pointer      :: rainp(:,:)         => null()  !<
    real (kind=kind_phys), pointer      :: rb(:)              => null()  !<
    real (kind=kind_phys), pointer      :: rb_ice(:)          => null()  !<
    real (kind=kind_phys), pointer      :: rb_land(:)         => null()  !<
    real (kind=kind_phys), pointer      :: rb_water(:)        => null()  !<
    logical                             :: max_hourly_reset              !<
    logical                             :: ext_diag_thompson_reset       !<
    real (kind=kind_phys), pointer      :: rhc(:,:)           => null()  !<
    real (kind=kind_phys), pointer      :: runoff(:)          => null()  !<
    real (kind=kind_phys), pointer      :: save_q(:,:,:)      => null()  !<
    real (kind=kind_phys), pointer      :: save_t(:,:)        => null()  !<
    real (kind=kind_phys), pointer      :: save_tcp(:,:)      => null()  !<
    real (kind=kind_phys), pointer      :: save_u(:,:)        => null()  !<
    real (kind=kind_phys), pointer      :: save_v(:,:)        => null()  !<
    real (kind=kind_phys), pointer      :: sbsno(:)           => null()  !<
    type (cmpfsw_type),    pointer      :: scmpsw(:)          => null()  !<
    real (kind=kind_phys), pointer      :: sfcalb(:,:)        => null()  !<
    real (kind=kind_phys), pointer      :: sigma(:)           => null()  !<
    real (kind=kind_phys), pointer      :: sigmaf(:)          => null()  !<
    real (kind=kind_phys), pointer      :: sigmafrac(:,:)     => null()  !<
    real (kind=kind_phys), pointer      :: sigmatot(:,:)      => null()  !<
    logical                             :: skip_macro                    !<
    real (kind=kind_phys), pointer      :: snowc(:)           => null()  !<
    real (kind=kind_phys), pointer      :: snohf(:)           => null()  !<
    real (kind=kind_phys), pointer      :: snowmp(:)          => null()  !<
    real (kind=kind_phys), pointer      :: snowmt(:)          => null()  !<
    real (kind=kind_phys), pointer      :: stress(:)          => null()  !<
    real (kind=kind_phys), pointer      :: stress_ice(:)      => null()  !<
    real (kind=kind_phys), pointer      :: stress_land(:)     => null()  !<
    real (kind=kind_phys), pointer      :: stress_water(:)    => null()  !<
    real (kind=kind_phys), pointer      :: t2mmp(:)           => null()  !<
    real (kind=kind_phys), pointer      :: ten_q(:,:,:)       => null()
    real (kind=kind_phys), pointer      :: ten_t(:,:)         => null()
    real (kind=kind_phys), pointer      :: ten_u(:,:)         => null()
    real (kind=kind_phys), pointer      :: ten_v(:,:)         => null()
    real (kind=kind_phys), pointer      :: theta(:)           => null()  !<
    real (kind=kind_phys), pointer      :: tlvl(:,:)          => null()  !<
    real (kind=kind_phys), pointer      :: tkeh(:,:)          => null()  !< vertical turbulent kinetic energy (m2/s2) at the model layer interfaces
    real (kind=kind_phys), pointer      :: tlyr(:,:)          => null()  !<
    real (kind=kind_phys), pointer      :: tprcp_ice(:)       => null()  !<
    real (kind=kind_phys), pointer      :: tprcp_land(:)      => null()  !<
    real (kind=kind_phys), pointer      :: tprcp_water(:)     => null()  !<
    integer                             :: tracers_start_index           !<
    integer                             :: tracers_total                 !<
    integer                             :: tracers_water                 !<
    logical                             :: trans_aero                    !<
    real (kind=kind_phys), pointer      :: trans(:)           => null()  !<
    real (kind=kind_phys), pointer      :: tseal(:)           => null()  !<
    real (kind=kind_phys), pointer      :: tsfa(:)            => null()  !<
    real (kind=kind_phys), pointer      :: tsfc_water(:)      => null()  !<
    real (kind=kind_phys), pointer      :: tsfg(:)            => null()  !<
    real (kind=kind_phys), pointer      :: tsurf_ice(:)       => null()  !<
    real (kind=kind_phys), pointer      :: tsurf_land(:)      => null()  !<
    real (kind=kind_phys), pointer      :: tsurf_water(:)     => null()  !<
    real (kind=kind_phys), pointer      :: ud_mf(:,:)         => null()  !<
    real (kind=kind_phys), pointer      :: uustar_ice(:)      => null()  !<
    real (kind=kind_phys), pointer      :: uustar_land(:)     => null()  !<
    real (kind=kind_phys), pointer      :: uustar_water(:)    => null()  !<
    real (kind=kind_phys), pointer      :: vdftra(:,:,:)      => null()  !<
    real (kind=kind_phys), pointer      :: vegf1d(:)          => null()  !<
    real (kind=kind_phys)               :: lndp_vgf                      !<

    real (kind=kind_phys), pointer      :: w_upi(:,:)         => null()  !<
    real (kind=kind_phys), pointer      :: wcbmax(:)          => null()  !<
    real (kind=kind_phys), pointer      :: wind(:)            => null()  !<
    real (kind=kind_phys), pointer      :: work1(:)           => null()  !<
    real (kind=kind_phys), pointer      :: work2(:)           => null()  !<
    real (kind=kind_phys), pointer      :: work3(:)           => null()  !<
    real (kind=kind_phys), pointer      :: xcosz(:)           => null()  !<
    real (kind=kind_phys), pointer      :: xlai1d(:)          => null()  !<
    real (kind=kind_phys), pointer      :: xmu(:)             => null()  !<
    real (kind=kind_phys), pointer      :: z01d(:)            => null()  !<
    real (kind=kind_phys), pointer      :: zt1d(:)            => null()  !<
    real (kind=kind_phys), pointer      :: ztmax_ice(:)       => null()  !<
    real (kind=kind_phys), pointer      :: ztmax_land(:)      => null()  !<
    real (kind=kind_phys), pointer      :: ztmax_water(:)     => null()  !<
!==================================================================================================
! UGWP - five mechnanisms of momentum deposition due to various types of GWs
! (oss, ofd, obl, ogw) + ngw = sum( sso + ngw)
!==================================================================================================
! nGWs
    real (kind=kind_phys), pointer      :: dudt_ngw(:,:)      => null()  !<
    real (kind=kind_phys), pointer      :: dvdt_ngw(:,:)      => null()  !<
    real (kind=kind_phys), pointer      :: dtdt_ngw(:,:)      => null()  !<
    real (kind=kind_phys), pointer      :: kdis_ngw(:,:)      => null()  !<

    real (kind=kind_phys), pointer      :: tau_oss(: )        => null()  !< instantaneous momentum flux due to OSS
    real (kind=kind_phys), pointer      :: tau_tofd(:)        => null()  !< instantaneous momentum flux due to TOFD
    real (kind=kind_phys), pointer      :: tau_mtb(:)         => null()  !< instantaneous momentum of mountain blocking drag
    real (kind=kind_phys), pointer      :: tau_ogw(:)         => null()  !< instantaneous momentum flux of OGWs
    real (kind=kind_phys), pointer      :: tau_ngw(:)         => null()  !< instantaneous momentum flux of NGWs

    real (kind=kind_phys), pointer      :: zngw(:)            => null()  !< launch levels of NGWs
    real (kind=kind_phys), pointer      :: zmtb(:)            => null()  !< mountain blocking height
    real (kind=kind_phys), pointer      :: zlwb(:)            => null()  !< low level wave breaking height
    real (kind=kind_phys), pointer      :: zogw(:)            => null()  !< height of OGW-launch

    real (kind=kind_phys), pointer      :: dudt_mtb(:,:)      => null()  !< daily aver u-wind tend due to mountain blocking
    real (kind=kind_phys), pointer      :: dudt_tms(:,:)      => null()  !< daily aver u-wind tend due to TMS

    ! RRTMGP
    real (kind=kind_phys), pointer      :: p_lay(:,:)                => null()  !<
    real (kind=kind_phys), pointer      :: p_lev(:,:)                => null()  !<
    real (kind=kind_phys), pointer      :: t_lev(:,:)                => null()  !<
    real (kind=kind_phys), pointer      :: t_lay(:,:)                => null()  !<
    real (kind=kind_phys), pointer      :: relhum(:,:)               => null()  !<
    real (kind=kind_phys), pointer      :: tv_lay(:,:)               => null()  !<
    real (kind=kind_phys), pointer      :: qs_lay(:,:)               => null()  !<
    real (kind=kind_phys), pointer      :: q_lay(:,:)                => null()  !<
    real (kind=kind_phys), pointer      :: deltaZ(:,:)               => null()  !<
    real (kind=kind_phys), pointer      :: deltaZc(:,:)              => null()  !< 
    real (kind=kind_phys), pointer      :: deltaP(:,:)               => null()  !< 
    real (kind=kind_phys), pointer      :: cloud_overlap_param(:,:)  => null()  !< Cloud overlap parameter
    real (kind=kind_phys), pointer      :: cnv_cloud_overlap_param(:,:) => null()  !< Convective cloud overlap parameter
    real (kind=kind_phys), pointer      :: precip_overlap_param(:,:) => null()  !< Precipitation overlap parameter
    real (kind=kind_phys), pointer      :: tracer(:,:,:)             => null()  !<
    real (kind=kind_phys), pointer      :: aerosolslw(:,:,:,:)       => null()  !< Aerosol radiative properties in each LW band.
    real (kind=kind_phys), pointer      :: aerosolssw(:,:,:,:)       => null()  !< Aerosol radiative properties in each SW band.
    real (kind=kind_phys), pointer      :: precip_frac(:,:)          => null()  !< Precipitation fraction
    real (kind=kind_phys), pointer      :: cld_cnv_frac(:,:)         => null()  !< SGS convective cloud fraction 
    real (kind=kind_phys), pointer      :: cld_cnv_lwp(:,:)          => null()  !< SGS convective cloud liquid water path
    real (kind=kind_phys), pointer      :: cld_cnv_reliq(:,:)        => null()  !< SGS convective cloud liquid effective radius
    real (kind=kind_phys), pointer      :: cld_cnv_iwp(:,:)          => null()  !< SGS convective cloud ice water path
    real (kind=kind_phys), pointer      :: cld_cnv_reice(:,:)        => null()  !< SGS convective cloud ice effecive radius
    real (kind=kind_phys), pointer      :: cld_pbl_lwp(:,:)          => null()  !< SGS PBL        cloud liquid water path 
    real (kind=kind_phys), pointer      :: cld_pbl_reliq(:,:)        => null()  !< SGS PBL        cloud liquid effective radius
    real (kind=kind_phys), pointer      :: cld_pbl_iwp(:,:)          => null()  !< SGS PBL        cloud ice water path
    real (kind=kind_phys), pointer      :: cld_pbl_reice(:,:)        => null()  !< SGS PBL        cloud ice effecive radius
    real (kind=kind_phys), pointer      :: fluxlwUP_allsky(:,:)      => null()  !< RRTMGP upward   longwave  all-sky flux profile
    real (kind=kind_phys), pointer      :: fluxlwDOWN_allsky(:,:)    => null()  !< RRTMGP downward longwave  all-sky flux profile
    real (kind=kind_phys), pointer      :: fluxlwUP_clrsky(:,:)      => null()  !< RRTMGP upward   longwave  clr-sky flux profile
    real (kind=kind_phys), pointer      :: fluxlwDOWN_clrsky(:,:)    => null()  !< RRTMGP downward longwave  clr-sky flux profile
    real (kind=kind_phys), pointer      :: fluxswUP_allsky(:,:)      => null()  !< RRTMGP upward   shortwave all-sky flux profile
    real (kind=kind_phys), pointer      :: fluxswDOWN_allsky(:,:)    => null()  !< RRTMGP downward shortwave all-sky flux profile
    real (kind=kind_phys), pointer      :: fluxswUP_clrsky(:,:)      => null()  !< RRTMGP upward   shortwave clr-sky flux profile
    real (kind=kind_phys), pointer      :: fluxswDOWN_clrsky(:,:)    => null()  !< RRTMGP downward shortwave clr-sky flux profile
    real (kind=kind_phys), pointer      :: sfc_emiss_byband(:,:)     => null()  !< 
    real (kind=kind_phys), pointer      :: sec_diff_byband(:,:)      => null()  !<
    real (kind=kind_phys), pointer      :: sfc_alb_nir_dir(:,:)      => null()  !<
    real (kind=kind_phys), pointer      :: sfc_alb_nir_dif(:,:)      => null()  !<
    real (kind=kind_phys), pointer      :: sfc_alb_uvvis_dir(:,:)    => null()  !<
    real (kind=kind_phys), pointer      :: sfc_alb_uvvis_dif(:,:)    => null()  !<
    real (kind=kind_phys), pointer      :: toa_src_lw(:,:)           => null()  !<
    real (kind=kind_phys), pointer      :: toa_src_sw(:,:)           => null()  !<
    type(proflw_type), pointer          :: flxprf_lw(:,:)            => null()  !< DDT containing RRTMGP longwave fluxes
    type(profsw_type), pointer          :: flxprf_sw(:,:)            => null()  !< DDT containing RRTMGP shortwave fluxes
    real (kind=kind_phys), pointer      :: vmr_o2(:,:)               => null()  !<
    real (kind=kind_phys), pointer      :: vmr_h2o(:,:)              => null()  !<
    real (kind=kind_phys), pointer      :: vmr_o3(:,:)               => null()  !<
    real (kind=kind_phys), pointer      :: vmr_ch4(:,:)              => null()  !<
    real (kind=kind_phys), pointer      :: vmr_n2o(:,:)              => null()  !<
    real (kind=kind_phys), pointer      :: vmr_co2(:,:)              => null()  !<

    !-- GSL drag suite
    real (kind=kind_phys), pointer      :: varss(:)           => null()  !<
    real (kind=kind_phys), pointer      :: ocss(:)            => null()  !<
    real (kind=kind_phys), pointer      :: oa4ss(:,:)         => null()  !<
    real (kind=kind_phys), pointer      :: clxss(:,:)         => null()  !<

    !-- 3D diagnostics
    integer :: rtg_ozone_index, rtg_tke_index

    contains

      procedure :: create      => gfs_interstitial_create     !<   allocate array data
      procedure :: destroy     => gfs_interstitial_destroy    !<   deallocate array data
      procedure :: reset       => gfs_interstitial_reset      !<   reset array data

  end type GFS_interstitial_type

!! \section arg_table_GFDL_interstitial_type Argument Table
!! \htmlinclude GFDL_interstitial_type.html
!!
  type GFDL_interstitial_type

     real(kind_dyn)                      :: akap
     real(kind_dyn)                      :: bdt
     real(kind_dyn), pointer             :: cappa(:,:,:)
     logical                             :: do_qa
     real(kind_dyn), pointer             :: dtdt(:,:,:)
     logical                             :: fast_mp_consv
     integer                             :: kmp
     logical                             :: last_step
     real(kind_dyn)                      :: mdt
     logical                             :: out_dt
     real(kind_dyn), pointer             :: pfull(:)
     real(kind_dyn), pointer             :: te0_2d(:,:) ! called te_2d in fv_dynamics, te0_2d in Lagrangian_to_Eulerian, te0_2d in fv_sat_adj
     real(kind_dyn), pointer             :: te0(:,:,:)  ! called dp1 in fv_dynamics, te in Lagrangian_to_Eulerian, te0 in fv_sat_adj
     real(kind_dyn)                      :: zvir
     logical                             :: do_sat_adj
     integer                             :: is
     integer                             :: ie
     integer                             :: isd
     integer                             :: ied
     integer                             :: isc1
     integer                             :: iec1
     integer                             :: isc2
     integer                             :: iec2
     integer                             :: js
     integer                             :: je
     integer                             :: jsd
     integer                             :: jed
     integer                             :: jsc1
     integer                             :: jec1
     integer                             :: jsc2
     integer                             :: jec2
     integer                             :: ng
     integer                             :: npz
     integer                             :: npzp1
     integer                             :: npzdelz
     integer                             :: npzq_con
     integer                             :: npzcappa
     real(kind_dyn),  pointer            :: delp(:,:,:)
     real(kind_dyn),  pointer            :: delz(:,:,:)
     real(kind_grid), pointer            :: area(:,:)
     real(kind_dyn),  pointer            :: peln(:,:,:)
     real(kind_dyn),  pointer            :: phis(:,:)
     real(kind_dyn),  pointer            :: pkz(:,:,:)
     real(kind_dyn),  pointer            :: pt(:,:,:)
     real(kind_dyn),  pointer            :: qvi(:,:,:,:)
     real(kind_dyn),  pointer            :: qv(:,:,:)
     real(kind_dyn),  pointer            :: ql(:,:,:)
     real(kind_dyn),  pointer            :: qi(:,:,:)
     real(kind_dyn),  pointer            :: qr(:,:,:)
     real(kind_dyn),  pointer            :: qs(:,:,:)
     real(kind_dyn),  pointer            :: qg(:,:,:)
     real(kind_dyn),  pointer            :: qc(:,:,:)
     real(kind_dyn),  pointer            :: q_con(:,:,:)
     integer                             :: nthreads
     logical                             :: hydrostatic
     integer                             :: nwat
     integer                             :: ngas
     real(kind_dyn),  pointer            :: rilist(:)
     real(kind_dyn),  pointer            :: cpilist(:)
     integer                             :: mpirank
     integer                             :: mpiroot

  contains

    procedure :: create  => gfdl_interstitial_create     !<   allocate array data
    procedure :: reset   => gfdl_interstitial_reset      !<   reset array data
    procedure :: mprint  => gfdl_interstitial_print      !<   print array data

  end type GFDL_interstitial_type

contains

!----------------------
! GFS_interstitial_type
!----------------------

  subroutine gfs_interstitial_create (Interstitial, ixs, ixe, Model)
    !
    implicit none
    !
    class(GFS_interstitial_type)       :: Interstitial
    integer,                intent(in) :: ixs, ixe
    type(GFS_control_type), intent(in) :: Model
    integer                            :: iGas
    !
    allocate (Interstitial%otspt      (Model%ntracp1,2))
    allocate (Interstitial%otsptflag  (Model%ntrac))
    ! Set up numbers of tracers for PBL, convection, etc: sets
    ! Interstitial%{nvdiff,mg3_as_mg2,nn,tracers_total,ntcwx,ntiwx,ntk,ntkev,otspt,nsamftrac,ncstrac,nscav}
    call gfs_interstitial_setup_tracers(Interstitial, Model)
    ! Allocate arrays
    allocate (Interstitial%adjsfculw_land  (ixs:ixe))
    allocate (Interstitial%adjsfculw_ice   (ixs:ixe))
    allocate (Interstitial%adjsfculw_water (ixs:ixe))
    allocate (Interstitial%adjnirbmd       (ixs:ixe))
    allocate (Interstitial%adjnirbmu       (ixs:ixe))
    allocate (Interstitial%adjnirdfd       (ixs:ixe))
    allocate (Interstitial%adjnirdfu       (ixs:ixe))
    allocate (Interstitial%adjvisbmd       (ixs:ixe))
    allocate (Interstitial%adjvisbmu       (ixs:ixe))
    allocate (Interstitial%adjvisdfu       (ixs:ixe))
    allocate (Interstitial%adjvisdfd       (ixs:ixe))
    allocate (Interstitial%aerodp          (ixs:ixe,NSPC1))
    allocate (Interstitial%alb1d           (ixs:ixe))
    if (.not. Model%do_RRTMGP) then
      ! RRTMGP uses its own cloud_overlap_param
      allocate (Interstitial%alpha         (ixs:ixe,Model%levr+LTP))
    end if
    allocate (Interstitial%bexp1d          (ixs:ixe))
    allocate (Interstitial%cd              (ixs:ixe))
    allocate (Interstitial%cd_ice          (ixs:ixe))
    allocate (Interstitial%cd_land         (ixs:ixe))
    allocate (Interstitial%cd_water        (ixs:ixe))
    allocate (Interstitial%cdq             (ixs:ixe))
    allocate (Interstitial%cdq_ice         (ixs:ixe))
    allocate (Interstitial%cdq_land        (ixs:ixe))
    allocate (Interstitial%cdq_water       (ixs:ixe))
    allocate (Interstitial%chh_ice         (ixs:ixe))
    allocate (Interstitial%chh_land        (ixs:ixe))
    allocate (Interstitial%chh_water       (ixs:ixe))
    allocate (Interstitial%cldf            (ixs:ixe))
    allocate (Interstitial%cldsa           (ixs:ixe,5))
    allocate (Interstitial%cldtaulw        (ixs:ixe,Model%levr+LTP))
    allocate (Interstitial%cldtausw        (ixs:ixe,Model%levr+LTP))
    allocate (Interstitial%cld1d           (ixs:ixe))
    allocate (Interstitial%clouds          (ixs:ixe,Model%levr+LTP,NF_CLDS))
    allocate (Interstitial%clw             (ixs:ixe,Model%levs,Interstitial%nn))
    allocate (Interstitial%dclw            (ixs:ixe,Model%levs,Interstitial%nn))
    allocate (Interstitial%clx             (ixs:ixe,4))
    allocate (Interstitial%cmm_ice         (ixs:ixe))
    allocate (Interstitial%cmm_land        (ixs:ixe))
    allocate (Interstitial%cmm_water       (ixs:ixe))
    allocate (Interstitial%cnvc            (ixs:ixe,Model%levs))
    allocate (Interstitial%ctei_r          (ixs:ixe))
    allocate (Interstitial%ctei_rml        (ixs:ixe))
    allocate (Interstitial%cumabs          (ixs:ixe))
    allocate (Interstitial%dd_mf           (ixs:ixe,Model%levs))
    allocate (Interstitial%de_lgth         (ixs:ixe))
    allocate (Interstitial%del             (ixs:ixe,Model%levs))
    allocate (Interstitial%del_gz          (ixs:ixe,Model%levs+1))
    allocate (Interstitial%delr            (ixs:ixe,Model%levr+LTP))
    allocate (Interstitial%dlength         (ixs:ixe))
    allocate (Interstitial%dqdt            (ixs:ixe,Model%levs,Model%ntrac))
    allocate (Interstitial%dqsfc1          (ixs:ixe))
    allocate (Interstitial%drain           (ixs:ixe))
    allocate (Interstitial%dtdt            (ixs:ixe,Model%levs))
    allocate (Interstitial%dtsfc1          (ixs:ixe))
    allocate (Interstitial%dt_mf           (ixs:ixe,Model%levs))
    allocate (Interstitial%dtzm            (ixs:ixe))
    allocate (Interstitial%dudt            (ixs:ixe,Model%levs))
    allocate (Interstitial%dusfcg          (ixs:ixe))
    allocate (Interstitial%dusfc1          (ixs:ixe))
    allocate (Interstitial%dvdt            (ixs:ixe,Model%levs))
    allocate (Interstitial%dvsfcg          (ixs:ixe))
    allocate (Interstitial%dvsfc1          (ixs:ixe))
    allocate (Interstitial%dvdftra         (ixs:ixe,Model%levs,Interstitial%nvdiff))
    allocate (Interstitial%dzlyr           (ixs:ixe,Model%levr+LTP))
    allocate (Interstitial%elvmax          (ixs:ixe))
    allocate (Interstitial%ep1d            (ixs:ixe))
    allocate (Interstitial%ep1d_ice        (ixs:ixe))
    allocate (Interstitial%ep1d_land       (ixs:ixe))
    allocate (Interstitial%ep1d_water      (ixs:ixe))
    allocate (Interstitial%evap_ice        (ixs:ixe))
    allocate (Interstitial%evap_land       (ixs:ixe))
    allocate (Interstitial%evap_water      (ixs:ixe))
    allocate (Interstitial%evbs            (ixs:ixe))
    allocate (Interstitial%evcw            (ixs:ixe))
    allocate (Interstitial%pah             (ixs:ixe))
    allocate (Interstitial%ecan            (ixs:ixe))
    allocate (Interstitial%etran           (ixs:ixe))
    allocate (Interstitial%edir            (ixs:ixe))
    allocate (Interstitial%faerlw          (ixs:ixe,Model%levr+LTP,NBDLW,NF_AELW))
    allocate (Interstitial%faersw          (ixs:ixe,Model%levr+LTP,NBDSW,NF_AESW))
    allocate (Interstitial%ffhh_ice        (ixs:ixe))
    allocate (Interstitial%ffhh_land       (ixs:ixe))
    allocate (Interstitial%ffhh_water      (ixs:ixe))
    allocate (Interstitial%fh2             (ixs:ixe))
    allocate (Interstitial%fh2_ice         (ixs:ixe))
    allocate (Interstitial%fh2_land        (ixs:ixe))
    allocate (Interstitial%fh2_water       (ixs:ixe))
    allocate (Interstitial%flag_cice       (ixs:ixe))
    allocate (Interstitial%flag_guess      (ixs:ixe))
    allocate (Interstitial%flag_iter       (ixs:ixe))
    allocate (Interstitial%flag_lakefreeze (ixs:ixe))
    allocate (Interstitial%ffmm_ice        (ixs:ixe))
    allocate (Interstitial%ffmm_land       (ixs:ixe))
    allocate (Interstitial%ffmm_water      (ixs:ixe))
    allocate (Interstitial%fm10            (ixs:ixe))
    allocate (Interstitial%fm10_ice        (ixs:ixe))
    allocate (Interstitial%fm10_land       (ixs:ixe))
    allocate (Interstitial%fm10_water      (ixs:ixe))
    allocate (Interstitial%frland          (ixs:ixe))
    allocate (Interstitial%fscav           (Interstitial%nscav))
    allocate (Interstitial%fswtr           (Interstitial%nscav))
    allocate (Interstitial%gabsbdlw        (ixs:ixe))
    allocate (Interstitial%gabsbdlw_ice    (ixs:ixe))
    allocate (Interstitial%gabsbdlw_land   (ixs:ixe))
    allocate (Interstitial%gabsbdlw_water  (ixs:ixe))
    allocate (Interstitial%gamma           (ixs:ixe))
    allocate (Interstitial%gamq            (ixs:ixe))
    allocate (Interstitial%gamt            (ixs:ixe))
    allocate (Interstitial%gasvmr          (ixs:ixe,Model%levr+LTP,NF_VGAS))
    allocate (Interstitial%gflx            (ixs:ixe))
    allocate (Interstitial%gflx_ice        (ixs:ixe))
    allocate (Interstitial%gflx_land       (ixs:ixe))
    allocate (Interstitial%gflx_water      (ixs:ixe))
    allocate (Interstitial%gwdcu           (ixs:ixe,Model%levs))
    allocate (Interstitial%gwdcv           (ixs:ixe,Model%levs))
    allocate (Interstitial%zvfun           (ixs:ixe))
    allocate (Interstitial%hffac           (ixs:ixe))
    allocate (Interstitial%hflxq           (ixs:ixe))
    allocate (Interstitial%hflx_ice        (ixs:ixe))
    allocate (Interstitial%hflx_land       (ixs:ixe))
    allocate (Interstitial%hflx_water      (ixs:ixe))
    allocate (Interstitial%htlwc           (ixs:ixe,Model%levr+LTP))
    allocate (Interstitial%htlw0           (ixs:ixe,Model%levr+LTP))
    allocate (Interstitial%htswc           (ixs:ixe,Model%levr+LTP))
    allocate (Interstitial%htsw0           (ixs:ixe,Model%levr+LTP))
    allocate (Interstitial%dry             (ixs:ixe))
    allocate (Interstitial%idxday          (ixs:ixe))
    allocate (Interstitial%icy             (ixs:ixe))
    allocate (Interstitial%lake            (ixs:ixe))
    allocate (Interstitial%ocean           (ixs:ixe))
    allocate (Interstitial%islmsk          (ixs:ixe))
    allocate (Interstitial%islmsk_cice     (ixs:ixe))
    allocate (Interstitial%wet             (ixs:ixe))
    allocate (Interstitial%kbot            (ixs:ixe))
    allocate (Interstitial%kcnv            (ixs:ixe))
    allocate (Interstitial%kinver          (ixs:ixe))
    allocate (Interstitial%kpbl            (ixs:ixe))
    allocate (Interstitial%ktop            (ixs:ixe))
    allocate (Interstitial%mbota           (ixs:ixe,3))
    allocate (Interstitial%mtopa           (ixs:ixe,3))
    allocate (Interstitial%oa4             (ixs:ixe,4))
    allocate (Interstitial%oc              (ixs:ixe))
    allocate (Interstitial%olyr            (ixs:ixe,Model%levr+LTP))
    allocate (Interstitial%plvl            (ixs:ixe,Model%levr+1+LTP))
    allocate (Interstitial%plyr            (ixs:ixe,Model%levr+LTP))
    allocate (Interstitial%prnum           (ixs:ixe,Model%levs))
    allocate (Interstitial%qlyr            (ixs:ixe,Model%levr+LTP))
    allocate (Interstitial%prcpmp          (ixs:ixe))
    allocate (Interstitial%qss_ice         (ixs:ixe))
    allocate (Interstitial%qss_land        (ixs:ixe))
    allocate (Interstitial%qss_water       (ixs:ixe))
    allocate (Interstitial%raincd          (ixs:ixe))
    allocate (Interstitial%raincs          (ixs:ixe))
    allocate (Interstitial%rainmcadj       (ixs:ixe))
    allocate (Interstitial%rainp           (ixs:ixe,Model%levs))
    allocate (Interstitial%rb              (ixs:ixe))
    allocate (Interstitial%rb_ice          (ixs:ixe))
    allocate (Interstitial%rb_land         (ixs:ixe))
    allocate (Interstitial%rb_water        (ixs:ixe))
    allocate (Interstitial%rhc             (ixs:ixe,Model%levs))
    allocate (Interstitial%runoff          (ixs:ixe))
    allocate (Interstitial%save_q          (ixs:ixe,Model%levs,Model%ntrac))
    allocate (Interstitial%save_t          (ixs:ixe,Model%levs))
    allocate (Interstitial%save_tcp        (ixs:ixe,Model%levs))
    allocate (Interstitial%save_u          (ixs:ixe,Model%levs))
    allocate (Interstitial%save_v          (ixs:ixe,Model%levs))
    allocate (Interstitial%sbsno           (ixs:ixe))
    allocate (Interstitial%scmpsw          (ixs:ixe))
    allocate (Interstitial%sfcalb          (ixs:ixe,NF_ALBD))
    allocate (Interstitial%sigma           (ixs:ixe))
    allocate (Interstitial%sigmaf          (ixs:ixe))
    allocate (Interstitial%sigmafrac       (ixs:ixe,Model%levs))
    allocate (Interstitial%sigmatot        (ixs:ixe,Model%levs+1))
    allocate (Interstitial%snowc           (ixs:ixe))
    allocate (Interstitial%snohf           (ixs:ixe))
    allocate (Interstitial%snowmt          (ixs:ixe))
    allocate (Interstitial%stress          (ixs:ixe))
    allocate (Interstitial%stress_ice      (ixs:ixe))
    allocate (Interstitial%stress_land     (ixs:ixe))
    allocate (Interstitial%stress_water    (ixs:ixe))
    allocate (Interstitial%ten_q           (ixs:ixe,Model%levs,Model%ntrac))
    allocate (Interstitial%ten_t           (ixs:ixe,Model%levs))
    allocate (Interstitial%ten_u           (ixs:ixe,Model%levs))
    allocate (Interstitial%ten_v           (ixs:ixe,Model%levs))
    allocate (Interstitial%theta           (ixs:ixe))
    allocate (Interstitial%tkeh            (ixs:ixe,Model%levs+1)) !Vertical turbulent kinetic energy at model layer interfaces
    allocate (Interstitial%tlvl            (ixs:ixe,Model%levr+1+LTP))
    allocate (Interstitial%tlyr            (ixs:ixe,Model%levr+LTP))
    allocate (Interstitial%tprcp_ice       (ixs:ixe))
    allocate (Interstitial%tprcp_land      (ixs:ixe))
    allocate (Interstitial%tprcp_water     (ixs:ixe))
    allocate (Interstitial%trans           (ixs:ixe))
    allocate (Interstitial%tseal           (ixs:ixe))
    allocate (Interstitial%tsfa            (ixs:ixe))
    allocate (Interstitial%tsfc_water      (ixs:ixe))
    allocate (Interstitial%tsfg            (ixs:ixe))
    allocate (Interstitial%tsurf_ice       (ixs:ixe))
    allocate (Interstitial%tsurf_land      (ixs:ixe))
    allocate (Interstitial%tsurf_water     (ixs:ixe))
    allocate (Interstitial%ud_mf           (ixs:ixe,Model%levs))
    allocate (Interstitial%uustar_ice      (ixs:ixe))
    allocate (Interstitial%uustar_land     (ixs:ixe))
    allocate (Interstitial%uustar_water    (ixs:ixe))
    allocate (Interstitial%vdftra          (ixs:ixe,Model%levs,Interstitial%nvdiff))  !GJF first dimension was set as 'IX' in GFS_physics_driver
    allocate (Interstitial%vegf1d          (ixs:ixe))
    allocate (Interstitial%wcbmax          (ixs:ixe))
    allocate (Interstitial%wind            (ixs:ixe))
    allocate (Interstitial%work1           (ixs:ixe))
    allocate (Interstitial%work2           (ixs:ixe))
    allocate (Interstitial%work3           (ixs:ixe))
    allocate (Interstitial%xcosz           (ixs:ixe))
    allocate (Interstitial%xlai1d          (ixs:ixe))
    allocate (Interstitial%xmu             (ixs:ixe))
    allocate (Interstitial%z01d            (ixs:ixe))
    allocate (Interstitial%zt1d            (ixs:ixe))
    allocate (Interstitial%ztmax_ice       (ixs:ixe))
    allocate (Interstitial%ztmax_land      (ixs:ixe))
    allocate (Interstitial%ztmax_water     (ixs:ixe))

    ! RRTMGP
    if (Model%do_RRTMGP) then
       allocate (Interstitial%tracer               (ixs:ixe, Model%levs,Model%ntrac))
       allocate (Interstitial%tv_lay               (ixs:ixe, Model%levs))
       allocate (Interstitial%relhum               (ixs:ixe, Model%levs))
       allocate (Interstitial%qs_lay               (ixs:ixe, Model%levs))
       allocate (Interstitial%q_lay                (ixs:ixe, Model%levs))
       allocate (Interstitial%deltaZ               (ixs:ixe, Model%levs))
       allocate (Interstitial%deltaZc              (ixs:ixe, Model%levs))
       allocate (Interstitial%deltaP               (ixs:ixe, Model%levs))
       allocate (Interstitial%p_lev                (ixs:ixe, Model%levs+1))
       allocate (Interstitial%p_lay                (ixs:ixe, Model%levs))
       allocate (Interstitial%t_lev                (ixs:ixe, Model%levs+1))
       allocate (Interstitial%t_lay                (ixs:ixe, Model%levs))
       allocate (Interstitial%cloud_overlap_param  (ixs:ixe, Model%levs))
       allocate (Interstitial%precip_overlap_param (ixs:ixe, Model%levs))
       allocate (Interstitial%fluxlwUP_allsky      (ixs:ixe, Model%levs+1))
       allocate (Interstitial%fluxlwDOWN_allsky    (ixs:ixe, Model%levs+1))
       allocate (Interstitial%fluxlwUP_clrsky      (ixs:ixe, Model%levs+1))
       allocate (Interstitial%fluxlwDOWN_clrsky    (ixs:ixe, Model%levs+1))
       allocate (Interstitial%fluxswUP_allsky      (ixs:ixe, Model%levs+1))
       allocate (Interstitial%fluxswDOWN_allsky    (ixs:ixe, Model%levs+1))
       allocate (Interstitial%fluxswUP_clrsky      (ixs:ixe, Model%levs+1))
       allocate (Interstitial%fluxswDOWN_clrsky    (ixs:ixe, Model%levs+1))
       allocate (Interstitial%aerosolslw           (ixs:ixe, Model%levs, Model%rrtmgp_nBandsLW, NF_AELW))
       allocate (Interstitial%aerosolssw           (ixs:ixe, Model%levs, Model%rrtmgp_nBandsSW, NF_AESW))
       allocate (Interstitial%precip_frac          (ixs:ixe, Model%levs))
       allocate (Interstitial%cld_cnv_frac         (ixs:ixe, Model%levs))
       allocate (Interstitial%cnv_cloud_overlap_param(ixs:ixe, Model%levs))
       allocate (Interstitial%cld_cnv_lwp          (ixs:ixe, Model%levs))
       allocate (Interstitial%cld_cnv_reliq        (ixs:ixe, Model%levs))
       allocate (Interstitial%cld_cnv_iwp          (ixs:ixe, Model%levs))
       allocate (Interstitial%cld_cnv_reice        (ixs:ixe, Model%levs))
       allocate (Interstitial%cld_pbl_lwp          (ixs:ixe, Model%levs))
       allocate (Interstitial%cld_pbl_reliq        (ixs:ixe, Model%levs))
       allocate (Interstitial%cld_pbl_iwp          (ixs:ixe, Model%levs))
       allocate (Interstitial%cld_pbl_reice        (ixs:ixe, Model%levs))
       allocate (Interstitial%flxprf_lw            (ixs:ixe, Model%levs+1))
       allocate (Interstitial%flxprf_sw            (ixs:ixe, Model%levs+1))
       allocate (Interstitial%sfc_emiss_byband     (Model%rrtmgp_nBandsLW,ixs:ixe))
       allocate (Interstitial%sec_diff_byband      (Model%rrtmgp_nBandsLW,ixs:ixe))
       allocate (Interstitial%sfc_alb_nir_dir      (Model%rrtmgp_nBandsSW,ixs:ixe))
       allocate (Interstitial%sfc_alb_nir_dif      (Model%rrtmgp_nBandsSW,ixs:ixe))
       allocate (Interstitial%sfc_alb_uvvis_dir    (Model%rrtmgp_nBandsSW,ixs:ixe))
       allocate (Interstitial%sfc_alb_uvvis_dif    (Model%rrtmgp_nBandsSW,ixs:ixe))
       allocate (Interstitial%toa_src_sw           (ixs:ixe,Model%rrtmgp_nGptsSW))
       allocate (Interstitial%toa_src_lw           (ixs:ixe,Model%rrtmgp_nGptsLW))
       allocate (Interstitial%vmr_o2               (ixs:ixe, Model%levs))
       allocate (Interstitial%vmr_h2o              (ixs:ixe, Model%levs))
       allocate (Interstitial%vmr_o3               (ixs:ixe, Model%levs))
       allocate (Interstitial%vmr_ch4              (ixs:ixe, Model%levs))
       allocate (Interstitial%vmr_n2o              (ixs:ixe, Model%levs))
       allocate (Interstitial%vmr_co2              (ixs:ixe, Model%levs))

    end if

! UGWP common
    allocate (Interstitial%tau_mtb         (ixs:ixe))
    allocate (Interstitial%tau_ogw         (ixs:ixe))
    allocate (Interstitial%tau_tofd        (ixs:ixe))
    allocate (Interstitial%tau_ngw         (ixs:ixe))
    allocate (Interstitial%tau_oss         (ixs:ixe))
    allocate (Interstitial%dudt_mtb        (ixs:ixe,Model%levs))
    allocate (Interstitial%dudt_tms        (ixs:ixe,Model%levs))
    allocate (Interstitial%zmtb            (ixs:ixe)           )
    allocate (Interstitial%zlwb            (ixs:ixe)           )
    allocate (Interstitial%zogw            (ixs:ixe)           )
    allocate (Interstitial%zngw            (ixs:ixe)           )

! CIRES UGWP v1
    if (Model%ldiag_ugwp .or. Model%do_ugwp_v0 .or. Model%do_ugwp_v0_nst_only &
        .or. Model%do_ugwp_v1) then
      allocate (Interstitial%dudt_ngw        (ixs:ixe,Model%levs))
      allocate (Interstitial%dvdt_ngw        (ixs:ixe,Model%levs))
      allocate (Interstitial%dtdt_ngw        (ixs:ixe,Model%levs))
      allocate (Interstitial%kdis_ngw        (ixs:ixe,Model%levs))
    end if

!-- GSL drag suite
    if (Model%gwd_opt==3 .or. Model%gwd_opt==33 .or. &
        Model%gwd_opt==2 .or. Model%gwd_opt==22 ) then
       allocate (Interstitial%varss           (ixs:ixe))
       allocate (Interstitial%ocss            (ixs:ixe))
       allocate (Interstitial%oa4ss           (ixs:ixe,4))
       allocate (Interstitial%clxss           (ixs:ixe,4))
    end if
!
    ! Allocate arrays that are conditional on physics choices
    if (Model%imp_physics == Model%imp_physics_gfdl .or. Model%imp_physics == Model%imp_physics_thompson &
         .or. Model%imp_physics == Model%imp_physics_tempo .or. Model%imp_physics == Model%imp_physics_nssl &
        ) then
       allocate (Interstitial%graupelmp  (ixs:ixe))
       allocate (Interstitial%icemp      (ixs:ixe))
       allocate (Interstitial%rainmp     (ixs:ixe))
       allocate (Interstitial%snowmp     (ixs:ixe))
    else if (Model%imp_physics == Model%imp_physics_mg) then
       allocate (Interstitial%ncgl       (ixs:ixe,Model%levs))
       allocate (Interstitial%ncpr       (ixs:ixe,Model%levs))
       allocate (Interstitial%ncps       (ixs:ixe,Model%levs))
       allocate (Interstitial%qgl        (ixs:ixe,Model%levs))
       allocate (Interstitial%qrn        (ixs:ixe,Model%levs))
       allocate (Interstitial%qsnw       (ixs:ixe,Model%levs))
       allocate (Interstitial%qlcn       (ixs:ixe,Model%levs))
       allocate (Interstitial%qicn       (ixs:ixe,Model%levs))
       allocate (Interstitial%w_upi      (ixs:ixe,Model%levs))
       allocate (Interstitial%cf_upi     (ixs:ixe,Model%levs))
       allocate (Interstitial%cnv_mfd    (ixs:ixe,Model%levs))
       allocate (Interstitial%cnv_dqldt  (ixs:ixe,Model%levs))
       allocate (Interstitial%clcn       (ixs:ixe,Model%levs))
       allocate (Interstitial%cnv_fice   (ixs:ixe,Model%levs))
       allocate (Interstitial%cnv_ndrop  (ixs:ixe,Model%levs))
       allocate (Interstitial%cnv_nice   (ixs:ixe,Model%levs))
    end if
    if (Model%lsm == Model%lsm_noahmp) then
       allocate (Interstitial%t2mmp (ixs:ixe))
       allocate (Interstitial%q2mp  (ixs:ixe))
    end if
    !
    ! Set components that do not change
    Interstitial%frain            = Model%dtf/Model%dtp
    Interstitial%ipr              = min(ixe-ixs+1,10)
    Interstitial%latidxprnt       = 1
    Interstitial%levi             = Model%levs+1
    Interstitial%lmk              = Model%levr+LTP
    Interstitial%lmp              = Model%levr+1+LTP
    Interstitial%nbdlw            = NBDLW
    Interstitial%nbdsw            = NBDSW
    Interstitial%nf_aelw          = NF_AELW
    Interstitial%nf_aesw          = NF_AESW
    Interstitial%nf_albd          = NF_ALBD
    Interstitial%nspc1            = NSPC1
    if (Model%oz_phys .or. Model%oz_phys_2015) then
      Interstitial%oz_coeffp5     = Model%oz_coeff+5
    else
      Interstitial%oz_coeffp5     = 5
    endif
    !
    Interstitial%skip_macro       = .false.
    ! The value phys_hydrostatic from dynamics does not match the
    ! hardcoded value for calling GFDL MP in GFS_physics_driver.F90,
    ! which is set to .true.
    Interstitial%phys_hydrostatic = .true.
    !
    ! Reset all other variables
    call Interstitial%reset (Model)
    !
  end subroutine gfs_interstitial_create

  subroutine gfs_interstitial_destroy (Interstitial, Model)
    !
    implicit none
    !
    class(GFS_interstitial_type)       :: Interstitial
    type(GFS_control_type), intent(in) :: Model

    deallocate (Interstitial%otspt)
    deallocate (Interstitial%otsptflag)
    ! Allocate arrays
    deallocate (Interstitial%adjsfculw_land)
    deallocate (Interstitial%adjsfculw_ice)
    deallocate (Interstitial%adjsfculw_water)
    deallocate (Interstitial%adjnirbmd)
    deallocate (Interstitial%adjnirbmu)
    deallocate (Interstitial%adjnirdfd)
    deallocate (Interstitial%adjnirdfu)
    deallocate (Interstitial%adjvisbmd)
    deallocate (Interstitial%adjvisbmu)
    deallocate (Interstitial%adjvisdfu)
    deallocate (Interstitial%adjvisdfd)
    deallocate (Interstitial%aerodp)
    deallocate (Interstitial%alb1d)
    if (.not. Model%do_RRTMGP) then
      deallocate (Interstitial%alpha)
    end if
    deallocate (Interstitial%bexp1d)
    deallocate (Interstitial%cd)
    deallocate (Interstitial%cd_ice)
    deallocate (Interstitial%cd_land)
    deallocate (Interstitial%cd_water)
    deallocate (Interstitial%cdq)
    deallocate (Interstitial%cdq_ice)
    deallocate (Interstitial%cdq_land)
    deallocate (Interstitial%cdq_water)
    deallocate (Interstitial%chh_ice)
    deallocate (Interstitial%chh_land)
    deallocate (Interstitial%chh_water)
    deallocate (Interstitial%cldf)
    deallocate (Interstitial%cldsa)
    deallocate (Interstitial%cldtaulw)
    deallocate (Interstitial%cldtausw)
    deallocate (Interstitial%cld1d)
    deallocate (Interstitial%clouds)
    deallocate (Interstitial%clw)
    deallocate (Interstitial%dclw)
    deallocate (Interstitial%clx)
    deallocate (Interstitial%cmm_ice)
    deallocate (Interstitial%cmm_land)
    deallocate (Interstitial%cmm_water)
    deallocate (Interstitial%cnvc)
    deallocate (Interstitial%ctei_r)
    deallocate (Interstitial%ctei_rml)
    deallocate (Interstitial%cumabs)
    deallocate (Interstitial%dd_mf)
    deallocate (Interstitial%de_lgth)
    deallocate (Interstitial%del)
    deallocate (Interstitial%del_gz)
    deallocate (Interstitial%delr)
    deallocate (Interstitial%dlength)
    deallocate (Interstitial%dqdt)
    deallocate (Interstitial%dqsfc1)
    deallocate (Interstitial%drain)
    deallocate (Interstitial%dtdt)
    deallocate (Interstitial%dtsfc1)
    deallocate (Interstitial%dt_mf)
    deallocate (Interstitial%dtzm)
    deallocate (Interstitial%dudt)
    deallocate (Interstitial%dusfcg)
    deallocate (Interstitial%dusfc1)
    deallocate (Interstitial%dvdt)
    deallocate (Interstitial%dvsfcg)
    deallocate (Interstitial%dvsfc1)
    deallocate (Interstitial%dvdftra)
    deallocate (Interstitial%dzlyr)
    deallocate (Interstitial%elvmax)
    deallocate (Interstitial%ep1d)
    deallocate (Interstitial%ep1d_ice)
    deallocate (Interstitial%ep1d_land)
    deallocate (Interstitial%ep1d_water)
    deallocate (Interstitial%evap_ice)
    deallocate (Interstitial%evap_land)
    deallocate (Interstitial%evap_water)
    deallocate (Interstitial%evbs)
    deallocate (Interstitial%evcw)
    deallocate (Interstitial%pah)
    deallocate (Interstitial%ecan)
    deallocate (Interstitial%etran)
    deallocate (Interstitial%edir)
    deallocate (Interstitial%faerlw)
    deallocate (Interstitial%faersw)
    deallocate (Interstitial%ffhh_ice)
    deallocate (Interstitial%ffhh_land)
    deallocate (Interstitial%ffhh_water)
    deallocate (Interstitial%fh2)
    deallocate (Interstitial%fh2_ice)
    deallocate (Interstitial%fh2_land)
    deallocate (Interstitial%fh2_water)
    deallocate (Interstitial%flag_cice)
    deallocate (Interstitial%flag_guess)
    deallocate (Interstitial%flag_iter)
    deallocate (Interstitial%flag_lakefreeze)
    deallocate (Interstitial%ffmm_ice)
    deallocate (Interstitial%ffmm_land)
    deallocate (Interstitial%ffmm_water)
    deallocate (Interstitial%fm10)
    deallocate (Interstitial%fm10_ice)
    deallocate (Interstitial%fm10_land)
    deallocate (Interstitial%fm10_water)
    deallocate (Interstitial%frland)
    deallocate (Interstitial%fscav)
    deallocate (Interstitial%fswtr)
    deallocate (Interstitial%gabsbdlw)
    deallocate (Interstitial%gabsbdlw_ice)
    deallocate (Interstitial%gabsbdlw_land)
    deallocate (Interstitial%gabsbdlw_water)
    deallocate (Interstitial%gamma)
    deallocate (Interstitial%gamq)
    deallocate (Interstitial%gamt)
    deallocate (Interstitial%gasvmr)
    deallocate (Interstitial%gflx)
    deallocate (Interstitial%gflx_ice)
    deallocate (Interstitial%gflx_land)
    deallocate (Interstitial%gflx_water)
    deallocate (Interstitial%gwdcu)
    deallocate (Interstitial%gwdcv)
    deallocate (Interstitial%zvfun)
    deallocate (Interstitial%hffac)
    deallocate (Interstitial%hflxq)
    deallocate (Interstitial%hflx_ice)
    deallocate (Interstitial%hflx_land)
    deallocate (Interstitial%hflx_water)
    deallocate (Interstitial%htlwc)
    deallocate (Interstitial%htlw0)
    deallocate (Interstitial%htswc)
    deallocate (Interstitial%htsw0)
    deallocate (Interstitial%dry)
    deallocate (Interstitial%idxday)
    deallocate (Interstitial%icy)
    deallocate (Interstitial%lake)
    deallocate (Interstitial%ocean)
    deallocate (Interstitial%islmsk)
    deallocate (Interstitial%islmsk_cice)
    deallocate (Interstitial%wet)
    deallocate (Interstitial%kbot)
    deallocate (Interstitial%kcnv)
    deallocate (Interstitial%kinver)
    deallocate (Interstitial%kpbl)
    deallocate (Interstitial%ktop)
    deallocate (Interstitial%mbota)
    deallocate (Interstitial%mtopa)
    deallocate (Interstitial%oa4)
    deallocate (Interstitial%oc)
    deallocate (Interstitial%olyr)
    deallocate (Interstitial%plvl)
    deallocate (Interstitial%plyr)
    deallocate (Interstitial%prnum)
    deallocate (Interstitial%qlyr)
    deallocate (Interstitial%prcpmp)
    deallocate (Interstitial%qss_ice)
    deallocate (Interstitial%qss_land)
    deallocate (Interstitial%qss_water)
    deallocate (Interstitial%raincd)
    deallocate (Interstitial%raincs)
    deallocate (Interstitial%rainmcadj)
    deallocate (Interstitial%rainp)
    deallocate (Interstitial%rb)
    deallocate (Interstitial%rb_ice)
    deallocate (Interstitial%rb_land)
    deallocate (Interstitial%rb_water)
    deallocate (Interstitial%rhc)
    deallocate (Interstitial%runoff)
    deallocate (Interstitial%save_q)
    deallocate (Interstitial%save_t)
    deallocate (Interstitial%save_tcp)
    deallocate (Interstitial%save_u)
    deallocate (Interstitial%save_v)
    deallocate (Interstitial%sbsno)
    deallocate (Interstitial%scmpsw)
    deallocate (Interstitial%sfcalb)
    deallocate (Interstitial%sigma)
    deallocate (Interstitial%sigmaf)
    deallocate (Interstitial%sigmafrac)
    deallocate (Interstitial%sigmatot)
    deallocate (Interstitial%snowc)
    deallocate (Interstitial%snohf)
    deallocate (Interstitial%snowmt)
    deallocate (Interstitial%stress)
    deallocate (Interstitial%stress_ice)
    deallocate (Interstitial%stress_land)
    deallocate (Interstitial%stress_water)
    deallocate (Interstitial%ten_q)
    deallocate (Interstitial%ten_t)
    deallocate (Interstitial%ten_u)
    deallocate (Interstitial%ten_v)
    deallocate (Interstitial%theta)
    deallocate (Interstitial%tkeh)
    deallocate (Interstitial%tlvl)
    deallocate (Interstitial%tlyr)
    deallocate (Interstitial%tprcp_ice)
    deallocate (Interstitial%tprcp_land)
    deallocate (Interstitial%tprcp_water)
    deallocate (Interstitial%trans)
    deallocate (Interstitial%tseal)
    deallocate (Interstitial%tsfa)
    deallocate (Interstitial%tsfc_water)
    deallocate (Interstitial%tsfg)
    deallocate (Interstitial%tsurf_ice)
    deallocate (Interstitial%tsurf_land)
    deallocate (Interstitial%tsurf_water)
    deallocate (Interstitial%ud_mf)
    deallocate (Interstitial%uustar_ice)
    deallocate (Interstitial%uustar_land)
    deallocate (Interstitial%uustar_water)
    deallocate (Interstitial%vdftra)
    deallocate (Interstitial%vegf1d)
    deallocate (Interstitial%wcbmax)
    deallocate (Interstitial%wind)
    deallocate (Interstitial%work1)
    deallocate (Interstitial%work2)
    deallocate (Interstitial%work3)
    deallocate (Interstitial%xcosz)
    deallocate (Interstitial%xlai1d)
    deallocate (Interstitial%xmu)
    deallocate (Interstitial%z01d)
    deallocate (Interstitial%zt1d)
    deallocate (Interstitial%ztmax_ice)
    deallocate (Interstitial%ztmax_land)
    deallocate (Interstitial%ztmax_water)

    ! RRTMGP
    if (Model%do_RRTMGP) then
       deallocate (Interstitial%tracer)
       deallocate (Interstitial%tv_lay)
       deallocate (Interstitial%relhum)
       deallocate (Interstitial%qs_lay)
       deallocate (Interstitial%q_lay)
       deallocate (Interstitial%deltaZ)
       deallocate (Interstitial%deltaZc)
       deallocate (Interstitial%deltaP)
       deallocate (Interstitial%p_lev)
       deallocate (Interstitial%p_lay)
       deallocate (Interstitial%t_lev)
       deallocate (Interstitial%t_lay)
       deallocate (Interstitial%cloud_overlap_param)
       deallocate (Interstitial%precip_overlap_param)
       deallocate (Interstitial%fluxlwUP_allsky)
       deallocate (Interstitial%fluxlwDOWN_allsky)
       deallocate (Interstitial%fluxlwUP_clrsky)
       deallocate (Interstitial%fluxlwDOWN_clrsky)
       deallocate (Interstitial%fluxswUP_allsky)
       deallocate (Interstitial%fluxswDOWN_allsky)
       deallocate (Interstitial%fluxswUP_clrsky)
       deallocate (Interstitial%fluxswDOWN_clrsky)
       deallocate (Interstitial%aerosolslw)
       deallocate (Interstitial%aerosolssw)
       deallocate (Interstitial%precip_frac)
       deallocate (Interstitial%cld_cnv_frac)
       deallocate (Interstitial%cnv_cloud_overlap_param)
       deallocate (Interstitial%cld_cnv_lwp)
       deallocate (Interstitial%cld_cnv_reliq)
       deallocate (Interstitial%cld_cnv_iwp)
       deallocate (Interstitial%cld_cnv_reice)
       deallocate (Interstitial%cld_pbl_lwp)
       deallocate (Interstitial%cld_pbl_reliq)
       deallocate (Interstitial%cld_pbl_iwp)
       deallocate (Interstitial%cld_pbl_reice)
       deallocate (Interstitial%flxprf_lw)
       deallocate (Interstitial%flxprf_sw)
       deallocate (Interstitial%sfc_emiss_byband)
       deallocate (Interstitial%sec_diff_byband)
       deallocate (Interstitial%sfc_alb_nir_dir)
       deallocate (Interstitial%sfc_alb_nir_dif)
       deallocate (Interstitial%sfc_alb_uvvis_dir)
       deallocate (Interstitial%sfc_alb_uvvis_dif)
       deallocate (Interstitial%toa_src_sw)
       deallocate (Interstitial%toa_src_lw)
       deallocate (Interstitial%vmr_o2)
       deallocate (Interstitial%vmr_h2o)
       deallocate (Interstitial%vmr_o3)
       deallocate (Interstitial%vmr_ch4)
       deallocate (Interstitial%vmr_n2o)
       deallocate (Interstitial%vmr_co2)
    end if

    ! UGWP common
    deallocate (Interstitial%tau_mtb)
    deallocate (Interstitial%tau_ogw)
    deallocate (Interstitial%tau_tofd)
    deallocate (Interstitial%tau_ngw)
    deallocate (Interstitial%tau_oss)
    deallocate (Interstitial%dudt_mtb)
    deallocate (Interstitial%dudt_tms)
    deallocate (Interstitial%zmtb)
    deallocate (Interstitial%zlwb)
    deallocate (Interstitial%zogw)
    deallocate (Interstitial%zngw)

    ! CIRES UGWP v1
    if (Model%ldiag_ugwp .or. Model%do_ugwp_v0 .or. Model%do_ugwp_v0_nst_only &
        .or. Model%do_ugwp_v1) then
      deallocate (Interstitial%dudt_ngw)
      deallocate (Interstitial%dvdt_ngw)
      deallocate (Interstitial%dtdt_ngw)
      deallocate (Interstitial%kdis_ngw)
    end if

    !-- GSL drag suite
    if (Model%gwd_opt==3 .or. Model%gwd_opt==33 .or. &
        Model%gwd_opt==2 .or. Model%gwd_opt==22 ) then
       deallocate (Interstitial%varss)
       deallocate (Interstitial%ocss)
       deallocate (Interstitial%oa4ss)
       deallocate (Interstitial%clxss)
    end if

    ! Allocate arrays that are conditional on physics choices
    if (Model%imp_physics == Model%imp_physics_gfdl .or. Model%imp_physics == Model%imp_physics_thompson &
         .or. Model%imp_physics == Model%imp_physics_tempo .or. Model%imp_physics == Model%imp_physics_nssl &
        ) then
       deallocate (Interstitial%graupelmp)
       deallocate (Interstitial%icemp)
       deallocate (Interstitial%rainmp)
       deallocate (Interstitial%snowmp)
    else if (Model%imp_physics == Model%imp_physics_mg) then
       deallocate (Interstitial%ncgl)
       deallocate (Interstitial%ncpr)
       deallocate (Interstitial%ncps)
       deallocate (Interstitial%qgl)
       deallocate (Interstitial%qrn)
       deallocate (Interstitial%qsnw)
       deallocate (Interstitial%qlcn)
       deallocate (Interstitial%qicn)
       deallocate (Interstitial%w_upi)
       deallocate (Interstitial%cf_upi)
       deallocate (Interstitial%cnv_mfd)
       deallocate (Interstitial%cnv_dqldt)
       deallocate (Interstitial%clcn)
       deallocate (Interstitial%cnv_fice)
       deallocate (Interstitial%cnv_ndrop)
       deallocate (Interstitial%cnv_nice)
    end if
    if (Model%lsm == Model%lsm_noahmp) then
       deallocate (Interstitial%t2mmp)
       deallocate (Interstitial%q2mp)
    end if
    
  end subroutine gfs_interstitial_destroy

  subroutine gfs_interstitial_setup_tracers(Interstitial, Model)
    !
    implicit none
    !
    class(GFS_interstitial_type)       :: Interstitial
    type(GFS_control_type), intent(in) :: Model
    integer :: n, tracers
    logical :: ltest

    !first, initialize the values (in case the values don't get initialized within if statements below)
    Interstitial%nvdiff           = Model%ntrac
    Interstitial%mg3_as_mg2       = .false.
    Interstitial%nn               = Model%ntrac + 1
    Interstitial%itc              = 0
    Interstitial%ntk              = 0
    Interstitial%ntkev            = 0
    Interstitial%tracers_total    = 0
    Interstitial%otspt(:,:)       = .true.
    Interstitial%otsptflag(:)     = .true.
    Interstitial%nsamftrac        = 0
    Interstitial%ncstrac          = 0
    Interstitial%ntcwx            = 0
    Interstitial%ntiwx            = 0
    Interstitial%ntrwx            = 0

    ! perform aerosol convective transport and PBL diffusion
    Interstitial%trans_aero = Model%cplchm .and. Model%trans_trac

    if (Model%imp_physics == Model%imp_physics_thompson .or. &
         Model%imp_physics == Model%imp_physics_tempo) then
      if (Model%ltaerosol) then
        Interstitial%nvdiff = 12
     else if (Model%mraerosol) then
        Interstitial%nvdiff = 10
      else
        Interstitial%nvdiff = 9
      endif
      if (Model%satmedmf) Interstitial%nvdiff = Interstitial%nvdiff + 1
    elseif ( Model%imp_physics == Model%imp_physics_nssl ) then
      if (Model%me == Model%master)  write(0,*) 'nssl_settings1: nvdiff,ntrac = ', Interstitial%nvdiff, Model%ntrac

      IF ( Model%nssl_hail_on ) THEN
        Interstitial%nvdiff = 16 !  Model%ntrac ! 17
      ELSE
        Interstitial%nvdiff = 13 ! turn off hail q,N, and volume
      ENDIF
      ! write(*,*) 'NSSL: nvdiff, ntrac = ',Interstitial%nvdiff, Model%ntrac
      if (Model%satmedmf) Interstitial%nvdiff = Interstitial%nvdiff + 1
      IF ( Model%nssl_ccn_on ) THEN
        Interstitial%nvdiff = Interstitial%nvdiff + 1
      ENDIF
      if (Model%me == Model%master)  write(0,*) 'nssl_settings2: nvdiff,ntrac = ', Interstitial%nvdiff, Model%ntrac
    
    elseif (Model%imp_physics == Model%imp_physics_wsm6) then
      Interstitial%nvdiff = Model%ntrac -3
      if (Model%satmedmf) Interstitial%nvdiff = Interstitial%nvdiff + 1
    elseif (Model%ntclamt > 0) then             ! for GFDL MP don't diffuse cloud amount
      Interstitial%nvdiff = Model%ntrac - 1
    endif

    if (Model%imp_physics == Model%imp_physics_mg) then
      if (abs(Model%fprcp) == 1) then
        Interstitial%mg3_as_mg2 = .false.
      elseif (Model%fprcp >= 2) then
        if(Model%ntgl > 0 .and. (Model%mg_do_graupel .or. Model%mg_do_hail)) then
          Interstitial%mg3_as_mg2 = .false.
        else                              ! MG3 code run without graupel/hail i.e. as MG2
          Interstitial%mg3_as_mg2 = .true.
        endif
      endif
    endif

    Interstitial%nscav = Model%ntrac - Model%ncnd + 2

    if (Interstitial%nvdiff == Model%ntrac) then
      Interstitial%ntcwx = Model%ntcw
      Interstitial%ntiwx = Model%ntiw
      Interstitial%ntrwx = Model%ntrw
    else
      if (Model%imp_physics == Model%imp_physics_wsm6) then
        Interstitial%ntcwx = 2
        Interstitial%ntiwx = 3
     elseif (Model%imp_physics == Model%imp_physics_thompson .or. &
          Model%imp_physics == Model%imp_physics_tempo) then
        Interstitial%ntcwx = 2
        Interstitial%ntiwx = 3
        Interstitial%ntrwx = 4
      elseif (Model%imp_physics == Model%imp_physics_nssl) then
        Interstitial%ntcwx = 2
        Interstitial%ntiwx = 3
        Interstitial%ntrwx = 4
      elseif (Model%imp_physics == Model%imp_physics_gfdl) then
        Interstitial%ntcwx = 2
        Interstitial%ntiwx = 3
        Interstitial%ntrwx = 4
      ! F-A MP scheme
      elseif (Model%imp_physics == Model%imp_physics_fer_hires) then
        Interstitial%ntcwx = 2
        Interstitial%ntiwx = 3
        Interstitial%ntrwx = 4
      elseif (Model%imp_physics == Model%imp_physics_mg) then
        Interstitial%ntcwx = 2
        Interstitial%ntiwx = 3
        Interstitial%ntrwx = 4
      elseif (Model%imp_physics == Model%imp_physics_zhao_carr) then
        Interstitial%ntcwx = 2
      endif
    endif

    if (Model%cplchm) then
      ! Only the following microphysics schemes are supported with coupled chemistry
      if (Model%imp_physics == Model%imp_physics_zhao_carr) then
        Interstitial%nvdiff = 3
      elseif (Model%imp_physics == Model%imp_physics_mg) then
        if (Model%ntgl > 0) then
          Interstitial%nvdiff = 12
        else
          Interstitial%nvdiff = 10
        endif
      elseif (Model%imp_physics == Model%imp_physics_gfdl) then
        Interstitial%nvdiff = 7
     elseif (Model%imp_physics == Model%imp_physics_thompson .or. &
          Model%imp_physics == Model%imp_physics_tempo) then
        if (Model%ltaerosol) then
          Interstitial%nvdiff = 12
        else if (Model%mraerosol) then
          Interstitial%nvdiff = 10
        else
          Interstitial%nvdiff = 9
        endif
      else
        write(0,*) "Selected microphysics scheme is not supported when coupling with chemistry"
        stop
      endif
      if (Interstitial%trans_aero) Interstitial%nvdiff = Interstitial%nvdiff + Model%ntchm
      if (Model%ntke > 0) Interstitial%nvdiff = Interstitial%nvdiff + 1    !  adding tke to the list
    endif

    if (Model%ntke > 0) Interstitial%ntkev = Interstitial%nvdiff

    if (Model%ntiw > 0) then
        if (Model%ntclamt > 0 .and. Model%ntsigma > 0 .and. Model%ntomega > 0) then
           Interstitial%nn = Model%ntrac - 4
        elseif (Model%ntclamt > 0 .and. Model%ntsigma > 0 .and. Model%ntomega <= 0) then
           Interstitial%nn = Model%ntrac - 3
        elseif (Model%ntclamt > 0 .and. Model%ntsigma <= 0 .and. Model%ntomega > 0) then
           Interstitial%nn = Model%ntrac - 3
        elseif (Model%ntclamt > 0 .and. Model%ntsigma <= 0 .and. Model%ntomega <= 0) then
           Interstitial%nn = Model%ntrac - 2
        elseif (Model%ntclamt <= 0 .and. Model%ntsigma > 0 .and. Model%ntomega > 0) then
           Interstitial%nn = Model%ntrac - 3
        elseif (Model%ntclamt <= 0 .and. Model%ntsigma > 0 .and. Model%ntomega <= 0) then
           Interstitial%nn = Model%ntrac - 2
        elseif (Model%ntclamt <= 0 .and. Model%ntsigma <= 0 .and. Model%ntomega > 0) then
           Interstitial%nn = Model%ntrac - 2
        else
           Interstitial%nn = Model%ntrac - 1
        endif
    elseif (Model%ntcw > 0) then
      Interstitial%nn = Model%ntrac
    else
      Interstitial%nn = Model%ntrac + 1
    endif

    if (Model%cscnv .or. Model%satmedmf .or. Model%trans_trac ) then
      Interstitial%otspt(:,:)   = .true.     ! otspt is used only for cscnv
      Interstitial%otspt(1:3,:) = .false.    ! this is for sp.hum, ice and liquid water
      Interstitial%otsptflag(:) = .true.
      tracers = 2
      do n=2,Model%ntrac
        ltest = ( n /= Model%ntcw  .and. n /= Model%ntiw  .and. n /= Model%ntclamt .and. &
                  n /= Model%ntrw  .and. n /= Model%ntsw  .and. n /= Model%ntrnc   .and. &
                  n /= Model%ntlnc .and. n /= Model%ntinc                          .and. &
                  n /= Model%ntsnc .and. n /= Model%ntgl  .and. n /= Model%ntgnc   .and. &
                  n /= Model%nthl  .and. n /= Model%nthnc .and. n /= Model%ntgv    .and. &
                  n /= Model%nthv  .and. n /= Model%ntccn .and. n /= Model%ntccna  .and. &
                  n /= Model%ntrz  .and. n /= Model%ntgz  .and. n /= Model%nthz    .and. &
                  n /= Model%ntsigma .and.  n /= Model%ntomega)
        Interstitial%otsptflag(n) = ltest
        if ( ltest ) then
          tracers = tracers + 1
          if (Model%ntke  == n ) then
            Interstitial%otspt(tracers+1,1) = .false.
            Interstitial%ntk = tracers
          endif
          if (Model%ntlnc == n .or. Model%ntinc == n .or. Model%ntrnc == n .or. Model%ntsnc == n .or. Model%ntgnc == n)    &
!           if (ntlnc == n .or. ntinc == n .or. ntrnc == n .or. ntsnc == n .or.&
!               ntrw  == n .or. ntsw  == n .or. ntgl  == n)                    &
                  Interstitial%otspt(tracers+1,1) = .false.
          if (Interstitial%trans_aero .and. Model%ntchs == n) Interstitial%itc = tracers
        endif
      enddo
      Interstitial%tracers_total = tracers - 2
    endif   ! end if_ras or cfscnv or samf
    if (.not. Model%satmedmf .and. .not. Model%trans_trac .and. &
        .not. Model%ras      .and. .not. Model%do_shoc) then
       Interstitial%nsamftrac = 0
    else
       Interstitial%nsamftrac = Interstitial%tracers_total
    endif
    Interstitial%ncstrac = Interstitial%tracers_total + 3

  end subroutine gfs_interstitial_setup_tracers

  subroutine gfs_interstitial_reset (Interstitial, Model)
    !
    implicit none
    !
    class(GFS_interstitial_type) :: Interstitial
    type(GFS_control_type), intent(in) :: Model
    !
    Interstitial%adjsfculw_land  = clear_val
    Interstitial%adjsfculw_ice   = clear_val
    Interstitial%adjsfculw_water = clear_val
    Interstitial%adjnirbmd       = clear_val
    Interstitial%adjnirbmu       = clear_val
    Interstitial%adjnirdfd       = clear_val
    Interstitial%adjnirdfu       = clear_val
    Interstitial%adjvisbmd       = clear_val
    Interstitial%adjvisbmu       = clear_val
    Interstitial%adjvisdfu       = clear_val
    Interstitial%adjvisdfd       = clear_val
    Interstitial%aerodp          = clear_val
    Interstitial%alb1d           = clear_val
    if (.not. Model%do_RRTMGP) then
      Interstitial%alpha         = clear_val
    end if
    Interstitial%bexp1d          = clear_val
    Interstitial%cd              = clear_val
    Interstitial%cd_ice          = Model%huge
    Interstitial%cd_land         = Model%huge
    Interstitial%cd_water        = Model%huge
    Interstitial%cdq             = clear_val
    Interstitial%cdq_ice         = Model%huge
    Interstitial%cdq_land        = Model%huge
    Interstitial%cdq_water       = Model%huge
    Interstitial%chh_ice         = Model%huge
    Interstitial%chh_land        = Model%huge
    Interstitial%chh_water       = Model%huge
    Interstitial%cldf            = clear_val
    Interstitial%cldsa           = clear_val
    Interstitial%cldtaulw        = clear_val
    Interstitial%cldtausw        = clear_val
    Interstitial%cld1d           = clear_val
    Interstitial%clouds          = clear_val
    Interstitial%clw             = clear_val
    Interstitial%clw(:,:,2)      = -999.9
    Interstitial%dclw            = clear_val
    Interstitial%clx             = clear_val
    Interstitial%cmm_ice         = Model%huge
    Interstitial%cmm_land        = Model%huge
    Interstitial%cmm_water       = Model%huge
    Interstitial%cnvc            = clear_val
    Interstitial%ctei_r          = clear_val
    Interstitial%ctei_rml        = clear_val
    Interstitial%cumabs          = clear_val
    Interstitial%dd_mf           = clear_val
    Interstitial%de_lgth         = clear_val
    Interstitial%del             = clear_val
    Interstitial%del_gz          = clear_val
    Interstitial%delr            = clear_val
    Interstitial%dlength         = clear_val
    Interstitial%dqdt            = clear_val
    Interstitial%dqsfc1          = clear_val
    Interstitial%drain           = clear_val
    Interstitial%dtdt            = clear_val
    Interstitial%dtsfc1          = clear_val
    Interstitial%dt_mf           = clear_val
    Interstitial%dtzm            = clear_val
    Interstitial%dudt            = clear_val
    Interstitial%dusfcg          = clear_val
    Interstitial%dusfc1          = clear_val
    Interstitial%dvdt            = clear_val
    Interstitial%dvsfcg          = clear_val
    Interstitial%dvsfc1          = clear_val
    Interstitial%dvdftra         = clear_val
    Interstitial%dzlyr           = clear_val
    Interstitial%elvmax          = clear_val
    Interstitial%ep1d            = clear_val
    Interstitial%ep1d_ice        = Model%huge
    Interstitial%ep1d_land       = Model%huge
    Interstitial%ep1d_water      = Model%huge
    Interstitial%evap_ice        = Model%huge
    Interstitial%evap_land       = Model%huge
    Interstitial%evap_water      = Model%huge
    Interstitial%evbs            = clear_val
    Interstitial%evcw            = clear_val
    Interstitial%pah             = clear_val
    Interstitial%ecan            = clear_val
    Interstitial%etran           = clear_val
    Interstitial%edir            = clear_val
    Interstitial%faerlw          = clear_val
    Interstitial%faersw          = clear_val
    Interstitial%ffhh_ice        = Model%huge
    Interstitial%ffhh_land       = Model%huge
    Interstitial%ffhh_water      = Model%huge
    Interstitial%fh2             = clear_val
    Interstitial%fh2_ice         = Model%huge
    Interstitial%fh2_land        = Model%huge
    Interstitial%fh2_water       = Model%huge
    Interstitial%flag_cice       = .false.
    Interstitial%flag_guess      = .false.
    Interstitial%flag_iter       = .true.
    Interstitial%flag_lakefreeze = .false.
    Interstitial%ffmm_ice        = Model%huge
    Interstitial%ffmm_land       = Model%huge
    Interstitial%ffmm_water      = Model%huge
    Interstitial%fm10            = clear_val
    Interstitial%fm10_ice        = Model%huge
    Interstitial%fm10_land       = Model%huge
    Interstitial%fm10_water      = Model%huge
    Interstitial%frland          = clear_val
    Interstitial%fscav           = clear_val
    Interstitial%fswtr           = clear_val
    Interstitial%gabsbdlw        = clear_val
    Interstitial%gabsbdlw_ice    = clear_val
    Interstitial%gabsbdlw_land   = clear_val
    Interstitial%gabsbdlw_water  = clear_val
    Interstitial%gamma           = clear_val
    Interstitial%gamq            = clear_val
    Interstitial%gamt            = clear_val
    Interstitial%gasvmr          = clear_val
    Interstitial%gflx            = clear_val
    Interstitial%gflx_ice        = clear_val
    Interstitial%gflx_land       = clear_val
    Interstitial%gflx_water      = clear_val
    Interstitial%gwdcu           = clear_val
    Interstitial%gwdcv           = clear_val
    Interstitial%zvfun           = clear_val
    Interstitial%hffac           = clear_val
    Interstitial%hflxq           = clear_val
    Interstitial%hflx_ice        = Model%huge
    Interstitial%hflx_land       = Model%huge
    Interstitial%hflx_water      = Model%huge
    Interstitial%htlwc           = clear_val
    Interstitial%htlw0           = clear_val
    Interstitial%htswc           = clear_val
    Interstitial%htsw0           = clear_val
    Interstitial%dry             = .false.
    Interstitial%idxday          = 0
    Interstitial%icy             = .false.
    Interstitial%lake            = .false.
    Interstitial%lndp_vgf        = clear_val
    Interstitial%ocean           = .false.
    Interstitial%islmsk          = 0
    Interstitial%islmsk_cice     = 0
    Interstitial%wet             = .false.
    Interstitial%kb              = 0
    Interstitial%kbot            = Model%levs
    Interstitial%kcnv            = 0
    Interstitial%kd              = 0
    Interstitial%kinver          = Model%levs
    Interstitial%kpbl            = 0
    Interstitial%kt              = 0
    Interstitial%ktop            = 1
    Interstitial%mbota           = 0
    Interstitial%mtopa           = 0
    Interstitial%nday            = 0
    Interstitial%oa4             = clear_val
    Interstitial%oc              = clear_val
    Interstitial%olyr            = clear_val
    Interstitial%plvl            = clear_val
    Interstitial%plyr            = clear_val
    Interstitial%prnum           = clear_val
    Interstitial%qlyr            = clear_val
    Interstitial%prcpmp          = clear_val
    Interstitial%qss_ice         = Model%huge
    Interstitial%qss_land        = Model%huge
    Interstitial%qss_water       = Model%huge
    Interstitial%raddt           = clear_val
    Interstitial%raincd          = clear_val
    Interstitial%raincs          = clear_val
    Interstitial%rainmcadj       = clear_val
    Interstitial%rainp           = clear_val
    Interstitial%rb              = clear_val
    Interstitial%rb_ice          = Model%huge
    Interstitial%rb_land         = Model%huge
    Interstitial%rb_water        = Model%huge
    Interstitial%rhc             = clear_val
    Interstitial%runoff          = clear_val
    Interstitial%save_q          = clear_val
    Interstitial%save_t          = clear_val
    Interstitial%save_tcp        = clear_val
    Interstitial%save_u          = clear_val
    Interstitial%save_v          = clear_val
    Interstitial%sbsno           = clear_val
    Interstitial%scmpsw%uvbfc    = clear_val
    Interstitial%scmpsw%uvbf0    = clear_val
    Interstitial%scmpsw%nirbm    = clear_val
    Interstitial%scmpsw%nirdf    = clear_val
    Interstitial%scmpsw%visbm    = clear_val
    Interstitial%scmpsw%visdf    = clear_val
    Interstitial%sfcalb          = clear_val
    Interstitial%sigma           = clear_val
    Interstitial%sigmaf          = clear_val
    Interstitial%sigmafrac       = clear_val
    Interstitial%sigmatot        = clear_val
    Interstitial%snowc           = clear_val
    Interstitial%snohf           = clear_val
    Interstitial%snowmt          = clear_val
    Interstitial%stress          = clear_val
    Interstitial%stress_ice      = Model%huge
    Interstitial%stress_land     = Model%huge
    Interstitial%stress_water    = Model%huge
    Interstitial%ten_q           = clear_val
    Interstitial%ten_t           = clear_val
    Interstitial%ten_u           = clear_val
    Interstitial%ten_v           = clear_val
    Interstitial%theta           = clear_val
    Interstitial%tkeh            = 0
    Interstitial%tlvl            = clear_val
    Interstitial%tlyr            = clear_val
    Interstitial%tprcp_ice       = Model%huge
    Interstitial%tprcp_land      = Model%huge
    Interstitial%tprcp_water     = Model%huge
    Interstitial%trans           = clear_val
    Interstitial%tseal           = clear_val
    Interstitial%tsfa            = clear_val
    Interstitial%tsfc_water      = Model%huge
    Interstitial%tsfg            = clear_val
    Interstitial%tsurf_ice       = Model%huge
    Interstitial%tsurf_land      = Model%huge
    Interstitial%tsurf_water     = Model%huge
    Interstitial%ud_mf           = clear_val
    Interstitial%uustar_ice      = Model%huge
    Interstitial%uustar_land     = Model%huge
    Interstitial%uustar_water    = Model%huge
    Interstitial%vdftra          = clear_val
    Interstitial%vegf1d          = clear_val
    Interstitial%wcbmax          = clear_val
    Interstitial%wind            = Model%huge
    Interstitial%work1           = clear_val
    Interstitial%work2           = clear_val
    Interstitial%work3           = clear_val
    Interstitial%xcosz           = clear_val
    Interstitial%xlai1d          = clear_val
    Interstitial%xmu             = clear_val
    Interstitial%z01d            = clear_val
    Interstitial%zt1d            = clear_val
    Interstitial%ztmax_ice       = clear_val
    Interstitial%ztmax_land      = clear_val
    Interstitial%ztmax_water     = clear_val

    ! RRTMGP
    if (Model%do_RRTMGP) then
       Interstitial%tracer                      = clear_val
       Interstitial%tv_lay                      = clear_val
       Interstitial%relhum                      = clear_val
       Interstitial%qs_lay                      = clear_val
       Interstitial%q_lay                       = clear_val
       Interstitial%deltaZ                      = clear_val
       Interstitial%deltaZc                     = clear_val
       Interstitial%deltaP                      = clear_val
       Interstitial%p_lev                       = clear_val
       Interstitial%p_lay                       = clear_val
       Interstitial%t_lev                       = clear_val
       Interstitial%t_lay                       = clear_val
       Interstitial%cloud_overlap_param         = clear_val
       Interstitial%precip_overlap_param        = clear_val
       Interstitial%fluxlwUP_allsky             = clear_val
       Interstitial%fluxlwDOWN_allsky           = clear_val
       Interstitial%fluxlwUP_clrsky             = clear_val
       Interstitial%fluxlwDOWN_clrsky           = clear_val
       Interstitial%fluxswUP_allsky             = clear_val
       Interstitial%fluxswDOWN_allsky           = clear_val
       Interstitial%fluxswUP_clrsky             = clear_val
       Interstitial%fluxswDOWN_clrsky           = clear_val
       Interstitial%aerosolslw                  = clear_val
       Interstitial%aerosolssw                  = clear_val
       Interstitial%precip_frac                 = clear_val
       Interstitial%cld_cnv_frac                = clear_val
       Interstitial%cnv_cloud_overlap_param     = clear_val
       Interstitial%cld_cnv_lwp                 = clear_val
       Interstitial%cld_cnv_reliq               = clear_val
       Interstitial%cld_cnv_iwp                 = clear_val
       Interstitial%cld_cnv_reice               = clear_val
       Interstitial%cld_pbl_lwp                 = clear_val
       Interstitial%cld_pbl_reliq               = clear_val
       Interstitial%cld_pbl_iwp                 = clear_val
       Interstitial%cld_pbl_reice               = clear_val
       Interstitial%flxprf_lw%upfxc             = clear_val
       Interstitial%flxprf_lw%dnfxc             = clear_val
       Interstitial%flxprf_lw%upfx0             = clear_val
       Interstitial%flxprf_lw%dnfx0             = clear_val
       Interstitial%flxprf_sw%upfxc             = clear_val
       Interstitial%flxprf_sw%dnfxc             = clear_val
       Interstitial%flxprf_sw%upfx0             = clear_val
       Interstitial%flxprf_sw%dnfx0             = clear_val
       Interstitial%sfc_emiss_byband            = clear_val
       Interstitial%sec_diff_byband             = clear_val
       Interstitial%sfc_alb_nir_dir             = clear_val
       Interstitial%sfc_alb_nir_dif             = clear_val
       Interstitial%sfc_alb_uvvis_dir           = clear_val
       Interstitial%sfc_alb_uvvis_dif           = clear_val
       Interstitial%toa_src_sw                  = clear_val
       Interstitial%toa_src_lw                  = clear_val
       Interstitial%vmr_o2                      = clear_val
       Interstitial%vmr_h2o                     = clear_val
       Interstitial%vmr_o3                      = clear_val
       Interstitial%vmr_ch4                     = clear_val
       Interstitial%vmr_n2o                     = clear_val
       Interstitial%vmr_co2                     = clear_val
    end if

    ! UGWP common
    Interstitial%tau_mtb         = clear_val
    Interstitial%tau_ogw         = clear_val
    Interstitial%tau_tofd        = clear_val
    Interstitial%tau_ngw         = clear_val
    Interstitial%tau_oss         = clear_val
    Interstitial%dudt_mtb        = clear_val
    Interstitial%dudt_tms        = clear_val
    Interstitial%zmtb            = clear_val
    Interstitial%zlwb            = clear_val
    Interstitial%zogw            = clear_val
    Interstitial%zngw            = clear_val

! CIRES UGWP v1
    if (Model%ldiag_ugwp .or. Model%do_ugwp_v0 .or. Model%do_ugwp_v0_nst_only &
        .or. Model%do_ugwp_v1) then
      Interstitial%dudt_ngw      = clear_val
      Interstitial%dvdt_ngw      = clear_val
      Interstitial%dtdt_ngw      = clear_val
      Interstitial%kdis_ngw      = clear_val
    end if

!-- GSL drag suite
    if (Model%gwd_opt==3 .or. Model%gwd_opt==33 .or. &
        Model%gwd_opt==2 .or. Model%gwd_opt==22 ) then
       Interstitial%varss        = clear_val
       Interstitial%ocss         = clear_val
       Interstitial%oa4ss        = clear_val
       Interstitial%clxss        = clear_val
    end if
!
    ! Allocate arrays that are conditional on physics choices
    if (Model%imp_physics == Model%imp_physics_gfdl .or. Model%imp_physics == Model%imp_physics_thompson &
         .or. Model%imp_physics == Model%imp_physics_tempo .or. Model%imp_physics == Model%imp_physics_nssl &
        ) then
       Interstitial%graupelmp    = clear_val
       Interstitial%icemp        = clear_val
       Interstitial%rainmp       = clear_val
       Interstitial%snowmp       = clear_val
    else if (Model%imp_physics == Model%imp_physics_mg) then
       Interstitial%ncgl         = clear_val
       Interstitial%ncpr         = clear_val
       Interstitial%ncps         = clear_val
       Interstitial%qgl          = clear_val
       Interstitial%qrn          = clear_val
       Interstitial%qsnw         = clear_val
       Interstitial%qlcn         = clear_val
       Interstitial%qicn         = clear_val
       Interstitial%w_upi        = clear_val
       Interstitial%cf_upi       = clear_val
       Interstitial%cnv_mfd      = clear_val
       Interstitial%cnv_dqldt    = clear_val
       Interstitial%clcn         = clear_val
       Interstitial%cnv_fice     = clear_val
       Interstitial%cnv_ndrop    = clear_val
       Interstitial%cnv_nice     = clear_val
    end if
    if (Model%lsm == Model%lsm_noahmp) then
       Interstitial%t2mmp        = clear_val
       Interstitial%q2mp         = clear_val
    end if

    ! Set flag for resetting maximum hourly output fields
    Interstitial%max_hourly_reset = mod(Model%kdt-1, nint(Model%avg_max_length/Model%dtp)) == 0
    ! Use same logic in UFS to reset Thompson extended diagnostics
    Interstitial%ext_diag_thompson_reset = Interstitial%max_hourly_reset

    ! Frequency flag for computing the full radar reflectivity (water coated ice) 
    if (Model%nsfullradar_diag<0) then
      Interstitial%fullradar_diag = .true.
    else
      Interstitial%fullradar_diag = (Model%kdt == 1 .or. mod(Model%kdt, nint(Model%nsfullradar_diag/Model%dtp)) == 0) 
    end if
  end subroutine gfs_interstitial_reset

!-----------------------
! GFDL_interstitial_type
!-----------------------

  subroutine gfdl_interstitial_create (Interstitial, is, ie, isd, ied, js, je, jsd, jed, npz, ng, &
                                       dt_atmos, p_split, k_split, zvir, p_ref, ak, bk,           &
                                       do_ql, do_qi, do_qr, do_qs, do_qg, do_qa,                  &
                                       kappa, hydrostatic, do_sat_adj,                            &
                                       delp, delz, area, peln, phis, pkz, pt,                     &
                                       qvi, qv, ql, qi, qr, qs, qg, qc, q_con,                    &
                                       nthreads, nwat, ngas, rilist, cpilist, mpirank, mpiroot)
    !
    implicit none
    !
    class(GFDL_interstitial_type) :: Interstitial
    integer, intent(in) :: is
    integer, intent(in) :: ie
    integer, intent(in) :: isd
    integer, intent(in) :: ied
    integer, intent(in) :: js
    integer, intent(in) :: je
    integer, intent(in) :: jsd
    integer, intent(in) :: jed
    integer, intent(in) :: npz
    integer, intent(in) :: ng
    real(kind_dyn),    intent(in) :: dt_atmos
    integer, intent(in) :: p_split
    integer, intent(in) :: k_split
    real(kind_dyn),    intent(in) :: zvir
    real(kind_dyn),    intent(in) :: p_ref
    real(kind_dyn),    intent(in) :: ak(:)
    real(kind_dyn),    intent(in) :: bk(:)
    logical, intent(in) :: do_ql
    logical, intent(in) :: do_qi
    logical, intent(in) :: do_qr
    logical, intent(in) :: do_qs
    logical, intent(in) :: do_qg
    logical, intent(in) :: do_qa
    real(kind_dyn),    intent(in) :: kappa
    logical, intent(in) :: hydrostatic
    logical, intent(in) :: do_sat_adj
    real(kind_dyn),  target, intent(in) :: delp(:,:,:)
    real(kind_dyn),  target, intent(in) :: delz(:,:,:)
    real(kind_grid), target, intent(in) :: area(:,:)
    real(kind_dyn),  target, intent(in) :: peln(:,:,:)
    real(kind_dyn),  target, intent(in) :: phis(:,:)
    real(kind_dyn),  target, intent(in) :: pkz(:,:,:)
    real(kind_dyn),  target, intent(in) :: pt(:,:,:)
    real(kind_dyn),  target, intent(in) :: qvi(:,:,:,:)
    real(kind_dyn),  target, intent(in) :: qv(:,:,:)
    real(kind_dyn),  target, intent(in) :: ql(:,:,:)
    real(kind_dyn),  target, intent(in) :: qi(:,:,:)
    real(kind_dyn),  target, intent(in) :: qr(:,:,:)
    real(kind_dyn),  target, intent(in) :: qs(:,:,:)
    real(kind_dyn),  target, intent(in) :: qg(:,:,:)
    real(kind_dyn),  target, intent(in) :: qc(:,:,:)
    real(kind_dyn),  target, intent(in) :: q_con(:,:,:)
    integer, intent(in) :: nthreads
    ! For multi-gases physics
    integer,        intent(in)           :: nwat
    integer,        intent(in), optional :: ngas
    real(kind_dyn), intent(in), optional :: rilist(0:)
    real(kind_dyn), intent(in), optional :: cpilist(0:)
    integer,        intent(in)           :: mpirank
    integer,        intent(in)           :: mpiroot
    !
    integer :: isc1, jsc1, iec1, jec1
    integer :: isc2, jsc2, iec2, jec2
    !
    isc1 = lbound(delp, dim=1)
    jsc1 = lbound(delp, dim=2)
    iec1 = ubound(delp, dim=1)
    jec1 = ubound(delp, dim=2)
    isc2 = lbound(delz, dim=1)
    jsc2 = lbound(delz, dim=2)
    iec2 = ubound(delz, dim=1)
    jec2 = ubound(delz, dim=2)
    !
#ifdef MOIST_CAPPA
    Interstitial%npzcappa = npz
    allocate (Interstitial%cappa  (isd:ied, jsd:jed, 1:npz) )
#else
    Interstitial%npzcappa = 1
    allocate (Interstitial%cappa  (isd:isd, jsd:jsd, 1)     )
#endif
    allocate (Interstitial%dtdt   (is:ie, js:je, 1:npz)     )
    allocate (Interstitial%pfull  (1:npz)                   )
    allocate (Interstitial%te0_2d (is:ie, js:je)            )
    allocate (Interstitial%te0    (isd:ied, jsd:jed, 1:npz) )
    !
    ! Initialize variables to default values
#ifdef SW_DYNAMICS
    Interstitial%akap = 1.
#else
    Interstitial%akap = kappa
#endif
    Interstitial%bdt       = dt_atmos/real(abs(p_split))
    Interstitial%do_qa     = do_qa
    Interstitial%mdt       = Interstitial%bdt/real(k_split)
    !
    ! Flag for hydrostatic/non-hydrostatic physics
    Interstitial%hydrostatic = hydrostatic
    if (hydrostatic) then
       Interstitial%npzdelz = 1
    else
       Interstitial%npzdelz = npz
    end if
    !
    Interstitial%zvir      = zvir
    !
    Interstitial%do_sat_adj =  do_sat_adj
    Interstitial%is         =  is
    Interstitial%ie         =  ie
    Interstitial%isd        =  isd
    Interstitial%ied        =  ied
    Interstitial%isc1       =  isc1
    Interstitial%iec1       =  iec1
    Interstitial%isc2       =  isc2
    Interstitial%iec2       =  iec2
    Interstitial%js         =  js
    Interstitial%je         =  je
    Interstitial%jsd        =  jsd
    Interstitial%jed        =  jed
    Interstitial%jsc1       =  jsc1
    Interstitial%jec1       =  jec1
    Interstitial%jsc2       =  jsc2
    Interstitial%jec2       =  jec2
    Interstitial%ng         =  ng
    Interstitial%npz        =  npz
    Interstitial%npzp1      =  npz+1
    !
    ! Set up links from GFDL_interstitial DDT to ATM DDT
    Interstitial%delp       => delp
    Interstitial%delz       => delz
    Interstitial%area       => area
    Interstitial%peln       => peln
    Interstitial%phis       => phis
    Interstitial%pkz        => pkz
    Interstitial%pt         => pt
    Interstitial%qvi        => qvi
    Interstitial%qv         => qv
    if (do_ql) Interstitial%ql => ql
    if (do_qi) Interstitial%qi => qi
    if (do_qr) Interstitial%qr => qr
    if (do_qs) Interstitial%qs => qs
    if (do_qg) Interstitial%qg => qg
    if (do_qa) Interstitial%qc => qc
    !
#ifdef USE_COND
    Interstitial%npzq_con = npz
#else
    Interstitial%npzq_con = 1
#endif
    Interstitial%q_con      => q_con
    !
    ! Number of OpenMP threads available for schemes
    Interstitial%nthreads   = nthreads
    !
    ! For multi-gases physics
    Interstitial%nwat  = nwat
    ! If ngas, rilist and cpilist are present, then
    ! multi-gases physics are used. If not, set ngas=0
    ! (safe value), allocate rilist/cpilist and set to zero
    if(present(ngas)) then
      Interstitial%ngas  = ngas
    else
      Interstitial%ngas  = 0
    end if
    allocate(Interstitial%rilist(0:Interstitial%ngas))
    allocate(Interstitial%cpilist(0:Interstitial%ngas))
    if (present(rilist)) then
      Interstitial%rilist  = rilist(0:Interstitial%ngas)
      Interstitial%cpilist = cpilist(0:Interstitial%ngas)
    else
      Interstitial%rilist  = 0.0
      Interstitial%cpilist = 0.0
    end if
    !
    Interstitial%mpirank = mpirank
    Interstitial%mpiroot = mpiroot
    !
    ! Calculate vertical pressure levels
    call gfdl_interstitital_calculate_pressure_levels(Interstitial, npz, p_ref, ak, bk)
    !
    ! Reset all other variables
    call Interstitial%reset()
    !
  end subroutine gfdl_interstitial_create

  subroutine gfdl_interstitital_calculate_pressure_levels(Interstitial, npz, p_ref, ak, bk)

     implicit none

     class(GFDL_interstitial_type)       :: Interstitial
     integer, intent(in) :: npz
     real(kind_dyn), intent(in) :: p_ref
     real(kind_dyn), intent(in) :: ak(:)
     real(kind_dyn), intent(in) :: bk(:)

     real(kind_dyn) :: ph1
     real(kind_dyn) :: ph2
     integer :: k

#ifdef SW_DYNAMICS
      Interstitial%pfull(1) = 0.5*p_ref
#else
      do k=1,npz
         ph1 = ak(k  ) + bk(k  )*p_ref
         ph2 = ak(k+1) + bk(k+1)*p_ref
         Interstitial%pfull(k) = (ph2 - ph1) / log(ph2/ph1)
      enddo
#endif
      ! DH* This is copied from fv_mapz.F90, does it work with SW_DYNAMICS?
      do k=1,npz
         Interstitial%kmp = k
         if ( Interstitial%pfull(k) > 10.E2 ) exit
      enddo
  end subroutine gfdl_interstitital_calculate_pressure_levels

  subroutine gfdl_interstitial_reset (Interstitial)
    !
    implicit none
    !
    class(GFDL_interstitial_type) :: Interstitial
    !
    Interstitial%cappa         = 0.0
    Interstitial%dtdt          = 0.0
    Interstitial%fast_mp_consv = .false.
    Interstitial%last_step     = .false.
    Interstitial%out_dt        = .false.
    Interstitial%te0_2d        = 0.0
    Interstitial%te0           = 0.0
    !
  end subroutine gfdl_interstitial_reset

  subroutine gfdl_interstitial_print(Interstitial)
    !
    implicit none
    !
    class(GFDL_interstitial_type) :: Interstitial
    !
    ! Print static variables
    write (0,'(a)') 'Interstitial_print'
    write (0,*) 'Interstitial_print: values that do not change'
    write (0,*) 'Interstitial%akap              = ', Interstitial%akap
    write (0,*) 'Interstitial%bdt               = ', Interstitial%bdt
    write (0,*) 'Interstitial%kmp               = ', Interstitial%kmp
    write (0,*) 'Interstitial%mdt               = ', Interstitial%mdt
    write (0,*) 'sum(Interstitial%pfull)        = ', sum(Interstitial%pfull)
    write (0,*) 'Interstitial%zvir              = ', Interstitial%zvir
    write (0,*) 'Interstitial%do_qa             = ', Interstitial%do_qa
    write (0,*) 'Interstitial%do_sat_adj        = ', Interstitial%do_sat_adj
    write (0,*) 'Interstitial%is                = ', Interstitial%is
    write (0,*) 'Interstitial%ie                = ', Interstitial%ie
    write (0,*) 'Interstitial%isd               = ', Interstitial%isd
    write (0,*) 'Interstitial%ied               = ', Interstitial%ied
    write (0,*) 'Interstitial%js                = ', Interstitial%js
    write (0,*) 'Interstitial%je                = ', Interstitial%je
    write (0,*) 'Interstitial%jsd               = ', Interstitial%jsd
    write (0,*) 'Interstitial%jed               = ', Interstitial%jed
    write (0,*) 'sum(Interstitial%area)         = ', sum(Interstitial%area)
    write (0,*) 'Interstitial%ng                = ', Interstitial%ng
    write (0,*) 'Interstitial%npz               = ', Interstitial%npz
    write (0,*) 'Interstitial%npzp1             = ', Interstitial%npzp1
    write (0,*) 'Interstitial%npzdelz           = ', Interstitial%npzdelz
    write (0,*) 'Interstitial%npzq_con          = ', Interstitial%npzq_con
    write (0,*) 'Interstitial%npzcappa          = ', Interstitial%npzcappa
    ! Print all other variables
    write (0,*) 'Interstitial_print: values that change'
    write (0,*) 'sum(Interstitial%cappa)        = ', sum(Interstitial%cappa)
    write (0,*) 'sum(Interstitial%dtdt)         = ', sum(Interstitial%dtdt)
    write (0,*) 'Interstitial%fast_mp_consv     = ', Interstitial%fast_mp_consv
    write (0,*) 'Interstitial%last_step         = ', Interstitial%last_step
    write (0,*) 'Interstitial%out_dt            = ', Interstitial%out_dt
    write (0,*) 'sum(Interstitial%te0_2d)       = ', sum(Interstitial%te0_2d)
    write (0,*) 'sum(Interstitial%te0)          = ', sum(Interstitial%te0)
    write (0,*) 'sum(Interstitial%delp)         = ', Interstitial%delp
    write (0,*) 'sum(Interstitial%delz)         = ', Interstitial%delz
    write (0,*) 'sum(Interstitial%peln)         = ', Interstitial%peln
    write (0,*) 'sum(Interstitial%phis)         = ', Interstitial%phis
    write (0,*) 'sum(Interstitial%pkz)          = ', Interstitial%pkz
    write (0,*) 'sum(Interstitial%pt)           = ', Interstitial%pt
    write (0,*) 'sum(Interstitial%qvi)          = ', Interstitial%qvi
    write (0,*) 'sum(Interstitial%qv)           = ', Interstitial%qv
    if (associated(Interstitial%ql)) write (0,*) 'sum(Interstitial%ql)           = ', Interstitial%ql
    if (associated(Interstitial%qi)) write (0,*) 'sum(Interstitial%qi)           = ', Interstitial%qi
    if (associated(Interstitial%qr)) write (0,*) 'sum(Interstitial%qr)           = ', Interstitial%qr
    if (associated(Interstitial%qs)) write (0,*) 'sum(Interstitial%qs)           = ', Interstitial%qs
    if (associated(Interstitial%qg)) write (0,*) 'sum(Interstitial%qg)           = ', Interstitial%qg
    if (associated(Interstitial%qc)) write (0,*) 'sum(Interstitial%qc)           = ', Interstitial%qc
    write (0,*) 'sum(Interstitial%q_con)        = ', Interstitial%q_con
    write (0,*) 'Interstitial%hydrostatic       = ', Interstitial%hydrostatic
    write (0,*) 'Interstitial%nwat              = ', Interstitial%nwat
    write (0,*) 'Interstitial%ngas              = ', Interstitial%ngas
    write (0,*) 'Interstitial%rilist            = ', Interstitial%rilist
    write (0,*) 'Interstitial%cpilist           = ', Interstitial%cpilist
    write (0,*) 'Interstitial%mpirank           = ', Interstitial%mpirank
    write (0,*) 'Interstitial%mpiroot           = ', Interstitial%mpiroot
    write (0,*) 'Interstitial%nthreads          = ', Interstitial%nthreads
    write (0,*) 'Interstitial_print: end'
    !
  end subroutine gfdl_interstitial_print

end module CCPP_typedefs

