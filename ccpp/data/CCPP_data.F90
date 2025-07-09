module CCPP_data

!! \section arg_table_CCPP_data Argument Table
!! \htmlinclude CCPP_data.html
!!

    use ccpp_types,    only: ccpp_t
    use CCPP_typedefs, only: GFS_interstitial_type,  &
                             GFDL_interstitial_type
    use GFS_typedefs,  only: GFS_control_type,       &
                             GFS_statein_type,       &
                             GFS_stateout_type,      &
                             GFS_grid_type,          &
                             GFS_tbd_type,           &
                             GFS_cldprop_type,       &
                             GFS_sfcprop_type,       &
                             GFS_radtend_type,       &
                             GFS_coupling_type,      &
                             GFS_diag_type

    implicit none

    private

    public cdata_tile,             &
           cdata_domain,           &
           cdata_block,            &
           ccpp_suite,             &
           GFDL_interstitial,      &
           GFS_control,            &
           GFS_statein,            &
           GFS_stateout,           &
           GFS_grid,               &
           GFS_tbd,                &
           GFS_cldprop,            &
           GFS_sfcprop,            &
           GFS_radtend,            &
           GFS_coupling,           &
           GFS_intdiag,            &
           GFS_interstitial

    !-------------------------------------------------------!
    !  GFS data containers;                                 !
    !  GFS_Interstitial has dimension nthreads              !
    !-------------------------------------------------------!
    type(GFS_control_type),                                    save, target :: GFS_control
    type(GFS_statein_type),                                    save, target :: GFS_statein
    type(GFS_stateout_type),                                   save, target :: GFS_stateout
    type(GFS_grid_type),                                       save, target :: GFS_grid
    type(GFS_tbd_type),                                        save, target :: GFS_tbd
    type(GFS_cldprop_type),                                    save, target :: GFS_cldprop
    type(GFS_sfcprop_type),                                    save, target :: GFS_sfcprop
    type(GFS_radtend_type),                                    save, target :: GFS_radtend
    type(GFS_coupling_type),                                   save, target :: GFS_coupling
    type(GFS_diag_type),                                       save, target :: GFS_intdiag
    type(GFS_interstitial_type),  dimension(:),   allocatable, save, target :: GFS_interstitial

    !------------------------------------------------------!
    !  CCPP data containers for dynamics (fast physics)    !
    !------------------------------------------------------!
    type(GFDL_interstitial_type),                              save, target :: GFDL_interstitial

    !------------------------------------------------------!
    !  CCPP containers for the six tiles used in dynamics, !
    !  for the entire domain and for the individual blocks !
    !  with dimensions nblocks and nthreads                !
    !------------------------------------------------------!
    type(ccpp_t),                                              save, target :: cdata_tile
    type(ccpp_t),                                              save, target :: cdata_domain
    type(ccpp_t),                 dimension(:,:), allocatable, save, target :: cdata_block

    !------------------------------------------------------!
    !  CCPP suite name                                     !
    !------------------------------------------------------!
    character(len=256)      :: ccpp_suite='undefined'

end module CCPP_data
