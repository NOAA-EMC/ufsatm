module GFS_restart

  use machine,          only: kind_phys
  use GFS_typedefs,     only: GFS_control_type,  GFS_statein_type,  &
                              GFS_stateout_type, GFS_sfcprop_type,  &
                              GFS_coupling_type, GFS_grid_type,     &
                              GFS_tbd_type,      GFS_cldprop_type,  &
                              GFS_radtend_type,  GFS_diag_type,     &
                              GFS_init_type
  use GFS_diagnostics,  only: GFS_externaldiag_type

  type var_subtype
     real(kind=kind_phys), dimension(:),   pointer :: var2 => null()
     real(kind=kind_phys), dimension(:,:), pointer :: var3 => null()
  end type var_subtype

  type GFS_restart_type
     integer           :: axes  !< Rank of data (2D or 3D).
     logical           :: diag  !< True for diagnostic field.
     logical           :: reset !< If true, zero out diagnostic field.
     character(len=32) :: name  !< variable name as it will appear in the restart file.
     type(var_subtype) :: data  !< Holds pointers to contiguous data.
  end type GFS_restart_type

  public GFS_restart_type, GFS_restart_populate

  CONTAINS
!*******************************************************************************************

!---------------------
! GFS_restart_populate
!---------------------
  subroutine GFS_restart_populate (Restart, Model, Statein, Stateout, Sfcprop,     &
                                   Coupling, Grid, Tbd, Cldprop, Radtend, IntDiag, &
                                   Init_parm, ExtDiag)
!----------------------------------------------------------------------------------------!
!   RESTART_METADATA                                                                     !
!     Restart%axes           [int*4  ]  Number of axes (rank) of variable                !
!     Restart%diag           [logical]  Flag to indicate diagnostic variable             !
!     Restart%reset          [logical]  Flag to indicate diagnostics need to be reset    !
!     Restart%name           [char=32]  Variable name in restart file                    !
!     Restart%data%var2(:)   [real*8 ]  pointer to 2D data (im)                          !
!     Restart%data%var3(:,:) [real*8 ]  pointer to 3D data (im,levs)                     !
!----------------------------------------------------------------------------------------!
    type(GFS_restart_type),     intent(inout), allocatable :: Restart(:)
    type(GFS_control_type),     intent(in)    :: Model
    type(GFS_statein_type),     intent(in)    :: Statein
    type(GFS_stateout_type),    intent(in)    :: Stateout
    type(GFS_sfcprop_type),     intent(in)    :: Sfcprop
    type(GFS_coupling_type),    intent(in)    :: Coupling
    type(GFS_grid_type),        intent(in)    :: Grid
    type(GFS_tbd_type),         intent(in)    :: Tbd
    type(GFS_cldprop_type),     intent(in)    :: Cldprop
    type(GFS_radtend_type),     intent(in)    :: Radtend
    type(GFS_diag_type),        intent(in)    :: IntDiag
    type(GFS_init_type),        intent(in)    :: Init_parm
    type(GFS_externaldiag_type),intent(in)    :: ExtDiag(:)

    !--- local variables
    integer :: idx, ndiag_rst
    integer :: ndiag_idx(20), itime
    integer ::  num, offset
    character(len=2) :: c2 = ''
    logical :: surface_layer_saves_rainprev
    integer :: num2d, num3d

    !--- check if continuous accumulated total precip and total cnvc precip are
    !    requested in output. If so, store location into Diagnsotic type.
    ndiag_rst = 0
    ndiag_idx(1:20) = 0
    do idx=1, size(ExtDiag)
      if( ExtDiag(idx)%id > 0) then
        if( trim(ExtDiag(idx)%name) == 'totprcp_ave') then
          ndiag_rst = ndiag_rst +1
          ndiag_idx(ndiag_rst) = idx
        else if( trim(ExtDiag(idx)%name) == 'cnvprcp_ave') then
          ndiag_rst = ndiag_rst +1
          ndiag_idx(ndiag_rst) = idx
        else if( trim(ExtDiag(idx)%name) == 'totice_ave') then
          ndiag_rst = ndiag_rst +1
          ndiag_idx(ndiag_rst) = idx
        else if( trim(ExtDiag(idx)%name) == 'totsnw_ave') then
          ndiag_rst = ndiag_rst +1
          ndiag_idx(ndiag_rst) = idx
        else if( trim(ExtDiag(idx)%name) == 'totgrp_ave') then
          ndiag_rst = ndiag_rst +1
          ndiag_idx(ndiag_rst) = idx
        else if( trim(ExtDiag(idx)%name) == 'tsnowp') then
          ndiag_rst = ndiag_rst +1
          ndiag_idx(ndiag_rst) = idx
        else if( trim(ExtDiag(idx)%name) == 'frozr') then
          ndiag_rst = ndiag_rst +1
          ndiag_idx(ndiag_rst) = idx
        else if( trim(ExtDiag(idx)%name) == 'frzr') then
          ndiag_rst = ndiag_rst +1
          ndiag_idx(ndiag_rst) = idx
        endif
      endif
    enddo
  
    ! Number of required 2D restart variables.
    num2d = 3 + Model%ntot2d + Model%nctp + ndiag_rst

    ! The CLM Lake Model needs raincprev and rainncprv, which some
    ! surface layer schemes save, and some don't. If the surface layer
    ! scheme does not save that variable, then it'll be saved
    ! separately for clm_lake.
    surface_layer_saves_rainprev = .false.

    ! Do we have any 2D restart varaibles dependent on physics scheme?
    ! GF
    if (Model%imfdeepcnv == Model%imfdeepcnv_gf) then
       num2d = num2d + 3
    endif
    ! Unified convection
    if (Model%imfdeepcnv == Model%imfdeepcnv_c3) then
       num2d = num2d + 3
    endif
    ! CA
    if (Model%imfdeepcnv == 2 .and. Model%do_ca) then
       num2d = num2d + 1
    endif
    ! NoahMP
    if (Model%lsm == Model%lsm_noahmp) then
       num2d = num2d + 10
       surface_layer_saves_rainprev = .true.
    endif
    ! RUC 
    if (Model%lsm == Model%lsm_ruc) then
       num2d = num2d + 5
       surface_layer_saves_rainprev = .true.
    endif
    ! MYNN SFC
    if (Model%do_mynnsfclay) then
       num2d = num2d + 13
    endif
    ! Save rain prev for lake if surface layer doesn't.
    if (Model%lkm>0 .and. Model%iopt_lake==Model%iopt_lake_clm .and. &
         .not.surface_layer_saves_rainprev) then
       num2d = num2d + 2
    endif
    ! Thompson aerosol-aware
    if ((Model%imp_physics == Model%imp_physics_thompson .or. &
         Model%imp_physics == Model%imp_physics_tempo) .and. (Model%ltaerosol)) then
       num2d = num2d + 2
    endif
    if (Model%do_cap_suppress .and. Model%num_dfi_radar>0) then
       num2d = num2d + Model%num_dfi_radar
    endif
    if (Model%rrfs_sd) then
       num2d = num2d + 6
    endif

    ! Number of required 3D restart variables.
    num3d = Model%ntot3d
    
    ! Do we have any 3D restart varaibles dependent on physics scheme?
    if (Model%num_dfi_radar>0) then
       num3d = num3d + Model%num_dfi_radar
    endif
    if(Model%lrefres) then
       num3d = Model%ntot3d+1
    endif
    ! General Convection
    if (Model%imfdeepcnv == Model%imfdeepcnv_gf) then
       num3d = num3d + 1
    endif
    ! GF
    if (Model%imfdeepcnv == 3) then
       num3d = num3d + 3
    endif
    ! Unified convection
    if (Model%imfdeepcnv == 5) then
       num3d = num3d + 4
    endif
    ! MYNN PBL
    if (Model%do_mynnedmf) then
       num3d = num3d + 9
    endif
    if (Model%rrfs_sd) then
       num3d = num3d + 4
    endif
    !Prognostic area fraction
    if (Model%progsigma) then
       num3d = num3d + 2
    endif

    if (Model%num_dfi_radar > 0) then
      do itime=1,Model%dfi_radar_max_intervals
        if(Model%ix_dfi_radar(itime)>0) then
           num3d = num3d + 1
        endif
      enddo
    endif

    !--- Allocate Restart data type.
    allocate (Restart(num2d+num3d))
    Restart(:)%diag  = .false.
    Restart(:)%reset = .false.
    idx = 0

    !--- Cldprop variables
    idx = idx + 1
    Restart(idx)%name = 'cv'
    Restart(idx)%axes = 2
    Restart(idx)%data%var2 => Cldprop%cv(:)
    
    idx = idx + 1
    Restart(idx)%name = 'cvt'
    Restart(idx)%axes = 2
    Restart(idx)%data%var2 => Cldprop%cvt(:)

    idx = idx + 1
    Restart(idx)%name = 'cvb'
    Restart(idx)%axes = 2
    Restart(idx)%data%var2 => Cldprop%cvb(:)

    !--- phy_f2d variables
    do num = 1,Model%ntot2d
      idx = idx + 1
       !--- set the variable name
      write(c2,'(i2.2)') num
      Restart(idx)%name = 'phy_f2d_'//c2
      Restart(idx)%axes = 2
      Restart(idx)%data%var2 => Tbd%phy_f2d(:,num)
    enddo

    !--- phy_fctd variables
    if (Model%nctp > 0) then
      do num = 1, Model%nctp
        idx = idx + 1
        !--- set the variable name
        write(c2,'(i2.2)') num
        Restart(idx)%name = 'phy_fctd_'//c2
        Restart(idx)%axes = 2
        Restart(idx)%data%var2 => Tbd%phy_fctd(:,num)
      enddo
    endif

    !--- Diagnostic variables
    do num = 1,ndiag_rst
      if( ndiag_idx(num) > 0 ) then
        idx = idx + 1
        Restart(idx)%name  = trim(ExtDiag(ndiag_idx(num))%name)
        Restart(idx)%axes  = 2
        Restart(idx)%diag  = .true.
        Restart(idx)%reset = .true.
        Restart(idx)%data%var2 => ExtDiag(ndiag_idx(num))%data%var2(:)
      endif
    enddo

    !--- Celluluar Automaton, 2D
    !CA
    if (Model%imfdeepcnv == 2 .and. Model%do_ca) then
      idx = idx + 1
      Restart(idx)%name = 'ca_condition'
      Restart(idx)%axes = 2
      Restart(idx)%data%var2 => Coupling%condition(:)
    endif
    ! Unified convection
    if (Model%imfdeepcnv == Model%imfdeepcnv_c3) then
      idx = idx + 1
      Restart(idx)%name = 'gf_2d_conv_act'
      Restart(idx)%axes = 2
      Restart(idx)%data%var2 => Sfcprop%conv_act(:)
      idx = idx + 1
      Restart(idx)%name = 'gf_2d_conv_act_m'
      Restart(idx)%axes = 2
      Restart(idx)%data%var2 => Sfcprop%conv_act_m(:)
      idx = idx + 1
      Restart(idx)%name = 'aod_gf'
      Restart(idx)%axes = 2
      Restart(idx)%data%var2 => Tbd%aod_gf(:)
    endif
    !--- RAP/HRRR-specific variables, 2D
    ! GF
    if (Model%imfdeepcnv == Model%imfdeepcnv_gf) then
      idx = idx + 1
      Restart(idx)%name = 'gf_2d_conv_act'
      Restart(idx)%axes = 2
      Restart(idx)%data%var2 => Sfcprop%conv_act(:)
      idx = idx + 1
      Restart(idx)%name = 'gf_2d_conv_act_m'
      Restart(idx)%axes = 2
      Restart(idx)%data%var2 => Sfcprop%conv_act_m(:)
      idx = idx + 1
      Restart(idx)%name = 'aod_gf'
      Restart(idx)%axes = 2
      Restart(idx)%data%var2 => Tbd%aod_gf(:)
    endif
    ! NoahMP
    if (Model%lsm == Model%lsm_noahmp) then
      idx = idx + 1
      Restart(idx)%name = 'noahmp_2d_raincprv'
      Restart(idx)%axes = 2
      Restart(idx)%data%var2 => Sfcprop%raincprv(:)
      idx = idx + 1
      Restart(idx)%name = 'noahmp_2d_rainncprv'
      Restart(idx)%axes = 2
      Restart(idx)%data%var2 => Sfcprop%rainncprv(:)
      idx = idx + 1
      Restart(idx)%name = 'noahmp_2d_iceprv'
      Restart(idx)%axes = 2
      Restart(idx)%data%var2 => Sfcprop%iceprv(:)
      idx = idx + 1
      Restart(idx)%name = 'noahmp_2d_snowprv'
      Restart(idx)%axes = 2
      Restart(idx)%data%var2 => Sfcprop%snowprv(:)
      idx = idx + 1
      Restart(idx)%name = 'noahmp_2d_graupelprv'
      Restart(idx)%axes = 2
      Restart(idx)%data%var2 => Sfcprop%graupelprv(:)
      idx = idx + 1
      Restart(idx)%name = 'noahmp_2d_draincprv'
      Restart(idx)%axes = 2
      Restart(idx)%data%var2 => Sfcprop%draincprv(:)
      idx = idx + 1
      Restart(idx)%name = 'noahmp_2d_drainncprv'
      Restart(idx)%axes = 2
      Restart(idx)%data%var2 => Sfcprop%drainncprv(:)
      idx = idx + 1
      Restart(idx)%name = 'noahmp_2d_diceprv'
      Restart(idx)%axes = 2
      Restart(idx)%data%var2 => Sfcprop%diceprv(:)
      idx = idx + 1
      Restart(idx)%name = 'noahmp_2d_dsnowprv'
      Restart(idx)%axes = 2
      Restart(idx)%data%var2 => Sfcprop%dsnowprv(:)
      idx = idx + 1
      Restart(idx)%name = 'noahmp_2d_dgraupelprv'
      Restart(idx)%axes = 2
      Restart(idx)%data%var2 => Sfcprop%dgraupelprv(:)
    endif
    ! RUC 
    if (Model%lsm == Model%lsm_ruc) then
      idx = idx + 1
      Restart(idx)%name = 'ruc_2d_raincprv'
      Restart(idx)%axes = 2
      Restart(idx)%data%var2 => Sfcprop%raincprv(:)
      idx = idx + 1
      Restart(idx)%name = 'ruc_2d_rainncprv'
      Restart(idx)%axes = 2
      Restart(idx)%data%var2 => Sfcprop%rainncprv(:)
      idx = idx + 1
      Restart(idx)%name = 'ruc_2d_iceprv'
      Restart(idx)%axes = 2
      Restart(idx)%data%var2 => Sfcprop%iceprv(:)
      idx = idx + 1
      Restart(idx)%name = 'ruc_2d_snowprv'
      Restart(idx)%axes = 2
      Restart(idx)%data%var2 => Sfcprop%snowprv(:)
      idx = idx + 1
      Restart(idx)%name = 'ruc_2d_graupelprv'
      Restart(idx)%axes = 2
      Restart(idx)%data%var2 => Sfcprop%graupelprv(:)
    endif
    ! MYNN SFC
    if (Model%do_mynnsfclay) then
        idx = idx + 1
        Restart(idx)%name = 'mynn_2d_uustar'
        Restart(idx)%axes = 2
        Restart(idx)%data%var2 => Sfcprop%uustar(:)
        idx = idx + 1
        Restart(idx)%name = 'mynn_2d_hpbl'
        Restart(idx)%axes = 2
        Restart(idx)%data%var2 => Tbd%hpbl(:)
        idx = idx + 1
        Restart(idx)%name = 'mynn_2d_ustm'
        Restart(idx)%axes = 2
        Restart(idx)%data%var2 => Sfcprop%ustm(:)
        idx = idx + 1
        Restart(idx)%name = 'mynn_2d_zol'
        Restart(idx)%axes = 2
        Restart(idx)%data%var2 => Sfcprop%zol(:)
        idx = idx + 1
        Restart(idx)%name = 'mynn_2d_mol'
        Restart(idx)%axes = 2
        Restart(idx)%data%var2 => Sfcprop%mol(:)
        idx = idx + 1
        Restart(idx)%name = 'mynn_2d_flhc'
        Restart(idx)%axes = 2
        Restart(idx)%data%var2 => Sfcprop%flhc(:)
        idx = idx + 1
        Restart(idx)%name = 'mynn_2d_flqc'
        Restart(idx)%axes = 2
        Restart(idx)%data%var2 => Sfcprop%flqc(:)
        idx = idx + 1
        Restart(idx)%name = 'mynn_2d_chs2'
        Restart(idx)%axes = 2
        Restart(idx)%data%var2 => Sfcprop%chs2(:)
        idx = idx + 1
        Restart(idx)%name = 'mynn_2d_cqs2'
        Restart(idx)%axes = 2
        Restart(idx)%data%var2 => Sfcprop%cqs2(:)
        idx = idx + 1
        Restart(idx)%name = 'mynn_2d_lh'
        Restart(idx)%axes = 2
        Restart(idx)%data%var2 => Sfcprop%lh(:)
        idx = idx + 1
        Restart(idx)%name = 'mynn_2d_hflx'
        Restart(idx)%axes = 2
        Restart(idx)%data%var2 => Sfcprop%hflx(:)
        idx = idx + 1
        Restart(idx)%name = 'mynn_2d_evap'
        Restart(idx)%axes = 2
        Restart(idx)%data%var2 => Sfcprop%evap(:)
        idx = idx + 1
        Restart(idx)%name = 'mynn_2d_qss'
        Restart(idx)%axes = 2
        Restart(idx)%data%var2 => Sfcprop%qss(:)
    endif
    ! Save rain prev for lake if surface layer doesn't.
    if (Model%lkm>0 .and. Model%iopt_lake==Model%iopt_lake_clm .and. &
         .not.surface_layer_saves_rainprev) then
      idx = idx + 1
      Restart(idx)%name = 'raincprv'
      Restart(idx)%axes = 2
      Restart(idx)%data%var2 => Sfcprop%raincprv(:)
      idx = idx + 1
      Restart(idx)%name = 'rainncprv'
      Restart(idx)%axes = 2
      Restart(idx)%data%var2 => Sfcprop%rainncprv(:)
    endif
    ! Thompson aerosol-aware
    if ((Model%imp_physics == Model%imp_physics_thompson .or. &
         Model%imp_physics == Model%imp_physics_tempo) .and. Model%ltaerosol) then
      idx = idx + 1
      Restart(idx)%name = 'thompson_2d_nwfa2d'
      Restart(idx)%axes = 2
      Restart(idx)%data%var2 => Coupling%nwfa2d(:)
      idx = idx + 1
      Restart(idx)%name = 'thompson_2d_nifa2d'
      Restart(idx)%axes = 2
      Restart(idx)%data%var2 => Coupling%nifa2d(:)
    endif

    ! Convection suppression
    if (Model%do_cap_suppress .and. Model%num_dfi_radar > 0) then
      do itime=1,Model%dfi_radar_max_intervals
        if(Model%ix_dfi_radar(itime)>0) then
          idx = idx + 1
          if(itime==1) then
            Restart(idx)%name = 'cap_suppress'
          else
            write(Restart(idx)%name,'("cap_suppress_",I0)') itime
          endif
          Restart(idx)%axes = 2
          Restart(idx)%data%var2 => Tbd%cap_suppress(:,Model%ix_dfi_radar(itime))
        endif
      enddo
    endif

    ! RRFS-SD
    if (Model%rrfs_sd) then
      idx = idx + 1
      Restart(idx)%name = 'ddvel_1'
      Restart(idx)%axes = 2
      Restart(idx)%data%var2 => Coupling%ddvel(:,1)
      idx = idx + 1
      Restart(idx)%name = 'ddvel_2'
      Restart(idx)%axes = 2
      Restart(idx)%data%var2 => Coupling%ddvel(:,2)
      idx = idx + 1
      Restart(idx)%name = 'min_fplume'
      Restart(idx)%axes = 2
      Restart(idx)%data%var2 => Coupling%min_fplume(:)
      idx = idx + 1
      Restart(idx)%name = 'max_fplume'
      Restart(idx)%axes = 2
      Restart(idx)%data%var2 => Coupling%max_fplume(:)
      idx = idx + 1
      Restart(idx)%name = 'rrfs_hwp'
      Restart(idx)%axes = 2
      Restart(idx)%data%var2 => Coupling%rrfs_hwp(:)
      idx = idx + 1
      Restart(idx)%name = 'rrfs_hwp_ave'
      Restart(idx)%axes = 2
      Restart(idx)%data%var2 => Coupling%rrfs_hwp_ave(:)
    endif

    !--- phy_f3d variables
    do num = 1,Model%ntot3d
      idx = idx + 1
      !--- set the variable name
      write(c2,'(i2.2)') num
      Restart(idx)%name = 'phy_f3d_'//c2
      Restart(idx)%axes = 3
      Restart(idx)%data%var3 => Tbd%phy_f3d(:,:,num)
   enddo

   if (Model%lrefres) then
      idx = idx + 1
      Restart(idx)%name = 'ref_f3d'
      Restart(idx)%axes = 3
      Restart(idx)%data%var3 => IntDiag%refl_10cm(:,:)
   endif

    !Prognostic closure
    if(Model%progsigma)then
       idx = idx + 1
       Restart(idx)%name = 'sas_3d_qgrs_dsave'
       Restart(idx)%axes = 3
       Restart(idx)%data%var3 => Tbd%prevsq(:,:)
       idx = idx + 1
       Restart(idx)%name = 'sas_3d_dqdt_qmicro'
       Restart(idx)%axes = 3
       Restart(idx)%data%var3 => Coupling%dqdt_qmicro(:,:)
    endif

    !--Convection variable used in CB cloud fraction. Presently this
    !--is only needed in sgscloud_radpre for imfdeepcnv == imfdeepcnv_gf.
    if (Model%imfdeepcnv == Model%imfdeepcnv_gf .or. Model%imfdeepcnv == Model%imfdeepcnv_c3) then
      idx = idx + 1
      Restart(idx)%name = 'cnv_3d_ud_mf'
      Restart(idx)%axes = 3
      Restart(idx)%data%var3 => Tbd%ud_mf(:,:)
    endif

    !Unified convection scheme                                                                                                                                                                    
    if (Model%imfdeepcnv == Model%imfdeepcnv_c3) then
      idx = idx + 1
      Restart(idx)%name = 'gf_3d_prevst'
      Restart(idx)%axes = 3
      Restart(idx)%data%var3 => Tbd%prevst(:,:)
      idx = idx + 1
      Restart(idx)%name = 'gf_3d_prevsq'
      Restart(idx)%axes = 3
      Restart(idx)%data%var3 => Tbd%prevsq(:,:)
      idx = idx + 1
      Restart(idx)%name = 'gf_3d_qci_conv'
      Restart(idx)%axes = 3
      Restart(idx)%data%var3 => Coupling%qci_conv(:,:)
    endif

    !--- RAP/HRRR-specific variables, 3D
    ! GF
    if (Model%imfdeepcnv == Model%imfdeepcnv_gf) then
      idx = idx + 1
      Restart(idx)%name = 'gf_3d_prevst'
      Restart(idx)%axes = 3
      Restart(idx)%data%var3 => Tbd%prevst(:,:)
      idx = idx + 1
      Restart(idx)%name = 'gf_3d_prevsq'
      Restart(idx)%axes = 3
      Restart(idx)%data%var3 => Tbd%prevsq(:,:)
      idx = idx + 1
      Restart(idx)%name = 'gf_3d_qci_conv'
      Restart(idx)%axes = 3
      Restart(idx)%data%var3 => Coupling%qci_conv(:,:)
    endif
    ! MYNN PBL
    if (Model%do_mynnedmf) then
      idx = idx + 1
      Restart(idx)%name = 'mynn_3d_cldfra_bl'
      Restart(idx)%axes = 3
      Restart(idx)%data%var3 => Tbd%cldfra_bl(:,:)
      idx = idx + 1
      Restart(idx)%name = 'mynn_3d_qc_bl'
      Restart(idx)%axes = 3
      Restart(idx)%data%var3 => Tbd%qc_bl(:,:)
      idx = idx + 1
      Restart(idx)%name = 'mynn_3d_qi_bl'
      Restart(idx)%axes = 3
      Restart(idx)%data%var3 => Tbd%qi_bl(:,:)
      idx = idx + 1
      Restart(idx)%name = 'mynn_3d_el_pbl'
      Restart(idx)%axes = 3
      Restart(idx)%data%var3 => Tbd%el_pbl(:,:)
      idx = idx + 1
      Restart(idx)%name = 'mynn_3d_sh3d'
      Restart(idx)%axes = 3
      Restart(idx)%data%var3 => Tbd%sh3d(:,:)
      idx = idx + 1
      Restart(idx)%name = 'mynn_3d_qke'
      Restart(idx)%axes = 3
      Restart(idx)%data%var3 => Tbd%qke(:,:)
      idx = idx + 1
      Restart(idx)%name = 'mynn_3d_tsq'
      Restart(idx)%axes = 3
      Restart(idx)%data%var3 => Tbd%tsq(:,:)
      idx = idx + 1
      Restart(idx)%name = 'mynn_3d_qsq'
      Restart(idx)%axes = 3
      Restart(idx)%data%var3 => Tbd%qsq(:,:)
      idx = idx + 1
      Restart(idx)%name = 'mynn_3d_cov'
      Restart(idx)%axes = 3
      Restart(idx)%data%var3 => Tbd%cov(:,:)
    endif

    ! Radar-derived microphysics temperature tendencies
    if (Model%num_dfi_radar > 0) then
      do itime=1,Model%dfi_radar_max_intervals
        if(Model%ix_dfi_radar(itime)>0) then
          idx = idx + 1
          if(itime==1) then
            Restart(idx)%name = 'radar_tten'
          else
            write(Restart(idx)%name,'("radar_tten_",I0)') itime
          endif
          Restart(idx)%axes = 3
          Restart(idx)%data%var3 => Tbd%dfi_radar_tten(:,:,Model%ix_dfi_radar(itime))
        endif
      enddo
    endif

    if(Model%rrfs_sd) then
      idx = idx + 1
      Restart(idx)%name = 'chem3d_1'
      Restart(idx)%axes = 3
      Restart(idx)%data%var3 => Coupling%chem3d(:,:,1)
      idx = idx + 1
      Restart(idx)%name = 'chem3d_2'
      Restart(idx)%axes = 3
      Restart(idx)%data%var3 => Coupling%chem3d(:,:,2)
      idx = idx + 1
      Restart(idx)%name = 'chem3d_3'
      Restart(idx)%axes = 3
      Restart(idx)%data%var3 => Coupling%chem3d(:,:,3)
      idx = idx + 1
      Restart(idx)%name = 'ext550'
      Restart(idx)%axes = 3
      Restart(idx)%data%var3 => Radtend%ext550(:,:)
    endif

  end subroutine GFS_restart_populate

end module GFS_restart
