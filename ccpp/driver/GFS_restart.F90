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
     real(kind=kind_phys), dimension(:),   pointer :: var2(:)   => null()
     real(kind=kind_phys), dimension(:,:), pointer :: var3(:,:) => null()
  end type var_subtype

  type GFS_restart_type
     integer           :: axes     !< Rank of data (2D or 3D).
     logical           :: diag     !< True for diagnostic field.
     logical           :: reset    !< If true, zero out diagnostic field.
     character(len=32) :: name     !< variable name as it will appear in the restart file.
     type(var_subtype) :: data(1)  !< Holds pointers to contiguous data.
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
!     Restart%data(1)%var2(:)   [real*8 ]  pointer to 2D data (im)                       !
!     Restart%data(1)%var3(:,:) [real*8 ]  pointer to 3D data (im,levs)                  !
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
    if (Model%imp_physics == Model%imp_physics_thompson .and. Model%ltaerosol) then
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
    print*,'SWALES GFS_restart_populate num2d+num3d = ',num2d+num3d
    !--- Allocate Restart data type.
    allocate (Restart(num2d+num3d))
    Restart(:)%diag  = .false.
    Restart(:)%reset = .false.

    !--- Cldprop variables
    Restart(1)%name = 'cv'
    Restart(1)%axes = 2
    Restart(1)%data(1)%var2 => Cldprop%cv(:)
    Restart(2)%name = 'cvt'
    Restart(2)%axes = 2
    Restart(2)%data(1)%var2 => Cldprop%cvt(:)
    Restart(3)%name = 'cvb'
    Restart(3)%axes = 2
    Restart(3)%data(1)%var2 => Cldprop%cvb(:)

    !--- phy_f2d variables
    offset = 3
    do num = 1,Model%ntot2d
       !--- set the variable name
      write(c2,'(i2.2)') num
      Restart(num+offset)%name = 'phy_f2d_'//c2
      Restart(num+offset)%axes = 2
      Restart(num+offset)%data(1)%var2 => Tbd%phy_f2d(:,num)
    enddo
    offset = offset + Model%ntot2d

    !--- phy_fctd variables
    if (Model%nctp > 0) then
      do num = 1, Model%nctp
       !--- set the variable name
        write(c2,'(i2.2)') num
        Restart(num+offset)%name = 'phy_fctd_'//c2
        Restart(num+offset)%axes = 2
        Restart(num+offset)%data(1)%var2 => Tbd%phy_fctd(:,num)
      enddo
      offset = offset + Model%nctp
    endif

    !--- Diagnostic variables
    do idx = 1,ndiag_rst
      if( ndiag_idx(idx) > 0 ) then
        Restart(offset+idx)%name  = trim(ExtDiag(ndiag_idx(idx))%name)
        Restart(offset+idx)%axes  = 2
        Restart(offset+idx)%diag  = .true.
        Restart(offset+idx)%reset = .true.
        Restart(offset+idx)%data(1)%var2 => ExtDiag(ndiag_idx(idx))%data(1)%var2(:)
      endif
    enddo

    num = offset + ndiag_rst
    !--- Celluluar Automaton, 2D
    !CA
    if (Model%imfdeepcnv == 2 .and. Model%do_ca) then
      num = num + 1
      Restart(num)%name = 'ca_condition'
      Restart(num)%axes = 2
      Restart(num)%data(1)%var2 => Coupling%condition(:)
    endif
    ! Unified convection
    if (Model%imfdeepcnv == Model%imfdeepcnv_c3) then
      num = num + 1
      Restart(num)%name = 'gf_2d_conv_act'
      Restart(num)%axes = 2
      Restart(num)%data(1)%var2 => Sfcprop%conv_act(:)
      num = num + 1
      Restart(num)%name = 'gf_2d_conv_act_m'
      Restart(num)%axes = 2
      Restart(num)%data(1)%var2 => Sfcprop%conv_act_m(:)
      num = num + 1
      Restart(num)%name = 'aod_gf'
      Restart(num)%axes = 2
      Restart(num)%data(1)%var2 => Tbd%aod_gf(:)
    endif
    !--- RAP/HRRR-specific variables, 2D
    ! GF
    if (Model%imfdeepcnv == Model%imfdeepcnv_gf) then
      num = num + 1
      Restart(num)%name = 'gf_2d_conv_act'
      Restart(num)%axes = 2
      Restart(num)%data(1)%var2 => Sfcprop%conv_act(:)
      num = num + 1
      Restart(num)%name = 'gf_2d_conv_act_m'
      Restart(num)%axes = 2
      Restart(num)%data(1)%var2 => Sfcprop%conv_act_m(:)
      num = num + 1
      Restart(num)%name = 'aod_gf'
      Restart(num)%axes = 2
      Restart(num)%data(1)%var2 => Tbd%aod_gf(:)
    endif
    ! NoahMP
    if (Model%lsm == Model%lsm_noahmp) then
      num = num + 1
      Restart(num)%name = 'noahmp_2d_raincprv'
      Restart(num)%axes = 2
      Restart(num)%data(1)%var2 => Sfcprop%raincprv(:)
      num = num + 1
      Restart(num)%name = 'noahmp_2d_rainncprv'
      Restart(num)%axes = 2
      Restart(num)%data(1)%var2 => Sfcprop%rainncprv(:)
      num = num + 1
      Restart(num)%name = 'noahmp_2d_iceprv'
      Restart(num)%axes = 2
      Restart(num)%data(1)%var2 => Sfcprop%iceprv(:)
      num = num + 1
      Restart(num)%name = 'noahmp_2d_snowprv'
      Restart(num)%axes = 2
      Restart(num)%data(1)%var2 => Sfcprop%snowprv(:)
      num = num + 1
      Restart(num)%name = 'noahmp_2d_graupelprv'
      Restart(num)%axes = 2
      Restart(num)%data(1)%var2 => Sfcprop%graupelprv(:)
      num = num + 1
      Restart(num)%name = 'noahmp_2d_draincprv'
      Restart(num)%axes = 2
      Restart(num)%data(1)%var2 => Sfcprop%draincprv(:)
      num = num + 1
      Restart(num)%name = 'noahmp_2d_drainncprv'
      Restart(num)%axes = 2
      Restart(num)%data(1)%var2 => Sfcprop%drainncprv(:)
      num = num + 1
      Restart(num)%name = 'noahmp_2d_diceprv'
      Restart(num)%axes = 2
      Restart(num)%data(1)%var2 => Sfcprop%diceprv(:)
      num = num + 1
      Restart(num)%name = 'noahmp_2d_dsnowprv'
      Restart(num)%axes = 2
      Restart(num)%data(1)%var2 => Sfcprop%dsnowprv(:)
      num = num + 1
      Restart(num)%name = 'noahmp_2d_dgraupelprv'
      Restart(num)%axes = 2
      Restart(num)%data(1)%var2 => Sfcprop%dgraupelprv(:)
    endif
    ! RUC 
    if (Model%lsm == Model%lsm_ruc) then
      num = num + 1
      Restart(num)%name = 'ruc_2d_raincprv'
      Restart(num)%axes = 2
      Restart(num)%data(1)%var2 => Sfcprop%raincprv(:)
      num = num + 1
      Restart(num)%name = 'ruc_2d_rainncprv'
      Restart(num)%axes = 2
      Restart(num)%data(1)%var2 => Sfcprop%rainncprv(:)
      num = num + 1
      Restart(num)%name = 'ruc_2d_iceprv'
      Restart(num)%axes = 2
      Restart(num)%data(1)%var2 => Sfcprop%iceprv(:)
      num = num + 1
      Restart(num)%name = 'ruc_2d_snowprv'
      Restart(num)%axes = 2
      Restart(num)%data(1)%var2 => Sfcprop%snowprv(:)
      num = num + 1
      Restart(num)%name = 'ruc_2d_graupelprv'
      Restart(num)%axes = 2
      Restart(num)%data(1)%var2 => Sfcprop%graupelprv(:)
    endif
    ! MYNN SFC
    if (Model%do_mynnsfclay) then
        num = num + 1
        Restart(num)%name = 'mynn_2d_uustar'
        Restart(num)%axes = 2
        Restart(num)%data(1)%var2 => Sfcprop%uustar(:)
        num = num + 1
        Restart(num)%name = 'mynn_2d_hpbl'
        Restart(num)%axes = 2
        Restart(num)%data(1)%var2 => Tbd%hpbl(:)
        num = num + 1
        Restart(num)%name = 'mynn_2d_ustm'
        Restart(num)%axes = 2
        Restart(num)%data(1)%var2 => Sfcprop%ustm(:)
        num = num + 1
        Restart(num)%name = 'mynn_2d_zol'
        Restart(num)%axes = 2
        Restart(num)%data(1)%var2 => Sfcprop%zol(:)
        num = num + 1
        Restart(num)%name = 'mynn_2d_mol'
        Restart(num)%axes = 2
        Restart(num)%data(1)%var2 => Sfcprop%mol(:)
        num = num + 1
        Restart(num)%name = 'mynn_2d_flhc'
        Restart(num)%axes = 2
        Restart(num)%data(1)%var2 => Sfcprop%flhc(:)
        num = num + 1
        Restart(num)%name = 'mynn_2d_flqc'
        Restart(num)%axes = 2
        Restart(num)%data(1)%var2 => Sfcprop%flqc(:)
        num = num + 1
        Restart(num)%name = 'mynn_2d_chs2'
        Restart(num)%axes = 2
        Restart(num)%data(1)%var2 => Sfcprop%chs2(:)
        num = num + 1
        Restart(num)%name = 'mynn_2d_cqs2'
        Restart(num)%axes = 2
        Restart(num)%data(1)%var2 => Sfcprop%cqs2(:)
        num = num + 1
        Restart(num)%name = 'mynn_2d_lh'
        Restart(num)%axes = 2
        Restart(num)%data(1)%var2 => Sfcprop%lh(:)
        num = num + 1
        Restart(num)%name = 'mynn_2d_hflx'
        Restart(num)%axes = 2
        Restart(num)%data(1)%var2 => Sfcprop%hflx(:)
        num = num + 1
        Restart(num)%name = 'mynn_2d_evap'
        Restart(num)%axes = 2
        Restart(num)%data(1)%var2 => Sfcprop%evap(:)
        num = num + 1
        Restart(num)%name = 'mynn_2d_qss'
        Restart(num)%axes = 2
        Restart(num)%data(1)%var2 => Sfcprop%qss(:)
    endif
    ! Save rain prev for lake if surface layer doesn't.
    if (Model%lkm>0 .and. Model%iopt_lake==Model%iopt_lake_clm .and. &
         .not.surface_layer_saves_rainprev) then
      num = num + 1
      Restart(num)%name = 'raincprv'
      Restart(num)%axes = 2
      Restart(num)%data(1)%var2 => Sfcprop%raincprv(:)
      num = num + 1
      Restart(num)%name = 'rainncprv'
      Restart(num)%axes = 2
      Restart(num)%data(1)%var2 => Sfcprop%rainncprv(:)
    endif
    ! Thompson aerosol-aware
    if (Model%imp_physics == Model%imp_physics_thompson .and. Model%ltaerosol) then
      num = num + 1
      Restart(num)%name = 'thompson_2d_nwfa2d'
      Restart(num)%axes = 2
      Restart(num)%data(1)%var2 => Coupling%nwfa2d(:)
      num = num + 1
      Restart(num)%name = 'thompson_2d_nifa2d'
      Restart(num)%axes = 2
      Restart(num)%data(1)%var2 => Coupling%nifa2d(:)
    endif

    ! Convection suppression
    if (Model%do_cap_suppress .and. Model%num_dfi_radar > 0) then
      do itime=1,Model%dfi_radar_max_intervals
        if(Model%ix_dfi_radar(itime)>0) then
          num = num + 1
          if(itime==1) then
            Restart(num)%name = 'cap_suppress'
          else
            write(Restart(num)%name,'("cap_suppress_",I0)') itime
          endif
          Restart(num)%axes = 2
          Restart(num)%data(1)%var2 => Tbd%cap_suppress(:,Model%ix_dfi_radar(itime))
        endif
      enddo
    endif

    ! RRFS-SD
    if (Model%rrfs_sd) then
      num = num + 1
      Restart(num)%name = 'ddvel_1'
      Restart(num)%axes = 2
      Restart(num)%data(1)%var2 => Coupling%ddvel(:,1)
      num = num + 1
      Restart(num)%name = 'ddvel_2'
      Restart(num)%axes = 2
      Restart(num)%data(1)%var2 => Coupling%ddvel(:,2)
      num = num + 1
      Restart(num)%name = 'min_fplume'
      Restart(num)%axes = 2
      Restart(num)%data(1)%var2 => Coupling%min_fplume(:)
      num = num + 1
      Restart(num)%name = 'max_fplume'
      Restart(num)%axes = 2
      Restart(num)%data(1)%var2 => Coupling%max_fplume(:)
      num = num + 1
      Restart(num)%name = 'rrfs_hwp'
      Restart(num)%axes = 2
      Restart(num)%data(1)%var2 => Coupling%rrfs_hwp(:)
      num = num + 1
      Restart(num)%name = 'rrfs_hwp_ave'
      Restart(num)%axes = 2
      Restart(num)%data(1)%var2 => Coupling%rrfs_hwp_ave(:)
    endif

    !--- phy_f3d variables
    do num = 1,Model%ntot3d
       !--- set the variable name
      write(c2,'(i2.2)') num
      Restart(num)%name = 'phy_f3d_'//c2
      Restart(num)%axes = 3
      Restart(num)%data(1)%var3 => Tbd%phy_f3d(:,:,num)
    enddo
    if (Model%lrefres) then
      num = Model%ntot3d+1
      Restart(num)%name = 'ref_f3d'
      Restart(num)%axes = 3
      Restart(num)%data(1)%var3 => IntDiag%refl_10cm(:,:)
    endif
    if (Model%lrefres) then
       num = Model%ntot3d+1
    else
       num = Model%ntot3d
    endif

    !Prognostic closure
    if(Model%progsigma)then
       num = num + 1
       Restart(num)%name = 'sas_3d_qgrs_dsave'
       Restart(num)%axes = 3
       Restart(num)%data(1)%var3 => Tbd%prevsq(:,:)
       num = num + 1
       Restart(num)%name = 'sas_3d_dqdt_qmicro'
       Restart(num)%axes = 3
       Restart(num)%data(1)%var3 => Coupling%dqdt_qmicro(:,:)
    endif

    !--Convection variable used in CB cloud fraction. Presently this
    !--is only needed in sgscloud_radpre for imfdeepcnv == imfdeepcnv_gf.
    if (Model%imfdeepcnv == Model%imfdeepcnv_gf .or. Model%imfdeepcnv == Model%imfdeepcnv_c3) then
      num = num + 1
      Restart(num)%name = 'cnv_3d_ud_mf'
      Restart(num)%axes = 3
      Restart(num)%data(1)%var3 => Tbd%ud_mf(:,:)
    endif

    !Unified convection scheme                                                                                                                                                                    
    if (Model%imfdeepcnv == Model%imfdeepcnv_c3) then
      num = num + 1
      Restart(num)%name = 'gf_3d_prevst'
      Restart(num)%axes = 3
      Restart(num)%data(1)%var3 => Tbd%prevst(:,:)
      num = num + 1
      Restart(num)%name = 'gf_3d_prevsq'
      Restart(num)%axes = 3
      Restart(num)%data(1)%var3 => Tbd%prevsq(:,:)
      num = num + 1
      Restart(num)%name = 'gf_3d_qci_conv'
      Restart(num)%axes = 3
      Restart(num)%data(1)%var3 => Coupling%qci_conv(:,:)
    endif

    !--- RAP/HRRR-specific variables, 3D
    ! GF
    if (Model%imfdeepcnv == Model%imfdeepcnv_gf) then
      num = num + 1
      Restart(num)%name = 'gf_3d_prevst'
      Restart(num)%axes = 3
      Restart(num)%data(1)%var3 => Tbd%prevst(:,:)
      num = num + 1
      Restart(num)%name = 'gf_3d_prevsq'
      Restart(num)%axes = 3
      Restart(num)%data(1)%var3 => Tbd%prevsq(:,:)
      num = num + 1
      Restart(num)%name = 'gf_3d_qci_conv'
      Restart(num)%axes = 3
      Restart(num)%data(1)%var3 => Coupling%qci_conv(:,:)
    endif
    ! MYNN PBL
    if (Model%do_mynnedmf) then
      num = num + 1
      Restart(num)%name = 'mynn_3d_cldfra_bl'
      Restart(num)%axes = 3
      Restart(num)%data(1)%var3 => Tbd%cldfra_bl(:,:)
      num = num + 1
      Restart(num)%name = 'mynn_3d_qc_bl'
      Restart(num)%axes = 3
      Restart(num)%data(1)%var3 => Tbd%qc_bl(:,:)
      num = num + 1
      Restart(num)%name = 'mynn_3d_qi_bl'
      Restart(num)%axes = 3
      Restart(num)%data(1)%var3 => Tbd%qi_bl(:,:)
      num = num + 1
      Restart(num)%name = 'mynn_3d_el_pbl'
      Restart(num)%axes = 3
      Restart(num)%data(1)%var3 => Tbd%el_pbl(:,:)
      num = num + 1
      Restart(num)%name = 'mynn_3d_sh3d'
      Restart(num)%axes = 3
      Restart(num)%data(1)%var3 => Tbd%sh3d(:,:)
      num = num + 1
      Restart(num)%name = 'mynn_3d_qke'
      Restart(num)%axes = 3
      Restart(num)%data(1)%var3 => Tbd%qke(:,:)
      num = num + 1
      Restart(num)%name = 'mynn_3d_tsq'
      Restart(num)%axes = 3
      Restart(num)%data(1)%var3 => Tbd%tsq(:,:)
      num = num + 1
      Restart(num)%name = 'mynn_3d_qsq'
      Restart(num)%axes = 3
      Restart(num)%data(1)%var3 => Tbd%qsq(:,:)
      num = num + 1
      Restart(num)%name = 'mynn_3d_cov'
      Restart(num)%axes = 3
      Restart(num)%data(1)%var3 => Tbd%cov(:,:)
    endif

    ! Radar-derived microphysics temperature tendencies
    if (Model%num_dfi_radar > 0) then
      do itime=1,Model%dfi_radar_max_intervals
        if(Model%ix_dfi_radar(itime)>0) then
          num = num + 1
          if(itime==1) then
            Restart(num)%name = 'radar_tten'
          else
            write(Restart(num)%name,'("radar_tten_",I0)') itime
          endif
          Restart(num)%axes = 3
          Restart(num)%data(1)%var3 => Tbd%dfi_radar_tten(:,:,Model%ix_dfi_radar(itime))
        endif
      enddo
    endif

    if(Model%rrfs_sd) then
      num = num + 1
      Restart(num)%name = 'chem3d_1'
      Restart(num)%axes = 3
      Restart(num)%data(1)%var3 => Coupling%chem3d(:,:,1)
      num = num + 1
      Restart(num)%name = 'chem3d_2'
      Restart(num)%axes = 3
      Restart(num)%data(1)%var3 => Coupling%chem3d(:,:,2)
      num = num + 1
      Restart(num)%name = 'chem3d_3'
      Restart(num)%axes = 3
      Restart(num)%data(1)%var3 => Coupling%chem3d(:,:,3)
      num = num + 1
      Restart(num)%name = 'ext550'
      Restart(num)%axes = 3
      Restart(num)%data(1)%var3 => Radtend%ext550(:,:)
    endif

  end subroutine GFS_restart_populate

end module GFS_restart
