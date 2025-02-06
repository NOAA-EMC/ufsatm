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
    real(kind=kind_phys), pointer :: var2p(:)   => null()  !< 2D data saved in packed format [dim(ix)]
    real(kind=kind_phys), pointer :: var3p(:,:) => null()  !< 3D data saved in packed format [dim(ix,levs)]
  end type var_subtype

  type GFS_restart_type
    integer           :: num2d                    !< current number of registered 2D restart variables
    integer           :: num3d                    !< current number of registered 3D restart variables
    integer           :: fdiag                    !< index of first diagnostic field in restart file
    integer           :: ldiag                    !< index of last diagnostic field in restart file

    character(len=32), allocatable :: name2d(:)   !< variable name as it will appear in the restart file
    character(len=32), allocatable :: name3d(:)   !< variable name as it will appear in the restart file
    type(var_subtype), allocatable :: data(:,:)   !< holds pointers to data in packed format (allocated to (nblks,max(2d/3dfields))
  end type GFS_restart_type

  public GFS_restart_type, GFS_restart_populate

  CONTAINS
!*******************************************************************************************

!---------------------
! GFS_restart_populate
!---------------------
  subroutine GFS_restart_populate (Restart, Model, Statein, Stateout, Sfcprop, &
                                   Coupling, Grid, Tbd, Cldprop, Radtend, IntDiag, Init_parm, ExtDiag)
!----------------------------------------------------------------------------------------!
!   RESTART_METADATA                                                                         !
!     Restart%num2d          [int*4  ]  number of 2D variables to output             !
!     Restart%num3d          [int*4  ]  number of 3D variables to output             !
!     Restart%name2d         [char=32]  variable name in restart file                !
!     Restart%name3d         [char=32]  variable name in restart file                !
!     Restart%fld2d(:,:,:)   [real*8 ]  pointer to 2D data (im,nblks,MAX_RSTRT)      !
!     Restart%fld3d(:,:,:,:) [real*8 ]  pointer to 3D data (im,levs,nblks,MAX_RSTRT) !
!----------------------------------------------------------------------------------------!
    type(GFS_restart_type),     intent(inout) :: Restart
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
    integer :: nblks, num, nb, max_rstrt, offset 
    character(len=2) :: c2 = ''
    logical :: surface_layer_saves_rainprev
    
    nblks = size(Init_parm%blksz)
    max_rstrt = size(Restart%name2d)

    !--- check if continuous accumulated total precip and total cnvc precip are
    !requested in output 
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

    ! Store first and last index of diagnostic fields:
    Restart%fdiag = 3 + Model%ntot2d + Model%nctp + 1
    Restart%ldiag = 3 + Model%ntot2d + Model%nctp + ndiag_rst
    Restart%num2d = 3 + Model%ntot2d + Model%nctp + ndiag_rst

    ! The CLM Lake Model needs raincprev and rainncprv, which some
    ! surface layer schemes save, and some don't. If the surface layer
    ! scheme does not save that variable, then it'll be saved
    ! separately for clm_lake.
    surface_layer_saves_rainprev = .false.

    ! GF
    if (Model%imfdeepcnv == Model%imfdeepcnv_gf) then
      Restart%num2d = Restart%num2d + 3
    endif
    ! Unified convection
    if (Model%imfdeepcnv == Model%imfdeepcnv_c3) then
      Restart%num2d = Restart%num2d + 3
    endif
    ! CA
    if (Model%imfdeepcnv == 2 .and. Model%do_ca) then
      Restart%num2d = Restart%num2d + 1
    endif
    ! NoahMP
    if (Model%lsm == Model%lsm_noahmp) then
      Restart%num2d = Restart%num2d + 10
      surface_layer_saves_rainprev = .true.
    endif
    ! RUC 
    if (Model%lsm == Model%lsm_ruc) then
      Restart%num2d = Restart%num2d + 5
      surface_layer_saves_rainprev = .true.
    endif
    ! MYNN SFC
    if (Model%do_mynnsfclay) then
      Restart%num2d = Restart%num2d + 13
    endif
    ! Save rain prev for lake if surface layer doesn't.
    if (Model%lkm>0 .and. Model%iopt_lake==Model%iopt_lake_clm .and. &
         .not.surface_layer_saves_rainprev) then
      Restart%num2d = Restart%num2d + 2
    endif
    ! Thompson aerosol-aware
    if (Model%imp_physics == Model%imp_physics_thompson .and. Model%ltaerosol) then
      Restart%num2d = Restart%num2d + 2
    endif
    if (Model%do_cap_suppress .and. Model%num_dfi_radar>0) then
      Restart%num2d = Restart%num2d + Model%num_dfi_radar
    endif
    if (Model%rrfs_sd) then
      Restart%num2d = Restart%num2d + 6
    endif

    Restart%num3d = Model%ntot3d
    if (Model%num_dfi_radar>0) then
      Restart%num3d = Restart%num3d + Model%num_dfi_radar
    endif
    if(Model%lrefres) then
       Restart%num3d = Model%ntot3d+1
    endif
    ! General Convection
    if (Model%imfdeepcnv == Model%imfdeepcnv_gf) then
      Restart%num3d = Restart%num3d + 1
    endif
    ! GF
    if (Model%imfdeepcnv == 3) then
      Restart%num3d = Restart%num3d + 3
    endif
    ! Unified convection
    if (Model%imfdeepcnv == 5) then
      Restart%num3d = Restart%num3d + 4
    endif
    ! MYNN PBL
    if (Model%do_mynnedmf) then
      Restart%num3d = Restart%num3d + 9
    endif
    if (Model%rrfs_sd) then
      Restart%num3d = Restart%num3d + 4
    endif
    !Prognostic area fraction
    if (Model%progsigma) then
       Restart%num3d = Restart%num3d + 2
    endif

    if (Model%num_dfi_radar > 0) then
      do itime=1,Model%dfi_radar_max_intervals
        if(Model%ix_dfi_radar(itime)>0) then
          Restart%num3d = Restart%num3d + 1
        endif
      enddo
    endif

    allocate (Restart%name2d(Restart%num2d))
    allocate (Restart%name3d(Restart%num3d))
    allocate (Restart%data(nblks,max(Restart%num2d,Restart%num3d)))

    Restart%name2d(:) = ' '
    Restart%name3d(:) = ' '

    !--- Cldprop variables
    Restart%name2d(1) = 'cv'
    Restart%name2d(2) = 'cvt'
    Restart%name2d(3) = 'cvb'
    do nb = 1,nblks
      Restart%data(nb,1)%var2p => Cldprop%cv(Model%chunk_begin(nb):Model%chunk_end(nb))
      Restart%data(nb,2)%var2p => Cldprop%cvt(Model%chunk_begin(nb):Model%chunk_end(nb))
      Restart%data(nb,3)%var2p => Cldprop%cvb(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    !--- phy_f2d variables
    offset = 3
    do num = 1,Model%ntot2d
       !--- set the variable name
      write(c2,'(i2.2)') num
      Restart%name2d(num+offset) = 'phy_f2d_'//c2
      do nb = 1,nblks
        Restart%data(nb,num+offset)%var2p => Tbd%phy_f2d(Model%chunk_begin(nb):Model%chunk_end(nb),num)
      enddo
    enddo
    offset = offset + Model%ntot2d

    !--- phy_fctd variables
    if (Model%nctp > 0) then
      do num = 1, Model%nctp
       !--- set the variable name
        write(c2,'(i2.2)') num
        Restart%name2d(num+offset) = 'phy_fctd_'//c2
        do nb = 1,nblks
          Restart%data(nb,num+offset)%var2p => Tbd%phy_fctd(Model%chunk_begin(nb):Model%chunk_end(nb),num)
        enddo
      enddo
      offset = offset + Model%nctp
    endif

    !--- Diagnostic variables
    do idx = 1,ndiag_rst
      if( ndiag_idx(idx) > 0 ) then
        Restart%name2d(offset+idx) = trim(ExtDiag(ndiag_idx(idx))%name)
        do nb = 1,nblks
          Restart%data(nb,offset+idx)%var2p => ExtDiag(ndiag_idx(idx))%data(nb)%var2
        enddo
      endif
!      print *,'in restart 2d field, Restart%name2d(',offset+idx,')=',trim(Restart%name2d(offset+idx))
    enddo

    num = offset + ndiag_rst
    !--- Celluluar Automaton, 2D
    !CA
    if (Model%imfdeepcnv == 2 .and. Model%do_ca) then
      num = num + 1
      Restart%name2d(num) = 'ca_condition'
      do nb = 1,nblks
        Restart%data(nb,num)%var2p => Coupling%condition(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo
    endif
    ! Unified convection
    if (Model%imfdeepcnv == Model%imfdeepcnv_c3) then
      num = num + 1
      Restart%name2d(num) = 'gf_2d_conv_act'
      do nb = 1,nblks
        Restart%data(nb,num)%var2p => Sfcprop%conv_act(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo
      num = num + 1
      Restart%name2d(num) = 'gf_2d_conv_act_m'
      do nb = 1,nblks
        Restart%data(nb,num)%var2p => Sfcprop%conv_act_m(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo
      num = num + 1
      Restart%name2d(num) = 'aod_gf'
      do nb = 1,nblks
        Restart%data(nb,num)%var2p => Tbd%aod_gf(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo
    endif
    !--- RAP/HRRR-specific variables, 2D
    ! GF
    if (Model%imfdeepcnv == Model%imfdeepcnv_gf) then
      num = num + 1
      Restart%name2d(num) = 'gf_2d_conv_act'
      do nb = 1,nblks
        Restart%data(nb,num)%var2p => Sfcprop%conv_act(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo
      num = num + 1
      Restart%name2d(num) = 'gf_2d_conv_act_m'
      do nb = 1,nblks
        Restart%data(nb,num)%var2p => Sfcprop%conv_act_m(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo
      num = num + 1
      Restart%name2d(num) = 'aod_gf'
      do nb = 1,nblks
        Restart%data(nb,num)%var2p => Tbd%aod_gf(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo
    endif
    ! NoahMP
    if (Model%lsm == Model%lsm_noahmp) then
      num = num + 1
      Restart%name2d(num) = 'noahmp_2d_raincprv'
      do nb = 1,nblks
        Restart%data(nb,num)%var2p => Sfcprop%raincprv(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo
      num = num + 1
      Restart%name2d(num) = 'noahmp_2d_rainncprv'
      do nb = 1,nblks
        Restart%data(nb,num)%var2p => Sfcprop%rainncprv(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo
      num = num + 1
      Restart%name2d(num) = 'noahmp_2d_iceprv'
      do nb = 1,nblks
        Restart%data(nb,num)%var2p => Sfcprop%iceprv(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo
      num = num + 1
      Restart%name2d(num) = 'noahmp_2d_snowprv'
      do nb = 1,nblks
        Restart%data(nb,num)%var2p => Sfcprop%snowprv(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo
      num = num + 1
      Restart%name2d(num) = 'noahmp_2d_graupelprv'
      do nb = 1,nblks
        Restart%data(nb,num)%var2p => Sfcprop%graupelprv(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo
      num = num + 1
      Restart%name2d(num) = 'noahmp_2d_draincprv'
      do nb = 1,nblks
        Restart%data(nb,num)%var2p => Sfcprop%draincprv(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo
      num = num + 1
      Restart%name2d(num) = 'noahmp_2d_drainncprv'
      do nb = 1,nblks
        Restart%data(nb,num)%var2p => Sfcprop%drainncprv(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo
      num = num + 1
      Restart%name2d(num) = 'noahmp_2d_diceprv'
      do nb = 1,nblks
        Restart%data(nb,num)%var2p => Sfcprop%diceprv(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo
      num = num + 1
      Restart%name2d(num) = 'noahmp_2d_dsnowprv'
      do nb = 1,nblks
        Restart%data(nb,num)%var2p => Sfcprop%dsnowprv(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo
      num = num + 1
      Restart%name2d(num) = 'noahmp_2d_dgraupelprv'
      do nb = 1,nblks
        Restart%data(nb,num)%var2p => Sfcprop%dgraupelprv(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo
    endif
    ! RUC 
    if (Model%lsm == Model%lsm_ruc) then
      num = num + 1
      Restart%name2d(num) = 'ruc_2d_raincprv'
      do nb = 1,nblks
        Restart%data(nb,num)%var2p => Sfcprop%raincprv(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo
      num = num + 1
      Restart%name2d(num) = 'ruc_2d_rainncprv'
      do nb = 1,nblks
        Restart%data(nb,num)%var2p => Sfcprop%rainncprv(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo
      num = num + 1
      Restart%name2d(num) = 'ruc_2d_iceprv'
      do nb = 1,nblks
        Restart%data(nb,num)%var2p => Sfcprop%iceprv(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo
      num = num + 1
      Restart%name2d(num) = 'ruc_2d_snowprv'
      do nb = 1,nblks
        Restart%data(nb,num)%var2p => Sfcprop%snowprv(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo
      num = num + 1
      Restart%name2d(num) = 'ruc_2d_graupelprv'
      do nb = 1,nblks
        Restart%data(nb,num)%var2p => Sfcprop%graupelprv(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo
    endif
    ! MYNN SFC
    if (Model%do_mynnsfclay) then
        num = num + 1
        Restart%name2d(num) = 'mynn_2d_uustar'
        do nb = 1,nblks
          Restart%data(nb,num)%var2p => Sfcprop%uustar(Model%chunk_begin(nb):Model%chunk_end(nb))
        enddo
        num = num + 1
        Restart%name2d(num) = 'mynn_2d_hpbl'
        do nb = 1,nblks
          Restart%data(nb,num)%var2p => Tbd%hpbl(Model%chunk_begin(nb):Model%chunk_end(nb))
        enddo
        num = num + 1
        Restart%name2d(num) = 'mynn_2d_ustm'
        do nb = 1,nblks
          Restart%data(nb,num)%var2p => Sfcprop%ustm(Model%chunk_begin(nb):Model%chunk_end(nb))
        enddo
        num = num + 1
        Restart%name2d(num) = 'mynn_2d_zol'
        do nb = 1,nblks
          Restart%data(nb,num)%var2p => Sfcprop%zol(Model%chunk_begin(nb):Model%chunk_end(nb))
        enddo
        num = num + 1
        Restart%name2d(num) = 'mynn_2d_mol'
        do nb = 1,nblks
          Restart%data(nb,num)%var2p => Sfcprop%mol(Model%chunk_begin(nb):Model%chunk_end(nb))
        enddo
        num = num + 1
        Restart%name2d(num) = 'mynn_2d_flhc'
        do nb = 1,nblks
          Restart%data(nb,num)%var2p => Sfcprop%flhc(Model%chunk_begin(nb):Model%chunk_end(nb))
        enddo
        num = num + 1
        Restart%name2d(num) = 'mynn_2d_flqc'
        do nb = 1,nblks
          Restart%data(nb,num)%var2p => Sfcprop%flqc(Model%chunk_begin(nb):Model%chunk_end(nb))
        enddo
        num = num + 1
        Restart%name2d(num) = 'mynn_2d_chs2'
        do nb = 1,nblks
          Restart%data(nb,num)%var2p => Sfcprop%chs2(Model%chunk_begin(nb):Model%chunk_end(nb))
        enddo
        num = num + 1
        Restart%name2d(num) = 'mynn_2d_cqs2'
        do nb = 1,nblks
          Restart%data(nb,num)%var2p => Sfcprop%cqs2(Model%chunk_begin(nb):Model%chunk_end(nb))
        enddo
        num = num + 1
        Restart%name2d(num) = 'mynn_2d_lh'
        do nb = 1,nblks
          Restart%data(nb,num)%var2p => Sfcprop%lh(Model%chunk_begin(nb):Model%chunk_end(nb))
        enddo
        num = num + 1
        Restart%name2d(num) = 'mynn_2d_hflx'
        do nb = 1,nblks
          Restart%data(nb,num)%var2p => Sfcprop%hflx(Model%chunk_begin(nb):Model%chunk_end(nb))
        enddo
        num = num + 1
        Restart%name2d(num) = 'mynn_2d_evap'
        do nb = 1,nblks
          Restart%data(nb,num)%var2p => Sfcprop%evap(Model%chunk_begin(nb):Model%chunk_end(nb))
        enddo
        num = num + 1
        Restart%name2d(num) = 'mynn_2d_qss'
        do nb = 1,nblks
          Restart%data(nb,num)%var2p => Sfcprop%qss(Model%chunk_begin(nb):Model%chunk_end(nb))
        enddo
    endif
    ! Save rain prev for lake if surface layer doesn't.
    if (Model%lkm>0 .and. Model%iopt_lake==Model%iopt_lake_clm .and. &
         .not.surface_layer_saves_rainprev) then
      num = num + 1
      Restart%name2d(num) = 'raincprv'
      do nb = 1,nblks
        Restart%data(nb,num)%var2p => Sfcprop%raincprv(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo
      num = num + 1
      Restart%name2d(num) = 'rainncprv'
      do nb = 1,nblks
        Restart%data(nb,num)%var2p => Sfcprop%rainncprv(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo
    endif
    ! Thompson aerosol-aware
    if (Model%imp_physics == Model%imp_physics_thompson .and. Model%ltaerosol) then
      num = num + 1
      Restart%name2d(num) = 'thompson_2d_nwfa2d'
      do nb = 1,nblks
        Restart%data(nb,num)%var2p => Coupling%nwfa2d(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo
      num = num + 1
      Restart%name2d(num) = 'thompson_2d_nifa2d'
      do nb = 1,nblks
        Restart%data(nb,num)%var2p => Coupling%nifa2d(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo
    endif

    ! Convection suppression
    if (Model%do_cap_suppress .and. Model%num_dfi_radar > 0) then
      do itime=1,Model%dfi_radar_max_intervals
        if(Model%ix_dfi_radar(itime)>0) then
          num = num + 1
          if(itime==1) then
            Restart%name2d(num) = 'cap_suppress'
          else
            write(Restart%name2d(num),'("cap_suppress_",I0)') itime
          endif
          do nb = 1,nblks
            Restart%data(nb,num)%var2p => Tbd%cap_suppress(Model%chunk_begin(nb):Model%chunk_end(nb),Model%ix_dfi_radar(itime))
          enddo
        endif
      enddo
    endif

    ! RRFS-SD
    if (Model%rrfs_sd) then
      num = num + 1
      Restart%name2d(num) = 'ddvel_1'
      do nb = 1,nblks
        Restart%data(nb,num)%var2p => Coupling%ddvel(Model%chunk_begin(nb):Model%chunk_end(nb),1)
      enddo
      num = num + 1
      Restart%name2d(num) = 'ddvel_2'
      do nb = 1,nblks
        Restart%data(nb,num)%var2p => Coupling%ddvel(Model%chunk_begin(nb):Model%chunk_end(nb),2)
      enddo
      num = num + 1
      Restart%name2d(num) = 'min_fplume'
      do nb = 1,nblks
        Restart%data(nb,num)%var2p => Coupling%min_fplume(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo
      num = num + 1
      Restart%name2d(num) = 'max_fplume'
      do nb = 1,nblks
        Restart%data(nb,num)%var2p => Coupling%max_fplume(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo
      num = num + 1
      Restart%name2d(num) = 'rrfs_hwp'
      do nb = 1,nblks
        Restart%data(nb,num)%var2p => Coupling%rrfs_hwp(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo
      num = num + 1
      Restart%name2d(num) = 'rrfs_hwp_ave'
      do nb = 1,nblks
        Restart%data(nb,num)%var2p => Coupling%rrfs_hwp_ave(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo
    endif

    !--- phy_f3d variables
    do num = 1,Model%ntot3d
       !--- set the variable name
      write(c2,'(i2.2)') num
      Restart%name3d(num) = 'phy_f3d_'//c2
      do nb = 1,nblks
        Restart%data(nb,num)%var3p => Tbd%phy_f3d(Model%chunk_begin(nb):Model%chunk_end(nb),:,num)
      enddo
    enddo
    if (Model%lrefres) then
      num = Model%ntot3d+1
      restart%name3d(num) = 'ref_f3d'
      do nb = 1,nblks
        Restart%data(nb,num)%var3p => IntDiag%refl_10cm(Model%chunk_begin(nb):Model%chunk_end(nb),:)
      enddo
    endif
    if (Model%lrefres) then
       num = Model%ntot3d+1
    else
       num = Model%ntot3d
    endif

    !Prognostic closure
    if(Model%progsigma)then
       num = num + 1
       Restart%name3d(num) = 'sas_3d_qgrs_dsave'
      do nb = 1,nblks
        Restart%data(nb,num)%var3p => Tbd%prevsq(Model%chunk_begin(nb):Model%chunk_end(nb),:)
      enddo
      num = num + 1
      Restart%name3d(num) = 'sas_3d_dqdt_qmicro'
      do nb = 1,nblks
        Restart%data(nb,num)%var3p => Coupling%dqdt_qmicro(Model%chunk_begin(nb):Model%chunk_end(nb),:)
      enddo
    endif

    !--Convection variable used in CB cloud fraction. Presently this
    !--is only needed in sgscloud_radpre for imfdeepcnv == imfdeepcnv_gf.
    if (Model%imfdeepcnv == Model%imfdeepcnv_gf .or. Model%imfdeepcnv == Model%imfdeepcnv_c3) then
      num = num + 1
      Restart%name3d(num) = 'cnv_3d_ud_mf'
      do nb = 1,nblks
        Restart%data(nb,num)%var3p => Tbd%ud_mf(Model%chunk_begin(nb):Model%chunk_end(nb),:)
      enddo
    endif

    !Unified convection scheme                                                                                                                                                                    
    if (Model%imfdeepcnv == Model%imfdeepcnv_c3) then
      num = num + 1
      Restart%name3d(num) = 'gf_3d_prevst'
      do nb = 1,nblks
        Restart%data(nb,num)%var3p => Tbd%prevst(Model%chunk_begin(nb):Model%chunk_end(nb),:)
      enddo
      num = num + 1
      Restart%name3d(num) = 'gf_3d_prevsq'
      do nb = 1,nblks
        Restart%data(nb,num)%var3p => Tbd%prevsq(Model%chunk_begin(nb):Model%chunk_end(nb),:)
      enddo
      num = num + 1
      Restart%name3d(num) = 'gf_3d_qci_conv'
      do nb = 1,nblks
        Restart%data(nb,num)%var3p => Coupling%qci_conv(Model%chunk_begin(nb):Model%chunk_end(nb),:)
      enddo
    endif

    !--- RAP/HRRR-specific variables, 3D
    ! GF
    if (Model%imfdeepcnv == Model%imfdeepcnv_gf) then
      num = num + 1
      Restart%name3d(num) = 'gf_3d_prevst'
      do nb = 1,nblks
        Restart%data(nb,num)%var3p => Tbd%prevst(Model%chunk_begin(nb):Model%chunk_end(nb),:)
      enddo
      num = num + 1
      Restart%name3d(num) = 'gf_3d_prevsq'
      do nb = 1,nblks
        Restart%data(nb,num)%var3p => Tbd%prevsq(Model%chunk_begin(nb):Model%chunk_end(nb),:)
      enddo
      num = num + 1
      Restart%name3d(num) = 'gf_3d_qci_conv'
      do nb = 1,nblks
        Restart%data(nb,num)%var3p => Coupling%qci_conv(Model%chunk_begin(nb):Model%chunk_end(nb),:)
      enddo
    endif
    ! MYNN PBL
    if (Model%do_mynnedmf) then
      num = num + 1
      Restart%name3d(num) = 'mynn_3d_cldfra_bl'
      do nb = 1,nblks
        Restart%data(nb,num)%var3p => Tbd%cldfra_bl(Model%chunk_begin(nb):Model%chunk_end(nb),:)
      enddo
      num = num + 1
      Restart%name3d(num) = 'mynn_3d_qc_bl'
      do nb = 1,nblks
        Restart%data(nb,num)%var3p => Tbd%qc_bl(Model%chunk_begin(nb):Model%chunk_end(nb),:)
      enddo
      num = num + 1
      Restart%name3d(num) = 'mynn_3d_qi_bl'
      do nb = 1,nblks
        Restart%data(nb,num)%var3p => Tbd%qi_bl(Model%chunk_begin(nb):Model%chunk_end(nb),:)
      enddo
      num = num + 1
      Restart%name3d(num) = 'mynn_3d_el_pbl'
      do nb = 1,nblks
        Restart%data(nb,num)%var3p => Tbd%el_pbl(Model%chunk_begin(nb):Model%chunk_end(nb),:)
      enddo
      num = num + 1
      Restart%name3d(num) = 'mynn_3d_sh3d'
      do nb = 1,nblks
        Restart%data(nb,num)%var3p => Tbd%sh3d(Model%chunk_begin(nb):Model%chunk_end(nb),:)
      enddo
      num = num + 1
      Restart%name3d(num) = 'mynn_3d_qke'
      do nb = 1,nblks
        Restart%data(nb,num)%var3p => Tbd%qke(Model%chunk_begin(nb):Model%chunk_end(nb),:)
      enddo
      num = num + 1
      Restart%name3d(num) = 'mynn_3d_tsq'
      do nb = 1,nblks
        Restart%data(nb,num)%var3p => Tbd%tsq(Model%chunk_begin(nb):Model%chunk_end(nb),:)
      enddo
      num = num + 1
      Restart%name3d(num) = 'mynn_3d_qsq'
      do nb = 1,nblks
        Restart%data(nb,num)%var3p => Tbd%qsq(Model%chunk_begin(nb):Model%chunk_end(nb),:)
      enddo
      num = num + 1
      Restart%name3d(num) = 'mynn_3d_cov'
      do nb = 1,nblks
        Restart%data(nb,num)%var3p => Tbd%cov(Model%chunk_begin(nb):Model%chunk_end(nb),:)
      enddo
    endif

    ! Radar-derived microphysics temperature tendencies
    if (Model%num_dfi_radar > 0) then
      do itime=1,Model%dfi_radar_max_intervals
        if(Model%ix_dfi_radar(itime)>0) then
          num = num + 1
          if(itime==1) then
            Restart%name3d(num) = 'radar_tten'
          else
            write(Restart%name3d(num),'("radar_tten_",I0)') itime
          endif
          do nb = 1,nblks
            Restart%data(nb,num)%var3p => Tbd%dfi_radar_tten( &
              Model%chunk_begin(nb):Model%chunk_end(nb),:,Model%ix_dfi_radar(itime))
          enddo
        endif
      enddo
    endif

    if(Model%rrfs_sd) then
      num = num + 1
      Restart%name3d(num) = 'chem3d_1'
      do nb = 1,nblks
        Restart%data(nb,num)%var3p => Coupling%chem3d(Model%chunk_begin(nb):Model%chunk_end(nb),:,1)
      enddo
      num = num + 1
      Restart%name3d(num) = 'chem3d_2'
      do nb = 1,nblks
        Restart%data(nb,num)%var3p => Coupling%chem3d(Model%chunk_begin(nb):Model%chunk_end(nb),:,2)
      enddo
      num = num + 1
      Restart%name3d(num) = 'chem3d_3'
      do nb = 1,nblks
        Restart%data(nb,num)%var3p => Coupling%chem3d(Model%chunk_begin(nb):Model%chunk_end(nb),:,3)
      enddo
      num = num + 1
      Restart%name3d(num) = 'ext550'
      do nb = 1,nblks
        Restart%data(nb,num)%var3p => Radtend%ext550(Model%chunk_begin(nb):Model%chunk_end(nb),:)
      enddo
    endif

  end subroutine GFS_restart_populate

end module GFS_restart
