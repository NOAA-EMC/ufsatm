!> \file fv3atm_restart_io.F90
!! This file contains the restart reading and writing code, for quilt and non-quilt
!! of the Sfcprop and physics data.

module fv3atm_restart_io_mod

  use block_control_mod,  only: block_control_type
  use mpp_mod,            only: mpp_error, mpp_chksum, NOTE,   FATAL
  use GFS_typedefs,       only: GFS_sfcprop_type, GFS_control_type, kind_phys, GFS_data_type
  use GFS_restart,        only: GFS_restart_type
  use fms_mod,            only: stdout
  use fms2_io_mod,        only: FmsNetcdfDomainFile_t, unlimited,      &
                                open_file, close_file,                 &
                                register_axis, register_restart_field, &
                                register_variable_attribute, register_field, &
                                read_restart, write_restart, write_data,     &
                                get_global_io_domain_indices, get_dimension_size
  use mpp_domains_mod,    only: domain2d
  use fv3atm_common_io,   only: create_2d_field_and_add_to_bundle, &
       create_3d_field_and_add_to_bundle, copy_from_gfs_data, axis_type
  use fv3atm_sfc_io
  use fv3atm_rrfs_sd_io
  use fv3atm_clm_lake_io
  use fv3atm_oro_io

  implicit none
  private

  public fv3atm_checksum
  public fv3atm_restart_read
  public fv3atm_restart_write
  public fv3atm_restart_register
  public fv_phy_restart_output
  public fv_phy_restart_bundle_setup
  public fv_sfc_restart_output
  public fv_sfc_restart_bundle_setup

  !>\defgroup fv3atm_restart_io_mod module
  !> @{

  !>@Internal storage for reading and writing physics restart files.
  type phy_data_type
    real(kind=kind_phys), pointer, dimension(:,:,:)   :: var2 => null()
    real(kind=kind_phys), pointer, dimension(:,:,:,:) :: var3 => null()
    character(len=32),dimension(:),pointer :: var2_names => null()
    character(len=32),dimension(:),pointer :: var3_names => null()
    integer :: nvar2d = 0, nvar3d = 0, npz = 0
  contains
    procedure :: alloc => phy_data_alloc
    procedure :: transfer_data => phy_data_transfer_data
    final phy_data_final
  end type phy_data_type

  !--- GFDL filenames

  !>@ Filename template for orography data. FMS may add grid and tile information to the name
  character(len=32), parameter  :: fn_oro    = 'oro_data.nc'

  !>@ Filename template for gravity wave drag large-scale orography data. FMS may add grid and tile information to the name
  character(len=32), parameter  :: fn_oro_ls = 'oro_data_ls.nc'

  !>@ Filename template for gravity wave drag small-scale orography data. FMS may add grid and tile information to the name
  character(len=32), parameter  :: fn_oro_ss = 'oro_data_ss.nc'

  !>@ Filename template for surface data that doesn't fall under other categories. FMS may add grid and tile information to the name
  character(len=32), parameter  :: fn_srf    = 'sfc_data.nc'

  !>@ Filename template for physics diagnostic data. FMS may add grid and tile information to the name
  character(len=32), parameter  :: fn_phy    = 'phy_data.nc'

  !>@ Filename template for monthly dust data for RRFS_SD. FMS may add grid and tile information to the name
  character(len=32), parameter  :: fn_dust12m= 'dust12m_data.nc'

  !>@ Filename template for RRFS-SD emissions data. FMS may add grid and tile information to the name
  character(len=32), parameter  :: fn_emi    = 'emi_data.nc'

  !>@ Filename template for RRFS-SD smoke data. FMS may add grid and tile information to the name
  character(len=32), parameter  :: fn_rrfssd = 'SMOKE_RRFS_data.nc'

  real(kind_phys), parameter:: zero = 0.0, one = 1.0

  !>@ Instance of phy_data_type for quilt output of physics diagnostic data
  type(phy_data_type) :: phy_quilt

  !>@ Instance of clm_lake_data_type for quilt output of CLM Lake model restart data
  type(clm_lake_data_type) :: clm_lake_quilt

  !>@ Instance of Sfc_io_data_type for quilt output of surface restart data
  type(Sfc_io_data_type) :: sfc_quilt

  !>@ Instance of rrfs_sd_state_type for quilt output of RRFS-SD scheme restart data
  type(rrfs_sd_state_type) :: rrfs_sd_quilt

contains

  !>@brief Reads physics and surface fields.
  !> \section fv3atm_restart_read subroutine
  !! Calls sfc_prop_restart_read and phys_restart_read to read all surface and physics restart files.
  subroutine fv3atm_restart_read (GFS_Data, GFS_Restart, Atm_block, Model, fv_domain, warm_start, ignore_rst_cksum)
    implicit none
    type(GFS_data_type),      intent(inout) :: GFS_Data(:)
    type(GFS_restart_type),   intent(inout) :: GFS_Restart
    type(block_control_type), intent(in)    :: Atm_block
    type(GFS_control_type),   intent(inout) :: Model
    type(domain2d),           intent(in)    :: fv_domain
    logical,                  intent(in)    :: warm_start
    logical,                  intent(in)    :: ignore_rst_cksum

    !--- read in surface data from chgres
    call sfc_prop_restart_read (GFS_Data%Sfcprop, Atm_block, Model, fv_domain, warm_start, ignore_rst_cksum)

    !--- read in physics restart data
    call phys_restart_read (GFS_Restart, Atm_block, Model, fv_domain, ignore_rst_cksum)

  end subroutine fv3atm_restart_read

  !>@brief Writes surface and physics restart fields without using the write component (quilt).
  !> \section fv3atm_restart_write subroutine
  !! Calls sfc_prop_restart_write and phys_restart_write to write
  !! surface and physics restart fields. This pauses the model to
  !! write; it does not use the write component (quilt).
  subroutine fv3atm_restart_write (GFS_Data, GFS_Restart, Atm_block, Model, fv_domain, timestamp)
    implicit none
    type(GFS_data_type),         intent(inout) :: GFS_Data(:)
    type(GFS_restart_type),      intent(inout) :: GFS_Restart
    type(block_control_type),    intent(in)    :: Atm_block
    type(GFS_control_type),      intent(in)    :: Model
    type(domain2d),              intent(in)    :: fv_domain
    character(len=32), optional, intent(in)    :: timestamp

    !--- write surface data from chgres
    call sfc_prop_restart_write (GFS_Data%Sfcprop, Atm_block, Model, fv_domain, timestamp)

    !--- write physics restart data
    call phys_restart_write (GFS_Restart, Atm_block, Model, fv_domain, timestamp)

  end subroutine fv3atm_restart_write

  !----------------
  ! fv3atm_checksum
  !----------------
  subroutine fv3atm_checksum (Model, GFS_Data, Atm_block)
    implicit none
    !--- interface variables
    type(GFS_control_type),    intent(in) :: Model
    type(GFS_data_type),       intent(in) :: GFS_Data(:)
    type (block_control_type), intent(in) :: Atm_block
    !--- local variables
    integer :: outunit, i, ix, nb, isc, iec, jsc, jec, lev, ntr, k
    integer :: nsfcprop2d, nt
    real(kind=kind_phys), allocatable :: temp2d(:,:,:)
    real(kind=kind_phys), allocatable :: temp3d(:,:,:,:)
    real(kind=kind_phys), allocatable :: temp3dlevsp1(:,:,:,:)
    integer, allocatable :: ii1(:), jj1(:)
    character(len=32) :: name

    isc = Model%isc
    iec = Model%isc+Model%nx-1
    jsc = Model%jsc
    jec = Model%jsc+Model%ny-1
    lev = Model%levs

    ntr = size(GFS_Data(1)%Statein%qgrs,3)

    nsfcprop2d = 94
    if (Model%lsm == Model%lsm_noahmp) then
      nsfcprop2d = nsfcprop2d + 49
      if (Model%use_cice_alb) then
        nsfcprop2d = nsfcprop2d + 4
      endif
    elseif (Model%lsm == Model%lsm_ruc) then
      nsfcprop2d = nsfcprop2d + 4 + 12
      if (Model%rdlai) then
        nsfcprop2d = nsfcprop2d + 1
      endif
    else
      if (Model%use_cice_alb) then
        nsfcprop2d = nsfcprop2d + 4
      endif
    endif

    if (Model%nstf_name(1) > 0) then
      nsfcprop2d = nsfcprop2d + 16
    endif

    if(Model%lkm>0 .and. Model%iopt_lake==Model%iopt_lake_flake) then
      nsfcprop2d = nsfcprop2d + 10
    endif

    allocate (temp2d(isc:iec,jsc:jec,nsfcprop2d+Model%ntot2d+Model%nctp))
    allocate (temp3d(isc:iec,jsc:jec,1:lev,14+Model%ntot3d+2*ntr))
    allocate (temp3dlevsp1(isc:iec,jsc:jec,1:lev+1,3))

    temp2d = zero
    temp3d = zero
    temp3dlevsp1 = zero

    !$omp parallel do default(shared) private(i, k, nb, ix, nt, ii1, jj1)
    block_loop: do nb = 1, Atm_block%nblks
      allocate(ii1(Atm_block%blksz(nb)))
      allocate(jj1(Atm_block%blksz(nb)))
      ii1=Atm_block%index(nb)%ii - isc + 1
      jj1=Atm_block%index(nb)%jj - jsc + 1

      ! Copy into temp2d
      nt=0

      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Statein%pgr)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%slmsk)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%tsfc)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%tisfc)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%snowd)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%zorl)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%fice)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%hprime(:,1))
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%sncovr)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%snoalb)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%alvsf)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%alnsf)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%alvwf)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%alnwf)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%facsf)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%facwf)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%slope)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%shdmin)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%shdmax)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%tg3)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%vfrac)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%vtype)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%stype)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%scolor)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%uustar)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%oro)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%oro_uf)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%hice)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%weasd)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%canopy)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%ffmm)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%ffhh)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%f10m)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%tprcp)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%srflag)
      lsm_choice: if (Model%lsm == Model%lsm_noah .or. Model%lsm == Model%lsm_noahmp) then
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%slc)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%smc)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%stc)
      elseif (Model%lsm == Model%lsm_ruc) then
        do k=1,3
          call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%sh2o(:,k))
        enddo
        ! Combine levels 4 to lsoil_lsm (9 for RUC) into one
        nt=nt+1
        do ix=1,Atm_block%blksz(nb)
          temp2d(ii1(ix),jj1(ix),nt) = sum(GFS_Data(nb)%Sfcprop%sh2o(ix,4:Model%lsoil_lsm))
        enddo
        do k=1,3
          call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%smois(:,k))
        enddo
        ! Combine levels 4 to lsoil_lsm (9 for RUC) into one
        nt=nt+1
        do ix=1,Atm_block%blksz(nb)
          temp2d(ii1(ix),jj1(ix),nt) = sum(GFS_Data(nb)%Sfcprop%smois(ix,4:Model%lsoil_lsm))
        enddo
        do k=1,3
          call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%tslb(:,k))
        enddo
        ! Combine levels 4 to lsoil_lsm (9 for RUC) into one
        nt=nt+1
        do ix=1,Atm_block%blksz(nb)
          temp2d(ii1(ix),jj1(ix),nt) = sum(GFS_Data(nb)%Sfcprop%tslb(ix,4:Model%lsoil_lsm))
        enddo
      endif lsm_choice

      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%t2m)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%q2m)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Coupling%nirbmdi)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Coupling%nirdfdi)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Coupling%visbmdi)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Coupling%visdfdi)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Coupling%nirbmui)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Coupling%nirdfui)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Coupling%visbmui)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Coupling%visdfui)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Coupling%sfcdsw)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Coupling%sfcnsw)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Coupling%sfcdlw)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Grid%xlon)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Grid%xlat)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Grid%xlat_d)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Grid%sinlat)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Grid%coslat)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Grid%area)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Grid%dx)
      if (Model%ntoz > 0) then
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Grid%ddy_o3)
      endif
      if (Model%h2o_phys) then
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Grid%ddy_h)
      endif
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Cldprop%cv)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Cldprop%cvt)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Cldprop%cvb)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Radtend%sfalb)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Radtend%coszen)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Radtend%tsflw)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Radtend%semis)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Radtend%coszdg)

      ! Radtend%sfcfsw is an array of derived type, so we copy all
      ! eight elements of the type in one loop
      do ix=1,Atm_block%blksz(nb)
        temp2d(ii1(ix),jj1(ix),nt+1) = GFS_Data(nb)%Radtend%sfcfsw(ix)%upfxc
        temp2d(ii1(ix),jj1(ix),nt+2) = GFS_Data(nb)%Radtend%sfcfsw(ix)%upfx0
        temp2d(ii1(ix),jj1(ix),nt+3) = GFS_Data(nb)%Radtend%sfcfsw(ix)%dnfxc
        temp2d(ii1(ix),jj1(ix),nt+4) = GFS_Data(nb)%Radtend%sfcfsw(ix)%dnfx0
        temp2d(ii1(ix),jj1(ix),nt+5) = GFS_Data(nb)%Radtend%sfcflw(ix)%upfxc
        temp2d(ii1(ix),jj1(ix),nt+6) = GFS_Data(nb)%Radtend%sfcflw(ix)%upfx0
        temp2d(ii1(ix),jj1(ix),nt+7) = GFS_Data(nb)%Radtend%sfcflw(ix)%dnfxc
        temp2d(ii1(ix),jj1(ix),nt+8) = GFS_Data(nb)%Radtend%sfcflw(ix)%dnfx0
      enddo
      nt = nt + 8

      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%tiice(:,1))
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%tiice(:,2))
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%albdirvis_lnd)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%albdirnir_lnd)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%albdifvis_lnd)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%albdifnir_lnd)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%emis_lnd)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%emis_ice)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%sncovr_ice)

      if (Model%use_cice_alb .or. Model%lsm == Model%lsm_ruc) then
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%albdirvis_ice)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%albdirnir_ice)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%albdifvis_ice)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%albdifnir_ice)
      endif

      lsm_choice_2: if (Model%lsm == Model%lsm_noahmp) then
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%snowxy)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%tvxy)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%tgxy)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%canicexy)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%canliqxy)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%eahxy)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%tahxy)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%cmxy)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%chxy)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%fwetxy)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%sneqvoxy)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%alboldxy)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%qsnowxy)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%wslakexy)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%zwtxy)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%waxy)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%wtxy)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%lfmassxy)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%rtmassxy)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%stmassxy)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%woodxy)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%stblcpxy)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%fastcpxy)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%xsaixy)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%xlaixy)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%taussxy)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%smcwtdxy)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%deeprechxy)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%rechxy)

        ! These five arrays use bizarre indexing, so we use loops:
        do k=-2,0
          call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%snicexy(:,k))
        enddo

        do k=-2,0
          call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%snliqxy(:,k))
        enddo

        do k=-2,0
          call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%tsnoxy(:,k))
        enddo

        do k=1,4
          call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%smoiseq(:,k))
        enddo

        do k=-2,4
          call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%zsnsoxy(:,k))
        enddo
      elseif (Model%lsm == Model%lsm_ruc) then
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%wetness)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%clw_surf_land)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%clw_surf_ice)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%qwv_surf_land)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%qwv_surf_ice)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%tsnow_land)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%tsnow_ice)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%snowfallac_land)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%snowfallac_ice)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%sfalb_lnd)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%sfalb_lnd_bck)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%sfalb_ice)
        if (Model%rdlai) then
          call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%xlaixy)
        endif
      endif lsm_choice_2

      nstf_name_choice: if (Model%nstf_name(1) > 0) then
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%tref)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%z_c)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%c_0)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%c_d)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%w_0)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%w_d)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%xt)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%xs)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%xu)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%xz)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%zm)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%xtts)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%xzts)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%ifd)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%dt_cool)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%qrain)
      endif nstf_name_choice

      ! Flake
      if (Model%lkm > 0 .and. Model%iopt_lake==Model%iopt_lake_flake) then
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%T_snow)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%T_ice)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%h_ML)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%t_ML)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%t_mnw)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%h_talb)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%t_talb)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%t_bot1)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%t_bot2)
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Sfcprop%c_t)
      endif

      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Tbd%phy_f2d)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp2d,GFS_Data(nb)%Tbd%phy_fctd)

      ! Copy to temp3dlevsp1
      nt=0

      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp3dlevsp1, GFS_Data(nb)%Statein%phii)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp3dlevsp1, GFS_Data(nb)%Statein%prsi)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp3dlevsp1, GFS_Data(nb)%Statein%prsik)

      ! Copy to temp3d
      nt=0

      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp3d,GFS_Data(nb)%Statein%phil)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp3d,GFS_Data(nb)%Statein%prsl)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp3d,GFS_Data(nb)%Statein%prslk)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp3d,GFS_Data(nb)%Statein%ugrs)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp3d,GFS_Data(nb)%Statein%vgrs)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp3d,GFS_Data(nb)%Statein%vvl)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp3d,GFS_Data(nb)%Statein%tgrs)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp3d,GFS_Data(nb)%Stateout%gu0)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp3d,GFS_Data(nb)%Stateout%gv0)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp3d,GFS_Data(nb)%Stateout%gt0)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp3d,GFS_Data(nb)%Radtend%htrsw)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp3d,GFS_Data(nb)%Radtend%htrlw)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp3d,GFS_Data(nb)%Radtend%swhc)
      call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp3d,GFS_Data(nb)%Radtend%lwhc)
      do k = 1,Model%ntot3d
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp3d,GFS_Data(nb)%Tbd%phy_f3d(:,:,k))
      enddo
      do k = 1,ntr
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp3d,GFS_Data(nb)%Statein%qgrs(:,:,k))
        call copy_from_GFS_Data(ii1,jj1,isc,jsc,nt,temp3d,GFS_Data(nb)%Stateout%gq0(:,:,k))
      enddo
    enddo block_loop


    outunit = stdout()
    do i = 1,nsfcprop2d+Model%ntot2d+Model%nctp
      write (name, '(i3.3,3x,4a)') i, ' 2d '
      write(outunit,100) name, mpp_chksum(temp2d(:,:,i:i))
    enddo
    do i = 1,3
      write (name, '(i2.2,3x,4a)') i, ' 3d levsp1'
      write(outunit,100) name, mpp_chksum(temp3dlevsp1(:,:,:,i:i))
    enddo
    do i = 1,14+Model%ntot3d+2*ntr
      write (name, '(i2.2,3x,4a)') i, ' 3d levs'
      write(outunit,100) name, mpp_chksum(temp3d(:,:,:,i:i))
    enddo
100 format("CHECKSUM::",A32," = ",Z20)

    deallocate(temp2d)
    deallocate(temp3d)
    deallocate(temp3dlevsp1)
  end subroutine fv3atm_checksum

  !>@brief Reads surface, orography, CLM Lake, and RRFS-SD data.
  !> \section sfc_prop_restart_read subroutine
  !!  Creates and populates a data type which is then used to "register"
  !!  restart variables with the FMS restart subsystem.
  !!  Calls an FMS routine to restore the data from a restart file.
  !!  Also calculates sncovr if it is not present in the restart file.
  subroutine sfc_prop_restart_read (Sfcprop, Atm_block, Model, fv_domain, warm_start, ignore_rst_cksum)
    use fv3atm_rrfs_sd_io
    implicit none
    !--- interface variable definitions
    type(GFS_sfcprop_type),    intent(inout) :: Sfcprop(:)
    type (block_control_type), intent(in)    :: Atm_block
    type(GFS_control_type),    intent(inout) :: Model
    type (domain2d),           intent(in)    :: fv_domain
    logical,                   intent(in)    :: warm_start
    logical,                   intent(in)    :: ignore_rst_cksum
    !--- directory of the input files
    character(5)  :: indir='INPUT'
    character(37) :: infile
    !--- fms2_io file open logic
    logical :: amiopen
    logical :: override_frac_grid

    type(clm_lake_data_type) :: clm_lake
    type(rrfs_sd_state_type) :: rrfs_sd_state
    type(rrfs_sd_emissions_type) :: rrfs_sd_emis
    type(Oro_scale_io_data_type) :: oro_ss
    type(Oro_scale_io_data_type) :: oro_ls
    type(Sfc_io_data_type) :: sfc
    type(Oro_io_data_type) :: oro

    type(FmsNetcdfDomainFile_t) :: Oro_restart, Sfc_restart, dust12m_restart, emi_restart, rrfssd_restart
    type(FmsNetcdfDomainFile_t) :: Oro_ls_restart, Oro_ss_restart

    !--- OROGRAPHY FILE

    !--- open file
    infile=trim(indir)//'/'//trim(fn_oro)
    amiopen=open_file(Oro_restart, trim(infile), 'read', domain=fv_domain, is_restart=.true., dont_add_res_to_filename=.true.)
    if (.not.amiopen) call mpp_error( FATAL, 'Error with opening file '//trim(infile) )

    call oro%register(Model,Oro_restart,Atm_block)

    !--- read the orography restart/data
    call mpp_error(NOTE,'reading topographic/orographic information from INPUT/oro_data.tile*.nc')
    call read_restart(Oro_restart, ignore_checksum=ignore_rst_cksum)
    call close_file(Oro_restart)

    !--- copy data into GFS containers
    call oro%copy(Model, Sfcprop, Atm_block)

    if_smoke: if(Model%rrfs_sd) then  ! for RRFS-SD

      !--- Dust input FILE
      !--- open file
      infile=trim(indir)//'/'//trim(fn_dust12m)
      amiopen=open_file(dust12m_restart, trim(infile), 'read', domain=fv_domain, is_restart=.true., dont_add_res_to_filename=.true.)
      if (.not.amiopen) call mpp_error( FATAL, 'Error with opening file'//trim(infile) )

      !--- Register axes and variables, allocate memory:
      call rrfs_sd_emis%register_dust12m(dust12m_restart, Atm_block)

      !--- read new GSL created dust12m restart/data
      call mpp_error(NOTE,'reading dust12m information from INPUT/dust12m_data.tile*.nc')
      call read_restart(dust12m_restart)
      call close_file(dust12m_restart)

      !--- Copy to Sfcprop and free temporary arrays:
      call rrfs_sd_emis%copy_dust12m(Sfcprop, Atm_block)

      !----------------------------------------------

      !--- open anthropogenic emission file
      infile=trim(indir)//'/'//trim(fn_emi)
      amiopen=open_file(emi_restart, trim(infile), 'read', domain=fv_domain, is_restart=.true., dont_add_res_to_filename=.true.)
      if (.not.amiopen) call mpp_error( FATAL, 'Error with opening file'//trim(infile) )

      ! Register axes and variables, allocate memory
      call rrfs_sd_emis%register_emi(emi_restart, Atm_block)

      !--- read anthropogenic emi restart/data
      call mpp_error(NOTE,'reading emi information from INPUT/emi_data.tile*.nc')
      call read_restart(emi_restart)
      call close_file(emi_restart)

      !--- Copy to Sfcprop and free temporary arrays:
      call rrfs_sd_emis%copy_emi(Sfcprop, Atm_block)

      !----------------------------------------------

      !--- Dust input FILE
      !--- open file
      infile=trim(indir)//'/'//trim(fn_rrfssd)
      amiopen=open_file(rrfssd_restart, trim(infile), 'read', domain=fv_domain, is_restart=.true., dont_add_res_to_filename=.true.)
      if (.not.amiopen) call mpp_error( FATAL, 'Error with opening file'//trim(infile) )

      ! Register axes and variables, allocate memory
      call rrfs_sd_emis%register_fire(Model, rrfssd_restart, Atm_block)

      !--- read new GSL created rrfssd restart/data
      call mpp_error(NOTE,'reading rrfssd information from INPUT/SMOKE_RRFS_data.nc')
      call read_restart(rrfssd_restart)
      call close_file(rrfssd_restart)

      !--- Copy to Sfcprop and free temporary arrays:
      call rrfs_sd_emis%copy_fire(Model, Sfcprop, Atm_block)

    endif if_smoke  ! RRFS_SD

    !--- Modify/read-in additional orographic static fields for GSL drag suite
    if (Model%gwd_opt==3 .or. Model%gwd_opt==33 .or. &
         Model%gwd_opt==2 .or. Model%gwd_opt==22 ) then

      if ( (Model%gwd_opt==3 .or. Model%gwd_opt==33) .or.    &
           ( (Model%gwd_opt==2 .or. Model%gwd_opt==22) .and. &
           Model%do_gsl_drag_ls_bl ) ) then
        !--- open restart file
        infile=trim(indir)//'/'//trim(fn_oro_ls)
        amiopen=open_file(Oro_ls_restart, trim(infile), 'read', domain=fv_domain, is_restart=.true., dont_add_res_to_filename=.true.)
        if( .not.amiopen ) call mpp_error( FATAL, 'Error with opening file '//trim(infile) )
        call oro_ls%register(Model,Oro_ls_restart,Atm_block)
        !--- read new GSL created orography restart/data
        call mpp_error(NOTE,'reading topographic/orographic information from &
             &INPUT/oro_data_ls.tile*.nc')
        call read_restart(Oro_ls_restart, ignore_checksum=ignore_rst_cksum)
        call close_file(Oro_ls_restart)
        call oro_ls%copy(Sfcprop,Atm_block,1)
      endif

      !--- open restart file
      infile=trim(indir)//'/'//trim(fn_oro_ss)
      amiopen=open_file(Oro_ss_restart, trim(infile), 'read', domain=fv_domain, is_restart=.true., dont_add_res_to_filename=.true.)
      if( .not.amiopen ) call mpp_error( FATAL, 'Error with opening file '//trim(infile) )
      call oro_ss%register(Model,Oro_ss_restart,Atm_block)
      call mpp_error(NOTE,'reading topographic/orographic information from &
           &INPUT/oro_data_ss.tile*.nc')
      call read_restart(Oro_ss_restart, ignore_checksum=ignore_rst_cksum)
      call close_file(Oro_ss_restart)
      call oro_ss%copy(Sfcprop,Atm_block,15)
    end if

    !--- SURFACE FILE

    !--- open file
    infile=trim(indir)//'/'//trim(fn_srf)
    amiopen=open_file(Sfc_restart, trim(infile), "read", domain=fv_domain, is_restart=.true., dont_add_res_to_filename=.true.)
    if( .not.amiopen ) call mpp_error(FATAL, 'Error opening file'//trim(infile))

    if(sfc%allocate_arrays(Model, Atm_block, .true., warm_start)) then
      call sfc%fill_2d_names(Model, warm_start)
      call sfc%register_axes(Model, Sfc_restart, .true., warm_start)

      ! Tell CLM Lake to allocate data, and register its axes and fields
      if(Model%lkm>0 .and. Model%iopt_lake==Model%iopt_lake_clm) then
        call clm_lake%allocate_data(Model)
        call clm_lake%fill_data(Model,Atm_block,Sfcprop)
        call clm_lake%copy_from_grid(Model,Atm_block,Sfcprop)
        call clm_lake%register_axes(Model, Sfc_restart)
        call clm_lake%register_fields(Sfc_restart)
      endif

      if(Model%rrfs_sd) then
        call rrfs_sd_state%allocate_data(Model)
        call rrfs_sd_state%fill_data(Model, Atm_block, Sfcprop)
        call rrfs_sd_state%register_axis(Model, Sfc_restart)
        call rrfs_sd_state%register_fields(Sfc_restart)
      endif

      call sfc%register_2d_fields(Model,Sfc_restart,.true.,warm_start)
    endif  ! if not allocated

    call sfc%fill_3d_names(Model,warm_start)
    call sfc%register_3d_fields(Model,Sfc_restart,.true.,warm_start)
    call sfc%init_fields(Model)

    !--- read the surface restart/data
    call mpp_error(NOTE,'reading surface properties data from INPUT/sfc_data.tile*.nc')
    call read_restart(Sfc_restart, ignore_checksum=ignore_rst_cksum)
    call close_file(Sfc_restart)

    ! Tell clm_lake to copy data to temporary arrays
    if(Model%lkm>0 .and. Model%iopt_lake==Model%iopt_lake_clm) then
      call clm_lake%copy_to_grid(Model,Atm_block,Sfcprop)
    endif

    if(Model%rrfs_sd) then
      call rrfs_sd_state%copy_to_grid(Model,Atm_block,Sfcprop)
    end if

    !   write(0,*)' stype read in min,max=',minval(sfc%var2(:,:,35)),maxval(sfc%var2(:,:,35)),' sfc%name2=',sfc%name2(35)
    !   write(0,*)' stype read in min,max=',minval(sfc%var2(:,:,18)),maxval(sfc%var2(:,:,18))
    !   write(0,*)' sfc%var2=',sfc%var2(:,:,12)

    !--- place the data into the block GFS containers
    override_frac_grid=Model%frac_grid
    call sfc%copy_to_grid(Model, Atm_block, Sfcprop, warm_start, override_frac_grid)
    Model%frac_grid=override_frac_grid

    call mpp_error(NOTE, 'gfs_driver:: - after put to container ')

    call sfc%apply_safeguards(Model, Atm_block, Sfcprop)

    ! A standard-compliant Fortran 2003 compiler will call clm_lake_final and rrfs_sd_final here.

  end subroutine sfc_prop_restart_read

  !>@brief Writes surface restart data without using the write component.
  !> \section sfc_prop_restart_write procedure
  !! Routine to write out GFS surface restarts via the FMS restart
  !! subsystem. Takes an optional argument to append timestamps for intermediate
  !! restarts.
  subroutine sfc_prop_restart_write (Sfcprop, Atm_block, Model, fv_domain, timestamp)
    use fv3atm_rrfs_sd_io
    implicit none
    !--- interface variable definitions
    type(GFS_sfcprop_type),      intent(in) :: Sfcprop(:)
    type(block_control_type),    intent(in) :: Atm_block
    type(GFS_control_type),      intent(in) :: Model
    type(domain2d),              intent(in) :: fv_domain
    character(len=32), optional, intent(in) :: timestamp
    !--- directory of the input files
    character(7)  :: indir='RESTART'
    character(72) :: infile
    !--- fms2_io file open logic
    logical :: amiopen
    !--- variables used for fms2_io register axis

    type(clm_lake_data_type), target :: clm_lake
    type(rrfs_sd_state_type) :: rrfs_sd_state
    type(Sfc_io_data_type) :: sfc
    type(FmsNetcdfDomainFile_t) :: Sfc_restart

    !--- set filename
    infile=trim(indir)//'/'//trim(fn_srf)
    if( present(timestamp) ) infile=trim(indir)//'/'//trim(timestamp)//'.'//trim(fn_srf)

    !--- register axis
    amiopen=open_file(Sfc_restart, trim(infile), 'overwrite', domain=fv_domain, is_restart=.true., dont_add_res_to_filename=.true.)
    if_amiopen: if( amiopen ) then
      call sfc%register_axes(Model, Sfc_restart, .false., .true.)
      call sfc%write_axes(Model, Sfc_restart)
    else
      call mpp_error(FATAL, 'Error in opening file'//trim(infile) )
    end if if_amiopen

    ! Tell clm_lake to allocate data, register its axes, and call write_data for each axis's variable
    if(Model%lkm>0 .and. Model%iopt_lake==Model%iopt_lake_clm) then
      call clm_lake%allocate_data(Model)
      call clm_lake%register_axes(Model, Sfc_restart)
      call clm_lake%write_axes(Model, Sfc_restart)
    endif

    if(Model%rrfs_sd) then
      call rrfs_sd_state%allocate_data(Model)
      call rrfs_sd_state%register_axis(Model,Sfc_restart)
      call rrfs_sd_state%write_axis(Model,Sfc_restart)
    end if

    if (sfc%allocate_arrays(Model, Atm_block, .false., .true.)) then
      call sfc%fill_2d_names(Model,.true.)
    end if

    if(Model%lkm>0 .and. Model%iopt_lake==Model%iopt_lake_clm) then
      ! Tell clm_lake to register all of its fields
      call clm_lake%register_fields(Sfc_restart)
    endif

    if(Model%rrfs_sd) then
      call rrfs_sd_state%register_fields(Sfc_restart)
    endif

    ! Register 2D surface property fields (except lake, smoke, and dust)
    call sfc%register_2d_fields(Model, Sfc_restart, .false., .true.)

    ! Determine list of 3D surface property fields names:
    call sfc%fill_3d_names(Model, .true.)

    ! Register 3D surface property fields (except lake, smoke, and dust)
    call sfc%register_3d_fields(Model, Sfc_restart, .false., .true.)

    ! Tell clm_lake to copy Sfcprop data to its internal temporary arrays.
    if(Model%lkm>0 .and. Model%iopt_lake==Model%iopt_lake_clm) then
      call clm_lake%copy_from_grid(Model,Atm_block,Sfcprop)
    endif

    if(Model%rrfs_sd) then
      call rrfs_sd_state%copy_from_grid(Model,Atm_block,Sfcprop)
    endif

    call sfc%copy_from_grid(Model, Atm_block, Sfcprop)

    call write_restart(Sfc_restart)
    call close_file(Sfc_restart)

    ! A standard-compliant Fortran 2003 compiler will call rrfs_sd_final and clm_lake_final here

  end subroutine sfc_prop_restart_write

  !>@brief Reads the physics restart data.
  !> \section phys_restart_read subroutine
  !! Creates and populates a data type which is then used to "register"
  !! restart variables with the GFDL FMS restart subsystem.
  !! Calls a GFDL FMS routine to restore the data from a restart file.
  subroutine phys_restart_read (GFS_Restart, Atm_block, Model, fv_domain, ignore_rst_cksum)
    implicit none
    !--- interface variable definitions
    type(GFS_restart_type),      intent(in) :: GFS_Restart
    type(block_control_type),    intent(in) :: Atm_block
    type(GFS_control_type),      intent(in) :: Model
    type(domain2d),              intent(in) :: fv_domain
    logical,                     intent(in) :: ignore_rst_cksum
    !--- local variables
    integer :: i, j, k, nb, ix, num
    integer :: isc, iec, jsc, jec, nx, ny
    character(len=64) :: fname
    real(kind=kind_phys), pointer, dimension(:,:)   :: var2_p => NULL()
    real(kind=kind_phys), pointer, dimension(:,:,:) :: var3_p => NULL()
    !--- directory of the input files
    character(5)  :: indir='INPUT'
    logical :: amiopen, was_allocated

    type(phy_data_type) :: phy
    type(FmsNetcdfDomainFile_t) :: Phy_restart

    isc = Atm_block%isc
    iec = Atm_block%iec
    jsc = Atm_block%jsc
    jec = Atm_block%jec
    nx  = (iec - isc + 1)
    ny  = (jec - jsc + 1)

    was_allocated = phy%alloc(GFS_Restart, Atm_block)

    !--- open restart file and register axes
    fname = trim(indir)//'/'//trim(fn_phy)
    amiopen=open_file(Phy_restart, trim(fname), 'read', domain=fv_domain, is_restart=.true., dont_add_res_to_filename=.true.)
    if( amiopen ) then
      call register_axis(Phy_restart, 'xaxis_1', 'X')
      call register_axis(Phy_restart, 'yaxis_1', 'Y')
      call register_axis(Phy_restart, 'zaxis_1', phy%npz)
      call register_axis(Phy_restart, 'Time', unlimited)
    else
      call mpp_error(NOTE,'No physics restarts - cold starting physical parameterizations')
      return
    endif

    !--- register the restart fields
    if(was_allocated) then

      do num = 1,phy%nvar2d
        var2_p => phy%var2(:,:,num)
        call register_restart_field(Phy_restart, trim(GFS_Restart%name2d(num)), var2_p, dimensions=(/'xaxis_1','yaxis_1','Time   '/),&
             &is_optional=.true.)
      enddo
      do num = 1,phy%nvar3d
        var3_p => phy%var3(:,:,:,num)
        call register_restart_field(Phy_restart, trim(GFS_restart%name3d(num)), var3_p, dimensions=(/'xaxis_1','yaxis_1','zaxis_1','Time   '/), is_optional=.true.)
      enddo
      nullify(var2_p)
      nullify(var3_p)
    endif

    !--- read the surface restart/data
    call mpp_error(NOTE,'reading physics restart data from INPUT/phy_data.tile*.nc')
    call read_restart(Phy_restart, ignore_checksum=ignore_rst_cksum)
    call close_file(Phy_restart)

    call phy%transfer_data(.true., GFS_Restart, Atm_block, Model)

  end subroutine phys_restart_read

  !>@brief Writes the physics restart file without using the write component
  !> \section phys_restart_write subroutine
  !! Routine to write out GFS surface restarts via the FMS restart
  !! subsystem. Takes an optional argument to append timestamps for intermediate
  !! restarts.
  subroutine phys_restart_write (GFS_Restart, Atm_block, Model, fv_domain, timestamp)
    implicit none
    !--- interface variable definitions
    type(GFS_restart_type),      intent(in) :: GFS_Restart
    type(block_control_type),    intent(in) :: Atm_block
    type(GFS_control_type),      intent(in) :: Model
    type(domain2d),              intent(in) :: fv_domain
    character(len=32), optional, intent(in) :: timestamp
    !--- local variables
    integer :: i, j, k, nb, ix, num
    integer :: isc, iec, jsc, jec, nx, ny
    real(kind=kind_phys), pointer, dimension(:,:)   :: var2_p => NULL()
    real(kind=kind_phys), pointer, dimension(:,:,:) :: var3_p => NULL()
    !--- used for axis data for fms2_io
    integer :: is, ie
    integer, allocatable, dimension(:) :: buffer
    character(7) :: indir='RESTART'
    character(72) :: infile
    logical :: amiopen, allocated_something
    integer :: xaxis_1_chunk, yaxis_1_chunk

    type(phy_data_type) :: phy
    type(FmsNetcdfDomainFile_t) :: Phy_restart

    isc = Atm_block%isc
    iec = Atm_block%iec
    jsc = Atm_block%jsc
    jec = Atm_block%jec
    nx  = (iec - isc + 1)
    ny  = (jec - jsc + 1)

    !--- register the restart fields
    allocated_something = phy%alloc(GFS_Restart, Atm_block)

    !--- set file name
    infile=trim(indir)//'/'//trim(fn_phy)
    if( present(timestamp) ) infile=trim(indir)//'/'//trim(timestamp)//'.'//trim(fn_phy)
    !--- register axis
    amiopen=open_file(Phy_restart, trim(infile), 'overwrite', domain=fv_domain, is_restart=.true., dont_add_res_to_filename=.true.)
    if( amiopen ) then
      call register_axis(Phy_restart, 'xaxis_1', 'X')
      call register_field(Phy_restart, 'xaxis_1', axis_type, (/'xaxis_1'/))
      call register_variable_attribute(Phy_restart, 'xaxis_1', 'cartesian_axis', 'X', str_len=1)
      call get_global_io_domain_indices(Phy_restart, 'xaxis_1', is, ie, indices=buffer)
      call write_data(Phy_restart, "xaxis_1", buffer)
      deallocate(buffer)
      call get_dimension_size(Phy_restart, 'xaxis_1', xaxis_1_chunk)

      call register_axis(Phy_restart, 'yaxis_1', 'Y')
      call register_field(Phy_restart, 'yaxis_1', axis_type, (/'yaxis_1'/))
      call register_variable_attribute(Phy_restart, 'yaxis_1', 'cartesian_axis', 'Y', str_len=1)
      call get_global_io_domain_indices(Phy_restart, 'yaxis_1', is, ie, indices=buffer)
      call write_data(Phy_restart, "yaxis_1", buffer)
      deallocate(buffer)
      call get_dimension_size(Phy_restart, 'yaxis_1', yaxis_1_chunk)

      call register_axis(Phy_restart, 'zaxis_1', phy%npz)
      call register_field(Phy_restart, 'zaxis_1', axis_type, (/'zaxis_1'/))
      call register_variable_attribute(Phy_restart, 'zaxis_1', 'cartesian_axis', 'Z', str_len=1)
      allocate( buffer(phy%npz) )
      do i=1, phy%npz
        buffer(i)=i
      end do
      call write_data(Phy_restart, "zaxis_1", buffer)
      deallocate(buffer)

      call register_axis(Phy_restart, 'Time', unlimited)
      call register_field(Phy_restart, 'Time', axis_type, (/'Time'/))
      call register_variable_attribute(Phy_restart, 'Time', 'cartesian_axis', 'T', str_len=1)
      call write_data(Phy_restart, "Time", 1)
    else
      call mpp_error(FATAL, 'Error opening file '//trim(infile))
    end if

    do num = 1,phy%nvar2d
      var2_p => phy%var2(:,:,num)
      call register_restart_field(Phy_restart, trim(GFS_Restart%name2d(num)), var2_p, dimensions=(/'xaxis_1','yaxis_1','Time   '/),&
           & chunksizes=(/xaxis_1_chunk,yaxis_1_chunk,1/), is_optional=.true.)
    enddo
    do num = 1,phy%nvar3d
      var3_p => phy%var3(:,:,:,num)
      call register_restart_field(Phy_restart, trim(GFS_Restart%name3d(num)), var3_p, dimensions=(/'xaxis_1','yaxis_1','zaxis_1','Time   '/),&
           & chunksizes=(/xaxis_1_chunk,yaxis_1_chunk,1,1/), is_optional=.true.)
    enddo
    nullify(var2_p)
    nullify(var3_p)

    call phy%transfer_data(.false., GFS_Restart, Atm_block, Model)

    call write_restart(Phy_restart)
    call close_file(Phy_restart)

  end subroutine phys_restart_write

  !>@brief Allocates buffers and registers fields for a quilting (write component) restart.
  !> \section fv3atm_restart_register subroutine
  !! Allocates all data buffers and sets variable names for surface and physics restarts.
  subroutine fv3atm_restart_register (Sfcprop, GFS_restart, Atm_block, Model)
    implicit none

    type(GFS_sfcprop_type),      intent(in) :: Sfcprop(:)
    type(GFS_restart_type),      intent(in) :: GFS_Restart
    type(block_control_type),    intent(in) :: Atm_block
    type(GFS_control_type),      intent(in) :: Model

    logical was_changed

    !--------------- phy
    was_changed = phy_quilt%alloc(GFS_Restart, Atm_block)

    !--------------- sfc
    was_changed = sfc_quilt%allocate_arrays(Model, Atm_block, .false., .true.)
    call sfc_quilt%fill_2d_names(Model, .true.)
    call sfc_quilt%fill_3d_names(Model, .true.)

    if(Model%iopt_lake == 2 .and. Model%lkm > 0) then
      call clm_lake_quilt%allocate_data(Model)
      call clm_lake_quilt%fill_data(Model, Atm_block, Sfcprop)
    endif

    if(Model%rrfs_sd) then
      call rrfs_sd_quilt%allocate_data(Model)
      call rrfs_sd_quilt%fill_data(Model, Atm_block, Sfcprop)
    endif

  end subroutine fv3atm_restart_register

  !>@Copies physics restart fields from write component data structures to the model grid.
  subroutine fv_phy_restart_output(GFS_Restart, Atm_block)

    implicit none

    type(GFS_restart_type),      intent(in) :: GFS_Restart
    type(block_control_type),    intent(in) :: Atm_block

    call phy_quilt%transfer_data(.false., GFS_Restart, Atm_block)

  end subroutine fv_phy_restart_output

  !>@Copies physics restart fields from the model grid to write component data structures
  subroutine fv_sfc_restart_output(Sfcprop, Atm_block, Model)
    !--- interface variable definitions
    implicit none

    type(GFS_sfcprop_type),      intent(in) :: Sfcprop(:)
    type(block_control_type),    intent(in) :: Atm_block
    type(GFS_control_type),      intent(in) :: Model

    call sfc_quilt%copy_from_grid(Model, Atm_block, Sfcprop)
    if(Model%iopt_lake == 2 .and. Model%lkm > 0) then
      call clm_lake_quilt%copy_from_grid(Model, Atm_block, Sfcprop)
    endif
    if(Model%rrfs_sd) then
      call rrfs_sd_quilt%copy_from_grid(Model, Atm_block, Sfcprop)
    endif

  end subroutine fv_sfc_restart_output

  !>@ Creates the ESMF bundle for physics restart data
  subroutine fv_phy_restart_bundle_setup(bundle, grid, rc)
    use esmf

    implicit none

    type(ESMF_FieldBundle),intent(inout)        :: bundle
    type(ESMF_Grid),intent(inout)               :: grid
    integer,intent(out)                         :: rc

    !*** local variables
    integer i
    character(128)    :: bdl_name
    character(128)    :: outputfile
    real(kind_phys),dimension(:,:),pointer   :: temp_r2d
    real(kind_phys),dimension(:,:,:),pointer   :: temp_r3d
    integer :: num
    real(kind_phys), allocatable :: axis_values(:)

    if (.not. associated(phy_quilt%var2)) then
      write(0,*)'ERROR phy_quilt%var2, NOT allocated'
    endif
    if (.not. associated(phy_quilt%var3)) then
      write(0,*)'ERROR phy_quilt%var3 NOT allocated'
    endif

    call ESMF_FieldBundleGet(bundle, name=bdl_name,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    outputfile = trim(bdl_name)

    !*** add esmf fields

    do num = 1,phy_quilt%nvar2d
      temp_r2d => phy_quilt%var2(:,:,num)
      call create_2d_field_and_add_to_bundle(temp_r2d, trim(phy_quilt%var2_names(num)), trim(outputfile), grid, bundle)
    enddo

    allocate(axis_values(phy_quilt%npz))
    axis_values = (/ (i, i=1,phy_quilt%npz) /)

    do num = 1,phy_quilt%nvar3d
      temp_r3d => phy_quilt%var3(:,:,:,num)
      call create_3d_field_and_add_to_bundle(temp_r3d, trim(phy_quilt%var3_names(num)), "zaxis_1", axis_values, trim(outputfile), grid, bundle)
    enddo

    deallocate(axis_values)

  end subroutine fv_phy_restart_bundle_setup

  !>@ Creates the ESMF bundle for surface restart data
  subroutine fv_sfc_restart_bundle_setup(bundle, grid, Model, rc)
    use esmf

    implicit none

    type(ESMF_FieldBundle),intent(inout)        :: bundle
    type(ESMF_Grid),intent(inout)               :: grid
    type(GFS_control_type),          intent(in) :: Model
    integer,intent(out)                         :: rc

    !*** local variables
    character(128)    :: sfcbdl_name
    character(128)    :: outputfile

    call ESMF_FieldBundleGet(bundle, name=sfcbdl_name,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    outputfile = trim(sfcbdl_name)

    !*** add esmf fields

    call sfc_quilt%bundle_2d_fields(bundle, grid, Model, outputfile)
    call sfc_quilt%bundle_3d_fields(bundle, grid, Model, outputfile)

    if(Model%iopt_lake == 2 .and. Model%lkm > 0) then
      call clm_lake_quilt%bundle_fields(bundle, grid, Model, outputfile)
    endif
    if(Model%rrfs_sd) then
      call rrfs_sd_quilt%bundle_fields(bundle, grid, Model, outputfile)
    endif

  end subroutine fv_sfc_restart_bundle_setup

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !
  !                     PRIVATE SUBROUTINES
  !
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !>@brief Allocates and fills internal data structures for quilt or non-quilt physics restart I/O
  !> \section phy_data_type%alloc procedure
  !! Allocates the variable and variable name data structures in the phy_data_type.
  !! Also, copies the GFS_Restart names to the phy_data_type arrays.
  !! Do not call from outside this module; it is part of the internal implementation.
  logical function phy_data_alloc(phy, GFS_Restart, Atm_block)
    use fv3atm_common_io, only: get_nx_ny_from_atm
    implicit none
    class(phy_data_type) :: phy
    type(GFS_restart_type),      intent(in) :: GFS_Restart
    type(block_control_type),    intent(in) :: Atm_block

    integer :: nx, ny, num

    phy_data_alloc = .false.

    if(associated(phy%var2)) return

    call get_nx_ny_from_atm(Atm_block, nx, ny)

    phy%npz = Atm_block%npz
    phy%nvar2d = GFS_Restart%num2d
    phy%nvar3d = GFS_Restart%num3d

    allocate (phy%var2(nx,ny,phy%nvar2d), phy%var2_names(phy%nvar2d))
    allocate (phy%var3(nx,ny,phy%npz,phy%nvar3d), phy%var3_names(phy%nvar3d))
    phy%var2 = zero
    phy%var3 = zero
    do num = 1,phy%nvar2d
      phy%var2_names(num) = trim(GFS_Restart%name2d(num))
    enddo
    do num = 1,phy%nvar3d
      phy%var3_names(num) = trim(GFS_Restart%name3d(num))
    enddo

    phy_data_alloc = .true.
  end function phy_data_alloc

  !>@brief Copies data between the internal physics restart data structures and the model grid
  !> \section phy_data_type%transfer_data procedure
  !! Restart I/O stores data in temporary arrays while interfacing with ESMF or FMS. This procedure
  !! copies between the temporary arrays and the model grid. The "reading" flag controls the
  !! direction of the copy. For reading=.true., data is copied from the temporary arrays to the
  !! model grid (during restart read). For reading=.false., data is copied from the model grid to
  !! temporary arrays (for writing the restart).
  subroutine phy_data_transfer_data(phy, reading, GFS_Restart, Atm_block, Model)
    use mpp_mod,            only: FATAL, mpp_error
    implicit none
    class(phy_data_type) :: phy
    logical, intent(in) :: reading
    type(GFS_restart_type) :: GFS_Restart
    type(block_control_type) :: Atm_block
    type(GFS_control_type), optional, intent(in) :: Model

    integer :: i, j, k, num, nb, ix

    !--- register the restart fields
    if (.not. associated(phy%var2)) then
      call mpp_error(FATAL,'phy%var2 must be allocated')
      return ! should never get here
    endif
    if (.not. associated(phy%var3)) then
      call mpp_error(FATAL,'phy%var3 must be allocated')
      return ! should never get here
    endif

    ! Copy 2D Vars

    if(reading) then
      !--- place the data into the block GFS containers
      !--- phy%var* variables
      do num = 1,phy%nvar2d
        !$omp parallel do default(shared) private(i, j, nb, ix)
        do nb = 1,Atm_block%nblks
          do ix = 1, Atm_block%blksz(nb)
            i = Atm_block%index(nb)%ii(ix) - Atm_block%isc + 1
            j = Atm_block%index(nb)%jj(ix) - Atm_block%jsc + 1
            GFS_Restart%data(nb,num)%var2p(ix) = phy%var2(i,j,num)
          enddo
        enddo
      enddo
    else
      !--- 2D variables
      do num = 1,phy%nvar2d
        !$omp parallel do default(shared) private(i, j, nb, ix)
        do nb = 1,Atm_block%nblks
          do ix = 1, Atm_block%blksz(nb)
            i = Atm_block%index(nb)%ii(ix) - Atm_block%isc + 1
            j = Atm_block%index(nb)%jj(ix) - Atm_block%jsc + 1
            phy%var2(i,j,num) = GFS_Restart%data(nb,num)%var2p(ix)
          enddo
        enddo
      enddo
    endif

    !-- if restart from init time, reset accumulated diag fields

    if(reading .and. present(Model)) then
      if(Model%phour < 1.e-7) then
        do num = GFS_Restart%fdiag,GFS_Restart%ldiag
          !$omp parallel do default(shared) private(i, j, nb, ix)
          do nb = 1,Atm_block%nblks
            do ix = 1, Atm_block%blksz(nb)
              i = Atm_block%index(nb)%ii(ix) - Atm_block%isc + 1
              j = Atm_block%index(nb)%jj(ix) - Atm_block%jsc + 1
              GFS_Restart%data(nb,num)%var2p(ix) = zero
            enddo
          enddo
        enddo
      endif
    endif

    ! Copy 3D Vars

    if(reading) then
      do num = 1,phy%nvar3d
        !$omp parallel do default(shared) private(i, j, k, nb, ix)
        do nb = 1,Atm_block%nblks
          do k=1,phy%npz
            do ix = 1, Atm_block%blksz(nb)
              i = Atm_block%index(nb)%ii(ix) - Atm_block%isc + 1
              j = Atm_block%index(nb)%jj(ix) - Atm_block%jsc + 1
              GFS_Restart%data(nb,num)%var3p(ix,k) = phy%var3(i,j,k,num)
            enddo
          enddo
        enddo
      enddo
    else
      !--- 3D variables
      do num = 1,phy%nvar3d
        !$omp parallel do default(shared) private(i, j, k, nb, ix)
        do nb = 1,Atm_block%nblks
          do k=1,phy%npz
            do ix = 1, Atm_block%blksz(nb)
              i = Atm_block%index(nb)%ii(ix) - Atm_block%isc + 1
              j = Atm_block%index(nb)%jj(ix) - Atm_block%jsc + 1
              phy%var3(i,j,k,num) = GFS_Restart%data(nb,num)%var3p(ix,k)
            enddo
          enddo
        enddo
      enddo
    endif

  end subroutine phy_data_transfer_data

  !>@ Destructor for phy_data_type
  subroutine phy_data_final(phy)
    implicit none
    type(phy_data_type) :: phy

    ! This #define reduces code length by a lot
#define IF_ASSOC_DEALLOC_NULL(var) \
    if(associated(phy%var)) then ; \
      deallocate(phy%var) ; \
      nullify(phy%var) ; \
    endif

    IF_ASSOC_DEALLOC_NULL(var2)
    IF_ASSOC_DEALLOC_NULL(var3)
    IF_ASSOC_DEALLOC_NULL(var2_names)
    IF_ASSOC_DEALLOC_NULL(var3_names)

#undef IF_ASSOC_DEALLOC_NULL
  end subroutine phy_data_final

end module fv3atm_restart_io_mod
!> @}
