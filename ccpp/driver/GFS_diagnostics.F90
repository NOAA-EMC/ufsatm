module GFS_diagnostics

!-----------------------------------------------------------------------
!    GFS_diagnostics_mod defines a data type and contains the routine
!    to populate said type with diagnostics from the GFS physics for
!    use by the modeling system for output
!-----------------------------------------------------------------------

  use machine,            only: kind_phys

  !--- GFS_typedefs ---
  use GFS_typedefs,       only: GFS_control_type,  GFS_statein_type,  &
                                GFS_stateout_type, GFS_sfcprop_type,  &
                                GFS_coupling_type, GFS_grid_type,     &
                                GFS_tbd_type,      GFS_cldprop_type,  &
                                GFS_radtend_type,  GFS_diag_type,     &
                                GFS_init_type
  implicit none
  private

  !--- private data type definition ---
  type data_subtype
    integer,              dimension(:),   pointer :: int2  => NULL()
    real(kind=kind_phys), dimension(:),   pointer :: var2  => NULL()
    real(kind=kind_phys), dimension(:),   pointer :: var21 => NULL()
    real(kind=kind_phys), dimension(:,:), pointer :: var3  => NULL()
  end type data_subtype

  !--- data type definition for use with GFDL FMS diagnostic manager until write component is working
  type GFS_externaldiag_type
    integer :: id
    integer :: axes
    logical :: time_avg
    character(len=64)    :: time_avg_kind
    character(len=64)    :: mod_name
    character(len=64)    :: name
    character(len=128)   :: desc
    character(len=64)    :: unit
    character(len=64)    :: mask
    character(len=64)    :: intpl_method
    real(kind=kind_phys) :: cnvfac
    type(data_subtype), dimension(:), allocatable :: data
   end type GFS_externaldiag_type

  !--- public data type ---
  public  GFS_externaldiag_type

  !--- public interfaces ---
  public  GFS_externaldiag_populate

  CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 ! Helper function for GFS_externaldiag_populate to handle the massive dtend(:,:,dtidx(:,:)) array
    subroutine add_dtend(Model,ExtDiag,IntDiag,idx,nblks,itrac,iprocess,desc,unit)
    implicit none
    type(GFS_control_type),       intent(in)    :: Model
    type(GFS_externaldiag_type),  intent(inout) :: ExtDiag(:)
    type(GFS_diag_type),          intent(in)    :: IntDiag
    integer, intent(in) :: nblks, itrac, iprocess
    integer, intent(inout) :: idx
    real(kind=kind_phys), pointer :: dtend(:,:,:) ! Assumption: dtend is null iff all(dtidx <= 1)
    character(len=*), intent(in), optional :: desc, unit

    integer :: idtend, nb

    idtend = Model%dtidx(itrac,iprocess)
    if(idtend>=1) then
       idx = idx + 1
       ExtDiag(idx)%axes = 3
       ExtDiag(idx)%name = 'dtend_'//trim(Model%dtend_var_labels(itrac)%name)//'_'//trim(Model%dtend_process_labels(iprocess)%name)
       ExtDiag(idx)%mod_name = Model%dtend_process_labels(iprocess)%mod_name
       ExtDiag(idx)%time_avg = Model%dtend_process_labels(iprocess)%time_avg
       if(present(desc)) then
          ExtDiag(idx)%desc = desc
       else
          ExtDiag(idx)%desc = trim(Model%dtend_var_labels(itrac)%desc)//' '//trim(Model%dtend_process_labels(iprocess)%desc)
       endif
       if(present(unit)) then
          ExtDiag(idx)%unit = trim(unit)
       else
          ExtDiag(idx)%unit = trim(Model%dtend_var_labels(itrac)%unit)
       endif
       allocate (ExtDiag(idx)%data(nblks))
       do nb = 1,nblks
          ExtDiag(idx)%data(nb)%var3 => IntDiag%dtend(Model%chunk_begin(nb):Model%chunk_end(nb),:,idtend)
       enddo
    endif
  end subroutine add_dtend

!-------------------------------------------------------------------------
!--- GFS_externaldiag_populate ---
!-------------------------------------------------------------------------
!    creates and populates a data type with GFS physics diagnostic
!    variables which is then handed off to the IPD for use by the model
!    infrastructure layer to output as needed.  The data type includes
!    names, units, conversion factors, etc.  There is no copying of data,
!    but instead pointers are associated to the internal representation
!    of each individual physics diagnostic.
!-------------------------------------------------------------------------
  subroutine GFS_externaldiag_populate (ExtDiag, Model, Statein, Stateout, Sfcprop, Coupling,  &
                                        Grid, Tbd, Cldprop, Radtend, IntDiag, Init_parm)
!---------------------------------------------------------------------------------------------!
!   DIAGNOSTIC_METADATA                                                                       !
!     ExtDiag%id                   [integer ]   switch to turn on/off variable output         !
!     ExtDiag%axes                 [integer ]   dimensionality of variable (2 or 3)           !
!     ExtDiag%time_avg             [logical ]   bucketed accumulation time average            !
!     ExtDiag%time_avg_kind        [char*64 ]   time average period                           !
!     ExtDiag%mod_name             [char*64 ]   classification of the variable                !
!     ExtDiag%name                 [char*64 ]   output name for variable                      !
!     ExtDiag%desc                 [char*128]   long description of field                     !
!     ExtDiag%unit                 [char*64 ]   units associated with field                   !
!     ExtDiag%mask                 [char*64 ]   description of mask-type                      !
!     ExtDiag%intpl_method         [char*64 ]   method to use for interpolation               !
!     ExtDiag%cnvfac               [real*8  ]   conversion factor to output specified units   !
!     ExtDiag%data(nb)%int2(:)     [integer ]   pointer to 2D data [=> null() for a 3D field] !
!     ExtDiag%data(nb)%var2(:)     [real*8  ]   pointer to 2D data [=> null() for a 3D field] !
!     ExtDiag%data(nb)%var21(:)    [real*8  ]   pointer to 2D data for ratios                 !
!     ExtDiag%data(nb)%var3(:,:)   [real*8  ]   pointer to 3D data [=> null() for a 2D field] !
!---------------------------------------------------------------------------------------------!

    use parse_tracers,    only: get_tracer_index
    implicit none
!
!  ---  interface variables
    type(GFS_externaldiag_type),  intent(inout) :: ExtDiag(:)
    type(GFS_control_type),       intent(in)    :: Model
    type(GFS_statein_type),       intent(in)    :: Statein
    type(GFS_stateout_type),      intent(in)    :: Stateout
    type(GFS_sfcprop_type),       intent(in)    :: Sfcprop
    type(GFS_coupling_type),      intent(in)    :: Coupling
    type(GFS_grid_type),          intent(in)    :: Grid
    type(GFS_tbd_type),           intent(in)    :: Tbd
    type(GFS_cldprop_type),       intent(in)    :: Cldprop
    type(GFS_radtend_type),       intent(in)    :: Radtend
    type(GFS_diag_type),          intent(in)    :: IntDiag
    type(GFS_init_type),          intent(in)    :: Init_parm

!--- local variables
    integer :: idt, idx, num, nb, nblks, NFXR, idtend, ichem, itrac, iprocess, i
    character(len=2) :: xtra
    real(kind=kind_phys), parameter :: cn_one = 1._kind_phys
    real(kind=kind_phys), parameter :: cn_100 = 100._kind_phys
    real(kind=kind_phys), parameter :: cn_th  = 1000._kind_phys
    real(kind=kind_phys), parameter :: cn_hr  = 3600._kind_phys
    character(len=30) :: namestr, descstr

    NFXR = Model%NFXR
    nblks = Model%nchunks

    ExtDiag(:)%id = -99
    ExtDiag(:)%axes = -99
    ExtDiag(:)%cnvfac = cn_one
    ExtDiag(:)%time_avg = .FALSE.
    ExtDiag(:)%time_avg_kind = ''
    ExtDiag(:)%mask = ''
    ExtDiag(:)%name = ''
    ExtDiag(:)%intpl_method = 'nearest_stod'

    idx = 0

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'cldfra2d'
    ExtDiag(idx)%desc = 'instantaneous 2D (max-in-column) cloud fraction'
    ExtDiag(idx)%unit = 'frac'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%cldfra2d(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'total_albedo'
    ExtDiag(idx)%desc = 'total sky albedo at toa'
    ExtDiag(idx)%unit = 'frac'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%total_albedo(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'lwp_ex'
    ExtDiag(idx)%desc = 'total liquid water path from explicit microphysics'
    ExtDiag(idx)%unit = 'kg m-2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%lwp_ex(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'iwp_ex'
    ExtDiag(idx)%desc = 'total ice water path from explicit microphysics'
    ExtDiag(idx)%unit = 'kg m-2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%iwp_ex(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'lwp_fc'
    ExtDiag(idx)%desc = 'total liquid water path from cloud fraction scheme'
    ExtDiag(idx)%unit = 'kg m-2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%lwp_fc(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'iwp_fc'
    ExtDiag(idx)%desc = 'total ice water path from cloud fraction scheme'
    ExtDiag(idx)%unit = 'kg m-2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%iwp_fc(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'ALBDO_ave'
    ExtDiag(idx)%desc = 'surface albedo'
    ExtDiag(idx)%unit = '%'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%cnvfac = cn_100
    ExtDiag(idx)%mask = 'positive_flux'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2  => IntDiag%fluxr(Model%chunk_begin(nb):Model%chunk_end(nb),3)
      ExtDiag(idx)%data(nb)%var21 => IntDiag%fluxr(Model%chunk_begin(nb):Model%chunk_end(nb),4)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'DLWRF'
    ExtDiag(idx)%desc = 'surface downward longwave flux'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%cnvfac = cn_one
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%dlwsfc(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'DLWRFI'
    ExtDiag(idx)%desc = 'instantaneous surface downward longwave flux'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%dlwsfci(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'ULWRF'
    ExtDiag(idx)%desc = 'surface upward longwave flux'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%cnvfac = cn_one
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%ulwsfc(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'DSWRFItoa'
    ExtDiag(idx)%desc = 'instantaneous top of atmos downward shortwave flux'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%cnvfac = cn_one
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%fluxr(Model%chunk_begin(nb):Model%chunk_end(nb),23)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'USWRFItoa'
    ExtDiag(idx)%desc = 'instantaneous top of atmos upward shortwave flux'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%cnvfac = cn_one
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%fluxr(Model%chunk_begin(nb):Model%chunk_end(nb),2)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'ULWRFItoa'
    ExtDiag(idx)%desc = 'instantaneous top of atmos upward longwave flux'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%cnvfac = cn_one
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%fluxr(Model%chunk_begin(nb):Model%chunk_end(nb),1)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'ULWRFI'
    ExtDiag(idx)%desc = 'instantaneous surface upward longwave flux'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%ulwsfci(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'DSWRF'
    ExtDiag(idx)%desc = 'averaged surface downward shortwave flux'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%cnvfac = cn_one
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%time_avg_kind = 'rad_sw'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%fluxr(Model%chunk_begin(nb):Model%chunk_end(nb),4)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'DSWRFI'
    ExtDiag(idx)%desc = 'instantaneous surface downward shortwave flux'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%dswsfci(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'USWRF'
    ExtDiag(idx)%desc = 'averaged surface upward shortwave flux'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%cnvfac = cn_one
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%time_avg_kind = 'rad_sw'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%fluxr(Model%chunk_begin(nb):Model%chunk_end(nb),3)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'USWRFI'
    ExtDiag(idx)%desc = 'instantaneous surface upward shortwave flux'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%uswsfci(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'duvb_ave'
    ExtDiag(idx)%desc = 'UV-B Downward Solar Flux'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%time_avg_kind = 'rad_sw'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%fluxr(Model%chunk_begin(nb):Model%chunk_end(nb),21)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'cduvb_ave'
    ExtDiag(idx)%desc = 'Clear sky UV-B Downward Solar Flux'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%time_avg_kind = 'rad_sw'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%fluxr(Model%chunk_begin(nb):Model%chunk_end(nb),22)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'vbdsf_ave'
    ExtDiag(idx)%desc = 'Visible Beam Downward Solar Flux'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%time_avg_kind = 'rad_sw'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%fluxr(Model%chunk_begin(nb):Model%chunk_end(nb),24)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'vddsf_ave'
    ExtDiag(idx)%desc = 'Visible Diffuse Downward Solar Flux'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%time_avg_kind = 'rad_sw'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%fluxr(Model%chunk_begin(nb):Model%chunk_end(nb),25)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'nbdsf_ave'
    ExtDiag(idx)%desc = 'Near IR Beam Downward Solar Flux'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%time_avg_kind = 'rad_sw'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%fluxr(Model%chunk_begin(nb):Model%chunk_end(nb),26)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'nddsf_ave'
    ExtDiag(idx)%desc = 'Near IR Diffuse Downward Solar Flux'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%time_avg_kind = 'rad_sw'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%fluxr(Model%chunk_begin(nb):Model%chunk_end(nb),27)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'csulf_avetoa'
    ExtDiag(idx)%desc = 'Clear Sky Upward Long Wave Flux at toa'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%time_avg_kind = 'rad_lw'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%fluxr(Model%chunk_begin(nb):Model%chunk_end(nb),28)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'csusf_avetoa'
    ExtDiag(idx)%desc = 'Clear Sky Upward Short Wave Flux at toa'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%time_avg_kind = 'rad_sw'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%fluxr(Model%chunk_begin(nb):Model%chunk_end(nb),29)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'csdlf_ave'
    ExtDiag(idx)%desc = 'Clear Sky Downward Long Wave Flux'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%time_avg_kind = 'rad_lw'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%fluxr(Model%chunk_begin(nb):Model%chunk_end(nb),30)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'csusf_ave'
    ExtDiag(idx)%desc = 'Clear Sky Upward Short Wave Flux'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%time_avg_kind = 'rad_sw'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%fluxr(Model%chunk_begin(nb):Model%chunk_end(nb),31)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'csdsf_ave'
    ExtDiag(idx)%desc = 'Clear Sky Downward Short Wave Flux'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%time_avg_kind = 'rad_sw'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%fluxr(Model%chunk_begin(nb):Model%chunk_end(nb),32)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'csdsf'
    ExtDiag(idx)%desc = 'Clear Sky Instantateous Downward Short Wave Flux'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%fluxr(Model%chunk_begin(nb):Model%chunk_end(nb),32)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'csulf_ave'
    ExtDiag(idx)%desc = 'Clear Sky Upward Long Wave Flux'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%time_avg_kind = 'rad_lw'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%fluxr(Model%chunk_begin(nb):Model%chunk_end(nb),33)
    enddo


    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'DSWRFtoa'
    ExtDiag(idx)%desc = 'top of atmos downward shortwave flux'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%cnvfac = cn_one
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%time_avg_kind = 'rad_sw'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%fluxr(Model%chunk_begin(nb):Model%chunk_end(nb),23)
    enddo


    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'USWRFtoa'
    ExtDiag(idx)%desc = 'top of atmos upward shortwave flux'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%cnvfac = cn_one
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%time_avg_kind = 'rad_sw'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%fluxr(Model%chunk_begin(nb):Model%chunk_end(nb),2)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'ULWRFtoa'
    ExtDiag(idx)%desc = 'top of atmos upward longwave flux'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%cnvfac = cn_one
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%time_avg_kind = 'rad_lw'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%fluxr(Model%chunk_begin(nb):Model%chunk_end(nb),1)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'TCDC_aveclm'
    ExtDiag(idx)%desc = 'atmos column total cloud cover'
    ExtDiag(idx)%unit = '%'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%cnvfac = cn_100
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%time_avg_kind = 'rad_swlw_min'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%fluxr(Model%chunk_begin(nb):Model%chunk_end(nb),17)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'TCDC_avebndcl'
    ExtDiag(idx)%desc = 'boundary layer cloud layer total cloud cover'
    ExtDiag(idx)%unit = '%'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%cnvfac = cn_100
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%time_avg_kind = 'rad_swlw_min'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%fluxr(Model%chunk_begin(nb):Model%chunk_end(nb),18)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'TCDCcnvcl'
    ExtDiag(idx)%desc = 'convective cloud layer total cloud cover'
    ExtDiag(idx)%unit = '%'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%cnvfac = cn_100
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => Cldprop%cv(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'PREScnvclt'
    ExtDiag(idx)%desc = 'pressure at convective cloud top level'
    ExtDiag(idx)%unit = 'pa'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%mask = 'cldmask'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => Cldprop%cvt(Model%chunk_begin(nb):Model%chunk_end(nb))
      ExtDiag(idx)%data(nb)%var21 => Cldprop%cv(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'PREScnvclb'
    ExtDiag(idx)%desc = 'pressure at convective cloud bottom level'
    ExtDiag(idx)%unit = 'pa'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%mask = 'cldmask'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => Cldprop%cvb(Model%chunk_begin(nb):Model%chunk_end(nb))
      ExtDiag(idx)%data(nb)%var21 => Cldprop%cv(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo
!    if(mpp_pe()==mpp_root_pe())print *,'in gfdl_diag_register,af PREScnvclb,idx=',idx


    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'TCDC_avehcl'
    ExtDiag(idx)%desc = 'high cloud level total cloud cover'
    ExtDiag(idx)%unit = '%'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%cnvfac = cn_100
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%time_avg_kind = 'rad_swlw_min'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%fluxr(Model%chunk_begin(nb):Model%chunk_end(nb),5)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'PRES_avehct'
    ExtDiag(idx)%desc = 'pressure high cloud top level'
    ExtDiag(idx)%unit = '%'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%time_avg_kind = 'rad_swlw_min'
    ExtDiag(idx)%mask = "cldmask_ratio"
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%fluxr(Model%chunk_begin(nb):Model%chunk_end(nb),8)
      ExtDiag(idx)%data(nb)%var21 => IntDiag%fluxr(Model%chunk_begin(nb):Model%chunk_end(nb),5)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'PRES_avehcb'
    ExtDiag(idx)%desc = 'pressure high cloud bottom level'
    ExtDiag(idx)%unit = '%'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%time_avg_kind = 'rad_swlw_min'
    ExtDiag(idx)%mask = "cldmask_ratio"
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%fluxr(Model%chunk_begin(nb):Model%chunk_end(nb),11)
      ExtDiag(idx)%data(nb)%var21 => IntDiag%fluxr(Model%chunk_begin(nb):Model%chunk_end(nb),5)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'TEMP_avehct'
    ExtDiag(idx)%desc = 'temperature high cloud top level'
    ExtDiag(idx)%unit = 'K'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%time_avg_kind = 'rad_swlw_min'
    ExtDiag(idx)%mask = "cldmask_ratio"
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%fluxr(Model%chunk_begin(nb):Model%chunk_end(nb),14)
      ExtDiag(idx)%data(nb)%var21 => IntDiag%fluxr(Model%chunk_begin(nb):Model%chunk_end(nb),5)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'TCDC_avemcl'
    ExtDiag(idx)%desc = 'mid cloud level total cloud cover'
    ExtDiag(idx)%unit = '%'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%cnvfac = cn_100
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%time_avg_kind = 'rad_swlw_min'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%fluxr(Model%chunk_begin(nb):Model%chunk_end(nb),6)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'PRES_avemct'
    ExtDiag(idx)%desc = 'pressure middle cloud top level'
    ExtDiag(idx)%unit = '%'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%time_avg_kind = 'rad_swlw_min'
    ExtDiag(idx)%mask = "cldmask_ratio"
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%fluxr(Model%chunk_begin(nb):Model%chunk_end(nb),9)
      ExtDiag(idx)%data(nb)%var21 => IntDiag%fluxr(Model%chunk_begin(nb):Model%chunk_end(nb),6)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'PRES_avemcb'
    ExtDiag(idx)%desc = 'pressure middle cloud bottom level'
    ExtDiag(idx)%unit = '%'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%time_avg_kind = 'rad_swlw_min'
    ExtDiag(idx)%mask = "cldmask_ratio"
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%fluxr(Model%chunk_begin(nb):Model%chunk_end(nb),12)
      ExtDiag(idx)%data(nb)%var21 => IntDiag%fluxr(Model%chunk_begin(nb):Model%chunk_end(nb),6)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'TEMP_avemct'
    ExtDiag(idx)%desc = 'temperature middle cloud top level'
    ExtDiag(idx)%unit = 'K'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%time_avg_kind = 'rad_swlw_min'
    ExtDiag(idx)%mask = "cldmask_ratio"
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%fluxr(Model%chunk_begin(nb):Model%chunk_end(nb),15)
      ExtDiag(idx)%data(nb)%var21 => IntDiag%fluxr(Model%chunk_begin(nb):Model%chunk_end(nb),6)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'TCDC_avelcl'
    ExtDiag(idx)%desc = 'low cloud level total cloud cover'
    ExtDiag(idx)%unit = '%'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%cnvfac = cn_100
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%time_avg_kind = 'rad_swlw_min'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%fluxr(Model%chunk_begin(nb):Model%chunk_end(nb),7)
    enddo

   idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'PRES_avelct'
    ExtDiag(idx)%desc = 'pressure low cloud top level'
    ExtDiag(idx)%unit = '%'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%time_avg_kind = 'rad_swlw_min'
    ExtDiag(idx)%mask = "cldmask_ratio"
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%fluxr(Model%chunk_begin(nb):Model%chunk_end(nb),10)
      ExtDiag(idx)%data(nb)%var21 => IntDiag%fluxr(Model%chunk_begin(nb):Model%chunk_end(nb),7)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'PRES_avelcb'
    ExtDiag(idx)%desc = 'pressure low cloud bottom level'
    ExtDiag(idx)%unit = '%'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%time_avg_kind = 'rad_swlw_min'
    ExtDiag(idx)%mask = "cldmask_ratio"
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%fluxr(Model%chunk_begin(nb):Model%chunk_end(nb),13)
      ExtDiag(idx)%data(nb)%var21 => IntDiag%fluxr(Model%chunk_begin(nb):Model%chunk_end(nb),7)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'TEMP_avelct'
    ExtDiag(idx)%desc = 'temperature low cloud top level'
    ExtDiag(idx)%unit = 'K'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%time_avg_kind = 'rad_swlw_min'
    ExtDiag(idx)%mask = "cldmask_ratio"
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%fluxr(Model%chunk_begin(nb):Model%chunk_end(nb),16)
      ExtDiag(idx)%data(nb)%var21 => IntDiag%fluxr(Model%chunk_begin(nb):Model%chunk_end(nb),7)
    enddo
!--- aerosol diagnostics ---
    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'AOD_550'
    ExtDiag(idx)%desc = 'total aerosol optical depth at 550 nm'
    ExtDiag(idx)%unit = 'numerical'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%fluxr(Model%chunk_begin(nb):Model%chunk_end(nb),34)
    enddo

!--- aerosol diagnostics ---
    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'DU_AOD_550'
    ExtDiag(idx)%desc = 'dust aerosol optical depth at 550 nm'
    ExtDiag(idx)%unit = 'numerical'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%fluxr(Model%chunk_begin(nb):Model%chunk_end(nb),35)
    enddo

!--- aerosol diagnostics ---
    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'BC_AOD_550'
    ExtDiag(idx)%desc = 'soot aerosol optical depth at 550 nm'
    ExtDiag(idx)%unit = 'numerical'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%fluxr(Model%chunk_begin(nb):Model%chunk_end(nb),36)
    enddo

!--- aerosol diagnostics ---
    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'OC_AOD_550'
    ExtDiag(idx)%desc = 'waso aerosol optical depth at 550 nm'
    ExtDiag(idx)%unit = 'numerical'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%fluxr(Model%chunk_begin(nb):Model%chunk_end(nb),37)
    enddo

!--- aerosol diagnostics ---
    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'SU_AOD_550'
    ExtDiag(idx)%desc = 'suso aerosol optical depth at 550 nm'
    ExtDiag(idx)%unit = 'numerical'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%fluxr(Model%chunk_begin(nb):Model%chunk_end(nb),38)
    enddo

!--- aerosol diagnostics ---
    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'SS_AOD_550'
    ExtDiag(idx)%desc = 'salt aerosol optical depth at 550 nm'
    ExtDiag(idx)%unit = 'numerical'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%fluxr(Model%chunk_begin(nb):Model%chunk_end(nb),39)
    enddo

!--- air quality diagnostics ---
  if (Model%cplaqm) then
    if (associated(IntDiag%aod)) then
      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'aod'
      ExtDiag(idx)%desc = 'total aerosol optical depth at 550 nm'
      ExtDiag(idx)%unit = 'numerical'
      ExtDiag(idx)%mod_name = 'gfs_phys'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => IntDiag%aod(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo
    endif
  endif

!
!
!--- accumulated diagnostics ---
    do num = 1,NFXR
      write (xtra,'(I2.2)') num
      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'fluxr_'//trim(xtra)
      ExtDiag(idx)%desc = 'fluxr diagnostic '//trim(xtra)//' - GFS radiation'
      ExtDiag(idx)%unit = 'XXX'
      ExtDiag(idx)%mod_name = 'gfs_phys'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => IntDiag%fluxr(Model%chunk_begin(nb):Model%chunk_end(nb),num)
      enddo
    enddo

!--- the next two appear to be appear to be coupling fields in gloopr
!--- each has four elements
!rab    do num = 1,4
!rab      write (xtra,'(I1)') num
!rab      idx = idx + 1
!rab      ExtDiag(idx)%axes = 2
!rab      ExtDiag(idx)%name = 'dswcmp_'//trim(xtra)
!rab      ExtDiag(idx)%desc = 'dswcmp dagnostic '//trim(xtra)//' - GFS radiation'
!rab      ExtDiag(idx)%unit = 'XXX'
!rab      ExtDiag(idx)%mod_name = 'gfs_phys'
!rab      do nb = 1,nblks
!rab        ExtDiag(idx)%data(nb)%var2 => IntDiag%dswcmp(Model%chunk_begin(nb):Model%chunk_end(nb),num)
!rab      enddo
!rab    enddo
!rab
!rab    do num = 1,4
!rab      write (xtra,'(I1)') num
!rab      idx = idx + 1
!rab      ExtDiag(idx)%axes = 2
!rab      ExtDiag(idx)%name = 'uswcmp_'//trim(xtra)
!rab      ExtDiag(idx)%desc = 'uswcmp dagnostic '//trim(xtra)//' - GFS radiation'
!rab      ExtDiag(idx)%unit = 'XXX'
!rab      ExtDiag(idx)%mod_name = 'gfs_phys'
!rab      do nb = 1,nblks
!rab        ExtDiag(idx)%data(nb)%var2 => IntDiag%uswcmp(Model%chunk_begin(nb):Model%chunk_end(nb),num)
!rab      enddo
!rab    enddo

! DH gfortran cannot point to members of arrays of derived types such
! as IntDiag(nb)%topfsw(:)%upfxc (the compilation succeeds, but the
! pointers do not reference the correct data and the output either
! contains garbage (Inf, NaN), or the netCDF I/O layer crashes.
#ifndef __GFORTRAN__
    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'sw_upfxc'
    ExtDiag(idx)%desc = 'total sky upward sw flux at toa - GFS radiation'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%topfsw(Model%chunk_begin(nb):Model%chunk_end(nb))%upfxc
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'sw_dnfxc'
    ExtDiag(idx)%desc = 'total sky downward sw flux at toa - GFS radiation'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%topfsw(Model%chunk_begin(nb):Model%chunk_end(nb))%dnfxc
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'sw_upfx0'
    ExtDiag(idx)%desc = 'clear sky upward sw flux at toa - GFS radiation'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%topfsw(Model%chunk_begin(nb):Model%chunk_end(nb))%upfx0
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'lw_upfxc'
    ExtDiag(idx)%desc = 'total sky upward lw flux at toa - GFS radiation'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%topflw(Model%chunk_begin(nb):Model%chunk_end(nb))%upfxc
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'lw_upfx0'
    ExtDiag(idx)%desc = 'clear sky upward lw flux at toa - GFS radiation'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%topflw(Model%chunk_begin(nb):Model%chunk_end(nb))%upfx0
    enddo
#endif

!--- physics accumulated diagnostics ---
    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'ssrun_acc'
    ExtDiag(idx)%desc = 'Accumulated surface storm water runoff'
    ExtDiag(idx)%unit = 'kg/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%srunoff(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'evbs_ave'
    ExtDiag(idx)%desc = 'Direct Evaporation from Bare Soil'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%evbsa(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'evcw_ave'
    ExtDiag(idx)%desc = 'Canopy water evaporation'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%evcwa(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'snohf'
    ExtDiag(idx)%desc = 'Snow Phase Change Heat Flux'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%snohfa(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    if (Model%lsm == Model%lsm_noahmp) then
     idx = idx + 1
     ExtDiag(idx)%axes = 2
     ExtDiag(idx)%name = 'pah_ave'
     ExtDiag(idx)%desc = ' Total Precipitation Advected Heat'
     ExtDiag(idx)%unit = 'W/m**2'
     ExtDiag(idx)%mod_name = 'gfs_phys'
     ExtDiag(idx)%time_avg = .TRUE.
     allocate (ExtDiag(idx)%data(nblks))
     do nb = 1,nblks
       ExtDiag(idx)%data(nb)%var2 => IntDiag%paha(Model%chunk_begin(nb):Model%chunk_end(nb))
     enddo
    endif

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'trans_ave'
    ExtDiag(idx)%desc = 'transpiration'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%transa(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'sbsno_ave'
    ExtDiag(idx)%desc = 'Sublimation (evaporation from snow)'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%sbsnoa(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'snowc_ave'
    ExtDiag(idx)%desc = 'snow cover - GFS lsm'
    ExtDiag(idx)%unit = '%'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%cnvfac = cn_100
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%snowca(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'snowc'
    ExtDiag(idx)%desc = 'snow cover '
    ExtDiag(idx)%unit = 'fraction'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => Sfcprop%sncovr(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'soilm'
    ExtDiag(idx)%desc = 'total column soil moisture content'
    ExtDiag(idx)%unit = 'kg/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%mask = "land_only"
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2  => IntDiag%soilm(Model%chunk_begin(nb):Model%chunk_end(nb))
      ExtDiag(idx)%data(nb)%var21 => Sfcprop%slmsk(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'tmpmin2m'
    ExtDiag(idx)%desc = 'min temperature at 2m height'
    ExtDiag(idx)%unit = 'K'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%tmpmin(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'tmpmax2m'
    ExtDiag(idx)%desc = 'max temperature at 2m height'
    ExtDiag(idx)%unit = 'K'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%tmpmax(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'dusfc'
    ExtDiag(idx)%desc = 'surface zonal momentum flux'
    ExtDiag(idx)%unit = 'N/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%cnvfac = cn_one
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%dusfc(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'dvsfc'
    ExtDiag(idx)%desc = 'surface meridional momentum flux'
    ExtDiag(idx)%unit = 'N/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%cnvfac = cn_one
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%dvsfc(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'shtfl_ave'
    ExtDiag(idx)%desc = 'surface sensible heat flux'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%cnvfac = cn_one
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%dtsfc(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'lhtfl_ave'
    ExtDiag(idx)%desc = 'surface latent heat flux'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%cnvfac = cn_one
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%dqsfc(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'totprcp_ave'
    ExtDiag(idx)%desc = 'surface precipitation rate'
    ExtDiag(idx)%unit = 'kg/m**2/s'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%cnvfac = cn_th
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%time_avg_kind = 'full'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%totprcp(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'totprcpb_ave'
    ExtDiag(idx)%desc = 'bucket surface precipitation rate'
    ExtDiag(idx)%unit = 'kg/m**2/s'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%cnvfac = cn_th
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%totprcpb(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'gflux_ave'
    ExtDiag(idx)%desc = 'surface ground heat flux'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%cnvfac = cn_one
    ExtDiag(idx)%time_avg = .TRUE.
!    ExtDiag(idx)%mask = "land_ice_only"
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2  => IntDiag%gflux(Model%chunk_begin(nb):Model%chunk_end(nb))
!      ExtDiag(idx)%data(nb)%var21 => Sfcprop%slmsk(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'dlwsfc'
    ExtDiag(idx)%desc = 'time accumulated downward lw flux at surface'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%dlwsfc(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'ulwsfc'
    ExtDiag(idx)%desc = 'time accumulated upward lw flux at surface'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%ulwsfc(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'sunsd_acc'
    ExtDiag(idx)%desc = 'Sunshine Duration'
    ExtDiag(idx)%unit = 's'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%suntim(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'watr_acc'
    ExtDiag(idx)%desc = 'total water runoff'
    ExtDiag(idx)%unit = 'kg/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%runoff(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'ecan_acc'
    ExtDiag(idx)%desc = 'total evaporation of intercepted water'
    ExtDiag(idx)%unit = 'kg/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%tecan(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'etran_acc'
    ExtDiag(idx)%desc = 'total plant transpiration'
    ExtDiag(idx)%unit = 'kg/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%tetran(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'edir_acc'
    ExtDiag(idx)%desc = 'total soil surface evaporation'
    ExtDiag(idx)%unit = 'kg/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%tedir(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    if (Model%lsm == Model%lsm_noahmp) then
     idx = idx + 1
     ExtDiag(idx)%axes = 2
     ExtDiag(idx)%name = 'wa_acc'
     ExtDiag(idx)%desc = 'total water storage in aquifer'
     ExtDiag(idx)%unit = 'kg/m**2'
     ExtDiag(idx)%mod_name = 'gfs_phys'
     allocate (ExtDiag(idx)%data(nblks))
     do nb = 1,nblks
       ExtDiag(idx)%data(nb)%var2 => IntDiag%twa(Model%chunk_begin(nb):Model%chunk_end(nb))
     enddo
    endif

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'pevpr_ave'
    ExtDiag(idx)%desc = 'averaged potential evaporation rate'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%ep(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'cwork_ave'
    ExtDiag(idx)%desc = 'cloud work function (valid only with sas)'
    ExtDiag(idx)%unit = 'J/kg'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%cldwrk(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo


    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'u-gwd_ave'
    ExtDiag(idx)%desc = 'surface zonal gravity wave stress'
    ExtDiag(idx)%unit = 'N/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%cnvfac = cn_one
    ExtDiag(idx)%time_avg = .TRUE.
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%dugwd(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'v-gwd_ave'
    ExtDiag(idx)%desc = 'surface meridional gravity wave stress'
    ExtDiag(idx)%unit = 'N/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%cnvfac = cn_one
    ExtDiag(idx)%time_avg = .TRUE.
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%dvgwd(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'psmean'
    ExtDiag(idx)%desc = 'surface pressure'
    ExtDiag(idx)%unit = 'kPa'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%psmean(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'cnvprcp_ave'
    ExtDiag(idx)%desc = 'averaged surface convective precipitation rate'
    ExtDiag(idx)%unit = 'kg/m**2/s'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%cnvfac = cn_th
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%time_avg_kind = 'full'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%cnvprcp(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'cnvprcpb_ave'
    ExtDiag(idx)%desc = 'averaged bucket surface convective precipitation rate'
    ExtDiag(idx)%unit = 'kg/m**2/s'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%cnvfac = cn_th
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%cnvprcpb(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'cnvprcp'
    ExtDiag(idx)%desc = 'surface convective precipitation rate'
    ExtDiag(idx)%unit = 'kg/m**2/s'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%cnvprcp(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'spfhmin2m'
    ExtDiag(idx)%desc = 'minimum specific humidity'
    ExtDiag(idx)%unit = 'kg/kg'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%spfhmin(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'spfhmax2m'
    ExtDiag(idx)%desc = 'maximum specific humidity'
    ExtDiag(idx)%unit = 'kg/kg'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%spfhmax(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'u10mmax'
    ExtDiag(idx)%desc = 'maximum (magnitude) u-wind'
    ExtDiag(idx)%unit = 'm/s'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'vector_bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%u10mmax(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'v10mmax'
    ExtDiag(idx)%desc = 'maximum (magnitude) v-wind'
    ExtDiag(idx)%unit = 'm/s'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'vector_bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%v10mmax(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'wind10mmax'
    ExtDiag(idx)%desc = 'maximum wind speed'
    ExtDiag(idx)%unit = 'm/s'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%wind10mmax(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo
    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'u10max'
    ExtDiag(idx)%desc = 'hourly maximum (magnitude) u-wind'
    ExtDiag(idx)%unit = 'm/s'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'vector_bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%u10max(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'v10max'
    ExtDiag(idx)%desc = 'hourly maximum (magnitude) v-wind'
    ExtDiag(idx)%unit = 'm/s'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'vector_bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%v10max(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'spd10max'
    ExtDiag(idx)%desc = 'hourly maximum wind speed'
    ExtDiag(idx)%unit = 'm/s'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%spd10max(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 't02max'
    ExtDiag(idx)%desc = 'max hourly 2m Temperature'
    ExtDiag(idx)%unit = 'K'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%t02max(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 't02min'
    ExtDiag(idx)%desc = 'min hourly 2m Temperature'
    ExtDiag(idx)%unit = 'K'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%t02min(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo
    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'rh02max'
    ExtDiag(idx)%desc = 'max hourly 2m RH'
    ExtDiag(idx)%unit = '%'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%rh02max(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'rh02min'
    ExtDiag(idx)%desc = 'min hourly 2m RH'
    ExtDiag(idx)%unit = '%'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%rh02min(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'pratemax'
    ExtDiag(idx)%desc = 'max hourly precipitation rate'
    ExtDiag(idx)%unit = 'mm h-1'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%pratemax(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'frzr'
    ExtDiag(idx)%desc = 'accumulated surface freezing rain'
    ExtDiag(idx)%unit = 'kg/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%cnvfac = cn_th
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%frzr(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'frzrb'
    ExtDiag(idx)%desc = 'accumulated surface freezing rain in bucket'
    ExtDiag(idx)%unit = 'kg/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%cnvfac = cn_th
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%frzrb(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'frozr'
    ExtDiag(idx)%desc = 'accumulated surface graupel'
    ExtDiag(idx)%unit = 'kg/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%cnvfac = cn_th
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%frozr(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'frozrb'
    ExtDiag(idx)%desc = 'accumulated surface graupel in bucket'
    ExtDiag(idx)%unit = 'kg/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%cnvfac = cn_th
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%frozrb(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'tsnowp'
    ExtDiag(idx)%desc = 'accumulated surface snow'
    ExtDiag(idx)%unit = 'kg/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%cnvfac = cn_th
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%tsnowp(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'tsnowpb'
    ExtDiag(idx)%desc = 'accumulated surface snow in bucket'
    ExtDiag(idx)%unit = 'kg/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%cnvfac = cn_th
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%tsnowpb(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'rhonewsn'
    ExtDiag(idx)%desc = 'precipitation ice density'
    ExtDiag(idx)%unit = 'kg m^-3'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%rhonewsn1(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'rain'
    ExtDiag(idx)%desc = 'total rain at this time step'
    ExtDiag(idx)%unit = 'XXX'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%rain(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'rainc'
    ExtDiag(idx)%desc = 'convective rain at this time step'
    ExtDiag(idx)%unit = 'XXX'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%rainc(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'ice'
    ExtDiag(idx)%desc = 'ice fall at this time step'
    ExtDiag(idx)%unit = 'XXX'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%ice(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'snow'
    ExtDiag(idx)%desc = 'snow fall at this time step'
    ExtDiag(idx)%unit = 'XXX'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%snow(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'graupel'
    ExtDiag(idx)%desc = 'graupel fall at this time step'
    ExtDiag(idx)%unit = 'XXX'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%graupel(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'totice_ave'
    ExtDiag(idx)%desc = 'surface ice precipitation rate'
    ExtDiag(idx)%unit = 'kg/m**2/s'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%cnvfac = cn_th
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%time_avg_kind = 'full'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%totice(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'toticeb_ave'
    ExtDiag(idx)%desc = 'bucket surface ice precipitation rate'
    ExtDiag(idx)%unit = 'kg/m**2/s'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%cnvfac = cn_th
    ExtDiag(idx)%time_avg = .TRUE.
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%toticeb(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'totsnw_ave'
    ExtDiag(idx)%desc = 'surface snow precipitation rate'
    ExtDiag(idx)%unit = 'kg/m**2/s'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%cnvfac = cn_th
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%time_avg_kind = 'full'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%totsnw(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'totsnwb_ave'
    ExtDiag(idx)%desc = 'bucket surface snow precipitation rate'
    ExtDiag(idx)%unit = 'kg/m**2/s'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%cnvfac = cn_th
    ExtDiag(idx)%time_avg = .TRUE.
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%totsnwb(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'totgrp_ave'
    ExtDiag(idx)%desc = 'surface graupel precipitation rate'
    ExtDiag(idx)%unit = 'kg/m**2/s'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%cnvfac = cn_th
    ExtDiag(idx)%time_avg = .TRUE.
    ExtDiag(idx)%time_avg_kind = 'full'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%totgrp(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'totgrpb_ave'
    ExtDiag(idx)%desc = 'bucket surface graupel precipitation rate'
    ExtDiag(idx)%unit = 'kg/m**2/s'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%cnvfac = cn_th
    ExtDiag(idx)%time_avg = .TRUE.
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%totgrpb(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

!    if(mpp_pe()==mpp_root_pe())print *,'in gfdl_diag_register,af totgrp,idx=',idx

    if(associated(Coupling%sfcdlw)) then
    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'sfcdlw'
    ExtDiag(idx)%desc = 'sfcdlw'
    ExtDiag(idx)%unit = 'W m-2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => Coupling%sfcdlw(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo
    endif

    if(associated(Coupling%htrlw)) then
    idx = idx + 1
    ExtDiag(idx)%axes = 3
    ExtDiag(idx)%name = 'htrlw'
    ExtDiag(idx)%desc = 'htrlw'
    ExtDiag(idx)%unit = 'W m-2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var3 => Coupling%htrlw(Model%chunk_begin(nb):Model%chunk_end(nb),:)
    enddo
    endif

    if(associated(Radtend%lwhc)) then
    idx = idx + 1
    ExtDiag(idx)%axes = 3
    ExtDiag(idx)%name = 'lwhc'
    ExtDiag(idx)%desc = 'lwhc'
    ExtDiag(idx)%unit = 'K s-1'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var3 => Radtend%lwhc(Model%chunk_begin(nb):Model%chunk_end(nb),:)
    enddo
    endif

!--- physics instantaneous diagnostics ---
    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'u10m'
    ExtDiag(idx)%desc = '10 meter u wind'
    ExtDiag(idx)%unit = 'm/s'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'vector_bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%u10m(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'v10m'
    ExtDiag(idx)%desc = '10 meter v wind'
    ExtDiag(idx)%unit = 'm/s'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'vector_bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%v10m(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'dpt2m'
    ExtDiag(idx)%desc = '2 meter dew point temperature'
    ExtDiag(idx)%unit = 'K'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%dpt2m(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'hgt_hyblev1'
    ExtDiag(idx)%desc = 'layer 1 height'
    ExtDiag(idx)%unit = 'm'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%zlvl(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'psurf'
    ExtDiag(idx)%desc = 'surface pressure'
    ExtDiag(idx)%unit = 'Pa'
    ExtDiag(idx)%mask = 'pseudo_ps'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%psurf(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'hpbl'
    ExtDiag(idx)%desc = 'surface planetary boundary layer height'
    ExtDiag(idx)%unit = 'm'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => Tbd%hpbl(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'pwat'
    ExtDiag(idx)%desc = 'atmos column precipitable water'
    ExtDiag(idx)%unit = 'kg/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%pwat(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'tmp_hyblev1'
    ExtDiag(idx)%desc = 'layer 1 temperature'
    ExtDiag(idx)%unit = 'K'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%t1(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'spfh_hyblev1'
    ExtDiag(idx)%desc = 'layer 1 specific humidity'
    ExtDiag(idx)%unit = 'kg/kg'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%q1(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'ugrd_hyblev1'
    ExtDiag(idx)%desc = 'layer 1 zonal wind'
    ExtDiag(idx)%unit = 'm/s'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'vector_bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%u1(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'vgrd_hyblev1'
    ExtDiag(idx)%desc = 'layer 1 meridional wind'
    ExtDiag(idx)%unit = 'm/s'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'vector_bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%v1(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'sfexc'
    ExtDiag(idx)%desc = 'Exchange Coefficient'
    ExtDiag(idx)%unit = 'kg/m2/s'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%chh(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'acond'
    ExtDiag(idx)%desc = 'Aerodynamic conductance'
    ExtDiag(idx)%unit = 'm/s'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%cmm(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'dlwsfci'
    ExtDiag(idx)%desc = 'instantaneous sfc downward lw flux'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%dlwsfci(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'ulwsfci'
    ExtDiag(idx)%desc = 'instantaneous sfc upward lw flux'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%ulwsfci(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'dswsfci'
    ExtDiag(idx)%desc = 'instantaneous sfc downward sw flux'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%dswsfci(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'uswsfci'
    ExtDiag(idx)%desc = 'instantaneous sfc upward sw flux'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%uswsfci(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'dusfci'
    ExtDiag(idx)%desc = 'instantaneous u component of surface stress'
    ExtDiag(idx)%unit = 'XXX'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%dusfci(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'dvsfci'
    ExtDiag(idx)%desc = 'instantaneous v component of surface stress'
    ExtDiag(idx)%unit = 'XXX'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%dvsfci(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'shtfl'
    ExtDiag(idx)%desc = 'instantaneous surface sensible heat net flux'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%dtsfci(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'lhtfl'
    ExtDiag(idx)%desc = 'instantaneous surface latent heat net flux'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%dqsfci(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'gfluxi'
    ExtDiag(idx)%desc = 'instantaneous surface ground heat flux'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%gfluxi(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'wilt'
    ExtDiag(idx)%desc = 'wiltimg point (volumetric)'
    ExtDiag(idx)%unit = 'Proportion'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%smcwlt2(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'fldcp'
    ExtDiag(idx)%desc = 'Field Capacity (volumetric)'
    ExtDiag(idx)%unit = 'fraction'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%smcref2(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    if (Model%lsm == Model%lsm_noahmp) then
     idx = idx + 1
     ExtDiag(idx)%axes = 2
     ExtDiag(idx)%name = 'pahi'
     ExtDiag(idx)%desc = 'instantaneous precipitation advected heat flux'
     ExtDiag(idx)%unit = 'W/m**2'
     ExtDiag(idx)%mod_name = 'gfs_phys'
     allocate (ExtDiag(idx)%data(nblks))
     do nb = 1,nblks
       ExtDiag(idx)%data(nb)%var2 => IntDiag%pahi(Model%chunk_begin(nb):Model%chunk_end(nb))
     enddo
    endif

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'pevpr'
    ExtDiag(idx)%desc = 'instantaneous surface potential evaporation'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%epi(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'wet1'
    ExtDiag(idx)%desc = 'normalized soil wetness'
    ExtDiag(idx)%unit = 'XXX'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    if (Model%lsm==Model%lsm_ruc) then
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => Sfcprop%wetness(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo
    else
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => IntDiag%wet1(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo
    endif

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'cpofp'
    ExtDiag(idx)%desc = 'Percent frozen precipitation'
    ExtDiag(idx)%unit = 'fraction'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%sr(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'crain_ave'
    ExtDiag(idx)%desc = 'averaged categorical rain'
    ExtDiag(idx)%unit = 'number'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    ExtDiag(idx)%cnvfac = cn_one
    ExtDiag(idx)%time_avg = .TRUE.
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%tdomr(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'csnow_ave'
    ExtDiag(idx)%desc = 'averaged categorical snow'
    ExtDiag(idx)%unit = 'number'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    ExtDiag(idx)%cnvfac = cn_one
    ExtDiag(idx)%time_avg = .TRUE.
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%tdoms(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'cfrzr_ave'
    ExtDiag(idx)%desc = 'averaged categorical freezing rain'
    ExtDiag(idx)%unit = 'number'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    ExtDiag(idx)%cnvfac = cn_one
    ExtDiag(idx)%time_avg = .TRUE.
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%tdomzr(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'cicep_ave'
    ExtDiag(idx)%desc = 'averaged categorical sleet'
    ExtDiag(idx)%unit = 'number'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    ExtDiag(idx)%cnvfac = cn_one
    ExtDiag(idx)%time_avg = .TRUE.
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%tdomip(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 3
    ExtDiag(idx)%name = 'refl_10cm'
    ExtDiag(idx)%desc = 'Radar reflectivity'
    ExtDiag(idx)%unit = 'dBz'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var3 => IntDiag%refl_10cm(Model%chunk_begin(nb):Model%chunk_end(nb),:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'max_hail_diam_sfc'
    ExtDiag(idx)%desc = 'Maximum hail diameter at lowest model level'
    ExtDiag(idx)%unit = 'm'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%max_hail_diam_sfc(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 3
    ExtDiag(idx)%name = 'dkt'
    ExtDiag(idx)%desc = 'Atmospheric heat diffusivity'
    ExtDiag(idx)%unit = 'm2s-1'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var3 => IntDiag%dkt(Model%chunk_begin(nb):Model%chunk_end(nb),:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 3
    ExtDiag(idx)%name = 'dku'
    ExtDiag(idx)%desc = 'Atmospheric momentum diffusivity'
    ExtDiag(idx)%unit = 'm2s-1'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var3 => IntDiag%dku(Model%chunk_begin(nb):Model%chunk_end(nb),:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 3
    ExtDiag(idx)%name = 'cldfra'
    ExtDiag(idx)%desc = 'Instantaneous 3D Cloud Fraction'
    ExtDiag(idx)%unit = 'frac'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var3 => IntDiag%cldfra(Model%chunk_begin(nb):Model%chunk_end(nb),:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 3
    ExtDiag(idx)%name = 'cnvw'
    ExtDiag(idx)%desc = 'subgrid scale convective cloud water'
    ExtDiag(idx)%unit = 'kg/kg'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    if( Model%ncnvw > 0 ) then
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var3 => Tbd%phy_f3d(Model%chunk_begin(nb):Model%chunk_end(nb),:,Model%ncnvw)
      enddo
    endif

    if (Model%do_skeb) then
      idx = idx + 1
      ExtDiag(idx)%axes = 3
      ExtDiag(idx)%name = 'skebu_wts'
      ExtDiag(idx)%desc = 'perturbation velocity'
      ExtDiag(idx)%unit = 'm/s'
      ExtDiag(idx)%mod_name = 'gfs_phys'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var3 => Coupling%skebu_wts(Model%chunk_begin(nb):Model%chunk_end(nb),:)
      enddo

      idx = idx + 1
      ExtDiag(idx)%axes = 3
      ExtDiag(idx)%name = 'skebv_wts'
      ExtDiag(idx)%desc = 'perturbation velocity'
      ExtDiag(idx)%unit = 'm/s'
      ExtDiag(idx)%mod_name = 'gfs_phys'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var3 => Coupling%skebv_wts(Model%chunk_begin(nb):Model%chunk_end(nb),:)
      enddo
    endif

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'zmtnblck'
    ExtDiag(idx)%desc = 'level of dividing streamline'
    ExtDiag(idx)%unit = 'm/s'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%zmtnblck(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'refdmax'
    ExtDiag(idx)%desc = 'max hourly 1-km agl reflectivity'
    ExtDiag(idx)%unit = 'dBZ'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%refdmax(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo
    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'refdmax263k'
    ExtDiag(idx)%desc = 'max hourly -10C reflectivity'
    ExtDiag(idx)%unit = 'dBZ'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%refdmax263k(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    if (Model%do_sppt .or. Model%ca_global) then
      idx = idx + 1
      ExtDiag(idx)%axes = 3
      ExtDiag(idx)%name = 'sppt_wts'
      ExtDiag(idx)%desc = 'perturbation velocity'
      ExtDiag(idx)%unit = 'm/s'
      ExtDiag(idx)%mod_name = 'gfs_phys'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var3 => Coupling%sppt_wts(Model%chunk_begin(nb):Model%chunk_end(nb),:)
      enddo
    endif

    if (Model%do_shum) then
      idx = idx + 1
      ExtDiag(idx)%axes = 3
      ExtDiag(idx)%name = 'shum_wts'
      ExtDiag(idx)%desc = 'perturbation velocity'
      ExtDiag(idx)%unit = 'm/s'
      ExtDiag(idx)%mod_name = 'gfs_phys'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var3 => Coupling%shum_wts(Model%chunk_begin(nb):Model%chunk_end(nb),:)
      enddo
    endif

    if (Model%do_spp) then
      idx = idx + 1
      ExtDiag(idx)%axes = 3
      ExtDiag(idx)%name = 'spp_wts_pbl'
      ExtDiag(idx)%desc = 'spp pbl perturbation wts'
      ExtDiag(idx)%unit = 'm/s'
      ExtDiag(idx)%mod_name = 'gfs_phys'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var3 => Coupling%spp_wts_pbl(Model%chunk_begin(nb):Model%chunk_end(nb),:)
      enddo
    endif

    if (Model%do_spp) then
      idx = idx + 1
      ExtDiag(idx)%axes = 3
      ExtDiag(idx)%name = 'spp_wts_sfc'
      ExtDiag(idx)%desc = 'spp sfc perturbation wts'
      ExtDiag(idx)%unit = 'm/s'
      ExtDiag(idx)%mod_name = 'gfs_phys'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var3 => Coupling%spp_wts_sfc(Model%chunk_begin(nb):Model%chunk_end(nb),:)
      enddo
    endif

    if (Model%do_spp) then
      idx = idx + 1
      ExtDiag(idx)%axes = 3
      ExtDiag(idx)%name = 'spp_wts_mp'
      ExtDiag(idx)%desc = 'spp mp perturbation wts'
      ExtDiag(idx)%unit = 'm/s'
      ExtDiag(idx)%mod_name = 'gfs_phys'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var3 => Coupling%spp_wts_mp(Model%chunk_begin(nb):Model%chunk_end(nb),:)
      enddo
    endif

    if (Model%do_spp) then
      idx = idx + 1
      ExtDiag(idx)%axes = 3
      ExtDiag(idx)%name = 'spp_wts_gwd'
      ExtDiag(idx)%desc = 'spp gwd perturbation wts'
      ExtDiag(idx)%unit = 'm/s'
      ExtDiag(idx)%mod_name = 'gfs_phys'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var3 => Coupling%spp_wts_gwd(Model%chunk_begin(nb):Model%chunk_end(nb),:)
      enddo
    endif

    if (Model%do_spp) then
      idx = idx + 1
      ExtDiag(idx)%axes = 3
      ExtDiag(idx)%name = 'spp_wts_rad'
      ExtDiag(idx)%desc = 'spp rad perturbation wts'
      ExtDiag(idx)%unit = 'm/s'
      ExtDiag(idx)%mod_name = 'gfs_phys'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var3 => Coupling%spp_wts_rad(Model%chunk_begin(nb):Model%chunk_end(nb),:)
      enddo
    endif

    if (Model%do_spp) then
      idx = idx + 1
      ExtDiag(idx)%axes = 3
      ExtDiag(idx)%name = 'spp_wts_cu_deep'
      ExtDiag(idx)%desc = 'spp cu deep perturbation wts'
      ExtDiag(idx)%unit = 'm/s'
      ExtDiag(idx)%mod_name = 'gfs_phys'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var3 => Coupling%spp_wts_cu_deep(Model%chunk_begin(nb):Model%chunk_end(nb),:)
      enddo
    endif

    if (Model%lndp_type /= 0) then
      idx = idx + 1
      ExtDiag(idx)%axes = 3
      ExtDiag(idx)%name = 'sfc_wts'
      ExtDiag(idx)%desc = 'perturbation amplitude'
      ExtDiag(idx)%unit = 'none'
      ExtDiag(idx)%mod_name = 'gfs_phys'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var3 => Coupling%sfc_wts(Model%chunk_begin(nb):Model%chunk_end(nb),:)
      enddo

    endif

    if (Model%do_ca) then
      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'ca1'
      ExtDiag(idx)%desc = 'Cellular Automata'
      ExtDiag(idx)%unit = '%'
      ExtDiag(idx)%mod_name = 'gfs_phys'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => Coupling%ca1(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo

      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'ca_deep'
      ExtDiag(idx)%desc = 'CA deep conv'
      ExtDiag(idx)%unit = '%'
      ExtDiag(idx)%mod_name = 'gfs_phys'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => Coupling%ca_deep(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo

      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'ca_turb'
      ExtDiag(idx)%desc = 'CA turbulence'
      ExtDiag(idx)%unit = '%'
      ExtDiag(idx)%mod_name = 'gfs_phys'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => Coupling%ca_turb(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo

      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'ca_shal'
      ExtDiag(idx)%desc = 'CA shallow conv'
      ExtDiag(idx)%unit = '%'
      ExtDiag(idx)%mod_name = 'gfs_phys'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => Coupling%ca_shal(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo

      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'ca_rad'
      ExtDiag(idx)%desc = 'CA radiation'
      ExtDiag(idx)%unit = '%'
      ExtDiag(idx)%mod_name = 'gfs_phys'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => Coupling%ca_rad(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo

      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'ca_micro'
      ExtDiag(idx)%desc = 'CA microphys'
      ExtDiag(idx)%unit = '%'
      ExtDiag(idx)%mod_name = 'gfs_phys'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => Coupling%ca_micro(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo
    endif

  if (Model%lkm/=0) then

      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'lakefrac'
      ExtDiag(idx)%desc = 'Lake Fraction'
      ExtDiag(idx)%unit = 'fraction'
      ExtDiag(idx)%mod_name = 'gfs_sfc'
      ExtDiag(idx)%intpl_method = 'nearest_stod'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => Sfcprop%lakefrac(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo

      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'lakedepth'
      ExtDiag(idx)%desc = 'Lake Depth'
      ExtDiag(idx)%unit = 'm'
      ExtDiag(idx)%mod_name = 'gfs_sfc'
      ExtDiag(idx)%intpl_method = 'nearest_stod'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => Sfcprop%lakedepth(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo

      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'T_snow'
      ExtDiag(idx)%desc = 'Temperature of snow on a lake'
      ExtDiag(idx)%unit = 'K'
      ExtDiag(idx)%mod_name = 'gfs_sfc'
      ExtDiag(idx)%intpl_method = 'nearest_stod'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => Sfcprop%T_snow(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo

      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'T_ice'
      ExtDiag(idx)%desc = 'Temperature of ice on a lake'
      ExtDiag(idx)%unit = 'K'
      ExtDiag(idx)%mod_name = 'gfs_sfc'
      ExtDiag(idx)%intpl_method = 'nearest_stod'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => Sfcprop%T_ice(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo

      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'use_lake_model'
      ExtDiag(idx)%desc = 'Lake Model Flag'
      ExtDiag(idx)%unit = 'flag'
      ExtDiag(idx)%mod_name = 'gfs_sfc'
      ExtDiag(idx)%intpl_method = 'nearest_stod'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%int2 => Sfcprop%use_lake_model(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo

      if(Model%iopt_lake==Model%iopt_lake_clm) then

        ! Populate the 3D arrays separately since the code is complicated:
        call clm_lake_externaldiag_populate(ExtDiag, Model, Sfcprop, idx, cn_one, nblks)

        idx = idx + 1
        ExtDiag(idx)%axes = 2
        ExtDiag(idx)%name = 'lake_is_salty'
        ExtDiag(idx)%desc = 'lake point is considered salty by clm lake model'
        ExtDiag(idx)%unit = '1'
        ExtDiag(idx)%mod_name = 'gfs_sfc'
        allocate (ExtDiag(idx)%data(nblks))
        do nb = 1,nblks
           ExtDiag(idx)%data(nb)%int2 => Sfcprop%lake_is_salty(Model%chunk_begin(nb):Model%chunk_end(nb))
        enddo

        idx = idx + 1
        ExtDiag(idx)%axes = 2
        ExtDiag(idx)%name = 'lake_cannot_freeze'
        ExtDiag(idx)%desc = 'clm lake model considers the point to be so salty it cannot freeze'
        ExtDiag(idx)%unit = '1'
        ExtDiag(idx)%mod_name = 'gfs_sfc'
        allocate (ExtDiag(idx)%data(nblks))
        do nb = 1,nblks
           ExtDiag(idx)%data(nb)%int2 => Sfcprop%lake_cannot_freeze(Model%chunk_begin(nb):Model%chunk_end(nb))
        enddo

        idx = idx + 1
        ExtDiag(idx)%axes = 2
        ExtDiag(idx)%name = 'lake_t2m'
        ExtDiag(idx)%desc = 'Temperature at 2 m from Lake Model'
        ExtDiag(idx)%unit = 'K'
        ExtDiag(idx)%intpl_method = 'nearest_stod'
        ExtDiag(idx)%mod_name = 'gfs_sfc'
        allocate (ExtDiag(idx)%data(nblks))
        do nb = 1,nblks
          ExtDiag(idx)%data(nb)%var2 => Sfcprop%lake_t2m(Model%chunk_begin(nb):Model%chunk_end(nb))
        enddo

        idx = idx + 1
        ExtDiag(idx)%axes = 2
        ExtDiag(idx)%name = 'lake_q2m'
        ExtDiag(idx)%desc = '2m specific humidity from Lake Model'
        ExtDiag(idx)%unit = 'kg/kg'
        ExtDiag(idx)%mod_name = 'gfs_sfc'
        ExtDiag(idx)%intpl_method = 'nearest_stod'
        allocate (ExtDiag(idx)%data(nblks))
        do nb = 1,nblks
          ExtDiag(idx)%data(nb)%var2 => Sfcprop%lake_q2m(Model%chunk_begin(nb):Model%chunk_end(nb))
        enddo

        idx = idx + 1
        ExtDiag(idx)%axes = 2
        ExtDiag(idx)%name = 'lake_albedo'
        ExtDiag(idx)%desc = 'mid day surface albedo over lake'
        ExtDiag(idx)%unit = 'fraction'
        ExtDiag(idx)%mod_name = 'gfs_sfc'
        ExtDiag(idx)%intpl_method = 'nearest_stod'
        allocate (ExtDiag(idx)%data(nblks))
        do nb = 1,nblks
          ExtDiag(idx)%data(nb)%var2 => Sfcprop%lake_albedo(Model%chunk_begin(nb):Model%chunk_end(nb))
        enddo

        idx = idx + 1
        ExtDiag(idx)%axes = 2
        ExtDiag(idx)%name = 'lake_h2osno2d'
        ExtDiag(idx)%desc = 'water equiv of acc snow depth over lake'
        ExtDiag(idx)%unit = 'mm'
        ExtDiag(idx)%mod_name = 'gfs_sfc'
        ExtDiag(idx)%intpl_method = 'nearest_stod'
        allocate (ExtDiag(idx)%data(nblks))
        do nb = 1,nblks
          ExtDiag(idx)%data(nb)%var2 => Sfcprop%lake_h2osno2d(Model%chunk_begin(nb):Model%chunk_end(nb))
        enddo

        idx = idx + 1
        ExtDiag(idx)%axes = 2
        ExtDiag(idx)%name = 'lake_sndpth2d'
        ExtDiag(idx)%desc = 'actual acc snow depth over lake in clm lake model'
        ExtDiag(idx)%unit = 'm'
        ExtDiag(idx)%mod_name = 'gfs_sfc'
        ExtDiag(idx)%intpl_method = 'nearest_stod'
        allocate (ExtDiag(idx)%data(nblks))
        do nb = 1,nblks
          ExtDiag(idx)%data(nb)%var2 => Sfcprop%lake_sndpth2d(Model%chunk_begin(nb):Model%chunk_end(nb))
        enddo

        idx = idx + 1
        ExtDiag(idx)%axes = 2
        ExtDiag(idx)%name = 'lake_snl2d'
        ExtDiag(idx)%desc = 'snow layers in clm lake model (treated as integer)'
        ExtDiag(idx)%unit = 'count'
        ExtDiag(idx)%mod_name = 'gfs_sfc'
        ExtDiag(idx)%intpl_method = 'nearest_stod'
        allocate (ExtDiag(idx)%data(nblks))
        do nb = 1,nblks
          ExtDiag(idx)%data(nb)%var2 => Sfcprop%lake_snl2d(Model%chunk_begin(nb):Model%chunk_end(nb))
        enddo

        idx = idx + 1
        ExtDiag(idx)%axes = 2
        ExtDiag(idx)%name = 'lake_tsfc'
        ExtDiag(idx)%desc = 'skin temperature from clm lake model'
        ExtDiag(idx)%unit = 'K'
        ExtDiag(idx)%mod_name = 'gfs_sfc'
        ExtDiag(idx)%intpl_method = 'nearest_stod'
        allocate (ExtDiag(idx)%data(nblks))
        do nb = 1,nblks
          ExtDiag(idx)%data(nb)%var2 => Sfcprop%lake_tsfc(Model%chunk_begin(nb):Model%chunk_end(nb))
        enddo

        idx = idx + 1
        ExtDiag(idx)%axes = 2
        ExtDiag(idx)%name = 'lake_savedtke12d'
        ExtDiag(idx)%desc = 'top level eddy conductivity from previous timestep in clm lake model'
        ExtDiag(idx)%unit = 'kg m-3'
        ExtDiag(idx)%mod_name = 'gfs_sfc'
        ExtDiag(idx)%intpl_method = 'nearest_stod'
        allocate (ExtDiag(idx)%data(nblks))
        do nb = 1,nblks
          ExtDiag(idx)%data(nb)%var2 => Sfcprop%lake_savedtke12d(Model%chunk_begin(nb):Model%chunk_end(nb))
        enddo

        idx = idx + 1
        ExtDiag(idx)%axes = 2
        ExtDiag(idx)%name = 'lake_ht'
        ExtDiag(idx)%desc = 'lake_ht'
        ExtDiag(idx)%unit = 'unitless'
        ExtDiag(idx)%mod_name = 'gfs_sfc'
        ExtDiag(idx)%intpl_method = 'nearest_stod'
        allocate (ExtDiag(idx)%data(nblks))
        do nb = 1,nblks
          ExtDiag(idx)%data(nb)%var2 => Sfcprop%lake_ht(Model%chunk_begin(nb):Model%chunk_end(nb))
        enddo

      endif

  endif

  if (Model%ldiag_ugwp) THEN

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'zmtb'
    ExtDiag(idx)%desc = 'height of dividing streamline'
    ExtDiag(idx)%unit = 'm'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%zmtb(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'zogw'
    ExtDiag(idx)%desc = 'height of OGW-launch'
    ExtDiag(idx)%unit = 'm'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%zogw(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'zlwb'
    ExtDiag(idx)%desc = 'height of LWB-level'
    ExtDiag(idx)%unit = 'm'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%zlwb(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'tau_ogw'
    ExtDiag(idx)%desc = ' OGW vertical MF at launch level'
    ExtDiag(idx)%unit = 'N/m2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%tau_ogw(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'tau_mtb'
    ExtDiag(idx)%desc = ' ORO-MTB integrated flux from surface'
    ExtDiag(idx)%unit = 'N/m2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%tau_mtb(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'tau_tofd'
    ExtDiag(idx)%desc = ' ORO-TOFD integrated flux from surface'
    ExtDiag(idx)%unit = 'N/m2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%tau_tofd(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'tau_ngw'
    ExtDiag(idx)%desc = ' NGW momentum flux at launch level '
    ExtDiag(idx)%unit = 'N/m2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%tau_ngw(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 3
    ExtDiag(idx)%name = 'du3dt_ogw'
    ExtDiag(idx)%desc = 'axz_oro averaged E-W OROGW-tendency'
    ExtDiag(idx)%unit = 'm/s/s'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var3 => IntDiag%du3dt_ogw(Model%chunk_begin(nb):Model%chunk_end(nb),:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 3
    ExtDiag(idx)%name = 'du3dt_ngw'
    ExtDiag(idx)%desc = 'axz_oro averaged E-W GWALL-tendency'
    ExtDiag(idx)%unit = 'm/s/s'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var3 => IntDiag%du3dt_ngw(Model%chunk_begin(nb):Model%chunk_end(nb),:)
    enddo
!
!
    idx = idx + 1
    ExtDiag(idx)%axes = 3
    ExtDiag(idx)%name = 'du3dt_mtb'
    ExtDiag(idx)%desc = 'axz_oro averaged E-W MTB-tendency'
    ExtDiag(idx)%unit = 'm/s/s'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var3 => IntDiag%du3dt_mtb(Model%chunk_begin(nb):Model%chunk_end(nb),:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 3
    ExtDiag(idx)%name = 'du3dt_tms'
    ExtDiag(idx)%desc = 'axz_oro averaged E-W TMS-tendency'
    ExtDiag(idx)%unit = 'm/s/s'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var3 => IntDiag%du3dt_tms(Model%chunk_begin(nb):Model%chunk_end(nb),:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 3
    ExtDiag(idx)%name = 'dudt_ogw'
    ExtDiag(idx)%desc = 'x wind tendency from mesoscale OGWD'
    ExtDiag(idx)%unit = 'm s-2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var3 => IntDiag%dudt_ogw(Model%chunk_begin(nb):Model%chunk_end(nb),:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 3
    ExtDiag(idx)%name = 'dvdt_ogw'
    ExtDiag(idx)%desc = 'y wind tendency from mesoscale OGWD'
    ExtDiag(idx)%unit = 'm s-2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var3 => IntDiag%dvdt_ogw(Model%chunk_begin(nb):Model%chunk_end(nb),:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 3
    ExtDiag(idx)%name = 'dudt_obl'
    ExtDiag(idx)%desc = 'x wind tendency from blocking drag'
    ExtDiag(idx)%unit = 'm s-2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var3 => IntDiag%dudt_obl(Model%chunk_begin(nb):Model%chunk_end(nb),:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 3
    ExtDiag(idx)%name = 'dvdt_obl'
    ExtDiag(idx)%desc = 'y wind tendency from blocking drag'
    ExtDiag(idx)%unit = 'm s-2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var3 => IntDiag%dvdt_obl(Model%chunk_begin(nb):Model%chunk_end(nb),:)
    enddo

    ! 2D variables

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'du_ogwcol'
    ExtDiag(idx)%desc = 'integrated x momentum flux from meso scale ogw'
    ExtDiag(idx)%unit = 'Pa'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%du_ogwcol(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'dv_ogwcol'
    ExtDiag(idx)%desc = 'integrated y momentum flux from meso scale ogw'
    ExtDiag(idx)%unit = 'Pa'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%dv_ogwcol(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'du_oblcol'
    ExtDiag(idx)%desc = 'integrated x momentum flux from blocking drag'
    ExtDiag(idx)%unit = 'Pa'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%du_oblcol(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'dv_oblcol'
    ExtDiag(idx)%desc = 'integrated y momentum flux from blocking drag'
    ExtDiag(idx)%unit = 'Pa'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%dv_oblcol(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 3
    ExtDiag(idx)%name = 'dws3dt_ogw'
    ExtDiag(idx)%desc = 'averaged wind speed tendency due to mesoscale gravity wave drag'
    ExtDiag(idx)%unit = 'm s-2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var3 => IntDiag%dws3dt_ogw(Model%chunk_begin(nb):Model%chunk_end(nb),:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 3
    ExtDiag(idx)%name = 'dws3dt_obl'
    ExtDiag(idx)%desc = 'averaged wind speed tendency due to blocking drag'
    ExtDiag(idx)%unit = 'm s-2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var3 => IntDiag%dws3dt_obl(Model%chunk_begin(nb):Model%chunk_end(nb),:)
    enddo

    ! Variables for GSL drag suite

    idx = idx + 1
    ExtDiag(idx)%axes = 3
    ExtDiag(idx)%name = 'dudt_oss'
    ExtDiag(idx)%desc = 'x wind tendency from small scale GWD'
    ExtDiag(idx)%unit = 'm s-2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var3 => IntDiag%dudt_oss(Model%chunk_begin(nb):Model%chunk_end(nb),:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 3
    ExtDiag(idx)%name = 'dvdt_oss'
    ExtDiag(idx)%desc = 'y wind tendency from small scale GWD'
    ExtDiag(idx)%unit = 'm s-2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var3 => IntDiag%dvdt_oss(Model%chunk_begin(nb):Model%chunk_end(nb),:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 3
    ExtDiag(idx)%name = 'dudt_ofd'
    ExtDiag(idx)%desc = 'x wind tendency from form drag'
    ExtDiag(idx)%unit = 'm s-2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var3 => IntDiag%dudt_ofd(Model%chunk_begin(nb):Model%chunk_end(nb),:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 3
    ExtDiag(idx)%name = 'dvdt_ofd'
    ExtDiag(idx)%desc = 'y wind tendency from form drag'
    ExtDiag(idx)%unit = 'm s-2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var3 => IntDiag%dvdt_ofd(Model%chunk_begin(nb):Model%chunk_end(nb),:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 3
    ExtDiag(idx)%name = 'dws3dt_oss'
    ExtDiag(idx)%desc = 'averaged wind speed tendency due to small-scale gravity wave drag'
    ExtDiag(idx)%unit = 'm s-2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var3 => IntDiag%dws3dt_oss(Model%chunk_begin(nb):Model%chunk_end(nb),:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 3
    ExtDiag(idx)%name = 'dws3dt_ofd'
    ExtDiag(idx)%desc = 'averaged wind speed tendency due to turbulent orographic form drag'
    ExtDiag(idx)%unit = 'm s-2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var3 => IntDiag%dws3dt_ofd(Model%chunk_begin(nb):Model%chunk_end(nb),:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 3
    ExtDiag(idx)%name = 'ldu3dt_ogw'
    ExtDiag(idx)%desc = 'averaged x wind tendency due to mesoscale orographic gravity wave drag'
    ExtDiag(idx)%unit = 'm s-2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var3 => IntDiag%ldu3dt_ogw(Model%chunk_begin(nb):Model%chunk_end(nb),:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 3
    ExtDiag(idx)%name = 'ldu3dt_obl'
    ExtDiag(idx)%desc = 'averaged x wind tendency due to blocking drag'
    ExtDiag(idx)%unit = 'm s-2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var3 => IntDiag%ldu3dt_obl(Model%chunk_begin(nb):Model%chunk_end(nb),:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 3
    ExtDiag(idx)%name = 'ldu3dt_ofd'
    ExtDiag(idx)%desc = 'averaged x wind tendency due to form drag'
    ExtDiag(idx)%unit = 'm s-2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var3 => IntDiag%ldu3dt_ofd(Model%chunk_begin(nb):Model%chunk_end(nb),:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 3
    ExtDiag(idx)%name = 'ldu3dt_oss'
    ExtDiag(idx)%desc = 'averaged x wind tendency due to small scale gravity wave drag'
    ExtDiag(idx)%unit = 'm s-2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var3 => IntDiag%ldu3dt_oss(Model%chunk_begin(nb):Model%chunk_end(nb),:)
    enddo

    ! 2D variables

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'du_osscol'
    ExtDiag(idx)%desc = 'integrated x momentum flux from small scale gwd'
    ExtDiag(idx)%unit = 'Pa'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%du_osscol(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'dv_osscol'
    ExtDiag(idx)%desc = 'integrated y momentum flux from small scale gwd'
    ExtDiag(idx)%unit = 'Pa'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%dv_osscol(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'du_ofdcol'
    ExtDiag(idx)%desc = 'integrated x momentum flux from form drag'
    ExtDiag(idx)%unit = 'Pa'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%du_ofdcol(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'dv_ofdcol'
    ExtDiag(idx)%desc = 'integrated y momentum flux from form drag'
    ExtDiag(idx)%unit = 'Pa'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%dv_ofdcol(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'du3_ogwcol'
    ExtDiag(idx)%desc = 'time averaged surface x momentum flux from mesoscale orographic gravity wave drag'
    ExtDiag(idx)%unit = 'Pa'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%du3_ogwcol(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'dv3_ogwcol'
    ExtDiag(idx)%desc = 'time averaged surface y momentum flux from mesoscale orographic gravity wave drag'
    ExtDiag(idx)%unit = 'Pa'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%dv3_ogwcol(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'du3_oblcol'
    ExtDiag(idx)%desc = 'time averaged surface x momentum flux from blocking drag'
    ExtDiag(idx)%unit = 'Pa'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%du3_oblcol(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'dv3_oblcol'
    ExtDiag(idx)%desc = 'time averaged surface y momentum flux from blocking drag'
    ExtDiag(idx)%unit = 'Pa'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%dv3_oblcol(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'du3_osscol'
    ExtDiag(idx)%desc = 'time averaged surface x momentum flux from small scale gravity wave drag'
    ExtDiag(idx)%unit = 'Pa'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%du3_osscol(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'dv3_osscol'
    ExtDiag(idx)%desc = 'time averaged surface y momentum flux from small scale gravity wave drag'
    ExtDiag(idx)%unit = 'Pa'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%dv3_osscol(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'du3_ofdcol'
    ExtDiag(idx)%desc = 'time averaged surface x momentum flux from form drag'
    ExtDiag(idx)%unit = 'Pa'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%du3_ofdcol(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'dv3_ofdcol'
    ExtDiag(idx)%desc = 'time averaged surface y momentum flux from form drag'
    ExtDiag(idx)%unit = 'Pa'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%dv3_ofdcol(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    ! UGWP non-stationary GWD outputs

    idx = idx + 1
    ExtDiag(idx)%axes = 3
    ExtDiag(idx)%name = 'ldu3dt_ngw'
    ExtDiag(idx)%desc = 'time averaged u momentum tendency due to non-stationary gravity wave drag'
    ExtDiag(idx)%unit = 'm s-2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var3 => IntDiag%ldu3dt_ngw(Model%chunk_begin(nb):Model%chunk_end(nb),:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 3
    ExtDiag(idx)%name = 'ldv3dt_ngw'
    ExtDiag(idx)%desc = 'time averaged v momentum tendency due to non-stationary gravity wave drag'
    ExtDiag(idx)%unit = 'm s-2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var3 => IntDiag%ldv3dt_ngw(Model%chunk_begin(nb):Model%chunk_end(nb),:)
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 3
    ExtDiag(idx)%name = 'ldt3dt_ngw'
    ExtDiag(idx)%desc = 'time averaged temperature tendency due to non-stationary gravity wave drag'
    ExtDiag(idx)%unit = 'K s-1'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%time_avg = .TRUE.
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var3 => IntDiag%ldt3dt_ngw(Model%chunk_begin(nb):Model%chunk_end(nb),:)
    enddo

  ENDIF  ! if (Model%ldiag_ugwp)


!    if(mpp_pe()==mpp_root_pe())print *,'in gfdl_diag_register,af shum_wts,idx=',idx

!--- Three-dimensional diagnostic tendencies stored in a 4D sparse
!--- array need special handling:
    if_ldiag3d: if(Model%ldiag3d) then
      do iprocess=1,Model%nprocess
        do itrac=1,Model%ntracp100
          if(Model%dtidx(itrac,iprocess)>=1) then
            call add_dtend(Model,ExtDiag,IntDiag,idx,nblks,itrac,iprocess)
          endif
        enddo
      enddo

      if_qdiag3d: if(Model%qdiag3d) then

        idx = idx + 1
        ExtDiag(idx)%axes = 3
        ExtDiag(idx)%name = 'upd_mf'
        ExtDiag(idx)%desc = 'updraft convective mass flux'
        ExtDiag(idx)%unit = 'kg m-1 s-3'
        ExtDiag(idx)%mod_name = 'gfs_phys'
        ExtDiag(idx)%time_avg = .TRUE.
        allocate (ExtDiag(idx)%data(nblks))
        do nb = 1,nblks
          ExtDiag(idx)%data(nb)%var3 => IntDiag%upd_mf(Model%chunk_begin(nb):Model%chunk_end(nb),:)
        enddo

        idx = idx + 1
        ExtDiag(idx)%axes = 3
        ExtDiag(idx)%name = 'dwn_mf'
        ExtDiag(idx)%desc = 'downdraft convective mass flux'
        ExtDiag(idx)%unit = 'kg m-1 s-3'
        ExtDiag(idx)%mod_name = 'gfs_phys'
        ExtDiag(idx)%time_avg = .TRUE.
        allocate (ExtDiag(idx)%data(nblks))
        do nb = 1,nblks
          ExtDiag(idx)%data(nb)%var3 => IntDiag%dwn_mf(Model%chunk_begin(nb):Model%chunk_end(nb),:)
        enddo

        idx = idx + 1
        ExtDiag(idx)%axes = 3
        ExtDiag(idx)%name = 'det_mf'
        ExtDiag(idx)%desc = 'detrainment convective mass flux'
        ExtDiag(idx)%unit = 'kg m-1 s-3'
        ExtDiag(idx)%mod_name = 'gfs_phys'
        ExtDiag(idx)%time_avg = .TRUE.
        allocate (ExtDiag(idx)%data(nblks))
        do nb = 1,nblks
          ExtDiag(idx)%data(nb)%var3 => IntDiag%det_mf(Model%chunk_begin(nb):Model%chunk_end(nb),:)
        enddo

      end if if_qdiag3d

    end if if_ldiag3d

!rab
!rab    do num = 1,5+Mdl_parms%pl_coeff
!rab      write (xtra,'(I1)') num
!rab      idx = idx + 1
!rab      ExtDiag(idx)%axes = 3
!rab      ExtDiag(idx)%name = 'dtend_'//trim(xtra)
!rab      ExtDiag(idx)%desc = 'moisture change due to physics '//trim(xtra)//''
!rab      ExtDiag(idx)%unit = 'XXX'
!rab      ExtDiag(idx)%mod_name = 'gfs_phys'
!rab    enddo
!rab
!rab    idx = idx + 1
!rab    ExtDiag(idx)%axes = 3
!rab    !Requires lgocart = .T.
!rab    ExtDiag(idx)%name = 'dqdt_v'
!rab    ExtDiag(idx)%desc = 'instantaneous total moisture tendency'
!rab    ExtDiag(idx)%unit = 'XXX'
!rab    ExtDiag(idx)%mod_name = 'gfs_phys'

!--- Surface diagnostics in gfs_sfc
    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'alnsf'
    ExtDiag(idx)%desc = 'mean nir albedo with strong cosz dependency'
    ExtDiag(idx)%unit = '%'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => Sfcprop%alnsf(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'alnwf'
    ExtDiag(idx)%desc = 'mean nir albedo with weak cosz dependency'
    ExtDiag(idx)%unit = '%'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => Sfcprop%alnwf(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'alvsf'
    ExtDiag(idx)%desc = 'mean vis albedo with strong cosz dependency'
    ExtDiag(idx)%unit = '%'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => Sfcprop%alvsf(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'alvwf'
    ExtDiag(idx)%desc = 'mean vis albedo with weak cosz dependency'
    ExtDiag(idx)%unit = '%'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => Sfcprop%alvwf(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'canopy'
    ExtDiag(idx)%desc = 'canopy water (cnwat in gfs data)'
    ExtDiag(idx)%unit = 'mm'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => Sfcprop%canopy(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'f10m'
    ExtDiag(idx)%desc = '10-meter wind speed divided by lowest model wind speed'
    ExtDiag(idx)%unit = 'N/A'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => Sfcprop%f10m(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'facsf'
    ExtDiag(idx)%desc = 'fractional coverage with strong cosz dependency'
    ExtDiag(idx)%unit = 'XXX'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => Sfcprop%facsf(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'facwf'
    ExtDiag(idx)%desc = 'fractional coverage with weak cosz dependency'
    ExtDiag(idx)%unit = 'XXX'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => Sfcprop%facwf(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'ffhh'
    ExtDiag(idx)%desc = 'fh parameter from PBL scheme'
    ExtDiag(idx)%unit = 'XXX'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => Sfcprop%ffhh(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'ffmm'
    ExtDiag(idx)%desc = 'fm parameter from PBL scheme'
    ExtDiag(idx)%unit = 'XXX'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => Sfcprop%ffmm(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'uustar'
    ExtDiag(idx)%desc = 'uustar surface frictional wind'
    ExtDiag(idx)%unit = 'XXX'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => Sfcprop%uustar(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'slope'
    ExtDiag(idx)%desc = 'surface slope type'
    ExtDiag(idx)%unit = 'XXX'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%int2 => Sfcprop%slope(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'fice'
    ExtDiag(idx)%desc = 'surface ice concentration (ice=1; no ice=0)'
    ExtDiag(idx)%unit = 'fraction'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => Sfcprop%fice(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'hice'
    ExtDiag(idx)%desc = 'sea ice thickness (icetk in gfs_data)'
    ExtDiag(idx)%unit = 'm'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => Sfcprop%hice(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'snoalb'
    ExtDiag(idx)%desc = 'maximum snow albedo in fraction'
    ExtDiag(idx)%unit = 'XXX'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => Sfcprop%snoalb(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'shdmax'
    ExtDiag(idx)%desc = 'maximum fractional coverage of green vegetation'
    ExtDiag(idx)%unit = 'XXX'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => Sfcprop%shdmax(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'shdmin'
    ExtDiag(idx)%desc = 'minimum fractional coverage of green vegetation'
    ExtDiag(idx)%unit = 'XXX'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => Sfcprop%shdmin(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'snowd'
    ExtDiag(idx)%desc = 'surface snow depth'
    ExtDiag(idx)%unit = 'm'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    ExtDiag(idx)%cnvfac = cn_one/cn_th
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => Sfcprop%snowd(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'sbsno'
    ExtDiag(idx)%desc = 'instantaneous sublimation (evaporation from snow)'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%sbsno(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'evbs'
    ExtDiag(idx)%desc = 'instantaneous direct evaporation over land'
    ExtDiag(idx)%unit = 'W m-2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%evbs(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'evcw'
    ExtDiag(idx)%desc = 'instantaneous canopy evaporation'
    ExtDiag(idx)%unit = 'W m-2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%evcw(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'trans'
    ExtDiag(idx)%desc = 'instantaneous transpiration'
    ExtDiag(idx)%unit = 'W m-2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%trans(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    if (Model%lsm == Model%lsm_ruc) then

      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'sfalb'
      ExtDiag(idx)%desc = 'surface albedo over land'
      ExtDiag(idx)%unit = 'fraction'
      ExtDiag(idx)%mod_name = 'gfs_sfc'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => Sfcprop%sfalb_lnd(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo

      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'rhofr'
      ExtDiag(idx)%desc = 'density of frozen precipitation'
      ExtDiag(idx)%unit = 'kg m-3'
      ExtDiag(idx)%mod_name = 'gfs_sfc'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => Sfcprop%rhofr(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo

      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'snowfall_acc_land'
      ExtDiag(idx)%desc = 'total accumulated frozen precipitation over land'
      ExtDiag(idx)%unit = 'kg m-2'
      ExtDiag(idx)%mod_name = 'gfs_sfc'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => Sfcprop%snowfallac_land(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo

      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'acsnow_land'
      ExtDiag(idx)%desc = 'total accumulated SWE of frozen precipitation over land'
      ExtDiag(idx)%unit = 'kg m-2'
      ExtDiag(idx)%mod_name = 'gfs_sfc'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => Sfcprop%acsnow_land(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo

      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'snowmt_land'
      ExtDiag(idx)%desc = 'accumulated snow melt over land'
      ExtDiag(idx)%unit = 'kg m-2'
      ExtDiag(idx)%mod_name = 'gfs_phys'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => IntDiag%snowmt_land(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo

      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'snowfall_acc_ice'
      ExtDiag(idx)%desc = 'total accumulated frozen precipitation over ice'
      ExtDiag(idx)%unit = 'kg m-2'
      ExtDiag(idx)%mod_name = 'gfs_sfc'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => Sfcprop%snowfallac_ice(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo

      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'acsnow_ice'
      ExtDiag(idx)%desc = 'total accumulated SWE of frozen precipitation over ice'
      ExtDiag(idx)%unit = 'kg m-2'
      ExtDiag(idx)%mod_name = 'gfs_sfc'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => Sfcprop%acsnow_ice(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo

      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'snowmt_ice'
      ExtDiag(idx)%desc = 'accumulated snow melt over ice'
      ExtDiag(idx)%unit = 'kg m-2'
      ExtDiag(idx)%mod_name = 'gfs_phys'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => IntDiag%snowmt_ice(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo
    endif ! RUC lsm

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'crain'
    ExtDiag(idx)%desc = 'instantaneous categorical rain'
    ExtDiag(idx)%unit = 'number'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    ExtDiag(idx)%cnvfac = cn_one
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => Sfcprop%srflag(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'stype'
    ExtDiag(idx)%desc = 'soil type in integer 1-9'
    ExtDiag(idx)%unit = 'number'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%int2 => Sfcprop%stype(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'scolor'
    ExtDiag(idx)%desc = 'soil color in integer 1-20'
    ExtDiag(idx)%unit = 'number'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%int2 => Sfcprop%scolor(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo


    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'lfrac'
    ExtDiag(idx)%desc = 'land fraction'
    ExtDiag(idx)%unit = 'fraction [0:1]'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => Sfcprop%landfrac(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'q2m'
    ExtDiag(idx)%desc = '2m specific humidity'
    ExtDiag(idx)%unit = 'kg/kg'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => Sfcprop%q2m(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 't2m'
    ExtDiag(idx)%desc = '2m temperature'
    ExtDiag(idx)%unit = 'K'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    ExtDiag(idx)%intpl_method = 'bilinear'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => Sfcprop%t2m(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'tsfc'
    ExtDiag(idx)%desc = 'surface temperature'
    ExtDiag(idx)%unit = 'K'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => Sfcprop%tsfc(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'usfco'
    ExtDiag(idx)%desc = 'surface zonal current'
    ExtDiag(idx)%unit = 'm/s'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => Sfcprop%usfco(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'vsfco'
    ExtDiag(idx)%desc = 'surface meridional current'
    ExtDiag(idx)%unit = 'm/s'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => Sfcprop%vsfco(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    if (Model%frac_grid) then
      do num = 1,Model%kice
        write (xtra,'(i1)') num
        idx = idx + 1
        ExtDiag(idx)%axes = 2
        ExtDiag(idx)%name = 'tiice'//trim(xtra)
        ExtDiag(idx)%desc = 'internal ice temperature layer ' // trim(xtra)
        ExtDiag(idx)%unit = 'K'
        ExtDiag(idx)%mod_name = 'gfs_sfc'
        allocate (ExtDiag(idx)%data(nblks))
        do nb = 1,nblks
          ExtDiag(idx)%data(nb)%var2 => Sfcprop%tiice(Model%chunk_begin(nb):Model%chunk_end(nb),num)
        enddo
      enddo
    end if

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'tg3'
    ExtDiag(idx)%desc = 'deep soil temperature'
    ExtDiag(idx)%unit = 'K'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => Sfcprop%tg3(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'tisfc'
    ExtDiag(idx)%desc = 'surface temperature over ice fraction'
    ExtDiag(idx)%unit = 'K'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => Sfcprop%tisfc(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'tprcp'
    ExtDiag(idx)%desc = 'total time-step precipitation'
    ExtDiag(idx)%unit = 'm'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => Sfcprop%tprcp(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'vtype'
    ExtDiag(idx)%desc = 'vegetation type in integer'
    ExtDiag(idx)%unit = 'number'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%int2 => sfcprop%vtype(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'weasd'
    ExtDiag(idx)%desc = 'surface snow water equivalent'
    ExtDiag(idx)%unit = 'kg/m**2'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => sfcprop%weasd(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'weasdi'
    ExtDiag(idx)%desc = 'surface snow water equivalent over ice'
    ExtDiag(idx)%unit = 'kg/m**2'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => sfcprop%weasdi(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'snodi'
    ExtDiag(idx)%desc = 'snow depth over ice'
    ExtDiag(idx)%unit = 'mm'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => sfcprop%snodi(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'hgtsfc'
    ExtDiag(idx)%desc = 'surface geopotential height'
    ExtDiag(idx)%unit = 'gpm'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    ExtDiag(idx)%cnvfac = cn_one
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => sfcprop%oro(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'slmsksfc'
    ExtDiag(idx)%desc = 'sea-land-ice mask (0-sea, 1-land, 2-ice)'
    ExtDiag(idx)%unit = 'numerical'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => sfcprop%slmsk(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'zorlsfc'
    ExtDiag(idx)%desc = 'surface roughness'
    ExtDiag(idx)%unit = 'm'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    ExtDiag(idx)%cnvfac = cn_one/cn_100
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => sfcprop%zorl(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'vfracsfc'
    ExtDiag(idx)%desc = 'vegetation fraction'
    ExtDiag(idx)%unit = 'fraction'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    ExtDiag(idx)%cnvfac = cn_100
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => sfcprop%vfrac(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    if (Model%lsm==Model%lsm_ruc) then
      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'wetness'
      ExtDiag(idx)%desc = 'soil moisture availability in top soil layer'
      ExtDiag(idx)%unit = 'fraction'
      ExtDiag(idx)%mod_name = 'gfs_sfc'
      ExtDiag(idx)%cnvfac = cn_100
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => sfcprop%wetness(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo
    end if

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'nirbmdi'
    ExtDiag(idx)%desc = 'sfc nir beam sw downward flux'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => Coupling%nirbmdi(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'nirdfdi'
    ExtDiag(idx)%desc = 'sfc nir diff sw downward flux'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => Coupling%nirdfdi(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'visbmdi'
    ExtDiag(idx)%desc = 'sfc uv+vis beam sw downward flux'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => Coupling%visbmdi(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'visdfdi'
    ExtDiag(idx)%desc = ' sfc uv+vis diff sw downward flux'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => Coupling%visdfdi(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'xlaixy'
    ExtDiag(idx)%desc = 'leaf area index'
    ExtDiag(idx)%unit = 'number'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => sfcprop%xlaixy(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    do num = 1,Model%nvegcat
      write (xtra,'(i2)') num
      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'vfrac_'//trim(xtra)
      ExtDiag(idx)%desc = 'fraction of vegetation category'//trim(xtra)
      ExtDiag(idx)%unit = 'frac'
      ExtDiag(idx)%mod_name = 'gfs_sfc'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => sfcprop%vegtype_frac(Model%chunk_begin(nb):Model%chunk_end(nb),num)
      enddo
    enddo

    do num = 1,Model%nsoilcat
      write (xtra,'(i2)') num
      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'sfrac_'//trim(xtra)
      ExtDiag(idx)%desc = 'fraction of soil category'//trim(xtra)
      ExtDiag(idx)%unit = 'frac'
      ExtDiag(idx)%mod_name = 'gfs_sfc'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => sfcprop%soiltype_frac(Model%chunk_begin(nb):Model%chunk_end(nb),num)
      enddo
    enddo

  if (Model%lsm == Model%lsm_ruc) then
    do num = 1,Model%lsoil_lsm
      write (xtra,'(i1)') num
      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'slc_'//trim(xtra)
      ExtDiag(idx)%desc = 'liquid soil moisture ' // trim(soil_layer_depth(Model%lsm, Model%lsm_ruc, Model%lsm_noah, num))
      ExtDiag(idx)%unit = 'm**3/m**3'
      ExtDiag(idx)%mod_name = 'gfs_sfc'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => sfcprop%sh2o(Model%chunk_begin(nb):Model%chunk_end(nb),num)
      enddo
    enddo
    idx = idx + 1
    ExtDiag(idx)%axes = 3
    ExtDiag(idx)%name = 'soill'
    ExtDiag(idx)%desc = 'liquid soil moisture'
    ExtDiag(idx)%unit = 'm**3/m**3'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var3 => sfcprop%sh2o(Model%chunk_begin(nb):Model%chunk_end(nb),:)
    enddo
  else
    do num = 1,Model%lsoil_lsm
      write (xtra,'(i1)') num
      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'slc_'//trim(xtra)
! DH* Can't use correct unit/description because of the way
! bit for bit tests are conducted (using cmp -> test fails)
#if 0
      ExtDiag(idx)%desc = 'liquid soil moisture ' // trim(soil_layer_depth(Model%lsm, Model%lsm_ruc, Model%lsm_noah, num))
      ExtDiag(idx)%unit = 'm**3/m**3'
#else
      ExtDiag(idx)%desc = 'liquid soil mositure at layer-'//trim(xtra)
      ExtDiag(idx)%unit = 'xxx'
#endif
! *DH
      ExtDiag(idx)%mod_name = 'gfs_sfc'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => sfcprop%slc(Model%chunk_begin(nb):Model%chunk_end(nb),num)
      enddo
    enddo
    idx = idx + 1
    ExtDiag(idx)%axes = 3
    ExtDiag(idx)%name = 'soill'
    ExtDiag(idx)%desc = 'liquid soil moisture'
    ExtDiag(idx)%unit = 'm**3/m**3'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var3 => sfcprop%slc(Model%chunk_begin(nb):Model%chunk_end(nb),:)
    enddo
  endif

  if (Model%lsm == Model%lsm_ruc) then
    do num = 1,Model%lsoil_lsm
      write (xtra,'(i1)') num
      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'soilw'//trim(xtra)
      ExtDiag(idx)%desc = 'volumetric soil moisture ' // trim(soil_layer_depth(Model%lsm, Model%lsm_ruc, Model%lsm_noah, num))
      ExtDiag(idx)%unit = 'fraction'
      ExtDiag(idx)%mod_name = 'gfs_sfc'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => sfcprop%smois(Model%chunk_begin(nb):Model%chunk_end(nb),num)
      enddo
    enddo
    idx = idx + 1
    ExtDiag(idx)%axes = 3
    ExtDiag(idx)%name = 'soilw'
    ExtDiag(idx)%desc = 'volumetric soil moisture'
    ExtDiag(idx)%unit = 'fraction'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var3 => sfcprop%smois(Model%chunk_begin(nb):Model%chunk_end(nb),:)
    enddo
  else
    do num = 1,Model%lsoil_lsm
      write (xtra,'(i1)') num
      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'soilw'//trim(xtra)
      ExtDiag(idx)%desc = 'volumetric soil moisture ' // trim(soil_layer_depth(Model%lsm, Model%lsm_ruc, Model%lsm_noah, num))
      ExtDiag(idx)%unit = 'fraction'
      ExtDiag(idx)%mod_name = 'gfs_sfc'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => sfcprop%smc(Model%chunk_begin(nb):Model%chunk_end(nb),num)
      enddo
    enddo
    idx = idx + 1
    ExtDiag(idx)%axes = 3
    ExtDiag(idx)%name = 'soilw'
    ExtDiag(idx)%desc = 'volumetric soil moisture'
    ExtDiag(idx)%unit = 'fraction'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var3 => sfcprop%smc(Model%chunk_begin(nb):Model%chunk_end(nb),:)
    enddo
  endif

    if (Model%lsm == Model%lsm_ruc) then
      do num = 1,Model%lsoil_lsm
        write (xtra,'(i1)') num
        idx = idx + 1
        ExtDiag(idx)%axes = 2
        ExtDiag(idx)%name = 'soilt'//trim(xtra)
        ExtDiag(idx)%desc = 'soil temperature ' // trim(soil_layer_depth(Model%lsm, Model%lsm_ruc, Model%lsm_noah, num))
        ExtDiag(idx)%unit = 'K'
        ExtDiag(idx)%mod_name = 'gfs_sfc'
        allocate (ExtDiag(idx)%data(nblks))
        do nb = 1,nblks
          ExtDiag(idx)%data(nb)%var2 => sfcprop%tslb(Model%chunk_begin(nb):Model%chunk_end(nb),num)
        enddo
      enddo
      idx = idx + 1
      ExtDiag(idx)%axes = 3
      ExtDiag(idx)%name = 'soilt'
      ExtDiag(idx)%desc = 'soil temperature'
      ExtDiag(idx)%unit = 'K'
      ExtDiag(idx)%mod_name = 'gfs_sfc'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var3 => sfcprop%tslb(Model%chunk_begin(nb):Model%chunk_end(nb),:)
      enddo
    else
      do num = 1,Model%lsoil_lsm
        write (xtra,'(i1)') num
        idx = idx + 1
        ExtDiag(idx)%axes = 2
        ExtDiag(idx)%name = 'soilt'//trim(xtra)
        ExtDiag(idx)%desc = 'soil temperature ' // trim(soil_layer_depth(Model%lsm, Model%lsm_ruc, Model%lsm_noah, num))
        ExtDiag(idx)%unit = 'K'
        ExtDiag(idx)%mod_name = 'gfs_sfc'
        allocate (ExtDiag(idx)%data(nblks))
        do nb = 1,nblks
          ExtDiag(idx)%data(nb)%var2 => sfcprop%stc(Model%chunk_begin(nb):Model%chunk_end(nb),num)
        enddo
      enddo
      idx = idx + 1
      ExtDiag(idx)%axes = 3
      ExtDiag(idx)%name = 'soilt'
      ExtDiag(idx)%desc = 'soil temperature'
      ExtDiag(idx)%unit = 'K'
      ExtDiag(idx)%mod_name = 'gfs_sfc'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var3 => sfcprop%stc(Model%chunk_begin(nb):Model%chunk_end(nb),:)
      enddo
    endif

!--------------------------nsst variables
  if (model%nstf_name(1) > 0) then
!--------------------------nsst variables

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'tref'
    ExtDiag(idx)%desc = 'nsst reference or foundation temperature'
    ExtDiag(idx)%unit = 'K'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => sfcprop%tref(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'z_c'
    ExtDiag(idx)%desc = 'nsst sub-layer cooling thickness'
    ExtDiag(idx)%unit = 'm'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => sfcprop%z_c(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'c_0'
    ExtDiag(idx)%desc = 'nsst coefficient1 to calculate d(tz)/d(ts)'
    ExtDiag(idx)%unit = 'numerical'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => sfcprop%c_0(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'c_d'
    ExtDiag(idx)%desc = 'nsst coefficient2 to calculate d(tz)/d(ts)'
    ExtDiag(idx)%unit = 'n/a'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => sfcprop%c_d(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'w_0'
    ExtDiag(idx)%desc = 'nsst coefficient3 to calculate d(tz)/d(ts)'
    ExtDiag(idx)%unit = 'n/a'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => sfcprop%w_0(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'w_d'
    ExtDiag(idx)%desc = 'nsst coefficient4 to calculate d(tz)/d(ts)'
    ExtDiag(idx)%unit = 'n/a'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => sfcprop%w_d(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'xt'
    ExtDiag(idx)%desc = 'nsst heat content in diurnal thermocline layer'
    ExtDiag(idx)%unit = 'k*m'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => sfcprop%xt(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'xs'
    ExtDiag(idx)%desc = 'nsst salinity content in diurnal thermocline layer'
    ExtDiag(idx)%unit = 'n/a'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => sfcprop%xs(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'xu'
    ExtDiag(idx)%desc = 'nsst u-current content in diurnal thermocline layer'
    ExtDiag(idx)%unit = 'm2/s'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => sfcprop%xu(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'xv'
    ExtDiag(idx)%desc = 'nsst v-current content in diurnal thermocline layer'
    ExtDiag(idx)%unit = 'm2/s'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => sfcprop%xv(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'xz'
    ExtDiag(idx)%desc = 'nsst diurnal thermocline layer thickness'
    ExtDiag(idx)%unit = 'm'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => sfcprop%xz(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'zm'
    ExtDiag(idx)%desc = 'nsst mixed layer thickness'
    ExtDiag(idx)%unit = 'm'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => sfcprop%zm(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'xtts'
    ExtDiag(idx)%desc = 'nsst d(xt)/d(ts)'
    ExtDiag(idx)%unit = 'm'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => sfcprop%xtts(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'xzts'
    ExtDiag(idx)%desc = 'nsst d(xt)/d(ts)'
    ExtDiag(idx)%unit = 'm/k'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => sfcprop%xzts(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'd_conv'
    ExtDiag(idx)%desc = 'nsst thickness of free convection layer'
    ExtDiag(idx)%unit = 'm'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => sfcprop%d_conv(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'ifd'
    ExtDiag(idx)%desc = 'nsst index to start dtlm run or not'
    ExtDiag(idx)%unit = 'n/a'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => sfcprop%ifd(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'dt_cool'
    ExtDiag(idx)%desc = 'nsst sub-layer cooling amount'
    ExtDiag(idx)%unit = 'k'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => sfcprop%dt_cool(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo

    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'qrain'
    ExtDiag(idx)%desc = 'nsst sensible heat flux due to rainfall'
    ExtDiag(idx)%unit = 'W/m**2'
    ExtDiag(idx)%mod_name = 'gfs_sfc'
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => sfcprop%qrain(Model%chunk_begin(nb):Model%chunk_end(nb))
    enddo
!--------------------------nsst variables
  endif

!--------------------------aerosols
    if (Model%ntwa>0) then
      if (Model%ltaerosol) then
        idx = idx + 1
        ExtDiag(idx)%axes = 3
        ExtDiag(idx)%name = 'nwfa'
        ExtDiag(idx)%desc = 'number concentration of water-friendly aerosols'
        ExtDiag(idx)%unit = 'kg-1'
        ExtDiag(idx)%mod_name = 'gfs_phys'
        allocate (ExtDiag(idx)%data(nblks))
        do nb = 1,nblks
          ExtDiag(idx)%data(nb)%var3 => Statein%qgrs(Model%chunk_begin(nb):Model%chunk_end(nb),:,Model%ntwa)
        enddo

        idx = idx + 1
        ExtDiag(idx)%axes = 2
        ExtDiag(idx)%name = 'nwfa2d'
        ExtDiag(idx)%desc = 'water-friendly surface aerosol source'
        ExtDiag(idx)%unit = 'kg-1 s-1'
        ExtDiag(idx)%mod_name = 'gfs_sfc'
        allocate (ExtDiag(idx)%data(nblks))
        do nb = 1,nblks
          ExtDiag(idx)%data(nb)%var2 => Coupling%nwfa2d(Model%chunk_begin(nb):Model%chunk_end(nb))
        enddo
      elseif (Model%mraerosol) then
        idx = idx + 1
        ExtDiag(idx)%axes = 3
        ExtDiag(idx)%name = 'nwfa'
        ExtDiag(idx)%desc = 'number concentration of water-friendly aerosols'
        ExtDiag(idx)%unit = 'kg-1'
        ExtDiag(idx)%mod_name = 'gfs_phys'
        allocate (ExtDiag(idx)%data(nblks))
        do nb = 1,nblks
          ExtDiag(idx)%data(nb)%var3 => Stateout%gq0(Model%chunk_begin(nb):Model%chunk_end(nb),:,Model%ntwa)
        enddo
      endif
    endif

    if (Model%ntia>0) then
      if (Model%ltaerosol) then
        idx = idx + 1
        ExtDiag(idx)%axes = 3
        ExtDiag(idx)%name = 'nifa'
        ExtDiag(idx)%desc = 'number concentration of ice-friendly aerosols'
        ExtDiag(idx)%unit = 'kg-1'
        ExtDiag(idx)%mod_name = 'gfs_phys'
        allocate (ExtDiag(idx)%data(nblks))
        do nb = 1,nblks
          ExtDiag(idx)%data(nb)%var3 => Statein%qgrs(Model%chunk_begin(nb):Model%chunk_end(nb),:,Model%ntia)
        enddo

        idx = idx + 1
        ExtDiag(idx)%axes = 2
        ExtDiag(idx)%name = 'nifa2d'
        ExtDiag(idx)%desc = 'ice-friendly surface aerosol source'
        ExtDiag(idx)%unit = 'kg-1 s-1'
        ExtDiag(idx)%mod_name = 'gfs_sfc'
        allocate (ExtDiag(idx)%data(nblks))
        do nb = 1,nblks
          ExtDiag(idx)%data(nb)%var2 => Coupling%nifa2d(Model%chunk_begin(nb):Model%chunk_end(nb))
        enddo
      else if (Model%mraerosol) then
        idx = idx + 1
        ExtDiag(idx)%axes = 3
        ExtDiag(idx)%name = 'nifa'
        ExtDiag(idx)%desc = 'number concentration of ice-friendly aerosols'
        ExtDiag(idx)%unit = 'kg-1'
        ExtDiag(idx)%mod_name = 'gfs_phys'
        allocate (ExtDiag(idx)%data(nblks))
        do nb = 1,nblks
          ExtDiag(idx)%data(nb)%var3 => Stateout%gq0(Model%chunk_begin(nb):Model%chunk_end(nb),:,Model%ntia)
        enddo
      end if
    endif

    ! Extended diagnostics from Thompson MP
    thompson_extended_diagnostics: if (Model%ext_diag_thompson) then
      do num=1,Model%thompson_ext_ndiag3d
        idx = idx + 1
        ExtDiag(idx)%axes = 3
        select case (num)
          ! This is the place to add specific names, descriptions,
          ! and units if so desired
          !case (1)
          ! ...
          case default
            write (xtra,'(I2.2)') num
            ExtDiag(idx)%name = 'thompson_diag3d_' // trim(xtra)
            ExtDiag(idx)%desc = 'Thompson extended diagnostics array ' // trim(xtra)
            ExtDiag(idx)%unit = 'unknown'
        end select
        ExtDiag(idx)%mod_name = 'gfs_phys'
        ExtDiag(idx)%intpl_method = 'bilinear'
        ExtDiag(idx)%time_avg = .false.
        allocate (ExtDiag(idx)%data(nblks))
        do nb = 1,nblks
          ExtDiag(idx)%data(nb)%var3 => IntDiag%thompson_ext_diag3d(Model%chunk_begin(nb):Model%chunk_end(nb),:,num)
        enddo
      enddo
    end if thompson_extended_diagnostics

    if (Model%cpl_fire .and. Model%ntfsmoke>0) then
      idx = idx + 1
      ExtDiag(idx)%axes = 3
      ExtDiag(idx)%name = 'fsmoke'
      ExtDiag(idx)%desc = 'smoke concentration'
      ExtDiag(idx)%unit = 'kg kg-1'
      ExtDiag(idx)%mod_name = 'gfs_phys'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var3 => Statein%qgrs(Model%chunk_begin(nb):Model%chunk_end(nb),:,Model%ntfsmoke)
      enddo
    endif

    if (Model%rrfs_sd .and. Model%ntsmoke>0) then

      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'fire_heat'
      ExtDiag(idx)%desc = 'surface fire heat flux'
      ExtDiag(idx)%unit = 'W m-2'
      ExtDiag(idx)%mod_name = 'gfs_sfc'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => Sfcprop%fire_heat_flux(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo

      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'burned'
      ExtDiag(idx)%desc = 'ration of the burnt area to the grid cell area'
      ExtDiag(idx)%unit = 'frac'
      ExtDiag(idx)%mod_name = 'gfs_sfc'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => Sfcprop%frac_grid_burned(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo

      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'emdust'
      ExtDiag(idx)%desc = 'emission of fine dust for smoke'
      ExtDiag(idx)%unit = 'ug m-2 s-1'
      ExtDiag(idx)%mod_name = 'gfs_sfc'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => Sfcprop%emdust(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo

      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'emseas'
      ExtDiag(idx)%desc = 'emission of seas for smoke'
      ExtDiag(idx)%unit = 'ug m-2 s-1'
      ExtDiag(idx)%mod_name = 'gfs_sfc'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => Sfcprop%emseas(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo

      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'emanoc'
      ExtDiag(idx)%desc = 'emission of anoc for thompson mp'
      ExtDiag(idx)%unit = 'ug m-2 s-1'
      ExtDiag(idx)%mod_name = 'gfs_sfc'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => Sfcprop%emanoc(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo

      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'coef_bb_dc'
      ExtDiag(idx)%desc = 'coeff bb for smoke'
      ExtDiag(idx)%unit = ''
      ExtDiag(idx)%mod_name = 'gfs_sfc'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => Sfcprop%coef_bb_dc(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo

      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'min_fplume'
      ExtDiag(idx)%desc = 'minimum smoke plume height'
      ExtDiag(idx)%unit = ''
      ExtDiag(idx)%mod_name = 'gfs_sfc'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => Coupling%min_fplume(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo

      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'max_fplume'
      ExtDiag(idx)%desc = 'maximum smoke plume height'
      ExtDiag(idx)%unit = ''
      ExtDiag(idx)%mod_name = 'gfs_sfc'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => Coupling%max_fplume(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo

      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'HWP'
      ExtDiag(idx)%desc = 'hourly fire weather potential'
      ExtDiag(idx)%unit = ' '
      ExtDiag(idx)%mod_name = 'gfs_sfc'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => Coupling%rrfs_hwp(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo
      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'HWP_ave'
      ExtDiag(idx)%desc = 'averaged fire weather potential'
      ExtDiag(idx)%unit = ' '
      ExtDiag(idx)%mod_name = 'gfs_sfc'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => Coupling%rrfs_hwp_ave(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo

      extended_smoke_dust_diagnostics: if ( Model%extended_sd_diags ) then

      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'uspdavg'
      ExtDiag(idx)%desc = 'BL average wind speed'
      ExtDiag(idx)%unit = ''
      ExtDiag(idx)%mod_name = 'gfs_sfc'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => Coupling%uspdavg(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo

      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'hpbl_thetav'
      ExtDiag(idx)%desc = 'BL depth modified parcel method'
      ExtDiag(idx)%unit = ''
      ExtDiag(idx)%mod_name = 'gfs_sfc'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => Coupling%hpbl_thetav(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo

      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'drydep_smoke'
      ExtDiag(idx)%desc = 'dry deposition smoke'
      ExtDiag(idx)%unit = ' '
      ExtDiag(idx)%mod_name = 'gfs_sfc'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => Coupling%drydep_flux(Model%chunk_begin(nb):Model%chunk_end(nb),1)
      enddo

      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'drydep_dust'
      ExtDiag(idx)%desc = 'dry deposition dust'
      ExtDiag(idx)%unit = ' '
      ExtDiag(idx)%mod_name = 'gfs_sfc'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => Coupling%drydep_flux(Model%chunk_begin(nb):Model%chunk_end(nb),2)
      enddo

      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'drydep_coarsepm'
      ExtDiag(idx)%desc = 'dry deposition coarsepm'
      ExtDiag(idx)%unit = ' '
      ExtDiag(idx)%mod_name = 'gfs_sfc'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => Coupling%drydep_flux(Model%chunk_begin(nb):Model%chunk_end(nb),3)
      enddo

      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'wetdpr_smoke'
      ExtDiag(idx)%desc = 'resolved wet deposition smoke'
      ExtDiag(idx)%unit = ' '
      ExtDiag(idx)%mod_name = 'gfs_sfc'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => Coupling%wetdpr_flux(Model%chunk_begin(nb):Model%chunk_end(nb),1)
      enddo

      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'wetdpr_dust'
      ExtDiag(idx)%desc = 'resolved wet deposition dust'
      ExtDiag(idx)%unit = ' '
      ExtDiag(idx)%mod_name = 'gfs_sfc'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => Coupling%wetdpr_flux(Model%chunk_begin(nb):Model%chunk_end(nb),2)
      enddo

      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'wetdpr_coarsepm'
      ExtDiag(idx)%desc = 'resolved wet deposition coarsepm'
      ExtDiag(idx)%unit = ' '
      ExtDiag(idx)%mod_name = 'gfs_sfc'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => Coupling%wetdpr_flux(Model%chunk_begin(nb):Model%chunk_end(nb),3)
      enddo

      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'wetdpc_smoke'
      ExtDiag(idx)%desc = 'convective wet deposition smoke'
      ExtDiag(idx)%unit = ' '
      ExtDiag(idx)%mod_name = 'gfs_sfc'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => Coupling%wetdpc_flux(Model%chunk_begin(nb):Model%chunk_end(nb),1)
      enddo

      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'wetdpc_dust'
      ExtDiag(idx)%desc = 'convective wet deposition dust'
      ExtDiag(idx)%unit = ' '
      ExtDiag(idx)%mod_name = 'gfs_sfc'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => Coupling%wetdpc_flux(Model%chunk_begin(nb):Model%chunk_end(nb),2)
      enddo

      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'wetdpc_coarsepm'
      ExtDiag(idx)%desc = 'convective wet deposition coarsepm'
      ExtDiag(idx)%unit = ' '
      ExtDiag(idx)%mod_name = 'gfs_sfc'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => Coupling%wetdpc_flux(Model%chunk_begin(nb):Model%chunk_end(nb),3)
      enddo

      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'peak_hr'
      ExtDiag(idx)%desc = 'hour of peak smoke emissions'
      ExtDiag(idx)%unit = ' '
      ExtDiag(idx)%mod_name = 'gfs_sfc'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
       ExtDiag(idx)%data(nb)%var2 => Sfcprop%peak_hr(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo

      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'fire_type'
      ExtDiag(idx)%desc = 'fire type'
      ExtDiag(idx)%unit = ''
      ExtDiag(idx)%mod_name = 'gfs_sfc'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%int2 => Sfcprop%fire_type(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo

      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'lu_nofire'
      ExtDiag(idx)%desc = 'lu nofire pixes'
      ExtDiag(idx)%unit = ''
      ExtDiag(idx)%mod_name = 'gfs_sfc'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => Sfcprop%lu_nofire(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo

      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'lu_qfire'
      ExtDiag(idx)%desc = 'lu qfire pixes'
      ExtDiag(idx)%unit = ''
      ExtDiag(idx)%mod_name = 'gfs_sfc'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => Sfcprop%lu_qfire(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo

      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'fhist'
      ExtDiag(idx)%desc = 'coefficient to scale the fire activity depending on the fire duration'
      ExtDiag(idx)%unit = ' '
      ExtDiag(idx)%mod_name = 'gfs_sfc'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
       ExtDiag(idx)%data(nb)%var2 => Sfcprop%fhist(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo

      if (Model%ebb_dcycle == 2 ) then

      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'fire_end_hr'
      ExtDiag(idx)%desc = 'Hours since fire was last detected'
      ExtDiag(idx)%unit = 'hrs'
      ExtDiag(idx)%mod_name = 'gfs_sfc'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
       ExtDiag(idx)%data(nb)%var2 => Sfcprop%smoke2d_RRFS(Model%chunk_begin(nb):Model%chunk_end(nb),3)
      enddo

      endif

      endif  extended_smoke_dust_diagnostics


      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'ebb_smoke_in'
      ExtDiag(idx)%desc = 'input smoke emission'
      ExtDiag(idx)%unit = 'ug m-2 s-1'
      ExtDiag(idx)%mod_name = 'gfs_sfc'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
       ExtDiag(idx)%data(nb)%var2 => Sfcprop%ebb_smoke_in(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo


      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'frp_output'
      ExtDiag(idx)%desc = 'output frp'
      ExtDiag(idx)%unit = 'mw'
      ExtDiag(idx)%mod_name = 'gfs_sfc'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
       ExtDiag(idx)%data(nb)%var2 => Sfcprop%frp_output(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo

      smoke_forecast_mode: if (Model%ebb_dcycle == 2 ) then

      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'ebb_rate'
      ExtDiag(idx)%desc = 'Total EBB Emissions'
      ExtDiag(idx)%unit = 'ug m-2 s-1'
      ExtDiag(idx)%mod_name = 'gfs_sfc'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
       ExtDiag(idx)%data(nb)%var2 => Sfcprop%smoke2d_RRFS(Model%chunk_begin(nb):Model%chunk_end(nb),1)
      enddo

      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'frp_davg'
      ExtDiag(idx)%desc = 'Daily mean Fire Radiative Power'
      ExtDiag(idx)%unit = 'mw'
      ExtDiag(idx)%mod_name = 'gfs_sfc'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
       ExtDiag(idx)%data(nb)%var2 => Sfcprop%smoke2d_RRFS(Model%chunk_begin(nb):Model%chunk_end(nb),2)
      enddo


      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'hwp_davg'
      ExtDiag(idx)%desc = 'Daily mean Hourly Wildfire Potential'
      ExtDiag(idx)%unit = ' '
      ExtDiag(idx)%mod_name = 'gfs_sfc'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
       ExtDiag(idx)%data(nb)%var2 => Sfcprop%smoke2d_RRFS(Model%chunk_begin(nb):Model%chunk_end(nb),4)
      enddo

      endif smoke_forecast_mode

      idx = idx + 1
      ExtDiag(idx)%axes = 3
      ExtDiag(idx)%name = 'ebu_smoke'
      ExtDiag(idx)%desc = 'smoke emission'
      ExtDiag(idx)%unit = 'ug/m2/s'
      ExtDiag(idx)%mod_name = 'gfs_phys'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
       ExtDiag(idx)%data(nb)%var3 => Coupling%ebu_smoke(Model%chunk_begin(nb):Model%chunk_end(nb),:)
      enddo

      idx = idx + 1
      ExtDiag(idx)%axes = 3
      ExtDiag(idx)%name = 'ext550'
      ExtDiag(idx)%desc = '3d total extinction at 550nm'
      ExtDiag(idx)%unit = ' '
      ExtDiag(idx)%mod_name = 'gfs_phys'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
       ExtDiag(idx)%data(nb)%var3 => Radtend%ext550(Model%chunk_begin(nb):Model%chunk_end(nb),:)
      enddo
    endif

    do i=1,Model%num_dfi_radar
       idx = idx + 1
       ExtDiag(idx)%axes = 3
       if(i>1) then
          write(ExtDiag(idx)%name,'(A,I0)') 'radar_tten_',i
       else
          ExtDiag(idx)%name = 'radar_tten'
       endif
       write(ExtDiag(idx)%desc,'(A,I0,A,I0)') 'temperature tendency due to dfi radar tendencies ',i,' of ',Model%num_dfi_radar
       ExtDiag(idx)%unit = 'K s-1'
       ExtDiag(idx)%mod_name = 'gfs_phys'
       ExtDiag(idx)%time_avg = .FALSE.

       allocate (ExtDiag(idx)%data(nblks))
       do nb = 1,nblks
          ExtDiag(idx)%data(nb)%var3 => Tbd%dfi_radar_tten(Model%chunk_begin(nb):Model%chunk_end(nb),:,i)
       enddo
    enddo

    if(Model%lightning_threat) then
      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'ltg1_max'
      ExtDiag(idx)%desc = 'Max Lightning Threat 1'
      ExtDiag(idx)%unit = 'flashes/(5 min)'
      ExtDiag(idx)%mod_name = 'gfs_sfc'
      ! CCPP physics units are flashes per minute
      ExtDiag(idx)%cnvfac = 5.0_kind_phys
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => IntDiag%ltg1_max(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo

      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'ltg2_max'
      ExtDiag(idx)%desc = 'Max Lightning Threat 2'
      ExtDiag(idx)%unit = 'flashes/(5 min)'
      ExtDiag(idx)%mod_name = 'gfs_sfc'
      ! CCPP physics units are flashes per minute
      ExtDiag(idx)%cnvfac = 5.0_kind_phys
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => IntDiag%ltg2_max(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo

      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'ltg3_max'
      ExtDiag(idx)%desc = 'Max Lightning Threat 3'
      ExtDiag(idx)%unit = 'flashes/(5 min)'
      ExtDiag(idx)%mod_name = 'gfs_sfc'
      ! CCPP physics units are flashes per minute
      ExtDiag(idx)%cnvfac = 5.0_kind_phys
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => IntDiag%ltg3_max(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo
    endif

    ! Cloud effective radii from Microphysics
    if (Model%imp_physics == Model%imp_physics_thompson .or. Model%imp_physics == Model%imp_physics_fer_hires .or. &
        Model%imp_physics == Model%imp_physics_nssl  ) then
      idx = idx + 1
      ExtDiag(idx)%axes = 3
      ExtDiag(idx)%name = 'cleffr'
      ExtDiag(idx)%desc = 'effective radius of cloud liquid water particle'
      ExtDiag(idx)%unit = 'um'
      ExtDiag(idx)%mod_name = 'gfs_phys'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var3 => Tbd%phy_f3d(Model%chunk_begin(nb):Model%chunk_end(nb),:,Model%nleffr)
      enddo
      idx = idx + 1
      ExtDiag(idx)%axes = 3
      ExtDiag(idx)%name = 'cieffr'
      ExtDiag(idx)%desc = 'effective radius of stratiform cloud ice particle in um'
      ExtDiag(idx)%unit = 'um'
      ExtDiag(idx)%mod_name = 'gfs_phys'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var3 => Tbd%phy_f3d(Model%chunk_begin(nb):Model%chunk_end(nb),:,Model%nieffr)
      enddo
      idx = idx + 1
      ExtDiag(idx)%axes = 3
      ExtDiag(idx)%name = 'cseffr'
      ExtDiag(idx)%desc = 'effective radius of stratiform cloud snow particle in um'
      ExtDiag(idx)%unit = 'um'
      ExtDiag(idx)%mod_name = 'gfs_phys'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var3 => Tbd%phy_f3d(Model%chunk_begin(nb):Model%chunk_end(nb),:,Model%nseffr)
      enddo
    endif

    !MYNN
    if (Model%do_mynnedmf) then

      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'ztop_plume'
      ExtDiag(idx)%desc = 'height of highest plume'
      ExtDiag(idx)%unit = 'm'
      ExtDiag(idx)%mod_name = 'gfs_sfc'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => IntDiag%ztop_plume(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo

      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'maxmf'
      ExtDiag(idx)%desc = 'maximum mass-flux in column'
      ExtDiag(idx)%unit = 'm s-1'
      ExtDiag(idx)%mod_name = 'gfs_sfc'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => IntDiag%maxmf(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo

      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'maxwidth'
      ExtDiag(idx)%desc = 'maximum width of plumes in grid column'
      ExtDiag(idx)%unit = 'm'
      ExtDiag(idx)%mod_name = 'gfs_sfc'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => IntDiag%maxwidth(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo
    endif

    if (Model%do_mynnsfclay) then
      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'zol'
      ExtDiag(idx)%desc = 'monin obukhov surface stability parameter'
      ExtDiag(idx)%unit = 'n/a'
      ExtDiag(idx)%mod_name = 'gfs_sfc'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => Sfcprop%zol(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo

      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'flhc'
      ExtDiag(idx)%desc = 'surface exchange coefficient for heat'
      ExtDiag(idx)%unit = 'W m-2 K-1'
      ExtDiag(idx)%mod_name = 'gfs_sfc'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => Sfcprop%flhc(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo

      idx = idx + 1
      ExtDiag(idx)%axes = 2
      ExtDiag(idx)%name = 'flqc'
      ExtDiag(idx)%desc = 'surface exchange coefficient for moisture'
      ExtDiag(idx)%unit = 'kg m-2 s-1'
      ExtDiag(idx)%mod_name = 'gfs_sfc'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var2 => Sfcprop%flqc(Model%chunk_begin(nb):Model%chunk_end(nb))
      enddo
    endif

    if (Model%do_mynnedmf) then
      idx = idx + 1
      ExtDiag(idx)%axes = 3
      ExtDiag(idx)%name = 'CLDFRA_BL'
      ExtDiag(idx)%desc = 'subgrid cloud fraction'
      ExtDiag(idx)%unit = 'frac'
      ExtDiag(idx)%mod_name = 'gfs_phys'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var3 => Tbd%CLDFRA_BL(Model%chunk_begin(nb):Model%chunk_end(nb),:)
      enddo

      idx = idx + 1
      ExtDiag(idx)%axes = 3
      ExtDiag(idx)%name = 'QC_BL'
      ExtDiag(idx)%desc = 'subgrid cloud mixing ratio'
      ExtDiag(idx)%unit = 'frac'
      ExtDiag(idx)%mod_name = 'gfs_phys'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var3 => Tbd%QC_BL(Model%chunk_begin(nb):Model%chunk_end(nb),:)
      enddo

      idx = idx + 1
      ExtDiag(idx)%axes = 3
      ExtDiag(idx)%name = 'EL_PBL'
      ExtDiag(idx)%desc = 'turbulent mixing length'
      ExtDiag(idx)%unit = 'm'
      ExtDiag(idx)%mod_name = 'gfs_phys'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var3 => Tbd%el_pbl(Model%chunk_begin(nb):Model%chunk_end(nb),:)
      enddo

      idx = idx + 1
      ExtDiag(idx)%axes = 3
      ExtDiag(idx)%name = 'QKE'
      ExtDiag(idx)%desc = '2 X TKE (from mynn)'
      ExtDiag(idx)%unit = 'm2 s-2'
      ExtDiag(idx)%mod_name = 'gfs_phys'
      allocate (ExtDiag(idx)%data(nblks))
      do nb = 1,nblks
        ExtDiag(idx)%data(nb)%var3 => Tbd%QKE(Model%chunk_begin(nb):Model%chunk_end(nb),:)
      enddo

      if (Model%bl_mynn_output > 0) then

        idx = idx + 1
        ExtDiag(idx)%axes = 3
        ExtDiag(idx)%name = 'edmf_a'
        ExtDiag(idx)%desc = 'updraft area fraction (from mynn)'
        ExtDiag(idx)%unit = '-'
        ExtDiag(idx)%mod_name = 'gfs_phys'
        allocate (ExtDiag(idx)%data(nblks))
        do nb = 1,nblks
          ExtDiag(idx)%data(nb)%var3 => IntDiag%edmf_a(Model%chunk_begin(nb):Model%chunk_end(nb),:)
        enddo

        idx = idx + 1
        ExtDiag(idx)%axes = 3
        ExtDiag(idx)%name = 'edmf_w'
        ExtDiag(idx)%desc = 'mean updraft vertical veloctity (mynn)'
        ExtDiag(idx)%unit = 'm s-1'
        ExtDiag(idx)%mod_name = 'gfs_phys'
        allocate (ExtDiag(idx)%data(nblks))
        do nb = 1,nblks
          ExtDiag(idx)%data(nb)%var3 => IntDiag%edmf_w(Model%chunk_begin(nb):Model%chunk_end(nb),:)
        enddo

        idx = idx + 1
        ExtDiag(idx)%axes = 3
        ExtDiag(idx)%name = 'edmf_qt'
        ExtDiag(idx)%desc = 'updraft total water (from mynn)'
        ExtDiag(idx)%unit = 'kg kg-1'
        ExtDiag(idx)%mod_name = 'gfs_phys'
        allocate (ExtDiag(idx)%data(nblks))
        do nb = 1,nblks
          ExtDiag(idx)%data(nb)%var3 => IntDiag%edmf_qt(Model%chunk_begin(nb):Model%chunk_end(nb),:)
        enddo

        idx = idx + 1
        ExtDiag(idx)%axes = 3
        ExtDiag(idx)%name = 'edmf_thl'
        ExtDiag(idx)%desc = 'mean liquid potential temperature (mynn)'
        ExtDiag(idx)%unit = 'K'
        ExtDiag(idx)%mod_name = 'gfs_phys'
        allocate (ExtDiag(idx)%data(nblks))
        do nb = 1,nblks
          ExtDiag(idx)%data(nb)%var3 => IntDiag%edmf_thl(Model%chunk_begin(nb):Model%chunk_end(nb),:)
        enddo

        idx = idx + 1
        ExtDiag(idx)%axes = 3
        ExtDiag(idx)%name = 'edmf_ent'
        ExtDiag(idx)%desc = 'updraft entrainment rate (from mynn)'
        ExtDiag(idx)%unit = 'm-1'
        ExtDiag(idx)%mod_name = 'gfs_phys'
        allocate (ExtDiag(idx)%data(nblks))
        do nb = 1,nblks
          ExtDiag(idx)%data(nb)%var3 => IntDiag%edmf_ent(Model%chunk_begin(nb):Model%chunk_end(nb),:)
        enddo

        idx = idx + 1
        ExtDiag(idx)%axes = 3
        ExtDiag(idx)%name = 'edmf_qc'
        ExtDiag(idx)%desc = 'mean updraft liquid water (mynn)'
        ExtDiag(idx)%unit = 'kg kg-1'
        ExtDiag(idx)%mod_name = 'gfs_phys'
        allocate (ExtDiag(idx)%data(nblks))
        do nb = 1,nblks
          ExtDiag(idx)%data(nb)%var3 => IntDiag%edmf_qc(Model%chunk_begin(nb):Model%chunk_end(nb),:)
        enddo

        idx = idx + 1
        ExtDiag(idx)%axes = 3
        ExtDiag(idx)%name = 'sub_thl'
        ExtDiag(idx)%desc = 'subsidence temperature tendency (from mynn)'
        ExtDiag(idx)%unit = 'K s-1'
        ExtDiag(idx)%mod_name = 'gfs_phys'
        allocate (ExtDiag(idx)%data(nblks))
        do nb = 1,nblks
          ExtDiag(idx)%data(nb)%var3 => IntDiag%sub_thl(Model%chunk_begin(nb):Model%chunk_end(nb),:)
        enddo

        idx = idx + 1
        ExtDiag(idx)%axes = 3
        ExtDiag(idx)%name = 'sub_sqv'
        ExtDiag(idx)%desc = 'subsidence water vapor tendency (mynn)'
        ExtDiag(idx)%unit = 'kg kg-1 s-1'
        ExtDiag(idx)%mod_name = 'gfs_phys'
        allocate (ExtDiag(idx)%data(nblks))
        do nb = 1,nblks
          ExtDiag(idx)%data(nb)%var3 => IntDiag%sub_sqv(Model%chunk_begin(nb):Model%chunk_end(nb),:)
        enddo

        idx = idx + 1
        ExtDiag(idx)%axes = 3
        ExtDiag(idx)%name = 'det_thl'
        ExtDiag(idx)%desc = 'detrainment temperature tendency (from mynn)'
        ExtDiag(idx)%unit = 'K s-1'
        ExtDiag(idx)%mod_name = 'gfs_phys'
        allocate (ExtDiag(idx)%data(nblks))
        do nb = 1,nblks
          ExtDiag(idx)%data(nb)%var3 => IntDiag%det_thl(Model%chunk_begin(nb):Model%chunk_end(nb),:)
        enddo

        idx = idx + 1
        ExtDiag(idx)%axes = 3
        ExtDiag(idx)%name = 'det_sqv'
        ExtDiag(idx)%desc = 'detrainment water vapor tendency (mynn)'
        ExtDiag(idx)%unit = 'kg kg-1 s-1'
        ExtDiag(idx)%mod_name = 'gfs_phys'
        allocate (ExtDiag(idx)%data(nblks))
        do nb = 1,nblks
          ExtDiag(idx)%data(nb)%var3 => IntDiag%det_sqv(Model%chunk_begin(nb):Model%chunk_end(nb),:)
        enddo

      endif
    endif

!  print *,'in gfdl_diag_register,af all extdiag, idx=',idx

!--- prognostic variable tendencies (t, u, v, sph, clwmr, o3)
!rab    idx = idx + 1
!rab    ExtDiag(idx)%axes = 3
!rab    ExtDiag(idx)%name = 'dtemp_dt'
!rab    ExtDiag(idx)%desc = 'gfs radiation/physics temperature tendency'
!rab    ExtDiag(idx)%unit = 'k/s'
!rab    ExtDiag(idx)%mod_name = 'gfs_phys'
!rab
!rab    idx = idx + 1
!rab    ExtDiag(idx)%axes = 3
!rab    ExtDiag(idx)%name = 'du_dt'
!rab    ExtDiag(idx)%desc = 'gfs radiation/physics horizontal wind component tendency'
!rab    ExtDiag(idx)%unit = 'm/s/s'
!rab    ExtDiag(idx)%mod_name = 'gfs_phys'
!rab
!rab    idx = idx + 1
!rab    ExtDiag(idx)%axes = 3
!rab    ExtDiag(idx)%name = 'dv_dt'
!rab    ExtDiag(idx)%desc = 'gfs radiation/physics meridional wind component tendency'
!rab    ExtDiag(idx)%unit = 'm/s/s'
!rab    ExtDiag(idx)%mod_name = 'gfs_phys'
!rab
!rab    idx = idx + 1
!rab    ExtDiag(idx)%axes = 3
!rab    ExtDiag(idx)%name = 'dsphum_dt'
!rab    ExtDiag(idx)%desc = 'gfs radiation/physics specific humidity tendency'
!rab    ExtDiag(idx)%unit = 'kg/kg/s'
!rab    ExtDiag(idx)%mod_name = 'gfs_phys'
!rab
!rab    idx = idx + 1
!rab    ExtDiag(idx)%axes = 3
!rab    ExtDiag(idx)%name = 'dclwmr_dt'
!rab    ExtDiag(idx)%desc = 'gfs radiation/radiation cloud water mixing ratio tendency'
!rab    ExtDiag(idx)%unit = 'kg/kg/s'
!rab    ExtDiag(idx)%mod_name = 'gfs_phys'
!rab
!rab    idx = idx + 1
!rab    ExtDiag(idx)%axes = 3
!rab    ExtDiag(idx)%name = 'do3mr_dt'
!rab    ExtDiag(idx)%desc = 'gfs radiation/radiation ozone mixing ratio tendency'
!rab    ExtDiag(idx)%unit = 'kg/kg/s'
!rab    ExtDiag(idx)%mod_name = 'gfs_phys'

  ! Auxiliary 2d arrays to output (for debugging)
  do num=1,Model%naux2d
    write (xtra,'(I2.2)') num
    idx = idx + 1
    ExtDiag(idx)%axes = 2
    ExtDiag(idx)%name = 'aux2d_'//trim(xtra)
    ExtDiag(idx)%desc = 'auxiliary 2d array '//trim(xtra)
    ExtDiag(idx)%unit = 'unknown'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    ExtDiag(idx)%time_avg = Model%aux2d_time_avg(num)
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var2 => IntDiag%aux2d(Model%chunk_begin(nb):Model%chunk_end(nb),num)
    enddo
  enddo

  ! Auxiliary 3d arrays to output (for debugging)
  do num=1,Model%naux3d
    write (xtra,'(I2.2)') num
    idx = idx + 1
    ExtDiag(idx)%axes = 3
    ExtDiag(idx)%name = 'aux3d_'//trim(xtra)
    ExtDiag(idx)%desc = 'auxiliary 3d array '//trim(xtra)
    ExtDiag(idx)%unit = 'unknown'
    ExtDiag(idx)%mod_name = 'gfs_phys'
    ExtDiag(idx)%intpl_method = 'bilinear'
    ExtDiag(idx)%time_avg = Model%aux3d_time_avg(num)
    allocate (ExtDiag(idx)%data(nblks))
    do nb = 1,nblks
      ExtDiag(idx)%data(nb)%var3 => IntDiag%aux3d(Model%chunk_begin(nb):Model%chunk_end(nb),:,num)
    enddo
  enddo

  end subroutine GFS_externaldiag_populate

  subroutine clm_lake_externaldiag_populate(ExtDiag, Model, Sfcprop, idx, cn_one, nblks)
    implicit none
    type(GFS_externaldiag_type),  intent(inout) :: ExtDiag(:)
    type(GFS_control_type),       intent(in)    :: Model
    type(GFS_sfcprop_type),       intent(in)    :: Sfcprop
    integer,                      intent(inout) :: idx
    integer,                      intent(in)    :: nblks
    real(kind=kind_phys),         intent(in)    :: cn_one
    character(:), allocatable :: fullname

    integer :: nk, idx0, iblk

    do iblk=1,nblks
      call link_all_levels(Sfcprop%lake_snow_z3d(Model%chunk_begin(iblk):Model%chunk_end(iblk),:), 'lake_snow_z3d', 'lake snow level depth', 'm')
    enddo

    do iblk=1,nblks
      call link_all_levels(Sfcprop%lake_snow_dz3d(Model%chunk_begin(iblk):Model%chunk_end(iblk),:), 'lake_snow_dz3d', 'lake snow level thickness', 'm')
    enddo

    do iblk=1,nblks
      call link_all_levels(Sfcprop%lake_snow_zi3d(Model%chunk_begin(iblk):Model%chunk_end(iblk),:), 'lake_snow_zi3d', 'lake snow interface depth', 'm')
    enddo

    do iblk=1,nblks
      call link_all_levels(Sfcprop%lake_h2osoi_vol3d(Model%chunk_begin(iblk):Model%chunk_end(iblk),:), 'lake_h2osoi_vol3d', 'volumetric soil water', 'm3 m-3')
    enddo

    do iblk=1,nblks
      call link_all_levels(Sfcprop%lake_h2osoi_liq3d(Model%chunk_begin(iblk):Model%chunk_end(iblk),:), 'lake_h2osoi_liq3d', 'soil liquid water content', 'kg m-2')
    enddo

    do iblk=1,nblks
      call link_all_levels(Sfcprop%lake_h2osoi_ice3d(Model%chunk_begin(iblk):Model%chunk_end(iblk),:), 'lake_h2osoi_ice3d', 'soil ice water content', 'kg m-2')
    enddo

    do iblk=1,nblks
      call link_all_levels(Sfcprop%lake_t_soisno3d(Model%chunk_begin(iblk):Model%chunk_end(iblk),:), 'lake_t_soisno3d', 'snow or soil level temperature', 'K')
    enddo

    do iblk=1,nblks
      call link_all_levels(Sfcprop%lake_t_lake3d(Model%chunk_begin(iblk):Model%chunk_end(iblk),:), 'lake_t_lake3d', 'lake layer temperature', 'K')
    enddo

    do iblk=1,nblks
      call link_all_levels(Sfcprop%lake_icefrac3d(Model%chunk_begin(iblk):Model%chunk_end(iblk),:), 'lake_icefrac3d', 'lake fractional ice cover', 'fraction')
    enddo

  contains

    subroutine link_all_levels(var3d, varname, levelname, unit)
      implicit none
      real(kind=kind_phys), target :: var3d(:,:)
      character(len=*), intent(in) :: varname, levelname, unit
      integer k, b, namelen

      if(iblk==1) then
        namelen = 30+max(len(varname),len(levelname))
        allocate(character(namelen) :: fullname)
        idx0 = idx
      endif

      var_z_loop: do k=1,size(var3d,2)
        idx = idx0 + k
        if(iblk==1) then
          ExtDiag(idx)%axes = 2
          write(fullname,"(A,'_',I0)") trim(varname),k
          ExtDiag(idx)%name = trim(fullname)
          write(fullname,"(A,' level ',I0,' of ',I0)") trim(levelname),k,size(var3d,2)
          ExtDiag(idx)%desc = trim(fullname)
          ExtDiag(idx)%unit = trim(unit)
          ExtDiag(idx)%mod_name = 'gfs_sfc'
          ExtDiag(idx)%intpl_method = 'nearest_stod'

          allocate (ExtDiag(idx)%data(nblks))
          do b=1,nblks
            nullify(ExtDiag(idx)%data(b)%var2)
          enddo
        endif

        ExtDiag(idx)%data(iblk)%var2 => var3d(:,k)
      enddo var_z_loop

      if(iblk==nblks) then
        deallocate(fullname)
      endif
    end subroutine link_all_levels
  end subroutine clm_lake_externaldiag_populate

  function soil_layer_depth(lsm, lsm_ruc, lsm_noah, layer) result(layer_depth)
     character(len=30)   :: layer_depth
     integer, intent(in) :: lsm, lsm_ruc, lsm_noah, layer
     !
     continue
     !
     if (lsm==lsm_ruc) then
        select case (layer)
           case (1)
              layer_depth = 'at 0 cm depth'
           case (2)
              layer_depth = 'at 5 cm depth'
           case (3)
              layer_depth = 'at 20 cm depth'
           case (4)
              layer_depth = 'at 40 cm depth'
           case (5)
              layer_depth = 'at 60 cm depth'
           case (6)
              layer_depth = 'at 100 cm depth'
           case (7)
              layer_depth = 'at 160 cm depth'
           case (8)
              layer_depth = 'at 220 cm depth'
           case (9)
              layer_depth = 'at 300 cm depth'
           case default
              write (layer_depth,'(a,i0)') 'invalid layer ', layer
        end select
     else if (lsm==lsm_noah) then
        select case (layer)
           case (1)
              layer_depth = '0-10cm'
           case (2)
              layer_depth = '10-40cm'
           case (3)
              layer_depth = '40-100cm'
           case (4)
              layer_depth = '100-200cm'
           case default
              write (layer_depth,'(a,i0)') 'invalid layer ', layer
        end select
     else
        write (layer_depth,'(a,i0)') 'unknown layer ', layer
     end if
     !
     return
     !
  end function soil_layer_depth

!-------------------------------------------------------------------------

end module GFS_diagnostics
