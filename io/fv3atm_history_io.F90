!> \file fv3atm_history_io.F90
!! This file defines routines used to output atmosphere diagnostic
!! (history) data from the physics and surface fields, both for quilt
!! and non-quilt output.
module fv3atm_history_io_mod

  !
  !--- FMS/GFDL modules
  use block_control_mod,  only: block_control_type
  use mpp_mod,            only: mpp_error,  mpp_pe, mpp_root_pe, FATAL
  use mpp_domains_mod,    only: domain1d, domainUG
  use time_manager_mod,   only: time_type
  use diag_manager_mod,   only: register_diag_field, send_data
  use diag_axis_mod,      only: get_axis_global_length, get_diag_axis, &
       get_diag_axis_name
  use diag_data_mod,      only: output_fields, max_output_fields
  use diag_util_mod,      only: find_input_field
  use constants_mod,      only: grav, rdgas
  !
  !--- GFS_typedefs
  use GFS_typedefs,       only: GFS_control_type, kind_phys
  use GFS_diagnostics,    only: GFS_externaldiag_type

  !
  !-----------------------------------------------------------------------
  implicit none
  private

  !--- public interfaces ---
  public  fv3atm_diag_register, fv3atm_diag_output
#ifdef use_WRTCOMP
  public  fv_phys_bundle_setup
#endif

  !>\defgroup fv3atm_history_io_mod FV3ATM History I/O Module
  !> @{

  !>@ The maximum allowed number of diagnostic fields that can be defined in any given model run.
  !! This does not include rrfs-sd or clm lake, which have their own data structures.
  integer, parameter, public :: DIAG_SIZE = 800

  real, parameter :: missing_value = 9.99e20_kind_phys
  real, parameter :: stndrd_atmos_ps = 101325.0_kind_phys
  real, parameter :: stndrd_atmos_lapse = 0.0065_kind_phys
  real, parameter :: drythresh = 1.e-4_kind_phys
  real, parameter :: zero = 0.0_kind_phys, one = 1.0_kind_phys

  !>@ Storage type for temporary data during output of diagnostic (history) files
  type history_type
    integer :: tot_diag_idx = 0

    integer :: isco=0,ieco=0,jsco=0,jeco=0,num_axes_phys=0
    integer :: ncld=0, nsoil=0, nsoil_lsm=0, imp_physics=0, landsfcmdl=0
    real(4) :: fhzero=0.
    real(4) :: dtp=0
    integer,dimension(:),        pointer         :: levo => null()
    integer,dimension(:),        pointer         :: nstt => null()
    integer,dimension(:),        pointer         :: nstt_vctbl => null()
    integer,dimension(:),        pointer         :: all_axes => null()
    character(20),dimension(:),  pointer         :: axis_name => null()
    real(4), dimension(:,:,:),   pointer         :: buffer_phys_bl => null()
    real(4), dimension(:,:,:),   pointer         :: buffer_phys_nb => null()
    real(4), dimension(:,:,:,:), pointer         :: buffer_phys_windvect => null()
    real(kind=kind_phys),dimension(:,:),pointer  :: lon => null()
    real(kind=kind_phys),dimension(:,:),pointer  :: lat => null()
    real(kind=kind_phys),dimension(:,:),pointer  :: uwork => null()
    real(kind=kind_phys),dimension(:,:,:),pointer:: uwork3d => null()
    logical                    :: uwork_set = .false.
    character(128)             :: uwindname = "(noname)"

    !--- miscellaneous other variables
    logical :: use_wrtgridcomp_output = .FALSE.
  contains
    procedure :: register => history_type_register
    procedure :: output => history_type_output
    procedure :: store_data => history_type_store_data
    procedure :: store_data3D => history_type_store_data3D
#ifdef use_WRTCOMP
    procedure :: bundle_setup => history_type_bundle_setup
    procedure :: add_field_to_phybundle => history_type_add_field_to_phybundle
    procedure :: find_output_name => history_type_find_output_name
#endif
  end type history_type

  !>@ This shared_history_data instance of history_type is shared between all calls to public module subroutines.
  type(history_type) :: shared_history_data

CONTAINS

  !>@brief Registers diagnostic variables with the FMS diagnostic manager.
  !> \section fv3atm_diag_register subroutine
  !! Creates and populates a data type which is then used to "register"
  !! diagnostic variables with the GFDL FMS diagnostic manager.
  !! includes short & long names, units, conversion factors, etc.
  !! there is no copying of data, but instead a clever use of pointers.
  !! calls a GFDL FMS routine to register diagnositcs and compare against
  !! the diag_table to determine what variables are to be output.
  subroutine fv3atm_diag_register(Diag, Time, Atm_block, Model, xlon, xlat, axes)
    use physcons,  only: con_g
    implicit none
    !--- subroutine interface variable definitions
    type(GFS_externaldiag_type),       intent(inout) :: Diag(:)
    type(time_type),           intent(in)    :: Time
    type (block_control_type), intent(in)    :: Atm_block
    type(GFS_control_type),    intent(in)    :: Model
    real(kind=kind_phys),      intent(in)    :: xlon(:,:)
    real(kind=kind_phys),      intent(in)    :: xlat(:,:)
    integer, dimension(4),     intent(in)    :: axes

    call shared_history_data%register(Diag, Time, Atm_block, Model, xlon, xlat, axes)
  end subroutine fv3atm_diag_register

  !>@brief Transfers diagnostic data to the FMS diagnostic manager
  !> \section fv3atm_diag_output subroutine
  !! This routine transfers diagnostic data to the FMS diagnostic
  !!  manager for eventual output to the history files.
  subroutine fv3atm_diag_output(time, diag, atm_block, nx, ny, levs, ntcw, ntoz, &
       dt, time_int, time_intfull, time_radsw, time_radlw)
    !--- subroutine interface variable definitions
    type(time_type),           intent(in) :: time
    type(GFS_externaldiag_type),       intent(in) :: diag(:)
    type (block_control_type), intent(in) :: atm_block
    integer,                   intent(in) :: nx, ny, levs, ntcw, ntoz
    real(kind=kind_phys),      intent(in) :: dt
    real(kind=kind_phys),      intent(in) :: time_int
    real(kind=kind_phys),      intent(in) :: time_intfull
    real(kind=kind_phys),      intent(in) :: time_radsw
    real(kind=kind_phys),      intent(in) :: time_radlw

    call shared_history_data%output(time, diag, atm_block, nx, ny, levs, ntcw, ntoz, &
         dt, time_int, time_intfull, time_radsw, time_radlw)

  end subroutine fv3atm_diag_output

#ifdef use_WRTCOMP
  !>@brief Sets up the ESMF bundle to use for quilt diagnostic output
  !> \section fv_phys_bundle_setup subroutine
  !! This part of the write component (quilt) sets up the ESMF bundles
  !! to use for writing diagnostic output. It is only defined when the
  !! write component is enabled at compile time.
  subroutine fv_phys_bundle_setup(Diag, axes, phys_bundle, fcst_grid, quilting, nbdlphys, rc)
    !
    !-------------------------------------------------------------
    !*** set esmf bundle for phys output fields
    !------------------------------------------------------------
    !
    use esmf
    use diag_data_mod, ONLY:  diag_atttype
    !
    implicit none
    !
    type(GFS_externaldiag_type),intent(in)      :: Diag(:)
    integer, intent(in)                         :: axes(:)
    type(ESMF_FieldBundle),intent(inout)        :: phys_bundle(:)
    type(ESMF_Grid),intent(inout)               :: fcst_grid
    logical,intent(in)                          :: quilting
    integer, intent(in)                         :: nbdlphys
    integer,intent(out)                         :: rc

    call shared_history_data%bundle_setup(Diag, axes, phys_bundle, fcst_grid, quilting, nbdlphys, rc)
  end subroutine fv_phys_bundle_setup
#endif

  !>@brief Private implementation of fv3atm_diag_register. Do not call directly.
  !> \section history_type%register procedure
  !! This is the history_type%register procedure, which provides the internal
  !! implementation of fv3atm_diag_register. Do not call this directly.
  subroutine history_type_register(hist, Diag, Time, Atm_block, Model, xlon, xlat, axes)
    use physcons,  only: con_g
    implicit none
    !--- subroutine interface variable definitions
    class(history_type)                      :: hist
    type(GFS_externaldiag_type),       intent(inout) :: Diag(:)
    type(time_type),           intent(in)    :: Time
    type (block_control_type), intent(in)    :: Atm_block
    type(GFS_control_type),    intent(in)    :: Model
    real(kind=kind_phys),      intent(in)    :: xlon(:,:)
    real(kind=kind_phys),      intent(in)    :: xlat(:,:)
    integer, dimension(4),     intent(in)    :: axes
    !--- local variables
    integer :: idx, nrgst_bl, nrgst_nb, nrgst_vctbl

    hist%isco   = Atm_block%isc
    hist%ieco   = Atm_block%iec
    hist%jsco   = Atm_block%jsc
    hist%jeco   = Atm_block%jec
    hist%fhzero = Model%fhzero
    !   hist%ncld   = Model%ncld
    hist%ncld   = Model%imp_physics
    hist%nsoil  = Model%lsoil
    hist%nsoil_lsm  = Model%lsoil_lsm
    hist%dtp    = Model%dtp
    hist%imp_physics  = Model%imp_physics
    hist%landsfcmdl  = Model%lsm
    !    print *,'in fv3atm_diag_register,hist%ncld=',Model%ncld,Model%lsoil,Model%imp_physics, &
    !      ' hist%dtp=',hist%dtp,' hist%landsfcmdl=',Model%lsm
    !
    !save lon/lat for vector interpolation
    allocate(hist%lon(hist%isco:hist%ieco,hist%jsco:hist%jeco))
    allocate(hist%lat(hist%isco:hist%ieco,hist%jsco:hist%jeco))
    hist%lon = xlon
    hist%lat = xlat

    do idx = 1,DIAG_SIZE
      if (trim(Diag(idx)%name) == '') exit
      hist%tot_diag_idx = idx
    enddo

    if (hist%tot_diag_idx == DIAG_SIZE) then
      call mpp_error(fatal, 'fv3atm_io::fv3atm_diag_register - need to increase parameter DIAG_SIZE')
    endif

    allocate(hist%levo(hist%tot_diag_idx))
    allocate(hist%nstt(hist%tot_diag_idx), hist%nstt_vctbl(hist%tot_diag_idx))
    hist%levo          = 0
    hist%nstt          = 0
    hist%nstt_vctbl    = 0
    nrgst_bl      = 0
    nrgst_nb      = 0
    nrgst_vctbl   = 0
    hist%num_axes_phys = 2
    do idx = 1,hist%tot_diag_idx
      if (diag(idx)%axes == -99) then
        call mpp_error(fatal, 'gfs_driver::gfs_diag_register - attempt to register an undefined variable')
      endif
      Diag(idx)%id = register_diag_field (trim(Diag(idx)%mod_name), trim(Diag(idx)%name),  &
           axes(1:Diag(idx)%axes), Time, trim(Diag(idx)%desc), &
           trim(Diag(idx)%unit), missing_value=real(missing_value))
      if(Diag(idx)%id > 0) then
        if (Diag(idx)%axes == 2) then
          hist%levo(idx) = 1
          if( index(trim(Diag(idx)%intpl_method),'bilinear') > 0 ) then
            nrgst_bl = nrgst_bl + 1
            hist%nstt(idx) = nrgst_bl
          else if (trim(Diag(idx)%intpl_method) == 'nearest_stod' ) then
            nrgst_nb = nrgst_nb + 1
            hist%nstt(idx) = nrgst_nb
          endif
          if(trim(Diag(idx)%intpl_method) == 'vector_bilinear') then
            if(Diag(idx)%name(1:1) == 'v' .or. Diag(idx)%name(1:1) == 'V') then
              nrgst_vctbl = nrgst_vctbl + 1
              hist%nstt_vctbl(idx) = nrgst_vctbl
              !             print *,'in phy_setup, vector_bilinear, name=', trim(Diag(idx)%name),' nstt_vctbl=', hist%nstt_vctbl(idx), 'idx=',idx
            endif
          endif
        else if (diag(idx)%axes == 3) then
          hist%levo(idx) = size(Diag(idx)%data(1)%var3, dim=2)
          if( index(trim(diag(idx)%intpl_method),'bilinear') > 0 ) then
            hist%nstt(idx) = nrgst_bl + 1
            nrgst_bl  = nrgst_bl + hist%levo(idx)
          else if (trim(diag(idx)%intpl_method) == 'nearest_stod' ) then
            hist%nstt(idx) = nrgst_nb + 1
            nrgst_nb  = nrgst_nb + hist%levo(idx)
          endif
          if(trim(diag(idx)%intpl_method) == 'vector_bilinear') then
            if(diag(idx)%name(1:1) == 'v' .or. diag(idx)%name(1:1) == 'V') then
              hist%nstt_vctbl(idx) = nrgst_vctbl + 1
              nrgst_vctbl = nrgst_vctbl + hist%levo(idx)
              !             print *,'in phy_setup, vector_bilinear, name=', trim(diag(idx)%name),' nstt_vctbl=', hist%nstt_vctbl(idx), 'idx=',idx
            endif
          endif
          hist%num_axes_phys = 3
        endif
      endif

    enddo

    allocate(hist%buffer_phys_bl(hist%isco:hist%ieco,hist%jsco:hist%jeco,nrgst_bl))
    allocate(hist%buffer_phys_nb(hist%isco:hist%ieco,hist%jsco:hist%jeco,nrgst_nb))
    allocate(hist%buffer_phys_windvect(3,hist%isco:hist%ieco,hist%jsco:hist%jeco,nrgst_vctbl))
    hist%buffer_phys_bl = zero
    hist%buffer_phys_nb = zero
    hist%buffer_phys_windvect = zero
    if(mpp_pe() == mpp_root_pe()) print *,'in fv3atm_diag_register, nrgst_bl=',nrgst_bl,' nrgst_nb=',nrgst_nb, &
         ' nrgst_vctbl=',nrgst_vctbl, 'hist%isco=',hist%isco,hist%ieco,'hist%jsco=',hist%jsco,hist%jeco,' hist%num_axes_phys=', hist%num_axes_phys

  end subroutine history_type_register

  !>@brief Internal implementation of fv3atm_diag_output
  !> \section history_type%output procedure
  !! This is history_type%output, which provides the internal
  !! implementation of the public fv3atm_diag_output routine. Never
  !! call this directly.
  subroutine history_type_output(hist, time, diag, atm_block, nx, ny, levs, ntcw, ntoz, &
       dt, time_int, time_intfull, time_radsw, time_radlw)
    !--- subroutine interface variable definitions
    class(history_type)                   :: hist
    type(time_type),           intent(in) :: time
    type(GFS_externaldiag_type),       intent(in) :: diag(:)
    type (block_control_type), intent(in) :: atm_block
    integer,                   intent(in) :: nx, ny, levs, ntcw, ntoz
    real(kind=kind_phys),      intent(in) :: dt
    real(kind=kind_phys),      intent(in) :: time_int
    real(kind=kind_phys),      intent(in) :: time_intfull
    real(kind=kind_phys),      intent(in) :: time_radsw
    real(kind=kind_phys),      intent(in) :: time_radlw
    !--- local variables
    integer :: i, j, k, idx, nb, ix, ii, jj, levo_3d
    character(len=2) :: xtra
#ifdef CCPP_32BIT
    real, dimension(nx,ny)      :: var2
    real, dimension(:,:,:), allocatable :: var3
#else
    real(kind=kind_phys), dimension(nx,ny)      :: var2
    real(kind=kind_phys), dimension(:,:,:), allocatable :: var3
#endif
    real(kind=kind_phys) :: rtime_int, rtime_intfull, lcnvfac
    real(kind=kind_phys) :: rtime_radsw, rtime_radlw

    rtime_int     = one/time_int
    rtime_intfull = one/time_intfull
    rtime_radsw   = one/time_radsw
    rtime_radlw   = one/time_radlw

    !     if(mpp_pe()==mpp_root_pe())print *,'in,fv3atm_io. time avg, time_int=',time_int
    history_loop: do idx = 1,hist%tot_diag_idx
      has_id: if (diag(idx)%id > 0) then
        lcnvfac = diag(idx)%cnvfac
        if (diag(idx)%time_avg) then
          if ( trim(diag(idx)%time_avg_kind) == 'full' ) then
            lcnvfac = lcnvfac*rtime_intfull
            !             if(mpp_pe()==mpp_root_pe())print *,'in,fv3atm_io. full time avg, field=',trim(Diag(idx)%name),' time=',time_intfull
          else if ( trim(diag(idx)%time_avg_kind) == 'rad_lw' ) then
            lcnvfac = lcnvfac*min(rtime_radlw,rtime_int)
            !             if(mpp_pe()==mpp_root_pe())print *,'in,fv3atm_io. rad longwave avg, field=',trim(Diag(idx)%name),' time=',time_radlw
          else if ( trim(diag(idx)%time_avg_kind) == 'rad_sw' ) then
            lcnvfac = lcnvfac*min(rtime_radsw,rtime_int)
            !             if(mpp_pe()==mpp_root_pe())print *,'in,fv3atm_io. rad shortwave avg, field=',trim(Diag(idx)%name),' time=',time_radsw
          else if ( trim(diag(idx)%time_avg_kind) == 'rad_swlw_min' ) then
            lcnvfac = lcnvfac*min(max(rtime_radsw,rtime_radlw),rtime_int)
            !             if(mpp_pe()==mpp_root_pe())print *,'in,fv3atm_io. rad swlw min avg, field=',trim(Diag(idx)%name),' time=',time_radlw,time_radsw,time_int
          else
            lcnvfac = lcnvfac*rtime_int
          endif
        endif
        if_2d: if (diag(idx)%axes == 2) then
          ! Integer data
          int_or_real: if (associated(Diag(idx)%data(1)%int2)) then
            if (trim(Diag(idx)%intpl_method) == 'nearest_stod') then
              var2(1:nx,1:ny) = 0._kind_phys
              do j = 1, ny
                jj = j + Atm_block%jsc -1
                do i = 1, nx
                  ii = i + Atm_block%isc -1
                  nb = Atm_block%blkno(ii,jj)
                  ix = Atm_block%ixp(ii,jj)
                  var2(i,j) = real(Diag(idx)%data(nb)%int2(ix), kind=kind_phys)
                enddo
              enddo
              call hist%store_data(Diag(idx)%id, var2, Time, idx, Diag(idx)%intpl_method, Diag(idx)%name)
            else
              call mpp_error(FATAL, 'Interpolation method ' // trim(Diag(idx)%intpl_method) // ' for integer variable ' &
                   // trim(Diag(idx)%name) // ' not supported.')
            endif
            ! Real data
          else ! int_or_real
            if_mask: if (trim(diag(idx)%mask) == 'positive_flux') then
              !--- albedos are actually a ratio of two radiation surface properties
              var2(1:nx,1:ny) = 0._kind_phys
              do j = 1, ny
                jj = j + Atm_block%jsc -1
                do i = 1, nx
                  ii = i + Atm_block%isc -1
                  nb = Atm_block%blkno(ii,jj)
                  ix = Atm_block%ixp(ii,jj)
                  if (Diag(idx)%data(nb)%var21(ix) > 0._kind_phys) &
                       var2(i,j) = max(0._kind_phys,min(1._kind_phys,Diag(idx)%data(nb)%var2(ix)/Diag(idx)%data(nb)%var21(ix)))*lcnvfac
                enddo
              enddo
            elseif (trim(Diag(idx)%mask) == 'land_ice_only') then
              !--- need to "mask" gflux to output valid data over land/ice only
              var2(1:nx,1:ny) = missing_value
              do j = 1, ny
                jj = j + Atm_block%jsc -1
                do i = 1, nx
                  ii = i + Atm_block%isc -1
                  nb = Atm_block%blkno(ii,jj)
                  ix = Atm_block%ixp(ii,jj)
                  if (Diag(idx)%data(nb)%var21(ix) /= 0) var2(i,j) = Diag(idx)%data(nb)%var2(ix)*lcnvfac
                enddo
              enddo
            elseif (trim(Diag(idx)%mask) == 'land_only') then
              !--- need to "mask" soilm to have value only over land
              var2(1:nx,1:ny) = missing_value
              do j = 1, ny
                jj = j + Atm_block%jsc -1
                do i = 1, nx
                  ii = i + Atm_block%isc -1
                  nb = Atm_block%blkno(ii,jj)
                  ix = Atm_block%ixp(ii,jj)
                  if (Diag(idx)%data(nb)%var21(ix) == 1) var2(i,j) = Diag(idx)%data(nb)%var2(ix)*lcnvfac
                enddo
              enddo
            elseif (trim(Diag(idx)%mask) == 'cldmask') then
              !--- need to "mask" soilm to have value only over land
              var2(1:nx,1:ny) = missing_value
              do j = 1, ny
                jj = j + Atm_block%jsc -1
                do i = 1, nx
                  ii = i + Atm_block%isc -1
                  nb = Atm_block%blkno(ii,jj)
                  ix = Atm_block%ixp(ii,jj)
                  if (Diag(idx)%data(nb)%var21(ix)*100. > 0.5) var2(i,j) = Diag(idx)%data(nb)%var2(ix)*lcnvfac
                enddo
              enddo
            elseif (trim(Diag(idx)%mask) == 'cldmask_ratio') then
              !--- need to "mask" soilm to have value only over land
              var2(1:nx,1:ny) = missing_value
              do j = 1, ny
                jj = j + Atm_block%jsc -1
                do i = 1, nx
                  ii = i + Atm_block%isc -1
                  nb = Atm_block%blkno(ii,jj)
                  ix = Atm_block%ixp(ii,jj)
                  if (Diag(idx)%data(nb)%var21(ix)*100.*lcnvfac > 0.5) var2(i,j) = Diag(idx)%data(nb)%var2(ix)/ &
                       Diag(idx)%data(nb)%var21(ix)
                enddo
              enddo
            elseif (trim(Diag(idx)%mask) == 'pseudo_ps') then
              if ( hist%use_wrtgridcomp_output ) then
                do j = 1, ny
                  jj = j + Atm_block%jsc -1
                  do i = 1, nx
                    ii = i + Atm_block%isc -1
                    nb = Atm_block%blkno(ii,jj)
                    ix = Atm_block%ixp(ii,jj)
                    var2(i,j) = (Diag(idx)%data(nb)%var2(ix)/stndrd_atmos_ps)**(rdgas/grav*stndrd_atmos_lapse)
                  enddo
                enddo
              else
                do j = 1, ny
                  jj = j + Atm_block%jsc -1
                  do i = 1, nx
                    ii = i + Atm_block%isc -1
                    nb = Atm_block%blkno(ii,jj)
                    ix = Atm_block%ixp(ii,jj)
                    var2(i,j) = Diag(idx)%data(nb)%var2(ix)
                  enddo
                enddo
              endif
            elseif (trim(Diag(idx)%mask) == '') then
              do j = 1, ny
                jj = j + Atm_block%jsc -1
                do i = 1, nx
                  ii = i + Atm_block%isc -1
                  nb = Atm_block%blkno(ii,jj)
                  ix = Atm_block%ixp(ii,jj)
                  var2(i,j) = Diag(idx)%data(nb)%var2(ix)*lcnvfac
                enddo
              enddo
            endif if_mask
          endif int_or_real

          call hist%store_data(Diag(idx)%id, var2, Time, idx, Diag(idx)%intpl_method, Diag(idx)%name)

        elseif (Diag(idx)%axes == 3) then
          !---
          !--- skipping other 3D variables with the following else statement
          !---

          levo_3d = hist%levo(idx)
          allocate(var3(nx,ny,levo_3d))

          do k=1, levo_3d
            do j = 1, ny
              jj = j + Atm_block%jsc -1
              do i = 1, nx
                ii = i + Atm_block%isc -1
                nb = Atm_block%blkno(ii,jj)
                ix = Atm_block%ixp(ii,jj)
                ! flip only 3d variables with vertical dimension == levs (atm model levels)
                if (levo_3d == levs) then
                  var3(i,j,k) = Diag(idx)%data(nb)%var3(ix,levo_3d-k+1)*lcnvfac
                else
                  var3(i,j,k) = Diag(idx)%data(nb)%var3(ix,        k  )*lcnvfac
                endif
              enddo
            enddo
          enddo

          call hist%store_data3D(Diag(idx)%id, var3, Time, idx, Diag(idx)%intpl_method, Diag(idx)%name)
          deallocate(var3)

        endif if_2d
      endif has_id
    end do history_loop
  end subroutine history_type_output

  !>@brief Part of the internal implementation of history_type_output (history_type%output)
  !> \section history_type%store_data procedure
  !! This routine copies data from an x-y array to internal buffers for later output.
  !! Never call this subroutine directly; call fv3atm_diag_output instead.
  subroutine history_type_store_data(hist,id, work, Time, idx, intpl_method, fldname)
    implicit none
    class(history_type)                 :: hist
    integer, intent(in)                 :: id
    integer, intent(in)                 :: idx
#ifdef CCPP_32BIT
    real, intent(in)                    :: work(:,:)
#else
    real(kind=kind_phys), intent(in)    :: work(hist%ieco-hist%isco+1,hist%jeco-hist%jsco+1)
#endif
    type(time_type), intent(in)         :: Time
    character(*), intent(in)            :: intpl_method
    character(*), intent(in)            :: fldname
    !
    real(kind=kind_phys)                :: sinlat, sinlon, coslon
    integer j,i,nv,i1,j1
    logical used
    !
    if_has_id: if( id > 0 ) then
      if_gridcomp: if( hist%use_wrtgridcomp_output ) then
        if_interp: if( trim(intpl_method) == 'bilinear') then
          !$omp parallel do default(shared) private(i,j)
          do j= hist%jsco,hist%jeco
            do i= hist%isco,hist%ieco
              hist%buffer_phys_bl(i,j,hist%nstt(idx)) = work(i-hist%isco+1,j-hist%jsco+1)
            enddo
          enddo
        else if(trim(intpl_method) == 'nearest_stod') then
          !$omp parallel do default(shared) private(i,j)
          do j= hist%jsco,hist%jeco
            do i= hist%isco,hist%ieco
              hist%buffer_phys_nb(i,j,hist%nstt(idx)) = work(i-hist%isco+1,j-hist%jsco+1)
            enddo
          enddo
        else if(trim(intpl_method) == 'vector_bilinear') then
          !first save the data
          !$omp parallel do default(shared) private(i,j)
          do j= hist%jsco,hist%jeco
            do i= hist%isco,hist%ieco
              hist%buffer_phys_bl(i,j,hist%nstt(idx)) = work(i-hist%isco+1,j-hist%jsco+1)
            enddo
          enddo
          if_u_wind: if( fldname(1:1) == 'u' .or. fldname(1:1) == 'U') then
            if(.not.associated(hist%uwork)) allocate(hist%uwork(hist%isco:hist%ieco,hist%jsco:hist%jeco))
            !$omp parallel do default(shared) private(i,j)
            do j= hist%jsco,hist%jeco
              do i= hist%isco,hist%ieco
                hist%uwork(i,j) = work(i-hist%isco+1,j-hist%jsco+1)
              enddo
            enddo
            hist%uwindname = fldname
            hist%uwork_set = .true.
          endif if_u_wind
          if_v_wind: if( fldname(1:1) == 'v' .or. fldname(1:1) == 'V') then
            !set up wind vector
            if( hist%uwork_set .and. trim(hist%uwindname(2:)) == trim(fldname(2:))) then
              nv = hist%nstt_vctbl(idx)
              !$omp parallel do default(shared) private(i,j,i1,j1,sinlat,sinlon,coslon)
              do j= hist%jsco,hist%jeco
                j1 = j-hist%jsco+1
                do i= hist%isco,hist%ieco
                  i1 = i-hist%isco+1
                  sinlat = sin(hist%lat(i,j))
                  sinlon = sin(hist%lon(i,j))
                  coslon = cos(hist%lon(i,j))
                  hist%buffer_phys_windvect(1,i,j,nv) = hist%uwork(i,j)*coslon - work(i1,j1)*sinlat*sinlon
                  hist%buffer_phys_windvect(2,i,j,nv) = hist%uwork(i,j)*sinlon + work(i1,j1)*sinlat*coslon
                  hist%buffer_phys_windvect(3,i,j,nv) =                     work(i1,j1)*cos(hist%lat(i,j))
                enddo
              enddo
            endif
            hist%uwork     = zero
            hist%uwindname = ''
            hist%uwork_set = .false.
          endif if_v_wind

        endif if_interp
      else
        used = send_data(id, work, Time)
      endif if_gridcomp
    endif if_has_id
    !
  end subroutine history_type_store_data

  !>@brief Part of the internal implementation of history_type_output (history_type%output)
  !> \section history_type%store_data3D procedure
  !! This routine copies data from an x-y-z array to internal buffers for later output.
  !! Never call this subroutine directly; call fv3atm_diag_output instead.
  subroutine history_type_store_data3D(hist, id, work, Time, idx, intpl_method, fldname)
    implicit none
    class(history_type)                 :: hist
    integer, intent(in)                 :: id
    integer, intent(in)                 :: idx
#ifdef CCPP_32BIT
    real, intent(in)                    :: work(:,:,:)
#else
    real(kind=kind_phys), intent(in)    :: work(hist%ieco-hist%isco+1,hist%jeco-hist%jsco+1,hist%levo(idx))
#endif
    type(time_type), intent(in)         :: Time
    character(*), intent(in)            :: intpl_method
    character(*), intent(in)            :: fldname
    !
    real(kind=kind_phys), allocatable, dimension(:,:) :: sinlon, coslon, sinlat, coslat
    integer k,j,i,nv,i1,j1
    logical used
    !
    if( id > 0 ) then
      if( hist%use_wrtgridcomp_output ) then
        if( trim(intpl_method) == 'bilinear') then
          !$omp parallel do default(shared) private(i,j,k)
          do k= 1,hist%levo(idx)
            do j= hist%jsco,hist%jeco
              do i= hist%isco,hist%ieco
                hist%buffer_phys_bl(i,j,hist%nstt(idx)+k-1) = work(i-hist%isco+1,j-hist%jsco+1,k)
              enddo
            enddo
          enddo
        else if(trim(intpl_method) == 'nearest_stod') then
          !$omp parallel do default(shared) private(i,j,k)
          do k= 1,hist%levo(idx)
            do j= hist%jsco,hist%jeco
              do i= hist%isco,hist%ieco
                hist%buffer_phys_nb(i,j,hist%nstt(idx)+k-1) = work(i-hist%isco+1,j-hist%jsco+1,k)
              enddo
            enddo
          enddo
        else if(trim(intpl_method) == 'vector_bilinear') then
          !first save the data
          !$omp parallel do default(shared) private(i,j,k)
          do k= 1,hist%levo(idx)
            do j= hist%jsco,hist%jeco
              do i= hist%isco,hist%ieco
                hist%buffer_phys_bl(i,j,hist%nstt(idx)+k-1) = work(i-hist%isco+1,j-hist%jsco+1,k)
              enddo
            enddo
          enddo
          if( fldname(1:1) == 'u' .or. fldname(1:1) == 'U') then
            if(.not.associated(hist%uwork3d)) allocate(hist%uwork3d(hist%isco:hist%ieco,hist%jsco:hist%jeco,hist%levo(idx)))
            !$omp parallel do default(shared) private(i,j,k)
            do k= 1, hist%levo(idx)
              do j= hist%jsco,hist%jeco
                do i= hist%isco,hist%ieco
                  hist%uwork3d(i,j,k) = work(i-hist%isco+1,j-hist%jsco+1,k)
                enddo
              enddo
            enddo
            hist%uwindname = fldname
            hist%uwork_set = .true.
          endif
          if( fldname(1:1) == 'v' .or. fldname(1:1) == 'V') then
            !set up wind vector
            if( hist%uwork_set .and. trim(hist%uwindname(2:)) == trim(fldname(2:))) then
              allocate (sinlon(hist%isco:hist%ieco,hist%jsco:hist%jeco), coslon(hist%isco:hist%ieco,hist%jsco:hist%jeco), &
                   sinlat(hist%isco:hist%ieco,hist%jsco:hist%jeco), coslat(hist%isco:hist%ieco,hist%jsco:hist%jeco))
              !$omp parallel do default(shared) private(i,j)
              do j= hist%jsco,hist%jeco
                do i= hist%isco,hist%ieco
                  sinlon(i,j) = sin(hist%lon(i,j))
                  coslon(i,j) = cos(hist%lon(i,j))
                  sinlat(i,j) = sin(hist%lat(i,j))
                  coslat(i,j) = cos(hist%lat(i,j))
                enddo
              enddo
              !$omp parallel do default(shared) private(i,j,k,nv,i1,j1)
              do k= 1, hist%levo(idx)
                nv = hist%nstt_vctbl(idx)+k-1
                do j= hist%jsco,hist%jeco
                  j1 = j-hist%jsco+1
                  do i= hist%isco,hist%ieco
                    i1 = i-hist%isco+1
                    hist%buffer_phys_windvect(1,i,j,nv) = hist%uwork3d(i,j,k)*coslon(i,j) &
                         - work(i1,j1,k)*sinlat(i,j)*sinlon(i,j)
                    hist%buffer_phys_windvect(2,i,j,nv) = hist%uwork3d(i,j,k)*sinlon(i,j) &
                         + work(i1,j1,k)*sinlat(i,j)*coslon(i,j)
                    hist%buffer_phys_windvect(3,i,j,nv) = work(i1,j1,k)*coslat(i,j)
                  enddo
                enddo
              enddo
              deallocate (sinlon, coslon, sinlat, coslat)
            endif
            hist%uwork3d   = zero
            hist%uwindname = ''
            hist%uwork_set = .false.
          endif

        endif
      else
        used = send_data(id, work, Time)
      endif
    endif
    !
  end subroutine history_type_store_data3D

#ifdef use_WRTCOMP
  !>@brief Sets up the ESMF bundle to use for quilt diagnostic output
  !> \section history_type%bundle_setup procedure
  !! This part of the write component (quilt) sets up the ESMF bundles
  !! to use for writing diagnostic output. It is only defined when the
  !! write component is enabled at compile time.

  subroutine history_type_bundle_setup(hist, Diag, axes, phys_bundle, fcst_grid, quilting, nbdlphys, rc)
    ! set esmf bundle for phys output fields
    use esmf
    use diag_data_mod, ONLY:  diag_atttype
    !
    implicit none
    !
    class(history_type)                         :: hist
    type(GFS_externaldiag_type),intent(in)      :: Diag(:)
    integer, intent(in)                         :: axes(:)
    type(ESMF_FieldBundle),intent(inout)        :: phys_bundle(:)
    type(ESMF_Grid),intent(inout)               :: fcst_grid
    logical,intent(in)                          :: quilting
    integer, intent(in)                         :: nbdlphys
    integer,intent(out)                         :: rc

    !
    !*** local variables
    integer i, idx, ibdl
    integer id, axis_length, direction, edges
    integer num_attributes
    character(255)    :: units, long_name, cart_name, axis_direct, edgesS
    character(128)    :: output_name, physbdl_name, outputfile1
    logical           :: lput2physbdl, loutputfile, l2dvector
    type(domain1d)    :: Domain
    type(domainUG)    :: DomainU
    real,dimension(:),allocatable               :: axis_data
    character(128),dimension(:), allocatable    :: bdl_intplmethod, outputfile
    type(diag_atttype),dimension(:),allocatable :: attributes
    !
    logical isPresent
    integer udimCount
    character(80),dimension(:),allocatable :: udimList
    character(20),dimension(:),  allocatable         :: axis_name_vert
    !
    !------------------------------------------------------------
    !--- use wrte grid component for output
    hist%use_wrtgridcomp_output = quilting
    !   if(mpp_pe()==mpp_root_pe())print *,'in fv_phys bundle,use_wrtgridcomp_output=',hist%use_wrtgridcomp_output, &
    !   print *,'in fv_phys bundle,use_wrtgridcomp_output=',hist%use_wrtgridcomp_output, &
    !       'hist%isco=',hist%isco,hist%ieco,'hist%jsco=',hist%jsco,hist%jeco,'hist%tot_diag_idx=',hist%tot_diag_idx
    !
    !------------------------------------------------------------
    !*** add attributes to the bundle such as subdomain limtis,
    !*** axes, output time, etc
    !------------------------------------------------------------
    !
    allocate(bdl_intplmethod(nbdlphys), outputfile(nbdlphys))
    if(mpp_pe()==mpp_root_pe()) print *,'in fv_phys bundle,nbdl=',nbdlphys
    do ibdl = 1, nbdlphys
      loutputfile = .false.
      call ESMF_FieldBundleGet(phys_bundle(ibdl), name=physbdl_name,rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      idx = index(physbdl_name,'_bilinear')
      if(idx > 0) then
        outputfile(ibdl)      = physbdl_name(1:idx-1)
        bdl_intplmethod(ibdl) = 'bilinear'
        loutputfile           = .true.
      endif
      idx = index(physbdl_name,'_nearest_stod')
      if(idx > 0) then
        outputfile(ibdl)      = physbdl_name(1:idx-1)
        bdl_intplmethod(ibdl) = 'nearest_stod'
        loutputfile           = .true.
      endif
      if( .not. loutputfile) then
        outputfile(ibdl)      = 'phy'
        bdl_intplmethod(ibdl) = 'nearest_stod'
      endif
      !    print *,'in fv_phys bundle,i=',ibdl,'outputfile=',trim(outputfile(ibdl)), &
      !      'bdl_intplmethod=',trim(bdl_intplmethod(ibdl))

      call ESMF_AttributeAdd(phys_bundle(ibdl), convention="NetCDF", purpose="FV3", &
           attrList=(/"fhzero     ", "ncld       ", "nsoil      ",&
           "imp_physics", "dtp        ", "landsfcmdl "/), rc=rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call ESMF_AttributeSet(phys_bundle(ibdl), convention="NetCDF", purpose="FV3", &
           name="fhzero", value=hist%fhzero, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call ESMF_AttributeSet(phys_bundle(ibdl), convention="NetCDF", purpose="FV3", &
           name="ncld", value=hist%ncld, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call ESMF_AttributeSet(phys_bundle(ibdl), convention="NetCDF", purpose="FV3", &
           name="nsoil", value=hist%nsoil, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call ESMF_AttributeSet(phys_bundle(ibdl), convention="NetCDF", purpose="FV3", &
           name="imp_physics", value=hist%imp_physics, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call ESMF_AttributeSet(phys_bundle(ibdl), convention="NetCDF", purpose="FV3", &
           name="dtp", value=hist%dtp, rc=rc)
      !     print *,'in fcst gfdl diag, hist%dtp=',hist%dtp,' ibdl=',ibdl
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call ESMF_AttributeSet(phys_bundle(ibdl), convention="NetCDF", purpose="FV3", &
           name="landsfcmdl", value=hist%landsfcmdl, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      !end ibdl
    enddo
    !
    !*** get axis names
    allocate(hist%axis_name(hist%num_axes_phys))
    do id = 1,hist%num_axes_phys
      call get_diag_axis_name( axes(id), hist%axis_name(id))
    enddo
    isPresent = .false.
    if( hist%num_axes_phys>2 ) then
      allocate(axis_name_vert(hist%num_axes_phys-2))
      do id=3,hist%num_axes_phys
        axis_name_vert(id-2) = hist%axis_name(id)
      enddo
      !
      call ESMF_AttributeGet(fcst_grid, convention="NetCDF", purpose="FV3", &
           name="vertical_dim_labels", isPresent=isPresent, &
           itemCount=udimCount, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      if (isPresent .and. (udimCount>hist%num_axes_phys-2) ) then
        allocate(udimList(udimCount))
        call ESMF_AttributeGet(fcst_grid, convention="NetCDF", purpose="FV3", &
             name="vertical_dim_labels", valueList=udimList, rc=rc)
        !       if(mpp_pe()==mpp_root_pe()) print *,'in fv3atmio, vertical
        !       list=',udimList(1:udimCount),'rc=',rc

        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      else

        if(mpp_pe()==mpp_root_pe()) print *,'in fv_dyn bundle,axis_name_vert=',axis_name_vert
        call ESMF_AttributeAdd(fcst_grid, convention="NetCDF", purpose="FV3",  &
             attrList=(/"vertical_dim_labels"/), rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
        call ESMF_AttributeSet(fcst_grid, convention="NetCDF", purpose="FV3", &
             name="vertical_dim_labels", valueList=axis_name_vert, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      endif
      deallocate(axis_name_vert)
    endif

    !*** add attributes
    if(associated(hist%all_axes)) then
      deallocate(hist%all_axes)
      nullify(hist%all_axes)
    endif
    allocate(hist%all_axes(hist%num_axes_phys))
    hist%all_axes(1:hist%num_axes_phys) = axes(1:hist%num_axes_phys)
    if (.not. isPresent .or. (udimCount<hist%num_axes_phys-2) ) then
      do id = 1,hist%num_axes_phys
        axis_length =  get_axis_global_length(axes(id))
        allocate(axis_data(axis_length))
        call get_diag_axis( axes(id), hist%axis_name(id), units, long_name, cart_name, &
             direction, edges, Domain, DomainU, axis_data,           &
             num_attributes=num_attributes, attributes=attributes)
        !
        edgesS = ''
        do i = 1,hist%num_axes_phys
          if(axes(i) == edges) edgesS=hist%axis_name(i)
        enddo
        ! Add vertical dimension Attributes to Grid
        if( id>2 ) then
          !      if(mpp_pe()==mpp_root_pe()) print *,' in dyn add grid, axis_name=',     &
          !         trim(hist%axis_name(id)),'axis_data=',axis_data
          if(trim(edgesS)/='') then
            call ESMF_AttributeAdd(fcst_grid, convention="NetCDF", purpose="FV3",  &
                 attrList=(/trim(hist%axis_name(id)),trim(hist%axis_name(id))//":long_name",    &
                 trim(hist%axis_name(id))//":units", trim(hist%axis_name(id))//":cartesian_axis", &
                 trim(hist%axis_name(id))//":positive", trim(hist%axis_name(id))//":edges"/), rc=rc)
          else
            call ESMF_AttributeAdd(fcst_grid, convention="NetCDF", purpose="FV3",  &
                 attrList=(/trim(hist%axis_name(id)),trim(hist%axis_name(id))//":long_name",    &
                 trim(hist%axis_name(id))//":units", trim(hist%axis_name(id))//":cartesian_axis", &
                 trim(hist%axis_name(id))//":positive"/), rc=rc)
          endif
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          call ESMF_AttributeSet(fcst_grid, convention="NetCDF", purpose="FV3", &
               name=trim(hist%axis_name(id)), valueList=axis_data, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          call ESMF_AttributeSet(fcst_grid, convention="NetCDF", purpose="FV3", &
               name=trim(hist%axis_name(id))//":long_name", value=trim(long_name), rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          call ESMF_AttributeSet(fcst_grid, convention="NetCDF", purpose="FV3", &
               name=trim(hist%axis_name(id))//":units", value=trim(units), rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          call ESMF_AttributeSet(fcst_grid, convention="NetCDF", purpose="FV3", &
               name=trim(hist%axis_name(id))//":cartesian_axis", value=trim(cart_name), rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          if(direction > 0) then
            axis_direct = "up"
          else
            axis_direct = "down"
          endif
          call ESMF_AttributeSet(fcst_grid, convention="NetCDF", purpose="FV3", &
               name=trim(hist%axis_name(id))//":positive", value=trim(axis_direct), rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          if(trim(edgesS)/='') then
            call ESMF_AttributeSet(fcst_grid, convention="NetCDF", purpose="FV3", &
                 name=trim(hist%axis_name(id))//":edges", value=trim(edgesS), rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
          endif

        endif
        !
        deallocate(axis_data)
      enddo
    endif

    ! add zsoil axis
    call ESMF_AttributeAdd(fcst_grid, convention="NetCDF", purpose="FV3",  &
         attrList=(/"zsoil               ", &
                    "zsoil:long_name     ", &
                    "zsoil:units         ", &
                    "zsoil:cartesian_axis", &
                    "zsoil:positive      "/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_AttributeSet(fcst_grid, convention="NetCDF", purpose="FV3", &
         name="zsoil", valueList=(/ (i, i=1,hist%nsoil_lsm) /), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_AttributeSet(fcst_grid, convention="NetCDF", purpose="FV3", &
         name="zsoil:long_name", value="soil level", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_AttributeSet(fcst_grid, convention="NetCDF", purpose="FV3", &
         name="zsoil:units", value="level", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_AttributeSet(fcst_grid, convention="NetCDF", purpose="FV3", &
         name="zsoil:cartesian_axis", value="Z", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_AttributeSet(fcst_grid, convention="NetCDF", purpose="FV3", &
         name="zsoil:positive", value="down", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    !   print *,'in setup fieldbundle_phys, hist%num_axes_phys=',hist%num_axes_phys,'hist%tot_diag_idx=',hist%tot_diag_idx, &
    !       'nbdlphys=',nbdlphys
    !
    !-----------------------------------------------------------------------------------------
    !*** add esmf fields
    !
    do idx= 1,hist%tot_diag_idx

      lput2physbdl = .false.
      do ibdl = 1, nbdlphys

        if( index(trim(Diag(idx)%intpl_method),trim(bdl_intplmethod(ibdl))) > 0) then
          lput2physbdl = .true.
          if( Diag(idx)%id > 0 ) then
            call hist%find_output_name(trim(Diag(idx)%mod_name),trim(Diag(idx)%name),output_name)

            !add origin field
            call hist%add_field_to_phybundle(trim(output_name),trim(Diag(idx)%desc),trim(Diag(idx)%unit), "time: point",         &
                 axes(1:Diag(idx)%axes), fcst_grid, hist%nstt(idx), hist%levo(idx), phys_bundle(ibdl), outputfile(ibdl),  &
                 bdl_intplmethod(ibdl), rcd=rc)
            !           if( mpp_pe() == mpp_root_pe()) print *,'phys, add field,',trim(Diag(idx)%name),'idx=',idx,'ibdl=',ibdl
            !
            if( index(trim(Diag(idx)%intpl_method), "vector") > 0) then
              l2dvector = .true.
              if (hist%nstt_vctbl(idx) > 0) then
                output_name = 'wind'//trim(output_name)//'vector'
                outputfile1 = 'none'
                call hist%add_field_to_phybundle(trim(output_name),trim(Diag(idx)%desc),trim(Diag(idx)%unit), "time: point",       &
                     axes(1:Diag(idx)%axes), fcst_grid, hist%nstt_vctbl(idx), hist%levo(idx), phys_bundle(ibdl), outputfile1, &
                     bdl_intplmethod(ibdl),l2dvector=l2dvector,  rcd=rc)
                !               if( mpp_pe() == mpp_root_pe()) print *,'in phys, add vector field,',trim(Diag(idx)%name),' idx=',idx,' ibdl=',ibdl
              endif
            endif

          endif
        endif
      enddo
      if( .not. lput2physbdl ) then
        if( mpp_pe() == mpp_root_pe()) print *,'WARNING: not matching interpolation method, field ',trim(Diag(idx)%name), &
             ' is not added to phys bundle '
      endif

    enddo
    deallocate(hist%axis_name)
    deallocate(hist%all_axes)
    nullify(hist%axis_name)
    nullify(hist%all_axes)

  end subroutine history_type_bundle_setup

  !>@brief Adds one field to an ESMF field bundle for later output. Internal subroutine; do not call this directly.
  !> \section history_type%add_field_to_phybundle procedure
  !! This is part of the internal implementation of history_type_bundle_setup (history_type%bundle_setup).
  !! It sets attributes for and logs information about a single ESMF field. Do not call this subroutine directly.
  !! Call fv_phys_bundle_setup instead.
  subroutine history_type_add_field_to_phybundle(hist,var_name,long_name,units,cell_methods, axes,phys_grid, &
       kstt,levo,phys_bundle,output_file,intpl_method,range,l2dvector,rcd)
    !
    use esmf
    !
    implicit none
    class(history_type)                  :: hist
    character(*), intent(in)             :: var_name, long_name, units, cell_methods
    character(*), intent(in)             :: output_file, intpl_method
    integer, intent(in)                  :: axes(:)
    type(esmf_grid), intent(in)          :: phys_grid
    integer, intent(in)                  :: kstt, levo
    type(esmf_fieldbundle),intent(inout) :: phys_bundle
    real, intent(in), optional           :: range(2)
    logical, intent(in), optional        :: l2dvector
    integer, intent(out), optional       :: rcd
    !
    !*** local variable
    type(ESMF_Field)         :: field
    type(ESMF_DataCopy_Flag) :: copyflag=ESMF_DATACOPY_REFERENCE
    integer rc, i, j, idx
    real(4),dimension(:,:),pointer   :: temp_r2d
    real(4),dimension(:,:,:),pointer :: temp_r3d
    logical :: l2dvector_local
    !
    ! fix for non-standard compilers (e.g. PGI)
    l2dvector_local = .false.
    if (present(l2dvector)) then
      if (l2dvector) then
        l2dvector_local = .true.
      end if
    end if
    !
    !*** create esmf field
    if (l2dvector_local .and. size(axes)==2) then
      temp_r3d => hist%buffer_phys_windvect(1:3,hist%isco:hist%ieco,hist%jsco:hist%jeco,kstt)
      !     if( mpp_root_pe() == 0) print *,'phys, create wind vector esmf field'
      call ESMF_LogWrite('bf create winde vector esmf field '//trim(var_name), ESMF_LOGMSG_INFO, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      !datacopyflag=ESMF_DATACOPY_VALUE, &
      field = ESMF_FieldCreate(phys_grid, temp_r3d, datacopyflag=ESMF_DATACOPY_REFERENCE,          &
           gridToFieldMap=(/2,3/), ungriddedLBound=(/1/), ungriddedUBound=(/3/), &
           name=var_name, indexFlag=ESMF_INDEX_DELOCAL, rc=rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      call ESMF_LogWrite('af winde vector esmf field create '//trim(var_name), ESMF_LOGMSG_INFO, rc=rc)

      call ESMF_AttributeAdd(field, convention="NetCDF", purpose="FV3", &
           attrList=(/"output_file"/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
           line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

      call ESMF_AttributeSet(field, convention="NetCDF", purpose="FV3", &
           name='output_file',value=trim(output_file),rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
           line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

      call ESMF_LogWrite('before winde vector esmf field add output_file', ESMF_LOGMSG_INFO, rc=rc)

      !     if( mpp_root_pe() == 0)print *,'phys, aftercreate wind vector esmf field'
      call ESMF_FieldBundleAdd(phys_bundle,(/field/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
           line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

      if( present(rcd)) rcd=rc
      call ESMF_LogWrite('aft winde vector esmf field add to fieldbundle'//trim(var_name), ESMF_LOGMSG_INFO, rc=rc)
      return
    else if( trim(intpl_method) == 'nearest_stod' ) then
      if(size(axes) == 2) then
        temp_r2d => hist%buffer_phys_nb(hist%isco:hist%ieco,hist%jsco:hist%jeco,kstt)
        field = ESMF_FieldCreate(phys_grid, temp_r2d, datacopyflag=copyflag, &
             name=var_name, indexFlag=ESMF_INDEX_DELOCAL, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,       &
             line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

      else if(size(axes) == 3) then
        temp_r3d => hist%buffer_phys_nb(hist%isco:hist%ieco,hist%jsco:hist%jeco,kstt:kstt+levo-1)
        field = ESMF_FieldCreate(phys_grid, temp_r3d, datacopyflag=copyflag, &
             name=var_name, indexFlag=ESMF_INDEX_DELOCAL, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,       &
             line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

        if( mpp_pe() == mpp_root_pe()) print *,'add 3D field to after nearest_stod, fld=', trim(var_name)
      endif
    else if( trim(intpl_method) == 'bilinear' ) then
      if(size(axes) == 2) then
        temp_r2d => hist%buffer_phys_bl(hist%isco:hist%ieco,hist%jsco:hist%jeco,kstt)
        field = ESMF_FieldCreate(phys_grid, temp_r2d, datacopyflag=copyflag, &
             name=var_name, indexFlag=ESMF_INDEX_DELOCAL, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,       &
             line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
      else if(size(axes) == 3) then
        temp_r3d => hist%buffer_phys_bl(hist%isco:hist%ieco,hist%jsco:hist%jeco,kstt:kstt+levo-1)
        field = ESMF_FieldCreate(phys_grid, temp_r3d, datacopyflag=copyflag, &
             name=var_name, indexFlag=ESMF_INDEX_DELOCAL, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,       &
             line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
        if( mpp_pe() == mpp_root_pe()) print *,'add field to after bilinear, fld=', trim(var_name)
      endif
    endif
    !
    !*** add field attributes
    call ESMF_AttributeAdd(field, convention="NetCDF", purpose="FV3", &
         attrList=(/"long_name"/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
         line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

    call ESMF_AttributeSet(field, convention="NetCDF", purpose="FV3", &
         name='long_name',value=trim(long_name),rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
         line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

    call ESMF_AttributeAdd(field, convention="NetCDF", purpose="FV3", &
         attrList=(/"units"/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
         line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

    call ESMF_AttributeSet(field, convention="NetCDF", purpose="FV3", &
         name='units',value=trim(units),rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
         line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

    call ESMF_AttributeAdd(field, convention="NetCDF", purpose="FV3", &
         attrList=(/"missing_value"/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
         line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

    call ESMF_AttributeSet(field, convention="NetCDF", purpose="FV3", &
         name='missing_value',value=real(missing_value,kind=4),rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
         line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

    call ESMF_AttributeAdd(field, convention="NetCDF", purpose="FV3", &
         attrList=(/"_FillValue"/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
         line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

    call ESMF_AttributeSet(field, convention="NetCDF", purpose="FV3", &
         name='_FillValue',value=real(missing_value,kind=4),rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
         line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

    call ESMF_AttributeAdd(field, convention="NetCDF", purpose="FV3", &
         attrList=(/"cell_methods"/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
         line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

    call ESMF_AttributeSet(field, convention="NetCDF", purpose="FV3", &
         name='cell_methods',value=trim(cell_methods),rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
         line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
    !
    call ESMF_AttributeAdd(field, convention="NetCDF", purpose="FV3", &
         attrList=(/"output_file"/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
         line=__LINE__, file=__FILE__))  call ESMF_Finalize(endflag=ESMF_END_ABORT)

    call ESMF_AttributeSet(field, convention="NetCDF", purpose="FV3", &
         name='output_file',value=trim(output_file),rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
         line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

    !
    !*** add vertical coord attribute:
    if( size(axes) > 2) then
      do i=3,size(axes)
        idx=0
        do j=1,size(hist%all_axes)
          if (axes(i)==hist%all_axes(j)) then
            idx=j
            exit
          endif
        enddo
        if (idx>0) then
          call ESMF_AttributeAdd(field, convention="NetCDF", purpose="FV3", &
               attrList=(/"ESMF:ungridded_dim_labels"/), rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          if (levo == hist%nsoil_lsm) then
            call ESMF_AttributeSet(field, convention="NetCDF", purpose="FV3", &
                 name="ESMF:ungridded_dim_labels", valueList=(/"zsoil"/), rc=rc)
          else
            call ESMF_AttributeSet(field, convention="NetCDF", purpose="FV3", &
                 name="ESMF:ungridded_dim_labels", valueList=(/trim(hist%axis_name(idx))/), rc=rc)
          endif
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
        endif
      enddo
    endif

    !*** add field into bundle
    call ESMF_FieldBundleAdd(phys_bundle,(/field/), rc=rc)
    if( present(rcd)) rcd=rc
    !
    call ESMF_LogWrite('phys field add to fieldbundle '//trim(var_name), ESMF_LOGMSG_INFO, rc=rc)

  end subroutine history_type_add_field_to_phybundle

  !>@brief Private subroutine to search a field list for a specific name.
  !> \section history_type%find_output_name procedure
  !! Searches the GFS_Diagnostic-generated field list for a
  !! specific name and retrieves the name that should be used for
  !! outputting the variable. This is part of the internal
  !! implementation of history_type_bundle_setup
  !! (history_type%bundle_setup) and should not be called
  !! directly. Call fv_phys_bundle_setup instead.
  subroutine history_type_find_output_name(hist,module_name,field_name,output_name)
    implicit none
    class(history_type)          :: hist
    character(*), intent(in)     :: module_name
    character(*), intent(in)     :: field_name
    character(*), intent(out)    :: output_name
    !
    integer i,in_num
    integer tile_count
    !
    tile_count = 1
    in_num = find_input_field(module_name, field_name, tile_count)
    !
    output_name = ''
    do i=1, max_output_fields
      if(output_fields(i)%input_field == in_num) then
        output_name = output_fields(i)%output_name
        exit
      endif
    enddo
    if(output_name == '') then
19    format("Error: can't find output name for model field ",'"',A,'"')
      print 19,trim(field_name)
    endif

  end subroutine history_type_find_output_name
#endif
  !-------------------------------------------------------------------------

end module fv3atm_history_io_mod
!> @}
