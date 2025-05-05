!> \file fv3atm_sfc_io.F90
!! This file contains a derived type and subroutines to read and write restart files for
!! most FV3ATM surface fields. It works both for quilt (via ESMF) and non-quilt (via FMS)
!! restarts. Certain fields are handled by other files: fv3atm_oro_io.F90, fv3atm_rrfs_sd_io.F90,
!! and fv3atm_clm_lake_io.F90.
module fv3atm_sfc_io

  use block_control_mod,  only: block_control_type
  use fms2_io_mod,        only: FmsNetcdfDomainFile_t, unlimited, write_data,&
                                register_axis, register_restart_field,       &
                                register_variable_attribute, register_field, &
                                get_global_io_domain_indices, variable_exists, &
                                get_dimension_size
  use fv3atm_common_io,   only: GFS_Data_transfer, axis_type, &
       create_2d_field_and_add_to_bundle, create_3d_field_and_add_to_bundle
  use GFS_typedefs,       only: GFS_sfcprop_type, GFS_control_type, kind_phys
  use mpp_mod,            only: mpp_error,  NOTE
  use physcons,           only: con_tice          !saltwater freezing temp (K)

  implicit none
  private

  public :: Sfc_io_data_type
  public :: Sfc_io_fill_2d_names, Sfc_io_fill_2d_names_v2, &
       Sfc_io_fill_3d_names, Sfc_io_allocate_arrays, &
       Sfc_io_register_axes, Sfc_io_write_axes, Sfc_io_register_2d_fields, &
       Sfc_io_register_3d_fields, Sfc_io_copy_to_grid, Sfc_io_copy_from_grid, &
       Sfc_io_apply_safeguards, Sfc_io_transfer, Sfc_io_final

  !> \defgroup fv3atm_sfc_io module
  !> @{

  !>@ Minimum temperature allowed for snow/ice
  real(kind=kind_phys), parameter :: timin = 173.0_kind_phys

  real(kind_phys), parameter:: min_lake_orog = 200.0_kind_phys
  real(kind_phys), parameter:: zero = 0, one = 1

  !> Internal data storage type for reading and writing surface restart files
  type Sfc_io_data_type
    integer, public :: nvar2o = 0
    integer, public :: nvar3 = 0
    integer, public :: nvar2r = 0
    integer, public :: nvar2mp = 0
    integer, public :: nvar3mp = 0
    integer, public :: nvar2l = 0
    integer, public :: nvar2m = 0
    integer, public :: nvar_before_lake = 0

    ! The lsoil flag is only meaningful when reading:;
    logical, public :: is_lsoil = .false.

    logical, public :: is_v2_file = .false.

    ! SYNONYMS: Some nvar variables had two names in fv3atm_io.F90. They have
    ! only one name here. The "_s" is redundant because this file only has
    ! surface restart variables.
    !
    !  - nvar2m = nvar_s2m
    !  - nvar2o = nvar_s2o
    !  - nvar2r = nvar_s2r
    !  - nvar2mp = nvar_s2mp
    !  - nvar3mp = nvar_s3mp

    real(kind=kind_phys), pointer, dimension(:,:,:), public :: var2 => null()
    real(kind=kind_phys), pointer, dimension(:,:,:), public :: var3ice => null()
    real(kind=kind_phys), pointer, dimension(:,:,:,:), public :: var3 => null()
    real(kind=kind_phys), pointer, dimension(:,:,:,:), public :: var3sn => null()
    real(kind=kind_phys), pointer, dimension(:,:,:,:), public :: var3eq => null()
    real(kind=kind_phys), pointer, dimension(:,:,:,:), public :: var3zn => null()

    character(len=32), pointer, dimension(:), public :: name2 => null()
    character(len=32), pointer, dimension(:), public :: name3 => null()

  contains

    procedure, public :: allocate_arrays => Sfc_io_allocate_arrays
    procedure, public :: register_axes => Sfc_io_register_axes
    procedure, public :: write_axes => Sfc_io_write_axes
    procedure, public :: register_2d_fields => Sfc_io_register_2d_fields
    procedure, public :: register_3d_fields => Sfc_io_register_3d_fields
    procedure, public :: bundle_2d_fields => Sfc_io_bundle_2d_fields
    procedure, public :: bundle_3d_fields => Sfc_io_bundle_3d_fields
    procedure, public :: fill_2d_names => Sfc_io_fill_2d_names
    procedure, public :: fill_2d_names_v2 => Sfc_io_fill_2d_names_v2
    procedure, public :: fill_3d_names => Sfc_io_fill_3d_names
    procedure, public :: init_fields => Sfc_io_init_fields
    procedure, public :: transfer => Sfc_io_transfer
    procedure, public :: copy_to_grid => Sfc_io_copy_to_grid
    procedure, public :: copy_from_grid => Sfc_io_copy_from_grid
    procedure, public :: apply_safeguards => Sfc_io_apply_safeguards

    procedure, private :: calculate_indices => Sfc_io_calculate_indices

    final :: Sfc_io_final
  end type Sfc_io_data_type

contains

  !>@brief Calculates all nvar indices in the Sfc_io_data_type
  !> \section Sfc_io_data_type%calculate_indices() procedure
  !! Calculates all nvar counts, which record the number of fields
  !! of various types. These determine array sizes.
  !! Returns .true. if any nvar counts changed, or .false. otherwise.
  function Sfc_io_calculate_indices(sfc, Model, reading, warm_start)
    implicit none
    class(Sfc_io_data_type)             :: sfc
    type(GFS_control_type),      intent(in) :: Model
    logical :: Sfc_io_calculate_indices
    logical, intent(in) :: reading, warm_start

    integer :: nvar2m, nvar2o, nvar3, nvar2r, nvar2mp, nvar3mp, nvar2l
    integer :: nvar_before_lake

    nvar2m = 49
    if (Model%use_cice_alb .or. Model%lsm == Model%lsm_ruc) then
      nvar2m = nvar2m + 4
      !nvar2m = nvar2m + 5
    endif
    if (Model%cplwav) then
      nvar2m = nvar2m + 1
    endif
    if (Model%nstf_name(1) > 0) then
      nvar2o = 18
    else
      nvar2o = 0
    endif
    if (Model%lsm == Model%lsm_ruc .and. warm_start) then
      if (Model%rdlai) then
        nvar2r = 13
      else
        nvar2r = 12
      endif
      nvar3  = 5
    else
      if(.not.reading .and. Model%rdlai) then
        nvar2r = 1
      else
        nvar2r = 0
      endif
      nvar3  = 3
    endif
    if (Model%lsm == Model%lsm_noahmp) then
      nvar2mp = 29
      nvar3mp = 5
    else
      nvar2mp = 0
      nvar3mp = 0
    endif
    !CLM Lake and Flake
    if (Model%lkm > 0 .and. Model%iopt_lake==Model%iopt_lake_flake) then
      nvar2l = 10
    else
      nvar2l = 0
    endif

    nvar_before_lake=nvar2m+nvar2o+nvar2r+nvar2mp

    Sfc_io_calculate_indices = &
         nvar2m /= sfc%nvar2m .or. &
         nvar2o /= sfc%nvar2o .or. &
         nvar3 /= sfc%nvar3 .or. &
         nvar2r /= sfc%nvar2r .or. &
         nvar2mp /= sfc%nvar2mp .or. &
         nvar3mp /= sfc%nvar3mp .or. &
         nvar2l /= sfc%nvar2l .or. &
         nvar2m /= sfc%nvar2m .or. &
         nvar_before_lake /= sfc%nvar_before_lake

    sfc%nvar2m = nvar2m
    sfc%nvar2o = nvar2o
    sfc%nvar3 = nvar3
    sfc%nvar2r = nvar2r
    sfc%nvar2mp = nvar2mp
    sfc%nvar3mp = nvar3mp
    sfc%nvar2l = nvar2l
    sfc%nvar2m = nvar2m
    sfc%nvar_before_lake = nvar_before_lake

  end function Sfc_io_calculate_indices

  !>@brief Allocates internal Sfc_io_data_type arrays if array sizes should change.
  !> \section Sfc_io_data_type%allocate_arrays() procedure
  !! Calls calculate_arrays() to determine if any nvar counts have changed, based
  !! on the new arguments. If they have changed, then arrays are reallocated.
  !! The arrays will need to be filled with new data at that point, as the contents
  !! will be unknown. Returns .true. if arrays were reallocated, and .false. otherwise.
  function Sfc_io_allocate_arrays(sfc, Model, Atm_block, reading, warm_start)
    implicit none
    class(Sfc_io_data_type)             :: sfc
    type(block_control_type),    intent(in) :: Atm_block
    type(GFS_control_type),      intent(in) :: Model
    logical :: Sfc_io_allocate_arrays
    logical, intent(in) :: reading, warm_start

    integer :: isc, iec, jsc, jec, npz, nx, ny

    isc = Atm_block%isc
    iec = Atm_block%iec
    jsc = Atm_block%jsc
    jec = Atm_block%jec
    npz = Atm_block%npz
    nx  = (iec - isc + 1)
    ny  = (jec - jsc + 1)

    Sfc_io_allocate_arrays = sfc%calculate_indices(Model, reading, warm_start)
    Sfc_io_allocate_arrays = Sfc_io_allocate_arrays .or. .not. associated(sfc%name2)

    if(Sfc_io_allocate_arrays) then
      !--- allocate the various containers needed for restarts
      allocate(sfc%name2(sfc%nvar2m+sfc%nvar2o+sfc%nvar2mp+sfc%nvar2r+sfc%nvar2l))
      allocate(sfc%name3(0:sfc%nvar3+sfc%nvar3mp))
      allocate(sfc%var2(nx,ny,sfc%nvar2m+sfc%nvar2o+sfc%nvar2mp+sfc%nvar2r+sfc%nvar2l))

      ! Note that this may cause problems with RUC LSM for coldstart runs from GFS data
      ! if the initial conditions do contain this variable, because Model%kice is 9 for
      ! RUC LSM, but tiice in the initial conditions will only have two vertical layers
      allocate(sfc%var3ice(nx,ny,Model%kice))

      if (Model%lsm == Model%lsm_noah .or. Model%lsm == Model%lsm_noahmp .or. (.not.warm_start)) then
        allocate(sfc%var3(nx,ny,Model%lsoil,sfc%nvar3))
      elseif (Model%lsm == Model%lsm_ruc) then
        allocate(sfc%var3(nx,ny,Model%lsoil_lsm,sfc%nvar3))
      endif

      sfc%var2   = -9999.0_kind_phys
      sfc%var3   = -9999.0_kind_phys
      sfc%var3ice= -9999.0_kind_phys

      if (Model%lsm == Model%lsm_noahmp) then
        allocate(sfc%var3sn(nx,ny,-2:0,4:6))
        allocate(sfc%var3eq(nx,ny,1:4,7:7))
        allocate(sfc%var3zn(nx,ny,-2:4,8:8))

        sfc%var3sn = -9999.0_kind_phys
        sfc%var3eq = -9999.0_kind_phys
        sfc%var3zn = -9999.0_kind_phys
      endif
    endif

    if(Model%lkm>0 .and. Model%iopt_lake==Model%iopt_lake_flake .and. Model%me==0) then
      if(size(sfc%name2)/=sfc%nvar_before_lake+10) then
3814    format("ERROR: size mismatch size(sfc%name2)=",I0," /= nvar_before_lake+10=",I0)
        write(0,3814) size(sfc%name2),sfc%nvar_before_lake+10
      endif
    endif
  end function Sfc_io_allocate_arrays

  !>@ Registers all axes for reading or writing restarts using FMS (non-quilt)
  subroutine Sfc_io_register_axes(sfc, Model, Sfc_restart, reading, warm_start)
    implicit none
    class(Sfc_io_data_type)             :: sfc
    type(GFS_control_type),      intent(in) :: Model
    type(FmsNetcdfDomainFile_t) :: Sfc_restart
    logical, intent(in) :: reading, warm_start

    if(reading) then
      sfc%is_lsoil = .false.
    endif

    if(.not.warm_start .and. reading) then
      if( variable_exists(Sfc_restart,"lsoil") ) then
        if(reading) then
          sfc%is_lsoil=.true.
        endif
        call register_axis(Sfc_restart, 'lon', 'X')
        call register_axis(Sfc_restart, 'lat', 'Y')
        call register_axis(Sfc_restart, 'lsoil', dimension_length=Model%lsoil)
      else
        call register_axis(Sfc_restart, 'xaxis_1', 'X')
        call register_axis(Sfc_restart, 'yaxis_1', 'Y')
        call register_axis(Sfc_restart, 'zaxis_1', dimension_length=4)
        call register_axis(Sfc_restart, 'Time', 1)
      end if
    else
      call register_axis(Sfc_restart, 'xaxis_1', 'X')
      call register_axis(Sfc_restart, 'yaxis_1', 'Y')
      call register_axis(Sfc_restart, 'zaxis_1', dimension_length=Model%kice)
      if (Model%lsm == Model%lsm_noah .or. Model%lsm == Model%lsm_noahmp) then
        call register_axis(Sfc_restart, 'zaxis_2', dimension_length=Model%lsoil)
      else if(Model%lsm == Model%lsm_ruc .and. reading) then
        call register_axis(Sfc_restart, 'zaxis_2', dimension_length=Model%lsoil_lsm)
        ! The RUC only ever writes zaxis_1, which is combined soil/ice
        ! vertical dimension, lsoil_lsm/kice, which is 9.  Other LSMs read and
        ! write zaxis_2, which is lsoil for them, and that's always 4.
        ! Defining zaxis_2 here lets RUC LSM read from a different soil
        ! vertical coordinate (lsoil_lsm). It is needed for restart of RUC LSM
        ! from RUC LSM. This capability exists for historical reasons, because
        ! there are two sets of soil state variables: one set has lsoil=4
        ! vertical layers (Noah LSM. NoahMP LSM), and another set has
        ! lsoil_lsm=9 vertical levels (RUC LSM). Ideally there should be just
        ! one set of soil variables that could have different vertical
        ! dimension depending on the choice of LSM.  For now: just make sure
        ! you only restart RUC LSM off of RUC LSM, and always have kice =
        ! lsoil = lsoil_lsm = 9 and everything will be fine.
      endif
      if(Model%lsm == Model%lsm_noahmp) then
        call register_axis(Sfc_restart, 'zaxis_3', dimension_length=3)
        call register_axis(Sfc_restart, 'zaxis_4', dimension_length=7)
      end if
      call register_axis(Sfc_restart, 'Time', unlimited)
    endif
  end subroutine Sfc_io_register_axes

  !>@ Writes axis index variables and related metadata for all axes when writing FMS (non-quilt) restarts
  subroutine Sfc_io_write_axes(sfc, Model, Sfc_restart)
    implicit none
    class(Sfc_io_data_type)             :: sfc
    type(GFS_control_type),      intent(in) :: Model
    type(FmsNetcdfDomainFile_t) :: Sfc_restart

    integer, allocatable :: buffer(:)
    integer :: i, is, ie
    logical :: mand

    call register_field(Sfc_restart, 'xaxis_1', axis_type, (/'xaxis_1'/))
    call register_variable_attribute(Sfc_restart, 'xaxis_1', 'cartesian_axis', 'X', str_len=1)
    call get_global_io_domain_indices(Sfc_restart, 'xaxis_1', is, ie, indices=buffer)
    call write_data(Sfc_restart, "xaxis_1", buffer)
    deallocate(buffer)

    call register_field(Sfc_restart, 'yaxis_1', axis_type, (/'yaxis_1'/))
    call register_variable_attribute(Sfc_restart, 'yaxis_1', 'cartesian_axis', 'Y', str_len=1)
    call get_global_io_domain_indices(Sfc_restart, 'yaxis_1', is, ie, indices=buffer)
    call write_data(Sfc_restart, "yaxis_1", buffer)
    deallocate(buffer)

    call register_field(Sfc_restart, 'zaxis_1', axis_type, (/'zaxis_1'/))
    call register_variable_attribute(Sfc_restart, 'zaxis_1', 'cartesian_axis', 'Z', str_len=1)
    allocate( buffer(Model%kice) )
    do i=1, Model%kice
      buffer(i) = i
    end do
    call write_data(Sfc_restart, 'zaxis_1', buffer)
    deallocate(buffer)

    if (Model%lsm == Model%lsm_noah .or. Model%lsm == Model%lsm_noahmp) then
      call register_field(Sfc_restart, 'zaxis_2', axis_type, (/'zaxis_2'/))
      call register_variable_attribute(Sfc_restart, 'zaxis_2', 'cartesian_axis', 'Z', str_len=1)
      allocate( buffer(Model%lsoil) )
      do i=1, Model%lsoil
        buffer(i)=i
      end do
      call write_data(Sfc_restart, 'zaxis_2', buffer)
      deallocate(buffer)
    endif

    if(Model%lsm == Model%lsm_noahmp) then
      call register_field(Sfc_restart, 'zaxis_3', axis_type, (/'zaxis_3'/))
      call register_variable_attribute(Sfc_restart, 'zaxis_3', 'cartesian_axis', 'Z', str_len=1)
      allocate(buffer(3))
      do i=1, 3
        buffer(i) = i
      end do
      call write_data(Sfc_restart, 'zaxis_3', buffer)
      deallocate(buffer)

      call register_field(Sfc_restart, 'zaxis_4', axis_type, (/'zaxis_4'/))
      call register_variable_attribute(Sfc_restart, 'zaxis_4', 'cartesian_axis' ,'Z', str_len=1)
      allocate(buffer(7))
      do i=1, 7
        buffer(i)=i
      end do
      call write_data(Sfc_restart, 'zaxis_4', buffer)
      deallocate(buffer)
    end if
    call register_field(Sfc_restart, 'Time', axis_type, (/'Time'/))
    call register_variable_attribute(Sfc_restart, 'Time', 'cartesian_axis', 'T', str_len=1)
    call write_data( Sfc_restart, 'Time', 1)
  end subroutine Sfc_io_write_axes

  !>@ Fills the name3d array with all surface 3D field names.
  subroutine Sfc_io_fill_3d_names(sfc,Model,warm_start)
    implicit none
    class(Sfc_io_data_type)           :: sfc
    type(GFS_control_type),    intent(in) :: Model
    logical, intent(in) :: warm_start
    integer :: nt

    !--- names of the 3d variables to save
    if (Model%lsm == Model%lsm_noah .or. Model%lsm == Model%lsm_noahmp .or. (.not.warm_start)) then
      !--- names of the 3D variables to save
      sfc%name3(1) = 'stc'
      sfc%name3(2) = 'smc'
      sfc%name3(3) = 'slc'
      if (Model%lsm == Model%lsm_noahmp) then
        sfc%name3(4) = 'snicexy'
        sfc%name3(5) = 'snliqxy'
        sfc%name3(6) = 'tsnoxy'
        sfc%name3(7) = 'smoiseq'
        sfc%name3(8) = 'zsnsoxy'
      endif
    else if (Model%lsm == Model%lsm_ruc) then
      !--- names of the 3D variables to save
      sfc%name3(1) = 'tslb'
      sfc%name3(2) = 'smois'
      sfc%name3(3) = 'sh2o'
      sfc%name3(4) = 'smfr'
      sfc%name3(5) = 'flfr'
    end if
    sfc%name3(0) = 'tiice'
  end subroutine Sfc_io_fill_3d_names

  !>@ Fills the name2d array with all surface 2D field names. Updates nvar2m if needed.
  subroutine Sfc_io_fill_2d_names(sfc,Model,warm_start)
    implicit none
    class(Sfc_io_data_type)           :: sfc
    type(GFS_control_type),    intent(in) :: Model
    logical, intent(in) :: warm_start
    integer :: nt

    !--- names of the 2D variables to save
    nt=0
    nt=nt+1 ; sfc%name2(nt) = 'slmsk'
    nt=nt+1 ; sfc%name2(nt) = 'tsea'    !tsfc
    nt=nt+1 ; sfc%name2(nt) = 'sheleg'  !weasd
    nt=nt+1 ; sfc%name2(nt) = 'tg3'
    nt=nt+1 ; sfc%name2(nt) = 'zorl'
    nt=nt+1 ; sfc%name2(nt) = 'alvsf'
    nt=nt+1 ; sfc%name2(nt) = 'alvwf'
    nt=nt+1 ; sfc%name2(nt) = 'alnsf'
    nt=nt+1 ; sfc%name2(nt) = 'alnwf'
    nt=nt+1 ; sfc%name2(nt) = 'facsf'
    nt=nt+1 ; sfc%name2(nt) = 'facwf'
    nt=nt+1 ; sfc%name2(nt) = 'vfrac'
    nt=nt+1 ; sfc%name2(nt) = 'canopy'
    nt=nt+1 ; sfc%name2(nt) = 'f10m'
    nt=nt+1 ; sfc%name2(nt) = 't2m'
    nt=nt+1 ; sfc%name2(nt) = 'q2m'
    nt=nt+1 ; sfc%name2(nt) = 'vtype'
    nt=nt+1 ; sfc%name2(nt) = 'stype'
    nt=nt+1 ; sfc%name2(nt) = 'uustar'
    nt=nt+1 ; sfc%name2(nt) = 'ffmm'
    nt=nt+1 ; sfc%name2(nt) = 'ffhh'
    nt=nt+1 ; sfc%name2(nt) = 'hice'
    nt=nt+1 ; sfc%name2(nt) = 'fice'
    nt=nt+1 ; sfc%name2(nt) = 'tisfc'
    nt=nt+1 ; sfc%name2(nt) = 'tprcp'
    nt=nt+1 ; sfc%name2(nt) = 'srflag'
    nt=nt+1 ; sfc%name2(nt) = 'snwdph'  !snowd
    nt=nt+1 ; sfc%name2(nt) = 'shdmin'
    nt=nt+1 ; sfc%name2(nt) = 'shdmax'
    nt=nt+1 ; sfc%name2(nt) = 'slope'
    nt=nt+1 ; sfc%name2(nt) = 'snoalb'
    !--- variables below here are optional
    nt=nt+1 ; sfc%name2(nt) = 'scolor'
    nt=nt+1 ; sfc%name2(nt) = 'sncovr'
    nt=nt+1 ; sfc%name2(nt) = 'snodl' !snowd on land portion of a cell
    nt=nt+1 ; sfc%name2(nt) = 'weasdl'!weasd on land portion of a cell
    nt=nt+1 ; sfc%name2(nt) = 'tsfc'  !tsfc composite
    nt=nt+1 ; sfc%name2(nt) = 'tsfcl' !temp on land portion of a cell
    nt=nt+1 ; sfc%name2(nt) = 'zorlw' !zorl on water portion of a cell
    nt=nt+1 ; sfc%name2(nt) = 'zorll' !zorl on land portion of a cell
    nt=nt+1 ; sfc%name2(nt) = 'zorli' !zorl on ice portion of a cell
    nt=nt+1 ; sfc%name2(nt) = 'albdirvis_lnd'
    nt=nt+1 ; sfc%name2(nt) = 'albdirnir_lnd'
    nt=nt+1 ; sfc%name2(nt) = 'albdifvis_lnd'
    nt=nt+1 ; sfc%name2(nt) = 'albdifnir_lnd'
    nt=nt+1 ; sfc%name2(nt) = 'emis_lnd'
    nt=nt+1 ; sfc%name2(nt) = 'emis_ice'
    nt=nt+1 ; sfc%name2(nt) = 'sncovr_ice'
    nt=nt+1 ; sfc%name2(nt) = 'snodi' ! snowd on ice portion of a cell
    nt=nt+1 ; sfc%name2(nt) = 'weasdi'! weasd on ice portion of a cell

    if (Model%use_cice_alb .or. Model%lsm == Model%lsm_ruc) then
      nt=nt+1 ; sfc%name2(nt) = 'albdirvis_ice'
      nt=nt+1 ; sfc%name2(nt) = 'albdifvis_ice'
      nt=nt+1 ; sfc%name2(nt) = 'albdirnir_ice'
      nt=nt+1 ; sfc%name2(nt) = 'albdifnir_ice'
    endif

    if(Model%cplwav) then
      nt=nt+1 ; sfc%name2(nt) = 'zorlwav' !zorl from wave component
      sfc%nvar2m = nt
    endif

    if (Model%nstf_name(1) > 0) then
      !--- NSSTM inputs only needed when (nstf_name(1) > 0) .and. (nstf_name(2)) == 0)
      nt=nt+1 ; sfc%name2(nt) = 'tref'
      nt=nt+1 ; sfc%name2(nt) = 'z_c'
      nt=nt+1 ; sfc%name2(nt) = 'c_0'
      nt=nt+1 ; sfc%name2(nt) = 'c_d'
      nt=nt+1 ; sfc%name2(nt) = 'w_0'
      nt=nt+1 ; sfc%name2(nt) = 'w_d'
      nt=nt+1 ; sfc%name2(nt) = 'xt'
      nt=nt+1 ; sfc%name2(nt) = 'xs'
      nt=nt+1 ; sfc%name2(nt) = 'xu'
      nt=nt+1 ; sfc%name2(nt) = 'xv'
      nt=nt+1 ; sfc%name2(nt) = 'xz'
      nt=nt+1 ; sfc%name2(nt) = 'zm'
      nt=nt+1 ; sfc%name2(nt) = 'xtts'
      nt=nt+1 ; sfc%name2(nt) = 'xzts'
      nt=nt+1 ; sfc%name2(nt) = 'd_conv'
      nt=nt+1 ; sfc%name2(nt) = 'ifd'
      nt=nt+1 ; sfc%name2(nt) = 'dt_cool'
      nt=nt+1 ; sfc%name2(nt) = 'qrain'
    endif
    !
    ! Only needed when Noah MP LSM is used - 29 2D
    !
    if (Model%lsm == Model%lsm_noahmp) then
      nt=nt+1 ; sfc%name2(nt) = 'snowxy'
      nt=nt+1 ; sfc%name2(nt) = 'tvxy'
      nt=nt+1 ; sfc%name2(nt) = 'tgxy'
      nt=nt+1 ; sfc%name2(nt) = 'canicexy'
      nt=nt+1 ; sfc%name2(nt) = 'canliqxy'
      nt=nt+1 ; sfc%name2(nt) = 'eahxy'
      nt=nt+1 ; sfc%name2(nt) = 'tahxy'
      nt=nt+1 ; sfc%name2(nt) = 'cmxy'
      nt=nt+1 ; sfc%name2(nt) = 'chxy'
      nt=nt+1 ; sfc%name2(nt) = 'fwetxy'
      nt=nt+1 ; sfc%name2(nt) = 'sneqvoxy'
      nt=nt+1 ; sfc%name2(nt) = 'alboldxy'
      nt=nt+1 ; sfc%name2(nt) = 'qsnowxy'
      nt=nt+1 ; sfc%name2(nt) = 'wslakexy'
      nt=nt+1 ; sfc%name2(nt) = 'zwtxy'
      nt=nt+1 ; sfc%name2(nt) = 'waxy'
      nt=nt+1 ; sfc%name2(nt) = 'wtxy'
      nt=nt+1 ; sfc%name2(nt) = 'lfmassxy'
      nt=nt+1 ; sfc%name2(nt) = 'rtmassxy'
      nt=nt+1 ; sfc%name2(nt) = 'stmassxy'
      nt=nt+1 ; sfc%name2(nt) = 'woodxy'
      nt=nt+1 ; sfc%name2(nt) = 'stblcpxy'
      nt=nt+1 ; sfc%name2(nt) = 'fastcpxy'
      nt=nt+1 ; sfc%name2(nt) = 'xsaixy'
      nt=nt+1 ; sfc%name2(nt) = 'xlaixy'
      nt=nt+1 ; sfc%name2(nt) = 'taussxy'
      nt=nt+1 ; sfc%name2(nt) = 'smcwtdxy'
      nt=nt+1 ; sfc%name2(nt) = 'deeprechxy'
      nt=nt+1 ; sfc%name2(nt) = 'rechxy'
    else if (Model%lsm == Model%lsm_ruc .and. warm_start) then
      nt=nt+1 ; sfc%name2(nt) = 'wetness'
      nt=nt+1 ; sfc%name2(nt) = 'clw_surf_land'
      nt=nt+1 ; sfc%name2(nt) = 'clw_surf_ice'
      nt=nt+1 ; sfc%name2(nt) = 'qwv_surf_land'
      nt=nt+1 ; sfc%name2(nt) = 'qwv_surf_ice'
      nt=nt+1 ; sfc%name2(nt) = 'tsnow_land'
      nt=nt+1 ; sfc%name2(nt) = 'tsnow_ice'
      nt=nt+1 ; sfc%name2(nt) = 'snowfall_acc_land'
      nt=nt+1 ; sfc%name2(nt) = 'snowfall_acc_ice'
      nt=nt+1 ; sfc%name2(nt) = 'sfalb_lnd'
      nt=nt+1 ; sfc%name2(nt) = 'sfalb_lnd_bck'
      nt=nt+1 ; sfc%name2(nt) = 'sfalb_ice'
      if (Model%rdlai) then
        nt=nt+1 ; sfc%name2(nt) = 'lai'
      endif
    else if (Model%lsm == Model%lsm_ruc .and. Model%rdlai) then
      nt=nt+1 ; sfc%name2(nt) = 'lai'
    endif

    if (Model%lkm > 0 .and. Model%iopt_lake==Model%iopt_lake_flake) then
      nt=nt+1 ; sfc%name2(nt) = 'T_snow'
      nt=nt+1 ; sfc%name2(nt) = 'T_ice'
      nt=nt+1 ; sfc%name2(nt) = 'h_ML'
      nt=nt+1 ; sfc%name2(nt) = 't_ML'
      nt=nt+1 ; sfc%name2(nt) = 't_mnw'
      nt=nt+1 ; sfc%name2(nt) = 'h_talb'
      nt=nt+1 ; sfc%name2(nt) = 't_talb'
      nt=nt+1 ; sfc%name2(nt) = 't_bot1'
      nt=nt+1 ; sfc%name2(nt) = 't_bot2'
      nt=nt+1 ; sfc%name2(nt) = 'c_t'
    endif
  end subroutine Sfc_io_fill_2d_names

  !>@ Fills the name2d array with all surface 2D field names. Updates nvar2m if needed.
  !!  This routine is for v2 coldstart files.
  subroutine Sfc_io_fill_2d_names_v2(sfc,Model,warm_start)
    implicit none
    class(Sfc_io_data_type)           :: sfc
    type(GFS_control_type),    intent(in) :: Model
    logical, intent(in) :: warm_start
    integer :: nt

    !--- names of the 2D variables to save
    nt=0
    nt=nt+1 ; sfc%name2(nt) = 'slmsk'
    nt=nt+1 ; sfc%name2(nt) = 'tsea'    ! tsfc
    nt=nt+1 ; sfc%name2(nt) = 'sheleg'  ! weasd in file. Optional for cold starts.
    nt=nt+1 ; sfc%name2(nt) = 'tg3'
    nt=nt+1 ; sfc%name2(nt) = 'zorl'    ! Optional for cold starts.
    nt=nt+1 ; sfc%name2(nt) = 'alvsf'
    nt=nt+1 ; sfc%name2(nt) = 'alvwf'
    nt=nt+1 ; sfc%name2(nt) = 'alnsf'
    nt=nt+1 ; sfc%name2(nt) = 'alnwf'
    nt=nt+1 ; sfc%name2(nt) = 'facsf'
    nt=nt+1 ; sfc%name2(nt) = 'facwf'
    nt=nt+1 ; sfc%name2(nt) = 'vfrac'
    nt=nt+1 ; sfc%name2(nt) = 'canopy'
    nt=nt+1 ; sfc%name2(nt) = 'f10m'
    nt=nt+1 ; sfc%name2(nt) = 't2m'
    nt=nt+1 ; sfc%name2(nt) = 'q2m'
    nt=nt+1 ; sfc%name2(nt) = 'vtype'
    nt=nt+1 ; sfc%name2(nt) = 'stype'
    nt=nt+1 ; sfc%name2(nt) = 'uustar'
    nt=nt+1 ; sfc%name2(nt) = 'ffmm'
    nt=nt+1 ; sfc%name2(nt) = 'ffhh'
    nt=nt+1 ; sfc%name2(nt) = 'hice'
    nt=nt+1 ; sfc%name2(nt) = 'fice'
    nt=nt+1 ; sfc%name2(nt) = 'tisfc'
    nt=nt+1 ; sfc%name2(nt) = 'tprcp'
    nt=nt+1 ; sfc%name2(nt) = 'srflag'
    nt=nt+1 ; sfc%name2(nt) = 'snwdph'  ! snowd in file. Optional for cold starts.
    nt=nt+1 ; sfc%name2(nt) = 'shdmin'
    nt=nt+1 ; sfc%name2(nt) = 'shdmax'
    nt=nt+1 ; sfc%name2(nt) = 'slope'
    nt=nt+1 ; sfc%name2(nt) = 'snoalb'
    !--- variables below here are optional, unless indicated.
    nt=nt+1 ; sfc%name2(nt) = 'scolor'
    nt=nt+1 ; sfc%name2(nt) = 'sncovr'
    nt=nt+1 ; sfc%name2(nt) = 'snodl' ! snowd on land portion of a cell.
                                      ! Mandatory for cold starts.
    nt=nt+1 ; sfc%name2(nt) = 'weasdl'! weasd on land portion of a cell.
                                      ! Mandatory for cold starts.
    nt=nt+1 ; sfc%name2(nt) = 'tsfc'  ! tsfc composite
    nt=nt+1 ; sfc%name2(nt) = 'tsfcl' ! temp on land portion of a cell
    nt=nt+1 ; sfc%name2(nt) = 'zorlw' ! zorl on water portion of a cell
                                      ! Mandatory for cold starts.
    nt=nt+1 ; sfc%name2(nt) = 'zorll' ! zorl on land portion of a cell
    nt=nt+1 ; sfc%name2(nt) = 'zorli' ! zorl on ice portion of a cell
                                      ! Mandatory for cold starts.
    nt=nt+1 ; sfc%name2(nt) = 'albdirvis_lnd'
    nt=nt+1 ; sfc%name2(nt) = 'albdirnir_lnd'
    nt=nt+1 ; sfc%name2(nt) = 'albdifvis_lnd'
    nt=nt+1 ; sfc%name2(nt) = 'albdifnir_lnd'
    nt=nt+1 ; sfc%name2(nt) = 'emis_lnd'
    nt=nt+1 ; sfc%name2(nt) = 'emis_ice'
    nt=nt+1 ; sfc%name2(nt) = 'sncovr_ice'
    nt=nt+1 ; sfc%name2(nt) = 'snodi' ! snowd on ice portion of a cell.
                                      ! Mandatory for cold starts.
    nt=nt+1 ; sfc%name2(nt) = 'weasdi'! weasd on ice portion of a cell
                                      ! Mandatory for cold starts.

    if (Model%use_cice_alb .or. Model%lsm == Model%lsm_ruc) then
      nt=nt+1 ; sfc%name2(nt) = 'albdirvis_ice'
      nt=nt+1 ; sfc%name2(nt) = 'albdifvis_ice'
      nt=nt+1 ; sfc%name2(nt) = 'albdirnir_ice'
      nt=nt+1 ; sfc%name2(nt) = 'albdifnir_ice'
    endif

    if(Model%cplwav) then
      nt=nt+1 ; sfc%name2(nt) = 'zorlwav' !zorl from wave component
      sfc%nvar2m = nt
    endif

    if (Model%nstf_name(1) > 0) then
      !--- NSSTM inputs only needed when (nstf_name(1) > 0) .and. (nstf_name(2)) == 0)
      nt=nt+1 ; sfc%name2(nt) = 'tref'
      nt=nt+1 ; sfc%name2(nt) = 'z_c'
      nt=nt+1 ; sfc%name2(nt) = 'c_0'
      nt=nt+1 ; sfc%name2(nt) = 'c_d'
      nt=nt+1 ; sfc%name2(nt) = 'w_0'
      nt=nt+1 ; sfc%name2(nt) = 'w_d'
      nt=nt+1 ; sfc%name2(nt) = 'xt'
      nt=nt+1 ; sfc%name2(nt) = 'xs'
      nt=nt+1 ; sfc%name2(nt) = 'xu'
      nt=nt+1 ; sfc%name2(nt) = 'xv'
      nt=nt+1 ; sfc%name2(nt) = 'xz'
      nt=nt+1 ; sfc%name2(nt) = 'zm'
      nt=nt+1 ; sfc%name2(nt) = 'xtts'
      nt=nt+1 ; sfc%name2(nt) = 'xzts'
      nt=nt+1 ; sfc%name2(nt) = 'd_conv'
      nt=nt+1 ; sfc%name2(nt) = 'ifd'
      nt=nt+1 ; sfc%name2(nt) = 'dt_cool'
      nt=nt+1 ; sfc%name2(nt) = 'qrain'
    endif
    !
    ! Only needed when Noah MP LSM is used - 29 2D
    !
    if (Model%lsm == Model%lsm_noahmp) then
      nt=nt+1 ; sfc%name2(nt) = 'snowxy'
      nt=nt+1 ; sfc%name2(nt) = 'tvxy'
      nt=nt+1 ; sfc%name2(nt) = 'tgxy'
      nt=nt+1 ; sfc%name2(nt) = 'canicexy'
      nt=nt+1 ; sfc%name2(nt) = 'canliqxy'
      nt=nt+1 ; sfc%name2(nt) = 'eahxy'
      nt=nt+1 ; sfc%name2(nt) = 'tahxy'
      nt=nt+1 ; sfc%name2(nt) = 'cmxy'
      nt=nt+1 ; sfc%name2(nt) = 'chxy'
      nt=nt+1 ; sfc%name2(nt) = 'fwetxy'
      nt=nt+1 ; sfc%name2(nt) = 'sneqvoxy'
      nt=nt+1 ; sfc%name2(nt) = 'alboldxy'
      nt=nt+1 ; sfc%name2(nt) = 'qsnowxy'
      nt=nt+1 ; sfc%name2(nt) = 'wslakexy'
      nt=nt+1 ; sfc%name2(nt) = 'zwtxy'
      nt=nt+1 ; sfc%name2(nt) = 'waxy'
      nt=nt+1 ; sfc%name2(nt) = 'wtxy'
      nt=nt+1 ; sfc%name2(nt) = 'lfmassxy'
      nt=nt+1 ; sfc%name2(nt) = 'rtmassxy'
      nt=nt+1 ; sfc%name2(nt) = 'stmassxy'
      nt=nt+1 ; sfc%name2(nt) = 'woodxy'
      nt=nt+1 ; sfc%name2(nt) = 'stblcpxy'
      nt=nt+1 ; sfc%name2(nt) = 'fastcpxy'
      nt=nt+1 ; sfc%name2(nt) = 'xsaixy'
      nt=nt+1 ; sfc%name2(nt) = 'xlaixy'
      nt=nt+1 ; sfc%name2(nt) = 'taussxy'
      nt=nt+1 ; sfc%name2(nt) = 'smcwtdxy'
      nt=nt+1 ; sfc%name2(nt) = 'deeprechxy'
      nt=nt+1 ; sfc%name2(nt) = 'rechxy'
    else if (Model%lsm == Model%lsm_ruc .and. warm_start) then
      nt=nt+1 ; sfc%name2(nt) = 'wetness'
      nt=nt+1 ; sfc%name2(nt) = 'clw_surf_land'
      nt=nt+1 ; sfc%name2(nt) = 'clw_surf_ice'
      nt=nt+1 ; sfc%name2(nt) = 'qwv_surf_land'
      nt=nt+1 ; sfc%name2(nt) = 'qwv_surf_ice'
      nt=nt+1 ; sfc%name2(nt) = 'tsnow_land'
      nt=nt+1 ; sfc%name2(nt) = 'tsnow_ice'
      nt=nt+1 ; sfc%name2(nt) = 'snowfall_acc_land'
      nt=nt+1 ; sfc%name2(nt) = 'snowfall_acc_ice'
      nt=nt+1 ; sfc%name2(nt) = 'sfalb_lnd'
      nt=nt+1 ; sfc%name2(nt) = 'sfalb_lnd_bck'
      nt=nt+1 ; sfc%name2(nt) = 'sfalb_ice'
      if (Model%rdlai) then
        nt=nt+1 ; sfc%name2(nt) = 'lai'
      endif
    else if (Model%lsm == Model%lsm_ruc .and. Model%rdlai) then
      nt=nt+1 ; sfc%name2(nt) = 'lai'
    endif

    if (Model%lkm > 0 .and. Model%iopt_lake==Model%iopt_lake_flake) then
      nt=nt+1 ; sfc%name2(nt) = 'T_snow'
      nt=nt+1 ; sfc%name2(nt) = 'T_ice'
      nt=nt+1 ; sfc%name2(nt) = 'h_ML'
      nt=nt+1 ; sfc%name2(nt) = 't_ML'
      nt=nt+1 ; sfc%name2(nt) = 't_mnw'
      nt=nt+1 ; sfc%name2(nt) = 'h_talb'
      nt=nt+1 ; sfc%name2(nt) = 't_talb'
      nt=nt+1 ; sfc%name2(nt) = 't_bot1'
      nt=nt+1 ; sfc%name2(nt) = 't_bot2'
      nt=nt+1 ; sfc%name2(nt) = 'c_t'
    endif
  end subroutine Sfc_io_fill_2d_names_v2

  !>@ Registers 2D fields with FMS for reading or writing non-quilt restart files
  subroutine Sfc_io_register_2d_fields(sfc,Model,Sfc_restart,reading,warm_start)
    implicit none
    class(Sfc_io_data_type)             :: sfc
    type(GFS_control_type),      intent(in) :: Model
    type(FmsNetcdfDomainFile_t) :: Sfc_restart
    logical, intent(in) :: reading, warm_start

    real(kind=kind_phys), pointer, dimension(:,:)   :: var2_p  => NULL()
    real(kind=kind_phys), pointer, dimension(:,:,:) :: var3_p  => NULL()
    real(kind=kind_phys), pointer, dimension(:,:,:) :: var3_p1 => NULL()
    real(kind=kind_phys), pointer, dimension(:,:,:) :: var3_p2 => NULL()
    real(kind=kind_phys), pointer, dimension(:,:,:) :: var3_p3 => NULL()
    integer :: num
    logical :: mand

    character(len=7) :: time2d(3)

    integer :: xaxis_1_chunk, yaxis_1_chunk
    integer :: chunksizes2d(3)

    call get_dimension_size(Sfc_restart, 'xaxis_1', xaxis_1_chunk)
    call get_dimension_size(Sfc_restart, 'yaxis_1', yaxis_1_chunk)

    if(.not.reading) then
      time2d=(/'xaxis_1','yaxis_1','Time   '/)
      chunksizes2d=(/xaxis_1_chunk, yaxis_1_chunk, 1/)
    else
      time2d=(/'Time   ','yaxis_1','xaxis_1'/)
    endif

    !--- register the 2D fields
    if (sfc%is_v2_file) then
      do num = 1,sfc%nvar2m
        var2_p => sfc%var2(:,:,num)
        if (trim(sfc%name2(num)) == 'sncovr' .or. trim(sfc%name2(num)) == 'zorll' &
           .or. trim(sfc%name2(num)) == 'zorl' .or. trim(sfc%name2(num)) == 'zorlwav' &
           .or. trim(sfc%name2(num)) == 'snwdph' .or. trim(sfc%name2(num)) == 'sheleg'  &
           .or. trim(sfc%name2(num)) == 'tsfc' &
           .or. trim(sfc%name2(num)) == 'albdirvis_lnd' .or. trim(sfc%name2(num)) == 'albdirnir_lnd' &
           .or. trim(sfc%name2(num)) == 'albdifvis_lnd' .or. trim(sfc%name2(num)) == 'albdifnir_lnd' &
           .or. trim(sfc%name2(num)) == 'albdirvis_ice' .or. trim(sfc%name2(num)) == 'albdirnir_ice' &
           .or. trim(sfc%name2(num)) == 'albdifvis_ice' .or. trim(sfc%name2(num)) == 'albdifnir_ice' &
           .or. trim(sfc%name2(num)) == 'emis_lnd'      .or. trim(sfc%name2(num)) == 'emis_ice'      &
           .or. trim(sfc%name2(num)) == 'sncovr_ice'    .or. trim(sfc%name2(num)) == 'scolor') then
          if(reading .and. sfc%is_lsoil) then
            call register_restart_field(Sfc_restart, sfc%name2(num), var2_p, dimensions=(/'lat','lon'/), is_optional=.true.)
          else
            call register_restart_field(Sfc_restart, sfc%name2(num), var2_p, dimensions=time2d,&
               & chunksizes=chunksizes2d, is_optional=.true.)
          end if
        else
          if(reading .and. sfc%is_lsoil) then
            call register_restart_field(Sfc_restart,sfc%name2(num),var2_p, dimensions=(/'lat','lon'/))
          else
            call register_restart_field(Sfc_restart,sfc%name2(num),var2_p, dimensions=time2d, chunksizes=chunksizes2d)
          end if
        endif
      enddo
    else
      do num = 1,sfc%nvar2m
        var2_p => sfc%var2(:,:,num)
        if (trim(sfc%name2(num)) == 'sncovr' .or. trim(sfc%name2(num)) == 'tsfcl' .or. trim(sfc%name2(num)) == 'zorll' &
             .or. trim(sfc%name2(num)) == 'zorli' .or. trim(sfc%name2(num)) == 'zorlwav' &
             .or. trim(sfc%name2(num)) == 'snodl' .or. trim(sfc%name2(num)) == 'weasdl'  &
             .or. trim(sfc%name2(num)) == 'snodi' .or. trim(sfc%name2(num)) == 'weasdi'  &
             .or. trim(sfc%name2(num)) == 'tsfc'  .or. trim(sfc%name2(num)) ==  'zorlw'  &
             .or. trim(sfc%name2(num)) == 'albdirvis_lnd' .or. trim(sfc%name2(num)) == 'albdirnir_lnd' &
             .or. trim(sfc%name2(num)) == 'albdifvis_lnd' .or. trim(sfc%name2(num)) == 'albdifnir_lnd' &
             .or. trim(sfc%name2(num)) == 'albdirvis_ice' .or. trim(sfc%name2(num)) == 'albdirnir_ice' &
             .or. trim(sfc%name2(num)) == 'albdifvis_ice' .or. trim(sfc%name2(num)) == 'albdifnir_ice' &
             .or. trim(sfc%name2(num)) == 'emis_lnd'      .or. trim(sfc%name2(num)) == 'emis_ice'      &
             .or. trim(sfc%name2(num)) == 'sncovr_ice'    .or. trim(sfc%name2(num)) == 'scolor') then
          if(reading .and. sfc%is_lsoil) then
            call register_restart_field(Sfc_restart, sfc%name2(num), var2_p, dimensions=(/'lat','lon'/), is_optional=.true.)
          else
            call register_restart_field(Sfc_restart, sfc%name2(num), var2_p, dimensions=time2d,&
               & chunksizes=chunksizes2d, is_optional=.true.)
          end if
        else
          if(reading .and. sfc%is_lsoil) then
            call register_restart_field(Sfc_restart,sfc%name2(num),var2_p, dimensions=(/'lat','lon'/))
          else
            call register_restart_field(Sfc_restart,sfc%name2(num),var2_p, dimensions=time2d, chunksizes=chunksizes2d)
          end if
        endif
      enddo
    endif

    if (Model%nstf_name(1) > 0) then
      mand = .false.
      if (Model%nstf_name(2) == 0) mand = .true.
      do num = sfc%nvar2m+1,sfc%nvar2m+sfc%nvar2o
        var2_p => sfc%var2(:,:,num)
        if(sfc%is_lsoil) then
          call register_restart_field(Sfc_restart, sfc%name2(num), var2_p, dimensions=(/'lat','lon'/), is_optional=.not.mand)
        else
          call register_restart_field(Sfc_restart, sfc%name2(num), var2_p, dimensions=time2d, chunksizes=chunksizes2d, is_optional=.not.mand)
        endif
      enddo
    endif

    if (Model%lsm == Model%lsm_ruc) then ! sfc%nvar2mp = 0
      do num = sfc%nvar2m+sfc%nvar2o+1, sfc%nvar2m+sfc%nvar2o+sfc%nvar2r
        var2_p => sfc%var2(:,:,num)
        if(sfc%is_lsoil) then
          call register_restart_field(Sfc_restart, sfc%name2(num), var2_p, dimensions=(/'lat','lon'/) )
        else
          call register_restart_field(Sfc_restart, sfc%name2(num), var2_p, dimensions=time2d, chunksizes=chunksizes2d)
        end if
      enddo
    endif ! mp/ruc

    ! Noah MP register only necessary only lsm = 2, not necessary has values
    if ( (.not.reading .and. Model%lsm == Model%lsm_noahmp) &
         .or. (reading .and. sfc%nvar2mp > 0) ) then
      mand = .not.reading
      do num = sfc%nvar2m+sfc%nvar2o+1,sfc%nvar2m+sfc%nvar2o+sfc%nvar2mp
        var2_p => sfc%var2(:,:,num)
        if(sfc%is_lsoil) then
          call register_restart_field(Sfc_restart, sfc%name2(num), var2_p, dimensions=(/'lat','lon'/), is_optional=.not.mand)
        else
          call register_restart_field(Sfc_restart, sfc%name2(num), var2_p, dimensions=time2d, chunksizes=chunksizes2d, is_optional=.not.mand)
        end if
      enddo
    endif ! noahmp

    ! Flake
    if (Model%lkm > 0 .and. Model%iopt_lake==Model%iopt_lake_flake) then
      mand = .not.reading
      do num = sfc%nvar_before_lake+1,sfc%nvar_before_lake+sfc%nvar2l
        var2_p => sfc%var2(:,:,num)
        if(sfc%is_lsoil) then
          call register_restart_field(Sfc_restart, sfc%name2(num),var2_p,dimensions=(/'lat','lon'/), is_optional=.not.mand)
        else
          call register_restart_field(Sfc_restart, sfc%name2(num),var2_p,dimensions=time2d, chunksizes=chunksizes2d, is_optional=.not.mand)
        endif
      enddo
    endif

  end subroutine Sfc_io_register_2d_fields

  !>@ Registers 3D fields with FMS for reading or writing non-quilt restart files
  subroutine Sfc_io_register_3d_fields(sfc,Model,Sfc_restart,reading,warm_start)
    implicit none
    class(Sfc_io_data_type)             :: sfc
    type(GFS_control_type),      intent(in) :: Model
    type(FmsNetcdfDomainFile_t) :: Sfc_restart
    logical, intent(in) :: reading, warm_start

    real(kind=kind_phys), pointer, dimension(:,:,:) :: var3_p  => NULL()
    real(kind=kind_phys), pointer, dimension(:,:,:) :: var3_p1 => NULL()
    real(kind=kind_phys), pointer, dimension(:,:,:) :: var3_p2 => NULL()
    real(kind=kind_phys), pointer, dimension(:,:,:) :: var3_p3 => NULL()
    real(kind=kind_phys), pointer, dimension(:,:,:) :: var3_fr => NULL()
    integer :: num
    logical :: mand

    character(len=7), parameter :: xyz1_time(4) = (/'xaxis_1', 'yaxis_1', 'zaxis_1', 'Time   '/)
    character(len=7), parameter :: xyz2_time(4) = (/'xaxis_1', 'yaxis_1', 'zaxis_2', 'Time   '/)
    character(len=7), parameter :: xyz3_time(4) = (/'xaxis_1', 'yaxis_1', 'zaxis_3', 'Time   '/)
    character(len=7), parameter :: xyz4_time(4) = (/'xaxis_1', 'yaxis_1', 'zaxis_4', 'Time   '/)

    integer :: xaxis_1_chunk, yaxis_1_chunk
    integer :: chunksizes3d(4)

    call get_dimension_size(Sfc_restart, 'xaxis_1', xaxis_1_chunk)
    call get_dimension_size(Sfc_restart, 'yaxis_1', yaxis_1_chunk)

    chunksizes3d = (/xaxis_1_chunk, yaxis_1_chunk, 1, 1/)

    !--- register the 3D fields
    var3_p => sfc%var3ice(:,:,:)
    if (sfc%is_v2_file) then
      call register_restart_field(Sfc_restart, sfc%name3(0), var3_p, dimensions=xyz1_time, chunksizes=chunksizes3d, is_optional=.false.)
    else
      call register_restart_field(Sfc_restart, sfc%name3(0), var3_p, dimensions=xyz1_time, chunksizes=chunksizes3d, is_optional=.true.)
    endif

    if(reading) then
      do num = 1,sfc%nvar3
        var3_p => sfc%var3(:,:,:,num)
        if ( warm_start ) then
          call register_restart_field(Sfc_restart, sfc%name3(num), var3_p, dimensions=(/'xaxis_1', 'yaxis_1', 'lsoil  ', 'Time   '/),&
               &is_optional=.true.)
        else
          if(sfc%is_lsoil) then
            call register_restart_field(Sfc_restart, sfc%name3(num), var3_p, dimensions=(/'lat  ', 'lon  ', 'lsoil'/), is_optional=.true.)
          else
            call register_restart_field(Sfc_restart, sfc%name3(num), var3_p, dimensions=xyz2_time,&
                 &is_optional=.true.)
          end if
        end if
      enddo
    elseif(Model%lsm == Model%lsm_ruc) then
      do num = 1,sfc%nvar3
        var3_p => sfc%var3(:,:,:,num)
        call register_restart_field(Sfc_restart, sfc%name3(num), var3_p, dimensions=xyz1_time, chunksizes=chunksizes3d)
      enddo
      nullify(var3_p)
    else ! writing something other than ruc
      do num = 1,sfc%nvar3
        var3_p => sfc%var3(:,:,:,num)
        call register_restart_field(Sfc_restart, sfc%name3(num), var3_p, dimensions=xyz2_time, chunksizes=chunksizes3d)
      enddo
      nullify(var3_p)
    endif

    if (Model%lsm == Model%lsm_noahmp) then
      mand = .not.reading
      do num = sfc%nvar3+1,sfc%nvar3+3
        var3_p1 => sfc%var3sn(:,:,:,num)
        call register_restart_field(Sfc_restart, sfc%name3(num), var3_p1, dimensions=xyz3_time, chunksizes=chunksizes3d, is_optional=.not.mand)
      enddo

      var3_p2 => sfc%var3eq(:,:,:,7)
      call register_restart_field(Sfc_restart, sfc%name3(7), var3_p2, dimensions=xyz2_time, chunksizes=chunksizes3d, is_optional=.not.mand)

      var3_p3 => sfc%var3zn(:,:,:,8)
      call register_restart_field(Sfc_restart, sfc%name3(8), var3_p3, dimensions=xyz4_time, chunksizes=chunksizes3d, is_optional=.not.mand)
    endif   !mp

  end subroutine Sfc_io_register_3d_fields

  !>@ Initializes some surface fields with reasonable defaults
  subroutine Sfc_io_init_fields(sfc,Model)
    implicit none
    class(Sfc_io_data_type)           :: sfc
    type(GFS_control_type),    intent(in) :: Model

    !--- Noah MP define arbitrary value (number layers of snow) to indicate
    !coldstart(sfcfile doesn't include noah mp fields) or not

    if (Model%lsm == Model%lsm_noahmp) then
      sfc%var2(1,1,sfc%nvar2m+19) = -66666.0_kind_phys
    endif
  end subroutine Sfc_io_init_fields

  !>@ Copies data to the model grid (reading=true) or from the model grid (reading=false)
  !> \section Sfc_io_data_type%transfer
  !! Called to transfer data between the model grid and Sfc_io_data_type temporary arrays.
  !! The FMS and ESMF restarts use the temporary arrays, not the model grid arrays. This
  !! transfer routine copies to the model grid if reading=.true. or from the model grid
  !! if reading=.false. This is mostly loops around GFS_data_transfer() interface calls.
  !!
  !! In addition, if override_frac_grid is provided, it will be set to Model%frac_grid.
  subroutine Sfc_io_transfer(sfc, reading, Model, Atm_block, Sfcprop, warm_start, override_frac_grid)
    !--- interface variable definitions
    implicit none

    class(Sfc_io_data_type)                 :: sfc
    logical, intent(in)                     :: reading
    type(GFS_sfcprop_type)                  :: Sfcprop
    type(block_control_type),    intent(in) :: Atm_block
    type(GFS_control_type),      intent(in) :: Model
    logical, intent(in) :: warm_start
    logical, intent(out), optional :: override_frac_grid

    integer :: i, j, k, nb, ix, lsoil, num, nt, ixs, ixe, im
    integer :: isc, iec, jsc, jec, npz, nx, ny
    integer, allocatable :: ii1(:), jj1(:)
    real(kind_phys) :: ice

    ! "To" variable:
    !   to=.TRUE.   means       transfer sfc data  TO  Sfcprop grid
    !   to=.FALSE.  means  transfer into sfc data FROM Sfcprop grid

    isc = Atm_block%isc
    iec = Atm_block%iec
    jsc = Atm_block%jsc
    jec = Atm_block%jec
    npz = Atm_block%npz
    nx  = (iec - isc + 1)
    ny  = (jec - jsc + 1)

    !   write(0,*)' stype read in min,max=',minval(sfc_var2(:,:,35)),maxval(sfc_var2(:,:,35)),' sfc_name2=',sfc_name2(35)
    !   write(0,*)' stype read in min,max=',minval(sfc_var2(:,:,18)),maxval(sfc_var2(:,:,18))
    !   write(0,*)' sfc_var2=',sfc_var2(:,:,12)

    !$omp parallel do default(shared) private(i, j, nb, ix, nt, ii1, jj1, lsoil, ixs, ixe, im)
    block_loop: do nb = 1, Atm_block%nblks
      allocate(ii1(Atm_block%blksz(nb)))
      allocate(jj1(Atm_block%blksz(nb)))
      ii1=Atm_block%index(nb)%ii - isc + 1
      jj1=Atm_block%index(nb)%jj - jsc + 1
      ixs = Model%chunk_begin(nb)
      ixe = Model%chunk_end(nb)

      nt=0

      !--- 2D variables
      !    ------------
      call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%slmsk(ixs:ixe))   !--- slmsk
      call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%tsfco(ixs:ixe))   !--- tsfc (tsea in sfc file)
      call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%weasd(ixs:ixe))   !--- weasd (sheleg in sfc file)
      call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%tg3(ixs:ixe))     !--- tg3
      call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%zorl(ixs:ixe))    !--- zorl composite
      call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%alvsf(ixs:ixe))   !--- alvsf
      call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%alvwf(ixs:ixe))   !--- alvwf
      call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%alnsf(ixs:ixe))   !--- alnsf
      call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%alnwf(ixs:ixe))   !--- alnwf
      call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%facsf(ixs:ixe))   !--- facsf
      call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%facwf(ixs:ixe))   !--- facwf
      call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%vfrac(ixs:ixe))   !--- vfrac
      call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%canopy(ixs:ixe))  !--- canopy
      call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%f10m(ixs:ixe))    !--- f10m
      call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%t2m(ixs:ixe))     !--- t2m
      call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%q2m(ixs:ixe))     !--- q2m
      call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%vtype(ixs:ixe))   !--- vtype
      call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%stype(ixs:ixe))   !--- stype
      call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%uustar(ixs:ixe))  !--- uustar
      call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%ffmm(ixs:ixe))    !--- ffmm
      call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%ffhh(ixs:ixe))    !--- ffhh
      call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%hice(ixs:ixe))    !--- hice
      call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%fice(ixs:ixe))    !--- fice
      call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%tisfc(ixs:ixe))   !--- tisfc
      call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%tprcp(ixs:ixe))   !--- tprcp
      call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%srflag(ixs:ixe))  !--- srflag
      call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%snowd(ixs:ixe))   !--- snowd (snwdph in the file)
      call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%shdmin(ixs:ixe))  !--- shdmin
      call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%shdmax(ixs:ixe))  !--- shdmax
      call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%slope(ixs:ixe))   !--- slope
      call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%snoalb(ixs:ixe))  !--- snoalb
      call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%scolor(ixs:ixe))  !--- scolor
      call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%sncovr(ixs:ixe))  !--- sncovr
      call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%snodl(ixs:ixe))   !--- snodl (snowd on land  portion of a cell)
      call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%weasdl(ixs:ixe))  !--- weasdl (weasd on land  portion of a cell)
      call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%tsfc(ixs:ixe))    !--- tsfc composite
      call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%tsfcl(ixs:ixe))   !--- tsfcl  (temp on land portion of a cell)
      call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%zorlw(ixs:ixe))   !--- zorlw (zorl on water portion of a cell)
      call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%zorll(ixs:ixe))   !--- zorll (zorl on land portion of a cell)
      call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%zorli(ixs:ixe))   !--- zorli (zorl on ice  portion of a cell)
      call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%albdirvis_lnd(ixs:ixe))
      call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%albdirnir_lnd(ixs:ixe))
      call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%albdifvis_lnd(ixs:ixe))
      call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%albdifnir_lnd(ixs:ixe))
      call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%emis_lnd(ixs:ixe))
      call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%emis_ice(ixs:ixe))
      call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%sncovr_ice(ixs:ixe))
      call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%snodi(ixs:ixe))   !--- snodi (snowd on ice  portion of a cell)
      call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%weasdi(ixs:ixe))  !--- weasdi (weasd on ice  portion of a cell)
      if (Model%use_cice_alb .or. Model%lsm == Model%lsm_ruc) then
        call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%albdirvis_ice(ixs:ixe))
        call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%albdifvis_ice(ixs:ixe))
        call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%albdirnir_ice(ixs:ixe))
        call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%albdifnir_ice(ixs:ixe))
        !         call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%sfalb_ice(ixs:ixe))
      endif
      if(Model%cplwav) then
        !tgs - the following line is a bug. It should be nt = nt
        !nt = sfc%nvar2m-1 ! Next item will be at sfc%nvar2m
        call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%zorlwav(ixs:ixe)) !--- (zorl from wave model)
      else if(reading) then
        Sfcprop%zorlwav(ixs:ixe)  = Sfcprop%zorlw(ixs:ixe)
      endif

      if(present(override_frac_grid)) then
        override_frac_grid=Model%frac_grid
      endif

      if(reading) then
        do_lsi_fractions: do ix = 1, Atm_block%blksz(nb)
          im = Model%chunk_begin(nb)+ix-1
          if (Sfcprop%stype(im) == 14 .or. Sfcprop%stype(im) <= 0) then
            Sfcprop%landfrac(im) = zero
            Sfcprop%stype(im) = 0
            if (Sfcprop%lakefrac(im) > zero) then
              Sfcprop%lakefrac(im) = one
            endif
          endif

          if_frac_grid: if (Model%frac_grid) then
            if (Sfcprop%landfrac(im) > -999.0_kind_phys) then
              Sfcprop%slmsk(im) = ceiling(Sfcprop%landfrac(im)-1.0e-6)
              if (Sfcprop%slmsk(im) == 1 .and. Sfcprop%stype(im) == 14) &
                   Sfcprop%slmsk(im) = 0
              if (Sfcprop%lakefrac(im) > zero) then
                Sfcprop%oceanfrac(im) = zero ! lake & ocean don't coexist in a cell
                if (nint(Sfcprop%slmsk(im)) /= 1) then
                  if(Sfcprop%fice(im) >= Model%min_lakeice) then
                    Sfcprop%slmsk(im) = 2
                  else
                    Sfcprop%slmsk(im) = 0
                  endif
                endif
              else
                Sfcprop%lakefrac(im)  = zero
                Sfcprop%oceanfrac(im) = one - Sfcprop%landfrac(im)
                if (nint(Sfcprop%slmsk(im)) /= 1) then
                  if (Sfcprop%fice(im) >= Model%min_seaice) then
                    Sfcprop%slmsk(im) = 2
                  else
                    Sfcprop%slmsk(im) = 0
                  endif
                endif
              endif
            else
              if(present(override_frac_grid)) then
                override_frac_grid = .false.
              endif
              if (nint(Sfcprop%slmsk(im)) == 1) then
                Sfcprop%landfrac(im)  = one
                Sfcprop%lakefrac(im)  = zero
                Sfcprop%oceanfrac(im) = zero
              else
                if (Sfcprop%slmsk(im) < 0.1_kind_phys .or. Sfcprop%slmsk(im) > 1.9_kind_phys) then
                  Sfcprop%landfrac(im) = zero
                  if (Sfcprop%oro_uf(im) > min_lake_orog) then   ! lakes
                    Sfcprop%lakefrac(im)  = one
                    Sfcprop%oceanfrac(im) = zero
                  else                                               ! ocean
                    Sfcprop%lakefrac(im)  = zero
                    Sfcprop%oceanfrac(im) = one
                  endif
                endif
              endif
            endif
          else                                             ! not a fractional grid
            if (Sfcprop%landfrac(im) > -999.0_kind_phys) then
              if (Sfcprop%lakefrac(im) > zero) then
                Sfcprop%oceanfrac(im) = zero
                Sfcprop%landfrac(im)  = zero
                Sfcprop%lakefrac(im)  = one
                Sfcprop%slmsk(im)     = zero
                if (Sfcprop%fice(im) >= Model%min_lakeice) Sfcprop%slmsk(im) = 2.0
              else
                Sfcprop%slmsk(im) = nint(Sfcprop%landfrac(im))
                if (Sfcprop%stype(im) <= 0 .or. Sfcprop%stype(im) == 14) &
                     Sfcprop%slmsk(im) = zero
                if (nint(Sfcprop%slmsk(im)) == 0) then
                  Sfcprop%oceanfrac(im) = one
                  Sfcprop%landfrac(im)  = zero
                  Sfcprop%lakefrac(im)  = zero
                  if (Sfcprop%fice(im) >= Model%min_seaice) Sfcprop%slmsk(im) = 2.0
                else
                  Sfcprop%landfrac(im)  = one
                  Sfcprop%lakefrac(im)  = zero
                  Sfcprop%oceanfrac(im) = zero
                endif
              endif
            else
              if (nint(Sfcprop%slmsk(im)) == 1 .and. Sfcprop%stype(im) > 0      &
                   .and. Sfcprop%stype(im) /= 14) then
                Sfcprop%landfrac(im)  = one
                Sfcprop%lakefrac(im)  = zero
                Sfcprop%oceanfrac(im) = zero
              else
                Sfcprop%slmsk(im)    = zero
                Sfcprop%landfrac(im) = zero
                if (Sfcprop%oro_uf(im) > min_lake_orog) then   ! lakes
                  Sfcprop%lakefrac(im) = one
                  Sfcprop%oceanfrac(im) = zero
                  if (Sfcprop%fice(im) > Model%min_lakeice) Sfcprop%slmsk(im) = 2.0
                else                                       ! ocean
                  Sfcprop%lakefrac(im)  = zero
                  Sfcprop%oceanfrac(im) = one
                  if (Sfcprop%fice(im) > Model%min_seaice) Sfcprop%slmsk(im) = 2.0
                endif
              endif
            endif
          endif if_frac_grid
        enddo do_lsi_fractions
      endif

      if (reading .and. warm_start .and. Model%kdt > 1) then
        do ix = 1, Atm_block%blksz(nb)
          im = Model%chunk_begin(nb)+ix-1
          Sfcprop%slmsk(im)  = sfc%var2(ii1(ix),jj1(ix),1)    !--- slmsk
        enddo
      endif

      !
      !--- NSSTM variables
      !tgs - the following line is a bug that will show if(Model%cplwav) = true
      !nt = sfc%nvar2m
      if (Model%nstf_name(1) > 0) then
        if (reading .and. Model%nstf_name(2) == 1) then             ! nsst spinup
          !--- nsstm tref
          nt = nt + 18
          Sfcprop%tref(ixs:ixe)    = Sfcprop%tsfco(ixs:ixe)
          Sfcprop%z_c(ixs:ixe)     = zero
          Sfcprop%c_0(ixs:ixe)     = zero
          Sfcprop%c_d(ixs:ixe)     = zero
          Sfcprop%w_0(ixs:ixe)     = zero
          Sfcprop%w_d(ixs:ixe)     = zero
          Sfcprop%xt(ixs:ixe)      = zero
          Sfcprop%xs(ixs:ixe)      = zero
          Sfcprop%xu(ixs:ixe)      = zero
          Sfcprop%xv(ixs:ixe)      = zero
          Sfcprop%xz(ixs:ixe)      = 20.0_kind_phys
          Sfcprop%zm(ixs:ixe)      = zero
          Sfcprop%xtts(ixs:ixe)    = zero
          Sfcprop%xzts(ixs:ixe)    = zero
          Sfcprop%d_conv(ixs:ixe)  = zero
          Sfcprop%ifd(ixs:ixe)     = zero
          Sfcprop%dt_cool(ixs:ixe) = zero
          Sfcprop%qrain(ixs:ixe)   = zero
        elseif (.not.reading .or. Model%nstf_name(2) == 0) then         ! nsst restart
          call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%tref(ixs:ixe))  !--- nsstm tref
          call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%z_c(ixs:ixe))  !--- nsstm z_c
          call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%c_0(ixs:ixe))  !--- nsstm c_0
          call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%c_d(ixs:ixe))  !--- nsstm c_d
          call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%w_0(ixs:ixe))  !--- nsstm w_0
          call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%w_d(ixs:ixe))  !--- nsstm w_d
          call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%xt(ixs:ixe))  !--- nsstm xt
          call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%xs(ixs:ixe))  !--- nsstm xs
          call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%xu(ixs:ixe))  !--- nsstm xu
          call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%xv(ixs:ixe)) !--- nsstm xv
          call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%xz(ixs:ixe)) !--- nsstm xz
          call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%zm(ixs:ixe)) !--- nsstm zm
          call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%xtts(ixs:ixe)) !--- nsstm xtts
          call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%xzts(ixs:ixe)) !--- nsstm xzts
          call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%d_conv(ixs:ixe)) !--- nsstm d_conv
          call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%ifd(ixs:ixe)) !--- nsstm ifd
          call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%dt_cool(ixs:ixe)) !--- nsstm dt_cool
          call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%qrain(ixs:ixe)) !--- nsstm qrain
        endif
      endif

      if (Model%lsm == Model%lsm_ruc .and. (warm_start .or. .not. reading)) then
        !--- Extra RUC variables
        call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%wetness(ixs:ixe))
        call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%clw_surf_land(ixs:ixe))
        call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%clw_surf_ice(ixs:ixe))
        call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%qwv_surf_land(ixs:ixe))
        call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%qwv_surf_ice(ixs:ixe))
        call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%tsnow_land(ixs:ixe))
        call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%tsnow_ice(ixs:ixe))
        call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%snowfallac_land(ixs:ixe))
        call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%snowfallac_ice(ixs:ixe))
        call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%sfalb_lnd(ixs:ixe))
        call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%sfalb_lnd_bck(ixs:ixe))
        call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%sfalb_ice(ixs:ixe))
        if (Model%rdlai) then
          call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%xlaixy(ixs:ixe))
        endif
      else if (reading .and. Model%lsm == Model%lsm_ruc) then
        ! Initialize RUC snow cover on ice from snow cover
        Sfcprop%sncovr_ice(ixs:ixe) = Sfcprop%sncovr(ixs:ixe)
        if (Model%rdlai) then
          call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%xlaixy(ixs:ixe))
        end if
      elseif (Model%lsm == Model%lsm_noahmp) then
        !--- Extra Noah MP variables
        call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%snowxy(ixs:ixe))
        call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%tvxy(ixs:ixe))
        call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%tgxy(ixs:ixe))
        call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%canicexy(ixs:ixe))
        call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%canliqxy(ixs:ixe))
        call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%eahxy(ixs:ixe))
        call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%tahxy(ixs:ixe))
        call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%cmxy(ixs:ixe))
        call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%chxy(ixs:ixe))
        call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%fwetxy(ixs:ixe))
        call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%sneqvoxy(ixs:ixe))
        call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%alboldxy(ixs:ixe))
        call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%qsnowxy(ixs:ixe))
        call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%wslakexy(ixs:ixe))
        call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%zwtxy(ixs:ixe))
        call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%waxy(ixs:ixe))
        call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%wtxy(ixs:ixe))
        call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%lfmassxy(ixs:ixe))
        call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%rtmassxy(ixs:ixe))
        call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%stmassxy(ixs:ixe))
        call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%woodxy(ixs:ixe))
        call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%stblcpxy(ixs:ixe))
        call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%fastcpxy(ixs:ixe))
        call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%xsaixy(ixs:ixe))
        call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%xlaixy(ixs:ixe))
        call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%taussxy(ixs:ixe))
        call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%smcwtdxy(ixs:ixe))
        call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%deeprechxy(ixs:ixe))
        call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%rechxy(ixs:ixe))
      endif
      if (Model%lkm > 0 .and. Model%iopt_lake==Model%iopt_lake_flake) then
        call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%T_snow(ixs:ixe))
        call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%T_ice(ixs:ixe))
        call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%h_ML(ixs:ixe))
        call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%t_ML(ixs:ixe))
        call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%t_mnw(ixs:ixe))
        call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%h_talb(ixs:ixe))
        call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%t_talb(ixs:ixe))
        call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%t_bot1(ixs:ixe))
        call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%t_bot2(ixs:ixe))
        call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,sfc%var2,Sfcprop%c_t(ixs:ixe))
      endif
      if(.not.reading) then
        do k = 1,Model%kice
          do ix = 1, Atm_block%blksz(nb)
            im = Model%chunk_begin(nb)+ix-1
            ice=Sfcprop%tiice(im,k)
            if(ice<one) then
              sfc%var3ice(ii1(ix),jj1(ix),k) = zero
            else
              sfc%var3ice(ii1(ix),jj1(ix),k) = ice
            endif
          enddo
        enddo
      endif

      if (Model%lsm == Model%lsm_noah .or. Model%lsm == Model%lsm_noahmp .or. (reading .and. .not.warm_start)) then
        !--- 3D variables
        nt=0
        call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,1,Model%lsoil,sfc%var3,Sfcprop%stc(ixs:ixe,:))
        call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,1,Model%lsoil,sfc%var3,Sfcprop%smc(ixs:ixe,:))
        call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,1,Model%lsoil,sfc%var3,Sfcprop%slc(ixs:ixe,:))

        if (Model%lsm == Model%lsm_noahmp) then

          ! These use weird indexing which is lost during a Fortran
          ! subroutine call, so we use loops instead:

          if(reading) then
            nt=nt+1
            do lsoil = -2, 0
              do ix = 1, Atm_block%blksz(nb)
                im = Model%chunk_begin(nb)+ix-1
                Sfcprop%snicexy(im,lsoil) = sfc%var3sn(ii1(ix),jj1(ix),lsoil,nt)
              enddo
            enddo

            nt=nt+1
            do lsoil = -2, 0
              do ix = 1, Atm_block%blksz(nb)
                im = Model%chunk_begin(nb)+ix-1
                Sfcprop%snliqxy(im,lsoil) = sfc%var3sn(ii1(ix),jj1(ix),lsoil,nt)
              enddo
            enddo

            nt=nt+1
            do lsoil = -2, 0
              do ix = 1, Atm_block%blksz(nb)
                im = Model%chunk_begin(nb)+ix-1
                Sfcprop%tsnoxy(im,lsoil)  = sfc%var3sn(ii1(ix),jj1(ix),lsoil,nt)
              enddo
            enddo

            nt=nt+1
            do lsoil = 1, 4
              do ix = 1, Atm_block%blksz(nb)
                im = Model%chunk_begin(nb)+ix-1
                Sfcprop%smoiseq(im,lsoil)  = sfc%var3eq(ii1(ix),jj1(ix),lsoil,nt)
              enddo
            enddo

            nt=nt+1
            do lsoil = -2, 4
              do ix = 1, Atm_block%blksz(nb)
                im = Model%chunk_begin(nb)+ix-1
                Sfcprop%zsnsoxy(im,lsoil)  = sfc%var3zn(ii1(ix),jj1(ix),lsoil,nt)
              enddo
            enddo

          else

            nt=nt+1
            do lsoil = -2,0
              do ix = 1, Atm_block%blksz(nb)
                im = Model%chunk_begin(nb)+ix-1
                sfc%var3sn(ii1(ix),jj1(ix),lsoil,nt) = Sfcprop%snicexy(im,lsoil)
              enddo
            enddo

            nt=nt+1
            do lsoil = -2,0
              do ix = 1, Atm_block%blksz(nb)
                im = Model%chunk_begin(nb)+ix-1
                sfc%var3sn(ii1(ix),jj1(ix),lsoil,nt) = Sfcprop%snliqxy(im,lsoil)
              enddo
            enddo

            nt=nt+1
            do lsoil = -2,0
              do ix = 1, Atm_block%blksz(nb)
                im = Model%chunk_begin(nb)+ix-1
                sfc%var3sn(ii1(ix),jj1(ix),lsoil,nt) = Sfcprop%tsnoxy(im,lsoil)
              enddo
            enddo

            nt=nt+1
            do lsoil = 1,Model%lsoil
              do ix = 1, Atm_block%blksz(nb)
                im = Model%chunk_begin(nb)+ix-1
                sfc%var3eq(ii1(ix),jj1(ix),lsoil,nt)  = Sfcprop%smoiseq(im,lsoil)
              enddo
            enddo

            nt=nt+1
            do lsoil = -2,4
              do ix = 1, Atm_block%blksz(nb)
                im = Model%chunk_begin(nb)+ix-1
                sfc%var3zn(ii1(ix),jj1(ix),lsoil,nt)  = Sfcprop%zsnsoxy(im,lsoil)
              enddo
            enddo
          endif
        endif
      else if (Model%lsm == Model%lsm_ruc) then
        !--- 3D variables
        nt=0
        call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,1,Model%lsoil_lsm,sfc%var3,Sfcprop%tslb(ixs:ixe,:))
        call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,1,Model%lsoil_lsm,sfc%var3,Sfcprop%smois(ixs:ixe,:))
        call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,1,Model%lsoil_lsm,sfc%var3,Sfcprop%sh2o(ixs:ixe,:))
        call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,1,Model%lsoil_lsm,sfc%var3,Sfcprop%keepsmfr(ixs:ixe,:))
        call GFS_Data_transfer(reading,ii1,jj1,isc,jsc,nt,1,Model%lsoil_lsm,sfc%var3,Sfcprop%flag_frsoil(ixs:ixe,:))
      endif

      if(reading) then
        do k = 1,Model%kice
          do ix = 1, Atm_block%blksz(nb)
            im = Model%chunk_begin(nb)+ix-1
            Sfcprop%tiice(im,k) = sfc%var3ice(ii1(ix),jj1(ix),k)   !--- internal ice temp
          enddo
        enddo
      endif

      deallocate(ii1,jj1)

    end do block_loop
  end subroutine Sfc_io_transfer

  !>@ Copies from Sfc_io_data_type internal arrays to the model grid by calling transfer() with reading=.true.
  subroutine Sfc_io_copy_to_grid(sfc, Model, Atm_block, Sfcprop, warm_start, override_frac_grid)
    !--- interface variable definitions
    implicit none

    class(Sfc_io_data_type)             :: sfc
    type(GFS_sfcprop_type),      intent(in) :: Sfcprop
    type(block_control_type),    intent(in) :: Atm_block
    type(GFS_control_type),      intent(in) :: Model
    logical, intent(in) :: warm_start
    logical, intent(out), optional :: override_frac_grid

    call sfc%transfer(.true.,Model, Atm_block, Sfcprop, warm_start, override_frac_grid)

  end subroutine Sfc_io_copy_to_grid

  !>@ Copies from the model grid to Sfc_io_data_type internal arrays by calling transfer() with reading=.false.
  subroutine Sfc_io_copy_from_grid(sfc, Model, Atm_block, Sfcprop)
    !--- interface variable definitions
    implicit none

    class(Sfc_io_data_type)             :: sfc
    type(GFS_sfcprop_type),      intent(in) :: Sfcprop
    type(block_control_type),    intent(in) :: Atm_block
    type(GFS_control_type),      intent(in) :: Model

    call sfc%transfer(.false., Model, Atm_block, Sfcprop, warm_start=.false.)

  end subroutine Sfc_io_copy_from_grid

  !>@ Calculates values and applies safeguards after reading restart data.
  subroutine Sfc_io_apply_safeguards(sfc, Model, Atm_block, Sfcprop)
    !--- interface variable definitions
    implicit none

    class(Sfc_io_data_type)             :: sfc
    type(GFS_sfcprop_type),      intent(in) :: Sfcprop
    type(block_control_type),    intent(in) :: Atm_block
    type(GFS_control_type),      intent(in) :: Model

    integer :: i, j, k, nb, ix, lsoil, num, nt, im
    integer :: isc, iec, jsc, jec, npz, nx, ny
    real(kind_phys) :: ice, tem, tem1

    isc = Atm_block%isc
    iec = Atm_block%iec
    jsc = Atm_block%jsc
    jec = Atm_block%jec
    npz = Atm_block%npz
    nx  = (iec - isc + 1)
    ny  = (jec - jsc + 1)

    ! so far: At cold start everything is 9999.0, warm start snowxy has values
    !         but the 3D of snow fields are not available because not allocated yet.
    !         ix,nb loops may be consolidate with the Noah MP isnowxy init
    !         restore traditional vars first,we need some of them to init snow fields
    !         snow depth to actual snow layers; so we can allocate and register
    !         note zsnsoxy is from -2:4 - isnowxy is from 0:-2, but we need
    !         exact snow layers to pass 3D fields correctly, snow layers are
    !         different fro grid to grid, we have to init point by point/grid.
    !         It has to be done after the weasd is available
    !         sfc%var2(1,1,32) is the first; we need this to allocate snow related fields

    i = Atm_block%index(1)%ii(1) - isc + 1
    j = Atm_block%index(1)%jj(1) - jsc + 1

    if (sfc%var2(i,j,32) < -9990.0_kind_phys) then
      if (Model%me == Model%master ) call mpp_error(NOTE, 'gfs_driver::surface_props_input - set init soil color')
      !$omp parallel do default(shared) private(nb, ix, im)
      do nb = 1, Atm_block%nblks
        do ix = 1, Atm_block%blksz(nb)
          im = Model%chunk_begin(nb)+ix-1
          if ( nint (Sfcprop%slmsk(im)) == 1 ) then  !including glacier
            Sfcprop%scolor(im)  = 4
          else
            Sfcprop%scolor(im)  = zero
          endif
        enddo
      enddo
    endif

    if (sfc%is_v2_file) then

    if (sfc%var2(i,j,27) < -9990.0_kind_phys) then
      if (Model%me == Model%master ) call mpp_error(NOTE, 'gfs_driver::surface_props_input - computing snowd')
      !$omp parallel do default(shared) private(nb, ix, im)
      do nb = 1, Atm_block%nblks
        do ix = 1, Atm_block%blksz(nb)
          im = Model%chunk_begin(nb)+ix-1
          if (Sfcprop%fice(im) > zero) then
            Sfcprop%snowd(im)  = Sfcprop%snodi(im)
          elseif (Sfcprop%landfrac(im) > zero) then
            Sfcprop%snowd(im)  = Sfcprop%snodl(im)
          else
            Sfcprop%snowd(im)  = zero
          endif
        enddo
      enddo
    endif

    if (sfc%var2(i,j,3) < -9990.0_kind_phys) then
      if (Model%me == Model%master ) call mpp_error(NOTE, 'gfs_driver::surface_props_input - computing weasd')
      !$omp parallel do default(shared) private(nb, ix, im)
      do nb = 1, Atm_block%nblks
        do ix = 1, Atm_block%blksz(nb)
          im = Model%chunk_begin(nb)+ix-1
          if (Sfcprop%fice(im) > zero) then
            Sfcprop%weasd(im) = Sfcprop%weasdi(im)
          elseif (Sfcprop%landfrac(im) > zero) then
            Sfcprop%weasd(im) = Sfcprop%weasdl(im)
          else
            Sfcprop%weasd(im) = zero
          endif
        enddo
      enddo
    endif

! Needed for first time step in radiation before Noah/NoahMP sets it from look up table.
! Just use a nominal value.
    if (sfc%var2(i,j,39) < -9990.0_kind_phys) then
      if (Model%me == Model%master ) call mpp_error(NOTE, 'gfs_driver::surface_props_input - computing zorll')
      !$omp parallel do default(shared) private(nb, ix, im)
      do nb = 1, Atm_block%nblks
        do ix = 1, Atm_block%blksz(nb)
          im = Model%chunk_begin(nb)+ix-1
          if (Sfcprop%landfrac(im) > zero) then
            Sfcprop%zorll(im) = 25.0
          endif
        enddo
      enddo
    endif

    if (sfc%var2(i,j,5) < -9990.0_kind_phys) then
      if (Model%me == Model%master ) call mpp_error(NOTE, 'gfs_driver::surface_props_input - computing zorl')
      !$omp parallel do default(shared) private(nb, ix, im)
      do nb = 1, Atm_block%nblks
        do ix = 1, Atm_block%blksz(nb)
          im = Model%chunk_begin(nb)+ix-1
          if (Sfcprop%fice(im) > zero) then
            Sfcprop%zorl(im) = Sfcprop%zorli(im)
          elseif (Sfcprop%landfrac(im) > zero) then
            Sfcprop%zorl(im) = Sfcprop%zorll(im)
          else
            Sfcprop%zorl(im) = Sfcprop%zorlw(im)
          endif
        enddo
      enddo
    endif

    if (sfc%var2(i,j,46) < -9990.0_kind_phys) then
      if (Model%me == Model%master ) call mpp_error(NOTE, 'gfs_driver::surface_props_input - computing emis_ice')
      !$omp parallel do default(shared) private(nb, ix, im)
      do nb = 1, Atm_block%nblks
        do ix = 1, Atm_block%blksz(nb)
          im = Model%chunk_begin(nb)+ix-1
          Sfcprop%emis_ice(im) = 0.96
        enddo
      enddo
    endif

    if (sfc%var2(i,j,47) < -9990.0_kind_phys .and. Model%lsm /= Model%lsm_ruc) then
      if (Model%me == Model%master ) call mpp_error(NOTE, 'gfs_driver::surface_props_input - computing sncovr_ice')
      !$omp parallel do default(shared) private(nb, ix, im)
      do nb = 1, Atm_block%nblks
        do ix = 1, Atm_block%blksz(nb)
          im = Model%chunk_begin(nb)+ix-1
          !Sfcprop%sncovr_ice(im) = Sfcprop%sncovr(im)
          Sfcprop%sncovr_ice(im) = zero
        enddo
      enddo
    endif

    if (Model%use_cice_alb) then
      if (sfc%var2(i,j,50) < -9990.0_kind_phys) then
        !$omp parallel do default(shared) private(nb, ix, im)
        do nb = 1, Atm_block%nblks
          do ix = 1, Atm_block%blksz(nb)
            im = Model%chunk_begin(nb)+ix-1
            if (Sfcprop%oceanfrac(im) > zero .and. &
                 Sfcprop%fice(im) >= Model%min_seaice) then
              Sfcprop%albdirvis_ice(im) = 0.6_kind_phys
              Sfcprop%albdifvis_ice(im) = 0.6_kind_phys
              Sfcprop%albdirnir_ice(im) = 0.6_kind_phys
              Sfcprop%albdifnir_ice(im) = 0.6_kind_phys
            endif
          enddo
        enddo
      endif

    endif

    ! Fill in composite tsfc for coldstart runs - must happen after tsfcl is computed
    compute_tsfc_for_coldstart: if (sfc%var2(i,j,36) < -9990.0_kind_phys) then
      if (Model%me == Model%master ) call mpp_error(NOTE, 'gfs_driver::surface_props_input - computing composite tsfc')
      if(Model%frac_grid) then ! 3-way composite
        !$omp parallel do default(shared) private(nb, ix, im, tem, tem1)
        do nb = 1, Atm_block%nblks
          do ix = 1, Atm_block%blksz(nb)
            im = Model%chunk_begin(nb)+ix-1
            Sfcprop%tsfco(im) = max(con_tice, Sfcprop%tsfco(im)) ! this may break restart reproducibility
            tem1 = one - Sfcprop%landfrac(im)
            tem  = tem1 * Sfcprop%fice(im) ! tem = ice fraction wrt whole cell
            Sfcprop%tsfc(im) = Sfcprop%tsfcl(im) * Sfcprop%landfrac(im) &
                 + Sfcprop%tisfc(im) * tem                      &
                 + Sfcprop%tsfco(im) * (tem1-tem)
          enddo
        enddo
      else
        !$omp parallel do default(shared) private(nb, ix, im, tem)
        do nb = 1, Atm_block%nblks
          do ix = 1, Atm_block%blksz(nb)
            im = Model%chunk_begin(nb)+ix-1
            if (Sfcprop%slmsk(im) == 1) then
              Sfcprop%tsfc(im) = Sfcprop%tsfcl(im)
              if (Sfcprop%tsfc(im) < -99 .or. Sfcprop%tsfc(im) > 999.) print*,'bad tsfc land ',nb,ix,Sfcprop%tsfcl(im)
            elseif(Sfcprop%fice(im) > 0.0)then
              Sfcprop%tsfc(im) = Sfcprop%tisfc(im)
              if (Sfcprop%tsfc(im) < -99 .or. Sfcprop%tsfc(im) > 999.) print*,'bad tsfc ice  ',nb,ix,Sfcprop%tisfc(im)
            else
              Sfcprop%tsfc(im) = Sfcprop%tsfco(im)
              if (Sfcprop%tsfc(im) < -99 .or. Sfcprop%tsfc(im) > 999.) print*,'bad tsfc water ',nb,ix,Sfcprop%tsfco(im)
            endif
          enddo
        enddo
      endif
    endif compute_tsfc_for_coldstart

    if (sfc%var2(i,j,sfc%nvar2m) < -9990.0_kind_phys) then
      if (Model%me == Model%master ) call mpp_error(NOTE, 'gfs_driver::surface_props_input - computing zorlwav')
      !$omp parallel do default(shared) private(nb, ix, im)
      do nb = 1, Atm_block%nblks
        do ix = 1, Atm_block%blksz(nb)
          im = Model%chunk_begin(nb)+ix-1
          Sfcprop%zorlwav(im) = Sfcprop%zorl(im) !--- compute zorlwav from existing variables
        enddo
      enddo
    endif


    else ! old verion of coldstart file


    if (sfc%var2(i,j,34) < -9990.0_kind_phys) then
      if (Model%me == Model%master ) call mpp_error(NOTE, 'gfs_driver::surface_props_input - computing snodl')
      !$omp parallel do default(shared) private(nb, ix, im, tem)
      do nb = 1, Atm_block%nblks
        do ix = 1, Atm_block%blksz(nb)
          im = Model%chunk_begin(nb)+ix-1
          if (Sfcprop%landfrac(im) > zero) then
            tem = one / (Sfcprop%fice(im)*(one-Sfcprop%landfrac(im))+Sfcprop%landfrac(im))
            Sfcprop%snodl(im)  = Sfcprop%snowd(im) * tem
          else
            Sfcprop%snodl(im)  = zero
          endif
        enddo
      enddo
    endif

    if (sfc%var2(i,j,35) < -9990.0_kind_phys) then
      if (Model%me == Model%master ) call mpp_error(NOTE, 'gfs_driver::surface_props_input - computing weasdl')
      !$omp parallel do default(shared) private(nb, ix, im, tem)
      do nb = 1, Atm_block%nblks
        do ix = 1, Atm_block%blksz(nb)
          im = Model%chunk_begin(nb)+ix-1
          if (Sfcprop%landfrac(im) > zero) then
            tem = one / (Sfcprop%fice(im)*(one-Sfcprop%landfrac(im))+Sfcprop%landfrac(im))
            Sfcprop%weasdl(im) = Sfcprop%weasd(im) * tem
          else
            Sfcprop%weasdl(im) = zero
          endif
        enddo
      enddo
    endif

    if (sfc%var2(i,j,37) < -9990.0_kind_phys) then
      if (Model%me == Model%master ) call mpp_error(NOTE, 'gfs_driver::surface_props_input - computing tsfcl')
      !$omp parallel do default(shared) private(nb, ix, im)
      do nb = 1, Atm_block%nblks
        do ix = 1, Atm_block%blksz(nb)
          im = Model%chunk_begin(nb)+ix-1
          Sfcprop%tsfcl(im) = Sfcprop%tsfco(im) !--- compute tsfcl from existing variables
        enddo
      enddo
    endif

    if (sfc%var2(i,j,38) < -9990.0_kind_phys) then
      if (Model%me == Model%master ) call mpp_error(NOTE, 'gfs_driver::surface_props_input - computing zorlw')
      !$omp parallel do default(shared) private(nb, ix, im)
      do nb = 1, Atm_block%nblks
        do ix = 1, Atm_block%blksz(nb)
          im = Model%chunk_begin(nb)+ix-1
          if (Sfcprop%landfrac(im) < one .and. Sfcprop%fice(im) < one) then
            Sfcprop%zorlw(im) = min(Sfcprop%zorl(im), 0.317)
          endif
        enddo
      enddo
    endif

    if (sfc%var2(i,j,39) < -9990.0_kind_phys) then
      if (Model%me == Model%master ) call mpp_error(NOTE, 'gfs_driver::surface_props_input - computing zorll')
      !$omp parallel do default(shared) private(nb, ix, im)
      do nb = 1, Atm_block%nblks
        do ix = 1, Atm_block%blksz(nb)
          im = Model%chunk_begin(nb)+ix-1
          Sfcprop%zorll(im) = Sfcprop%zorl(im) !--- compute zorll from existing variables
        enddo
      enddo
    endif

    if (sfc%var2(i,j,40) < -9990.0_kind_phys) then
      if (Model%me == Model%master ) call mpp_error(NOTE, 'gfs_driver::surface_props_input - computing zorli')
      !$omp parallel do default(shared) private(nb, ix, im)
      do nb = 1, Atm_block%nblks
        do ix = 1, Atm_block%blksz(nb)
          im = Model%chunk_begin(nb)+ix-1
          if (Sfcprop%fice(im)*(one-Sfcprop%landfrac(im)) > zero) then
            Sfcprop%zorli(im) = one
          endif
        enddo
      enddo
    endif

    if (sfc%var2(i,j,46) < -9990.0_kind_phys) then
      if (Model%me == Model%master ) call mpp_error(NOTE, 'gfs_driver::surface_props_input - computing emis_ice')
      !$omp parallel do default(shared) private(nb, ix, im)
      do nb = 1, Atm_block%nblks
        do ix = 1, Atm_block%blksz(nb)
          im = Model%chunk_begin(nb)+ix-1
          Sfcprop%emis_ice(im) = 0.96
        enddo
      enddo
    endif

    if (sfc%var2(i,j,47) < -9990.0_kind_phys .and. Model%lsm /= Model%lsm_ruc) then
      if (Model%me == Model%master ) call mpp_error(NOTE, 'gfs_driver::surface_props_input - computing sncovr_ice')
      !$omp parallel do default(shared) private(nb, ix, im)
      do nb = 1, Atm_block%nblks
        do ix = 1, Atm_block%blksz(nb)
          im = Model%chunk_begin(nb)+ix-1
          !Sfcprop%sncovr_ice(im) = Sfcprop%sncovr(im)
          Sfcprop%sncovr_ice(im) = zero
        enddo
      enddo
    endif

    if (sfc%var2(i,j,48) < -9990.0_kind_phys) then
      if (Model%me == Model%master ) call mpp_error(NOTE, 'gfs_driver::surface_props_input - computing snodi')
      !$omp parallel do default(shared) private(nb, ix, im, tem)
      do nb = 1, Atm_block%nblks
        do ix = 1, Atm_block%blksz(nb)
          im = Model%chunk_begin(nb)+ix-1
          if (Sfcprop%fice(im) > zero) then
            tem = one / (Sfcprop%fice(im)*(one-Sfcprop%landfrac(im))+Sfcprop%landfrac(im))
            Sfcprop%snodi(im)  = min(Sfcprop%snowd(im) * tem, 3.0)
          else
            Sfcprop%snodi(im)  = zero
          endif
        enddo
      enddo
    endif

    if (sfc%var2(i,j,49) < -9990.0_kind_phys) then
      if (Model%me == Model%master ) call mpp_error(NOTE, 'gfs_driver::surface_props_input - computing weasdi')
      !$omp parallel do default(shared) private(nb, ix, im, tem)
      do nb = 1, Atm_block%nblks
        do ix = 1, Atm_block%blksz(nb)
          im = Model%chunk_begin(nb)+ix-1
          if (Sfcprop%fice(im) > zero) then
            tem = one / (Sfcprop%fice(im)*(one-Sfcprop%landfrac(im))+Sfcprop%landfrac(im))
            Sfcprop%weasdi(im)  = Sfcprop%weasd(im)*tem
          else
            Sfcprop%weasdi(im)  = zero
          endif
        enddo
      enddo
    endif

    if (Model%use_cice_alb) then
      if (sfc%var2(i,j,50) < -9990.0_kind_phys) then
        !$omp parallel do default(shared) private(nb, ix, im)
        do nb = 1, Atm_block%nblks
          do ix = 1, Atm_block%blksz(nb)
            im = Model%chunk_begin(nb)+ix-1
            if (Sfcprop%oceanfrac(im) > zero .and. &
                 Sfcprop%fice(im) >= Model%min_seaice) then
              Sfcprop%albdirvis_ice(im) = 0.6_kind_phys
              Sfcprop%albdifvis_ice(im) = 0.6_kind_phys
              Sfcprop%albdirnir_ice(im) = 0.6_kind_phys
              Sfcprop%albdifnir_ice(im) = 0.6_kind_phys
            endif
          enddo
        enddo
      endif

    endif

    ! Fill in composite tsfc for coldstart runs - must happen after tsfcl is computed
    compute_tsfc_for_colstart: if (sfc%var2(i,j,35) < -9990.0_kind_phys) then
      if (Model%me == Model%master ) call mpp_error(NOTE, 'gfs_driver::surface_props_input - computing composite tsfc')
      if(Model%frac_grid) then ! 3-way composite
        !$omp parallel do default(shared) private(nb, ix, im, tem, tem1)
        do nb = 1, Atm_block%nblks
          do ix = 1, Atm_block%blksz(nb)
            im = Model%chunk_begin(nb)+ix-1
            Sfcprop%tsfco(im) = max(con_tice, Sfcprop%tsfco(im)) ! this may break restart reproducibility
            tem1 = one - Sfcprop%landfrac(im)
            tem  = tem1 * Sfcprop%fice(im) ! tem = ice fraction wrt whole cell
            Sfcprop%tsfc(im) = Sfcprop%tsfcl(im) * Sfcprop%landfrac(im) &
                 + Sfcprop%tisfc(im) * tem                      &
                 + Sfcprop%tsfco(im) * (tem1-tem)
          enddo
        enddo
      else
        !$omp parallel do default(shared) private(nb, ix, im, tem)
        do nb = 1, Atm_block%nblks
          do ix = 1, Atm_block%blksz(nb)
            im = Model%chunk_begin(nb)+ix-1
            if (Sfcprop%slmsk(im) == 1) then
              Sfcprop%tsfc(im) = Sfcprop%tsfcl(im)
            else
              tem = one - Sfcprop%fice(im)
              Sfcprop%tsfc(im) = Sfcprop%tisfc(im) * Sfcprop%fice(im) &
                   + Sfcprop%tsfco(im) * tem
            endif
          enddo
        enddo
      endif
    endif compute_tsfc_for_colstart

    if (sfc%var2(i,j,sfc%nvar2m) < -9990.0_kind_phys) then
      if (Model%me == Model%master ) call mpp_error(NOTE, 'gfs_driver::surface_props_input - computing zorlwav')
      !$omp parallel do default(shared) private(nb, ix, im)
      do nb = 1, Atm_block%nblks
        do ix = 1, Atm_block%blksz(nb)
          im = Model%chunk_begin(nb)+ix-1
          Sfcprop%zorlwav(im) = Sfcprop%zorl(im) !--- compute zorlwav from existing variables
        enddo
      enddo
    endif

    if (nint(sfc%var3ice(1,1,1)) == -9999) then    !--- initialize internal ice temp from layer 1 and 2 soil temp
      if (Model%me == Model%master ) call mpp_error(NOTE, 'gfs_driver::surface_props_input - computing tiice')
      do nb = 1, Atm_block%nblks
        do ix = 1, Atm_block%blksz(nb)
          im = Model%chunk_begin(nb)+ix-1
          Sfcprop%tiice(im,1) = max(timin, min(con_tice, Sfcprop%stc(im,1)))
          Sfcprop%tiice(im,2) = max(timin, min(con_tice, Sfcprop%stc(im,2)))
        enddo
      enddo
    endif

    endif ! check on which version of the surface file.

  end subroutine Sfc_io_apply_safeguards

  !>@ destructor for Sfc_io_data_type
  subroutine Sfc_io_final(sfc)
    implicit none
    type(Sfc_io_data_type)             :: sfc

    sfc%nvar2m=0
    sfc%nvar2o=0
    sfc%nvar2l=0
    sfc%nvar3=0
    sfc%nvar2r=0
    sfc%nvar2mp=0
    sfc%nvar3mp=0
    sfc%nvar2m=0
    sfc%nvar_before_lake=0
    sfc%is_lsoil=.false.

    ! This #define reduces code length by a lot
#define IF_ASSOC_DEALLOC_NULL(var) \
    if(associated(sfc%var)) then ; \
      deallocate(sfc%var) ; \
      nullify(sfc%var) ; \
    endif

    IF_ASSOC_DEALLOC_NULL(var2)
    IF_ASSOC_DEALLOC_NULL(var3ice)
    IF_ASSOC_DEALLOC_NULL(var3)
    IF_ASSOC_DEALLOC_NULL(var3sn)
    IF_ASSOC_DEALLOC_NULL(var3eq)
    IF_ASSOC_DEALLOC_NULL(var3zn)
    IF_ASSOC_DEALLOC_NULL(name2)
    IF_ASSOC_DEALLOC_NULL(name3)

#undef IF_ASSOC_DEALLOC_NULL

  end subroutine Sfc_io_final

  !>@ Creates ESMF bundles for 2D fields, for writing surface restart files using the write component (quilt)
  subroutine Sfc_io_bundle_2d_fields(sfc, bundle, grid, Model, outputfile)
    use esmf
    use GFS_typedefs, only: GFS_control_type
    implicit none
    class(Sfc_io_data_type)             :: sfc
    type(ESMF_FieldBundle),intent(inout)        :: bundle
    type(ESMF_Grid),intent(inout)               :: grid
    type(GFS_control_type),          intent(in) :: Model
    character(*), intent(in)                    :: outputfile

    real(kind_phys),dimension(:,:),pointer :: temp_r2d
    integer :: num

    if (.not. associated(sfc%var2)) then
      write(0,*)'ERROR sfc%var2, NOT associated'
      return
    endif
    if (.not. associated(sfc%name2)) then
      write(0,*)'ERROR sfc%name2 NOT associated'
      return
    endif

    do num = 1,sfc%nvar2m
      temp_r2d => sfc%var2(:,:,num)
      call create_2d_field_and_add_to_bundle(temp_r2d, trim(sfc%name2(num)), outputfile, grid, bundle)
    enddo

    if (Model%nstf_name(1) > 0) then
      do num = sfc%nvar2m+1,sfc%nvar2m+sfc%nvar2o
        temp_r2d => sfc%var2(:,:,num)
        call create_2d_field_and_add_to_bundle(temp_r2d, trim(sfc%name2(num)), outputfile, grid, bundle)
      enddo
    endif

    if (Model%lsm == Model%lsm_ruc) then ! sfc%nvar2mp =0
      do num = sfc%nvar2m+sfc%nvar2o+1, sfc%nvar2m+sfc%nvar2o+sfc%nvar2r
        temp_r2d => sfc%var2(:,:,num)
        call create_2d_field_and_add_to_bundle(temp_r2d, trim(sfc%name2(num)), outputfile, grid, bundle)
      enddo
    else if (Model%lsm == Model%lsm_noahmp) then ! sfc%nvar2r =0
      do num = sfc%nvar2m+sfc%nvar2o+1,sfc%nvar2m+sfc%nvar2o+sfc%nvar2mp
        temp_r2d => sfc%var2(:,:,num)
        call create_2d_field_and_add_to_bundle(temp_r2d, trim(sfc%name2(num)), outputfile, grid, bundle)
      enddo
    endif
  end subroutine Sfc_io_bundle_2d_fields

  !>@ Creates ESMF bundles for 3D fields, for writing surface restart files using the write component (quilt)
  subroutine Sfc_io_bundle_3d_fields(sfc, bundle, grid, Model, outputfile)
    use esmf
    use GFS_typedefs, only: GFS_control_type
    implicit none
    class(Sfc_io_data_type)             :: sfc
    type(ESMF_FieldBundle),intent(inout)        :: bundle
    type(ESMF_Grid),intent(inout)               :: grid
    type(GFS_control_type),          intent(in) :: Model
    character(*), intent(in)                    :: outputfile

    real(kind_phys),dimension(:,:,:),pointer :: temp_r3d
    integer :: num, i
    real(kind_phys), dimension(:), allocatable :: zaxis_1, zaxis_2, zaxis_3, zaxis_4

    allocate(zaxis_1(Model%kice))
    zaxis_1 = (/ (i, i=1,Model%kice) /)

    temp_r3d => sfc%var3ice(:,:,:)
    call create_3d_field_and_add_to_bundle(temp_r3d, trim(sfc%name3(0)), "zaxis_1", zaxis_1, trim(outputfile), grid, bundle)

    if(Model%lsm == Model%lsm_ruc) then
      do num = 1,sfc%nvar3
        temp_r3d => sfc%var3(:,:,:,num)
        call create_3d_field_and_add_to_bundle(temp_r3d, trim(sfc%name3(num)), "zaxis_1", zaxis_1, trim(outputfile), grid, bundle)
      enddo
    else
      allocate(zaxis_2(Model%lsoil))
      zaxis_2 = (/ (i, i=1,Model%lsoil) /)
      do num = 1,sfc%nvar3
        temp_r3d => sfc%var3(:,:,:,num)
        call create_3d_field_and_add_to_bundle(temp_r3d, trim(sfc%name3(num)), "zaxis_2", zaxis_2, trim(outputfile), grid, bundle)
      enddo
      deallocate(zaxis_2)
    endif

    if (Model%lsm == Model%lsm_noahmp) then
      allocate(zaxis_3(3))
      zaxis_3 = (/ (i, i=1,3) /)

      do num = sfc%nvar3+1,sfc%nvar3+3
        temp_r3d => sfc%var3sn(:,:,:,num)
        call create_3d_field_and_add_to_bundle(temp_r3d, trim(sfc%name3(num)), "zaxis_3", zaxis_3, trim(outputfile), grid, bundle)
      enddo

      allocate(zaxis_2(Model%lsoil))
      zaxis_2 = (/ (i, i=1,Model%lsoil) /)

      temp_r3d => sfc%var3eq(:,:,:,7)
      call create_3d_field_and_add_to_bundle(temp_r3d, trim(sfc%name3(7)), "zaxis_2", zaxis_2, trim(outputfile), grid, bundle)

      allocate(zaxis_4(7))
      zaxis_4 = (/ (i, i=1,7) /)

      temp_r3d => sfc%var3zn(:,:,:,8)
      call create_3d_field_and_add_to_bundle(temp_r3d, trim(sfc%name3(8)), "zaxis_4", zaxis_4, trim(outputfile), grid, bundle)
    endif ! lsm = lsm_noahmp

    if(allocated(zaxis_1)) deallocate(zaxis_1)
    if(allocated(zaxis_2)) deallocate(zaxis_2)
    if(allocated(zaxis_3)) deallocate(zaxis_3)
    if(allocated(zaxis_4)) deallocate(zaxis_4)

  end subroutine Sfc_io_bundle_3d_fields
end module fv3atm_sfc_io
!> @}
