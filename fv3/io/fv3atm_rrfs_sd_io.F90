!> \file fv3atm_rrfs_sd_io.F90
!! This file contains derived types and subroutines for RRFS-SD scheme I/O.
!! They read and write restart files, and read emissions data.

module fv3atm_rrfs_sd_io
  use block_control_mod,  only: block_control_type
  use fms2_io_mod,        only: FmsNetcdfDomainFile_t, write_data, &
                                register_axis, register_restart_field, &
                                register_variable_attribute, register_field, &
                                get_dimension_size
  use GFS_typedefs,       only: GFS_sfcprop_type, GFS_control_type, kind_phys
  use fv3atm_common_io,   only: get_nx_ny_from_atm, create_2d_field_and_add_to_bundle, &
                                create_3d_field_and_add_to_bundle, axis_type

  implicit none

  private

  public :: rrfs_sd_state_type, rrfs_sd_state_register_axis, rrfs_sd_state_write_axis, &
       rrfs_sd_state_fill_data, rrfs_sd_state_register_fields, rrfs_sd_state_deallocate_data, &
       rrfs_sd_state_copy_from_grid, rrfs_sd_state_copy_to_grid, &
       rrfs_sd_state_final

  public :: rrfs_sd_emissions_type, rrfs_sd_emissions_final, &
       rrfs_sd_emissions_register_dust12m, rrfs_sd_emissions_copy_dust12m, &
       rrfs_sd_emissions_register_emi, rrfs_sd_emissions_copy_emi, &
       rrfs_sd_emissions_register_fire, rrfs_sd_emissions_copy_fire

  !>\defgroup fv3atm_rrfs_sd_io module
  !> @{

  !>@ Temporary data storage for reading and writing restart data for the RRFS-SD scheme.
  type rrfs_sd_state_type
    ! The rrfs_sd_state_type stores temporary arrays used to read or
    ! write RRFS-SD restart and axis variables.

    real(kind_phys), pointer, private, dimension(:,:) :: & ! i,j variables
         emdust=>null(), emseas=>null(), emanoc=>null(), fhist=>null(), coef_bb_dc=>null()

    real(kind_phys), pointer, private, dimension(:,:,:) :: &
         fire_in=>null() ! i, j, fire_aux_data_levels

    real(kind_phys), pointer, private, dimension(:) :: &
         fire_aux_data_levels=>null() ! 1:Model%fire_aux_data_levels index array for metadata write

  contains
    procedure, public :: register_axis => rrfs_sd_state_register_axis ! register fire_aux_data_levels axis
    procedure, public :: write_axis => rrfs_sd_state_write_axis ! write fire_aux_data_levels variable
    procedure, public :: allocate_data => rrfs_sd_state_allocate_data ! allocate all pointers
    procedure, public :: fill_data => rrfs_sd_state_fill_data ! fill data with default values
    procedure, public :: register_fields => rrfs_sd_state_register_fields ! register rrfs_sd fields
    procedure, public :: deallocate_data => rrfs_sd_state_deallocate_data ! deallocate pointers
    procedure, public :: copy_from_grid => rrfs_sd_state_copy_from_grid ! Copy Sfcprop to arrays
    procedure, public :: copy_to_grid => rrfs_sd_state_copy_to_grid ! Copy arrays to Sfcprop
    procedure, public :: bundle_fields => rrfs_sd_bundle_fields ! Point esmf bundles to arrays
    final :: rrfs_sd_state_final ! Destructor; calls deallocate_data
  end type rrfs_sd_state_type

  ! --------------------------------------------------------------------

  !>@ Temporary data storage for reading RRFS-SD emissions data
  type rrfs_sd_emissions_type
    integer, private :: nvar_dust12m = 5
    integer, private :: nvar_emi = 1
    integer, private :: nvar_fire = 2
    integer, private :: nvar_fire2d = 5

    character(len=32), pointer, dimension(:), private :: dust12m_name => null()
    character(len=32), pointer, dimension(:), private :: emi_name => null()
    character(len=32), pointer, dimension(:), private :: fire_name => null()
    character(len=32), pointer, dimension(:), private :: fire_name2d => null()

    real(kind=kind_phys), pointer, dimension(:,:,:,:), private :: dust12m_var => null()
    real(kind=kind_phys), pointer, dimension(:,:,:,:), private :: emi_var => null()
    real(kind=kind_phys), pointer, dimension(:,:,:,:), private :: fire_var => null()
    real(kind=kind_phys), pointer, dimension(:,:,:  ), private :: fire_var2d => null()

  contains

    procedure, public :: register_dust12m => rrfs_sd_emissions_register_dust12m
    procedure, public :: copy_dust12m => rrfs_sd_emissions_copy_dust12m

    procedure, public :: register_emi => rrfs_sd_emissions_register_emi
    procedure, public :: copy_emi => rrfs_sd_emissions_copy_emi

    procedure, public :: register_fire => rrfs_sd_emissions_register_fire
    procedure, public :: copy_fire => rrfs_sd_emissions_copy_fire

    final :: rrfs_sd_emissions_final
  end type rrfs_sd_emissions_type

  ! --------------------------------------------------------------------

contains


  ! --------------------------------------------------------------------
  ! -- RRFS_SD_STATE IMPLEMENTATION ------------------------------------
  ! --------------------------------------------------------------------

  !>@ Registers the fire_aux_data_levels axis for restart I/O
  subroutine rrfs_sd_state_register_axis(data,Model,Sfc_restart)
    implicit none
    class(rrfs_sd_state_type) :: data
    type(FmsNetcdfDomainFile_t) :: Sfc_restart
    type(GFS_control_type),      intent(in) :: Model
    call register_axis(Sfc_restart, 'fire_aux_data_levels', &
         dimension_length=Model%fire_aux_data_levels)
  end subroutine rrfs_sd_state_register_axis

  ! --------------------------------------------------------------------

  !>@ Registers and writes the axis indices for the fire_aux_data_levels axis
  subroutine rrfs_sd_state_write_axis(data,Model,Sfc_restart)
    implicit none
    class(rrfs_sd_state_type) :: data
    type(FmsNetcdfDomainFile_t) :: Sfc_restart
    type(GFS_control_type),      intent(in) :: Model

    call register_field(Sfc_restart, 'fire_aux_data_levels', axis_type, (/'fire_aux_data_levels'/))
    call register_variable_attribute(Sfc_restart, 'fire_aux_data_levels', 'cartesian_axis' ,'Z', str_len=1)
    call write_data(Sfc_restart, 'fire_aux_data_levels', data%fire_aux_data_levels)
  end subroutine rrfs_sd_state_write_axis

  ! --------------------------------------------------------------------

  !>@ Allocates temporary arrays for RRFS-SD scheme I/O and stores fire_aux_data_levels axis indices
  subroutine rrfs_sd_state_allocate_data(data,Model)
    implicit none
    class(rrfs_sd_state_type) :: data
    type(GFS_control_type),   intent(in) :: Model
    integer :: nx, ny, i

    call data%deallocate_data

    nx=Model%nx
    ny=Model%ny

    allocate(data%emdust(nx,ny))
    allocate(data%emseas(nx,ny))
    allocate(data%emanoc(nx,ny))
    allocate(data%fhist(nx,ny))
    allocate(data%coef_bb_dc(nx,ny))
    allocate(data%fire_aux_data_levels(Model%fire_aux_data_levels))
    allocate(data%fire_in(nx,ny,Model%fire_aux_data_levels))

    do i=1,Model%fire_aux_data_levels
      data%fire_aux_data_levels(i) = i
    enddo

  end subroutine rrfs_sd_state_allocate_data

  ! --------------------------------------------------------------------

  !>@brief Fills RRFS-SD temporary arrays with reasonable defaults.
  !> \section rrfs_sd_state_type%fill_data() procedure
  !! Fills all temporary variables with default values.
  !! Terrible things will happen if you don't call data%allocate_data first.
  subroutine rrfs_sd_state_fill_data(data, Model, Atm_block, Sfcprop)
    implicit none
    class(rrfs_sd_state_type) :: data
    type(GFS_sfcprop_type),   intent(in) :: Sfcprop
    type(GFS_control_type),   intent(in) :: Model
    type(block_control_type), intent(in) :: Atm_block

    integer :: nb, ix, isc, jsc, i, j

    isc = Model%isc
    jsc = Model%jsc

    !$omp parallel do default(shared) private(i, j, nb, ix)
    do nb = 1, Atm_block%nblks
      do ix = 1, Atm_block%blksz(nb)
        i = Atm_block%index(nb)%ii(ix) - isc + 1
        j = Atm_block%index(nb)%jj(ix) - jsc + 1

        data%emdust(i,j) = 0
        data%emseas(i,j) = 0
        data%emanoc(i,j) = 0
        data%fhist(i,j) = 1.
        data%coef_bb_dc(i,j) = 0

        data%fire_in(i,j,:) = 0
      end do
    end do
  end subroutine rrfs_sd_state_fill_data

  ! --------------------------------------------------------------------

  !>@brief Registers RRFS-SD restart variables (for read or write)
  !> \section rrfs_sd_state_type%register_fields() procedure
  !! Registers all restart fields needed by the RRFS-SD
  !! Terrible things will happen if you don't call data%allocate_data
  !! and data%register_axes first.
  subroutine rrfs_sd_state_register_fields(data,Sfc_restart)
    implicit none
    class(rrfs_sd_state_type) :: data
    type(FmsNetcdfDomainFile_t) :: Sfc_restart

    integer :: xaxis_1_chunk, yaxis_1_chunk
    integer :: chunksizes2d(3), chunksizes3d(4)

    call get_dimension_size(Sfc_restart, 'xaxis_1', xaxis_1_chunk)
    call get_dimension_size(Sfc_restart, 'yaxis_1', yaxis_1_chunk)
    chunksizes2d = (/xaxis_1_chunk, yaxis_1_chunk, 1/)
    chunksizes3d = (/xaxis_1_chunk, yaxis_1_chunk, 1, 1/)

    ! Register 2D fields
    call register_restart_field(Sfc_restart, 'emdust', data%emdust, &
         dimensions=(/'xaxis_1', 'yaxis_1', 'Time   '/), chunksizes=chunksizes2d, is_optional=.true.)
    call register_restart_field(Sfc_restart, 'emseas', data%emseas, &
         dimensions=(/'xaxis_1', 'yaxis_1', 'Time   '/), chunksizes=chunksizes2d, is_optional=.true.)
    call register_restart_field(Sfc_restart, 'emanoc', data%emanoc, &
         dimensions=(/'xaxis_1', 'yaxis_1', 'Time   '/), chunksizes=chunksizes2d, is_optional=.true.)
    call register_restart_field(Sfc_restart, 'fhist', data%fhist, &
         dimensions=(/'xaxis_1', 'yaxis_1', 'Time   '/), chunksizes=chunksizes2d, is_optional=.true.)
    call register_restart_field(Sfc_restart, 'coef_bb_dc', data%coef_bb_dc, &
         dimensions=(/'xaxis_1', 'yaxis_1', 'Time   '/), chunksizes=chunksizes2d, is_optional=.true.)

    ! Register 3D field
    call register_restart_field(Sfc_restart, 'fire_in', data%fire_in, &
         dimensions=(/'xaxis_1             ', 'yaxis_1             ', &
         'fire_aux_data_levels', 'Time                '/), &
         chunksizes=chunksizes3d, is_optional=.true.)
  end subroutine rrfs_sd_state_register_fields

  ! --------------------------------------------------------------------

  !>@brief Creates ESMF bundles for writing RRFS-SD restarts via the write component (quilt)
  !> \section rrfs_sd_state_type%bundle_fields() procedure
  !! Registers all restart fields needed by the RRFS-SD
  !! Terrible things will happen if you don't call data%allocate_data
  !! and data%register_axes first.
  subroutine rrfs_sd_bundle_fields(data, bundle, grid, Model, outputfile)
    use esmf
    use GFS_typedefs, only: GFS_control_type
    implicit none
    class(rrfs_sd_state_type) :: data
    type(ESMF_FieldBundle),intent(inout)        :: bundle
    type(ESMF_Grid),intent(inout)               :: grid
    type(GFS_control_type),          intent(in) :: Model
    character(*), intent(in)                    :: outputfile

    ! Register 2D fields
    call create_2d_field_and_add_to_bundle(data%emdust, "emdust", trim(outputfile), grid, bundle)
    call create_2d_field_and_add_to_bundle(data%emseas, "emseas", trim(outputfile), grid, bundle)
    call create_2d_field_and_add_to_bundle(data%emanoc, "emanoc", trim(outputfile), grid, bundle)
    call create_2d_field_and_add_to_bundle(data%fhist, "fhist", trim(outputfile), grid, bundle)
    call create_2d_field_and_add_to_bundle(data%coef_bb_dc, "coef_bb_dc", trim(outputfile), grid, bundle)

    ! Register 3D field
    call create_3d_field_and_add_to_bundle(data%fire_in, 'fire_in', 'fire_aux_data_levels', &
         data%fire_aux_data_levels, trim(outputfile), grid, bundle)
  end subroutine rrfs_sd_bundle_fields

  ! --------------------------------------------------------------------

  !>@brief Destructor for the rrfs_sd_state_type
  !> \section rrfs_sd_state_type destructor() procedure
  !! Final routine for rrfs_sd_state_type, called automatically when
  !! an object of that type goes out of scope.  This is a wrapper
  !! around data%deallocate_data() with necessary syntactic
  !! differences.
  subroutine rrfs_sd_state_final(data)
    implicit none
    type(rrfs_sd_state_type) :: data
    call rrfs_sd_state_deallocate_data(data)
  end subroutine rrfs_sd_state_final

  ! --------------------------------------------------------------------

  !>@brief Deallocates internal arrays in an rrfs_sd_state_type
  !> \section rrfs_sd_state_type%deallocate_data() procedure
  !! Deallocates all data used, and nullifies the pointers. The data
  !! object can safely be used again after this call. This is also
  !! the implementation of the rrfs_sd_state_deallocate_data final routine.
  subroutine rrfs_sd_state_deallocate_data(data)
    implicit none
    class(rrfs_sd_state_type) :: data

    ! This #define reduces code length by a lot
#define IF_ASSOC_DEALLOC_NULL(var) \
    if(associated(data%var)) then ; \
      deallocate(data%var) ; \
      nullify(data%var) ; \
    endif

    IF_ASSOC_DEALLOC_NULL(emdust)
    IF_ASSOC_DEALLOC_NULL(emseas)
    IF_ASSOC_DEALLOC_NULL(emanoc)
    IF_ASSOC_DEALLOC_NULL(fhist)
    IF_ASSOC_DEALLOC_NULL(coef_bb_dc)

    IF_ASSOC_DEALLOC_NULL(fire_in)

    ! Undefine this to avoid cluttering the cpp scope:
#undef IF_ASSOC_DEALLOC_NULL
  end subroutine rrfs_sd_state_deallocate_data

  ! --------------------------------------------------------------------

  !>@brief Copies from rrfs_sd_state_type internal arrays to the model grid.
  !> \section rrfs_sd_state_type%copy_to_grid() procedure
  !! This procedure is called after reading a restart, to copy restart data
  !! from the rrfs_sd_state_type to the model grid.
  subroutine rrfs_sd_state_copy_to_grid(data, Model, Atm_block, Sfcprop)
    implicit none
    class(rrfs_sd_state_type) :: data
    type(GFS_sfcprop_type),   intent(in) :: Sfcprop
    type(GFS_control_type),   intent(in) :: Model
    type(block_control_type), intent(in) :: Atm_block

    integer :: nb, ix, i, j, im

    !$omp parallel do default(shared) private(i, j, nb, ix, im)
    do nb = 1, Atm_block%nblks
      do ix = 1, Atm_block%blksz(nb)
        i = Atm_block%index(nb)%ii(ix) - Atm_block%isc + 1
        j = Atm_block%index(nb)%jj(ix) - Atm_block%jsc + 1
        im = Model%chunk_begin(nb)+ix-1

        Sfcprop%emdust(im) = data%emdust(i,j)
        Sfcprop%emseas(im) = data%emseas(i,j)
        Sfcprop%emanoc(im) = data%emanoc(i,j)
        Sfcprop%fhist(im) = data%fhist(i,j)
        Sfcprop%coef_bb_dc(im) = data%coef_bb_dc(i,j)

        Sfcprop%fire_in(im,:) = data%fire_in(i,j,:)
      enddo
    enddo
  end subroutine rrfs_sd_state_copy_to_grid

  ! --------------------------------------------------------------------

  !>@brief Copies from the model grid to rrfs_sd_state_type internal arrays
  !> \section rrfs_sd_state_type%copy_from_grid() procedure
  !! This procedure is called before writing the restart, to copy data from
  !! the model grid to rrfs_sd_state_type internal arrays. The ESMF or FMS
  !! restart code will write data from those arrays, not the model grid.
  subroutine rrfs_sd_state_copy_from_grid(data, Model, Atm_block, Sfcprop)
    implicit none
    class(rrfs_sd_state_type) :: data
    type(GFS_sfcprop_type),   intent(in) :: Sfcprop
    type(GFS_control_type),   intent(in) :: Model
    type(block_control_type), intent(in) :: Atm_block

    integer :: nb, ix, i, j, im

    !$omp parallel do default(shared) private(i, j, nb, ix, im)
    do nb = 1, Atm_block%nblks
      do ix = 1, Atm_block%blksz(nb)
        i = Atm_block%index(nb)%ii(ix) - Atm_block%isc + 1
        j = Atm_block%index(nb)%jj(ix) - Atm_block%jsc + 1
        im = Model%chunk_begin(nb)+ix-1

        data%emdust(i,j) = Sfcprop%emdust(im)
        data%emseas(i,j) = Sfcprop%emseas(im)
        data%emanoc(i,j) = Sfcprop%emanoc(im)
        data%fhist(i,j) = Sfcprop%fhist(im)
        data%coef_bb_dc(i,j) = Sfcprop%coef_bb_dc(im)

        data%fire_in(i,j,:) = Sfcprop%fire_in(im,:)
      enddo
    enddo
  end subroutine rrfs_sd_state_copy_from_grid

  ! --------------------------------------------------------------------
  ! -- RRFS_SD_EMISSIONS IMPLEMENTATION --------------------------------
  ! --------------------------------------------------------------------

  !>@ Allocates temporary arrays and registers variables for reading the dust12m file.
  subroutine rrfs_sd_emissions_register_dust12m(data, restart, Atm_block)
    implicit none
    class(rrfs_sd_emissions_type) :: data
    type(FmsNetcdfDomainFile_t) :: restart
    type(block_control_type), intent(in) :: Atm_block

    real(kind=kind_phys), pointer, dimension(:,:,:) :: var3_p2 => NULL()
    integer :: num, nx, ny

    if(associated(data%dust12m_name)) then
      deallocate(data%dust12m_name)
      nullify(data%dust12m_name)
    endif
    if(associated(data%dust12m_var)) then
      deallocate(data%dust12m_var)
      nullify(data%dust12m_var)
    endif

    call get_nx_ny_from_atm(Atm_block, nx, ny)
    allocate(data%dust12m_name(data%nvar_dust12m))
    allocate(data%dust12m_var(nx,ny,12,data%nvar_dust12m))

    data%dust12m_name(1)  = 'clay'
    data%dust12m_name(2)  = 'rdrag'
    data%dust12m_name(3)  = 'sand'
    data%dust12m_name(4)  = 'ssm'
    data%dust12m_name(5)  = 'uthr'

    !--- register axis
    call register_axis(restart, 'lon', 'X')
    call register_axis(restart, 'lat', 'Y')
    call register_axis(restart, 'time', 12)
    !--- register the 3D fields
    do num = 1,data%nvar_dust12m
      var3_p2 => data%dust12m_var(:,:,:,num)
      call register_restart_field(restart, data%dust12m_name(num), var3_p2, &
           dimensions=(/'time', 'lat ', 'lon '/),&
           &is_optional=.true.)
      ! That was "is_optional=.not.mand" in the original, but mand was never initialized.
    enddo
  end subroutine rrfs_sd_emissions_register_dust12m

  ! --------------------------------------------------------------------

  !>@ Called after register_dust12m() to copy data from internal arrays to the model grid and deallocate arrays
  subroutine rrfs_sd_emissions_copy_dust12m(data, Model, Sfcprop, Atm_block)
    implicit none
    type(GFS_control_type),    intent(in) :: Model
    type(GFS_sfcprop_type),    intent(inout) :: Sfcprop
    class(rrfs_sd_emissions_type) :: data
    type(block_control_type), intent(in) :: Atm_block

    integer :: num, nb, i, j, ix, k, im

    if(.not.associated(data%dust12m_name) .or. .not.associated(data%dust12m_var)) then
      write(0,*) 'ERROR: Called copy_dust12m before register_dust12m'
      return
    endif

    !$omp parallel do default(shared) private(i, j, nb, ix, k, im)
    do nb = 1, Atm_block%nblks
      !--- 3D variables
      do ix = 1, Atm_block%blksz(nb)
        i = Atm_block%index(nb)%ii(ix) - Atm_block%isc + 1
        j = Atm_block%index(nb)%jj(ix) - Atm_block%jsc + 1
        im = Model%chunk_begin(nb)+ix-1
        do k = 1, 12
          Sfcprop%dust12m_in(im,k,1)  = data%dust12m_var(i,j,k,1)
          Sfcprop%dust12m_in(im,k,2)  = data%dust12m_var(i,j,k,2)
          Sfcprop%dust12m_in(im,k,3)  = data%dust12m_var(i,j,k,3)
          Sfcprop%dust12m_in(im,k,4)  = data%dust12m_var(i,j,k,4)
          Sfcprop%dust12m_in(im,k,5)  = data%dust12m_var(i,j,k,5)
        enddo
      enddo
    enddo

    deallocate(data%dust12m_name)
    nullify(data%dust12m_name)
    deallocate(data%dust12m_var)
    nullify(data%dust12m_var)
  end subroutine rrfs_sd_emissions_copy_dust12m

  ! --------------------------------------------------------------------

  !>@ Allocates temporary arrays and registers variables for reading the emissions file.
  subroutine rrfs_sd_emissions_register_emi(data, restart, Atm_block)
    implicit none
    class(rrfs_sd_emissions_type) :: data
    type(FmsNetcdfDomainFile_t) :: restart
    type(block_control_type), intent(in) :: Atm_block

    real(kind=kind_phys), pointer, dimension(:,:,:) :: var3_p2 => NULL()
    integer :: num, nx, ny

    if(associated(data%emi_name)) then
      deallocate(data%emi_name)
      nullify(data%emi_name)
    endif

    if(associated(data%emi_var)) then
      deallocate(data%emi_var)
      nullify(data%emi_var)
    endif

    call get_nx_ny_from_atm(Atm_block, nx, ny)
    allocate(data%emi_name(data%nvar_emi))
    allocate(data%emi_var(nx,ny,1,data%nvar_emi))

    data%emi_name(1)  = 'e_oc'
    !--- register axis
    call register_axis( restart, 'time', 1) ! only read first time level, even if multiple are present
    call register_axis( restart, "grid_xt", 'X' )
    call register_axis( restart, "grid_yt", 'Y' )
    !--- register the 2D fields
    do num = 1,data%nvar_emi
      var3_p2 => data%emi_var(:,:,:,num)
      call register_restart_field(restart, data%emi_name(num), var3_p2, &
           dimensions=(/'time   ','grid_yt','grid_xt'/))
    enddo
  end subroutine rrfs_sd_emissions_register_emi

  ! --------------------------------------------------------------------

  !>@ Called after register_emi() to copy data from internal arrays to the model grid and deallocate arrays
  subroutine rrfs_sd_emissions_copy_emi(data, Model, Sfcprop, Atm_block)
    implicit none
    type(GFS_control_type),    intent(in) :: Model
    type(GFS_sfcprop_type),    intent(inout) :: Sfcprop
    class(rrfs_sd_emissions_type) :: data
    type(block_control_type), intent(in) :: Atm_block

    integer :: num, nb, i, j, ix, im

    if(.not.associated(data%emi_name) .or. .not.associated(data%emi_var)) then
      write(0,*) 'ERROR: Called copy_emi before register_emi'
      return
    endif

    do num=1,data%nvar_emi
      !$omp parallel do default(shared) private(i, j, nb, ix, im)
      do nb = 1, Atm_block%nblks
        !--- 2D variables
        do ix = 1, Atm_block%blksz(nb)
          i = Atm_block%index(nb)%ii(ix) - Atm_block%isc + 1
          j = Atm_block%index(nb)%jj(ix) - Atm_block%jsc + 1
          im = Model%chunk_begin(nb)+ix-1
          Sfcprop%emi_in(im,num)  = data%emi_var(i,j,1,num)
        enddo
      enddo
    enddo

    deallocate(data%emi_name)
    nullify(data%emi_name)
    deallocate(data%emi_var)
    nullify(data%emi_var)
  end subroutine rrfs_sd_emissions_copy_emi

  ! --------------------------------------------------------------------

  !>@ Allocates temporary arrays and registers variables for reading the fire data file.
  subroutine rrfs_sd_emissions_register_fire(data, Model, restart, Atm_block)
    implicit none
    class(rrfs_sd_emissions_type) :: data
    type(GFS_control_type),   intent(in) :: Model
    type(FmsNetcdfDomainFile_t) :: restart
    type(block_control_type), intent(in) :: Atm_block

    real(kind=kind_phys), pointer, dimension(:,:) :: var_p2 => NULL()
    real(kind=kind_phys), pointer, dimension(:,:,:) :: var3_p2 => NULL()
    integer :: num, nx, ny
    integer :: ebb_dcycle

    ebb_dcycle=Model%ebb_dcycle

    if(associated(data%fire_name)) then
      deallocate(data%fire_name)
      nullify(data%fire_name)
    endif

    if(associated(data%fire_name2d)) then
      deallocate(data%fire_name2d)
      nullify(data%fire_name2d)
    endif

    if(associated(data%fire_var)) then
      deallocate(data%fire_var)
      nullify(data%fire_var)
    endif

    if(associated(data%fire_var2d)) then
      deallocate(data%fire_var2d)
      nullify(data%fire_var2d)
    endif

    !--- allocate the various containers needed for rrfssd fire data
    call get_nx_ny_from_atm(Atm_block, nx, ny)
    allocate(data%fire_name(data%nvar_fire))
    allocate(data%fire_name2d(data%nvar_fire2d))
    allocate(data%fire_var(nx,ny,24,data%nvar_fire))
    allocate(data%fire_var2d(nx,ny,data%nvar_fire2d))

    data%fire_name(1)  = 'ebb_smoke_hr'  ! 2d x 24 hours
    data%fire_name(2)  = 'frp_avg_hr'    ! 2d x 24 hours

    ! For the operational system
    data%fire_name2d(1)  = 'ebb_rate'  ! 2d
    data%fire_name2d(2)  = 'frp_davg'
    data%fire_name2d(3)  = 'fire_end_hr'
    data%fire_name2d(4)  = 'hwp_davg'
    data%fire_name2d(5)  = 'totprcp_24hrs'

    !--- register axis
    call register_axis(restart, 'lon', 'X')
    call register_axis(restart, 'lat', 'Y')
    if (ebb_dcycle==1) then ! -- retro mode
     !--- register the 3D fields
     call register_axis(restart, 't', 24)
     do num = 1,data%nvar_fire
      var3_p2 => data%fire_var(:,:,:,num)
      call register_restart_field(restart, data%fire_name(num), var3_p2, &
           dimensions=(/'t  ', 'lat', 'lon'/), is_optional=.true.)
     enddo
    elseif (ebb_dcycle==2) then ! -- forecast mode
     !--- register the 2D fields
     call register_axis(restart, 't', 1)
     do num = 1,data%nvar_fire2d
      var_p2 => data%fire_var2d(:,:,num)
      call register_restart_field(restart, data%fire_name2d(num), var_p2, &
           dimensions=(/'lat', 'lon'/), is_optional=.true.)
     enddo
    else
     ! -- user define their own fire emission
    endif

  end subroutine rrfs_sd_emissions_register_fire

  ! --------------------------------------------------------------------

  !>@ Called after register_fire() to copy data from internal arrays to the model grid and deallocate arrays
  subroutine rrfs_sd_emissions_copy_fire(data, Model, Sfcprop, Atm_block)
    implicit none
    class(rrfs_sd_emissions_type) :: data
    type(GFS_control_type),   intent(in) :: Model
    type(GFS_sfcprop_type),    intent(inout) :: Sfcprop
    type(block_control_type), intent(in) :: Atm_block

    integer :: nb, ix, k, i, j, im
    integer :: ebb_dcycle

    ebb_dcycle=Model%ebb_dcycle

    !$omp parallel do default(shared) private(i, j, nb, ix, k, im)
    do nb = 1, Atm_block%nblks
      do ix = 1, Atm_block%blksz(nb)
        i = Atm_block%index(nb)%ii(ix) - Atm_block%isc + 1
        j = Atm_block%index(nb)%jj(ix) - Atm_block%jsc + 1
        im = Model%chunk_begin(nb)+ix-1
        if (ebb_dcycle==1) then ! -- retro mode
        !--- 3D variables
         do k = 1, 24
          Sfcprop%smoke_RRFS(im,k,1)  = data%fire_var(i,j,k,1)
          Sfcprop%smoke_RRFS(im,k,2)  = data%fire_var(i,j,k,2)
         enddo
        elseif (ebb_dcycle==2) then ! -- forecast mode
        !--- 2D variables
          Sfcprop%smoke2d_RRFS(im,1)  = data%fire_var2d(i,j,1)
          Sfcprop%smoke2d_RRFS(im,2)  = data%fire_var2d(i,j,2)
          Sfcprop%smoke2d_RRFS(im,3)  = data%fire_var2d(i,j,3)
          Sfcprop%smoke2d_RRFS(im,4)  = data%fire_var2d(i,j,4)
          Sfcprop%smoke2d_RRFS(im,5)  = data%fire_var2d(i,j,5)
        else
         ! -- user define their own fire emission
        endif
      enddo
    enddo
  end subroutine rrfs_sd_emissions_copy_fire

  !>@ Destructor for rrfs_sd_emissions_type
  subroutine rrfs_sd_emissions_final(data)
    implicit none
    type(rrfs_sd_emissions_type) :: data

    ! This #define reduces code length by a lot
#define IF_ASSOC_DEALLOC_NULL(var) \
    if(associated(data%var)) then ; \
      deallocate(data%var) ; \
      nullify(data%var) ; \
    endif

    IF_ASSOC_DEALLOC_NULL(dust12m_name)
    IF_ASSOC_DEALLOC_NULL(emi_name)
    IF_ASSOC_DEALLOC_NULL(fire_name)
    IF_ASSOC_DEALLOC_NULL(dust12m_var)
    IF_ASSOC_DEALLOC_NULL(emi_var)
    IF_ASSOC_DEALLOC_NULL(fire_var)

    ! Undefine this to avoid cluttering the cpp scope:
#undef IF_ASSOC_DEALLOC_NULL
  end subroutine rrfs_sd_emissions_final

end module fv3atm_rrfs_sd_io

!> @}
