!> \file fv3atm_clm_lake_io.F90
!! This code reads and writes restart files for the CLM Lake Model. The source code of
!! that model can be found in CCPP. Only the fv3atm_restart_io.F90 should ever access
!! these routines.
!!
!! The CLM Lake Model has its own restart code due to its five alternative vertical
!! levels, which don't match the five found in the other surface fields. For the sake
!! of code simplicity, a dedicated file was a better implementation.

module fv3atm_clm_lake_io
  use GFS_typedefs,       only: GFS_sfcprop_type, GFS_control_type, kind_phys
  use block_control_mod,  only: block_control_type
  use fms2_io_mod,        only: FmsNetcdfDomainFile_t, register_axis, &
       register_restart_field, write_data, &
       register_variable_attribute, register_field, get_dimension_size
  use fv3atm_common_io,   only: create_2d_field_and_add_to_bundle, &
       create_3d_field_and_add_to_bundle, axis_type

  implicit none

  private
  public :: clm_lake_data_type, clm_lake_register_axes, clm_lake_allocate_data, &
       clm_lake_register_fields, clm_lake_deallocate_data, clm_lake_write_axes, &
       clm_lake_copy_from_grid, clm_lake_copy_to_grid, clm_lake_bundle_fields, &
       clm_lake_final, clm_lake_fill_data

  !>\defgroup CLM Lake Model restart public interface
  !>  @{

  !>@ The clm_lake_data_type derived type is a class that stores
  !!  temporary arrays used to read or write CLM Lake model restart
  !!  and axis variables. It can safely be declared and unused, but
  !!  you should only call these routines if the CLM Lake Model was
  !!  (or will be) used by this execution of the FV3. It is the
  !!  responsibility of the caller to ensure the necessary data is in
  !!  Sfc_restart, Sfcprop, and Model.
  type clm_lake_data_type
    ! All 2D variables needed for a restart
    real(kind_phys), pointer, private, dimension(:,:) :: &
         T_snow=>null(), T_ice=>null(), &
         lake_snl2d=>null(), lake_h2osno2d=>null(), lake_tsfc=>null(), clm_lakedepth=>null(), &
         lake_savedtke12d=>null(), lake_sndpth2d=>null(), clm_lake_initialized=>null(), &
         input_lakedepth=>null()

    ! All 3D variables needed for a restart
    real(kind_phys), pointer, private, dimension(:,:,:) :: &
         lake_snow_z3d=>null(), lake_snow_dz3d=>null(), &
         lake_snow_zi3d=>null(), lake_h2osoi_vol3d=>null(), lake_h2osoi_liq3d=>null(), &
         lake_h2osoi_ice3d=>null(), lake_t_soisno3d=>null(), lake_t_lake3d=>null(), &
         lake_icefrac3d=>null()

    ! Axis indices in 1-based array, containing non-1-based indices
    real(kind_phys), pointer, private, dimension(:) :: &
         levlake_clm_lake, levsnowsoil_clm_lake, levsnowsoil1_clm_lake
  contains

    ! register_axes calls registers_axis on Sfc_restart for all required axes
    procedure, public :: register_axes => clm_lake_register_axes

    ! allocate_data allocates all of the pointers in this object
    procedure, public :: allocate_data => clm_lake_allocate_data

    ! register_fields calls register_field on Sfc_restart for all CLM Lake model restart variables
    procedure, public :: register_fields => clm_lake_register_fields

    ! deallocate_data deallocates all pointers, allowing this object to be used repeatedly.
    ! It is safe to call deallocate_data if no data has been allocated.
    procedure, public :: deallocate_data => clm_lake_deallocate_data

    ! write_axes writes variables to Sfc_restart, with the name of
    ! each axis, containing the appropriate information
    procedure, public :: write_axes => clm_lake_write_axes

    ! fills internal arrays with zero:
    procedure, public :: fill_data => clm_lake_fill_data

    ! copy_from_grid copies from Sfcprop to internal pointers (declared above)
    procedure, public :: copy_from_grid => clm_lake_copy_from_grid

    ! copy_from_grid copies from internal pointers (declared above) to Sfcprop
    procedure, public :: copy_to_grid => clm_lake_copy_to_grid

    ! send field bundles in restart quilt server
    procedure, public :: bundle_fields => clm_lake_bundle_fields

    ! A fortran 2003 compliant compiler will call clm_lake_final
    ! automatically when an object of this type goes out of
    ! scope. This will deallocate any arrays via a call to
    ! deallocate_data. It is safe to call this routine if no data has
    ! been allocated.
    final :: clm_lake_final
  end type clm_lake_data_type

CONTAINS

  !>@ This subroutine is clm_lake%alocate_data. It deallocates all
  !!  data, and reallocate to the size specified in Model
  subroutine clm_lake_allocate_data(clm_lake,Model)
    implicit none
    class(clm_lake_data_type) :: clm_lake
    type(GFS_control_type),   intent(in) :: Model

    integer :: nx, ny, i
    call clm_lake%deallocate_data

    nx=Model%nx
    ny=Model%ny

    allocate(clm_lake%T_snow(nx,ny))
    allocate(clm_lake%T_ice(nx,ny))
    allocate(clm_lake%lake_snl2d(nx,ny))
    allocate(clm_lake%lake_h2osno2d(nx,ny))
    allocate(clm_lake%lake_tsfc(nx,ny))
    allocate(clm_lake%lake_savedtke12d(nx,ny))
    allocate(clm_lake%lake_sndpth2d(nx,ny))
    allocate(clm_lake%clm_lakedepth(nx,ny))
    allocate(clm_lake%clm_lake_initialized(nx,ny))
    allocate(clm_lake%input_lakedepth(nx,ny))

    allocate(clm_lake%lake_snow_z3d(nx,ny,Model%nlevsnowsoil1_clm_lake))
    allocate(clm_lake%lake_snow_dz3d(nx,ny,Model%nlevsnowsoil1_clm_lake))
    allocate(clm_lake%lake_snow_zi3d(nx,ny,Model%nlevsnowsoil_clm_lake))
    allocate(clm_lake%lake_h2osoi_vol3d(nx,ny,Model%nlevsnowsoil1_clm_lake))
    allocate(clm_lake%lake_h2osoi_liq3d(nx,ny,Model%nlevsnowsoil1_clm_lake))
    allocate(clm_lake%lake_h2osoi_ice3d(nx,ny,Model%nlevsnowsoil1_clm_lake))
    allocate(clm_lake%lake_t_soisno3d(nx,ny,Model%nlevsnowsoil1_clm_lake))
    allocate(clm_lake%lake_t_lake3d(nx,ny,Model%nlevlake_clm_lake))
    allocate(clm_lake%lake_icefrac3d(nx,ny,Model%nlevlake_clm_lake))

    allocate(clm_lake%levlake_clm_lake(Model%nlevlake_clm_lake))
    allocate(clm_lake%levsnowsoil_clm_lake(Model%nlevsnowsoil_clm_lake))
    allocate(clm_lake%levsnowsoil1_clm_lake(Model%nlevsnowsoil1_clm_lake))

    do i=1,Model%nlevlake_clm_lake
      clm_lake%levlake_clm_lake(i) = i
    enddo
    do i=-Model%nlevsnow_clm_lake,Model%nlevsoil_clm_lake
      clm_lake%levsnowsoil_clm_lake(i+Model%nlevsnow_clm_lake+1) = i
    enddo
    do i=-Model%nlevsnow_clm_lake+1,Model%nlevsoil_clm_lake
      clm_lake%levsnowsoil1_clm_lake(i+Model%nlevsnow_clm_lake) = i
    enddo
  end subroutine clm_lake_allocate_data

  !>@ This is clm_lake%register_axes. It registers all five axes needed
  !!  by CLM Lake restart data.
  subroutine clm_lake_register_axes(clm_lake,Model,Sfc_restart)
    implicit none
    class(clm_lake_data_type) :: clm_lake
    type(GFS_control_type),      intent(in) :: Model
    type(FmsNetcdfDomainFile_t) :: Sfc_restart

    call register_axis(Sfc_restart, 'levlake_clm_lake', dimension_length=Model%nlevlake_clm_lake)
    call register_axis(Sfc_restart, 'levsnowsoil_clm_lake', dimension_length=Model%nlevsnowsoil_clm_lake)
    call register_axis(Sfc_restart, 'levsnowsoil1_clm_lake', dimension_length=Model%nlevsnowsoil1_clm_lake)
  end subroutine clm_lake_register_axes

  !>@ This is clm_lake%write_axes. It creates variables with the name
  !!  name as each clm_lake axis, and fills the variable with the
  !!  appropriate indices
  subroutine clm_lake_write_axes(clm_lake, Model, Sfc_restart)
    implicit none
    class(clm_lake_data_type) :: clm_lake
    type(GFS_control_type),      intent(in) :: Model
    type(FmsNetcdfDomainFile_t) :: Sfc_restart
    integer :: i
    call register_field(Sfc_restart, 'levlake_clm_lake', axis_type, (/'levlake_clm_lake'/))
    call register_variable_attribute(Sfc_restart, 'levlake_clm_lake', 'cartesian_axis' ,'Z', str_len=1)

    call register_field(Sfc_restart, 'levsnowsoil_clm_lake', axis_type, (/'levsnowsoil_clm_lake'/))
    call register_variable_attribute(Sfc_restart, 'levsnowsoil_clm_lake', 'cartesian_axis' ,'Z', str_len=1)

    call register_field(Sfc_restart, 'levsnowsoil1_clm_lake', axis_type, (/'levsnowsoil1_clm_lake'/))
    call register_variable_attribute(Sfc_restart, 'levsnowsoil1_clm_lake', 'cartesian_axis' ,'Z', str_len=1)

    call write_data(Sfc_restart, 'levlake_clm_lake', clm_lake%levlake_clm_lake)
    call write_data(Sfc_restart, 'levsnowsoil_clm_lake', clm_lake%levsnowsoil_clm_lake)
    call write_data(Sfc_restart, 'levsnowsoil1_clm_lake', clm_lake%levsnowsoil1_clm_lake)
  end subroutine clm_lake_write_axes

  !>@ This is clm_lake%fill_data. It fills internal arrays with zero
  !!  Terrible things will happen if you don't call
  !!  clm_lake%allocate_data first.
  subroutine clm_lake_fill_data(clm_lake, Model, Atm_block, Sfcprop)
    implicit none
    class(clm_lake_data_type) :: clm_lake
    type(GFS_sfcprop_type),   intent(in) :: Sfcprop
    type(GFS_control_type),   intent(in) :: Model
    type(block_control_type), intent(in) :: Atm_block

    real(kind_phys), parameter :: zero = 0
    integer :: nb, ix, isc, jsc, i, j
    isc = Model%isc
    jsc = Model%jsc

    ! Copy data to temporary arrays

    !$omp parallel do default(shared) private(i, j, nb, ix)
    do nb = 1, Atm_block%nblks
      do ix = 1, Atm_block%blksz(nb)
        i = Atm_block%index(nb)%ii(ix) - isc + 1
        j = Atm_block%index(nb)%jj(ix) - jsc + 1

        clm_lake%T_snow(i,j) = zero
        clm_lake%T_ice(i,j) = zero
        clm_lake%lake_snl2d(i,j) = zero
        clm_lake%lake_h2osno2d(i,j) = zero
        clm_lake%lake_tsfc(i,j) = zero
        clm_lake%lake_savedtke12d(i,j) = zero
        clm_lake%lake_sndpth2d(i,j) = zero
        clm_lake%clm_lakedepth(i,j) = zero
        clm_lake%clm_lake_initialized(i,j) = zero
        clm_lake%input_lakedepth(i,j) = zero

        clm_lake%lake_snow_z3d(i,j,:) = zero
        clm_lake%lake_snow_dz3d(i,j,:) = zero
        clm_lake%lake_snow_zi3d(i,j,:) = zero
        clm_lake%lake_h2osoi_vol3d(i,j,:) = zero
        clm_lake%lake_h2osoi_liq3d(i,j,:) = zero
        clm_lake%lake_h2osoi_ice3d(i,j,:) = zero
        clm_lake%lake_t_soisno3d(i,j,:) = zero
        clm_lake%lake_t_lake3d(i,j,:) = zero
        clm_lake%lake_icefrac3d(i,j,:) = zero
      enddo
    enddo
  end subroutine clm_lake_fill_data

  !>@ This is clm_lake%copy_from_grid. It copies from Sfcprop
  !!  variables to the corresponding data temporary variables.
  !!  Terrible things will happen if you don't call
  !!  clm_lake%allocate_data first.
  subroutine clm_lake_copy_from_grid(clm_lake, Model, Atm_block, Sfcprop)
    implicit none
    class(clm_lake_data_type) :: clm_lake
    type(GFS_sfcprop_type),   intent(in) :: Sfcprop
    type(GFS_control_type),   intent(in) :: Model
    type(block_control_type), intent(in) :: Atm_block

    integer :: nb, ix, isc, jsc, i, j, im
    isc = Model%isc
    jsc = Model%jsc

    ! Copy data to temporary arrays

    !$omp parallel do default(shared) private(i, j, nb, ix, im)
    do nb = 1, Atm_block%nblks
      do ix = 1, Atm_block%blksz(nb)
        i = Atm_block%index(nb)%ii(ix) - isc + 1
        j = Atm_block%index(nb)%jj(ix) - jsc + 1
        im = Model%chunk_begin(nb)+ix-1

        clm_lake%T_snow(i,j) = Sfcprop%T_snow(im)
        clm_lake%T_ice(i,j) = Sfcprop%T_ice(im)
        clm_lake%lake_snl2d(i,j) = Sfcprop%lake_snl2d(im)
        clm_lake%lake_h2osno2d(i,j) = Sfcprop%lake_h2osno2d(im)
        clm_lake%lake_tsfc(i,j) = Sfcprop%lake_tsfc(im)
        clm_lake%lake_savedtke12d(i,j) = Sfcprop%lake_savedtke12d(im)
        clm_lake%lake_sndpth2d(i,j) = Sfcprop%lake_sndpth2d(im)
        clm_lake%clm_lakedepth(i,j) = Sfcprop%clm_lakedepth(im)
        clm_lake%clm_lake_initialized(i,j) = Sfcprop%clm_lake_initialized(im)
        clm_lake%input_lakedepth(i,j) = Sfcprop%input_lakedepth(im)

        clm_lake%lake_snow_z3d(i,j,:) = Sfcprop%lake_snow_z3d(im,:)
        clm_lake%lake_snow_dz3d(i,j,:) = Sfcprop%lake_snow_dz3d(im,:)
        clm_lake%lake_snow_zi3d(i,j,:) = Sfcprop%lake_snow_zi3d(im,:)
        clm_lake%lake_h2osoi_vol3d(i,j,:) = Sfcprop%lake_h2osoi_vol3d(im,:)
        clm_lake%lake_h2osoi_liq3d(i,j,:) = Sfcprop%lake_h2osoi_liq3d(im,:)
        clm_lake%lake_h2osoi_ice3d(i,j,:) = Sfcprop%lake_h2osoi_ice3d(im,:)
        clm_lake%lake_t_soisno3d(i,j,:) = Sfcprop%lake_t_soisno3d(im,:)
        clm_lake%lake_t_lake3d(i,j,:) = Sfcprop%lake_t_lake3d(im,:)
        clm_lake%lake_icefrac3d(i,j,:) = Sfcprop%lake_icefrac3d(im,:)
      enddo
    enddo
  end subroutine clm_lake_copy_from_grid

  !>@ This is clm_lake%copy_to_grid. It copies from data temporary
  !!  variables to the corresponding Sfcprop variables.  Terrible
  !!  things will happen if you don't call data%allocate_data first.
  subroutine clm_lake_copy_to_grid(clm_lake, Model, Atm_block, Sfcprop)
    implicit none
    class(clm_lake_data_type) :: clm_lake
    type(GFS_sfcprop_type),   intent(in) :: Sfcprop
    type(GFS_control_type),   intent(in) :: Model
    type(block_control_type), intent(in) :: Atm_block

    integer :: nb, ix, isc, jsc, i, j, im
    isc = Model%isc
    jsc = Model%jsc

    ! Copy data to temporary arrays

    !$omp parallel do default(shared) private(i, j, nb, ix, im)
    do nb = 1, Atm_block%nblks
      do ix = 1, Atm_block%blksz(nb)
        i = Atm_block%index(nb)%ii(ix) - isc + 1
        j = Atm_block%index(nb)%jj(ix) - jsc + 1
        im = Model%chunk_begin(nb)+ix-1

        Sfcprop%T_snow(im) = clm_lake%T_snow(i,j)
        Sfcprop%T_ice(im) = clm_lake%T_ice(i,j)
        Sfcprop%lake_snl2d(im) = clm_lake%lake_snl2d(i,j)
        Sfcprop%lake_h2osno2d(im) = clm_lake%lake_h2osno2d(i,j)
        Sfcprop%lake_tsfc(im) = clm_lake%lake_tsfc(i,j)
        Sfcprop%lake_savedtke12d(im) = clm_lake%lake_savedtke12d(i,j)
        Sfcprop%lake_sndpth2d(im) = clm_lake%lake_sndpth2d(i,j)
        Sfcprop%clm_lakedepth(im) = clm_lake%clm_lakedepth(i,j)
        Sfcprop%clm_lake_initialized(im) = clm_lake%clm_lake_initialized(i,j)
        Sfcprop%input_lakedepth(im) = clm_lake%input_lakedepth(i,j)

        Sfcprop%lake_snow_z3d(im,:) = clm_lake%lake_snow_z3d(i,j,:)
        Sfcprop%lake_snow_dz3d(im,:) = clm_lake%lake_snow_dz3d(i,j,:)
        Sfcprop%lake_snow_zi3d(im,:) = clm_lake%lake_snow_zi3d(i,j,:)
        Sfcprop%lake_h2osoi_vol3d(im,:) = clm_lake%lake_h2osoi_vol3d(i,j,:)
        Sfcprop%lake_h2osoi_liq3d(im,:) = clm_lake%lake_h2osoi_liq3d(i,j,:)
        Sfcprop%lake_h2osoi_ice3d(im,:) = clm_lake%lake_h2osoi_ice3d(i,j,:)
        Sfcprop%lake_t_soisno3d(im,:) = clm_lake%lake_t_soisno3d(i,j,:)
        Sfcprop%lake_t_lake3d(im,:) = clm_lake%lake_t_lake3d(i,j,:)
        Sfcprop%lake_icefrac3d(im,:) = clm_lake%lake_icefrac3d(i,j,:)
      enddo
    enddo
  end subroutine clm_lake_copy_to_grid

  !>@ This is clm_lake%register_fields, and it is only used in the
  !!  non-quilt restart. It registers all restart fields needed by the
  !!  CLM Lake Model.  Terrible things will happen if you don't call
  !!  clm_lake%allocate_data and clm_lake%register_axes first.
  subroutine clm_lake_register_fields(clm_lake, Sfc_restart)
    implicit none
    class(clm_lake_data_type) :: clm_lake
    type(FmsNetcdfDomainFile_t) :: Sfc_restart

    integer :: xaxis_1_chunk, yaxis_1_chunk
    integer :: chunksizes2d(3), chunksizes3d(4)

    call get_dimension_size(Sfc_restart, 'xaxis_1', xaxis_1_chunk)
    call get_dimension_size(Sfc_restart, 'yaxis_1', yaxis_1_chunk)
    chunksizes2d = (/xaxis_1_chunk, yaxis_1_chunk, 1/)
    chunksizes3d = (/xaxis_1_chunk, yaxis_1_chunk, 1, 1/)

    ! Register 2D fields
    call register_restart_field(Sfc_restart, 'T_snow', clm_lake%T_snow, &
         dimensions=(/'xaxis_1', 'yaxis_1', 'Time   '/), chunksizes=chunksizes2d, is_optional=.true.)
    call register_restart_field(Sfc_restart, 'T_ice', clm_lake%T_ice, &
         dimensions=(/'xaxis_1', 'yaxis_1', 'Time   '/), chunksizes=chunksizes2d, is_optional=.true.)
    call register_restart_field(Sfc_restart, 'lake_snl2d', clm_lake%lake_snl2d, &
         dimensions=(/'xaxis_1', 'yaxis_1', 'Time   '/), chunksizes=chunksizes2d, is_optional=.true.)
    call register_restart_field(Sfc_restart, 'lake_h2osno2d', clm_lake%lake_h2osno2d, &
         dimensions=(/'xaxis_1', 'yaxis_1', 'Time   '/), chunksizes=chunksizes2d, is_optional=.true.)
    call register_restart_field(Sfc_restart, 'lake_tsfc', clm_lake%lake_tsfc, &
         dimensions=(/'xaxis_1', 'yaxis_1', 'Time   '/), chunksizes=chunksizes2d, is_optional=.true.)
    call register_restart_field(Sfc_restart, 'lake_savedtke12d', clm_lake%lake_savedtke12d, &
         dimensions=(/'xaxis_1', 'yaxis_1', 'Time   '/), chunksizes=chunksizes2d, is_optional=.true.)
    call register_restart_field(Sfc_restart, 'lake_sndpth2d', clm_lake%lake_sndpth2d, &
         dimensions=(/'xaxis_1', 'yaxis_1', 'Time   '/), chunksizes=chunksizes2d, is_optional=.true.)
    call register_restart_field(Sfc_restart, 'clm_lakedepth', clm_lake%clm_lakedepth, &
         dimensions=(/'xaxis_1', 'yaxis_1', 'Time   '/), chunksizes=chunksizes2d, is_optional=.true.)
    call register_restart_field(Sfc_restart, 'clm_lake_initialized', clm_lake%clm_lake_initialized, &
         dimensions=(/'xaxis_1', 'yaxis_1', 'Time   '/), chunksizes=chunksizes2d, is_optional=.true.)
    call register_restart_field(Sfc_restart, 'input_lakedepth', clm_lake%input_lakedepth, &
         dimensions=(/'xaxis_1', 'yaxis_1', 'Time   '/), chunksizes=chunksizes2d, is_optional=.true.)

    ! Register 3D fields
    call register_restart_field(Sfc_restart,'lake_snow_z3d', clm_lake%lake_snow_z3d, &
         dimensions=(/'xaxis_1              ', 'yaxis_1              ', &
         'levsnowsoil1_clm_lake', 'Time                 '/), chunksizes=chunksizes3d, is_optional=.true.)
    call register_restart_field(Sfc_restart,'lake_snow_dz3d', clm_lake%lake_snow_dz3d, &
         dimensions=(/'xaxis_1              ', 'yaxis_1              ', &
         'levsnowsoil1_clm_lake', 'Time                 '/), chunksizes=chunksizes3d, is_optional=.true.)
    call register_restart_field(Sfc_restart,'lake_snow_zi3d', clm_lake%lake_snow_zi3d, &
         dimensions=(/'xaxis_1              ', 'yaxis_1              ', &
         'levsnowsoil_clm_lake ', 'Time                 '/), chunksizes=chunksizes3d, is_optional=.true.)
    call register_restart_field(Sfc_restart,'lake_h2osoi_vol3d', clm_lake%lake_h2osoi_vol3d, &
         dimensions=(/'xaxis_1              ', 'yaxis_1              ', &
         'levsnowsoil1_clm_lake', 'Time                 '/), chunksizes=chunksizes3d, is_optional=.true.)
    call register_restart_field(Sfc_restart,'lake_h2osoi_liq3d', clm_lake%lake_h2osoi_liq3d, &
         dimensions=(/'xaxis_1              ', 'yaxis_1              ', &
         'levsnowsoil1_clm_lake', 'Time                 '/), chunksizes=chunksizes3d, is_optional=.true.)
    call register_restart_field(Sfc_restart,'lake_h2osoi_ice3d', clm_lake%lake_h2osoi_ice3d, &
         dimensions=(/'xaxis_1              ', 'yaxis_1              ', &
         'levsnowsoil1_clm_lake', 'Time                 '/), chunksizes=chunksizes3d, is_optional=.true.)
    call register_restart_field(Sfc_restart,'lake_t_soisno3d', clm_lake%lake_t_soisno3d, &
         dimensions=(/'xaxis_1              ', 'yaxis_1              ', &
         'levsnowsoil1_clm_lake', 'Time                 '/), chunksizes=chunksizes3d, is_optional=.true.)
    call register_restart_field(Sfc_restart,'lake_t_lake3d', clm_lake%lake_t_lake3d, &
         dimensions=(/'xaxis_1              ', 'yaxis_1              ', &
         'levlake_clm_lake     ', 'Time                 '/), chunksizes=chunksizes3d, is_optional=.true.)
    call register_restart_field(Sfc_restart,'lake_icefrac3d', clm_lake%lake_icefrac3d, &
         dimensions=(/'xaxis_1              ', 'yaxis_1              ', &
         'levlake_clm_lake     ', 'Time                 '/), chunksizes=chunksizes3d, is_optional=.true.)
  end subroutine clm_lake_register_fields

  !>@ This is clm_lake%bundle_fields, and it is only used in the
  !!  quilt restart. It bundles all fields needed by the CLM Lake
  !!  Model, which makes them available to ESMF for restart I/O.
  !!  Terrible things will happen if you don't call
  !!  clm_lake%allocate_data and clm_lake%register_axes first.
  subroutine clm_lake_bundle_fields(clm_lake, bundle, grid, Model, outputfile)
    use esmf
    use GFS_typedefs, only: GFS_control_type
    implicit none
    class(Clm_lake_data_type)             :: clm_lake
    type(ESMF_FieldBundle),intent(inout)        :: bundle
    type(ESMF_Grid),intent(inout)               :: grid
    type(GFS_control_type),          intent(in) :: Model
    character(*), intent(in)                    :: outputfile

    real(kind_phys),dimension(:,:),pointer :: temp_r2d
    real(kind_phys),dimension(:,:,:),pointer :: temp_r3d
    integer :: num

    ! Register 2D fields
    call create_2d_field_and_add_to_bundle(clm_lake%T_snow, "T_snow", trim(outputfile), grid, bundle)
    call create_2d_field_and_add_to_bundle(clm_lake%T_ice, 'T_ice', trim(outputfile), grid, bundle)
    call create_2d_field_and_add_to_bundle(clm_lake%lake_snl2d, "lake_snl2d", trim(outputfile), grid, bundle)
    call create_2d_field_and_add_to_bundle(clm_lake%lake_h2osno2d, "lake_h2osno2d", trim(outputfile), grid, bundle)
    call create_2d_field_and_add_to_bundle(clm_lake%lake_tsfc, "lake_tsfc", trim(outputfile), grid, bundle)
    call create_2d_field_and_add_to_bundle(clm_lake%lake_savedtke12d, "lake_savedtke12d", trim(outputfile), grid, bundle)
    call create_2d_field_and_add_to_bundle(clm_lake%lake_sndpth2d, "lake_sndpth2d", trim(outputfile), grid, bundle)
    call create_2d_field_and_add_to_bundle(clm_lake%clm_lakedepth, "clm_lakedepth", trim(outputfile), grid, bundle)
    call create_2d_field_and_add_to_bundle(clm_lake%clm_lake_initialized, "clm_lake_initialized", trim(outputfile), grid, bundle)
    call create_2d_field_and_add_to_bundle(clm_lake%input_lakedepth, "input_lakedepth", trim(outputfile), grid, bundle)

    ! Register 3D fields
    call create_3d_field_and_add_to_bundle(clm_lake%lake_snow_z3d, 'lake_snow_z3d', 'levsnowsoil1_clm_lake', &
         clm_lake%levsnowsoil1_clm_lake, trim(outputfile), grid, bundle)
    call create_3d_field_and_add_to_bundle(clm_lake%lake_snow_dz3d, 'lake_snow_dz3d', 'levsnowsoil1_clm_lake', &
         clm_lake%levsnowsoil1_clm_lake, trim(outputfile), grid, bundle)
    call create_3d_field_and_add_to_bundle(clm_lake%lake_snow_zi3d, 'lake_snow_zi3d', 'levsnowsoil_clm_lake', &
         clm_lake%levsnowsoil_clm_lake, trim(outputfile), grid, bundle)
    call create_3d_field_and_add_to_bundle(clm_lake%lake_h2osoi_vol3d, 'lake_h2osoi_vol3d', 'levsnowsoil1_clm_lake', &
         clm_lake%levsnowsoil1_clm_lake, trim(outputfile), grid, bundle)
    call create_3d_field_and_add_to_bundle(clm_lake%lake_h2osoi_liq3d, 'lake_h2osoi_liq3d', 'levsnowsoil1_clm_lake', &
         clm_lake%levsnowsoil1_clm_lake, trim(outputfile), grid, bundle)
    call create_3d_field_and_add_to_bundle(clm_lake%lake_h2osoi_ice3d, 'lake_h2osoi_ice3d', 'levsnowsoil1_clm_lake', &
         clm_lake%levsnowsoil1_clm_lake, trim(outputfile), grid, bundle)
    call create_3d_field_and_add_to_bundle(clm_lake%lake_t_soisno3d, 'lake_t_soisno3d', 'levsnowsoil1_clm_lake', &
         clm_lake%levsnowsoil1_clm_lake, trim(outputfile), grid, bundle)
    call create_3d_field_and_add_to_bundle(clm_lake%lake_t_lake3d, 'lake_t_lake3d', 'levlake_clm_lake', &
         clm_lake%levlake_clm_lake, trim(outputfile), grid, bundle)
    call create_3d_field_and_add_to_bundle(clm_lake%lake_icefrac3d, 'lake_icefrac3d', 'levlake_clm_lake', &
         clm_lake%levlake_clm_lake, trim(outputfile), grid, bundle)

  end subroutine Clm_lake_bundle_fields

  !>@ Final routine (destructor) for the clm_lake_data_type, called
  !!  automatically when an object of that type goes out of scope.  This
  !!  is simply a wrapper around clm_lake%deallocate_data().
  subroutine clm_lake_final(clm_lake)
    implicit none
    type(clm_lake_data_type) :: clm_lake
    call clm_lake_deallocate_data(clm_lake)
  end subroutine clm_lake_final

  !>@ This is clm_lake%deallocate_data. It deallocates all data used,
  !!  and nullifies the pointers. The clm_lake object can safely be
  !!  used again after this call. This is also the implementation of
  !!  the clm_lake_data_type final routine.
  subroutine clm_lake_deallocate_data(clm_lake)
    implicit none
    class(clm_lake_data_type) :: clm_lake

    ! Deallocate and nullify any associated pointers

    ! This #define reduces code length by a lot
#define IF_ASSOC_DEALLOC_NULL(var) \
    if(associated(clm_lake%var)) then ; \
      deallocate(clm_lake%var) ; \
      nullify(clm_lake%var) ; \
    endif

    IF_ASSOC_DEALLOC_NULL(T_snow)
    IF_ASSOC_DEALLOC_NULL(T_ice)
    IF_ASSOC_DEALLOC_NULL(lake_snl2d)
    IF_ASSOC_DEALLOC_NULL(lake_h2osno2d)
    IF_ASSOC_DEALLOC_NULL(lake_tsfc)
    IF_ASSOC_DEALLOC_NULL(lake_savedtke12d)
    IF_ASSOC_DEALLOC_NULL(lake_sndpth2d)
    IF_ASSOC_DEALLOC_NULL(clm_lakedepth)
    IF_ASSOC_DEALLOC_NULL(clm_lake_initialized)
    IF_ASSOC_DEALLOC_NULL(input_lakedepth)

    IF_ASSOC_DEALLOC_NULL(lake_snow_z3d)
    IF_ASSOC_DEALLOC_NULL(lake_snow_dz3d)
    IF_ASSOC_DEALLOC_NULL(lake_snow_zi3d)
    IF_ASSOC_DEALLOC_NULL(lake_h2osoi_vol3d)
    IF_ASSOC_DEALLOC_NULL(lake_h2osoi_liq3d)
    IF_ASSOC_DEALLOC_NULL(lake_h2osoi_ice3d)
    IF_ASSOC_DEALLOC_NULL(lake_t_soisno3d)
    IF_ASSOC_DEALLOC_NULL(lake_t_lake3d)
    IF_ASSOC_DEALLOC_NULL(lake_icefrac3d)

#undef IF_ASSOC_DEALLOC_NULL
  end subroutine clm_lake_deallocate_data

end module fv3atm_clm_lake_io
!> @}
