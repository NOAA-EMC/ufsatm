! #########################################################################################
!
! MPAS configuration information
!
! #########################################################################################
module module_mpas_config
  use MPAS_typedefs, only: r8 => kind_dbl_prec, r4 => kind_sngl_prec
  use GFS_typedefs, only: pi => con_pi
  use mpi_f08
  use pio, only : iosystem_desc_t, file_desc_t, io_desc_t
  use esmf

  implicit none

  !> Atmosphere time step in seconds
  integer                  :: dt_atmos

  !> Number of MPAS dycore calls per ATMosphere time step.
  integer                  :: n_atmos

  !> MPI communicator for the forecast grid component
  type(MPI_Comm)           :: fcst_mpi_comm

  !> Total number of mpi tasks for the forecast grid components
  integer                  :: fcst_ntasks

  !> The first integration step
  integer                  :: first_kdt

  !> ID number for the coupled grids
  integer                  :: cpl_grid_id

  !> Flag to decide if model writes out coupled diagnostic fields
  logical                  :: cplprint_flag = .false.

  !> Flag to decide if write grid components is used
  logical                  :: quilting = .false.

  !> Flag to decide if write grid component writes out restart files
  logical                  :: quilting_restart = .false.

  !> Output frequency if this array has only two elements and the value of
  !! the second eletment is -1. Otherwise, it is the specific output forecast
  !! hours
  real,dimension(:),allocatable :: output_fh

  !> Calendar type
  character(17)            :: calendar='                 '

  !> MPAS Initial Condition file (via UFSATM NML)
  character(len=256) :: ic_filename

  !> MPAS Lateral Boundary Condition file (via UFSATM NML)
  character(len=256) :: lbc_filename

  !> MPAS output filenames
  character(len=256) :: output_filename = "output.mpas.nc"
  character(len=256) :: restart_filename = "restart.mpas.nc"

  !> PIO
  type(iosystem_desc_t), pointer :: pio_subsystem_ic
  type(iosystem_desc_t), pointer :: pio_subsystem_lbc
  type(iosystem_desc_t), pointer :: pio_subsystem_output
  type(file_desc_t), target :: pioid_ic
  type(file_desc_t), target :: pioid_lbc
  type(file_desc_t), target :: pioid_output
  type(io_desc_t) :: pio_iodesc
  integer :: pio_iotype
  integer :: pio_ioformat
  integer :: pio_stride
  integer :: pio_numiotasks
  logical :: pio_subsystem_output_file_created = .false.
  integer :: pio_subsystem_output_record = 1
  integer, parameter :: TIMELEVEL_NOW = 1 ! current time
  integer, parameter :: TIMELEVEL_NEXT = 2 ! updated/next time

  !> MPAS Grid information
  real(r8), target, allocatable :: zref(:)
  real(r8), target, allocatable :: zref_edge(:)
  real(r8), target, allocatable :: pref(:)
  real(r8), target, allocatable :: pref_edge(:)

  !> sphere_radius is a global attribute in the MPAS initial file.  It is needed to
  !> normalize the cell areas to a unit sphere.
  real(r8) :: sphere_radius

  integer :: maxNCells     ! maximum number of cells for any task (nCellsSolve <= maxNCells)
  integer :: maxEdges      ! maximum number of edges per cell
  integer :: nVertLevels   ! number of vertical layers (midpoints)

  integer, pointer :: &
       nCells,          & ! number of cells in task
       nCellsSolve,     & ! number of cells that a task solves
       nEdgesSolve,     & ! number of edges (velocity) that a task solves
       nVerticesSolve,  & ! number of vertices (vorticity) that a task solves
       nVertLevelsSolve

  real(r4), pointer :: latCell(:), lonCell(:)

  !> Global gridded data
  integer :: nCellsGlobal     ! global number of cells/columns
  integer :: nEdgesGlobal     ! global number of edges
  integer :: nVerticesGlobal  ! global number of vertices

  !> GridCell Longitue/Latitue/Area
  real(r4), allocatable :: latCellGlobal(:)
  real(r4), allocatable :: lonCellGlobal(:)
  real(r4), allocatable :: areaCellGlobal(:)

end module module_mpas_config
