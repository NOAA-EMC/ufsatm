! #########################################################################################
!
! MPAS configuration information
!
! #########################################################################################
module module_mpas_config
  use MPAS_typedefs, only: r8 => kind_dbl_prec
  use mpi_f08
  use pio, only : iosystem_desc_t, file_desc_t, io_desc_t
  use esmf

  implicit none

  !> Atmosphere time step in seconds
  integer                  :: dt_atmos

  !> MPI communicator for the forecast grid component
  type(MPI_Comm)           :: fcst_mpi_comm

  !> Total number of mpi tasks for the forecast grid components
  integer                  :: fcst_ntasks

  !> Output frequency if this array has only two elements and the value of
  !! the second eletment is -1. Otherwise, it is the specific output forecast
  !! hours
  real,dimension(:),allocatable :: output_fh

  !> Calendar type
  character(17)            :: calendar='                 '

  !> Files (Should come from namelist. ToDo)
  !character(len=256) :: mesh_filename = "external mesh file"!"x1.40962.grid.nc" ! This is not actually used during INIT.
  character(len=256) :: ic_filename   = "mpas.init.nc"

  !> PIO
  type(iosystem_desc_t), pointer :: pio_subsystem
  integer :: pio_iotype
  integer :: pio_ioformat
  integer :: pio_stride
  integer :: pio_numiotasks
  type(file_desc_t), target :: pioid
  type(io_desc_t) :: pio_iodesc
  
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

  !> Global gridded data
  integer :: nCells_g      ! global number of cells/columns
  integer :: nEdges_g      ! global number of edges
  integer :: nVertices_g   ! global number of vertices
  
end module module_mpas_config
