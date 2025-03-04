module module_mpas_config

  use mpi_f08
  use pio, only : iosystem_desc_t, file_desc_t
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

  !> PIO
  type(iosystem_desc_t), dimension(1), target, public :: pio_subsystems
  type(file_desc_t),     dimension(1), target, public :: fh_init
  
end module module_mpas_config
