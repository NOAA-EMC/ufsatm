!> @file
!> @brief fv3 configure variables from model_configure.
!> @author Jun Wang @date 01/2017

!> @brief fv3 configure variables from model_configure.
!>
!> @author Jun Wang @date 01/2017
  module module_fv3_config
  use esmf

  implicit none


  !> Atmosphere time step in seconds
  integer                  :: dt_atmos

  !> The first integration step
  integer                  :: first_kdt

  !> MPI communicator for the forecast grid component
  integer                  :: fcst_mpi_comm

  !> Total number of mpi tasks for the forecast grid components
  integer                  :: fcst_ntasks


  !> ID number for the coupled grids
  integer                  :: cpl_grid_id

  !> Flag to decide if model writes out coupled diagnostic fields
  logical                  :: cplprint_flag

  !> Flag to decide if write grid components is used
  logical                  :: quilting

  !> Flag to decide if write grid component writes out restart files
  logical                  :: quilting_restart


  !> Output frequency if this array has only two elements and the value of
  !! the second eletment is -1. Otherwise, it is the specific output forecast
  !! hours
  real,dimension(:),allocatable                   :: output_fh

  !> Calendar type
  character(17)            :: calendar='                 '

  end module module_fv3_config
