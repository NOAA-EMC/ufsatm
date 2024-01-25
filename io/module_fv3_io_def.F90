!> @file
!> @brief fv3 I/O related configration variables.
!> @author Jun Wang @date 01/2017

!> @brief fv3 I/O related configration variables.
!>
!> @author Jun Wang @date 01/2017
module module_fv3_io_def
  use esmf, only     : esmf_maxstr
  implicit none

  !> Number of processors used in the forecast grid component
  integer           :: num_pes_fcst

  !> Number of write tasks per write group.  
  integer           :: wrttasks_per_group

  !> Number of the write groups
  integer           :: write_groups 

  !> Current write group 
  integer           :: n_group 

  !> Number of history files
  integer           :: num_files 

  !> Number of the ESMF field bundles for physics fields
  integer           :: nbdlphys 

  !> IAU running window length
  integer           :: iau_offset 

  !> Logical variable to decide if full time (HH.MM.SS) is used in the history
  !! file names
  logical           :: lflname_fulltime 

  !> Logical variable to decide if unlimited time dimension is used
  logical           :: time_unlimited 


  !> Base names for model history output files
  character(len=esmf_maxstr),dimension(:),allocatable :: filename_base  

  !> Output file format 
  character(len=esmf_maxstr),dimension(:),allocatable :: output_file 


  !> The first write task in a write group
  integer,dimension(:),allocatable     :: lead_wrttask 

  !> The last write task in a write group
  integer,dimension(:),allocatable     :: last_wrttask 

  !> Output grid type, e.g. "gaussian_grid"
  character(len=esmf_maxstr),dimension(:),allocatable :: output_grid 

  !> The i-dimension in the output grid
  integer,dimension(:),allocatable  :: imo 

  !> The j-dimension in the output grid
  integer,dimension(:),allocatable  :: jmo 

  !> Longitude of the center point in the output grid
  real,dimension(:),allocatable     :: cen_lon 

  !> Latitude of the center pointer in the output grid
  real,dimension(:),allocatable     :: cen_lat 

  !> Longitude of the first grid point in the output grid
  real,dimension(:),allocatable     :: lon1 

  !> Latitude of the first pointer in the output grid
  real,dimension(:),allocatable     :: lat1 

  !> Longitude of the last grid point in the output grid
  real,dimension(:),allocatable     :: lon2 

  !> Latitude of the last pointer in the output grid
  real,dimension(:),allocatable     :: lat2 

  !> Longitude increment
  real,dimension(:),allocatable     :: dlon 

  !> Latitude increment
  real,dimension(:),allocatable     :: dlat 

  !> The first latitude from the pole at which the secant cone cuts the sphere
  real,dimension(:),allocatable     :: stdlat1 

  !> The second latitude from the pole at which the secant cone cuts the sphere
  real,dimension(:),allocatable     :: stdlat2 

  !> x-direction grid length
  real,dimension(:),allocatable     :: dx 

  !> y-direction grid length
  real,dimension(:),allocatable     :: dy 

  !> Deflate level to use, 0 means no deflate.
  integer,dimension(:),allocatable  :: ideflate

  !> Number of significant digits for lossy compression.
  integer,dimension(:),allocatable  :: quantize_nsd

  !> Zstandard compression level, 0 means no zstandard compression.
  integer,dimension(:),allocatable  :: zstandard_level

  !> Quantize mode to use for lossy compression.
  character(len=esmf_maxstr),dimension(:),allocatable :: quantize_mode

  !> Chunk size in i dimension for 2D data.
  integer,dimension(:),allocatable  :: ichunk2d

  !> Chunk size in j dimension for 2D data.
  integer,dimension(:),allocatable  :: jchunk2d

  !> Chunk size in i dimension for 3D data.
  integer,dimension(:),allocatable  :: ichunk3d

  !> Chunk size in j dimension for 3D data.
  integer,dimension(:),allocatable  :: jchunk3d

  !> Chunk size in k dimension for 3D data.
  integer,dimension(:),allocatable  :: kchunk3d

end module module_fv3_io_def
