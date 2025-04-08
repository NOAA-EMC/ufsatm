! ###########################################################################################
!> \file MPAS_typedefs.F90
! ########################################################################################### 
module MPAS_typedefs

  use mpi_f08
  use machine, only: kind_phys, kind_dbl_prec, kind_sngl_prec
  implicit none

!> \section arg_table_MPAS_typedefs
!! \htmlinclude MPAS_typedefs.html
!!

  ! MPAS_init_type           !<
  ! MPAS_statein_type        !< prognostic state data in from dycore
  ! MPAS_stateint_type       !< prognostic state or tendencies before call MP.
  ! MPAS_stateout_type       !< prognostic state or tendencies return to dycore

  ! #########################################################################################
  ! MPAS_init_type
  ! #########################################################################################
!! \section arg_table_MPAS_init_type
!! \htmlinclude MPAS_init_type.html
!!
  type MPAS_init_type
     integer :: me                    !< my MPI-rank
     integer :: master                !< master MPI-rank
     type(MPI_Comm) :: mpi_comm       !< forecast tasks mpi communicator
     integer :: fcst_ntasks           !< total number of forecast tasks
     integer :: nlunit                !< fortran unit number for file opens
     integer :: logunit               !< fortran unit number for writing logfile
     integer :: bdat(8)               !< model begin date in GFS format   (same as idat)
     integer :: cdat(8)               !< model current date in GFS format (same as jdat)

     real(kind_phys) :: dt_dycore     !< dynamics time step in seconds
     real(kind_phys) :: dt_phys       !< physics  time step in seconds
     integer, pointer :: blksz(:)     !< for explicit data blocking
                                      !< default blksz(1)=[nx*ny]
     integer :: levs                  !< number of vertical levels

     integer                    :: nConstituents   !< Number of constituents (tracers).
     integer                    :: nwat            !< number of hydrometeors in dcyore (including water vapor)
     character(len=32), pointer :: tracer_names(:) !< tracers names to dereference tracer id
     integer,           pointer :: tracer_types(:) !< tracers types: 0=generic, 1=chem,prog, 2=chem,diag
     character(len=64) :: fn_nml                   !< namelist filename
     character(len=:), pointer, dimension(:) :: input_nml_file => null() !< character string containing full namelist
                                                                         !< for use with internal file reads

     ! NOT NEEDED FOR MPAS, BUT NEEDED FOR CONTROL_INITIALZE. NEED TO PARTITION
     real(kind=kind_phys), pointer :: ak(:)       !< from surface (k=1) to TOA (k=levs)
     real(kind=kind_phys), pointer :: bk(:)       !< from surface (k=1) to TOA (k=levs)
     integer :: isc                               !< starting i-index for this MPI-domain
     integer :: jsc                               !< starting j-index for this MPI-domain
     integer :: nx                                !< number of points in i-dir for this MPI rank
     integer :: ny                                !< number of points in j-dir for this MPI rank
     integer :: cnx                               !< number of points in i-dir for this cubed-sphere face
                                                  !< equal to gnx for lat-lon grids
     integer :: cny                               !< number of points in j-dir for this cubed-sphere face
                                                  !< equal to gny for lat-lon grids
     integer :: gnx                               !< number of global points in x-dir (i) along the equator
     integer :: gny                               !< number of global points in y-dir (j) along any meridian
     integer :: iau_offset                        !< iau running window length
     integer :: tile_num                          !< tile number for this MPI rank
     logical :: restart                           !< flag whether this is a coldstart (.false.) or a warmstart/restart (.true.)
     logical :: hydrostatic                       !< flag whether this is a hydrostatic or non-hydrostatic run
  end type MPAS_init_type

  ! #########################################################################################
  ! MPAS_statein_type
  !  Prognostic state variables with layer and level specific data from dycore.
  ! #########################################################################################
!! \section arg_table_MPAS_statein_type
!! \htmlinclude MPAS_statein_type.html
!!
  type MPAS_statein_type
     real (kind_phys), pointer :: theta (:,:) => null() !< Potential temperature
     real (kind_phys), pointer :: u (:,:)     => null() !< Zonal wind
     real (kind_phys), pointer :: v (:,:)     => null() !< Meridional wind
     real (kind_phys), pointer :: rho (:,:)   => null() !< Dry air density
     real (kind_phys), pointer :: q (:,:,:)   => null() !< Tracers
   contains
     procedure :: create => statein_create
  end type MPAS_statein_type
  
  ! #########################################################################################
  ! MPAS_stateint_type
  !  Prognostic state or tendencies BEFORE calling microphysics.
  ! #########################################################################################
!! \section arg_table_MPAS_stateint_type
!! \htmlinclude MPAS_stateint_type.html
!!
  type MPAS_stateint_type
     real (kind_phys), pointer :: u (:,:)    => null()  !< Zonal wind
     real (kind_phys), pointer :: v (:,:)    => null()  !< Meridional wind
     real (kind_phys), pointer :: temp (:,:) => null()  !< Temperature
     real (kind_phys), pointer :: q (:,:,:)  => null()  !< Tracers
   contains
     procedure :: create  => stateint_create
  end type MPAS_stateint_type
  
  ! #########################################################################################
  ! MPAS_stateout_type
  !  Prognostic state or tendencies after ALL physical parameterizations.
  ! #########################################################################################
!! \section arg_table_MPAS_stateout_type
!! \htmlinclude MPAS_stateout_type.html
!!
  type MPAS_stateout_type
     real (kind_phys), pointer :: u (:,:)    => null()  !< Zonal wind
     real (kind_phys), pointer :: v (:,:)    => null()  !< Meridional wind
     real (kind_phys), pointer :: temp (:,:) => null()  !< Temperature
     real (kind_phys), pointer :: q (:,:,:)  => null()  !< Tracers
   contains
     procedure :: create  => stateout_create
  end type MPAS_stateout_type

  public MPAS_init_type, MPAS_statein_type, MPAS_stateint_type, MPAS_stateout_type

contains
  ! #########################################################################################
  ! Allocation pricedures (type-bound)
  ! #########################################################################################
  subroutine statein_create(Statein, Model)
    use GFS_typedefs, only: GFS_control_type
    implicit none
    class(MPAS_statein_type)           :: Statein
    type(GFS_control_type), intent(in) :: Model
    allocate (Statein%theta (Model%ncols,Model%levs))
    allocate (Statein%u     (Model%ncols,Model%levs))
    allocate (Statein%v     (Model%ncols,Model%levs))
    allocate (Statein%rho   (Model%ncols,Model%levs))
    allocate (Statein%q     (Model%ncols,Model%levs,Model%ntrac))
    !
    Statein%theta = 0.0_kind_phys
    Statein%u     = 0.0_kind_phys
    Statein%v     = 0.0_kind_phys
    Statein%rho   = 0.0_kind_phys
    Statein%q     = 0.0_kind_phys

  end subroutine statein_create

  ! #########################################################################################
  !
  ! #########################################################################################
  subroutine stateout_create (Stateout, Model)
    use GFS_typedefs, only: GFS_control_type
    implicit none
    class(MPAS_stateout_type)          :: Stateout
    type(GFS_control_type), intent(in) :: Model

    allocate (Stateout%u    (Model%ncols,Model%levs))
    allocate (Stateout%v    (Model%ncols,Model%levs))
    allocate (Stateout%temp (Model%ncols,Model%levs))
    allocate (Stateout%q    (Model%ncols,Model%levs,Model%ntrac))
    !
    Stateout%u    = 0.0_kind_phys
    Stateout%v    = 0.0_kind_phys
    Stateout%temp = 0.0_kind_phys
    Stateout%q    = 0.0_kind_phys

  end subroutine stateout_create

  ! #########################################################################################
  !
  ! #########################################################################################
  subroutine stateint_create (Stateint, Model)
    use GFS_typedefs, only: GFS_control_type
    implicit none
    class(MPAS_stateint_type)          :: Stateint
    type(GFS_control_type), intent(in) :: Model

    allocate (Stateint%u    (Model%ncols,Model%levs))
    allocate (Stateint%v    (Model%ncols,Model%levs))
    allocate (Stateint%temp (Model%ncols,Model%levs))
    allocate (Stateint%q    (Model%ncols,Model%levs,Model%ntrac))
    !
    Stateint%u    = 0.0_kind_phys
    Stateint%v    = 0.0_kind_phys
    Stateint%temp = 0.0_kind_phys
    Stateint%q    = 0.0_kind_phys

  end subroutine stateint_create

end module MPAS_typedefs
