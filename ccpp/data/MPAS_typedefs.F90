! ###########################################################################################
!> \file MPAS_typedefs.F90
! ########################################################################################### 
module MPAS_typedefs
  use mpi_f08
  use machine, only: kind_phys, kind_dbl_prec, kind_sngl_prec
  use mpas_kind_types, only : mpas_kind => RKIND
  implicit none

!> \section arg_table_MPAS_typedefs
!! \htmlinclude MPAS_typedefs.html
!!

  ! #########################################################################################
  ! MPAS_control_type
  ! #########################################################################################
!! \section arg_table_MPAS_control_type
!! \htmlinclude MPAS_control_type.html
!!
  type MPAS_control_type

     ! Namelist filename
     character(len=64) :: fn_nml

     ! Full namelist for use with internal file reads
     character(len=:), pointer, dimension(:) :: input_nml_file => null()

     ! MPI Bookkeeping
     integer        :: me             !< current MPI-rank
     integer        :: master         !< master MPI-rank
     type(MPI_Comm) :: mpi_comm       !< forecast tasks mpi communicator

     ! ESMF 
     integer          :: fcst_ntasks  !< total number of forecast tasks

     ! Log file identifier
     integer          :: nlunit       !< fortran unit number for file opens
     integer          :: logunit      !< fortran unit number for writing logfile

     ! UFS date(s) for model time.
     integer          :: bdat(8)      !< model begin date in GFS format   (same as idat)
     integer          :: cdat(8)      !< model current date in GFS format (same as jdat)

     ! Spatial/Temporal parameters for physics/dynamics coupling.
     real(kind_phys)  :: dt_dycore    !< dynamics time step in seconds
     real(kind_phys)  :: dt_phys      !< physics  time step in seconds
     integer          :: nblks        !< Number of data (physics) blocks.
     integer, pointer :: blksz(:)     !< Block size for  data blocking (default blksz(1)=[nCells])
     integer          :: levs         !< number of vertical levels

     ! Tracer information
     integer                    :: nConstituents   !< Number of constituents (tracers).
     integer                    :: nwat            !< number of hydrometeors in dcyore (including water vapor)
     character(len=32), pointer :: tracer_names(:) !< tracers names to dereference tracer id
     integer,           pointer :: tracer_types(:) !< tracers types: 0=generic, 1=chem,prog, 2=chem,diag

     ! #####################################################################################
     ! These are fields are not needed for MPAS, but instead are needed to call
     ! control_initialize (i.e. GFS physics initialization) from  within MPAS_INIT().
     !
     ! Need to partition GFS_CONTROL_TYPE, in GFS_TYPEDEFS.F90, into FV3 and MPAS pieces.
     !
     ! GFS_CONTROL_TYPE => GFS_CONTROL_TYPE       - Dycore agnostic physics configuration
     !                 \                           (Most stuff in gfs_control_type will remain in here)
     !                  +  FV3_CONTROL_TYPE       - FV3 dycore specific physics configuration
     !                                             (all this stuff below and any derived fields
     !                                              in control_initialize shall reside here)
     !                  +  MPAS_CONTROL_TYPE      - MPAS dycore specific physics configuration
     !                                             (TBD)
     ! ##################################################################################### 
     !real(kind=kind_phys), pointer :: ak(:)       !< from surface (k=1) to TOA (k=levs)
     !real(kind=kind_phys), pointer :: bk(:)       !< from surface (k=1) to TOA (k=levs)
     !integer :: isc                               !< starting i-index for this MPI-domain
     !integer :: jsc                               !< starting j-index for this MPI-domain
     !integer :: nx                                !< number of points in i-dir for this MPI rank
     !integer :: ny                                !< number of points in j-dir for this MPI rank
     !integer :: cnx                               !< number of points in i-dir for this cubed-sphere face
     !                                             !< equal to gnx for lat-lon grids
     !integer :: cny                               !< number of points in j-dir for this cubed-sphere face
     !                                             !< equal to gny for lat-lon grids
     !integer :: gnx                               !< number of global points in x-dir (i) along the equator
     !integer :: gny                               !< number of global points in y-dir (j) along any meridian
     integer :: iau_offset                        !< iau running window length
     !integer :: tile_num                          !< tile number for this MPI rank
     logical :: restart                           !< flag whether this is a coldstart (.false.) or a warmstart/restart (.true.)
     !logical :: hydrostatic                       !< flag whether this is a hydrostatic or non-hydrostatic run
  end type MPAS_control_type

  ! #########################################################################################
  ! MPAS_statein_type
  !  Prognostic state variables INTO dycore.
  ! #########################################################################################
!! \section arg_table_MPAS_statein_type
!! \htmlinclude MPAS_statein_type.html
!!
  type MPAS_statein_type
     ! Dimensions
     integer :: nCells                            ! Number of cells, including halo cells
     integer :: nEdges                            ! Number of edges, including halo edges
     integer :: nVertices                         ! Number of vertices, including halo vertices
     integer :: nVertLevels                       ! Number of vertical layers
     !
     integer :: nCellsSolve                       ! Number of cells, excluding halo cells
     integer :: nEdgesSolve                       ! Number of edges, excluding halo edges
     integer :: nVerticesSolve                    ! Number of vertices, excluding halo vertices

     ! MPAS vertical coordiante (invariant)
     real(mpas_kind), pointer :: zint(:,:)        ! Geometric height [m]  at layer interfaces (nlev+1,ncol)
     real(mpas_kind), pointer :: zz(:,:)          ! Vertical coordinate metric [1] at layer
                                                  ! midpoints (nlev,ncol)
     real(mpas_kind), pointer :: fzm(:)           ! Interp weight from k layer midpoint to k
                                                  ! layer interface [1] (nlev)
     real(mpas_kind), pointer :: fzp(:)           ! Interp weight from k-1 layer midpoint to k
                                                  ! layer interface [dimensionless] (nlev)
     ! Cell area (invariant)
     real(mpas_kind), pointer :: areaCell(:)      ! cell area [m^2]

     ! For edge-normal velocity calculations (invariant)
     real(mpas_kind), pointer :: east(:,:)        ! Cartesian components of unit east vector
                                                  ! at cell centers [dimensionless]       (3,ncol)
     real(mpas_kind), pointer :: north(:,:)       ! Cartesian components of unit north vector
                                                  ! at cell centers [dimensionless]       (3,ncol)
     real(mpas_kind), pointer :: normal(:,:)      ! Cartesian components of the vector normal
                                                  ! to an edge and tangential to the surface
                                                  ! of the sphere [dimensionless]         (3,ncol)
     integer, pointer :: cellsOnEdge(:,:)         ! Indices of cells separated by an edge (2,nedge)

     ! Indices for tracer (scalar) indices
     integer  :: index_qv                         ! Tracer index for water-vapor mixing-ratio
     integer  :: index_qc                         ! Tracer index for cloud-water mixing-ratio
     integer  :: index_qr                         ! Tracer index for rain-water mixing-ratio
     integer  :: index_qs                         ! Tracer index for snow mixing-ratio
     integer  :: index_qi                         ! Tracer index for ice mixing ratio
     integer  :: index_qh                         ! Tracer index for hail mixing ratio

     ! Base state variables
     real(mpas_kind), pointer :: rho_base(:,:)    ! Base-state dry air density [kg/m^3]  (nlev,ncol)
     real(mpas_kind), pointer :: theta_base(:,:)  ! Base-state potential temperature [K] (nlev,ncol)

     ! State that is directly prognosed by the dycore
     real(mpas_kind), pointer :: uperp(:,:)       ! Normal velocity at edges [m/s]  (nlev  ,nedge)
     real(mpas_kind), pointer :: w(:,:)           ! Vertical velocity [m/s]         (nlev+1,ncol)
     real(mpas_kind), pointer :: theta_m(:,:)     ! Moist potential temperature [K] (nlev  ,ncol)
     real(mpas_kind), pointer :: rho_zz(:,:)      ! Dry density [kg/m^3]
                                                  ! divided by d(zeta)/dz            (nlev ,ncol)
     real(mpas_kind), pointer :: tracers(:,:,:)   ! Tracers [kg/kg dry air]       (nq,nlev ,ncol)


     ! State that may be directly derived from dycore prognostic state
     real(mpas_kind), pointer :: theta(:,:)       ! Potential temperature [K]        (nlev,ncol)
     real(mpas_kind), pointer :: exner(:,:)       ! Exner function [-]               (nlev,ncol)
     real(mpas_kind), pointer :: rho(:,:)         ! Dry density [kg/m^3]             (nlev,ncol)
     real(mpas_kind), pointer :: ux(:,:)          ! Zonal veloc at center [m/s]      (nlev,ncol)
     real(mpas_kind), pointer :: uy(:,:)          ! Meridional veloc at center [m/s] (nlev,ncol)

     ! Tendencies from physics
     real(mpas_kind), pointer :: ru_tend(:,:)     ! Normal horizontal momentum tendency
                                                  ! from physics [kg/m^2/s]          (nlev,nedge)
     real(mpas_kind), pointer :: rtheta_tend(:,:) ! Tendency of rho*theta/zz
                                                  ! from physics [kg K/m^3/s]        (nlev,ncol)
     real(mpas_kind), pointer :: rho_tend(:,:)    ! Dry air density tendency
                                                  ! from physics [kg/m^3/s]          (nlev,ncol)

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
  end type MPAS_stateint_type
  
  ! #########################################################################################
  ! MPAS_stateout_type
  !  Prognostic state or tendencies FROM MPAS dycore.
  ! #########################################################################################
!! \section arg_table_MPAS_stateout_type
!! \htmlinclude MPAS_stateout_type.html
!!
  type MPAS_stateout_type
     ! Dimensions
     integer :: nCells                            ! Number of cells, including halo cells
     integer :: nEdges                            ! Number of edges, including halo edges
     integer :: nVertices                         ! Number of vertices, including halo vertices
     integer :: nVertLevels                       ! Number of vertical layers
     !
     integer :: nCellsSolve                       ! Number of cells, excluding halo cells
     integer :: nEdgesSolve                       ! Number of edges, excluding halo edges
     integer :: nVerticesSolve                    ! Number of vertices, excluding halo vertices
     
     ! MPAS vertical coordiante (invariant)
     real(mpas_kind), pointer :: zint(:,:)        ! Geometric height [m]  at layer interfaces (nlev+1,ncol)
     real(mpas_kind), pointer :: zz(:,:)          ! Vertical coordinate metric [1] at layer
                                                  ! midpoints (nlev,ncol)
     real(mpas_kind), pointer :: fzm(:)           ! Interp weight from k layer midpoint to k
                                                  ! layer interface [1] (nlev)
     real(mpas_kind), pointer :: fzp(:)           ! Interp weight from k-1 layer midpoint to k
                                                  ! layer interface [dimensionless] (nlev)

     ! Indices for tracer (scalar) indices
     integer  :: index_qv                         ! Tracer index for water-vapor mixing-ratio
     integer  :: index_qc                         ! Tracer index for cloud-water mixing-ratio
     integer  :: index_qr                         ! Tracer index for rain-water mixing-ratio
     integer  :: index_qs                         ! Tracer index for snow mixing-ratio
     integer  :: index_qi                         ! Tracer index for ice mixing ratio
     integer  :: index_qh                         ! Tracer index for hail mixing ratio

     ! State that is directly prognosed by the dycore
     real(mpas_kind), pointer :: uperp(:,:)       ! Normal velocity at edges [m/s]  (nlev  ,nedge)
     real(mpas_kind), pointer :: w(:,:)           ! Vertical velocity [m/s]         (nlev+1,ncol)
     real(mpas_kind), pointer :: theta_m(:,:)     ! Moist potential temperature [K] (nlev  ,ncol)
     real(mpas_kind), pointer :: rho_zz(:,:)      ! Dry density [kg/m^3]
                                                  ! divided by d(zeta)/dz            (nlev ,ncol)
     real(mpas_kind), pointer :: tracers(:,:,:)   ! Tracers [kg/kg dry air]       (nq,nlev ,ncol)

     ! State that may be directly derived from dycore prognostic state (ToDo)
     real(mpas_kind), pointer :: theta(:,:)       ! Potential temperature [K]        (nlev,ncol)
     real(mpas_kind), pointer :: exner(:,:)       ! Exner function [-]               (nlev,ncol)
     real(mpas_kind), pointer :: rho(:,:)         ! Dry density [kg/m^3]             (nlev,ncol)
     real(mpas_kind), pointer :: ux(:,:)          ! Zonal veloc at center [m/s]      (nlev,ncol)
     real(mpas_kind), pointer :: uy(:,:)          ! Meridional veloc at center [m/s] (nlev,ncol)
     real(mpas_kind), pointer :: pmiddry(:,:)     ! Dry hydrostatic pressure [Pa]
                                                  ! at layer midpoints               (nlev,ncol)
     real(mpas_kind), pointer :: pintdry(:,:)     ! Dry hydrostatic pressure [Pa]
                                                  ! at layer interfaces            (nlev+1,ncol)
     real(mpas_kind), pointer :: vorticity(:,:)   ! Relative vertical vorticity [s^-1]
                                                  !                                  (nlev,nvtx)
     real(mpas_kind), pointer :: divergence(:,:)  ! Horizontal velocity divergence [s^-1]
                                                  !                                  (nlev,ncol)
  end type MPAS_stateout_type

  public MPAS_control_type, MPAS_statein_type, MPAS_stateint_type, MPAS_stateout_type

end module MPAS_typedefs
