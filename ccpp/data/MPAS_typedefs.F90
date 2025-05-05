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

     ! Index for h2o mixing ratio.
     integer  :: index_qv

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
     ! Hydrostatic Pressure
     real(mpas_kind), pointer :: pmiddry(:,:)     ! Pressure at layer centers        (nlev,ncol)
     real(mpas_kind), pointer :: pintdry(:,:)     ! Pressure at layer interfaces     (nlev+1,ncol)

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
     integer :: nCells           ! Number of cells, including halo cells
     integer :: nEdges           ! Number of edges, including halo edges
     integer :: nVertices        ! Number of vertices, including halo vertices
     integer :: nVertLevels      ! Number of vertical layers
     !
     integer :: nCellsSolve      ! Number of cells, excluding halo cells
     integer :: nEdgesSolve      ! Number of edges, excluding halo edges
     integer :: nVerticesSolve   ! Number of vertices, excluding halo vertices

     ! MPAS vertical coordiante (invariant)
     real(mpas_kind), pointer :: zint(:,:)        ! Geometric height [m]  at layer interfaces (nlev+1,ncol)
     real(mpas_kind), pointer :: zz(:,:)          ! Vertical coordinate metric [1] at layer
                                                  ! midpoints (nlev,ncol)
     real(mpas_kind), pointer :: fzm(:)           ! Interp weight from k layer midpoint to k
                                                  ! layer interface [1] (nlev)
     real(mpas_kind), pointer :: fzp(:)           ! Interp weight from k-1 layer midpoint to k
                                                  ! layer interface [dimensionless] (nlev)

     ! Index for h2o mixing ratio.
     integer  :: index_qv

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
     real(mpas_kind), pointer :: pmiddry(:,:)     ! Dry hydrostatic pressure [Pa]
                                                  ! at layer midpoints               (nlev,ncol)
     real(mpas_kind), pointer :: pintdry(:,:)     ! Dry hydrostatic pressure [Pa]
                                                  ! at layer interfaces            (nlev+1,ncol)
     real(mpas_kind), pointer :: vorticity(:,:)   ! Relative vertical vorticity [s^-1]
                                                  !                                  (nlev,nvtx)
     real(mpas_kind), pointer :: divergence(:,:)  ! Horizontal velocity divergence [s^-1]
                                                  !                                  (nlev,ncol)
  end type MPAS_stateout_type

  public MPAS_init_type, MPAS_statein_type, MPAS_stateint_type, MPAS_stateout_type

end module MPAS_typedefs
