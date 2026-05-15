! ###########################################################################################
!> \file atmos_coupling.F90
!> Procedures for coupling the MPAS dynamical core to the CCPP Physics.
!>
! ###########################################################################################
module atmos_coupling_mod
  use mpas_kind_types, only : mpas_kind => RKIND
  use ufs_mpas_io,     only : domain_ptr
  
  implicit none
  public :: MPAS_statein_type
  public :: MPAS_stateout_type
  public :: ufs_mpas_to_physics
  public :: ufs_microphysics_to_mpas
  public :: ufs_mpas_to_microphysics
  public :: ufs_mpas_grid_to_physics

  !> #######################################################################################
  !> MPAS_statein_type
  !> Fields needed by the MPAS dynamical core for forward integration.
  !>
  !> #######################################################################################
  type MPAS_statein_type
     ! Dimensions
     integer, pointer :: nCells                   ! Number of cells, including halo cells
     integer, pointer :: nEdges                   ! Number of edges, including halo edges
     integer, pointer :: nVertices                ! Number of vertices, including halo vertices
     integer, pointer :: nVertLevels              ! Number of vertical layers
     !
     integer, pointer :: nCellsSolve              ! Number of cells, excluding halo cells
     integer, pointer :: nEdgesSolve              ! Number of edges, excluding halo edges
     integer, pointer :: nVerticesSolve           ! Number of vertices, excluding halo vertices

     ! MPAS vertical coordiante (invariant)
     real(mpas_kind), pointer :: zgrid(:,:)       ! Geometric height [m]  at layer interfaces (nlev+1,ncol)
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
     integer, pointer  :: index_qv                ! Tracer index for water-vapor mixing-ratio

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

     ! Diagnostics
     real(mpas_kind), pointer :: pressure_b(:,:)
     real(mpas_kind), pointer :: pressure_p(:,:)
     real(mpas_kind), pointer :: surface_pressure(:)

  end type MPAS_statein_type

  !> #######################################################################################
  !> MPAS_stateout_type
  !> Fields prognosed (or diagnosed) by the MPAS dynamical core.
  !>
  !> #######################################################################################
    type MPAS_stateout_type
     ! Dimensions
     integer, pointer :: nCells                   ! Number of cells, including halo cells
     integer, pointer :: nEdges                   ! Number of edges, including halo edges
     integer, pointer :: nVertices                ! Number of vertices, including halo vertices
     integer, pointer :: nVertLevels              ! Number of vertical layers
     !
     integer, pointer :: nCellsSolve              ! Number of cells, excluding halo cells
     integer, pointer :: nEdgesSolve              ! Number of edges, excluding halo edges
     integer, pointer :: nVerticesSolve           ! Number of vertices, excluding halo vertices
     
     ! MPAS vertical coordiante (invariant)
     real(mpas_kind), pointer :: zgrid(:,:)       ! Geometric height [m]  at layer interfaces (nlev+1,ncol)
     real(mpas_kind), pointer :: zz(:,:)          ! Vertical coordinate metric [1] at layer
                                                  ! midpoints (nlev,ncol)
     real(mpas_kind), pointer :: fzm(:)           ! Interp weight from k layer midpoint to k
                                                  ! layer interface [1] (nlev)
     real(mpas_kind), pointer :: fzp(:)           ! Interp weight from k-1 layer midpoint to k
                                                  ! layer interface [dimensionless] (nlev)

     ! Indices for tracer (scalar) indices
     integer, pointer  :: index_qv                ! Tracer index for water-vapor mixing-ratio
     
     ! State that is directly prognosed by the dycore
     real(mpas_kind), pointer :: uperp(:,:)       ! Normal velocity at edges [m/s]  (nlev  ,nedge)
     real(mpas_kind), pointer :: w(:,:)           ! Vertical velocity [m/s]         (nlev+1,ncol)
     real(mpas_kind), pointer :: theta_m(:,:)     ! Moist potential temperature [K] (nlev  ,ncol)
     real(mpas_kind), pointer :: rho_zz(:,:)      ! Dry density [kg/m^3]
                                                  ! divided by d(zeta)/dz            (nlev ,ncol)
     real(mpas_kind), pointer :: tracers(:,:,:)   ! Tracers [kg/kg dry air]       (nq,nlev ,ncol)

     ! State that may be directly derived from dycore prognostic state.
     real(mpas_kind), pointer :: theta(:,:)       ! Potential temperature [K]        (nlev,ncol)
     real(mpas_kind), pointer :: exner(:,:)       ! Exner function [-]               (nlev,ncol)
     real(mpas_kind), pointer :: rho(:,:)         ! Dry density [kg/m^3]             (nlev,ncol)
     real(mpas_kind), pointer :: ux(:,:)          ! Zonal veloc at center [m/s]      (nlev,ncol)
     real(mpas_kind), pointer :: uy(:,:)          ! Meridional veloc at center [m/s] (nlev,ncol)
     real(mpas_kind), pointer :: pmiddry(:,:)     ! Dry hydrostatic pressure [Pa]
                                                  ! at layer midpoints               (nlev,ncol)
     real(mpas_kind), pointer :: pintdry(:,:)     ! Dry hydrostatic pressure [Pa]
                                                  ! at layer interfaces            (nlev+1,ncol)
     real(mpas_kind), pointer :: pmid(:,:)        ! Pressure at layer midpoints      (nlev,ncol)
     real(mpas_kind), pointer :: vorticity(:,:)   ! Relative vertical vorticity [s^-1]
                                                  !                                  (nlev,nvtx)
     real(mpas_kind), pointer :: divergence(:,:)  ! Horizontal velocity divergence [s^-1]
                                                  !                                  (nlev,ncol)
     ! Diagnostics
     real(mpas_kind), pointer :: pressure_b(:,:)
     real(mpas_kind), pointer :: pressure_p(:,:)
     real(mpas_kind), pointer :: surface_pressure(:)

  end type MPAS_stateout_type
  
contains
  !> #########################################################################################
  !> Procedure to convert input "MPAS" variables to "CCPP" variables.
  !> Called prior to MPAS dynamical core (initial-step only).
  !>
  !> Analogous to MPAS_to_physics in src/core_atmosphere/physics/mpas_atmphys_interface.F
  !>
  !> This procedure accesses MPAS data using MPAS native procedures and stores the data
  !> locally in the data-containers defined above. The MPAS "state" is then translated to the
  !> CCPP "state" needed by the physics.
  !>
  !> #########################################################################################
  subroutine ufs_mpas_to_physics(physics_state, surface_state)
    use GFS_typedefs,         only : GFS_statein_type, GFS_sfcprop_type
    use mpas_derived_types,   only : mpas_pool_type
    use mpas_pool_routines,   only : mpas_pool_get_subpool, mpas_pool_get_array, mpas_pool_get_dimension
    use atm_core,             only : atm_compute_output_diagnostics
    use mpas_kind_types,      only : RKIND
    use mpas_constants,       only : gravity

    ! Arguments
    type(GFS_statein_type),   intent(inout) :: physics_state
    type(GFS_sfcprop_type),   intent(inout) :: surface_state
    ! Locals
    type(mpas_stateout_type) :: mpas_state
    type(mpas_pool_type), pointer :: state_pool
    type(mpas_pool_type), pointer :: diag_pool
    type(mpas_pool_type), pointer :: mesh_pool
    type(mpas_pool_type), pointer :: sfc_pool
    integer :: iCol, iLay, iTracer
    integer, pointer :: nCellsSolve, num_scalars, nwat, nVertLevels, index_qv
    integer, dimension(:), pointer :: isltyp
    real(kind=RKIND) :: rho1, rho2, tem1, tem2
    real(kind=RKIND),dimension(:,:),pointer  :: qv, qc, qr, qi, qs, qg

    ! Access MPAS data pools.
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'state',     state_pool)
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'diag',      diag_pool)
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh',      mesh_pool)

    ! DJS to GFS: Sanity check to ensure data is in "sfc_pool" to pass to physics types.
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'sfc_input', sfc_pool)
    call mpas_pool_get_array(sfc_pool, 'isltyp', isltyp, 1)

    ! Get MPAS dimensions
    call mpas_pool_get_dimension(mesh_pool,  'nCellsSolve', nCellsSolve)
    call mpas_pool_get_dimension(state_pool, 'num_scalars', num_scalars)
    call mpas_pool_get_dimension(state_pool, 'index_qv',    index_qv)
    call mpas_pool_get_dimension(state_pool, 'moist_end',   nwat)
    call mpas_pool_get_dimension(mesh_pool,  'nVertLevels', nVertLevels)

    ! Grab fields from MPAS pools
    call mpas_pool_get_array(diag_pool,  'theta',                  MPAS_state % theta)
    call mpas_pool_get_array(diag_pool,  'uReconstructZonal',      MPAS_state % ux)
    call mpas_pool_get_array(diag_pool,  'uReconstructMeridional', MPAS_state % uy)
    call mpas_pool_get_array(state_pool, 'scalars',                MPAS_state % tracers, timeLevel=1)
    call mpas_pool_get_array(state_pool, 'w',                      MPAS_state % w, timeLevel=1)
    call mpas_pool_get_array(diag_pool,  'exner',                  MPAS_state % exner)
    call mpas_pool_get_array(mesh_pool,  'zgrid',                  MPAS_state % zgrid)
    call mpas_pool_get_array(mesh_pool,  'zz',                     MPAS_state % zz)
    call mpas_pool_get_array(state_pool, 'theta_m',                MPAS_state % theta_m, timeLevel=1)
    call mpas_pool_get_array(state_pool, 'rho_zz',                 MPAS_state % rho_zz,  timeLevel=1)
    call mpas_pool_get_array(diag_pool,  'pressure_base',          MPAS_state % pressure_b)
    call mpas_pool_get_array(diag_pool,  'pressure_p',             MPAS_state % pressure_p)

    ! Copy fields from MPAS data containers to physics data containers.
    ! [k, i] -> [i, k]
    ! Retain bottom-up convention
    do iCol = 1, nCellsSolve
       physics_state % tgrs(iCol,:)   = MPAS_state % theta(:,iCol)*MPAS_state % exner(:,iCol)
       physics_state % ugrs(iCol,:)   = MPAS_state % ux(:,iCol)
       physics_state % vgrs(iCol,:)   = MPAS_state % uy(:,iCol)
       physics_state % phil(iCol,:)   = MPAS_state % zz(:,iCol)
       physics_state % phii(iCol,:)   = MPAS_state % zgrid(:,iCol)
       physics_state % prslk(iCol,:)  = MPAS_state % exner(:,iCol)
       ! MPAS provides vertical velocity at interfaces, compute layer mean.
       do iLay=1,nVertLevels
          physics_state % vvl(iCol,iLay) = 0.5*(MPAS_state % w(iLay,iCol) + MPAS_state % w(iLay+1,iCol))
       enddo
       do iTracer = 1,num_scalars
          physics_state % qgrs(iCol,:,iTracer) = MPAS_state % tracers(iTracer,:,iCol)
       enddo
    enddo    

    ! Set surface temperature to lowest level temperature (revisit for coupling)
    do iCol = 1, nCellsSolve
       surface_state % tsfc(iCol) = MPAS_state % theta(1,iCol)*MPAS_state % exner(1,iCol)
    enddo

    ! Calculation of the surface pressure using hydrostatic assumption down to the surface.
    ! (from mpas_atmphys_interface.F:MPAS_to_physics())
    call mpas_pool_get_array(diag_pool,   'surface_pressure'      ,MPAS_state % surface_pressure)
    do iCol = 1, nCellsSolve
       tem1 = MPAS_state % zgrid(2,iCol) - MPAS_state % zgrid(1,iCol)
       tem2 = MPAS_state % zgrid(3,iCol) - MPAS_state % zgrid(2,iCol)
       rho1 = MPAS_state % rho_zz(1,iCol) * MPAS_state % zz(1,iCol) * (1. + MPAS_state % tracers(index_qv,1,iCol))
       rho2 = MPAS_state % rho_zz(2,iCol) * MPAS_state % zz(2,iCol) * (1. + MPAS_state % tracers(index_qv,2,iCol))
       MPAS_state % surface_pressure(iCol) = 0.5*gravity*(MPAS_state % zgrid(2,iCol) - MPAS_state % zgrid(1,iCol)) &
            * (rho1 - 0.5*(rho2-rho1)*tem1/(tem1+tem2))
       MPAS_state % surface_pressure(iCol) = MPAS_state % surface_pressure(iCol) + &
                                             MPAS_state % pressure_p(1,iCol) + &
                                             MPAS_state % pressure_b(1,iCol)
    enddo


    ! Compute hydrostatic pressures
    allocate(MPAS_state % pmid(   nVertLevels,   nCellsSolve))
    allocate(MPAS_state % pmiddry(nVertLevels,   nCellsSolve))
    allocate(MPAS_state % pintdry(nVertLevels+1, nCellsSolve))
    call hydrostatic_pressure(nCellsSolve, nVertLevels, nwat, index_qv, MPAS_state % zz,    &
         MPAS_state % zgrid, MPAS_state % rho_zz, MPAS_state % theta_m, MPAS_state % exner,  &
         MPAS_state % tracers, MPAS_state % pmiddry, MPAS_state % pintdry, MPAS_state % pmid)

    ! Copy MPAS pressures into physics data containers.
    ! [k, i] -> [i, k]
    ! Retain bottom-up convention
    do iCol = 1, nCellsSolve
       physics_state % pgr(iCol)    = MPAS_state % pintdry(nVertLevels+1,iCol)
       physics_state % prsl(iCol,:) = MPAS_state % pmiddry(:,iCol)
       physics_state % prsi(iCol,:) = MPAS_state % pintdry(:,iCol)
    enddo
    ! Housekeeping
    nullify (mesh_pool)
    nullify (state_pool)
    nullify (diag_pool)

  end subroutine ufs_mpas_to_physics

  !> #########################################################################################
  !> Procedure to update state with physics tendencies prior to calling MPAS dynamical core.
  !>
  !> Analogous to phys_get_tend in physics/mpas_atmphys_todynamics.F
  !> Instead of updating the state with physics tendencies from the MPAS "tend_pool", we
  !> will use tendencies from the CCPP Physics.
  !>
  !> #########################################################################################
  subroutine ufs_physics_to_mpas(control, diag)
    use GFS_typedefs,            only : GFS_diag_type, GFS_control_type
    use mpas_derived_types,      only : mpas_pool_type
    use mpas_pool_routines,      only : mpas_pool_get_subpool, mpas_pool_get_array, mpas_pool_get_config
    use mpas_stochastic_physics, only : stochastic_physics_pattern_apply
    use mpas_kind_types,         only : RKIND

    ! Arguments
    type(GFS_diag_type),    intent(in) :: diag
    type(GFS_control_type), intent(in) :: control
    ! Locals
    character(len=32)  :: tend_names(4)
    type(mpas_pool_type), pointer :: tend_physics_pool
    real(kind=RKIND), pointer     :: tend_array(:,:)
    logical, pointer :: do_sppt
    integer :: ierr
    
    call mpas_pool_get_config(domain_ptr % blocklist % configs, 'do_sppt', do_sppt)
    
    ! Call stochastic physics (SPPT) for (hydrostatic) physics tendencies.
    if (do_sppt) then
       ! Grab tend_physics pool from MPAS, populate with CCPP Physics data.
       call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'tend_physics', tend_physics_pool)
       !
       call mpas_pool_get_array(tend_physics_pool, 'rucuten', tend_array)
       tend_array(:,:) = &
            diag%dtend(:,:,control%dtidx(control%index_of_x_wind, control%index_of_process_dcnv)) + &
            diag%dtend(:,:,control%dtidx(control%index_of_x_wind, control%index_of_process_scnv))
       call mpas_pool_get_array(tend_physics_pool, 'rvcuten', tend_array)
       tend_array(:,:) = &
            diag%dtend(:,:,control%dtidx(control%index_of_y_wind, control%index_of_process_dcnv)) + &
            diag%dtend(:,:,control%dtidx(control%index_of_y_wind, control%index_of_process_scnv))
       call mpas_pool_get_array(tend_physics_pool, 'rublten', tend_array)
       tend_array(:,:) = &
            diag%dtend(:,:, control%dtidx(control%index_of_x_wind, control%index_of_process_pbl))
       call mpas_pool_get_array(tend_physics_pool, 'rvblten', tend_array)
       tend_array(:,:) = &
            diag%dtend(:,:, control%dtidx(control%index_of_y_wind, control%index_of_process_pbl))

       ! Apply pattern to tendencies
       !call stochastic_physics_pattern_apply(domain_ptr, 'phys', ierr)
       tend_names(1) = "rucuten"
       tend_names(2) = "rvcuten"
       tend_names(3) = "rublten"
       tend_names(4) = "rvblten"
       call stochastic_physics_pattern_apply(domain_ptr, 4, tend_names, ierr)
    end if




    !
    ! ACCUMMULATE PHYSICS TENDENCIES, CONVERT TO (NON-HYDROSTATIC) DYCORE PROGNOSIC STATE VARIABLES
    !









    
    ! Call stochastic physics (SPPT) for (non-hydrostatic) physics tendencies. 
    if (do_sppt) then

       call mpas_pool_get_array(tend_physics_pool, 'tend_rtheta_physics', tend_array)
       !call mpas_pool_get_array(tend_physics_pool, 'tend_rho_physics', tend_array)
       tend_names(1) = "tend_rtheta_physics"
       !tend_names(2) = "tend_rho_physics"

       ! Apply pattern to tendencies
       !call stochastic_physics_pattern_apply(domain_ptr, 'prog', ierr)
       call stochastic_physics_pattern_apply(domain_ptr, 1, tend_names, ierr)
       
    end if

  end subroutine ufs_physics_to_mpas

  !> #########################################################################################
  !> Procedure to convert of output "CCPP" variables to "MPAS" variables
  !> Called prior to MPAS dynamical core (integration)
  !>
  !> This procedure updates the MPAS "state" using prognosed physics/microphysics variables.
  !>
  !> Analogous to microphysics_to_MPAS in src/core_atmosphere/physics/mpas_atmphys_interface.F
  !>
  !> #########################################################################################
  subroutine ufs_microphysics_to_mpas(physics_state)
    use GFS_typedefs,       only : GFS_stateout_type
    use mpas_derived_types, only : mpas_pool_type
    use mpas_pool_routines, only : mpas_pool_get_subpool, mpas_pool_get_array, mpas_pool_get_dimension
    use mpas_constants,     only : gravity
    use mpas_kind_types,    only : RKIND

    ! Arguments
    type(GFS_stateout_type), intent(in   ) :: physics_state
    ! Locals
    type(mpas_statein_type) :: mpas_state
    type(mpas_pool_type), pointer :: diag_pool
    type(mpas_pool_type), pointer :: mesh_pool
    type(mpas_pool_type), pointer :: state_pool
    integer, pointer :: nCellsSolve, index_qv
    integer :: iCol
    real(kind=RKIND) :: rho1, rho2, tem1, tem2

    ! Access MPAS data pools
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'state',     state_pool)
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'diag',      diag_pool)
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh',      mesh_pool)

    ! Get MPAS dimensions
    call mpas_pool_get_dimension(mesh_pool,  'nCellsSolve', nCellsSolve)
    call mpas_pool_get_dimension(state_pool, 'index_qv',    index_qv)

    ! Grab fields from MPAS pools
    call mpas_pool_get_array(state_pool, 'scalars',                MPAS_state % tracers, timeLevel=1)
    call mpas_pool_get_array(mesh_pool,  'zgrid',                  MPAS_state % zgrid)
    call mpas_pool_get_array(mesh_pool,  'zz',                     MPAS_state % zz)
    call mpas_pool_get_array(state_pool, 'rho_zz',                 MPAS_state % rho_zz,  timeLevel=1)
    call mpas_pool_get_array(diag_pool,  'pressure_base',          MPAS_state % pressure_b)
    call mpas_pool_get_array(diag_pool,  'pressure_p',             MPAS_state % pressure_p)

    ! [i, k] -> [k, i]
    ! top-down -> bottom-up ordering convention
    ! Thermodynamic conversions from moist (CCPP) to dry (MPAS)

    ! Calculation of the surface pressure using hydrostatic assumption down to the surface.
    ! (from mpas_atmphys_interface.F:MPAS_to_physics())
    call mpas_pool_get_array(diag_pool,   'surface_pressure'      ,MPAS_state % surface_pressure)
    do iCol = 1, nCellsSolve
       tem1 = MPAS_state % zgrid(2,iCol) - MPAS_state % zgrid(1,iCol)
       tem2 = MPAS_state % zgrid(3,iCol) - MPAS_state % zgrid(2,iCol)
       rho1 = MPAS_state % rho_zz(1,iCol) * MPAS_state % zz(1,iCol) * (1. + MPAS_state % tracers(index_qv,1,iCol))
       rho2 = MPAS_state % rho_zz(2,iCol) * MPAS_state % zz(2,iCol) * (1. + MPAS_state % tracers(index_qv,2,iCol))
       MPAS_state % surface_pressure(iCol) = 0.5*gravity*(MPAS_state % zgrid(2,iCol) - MPAS_state % zgrid(1,iCol)) &
            * (rho1 - 0.5*(rho2-rho1)*tem1/(tem1+tem2))
       MPAS_state % surface_pressure(iCol) = MPAS_state % surface_pressure(iCol) + &
                                             MPAS_state % pressure_p(1,iCol) + &
                                             MPAS_state % pressure_b(1,iCol)
    enddo

    ! Housekeeping
    nullify (state_pool)
    nullify (mesh_pool)
    nullify (diag_pool)

  end subroutine ufs_microphysics_to_mpas

  !> #########################################################################################
  !> Procedure to convert of "MPAS" variables to "CCPP" variables.
  !> Called prior to CCPP Microphysics Group.
  !> 
  !> Analogous to microphysics_from_MPAS in src/core_atmosphere/physics/mpas_atmphys_interface.F
  !>
  !> This procedure accesses MPAS data using MPAS native procedures and stores the data
  !> locally in the data-containers defined above. The MPAS "state" is then translated to the
  !> CCPP "state" needed by the microphysics.
  !>
  !> #########################################################################################
  subroutine ufs_mpas_to_microphysics(physics_state)
    use GFS_typedefs,         only : GFS_statein_type
    ! Arguments
    type(GFS_statein_type),   intent(inout) :: physics_state
 
  end subroutine ufs_mpas_to_microphysics

  !> #########################################################################################
  !> Procedure to compute dry hydrostatic pressure at layer interfaces and midpoints.
  !>
  !> Given arrays of zz, zgrid, rho_zz, and theta_m from the MPAS-A prognostic state, compute
  !> dry hydrostatic pressure at layer interfaces and midpoints.
  !> The vertical dimension for 3-d arrays is innermost, and k=1 represents the lowest layer
  !> or level in the fields.
  !>
  !> \update: Dustin Swales April 2025 - Modified for use in UWM
  !>
  !> DJS to GJF: We shouldn't need this once you port the MPAS_to_physics/MPAS_to_microphysics
  !> routines from MPAS.
  !>
  !> ######################################################################################### 
  subroutine hydrostatic_pressure(nCells, nVertLevels, qsize, index_qv, zz, zgrid, rho_zz,   &
       theta_m, exner, q, pmiddry, pintdry,pmid) 
    use mpas_constants,  only: cp, rgas, cv, gravity, p0, Rv_over_Rd => rvord
    use mpas_kind_types, only: RKIND
    ! Arguments
    integer, intent(in) :: nCells
    integer, intent(in) :: nVertLevels
    integer, intent(in) :: qsize
    integer, intent(in) :: index_qv
    real(RKIND), dimension(nVertLevels, nCells),       intent(in) :: zz      ! d(zeta)/dz [-]
    real(RKIND), dimension(nVertLevels+1, nCells),     intent(in) :: zgrid   ! geometric heights of layer interfaces [m]
    real(RKIND), dimension(nVertLevels, nCells),       intent(in) :: rho_zz  ! dry density / zz [kg m^-3]
    real(RKIND), dimension(nVertLevels, nCells),       intent(in) :: theta_m ! modified potential temperature
    real(RKIND), dimension(nVertLevels, nCells),       intent(in) :: exner   ! Exner function
    real(RKIND), dimension(qsize,nVertLevels, nCells), intent(in) :: q       ! water vapor dry mixing ratio
    real(RKIND), dimension(nVertLevels, nCells),       intent(out):: pmiddry ! layer midpoint dry hydrostatic pressure [Pa]
    real(RKIND), dimension(nVertLevels+1, nCells),     intent(out):: pintdry ! layer interface dry hydrostatic pressure [Pa]
    real(RKIND), dimension(nVertLevels, nCells),       intent(out):: pmid    ! layer midpoint hydrostatic pressure [Pa]

    ! Local variables
    integer :: iCell, k, idx
    real(RKIND), dimension(nVertLevels)          :: dz       ! Geometric layer thickness in column
    real(RKIND), dimension(nVertLevels)          :: dp,dpdry ! Pressure thickness
    real(RKIND), dimension(nVertLevels+1,nCells) :: pint  ! hydrostatic pressure at interface
    real(RKIND) :: sum_water
    real(RKIND) :: pk,rhok,rhodryk,thetavk,kap1,kap2,tvk,tk
    real(RKIND), parameter :: epsilon = 0.05_RKIND
    real(RKIND) :: dp_epsilon, dpdry_epsilon

    !
    ! For each column, integrate downward from model top to compute dry hydrostatic pressure at layer
    ! midpoints and interfaces. The pressure averaged to layer midpoints should be consistent with
    ! the ideal gas law using the rho_zz and theta values prognosed by MPAS at layer midpoints.
    !
    do iCell = 1, nCells
       dz(:) = zgrid(2:nVertLevels+1,iCell) - zgrid(1:nVertLevels,iCell)
       do k = nVertLevels, 1, -1
          rhodryk  = zz(k,iCell)* rho_zz(k,iCell) !full CAM physics density
          rhok = 1.0_RKIND
          do idx=2,qsize!dry_air_species_num+1,thermodynamic_active_species_num
             rhok = rhok+q(idx,k,iCell)
          end do
          rhok     = rhok*rhodryk
          dp(k)    = gravity*dz(k)*rhok
          dpdry(k) = gravity*dz(k)*rhodryk
       end do

       k = nVertLevels
       sum_water = 1.0_RKIND
       do idx=2,qsize!dry_air_species_num+1,thermodynamic_active_species_num
          sum_water = sum_water+q(idx,k,iCell)
       end do
       rhok     = sum_water*zz(k,iCell) * rho_zz(k,iCell)
       thetavk  = theta_m(k,iCell)/sum_water
       tvk      = thetavk*exner(k,iCell)
       pk       = dp(k)*rgas*tvk/(gravity*dz(k))
       !
       ! model top pressure consistently diagnosed using the assumption that the mid level
       ! is at height z(nVertLevels-1)+0.5*dz
       !
       pintdry(nVertLevels+1,iCell) = pk-0.5_RKIND*dz(nVertLevels)*rhok*gravity  !hydrostatic
       pint   (nVertLevels+1,iCell) = pintdry(nVertLevels+1,iCell)
       do k = nVertLevels, 1, -1
          !
          ! compute hydrostatic dry interface pressure so that (pintdry(k+1)-pintdry(k))/g is pseudo density
          !
          sum_water = 1.0_RKIND
          do idx=2,qsize!dry_air_species_num+1,thermodynamic_active_species_num
             sum_water = sum_water+q(idx,k,iCell)
          end do
          thetavk = theta_m(k,iCell)/sum_water!convert modified theta to virtual theta
          tvk     = thetavk*exner(k,iCell)
          tk      = tvk*sum_water/(1.0_RKIND+Rv_over_Rd*q(index_qv,k,iCell))
          pint   (k,iCell) = pint   (k+1,iCell)+dp(k)
          pintdry(k,iCell) = pintdry(k+1,iCell)+dpdry(k)
          pmid(k,iCell)    = dp(k)   *rgas*tvk/(gravity*dz(k))
          pmiddry(k,iCell) = dpdry(k)*rgas*tk /(gravity*dz(k))
          !
          ! PMID is not necessarily bounded by the hydrostatic interface pressure.
          ! (has been found to be an issue at ~3.75km resolution in surface layer)
          !
          dp_epsilon = dp(k) * epsilon
          dpdry_epsilon = dpdry(k)*epsilon
          pmid   (k, iCell) = max(min(pmid   (k, iCell), pint   (k, iCell) - dp_epsilon), pint   (k + 1, iCell) + dp_epsilon)
          pmiddry(k, iCell) = max(min(pmiddry(k, iCell), pintdry(k, iCell) - dpdry_epsilon), pintdry(k + 1, iCell) + dpdry_epsilon)
       end do
    end do
  end subroutine hydrostatic_pressure

!> #########################################################################################
!> Procedure to transfer MPAS grid information to physics DDTs.
!>
!> #########################################################################################
  subroutine ufs_mpas_grid_to_physics(physics_grid)
    use GFS_typedefs,         only : GFS_grid_type
    use mpas_derived_types,   only : mpas_pool_type
    use mpas_pool_routines,   only : mpas_pool_get_subpool, mpas_pool_get_dimension, mpas_pool_get_array, mpas_pool_get_config
    use mpas_kind_types,      only : RKIND
    use mpas_constants,       only : pii
    use mpas_log,             only : mpas_log_write
    use mpas_derived_types,   only : MPAS_LOG_ERR, MPAS_LOG_WARN
    use mpp_mod,              only : mpp_error, FATAL
    ! Arguments
    type(GFS_grid_type),      intent(inout) :: physics_grid
    ! Locals
    type(mpas_pool_type), pointer :: mesh_pool
    integer :: i, ierr
    integer, pointer :: nCellsSolve
    real(RKIND), pointer :: lat(:), lon(:), area(:), meshDensity(:)
    
    real(RKIND), pointer :: nominalMinDc
    real(RKIND), pointer :: config_len_disp
    real(RKIND)          :: rad2deg
    
    ierr = 0
    rad2deg = 180.0_RKIND/pii
    
    ! Access MPAS data pools.
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh',  mesh_pool)
    
    ! Get MPAS dimensions
    call mpas_pool_get_dimension(mesh_pool,  'nCellsSolve', nCellsSolve)
    
    call mpas_pool_get_array(mesh_pool,  'latCell',                lat)
    call mpas_pool_get_array(mesh_pool,  'lonCell',                lon)
    call mpas_pool_get_array(mesh_pool,  'areaCell',               area)
    call mpas_pool_get_array(mesh_pool,  'meshDensity',            meshDensity)
    
    ! (from mpas_atm_core.F/atm_core_init Determine horizontal length scale used by horizontal diffusion and 3-d divergence damping
    nullify(nominalMinDc)
    call mpas_pool_get_array(mesh_pool, 'nominalMinDc', nominalMinDc)

    nullify(config_len_disp)
    call mpas_pool_get_config(domain_ptr % blocklist % configs, 'config_len_disp', config_len_disp)

    ! If config_len_disp was specified as a valid value, use that
    if (config_len_disp > 0.0_RKIND) then
      ! But if nominalMinDc was available in the input file and is different, print a warning
      if (nominalMinDc > 0.0_RKIND .and. abs(nominalMinDc - config_len_disp) > 1.0e-6_RKIND * config_len_disp) then
        call mpas_log_write('nominalMinDc was read from input file as a positive value ($r) that differs', &
                                realArgs=[nominalMinDc], messageType=MPAS_LOG_WARN)
        call mpas_log_write('from the specified config_len_disp value ($r)', &
                                realArgs=[config_len_disp], messageType=MPAS_LOG_WARN)
      end if
      nominalMinDc = config_len_disp
    ! Otherwise, try to use nominalMinDc
    else
      if (nominalMinDc > 0.0_RKIND) then
        call mpas_log_write('Setting config_len_disp to $r based on nominalMinDc value in input file', realArgs=[nominalMinDc])
          config_len_disp = nominalMinDc
      else
        call mpas_log_write('Both config_len_disp and nominalMinDc are <= 0.0.', messageType=MPAS_LOG_ERR)
        call mpas_log_write('Please either specify config_len_disp in the &nhyd_model namelist group,', &
                                messageType=MPAS_LOG_ERR)
        call mpas_log_write('or use an input file that provides a valid value for the nominalMinDc variable.', &
                                messageType=MPAS_LOG_ERR)
        ierr = 1
      end if
    end if
    if (ierr/=0)  call mpp_error(FATAL, 'Call to ufs_mpas_grid_to_physics() failed')  

    do i=1, nCellsSolve
      physics_grid % xlat(i)   = lat(i)
      physics_grid % xlon(i)   = lon(i)
      physics_grid % xlat_d(i) = physics_grid % xlat(i) * rad2deg
      physics_grid % xlon_d(i) = physics_grid % xlon(i) * rad2deg
      physics_grid % sinlat(i) = sin(physics_grid % xlat(i))
      physics_grid % coslat(i) = sqrt(1.0_RKIND - physics_grid % sinlat(i) * physics_grid % sinlat(i))
      physics_grid % area(i)   = area(i)
      !formula for dx comes from mpas_atmphys_driver_gwdo.F instead of sqrt(area) as in FV3
      physics_grid % dx(i)     = config_len_disp / meshDensity(i)**0.25
    end do
    
  end subroutine ufs_mpas_grid_to_physics
  
end module atmos_coupling_mod
