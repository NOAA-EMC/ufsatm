! ###########################################################################################
!> \file atmos_coupling.F90
!> Procedures for coupling the MPAS dynamical core to the CCPP Physics.
!>
! ###########################################################################################
module atmos_coupling_mod
  use mpas_kind_types, only : strKIND, RKIND, R8KIND
  use ufs_mpas_io,     only : domain_ptr,dyn_mpas_exchange_halo
  use ufs_mpas_io,     only : lucats,lumatch,luseas,lutype,li,albd,slmo,sfem,sfz0,therin,scfx,sfhc

  implicit none

  public :: ufs_physics_to_mpas
  public :: ufs_mpas_to_physics
  public :: ufs_microphysics_to_mpas
  public :: ufs_mpas_to_microphysics
  public :: ufs_mpas_grid_to_physics
  public :: ufs_mpas_sfc_to_physics
  public :: ufs_mpas_sfc_rad_update
  public :: ufs_mpas_landuse_update
  public :: ufs_mpas_gwd_to_physics
  public :: ufs_mpas_reference_pressure
  public :: ufs_mpas_surface_update
  public :: ufs_mpas_ozone_init

contains
  !> #########################################################################################
  !> Procedure to populate CCPP data containers with MPAS pool data.
  !> Called BEFORE CCPP Radiation and Physics Groups.
  !>
  !> Analogous to MPAS_to_physics in src/core_atmosphere/physics/mpas_atmphys_interface.F
  !>
  !> This procedure accesses MPAS data using MPAS native procedures and stores the data
  !> locally in the data-containers defined above. The MPAS "state" is then translated to the
  !> CCPP "state" needed by the physics.
  !>
  !> #########################################################################################
  subroutine ufs_mpas_to_physics(statein, stateout, surface, radiation, mpas_from_ufs_cnst)
    use GFS_typedefs,         only : GFS_sfcprop_type, GFS_radtend_type
    use GFS_typedefs,         only : GFS_statein_type, GFS_stateout_type
    use mpas_derived_types,   only : mpas_pool_type
    use mpas_pool_routines,   only : mpas_pool_get_subpool, mpas_pool_get_array, mpas_pool_get_dimension
    use atm_core,             only : atm_compute_output_diagnostics
    use mpas_constants,       only : gravity, rvord

    ! Arguments
    type(GFS_statein_type),   intent(inout) :: statein
    type(GFS_stateout_type),  intent(inout) :: stateout
    type(GFS_sfcprop_type),   intent(inout) :: surface
    type(GFS_radtend_type),   intent(inout) :: radiation
    integer, pointer,         intent(in   ) :: mpas_from_ufs_cnst(:)
    ! Locals
    type(mpas_pool_type), pointer :: state_pool, diag_pool, mesh_pool, diag_phys, sfc_input
    integer :: iCol, iLay, iTracer, ithread
    integer, pointer :: nCellsSolve, num_scalars, nVertLevels, index_qv
    integer, pointer :: nThreads, cellSolveThreadStart(:), cellSolveThreadEnd(:)
    real(kind=RKIND) :: rho1, rho2, tem1, tem2, theta
    real(kind=RKIND), pointer :: ux(:,:), uy(:,:), theta_m(:,:), rho_zz(:,:), zgrid(:,:), zz(:,:)
    real(kind=RKIND), pointer :: exner(:,:), tracers(:,:,:), pressure_b(:,:), pressure_p(:,:)
    real(kind=RKIND), pointer :: w(:,:), surface_pressure(:), rho(:,:)
    real(RKIND), pointer :: sfc_albedo(:),sfc_emiss(:)
    logical, allocatable :: per_kg(:)
    real(kind=RKIND)     :: rv0
    character(len=*), parameter :: subname = 'atmos_coupling::ufs_mpas_to_physics'

    ! Get openMP information
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'nThreads',             nThreads)
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'cellSolveThreadStart', cellSolveThreadStart)
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'cellSolveThreadEnd',   cellSolveThreadEnd)

    ! Access MPAS data pools.
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'state', state_pool)
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'diag',  diag_pool)
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh',  mesh_pool)
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'diag_physics',  diag_phys)
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'sfc_input', sfc_input)

    ! Get MPAS dimensions
    call mpas_pool_get_dimension(mesh_pool,  'nCellsSolve', nCellsSolve)
    call mpas_pool_get_dimension(mesh_pool,  'nVertLevels', nVertLevels)
    call mpas_pool_get_dimension(state_pool, 'num_scalars', num_scalars)
    call mpas_pool_get_dimension(state_pool, 'index_qv',    index_qv)

    ! Identify which scalars need conversion to dry mass or number mixing ratios.
    allocate(per_kg(num_scalars))
    call ufs_mpas_per_kg_air_mask(state_pool, num_scalars, per_kg)

    ! Grab fields from MPAS pools
    call mpas_pool_get_array(diag_pool,  'uReconstructZonal',      ux)
    call mpas_pool_get_array(diag_pool,  'uReconstructMeridional', uy)
    call mpas_pool_get_array(state_pool, 'scalars',                tracers, timeLevel=1)
    call mpas_pool_get_array(state_pool, 'w',                      w, timeLevel=1)
    call mpas_pool_get_array(diag_pool,  'exner',                  exner)
    call mpas_pool_get_array(mesh_pool,  'zgrid',                  zgrid)
    call mpas_pool_get_array(mesh_pool,  'zz',                     zz)
    call mpas_pool_get_array(state_pool, 'theta_m',                theta_m, timeLevel=1)
    call mpas_pool_get_array(state_pool, 'rho_zz',                 rho_zz,  timeLevel=1)
    call mpas_pool_get_array(diag_pool,  'pressure_base',          pressure_b)
    call mpas_pool_get_array(diag_pool,  'pressure_p',             pressure_p)
    call mpas_pool_get_array(diag_pool,  'surface_pressure',       surface_pressure)
    call mpas_pool_get_array(diag_phys,  'sfc_albedo',             sfc_albedo)
    call mpas_pool_get_array(diag_phys,  'sfc_emiss' ,             sfc_emiss)

    ! Local variables
    allocate(rho( nCellsSolve, nVertLevels))

    ! Copy fields from MPAS data containers to physics data containers.
    ! [k, i] -> [i, k]
    ! Retain bottom-up convention
    do ithread = 1,nThreads
       do iCol = cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread)
          do iLay = 1,nVertLevels
             ! Scalars (tracer,layer,col) -> (col,layer,tracer), per kg dry air -> per kg moist air
             rv0 = max(0._RKIND, tracers(index_qv,iLay,iCol))
             do iTracer = 1,num_scalars
                if (per_kg(iTracer)) then
                   statein % qgrs(iCol,iLay,mpas_from_ufs_cnst(iTracer)) = &
                        max(0._RKIND, tracers(iTracer,iLay,iCol)) / (1._RKIND + rv0)
                else
                   statein % qgrs(iCol,iLay,mpas_from_ufs_cnst(iTracer)) = max(0._RKIND, tracers(iTracer,iLay,iCol))
                end if
             end do
             
             ! Air denisty (rho) (TODO: Pass to CCPP Physics)
             rho(iCol,iLay) = zz(iLay,iCol) * rho_zz(iLay,iCol)

             ! Potential temperature (theta_m -> theta)
             theta = theta_m(iLay,iCol) / (1._RKIND + rvord * tracers(index_qv,iLay,iCol))

             ! Air temperature (theta -> temp)
             statein % tgrs(iCol,iLay)   = theta*exner(iLay,iCol)

             ! Winds at grid center
             statein % ugrs(iCol,iLay)   = ux(iLay,iCol)
             statein % vgrs(iCol,iLay)   = uy(iLay,iCol)

             ! Layer geopotential and height. CCPP geopotential is relative to the
             ! model surface (FV3: get_phi_fv3 sets phii(:,1) = 0), so subtract the
             ! terrain height zgrid(1,:). Schemes use phil/g as height above ground
             ! (e.g. drag_suite TOFD: zl**(-1.2) -> NaN where the absolute height
             ! is <= 0, i.e. land below sea level).
             statein % phil(iCol,iLay)   = (0.5*(zgrid(iLay+1,iCol)+zgrid(iLay,iCol)) - zgrid(1,iCol))*gravity !(m -> m2/s2)
             statein % zgrid(iCol,iLay)  = 0.5*(zgrid(iLay+1,iCol)+zgrid(iLay,iCol))

             ! Level geopotential (surface-relative) and height.
             statein % phii(iCol,iLay)   = (zgrid(iLay,iCol) - zgrid(1,iCol))*gravity !(m -> m2/s2)
             statein % zigrid(iCol,iLay) = zgrid(iLay,iCol)

             ! Layer thickness.
             statein % dzgrid(iCol,iLay) = zgrid(iLay+1,iCol) - zgrid(iLay,iCol)

             ! Exner funciton
             statein % prslk(iCol,iLay)  = exner(iLay,iCol)

             ! MPAS provides vertical velocity at interfaces, compute layer mean, and convert from w -> omega
             statein % vvl(iCol,iLay) =  -0.5*(w(iLay,iCol) + w(iLay+1,iCol))*rho(iCol,iLay)*gravity
          end do
          do iLay = nVertLevels,nVertLevels+1
             statein % phii(iCol,iLay)     = (zgrid(iLay,iCol) - zgrid(1,iCol))*gravity !(m -> m2/s2)
             statein % zigrid(iCol,iLay)   = zgrid(iLay,iCol)
             statein % dzgrid(iCol,iLay-1) = zgrid(iLay,iCol) - zgrid(iLay-1,iCol)
          end do

          ! Surface radiative properties
          radiation % salb(iCol)  = sfc_albedo(iCol)
          radiation % semis(iCol) = sfc_emiss(iCol)
       end do
    end do

    ! Calculation of the surface pressure using hydrostatic assumption down to the surface.
    ! (from mpas_atmphys_interface.F:MPAS_to_physics())
    do ithread = 1,nThreads
       do iCol = cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread)
          tem1 = zgrid(2,iCol) - zgrid(1,iCol)
          tem2 = zgrid(3,iCol) - zgrid(2,iCol)
          rho1 = rho_zz(1,iCol) * zz(1,iCol) * (1. + tracers(index_qv,1,iCol))
          rho2 = rho_zz(2,iCol) * zz(2,iCol) * (1. + tracers(index_qv,2,iCol))
          surface_pressure(iCol) = 0.5*gravity*(zgrid(2,iCol) - zgrid(1,iCol)) &
               * (rho1 - 0.5*(rho2-rho1)*tem1/(tem1+tem2))
          surface_pressure(iCol) = surface_pressure(iCol) + pressure_p(1,iCol) + pressure_b(1,iCol)
       end do
    end do

    ! Compute hydrostatic pressure
    call ufs_mpas_hydrostatic_pressure(statein, tracers(index_qv,:,:),'physics')

    ! DJS: Why do some process-split schemes use one temperature (air_temperature) stored
    !      in GFS_stateout, while others use one (physics_timestep_initial_air_temperature)
    !      stored in GFS_statein?
    !      Shouldn't all schemes use "air_temperature"?
    !      "air_temperature" = "physics_timestep_initial_air_temperature" until the state is
    !      updated?
    stateout % gt0(:,:)   = statein % tgrs(:,:)
    stateout % gq0(:,:,:) = statein % qgrs(:,:,:)
    stateout % gu0(:,:)   = statein % ugrs(:,:)
    stateout % gv0(:,:)   = statein % vgrs(:,:)

    ! Reset tendencies for clean accumulation of the (process-split) physics group tendencies.
    stateout % dtdt(:,:)   = 0._RKIND
    stateout % dqdt(:,:,:) = 0._RKIND
    stateout % dudt(:,:)   = 0._RKIND
    stateout % dvdt(:,:)   = 0._RKIND

    ! Housekeeping
    deallocate (rho)
    deallocate (per_kg)
    nullify (mesh_pool)
    nullify (state_pool)
    nullify (diag_pool)

  end subroutine ufs_mpas_to_physics

  !> #########################################################################################
  !> Procedure to update MPAS prognostic tendencies with physics (CCPP) tendencies.
  !> Called AFTER physics, BEFORE calling dynamics (current timestep).
  !>
  !> - Update scalar tendencies "scalars_tend" from MPAS "tend" pool
  !> - Convert from potential temerature (theta) to modified potential temperature (theta_m).
  !> - Update "theta_m" tendency from MPAS "tend" pool
  !> - Compute winds on grid edges
  !> - Update "ru" tendency from MPAS "tend" pool
  !>
  !> Analogous to phys_get_tend in physics/mpas_atmphys_todynamics.F
  !> Here, instead of updating the state with physics tendencies from the MPAS "tend_phys"
  !> pool, we use tendencies from the CCPP Physics data containers.
  !>
  !> #########################################################################################
  subroutine ufs_physics_to_mpas(stateout, mpas_from_ufs_cnst)
    use GFS_typedefs,       only : GFS_stateout_type
    use mpas_derived_types, only : mpas_pool_type
    use mpas_pool_routines, only : mpas_pool_get_subpool, mpas_pool_get_array, mpas_pool_get_dimension, mpas_pool_get_config
    use mpas_constants,     only : rv, rgas, gravity

    ! Arguments
    type(GFS_stateout_type), intent(in) :: stateout
    integer, pointer,        intent(in) :: mpas_from_ufs_cnst(:)

    ! Locals
    type(mpas_pool_type),  pointer :: state_pool, mesh_pool, tend_pool, diag_pool, tend_phys
    real(kind=RKIND), pointer :: mass(:,:), mass_edge(:,:), exner(:,:), theta_m(:,:), zgrid(:,:), zz(:,:)
    real(kind=RKIND), pointer :: pressure_b(:,:), pressure_p(:,:)
    real(kind=RKIND), pointer :: tend_u_phys(:,:)
    real(kind=RKIND), pointer :: tend_uzonal(:,:), tend_umerid(:,:), tend_rtheta_phys(:,:), pattern(:,:)
    real(kind=RKIND), pointer :: tend_rho_phys(:,:), tend_ru_phys(:,:)
    real(kind=RKIND), pointer :: scalars(:,:,:), tend_scalars(:,:,:)
    real(kind=RKIND), pointer :: surface_pressure(:)
    integer, pointer :: nCells, nCellsSolve, num_scalars, nVertLevels, nEdges, nEdgesSolve
    integer, pointer :: index_qv => null()
    integer, pointer :: index_qc => null()
    integer, pointer :: index_qi => null()
    integer, pointer :: index_qr => null()
    integer, pointer :: index_qs => null()
    integer, pointer :: index_qg => null()
    integer, pointer :: index_nc => null()
    integer, pointer :: index_ni => null()
    integer, pointer :: index_nifa => null()
    integer, pointer :: index_nwfa => null()
    integer, pointer :: index_sgs_tke => null()

    integer, pointer :: nThreads, cellSolveThreadStart(:), cellSolveThreadEnd(:)
    integer :: iCol,iLay,ithread,iTracer
    real(kind=RKIND):: coeff, tem1, tem2, rho1, rho2, tend_th_phys
    logical, allocatable :: per_kg(:)
    logical, pointer :: do_sppt
    real(kind=RKIND), pointer :: config_dt
    integer :: iq_v
    real(kind=8) :: rv0, qv0, qv1, dqv, dq, qx0, den, drdt
    character(len=*), parameter :: subname = 'atmos_coupling::ufs_mpas_physics_to_mpas'

    ! Get openMP information
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'nThreads',             nThreads)
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'cellSolveThreadStart', cellSolveThreadStart)
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'cellSolveThreadEnd',   cellSolveThreadEnd)
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'nEdges',               nEdges)
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'nEdgesSolve',          nEdgesSolve)

    ! Access MPAS data pools
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'state',        state_pool)
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh',         mesh_pool)
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'tend',         tend_pool)
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'diag',         diag_pool)
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'tend_physics', tend_phys)

    ! Get MPAS dimensions
    call mpas_pool_get_dimension(mesh_pool,  'nCellsSolve', nCellsSolve)
    call mpas_pool_get_dimension(mesh_pool,  'nVertLevels', nVertLevels)
    call mpas_pool_get_dimension(mesh_pool,  'nCells',      nCells)
    call mpas_pool_get_dimension(state_pool, 'num_scalars', num_scalars)
    call mpas_pool_get_dimension(state_pool, 'index_qv',    index_qv)
    ! DJS: In DEBUG mode, MPAS will report a warning for any of these fields that haven't
    !      been set. This can be ignored. Below we check to see if any of the fields are
    !      associated before trying to access them.
    call mpas_pool_get_dimension(state_pool, 'index_qc',    index_qc)
    call mpas_pool_get_dimension(state_pool, 'index_qi',    index_qi)
    call mpas_pool_get_dimension(state_pool, 'index_qr',    index_qr)
    call mpas_pool_get_dimension(state_pool, 'index_qs',    index_qs)
    call mpas_pool_get_dimension(state_pool, 'index_qg',    index_qg)
    call mpas_pool_get_dimension(state_pool, 'index_nc',    index_nc)
    call mpas_pool_get_dimension(state_pool, 'index_ni',    index_ni)
    call mpas_pool_get_dimension(state_pool, 'index_nifa',  index_nifa)
    call mpas_pool_get_dimension(state_pool, 'index_nwfa',  index_nwfa)
    call mpas_pool_get_dimension(state_pool, 'index_sgs_tke', index_sgs_tke)

    call mpas_pool_get_config(domain_ptr % blocklist % configs, 'config_dt', config_dt)

    ! Identify which scalars need conversion to dry mass or number mixing ratios.
    allocate(per_kg(num_scalars))
    call ufs_mpas_per_kg_air_mask(state_pool, num_scalars, per_kg)
    iq_v = mpas_from_ufs_cnst(index_qv)

    ! Grab fields from MPAS pools
    call mpas_pool_get_array(state_pool,'theta_m',          theta_m,1)
    call mpas_pool_get_array(state_pool,'scalars',          scalars,1)
    call mpas_pool_get_array(state_pool,'rho_zz',           mass,   1)
    call mpas_pool_get_array(mesh_pool, 'zz',               zz)
    call mpas_pool_get_array(mesh_pool, 'zgrid',            zgrid)
    call mpas_pool_get_array(diag_pool, 'rho_edge',         mass_edge)
    call mpas_pool_get_array(diag_pool, 'surface_pressure', surface_pressure)
    call mpas_pool_get_array(diag_pool, 'pressure_base',    pressure_b)
    call mpas_pool_get_array(diag_pool, 'pressure_p',       pressure_p)
    call mpas_pool_get_array(diag_pool, 'exner',            exner)
    ! Physics tendencies used by MPAS during time integration.
    call mpas_pool_get_array(tend_phys, 'tend_uzonal',         tend_uzonal)
    call mpas_pool_get_array(tend_phys, 'tend_umerid',         tend_umerid)
    call mpas_pool_get_array(tend_phys, 'tend_rtheta_physics', tend_rtheta_phys)
    call mpas_pool_get_array(tend_phys, 'tend_rho_physics',    tend_rho_phys)
    call mpas_pool_get_array(tend_phys, 'tend_ru_physics',     tend_ru_phys)
    call mpas_pool_get_array(tend_pool, 'scalars_tend',        tend_scalars)

    ! Allocate/initialize local variables
    allocate(tend_u_phys(nVertLevels,nEdges+1))
    tend_u_phys(:,:) = 0._RKIND

    ! Initalize physics tendencies in MPAS.
    tend_scalars(:,:,:)   = 0._RKIND
    tend_rtheta_phys(:,:) = 0._RKIND
    tend_rho_phys(:,:)    = 0._RKIND
    tend_ru_phys(:,:)     = 0._RKIND

    !> #####################################################################################
    !> 0) Stochastic physics (physics processes)
    !> #####################################################################################
    call mpas_pool_get_config(domain_ptr % blocklist % configs, 'do_sppt', do_sppt)

    ! Call stochastic physics (SPPT) for (hydrostatic) physics tendencies.
    if (do_sppt) then
       ! Grab stochastic physics pattern.
       call mpas_pool_get_array(tend_phys,'stoch_pattern_sppt',pattern)

       ! Apply pattern to scheme tendencies
       do ithread = 1,nThreads
          do iCol = cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread)
             do iLay = 1,nVertLevels
                stateout % ten_u_conv(iCol,iLay) = stateout % ten_u_conv(iCol,iLay) * pattern(ilay,iCol)
                stateout % ten_v_conv(iCol,iLay) = stateout % ten_v_conv(iCol,iLay) * pattern(ilay,iCol)
                stateout % ten_u_pbl(iCol,iLay)  = stateout % ten_u_pbl(iCol,iLay)  * pattern(ilay,iCol)
                stateout % ten_v_pbl(iCol,iLay)  = stateout % ten_v_pbl(iCol,iLay)  * pattern(ilay,iCol)
             end do
          end do
       end do

    end if

    !> #####################################################################################
    !> 1) Update MPAS tendency "tend_scalars".
    !> #####################################################################################

    do ithread=1,nThreads
      do iCol=cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread)
        do iLay = 1,nVertLevels
          ! State the physics saw (per kg moist air) and its vapor tendency
          rv0 = max(0._RKIND, scalars(index_qv,iLay,iCol))
          qv0 = rv0 / (1._R8KIND + rv0)
          dqv = stateout % dqdt(iCol,iLay,iq_v)
          qv1 = qv0 + config_dt*dqv
          den = (1._R8KIND - qv1) * (1._R8KIND - qv0)
          do iTracer = 1, num_scalars
            dq = stateout % dqdt(iCol,iLay,mpas_from_ufs_cnst(iTracer))
            if (per_kg(iTracer)) then
              qx0  = max(0._RKIND, scalars(iTracer,iLay,iCol)) / (1._R8KIND + rv0)
              drdt = (dq*(1._R8KIND - qv0) + qx0*dqv) / den      ! d(mixing ratio)/dt
            else
              drdt = dq                                         ! e.g. sgs_tke
              ! TKE is updated in place by the PBL scheme (no tendency), so form one here.
              ! DJS2026: MYNN-EDMF outputs the updated SGS TKE scaled by 2x. Is this needed
              ! in MYNN?
              if (associated(index_sgs_tke)) then
                 if (iTracer == index_sgs_tke) then
                    drdt = (stateout % gq0(iCol,iLay,mpas_from_ufs_cnst(iTracer)) - &
                            max(0._RKIND, scalars(iTracer,iLay,iCol))) / config_dt
                 end if
              end if
            end if
            tend_scalars(iTracer,iLay,iCol) = tend_scalars(iTracer,iLay,iCol) + &
                 real(drdt, RKIND) * mass(iLay,iCol)
          end do
        end do
      end do
    end do
    
    ! Update halo points.
    call dyn_mpas_exchange_halo('scalars_tend',.true.)

    !> #####################################################################################
    !> 2) Update MPAS tendency "tend_rtheta_phys".
    !> #####################################################################################
    !
    do ithread = 1,nThreads
       do iCol = cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread)
          do iLay = 1, nVertLevels
             ! Add accumulated temperature tendencies from the physics group, convert from temperature
             ! to rho*potential temperature (T -> rtheta). 
             tend_th_phys = (stateout % dtdt(iCol,iLay)/exner(iLay,iCol))*mass(iLay,iCol)
             coeff = (1. + rv/rgas * scalars(index_qv,iLay,iCol))
             ! Convert from the tendency of potential temperature to the tendency of the  modified
             ! potential temperature (rtheta -> rtheta_m). 
             tend_rtheta_phys(iLay,iCol) = coeff * tend_th_phys + rv/rgas * theta_m(iLay,iCol) * tend_scalars(index_qv,iLay,iCol) / coeff
          end do
       end do
    end do

    ! Update halo points.
    call dyn_mpas_exchange_halo('tend_rtheta_physics',.true.)

    !> #####################################################################################
    !> 3) Update MPAS tendency "tend_ru_phys".
    !> #####################################################################################

    ! First, need to regrid from grid-centers (physics) to grid-edges (dynamics)...
    ! Copy CCPP accumulated physics tendencies to MPAS scratch arrays.
    do ithread = 1,nThreads
       do iCol = cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread)
           do iLay = 1, nVertLevels
              tend_uzonal(iLay,iCol) = stateout % dudt(iCol,iLay)
              tend_umerid(iLay,iCol) = stateout % dvdt(iCol,iLay)
           end do
       end do
    end do

    ! Next, update the halo points of the scratch arrays.
    call dyn_mpas_exchange_halo('tend_uzonal',.true.)
    call dyn_mpas_exchange_halo('tend_umerid',.true.)

    ! Finally, compute wind tendency at grid-edges.
    call tend_toEdges(mesh_pool, tend_uzonal, tend_umerid, tend_u_phys)

    ! Update MPAS tendency (ru)
    do iCol = 1,nEdgesSolve
       do iLay = 1, nVertLevels
          tend_ru_phys(iLay,iCol) = tend_u_phys(iLay,iCol)*mass_edge(iLay,iCol)
       end do
    end do

    ! Update halo points.
    call dyn_mpas_exchange_halo('tend_ru_physics',.true.)

    !> #####################################################################################
    !> 4) Stochastic physics (non-hydrostatic physics)
    !> #####################################################################################
    ! Call stochastic physics (SPPT) for (non-hydrostatic) physics tendencies.
    if (do_sppt) then
       ! Grab stochastic physics pattern.
       call mpas_pool_get_array(tend_phys,'stoch_pattern_sppt',pattern)

       ! Apply pattern to accumulated tendencies.
    end if
    
    !> #####################################################################################
    !> Diagnostics
    !> #####################################################################################
    ! Calculation of the surface pressure using hydrostatic assumption down to the surface.
    ! (from mpas_atmphys_interface.F:MPAS_to_physics())
    do ithread = 1,nThreads
       do iCol = cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread)
          tem1 = zgrid(2,iCol) - zgrid(1,iCol)
          tem2 = zgrid(3,iCol) - zgrid(2,iCol)
          rho1 = mass(1,iCol) * zz(1,iCol) * (1. + scalars(index_qv,1,iCol))
          rho2 = mass(2,iCol) * zz(2,iCol) * (1. + scalars(index_qv,2,iCol))
          surface_pressure(iCol) = 0.5*gravity*( zgrid(2,iCol) -  zgrid(1,iCol)) &
               * (rho1 - 0.5*(rho2-rho1)*tem1/(tem1+tem2))
          surface_pressure(iCol) = surface_pressure(iCol) + pressure_p(1,iCol) + &
               pressure_b(1,iCol)
       enddo
    enddo

    ! Housekeeping
    deallocate(per_kg)
    deallocate(tend_u_phys)
    nullify (state_pool)
    nullify (mesh_pool)
    nullify (diag_pool)
    nullify (tend_phys)

  end subroutine ufs_physics_to_mpas

  !> #########################################################################################
  !> Procedure to compute diabatic heating tendency from microphysics and store in MPAS pool.
  !> Called AFTER microphysics, BEFORE calling dynamics (next timestep)
  !>
  !> This routine is where we compute the diabatic heating rate, "rt_diabatic_tend", due to the
  !> microphysics.
  !> Follows microphysics_to_MPAS from src/core_atmosphere/physics/mpas_atmphys_interface.F.
  !> compute the diabatic heating rate, "rt_diabatic_tend" and save in MPAS memory for use by
  !> the dynamics at the subsequent time step in dynamics/mpas_atm_time_integration.F.
  !> Additionally, update any other fields needed by the dynamics (e.g., theta_m, rtheta_p)
  !>
  !> #########################################################################################
  subroutine ufs_microphysics_to_mpas(stateout, mpas_from_ufs_cnst)
    use GFS_typedefs,       only : GFS_stateout_type
    use mpas_derived_types, only : mpas_pool_type
    use mpas_pool_routines, only : mpas_pool_get_subpool, mpas_pool_get_array
    use mpas_pool_routines, only : mpas_pool_get_dimension, mpas_pool_get_config
    use mpas_constants,     only : gravity, rvord, rv, rgas, p0, cp
    use mpas_atm_boundaries, only : nRelaxZone

    ! Arguments
    type(GFS_stateout_type), intent(in) :: stateout
    integer, pointer,        intent(in) :: mpas_from_ufs_cnst(:)

    ! Locals
    type(mpas_pool_type), pointer :: diag_pool, mesh_pool, state_pool, tend_pool
    integer :: iCol, ithread, iLay, iTracer
    integer, pointer :: nCellsSolve
    integer, pointer :: index_qv => null()
    integer, pointer :: num_scalars, nVertLevels, bdyMaskCell(:)
    integer, pointer :: nThreads, cellSolveThreadStart(:), cellSolveThreadEnd(:)
    real(kind=RKIND) :: rho1, rho2, tem1, tem2, coeff, rcv, theta_dyn, qv_old
    real(kind=RKIND), pointer :: config_dt
    real(kind=RKIND), pointer :: tracers(:,:,:), rt_diabatic_tend(:,:), rho_zz(:,:), theta_m(:,:)
    real(kind=RKIND), pointer :: zz(:,:), zgrid(:,:), exner(:,:), exner_b(:,:), rtheta_b(:,:), theta(:,:)
    real(kind=RKIND), pointer :: rtheta_p(:,:), pressure_b(:,:), pressure_p(:,:), surface_pressure(:), dtheta_dt_mp(:,:)
    logical, allocatable :: per_kg(:)
    integer :: iq_v
    real(kind=8) :: rv0, qv0, qv1, dqv, dq, qx0, den, drdt
    character(len=*), parameter :: subname = 'atmos_coupling::ufs_microphysics_to_mpas'

    ! Get openMP information
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'nThreads',             nThreads)
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'cellSolveThreadStart', cellSolveThreadStart)
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'cellSolveThreadEnd',   cellSolveThreadEnd)

    ! Access MPAS data pools
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'state', state_pool)
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'diag',  diag_pool)
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh',  mesh_pool)
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'tend',  tend_pool)

    ! Model timestep
    call mpas_pool_get_config( domain_ptr % blocklist % configs, 'config_dt', config_dt)

    ! MPAS dimensions
    call mpas_pool_get_dimension(mesh_pool,  'nCellsSolve', nCellsSolve)
    call mpas_pool_get_dimension(state_pool, 'index_qv',    index_qv)
    call mpas_pool_get_dimension(state_pool, 'num_scalars', num_scalars)
    call mpas_pool_get_dimension(mesh_pool,  'nVertLevels', nVertLevels)

    ! Identify which scalars need conversion to dry mass or number mixing ratios.
    allocate(per_kg(num_scalars))
    call ufs_mpas_per_kg_air_mask(state_pool, num_scalars, per_kg)
    iq_v = mpas_from_ufs_cnst(index_qv)

    ! Grab fields from MPAS pools
    call mpas_pool_get_array(state_pool, 'scalars',          tracers, timeLevel=1)
    call mpas_pool_get_array(state_pool, 'rho_zz',           rho_zz,  timeLevel=1)
    call mpas_pool_get_array(state_pool, 'theta_m',          theta_m, timeLevel=1)
    call mpas_pool_get_array(mesh_pool,  'zgrid',            zgrid)
    call mpas_pool_get_array(mesh_pool,  'zz',               zz)
    call mpas_pool_get_array(mesh_pool,  'bdyMaskCell',      bdyMaskCell)
    call mpas_pool_get_array(diag_pool,  'pressure_base',    pressure_b)
    call mpas_pool_get_array(diag_pool,  'pressure_p',       pressure_p)
    call mpas_pool_get_array(diag_pool,  'surface_pressure', surface_pressure)
    call mpas_pool_get_array(diag_pool,  'exner',            exner)
    call mpas_pool_get_array(diag_pool,  'exner_base',       exner_b)
    call mpas_pool_get_array(diag_pool,  'rtheta_p',         rtheta_p)
    call mpas_pool_get_array(diag_pool,  'rtheta_base',      rtheta_b)
    call mpas_pool_get_array(tend_pool,  'rt_diabatic_tend', rt_diabatic_tend)
    call mpas_pool_get_array(diag_pool,  'dtheta_dt_mp',     dtheta_dt_mp)

    ! The MPAS version of microphysics schemes update the state within the dynamics;
    ! for CCPP/UFS, we will need to update the state variables here for use by the dynamics.
    ! Also, Update water vapor, cloud liquid water, rain mixing ratios, modified potential temperature,
    ! and potential temperature heating rate from microphysics
    allocate(theta(nVertLevels, nCellsSolve))
    rcv = rgas/(cp-rgas)
    do ithread=1,nThreads
       do iCol=cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread)
        ! Limited-area: cells in the specified zone are set from the LBC data by the dycore.
         ! Leave them untouched by microphysics
         !(standalone MPAS resets them after microphysics, from the dycore.).
        if (bdyMaskCell(iCol) > nRelaxZone) then
           rt_diabatic_tend(:,iCol) = 0._RKIND
           dtheta_dt_mp(:,iCol)     = 0._RKIND
           cycle
        end if
        do iLay = 1,nVertLevels

          ! Initialize diabatic heating tendency (initially fill with theta_m before updating)
          rt_diabatic_tend(iLay,iCol) = theta_m(iLay,iCol)

          ! Update potential temperature (theta) with microphysics tendency
          coeff = (1._RKIND + rvord * tracers(index_qv,iLay,iCol))
          theta_dyn = theta_m(ilay,iCol)/coeff
          theta(iLay,iCol) = theta_dyn + config_dt * (stateout % dtdt(iCol,iLay) / exner(iLay,iCol))
          qv_old = tracers(index_qv,iLay,iCol)
          ! Old state as the microphysics saw it; must be taken BEFORE any tracer is updated
          rv0 = max(0._RKIND, tracers(index_qv,iLay,iCol))
          qv0 = rv0 / (1._R8KIND + rv0)
          dqv = stateout % dqdt(iCol,iLay,iq_v)
          qv1 = qv0 + config_dt*dqv
          den = (1._R8KIND - qv1) * (1._R8KIND - qv0)
          do iTracer = 1,num_scalars
            dq = stateout % dqdt(iCol,iLay,mpas_from_ufs_cnst(iTracer))
            if (per_kg(iTracer)) then
              qx0  = max(0._RKIND, tracers(iTracer,iLay,iCol)) / (1._R8KIND + rv0)
              drdt = (dq*(1._R8KIND - qv0) + qx0*dqv) / den
            else
              drdt = dq
            end if
            tracers(iTracer,iLay,iCol) = max(0._RKIND, tracers(iTracer,iLay,iCol) + real(config_dt*drdt, RKIND))
          end do

          ! Modified potential temperature (theta ->theta_m)
          coeff = (1._RKIND + rvord * tracers(index_qv,iLay,iCol))
          theta_m(iLay,iCol) = theta(iLay,iCol)*coeff

          ! Now compute diabatic heating due to microphsyics, save for next time step
          rt_diabatic_tend(iLay,iCol) = (theta_m(iLay,iCol) - rt_diabatic_tend(iLay,iCol)) / config_dt

          ! Save the straight theta tendency due to microphysics
          dtheta_dt_mp(iLay,iCol) =  (theta(iLay,iCol) - theta_dyn) / config_dt

          ! Density weighted perturbation potential temperature
          rtheta_p(iLay,iCol) = rho_zz(iLay,iCol) * theta_m(iLay,iCol) - rtheta_b(iLay,iCol)

          ! Exner function
          exner(iLay,iCol) = (zz(iLay,iCol)*(rgas/P0)*(rtheta_p(iLay,iCol)+rtheta_b(iLay,iCol)))**rcv

          ! Perturbation pressure
          pressure_p(iLay,iCol) = zz(iLay,iCol)*rgas*(exner(iLay,iCol)*rtheta_p(iLay,iCol) + &
                                    (exner(iLay,iCol)-exner_b(iLay,iCol))*rtheta_b(iLay,iCol))
        end do
      end do
    end do

    ! Calculation of the surface pressure using hydrostatic assumption down to the surface.
    ! (from mpas_atmphys_interface.F:MPAS_to_physics())
    do ithread = 1,nThreads
       do iCol = cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread)
          tem1 = zgrid(2,iCol) - zgrid(1,iCol)
          tem2 = zgrid(3,iCol) - zgrid(2,iCol)
          rho1 = rho_zz(1,iCol) * zz(1,iCol) * (1. + tracers(index_qv,1,iCol))
          rho2 = rho_zz(2,iCol) * zz(2,iCol) * (1. + tracers(index_qv,2,iCol))
          surface_pressure(iCol) = 0.5*gravity*(zgrid(2,iCol) - zgrid(1,iCol)) &
               * (rho1 - 0.5*(rho2-rho1)*tem1/(tem1+tem2))
          surface_pressure(iCol) = surface_pressure(iCol) + pressure_p(1,iCol) + pressure_b(1,iCol)
       end do
    end do

    ! Housekeeping
    deallocate (per_kg)
    nullify (state_pool)
    nullify (mesh_pool)
    nullify (diag_pool)

  end subroutine ufs_microphysics_to_mpas

  !> #########################################################################################
  !> Procedure to update physics (CCPP) state using updated MPAS state.
  !> Called AFTER dynamics, BEFORE microphysics.
  !>
  !> Analogous to microphysics_from_MPAS in src/core_atmosphere/physics/mpas_atmphys_interface.F
  !>
  !> #########################################################################################
  subroutine ufs_mpas_to_microphysics(stateout, statein, mpas_from_ufs_cnst)
    use GFS_typedefs,         only : GFS_stateout_type, GFS_statein_type
    use mpas_derived_types,   only : mpas_pool_type
    use mpas_pool_routines,   only : mpas_pool_get_subpool, mpas_pool_get_array, mpas_pool_get_dimension
    use mpas_constants,       only : rvord, gravity

    ! Arguments
    type(GFS_stateout_type), intent(inout) :: stateout
    type(GFS_statein_type),  intent(inout) :: statein
    integer, pointer,        intent(in   ) :: mpas_from_ufs_cnst(:)

    ! Locals
    type(mpas_pool_type), pointer :: state_pool, diag_pool, mesh_pool
    integer :: ithread, iCol, iTracer, iLay
    integer, pointer :: num_scalars, nVertLevels, nCellsSolve, index_qv
    integer, pointer :: nThreads, cellSolveThreadStart(:), cellSolveThreadEnd(:)
    real(kind=RKIND), pointer :: rho_zz(:,:), theta_m(:,:), zz(:,:), zgrid(:,:), exner(:,:)
    real(kind=RKIND), pointer :: tracers(:,:,:), w(:,:), pressure_b(:,:), pressure_p(:,:)
    real(kind=RKIND) :: theta, rho
    logical, allocatable :: per_kg(:)
    real(kind=RKIND)     :: rv0
    character(len=*), parameter :: subname = 'atmos_coupling::ufs_mpas_to_microphysics'

    ! Get openMP information
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'nThreads',             nThreads)
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'cellSolveThreadStart', cellSolveThreadStart)
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'cellSolveThreadEnd',   cellSolveThreadEnd)

    ! Access MPAS data pools
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'state', state_pool)
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'diag',  diag_pool)
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh',  mesh_pool)

    ! MPAS dimensions
    call mpas_pool_get_dimension(mesh_pool,  'nVertLevels', nVertLevels)
    call mpas_pool_get_dimension(mesh_pool,  'nCellsSolve', nCellsSolve)
    call mpas_pool_get_dimension(state_pool, 'num_scalars', num_scalars)
    call mpas_pool_get_dimension(state_pool, 'index_qv',    index_qv)

    ! Identify which scalars need conversion to dry mass or number mixing ratios.
    allocate(per_kg(num_scalars))
    call ufs_mpas_per_kg_air_mask(state_pool, num_scalars, per_kg)

    ! Grab fields from MPAS pools
    call mpas_pool_get_array(state_pool, 'rho_zz',  rho_zz,  timeLevel=1)
    call mpas_pool_get_array(state_pool, 'theta_m', theta_m, timeLevel=1)
    call mpas_pool_get_array(state_pool, 'scalars', tracers, timeLevel=1)
    call mpas_pool_get_array(state_pool, 'w',       w,       timeLevel=1)
    call mpas_pool_get_array(mesh_pool,  'zz',      zz)
    call mpas_pool_get_array(mesh_pool,  'zgrid',   zgrid)
    call mpas_pool_get_array(diag_pool,  'exner',   exner)
    call mpas_pool_get_array(diag_pool,  'pressure_base',    pressure_b)
    call mpas_pool_get_array(diag_pool,  'pressure_p',       pressure_p)

    ! Update fields needed by microphysics...
    do ithread = 1,nThreads
       do iCol = cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread)
          do iLay = 1,nVertLevels
             ! Scalars (tracer,layer,col) -> (col,layer,tracer)
             rv0 = max(0._RKIND, tracers(index_qv,iLay,iCol))
             do iTracer = 1,num_scalars
                if (per_kg(iTracer)) then
                   stateout % gq0(iCol,iLay,mpas_from_ufs_cnst(iTracer)) = &
                        max(0._RKIND, tracers(iTracer,iLay,iCol)) / (1._RKIND + rv0)
                else
                   stateout % gq0(iCol,iLay,mpas_from_ufs_cnst(iTracer)) = max(0._RKIND, tracers(iTracer,iLay,iCol))
                end if
             end do

             ! Potential temperature (theta_m -> theta)
             theta = theta_m(iLay,iCol) / (1._RKIND + rvord * max(0._RKIND,tracers(index_qv,iLay,iCol)))

             ! Air temperature (theta -> t)
             stateout % gt0(iCol,iLay) = theta*exner(iLay,iCol)

             ! Vertical velocity (w -> omega)
             rho = zz(iLay,iCol) * rho_zz(iLay,iCol)
             statein % vvl(iCol,iLay) = -0.5*(w(iLay,iCol) + w(iLay+1,iCol))*rho*gravity

             ! Non-hydrostatic pressure
             statein % prsl(iCol,iLay) = pressure_p(iLay,iCol) + pressure_b(iLay,iCol)
             statein % prsi(iCol,iLay) = exner(iLay,iCol)
          end do
       end do
    end do

    ! Update hydrostatic pressure. NO. MP in MPAS uses full non-hydrostatic pressure.    
    call ufs_mpas_hydrostatic_pressure(statein, tracers(index_qv,:,:),'microphysics')

    ! Reset tendencies for clean accumulation of ONLY the microphysics tendencies.
    stateout % dtdt(:,:)   = 0._RKIND
    stateout % dqdt(:,:,:) = 0._RKIND
    stateout % dudt(:,:)   = 0._RKIND
    stateout % dvdt(:,:)   = 0._RKIND

    ! Houzekeeping
    deallocate(per_kg)
    nullify(diag_pool)
    nullify(mesh_pool)
    nullify(state_pool)

  end subroutine ufs_mpas_to_microphysics

!> #########################################################################################
!> Procedure to transfer MPAS grid information to physics DDTs.
!>
!> #########################################################################################
  subroutine ufs_mpas_grid_to_physics(physics_grid)
    use GFS_typedefs,         only : GFS_grid_type
    use mpas_derived_types,   only : mpas_pool_type
    use mpas_pool_routines,   only : mpas_pool_get_subpool, mpas_pool_get_dimension
    use mpas_pool_routines,   only : mpas_pool_get_array, mpas_pool_get_config
    use mpas_constants,       only : pii
    use mpas_log,             only : mpas_log_write
    use mpas_derived_types,   only : MPAS_LOG_ERR, MPAS_LOG_WARN, MPAS_LOG_CRIT
    ! Arguments
    type(GFS_grid_type),      intent(inout) :: physics_grid

    ! Locals
    type(mpas_pool_type), pointer :: mesh_pool
    integer :: i, ierr, ithread
    integer, pointer :: nThreads, cellSolveThreadStart(:), cellSolveThreadEnd(:)
    real(RKIND), pointer :: lat(:), lon(:), area(:), meshDensity(:)
    real(RKIND), pointer :: nominalMinDc, config_len_disp
    real(RKIND)          :: rad2deg
    character(len=*), parameter :: subname = 'atmos_coupling::ufs_mpas_grid_to_physics'

    ierr = 0
    rad2deg = 180.0_RKIND/pii

    ! Get openMP information
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'nThreads',             nThreads)
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'cellSolveThreadStart', cellSolveThreadStart)
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'cellSolveThreadEnd',   cellSolveThreadEnd)

    ! Access MPAS data pools.
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh',  mesh_pool)

    ! Grab fields from MPAS pools
    call mpas_pool_get_array(mesh_pool,  'latCell',     lat)
    call mpas_pool_get_array(mesh_pool,  'lonCell',     lon)
    call mpas_pool_get_array(mesh_pool,  'areaCell',    area)
    call mpas_pool_get_array(mesh_pool,  'meshDensity', meshDensity)

    ! (from mpas_atm_core.F/atm_core_init Determine horizontal length scale used by horizontal diffusion and 3-d divergence damping
    nullify(nominalMinDc)
    call mpas_pool_get_array(mesh_pool, 'nominalMinDc', nominalMinDc)

    nullify(config_len_disp)
    call mpas_pool_get_config(domain_ptr % blocklist % configs, 'config_len_disp', config_len_disp)

    ! If config_len_disp was specified as a valid value, use that
    if (config_len_disp > 0.0_RKIND) then
      ! But if nominalMinDc was available in the input file and is different, print a warning
      if (nominalMinDc > 0.0_RKIND .and. abs(nominalMinDc - config_len_disp) > 1.0e-6_RKIND * config_len_disp) then
        call mpas_log_write(subname // ' WARNING: nominalMinDc was read from input file as a positive value ($r) that differs', &
                            realArgs=[nominalMinDc], messageType=MPAS_LOG_WARN)
        call mpas_log_write(subname // ' WARNING: from the specified config_len_disp value ($r)', &
                            realArgs=[config_len_disp], messageType=MPAS_LOG_WARN)
      end if
      nominalMinDc = config_len_disp
    ! Otherwise, try to use nominalMinDc
    else
      if (nominalMinDc > 0.0_RKIND) then
        call mpas_log_write('Setting config_len_disp to $r based on nominalMinDc value in input file', realArgs=[nominalMinDc])
          config_len_disp = nominalMinDc
      else
         call mpas_log_write(subname // ' ERROR: Both config_len_disp and nominalMinDc are <= 0.0.', &
                             messageType=MPAS_LOG_ERR)
         call mpas_log_write(subname // ' ERROR: Please either specify config_len_disp in the &nhyd_model namelist group,', &
                             messageType=MPAS_LOG_ERR)
         call mpas_log_write(subname // ' ERROR: or use an input file that provides a valid value for the nominalMinDc variable.', &
                             messageType=MPAS_LOG_ERR)
        ierr = 1
      end if
    end if
    if (ierr/=0)  call mpas_log_write(subname // ' ERROR: Call to ufs_mpas_grid_to_physics() failed', messageType=MPAS_LOG_CRIT)

    do ithread = 1,nThreads
       do i = cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread)
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
    end do

  end subroutine ufs_mpas_grid_to_physics

!> #########################################################################################
!> Procedure to transfer MPAS information to physics srfprop DDT
!>
!> #########################################################################################
  subroutine ufs_mpas_sfc_to_physics(surface, control)
    use GFS_typedefs,         only : GFS_sfcprop_type, GFS_control_type
    use mpas_derived_types,   only : mpas_pool_type
    use mpas_pool_routines,   only : mpas_pool_get_subpool, mpas_pool_get_dimension, mpas_pool_get_array
    use mpas_kind_types,      only : RKIND, StrKIND

    ! Arguments
    type(GFS_control_type),      intent(in) :: control
    type(GFS_sfcprop_type),      intent(inout) :: surface
    ! Locals
    type(mpas_pool_type), pointer :: sfc_input, mesh, diag_phys
    integer :: iCol, iLev, ithread
    integer, pointer :: nThreads, cellSolveThreadStart(:), cellSolveThreadEnd(:)
    integer, pointer :: isltyp(:), ivgtyp(:), landmask(:), nSoilLevels
    real(RKIND), pointer :: dzs(:,:), sh2o(:,:), smois(:,:), tslb(:,:)
    real(RKIND), pointer :: albbck(:), skintemp(:), snow(:), snowc(:), snowh(:)
    real(RKIND), pointer :: sst(:), tmn(:), vegfra(:), seaice(:), xice(:), xland(:), znt(:), sfc_albedo(:), canwat(:)
    real(RKIND), pointer :: greenfrac(:,:), albedo12m(:,:), landusef(:,:), soilf(:,:), sfc_emiss(:)
    real(RKIND), pointer :: ter(:), shdmin(:), shdmax(:), snoalb(:)
    character(len=StrKIND), pointer :: mminlu
    character(len=*), parameter :: subname = 'atmos_coupling::ufs_mpas_sfc_to_physics'

    ! Get openMP information
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'nThreads',             nThreads)
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'cellSolveThreadStart', cellSolveThreadStart)
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'cellSolveThreadEnd',   cellSolveThreadEnd)
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'nSoilLevels',          nSoilLevels)

    ! Access MPAS data pools.
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'sfc_input', sfc_input)
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh', mesh)
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'diag_physics', diag_phys)

    ! Grab fields from MPAS pools    
    call mpas_pool_get_array(sfc_input, 'isltyp',    isltyp)    !dim (nCells); dominant soil category
    call mpas_pool_get_array(sfc_input, 'ivgtyp',    ivgtyp)    !dim (nCells); dominant vegetation category
    call mpas_pool_get_array(sfc_input, 'landmask',  landmask)  !dim (nCells); land-ocean mask (1=land ; 0=ocean)
    call mpas_pool_get_array(sfc_input, 'mminlu',    mminlu)    ! (string) land use classification
    call mpas_pool_get_array(sfc_input, 'sfc_albbck',albbck)    !dim (nCells Time); background surface albedo
    call mpas_pool_get_array(sfc_input, 'sh2o',      sh2o)      !dim (nSoilLevels nCells Time); soil equivalent liquid water (m3 m^{-3})
    call mpas_pool_get_array(sfc_input, 'smois',     smois)     !dim (nSoilLevels nCells Time); soil moisture (m3 m^{-3})
    call mpas_pool_get_array(sfc_input, 'skintemp',  skintemp)  !dim (nCells Time); ground or water surface temperature (K)
    call mpas_pool_get_array(sfc_input, 'snow',      snow)      !dim (nCells Time); snow water equivalent (kg m^-2)
    call mpas_pool_get_array(sfc_input, 'snowc',     snowc)     !dim (nCells Time); flag for snow on ground (=0 no snow; =1,otherwise
    call mpas_pool_get_array(sfc_input, 'snowh',     snowh)     !dim (nCells Time); physical snow depth (m)
    call mpas_pool_get_array(sfc_input, 'sst',       sst)       !dim (nCells Time); sea-surface temperature (K)
    call mpas_pool_get_array(sfc_input, 'ter',       ter)       !dim (nCells); terrain height (m)
    call mpas_pool_get_array(sfc_input, 'tmn' ,      tmn)       !dim (nCells Time); deep soil temperature (K)
    call mpas_pool_get_array(sfc_input, 'tslb',      tslb)      !dim (nSoilLevels nCells Time); soil layer temperature (K)
    call mpas_pool_get_array(sfc_input, 'vegfra',    vegfra)    !dim (nCells Time); vegetation fraction (percent)
    call mpas_pool_get_array(sfc_input, 'seaice',    seaice)    !dim (nCells Time); sea-ice flag (0=no seaice; =1 otherwise)
    call mpas_pool_get_array(sfc_input, 'xice',      xice)      !dim (nCells Time); fractional area coverage of sea-ice
    call mpas_pool_get_array(sfc_input, 'xland',     xland)     !dim (nCells Time); land-ocean mask (1=land including sea-ice ; 2=ocean)
    call mpas_pool_get_array(sfc_input, 'shdmin',    shdmin)    !dim (nCells); minimum fractional coverage of annual green vegetation fraction (percent)
    call mpas_pool_get_array(sfc_input, 'shdmax',    shdmax)    !dim (nCells); maximum fractional coverage of annual green vegetation fraction (percent)
    call mpas_pool_get_array(sfc_input, 'snoalb',    snoalb)    !dim (nCells); annual maximum snow albedo
    call mpas_pool_get_array(sfc_input, 'greenfrac', greenfrac) !dim (nMonths nCells); monthly-mean climatological greenness fraction (percent)
    call mpas_pool_get_array(sfc_input, 'albedo12m', albedo12m) !dim (nMonhts nCells); monthly-mean climatological surface albedo (percent)
    call mpas_pool_get_array(sfc_input, 'canwat',    canwat)    !dim (nCells); water in canopy (kg m^-2)
    call mpas_pool_get_array(sfc_input, 'landusef',  landusef)
    call mpas_pool_get_array(sfc_input, 'soilf',     soilf)
    call mpas_pool_get_array(diag_phys, 'znt',       znt)       !dim (nCells); roughness length (m)
    call mpas_pool_get_array(diag_phys, 'sfc_albedo',sfc_albedo)!dim (nCells); surface albedo (fraction)
    call mpas_pool_get_array(diag_phys, 'sfc_emiss', sfc_emiss)

    do ithread = 1,nThreads
      do iCol = cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread)
        surface % tsfco(iCol)           = sst(iCol)
        surface % tg3(iCol)             = tmn(iCol)
        surface % zorl(iCol)            = znt(iCol)*100.0_RKIND
        surface % alvsf(iCol)           = albbck(iCol)
        surface % alvwf(iCol)           = albbck(iCol)
        surface % alnsf(iCol)           = albbck(iCol)
        surface % alnwf(iCol)           = albbck(iCol)
        surface % facsf(iCol)           = 0.5!? - from gcycle?
        surface % facwf(iCol)           = 0.5!? - from gcycle?
        surface % vfrac(iCol)           = vegfra(iCol)*0.01_RKIND !conversion to decimal from percent
        surface % canopy(iCol)          = canwat(iCol)
        surface % vtype(iCol)           = ivgtyp(iCol)
        surface % stype(iCol)           = isltyp(iCol)
        surface % vegtype_frac(iCol,:)  = landusef(:,iCol)
        surface % soiltype_frac(iCol,:) = soilf(:,iCol)
        surface % fice(iCol)            = xice(iCol) 
        surface % tisfc(iCol)           = skintemp(iCol)
        surface % srflag(iCol)          = 0.0_RKIND
        surface % shdmin(iCol)          = shdmin(iCol) * 0.01_RKIND   ! percent -> fraction
        surface % shdmax(iCol)          = shdmax(iCol) * 0.01_RKIND   ! percent -> fraction
        surface % snoalb(iCol)          = snoalb(iCol)
        surface % tsfc(iCol)            = skintemp(iCol)
        surface % tsfcl(iCol)           = skintemp(iCol)
        surface % zorlw(iCol)           = znt(iCol)*100.0_RKIND
        surface % zorll(iCol)           = znt(iCol)*100.0_RKIND
        surface % zorli(iCol)           = znt(iCol)*100.0_RKIND
        surface % albdirvis_lnd(iCol)   = sfc_albedo(iCol)
        surface % albdirnir_lnd(iCol)   = sfc_albedo(iCol)
        surface % albdifvis_lnd(iCol)   = sfc_albedo(iCol)
        surface % albdifnir_lnd(iCol)   = sfc_albedo(iCol)
        surface % emis_lnd(iCol)        = sfc_emiss(iCol)
        surface % emis_ice(iCol)        = sfc_emiss(iCol)
        surface % weasd(iCol)           = snow(iCol)                  ! SWE   [kg m-2 = mm]
        surface % snowd(iCol)           = snowh(iCol) * 1000._RKIND   ! depth [m -> mm]
        surface % weasdl(iCol)          = snow(iCol)                  ! land SWE   [mm]
        surface % snodl(iCol)           = snowh(iCol) * 1000._RKIND   ! land depth [mm]
        surface % weasdi(iCol)          = snow(iCol)                  ! ice SWE    [mm]
        surface % snodi(iCol)           = snowh(iCol) * 1000._RKIND   ! ice depth  [mm]
        surface % sncovr(iCol)          = snowc(iCol)
        surface % sncovr_ice(iCol)      = snowc(iCol)
        !surface % slope(iCol)           = Not needed for RUC, internal to scheme.

        ! Land/sea masking.
        if (landmask(iCol) == 1) then
           surface % slmsk(iCol) = 1.0_RKIND
           surface % landfrac(iCol) = 1.0_RKIND
        else
           surface % oceanfrac(iCol) = 1.0_RKIND
           if (seaice(iCol) > 0.0_RKIND) then
              surface % slmsk(iCol) = 2.0_RKIND
           end if
        endif

        ! Soild colots (depends on land/sea masking)
        if (nint (surface % slmsk(iCol)) == 1) then
          surface % scolor(iCol) = 4
        else
          surface % scolor(iCol) = 0
        endif

        ! RUC specific fields (allocated wlen lsm=lsm_ruc)
        if (control % lsm == control % lsm_ruc) then
          do iLev=1, nSoilLevels
            surface % sh2o(iCol, iLev)  = sh2o(iLev, iCol)
            surface % smois(iCol, iLev) = smois(iLev, iCol)
            surface % tslb(iCol, iLev)  = tslb(iLev, iCol)
          end do
        end if
        ! smc/stc are allocated for every LSM, but only with Model%lsoil levels
        ! (4 for Noah/Noah-MP) -- not nSoilLevels, which is the MPAS IC's own soil
        ! discretization (9 here, a RUC-style IC) and can exceed Model%lsoil, so
        ! cap the loop or this overruns the smc/stc arrays.
        do iLev=1, min(nSoilLevels, control % lsoil)
          surface % smc(iCol, iLev) = smois(iLev, iCol)
          surface % stc(iCol, iLev) = tslb(iLev, iCol)
       end do
      end do
    end do

    !dzs is defined for every cell in the input data, but the variable in GFS_typedefs only has soil depth as a dimension (is
    !uniform for every grid cell)
    !control % dzs(1) = dzs(1,1)
    !control % zs(1) = -1*control % dzs(1)
    !do iLev = 2, nSoilLevels
    !  control % dzs(iLev) = dzs(iLev,1)
    !  control % zs(iLev) = control % zs(iLev-1) - control % dzs(iLev)
    !end do


  end subroutine ufs_mpas_sfc_to_physics

  !> #########################################################################################
  !> Feed RUC's snow/ice-aware albedo and emissivity back to MPAS sfc_albedo / sfc_emiss,
  !> so the next radiation call sees them (same as stock MPAS, where RUC writes
  !> sfc_albedo/sfc_emiss in driver_lsm and RRTMG reads them on the following call).
  !>
  !> #########################################################################################
  subroutine ufs_mpas_sfc_rad_update(surface, control)
    use GFS_typedefs,       only : GFS_sfcprop_type, GFS_control_type
    use mpas_derived_types, only : mpas_pool_type
    use mpas_pool_routines, only : mpas_pool_get_subpool, mpas_pool_get_array, &
                                   mpas_pool_get_dimension, mpas_pool_get_config
 
    type(GFS_sfcprop_type), intent(in) :: surface
    type(GFS_control_type), intent(in) :: control
 
    type(mpas_pool_type), pointer :: diag_phys
    integer,  pointer :: nThreads, cellSolveThreadStart(:), cellSolveThreadEnd(:)
    logical,  pointer :: config_frac_seaice
    real(RKIND), pointer :: sfc_albedo(:), sfc_emiss(:)
    real(RKIND) :: fi
    integer :: iCol, ithread
    real(RKIND), parameter :: alb_wat = 0.08_RKIND, ems_wat = 0.98_RKIND   ! MPAS open-water values
 
    ! Only meaningful when RUC provides the land/ice albedo and emissivity.
    if (control%lsm /= control%lsm_ruc) return

    ! Get dimension
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions, 'nThreads',             nThreads)
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions, 'cellSolveThreadStart', cellSolveThreadStart)
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions, 'cellSolveThreadEnd',   cellSolveThreadEnd)

    ! MPAS configuration
    call mpas_pool_get_config(domain_ptr % blocklist % configs, 'config_frac_seaice', config_frac_seaice)

    ! Access MPAS data pools.
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'diag_physics', diag_phys)

    ! Grab fields from MPAS pools  
    call mpas_pool_get_array(diag_phys, 'sfc_albedo', sfc_albedo)
    call mpas_pool_get_array(diag_phys, 'sfc_emiss',  sfc_emiss)
    
    do ithread = 1, nThreads
       do iCol = cellSolveThreadStart(ithread), cellSolveThreadEnd(ithread)
          select case (nint(surface%slmsk(iCol)))
          case (1)                                   ! land: RUC albedo/emissivity incl. snow
             sfc_albedo(iCol) = surface % sfalb_lnd(iCol)
             sfc_emiss(iCol)  = surface % emis_lnd(iCol)
          case (2)                                   ! sea ice: RUC ice albedo/emissivity incl. snow on ice
             if (config_frac_seaice) then
                fi = min(1._RKIND, max(0._RKIND, real(surface % fice(iCol), RKIND)))
                sfc_albedo(iCol) = fi*surface % sfalb_ice(iCol) + (1._RKIND-fi)*alb_wat
                sfc_emiss(iCol)  = fi*surface % emis_ice(iCol)  + (1._RKIND-fi)*ems_wat
             else
                sfc_albedo(iCol) = surface % sfalb_ice(iCol)
                sfc_emiss(iCol)  = surface % emis_ice(iCol)
             end if
          case default                               ! open water: keep MPAS values (0.08 / 0.98)
          end select
       end do
    end do
 
  end subroutine ufs_mpas_sfc_rad_update

  !> #########################################################################################
  !> Procedure to populate CCPP data container with MPAS pool data.
  !> Initial condition fields needed by the CCPP GWD parameterization.
  !>
  !> #########################################################################################
  subroutine ufs_mpas_gwd_to_physics(control,surface)
    use GFS_typedefs,         only : GFS_control_type
    use GFS_typedefs,         only : GFS_sfcprop_type
    use mpas_derived_types,   only : mpas_pool_type
    use mpas_pool_routines,   only : mpas_pool_get_subpool, mpas_pool_get_dimension, mpas_pool_get_array
    ! Arguments
    type(GFS_control_type), intent(in   ) :: control
    type(GFS_sfcprop_type), intent(inout) :: surface
    ! Locals
    type(mpas_pool_type), pointer :: mesh_pool, sfc_input
    integer :: iCol, ithread
    integer, pointer :: nThreads, cellSolveThreadStart(:), cellSolveThreadEnd(:)
    real(RKIND), pointer :: var2d(:), con(:), oa1(:), oa2(:), oa3(:), oa4(:)
    real(RKIND), pointer :: ol1(:), ol2(:), ol3(:), ol4(:)
    real(RKIND), pointer :: var2dss(:), conss(:), oa1ss(:), oa2ss(:), oa3ss(:), oa4ss(:)
    real(RKIND), pointer :: ol1ss(:), ol2ss(:), ol3ss(:), ol4ss(:)
    real(RKIND), pointer :: var2dls(:), conls(:), oa1ls(:), oa2ls(:), oa3ls(:), oa4ls(:)
    real(RKIND), pointer :: ol1ls(:), ol2ls(:), ol3ls(:), ol4ls(:)
    character(len=*), parameter :: subname = 'atmos_coupling::ufs_mpas_gwd_to_physics'

    ! Get openMP information
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'nThreads',             nThreads)
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'cellSolveThreadStart', cellSolveThreadStart)
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'cellSolveThreadEnd',   cellSolveThreadEnd)

    ! Access MPAS data pools
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh',      mesh_pool)
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'sfc_input', sfc_input)

    ! DJS: Do we need these? Only used when coupling to YSU scheme?
    call mpas_pool_get_array(sfc_input, 'var2d', var2d)
    call mpas_pool_get_array(sfc_input, 'con',   con)
    call mpas_pool_get_array(sfc_input, 'oa1',   oa1)
    call mpas_pool_get_array(sfc_input,	'oa2',   oa2)
    call mpas_pool_get_array(sfc_input,	'oa3',   oa3)
    call mpas_pool_get_array(sfc_input,	'oa4',   oa4)
    call mpas_pool_get_array(sfc_input, 'ol1',   ol1)
    call mpas_pool_get_array(sfc_input, 'ol2',   ol2)
    call mpas_pool_get_array(sfc_input, 'ol3',   ol3)
    call mpas_pool_get_array(sfc_input, 'ol4',   ol4)
    !
    if (control%gwd_opt==3 .or. control%gwd_opt==33 .or. &
        control%gwd_opt==2 .or. control%gwd_opt==22 ) then
       call mpas_pool_get_array(sfc_input, 'var2dss', var2dss)
       call mpas_pool_get_array(sfc_input, 'conss',   conss)
       call mpas_pool_get_array(sfc_input, 'oa1ss',   oa1ss)
       call mpas_pool_get_array(sfc_input, 'oa2ss',   oa2ss)
       call mpas_pool_get_array(sfc_input, 'oa3ss',   oa3ss)
       call mpas_pool_get_array(sfc_input, 'oa4ss',   oa4ss)
       call mpas_pool_get_array(sfc_input, 'ol1ss',   ol1ss)
       call mpas_pool_get_array(sfc_input, 'ol2ss',   ol2ss)
       call mpas_pool_get_array(sfc_input, 'ol3ss',   ol3ss)
       call mpas_pool_get_array(sfc_input, 'ol4ss',   ol4ss)
       call mpas_pool_get_array(sfc_input, 'var2dls', var2dls)
       call mpas_pool_get_array(sfc_input, 'conls',   conls)
       call mpas_pool_get_array(sfc_input, 'oa1ls',   oa1ls)
       call mpas_pool_get_array(sfc_input, 'oa2ls',   oa2ls)
       call mpas_pool_get_array(sfc_input, 'oa3ls',   oa3ls)
       call mpas_pool_get_array(sfc_input, 'oa4ls',   oa4ls)
       call mpas_pool_get_array(sfc_input, 'ol1ls',   ol1ls)
       call mpas_pool_get_array(sfc_input, 'ol2ls',   ol2ls)
       call mpas_pool_get_array(sfc_input, 'ol3ls',   ol3ls)
       call mpas_pool_get_array(sfc_input, 'ol4ls',   ol4ls)
    end if

    do ithread = 1,nThreads
       do iCol = cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread)
          ! DJS: Do we need these? Only used when coupling to YSU scheme?
          surface % hprime(iCol,1) = var2d(iCol)
          surface % oc(iCol)       = con(iCol)
          surface % oa4(iCol,1)    = oa1(iCol)
          surface % oa4(iCol,2)    = oa2(iCol)
          surface % oa4(iCol,3)    = oa3(iCol)
          surface % oa4(iCol,4)    = oa4(iCol)
          surface % clx(iCol,1)    = ol1(iCol)
          surface % clx(iCol,2)    = ol2(iCol)
          surface % clx(iCol,3)    = ol3(iCol)
          surface % clx(iCol,4)    = ol4(iCol)
          ! DJS: Where are these set? Do we need them?
          !      drag_suite_run doesn't use these...
          !      GFS_GWD_generic_pre sets them to zero, since mntvar(11:14) is initialized to zero.
          !      Not using GFS_GWD_generic_pre in UFS-MPAS, instead these are initializaed to zero elsewhere.

          ! UGWPv1 requires these fields, TODO: do the other options need them?
          if (control%gwd_opt==2) then
             surface % gamma(iCol) = 0.0
             surface % sigma(iCol) = 0.0
             surface % theta(iCol) = 0.0
             surface % elvmax(iCol) = 0.0
          end if
          if (control%gwd_opt==3 .or. control%gwd_opt==33 .or. &
              control%gwd_opt==2 .or. control%gwd_opt==22 ) then
             surface % hprime(iCol,1) = var2dls(iCol)
             surface % oc(iCol)       = conls(iCol)
             surface % oa4(iCol,1)    = oa1ls(iCol)
             surface % oa4(iCol,2)    = oa2ls(iCol)
             surface % oa4(iCol,3)    = oa3ls(iCol)
             surface % oa4(iCol,4)    = oa4ls(iCol)
             surface % clx(iCol,1)    = ol1ls(iCol)
             surface % clx(iCol,2)    = ol2ls(iCol)
             surface % clx(iCol,3)    = ol3ls(iCol)
             surface % clx(iCol,4)    = ol4ls(iCol)
             surface % varss(iCol)    = var2dss(iCol)
             surface % ocss(iCol)     = conss(iCol)
             surface % oa4ss(iCol,1)  = oa1ss(iCol)
             surface % oa4ss(iCol,2)  = oa2ss(iCol)
             surface % oa4ss(iCol,3)  = oa3ss(iCol)
             surface % oa4ss(iCol,4)  = oa4ss(iCol)
             surface % clxss(iCol,1)  = ol1ss(iCol)
             surface % clxss(iCol,2)  = ol2ss(iCol)
             surface % clxss(iCol,3)  = ol3ss(iCol)
             surface % clxss(iCol,4)  = ol4ss(iCol)
          end if
       end do
    end do

  end subroutine ufs_mpas_gwd_to_physics

  !> #########################################################################################
  !> Build the dycore-neutral reference pressure profile that UGWPv1 expects in ak/bk.
  !>
  !> FV3 passes its hybrid sigma-pressure coefficients to control_initialize(); MPAS is on a
  !> terrain-following height coordinate and has no such coefficients, so cires_ugwpv1_init
  !> would otherwise dereference null pointers at cires_ugwpv1_module.F90:252:
  !>
  !>     pmb(k) = ak(k) + pref*bk(k)        ! pref = 1.e5
  !>     zkm(k) = -hpskm*alog(pmb(k)/pref)
  !>
  !> Setting bk = 0 and ak = a reference layer pressure makes that evaluate pmb(k) = ak(k).
  !> The reference used here is the global mean of the MPAS base-state pressure, which is
  !> MPAS's own notion of a reference profile rather than an invented one.
  !>
  !> The mean is taken with a global reduction on purpose. pressure_base is (nVertLevels,
  !> nCells) because the terrain-following coordinate makes the base state terrain-dependent,
  !> so a per-rank mean or a single arbitrary column would change with the MPI decomposition
  !> and break decomposition-independence. Terrain spread in layer height is under 2% above
  !> 15 km, so the choice only perturbs UGWP's launch level by a level or two.
  !>
  !> Ordered surface -> TOA, matching both the MPAS bottom-up convention and the
  !> "ak -pa bk-dimensionless from surf" comment in cires_ugwpv1_module.F90.
  !> Called once, before MPAS_initialize().
  !> #########################################################################################
  subroutine ufs_mpas_reference_pressure(levs, p_ref)
    use mpas_derived_types,   only : mpas_pool_type
    use mpas_derived_types,   only : MPAS_LOG_CRIT
    use mpas_pool_routines,   only : mpas_pool_get_subpool, mpas_pool_get_dimension, mpas_pool_get_array
    use mpas_dmpar,           only : mpas_dmpar_sum_real_array, mpas_dmpar_sum_int
    use mpas_log,             only : mpas_log_write

    ! Arguments
    integer,          intent(in)  :: levs
    real(kind=RKIND), intent(out) :: p_ref(levs+1)   !< reference pressure at layer centres (Pa).
    !< MPAS has no hybrid sigma-pressure coordinate, so unlike FV3 this is not derived
    !< from ak/bk coefficients -- it is computed directly as the global-mean base-state
    !< pressure profile, and handed to the host as the dycore-neutral p_ref quantity
    !< (see GFS_typedefs.F90::control_initialize and cires_ugwpv1_module.F90).

    ! Locals
    type(mpas_pool_type), pointer :: mesh_pool, diag_pool
    integer, pointer :: nCellsSolve, nVertLevels
    real(kind=RKIND), pointer :: pressure_base(:,:)
    real(kind=RKIND), allocatable :: localSum(:), globalSum(:)
    integer :: iCol, iLay, nCellsGlobal
    character(len=*), parameter :: subname = 'atmos_coupling::ufs_mpas_reference_pressure'

    ! Access MPAS data pools
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh', mesh_pool)
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'diag', diag_pool)
    call mpas_pool_get_dimension(mesh_pool, 'nCellsSolve', nCellsSolve)
    call mpas_pool_get_dimension(mesh_pool, 'nVertLevels', nVertLevels)
    call mpas_pool_get_array(diag_pool, 'pressure_base', pressure_base)

    if (nVertLevels /= levs) then
       call mpas_log_write(subname // ' ERROR: nVertLevels does not match levs', &
                           messageType=MPAS_LOG_CRIT)
    end if
    if (.not. associated(pressure_base)) then
       call mpas_log_write(subname // ' ERROR: diag pressure_base is not associated; ' // &
                           'MPAS base state must be initialized before this call', &
                           messageType=MPAS_LOG_CRIT)
    end if

    allocate(localSum(levs), globalSum(levs))
    localSum(:) = 0.0_RKIND

    do iCol = 1, nCellsSolve
       do iLay = 1, levs
          localSum(iLay) = localSum(iLay) + pressure_base(iLay,iCol)
       end do
    end do

    call mpas_dmpar_sum_real_array(domain_ptr % dminfo, levs, localSum, globalSum)
    call mpas_dmpar_sum_int(domain_ptr % dminfo, nCellsSolve, nCellsGlobal)

    if (nCellsGlobal <= 0) then
       call mpas_log_write(subname // ' ERROR: global cell count is not positive', &
                           messageType=MPAS_LOG_CRIT)
    end if

    p_ref(1:levs) = globalSum(1:levs) / real(nCellsGlobal, RKIND)
    ! Model lid. cires_ugwpv1_module.F90:248 documents the top of the profile as zero
    ! pressure; UGWPv1 reads only 1:levs, so this element exists to satisfy the p_ref(levs+1)
    ! declaration in its interface.
    p_ref(levs+1) = 0.0_RKIND

    call mpas_log_write(subname // ': MPAS reference pressure profile from base state, ' // &
                        '$i global cells', intArgs=[nCellsGlobal])
    call mpas_log_write(subname // ': surface p_ref = $r Pa, model-top layer p_ref = $r Pa', &
                        realArgs=[p_ref(1), p_ref(levs)])

    deallocate(localSum, globalSum)
    nullify (mesh_pool)
    nullify (diag_pool)

  end subroutine ufs_mpas_reference_pressure

  !> ########################################################################################
  !> This routine computes physics surface properties using data from the LANDUSE.TBL (read in
  !> during model initialization), and MPAS grid information.
  !> Called once per calendar day.
  !> ########################################################################################
  subroutine ufs_mpas_landuse_update(julday)
    use mpas_derived_types,   only : mpas_pool_type
    use mpas_pool_routines,   only : mpas_pool_get_subpool, mpas_pool_get_dimension
    use mpas_log,             only : mpas_log_write
    use mpas_derived_types,   only : MPAS_LOG_ERR, MPAS_LOG_WARN, MPAS_LOG_CRIT
    use mpas_pool_routines,   only : mpas_pool_get_array, mpas_pool_get_config
    implicit none
    ! Arguments
    integer, intent(in) :: julday
    ! Locals
    type(mpas_pool_type), pointer :: sfc_input, diag_phys, mesh
    real(kind=RKIND), pointer :: latCell(:),snoalb(:),snowc(:),xice(:),albbck(:),embck(:),&
         xicem(:),xland(:),z0(:),mavail(:),sfc_albedo(:),sfc_emiss(:),thc(:),ust(:),znt(:)
    integer, pointer :: nCells,isice,iswater,ivgtyp(:),landmask(:)
    logical,pointer :: config_frac_seaice, config_sfc_albedo
    integer :: iCell,isn,is
    real(kind=RKIND) :: xice_threshold
    character(len=*), parameter :: subname = 'atmos_coupling::ufs_mpas_landuse_update'

    call mpas_log_write(subname //'   Updating MPAS surface properting', messageType=MPAS_LOG_WARN)

    ! Are these really physics flags? (i.e. Do the scheme need theses? Or just here for the surface init?)
    call mpas_pool_get_config( domain_ptr % blocklist % configs, 'config_frac_seaice', config_frac_seaice)
    call mpas_pool_get_config( domain_ptr % blocklist % configs, 'config_sfc_albedo' , config_sfc_albedo)

    ! Access MPAS data pools.
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'diag_physics', diag_phys)
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh',         mesh)
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'sfc_input',    sfc_input)

    ! Get dimensions
    call mpas_pool_get_dimension(mesh,'nCells',nCells)

    ! Access MPAS data pool
    call mpas_pool_get_array(mesh,     'latCell'   ,latCell    )
    call mpas_pool_get_array(sfc_input,'isice_lu'  ,isice      )
    call mpas_pool_get_array(sfc_input,'iswater_lu',iswater    )
    call mpas_pool_get_array(sfc_input,'landmask'  ,landmask   )
    call mpas_pool_get_array(sfc_input,'ivgtyp'    ,ivgtyp     )
    call mpas_pool_get_array(sfc_input,'snoalb'    ,snoalb     )
    call mpas_pool_get_array(sfc_input,'snowc'     ,snowc      )
    call mpas_pool_get_array(sfc_input,'xice'      ,xice       )
    call mpas_pool_get_array(sfc_input,'xland'     ,xland      )
    call mpas_pool_get_array(sfc_input,'sfc_albbck',albbck     )
    call mpas_pool_get_array(diag_phys,'sfc_emibck',embck      )
    call mpas_pool_get_array(diag_phys,'mavail'    ,mavail     )
    call mpas_pool_get_array(diag_phys,'sfc_albedo',sfc_albedo )
    call mpas_pool_get_array(diag_phys,'sfc_emiss' ,sfc_emiss  )
    call mpas_pool_get_array(diag_phys,'thc'       ,thc        )
    call mpas_pool_get_array(diag_phys,'ust'       ,ust        )
    call mpas_pool_get_array(diag_phys,'xicem'     ,xicem      )
    call mpas_pool_get_array(diag_phys,'z0'        ,z0         )
    call mpas_pool_get_array(diag_phys,'znt'       ,znt        )

    if(.not. config_frac_seaice) then
       xice_threshold = 0.5_RKIND
    elseif(config_frac_seaice) then
       xice_threshold = 0.02
    endif

    do iCell = 1, nCells
       !finds the season as function of julian day (summer=1, winter=2): summer in the Northern
       !Hemisphere is defined between March 15th and September 15th (winter otherwise).
       isn = 1
       if(julday.lt.105 .or. julday.ge.288) isn=2
       if(latCell(iCell) .lt. 0.) isn=3-isn

       is = ivgtyp(iCell)

       !set no data points to water:
       if(is.eq.0) is = iswater
       if(.not. config_sfc_albedo) albbck(iCell) = albd(is,isn)/100.
       sfc_albedo(iCell) = albbck(iCell)

       if(snowc(iCell) .gt. 0.5) then
          albbck(iCell) = albd(isice,isn) / 100.
          if(config_sfc_albedo) then
             sfc_albedo(iCell) = snoalb(iCell)
          else
             sfc_albedo(iCell) = albbck(iCell) / (1+scfx(is,isn))
          endif
       endif
       thc(iCell)    = therin(is,isn) / 100.
       z0(iCell)     = sfz0(is,isn) / 100.
       znt(iCell)    = z0(icell)
       embck(iCell)  = sfem(is,isn)
       sfc_emiss(iCell) = embck(iCell)

       if(associated(mavail)) mavail(iCell) = slmo(is,isn)
       if(associated(ust)) ust(iCell) = 0.0001

       !set sea-ice points to land with ice/snow surface properties:
       if(xice(iCell) .ge. xice_threshold) then
          albbck(iCell) = albd(isice,isn) / 100.
          embck(iCell)  = sfem(isice,isn)
          if(config_frac_seaice) then
             !0.08 is the albedo over open water.
             !0.98 is the emissivity over open water.
             sfc_albedo(iCell) = xice(iCell)*albbck(iCell) + (1-xice(iCell))*0.08
             sfc_emiss(iCell)  = xice(iCell)*embck(iCell)  + (1-xice(iCell))*0.98
          else
             sfc_albedo(iCell) = albbck(iCell)
             sfc_emiss(iCell)  = embck(iCell)
          endif
          thc(iCell) = therin(isice,isn) / 100.
          z0(icell)  = sfz0(isice,isn) / 100.
          znt(iCell) = z0(iCell)

          if(associated(mavail)) mavail(iCell) = slmo(isice,isn)
       endif
    end do ! nCells

    call mpas_log_write(subname //'   Finished updating MPAS surface properties ', messageType=MPAS_LOG_WARN)
  end subroutine ufs_mpas_landuse_update

  !> ########################################################################################
  !> Interpolate wind-tendencies from centers to edges of grid-cells.
  !>
  !> See MPAS-Model/src/core_atmosphere/physics/mpas_atmphys_todynamics.F: tend_toEdges
  !> ########################################################################################
  subroutine tend_toEdges(mesh,Ux_tend,Uy_tend,U_tend)
    use mpas_derived_types,   only : mpas_pool_type
    use mpas_pool_routines,   only : mpas_pool_get_dimension, mpas_pool_get_array

    type(mpas_pool_type),intent(in):: mesh
    real(kind=RKIND),intent(in),dimension(:,:),target:: Ux_tend,Uy_tend
    real(kind=RKIND),intent(out),dimension(:,:):: U_tend

    ! locals
    integer:: iEdge
    integer:: cell1, cell2
    integer,pointer:: nCells,nCellsSolve,nEdges
    integer,dimension(:,:),pointer:: cellsOnEdge
    real(kind=RKIND),dimension(:,:),pointer:: east,north,edgeNormalVectors

    ! Grab dimensions
    call mpas_pool_get_dimension(mesh,'nCells',nCells)
    call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
    call mpas_pool_get_dimension(mesh,'nEdges',nEdges)

    ! Grab arrays
    call mpas_pool_get_array(mesh,'east',east)
    call mpas_pool_get_array(mesh,'north',north)
    call mpas_pool_get_array(mesh,'edgeNormalVectors',edgeNormalVectors)
    call mpas_pool_get_array(mesh,'cellsOnEdge',cellsOnEdge)

    do iEdge = 1, nEdges
       cell1 = cellsOnEdge(1,iEdge)
       cell2 = cellsOnEdge(2,iEdge)
       U_tend(:,iEdge) =  Ux_tend(:,cell1) * 0.5 * (edgeNormalVectors(1,iEdge) * east(1,cell1)   &
                                                 +  edgeNormalVectors(2,iEdge) * east(2,cell1)   &
                                                 +  edgeNormalVectors(3,iEdge) * east(3,cell1))  &
                        + Uy_tend(:,cell1) * 0.5 * (edgeNormalVectors(1,iEdge) * north(1,cell1)  &
                                                 +  edgeNormalVectors(2,iEdge) * north(2,cell1)  &
                                                 +  edgeNormalVectors(3,iEdge) * north(3,cell1)) &
                        + Ux_tend(:,cell2) * 0.5 * (edgeNormalVectors(1,iEdge) * east(1,cell2)   &
                                                 +  edgeNormalVectors(2,iEdge) * east(2,cell2)   &
                                                 +  edgeNormalVectors(3,iEdge) * east(3,cell2))  &
                        + Uy_tend(:,cell2) * 0.5 * (edgeNormalVectors(1,iEdge) * north(1,cell2)  &
                                                 +  edgeNormalVectors(2,iEdge) * north(2,cell2)  &
                                                 +  edgeNormalVectors(3,iEdge) * north(3,cell2))
    end do
  end subroutine tend_toEdges

  !> ########################################################################################
  !> Procedure to compute hydrostatic pressure for the physics from the MPAS state.
  !>
  !> ########################################################################################
  subroutine ufs_mpas_hydrostatic_pressure(statein, qv, when)
    use mpas_log,             only : mpas_log_write
    use mpas_derived_types,   only : MPAS_LOG_ERR, MPAS_LOG_WARN, MPAS_LOG_CRIT
    use mpas_derived_types,   only : mpas_pool_type
    use mpas_pool_routines,   only : mpas_pool_get_dimension, mpas_pool_get_array, mpas_pool_get_subpool
    use GFS_typedefs,         only : GFS_statein_type
    use mpas_constants,       only : gravity

    type(GFS_statein_type), intent(inout) :: statein
    real(kind=RKIND), intent(in) :: qv(:,:)
    character(len=*), intent(in) :: when
    !
    type(mpas_pool_type), pointer :: diag_pool, mesh_pool, state_pool
    integer, pointer :: nThreads,  nCellsSolve, nVertLevels
    integer, pointer :: cellSolveThreadStart(:), cellSolveThreadEnd(:)
    integer :: iCol, ithread, iLay
    real(kind=RKIND) :: tem1, fzm_p, fzp_p, z0, z1, z2, w1, w2, rho_moist
    real(kind=RKIND), pointer :: prsi(:,:), prsl(:,:), rho(:,:), exner(:,:), pressure_b(:,:)
    real(kind=RKIND), pointer :: pressure_p(:,:), zgrid(:,:), mass(:,:), zz(:,:)
    character(len=*), parameter :: subname = 'atmos_coupling::ufs_mpas_hydrostatic_pressure'

    call mpas_log_write(subname //'   Computing hydrostatic pressure for '//trim(when), messageType=MPAS_LOG_WARN)

    ! Get openMP information
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'nThreads',             nThreads)
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'cellSolveThreadStart', cellSolveThreadStart)
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'cellSolveThreadEnd',   cellSolveThreadEnd)

    ! Access MPAS data pool
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'diag',  diag_pool)
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh',  mesh_pool)
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'state', state_pool)

    ! Get MPAS dimensions
    call mpas_pool_get_dimension(mesh_pool,  'nCellsSolve', nCellsSolve)
    call mpas_pool_get_dimension(mesh_pool,  'nVertLevels', nVertLevels)

    ! Grab fields from MPAS pools
    call mpas_pool_get_array(diag_pool,  'exner',         exner)
    call mpas_pool_get_array(diag_pool,  'pressure_base', pressure_b)
    call mpas_pool_get_array(diag_pool,  'pressure_p',    pressure_p)
    call mpas_pool_get_array(mesh_pool,  'zgrid',         zgrid)
    call mpas_pool_get_array(state_pool, 'rho_zz',        mass,   1)
    call mpas_pool_get_array(mesh_pool,  'zz',            zz)

    allocate(prsl(nCellsSolve, nVertLevels))
    allocate(prsi(nCellsSolve, nVertLevels + 1))
    allocate(rho( nCellsSolve, nVertLevels))

    ! Pressure (non-hydrostatic)
    do ithread = 1,nThreads
       do iCol = cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread)
          do iLay = 1,nVertLevels
             rho(iCol,iLay) = zz(iLay,iCol) * mass(iLay,iCol)
             prsl(iCol,iLay) = pressure_p(iLay,iCol) + pressure_b(iLay,iCol)
          end do
       end do
    end do

    ! Interpolation of pressure and temperature from layer-center to layer-interface
    do ithread = 1,nThreads
       do iCol = cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread)
          do iLay = 2,nVertLevels
             tem1 = 1./(zgrid(iLay+1,iCol)-zgrid(iLay-1,iCol))
             fzm_p = (zgrid(iLay,  iCol)-zgrid(iLay-1,iCol)) * tem1
             fzp_p = (zgrid(iLay+1,iCol)-zgrid(iLay,  iCol)) * tem1
             statein % tgri(iCol,iLay) = fzm_p*statein % tgrs(iCol,iLay) + fzp_p*statein % tgrs(iCol,iLay-1)
             prsi(iCol,iLay) = fzm_p*prsl(iCol,iLay) + fzp_p*prsl(iCol,iLay-1)
          enddo
       enddo
    enddo

    ! Interpolation of pressure and temperature to the top-of-the-model
    iLay = nVertLevels + 1
    do ithread = 1,nThreads
       do iCol = cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread)
          z0 = zgrid(iLay,iCol)
          z1 = 0.5*(zgrid(iLay  ,iCol)+zgrid(iLay-1,iCol))
          z2 = 0.5*(zgrid(iLay-1,iCol)+zgrid(iLay-2,iCol))
          w1 = (z0-z2)/(z1-z2)
          w2 = 1.-w1
          statein % tgri(iCol,iLay) = w1*statein % tgrs(iCol,iLay-1) + w2*statein % tgrs(iCol,iLay-2)
          prsi(iCol,iLay) = exp(w1*log(prsl(iCol,iLay-1))+w2*log(prsl(iCol,iLay-2)))
       end do
    end do

    ! Recalculate the pressure and temperature  at the surface as an extrapolation of
    ! the pressures in the 2 layers above the surface.
    iLay = 1
    do ithread = 1,nThreads
       do iCol = cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread)
          z0 = zgrid(iLay,iCol)
          z1 = 0.5*(zgrid(iLay  ,iCol)+zgrid(iLay+1,iCol))
          z2 = 0.5*(zgrid(iLay+1,iCol)+zgrid(iLay+2,iCol))
          w1 = (z0-z2)/(z1-z2)
          w2 = 1.-w1
          statein % tgri(iCol,iLay) = w1*statein % tgrs(iCol,iLay) + w2*statein % tgrs(iCol,iLay+1)
          prsi(iCol,iLay) = w1*prsl(iCol,iLay)+w2*prsl(iCol,iLay+1)
       end do
    end do

    ! Calculation of the hydrostatic pressure
    do ithread = 1,nThreads
       do iCol = cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread)
          ! Pressure at layer-interfaces
          iLay = nVertLevels + 1
          statein % prsi(iCol,iLay) = prsi(iCol,iLay)
          do iLay = nVertLevels,1,-1
             rho_moist = rho(iCol,iLay) * (1.+qv(iLay,iCol)) ! Moist-density
             statein % prsi(iCol,iLay)  = statein % prsi(iCol,iLay+1) + &
                  gravity*rho_moist*statein % dzgrid(iCol,iLay)
          end do
          ! Pressure at layer-centers
          do iLay = nVertLevels,1,-1
             statein % prsl(iCol,iLay) = 0.5*(statein % prsi(iCol,iLay+1) + statein % prsi(iCol,iLay) )
          end do
          ! Pressure layer thickness
          do iLay = 1,nVertLevels
             statein % dp(iCol,iLay) = statein % prsi(iCol,iLay) - statein % prsi(iCol,iLay+1)
          end do
          ! Pressure difference across layer-centers
          statein % dp(iCol,1) = statein % prsi(iCol,1) - statein % prsi(iCol,2)
          do iLay = 2,nVertLevels
             statein % dp(iCol,iLay) = statein % prsi(iCol,iLay) - statein % prsi(iCol,iLay+1)
          end do
          do iLay = 1, nVertLevels
             statein % dpc(iCol,iLay) = statein % prsi(iCol,iLay) - statein % prsi(iCol,iLay+1)
          end do
          ! Surface pressure
          statein % pgr(iCol) = statein % prsi(iCol,1)
       end do
    end do

    call mpas_log_write(subname //'   Finished computing hydrostatic pressure for '//trim(when), messageType=MPAS_LOG_WARN)

  end subroutine ufs_mpas_hydrostatic_pressure

  !> ########################################################################################
  !> Procedure to interpolate monthly surface data for current day.
  !>
  !> ########################################################################################
  subroutine ufs_mpas_surface_update(day, month, control, surface)
    use mpas_log,             only : mpas_log_write
    use mpas_derived_types,   only : MPAS_LOG_ERR, MPAS_LOG_WARN, MPAS_LOG_CRIT
    use mpas_derived_types,   only : mpas_pool_type
    use mpas_pool_routines,   only : mpas_pool_get_dimension, mpas_pool_get_array, mpas_pool_get_subpool
    use GFS_typedefs,         only : GFS_sfcprop_type, GFS_control_type

    type(GFS_control_type), intent(in   ) :: control
    type(GFS_sfcprop_type), intent(inout) :: surface
    integer, intent(in) :: day, month
    !
    integer, pointer :: nCellsSolve, landmask(:)
    real(RKIND), pointer :: sfc_albbck(:), albedo12m(:,:), sfc_albedo(:)
    real(RKIND) :: wt1, wt2
    type(mpas_pool_type), pointer :: sfc_input, diag_phys, mesh
    integer, parameter, dimension(12) :: mmd = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    integer ::  iCell, m1, m2
    character(len=*), parameter :: subname = 'atmos_coupling::ufs_mpas_surface_update'

    !
    call mpas_log_write(subname //'   Updating MPAS surface radiative properties ', messageType=MPAS_LOG_WARN)

    ! Determine temporal interpoaltion indices.
    wt1 = 0.5
    wt2 = 0.5
    if (day > 15) then
       m1 = month
       m2 = month + 1
       if (month == 12) m2 = 1
       wt1 = 1. - (day - 15.)/real(mmd(m2),kind=RKIND)
       wt2 = 1. - wt1
    else
       m1 = month - 1
       m2 = month
       if (month == 1) m1 = 12
       wt1 = (15. - day)/real(mmd(m1),kind=RKIND)
       wt2 = 1. - wt1
    endif

    ! Access MPAS data pools.
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'sfc_input', sfc_input)
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh', mesh)
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'diag_physics', diag_phys)

    ! Dimensions
    call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)

    ! Arrays
    call mpas_pool_get_array(sfc_input,'landmask'  , landmask  )
    call mpas_pool_get_array(sfc_input,'albedo12m' , albedo12m )
    call mpas_pool_get_array(sfc_input,'sfc_albbck', sfc_albbck)
    call mpas_pool_get_array(diag_phys,'sfc_albedo', sfc_albedo)

    ! Update background surface albedo
    do iCell = 1, nCellsSolve
       sfc_albbck(iCell) = wt1*albedo12m(m1,iCell) + wt2*albedo12m(m2,iCell)
       sfc_albbck(iCell) = sfc_albbck(iCell) / 100.
       if(landmask(iCell) .eq. 0) sfc_albbck(iCell) = 0.08
    enddo

    surface % alvsf(1:nCellsSolve) = sfc_albbck(1:nCellsSolve)
    surface % alvwf(1:nCellsSolve) = sfc_albbck(1:nCellsSolve)
    surface % alnsf(1:nCellsSolve) = sfc_albbck(1:nCellsSolve)
    surface % alnwf(1:nCellsSolve) = sfc_albbck(1:nCellsSolve)
    ! RUC only turns alvsf/alnsf/... into sfalb_lnd_bck inside lsm_ruc_init;
    ! lsm_ruc_run takes sfalb_lnd_bck as intent(in). Updating alvsf daily has
    ! no effect. Update the RUC background directly (land points only):
    ! RUC specific fields (allocated wlen lsm=lsm_ruc)
    if (control % lsm == control % lsm_ruc) then
       do iCell = 1, nCellsSolve
          if (landmask(iCell) == 1) surface % sfalb_lnd_bck(iCell) = sfc_albbck(iCell)
       enddo
    endif
    !
    call mpas_log_write(subname //'   Finished updating MPAS surface radiative properties ', messageType=MPAS_LOG_WARN)
 end subroutine ufs_mpas_surface_update

 !> ########################################################################################
 !> Initialize the prognostic ozone tracer (o3mr) from the ozone data that the 2015 ozone
 !> photochemistry already uses, when the MPAS initial conditions carry no ozone.
 !>
 !>   It fills (a) the MPAS scalar, so the dycore advects it and ufs_mpas_to_physics picks it
 !>   up on later steps, and (b) Statein%qgrs / Stateout%gq0, because ufs_mpas_to_physics
 !>   has already copied the (zero) tracer for this first step.
 !>   Only columns whose ozone is essentially zero are filled, so ICs that do contain ozone
 !>   (e.g. interpolated from GFS) are left alone.
 !>
 !> ######################################################################################## 
  subroutine ufs_mpas_ozone_init(control, tbd, statein, stateout, mpas_from_ufs_cnst)
    use GFS_typedefs,       only : GFS_control_type, GFS_tbd_type, GFS_statein_type, GFS_stateout_type
    use mpas_derived_types, only : mpas_pool_type, MPAS_LOG_WARN, MPAS_LOG_CRIT
    use mpas_pool_routines, only : mpas_pool_get_subpool, mpas_pool_get_array, mpas_pool_get_dimension
    use mpas_log,           only : mpas_log_write

    type(GFS_control_type),  intent(in)    :: control
    type(GFS_tbd_type),      intent(in)    :: tbd
    type(GFS_statein_type),  intent(inout) :: statein
    type(GFS_stateout_type), intent(inout) :: stateout
    integer, pointer,        intent(in)    :: mpas_from_ufs_cnst(:)

    type(mpas_pool_type), pointer :: state_pool, mesh_pool
    real(kind=RKIND),     pointer :: scalars(:,:,:)
    integer,              pointer :: nCellsSolve, nVertLevels, num_scalars
    integer  :: iCol, iLay, k, kmpas, nfill, nlevoz
    real(RKIND) :: lnp, w, oz
    real(RKIND), parameter :: oz_missing = 1.0e-9_RKIND   ! column max below this = no ozone
    character(len=*), parameter :: subname = 'atmos_coupling::ufs_mpas_ozone_init'

    if (control%ntoz <= 0) return
    if (control%oz_coeff < 6) then
       call mpas_log_write(subname//': ozone forcing data has no climatology coefficient '// &
            '(oz_coeff < 6, 2006 scheme); o3mr not initialized', messageType=MPAS_LOG_WARN)
       return
    end if

    ! Access MPAS data pools
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'state', state_pool)
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh',  mesh_pool)

    ! Get dimensions
    call mpas_pool_get_dimension(state_pool, 'num_scalars', num_scalars)
    call mpas_pool_get_dimension(mesh_pool,  'nCellsSolve', nCellsSolve)
    call mpas_pool_get_dimension(mesh_pool,  'nVertLevels', nVertLevels)

    ! Grab fields from MPAS pools
    call mpas_pool_get_array(state_pool, 'scalars', scalars, timeLevel=1)

    ! MPAS scalar index that holds the UFS ozone tracer
    kmpas = 0
    do k = 1, num_scalars
       if (mpas_from_ufs_cnst(k) == control%ntoz) kmpas = k
    end do
    if (kmpas == 0) call mpas_log_write(subname//': no MPAS scalar maps to ntoz', messageType=MPAS_LOG_CRIT)

    nlevoz = control%levozp
    nfill  = 0
    do iCol = 1, nCellsSolve
       if (maxval(statein%qgrs(iCol,:,control%ntoz)) > oz_missing) cycle   ! IC already has ozone
       nfill = nfill + 1
       do iLay = 1, nVertLevels
          ! same log-pressure interpolation as run_o3prog_2015 (po3 decreases with index)
          lnp = log(statein%prsl(iCol,iLay))
          if (lnp >= control%ozphys%po3(1)) then
             oz = tbd%ozpl(iCol,1,6)
          else if (lnp < control%ozphys%po3(nlevoz)) then
             oz = tbd%ozpl(iCol,nlevoz,6)
          else
             do k = 1, nlevoz-1
                if (lnp < control%ozphys%po3(k) .and. lnp >= control%ozphys%po3(k+1)) then
                   w  = (lnp - control%ozphys%po3(k+1)) / (control%ozphys%po3(k) - control%ozphys%po3(k+1))
                   oz = w*tbd%ozpl(iCol,k,6) + (1.0_RKIND - w)*tbd%ozpl(iCol,k+1,6)
                   exit
                end if
             end do
          end if
          oz = max(oz, 0.0_RKIND)
          scalars(kmpas,iLay,iCol)            = real(oz, RKIND)
          statein%qgrs(iCol,iLay,control%ntoz) = oz
          stateout%gq0(iCol,iLay,control%ntoz) = oz
       end do
    end do

    ! Update halo points.
    call dyn_mpas_exchange_halo('scalars',.true.)
    
    call mpas_log_write(subname//': initialized o3mr from ozone climatology in $i of $i columns', &
         intArgs=(/nfill, nCellsSolve/), messageType=MPAS_LOG_WARN)

  end subroutine ufs_mpas_ozone_init
 
 !> ########################################################################################
 !> Procedure to update surface boundary conditions with input SST and fractional sea-ice
 !> coverage.
 !> NOT YET IMPLEMENTED
 !> ########################################################################################
 subroutine ufs_mpas_sst_update()
   use mpas_derived_types,   only : mpas_pool_type
   use mpas_pool_routines,   only : mpas_pool_get_subpool, mpas_pool_get_dimension
   use mpas_log,             only : mpas_log_write
   use mpas_derived_types,   only : MPAS_LOG_ERR, MPAS_LOG_WARN, MPAS_LOG_CRIT
   use mpas_pool_routines,   only : mpas_pool_get_array, mpas_pool_get_config
   implicit none
   logical,pointer :: config_frac_seaice
   integer :: iCell
   integer,pointer :: nCells
   type(mpas_pool_type), pointer :: sfc_input, diag_phys, mesh
   real(kind=RKIND), pointer :: sfc_albedo(:), sfc_emiss(:), xice(:), xicem(:)
   real(kind=RKIND) :: xice_threshold
   character(len=*), parameter :: subname = 'atmos_coupling::ufs_mpas_sst_update'

   call mpas_log_write(subname //'   Updating MPAS surface(sst) properties', messageType=MPAS_LOG_WARN)

   ! Access MPAS data pools.
   call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'diag_physics', diag_phys)
   call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh',         mesh)
   call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'sfc_input',    sfc_input)

   ! Get dimensions
   call mpas_pool_get_dimension(mesh,'nCells',nCells)
   call mpas_pool_get_config( domain_ptr % blocklist % configs, 'config_frac_seaice', config_frac_seaice)

   ! Grab fields from MPAS pools
   call mpas_pool_get_array(sfc_input,'sfc_albedo', sfc_albedo)
   call mpas_pool_get_array(sfc_input,'sfc_emiss',  sfc_emiss)
   call mpas_pool_get_array(diag_phys,'xicem',      xicem)
   call mpas_pool_get_array(sfc_input,'xice',       xice)

   if(.not. config_frac_seaice) then
      xice_threshold = 0.5_RKIND
   elseif(config_frac_seaice) then
      xice_threshold = 0.02
   endif

   if(config_frac_seaice) then
      do iCell = 1,nCells
         if(xice(iCell) < xice_threshold) xice(iCell) = 0._RKIND
      enddo
    elseif(.not.config_frac_seaice) then
       do iCell = 1,nCells
          if(xice(iCell) >= xice_threshold) then
             xice(iCell) = 1._RKIND
          else
             xice(iCell) = 0._RKIND
          endif
       enddo
    endif

    !update the surface albedo and surface emissivity. before updating xice, sfc_albedo and sfc_emiss
    !are valid according to the earlier value of xice, xicem. now that xice has been updated, we also
    !update sfc_albedo and sfc_emiss accordingly:
    if(config_frac_seaice) then
       do iCell = 1, nCells
          if(xice(iCell) /= xicem(iCell) .and. xicem(iCell) >= xice_threshold) then
             sfc_albedo(iCell) = 0.08 + (sfc_albedo(iCell) -0.08) * xice(iCell)/xicem(iCell)
             sfc_emiss(iCell)  = 0.98 + (sfc_emiss(iCell)-0.98  ) * xice(iCell)/xicem(iCell)
          endif
       enddo
    endif

  end subroutine ufs_mpas_sst_update

  !> ########################################################################################
  !> Mask of the MPAS scalars that are "per kg of air" (mass or number mixing ratios) and
  !> therefore need converting between per kg dry air (MPAS) and per kg moist air (CCPP).
  !> ########################################################################################
  subroutine ufs_mpas_per_kg_air_mask(state_pool, num_scalars, mask)
    use mpas_derived_types, only : mpas_pool_type
    use mpas_pool_routines, only : mpas_pool_get_dimension
    type(mpas_pool_type), pointer     :: state_pool
    integer,              intent(in)  :: num_scalars
    logical,              intent(out) :: mask(num_scalars)
    ! Scalars that are not per kg of air. Add to this list if you register others.
    character(len=16), parameter :: not_per_kg(3) = &
         [character(len=16) :: 'index_sgs_tke', 'index_sigmab', 'index_cld_amt']
    integer, pointer :: idx
    integer :: n
 
    mask(:) = .true.
    do n = 1, size(not_per_kg)
       nullify(idx)
       call mpas_pool_get_dimension(state_pool, trim(not_per_kg(n)), idx)   ! null if not registered
       if (associated(idx)) then
          if (idx >= 1 .and. idx <= num_scalars) mask(idx) = .false.
       end if
    end do
  end subroutine ufs_mpas_per_kg_air_mask
end module atmos_coupling_mod
