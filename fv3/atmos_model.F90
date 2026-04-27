!> @file
!> @brief Driver for the atmospheric model, contains routines to advance the
!>   atmospheric model state by one time step.
!>
!> @details This version of atmos_model_mod has been designed around the implicit
!>    version diffusion scheme of the GCM. It requires two routines to advance
!>    the atmospheric model one time step into the future. These two routines
!>    correspond to the down and up sweeps of the standard tridiagonal solver.
!>    Most atmospheric processes (dynamics,radiation,etc.) are performed
!>    in the down routine. The up routine finishes the vertical diffusion
!>    and computes moisture related terms (convection,large-scale condensation,
!>    and precipitation).
!>    The boundary variables needed by other component models for coupling
!>    are contained in a derived data type. A variable of this derived type
!>    is returned when initializing the atmospheric model. It is used by other
!>    routines in this module and by coupling routines. The contents of
!>    this derived type should only be modified by the atmospheric model.
!>
!> @author

!***********************************************************************
!*                   GNU General Public License                        *
!* This file is a part of fvGFS.                                       *
!*                                                                     *
!* fvGFS is free software; you can redistribute it and/or modify it    *
!* and are expected to follow the terms of the GNU General Public      *
!* License as published by the Free Software Foundation; either        *
!* version 2 of the License, or (at your option) any later version.    *
!*                                                                     *
!* fvGFS is distributed in the hope that it will be useful, but        *
!* WITHOUT ANY WARRANTY; without even the implied warranty of          *
!* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU   *
!* General Public License for more details.                            *
!*                                                                     *
!* For the full text of the GNU General Public License,                *
!* write to: Free Software Foundation, Inc.,                           *
!*           675 Mass Ave, Cambridge, MA 02139, USA.                   *
!* or see:   http://www.gnu.org/licenses/gpl.html                      *
!***********************************************************************
module atmos_model_mod

use mpp_mod,            only: mpp_pe, mpp_root_pe, mpp_clock_id, mpp_clock_begin
use mpp_mod,            only: mpp_clock_end, CLOCK_COMPONENT, MPP_CLOCK_SYNC
use mpp_mod,            only: FATAL, mpp_min, mpp_max, mpp_error, mpp_chksum
use mpp_domains_mod,    only: domain2d
use mpp_mod,            only: mpp_get_current_pelist_name
use mpp_mod,            only: input_nml_file
use fms2_io_mod,        only: file_exists
use fms_mod,            only: write_version_number, stdlog, stdout
use fms_mod,            only: clock_flag_default
use fms_mod,            only: check_nml_error
use diag_manager_mod,   only: diag_send_complete_instant
use time_manager_mod,   only: time_type, get_time, get_date, &
                              operator(+), operator(-), real_to_time_type
use field_manager_mod,  only: MODEL_ATMOS
use tracer_manager_mod, only: get_number_tracers, get_tracer_names, &
                              get_tracer_index, NO_TRACER
use xgrid_mod,          only: grid_box_type
use atmosphere_mod,     only: atmosphere_init
use atmosphere_mod,     only: atmosphere_restart
use atmosphere_mod,     only: atmosphere_end
use atmosphere_mod,     only: atmosphere_state_update
use atmosphere_mod,     only: atmosphere_fill_nest_cpl
use atmosphere_mod,     only: atmos_phys_driver_statein
use atmosphere_mod,     only: atmosphere_control_data
use atmosphere_mod,     only: atmosphere_resolution, atmosphere_domain
use atmosphere_mod,     only: atmosphere_grid_bdry, atmosphere_grid_ctr
use atmosphere_mod,     only: atmosphere_dynamics, atmosphere_diag_axes
use atmosphere_mod,     only: atmosphere_etalvls, atmosphere_hgt
!rab use atmosphere_mod,     only: atmosphere_tracer_postinit
use atmosphere_mod,     only: atmosphere_diss_est, atmosphere_nggps_diag
use atmosphere_mod,     only: atmosphere_scalar_field_halo
use atmosphere_mod,     only: atmosphere_get_bottom_layer
use atmosphere_mod,     only: set_atmosphere_pelist
use atmosphere_mod,     only: Atm, mygrid, get_nth_domain_info
use block_control_mod,  only: block_control_type, define_blocks_packed
use DYCORE_typedefs,    only: DYCORE_data_type, DYCORE_diag_type

use GFS_typedefs,       only: GFS_init_type, GFS_kind_phys => kind_phys
use GFS_restart,        only: GFS_restart_type, GFS_restart_populate
use GFS_diagnostics,    only: GFS_externaldiag_type, &
                              GFS_externaldiag_populate
use CCPP_data,          only: ccpp_suite, GFS_control, &
                              GFS_statein, GFS_stateout, &
                              GFS_grid, GFS_tbd, GFS_cldprop, &
                              GFS_sfcprop, GFS_radtend, &
                              GFS_coupling, GFS_intdiag, &
                              GFS_interstitial
use GFS_init,           only: GFS_initialize
use CCPP_driver,        only: CCPP_step
use mod_ufsatm_util,    only: get_atmos_tracer_types
use stochastic_physics_wrapper_mod, only: stochastic_physics_wrapper,stochastic_physics_wrapper_end

use fv3atm_history_io_mod,    only: fv3atm_diag_register, fv3atm_diag_output,  &
                              DIAG_SIZE
use fv3atm_restart_io_mod,    only: fv3atm_restart_register, &
                                    fv3atm_checksum, &
                                    fv_phy_restart_output, &
                                    fv_sfc_restart_output, &
                                    fv3atm_restart_read, &
                                    fv3atm_restart_write
use fv_ufs_restart_io_mod,    only: fv_dyn_restart_register, &
                                    fv_dyn_restart_output
use fv_iau_mod,         only: iau_external_data_type,getiauforcing,iau_initialize
use module_fv3_config,  only: first_kdt, output_fh,                      &
                              fcst_mpi_comm, fcst_ntasks,                &
                              quilting_restart
use module_block_data,  only: block_atmos_copy, block_data_copy,         &
                              block_data_copy_or_fill,                   &
                              block_data_combine_fractions

#ifdef MOVING_NEST
use fv_moving_nest_main_mod,  only: update_moving_nest, dump_moving_nest
use fv_moving_nest_main_mod,  only: nest_tracker_init
use fv_moving_nest_main_mod,  only: moving_nest_end, nest_tracker_end
use fv_moving_nest_types_mod, only: fv_moving_nest_init
use fv_tracker_mod,           only: check_is_moving_nest, execute_tracker
#endif
!-----------------------------------------------------------------------

implicit none
private

public update_atmos_radiation_physics
public update_atmos_model_state
public update_atmos_model_dynamics
public atmos_model_init, atmos_model_end, atmos_data_type
public atmos_model_exchange_phase_1, atmos_model_exchange_phase_2
public atmos_model_restart
public get_atmos_model_ungridded_dim
public atmos_model_get_nth_domain_info
public addLsmask2grid
public setup_exportdata
public setup_inlinedata
public set_fhzero_loop, InitTimeFromIAUOffset
public get_atmos_tracer_types
public copy2block

interface merge_importfield
  module procedure merge_importfield_with_field
  module procedure merge_importfield_with_scalar
end interface merge_importfield

public merge_importfield
!-----------------------------------------------------------------------

!<PUBLICTYPE >
!> Calculate gradient on cubic sphere grid.
 type atmos_data_type
     integer                       :: axes(4)            !< axis indices (returned by diag_manager) for the atmospheric grid
                                                         !< (they correspond to the x, y, pfull, phalf axes)
     integer, pointer              :: pelist(:) =>null() !< pelist where atmosphere is running.
     integer                       :: layout(2)          !< computer task laytout
     integer                       :: grid_type
     logical                       :: regional           !< true if domain is regional
     logical                       :: nested             !< true if there is a nest
     logical                       :: moving_nest_parent !< true if this grid has a moving nest child
     logical                       :: is_moving_nest     !< true if this is a moving nest grid
     logical                       :: isAtCapTime        !< true if currTime is at the cap driverClock's currTime
     integer                       :: ngrids             !< number of grids
     integer                       :: mygrid             !< current grid
     integer                       :: mlon, mlat         !< longitude and latitude
     integer                       :: iau_offset         !< iau running window length
     logical                       :: pe                 !< current pe.
     real(kind=GFS_kind_phys), pointer, dimension(:)     :: ak, bk
     real(kind=GFS_kind_phys), pointer, dimension(:,:)   :: lon_bnd  => null() !< local longitude axis grid box corners in radians.
     real(kind=GFS_kind_phys), pointer, dimension(:,:)   :: lat_bnd  => null() !< local latitude axis grid box corners in radians.
     real(kind=GFS_kind_phys), pointer, dimension(:,:)   :: lon      => null() !< local longitude axis grid box centers in radians.
     real(kind=GFS_kind_phys), pointer, dimension(:,:)   :: lat      => null() !< local latitude axis grid box centers in radians.
     real(kind=GFS_kind_phys), pointer, dimension(:,:)   :: dx, dy
     real(kind=GFS_kind_phys), pointer, dimension(:,:)   :: area
     real(kind=GFS_kind_phys), pointer, dimension(:,:,:) :: layer_hgt, level_hgt
     type(domain2d)                :: domain             !< domain decomposition
     type(domain2d)                :: domain_for_read    !< domain decomposition
     type(time_type)               :: Time               !< current time
     type(time_type)               :: Time_step          !< atmospheric time step.
     type(time_type)               :: Time_init          !< reference time.
     type(grid_box_type)           :: grid               !< hold grid information needed for 2nd order conservative flux exchange
     type(GFS_externaldiag_type), pointer, dimension(:) :: Diag
 end type atmos_data_type
!</PUBLICTYPE >

!> these two arrays, lon_bnd_work and lat_bnd_work are 'working' arrays, always allocated
!> as (nlon+1, nlat+1) and are used to get the corner lat/lon values from the dycore.
!> these values are then copied to Atmos%lon_bnd, Atmos%lat_bnd which are allocated with
!> sizes that correspond to the corner coordinates distgrid in fcstGrid
real(kind=GFS_kind_phys), pointer, dimension(:,:), save :: lon_bnd_work  => null()
real(kind=GFS_kind_phys), pointer, dimension(:,:), save :: lat_bnd_work  => null()
integer, save :: i_bnd_size, j_bnd_size !< Boundary array size
!> Timing clocks
integer :: fv3Clock, getClock, updClock, setupClock, radClock, physClock

!-----------------------------------------------------------------------
integer :: blocksize    = 1       !< Number of grid points in a block
logical :: chksum_debug = .false. !< Logical for checksum debugging
logical :: dycore_only  = .false. !< Logical for running only dynamical core
logical :: debug        = .false. !< Logical for running debug mode
!logical :: debug        = .true.
logical :: sync         = .false. !< Logical to enable sync for timing
real    :: avg_max_length=3600.   !< Maximum length for time averaging
logical :: ignore_rst_cksum = .false. !< Logical to ignore restart file checksum
namelist /atmos_model_nml/ blocksize, chksum_debug, dycore_only, debug, sync, ccpp_suite, avg_max_length, &
                           ignore_rst_cksum

type (time_type) :: diag_time, diag_time_fhzero !< Time diagnostic and forecast hour zero time diagnostic

!--- concurrent and decoupled radiation and physics variables
!-------------------
!  DYCORE containers
!-------------------
type(DYCORE_data_type),    allocatable :: DYCORE_Data(:)  !< number of blocks

!----------------
!  GFS containers
!----------------
type(GFS_externaldiag_type), target :: GFS_Diag(DIAG_SIZE) !< Contains external diagnostic data
type(GFS_restart_type)     , allocatable, target :: GFS_restart_var(:) !< Contains restart variables

!--------------
! IAU container
!--------------
type(iau_external_data_type)        :: IAU_Data !< number of blocks

!-----------------
!  Block container
!-----------------
type (block_control_type), target   :: Atm_block

!-----------------------------------------------------------------------

character(len=128) :: version = '$Id$'   !< Version control string
character(len=128) :: tagname = '$Name$' !< Version control tag string

#ifdef NAM_phys
  logical,parameter :: flip_vc = .false.
#else
  logical,parameter :: flip_vc = .true.
#endif
  !> Setting constant parameters
  real(kind=GFS_kind_phys), parameter :: zero    = 0.0_GFS_kind_phys,     &
                                         one     = 1.0_GFS_kind_phys,     &
                                         epsln   = 1.0e-10_GFS_kind_phys, &
                                         zorlmin = 1.0e-7_GFS_kind_phys

contains

!#######################################################################
!> @brief Update radiation physics in Atmos
!
!> @details Called every time step as the atmospheric driver to compute the
!>   atmospheric tendencies for dynamics, radiation, vertical diffusion of
!>   momentum, tracers, and heat/moisture.  For heat/moisture only the
!>   downward sweep of the tridiagonal elimination is performed, hence
!>   the name "_down".
!>
!>  @param[in,out] Atmos  Derived-type variable that contains fields needed
!>    by the flux exchange module. These fields describe the atmospheric grid
!>    and are needed to compute/exchange fluxes with other component models.
!>    All fields in this variable type are allocated for the global grid
!>    (without halo regions).
subroutine update_atmos_radiation_physics (Atmos)
!-----------------------------------------------------------------------
  implicit none
  type (atmos_data_type), intent(in) :: Atmos
!--- local variables---
    integer :: idtend, itrac
    integer :: nb, jdat(8), rc, ierr

    if (mpp_pe() == mpp_root_pe() .and. debug) write(6,*) "statein driver"
!--- get atmospheric state from the dynamic core
    call set_atmosphere_pelist()
    call mpp_clock_begin(getClock)
    if (GFS_control%do_skeb) call atmosphere_diss_est (GFS_control%skeb_npass) !  do smoothing for SKEB
    ! SA-3D-TKE added GFS_Tbd (kyf)
    call atmos_phys_driver_statein (GFS_Control, GFS_Statein, GFS_Tbd, Atm_block, flip_vc)
    call mpp_clock_end(getClock)

!--- if dycore only run, set up the dummy physics output state as the input state
    if (dycore_only) then
        GFS_Stateout%gu0 = GFS_Statein%ugrs
        GFS_Stateout%gv0 = GFS_Statein%vgrs
        GFS_Stateout%gt0 = GFS_Statein%tgrs
        GFS_Stateout%gq0 = GFS_Statein%qgrs
    else
      if (mpp_pe() == mpp_root_pe() .and. debug) write(6,*) "setup step"

!--- update GFS_control%jdat(8)
      jdat(:) = 0
      call get_date (Atmos%Time, jdat(1), jdat(2), jdat(3),  &
                                 jdat(5), jdat(6), jdat(7))
      GFS_control%jdat(:) = jdat(:)

!--- execute the atmospheric setup step
      call mpp_clock_begin(setupClock)
      call CCPP_step (step="timestep_init", nblks=Atm_block%nblks, ierr=ierr, dycore='fv3')
      if (ierr/=0)  call mpp_error(FATAL, 'Call to CCPP timestep_init step failed')

      if (GFS_Control%do_sppt .or. GFS_Control%do_shum .or. GFS_Control%do_skeb .or. &
          GFS_Control%lndp_type > 0  .or. GFS_Control%do_ca .or. GFS_Control%do_spp) then
!--- call stochastic physics pattern generation / cellular automata
        call stochastic_physics_wrapper(GFS_Control, GFS_Statein, GFS_Grid, GFS_Sfcprop, GFS_Radtend, GFS_Coupling, Atm_block, ierr)
        if (ierr/=0)  call mpp_error(FATAL, 'Call to stochastic_physics_wrapper failed')
      endif

!--- if coupled, assign coupled fields
      call assign_importdata(Atmos%Time,Atmos%Time_step,Atmos%regional,Atmos%ngrids,rc)
      if (rc/=0)  call mpp_error(FATAL, 'Call to assign_importdata failed')

      ! Currently for FV3ATM, it is only enabled for parent domain coupling
      ! with other model components. In this case, only the parent domain
      ! receives coupled fields through the above assign_importdata step. Thus,
      ! an extra step is needed to fill the coupling variables in the nest,
      ! by downscaling the coupling variables from its parent.
      if (Atmos%isAtCapTime .and. Atmos%ngrids > 1) then
        if (GFS_control%cplocn2atm .or. GFS_control%cplwav2atm) then
          call atmosphere_fill_nest_cpl(Atm_block, GFS_control, GFS_sfcprop)
        endif
      endif

      ! Calculate total non-physics tendencies by substracting old GFS Stateout
      ! variables from new/updated GFS Statein variables (gives the tendencies
      ! due to anything else than physics)
      if (GFS_Control%ldiag3d) then
        idtend = GFS_Control%dtidx(GFS_Control%index_of_x_wind,GFS_Control%index_of_process_non_physics)
        if(idtend>=1) then
          do nb = 1,Atm_block%nblks
            GFS_Intdiag%dtend(GFS_Control%chunk_begin(nb):GFS_Control%chunk_end(nb),:,idtend) = GFS_Intdiag%dtend(GFS_Control%chunk_begin(nb):GFS_Control%chunk_end(nb),:,idtend) &
                 + (GFS_Statein%ugrs(GFS_Control%chunk_begin(nb):GFS_Control%chunk_end(nb),:) - GFS_Stateout%gu0(GFS_Control%chunk_begin(nb):GFS_Control%chunk_end(nb),:))
          enddo
        endif

        idtend = GFS_Control%dtidx(GFS_Control%index_of_y_wind,GFS_Control%index_of_process_non_physics)
        if(idtend>=1) then
          do nb = 1,Atm_block%nblks
            GFS_Intdiag%dtend(GFS_Control%chunk_begin(nb):GFS_Control%chunk_end(nb),:,idtend) = GFS_Intdiag%dtend(GFS_Control%chunk_begin(nb):GFS_Control%chunk_end(nb),:,idtend) &
                 + (GFS_Statein%vgrs(GFS_Control%chunk_begin(nb):GFS_Control%chunk_end(nb),:) - GFS_Stateout%gv0(GFS_Control%chunk_begin(nb):GFS_Control%chunk_end(nb),:))
          enddo
        endif

        idtend = GFS_Control%dtidx(GFS_Control%index_of_temperature,GFS_Control%index_of_process_non_physics)
        if(idtend>=1) then
          do nb = 1,Atm_block%nblks
            GFS_Intdiag%dtend(GFS_Control%chunk_begin(nb):GFS_Control%chunk_end(nb),:,idtend) = GFS_Intdiag%dtend(GFS_Control%chunk_begin(nb):GFS_Control%chunk_end(nb),:,idtend) &
                 + (GFS_Statein%tgrs(GFS_Control%chunk_begin(nb):GFS_Control%chunk_end(nb),:) - GFS_Stateout%gt0(GFS_Control%chunk_begin(nb):GFS_Control%chunk_end(nb),:))
          enddo
        endif

        if (GFS_Control%qdiag3d) then
          do itrac=1,GFS_Control%ntrac
            idtend = GFS_Control%dtidx(itrac+100,GFS_Control%index_of_process_non_physics)
            if(idtend>=1) then
              do nb = 1,Atm_block%nblks
                GFS_Intdiag%dtend(GFS_Control%chunk_begin(nb):GFS_Control%chunk_end(nb),:,idtend) = GFS_Intdiag%dtend(GFS_Control%chunk_begin(nb):GFS_Control%chunk_end(nb),:,idtend) &
                     + (GFS_Statein%qgrs(GFS_Control%chunk_begin(nb):GFS_Control%chunk_end(nb),:,itrac) - GFS_Stateout%gq0(GFS_Control%chunk_begin(nb):GFS_Control%chunk_end(nb),:,itrac))
              enddo
            endif
          enddo
        endif
      endif

      call mpp_clock_end(setupClock)

      if (mpp_pe() == mpp_root_pe() .and. debug) write(6,*) "radiation driver"

!--- execute the atmospheric radiation subcomponent (RRTM)

      call mpp_clock_begin(radClock)
      ! Performance improvement. Only enter if it is time to call the radiation physics.
      if (GFS_control%lsswr .or. GFS_control%lslwr) then
        call CCPP_step (step="radiation", nblks=Atm_block%nblks, ierr=ierr, dycore='fv3')
        if (ierr/=0)  call mpp_error(FATAL, 'Call to CCPP radiation step failed')
      endif
      call mpp_clock_end(radClock)

      if (chksum_debug) then
        if (mpp_pe() == mpp_root_pe()) print *,'RADIATION STEP  ', GFS_control%kdt, GFS_control%fhour
        call fv3atm_checksum(GFS_control, GFS_Statein, GFS_Stateout, GFS_Grid, GFS_Tbd, GFS_Cldprop, GFS_Sfcprop, GFS_Radtend, GFS_Coupling, Atm_block)
      endif

      if (mpp_pe() == mpp_root_pe() .and. debug) write(6,*) "physics driver"

!--- execute the atmospheric physics step1 subcomponent (main physics driver)

      call mpp_clock_begin(physClock)
      call CCPP_step (step="physics", nblks=Atm_block%nblks, ierr=ierr, dycore='fv3')
      if (ierr/=0)  call mpp_error(FATAL, 'Call to CCPP physics step failed')
      call mpp_clock_end(physClock)

      if (chksum_debug) then
        if (mpp_pe() == mpp_root_pe()) print *,'PHYSICS STEP1   ', GFS_control%kdt, GFS_control%fhour
        call fv3atm_checksum(GFS_control, GFS_Statein, GFS_Stateout, GFS_Grid, GFS_Tbd, GFS_Cldprop, GFS_Sfcprop, GFS_Radtend, GFS_Coupling, Atm_block)
      endif

      if (GFS_Control%do_sppt .or. GFS_Control%do_shum .or. GFS_Control%do_skeb .or. &
          GFS_Control%lndp_type > 0  .or. GFS_Control%do_ca ) then

        if (mpp_pe() == mpp_root_pe() .and. debug) write(6,*) "stochastic physics driver"

!--- execute the atmospheric physics step2 subcomponent (stochastic physics driver)

        call mpp_clock_begin(physClock)
        call CCPP_step (step="stochastics", nblks=Atm_block%nblks, ierr=ierr, dycore='fv3')
        if (ierr/=0)  call mpp_error(FATAL, 'Call to CCPP stochastics step failed')
        call mpp_clock_end(physClock)

      endif

      if (chksum_debug) then
        if (mpp_pe() == mpp_root_pe()) print *,'PHYSICS STEP2   ', GFS_control%kdt, GFS_control%fhour
        call fv3atm_checksum(GFS_control, GFS_Statein, GFS_Stateout, GFS_Grid, GFS_Tbd, GFS_Cldprop, GFS_Sfcprop, GFS_Radtend, GFS_Coupling, Atm_block)
      endif
      call getiauforcing(GFS_control,IAU_data,Atm(mygrid))
      if (mpp_pe() == mpp_root_pe() .and. debug) write(6,*) "end of radiation and physics step"

!--- execute the atmospheric timestep finalize step
      call mpp_clock_begin(setupClock)
      call CCPP_step (step="timestep_finalize", nblks=Atm_block%nblks, ierr=ierr, dycore='fv3')
      if (ierr/=0)  call mpp_error(FATAL, 'Call to CCPP timestep_finalize step failed')
      call mpp_clock_end(setupClock)

    endif

    ! Per-timestep diagnostics must be after physics but before
    ! flagging the first timestep.
    if(GFS_control%print_diff_pgr) then
      call atmos_timestep_diagnostics(Atmos)
    endif

    ! Update flag for first time step of time integration
    GFS_control%first_time_step = .false.

!-----------------------------------------------------------------------
 end subroutine update_atmos_radiation_physics
!> @brief Calculate diagnositc information every time step
!> @details Calculates per-timestep, domain-wide, diagnostic, information and
!>   prints to stdout from master rank. Must be called after physics
!>   update but before first_time_step flag is cleared.
!>
!> @param[inout] Atmos Derived-type variable that contains fields needed by the
!>   flux exchange module. These fields describe the atmospheric grid and are needed
!>   to compute/exchange fluxes with other component models.  All fields in this
!>   variable type are allocated for the global grid (without halo regions).
subroutine atmos_timestep_diagnostics(Atmos)
  use mpi_f08
  implicit none
  type (atmos_data_type), intent(in) :: Atmos
!--- local variables---
    integer :: i, nb, count, ierror
    integer :: j
    ! double precision ensures ranks and sums are not truncated
    ! regardless of compilation settings
    double precision :: pdiff, psum, pcount, maxabs, pmaxloc(7), adiff
    double precision :: sendbuf(2), recvbuf(2), global_average

    if(GFS_control%print_diff_pgr) then
      if(.not. GFS_control%first_time_step) then
        pmaxloc = 0.0d0
        recvbuf = 0.0d0
        psum    = 0.0d0
        pcount  = 0.0d0
        maxabs  = 0.0d0

        ! Put pgr stats in pmaxloc, psum, and pcount:
        pmaxloc(1) = GFS_Control%tile_num
        j = 0
        do nb = 1,ATM_block%nblks
          count = size(GFS_Statein%pgr(GFS_Control%chunk_begin(nb):GFS_Control%chunk_end(nb)))
          do i=1,count
            j = j+1
            pdiff = GFS_Statein%pgr(j)-GFS_Intdiag%old_pgr(j)
            adiff = abs(pdiff)
            psum  = psum + adiff
            if(adiff>=maxabs) then
              maxabs=adiff
              pmaxloc(2:3) = (/ dble(ATM_block%index(nb)%ii(i)), dble(ATM_block%index(nb)%jj(i)) /)
              pmaxloc(4:7) = (/ dble(pdiff), dble(GFS_Statein%pgr(j)), &
                   dble(GFS_Grid%xlat(j)), dble(GFS_Grid%xlon(j)) /)
            endif
          enddo
          pcount = pcount+count
        enddo

        ! Sum pgr stats from psum/pcount and convert to hPa/hour global avg:
        sendbuf(1:2) = (/ psum, pcount /)
        call MPI_Allreduce(sendbuf,recvbuf,2,MPI_DOUBLE_PRECISION,MPI_SUM,GFS_Control%communicator,ierror)
        global_average = recvbuf(1)/recvbuf(2) * 36.0d0/GFS_control%dtp

        ! Get the pmaxloc for the global maximum:
        sendbuf(1:2) = (/ maxabs, dble(GFS_Control%me) /)
        call MPI_Allreduce(sendbuf,recvbuf,1,MPI_2DOUBLE_PRECISION,MPI_MAXLOC,GFS_Control%communicator,ierror)
        call MPI_Bcast(pmaxloc,size(pmaxloc),MPI_DOUBLE_PRECISION,nint(recvbuf(2)),GFS_Control%communicator,ierror)

        if(GFS_Control%me == GFS_Control%master) then
2933      format('At forecast hour ',F9.3,' mean abs pgr change is ',F16.8,' hPa/hr')
2934      format('  max abs change   ',F15.10,' bar  at  tile=',I0,' i=',I0,' j=',I0)
2935      format('  pgr at that point',F15.10,' bar      lat=',F12.6,' lon=',F12.6)
          print 2933, GFS_control%fhour, global_average
          print 2934, pmaxloc(4)*1d-5, nint(pmaxloc(1:3))
          print 2935, pmaxloc(5)*1d-5, pmaxloc(6:7)*57.29577951308232d0 ! 180/pi
        endif
      endif
      ! old_pgr is updated every timestep, including the first one where stats aren't printed:
      GFS_Intdiag%old_pgr = GFS_Statein%pgr
    endif

!-----------------------------------------------------------------------
end subroutine atmos_timestep_diagnostics
!> @brief Routine to initialize the atmospheric model
!>
!> @param[inout] Atmos Derived-type variable describing atmospheric grid
!> @param[in] Time_init Reference time
!> @param[in] Time Current time
!> @param[in] Time_step Atmospheric time step
subroutine atmos_model_init (Atmos, Time_init, Time, Time_step)

#ifdef _OPENMP
  use omp_lib
#endif
  use update_ca, only: read_ca_restart

  type (atmos_data_type), intent(inout) :: Atmos
  type (time_type), intent(in) :: Time_init, Time, Time_step
!--- local variables ---
  integer :: unit, i
  ! NEEDED? integer :: j, ix
  integer :: mlon, mlat, nlon, nlat, nlev, sec, sec_lastfhzerofh
  integer :: ierr, io, logunit
  integer :: tile_num
  integer :: isc, iec, jsc, jec
  real(kind=GFS_kind_phys) :: dt_phys
  logical              :: p_hydro, hydro, tmpflag_fhzero
  logical, save        :: block_message = .true.
  type(GFS_init_type)  :: Init_parm
  integer              :: bdat(8), cdat(8)
  integer              :: ntracers
  character(len=32), allocatable, target :: tracer_names(:)
  integer,           allocatable, target :: tracer_types(:)
  integer :: nthrds, nb

!-----------------------------------------------------------------------

!---- set the atmospheric model time ------

   Atmos % isAtCapTime = .false.
   Atmos % Time_init = Time_init
   Atmos % Time      = Time
   Atmos % Time_step = Time_step
   call get_time (Atmos % Time_step, sec)
   dt_phys = real(sec)      ! integer seconds

   logunit = stdlog()

!---------- initialize atmospheric dynamics after reading the namelist -------
!---------- (need name of CCPP suite definition file from input.nml) ---------
   call atmosphere_init (Atmos%Time_init, Atmos%Time, Atmos%Time_step,&
                         Atmos%grid, Atmos%area)
#ifdef MOVING_NEST
   call fv_moving_nest_init(Atm, mygrid)
   call nest_tracker_init()
#endif
!-----------------------------------------------------------------------
   call atmosphere_resolution (nlon, nlat, global=.false.)
   call atmosphere_resolution (mlon, mlat, global=.true.)
   call atmosphere_domain (Atmos%domain, Atmos%domain_for_read, Atmos%layout, &
                           Atmos%regional, Atmos%nested, &
                           Atmos%ngrids, Atmos%mygrid, Atmos%pelist, Atmos%grid_type)
   Atmos%moving_nest_parent = .false.
   Atmos%is_moving_nest = .false.
#ifdef MOVING_NEST
   call check_is_moving_nest(Atm, Atmos%mygrid, Atmos%ngrids, Atmos%is_moving_nest, Atmos%moving_nest_parent)
#endif
   call atmosphere_diag_axes (Atmos%axes)
   call atmosphere_etalvls (Atmos%ak, Atmos%bk, flip=flip_vc)

   tile_num=-1
   call atmosphere_control_data (isc, iec, jsc, jec, nlev, p_hydro, hydro, global_tile_num=tile_num)

   allocate (Atmos%lon(nlon,nlat), Atmos%lat(nlon,nlat))
   call atmosphere_grid_ctr (Atmos%lon, Atmos%lat)

   i_bnd_size = nlon
   j_bnd_size = nlat
   if (iec == mlon) then
      ! we are on task at the 'east' edge of the cubed sphere face or regional domain
      ! corner arrays should have one extra element in 'i' direction
      i_bnd_size = nlon + 1
   end if
   if (jec == mlat) then
      ! we are on task at the 'north' edge of the cubed sphere face or regional domain
      ! corner arrays should have one extra element in 'j' direction
      j_bnd_size = nlat + 1
   end if
   allocate (Atmos%lon_bnd(i_bnd_size,j_bnd_size), Atmos%lat_bnd(i_bnd_size,j_bnd_size))
   allocate (lon_bnd_work(nlon+1,nlat+1), lat_bnd_work(nlon+1,nlat+1))
   call atmosphere_grid_bdry (lon_bnd_work, lat_bnd_work)
   Atmos%lon_bnd(1:i_bnd_size,1:j_bnd_size) = lon_bnd_work(1:i_bnd_size,1:j_bnd_size)
   Atmos%lat_bnd(1:i_bnd_size,1:j_bnd_size) = lat_bnd_work(1:i_bnd_size,1:j_bnd_size)

   call atmosphere_hgt (Atmos%layer_hgt, 'layer', relative=.false., flip=flip_vc)
   call atmosphere_hgt (Atmos%level_hgt, 'level', relative=.false., flip=flip_vc)

   Atmos%mlon = mlon
   Atmos%mlat = mlat

!----------------------------------------------------------------------------------------------
! initialize atmospheric model - must happen AFTER atmosphere_init so that nests work correctly

   if (file_exists('input.nml')) then
      read(input_nml_file, nml=atmos_model_nml, iostat=io)
      ierr = check_nml_error(io, 'atmos_model_nml')
   endif

!-----------------------------------------------------------------------
!--- before going any further check definitions for 'blocks'
!-----------------------------------------------------------------------
   call define_blocks_packed ('atmos_model', Atm_block, isc, iec, jsc, jec, nlev, &
                              blocksize, block_message)

   allocate(DYCORE_Data(Atm_block%nblks))

#ifdef _OPENMP
   nthrds = omp_get_max_threads()
#else
   nthrds = 1
#endif
   allocate(GFS_interstitial(nthrds+1))

!--- update GFS_control%jdat(8)
   bdat(:) = 0
   call get_date (Time_init, bdat(1), bdat(2), bdat(3),  &
                             bdat(5), bdat(6), bdat(7))
   cdat(:) = 0
   call get_date (Time,      cdat(1), cdat(2), cdat(3),  &
                             cdat(5), cdat(6), cdat(7))
   call get_number_tracers(MODEL_ATMOS, num_tracers=ntracers)
   allocate (tracer_names(ntracers), tracer_types(ntracers))
   do i = 1, ntracers
     call get_tracer_names(MODEL_ATMOS, i, tracer_names(i))
   enddo
   call get_atmos_tracer_types(tracer_types)
!--- setup Init_parm
   Init_parm%me              =  mpp_pe()
   Init_parm%master          =  mpp_root_pe()
   Init_parm%fcst_mpi_comm   =  fcst_mpi_comm
   Init_parm%fcst_ntasks     =  fcst_ntasks
   Init_parm%tile_num        =  tile_num
   Init_parm%isc             =  isc
   Init_parm%jsc             =  jsc
   Init_parm%nx              =  nlon
   Init_parm%ny              =  nlat
   Init_parm%levs            =  nlev
   Init_parm%cnx             =  mlon
   Init_parm%cny             =  mlat
   Init_parm%gnx             =  Init_parm%cnx*4
   Init_parm%gny             =  Init_parm%cny*2
   Init_parm%nlunit          =  9999
   Init_parm%logunit         =  logunit
   Init_parm%bdat(:)         =  bdat(:)
   Init_parm%cdat(:)         =  cdat(:)
   Init_parm%dt_dycore       =  dt_phys
   Init_parm%dt_phys         =  dt_phys
   Init_parm%iau_offset      =  Atmos%iau_offset
   Init_parm%blksz           => Atm_block%blksz
   Init_parm%ak              => Atmos%ak
   Init_parm%bk              => Atmos%bk
   Init_parm%xlon            => Atmos%lon
   Init_parm%xlat            => Atmos%lat
   Init_parm%area            => Atmos%area
   Init_parm%nwat            = Atm(mygrid)%flagstruct%nwat
   Init_parm%tracer_names    => tracer_names
   Init_parm%tracer_types    => tracer_types
   Init_parm%restart         = Atm(mygrid)%flagstruct%warm_start
   Init_parm%hydrostatic     = Atm(mygrid)%flagstruct%hydrostatic

   allocate(Init_parm%input_nml_file, mold=input_nml_file)
   Init_parm%input_nml_file = input_nml_file
   Init_parm%fn_nml='using internal file'

   call GFS_initialize (GFS_control, GFS_Statein, GFS_Stateout, GFS_Sfcprop, &
                        GFS_Coupling, GFS_Grid, GFS_Tbd, GFS_Cldprop, GFS_Radtend, &
                        GFS_Intdiag, Init_parm)

   !--- populate/associate the Diag container elements
   call GFS_externaldiag_populate (GFS_Diag, GFS_Control, GFS_Statein, GFS_Stateout,   &
                                             GFS_Sfcprop, GFS_Coupling, GFS_Grid,      &
                                             GFS_Tbd, GFS_Cldprop, GFS_Radtend,        &
                                             GFS_Intdiag, Init_parm)

   Atmos%Diag => GFS_Diag

   Atm(mygrid)%flagstruct%do_skeb = GFS_control%do_skeb

!  initialize the IAU module
   call iau_initialize (GFS_control,IAU_data,Init_parm,Atm(mygrid))

   Init_parm%blksz           => null()
   Init_parm%ak              => null()
   Init_parm%bk              => null()
   Init_parm%xlon            => null()
   Init_parm%xlat            => null()
   Init_parm%area            => null()
   Init_parm%tracer_names    => null()
   deallocate (tracer_names)
   deallocate (tracer_types)

   call atmosphere_nggps_diag (Time, init=.true.)
   call fv3atm_diag_register (GFS_Diag, Time, Atm_block, GFS_control, Atmos%lon, Atmos%lat, Atmos%axes)
   call GFS_restart_populate (GFS_restart_var, GFS_control, GFS_statein, GFS_stateout, GFS_sfcprop, &
                              GFS_coupling, GFS_grid, GFS_tbd, GFS_cldprop,  GFS_Radtend, &
                              GFS_IntDiag, Init_parm, GFS_Diag)
   if (quilting_restart) then
      call fv_dyn_restart_register (Atm(mygrid))
      call fv3atm_restart_register (GFS_Sfcprop, GFS_restart_var, Atm_block, GFS_control)
   endif
   call fv3atm_restart_read (GFS_sfcprop, GFS_restart_var, Atm_block, GFS_control, Atmos%domain_for_read, &
                             Atm(mygrid)%flagstruct%warm_start, ignore_rst_cksum)
  if(GFS_control%do_ca .and. Atm(mygrid)%flagstruct%warm_start)then
    call read_ca_restart (Atmos%domain,3,GFS_control%ncells,GFS_control%nca,GFS_control%ncells_g,GFS_control%nca_g)
  endif
   ! Populate the GFS_Statein container with the prognostic state
   ! in Atm_block, which contains the initial conditions/restart data.
   ! SA-3D-TKE added GFS_Tbd (kyf)
   call atmos_phys_driver_statein (GFS_control, GFS_statein, GFS_Tbd, Atm_block, flip_vc)

   ! When asked to calculate 3-dim. tendencies, set Stateout variables to
   ! Statein variables here in order to capture the first call to dycore
    if (GFS_control%ldiag3d) then
        GFS_Stateout%gu0 = GFS_Statein%ugrs
        GFS_Stateout%gv0 = GFS_Statein%vgrs
        GFS_Stateout%gt0 = GFS_Statein%tgrs
        GFS_Stateout%gq0 = GFS_Statein%qgrs
    endif

   ! Initialize the CCPP framework
   call CCPP_step (step="init", nblks=Atm_block%nblks, ierr=ierr, dycore='fv3')
   if (ierr/=0)  call mpp_error(FATAL, 'Call to CCPP init step failed')
   ! Initialize the CCPP physics
   call CCPP_step (step="physics_init", nblks=Atm_block%nblks, ierr=ierr, dycore='fv3')
   if (ierr/=0)  call mpp_error(FATAL, 'Call to CCPP physics_init step failed')

   if (GFS_Control%do_sppt .or. GFS_Control%do_shum .or. GFS_Control%do_skeb .or. &
       GFS_Control%lndp_type > 0  .or. GFS_Control%do_ca .or. GFS_Control%do_spp) then

!--- Initialize stochastic physics pattern generation / cellular automata for first time step
     call stochastic_physics_wrapper(GFS_control, GFS_statein, GFS_grid, GFS_sfcprop, GFS_Radtend, GFS_Coupling, Atm_block, ierr)
     if (ierr/=0)  call mpp_error(FATAL, 'Call to stochastic_physics_wrapper failed')

   endif

   !--- set the initial diagnostic timestamp
   diag_time = Time
   diag_time_fhzero = Atmos%Time_init
   call get_time (Atmos%Time - Atmos%Time_init, sec)
   call set_fhzero_loop(sec, sec_lastfhzerofh)
   if (mpp_pe() == mpp_root_pe()) print *,'in atmos_model, fhzero=',GFS_Control%fhzero, 'fhour=',sec/3600.,sec_lastfhzerofh/3600

   if (mod((sec-sec_lastfhzerofh),int(GFS_Control%fhzero*3600.)) /= 0) then
     diag_time = Time - real_to_time_type(real(mod(int((GFS_Control%kdt - 1)*dt_phys-sec_lastfhzerofh),int(GFS_Control%fhzero*3600.0))))
     if (mpp_pe() == mpp_root_pe()) print *,'Warning: in atmos_init,start at non multiple of fhzero'
   endif
   if (Atmos%iau_offset > zero) then
     call get_time (Atmos%Time - Atmos%Time_init, sec)
     if (sec < Atmos%iau_offset*3600) then
       diag_time = Atmos%Time_init
       diag_time_fhzero = Atmos%Time
     endif
   endif

   !---- print version number to logfile ----

   call write_version_number ( version, tagname )
   !--- write the namelist to a log file
   if (mpp_pe() == mpp_root_pe()) then
      unit = stdlog( )
      write (unit, nml=atmos_model_nml)
      close (unit)
   endif

   !--- set up clock time

   setupClock = mpp_clock_id( 'GFS Step Setup        ', flags=clock_flag_default, grain=CLOCK_COMPONENT )
   radClock   = mpp_clock_id( 'GFS Radiation         ', flags=clock_flag_default, grain=CLOCK_COMPONENT )
   physClock  = mpp_clock_id( 'GFS Physics           ', flags=clock_flag_default, grain=CLOCK_COMPONENT )
   getClock   = mpp_clock_id( 'Dynamics get state    ', flags=clock_flag_default, grain=CLOCK_COMPONENT )
   updClock   = mpp_clock_id( 'Dynamics update state ', flags=clock_flag_default, grain=CLOCK_COMPONENT )
   if (sync) then
     fv3Clock = mpp_clock_id( 'FV3 Dycore            ', flags=clock_flag_default+MPP_CLOCK_SYNC, grain=CLOCK_COMPONENT )
   else
     fv3Clock = mpp_clock_id( 'FV3 Dycore            ', flags=clock_flag_default, grain=CLOCK_COMPONENT )
   endif

!--- get bottom layer data from dynamical core for coupling
   call atmosphere_get_bottom_layer (Atm_block, DYCORE_Data)

   ! Set flag for first time step of time integration
   GFS_control%first_time_step = .true.

!-----------------------------------------------------------------------
end subroutine atmos_model_init
! </SUBROUTINE>

  !> This will set the forecast hour based on the fhzero_array
  !> input. This should handle both an input of a full fhzero
  !> array and one that uses increments.
  !>
  !> @param tmpflag_fhzero logical if current timestep is in between output hours
  !> @param[inout] sec time since model initialization, in sec
  !> @param[inout] sec_lastfhzerofh time since last fhzero time, in sec
  !>
  !> @author Daniel Sarmiento @date May 16, 2025
subroutine set_fhzero_loop(sec, sec_lastfhzerofh)

   logical                      :: tmpflag_fhzero
   integer                      :: i
   integer, intent(inout)       :: sec, sec_lastfhzerofh

   !--- Model should restart at the forecast hours that are multiples of fhzero.
   !--- WARNING: For special cases that model needs to restart at non-multiple of fhzero
   !--- the fields in first output files are not accumulated from the beginning of
   !--- the bucket, but the restart time.
   if( GFS_Control%fhzero_array(1) > 0. ) then
     fhzero_loop: do i=1,size(GFS_Control%fhzero_array)
       tmpflag_fhzero= .false.
       if( GFS_Control%fhzero_array(i) > 0.) then
         if( i == 1 ) then
           if( sec <= GFS_Control%fhzero_fhour(i)*3600. ) tmpflag_fhzero = .true.
         else if( i > 1 ) then
           if( sec > GFS_Control%fhzero_fhour(i-1)*3600. .and. sec <=GFS_Control%fhzero_fhour(i)*3600. ) &
             tmpflag_fhzero = .true.
         endif
         if( tmpflag_fhzero ) then
           GFS_Control%fhzero = GFS_Control%fhzero_array(i)
           if( GFS_Control%fhzero > 0) then
             sec_lastfhzerofh = (int(sec/3600.)/int(GFS_Control%fhzero))*int(GFS_Control%fhzero)*3600
           else
             sec_lastfhzerofh = 0
           endif
         endif
       endif
     enddo fhzero_loop
   else
     sec_lastfhzerofh = 0
   endif

end subroutine set_fhzero_loop
!> @brief Run the atmospheric dynamics to advect the properties
!>
!> @param[inout] Atmos Derived-type variable describing atmospheric grid
subroutine update_atmos_model_dynamics (Atmos)
! run the atmospheric dynamics to advect the properties
  type (atmos_data_type), intent(in) :: Atmos

    call set_atmosphere_pelist()
#ifdef MOVING_NEST
    ! W. Ramstrom, AOML/HRD -- May 28, 2021
    ! Evaluates whether to move nest, then performs move if needed
    if (Atmos%moving_nest_parent .or. Atmos%is_moving_nest ) then
      call update_moving_nest (Atm_block, GFS_control, GFS_sfcprop, GFS_tbd, &
                               GFS_cldprop, GFS_intdiag, GFS_grid, Atmos%Time)
    endif
#endif
    call mpp_clock_begin(fv3Clock)
    call atmosphere_dynamics (Atmos%Time)
#ifdef MOVING_NEST
    ! W. Ramstrom, AOML/HRD -- June 9, 2021
    ! Debugging output of moving nest code.  Called from this level to access needed input variables.
    if (Atmos%moving_nest_parent .or. Atmos%is_moving_nest ) then
      call dump_moving_nest (Atm_block, GFS_control, GFS_sfcprop, GFS_tbd, Atmos%Time)
    endif
#endif

    call mpp_clock_end(fv3Clock)

end subroutine update_atmos_model_dynamics
!> @brief Perform data exchange with coupled components in run phase 1
!>
!> @details This subroutine currently exports atmospheric fields and tracers
!>   to the chemistry component during the model's run phase 1, i.e.
!>   before chemistry is run.
!>
!> @param[inout] Atmos Derived-type variable describing atmospheric grid
!> @param[out] rc Return code
subroutine atmos_model_exchange_phase_1 (Atmos, rc)

  use ESMF

  type (atmos_data_type), intent(inout) :: Atmos
  integer, optional,      intent(out)   :: rc
!--- local variables
  integer :: localrc

    !--- begin
    if (present(rc)) rc = ESMF_SUCCESS

    !--- if coupled, exchange coupled fields
    if( GFS_control%cplchm ) then
      ! -- export fields to chemistry
      call update_atmos_chemistry('export', rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
    endif

 end subroutine atmos_model_exchange_phase_1
!> @brief Perform data exchange with coupled components in run phase 2
!>
!> @details This subroutine currently imports fields updated by the coupled
!>  chemistry component back into the atmospheric model during run phase 2.
!>
!> @param[inout] Atmos Derived-type variable describing atmospheric grid
!> @param[out] rc Return code
subroutine atmos_model_exchange_phase_2 (Atmos, rc)

  use ESMF

  type (atmos_data_type), intent(inout) :: Atmos
  integer, optional,      intent(out)   :: rc
!--- local variables
  integer :: localrc

    !--- begin
    if (present(rc)) rc = ESMF_SUCCESS

    !--- if coupled, exchange coupled fields
    if( GFS_control%cplchm ) then
      ! -- import fields from chemistry
      call update_atmos_chemistry('import', rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
    endif

end subroutine atmos_model_exchange_phase_2
!> @brief Update the model state after all concurrency is completed
!>
!> @param[inout] Atmos Derived-type variable describing atmospheric grid
!> @param[out] rc Return code
subroutine update_atmos_model_state (Atmos, rc)
! to update the model state after all concurrency is completed
  use ESMF
  type (atmos_data_type), intent(inout) :: Atmos
  integer, optional,      intent(out)   :: rc
!--- local variables
  integer :: i, localrc, sec_lastfhzerofh
  integer :: isec, seconds, isec_fhzero
  integer :: dtatm_temp
  logical :: tmpflag_fhzero
  real(kind=GFS_kind_phys) :: time_int, time_intfull
!
    if (present(rc)) rc = ESMF_SUCCESS

    call set_atmosphere_pelist()
    call mpp_clock_begin(fv3Clock)
    call mpp_clock_begin(updClock)
    call atmosphere_state_update(Atmos%Time, GFS_control, GFS_statein, GFS_stateout, IAU_Data, Atm_block, flip_vc)
#ifdef MOVING_NEST
    call execute_tracker(Atm, mygrid, Atmos%Time, Atmos%Time_step)
#endif
    call mpp_clock_end(updClock)
    call mpp_clock_end(fv3Clock)

    if (chksum_debug) then
      if (mpp_pe() == mpp_root_pe()) print *,'UPDATE STATE    ', GFS_control%kdt, GFS_control%fhour
      call fv3atm_checksum(GFS_control, GFS_statein, GFS_stateout, GFS_grid, GFS_tbd, GFS_cldprop, GFS_sfcprop, GFS_Radtend, GFS_Coupling, Atm_block)
    endif

    !--- advance time ---
    Atmos % Time = Atmos % Time + Atmos % Time_step

    call get_time (Atmos%Time - diag_time, isec)
    call get_time (Atmos%Time - Atmos%Time_init, seconds)
    call get_time (Atmos%Time - diag_time_fhzero, isec_fhzero)
    call atmosphere_nggps_diag(Atmos%Time,ltavg=.true.,avg_max_length=avg_max_length)
    if (ANY(nint(output_fh(:)*3600.0) == seconds) .or. (GFS_control%kdt == first_kdt)) then
      if (mpp_pe() == mpp_root_pe()) write(6,*) "---isec,seconds",isec,seconds
      time_int = real(isec)
      time_intfull = real(seconds)
      call InitTimeFromIAUOffset(Atmos, time_int, time_intfull, seconds, isec_fhzero)
      if (mpp_pe() == mpp_root_pe()) write(6,*) 'gfs diags time since last bucket empty: ',time_int,' time_intfull=', &
         time_intfull,' kdt=',GFS_control%kdt
      call atmosphere_nggps_diag(Atmos%Time)
      call get_time ( Atmos%Time_step, dtatm_temp)
      call fv3atm_diag_output(Atmos%Time, GFS_Diag, Atm_block, GFS_control%nx, GFS_control%ny, &
                            GFS_control%levs, 1, 1, 1.0_GFS_kind_phys, time_int, time_intfull, &
                            GFS_control%fhswr, GFS_control%fhlwr, GFS_control, dtatm_temp)
    endif

    !---  find current fhzero
    call set_fhzero_loop(seconds,sec_lastfhzerofh)
    if (mpp_pe() == mpp_root_pe()) print *,'in atmos_model update, fhzero=',GFS_Control%fhzero, 'fhour=',seconds/3600.,sec_lastfhzerofh/3600.

    if (nint(GFS_Control%fhzero) > 0) then
      if (mod(isec - sec_lastfhzerofh,nint(GFS_Control%fhzero*3600.)) == 0) diag_time = Atmos%Time
!    if (mpp_pe() == mpp_root_pe()) print *,'in atmos_model update time=',isec/3600.,'last fhzeo=',sec_lastfhzerofh
    endif
    call diag_send_complete_instant (Atmos%Time)

    !--- get bottom layer data from dynamical core for coupling
    call atmosphere_get_bottom_layer (Atm_block, DYCORE_Data)

    !--- if in coupled mode, set up coupled fields
    call setup_exportdata(rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__, rcToReturn=rc)) return

    !--- conditionally update the coordinate arrays for moving domains
    if (Atmos%is_moving_nest) then
      call atmosphere_grid_ctr (Atmos%lon, Atmos%lat)
      call atmosphere_grid_bdry (lon_bnd_work, lat_bnd_work, global=.false.)
      Atmos%lon_bnd(1:i_bnd_size,1:j_bnd_size) = lon_bnd_work(1:i_bnd_size,1:j_bnd_size)
      Atmos%lat_bnd(1:i_bnd_size,1:j_bnd_size) = lat_bnd_work(1:i_bnd_size,1:j_bnd_size)
    endif

end subroutine update_atmos_model_state
! </SUBROUTINE>

  !> This will calculate time if an IAU offest has been defined
  !> in the model configuration.
  !>
  !> @param[inout] atmos the main atmos model configurations
  !> @param[inout] time_init model initialization time
  !> @param[inout] time_intfull model time remaining
  !> @param[in]    seconds runtime from model initialization
  !> @param[in]    isec_fhzero model time delta from init to forecast hour 00
  !>
  !> @author Daniel Sarmiento @date May 16, 2025
  !> @date Mar 24, 2026 :: Fix average fields at f000 when using IAU
subroutine InitTimeFromIAUOffset(Atmos, time_int, time_intfull, seconds, isec_fhzero)

  type (atmos_data_type),   intent(inout)  :: Atmos
  real(kind=GFS_kind_phys), intent(inout)  :: time_int, time_intfull
  integer,                  intent(in)     :: seconds, isec_fhzero

  if(Atmos%iau_offset > zero) then
    if( time_int - Atmos%iau_offset*3600. > zero ) then
      time_int = time_int - Atmos%iau_offset*3600.
    else if (seconds == nint(Atmos%iau_offset*3600.)) then
      time_int = real(isec_fhzero)
    endif
    if( time_intfull - Atmos%iau_offset*3600. > zero) then
      time_intfull = time_intfull - Atmos%iau_offset*3600.
    endif
  endif

end subroutine InitTimeFromIAUOffset

!> @brief Terminate routine for atmospheric model
!> @details Call once to terminate this module and any other modules used.
!>   This routine writes a restart file and deallocates storage
!>   used by the derived-type variable atmos_boundary_data_type.
!>
!> @param[inout] Atmos Derived-type variable describing atmospheric grid
subroutine atmos_model_end (Atmos)
  use get_stochy_pattern_mod, only: write_stoch_restart_atm
  use update_ca, only: write_ca_restart
  type (atmos_data_type), intent(inout) :: Atmos
!---local variables
  integer :: ierr

!-----------------------------------------------------------------------
!---- termination routine for atmospheric model ----

#ifdef MOVING_NEST
    !  Call this before atmosphere_end(), because that deallocates Atm
    if (Atmos%is_moving_nest) then
      call moving_nest_end()
      call nest_tracker_end()
    endif
#endif

    call atmosphere_end (Atmos % Time, Atmos%grid, .false.)

    if (GFS_Control%do_sppt .or. GFS_Control%do_shum .or. GFS_Control%do_skeb .or. &
        GFS_Control%lndp_type > 0  .or. GFS_Control%do_ca .or. GFS_Control%do_spp) then
      call stochastic_physics_wrapper_end(GFS_control)
    endif

!   Fast physics (from dynamics) are finalized in atmosphere_end above;
!   standard/slow physics (from CCPP) are finalized in CCPP_step 'physics_finalize'.
    call CCPP_step (step="physics_finalize", nblks=Atm_block%nblks, ierr=ierr, dycore='fv3')
    if (ierr/=0)  call mpp_error(FATAL, 'Call to CCPP physics_finalize step failed')

!   The CCPP framework for all cdata structures is finalized in CCPP_step 'finalize'.
    call CCPP_step (step="finalize", nblks=Atm_block%nblks, ierr=ierr, dycore='fv3')
    if (ierr/=0)  call mpp_error(FATAL, 'Call to CCPP finalize step failed')

    deallocate (Atmos%lon, Atmos%lat)
    deallocate (Atmos%lon_bnd, Atmos%lat_bnd)
    deallocate (lon_bnd_work, lat_bnd_work)

end subroutine atmos_model_end

!> @brief Write out restart files registered through register_restart_file
!>
!> @param[inout] Atmos Derived-type variable describing atmospheric grid
!> @param[in] timestamp Model timestamp
subroutine atmos_model_restart(Atmos, timestamp)
  use update_ca, only: write_ca_restart
  type (atmos_data_type),   intent(inout) :: Atmos
  character(len=*),  intent(in)           :: timestamp

    if (quilting_restart) then
       call fv_sfc_restart_output(GFS_sfcprop, Atm_block, GFS_control)
       call fv_phy_restart_output(GFS_restart_var, Atm_block, GFS_Control)
       call fv_dyn_restart_output(Atm(mygrid), timestamp)
    else
       call atmosphere_restart(timestamp)
       call fv3atm_restart_write (GFS_sfcprop, GFS_restart_var, Atm_block, &
                                  GFS_control, Atmos%domain, timestamp)
    endif
    if(GFS_control%do_ca)then
       call write_ca_restart(timestamp)
    endif
end subroutine atmos_model_restart
!> @brief Retrieve ungridded dimensions of atmospheric model arrays
!>
!> @param[out] nlev Number of atmospheric levels
!> @param[out] nsoillev Number of soil levels
!> @param[out] ntracers Number of atmospheric tracers
subroutine get_atmos_model_ungridded_dim(nlev, nsoillev, ntracers)

  integer, optional, intent(out) :: nlev, nsoillev, ntracers

  !--- number of atmospheric vertical levels
  if (present(nlev)) nlev = Atm_block%npz

  !--- number of soil levels
  if (present(nsoillev)) then
    nsoillev = 0
    if (associated(GFS_Sfcprop%slc)) nsoillev = size(GFS_Sfcprop%slc, dim=2)
  end if

  !--- total number of atmospheric tracers
  if (present(ntracers)) call get_number_tracers(MODEL_ATMOS, num_tracers=ntracers)

end subroutine get_atmos_model_ungridded_dim

!> @brief Populate exported chemistry fields with current atmospheric state
!> @details Update tracer concentrations for atmospheric chemistry with values
!>   from chemistry component (state='import'). Fields should be exported/imported
!>   from/to the atmospheric state after physics calculations.
!>   NOTE: It is assumed that all the chemical tracers follow the standard
!>   atmospheric tracers, which end with ozone. The order of the chemical
!>   tracers must match their order in the chemistry component.
!>
!> @param[in] state Defines wheter field should be imported or exported
!> @param[out] rc Return code
subroutine update_atmos_chemistry(state, rc)

  use ESMF
  use module_cplfields,   only: cplFieldGet

  character(len=*),  intent(in)  :: state
  integer, optional, intent(out) :: rc

  !--- local variables
  integer :: localrc
  integer :: ni, nj, nk, nt, ntb, nte
  integer :: nb, ix, i, j, k, k1, it
  integer :: ib, jb, im

  real(ESMF_KIND_R8), dimension(:,:,:),   pointer :: cldfra,       &
                                                     pfils, pflls, &
                                                     phii,  phil,  &
                                                     prsi,  prsl,  &
                                                     slc,   smc,   &
                                                     stc,   temp,  &
                                                     ua,    va

  real(ESMF_KIND_R8), dimension(:,:,:,:), pointer :: q

!IVAI: add coszens, jo3o1d, jno2, claie, cfch, cfrt, cclu, cpopu
  real(ESMF_KIND_R8), dimension(:,:), pointer :: aod, area, canopy, cmm,  &
    claie, cfch, cfrt, cclu, cpopu, & !IVAI
    dqsfc, dtsfc, fice, flake, focn, fsnow, hpbl, &
    coszens, jo3o1d, jno2, &  !IVAI
    nswsfc, oro, psfc, &
    q2m, rain, rainc, rca, shfsfc, slmsk, stype, swet, t2m, tsfc,    &
    u10m, uustar, v10m, vfrac, xlai, zorl, vtype

! logical, parameter :: diag = .true.

  ! -- begin
  if (present(rc)) rc = ESMF_SUCCESS

  ni  = Atm_block%iec - Atm_block%isc + 1
  nj  = Atm_block%jec - Atm_block%jsc + 1
  nk  = Atm_block%npz

  !--- get total number of tracers
  call get_number_tracers(MODEL_ATMOS, num_tracers=nt)

  select case (trim(state))
    case ('import')
      !--- retrieve references to allocated memory for each field
      call cplFieldGet(state,'inst_tracer_mass_frac', farrayPtr4d=q, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      if (GFS_control%cplaqm) then
        call cplFieldGet(state,'inst_tracer_diag_aod', farrayPtr2d=aod, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__, rcToReturn=rc)) return

!IVAI: case ('import') canopy arrays read in via 'aqm_emis_read'

        if (GFS_control%do_canopy) then
          call cplFieldGet(state,'inst_tracer_diag_claie', farrayPtr2d=claie, rc=localrc)
          if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__, rcToReturn=rc)) return

          call cplFieldGet(state,'inst_tracer_diag_cfch', farrayPtr2d=cfch, rc=localrc)
          if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__, rcToReturn=rc)) return

          call cplFieldGet(state,'inst_tracer_diag_cfrt', farrayPtr2d=cfrt, rc=localrc)
          if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__, rcToReturn=rc)) return

          call cplFieldGet(state,'inst_tracer_diag_cclu', farrayPtr2d=cclu, rc=localrc)
          if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__, rcToReturn=rc)) return

          call cplFieldGet(state,'inst_tracer_diag_cpopu', farrayPtr2d=cpopu, rc=localrc)
          if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__, rcToReturn=rc)) return
        end if

!IVAI: case ('import') photdiag arrays
        call cplFieldGet(state,'inst_tracer_diag_coszens', farrayPtr2d=coszens, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__, rcToReturn=rc)) return

        call cplFieldGet(state,'inst_tracer_diag_jo3o1d', farrayPtr2d=jo3o1d, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__, rcToReturn=rc)) return

        call cplFieldGet(state,'inst_tracer_diag_jno2', farrayPtr2d=jno2, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__, rcToReturn=rc)) return
!IVAI
      end if

      !--- do not import tracer concentrations by default
      ntb = nt + 1
      nte = nt

      !--- if chemical tracers are present, set bounds appropriately
      if (GFS_control%ntchm > 0) then
        ntb = GFS_control%ntchs
        nte = GFS_control%ntche
      end if

      !--- prognostic tracer concentrations
      do it = ntb, nte
!$OMP parallel do default (none) &
!$OMP             shared  (it, nk, nj, ni, Atm_block, GFS_Control, GFS_Stateout, q)  &
!$OMP             private (k, j, jb, i, ib, nb, ix, im)
        do k = 1, nk
          do j = 1, nj
            jb = j + Atm_block%jsc - 1
            do i = 1, ni
              ib = i + Atm_block%isc - 1
              nb = Atm_block%blkno(ib,jb)
              ix = Atm_block%ixp(ib,jb)
              im = GFS_Control%chunk_begin(nb)+ix-1
              GFS_Stateout%gq0(im,k,it) = q(i,j,k,it)
            enddo
          enddo
        enddo
      enddo

      !--- diagnostic tracers
      !--- set tracer concentrations in the atmospheric state directly
      !--- since the atmosphere's driver cannot perform this step while
      !--- updating the state
      if (GFS_control%ndchm > 0) then
        ntb = GFS_control%ndchs
        nte = GFS_control%ndche
!$OMP parallel do default (none) &
!$OMP             shared  (mygrid, nk, ntb, nte, Atm, Atm_block, q) &
!$OMP             private (i, ib, ix, j, jb, k, k1, nb)
        do nb = 1, Atm_block%nblks
          do k = 1, nk
            if(flip_vc) then
              k1 = nk+1-k !reverse the k direction
            else
              k1 = k
            endif
            do ix = 1, Atm_block%blksz(nb)
              ib = Atm_block%index(nb)%ii(ix)
              jb = Atm_block%index(nb)%jj(ix)
              i = ib - Atm_block%isc + 1
              j = jb - Atm_block%jsc + 1
              Atm(mygrid)%q(ib,jb,k1,ntb:nte) = q(i,j,k,ntb:nte)
            enddo
          end do
        end do
      end if

      if (GFS_control%cplaqm) then
        !--- other diagnostics
!$OMP   parallel do default (none) &
!$OMP               shared  (nj, ni, Atm_block, GFS_Control, GFS_Intdiag, aod) &
!$OMP               private (j, jb, i, ib, nb, ix, im)
        do j = 1, nj
          jb = j + Atm_block%jsc - 1
          do i = 1, ni
            ib = i + Atm_block%isc - 1
            nb = Atm_block%blkno(ib,jb)
            ix = Atm_block%ixp(ib,jb)
            im = GFS_Control%chunk_begin(nb)+ix-1
            GFS_IntDiag%aod(im) = aod(i,j)
          enddo
        enddo

        if (GFS_control%do_canopy) then
!IVAI: case ('import') canopy arrays read in via aqm_emis_read
!$OMP   parallel do default (none) &
!$OMP               shared  (nj, ni, Atm_block, GFS_Control, GFS_Intdiag, claie) &
!$OMP               private (j, jb, i, ib, nb, ix, im)
          do j = 1, nj
            jb = j + Atm_block%jsc - 1
            do i = 1, ni
              ib = i + Atm_block%isc - 1
              nb = Atm_block%blkno(ib,jb)
              ix = Atm_block%ixp(ib,jb)
              im = GFS_Control%chunk_begin(nb)+ix-1
              GFS_IntDiag%claie(im) = claie(i,j)
            enddo
          enddo

!$OMP   parallel do default (none) &
!$OMP               shared  (nj, ni, Atm_block, GFS_Control, GFS_Intdiag, cfch) &
!$OMP               private (j, jb, i, ib, nb, ix, im)
          do j = 1, nj
            jb = j + Atm_block%jsc - 1
            do i = 1, ni
              ib = i + Atm_block%isc - 1
              nb = Atm_block%blkno(ib,jb)
              ix = Atm_block%ixp(ib,jb)
              im = GFS_Control%chunk_begin(nb)+ix-1
              GFS_IntDiag%cfch(im) = cfch(i,j)
            enddo
          enddo

!$OMP   parallel do default (none) &
!$OMP               shared  (nj, ni, Atm_block, GFS_Control, GFS_Intdiag, cfrt) &
!$OMP               private (j, jb, i, ib, nb, ix, im)
          do j = 1, nj
            jb = j + Atm_block%jsc - 1
            do i = 1, ni
              ib = i + Atm_block%isc - 1
              nb = Atm_block%blkno(ib,jb)
              ix = Atm_block%ixp(ib,jb)
              im = GFS_Control%chunk_begin(nb)+ix-1
              GFS_IntDiag%cfrt(im) = cfrt(i,j)
            enddo
          enddo

!$OMP   parallel do default (none) &
!$OMP               shared  (nj, ni, Atm_block, GFS_Control, GFS_Intdiag, cclu) &
!$OMP               private (j, jb, i, ib, nb, ix, im)
          do j = 1, nj
            jb = j + Atm_block%jsc - 1
            do i = 1, ni
              ib = i + Atm_block%isc - 1
              nb = Atm_block%blkno(ib,jb)
              ix = Atm_block%ixp(ib,jb)
              im = GFS_Control%chunk_begin(nb)+ix-1
              GFS_IntDiag%cclu(im) = cclu(i,j)
            enddo
          enddo

!$OMP   parallel do default (none) &
!$OMP               shared  (nj, ni, Atm_block, GFS_Control, GFS_Intdiag, cpopu) &
!$OMP               private (j, jb, i, ib, nb, ix, im)
          do j = 1, nj
            jb = j + Atm_block%jsc - 1
            do i = 1, ni
              ib = i + Atm_block%isc - 1
              nb = Atm_block%blkno(ib,jb)
              ix = Atm_block%ixp(ib,jb)
              im = GFS_Control%chunk_begin(nb)+ix-1
              GFS_IntDiag%cpopu(im) = cpopu(i,j)
            enddo
          enddo
        endif ! GFS_control%do_canopy

!IVAI: case ('import') photdiag arrays
!$OMP   parallel do default (none) &
!$OMP               shared  (nj, ni, Atm_block, GFS_Control, GFS_Intdiag, coszens) &
!$OMP               private (j, jb, i, ib, nb, ix, im)
        do j = 1, nj
          jb = j + Atm_block%jsc - 1
          do i = 1, ni
            ib = i + Atm_block%isc - 1
            nb = Atm_block%blkno(ib,jb)
            ix = Atm_block%ixp(ib,jb)
            im = GFS_Control%chunk_begin(nb)+ix-1
            GFS_IntDiag%coszens(im) = coszens(i,j)
          enddo
        enddo

!$OMP   parallel do default (none) &
!$OMP               shared  (nj, ni, Atm_block, GFS_Control, GFS_Intdiag, jo3o1d) &
!$OMP               private (j, jb, i, ib, nb, ix, im)
        do j = 1, nj
          jb = j + Atm_block%jsc - 1
          do i = 1, ni
            ib = i + Atm_block%isc - 1
            nb = Atm_block%blkno(ib,jb)
            ix = Atm_block%ixp(ib,jb)
            im = GFS_Control%chunk_begin(nb)+ix-1
            GFS_IntDiag%jo3o1d(im) = jo3o1d(i,j)
          enddo
        enddo

!$OMP   parallel do default (none) &
!$OMP               shared  (nj, ni, Atm_block, GFS_Control, GFS_Intdiag, jno2) &
!$OMP               private (j, jb, i, ib, nb, ix, im)
        do j = 1, nj
          jb = j + Atm_block%jsc - 1
          do i = 1, ni
            ib = i + Atm_block%isc - 1
            nb = Atm_block%blkno(ib,jb)
            ix = Atm_block%ixp(ib,jb)
            im = GFS_Control%chunk_begin(nb)+ix-1
            GFS_IntDiag%jno2(im) = jno2(i,j)
          enddo
        enddo
!IVAI
      end if

      if (GFS_control%debug) then
        write(6,'("update_atmos: ",a,": qgrs - min/max/avg",3g16.6)') &
          trim(state), minval(q), maxval(q), sum(q)/size(q)
        if (GFS_control%cplaqm) &
          write(6,'("update_atmos: ",a,": aod  - min/max    ",3g16.6)') &
            trim(state), minval(aod), maxval(aod)
!IVAI: case ('import') canopy arrays read via aqm_emis_read
        if (GFS_control%cplaqm .and. GFS_control%do_canopy) &
          write(6,'("update_atmos: ",a,": claie - min/max    ",3g16.6)') &
            trim(state), minval(claie), maxval(claie)
        if (GFS_control%cplaqm .and. GFS_control%do_canopy) &
          write(6,'("update_atmos: ",a,": cfch  - min/max    ",3g16.6)') &
            trim(state), minval(cfch), maxval(cfch)
        if (GFS_control%cplaqm .and. GFS_control%do_canopy) &
          write(6,'("update_atmos: ",a,": cfrt  - min/max    ",3g16.6)') &
            trim(state), minval(cfrt), maxval(cfrt)
        if (GFS_control%cplaqm .and. GFS_control%do_canopy) &
          write(6,'("update_atmos: ",a,": cclu  - min/max    ",3g16.6)') &
            trim(state), minval(cclu), maxval(cclu)
        if (GFS_control%cplaqm .and. GFS_control%do_canopy) &
          write(6,'("update_atmos: ",a,": cpopu - min/max    ",3g16.6)') &
            trim(state), minval(cpopu), maxval(cpopu)
!IVAI: case ('import') photdiag arrays
        if (GFS_control%cplaqm) &
          write(6,'("update_atmos: ",a,": coszens - min/max    ",3g16.6)') &
            trim(state), minval(coszens), maxval(coszens)
        if (GFS_control%cplaqm) &
          write(6,'("update_atmos: ",a,": jo3o1d  - min/max    ",3g16.6)') &
            trim(state), minval(jo3o1d), maxval(jo3o1d)
        if (GFS_control%cplaqm) &
          write(6,'("update_atmos: ",a,": jno2    - min/max    ",3g16.6)') &
            trim(state), minval(jno2), maxval(jno2)
!IVAI
      end if

    case ('export')
      !--- retrieve references to allocated memory for each field
      call cplFieldGet(state,'inst_pres_levels', farrayPtr3d=prsl, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'inst_geop_levels', farrayPtr3d=phil, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'inst_geop_interface', farrayPtr3d=phii, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'inst_temp_levels', farrayPtr3d=temp, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'inst_zonal_wind_levels', farrayPtr3d=ua, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'inst_merid_wind_levels', farrayPtr3d=va, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'inst_tracer_mass_frac', farrayPtr4d=q, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'inst_pbl_height', farrayPtr2d=hpbl, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'surface_cell_area', farrayPtr2d=area, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'inst_convective_rainfall_amount', &
                       farrayPtr2d=rainc, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'inst_friction_velocity', farrayPtr2d=uustar, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'inst_rainfall_amount', farrayPtr2d=rain, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'inst_land_sea_mask', farrayPtr2d=slmsk, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'inst_temp_height_surface', farrayPtr2d=tsfc, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'inst_surface_roughness', farrayPtr2d=zorl, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'inst_cloud_frac_levels', farrayPtr3d=cldfra, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'inst_zonal_wind_height10m', farrayPtr2d=u10m, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'inst_merid_wind_height10m', farrayPtr2d=v10m, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'ice_fraction_in_atm', farrayPtr2d=fice, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'surface_snow_area_fraction', farrayPtr2d=fsnow, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      call cplFieldGet(state,'inst_pres_interface', farrayPtr3d=prsi, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      if (GFS_Control%cplaqm) then

        call cplFieldGet(state,'canopy_moisture_storage', farrayPtr2d=canopy, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__, rcToReturn=rc)) return

        call cplFieldGet(state,'inst_aerodynamic_conductance', farrayPtr2d=cmm, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__, rcToReturn=rc)) return

        call cplFieldGet(state,'inst_laten_heat_flx', farrayPtr2d=dqsfc, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__, rcToReturn=rc)) return

        call cplFieldGet(state,'inst_sensi_heat_flx', farrayPtr2d=dtsfc, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__, rcToReturn=rc)) return

        call cplFieldGet(state,'inst_net_sw_flx', farrayPtr2d=nswsfc, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__, rcToReturn=rc)) return

        call cplFieldGet(state,'height', farrayPtr2d=oro, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__, rcToReturn=rc)) return

        call cplFieldGet(state,'inst_pres_height_surface', farrayPtr2d=psfc, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__, rcToReturn=rc)) return

        call cplFieldGet(state,'inst_spec_humid_height2m', farrayPtr2d=q2m, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__, rcToReturn=rc)) return

        call cplFieldGet(state,'inst_canopy_resistance', farrayPtr2d=rca, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__, rcToReturn=rc)) return

        call cplFieldGet(state,'inst_soil_moisture_content', farrayPtr3d=smc, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__, rcToReturn=rc)) return

        call cplFieldGet(state,'temperature_of_soil_layer', farrayPtr3d=stc, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__, rcToReturn=rc)) return

        call cplFieldGet(state,'inst_temp_height2m', farrayPtr2d=t2m, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__, rcToReturn=rc)) return

        call cplFieldGet(state,'inst_vegetation_area_frac', farrayPtr2d=vfrac, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__, rcToReturn=rc)) return

        call cplFieldGet(state,'leaf_area_index', farrayPtr2d=xlai, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__, rcToReturn=rc)) return

        call cplFieldGet(state,'soil_type', farrayPtr2d=stype, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__, rcToReturn=rc)) return

        call cplFieldGet(state,'vegetation_type', farrayPtr2d=vtype, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      else

        call cplFieldGet(state,'inst_liq_nonconv_tendency_levels', &
                         farrayPtr3d=pflls, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__, rcToReturn=rc)) return

        call cplFieldGet(state,'inst_ice_nonconv_tendency_levels', &
                         farrayPtr3d=pfils, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__, rcToReturn=rc)) return

        call cplFieldGet(state,'lake_fraction', farrayPtr2d=flake, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__, rcToReturn=rc)) return

        call cplFieldGet(state,'ocean_fraction', farrayPtr2d=focn, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__, rcToReturn=rc)) return

        call cplFieldGet(state,'inst_up_sensi_heat_flx', farrayPtr2d=shfsfc, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__, rcToReturn=rc)) return

        call cplFieldGet(state,'inst_soil_moisture_content', farrayPtr3d=slc, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__, rcToReturn=rc)) return

        call cplFieldGet(state,'inst_surface_soil_wetness', farrayPtr2d=swet, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      end if

      !--- handle all three-dimensional variables

      !--- interface values
      phii = reshape(GFS_Statein%phii, shape(phii))
      prsi = reshape(GFS_Statein%prsi, shape(prsi))
      !--- layer values
      prsl = reshape(GFS_Statein%prsl, shape(prsl))
      phil = reshape(GFS_Statein%phil, shape(phil))
      temp = reshape(GFS_Stateout%gt0, shape(temp))
      ua   = reshape(GFS_Stateout%gu0, shape(ua  ))
      va   = reshape(GFS_Stateout%gv0, shape(va  ))
      cldfra = reshape(GFS_IntDiag%cldfra, shape(cldfra))

      if (.not.GFS_Control%cplaqm) then
        !--- layer values
        pfils = reshape(GFS_Coupling%pfi_lsan, shape(pfils))
        pflls = reshape(GFS_Coupling%pfl_lsan, shape(pflls))
      end if

      !--- top interface values
      k = nk+1
!$omp parallel do default(shared) private(i,j,jb,ib,nb,ix)
      do j = 1, nj
        jb = j + Atm_block%jsc - 1
        do i = 1, ni
          ib = i + Atm_block%isc - 1
          nb = Atm_block%blkno(ib,jb)
          ix = Atm_block%ixp(ib,jb)
          im = GFS_Control%chunk_begin(nb)+ix-1
          phii(i,j,k) = GFS_Statein%phii(im,k)
          prsi(i,j,k) = GFS_Statein%prsi(im,k)
        enddo
      enddo

      !--- tracers quantities
      do it = 1, nt
!$OMP parallel do default (none) &
!$OMP             shared  (it, nk, nj, ni, Atm_block, GFS_Control, GFS_Stateout, q)  &
!$OMP             private (k, j, jb, i, ib, nb, ix, im)
        do k = 1, nk
          do j = 1, nj
            jb = j + Atm_block%jsc - 1
            do i = 1, ni
              ib = i + Atm_block%isc - 1
              nb = Atm_block%blkno(ib,jb)
              ix = Atm_block%ixp(ib,jb)
              im = GFS_Control%chunk_begin(nb)+ix-1
              q(i,j,k,it) = GFS_Stateout%gq0(im,k,it)
            enddo
          enddo
        enddo
      enddo

      area = reshape(GFS_Grid%area, shape(area))
      hpbl = reshape(GFS_Tbd%hpbl, shape(hpbl))
      uustar = reshape(GFS_Sfcprop%uustar, shape(uustar))
      slmsk = reshape(GFS_Sfcprop%slmsk, shape(slmsk))
      zorl = reshape(GFS_Sfcprop%zorl, shape(zorl))
      fice = reshape(GFS_Sfcprop%fice, shape(fice))
      fsnow = reshape(GFS_Sfcprop%sncovr, shape(fsnow))

      rainc = reshape(GFS_Coupling%rainc_cpl, shape(rainc))
      rain = reshape(GFS_Coupling%rain_cpl, shape(rain))  &
           + reshape(GFS_Coupling%snow_cpl, shape(rain))
      tsfc = reshape(GFS_Coupling%tsfci_cpl, shape(tsfc))
      u10m = reshape(GFS_Coupling%u10mi_cpl, shape(u10m))
      v10m = reshape(GFS_Coupling%v10mi_cpl, shape(v10m))

      if (GFS_Control%cplaqm) then
        cmm = reshape(GFS_IntDiag%cmm, shape(cmm))
        canopy = reshape(GFS_Sfcprop%canopy, shape(canopy))
        !oro(i,j)    = max(0.d0, GFS_Data(nb)%Sfcprop%oro(ix))
        oro = reshape(GFS_Sfcprop%oro, shape(oro))
        where (oro < 0.d0) oro = 0.d0
        rca = reshape(GFS_Sfcprop%rca, shape(rca))
        !smc(i,j,:)  = GFS_Data(nb)%Sfcprop%smc(ix,:)
        !stc(i,j,:)  = GFS_Data(nb)%Sfcprop%stc(ix,:)
        smc = reshape(GFS_Sfcprop%smc, shape(smc))
        stc = reshape(GFS_Sfcprop%stc, shape(stc))
        vfrac = reshape(GFS_Sfcprop%vfrac, shape(vfrac))
        xlai = reshape(GFS_Sfcprop%xlaixy, shape(xlai))
        !if (nint(slmsk(i,j)) == 2) then
        !  if (GFS_Control%isot == 1) then
        !    stype(i,j) = 16._ESMF_KIND_R8
        !  else
        !    stype(i,j) = 9._ESMF_KIND_R8
        !  endif
        !else
        !  stype(i,j) = real(int( GFS_Data(nb)%Sfcprop%stype(ix)+0.5 ), kind=ESMF_KIND_R8)
        !endif
        stype = real(int(reshape(GFS_Sfcprop%stype, shape(stype))+0.5), kind=ESMF_KIND_R8)
        vtype = real(int(reshape(GFS_Sfcprop%vtype, shape(vtype))+0.5), kind=ESMF_KIND_R8)
        if (GFS_Control%isot == 1) then
          where (slmsk == 2) stype = 16._ESMF_KIND_R8
        else
          where (slmsk == 2) stype = 9._ESMF_KIND_R8
        end if
        dqsfc  = reshape(GFS_Coupling%dqsfci_cpl, shape(dqsfc))
        dtsfc  = reshape(GFS_Coupling%dtsfci_cpl, shape(dtsfc))
        nswsfc = reshape(GFS_Coupling%nswsfci_cpl, shape(nswsfc))
        psfc   = reshape(GFS_Coupling%psurfi_cpl, shape(psfc))
        q2m    = reshape(GFS_Coupling%q2mi_cpl, shape(q2m))
        t2m    = reshape(GFS_Coupling%t2mi_cpl, shape(t2m))
      else
        !flake(i,j)  = max(zero, GFS_Data(nb)%Sfcprop%lakefrac(ix))
        flake = reshape(GFS_Sfcprop%lakefrac, shape(flake))
        where (flake<zero) flake = zero
        focn = reshape(GFS_Sfcprop%oceanfrac, shape(focn))
        slc = reshape(GFS_Sfcprop%slc, shape(slc))
        shfsfc = reshape(GFS_Coupling%ushfsfci, shape(shfsfc))
        if (GFS_Control%lsm == GFS_Control%lsm_ruc) then
          swet = reshape(GFS_Sfcprop%wetness, shape(swet))
        else
          swet = reshape(GFS_IntDiag%wet1, shape(swet))
        end if
      end if

      ! -- zero out accumulated fields
      if (.not. GFS_control%cplflx .and. .not. GFS_control%cpllnd) then
! DH* I think we can make this a lot simpler?
!$OMP parallel do default (none) &
!$OMP             shared  (nj, ni, Atm_block, GFS_control, GFS_coupling) &
!$OMP             private (j, jb, i, ib, nb, ix, im)
        do j = 1, nj
          jb = j + Atm_block%jsc - 1
          do i = 1, ni
            ib = i + Atm_block%isc - 1
            nb = Atm_block%blkno(ib,jb)
            ix = Atm_block%ixp(ib,jb)
            im = GFS_Control%chunk_begin(nb)+ix-1
            GFS_coupling%rainc_cpl(im) = zero
            GFS_coupling%rain_cpl(im)  = zero
            GFS_coupling%snow_cpl(im)  = zero
          enddo
        enddo
! *DH
      end if

      if (GFS_control%debug) then
        ! -- diagnostics
        write(6,'("update_atmos: prsi   - min/max/avg",3g16.6)') minval(prsi),   maxval(prsi),   sum(prsi)/size(prsi)
        write(6,'("update_atmos: phii   - min/max/avg",3g16.6)') minval(phii),   maxval(phii),   sum(phii)/size(phii)
        write(6,'("update_atmos: prsl   - min/max/avg",3g16.6)') minval(prsl),   maxval(prsl),   sum(prsl)/size(prsl)
        write(6,'("update_atmos: phil   - min/max/avg",3g16.6)') minval(phil),   maxval(phil),   sum(phil)/size(phil)
        write(6,'("update_atmos: tgrs   - min/max/avg",3g16.6)') minval(temp),   maxval(temp),   sum(temp)/size(temp)
        write(6,'("update_atmos: ugrs   - min/max/avg",3g16.6)') minval(ua),     maxval(ua),     sum(ua)/size(ua)
        write(6,'("update_atmos: vgrs   - min/max/avg",3g16.6)') minval(va),     maxval(va),     sum(va)/size(va)
        write(6,'("update_atmos: qgrs   - min/max/avg",3g16.6)') minval(q),      maxval(q),      sum(q)/size(q)

        write(6,'("update_atmos: hpbl   - min/max/avg",3g16.6)') minval(hpbl),   maxval(hpbl),   sum(hpbl)/size(hpbl)
        write(6,'("update_atmos: rainc  - min/max/avg",3g16.6)') minval(rainc),  maxval(rainc),  sum(rainc)/size(rainc)
        write(6,'("update_atmos: rain   - min/max/avg",3g16.6)') minval(rain),   maxval(rain),   sum(rain)/size(rain)
        write(6,'("update_atmos: slmsk  - min/max/avg",3g16.6)') minval(slmsk),  maxval(slmsk),  sum(slmsk)/size(slmsk)
        write(6,'("update_atmos: tsfc   - min/max/avg",3g16.6)') minval(tsfc),   maxval(tsfc),   sum(tsfc)/size(tsfc)
        write(6,'("update_atmos: area   - min/max/avg",3g16.6)') minval(area),   maxval(area),   sum(area)/size(area)
        write(6,'("update_atmos: zorl   - min/max/avg",3g16.6)') minval(zorl),   maxval(zorl),   sum(zorl)/size(zorl)
        write(6,'("update_atmos: cldfra - min/max/avg",3g16.6)') minval(cldfra), maxval(cldfra), sum(cldfra)/size(cldfra)
        write(6,'("update_atmos: fice   - min/max/avg",3g16.6)') minval(fice),   maxval(fice),   sum(fice)/size(fice)
        write(6,'("update_atmos: pfils  - min/max/avg",3g16.6)') minval(pfils),  maxval(pfils),  sum(pfils)/size(pfils)
        write(6,'("update_atmos: pflls  - min/max/avg",3g16.6)') minval(pflls),  maxval(pflls),  sum(pflls)/size(pflls)
        write(6,'("update_atmos: u10m   - min/max/avg",3g16.6)') minval(u10m),   maxval(u10m),   sum(u10m)/size(u10m)
        write(6,'("update_atmos: v10m   - min/max/avg",3g16.6)') minval(v10m),   maxval(v10m),   sum(v10m)/size(v10m)
        if (GFS_Control%cplaqm) then
          write(6,'("update_atmos: canopy - min/max/avg",3g16.6)') minval(canopy), maxval(canopy), sum(canopy)/size(canopy)
          write(6,'("update_atmos: cmm    - min/max/avg",3g16.6)') minval(cmm),    maxval(cmm),    sum(cmm)/size(cmm)
          write(6,'("update_atmos: dqsfc  - min/max/avg",3g16.6)') minval(dqsfc),  maxval(dqsfc),  sum(dqsfc)/size(dqsfc)
          write(6,'("update_atmos: dtsfc  - min/max/avg",3g16.6)') minval(dtsfc),  maxval(dtsfc),  sum(dtsfc)/size(dtsfc)
          write(6,'("update_atmos: nswsfc - min/max/avg",3g16.6)') minval(nswsfc), maxval(nswsfc), sum(nswsfc)/size(nswsfc)
          write(6,'("update_atmos: oro    - min/max/avg",3g16.6)') minval(oro),    maxval(oro),    sum(oro)/size(oro)
          write(6,'("update_atmos: psfc   - min/max/avg",3g16.6)') minval(psfc),   maxval(psfc),   sum(psfc)/size(psfc)
          write(6,'("update_atmos: q2m    - min/max/avg",3g16.6)') minval(q2m),    maxval(q2m),    sum(q2m)/size(q2m)
          write(6,'("update_atmos: rca    - min/max/avg",3g16.6)') minval(rca),    maxval(rca),    sum(rca)/size(rca)
          write(6,'("update_atmos: smc    - min/max/avg",3g16.6)') minval(smc),    maxval(smc),    sum(smc)/size(smc)
          write(6,'("update_atmos: stc    - min/max/avg",3g16.6)') minval(stc),    maxval(stc),    sum(stc)/size(stc)
          write(6,'("update_atmos: t2m    - min/max/avg",3g16.6)') minval(t2m),    maxval(t2m),    sum(t2m)/size(t2m)
          write(6,'("update_atmos: vfrac  - min/max/avg",3g16.6)') minval(vfrac),  maxval(vfrac),  sum(vfrac)/size(vfrac)
          write(6,'("update_atmos: xlai   - min/max/avg",3g16.6)') minval(xlai),   maxval(xlai),   sum(xlai)/size(xlai)
          write(6,'("update_atmos: stype  - min/max/avg",3g16.6)') minval(stype),  maxval(stype),  sum(stype)/size(stype)
          write(6,'("update_atmos: vtype  - min/max/avg",3g16.6)') minval(vtype),  maxval(vtype),  sum(vtype)/size(vtype)
        else
          write(6,'("update_atmos: flake  - min/max/avg",3g16.6)') minval(flake),  maxval(flake),  sum(flake)/size(flake)
          write(6,'("update_atmos: focn   - min/max/avg",3g16.6)') minval(focn),   maxval(focn),   sum(focn)/size(focn)
          write(6,'("update_atmos: shfsfc - min/max/avg",3g16.6)') minval(shfsfc), maxval(shfsfc), sum(shfsfc)/size(shfsfc)
          write(6,'("update_atmos: slc    - min/max/avg",3g16.6)') minval(slc),    maxval(slc),    sum(slc)/size(slc)
          write(6,'("update_atmos: swet   - min/max/avg",3g16.6)') minval(swet),   maxval(swet),   sum(swet)/size(swet)
        end if
      end if

    case default
      ! -- do nothing
  end select

end subroutine update_atmos_chemistry

!> @brief Assigns imported data from coupled components to atmospheric model variables
!>
!> @param[in] atmtime      Current model time
!> @param[in] atmtimestep  Model timestep
!> @param[in] isregional   Flag for regional configuration
!> @param[in] ngrids       The number of grids
!> @param[out] rc          Return code
subroutine assign_importdata(atmtime,atmtimestep,isregional,ngrids,rc)

  use ESMF
  use module_cplfields,  only: importFields, nImportFields, queryImportFields, importFieldsValid

  implicit none
  type(time_type), intent(in)  :: atmtime, atmtimestep
  logical,         intent(in)  :: isregional
  integer,         intent(in)  :: ngrids
  integer ,        intent(out) :: rc

  !--- local variables
  integer :: n, j, i, k, ix, nb, im, isc, iec, jsc, jec, nk, dimCount
  character(len=128) :: impfield_name, fldname
  type(ESMF_TypeKind_Flag)                           :: datatype
  real(kind=ESMF_KIND_R8),  dimension(:,:), pointer  :: datar82d
  real(kind=ESMF_KIND_R8),  dimension(:,:,:), pointer:: datar83d
  real(kind=GFS_kind_phys), dimension(:,:), pointer  :: dataptr
  logical,                  dimension(:,:), pointer  :: mergeflg
  real(kind=GFS_kind_phys)                           :: tem, ofrac
  logical :: lcpl_fice
  real(ESMF_KIND_R8), parameter :: missing_value = 9.99e20_ESMF_KIND_R8

  type(ESMF_Grid)               :: grid
  type(ESMF_FieldBundle)        :: FBcpl2phys
  type(ESMF_Field)              :: dbgField
  logical                       :: add2FB
  character(len=128)            :: fname
  character(19)                 :: timestring
  character(len=:), allocatable :: fieldlist(:)
  integer                       :: nfields
  integer                       :: iyear, imonth, iday, ihour, iminute, isecond

  real(kind=GFS_kind_phys), parameter :: z0ice=1.0        !< ice roughness (cm)
  real(kind=GFS_kind_phys), parameter :: himax = 1.0e12   !< maximum ice thickness allowed
  real(kind=GFS_kind_phys), parameter :: hsmax = 1.0e12   !< maximum snow depth (m) allowed
  real(kind=GFS_kind_phys), parameter :: con_sbc = 5.670400e-8_GFS_kind_phys !< stefan-boltzmann
  !------------------------------------------------------------------------------

  rc  = -999
  ! configurations with nests cannot create debug FBs in this routine
  if (ngrids > 1 .and. GFS_control%cpl_imp_dbg) then
    print '(A)','cpl_imp_dbg=.T. is incompatible with ngrids>1'
    return
  endif

  ! set up local dimension
  isc = GFS_control%isc
  iec = GFS_control%isc+GFS_control%nx-1
  jsc = GFS_control%jsc
  jec = GFS_control%jsc+GFS_control%ny-1
  nk  = Atm_block%npz
  lcpl_fice = .false.

  allocate(dataptr(isc:iec,jsc:jec))
  allocate(mergeflg(isc:iec,jsc:jec))

  if (GFS_control%cpl_imp_dbg) then
    FBcpl2phys = ESMF_FieldBundleCreate(rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    do n = 1,nImportFields
      if (.not. ESMF_FieldIsCreated(importFields(n), rc=rc)) then
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
        cycle
      endif
      call ESMF_FieldGet(importFields(n), grid=grid, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      exit
    end do
  endif

  do n=1,nImportFields ! Each import field is only available if it was connected in the import state.

    if (.not. ESMF_FieldIsCreated(importFields(n), rc=rc)) then
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      cycle
    endif

    ! put the data from local cubed sphere grid to column grid for phys
    add2FB = .false.
    dataptr = -99999.0
    mergeflg = .false.
    call ESMF_FieldGet(importFields(n), dimCount=dimCount ,typekind=datatype, name=impfield_name, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    if ( dimCount == 2) then
      if ( datatype == ESMF_TYPEKIND_R8) then
        call ESMF_FieldGet(importFields(n),farrayPtr=datar82d,localDE=0, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
        dataptr = datar82d
        if (GFS_control%cpl_imp_mrg) then
          mergeflg(:,:) = datar82d(:,:).eq.missing_value
        endif
        if (mpp_pe() == mpp_root_pe() .and. debug) print '(A,3g16.7)','in cplIMP,atmos gets '//trim(impfield_name) &
             //' dataptr= ', dataptr(isc,jsc), maxval(dataptr), minval(dataptr)
      endif
    else if( dimCount == 3) then
      if ( datatype == ESMF_TYPEKIND_R8) then
        call ESMF_FieldGet(importFields(n),farrayPtr=datar83d,localDE=0, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      endif
    endif

    if (dataptr(isc,jsc) > -99998.0) then

      if(GFS_control%cplwav2atm) then
        ! get sea-state dependent surface roughness
        !----------------------------
        fldname = 'wave_z0_roughness_length'
        if (trim(impfield_name) == trim(fldname)) then
          if (importFieldsValid(queryImportFields(fldname))) then
            add2FB = .true.
            call copy2block(GFS_Sfcprop%zorlwav, dataptr, mask=GFS_Sfcprop%oceanfrac)
            if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'get wave roughness from mediator'
          endif
        endif
      endif ! GFS_control%cplwav2atm

      if (GFS_control%cplocn2atm) then
        ! get sst:  sst needs to be adjusted by land sea mask before passing to fv3
        !--------------------------------------------------------------------------
        fldname = 'sea_surface_temperature'
        if (trim(impfield_name) == trim(fldname)) then
          if (importFieldsValid(queryImportFields(fldname))) then
            add2FB = .true.
            call copy2block(GFS_Sfcprop%tsfco, dataptr, mask=GFS_Sfcprop%oceanfrac, validmin=150.0_GFS_kind_phys)
            if (GFS_control%cpl_imp_mrg) then
              call merge_importfield(GFS_Sfcprop%tsfco, GFS_Sfcprop%tsfc, mergeflg, mask=GFS_Sfcprop%oceanfrac)
            end if
            if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'get sst from mediator'
          endif
        end if
        ! get zonal ocean current:
        !--------------------------------------------------------------------------
        fldname = 'ocn_current_zonal'
        if (trim(impfield_name) == trim(fldname)) then
          if (importFieldsValid(queryImportFields(fldname))) then
            add2FB = .true.
            call copy2block(GFS_Sfcprop%usfco, dataptr, mask=GFS_Sfcprop%oceanfrac)
            if (GFS_control%cpl_imp_mrg) then
              call merge_importfield(GFS_Sfcprop%usfco, zero, mergeflg, mask=GFS_Sfcprop%oceanfrac)
            end if
            if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'get usfco from mediator'
          end if
        end if
        ! get meridional ocean current:
        !--------------------------------------------------------------------------
        fldname = 'ocn_current_merid'
        if (trim(impfield_name) == trim(fldname)) then
          if (importFieldsValid(queryImportFields(fldname))) then
            add2FB = .true.
            call copy2block(GFS_Sfcprop%vsfco, dataptr, mask=GFS_Sfcprop%oceanfrac)
            if (GFS_control%cpl_imp_mrg) then
              call merge_importfield(GFS_Sfcprop%vsfco, zero, mergeflg, mask=GFS_Sfcprop%oceanfrac)
            end if
            if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'get vsfco from mediator'
          end if
        end if
      end if ! GFS_control%cplocn2atm

      ! get sea ice surface temperature
      !--------------------------------
      fldname = 'sea_ice_surface_temperature'
      if (trim(impfield_name) == trim(fldname)) then
        if (importFieldsValid(queryImportFields(fldname))) then
          add2FB = .true.
          call copy2block(GFS_Sfcprop%tisfc, dataptr, mask=GFS_Sfcprop%oceanfrac, validmin=150.0_GFS_kind_phys)
          if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'get sea ice surface temperature from mediator'
        endif
      endif
      ! get sea ice fraction:  fice or sea ice concentration from the mediator
      !-----------------------------------------------------------------------
      fldname = 'ice_fraction'
      if (trim(impfield_name) == trim(fldname)) then
        if (importFieldsValid(queryImportFields(fldname))) then
          lcpl_fice = .true.
          add2FB = .true.
          call copy2block(GFS_Sfcprop%fice, dataptr, mask=GFS_Sfcprop%oceanfrac)
          if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get fice from mediator'
        endif
      endif
      ! get upward LW flux: for sea ice covered area
      !----------------------------------------------
      fldname = 'lwup_flx_ice'
      if (trim(impfield_name) == trim(fldname)) then
        if (importFieldsValid(queryImportFields(fldname))) then
          add2FB = .true.
          call copy2block(GFS_Coupling%ulwsfcin_cpl, dataptr, mask=GFS_Sfcprop%oceanfrac, flipsign=.true.)
          if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get lwflx from mediator'
        endif
      endif
      ! get latent heat flux: for sea ice covered area
      !------------------------------------------------
      fldname = 'laten_heat_flx_atm_into_ice'
      if (trim(impfield_name) == trim(fldname)) then
        if (importFieldsValid(queryImportFields(fldname))) then
          call copy2block(GFS_Coupling%dqsfcin_cpl, dataptr, mask=GFS_Sfcprop%oceanfrac, flipsign=.true.)
          if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get laten_heat from mediator'
        endif
      endif
      ! get sensible heat flux: for sea ice covered area
      !--------------------------------------------------
      fldname = 'sensi_heat_flx_atm_into_ice'
      if (trim(impfield_name) == trim(fldname)) then
        if (importFieldsValid(queryImportFields(fldname))) then
          call copy2block(GFS_Coupling%dtsfcin_cpl, dataptr, mask=GFS_Sfcprop%oceanfrac, flipsign=.true.)
          if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get sensi_heat from mediator'
        endif
      endif
      ! get zonal compt of momentum flux: for sea ice covered area
      !------------------------------------------------------------
      fldname = 'stress_on_air_ice_zonal'
      if (trim(impfield_name) == trim(fldname)) then
        if (importFieldsValid(queryImportFields(fldname))) then
          call copy2block(GFS_Coupling%dusfcin_cpl, dataptr, mask=GFS_Sfcprop%oceanfrac, flipsign=.true.)
          if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get zonal_moment_flx from mediator'
        endif
      endif
      ! get meridional compt of momentum flux: for sea ice covered area
      !-----------------------------------------------------------------
      fldname = 'stress_on_air_ice_merid'
      if (trim(impfield_name) == trim(fldname)) then
        if (importFieldsValid(queryImportFields(fldname))) then
          call copy2block(GFS_Coupling%dvsfcin_cpl, dataptr, mask=GFS_Sfcprop%oceanfrac, flipsign=.true.)
          if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get merid_moment_flx from mediator'
        endif
      endif
      ! get sea ice volume: for sea ice covered area
      !----------------------------------------------
      fldname = 'sea_ice_volume'
      if (trim(impfield_name) == trim(fldname)) then
        if (importFieldsValid(queryImportFields(fldname))) then
          add2FB = .true.
          call copy2block(GFS_Sfcprop%hice, dataptr, mask=GFS_Sfcprop%oceanfrac, validmax=himax)
          if (mpp_pe() == mpp_root_pe() .and. debug) print *,'fv3 assign_import: get ice_volume from mediator'
        endif
      endif
      ! get snow volume: for sea ice covered area
      !-------------------------------------------
      fldname = 'snow_volume_on_sea_ice'
      if (trim(impfield_name) == trim(fldname)) then
        if (importFieldsValid(queryImportFields(fldname))) then
          add2FB = .true.
          call copy2block(GFS_Coupling%hsnoin_cpl, dataptr, mask=GFS_Sfcprop%oceanfrac)
          if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get snow_volume from mediator'
        endif
      endif

      if (GFS_control%use_cice_alb) then
        ! get instantaneous near IR albedo for diffuse radiation: for sea ice covered area
        !---------------------------------------------------------------------------------
        fldname = 'inst_ice_ir_dif_albedo'
        if (trim(impfield_name) == trim(fldname)) then
          if (importFieldsValid(queryImportFields(fldname))) then
            call copy2block(GFS_Sfcprop%albdifnir_ice, dataptr, mask=GFS_Sfcprop%oceanfrac)
            if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get sfc_alb_nir_dif_cpl from mediator'
          endif
        endif
        ! get instantaneous near IR albedo for direct radiation: for sea ice covered area
        !---------------------------------------------------------------------------------
        fldname = 'inst_ice_ir_dir_albedo'
        if (trim(impfield_name) == trim(fldname)) then
          if (importFieldsValid(queryImportFields(fldname))) then
            call copy2block(GFS_Sfcprop%albdirnir_ice, dataptr, mask=GFS_Sfcprop%oceanfrac)
            if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get sfc_alb_nir_dir_cpl from mediator'
          endif
        endif
        ! get instantaneous visible albedo for diffuse radiation: for sea ice covered area
        !---------------------------------------------------------------------------------
        fldname = 'inst_ice_vis_dif_albedo'
        if (trim(impfield_name) == trim(fldname)) then
          if (importFieldsValid(queryImportFields(fldname))) then
            call copy2block(GFS_Sfcprop%albdifvis_ice, dataptr, mask=GFS_Sfcprop%oceanfrac)
            if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get sfc_alb_vis_dif_cpl from mediator'
          endif
        endif
        ! get instantaneous visible IR albedo for direct radiation: for sea ice covered area
        !---------------------------------------------------------------------------------
        fldname = 'inst_ice_vis_dir_albedo'
        if (trim(impfield_name) == trim(fldname)) then
          if (importFieldsValid(queryImportFields(fldname))) then
            call copy2block(GFS_Sfcprop%albdirvis_ice, dataptr, mask=GFS_Sfcprop%oceanfrac)
            if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get inst_ice_vis_dir_albedo from mediator'
          endif
        endif
      endif ! GFS_control%use_cice_alb

      if (GFS_control%use_med_flux) then
        ! get upward LW flux: for open ocean
        !----------------------------------------------
        fldname = 'lwup_flx_ocn'
        if (trim(impfield_name) == trim(fldname)) then
          if (importFieldsValid(queryImportFields(fldname))) then
            call copy2block(GFS_Coupling%ulwsfcin_med, dataptr, mask=GFS_Sfcprop%oceanfrac, flipsign=.true.)
            if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get lwflx for open ocean from mediator'
          endif
        endif
        ! get latent heat flux: for open ocean
        !------------------------------------------------
        fldname = 'laten_heat_flx_atm_into_ocn'
        if (trim(impfield_name) == trim(fldname)) then
          if (importFieldsValid(queryImportFields(fldname))) then
            call copy2block(GFS_Coupling%dqsfcin_med, dataptr, mask=GFS_Sfcprop%oceanfrac, flipsign=.true.)
            if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get laten_heat for open ocean from mediator'
          endif
        endif
        ! get sensible heat flux: for open ocean
        !--------------------------------------------------
        fldname = 'sensi_heat_flx_atm_into_ocn'
        if (trim(impfield_name) == trim(fldname)) then
          if (importFieldsValid(queryImportFields(fldname))) then
            call copy2block(GFS_Coupling%dtsfcin_med, dataptr, mask=GFS_Sfcprop%oceanfrac, flipsign=.true.)
            if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get sensi_heat for open ocean from mediator'
          endif
        endif
        ! get zonal compt of momentum flux: for open ocean
        !------------------------------------------------------------
        fldname = 'stress_on_air_ocn_zonal'
        if (trim(impfield_name) == trim(fldname)) then
          if (importFieldsValid(queryImportFields(fldname))) then
            call copy2block(GFS_Coupling%dusfcin_med, dataptr, mask=GFS_Sfcprop%oceanfrac, flipsign=.true.)
            if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get zonal_moment_flx for open ocean from mediator'
          endif
        endif
        ! get meridional compt of momentum flux: for open ocean
        !-----------------------------------------------------------------
        fldname = 'stress_on_air_ocn_merid'
        if (trim(impfield_name) == trim(fldname)) then
          if (importFieldsValid(queryImportFields(fldname))) then
            call copy2block(GFS_Coupling%dvsfcin_med, dataptr, mask=GFS_Sfcprop%oceanfrac, flipsign=.true.)
            if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get merid_moment_flx for open ocean from mediator'
          endif
        endif
      end if ! GFS_control%use_med_flux

      if (GFS_control%cpllnd .and. GFS_control%cpllnd2atm) then
        ! get surface snow area fraction: over land
        !------------------------------------------------
        fldname = 'inst_snow_area_fraction_lnd'
        if (trim(impfield_name) == trim(fldname)) then
          if (importFieldsValid(queryImportFields(fldname))) then
            call copy2block(GFS_Coupling%sncovr1_lnd, dataptr, mask=GFS_Sfcprop%landfrac)
            if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get snow area fraction from land'
          endif
        endif
        ! get latent heat flux: over land
        !------------------------------------------------
        fldname = 'inst_laten_heat_flx_lnd'
        if (trim(impfield_name) == trim(fldname)) then
          if (importFieldsValid(queryImportFields(fldname))) then
            call copy2block(GFS_Coupling%evap_lnd, dataptr, mask=GFS_Sfcprop%landfrac)
            if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get latent heat flux from land'
          endif
        endif
        ! get sensible heat flux: over land
        !--------------------------------------------------
        fldname = 'inst_sensi_heat_flx_lnd'
        if (trim(impfield_name) == trim(fldname)) then
          if (importFieldsValid(queryImportFields(fldname))) then
            call copy2block(GFS_Coupling%hflx_lnd, dataptr, mask=GFS_Sfcprop%landfrac)
            if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get sensible heat flux from land'
          endif
        endif
        ! get surface upward potential latent heat flux: over land
        !------------------------------------------------
        fldname = 'inst_potential_laten_heat_flx_lnd'
        if (trim(impfield_name) == trim(fldname)) then
          if (importFieldsValid(queryImportFields(fldname))) then
            call copy2block(GFS_Coupling%ep_lnd, dataptr, mask=GFS_Sfcprop%landfrac)
            if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get potential latent heat flux from land'
          endif
        endif
        ! get 2m air temperature: over land
        !------------------------------------------------
        fldname = 'inst_temp_height2m_lnd'
        if (trim(impfield_name) == trim(fldname)) then
          if (importFieldsValid(queryImportFields(fldname))) then
            call copy2block(GFS_Coupling%t2mmp_lnd, dataptr, mask=GFS_Sfcprop%landfrac)
            if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get temperature at 2m from land'
          endif
        endif
        ! get 2m specific humidity: over land
        !------------------------------------------------
        fldname = 'inst_spec_humid_height2m_lnd'
        if (trim(impfield_name) == trim(fldname)) then
          if (importFieldsValid(queryImportFields(fldname))) then
            call copy2block(GFS_Coupling%q2mp_lnd, dataptr, mask=GFS_Sfcprop%landfrac)
            if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get specific humidity at 2m from land'
          endif
        endif
        ! get specific humidity: over land
        !------------------------------------------------
        fldname = 'inst_spec_humid_lnd'
        if (trim(impfield_name) == trim(fldname)) then
          if (importFieldsValid(queryImportFields(fldname))) then
            call copy2block(GFS_Coupling%qsurf_lnd, dataptr, mask=GFS_Sfcprop%landfrac)
            if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get specific humidity from land'
          endif
        endif
        ! get upward heat flux in soil
        !------------------------------------------------
        fldname = 'inst_upward_heat_flux_lnd'
        if (trim(impfield_name) == trim(fldname)) then
          if (importFieldsValid(queryImportFields(fldname))) then
            call copy2block(GFS_Coupling%gflux_lnd, dataptr, mask=GFS_Sfcprop%landfrac)
            if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get upward heat flux from land'
          endif
        endif
        ! get surface runoff in soil
        !------------------------------------------------
        fldname = 'inst_runoff_rate_lnd'
        if (trim(impfield_name) == trim(fldname)) then
          if (importFieldsValid(queryImportFields(fldname))) then
            call copy2block(GFS_Coupling%runoff_lnd, dataptr, mask=GFS_Sfcprop%landfrac)
            if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get surface runoff from land'
          endif
        endif
        ! get subsurface runoff in soil
        !------------------------------------------------
        fldname = 'inst_subsurface_runoff_rate_lnd'
        if (trim(impfield_name) == trim(fldname)) then
          if (importFieldsValid(queryImportFields(fldname))) then
            call copy2block(GFS_Coupling%drain_lnd, dataptr, mask=GFS_Sfcprop%landfrac)
            if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get subsurface runoff from land'
          endif
        endif
        ! get momentum exchange coefficient
        !------------------------------------------------
        fldname = 'inst_drag_wind_speed_for_momentum'
        if (trim(impfield_name) == trim(fldname)) then
          if (importFieldsValid(queryImportFields(fldname))) then
            call copy2block(GFS_Coupling%cmm_lnd, dataptr, mask=GFS_Sfcprop%landfrac)
            if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get drag wind speed for momentum from land'
          endif
        endif
        ! get thermal exchange coefficient
        !------------------------------------------------
        fldname = 'inst_drag_mass_flux_for_heat_and_moisture'
        if (trim(impfield_name) == trim(fldname)) then
          if (importFieldsValid(queryImportFields(fldname))) then
            call copy2block(GFS_Coupling%chh_lnd, dataptr, mask=GFS_Sfcprop%landfrac)
            if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get thermal exchange coefficient form land'
          endif
        endif
        ! get function of surface roughness length and green vegetation fraction
        !------------------------------------------------
        fldname = 'inst_func_of_roughness_length_and_vfrac'
        if (trim(impfield_name) == trim(fldname)) then
          if (importFieldsValid(queryImportFields(fldname))) then
            call copy2block(GFS_Coupling%zvfun_lnd, dataptr, mask=GFS_Sfcprop%landfrac)
            if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get func. of roughness length and vfrac form land'
          endif
        endif
      endif ! GFS_control%cpllnd .and. GFS_control%cpllnd2atm

    endif ! if (dataptr(isc,jsc) > -99999.0) then

    if (GFS_control%cpl_fire) then
      ! get kinematic surface upward sensible heat flux of fire from Fire Behaviour model
      !------------------------------------------------
      fldname = 'hflx_fire'
      if (trim(impfield_name) == trim(fldname)) then
        if (importFieldsValid(queryImportFields(fldname))) then
          call copy2block(GFS_sfcprop%hflx_fire, dataptr)
          if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get hflx_fire from FBH model'
        endif
      endif
      ! get kinematic surface upward latent heat flux of fire from Fire Behaviour model
      !------------------------------------------------
      fldname = 'evap_fire'
      if (trim(impfield_name) == trim(fldname)) then
        if (importFieldsValid(queryImportFields(fldname))) then
          call copy2block(GFS_sfcprop%evap_fire, dataptr)
          if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get evap_fire from FBH model'
        endif
      endif
      ! get smoke_fire from Fire Behaviour model
      !------------------------------------------------
      fldname = 'smoke_fire'
      if (trim(impfield_name) == trim(fldname)) then
        if (importFieldsValid(queryImportFields(fldname))) then
          call copy2block(GFS_sfcprop%smoke_fire, dataptr)
          if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get smoke_fire from FBH model'
        endif
      endif
    endif ! (GFS_control%cpl_fire)

    if (GFS_control%cpl_imp_dbg .and. add2FB) then
      dbgField = ESMF_FieldCreate(grid=grid, typekind=ESMF_TYPEKIND_R8, name=trim(impfield_name), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      call ESMF_FieldBundleAdd(FBcpl2phys, (/dbgField/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      call ESMF_LogWrite('field '//trim(impfield_name)//' added to FBcpl2phys', ESMF_LOGMSG_INFO)
    endif
  enddo
  deallocate(mergeflg)
  deallocate(dataptr)

  !add fields not present in importstate to debug FB
  if (GFS_control%cpl_imp_dbg) then
    allocate(fieldlist(4), source=[character(len=14) :: 'ocean_fraction', 'slimskin_cpl', 'slmsk', 'zorlw'])
    do n = 1,4
      dbgField = ESMF_FieldCreate(grid=grid, typekind=ESMF_TYPEKIND_R8, name=trim(fieldlist(n)), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      call ESMF_FieldBundleAdd(FBcpl2phys, (/dbgField/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      call ESMF_LogWrite('field '//trim(fieldlist(n))//' added to FBcpl2phys', ESMF_LOGMSG_INFO)
    enddo
    deallocate(fieldlist)
  endif

  !$omp parallel do default(shared) private(i,j,nb,ix,tem,im,ofrac)
  do j=jsc,jec
    do i=isc,iec
      nb = Atm_block%blkno(i,j)
      ix = Atm_block%ixp(i,j)
      im = GFS_control%chunk_begin(nb)+ix-1

      if (GFS_control%cplwav2atm) then
        if (GFS_Sfcprop%oceanfrac(im) > zero .and. GFS_Sfcprop%zorlwav(im) > zorlmin) then
          tem = 100.0_GFS_kind_phys * min(0.1_GFS_kind_phys, GFS_Sfcprop%zorlwav(im))
          GFS_Sfcprop%zorlwav(im)      = tem
          GFS_Sfcprop%zorlw(im)        = tem
        else
          GFS_Sfcprop%zorlwav(im) = -999.0_GFS_kind_phys
        endif
      endif

      if (lcpl_fice) then
        GFS_Coupling%slimskin_cpl(im) = GFS_Sfcprop%slmsk(im)
        ofrac = GFS_Sfcprop%oceanfrac(im)

        if (ofrac > zero) then
          GFS_Sfcprop%fice(im) = max(zero, min(one, GFS_Sfcprop%fice(im)/ofrac)) !LHS: ice frac wrt water area
          if (GFS_Sfcprop%fice(im) >= GFS_control%min_seaice) then
            if (GFS_Sfcprop%fice(im) > one-epsln) GFS_Sfcprop%fice(im) = one

            if (abs(one-ofrac) < epsln) GFS_Sfcprop%slmsk(im) = 2.0_GFS_kind_phys !slmsk=2 crashes in gcycle on partial land points
            GFS_Coupling%slimskin_cpl(im) = 4.0_GFS_kind_phys

            GFS_Coupling%hsnoin_cpl(im) = min(hsmax, GFS_Coupling%hsnoin_cpl(im) / GFS_Sfcprop%fice(im))
            GFS_Sfcprop%zorli(im)       = z0ice

            tem = GFS_Sfcprop%tisfc(im) * GFS_Sfcprop%tisfc(im)
            tem = con_sbc * tem * tem
            if (GFS_Coupling%ulwsfcin_cpl(im) > zero) then
              GFS_Sfcprop%emis_ice(im) = GFS_Coupling%ulwsfcin_cpl(im) / tem
              GFS_Sfcprop%emis_ice(im) = max(0.9, min(one, GFS_Sfcprop%emis_ice(im)))
            else
              GFS_Sfcprop%emis_ice(im) = 0.96
            endif
            GFS_Coupling%ulwsfcin_cpl(im) = tem * GFS_Sfcprop%emis_ice(im)
          else
            GFS_Sfcprop%tisfc(im)       = GFS_Sfcprop%tsfco(im)
            GFS_Sfcprop%fice(im)        = zero
            GFS_Sfcprop%hice(im)        = zero
            GFS_Coupling%hsnoin_cpl(im) = zero
            !
            GFS_Coupling%dtsfcin_cpl(im)  = -99999.0 ! over open water - should not be used in ATM
            GFS_Coupling%dqsfcin_cpl(im)  = -99999.0 !                 ,,
            GFS_Coupling%dusfcin_cpl(im)  = -99999.0 !                 ,,
            GFS_Coupling%dvsfcin_cpl(im)  = -99999.0 !                 ,,
            GFS_Coupling%dtsfcin_cpl(im)  = -99999.0 !                 ,,
            GFS_Coupling%ulwsfcin_cpl(im) = -99999.0 !                 ,,
            if (abs(one-GFS_Sfcprop%oceanfrac(im)) < epsln) then !  100% open water
              GFS_Coupling%slimskin_cpl(im) = zero
              GFS_Sfcprop%slmsk(im)         = zero
            endif
          endif ! GFS_Sfcprop%fice(im) >= GFS_control%min_seaice
        endif ! GFS_Sfcprop%oceanfrac(im) > zero
      endif ! lcpl_fice
    enddo
  enddo

  if (GFS_control%cpl_imp_dbg) then
    call get_date(atmtime+atmtimestep,iyear,imonth,iday,ihour,iminute,isecond)
    write(timestring, "(I4.4,'-',I2.2,'-',I2.2,'T',I2.2,':',I2.2,':',I2.2)") iyear,imonth,iday,ihour,iminute,isecond
    if (isregional) then
      fname = 'fv3_merge_'//trim(timestring)//'.nc'
    else
      fname = 'fv3_merge_'//trim(timestring)//'.tile*.nc'
    end if

    call ESMF_FieldBundleGet(FBcpl2phys, fieldCount=nfields, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    if (nfields > 0) then
      call write_FB(FBcpl2phys, trim(fname), nfields, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      call ESMF_FieldBundleDestroy(FBcpl2phys, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    endif
  endif
  rc=0

end subroutine assign_importdata
!
subroutine setup_inlinedata(fieldName, datar82d, logunit)

  use ESMF, only: ESMF_KIND_R8, ESMF_LogFoundError, ESMF_LOGERR_PASSTHRU

  !--- arguments
  character(len=*), intent(in) :: fieldName
  real(kind=ESMF_KIND_R8), dimension(:,:), target, intent(in) :: datar82d
  integer, intent(in) :: logunit

  !--- local variables
  real(kind=GFS_kind_phys), dimension(:,:), pointer  :: dataptr

  allocate(dataptr(size(datar82d,1),size(datar82d,2)))
  dataptr = datar82d

  ! fill variables
  select case(trim(fieldName))
  case ('Si_ifrac')
    call copy2block(GFS_Coupling%fice_dat, dataptr)
  case ('Si_thick')
    call copy2block(GFS_Coupling%hice_dat, dataptr)
  case ('So_omask')
    call copy2block(GFS_Coupling%mask_dat, dataptr)
  case ('So_t')
    call copy2block(GFS_Coupling%tsfco_dat, dataptr)
  case ('Si_t')
    call copy2block(GFS_Coupling%tice_dat, dataptr)
  case default
    write(logunit,*) trim(fieldName)//' can not be used by cdeps inline! Skipping field ...'
  end select

end subroutine setup_inlinedata
!
  subroutine setup_exportdata(rc)

    use ESMF

    use module_cplfields,  only: exportFields, chemistryFieldNames
    use module_cplscalars, only: flds_scalar_name

    !--- arguments
    integer, optional, intent(out) :: rc

    !--- local variables
    integer                :: i, j, ix, im
    integer                :: isc, iec, jsc, jec
    integer                :: nb, nk
    real(GFS_kind_phys)    :: rtime, rtimek, spval

    integer                                     :: localrc
    integer                                     :: n,rank
    logical                                     :: isFound
    type(ESMF_TypeKind_Flag)                    :: datatype
    character(len=ESMF_MAXSTR)                  :: fieldName
    real(kind=ESMF_KIND_R4), dimension(:,:), pointer   :: datar42d
    real(kind=ESMF_KIND_R8), dimension(:,:), pointer   :: datar82d
    real(kind=ESMF_KIND_R8), dimension(:,:,:), pointer :: datar83d

    !--- local parameters
    real(kind=ESMF_KIND_R8), parameter :: zeror8 = 0._ESMF_KIND_R8
    real(GFS_kind_phys),     parameter :: revap  = one/2.501E+06_GFS_kind_phys ! reciprocal of specific
                                                                               ! heat of vaporization J/kg

    !--- begin
    if (present(rc)) rc = ESMF_SUCCESS

    isc = Atm_block%isc
    iec = Atm_block%iec
    jsc = Atm_block%jsc
    jec = Atm_block%jec
    nk  = Atm_block%npz

    rtime  = one / GFS_control%dtp
    rtimek = GFS_control%rho_h2o * rtime
    spval  = GFS_control%huge

    do n=1, size(exportFields)

      datar42d => null()
      datar82d => null()
      datar83d => null()

      isFound = ESMF_FieldIsCreated(exportFields(n), rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      if (isFound) then
        call ESMF_FieldGet(exportFields(n), name=fieldname, rank=rank, typekind=datatype, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
        if (trim(fieldname) == trim(flds_scalar_name)) then
          isFound = .false.
        else
          if (datatype == ESMF_TYPEKIND_R8) then
            select case (rank)
            case (2)
              call ESMF_FieldGet(exportFields(n),farrayPtr=datar82d,localDE=0, rc=localrc)
              if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
            case (3)
              call ESMF_FieldGet(exportFields(n),farrayPtr=datar83d,localDE=0, rc=localrc)
              if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
            case default
              !--- skip field
              isFound = .false.
            end select
          else if (datatype == ESMF_TYPEKIND_R4) then
            select case (rank)
            case (2)
              call ESMF_FieldGet(exportFields(n),farrayPtr=datar42d,localDE=0, rc=localrc)
              if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
            case default
              !--- skip field
              isFound = .false.
            end select
          else
            !--- skip field
            isFound = .false.
          end if
        end if
      end if
      !--- skip field if only required for chemistry
      if (isFound .and. GFS_control%cplchm) isFound = .not.any(trim(fieldname) == chemistryFieldNames)

      if (isFound) then
!$omp parallel do default(shared) private(nb) reduction(max:localrc)
        do nb = 1, Atm_block%nblks
          select case (trim(fieldname))
            !--- Instantaneous quantities
            ! Instantaneous mean layer pressure (Pa)
            case ('inst_pres_levels')
              call block_data_copy_or_fill(datar83d, GFS_statein%prsl, zeror8, Atm_block, nb, offset=GFS_Control%chunk_begin(nb), rc=localrc)
            ! Instantaneous geopotential at model layer centers (m2 s-2)
            case ('inst_geop_levels')
              call block_data_copy_or_fill(datar83d, GFS_statein%phil, zeror8, Atm_block, nb, offset=GFS_Control%chunk_begin(nb), rc=localrc)
            ! Instantaneous zonal wind (m s-1)
            case ('inst_zonal_wind_levels')
              call block_data_copy_or_fill(datar83d, GFS_statein%ugrs, zeror8, Atm_block, nb, offset=GFS_Control%chunk_begin(nb), rc=localrc)
            ! Instantaneous meridional wind (m s-1)
            case ('inst_merid_wind_levels')
              call block_data_copy_or_fill(datar83d, GFS_statein%vgrs, zeror8, Atm_block, nb, offset=GFS_Control%chunk_begin(nb), rc=localrc)
            ! Instantaneous surface roughness length (cm)
            case ('inst_surface_roughness')
              call block_data_copy(datar82d, GFS_sfcprop%zorl, Atm_block, nb, offset=GFS_Control%chunk_begin(nb), rc=localrc)
            ! Instantaneous u wind (m/s) 10 m above ground
            case ('inst_zonal_wind_height10m')
              call block_data_copy(datar82d, GFS_coupling%u10mi_cpl, Atm_block, nb, offset=GFS_Control%chunk_begin(nb), rc=localrc)
            ! Instantaneous v wind (m/s) 10 m above ground
            case ('inst_merid_wind_height10m')
              call block_data_copy(datar82d, GFS_coupling%v10mi_cpl, Atm_block, nb, offset=GFS_Control%chunk_begin(nb), rc=localrc)
            ! Instantaneous Zonal compt of momentum flux (N/m**2)
            case ('inst_zonal_moment_flx')
              call block_data_copy(datar82d, GFS_coupling%dusfci_cpl, Atm_block, nb, -one, spval, offset=GFS_Control%chunk_begin(nb), rc=localrc)
            ! Instantaneous Merid compt of momentum flux (N/m**2)
            case ('inst_merid_moment_flx')
              call block_data_copy(datar82d, GFS_coupling%dvsfci_cpl, Atm_block, nb, -one, spval, offset=GFS_Control%chunk_begin(nb), rc=localrc)
            ! Instantaneous Sensible heat flux (W/m**2)
            case ('inst_sensi_heat_flx')
              call block_data_copy(datar82d, GFS_coupling%dtsfci_cpl, Atm_block, nb, -one, spval, offset=GFS_Control%chunk_begin(nb), rc=localrc)
            ! Instantaneous Latent heat flux (W/m**2)
            case ('inst_laten_heat_flx')
              call block_data_copy(datar82d, GFS_coupling%dqsfci_cpl, Atm_block, nb, -one, spval, offset=GFS_Control%chunk_begin(nb), rc=localrc)
            ! Instantaneous Evap flux (kg/m**2/s)
            case ('inst_evap_rate')
              call block_data_copy(datar82d, GFS_coupling%dqsfci_cpl, Atm_block, nb, -revap, spval, offset=GFS_Control%chunk_begin(nb), rc=localrc)
            ! Instantaneous precipitation rate (kg/m2/s)
            case ('inst_prec_rate')
              call block_data_copy(datar82d, GFS_coupling%rain_cpl, Atm_block, nb, rtimek, spval, offset=GFS_Control%chunk_begin(nb), rc=localrc)
            ! Instantaneous convective precipitation rate (kg/m2/s)
            case ('inst_prec_rate_conv')
              call block_data_copy(datar82d, GFS_coupling%rainc_cpl, Atm_block, nb, rtimek, spval, offset=GFS_Control%chunk_begin(nb), rc=localrc)
            ! Instaneous snow precipitation rate (kg/m2/s)
            case ('inst_fprec_rate')
              call block_data_copy(datar82d, GFS_coupling%snow_cpl, Atm_block, nb, rtimek, spval, offset=GFS_Control%chunk_begin(nb), rc=localrc)
            ! Instantaneous Downward long wave radiation flux (W/m**2)
            case ('inst_down_lw_flx')
              call block_data_copy(datar82d, GFS_coupling%dlwsfci_cpl, Atm_block, nb, offset=GFS_Control%chunk_begin(nb), rc=localrc)
            ! Instantaneous Downward solar radiation flux (W/m**2)
            case ('inst_down_sw_flx')
              call block_data_copy(datar82d, GFS_coupling%dswsfci_cpl, Atm_block, nb, offset=GFS_Control%chunk_begin(nb), rc=localrc)
            ! Instantaneous Temperature (K) 2 m above ground
            case ('inst_temp_height2m')
              call block_data_copy(datar82d, GFS_coupling%t2mi_cpl, Atm_block, nb, offset=GFS_Control%chunk_begin(nb), rc=localrc)
            ! Instantaneous Specific humidity (kg/kg) 2 m above ground
            case ('inst_spec_humid_height2m')
              call block_data_copy(datar82d, GFS_coupling%q2mi_cpl, Atm_block, nb, offset=GFS_Control%chunk_begin(nb), rc=localrc)
            ! Instantaneous Temperature (K) at surface
            case ('inst_temp_height_surface')
              call block_data_copy(datar82d, GFS_coupling%tsfci_cpl, Atm_block, nb, offset=GFS_Control%chunk_begin(nb), rc=localrc)
            ! Instantaneous Pressure (Pa) land and sea surface
            case ('inst_pres_height_surface')
              call block_data_copy(datar82d, GFS_coupling%psurfi_cpl, Atm_block, nb, offset=GFS_Control%chunk_begin(nb), rc=localrc)
            ! Instantaneous Surface height (m)
            case ('inst_surface_height')
              call block_data_copy(datar82d, GFS_coupling%oro_cpl, Atm_block, nb, offset=GFS_Control%chunk_begin(nb), rc=localrc)
            ! Instantaneous NET long wave radiation flux (W/m**2)
            case ('inst_net_lw_flx')
              call block_data_copy(datar82d, GFS_coupling%nlwsfci_cpl, Atm_block, nb, offset=GFS_Control%chunk_begin(nb), rc=localrc)
            ! Instantaneous NET solar radiation flux over the ocean (W/m**2)
            case ('inst_net_sw_flx')
              call block_data_copy(datar82d, GFS_coupling%nswsfci_cpl, Atm_block, nb, offset=GFS_Control%chunk_begin(nb), rc=localrc)
            ! Instantaneous sfc downward nir direct flux (W/m**2)
            case ('inst_down_sw_ir_dir_flx')
              call block_data_copy(datar82d, GFS_coupling%dnirbmi_cpl, Atm_block, nb, offset=GFS_Control%chunk_begin(nb), rc=localrc)
            ! Instantaneous sfc downward nir diffused flux (W/m**2)
            case ('inst_down_sw_ir_dif_flx')
              call block_data_copy(datar82d, GFS_coupling%dnirdfi_cpl, Atm_block, nb, offset=GFS_Control%chunk_begin(nb), rc=localrc)
            ! Instantaneous sfc downward uv+vis direct flux (W/m**2)
            case ('inst_down_sw_vis_dir_flx')
              call block_data_copy(datar82d, GFS_coupling%dvisbmi_cpl, Atm_block, nb, offset=GFS_Control%chunk_begin(nb), rc=localrc)
            ! Instantaneous sfc downward uv+vis diffused flux (W/m**2)
            case ('inst_down_sw_vis_dif_flx')
              call block_data_copy(datar82d, GFS_coupling%dvisdfi_cpl, Atm_block, nb, offset=GFS_Control%chunk_begin(nb), rc=localrc)
            ! Instantaneous net sfc nir direct flux (W/m**2)
            case ('inst_net_sw_ir_dir_flx')
              call block_data_copy(datar82d, GFS_coupling%nnirbmi_cpl, Atm_block, nb, offset=GFS_Control%chunk_begin(nb), rc=localrc)
            ! Instantaneous net sfc nir diffused flux (W/m**2)
            case ('inst_net_sw_ir_dif_flx')
              call block_data_copy(datar82d, GFS_coupling%nnirdfi_cpl, Atm_block, nb, offset=GFS_Control%chunk_begin(nb), rc=localrc)
            ! Instantaneous net sfc uv+vis direct flux (W/m**2)
            case ('inst_net_sw_vis_dir_flx')
              call block_data_copy(datar82d, GFS_coupling%nvisbmi_cpl, Atm_block, nb, offset=GFS_Control%chunk_begin(nb), rc=localrc)
            ! Instantaneous net sfc uv+vis diffused flux (W/m**2)
            case ('inst_net_sw_vis_dif_flx')
              call block_data_copy(datar82d, GFS_coupling%nvisdfi_cpl, Atm_block, nb, offset=GFS_Control%chunk_begin(nb), rc=localrc)
            ! Land/Sea mask (sea:0,land:1)
            case ('inst_land_sea_mask', 'slmsk')
              call block_data_copy(datar82d, GFS_sfcprop%slmsk, Atm_block, nb, offset=GFS_Control%chunk_begin(nb), rc=localrc)
            ! Total precipitation amount in each time step
            case ('inst_rainfall_amount')
              call block_data_copy(datar82d, GFS_sfcprop%tprcp, Atm_block, nb, offset=GFS_Control%chunk_begin(nb), rc=localrc)
            !--- Mean quantities
            ! MEAN Zonal compt of momentum flux (N/m**2)
            case ('mean_zonal_moment_flx_atm')
              call block_data_copy(datar82d, GFS_coupling%dusfc_cpl, Atm_block, nb, -rtime, spval, offset=GFS_Control%chunk_begin(nb), rc=localrc)
            ! MEAN Merid compt of momentum flux (N/m**2)
            case ('mean_merid_moment_flx_atm')
              call block_data_copy(datar82d, GFS_coupling%dvsfc_cpl, Atm_block, nb, -rtime, spval, offset=GFS_Control%chunk_begin(nb), rc=localrc)
            ! MEAN Sensible heat flux (W/m**2)
            case ('mean_sensi_heat_flx')
              call block_data_copy(datar82d, GFS_coupling%dtsfc_cpl, Atm_block, nb, -rtime, spval, offset=GFS_Control%chunk_begin(nb), rc=localrc)
            ! MEAN Latent heat flux (W/m**2)
            case ('mean_laten_heat_flx')
              call block_data_copy(datar82d, GFS_coupling%dqsfc_cpl, Atm_block, nb, -rtime, spval, offset=GFS_Control%chunk_begin(nb), rc=localrc)
            ! MEAN Evap rate (kg/m**2/s)
            case ('mean_evap_rate')
              call block_data_copy(datar82d, GFS_coupling%dqsfc_cpl, Atm_block, nb, -rtime*revap, spval, offset=GFS_Control%chunk_begin(nb), rc=localrc)
            ! MEAN Downward LW heat flux (W/m**2)
            case ('mean_down_lw_flx')
              call block_data_copy(datar82d, GFS_coupling%dlwsfc_cpl, Atm_block, nb, rtime, spval, offset=GFS_Control%chunk_begin(nb), rc=localrc)
            ! MEAN Downward SW heat flux (W/m**2)
            case ('mean_down_sw_flx')
              call block_data_copy(datar82d, GFS_coupling%dswsfc_cpl, Atm_block, nb, rtime, spval, offset=GFS_Control%chunk_begin(nb), rc=localrc)
            ! MEAN NET long wave radiation flux (W/m**2)
            case ('mean_net_lw_flx')
              call block_data_copy(datar82d, GFS_coupling%nlwsfc_cpl, Atm_block, nb, rtime, spval, offset=GFS_Control%chunk_begin(nb), rc=localrc)
            ! MEAN NET solar radiation flux over the ocean (W/m**2)
            case ('mean_net_sw_flx')
              call block_data_copy(datar82d, GFS_coupling%nswsfc_cpl, Atm_block, nb, rtime, spval, offset=GFS_Control%chunk_begin(nb), rc=localrc)
            ! MEAN sfc downward nir direct flux (W/m**2)
            case ('mean_down_sw_ir_dir_flx')
              call block_data_copy(datar82d, GFS_coupling%dnirbm_cpl, Atm_block, nb, rtime, spval, offset=GFS_Control%chunk_begin(nb), rc=localrc)
            ! MEAN sfc downward nir diffused flux (W/m**2)
            case ('mean_down_sw_ir_dif_flx')
              call block_data_copy(datar82d, GFS_coupling%dnirdf_cpl, Atm_block, nb, rtime, spval, offset=GFS_Control%chunk_begin(nb), rc=localrc)
            ! MEAN sfc downward uv+vis direct flux (W/m**2)
            case ('mean_down_sw_vis_dir_flx')
              call block_data_copy(datar82d, GFS_coupling%dvisbm_cpl, Atm_block, nb, rtime, spval, offset=GFS_Control%chunk_begin(nb), rc=localrc)
            ! MEAN sfc downward uv+vis diffused flux (W/m**2)
            case ('mean_down_sw_vis_dif_flx')
              call block_data_copy(datar82d, GFS_coupling%dvisdf_cpl, Atm_block, nb, rtime, spval, offset=GFS_Control%chunk_begin(nb), rc=localrc)
            ! MEAN NET sfc nir direct flux (W/m**2)
            case ('mean_net_sw_ir_dir_flx')
              call block_data_copy(datar82d, GFS_coupling%nnirbm_cpl, Atm_block, nb, rtime, spval, offset=GFS_Control%chunk_begin(nb), rc=localrc)
            ! MEAN NET sfc nir diffused flux (W/m**2)
            case ('mean_net_sw_ir_dif_flx')
              call block_data_copy(datar82d, GFS_coupling%nnirdf_cpl, Atm_block, nb, rtime, spval, offset=GFS_Control%chunk_begin(nb), rc=localrc)
            ! MEAN NET sfc uv+vis direct flux (W/m**2)
            case ('mean_net_sw_vis_dir_flx')
              call block_data_copy(datar82d, GFS_coupling%nvisbm_cpl, Atm_block, nb, rtime, spval, offset=GFS_Control%chunk_begin(nb), rc=localrc)
            ! MEAN NET sfc uv+vis diffused flux (W/m**2)
            case ('mean_net_sw_vis_dif_flx')
              call block_data_copy(datar82d, GFS_coupling%nvisdf_cpl, Atm_block, nb, rtime, spval, offset=GFS_Control%chunk_begin(nb), rc=localrc)
            ! MEAN precipitation rate (kg/m2/s)
            case ('mean_prec_rate')
              call block_data_copy(datar82d, GFS_sfcprop%tprcp, Atm_block, nb, rtimek, spval, offset=GFS_Control%chunk_begin(nb), rc=localrc)
            ! MEAN convective precipitation rate (kg/m2/s)
            case ('mean_prec_rate_conv')
              call block_data_copy(datar82d, GFS_coupling%rainc_cpl, Atm_block, nb, rtimek, spval, offset=GFS_Control%chunk_begin(nb), rc=localrc)
            ! MEAN snow precipitation rate (kg/m2/s)
            case ('mean_fprec_rate')
              call block_data_copy(datar82d, GFS_coupling%snow_cpl, Atm_block, nb, rtimek, spval, offset=GFS_Control%chunk_begin(nb), rc=localrc)
            ! oceanfrac used by atm to calculate fluxes
            case ('openwater_frac_in_atm')
              call block_data_combine_fractions(datar82d, GFS_sfcprop%oceanfrac, GFS_sfcprop%fice, Atm_block, nb, offset=GFS_Control%chunk_begin(nb), rc=localrc)
            !--- Dycore quantities
            ! bottom layer temperature (t)
            case('inst_temp_height_lowest')
              call block_data_copy_or_fill(datar82d, DYCORE_data(nb)%coupling%t_bot, zeror8, Atm_block, nb, offset=1, rc=localrc)
            case('inst_temp_height_lowest_from_phys')
              call block_data_copy_or_fill(datar82d, GFS_Statein%tgrs, 1, zeror8, Atm_block, nb, offset=GFS_Control%chunk_begin(nb), rc=localrc)
            ! bottom layer specific humidity (q)
            !    !    ! CHECK if tracer 1 is for specific humidity     !    !    !
            case('inst_spec_humid_height_lowest')
              call block_data_copy_or_fill(datar82d, DYCORE_data(nb)%coupling%tr_bot, 1, zeror8, Atm_block, nb, offset=1, rc=localrc)
            case('inst_spec_humid_height_lowest_from_phys')
              call block_data_copy_or_fill(datar82d, GFS_Statein%qgrs, 1, GFS_Control%ntqv, zeror8, Atm_block, nb, offset=GFS_Control%chunk_begin(nb), rc=localrc)
            ! bottom layer zonal wind (u)
            case('inst_zonal_wind_height_lowest')
              call block_data_copy_or_fill(datar82d, DYCORE_data(nb)%coupling%u_bot, zeror8, Atm_block, nb, offset=1, rc=localrc)
            ! bottom layer meridional wind (v)
            case('inst_merid_wind_height_lowest')
              call block_data_copy_or_fill(datar82d, DYCORE_data(nb)%coupling%v_bot, zeror8, Atm_block, nb, offset=1, rc=localrc)
            ! surface friction velocity
            case('surface_friction_velocity')
              call block_data_copy_or_fill(datar82d, GFS_Sfcprop%uustar, zeror8, Atm_block, nb, offset=GFS_Control%chunk_begin(nb), rc=localrc)
            ! bottom layer pressure (p)
            case('inst_pres_height_lowest')
              call block_data_copy_or_fill(datar82d, DYCORE_data(nb)%coupling%p_bot, zeror8, Atm_block, nb, offset=1, rc=localrc)
            ! bottom layer pressure (p) from physics
            case('inst_pres_height_lowest_from_phys')
              call block_data_copy_or_fill(datar82d, GFS_Statein%prsl, 1, zeror8, Atm_block, nb, offset=GFS_Control%chunk_begin(nb), rc=localrc)
            ! dimensionless exner function at surface adjacent layer
            case('inst_exner_function_height_lowest')
              call block_data_copy_or_fill(datar82d, GFS_Statein%prslk, 1, zeror8, Atm_block, nb, offset=GFS_Control%chunk_begin(nb), rc=localrc)
            ! bottom layer height (z)
            case('inst_height_lowest')
              call block_data_copy_or_fill(datar82d, DYCORE_data(nb)%coupling%z_bot, zeror8, Atm_block, nb, offset=1, rc=localrc)
            case ('vfrac')
              call block_data_copy(datar82d, GFS_Sfcprop%vfrac, Atm_block, nb, offset=GFS_Control%chunk_begin(nb), rc=localrc)
            case ('zorl')
              call block_data_copy(datar82d, GFS_Sfcprop%zorl, Atm_block, nb, offset=GFS_Control%chunk_begin(nb), rc=localrc)
            case default
              localrc = ESMF_RC_NOT_FOUND
          end select
        enddo
        if (ESMF_LogFoundError(rcToCheck=localrc, msg="Failure to populate exported field: "//trim(fieldname), &
          line=__LINE__, file=__FILE__, rcToReturn=rc)) return
      endif
    enddo ! exportFields

!---
    if (GFS_control%cplflx) then
! zero out accumulated fields
!$omp parallel do default(shared) private(i,j,nb,ix,im)
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          im = GFS_control%chunk_begin(nb)+ix-1
          GFS_coupling%dusfc_cpl(im)  = zero
          GFS_coupling%dvsfc_cpl(im)  = zero
          GFS_coupling%dtsfc_cpl(im)  = zero
          GFS_coupling%dqsfc_cpl(im)  = zero
          GFS_coupling%nlwsfc_cpl(im) = zero
          GFS_coupling%dnirbm_cpl(im) = zero
          GFS_coupling%dnirdf_cpl(im) = zero
          GFS_coupling%dvisbm_cpl(im) = zero
          GFS_coupling%dvisdf_cpl(im) = zero
        enddo
      enddo
      if (mpp_pe() == mpp_root_pe()) print *,'zeroing coupling accumulated fields at kdt= ',GFS_control%kdt
    endif !cplflx
!---
    if (GFS_control%cplflx .or. GFS_control%cpllnd) then
! zero out accumulated fields
!$omp parallel do default(shared) private(i,j,nb,ix,im)
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          im = GFS_control%chunk_begin(nb)+ix-1
          GFS_coupling%dlwsfc_cpl(im) = zero
          GFS_coupling%dswsfc_cpl(im) = zero
          GFS_coupling%rain_cpl(im)   = zero
          GFS_coupling%rainc_cpl(im)  = zero
          GFS_coupling%snow_cpl(im)   = zero
          GFS_coupling%nswsfc_cpl(im) = zero
          GFS_coupling%nnirbm_cpl(im) = zero
          GFS_coupling%nnirdf_cpl(im) = zero
          GFS_coupling%nvisbm_cpl(im) = zero
          GFS_coupling%nvisdf_cpl(im) = zero
        enddo
      enddo
      if (mpp_pe() == mpp_root_pe()) print *,'zeroing coupling accumulated fields at kdt= ',GFS_control%kdt
    endif !cplflx or cpllnd

  end subroutine setup_exportdata

!> @brief Adds land-sea mask information to the grid object
!>
!> @param[in] fcstGrid Grid object
!> @param[out] rc Return code
  subroutine addLsmask2grid(fcstGrid, rc)

    use ESMF
!
    implicit none
    type(ESMF_Grid)      :: fcstGrid
    integer, optional, intent(out) :: rc
!
!  local vars
    integer isc, iec, jsc, jec
    integer i, j, nb, ix, im
!    integer CLbnd(2), CUbnd(2), CCount(2), TLbnd(2), TUbnd(2), TCount(2)
    integer, allocatable  :: lsmask(:,:)
    integer(kind=ESMF_KIND_I4), pointer  :: maskPtr(:,:)
!
    isc = GFS_control%isc
    iec = GFS_control%isc+GFS_control%nx-1
    jsc = GFS_control%jsc
    jec = GFS_control%jsc+GFS_control%ny-1
    allocate(lsmask(isc:iec,jsc:jec))
!
!$omp parallel do default(shared) private(i,j,nb,ix,im)
    do j=jsc,jec
      do i=isc,iec
        nb = Atm_block%blkno(i,j)
        ix = Atm_block%ixp(i,j)
        im = GFS_control%chunk_begin(nb)+ix-1
! use land sea mask: land:1, ocean:0
        lsmask(i,j) = floor(one + epsln - GFS_sfcprop%oceanfrac(im))
      enddo
    enddo
!
! Get mask
    call ESMF_GridAddItem(fcstGrid, itemflag=ESMF_GRIDITEM_MASK,   &
                          staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

!    call ESMF_GridGetItemBounds(fcstGrid, itemflag=ESMF_GRIDITEM_MASK,   &
!         staggerloc=ESMF_STAGGERLOC_CENTER, computationalLBound=ClBnd,  &
!         computationalUBound=CUbnd, computationalCount=Ccount,  &
!         totalLBound=TLbnd, totalUBound=TUbnd, totalCount=Tcount, rc=rc)
!    print *,'in set up grid, aft add esmfgridadd item mask, rc=',rc, &
!     'ClBnd=',ClBnd,'CUbnd=',CUbnd,'Ccount=',Ccount, &
!     'TlBnd=',TlBnd,'TUbnd=',TUbnd,'Tcount=',Tcount
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_GridGetItem(fcstGrid, itemflag=ESMF_GRIDITEM_MASK,   &
                          staggerloc=ESMF_STAGGERLOC_CENTER,farrayPtr=maskPtr, rc=rc)
!    print *,'in set up grid, aft get maskptr, rc=',rc, 'size=',size(maskPtr,1),size(maskPtr,2), &
!      'bound(maskPtr)=', LBOUND(maskPtr,1),LBOUND(maskPtr,2),UBOUND(maskPtr,1),UBOUND(maskPtr,2)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
!
!$omp parallel do default(shared) private(i,j)
    do j=jsc,jec
      do i=isc,iec
        maskPtr(i-isc+1,j-jsc+1) = lsmask(i,j)
      enddo
    enddo
!      print *,'in set set lsmask, maskPtr=', maxval(maskPtr), minval(maskPtr)
!
    deallocate(lsmask)

  end subroutine addLsmask2grid

!> @brief Retrieves domain information from grid
!>
!> @param[in] n Grid number
!> @param[out] layout Processor grid layout
!> @param[out] nx Grid points in x-direction
!> @param[out] ny Grid points in y-direction
!> @param[out] pelist List of processor IDs
subroutine atmos_model_get_nth_domain_info(n, layout, nx, ny, pelist)
  integer, intent(in)  :: n
  integer, intent(out) :: layout(2)
  integer, intent(out) :: nx, ny
  integer, pointer, intent(out) :: pelist(:)

  call get_nth_domain_info(n, layout, nx, ny, pelist)

end subroutine atmos_model_get_nth_domain_info
!> @brief Copy data from a 2D source array to a 1D block-distributed destination array.
!>
!> @param[inout] destin_ptr 1D destination array distributed by block
!> @param[in]    source_ptr 2D source array in local coordinates
!> @param[in]    mask       Optional 1D mask array; elements with mask <= 0 are skipped
!> @param[in]    validmin   Optional minimum acceptable value; values <= validmin are rejected
!> @param[in]    validmax   Optional maximum acceptable value; values >= validmax are rejected
!> @param[in]    flipsign   Optional flag to negate values before assignment
!> @param[in]    block      Optional block control structure; uses Atm_block if not provided
!>
!> @author Denise.Worthen@noaa.gov
subroutine copy2block(destin_ptr, source_ptr, mask, validmin, validmax, flipsign, block)

  real(kind=GFS_kind_phys), intent(inout), target :: destin_ptr(:)
  real(kind=GFS_kind_phys), intent(in),    target :: source_ptr(:,:)
  real(kind=GFS_kind_phys), intent(in),  target, optional :: mask(:)
  type(block_control_type), intent(in),  target, optional :: block
  real(kind=GFS_kind_phys), intent(in), optional :: validmin
  real(kind=GFS_kind_phys), intent(in), optional :: validmax
  logical,                  intent(in), optional :: flipsign

  integer :: isc, jsc, iec, jec
  integer :: i, j, nb, ix, im
  type(block_control_type), pointer :: active_block

  real(kind=GFS_kind_phys) :: fval, spval
  real(kind=GFS_kind_phys) :: lvmin, lvmax
  logical :: lflip

  if (present(block)) then
    active_block => block
  else
    active_block => Atm_block
  end if

  isc = GFS_control%isc
  iec = GFS_control%isc + GFS_control%nx - 1
  jsc = GFS_control%jsc
  jec = GFS_control%jsc + GFS_control%ny - 1
  spval = GFS_control%huge

  lvmin = -spval
  lvmax = spval
  lflip = .false.
  if(present(validmin)) lvmin = validmin
  if(present(validmax)) lvmax = validmax
  if(present(flipsign)) lflip = flipsign

  !$omp parallel do default(shared) private(i,j,nb,ix,im,fval)
  do j = jsc, jec
    do i = isc, iec
      nb = active_block%blkno(i,j)
      ix = active_block%ixp(i,j)
      im = GFS_control%chunk_begin(nb) + ix - 1
      if (present(mask)) then
        if (mask(im) <= zero) cycle
      end if
      fval = source_ptr(i-isc+1, j-jsc+1)
      if (lvmin /= -spval .and. fval <= lvmin) cycle
      if (lvmax /=  spval .and. fval >= lvmax) cycle
      if (lflip) then
        destin_ptr(im) = -fval
      else
        destin_ptr(im) = fval
      end if
    end do
  end do
end subroutine copy2block
!> @brief Merge values from a source field into a destination array based on a merge flag.
!>
!> @param[inout] destin_ptr   1D destination array distributed by block
!> @param[in]    source_ptr   1D source array distributed by block
!> @param[in]    mergeflg     2D logical merge flag array in local coordinates
!> @param[in]    mask         Optional 1D mask array; elements with mask <= 0 are skipped
!> @param[in]    block        Optional block control structure; uses Atm_block if not provided
!>
!> @author Denise.Worthen@noaa.gov
subroutine merge_importfield_with_field(destin_ptr, source_ptr, mergeflg, mask, block)

  real(kind=GFS_kind_phys), intent(inout), target :: destin_ptr(:)
  real(kind=GFS_kind_phys), intent(in),    target :: source_ptr(:)
  logical,                  intent(in),    target :: mergeflg(:,:)
  real(kind=GFS_kind_phys), intent(in),    target, optional :: mask(:)
  type(block_control_type), intent(in),    target, optional :: block

  real(kind=GFS_kind_phys) :: fval

  integer :: isc, jsc, iec, jec
  integer :: i, j, nb, ix, im
  type(block_control_type), pointer :: active_block

  if (present(block)) then
    active_block => block
  else
    active_block => Atm_block
  end if

  isc = GFS_control%isc
  iec = GFS_control%isc + GFS_control%nx - 1
  jsc = GFS_control%jsc
  jec = GFS_control%jsc + GFS_control%ny - 1

  !$omp parallel do default(shared) private(i,j,nb,ix,im,fval)
  do j = jsc, jec
    do i = isc, iec
      nb = active_block%blkno(i,j)
      ix = active_block%ixp(i,j)
      im = GFS_control%chunk_begin(nb) + ix - 1
      fval = source_ptr(im)
      if (present(mask)) then
        if (mask(im) <= zero) cycle
      end if
      if (mergeflg(i-isc+1, j-jsc+1)) then
        destin_ptr(im) = fval
      end if
    end do
  end do
end subroutine merge_importfield_with_field
!> @brief Merge scalar values into a destination array based on a merge flag.
!>
!> @param[inout] destin_ptr   1D destination array distributed by block
!> @param[in]    scalarfill   Scalar value to assign where merge flag is true
!> @param[in]    mergeflg     2D logical merge flag array in local coordinates
!> @param[in]    mask         Optional 1D mask array; elements with mask <= 0 are skipped
!> @param[in]    block        Optional block control structure; uses Atm_block if not provided
!>
!> @author Denise.Worthen@noaa.gov
subroutine merge_importfield_with_scalar(destin_ptr, scalarfill, mergeflg, mask, block)

  real(kind=GFS_kind_phys), intent(inout), target :: destin_ptr(:)
  real(kind=GFS_kind_phys), intent(in)            :: scalarfill
  logical,                  intent(in),    target :: mergeflg(:,:)
  real(kind=GFS_kind_phys), intent(in),    target, optional :: mask(:)
  type(block_control_type), intent(in),    target, optional :: block

  integer :: isc, jsc, iec, jec
  integer :: i, j, nb, ix, im
  type(block_control_type), pointer :: active_block

  if (present(block)) then
    active_block => block
  else
    active_block => Atm_block
  end if

  isc = GFS_control%isc
  iec = GFS_control%isc + GFS_control%nx - 1
  jsc = GFS_control%jsc
  jec = GFS_control%jsc + GFS_control%ny - 1

  !$omp parallel do default(shared) private(i,j,nb,ix,im)
  do j = jsc, jec
    do i = isc, iec
      nb = active_block%blkno(i,j)
      ix = active_block%ixp(i,j)
      im = GFS_control%chunk_begin(nb) + ix - 1
      if (present(mask)) then
        if (mask(im) <= zero) cycle
      end if
      if (mergeflg(i-isc+1, j-jsc+1)) then
        destin_ptr(im) = scalarfill
      end if
    end do
  end do
end subroutine merge_importfield_with_scalar
!> @brief Write a FB for debugging the transfer of field data from importState to physics
!>
!> @param[inout] FBcpl2phys   an ESMF FieldBundle
!> @param[in]    fname        the name of the file where the FB is written
!> @param[in]    nfields      the number of fields in the FB
!> @param[out]   rc           return code
!>
!> @author Denise.Worthen@noaa.gov
subroutine write_FB(FBcpl2phys,fname,nfields,rc)

  use ESMF
  type(ESMF_FieldBundle), intent(inout) :: FBcpl2phys
  character(len=*),       intent(in)    :: fname
  integer,                intent(in)    :: nfields
  integer,                intent(out)   :: rc
  !--- local variables
  type(ESMF_Field)                 :: dbgField
  real(kind=ESMF_KIND_R8), pointer :: dbgptr(:,:)
  character(len=128), allocatable  :: fieldlist(:)
  integer :: n,nb,localrc

  rc = ESMF_SUCCESS
  nullify(dbgptr)
  allocate(character(len=128) :: fieldlist(1:nfields))

  call ESMF_FieldBundleGet(FBcpl2phys, fieldNameList=fieldlist, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
  do n = 1,nfields
    call ESMF_FieldBundleGet(FBcpl2phys, fieldName=trim(fieldlist(n)), field=dbgField, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_FieldGet(dbgField, farrayPtr=dbgptr, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_LogWrite('field '//trim(fieldlist(n))//' retrieved from FBcpl2phys', ESMF_LOGMSG_INFO)

    dbgptr = -GFS_control%huge
    localrc = ESMF_SUCCESS
    select case(trim(fieldlist(n)))
    case ('wave_z0_roughness_length')
      !$omp parallel do default(shared) private(nb) reduction(max:localrc)
      do nb = 1, Atm_block%nblks
        call block_data_copy(dbgptr, GFS_Sfcprop%zorlwav, Atm_block, nb, offset=GFS_Control%chunk_begin(nb), rc=localrc)
      enddo
    case ('sea_ice_surface_temperature')
      !$omp parallel do default(shared) private(nb) reduction(max:localrc)
      do nb = 1, Atm_block%nblks
        call block_data_copy(dbgptr, GFS_Sfcprop%tisfc, Atm_block, nb, offset=GFS_Control%chunk_begin(nb), rc=localrc)
      enddo
    case ('sea_surface_temperature')
      !$omp parallel do default(shared) private(nb) reduction(max:localrc)
      do nb = 1, Atm_block%nblks
        call block_data_copy(dbgptr, GFS_Sfcprop%tsfco, Atm_block, nb, offset=GFS_Control%chunk_begin(nb), rc=localrc)
      enddo
    case ('ocn_current_zonal')
      !$omp parallel do default(shared) private(nb) reduction(max:localrc)
      do nb = 1, Atm_block%nblks
        call block_data_copy(dbgptr, GFS_Sfcprop%usfco, Atm_block, nb, offset=GFS_Control%chunk_begin(nb), rc=localrc)
      enddo
    case ('ocn_current_merid')
      !$omp parallel do default(shared) private(nb) reduction(max:localrc)
      do nb = 1, Atm_block%nblks
        call block_data_copy(dbgptr, GFS_Sfcprop%vsfco, Atm_block, nb, offset=GFS_Control%chunk_begin(nb), rc=localrc)
      enddo
    case ('ice_fraction')
      !$omp parallel do default(shared) private(nb) reduction(max:localrc)
      do nb = 1, Atm_block%nblks
        call block_data_copy(dbgptr, GFS_Sfcprop%fice, Atm_block, nb, offset=GFS_Control%chunk_begin(nb), rc=localrc)
      enddo
    case ('lwup_flx_ice')
      !$omp parallel do default(shared) private(nb) reduction(max:localrc)
      do nb = 1, Atm_block%nblks
        call block_data_copy(dbgptr, GFS_Coupling%ulwsfcin_cpl, Atm_block, nb, offset=GFS_Control%chunk_begin(nb), rc=localrc)
      enddo
    case ('sea_ice_volume')
      !$omp parallel do default(shared) private(nb) reduction(max:localrc)
      do nb = 1, Atm_block%nblks
        call block_data_copy(dbgptr, GFS_Sfcprop%hice, Atm_block, nb, offset=GFS_Control%chunk_begin(nb), rc=localrc)
      enddo
    case ('snow_volume_on_sea_ice')
      !$omp parallel do default(shared) private(nb) reduction(max:localrc)
      do nb = 1, Atm_block%nblks
        call block_data_copy(dbgptr, GFS_Coupling%hsnoin_cpl, Atm_block, nb, offset=GFS_Control%chunk_begin(nb), rc=localrc)
      enddo
    case ('ocean_fraction')
      !$omp parallel do default(shared) private(nb) reduction(max:localrc)
      do nb = 1, Atm_block%nblks
        call block_data_copy(dbgptr, GFS_Sfcprop%oceanfrac, Atm_block, nb, offset=GFS_Control%chunk_begin(nb), rc=localrc)
      enddo
    case ('slimskin_cpl')
      !$omp parallel do default(shared) private(nb) reduction(max:localrc)
      do nb = 1, Atm_block%nblks
        call block_data_copy(dbgptr, GFS_Coupling%slimskin_cpl, Atm_block, nb, offset=GFS_Control%chunk_begin(nb), rc=localrc)
      enddo
    case ('slmsk')
      !$omp parallel do default(shared) private(nb) reduction(max:localrc)
      do nb = 1, Atm_block%nblks
        call block_data_copy(dbgptr, GFS_Sfcprop%slmsk, Atm_block, nb, offset=GFS_Control%chunk_begin(nb), rc=localrc)
      enddo
    case ('zorlw')
      !$omp parallel do default(shared) private(nb) reduction(max:localrc)
      do nb = 1, Atm_block%nblks
        call block_data_copy(dbgptr, GFS_Sfcprop%zorlw, Atm_block, nb, offset=GFS_Control%chunk_begin(nb), rc=localrc)
      enddo
    case default
      localrc = ESMF_RC_NOT_FOUND
    end select
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
  enddo
  rc = localrc

  call ESMF_FieldBundleWrite(FBcpl2phys, fileName=trim(fname), timeslice=1, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
  ! clean up
  do n = 1,nfields
    call ESMF_FieldBundleGet(FBcpl2phys, fieldName=trim(fieldlist(n)), field=dbgField, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_FieldDestroy(dbgField, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
  enddo
  nullify(dbgptr)
  deallocate(fieldlist)
end subroutine write_FB

end module atmos_model_mod
