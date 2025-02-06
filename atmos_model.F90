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
!-----------------------------------------------------------------------
!<OVERVIEW>
!  Driver for the atmospheric model, contains routines to advance the
!  atmospheric model state by one time step.
!</OVERVIEW>

!<DESCRIPTION>
!     This version of atmos_model_mod has been designed around the implicit
!     version diffusion scheme of the GCM. It requires two routines to advance
!     the atmospheric model one time step into the future. These two routines
!     correspond to the down and up sweeps of the standard tridiagonal solver.
!     Most atmospheric processes (dynamics,radiation,etc.) are performed
!     in the down routine. The up routine finishes the vertical diffusion
!     and computes moisture related terms (convection,large-scale condensation,
!     and precipitation).

!     The boundary variables needed by other component models for coupling
!     are contained in a derived data type. A variable of this derived type
!     is returned when initializing the atmospheric model. It is used by other
!     routines in this module and by coupling routines. The contents of
!     this derived type should only be modified by the atmospheric model.

!</DESCRIPTION>

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
use CCPP_driver,        only: CCPP_step, non_uniform_blocks

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
!-----------------------------------------------------------------------

!<PUBLICTYPE >
 type atmos_data_type
     integer                       :: axes(4)            ! axis indices (returned by diag_manager) for the atmospheric grid
                                                         ! (they correspond to the x, y, pfull, phalf axes)
     integer, pointer              :: pelist(:) =>null() ! pelist where atmosphere is running.
     integer                       :: layout(2)          ! computer task laytout
     logical                       :: regional           ! true if domain is regional
     logical                       :: nested             ! true if there is a nest
     logical                       :: moving_nest_parent ! true if this grid has a moving nest child
     logical                       :: is_moving_nest     ! true if this is a moving nest grid
     logical                       :: isAtCapTime        ! true if currTime is at the cap driverClock's currTime
     integer                       :: ngrids             !
     integer                       :: mygrid             !
     integer                       :: mlon, mlat
     integer                       :: iau_offset         ! iau running window length
     logical                       :: pe                 ! current pe.
     real(kind=GFS_kind_phys), pointer, dimension(:)     :: ak, bk
     real(kind=GFS_kind_phys), pointer, dimension(:,:)   :: lon_bnd  => null() ! local longitude axis grid box corners in radians.
     real(kind=GFS_kind_phys), pointer, dimension(:,:)   :: lat_bnd  => null() ! local latitude axis grid box corners in radians.
     real(kind=GFS_kind_phys), pointer, dimension(:,:)   :: lon      => null() ! local longitude axis grid box centers in radians.
     real(kind=GFS_kind_phys), pointer, dimension(:,:)   :: lat      => null() ! local latitude axis grid box centers in radians.
     real(kind=GFS_kind_phys), pointer, dimension(:,:)   :: dx, dy
     real(kind=GFS_kind_phys), pointer, dimension(:,:)   :: area
     real(kind=GFS_kind_phys), pointer, dimension(:,:,:) :: layer_hgt, level_hgt
     type(domain2d)                :: domain             ! domain decomposition
     type(domain2d)                :: domain_for_read    ! domain decomposition
     type(time_type)               :: Time               ! current time
     type(time_type)               :: Time_step          ! atmospheric time step.
     type(time_type)               :: Time_init          ! reference time.
     type(grid_box_type)           :: grid               ! hold grid information needed for 2nd order conservative flux exchange
     type(GFS_externaldiag_type), pointer, dimension(:) :: Diag
 end type atmos_data_type
                                                         ! to calculate gradient on cubic sphere grid.
!</PUBLICTYPE >

! these two arrays, lon_bnd_work and lat_bnd_work are 'working' arrays, always allocated
! as (nlon+1, nlat+1) and are used to get the corner lat/lon values from the dycore.
! these values are then copied to Atmos%lon_bnd, Atmos%lat_bnd which are allocated with
! sizes that correspond to the corner coordinates distgrid in fcstGrid
real(kind=GFS_kind_phys), pointer, dimension(:,:), save :: lon_bnd_work  => null()
real(kind=GFS_kind_phys), pointer, dimension(:,:), save :: lat_bnd_work  => null()
integer, save :: i_bnd_size, j_bnd_size

integer :: fv3Clock, getClock, updClock, setupClock, radClock, physClock

!-----------------------------------------------------------------------
integer :: blocksize    = 1
logical :: chksum_debug = .false.
logical :: dycore_only  = .false.
logical :: debug        = .false.
!logical :: debug        = .true.
logical :: sync         = .false.
real    :: avg_max_length=3600.
logical :: ignore_rst_cksum = .false.
namelist /atmos_model_nml/ blocksize, chksum_debug, dycore_only, debug, sync, ccpp_suite, avg_max_length, &
                           ignore_rst_cksum

type (time_type) :: diag_time, diag_time_fhzero

!--- concurrent and decoupled radiation and physics variables
!-------------------
!  DYCORE containers
!-------------------
type(DYCORE_data_type),    allocatable :: DYCORE_Data(:)  ! number of blocks

!----------------
!  GFS containers
!----------------
type(GFS_externaldiag_type), target :: GFS_Diag(DIAG_SIZE)
type(GFS_restart_type)              :: GFS_restart_var

!--------------
! IAU container
!--------------
type(iau_external_data_type)        :: IAU_Data ! number of blocks

!-----------------
!  Block container
!-----------------
type (block_control_type), target   :: Atm_block

!-----------------------------------------------------------------------

character(len=128) :: version = '$Id$'
character(len=128) :: tagname = '$Name$'

#ifdef NAM_phys
  logical,parameter :: flip_vc = .false.
#else
  logical,parameter :: flip_vc = .true.
#endif

  real(kind=GFS_kind_phys), parameter :: zero    = 0.0_GFS_kind_phys,     &
                                         one     = 1.0_GFS_kind_phys,     &
                                         epsln   = 1.0e-10_GFS_kind_phys, &
                                         zorlmin = 1.0e-7_GFS_kind_phys

contains

!#######################################################################
! <SUBROUTINE NAME="update_atmos_radiation_physics">
!
!<DESCRIPTION>
!   Called every time step as the atmospheric driver to compute the
!   atmospheric tendencies for dynamics, radiation, vertical diffusion of
!   momentum, tracers, and heat/moisture.  For heat/moisture only the
!   downward sweep of the tridiagonal elimination is performed, hence
!   the name "_down".
!</DESCRIPTION>

!   <TEMPLATE>
!     call  update_atmos_radiation_physics (Atmos)
!   </TEMPLATE>

! <INOUT NAME="Atmos" TYPE="type(atmos_data_type)">
!   Derived-type variable that contains fields needed by the flux exchange module.
!   These fields describe the atmospheric grid and are needed to
!   compute/exchange fluxes with other component models.  All fields in this
!   variable type are allocated for the global grid (without halo regions).
! </INOUT>

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
    call atmos_phys_driver_statein (GFS_Control, GFS_Statein, Atm_block, flip_vc)
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
      call CCPP_step (step="timestep_init", nblks=Atm_block%nblks, ierr=ierr)
      if (ierr/=0)  call mpp_error(FATAL, 'Call to CCPP timestep_init step failed')

      if (GFS_Control%do_sppt .or. GFS_Control%do_shum .or. GFS_Control%do_skeb .or. &
          GFS_Control%lndp_type > 0  .or. GFS_Control%do_ca .or. GFS_Control%do_spp) then
!--- call stochastic physics pattern generation / cellular automata
        call stochastic_physics_wrapper(GFS_Control, GFS_Statein, GFS_Grid, GFS_Sfcprop, GFS_Radtend, GFS_Coupling, Atm_block, ierr)
        if (ierr/=0)  call mpp_error(FATAL, 'Call to stochastic_physics_wrapper failed')
      endif

!--- if coupled, assign coupled fields
      call assign_importdata(jdat(:),rc)
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
        call CCPP_step (step="radiation", nblks=Atm_block%nblks, ierr=ierr)
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
      call CCPP_step (step="physics", nblks=Atm_block%nblks, ierr=ierr)
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
        call CCPP_step (step="stochastics", nblks=Atm_block%nblks, ierr=ierr)
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
      call CCPP_step (step="timestep_finalize", nblks=Atm_block%nblks, ierr=ierr)
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
! </SUBROUTINE>


!#######################################################################
! <SUBROUTINE NAME="atmos_timestep_diagnostics">
!
! <OVERVIEW>
! Calculates per-timestep, domain-wide, diagnostic, information and
! prints to stdout from master rank. Must be called after physics
! update but before first_time_step flag is cleared.
! </OVERVIEW>

!   <TEMPLATE>
!     call  atmos_timestep_diagnostics (Atmos)
!   </TEMPLATE>

! <INOUT NAME="Atmos" TYPE="type(atmos_data_type)">
!   Derived-type variable that contains fields needed by the flux exchange module.
!   These fields describe the atmospheric grid and are needed to
!   compute/exchange fluxes with other component models.  All fields in this
!   variable type are allocated for the global grid (without halo regions).
! </INOUT>
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
! </SUBROUTINE>

!#######################################################################
! <SUBROUTINE NAME="atmos_model_init">
!
! <OVERVIEW>
! Routine to initialize the atmospheric model
! </OVERVIEW>

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
                           Atmos%ngrids, Atmos%mygrid, Atmos%pelist)
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

   ! This logic deals with non-uniform block sizes for CCPP.
   ! When non-uniform block sizes are used, it is required
   ! that only the last block has a different (smaller)
   ! size than all other blocks. This is the standard in
   ! FV3. If this is the case, set non_uniform_blocks (a
   ! variable imported from CCPP_driver) to .true. and
   ! allocate nthreads+1 elements of the interstitial array.
   ! The extra element will be used by the thread that
   ! runs over the last, smaller block.
   if (minval(Atm_block%blksz)==maxval(Atm_block%blksz)) then
      non_uniform_blocks = .false.
      allocate(GFS_interstitial(nthrds))
   else if (all(minloc(Atm_block%blksz)==(/size(Atm_block%blksz)/))) then
      non_uniform_blocks = .true.
      allocate(GFS_interstitial(nthrds+1))
   else
      call mpp_error(FATAL, 'For non-uniform blocksizes, only the last element ' // &
                            'in Atm_block%blksz can be different from the others')
   end if

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

   ! allocate required to work around GNU compiler bug 100886 https://gcc.gnu.org/bugzilla/show_bug.cgi?id=100886
   allocate(Init_parm%input_nml_file, mold=input_nml_file)
   Init_parm%input_nml_file  => input_nml_file
   Init_parm%fn_nml='using internal file'

   call GFS_initialize (GFS_control, GFS_Statein, GFS_Stateout, GFS_Sfcprop, &
                        GFS_Coupling, GFS_Grid, GFS_Tbd, GFS_Cldprop, GFS_Radtend, &
                        GFS_Intdiag, GFS_interstitial, Init_parm)

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
   call atmos_phys_driver_statein (GFS_control, GFS_statein, Atm_block, flip_vc)

   ! When asked to calculate 3-dim. tendencies, set Stateout variables to
   ! Statein variables here in order to capture the first call to dycore
    if (GFS_control%ldiag3d) then
        GFS_Stateout%gu0 = GFS_Statein%ugrs
        GFS_Stateout%gv0 = GFS_Statein%vgrs
        GFS_Stateout%gt0 = GFS_Statein%tgrs
        GFS_Stateout%gq0 = GFS_Statein%qgrs
    endif

   ! Initialize the CCPP framework
   call CCPP_step (step="init", nblks=Atm_block%nblks, ierr=ierr)
   if (ierr/=0)  call mpp_error(FATAL, 'Call to CCPP init step failed')
   ! Initialize the CCPP physics
   call CCPP_step (step="physics_init", nblks=Atm_block%nblks, ierr=ierr)
   if (ierr/=0)  call mpp_error(FATAL, 'Call to CCPP physics_init step failed')

   if (GFS_Control%do_sppt .or. GFS_Control%do_shum .or. GFS_Control%do_skeb .or. &
       GFS_Control%lndp_type > 0  .or. GFS_Control%do_ca .or. GFS_Control%do_spp) then

!--- Initialize stochastic physics pattern generation / cellular automata for first time step
     call stochastic_physics_wrapper(GFS_control, GFS_statein, GFS_grid, GFS_sfcprop, GFS_Radtend, GFS_Coupling, Atm_block, ierr)
     if (ierr/=0)  call mpp_error(FATAL, 'Call to stochastic_physics_wrapper failed')

   endif

   !--- set the initial diagnostic timestamp
   diag_time = Time
   call get_time (Atmos%Time - Atmos%Time_init, sec)
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


!#######################################################################
! <SUBROUTINE NAME="update_atmos_model_dynamics"
!
! <OVERVIEW>
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
! </SUBROUTINE>


!#######################################################################
! <SUBROUTINE NAME="atmos_model_exchange_phase_1"
!
! <OVERVIEW>
!   Perform data exchange with coupled components in run phase 1
! </OVERVIEW>
!
! <DESCRIPTION>
!  This subroutine currently exports atmospheric fields and tracers
!  to the chemistry component during the model's run phase 1, i.e.
!  before chemistry is run.
! </DESCRIPTION>

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
! </SUBROUTINE>


!#######################################################################
! <SUBROUTINE NAME="atmos_model_exchange_phase_2"
!
! <OVERVIEW>
!   Perform data exchange with coupled components in run phase 2
! </OVERVIEW>
!
! <DESCRIPTION>
!  This subroutine currently imports fields updated by the coupled
!  chemistry component back into the atmospheric model during run
!  phase 2.
! </DESCRIPTION>

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
! </SUBROUTINE>


!#######################################################################
! <SUBROUTINE NAME="update_atmos_model_state"
!
! <OVERVIEW>
subroutine update_atmos_model_state (Atmos, rc)
! to update the model state after all concurrency is completed
  use ESMF
  type (atmos_data_type), intent(inout) :: Atmos
  integer, optional,      intent(out)   :: rc
!--- local variables
  integer :: i, localrc, sec_lastfhzerofh
  integer :: isec, seconds, isec_fhzero
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
    call atmosphere_nggps_diag(Atmos%Time,ltavg=.true.,avg_max_length=avg_max_length)
    if (ANY(nint(output_fh(:)*3600.0) == seconds) .or. (GFS_control%kdt == first_kdt)) then
      if (mpp_pe() == mpp_root_pe()) write(6,*) "---isec,seconds",isec,seconds
      time_int = real(isec)
      if(Atmos%iau_offset > zero) then
        if( time_int - Atmos%iau_offset*3600. > zero ) then
          time_int = time_int - Atmos%iau_offset*3600.
        else if(seconds == Atmos%iau_offset*3600) then
          call get_time (Atmos%Time - diag_time_fhzero, isec_fhzero)
          time_int = real(isec_fhzero)
          if (mpp_pe() == mpp_root_pe()) write(6,*) "---iseczero",isec_fhzero
        endif
      endif
      time_intfull = real(seconds)
      if(Atmos%iau_offset > zero) then
        if( time_intfull - Atmos%iau_offset*3600. > zero) then
          time_intfull = time_intfull - Atmos%iau_offset*3600.
        endif
      endif
      if (mpp_pe() == mpp_root_pe()) write(6,*) ' gfs diags time since last bucket empty: ',time_int/3600.,'hrs'
      call atmosphere_nggps_diag(Atmos%Time)
      call fv3atm_diag_output(Atmos%Time, GFS_Diag, Atm_block, GFS_control%nx, GFS_control%ny, &
                            GFS_control%levs, 1, 1, 1.0_GFS_kind_phys, time_int, time_intfull, &
                            GFS_control%fhswr, GFS_control%fhlwr)
    endif

    !---  find current fhzero
    if( GFS_Control%fhzero_array(1) > 0. ) then
      fhzero_loop: do i=1,size(GFS_Control%fhzero_array)
        tmpflag_fhzero = .false.
        if( GFS_Control%fhzero_array(i) > 0.) then
          if( i == 1 ) then
            if( seconds <= GFS_Control%fhzero_fhour(i)*3600. ) tmpflag_fhzero = .true.
          else if( i > 1 ) then
            if( seconds > GFS_Control%fhzero_fhour(i-1)*3600. .and. seconds <= GFS_Control%fhzero_fhour(i)*3600. ) &
              tmpflag_fhzero = .true.
          endif
          if( tmpflag_fhzero) then
            GFS_Control%fhzero = GFS_Control%fhzero_array(i)
            if( GFS_Control%fhzero > 0) then
              sec_lastfhzerofh = (int(seconds/3600.)/int(GFS_Control%fhzero))*int(GFS_Control%fhzero)*3600
            else
              sec_lastfhzerofh = 0
            endif
          endif
        endif
      enddo fhzero_loop
    else
      sec_lastfhzerofh = 0
    endif
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



!#######################################################################
! <SUBROUTINE NAME="atmos_model_end">
!
! <OVERVIEW>
!  termination routine for atmospheric model
! </OVERVIEW>

! <DESCRIPTION>
!  Call once to terminate this module and any other modules used.
!  This routine writes a restart file and deallocates storage
!  used by the derived-type variable atmos_boundary_data_type.
! </DESCRIPTION>

! <TEMPLATE>
!   call atmos_model_end (Atmos)
! </TEMPLATE>

! <INOUT NAME="Atmos" TYPE="type(atmos_data_type)">
!   Derived-type variable that contains fields needed by the flux exchange module.
! </INOUT>

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
    call CCPP_step (step="physics_finalize", nblks=Atm_block%nblks, ierr=ierr)
    if (ierr/=0)  call mpp_error(FATAL, 'Call to CCPP physics_finalize step failed')

!   The CCPP framework for all cdata structures is finalized in CCPP_step 'finalize'.
    call CCPP_step (step="finalize", nblks=Atm_block%nblks, ierr=ierr)
    if (ierr/=0)  call mpp_error(FATAL, 'Call to CCPP finalize step failed')

    deallocate (Atmos%lon, Atmos%lat)
    deallocate (Atmos%lon_bnd, Atmos%lat_bnd)
    deallocate (lon_bnd_work, lat_bnd_work)

end subroutine atmos_model_end

! </SUBROUTINE>
!#######################################################################
! <SUBROUTINE NAME="atmos_model_restart">
! <DESCRIPTION>
!  Write out restart files registered through register_restart_file
! </DESCRIPTION>
subroutine atmos_model_restart(Atmos, timestamp)
  use update_ca, only: write_ca_restart
  type (atmos_data_type),   intent(inout) :: Atmos
  character(len=*),  intent(in)           :: timestamp

    if (quilting_restart) then
       call fv_sfc_restart_output(GFS_sfcprop, Atm_block, GFS_control)
       call fv_phy_restart_output(GFS_restart_var, Atm_block)
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
! </SUBROUTINE>

!#######################################################################
! <SUBROUTINE NAME="get_atmos_model_ungridded_dim">
!
! <DESCRIPTION>
!  Retrieve ungridded dimensions of atmospheric model arrays
! </DESCRIPTION>

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
! </SUBROUTINE>

!#######################################################################
! <SUBROUTINE NAME="get_atmos_tracer_types">
! <DESCRIPTION>
!  Identify and return usage and type id of atmospheric tracers.
!  Ids are defined as:
!    0 = generic tracer
!    1 = chemistry - prognostic
!    2 = chemistry - diagnostic
!
!  Tracers are identified via the additional 'tracer_usage' keyword and
!  their optional 'type' qualifier. A tracer is assumed prognostic if
!  'type' is not provided. See examples from the field_table file below:
!
!  Prognostic tracer:
!  ------------------
!  "TRACER", "atmos_mod",    "so2"
!            "longname",     "so2 mixing ratio"
!            "units",        "ppm"
!            "tracer_usage", "chemistry"
!            "profile_type", "fixed", "surface_value=5.e-6" /
!
!  Diagnostic tracer:
!  ------------------
!  "TRACER", "atmos_mod",    "pm25"
!            "longname",     "PM2.5"
!            "units",        "ug/m3"
!            "tracer_usage", "chemistry", "type=diagnostic"
!            "profile_type", "fixed", "surface_value=5.e-6" /
!
!  For atmospheric chemistry, the order of both prognostic and diagnostic
!  tracers is validated against the model's internal assumptions.
!
! </DESCRIPTION>
subroutine get_atmos_tracer_types(tracer_types)

  use field_manager_mod,  only: parse
  use tracer_manager_mod, only: query_method

  integer, intent(out) :: tracer_types(:)

  !--- local variables
  logical :: found
  integer :: n, num_tracers, num_types
  integer :: id_max, id_min, id_num, ip_max, ip_min, ip_num
  character(len=32)  :: tracer_usage
  character(len=128) :: control, tracer_type

  !--- begin

  !--- validate array size
  call get_number_tracers(MODEL_ATMOS, num_tracers=num_tracers)

  if (size(tracer_types) < num_tracers) &
    call mpp_error(FATAL, 'insufficient size of tracer type array')

  !--- initialize tracer indices
  id_min = num_tracers + 1
  id_max = -id_min
  ip_min = id_min
  ip_max = id_max
  id_num = 0
  ip_num = 0

  do n = 1, num_tracers
    tracer_types(n) = 0
    found = query_method('tracer_usage',MODEL_ATMOS,n,tracer_usage,control)
    if (found) then
      if (trim(tracer_usage) == 'chemistry') then
        !--- set default to prognostic
        tracer_type = 'prognostic'
        num_types = parse(control, 'type', tracer_type)
        select case (trim(tracer_type))
          case ('diagnostic')
            tracer_types(n) = 2
            id_num = id_num + 1
            id_max = n
            if (id_num == 1) id_min = n
          case ('prognostic')
            tracer_types(n) = 1
            ip_num = ip_num + 1
            ip_max = n
            if (ip_num == 1) ip_min = n
        end select
      end if
    end if
  end do

  if (ip_num > 0) then
    !--- check if prognostic tracers are contiguous
    if (ip_num > ip_max - ip_min + 1) &
      call mpp_error(FATAL, 'prognostic chemistry tracers must be contiguous')
  end if

  if (id_num > 0) then
    !--- check if diagnostic tracers are contiguous
    if (id_num > id_max - id_min + 1) &
      call mpp_error(FATAL, 'diagnostic chemistry tracers must be contiguous')
  end if

  !--- prognostic tracers must precede diagnostic ones
  if (ip_max > id_min) &
    call mpp_error(FATAL, 'diagnostic chemistry tracers must follow prognostic ones')

end subroutine get_atmos_tracer_types
! </SUBROUTINE>

!#######################################################################
! <SUBROUTINE NAME="update_atmos_chemistry">
! <DESCRIPTION>
!  Populate exported chemistry fields with current atmospheric state
!  data (state='export'). Update tracer concentrations for atmospheric
!  chemistry with values from chemistry component (state='import').
!  Fields should be exported/imported from/to the atmospheric state
!  after physics calculations.
!
!  NOTE: It is assumed that all the chemical tracers follow the standard
!  atmospheric tracers, which end with ozone. The order of the chemical
!  tracers must match their order in the chemistry component.
!
!  Requires:
!         GFS data types
!         Atm_block
! </DESCRIPTION>
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

  real(ESMF_KIND_R8), dimension(:,:), pointer :: aod, area, canopy, cmm,  &
    dqsfc, dtsfc, fice, flake, focn, fsnow, hpbl, nswsfc, oro, psfc, &
    q2m, rain, rainc, rca, shfsfc, slmsk, stype, swet, t2m, tsfc,    &
    u10m, uustar, v10m, vfrac, xlai, zorl

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
      end if

      if (GFS_control%debug) then
        write(6,'("update_atmos: ",a,": qgrs - min/max/avg",3g16.6)') &
          trim(state), minval(q), maxval(q), sum(q)/size(q)
        if (GFS_control%cplaqm) &
          write(6,'("update_atmos: ",a,": aod  - min/max    ",3g16.6)') &
            trim(state), minval(aod), maxval(aod)
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
! </SUBROUTINE>

  subroutine assign_importdata(jdat, rc)

    use module_cplfields,  only: importFields, nImportFields, queryImportFields, &
                                 importFieldsValid
    use ESMF
!
    implicit none
    integer, intent(in)  :: jdat(8)
    integer, intent(out) :: rc

    !--- local variables
    integer :: n, j, i, k, ix, nb, im, isc, iec, jsc, jec, nk, dimCount, findex
    integer :: sphum, liq_wat, ice_wat, o3mr
    character(len=128) :: impfield_name, fldname
    type(ESMF_TypeKind_Flag)                           :: datatype
    real(kind=ESMF_KIND_R8),  dimension(:,:), pointer  :: datar82d
    real(kind=ESMF_KIND_R8),  dimension(:,:,:), pointer:: datar83d
    real(kind=GFS_kind_phys), dimension(:,:), pointer  :: datar8
    logical,                  dimension(:,:), pointer  :: mergeflg
    real(kind=GFS_kind_phys)                           :: tem, ofrac
    logical found, isFieldCreated, lcpl_fice
    real(ESMF_KIND_R8), parameter :: missing_value = 9.99e20_ESMF_KIND_R8
    type(ESMF_Grid)  :: grid
    type(ESMF_Field) :: dbgField
    character(19)    :: currtimestring
    real (kind=GFS_kind_phys), parameter :: z0ice=1.0    !  (in cm)

!
!     real(kind=GFS_kind_phys), parameter :: himax = 8.0      !< maximum ice thickness allowed
!     real(kind=GFS_kind_phys), parameter :: himin = 0.1      !< minimum ice thickness required
!     real(kind=GFS_kind_phys), parameter :: hsmax = 100.0    !< maximum snow depth (m) allowed
      real(kind=GFS_kind_phys), parameter :: himax = 1.0e12   !< maximum ice thickness allowed
      real(kind=GFS_kind_phys), parameter :: hsmax = 1.0e12   !< maximum snow depth (m) allowed
      real(kind=GFS_kind_phys), parameter :: con_sbc = 5.670400e-8_GFS_kind_phys !< stefan-boltzmann
!
!------------------------------------------------------------------------------
!
    rc  = -999

! set up local dimension
    isc = GFS_control%isc
    iec = GFS_control%isc+GFS_control%nx-1
    jsc = GFS_control%jsc
    jec = GFS_control%jsc+GFS_control%ny-1
    nk  = Atm_block%npz
    lcpl_fice = .false.

    allocate(datar8(isc:iec,jsc:jec))
    allocate(mergeflg(isc:iec,jsc:jec))

!   if (mpp_pe() == mpp_root_pe() .and. debug) print *,'in cplImp,dim=',isc,iec,jsc,jec
!   if (mpp_pe() == mpp_root_pe() .and. debug) print *,'in cplImp,GFS_data, size', size(GFS_data)
!   if (mpp_pe() == mpp_root_pe() .and. debug) print *,'in cplImp,tsfc, size', size(GFS_data(1)%sfcprop%tsfc)
!   if (mpp_pe() == mpp_root_pe() .and. debug) print *,'in cplImp,tsfc, min_seaice', GFS_control%min_seaice

    do n=1,nImportFields ! Each import field is only available if it was connected in the import state.

      found = .false.

      isFieldCreated = ESMF_FieldIsCreated(importFields(n), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      if (isFieldCreated) then ! put the data from local cubed sphere grid to column grid for phys

        datar8 = -99999.0
        mergeflg = .false.
        call ESMF_FieldGet(importFields(n), dimCount=dimCount ,typekind=datatype, &
                           name=impfield_name, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        if ( dimCount == 2) then
          if ( datatype == ESMF_TYPEKIND_R8) then
            call ESMF_FieldGet(importFields(n),farrayPtr=datar82d,localDE=0, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
            datar8 = datar82d
            if (GFS_control%cpl_imp_mrg) then
              mergeflg(:,:) = datar82d(:,:).eq.missing_value
            endif
            if (mpp_pe() == mpp_root_pe() .and. debug) print *,'in cplIMP,atmos gets ',trim(impfield_name),' datar8=', &
                                                               datar8(isc,jsc), maxval(datar8), minval(datar8)
            found = .true.
          endif

        else if( dimCount == 3) then
          if ( datatype == ESMF_TYPEKIND_R8) then
            call ESMF_FieldGet(importFields(n),farrayPtr=datar83d,localDE=0, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
            found = .true.
          endif
        endif
!
        if (found) then
         if (datar8(isc,jsc) > -99998.0) then
!
        ! get sea land mask: in order to update the coupling fields over the ocean/ice
!        fldname = 'land_mask'
!        if (trim(impfield_name) == trim(fldname)) then
!          findex = queryImportFields(fldname)
!          if (importFieldsValid(findex)) then
!!$omp parallel do default(shared) private(i,j,nb,ix)
!            do j=jsc,jec
!              do i=isc,iec
!                nb = Atm_block%blkno(i,j)
!                ix = Atm_block%ixp(i,j)
!                GFS_data(nb)%Coupling%slimskin_cpl(ix) = datar8(i,j)
!              enddo
!            enddo
!            if( mpp_pe()==mpp_root_pe()) print *,'get land mask from mediator'
!          endif
!        endif


! get sea-state dependent surface roughness (if cplwav2atm=true)
!----------------------------
          fldname = 'wave_z0_roughness_length'
          if (trim(impfield_name) == trim(fldname)) then
            findex = queryImportFields(fldname)
            if (importFieldsValid(findex) .and. GFS_control%cplwav2atm) then
!$omp parallel do default(shared) private(i,j,nb,ix,im,tem)
              do j=jsc,jec
                do i=isc,iec
                  nb = Atm_block%blkno(i,j)
                  ix = Atm_block%ixp(i,j)
                  im = GFS_control%chunk_begin(nb)+ix-1
                  if (GFS_Sfcprop%oceanfrac(im) > zero .and.  datar8(i,j) > zorlmin) then
                    tem = 100.0_GFS_kind_phys * min(0.1_GFS_kind_phys, datar8(i,j))
!                   GFS_Coupling%zorlwav_cpl(im) = tem
                    GFS_Sfcprop%zorlwav(im)      = tem
                    GFS_Sfcprop%zorlw(im)        = tem
                  else
                    GFS_Sfcprop%zorlwav(im) = -999.0_GFS_kind_phys
                  endif
                enddo
              enddo
            endif
          endif

! get sea ice surface temperature
!--------------------------------
          fldname = 'sea_ice_surface_temperature'
          if (trim(impfield_name) == trim(fldname)) then
            findex  = queryImportFields(fldname)
            if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,nb,ix,im)
              do j=jsc,jec
                do i=isc,iec
                  nb = Atm_block%blkno(i,j)
                  ix = Atm_block%ixp(i,j)
                  im = GFS_control%chunk_begin(nb)+ix-1
                  if (GFS_Sfcprop%oceanfrac(im) > zero .and.  datar8(i,j) > 150.0) then
!                   GFS_Coupling%tisfcin_cpl(im) = datar8(i,j)
                    GFS_Sfcprop%tisfc(im)       = datar8(i,j)
                  endif
                enddo
              enddo
            endif
          endif

! get sst:  sst needs to be adjusted by land sea mask before passing to fv3
!--------------------------------------------------------------------------
          fldname = 'sea_surface_temperature'
          if (trim(impfield_name) == trim(fldname)) then
            findex  = queryImportFields(fldname)
            if (importFieldsValid(findex) .and. GFS_control%cplocn2atm) then
!$omp parallel do default(shared) private(i,j,nb,ix,im)
              do j=jsc,jec
                do i=isc,iec
                  nb = Atm_block%blkno(i,j)
                  ix = Atm_block%ixp(i,j)
                  im = GFS_control%chunk_begin(nb)+ix-1
                  if (GFS_Sfcprop%oceanfrac(im) > zero .and. datar8(i,j) > 150.0) then
                    if(mergeflg(i,j)) then
!                     GFS_Coupling%tseain_cpl(im) = GFS_Sfcprop%tsfc(im)
                      GFS_Sfcprop%tsfco(im) = GFS_Sfcprop%tsfc(im)
                      datar8(i,j) = GFS_Sfcprop%tsfc(im)
                    else
!                     GFS_Coupling%tseain_cpl(im) = datar8(i,j)
                      GFS_Sfcprop%tsfco(im)       = datar8(i,j)
                    endif
                  endif
                enddo
              enddo
              if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'get sst from mediator'
            endif
          endif

! get zonal ocean current:
!--------------------------------------------------------------------------
          fldname = 'ocn_current_zonal'
          if (trim(impfield_name) == trim(fldname)) then
            findex  = queryImportFields(fldname)
            if (importFieldsValid(findex) .and. GFS_control%cplocn2atm) then
!$omp parallel do default(shared) private(i,j,nb,ix,im)
              do j=jsc,jec
                do i=isc,iec
                  nb = Atm_block%blkno(i,j)
                  ix = Atm_block%ixp(i,j)
                  im = GFS_control%chunk_begin(nb)+ix-1
                  GFS_Sfcprop%usfco(im) = zero
                  if (GFS_Sfcprop%oceanfrac(im) > zero) then  ! ocean points
                    if(mergeflg(i,j)) then
                      GFS_Sfcprop%usfco(im)       =  zero
                      datar8(i,j) = zero
                    else
                      GFS_Sfcprop%usfco(im)       = datar8(i,j)
                    endif
                  endif
                enddo
              enddo
              if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'get usfco from mediator'
            endif
          endif

! get meridional ocean current:
!--------------------------------------------------------------------------
          fldname = 'ocn_current_merid'
          if (trim(impfield_name) == trim(fldname)) then
            findex  = queryImportFields(fldname)
            if (importFieldsValid(findex) .and. GFS_control%cplocn2atm) then
!$omp parallel do default(shared) private(i,j,nb,ix,im)
              do j=jsc,jec
                do i=isc,iec
                  nb = Atm_block%blkno(i,j)
                  ix = Atm_block%ixp(i,j)
                  im = GFS_control%chunk_begin(nb)+ix-1
                  GFS_Sfcprop%vsfco(im) = zero
                  if (GFS_Sfcprop%oceanfrac(im) > zero) then  ! ocean points
                    if(mergeflg(i,j)) then
                      GFS_Sfcprop%vsfco(im)       =  zero
                      datar8(i,j) = zero
                    else
                      GFS_Sfcprop%vsfco(im)       = datar8(i,j)
                    endif
                  endif
                enddo
              enddo
              if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'get vsfco from mediator'
            endif
          endif

! get sea ice fraction:  fice or sea ice concentration from the mediator
!-----------------------------------------------------------------------
          fldname = 'ice_fraction'
          if (trim(impfield_name) == trim(fldname)) then
            findex  = queryImportFields(fldname)
            if (importFieldsValid(findex)) then
              lcpl_fice = .true.
!$omp parallel do default(shared) private(i,j,nb,ix,im,ofrac)
              do j=jsc,jec
                do i=isc,iec
                  nb = Atm_block%blkno(i,j)
                  ix = Atm_block%ixp(i,j)
                  im = GFS_control%chunk_begin(nb)+ix-1
                  GFS_Coupling%slimskin_cpl(im) = GFS_Sfcprop%slmsk(im)
                  ofrac = GFS_Sfcprop%oceanfrac(im)
                  if (ofrac > zero) then
                    GFS_Sfcprop%fice(im) = max(zero, min(one, datar8(i,j)/ofrac)) !LHS: ice frac wrt water area
                    if (GFS_Sfcprop%fice(im) >= GFS_control%min_seaice) then
                      if (GFS_Sfcprop%fice(im) > one-epsln) GFS_Sfcprop%fice(im) = one
                      if (abs(one-ofrac) < epsln) GFS_Sfcprop%slmsk(im) = 2.0_GFS_kind_phys !slmsk=2 crashes in gcycle on partial land points
!                     GFS_Sfcprop%slmsk(im)         = 2.0_GFS_kind_phys
                      GFS_Coupling%slimskin_cpl(im) = 4.0_GFS_kind_phys
                    else
                      GFS_Sfcprop%fice(im) = zero
                      if (abs(one-ofrac) < epsln) then
                        GFS_Sfcprop%slmsk(im)         = zero
                        GFS_Coupling%slimskin_cpl(im) = zero
                      endif
                    endif
                  endif
                enddo
              enddo
              if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get fice from mediator'
            endif
          endif

! get upward LW flux:  for sea ice covered area
!----------------------------------------------
          fldname = 'lwup_flx_ice'
          if (trim(impfield_name) == trim(fldname)) then
            findex  = queryImportFields(fldname)
            if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,nb,ix,im)
              do j=jsc,jec
!               do i=isc,iec
!                 nb = Atm_block%blkno(i,j)
!                 ix = Atm_block%ixp(i,j)
!                if (GFS_Sfcprop%slmsk(im) < 0.1 .or. GFS_Sfcprop%slmsk(im) > 1.9) then
!                   GFS_Coupling%ulwsfcin_cpl(im) = -datar8(i,j)
!                 endif
!               enddo
                do i=isc,iec
                  nb = Atm_block%blkno(i,j)
                  ix = Atm_block%ixp(i,j)
                  im = GFS_control%chunk_begin(nb)+ix-1
                  if (GFS_Sfcprop%oceanfrac(im) > zero) then
                    GFS_Coupling%ulwsfcin_cpl(im) = -datar8(i,j)
                  endif
                enddo
              enddo
              if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get lwflx from mediator'
            endif
          endif

! get latent heat flux:  for sea ice covered area
!------------------------------------------------
          fldname = 'laten_heat_flx_atm_into_ice'
          if (trim(impfield_name) == trim(fldname)) then
            findex  = queryImportFields(fldname)
            if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,nb,ix,im)
              do j=jsc,jec
                do i=isc,iec
                  nb = Atm_block%blkno(i,j)
                  ix = Atm_block%ixp(i,j)
                  im = GFS_control%chunk_begin(nb)+ix-1
                  if (GFS_Sfcprop%oceanfrac(im) > zero) then
                    GFS_Coupling%dqsfcin_cpl(im) = -datar8(i,j)
                  endif
                enddo
              enddo
              if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get laten_heat from mediator'
            endif
          endif

! get sensible heat flux:  for sea ice covered area
!--------------------------------------------------
          fldname = 'sensi_heat_flx_atm_into_ice'
          if (trim(impfield_name) == trim(fldname)) then
            findex  = queryImportFields(fldname)
            if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,nb,ix,im)
              do j=jsc,jec
                do i=isc,iec
                  nb = Atm_block%blkno(i,j)
                  ix = Atm_block%ixp(i,j)
                  im = GFS_control%chunk_begin(nb)+ix-1
                  if (GFS_Sfcprop%oceanfrac(im) > zero) then
                    GFS_Coupling%dtsfcin_cpl(im) = -datar8(i,j)
                  endif
                enddo
              enddo
              if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get sensi_heat from mediator'
            endif
          endif

! get zonal compt of momentum flux:  for sea ice covered area
!------------------------------------------------------------
          fldname = 'stress_on_air_ice_zonal'
          if (trim(impfield_name) == trim(fldname)) then
            findex  = queryImportFields(fldname)
            if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,nb,ix,im)
              do j=jsc,jec
                do i=isc,iec
                  nb = Atm_block%blkno(i,j)
                  ix = Atm_block%ixp(i,j)
                  im = GFS_control%chunk_begin(nb)+ix-1
                  if (GFS_Sfcprop%oceanfrac(im) > zero) then
                    GFS_Coupling%dusfcin_cpl(im) = -datar8(i,j)
                  endif
                enddo
              enddo
              if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get zonal_moment_flx from mediator'
            endif
          endif

! get meridional compt of momentum flux:  for sea ice covered area
!-----------------------------------------------------------------
          fldname = 'stress_on_air_ice_merid'
          if (trim(impfield_name) == trim(fldname)) then
            findex  = queryImportFields(fldname)
            if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,nb,ix,im)
              do j=jsc,jec
                do i=isc,iec
                  nb = Atm_block%blkno(i,j)
                  ix = Atm_block%ixp(i,j)
                  im = GFS_control%chunk_begin(nb)+ix-1
                  if (GFS_Sfcprop%oceanfrac(im) > zero) then
                    GFS_Coupling%dvsfcin_cpl(im) = -datar8(i,j)
                  endif
                enddo
              enddo
              if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get merid_moment_flx from mediator'
            endif
          endif

! get sea ice volume:  for sea ice covered area
!----------------------------------------------
          fldname = 'sea_ice_volume'
          if (trim(impfield_name) == trim(fldname)) then
            findex  = queryImportFields(fldname)
            if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,nb,ix,im)
              do j=jsc,jec
                do i=isc,iec
                  nb = Atm_block%blkno(i,j)
                  ix = Atm_block%ixp(i,j)
                  im = GFS_control%chunk_begin(nb)+ix-1
                  if (GFS_Sfcprop%oceanfrac(im) > zero) then
!                   GFS_Coupling%hicein_cpl(im) = datar8(i,j)
                    GFS_Sfcprop%hice(im)        = min(datar8(i,j), himax)
                  endif
                enddo
              enddo
              if (mpp_pe() == mpp_root_pe() .and. debug) print *,'fv3 assign_import: get ice_volume from mediator'
            endif
          endif

! get snow volume:  for sea ice covered area
!-------------------------------------------
          fldname = 'snow_volume_on_sea_ice'
          if (trim(impfield_name) == trim(fldname)) then
            findex  = queryImportFields(fldname)
            if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,nb,ix,im)
              do j=jsc,jec
                do i=isc,iec
                  nb = Atm_block%blkno(i,j)
                  ix = Atm_block%ixp(i,j)
                  im = GFS_control%chunk_begin(nb)+ix-1
                  if (GFS_Sfcprop%oceanfrac(im) > zero) then
                    GFS_Coupling%hsnoin_cpl(im) = datar8(i,j)
                  endif
                enddo
              enddo
              if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get snow_volume from mediator'
            endif
          endif

          if (GFS_control%use_cice_alb) then
!
! get instantaneous near IR albedo for diffuse radiation: for sea ice covered area
!---------------------------------------------------------------------------------
            fldname = 'inst_ice_ir_dif_albedo'
            if (trim(impfield_name) == trim(fldname)) then
              findex  = queryImportFields(fldname)
              if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,nb,ix,im)
                do j=jsc,jec
                  do i=isc,iec
                    nb = Atm_block%blkno(i,j)
                    ix = Atm_block%ixp(i,j)
                    im = GFS_control%chunk_begin(nb)+ix-1
                    if (GFS_Sfcprop%oceanfrac(im) > zero) then
!                     GFS_Coupling%sfc_alb_nir_dif_cpl(im) = datar8(i,j)
                      GFS_Sfcprop%albdifnir_ice(im) = datar8(i,j)
                    endif
                  enddo
                enddo
                if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get sfc_alb_nir_dif_cpl from mediator'
              endif
            endif
!
! get instantaneous near IR albedo for direct radiation: for sea ice covered area
!---------------------------------------------------------------------------------
            fldname = 'inst_ice_ir_dir_albedo'
            if (trim(impfield_name) == trim(fldname)) then
              findex  = queryImportFields(fldname)
              if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,nb,ix,im)
                do j=jsc,jec
                  do i=isc,iec
                    nb = Atm_block%blkno(i,j)
                    ix = Atm_block%ixp(i,j)
                    im = GFS_control%chunk_begin(nb)+ix-1
                    if (GFS_Sfcprop%oceanfrac(im) > zero) then
!                     GFS_Coupling%sfc_alb_nir_dir_cpl(im) = datar8(i,j)
                      GFS_Sfcprop%albdirnir_ice(im) = datar8(i,j)
                    endif
                  enddo
                enddo
                if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get sfc_alb_nir_dir_cpl from mediator'
              endif
            endif
!
! get instantaneous visible albedo for diffuse radiation: for sea ice covered area
!---------------------------------------------------------------------------------
            fldname = 'inst_ice_vis_dif_albedo'
            if (trim(impfield_name) == trim(fldname)) then
              findex  = queryImportFields(fldname)
              if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,nb,ix,im)
                do j=jsc,jec
                  do i=isc,iec
                    nb = Atm_block%blkno(i,j)
                    ix = Atm_block%ixp(i,j)
                    im = GFS_control%chunk_begin(nb)+ix-1
                    if (GFS_Sfcprop%oceanfrac(im) > zero) then
!                     GFS_Coupling%sfc_alb_vis_dif_cpl(im) = datar8(i,j)
                      GFS_Sfcprop%albdifvis_ice(im) = datar8(i,j)
                    endif
                  enddo
                enddo
                if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get sfc_alb_vis_dif_cpl from mediator'
              endif
            endif

!
! get instantaneous visible IR albedo for direct radiation: for sea ice covered area
!---------------------------------------------------------------------------------
            fldname = 'inst_ice_vis_dir_albedo'
            if (trim(impfield_name) == trim(fldname)) then
              findex  = queryImportFields(fldname)
              if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,nb,ix,im)
                do j=jsc,jec
                  do i=isc,iec
                    nb = Atm_block%blkno(i,j)
                    ix = Atm_block%ixp(i,j)
                    im = GFS_control%chunk_begin(nb)+ix-1
                    if (GFS_Sfcprop%oceanfrac(im) > zero) then
!                     GFS_Coupling%sfc_alb_vis_dir_cpl(im) = datar8(i,j)
                      GFS_Sfcprop%albdirvis_ice(im) = datar8(i,j)
                    endif
                  enddo
                enddo
                if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get inst_ice_vis_dir_albedo from mediator'
              endif
            endif
          endif

! get upward LW flux:  for open ocean
!----------------------------------------------
          fldname = 'lwup_flx_ocn'
          if (trim(impfield_name) == trim(fldname) .and. GFS_control%use_med_flux) then
            findex  = queryImportFields(fldname)
            if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,nb,ix,im)
              do j=jsc,jec
                do i=isc,iec
                  nb = Atm_block%blkno(i,j)
                  ix = Atm_block%ixp(i,j)
                  im = GFS_control%chunk_begin(nb)+ix-1
                  if (GFS_Sfcprop%oceanfrac(im) > zero) then
                    GFS_Coupling%ulwsfcin_med(im) = -datar8(i,j)
                  endif
                enddo
              enddo
              if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get lwflx for open ocean from mediator'
            endif
          endif

! get latent heat flux:  for open ocean
!------------------------------------------------
          fldname = 'laten_heat_flx_atm_into_ocn'
          if (trim(impfield_name) == trim(fldname) .and. GFS_control%use_med_flux) then
            findex  = queryImportFields(fldname)
            if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,nb,ix,im)
              do j=jsc,jec
                do i=isc,iec
                  nb = Atm_block%blkno(i,j)
                  ix = Atm_block%ixp(i,j)
                  im = GFS_control%chunk_begin(nb)+ix-1
                  if (GFS_Sfcprop%oceanfrac(im) > zero) then
                    GFS_Coupling%dqsfcin_med(im) = -datar8(i,j)
                  endif
                enddo
              enddo
              if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get laten_heat for open ocean from mediator'
            endif
          endif

! get sensible heat flux:  for open ocean
!--------------------------------------------------
          fldname = 'sensi_heat_flx_atm_into_ocn'
          if (trim(impfield_name) == trim(fldname) .and. GFS_control%use_med_flux) then
            findex  = queryImportFields(fldname)
            if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,nb,ix,im)
              do j=jsc,jec
                do i=isc,iec
                  nb = Atm_block%blkno(i,j)
                  ix = Atm_block%ixp(i,j)
                  im = GFS_control%chunk_begin(nb)+ix-1
                  if (GFS_Sfcprop%oceanfrac(im) > zero) then
                    GFS_Coupling%dtsfcin_med(im) = -datar8(i,j)
                  endif
                enddo
              enddo
              if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get sensi_heat for open ocean from mediator'
            endif
          endif

! get zonal compt of momentum flux:  for open ocean
!------------------------------------------------------------
          fldname = 'stress_on_air_ocn_zonal'
          if (trim(impfield_name) == trim(fldname) .and. GFS_control%use_med_flux) then
            findex  = queryImportFields(fldname)
            if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,nb,ix,im)
              do j=jsc,jec
                do i=isc,iec
                  nb = Atm_block%blkno(i,j)
                  ix = Atm_block%ixp(i,j)
                  im = GFS_control%chunk_begin(nb)+ix-1
                  if (GFS_Sfcprop%oceanfrac(im) > zero) then
                    GFS_Coupling%dusfcin_med(im) = -datar8(i,j)
                  endif
                enddo
              enddo
              if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get zonal_moment_flx for open ocean from mediator'
            endif
          endif

! get meridional compt of momentum flux:  for open ocean
!-----------------------------------------------------------------
          fldname = 'stress_on_air_ocn_merid'
          if (trim(impfield_name) == trim(fldname) .and. GFS_control%use_med_flux) then
            findex  = queryImportFields(fldname)
            if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,nb,ix,im)
              do j=jsc,jec
                do i=isc,iec
                  nb = Atm_block%blkno(i,j)
                  ix = Atm_block%ixp(i,j)
                  im = GFS_control%chunk_begin(nb)+ix-1
                  if (GFS_Sfcprop%oceanfrac(im) > zero) then
                    GFS_Coupling%dvsfcin_med(im) = -datar8(i,j)
                  endif
                enddo
              enddo
              if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get merid_moment_flx for open ocean from mediator'
            endif
          endif

! get surface snow area fraction: over land (if cpllnd=true and cpllnd2atm=true)
!------------------------------------------------
          fldname = 'inst_snow_area_fraction_lnd'
          if (trim(impfield_name) == trim(fldname)) then
            findex  = queryImportFields(fldname)
            if (importFieldsValid(findex) .and. GFS_control%cpllnd .and. GFS_control%cpllnd2atm) then
!$omp parallel do default(shared) private(i,j,nb,ix,im)
              do j=jsc,jec
                do i=isc,iec
                  nb = Atm_block%blkno(i,j)
                  ix = Atm_block%ixp(i,j)
                  im = GFS_control%chunk_begin(nb)+ix-1
                  if (GFS_Sfcprop%landfrac(im) > zero) then
                    GFS_Coupling%sncovr1_lnd(im) = datar8(i,j)
                  endif
                enddo
              enddo
              if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get snow area fraction from land'
            endif
          endif

! get latent heat flux: over land (if cpllnd=true and cpllnd2atm=true)
!------------------------------------------------
          fldname = 'inst_laten_heat_flx_lnd'
          if (trim(impfield_name) == trim(fldname)) then
            findex  = queryImportFields(fldname)
            if (importFieldsValid(findex) .and. GFS_control%cpllnd .and. GFS_control%cpllnd2atm) then
!$omp parallel do default(shared) private(i,j,nb,ix,im)
              do j=jsc,jec
                do i=isc,iec
                  nb = Atm_block%blkno(i,j)
                  ix = Atm_block%ixp(i,j)
                  im = GFS_control%chunk_begin(nb)+ix-1
                  if (GFS_Sfcprop%landfrac(im) > zero) then
                    GFS_Coupling%evap_lnd(im) = datar8(i,j)
                  endif
                enddo
              enddo
              if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get latent heat flux from land'
            endif
          endif

! get sensible heat flux: over land (if cpllnd=true and cpllnd2atm=true)
!--------------------------------------------------
          fldname = 'inst_sensi_heat_flx_lnd'
          if (trim(impfield_name) == trim(fldname)) then
            findex  = queryImportFields(fldname)
            if (importFieldsValid(findex) .and. GFS_control%cpllnd .and. GFS_control%cpllnd2atm) then
!$omp parallel do default(shared) private(i,j,nb,ix,im)
              do j=jsc,jec
                do i=isc,iec
                  nb = Atm_block%blkno(i,j)
                  ix = Atm_block%ixp(i,j)
                  im = GFS_control%chunk_begin(nb)+ix-1
                  if (GFS_Sfcprop%landfrac(im) > zero) then
                    GFS_Coupling%hflx_lnd(im) = datar8(i,j)
                  endif
                enddo
              enddo
              if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get sensible heat flux from land'
            endif
          endif

! get surface upward potential latent heat flux: over land (if cpllnd=true and cpllnd2atm=true)
!------------------------------------------------
          fldname = 'inst_potential_laten_heat_flx_lnd'
          if (trim(impfield_name) == trim(fldname)) then
            findex  = queryImportFields(fldname)
            if (importFieldsValid(findex) .and. GFS_control%cpllnd .and. GFS_control%cpllnd2atm) then
!$omp parallel do default(shared) private(i,j,nb,ix,im)
              do j=jsc,jec
                do i=isc,iec
                  nb = Atm_block%blkno(i,j)
                  ix = Atm_block%ixp(i,j)
                  im = GFS_control%chunk_begin(nb)+ix-1
                  if (GFS_Sfcprop%landfrac(im) > zero) then
                    GFS_Coupling%ep_lnd(im) = datar8(i,j)
                  endif
                enddo
              enddo
              if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get potential latent heat flux from land'
            endif
          endif

! get 2m air temperature: over land (if cpllnd=true and cpllnd2atm=true)
!------------------------------------------------
          fldname = 'inst_temp_height2m_lnd'
          if (trim(impfield_name) == trim(fldname)) then
            findex  = queryImportFields(fldname)
            if (importFieldsValid(findex) .and. GFS_control%cpllnd .and. GFS_control%cpllnd2atm) then
!$omp parallel do default(shared) private(i,j,nb,ix,im)
              do j=jsc,jec
                do i=isc,iec
                  nb = Atm_block%blkno(i,j)
                  ix = Atm_block%ixp(i,j)
                  im = GFS_control%chunk_begin(nb)+ix-1
                  if (GFS_Sfcprop%landfrac(im) > zero) then
                    GFS_Coupling%t2mmp_lnd(im) = datar8(i,j)
                  endif
                enddo
              enddo
              if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get temperature at 2m from land'
            endif
          endif

! get 2m specific humidity: over land (if cpllnd=true and cpllnd2atm=true)
!------------------------------------------------
          fldname = 'inst_spec_humid_height2m_lnd'
          if (trim(impfield_name) == trim(fldname)) then
            findex  = queryImportFields(fldname)
            if (importFieldsValid(findex) .and. GFS_control%cpllnd .and. GFS_control%cpllnd2atm) then
!$omp parallel do default(shared) private(i,j,nb,ix,im)
              do j=jsc,jec
                do i=isc,iec
                  nb = Atm_block%blkno(i,j)
                  ix = Atm_block%ixp(i,j)
                  im = GFS_control%chunk_begin(nb)+ix-1
                  if (GFS_Sfcprop%landfrac(im) > zero) then
                    GFS_Coupling%q2mp_lnd(im) = datar8(i,j)
                  endif
                enddo
              enddo
              if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get specific humidity at 2m from land'
            endif
          endif

! get specific humidity: over land (if cpllnd=true and cpllnd2atm=true)
!------------------------------------------------
          fldname = 'inst_spec_humid_lnd'
          if (trim(impfield_name) == trim(fldname)) then
            findex  = queryImportFields(fldname)
            if (importFieldsValid(findex) .and. GFS_control%cpllnd .and. GFS_control%cpllnd2atm) then
!$omp parallel do default(shared) private(i,j,nb,ix,im)
              do j=jsc,jec
                do i=isc,iec
                  nb = Atm_block%blkno(i,j)
                  ix = Atm_block%ixp(i,j)
                  im = GFS_control%chunk_begin(nb)+ix-1
                  if (GFS_Sfcprop%landfrac(im) > zero) then
                    GFS_Coupling%qsurf_lnd(im) = datar8(i,j)
                  endif
                enddo
              enddo
              if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get specific humidity from land'
            endif
          endif

! get upward heat flux in soil (if cpllnd=true and cpllnd2atm=true)
!------------------------------------------------
          fldname = 'inst_upward_heat_flux_lnd'
          if (trim(impfield_name) == trim(fldname)) then
            findex  = queryImportFields(fldname)
            if (importFieldsValid(findex) .and. GFS_control%cpllnd .and. GFS_control%cpllnd2atm) then
!$omp parallel do default(shared) private(i,j,nb,ix,im)
              do j=jsc,jec
                do i=isc,iec
                  nb = Atm_block%blkno(i,j)
                  ix = Atm_block%ixp(i,j)
                  im = GFS_control%chunk_begin(nb)+ix-1
                  if (GFS_Sfcprop%landfrac(im) > zero) then
                    GFS_Coupling%gflux_lnd(im) = datar8(i,j)
                  endif
                enddo
              enddo
              if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get upward heat flux from land'
            endif
          endif

! get surface runoff in soil (if cpllnd=true and cpllnd2atm=true)
!------------------------------------------------
          fldname = 'inst_runoff_rate_lnd'
          if (trim(impfield_name) == trim(fldname)) then
            findex  = queryImportFields(fldname)
            if (importFieldsValid(findex) .and. GFS_control%cpllnd .and. GFS_control%cpllnd2atm) then
!$omp parallel do default(shared) private(i,j,nb,ix,im)
              do j=jsc,jec
                do i=isc,iec
                  nb = Atm_block%blkno(i,j)
                  ix = Atm_block%ixp(i,j)
                  im = GFS_control%chunk_begin(nb)+ix-1
                  if (GFS_Sfcprop%landfrac(im) > zero) then
                    GFS_Coupling%runoff_lnd(im) = datar8(i,j)
                  endif
                enddo
              enddo
              if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get surface runoff from land'
            endif
          endif

! get subsurface runoff in soil (if cpllnd=true and cpllnd2atm=true)
!------------------------------------------------
          fldname = 'inst_subsurface_runoff_rate_lnd'
          if (trim(impfield_name) == trim(fldname)) then
            findex  = queryImportFields(fldname)
            if (importFieldsValid(findex) .and. GFS_control%cpllnd .and. GFS_control%cpllnd2atm) then
!$omp parallel do default(shared) private(i,j,nb,ix,im)
              do j=jsc,jec
                do i=isc,iec
                  nb = Atm_block%blkno(i,j)
                  ix = Atm_block%ixp(i,j)
                  im = GFS_control%chunk_begin(nb)+ix-1
                  if (GFS_Sfcprop%landfrac(im) > zero) then
                    GFS_Coupling%drain_lnd(im) = datar8(i,j)
                  endif
                enddo
              enddo
              if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get subsurface runoff from land'
            endif
          endif

! get momentum exchange coefficient (if cpllnd=true and cpllnd2atm=true)
!------------------------------------------------
          fldname = 'inst_drag_wind_speed_for_momentum'
          if (trim(impfield_name) == trim(fldname)) then
            findex  = queryImportFields(fldname)
            if (importFieldsValid(findex) .and. GFS_control%cpllnd .and. GFS_control%cpllnd2atm) then
!$omp parallel do default(shared) private(i,j,nb,ix,im)
              do j=jsc,jec
                do i=isc,iec
                  nb = Atm_block%blkno(i,j)
                  ix = Atm_block%ixp(i,j)
                  im = GFS_control%chunk_begin(nb)+ix-1
                  if (GFS_Sfcprop%landfrac(im) > zero) then
                    GFS_Coupling%cmm_lnd(im) = datar8(i,j)
                  endif
                enddo
              enddo
              if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get drag wind speed for momentum from land'
            endif
          endif

! get thermal exchange coefficient (if cpllnd=true and cpllnd2atm=true)
!------------------------------------------------
          fldname = 'inst_drag_mass_flux_for_heat_and_moisture'
          if (trim(impfield_name) == trim(fldname)) then
            findex  = queryImportFields(fldname)
            if (importFieldsValid(findex) .and. GFS_control%cpllnd .and. GFS_control%cpllnd2atm) then
!$omp parallel do default(shared) private(i,j,nb,ix,im)
              do j=jsc,jec
                do i=isc,iec
                  nb = Atm_block%blkno(i,j)
                  ix = Atm_block%ixp(i,j)
                  im = GFS_control%chunk_begin(nb)+ix-1
                  if (GFS_Sfcprop%landfrac(im) > zero) then
                    GFS_Coupling%chh_lnd(im) = datar8(i,j)
                  endif
                enddo
              enddo
              if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get thermal exchange coefficient form land'
            endif
          endif

! get function of surface roughness length and green vegetation fraction (if cpllnd=true and cpllnd2atm=true)
!------------------------------------------------
          fldname = 'inst_func_of_roughness_length_and_vfrac'
          if (trim(impfield_name) == trim(fldname)) then
            findex  = queryImportFields(fldname)
            if (importFieldsValid(findex) .and. GFS_control%cpllnd .and. GFS_control%cpllnd2atm) then
!$omp parallel do default(shared) private(i,j,nb,ix,im)
              do j=jsc,jec
                do i=isc,iec
                  nb = Atm_block%blkno(i,j)
                  ix = Atm_block%ixp(i,j)
                  im = GFS_control%chunk_begin(nb)+ix-1
                  if (GFS_Sfcprop%landfrac(im) > zero) then
                    GFS_Coupling%zvfun_lnd(im) = datar8(i,j)
                  endif
                enddo
              enddo
              if (mpp_pe() == mpp_root_pe() .and. debug)  print *,'fv3 assign_import: get func. of roughness length and vfrac form land'
            endif
          endif

        endif ! if (datar8(isc,jsc) > -99999.0) then

!-------------------------------------------------------

       ! For JEDI

        sphum   = get_tracer_index(MODEL_ATMOS, 'sphum')
        liq_wat = get_tracer_index(MODEL_ATMOS, 'liq_wat')
        ice_wat = get_tracer_index(MODEL_ATMOS, 'ice_wat')
        o3mr    = get_tracer_index(MODEL_ATMOS, 'o3mr')

        fldname = 'u'
        if (trim(impfield_name) == trim(fldname)) then
          findex  = queryImportFields(fldname)
          if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,k)
            do k=1,nk
            do j=jsc,jec
              do i=isc,iec
                Atm(mygrid)%u(i,j,k) = datar83d(i-isc+1,j-jsc+1,k)
              enddo
            enddo
            enddo
          endif
        endif

        fldname = 'v'
        if (trim(impfield_name) == trim(fldname)) then
          findex  = queryImportFields(fldname)
          if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,k)
            do k=1,nk
            do j=jsc,jec
              do i=isc,iec
                Atm(mygrid)%v(i,j,k) = datar83d(i-isc+1,j-jsc+1,k)
              enddo
            enddo
            enddo
          endif
        endif

        fldname = 'ua'
        if (trim(impfield_name) == trim(fldname)) then
          findex  = queryImportFields(fldname)
          if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,k)
            do k=1,nk
            do j=jsc,jec
              do i=isc,iec
                Atm(mygrid)%ua(i,j,k) = datar83d(i-isc+1,j-jsc+1,k)
              enddo
            enddo
            enddo
          endif
        endif

        fldname = 'va'
        if (trim(impfield_name) == trim(fldname)) then
          findex  = queryImportFields(fldname)
          if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,k)
            do k=1,nk
            do j=jsc,jec
              do i=isc,iec
                Atm(mygrid)%va(i,j,k) = datar83d(i-isc+1,j-jsc+1,k)
              enddo
            enddo
            enddo
          endif
        endif

        fldname = 't'
        if (trim(impfield_name) == trim(fldname)) then
          findex  = queryImportFields(fldname)
          if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,k)
            do k=1,nk
            do j=jsc,jec
              do i=isc,iec
                Atm(mygrid)%pt(i,j,k) = datar83d(i-isc+1,j-jsc+1,k)
              enddo
            enddo
            enddo
          endif
        endif

        fldname = 'delp'
        if (trim(impfield_name) == trim(fldname)) then
          findex  = queryImportFields(fldname)
          if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,k)
            do k=1,nk
            do j=jsc,jec
              do i=isc,iec
                Atm(mygrid)%delp(i,j,k) = datar83d(i-isc+1,j-jsc+1,k)
              enddo
            enddo
            enddo
          endif
        endif

        fldname = 'sphum'
        if (trim(impfield_name) == trim(fldname) .and. sphum > 0) then
          findex  = queryImportFields(fldname)
          if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,k)
            do k=1,nk
            do j=jsc,jec
              do i=isc,iec
                Atm(mygrid)%q(i,j,k,sphum) = datar83d(i-isc+1,j-jsc+1,k)
              enddo
            enddo
            enddo
          endif
        endif

        fldname = 'ice_wat'
        if (trim(impfield_name) == trim(fldname) .and. ice_wat > 0) then
          findex  = queryImportFields(fldname)
          if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,k)
            do k=1,nk
            do j=jsc,jec
              do i=isc,iec
                Atm(mygrid)%q(i,j,k,ice_wat) = datar83d(i-isc+1,j-jsc+1,k)
              enddo
            enddo
            enddo
          endif
        endif

        fldname = 'liq_wat'
        if (trim(impfield_name) == trim(fldname) .and. liq_wat > 0) then
          findex  = queryImportFields(fldname)
          if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,k)
            do k=1,nk
            do j=jsc,jec
              do i=isc,iec
                Atm(mygrid)%q(i,j,k,sphum) = datar83d(i-isc+1,j-jsc+1,k)
              enddo
            enddo
            enddo
          endif
        endif

        fldname = 'o3mr'
        if (trim(impfield_name) == trim(fldname) .and. o3mr > 0) then
          findex  = queryImportFields(fldname)
          if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,k)
            do k=1,nk
            do j=jsc,jec
              do i=isc,iec
                Atm(mygrid)%q(i,j,k,o3mr) = datar83d(i-isc+1,j-jsc+1,k)
              enddo
            enddo
            enddo
          endif
        endif

        fldname = 'phis'
        if (trim(impfield_name) == trim(fldname)) then
          findex  = queryImportFields(fldname)
          if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j)
            do j=jsc,jec
              do i=isc,iec
                Atm(mygrid)%phis(i,j) = datar82d(i-isc+1,j-jsc+1)
              enddo
            enddo
          endif
        endif

        fldname = 'u_srf'
        if (trim(impfield_name) == trim(fldname)) then
          findex  = queryImportFields(fldname)
          if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j)
            do j=jsc,jec
              do i=isc,iec
                Atm(mygrid)%u_srf(i,j) = datar82d(i-isc+1,j-jsc+1)
              enddo
            enddo
          endif
        endif

        fldname = 'v_srf'
        if (trim(impfield_name) == trim(fldname)) then
          findex  = queryImportFields(fldname)
          if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j)
            do j=jsc,jec
              do i=isc,iec
                Atm(mygrid)%v_srf(i,j) = datar82d(i-isc+1,j-jsc+1)
              enddo
            enddo
          endif
        endif

        ! physics
        fldname = 'slmsk'
        if (trim(impfield_name) == trim(fldname)) then
          findex  = queryImportFields(fldname)
          if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,nb,ix,im)
            do j=jsc,jec
              do i=isc,iec
                nb = Atm_block%blkno(i,j)
                ix = Atm_block%ixp(i,j)
                im = GFS_control%chunk_begin(nb)+ix-1
                GFS_sfcprop%slmsk(im) = datar82d(i-isc+1,j-jsc+1)
              enddo
            enddo
          endif
        endif

        fldname = 'weasd'
        if (trim(impfield_name) == trim(fldname)) then
          findex  = queryImportFields(fldname)
          if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,nb,ix,im)
            do j=jsc,jec
              do i=isc,iec
                nb = Atm_block%blkno(i,j)
                ix = Atm_block%ixp(i,j)
                im = GFS_control%chunk_begin(nb)+ix-1
                GFS_sfcprop%weasd(im) = datar82d(i-isc+1,j-jsc+1)
              enddo
            enddo
          endif
        endif

        fldname = 'tsea'
        if (trim(impfield_name) == trim(fldname)) then
          findex  = queryImportFields(fldname)
          if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,nb,ix,im)
            do j=jsc,jec
              do i=isc,iec
                nb = Atm_block%blkno(i,j)
                ix = Atm_block%ixp(i,j)
                im = GFS_control%chunk_begin(nb)+ix-1
                GFS_sfcprop%tsfco(im) = datar82d(i-isc+1,j-jsc+1)
              enddo
            enddo
          endif
        endif

        fldname = 'vtype'
        if (trim(impfield_name) == trim(fldname)) then
          findex  = queryImportFields(fldname)
          if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,nb,ix,im)
            do j=jsc,jec
              do i=isc,iec
                nb = Atm_block%blkno(i,j)
                ix = Atm_block%ixp(i,j)
                im = GFS_control%chunk_begin(nb)+ix-1
                GFS_sfcprop%vtype(im) = int(datar82d(i-isc+1,j-jsc+1))
              enddo
            enddo
          endif
        endif

        fldname = 'stype'
        if (trim(impfield_name) == trim(fldname)) then
          findex  = queryImportFields(fldname)
          if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,nb,ix,im)
            do j=jsc,jec
              do i=isc,iec
                nb = Atm_block%blkno(i,j)
                ix = Atm_block%ixp(i,j)
                im = GFS_control%chunk_begin(nb)+ix-1
                GFS_sfcprop%stype(im) = int(datar82d(i-isc+1,j-jsc+1))
              enddo
            enddo
          endif
        endif

        fldname = 'vfrac'
        if (trim(impfield_name) == trim(fldname)) then
          findex  = queryImportFields(fldname)
          if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,nb,ix,im)
            do j=jsc,jec
              do i=isc,iec
                nb = Atm_block%blkno(i,j)
                ix = Atm_block%ixp(i,j)
                im = GFS_control%chunk_begin(nb)+ix-1
                GFS_sfcprop%vfrac(im) = datar82d(i-isc+1,j-jsc+1)
              enddo
            enddo
          endif
        endif

        fldname = 'stc'
        if (trim(impfield_name) == trim(fldname)) then
          findex  = queryImportFields(fldname)
          if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,nb,ix,im)
            do j=jsc,jec
              do i=isc,iec
                nb = Atm_block%blkno(i,j)
                ix = Atm_block%ixp(i,j)
                im = GFS_control%chunk_begin(nb)+ix-1
                GFS_sfcprop%stc(im,:) = datar83d(i-isc+1,j-jsc+1,:)
              enddo
            enddo
          endif
        endif

        fldname = 'smc'
        if (trim(impfield_name) == trim(fldname)) then
          findex  = queryImportFields(fldname)
          if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,nb,ix,im)
            do j=jsc,jec
              do i=isc,iec
                nb = Atm_block%blkno(i,j)
                ix = Atm_block%ixp(i,j)
                im = GFS_control%chunk_begin(nb)+ix-1
                GFS_sfcprop%smc(im,:) = datar83d(i-isc+1,j-jsc+1,:)
              enddo
            enddo
          endif
        endif

        fldname = 'snwdph'
        if (trim(impfield_name) == trim(fldname)) then
          findex  = queryImportFields(fldname)
          if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,nb,ix,im)
            do j=jsc,jec
              do i=isc,iec
                nb = Atm_block%blkno(i,j)
                ix = Atm_block%ixp(i,j)
                im = GFS_control%chunk_begin(nb)+ix-1
                GFS_sfcprop%snowd(im) = datar82d(i-isc+1,j-jsc+1)
              enddo
            enddo
          endif
        endif

        fldname = 'f10m'
        if (trim(impfield_name) == trim(fldname)) then
          findex  = queryImportFields(fldname)
          if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,nb,ix,im)
            do j=jsc,jec
              do i=isc,iec
                nb = Atm_block%blkno(i,j)
                ix = Atm_block%ixp(i,j)
                im = GFS_control%chunk_begin(nb)+ix-1
                GFS_sfcprop%f10m(im) = datar82d(i-isc+1,j-jsc+1)
              enddo
            enddo
          endif
        endif

        fldname = 'zorl'
        if (trim(impfield_name) == trim(fldname)) then
          findex  = queryImportFields(fldname)
          if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,nb,ix,im)
            do j=jsc,jec
              do i=isc,iec
                nb = Atm_block%blkno(i,j)
                ix = Atm_block%ixp(i,j)
                im = GFS_control%chunk_begin(nb)+ix-1
                GFS_sfcprop%zorl(im) = datar82d(i-isc+1,j-jsc+1)
              enddo
            enddo
          endif
        endif

        fldname = 't2m'
        if (trim(impfield_name) == trim(fldname)) then
          findex  = queryImportFields(fldname)
          if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,nb,ix,im)
            do j=jsc,jec
              do i=isc,iec
                nb = Atm_block%blkno(i,j)
                ix = Atm_block%ixp(i,j)
                im = GFS_control%chunk_begin(nb)+ix-1
                GFS_sfcprop%t2m(im) = datar82d(i-isc+1,j-jsc+1)
              enddo
            enddo
          endif
        endif

        fldname = 'hflx_fire'
        if (trim(impfield_name) == trim(fldname)) then
          findex  = queryImportFields(fldname)
          if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,nb,ix)
            do j=jsc,jec
              do i=isc,iec
                nb = Atm_block%blkno(i,j)
                ix = Atm_block%ixp(i,j)
                im = GFS_control%chunk_begin(nb)+ix-1
                GFS_sfcprop%hflx_fire(im) = datar82d(i-isc+1,j-jsc+1)
              enddo
            enddo
          endif
        endif

        fldname = 'evap_fire'
        if (trim(impfield_name) == trim(fldname)) then
          findex  = queryImportFields(fldname)
          if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,nb,ix)
            do j=jsc,jec
              do i=isc,iec
                nb = Atm_block%blkno(i,j)
                ix = Atm_block%ixp(i,j)
                im = GFS_control%chunk_begin(nb)+ix-1
                GFS_sfcprop%evap_fire(im) = datar82d(i-isc+1,j-jsc+1)
              enddo
            enddo
          endif
        endif

        fldname = 'smoke_fire'
        if (trim(impfield_name) == trim(fldname)) then
          findex  = queryImportFields(fldname)
          if (importFieldsValid(findex)) then
!$omp parallel do default(shared) private(i,j,nb,ix)
            do j=jsc,jec
              do i=isc,iec
                nb = Atm_block%blkno(i,j)
                ix = Atm_block%ixp(i,j)
                im = GFS_control%chunk_begin(nb)+ix-1
                GFS_sfcprop%smoke_fire(im) = datar82d(i-isc+1,j-jsc+1)
              enddo
            enddo
          endif
        endif

          ! write post merge import data to NetCDF file.
          if (GFS_control%cpl_imp_dbg) then
            call ESMF_FieldGet(importFields(n), grid=grid, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

            dbgField = ESMF_FieldCreate(grid=grid, farrayPtr=datar8, name=impfield_name, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

            write (currtimestring, "(I4.4,'-',I2.2,'-',I2.2,'T',I2.2,':',I2.2,':',I2.2)") &
                                   jdat(1), jdat(2), jdat(3), jdat(5), jdat(6), jdat(7)
            call ESMF_FieldWrite(dbgField, fileName='fv3_merge_'//trim(impfield_name)//'_'// &
                                 trim(currtimestring)//'.nc', rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

            call ESMF_FieldDestroy(dbgField, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
          endif

        endif ! if (found) then
      endif   ! if (isFieldCreated) then
    enddo
!
    deallocate(mergeflg)
    deallocate(datar8)

! update sea ice related fields:
    if( lcpl_fice ) then
!$omp parallel do default(shared) private(i,j,nb,ix,tem,im)
      do j=jsc,jec
        do i=isc,iec
          nb = Atm_block%blkno(i,j)
          ix = Atm_block%ixp(i,j)
          im = GFS_control%chunk_begin(nb)+ix-1
          if (GFS_Sfcprop%oceanfrac(im) > zero) then
            if (GFS_Sfcprop%fice(im) >= GFS_control%min_seaice) then

              GFS_Coupling%hsnoin_cpl(im) = min(hsmax, GFS_Coupling%hsnoin_cpl(im) &
                                                              / GFS_Sfcprop%fice(im))
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
!             GFS_Sfcprop%albdirvis_ice(im) = -9999.0  !                 ,,
!             GFS_Sfcprop%albdirnir_ice(im) = -9999.0  !                 ,,
!             GFS_Sfcprop%albdifvis_ice(im) = -9999.0  !                 ,,
!             GFS_Sfcprop%albdifnir_ice(im) = -9999.0  !                 ,,
              if (abs(one-GFS_Sfcprop%oceanfrac(im)) < epsln) then !  100% open water
                GFS_Coupling%slimskin_cpl(im) = zero
                GFS_Sfcprop%slmsk(im)         = zero
              endif
            endif
          endif
        enddo
      enddo
    endif

    rc=0
!
  end subroutine assign_importdata

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
    integer                :: sphum, liq_wat, ice_wat, o3mr
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
            !--- JEDI fields
            case ('u')
              call block_atmos_copy(datar83d, Atm(mygrid)%u(isc:iec,jsc:jec,:), Atm_block, nb, offset=1, rc=localrc)
            case ('v')
              call block_atmos_copy(datar83d, Atm(mygrid)%v(isc:iec,jsc:jec,:), Atm_block, nb, offset=1, rc=localrc)
            case ('ua')
              call block_atmos_copy(datar83d, Atm(mygrid)%ua(isc:iec,jsc:jec,:),Atm_block, nb, offset=1, rc=localrc)
            case ('va')
              call block_atmos_copy(datar83d, Atm(mygrid)%va(isc:iec,jsc:jec,:), Atm_block, nb, offset=1, rc=localrc)
            case ('t')
              call block_atmos_copy(datar83d, Atm(mygrid)%pt(isc:iec,jsc:jec,:), Atm_block, nb, offset=1, rc=localrc)
            case ('delp')
              call block_atmos_copy(datar83d, Atm(mygrid)%delp(isc:iec,jsc:jec,:), Atm_block, nb, offset=1, rc=localrc)
            case ('sphum')
              sphum = get_tracer_index(MODEL_ATMOS, 'sphum')
              call block_atmos_copy(datar83d, Atm(mygrid)%q(isc:iec,jsc:jec,:,:), sphum, Atm_block, nb, offset=1, rc=localrc)
            case ('ice_wat')
              ice_wat = get_tracer_index(MODEL_ATMOS, 'ice_wat')
              call block_atmos_copy(datar83d, Atm(mygrid)%q(isc:iec,jsc:jec,:,:), ice_wat, Atm_block, nb, offset=1, rc=localrc)
            case ('liq_wat')
              liq_wat = get_tracer_index(MODEL_ATMOS, 'liq_wat')
              call block_atmos_copy(datar83d, Atm(mygrid)%q(isc:iec,jsc:jec,:,:), liq_wat, Atm_block, nb, offset=1, rc=localrc)
            case ('o3mr')
              o3mr = get_tracer_index(MODEL_ATMOS, 'o3mr')
              call block_atmos_copy(datar83d, Atm(mygrid)%q(isc:iec,jsc:jec,:,:), o3mr, Atm_block, nb, offset=1, rc=localrc)
            case ('phis')
              call block_atmos_copy(datar82d, Atm(mygrid)%phis(isc:iec,jsc:jec), Atm_block, nb, offset=1, rc=localrc)
            case ('u_srf')
              call block_atmos_copy(datar82d, Atm(mygrid)%u_srf(isc:iec,jsc:jec), Atm_block, nb, offset=1, rc=localrc)
            case ('v_srf')
              call block_atmos_copy(datar82d, Atm(mygrid)%v_srf(isc:iec,jsc:jec), Atm_block, nb, offset=1, rc=localrc)
            case ('weasd')
              call block_data_copy(datar82d, GFS_Sfcprop%weasd, Atm_block, nb, offset=GFS_Control%chunk_begin(nb), rc=localrc)
            case ('tsea')
              call block_data_copy(datar82d, GFS_Sfcprop%tsfco, Atm_block, nb, offset=GFS_Control%chunk_begin(nb), rc=localrc)
            case ('vtype')
              call block_data_copy(datar82d, GFS_Sfcprop%vtype, Atm_block, nb, offset=GFS_Control%chunk_begin(nb), rc=localrc)
            case ('stype')
              call block_data_copy(datar82d, GFS_Sfcprop%stype, Atm_block, nb, offset=GFS_Control%chunk_begin(nb), rc=localrc)
            case ('vfrac')
              call block_data_copy(datar82d, GFS_Sfcprop%vfrac, Atm_block, nb, offset=GFS_Control%chunk_begin(nb), rc=localrc)
            case ('stc')
              call block_data_copy(datar83d, GFS_Sfcprop%stc, Atm_block, nb, offset=GFS_Control%chunk_begin(nb), rc=localrc)
            case ('smc')
              call block_data_copy(datar83d, GFS_Sfcprop%smc, Atm_block, nb, offset=GFS_Control%chunk_begin(nb), rc=localrc)
            case ('snwdph')
              call block_data_copy(datar82d, GFS_Sfcprop%snowd, Atm_block, nb, offset=GFS_Control%chunk_begin(nb), rc=localrc)
            case ('f10m')
              call block_data_copy(datar82d, GFS_Sfcprop%f10m, Atm_block, nb, offset=GFS_Control%chunk_begin(nb), rc=localrc)
            case ('zorl')
              call block_data_copy(datar82d, GFS_Sfcprop%zorl, Atm_block, nb, offset=GFS_Control%chunk_begin(nb), rc=localrc)
            case ('t2m')
              call block_data_copy(datar82d, GFS_Sfcprop%t2m, Atm_block, nb, offset=GFS_Control%chunk_begin(nb), rc=localrc)
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
!------------------------------------------------------------------------------
  subroutine atmos_model_get_nth_domain_info(n, layout, nx, ny, pelist)
   integer, intent(in)  :: n
   integer, intent(out) :: layout(2)
   integer, intent(out) :: nx, ny
   integer, pointer, intent(out) :: pelist(:)

   call get_nth_domain_info(n, layout, nx, ny, pelist)

  end subroutine atmos_model_get_nth_domain_info

end module atmos_model_mod
