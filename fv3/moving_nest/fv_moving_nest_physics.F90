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

!***********************************************************************
!> @file
!! @brief Provides Moving Nest functionality for physics and surface variables
!! @author W. Ramstrom.  Collaboration with Bin Liu and Chunxi Zhang, EMC
!! @email William.Ramstrom@noaa.gov
! =======================================================================!


! =======================================================================!
!
! Notes
!
!------------------------------------------------------------------------
! Moving Nest Subroutine Naming Convention
!-----------------------------------------------------------------------
!
! mn_meta_* subroutines perform moving nest operations for FV3 metadata.
!               These routines will run only once per nest move.
!
! mn_var_*  subroutines perform moving nest operations for an individual FV3 variable.
!               These routines will run many times per nest move.
!
! mn_prog_* subroutines perform moving nest operations for the list of prognostic fields.
!               These routines will run only once per nest move.
!
! mn_phys_* subroutines perform moving nest operations for the list of physics fields.
!               These routines will run only once per nest move.
!
! =======================================================================!

module fv_moving_nest_physics_mod

  use block_control_mod,      only: block_control_type
  use mpp_mod,                only: mpp_pe, mpp_sync, mpp_sync_self, mpp_send, mpp_error, NOTE, FATAL
  use mpp_domains_mod,        only: mpp_update_domains, mpp_get_data_domain, mpp_get_global_domain
  use mpp_domains_mod,        only: mpp_define_nest_domains, mpp_shift_nest_domains, nest_domain_type, domain2d
  use mpp_domains_mod,        only: mpp_get_C2F_index, mpp_update_nest_fine
  use mpp_domains_mod,        only: mpp_get_F2C_index, mpp_update_nest_coarse
  use mpp_domains_mod,        only: NORTH, SOUTH, EAST, WEST, CORNER, CENTER
  use mpp_domains_mod,        only: NUPDATE, SUPDATE, EUPDATE, WUPDATE, DGRID_NE

  use GFS_typedefs,           only: GFS_sfcprop_type, GFS_tbd_type, GFS_cldprop_type, &
                                    GFS_grid_type, GFS_diag_type, GFS_control_type, kind_phys
  use GFS_init,               only: GFS_grid_populate

  use boundary_mod,           only: update_coarse_grid, update_coarse_grid_mpp
#ifdef OVERLOAD_R4
  use constantsR4_mod,        only: cp_air, rdgas, grav, rvgas, kappa, pstd_mks, hlv
#else
  use constants_mod,          only: cp_air, rdgas, grav, rvgas, kappa, pstd_mks, hlv
#endif
  use field_manager_mod,      only: MODEL_ATMOS
  use fv_arrays_mod,          only: fv_atmos_type, fv_nest_type, fv_grid_type, R_GRID
  use fv_moving_nest_types_mod,   only: fv_moving_nest_prog_type, fv_moving_nest_physics_type, mn_surface_grids, fv_moving_nest_type
  use fv_arrays_mod,          only: allocate_fv_nest_bc_type, deallocate_fv_nest_bc_type
  use fv_grid_tools_mod,      only: init_grid
  use fv_grid_utils_mod,      only: grid_utils_init, ptop_min, dist2side_latlon
  use fv_mapz_mod,            only: Lagrangian_to_Eulerian, moist_cv, compute_total_energy
  use fv_nesting_mod,         only: dealloc_nested_buffers
  use fv_nwp_nudge_mod,       only: do_adiabatic_init
  use init_hydro_mod,         only: p_var
  use tracer_manager_mod,     only: get_tracer_index, get_tracer_names
  use fv_moving_nest_utils_mod,  only: alloc_halo_buffer, grid_geometry, output_grid_to_nc
  use fv_moving_nest_utils_mod,  only: fill_nest_from_buffer, fill_nest_from_buffer_cell_center, fill_nest_from_buffer_nearest_neighbor
  use fv_moving_nest_utils_mod,  only: fill_nest_halos_from_parent, fill_grid_from_supergrid, fill_weight_grid
  use fv_moving_nest_utils_mod,  only: alloc_read_data
  use fv_moving_nest_utils_mod,  only: fill_nest_from_buffer_cell_center_masked
  use fv_moving_nest_utils_mod,  only: fill_nest_halos_from_parent_masked

  use fv_moving_nest_mod,     only: mn_var_fill_intern_nest_halos, mn_var_dump_to_netcdf, mn_var_shift_data, calc_nest_alignment
  use fv_moving_nest_types_mod, only: Moving_nest
  implicit none

#ifdef NO_QUAD_PRECISION
  ! 64-bit precision (kind=8)
  integer, parameter:: f_p = selected_real_kind(15)
#else
  ! Higher precision (kind=16) for grid geometrical factors:
  integer, parameter:: f_p = selected_real_kind(20)
#endif

#ifdef OVERLOAD_R4
  real, parameter:: real_snan=x'FFBFFFFF'
#else
  real, parameter:: real_snan=x'FFF7FFFFFFFFFFFF'
#endif

  logical :: debug_log = .false.
  logical :: move_physics = .true.       ! Always true, unless developer sets move_physics to .False. here for debugging.
  logical :: move_nsst = .true.          ! Value is reset in fv_moving_nest_main.F90 from namelist options

  !! Persistent variables to enable debug printing after range warnings.
  !type (fv_atmos_type), pointer                 :: save_Atm_n
  !type (block_control_type), pointer            :: save_Atm_block
  !type(GFS_control_type), pointer               :: save_GFS_control
  !type(GFS_sfcprop_type), pointer               :: save_GFS_sfcprop

#include <fms_platform.h>

contains


  subroutine mn_phys_set_slmsk(Atm, n, mn_static, ioffset, joffset, refine)
    type(fv_atmos_type), intent(inout),allocatable   :: Atm(:)              !< Array of atmospheric data
    integer, intent(in)                              :: n                   !< Current grid number
    type(mn_surface_grids), intent(in)               :: mn_static           !< Static surface data
    integer, intent(in)                              :: ioffset, joffset    !< Current nest offset in i,j direction
    integer, intent(in)                              :: refine              !< Nest refinement ratio

    integer                 :: i_pe, j_pe, i_idx, j_idx

    !print '("[INFO] MASK inside mn_phys_set_slmsk npe=",I0)', mpp_pe()
    ! Setup local land sea mask grid for masked interpolations
    do i_pe = Atm(n)%bd%isd, Atm(n)%bd%ied
      do j_pe = Atm(n)%bd%jsd, Atm(n)%bd%jed
        i_idx = (ioffset-1)*refine + i_pe
        j_idx = (joffset-1)*refine + j_pe

        Moving_nest(n)%mn_phys%slmsk(i_pe, j_pe) = mn_static%fp_ls%ls_mask_grid(i_idx, j_idx)
      enddo
    enddo
  end subroutine mn_phys_set_slmsk


  !>@brief The subroutine 'mn_phys_reset_sfc_props' sets the static surface parameters from the high-resolution input file data
  !>@details This subroutine relies on earlier code reading the data from files into the mn_static data structure
  !!  This subroutine does not yet handle ice points or frac_grid - fractional landfrac/oceanfrac values
  subroutine mn_phys_reset_sfc_props(Atm, n, mn_static, Atm_block, GFS_Sfcprop, ioffset, joffset, refine)
    type(fv_atmos_type), intent(inout),allocatable   :: Atm(:)              !< Array of atmospheric data
    integer, intent(in)                              :: n                   !< Current grid number
    type(mn_surface_grids), intent(in)               :: mn_static           !< Static surface data
    type(block_control_type), intent(in)             :: Atm_block           !< Physics block layout
    type(GFS_sfcprop_type), intent(inout)            :: GFS_Sfcprop         !< Physics variable data
    integer, intent(in)                              :: ioffset, joffset    !< Current nest offset in i,j direction
    integer, intent(in)                              :: refine              !< Nest refinement ratio

    ! For iterating through physics/surface vector data
    integer                 :: nb, blen, ix, i_pe, j_pe, i_idx, j_idx, im
    real(kind=kind_phys)    :: phys_oro

    !print '("[INFO] MASK inside mn_phys_reset_sfc_props npe=",I0)', mpp_pe()
    call mn_phys_set_slmsk(Atm, n, mn_static, ioffset, joffset, refine)

    !  Reset the variables from the fix_sfc files
    im = 0
    do nb = 1,Atm_block%nblks
      blen = Atm_block%blksz(nb)
      do ix = 1, blen
        i_pe = Atm_block%index(nb)%ii(ix)
        j_pe = Atm_block%index(nb)%jj(ix)

        i_idx = (ioffset-1)*refine + i_pe
        j_idx = (joffset-1)*refine + j_pe

        im = im + 1

        ! Reset the land sea mask from the hires parent data
        GFS_Sfcprop%slmsk(im) = mn_static%fp_ls%ls_mask_grid(i_idx, j_idx)

        !  IFD values are 0 for land, and 1 for oceans/lakes -- reverse of the land sea mask
        !  Land Sea Mask has values of 0 for oceans/lakes, 1 for land, 2 for sea ice
        !  TODO figure out what ifd should be for sea ice
        if (mn_static%fp_ls%ls_mask_grid(i_idx, j_idx) .eq. 1 ) then
          if (move_nsst) GFS_Sfcprop%ifd(im) = 0         ! Land
          GFS_Sfcprop%oceanfrac(im) = 0   ! Land -- TODO permit fractions
          GFS_Sfcprop%landfrac(im) = 1    ! Land -- TODO permit fractions
        else
          if (move_nsst) GFS_Sfcprop%ifd(im) = 1         ! Ocean
          GFS_Sfcprop%oceanfrac(im) = 1   ! Ocean -- TODO permit fractions
          GFS_Sfcprop%landfrac(im) = 0    ! Ocean -- TODO permit fractions
        endif

        GFS_Sfcprop%tg3(im) = mn_static%fp_fix%deep_soil_temp_grid(i_idx, j_idx)

        ! Follow logic from FV3/io/FV3GFS_io.F90 line 1187
        ! TODO this will need to be more complicated if we support frac_grid
        !if (nint(mn_static%soil_type_grid(i_idx, j_idx)) == 14 .or. int(mn_static%soil_type_grid(i_idx, j_idx)+0.5) <= 0) then
        !if (nint(mn_static%soil_type_grid(i_idx, j_idx)) == 14 .or.

        !if ( (mn_static%ls_mask_grid(i_idx, j_idx) .eq. 1 .and. nint(mn_static%land_frac_grid(i_idx, j_idx)) == 0) .or. &
        !    mn_static%soil_type_grid(i_idx, j_idx) < 0.5) then
        if (mn_static%fp_ls%ls_mask_grid(i_idx, j_idx) .eq. 1 .and. nint(mn_static%fp_ls%land_frac_grid(i_idx, j_idx)) == 0 ) then
          ! Water soil type == lake, etc. -- override the other variables and make this water

          if (move_nsst) GFS_Sfcprop%ifd(im) = 1         ! Ocean
          GFS_Sfcprop%oceanfrac(im) = 1   ! Ocean -- TODO permit fractions
          GFS_Sfcprop%landfrac(im) = 0    ! Ocean -- TODO permit fractions

          GFS_Sfcprop%stype(im) = 0
          GFS_Sfcprop%slmsk(im) = 0
        else
          GFS_Sfcprop%stype(im) = nint(mn_static%fp_ls%soil_type_grid(i_idx, j_idx))
        endif

        !GFS_Sfcprop%vfrac(im) = mn_static%veg_frac_grid(i_idx, j_idx)
        GFS_Sfcprop%vtype(im) = nint(mn_static%fp_fix%veg_type_grid(i_idx, j_idx))
        GFS_Sfcprop%slope(im) = nint(mn_static%fp_fix%slope_type_grid(i_idx, j_idx))
        GFS_Sfcprop%snoalb(im) = mn_static%fp_fix%max_snow_alb_grid(i_idx, j_idx)

        GFS_Sfcprop%facsf(im) = mn_static%fp_fix%facsf_grid(i_idx, j_idx)
        GFS_Sfcprop%facwf(im) = mn_static%fp_fix%facwf_grid(i_idx, j_idx)

        GFS_Sfcprop%alvsf(im) = mn_static%fp_fix%alvsf_grid(i_idx, j_idx)
        GFS_Sfcprop%alvwf(im) = mn_static%fp_fix%alvwf_grid(i_idx, j_idx)
        GFS_Sfcprop%alnsf(im) = mn_static%fp_fix%alnsf_grid(i_idx, j_idx)
        GFS_Sfcprop%alnwf(im) = mn_static%fp_fix%alnwf_grid(i_idx, j_idx)

        ! Reset the orography in the physics arrays, using the smoothed values from above
        phys_oro =  Atm(n)%phis(i_pe, j_pe) / grav
        GFS_Sfcprop%oro(im) = phys_oro
        GFS_Sfcprop%oro_uf(im) = phys_oro

      enddo
    enddo

  end subroutine mn_phys_reset_sfc_props

  !>@brief The subroutine 'mn_phys_reset_phys_latlon' sets the lat/lons from the high-resolution input file data
  !>@details This subroutine sets lat/lons of the moved nest, then recalculates all the derived quantities (dx,dy,etc.)
  subroutine mn_reset_phys_latlon(Atm, n, tile_geo, fp_super_tile_geo, Atm_block, GFS_control, GFS_grid)
    type(fv_atmos_type), allocatable, intent(in)      :: Atm(:)               !< Array of atmospheric data
    integer, intent(in)                  :: n                    !< Current grid number
    type(grid_geometry), intent(in)      :: tile_geo             !< Bounds of this grid
    type(grid_geometry), intent(in)      :: fp_super_tile_geo    !< Bounds of high-resolution parent grid
    type(block_control_type), intent(in) :: Atm_block            !< Physics block layout
    type(GFS_control_type), intent(in)   :: GFS_control          !< Physics metadata
    type(GFS_grid_type), intent(inout)   :: GFS_grid             !< Physics variable data

    integer :: isc, jsc, iec, jec
    integer :: x, y, fp_i, fp_j
    integer :: nest_x, nest_y, parent_x, parent_y
    integer :: this_pe

    real(kind=kind_phys), allocatable :: lats(:,:), lons(:,:), area(:,:)

    this_pe = mpp_pe()

    isc = Atm(n)%bd%isc
    jsc = Atm(n)%bd%jsc
    iec = Atm(n)%bd%iec
    jec = Atm(n)%bd%jec

    allocate(lats(isc:iec, jsc:jec))
    allocate(lons(isc:iec, jsc:jec))
    allocate(area(isc:iec, jsc:jec))

    call calc_nest_alignment(Atm, n, nest_x, nest_y, parent_x, parent_y)

    do x = isc, iec
      do y = jsc, jec
        fp_i = (x - nest_x) * 2 + parent_x
        fp_j = (y - nest_y) * 2 + parent_y

        lons(x,y) = fp_super_tile_geo%lons(fp_i, fp_j)
        lats(x,y) = fp_super_tile_geo%lats(fp_i, fp_j)

        ! Need to add the areas from 4 squares, because the netCDF file has areas calculated for the supergrid cells
        !  We need the area of the whole center of the cell.
        !  Example dimensions for C288_grid.tile6.nc
        !   longitude -- x(577,577)
        !   latitude  -- y(577,577)
        !   area      -- x(576,576)

        !  Extracting lat/lon/area from Supergrid
        !
        !   1,1----2,1----3,1
        !    |      |      |
        !    | a1,1 | a2,1 |
        !    |      |      |
        !   1,2----2,2----3,2
        !    |      |      |
        !    | a1,2 | a2,2 |
        !    |      |      |
        !   1,3----2,3----3,3
        !
        !  The model A-grid cell 1,1 is centered at supergrid location 2,2
        !    The area of the A-grid cell is the sum of the 4 supergrid areas   A = a(1,1) + a(1,2) + a(2,1) + a(2,2)

        area(x,y) = fp_super_tile_geo%area(fp_i - 1, fp_j - 1) + fp_super_tile_geo%area(fp_i - 1, fp_j) + &
            fp_super_tile_geo%area(fp_i, fp_j - 1) + fp_super_tile_geo%area(fp_i, fp_j)   ! TODO make sure these offsets are correct.
      enddo
    enddo

    call GFS_grid_populate(GFS_Grid, lons, lats, area)

    deallocate(lats)
    deallocate(lons)
    deallocate(area)

  end subroutine mn_reset_phys_latlon

  !>@brief The subroutine 'mn_phys_fill_temp_variables' extracts 1D physics data into a 2D array for nest motion
  !>@details This subroutine fills in the mn_phys structure on the Atm object with 2D arrays of physics/surface variables.
  !!  Note that ice variables are not yet handled.
  subroutine mn_phys_fill_temp_variables(Atm, Atm_block, GFS_control, GFS_sfcprop, GFS_tbd, GFS_cldprop, GFS_intdiag, n, child_grid_num, is_fine_pe, npz)
    type(fv_atmos_type), allocatable, target, intent(inout)  :: Atm(:)            !< Array of atmospheric data
    type (block_control_type), target, intent(in)            :: Atm_block         !< Physics block layout
    type(GFS_control_type), target, intent(in)               :: GFS_control       !< Physics metadata
    type(GFS_sfcprop_type), target, intent(in)               :: GFS_sfcprop       !< Physics variable data (surface)
    type(GFS_tbd_type), target, intent(in)                   :: GFS_tbd           !< Physics variable data (tbd)
    type(GFS_cldprop_type), target, intent(in)               :: GFS_cldprop       !< Physics variable data (clouds)
    type(GFS_diag_type), target, intent(in)                  :: GFS_intdiag       !< Physics variable data (clouds)
    integer, intent(in)                                      :: n, child_grid_num !< Current grid number, child grid number
    logical, intent(in)                                      :: is_fine_pe        !< Is this a nest PE?
    integer, intent(in)                                      :: npz               !< Number of vertical levels

    integer :: isd, ied, jsd, jed
    integer :: is, ie, js, je
    integer :: this_pe

    integer :: nb, blen, i, j, k, ix, nv, im
    type(fv_moving_nest_physics_type), pointer       :: mn_phys
    integer :: err_field = 0

    this_pe = mpp_pe()

    !save_Atm_n => Atm(n)
    !save_Atm_block => Atm_block
    !save_GFS_control => GFS_control
    !save_GFS_sfcprop => GFS_sfcprop

    isd = Atm(n)%bd%isd
    ied = Atm(n)%bd%ied
    jsd = Atm(n)%bd%jsd
    jed = Atm(n)%bd%jed

    !if (is_fine_pe) call dump_surface_physics(isd+8, jsd+8, npz-1)

    is = Atm(n)%bd%is
    ie = Atm(n)%bd%ie
    js = Atm(n)%bd%js
    je = Atm(n)%bd%je

    mn_phys => Moving_nest(n)%mn_phys

    mn_phys%ts(is:ie, js:je) =  Atm(n)%ts(is:ie, js:je)

    im = 0
    do nb = 1,Atm_block%nblks
      blen = Atm_block%blksz(nb)
      do ix = 1, blen
        ! Get the indices only once, before iterating through vertical levels or number of variables
        i = Atm_block%index(nb)%ii(ix)
        j = Atm_block%index(nb)%jj(ix)
        im = im + 1

        if (move_physics) then
          do k = 1, GFS_control%lsoil
            mn_phys%smc(i,j,k) = GFS_sfcprop%smc(im,k)
            mn_phys%stc(i,j,k) = GFS_sfcprop%stc(im,k)
            mn_phys%slc(i,j,k) = GFS_sfcprop%slc(im,k)
          enddo

          mn_phys%emis_lnd(i,j)      = GFS_sfcprop%emis_lnd(im)
          mn_phys%emis_ice(i,j)      = GFS_sfcprop%emis_ice(im)
          mn_phys%emis_wat(i,j)      = GFS_sfcprop%emis_wat(im)

          !mn_phys%sfalb_lnd(i,j)     = GFS_sfcprop%sfalb_lnd(im)
          !mn_phys%sfalb_lnd_bck(i,j) = GFS_sfcprop%sfalb_lnd_bck(im)
          !mn_phys%semis(i,j)      = GFS_Radtend%semis(im)
          !mn_phys%semisbase(i,j)      = GFS_sfcprop%semisbase(im)
          !mn_phys%sfalb(i,j)      = GFS_Radtend%sfalb(im)

          mn_phys%albdirvis_lnd(i,j) = GFS_sfcprop%albdirvis_lnd(im)
          mn_phys%albdirnir_lnd(i,j) = GFS_sfcprop%albdirnir_lnd(im)
          mn_phys%albdifvis_lnd(i,j) = GFS_sfcprop%albdifvis_lnd(im)
          mn_phys%albdifnir_lnd(i,j) = GFS_sfcprop%albdifnir_lnd(im)

          mn_phys%u10m(i,j)  = GFS_intdiag%u10m(im)
          mn_phys%v10m(i,j)  = GFS_intdiag%v10m(im)
          mn_phys%tprcp(i,j) = GFS_sfcprop%tprcp(im)

          do k = 1, GFS_control%nmtvr
            mn_phys%hprime(i,j,k)  = GFS_sfcprop%hprime(im,k)
          enddo

          mn_phys%lakefrac(i,j) = GFS_Sfcprop%lakefrac(im)
          mn_phys%lakedepth(i,j) = GFS_Sfcprop%lakedepth(im)

          mn_phys%canopy(i,j) = GFS_Sfcprop%canopy(im)
          mn_phys%vegfrac(i,j)= GFS_Sfcprop%vfrac(im)
          mn_phys%uustar(i,j) = GFS_Sfcprop%uustar(im)
          mn_phys%shdmin(i,j) = GFS_Sfcprop%shdmin(im)
          mn_phys%shdmax(i,j) = GFS_Sfcprop%shdmax(im)
          mn_phys%zorl(i,j)   = GFS_Sfcprop%zorl(im)
          mn_phys%zorll(i,j)  = GFS_Sfcprop%zorll(im)
          mn_phys%zorlwav(i,j)= GFS_Sfcprop%zorlwav(im)
          mn_phys%zorlw(i,j)  = GFS_Sfcprop%zorlw(im)
          mn_phys%usfco(i,j)  = GFS_Sfcprop%usfco(im)
          mn_phys%vsfco(i,j)  = GFS_Sfcprop%vsfco(im)
          mn_phys%tsfco(i,j)  = GFS_Sfcprop%tsfco(im)
          mn_phys%tsfcl(i,j)  = GFS_Sfcprop%tsfcl(im)
          mn_phys%tsfc(i,j)   = GFS_Sfcprop%tsfc(im)

          mn_phys%albdirvis_lnd(i,j)   = GFS_Sfcprop%albdirvis_lnd(im)
          mn_phys%albdirnir_lnd(i,j)   = GFS_Sfcprop%albdirnir_lnd(im)
          mn_phys%albdifvis_lnd(i,j)   = GFS_Sfcprop%albdifvis_lnd(im)
          mn_phys%albdifnir_lnd(i,j)   = GFS_Sfcprop%albdifnir_lnd(im)

          do nv = 1, GFS_Control%ntot2d
            mn_phys%phy_f2d(i,j,nv) = GFS_tbd%phy_f2d(im, nv)
          enddo

          do k = 1, GFS_control%levs
            do nv = 1, GFS_control%ntot3d
              mn_phys%phy_f3d(i,j,k,nv) = GFS_tbd%phy_f3d(im, k, nv)
            enddo
          enddo

          ! Cloud prop data has x,y dimensions
          mn_phys%cv(i,j)  = GFS_cldprop%cv(im)
          mn_phys%cvt(i,j) = GFS_cldprop%cvt(im)
          mn_phys%cvb(i,j) = GFS_cldprop%cvb(im)
        endif

        if (move_nsst) then
          mn_phys%tref(i,j)   = GFS_sfcprop%tref(im)
          mn_phys%z_c(i,j)    = GFS_sfcprop%z_c(im)
          mn_phys%c_0(i,j)    = GFS_sfcprop%c_0(im)
          mn_phys%c_d(i,j)    = GFS_sfcprop%c_d(im)
          mn_phys%w_0(i,j)    = GFS_sfcprop%w_0(im)
          mn_phys%w_d(i,j)    = GFS_sfcprop%w_d(im)
          mn_phys%xt(i,j)     = GFS_sfcprop%xt(im)
          mn_phys%xs(i,j)     = GFS_sfcprop%xs(im)
          mn_phys%xu(i,j)     = GFS_sfcprop%xu(im)
          mn_phys%xv(i,j)     = GFS_sfcprop%xv(im)
          mn_phys%xz(i,j)     = GFS_sfcprop%xz(im)
          mn_phys%zm(i,j)     = GFS_sfcprop%zm(im)
          mn_phys%xtts(i,j)   = GFS_sfcprop%xtts(im)
          mn_phys%xzts(i,j)   = GFS_sfcprop%xzts(im)
          mn_phys%d_conv(i,j) = GFS_sfcprop%d_conv(im)
          mn_phys%dt_cool(i,j)= GFS_sfcprop%dt_cool(im)
          mn_phys%qrain(i,j)  = GFS_sfcprop%qrain(im)
        endif

        if (GFS_control%lsm == GFS_control%lsm_noahmp) then
          mn_phys%soilcolor(i,j)  = GFS_sfcprop%scolor(im)
          mn_phys%snowxy(i,j)     = GFS_sfcprop%snowxy(im)
          !if (i .eq. 149 .and. j .eq. 169) print '("[INFO] WDR SNOWXY MASK2D npe=",I0," i=",I0," j=",I0," snowxy=",E10.5)', this_pe, i, j, mn_phys%snowxy(i,j)


          mn_phys%tvxy(i,j)       = GFS_sfcprop%tvxy(im)
          mn_phys%tgxy(i,j)       = GFS_sfcprop%tgxy(im)
          mn_phys%canicexy(i,j)   = GFS_sfcprop%canicexy(im)
          mn_phys%canliqxy(i,j)   = GFS_sfcprop%canliqxy(im)
          mn_phys%eahxy(i,j)      = GFS_sfcprop%eahxy(im)
          mn_phys%tahxy(i,j)      = GFS_sfcprop%tahxy(im)
          mn_phys%cmxy(i,j)       = GFS_sfcprop%cmxy(im)
          mn_phys%chxy(i,j)       = GFS_sfcprop%chxy(im)
          mn_phys%fwetxy(i,j)     = GFS_sfcprop%fwetxy(im)
          mn_phys%sneqvoxy(i,j)   = GFS_sfcprop%sneqvoxy(im)
          mn_phys%alboldxy(i,j)   = GFS_sfcprop%alboldxy(im)
          mn_phys%qsnowxy(i,j)    = GFS_sfcprop%qsnowxy(im)
          mn_phys%wslakexy(i,j)   = GFS_sfcprop%wslakexy(im)
          mn_phys%zwtxy(i,j)      = GFS_sfcprop%zwtxy(im)
          mn_phys%waxy(i,j)       = GFS_sfcprop%waxy(im)
          mn_phys%wtxy(i,j)       = GFS_sfcprop%wtxy(im)
          mn_phys%lfmassxy(i,j)   = GFS_sfcprop%lfmassxy(im)
          mn_phys%rtmassxy(i,j)   = GFS_sfcprop%rtmassxy(im)
          mn_phys%stmassxy(i,j)   = GFS_sfcprop%stmassxy(im)
          mn_phys%woodxy(i,j)     = GFS_sfcprop%woodxy(im)
          mn_phys%stblcpxy(i,j)   = GFS_sfcprop%stblcpxy(im)
          mn_phys%fastcpxy(i,j)   = GFS_sfcprop%fastcpxy(im)
          mn_phys%xsaixy(i,j)     = GFS_sfcprop%xsaixy(im)
          mn_phys%xlaixy(i,j)     = GFS_sfcprop%xlaixy(im)
          mn_phys%taussxy(i,j)    = GFS_sfcprop%taussxy(im)
          mn_phys%smcwtdxy(i,j)   = GFS_sfcprop%smcwtdxy(im)
          mn_phys%deeprechxy(i,j) = GFS_sfcprop%deeprechxy(im)
          mn_phys%rechxy(i,j)     = GFS_sfcprop%rechxy(im)

          do k = 1, GFS_control%lsoil
             mn_phys%smoiseq(i,j,k) = GFS_sfcprop%smoiseq(im,k)
          enddo

          ! lsnow_lsm_lbound is a negative value, lsnow_ubound is usually 0
          do k = GFS_control%lsnow_lsm_lbound, GFS_control%lsnow_lsm_ubound
            mn_phys%snicexy(i,j,k)    = GFS_sfcprop%snicexy(im,k)
            mn_phys%snliqxy(i,j,k)    = GFS_sfcprop%snliqxy(im,k)
            mn_phys%tsnoxy(i,j,k)     = GFS_sfcprop%tsnoxy(im,k)
          enddo

          mn_phys%snowd(i,j)      = GFS_sfcprop%snowd(im)
          mn_phys%weasd(i,j)      = GFS_sfcprop%weasd(im)

          do k = GFS_control%lsnow_lsm_lbound, GFS_control%lsoil
            mn_phys%zsnsoxy(i,j,k)    = GFS_sfcprop%zsnsoxy(im,k)
          enddo
        endif
      enddo
    enddo

  end subroutine mn_phys_fill_temp_variables

  !>@brief The subroutine 'mn_phys_apply_temp_variables' copies moved 2D data back into 1D physics arryas for nest motion
  !>@details This subroutine fills the 1D physics arrays from the mn_phys structure on the Atm object
  !!  Note that ice variables are not yet handled.
  subroutine mn_phys_apply_temp_variables(Atm, Atm_block, GFS_control, GFS_sfcprop, GFS_tbd, GFS_cldprop, GFS_intdiag, n, child_grid_num, is_fine_pe, npz)
    type(fv_atmos_type), allocatable, target, intent(inout)  :: Atm(:)            !< Array of atmospheric data
    type (block_control_type), intent(in)                    :: Atm_block         !< Physics block layout
    type(GFS_control_type), intent(in)                       :: GFS_control       !< Physics metadata
    type(GFS_sfcprop_type), intent(inout)                    :: GFS_sfcprop       !< Physics variable data (surface)
    type(GFS_tbd_type), intent(inout)                        :: GFS_tbd           !< Physics variable data (tbd)
    type(GFS_cldprop_type), intent(inout)                    :: GFS_cldprop       !< Physics variable data (clouds)
    type(GFS_diag_type), intent(inout)                       :: GFS_intdiag       !< Physics variable data (diagnostic)
    integer, intent(in)                                      :: n, child_grid_num !< Current grid number, child grid number
    logical, intent(in)                                      :: is_fine_pe        !< Is this a nest PE?
    integer, intent(in)                                      :: npz               !< Number of vertical levels
    integer, intent(in)                                      :: a_step            !< Timestep for logging

    integer :: is, ie, js, je
    integer :: this_pe
    integer :: nb, blen, i, j ,k, ix, nv, im
    type(fv_moving_nest_physics_type), pointer       :: mn_phys

    this_pe = mpp_pe()
    mn_phys => Moving_nest(n)%mn_phys

    !  Needed to fill the local grids for parent and nest PEs in order to transmit/interpolate data from parent to nest
    !  But only the nest PE's have changed the values with nest motion, so they are the only ones that need to update the original arrays
    if (is_fine_pe) then
      is = Atm(n)%bd%is
      ie = Atm(n)%bd%ie
      js = Atm(n)%bd%js
      je = Atm(n)%bd%je

      ! SST directly in Atm structure
      Atm(n)%ts(is:ie, js:je) =  mn_phys%ts(is:ie, js:je)

      im = 0
      do nb = 1,Atm_block%nblks
        blen = Atm_block%blksz(nb)
        do ix = 1, blen
          i = Atm_block%index(nb)%ii(ix)
          j = Atm_block%index(nb)%jj(ix)
          im = im +1

          if (move_physics) then
            ! Surface properties
            do k = 1, GFS_control%lsoil
              GFS_sfcprop%smc(im,k) = mn_phys%smc(i,j,k)
              GFS_sfcprop%stc(im,k) = mn_phys%stc(i,j,k)
              GFS_sfcprop%slc(im,k) = mn_phys%slc(i,j,k)
            enddo

            ! EMIS PATCH - Force to positive at all locations.
            if (mn_phys%emis_lnd(i,j) .ge. 0.0) then
              GFS_sfcprop%emis_lnd(im) = mn_phys%emis_lnd(i,j)
            else
              GFS_sfcprop%emis_lnd(im) = 0.5
            endif
            if (mn_phys%emis_ice(i,j) .ge. 0.0) then
              GFS_sfcprop%emis_ice(im) = mn_phys%emis_ice(i,j)
            else
              GFS_sfcprop%emis_ice(im) = 0.5
            endif
            if (mn_phys%emis_wat(i,j) .ge. 0.0) then
              GFS_sfcprop%emis_wat(im) = mn_phys%emis_wat(i,j)
            else
              GFS_sfcprop%emis_wat(im) = 0.5
            endif

            !GFS_sfcprop%sfalb_lnd(im) = mn_phys%sfalb_lnd(i,j)
            !GFS_sfcprop%sfalb_lnd_bck(im) = mn_phys%sfalb_lnd_bck(i,j)
            !GFS_radtend%semis(im) = mn_phys%semis(i,j)
            !GFS_sfcprop%semisbase(im) = mn_phys%semisbase(i,j)
            !GFS_radtend%sfalb(im) = mn_phys%sfalb(i,j)

            GFS_intdiag%u10m(im) = mn_phys%u10m(i,j)
            GFS_intdiag%v10m(im) = mn_phys%v10m(i,j)
            GFS_sfcprop%tprcp(im) = mn_phys%tprcp(i,j)

            do k = 1, GFS_control%nmtvr
              GFS_sfcprop%hprime(im,k) = mn_phys%hprime(i,j,k)
            enddo

            GFS_sfcprop%lakefrac(im) = mn_phys%lakefrac(i,j)
            GFS_sfcprop%lakedepth(im) = mn_phys%lakedepth(i,j)

            GFS_sfcprop%canopy(im) = mn_phys%canopy(i,j)
            GFS_sfcprop%vfrac(im)  = mn_phys%vegfrac(i,j)
            GFS_sfcprop%uustar(im) = mn_phys%uustar(i,j)
            GFS_sfcprop%shdmin(im) = mn_phys%shdmin(i,j)
            GFS_sfcprop%shdmax(im) = mn_phys%shdmax(i,j)

            ! Set roughness lengths to physically reasonable values if they have fill value (possible at coastline)
            ! sea/land mask array (sea:0,land:1,sea-ice:2)
            if (nint(GFS_sfcprop%slmsk(im)) .eq. 1 .and. mn_phys%zorll(i,j) .gt. 1e6) then
              GFS_sfcprop%zorll(im)  = 82.0   !
            else
              GFS_sfcprop%zorll(im)  = mn_phys%zorll(i,j)
            endif

            if (nint(GFS_sfcprop%slmsk(im)) .eq. 0 .and. mn_phys%zorlw(i,j) .gt. 1e6) then
              GFS_sfcprop%zorlw(im)  = 83.0   !
            else
              GFS_sfcprop%zorlw(im)  = mn_phys%zorlw(i,j)
            endif

            if (nint(GFS_sfcprop%slmsk(im)) .eq. 0 .and. mn_phys%zorlwav(i,j) .gt. 1e6) then
              GFS_sfcprop%zorlwav(im)  = 84.0   !
            else
              GFS_sfcprop%zorlwav(im)  = mn_phys%zorlwav(i,j)
            endif

            if (mn_phys%zorl(i,j) .gt. 1e6) then
              GFS_sfcprop%zorl(im)   = 85.0
            else
              GFS_sfcprop%zorl(im)   = mn_phys%zorl(i,j)
            endif

            if (nint(GFS_sfcprop%slmsk(im)) .eq. 0 .and. mn_phys%usfco(i,j) .gt. 1e6) then
              GFS_sfcprop%usfco(im)  = 0.0
            else
              GFS_sfcprop%usfco(im)  = mn_phys%usfco(i,j)
            endif
            if (nint(GFS_sfcprop%slmsk(im)) .eq. 0 .and. mn_phys%vsfco(i,j) .gt. 1e6) then
              GFS_sfcprop%vsfco(im)  = 0.0
            else
              GFS_sfcprop%vsfco(im)  = mn_phys%vsfco(i,j)
            endif

            GFS_sfcprop%tsfco(im)  = mn_phys%tsfco(i,j)
            GFS_sfcprop%tsfcl(im)  = mn_phys%tsfcl(i,j)
            GFS_sfcprop%tsfc(im)   = mn_phys%tsfc(i,j)

            ! Set albedo values to physically reasonable values if they have negative fill values.
            if (mn_phys%albdirvis_lnd (i,j) .ge. 0.0) then
              GFS_sfcprop%albdirvis_lnd (im)   = mn_phys%albdirvis_lnd (i,j)
            else
              GFS_sfcprop%albdirvis_lnd (im)   = 0.5
            endif

            if (mn_phys%albdirnir_lnd (i,j) .ge. 0.0) then
              GFS_sfcprop%albdirnir_lnd (im)   = mn_phys%albdirnir_lnd (i,j)
            else
              GFS_sfcprop%albdirnir_lnd (im)   = 0.5
            endif

            if (mn_phys%albdifvis_lnd (i,j) .ge. 0.0) then
              GFS_sfcprop%albdifvis_lnd (im)   = mn_phys%albdifvis_lnd (i,j)
            else
              GFS_sfcprop%albdifvis_lnd (im)   = 0.5
            endif

            if (mn_phys%albdifnir_lnd (i,j) .ge. 0.0) then
              GFS_sfcprop%albdifnir_lnd (im)   = mn_phys%albdifnir_lnd (i,j)
            else
              GFS_sfcprop%albdifnir_lnd (im)   = 0.5
            endif

            ! Cloud properties
            GFS_cldprop%cv(im) = mn_phys%cv(i,j)
            GFS_cldprop%cvt(im) = mn_phys%cvt(i,j)
            GFS_cldprop%cvb(im) = mn_phys%cvb(i,j)

            do nv = 1, GFS_control%ntot2d
              GFS_tbd%phy_f2d(im, nv) = mn_phys%phy_f2d(i,j,nv)
            enddo

            do k = 1, GFS_control%levs
              do nv = 1, GFS_control%ntot3d
                GFS_tbd%phy_f3d(im, k, nv) = mn_phys%phy_f3d(i,j,k,nv)
              enddo
            enddo
          endif

          if (move_nsst) then
            GFS_sfcprop%tref(im)    = mn_phys%tref(i,j)
            GFS_sfcprop%z_c(im)     = mn_phys%z_c(i,j)
            GFS_sfcprop%c_0(im)     = mn_phys%c_0(i,j)
            GFS_sfcprop%c_d(im)     = mn_phys%c_d(i,j)
            GFS_sfcprop%w_0(im)     = mn_phys%w_0(i,j)
            GFS_sfcprop%w_d(im)     = mn_phys%w_d(i,j)
            GFS_sfcprop%xt(im)      = mn_phys%xt(i,j)
            GFS_sfcprop%xs(im)      = mn_phys%xs(i,j)
            GFS_sfcprop%xu(im)      = mn_phys%xu(i,j)
            GFS_sfcprop%xv(im)      = mn_phys%xv(i,j)
            GFS_sfcprop%xz(im)      = mn_phys%xz(i,j)
            GFS_sfcprop%zm(im)      = mn_phys%zm(i,j)
            GFS_sfcprop%xtts(im)    = mn_phys%xtts(i,j)
            GFS_sfcprop%xzts(im)    = mn_phys%xzts(i,j)
            GFS_sfcprop%d_conv(im)  = mn_phys%d_conv(i,j)
            GFS_sfcprop%dt_cool(im) = mn_phys%dt_cool(i,j)
            GFS_sfcprop%qrain(im)   = mn_phys%qrain(i,j)
          endif

          if (GFS_control%lsm == GFS_control%lsm_noahmp) then

            GFS_sfcprop%scolor(im)  = mn_phys%soilcolor(i,j)
            GFS_sfcprop%snowxy(im)     = mn_phys%snowxy(i,j)
            GFS_sfcprop%tvxy(im)       = mn_phys%tvxy(i,j)
            GFS_sfcprop%tgxy(im)       = mn_phys%tgxy(i,j)
            GFS_sfcprop%canicexy(im)   = mn_phys%canicexy(i,j)
            GFS_sfcprop%canliqxy(im)   = mn_phys%canliqxy(i,j)
            GFS_sfcprop%eahxy(im)      = mn_phys%eahxy(i,j)
            GFS_sfcprop%tahxy(im)      = mn_phys%tahxy(i,j)
            GFS_sfcprop%cmxy(im)       = mn_phys%cmxy(i,j)
            GFS_sfcprop%chxy(im)       = mn_phys%chxy(i,j)
            GFS_sfcprop%fwetxy(im)     = mn_phys%fwetxy(i,j)
            GFS_sfcprop%sneqvoxy(im)   = mn_phys%sneqvoxy(i,j)
            GFS_sfcprop%alboldxy(im)   = mn_phys%alboldxy(i,j)
            GFS_sfcprop%qsnowxy(im)    = mn_phys%qsnowxy(i,j)
            GFS_sfcprop%wslakexy(im)   = mn_phys%wslakexy(i,j)
            GFS_sfcprop%zwtxy(im)      = mn_phys%zwtxy(i,j)
            GFS_sfcprop%waxy(im)       = mn_phys%waxy(i,j)
            GFS_sfcprop%wtxy(im)       = mn_phys%wtxy(i,j)
            GFS_sfcprop%lfmassxy(im)   = mn_phys%lfmassxy(i,j)
            GFS_sfcprop%rtmassxy(im)   = mn_phys%rtmassxy(i,j)
            GFS_sfcprop%stmassxy(im)   = mn_phys%stmassxy(i,j)
            GFS_sfcprop%woodxy(im)     = mn_phys%woodxy(i,j)
            GFS_sfcprop%stblcpxy(im)   = mn_phys%stblcpxy(i,j)
            GFS_sfcprop%fastcpxy(im)   = mn_phys%fastcpxy(i,j)
            GFS_sfcprop%xsaixy(im)     = mn_phys%xsaixy(i,j)
            GFS_sfcprop%xlaixy(im)     = mn_phys%xlaixy(i,j)
            GFS_sfcprop%taussxy(im)    = mn_phys%taussxy(i,j)
            GFS_sfcprop%smcwtdxy(im)   = mn_phys%smcwtdxy(i,j)
            GFS_sfcprop%deeprechxy(im) = mn_phys%deeprechxy(i,j)
            GFS_sfcprop%rechxy(im)     = mn_phys%rechxy(i,j)

            do k = 1, GFS_control%lsoil
               GFS_sfcprop%smoiseq(im,k)    = mn_phys%smoiseq(i,j,k)
            enddo


            do k = GFS_control%lsnow_lsm_lbound, GFS_control%lsnow_lsm_ubound
              GFS_sfcprop%snicexy(im,k)    = mn_phys%snicexy(i,j,k)
              GFS_sfcprop%snliqxy(im,k)    = mn_phys%snliqxy(i,j,k)
              GFS_sfcprop%tsnoxy(im,k)     = mn_phys%tsnoxy(i,j,k)
            enddo

            GFS_sfcprop%snowd(im)      = mn_phys%snowd(i,j)
            GFS_sfcprop%weasd(im)      = mn_phys%weasd(i,j)

            do k = GFS_control%lsnow_lsm_lbound, GFS_control%lsoil
              GFS_sfcprop%zsnsoxy(im,k)    = mn_phys%zsnsoxy(i,j,k)
            enddo

          endif

          ! Check if stype and vtype are properly set for land points.  Set to reasonable values if they have fill values.
          if ( (int(GFS_sfcprop%slmsk(im)) .eq. 1) )  then

            if (GFS_sfcprop%vtype(im) .lt. 0.5) then
              GFS_sfcprop%vtype(im) = 7    ! Force to grassland
            endif

            if (GFS_sfcprop%stype(im) .lt. 0.5) then
              GFS_sfcprop%stype(im) = 3    ! Force to sandy loam
            endif

            if (GFS_sfcprop%vtype_save(im) .lt. 0.5) then
              GFS_sfcprop%vtype_save(im) = 7    ! Force to grassland
            endif
            if (GFS_sfcprop%stype_save(im) .lt. 0.5) then
              GFS_sfcprop%stype_save(im) = 3    ! Force to sandy loam
            endif

          endif
        enddo
      enddo
    endif

  end subroutine mn_phys_apply_temp_variables


  !>@brief The subroutine 'mn_physfill_nest_halos_from_parent' transfers data from the coarse grid to the nest edge
  !>@details This subroutine must run on parent and nest PEs to complete the data transfers
  subroutine mn_phys_fill_nest_halos_from_parent(Atm, GFS_control, mn_static, n, child_grid_num, is_fine_pe, nest_domain, nz)
    type(fv_atmos_type), allocatable, target, intent(inout)  :: Atm(:)            !< Array of atmospheric data
    type(GFS_control_type), intent(in)                       :: GFS_control       !< Physics metadata
    type(mn_surface_grids), intent(in)                       :: mn_static         !< Static data
    integer, intent(in)                                      :: n, child_grid_num !< Current grid number, child grid number
    logical, intent(in)                                      :: is_fine_pe        !< Is this a nest PE?
    type(nest_domain_type), intent(inout)                    :: nest_domain       !< Nest domain for FMS
    integer, intent(in)                                      :: nz                !< Number of vertical levels

    integer  :: position, position_u, position_v
    integer  :: interp_type, interp_type_u, interp_type_v, interp_type_lmask
    integer  :: x_refine, y_refine
    type(fv_moving_nest_physics_type), pointer :: mn_phys

    integer, parameter :: M_WATER = 0, M_LAND = 1, M_SEAICE = 2
    !! For NOAHMP
    ! (/0.0, 0.0, 0.0,  0.1,0.4,1.0,2.0/) -- 3 snow levels, 4 soil levels
    ! TODO make this more flexible for number of snow and soil levels
      !do k = GFS_control%lsnow_lsm_lbound, GFS_control%lsoil
    real(kind=kind_phys) :: zsns_default(-2:4)

    if (GFS_control%lsm == GFS_control%lsm_noahmp) then
      zsns_default = [0.0, 0.0, 0.0,  -0.1,-0.4,-1.0,-2.0 ]
    else
      ! Expect that zsns_default is not used in this case, but just to be safe, set to 0
      zsns_default = 0.0
    endif


    interp_type = 1        ! cell-centered A-grid
    interp_type_u = 4      ! D-grid
    interp_type_v = 4      ! D-grid
    interp_type_lmask = 7  ! land mask, cell-centered A-grid

    position = CENTER
    position_u = NORTH
    position_v = EAST

    x_refine = Atm(child_grid_num)%neststruct%refinement
    y_refine = x_refine

    mn_phys => Moving_nest(n)%mn_phys

    !  Fill centered-grid variables

    call fill_nest_halos_from_parent("ts", mn_phys%ts, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
        Atm(child_grid_num)%neststruct%ind_h, &
        x_refine, y_refine, &
        is_fine_pe, nest_domain, position)

    if (move_physics) then
      ! Default - Arbitrary value 0.3
      call fill_nest_halos_from_parent_masked("smc", mn_phys%smc, interp_type_lmask, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position, 1, GFS_Control%lsoil, mn_phys%slmsk, mn_static%parent_ls%ls_mask_grid, M_LAND, 0.3D0)
      ! Defaults - use surface temperature to set soil temperature at each level
      call fill_nest_halos_from_parent_masked("stc", mn_phys%stc, interp_type_lmask, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position, 1, GFS_Control%lsoil, mn_phys%slmsk, mn_static%parent_ls%ls_mask_grid, M_LAND, mn_phys%ts)
      ! Default - Arbitrary value 0.3
      call fill_nest_halos_from_parent_masked("slc", mn_phys%slc, interp_type_lmask, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position, 1, GFS_Control%lsoil, mn_phys%slmsk, mn_static%parent_ls%ls_mask_grid, M_LAND, 0.3D0)

      call fill_nest_halos_from_parent("phy_f2d", mn_phys%phy_f2d, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position, GFS_control%ntot2d)

      call fill_nest_halos_from_parent("phy_f3d", mn_phys%phy_f3d, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position, GFS_control%levs)

      !!  Surface variables

      !call fill_nest_halos_from_parent("sfalb_lnd", mn_phys%sfalb_lnd, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
      !     Atm(child_grid_num)%neststruct%ind_h, &
      !     x_refine, y_refine, &
      !     is_fine_pe, nest_domain, position)

      ! sea/land mask array (sea:0,land:1,sea-ice:2)
      !integer, parameter :: M_WATER = 0, M_LAND = 1, M_SEAICE = 2

      call fill_nest_halos_from_parent_masked("emis_lnd", mn_phys%emis_lnd, interp_type_lmask, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position, mn_phys%slmsk, mn_static%parent_ls%ls_mask_grid, M_LAND, 0.5D0)

      call fill_nest_halos_from_parent_masked("emis_ice", mn_phys%emis_ice, interp_type_lmask, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position, mn_phys%slmsk, mn_static%parent_ls%ls_mask_grid, M_SEAICE, 0.5D0)

      call fill_nest_halos_from_parent_masked("emis_wat", mn_phys%emis_wat, interp_type_lmask, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position, mn_phys%slmsk, mn_static%parent_ls%ls_mask_grid, M_WATER, 0.5D0)

      !call fill_nest_halos_from_parent("sfalb_lnd_bck", mn_phys%sfalb_lnd_bck, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
      !     Atm(child_grid_num)%neststruct%ind_h, &
      !     x_refine, y_refine, &
      !     is_fine_pe, nest_domain, position)


      !call fill_nest_halos_from_parent("semis", mn_phys%semis, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
      !     Atm(child_grid_num)%neststruct%ind_h, &
      !     x_refine, y_refine, &
      !     is_fine_pe, nest_domain, position)
      !call fill_nest_halos_from_parent("semisbase", mn_phys%semisbase, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
      !     Atm(child_grid_num)%neststruct%ind_h, &
      !     x_refine, y_refine, &
      !     is_fine_pe, nest_domain, position)
      !call fill_nest_halos_from_parent("sfalb", mn_phys%sfalb, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
      !     Atm(child_grid_num)%neststruct%ind_h, &
      !     x_refine, y_refine, &
      !     is_fine_pe, nest_domain, position)


      call fill_nest_halos_from_parent("u10m", mn_phys%u10m, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position)
      call fill_nest_halos_from_parent("v10m", mn_phys%v10m, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position)
      call fill_nest_halos_from_parent("tprcp", mn_phys%tprcp, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position)

      call fill_nest_halos_from_parent("hprime", mn_phys%hprime, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position, GFS_control%nmtvr)

      call fill_nest_halos_from_parent("lakefrac", mn_phys%lakefrac, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position)
      call fill_nest_halos_from_parent("lakedepth", mn_phys%lakedepth, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position)

!      call fill_nest_halos_from_parent("canopy", mn_phys%canopy, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
!          Atm(child_grid_num)%neststruct%ind_h, &
!          x_refine, y_refine, &
!          is_fine_pe, nest_domain, position)
!      call fill_nest_halos_from_parent("vegfrac", mn_phys%vegfrac, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
!          Atm(child_grid_num)%neststruct%ind_h, &
!          x_refine, y_refine, &
!          is_fine_pe, nest_domain, position)
      call fill_nest_halos_from_parent_masked("canopy", mn_phys%canopy, interp_type_lmask, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, x_refine, y_refine, &
          is_fine_pe, nest_domain, position, mn_phys%slmsk, mn_static%parent_ls%ls_mask_grid, M_LAND, 0.0D0)
      call fill_nest_halos_from_parent_masked("vegfrac", mn_phys%vegfrac, interp_type_lmask, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, x_refine, y_refine, &
          is_fine_pe, nest_domain, position, mn_phys%slmsk, mn_static%parent_ls%ls_mask_grid, M_LAND, 0.50D0)



      call fill_nest_halos_from_parent("uustar", mn_phys%uustar, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position)
      call fill_nest_halos_from_parent("shdmin", mn_phys%shdmin, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position)
      call fill_nest_halos_from_parent("shdmax", mn_phys%shdmax, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position)
      call fill_nest_halos_from_parent("zorl", mn_phys%zorl, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position)

      call fill_nest_halos_from_parent_masked("zorll", mn_phys%zorll, interp_type_lmask, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position, mn_phys%slmsk, mn_static%parent_ls%ls_mask_grid, 1, 86.0D0)
      call fill_nest_halos_from_parent_masked("zorlwav", mn_phys%zorlwav, interp_type_lmask, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position, mn_phys%slmsk, mn_static%parent_ls%ls_mask_grid, 0, 77.0D0)
      call fill_nest_halos_from_parent_masked("zorlw", mn_phys%zorlw, interp_type_lmask, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position, mn_phys%slmsk, mn_static%parent_ls%ls_mask_grid, 0, 78.0D0)

      call fill_nest_halos_from_parent_masked("usfco", mn_phys%usfco, interp_type_lmask, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position, mn_phys%slmsk, mn_static%parent_ls%ls_mask_grid, 0, 0.0D0)
      call fill_nest_halos_from_parent_masked("vsfco", mn_phys%vsfco, interp_type_lmask, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position, mn_phys%slmsk, mn_static%parent_ls%ls_mask_grid, 0, 0.0D0)

      call fill_nest_halos_from_parent("tsfco", mn_phys%tsfco, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position)
      call fill_nest_halos_from_parent("tsfcl", mn_phys%tsfcl, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position)
      call fill_nest_halos_from_parent("tsfc", mn_phys%tsfc, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position)

      call fill_nest_halos_from_parent_masked("albdirvis_lnd", mn_phys%albdirvis_lnd, interp_type_lmask, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position, mn_phys%slmsk, mn_static%parent_ls%ls_mask_grid, M_LAND, 0.5D0)
      call fill_nest_halos_from_parent_masked("albdirnir_lnd", mn_phys%albdirnir_lnd, interp_type_lmask, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position, mn_phys%slmsk, mn_static%parent_ls%ls_mask_grid, M_LAND, 0.5D0)
      call fill_nest_halos_from_parent_masked("albdifvis_lnd", mn_phys%albdifvis_lnd, interp_type_lmask, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position, mn_phys%slmsk, mn_static%parent_ls%ls_mask_grid, M_LAND, 0.5D0)
      call fill_nest_halos_from_parent_masked("albdifnir_lnd", mn_phys%albdifnir_lnd, interp_type_lmask, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position, mn_phys%slmsk, mn_static%parent_ls%ls_mask_grid, M_LAND, 0.5D0)



      call fill_nest_halos_from_parent("cv", mn_phys%cv, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position)
      call fill_nest_halos_from_parent("cvt", mn_phys%cvt, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position)
      call fill_nest_halos_from_parent("cvb", mn_phys%cvb, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position)
    endif

    if (move_nsst) then

      call fill_nest_halos_from_parent("tref", mn_phys%tref, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position)
      call fill_nest_halos_from_parent("z_c", mn_phys%z_c, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position)
      call fill_nest_halos_from_parent("c_0", mn_phys%c_0, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position)
      call fill_nest_halos_from_parent("c_d", mn_phys%c_d, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position)
      call fill_nest_halos_from_parent("w_0", mn_phys%w_0, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position)
      call fill_nest_halos_from_parent("w_d", mn_phys%w_d, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position)
      call fill_nest_halos_from_parent("xt", mn_phys%xt, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position)
      call fill_nest_halos_from_parent("xs", mn_phys%xs, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position)
      call fill_nest_halos_from_parent("xu", mn_phys%xu, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position)
      call fill_nest_halos_from_parent("xv", mn_phys%xv, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position)
      call fill_nest_halos_from_parent("xz", mn_phys%xz, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position)
      call fill_nest_halos_from_parent("zm", mn_phys%zm, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position)
      call fill_nest_halos_from_parent("xtts", mn_phys%xtts, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position)
      call fill_nest_halos_from_parent("xzts", mn_phys%xzts, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position)
      call fill_nest_halos_from_parent("d_conv", mn_phys%d_conv, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position)
      call fill_nest_halos_from_parent("dt_cool", mn_phys%dt_cool, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position)
      call fill_nest_halos_from_parent("qrain", mn_phys%qrain, interp_type, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, &
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position)

    endif

    if (move_physics .and. GFS_control%lsm == GFS_control%lsm_noahmp) then

      !integer, parameter :: M_WATER = 0, M_LAND = 1, M_SEAICE = 2

      !  Land Sea Mask has values of 0 for oceans/lakes, 1 for land, 2 for sea ice

      ! Soil color.  Default is set to sandy soil/desert 1, which seems appropriate for isolated islands
      !  Reference: https://www.jsg.utexas.edu/noah-mp/files/Users_Guide_v0.pdf
      !  Default changed to 10 based on suggestion from Mike Barlage; more middle of the spectrum value.
      call fill_nest_halos_from_parent_masked("soilcol", mn_phys%soilcolor, interp_type_lmask, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, x_refine, y_refine, &
          is_fine_pe, nest_domain, position, mn_phys%slmsk, mn_static%parent_ls%ls_mask_grid, M_LAND, 10.0D0)

      call fill_nest_halos_from_parent_masked("snowxy", mn_phys%snowxy, interp_type_lmask, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, x_refine, y_refine, &
          is_fine_pe, nest_domain, position, mn_phys%slmsk, mn_static%parent_ls%ls_mask_grid, M_LAND, 0.0D0)
      call fill_nest_halos_from_parent_masked("tvxy", mn_phys%tvxy, interp_type_lmask, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, x_refine, y_refine, &
          is_fine_pe, nest_domain, position, mn_phys%slmsk, mn_static%parent_ls%ls_mask_grid, M_LAND, mn_phys%ts)
      call fill_nest_halos_from_parent_masked("tgxy", mn_phys%tgxy, interp_type_lmask, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, x_refine, y_refine, &
          is_fine_pe, nest_domain, position, mn_phys%slmsk, mn_static%parent_ls%ls_mask_grid, M_LAND, mn_phys%ts)

      call fill_nest_halos_from_parent_masked("canicexy", mn_phys%canicexy, interp_type_lmask, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, x_refine, y_refine, &
          is_fine_pe, nest_domain, position, mn_phys%slmsk, mn_static%parent_ls%ls_mask_grid, M_LAND, 0.0D0)
      call fill_nest_halos_from_parent_masked("canliqxy", mn_phys%canliqxy, interp_type_lmask, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, x_refine, y_refine, &
          is_fine_pe, nest_domain, position, mn_phys%slmsk, mn_static%parent_ls%ls_mask_grid, M_LAND, 0.0D0)

      call fill_nest_halos_from_parent_masked("eahxy", mn_phys%eahxy, interp_type_lmask, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, x_refine, y_refine, &
          is_fine_pe, nest_domain, position, mn_phys%slmsk, mn_static%parent_ls%ls_mask_grid, M_LAND, 2000.0D0)

      call fill_nest_halos_from_parent_masked("tahxy", mn_phys%tahxy, interp_type_lmask, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, x_refine, y_refine, &
          is_fine_pe, nest_domain, position, mn_phys%slmsk, mn_static%parent_ls%ls_mask_grid, M_LAND, mn_phys%ts)

      ! TODO get realistic default value here  -- bulk momentum drag coefficient
      call fill_nest_halos_from_parent_masked("cmxy", mn_phys%cmxy, interp_type_lmask, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, x_refine, y_refine, &
          is_fine_pe, nest_domain, position, mn_phys%slmsk, mn_static%parent_ls%ls_mask_grid, M_LAND, 2.4D-3)

      ! TODO get realistic default value here  -- bulk sensible heat drag coefficient
      call fill_nest_halos_from_parent_masked("chxy", mn_phys%chxy, interp_type_lmask, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, x_refine, y_refine, &
          is_fine_pe, nest_domain, position, mn_phys%slmsk, mn_static%parent_ls%ls_mask_grid, M_LAND, 2.4D-3)

      ! wetted or snowed fraction of the canopy
      call fill_nest_halos_from_parent_masked("fwetxy", mn_phys%fwetxy, interp_type_lmask, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, x_refine, y_refine, &
          is_fine_pe, nest_domain, position, mn_phys%slmsk, mn_static%parent_ls%ls_mask_grid, M_LAND, 0.0D0)

      ! snow mass at last time step[mm h2o]
      call fill_nest_halos_from_parent_masked("sneqvoxy", mn_phys%sneqvoxy, interp_type_lmask, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, x_refine, y_refine, &
          is_fine_pe, nest_domain, position, mn_phys%slmsk, mn_static%parent_ls%ls_mask_grid, M_LAND, 0.0D0)

      ! Albedo assuming deep snow on prev timestep - default to 0.65
      call fill_nest_halos_from_parent_masked("alboldxy", mn_phys%alboldxy, interp_type_lmask, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, x_refine, y_refine, &
          is_fine_pe, nest_domain, position, mn_phys%slmsk, mn_static%parent_ls%ls_mask_grid, M_LAND, 0.65D0)

      ! Liquid equivalent snow - default to 0
      call fill_nest_halos_from_parent_masked("qsnowxy", mn_phys%qsnowxy, interp_type_lmask, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, x_refine, y_refine, &
          is_fine_pe, nest_domain, position, mn_phys%slmsk, mn_static%parent_ls%ls_mask_grid, M_LAND, 0.0D0)

      ! Lake water storage [mm] -- TODO find better default
      call fill_nest_halos_from_parent_masked("wslakexy", mn_phys%wslakexy, interp_type_lmask, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, x_refine, y_refine, &
          is_fine_pe, nest_domain, position, mn_phys%slmsk, mn_static%parent_ls%ls_mask_grid, M_LAND, 0.0D0)

      ! Water table depth - set to 2.5, cold start value
      call fill_nest_halos_from_parent_masked("zwtxy", mn_phys%zwtxy, interp_type_lmask, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, x_refine, y_refine, &
          is_fine_pe, nest_domain, position, mn_phys%slmsk, mn_static%parent_ls%ls_mask_grid, M_LAND, 2.5D0)

      ! Water storage in aquifer - set to 4900.0, cold start value
      call fill_nest_halos_from_parent_masked("waxy", mn_phys%waxy, interp_type_lmask, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, x_refine, y_refine, &
          is_fine_pe, nest_domain, position, mn_phys%slmsk, mn_static%parent_ls%ls_mask_grid, M_LAND, 4900.0D0)
      ! Water storage in aquifer and saturated soil - set to 4900.0, cold start value
      call fill_nest_halos_from_parent_masked("wtxy", mn_phys%wtxy, interp_type_lmask, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, x_refine, y_refine, &
          is_fine_pe, nest_domain, position, mn_phys%slmsk, mn_static%parent_ls%ls_mask_grid, M_LAND, 4900.0D0)


      ! Leaf mass [g/m2] -- TODO find better default
      call fill_nest_halos_from_parent_masked("lfmassxy", mn_phys%lfmassxy, interp_type_lmask, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, x_refine, y_refine, &
          is_fine_pe, nest_domain, position, mn_phys%slmsk, mn_static%parent_ls%ls_mask_grid, M_LAND, 0.0D0)
      ! Fine root mass [g/m2] -- TODO find better default
      call fill_nest_halos_from_parent_masked("rtmassxy", mn_phys%rtmassxy, interp_type_lmask, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, x_refine, y_refine, &
          is_fine_pe, nest_domain, position, mn_phys%slmsk, mn_static%parent_ls%ls_mask_grid, M_LAND, 0.0D0)
      ! Stem mass [g/m2] -- TODO find better default
      call fill_nest_halos_from_parent_masked("stmassxy", mn_phys%stmassxy, interp_type_lmask, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, x_refine, y_refine, &
          is_fine_pe, nest_domain, position, mn_phys%slmsk, mn_static%parent_ls%ls_mask_grid, M_LAND, 0.0D0)
      ! Wood mass [g/m2] -- TODO find better default
      call fill_nest_halos_from_parent_masked("woodxy", mn_phys%woodxy, interp_type_lmask, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, x_refine, y_refine, &
          is_fine_pe, nest_domain, position, mn_phys%slmsk, mn_static%parent_ls%ls_mask_grid, M_LAND, 0.0D0)

      ! stable carbon in deep soil [g/m2] -- TODO find a better default
      call fill_nest_halos_from_parent_masked("stblcpxy", mn_phys%stblcpxy, interp_type_lmask, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, x_refine, y_refine, &
          is_fine_pe, nest_domain, position, mn_phys%slmsk, mn_static%parent_ls%ls_mask_grid, M_LAND, 0.0D0)
      ! short-lived carbon, shallow soil [g/m2] -- TODO find a better default
      call fill_nest_halos_from_parent_masked("fastcpxy", mn_phys%fastcpxy, interp_type_lmask, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, x_refine, y_refine, &
          is_fine_pe, nest_domain, position, mn_phys%slmsk, mn_static%parent_ls%ls_mask_grid, M_LAND, 0.0D0)

      ! stem area index [m2/m2] -- TODO find a better default
      call fill_nest_halos_from_parent_masked("xsaixy", mn_phys%xsaixy, interp_type_lmask, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, x_refine, y_refine, &
          is_fine_pe, nest_domain, position, mn_phys%slmsk, mn_static%parent_ls%ls_mask_grid, M_LAND, 0.0D0)
      ! leaf area index [m2/m2] -- TODO find a better default
      call fill_nest_halos_from_parent_masked("xlaixy", mn_phys%xlaixy, interp_type_lmask, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, x_refine, y_refine, &
          is_fine_pe, nest_domain, position, mn_phys%slmsk, mn_static%parent_ls%ls_mask_grid, M_LAND, 0.0D0)

      ! snow age factor [-] -- TODO find a better default
      call fill_nest_halos_from_parent_masked("taussxy", mn_phys%taussxy, interp_type_lmask, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, x_refine, y_refine, &
          is_fine_pe, nest_domain, position, mn_phys%slmsk, mn_static%parent_ls%ls_mask_grid, M_LAND, 0.0D0)

      ! soil moisture content in the layer to the water table when deep -- TODO find a better default
      call fill_nest_halos_from_parent_masked("smcwtdxy", mn_phys%smcwtdxy, interp_type_lmask, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, x_refine, y_refine, &
          is_fine_pe, nest_domain, position, mn_phys%slmsk, mn_static%parent_ls%ls_mask_grid, M_LAND, 0.0D0)

      ! recharge to the water table when deep -- TODO find a better default
      call fill_nest_halos_from_parent_masked("deeprechxy", mn_phys%deeprechxy, interp_type_lmask, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, x_refine, y_refine, &
          is_fine_pe, nest_domain, position, mn_phys%slmsk, mn_static%parent_ls%ls_mask_grid, M_LAND, 0.0D0)
      ! recharge to the water table  -- TODO find a better default
      call fill_nest_halos_from_parent_masked("rechxy", mn_phys%rechxy, interp_type_lmask, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, x_refine, y_refine, &
          is_fine_pe, nest_domain, position, mn_phys%slmsk, mn_static%parent_ls%ls_mask_grid, M_LAND, 0.0D0)

      call fill_nest_halos_from_parent_masked("snicexy", mn_phys%snicexy, interp_type_lmask, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, x_refine, y_refine, &
          is_fine_pe, nest_domain, position, GFS_control%lsnow_lsm_lbound, GFS_control%lsnow_lsm_ubound, &
          mn_phys%slmsk, mn_static%parent_ls%ls_mask_grid, M_LAND, 0.0D0)

      call fill_nest_halos_from_parent_masked("snliqxy", mn_phys%snliqxy, interp_type_lmask, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, x_refine, y_refine, &
          is_fine_pe, nest_domain, position, GFS_control%lsnow_lsm_lbound,  GFS_control%lsnow_lsm_ubound, &
          mn_phys%slmsk, mn_static%parent_ls%ls_mask_grid, M_LAND, 0.0D0)

      ! surface snow thickness water equivalent over land - - default to 0
      call fill_nest_halos_from_parent_masked("snowd", mn_phys%snowd, interp_type_lmask, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, x_refine, y_refine, &
          is_fine_pe, nest_domain, position, mn_phys%slmsk, mn_static%parent_ls%ls_mask_grid, M_LAND, 0.0D0)

      ! Temperature in surface snow -- TODO notes say default to 0, but I will put 273.15K
      call fill_nest_halos_from_parent_masked("tsnoxy", mn_phys%tsnoxy, interp_type_lmask, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, x_refine, y_refine, &
          is_fine_pe, nest_domain, position, GFS_control%lsnow_lsm_lbound, GFS_control%lsnow_lsm_ubound, &
          !mn_phys%slmsk, mn_static%parent_ls%ls_mask_grid, M_LAND, 273.15D0)
          mn_phys%slmsk, mn_static%parent_ls%ls_mask_grid, M_LAND, 0.0D0)

      ! water equivalent accumulated snow depth over land - - default to 0
      call fill_nest_halos_from_parent_masked("weasd", mn_phys%weasd, interp_type_lmask, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, x_refine, y_refine, &
          is_fine_pe, nest_domain, position, mn_phys%slmsk, mn_static%parent_ls%ls_mask_grid, M_LAND, 0.0D0)

      call fill_nest_halos_from_parent_masked("smoiseq", mn_phys%smoiseq, interp_type_lmask, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, & 
          x_refine, y_refine, &
          is_fine_pe, nest_domain, position, 1, GFS_control%lsoil, mn_phys%slmsk, mn_static%parent_ls%ls_mask_grid, M_LAND, 0.3D0)


      call fill_nest_halos_from_parent_masked("zsnsoxy", mn_phys%zsnsoxy, interp_type_lmask, Atm(child_grid_num)%neststruct%wt_h, &
          Atm(child_grid_num)%neststruct%ind_h, x_refine, y_refine, &
          is_fine_pe, nest_domain, position, GFS_control%lsnow_lsm_lbound, GFS_control%lsoil, &
          mn_phys%slmsk, mn_static%parent_ls%ls_mask_grid, M_LAND, zsns_default)

    endif

  end subroutine mn_phys_fill_nest_halos_from_parent

  !>@brief The subroutine 'mn_phys_fill_intern_nest_halos' fills the intenal nest halos for the physics variables
  !>@details This subroutine is only called for the nest PEs.
  subroutine mn_phys_fill_intern_nest_halos(moving_nest, GFS_control, domain_fine, is_fine_pe)
    type(fv_moving_nest_type), target, intent(inout) :: moving_nest         !< Single instance of moving nest data
    type(GFS_control_type), intent(in)               :: GFS_control         !< Physics metadata
    type(domain2d), intent(inout)                    :: domain_fine         !< Domain structure for this nest
    logical, intent(in)                              :: is_fine_pe          !< Is nest PE - should be True.  Argument is redundant.

    type(fv_moving_nest_physics_type), pointer :: mn_phys

    mn_phys => moving_nest%mn_phys

    call mn_var_fill_intern_nest_halos(mn_phys%ts, domain_fine, is_fine_pe)   !! Skin Temp/SST
    if (move_physics) then
      call mn_var_fill_intern_nest_halos(mn_phys%smc, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%stc, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%slc, domain_fine, is_fine_pe)

      call mn_var_fill_intern_nest_halos(mn_phys%phy_f2d, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%phy_f3d, domain_fine, is_fine_pe)

      call mn_var_fill_intern_nest_halos(mn_phys%emis_lnd, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%emis_ice, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%emis_wat, domain_fine, is_fine_pe)

      !call mn_var_fill_intern_nest_halos(mn_phys%sfalb_lnd, domain_fine, is_fine_pe)
      !call mn_var_fill_intern_nest_halos(mn_phys%sfalb_lnd_bck, domain_fine, is_fine_pe)
      !call mn_var_fill_intern_nest_halos(mn_phys%semis, domain_fine, is_fine_pe)
      !call mn_var_fill_intern_nest_halos(mn_phys%semisbase, domain_fine, is_fine_pe)
      !call mn_var_fill_intern_nest_halos(mn_phys%sfalb, domain_fine, is_fine_pe)

      call mn_var_fill_intern_nest_halos(mn_phys%u10m, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%v10m, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%tprcp, domain_fine, is_fine_pe)

      call mn_var_fill_intern_nest_halos(mn_phys%hprime, domain_fine, is_fine_pe)

      call mn_var_fill_intern_nest_halos(mn_phys%lakefrac, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%lakedepth, domain_fine, is_fine_pe)

      call mn_var_fill_intern_nest_halos(mn_phys%canopy, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%vegfrac, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%uustar, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%shdmin, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%shdmax, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%zorl, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%zorll, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%zorlwav, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%zorlw, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%usfco, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%vsfco, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%tsfco, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%tsfcl, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%tsfc, domain_fine, is_fine_pe)

      call mn_var_fill_intern_nest_halos(mn_phys%albdirvis_lnd, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%albdirnir_lnd, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%albdifvis_lnd, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%albdifnir_lnd, domain_fine, is_fine_pe)

      call mn_var_fill_intern_nest_halos(mn_phys%cv, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%cvt, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%cvb, domain_fine, is_fine_pe)
    endif

    if (move_nsst) then
      call mn_var_fill_intern_nest_halos(mn_phys%tref, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%z_c, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%c_0, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%c_d, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%w_0, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%w_d, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%xt, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%xs, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%xu, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%xv, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%xz, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%zm, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%xtts, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%xzts, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%d_conv, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%dt_cool, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%qrain, domain_fine, is_fine_pe)
    endif

    if (move_physics .and. GFS_control%lsm == GFS_control%lsm_noahmp) then
      call mn_var_fill_intern_nest_halos(mn_phys%soilcolor, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%snowxy, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%tvxy, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%tgxy, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%canicexy, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%canliqxy, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%eahxy, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%tahxy, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%cmxy, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%chxy, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%fwetxy, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%sneqvoxy, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%alboldxy, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%qsnowxy, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%wslakexy, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%zwtxy, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%waxy, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%wtxy, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%lfmassxy, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%rtmassxy, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%stmassxy, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%woodxy, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%stblcpxy, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%fastcpxy, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%xsaixy, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%xlaixy, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%taussxy, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%smcwtdxy, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%deeprechxy, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%rechxy, domain_fine, is_fine_pe)

      call mn_var_fill_intern_nest_halos(mn_phys%snicexy, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%snliqxy, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%snowd, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%tsnoxy, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%weasd, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%smoiseq, domain_fine, is_fine_pe)
      call mn_var_fill_intern_nest_halos(mn_phys%zsnsoxy, domain_fine, is_fine_pe)

    endif

  end subroutine mn_phys_fill_intern_nest_halos

  !>@brief The subroutine 'mn_phys_shift_data' shifts the variable in the nest, including interpolating at the leading edge
  !>@details This subroutine is called for the nest and parent PEs.
  subroutine mn_phys_shift_data(Atm, GFS_control, n, child_grid_num, wt_h, wt_u, wt_v, &
      delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, nz)
    type(fv_atmos_type), allocatable, target, intent(inout)  :: Atm(:)                  !< Array of atmospheric data
    type(GFS_control_type), intent(in)                       :: GFS_control             !< Physics metadata
    integer, intent(in)                                      :: n, child_grid_num       !< Current grid number, child grid number
    real, allocatable, intent(in)                            :: wt_h(:,:,:), wt_u(:,:,:), wt_v(:,:,:) !< Interpolation weights
    integer, intent(in)                                      :: delta_i_c, delta_j_c    !< Nest motion in i,j direction
    integer, intent(in)                                      :: x_refine, y_refine      !< Nest refinement
    logical, intent(in)                                      :: is_fine_pe              !< Is this the nest PE?
    type(nest_domain_type), intent(inout)                    :: nest_domain             !< Nest domain structure
    integer, intent(in)                                      :: nz                      !< Number of vertical levels

    ! Constants for mpp calls
    integer  :: interp_type   = 1    ! cell-centered A-grid
    integer  :: interp_type_u = 4    ! D-grid
    integer  :: interp_type_v = 4    ! D-grid
    integer  :: position      = CENTER
    integer  :: position_u    = NORTH
    integer  :: position_v    = EAST
    type(fv_moving_nest_physics_type), pointer :: mn_phys

    mn_phys => Moving_nest(n)%mn_phys

    !! Skin temp/SST
    call mn_var_shift_data(mn_phys%ts, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
        delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)

    if (move_physics) then
      !! Soil variables
      call mn_var_shift_data(mn_phys%smc, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position, GFS_control%lsoil)
      call mn_var_shift_data(mn_phys%stc, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position, GFS_control%lsoil)
      call mn_var_shift_data(mn_phys%slc, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position, GFS_control%lsoil)

      !! Physics arrays
      call mn_var_shift_data(mn_phys%phy_f2d, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position, GFS_control%ntot2d)

      call mn_var_shift_data(mn_phys%phy_f3d, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position, GFS_control%levs)

      ! Surface variables

      call mn_var_shift_data(mn_phys%emis_lnd, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%emis_ice, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%emis_wat, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)


      !call mn_var_shift_data(mn_phys%sfalb_lnd, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
      !  delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      !call mn_var_shift_data(mn_phys%sfalb_lnd_bck, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
      !  delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      !call mn_var_shift_data(mn_phys%semis, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
      !  delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      !call mn_var_shift_data(mn_phys%semisbase, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
      !  delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      !call mn_var_shift_data(mn_phys%sfalb, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
      !  delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)

      call mn_var_shift_data(mn_phys%u10m, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%v10m, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%tprcp, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%hprime, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position, GFS_control%nmtvr)
      call mn_var_shift_data(mn_phys%lakefrac, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%lakedepth, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%canopy, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%vegfrac, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%uustar, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%shdmin, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%shdmax, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%zorl, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%zorll, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%zorlwav, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%zorlw, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%usfco, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%vsfco, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%tsfco, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%tsfcl, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%tsfc, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%albdirvis_lnd, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%albdirnir_lnd, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%albdifvis_lnd, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%albdifnir_lnd, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%cv, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%cvt, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%cvb, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
    endif

    if (move_nsst) then
      call mn_var_shift_data(mn_phys%tref, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%z_c, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%c_0, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%c_d, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%w_0, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%w_d, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%xt, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%xs, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%xu, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%xv, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%xz, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%zm, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%xtts, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%xzts, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%d_conv, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%dt_cool, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%qrain, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
    endif

    if (move_physics .and. GFS_control%lsm == GFS_control%lsm_noahmp) then
      call mn_var_shift_data(mn_phys%soilcolor, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%snowxy, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%tvxy, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%tgxy, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%canicexy, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%canliqxy, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%eahxy, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%tahxy, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%cmxy, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%chxy, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%fwetxy, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%sneqvoxy, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%alboldxy, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%qsnowxy, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%wslakexy, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%zwtxy, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%waxy, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%wtxy, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%lfmassxy, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%rtmassxy, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%stmassxy, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%woodxy, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%stblcpxy, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%fastcpxy, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%xsaixy, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%xlaixy, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%taussxy, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%smcwtdxy, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%deeprechxy, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%rechxy, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%smoiseq, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position, GFS_control%lsoil)


      call mn_var_shift_data(mn_phys%snicexy, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position, GFS_control%lsnow_lsm_lbound, GFS_control%lsnow_lsm_ubound)
      call mn_var_shift_data(mn_phys%snliqxy, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position, GFS_control%lsnow_lsm_lbound, GFS_control%lsnow_lsm_ubound)
      call mn_var_shift_data(mn_phys%snowd, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%tsnoxy, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position, GFS_control%lsnow_lsm_lbound, GFS_control%lsnow_lsm_ubound)
      call mn_var_shift_data(mn_phys%weasd, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position)
      call mn_var_shift_data(mn_phys%zsnsoxy, interp_type, wt_h, Atm(child_grid_num)%neststruct%ind_h, &
          delta_i_c, delta_j_c, x_refine, y_refine, is_fine_pe, nest_domain, position, GFS_control%lsnow_lsm_lbound, GFS_control%lsoil)

    endif


  end subroutine mn_phys_shift_data

  !>@brief The subroutine 'mn_phys_dump_to_netcdf' dumps physics variables to debugging netCDF files
  !>@details This subroutine is called for the nest and parent PEs.
  subroutine mn_phys_dump_to_netcdf(Atm, Atm_block, GFS_control, GFS_sfcprop, GFS_tbd, time_val, file_prefix, is_fine_pe, domain_coarse, domain_fine, nz)
    type(fv_atmos_type), intent(in)            :: Atm                           !< Single instance of atmospheric data
    type (block_control_type), intent(in)      :: Atm_block                     !< Physics block layout
    type(GFS_control_type), intent(in)         :: GFS_control                   !< Physics metadata
    type(GFS_sfcprop_type), intent(in)         :: GFS_sfcprop                   !< Physics variable data (surface)
    type(GFS_tbd_type), intent(in)             :: GFS_tbd                       !< Physics variable data (tbd)
    integer, intent(in)                        :: time_val                      !< Timestep number for filename
    character(len=*), intent(in)               :: file_prefix                   !< Prefix for output netCDF filenames
    logical, intent(in)                        :: is_fine_pe                    !< Is this the nest PE?
    type(domain2d), intent(in)                 :: domain_coarse, domain_fine    !< Domain structures for parent and nest
    integer, intent(in)                        :: nz                            !< Number of vertical levels

    integer :: is, ie, js, je
    integer :: nb, blen, i, j, k, ix, nv, im
    integer :: this_pe

    integer            :: n_moist
    character(len=16)  :: out_var_name, phys_var_name
    integer            :: position = CENTER

    ! Coerce the double precision variables from physics into single precision for debugging netCDF output
    ! Does not affect values used in calculations.
    ! TODO do we want to dump these as double precision??
    real, allocatable :: smc_pr_local (:,:,:)  !< soil moisture content
    real, allocatable :: stc_pr_local (:,:,:)  !< soil temperature
    real, allocatable :: slc_pr_local (:,:,:)  !< soil liquid water content
    real, allocatable, dimension(:,:) :: sealand_pr_local, deep_soil_t_pr_local, soil_type_pr_local, veg_type_pr_local, slope_type_pr_local, max_snow_alb_pr_local
    real, allocatable, dimension(:,:) :: tsfco_pr_local, tsfcl_pr_local, tsfc_pr_local, vegfrac_pr_local
    real, allocatable, dimension(:,:) :: tref_pr_local, c_0_pr_local, xt_pr_local,  xu_pr_local,  xv_pr_local, ifd_pr_local
    real, allocatable, dimension(:,:) :: facsf_pr_local, facwf_pr_local
    real, allocatable, dimension(:,:) :: alvsf_pr_local, alvwf_pr_local, alnsf_pr_local, alnwf_pr_local
    real, allocatable, dimension(:,:) :: zorl_pr_local, zorll_pr_local, zorlw_pr_local, zorli_pr_local
    real, allocatable, dimension(:,:) :: usfco_pr_local, vsfco_pr_local
    real, allocatable :: phy_f2d_pr_local (:,:,:)
    real, allocatable :: phy_f3d_pr_local (:,:,:,:)
    real, allocatable :: lakefrac_pr_local (:,:)  !< lake fraction
    real, allocatable :: landfrac_pr_local (:,:)  !< land fraction
    real, allocatable :: emis_lnd_pr_local (:,:)  !< emissivity land
    real, allocatable :: snowxy_pr_local (:,:)     !< number of snow layers

    logical :: move_noahmp
    move_noahmp = .True.

    this_pe = mpp_pe()

    !  Skin temp/SST
    call mn_var_dump_to_netcdf(Atm%ts, is_fine_pe, domain_coarse, domain_fine, position, time_val, Atm%global_tile, file_prefix, "SSTK")
    !  Terrain height == phis / grav
    call mn_var_dump_to_netcdf(Atm%phis / grav, is_fine_pe, domain_coarse, domain_fine, position, time_val, Atm%global_tile, file_prefix, "orog")

    ! sgh and oro were only fully allocated if fv_land is True
    !      if false, oro is (1,1), and sgh is not allocated
    if ( Atm%flagstruct%fv_land ) then
      ! land frac --  called oro in fv_array.F90
      call mn_var_dump_to_netcdf(Atm%oro, is_fine_pe, domain_coarse, domain_fine, position, time_val, Atm%global_tile, file_prefix, "LFRAC")
      ! terrain standard deviation --  called sgh in fv_array.F90
      call mn_var_dump_to_netcdf(Atm%sgh, is_fine_pe, domain_coarse, domain_fine, position, time_val, Atm%global_tile, file_prefix, "STDDEV")
    endif

    is = Atm%bd%is
    ie = Atm%bd%ie
    js = Atm%bd%js
    je = Atm%bd%je

    ! Just allocate compute domain size here for outputs;  the nest moving code also has halos added, but we don't need them here.
    if (move_physics) then
      allocate ( smc_pr_local(is:ie, js:je, GFS_control%lsoil) )
      allocate ( stc_pr_local(is:ie, js:je, GFS_control%lsoil) )
      allocate ( slc_pr_local(is:ie, js:je, GFS_control%lsoil) )
      allocate ( sealand_pr_local(is:ie, js:je) )
      allocate ( lakefrac_pr_local(is:ie, js:je) )
      allocate ( landfrac_pr_local(is:ie, js:je) )
      allocate ( emis_lnd_pr_local(is:ie, js:je) )
      allocate ( phy_f2d_pr_local(is:ie, js:je, GFS_control%ntot2d) )
      allocate ( phy_f3d_pr_local(is:ie, js:je, GFS_control%levs, GFS_control%ntot3d) )
      allocate ( tsfco_pr_local(is:ie, js:je) )
      allocate ( tsfcl_pr_local(is:ie, js:je) )
      allocate ( tsfc_pr_local(is:ie, js:je) )
      allocate ( vegfrac_pr_local(is:ie, js:je) )
      allocate ( alvsf_pr_local(is:ie, js:je) )
      allocate ( alvwf_pr_local(is:ie, js:je) )
      allocate ( alnsf_pr_local(is:ie, js:je) )
      allocate ( alnwf_pr_local(is:ie, js:je) )
      allocate ( deep_soil_t_pr_local(is:ie, js:je) )
      allocate ( soil_type_pr_local(is:ie, js:je) )
      !allocate ( veg_frac_pr_local(is:ie, js:je) )
      allocate ( veg_type_pr_local(is:ie, js:je) )
      allocate ( slope_type_pr_local(is:ie, js:je) )
      allocate ( max_snow_alb_pr_local(is:ie, js:je) )
      allocate ( facsf_pr_local(is:ie, js:je) )
      allocate ( facwf_pr_local(is:ie, js:je) )
      allocate ( zorl_pr_local(is:ie, js:je) )
      allocate ( zorll_pr_local(is:ie, js:je) )
      allocate ( zorlw_pr_local(is:ie, js:je) )
      allocate ( zorli_pr_local(is:ie, js:je) )
      allocate ( usfco_pr_local(is:ie, js:je) )
      allocate ( vsfco_pr_local(is:ie, js:je) )
    endif

    if (move_nsst) then
      allocate ( tref_pr_local(is:ie, js:je) )
      allocate ( c_0_pr_local(is:ie, js:je) )
      allocate ( xt_pr_local(is:ie, js:je) )
      allocate ( xu_pr_local(is:ie, js:je) )
      allocate ( xv_pr_local(is:ie, js:je) )
      allocate ( ifd_pr_local(is:ie, js:je) )
    endif

    if (move_noahmp) then
      allocate ( snowxy_pr_local(is:ie, js:je) )
    endif

    if (move_physics) then
      smc_pr_local = +99999.9
      stc_pr_local = +99999.9
      slc_pr_local = +99999.9
      sealand_pr_local = +99999.9
      lakefrac_pr_local = +99999.9
      landfrac_pr_local = +99999.9
      emis_lnd_pr_local = +99999.9
      phy_f2d_pr_local = +99999.9
      phy_f3d_pr_local = +99999.9
      tsfco_pr_local = +99999.9
      tsfcl_pr_local = +99999.9
      tsfc_pr_local = +99999.9
      vegfrac_pr_local = +99999.9
      alvsf_pr_local = +99999.9
      alvwf_pr_local = +99999.9
      alnsf_pr_local = +99999.9
      alnwf_pr_local = +99999.9
    endif
    if (move_nsst) then
      tref_pr_local = +99999.9
      c_0_pr_local = +99999.9
      xt_pr_local = +99999.9
      xu_pr_local = +99999.9
      xv_pr_local = +99999.9
      ifd_pr_local = +99999.9
    endif
    if (move_nsst) then
      snowxy_pr_local = +99999.9
    endif

    im = 0
    do nb = 1,Atm_block%nblks
      blen = Atm_block%blksz(nb)
      do ix = 1, blen
        i = Atm_block%index(nb)%ii(ix)
        j = Atm_block%index(nb)%jj(ix)
        im = im + 1

        if (move_physics) then
          do k = 1, GFS_control%lsoil
            ! Use real() to lower the precision
            smc_pr_local(i,j,k) = real(GFS_sfcprop%smc(im,k))
            stc_pr_local(i,j,k) = real(GFS_sfcprop%stc(im,k))
            slc_pr_local(i,j,k) = real(GFS_sfcprop%slc(im,k))
          enddo

          sealand_pr_local(i,j) = real(GFS_sfcprop%slmsk(im))
          lakefrac_pr_local(i,j) = real(GFS_sfcprop%lakefrac(im))
          landfrac_pr_local(i,j) = real(GFS_sfcprop%landfrac(im))
          emis_lnd_pr_local(i,j) = real(GFS_sfcprop%emis_lnd(im))
          deep_soil_t_pr_local(i, j) = GFS_sfcprop%tg3(im)
          soil_type_pr_local(i, j) = GFS_sfcprop%stype(im)
          !veg_frac_pr_local(i, j) = GFS_sfcprop%vfrac(im)
          veg_type_pr_local(i, j) = GFS_sfcprop%vtype(im)
          slope_type_pr_local(i, j) = GFS_sfcprop%slope(im)
          facsf_pr_local(i, j) = GFS_sfcprop%facsf(im)
          facwf_pr_local(i, j) = GFS_sfcprop%facwf(im)
          zorl_pr_local(i, j) = GFS_sfcprop%zorl(im)
          zorlw_pr_local(i, j) = GFS_sfcprop%zorlw(im)
          zorll_pr_local(i, j) = GFS_sfcprop%zorll(im)
          zorli_pr_local(i, j) = GFS_sfcprop%zorli(im)
          usfco_pr_local(i, j) = GFS_sfcprop%usfco(im)
          vsfco_pr_local(i, j) = GFS_sfcprop%vsfco(im)
          max_snow_alb_pr_local(i, j) = GFS_sfcprop%snoalb(im)
          tsfco_pr_local(i, j) = GFS_sfcprop%tsfco(im)
          tsfcl_pr_local(i, j) = GFS_sfcprop%tsfcl(im)
          tsfc_pr_local(i, j)  = GFS_sfcprop%tsfc(im)
          vegfrac_pr_local(i, j) = GFS_sfcprop%vfrac(im)
          alvsf_pr_local(i, j) = GFS_sfcprop%alvsf(im)
          alvwf_pr_local(i, j) = GFS_sfcprop%alvwf(im)
          alnsf_pr_local(i, j) = GFS_sfcprop%alnsf(im)
          alnwf_pr_local(i, j) = GFS_sfcprop%alnwf(im)

          do nv = 1, GFS_Control%ntot2d
            ! Use real() to lower the precision
            phy_f2d_pr_local(i,j,nv) = real(GFS_tbd%phy_f2d(im, nv))
          enddo

          do k = 1, GFS_control%levs
            do nv = 1, GFS_control%ntot3d
              ! Use real() to lower the precision
              phy_f3d_pr_local(i,j,k,nv) = real(GFS_tbd%phy_f3d(im, k, nv))
            enddo
          enddo
        endif

        if (move_nsst) then
          tref_pr_local(i,j) = GFS_sfcprop%tref(im)
          c_0_pr_local(i,j) = GFS_sfcprop%c_0(im)
          xt_pr_local(i,j) = GFS_sfcprop%xt(im)
          xu_pr_local(i,j) = GFS_sfcprop%xu(im)
          xv_pr_local(i,j) = GFS_sfcprop%xv(im)
          ifd_pr_local(i,j) = GFS_sfcprop%ifd(im)
        endif

        if (move_noahmp) then
          snowxy_pr_local(i,j) = GFS_sfcprop%snowxy(im)
        endif

      enddo
    enddo

    if (move_physics) then
      !call mn_var_dump_to_netcdf(stc_pr_local, is_fine_pe, domain_coarse, domain_fine, position, GFS_control%lsoil, time_val, Atm%global_tile, file_prefix, "SOILT")
      !call mn_var_dump_to_netcdf(smc_pr_local, is_fine_pe, domain_coarse, domain_fine, position, GFS_control%lsoil, time_val, Atm%global_tile, file_prefix, "SOILM")
      !call mn_var_dump_to_netcdf(slc_pr_local, is_fine_pe, domain_coarse, domain_fine, position, GFS_control%lsoil, time_val, Atm%global_tile, file_prefix, "SOILL")
      call mn_var_dump_to_netcdf(sealand_pr_local, is_fine_pe, domain_coarse, domain_fine, position, time_val, Atm%global_tile, file_prefix, "LMASK")
      call mn_var_dump_to_netcdf(lakefrac_pr_local, is_fine_pe, domain_coarse, domain_fine, position, time_val, Atm%global_tile, file_prefix, "LAKEFRAC")
      call mn_var_dump_to_netcdf(landfrac_pr_local, is_fine_pe, domain_coarse, domain_fine, position, time_val, Atm%global_tile, file_prefix, "LANDFRAC")
      call mn_var_dump_to_netcdf(emis_lnd_pr_local, is_fine_pe, domain_coarse, domain_fine, position, time_val, Atm%global_tile, file_prefix, "EMISLAND")
      call mn_var_dump_to_netcdf(deep_soil_t_pr_local, is_fine_pe, domain_coarse, domain_fine, position, time_val, Atm%global_tile, file_prefix, "DEEPSOIL")
      call mn_var_dump_to_netcdf(soil_type_pr_local, is_fine_pe, domain_coarse, domain_fine, position, time_val, Atm%global_tile, file_prefix, "SOILTP")
      !call mn_var_dump_to_netcdf(veg_frac_pr_local, is_fine_pe, domain_coarse, domain_fine, position, time_val, Atm%global_tile, file_prefix, "VEGFRAC")
      call mn_var_dump_to_netcdf(veg_type_pr_local, is_fine_pe, domain_coarse, domain_fine, position, time_val, Atm%global_tile, file_prefix, "VEGTYPE")
      call mn_var_dump_to_netcdf(slope_type_pr_local, is_fine_pe, domain_coarse, domain_fine, position, time_val, Atm%global_tile, file_prefix, "SLOPE")
      call mn_var_dump_to_netcdf(max_snow_alb_pr_local, is_fine_pe, domain_coarse, domain_fine, position, time_val, Atm%global_tile, file_prefix, "SNOWALB")
      call mn_var_dump_to_netcdf(tsfco_pr_local, is_fine_pe, domain_coarse, domain_fine, position, time_val, Atm%global_tile, file_prefix, "TSFCO")
      call mn_var_dump_to_netcdf(tsfcl_pr_local, is_fine_pe, domain_coarse, domain_fine, position, time_val, Atm%global_tile, file_prefix, "TSFCL")
      call mn_var_dump_to_netcdf(tsfc_pr_local, is_fine_pe, domain_coarse, domain_fine, position, time_val, Atm%global_tile, file_prefix, "TSFC")
      call mn_var_dump_to_netcdf(vegfrac_pr_local, is_fine_pe, domain_coarse, domain_fine, position, time_val, Atm%global_tile, file_prefix, "VEGFRAC")
      call mn_var_dump_to_netcdf(alvsf_pr_local, is_fine_pe, domain_coarse, domain_fine, position, time_val, Atm%global_tile, file_prefix, "ALVSF")
      call mn_var_dump_to_netcdf(alvwf_pr_local, is_fine_pe, domain_coarse, domain_fine, position, time_val, Atm%global_tile, file_prefix, "ALVWF")
      call mn_var_dump_to_netcdf(alnsf_pr_local, is_fine_pe, domain_coarse, domain_fine, position, time_val, Atm%global_tile, file_prefix, "ALNSF")
      call mn_var_dump_to_netcdf(alnwf_pr_local, is_fine_pe, domain_coarse, domain_fine, position, time_val, Atm%global_tile, file_prefix, "ALNWF")
      call mn_var_dump_to_netcdf(facsf_pr_local, is_fine_pe, domain_coarse, domain_fine, position, time_val, Atm%global_tile, file_prefix, "FACSF")
      call mn_var_dump_to_netcdf(facwf_pr_local, is_fine_pe, domain_coarse, domain_fine, position, time_val, Atm%global_tile, file_prefix, "FACWF")
      call mn_var_dump_to_netcdf(zorl_pr_local, is_fine_pe, domain_coarse, domain_fine, position, time_val, Atm%global_tile, file_prefix, "ZORL")
      call mn_var_dump_to_netcdf(zorlw_pr_local, is_fine_pe, domain_coarse, domain_fine, position, time_val, Atm%global_tile, file_prefix, "ZORLW")
      call mn_var_dump_to_netcdf(zorll_pr_local, is_fine_pe, domain_coarse, domain_fine, position, time_val, Atm%global_tile, file_prefix, "ZORLL")
      call mn_var_dump_to_netcdf(zorli_pr_local, is_fine_pe, domain_coarse, domain_fine, position, time_val, Atm%global_tile, file_prefix, "ZORLI")
      call mn_var_dump_to_netcdf(usfco_pr_local, is_fine_pe, domain_coarse, domain_fine, position, time_val, Atm%global_tile, file_prefix, "SSU")
      call mn_var_dump_to_netcdf(vsfco_pr_local, is_fine_pe, domain_coarse, domain_fine, position, time_val, Atm%global_tile, file_prefix, "SSV")

      do nv = 1, GFS_control%ntot2d
        write (phys_var_name, "(A4,I0.3)")  'PH2D', nv
        !call mn_var_dump_to_netcdf(phy_f2d_pr_local(:,:,nv), is_fine_pe, domain_coarse, domain_fine, position, 1, &
        !    time_val, Atm%global_tile, file_prefix, phys_var_name)
      enddo

      do nv = 1, GFS_control%ntot3d
        write (phys_var_name, "(A4,I0.3)")  'PH3D', nv
        !call mn_var_dump_to_netcdf(phy_f3d_pr_local(:,:,:,nv), is_fine_pe, domain_coarse, domain_fine, position, GFS_control%levs, &
        !    time_val, Atm%global_tile, file_prefix, phys_var_name)
      enddo
    endif

    if (move_nsst) then
      call mn_var_dump_to_netcdf(tref_pr_local, is_fine_pe, domain_coarse, domain_fine, position, time_val, Atm%global_tile, file_prefix, "TREF")
      call mn_var_dump_to_netcdf(c_0_pr_local, is_fine_pe, domain_coarse, domain_fine, position, time_val, Atm%global_tile, file_prefix, "C_0")
      call mn_var_dump_to_netcdf(xt_pr_local, is_fine_pe, domain_coarse, domain_fine, position, time_val, Atm%global_tile, file_prefix, "XT")
      call mn_var_dump_to_netcdf(xu_pr_local, is_fine_pe, domain_coarse, domain_fine, position, time_val, Atm%global_tile, file_prefix, "XU")
      call mn_var_dump_to_netcdf(xv_pr_local, is_fine_pe, domain_coarse, domain_fine, position, time_val, Atm%global_tile, file_prefix, "XV")
      call mn_var_dump_to_netcdf(ifd_pr_local, is_fine_pe, domain_coarse, domain_fine, position, time_val, Atm%global_tile, file_prefix, "IFD")
    endif

    if (move_noahmp) then
      call mn_var_dump_to_netcdf(snowxy_pr_local, is_fine_pe, domain_coarse, domain_fine, position, time_val, Atm%global_tile, file_prefix, "SNOWXY")
    endif

    if (move_physics) then
      deallocate(smc_pr_local)
      deallocate(stc_pr_local)
      deallocate(slc_pr_local)
      deallocate(lakefrac_pr_local)
      deallocate(landfrac_pr_local)
      deallocate(emis_lnd_pr_local)
      deallocate(sealand_pr_local, deep_soil_t_pr_local, soil_type_pr_local, veg_type_pr_local, max_snow_alb_pr_local)
      deallocate(tsfco_pr_local, tsfcl_pr_local, tsfc_pr_local, vegfrac_pr_local)
      deallocate(alvsf_pr_local, alvwf_pr_local, alnsf_pr_local, alnwf_pr_local)
      deallocate(facsf_pr_local, facwf_pr_local)
      deallocate(zorl_pr_local, zorlw_pr_local, zorll_pr_local, zorli_pr_local)
      deallocate(usfco_pr_local, vsfco_pr_local)
      deallocate(phy_f2d_pr_local)
      deallocate(phy_f3d_pr_local)
    endif

    if (move_nsst) deallocate(tref_pr_local, c_0_pr_local, xt_pr_local,  xu_pr_local,  xv_pr_local, ifd_pr_local)

    if (move_noahmp) deallocate(snowxy_pr_local)


  end subroutine mn_phys_dump_to_netcdf

end module fv_moving_nest_physics_mod
