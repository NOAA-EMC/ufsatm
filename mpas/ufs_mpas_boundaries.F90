!> ###########################################################################################
!> \file ufs_mpas_boundaries.F90
!>
!> Routines adopted from MPAS src/core_atmosphere/dynamics/mpas_atm_boundaries.F for use in
!> the UFS Weather Model.
!>
!> ###########################################################################################
module ufs_mpas_boundaries
  use mpas_atm_boundaries, only : LBC_intv_end
  use ufs_mpas_io

  implicit none

  public :: ufs_mpas_atm_update_bdy_tend, ufs_mpas_atm_bdy_checks
  
contains

  !> #########################################################################################
  !>
  !> routine ufs_mpas_atm_update_bdy_tend
  !>
  !> \brief   Reads new boundary data and updates the LBC tendencies
  !> \author  Michael Duda
  !> \date    27 September 2016
  !> \details
  !>  This routine reads from the 'lbc_in' stream all variables in the 'lbc'
  !>  pool. When called with firstCall=.true., the latest time before the
  !>  present is read into time level 2 of the lbc pool; otherwise, the
  !>  contents of time level 2 are shifted to time level 1, the earliest
  !>  time strictly later than the present is read into time level 2, and
  !>  the tendencies for all fields in the lbc pool are computed and stored
  !>  in time level 1.
  !>
  !> \update: Dustin Swales September 2025 - Modified for use in UWM
  !>
  !> #########################################################################################
  subroutine ufs_mpas_atm_update_bdy_tend(clock, block, firstCall, nRecord, ierr, debug)
    use mpas_constants,      only : rvord
    use mpas_log,            only : mpas_log_write
    use mpas_derived_types,  only : MPAS_STREAM_MGR_NOERR, MPAS_LOG_ERR
    use mpas_derived_types,  only : mpas_pool_type, mpas_Clock_type, block_type
    use mpas_derived_types,  only : MPAS_TimeInterval_type, MPAS_Time_Type
    use mpas_timekeeping,    only : mpas_set_time
    use mpas_kind_types,     only : StrKIND, RKIND
    use mpas_derived_types,  only : MPAS_STREAM_LATEST_BEFORE
    use mpas_derived_types,  only : MPAS_STREAM_EARLIEST_STRICTLY_AFTER
    use mpas_timekeeping,    only : mpas_get_timeInterval, mpas_get_time, operator(-)
    use mpas_timekeeping,    only : mpas_get_clock_time, MPAS_NOW
    use mpas_pool_routines,  only : mpas_pool_get_config, mpas_pool_get_subpool
    use mpas_pool_routines,  only : mpas_pool_shift_time_levels, mpas_pool_get_array
    use mpas_pool_routines,  only : mpas_pool_get_dimension
    use module_mpas_config,  only : lbc_filename, pioid_lbc, pio_subsystem_lbc

    implicit none

    type (mpas_clock_type), intent(in) :: clock
    type (block_type), intent(inout) :: block
    logical, intent(in) :: firstCall
    integer, intent(in) :: nRecord
    integer, intent(out) :: ierr
    logical, intent(in) :: debug

    character(len=StrKIND) :: lbc_intv_start_string
    character(len=StrKIND) :: lbc_intv_end_string

    type (mpas_pool_type), pointer :: mesh
    type (mpas_pool_type), pointer :: state
    type (mpas_pool_type), pointer :: lbc
    real (kind=RKIND) :: dt

    integer, pointer :: nCells_ptr
    integer, pointer :: nEdges_ptr
    integer, pointer :: nVertLevels_ptr
    integer, pointer :: index_qv_ptr
    integer, pointer :: nScalars_ptr
    integer :: nCells, nEdges, nVertLevels, index_qv, nScalars

    real (kind=RKIND), dimension(:,:), pointer :: u
    real (kind=RKIND), dimension(:,:), pointer :: ru
    real (kind=RKIND), dimension(:,:), pointer :: rho_edge
    real (kind=RKIND), dimension(:,:), pointer :: w
    real (kind=RKIND), dimension(:,:), pointer :: theta
    real (kind=RKIND), dimension(:,:), pointer :: rtheta_m
    real (kind=RKIND), dimension(:,:), pointer :: rho_zz
    real (kind=RKIND), dimension(:,:), pointer :: rho
    real (kind=RKIND), dimension(:,:,:), pointer :: scalars
    real (kind=RKIND), dimension(:,:), pointer :: lbc_tend_u
    real (kind=RKIND), dimension(:,:), pointer :: lbc_tend_ru
    real (kind=RKIND), dimension(:,:), pointer :: lbc_tend_rho_edge
    real (kind=RKIND), dimension(:,:), pointer :: lbc_tend_w
    real (kind=RKIND), dimension(:,:), pointer :: lbc_tend_theta
    real (kind=RKIND), dimension(:,:), pointer :: lbc_tend_rtheta_m
    real (kind=RKIND), dimension(:,:), pointer :: lbc_tend_rho_zz
    real (kind=RKIND), dimension(:,:), pointer :: lbc_tend_rho
    real (kind=RKIND), dimension(:,:,:), pointer :: lbc_tend_scalars

    integer, dimension(:,:), pointer :: cellsOnEdge
    real (kind=RKIND), dimension(:,:), pointer :: zz

    integer :: dd_intv, s_intv, sn_intv, sd_intv
    type (MPAS_Time_Type) :: currTime
    type (MPAS_TimeInterval_Type) :: lbc_interval
    character(len=StrKIND) :: read_time
    integer :: iEdge, iCell, k, j
    integer :: cell1, cell2


    ierr = 0

    call mpas_pool_get_subpool(block % structs, 'mesh', mesh)
    call mpas_pool_get_subpool(block % structs, 'state', state)
    call mpas_pool_get_subpool(block % structs, 'lbc', lbc)

    if (firstCall) then
       call dyn_mpas_read_write_stream(clock, 'r', 'lbc_in', pio_file_desc=pioid_lbc, ierr=ierr, timeLevel=2, &
            whence = MPAS_STREAM_LATEST_BEFORE, actualWhen=read_time, nRecord=nRecord, debug=debug)
       if (ierr /= MPAS_STREAM_MGR_NOERR) then
          call mpas_log_write('Could not read from ''lbc_in'' stream on or before the current date '// &
                              'to update lateral boundary tendencies', messageType=MPAS_LOG_ERR)
          ierr = 1
       end if
    else
       call mpas_pool_shift_time_levels(lbc)
       call dyn_mpas_read_write_stream(clock, 'r', 'lbc_in', pio_file_desc=pioid_lbc, ierr=ierr, timeLevel=2,  &
            whence = MPAS_STREAM_EARLIEST_STRICTLY_AFTER, actualWhen=read_time, nRecord=nRecord, debug=debug)
       if (ierr /= MPAS_STREAM_MGR_NOERR) then
          call mpas_log_write('Could not read from ''lbc_in'' stream after the current date '// &
                              'to update lateral boundary tendencies', messageType=MPAS_LOG_ERR)
          ierr = 1
       end if
    end if
    if (ierr /= 0) then
       return
    end if

    call mpas_set_time(currTime, dateTimeString=trim(read_time))
    
    !
    ! Compute any derived fields from those that were read from the lbc_in stream
    !
    call mpas_pool_get_array(lbc, 'lbc_u', u, 2)
    call mpas_pool_get_array(lbc, 'lbc_ru', ru, 2)
    call mpas_pool_get_array(lbc, 'lbc_rho_edge', rho_edge, 2)
    call mpas_pool_get_array(lbc, 'lbc_w', w, 2)
    call mpas_pool_get_array(lbc, 'lbc_theta', theta, 2)
    call mpas_pool_get_array(lbc, 'lbc_rtheta_m', rtheta_m, 2)
    call mpas_pool_get_array(lbc, 'lbc_rho_zz', rho_zz, 2)
    call mpas_pool_get_array(lbc, 'lbc_rho', rho, 2)
    call mpas_pool_get_array(lbc, 'lbc_scalars', scalars, 2)

    call mpas_pool_get_array(mesh, 'cellsOnEdge', cellsOnEdge)
    call mpas_pool_get_dimension(mesh, 'nCells', nCells_ptr)
    call mpas_pool_get_dimension(mesh, 'nEdges', nEdges_ptr)
    call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels_ptr)
    call mpas_pool_get_dimension(state, 'num_scalars', nScalars_ptr)
    call mpas_pool_get_dimension(lbc, 'index_qv', index_qv_ptr)
    call mpas_pool_get_array(mesh, 'zz', zz)

    if (.not. firstCall) then
       call mpas_pool_get_array(lbc, 'lbc_u', lbc_tend_u, 1)
       call mpas_pool_get_array(lbc, 'lbc_ru', lbc_tend_ru, 1)
       call mpas_pool_get_array(lbc, 'lbc_rho_edge', lbc_tend_rho_edge, 1)
       call mpas_pool_get_array(lbc, 'lbc_w', lbc_tend_w, 1)
       call mpas_pool_get_array(lbc, 'lbc_theta', lbc_tend_theta, 1)
       call mpas_pool_get_array(lbc, 'lbc_rtheta_m', lbc_tend_rtheta_m, 1)
       call mpas_pool_get_array(lbc, 'lbc_rho_zz', lbc_tend_rho_zz, 1)
       call mpas_pool_get_array(lbc, 'lbc_rho', lbc_tend_rho, 1)
       call mpas_pool_get_array(lbc, 'lbc_scalars', lbc_tend_scalars, 1)
    endif

    ! Dereference the pointers to avoid non-array pointer for OpenACC
    nCells = nCells_ptr
    nEdges = nEdges_ptr
    nVertLevels = nVertLevels_ptr
    nScalars = nScalars_ptr
    index_qv = index_qv_ptr

    ! Compute lbc_rho_zz
    do k=1,nVertLevels
       zz(k,nCells+1) = 1.0_RKIND          ! Avoid potential division by zero in the following line
    end do

    do iCell=1,nCells+1
       do k=1,nVertLevels
          rho_zz(k,iCell) = rho(k,iCell) / zz(k,iCell)
       end do
    end do

    ! Average lbc_rho_zz to edges
    do iEdge=1,nEdges
       cell1 = cellsOnEdge(1,iEdge)
       cell2 = cellsOnEdge(2,iEdge)
       if (cell1 > 0 .and. cell2 > 0) then
          do k = 1, nVertLevels
             rho_edge(k,iEdge) = 0.5_RKIND * (rho_zz(k,cell1) + rho_zz(k,cell2))
          end do
       end if
    end do

    do iEdge=1,nEdges+1
       do k=1,nVertLevels
          ru(k,iEdge) = u(k,iEdge) * rho_edge(k,iEdge)
       end do
    end do
    
    do iCell=1,nCells+1
       do k=1,nVertLevels
          rtheta_m(k,iCell) = theta(k,iCell) * rho_zz(k,iCell) * (1.0_RKIND + rvord * scalars(index_qv,k,iCell))
       end do
    end do

    if (.not. firstCall) then

       lbc_interval = currTime - LBC_intv_end
       call mpas_get_timeInterval(interval=lbc_interval, DD=dd_intv, S=s_intv, S_n=sn_intv, S_d=sd_intv, ierr=ierr)
       dt = 86400.0_RKIND * real(dd_intv, kind=RKIND) + real(s_intv, kind=RKIND) &
            + (real(sn_intv, kind=RKIND) / real(sd_intv, kind=RKIND))

       dt = 1.0_RKIND / dt

       do iEdge=1,nEdges+1
          do k=1,nVertLevels
             lbc_tend_u(k,iEdge) = (u(k,iEdge) - lbc_tend_u(k,iEdge)) * dt
             lbc_tend_ru(k,iEdge) = (ru(k,iEdge) - lbc_tend_ru(k,iEdge)) * dt
             lbc_tend_rho_edge(k,iEdge) = (rho_edge(k,iEdge) - lbc_tend_rho_edge(k,iEdge)) * dt
          end do
       end do

       do iCell=1,nCells+1
          do k=1,nVertLevels+1
             lbc_tend_w(k,iCell) = (w(k,iCell) - lbc_tend_w(k,iCell)) * dt
          end do
       end do

       do iCell=1,nCells+1
          do k=1,nVertLevels
             lbc_tend_theta(k,iCell) = (theta(k,iCell) - lbc_tend_theta(k,iCell)) * dt
             lbc_tend_rtheta_m(k,iCell) = (rtheta_m(k,iCell) - lbc_tend_rtheta_m(k,iCell)) * dt
             lbc_tend_rho_zz(k,iCell) = (rho_zz(k,iCell) - lbc_tend_rho_zz(k,iCell)) * dt
             lbc_tend_rho(k,iCell) = (rho(k,iCell) - lbc_tend_rho(k,iCell)) * dt
          end do
       end do

       do iCell=1,nCells+1
          do k=1,nVertLevels
             do j = 1,nScalars
                lbc_tend_scalars(j,k,iCell) = (scalars(j,k,iCell) - lbc_tend_scalars(j,k,iCell)) * dt
             end do
          end do
       end do

       !
       ! Logging the lbc start and end times appears to be backwards, but
       ! until the end of this function, LBC_intv_end == the last interval
       ! time and currTime == the next interval time.
       !
       call mpas_get_time(LBC_intv_end, dateTimeString=lbc_intv_start_string)
       call mpas_get_time(currTime, dateTimeString=lbc_intv_end_string)
       call mpas_log_write('----------------------------------------------------------------------')
       call mpas_log_write('Updated lateral boundary conditions. LBCs are now valid')
       call mpas_log_write('from '//trim(lbc_intv_start_string)//' to '//trim(lbc_intv_end_string))
       call mpas_log_write('----------------------------------------------------------------------')

    end if
    LBC_intv_end = currTime
    
  end subroutine ufs_mpas_atm_update_bdy_tend

  !> ########################################################################################
  !
  !  routine ufs_mpas_atm_bdy_checks
  !
  !> \brief   Checks compatibility of limited-area settings
  !> \author  Michael Duda
  !> \date    12 May 2019
  !> \details
  !>  This routine checks that settings related to limited-area simulations
  !>  are compatible. Specifically, the following are checked by this routine:
  !>
  !>  1) If config_apply_lbcs = true, the bdyMaskCell field must have non-zero elements
  !>  2) If config_apply_lbcs = false, the bdyMaskCell field must not have non-zero elements
  !>
  !>  If any of the above are not true, this routine prints an error message and
  !>  returns a non-zero value in ierr; otherwise, a value of 0 is returned.
  !>
  !> \update: Dustin Swales March 2026 - Modified for use in UWM
  !>
  !> ########################################################################################
  subroutine ufs_mpas_atm_bdy_checks(dminfo, blockList, ierr)
    use mpas_log,           only : mpas_log_write
    use mpas_kind_types,    only : StrKIND
    use mpas_derived_types, only : dm_info, block_type, mpas_pool_type, MPAS_LOG_ERR
    use mpas_pool_routines, only : mpas_pool_get_config, mpas_pool_get_dimension
    use mpas_pool_routines, only : mpas_pool_get_subpool, mpas_pool_get_array
    use mpas_dmpar,         only : mpas_dmpar_max_int

    implicit none

    type (dm_info), pointer :: dminfo
    type (block_type), pointer :: blockList
    integer, intent(out) :: ierr

    character(len=StrKIND) :: input_interval
    logical, pointer :: config_apply_lbcs => null()
    integer, pointer :: nCellsSolve => null()
    type (mpas_pool_type), pointer :: meshPool => null()
    type (block_type), pointer :: block => null()
    integer, dimension(:), pointer :: bdyMaskCell => null()
    integer :: maxvar2d_local, maxvar2d_global

    call mpas_pool_get_config(blocklist % configs, 'config_apply_lbcs', config_apply_lbcs)

    call mpas_log_write('')
    call mpas_log_write('Checking consistency of limited-area settings...')
    call mpas_log_write(' - config_apply_lbcs = $l', logicArgs=(/config_apply_lbcs/))

    !
    ! Check whether any elements of bdyMaskCell have non-zero values
    !
    maxvar2d_local = -huge(maxvar2d_local)
    block => blockList
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
         call mpas_pool_get_array(meshPool, 'bdyMaskCell', bdyMaskCell)

         maxvar2d_local = max(maxvar2d_local, maxval(bdyMaskCell(1:nCellsSolve)))

         block => block % next
      end do

      call mpas_dmpar_max_int(dminfo, maxvar2d_local, maxvar2d_global)
      call mpas_log_write(' - Maximum value in bdyMaskCell = $i', intArgs=(/maxvar2d_global/))

      !
      ! If there are boundary cells, config_apply_lbcs must be set to true
      !
      if (.not. config_apply_lbcs .and. maxvar2d_global > 0) then
         call mpas_log_write('Boundary cells found in the bdyMaskCell field, but config_apply_lbcs = false.', &
                              messageType=MPAS_LOG_ERR)
         call mpas_log_write('Please ensure that config_apply_lbcs = true for limited-area simulations.', &
                              messageType=MPAS_LOG_ERR)
         ierr = 1
         return
      end if

      !
      ! If there are no boundary cells, config_apply_lbcs must be set to false
      !
      if (config_apply_lbcs .and. maxvar2d_global == 0) then
         call mpas_log_write('config_apply_lbcs = true, but no boundary cells found in the bdyMaskCell field.', &
                              messageType=MPAS_LOG_ERR)
         call mpas_log_write('Please ensure that config_apply_lbcs = false for global simulations.', &
                              messageType=MPAS_LOG_ERR)
         ierr = 1
         return
      end if

      call mpas_log_write(' ----- done checking limited-area settings -----')
      call mpas_log_write('')
      ierr = 0

  end subroutine ufs_mpas_atm_bdy_checks
end module ufs_mpas_boundaries
