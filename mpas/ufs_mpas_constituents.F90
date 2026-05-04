!> ###########################################################################################
!> \file ufs_mpas_constituents.F90
!>
!> This module contains the interface between MPAS constituents and the UFS Weather Model.
!> 
!> ###########################################################################################
module ufs_mpas_constituents
  use mpas_kind_types,  only : StrKIND
  use ufs_mpas_io,      only : domain_ptr
  implicit none

  public

  ! These are setup during ATM initialization.
  character(StrKIND), allocatable :: constituent_name(:)
  integer, allocatable :: index_constituent_to_mpas_scalar(:)
  integer, allocatable :: index_mpas_scalar_to_constituent(:)
  logical, allocatable :: is_water_species(:)
contains
  !> #########################################################################################
  !>
  !> \brief  Define the names of constituents at run-time
  !> \author Michael Duda
  !> \date   21 May 2020
  !> \details
  !>  Given an array of constituent names, which must have size equal to the number
  !>  of scalars that were set in the call to ufs_mpas_init_phase1, and given
  !>  a function to identify which scalars are moisture species, this routine defines
  !>  scalar constituents for the MPAS-A dycore.
  !>  Because the MPAS-A dycore expects all moisture constituents to appear in
  !>  a contiguous range of constituent indices, this routine may in general need
  !>  to reorder the constituents; to allow for mapping of indices between UFS
  !>  physics and the MPAS-A dycore, this routine returns index mapping arrays
  !>  mpas_from_ufs_cnst and ufs_from_mpas_cnst.
  !>
  !> \update: Dustin Swales April 2025 - Modified for use in UWM  
  !>
  !> #########################################################################################
  subroutine ufs_mpas_define_scalars(mpas_from_ufs_cnst, ufs_from_mpas_cnst, ierr)
    use mpas_derived_types, only : mpas_pool_type, field3dReal
    use mpas_pool_routines, only : mpas_pool_get_subpool, mpas_pool_get_field, &
                                   mpas_pool_get_dimension, mpas_pool_add_dimension
    use mpas_attlist,       only : mpas_add_att
    use mpas_log,           only : mpas_log_write
    use mpas_derived_types, only : MPAS_LOG_ERR
    ! FMS
    use mpp_mod,            only : FATAL, mpp_error

    ! Arguments
    integer, dimension(:), pointer :: mpas_from_ufs_cnst, ufs_from_mpas_cnst
    integer, intent(out) :: ierr

    ! Local variables
    character(len=*), parameter :: subname = 'ufs_mpas_subdriver::ufs_mpas_define_scalars'
    integer :: i, j, timeLevs
    integer, pointer :: num_scalars
    integer :: num_moist
    integer :: idx_passive
    type (mpas_pool_type), pointer :: statePool
    type (mpas_pool_type), pointer :: tendPool
    type (field3dReal), pointer :: scalarsField
    character(len=128) :: tempstr
    character :: moisture_char

    ierr = 0

    !
    ! Define scalars
    !
    nullify(statePool)
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'state', statePool)

    if (.not. associated(statePool)) then
       call mpas_log_write(trim(subname)//': ERROR: The ''state'' pool was not found.', &
                           messageType=MPAS_LOG_ERR)
       ierr = 1
       return
    end if

    nullify(num_scalars)
    call mpas_pool_get_dimension(statePool, 'num_scalars', num_scalars)

    !
    ! The num_scalars dimension should have been defined by atm_core_interface::atm_allocate_scalars, and
    ! if this dimension does not exist, something has gone wrong
    !
    if (.not. associated(num_scalars)) then
       call mpas_log_write(trim(subname)//': ERROR: The ''num_scalars'' dimension does not exist in the ''state'' pool.', &
                           messageType=MPAS_LOG_ERR)
       ierr = 1
       return
    end if

    !
    ! If at runtime there are not num_scalars names in the array of constituent names provided by UFS,
    ! something has gone wrong
    !
    if (size(constituent_name) /= num_scalars) then
       call mpas_log_write(trim(subname)//': ERROR: The number of constituent names is not equal to the num_scalars dimension', &
                           messageType=MPAS_LOG_ERR)
       call mpas_log_write('size(constituent_name) = $i, num_scalars = $i', intArgs=[size(constituent_name), num_scalars], &
                           messageType=MPAS_LOG_ERR)
       ierr = 1
       return
    end if

    !
    ! In UFS, the first scalar (if there are any) is always qv (specific humidity); if this is not
    ! the case, something has gone wrong
    !
    if (size(constituent_name) > 0) then
       if (trim(constituent_name(1)) /= 'qv') then
          call mpas_log_write(trim(subname)//': ERROR: The first constituent is not qv', messageType=MPAS_LOG_ERR)
          ierr = 1
          return
       end if
    end if

    !
    ! Determine which of the constituents are moisture species
    !
    allocate(mpas_from_ufs_cnst(num_scalars), stat=ierr)
    if( ierr /= 0 ) call mpp_error(FATAL,subname//':failed to allocate mpas_from_ufs_cnst array')
    mpas_from_ufs_cnst(:) = 0
    num_moist = 0
    do i = 1, size(constituent_name)
       if (is_water_species(i)) then
          num_moist = num_moist + 1
          mpas_from_ufs_cnst(num_moist) = i
       end if
    end do

    !
    ! If UFS has no scalars, let the only scalar in MPAS be 'qv' (a moisture species)
    !
    if (num_scalars == 1 .and. size(constituent_name) == 0) then
       num_moist = 1
    end if

    !
    ! Assign non-moisture constituents to mpas_from_ufs_cnst(num_moist+1:size(constituent_name))
    !
    idx_passive = num_moist + 1
    do i = 1, size(constituent_name)
       ! If UFS constituent i is not already mapped as a moist constituent
       if (.not. is_water_species(i)) then
          mpas_from_ufs_cnst(idx_passive) = i
          idx_passive = idx_passive + 1
       end if
    end do

    !
    ! Create inverse map, ufs_from_mpas_cnst
    !
    allocate(ufs_from_mpas_cnst(num_scalars), stat=ierr)
    if( ierr /= 0 ) call mpp_error(FATAL,subname//':failed to allocate ufs_from_mpas_cnst array')
    ufs_from_mpas_cnst(:) = 0

    do i = 1, size(constituent_name)
       ufs_from_mpas_cnst(mpas_from_ufs_cnst(i)) = i
    end do

    timeLevs = 2

    do i = 1, timeLevs
       nullify(scalarsField)
       call mpas_pool_get_field(statePool, 'scalars', scalarsField, timeLevel=i)

       if (.not. associated(scalarsField)) then
          call mpas_log_write(trim(subname)//': ERROR: The ''scalars'' field was not found in the ''state'' pool', &
                              messageType=MPAS_LOG_ERR)
          ierr = 1
          return
       end if

       if (i == 1) call mpas_pool_add_dimension(statePool, 'index_qv', 1)
       scalarsField % constituentNames(1) = 'qv'
       call mpas_add_att(scalarsField % attLists(1) % attList, 'units', 'kg kg^{-1}')
       call mpas_add_att(scalarsField % attLists(1) % attList, 'long_name', 'Water vapor mixing ratio')

       do j = 2, size(constituent_name)
          scalarsField % constituentNames(j) = trim(constituent_name(mpas_from_ufs_cnst(j)))
       end do

    end do

    call mpas_pool_add_dimension(statePool, 'moist_start', 1)
    call mpas_pool_add_dimension(statePool, 'moist_end', num_moist)

    !
    ! Print a tabular summary of the mapping between constituent indices
    !
    call mpas_log_write('')
    call mpas_log_write('  i MPAS constituent mpas_from_ufs_cnst(i)       i UFS constituent  ufs_from_mpas_cnst(i)')
    call mpas_log_write('------------------------------------------     ------------------------------------------')
    do i = 1, min(num_scalars, size(constituent_name))
       if (i <= num_moist) then
          moisture_char = '*'
       else
          moisture_char = ' '
       end if
       write(tempstr, '(i3,1x,a16,1x,i18,8x,i3,1x,a16,1x,i18)') i, trim(scalarsField % constituentNames(i))//moisture_char, &
                                                                mpas_from_ufs_cnst(i), &
                                                                i, trim(constituent_name(i)), &
                                                                ufs_from_mpas_cnst(i)
       call mpas_log_write(trim(tempstr))
    end do
    call mpas_log_write('------------------------------------------     ------------------------------------------')
    call mpas_log_write('* = constituent used as a moisture species in MPAS-A dycore')
    call mpas_log_write('')

    !
    ! Define scalars_tend
    !
    nullify(tendPool)
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'tend', tendPool)

    if (.not. associated(tendPool)) then
       call mpas_log_write(trim(subname)//': ERROR: The ''tend'' pool was not found.', &
                           messageType=MPAS_LOG_ERR)
       ierr = 1
       return
    end if

    timeLevs = 1

    do i = 1, timeLevs
       nullify(scalarsField)
       call mpas_pool_get_field(tendPool, 'scalars_tend', scalarsField, timeLevel=i)

       if (.not. associated(scalarsField)) then
          call mpas_log_write(trim(subname)//': ERROR: The ''scalars_tend'' field was not found in the ''tend'' pool', &
                              messageType=MPAS_LOG_ERR)
          ierr = 1
          return
       end if

       if (i == 1) call mpas_pool_add_dimension(tendPool, 'index_qv', 1)
       scalarsField % constituentNames(1) = 'tend_qv'
       call mpas_add_att(scalarsField % attLists(1) % attList, 'units', 'kg m^{-3} s^{-1}')
       call mpas_add_att(scalarsField % attLists(1) % attList, 'long_name', 'Tendency of water vapor mixing ratio')

       do j = 2, size(constituent_name)
          scalarsField % constituentNames(j) = 'tend_'//trim(constituent_name(mpas_from_ufs_cnst(j)))
       end do
    end do

    call mpas_pool_add_dimension(tendPool, 'moist_start', 1)
    call mpas_pool_add_dimension(tendPool, 'moist_end', num_moist)

  end subroutine ufs_mpas_define_scalars

  !> #########################################################################################
  !>
  !> \brief  Define the names of lateral-boundary condition constituents at run-time.
  !> \author Dustin Swales
  !> \date   01 March 2026
  !> \details
  !>  Follows ufs_mpas_define_scalars, but for scalars in the LBC pool.
  !>
  !> #########################################################################################
  subroutine ufs_mpas_define_lbc_scalars(mpas_from_ufs_cnst, ufs_from_mpas_cnst, ierr)
    use mpas_derived_types, only : mpas_pool_type, field3dReal, MPAS_LOG_ERR
    use mpas_pool_routines, only : mpas_pool_get_subpool, mpas_pool_get_field
    use mpas_pool_routines, only : mpas_pool_get_dimension, mpas_pool_add_dimension
    use mpas_attlist,       only : mpas_add_att
    use mpas_log,           only : mpas_log_write

    ! Arguments
    integer, dimension(:), pointer :: mpas_from_ufs_cnst, ufs_from_mpas_cnst
    integer, intent(out) :: ierr

    ! Local variables
    character(len=*), parameter :: subname = 'ufs_mpas_subdriver::ufs_mpas_define_lbc_scalars'
    type (mpas_pool_type), pointer :: lbcPool
    integer, pointer :: num_scalars
    integer :: i, j, timeLevs, num_moist
    type (field3dReal), pointer :: scalarsField

    ierr = 0

    !
    ! Define lbc_scalars
    !
    nullify(lbcPool)
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'lbc', lbcPool)

    if (.not. associated(lbcPool)) then
       call mpas_log_write(trim(subname)//': ERROR: The ''lbc'' pool was not found.', &
                           messageType=MPAS_LOG_ERR)
       ierr = 1
       return
    end if

    nullify(num_scalars)
    call mpas_pool_get_dimension(lbcPool, 'num_scalars', num_scalars)

    !
    ! The num_scalars dimension should have been defined by atm_core_interface::atm_allocate_lbc_scalars, and
    ! if this dimension does not exist, something has gone wrong.
    !
    if (.not. associated(num_scalars)) then
       call mpas_log_write(trim(subname)//': ERROR: The ''num_scalars'' dimension does not exist in the ''lbc'' pool.', &
                           messageType=MPAS_LOG_ERR)
       ierr = 1
       return
    end if

    timeLevs = 2

    do i = 1, timeLevs
       nullify(scalarsField)
       call mpas_pool_get_field(lbcPool, 'lbc_scalars', scalarsField, timeLevel=i)

       if (.not. associated(scalarsField)) then
          call mpas_log_write(trim(subname)//': ERROR: The ''lbc_scalars'' field was not found in the ''lbc'' pool', &
                              messageType=MPAS_LOG_ERR)
          ierr = 1
          return
       end if

       if (i == 1) call mpas_pool_add_dimension(lbcPool, 'index_qv', 1)
       scalarsField % constituentNames(1) = 'lbc_qv'
       call mpas_add_att(scalarsField % attLists(1) % attList, 'units', 'kg kg^{-1}')
       call mpas_add_att(scalarsField % attLists(1) % attList, 'long_name', 'Water vapor mixing ratio')

       do j = 2, size(constituent_name)
          scalarsField % constituentNames(j) = 'lbc_'//trim(constituent_name(mpas_from_ufs_cnst(j)))
       end do

    end do

    ! Define lbc_scalars_tend
    ! DJS: No need to do this for LBCs. Tendency/State for LBC stored in LBC pool created
    !      in ufs_mpas_update_bdy_tend()

    call mpas_pool_add_dimension(lbcPool, 'moist_start', 1)
    call mpas_pool_add_dimension(lbcPool, 'moist_end', num_moist)

  end subroutine ufs_mpas_define_lbc_scalars
end module ufs_mpas_constituents
