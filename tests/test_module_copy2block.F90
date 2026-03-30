program test_module_copy2block
  !! Unit test driver for copy2block
  !! Tests copy2block across multiple grid decompositions and block sizes

  use ESMF,              only: ESMF_KIND_R8, ESMF_SUCCESS
  use block_control_mod, only: block_control_type, create_block_control
  use CCPP_data,         only: GFS_control
  use GFS_typedefs,      only: GFS_kind_phys => kind_phys
  use atmos_model_mod,   only: copy2block, atmos_model_set_copy2block_test_state

  implicit none

  ! Define test configurations
  integer, parameter :: num_configs = 8
  integer, parameter :: max_grid_size = 16

  type :: test_config_type
     integer :: nx, ny                  !! Grid dimensions
     integer :: inpes, jnpes            !! Decomposition
     integer :: blocksize               !! Block size
     character(len=64) :: description   !! Configuration description
  end type test_config_type

  type(test_config_type) :: configs(num_configs)
  type(block_control_type) :: block_control

  integer :: config_idx, rc, test_count, test_passed
  integer :: current_test_count, current_test_passed

  ! Initialize test counters
  test_count = 0
  test_passed = 0

  print *, "=========================================="
  print *, "Unit Tests: module_block_data"
  print *, "Testing copy2block Across Multiple Decompositions"
  print *, "=========================================="
  print *, " "

  ! Define test configurations
  call setup_test_configurations(configs)

  ! Run tests for each configuration
  do config_idx = 1, num_configs
     print *, "=========================================="
     print *, "Configuration ", config_idx, " of ", num_configs
     print *, "=========================================="
     print *, "Grid Configuration:"
     print *, "  Grid Size: ", configs(config_idx)%nx, " x ", configs(config_idx)%ny
     print *, "  Decomposition (inpes x jnpes): ", configs(config_idx)%inpes, " x ", configs(config_idx)%jnpes
     print *, "  Block Size: ", configs(config_idx)%blocksize
     print *, "  Description: ", trim(configs(config_idx)%description)
     print *, " "

     ! Initialize block control structure
     call initialize_block_control(block_control, &
          configs(config_idx)%nx, &
          configs(config_idx)%ny, &
          configs(config_idx)%inpes, &
          configs(config_idx)%jnpes, &
          configs(config_idx)%blocksize)

     ! Initialize per-configuration counters
     current_test_count = 0
     current_test_passed = 0

     ! Test 1: Block structure initialization
     current_test_count = current_test_count + 1
     call test_block_initialization(block_control, current_test_count, current_test_passed, config_idx)

     ! Test 2: copy2block basic mapping
     current_test_count = current_test_count + 1
     call test_copy2block_basic(block_control, current_test_count, current_test_passed, config_idx)

     ! Test 3: copy2block factor and bounds
     current_test_count = current_test_count + 1
     call test_copy2block_factor_and_bounds(block_control, current_test_count, current_test_passed, config_idx)

     ! Test 4: copy2block mask handling
     current_test_count = current_test_count + 1
     call test_copy2block_mask(block_control, current_test_count, current_test_passed, config_idx)

     ! Test 5: copy2block without mask
     current_test_count = current_test_count + 1
     call test_copy2block_no_mask(block_control, current_test_count, current_test_passed, config_idx)

     ! Test 6: copy2block unassociated targets
     current_test_count = current_test_count + 1
     call test_copy2block_unassociated_targets(block_control, current_test_count, current_test_passed, config_idx)

     ! Test 7: copy2block validmin == validmax
     current_test_count = current_test_count + 1
     call test_copy2block_validmin_eq_validmax(block_control, current_test_count, current_test_passed, config_idx)

     ! Accumulate counts
     test_count = test_count + current_test_count
     test_passed = test_passed + current_test_passed

     ! Print configuration summary
     print *, "Config ", config_idx, " Results: ", current_test_passed, " / ", current_test_count, " passed"
     print *, " "

     ! Clean up block control
     call cleanup_block_control(block_control)
  end do

  ! Print overall summary
  print *, "=========================================="
  print *, "Overall Test Summary"
  print *, "=========================================="
  print *, "Total Tests Passed: ", test_passed, " / ", test_count
  print *, "=========================================="

  if (test_passed == test_count) then
     print *, "All tests passed!"
     stop 0
  else
     print *, "Some tests failed!"
     stop 1
  end if

contains

  !============================================================================
  ! Setup test configurations
  !============================================================================
  subroutine setup_test_configurations(configs)
    type(test_config_type), intent(out) :: configs(:)

    ! Configuration 1: Small 8x8 grid, 2x4 decomposition, blocksize 8
    configs(1)%nx = 8
    configs(1)%ny = 8
    configs(1)%inpes = 2
    configs(1)%jnpes = 4
    configs(1)%blocksize = 8
    configs(1)%description = "8x8, inpes=2, jnpes=4, bs=8"

    ! Configuration 2: 8x8 grid, 2x2 decomposition, blocksize 16
    configs(2)%nx = 8
    configs(2)%ny = 8
    configs(2)%inpes = 2
    configs(2)%jnpes = 2
    configs(2)%blocksize = 16
    configs(2)%description = "8x8, inpes=2, jnpes=2, bs=16"

    ! Configuration 3: 8x8 grid, 1x8 decomposition (linear, uneven X)
    configs(3)%nx = 8
    configs(3)%ny = 8
    configs(3)%inpes = 1
    configs(3)%jnpes = 8
    configs(3)%blocksize = 8
    configs(3)%description = "8x8, inpes=1, jnpes=8, bs=8 (linear in Y)"

    ! Configuration 4: 8x8 grid, 4x2 decomposition (uneven both ways)
    configs(4)%nx = 8
    configs(4)%ny = 8
    configs(4)%inpes = 4
    configs(4)%jnpes = 2
    configs(4)%blocksize = 8
    configs(4)%description = "8x8, inpes=4, jnpes=2, bs=8 (uneven both)"

    ! Configuration 5: 16x16 grid, 1x4 decomposition (linear, uneven Y)
    configs(5)%nx = 16
    configs(5)%ny = 16
    configs(5)%inpes = 1
    configs(5)%jnpes = 4
    configs(5)%blocksize = 16
    configs(5)%description = "16x16, inpes=1, jnpes=4, bs=16 (linear in Y)"

    ! Configuration 6: 16x16 grid, 3x3 decomposition (uneven both ways)
    configs(6)%nx = 16
    configs(6)%ny = 16
    configs(6)%inpes = 3
    configs(6)%jnpes = 3
    configs(6)%blocksize = 16
    configs(6)%description = "16x16, inpes=3, jnpes=3, bs=16 (uneven both)"

    ! Configuration 7: 8x8 grid, 2x6 decomposition (blocksize remainders in Y)
    configs(7)%nx = 8
    configs(7)%ny = 8
    configs(7)%inpes = 2
    configs(7)%jnpes = 6
    configs(7)%blocksize = 8
    configs(7)%description = "8x8, inpes=2, jnpes=6, bs=8 (uneven Y division)"

    ! Configuration 8: 16x16 grid, 2x8 decomposition (many blocks, uneven X)
    configs(8)%nx = 16
    configs(8)%ny = 16
    configs(8)%inpes = 2
    configs(8)%jnpes = 8
    configs(8)%blocksize = 16
    configs(8)%description = "16x16, inpes=2, jnpes=8, bs=16 (many blocks)"

  end subroutine setup_test_configurations

  !============================================================================
  ! Initialize block control structure
  !============================================================================
  subroutine initialize_block_control(block_ctl, nx, ny, inpes, jnpes, blocksize)

    type(block_control_type), intent(out) :: block_ctl
    integer, intent(in) :: nx, ny, inpes, jnpes, blocksize

    integer :: nblocks, iblock, jblock, block_id
    integer :: i, j, istart, jstart, iend, jend
    integer :: npts, ipt

    ! This subroutine simulates a SINGLE MPI task's view of the domain
    ! The full global domain is (nx*inpes) x (ny*jnpes)
    ! This task is responsible for nx x ny portion

    ! Calculate number of blocks within this MPI task
    nblocks = inpes * jnpes

    !Allocate block control arrays
    allocate(block_ctl%blksz(nblocks))
    allocate(block_ctl%index(nblocks))

    !Set local domain bounds for this MPI task
    block_ctl%isc = 1
    block_ctl%iec = nx
    block_ctl%jsc = 1
    block_ctl%jec = ny
    block_ctl%ni = nx
    block_ctl%nj = ny
    block_ctl%nblocks = nblocks

    allocate(block_ctl%blkno(block_ctl%isc:block_ctl%iec, block_ctl%jsc:block_ctl%jec))
    allocate(block_ctl%ixp(block_ctl%isc:block_ctl%iec, block_ctl%jsc:block_ctl%jec))
    block_ctl%blkno = 0
    block_ctl%ixp = 0

    ! Create blocks in a row-major order within this MPI task's domain
    block_id = 1
    do jblock = 1, jnpes
       do iblock = 1, inpes
          ! Calculate block bounds within this task's local domain
          istart = ((iblock - 1) * nx) / inpes + 1
          iend = (iblock * nx) / inpes
          jstart = ((jblock - 1) * ny) / jnpes + 1
          jend = (jblock * ny) / jnpes

          ! Number of points in this block
          npts = (iend - istart + 1) * (jend - jstart + 1)
          block_ctl%blksz(block_id) = npts

          ! Allocate index arrays for this block
          allocate(block_ctl%index(block_id)%ii(npts))
          allocate(block_ctl%index(block_id)%jj(npts))

          ! Populate local domain indices for this block
          ipt = 1
          do j = jstart, jend
             do i = istart, iend
                block_ctl%index(block_id)%ii(ipt) = i
                block_ctl%index(block_id)%jj(ipt) = j
                block_ctl%blkno(i, j) = block_id
                block_ctl%ixp(i, j) = ipt
                ipt = ipt + 1
             end do
          end do
          block_id = block_id + 1
       end do
    end do

    print *, "Block Control Initialized (MPI task local domain):"
    print *, "  Local domain: isc=", block_ctl%isc, " iec=", block_ctl%iec, " jsc=", block_ctl%jsc, " jec=", block_ctl%jec
    print *, "  Number of blocks: ", nblocks
    print *, "  Block sizes range from: 1 to ", maxval(block_ctl%blksz)

  end subroutine initialize_block_control

  !============================================================================
  ! Clean up block control structure
  !============================================================================
  subroutine cleanup_block_control(block_ctl)
    type(block_control_type), intent(inout) :: block_ctl
    integer :: i

    if (allocated(block_ctl%blksz)) deallocate(block_ctl%blksz)
    if (allocated(block_ctl%blkno)) deallocate(block_ctl%blkno)
    if (allocated(block_ctl%ixp)) deallocate(block_ctl%ixp)
    if (allocated(block_ctl%index)) then
       do i = 1, size(block_ctl%index)
          if (allocated(block_ctl%index(i)%ii)) deallocate(block_ctl%index(i)%ii)
          if (allocated(block_ctl%index(i)%jj)) deallocate(block_ctl%index(i)%jj)
       end do
       deallocate(block_ctl%index)
    end if
    if (associated(GFS_control%chunk_begin)) deallocate(GFS_control%chunk_begin)

  end subroutine cleanup_block_control

  subroutine setup_copy2block_state(block_ctl)
    type(block_control_type), intent(in) :: block_ctl

    integer :: block_id, offset

    GFS_control%isc = block_ctl%isc
    GFS_control%jsc = block_ctl%jsc
    GFS_control%nx = block_ctl%ni
    GFS_control%ny = block_ctl%nj
    GFS_control%huge = 9.9692099683868690E30_GFS_kind_phys

    if (associated(GFS_control%chunk_begin)) deallocate(GFS_control%chunk_begin)
    allocate(GFS_control%chunk_begin(block_ctl%nblocks))

    offset = 1
    do block_id = 1, block_ctl%nblocks
       GFS_control%chunk_begin(block_id) = offset
       offset = offset + block_ctl%blksz(block_id)
    end do

    call atmos_model_set_copy2block_test_state(block)

  end subroutine setup_copy2block_state

  integer function packed_index(block_ctl, i, j)
    type(block_control_type), intent(in) :: block_ctl
    integer, intent(in) :: i, j

    packed_index = GFS_control%chunk_begin(block_ctl%blkno(i, j)) + block_ctl%ixp(i, j) - 1

  end function packed_index

  !============================================================================
  ! TEST 1: Block initialization
  !============================================================================
  subroutine test_block_initialization(block, test_num, passed_count, config_idx)
    type(block_control_type), intent(in) :: block
    integer, intent(in) :: test_num, config_idx
    integer, intent(inout) :: passed_count

    integer :: total_pts, computed_total
    integer :: i

    print *, "  [Config ", config_idx, "] Test ", test_num, ": Block Structure Initialization"

    ! Check that block size array is properly initialized
    if (.not. allocated(block_ctl%blksz)) then
       print *, "    FAILED: blksz array not allocated"
       return
    end if

    ! Check grid dimensions
    if (block_ctl%ni /= block_ctl%iec - block_ctl%isc + 1 .or. block_ctl%nj /= block_ctl%jec - block_ctl%jsc + 1) then
       print *, "    FAILED: Grid dimensions incorrect"
       return
    end if

    ! Check that total points match grid size
    total_pts = (block_ctl%iec - block_ctl%isc + 1) * (block_ctl%jec - block_ctl%jsc + 1)
    computed_total = sum(block_ctl%blksz)
    if (computed_total /= total_pts) then
       print *, "    FAILED: Total points mismatch"
       print *, "      Expected: ", total_pts, " Computed: ", computed_total
       return
    end if

    ! Check that global indices are within domain bounds
    do i = 1, block_ctl%nblocks
       if (minval(block_ctl%index(i)%ii) < block_ctl%isc) then
          print *, "    FAILED: I index below isc in block ", i
          return
       end if
       if (maxval(block_ctl%index(i)%ii) > block_ctl%iec) then
          print *, "    FAILED: I index exceeds iec in block ", i
          return
       end if
       if (minval(block_ctl%index(i)%jj) < block_ctl%jsc) then
          print *, "    FAILED: J index below jsc in block ", i
          return
       end if
       if (maxval(block_ctl%index(i)%jj) > block_ctl%jec) then
          print *, "    FAILED: J index exceeds jec in block ", i
          return
       end if
    end do

    print *, "    PASSED"
    passed_count = passed_count + 1
  end subroutine test_block_initialization

  !============================================================================
  ! TEST: copy2block basic mapping
  !============================================================================
  subroutine test_copy2block_basic(block, test_num, passed_count, config_idx)
    type(block_control_type), intent(in) :: block
    integer, intent(in) :: test_num, config_idx
    integer, intent(inout) :: passed_count

    real(GFS_kind_phys), allocatable :: dest_1d(:), mask(:)
    real(ESMF_KIND_R8), allocatable :: source_2d(:,:)

    integer :: i, j, im, rc, total_pts
    logical :: test_pass

    print *, "  [Config ", config_idx, "] Test ", test_num, ": copy2block Basic Mapping"

    call setup_copy2block_state(block)
    total_pts = sum(block%blksz)
    allocate(dest_1d(total_pts), mask(total_pts))
    allocate(source_2d(block%ni, block%nj))
    dest_1d = -1.0_GFS_kind_phys
    mask = 1.0_GFS_kind_phys
    do j = 1, block%nj
       do i = 1, block%ni
          source_2d(i, j) = real(100 * j + i, ESMF_KIND_R8)
       end do
    end do

    call copy2block(dest_1d, source_2d, mask, rc=rc)

    test_pass = rc == ESMF_SUCCESS
    if (.not. test_pass) print *, "    FAILED: copy2block returned rc=", rc
    if (test_pass) then
       do j = block%jsc, block%jec
          do i = block%isc, block%iec
             im = packed_index(block, i, j)
             if (abs(dest_1d(im) - real(source_2d(i - block%isc + 1, j - block%jsc + 1), GFS_kind_phys)) > 1.0e-10_GFS_kind_phys) then
                print *, "    FAILED: Data value mismatch at i=", i, " j=", j
                test_pass = .false.
                exit
             end if
          end do
          if (.not. test_pass) exit
       end do
    end if
    if (test_pass) then
       print *, "    PASSED"
       passed_count = passed_count + 1
    else
       print *, "    FAILED: copy2block basic mapping verification failed"
    end if
    deallocate(dest_1d, mask, source_2d)
  end subroutine test_copy2block_basic

  !============================================================================
  ! TEST: copy2block factor and bounds
  !============================================================================
  subroutine test_copy2block_factor_and_bounds(block, test_num, passed_count, config_idx)
    type(block_control_type), intent(in) :: block
    integer, intent(in) :: test_num, config_idx
    integer, intent(inout) :: passed_count

    real(GFS_kind_phys), allocatable :: dest_1d(:), mask(:)
    real(ESMF_KIND_R8), allocatable :: source_2d(:,:)

    real(GFS_kind_phys), parameter :: factor = -2.0_GFS_kind_phys
    real(GFS_kind_phys), parameter :: validmin = 1.0_GFS_kind_phys
    real(GFS_kind_phys), parameter :: validmax = 5.0_GFS_kind_phys
    real(GFS_kind_phys) :: expected

    integer :: i, j, im, rc, total_pts
    logical :: test_pass

    print *, "  [Config ", config_idx, "] Test ", test_num, ": copy2block Factor and Bounds"

    call setup_copy2block_state(block)
    total_pts = sum(block%blksz)
    allocate(dest_1d(total_pts), mask(total_pts))
    allocate(source_2d(block%ni, block%nj))
    dest_1d = -1.0_GFS_kind_phys
    mask = 1.0_GFS_kind_phys
    do j = 1, block%nj
       do i = 1, block%ni
          source_2d(i, j) = real(i - 2 * j, ESMF_KIND_R8)
       end do
    end do

    call copy2block(dest_1d, source_2d, mask, validmin=validmin, validmax=validmax, factor=factor, rc=rc)

    test_pass = rc == ESMF_SUCCESS
    if (.not. test_pass) print *, "    FAILED: copy2block returned rc=", rc
    if (test_pass) then
       do j = block%jsc, block%jec
          do i = block%isc, block%iec
             im = packed_index(block, i, j)
             expected = real(source_2d(i - block%isc + 1, j - block%jsc + 1), GFS_kind_phys) * factor
             expected = max(validmin, expected)
             expected = min(validmax, expected)
             if (abs(dest_1d(im) - expected) > 1.0e-10_GFS_kind_phys) then
                print *, "    FAILED: Factor/bounds mismatch at i=", i, " j=", j
                test_pass = .false.
                exit
             end if
          end do
          if (.not. test_pass) exit
       end do
    end if
    if (test_pass) then
       print *, "    PASSED"
       passed_count = passed_count + 1
    else
       print *, "    FAILED: copy2block factor/bounds verification failed"
    end if
    deallocate(dest_1d, mask, source_2d)
  end subroutine test_copy2block_factor_and_bounds

  !============================================================================
  ! TEST: copy2block mask handling
  !============================================================================
  subroutine test_copy2block_mask(block, test_num, passed_count, config_idx)
    type(block_control_type), intent(in) :: block
    integer, intent(in) :: test_num, config_idx
    integer, intent(inout) :: passed_count

    real(GFS_kind_phys), allocatable :: dest_1d(:), mask(:)
    real(ESMF_KIND_R8), allocatable :: source_2d(:,:)
    real(GFS_kind_phys), parameter :: sentinel = -777.0_GFS_kind_phys

    integer :: i, j, im, rc, total_pts
    logical :: test_pass

    print *, "  [Config ", config_idx, "] Test ", test_num, ": copy2block Mask Handling"

    call setup_copy2block_state(block)
    total_pts = sum(block%blksz)
    allocate(dest_1d(total_pts), mask(total_pts))
    allocate(source_2d(block%ni, block%nj))
    dest_1d = sentinel
    mask = 0.0_GFS_kind_phys
    do j = 1, block%nj
       do i = 1, block%ni
          source_2d(i, j) = real(10 * j + i, ESMF_KIND_R8)
          if (mod(i + j, 2) == 0) then
             im = packed_index(block, i, j)
             mask(im) = 1.0_GFS_kind_phys
          end if
       end do
    end do

    call copy2block(dest_1d, source_2d, mask, rc=rc)

    test_pass = rc == ESMF_SUCCESS
    if (.not. test_pass) print *, "    FAILED: copy2block returned rc=", rc
    if (test_pass) then
       do j = block%jsc, block%jec
          do i = block%isc, block%iec
             im = packed_index(block, i, j)
             if (mask(im) > 0.0_GFS_kind_phys) then
                if (abs(dest_1d(im) - real(source_2d(i - block%isc + 1, j - block%jsc + 1), GFS_kind_phys)) > 1.0e-10_GFS_kind_phys) then
                   print *, "    FAILED: Masked point not copied at i=", i, " j=", j
                   test_pass = .false.
                   exit
                end if
             else if (abs(dest_1d(im) - sentinel) > 1.0e-10_GFS_kind_phys) then
                print *, "    FAILED: Unmasked point was modified at i=", i, " j=", j
                test_pass = .false.
                exit
             end if
          end do
          if (.not. test_pass) exit
       end do
    end if
    if (test_pass) then
       print *, "    PASSED"
       passed_count = passed_count + 1
    else
       print *, "    FAILED: copy2block mask verification failed"
    end if
    deallocate(dest_1d, mask, source_2d)
  end subroutine test_copy2block_mask

  !============================================================================
  ! TEST: copy2block without mask
  !============================================================================
  subroutine test_copy2block_no_mask(block, test_num, passed_count, config_idx)
    type(block_control_type), intent(in) :: block
    integer, intent(in) :: test_num, config_idx
    integer, intent(inout) :: passed_count

    real(GFS_kind_phys), allocatable :: dest_1d(:)
    real(ESMF_KIND_R8), allocatable :: source_2d(:,:)
    real(GFS_kind_phys), parameter :: sentinel = -999.0_GFS_kind_phys

    integer :: i, j, im, rc, total_pts
    logical :: test_pass

    print *, "  [Config ", config_idx, "] Test ", test_num, ": copy2block Without Mask"

    call setup_copy2block_state(block)
    total_pts = sum(block%blksz)
    allocate(dest_1d(total_pts))
    allocate(source_2d(block%ni, block%nj))

    dest_1d = sentinel
    do j = 1, block%nj
       do i = 1, block%ni
          source_2d(i, j) = real(1000 + 10 * j + i, ESMF_KIND_R8)
       end do
    end do

    call copy2block(dest_1d, source_2d, rc=rc)

    test_pass = rc == ESMF_SUCCESS
    if (.not. test_pass) print *, "    FAILED: copy2block returned rc=", rc

    if (test_pass) then
       do j = block%jsc, block%jec
          do i = block%isc, block%iec
             im = packed_index(block, i, j)
             if (abs(dest_1d(im) - real(source_2d(i - block%isc + 1, j - block%jsc + 1), GFS_kind_phys)) > 1.0e-10_GFS_kind_phys) then
                print *, "    FAILED: Value mismatch without mask at i=", i, " j=", j
                test_pass = .false.
                exit
             end if
          end do
          if (.not. test_pass) exit
       end do
    end if

    if (test_pass) then
       print *, "    PASSED"
       passed_count = passed_count + 1
    else
       print *, "    FAILED: copy2block without-mask verification failed"
    end if

    deallocate(dest_1d, source_2d)
  end subroutine test_copy2block_no_mask

  !============================================================================
  ! TEST: copy2block mismatched shapes
  !============================================================================
  subroutine test_copy2block_mismatched_shapes(block, test_num, passed_count, config_idx)
    type(block_control_type), intent(in) :: block
    integer, intent(in) :: test_num, config_idx
    integer, intent(inout) :: passed_count
    real(GFS_kind_phys), allocatable :: dest_1d(:), mask(:)
    real(ESMF_KIND_R8), allocatable :: source_2d(:,:)
    integer :: rc, total_pts
    logical :: test_pass
    print *, "  [Config ", config_idx, "] Test ", test_num, ": copy2block Mismatched Shapes"
    call setup_copy2block_state(block)
    total_pts = sum(block%blksz)
    allocate(dest_1d(total_pts), mask(total_pts))
    allocate(source_2d(block%ni+1, block%nj))  ! Deliberate shape mismatch
    dest_1d = -1.0_GFS_kind_phys
    mask = 1.0_GFS_kind_phys
    source_2d = 1.0
    rc = -999
    call copy2block(dest_1d, source_2d, mask, rc=rc)
    test_pass = rc /= ESMF_SUCCESS
    if (test_pass) then
       print *, "    PASSED (error detected as expected)"
       passed_count = passed_count + 1
    else
       print *, "    FAILED: copy2block did not detect shape mismatch"
    end if
    deallocate(dest_1d, mask, source_2d)
  end subroutine test_copy2block_mismatched_shapes

  !============================================================================
  ! TEST: copy2block invalid mask
  !============================================================================
  subroutine test_copy2block_invalid_mask(block, test_num, passed_count, config_idx)
    type(block_control_type), intent(in) :: block
    integer, intent(in) :: test_num, config_idx
    integer, intent(inout) :: passed_count

    real(GFS_kind_phys), allocatable :: dest_1d(:), mask(:)
    real(ESMF_KIND_R8), allocatable :: source_2d(:,:)

    integer :: rc, total_pts
    logical :: test_pass

    print *, "  [Config ", config_idx, "] Test ", test_num, ": copy2block Invalid Mask"

    call setup_copy2block_state(block)
    total_pts = sum(block%blksz)
    allocate(dest_1d(total_pts), mask(total_pts-1))  ! Deliberate mask size error
    allocate(source_2d(block%ni, block%nj))
    dest_1d = -1.0_GFS_kind_phys
    mask = 1.0_GFS_kind_phys
    source_2d = 1.0
    rc = -999

    call copy2block(dest_1d, source_2d, mask, rc=rc)
    test_pass = rc /= ESMF_SUCCESS
    if (test_pass) then
       print *, "    PASSED (error detected as expected)"
       passed_count = passed_count + 1
    else
       print *, "    FAILED: copy2block did not detect invalid mask"
    end if
    deallocate(dest_1d, mask, source_2d)
  end subroutine test_copy2block_invalid_mask

  !============================================================================
  ! TEST: copy2block validmin == validmax
  !============================================================================
  subroutine test_copy2block_validmin_eq_validmax(block, test_num, passed_count, config_idx)
    type(block_control_type), intent(in) :: block
    integer, intent(in) :: test_num, config_idx
    integer, intent(inout) :: passed_count

    real(GFS_kind_phys), allocatable :: dest_1d(:), mask(:)
    real(ESMF_KIND_R8), allocatable :: source_2d(:,:)
    real(GFS_kind_phys), parameter :: val = 3.14_GFS_kind_phys

    integer :: rc, total_pts
    logical :: test_pass

    print *, "  [Config ", config_idx, "] Test ", test_num, ": copy2block validmin == validmax"

    call setup_copy2block_state(block)

    total_pts = sum(block%blksz)
    allocate(dest_1d(total_pts), mask(total_pts))
    allocate(source_2d(block%ni, block%nj))
    dest_1d = -1.0_GFS_kind_phys
    mask = 1.0_GFS_kind_phys
    source_2d = 2.0
    rc = -999

    call copy2block(dest_1d, source_2d, mask, validmin=val, validmax=val, rc=rc)

    test_pass = rc /= ESMF_SUCCESS
    if (test_pass) then
       print *, "    PASSED (error detected as expected)"
       passed_count = passed_count + 1
    else
       print *, "    FAILED: copy2block did not detect validmin == validmax"
    end if
    deallocate(dest_1d, mask, source_2d)
  end subroutine test_copy2block_validmin_eq_validmax

  !============================================================================
  ! TEST: copy2block unassociated targets
  !============================================================================
  subroutine test_copy2block_unassociated_targets(block, test_num, passed_count, config_idx)
    type(block_control_type), intent(in) :: block
    integer, intent(in) :: test_num, config_idx
    integer, intent(inout) :: passed_count

    real(GFS_kind_phys), pointer :: dest_1d(:) => null()
    real(ESMF_KIND_R8), pointer :: source_2d(:,:) => null()
    real(GFS_kind_phys), pointer :: mask(:) => null()

    integer :: rc
    logical :: test_pass

    print *, "  [Config ", config_idx, "] Test ", test_num, ": copy2block Unassociated Targets"
    ! Do not allocate any arrays; pass as null pointers
    rc = -999

    call copy2block(dest_1d, source_2d, mask, rc=rc)
    test_pass = rc /= ESMF_SUCCESS
    if (test_pass) then
       print *, "    PASSED (error detected as expected)"
       passed_count = passed_count + 1
    else
       print *, "    FAILED: copy2block did not detect unassociated targets"
    end if
  end subroutine test_copy2block_unassociated_targets

end program test_module_copy2block
