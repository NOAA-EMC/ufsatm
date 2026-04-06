program test_copy2block
  !! Unit test driver for copy2block
  !! Tests copy2block across multiple grid decompositions and block sizes

  use ESMF,              only: ESMF_SUCCESS
  use block_control_mod, only: block_control_type
  use GFS_typedefs,      only: GFS_kind_phys => kind_phys
  use atmos_model_mod,   only: copy2block
  use module_block_harness, only: test_config_type, setup_test_configurations, &
       initialize_block_control, cleanup_block_control, &
       setup_copy2block_state, packed_index

  implicit none

  ! Define test configurations
  integer, parameter :: num_configs = 8

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

     ! Test 3: copy2block flip and bounds
     current_test_count = current_test_count + 1
     call test_copy2block_flip_and_bounds(block_control, current_test_count, current_test_passed, config_idx)

     ! Test 4: copy2block mask handling
     current_test_count = current_test_count + 1
     call test_copy2block_mask(block_control, current_test_count, current_test_passed, config_idx)

     ! Test 5: copy2block without mask
     current_test_count = current_test_count + 1
     call test_copy2block_no_mask(block_control, current_test_count, current_test_passed, config_idx)

     ! Test 6: copy2block source X shape validation
     current_test_count = current_test_count + 1
     call test_copy2block_mismatched_shapes_x_alt(block_control, current_test_count, current_test_passed, config_idx)

     ! Test 7: copy2block mismatched source shape
     current_test_count = current_test_count + 1
     call test_copy2block_mismatched_shapes(block_control, current_test_count, current_test_passed, config_idx)

     ! Test 8: copy2block source XY shape validation
     current_test_count = current_test_count + 1
     call test_copy2block_mismatched_shapes_xy(block_control, current_test_count, current_test_passed, config_idx)

     ! Test 9: copy2block validmin only
     current_test_count = current_test_count + 1
     call test_copy2block_validmin_only(block_control, current_test_count, current_test_passed, config_idx)

     ! Test 10: copy2block validmax only
     current_test_count = current_test_count + 1
     call test_copy2block_validmax_only(block_control, current_test_count, current_test_passed, config_idx)

     ! Test 11: copy2block mixed copied/skipped regression
     current_test_count = current_test_count + 1
     call test_copy2block_mixed_skip_regression(block_control, current_test_count, current_test_passed, config_idx)

     ! Test 12: copy2block flip only
     current_test_count = current_test_count + 1
     call test_copy2block_flip_only(block_control, current_test_count, current_test_passed, config_idx)

     ! Test 13: copy2block default flip
     current_test_count = current_test_count + 1
     call test_copy2block_default_flip(block_control, current_test_count, current_test_passed, config_idx)

     ! Test 14: copy2block negative mask values
     current_test_count = current_test_count + 1
     call test_copy2block_negative_mask_values(block_control, current_test_count, current_test_passed, config_idx)

     ! Test 15: copy2block mismatched source Y shape
     current_test_count = current_test_count + 1
     call test_copy2block_mismatched_shapes_y(block_control, current_test_count, current_test_passed, config_idx)

     ! Test 16: copy2block oversized source/mask acceptance
     current_test_count = current_test_count + 1
     call test_copy2block_oversized_inputs(block_control, current_test_count, current_test_passed, config_idx)

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
    if (.not. allocated(block%blksz)) then
       print *, "    FAILED: blksz array not allocated"
       return
    end if

    ! Check grid dimensions
    ! Check that total points match grid size
    total_pts = (block%iec - block%isc + 1) * (block%jec - block%jsc + 1)
    computed_total = sum(block%blksz)
    if (computed_total /= total_pts) then
       print *, "    FAILED: Total points mismatch"
       print *, "      Expected: ", total_pts, " Computed: ", computed_total
       return
    end if

    ! Check that global indices are within domain bounds
    do i = 1, block%nblks
       if (minval(block%index(i)%ii) < block%isc) then
          print *, "    FAILED: I index below isc in block ", i
          return
       end if
       if (maxval(block%index(i)%ii) > block%iec) then
          print *, "    FAILED: I index exceeds iec in block ", i
          return
       end if
       if (minval(block%index(i)%jj) < block%jsc) then
          print *, "    FAILED: J index below jsc in block ", i
          return
       end if
       if (maxval(block%index(i)%jj) > block%jec) then
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
    real(GFS_kind_phys), allocatable :: source_2d(:,:)

    integer :: i, j, im, rc, total_pts, nx_local, ny_local
    logical :: test_pass

    print *, "  [Config ", config_idx, "] Test ", test_num, ": copy2block Basic Mapping"

    call setup_copy2block_state(block)
    nx_local = block%iec - block%isc + 1
    ny_local = block%jec - block%jsc + 1
    total_pts = sum(block%blksz)

    allocate(dest_1d(total_pts), mask(total_pts))
    allocate(source_2d(nx_local, ny_local))

    dest_1d = -1.0_GFS_kind_phys
    mask = 1.0_GFS_kind_phys

    do j = 1, ny_local
       do i = 1, nx_local
          source_2d(i, j) = real(100 * j + i, GFS_kind_phys)
       end do
    end do

    call copy2block(dest_1d, source_2d, mask, block=block, rc=rc)

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
  ! TEST: copy2block flip and bounds
  !============================================================================
  subroutine test_copy2block_flip_and_bounds(block, test_num, passed_count, config_idx)
    type(block_control_type), intent(in) :: block
    integer, intent(in) :: test_num, config_idx
    integer, intent(inout) :: passed_count

    real(GFS_kind_phys), allocatable :: dest_1d(:), mask(:)
    real(GFS_kind_phys), allocatable :: source_2d(:,:)

    logical, parameter :: flip = .true.
    real(GFS_kind_phys), parameter :: validmin = 1.0_GFS_kind_phys
    real(GFS_kind_phys), parameter :: validmax = 5.0_GFS_kind_phys
    real(GFS_kind_phys), parameter :: sentinel = -333.0_GFS_kind_phys
    real(GFS_kind_phys) :: expected

    integer :: i, j, im, rc, total_pts, nx_local, ny_local
    logical :: test_pass

    print *, "  [Config ", config_idx, "] Test ", test_num, ": copy2block Flip and Bounds"

    call setup_copy2block_state(block)
    nx_local = block%iec - block%isc + 1
    ny_local = block%jec - block%jsc + 1
    total_pts = sum(block%blksz)

    allocate(dest_1d(total_pts), mask(total_pts))
    allocate(source_2d(nx_local, ny_local))

    dest_1d = sentinel
    mask = 1.0_GFS_kind_phys
    do j = 1, ny_local
       do i = 1, nx_local
          source_2d(i, j) = real(i - 2 * j, GFS_kind_phys)
       end do
    end do

    call copy2block(dest_1d, source_2d, mask, validmin=validmin, validmax=validmax, flipsign=flip, block=block, rc=rc)

    test_pass = rc == ESMF_SUCCESS
    if (.not. test_pass) print *, "    FAILED: copy2block returned rc=", rc
    if (test_pass) then
       do j = block%jsc, block%jec
          do i = block%isc, block%iec
             im = packed_index(block, i, j)
             expected = real(source_2d(i - block%isc + 1, j - block%jsc + 1), GFS_kind_phys)
             if (expected <= validmin .or. expected >= validmax) then
                if (abs(dest_1d(im) - sentinel) > 1.0e-10_GFS_kind_phys) then
                   print *, "    FAILED: Out-of-range point was modified at i=", i, " j=", j
                   test_pass = .false.
                   exit
                end if
             else
                expected = -expected
                if (abs(dest_1d(im) - expected) > 1.0e-10_GFS_kind_phys) then
                   print *, "    FAILED: Flip/bounds mismatch at i=", i, " j=", j
                   test_pass = .false.
                   exit
                end if
             end if
          end do
          if (.not. test_pass) exit
       end do
    end if
    if (test_pass) then
       print *, "    PASSED"
       passed_count = passed_count + 1
    else
       print *, "    FAILED: copy2block flip/bounds verification failed"
    end if
    deallocate(dest_1d, mask, source_2d)
  end subroutine test_copy2block_flip_and_bounds

  !============================================================================
  ! TEST: copy2block mask handling
  !============================================================================
  subroutine test_copy2block_mask(block, test_num, passed_count, config_idx)
    type(block_control_type), intent(in) :: block
    integer, intent(in) :: test_num, config_idx
    integer, intent(inout) :: passed_count

    real(GFS_kind_phys), allocatable :: dest_1d(:), mask(:)
    real(GFS_kind_phys), allocatable :: source_2d(:,:)
    real(GFS_kind_phys), parameter :: sentinel = -777.0_GFS_kind_phys

    integer :: i, j, im, rc, total_pts, nx_local, ny_local
    logical :: test_pass

    print *, "  [Config ", config_idx, "] Test ", test_num, ": copy2block Mask Handling"

    call setup_copy2block_state(block)
    nx_local = block%iec - block%isc + 1
    ny_local = block%jec - block%jsc + 1
    total_pts = sum(block%blksz)

    allocate(dest_1d(total_pts), mask(total_pts))
    allocate(source_2d(nx_local, ny_local))

    dest_1d = sentinel
    mask = 0.0_GFS_kind_phys
    do j = 1, ny_local
       do i = 1, nx_local
          source_2d(i, j) = real(10 * j + i, GFS_kind_phys)
          if (mod(i + j, 2) == 0) then
             im = packed_index(block, i, j)
             mask(im) = 1.0_GFS_kind_phys
          end if
       end do
    end do

    call copy2block(dest_1d, source_2d, mask, block=block, rc=rc)

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
    real(GFS_kind_phys), allocatable :: source_2d(:,:)
    real(GFS_kind_phys), parameter :: sentinel = -999.0_GFS_kind_phys

    integer :: i, j, im, rc, total_pts, nx_local, ny_local
    logical :: test_pass

    print *, "  [Config ", config_idx, "] Test ", test_num, ": copy2block Without Mask"

    call setup_copy2block_state(block)
    nx_local = block%iec - block%isc + 1
    ny_local = block%jec - block%jsc + 1
    total_pts = sum(block%blksz)

    allocate(dest_1d(total_pts))
    allocate(source_2d(nx_local, ny_local))

    dest_1d = sentinel
    do j = 1, ny_local
       do i = 1, nx_local
          source_2d(i, j) = real(1000 + 10 * j + i, GFS_kind_phys)
       end do
    end do

    call copy2block(dest_1d, source_2d, block=block, rc=rc)

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
  ! TEST: copy2block validmin only
  !============================================================================
  subroutine test_copy2block_validmin_only(block, test_num, passed_count, config_idx)
    type(block_control_type), intent(in) :: block
    integer, intent(in) :: test_num, config_idx
    integer, intent(inout) :: passed_count

    real(GFS_kind_phys), allocatable :: dest_1d(:), mask(:)
    real(GFS_kind_phys), allocatable :: source_2d(:,:)
    real(GFS_kind_phys), parameter :: validmin = -5.0_GFS_kind_phys
    real(GFS_kind_phys), parameter :: sentinel = -444.0_GFS_kind_phys
    real(GFS_kind_phys) :: expected

    integer :: i, j, im, rc, total_pts, nx_local, ny_local
    logical :: test_pass

    print *, "  [Config ", config_idx, "] Test ", test_num, ": copy2block validmin only"

    call setup_copy2block_state(block)
    nx_local = block%iec - block%isc + 1
    ny_local = block%jec - block%jsc + 1
    total_pts = sum(block%blksz)

    allocate(dest_1d(total_pts), mask(total_pts))
    allocate(source_2d(nx_local, ny_local))

    dest_1d = sentinel
    mask = 1.0_GFS_kind_phys

    do j = 1, ny_local
       do i = 1, nx_local
          source_2d(i, j) = real(i - 3 * j, GFS_kind_phys)
       end do
    end do

    call copy2block(dest_1d, source_2d, mask, validmin=validmin, block=block, rc=rc)

    test_pass = rc == ESMF_SUCCESS
    if (.not. test_pass) print *, "    FAILED: copy2block returned rc=", rc
    if (test_pass) then
       do j = block%jsc, block%jec
          do i = block%isc, block%iec
             im = packed_index(block, i, j)
             expected = real(source_2d(i - block%isc + 1, j - block%jsc + 1), GFS_kind_phys)
             if (expected <= validmin) then
                if (abs(dest_1d(im) - sentinel) > 1.0e-10_GFS_kind_phys) then
                   print *, "    FAILED: Below-validmin point was modified at i=", i, " j=", j
                   test_pass = .false.
                   exit
                end if
             else if (abs(dest_1d(im) - expected) > 1.0e-10_GFS_kind_phys) then
                print *, "    FAILED: validmin-only mismatch at i=", i, " j=", j
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
       print *, "    FAILED: copy2block validmin-only verification failed"
    end if
    deallocate(dest_1d, mask, source_2d)
  end subroutine test_copy2block_validmin_only

  !============================================================================
  ! TEST: copy2block validmax only
  !============================================================================
  subroutine test_copy2block_validmax_only(block, test_num, passed_count, config_idx)
    type(block_control_type), intent(in) :: block
    integer, intent(in) :: test_num, config_idx
    integer, intent(inout) :: passed_count

    real(GFS_kind_phys), allocatable :: dest_1d(:), mask(:)
    real(GFS_kind_phys), allocatable :: source_2d(:,:)
    real(GFS_kind_phys), parameter :: validmax = 2.5_GFS_kind_phys
    real(GFS_kind_phys), parameter :: sentinel = -555.0_GFS_kind_phys
    real(GFS_kind_phys) :: expected

    integer :: i, j, im, rc, total_pts, nx_local, ny_local
    logical :: test_pass

    print *, "  [Config ", config_idx, "] Test ", test_num, ": copy2block validmax only"

    call setup_copy2block_state(block)
    nx_local = block%iec - block%isc + 1
    ny_local = block%jec - block%jsc + 1
    total_pts = sum(block%blksz)

    allocate(dest_1d(total_pts), mask(total_pts))
    allocate(source_2d(nx_local, ny_local))

    dest_1d = sentinel
    mask = 1.0_GFS_kind_phys

    do j = 1, ny_local
       do i = 1, nx_local
          source_2d(i, j) = real(2 * i + j, GFS_kind_phys)
       end do
    end do

    call copy2block(dest_1d, source_2d, mask, validmax=validmax, block=block, rc=rc)

    test_pass = rc == ESMF_SUCCESS
    if (.not. test_pass) print *, "    FAILED: copy2block returned rc=", rc
    if (test_pass) then
       do j = block%jsc, block%jec
          do i = block%isc, block%iec
             im = packed_index(block, i, j)
             expected = real(source_2d(i - block%isc + 1, j - block%jsc + 1), GFS_kind_phys)
             if (expected >= validmax) then
                if (abs(dest_1d(im) - sentinel) > 1.0e-10_GFS_kind_phys) then
                   print *, "    FAILED: Above-validmax point was modified at i=", i, " j=", j
                   test_pass = .false.
                   exit
                end if
             else if (abs(dest_1d(im) - expected) > 1.0e-10_GFS_kind_phys) then
                print *, "    FAILED: validmax-only mismatch at i=", i, " j=", j
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
       print *, "    FAILED: copy2block validmax-only verification failed"
    end if
    deallocate(dest_1d, mask, source_2d)
  end subroutine test_copy2block_validmax_only

  !============================================================================
  ! TEST: copy2block mixed copied/skipped values preserve destination on skip
  !============================================================================
  subroutine test_copy2block_mixed_skip_regression(block, test_num, passed_count, config_idx)
    type(block_control_type), intent(in) :: block
    integer, intent(in) :: test_num, config_idx
    integer, intent(inout) :: passed_count

    real(GFS_kind_phys), allocatable :: dest_1d(:), mask(:)
    real(GFS_kind_phys), allocatable :: source_2d(:,:)
    real(GFS_kind_phys), parameter :: validmin = 0.0_GFS_kind_phys
    real(GFS_kind_phys), parameter :: validmax = 10.0_GFS_kind_phys
    real(GFS_kind_phys), parameter :: sentinel = -654.0_GFS_kind_phys
    real(GFS_kind_phys) :: expected

    integer :: i, j, im, rc, total_pts, nx_local, ny_local
    logical :: test_pass

    print *, "  [Config ", config_idx, "] Test ", test_num, ": copy2block Mixed Skip Regression"

    call setup_copy2block_state(block)
    nx_local = block%iec - block%isc + 1
    ny_local = block%jec - block%jsc + 1
    total_pts = sum(block%blksz)

    allocate(dest_1d(total_pts), mask(total_pts))
    allocate(source_2d(nx_local, ny_local))

    dest_1d = sentinel
    mask = 1.0_GFS_kind_phys

    do j = 1, ny_local
       do i = 1, nx_local
          select case (mod(i + j, 3))
          case (0)
             source_2d(i, j) = -2.0_GFS_kind_phys
          case (1)
             source_2d(i, j) = real(i + j, GFS_kind_phys)
          case default
             source_2d(i, j) = 12.0_GFS_kind_phys
          end select
       end do
    end do

    call copy2block(dest_1d, source_2d, mask, validmin=validmin, validmax=validmax, block=block, rc=rc)

    test_pass = rc == ESMF_SUCCESS
    if (.not. test_pass) print *, "    FAILED: copy2block returned rc=", rc
    if (test_pass) then
       do j = block%jsc, block%jec
          do i = block%isc, block%iec
             im = packed_index(block, i, j)
             expected = real(source_2d(i - block%isc + 1, j - block%jsc + 1), GFS_kind_phys)
             if (expected <= validmin .or. expected >= validmax) then
                if (abs(dest_1d(im) - sentinel) > 1.0e-10_GFS_kind_phys) then
                   print *, "    FAILED: Skipped point was modified at i=", i, " j=", j
                   test_pass = .false.
                   exit
                end if
             else if (abs(dest_1d(im) - expected) > 1.0e-10_GFS_kind_phys) then
                print *, "    FAILED: In-range point was not copied at i=", i, " j=", j
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
       print *, "    FAILED: copy2block mixed-skip regression verification failed"
    end if

    deallocate(dest_1d, mask, source_2d)
  end subroutine test_copy2block_mixed_skip_regression

  !============================================================================
  ! TEST: copy2block flip only
  !============================================================================
  subroutine test_copy2block_flip_only(block, test_num, passed_count, config_idx)
    type(block_control_type), intent(in) :: block
    integer, intent(in) :: test_num, config_idx
    integer, intent(inout) :: passed_count

    real(GFS_kind_phys), allocatable :: dest_1d(:), mask(:)
    real(GFS_kind_phys), allocatable :: source_2d(:,:)
    logical, parameter :: flip = .true.
    real(GFS_kind_phys) :: expected

    integer :: i, j, im, rc, total_pts, nx_local, ny_local
    logical :: test_pass

    print *, "  [Config ", config_idx, "] Test ", test_num, ": copy2block flip only"

    call setup_copy2block_state(block)
    nx_local = block%iec - block%isc + 1
    ny_local = block%jec - block%jsc + 1
    total_pts = sum(block%blksz)

    allocate(dest_1d(total_pts), mask(total_pts))
    allocate(source_2d(nx_local, ny_local))

    dest_1d = -1.0_GFS_kind_phys
    mask = 1.0_GFS_kind_phys

    do j = 1, ny_local
       do i = 1, nx_local
          source_2d(i, j) = real(7 * i - j, GFS_kind_phys)
       end do
    end do

    call copy2block(dest_1d, source_2d, mask, flipsign=flip, block=block, rc=rc)

    test_pass = rc == ESMF_SUCCESS
    if (.not. test_pass) print *, "    FAILED: copy2block returned rc=", rc
    if (test_pass) then
       do j = block%jsc, block%jec
          do i = block%isc, block%iec
             im = packed_index(block, i, j)
             expected = -real(source_2d(i - block%isc + 1, j - block%jsc + 1), GFS_kind_phys)
             if (abs(dest_1d(im) - expected) > 1.0e-10_GFS_kind_phys) then
                print *, "    FAILED: flip-only mismatch at i=", i, " j=", j
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
       print *, "    FAILED: copy2block flip-only verification failed"
    end if
    deallocate(dest_1d, mask, source_2d)
  end subroutine test_copy2block_flip_only

  !============================================================================
  ! TEST: copy2block default flip path
  !============================================================================
  subroutine test_copy2block_default_flip(block, test_num, passed_count, config_idx)
    type(block_control_type), intent(in) :: block
    integer, intent(in) :: test_num, config_idx
    integer, intent(inout) :: passed_count

    real(GFS_kind_phys), allocatable :: dest_1d(:), mask(:)
    real(GFS_kind_phys), allocatable :: source_2d(:,:)
    real(GFS_kind_phys) :: expected

    integer :: i, j, im, rc, total_pts, nx_local, ny_local
    logical :: test_pass

    print *, "  [Config ", config_idx, "] Test ", test_num, ": copy2block default flip path"

    call setup_copy2block_state(block)
    nx_local = block%iec - block%isc + 1
    ny_local = block%jec - block%jsc + 1
    total_pts = sum(block%blksz)

    allocate(dest_1d(total_pts), mask(total_pts))
    allocate(source_2d(nx_local, ny_local))

    dest_1d = -1.0_GFS_kind_phys
    mask = 1.0_GFS_kind_phys

    do j = 1, ny_local
       do i = 1, nx_local
          source_2d(i, j) = real(i, GFS_kind_phys) / 10.0_GFS_kind_phys - real(j, GFS_kind_phys) / 7.0_GFS_kind_phys
       end do
    end do

    call copy2block(dest_1d, source_2d, mask, block=block, rc=rc)

    test_pass = rc == ESMF_SUCCESS
    if (.not. test_pass) print *, "    FAILED: copy2block returned rc=", rc
    if (test_pass) then
       do j = block%jsc, block%jec
          do i = block%isc, block%iec
             im = packed_index(block, i, j)
             expected = real(source_2d(i - block%isc + 1, j - block%jsc + 1), GFS_kind_phys)
             if (abs(dest_1d(im) - expected) > 1.0e-10_GFS_kind_phys) then
                print *, "    FAILED: default-flip mismatch at i=", i, " j=", j
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
       print *, "    FAILED: copy2block default-flip verification failed"
    end if
    deallocate(dest_1d, mask, source_2d)
  end subroutine test_copy2block_default_flip

  !============================================================================
  ! TEST: copy2block negative mask values (<= 0 must be skipped)
  !============================================================================
  subroutine test_copy2block_negative_mask_values(block, test_num, passed_count, config_idx)
    type(block_control_type), intent(in) :: block
    integer, intent(in) :: test_num, config_idx
    integer, intent(inout) :: passed_count

    real(GFS_kind_phys), allocatable :: dest_1d(:), mask(:)
    real(GFS_kind_phys), allocatable :: source_2d(:,:)
    real(GFS_kind_phys), parameter :: sentinel = -456.0_GFS_kind_phys

    integer :: i, j, im, rc, total_pts, nx_local, ny_local
    logical :: test_pass

    print *, "  [Config ", config_idx, "] Test ", test_num, ": copy2block Negative Mask Values"

    call setup_copy2block_state(block)
    nx_local = block%iec - block%isc + 1
    ny_local = block%jec - block%jsc + 1
    total_pts = sum(block%blksz)

    allocate(dest_1d(total_pts), mask(total_pts))
    allocate(source_2d(nx_local, ny_local))

    dest_1d = sentinel
    mask = 0.0_GFS_kind_phys

    do j = 1, ny_local
       do i = 1, nx_local
          source_2d(i, j) = real(20 * j + i, GFS_kind_phys)
          im = packed_index(block, i, j)
          select case (mod(i + j, 3))
          case (0)
             mask(im) = 1.0_GFS_kind_phys
          case (1)
             mask(im) = 0.0_GFS_kind_phys
          case default
             mask(im) = -1.0_GFS_kind_phys
          end select
       end do
    end do

    call copy2block(dest_1d, source_2d, mask, block=block, rc=rc)

    test_pass = rc == ESMF_SUCCESS
    if (.not. test_pass) print *, "    FAILED: copy2block returned rc=", rc
    if (test_pass) then
       do j = block%jsc, block%jec
          do i = block%isc, block%iec
             im = packed_index(block, i, j)
             if (mask(im) > 0.0_GFS_kind_phys) then
                if (abs(dest_1d(im) - real(source_2d(i - block%isc + 1, j - block%jsc + 1), GFS_kind_phys)) > 1.0e-10_GFS_kind_phys) then
                   print *, "    FAILED: Positive-mask point not copied at i=", i, " j=", j
                   test_pass = .false.
                   exit
                end if
             else if (abs(dest_1d(im) - sentinel) > 1.0e-10_GFS_kind_phys) then
                print *, "    FAILED: Non-positive-mask point was modified at i=", i, " j=", j
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
       print *, "    FAILED: copy2block negative-mask verification failed"
    end if

    deallocate(dest_1d, mask, source_2d)
  end subroutine test_copy2block_negative_mask_values

  !============================================================================
  ! TEST: copy2block mismatched source Y dimension
  !============================================================================
  subroutine test_copy2block_mismatched_shapes_y(block, test_num, passed_count, config_idx)
    type(block_control_type), intent(in) :: block
    integer, intent(in) :: test_num, config_idx
    integer, intent(inout) :: passed_count

    real(GFS_kind_phys), allocatable :: dest_1d(:), mask(:)
    real(GFS_kind_phys), allocatable :: source_2d(:,:)

    integer :: rc, total_pts, nx_local, ny_local
    logical :: test_pass

    print *, "  [Config ", config_idx, "] Test ", test_num, ": copy2block Mismatched Source Y"

    call setup_copy2block_state(block)
    nx_local = block%iec - block%isc + 1
    ny_local = block%jec - block%jsc + 1
    total_pts = sum(block%blksz)

    allocate(dest_1d(total_pts), mask(total_pts))
    allocate(source_2d(nx_local, ny_local-1))  ! Deliberate Y mismatch

    dest_1d = -1.0_GFS_kind_phys
    mask = 1.0_GFS_kind_phys
    source_2d = 1.0
    rc = -999

    call copy2block(dest_1d, source_2d, mask, block=block, rc=rc)

    test_pass = rc /= ESMF_SUCCESS
    if (test_pass) then
       print *, "    PASSED (error detected as expected)"
       passed_count = passed_count + 1
    else
       print *, "    FAILED: copy2block did not detect source Y mismatch"
    end if
    deallocate(dest_1d, mask, source_2d)
  end subroutine test_copy2block_mismatched_shapes_y

  !============================================================================
  ! TEST: copy2block oversized inputs should be accepted
  !============================================================================
  subroutine test_copy2block_oversized_inputs(block, test_num, passed_count, config_idx)
    type(block_control_type), intent(in) :: block
    integer, intent(in) :: test_num, config_idx
    integer, intent(inout) :: passed_count

    real(GFS_kind_phys), allocatable :: dest_1d(:), mask(:)
    real(GFS_kind_phys), allocatable :: source_2d(:,:)
    real(GFS_kind_phys), parameter :: sentinel = -321.0_GFS_kind_phys

    integer :: i, j, im, rc, nx_local, ny_local, total_pts
    logical :: test_pass

    print *, "  [Config ", config_idx, "] Test ", test_num, ": copy2block Oversized Inputs"

    call setup_copy2block_state(block)
    nx_local = block%iec - block%isc + 1
    ny_local = block%jec - block%jsc + 1
    total_pts = sum(block%blksz)

    allocate(dest_1d(total_pts + 5), mask(total_pts + 7), source_2d(nx_local + 2, ny_local + 3))

    dest_1d = sentinel
    mask = 1.0_GFS_kind_phys
    source_2d = -999.0_GFS_kind_phys
    do j = 1, ny_local
       do i = 1, nx_local
          source_2d(i, j) = real(3 * i + 11 * j, GFS_kind_phys)
       end do
    end do

    call copy2block(dest_1d, source_2d, mask, block=block, rc=rc)

    test_pass = rc == ESMF_SUCCESS
    if (.not. test_pass) print *, "    FAILED: copy2block returned rc=", rc
    if (test_pass) then
       do j = block%jsc, block%jec
          do i = block%isc, block%iec
             im = packed_index(block, i, j)
             if (abs(dest_1d(im) - real(source_2d(i - block%isc + 1, j - block%jsc + 1), GFS_kind_phys)) > 1.0e-10_GFS_kind_phys) then
                print *, "    FAILED: Oversized-input mapping mismatch at i=", i, " j=", j
                test_pass = .false.
                exit
             end if
          end do
          if (.not. test_pass) exit
       end do
    end if

    if (test_pass) then
       do im = total_pts + 1, size(dest_1d)
          if (abs(dest_1d(im) - sentinel) > 1.0e-10_GFS_kind_phys) then
             print *, "    FAILED: Dest tail value unexpectedly modified at im=", im
             test_pass = .false.
             exit
          end if
       end do
    end if

    if (test_pass) then
       print *, "    PASSED"
       passed_count = passed_count + 1
    else
       print *, "    FAILED: copy2block oversized-input verification failed"
    end if

    deallocate(dest_1d, mask, source_2d)
  end subroutine test_copy2block_oversized_inputs

  !============================================================================
  ! TEST: copy2block mismatched source X shape
  !============================================================================
  subroutine test_copy2block_mismatched_shapes(block, test_num, passed_count, config_idx)

    type(block_control_type), intent(in) :: block
    integer, intent(in) :: test_num, config_idx
    integer, intent(inout) :: passed_count

    real(GFS_kind_phys), allocatable :: dest_1d(:), mask(:)
    real(GFS_kind_phys), allocatable :: source_2d(:,:)

    integer :: rc, total_pts, nx_local, ny_local, bad_nx
    logical :: test_pass

    print *, "  [Config ", config_idx, "] Test ", test_num, ": copy2block Mismatched Source X"

    call setup_copy2block_state(block)
    nx_local = block%iec - block%isc + 1
    ny_local = block%jec - block%jsc + 1
    total_pts = sum(block%blksz)

    allocate(dest_1d(total_pts), mask(total_pts))
    bad_nx = merge(nx_local - 1, nx_local + 1, nx_local > 1)
    allocate(source_2d(bad_nx, ny_local))

    dest_1d = -1.0_GFS_kind_phys
    mask = 1.0_GFS_kind_phys
    source_2d = 1.0
    rc = -999

    call copy2block(dest_1d, source_2d, mask, block=block, rc=rc)

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
  ! TEST: copy2block mismatched source X/Y shape
  !============================================================================
  subroutine test_copy2block_mismatched_shapes_xy(block, test_num, passed_count, config_idx)
    type(block_control_type), intent(in) :: block
    integer, intent(in) :: test_num, config_idx
    integer, intent(inout) :: passed_count

    real(GFS_kind_phys), allocatable :: dest_1d(:), mask(:)
    real(GFS_kind_phys), allocatable :: source_2d(:,:)

    integer :: rc, total_pts, nx_local, ny_local, bad_nx, bad_ny
    logical :: test_pass

    print *, "  [Config ", config_idx, "] Test ", test_num, ": copy2block Mismatched Source XY"

    call setup_copy2block_state(block)
    nx_local = block%iec - block%isc + 1
    ny_local = block%jec - block%jsc + 1
    total_pts = sum(block%blksz)

    allocate(dest_1d(total_pts), mask(total_pts))
    bad_nx = merge(nx_local - 1, nx_local + 1, nx_local > 1)
    bad_ny = merge(ny_local - 1, ny_local + 1, ny_local > 1)
    allocate(source_2d(bad_nx, bad_ny))

    dest_1d = -1.0_GFS_kind_phys
    mask = 1.0_GFS_kind_phys
    source_2d = 1.0
    rc = -999

    call copy2block(dest_1d, source_2d, mask, block=block, rc=rc)
    test_pass = rc /= ESMF_SUCCESS
    if (test_pass) then
       print *, "    PASSED (error detected as expected)"
       passed_count = passed_count + 1
    else
       print *, "    FAILED: copy2block did not detect source XY mismatch"
    end if
    deallocate(dest_1d, mask, source_2d)
  end subroutine test_copy2block_mismatched_shapes_xy


  !============================================================================
  ! TEST: copy2block mismatched source X shape (regression guard)
  !============================================================================
  subroutine test_copy2block_mismatched_shapes_x_alt(block, test_num, passed_count, config_idx)
    type(block_control_type), intent(in) :: block
    integer, intent(in) :: test_num, config_idx
    integer, intent(inout) :: passed_count

    real(GFS_kind_phys), allocatable :: dest_1d(:), mask(:)
    real(GFS_kind_phys), allocatable :: source_2d(:,:)

    integer :: rc, nx_local, ny_local, total_pts, bad_nx
    logical :: test_pass

    print *, "  [Config ", config_idx, "] Test ", test_num, ": copy2block Mismatched Source X (Alt)"
    call setup_copy2block_state(block)

    nx_local = block%iec - block%isc + 1
    ny_local = block%jec - block%jsc + 1
    total_pts = sum(block%blksz)

    allocate(dest_1d(total_pts), mask(total_pts))
    bad_nx = merge(nx_local - 1, nx_local + 1, nx_local > 1)
    allocate(source_2d(bad_nx, ny_local))
    mask = 1.0_GFS_kind_phys
    source_2d = 1.0

    rc = -999

    call copy2block(dest_1d, source_2d, mask, block=block, rc=rc)
    test_pass = rc /= ESMF_SUCCESS
    if (test_pass) then
       print *, "    PASSED (error detected as expected)"
       passed_count = passed_count + 1
    else
       print *, "    FAILED: copy2block did not detect source X mismatch"
    end if
    deallocate(dest_1d, mask, source_2d)
  end subroutine test_copy2block_mismatched_shapes_x_alt

end program test_copy2block
