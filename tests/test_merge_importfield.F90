program test_merge_importfield
  !! Unit test driver for merge_importfield

  use ESMF,                only: ESMF_SUCCESS
  use block_control_mod,   only: block_control_type
  use GFS_typedefs,        only: GFS_kind_phys => kind_phys
  use atmos_model_mod,     only: merge_importfield
  use module_block_harness, only: test_config_type, setup_test_configurations, &
       initialize_block_control, cleanup_block_control, &
       setup_copy2block_state, packed_index

  implicit none

  integer, parameter :: num_configs = 8

  type(test_config_type) :: configs(num_configs)
  type(block_control_type) :: block_control

  integer :: config_idx, test_count, test_passed
  integer :: current_test_count, current_test_passed

  test_count = 0
  test_passed = 0

  print *, "=========================================="
  print *, "Unit Tests: merge_importfield"
  print *, "=========================================="
  print *, " "

  call setup_test_configurations(configs)

  do config_idx = 1, num_configs
     print *, "=========================================="
     print *, "Configuration ", config_idx, " of ", num_configs
     print *, "=========================================="

     call initialize_block_control(block_control, &
          configs(config_idx)%nx, &
          configs(config_idx)%ny, &
          configs(config_idx)%inpes, &
          configs(config_idx)%jnpes, &
          configs(config_idx)%blocksize)

     current_test_count = 0
     current_test_passed = 0

     current_test_count = current_test_count + 1
     call test_merge_field_selected_points(block_control, current_test_count, current_test_passed, config_idx)

     current_test_count = current_test_count + 1
     call test_merge_scalar_selected_points(block_control, current_test_count, current_test_passed, config_idx)

     current_test_count = current_test_count + 1
     call test_merge_all_false_noop(block_control, current_test_count, current_test_passed, config_idx)

     current_test_count = current_test_count + 1
     call test_merge_mask_blocks_updates(block_control, current_test_count, current_test_passed, config_idx)

     test_count = test_count + current_test_count
     test_passed = test_passed + current_test_passed

     print *, "Config ", config_idx, " Results: ", current_test_passed, " / ", current_test_count, " passed"
     print *, " "

     call cleanup_block_control(block_control)
  end do

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

  subroutine test_merge_field_selected_points(block, test_num, passed_count, config_idx)
    type(block_control_type), intent(in) :: block
    integer, intent(in) :: test_num, config_idx
    integer, intent(inout) :: passed_count

    real(GFS_kind_phys), allocatable :: dest_1d(:), src_1d(:), mask(:)
    real(GFS_kind_phys), allocatable :: src_2d(:,:)
    logical, allocatable :: mergeflg(:,:)

    integer :: i, j, im, rc, total_pts
    logical :: ok
    real(GFS_kind_phys), parameter :: sentinel = -777.0_GFS_kind_phys
    real(GFS_kind_phys), parameter :: sentinel2d = -888.0_GFS_kind_phys

    print *, "  [Config ", config_idx, "] Test ", test_num, ": merge_importfield field path"

    call setup_copy2block_state(block)
    total_pts = sum(block%blksz)

    allocate(dest_1d(total_pts), src_1d(total_pts), mask(total_pts))
    allocate(src_2d(block%isc:block%iec, block%jsc:block%jec))
    allocate(mergeflg(block%isc:block%iec, block%jsc:block%jec))

    dest_1d = sentinel
    src_2d = sentinel2d
    mask = 1.0_GFS_kind_phys

    do im = 1, total_pts
       src_1d(im) = real(1000 + 3 * im, GFS_kind_phys)
    end do

    do j = block%jsc, block%jec
       do i = block%isc, block%iec
          mergeflg(i, j) = mod(i + j, 2) == 0
       end do
    end do

    rc = -999
    call merge_importfield(dest_1d, src_1d, mergeflg, src_2d, mask=mask, block=block, rc=rc)

    ok = rc == ESMF_SUCCESS
    if (.not. ok) print *, "    FAILED: merge_importfield returned rc=", rc

    if (ok) then
       do j = block%jsc, block%jec
          do i = block%isc, block%iec
             im = packed_index(block, i, j)
             if (mergeflg(i, j)) then
                if (abs(dest_1d(im) - src_1d(im)) > 1.0e-10_GFS_kind_phys) then
                   print *, "    FAILED: field merge did not update dest at i=", i, " j=", j
                   ok = .false.
                   exit
                end if
                if (abs(src_2d(i, j) - real(src_1d(im), GFS_kind_phys)) > 1.0e-10_GFS_kind_phys) then
                   print *, "    FAILED: field merge did not update src_2d at i=", i, " j=", j
                   ok = .false.
                   exit
                end if
             else
                if (abs(dest_1d(im) - sentinel) > 1.0e-10_GFS_kind_phys) then
                   print *, "    FAILED: field merge changed unflagged dest at i=", i, " j=", j
                   ok = .false.
                   exit
                end if
                if (abs(src_2d(i, j) - sentinel2d) > 1.0e-10_GFS_kind_phys) then
                   print *, "    FAILED: field merge changed unflagged src_2d at i=", i, " j=", j
                   ok = .false.
                   exit
                end if
             end if
          end do
          if (.not. ok) exit
       end do
    end if

    if (ok) then
       print *, "    PASSED"
       passed_count = passed_count + 1
    else
       print *, "    FAILED: merge_importfield field-path verification failed"
    end if

    deallocate(dest_1d, src_1d, mask, src_2d, mergeflg)
  end subroutine test_merge_field_selected_points

  subroutine test_merge_scalar_selected_points(block, test_num, passed_count, config_idx)
    type(block_control_type), intent(in) :: block
    integer, intent(in) :: test_num, config_idx
    integer, intent(inout) :: passed_count

    real(GFS_kind_phys), allocatable :: dest_1d(:)
    real(GFS_kind_phys), allocatable :: src_2d(:,:)
    logical, allocatable :: mergeflg(:,:)

    integer :: i, j, im, rc, total_pts
    logical :: ok
    real(GFS_kind_phys), parameter :: scalarfill = 12.5_GFS_kind_phys
    real(GFS_kind_phys), parameter :: sentinel = -222.0_GFS_kind_phys
    real(GFS_kind_phys), parameter :: sentinel2d = -333.0_GFS_kind_phys

    print *, "  [Config ", config_idx, "] Test ", test_num, ": merge_importfield scalar path"

    call setup_copy2block_state(block)
    total_pts = sum(block%blksz)

    allocate(dest_1d(total_pts))
    allocate(src_2d(block%isc:block%iec, block%jsc:block%jec))
    allocate(mergeflg(block%isc:block%iec, block%jsc:block%jec))

    dest_1d = sentinel
    src_2d = sentinel2d

    do j = block%jsc, block%jec
       do i = block%isc, block%iec
          mergeflg(i, j) = mod(i * 2 + j, 3) == 0
       end do
    end do

    rc = -999
    call merge_importfield(dest_1d, scalarfill, mergeflg, src_2d, block=block, rc=rc)

    ok = rc == ESMF_SUCCESS
    if (.not. ok) print *, "    FAILED: merge_importfield returned rc=", rc

    if (ok) then
       do j = block%jsc, block%jec
          do i = block%isc, block%iec
             im = packed_index(block, i, j)
             if (mergeflg(i, j)) then
                if (abs(dest_1d(im) - scalarfill) > 1.0e-10_GFS_kind_phys) then
                   print *, "    FAILED: scalar merge did not update dest at i=", i, " j=", j
                   ok = .false.
                   exit
                end if
                if (abs(src_2d(i, j) - real(scalarfill, GFS_kind_phys)) > 1.0e-10_GFS_kind_phys) then
                   print *, "    FAILED: scalar merge did not update src_2d at i=", i, " j=", j
                   ok = .false.
                   exit
                end if
             else
                if (abs(dest_1d(im) - sentinel) > 1.0e-10_GFS_kind_phys) then
                   print *, "    FAILED: scalar merge changed unflagged dest at i=", i, " j=", j
                   ok = .false.
                   exit
                end if
                if (abs(src_2d(i, j) - sentinel2d) > 1.0e-10_GFS_kind_phys) then
                   print *, "    FAILED: scalar merge changed unflagged src_2d at i=", i, " j=", j
                   ok = .false.
                   exit
                end if
             end if
          end do
          if (.not. ok) exit
       end do
    end if

    if (ok) then
       print *, "    PASSED"
       passed_count = passed_count + 1
    else
       print *, "    FAILED: merge_importfield scalar-path verification failed"
    end if

    deallocate(dest_1d, src_2d, mergeflg)
  end subroutine test_merge_scalar_selected_points

  subroutine test_merge_all_false_noop(block, test_num, passed_count, config_idx)
    type(block_control_type), intent(in) :: block
    integer, intent(in) :: test_num, config_idx
    integer, intent(inout) :: passed_count

    real(GFS_kind_phys), allocatable :: dest_1d(:), src_1d(:), mask(:)
    real(GFS_kind_phys), allocatable :: src_2d(:,:)
    logical, allocatable :: mergeflg(:,:)

    integer :: i, j, im, rc, total_pts
    logical :: ok
    real(GFS_kind_phys), parameter :: sentinel = -5.0_GFS_kind_phys
    real(GFS_kind_phys), parameter :: sentinel2d = -6.0_GFS_kind_phys

    print *, "  [Config ", config_idx, "] Test ", test_num, ": merge_importfield all-false noop"

    call setup_copy2block_state(block)
    total_pts = sum(block%blksz)

    allocate(dest_1d(total_pts), src_1d(total_pts), mask(total_pts))
    allocate(src_2d(block%isc:block%iec, block%jsc:block%jec))
    allocate(mergeflg(block%isc:block%iec, block%jsc:block%jec))

    dest_1d = sentinel
    src_2d = sentinel2d
    src_1d = 44.0_GFS_kind_phys
    mergeflg = .false.
    mask = 1.0_GFS_kind_phys

    rc = -999
    call merge_importfield(dest_1d, src_1d, mergeflg, src_2d, mask=mask, block=block, rc=rc)

    ok = rc == ESMF_SUCCESS
    if (.not. ok) print *, "    FAILED: merge_importfield returned rc=", rc

    if (ok) then
       do im = 1, total_pts
          if (abs(dest_1d(im) - sentinel) > 1.0e-10_GFS_kind_phys) then
             print *, "    FAILED: all-false merge changed destination at im=", im
             ok = .false.
             exit
          end if
       end do
    end if

    if (ok) then
       do j = block%jsc, block%jec
          do i = block%isc, block%iec
             if (abs(src_2d(i, j) - sentinel2d) > 1.0e-10_GFS_kind_phys) then
                print *, "    FAILED: all-false merge changed source2d at i=", i, " j=", j
                ok = .false.
                exit
             end if
          end do
          if (.not. ok) exit
       end do
    end if

    if (ok) then
       print *, "    PASSED"
       passed_count = passed_count + 1
    else
       print *, "    FAILED: merge_importfield all-false verification failed"
    end if

    deallocate(dest_1d, src_1d, mask, src_2d, mergeflg)
  end subroutine test_merge_all_false_noop

  subroutine test_merge_mask_blocks_updates(block, test_num, passed_count, config_idx)
    type(block_control_type), intent(in) :: block
    integer, intent(in) :: test_num, config_idx
    integer, intent(inout) :: passed_count

    real(GFS_kind_phys), allocatable :: dest_1d(:)
    real(GFS_kind_phys), allocatable :: src_2d(:,:)
    real(GFS_kind_phys), allocatable :: mask(:)
    logical, allocatable :: mergeflg(:,:)

    integer :: i, j, im, rc, total_pts
    logical :: ok
    real(GFS_kind_phys), parameter :: scalarfill = 99.0_GFS_kind_phys
    real(GFS_kind_phys), parameter :: sentinel = -101.0_GFS_kind_phys
    real(GFS_kind_phys), parameter :: sentinel2d = -202.0_GFS_kind_phys

    print *, "  [Config ", config_idx, "] Test ", test_num, ": merge_importfield mask gating"

    call setup_copy2block_state(block)
    total_pts = sum(block%blksz)

    allocate(dest_1d(total_pts), mask(total_pts))
    allocate(src_2d(block%isc:block%iec, block%jsc:block%jec))
    allocate(mergeflg(block%isc:block%iec, block%jsc:block%jec))

    dest_1d = sentinel
    src_2d = sentinel2d
    mergeflg = .true.

    do j = block%jsc, block%jec
       do i = block%isc, block%iec
          im = packed_index(block, i, j)
          if (mod(i + j, 3) == 0) then
             mask(im) = 0.0_GFS_kind_phys
          else
             mask(im) = 1.0_GFS_kind_phys
          end if
       end do
    end do

    rc = -999
    call merge_importfield(dest_1d, scalarfill, mergeflg, src_2d, mask=mask, block=block, rc=rc)

    ok = rc == ESMF_SUCCESS
    if (.not. ok) print *, "    FAILED: merge_importfield returned rc=", rc

    if (ok) then
       do j = block%jsc, block%jec
          do i = block%isc, block%iec
             im = packed_index(block, i, j)
             if (mask(im) > 0.0_GFS_kind_phys) then
                if (abs(dest_1d(im) - scalarfill) > 1.0e-10_GFS_kind_phys) then
                   print *, "    FAILED: mask-gated merge missed active point at i=", i, " j=", j
                   ok = .false.
                   exit
                end if
                if (abs(src_2d(i, j) - real(scalarfill, GFS_kind_phys)) > 1.0e-10_GFS_kind_phys) then
                   print *, "    FAILED: mask-gated merge missed active source2d at i=", i, " j=", j
                   ok = .false.
                   exit
                end if
             else
                if (abs(dest_1d(im) - sentinel) > 1.0e-10_GFS_kind_phys) then
                   print *, "    FAILED: mask-gated merge updated masked-out point at i=", i, " j=", j
                   ok = .false.
                   exit
                end if
                if (abs(src_2d(i, j) - sentinel2d) > 1.0e-10_GFS_kind_phys) then
                   print *, "    FAILED: mask-gated merge updated masked-out source2d at i=", i, " j=", j
                   ok = .false.
                   exit
                end if
             end if
          end do
          if (.not. ok) exit
       end do
    end if

    if (ok) then
       print *, "    PASSED"
       passed_count = passed_count + 1
    else
       print *, "    FAILED: merge_importfield mask-gating verification failed"
    end if

    deallocate(dest_1d, mask, src_2d, mergeflg)
  end subroutine test_merge_mask_blocks_updates

end program test_merge_importfield
