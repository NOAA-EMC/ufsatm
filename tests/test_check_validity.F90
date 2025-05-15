program test_check_validity
  use tracking_utils_mod, only: check_validity
  use mpp_mod, only: mpp_pe, mpp_init
  implicit none
  
#ifdef DEBUG
  integer, parameter :: NUM_TESTS = 10
  character(len=8), dimension(NUM_TESTS) :: test_params
  real, dimension(NUM_TESTS) :: test_values
  logical, dimension(NUM_TESTS) :: expect_valid
  integer :: i
  integer :: dummy_i = 5, dummy_j = 10
  
  ! Initialize mpp module
  call mpp_init()
  
  ! Initialize test data
  call init_test_data()
  
  ! Run tests
  do i = 1, NUM_TESTS
    call check_validity(test_params(i), test_values(i), dummy_i, dummy_j)
  end do
  
contains

  ! Initialize test cases and expected results
  subroutine init_test_data()
    ! Test case 1: Valid zeta value
    test_params(1) = "zeta"
    test_values(1) = 5.0E-3
    expect_valid(1) = .true.
    
    ! Test case 2: Invalid zeta value (too high)
    test_params(2) = "zeta"
    test_values(2) = 2.0E-2
    expect_valid(2) = .false.
    
    ! Test case 3: Invalid zeta value (too low)
    test_params(3) = "zeta"
    test_values(3) = -2.0E-2
    expect_valid(3) = .false.
    
    ! Test case 4: Valid hgt value
    test_params(4) = "hgt"
    test_values(4) = 5000.0
    expect_valid(4) = .true.
    
    ! Test case 5: Invalid hgt value (too high)
    test_params(5) = "hgt"
    test_values(5) = 2.0E4
    expect_valid(5) = .false.
    
    ! Test case 6: Valid slp value
    test_params(6) = "slp"
    test_values(6) = 1.0E5
    expect_valid(6) = .true.
    
    ! Test case 7: Invalid slp value (too low)
    test_params(7) = "slp"
    test_values(7) = 0.8E5
    expect_valid(7) = .false.
    
    ! Test case 8: Valid wind value
    test_params(8) = "wind"
    test_values(8) = 15.0
    expect_valid(8) = .true.
    
    ! Test case 9: Invalid wind value (too high)
    test_params(9) = "wind"
    test_values(9) = 250.0
    expect_valid(9) = .false.
    
    ! Test case 10: Unrecognized parameter
    test_params(10) = "unknown"
    test_values(10) = 1.0
    expect_valid(10) = .false.
  end subroutine init_test_data

#else
  ! In non-DEBUG mode, just print a message and exit
  print *, "Skipping test_check_validity - not in DEBUG mode"
#endif

end program test_check_validity
