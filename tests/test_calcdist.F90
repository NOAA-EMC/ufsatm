program test_calcdist
  use fv_tracker_mod, only: calcdist
  implicit none
  
  integer, parameter :: NUM_TESTS = 3
  real, dimension(NUM_TESTS) :: lon1, lat1, lon2, lat2
  real, dimension(NUM_TESTS) :: expected_dist, expected_degrees
  real :: xdist, degrees
  integer :: i
  real :: margin_error = 1e-3
  
  ! Initialize test data
  call init_test_data()
  
  ! Run tests
  do i = 1, NUM_TESTS
    degrees = 0.0
    xdist = 0.0
    
    call calcdist(lon1(i), lat1(i), lon2(i), lat2(i), xdist, degrees)
    
    ! Check results
    if (abs(xdist - expected_dist(i)) > margin_error) then
      stop i
    end if
  end do
  
contains

  subroutine init_test_data()
    
    ! Test 1: Zero distance
    lon1(1) = 100.0
    lat1(1) = 10.0
    lon2(1) = 100.0
    lat2(1) = 10.0
    expected_degrees(1) = 0.0
    expected_dist(1) = 0.0
    
    ! Test 2: Miami to Berlin
    lon1(2) = -80.19
    lat1(2) = 25.76
    lon2(2) = 13.41
    lat2(2) = 52.52
    expected_degrees(2) = 0.0
    expected_dist(2) = 7996.26 ! 8010.95km
    
    ! Test 3: Helsinki to Cape Town 
    lon1(3) = 24.94
    lat1(3) = 60.17
    lon2(3) = 18.42
    lat2(3) = -33.90
    expected_degrees(3) = 0.0
    expected_dist(3) = 10477.17 ! 10442.76km
    
  end subroutine init_test_data

end program test_calcdist
