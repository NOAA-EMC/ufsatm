program test_clean_lon_lat
  use fv_tracker
  implicit none
  
  integer, parameter :: NUM_TESTS = 11
  real, dimension(NUM_TESTS) :: test_lons, test_lats
  real, dimension(NUM_TESTS) :: expected_lons, expected_lats
  real :: lon, lat
  integer :: i
  
  ! Initialize test data
  call init_test_data()
  
  do i = 1, NUM_TESTS
    lon = test_lons(i)
    lat = test_lats(i)
    
    call clean_lon_lat(lon, lat)

    ! Check results    
    if (abs(lon - expected_lons(i)) > 10e-6 .or. abs(lat - expected_lats(i)) > 10e-6) then
      write(*,'(A,I2,A,F12.6,A,F12.6,A)') "Test ", i, ": clean_lon_lat(", lon, ", ", lat, ")" 
      stop i
    end if
      write(*,'(A,F12.6,A,F12.6,I2)') "Passed with ", lon, ", ", lat, ")", i
  end do
  
contains

  subroutine init_test_data()
    ! Test 1: No changes needed
    test_lons(1) = 45.0
    test_lats(1) = 30.0
    expected_lons(1) = 45.0
    expected_lats(1) = 30.0
    
    ! Test 2: Longitude is > 180
    test_lons(2) = 270.0
    test_lats(2) = 45.0
    expected_lons(2) = -90.0
    expected_lats(2) = 45.0
    
    ! Test 3: Longitude is < -180
    test_lons(3) = -270.0
    test_lats(3) = 45.0
    expected_lons(3) = 90.0
    expected_lats(3) = 45.0
    
    ! Test 4: Latitude is > 90
    test_lons(4) = 45.0
    test_lats(4) = 100.0
    expected_lons(4) = 45.0
    expected_lats(4) = -80.0    
    
    ! Test case 5: Latitude is < -90
    test_lons(5) = 45.0
    test_lats(5) = -100.0
    expected_lons(5) = 45.0
    expected_lats(5) = 80.0
    
    ! Test case 6: Latitude and longitude are > 90 and > 180
    test_lons(6) = 225.0
    test_lats(6) = 95.0
    expected_lons(6) = -135.0
    expected_lats(6) = -85.0
    
    ! Test case 7: Latitude and longitude are < 90 and < 180
    test_lons(7) = -200.0
    test_lats(7) = -105.0
    expected_lons(7) = 160.0
    expected_lats(7) = 75.0
    
    ! Test case 8: Test zeros
    test_lons(8) = 0.0
    test_lats(8) = 0.0      
    expected_lons(8) = 0.0
    expected_lats(8) = 0.0
    
    ! Test case 9: Test edges of range (positive)
    test_lons(9) = 180.0
    test_lats(9) = 90.0
    expected_lons(9) = 180.0
    expected_lats(9) = 90.0

    ! Test case 10: Test edges of range (negative)
    test_lons(10) = -180.0
    test_lats(10) = -90.0
    expected_lons(10) = -180.0
    expected_lats(10) = -90.0

    ! Test case 11: Large values
    test_lons(11) = 5000.0
    test_lats(11) = -990.0
    expected_lons(11) = 140.0
    expected_lats(11) = -90.0
  end subroutine init_test_data

end program test_clean_lon_lat
