program test_output_hours
  use fv3atm_cap_mod, only: OutputHours_FrequencyInput, OutputHours_ArrayInput
  use module_fv3_config, only: dt_atmos, output_fh
  use module_fv3_io_def, only: lflname_fulltime

  implicit none
  
  ! Test variables
  integer :: test_passed, test_failed
  integer :: i, expected_size
  logical :: test_result
  
  ! Variables for testing
  real :: nfhmax, output_startfh, outputfh2(2)
  integer :: noutput_fh
  
  dt_atmos = 1800
  
  call test_frequency_input()
  call test_array_input()
  
contains

  ! Test OutputHours_FrequencyInput subroutine
  subroutine test_frequency_input()
    
    !============================================
    ! Test 1: Basic frequency input with start at 0
    nfhmax = 24.0
    output_startfh = 0.0
    outputfh2(1) = 3.0
    outputfh2(2) = -1.0
    lflname_fulltime = .false.
    
    if (allocated(output_fh)) deallocate(output_fh)
    call OutputHours_FrequencyInput(nfhmax, output_startfh, outputfh2)
    
    expected_size = 9
    
    if (size(output_fh) /= expected_size) then
    !Expected size should be 9
      stop 1
    end if
    if (abs(output_fh(1) - dt_atmos/3600.0) > 1e-6) then
      !First output should be dt_atmos/3600
      stop 2
    end if
    ! lflname_fulltime should be false
    if (lflname_fulltime) then
      !Excluding first element, output_fh are all integers
      stop 3
    end if
    
    !============================================
    ! Test 2: Frequency input with non-zero start
    nfhmax = 48.0
    output_startfh = 6.0
    outputfh2(1) = 6.0
    outputfh2(2) = -1.0
    lflname_fulltime = .false.
    
    if (allocated(output_fh)) deallocate(output_fh)
    call OutputHours_FrequencyInput(nfhmax, output_startfh, outputfh2)
    
    expected_size = 8
    
    !Expected size should be 8
    if (size(output_fh) /= expected_size) then
      stop 4
    end if
    
    ! First value (should be 6.0)
    if (abs(output_fh(1) - 6.0) > 1e-6) then
      stop 5
    end if
      
    ! lflname_fulltime should be false
    if (lflname_fulltime) then
      !Excluding first element, output_fh are all integers
      stop 6
    end if
    
    !============================================
    ! Test 3: Frequency that creates non-integer hours
    ! Only checking lflname_fulltime since other aspects of
    !     array generation were already checked.
    nfhmax = 10.0
    output_startfh = 0.0
    outputfh2(1) = 2.5
    outputfh2(2) = -1.0
    lflname_fulltime = .false.
    
    if (allocated(output_fh)) deallocate(output_fh)
    call OutputHours_FrequencyInput(nfhmax, output_startfh, outputfh2)
 
    ! Check lflname_fulltime, should be True since non-integer values exist
    if (.not. lflname_fulltime) then
      stop 7
    end if

    !============================================
    ! Test 4: nfhmax equals output_startfh
    nfhmax = 6.0
    output_startfh = 6.0
    outputfh2(1) = 3.0
    outputfh2(2) = -1.0
    
    if (allocated(output_fh)) deallocate(output_fh)
    call OutputHours_FrequencyInput(nfhmax, output_startfh, outputfh2)

    ! output_fh should not allocate when nfhmax == output_startfh
    if (allocated(output_fh)) then
      stop 8
    end if
    
  end subroutine test_frequency_input
  
  subroutine test_array_input()
    
    !============================================
    ! Test 1: Basic array input with start at 0
    noutput_fh = 5
    output_startfh = 0.0
    lflname_fulltime = .false.
    
    if (allocated(output_fh)) deallocate(output_fh)
    allocate(output_fh(noutput_fh))
    output_fh = (/ 0.0, 3.0, 6.0, 9.0, 12.0 /)
    
    call OutputHours_ArrayInput(noutput_fh, output_startfh)
    
    ! Check first value. Should be dt_atmos/3600
    if (abs(output_fh(1) - dt_atmos/3600.0) > 1e-6) then
      stop 9
    end if
    
    ! Check lflname_fulltime, should be false
    !Excluding first element, output_fh are all integers
    if (lflname_fulltime) then
      stop 10
    end if
    
    !============================================
    ! Test 2: Array input with non-zero start
    noutput_fh = 4
    output_startfh = 6.0
    lflname_fulltime = .false.
    
    if (allocated(output_fh)) deallocate(output_fh)
    allocate(output_fh(noutput_fh))
    output_fh = (/ 0.0, 6.0, 12.0, 18.0 /)
    
    call OutputHours_ArrayInput(noutput_fh, output_startfh)
    
    ! Check values (should be shifted by 6: 6, 12, 18, 24)
    if (abs(output_fh(1) - 6.0) > 1e-6 .or. &
        abs(output_fh(2) - 12.0) > 1e-6 .or. &
        abs(output_fh(3) - 18.0) > 1e-6 .or. &
        abs(output_fh(4) - 24.0) > 1e-6) then
      stop 11
    end if
    
    ! Check lflname_fulltime (should be false)
    if (lflname_fulltime) then
      stop 12
    end if
    
    !============================================
    ! Test 3: Array with non-integer hours
    noutput_fh = 4
    output_startfh = 0.0
    lflname_fulltime = .false.
    
    if (allocated(output_fh)) deallocate(output_fh)
    allocate(output_fh(noutput_fh))
    output_fh = (/ 1.5, 3.0, 4.5, 6.0 /)
    
    call OutputHours_ArrayInput(noutput_fh, output_startfh)
    
    test_result = .true.
    print *, '  Output hours:', output_fh
    
    ! Check lflname_fulltime (should be true)
    if (.not. lflname_fulltime) then
      stop 13
    else
      print *, '  lflname_fulltime correctly set to true'
    end if
    
    !============================================
    ! Test 4: Array with fractional hours from non-zero start
    noutput_fh = 3
    output_startfh = 0.5
    lflname_fulltime = .false.
    
    if (allocated(output_fh)) deallocate(output_fh)
    allocate(output_fh(noutput_fh))
    output_fh = (/ 0.0, 3.0, 6.0 /)
    
    call OutputHours_ArrayInput(noutput_fh, output_startfh)
    
    ! Check values (should be 0.5, 3.5, 6.5)
    if (abs(output_fh(1) - 0.5) > 1e-6 .or. &
        abs(output_fh(2) - 3.5) > 1e-6 .or. &
        abs(output_fh(3) - 6.5) > 1e-6) then
      stop 14
    end if
    
    ! Check lflname_fulltime (should be true)
    if (.not. lflname_fulltime) then
      stop 15
    end if
    
  end subroutine test_array_input
  
end program test_output_hours
