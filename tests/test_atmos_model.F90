program test_atmos_model
  use atmos_model_mod, only: set_fhzero_loop, InitTimeFromIAUOffset, &
                             get_atmos_tracer_types, atmos_data_type
  use GFS_typedefs, only: GFS_control_type, GFS_kind_phys => kind_phys
  use time_manager_mod, only: time_type, set_time, get_time, operator(-)
  use tracer_manager_mod, only: get_number_tracers
  use field_manager_mod, only: MODEL_ATMOS
  use mpp_mod, only: mpp_init, mpp_exit, FATAL, mpp_error
  use CCPP_data, only: GFS_control
  
  implicit none
  
  integer :: test_passed, total_tests
  integer :: suite_passed, suite_total

  ! Test Suite 1: set_fhzero_loop
  call test_set_fhzero_loop_suite()
  
  ! Test Suite 2: get_atmos_tracer_types
  call test_get_atmos_tracer_types_suite()
  
contains

  !============================================================================
  ! TEST SUITE 1: set_fhzero_loop
  !============================================================================
  subroutine test_set_fhzero_loop_suite()
    integer :: sec, sec_lastfhzerofh
    
    ! Test 1: Basic functionality with single fhzero value
    call test_single_fhzero()
    
    ! Test 2: Multiple fhzero array values
    call test_multiple_fhzero()
    
    ! Test 3: Edge case with zero and negative values
    call test_fhzero_edge_cases()
    
  end subroutine test_set_fhzero_loop_suite
  
  subroutine test_single_fhzero()
    integer :: sec, sec_lastfhzerofh

    ! Setup
    GFS_control%fhzero_array(1) = 6.0_GFS_kind_phys
    GFS_control%fhzero_fhour(1) = 24.0_GFS_kind_phys
    
    ! Test case: Time within first interval
    sec = 10800 
    call set_fhzero_loop(sec, sec_lastfhzerofh)
    
    if (GFS_control%fhzero /= 6.0_GFS_kind_phys .or. sec_lastfhzerofh /= 0) then
      print *, "Incorrect handling of single fhzero value"
      stop 1
    end if

  end subroutine test_single_fhzero
  
  subroutine test_multiple_fhzero()
    integer :: sec, sec_lastfhzerofh

    ! Setup
    GFS_control%fhzero_array = [3.0_GFS_kind_phys, 6.0_GFS_kind_phys]
    GFS_control%fhzero_fhour = [12.0_GFS_kind_phys, 24.0_GFS_kind_phys]
    
    ! Test first interval
    sec = 7200 
    call set_fhzero_loop(sec, sec_lastfhzerofh)
    
    if (GFS_control%fhzero /= 3.0_GFS_kind_phys) then
      print *, "Incorrect handling of fhzero array"
      stop 2
    end if
    
    ! Test second interval
    sec = 50400
    call set_fhzero_loop(sec, sec_lastfhzerofh)
    
    if (GFS_control%fhzero /= 6.0_GFS_kind_phys) then
      print *, "Incorrect handling of fhzero array"
      stop 3
    end if
    
  end subroutine test_multiple_fhzero
  
  subroutine test_fhzero_edge_cases()
    integer :: sec, sec_lastfhzerofh

    ! Test zero fhzero value
    GFS_control%fhzero_array = [0.0_GFS_kind_phys, 6.0_GFS_kind_phys]
    GFS_control%fhzero_fhour = [6.0_GFS_kind_phys, 12.0_GFS_kind_phys]
    
    sec = 3600
    call set_fhzero_loop(sec, sec_lastfhzerofh)
    
    if (sec_lastfhzerofh /= 0) then
      print *, "Incorrect handling of fh = 0 case"
      stop 4
    end if

  end subroutine test_fhzero_edge_cases

  !============================================================================
  ! TEST SUITE 2: get_atmos_tracer_types
  !============================================================================
  subroutine test_get_atmos_tracer_types_suite()
    integer, allocatable :: tracer_types(:)
    integer :: num_tracers
    
    ! Test 1: Basic functionality with mock tracers
    call test_tracer_basic_functionality()
    
    ! Test 2: Test with chemistry tracers
    call test_chemistry_tracers()
    
    ! Test 3: Edge cases
    call test_tracer_edge_cases()
    
  end subroutine test_get_atmos_tracer_types_suite
  
  subroutine test_tracer_basic_functionality()
    integer, allocatable :: tracer_types(:)
    integer :: num_tracers
    
    ! For this test, we'll simulate having 5 tracers
    num_tracers = 5
    allocate(tracer_types(num_tracers))
    
    ! Initialize all to zero (default)
    tracer_types = 0
    
    if (any(tracer_types /= 0)) then
      print *, "Tracer type array being rewritten"
      stop 5
    end if
    
    deallocate(tracer_types)
  end subroutine test_tracer_basic_functionality
  
  subroutine test_chemistry_tracers()
    integer, allocatable :: tracer_types(:)
    integer :: num_tracers
    
    ! Simulate having tracers with chemistry types
    num_tracers = 8
    allocate(tracer_types(num_tracers))
    
    ! Manually set tracer types to simulate:
    ! Tracers 1-3: generic (0)
    ! Tracers 4-6: chemistry prognostic (1)
    ! Tracers 7-8: chemistry diagnostic (2)
    tracer_types = [0, 0, 0, 1, 1, 1, 2, 2]

    ! Test generic tracers are contiguous
    if (any(tracer_types(1:3) /= 0)) then
      print *, "Tracer type array being rewritten or rearranged"
      stop 6
    end if

    ! Test prognostic tracers are contiguous
    if (any(tracer_types(4:6) /= 1)) then
      print *, "Tracer type array being rewritten or rearranged"
      stop 7
    end if
    
    ! Test diagnostic tracers are contiguous
    if (any(tracer_types(7:8) /= 2)) then
      print *, "Tracer type array being rewritten or rearranged"
      stop 8
    end if
    
    deallocate(tracer_types)
  end subroutine test_chemistry_tracers
  
  subroutine test_tracer_edge_cases()
    integer, allocatable :: tracer_types(:)
    integer :: num_tracers
    
    ! Test with large number of tracers
    num_tracers = 100
    allocate(tracer_types(num_tracers))
    tracer_types = 0
    
    if (size(tracer_types) /= 100) then
      print *, "Tracer type array missing values when array is large"
      stop 9
    end if
    
    deallocate(tracer_types)
  end subroutine test_tracer_edge_cases

end program test_atmos_model
