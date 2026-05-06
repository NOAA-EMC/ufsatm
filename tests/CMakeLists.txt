# This file sets up CTest unit tests for fv3atm.
#
# Alex Richert, Jan 2024
#
# Add More unit tests
# Daniel Sarmiento, Jun 2025

find_package(PIO REQUIRED COMPONENTS C Fortran)

# Stage test data
execute_process(COMMAND cmake -E create_symlink
  "${CMAKE_CURRENT_SOURCE_DIR}/data"
  "${CMAKE_CURRENT_BINARY_DIR}/data"
  )

function(add_fv3atm_mpi_test TESTNAME)
  add_executable(${TESTNAME} ${TESTNAME}.F90)
  target_link_libraries(${TESTNAME} PRIVATE ufsatm_fv3 MPI::MPI_Fortran PIO::PIO_C)
  add_test(${TESTNAME} ${MPIEXEC_EXECUTABLE} -n 2 ${CMAKE_CURRENT_BINARY_DIR}/${TESTNAME})
endfunction()

function(add_fv3atm_serial_test TESTNAME)
  add_executable(${TESTNAME} ${TESTNAME}.F90)
  target_link_libraries(${TESTNAME} PRIVATE ufsatm_fv3)
  add_test(${TESTNAME} ${CMAKE_CURRENT_BINARY_DIR}/${TESTNAME})
endfunction()

# Build the shared test harness module once to avoid parallel .mod write races.
add_library(test_block_harness STATIC module_block_harness.F90)
set_target_properties(test_block_harness PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod/test_block_harness)
target_include_directories(test_block_harness PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}/mod/test_block_harness)
target_link_libraries(test_block_harness PRIVATE ufsatm_fv3 MPI::MPI_Fortran PIO::PIO_C)

foreach(testname test_atmos_model test_fv3_cap test_module_wrt_grid_comp test_copy2block)
  add_fv3atm_mpi_test(${testname})
endforeach()

target_link_libraries(test_copy2block PRIVATE test_block_harness)

add_executable(test_merge_importfield test_merge_importfield.F90)
target_link_libraries(test_merge_importfield PRIVATE ufsatm_fv3 MPI::MPI_Fortran PIO::PIO_C test_block_harness)
add_test(test_merge_importfield ${MPIEXEC_EXECUTABLE} -n 2 ${CMAKE_CURRENT_BINARY_DIR}/test_merge_importfield)

foreach(testname test_post_nems_routines)
  add_fv3atm_serial_test(${testname})
endforeach()
