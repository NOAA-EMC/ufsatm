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

foreach(testname test_atmos_model test_fv3_cap test_module_wrt_grid_comp test_module_copy2block)
  add_fv3atm_mpi_test(${testname})
endforeach()

foreach(testname test_post_nems_routines)
  add_fv3atm_serial_test(${testname})
endforeach()
