cmake_minimum_required(VERSION 3.0)

project(CCPP-FV3
        LANGUAGES C CXX Fortran)
set(PROJECT "CCPP-FV3")

#------------------------------------------------------------------------------
# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to 'Release' as none was specified.")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release")
endif()

#------------------------------------------------------------------------------
# Call to CCPP code generator
if(DEBUG)
  # Enable debugging features in auto-generated physics caps
  set(_ccpp_debug_arg "--debug")
  # Enable verbose output from ccpp_prebuild.py
  set(_ccpp_verbose_arg "--verbose")
endif()
if(DEFINED CCPP_SUITES)
  set(_ccpp_suites_arg "--suites=${CCPP_SUITES}")
  message("Calling CCPP code generator (ccpp_prebuild.py) for suites ${_ccpp_suites_arg} ...")
else()
  message("Calling CCPP code generator (ccpp_prebuild.py) for all available suites ...")
endif()
execute_process(COMMAND ${Python_EXECUTABLE}
                        "framework/scripts/ccpp_prebuild.py"
                        "--config=config/ccpp_prebuild_config.py"
                        "--builddir=${CMAKE_CURRENT_BINARY_DIR}" ${_ccpp_suites_arg} ${_ccpp_debug_arg} ${_ccpp_verbose_arg}
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                OUTPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/ccpp_prebuild.out
                ERROR_FILE  ${CMAKE_CURRENT_BINARY_DIR}/ccpp_prebuild.err
                RESULT_VARIABLE RC)
# Check return code from ccpp_prebuild.py
if(NOT RC EQUAL 0)
  message(FATAL_ERROR "An error occured while running ccpp_prebuild.py, check ${CMAKE_CURRENT_BINARY_DIR}/ccpp_prebuild.{out,err}")
endif()

#------------------------------------------------------------------------------
# Add host-model specific preprocessor flag (needed for some physics schemes)
add_definitions(-DFV3)

#------------------------------------------------------------------------------
# Set MPI flags for C/C++/Fortran preprocessor
if(MPI)
  add_definitions(-DMPI)
endif()

#------------------------------------------------------------------------------
# Set flag for 32bit dynamics build
if(32BIT)
  message(STATUS "Compile CCPP fast physics with 32-bit precision")
  add_definitions(-DOVERLOAD_R4)
  set(CMAKE_Fortran_FLAGS_DYNAMICS "")
else()
  message(STATUS "Compile CCPP fast physics with 64-bit precision")
  remove_definitions(-DOVERLOAD_R8)
  remove_definitions(-DOVERLOAD_R4)
endif()

if(CCPP_32BIT)
  message(STATUS "Compile CCPP slow physics with 32-bit precision")
  add_definitions(-DSINGLE_PREC)
  add_definitions(-DRTE_USE_SP)
  if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(CMAKE_Fortran_FLAGS_PHYSICS  "-real-size 32")
  elseif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(CMAKE_Fortran_FLAGS_PHYSICS  "-fno-default-real-8 -fdefault-double-8")
  endif()
else(CCPP_32BIT)
  message(STATUS "Compile CCPP slow physics with 64-bit precision")
  remove_definitions(-DSINGLE_PREC)
  remove_definitions(-DRTE_USE_SP)
  if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(CMAKE_Fortran_FLAGS_PHYSICS  "-real-size 64")
  elseif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(CMAKE_Fortran_FLAGS_PHYSICS  "-fdefault-real-8 -fdefault-double-8")
  endif()
endif(CCPP_32BIT)

#------------------------------------------------------------------------------
# Add model-specific flags for C/C++/Fortran preprocessor
if(NOT HYDRO)
  add_definitions(-DMOIST_CAPPA -DUSE_COND -DNEMS_GSM)
endif()
add_definitions(-DINTERNAL_FILE_NML -DNEMS_GSM)

if(MULTI_GASES)
  add_definitions(-DMULTI_GASES)
endif()

if(IDEA_PHYS)
  add_definitions(-DIDEA_PHYS)
endif()
#------------------------------------------------------------------------------
# Build CCPP framework and physics

add_subdirectory(framework)
add_subdirectory(physics)

#------------------------------------------------------------------------------
# Build fv3ccpp

# Can we move data/*.F90 and driver/*.F90 to this directory ???
add_library(
    fv3ccpp

    driver/GFS_diagnostics.F90
    driver/GFS_restart.F90
    driver/GFS_init.F90
    driver/CCPP_driver.F90

    ${CMAKE_CURRENT_BINARY_DIR}/physics/ccpp_static_api.F90
)

# Compile GFS_diagnostics.F90 without optimization, this leads to out of memory errors on wcoss_dell_p3
set_property(SOURCE driver/GFS_diagnostics.F90 APPEND_STRING PROPERTY COMPILE_FLAGS "-O0")

target_link_libraries(fv3ccpp PUBLIC ccpp_framework)
target_link_libraries(fv3ccpp PUBLIC ccpp_physics)

if(OPENMP)
  target_link_libraries(fv3ccpp PUBLIC OpenMP::OpenMP_Fortran)
endif()

set_target_properties(fv3ccpp PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
target_include_directories(fv3ccpp PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>)

###############################################################################
### Install
###############################################################################
install(
  TARGETS fv3ccpp ccpp_framework ccpp_physics
  EXPORT fv3ccpp-config
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod DESTINATION ${CMAKE_INSTALL_PREFIX})

install(EXPORT fv3ccpp-config
  DESTINATION lib/cmake)
