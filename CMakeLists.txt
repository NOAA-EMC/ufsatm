# This is the main CMake file for ufsatm.
#
# Dusan Jovic, Alex Richert, Dustin Swales
cmake_minimum_required(VERSION 3.19)

# Handle user build options.
option(ENABLE_FV3ATM_DOCS "Enable generation of doxygen-based documentation." OFF)

# Determine whether or not to generate documentation.
if(ENABLE_FV3ATM_DOCS)
  message(STATUS "We will build the fv3atm doxygen documentation.")
  find_package(Doxygen REQUIRED)
  add_subdirectory(docs)
endif()

# Enable CI build & unit testing:
if(BUILD_TESTING)
  if (FV3)
    project(ufsatm_fv3 VERSION 1.0 LANGUAGES C CXX Fortran)
  endif()
  if (MPAS)
    project(ufsatm_mpas VERSION 1.0 LANGUAGES C CXX Fortran)
  endif()
  include(ci/CMakeLists.txt)
endif()

# Set variables for all dycore build options in UFSATM.
set(DYCORE_TARGET_MPAS ufsatm_mpas)
set(DYCORE_TARGET_FV3  ufsatm_fv3)

###############################################################################
### All dynamical cores
###############################################################################
# CCPP
set(MPI ON)
add_subdirectory(ccpp)

if(CCPP_32BIT)
  add_definitions(-DCCPP_32BIT)
else()
  remove_definitions(-DCCPP_32BIT)
endif()

# Pnetcdf
if(NOT PARALLEL_NETCDF)
  list(APPEND _ufsatm_defs_private NO_PARALLEL_NETCDF)
endif()

if(BUILD_TESTING)
  include(CTest)
  add_subdirectory(tests)
endif()

###############################################################################
### UFSATM with FV3 dynamical core
###############################################################################
if (FV3)
  set(DYCORE_TARGET ${DYCORE_TARGET_FV3})

  # These ifdefs need to be turned ON in the dycore.
  set(use_WRTCOMP ON)
  set(GFS_PHYS    ON)
  set(GFS_TYPES   ON)
  set(USE_GFSL63  ON)
  if(MOVING_NEST)
    set(MOVING_NEST ON)
  endif()

  if(HYDRO)
    set(HYDRO ON)
  endif()
  add_subdirectory(fv3/atmos_cubed_sphere)
  target_compile_definitions(fv3 PRIVATE BYPASS_BREED_SLP_INLINE)

  if(INLINE_POST)
    set(BUILD_POSTEXEC OFF)
    add_subdirectory(upp)
    set(POST_SRC io/post_nems_routines.F90 io/post_fv3.F90)
    list(APPEND _ufsatm_defs_private INLINE_POST)
  endif()

  if(MOVING_NEST)
    list(APPEND _ufsatm_defs_private MOVING_NEST)
    if(NOT HYDRO)
      list(APPEND _ufsatm_defs_private MOIST_CAPPA USE_COND)
    endif()

    if(DEBUG)
      list(APPEND _ufsatm_defs_private DEBUG)
    endif()
    if(GFS_PHYS)
      list(APPEND _ufsatm_defs_private GFS_PHYS)
    endif()
    if(GFS_TYPES)
      list(APPEND _ufsatm_defs_private GFS_TYPES)
    endif()
    if(USE_GFSL63)
      list(APPEND _ufsatm_defs_private USE_GFSL63)
    endif()
    if(INTERNAL_FILE_NML)
      list(APPEND _ufsatm_defs_private INTERNAL_FILE_NML)
    endif()
    if(ENABLE_QUAD_PRECISION)
      list(APPEND _ufsatm_defs_private ENABLE_QUAD_PRECISION)
    endif()
    if(32BIT)
      list(APPEND _ufsatm_defs_private OVERLOAD_R4 OVERLOAD_R8)
    endif()
    if(MULTI_GASES)
      list(APPEND _ufsatm_defs_private MULTI_GASES)
    endif()

    list(APPEND moving_nest_srcs
      fv3/moving_nest/bounding_box.F90
      fv3/moving_nest/fv_tracker.F90
      fv3/moving_nest/fv_moving_nest.F90
      fv3/moving_nest/fv_moving_nest_main.F90
      fv3/moving_nest/fv_moving_nest_physics.F90
      fv3/moving_nest/fv_moving_nest_types.F90
      fv3/moving_nest/fv_moving_nest_utils.F90
      )
  else()
    list(APPEND moving_nest_srcs "")
  endif()

  # FV3 drivers and dependencies
  add_library(${DYCORE_TARGET}
    fv3/atmos_model.F90
    fv3/fv3_cap.F90
    fv3/module_fv3_config.F90
    fv3/module_fcst_grid_comp.F90
    stochastic_physics/stochastic_physics_wrapper.F90
    cpl/module_block_data.F90
    cpl/module_cplfields.F90
    cpl/module_cap_cpl.F90
    cpl/module_cplscalars.F90
    io/fv3atm_common_io.F90
    io/fv3atm_clm_lake_io.F90
    io/fv3atm_rrfs_sd_io.F90
    io/fv3atm_sfc_io.F90
    io/fv3atm_oro_io.F90
    io/fv3atm_history_io.F90
    io/fv3atm_restart_io.F90
    io/module_write_netcdf.F90
    io/module_write_restart_netcdf.F90
    io/module_fv3_io_def.F90
    io/module_write_internal_state.F90
    io/module_wrt_grid_comp.F90
    ${moving_nest_srcs}
    ${POST_SRC}
    )
  add_dependencies(${DYCORE_TARGET} fv3 fv3ccpp stochastic_physics)

  list(APPEND _ufsatm_defs_private GFS_PHYS
                                   INTERNAL_FILE_NML
                                   use_WRTCOMP)
endif()

###############################################################################
### UFSATM with MPAS dynamical core.
###############################################################################
if (MPAS)
  set(DYCORE_TARGET ${DYCORE_TARGET_MPAS})

  # Include MPAS Cmake tools.
  include(${CMAKE_CURRENT_SOURCE_DIR}/mpas/MPAS-Model/cmake/Functions/MPAS_Functions.cmake)

  # Set any pre-processor directive needed in MPAS dycore.
  get_mpas_version(MPAS_VERSION)
  set(MPAS_ALL_CORES atmosphere)
  set(MPAS_CORES atmosphere CACHE STRING "MPAS cores to build. Options: ${MPAS_ALL_CORES}")
  if(MPAS_CORES MATCHES " ") #Convert strings separated with spaces to CMake list separated with ';'
    string(REPLACE " " ";" MPAS_CORES ${MPAS_CORES})
    set(MPAS_CORES ${MPAS_CORES} CACHE STRING "MPAS cores to build. Options: ${MPAS_ALL_CORES}" FORCE)
  endif()
  #option(MPAS_CAM_DYCORE TRUE)
  add_definitions(-DMPAS_USE_MPI_F08)
  add_definitions(-DMPAS_PIO_SUPPORT)
  add_definitions(-DMPAS_CAM_DYCORE)
  add_definitions(-DMPAS_UFS_DYCORE)

  # Source files for MPAS dynamical core drivers.
  set(MPAS_MAIN_SRC      ${CMAKE_CURRENT_SOURCE_DIR}/mpas/MPAS-Model/src/driver/mpas.F)
  set(MPAS_SUBDRIVER_SRC ${CMAKE_CURRENT_SOURCE_DIR}/mpas/MPAS-Model/src/driver/mpas_subdriver.F)

  # MPAS dynamical core
  add_subdirectory(mpas)

  # MPAS drivers and dependencies
  add_library(${DYCORE_TARGET}
    mpas/atmos_model.F90
    mpas/mpas_cap.F90
    mpas/ufs_mpas_subdriver.F90
    mpas/module_mpas_config.F90
    mpas/module_fcst_grid_comp.F90
    ccpp/data/MPAS_typedefs.F90
    ccpp/driver/MPAS_init.F90
    )
  add_dependencies(${DYCORE_TARGET} mpas mpasccpp)

  if(NOT MPAS_GIT_VERSION)
    find_package(Git QUIET)
    if(GIT_FOUND)
      execute_process(COMMAND ${GIT_EXECUTABLE} describe --dirty
                      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/mpas/MPAS-Model"
		      OUTPUT_VARIABLE _mpas_git_version
		      ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
    else()
      set(_mpas_git_version "Unknown")
    endif()
    set(MPAS_GIT_VERSION ${_mpas_git_version} CACHE STRING "MPAS-Model git version")
    message(STATUS "Setting MPAS_GIT_VERSION ${_mpas_git_version}")
  endif()
endif()

###############################################################################
### Link libraries
###############################################################################
target_compile_definitions(${DYCORE_TARGET} PRIVATE "${_ufsatm_defs_private}")

set_target_properties(${DYCORE_TARGET} PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
target_include_directories(${DYCORE_TARGET} INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>
                                                  $<INSTALL_INTERFACE:mod>)

if (FV3)
  target_link_libraries(${DYCORE_TARGET} PUBLIC fv3
                                         fv3ccpp
					 stochastic_physics
					 fms)
  if(INLINE_POST)
    target_link_libraries(${DYCORE_TARGET} PUBLIC upp::upp)
  endif()				       
endif()

if (MPAS)
  target_link_libraries(${DYCORE_TARGET} PUBLIC mpas
                                         mpasccpp
					 fms)
endif()

# Always include EMC libraries in dycore install
target_link_libraries(${DYCORE_TARGET} PUBLIC w3emc::w3emc_d
                                              sp::sp_d
					      bacio::bacio_4
					      esmf)
# OpenMP
if(OPENMP)
  target_link_libraries(${DYCORE_TARGET} PUBLIC OpenMP::OpenMP_Fortran)
endif()

###############################################################################
### Install UFSATM
###############################################################################
install(
  TARGETS ${DYCORE_TARGET}
  EXPORT ufsatm-config
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod DESTINATION ${CMAKE_INSTALL_PREFIX})

install(EXPORT ufsatm-config
  DESTINATION lib/cmake)
