# This is the main CMake file for fv3atm.
#
# Dusan Jovic, Alex Richert
cmake_minimum_required(VERSION 3.19)

# Handle user build options.
option(ENABLE_DOCS "Enable generation of doxygen-based documentation." OFF)

# Determine whether or not to generate documentation.
if(ENABLE_DOCS)
  find_package(Doxygen REQUIRED)
  add_subdirectory(docs)  
endif()

# Enable CI build & unit testing:
if(BUILD_TESTING)
  project(fv3atm VERSION 1.0 LANGUAGES C CXX Fortran)
  include(ci/CMakeLists.txt)
endif()

###############################################################################
### CCPP
###############################################################################

add_subdirectory(ccpp)

###############################################################################
### fv3 dynamical core
###############################################################################

# These ifdefs need to be turned ON in the dycore.
set(use_WRTCOMP ON)
set(GFS_PHYS    ON)
set(GFS_TYPES   ON)
set(USE_GFSL63  ON)
if(MOVING_NEST)
  set(MOVING_NEST ON)
endif()
add_subdirectory(atmos_cubed_sphere)

###############################################################################
### fv3atm
###############################################################################
if(INLINE_POST)
  set(BUILD_POSTEXEC OFF)
  add_subdirectory(upp)
  set(POST_SRC io/post_nems_routines.F90 io/post_fv3.F90)
  list(APPEND _fv3atm_defs_private INLINE_POST)
endif()

if(CCPP_32BIT)
  add_definitions(-DCCPP_32BIT)
else()
  remove_definitions(-DCCPP_32BIT)
endif()

if(NOT PARALLEL_NETCDF)
  list(APPEND _fv3atm_defs_private NO_PARALLEL_NETCDF)
endif()

if(MOVING_NEST)
  list(APPEND _fv3atm_defs_private MOVING_NEST MOIST_CAPPA USE_COND)
  if(DEBUG)
    list(APPEND _fv3atm_defs_private DEBUG)
  endif()
  if(GFS_PHYS)
    list(APPEND _fv3atm_defs_private GFS_PHYS)
  endif()
  if(GFS_TYPES)
    list(APPEND _fv3atm_defs_private GFS_TYPES)
  endif()
  if(USE_GFSL63)
    list(APPEND _fv3atm_defs_private USE_GFSL63)
  endif()
  if(INTERNAL_FILE_NML)
    list(APPEND _fv3atm_defs_private INTERNAL_FILE_NML)
  endif()
  if(ENABLE_QUAD_PRECISION)
    list(APPEND _fv3atm_defs_private ENABLE_QUAD_PRECISION)
  endif()
  if(32BIT)
    list(APPEND _fv3atm_defs_private OVERLOAD_R4 OVERLOAD_R8)
  endif()

  list(APPEND moving_nest_srcs
    moving_nest/bounding_box.F90
    moving_nest/fv_tracker.F90
    moving_nest/fv_moving_nest.F90
    moving_nest/fv_moving_nest_main.F90
    moving_nest/fv_moving_nest_physics.F90
    moving_nest/fv_moving_nest_types.F90
    moving_nest/fv_moving_nest_utils.F90
  )
else()
   list(APPEND moving_nest_srcs "")
endif()

add_library(fv3atm
  atmos_model.F90
  fv3_cap.F90
  module_fv3_config.F90
  module_fcst_grid_comp.F90
  stochastic_physics/stochastic_physics_wrapper.F90
  cpl/module_block_data.F90
  cpl/module_cplfields.F90
  cpl/module_cap_cpl.F90
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

add_dependencies(fv3atm fv3 fv3ccpp stochastic_physics)

list(APPEND _fv3atm_defs_private GFS_PHYS
                                 INTERNAL_FILE_NML
                                 use_WRTCOMP)

target_compile_definitions(fv3atm PRIVATE "${_fv3atm_defs_private}")

set_target_properties(fv3atm PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
target_include_directories(fv3atm INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>
                                            $<INSTALL_INTERFACE:mod>)

if(MULTI_GASES)
  list(APPEND _fv3atm_defs_private MULTI_GASES)
endif()

target_link_libraries(fv3atm PUBLIC fv3
                                    fv3ccpp
                                    stochastic_physics
                                    fms)

target_link_libraries(fv3atm PUBLIC w3emc::w3emc_d
                                    sp::sp_d
                                    bacio::bacio_4
                                    esmf)

if(INLINE_POST)
  target_link_libraries(fv3atm PUBLIC upp::upp)
endif()

if(OPENMP)
  target_link_libraries(fv3atm PUBLIC OpenMP::OpenMP_Fortran)
endif()

###############################################################################
### Install
###############################################################################
install(
  TARGETS fv3atm
  EXPORT fv3atm-config
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod DESTINATION ${CMAKE_INSTALL_PREFIX})

install(EXPORT fv3atm-config
  DESTINATION lib/cmake)
