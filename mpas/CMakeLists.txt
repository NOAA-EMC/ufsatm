cmake_minimum_required(VERSION 3.19)

project(MPAS
  VERSION 1.0.0
  LANGUAGES Fortran)

include(${CMAKE_CURRENT_SOURCE_DIR}/MPAS-Model/cmake/Functions/MPAS_Functions.cmake)

list(INSERT CMAKE_MODULE_PATH 0 ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules)
set(CMAKE_DIRECTORY_LABELS ${PROJECT_NAME})
include(GNUInstallDirs)

# Build product output locations
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Set default build type to RelWithDebInfo
if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "Setting default build type to Release.  Specify CMAKE_BUILD_TYPE to override.")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "CMake Build type" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# Find C pre-processor
if(CMAKE_C_COMPILER_ID MATCHES GNU)
  find_program(CPP_EXECUTABLE NAMES cpp REQUIRED)
  set(CPP_EXTRA_FLAGS -traditional)
elseif(CMAKE_C_COMPILER_ID MATCHES "(Apple)?Clang" )
  find_program(CPP_EXECUTABLE NAMES cpp REQUIRED)
else()
  message(STATUS "Unknown compiler: ${CMAKE_C_COMPILER_ID}")
  set(CPP_EXECUTABLE ${CMAKE_C_COMPILER})
endif()

# Fortran module output directory for build interface
set(MPAS_MODULE_DIR ${PROJECT_NAME}/module/${CMAKE_Fortran_COMPILER_ID}/${CMAKE_Fortran_COMPILER_VERSION})

# Install Fortran module directory
install(DIRECTORY ${CMAKE_BINARY_DIR}/${MPAS_MODULE_DIR}/ DESTINATION ${CMAKE_INSTALL_LIBDIR}/${MPAS_MODULE_DIR}/)

###############################################################################
# Build MPAS libraries...
###############################################################################
add_definitions(-DMPAS_CAM_DYCORE)
add_definitions(-DMPAS_UFS_DYCORE)

# MPAS Utilities (Externals)
add_subdirectory(MPAS-Model/src/external/ezxml)

# ESMF libraries.
if(NOT ESMF_FOUND)
  find_package(ESMF REQUIRED)
endif()
add_definitions(-DMPAS_EXTERNAL_ESMF_LIB -DMPAS_NO_ESMF_INIT)
add_library(${PROJECT_NAME}::external::esmf ALIAS ESMF::ESMF)

# MPAS Namelist
add_subdirectory(MPAS-Model/src/tools/input_gen) # Targets: namelist_gen, streams_gen

# MPAS Registry
add_subdirectory(MPAS-Model/src/tools/registry)  # Targets: mpas_parse_<core_name>

# MPAS framework
add_subdirectory(MPAS-Model/src/framework)       # Target: MPAS::framework

# MPAS operators
add_subdirectory(MPAS-Model/src/operators)       # Target: MPAS::operators

# MPAS atmosphere
add_subdirectory(MPAS-Model/src/core_atmosphere) # Target: core_atmosphere
add_library(mpas ALIAS core_atmosphere)

############################################################################### 
# Package Configurations
############################################################################### 
include(CMakePackageConfigHelpers)

# Build-tree target exports
export(EXPORT ${PROJECT_NAME}ExportsExternal NAMESPACE ${PROJECT_NAME}::external:: FILE ${PROJECT_NAME}-targets-external.cmake)
export(EXPORT ${PROJECT_NAME}Exports NAMESPACE ${PROJECT_NAME}:: FILE ${PROJECT_NAME}-targets.cmake)
export(EXPORT ${PROJECT_NAME}ExportsCore NAMESPACE ${PROJECT_NAME}::core:: FILE ${PROJECT_NAME}-targets-core.cmake)

# CMake Config file install location
set(CONFIG_INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules DESTINATION ${CONFIG_INSTALL_DESTINATION})

####
set(BINDIR ${CMAKE_BINARY_DIR}/bin)
set(CORE_DATADIR_ROOT ${CMAKE_BINARY_DIR}/${PROJECT_NAME})
set(CMAKE_MODULE_INSTALL_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules)
string(TOLOWER ${PROJECT_NAME} PROJECT_NAME_LOWER)
configure_package_config_file(
  MPAS-Model/cmake/PackageConfig.cmake.in ${PROJECT_NAME_LOWER}-config.cmake
  INSTALL_DESTINATION .
  INSTALL_PREFIX ${CMAKE_CURRENT_BINARY_DIR}
  PATH_VARS BINDIR CORE_DATADIR_ROOT CMAKE_MODULE_INSTALL_PATH)

###
set(BINDIR ${CMAKE_INSTALL_BINDIR})
set(CORE_DATADIR_ROOT ${CMAKE_INSTALL_DATADIR}/${PROJECT_NAME})
set(CMAKE_MODULE_INSTALL_PATH ${CONFIG_INSTALL_DESTINATION}/Modules)
configure_package_config_file(
  MPAS-Model/cmake/PackageConfig.cmake.in install/${PROJECT_NAME_LOWER}-config.cmake
  INSTALL_DESTINATION ${CONFIG_INSTALL_DESTINATION}
  PATH_VARS BINDIR CORE_DATADIR_ROOT CMAKE_MODULE_INSTALL_PATH)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/install/${PROJECT_NAME_LOWER}-config.cmake
  DESTINATION ${CONFIG_INSTALL_DESTINATION})

###
write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME_LOWER}-config-version.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY AnyNewerVersion)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME_LOWER}-config-version.cmake
  DESTINATION ${CONFIG_INSTALL_DESTINATION})

###
install(EXPORT ${PROJECT_NAME}ExportsExternal
        NAMESPACE ${PROJECT_NAME}::external::
        FILE ${PROJECT_NAME_LOWER}-targets-external.cmake
        DESTINATION ${CONFIG_INSTALL_DESTINATION})
install(EXPORT ${PROJECT_NAME}Exports
        NAMESPACE ${PROJECT_NAME}::
        FILE ${PROJECT_NAME_LOWER}-targets.cmake
        DESTINATION ${CONFIG_INSTALL_DESTINATION})
install(EXPORT ${PROJECT_NAME}ExportsCore
	NAMESPACE ${PROJECT_NAME}::core::
        FILE ${PROJECT_NAME_LOWER}-targets-core.cmake
        DESTINATION ${CONFIG_INSTALL_DESTINATION})
