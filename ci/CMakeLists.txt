# This file is used by fv3atm's root CMakeLists.txt when BUILD_TESTING=ON. It is
# used for CI runs and is based on the ufs-weather-model root CMakeLists.txt. It
# cannot be built directly, and should not be used for compiling for general R&D
# or operations.
#
# Alex Richert, 6 Dec 2023

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMakeModules/Modules)

if(${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
  if(CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 9.0.0)
    message(FATAL_ERROR "GNU Compiler >= 9 is required")
  endif()
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ggdb -fbacktrace -cpp -fcray-pointer -ffree-line-length-none -fno-range-check")

  if(${CMAKE_Fortran_COMPILER_VERSION} VERSION_GREATER_EQUAL 10)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch -fallow-invalid-boz")
  endif()

  if(NOT 32BIT)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8 -fdefault-double-8")
  endif()
elseif(${CMAKE_Fortran_COMPILER_ID} MATCHES "Intel")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -traceback -fpp -fno-alias -auto -safe-cray-ptr -ftz -assume byterecl -nowarn -sox -align array64byte -qno-opt-dynamic-align")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -qno-opt-dynamic-align -sox -fp-model source")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -debug minimal -qoverride-limits")
  set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -fp-model consistent")
  set(CMAKE_C_FLAGS_RELEASE       "-O2 -debug minimal")

  if(NOT 32BIT)
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -real-size 64")
  endif()
endif()

set(32BIT           OFF CACHE BOOL "Enable 32BIT (single precision arithmetic in dycore and fast physics)")
set(CCPP_32BIT      OFF CACHE BOOL "Enable CCPP_32BIT (single precision arithmetic in slow physics)")
set(INLINE_POST     ON  CACHE BOOL "Enable inline post")
set(MULTI_GASES     OFF CACHE BOOL "Enable MULTI_GASES")
set(MOVING_NEST     OFF CACHE BOOL "Enable moving nest code")
set(OPENMP          ON  CACHE BOOL "Enable OpenMP threading")
set(PARALLEL_NETCDF OFF CACHE BOOL "Enable parallel NetCDF")

message("32BIT ............ ${32BIT}")
message("CCPP_32BIT ....... ${CCPP_32BIT}")
message("INLINE_POST ...... ${INLINE_POST}")
message("MULTI_GASES ...... ${MULTI_GASES}")
message("MOVING_NEST ...... ${MOVING_NEST}")
message("OPENMP ........... ${OPENMP}")
message("PARALLEL_NETCDF .. ${PARALLEL_NETCDF}")

find_package(MPI REQUIRED)
if(OPENMP)
  find_package(OpenMP REQUIRED)
endif()

find_package(NetCDF 4.7.4 REQUIRED C Fortran)
find_package(ESMF 8.3.0 MODULE REQUIRED)
find_package(FMS 2022.04 REQUIRED COMPONENTS R4 R8)
if(32BIT)
  add_library(fms ALIAS FMS::fms_r4)
else()
  add_library(fms ALIAS FMS::fms_r8)
endif()
find_package(bacio 2.4.0 REQUIRED)
find_package(sp 2.3.3 REQUIRED)
find_package(w3emc 2.9.2 REQUIRED)

find_package(Python 3.6 REQUIRED COMPONENTS Interpreter)

add_subdirectory(stochastic_physics_repo)
