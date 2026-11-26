# SPDX-License-Identifier: GPL-3.0-or-later
cmake_minimum_required(VERSION 3.24)
project(CPPMPAS VERSION 0.0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

# Working precision. The reference is MPAS-Atmosphere built with PRECISION=double (CONSTITUTION.md 1.1), so FP64 is
# the default. The FP32 build (MPAS real*4) has no verified port yet: asking for it is an error, not a silent switch.
option(CPPMPAS_SINGLE_PRECISION "real = float (MPAS PRECISION=single); not ported yet" OFF)
if(CPPMPAS_SINGLE_PRECISION)
  message(FATAL_ERROR "CPPMPAS_SINGLE_PRECISION: the FP32 (MPAS real*4) build has no verified port yet")
endif()
# GPU backend: the kernels are written CPPMPAS_HD-callable (src/core/backend.hpp) so the CUDA build can share them,
# but no CUDA code is ported yet.
option(CPPMPAS_ENABLE_CUDA "Build the CUDA backend; not ported yet" OFF)
if(CPPMPAS_ENABLE_CUDA)
  message(FATAL_ERROR "CPPMPAS_ENABLE_CUDA: there is no CUDA port yet")
endif()
option(CPPMPAS_OPENMP "CPU kernels run with OpenMP (src/core/backend.hpp parallel_for)" OFF)

# Strict IEEE, as CPPWRF (CMakeLists.txt of CPPWRF 7fac5bf): no FMA contraction (the gfortran reference on x86-64
# without -march emits none), no fast-math. Release keeps CMake's -O3 -DNDEBUG, without -march: the reference is
# built without -march. Flags the user gives (CMAKE_CXX_FLAGS, CMAKE_CXX_FLAGS_<CONFIG>) are kept, and checked here.
add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-ffp-contract=off>)
set(_cppmpas_configs ${CMAKE_CONFIGURATION_TYPES})
if(NOT _cppmpas_configs)
  set(_cppmpas_configs "${CMAKE_BUILD_TYPE}")
endif()
set(_cppmpas_x87 "-mfpmath=[^ ]*387|-m32|-mno-sse2?( |$)|-mx32|-fexcess-precision=fast")
foreach(_cfg IN LISTS _cppmpas_configs)
  string(TOUPPER "${_cfg}" _bt)
  if("${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${_bt}}" MATCHES
     "(-Ofast|-f(fast-math|unsafe-math-optimizations|associative-math|reciprocal-math|finite-math-only)|-ffp-contract=(fast|on)|${_cppmpas_x87})")
    message(FATAL_ERROR "CMAKE_CXX_FLAGS: '${CMAKE_MATCH_1}' changes the floating-point results (CPPMPAS needs strict IEEE)")
  endif()
endforeach()

# NetCDF-C from the system (nc-config / pkg-config); no MPI.
find_package(PkgConfig QUIET)
if(PKG_CONFIG_FOUND)
  pkg_check_modules(NETCDF QUIET netcdf)
endif()
find_path(NETCDF_INCDIR netcdf.h HINTS ${NETCDF_INCLUDE_DIRS} /usr/include /usr/local/include)
find_library(NETCDF_LIB NAMES netcdf HINTS ${NETCDF_LIBRARY_DIRS} /usr/lib /usr/lib/x86_64-linux-gnu /usr/local/lib)
if(NOT NETCDF_INCDIR OR NOT NETCDF_LIB)
  message(FATAL_ERROR "NetCDF-C not found (netcdf.h / libnetcdf)")
endif()

add_library(cppmpas STATIC
  src/core/decomp.cpp
  src/core/mesh.cpp
  src/core/namelist.cpp
  src/core_atmosphere/atm_dissipation.cpp
  src/core_atmosphere/atm_dycore.cpp
  src/core_atmosphere/atm_init.cpp
  src/core_atmosphere/atm_solve_diagnostics.cpp
  src/operators/geometry_utils.cpp
  src/operators/vector_operations.cpp
  src/operators/vector_reconstruct.cpp
  src/io/netcdf_file.cpp
  src/io/registry.cpp
  src/io/mpas_reader.cpp
  src/io/mpas_writer.cpp
)
target_include_directories(cppmpas PUBLIC ${CMAKE_SOURCE_DIR}/src ${NETCDF_INCDIR})
target_link_libraries(cppmpas PUBLIC ${NETCDF_LIB})
target_compile_options(cppmpas PRIVATE -Wall -Wextra)
if(CPPMPAS_OPENMP)
  find_package(OpenMP REQUIRED COMPONENTS CXX)
  target_link_libraries(cppmpas PUBLIC OpenMP::OpenMP_CXX)
endif()

add_executable(cppmpas_cli tools/cppmpas.cpp)
set_target_properties(cppmpas_cli PROPERTIES OUTPUT_NAME cppmpas)
target_link_libraries(cppmpas_cli PRIVATE cppmpas)
target_compile_options(cppmpas_cli PRIVATE -Wall -Wextra)

enable_testing()
add_subdirectory(tests)
