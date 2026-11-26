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
# GPU backend (src/core/backend.hpp): the same sources built a second time with nvcc as the CUDA library
# cppmpas_cuda and the executable cppmpas_gpu, next to the CPU library cppmpas and executable cppmpas (unchanged).
# CPPMPAS_CUDA is a synonym of CPPMPAS_ENABLE_CUDA.
option(CPPMPAS_ENABLE_CUDA "Also build the CUDA backend (cppmpas_gpu)" OFF)
if(CPPMPAS_CUDA)
  set(CPPMPAS_ENABLE_CUDA ON)
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

if(CPPMPAS_ENABLE_CUDA)
  if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
    set(CMAKE_CUDA_ARCHITECTURES native)  # the GPU of the build machine; pass e.g. 120 for an RTX 50xx
  endif()
  enable_language(CUDA)
  set(CMAKE_CUDA_STANDARD 20)
  set(CMAKE_CUDA_STANDARD_REQUIRED ON)
  find_package(CUDAToolkit REQUIRED)
  # kernels capture views of the whole mesh / state / dycore fields (core/views.hpp): kernel parameters above 4 KB
  # need CUDA 12.1
  if(CMAKE_CUDA_COMPILER_VERSION VERSION_LESS 12.1)
    message(FATAL_ERROR "the CUDA build needs CUDA 12.1 or later (kernel parameters above 4 KB); found ${CMAKE_CUDA_COMPILER_VERSION}")
  endif()
  # CUDA FP64 as the reference's arithmetic, bit for bit (CPPWRF 7fac5bf CMakeLists.txt:53-122): no fused
  # multiply-add contraction (--fmad=false; gfortran without -march emits none), IEEE division and square root,
  # no flush to zero, the host code of nvcc's sources without contraction either; glibc's libm / libmvec on the
  # device (src/core/glibc_math.hpp through src/core/libm.hpp), never CUDA's math library or intrinsics. These are
  # compile options of every CUDA source: they follow CMAKE_CUDA_FLAGS and CMAKE_CUDA_FLAGS_<CONFIG> on nvcc's
  # command line, so a flag the user gives cannot override them. Release keeps CMake's -O3 -DNDEBUG and the user's
  # CMAKE_CUDA_FLAGS_<CONFIG>, which are checked below.
  add_compile_options("$<$<COMPILE_LANGUAGE:CUDA>:--fmad=false;--ftz=false;--prec-div=true;--prec-sqrt=true;-Xcompiler=-ffp-contract=off>")
  add_compile_options("$<$<COMPILE_LANGUAGE:CUDA>:--extended-lambda;--expt-relaxed-constexpr>")
  # flags that change results, in any form nvcc takes them (--use_fast_math / -use_fast_math, --fmad / -fmad,
  # --ftz, --prec-div, --prec-sqrt, with '=' or a blank) and every host-compiler flag the CXX check refuses, however
  # it reaches the host compiler (-Xcompiler, --compiler-options, with '=', a blank or a comma list): nvcc has no
  # option of its own with these names, so any occurrence is refused
  foreach(_cfg IN LISTS _cppmpas_configs)
    string(TOUPPER "${_cfg}" _bt)
    if("${CMAKE_CUDA_FLAGS} ${CMAKE_CUDA_FLAGS_${_bt}}" MATCHES
       "(use_fast_math|fmad[= ]*(true|1)|ftz[= ]*(true|1)|prec-div[= ]*(false|0)|prec-sqrt[= ]*(false|0)|-Ofast|-f(fast-math|unsafe-math-optimizations|associative-math|reciprocal-math|finite-math-only)|-ffp-contract=(fast|on)|${_cppmpas_x87})")
      message(FATAL_ERROR "CMAKE_CUDA_FLAGS: '${CMAKE_MATCH_1}' changes the floating-point results (CPPMPAS needs strict IEEE)")
    endif()
  endforeach()
endif()

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

set(CPPMPAS_SOURCES
  src/core/backend.cpp
  src/core/decomp.cpp
  src/core/mesh.cpp
  src/core/namelist.cpp
  src/core_atmosphere/atm_advance_scalars.cpp
  src/core_atmosphere/atm_dissipation.cpp
  src/core_atmosphere/atm_dissipation_perturbation.cpp
  src/core_atmosphere/atm_dycore.cpp
  src/core_atmosphere/atm_microphysics.cpp
  src/core_atmosphere/atm_init.cpp
  src/core_atmosphere/atm_solve_diagnostics.cpp
  src/operators/geometry_utils.cpp
  src/operators/vector_operations.cpp
  src/operators/vector_reconstruct.cpp
  src/io/netcdf_file.cpp
  src/io/registry.cpp
  src/io/mpas_reader.cpp
  src/io/mpas_restart.cpp
  src/io/mpas_writer.cpp
)
add_library(cppmpas STATIC ${CPPMPAS_SOURCES})
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

# The CUDA build of the same sources: each source is compiled by nvcc through a one-line .cu that includes it
# (a source file has one language per directory in CMake, and the CPU library keeps g++). Separable compilation:
# kernels call the HD functions of other sources (geometry, libm wrappers) as in the CPU build.
if(CPPMPAS_ENABLE_CUDA)
  set(_cuda_srcs src/core/backend_cuda.cu)
  foreach(_s IN LISTS CPPMPAS_SOURCES ITEMS tools/cppmpas.cpp)
    set(_w ${CMAKE_BINARY_DIR}/cuda_src/${_s}.cu)
    file(CONFIGURE OUTPUT ${_w} CONTENT "// generated: the CUDA build of ${_s}\n#include \"${CMAKE_SOURCE_DIR}/${_s}\"\n")
    if(_s STREQUAL "tools/cppmpas.cpp")
      set(_cli_cu ${_w})
    else()
      list(APPEND _cuda_srcs ${_w})
    endif()
  endforeach()
  add_library(cppmpas_cuda STATIC ${_cuda_srcs})
  target_include_directories(cppmpas_cuda PUBLIC ${CMAKE_SOURCE_DIR}/src ${NETCDF_INCDIR})
  target_compile_definitions(cppmpas_cuda PUBLIC CPPMPAS_USE_CUDA)
  target_link_libraries(cppmpas_cuda PUBLIC ${NETCDF_LIB} CUDA::cudart)
  set_target_properties(cppmpas_cuda PROPERTIES CUDA_SEPARABLE_COMPILATION ON POSITION_INDEPENDENT_CODE ON)
  add_executable(cppmpas_gpu ${_cli_cu})
  target_link_libraries(cppmpas_gpu PRIVATE cppmpas_cuda)
  set_target_properties(cppmpas_gpu PROPERTIES CUDA_SEPARABLE_COMPILATION ON CUDA_RESOLVE_DEVICE_SYMBOLS ON)
endif()

enable_testing()
add_subdirectory(tests)
