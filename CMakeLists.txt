# SPDX-License-Identifier: GPL-3.0-or-later
cmake_minimum_required(VERSION 3.24)

# Flags that change floating-point results, refused wherever they appear. Each is matched as a whole flag: the
# character before it and the one after it must separate flags (blank, tab, ',', ';', a generator expression's ':' or
# '>', a quote; '=' before it for -Xcompiler=...), so -fno-fast-math, -fsigned-zeros, -mno-fma or --fmad=false never
# match. The host compiler's (also when nvcc passes them on through -Xcompiler / --compiler-options, '=' or a blank,
# comma lists): -Ofast, fast-math and the flags it implies that change results, -fno-signed-zeros, complex arithmetic
# without the range checks, FMA contraction (any -ffp-contract but off), FMA instructions (-mfma, any -march), x87
# arithmetic (-mfpmath with 387, -m32, -mx32, -mno-sse, -mno-sse2, -fexcess-precision=fast). nvcc's (and ptxas's
# through -Xptxas): --use_fast_math, --fmad=true, --ftz=true, --prec-div=false, --prec-sqrt=false, with '=' or a
# blank, one dash or two.
set(_cppmpas_sep "[ \t,;:=>\"']")
set(_cppmpas_end "[ \t,;>\"']")
set(_cppmpas_val "[^ \t,;>\"']")
set(_cppmpas_fp_bad "${_cppmpas_sep}(-Ofast|-ffast-math|-funsafe-math-optimizations|-fassociative-math|-freciprocal-math|-ffinite-math-only|-fno-signed-zeros|-fcx-limited-range|-ffp-contract=(fast|on)${_cppmpas_val}*|-mfma4?|-march=${_cppmpas_val}+|-mfpmath=${_cppmpas_val}*387${_cppmpas_val}*|-m32|-mx32|-mno-sse2?|-fexcess-precision=fast|-?-use_fast_math|-?-fmad[= ]+(true|1)|-?-ftz[= ]+(true|1)|-?-prec-div[= ]+(false|0)|-?-prec-sqrt[= ]+(false|0))${_cppmpas_end}")
# cppmpas_check_fp_flags(<where> <flags>): FATAL_ERROR when <flags> (a command-line string or a CMake list, generator
# expressions included) hold a refused flag
function(cppmpas_check_fp_flags where flags)
  if(" ${flags} " MATCHES "${_cppmpas_fp_bad}")
    message(FATAL_ERROR "${where}: '${CMAKE_MATCH_1}' changes the floating-point results (CPPMPAS needs strict IEEE)")
  endif()
endfunction()
# every CMAKE_<LANG>_FLAGS, CMAKE_<LANG>_FLAGS_<CONFIG> (any configuration, also one not in
# CMAKE_CONFIGURATION_TYPES) and CMAKE_<KIND>_LINKER_FLAGS[_<CONFIG>] (-ffast-math at the link adds crtfastmath.o,
# which sets flush-to-zero), whether CUDA is enabled or not, the environment flags CMake starts them from (CXXFLAGS,
# CFLAGS, CUDAFLAGS, FFLAGS, LDFLAGS) and nvcc's own (NVCC_PREPEND_FLAGS, NVCC_APPEND_FLAGS)
function(cppmpas_check_fp_variables)
  get_cmake_property(_vars VARIABLES)
  foreach(_v IN LISTS _vars)
    if(_v MATCHES "^CMAKE_[A-Za-z]+_FLAGS(_[A-Za-z0-9_]+)?$" OR _v MATCHES "^CMAKE_[A-Z]+_LINKER_FLAGS(_[A-Za-z0-9_]+)?$")
      cppmpas_check_fp_flags("${_v}" "${${_v}}")
    endif()
  endforeach()
  foreach(_e IN ITEMS CXXFLAGS CFLAGS CUDAFLAGS FFLAGS LDFLAGS NVCC_PREPEND_FLAGS NVCC_APPEND_FLAGS)
    cppmpas_check_fp_flags("environment ${_e}" "$ENV{${_e}}")
  endforeach()
endfunction()
# the compile and link options of every directory, every target (and the imported ones the top-level directory sees)
# and every source file of a target, from the top-level directory down
function(cppmpas_check_fp_directory dir)
  foreach(_p IN ITEMS COMPILE_OPTIONS LINK_OPTIONS)
    get_property(_x DIRECTORY "${dir}" PROPERTY ${_p})
    cppmpas_check_fp_flags("directory ${dir} ${_p}" "${_x}")
  endforeach()
  get_property(_targets DIRECTORY "${dir}" PROPERTY BUILDSYSTEM_TARGETS)
  get_property(_imported DIRECTORY "${dir}" PROPERTY IMPORTED_TARGETS)
  foreach(_t IN LISTS _targets _imported)
    if(NOT TARGET ${_t})  # an imported target that is not global is seen only in and below its own directory
      continue()
    endif()
    foreach(_p IN ITEMS COMPILE_OPTIONS INTERFACE_COMPILE_OPTIONS COMPILE_FLAGS LINK_OPTIONS INTERFACE_LINK_OPTIONS
                        LINK_FLAGS STATIC_LIBRARY_OPTIONS)
      get_property(_x TARGET ${_t} PROPERTY ${_p})
      cppmpas_check_fp_flags("target ${_t} ${_p}" "${_x}")
    endforeach()
    get_property(_is_imported TARGET ${_t} PROPERTY IMPORTED)
    if(_is_imported)
      continue()
    endif()
    get_property(_srcs TARGET ${_t} PROPERTY SOURCES)
    foreach(_s IN LISTS _srcs)
      foreach(_p IN ITEMS COMPILE_OPTIONS COMPILE_FLAGS)
        get_property(_x SOURCE "${_s}" TARGET_DIRECTORY ${_t} PROPERTY ${_p})
        cppmpas_check_fp_flags("source ${_s} of target ${_t} ${_p}" "${_x}")
      endforeach()
    endforeach()
  endforeach()
  get_property(_subdirs DIRECTORY "${dir}" PROPERTY SUBDIRECTORIES)
  foreach(_d IN LISTS _subdirs)
    cppmpas_check_fp_directory("${_d}")
  endforeach()
endfunction()
# before project(): a refused flag fails here, not in the compiler checks project() runs with it (-m32 without a
# 32-bit library); again after the languages are enabled and at the end of this file with every target
cppmpas_check_fp_variables()
project(CPPMPAS VERSION 0.0.1 LANGUAGES CXX)
cppmpas_check_fp_variables()

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
# built without -march. Flags the user gives are kept, and checked (cppmpas_check_fp_flags above).
add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-ffp-contract=off>)

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
  # the user's CMAKE_CUDA_FLAGS[_<CONFIG>] and nvcc's defaults are checked with every other flag
  # (cppmpas_check_fp_variables above and at the end of this file)
  cppmpas_check_fp_variables()
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
  src/physics/module_cu_ntiedtke.cpp
  src/physics/module_cu_gf.cpp
  src/physics/module_mp_wsm6.cpp
  src/physics/mpas_atmphys_microphysics.cpp
  src/physics/mpas_atmphys_microphysics_thompson.cpp
  src/physics/mp_thompson_init.cpp
  src/physics/mpas_atmphys_convection.cpp
  src/physics/mpas_atmphys_gwdo.cpp
  src/physics/mpas_atmphys_interface.cpp
  src/physics/mpas_atmphys_pbl.cpp
  src/physics/mpas_atmphys_sfclayer.cpp
  src/physics/mpas_atmphys_todynamics.cpp
  src/physics/module_ra_rrtmg.cpp
  src/physics/mpas_atmphys_radiation.cpp
  src/physics/rrtmg/rrtmg_lw_tables.cpp
  src/physics/rrtmg/rrtmg_sw_tables.cpp
  src/physics/module_sf_noahdrv.cpp
  src/physics/mpas_atmphys_driver_lsm.cpp
  src/physics/mpas_atmphys_lsm_noahinit.cpp
  src/physics/mpas_atmphys_control.cpp
  src/physics/mpas_atmphys_landuse.cpp
  src/core_atmosphere/diagnostics/mpas_atm_diagnostics_manager.cpp
  src/core_atmosphere/diagnostics/mpas_isobaric_diagnostics.cpp
  src/core_atmosphere/diagnostics/mpas_convective_diagnostics.cpp
  src/physics/mpas_atmphys_driver.cpp
  src/operators/vector_operations.cpp
  src/operators/vector_reconstruct.cpp
  src/io/netcdf_file.cpp
  src/io/registry.cpp
  src/io/registry_pools.cpp
  src/io/mpas_reader.cpp
  src/io/mpas_restart.cpp
  src/io/mpas_writer.cpp
  src/io/mpas_streams.cpp
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

# every flag of the finished configuration: the variables (also those enable_language set) and every directory,
# target and source file
cppmpas_check_fp_variables()
cppmpas_check_fp_directory("${CMAKE_SOURCE_DIR}")
