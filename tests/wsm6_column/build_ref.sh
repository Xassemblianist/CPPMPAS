#!/bin/bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Builds MPAS's own WSM6 path as a column program with the reference compiler and flags:
#   build_ref.sh <mpas v8.4.2 source tree> <extracted mmm-physics a4baf7f tree> <build dir>
# MPAS sources compiled unchanged: framework/mpas_kind_types.F, framework/mpas_constants.F,
# physics/ccpp_kind_types.F, mpas_atmphys_constants.F, mpas_atmphys_vars.F, mpas_atmphys_interface.F,
# physics_wrf/module_mp_wsm6.F and MMM-physics' module_libmassv.F90, mp_radar.F90, mp_wsm6.F90,
# mp_wsm6_effectRad.F90, and physics/mpas_atmphys_functions.F (RSLF / RSIF, which the reference's driver_microphysics
# calls). From mpas_atmphys_driver_microphysics.F the build takes, verbatim, lines 98-257
# (allocate / deallocate_microphysics), 317-500 without the Kessler / Thompson cases (371-436), and 503-850
# (precip_from_MPAS, precip_to_MPAS, compute_radar_reflectivity without its Thompson case 738-791, compute_relhum).
# Framework modules (pools, timers, log, derived types, physics_error_fatal) are the stand-ins of stubs.F90. The
# harness program (harness.F90) only moves data.
set -eu
mpas=$(realpath "$1"); mmm=$(realpath "$2"); out=$3; here=$(dirname "$(realpath "$0")")
FFLAGS="-O3 -ffree-line-length-none -fconvert=big-endian -ffree-form -fdefault-real-8 -fdefault-double-8"
CPP="-D_MPI -DCORE_ATMOSPHERE -DMPAS_NAMELIST_SUFFIX=atmosphere -DMPAS_EXE_NAME=atmosphere_model -DMPAS_NATIVE_TIMERS -DMPAS_GIT_VERSION=v8.4.2 -DMPAS_BUILD_TARGET=gfortran -DMPAS_SMIOL_SUPPORT -DMPAS_USE_MPI_F08"
gfortran --version | head -1
mkdir -p "$out"; cd "$out"
phys=$mpas/src/core_atmosphere/physics
drv=$phys/mpas_atmphys_driver_microphysics.F
{
  echo ' module drv_extract'
  sed -n '10,17p' "$drv" | grep -v 'module_mp_kessler\|module_mp_thompson\|mpas_atmphys_init_microphysics\|use mp_wsm6'
  echo ' use module_mp_wsm6,only: wsm6'
  echo ' use mp_wsm6,only: refl10cm_wsm6'
  echo ' use mpas_atmphys_functions,only: rslf,rsif'
  echo ' use mpas_atmphys_utilities,only: physics_error_fatal'
  echo ' implicit none'
  echo ' contains'
  sed -n '98,257p' "$drv"
  sed -n '317,370p' "$drv"
  sed -n '437,500p' "$drv"
  sed -n '503,737p' "$drv"
  sed -n '792,850p' "$drv"
  echo ' end module drv_extract'
} > drv_extract.F
c() { gfortran $FFLAGS "$@"; }
c -c $mpas/src/framework/mpas_kind_types.F
c $CPP -c $mpas/src/framework/mpas_constants.F
c -c "$here/stubs.F90"
c $CPP -Dmpas -c $phys/ccpp_kind_types.F
for f in module_libmassv mp_radar mp_wsm6 mp_wsm6_effectRad; do c $CPP -c $mmm/$f.F90; done
c $CPP -Dmpas -c $phys/physics_wrf/module_mp_wsm6.F
c $CPP -c $phys/mpas_atmphys_constants.F
c $CPP -c $phys/mpas_atmphys_vars.F
c $CPP -c $phys/mpas_atmphys_interface.F
c $CPP -c $phys/mpas_atmphys_functions.F
c $CPP -c drv_extract.F
gfortran -O1 -ffree-line-length-none -fdefault-real-8 -fdefault-double-8 -c "$here/harness.F90"
gfortran -o wsm6_column_ref harness.o drv_extract.o mpas_atmphys_functions.o mpas_atmphys_interface.o mpas_atmphys_vars.o \
  mpas_atmphys_constants.o module_mp_wsm6.o mp_wsm6_effectRad.o mp_wsm6.o mp_radar.o module_libmassv.o \
  ccpp_kind_types.o stubs.o mpas_constants.o mpas_kind_types.o -lm
