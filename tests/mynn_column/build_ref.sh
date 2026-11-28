#!/bin/bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Builds MPAS's own MYNN PBL, its MPAS driver, the physics interface and the physics tendencies as a column program
# with the reference compiler and flags (the build log of the mpas842 reference: gfortran 14.2 -O3
# -ffree-line-length-none -fconvert=big-endian -ffree-form -fdefault-real-8 -fdefault-double-8; -Dmpas for the
# physics_wrf and physics files as the log, not for the MMM-physics files):
#   build_ref.sh <mpas v8.4.2 source tree> <extracted mmm-physics a4baf7f tree> <build dir>
# MPAS sources compiled unchanged, whole files: framework/mpas_kind_types.F, framework/mpas_constants.F,
# physics/ccpp_kind_types.F, mpas_atmphys_constants.F, mpas_atmphys_vars.F, mpas_atmphys_interface.F,
# mpas_atmphys_driver_pbl.F, mpas_atmphys_todynamics.F, physics_wrf/module_bl_mynn.F, bl_mynn_pre.F, bl_mynn_post.F
# and module_bl_ysu.F (driver_pbl's module uses it); from MMM-physics mynn_shared.F90, bl_mynn_subroutines.F90,
# bl_mynn.F90 and bl_ysu.F90. Framework modules (pools, timers, log, derived types, inner dimensions) are the
# stand-ins of stubs.F90. The harness program (harness.F90) only moves data.
set -euo pipefail
mpas=$(realpath "$1"); mmm=$(realpath "$2"); out=$3; here=$(dirname "$(realpath "$0")")
FFLAGS="-O3 -ffree-line-length-none -fconvert=big-endian -ffree-form -fdefault-real-8 -fdefault-double-8"
CPP="-D_MPI -DCORE_ATMOSPHERE -DMPAS_NAMELIST_SUFFIX=atmosphere -DMPAS_EXE_NAME=atmosphere_model -DMPAS_NATIVE_TIMERS -DMPAS_GIT_VERSION=v8.4.2 -DMPAS_BUILD_TARGET=gfortran -DMPAS_SMIOL_SUPPORT -DMPAS_USE_MPI_F08"
gfortran --version | sed -n 1p  # sed reads all of it: no SIGPIPE under pipefail
# the reference compiler is gfortran 14.2 (the build log of the mpas842 reference); another version is refused
case "$(gfortran -dumpfullversion)" in 14.2.*) ;; *) echo "FAIL: gfortran $(gfortran -dumpfullversion), the reference is 14.2"; exit 1;; esac
mkdir -p "$out"; cd "$out"
phys=$mpas/src/core_atmosphere/physics
c() { gfortran $FFLAGS "$@"; }
c -c $mpas/src/framework/mpas_kind_types.F
c $CPP -c $mpas/src/framework/mpas_constants.F
c -c "$here/stubs.F90"
c $CPP -Dmpas -c $phys/ccpp_kind_types.F
for f in mynn_shared bl_mynn_subroutines bl_mynn bl_ysu; do c $CPP -c $mmm/$f.F90; done
for f in bl_mynn_post bl_mynn_pre module_bl_mynn module_bl_ysu; do c $CPP -Dmpas -c $phys/physics_wrf/$f.F; done
for f in mpas_atmphys_constants mpas_atmphys_vars mpas_atmphys_interface mpas_atmphys_driver_pbl; do
  c $CPP -Dmpas -c $phys/$f.F
done
c $CPP -Dmpas -I$mpas/src/framework -c $phys/mpas_atmphys_todynamics.F
gfortran -O1 -ffree-line-length-none -fdefault-real-8 -fdefault-double-8 -c "$here/harness.F90"
gfortran -o mynn_column_ref harness.o mpas_atmphys_todynamics.o mpas_atmphys_driver_pbl.o mpas_atmphys_interface.o \
  mpas_atmphys_vars.o mpas_atmphys_constants.o module_bl_ysu.o module_bl_mynn.o bl_mynn_pre.o bl_mynn_post.o \
  bl_ysu.o bl_mynn.o bl_mynn_subroutines.o mynn_shared.o ccpp_kind_types.o stubs.o mpas_constants.o \
  mpas_kind_types.o -lm
