#!/bin/bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Builds MPAS's own surface layer, YSU PBL, physics interface and physics tendency code as a column program with the
# reference compiler and flags (the build log of the mpas842 reference: gfortran 14.2 -O3 -ffree-line-length-none
# -fconvert=big-endian -ffree-form -fdefault-real-8 -fdefault-double-8, -Dmpas where the log has it):
#   build_ref.sh <mpas v8.4.2 source tree> <extracted mmm-physics a4baf7f tree> <build dir>
# MPAS sources compiled unchanged, whole files: framework/mpas_kind_types.F, framework/mpas_constants.F,
# physics/ccpp_kind_types.F, mpas_atmphys_constants.F, mpas_atmphys_vars.F, mpas_atmphys_interface.F,
# mpas_atmphys_driver_sfclayer.F, mpas_atmphys_driver_pbl.F, mpas_atmphys_todynamics.F, physics_wrf/
# module_sf_sfclayrev.F, sf_sfclayrev_pre.F, module_bl_ysu.F and the MYNN / MM5 modules the two drivers use
# (module_sf_sfclay.F, module_sf_mynn.F, sf_mynn_pre.F, module_bl_mynn.F, bl_mynn_pre.F, bl_mynn_post.F); from
# MMM-physics bl_ysu.F90, sf_sfclayrev.F90, mynn_shared.F90, bl_mynn_subroutines.F90, bl_mynn.F90, sf_mynn.F90.
# Framework modules (pools, timers, log, derived types, inner dimensions) are the stand-ins of stubs.F90. The harness
# program (harness.F90, with tables.c) only moves data.
set -eu
mpas=$(realpath "$1"); mmm=$(realpath "$2"); out=$3; here=$(dirname "$(realpath "$0")")
FFLAGS="-O3 -ffree-line-length-none -fconvert=big-endian -ffree-form -fdefault-real-8 -fdefault-double-8"
CPP="-D_MPI -DCORE_ATMOSPHERE -DMPAS_NAMELIST_SUFFIX=atmosphere -DMPAS_EXE_NAME=atmosphere_model -DMPAS_NATIVE_TIMERS -DMPAS_GIT_VERSION=v8.4.2 -DMPAS_BUILD_TARGET=gfortran -DMPAS_SMIOL_SUPPORT -DMPAS_USE_MPI_F08"
gfortran --version | head -1
mkdir -p "$out"; cd "$out"
phys=$mpas/src/core_atmosphere/physics
c() { gfortran $FFLAGS "$@"; }
c -c $mpas/src/framework/mpas_kind_types.F
c $CPP -c $mpas/src/framework/mpas_constants.F
c -c "$here/stubs.F90"
c $CPP -Dmpas -c $phys/ccpp_kind_types.F
for f in mynn_shared bl_mynn_subroutines bl_mynn sf_mynn bl_ysu sf_sfclayrev; do c $CPP -c $mmm/$f.F90; done
for f in bl_mynn_post bl_mynn_pre module_bl_mynn sf_mynn_pre module_sf_mynn module_sf_sfclay sf_sfclayrev_pre \
         module_sf_sfclayrev module_bl_ysu; do
  c $CPP -Dmpas -c $phys/physics_wrf/$f.F
done
for f in mpas_atmphys_constants mpas_atmphys_vars mpas_atmphys_interface mpas_atmphys_driver_sfclayer \
         mpas_atmphys_driver_pbl; do
  c $CPP -Dmpas -c $phys/$f.F
done
c $CPP -Dmpas -I$mpas/src/framework -c $phys/mpas_atmphys_todynamics.F
gfortran -O1 -ffree-line-length-none -fdefault-real-8 -fdefault-double-8 -c "$here/harness.F90"
gcc -O1 -c "$here/tables.c"
gfortran -o pbl_column_ref harness.o tables.o mpas_atmphys_todynamics.o mpas_atmphys_driver_pbl.o mpas_atmphys_driver_sfclayer.o \
  mpas_atmphys_interface.o mpas_atmphys_vars.o mpas_atmphys_constants.o module_bl_ysu.o module_sf_sfclayrev.o \
  sf_sfclayrev_pre.o module_sf_sfclay.o module_sf_mynn.o sf_mynn_pre.o module_bl_mynn.o bl_mynn_pre.o bl_mynn_post.o \
  sf_sfclayrev.o bl_ysu.o sf_mynn.o bl_mynn.o bl_mynn_subroutines.o mynn_shared.o ccpp_kind_types.o stubs.o \
  mpas_constants.o mpas_kind_types.o -lm
