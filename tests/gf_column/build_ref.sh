#!/bin/bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Builds MPAS's own Grell-Freitas convection with its driver, the physics interface and the physics tendencies as a
# column program with the reference compiler and flags (the build log of the mpas842 reference: gfortran 14.2 -O3
# -ffree-line-length-none -fconvert=big-endian -ffree-form -fdefault-real-8 -fdefault-double-8, -Dmpas where the log
# has it):
#   build_ref.sh <mpas v8.4.2 source tree> <extracted mmm-physics a4baf7f tree> <build dir>
# MPAS sources compiled unchanged, whole files: framework/mpas_kind_types.F, framework/mpas_constants.F,
# physics/ccpp_kind_types.F, mpas_atmphys_constants.F, mpas_atmphys_vars.F, mpas_atmphys_utilities.F,
# mpas_atmphys_interface.F, mpas_atmphys_driver_convection.F, mpas_atmphys_todynamics.F, physics_wrf/module_cu_gf.mpas.F
# and the other convection modules driver_convection uses (module_cu_ntiedtke.F, cu_ntiedtke_pre.F, cu_ntiedtke_post.F,
# module_cu_kfeta.F, module_cu_tiedtke.F, MMM-physics cu_ntiedtke.F90 without -Dmpas, as the build log). The framework
# modules are the stand-ins of tests/cu_gwdo_column/stubs.F90. The harness program (harness.F90) only moves data.
set -eu
mpas=$(realpath "$1"); mmm=$(realpath "$2"); out=$3; here=$(dirname "$(realpath "$0")")
FFLAGS="-O3 -ffree-line-length-none -fconvert=big-endian -ffree-form -fdefault-real-8 -fdefault-double-8 ${EXTRA_FFLAGS:-}"
CPP="-D_MPI -DCORE_ATMOSPHERE -DMPAS_NAMELIST_SUFFIX=atmosphere -DMPAS_EXE_NAME=atmosphere_model -DMPAS_NATIVE_TIMERS -DMPAS_GIT_VERSION=v8.4.2 -DMPAS_BUILD_TARGET=gfortran -DMPAS_SMIOL_SUPPORT -DMPAS_USE_MPI_F08"
gfortran --version | head -1
mkdir -p "$out"; cd "$out"
phys=$mpas/src/core_atmosphere/physics
c() { gfortran $FFLAGS "$@"; }
c -c $mpas/src/framework/mpas_kind_types.F
c $CPP -c $mpas/src/framework/mpas_constants.F
c -c "$here/../cu_gwdo_column/stubs.F90"
c $CPP -Dmpas -c $phys/ccpp_kind_types.F
c $CPP -c $mmm/cu_ntiedtke.F90
c $CPP -Dmpas -c $phys/mpas_atmphys_constants.F
for f in cu_ntiedtke_pre cu_ntiedtke_post module_cu_ntiedtke module_cu_gf.mpas module_cu_kfeta module_cu_tiedtke; do
  c $CPP -Dmpas -c $phys/physics_wrf/$f.F
done
for f in mpas_atmphys_vars mpas_atmphys_utilities mpas_atmphys_interface mpas_atmphys_driver_convection; do
  c $CPP -Dmpas -c $phys/$f.F
done
c $CPP -Dmpas -I$mpas/src/framework -c $phys/mpas_atmphys_todynamics.F
gfortran -O1 -ffree-line-length-none -fdefault-real-8 -fdefault-double-8 -c "$here/harness.F90"
gfortran -o gf_column_ref harness.o mpas_atmphys_todynamics.o mpas_atmphys_driver_convection.o \
  mpas_atmphys_interface.o mpas_atmphys_utilities.o mpas_atmphys_vars.o module_cu_tiedtke.o module_cu_kfeta.o \
  module_cu_gf.mpas.o module_cu_ntiedtke.o cu_ntiedtke_post.o cu_ntiedtke_pre.o mpas_atmphys_constants.o \
  cu_ntiedtke.o ccpp_kind_types.o stubs.o mpas_constants.o mpas_kind_types.o -lm
