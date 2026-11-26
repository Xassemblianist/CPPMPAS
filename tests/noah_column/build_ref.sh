#!/bin/bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Builds MPAS's own Noah path as a column program with the reference compiler and flags:
#   build_ref.sh <mpas v8.4.2 source tree> <build dir>
# MPAS sources compiled unchanged (flags and -D of the reference build log, which compiles every physics file with
# -Dmpas): framework/mpas_kind_types.F, framework/mpas_constants.F, physics/mpas_atmphys_constants.F,
# mpas_atmphys_utilities.F, mpas_atmphys_date_time.F, mpas_atmphys_vars.F, physics_wrf/module_sf_noahlsm.F,
# module_sf_noahlsm_glacial_only.F, module_sf_noah_seaice.F, module_sf_noah_seaice_drv.F, module_sf_sfcdiags.F,
# module_sf_urban.F, module_bep_bem_helper.F, module_sf_bem.F, module_sf_bep.F, module_sf_bep_bem.F,
# module_sf_noahdrv.F, physics/mpas_atmphys_lsm_noahinit.F, mpas_atmphys_driver_lsm.F, mpas_atmphys_lsm_shared.F,
# mpas_atmphys_driver_seaice.F, mpas_atmphys_sfc_diagnostics.F. The framework modules (pools, timers, log, dmpar,
# io units) and mpas_atmphys_landuse's isurban are the stand-ins of stubs.F90. The harness (harness.F90) only moves
# data.
set -eu
mpas=$(realpath "$1"); out=$2; here=$(dirname "$(realpath "$0")")
FFLAGS="-O3 -ffree-line-length-none -fconvert=big-endian -ffree-form -fdefault-real-8 -fdefault-double-8"
CPP="-D_MPI -DCORE_ATMOSPHERE -DMPAS_NAMELIST_SUFFIX=atmosphere -DMPAS_EXE_NAME=atmosphere_model -DMPAS_NATIVE_TIMERS -DMPAS_GIT_VERSION=v8.4.2 -DMPAS_BUILD_TARGET=gfortran -DMPAS_SMIOL_SUPPORT -DMPAS_USE_MPI_F08"
gfortran --version | head -1
mkdir -p "$out"; cd "$out"
phys=$mpas/src/core_atmosphere/physics
w=$phys/physics_wrf
c() { gfortran $FFLAGS "$@"; }
c -c $mpas/src/framework/mpas_kind_types.F
c $CPP -c $mpas/src/framework/mpas_constants.F
c -c "$here/stubs.F90"
for f in mpas_atmphys_constants mpas_atmphys_utilities mpas_atmphys_date_time mpas_atmphys_vars; do
  c $CPP -Dmpas -c $phys/$f.F
done
for f in module_sf_noahlsm module_sf_noahlsm_glacial_only module_sf_noah_seaice module_sf_noah_seaice_drv \
         module_sf_sfcdiags module_sf_urban module_bep_bem_helper module_sf_bem module_sf_bep module_sf_bep_bem \
         module_sf_noahdrv; do
  c $CPP -Dmpas -c $w/$f.F
done
for f in mpas_atmphys_lsm_noahinit mpas_atmphys_driver_lsm mpas_atmphys_lsm_shared mpas_atmphys_driver_seaice \
         mpas_atmphys_sfc_diagnostics; do
  c $CPP -Dmpas -c $phys/$f.F
done
gfortran -O1 -ffree-line-length-none -fdefault-real-8 -fdefault-double-8 -c "$here/harness.F90"
gfortran -o noah_column_ref harness.o mpas_atmphys_sfc_diagnostics.o mpas_atmphys_driver_seaice.o \
  mpas_atmphys_lsm_shared.o mpas_atmphys_driver_lsm.o mpas_atmphys_lsm_noahinit.o module_sf_noahdrv.o \
  module_sf_bep_bem.o module_sf_bep.o module_sf_bem.o module_bep_bem_helper.o module_sf_urban.o module_sf_sfcdiags.o \
  module_sf_noah_seaice_drv.o module_sf_noah_seaice.o module_sf_noahlsm_glacial_only.o module_sf_noahlsm.o \
  mpas_atmphys_vars.o mpas_atmphys_date_time.o mpas_atmphys_utilities.o mpas_atmphys_constants.o stubs.o \
  mpas_constants.o mpas_kind_types.o -lm
