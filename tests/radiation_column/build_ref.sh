#!/bin/bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Builds MPAS's own cloudiness, RRTMG LW / SW radiation drivers, RRTMG table initialisation and ozone climatology as a
# column program with the reference compiler and flags (the build log of the mpas842 reference: gfortran 14.2 -O3
# -ffree-line-length-none -fconvert=big-endian -ffree-form -fdefault-real-8 -fdefault-double-8, and -Dmpas):
#   build_ref.sh <mpas v8.4.2 source tree> <build dir>
# MPAS sources compiled unchanged, whole files: framework/mpas_kind_types.F, mpas_constants.F; physics/
# mpas_atmphys_constants.F, mpas_atmphys_utilities.F, mpas_atmphys_functions.F, mpas_atmphys_date_time.F,
# mpas_atmphys_vars.F, mpas_atmphys_camrad_init.F, mpas_atmphys_rrtmg_lwinit.F, mpas_atmphys_rrtmg_swinit.F,
# mpas_atmphys_o3climatology.F, mpas_atmphys_driver_cloudiness.F, mpas_atmphys_driver_radiation_sw.F,
# mpas_atmphys_driver_radiation_lw.F and the physics_wrf modules they use (module_ra_rrtmg_vinterp.F,
# module_ra_rrtmg_lw.F, module_ra_rrtmg_sw.F, module_ra_rrtmg_sw_aerosols.F, module_ra_cam_support.F,
# module_ra_cam.F, module_mp_thompson_cldfra3.F, module_mp_thompson_aerosols.F). Framework modules (pools, timers,
# log, dm_info, units) and the time variables of mpas_atmphys_manager are the stand-ins of stubs.F90. The harness
# program (harness.F90, dump_tables.F90) only moves data.
set -eu
mpas=$(realpath "$1"); out=$2; here=$(dirname "$(realpath "$0")")
FFLAGS="-O3 -ffree-line-length-none -fconvert=big-endian -ffree-form -fdefault-real-8 -fdefault-double-8"
CPP="-D_MPI -DCORE_ATMOSPHERE -DMPAS_NAMELIST_SUFFIX=atmosphere -DMPAS_EXE_NAME=atmosphere_model -DMPAS_NATIVE_TIMERS -DMPAS_GIT_VERSION=v8.4.2 -DMPAS_BUILD_TARGET=gfortran -DMPAS_SMIOL_SUPPORT -DMPAS_USE_MPI_F08 -Dmpas"
gfortran --version | head -1
mkdir -p "$out"; cd "$out"
phys=$mpas/src/core_atmosphere/physics
c() { gfortran $FFLAGS "$@"; }
c -c $mpas/src/framework/mpas_kind_types.F
c $CPP -c $mpas/src/framework/mpas_constants.F
c -c "$here/stubs.F90"
for f in mpas_atmphys_constants mpas_atmphys_utilities mpas_atmphys_functions mpas_atmphys_date_time mpas_atmphys_vars; do
  c $CPP -c $phys/$f.F
done
for f in module_ra_rrtmg_vinterp module_ra_rrtmg_lw module_ra_rrtmg_sw module_ra_rrtmg_sw_aerosols module_ra_cam_support \
         module_ra_cam module_mp_thompson_cldfra3 module_mp_thompson_aerosols; do
  c $CPP -c $phys/physics_wrf/$f.F
done
for f in mpas_atmphys_camrad_init mpas_atmphys_rrtmg_lwinit mpas_atmphys_rrtmg_swinit mpas_atmphys_o3climatology \
         mpas_atmphys_driver_cloudiness mpas_atmphys_driver_radiation_sw mpas_atmphys_driver_radiation_lw; do
  c $CPP -c $phys/$f.F
done
objs="mpas_atmphys_driver_radiation_lw.o mpas_atmphys_driver_radiation_sw.o mpas_atmphys_driver_cloudiness.o
  mpas_atmphys_o3climatology.o mpas_atmphys_rrtmg_swinit.o mpas_atmphys_rrtmg_lwinit.o mpas_atmphys_camrad_init.o
  module_mp_thompson_aerosols.o module_mp_thompson_cldfra3.o module_ra_cam.o module_ra_cam_support.o
  module_ra_rrtmg_sw_aerosols.o module_ra_rrtmg_sw.o module_ra_rrtmg_lw.o module_ra_rrtmg_vinterp.o mpas_atmphys_vars.o
  mpas_atmphys_date_time.o mpas_atmphys_functions.o mpas_atmphys_utilities.o mpas_atmphys_constants.o stubs.o
  mpas_constants.o mpas_kind_types.o"
gfortran -O1 -ffree-line-length-none -fconvert=big-endian -fdefault-real-8 -fdefault-double-8 -c "$here/harness.F90"
gfortran -o radiation_column_ref harness.o $objs -lm
# the RRTMG tables after rrtmg_initlw_forMPAS / rrtmg_initsw_forMPAS in CPPWRF's dump layout (big-endian unformatted
# reads of RRTMG_*_DATA.DBL as MPAS's build: this program is compiled with -fconvert=big-endian too)
gfortran -O1 -ffree-line-length-none -fconvert=big-endian -fdefault-real-8 -fdefault-double-8 -c "$here/dump_tables.F90"
gfortran -o radiation_dump_tables dump_tables.o $objs -lm
