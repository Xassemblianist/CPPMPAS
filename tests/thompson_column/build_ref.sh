#!/bin/bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Builds MPAS's own Thompson path as a column program with the reference compiler and flags:
#   build_ref.sh <mpas v8.4.2 source tree> <extracted mmm-physics a4baf7f tree> <build dir>
# MPAS sources compiled unchanged, with the reference build's flags and -D options (-Dmpas where the build log has it):
# framework/mpas_kind_types.F, framework/mpas_constants.F, physics/ccpp_kind_types.F, mpas_atmphys_constants.F,
# mpas_atmphys_functions.F, mpas_atmphys_vars.F, mpas_atmphys_interface.F, physics_wrf/module_mp_kessler.F,
# module_mp_thompson.F, module_mp_wsm6.F, mpas_atmphys_init_microphysics.F, mpas_atmphys_driver_microphysics.F, and
# MMM-physics' module_libmassv.F90, mp_radar.F90, mp_wsm6.F90, mp_wsm6_effectRad.F90 (the driver's other schemes).
# Framework modules (pools, dmpar, I/O units, timers, log, physics_error_fatal) are the stand-ins of stubs.F90; the
# harness program (harness.F90) and dump.c only move data.
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
for f in module_libmassv mp_radar mp_wsm6 mp_wsm6_effectRad; do c $CPP -c $mmm/$f.F90; done
c $CPP -Dmpas -c $phys/mpas_atmphys_functions.F
c $CPP -Dmpas -c $phys/mpas_atmphys_constants.F
c $CPP -Dmpas -c $phys/physics_wrf/module_mp_kessler.F
c $CPP -Dmpas -c $phys/physics_wrf/module_mp_thompson.F
c $CPP -Dmpas -c $phys/physics_wrf/module_mp_wsm6.F
c $CPP -Dmpas -c $phys/mpas_atmphys_vars.F
c $CPP -Dmpas -c $phys/mpas_atmphys_interface.F
c $CPP -Dmpas -c $phys/mpas_atmphys_init_microphysics.F
c $CPP -Dmpas -c $phys/mpas_atmphys_driver_microphysics.F
gfortran -O1 -ffree-line-length-none -fdefault-real-8 -fdefault-double-8 -fconvert=big-endian -c "$here/harness.F90"
gcc -O1 -c "$here/dump.c"
gfortran -o thompson_column_ref harness.o dump.o mpas_atmphys_driver_microphysics.o mpas_atmphys_init_microphysics.o \
  mpas_atmphys_interface.o mpas_atmphys_vars.o module_mp_wsm6.o module_mp_thompson.o module_mp_kessler.o \
  mpas_atmphys_constants.o mpas_atmphys_functions.o mp_wsm6_effectRad.o mp_wsm6.o mp_radar.o module_libmassv.o \
  ccpp_kind_types.o stubs.o mpas_constants.o mpas_kind_types.o -lm
