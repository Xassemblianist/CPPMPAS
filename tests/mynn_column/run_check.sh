#!/bin/bash
# SPDX-License-Identifier: GPL-3.0-or-later
# The MYNN PBL column check: MPAS's own MPAS_to_physics, init_pbl, driver_pbl (bl_mynn) and physics_get_tend
# (build_ref.sh, gfortran with the reference flags) and CPPMPAS's (test_mynn_column) on the same synthetic columns,
# compared bit for bit at every step.
#   run_check.sh <mpas src> <mmm-physics a4baf7f tree> <work dir> <test_mynn_column> [case ...]
# A case is nCellsSolve:nVertLevels:nsteps:seed:restart:dacycling:cu:radt:scaleaware:dheat_opt:len_disp
# (columns.py). The default list: 1, 2, 3, 7, 16 and 64 cells (the 1-trip loops, one pair, odd tails, long loops),
# 55 levels as RG, 41 and 56 levels (odd and even trip counts of the dissipative-heating lanes, k = 1..nVertLevels-1)
# and 100 (mynn::kMaxLevels; a first layer thin enough for DMP_mf's superadiabatic test above k = 1), the first
# step (initflag = 1: mym_initialize) and later ones, a restart (initflag = 0 from the first step), DA cycling, the
# convection branches of physics_get_tend (cu_grell_freitas as convection_permitting, off, cu_kain_fritsch),
# radiation on and off, config_len_disp 240 km (the scale-aware path at 240 km), 12 km, 3 km and 500 m
# (Psig_bl < 1), config_mynn_scaleaware 0 and config_mynn_dheat_opt 0.
set -u
[ $# -ge 4 ] || { echo "FAIL: usage: run_check.sh <mpas src> <mmm-physics> <work dir> <test_mynn_column> (CPPMPAS_MPAS_SRC / CPPMPAS_MMM_PHYSICS unset?)"; exit 1; }
mpas=$1; mmm=$2; work=$3; cpp=$4; shift 4
[ -f "$mpas/src/core_atmosphere/physics/mpas_atmphys_driver_pbl.F" ] || { echo "FAIL: no MPAS v8.4.2 source tree at '$mpas' (CPPMPAS_MPAS_SRC)"; exit 1; }
[ -f "$mmm/bl_mynn.F90" ] || { echo "FAIL: no MMM-physics a4baf7f tree at '$mmm' (CPPMPAS_MMM_PHYSICS)"; exit 1; }
here=$(dirname "$(realpath "$0")")
cases=${*:-"1:55:3:1:0:0:1:1:1:1:240000 2:55:3:2:0:0:1:1:1:1:240000 3:55:3:3:0:0:0:0:1:1:240000 7:55:3:4:0:0:1:1:1:1:240000 64:55:3:5:0:0:1:1:1:1:240000 64:41:3:6:0:0:2:1:1:1:240000 63:56:3:7:0:0:1:0:1:1:240000 9:55:2:8:1:0:1:1:1:1:240000 16:55:2:9:0:1:1:1:1:1:240000 24:55:3:10:0:0:1:1:1:1:3000 8:41:2:11:0:0:0:1:0:0:240000 33:56:2:12:1:0:1:1:1:1:12000 24:55:2:13:0:0:1:1:1:1:500 16:100:2:14:0:0:1:1:1:1:240000"}
mkdir -p "$work"
"$here/build_ref.sh" "$mpas" "$mmm" "$work/ref" > "$work/build_ref.log" 2>&1 || { cat "$work/build_ref.log"; echo "FAIL: reference build"; exit 1; }
head -1 "$work/build_ref.log"
rc=0
for c in $cases; do
  IFS=: read -r nc nv ns seed rs da cu rad sa dh ld <<< "$c"
  in=$work/in_$c.bin
  python3 -I "$here/columns.py" make "$in" "$nc" "$nv" "$ns" "$seed" "$rs" "$da" "$cu" "$rad" "$sa" "$dh" "$ld" || { echo "FAIL: columns $c"; rc=1; continue; }
  "$work/ref/mynn_column_ref" "$in" "$work/ref_$c.bin" > "$work/ref_$c.log" 2>&1 || { cat "$work/ref_$c.log"; echo "FAIL: reference run $c"; rc=1; continue; }
  "$cpp" "$in" "$work/cpp_$c.bin" || { echo "FAIL: CPPMPAS run $c"; rc=1; continue; }
  echo -n "case $c: "
  if python3 -I "$here/columns.py" compare "$work/ref_$c.bin" "$work/cpp_$c.bin"; then
    rm -f "$in" "$work/ref_$c.bin" "$work/cpp_$c.bin"  # kept when they differ
  else
    rc=1
  fi
done
[ $rc = 0 ] && echo "mynn column check: PASS" || echo "mynn column check: FAIL"
exit $rc
