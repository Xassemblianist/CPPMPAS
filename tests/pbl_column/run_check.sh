#!/bin/bash
# SPDX-License-Identifier: GPL-3.0-or-later
# The surface layer / YSU PBL column check: MPAS's own MPAS_to_physics, driver_sfclayer (sf_monin_obukhov_rev),
# driver_pbl (bl_ysu) and physics_get_tend (build_ref.sh, gfortran with the reference flags) and CPPMPAS's
# (test_pbl_column) on the same synthetic columns, compared bit for bit.
#   run_check.sh <mpas src> <mmm-physics a4baf7f tree> <work dir> <test_pbl_column> [case ...]
# A case is nCellsSolve:nVertLevels:nsteps:seed:frac_seaice:ysu_pblmix:cu:radt (columns.py). The default list has
# odd and even counts (1, 2, 3, 7, 14, 64 cells: the 1-trip loops, one lane pair, odd tails, long vector bodies),
# 55 levels as RG and 41, config_frac_seaice true (RG) and false, config_ysu_pblmix false (RG) and true, the
# convection branches of physics_get_tend (cu_ntiedtke as mesoscale_reference, off, cu_kain_fritsch,
# cu_grell_freitas) and radiation on and off.
set -u
mpas=$1; mmm=$2; work=$3; cpp=$4; shift 4
here=$(dirname "$(realpath "$0")")
cases=${*:-"1:55:2:1:1:0:1:1 2:55:2:2:1:0:1:1 3:55:3:3:1:1:1:1 7:55:3:4:1:0:1:1 14:55:3:5:1:1:0:0 64:55:3:6:1:0:1:1 64:41:2:7:0:1:2:1 9:55:2:8:0:0:3:0 63:55:2:9:1:1:1:1"}
mkdir -p "$work"
"$here/build_ref.sh" "$mpas" "$mmm" "$work/ref" > "$work/build_ref.log" 2>&1 || { cat "$work/build_ref.log"; echo "FAIL: reference build"; exit 1; }
head -1 "$work/build_ref.log"
rc=0
for c in $cases; do
  IFS=: read -r nc nv ns seed si mix cu rad <<< "$c"
  in=$work/in_$c.bin
  python3 -I "$here/columns.py" make "$in" "$nc" "$nv" "$ns" "$seed" "$si" "$mix" "$cu" "$rad" || { echo "FAIL: columns $c"; rc=1; continue; }
  "$work/ref/pbl_column_ref" "$in" "$work/ref_$c.bin" > "$work/ref_$c.log" 2>&1 || { cat "$work/ref_$c.log"; echo "FAIL: reference run $c"; rc=1; continue; }
  "$cpp" "$in" "$work/cpp_$c.bin" || { echo "FAIL: CPPMPAS run $c"; rc=1; continue; }
  echo -n "case $c: "
  if python3 -I "$here/columns.py" compare "$work/ref_$c.bin" "$work/cpp_$c.bin"; then
    rm -f "$in" "$work/ref_$c.bin" "$work/cpp_$c.bin"  # kept when they differ
  else
    rc=1
  fi
done
[ $rc = 0 ] && echo "pbl column check: PASS" || echo "pbl column check: FAIL"
exit $rc
