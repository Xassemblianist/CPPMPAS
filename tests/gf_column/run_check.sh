#!/bin/bash
# SPDX-License-Identifier: GPL-3.0-or-later
# The Grell-Freitas column check: MPAS's own MPAS_to_physics, update_convection_step1, driver_convection
# (cu_grell_freitas), update_convection_step2 and physics_get_tend (build_ref.sh, gfortran with the reference flags) and
# CPPMPAS's (test_gf_column) on the same synthetic columns, compared bit for bit.
#   run_check.sh <mpas src> <mmm-physics a4baf7f tree> <work dir> <test_gf_column> [case ...]
# A case is nCellsSolve:nVertLevels:nsteps:seed:len_disp:acrain:closure_deep:closure_shallow:dt (columns.py). The
# default list has odd and even cell and level counts (1, 2, 3, 7, 14, 64, 63 cells; 55 levels as RG, 56, 41), RG's
# len_disp 240 km and 3 km (the scale-aware branch), the Registry's closures (0, 8) and others (deep 1, 7, 10, 14,
# -1; shallow 0, 1, 4, 7, -1), the rainc bucket, and dt 600 s (RG's), 300, 90 and 45 s.
set -u
mpas=$1; mmm=$2; work=$3; cpp=$4; shift 4
here=$(dirname "$(realpath "$0")")
cases=${*:-"1:55:2:1:240000:0:0:8:600 2:55:2:2:240000:0:0:8:600 3:56:2:3:3000:0:0:8:600 7:55:2:4:240000:1:0:8:600 14:55:2:5:3000:0:0:8:300 64:55:2:6:240000:0:0:8:600 63:41:2:7:240000:1:0:8:90 64:56:2:8:3000:0:0:8:600 21:55:2:9:240000:0:1:0:600 21:55:2:10:240000:0:7:1:600 21:55:2:11:3000:0:10:4:600 21:55:2:12:240000:0:14:7:300 21:55:2:13:240000:0:-1:-1:45"}
mkdir -p "$work"
"$here/build_ref.sh" "$mpas" "$mmm" "$work/ref" > "$work/build_ref.log" 2>&1 || { cat "$work/build_ref.log"; echo "FAIL: reference build"; exit 1; }
head -1 "$work/build_ref.log"
rc=0
for c in $cases; do
  IFS=: read -r nc nv ns seed ld acr cd cs dt <<< "$c"
  in=$work/in_$c.bin
  python3 -I "$here/columns.py" make "$in" "$nc" "$nv" "$ns" "$seed" "$ld" "$acr" "$cd" "$cs" "$dt" || { echo "FAIL: columns $c"; rc=1; continue; }
  "$work/ref/gf_column_ref" "$in" "$work/ref_$c.bin" > "$work/ref_$c.log" 2>&1 || { cat "$work/ref_$c.log"; echo "FAIL: reference run $c"; rc=1; continue; }
  "$cpp" "$in" "$work/cpp_$c.bin" || { echo "FAIL: CPPMPAS run $c"; rc=1; continue; }
  echo -n "case $c: "
  python3 -I "$here/columns.py" stats "$work/ref_$c.bin" "$nc"
  if python3 -I "$here/columns.py" compare "$work/ref_$c.bin" "$work/cpp_$c.bin"; then
    rm -f "$in" "$work/ref_$c.bin" "$work/cpp_$c.bin"  # kept when they differ
  else
    rc=1
  fi
done
[ $rc = 0 ] && echo "gf column check: PASS" || echo "gf column check: FAIL"
exit $rc
