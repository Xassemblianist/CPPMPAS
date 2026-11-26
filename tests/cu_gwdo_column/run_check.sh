#!/bin/bash
# SPDX-License-Identifier: GPL-3.0-or-later
# The GWDO / New Tiedtke column check: MPAS's own MPAS_to_physics, driver_gwdo (bl_ysu_gwdo), update_convection_step1,
# driver_convection (cu_ntiedtke), update_convection_step2 and physics_get_tend (build_ref.sh, gfortran with the
# reference flags) and CPPMPAS's (test_cu_gwdo_column) on the same synthetic columns, compared bit for bit.
#   run_check.sh <mpas src> <mmm-physics a4baf7f tree> <work dir> <test_cu_gwdo_column> [case ...]
# A case is nCellsSolve:nVertLevels:nsteps:seed:len_disp:n_cu:acrain:restart:dt (columns.py). The default list has odd
# and even counts (1, 2, 3, 7, 12, 64, 63 cells: the 1-trip loops, one lane pair, odd tails, long vector bodies), 55
# levels as RG and 41, RG's len_disp 240 km and dt 600 s with dx on both sides of New Tiedtke's 15 km switch
# (len_disp 12 km), a convection step of two dynamics steps (n_cu = 2), the rainc bucket, a restart (initflag 0 at
# the first step) and dt 300 s; three steps each, so initflag 1 and 0 both run.
set -u
mpas=$1; mmm=$2; work=$3; cpp=$4; shift 4
here=$(dirname "$(realpath "$0")")
cases=${*:-"1:55:3:1:240000:1:0:0:600 2:55:3:2:240000:1:0:0:600 3:55:3:3:12000:1:0:0:600 7:55:3:4:240000:1:0:0:600 12:55:3:5:12000:2:1:0:600 64:55:3:6:240000:1:0:0:600 63:41:3:7:12000:1:1:1:300 13:55:2:8:240000:1:0:1:600"}
mkdir -p "$work"
"$here/build_ref.sh" "$mpas" "$mmm" "$work/ref" > "$work/build_ref.log" 2>&1 || { cat "$work/build_ref.log"; echo "FAIL: reference build"; exit 1; }
head -1 "$work/build_ref.log"
rc=0
for c in $cases; do
  IFS=: read -r nc nv ns seed ld ncu acr rst dt <<< "$c"
  in=$work/in_$c.bin
  python3 -I "$here/columns.py" make "$in" "$nc" "$nv" "$ns" "$seed" "$ld" "$ncu" "$acr" "$rst" "$dt" || { echo "FAIL: columns $c"; rc=1; continue; }
  "$work/ref/cu_gwdo_column_ref" "$in" "$work/ref_$c.bin" > "$work/ref_$c.log" 2>&1 || { cat "$work/ref_$c.log"; echo "FAIL: reference run $c"; rc=1; continue; }
  "$cpp" "$in" "$work/cpp_$c.bin" || { echo "FAIL: CPPMPAS run $c"; rc=1; continue; }
  echo -n "case $c: "
  python3 -I "$here/columns.py" stats "$work/ref_$c.bin" "$nc"
  if python3 -I "$here/columns.py" compare "$work/ref_$c.bin" "$work/cpp_$c.bin"; then
    rm -f "$in" "$work/ref_$c.bin" "$work/cpp_$c.bin"  # kept when they differ
  else
    rc=1
  fi
done
[ $rc = 0 ] && echo "cu_gwdo column check: PASS" || echo "cu_gwdo column check: FAIL"
exit $rc
