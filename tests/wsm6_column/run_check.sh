#!/bin/bash
# SPDX-License-Identifier: GPL-3.0-or-later
# The WSM6 column check: MPAS's own WSM6 path (build_ref.sh, gfortran with the reference flags) and CPPMPAS's
# (test_wsm6_column) on the same synthetic columns, compared bit for bit.
#   run_check.sh <mpas src> <mmm-physics a4baf7f tree> <work dir> <test_wsm6_column> [case ...]
# A case is nCells:nVertLevels:nsteps:mp_top_level:dt:seed[:d]; ':d' runs it with l_diags (the reflectivity and
# relhum of the step before a diagnostics alarm, harness.F90). The default list covers 1, 2, 3, 7 and 64 cells
# (1-trip loops, one pair, odd tails, long vector bodies), RG's dt = 600 s (5 sub-steps of 120 s) and dt = 90 s
# (one), 55 levels as RG, and mp_top_level < nVertLevels.
set -euo pipefail
mpas=$1; mmm=$2; work=$3; cpp=$4; shift 4
here=$(dirname "$(realpath "$0")")
cases=${*:-"1:55:2:55:600:1 2:55:2:55:600:2 3:55:3:55:600:3 7:55:3:55:600:4 64:55:3:55:600:5 64:55:2:55:90:6 9:41:2:33:600:7 1:55:2:55:600:11:d 2:55:2:55:600:12:d 7:55:3:55:600:14:d 64:55:3:55:600:15:d 64:56:2:56:90:16:d 9:41:2:33:600:17:d"}
mkdir -p "$work"
"$here/build_ref.sh" "$mpas" "$mmm" "$work/ref" > "$work/build_ref.log" 2>&1 || { cat "$work/build_ref.log"; echo "FAIL: reference build"; exit 1; }
head -1 "$work/build_ref.log"
rc=0
for c in $cases; do
  IFS=: read -r nc nv ns top dt seed d <<< "$c"
  diag=""; [ "${d:-}" = d ] && diag=l_diags
  in=$work/in_$c.bin
  rm -f "$in" "$work/ref_$c.bin" "$work/cpp_$c.bin"  # no file of an earlier run can stand in for this one's
  python3 -I "$here/columns.py" make "$in" "$nc" "$nv" "$ns" "$top" "$dt" "$seed" || { echo "FAIL: columns $c"; rc=1; continue; }
  "$work/ref/wsm6_column_ref" "$in" "$work/ref_$c.bin" $diag || { echo "FAIL: reference run $c"; rc=1; continue; }
  "$cpp" "$in" "$work/cpp_$c.bin" $diag || { echo "FAIL: CPPMPAS run $c"; rc=1; continue; }
  echo -n "case $c: "
  python3 -I "$here/columns.py" compare "$in" "$work/ref_$c.bin" "$work/cpp_$c.bin" $diag || rc=1
done
[ $rc = 0 ] && echo "wsm6 column check: PASS" || echo "wsm6 column check: FAIL"
exit $rc
