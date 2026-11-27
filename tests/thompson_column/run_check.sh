#!/bin/bash
# SPDX-License-Identifier: GPL-3.0-or-later
# The Thompson column check: MPAS's own Thompson path (build_ref.sh, gfortran with the reference flags) and CPPMPAS's
# (test_thompson_column) on the same synthetic columns, compared bit for bit, and the module variables and tables
# thompson_init leaves, compared record by record.
#   run_check.sh <mpas src> <mmm-physics a4baf7f tree> <Thompson tables dir> <work dir> <test_thompson_column> [case ...]
# The tables dir holds the reference run's MP_THOMPSON_{QRacrQG,QRacrQS,freezeH2O,QIautQS}_DATA.DBL, checked by sha256
# against tables.sha256 (the tables/ of
# mpas842-dbl-rgcp-001). A case is nCells:nVertLevels:nsteps:mp_top_level:dt:seed[:d]; ':d' runs it with l_diags
# (the reflectivity and relhum of the step before a diagnostics alarm). The default list covers 1, 2, 3, 7, 14 and 64
# cells, 55 levels as RG and odd / even level counts (the rain-slope and dBZ loops' libmvec pairs from kte), RG's
# dt = 600 s (n_microp = 7 sub-steps of 600/7 s), 90 s (one), 300 s (three) and 45 s, and mp_top_level < nVertLevels.
set -u
mpas=$1; mmm=$2; tables=$3; work=$4; cpp=$5; shift 5
here=$(dirname "$(realpath "$0")")
cases=${*:-"1:55:2:55:600:1 2:55:2:55:600:2 3:56:3:56:600:3 7:55:3:55:600:4 14:55:3:55:300:5 64:55:3:55:600:6 64:56:2:56:90:7 9:41:2:33:600:8 7:30:4:30:45:9 1:55:2:55:600:11:d 2:56:2:56:600:12:d 7:55:3:55:600:14:d 64:55:3:55:600:15:d 64:56:2:56:90:16:d 9:41:2:33:600:17:d 7:55:2:55:130:18 14:55:2:55:130:19:d 7:70:2:70:600:20 14:70:2:70:300:21:d"}
for t in QRacrQG QRacrQS freezeH2O QIautQS; do
  [ -f "$tables/MP_THOMPSON_${t}_DATA.DBL" ] || { echo "FAIL: no $tables/MP_THOMPSON_${t}_DATA.DBL"; exit 1; }
done
# the tables must be the reference run's (release mpas842-dbl-rgcp-001, mpas842-dbl-rgcp-001.tables.sha256)
(cd "$tables" && sha256sum --quiet -c "$here/tables.sha256") || { echo "FAIL: $tables are not mpas842-dbl-rgcp-001's tables"; exit 1; }
mkdir -p "$work"
"$here/build_ref.sh" "$mpas" "$mmm" "$work/ref" > "$work/build_ref.log" 2>&1 || { cat "$work/build_ref.log"; echo "FAIL: reference build"; exit 1; }
head -1 "$work/build_ref.log"
mkdir -p "$work/run"
for t in QRacrQG QRacrQS freezeH2O QIautQS; do ln -sf "$(realpath "$tables/MP_THOMPSON_${t}_DATA.DBL")" "$work/run/"; done
[ -n "${cases// /}" ] || { echo "FAIL: no case"; exit 1; }
rm -f "$work/ref_init.bin" "$work/cpp_init.bin"
rc=0
first=1
ncase=0
for c in $cases; do
  IFS=: read -r nc nv ns top dt seed d <<< "$c"
  diag=""; [ "${d:-}" = d ] && diag=l_diags
  in=$work/in_$c.bin; rm -f "$work/ref_$c.bin" "$work/cpp_$c.bin"
  python3 -I "$here/columns.py" make "$in" "$nc" "$nv" "$ns" "$top" "$dt" "$seed" || { echo "FAIL: columns $c"; rc=1; continue; }
  (cd "$work/run" && "$work/ref/thompson_column_ref" "$in" "$work/ref_$c.bin" "$work/ref_init.bin" $diag) || { echo "FAIL: reference run $c"; rc=1; continue; }
  "$cpp" "$in" "$work/cpp_$c.bin" "$work/cpp_init.bin" "$tables" $diag || { echo "FAIL: CPPMPAS run $c"; rc=1; continue; }
  if [ $first = 1 ]; then
    python3 -I "$here/columns.py" compare_init "$work/ref_init.bin" "$work/cpp_init.bin" || rc=1
    first=0
  fi
  ncase=$((ncase + 1))
  echo -n "case $c: "
  python3 -I "$here/columns.py" compare "$in" "$work/ref_$c.bin" "$work/cpp_$c.bin" $diag || rc=1
  rm -f "$work/ref_$c.bin" "$work/cpp_$c.bin"
done
[ $first = 0 ] || { echo "FAIL: thompson_init was not compared"; rc=1; }
[ $ncase -ge 1 ] || rc=1
[ $rc = 0 ] && echo "thompson column check: PASS" || echo "thompson column check: FAIL"
exit $rc
