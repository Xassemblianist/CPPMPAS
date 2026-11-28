#!/bin/bash
# SPDX-License-Identifier: GPL-3.0-or-later
# The Noah column check: MPAS's own Noah path (build_ref.sh, gfortran with the reference flags) and CPPMPAS's
# (test_noah_column) on the same synthetic points, compared bit for bit, tables included.
#   run_check.sh <mpas src> <table dir> <work dir> <test_noah_column> [case ...]
# <table dir> holds VEGPARM.TBL, SOILPARM.TBL and GENPARM.TBL of the reference build (the mpas-build asset's
# tables/, from MPAS-Data c57dbc7). A case is nCells:nsteps:dt:seed:landuse:flags (columns.py); the default list has
# 1, 2, 3, 12, 13, 24, 25, 65, 72 and 145 points (1-trip loops, pairs, odd tails, every regime, the 6 sea-ice fractions), RG's
# dt = 600 s and dt = 90 s, MODIS and USGS land use, and the options config_sfc_albedo = false,
# config_sfc_snowalbedo = false, config_frac_seaice = false and config_do_restart.
set -euo pipefail
mpas=$1; tables=$(realpath "$2"); work=$3; cpp=$4; shift 4
here=$(dirname "$(realpath "$0")")
cases=${*:-"1:2:600:1:MODIS:- 2:2:600:2:MODIS:- 3:3:600:3:MODIS:- 12:3:600:4:MODIS:- 13:4:600:5:MODIS:- 72:3:600:6:MODIS:- 145:4:600:12:MODIS:- 65:3:90:7:MODIS:- 25:3:600:8:USGS:- 24:3:600:9:MODIS:as 72:3:600:10:MODIS:f 24:2:600:11:MODIS:r"}
mkdir -p "$work"; work=$(realpath "$work")
"$here/build_ref.sh" "$mpas" "$work/ref" > "$work/build_ref.log" 2>&1 || { cat "$work/build_ref.log"; echo "FAIL: reference build"; exit 1; }
head -1 "$work/build_ref.log"
for t in VEGPARM.TBL SOILPARM.TBL GENPARM.TBL; do
  [ -f "$tables/$t" ] || { echo "FAIL: $tables/$t missing"; exit 1; }
  ln -sf "$tables/$t" "$work/$t"
done
rc=0
for c in $cases; do
  IFS=: read -r nc ns dt seed lu fl <<< "$c"
  in=$work/in_$c.bin
  rm -f "$in" "$work/ref_$c.bin" "$work/cpp_$c.bin"  # no file of an earlier run can stand in for this one's
  python3 -I "$here/columns.py" make "$in" "$nc" "$ns" "$dt" "$seed" "$lu" "$fl" || { echo "FAIL: columns $c"; rc=1; continue; }
  (cd "$work" && ./ref/noah_column_ref "$in" "$work/ref_$c.bin") || { echo "FAIL: reference run $c"; rc=1; continue; }
  "$cpp" "$in" "$work/cpp_$c.bin" "$tables" || { echo "FAIL: CPPMPAS run $c"; rc=1; continue; }
  echo -n "case $c: "
  python3 -I "$here/columns.py" compare "$in" "$work/ref_$c.bin" "$work/cpp_$c.bin" || rc=1
done
[ $rc = 0 ] && echo "noah column check: PASS" || echo "noah column check: FAIL"
exit $rc
