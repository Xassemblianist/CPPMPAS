#!/bin/bash
# SPDX-License-Identifier: GPL-3.0-or-later
# The MYNN surface layer column check: MPAS's own MPAS_to_physics, init_sfclayer and driver_sfclayer with
# config_sfclayer_scheme = 'sf_mynn' (build_ref.sh, gfortran with the reference flags) and CPPMPAS's
# (test_sfmynn_column) on the same synthetic columns, compared bit for bit at every step.
#   run_check.sh <mpas src> <mmm-physics a4baf7f tree> <work dir> <test_sfmynn_column> [case ...]
# A case is nCellsSolve:nVertLevels:nsteps:seed:frac_seaice:do_restart:len_disp (columns.py). The default list has
# 1, 2, 3, 7, 10, 21, 30, 63 and 64 cells (the 1-trip loop, one lane pair, odd tails of the vectorised loop), 55
# levels as RG and 41 / 56, config_frac_seaice true (the convection_permitting reference) and false,
# config_do_restart true (initflag 0 at the first step), dt 600 s, len_disp 240 km (RG) and 3 km (no sub-grid wind).
set -euo pipefail
mpas=$1; mmm=$2; work=$3; cpp=$4; shift 4
here=$(dirname "$(realpath "$0")")
[ -d "$mpas/src/core_atmosphere/physics" ] || { echo "FAIL: no MPAS v8.4.2 source tree at '$mpas' (CPPMPAS_MPAS_SRC)"; exit 1; }
[ -f "$mmm/sf_mynn.F90" ] || { echo "FAIL: no MMM-physics a4baf7f tree at '$mmm' (CPPMPAS_MMM_PHYSICS)"; exit 1; }
cases=${*:-"1:55:3:1:1:0:240000 2:55:3:2:1:0:240000 3:56:3:3:1:0:240000 7:55:4:4:1:0:240000 64:55:4:5:1:0:240000 64:41:3:6:0:0:240000 21:55:3:7:1:1:240000 30:55:3:8:1:0:3000 63:56:3:9:0:0:240000 10:55:2:10:0:1:3000"}
mkdir -p "$work"
"$here/build_ref.sh" "$mpas" "$mmm" "$work/ref" > "$work/build_ref.log" 2>&1 || { cat "$work/build_ref.log"; echo "FAIL: reference build"; exit 1; }
head -1 "$work/build_ref.log"
[ -n "${cases// /}" ] || { echo "FAIL: no case"; exit 1; }
rc=0
ncase=0
"$cpp" guards || { echo "FAIL: the sfmynn hard errors"; rc=1; }
for c in $cases; do
  IFS=: read -r nc nv ns seed si rs ld <<< "$c"
  in=$work/in_$c.bin
  rm -f "$in" "$work/ref_$c.bin" "$work/cpp_$c.bin" "$work/ref_$c.log"  # no file of an earlier run can stand in
  python3 -I "$here/columns.py" make "$in" "$nc" "$nv" "$ns" "$seed" "$si" "$rs" "$ld" || { echo "FAIL: columns $c"; rc=1; continue; }
  "$work/ref/sfmynn_column_ref" "$in" "$work/ref_$c.bin" > "$work/ref_$c.log" 2>&1 || { cat "$work/ref_$c.log"; echo "FAIL: reference run $c"; rc=1; continue; }
  "$cpp" "$in" "$work/cpp_$c.bin" || { echo "FAIL: CPPMPAS run $c"; rc=1; continue; }
  ncase=$((ncase + 1))
  echo -n "case $c: "
  if python3 -I "$here/columns.py" compare "$in" "$work/ref_$c.bin" "$work/cpp_$c.bin"; then
    rm -f "$in" "$work/ref_$c.bin" "$work/cpp_$c.bin"  # kept when they differ
  else
    rc=1
  fi
done
[ $ncase -ge 1 ] || { echo "FAIL: no case ran"; rc=1; }
[ $rc = 0 ] && echo "sfmynn column check: PASS" || echo "sfmynn column check: FAIL"
exit $rc
