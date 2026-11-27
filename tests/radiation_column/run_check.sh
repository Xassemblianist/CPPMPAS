#!/bin/bash
# SPDX-License-Identifier: GPL-3.0-or-later
# The radiation column check: MPAS's own init_o3climatology, rrtmg_initlw_forMPAS / rrtmg_initsw_forMPAS,
# o3climatology_from_MPAS, driver_cloudiness ('cld_fraction'), driver_radiation_sw ('rrtmg_sw') and
# driver_radiation_lw ('rrtmg_lw') (build_ref.sh, gfortran with the reference flags) and CPPMPAS's
# (test_radiation_column) on the same synthetic columns, compared bit for bit.
#   run_check.sh <mpas src> <MPAS build tables dir> <work dir> <test_radiation_column> [case ...]
# The tables dir holds the RRTMG_LW_DATA.DBL, RRTMG_SW_DATA.DBL and OZONE_*.TBL files of the reference build
# (cppmpas-data mpas-build asset, bin/../tables). The RRTMG tables MPAS's initialisation computes from them are dumped
# (radiation_dump_tables) and compared byte for byte with the shared CPPWRF dumps (CPPMPAS_CPPWRF_RRTMG, required:
# cppmpas run reads those, so without them the check fails); test_radiation_column reads the dumps. A case is nCellsSolve:nVertLevels:nsteps:seed:cldovrlp:dcorrlen:o3climatology
# (columns.py).
set -euo pipefail
mpas=$1; tables=$2; work=$(realpath -m "$3"); cpp=$(realpath "$4"); shift 4
here=$(dirname "$(realpath "$0")")
cases=${*:-"1:55:6:1:2:0:1 2:55:6:2:2:0:1 3:55:6:3:2:0:1 7:55:6:4:2:0:1 16:55:6:5:2:0:1 64:55:6:6:2:0:1 63:55:3:7:2:0:1 9:41:6:8:2:0:0 8:55:3:9:1:0:1 8:55:3:10:3:0:1 8:55:3:11:5:1:1 8:55:3:12:0:0:1"}
mkdir -p "$work"
"$here/build_ref.sh" "$mpas" "$work/ref" > "$work/build_ref.log" 2>&1 || { cat "$work/build_ref.log"; echo "FAIL: reference build"; exit 1; }
head -1 "$work/build_ref.log"
mkdir -p "$work/run"
for f in RRTMG_LW_DATA.DBL RRTMG_SW_DATA.DBL OZONE_DAT.TBL OZONE_LAT.TBL OZONE_PLEV.TBL; do
  [ -f "$tables/$f" ] || { echo "FAIL: no $tables/$f"; exit 1; }
  ln -sf "$(realpath "$tables/$f")" "$work/run/$f"
done
rc=0
rm -f "$work/run/lw_tables.bin" "$work/run/sw_tables.bin" "$work/run/dump.log"
(cd "$work/run" && "$work/ref/radiation_dump_tables" lw_tables.bin sw_tables.bin > dump.log 2>&1) || { cat "$work/run/dump.log"; echo "FAIL: table dump"; exit 1; }
for b in lw sw; do
  [ -s "$work/run/${b}_tables.bin" ] || { cat "$work/run/dump.log"; echo "FAIL: the table dump wrote no ${b}_tables.bin"; exit 1; }
done
if [ -z "${CPPMPAS_CPPWRF_RRTMG:-}" ]; then
  echo "FAIL: CPPMPAS_CPPWRF_RRTMG is not set: the RRTMG tables cppmpas run reads (CPPWRF's data/rrtmg) are not checked"
  rc=1
else
  for b in lw sw; do
    if cmp "$work/run/${b}_tables.bin" "$CPPMPAS_CPPWRF_RRTMG/rrtmg_${b}_tables.bin"; then
      echo "RRTMG ${b} tables: MPAS's initialisation = CPPWRF's rrtmg_${b}_tables.bin"
    else
      echo "FAIL: RRTMG ${b} tables differ from CPPWRF's (or are missing)"; rc=1
    fi
  done
fi
for c in $cases; do
  IFS=: read -r nc nv ns seed ov dc oz <<< "$c"
  in=$work/in_$c.bin
  rm -f "$in" "$work/ref_$c.bin" "$work/cpp_$c.bin" "$work/ref_$c.log"  # no file of an earlier run can stand in
  python3 -I "$here/columns.py" make "$in" "$nc" "$nv" "$ns" "$seed" "$ov" "$dc" "$oz" || { echo "FAIL: columns $c"; rc=1; continue; }
  (cd "$work/run" && "$work/ref/radiation_column_ref" "$in" "$work/ref_$c.bin" > "$work/ref_$c.log" 2>&1) || { cat "$work/ref_$c.log"; echo "FAIL: reference run $c"; rc=1; continue; }
  "$cpp" "$in" "$work/cpp_$c.bin" "$tables" "$work/run/lw_tables.bin" "$work/run/sw_tables.bin" || { echo "FAIL: CPPMPAS run $c"; rc=1; continue; }
  echo -n "case $c: "
  if python3 -I "$here/columns.py" compare "$in" "$work/ref_$c.bin" "$work/cpp_$c.bin"; then
    rm -f "$in" "$work/ref_$c.bin" "$work/cpp_$c.bin"  # kept when they differ
  else
    rc=1
  fi
done
[ $rc = 0 ] && echo "radiation column check: PASS" || echo "radiation column check: FAIL"
exit $rc
