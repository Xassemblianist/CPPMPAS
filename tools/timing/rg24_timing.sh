#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Time RG24 (mpas842-dbl-rg24-001: 24 h, x1.10242, mesoscale_reference) with CPPMPAS and with MPAS 8.4.2, one run at
# a time, and compare every run's 19 stream files with the reference at tolerance 0 (tests/streams_compare.py).
#
# usage: tools/timing/rg24_timing.sh <set>...   sets: gpu, cpu (CPPMPAS), mpas1, mpasN (MPAS on N MPI ranks)
# environment:
#   SRC      this source tree                         BUILD   the CPPMPAS build (cppmpas, cppmpas_gpu)
#   REF      .../mpas842-dbl-rg24-001/RG24            WORK    scratch: runs/, timing/, compare/
#   TABLES   the physics tables (as CPPMPAS_MPAS_TABLES)   RRTMG  CPPWRF's rrtmg tables (as CPPMPAS_CPPWRF_RRTMG)
#   PRELOAD  the Intel libmvec preload of the reference runs (cppmpas-data ci/mvec_intel)
#   MPAS_EXE atmosphere_model, built as cppmpas-data ci/build.sh builds it; MPAS_TABLES its *TBL *DATA* files
#   MPIEXEC  mpiexec; MESH_PARTS  the directory with x1.10242.graph.info.part.N (NCAR's x1.10242.tar.gz)
#   PYTHON   a python3 with netCDF4 and numpy;   NRUNS  runs per set (3)
# Each run writes WORK/timing/<set>-<i>.json (tools/timing/timed.py) and WORK/compare/<set>-<i>.log; the outputs are
# deleted after the comparison. MPAS runs as ci/ref.py runs it: mpiexec, the preload, OMPI_MCA_pml=ob1, the stack
# limit raised; one rank uses only the self transport, more ranks also the shared-memory one, bound to cores.
set -uo pipefail
: "${SRC:?}" "${BUILD:?}" "${REF:?}" "${WORK:?}" "${TABLES:?}" "${RRTMG:?}"
PY=${PYTHON:-python3}
N=${NRUNS:-3}
export LC_ALL=C TZ=UTC OMP_NUM_THREADS=1
ulimit -s unlimited 2>/dev/null || ulimit -s "$(ulimit -Hs)"
mkdir -p "$WORK/runs" "$WORK/timing" "$WORK/compare"
mpas_dir() {
  rm -rf "$1" && mkdir -p "$1"
  cp "$REF"/forecast/namelist.atmosphere "$REF"/forecast/streams.atmosphere "$REF"/forecast/stream_list.atmosphere.* "$1"/
  ln -s "$REF/init/x1.10242.init.nc" "$1/x1.10242.init.nc"
  for t in "${MPAS_TABLES:?}"/*; do ln -s "$t" "$1"/; done
  cp "${MESH_PARTS:?}"/x1.10242.graph.info.part.* "$1"/
}
for set in "$@"; do
  for i in $(seq 1 "$N"); do
    name=$set-$i R=$WORK/runs/$set-$i
    case $set in
      gpu|cpu)
        rm -rf "$R" && mkdir -p "$R"
        if [ "$set" = gpu ]; then X=$BUILD/cppmpas_gpu pre=; else X=$BUILD/cppmpas pre=${PRELOAD:-}; fi
        LD_PRELOAD=$pre "$PY" -I "$SRC/tools/timing/timed.py" "$name" "$R" "$WORK/timing" \
          "$X" run "$REF/init/x1.10242.init.nc" "$R" --namelist "$REF/forecast/namelist.atmosphere" \
          --streams "$REF/forecast/streams.atmosphere" --scalars qv,qc,qr,qi,qs,qg --every 0 --tables "$TABLES" \
          --rrtmg "$RRTMG" ;;
      mpas[0-9]*)
        n=${set#mpas}; mpas_dir "$R"
        btl=self; [ "$n" -gt 1 ] && btl=self,sm
        hw=; [ "$n" -gt "$(nproc --all)" ] && { echo "$n ranks > $(nproc --all) CPUs"; exit 2; }
        [ "$n" -gt "$(( $(nproc --all) / 2 ))" ] && hw=--use-hwthread-cpus
        OMPI_MCA_pml=ob1 OMPI_MCA_btl=$btl LD_PRELOAD=${PRELOAD:?} "$PY" -I "$SRC/tools/timing/timed.py" "$name" "$R" \
          "$WORK/timing" "${MPIEXEC:-mpiexec}" --allow-run-as-root -n "$n" $hw --bind-to core \
          -x LD_PRELOAD -x OMPI_MCA_pml -x OMPI_MCA_btl "${MPAS_EXE:?}" ;;
      *) echo "unknown set $set"; exit 2 ;;
    esac
    cp "$R/log.atmosphere.0000.out" "$WORK/timing/$name.log.atmosphere.0000.out" 2>/dev/null
    rm -f "$R/x1.10242.init.nc"
    "$PY" -I "$SRC/tests/streams_compare.py" "$SRC/tools/mpas_compare.py" "$REF/forecast" "$R" \
      "$WORK/compare/$name.json" 19 2168 > "$WORK/compare/$name.log" 2>&1
    echo "$name compare rc=$?"
    rm -f "$R"/*.nc
  done
done
