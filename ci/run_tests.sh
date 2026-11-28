#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Inside the suite's container: run one job's tests (a ci/plan.py list) in /w/build, under the Intel libmvec preload
# as on the LAN PC (glibc's own libmvec gives other vector log / log10 lanes on AMD than on Intel; the references
# were made with this preload, docs/REFERENCES.md).
#   run_tests.sh <tests file> [gpu]
# gpu: without the fixture setups (ctest -FA): every GPU test of this CPU build exits 77 (Skipped) before it reads
# anything. Writes /w/out/ctest.log, junit.xml, LastTest.log. Exit status: ctest's; 6 when the preload is missing,
# is not mapped into a process started with it, or ld.so reported it could not preload it (ld.so only warns and goes
# on, and on an Intel runner glibc's own lanes are the reference's, so a missing preload would pass there silently).
set -uo pipefail
list=$1; kind=${2:-cpu}
mkdir -p /w/out
flags=$(grep -m1 '^flags' /proc/cpuinfo)
for f in fma avx2; do
  # glibc's ifuncs must pick the FMA / AVX2 libm variants of the reference (src/core/libm.hpp)
  [[ " $flags " == *" $f "* ]] || { echo "this runner's CPU has no $f" >&2; exit 5; }
done
grep -m1 'model name' /proc/cpuinfo | tee /w/out/cpu.txt
getconf _NPROCESSORS_ONLN | tee -a /w/out/cpu.txt
so=/w/build/ci/libmvec_intel.so
[ -s "$so" ] || { echo "run_tests.sh: the preload $so is missing" >&2; exit 6; }
# a process started under the preload maps it (/proc/self/maps of that process), and ld.so says nothing
LD_PRELOAD=$so cat /proc/self/maps > /w/out/preload_maps.txt 2> /w/out/preload_err.txt
if ! grep -qF "$so" /w/out/preload_maps.txt || [ -s /w/out/preload_err.txt ]; then
  echo "run_tests.sh: the preload $so is not loaded:" >&2; cat /w/out/preload_err.txt >&2; exit 6
fi
echo "preload $so loaded ($(sha256sum "$so" | cut -c1-16))" | tee -a /w/out/cpu.txt
extra=()
[ "$kind" = gpu ] && extra=(-FA '.*')
cd /w/build
start=$(date +%s)
LD_PRELOAD=$so ctest -j "$(getconf _NPROCESSORS_ONLN)" --tests-from-file "$list" "${extra[@]}" \
  --output-on-failure --output-junit /w/out/junit.xml > /w/out/ctest.log 2>&1
rc=$?
echo "ctest exit status $rc, $(( $(date +%s) - start )) s" | tee -a /w/out/ctest.log
cp Testing/Temporary/LastTest.log /w/out/ 2>/dev/null
# LastTest.log holds every test's output, ctest.log the failing ones'
if grep -qs 'ld\.so: object .*LD_PRELOAD' /w/out/ctest.log /w/out/LastTest.log; then
  echo "run_tests.sh: ld.so could not load a preload in a test:" >&2
  grep -hs -m5 'ld\.so: object .*LD_PRELOAD' /w/out/ctest.log /w/out/LastTest.log >&2
  rc=6
fi
grep -E 'tests passed|tests failed|\((Failed|Not Run|Skipped|Timeout)\)|Total Test time' /w/out/ctest.log | head -80
exit $rc
