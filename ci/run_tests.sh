#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Inside the suite's container: run one job's tests (a ci/plan.py list) in /w/build, under the Intel libmvec preload
# as on the LAN PC (glibc's own libmvec gives other vector log / log10 lanes on AMD than on Intel; the references
# were made with this preload, docs/REFERENCES.md).
#   run_tests.sh <tests file> [gpu]
# gpu: without the fixture setups (ctest -FA): every GPU test of this CPU build exits 77 (Skipped) before it reads
# anything. Writes /w/out/ctest.log, junit.xml, LastTest.log. Exit status: ctest's.
set -uo pipefail
list=$1; kind=${2:-cpu}
mkdir -p /w/out
flags=$(grep -m1 '^flags' /proc/cpuinfo)
for f in fma avx2; do
  # glibc's ifuncs must pick the FMA / AVX2 libm variants of the reference (src/core/libm.hpp)
  [[ " $flags " == *" $f "* ]] || { echo "this runner's CPU has no $f" >&2; exit 5; }
done
grep -m1 'model name' /proc/cpuinfo | tee /w/out/cpu.txt
nproc | tee -a /w/out/cpu.txt
extra=()
[ "$kind" = gpu ] && extra=(-FA '.*')
cd /w/build
start=$(date +%s)
LD_PRELOAD=/w/build/ci/libmvec_intel.so ctest -j "$(nproc)" --tests-from-file "$list" "${extra[@]}" \
  --output-on-failure --output-junit /w/out/junit.xml > /w/out/ctest.log 2>&1
rc=$?
echo "ctest exit status $rc, $(( $(date +%s) - start )) s" | tee -a /w/out/ctest.log
cp Testing/Temporary/LastTest.log /w/out/ 2>/dev/null
grep -E 'tests passed|tests failed|\((Failed|Not Run|Skipped|Timeout)\)|Total Test time' /w/out/ctest.log | head -80
exit $rc
