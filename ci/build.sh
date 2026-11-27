#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Inside the suite's container (ci/Dockerfile): configure CPPMPAS (the CPU build, Release, as the LAN PC suite did)
# with every data location under /w/data (ci/data.sh extract), and optionally build it and the Intel libmvec preload.
#   build.sh configure     configure /w/src into /w/build, list the tests (/w/out/tests.json)
#   build.sh build         configure, build, and build the preload /w/build/ci/libmvec_intel.so from /w/mvec
#                          (Xassemblianist/cppmpas-data ci/mvec_intel, as its ci/ref.py build_preload() builds it)
# The paths are fixed: CTest stores them, and every test job mounts the same directories at the same places.
set -euo pipefail
mode=${1:-}
D=/w/data
mkdir -p /w/out
{ uname -m; grep -m1 'model name' /proc/cpuinfo; ldd --version | sed -n 1p; gcc --version | sed -n 1p
  gfortran --version | sed -n 1p; cmake --version | sed -n 1p; python3 --version; } > /w/out/versions.txt
cmake -S /w/src -B /w/build -DCMAKE_BUILD_TYPE=Release \
  -DCPPMPAS_CASES_DIR=$D/cases -DCPPMPAS_REFS_DIR=$D/refs \
  -DCPPMPAS_MPAS_REGISTRY=$D/mpas_src/src/core_atmosphere/Registry.xml -DCPPMPAS_MPAS_SRC=$D/mpas_src \
  -DCPPMPAS_MMM_PHYSICS=$D/mmm -DCPPMPAS_MPAS_TABLES=$D/mpasbuild/tables -DCPPMPAS_NOAH_TABLES=$D/mpasbuild/tables \
  -DCPPMPAS_CPPWRF_RRTMG=$D/rrtmg > /w/out/cmake.log 2>&1 || { tail -40 /w/out/cmake.log; exit 3; }
(cd /w/build && ctest --show-only=json-v1 > /w/out/tests.json)
case $mode in
configure) ;;
build)
  cmake --build /w/build -j "$(nproc)" > /w/out/build.log 2>&1 || { grep -m40 -E 'error|Error' /w/out/build.log; exit 4; }
  (cd /w/build && ctest --show-only=json-v1 > /w/out/tests.json)  # now with the test executables' paths
  (cd /w/mvec && sha256sum --quiet -c SHA256SUMS)
  mkdir -p /w/build/ci
  g++ -std=c++20 -O2 -ffp-contract=off -fPIC -shared -I /w/mvec /w/mvec/tests/phys/libmvec_intel.cpp \
    -o /w/build/ci/libmvec_intel.so
  sha256sum /w/build/ci/libmvec_intel.so >> /w/out/versions.txt ;;
*) sed -n 3,7p "$0"; exit 2 ;;
esac
cat /w/out/versions.txt
