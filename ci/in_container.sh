#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Run a command in the suite's container (image cppmpas-ci, ci/Dockerfile) as the runner's user, with the fixed
# paths CTest stores: /w/src (the CPPMPAS checkout), /w/build, /w/data, /w/out, /w/ci (these scripts), /w/plan,
# /w/mvec (cppmpas-data ci/mvec_intel).
#   in_container.sh <base dir: holds build/, data/, out/> <command...>
set -euo pipefail
base=$1; shift
ws=${GITHUB_WORKSPACE:-$PWD}
mkdir -p "$base"/build "$base"/data "$base"/out "$ws"/plan "$ws"/mvec-repo/ci/mvec_intel
exec docker run --rm --user "$(id -u):$(id -g)" -e HOME=/tmp \
  -v "$ws/src:/w/src" -v "$base/build:/w/build" -v "$base/data:/w/data" -v "$base/out:/w/out" \
  -v "$ws/tools/ci:/w/ci:ro" -v "$ws/plan:/w/plan:ro" -v "$ws/mvec-repo/ci/mvec_intel:/w/mvec:ro" \
  cppmpas-ci "$@"
