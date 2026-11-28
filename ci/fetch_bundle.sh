#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# The source of a commit that is not on GitHub (.github/workflows/cpu_suite.yml, input `bundle`): tools/ci/ci_submit.sh
# uploads a git bundle of the commits origin/main lacks to the `jobs` release of Xassemblianist/cppmpas-data. This
# downloads it, checks its SHA-256 (the one ci_submit.sh computed before the upload), verifies it against the clone
# (its prerequisites are commits of main, so the clone needs main's full history) and checks the commit out.
#   fetch_bundle.sh <clone of CPPMPAS with main's history> <asset name> <sha256> <commit>
set -euo pipefail
src=$1 asset=$2 sha=$3 commit=$4
DATA_REPO=Xassemblianist/cppmpas-data
[[ $asset =~ ^[A-Za-z0-9._-]+\.bundle$ ]] || { echo "fetch_bundle.sh: bad asset name '$asset'" >&2; exit 2; }
[[ $sha =~ ^[0-9a-f]{64}$ ]] || { echo "fetch_bundle.sh: bad sha256 '$sha'" >&2; exit 2; }
[[ $commit =~ ^[0-9a-f]{40}$ ]] || { echo "fetch_bundle.sh: the commit must be a full sha, not '$commit'" >&2; exit 2; }
dir=$(mktemp -d)
curl -fsSL --retry 5 -o "$dir/$asset" "https://github.com/$DATA_REPO/releases/download/jobs/$asset"
echo "$sha  $dir/$asset" | sha256sum --strict -c -
git -C "$src" bundle verify "$dir/$asset"
git -C "$src" fetch -q "$dir/$asset" "refs/heads/ci:refs/heads/ci-under-test"
git -C "$src" checkout -q --detach "$commit"
[ "$(git -C "$src" rev-parse HEAD)" = "$commit" ] || { echo "fetch_bundle.sh: HEAD is not $commit" >&2; exit 1; }
[ "$(git -C "$src" rev-parse refs/heads/ci-under-test)" = "$commit" ] ||
  { echo "fetch_bundle.sh: the bundle holds $(git -C "$src" rev-parse refs/heads/ci-under-test), not $commit" >&2; exit 1; }
rm -rf "$dir"
