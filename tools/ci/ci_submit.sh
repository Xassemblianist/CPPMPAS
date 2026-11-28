#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# The CPU suite (.github/workflows/cpu_suite.yml) of any local commit, without pushing it: a commit that is not on
# origin/main travels as a git bundle (only the objects origin/main lacks) uploaded under a unique name to the `jobs`
# release of Xassemblianist/cppmpas-data; the workflow's build job fetches it, checks its SHA-256 and checks the
# commit out (ci/fetch_bundle.sh). The asset is deleted again whatever the outcome (also on an error or Ctrl-C).
#
#   tools/ci/ci_submit.sh [COMMIT] [--labels short|long|all] [--no-wait] [--hours H] [--out DIR]
#       [--workflow-ref REF] [--max-parallel N]
#
# COMMIT: default HEAD. --labels: the tests to run (default all; short: without the label long). The script waits
# for the run (at most --hours, default 3), downloads its summary (the artifact "summary") to --out (default
# ~/cppwrf_scratch/ci/cppmpas/<run id>) and prints it (per label: tests, passed, skipped, failed, not run). Exit 0
# only when the run succeeded and its summary says PASSED. --no-wait: wait only until the build job has fetched the
# bundle (then delete it), print the run's URL and exit. --workflow-ref: the branch whose cpu_suite.yml and ci/ run
# (default main; the commit under test is only the source that is built and tested). --max-parallel: test jobs at
# once (default 10, at most 10; fewer when other jobs of the project run: its share is 20 runners).
# Needs gh, logged in with the rights to start CPPMPAS's workflows and upload to cppmpas-data's releases.
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
GH_REPO_SRC=${CPPMPAS_CI_REPO:-Xassemblianist/CPPMPAS}
DATA=Xassemblianist/cppmpas-data
COMMIT=HEAD LABELS=all WAIT=1 HOURS=3 OUT= WREF=main MAXP=10
while [ $# -gt 0 ]; do
  case $1 in
    --labels) LABELS=$2; shift 2 ;; --no-wait) WAIT=0; shift ;; --hours) HOURS=$2; shift 2 ;;
    --out) OUT=$2; shift 2 ;; --workflow-ref) WREF=$2; shift 2 ;; --max-parallel) MAXP=$2; shift 2 ;;
    -*) sed -n '3,19p' "$0"; exit 2 ;;
    *) COMMIT=$1; shift ;;
  esac
done
case $LABELS in all|short|long) ;; *) echo "[ci_submit] --labels: all, short or long"; exit 2 ;; esac
[[ $MAXP =~ ^[0-9]+$ ]] && [ "$MAXP" -ge 1 ] && [ "$MAXP" -le 10 ] || { echo "[ci_submit] --max-parallel: 1 to 10"; exit 2; }
SHA=$(git -C "$REPO" rev-parse --verify "$COMMIT^{commit}")
git -C "$REPO" fetch -q origin main 2> /dev/null || echo "[ci_submit] note: git fetch origin failed; origin/main as known locally"
BASE=$(git -C "$REPO" rev-parse origin/main)
TAG="ci-$(date +%Y%m%d-%H%M%S)-${SHA:0:8}-$(od -An -N2 -tx1 /dev/urandom | tr -d ' \n')"
BUNDLE= BSHA= TMPD=
cleanup() {
  [ -z "$TMPD" ] || rm -rf "$TMPD"
  if [ -n "$BUNDLE" ]; then
    if gh release delete-asset jobs "$BUNDLE" -R "$DATA" -y > /dev/null 2>&1; then
      echo "[ci_submit] bundle $BUNDLE deleted"
    else
      echo "[ci_submit] WARNING: could not delete $BUNDLE: gh release delete-asset jobs $BUNDLE -R $DATA -y"
    fi
    BUNDLE=
  fi
}
trap cleanup EXIT
trap 'exit 130' INT TERM
if ! git -C "$REPO" merge-base --is-ancestor "$SHA" "$BASE"; then
  # git bundle's own format (v2), without touching any ref: the boundary commits (on origin/main) as prerequisites,
  # the commit as refs/heads/ci, then a thin pack of the objects origin/main lacks
  TMPD=$(mktemp -d)
  F=$TMPD/$TAG.bundle
  { echo "# v2 git bundle"
    git -C "$REPO" rev-list --boundary "$SHA" "^$BASE" | sed -n 's/^-\(.*\)/-\1/p'
    echo "$SHA refs/heads/ci"; echo
    printf '%s\n^%s\n' "$SHA" "$BASE" | git -C "$REPO" pack-objects --revs --thin --stdout --delta-base-offset -q
  } > "$F"
  git -C "$REPO" bundle verify -q "$F" > /dev/null 2>&1 || { echo "[ci_submit] the bundle does not verify"; exit 1; }
  BSHA=$(sha256sum "$F" | cut -d' ' -f1)
  gh release view jobs -R "$DATA" > /dev/null 2>&1 ||
    gh release create jobs -R "$DATA" --title jobs --latest=false --prerelease \
      --notes "git bundles of CPPMPAS commits under test (tools/ci/ci_submit.sh; deleted after each run)" > /dev/null
  BUNDLE=$TAG.bundle  # from here on the trap deletes it
  gh release upload jobs "$F" -R "$DATA" > /dev/null
  echo "[ci_submit] bundle $BUNDLE ($(stat -c %s "$F") bytes, sha256 ${BSHA:0:16}, $(git -C "$REPO" rev-list --count "$SHA" "^$BASE") commits on origin/main ${BASE:0:8})"
fi
gh workflow run cpu_suite.yml -R "$GH_REPO_SRC" --ref "$WREF" -f ref="$SHA" -f labels="$LABELS" \
  -f bundle="${BUNDLE:-}" -f bundle_sha256="$BSHA" -f tag="$TAG" -f max_parallel="$MAXP" > /dev/null
RUN=
for _ in $(seq 30); do  # the run of this dispatch, by its name
  sleep 5
  RUN=$(gh run list -R "$GH_REPO_SRC" --workflow cpu_suite.yml --limit 30 --json databaseId,displayTitle \
        --jq ".[] | select(.displayTitle == \"$TAG\") | .databaseId" | head -1)
  [ -n "$RUN" ] && break
done
[ -n "$RUN" ] || { echo "[ci_submit] the run $TAG did not appear"; exit 1; }
URL=https://github.com/$GH_REPO_SRC/actions/runs/$RUN
echo "[ci_submit] commit ${SHA:0:8} ($LABELS), run $RUN: $URL"
t0=$SECONDS last=
while :; do  # every minute: the run's status, and the jobs that ended
  st=$(gh run view "$RUN" -R "$GH_REPO_SRC" --json status,conclusion --jq '.status + " " + .conclusion' 2> /dev/null || echo "unknown ")
  jobs=$(gh run view "$RUN" -R "$GH_REPO_SRC" --json jobs 2> /dev/null || echo '{"jobs":[]}')
  ndone=$(jq '[.jobs[] | select(.status == "completed")] | length' <<< "$jobs")
  line="$st, $ndone jobs done"
  [ "$line" = "$last" ] || { echo "[ci_submit] $(date +%H:%M) $line"; last=$line; }
  case $st in completed*) break ;; esac
  if [ "$WAIT" = 0 ]; then  # the bundle is needed until the build job has checked the commit out
    fetched=$(jq -r '[.jobs[] | select(.name == "build") | (.status == "completed" or
                     ([.steps[]? | select(.name == "source from the bundle" and .status == "completed")] | length > 0))]
                     | any' <<< "$jobs")
    [ -n "$BUNDLE" ] && [ "$fetched" != true ] || { cleanup; echo "[ci_submit] not waiting for the suite: $URL"; exit 0; }
  fi
  [ $((SECONDS - t0)) -lt $((HOURS * 3600)) ] || { echo "[ci_submit] still running after $HOURS h: $URL"; exit 1; }
  sleep 60
done
cleanup
OUT=${OUT:-$HOME/cppwrf_scratch/ci/cppmpas/$RUN}
mkdir -p "$OUT"
echo "$URL" > "$OUT/run.txt"
ok=0
for i in 1 2 3 4 5; do  # the artifact can take a moment to be listed after the run ends
  rm -rf "$OUT/summary"
  gh run download "$RUN" -R "$GH_REPO_SRC" -n summary -D "$OUT/summary" > /dev/null 2>&1 && ok=1 && break
  sleep 20
done
[ $ok = 1 ] || { echo "[ci_submit] no summary artifact (run ${st#completed }): $URL"; exit 1; }
cat "$OUT/summary/summary.md"
echo "[ci_submit] summary in $OUT/summary; run conclusion: ${st#completed }; $URL"
passed=$(jq -r .passed "$OUT/summary/summary.json")
[ "${st#completed }" = success ] && [ "$passed" = true ]
