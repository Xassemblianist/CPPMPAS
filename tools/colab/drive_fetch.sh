#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Fetch reference releases from the Google Drive mirror of Xassemblianist/cppmpas-data (tools/colab/drive_mirror.py,
# tools/colab/drive_manifest.tsv) into a refs directory laid out as CPPMPAS_REFS_DIR expects (tests/CMakeLists.txt):
#
#   drive_fetch.sh [--src SRC] [--manifest FILE] [--github] REFS TAG ...
#
#   SRC         where the mirror is: rclone:REMOTE:PATH (default rclone:cppwrf-drive:cppmpas-data; set RCLONE_CONFIG
#               to a file holding only that remote's section) or a directory, e.g. a Colab Drive mount
#               (/content/drive/MyDrive/cppmpas-data)
#   --manifest  default: drive_manifest.tsv next to this script
#   --github    an asset missing on the mirror, or whose copy fails its check, comes from the GitHub release instead
#               (the same check applies)
#   TAG         a reference tag (mpas842-dbl-*): its .tar.zst.partNN assets are checked against the manifest's
#               SHA-256 and the release's <tag>.sha256, joined (the joined archive checked against <tag>.sha256),
#               extracted into REFS/<tag>/ and the extracted files checked against the archive's SHA256SUMS, as
#               ci/data.sh extract does; other assets of the release (rgcp-001's Thompson tables) are copied to
#               REFS/<tag>/ and checked against their .sha256 list. A release without archive parts (inputs,
#               mpas-build) has its assets copied, checked, to REFS/<tag>/.
# A tag whose REFS/<tag>/.drive_fetch_ok names the manifest rows it was made from is left as it is. Needs
# sha256sum, zstd, tar and, for an rclone source, rclone. Exit 0 only when every tag is in place and checked.
set -euo pipefail
export LC_ALL=C
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
SRC=rclone:cppwrf-drive:cppmpas-data MAN=$here/drive_manifest.tsv GH=
while [ $# -gt 0 ]; do
  case $1 in
    --src) SRC=$2; shift 2 ;; --manifest) MAN=$2; shift 2 ;; --github) GH=1; shift ;;
    -*) sed -n 3,22p "$0"; exit 2 ;; *) break ;;
  esac
done
[ $# -ge 2 ] || { sed -n 3,22p "$0"; exit 2; }
REFS=$1; shift
mkdir -p "$REFS"
REFS=$(cd "$REFS" && pwd)
URL=https://github.com/Xassemblianist/cppmpas-data/releases/download

get() {  # get TAG ASSET DEST: the mirror's copy (or the release's), then its size and SHA-256 against the manifest
  local t=$1 a=$2 dest=$3 row bytes sha path
  row=$(awk -F'\t' -v t="$t" -v a="$a" '$1 == t && $2 == a' "$MAN")
  bytes=$(cut -f3 <<< "$row") sha=$(cut -f4 <<< "$row") path=$(cut -f6 <<< "$row")
  rm -f "$dest"
  case $SRC in
    rclone:*) rclone copyto --retries 3 --low-level-retries 20 --multi-thread-streams 8 --multi-thread-cutoff 64M "${SRC#rclone:}/${path#*/}" "$dest" || true ;;
    *) cp "$SRC/${path#*/}" "$dest" || true ;;
  esac
  if ! { [ -f "$dest" ] && [ "$(stat -c %s "$dest")" = "$bytes" ] && echo "$sha  $dest" | sha256sum --strict --quiet -c -; }; then
    [ -n "$GH" ] || { echo "[drive_fetch] $t/$a: missing on the mirror or not $bytes bytes / sha256 $sha" >&2; return 1; }
    echo "[drive_fetch] $t/$a: not usable on the mirror, from the GitHub release"
    curl -fsSL --retry 5 -o "$dest" "$URL/$t/$a"
    [ "$(stat -c %s "$dest")" = "$bytes" ] && echo "$sha  $dest" | sha256sum --strict --quiet -c - ||
      { echo "[drive_fetch] $t/$a: the release's copy fails the check too" >&2; return 1; }
  fi
}

rc=0
for t in "$@"; do
  rows=$(awk -F'\t' -v t="$t" '!/^#/ && $1 == t' "$MAN")
  [ -n "$rows" ] || { echo "[drive_fetch] $t is not in $MAN" >&2; rc=1; continue; }
  stamp=$(sha256sum <<< "$rows" | cut -d' ' -f1)
  if [ "$(cat "$REFS/$t/.drive_fetch_ok" 2>/dev/null)" = "$stamp" ]; then echo "[drive_fetch] $t: already in place"; continue; fi
  # (errexit holds in the subshell only outside a condition: its status is read afterwards)
  set +e
  (
    set -euo pipefail
    t0=$(date +%s)
    dl=$REFS/.drive_fetch/$t
    rm -rf "$REFS/$t" "$dl"; mkdir -p "$dl" "$REFS/$t"
    mapfile -t assets < <(cut -f2 <<< "$rows")
    parts=() lists=() other=()
    for a in "${assets[@]}"; do
      case $a in
        "$t".tar.zst.part[0-9][0-9]) parts+=("$a") ;;
        *.sha256) lists+=("$a") ;;
        *) other+=("$a") ;;
      esac
    done
    for a in "${lists[@]}" "${other[@]}"; do get "$t" "$a" "$dl/$a"; done
    for a in "${other[@]}"; do  # a loose asset named in a .sha256 list of the release is checked by it too
      for l in "${lists[@]}"; do
        d=$(awk -v a="$a" '{n = $2; sub(/^\*/, "", n)} n == a {print $1}' "$dl/$l")
        [ -z "$d" ] || echo "$d  $dl/$a" | sha256sum --strict --quiet -c -
      done
      mv "$dl/$a" "$REFS/$t/$a"
    done
    if [ ${#parts[@]} -gt 0 ]; then
      mapfile -t parts < <(printf '%s\n' "${parts[@]}" | sort)
      [ -f "$dl/$t.sha256" ] || { echo "[drive_fetch] $t: no $t.sha256 in the manifest" >&2; exit 1; }
      for a in "${parts[@]}"; do
        get "$t" "$a" "$dl/$a"
        (cd "$dl" && grep -E "  \*?$a\$" "$t.sha256" | sha256sum --strict --quiet -c -)
      done
      want=$(awk -v n="$t.tar.zst" '{f = $2; sub(/^\*/, "", f)} f == n {print $1}' "$dl/$t.sha256")
      [ -n "$want" ] || { echo "[drive_fetch] $t.sha256 does not list $t.tar.zst" >&2; exit 1; }
      got=$(cd "$dl" && cat "${parts[@]}" | sha256sum | cut -d' ' -f1)
      [ "$got" = "$want" ] || { echo "[drive_fetch] $t: the joined archive is sha256 $got, $t.sha256 says $want" >&2; exit 1; }
      (cd "$dl" && cat "${parts[@]}") | zstd -dc | tar -xf - -C "$REFS/$t"
      # the archive's own list, as ci/data.sh extract checks it
      (cd "$REFS/$t" && sha256sum --strict --quiet -c SHA256SUMS)
      for l in "${lists[@]}"; do [ "$l" = "$t.sha256" ] || cp "$dl/$l" "$REFS/$t/"; done
    else
      for l in "${lists[@]}"; do cp "$dl/$l" "$REFS/$t/"; done
    fi
    rm -rf "$dl"
    echo "$stamp" > "$REFS/$t/.drive_fetch_ok"
    echo "[drive_fetch] $t: $(du -sh "$REFS/$t" | cut -f1) in $REFS/$t, checked ($(( $(date +%s) - t0 )) s)"
  )
  st=$?
  set -e
  [ $st = 0 ] || { echo "[drive_fetch] FAILED $t (REFS/$t removed)" >&2; rm -rf "$REFS/.drive_fetch/$t" "${REFS:?}/$t"; rc=1; }
done
rmdir "$REFS/.drive_fetch" 2>/dev/null || true
exit $rc
