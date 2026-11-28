#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# The CPU suite's inputs (ci/data.tsv): download them, check every SHA-256, and extract them where the build's cache
# variables point (ci/build.sh). Runs on the GitHub runner, outside the container.
#   data.sh fetch <item> <dir>            download the files of <item> into <dir> and check them
#   data.sh extract <item> <dir> <data>   check the files of <item> in <dir> again, extract them under <data>
#                                         (common: cases/, mpas_src/, mmm/, mpasbuild/tables/, rrtmg/;
#                                         a reference tag: refs/<tag>/, then its own SHA256SUMS), delete <dir>
#   data.sh bytes <item>                  the compressed size of <item>
set -euo pipefail
export LC_ALL=C
here=$(cd "$(dirname "$0")" && pwd)
table=$here/data.tsv
DATA_REPO=Xassemblianist/cppmpas-data
cmd=${1:-}; item=${2:-}
[ -n "$cmd" ] && [ -n "$item" ] || { sed -n 3,9p "$0"; exit 2; }
rows=$(awk -F'\t' -v i="$item" '$1 == i' "$table")
[ -n "$rows" ] || { echo "data.sh: no item '$item' in ci/data.tsv" >&2; exit 2; }

check() {  # check <dir>: every file of the item against its SHA-256
  awk -F'\t' -v d="$1" '{print $3 "  " d "/" $2}' <<< "$rows" | sha256sum --strict -c -
}

case $cmd in
bytes)
  awk -F'\t' '{s += $4} END {print s}' <<< "$rows" ;;
fetch)
  dir=$3; mkdir -p "$dir"
  while IFS=$'\t' read -r _ file sha bytes source; do
    case $source in
      release:*)
        gh release download "${source#release:}" -R "$DATA_REPO" -p "$file" --dir "$dir" --clobber ;;
      parts:*)
        tag=${source#parts:}
        gh release download "$tag" -R "$DATA_REPO" -p "$file.part*" -p "$tag.sha256" --dir "$dir" --clobber
        # the release's own list names the joined archive with the SHA-256 of ci/data.tsv
        grep -qx "$sha  $file" "$dir/$tag.sha256" || { echo "data.sh: $tag.sha256 does not list $sha $file" >&2; exit 1; }
        (cd "$dir" && sha256sum --strict -c --ignore-missing "$tag.sha256")
        cat "$dir/$file".part* > "$dir/$file"
        rm -f "$dir/$file".part* "$dir/$tag.sha256" ;;
      url:*)
        curl -fsSL --retry 5 -o "$dir/$file" "${source#url:}" ;;
      *) echo "data.sh: unknown source '$source'" >&2; exit 2 ;;
    esac
    [ "$(stat -c %s "$dir/$file")" = "$bytes" ] || { echo "data.sh: $file is not $bytes bytes" >&2; exit 1; }
  done <<< "$rows"
  check "$dir"
  du -sb "$dir" ;;
extract)
  dir=$3; data=$4
  check "$dir"
  if [ "$item" = common ]; then
    mkdir -p "$data"/{cases,mpas_src,mmm,mpasbuild,rrtmg}
    for c in mountain_wave supercell jw_baroclinic_wave; do tar -xzf "$dir/$c.tar.gz" -C "$data/cases"; done
    tar -xzf "$dir"/mpas-v8.4.2-f34984b.tar.gz -C "$data/mpas_src" --strip-components=1
    tar -xzf "$dir"/mmm-physics-a4baf7f.tar.gz -C "$data/mmm" --strip-components=1
    # the reference build's tables (RRTMG_*_DATA.DBL, OZONE_*.TBL, VEGPARM.TBL, ...), checked with the asset's list
    zstd -dc "$dir"/mpas-8.4.2-double-*.tar.zst | tar -xf - -C "$data/mpasbuild" --strip-components=1 build/SHA256SUMS build/tables
    (cd "$data/mpasbuild" && grep -E '  \./tables/' SHA256SUMS | sha256sum --strict -c --quiet -)
    cp "$dir"/rrtmg_lw_tables.bin "$dir"/rrtmg_sw_tables.bin "$data/rrtmg/"
  else
    mkdir -p "$data/refs/$item"
    zstd -dc "$dir/$item.tar.zst" | tar -xf - -C "$data/refs/$item"
    (cd "$data/refs/$item" && sha256sum --strict -c --quiet SHA256SUMS)
  fi
  rm -rf "$dir"
  echo "extracted $item: $(du -sh "$data" | cut -f1) in $data" ;;
*) sed -n 3,9p "$0"; exit 2 ;;
esac
