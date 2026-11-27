#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""The references the tests mix are one recipe: tests run with a reference's steps/namelist.atmosphere and compare
with the files of its forecast/ run, and the SC long tests take the 30-min output.nc of mpas842-dbl-sc-002 with the
init file and namelists of mpas842-dbl-sc-001. This holds only if, in every reference, the two namelists differ in
config_run_duration alone (the steps run is the forecast run cut short), and if sc-001 and sc-002 have the same init
file and namelists, byte for byte.

usage: check_ref_recipes.py <CPPMPAS_REFS_DIR>
"""
import sys
from pathlib import Path

RUNS = ["mpas842-dbl-mw-001/MW", "mpas842-dbl-sc-001/SC", "mpas842-dbl-jw-001/JW", "mpas842-dbl-rg-004/RG",
        "mpas842-dbl-rg24-001/RG24"]
SAME = [("mpas842-dbl-sc-001/SC", "mpas842-dbl-sc-002/SC",
         ["init/supercell_init.nc", "steps/namelist.atmosphere", "forecast/namelist.atmosphere"])]
ALLOWED = {"config_run_duration"}


def options(path):
    """{option: line} of a namelist, with the group of each option"""
    out, group = {}, None
    for line in path.read_text().splitlines():
        s = line.strip()
        if s.startswith("&"):
            group = s[1:]
        elif "=" in s:
            out[(group, s.split("=", 1)[0].strip())] = s
    return out


def main():
    refs = Path(sys.argv[1])
    bad = 0
    for r in RUNS:
        a, b = refs / r / "steps" / "namelist.atmosphere", refs / r / "forecast" / "namelist.atmosphere"
        if not a.is_file() or not b.is_file():
            print(f"FAIL {r}: {a} or {b} is missing")
            bad += 1
            continue
        oa, ob = options(a), options(b)
        diff = sorted(k for k in oa.keys() | ob.keys() if oa.get(k) != ob.get(k))
        other = [k for k in diff if k[1] not in ALLOWED]
        if other or not oa:
            print(f"FAIL {r}: steps and forecast namelists differ in {other or 'everything (no option read)'}")
            bad += 1
        else:
            print(f"ok   {r}: steps and forecast namelists differ only in {[k[1] for k in diff]}")
    for x, y, files in SAME:
        for f in files:
            p, q = refs / x / f, refs / y / f
            if not p.is_file() or not q.is_file() or p.read_bytes() != q.read_bytes():
                print(f"FAIL {f}: {x} and {y} differ (or one is missing)")
                bad += 1
            else:
                print(f"ok   {f}: {x} = {y}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
