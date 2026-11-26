#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Run tools/mpas_compare.py on a reference and a candidate directory of step files and require that every variable
is identical (only the step.*.nc files of the reference directory: a reference's steps/ may also hold initdump.nc):
the report may hold no difference other than "global attributes" (the namelist options MPAS writes as
global attributes, which CPPMPAS's writer does not write yet), and exactly N variables must have been compared.

usage: dycore_compare.py <mpas_compare.py> <reference dir> <candidate dir> <report.json> <N>
"""
import json
import subprocess
import sys
import tempfile
from pathlib import Path


def main():
    tool, ref, cand, report, n = sys.argv[1:6]
    steps = sorted(Path(ref).glob("step.*.nc"))
    if not steps:
        print(f"{ref}: no step.*.nc files")
        return 1
    with tempfile.TemporaryDirectory(dir=Path(report).parent) as tmp:
        for f in steps:
            (Path(tmp) / f.name).symlink_to(f.resolve())
        r = subprocess.run([sys.executable, "-I", tool, tmp, cand, "--json", report], capture_output=True, text=True)
    print(r.stdout[-4000:], r.stderr[-2000:])
    if r.returncode not in (0, 1):
        print(f"mpas_compare.py failed with status {r.returncode}")
        return 1
    with open(report) as f:
        d = json.load(f)
    other = [x for x in d["differences"] if x.get("what") != "global attributes"]
    for x in other:
        print("DIFF", x)
    if d["compared_variables"] != int(n):
        print(f"compared {d['compared_variables']} variables, expected {n}")
        return 1
    return 1 if other else 0


if __name__ == "__main__":
    raise SystemExit(main())
