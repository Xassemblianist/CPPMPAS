#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Compare the output streams of a `cppmpas run --streams` with a reference forecast directory at tolerance 0.

Every NetCDF file of the reference directory must be in the candidate directory and the candidate may hold no other
NetCDF file; tools/mpas_compare.py must report no difference at all (values, types, dimensions, variable attributes
and the global attributes, the random file_id / parent_id excepted), over exactly the given numbers of files and
variables. A restart_timestamp file of the reference must be identical, byte for byte, in the candidate, and a
candidate restart_timestamp the reference does not have is a difference. A report file left by an earlier run is
deleted first.

usage: streams_compare.py <mpas_compare.py> <reference dir> <candidate dir> <report.json> <files> <variables>
"""
import json
import re
import subprocess
import sys
from pathlib import Path


def main():
    tool, ref, cand, report, nfiles, nvars = sys.argv[1:7]
    Path(report).unlink(missing_ok=True)
    r = subprocess.run([sys.executable, "-I", tool, ref, cand, "--json", report], capture_output=True, text=True)
    print(r.stdout[-6000:], r.stderr[-2000:])
    bad = r.returncode != 0
    m = re.search(r"SAME (\d+) file\(s\), (\d+) variable\(s\)", r.stdout)
    if not m or (m.group(1), m.group(2)) != (nfiles, nvars):
        print(f"expected SAME {nfiles} file(s), {nvars} variable(s)")
        bad = True
    if not Path(report).is_file():
        print("mpas_compare.py wrote no report")
        return 1
    with open(report) as f:
        d = json.load(f)
    if d["differences"]:
        bad = True
    ts, c = Path(ref) / "restart_timestamp", Path(cand) / "restart_timestamp"
    if ts.exists() != c.exists() or (ts.exists() and c.read_bytes() != ts.read_bytes()):
        print(f"restart_timestamp differs or is in one directory only (reference {ts.exists()}, "
              f"candidate {c.exists()})")
        bad = True
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
