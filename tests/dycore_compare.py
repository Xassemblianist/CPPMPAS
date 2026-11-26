#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Run tools/mpas_compare.py on a reference and a candidate directory of step files and require that every variable
is identical: the report may hold no difference other than "global attributes", and exactly N variables must have been
compared. The global attributes are then checked file by file here: every attribute the candidate writes has the
reference's value (file_id, a random identifier, excepted), and the reference's attributes the candidate lacks are
only the namelist options (config_*) and the physics attributes MPAS writes (NOT_WRITTEN), which CPPMPAS's writer
does not write yet.

usage: dycore_compare.py <mpas_compare.py> <reference dir> <candidate dir> <report.json> <N>
"""
import json
import subprocess
import sys
from pathlib import Path

import netCDF4
import numpy as np

NOT_WRITTEN = {"input_soil_data", "input_soil_temperature_lag", "num_soil_layers", "months", "noznlev", "naerlev",
               "camdim1"}


def global_attribute_problems(ref_dir, cand_dir):
    out = []
    names = sorted(p.name for p in Path(ref_dir).glob("*.nc"))
    if not names:
        return ["no reference files"]
    for n in names:
        with netCDF4.Dataset(Path(ref_dir) / n) as r, netCDF4.Dataset(Path(cand_dir) / n) as c:
            ra = {k: r.getncattr(k) for k in r.ncattrs()}
            ca = {k: c.getncattr(k) for k in c.ncattrs()}
        for k, v in ca.items():
            if k == "file_id":
                continue
            if k not in ra:
                out.append(f"{n}: attribute {k} not in the reference")
            elif type(v) is not type(ra[k]) or not np.array_equal(np.asarray(v), np.asarray(ra[k])):
                out.append(f"{n}: attribute {k} = {v!r}, reference {ra[k]!r}")
        for k in ra:
            if k not in ca and not k.startswith("config_") and k not in NOT_WRITTEN:
                out.append(f"{n}: attribute {k} missing")
    return out


def main():
    tool, ref, cand, report, n = sys.argv[1:6]
    r = subprocess.run([sys.executable, "-I", tool, ref, cand, "--json", report], capture_output=True, text=True)
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
    atts = global_attribute_problems(ref, cand)
    for x in atts:
        print("GLOBAL", x)
    return 1 if other or atts else 0


if __name__ == "__main__":
    raise SystemExit(main())
