#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Run tools/mpas_compare.py on a reference and a candidate directory of step files and require that every variable
is identical: the report may hold no difference other than "global attributes", and exactly N variables must have been
compared. The global attributes are then checked file by file here: the candidate has exactly the reference's
attributes, in the reference's order, each with the reference's value, the namelist options (config_*) and the
physics options MPAS writes (input_soil_data, ..., camdim1) included; only file_id's value (a random identifier) is
not compared. Only the reference's step.*.nc files are compared (a reference steps/ directory may also hold other
files, such as initdump.nc). Attribute values are compared as the raw header bytes with their NetCDF type (so -0.0
differs from +0.0, NC_FLOAT from NC_DOUBLE, and NUL bytes count). A report file left by an earlier run is deleted
first.

usage: dycore_compare.py <mpas_compare.py> <reference dir> <candidate dir> <report.json> <N>
"""
import json
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
import mpas_compare  # noqa: E402  (its raw header reader)


def global_attribute_problems(ref_dir, cand_dir):
    out = []
    names = sorted(p.name for p in Path(ref_dir).glob("step.*.nc"))
    if not names:
        return ["no reference files"]
    for n in names:
        if not (Path(cand_dir) / n).is_file():
            out.append(f"{n}: not in the candidate directory")
            continue
        hr = mpas_compare.classic_header(Path(ref_dir) / n)
        hc = mpas_compare.classic_header(Path(cand_dir) / n)
        if hr is None or hc is None:
            out.append(f"{n}: not a classic-format NetCDF file")
            continue
        ra = {k: (t, raw) for k, t, raw in hr["gatts"]}
        ca = {k: (t, raw) for k, t, raw in hc["gatts"]}
        if "file_id" not in ca:
            out.append(f"{n}: no file_id")
        for k, v in ca.items():
            if k == "file_id":
                continue
            if k not in ra:
                out.append(f"{n}: attribute {k} not in the reference")
            elif v != ra[k]:
                out.append(f"{n}: attribute {k} = {v!r}, reference {ra[k]!r}")
        for k in ra:
            if k not in ca:
                out.append(f"{n}: attribute {k} missing")
        order_c = [k for k, _, _ in hc["gatts"]]
        order_r = [k for k, _, _ in hr["gatts"]]
        if len(set(order_c)) != len(order_c):
            out.append(f"{n}: an attribute name is written twice")
        if order_c != order_r:
            out.append(f"{n}: attribute order {order_c}, reference {order_r}")
    return out


def main():
    tool, ref, cand, report, n = sys.argv[1:6]
    Path(report).unlink(missing_ok=True)
    with tempfile.TemporaryDirectory(dir=Path(report).parent) as steps:  # the reference's step files only
        for f in sorted(Path(ref).glob("step.*.nc")):
            (Path(steps) / f.name).symlink_to(f.resolve())
        r = subprocess.run([sys.executable, "-I", tool, steps, cand, "--json", report], capture_output=True, text=True)
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
