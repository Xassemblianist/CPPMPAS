#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Global attributes as MPAS's bootstrapping reads them (mpas_bootstrapping.F:313-339):

* a grid file without file_id / parent_id / mesh_spec gives MPAS's three warnings;
* parent_id = trim(file_id) // NEW_LINE // trim(parent_id), assigned to domain % parent_id, a
  character(len=StrKIND = 512) variable (mpas_domain_types.inc:33): the joined value is cut at 512 characters, and the
  history file carries that value (mpas_subdriver.F:527).

usage: check_global_atts.py <cppmpas> <grid without ids> <init with a long parent_id> <workdir>
"""
import subprocess
import sys
from pathlib import Path

import netCDF4 as nc


def main():
    exe, grid, init, work = sys.argv[1:5]
    fails = []
    r = subprocess.run([exe, "info", grid], capture_output=True, text=True)
    if r.returncode != 0:
        fails.append(f"info {grid} failed: {r.stderr.strip()}")
    for a in ("file_id", "parent_id", "mesh_spec"):
        if f"MPAS warning: Attribute {a} not found" not in r.stderr:
            fails.append(f"no warning for the missing {a}")
    with nc.Dataset(grid) as g:
        if any(a in g.ncattrs() for a in ("file_id", "parent_id", "mesh_spec")):
            fails.append(f"{grid} has file_id / parent_id / mesh_spec: not a test of the warnings")

    out = Path(work) / "long_parent.history.nc"
    subprocess.run([exe, "history", init, str(out), "--fields", "latCell", "--file-id", "0123456789"], check=True)
    with nc.Dataset(init) as s, nc.Dataset(out) as h:
        joined = s.getncattr("file_id").rstrip(" ") + "\n" + s.getncattr("parent_id").rstrip(" ")
        if len(joined) <= 512:
            fails.append(f"{init}: joined parent_id has {len(joined)} characters, not a test of the cut")
        want = joined[:512].rstrip(" ")
        got = h.getncattr("parent_id")
        if got != want:
            fails.append(f"parent_id has {len(got)} characters, MPAS writes {len(want)}")
    for f in fails:
        print("FAIL", f)
    if fails:
        return 1
    print("OK global attributes: warnings for missing ids, parent_id cut at StrKIND = 512")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
