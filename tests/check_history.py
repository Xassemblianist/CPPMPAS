#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""History writer check. The expected layout is derived here independently from MPAS's Registry.xml (not from
CPPMPAS's generated table): NetCDF type per Registry type, dimensions = Time first, Fortran dimensions reversed,
StrLen = 64 for text (mpas_io_streams.F:1518-1519); attributes = every XML attribute except those MPAS ignores,
'description' renamed 'long_name' (src/tools/registry/gen_inc.c:43-58); CDF-5 (SMIOL default, smiolf.F90:426-430);
global attributes of mpas_subdriver.F:513-530. Values must equal the source init file's bytes (file order).

usage: check_history.py <cppmpas> <init.nc> <Registry.xml> <workdir> <fields> [scalars]
"""
import re
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import netCDF4 as nc
import numpy as np

IGNORED = {"name", "type", "dimensions", "persistence", "packages", "time_levs", "name_in_code", "array_group",
           "default_value"}


def registry(path):
    text = open(path, encoding="utf-8").read()
    text = "\n".join(l for l in text.splitlines() if not re.match(r"\s*#", l))  # all #ifdef bodies kept
    root = ET.fromstring(text)
    out = {}
    for vs in root.iter("var_struct"):
        for c in vs:
            if c.tag == "var":
                out[c.get("name")] = (c.get("type"), c.get("dimensions", "").split(), c)
            elif c.tag == "var_array":
                for k in c.findall("var"):
                    out[k.get("name")] = (c.get("type"), c.get("dimensions", "").split(), k)
    return out


def main():
    exe, init, reg_path, work, fields = sys.argv[1:6]
    scalars = sys.argv[6] if len(sys.argv) > 6 else "qv"
    work = Path(work)
    work.mkdir(parents=True, exist_ok=True)
    out = work / (Path(init).stem + ".history.nc")
    subprocess.run([exe, "history", init, str(out), "--fields", fields, "--scalars", scalars,
                    "--file-id", "0123456789"], check=True)
    reg = registry(reg_path)
    names = fields.split(",")
    fails = []
    with nc.Dataset(out) as h, nc.Dataset(init) as s:
        for ds in (h, s):
            ds.set_auto_mask(False)
            ds.set_auto_chartostring(False)
        if h.data_model != "NETCDF3_64BIT_DATA":
            fails.append(f"format {h.data_model}")
        if sorted(h.variables) != sorted(names):
            fails.append(f"variables {sorted(set(h.variables) ^ set(names))}")
        for n in names:
            typ, dims, el = reg[n]
            want_dims = (["Time"] if dims and dims[-1] == "Time" else []) + \
                [d for d in reversed(dims) if d != "Time"] + (["StrLen"] if typ == "text" else [])
            want_type = {"real": "f8", "integer": "i4", "text": "S1"}[typ]
            v = h[n]
            if list(v.dimensions) != want_dims or v.dtype.str[1:] != want_type:
                fails.append(f"{n}: {v.dtype}{v.dimensions}, expected {want_type}{tuple(want_dims)}")
            want_atts = {("long_name" if k == "description" else k): val
                         for k, val in el.attrib.items() if k not in IGNORED}
            got_atts = {a: v.getncattr(a) for a in v.ncattrs()}
            if got_atts != want_atts:
                fails.append(f"{n}: attributes {got_atts} != {want_atts}")
            a = np.asarray(v[...])
            b = np.asarray(s[n][...])
            if a.tobytes() != b.tobytes():
                fails.append(f"{n}: values differ from the init file")
        for d in h.dimensions:
            want = {"Time": None, "StrLen": 64}.get(d, len(s.dimensions[d]) if d in s.dimensions else -1)
            if (want is None) != h.dimensions[d].isunlimited() or (want is not None and len(h.dimensions[d]) != want):
                fails.append(f"dimension {d}")
        want_g = {"model_name": "mpas", "core_name": "atmosphere", "version": "8.4.2", "source": "MPAS",
                  "Conventions": "MPAS", "git_version": "v8.4.2", "on_a_sphere": s.getncattr("on_a_sphere"),
                  "sphere_radius": s.getncattr("sphere_radius"), "is_periodic": "NO", "x_period": 0.0,
                  "y_period": 0.0, "history": "mpirun -n 1 ./atmosphere_model",
                  "parent_id": s.getncattr("file_id") + "\n" + s.getncattr("parent_id"),
                  "mesh_spec": s.getncattr("mesh_spec"), "file_id": "0123456789"}
        got_g = {a: h.getncattr(a) for a in h.ncattrs()}
        if got_g != want_g:
            fails.append(f"global attributes {got_g} != {want_g}")
    for f in fails:
        print("FAIL", f)
    if fails:
        return 1
    print(f"OK history: {len(names)} variables, layout from Registry.xml, values = init file")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
