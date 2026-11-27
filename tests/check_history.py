#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""History writer check. The expected layout is derived here independently from MPAS's Registry.xml (not from
CPPMPAS's generated table): NetCDF type per Registry type, dimensions = Time first, Fortran dimensions reversed,
StrLen = 64 for text (mpas_io_streams.F:1518-1519); dimensions defined in the order MPAS defines them, field by
field in Fortran order: StrLen, inner dimensions, outer dimension, Time (mpas_io_streams.F:1518-1519, 1700-1745);
attributes, in order, = every XML attribute except those MPAS ignores, 'description' renamed 'long_name'
(src/tools/registry/gen_inc.c:43-58), with the init file's value of each where the init file has it (the fields are
read from it: put_get_field_atts, mpas_io_streams.F:4315-4360); CDF-5 (SMIOL default, smiolf.F90:426-430); global attributes of
mpas_subdriver.F:513-530 in that order, then file_id (mpas_stream_manager.F:4306-4310): the given one, or a random
10-character [0-9a-z] id when none is given. Values must equal the source init file's bytes (file order).
The Registry is preprocessed as the atmosphere core's build does: DO_PHYSICS defined, other macros undefined,
#include expanded (core_atmosphere/Makefile:9-11, core_reg).

usage: check_history.py <cppmpas> <init.nc> <Registry.xml> <workdir> <fields> [scalars] [file_id]
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


def cpp(path, defined=("DO_PHYSICS",)):
    keep, lines = [True], []
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        t = line.strip()
        if t.startswith("#"):
            words = t[1:].split()
            if words[0] in ("ifdef", "ifndef"):
                keep.append(keep[-1] and ((words[1] in defined) == (words[0] == "ifdef")))
            elif words[0] == "else":
                keep[-1] = keep[-2] and not keep[-1]
            elif words[0] == "endif":
                keep.pop()
            elif words[0] == "include" and keep[-1]:
                inc = cpp(Path(path).parent / words[1].strip('"'), defined)
                lines.append(re.sub(r"<\?xml[^>]*\?>", "", inc))
            continue
        if keep[-1]:
            lines.append(line)
    return "\n".join(lines)


def registry(path):
    root = ET.fromstring(cpp(path))
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
    file_id = sys.argv[7] if len(sys.argv) > 7 else None
    work = Path(work)
    work.mkdir(parents=True, exist_ok=True)
    out = work / (Path(init).stem + ".history.nc")
    cmd = [exe, "history", init, str(out), "--fields", fields, "--scalars", scalars]
    subprocess.run(cmd + (["--file-id", file_id] if file_id else []), check=True)
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
            # a field read from the input file: put_get_field_atts (mpas_io_streams.F:4315-4360) replaces the value
            # of each of its Registry attributes by the input file's (a text attribute the file has)
            src_atts = {a: s[n].getncattr(a) for a in s[n].ncattrs()} if n in s.variables else {}
            want_atts = []
            for k, val in el.attrib.items():
                if k in IGNORED:
                    continue
                k = "long_name" if k == "description" else k
                fv = src_atts.get(k)
                want_atts.append((k, fv.rstrip(" ") if isinstance(fv, str) else val))
            got_atts = [(a, v.getncattr(a)) for a in v.ncattrs()]
            if got_atts != want_atts:
                fails.append(f"{n}: attributes {got_atts} != {want_atts}")
            a = np.asarray(v[...])
            b = np.asarray(s[n][...])
            if a.tobytes() != b.tobytes():
                fails.append(f"{n}: values differ from the init file")
        want_order = []
        for n in names:
            typ, dims, _ = reg[n]
            fortran = (["StrLen"] if typ == "text" else []) + [d for d in dims if d != "Time"] + \
                (["Time"] if dims and dims[-1] == "Time" else [])
            want_order += [d for d in fortran if d not in want_order]
        if list(h.dimensions) != want_order:
            fails.append(f"dimension order {list(h.dimensions)}, MPAS defines {want_order}")
        for d in h.dimensions:
            want = {"Time": None, "StrLen": 64}.get(d, len(s.dimensions[d]) if d in s.dimensions else -1)
            if (want is None) != h.dimensions[d].isunlimited() or (want is not None and len(h.dimensions[d]) != want):
                fails.append(f"dimension {d}")
        want_g = [("model_name", "mpas"), ("core_name", "atmosphere"), ("version", "8.4.2"), ("source", "MPAS"),
                  ("Conventions", "MPAS"), ("git_version", "v8.4.2"),
                  ("on_a_sphere", "YES" if "YES" in s.getncattr("on_a_sphere") else "NO"),
                  ("sphere_radius", s.getncattr("sphere_radius")), ("is_periodic", "NO"), ("x_period", 0.0),
                  ("y_period", 0.0), ("history", "mpirun -n 1 ./atmosphere_model"),
                  ("parent_id", s.getncattr("file_id").rstrip(" ") + "\n" + s.getncattr("parent_id").rstrip(" ")),
                  ("mesh_spec", s.getncattr("mesh_spec").rstrip(" "))]
        got_g = [(a, h.getncattr(a)) for a in h.ncattrs()]
        got_id = got_g.pop() if got_g and got_g[-1][0] == "file_id" else ("file_id", None)
        if got_g != want_g:
            fails.append(f"global attributes {got_g} != {want_g}")
        if file_id is not None and got_id[1] != file_id:
            fails.append(f"file_id {got_id[1]!r}, given {file_id!r}")
        if file_id is None and not re.fullmatch(r"[0-9a-z]{10}", got_id[1] or ""):
            fails.append(f"file_id {got_id[1]!r} is not a random 10-character [0-9a-z] id")
    for f in fails:
        print("FAIL", f)
    if fails:
        return 1
    print(f"OK history: {len(names)} variables, layout from Registry.xml, values = init file")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
