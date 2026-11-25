#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Exact comparison of MPAS NetCDF output: field by field, time record by time record.

Two files, or two directories whose NetCDF files are paired by name, are compared at tolerance 0 on the raw bytes
of every variable (no CF masking, scaling or character conversion). For a float variable that differs, the report
gives the first differing time record (with its xtime), the number of differing values, the first differing index,
both values and the largest distance in units in the last place. A run that compares nothing is a failure.

Exit status: 0 identical, 1 differences, 2 usage or input error.
"""
import argparse
import json
import sys
from pathlib import Path

import netCDF4 as nc
import numpy as np

# Global attributes that differ between two identical runs: wall-clock dates, and the random file_id MPAS seeds from
# time() (src/framework/random_id.c, mpas_stream_manager.F) with the parent_id that inherits it (mpas_bootstrapping.F).
WALL_CLOCK_ATTRS = frozenset({"creation_date", "creation_time", "date_created", "date_modified", "file_id", "parent_id"})
NC_SUFFIXES = (".nc", ".nc4")


class InputError(Exception):
    pass


def ordered_int(values):
    """Map floats to integers whose order and spacing follow the float line (1 apart per ulp)."""
    width = values.dtype.itemsize
    signed = values.view(np.int64 if width == 8 else np.int32).astype(np.int64)
    bias = np.int64(-(2 ** 63)) if width == 8 else np.int64(-(2 ** 31))
    # Negative floats: reflect so that -0.0 and +0.0 are adjacent (distance 1) and order is monotone.
    return np.where(signed < 0, bias - signed - 1, signed)


def ulp_distance(a, b):
    ia, ib = ordered_int(a), ordered_int(b)
    if a.dtype.itemsize == 8:
        # Avoid int64 overflow for opposite signs by working in Python ints only for the maximum.
        diff = np.abs(ia.astype(object) - ib.astype(object))
        return int(diff.max()) if diff.size else 0
    return int(np.abs(ia - ib).max()) if ia.size else 0


def attr_signature(value):
    array = np.asarray(value)
    if array.dtype.hasobject:
        return ("object", str(value))
    return (array.dtype.str, array.shape, array.tobytes())


def attrs(obj, excluded=frozenset()):
    return {name: attr_signature(obj.getncattr(name)) for name in obj.ncattrs() if name not in excluded}


def xtime_strings(ds):
    if "xtime" not in ds.variables:
        return None
    var = ds.variables["xtime"]
    var.set_auto_chartostring(False)
    raw = np.asarray(var[:])
    return [b"".join(row.tolist()).decode("ascii", "replace").strip("\x00 ") for row in raw]


def record_slices(var):
    """Yield (record index or None, index tuple) covering the variable one time record at a time."""
    dims = var.dimensions
    if dims and dims[0] == "Time":
        for t in range(var.shape[0]):
            yield t, (t,)
    else:
        yield None, ()


def describe_difference(name, rec, label, a, b):
    diff = a.reshape(-1).view(np.uint8).reshape(a.size, -1) != b.reshape(-1).view(np.uint8).reshape(b.size, -1)
    where = np.flatnonzero(diff.any(axis=1))
    first = np.unravel_index(int(where[0]), a.shape) if a.ndim else ()
    info = {
        "variable": name,
        "record": rec,
        "xtime": label,
        "count": int(where.size),
        "size": int(a.size),
        "first_index": [int(i) for i in first],
        "reference": repr(a[first].item() if a.ndim else a.item()),
        "candidate": repr(b[first].item() if b.ndim else b.item()),
    }
    if a.dtype.kind == "f":
        info["max_ulp"] = ulp_distance(a.reshape(-1)[where], b.reshape(-1)[where])
    return info


def compare_files(ref_path, cand_path, select=None, first_only=True):
    """Return (differences, compared_count). Each difference is a dict with at least 'variable' and 'what'."""
    out = []
    compared = 0
    with nc.Dataset(ref_path) as a, nc.Dataset(cand_path) as b:
        for ds in (a, b):
            ds.set_auto_maskandscale(False)
            ds.set_auto_chartostring(False)
            if ds.groups:
                raise InputError(f"{ds.filepath()}: NetCDF groups are not MPAS output")
        if a.data_model != b.data_model:
            out.append({"variable": None, "what": f"format {a.data_model} vs {b.data_model}"})
        da = {k: (len(v), v.isunlimited()) for k, v in a.dimensions.items()}
        db = {k: (len(v), v.isunlimited()) for k, v in b.dimensions.items()}
        if da != db:
            out.append({"variable": None, "what": f"dimensions {da} vs {db}"})
        if attrs(a, WALL_CLOCK_ATTRS) != attrs(b, WALL_CLOCK_ATTRS):
            out.append({"variable": None, "what": "global attributes"})
        labels = xtime_strings(a)

        names_a, names_b = set(a.variables), set(b.variables)
        if select is not None:
            missing = [n for n in select if n not in names_a]
            if missing:
                raise InputError(f"{ref_path}: selected variables not in reference: {', '.join(missing)}")
            names_a &= set(select)
            names_b &= set(select)
        for name in sorted(names_a - names_b):
            out.append({"variable": name, "what": "missing in candidate"})
        for name in sorted(names_b - names_a):
            out.append({"variable": name, "what": "only in candidate"})
        for name in sorted(names_a & names_b):
            x, y = a.variables[name], b.variables[name]
            for v in (x, y):
                v.set_auto_maskandscale(False)
                v.set_auto_chartostring(False)
            if x.datatype is str or y.datatype is str:
                raise InputError(f"{name}: variable-length strings are not MPAS output")
            if (np.dtype(x.dtype), x.shape, x.dimensions) != (np.dtype(y.dtype), y.shape, y.dimensions):
                out.append({"variable": name, "what": f"type/shape {x.dtype}{x.shape}{x.dimensions} vs "
                                                      f"{y.dtype}{y.shape}{y.dimensions}"})
                continue
            if attrs(x) != attrs(y):
                out.append({"variable": name, "what": "variable attributes"})
            compared += 1
            for rec, index in record_slices(x):
                ra = np.ascontiguousarray(x[index] if index else x[...])
                rb = np.ascontiguousarray(y[index] if index else y[...])
                if ra.tobytes() != rb.tobytes():
                    label = labels[rec] if (labels and rec is not None and rec < len(labels)) else None
                    info = describe_difference(name, rec, label, ra, rb)
                    info["what"] = "values"
                    out.append(info)
                    if first_only:
                        break
    return out, compared


def pair_inputs(ref, cand):
    if ref.is_file() and cand.is_file():
        return [(ref, cand)], []
    if ref.is_dir() and cand.is_dir():
        rf = {p.name: p for p in ref.iterdir() if p.suffix in NC_SUFFIXES}
        cf = {p.name: p for p in cand.iterdir() if p.suffix in NC_SUFFIXES}
        problems = [{"file": n, "variable": None, "what": "file missing in candidate"} for n in sorted(rf.keys() - cf.keys())]
        problems += [{"file": n, "variable": None, "what": "file only in candidate"} for n in sorted(cf.keys() - rf.keys())]
        return [(rf[n], cf[n]) for n in sorted(rf.keys() & cf.keys())], problems
    raise InputError("give two files or two directories")


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("reference", type=Path)
    p.add_argument("candidate", type=Path)
    p.add_argument("--vars", help="comma-separated variables to compare (default: all)")
    p.add_argument("--all-records", action="store_true", help="report every differing record, not the first")
    p.add_argument("--json", type=Path, help="write the full report as JSON")
    args = p.parse_args(argv)
    select = [v for v in args.vars.split(",") if v] if args.vars else None
    try:
        pairs, report = pair_inputs(args.reference, args.candidate)
        compared = 0
        for ra, ca in pairs:
            diffs, n = compare_files(ra, ca, select, first_only=not args.all_records)
            compared += n
            for d in diffs:
                d["file"] = ra.name
            report += diffs
    except (InputError, OSError) as exc:
        print(f"ERROR {exc}", file=sys.stderr)
        return 2
    if compared == 0:
        report.append({"file": None, "variable": None, "what": "nothing was compared"})
    for d in report:
        extra = ""
        if d.get("what") == "values":
            extra = (f" record={d['record']} xtime={d['xtime']} count={d['count']}/{d['size']}"
                     f" first={d['first_index']} ref={d['reference']} cand={d['candidate']}")
            if "max_ulp" in d:
                extra += f" max_ulp={d['max_ulp']}"
        print(f"DIFF {d.get('file')} {d.get('variable')} {d['what']}{extra}")
    if args.json:
        args.json.write_text(json.dumps({"compared_variables": compared, "differences": report}, indent=1) + "\n")
    if report:
        return 1
    print(f"SAME {len(pairs)} file(s), {compared} variable(s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
