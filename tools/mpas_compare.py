#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Exact comparison of MPAS NetCDF output: field by field, time record by time record.

Two files, or two directories whose NetCDF files are paired by name, are compared at tolerance 0 on the raw bytes
of every variable (no CF masking, scaling or character conversion). For a float variable that differs, the report
gives the first differing time record (with its xtime), the number of differing values, the first differing index,
both values and the largest distance in units in the last place. A run that compares nothing is a failure, and so is
a file whose time dimension has no record or that is shorter than its header says (a truncated file). The global and
variable attributes of a classic-format file are read raw from its header: name, type, every byte of the value (NUL
bytes included, which netCDF4-python drops) and their order. A file compared with itself is an input error.

Exit status: 0 identical, 1 differences, 2 usage or input error.
"""
import argparse
import json
import sys
from pathlib import Path

import netCDF4 as nc
import numpy as np

# The one global attribute that differs between two identical runs: the random file_id MPAS seeds from time()
# (src/framework/random_id.c, mpas_stream_manager.F). Its presence, type and position are still compared, its value
# not. MPAS writes no wall-clock date; parent_id is the input file's file_id and parent_id (mpas_bootstrapping.F), the
# same for every run that reads the same input, and is compared.
RANDOM_ATTRS = frozenset({"file_id"})
NC_SUFFIXES = (".nc", ".nc4")
# NetCDF external types (netcdf.h): numpy dtype and size
NC_TYPES = {1: ("i1", 1), 2: ("S1", 1), 3: (">i2", 2), 4: (">i4", 4), 5: (">f4", 4), 6: (">f8", 8), 7: ("u1", 1),
            8: (">u2", 2), 9: (">u4", 4), 10: (">i8", 8), 11: (">u8", 8)}


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
    """[(name, signature)] in file order, through netCDF4-python (files that are not classic format)"""
    return [(name, None if name in excluded else attr_signature(obj.getncattr(name))) for name in obj.ncattrs()]


class _Reader:
    def __init__(self, f, version):
        self.f, self.big = f, version == 5

    def take(self, n):
        b = self.f.read(n)
        if len(b) != n:
            raise InputError(f"{self.f.name}: the header ends early")
        return b

    def u32(self):
        return int.from_bytes(self.take(4), "big")

    def count(self):  # NON_NEG: 8 bytes in CDF-5, 4 before
        return int.from_bytes(self.take(8 if self.big else 4), "big")

    def name(self):
        n = self.count()
        raw = self.take((n + 3) // 4 * 4)[:n]
        return raw.decode("utf-8", "replace")

    def att_list(self):
        tag, n = self.u32(), self.count()
        if tag not in (0, 0x0C) or (tag == 0 and n):
            raise InputError(f"{self.f.name}: bad attribute list in the header")
        out = []
        for _ in range(n):
            name, typ = self.name(), self.u32()
            if typ not in NC_TYPES:
                raise InputError(f"{self.f.name}: attribute {name} has the unknown type {typ}")
            k = self.count()
            size = k * NC_TYPES[typ][1]
            out.append((name, typ, self.take((size + 3) // 4 * 4)[:size]))
        return out


def classic_header(path):
    """The header of a classic-format (CDF-1, CDF-2, CDF-5) file, read raw: numrecs, dimensions, the global attributes
    [(name, type, value bytes)] in file order and every variable's (dimension ids, attributes, type, vsize, begin).
    None for a file of another format (netCDF-4)."""
    with open(path, "rb") as f:
        magic = f.read(4)
        if magic[:3] != b"CDF" or len(magic) != 4:
            return None
        version = magic[3]
        if version not in (1, 2, 5):
            raise InputError(f"{path}: CDF version {version}")
        r = _Reader(f, version)
        numrecs = r.count()
        if numrecs == (2 ** 64 - 1 if version == 5 else 2 ** 32 - 1):
            raise InputError(f"{path}: the record count is indeterminate (streaming)")
        tag, n = r.u32(), r.count()
        dims = []
        for _ in range(n):
            dims.append((r.name(), r.count()))
        gatts = r.att_list()
        tag, n = r.u32(), r.count()
        variables = {}
        for _ in range(n):
            name = r.name()
            dimids = [r.count() for _ in range(r.count())]
            vatts = r.att_list()
            typ, vsize = r.u32(), r.count()
            begin = int.from_bytes(r.take(4 if version == 1 else 8), "big")
            variables[name] = (dimids, vatts, typ, vsize, begin)
    return {"numrecs": numrecs, "dims": dims, "gatts": gatts, "vars": variables}


def header_problems(path, h):
    """Data past the end of the file (a truncated file reads as zeros through the NetCDF library)."""
    out = []
    unlimited = [i for i, (_, n) in enumerate(h["dims"]) if n == 0]
    record_vars = {k: v for k, v in h["vars"].items() if v[0] and v[0][0] in unlimited}
    def nbytes(v, skip_first):
        n = NC_TYPES[v[2]][1]
        for d in v[0][1 if skip_first else 0:]:
            n *= h["dims"][d][1]
        return n
    recsize = sum(v[3] for v in record_vars.values())
    if len(record_vars) == 1:
        recsize = nbytes(next(iter(record_vars.values())), True)
    end = 0
    for name, v in h["vars"].items():
        if name in record_vars:
            if h["numrecs"]:
                end = max(end, v[4] + (h["numrecs"] - 1) * recsize + nbytes(v, True))
        else:
            end = max(end, v[4] + nbytes(v, False))
    size = Path(path).stat().st_size
    if size < end:
        out.append(f"the file has {size} bytes, its header places data up to byte {end} (truncated)")
    return out


def raw_attrs(atts, excluded=frozenset()):
    return [(name, typ, None if name in excluded else raw) for name, typ, raw in atts]


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
    ha, hb = classic_header(ref_path), classic_header(cand_path)
    for path, h in ((ref_path, ha), (cand_path, hb)):
        for problem in header_problems(path, h) if h else []:
            out.append({"variable": None, "what": f"{Path(path).name}: {problem}"})
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
        if ha and hb:
            global_a, global_b = raw_attrs(ha["gatts"], RANDOM_ATTRS), raw_attrs(hb["gatts"], RANDOM_ATTRS)
        else:
            global_a, global_b = attrs(a, RANDOM_ATTRS), attrs(b, RANDOM_ATTRS)
        if global_a != global_b:
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
            if ha and hb:
                vatts_x, vatts_y = raw_attrs(ha["vars"][name][1]), raw_attrs(hb["vars"][name][1])
            else:
                vatts_x, vatts_y = attrs(x), attrs(y)
            if vatts_x != vatts_y:
                out.append({"variable": name, "what": "variable attributes"})
            if x.dimensions and x.dimensions[0] == "Time" and x.shape[0] == 0:
                out.append({"variable": name, "what": "no time record"})
                continue
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
    if ref.resolve() == cand.resolve():
        raise InputError(f"{ref}: the reference and the candidate are the same")
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
    select = args.vars.split(",") if args.vars is not None else None  # an empty name is not in the reference
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
