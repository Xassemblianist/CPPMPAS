#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Exact comparison of named variables of two MPAS NetCDF files (e.g. a reference restart or history file and a file
written by `cppmpas init`): type, dimension names, shape and variable attributes must agree, and the raw bytes of the
values must be identical (tolerance 0). Only the first time record of a variable with a Time dimension is compared
unless --record selects another one in the reference. A run that compares nothing is a failure, and so is a
variable with no value, a reference record that does not exist, a file shorter than its header says (truncated) and
a file compared with itself. The attributes of a classic-format file are read raw from its header (tools/
mpas_compare.py): type, every byte (NUL bytes included, which netCDF4-python drops) and order.

usage: compare_arrays.py REFERENCE CANDIDATE VAR[,VAR...] [--record N]
Exit status: 0 identical, 1 differences, 2 usage or input error.
"""
import argparse
import sys
from pathlib import Path

import netCDF4 as nc
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
import mpas_compare  # noqa: E402  (its raw header reader)


def ulp_max(a, b):
    ia = a.view(np.int64).astype(object)
    ib = b.view(np.int64).astype(object)
    def ordered(v):
        return v if v >= 0 else -(2 ** 63) - v - 1
    return max(abs(ordered(x) - ordered(y)) for x, y in zip(ia, ib))


def attributes(v):
    """[(name, type, raw bytes)] of a variable's attributes in file order: an NC_INT 1 differs from an NC_CHAR "1",
    an NC_FLOAT from an NC_DOUBLE of the same printed value, and the same attributes in another order differ."""
    out = []
    for k in v.ncattrs():
        a = v.getncattr(k)
        if isinstance(a, str):
            out.append((k, "char", a.encode()))
        else:
            arr = np.asarray(a)
            out.append((k, arr.dtype.str, arr.tobytes()))
    return out


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("reference", type=Path)
    p.add_argument("candidate", type=Path)
    p.add_argument("vars")
    p.add_argument("--record", type=int, default=0, help="time record of the reference (default 0)")
    a = p.parse_args(argv)
    names = a.vars.split(",")
    if "" in names:
        print(f"ERROR empty name in the list {a.vars!r}")
        return 2
    bad = 0
    compared = 0
    if a.reference.resolve() == a.candidate.resolve():
        print(f"ERROR {a.reference}: the reference and the candidate are the same file")
        return 2
    try:
        hr, hc = mpas_compare.classic_header(a.reference), mpas_compare.classic_header(a.candidate)
    except mpas_compare.InputError as exc:
        print(f"ERROR {exc}")
        return 2
    for path, h in ((a.reference, hr), (a.candidate, hc)):
        for problem in mpas_compare.header_problems(path, h) if h else []:
            print(f"DIFF {path.name}: {problem}")
            bad += 1
    with nc.Dataset(a.reference) as r, nc.Dataset(a.candidate) as c:
        for ds in (r, c):
            ds.set_auto_maskandscale(False)
            ds.set_auto_chartostring(False)
        for n in names:
            if n not in r.variables or n not in c.variables:
                print(f"ERROR {n}: missing in {'reference' if n not in r.variables else 'candidate'}")
                return 2
            x, y = r.variables[n], c.variables[n]
            for v in (x, y):
                v.set_auto_maskandscale(False)
                v.set_auto_chartostring(False)
            if np.dtype(x.dtype) != np.dtype(y.dtype) or x.dimensions != y.dimensions:
                print(f"DIFF {n}: type/dims {x.dtype}{x.dimensions} vs {y.dtype}{y.dimensions}")
                bad += 1
                continue
            if hr and hc:
                atts_x, atts_y = mpas_compare.raw_attrs(hr["vars"][n][1]), mpas_compare.raw_attrs(hc["vars"][n][1])
            else:
                atts_x, atts_y = attributes(x), attributes(y)
            if atts_x != atts_y:
                print(f"DIFF {n}: variable attributes")
                bad += 1
            if x.dimensions and x.dimensions[0] == "Time":
                if y.shape[0] != 1:
                    print(f"DIFF {n}: candidate has {y.shape[0]} records")
                    bad += 1
                    continue
                if not 0 <= a.record < x.shape[0]:
                    print(f"ERROR {n}: the reference has {x.shape[0]} records, record {a.record} asked")
                    return 2
                ra, rb = np.ascontiguousarray(x[a.record]), np.ascontiguousarray(y[0])
            else:
                ra, rb = np.ascontiguousarray(x[...]), np.ascontiguousarray(y[...])
            if ra.shape != rb.shape:
                print(f"DIFF {n}: shape {ra.shape} vs {rb.shape}")
                bad += 1
                continue
            if ra.size == 0:
                print(f"DIFF {n}: no values")
                bad += 1
                continue
            compared += 1
            fa = ra.reshape(-1).view(np.uint8).reshape(ra.size, -1)
            fb = rb.reshape(-1).view(np.uint8).reshape(rb.size, -1)
            where = np.flatnonzero((fa != fb).any(axis=1))
            if where.size:
                i0 = np.unravel_index(int(where[0]), ra.shape)
                extra = ""
                if ra.dtype.kind == "f" and ra.dtype.itemsize == 8:
                    extra = f" max_ulp={ulp_max(ra.reshape(-1)[where], rb.reshape(-1)[where])}"
                print(f"DIFF {n}: {where.size}/{ra.size} values differ, first {list(map(int, i0))} "
                      f"ref={ra[i0]!r} cand={rb[i0]!r}{extra}")
                bad += 1
            else:
                print(f"SAME {n} ({ra.size} values)")
    if bad:
        return 1
    if compared == 0:
        print("ERROR nothing was compared")
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
