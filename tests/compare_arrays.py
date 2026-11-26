#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Exact comparison of named variables of two MPAS NetCDF files (e.g. a reference restart or history file and a file
written by `cppmpas init`): type, dimension names, shape and variable attributes must agree, and the raw bytes of the
values must be identical (tolerance 0). Only the first time record of a variable with a Time dimension is compared
unless --record selects another one in the reference. A run that compares nothing is a failure.

usage: compare_arrays.py REFERENCE CANDIDATE VAR[,VAR...] [--record N]
Exit status: 0 identical, 1 differences, 2 usage or input error.
"""
import argparse
import sys
from pathlib import Path

import netCDF4 as nc
import numpy as np


def ulp_max(a, b):
    ia = a.view(np.int64).astype(object)
    ib = b.view(np.int64).astype(object)
    def ordered(v):
        return v if v >= 0 else -(2 ** 63) - v - 1
    return max(abs(ordered(x) - ordered(y)) for x, y in zip(ia, ib))


def attributes(v):
    """{name: (type, raw bytes)} of a variable's attributes: an NC_INT 1 differs from an NC_CHAR "1", an NC_FLOAT
    from an NC_DOUBLE of the same printed value."""
    out = {}
    for k in v.ncattrs():
        a = v.getncattr(k)
        if isinstance(a, str):
            out[k] = ("char", a.encode())
        else:
            arr = np.asarray(a)
            out[k] = (arr.dtype.str, arr.tobytes())
    return out


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("reference", type=Path)
    p.add_argument("candidate", type=Path)
    p.add_argument("vars")
    p.add_argument("--record", type=int, default=0, help="time record of the reference (default 0)")
    a = p.parse_args(argv)
    names = [v for v in a.vars.split(",") if v]
    bad = 0
    compared = 0
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
            if attributes(x) != attributes(y):
                print(f"DIFF {n}: variable attributes")
                bad += 1
            if x.dimensions and x.dimensions[0] == "Time":
                if y.shape[0] != 1:
                    print(f"DIFF {n}: candidate has {y.shape[0]} records")
                    bad += 1
                    continue
                ra, rb = np.ascontiguousarray(x[a.record]), np.ascontiguousarray(y[0])
            else:
                ra, rb = np.ascontiguousarray(x[...]), np.ascontiguousarray(y[...])
            if ra.shape != rb.shape:
                print(f"DIFF {n}: shape {ra.shape} vs {rb.shape}")
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
