#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""tests/compare_arrays.py on small files: identical files compare equal; a difference in a value's bytes, a
variable's type, or an attribute's type (NC_INT 1 against NC_CHAR "1", NC_FLOAT against NC_DOUBLE) or value is
reported.

usage: test_compare_arrays.py <scratch dir>
"""
import contextlib
import io
import sys
from pathlib import Path

import netCDF4 as nc
import numpy as np

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import compare_arrays  # noqa: E402


def make(path, vtype="f8", att=("units", "m"), values=(1.0, 2.0, -0.0)):
    with nc.Dataset(path, "w", format="NETCDF3_64BIT_DATA") as d:
        d.createDimension("n", 3)
        v = d.createVariable("x", vtype, ("n",))
        v.setncattr(att[0], att[1])
        v[:] = np.array(values, dtype=vtype)


def run(a, b):
    with contextlib.redirect_stdout(io.StringIO()) as out:
        rc = compare_arrays.main([str(a), str(b), "x"])
    return rc, out.getvalue().strip()


def main():
    w = Path(sys.argv[1])
    w.mkdir(parents=True, exist_ok=True)
    ref = w / "ca_ref.nc"
    make(ref, att=("scale", np.int32(1)))
    cases = [
        ("identical", dict(att=("scale", np.int32(1))), 0),
        ("char attribute for int", dict(att=("scale", "1")), 1),
        ("int64 attribute for int", dict(att=("scale", np.int64(1))), 1),
        ("attribute value", dict(att=("scale", np.int32(2))), 1),
        ("signed zero", dict(att=("scale", np.int32(1)), values=(1.0, 2.0, 0.0)), 1),
        ("variable type", dict(vtype="f4", att=("scale", np.int32(1))), 1),
    ]
    ref2 = w / "ca_ref2.nc"
    make(ref2, att=("scale", np.float64(0.5)))
    bad = 0
    for name, kw, want in cases:
        cand = w / "ca_cand.nc"
        make(cand, **kw)
        rc, out = run(ref, cand)
        ok = rc == want
        print(f"{'ok  ' if ok else 'FAIL'} {name}: rc {rc} (want {want}) {out}")
        bad += not ok
    cand = w / "ca_cand.nc"
    make(cand, att=("scale", np.float32(0.5)))
    rc, out = run(ref2, cand)
    print(f"{'ok  ' if rc == 1 else 'FAIL'} float attribute for double: rc {rc} {out}")
    bad += rc != 1
    for p in (ref, ref2, cand):
        p.unlink()
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
