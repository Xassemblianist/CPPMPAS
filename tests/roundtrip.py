#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Round trip: cppmpas reads an MPAS file and writes it back; the result must be identical (tools/mpas_compare.py,
tolerance 0: format, dimensions, global and variable attributes, value bytes) to the source restricted to the
variables CPPMPAS handles, with one MPAS rule applied to the source: connectivity slots past nEdgesOnCell /
nEdgesOnEdge are written as 0 (prewrite_reindex, mpas_stream_manager.F:5053, 5074, 5095, 5152). file_id, whose
value mpas_compare.py skips, is compared here. The text attributes of the reference file are copied as the source's
raw header bytes through the NetCDF C library (netCDF4-python writes an empty text attribute as one NUL byte and drops
NUL bytes), so the reference holds exactly the source's attributes.

usage: roundtrip.py <cppmpas> <mpas_compare.py> <src.nc> <workdir> grid|init [scalars]
"""
import ctypes
import ctypes.util
import subprocess
import sys
from pathlib import Path

import netCDF4 as nc
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
import mpas_compare  # noqa: E402  (its raw header reader)

_lib = ctypes.CDLL(ctypes.util.find_library("netcdf"))
NC_CHAR, NC_GLOBAL, NC_WRITE = 2, -1, 1


def copy_text_attributes(src, dst, names):
    """Write every text attribute of src (global, and of the variables in names) into dst with its raw bytes."""
    h = mpas_compare.classic_header(src)
    ncid = ctypes.c_int()
    if _lib.nc_open(str(dst).encode(), NC_WRITE, ctypes.byref(ncid)) or _lib.nc_redef(ncid):
        raise SystemExit(f"cannot open {dst} for writing")
    targets = [(NC_GLOBAL, h["gatts"])]
    for n in names:
        varid = ctypes.c_int()
        if _lib.nc_inq_varid(ncid, n.encode(), ctypes.byref(varid)):
            raise SystemExit(f"{dst}: no variable {n}")
        targets.append((varid.value, h["vars"][n][1]))
    for varid, atts in targets:
        for name, typ, raw in atts:
            if typ == NC_CHAR and _lib.nc_put_att_text(ncid, varid, name.encode(), ctypes.c_size_t(len(raw)), raw):
                raise SystemExit(f"cannot write {name} to {dst}")
    if _lib.nc_close(ncid):
        raise SystemExit(f"cannot close {dst}")

MESH = ["latCell", "lonCell", "xCell", "yCell", "zCell", "indexToCellID", "latEdge", "lonEdge", "xEdge", "yEdge",
        "zEdge", "indexToEdgeID", "latVertex", "lonVertex", "xVertex", "yVertex", "zVertex", "indexToVertexID",
        "cellsOnEdge", "nEdgesOnCell", "nEdgesOnEdge", "edgesOnCell", "edgesOnEdge", "weightsOnEdge", "dvEdge",
        "dcEdge", "angleEdge", "areaCell", "areaTriangle", "cellsOnCell", "verticesOnCell", "verticesOnEdge",
        "edgesOnVertex", "cellsOnVertex", "kiteAreasOnVertex", "meshDensity"]
GRID_OPTIONAL = ["fEdge", "fVertex"]  # absent from the x1.40962 grid (init_atmosphere computes them)
INIT = ["fEdge", "fVertex", "edgeNormalVectors", "localVerticalUnitVectors", "cellTangentPlane", "nominalMinDc", "bdyMaskCell",
        "bdyMaskEdge", "bdyMaskVertex", "cf1", "cf2", "cf3", "zgrid", "rdzw", "dzu", "rdzu", "fzm", "fzp", "zxu",
        "zz", "zb", "zb3", "dss", "deriv_two", "coeffs_reconstruct", "initial_time", "xtime", "u_init", "v_init",
        "t_init", "u", "w", "rho", "theta", "rho_base", "theta_base", "surface_pressure", "Time"]
UNUSED_SLOTS = {"cellsOnCell": "nEdgesOnCell", "edgesOnCell": "nEdgesOnCell", "verticesOnCell": "nEdgesOnCell",
                "edgesOnEdge": "nEdgesOnEdge"}


def main():
    exe, compare, src, work, mode = sys.argv[1:6]
    scalars = sys.argv[6] if len(sys.argv) > 6 else "qv"
    work = Path(work)
    work.mkdir(parents=True, exist_ok=True)
    stem = Path(src).stem
    cand, ref = work / f"{stem}.rt.cand.nc", work / f"{stem}.rt.ref.nc"
    for p in (cand, ref):  # no file of an earlier run can stand in for this one's
        p.unlink(missing_ok=True)
    cmd = [exe, "roundtrip", src, str(cand)]
    names = list(MESH)
    if mode not in ("grid", "init"):
        print(f"FAIL mode {mode!r}: grid or init")
        return 2
    if mode == "init":
        cmd += ["--init", "--scalars", scalars]
        names += INIT + scalars.split(",")
    if mode == "grid":
        with nc.Dataset(src) as s:
            names += [n for n in GRID_OPTIONAL if n in s.variables]
    subprocess.run(cmd, check=True)
    zeroed = 0
    with nc.Dataset(src) as s:
        s.set_auto_maskandscale(False)
        s.set_auto_chartostring(False)
        with nc.Dataset(ref, "w", format=s.data_model) as d:
            d.set_fill_off()
            for n in names:
                for dn in s[n].dimensions:
                    if dn not in d.dimensions:
                        dim = s.dimensions[dn]
                        d.createDimension(dn, None if dim.isunlimited() else len(dim))
            for a in s.ncattrs():
                d.setncattr(a, s.getncattr(a))
            for n in names:
                v = s[n]
                o = d.createVariable(n, v.dtype, v.dimensions)
                for a in v.ncattrs():
                    o.setncattr(a, v.getncattr(a))
                o.set_auto_chartostring(False)
                vals = np.array(v[...])
                if n in UNUSED_SLOTS:
                    cnt = s[UNUSED_SLOTS[n]][:]
                    mask = np.arange(vals.shape[1])[None, :] >= cnt[:, None]
                    zeroed += int(np.count_nonzero(vals[mask]))
                    vals[mask] = 0
                o[...] = vals
    copy_text_attributes(src, ref, names)
    print(f"reference: {len(names)} variables, {zeroed} nonzero unused connectivity slots set to 0")
    r = subprocess.run([sys.executable, compare, str(ref), str(cand)], capture_output=True, text=True)
    print(r.stdout, end="")
    print(r.stderr, end="", file=sys.stderr)
    if r.returncode != 0 or "SAME" not in r.stdout:
        return 1
    with nc.Dataset(cand) as c:
        if sorted(c.variables) != sorted(names):
            print("FAIL candidate variables", sorted(set(c.variables) ^ set(names)))
            return 1
    # mpas_compare.py leaves out the value of the run-dependent file_id; a round trip must keep it (and parent_id)
    # exactly, byte for byte
    gs = {k: (t, raw) for k, t, raw in mpas_compare.classic_header(src)["gatts"]}
    gc = {k: (t, raw) for k, t, raw in mpas_compare.classic_header(cand)["gatts"]}
    for a in ("file_id", "parent_id"):
        if gc.get(a) != gs.get(a):
            print(f"FAIL global attribute {a}: {gc.get(a)!r}, source {gs.get(a)!r}")
            return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
