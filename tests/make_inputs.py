#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Build the test inputs derived from the NCAR test-case grids: synthetic init files and broken files.

No MPAS init file exists yet (the reference init_atmosphere runs are not done), so the init files here carry the
real horizontal mesh of the grid plus every variable of the atmosphere input stream that CPPMPAS reads, with the
dimensions, types and attributes init_atmosphere writes (core_init_atmosphere/Registry.xml:497-640: CDF-5,
Time records, StrLen = 64) and deterministic pseudo-random values. They test reading and writing, not values.

usage: make_inputs.py <cases_dir> <out_dir>
"""
import sys
from pathlib import Path

import netCDF4 as nc
import numpy as np

MESH = ["latCell", "lonCell", "xCell", "yCell", "zCell", "indexToCellID", "latEdge", "lonEdge", "xEdge", "yEdge",
        "zEdge", "indexToEdgeID", "latVertex", "lonVertex", "xVertex", "yVertex", "zVertex", "indexToVertexID",
        "cellsOnEdge", "nEdgesOnCell", "nEdgesOnEdge", "edgesOnCell", "edgesOnEdge", "weightsOnEdge", "dvEdge",
        "dcEdge", "angleEdge", "areaCell", "areaTriangle", "cellsOnCell", "verticesOnCell", "verticesOnEdge",
        "edgesOnVertex", "cellsOnVertex", "kiteAreasOnVertex", "fEdge", "fVertex", "meshDensity"]
# name: (C-order dims, NetCDF type)
INIT = {
    "edgeNormalVectors": (("nEdges", "R3"), "f8"), "localVerticalUnitVectors": (("nCells", "R3"), "f8"),
    "cellTangentPlane": (("nCells", "TWO", "R3"), "f8"), "nominalMinDc": ((), "f8"),
    "bdyMaskCell": (("nCells",), "i4"), "bdyMaskEdge": (("nEdges",), "i4"), "bdyMaskVertex": (("nVertices",), "i4"),
    "cf1": ((), "f8"), "cf2": ((), "f8"), "cf3": ((), "f8"), "zgrid": (("nCells", "nVertLevelsP1"), "f8"),
    "rdzw": (("nVertLevels",), "f8"), "dzu": (("nVertLevels",), "f8"), "rdzu": (("nVertLevels",), "f8"),
    "fzm": (("nVertLevels",), "f8"), "fzp": (("nVertLevels",), "f8"), "zxu": (("nEdges", "nVertLevels"), "f8"),
    "zz": (("nCells", "nVertLevels"), "f8"), "zb": (("nEdges", "TWO", "nVertLevelsP1"), "f8"),
    "zb3": (("nEdges", "TWO", "nVertLevelsP1"), "f8"), "dss": (("nCells", "nVertLevels"), "f8"),
    "deriv_two": (("nEdges", "TWO", "FIFTEEN"), "f8"), "coeffs_reconstruct": (("nCells", "maxEdges", "R3"), "f8"),
    "u_init": (("nVertLevels",), "f8"), "v_init": (("nVertLevels",), "f8"),
    "t_init": (("nCells", "nVertLevels"), "f8"),
    "u": (("Time", "nEdges", "nVertLevels"), "f8"), "w": (("Time", "nCells", "nVertLevelsP1"), "f8"),
    "rho": (("Time", "nCells", "nVertLevels"), "f8"), "theta": (("Time", "nCells", "nVertLevels"), "f8"),
    "rho_base": (("Time", "nCells", "nVertLevels"), "f8"), "theta_base": (("Time", "nCells", "nVertLevels"), "f8"),
    "surface_pressure": (("Time", "nCells"), "f8"),
    # read by MPAS but not used by the dycore (CPPMPAS skips them on purpose)
    "relhum": (("Time", "nCells", "nVertLevels"), "f8"), "Time": (("Time",), "f8"), "xland": (("nCells",), "f8"),
    # not in the atmosphere input stream at all (MPAS does not read them)
    "ter": (("nCells",), "f8"), "qv_init": (("nVertLevels",), "f8"),
}
SCALARS = ["qv", "qc", "qr", "qi", "qs", "qg"]


def copy_grid_var(src, dst, name, fill_unused=None):
    v = src[name]
    out = dst.createVariable(name, v.dtype, v.dimensions)
    for a in v.ncattrs():
        out.setncattr(a, v.getncattr(a))
    out[...] = v[...]


def write_init(grid_path, out_path, nvl, seed, extra=None, drop=(), modify=None, ntime=1):
    rng = np.random.default_rng(seed)
    with nc.Dataset(grid_path) as g, nc.Dataset(out_path, "w", format="NETCDF3_64BIT_DATA") as d:
        g.set_auto_mask(False)
        d.set_fill_off()
        for dn in ("nCells", "nEdges", "nVertices", "maxEdges", "maxEdges2", "TWO", "vertexDegree"):
            d.createDimension(dn, len(g.dimensions[dn]))
        d.createDimension("Time", None)
        d.createDimension("StrLen", 64)
        d.createDimension("R3", 3)
        d.createDimension("FIFTEEN", 15)
        d.createDimension("nVertLevels", nvl)
        d.createDimension("nVertLevelsP1", nvl + 1)
        d.setncattr("on_a_sphere", g.getncattr("on_a_sphere").strip())
        d.setncattr("sphere_radius", g.getncattr("sphere_radius"))
        d.setncattr("is_periodic", "NO")
        d.setncattr("x_period", 0.0)
        d.setncattr("y_period", 0.0)
        d.setncattr("file_id", "abcdefghij0123456789abcdefghij0123456789")
        d.setncattr("parent_id", "parent0123")
        d.setncattr("mesh_spec", "0.0")
        for n in MESH:
            if n in drop:
                continue
            copy_grid_var(g, d, n)
        times = ["2000-01-01_00:00:00", "2000-01-01_00:30:00"][:ntime]
        for n in ("xtime", "initial_time"):
            dims = ("Time", "StrLen") if n == "xtime" else ("StrLen",)
            v = d.createVariable(n, "S1", dims)
            v.setncattr("units", "YYYY-MM-DD_hh:mm:ss")
            if n == "xtime":
                for t, s in enumerate(times):
                    v[t, :] = nc.stringtochar(np.array([s.ljust(64)], "S64"))[0]
            else:
                v[:] = nc.stringtochar(np.array([times[0].ljust(64)], "S64"))[0]
        allvars = dict(INIT)
        for s in SCALARS:
            allvars[s] = (("Time", "nCells", "nVertLevels"), "f8")
        if extra:
            allvars.update(extra)
        for n, (dims, typ) in allvars.items():
            if n in drop:
                continue
            shape = [ntime if dn == "Time" else len(d.dimensions[dn]) for dn in dims]
            v = d.createVariable(n, typ, dims)
            v.setncattr("units", "test units")
            v.setncattr("long_name", f"synthetic {n}")
            if typ == "i4":
                v[...] = np.zeros(shape, "i4")
            else:
                v[...] = rng.standard_normal(shape).astype(typ)
        if modify:
            modify(d)


def main():
    cases, out = Path(sys.argv[1]), Path(sys.argv[2])
    out.mkdir(parents=True, exist_ok=True)
    mw = cases / "mountain_wave" / "mtn_wave_grid.nc"
    sc = cases / "supercell" / "supercell_grid.nc"
    write_init(mw, out / "mw_init.nc", 4, 1)
    write_init(sc, out / "sc_init.nc", 3, 2)

    # broken inputs (each must be rejected with a clear message)
    def grid_variant(name, fmt=None, modify=None, drop=()):
        with nc.Dataset(mw) as g:
            g.set_auto_mask(False)
            with nc.Dataset(out / name, "w", format=fmt or "NETCDF3_64BIT_OFFSET") as d:
                for dn, dim in g.dimensions.items():
                    d.createDimension(dn, None if dim.isunlimited() else len(dim))
                for a in g.ncattrs():
                    d.setncattr(a, g.getncattr(a))
                for n in MESH:
                    if n not in drop:
                        copy_grid_var(g, d, n)
                if modify:
                    modify(d)

    grid_variant("err_periodic.nc", modify=lambda d: d.setncattr("is_periodic", "YES"))
    grid_variant("err_missing_dvEdge.nc", drop=("dvEdge",))

    def float_dvedge(d):
        d.renameVariable("dvEdge", "dvEdge_old")
        v = d.createVariable("dvEdge", "f4", ("nEdges",))
        v[:] = d["dvEdge_old"][:].astype("f4")
    grid_variant("err_float_dvEdge.nc", modify=float_dvedge)
    grid_variant("err_netcdf4.nc", fmt="NETCDF4")

    def perm_ids(d):
        a = d["indexToEdgeID"][:]
        a[[0, 1]] = a[[1, 0]]
        d["indexToEdgeID"][:] = a
    grid_variant("err_index_ids.nc", modify=perm_ids)

    def boundary(d):
        a = d["cellsOnCell"][:]
        a[5, 2] = 0
        d["cellsOnCell"][:] = a
    grid_variant("err_boundary.nc", modify=boundary)

    def bdymask(d):
        a = d["bdyMaskCell"][:]
        a[7] = 1
        d["bdyMaskCell"][:] = a
    write_init(mw, out / "err_bdymask.nc", 2, 3, modify=bdymask)
    write_init(mw, out / "err_two_records.nc", 2, 4, ntime=2)
    write_init(mw, out / "err_missing_zz.nc", 2, 5, drop=("zz",))
    write_init(mw, out / "err_missing_qr.nc", 2, 6, drop=("qr",))
    write_init(mw, out / "err_bad_dims_u.nc", 2, 7, extra={"u": (("Time", "nVertLevels", "nEdges"), "f8")})
    print("inputs written to", out)


if __name__ == "__main__":
    main()
