#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""The refusals of `cppmpas run` with config_do_restart (src/io/mpas_restart.cpp) on variants of the reference MW
restart file, and the start-time check of a cold start: each variant must fail with the given message.

  * a restart-stream field the run reads (ru, a coupled diag field; ru_p, a dycore field; Time) is missing;
  * the file holds a restart-stream field of an active package with no port (rthblten, the PBL tendency a run with
    bl_ysu would write);
  * a cold start (config_do_restart = false) whose input file does not hold config_start_time.

usage: check_restart.py <cppmpas> <restart.nc> <restart namelist> <init.nc> <scratch dir>
"""
import re
import subprocess
import sys
from pathlib import Path

import netCDF4 as nc
import numpy as np


def variant(src, dst, drop=(), add=None):
    """Copy src to dst without the variables in drop; add = (name, dims) adds a zero double variable."""
    a = nc.Dataset(src)
    a.set_auto_maskandscale(False)
    b = nc.Dataset(dst, "w", format=a.data_model)
    for d, dim in a.dimensions.items():
        b.createDimension(d, None if dim.isunlimited() else len(dim))
    for n, v in a.variables.items():
        if n in drop:
            continue
        w = b.createVariable(n, v.datatype, v.dimensions)
        w.setncatts({k: v.getncattr(k) for k in v.ncattrs()})
        w.set_auto_maskandscale(False)
        w[:] = v[:]
    if add:
        w = b.createVariable(add[0], "f8", add[1])
        w[:] = np.zeros(w.shape)
    b.setncatts({k: a.getncattr(k) for k in a.ncattrs()})
    b.close()
    a.close()


def main():
    exe, rst, nml, init, work = sys.argv[1], Path(sys.argv[2]), Path(sys.argv[3]), Path(sys.argv[4]), Path(sys.argv[5])
    work.mkdir(parents=True, exist_ok=True)
    cold = work / "namelist.cold_at_restart_time"
    text = nml.read_text()
    text, n = re.subn(r"config_do_restart = true", "config_do_restart = false", text)
    if n != 1:
        raise SystemExit("test setup: config_do_restart = true not found once in the restart namelist")
    cold.write_text(text)
    cases = [
        ("missing_ru", dict(drop=("ru",)), nml, "variable ru of the restart stream is missing"),
        ("missing_ru_p", dict(drop=("ru_p",)), nml, "variable ru_p of the restart stream is missing"),
        ("missing_Time", dict(drop=("Time",)), nml, "variable Time of the restart stream is missing"),
        ("active_pbl_tendency", dict(add=("rthblten", ("Time", "nCells", "nVertLevels"))), nml,
         "MPAS reads variable rthblten from the restart stream; CPPMPAS has no port for it"),
        ("cold_start_time", None, cold, "does not contain the time 0000-01-01_00:30:06"),
    ]
    bad = 0
    for name, how, namelist, want in cases:
        f = init
        if how is not None:
            f = work / f"restart.{name}.nc"
            variant(rst, f, **how)
        out = work / f"out.{name}"
        out.mkdir(exist_ok=True)
        r = subprocess.run([exe, "run", str(f), str(out), "--namelist", str(namelist), "--steps", "0", "--every", "0"],
                           capture_output=True, text=True)
        msg = [l for l in r.stderr.splitlines() if l.startswith("ERROR:")]
        ok = r.returncode == 1 and len(msg) == 1 and want in msg[0]
        print(f"{'ok  ' if ok else 'FAIL'} {name}: {msg[0] if msg else r.stderr.strip()[-300:]}")
        bad += not ok
        if how is not None:
            f.unlink()
        for p in out.iterdir():
            p.unlink()
        out.rmdir()
    cold.unlink()
    print(f"{len(cases) - bad}/{len(cases)} restart cases as expected")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
