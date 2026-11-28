#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""A restart test (RGR, RGR24) runs from the forecast checkpoint of its restart reference (mpas842-dbl-rgr-001,
mpas842-dbl-rgr24-001) with that reference's restart/ namelist and streams. That reference is its base case's recipe
(RG of mpas842-dbl-rg-004, RG24 of mpas842-dbl-rg24-001) only if its forecast is the base case's: the init file by
tools/mpas_compare.py (the random file_id init_atmosphere writes excepted), the forecast namelist and streams byte for
byte, and any further files given (NetCDF by tools/mpas_compare.py, their parent_id left out: its first line is the
random file_id of that run's init file); and its restart namelist differs from the
forecast's in config_start_time, config_run_duration and config_do_restart only. The copy of the restart namelist in
tests/namelists (read by the namelist-table check) must be the reference's, byte for byte.

usage: check_restart_recipe.py <base case dir> <restart case dir> <namelist copy> [file ...]
"""
import contextlib
import io
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
import mpas_compare  # noqa: E402

SAME = ["init/x1.10242.init.nc", "forecast/namelist.atmosphere", "forecast/streams.atmosphere"]
RESTART_ALLOWED = {"config_run_duration", "config_start_time", "config_do_restart"}


def options(path):
    """{(group, option): line} of a namelist"""
    out, group = {}, None
    for line in path.read_text().splitlines():
        s = line.strip()
        if s.startswith("&"):
            group = s[1:]
        elif "=" in s:
            out[(group, s.split("=", 1)[0].strip())] = s
    return out


def same_file(p, q, random=frozenset()):
    """p and q identical; NetCDF by tools/mpas_compare.py, which also leaves out the global attributes in random"""
    if not p.is_file() or not q.is_file():
        return False
    if p.suffix == ".nc":
        kept = mpas_compare.RANDOM_ATTRS
        mpas_compare.RANDOM_ATTRS = kept | random
        try:
            with contextlib.redirect_stdout(io.StringIO()) as out, contextlib.redirect_stderr(out):
                return mpas_compare.main([str(p), str(q)]) == 0
        finally:
            mpas_compare.RANDOM_ATTRS = kept
    return p.read_bytes() == q.read_bytes()


def main():
    base, rst, copy = (Path(x) for x in sys.argv[1:4])
    RG, RGR = f"{base.parent.name}/{base.name}", f"{rst.parent.name}/{rst.name}"
    bad = 0
    # a further file is a model output: its parent_id starts with the random file_id of its own trial's init file
    for f in SAME + sys.argv[4:]:
        if same_file(base / f, rst / f, frozenset() if f in SAME else frozenset({"parent_id"})):
            print(f"ok   {f}: {RG} = {RGR}")
        else:
            print(f"FAIL {f}: {RG} and {RGR} differ (or one is missing)")
            bad += 1
    a, b = rst / "forecast" / "namelist.atmosphere", rst / "restart" / "namelist.atmosphere"
    if a.is_file() and b.is_file():
        oa, ob = options(a), options(b)
        diff = sorted(k for k in oa.keys() | ob.keys() if oa.get(k) != ob.get(k))
        other = [k for k in diff if k[1] not in RESTART_ALLOWED]
        if other or not oa or ("restart", "config_do_restart") not in diff:
            print(f"FAIL {RGR}: forecast and restart namelists differ in {other or diff}")
            bad += 1
        else:
            print(f"ok   {RGR}: forecast and restart namelists differ only in {[k[1] for k in diff]}")
    else:
        print(f"FAIL {RGR}: {a} or {b} is missing")
        bad += 1
    if same_file(copy, b):
        print(f"ok   {copy} = {RGR}/restart/namelist.atmosphere")
    else:
        print(f"FAIL {copy} is not {RGR}/restart/namelist.atmosphere")
        bad += 1
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
