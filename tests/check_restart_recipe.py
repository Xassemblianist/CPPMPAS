#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""The RGR restart test runs from the forecast checkpoint of mpas842-dbl-rgr-001 with that reference's restart/
namelist and streams. That reference is RG's recipe only if its forecast is RG's (the init file of
mpas842-dbl-rg-004 by tools/mpas_compare.py, the random file_id init_atmosphere writes excepted; the forecast
namelist and streams byte for byte), and its restart namelist differs from the forecast's in config_start_time,
config_run_duration and config_do_restart only. tests/namelists/RGR_restart/namelist.atmosphere (read by the
namelist-table check) must be the reference's restart namelist, byte for byte.

usage: check_restart_recipe.py <CPPMPAS_REFS_DIR> <source dir>
"""
import contextlib
import io
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
import mpas_compare  # noqa: E402

RG, RGR = "mpas842-dbl-rg-004/RG", "mpas842-dbl-rgr-001/RGR"
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


def same_file(p, q):
    if not p.is_file() or not q.is_file():
        return False
    if p.suffix == ".nc":
        with contextlib.redirect_stdout(io.StringIO()) as out, contextlib.redirect_stderr(out):
            return mpas_compare.main([str(p), str(q)]) == 0
    return p.read_bytes() == q.read_bytes()


def main():
    refs, src = Path(sys.argv[1]), Path(sys.argv[2])
    bad = 0
    for f in SAME:
        if same_file(refs / RG / f, refs / RGR / f):
            print(f"ok   {f}: {RG} = {RGR}")
        else:
            print(f"FAIL {f}: {RG} and {RGR} differ (or one is missing)")
            bad += 1
    a, b = refs / RGR / "forecast" / "namelist.atmosphere", refs / RGR / "restart" / "namelist.atmosphere"
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
    copy = src / "tests" / "namelists" / "RGR_restart" / "namelist.atmosphere"
    if same_file(copy, b):
        print(f"ok   {copy.relative_to(src)} = {RGR}/restart/namelist.atmosphere")
    else:
        print(f"FAIL {copy.relative_to(src)} is not {RGR}/restart/namelist.atmosphere")
        bad += 1
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
