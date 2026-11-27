#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Run one model run and record its timing.

usage: timed.py <name> <run dir> <timing dir> <command...>

The command runs in <run dir> with stdout and stderr in <run dir>/run.log. <timing dir>/<name>.json gets the wall
time (time.monotonic around the process), the children's user and system CPU time and maximum RSS, and, for each
NetCDF file the run wrote in <run dir>, the time after the start at which it was last modified. The output times
split a run into its initialisation (up to the first output) and the stretches between outputs; neither CPPMPAS
nor this script times anything inside the process. The exit status is the command's.
"""
import json
import resource
import subprocess
import sys
import time
from pathlib import Path


def main():
    name, rundir, timing, cmd = sys.argv[1], Path(sys.argv[2]), Path(sys.argv[3]), sys.argv[4:]
    timing.mkdir(parents=True, exist_ok=True)
    start, t0 = time.time(), time.monotonic()
    with open(rundir / "run.log", "wb") as log:
        rc = subprocess.run(cmd, cwd=rundir, stdout=log, stderr=subprocess.STDOUT).returncode
    wall = time.monotonic() - t0
    ru = resource.getrusage(resource.RUSAGE_CHILDREN)
    files = {p.name: round(p.stat().st_mtime - start, 2) for p in rundir.glob("*.nc") if not p.is_symlink()}
    d = {"name": name, "rc": rc, "wall_s": round(wall, 2), "user_s": round(ru.ru_utime, 1),
         "sys_s": round(ru.ru_stime, 1), "maxrss_mb": ru.ru_maxrss // 1024, "cmd": cmd, "start_unix": start,
         "file_mtime_s": dict(sorted(files.items(), key=lambda kv: kv[1]))}
    (timing / f"{name}.json").write_text(json.dumps(d, indent=1) + "\n")
    print(f"{name} rc={rc} wall={wall:.1f} s")
    return rc


if __name__ == "__main__":
    sys.exit(main())
