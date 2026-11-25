#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Mutation check for tools/mpas_compare.py: each one-line bug below must make the unit tests fail.

Run: python3 tests/compare/mutate.py   (exit 0 only if every mutant is killed and the original passes)
"""
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "tools" / "mpas_compare.py"

MUTANTS = [
    ("value bytes ignored", "if ra.tobytes() != rb.tobytes():", "if False:"),
    ("records stop early", "for t in range(var.shape[0]):", "for t in range(min(1, var.shape[0])):"),
    ("dimensions ignored", "if da != db:", "if False:"),
    ("global attributes ignored", "if attrs(a, WALL_CLOCK_ATTRS) != attrs(b, WALL_CLOCK_ATTRS):", "if False:"),
    ("all global attributes ignored", "if attrs(a, WALL_CLOCK_ATTRS) != attrs(b, WALL_CLOCK_ATTRS):",
     "if attrs(a, set(a.ncattrs())) != attrs(b, set(b.ncattrs())):"),
    ("variable attributes ignored", "if attrs(x) != attrs(y):", "if False:"),
    ("missing variables ignored", "for name in sorted(names_a - names_b):", "for name in []:"),
    ("compared count lost", "compared += 1", "pass"),
    ("nothing-compared accepted", "if compared == 0:", "if False:"),
    ("missing files ignored", 'problems = [{"file": n, "variable": None, "what": "file missing in candidate"}',
     'problems = [] and [{"file": n, "variable": None, "what": "file missing in candidate"}'),
    ("ulp off by sign fold", "return np.where(signed < 0, bias - signed - 1, signed)",
     "return np.where(signed < 0, bias - signed, signed)"),
    ("first index wrong", "first = np.unravel_index(int(where[0]), a.shape) if a.ndim else ()",
     "first = np.unravel_index(int(where[-1]), a.shape) if a.ndim else ()"),
    ("exit status lost", "    if report:\n        return 1", "    if report:\n        return 0"),
    ("selection error dropped", "if missing:", "if False:"),
    ("random file_id compared", '"file_id", "parent_id"})', '"parent_id"})'),
]


def run_tests(tools_dir, tmp):
    env = dict(os.environ, MPAS_COMPARE_DIR=str(tools_dir), MPAS_COMPARE_TMP=str(tmp), PYTHONDONTWRITEBYTECODE="1")
    proc = subprocess.run([sys.executable, "-m", "unittest", "discover", "-s", str(ROOT / "tests" / "compare")],
                          env=env, capture_output=True, text=True)
    return proc.returncode == 0


def main():
    tmp = Path(tempfile.mkdtemp(prefix="mpasmut", dir=os.environ.get("MPAS_COMPARE_TMP")))
    try:
        text = SOURCE.read_text()
        if not run_tests(SOURCE.parent, tmp):
            print("ORIGINAL FAILS")
            return 1
        survivors = 0
        for name, old, new in MUTANTS:
            if text.count(old) != 1:
                print(f"STALE   {name}: pattern not found exactly once")
                survivors += 1
                continue
            mdir = tmp / "mutant"
            mdir.mkdir(exist_ok=True)
            (mdir / "mpas_compare.py").write_text(text.replace(old, new))
            killed = not run_tests(mdir, tmp)
            print(f"{'KILLED ' if killed else 'SURVIVED'} {name}")
            survivors += not killed
        print(f"{len(MUTANTS) - survivors}/{len(MUTANTS)} mutants killed")
        return 1 if survivors else 0
    finally:
        shutil.rmtree(tmp)


if __name__ == "__main__":
    raise SystemExit(main())
