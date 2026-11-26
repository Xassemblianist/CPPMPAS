#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""tools/gen_namelist_table.py on variants of a reference namelist.atmosphere: a value it cannot read as MPAS's
namelist read does (set twice, a quoted number, a 32-bit integer overflow, a real that overflows) stops it with a
message; a '!' inside a quoted value is part of the value, not a comment.

usage: check_gen_namelist.py <gen_namelist_table.py> <Registry.xml> <reference namelist.atmosphere> <scratch dir>
"""
import re
import subprocess
import sys
from pathlib import Path


def sub(pattern, repl):
    def f(text):
        out, n = re.subn(pattern, repl, text, flags=re.M)
        if n != 1:
            raise SystemExit(f"test setup: pattern {pattern!r} matched {n} times")
        return out
    return f


OK = None
CASES = [
    ("reference", lambda t: t, OK, "CPPMPAS_NML_OPTION(\"damping\", \"config_xnutr\", \"real\", \"0.2\", \"0.1\")"),
    ("bang_in_quotes", sub(r"config_physics_suite = 'none'", "config_physics_suite = 'no!ne' ! comment"), OK,
     "\"config_physics_suite\", \"character\", \"mesoscale_reference\", \"no!ne\")"),
    ("set_twice", sub(r"^(\s*config_zd = 10500.0)$", r"\1\n    config_zd = 10500.0"), "config_zd is set twice", None),
    ("quoted_real", sub(r"config_xnutr = 0.1", "config_xnutr = '0.1'"), "a real value must not be quoted", None),
    ("int_range", sub(r"config_number_of_sub_steps = 6", "config_number_of_sub_steps = 2147483648"),
     "is outside a 32-bit integer", None),
    ("real_overflow", sub(r"config_xnutr = 0.1", "config_xnutr = 1e999"), "overflows a double", None),
]


def main():
    gen, registry, ref, work = sys.argv[1], sys.argv[2], Path(sys.argv[3]), Path(sys.argv[4])
    work.mkdir(parents=True, exist_ok=True)
    text = ref.read_text()
    bad = 0
    for name, f, err, want in CASES:
        nml = work / f"gen_namelist.{name}"
        nml.write_text(f(text))
        r = subprocess.run([sys.executable, "-I", gen, registry, str(nml)], capture_output=True, text=True)
        if err is None:
            ok = r.returncode == 0 and want in r.stdout
            msg = "table written" if ok else (r.stderr.strip() or "expected row missing")
        else:
            ok = r.returncode != 0 and err in r.stderr
            msg = r.stderr.strip() or f"exit {r.returncode}"
        print(f"{'ok  ' if ok else 'FAIL'} {name}: {msg}")
        bad += not ok
        nml.unlink()
    print(f"{len(CASES) - bad}/{len(CASES)} generator cases as expected")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
