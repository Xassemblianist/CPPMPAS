#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""The namelist reader of `cppmpas init` (src/core/namelist.cpp) on variants of a reference namelist.atmosphere.

Accepted variants (the same values written in other forms MPAS's namelist read takes) must pass the namelist and
stop at the missing init file; refused variants must fail with the given message. Every option of the Registry is
recognised; a value that is not a verified reference value, an unknown option, record or form is an error naming it.

usage: check_namelist.py <cppmpas> <reference namelist.atmosphere> <scratch dir>
"""
import re
import subprocess
import sys
from pathlib import Path

PASSED = "cannot open"  # the namelist was accepted: the init file (which does not exist) is read next


def sub(pattern, repl):
    def f(text):
        out, n = re.subn(pattern, repl, text, flags=re.M)
        if n != 1:
            raise SystemExit(f"test setup: pattern {pattern!r} matched {n} times")
        return out
    return f


def add(record, line):
    return sub(rf"^&{record}\s*$", f"&{record}\n    {line}")


CASES = [
    # accepted: the reference itself and other spellings of the same values
    ("reference", lambda t: t, PASSED),
    ("logical_forms", sub(r"config_do_restart = false", "config_do_restart = .F."), PASSED),
    ("real_exponent", sub(r"config_zd = 10500.0", "config_zd = 1.05d4"), PASSED),
    ("comment_commas", sub(r"config_xnutr = 0.1", "config_xnutr = 0.1,   ! a comment, with / and &"), PASSED),
    ("one_line_group", sub(r"&restart\n    config_do_restart = false\n/", "&restart config_do_restart = false /"),
     PASSED),
    ("two_per_line", sub(r"    config_zd = 10500.0\n    config_xnutr = 0.1",
                         "    config_zd = 10500.0, config_xnutr = 0.1"), PASSED),
    ("case_of_names", sub(r"config_h_ScaleWithMesh", "CONFIG_H_SCALEWITHMESH"), PASSED),
    ("default_written", add("restart", "config_do_DAcycling = false"), PASSED),
    # refused values
    ("do_restart", sub(r"config_do_restart = false", "config_do_restart = true"),
     "namelist option config_do_restart = true has no verified port"),
    ("do_restart_one_line", sub(r"&restart\n    config_do_restart = false\n/", "&restart config_do_restart = true /"),
     "namelist option config_do_restart = true has no verified port"),
    ("do_restart_dotted", sub(r"    config_do_restart = false", "    config_do_restart = .true. ,"),
     "namelist option config_do_restart = .true. has no verified port"),
    ("dacycling", add("restart", "config_do_DAcycling = true"),
     "namelist option config_do_DAcycling = true has no verified port"),
    ("xnutr", sub(r"config_xnutr = 0.1", "config_xnutr = 0.15"),
     "namelist option config_xnutr = 0.15 has no verified port"),
    ("default_refused", sub(r"    config_dt = 6.0\n", ""),
     "namelist option config_dt = 720.0 (Registry default) has no verified port"),
    ("physics_suite", sub(r"config_physics_suite = 'none'", "config_physics_suite = 'mesoscale_reference'"),
     "namelist option config_physics_suite = 'mesoscale_reference' has no verified port"),
    # refused forms
    ("unknown_key", add("damping", "config_not_an_option = 1"),
     "config_not_an_option is not an option of the namelist record &damping"),
    ("wrong_record", lambda t: add("restart", "config_xnutr = 0.1")(sub(r"    config_xnutr = 0.1\n", "")(t)),
     "config_xnutr is not an option of the namelist record &restart"),
    ("unknown_record", lambda t: t + "\n&no_such_record\n/\n", "&no_such_record is not a namelist record"),
    ("quoted_real", sub(r"config_xnutr = 0.1", "config_xnutr = '0.1'"), "config_xnutr = '0.1' is not a real"),
    ("unquoted_character", sub(r"config_physics_suite = 'none'", "config_physics_suite = none"),
     "config_physics_suite = none is not a quoted character value"),
    ("bad_logical", sub(r"config_do_restart = false", "config_do_restart = yes"),
     "config_do_restart = yes is not a logical"),
    ("bad_integer", sub(r"config_number_of_sub_steps = 6", "config_number_of_sub_steps = 6.0"),
     "config_number_of_sub_steps = 6.0 is not an integer"),
    ("unclosed", sub(r"config_do_restart = false\n/", "config_do_restart = false\n"),
     "group &restart is not closed with '/'"),
    ("unclosed_at_end", sub(r"config_physics_suite = 'none'\n/", "config_physics_suite = 'none'\n"),
     "group &physics is not closed with '/'"),
    ("bad_real", sub(r"config_xnutr = 0.1", "config_xnutr = 0.1x"), "config_xnutr = 0.1x is not a real"),
    ("set_twice", add("damping", "config_zd = 10500.0"), "option config_zd is set twice"),
    ("array", sub(r"config_zd = 10500.0", "config_zd = 10500.0 2.0"), "config_zd has more than one value"),
    ("repeat_count", sub(r"config_zd = 10500.0", "config_zd = 1*10500.0"), "repeat counts are not handled"),
    ("outside_group", lambda t: "config_zd = 1.0\n" + t, "text outside a namelist group"),
]


def main():
    exe, ref, work = sys.argv[1], Path(sys.argv[2]), Path(sys.argv[3])
    work.mkdir(parents=True, exist_ok=True)
    text = ref.read_text()
    bad = 0
    for name, f, want in CASES:
        nml = work / f"namelist.{name}"
        nml.write_text(f(text))
        r = subprocess.run([exe, "init", str(work / "no_such_init.nc"), str(work / "x.nc"), "--fields", "dss",
                            "--namelist", str(nml)], capture_output=True, text=True)
        msg = r.stderr.strip()
        ok = r.returncode == 1 and "ERROR:" in msg and want in msg
        if want != PASSED:
            ok = ok and PASSED not in msg
        print(f"{'ok  ' if ok else 'FAIL'} {name}: {msg}")
        bad += not ok
        nml.unlink()
    print(f"{len(CASES) - bad}/{len(CASES)} namelist cases as expected")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
