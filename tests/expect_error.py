#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Run a command that must fail with exit status 1 and an 'ERROR:' message containing the given text, or (--warning)
succeed with exit status 0 and an 'MPAS warning:' line containing it.

usage: expect_error.py [--warning] "<expected text>" <command> [args...]
"""
import subprocess
import sys


def main():
    warning = len(sys.argv) > 1 and sys.argv[1] == "--warning"
    args = sys.argv[2:] if warning else sys.argv[1:]
    expected, cmd = args[0], args[1:]
    r = subprocess.run(cmd, capture_output=True, text=True)
    msg = r.stderr.strip()
    print(msg)
    if warning:
        if r.returncode != 0 or not any(l.startswith("MPAS warning:") and expected in l for l in msg.splitlines()):
            print(f"FAIL: exit {r.returncode}, expected success with an MPAS warning containing {expected!r}")
            return 1
        return 0
    if r.returncode != 1 or "ERROR:" not in msg or expected not in msg:
        print(f"FAIL: exit {r.returncode}, expected an error containing {expected!r}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
