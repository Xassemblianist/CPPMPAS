#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Run a command that must fail with exit status 1 and an 'ERROR:' message containing the given text.

usage: expect_error.py "<expected text>" <command> [args...]
"""
import subprocess
import sys


def main():
    expected, cmd = sys.argv[1], sys.argv[2:]
    r = subprocess.run(cmd, capture_output=True, text=True)
    msg = r.stderr.strip()
    print(msg)
    if r.returncode != 1 or "ERROR:" not in msg or expected not in msg:
        print(f"FAIL: exit {r.returncode}, expected an error containing {expected!r}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
