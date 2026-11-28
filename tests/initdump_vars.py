#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""The variables of a reference initdump.nc must be exactly the given list (so that the per-field tests cover every
field MPAS dumped), with one time record.

usage: initdump_vars.py INITDUMP.nc VAR[,VAR...]      exit 0 equal, 1 different, 2 input error
"""
import sys

import netCDF4 as nc


def main(argv):
    if len(argv) != 3:
        print(__doc__)
        return 2
    want = argv[2].split(",")
    if "" in want or len(set(want)) != len(want):
        print(f"ERROR the list {argv[2]!r} has an empty or a repeated name")
        return 2
    with nc.Dataset(argv[1]) as d:
        have = list(d.variables)
        records = len(d.dimensions["Time"]) if "Time" in d.dimensions else -1
    missing = sorted(set(have) - set(want))
    extra = sorted(set(want) - set(have))
    for v in missing:
        print(f"DIFF {v}: in the reference, not tested")
    for v in extra:
        print(f"DIFF {v}: tested, not in the reference")
    if records != 1:
        print(f"DIFF Time: {records} records")
    ok = not missing and not extra and records == 1
    print(f"{'SAME' if ok else 'DIFF'} {len(have)} variables")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
