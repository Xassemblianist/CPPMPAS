#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""qv_init on a restart: MPAS reads it with the restart stream (Registry.xml:661, mesh pool) and no routine of the
atmosphere core sets it, so the next restart file holds the values the run read. RGR's checkpoint has qv_init = 0
(init_atmosphere's real-data case does not set it), so RGR alone cannot tell "read and written back" from "the Registry
default". This test gives the checkpoint a qv_init profile that is not 0 and checks the run's next restart file.

  check_qv_init.py make <checkpoint> <copy>
      copies the checkpoint and sets qv_init(k) = 1e-3 / k (the file is otherwise unchanged)
  check_qv_init.py check <copy> <reference restart> <candidate restart>
      qv_init of the candidate must be the copy's, byte for byte, and every other variable of the candidate must be
      the reference's (tools/mpas_compare.py): qv_init changes nothing else
"""
import shutil
import sys
from pathlib import Path

import netCDF4 as nc
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
import mpas_compare  # noqa: E402


def make(src, dst):
    Path(dst).unlink(missing_ok=True)
    shutil.copyfile(src, dst)
    with nc.Dataset(dst, "a") as d:
        v = d.variables["qv_init"]
        n = v.shape[-1]
        v[...] = (1e-3 / np.arange(1, n + 1, dtype=np.float64)).reshape(v.shape)
    return 0


def raw(path, name):
    with nc.Dataset(path) as d:
        d.set_auto_maskandscale(False)
        return np.asarray(d.variables[name][...]).tobytes()


def check(copy, ref, cand):
    bad = 0
    want, got = raw(copy, "qv_init"), raw(cand, "qv_init")
    if not np.frombuffer(want, dtype=np.float64).any():
        print(f"FAIL {copy}: qv_init is 0 (the copy was not made)")
        bad += 1
    if got != want:
        print(f"FAIL {cand}: qv_init is not the checkpoint's")
        bad += 1
    else:
        print("ok   qv_init of the next restart file is the checkpoint's")
    with nc.Dataset(ref) as d:
        others = [n for n in d.variables if n != "qv_init"]
    code = mpas_compare.main([str(ref), str(cand), "--vars", ",".join(others)])
    if code != 0:
        print(f"FAIL {cand}: the other {len(others)} variables differ from {ref}")
        bad += 1
    return 1 if bad else 0


if __name__ == "__main__":
    if len(sys.argv) == 4 and sys.argv[1] == "make":
        sys.exit(make(sys.argv[2], sys.argv[3]))
    if len(sys.argv) == 5 and sys.argv[1] == "check":
        sys.exit(check(*sys.argv[2:5]))
    sys.exit(__doc__)
