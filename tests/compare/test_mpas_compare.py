# SPDX-License-Identifier: GPL-3.0-or-later
"""Tests for tools/mpas_compare.py on synthetic MPAS-like files. Run: python3 -m unittest discover tests/compare"""
import contextlib
import io
import os
import shutil
import sys
import tempfile
import unittest
from pathlib import Path

import netCDF4 as nc
import numpy as np

TOOLS = Path(os.environ.get("MPAS_COMPARE_DIR", Path(__file__).resolve().parents[2] / "tools"))
sys.path.insert(0, str(TOOLS))
import mpas_compare  # noqa: E402

N_TIME, N_CELLS, N_LEV = 3, 7, 4


def write(path, mutate=None, n_time=N_TIME):
    rng = np.random.default_rng(42)
    with nc.Dataset(path, "w", format="NETCDF3_64BIT_DATA") as ds:
        ds.createDimension("Time", None)
        ds.createDimension("StrLen", 64)
        ds.createDimension("nCells", N_CELLS)
        ds.createDimension("nVertLevels", N_LEV)
        ds.model_name = "mpas"
        ds.creation_date = "2026-11-25 00:00:00"
        ds.file_id = "0123456789"
        ds.config_dt = 450.0
        x = ds.createVariable("xtime", "S1", ("Time", "StrLen"))
        theta = ds.createVariable("theta", "f8", ("Time", "nCells", "nVertLevels"))
        theta.units = "K"
        w = ds.createVariable("w32", "f4", ("Time", "nCells", "nVertLevels"))
        ids = ds.createVariable("indexToCellID", "i4", ("nCells",))
        area = ds.createVariable("areaCell", "f8", ("nCells",))
        for t in range(n_time):
            label = f"0000-01-01_00:{t:02d}:00".ljust(64, "\x00")
            x[t] = np.frombuffer(label.encode(), dtype="S1")
            theta[t] = 300.0 + rng.standard_normal((N_CELLS, N_LEV))
            w[t] = rng.standard_normal((N_CELLS, N_LEV)).astype(np.float32)
        ids[:] = np.arange(1, N_CELLS + 1, dtype=np.int32)
        area[:] = 1.0e9 + rng.standard_normal(N_CELLS)
        if mutate:
            mutate(ds)


def run(*args):
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(io.StringIO()):
        code = mpas_compare.main([str(a) for a in args])
    return code, buf.getvalue()


def bump(value, ulps):
    return np.nextafter(value, np.inf) if ulps > 0 else np.nextafter(value, -np.inf)


class CompareTests(unittest.TestCase):
    def setUp(self):
        self.dir = Path(tempfile.mkdtemp(prefix="mpascmp", dir=os.environ.get("MPAS_COMPARE_TMP")))
        self.ref = self.dir / "ref.nc"
        write(self.ref)

    def tearDown(self):
        shutil.rmtree(self.dir)

    def cand(self, mutate=None, **kw):
        path = self.dir / "cand.nc"
        write(path, mutate, **kw)
        return path

    def test_identical(self):
        code, out = run(self.ref, self.cand())
        self.assertEqual(code, 0, out)
        self.assertIn("SAME 1 file(s), 5 variable(s)", out)

    def test_one_ulp_f8_reported_at_record(self):
        def m(ds):
            v = ds["theta"]
            v[2, 3, 1] = bump(v[2, 3, 1], 1)
        code, out = run(self.ref, self.cand(m))
        self.assertEqual(code, 1)
        self.assertIn("theta values record=2 xtime=0000-01-01_00:02:00 count=1/28 first=[3, 1]", out)
        self.assertIn("max_ulp=1", out)

    def test_first_of_several_differences(self):
        def m(ds):
            ds["theta"][1, 2, 3] = 0.0
            ds["theta"][1, 5, 0] = 0.0
        code, out = run(self.ref, self.cand(m))
        self.assertEqual(code, 1)
        self.assertIn("record=1 xtime=0000-01-01_00:01:00 count=2/28 first=[2, 3]", out)

    def test_one_ulp_f4(self):
        def m(ds):
            v = ds["w32"]
            v[0, 0, 0] = bump(np.float32(v[0, 0, 0]), -1)
        code, out = run(self.ref, self.cand(m))
        self.assertEqual(code, 1)
        self.assertIn("w32 values record=0", out)
        self.assertIn("max_ulp=1", out)

    def test_signed_zero(self):
        def zero(sign):
            def m(ds):
                ds["areaCell"][4] = sign * 0.0
            return m
        ref = self.dir / "ref0.nc"
        write(ref, zero(1.0))
        code, out = run(ref, self.cand(zero(-1.0)))
        self.assertEqual(code, 1)
        self.assertIn("areaCell values", out)
        self.assertIn("max_ulp=1", out)

    def test_nan_payload(self):
        def nan(bits):
            def m(ds):
                ds["areaCell"][0] = np.array([bits], dtype=np.uint64).view(np.float64)[0]
            return m
        ref = self.dir / "refn.nc"
        write(ref, nan(0x7FF8000000000000))
        code, _ = run(ref, self.cand(nan(0x7FF8000000000001)))
        self.assertEqual(code, 1)

    def test_integer_change(self):
        code, out = run(self.ref, self.cand(lambda ds: ds["indexToCellID"].__setitem__(0, 99)))
        self.assertEqual(code, 1)
        self.assertIn("indexToCellID values", out)

    def test_missing_variable(self):
        def m(ds):
            ds.renameVariable("theta", "theta_renamed")
        code, out = run(self.ref, self.cand(m))
        self.assertEqual(code, 1)
        self.assertIn("theta missing in candidate", out)
        self.assertIn("theta_renamed only in candidate", out)

    def test_extra_time_record(self):
        code, out = run(self.ref, self.cand(n_time=N_TIME + 1))
        self.assertEqual(code, 1)
        self.assertIn("dimensions", out)

    def test_variable_attribute(self):
        code, out = run(self.ref, self.cand(lambda ds: ds["theta"].setncattr("units", "C")))
        self.assertEqual(code, 1)
        self.assertIn("theta variable attributes", out)

    def test_wall_clock_attribute_ignored(self):
        def m(ds):
            ds.setncattr("creation_date", "2030-01-01 00:00:00")
            ds.setncattr("file_id", "abcdefghij")
        code, out = run(self.ref, self.cand(m))
        self.assertEqual(code, 0, out)

    def test_other_global_attribute(self):
        code, out = run(self.ref, self.cand(lambda ds: ds.setncattr("config_dt", 300.0)))
        self.assertEqual(code, 1)
        self.assertIn("global attributes", out)

    def test_selected_variables_only(self):
        def m(ds):
            ds["w32"][1, 1, 1] = 5.0
        code, out = run(self.ref, self.cand(m), "--vars", "theta")
        self.assertEqual(code, 0, out)
        code, out = run(self.ref, self.cand(m), "--vars", "theta,w32")
        self.assertEqual(code, 1)

    def test_selected_variable_absent_is_error(self):
        code, _ = run(self.ref, self.cand(), "--vars", "nonexistent")
        self.assertEqual(code, 2)

    def test_all_records(self):
        def m(ds):
            ds["theta"][0, 0, 0] = 0.0
            ds["theta"][2, 0, 0] = 0.0
        code, out = run(self.ref, self.cand(m), "--all-records")
        self.assertEqual(code, 1)
        self.assertIn("record=0", out)
        self.assertIn("record=2", out)
        _, out = run(self.ref, self.cand(m))
        self.assertNotIn("record=2", out)

    def test_directories(self):
        a, b = self.dir / "a", self.dir / "b"
        a.mkdir()
        b.mkdir()
        write(a / "history.nc")
        write(b / "history.nc")
        write(a / "restart.nc")
        code, out = run(a, b)
        self.assertEqual(code, 1)
        self.assertIn("restart.nc None file missing in candidate", out)
        write(b / "restart.nc")
        code, out = run(a, b)
        self.assertEqual(code, 0, out)

    def test_nothing_compared_fails(self):
        a, b = self.dir / "ea", self.dir / "eb"
        a.mkdir()
        b.mkdir()
        code, out = run(a, b)
        self.assertEqual(code, 1)
        self.assertIn("nothing was compared", out)

    def test_ulp_distance_across_zero(self):
        tiny = np.array([5e-324], dtype=np.float64)
        self.assertEqual(mpas_compare.ulp_distance(tiny, -tiny), 3)
        self.assertEqual(mpas_compare.ulp_distance(np.array([1.0]), np.array([np.nextafter(1.0, 2.0)])), 1)


if __name__ == "__main__":
    unittest.main()
