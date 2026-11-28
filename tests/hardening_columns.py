#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Columns that reach a situation MPAS leaves undefined, through the column programs of the CPU or the GPU build
(tests/*_column/test_*_column.cpp): each must stop with the named Error (exit 1, the message on stderr), and the same
program on an ordinary column must run. With --vs <program of the other build>, the ordinary column's outputs of
the two programs must be the same bytes (the GPU build against the CPU build).

  hardening_columns.py wsm6_1km <test_wsm6_column> <work dir> [--vs <cpu program>]
      WSM6 with l_diags on a 55-level column 550 m deep: no level midpoint at or above 1 km AGL, where
      mpas_atmphys_driver_microphysics.F:722-724 read zp(kte+1) past the column (Error); a 30 km column runs.
  hardening_columns.py mcica <test_radiation_column> <work dir> <tables dir> <lw tables> <sw tables> [--vs ...]
      RRTMG on 2-level columns (maximum-random overlap): 3 radiation layers, where the KISS seed of MCICA_SUBCOL
      reads the bottom four (Error); 3 levels (4 layers) run.
  hardening_columns.py noah <test_noah_column> <work dir> <tables dir> [--vs ...]
      Noah (USGS) on a restart with a land point of vegetation category 0 (REDPRM reads its tables at 0: Error), a
      cold start with soil category 31 (above SOILPARM's 19: Error), a cold start without config_sfc_snowalbedo and
      vegetation category 0 (MAXALB(0): Error), a vegetation table whose cropland has NROOT = 0 on a vegetated point
      (CANRES reads ZSOIL(0): Error); the unmodified points run.
"""
import os
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PY = sys.executable


def run(cmd, cwd=None):
    r = subprocess.run([str(c) for c in cmd], cwd=cwd, capture_output=True, text=True)
    return r.returncode, r.stdout + r.stderr


def make(cmd):
    rc, out = run([PY, "-I"] + cmd)
    if rc != 0:
        sys.exit(f"FAIL: making the input: {out}")


fails = 0


def expect_error(prog_args, msg, what):
    global fails
    rc, out = run(prog_args)
    ok = rc == 1 and msg in out
    print(f"{'ok  ' if ok else 'FAIL'} {what}: exit {rc}, {out.strip().splitlines()[-1] if out.strip() else ''}")
    if not ok:
        fails += 1


def expect_run(prog_args, out_file, what, vs_args=None, vs_out=None):
    global fails
    rc, out = run(prog_args)
    ok = rc == 0 and Path(out_file).is_file() and Path(out_file).stat().st_size > 0
    print(f"{'ok  ' if ok else 'FAIL'} {what}: exit {rc}")
    if not ok:
        print(out)
        fails += 1
        return
    if vs_args:
        rc2, out2 = run(vs_args)
        same = rc2 == 0 and Path(out_file).read_bytes() == Path(vs_out).read_bytes()
        print(f"{'ok  ' if same else 'FAIL'} {what}: the same bytes as the other build's program")
        if not same:
            print(out2)
            fails += 1


def main():
    args = sys.argv[1:]
    vs = None
    if "--vs" in args:
        i = args.index("--vs")
        vs = args[i + 1]
        args = args[:i] + args[i + 2:]
    if len(args) < 3:
        sys.exit(__doc__)
    case, prog, work = args[0], args[1], Path(args[2])
    shutil.rmtree(work, ignore_errors=True)
    work.mkdir(parents=True)
    if case == "wsm6_1km":
        cols = HERE / "wsm6_column" / "columns.py"
        make([cols, "make", work / "shallow.bin", 3, 55, 1, 55, 60.0, 7, 550.0])
        expect_error([prog, work / "shallow.bin", work / "shallow.out", "l_diags"],
                     "refl10cm_1km: no level midpoint at or above 1 km AGL", "WSM6 l_diags, 55 levels 550 m deep")
        make([cols, "make", work / "deep.bin", 7, 55, 2, 55, 60.0, 8])
        expect_run([prog, work / "deep.bin", work / "deep.out", "l_diags"], work / "deep.out",
                   "WSM6 l_diags, 55 levels 30 km deep",
                   vs and [vs, work / "deep.bin", work / "deep_vs.out", "l_diags"], work / "deep_vs.out")
    elif case == "mcica":
        tables, lw, sw = args[3:6]
        cols = HERE / "radiation_column" / "columns.py"
        # 3 columns: the regimes 0-4 (clear, liquid, deep mixed-phase, cirrus, multi-layer), maximum-random
        make([cols, "make", work / "thin.bin", 3, 2, 1, 3, 2, 0, 1])
        expect_error([prog, work / "thin.bin", work / "thin.out", tables, lw, sw],
                     "MCICA_SUBCOL: fewer than 4 layers", "RRTMG on 2-level cloudy columns")
        make([cols, "make", work / "ok.bin", 3, 3, 1, 3, 2, 0, 1])
        expect_run([prog, work / "ok.bin", work / "ok.out", tables, lw, sw], work / "ok.out",
                   "RRTMG on 3-level cloudy columns",
                   vs and [vs, work / "ok.bin", work / "ok_vs.out", tables, lw, sw], work / "ok_vs.out")
    elif case == "noah":
        tables = Path(args[3])
        cols = HERE / "noah_column" / "columns.py"

        def point0(src, dst, ivgtyp=None, isltyp=None):
            b = bytearray(Path(src).read_bytes())
            nC = int.from_bytes(b[0:4], "little", signed=True)
            o = 40 + 8 + 64  # hdr[10], dt, mminlu
            if ivgtyp is not None:
                b[o:o + 4] = int(ivgtyp).to_bytes(4, "little", signed=True)
            if isltyp is not None:
                b[o + 4 * nC:o + 4 * nC + 4] = int(isltyp).to_bytes(4, "little", signed=True)
            Path(dst).write_bytes(bytes(b))

        make([cols, "make", work / "restart.bin", 12, 2, 900.0, 3, "USGS", "r"])
        point0(work / "restart.bin", work / "restart_veg0.bin", ivgtyp=0)
        expect_error([prog, work / "restart_veg0.bin", work / "o1.out", tables],
                     "a soil, vegetation or slope category below 1", "Noah restart, vegetation category 0")
        make([cols, "make", work / "cold.bin", 12, 2, 900.0, 4, "USGS", ""])
        point0(work / "cold.bin", work / "cold_soil31.bin", isltyp=31)
        expect_error([prog, work / "cold_soil31.bin", work / "o2.out", tables], "ISLTYP above the soil table",
                     "Noah cold start, soil category 31")
        make([cols, "make", work / "cold_s.bin", 12, 2, 900.0, 5, "USGS", "s"])
        point0(work / "cold_s.bin", work / "cold_s_veg0.bin", ivgtyp=0)
        expect_error([prog, work / "cold_s_veg0.bin", work / "o3.out", tables], "IVGTYP outside the vegetation table",
                     "Noah cold start without config_sfc_snowalbedo, vegetation category 0")
        t0 = work / "tables_nroot0"
        t0.mkdir()
        for t in ("SOILPARM.TBL", "GENPARM.TBL"):
            shutil.copy(tables / t, t0 / t)
        veg = (tables / "VEGPARM.TBL").read_text().split("\n")
        # the USGS section: category 2 (Dryland Cropland and Pasture, point 0's) with NROOT 3 -> 0
        usgs = veg.index("USGS")
        k = next(i for i in range(usgs, len(veg)) if veg[i].startswith("2,"))
        fields = veg[k].split(",")
        assert fields[2].strip() == "3", veg[k]
        fields[2] = fields[2].replace("3", "0")
        veg[k] = ",".join(fields)
        (t0 / "VEGPARM.TBL").write_text("\n".join(veg))
        expect_error([prog, work / "cold.bin", work / "o4.out", t0], "NROOT = 0 (CANRES reads ZSOIL(0))",
                     "Noah, cropland with NROOT = 0")
        # and with LAIMIN = LAIMAX = 0: XLAI = 0 skips CANRES, TRANSP is reached
        t1 = work / "tables_nroot0_lai0"
        shutil.copytree(t0, t1)
        fields[8] = fields[8].replace(fields[8].strip(), "0.00")
        fields[9] = fields[9].replace(fields[9].strip(), "0.00")
        veg[k] = ",".join(fields)
        (t1 / "VEGPARM.TBL").write_text("\n".join(veg))
        expect_error([prog, work / "cold.bin", work / "o5.out", t1], "NROOT = 0 (TRANSP divides by NROOT)",
                     "Noah, cropland with NROOT = 0 and no leaves")
        expect_run([prog, work / "cold.bin", work / "cold.out", tables], work / "cold.out", "Noah cold start",
                   vs and [vs, work / "cold.bin", work / "cold_vs.out", tables], work / "cold_vs.out")
        expect_run([prog, work / "restart.bin", work / "restart.out", tables], work / "restart.out", "Noah restart",
                   vs and [vs, work / "restart.bin", work / "restart_vs.out", tables], work / "restart_vs.out")
    else:
        sys.exit(f"unknown case {case}")
    print("PASS" if fails == 0 else f"FAIL ({fails})")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
