#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Run one GPU test: skip it (exit status 77, CTest's SKIP_RETURN_CODE, reported as "Skipped", never "Passed") when
there is no CUDA build (NONE) or no CUDA device, else run the command and return its status.

usage: gpu_test.py <cppmpas_gpu | NONE> [--clean PATH]... device           the device check: kernels run on the GPU
       gpu_test.py <cppmpas_gpu | NONE> [--clean PATH]... no-device        with no visible device, cppmpas_gpu skips (77)
       gpu_test.py <cppmpas_gpu | NONE> [--clean PATH]... -- <command...>  any other GPU test (a run or a comparison)

--clean PATH removes PATH (a file or a directory) first, before the device is probed, also when the test is then
skipped: a run test names its output with it, so a skipped run leaves no file of an earlier run that the tests that
compare its output could take for this run's.

A skip is only for a machine without a CUDA device: when cppmpas_gpu reports another CUDA error (a driver too old for
the runtime, an initialisation error) the test fails. After the device check passed, a command that exits 77 fails
too (77 from the command is not a missing device). With the environment variable CPPMPAS_GPU_REQUIRED=1 (CMake
option CPPMPAS_GPU_REQUIRED) a test that would be skipped fails instead.

What a "gpu" test proves: cppmpas_gpu runs the dycore, the physics drivers and schemes, the radiation and the
diagnostics as kernels on the device, but some setup and bookkeeping code is host code in that build too, so those
parts are proven on the CPU only: physics_core_init (init_dirs_forphys, the accumulators and the land-use and soil
initialisation: landuse_init_forMPAS, init_lsm's soil-water loop), physics_update_surface / monthly_interp_to_date
(when the greenness alarm rings), radconst (the solar declination and constant, scalars), the isobaric level arrays
t_iso_levels / z_iso_levels and the reset of refl10cm_1km_max in tools/cppmpas.cpp.
"""
import os
import shutil
import subprocess
import sys
from pathlib import Path

SKIP = 77
# what the CUDA runtime says when there is no device (cudaErrorNoDevice), and src/core/backend_cuda.cu when the
# device count is 0: the only reasons for a skip
NO_DEVICE = ("no CUDA-capable device is detected", "(no CUDA device)")


def device(exe):
    """(status, output) of `cppmpas_gpu device`"""
    r = subprocess.run([exe, "device"], capture_output=True, text=True)
    return r.returncode, r.stdout + r.stderr


def skip(why):
    print(f"SKIP: {why}")
    if os.environ.get("CPPMPAS_GPU_REQUIRED") == "1":
        print("CPPMPAS_GPU_REQUIRED=1: a GPU test that cannot run is a failure")
        return 1
    return SKIP


def main():
    args = sys.argv[1:]
    if len(args) < 2:
        print(__doc__)
        return 2
    exe, args = args[0], args[1:]
    while len(args) >= 2 and args[0] == "--clean":
        p = Path(args[1])
        if p.is_dir() and not p.is_symlink():
            shutil.rmtree(p)
        else:
            p.unlink(missing_ok=True)
        args = args[2:]
    if not args:
        print(__doc__)
        return 2
    mode = args[0]
    if exe == "NONE":
        return skip("this build has no CUDA backend (configure with -DCPPMPAS_ENABLE_CUDA=ON)")
    if mode == "no-device":
        env = dict(os.environ, CUDA_VISIBLE_DEVICES="")
        r = subprocess.run([exe, "device"], capture_output=True, text=True, env=env)
        out = r.stdout + r.stderr
        print(out)
        if r.returncode == SKIP and "SKIP: no usable CUDA device" in out and "kernels run" not in out:
            return 0
        print(f"with no visible device cppmpas_gpu exited {r.returncode}, expected the skip status {SKIP}")
        return 1
    rc, out = device(exe)
    print(out.strip())
    if rc == SKIP:
        if not any(x in out for x in NO_DEVICE):
            print("cppmpas_gpu found no usable device for another reason than a missing one (a CUDA error): FAIL")
            return 1
        return skip("no CUDA device")
    if rc != 0:
        print(f"cppmpas_gpu device failed with status {rc}")
        return 1
    if mode == "device":
        return 0 if "backend cuda:" in out and "kernels run on the GPU, sm_" in out else 1
    if mode != "--" or len(args) < 2:
        print(__doc__)
        return 2
    r = subprocess.run(args[1:])
    if r.returncode == SKIP:
        print(f"the command exited {SKIP} after the device check passed: a failure, not a skip")
        return 1
    return r.returncode


if __name__ == "__main__":
    raise SystemExit(main())
