#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Run one GPU test: skip it (exit status 77, CTest's SKIP_RETURN_CODE, reported as "Skipped", never "Passed") when
there is no CUDA build (NONE) or no usable CUDA device, else run the command and return its status.

usage: gpu_test.py <cppmpas_gpu | NONE> device               the device check: the kernels must run on the GPU
       gpu_test.py <cppmpas_gpu | NONE> no-device            with no visible device, cppmpas_gpu must skip (77)
       gpu_test.py <cppmpas_gpu | NONE> -- <command...>      any other GPU test (a run or a comparison of its output)
"""
import os
import subprocess
import sys

SKIP = 77


def device(exe):
    """(status, output) of `cppmpas_gpu device`"""
    r = subprocess.run([exe, "device"], capture_output=True, text=True)
    return r.returncode, r.stdout + r.stderr


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        return 2
    exe, mode = sys.argv[1], sys.argv[2]
    if exe == "NONE":
        print("SKIP: this build has no CUDA backend (configure with -DCPPMPAS_ENABLE_CUDA=ON)")
        return SKIP
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
        print("SKIP: no usable CUDA device")
        return SKIP
    if rc != 0:
        print(f"cppmpas_gpu device failed with status {rc}")
        return 1
    if mode == "device":
        return 0 if "backend cuda:" in out and "kernels run on the GPU, sm_" in out else 1
    if mode != "--" or len(sys.argv) < 4:
        print(__doc__)
        return 2
    r = subprocess.run(sys.argv[3:])
    # a command that itself reports a skip is a skip; any other status is the command's
    return r.returncode


if __name__ == "__main__":
    raise SystemExit(main())
