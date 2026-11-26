#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Show that every GPU test (label gpu, tests/CMakeLists.txt) can fail: apply one one-line bug at a time to the source
of the CUDA build, rebuild cppmpas_gpu and test_gpu_libm, run every GPU test, record the tests that failed (a test
that did not run, was skipped or passed does not count), restore the source. At the end every GPU test must have
been killed by at least one mutant. Run it through the gate on the machine with the GPU, in a CUDA build
(-DCPPMPAS_ENABLE_CUDA=ON) whose CPU fixtures (initdump_*, dycore_MW) the vs_cpu tests use.

usage: gpu_mutation_check.py <source dir> <build dir> [--long] [i,j,...]
  --long  include the label-long GPU tests (301 and 600 steps) in every mutant's run
The mutants: GPU-only bugs (the launch, the skip path, device math, FMA contraction) and the CPU harness's mutants
(tests/mutation_check.py, by index) of the arrays that no GPU-only bug reaches (xtime, Time, the zone masks, ...).
"""
import importlib.util
import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("mutation_check", HERE / "mutation_check.py")
mc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mc)

B, BC, X, L = "src/core/backend.hpp", "src/core/backend_cuda.cu", "tools/cppmpas.cpp", "src/core/libm.hpp"
GPU_MUTATIONS = [
    # a kernel skips its last element (one thread per element, CUDA only)
    (B, "  if (i <= long(last)) k(int(i));", "  if (i < long(last)) k(int(i));"),
    # a whole-array kernel skips its last element
    (B, "  if (i < n) k(i);", "  if (i + 1 < n) k(i);"),
    # the one-thread kernel does nothing
    (B, "template <class K> __global__ void kserial(K k) { k(); }", "template <class K> __global__ void kserial(K k) { (void)k; }"),
    # every launch reports a kernel failure
    (BC, "  if (g_err && g_err->set) {", "  if (g_err) {"),
    # no device: exit 0 instead of the skip status
    (X, "      return gpu::kNoDeviceExit;", "      return 0;"),
    # device math: CUDA's pow in place of the libmvec lane (the host keeps glibc's pow)
    (L, "  return ::cppwrf::glibc::mvec::pow(arg(x), arg(y));", "  return ::pow(x, y);"),
    # FMA contraction in the device code
    ("CMakeLists.txt", '" --fmad=false --ftz=false', '" --fmad=true --ftz=false'),
]
# tests/mutation_check.py mutants (by index) for what the GPU-only bugs do not reach: config_start_time (xtime of the
# steps), the epssm coefficients (serial kernel), the zone masks and nearestRelaxationCell (whole-array kernels), the
# xtime writer, the Time of the steps
CPU_MUTANTS = [84, 158, 160, 161, 162, 163, 178, 321]


def gpu_tests(build, long_too):
    args = ["ctest", "--test-dir", str(build), "-N", "-L", "gpu"] + ([] if long_too else ["-LE", "long"])
    out = subprocess.run(args, capture_output=True, text=True, check=True).stdout
    return re.findall(r"^\s*Test\s+#\d+: (\S+)$", out, re.M)


def main():
    argv = [a for a in sys.argv[1:] if a != "--long"]
    long_too = "--long" in sys.argv
    src, build = Path(argv[0]), Path(argv[1])
    muts = GPU_MUTATIONS + [mc.MUTATIONS[i][:3] for i in CPU_MUTANTS]
    if len(argv) > 2:
        muts = [muts[int(i)] for i in argv[2].split(",")]
    tests = gpu_tests(build, long_too)
    print(f"{len(tests)} GPU tests")
    covered, bad = set(), []
    targets = ["--target", "cppmpas_gpu", "--target", "test_gpu_libm", "--target", "cppmpas"]
    for path, old, new in muts:
        p = src / path
        orig = p.read_text()
        if orig.count(old) != 1:
            bad.append(f"{path}: fragment found {orig.count(old)} times: {old!r}")
            continue
        try:
            p.write_text(orig.replace(old, new))
            b = subprocess.run(["cmake", "--build", str(build), "-j14"] + targets, capture_output=True, text=True)
            if b.returncode != 0:
                bad.append(f"{path}: mutant does not build: {new.strip()}")
                print(b.stdout[-2000:], b.stderr[-2000:])
                continue
            args = ["ctest", "--test-dir", str(build), "-j8", "-L", "gpu"] + ([] if long_too else ["-LE", "long"])
            r = subprocess.run(args, capture_output=True, text=True)
            # the CPU fixtures the vs_cpu tests need run as fixtures (-L selects their dependents' setups too)
            killed = [t for t in tests if mc.killed(mc.test_status(r.stdout, t)) and
                      "Skipped" not in (mc.test_status(r.stdout, t) or "")]
            covered.update(killed)
            print(f"MUTANT {path}: {new.strip()[:80]}\n  killed {len(killed)} of {len(tests)}: "
                  + " ".join(sorted(killed)[:12]) + (" ..." if len(killed) > 12 else ""))
            if not killed:
                bad.append(f"{path}: {new.strip()} killed no GPU test")
        finally:
            p.write_text(orig)
    subprocess.run(["cmake", "--build", str(build), "-j14"] + targets, capture_output=True, text=True, check=True)
    survivors = [t for t in tests if t not in covered]
    for t in survivors:
        bad.append(f"{t} survived every mutant")
    for b in bad:
        print("PROBLEM", b)
    print(f"{len(covered)} of {len(tests)} GPU tests killed by at least one mutant")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
