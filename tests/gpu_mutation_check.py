#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Show that the GPU tests (label gpu, tests/CMakeLists.txt) can fail: apply one one-line bug at a time to the source
of the CUDA build, rebuild, run the GPU tests the mutant must break, require each of them to fail (a test that did not
run, was skipped, timed out or passed does not count), restore the source. The tests must pass on the unmutated
source first, and the tests that killed a mutant must pass again on the restored source. Run it on the machine with
the GPU, in a CUDA build (-DCPPMPAS_ENABLE_CUDA=ON) whose CPU fixtures (initdump_*, dycore_*) the vs_cpu tests use.

usage: gpu_mutation_check.py <source dir> <build dir> [i,j,...]
The mutants: GPU-only bugs (the launch wrappers, kernel_fail, the stack limit, the skip path, device math, FMA
contraction, the CUDA flag policy) and mutants of tests/mutation_check.py, found by their fragment, for the GPU tests
of the dycore paths the GPU-only bugs do not single out (JW's Smagorinsky viscosity and 4th-order filters, the
restart run).
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
CM = "CMakeLists.txt"
MW10 = [f"gpu_dycore_MW_step{s}_{f}" for s in (1, 10) for f in ("u", "theta_m")] + ["gpu_dycore_MW_vs_cpu"]
JW10 = [f"gpu_dycore_JW_step{s}_{f}" for s in (1, 10) for f in ("u", "theta_m")] + ["gpu_dycore_JW_vs_cpu"]
# (file, original fragment, mutated fragment, GPU tests that must fail)
GPU_MUTATIONS = [
    # a kernel skips the last element of its range
    (B, "  if (i <= long(last)) k(int(i));", "  if (i < long(last)) k(int(i));", ["gpu_device", "gpu_backend"]),
    # a whole-array kernel skips the last element (the garbage element of an MPAS array: only gpu_backend sees it)
    (B, "  if (i < n) k(i);", "  if (i + 1 < n) k(i);", ["gpu_backend"]),
    # the one-thread kernel does nothing
    (B, "template <class K> __global__ void kserial(K k) { k(); }",
     "template <class K> __global__ void kserial(K k) { (void)k; }", ["gpu_backend"] + MW10),
    # every launch reports a kernel failure
    (BC, "  if (g_err && g_err->set) {", "  if (g_err) {", ["gpu_device", "gpu_backend"]),
    # no device: exit 0 instead of the skip status
    (X, "      return gpu::kNoDeviceExit;", "      return 0;", ["gpu_skip_without_device"]),
    # device math: CUDA's pow in place of the libmvec lane on the device only (the host keeps glibc's)
    (L, "  return ::cppwrf::glibc::mvec::pow(arg(x), arg(y));",
     "#if defined(__CUDA_ARCH__)\n  return ::pow(x, y);\n#else\n  return ::cppwrf::glibc::mvec::pow(arg(x), arg(y));\n#endif",
     ["gpu_libm", "gpu_initdump_MW_exner", "gpu_initdump_JW_exner"] + MW10),
    # FMA contraction in the device code: nvcc's default --fmad=true (the strict flag left out)
    (CM, "$<$<COMPILE_LANGUAGE:CUDA>:--fmad=false;--ftz=false;", "$<$<COMPILE_LANGUAGE:CUDA>:--ftz=false;",
     ["cmake_cuda_strict_flags_last"] + MW10 + JW10),
    # kernel_fail does not end the thread: the statements after it run
    (B, '  detail::record_kernel_error(head, true, value, tail);\n  end_thread();',
     '  detail::record_kernel_error(head, true, value, tail);', ["gpu_backend"]),
    # a launch before anything initialised the runtime has no kernel error slot
    (B, "template <class Kern> void reserve_once(Kern* kernel) {\n  gpu::ensure_ready();",
     "template <class Kern> void reserve_once(Kern* kernel) {", ["gpu_backend"]),
    # the stack limit set to a kernel's frame even when it is smaller (lowered)
    (BC, "  if (a.localSizeBytes <= limit) {  // the limit is never lowered", "  if (false) {", ["gpu_backend"]),
    # the Release CUDA flags overwritten (the user's CMAKE_CUDA_FLAGS_RELEASE and CMake's -DNDEBUG dropped)
    (CM, "  enable_language(CUDA)\n", "  enable_language(CUDA)\n  set(CMAKE_CUDA_FLAGS_RELEASE \"-O3\")\n",
     ["cmake_cuda_strict_flags_last"]),
    # host fast-math through -Xcompiler not refused
    (CM, "|-Ofast|-f(fast-math|unsafe-math-optimizations|", "|-Ofast|-f(unsafe-math-optimizations|",
     ["cmake_rejects_cuda_6_CMAKE_CUDA_FLAGS", "cmake_rejects_cuda_6_CMAKE_CUDA_FLAGS_RELEASE"]),
]
# tests/mutation_check.py mutants by a piece of their original fragment, with the GPU tests they must break
CPU_MUTANTS = [
    # JW: smagorinsky_2d's eddy viscosity from u twice instead of u and v
    ("smagorinsky_2d(m, f.eddy_visc_horz, f.u2, d.v,", JW10),
    # JW: the del4 vorticity of u_dissipation_3d (DM:735-751)
    ("        const real edge_sign = m.invAreaTriangle(iVertex) * m.dcEdge(iEdge) * m.edgesOnVertex_sign(i, iVertex);",
     JW10),
    # the RBF reconstruction coefficients and the local vertical unit vectors (init block, not in the initdump set)
    ("    alpha = alpha / static_cast<real>(pointCount);", ["gpu_init_restart_MW", "gpu_init_restart_JW"]),
    ("      mpas_unit_vec_in_r3(v);\n      for (int d = 1; d <= 3; ++d) lvu(d, iCell) = v[d - 1];",
     ["gpu_init_restart_JW"]),
    # the restart run: the coupled fields of the restart stream recomputed instead of kept (mpas_atm_core.F:522)
    ("  if (!cfg.config_do_restart || (cfg.config_do_restart && cfg.config_do_DAcycling)) atm_in",
     ["gpu_dycore_restart0_MW_rho_zz", "gpu_dycore_restart0_MW_exner", "gpu_dycore_restart0_MW_ru"]),
]


def cpu_mutant(piece):
    hits = [m for m in mc.MUTATIONS if piece in m[1]]
    if len(hits) != 1:
        raise SystemExit(f"tests/mutation_check.py: {len(hits)} mutants with the fragment {piece!r}")
    return hits[0][:3]


def mutants():
    return [m for m in GPU_MUTATIONS] + [cpu_mutant(p) + (t,) for p, t in CPU_MUTANTS]


def build(b):
    return subprocess.run(["cmake", "--build", str(b), "-j4"], capture_output=True, text=True)


def main():
    if len(sys.argv) not in (3, 4):
        raise SystemExit(__doc__)
    src, b = Path(sys.argv[1]), Path(sys.argv[2])
    muts = mutants()
    if len(sys.argv) == 4:
        muts = [muts[int(i)] for i in sys.argv[3].split(",")]
    bad, covered = [], set()
    r0 = build(b)
    if r0.returncode != 0:
        raise SystemExit("the unmutated source does not build")
    needed = sorted({t for *_, ts in muts for t in ts})
    r = mc.ctest(b, mc.any_of(needed))
    broken = mc.not_passing(r.stdout, needed)
    for t in sorted(broken):
        bad.append(f"{t} does not pass on the unmutated source (status {mc.test_status(r.stdout, t)!r})")
    for path, old, new, tests in muts:
        p = src / path
        orig = p.read_text()
        if orig.count(old) != 1:
            bad.append(f"{path}: fragment found {orig.count(old)} times: {old!r}")
            continue
        try:
            p.write_text(orig.replace(old, new))
            rb = build(b)
            if rb.returncode != 0:
                bad.append(f"{path}: mutant does not build: {new.strip()[:80]}")
                print(rb.stdout[-1500:], rb.stderr[-1500:])
                continue
            r = mc.ctest(b, mc.any_of(tests))
            for t in tests:
                status = mc.test_status(r.stdout, t)
                k = mc.is_kill(r.returncode, status, t, broken)
                print(f"{'KILLED' if k else 'SURVIVED'} {t:40s} by {path}: {new.strip()[:70]!r}", flush=True)
                if k:
                    covered.add(t)
                else:
                    bad.append(f"{t} survived {path}: {new.strip()[:80]} (status {status!r})")
        finally:
            p.write_text(orig)
    if build(b).returncode != 0:
        bad.append("the restored source does not build")
    if covered:
        r = mc.ctest(b, mc.any_of(sorted(covered)))
        for t in sorted(mc.not_passing(r.stdout, covered)):
            bad.append(f"{t} does not pass on the restored source: its kills are not trusted")
    for x in bad:
        print("PROBLEM", x)
    print(f"{len(covered)} GPU tests killed by at least one mutant; {len(bad)} problems")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
