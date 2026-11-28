#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""The top-level CMakeLists.txt refuses every flag that changes floating-point results (cppmpas_check_fp_flags):
configure the project with each such flag, in each place a flag can come from, and require that the configure fails
with the refusal message naming that flag. Look-alike flags that do not change results (-fno-fast-math,
-fsigned-zeros, -mno-fma, --fmad=false, ...) must configure, so the check matches whole flags only.

usage: check_fp_flags.py <cmake> <source dir> <work dir> <generator> <C++ compiler>

Each configure runs in its own directory under <work dir>, which is deleted first and at the end. The CUDA flags are
checked also when the build has no CUDA (CMAKE_CUDA_FLAGS is checked whether CUDA is enabled or not), so this test
needs no nvcc.
"""
import os
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

MSG = "changes the floating-point results (CPPMPAS needs strict IEEE)"

# (cache variable, value, the flag the refusal must name)
REFUSED = [(v, f, f) for v in ("CMAKE_CXX_FLAGS",) for f in (
    "-fno-signed-zeros", "-mfma", "-mfma4", "-ffast-math", "-funsafe-math-optimizations", "-ffp-contract=fast",
    "-ffp-contract=on", "-ffp-contract=fast-honor-pragmas", "-freciprocal-math", "-fassociative-math", "-Ofast",
    "-ffinite-math-only", "-fcx-limited-range", "-march=native", "-march=x86-64-v3", "-mfpmath=387",
    "-mfpmath=sse+387", "-m32", "-mx32", "-mno-sse2", "-fexcess-precision=fast")]
REFUSED += [
    # inside a list of harmless flags
    ("CMAKE_CXX_FLAGS", "-O2 -g -fno-signed-zeros -Wall", "-fno-signed-zeros"),
    ("CMAKE_CXX_FLAGS", "-Wall\t-mfma", "-mfma"),
    # every configuration, also one that is not built and one CMake does not know
    ("CMAKE_CXX_FLAGS_RELEASE", "-ffast-math", "-ffast-math"),
    ("CMAKE_CXX_FLAGS_DEBUG", "-mfma", "-mfma"),
    ("CMAKE_CXX_FLAGS_RELWITHDEBINFO", "-freciprocal-math", "-freciprocal-math"),
    ("CMAKE_CXX_FLAGS_MINSIZEREL", "-fassociative-math", "-fassociative-math"),
    ("CMAKE_CXX_FLAGS_PROFILE", "-funsafe-math-optimizations", "-funsafe-math-optimizations"),
    # nvcc, one dash or two, '=' or a blank
    ("CMAKE_CUDA_FLAGS", "--use_fast_math", "--use_fast_math"),
    ("CMAKE_CUDA_FLAGS", "-use_fast_math", "-use_fast_math"),
    ("CMAKE_CUDA_FLAGS", "--fmad=true", "--fmad=true"),
    ("CMAKE_CUDA_FLAGS", "-fmad=true", "-fmad=true"),
    ("CMAKE_CUDA_FLAGS", "--fmad true", "--fmad true"),
    ("CMAKE_CUDA_FLAGS", "--fmad=1", "--fmad=1"),
    ("CMAKE_CUDA_FLAGS", "--ftz=true", "--ftz=true"),
    ("CMAKE_CUDA_FLAGS", "--prec-div=false", "--prec-div=false"),
    ("CMAKE_CUDA_FLAGS", "--prec-sqrt=false", "--prec-sqrt=false"),
    ("CMAKE_CUDA_FLAGS_RELEASE", "--use_fast_math", "--use_fast_math"),
    ("CMAKE_CUDA_FLAGS_DEBUG", "--fmad=true", "--fmad=true"),
    # host flags through nvcc, and ptxas's
    ("CMAKE_CUDA_FLAGS", "-Xcompiler=-ffast-math", "-ffast-math"),
    ("CMAKE_CUDA_FLAGS", "-Xcompiler -fno-signed-zeros", "-fno-signed-zeros"),
    ("CMAKE_CUDA_FLAGS", "--compiler-options=-O3,-mfma", "-mfma"),
    ("CMAKE_CUDA_FLAGS", '-Xcompiler="-ffp-contract=fast"', "-ffp-contract=fast"),
    ("CMAKE_CUDA_FLAGS", "-Xptxas --fmad=true", "--fmad=true"),
    # the link (-ffast-math links crtfastmath.o: flush to zero)
    ("CMAKE_EXE_LINKER_FLAGS", "-ffast-math", "-ffast-math"),
    ("CMAKE_SHARED_LINKER_FLAGS_RELEASE", "-Ofast", "-Ofast"),
]
# (environment variable, value, flag)
REFUSED_ENV = [("NVCC_APPEND_FLAGS", "--use_fast_math", "--use_fast_math"),
               ("NVCC_PREPEND_FLAGS", "--fmad=true", "--fmad=true"),
               ("CXXFLAGS", "-ffast-math", "-ffast-math")]
# CMake code run right after project() (CMAKE_PROJECT_INCLUDE): a directory, a target or a source file with a flag
REFUSED_INCLUDE = [
    ("add_compile_options(-fno-signed-zeros)", "-fno-signed-zeros"),
    ("add_link_options(-ffast-math)", "-ffast-math"),
    ("add_library(fp_probe INTERFACE)\ntarget_compile_options(fp_probe INTERFACE -mfma)", "-mfma"),
    ("add_library(fp_probe INTERFACE)\ntarget_link_options(fp_probe INTERFACE -Ofast)", "-Ofast"),
    ("add_library(fp_probe INTERFACE)\nset_property(TARGET fp_probe PROPERTY INTERFACE_COMPILE_OPTIONS "
     "\"$<$<COMPILE_LANGUAGE:CXX>:-ffast-math>\")", "-ffast-math"),
    ("set_source_files_properties(tools/cppmpas.cpp PROPERTIES COMPILE_OPTIONS -freciprocal-math)",
     "-freciprocal-math"),
    ("set_source_files_properties(src/io/mpas_writer.cpp PROPERTIES COMPILE_FLAGS \"-O3 -fassociative-math\")",
     "-fassociative-math"),
]
# flags that look like refused ones but do not change results: the whole project must configure with them. A refused
# flag at the end of a longer word (-I/opt/build-m32, -DTAG=x-mfma) and one at the start of a longer flag (-mno-sse3,
# -mno-sse2avx: -mno-sse and -mno-sse2 are refused) show that whole flags are matched, not substrings
ACCEPTED = [
    ("CMAKE_CXX_FLAGS", "-fno-fast-math -fsigned-zeros -mno-fma -ffp-contract=off -fno-reciprocal-math "
     "-fno-associative-math -fno-unsafe-math-optimizations -fno-finite-math-only -fexcess-precision=standard "
     "-mtune=generic -O3 -I/opt/build-m32 -I/opt/fftw-ffast-math -DTAG=x-mfma -mno-sse3 -mno-sse4 -mno-sse2avx"),
    ("CMAKE_CUDA_FLAGS", "--fmad=false -fmad=false --ftz=false --prec-div=true --prec-sqrt=true "
     "-Xcompiler=-fno-fast-math,-ffp-contract=off --compiler-options=-fsigned-zeros"),
]


def configure(cmake, src, wd, gen, cxx, defs, env=None, include=None):
    if wd.exists():
        shutil.rmtree(wd)
    wd.mkdir(parents=True)
    args = [cmake, "-S", src, "-B", str(wd / "b"), "-G", gen, f"-DCMAKE_CXX_COMPILER={cxx}"]
    if include is not None:
        (wd / "inject.cmake").write_text(include + "\n")
        args.append(f"-DCMAKE_PROJECT_INCLUDE={wd / 'inject.cmake'}")
    args += [f"-D{k}={v}" for k, v in defs]
    e = dict(os.environ)
    for k in ("CXXFLAGS", "CFLAGS", "CUDAFLAGS", "LDFLAGS", "NVCC_APPEND_FLAGS", "NVCC_PREPEND_FLAGS"):
        e.pop(k, None)
    e.update(env or {})
    r = subprocess.run(args, capture_output=True, text=True, env=e)
    shutil.rmtree(wd)
    return r.returncode, r.stdout + r.stderr


def main():
    if len(sys.argv) != 6:
        print(__doc__)
        return 2
    cmake, src, work, gen, cxx = sys.argv[1:]
    work = Path(work)
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    jobs = []  # (description, expected flag or None, kwargs)
    for v, val, flag in REFUSED:
        jobs.append((f"-D{v}='{val}'", flag, dict(defs=[(v, val)])))
    for k, val, flag in REFUSED_ENV:
        jobs.append((f"environment {k}='{val}'", flag, dict(defs=[], env={k: val})))
    for code, flag in REFUSED_INCLUDE:
        jobs.append((f"CMAKE_PROJECT_INCLUDE: {code!r}", flag, dict(defs=[], include=code)))
    jobs.append(("look-alike flags " + "; ".join(f"{k}='{v}'" for k, v in ACCEPTED), None, dict(defs=ACCEPTED)))

    def run(i):
        desc, flag, kw = jobs[i]
        rc, out = configure(cmake, src, work / f"c{i}", gen, cxx, **kw)
        if flag is None:
            ok = rc == 0 and "changes the floating-point" not in out
            return ok, desc, "configured" if ok else f"rc {rc}:\n{out[-3000:]}"
        # CMake wraps a long message at blanks: compared with the blanks normalised
        want = " ".join(f"'{flag}' {MSG}".split())
        ok = rc != 0 and want in " ".join(out.split())
        return ok, desc, f"refused '{flag}'" if ok else f"rc {rc}, no \"{want}\" in:\n{out[-3000:]}"

    with ThreadPoolExecutor(max_workers=4) as ex:
        results = list(ex.map(run, range(len(jobs))))
    bad = 0
    for ok, desc, what in results:
        print(("OK   " if ok else "FAIL ") + desc + ": " + what)
        bad += not ok
    shutil.rmtree(work, ignore_errors=True)
    if bad:
        print(f"{bad} of {len(jobs)} configures did not behave as required")
        return 1
    print(f"OK {len(jobs)} configures: {len(jobs) - 1} refused with the flag named, the look-alikes accepted")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
