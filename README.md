# CPPMPAS

A C++20/CUDA port of MPAS-Atmosphere, built to be **bit-identical** to the Fortran model.

CPPMPAS is the sister project of [CPPWRF](https://github.com/Xassemblianist/CPPWRF), a C++20/CUDA port of WRF-ARW 4.8.0
whose CPU and GPU builds reproduce stock WRF's output at tolerance 0: FP64 builds match WRF real\*8 and the FP32
identity build matches WRF real\*4, in every variable and output time of the verified configurations
([CPPWRF 1.4.0](https://github.com/Xassemblianist/CPPWRF/releases/tag/v1.4.0): 531 CPU checks and 423 CUDA FP64
checks against WRF, all at tolerance 0). CPPMPAS applies the same method to MPAS-Atmosphere.

## Status

On `main` (`3ca63da`), the CPU and CUDA FP64 builds are bit-identical (tolerance 0) to MPAS-Atmosphere 8.4.2 built
with gfortran 14.2, `-O3 -fdefault-real-8`, on glibc 2.41, run on one rank with the Intel `libmvec` preload.

The preload is part of the reference. Unmodified MPAS itself gives different results on Intel and AMD CPUs, because
glibc's vector `log` uses `rcpps`, whose result differs between the two vendors. The RG references
`mpas842-dbl-rg-001` and `-002` (AMD) differ from `-003` (Intel). With the preload, Intel and AMD runs are identical
(`rg-004` and `rg-006`). See [docs/REFERENCES.md](docs/REFERENCES.md).

- **RG24.** The 24 h real-data forecast on the x1.10242 mesh uses the full `mesoscale_reference` suite. Both
  builds match `mpas842-dbl-rg24-001` bit for bit across all 19 history, diagnostics and restart files and
  all 2168 variables. The comparison excludes the random `file_id` and `parent_id` global attributes.
- **GPU.** The RG24 run used a G4 NVIDIA RTX PRO 6000 Blackwell (`sm_120`, `nvcc` 13.0.88, CUDA runtime 13.1).
  Radiation ran on the device. Physics initialization, the greenness update, the clock and stream output ran
  on the host. The GPU run used `842125d`. Later changes through `3ca63da` were limited to the README and GPU
  mutation checks.
- **RG.** Both builds match `mpas842-dbl-rg-004` for the checked fields at steps 0 through 10 and all
  702 variables in the 5 forecast stream files.
- **Idealized cases.** MW (mountain wave), JW (Jablonowski–Williamson baroclinic wave) and SC (supercell with
  Kessler microphysics) match on the CPU and GPU for the 57 initdump fields per case, the tested initial
  restart-stream fields, and the checked fields at steps 0 through 10.
  - MW matches the checked output fields at 30 and 60 min and the restart fields at 00:30:06. A run that
    reads that restart matches the checked output fields at 01:00:00.
  - JW matches the checked output fields at 24 h and the restart fields at 24:07:30.
  - SC matches the checked output fields at 30 min and the restart fields at 00:30:03.
- **CPU suites.** With the Intel `libmvec` preload, `ctest -LE long` reports 2425/2425 and
  `ctest -L long` reports 306/306. The CPU build skips the GPU tests.

Open work includes restart runs with physics, the `convection_permitting` suite, variable-resolution meshes,
MPI and multi-GPU runs, and performance work.

## Build and tests

Build with `cmake -S . -B build && cmake --build build` (C++20, NetCDF-C, FP64). Configure with
`-DCPPMPAS_ENABLE_CUDA=ON` to build `cppmpas_gpu`.

For tests, set `CPPMPAS_CASES_DIR` to the extracted NCAR test cases, `CPPMPAS_MPAS_REGISTRY` to MPAS's
`src/core_atmosphere/Registry.xml`, and `CPPMPAS_REFS_DIR` to the directory holding the extracted
[reference releases](docs/REFERENCES.md). The tests compare with the in-tree `tools/mpas_compare.py`, whose own
tests run as `compare_unit` and `compare_mutants`. `unit_libm` needs glibc 2.41 on an x86-64 CPU with FMA and AVX2.

### CPU suite on GitHub runners

The workflow [`cpu_suite`](.github/workflows/cpu_suite.yml) runs the CPU suite in the Debian 13 image of the
references (glibc 2.41, gcc and gfortran 14.2). It builds the Release CPU build and the Intel `libmvec` preload, and
fetches every input from cppmpas-data. Each input is checked against its SHA-256 in [ci/data.tsv](ci/data.tsv) and
kept in the Actions cache. The tests are split into parallel jobs ([ci/plan.py](ci/plan.py)), and ctest runs under the
preload:

```bash
gh workflow run cpu_suite.yml -f ref=<branch or commit> -f labels=all   # or short (no label long), or long
gh run list --workflow cpu_suite.yml -L 1                              # then: gh run watch <id>
```

The run fails if any test fails or does not run. The `summary` artifact lists the counts per label and every
failing test. This is a CPU build, so the GPU tests exit 77; they are counted as skipped, never as passed. Pushes to
`main` run `all`.

A local commit that is not pushed goes through [tools/ci/ci_submit.sh](tools/ci/ci_submit.sh):

```bash
tools/ci/ci_submit.sh [COMMIT] [--labels short|long|all] [--no-wait]
```

It uploads a git bundle of the commits `origin/main` lacks to the `jobs` release of cppmpas-data under a unique
name, starts `cpu_suite` with it (the build job checks its SHA-256 and checks the commit out,
[ci/fetch_bundle.sh](ci/fetch_bundle.sh)), waits, prints the summary and deletes the bundle whatever the outcome.

Only the workflow file (and `ci/fetch_bundle.sh`) comes from the workflow's ref (`--workflow-ref`, default `main`).
Everything else in `ci/` comes from the commit under test: the image, the build, the plan with `ci/data.tsv` and
`ci/costs.tsv`, the data checks, the test runner and the summary. So a commit that adds a reference release to
`ci/data.tsv`, or changes the planner, is tested as it is. The workflow file keeps its own checks: the image's
SHA-256, and a verdict that reads every JUnit result itself (each selected CPU test passed, each GPU test skipped, at
least one CPU test, every job succeeded).

## Goal

- The same answers as MPAS-Atmosphere, bit for bit: a pinned MPAS release built with GNU Fortran is the reference, and
  every port is checked against its output at tolerance 0, field by field and output time by output time.
- Fast on NVIDIA GPUs, from one GPU to several (MPI), without giving up identity.
- One code base for CPU and GPU, readable C++ with the Fortran routine it reproduces named at each port.

## Plan

1. **Feasibility map.** Which MPAS-Atmosphere physics schemes are the WRF schemes CPPWRF already reproduces, which
   differ, and what the unstructured Voronoi-mesh dynamical core (data layout, halos, the split-explicit solver) needs.
   Measure the reference model's CPU and GPU (OpenACC) performance on standard meshes.
2. **Reference and harness.** A reproducible build of the pinned MPAS release, reference runs, and the comparison
   harness (CPPWRF's exact comparator and mutation-tested checks, adapted to MPAS's mesh files and streams).
3. **Mesh and I/O.** Reading MPAS meshes and initial conditions, writing output MPAS tools read.
4. **Dynamical core.** The atmosphere solver on the CPU in FP64, identical to MPAS real\*8, then on the GPU.
5. **Physics.** The physics suites, reusing CPPWRF's verified ports where MPAS calls the same code.
6. **FP32 identity, multi-GPU, performance.**

## Feasibility (MPAS-Atmosphere v8.4.2 against CPPWRF 1.4.0)

A first map of the work, counted from the two source trees:

- **Physics, `mesoscale_reference` suite** (WSM6, new Tiedtke, YSU + GWDO, RRTMG LW/SW, revised MM5 surface layer,
  Noah): every scheme has a CPPWRF port. The RRTMG cores and the Noah land model match WRF 4.8.0 almost line for line;
  the MPAS wrappers, and the external MMM-physics revision MPAS pins, differ in places, so each port is re-verified
  against MPAS, not assumed.
- **Physics, `convection_permitting` suite** (Thompson, Grell-Freitas, MYNN): MPAS calls older or different code than
  WRF 4.8.0 (Thompson based on WRF 4.1.4, the monolithic Grell-Freitas, MMM-physics MYNN), so these need new ports.
- **Constants:** MPAS uses g = 9.80616 m s⁻² where WRF uses 9.81, so physics shared with CPPWRF is parameterised, not
  copied.
- **Dynamical core:** about 11,700 lines of Fortran (CPPWRF's WRF dynamical core is about 9,000 lines of C++/CUDA): a
  split-explicit solver on the unstructured Voronoi C-grid, whose sums over the edges of each cell must keep MPAS's order
  for bit identity.
- **GPU in the reference model:** MPAS v8.4.2 carries 1,331 lines of OpenACC, all in the dynamical core and the halo
  and vector-reconstruction code; its physics runs on the CPU.
- **First milestone:** read MPAS's own `init.nc` / `lbc` files and mesh partitions directly (as CPPWRF reads WRF's
  `wrfinput` / `wrfbdy`), port the dynamical core on the CPU in FP64, and compare with MPAS real\*8 at tolerance 0.

## Method (from CPPWRF)

- The reference is the Fortran model as built by its default toolchain, including vectorised loops and the exact
  libm and libmvec routines it calls; GPU builds use device ports of the same glibc functions.
- Any difference is a bug until it is explained, and the explanation lands in the harness as a test.
- Ideas that would change results are kept for a separate fork.

## License

GPL-3.0-or-later (LICENSE). Files derived from MPAS keep MPAS's own copyright and license notices (MPAS is distributed
under a BSD-style license), and files derived from WRF keep WRF's notices, as in CPPWRF.

## Author

Ömer Kaan Yılmaz
