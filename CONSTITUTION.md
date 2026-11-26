# CPPMPAS constitution

The rules this project is built on. They come from CPPWRF, where they produced a full C++20/CUDA port of WRF-ARW
4.8.0 that matches WRF bit for bit on CPUs and GPUs. Changing a rule here needs the author's decision and a commit
that says why.

## 1. The goal

1. **Bit identity is the definition of correct.** The reference is MPAS-Atmosphere at one pinned release (v8.4.2),
   built with GNU Fortran in a pinned environment (Debian 13: gfortran 14.2, glibc 2.41), run in one process.
   CPPMPAS's FP64 CPU build must reproduce that reference in every compared field and output time at tolerance 0.
2. **A 1-ulp difference is a bug until it is explained.** Rounding differences compound in a nonlinear model and hide
   real errors. "Close enough" is never accepted for a build that claims identity. The explanation lands in the
   harness as a test.
3. **The GPU gives the same bits.** CUDA FP64 matches the same reference; GPU math uses device ports of the exact
   glibc routines the reference calls, no fused multiply-adds where the reference has none, and the reference's
   summation order.
4. **FP32 identity follows the reference's real\*4 build** once FP64 is done; a faster FP32 build that changes
   results is opt-in and is reported by its measured distance, never as identical.
5. **Ideas that change results belong to a separate fork.** CPPMPAS reproduces MPAS; it does not improve it.
6. **The reference is a recipe, not a version number.** It is pinned completely: compiler and version, the make
   target and its flags (`PRECISION=double`, optimisation and `-march` choices), the I/O layer (SMIOL or PIO) and
   its version, the MPI library and version (run on one process), the mesh and graph-partition files with their
   checksums, and the initial-condition generation (`init_atmosphere`) with its namelists and input data. The
   recipe is recorded with every reference. Any change to it makes a new reference, and comparisons against the
   old one are void until the references are regenerated.
7. **MPAS's own GPU (OpenACC) build is not a reference.** CPPMPAS is compared for correctness only with the pinned
   CPU real\*8 build (and later the real\*4 build). OpenACC output may be used to compare performance, never to
   judge correctness.

## 2. How porting is done

1. **Statement by statement.** Each C++ routine names the Fortran file and line it reproduces. Operation order,
   parenthesisation, constant folding, loop bounds and branches stay as the compiler built them, including
   integer powers, vectorised loops and the vector-math lanes (libmvec) the reference calls.
2. **One source for CPU and GPU.** Kernels are written once and run on the CPU (OpenMP) and the GPU (CUDA).
3. **Reductions keep the reference's order.** On the unstructured mesh every sum over the edges or vertices of a
   cell, and every global sum, follows MPAS's order on one process; decomposition (tiles, MPI ranks, GPU blocks)
   must not change it.
4. **No silent options.** Every namelist option or value without a verified port is rejected with a clear error,
   in the model and in any translator. Nothing runs "approximately".
5. **Shared physics is re-verified, never assumed.** Where MPAS calls the same scheme as WRF, CPPWRF's port is a
   starting point only: constants (g = 9.80616 m s⁻² in MPAS), wrappers and pinned physics revisions differ, so the
   scheme is checked against MPAS's own output.
6. **Code shared with CPPWRF never diverges silently.** Every scheme or routine taken from CPPWRF records its
   source (CPPWRF commit and file) in a header comment and in `SHARED_WITH_CPPWRF.md`. MPAS-specific differences
   (constants such as g, wrapper and revision differences) enter as parameters, never as edits of the shared
   logic. A bug fix found in either project is reported to the other and applied and verified there, and the
   ledger records both commits. A common library may replace the copies later; until then this rule holds.

## 3. How it is proven

1. **Every claim has a check.** A capability is listed as supported only after it passes against MPAS output at
   tolerance 0 in the automated suite. Unproven ports may exist in the code only behind an explicit
   "unverified" switch that is off by default.
2. **A test counts only if it can fail.** Each new check is shown to fail on a deliberate one-line bug before it is
   trusted. A check that did not run is a failure, not a pass. One deliberate bug per new check, shown once; full
   mutation sweeps only at release.
3. **Find differences by bisection.** By configuration (one option at a time), by time (output every step, first
   differing step and field), then by dumps of internal arrays in both codes at the same point.
4. **Proofs run on reproducible infrastructure.** Reference runs and the CPU suite run in the pinned environment on
   automated CI; GPU proofs run on recorded hardware. Every release states the commit, hardware and counts of its
   proofs.
5. **Honest reporting.** Numbers are reported with the command and commit that produced them; what was not verified
   is said plainly. Documentation never claims more than the checks show.

## 4. How the project is run

1. **Only `main` is published.** Work happens on branches that are merged after their proofs pass; history on
   `main` is never rewritten.
2. **Releases are tagged only after a full CPU suite and a GPU proof on the release commit** (or on a commit whose
   model code is identical, stated in the release notes).
3. **Respect other work.** Public text describes CPPMPAS by its own measured results and never disparages other
   ports or the MPAS team's work.
4. **Licences are kept.** CPPMPAS's own code is GPL-3.0-or-later; files derived from MPAS or WRF keep their
   upstream notices.

## 5. Order of work

1. Reference and harness: pinned MPAS build, reference runs of standard idealised and real cases, the comparator.
2. Mesh and I/O: MPAS meshes, `init.nc` / `lbc` files and partitions read directly; output MPAS tools can read.
3. Dynamical core on the CPU in FP64, then CUDA FP64.
4. Physics suites (`mesoscale_reference` first, then `convection_permitting`).
5. FP32 identity, multi-GPU, performance.
