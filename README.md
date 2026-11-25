# CPPMPAS

A C++20/CUDA port of MPAS-Atmosphere, built to be **bit-identical** to the Fortran model.

CPPMPAS is the sister project of [CPPWRF](https://github.com/Xassemblianist/CPPWRF), a C++20/CUDA port of WRF-ARW 4.8.0
whose CPU and GPU builds reproduce stock WRF's output at tolerance 0: FP64 builds match WRF real\*8 and the FP32
identity build matches WRF real\*4, in every variable and output time of the verified configurations
([CPPWRF 1.4.0](https://github.com/Xassemblianist/CPPWRF/releases/tag/v1.4.0): 531 CPU checks and 423 CUDA FP64
checks against WRF, all at tolerance 0). CPPMPAS applies the same method to MPAS-Atmosphere.

## Status

**Early work. There is no time step yet.** What exists, all on the CPU in FP64:

- The mesh and I/O layer: MPAS grid, init and restart files are read into arrays laid out as MPAS holds them in
  memory on one MPI rank (MPAS's local edge and vertex numbering, local connectivity, garbage elements), and written
  back in MPAS's file order. Checked on the reference's own init files (round trip identical, numbering against an
  independent re-implementation) for the JW, SC and MW cases.
- The init-time computations of `atm_mpas_init_block` before the first step (signs, inverses, deformation and
  advection coefficients, mesh scaling, damping coefficients, boundary masks, `mpas_initialize_vectors`, the RBF
  reconstruction coefficients, `atm_init_coupled_diagnostics`), with ports of the exact glibc 2.41 libm / libmvec
  entries the reference executable calls. **Verified at tolerance 0** against the MPAS v8.4.2 FP64 references
  (cppmpas-data releases `mpas842-dbl-{jw,sc,mw}-001`) for every array MPAS writes: `edgeNormalVectors`,
  `localVerticalUnitVectors`, `cellTangentPlane`, `coeffs_reconstruct`, `dss`, `exner_base`, `rtheta_base`,
  `pressure_base` (restart file) and `rho_zz`, `theta_m`, `pressure_p` (initial record of the every-step history), in
  all three cases. The other computed arrays are not written by MPAS and are **unverified** until the reference
  writes them (an `initdump` stream is prepared for the reference pipeline).

Build: `cmake -S . -B build && cmake --build build` (C++20, NetCDF-C; FP64). Tests need
`-DCPPMPAS_CASES_DIR=<extracted NCAR test-case tarballs> -DCPPMPAS_MPAS_COMPARE=<tools/mpas_compare.py>
-DCPPMPAS_MPAS_REGISTRY=<MPAS v8.4.2 src/core_atmosphere/Registry.xml> -DCPPMPAS_REFS_DIR=<directory with the
extracted mpas842-dbl-{mw,sc,jw}-001 releases>`; `unit_libm` needs glibc 2.41 on an x86-64 CPU with FMA and AVX2.

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
