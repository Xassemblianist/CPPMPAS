# CPPMPAS

A C++20/CUDA port of MPAS-Atmosphere, built to be **bit-identical** to the Fortran model.

CPPMPAS is the sister project of [CPPWRF](https://github.com/Xassemblianist/CPPWRF), a C++20/CUDA port of WRF-ARW 4.8.0
whose CPU and GPU builds reproduce stock WRF's output at tolerance 0: FP64 builds match WRF real\*8 and the FP32
identity build matches WRF real\*4, in every variable and output time of the verified configurations
([CPPWRF 1.4.0](https://github.com/Xassemblianist/CPPWRF/releases/tag/v1.4.0): 531 CPU checks and 423 CUDA FP64
checks against WRF, all at tolerance 0). CPPMPAS applies the same method to MPAS-Atmosphere.

## Status

**Early work.** What exists, all on the CPU in FP64:

- The mesh and I/O layer: MPAS grid, init and restart files are read into arrays laid out as MPAS holds them in
  memory on one MPI rank (MPAS's local edge and vertex numbering, local connectivity, garbage elements), and written
  back in MPAS's file order. Checked on the reference's own init files (round trip identical, numbering against an
  independent re-implementation) for the JW, SC and MW cases.
- All of `atm_mpas_init_block` before the first step: the mesh coefficients (signs, inverses, deformation and
  advection coefficients, mesh scaling, damping coefficients, boundary masks, `mpas_initialize_vectors`, the RBF
  reconstruction coefficients), `atm_init_coupled_diagnostics`, and the initial `atm_compute_solve_diagnostics` and
  `mpas_reconstruct`, with ports of the exact glibc 2.41 libm / libmvec entries the reference executable calls.
  **Verified at tolerance 0** against the MPAS v8.4.2 FP64 references in all three cases: every one of the 57 fields
  the init block computes (cppmpas-data releases `mpas842-dbl-mw-003`, `-sc-002`, `-jw-002`, whose `initdump.nc` holds
  them at the initial time), and the restart and first-step fields of `mpas842-dbl-{mw,sc,jw}-001`.
- The dynamical core for the mountain-wave (MW) and JW baroclinic-wave configurations (`cppmpas run --steps N`):
  `atm_srk3` with the acoustic steps, divergence damping, the 2nd-order horizontal mixing with a fixed (MW) or
  Smagorinsky (JW) eddy viscosity, the 4th-order filters (JW), the full vertical mixing (MW), Coriolis on the sphere
  and the per-step diagnostics. **Verified at tolerance 0** against `mpas842-dbl-mw-001` (steps 1-10 in every field
  of the step files, the restart at 30:06, the output at 30 and 60 min) and `mpas842-dbl-jw-001` (steps 1-10 in
  every field, the output at 24 h and the restart at 24:07:30).
- The supercell configuration on top of it: monotonic scalar transport of qv, qc, qr (`atm_advance_scalars` /
  `atm_advance_scalars_mono`), vertical mixing of the perturbation from the initial state, and Kessler
  microphysics with MPAS's microphysics interface (the Kessler kernel is shared with CPPWRF). **Verified at
  tolerance 0** against `mpas842-dbl-sc-001`: steps 1-10 in every field of the step files, and the output at
  30 min and the restart at 30:03 (cloud and rain present) in every field the port holds.
- Restart runs (`config_do_restart = true`): `cppmpas run` reads an MPAS restart file as MPAS reads its restart
  stream and continues from it. **Verified at tolerance 0** on MW: restarted from the reference's 00:30:06 restart
  file, the state at 01:00:00 is identical to the reference forecast's output in every field the port computes.
- The CUDA FP64 build (`-DCPPMPAS_ENABLE_CUDA=ON`, executable `cppmpas_gpu`, CUDA 12.1 or later): the same sources,
  every loop over cells, edges and vertices of the init block, the solve diagnostics and the dycore (the JW
  Smagorinsky viscosity and 4th-order filters, the SC scalar transport, perturbation mixing and Kessler microphysics
  included) a kernel with one thread per element (MPAS's neighbour and level order inside the thread), glibc's libm / libmvec ports on the device, `--fmad=false`. **Verified at
  tolerance 0** on an NVIDIA L4 (sm_89, CUDA 13.0, driver 580.82.07; `ctest -L gpu`, 706 GPU tests), against the
  references and against the CPU build: every initdump field of MW, SC and JW and the init block's restart-stream
  fields; MW steps 1-10, the 30:06 restart, the 30 and 60 min output and the restart run to 01:00:00; JW steps 1-10,
  the 24 h output and the 24:07:30 restart; SC steps 1-10, the 30-min output and the 30:03 restart. Without a CUDA build or device the GPU tests are reported as skipped.
- The mesoscale_reference physics on the GPU (`cppmpas_gpu run` with RG's namelists): the physics driver,
  MPAS_to_physics, the revised MM5 surface layer, Noah with the sea-ice driver and the surface diagnostics, YSU, GWDO,
  New Tiedtke, WSM6 with the radar reflectivity and relhum, the RRTMG LW / SW radiation with the cloudiness and the
  ozone climatology, physics_get_tend and the isobaric and convective diagnostics, every loop a kernel with one
  thread per column (MPAS's level order inside the thread, the reference's libmvec lanes). **Verified at tolerance
  0** on an NVIDIA RTX PRO 6000 Blackwell (sm_120, nvcc 13.0.88 with g++ 14 as the host compiler, driver
  580.82.07): RG steps 1-10, the 01:10 restart and the 1 h / 2 h histories against `mpas842-dbl-rg-004`, the
  rg-004 forecast's five stream files, and the 24 h forecast of `mpas842-dbl-rg24-001` (all 19 history, diag and
  restart files, 2168 variables); MW, JW and SC still pass on the GPU. The physics init (table reading, landuse,
  lsminit), the greenness update (physics_update_surface) and the clock and alarms run on the host, on the same
  managed arrays.
- Every other configuration (LES, LBCs, other physics, a Kessler restart, ...) is
  a hard error until it is ported and verified.
- The namelist: every option of MPAS-Atmosphere's Registry is known; a value that no verified reference runs with is a
  hard error naming the option.
- The MPAS reference itself is reproducible across machines: the MW reference built and run on an AMD EPYC and on an
  Intel Xeon runner is identical in every output file.

Build: `cmake -S . -B build && cmake --build build` (C++20, NetCDF-C; FP64). Tests need
`-DCPPMPAS_CASES_DIR=<extracted NCAR test-case tarballs>
-DCPPMPAS_MPAS_REGISTRY=<MPAS v8.4.2 src/core_atmosphere/Registry.xml> -DCPPMPAS_REFS_DIR=<directory with the
extracted mpas842-dbl-{mw,sc,jw}-001, mw-003, sc-002 and jw-002 releases>`; they compare with the in-tree `tools/mpas_compare.py`, whose own
tests run as `compare_unit` and `compare_mutants`. `unit_libm` needs glibc 2.41 on an x86-64 CPU with FMA and AVX2.

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
