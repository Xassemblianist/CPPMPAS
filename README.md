# CPPMPAS

A C++20/CUDA port of MPAS-Atmosphere, built to be **bit-identical** to the Fortran model.

CPPMPAS is the sister project of [CPPWRF](https://github.com/Xassemblianist/CPPWRF), a C++20/CUDA port of WRF-ARW 4.8.0
whose CPU and GPU builds reproduce stock WRF's output at tolerance 0: FP64 builds match WRF real\*8 and the FP32
identity build matches WRF real\*4, in every variable and output time of the verified configurations
([CPPWRF 1.4.0](https://github.com/Xassemblianist/CPPWRF/releases/tag/v1.4.0): 531 CPU checks and 423 CUDA FP64
checks against WRF, all at tolerance 0). CPPMPAS applies the same method to MPAS-Atmosphere.

## Status

**Planning. There is no model code here yet.** This repository records the plan; every capability will be listed in
this file only once it is verified against MPAS output.

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
