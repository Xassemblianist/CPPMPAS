# RG24 timing

`rg24_timing.sh` times RG24 (mpas842-dbl-rg24-001: 24 h on x1.10242, 240 km, 55 levels, `mesoscale_reference`,
144 steps of 600 s, 48 radiation calls, 19 stream files) with CPPMPAS and with MPAS 8.4.2. It runs one model at a
time and compares every run's outputs with the reference at tolerance 0. `timed.py` measures one run. The header of
the script lists its inputs.

Nothing inside CPPMPAS is timed. The breakdown comes from the output files' modification times:
- "init" is the time to the first output (reading the init file, the setup, the 00:00 history and diag files).
- "per step" is (last output - first output) / 144. It includes the physics, the radiation and the output
  writes of the 3-hourly files and the 12:10 restart.

MPAS's own timers (log.atmosphere.0000.out, MPAS_NATIVE_TIMERS) give the same split for MPAS.

## First measurement (2026-11-27, main 3ca63da)

Machine: one Colab G4 VM.
- CPU: AMD EPYC 9B45, 24 cores / 48 threads, 176 GB.
- GPU: NVIDIA RTX PRO 6000 Blackwell Server Edition (sm_120), driver 580.82.07.

CPPMPAS build: Release, conda g++ 14.4 (`-O3`), nvcc 13.0.88 (`--fmad=false`), and the 13.1 libcudart that the VM's
loader picks.

MPAS build: as cppmpas-data ci/build.sh builds it, with conda-forge's gfortran 14.2.0, Open MPI 5.0.11, netCDF 4.10.1,
netCDF-Fortran 4.6.4 and PnetCDF 1.14.1. The flags are `-O3 ... -fdefault-real-8 -fdefault-double-8`, with SMIOL. It
runs with the Intel libmvec preload.

Wall times, 3 runs each, in seconds:

| configuration | runs | median | init | per step | outputs vs rg24-001 |
|---|---|---|---|---|---|
| CPPMPAS, GPU | 173.6, 173.6, 173.4 | 173.6 | 4.2 | 1.18 | identical (19 files, 2168 variables) |
| MPAS, 1 rank | 700.8, 704.6, 708.0 | 704.6 | 4.1 | 4.84 | identical |
| MPAS, 24 ranks (1 per core) | 61.1, 61.2, 61.3 | 61.2 | 1.5 | 0.41 | differ from the one-rank run, the same way in all 3 runs |
| MPAS, 48 ranks (2 per core) | 65.1, 65.6, 65.2 | 65.2 | 2.0 | 0.44 | differ from the one-rank run, the same way in all 3 runs |
| CPPMPAS, CPU, 1 thread | 701.1, 700.9, 700.6 | 700.9 | 3.9 | 4.84 | identical |

How the GPU compares, from the medians:
- It is 4.06 times faster than MPAS on one rank.
- MPAS on 24 cores is 2.84 times faster than the GPU.

Caveats:
- This is one VM, and Colab VMs are noisy.
- A 10242-cell mesh is far too small to fill this GPU.
- CPPMPAS synchronises after every kernel launch and has not been tuned for speed.
- MPAS's partition files are NCAR's (x1.10242.tar.gz), and its multi-rank outputs are not bit-identical to the
  one-rank reference.
