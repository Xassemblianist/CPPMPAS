# Code shared with CPPWRF

Ledger required by CONSTITUTION.md §2.6. One row per scheme or routine taken from CPPWRF. A fix in either project
adds the commit in the other project once it is applied and verified there.

| CPPMPAS file | CPPWRF commit | CPPWRF file | MPAS-specific parameters | Fixes (CPPMPAS commit / CPPWRF commit) |
|---|---|---|---|---|
| src/core/backend.hpp (CPPMPAS_HD / CPPMPAS_DEV / CPPMPAS_FORCEINLINE / CPPMPAS_RESTRICT, GPU compiler detection) | 7fac5bf | src/core/gpu_runtime.hpp:20-30, src/core/backend.hpp:21-35 | macro prefix CPPMPAS_ instead of CPPWRF_; `parallel_for` over 1-based element ranges is new (CPPWRF's is over (i,j,k) boxes) | |
| src/core/precision.hpp (`real` typedef, no-excess-precision static_assert) | 7fac5bf | src/core/precision.hpp | FP64 only until the FP32 port exists (CPPMPAS_SINGLE_PRECISION is an error) | |
| CMakeLists.txt strict-IEEE flag policy (-ffp-contract=off, rejected fast-math / x87 flags) | 7fac5bf | CMakeLists.txt:54, 101-122 | -O3 without -march (the MPAS reference is built without -march) | |
