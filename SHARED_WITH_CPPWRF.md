# Code shared with CPPWRF

Ledger required by CONSTITUTION.md §2.6. One row per scheme or routine taken from CPPWRF. A fix in either project
adds the commit in the other project once it is applied and verified there.

| CPPMPAS file | CPPWRF commit | CPPWRF file | MPAS-specific parameters | Fixes (CPPMPAS commit / CPPWRF commit) |
|---|---|---|---|---|
| src/core/backend.hpp (CPPMPAS_HD / CPPMPAS_DEV / CPPMPAS_FORCEINLINE / CPPMPAS_RESTRICT, GPU compiler detection) | 7fac5bf | src/core/gpu_runtime.hpp:20-30, src/core/backend.hpp:21-35 | macro prefix CPPMPAS_ instead of CPPWRF_; `parallel_for` over 1-based element ranges is new (CPPWRF's is over (i,j,k) boxes) | |
| src/core/precision.hpp (`real` typedef, no-excess-precision static_assert) | 7fac5bf | src/core/precision.hpp | FP64 only until the FP32 port exists (CPPMPAS_SINGLE_PRECISION is an error) | |
| CMakeLists.txt strict-IEEE flag policy (-ffp-contract=off, rejected fast-math / x87 flags) | 7fac5bf | CMakeLists.txt:54, 101-122 | -O3 without -march (the MPAS reference is built without -march) | |
| src/core/glibc_math.hpp (glibc 2.41 libm FMA variants and libmvec SSE4 lanes: sin, cos, asin, acos, atan2, pow, mvec::sin/cos/pow and the rest of the file) | 7fac5bf | src/core/glibc_math.hpp | none (copied unchanged; CPPMPAS calls it through src/core/libm.hpp) | (1) asin / acos of \|x\| > 1: the shared asin / acos are the `__ieee754_*` kernels (0/0, -NaN), the public libm asin / acos are the SVID wrappers w_asin_compat.c:30-38 / w_acos_compat.c:30-38 (+NaN). Shared file unchanged; CPPMPAS applies the wrapper in libm.hpp (0f9cc98). Reported to CPPWRF 2026-11-26 (comms/cppmpas_to_cppwrf.md); CPPWRF: 1cd405d on branch fix/ASINNAN (src/core/glibc_math_compat.hpp, the same expression), not yet merged. (2) NaN arguments: once inlined, GCC may fold a caller's negation into the ports' NaN returns (x / x, (x - x) / (x - x), x + y), flipping the NaN's sign against glibc. Shared file unchanged; CPPMPAS's libm.hpp passes every argument through an empty asm (libm::arg) (acc6476). To be reported to CPPWRF, whose wrappers inline the same kernels. |
| src/core/glibc_math_tables.hpp (libm / libmvec data objects) | 7fac5bf | src/core/glibc_math_tables.hpp | none (copied unchanged) | |
| src/core/gpu_runtime.hpp (CPPWRF_GPUCC / CPPWRF_DEVICE_CODE macros glibc_math.hpp needs) | 7fac5bf | src/core/gpu_runtime.hpp | none (copied unchanged) | |
| src/core/libm.hpp `sincos` (__sincos_fma, s_sincos.c) — new, built on glibc_math.hpp's detail::do_sin / do_cos / reduce_sincos / branred | 7fac5bf | (not in CPPWRF; uses src/core/glibc_math.hpp) | to be offered to CPPWRF if WRF ever calls sincos | |
