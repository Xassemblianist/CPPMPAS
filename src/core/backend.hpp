// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// Single-source execution backend: kernels are trivially copyable functors whose operator() is CPPMPAS_HD, so the
// same kernel runs on the CPU (OpenMP) and, once the CUDA backend exists, on the GPU.
// Macro scheme taken from CPPWRF 7fac5bf src/core/gpu_runtime.hpp:20-30 and src/core/backend.hpp:21-35
// (SHARED_WITH_CPPWRF.md); renamed CPPWRF_ -> CPPMPAS_, logic unchanged.
#if defined(CPPMPAS_USE_CUDA) || defined(CPPMPAS_USE_HIP)
#define CPPMPAS_USE_GPU 1
#endif
#if defined(__CUDACC__) || defined(__HIP__)
#define CPPMPAS_GPUCC 1
#endif
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
#define CPPMPAS_DEVICE_CODE 1
#endif

#if defined(CPPMPAS_GPUCC)
#define CPPMPAS_HD __host__ __device__
#define CPPMPAS_DEV __device__
#define CPPMPAS_FORCEINLINE __host__ __device__ __forceinline__
#else
#define CPPMPAS_HD
#define CPPMPAS_DEV
#define CPPMPAS_FORCEINLINE inline
#endif

#if defined(__GNUC__) || defined(__clang__) || defined(CPPMPAS_GPUCC)
#define CPPMPAS_RESTRICT __restrict__
#else
#define CPPMPAS_RESTRICT
#endif

namespace cppmpas {

// k(i) for i = first..last (inclusive, MPAS 1-based element indices, e.g. cellStart..cellEnd), elements independent.
// MPAS's own loops with a cross-element accumulation are never written with this: they stay sequential loops.
template <class K> void parallel_for(int first, int last, const K& k) {
#if defined(_OPENMP)
#pragma omp parallel for schedule(static)
#endif
  for (int i = first; i <= last; ++i) k(i);
}

}  // namespace cppmpas
