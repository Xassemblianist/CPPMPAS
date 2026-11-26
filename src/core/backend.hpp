// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// Single-source execution backend: a kernel is a lambda over one MPAS element index (a cell, edge or vertex, MPAS's
// 1-based numbering) that captures only trivially copyable views (core/views.hpp), so the same source runs on the CPU
// (a loop, OpenMP when CPPMPAS_OPENMP) and on the GPU (CUDA, one thread per element).
// Macro scheme and launch design taken from CPPWRF 7fac5bf src/core/gpu_runtime.hpp:20-30 and src/core/backend.hpp:
// 21-35, 647-722 (k3d / parallel_for / kred_seq; SHARED_WITH_CPPWRF.md); renamed CPPWRF_ -> CPPMPAS_. MPAS
// parameters: 1-based element ranges [first, last] instead of CPPWRF's half-open (i, j, k) boxes, one thread per
// element (the k loop and every neighbour loop stay inside the thread, in MPAS's order), and a one-thread launch
// (serial) for the loops MPAS runs as one sequence.
//
// Writing a kernel (the statements inside are MPAS's, unchanged):
//   parallel_for(1, m.nEdges, [=, m = view(m), f = view(f)] CPPMPAS_LAMBDA(int iEdge) { ... });
// The init-captures make the names the body uses (m.cellsOnEdge, f.u2, ...) views of the same arrays; CPPMPAS_LAMBDA
// is __device__ in the CUDA build (CUDA allows init-captures only in device lambdas) and empty in the CPU build.
// A kernel never throws: kernel_fail(msg) throws Error(msg) on the CPU and, on the GPU, records msg for the launch,
// which throws the same Error when it returns.
#include <cstddef>
#include <exception>
#include <string>

#include "core/error.hpp"

#if defined(CPPMPAS_USE_CUDA) || defined(CPPMPAS_USE_HIP)
#define CPPMPAS_USE_GPU 1
#endif
#if defined(__CUDACC__) || defined(__HIP__)
#define CPPMPAS_GPUCC 1
#endif
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
#define CPPMPAS_DEVICE_CODE 1
#endif
#if defined(CPPMPAS_USE_HIP)
#error "CPPMPAS_USE_HIP: there is no HIP port"
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

// the annotation of a kernel lambda (see above)
#if defined(CPPMPAS_USE_CUDA) && defined(CPPMPAS_GPUCC)
#define CPPMPAS_LAMBDA __device__
#else
#define CPPMPAS_LAMBDA
#endif

// every source of a CUDA build is compiled by nvcc (CMakeLists.txt), as in CPPWRF
#if defined(CPPMPAS_USE_CUDA) && !defined(CPPMPAS_GPUCC)
#error "CPPMPAS_USE_CUDA: this source must be compiled by nvcc"
#endif

namespace cppmpas {

// ---- kernel errors ----
// The first failure of a launch: its message, with an optional integer between two parts (e.g. a column number).
struct KernelError {
  int set;
  int has_value;
  long value;
  char head[192];
  char tail[96];
};
namespace detail {
[[noreturn]] void throw_kernel_error(const char* head, bool has_value, long value, const char* tail);
#if defined(CPPMPAS_USE_CUDA)
extern __device__ KernelError* kernel_error_slot;  // managed memory, set once by the runtime (backend_cuda.cu)
#if defined(CPPMPAS_DEVICE_CODE)
__device__ inline void copy_message(char* dst, int cap, const char* src) {
  int i = 0;
  if (src)
    for (; i < cap - 1 && src[i] != '\0'; ++i) dst[i] = src[i];
  dst[i] = '\0';
}
__device__ inline void record_kernel_error(const char* head, bool has_value, long value, const char* tail) {
  KernelError* e = kernel_error_slot;
  if (e && atomicCAS(&e->set, 0, 1) == 0) {
    copy_message(e->head, int(sizeof e->head), head);
    copy_message(e->tail, int(sizeof e->tail), tail);
    e->has_value = has_value ? 1 : 0;
    e->value = value;
    __threadfence_system();
  }
}
#endif
#endif
}  // namespace detail

// Error(msg) on the CPU; on the GPU the launch throws it (the kernel returns right after the call).
CPPMPAS_HD inline void kernel_fail(const char* msg) {
#if defined(CPPMPAS_DEVICE_CODE)
  detail::record_kernel_error(msg, false, 0, "");
#else
  detail::throw_kernel_error(msg, false, 0, "");
#endif
}
CPPMPAS_HD inline void kernel_fail(const char* head, long value, const char* tail) {
#if defined(CPPMPAS_DEVICE_CODE)
  detail::record_kernel_error(head, true, value, tail);
#else
  detail::throw_kernel_error(head, true, value, tail);
#endif
}

// ---- device runtime (CUDA build; the CPU build has none) ----
#if defined(CPPMPAS_USE_CUDA)
namespace gpu {
// Exit status of a GPU executable when there is no usable CUDA device: CTest's SKIP_RETURN_CODE.
inline constexpr int kNoDeviceExit = 77;
// Selects device 0 and sets the kernel error slot. false, with the reason, when there is no usable device.
bool init(std::string& why_not);
// "NVIDIA GeForce ..., compute capability 12.0, driver 13.x (CUDA driver API), runtime 13.1, ..."
std::string describe();
void* managed_alloc(std::size_t bytes);  // cudaMallocManaged; Error when the device memory is short
void managed_free(void* p) noexcept;
// After every launch: wait for the device (managed memory is then the host's again: the GPUs this build targets
// may have no concurrent managed access), and throw Error for a launch failure or a kernel_fail.
void after_launch(const char* what);
// Before the first launch of a kernel: raise the stack limit to its local memory (CPPWRF backend.cpp
// reserve_local_memory), so a kernel with column work arrays does not fail at launch.
void reserve_local_memory(const void* kernel, std::size_t local_bytes);
}  // namespace gpu
#endif

#if defined(CPPMPAS_USE_CUDA) && defined(CPPMPAS_GPUCC)
namespace detail {
template <class K> __global__ void kfor(K k, int first, int last) {
  const long i = long(first) + long(blockIdx.x) * long(blockDim.x) + long(threadIdx.x);
  if (i <= long(last)) k(int(i));
}
template <class K> __global__ void kfor_n(K k, std::size_t n) {
  const std::size_t i = std::size_t(blockIdx.x) * blockDim.x + threadIdx.x;
  if (i < n) k(i);
}
template <class K> __global__ void kserial(K k) { k(); }
template <class Kern> void reserve_once(Kern* kernel) {
  static const bool done = [&] {
    cudaFuncAttributes a{};
    if (cudaFuncGetAttributes(&a, reinterpret_cast<const void*>(kernel)) == cudaSuccess)
      gpu::reserve_local_memory(reinterpret_cast<const void*>(kernel), a.localSizeBytes);
    return true;
  }();
  (void)done;
}
inline constexpr int kBlock = 128;
}  // namespace detail

// k(i) for i = first..last (inclusive, MPAS 1-based element indices), elements independent: one GPU thread each.
template <class K> void parallel_for(int first, int last, const K& k) {
  if (last < first) return;
  detail::reserve_once(&detail::kfor<K>);
  const long n = long(last) - long(first) + 1;
  detail::kfor<K><<<unsigned((n + detail::kBlock - 1) / detail::kBlock), detail::kBlock>>>(k, first, last);
  gpu::after_launch("parallel_for");
}
// k(i) for i = 0..n-1 over the elements of a whole array (garbage element included), elements independent.
template <class K> void parallel_for_n(std::size_t n, const K& k) {
  if (n == 0) return;
  detail::reserve_once(&detail::kfor_n<K>);
  detail::kfor_n<K><<<unsigned((n + detail::kBlock - 1) / detail::kBlock), detail::kBlock>>>(k, n);
  gpu::after_launch("parallel_for_n");
}
// k() once, on one GPU thread: a loop MPAS runs as one sequence (a carried dependence or a global order).
template <class K> void serial(const K& k) {
  detail::reserve_once(&detail::kserial<K>);
  detail::kserial<K><<<1, 1>>>(k);
  gpu::after_launch("serial");
}
#else
// k(i) for i = first..last (inclusive, MPAS 1-based element indices), elements independent.
// With OpenMP a kernel_fail (an exception) must not leave the parallel region: the first one is kept and rethrown.
template <class I, class K> void cpu_for(I first, I end, const K& k) {
#if defined(_OPENMP)
  std::exception_ptr err;
#pragma omp parallel for schedule(static)
  for (I i = first; i < end; ++i) {
    try {
      k(i);
    } catch (...) {
#pragma omp critical(cppmpas_kernel_fail)
      if (!err) err = std::current_exception();
    }
  }
  if (err) std::rethrow_exception(err);
#else
  for (I i = first; i < end; ++i) k(i);
#endif
}
template <class K> void parallel_for(int first, int last, const K& k) {
  if (last >= first) cpu_for<long>(long(first), long(last) + 1, [&](long i) { k(int(i)); });
}
template <class K> void parallel_for_n(std::size_t n, const K& k) { cpu_for<std::size_t>(0, n, k); }
template <class K> void serial(const K& k) { k(); }
#endif

// Fixed-capacity work array of a kernel, in place of MPAS's automatic arrays (std::vector on the CPU before):
// a[0..n-1] zero-filled as std::vector<T>(n) is, capacity N. The caller checks n <= N on the host first
// (fits_column); a larger n is a kernel_fail (the kernel's results are then not used: the launch throws).
template <class T, int N> struct LocalArray {
  T v[N];
  int n_;
  CPPMPAS_HD explicit LocalArray(int n) : n_(n < 0 ? 0 : n) {
    if (n_ > N) {
      kernel_fail("a kernel work array is larger than its capacity");
      n_ = N;
    }
    for (int i = 0; i < n_; ++i) v[i] = T{};
  }
  CPPMPAS_HD T& operator[](int i) { return v[i]; }
  CPPMPAS_HD const T& operator[](int i) const { return v[i]; }
  CPPMPAS_HD T* data() { return v; }
  CPPMPAS_HD T* begin() { return v; }
  CPPMPAS_HD T* end() { return v + n_; }
  CPPMPAS_HD int size() const { return n_; }
};
// the capacity of the vertical work arrays: nVertLevels + 2 elements (index 0 unused, 1..nVertLevels+1)
inline constexpr int kMaxColumn = 258;
template <class T> using Column = LocalArray<T, kMaxColumn>;
// Error naming the routine unless a column of nVertLevels levels (nVertLevels + 2 elements) fits a Column
void fits_column(int nVertLevels, const char* routine);

}  // namespace cppmpas
