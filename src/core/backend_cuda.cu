// SPDX-License-Identifier: GPL-3.0-or-later
// The CUDA runtime of core/backend.hpp (CUDA build only): device selection, managed memory, the checks after every
// launch and the kernel error slot. Launch checks and the local-memory reservation follow CPPWRF 7fac5bf
// src/core/backend.cpp (gpu_check_launch, reserve_local_memory; SHARED_WITH_CPPWRF.md).
#include <cuda_runtime.h>

#include <cstdio>
#include <mutex>
#include <string>

#include "core/backend.hpp"

namespace cppmpas {

namespace detail {
__device__ KernelError* kernel_error_slot = nullptr;
}  // namespace detail

namespace gpu {
namespace {
KernelError* g_err = nullptr;  // managed: written by the device, read by the host after the launch
bool g_ready = false;

[[noreturn]] void fail(const char* what, cudaError_t e) {
  throw Error(std::string("CUDA: ") + what + ": " + cudaGetErrorString(e));
}
void ck(cudaError_t e, const char* what) {
  if (e != cudaSuccess) fail(what, e);
}
}  // namespace

bool init(std::string& why_not) {
  if (g_ready) return true;
  int n = 0;
  const cudaError_t e = cudaGetDeviceCount(&n);
  if (e != cudaSuccess || n < 1) {
    why_not = e != cudaSuccess ? std::string(cudaGetErrorString(e)) : std::string("no CUDA device");
    (void)cudaGetLastError();
    return false;
  }
  ck(cudaSetDevice(0), "cudaSetDevice");
  void* p = nullptr;
  ck(cudaMallocManaged(&p, sizeof(KernelError)), "cudaMallocManaged (kernel error slot)");
  g_err = static_cast<KernelError*>(p);
  *g_err = KernelError{};
  ck(cudaMemcpyToSymbol(detail::kernel_error_slot, &g_err, sizeof g_err), "cudaMemcpyToSymbol (kernel error slot)");
  ck(cudaDeviceSynchronize(), "cudaDeviceSynchronize");
  g_ready = true;
  return true;
}

std::string describe() {
  int dev = 0, drv = 0, rt = 0;
  cudaDeviceProp p{};
  ck(cudaGetDevice(&dev), "cudaGetDevice");
  ck(cudaGetDeviceProperties(&p, dev), "cudaGetDeviceProperties");
  ck(cudaDriverGetVersion(&drv), "cudaDriverGetVersion");
  ck(cudaRuntimeGetVersion(&rt), "cudaRuntimeGetVersion");
  char buf[512];
  std::snprintf(buf, sizeof buf,
                "device %d: %s, compute capability %d.%d, %d SMs, %.0f MiB; CUDA driver API %d.%d, runtime %d.%d",
                dev, p.name, p.major, p.minor, p.multiProcessorCount, double(p.totalGlobalMem) / (1024.0 * 1024.0),
                drv / 1000, (drv % 1000) / 10, rt / 1000, (rt % 1000) / 10);
  return buf;
}

void* managed_alloc(std::size_t bytes) {
  if (!g_ready) {
    std::string why;
    if (!init(why)) throw Error("CUDA build: no usable CUDA device (" + why + ")");
  }
  if (bytes == 0) bytes = 1;
  void* p = nullptr;
  const cudaError_t e = cudaMallocManaged(&p, bytes);
  if (e != cudaSuccess) {
    (void)cudaGetLastError();
    std::size_t free_b = 0, total_b = 0;
    cudaMemGetInfo(&free_b, &total_b);
    char msg[320];
    std::snprintf(msg, sizeof msg,
                  "device memory too small: a managed allocation of %.0f MiB failed (%s); %.0f MiB of %.0f MiB free",
                  double(bytes) / 1048576.0, cudaGetErrorString(e), double(free_b) / 1048576.0,
                  double(total_b) / 1048576.0);
    throw Error(msg);
  }
  return p;
}
void managed_free(void* p) noexcept {
  if (p) cudaFree(p);
}

void after_launch(const char* what) {
  cudaError_t e = cudaGetLastError();
  if (e != cudaSuccess) fail(what, e);
  e = cudaDeviceSynchronize();
  if (e != cudaSuccess) fail(what, e);
  if (g_err && g_err->set) {
    const KernelError k = *g_err;
    *g_err = KernelError{};
    detail::throw_kernel_error(k.head, k.has_value != 0, k.value, k.tail);
  }
}

void reserve_local_memory(const void* /*kernel*/, std::size_t local_bytes) {
  static std::mutex mu;
  static std::size_t reserved = 0;
  std::lock_guard<std::mutex> lock(mu);
  if (local_bytes <= reserved) return;
  ck(cudaDeviceSynchronize(), "cudaDeviceSynchronize");
  ck(cudaDeviceSetLimit(cudaLimitStackSize, local_bytes), "cudaDeviceSetLimit (stack for a kernel's local memory)");
  reserved = local_bytes;
}

}  // namespace gpu
}  // namespace cppmpas
