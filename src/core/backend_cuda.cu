// SPDX-License-Identifier: GPL-3.0-or-later
// The CUDA runtime of core/backend.hpp (CUDA build only): device selection, managed memory, the checks after every
// launch and the kernel error slot. Launch checks and the local-memory reservation follow CPPWRF 7fac5bf
// src/core/backend.cpp (gpu_check_launch, reserve_local_memory; SHARED_WITH_CPPWRF.md).
#include <cuda_runtime.h>
#include <cxxabi.h>

#include <cstdio>
#include <cstdlib>
#include <map>
#include <mutex>
#include <string>

#include "core/backend.hpp"
#include "core/gpu_probe.hpp"

namespace cppmpas {

namespace detail {
__device__ KernelError* kernel_error_slot = nullptr;
}  // namespace detail

namespace gpu {
static_assert(static_cast<int>(cudaSuccess) == kCudaSuccess && static_cast<int>(cudaErrorNoDevice) == kCudaErrorNoDevice,
              "core/gpu_probe.hpp: the CUDA error codes differ from cuda_runtime.h");
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
  // core/gpu_probe.hpp: false (a skip) only when the machine has no device; any other CUDA error is an Error
  switch (classify_device_count(static_cast<int>(e), n)) {
    case DeviceProbe::ready:
      break;
    case DeviceProbe::no_device:
      why_not = e != cudaSuccess ? std::string(cudaGetErrorString(e)) : std::string("no CUDA device");
      (void)cudaGetLastError();
      return false;
    case DeviceProbe::error:
      (void)cudaGetLastError();
      fail("cudaGetDeviceCount", e);
  }
  ck(cudaSetDevice(0), "cudaSetDevice");
  void* p = nullptr;
  ck(cudaMallocManaged(&p, sizeof(KernelError)), "cudaMallocManaged (kernel error slot)");
  KernelError* slot = static_cast<KernelError*>(p);
  *slot = KernelError{};
  cudaError_t e2 = cudaMemcpyToSymbol(detail::kernel_error_slot, &slot, sizeof slot);
  if (e2 == cudaSuccess) e2 = cudaDeviceSynchronize();
  if (e2 != cudaSuccess) {  // the slot is not kept: a later init starts again
    (void)cudaFree(p);
    (void)cudaGetLastError();
    fail("setting the kernel error slot", e2);
  }
  g_err = slot;
  g_ready = true;
  return true;
}

void ensure_ready() {
  if (g_ready) return;
  std::string why;
  if (!init(why)) throw Error("CUDA build: no usable CUDA device (" + why + ")");
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
    if (cudaMemGetInfo(&free_b, &total_b) != cudaSuccess) {
      (void)cudaGetLastError();
      free_b = total_b = 0;
    }
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
  if (p && cudaFree(p) != cudaSuccess) (void)cudaGetLastError();
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

// CPPWRF 7fac5bf src/core/backend.cpp:384-428, without its stream capture (every launch here synchronises)
void reserve_local_memory(const void* kernel, const char* kernel_type) {
  cudaFuncAttributes a{};
  ck(cudaFuncGetAttributes(&a, kernel), "cudaFuncGetAttributes");
  static std::mutex m;
  static std::map<int, std::size_t> reserved_on;  // per device: the largest frame of the process's launches
  int cur = 0;
  ck(cudaGetDevice(&cur), "cudaGetDevice");
  std::lock_guard<std::mutex> lk(m);
  std::size_t& reserved = reserved_on[cur];
  if (a.localSizeBytes <= reserved) return;
  std::size_t limit = 0;
  ck(cudaDeviceGetLimit(&limit, cudaLimitStackSize), "cudaDeviceGetLimit");
  if (a.localSizeBytes <= limit) {  // the limit is never lowered
    reserved = a.localSizeBytes;
    return;
  }
  int sms = 0, tpsm = 0;
  ck(cudaDeviceGetAttribute(&sms, cudaDevAttrMultiProcessorCount, cur), "cudaDeviceGetAttribute");
  ck(cudaDeviceGetAttribute(&tpsm, cudaDevAttrMaxThreadsPerMultiProcessor, cur), "cudaDeviceGetAttribute");
  ck(cudaDeviceSynchronize(), "cudaDeviceSynchronize");
  std::size_t free_b = 0, total_b = 0;
  const cudaError_t mi = cudaMemGetInfo(&free_b, &total_b);
  if (mi != cudaSuccess) (void)cudaGetLastError();
  const cudaError_t e = cudaDeviceSetLimit(cudaLimitStackSize, a.localSizeBytes);
  if (e != cudaSuccess) {
    (void)cudaGetLastError();
    int st = 0;
    char* dm = abi::__cxa_demangle(kernel_type, nullptr, nullptr, &st);
    const std::string name = (st == 0 && dm) ? dm : kernel_type;
    std::free(dm);
    const double mib = 1024. * 1024.;
    char msg[640];
    std::snprintf(msg, sizeof msg,
                  "device memory too small: the kernel %.200s keeps %zu bytes per thread in local memory, which the GPU "
                  "reserves for each of its %d x %d resident threads (%.0f MiB in all, %.0f MiB of it reserved so far); "
                  "%.0f MiB of %.0f MiB are free (%s)",
                  name.c_str(), std::size_t(a.localSizeBytes), sms, tpsm, double(a.localSizeBytes) * sms * tpsm / mib,
                  double(limit) * sms * tpsm / mib, mi == cudaSuccess ? double(free_b) / mib : 0.,
                  mi == cudaSuccess ? double(total_b) / mib : 0., cudaGetErrorString(e));
    throw Error(msg);
  }
  reserved = a.localSizeBytes;
}

}  // namespace gpu
}  // namespace cppmpas
