// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// What the result of cudaGetDeviceCount means for a GPU executable (src/core/backend_cuda.cu gpu::init). Free of the
// CUDA headers, so that the CPU build's unit test (tests/test_gpu_probe.cpp) checks the same function the CUDA
// build calls; backend_cuda.cu static_asserts the two error codes below against cuda_runtime.h.
//   * cudaSuccess with at least one device: the device is used;
//   * cudaErrorNoDevice, or cudaSuccess with a count of 0: no device on this machine, the only reason for a skip
//     (exit status 77, CTest's SKIP_RETURN_CODE);
//   * any other error (a driver older than the runtime, an initialisation error, no memory, ...): a failure with the
//     runtime's message, never a skip.
namespace cppmpas::gpu {

inline constexpr int kCudaSuccess = 0;          // cudaSuccess (driver_types.h)
inline constexpr int kCudaErrorNoDevice = 100;  // cudaErrorNoDevice (driver_types.h)

enum class DeviceProbe { ready, no_device, error };

// err: cudaGetDeviceCount's return value as an int; count: the count it stored
constexpr DeviceProbe classify_device_count(int err, int count) {
  if (err == kCudaSuccess) return count >= 1 ? DeviceProbe::ready : DeviceProbe::no_device;
  if (err == kCudaErrorNoDevice) return DeviceProbe::no_device;
  return DeviceProbe::error;
}

}  // namespace cppmpas::gpu
