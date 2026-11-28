// SPDX-License-Identifier: GPL-3.0-or-later
// core/gpu_probe.hpp on the CPU: the classification gpu::init (src/core/backend_cuda.cu) gives the result of
// cudaGetDeviceCount. Only "no device" (cudaErrorNoDevice, or success with a count of 0) is a skip; every other CUDA
// error is a failure. The error codes are the values of CUDA's driver_types.h (backend_cuda.cu static_asserts
// kCudaSuccess and kCudaErrorNoDevice against it).
#include <cstdio>

#include "core/gpu_probe.hpp"

using cppmpas::gpu::classify_device_count;
using cppmpas::gpu::DeviceProbe;

namespace {
int failures = 0;
void expect(int err, int count, DeviceProbe want, const char* what) {
  const DeviceProbe got = classify_device_count(err, count);
  if (got != want) {
    std::printf("FAIL %s (error %d, count %d): got %d, expected %d\n", what, err, count, static_cast<int>(got),
                static_cast<int>(want));
    ++failures;
  }
}
}  // namespace

int main() {
  expect(0, 1, DeviceProbe::ready, "cudaSuccess, one device");
  expect(0, 4, DeviceProbe::ready, "cudaSuccess, four devices");
  expect(0, 0, DeviceProbe::no_device, "cudaSuccess, count 0");
  expect(100, 0, DeviceProbe::no_device, "cudaErrorNoDevice");
  // CUDA errors that are not a missing device (driver_types.h): each is a failure, never a skip
  expect(35, 0, DeviceProbe::error, "cudaErrorInsufficientDriver (driver older than the runtime)");
  expect(2, 0, DeviceProbe::error, "cudaErrorMemoryAllocation");
  expect(3, 0, DeviceProbe::error, "cudaErrorInitializationError");
  expect(719, 0, DeviceProbe::error, "cudaErrorLaunchFailure");
  expect(101, 0, DeviceProbe::error, "cudaErrorInvalidDevice");
  expect(34, 0, DeviceProbe::error, "cudaErrorStubLibrary");
  expect(999, 0, DeviceProbe::error, "cudaErrorUnknown");
  // an error with a count the runtime may have left nonzero is still an error
  expect(35, 1, DeviceProbe::error, "cudaErrorInsufficientDriver, count 1");
  expect(100, 1, DeviceProbe::no_device, "cudaErrorNoDevice, count 1");
  if (failures) return 1;
  std::printf("OK 13 classifications\n");
  return 0;
}
