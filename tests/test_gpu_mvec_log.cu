// SPDX-License-Identifier: GPL-3.0-or-later
// libm::log, log10 and the libmvec lanes mvec_log, mvec_log10 (src/core/libm.hpp) on the GPU against the same
// wrappers on the host, bit for bit, NaN results included (sign and payload): the arguments are passed as they are,
// never negated in the kernel (tests/test_gpu_libm.cu, case A). The inputs: the lanes outside the normal range,
// which libmvec computes with the scalar libm (+-0, every subnormal exponent +-4 ulps, random subnormal bits,
// negatives, +-Inf, quiet and signaling NaNs of either sign with payloads), the GPU audit's inputs (log10 of
// 0x1p-1023 and of the bits 0x0000080078dde6c4: CUDA's log10 gave ...5221 and ...debb there against glibc's
// ...5220 and ...debc), and random normal inputs. tests/test_libm.cpp checks the host wrappers against glibc.
//
// usage: test_gpu_mvec_log [n_random]    exit 0 ok, 1 mismatch, 77 no CUDA device
#include <cuda_runtime.h>

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

#include "core/backend.hpp"
#include "core/libm.hpp"

using namespace cppmpas;

namespace {

constexpr int kFuncs = 4;
const char* kName[kFuncs] = {"log", "log10", "mvec_log", "mvec_log10"};

CPPMPAS_HD void eval(double x, double* out) {
  out[0] = libm::log(x);
  out[1] = libm::log10(x);
  out[2] = libm::mvec_log(x);
  out[3] = libm::mvec_log10(x);
}

__global__ void kern(const double* xs, int n, double* out) {
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < n) eval(xs[i], out + std::size_t(i) * kFuncs);
}

std::uint64_t bits(double v) {
  std::uint64_t b;
  std::memcpy(&b, &v, 8);
  return b;
}
double from_bits(std::uint64_t b) {
  double v;
  std::memcpy(&v, &b, 8);
  return v;
}

}  // namespace

int main(int argc, char** argv) {
  int ndev = 0;
  const cudaError_t dev_err = cudaGetDeviceCount(&ndev);
  if (dev_err == cudaErrorNoDevice || (dev_err == cudaSuccess && ndev < 1)) {
    std::printf("SKIP: no usable CUDA device\n");
    return 77;
  }
  if (dev_err != cudaSuccess) {
    std::printf("FAIL cudaGetDeviceCount: %s\n", cudaGetErrorString(dev_err));
    return 1;
  }
  const int nrand = argc > 1 ? std::atoi(argv[1]) : 200000;
  std::vector<double> xs;
  for (std::uint64_t b : {0x0ull, 0x8000000000000000ull, 0x0000080078dde6c4ull, 0x0008000000000000ull,
                          0x000fffffffffffffull, 0x0000000000000001ull, 0x0010000000000000ull, 0x7fefffffffffffffull,
                          0x7ff0000000000000ull, 0xfff0000000000000ull, 0x7ff8000000000000ull, 0xfff8000000000000ull,
                          0x7ff8000000000123ull, 0xfff8000000000123ull, 0x7ff4000000000001ull, 0xfff4000000000001ull,
                          0x3ff0000000000000ull, 0xbff0000000000000ull, 0x4024000000000000ull})
    xs.push_back(from_bits(b));
  for (int e = 0; e < 52; ++e) {  // the subnormals around each power of two, both signs
    const std::uint64_t b = 1ull << e;
    for (std::uint64_t d = b > 4 ? b - 4 : 1; d <= b + 4; ++d) {
      xs.push_back(from_bits(d));
      xs.push_back(from_bits(d | 0x8000000000000000ull));
    }
  }
  std::uint64_t s = 0x9e3779b97f4a7c15ull;
  auto rnd = [&] {
    s ^= s << 13;
    s ^= s >> 7;
    s ^= s << 17;
    return s;
  };
  for (int i = 0; i < nrand; ++i) {
    const std::uint64_t r = rnd();
    switch (i % 4) {
      case 0: xs.push_back(from_bits(r & 0x000fffffffffffffull)); break;  // a subnormal
      case 1: xs.push_back(from_bits(r | 0x8000000000000000ull)); break;  // a negative (or -0, -Inf, -NaN)
      case 2: xs.push_back(from_bits((r & 0x7fffffffffffffffull) % 0x7ff0000000000000ull)); break;  // any finite >= 0
      default: xs.push_back(double(r >> 11) * 0x1.0p-53 * 2000.0 + 1.0); break;  // pressures, dBZ arguments
    }
  }
  const int n = int(xs.size());
  double *dx, *dout;
  cudaMalloc(&dx, std::size_t(n) * 8);
  cudaMalloc(&dout, std::size_t(n) * kFuncs * 8);
  cudaMemcpy(dx, xs.data(), std::size_t(n) * 8, cudaMemcpyHostToDevice);
  kern<<<(n + 127) / 128, 128>>>(dx, n, dout);
  std::vector<double> gpu(std::size_t(n) * kFuncs);
  if (cudaDeviceSynchronize() != cudaSuccess ||
      cudaMemcpy(gpu.data(), dout, gpu.size() * 8, cudaMemcpyDeviceToHost) != cudaSuccess) {
    std::printf("CUDA error: %s\n", cudaGetErrorString(cudaGetLastError()));
    return 1;
  }
  long bad[kFuncs] = {}, nan_results[kFuncs] = {};
  for (int i = 0; i < n; ++i) {
    double host[kFuncs];
    eval(xs[i], host);
    for (int f = 0; f < kFuncs; ++f) {
      const double h = host[f], g = gpu[std::size_t(i) * kFuncs + f];
      if (h != h) ++nan_results[f];
      if (bits(h) != bits(g) && ++bad[f] <= 5)
        std::printf("DIFF %s(%016llx): host %016llx gpu %016llx\n", kName[f], (unsigned long long)bits(xs[i]),
                    (unsigned long long)bits(h), (unsigned long long)bits(g));
    }
  }
  long total = 0;
  for (int f = 0; f < kFuncs; ++f) {
    std::printf("%-10s %d inputs (%ld NaN results), %ld differ\n", kName[f], n, nan_results[f], bad[f]);
    total += bad[f];
  }
  std::printf(total ? "FAIL\n" : "PASS\n");
  return total ? 1 : 0;
}
