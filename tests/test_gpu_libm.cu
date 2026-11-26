// SPDX-License-Identifier: GPL-3.0-or-later
// The libm / libmvec wrappers of src/core/libm.hpp on the GPU (CUDA build) against the same wrappers on the host:
// sin, cos, asin, acos, atan2, pow, sincos and the libmvec lanes mvec_sin, mvec_cos, mvec_pow, on special values
// (+-0, +-Inf, +-qNaN, +-sNaN, subnormals, |x| > 1 for asin / acos, ...) and on pseudo-random arguments, each also
// passed negated (-x computed in the kernel from a value in memory, as a caller's -x reaches the port).
// Rule: every non-NaN result is bit-identical (signed zero included); a NaN result is a NaN on both; a NaN result is
// bit-identical (sign and payload) too, except in the two cases measured on the sm_120 card (RTX 5070 Ti, driver
// 621.24, CUDA 13.1) and recorded here, which must occur exactly as often as recorded (a change fails the test):
//  (A) the argument is a NaN that the caller negated: on the device -x of a NaN is that NaN with its sign unchanged
//      (and an sNaN quieted), for -x and -1 * x alike, while the host flips the sign bit; the function then sees
//      the other sign and returns it. This is the hardware's / nvcc's negation, not the libm port: the result is the
//      same with and without libm::arg's device barrier (CPPMPAS_NO_DEVICE_ARG_BARRIER, measured);
//  (B) both arguments of atan2, pow or mvec_pow are NaN: glibc returns x + y, which is either argument's NaN
//      depending on which operand the compiler puts first (SSE addsd keeps its destination's); on the host that is a
//      register choice of the host compiler (measured: pow differs with g++ 14 on the sm_120 build, atan2 and
//      mvec_pow with g++ 13 on the L4 build), so at most one payload difference per function is allowed there.
// The sign differences (case A) are the device's and are counted exactly (the same on sm_120 and sm_89).
// No MPAS field compared by the GPU tests holds a NaN (they compare raw bytes at tolerance 0, unweakened).
//
// usage: test_gpu_libm [n_random]    exit 0 ok, 1 mismatch, 77 no CUDA device
#include <cuda_runtime.h>

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include "core/backend.hpp"
#include "core/libm.hpp"

using namespace cppmpas;

namespace {

constexpr int kFuncs = 12;
const char* kName[kFuncs] = {"sin",  "cos",      "asin",     "acos",     "atan2(y,x)", "pow(x,y)",
                             "sincos.sin", "sincos.cos", "mvec_sin", "mvec_cos", "mvec_pow(x,y)", "-x"};

// one evaluation of every function; neg: the arguments are negated first (in the caller, as -x)
CPPMPAS_HD void eval(double x, double y, bool neg, double* out) {
  if (neg) {
    x = -x;
    y = -y;
  }
  out[0] = libm::sin(x);
  out[1] = libm::cos(x);
  out[2] = libm::asin(x);
  out[3] = libm::acos(x);
  out[4] = libm::atan2(y, x);
  out[5] = libm::pow(x, y);
  libm::sincos(x, &out[6], &out[7]);
  out[8] = libm::mvec_sin(x);
  out[9] = libm::mvec_cos(x);
  out[10] = libm::mvec_pow(x, y);
  out[11] = x;  // the negated argument itself: what -x of a NaN gives
}

__global__ void kern(const double* xs, const double* ys, int n, double* out) {
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= 2 * n) return;
  const int j = i % n;
  eval(xs[j], ys[j], i >= n, out + std::size_t(i) * kFuncs);
}

// the three ways the device can negate a NaN read from memory
__global__ void negs(const double* x, double* out) {
  out[0] = -x[0];
  out[1] = 0.0 - x[0];
  out[2] = -1.0 * x[0];
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
  if (cudaGetDeviceCount(&ndev) != cudaSuccess || ndev < 1) {
    std::printf("SKIP: no usable CUDA device\n");
    return 77;
  }
  const int nrand = argc > 1 ? std::atoi(argv[1]) : 200000;
  std::vector<double> xs = {0.0, -0.0, 1.0, -1.0, 0.5, 2.0, 1.0000000000000002, 0.9999999999999999, 3.0e-310,
                            -3.0e-310, 1e300, -1e300, 0x1.921fb54442d18p0, 0x1.921fb54442d18p1, 1e22, 710.0,
                            from_bits(0x7ff0000000000000ull), from_bits(0xfff0000000000000ull),
                            from_bits(0x7ff8000000000000ull), from_bits(0xfff8000000000000ull),
                            from_bits(0x7ff4000000000001ull), from_bits(0x7ff8000000000123ull)};
  std::vector<double> ys = {2.0, 0.5, -1.0, 0.25, 3.0, -0.0, 0.0, 1e-300, 0.75, 1.5, 2.0, 0.25, 1.0, 0.0, 0.5, 0.1,
                            1.0, 2.0, 1.0, 0.5, 2.0, from_bits(0x7ff8000000000000ull)};
  // NaN y with ordinary x (atan2 / pow of a NaN second argument)
  for (double x : {0.5, -2.0, 1.0}) {
    xs.push_back(x);
    ys.push_back(from_bits(0xfff8000000000000ull));
  }
  std::uint64_t s = 0x9e3779b97f4a7c15ull;
  auto rnd = [&] {
    s ^= s << 13;
    s ^= s >> 7;
    s ^= s << 17;
    return s;
  };
  for (int i = 0; i < nrand; ++i) {
    const double u = double(rnd() >> 11) * 0x1.0p-53;
    const int range = int(rnd() % 4);
    const double x = range == 0 ? u * 2.0 - 1.0 : range == 1 ? (u - 0.5) * 20.0 : range == 2 ? (u - 0.5) * 2e6 : u * 3.0;
    const double y = double(rnd() >> 11) * 0x1.0p-53 * 4.0 - 1.0;
    xs.push_back(x);
    ys.push_back(range == 3 ? y * 0.5 + 0.2857142857142857 : y);  // pow near rcv = R/cv as the dycore calls it
  }
  const int n = int(xs.size());
  double *dx, *dy, *dout;
  cudaMalloc(&dx, n * 8);
  cudaMalloc(&dy, n * 8);
  cudaMalloc(&dout, std::size_t(2) * n * kFuncs * 8);
  cudaMemcpy(dx, xs.data(), n * 8, cudaMemcpyHostToDevice);
  cudaMemcpy(dy, ys.data(), n * 8, cudaMemcpyHostToDevice);
  kern<<<(2 * n + 127) / 128, 128>>>(dx, dy, n, dout);
  std::vector<double> gpu(std::size_t(2) * n * kFuncs);
  if (cudaDeviceSynchronize() != cudaSuccess ||
      cudaMemcpy(gpu.data(), dout, gpu.size() * 8, cudaMemcpyDeviceToHost) != cudaSuccess) {
    std::printf("CUDA error: %s\n", cudaGetErrorString(cudaGetLastError()));
    return 1;
  }
  long finite_bad = 0, nanness_bad = 0, compared = 0;
  long nan_total[kFuncs] = {}, nan_sign_diff[kFuncs] = {}, nan_payload_diff[kFuncs] = {}, unexplained = 0;
  for (int i = 0; i < 2 * n; ++i) {
    double host[kFuncs];
    eval(xs[i % n], ys[i % n], i >= n, host);
    for (int f = 0; f < kFuncs; ++f) {
      const double h = host[f], g = gpu[std::size_t(i) * kFuncs + f];
      ++compared;
      if (std::isnan(h) != std::isnan(g)) {
        if (++nanness_bad <= 10)
          std::printf("NaN-ness differs: %s x=%a y=%a neg=%d host %016llx gpu %016llx\n", kName[f], xs[i % n],
                      ys[i % n], i >= n, (unsigned long long)bits(h), (unsigned long long)bits(g));
      } else if (std::isnan(h)) {
        ++nan_total[f];
        const bool neg_nan = i >= n && (std::isnan(xs[i % n]) || std::isnan(ys[i % n]));          // case (A)
        const bool two_nan =
            std::isnan(xs[i % n]) && std::isnan(ys[i % n]) && (f == 4 || f == 5 || f == 10);  // case (B)
        if (bits(h) != bits(g) && !neg_nan && !(two_nan && (bits(h) >> 63) == (bits(g) >> 63))) {
          ++unexplained;
          std::printf("UNEXPLAINED NaN difference: %s x=%a y=%a neg=%d host %016llx gpu %016llx\n", kName[f],
                      xs[i % n], ys[i % n], i >= n, (unsigned long long)bits(h), (unsigned long long)bits(g));
        }
        if ((bits(h) >> 63) != (bits(g) >> 63)) {
          if (++nan_sign_diff[f] <= 3)
            std::printf("NaN sign differs: %s x=%a y=%a neg=%d host %016llx gpu %016llx\n", kName[f], xs[i % n],
                        ys[i % n], i >= n, (unsigned long long)bits(h), (unsigned long long)bits(g));
        } else if (bits(h) != bits(g)) {
          if (++nan_payload_diff[f] <= 3)
            std::printf("NaN payload differs: %s x=%a y=%a neg=%d host %016llx gpu %016llx\n", kName[f], xs[i % n],
                        ys[i % n], i >= n, (unsigned long long)bits(h), (unsigned long long)bits(g));
        }
      } else if (bits(h) != bits(g)) {
        if (++finite_bad <= 10)
          std::printf("DIFF %s x=%a y=%a neg=%d host %a gpu %a\n", kName[f], xs[i % n], ys[i % n], i >= n, h, g);
      }
    }
  }
  std::printf("%ld results compared (%d arguments, plain and negated, %d functions)\n", compared, n, kFuncs);
  std::printf("non-NaN results that differ (signed zero included): %ld\n", finite_bad);
  std::printf("results NaN on one side only: %ld\n", nanness_bad);
  long sign_total = 0;
  for (int f = 0; f < kFuncs; ++f) {
    std::printf("  %-14s NaN results %4ld, sign differs %4ld, payload differs %4ld\n", kName[f], nan_total[f],
                nan_sign_diff[f], nan_payload_diff[f]);
    sign_total += nan_sign_diff[f] + nan_payload_diff[f];
  }
  // negating a NaN on the device
  const double qn = from_bits(0x7ff8000000000000ull);
  double* dq;
  double* dn;
  cudaMalloc(&dq, 8);
  cudaMalloc(&dn, 3 * 8);
  cudaMemcpy(dq, &qn, 8, cudaMemcpyHostToDevice);
  negs<<<1, 1>>>(dq, dn);
  double ng[3];
  cudaMemcpy(ng, dn, 3 * 8, cudaMemcpyDeviceToHost);
  volatile double hq = qn;
  const double nh[3] = {-hq, 0.0 - hq, -1.0 * hq};
  const char* how[3] = {"-x", "0 - x", "-1 * x"};
  for (int k = 0; k < 3; ++k)
    std::printf("device %-7s of +qNaN: %016llx (host %016llx)\n", how[k], (unsigned long long)bits(ng[k]),
                (unsigned long long)bits(nh[k]));
  // the sign counts measured on the sm_120 and sm_89 cards (case A), exactly; at most one payload difference per
  // function of case B (the host compiler's operand order), none elsewhere
  const long want_sign[kFuncs] = {4, 4, 4, 4, 7, 6, 4, 4, 4, 4, 6, 4};
  const long max_payload[kFuncs] = {0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0};
  long recorded_bad = 0;
  for (int f = 0; f < kFuncs; ++f)
    if (nan_sign_diff[f] != want_sign[f] || nan_payload_diff[f] > max_payload[f]) {
      ++recorded_bad;
      std::printf("%s: NaN sign / payload differences %ld / %ld, recorded %ld / at most %ld\n", kName[f],
                  nan_sign_diff[f], nan_payload_diff[f], want_sign[f], max_payload[f]);
    }
  const bool dev_neg_keeps_sign = (bits(ng[0]) >> 63) == 0 && (bits(ng[2]) >> 63) == 0;
  std::printf("device -x of a NaN keeps its sign: %s (recorded: yes)\n", dev_neg_keeps_sign ? "yes" : "no");
  std::printf("NaN results with another sign or payload than the host's: %ld, all of cases A / B: %s\n", sign_total,
              unexplained == 0 ? "yes" : "no");
  if (unexplained || recorded_bad || !dev_neg_keeps_sign) return 1;
  if (finite_bad || nanness_bad) return 1;
  return 0;
}
