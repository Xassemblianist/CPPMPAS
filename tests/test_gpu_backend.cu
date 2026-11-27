// SPDX-License-Identifier: GPL-3.0-or-later
// The launch wrappers of src/core/backend.hpp on the GPU (CUDA build), each on a small array whose every element,
// the ones outside the range included, is checked afterwards:
//  1. a kernel_fail in the first launch of the process, before anything initialised the runtime, reaches the host
//     as its Error (the launch initialises the kernel error slot itself);
//  2. parallel_for(first, last) runs k(i) exactly once for each i in first..last, and for no other i;
//  3. parallel_for_n(n) runs k(i) exactly once for each i in 0..n-1 (the garbage element of an MPAS array included);
//  4. serial runs k once;
//  5. kernel_fail with a value: the launch throws "head<value>tail", and no statement after the call runs (the thread
//     ends there; a marker written after it stays untouched), also from a function the kernel calls;
//  6. a work array larger than its capacity is a kernel_fail before the array is used (no write past it);
//  7. after a failed launch the next launch runs normally;
//  8. the device's stack limit is never lowered by a launch (CPPWRF reserve_local_memory).
//
// usage: test_gpu_backend    exit 0 ok, 1 failure, 77 no CUDA device
#include <cuda_runtime.h>

#include <cstdio>
#include <string>

#include "core/backend.hpp"
#include "core/field.hpp"

using namespace cppmpas;

namespace {

int failures = 0;
void check(bool ok, const char* what) {
  std::printf("%s %s\n", ok ? "ok  " : "FAIL", what);
  if (!ok) ++failures;
}

// the Error message of f(), or "" when it does not throw
template <class F> std::string error_of(const F& f) {
  try {
    f();
  } catch (const Error& e) {
    return e.what();
  }
  return "";
}

// a function a kernel calls, which fails: the caller's statements after the call must not run either
CPPMPAS_HD int fails_deep(int i) {
  if (i == 3) kernel_fail("test_gpu_backend: element ", i, " failed");
  return i;
}

}  // namespace

int main() {
  int ndev = 0;
  if (cudaGetDeviceCount(&ndev) != cudaSuccess || ndev < 1) {
    std::printf("SKIP: no usable CUDA device\n");
    return 77;
  }
  std::size_t stack0 = 0;
  if (cudaDeviceGetLimit(&stack0, cudaLimitStackSize) != cudaSuccess) {
    std::printf("FAIL cudaDeviceGetLimit\n");
    return 1;
  }

  // 1. the first launch of the process fails: its Error reaches the host
  check(error_of([] { serial([=] CPPMPAS_LAMBDA() { kernel_fail("test_gpu_backend: first launch"); }); }) ==
            "test_gpu_backend: first launch",
        "a kernel_fail in the first launch of the process is the launch's Error");

  constexpr int n = 300;  // more than two blocks of threads
  FArray<int> a({n}, 0);
  const FView<int> av = a.view();

  // 2. parallel_for over 1-based 2..n-1: each element of the range once, elements 1 and n untouched
  parallel_for(2, n - 1, [=] CPPMPAS_LAMBDA(int i) { av(i) = av(i) + i; });
  bool ok = a(1) == 0 && a(n) == 0;
  for (int i = 2; i <= n - 1; ++i) ok = ok && a(i) == i;
  check(ok, "parallel_for(first, last): every element of first..last once, none outside");

  // 3. parallel_for_n over the whole array, its last element included
  parallel_for_n(a.size(), [=] CPPMPAS_LAMBDA(std::size_t i) { av[i] = av[i] + 1000; });
  ok = a(1) == 1000 && a(n) == 1000;
  for (int i = 2; i <= n - 1; ++i) ok = ok && a(i) == i + 1000;
  check(ok, "parallel_for_n(n): every element 0..n-1 once, the last one included");

  // 4. serial: once
  serial([=] CPPMPAS_LAMBDA() { av(1) = av(1) + 1; });
  check(a(1) == 1001, "serial: the kernel runs once");

  // 5. kernel_fail with a value: the message, and nothing after the call runs (also after a failing callee)
  FArray<int> marker({8}, 0);
  const FView<int> mv = marker.view();
  const std::string msg = error_of([=] {
    parallel_for(1, 8, [=] CPPMPAS_LAMBDA(int i) {
      const int r = fails_deep(i);
      mv(i) = r;  // element 3: never written
    });
  });
  check(msg == "test_gpu_backend: element 3 failed", "kernel_fail(head, value, tail) is the launch's Error");
  check(marker(3) == 0, "no statement after a kernel_fail runs (the thread ends at the call)");
  check(marker(1) == 1 && marker(8) == 8, "the other threads of the failed launch ran");

  // 6. a work array beyond its capacity fails before it is used
  FArray<double> guard({4}, 0.0);
  const FView<double> gv = guard.view();
  const std::string cap = error_of([=] {
    serial([=] CPPMPAS_LAMBDA() {
      LocalArray<double, 4> w(16);
      for (int i = 0; i < 16; ++i) w[i] = 1.0;  // never reached
      gv(1) = w[0];
    });
  });
  check(cap == "a kernel work array is larger than its capacity", "a work array beyond its capacity: Error");
  check(guard(1) == 0.0, "a work array beyond its capacity is never written");

  // 7. the next launch runs normally
  parallel_for(1, n, [=] CPPMPAS_LAMBDA(int i) { av(i) = -i; });
  ok = true;
  for (int i = 1; i <= n; ++i) ok = ok && a(i) == -i;
  check(ok, "a launch after a failed one runs normally");

  // 8. a kernel with a small local frame (a work array indexed at run time) leaves the stack limit as it was or
  // raises it
  const int idx = n;  // not known to the compiler: the work array lives in local memory
  parallel_for(1, 4, [=] CPPMPAS_LAMBDA(int i) {
    LocalArray<double, 4> w(4);
    w[(i + idx) % 4] = double(i);
    av(i) = int(w[(i + idx) % 4]);
  });
  std::size_t stack1 = 0;
  check(cudaDeviceGetLimit(&stack1, cudaLimitStackSize) == cudaSuccess && stack1 >= stack0,
        "the device's stack limit is never lowered by a launch");
  std::printf("stack limit: %zu bytes before, %zu after\n", stack0, stack1);

  std::printf("%s\n", failures ? "FAILED" : "all checks passed");
  return failures ? 1 : 0;
}
