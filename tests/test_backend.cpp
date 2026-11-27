// SPDX-License-Identifier: GPL-3.0-or-later
// The CPU build has none of the GPU's capacity limits (core/backend.hpp, operators/vector_operations.hpp): a kernel
// work array grows past its GPU capacity, a column of more than 256 levels is accepted, and mpas_legs /
// mpas_rbf_interp_func_3D_plane_vec_const_dir_comp_coeffs solve systems larger than the GPU's work arrays hold.
//
// usage: test_backend    exit 0 ok, 1 failure
#include <cmath>
#include <cstdio>
#include <vector>

#include "core/backend.hpp"
#include "operators/vector_operations.hpp"

using namespace cppmpas;

namespace {
int failures = 0;
void check(bool ok, const char* what) {
  std::printf("%s %s\n", ok ? "ok  " : "FAIL", what);
  if (!ok) ++failures;
}
}  // namespace

int main() {
  // a work array of 40 elements with a GPU capacity of 4: zero-filled, all writable, copied whole
  LocalArray<double, 4> w(40);
  bool ok = w.size() == 40;
  for (int i = 0; i < 40 && ok; ++i) ok = w[i] == 0.0;
  for (int i = 0; i < 40; ++i) w[i] = i;
  const LocalArray<double, 4> c = w;
  for (int i = 0; i < 40; ++i) ok = ok && c[i] == i && w[i] == i;
  check(ok, "a work array larger than its GPU capacity (CPU build): zero-filled, writable, copied");
  LocalArray<int, 8> small(3);
  small[2] = 7;
  const LocalArray<int, 8> sc = small;
  check(sc.size() == 3 && sc[0] == 0 && sc[2] == 7, "a work array within its capacity: copied");

  // a column of 1000 levels
  bool threw = false;
  try {
    fits_column(1000, "test_backend");
  } catch (const Error&) {
    threw = true;
  }
  check(!threw, "fits_column: 1000 levels accepted (no GPU column limit on the CPU)");

  // mpas_legs with N = 24 > kMaxRbfPoints + 2: A = 2 I, B = 2 x, solved exactly
  constexpr int N = 24;
  std::vector<real> A(N * N, 0.0), B(N), X(N, -1.0);
  std::vector<int> indx(N);
  for (int i = 0; i < N; ++i) {
    A[i + N * i] = 2.0;
    B[i] = 2.0 * (i + 1);
  }
  threw = false;
  try {
    mpas_legs(A.data(), N, B.data(), X.data(), indx.data());
  } catch (const Error&) {
    threw = true;
  }
  ok = !threw;
  for (int i = 0; i < N && ok; ++i) ok = X[i] == i + 1;
  check(ok, "mpas_legs: N = 24 solved (no GPU work-array limit on the CPU)");

  // the RBF coefficients of 20 points on a circle (more than kMaxRbfPoints)
  constexpr int pc = 20;
  std::vector<real> src(pc * 3), uv(pc * 3), coef(pc * 3, -7.0);
  for (int i = 0; i < pc; ++i) {
    const real t = 2.0 * 3.14159265358979 * i / pc;
    src[i] = std::cos(t);
    src[i + pc] = std::sin(t);
    src[i + 2 * pc] = 1.0;
    uv[i] = -std::sin(t);
    uv[i + pc] = std::cos(t);
    uv[i + 2 * pc] = 0.0;
  }
  const real dest[3] = {0.0, 0.0, 1.0};
  const real basis[6] = {1.0, 0.0, 0.0, 1.0, 0.0, 0.0};  // (2,3): rows (1,0,0) and (0,1,0)
  threw = false;
  try {
    mpas_rbf_interp_func_3D_plane_vec_const_dir_comp_coeffs(pc, src.data(), uv.data(), dest, 1.0, basis, coef.data());
  } catch (const Error&) {
    threw = true;
  }
  ok = !threw;
  for (int i = 0; i < pc * 3 && ok; ++i) ok = std::isfinite(coef[i]) && coef[i] != -7.0;
  check(ok, "RBF coefficients of 20 points computed (no GPU point limit on the CPU)");

  std::printf("%s\n", failures ? "FAILED" : "all checks passed");
  return failures ? 1 : 0;
}
