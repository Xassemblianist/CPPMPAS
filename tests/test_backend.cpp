// SPDX-License-Identifier: GPL-3.0-or-later
// The capacity limits of core/backend.hpp and operators/vector_operations.hpp on the CPU build, which refuses the
// inputs the GPU's work arrays cannot hold, as the GPU build does (one range of inputs for both builds): a column of
// 256 levels is accepted and one of 257 refused (fits_column), mpas_legs solves N = kMaxRbfPoints + 2 = 18 and
// refuses 19, the RBF coefficients of kMaxRbfPoints = 16 points are computed and 17 points refused, 8 scalars are
// accepted by advance_scalars and 9 refused. A kernel work array still grows past its GPU capacity on the CPU (the
// heap; no caller asks for that within the checked ranges).
//
// usage: test_backend    exit 0 ok, 1 failure
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "core/backend.hpp"
#include "core_atmosphere/atm_advance_scalars.hpp"
#include "operators/vector_operations.hpp"

using namespace cppmpas;

namespace {
int failures = 0;
void check(bool ok, const char* what) {
  std::printf("%s %s\n", ok ? "ok  " : "FAIL", what);
  if (!ok) ++failures;
}
template <class F> std::string error_of(const F& f) {
  try {
    f();
  } catch (const Error& e) {
    return e.what();
  }
  return "";
}

// mpas_legs on A = 2 I, B = 2 x: X = x exactly; "" when solved, else the Error
std::string legs(int N, bool* solved) {
  std::vector<real> A(N * N, 0.0), B(N), X(N, -1.0);
  std::vector<int> indx(N);
  for (int i = 0; i < N; ++i) {
    A[i + N * i] = 2.0;
    B[i] = 2.0 * (i + 1);
  }
  const std::string e = error_of([&] { mpas_legs(A.data(), N, B.data(), X.data(), indx.data()); });
  *solved = e.empty();
  for (int i = 0; i < N && *solved; ++i) *solved = X[i] == i + 1;
  return e;
}

// the RBF coefficients of pc points on a circle; "" when computed (all finite and written), else the Error
std::string rbf(int pc, bool* computed) {
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
  const std::string e = error_of([&] {
    mpas_rbf_interp_func_3D_plane_vec_const_dir_comp_coeffs(pc, src.data(), uv.data(), dest, 1.0, basis, coef.data());
  });
  *computed = e.empty();
  for (int i = 0; i < pc * 3 && *computed; ++i) *computed = std::isfinite(coef[i]) && coef[i] != -7.0;
  return e;
}
}  // namespace

int main() {
  // a work array of 40 elements with a GPU capacity of 4: zero-filled, all writable, copied whole
  LocalArray<double, 4> w(40);
  bool ok = w.size() == 40;
  for (int i = 0; i < 40 && ok; ++i) ok = w[i] == 0.0 && !std::signbit(w[i]);  // +0.0, not -0.0
  for (int i = 0; i < 40; ++i) w[i] = i;
  const LocalArray<double, 4> c = w;
  for (int i = 0; i < 40; ++i) ok = ok && c[i] == i && w[i] == i;
  check(ok, "a work array larger than its GPU capacity (CPU build): zero-filled, writable, copied");
  LocalArray<int, 8> small(3);
  small[2] = 7;
  const LocalArray<int, 8> sc = small;
  check(sc.size() == 3 && sc[0] == 0 && sc[2] == 7, "a work array within its capacity: copied");

  // columns: 256 levels fit, 257 do not (the GPU's kMaxColumn)
  check(error_of([] { fits_column(256, "test_backend"); }).empty(), "fits_column: 256 levels accepted");
  const std::string col = error_of([] { fits_column(257, "test_backend"); });
  check(col.find("test_backend: nVertLevels = 257 is more than the 256 levels") == 0,
        "fits_column: 257 levels refused on the CPU as on the GPU");

  // mpas_legs: N = 18 solved, N = 19 refused
  bool solved = false;
  check(legs(18, &solved).empty() && solved, "mpas_legs: N = 18 solved");
  check(legs(19, &solved) == "mpas_legs: N > kMaxRbfPoints + 2 (the work array of elgs)", "mpas_legs: N = 19 refused");

  // RBF coefficients: 16 points computed, 17 refused
  bool computed = false;
  check(rbf(16, &computed).empty() && computed, "RBF coefficients of 16 points computed");
  check(rbf(17, &computed) == "mpas_rbf_interp_func_3D_plane_vec_const_dir_comp_coeffs: pointCount > kMaxRbfPoints",
        "RBF coefficients of 17 points refused");

  // scalars: 8 accepted, 9 refused
  check(error_of([] { check_scalar_capacity(8); }).empty(), "advance_scalars: 8 scalars accepted");
  check(error_of([] { check_scalar_capacity(9); }) ==
            "advance_scalars: more than 8 scalars (the capacity of the per-thread work arrays)",
        "advance_scalars: 9 scalars refused");

  std::printf("%s\n", failures ? "FAILED" : "all checks passed");
  return failures ? 1 : 0;
}
