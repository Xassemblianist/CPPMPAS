// SPDX-License-Identifier: GPL-3.0-or-later
// cppmpas::libm (src/core/libm.hpp) against this machine's glibc libm / libmvec, bit for bit. Meaningful only where
// the reference's libraries are: glibc 2.41 (Debian 13) on x86-64 with FMA and AVX2 (libm's FMA variants, libmvec's
// SSE4 2-lane variants); anywhere else the test fails rather than pass without checking.
// A fixed set of special inputs comes first: +-0, subnormals, +-DBL_MIN / DBL_MAX, +-inf, NaN, +-4 ulps around the
// branch thresholds of s_sincos.c / s_sin.c / e_asin.c, |x| up to 2^1023 (__branred), the sincos strip
// hp0 -+ (2m + 1) / 256 (table-bucket midpoints of y = hp0 - |x|, where do_cos(a, da) and do_cos(y, hp1) differ),
// atan2 / pow on every pair of special values, and libmvec pow with distinct lanes. Results must have the same
// bits, NaN results of non-NaN inputs included (asin / acos of |x| > 1).
//   test_libm [N]   N random inputs per function (default 2'000'000)
#include <emmintrin.h>
#include <gnu/libc-version.h>

#include <cfloat>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <random>
#include <string>
#include <vector>

#include "core/libm.hpp"

extern "C" __m128d _ZGVbN2v_sin(__m128d);
extern "C" __m128d _ZGVbN2v_cos(__m128d);
extern "C" __m128d _ZGVbN2vv_pow(__m128d, __m128d);
extern "C" void sincos(double, double*, double*);

namespace {
double (*volatile c_sin)(double) = ::sin;
double (*volatile c_cos)(double) = ::cos;
double (*volatile c_asin)(double) = ::asin;
double (*volatile c_acos)(double) = ::acos;
double (*volatile c_atan2)(double, double) = ::atan2;
double (*volatile c_pow)(double, double) = ::pow;
void (*volatile c_sincos)(double, double*, double*) = ::sincos;

std::uint64_t bits(double x) {
  std::uint64_t b;
  std::memcpy(&b, &x, 8);
  return b;
}
bool same(double a, double b) { return bits(a) == bits(b); }
double from_bits(std::uint64_t b) {
  double x;
  std::memcpy(&x, &b, 8);
  return x;
}

struct Count {
  const char* name;
  long n = 0, bad = 0;
  void check(double x, double y, double got, double want) {
    ++n;
    // a NaN input propagates as a NaN whose sign the compiler need not keep (GCC folds (-x)/(-x) to x/x without
    // -fsignaling-nans): only NaN-ness is compared then; every other result must have the same bits
    const bool nan_in = std::isnan(x) || std::isnan(y);
    if (nan_in ? std::isnan(got) != std::isnan(want) || (!std::isnan(want) && !same(got, want)) : !same(got, want)) {
      if (bad < 5)
        std::printf("  %s(%a, %a): cppmpas %a (0x%016llx) glibc %a (0x%016llx)\n", name, x, y, got,
                    (unsigned long long)bits(got), want, (unsigned long long)bits(want));
      ++bad;
    }
  }
  bool report() const {
    std::printf("%-10s %ld inputs, %ld differ\n", name, n, bad);
    return bad == 0 && n > 0;
  }
};
}  // namespace

int main(int argc, char** argv) {
  const long N = argc > 1 ? std::atol(argv[1]) : 2000000;
  if (std::string(gnu_get_libc_version()) != "2.41") {
    std::printf("FAIL: glibc %s, the reference uses 2.41\n", gnu_get_libc_version());
    return 1;
  }
  __builtin_cpu_init();
  if (!__builtin_cpu_supports("fma") || !__builtin_cpu_supports("avx2") || !__builtin_cpu_supports("sse4.1")) {
    std::printf("FAIL: this CPU lacks FMA / AVX2 / SSE4.1: glibc would not select the reference's variants\n");
    return 1;
  }
  std::mt19937_64 rng(20261126);
  std::uniform_real_distribution<double> u01(0.0, 1.0);
  auto any_double = [&](int emin, int emax) {  // random sign, mantissa, exponent in [emin, emax]
    const double m = 1.0 + u01(rng);
    const int e = emin + int(rng() % std::uint64_t(emax - emin + 1));
    return ((rng() & 1) ? -1.0 : 1.0) * std::ldexp(m, e - 1);
  };
  Count sc{"sincos.sin"}, cc{"sincos.cos"}, s{"sin"}, co{"cos"}, as{"asin"}, ac{"acos"}, at{"atan2"}, pw{"pow"},
      vs{"mvec_sin"}, vc{"mvec_cos"}, vp{"mvec_pow"};
  auto unary = [&](double x) {
    double gs, gc, ps, pc;
    c_sincos(x, &gs, &gc);
    cppmpas::libm::sincos(x, &ps, &pc);
    sc.check(x, 0, ps, gs);
    cc.check(x, 0, pc, gc);
    s.check(x, 0, cppmpas::libm::sin(x), c_sin(x));
    co.check(x, 0, cppmpas::libm::cos(x), c_cos(x));
    as.check(x, 0, cppmpas::libm::asin(x), c_asin(x));
    ac.check(x, 0, cppmpas::libm::acos(x), c_acos(x));
    const __m128d vin = _mm_set_pd(-x, x);
    double outs[2], outc[2];
    _mm_storeu_pd(outs, _ZGVbN2v_sin(vin));
    _mm_storeu_pd(outc, _ZGVbN2v_cos(vin));
    vs.check(x, 0, cppmpas::libm::mvec_sin(x), outs[0]);
    vs.check(-x, 0, cppmpas::libm::mvec_sin(-x), outs[1]);
    vc.check(x, 0, cppmpas::libm::mvec_cos(x), outc[0]);
    vc.check(-x, 0, cppmpas::libm::mvec_cos(-x), outc[1]);
  };
  auto mvec_pow_pair = [&](double b1, double e1, double b2, double e2) {
    double outp[2];
    _mm_storeu_pd(outp, _ZGVbN2vv_pow(_mm_set_pd(b2, b1), _mm_set_pd(e2, e1)));
    vp.check(b1, e1, cppmpas::libm::mvec_pow(b1, e1), outp[0]);
    vp.check(b2, e2, cppmpas::libm::mvec_pow(b2, e2), outp[1]);
  };
  auto binary = [&](double y, double x) {
    at.check(y, x, cppmpas::libm::atan2(y, x), c_atan2(y, x));
    pw.check(y, x, cppmpas::libm::pow(y, x), c_pow(y, x));
  };

  // --- special inputs
  const double inf = HUGE_VAL, qnan = std::numeric_limits<double>::quiet_NaN();
  std::vector<double> sp = {0.0, DBL_TRUE_MIN, 2 * DBL_TRUE_MIN, DBL_MIN / 3, DBL_MIN, DBL_MAX, inf, qnan, 1.0, 0.5, 2.0,
                            M_PI / 2, M_PI, 2 * M_PI, 287.0 / 717.5, 0.25, 0.75, 3.0, 1e300, 1e-300, 1024.0, 1075.0};
  // high words of the branch thresholds: s_sincos.c / s_sin.c (2^-27, 2^-26, 0.855469, 2.426265, 105414350,
  // 2^1024), e_asin.c / e_acos (2^-26 .. 1), and powers of two up to 2^1023
  std::vector<std::uint32_t> hi = {0x3e400000u, 0x3e500000u, 0x3feb6000u, 0x400368fdu, 0x419921FBu, 0x7fefffffu,
                                   0x3c900000u, 0x3fc00000u, 0x3fe00000u, 0x3fe40000u, 0x3fe80000u, 0x3fec0000u,
                                   0x3fee8000u, 0x3fef0000u, 0x3ff00000u, 0x3ff00001u};
  for (std::uint32_t e = 0x3ff00000u; e <= 0x7fe00000u; e += 0x00100000u) hi.push_back(e);
  std::vector<double> xs;
  for (double x : sp) xs.push_back(x);
  for (std::uint32_t h : hi) {
    double x = from_bits(std::uint64_t(h) << 32);
    for (int k = 0; k < 4; ++k) x = std::nextafter(x, 0.0);
    for (int k = -4; k <= 4; ++k, x = std::nextafter(x, inf)) xs.push_back(x);
  }
  // asin / acos at +-1, +-(1 - 2^-53), +-1.0001
  xs.push_back(1.0 - 0x1p-53);
  xs.push_back(1.0001);
  // the sincos strip: y = hp0 - |x| = +-(2m + 1) / 256 exactly (hp0 - c is exact for these c), then +-4 ulps of x
  const double hp0 = 0x1.921fb54442d18p0;
  for (int m = 0; m <= 109; ++m) {
    const double c = (2 * m + 1) / 256.0;
    for (double x0 : {hp0 - c, hp0 + c}) {
      if (x0 < 0.855469 || x0 >= 2.426265) continue;
      double x = x0;
      for (int k = 0; k < 4; ++k) x = std::nextafter(x, 0.0);
      for (int k = -4; k <= 4; ++k, x = std::nextafter(x, inf)) xs.push_back(x);
    }
  }
  long nspecial = 0;
  for (double x : xs)
    for (double sx : {x, -x}) {
      unary(sx);
      ++nspecial;
    }
  // atan2 and pow on every pair of the special values (axis-aligned atan2 arguments included)
  for (double y : sp)
    for (double x : sp)
      for (double sy : {y, -y})
        for (double sx : {x, -x}) {
          binary(sy, sx);
          mvec_pow_pair(sy, sx, sx, sy);
        }
  std::printf("special inputs: %ld unary, %zu binary pairs\n", nspecial, 4 * sp.size() * sp.size());
  for (long i = 0; i < N; ++i) {
    // sincos over every branch of s_sincos.c (|x| from 2^-30 to 2^40) and the angle range MPAS passes (|x| < 2 pi)
    const double x = (i & 1) ? any_double(-30, 40) : (u01(rng) - 0.5) * 4.0 * M_PI;
    double gs, gc, ps, pc;
    c_sincos(x, &gs, &gc);
    cppmpas::libm::sincos(x, &ps, &pc);
    sc.check(x, 0, ps, gs);
    cc.check(x, 0, pc, gc);
    s.check(x, 0, cppmpas::libm::sin(x), c_sin(x));
    co.check(x, 0, cppmpas::libm::cos(x), c_cos(x));
    const double a = (i & 1) ? any_double(-30, 0) : 2.0 * u01(rng) - 1.0;
    as.check(a, 0, cppmpas::libm::asin(a), c_asin(a));
    ac.check(a, 0, cppmpas::libm::acos(a), c_acos(a));
    const double yy = any_double(-20, 20), xx = any_double(-20, 20);
    at.check(yy, xx, cppmpas::libm::atan2(yy, xx), c_atan2(yy, xx));
    // pow: positive bases as in MPAS (exner argument ~0.1 .. 3, meshDensity) with MPAS's exponents and random ones
    const double b = std::fabs(any_double(-12, 12));
    const double e = (i % 4 == 0) ? 287.0 / 717.5 : (i % 4 == 1) ? 0.25 : (i % 4 == 2) ? 0.75 : any_double(-4, 4);
    pw.check(b, e, cppmpas::libm::pow(b, e), c_pow(b, e));
    // libmvec lanes: two inputs per call
    const double v1 = (u01(rng) - 0.5) * 8.0 * M_PI, v2 = any_double(-30, 30);
    const __m128d vin = _mm_set_pd(v2, v1);
    double outs[2], outc[2];
    _mm_storeu_pd(outs, _ZGVbN2v_sin(vin));
    _mm_storeu_pd(outc, _ZGVbN2v_cos(vin));
    vs.check(v1, 0, cppmpas::libm::mvec_sin(v1), outs[0]);
    vs.check(v2, 0, cppmpas::libm::mvec_sin(v2), outs[1]);
    vc.check(v1, 0, cppmpas::libm::mvec_cos(v1), outc[0]);
    vc.check(v2, 0, cppmpas::libm::mvec_cos(v2), outc[1]);
    // libmvec pow: the same exponent in both lanes (as MPAS calls it), then distinct exponents
    const double b2 = std::fabs(any_double(-12, 12));
    mvec_pow_pair(b, e, b2, (i & 1) ? e : any_double(-4, 4));
  }
  bool ok = true;
  for (const Count* c : {&sc, &cc, &s, &co, &as, &ac, &at, &pw, &vs, &vc, &vp}) ok = c->report() && ok;
  std::printf(ok ? "PASS\n" : "FAIL\n");
  return ok ? 0 : 1;
}
