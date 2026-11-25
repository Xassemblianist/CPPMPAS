// SPDX-License-Identifier: GPL-3.0-or-later
// cppmpas::libm (src/core/libm.hpp) against this machine's glibc libm / libmvec, bit for bit. Meaningful only where
// the reference's libraries are: glibc 2.41 (Debian 13) on x86-64 with FMA and AVX2 (libm's FMA variants, libmvec's
// SSE4 2-lane variants); anywhere else the test fails rather than pass without checking.
//   test_libm [N]   N random inputs per function (default 2'000'000)
#include <emmintrin.h>
#include <gnu/libc-version.h>

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <string>

#include "core/libm.hpp"

extern "C" __m128d _ZGVbN2v_sin(__m128d);
extern "C" __m128d _ZGVbN2v_cos(__m128d);
extern "C" __m128d _ZGVbN2vv_pow(__m128d, __m128d);
extern "C" void sincos(double, double*, double*);

namespace {
double (*volatile c_sin)(double) = ::sin;
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
bool same(double a, double b) { return bits(a) == bits(b) || (std::isnan(a) && std::isnan(b)); }

struct Count {
  const char* name;
  long n = 0, bad = 0;
  void check(double x, double y, double got, double want) {
    ++n;
    if (!same(got, want)) {
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
  Count sc{"sincos.sin"}, cc{"sincos.cos"}, s{"sin"}, as{"asin"}, ac{"acos"}, at{"atan2"}, pw{"pow"},
      vs{"mvec_sin"}, vc{"mvec_cos"}, vp{"mvec_pow"};
  for (long i = 0; i < N; ++i) {
    // sincos over every branch of s_sincos.c (|x| from 2^-30 to 2^40) and the angle range MPAS passes (|x| < 2 pi)
    const double x = (i & 1) ? any_double(-30, 40) : (u01(rng) - 0.5) * 4.0 * M_PI;
    double gs, gc, ps, pc;
    c_sincos(x, &gs, &gc);
    cppmpas::libm::sincos(x, &ps, &pc);
    sc.check(x, 0, ps, gs);
    cc.check(x, 0, pc, gc);
    s.check(x, 0, cppmpas::libm::sin(x), c_sin(x));
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
    const double b2 = std::fabs(any_double(-12, 12));
    double outp[2];
    _mm_storeu_pd(outp, _ZGVbN2vv_pow(_mm_set_pd(b2, b), _mm_set_pd(e, e)));
    vp.check(b, e, cppmpas::libm::mvec_pow(b, e), outp[0]);
    vp.check(b2, e, cppmpas::libm::mvec_pow(b2, e), outp[1]);
  }
  bool ok = true;
  for (const Count* c : {&sc, &cc, &s, &as, &ac, &at, &pw, &vs, &vc, &vp}) ok = c->report() && ok;
  std::printf(ok ? "PASS\n" : "FAIL\n");
  return ok ? 0 : 1;
}
