// SPDX-License-Identifier: GPL-3.0-or-later AND LGPL-2.1-or-later
#pragma once
// The libm / libmvec entries the MPAS-Atmosphere reference calls, bit for bit.
//
// The reference (gfortran 14.2 -O3, no -march, glibc 2.41-12+deb13u5) calls libm through the PLT; on its CPU
// (AVX2 + FMA) glibc's ifuncs select the FMA variants (__sin_fma, __ieee754_pow_fma, ...), and in the loops gfortran
// vectorised it calls libmvec's 2-lane SSE4 functions (_ZGVbN2v_sin -> _ZGVbN2v_sin_sse4, ...). Which entry each MPAS
// statement calls was read from the reference executable (objdump -d of bin/atmosphere_model of the mpas-build
// asset 153c7206...; see the call sites named in src/core_atmosphere and src/operators). The functions come from
// CPPWRF's ports (core/glibc_math.hpp, shared with CPPWRF); sincos, which CPPWRF does not use, is ported here.
#include <cstdint>

#include "core/glibc_math.hpp"

namespace cppmpas::libm {

// scalar libm (FMA variants)
inline double sin(double x) { return ::cppwrf::glibc::sin(x); }        // __sin_fma
inline double cos(double x) { return ::cppwrf::glibc::cos(x); }        // __cos_fma
inline double asin(double x) { return ::cppwrf::glibc::asin(x); }      // __ieee754_asin_fma
inline double acos(double x) { return ::cppwrf::glibc::acos(x); }      // __ieee754_acos_fma
inline double atan2(double y, double x) { return ::cppwrf::glibc::atan2(y, x); }  // __ieee754_atan2_fma
inline double pow(double x, double y) { return ::cppwrf::glibc::pow(x, y); }      // __ieee754_pow_fma

// one lane of libmvec's 2-lane SSE4 functions (glibc 2.41 sysdeps/x86_64/fpu/multiarch/svml_d_*2_core_sse4.S)
inline double mvec_sin(double x) { return ::cppwrf::glibc::mvec::sin(x); }               // _ZGVbN2v_sin_sse4
inline double mvec_cos(double x) { return ::cppwrf::glibc::mvec::cos(x); }               // _ZGVbN2v_cos_sse4
inline double mvec_pow(double x, double y) { return ::cppwrf::glibc::mvec::pow(x, y); }  // _ZGVbN2vv_pow_sse4

// sincos: __sincos_fma, glibc 2.41 sysdeps/ieee754/dbl-64/s_sincos.c (built with FMA through
// sysdeps/x86_64/fpu/multiarch/s_sincos-fma.c), statement by statement, with s_sin.c's do_sin, do_cos and
// reduce_sincos (cppwrf::glibc::detail, the bodies of __sin_fma / __cos_fma) and __branred (branred.c).
// Source: glibc-2.41.tar.xz, sha256 a5a26b22f545d6b7d7b3dd828e11e428f24f4fac43c934fb071b6a7d0828e901
// (s_sincos.c sha256 c3e3cdbc..., s_sin.c e708cfa5...). The rounding-mode switch of SET_RESTORE_ROUND_53BIT and
// errno are left out (round to nearest throughout), as in CPPWRF's sin / cos.
// Unlike sin, sincos takes the 0.855469 <= |x| < 2.426265 branch through a + da = (hp0 - |x|) + hp1 split exactly
// (s_sincos.c:63-66), so its sine there is not always __sin's.
inline void sincos(double x, double* sinx, double* cosx) {
  namespace g = ::cppwrf::glibc::detail;
  constexpr double hp0 = 0x1.921fb54442d18p0, hp1 = 0x1.1a62633145c07p-54;  // usncs.h
  const std::uint32_t k = std::uint32_t(g::asu(x) >> 32) & 0x7fffffffu;    // s_sincos.c:44-45
  if (k < 0x400368fdu) {                                                    // :47
    if (k < 0x3e400000u) {                                                  // :51  |x| < 2^-27
      *sinx = x;
      *cosx = 1.0;
      return;
    }
    if (k < 0x3feb6000u) {  // :59  |x| < 0.855469
      *sinx = g::do_sin(x, 0);
      *cosx = g::do_cos(x, 0);
      return;
    }
    const double y = hp0 - g::fabs_(x);  // :66-71  |x| < 2.426265
    const double a = y + hp1;
    const double da = (y - a) + hp1;
    *sinx = g::copysign(g::do_cos(a, da), x);
    *cosx = g::do_sin(a, da);
    return;
  }
  if (k < 0x7ff00000u) {  // :74  |x| < 2^1024
    double a, da;
    unsigned n = k < 0x419921FBu ? unsigned(g::reduce_sincos(x, &a, &da)) : unsigned(g::branred(x, &a, &da));  // :80
    n = n & 3;
    if (n == 1 || n == 2) {  // :83-87
      a = -a;
      da = -da;
    }
    double* s = sinx;
    double* c = cosx;
    if (n & 1) {  // :89-94
      double* t = c;
      c = s;
      s = t;
    }
    *s = g::do_sin(a, da);  // :96-98
    const double xx = g::do_cos(a, da);
    *c = (n & 2) ? -xx : xx;
    return;
  }
  *sinx = *cosx = x / x;  // :104
}

}  // namespace cppmpas::libm
