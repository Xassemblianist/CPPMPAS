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

#include "core/backend.hpp"
#include "core/glibc_math.hpp"

namespace cppmpas::libm {

// The argument as the library receives it: a value in a register. glibc's NaN results of a NaN argument are that
// NaN quieted (s_sin.c x / x, e_asin.c (x - x) / (x - x), e_atan2.c / e_pow.c x + y), sign included. Once a port
// is inlined, GCC may fold the caller's negation into those expressions ((-y) / (-y) -> y / y, (-y) + x -> x - y),
// which is exact for every non-NaN value but flips the sign of a NaN result. The empty asm hides the argument's
// origin, so a NaN argument gives glibc's bits (tests/test_libm.cpp calls the ports with -x).
CPPMPAS_FORCEINLINE double arg(double x) {
#if (defined(__x86_64__) || defined(__i386__)) && !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__)
  asm("" : "+x"(x));
#endif
  return x;
}

// scalar libm (FMA variants); every entry is CPPMPAS_FORCEINLINE (host and device), as the CPPWRF ports are
CPPMPAS_FORCEINLINE double sin(double x) { return ::cppwrf::glibc::sin(arg(x)); }  // __sin_fma
CPPMPAS_FORCEINLINE double cos(double x) { return ::cppwrf::glibc::cos(arg(x)); }  // __cos_fma
// asin / acos are the wrappers of math/w_asin_compat.c:30-38 and w_acos_compat.c:30-38 (the default asin / acos of the
// x86-64 libm, LIBM_SVID_COMPAT): |x| > 1 returns __kernel_standard's NAN (+qNaN), not __ieee754_asin's 0/0
CPPMPAS_FORCEINLINE double asin(double x) {  // __ieee754_asin_fma
  x = arg(x);
  return ::cppwrf::glibc::detail::fabs_(x) > 1.0 ? __builtin_nan("") : ::cppwrf::glibc::asin(x);
}
CPPMPAS_FORCEINLINE double acos(double x) {  // __ieee754_acos_fma
  x = arg(x);
  return ::cppwrf::glibc::detail::fabs_(x) > 1.0 ? __builtin_nan("") : ::cppwrf::glibc::acos(x);
}
CPPMPAS_FORCEINLINE double atan2(double y, double x) {  // __ieee754_atan2_fma
  return ::cppwrf::glibc::atan2(arg(y), arg(x));
}
CPPMPAS_FORCEINLINE double pow(double x, double y) { return ::cppwrf::glibc::pow(arg(x), arg(y)); }  // __ieee754_pow_fma

// one lane of libmvec's 2-lane SSE4 functions (glibc 2.41 sysdeps/x86_64/fpu/multiarch/svml_d_*2_core_sse4.S)
CPPMPAS_FORCEINLINE double mvec_sin(double x) { return ::cppwrf::glibc::mvec::sin(arg(x)); }  // _ZGVbN2v_sin_sse4
CPPMPAS_FORCEINLINE double mvec_cos(double x) { return ::cppwrf::glibc::mvec::cos(arg(x)); }  // _ZGVbN2v_cos_sse4
CPPMPAS_FORCEINLINE double mvec_pow(double x, double y) {  // _ZGVbN2vv_pow_sse4
  return ::cppwrf::glibc::mvec::pow(arg(x), arg(y));
}

// sincos: __sincos_fma, glibc 2.41 sysdeps/ieee754/dbl-64/s_sincos.c (built with FMA through
// sysdeps/x86_64/fpu/multiarch/s_sincos-fma.c), statement by statement, with s_sin.c's do_sin, do_cos and
// reduce_sincos (cppwrf::glibc::detail, the bodies of __sin_fma / __cos_fma) and __branred (branred.c).
// Source: glibc-2.41.tar.xz, sha256 a5a26b22f545d6b7d7b3dd828e11e428f24f4fac43c934fb071b6a7d0828e901
// (s_sincos.c sha256 c3e3cdbc..., s_sin.c e708cfa5...). The rounding-mode switch of SET_RESTORE_ROUND_53BIT and
// errno are left out (round to nearest throughout), as in CPPWRF's sin / cos.
// Unlike sin, sincos takes the 0.855469 <= |x| < 2.426265 branch through a + da = (hp0 - |x|) + hp1 split exactly
// (s_sincos.c:67-69), so its sine there is not always __sin's.
CPPMPAS_FORCEINLINE void sincos(double x, double* sinx, double* cosx) {
  namespace g = ::cppwrf::glibc::detail;
  x = arg(x);
  constexpr double hp0 = 0x1.921fb54442d18p0, hp1 = 0x1.1a62633145c07p-54;  // usncs.h
  const std::uint32_t k = std::uint32_t(g::asu(x) >> 32) & 0x7fffffffu;    // s_sincos.c:43-44
  if (k < 0x400368fdu) {                                                    // :46
    if (k < 0x3e400000u) {                                                  // :50  |x| < 2^-27
      *sinx = x;
      *cosx = 1.0;
      return;
    }
    if (k < 0x3feb6000u) {  // :59-64  |x| < 0.855469
      *sinx = g::do_sin(x, 0);
      *cosx = g::do_cos(x, 0);
      return;
    }
    const double y = hp0 - g::fabs_(x);  // :66-72  |x| < 2.426265
    const double a = y + hp1;
    const double da = (y - a) + hp1;
    *sinx = g::copysign(g::do_cos(a, da), x);
    *cosx = g::do_sin(a, da);
    return;
  }
  if (k < 0x7ff00000u) {  // :74-75  |x| < 2^1024
    double a, da;
    unsigned n = k < 0x419921FBu ? unsigned(g::reduce_sincos(x, &a, &da)) : unsigned(g::branred(x, &a, &da));  // :81
    n = n & 3;
    if (n == 1 || n == 2) {  // :84-88
      a = -a;
      da = -da;
    }
    double* s = sinx;
    double* c = cosx;
    if (n & 1) {  // :90-95
      double* t = c;
      c = s;
      s = t;
    }
    *s = g::do_sin(a, da);  // :97-99
    const double xx = g::do_cos(a, da);
    *c = (n & 2) ? -xx : xx;
    return;
  }
  *sinx = *cosx = x / x;  // :106 (errno of :103-104 left out)
}

}  // namespace cppmpas::libm
