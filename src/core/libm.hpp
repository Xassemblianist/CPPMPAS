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

// The lanes of libmvec's 2-lane log / log10 outside the normal range (x < DBL_MIN: subnormals, +-0, negatives;
// +-Inf, NaN) call the scalar log / log10 of libm (svml_d_log2_core_sse4.S, svml_d_log102_core_sse4.S), the
// reference's __log_fma / __ieee754_log10. CPPWRF's core/libmvec_x86.hpp takes them from the glibc ports of
// core/glibc_math.hpp with CPPWRF_GLIBC_MATH and from std::log / std::log10 without: on the device that is CUDA's
// log, which rounds differently (log10 of 0x1p-1023: ...5221 on the GPU against glibc's ...5220). CPPMPAS defines it
// for every source, host and device, so both builds compute those lanes as glibc does.
#ifndef CPPWRF_GLIBC_MATH
#define CPPWRF_GLIBC_MATH 1
#endif

#include "core/backend.hpp"
#include "core/glibc_math.hpp"
#include "core/libmvec_x86.hpp"

namespace cppmpas::libm {

// The argument as the library receives it: a value in a register. glibc's NaN results of a NaN argument are that
// NaN quieted (s_sin.c x / x, e_asin.c (x - x) / (x - x), e_atan2.c / e_pow.c x + y), sign included. Once a port
// is inlined, GCC may fold the caller's negation into those expressions ((-y) / (-y) -> y / y, (-y) + x -> x - y),
// which is exact for every non-NaN value but flips the sign of a NaN result. The empty asm hides the argument's
// origin, so a NaN argument gives glibc's bits (tests/test_libm.cpp calls the ports with -x).
// On the device (CUDA) the same barrier is an identity PTX move, opaque to nvcc's optimiser (tests/test_gpu_libm.cu
// measures what NaN arguments give there; CPPMPAS_NO_DEVICE_ARG_BARRIER leaves it out, for that measurement only).
CPPMPAS_FORCEINLINE double arg(double x) {
#if (defined(__x86_64__) || defined(__i386__)) && !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__)
  asm("" : "+x"(x));
#elif defined(__CUDA_ARCH__) && !defined(CPPMPAS_NO_DEVICE_ARG_BARRIER)
  asm("mov.b64 %0, %0;" : "+d"(x));
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
CPPMPAS_FORCEINLINE double exp(double x) { return ::cppwrf::glibc::exp(arg(x)); }  // __exp_fma (Kessler)

// one lane of libmvec's 2-lane SSE4 functions (glibc 2.41 sysdeps/x86_64/fpu/multiarch/svml_d_*2_core_sse4.S)
CPPMPAS_FORCEINLINE double mvec_sin(double x) { return ::cppwrf::glibc::mvec::sin(arg(x)); }  // _ZGVbN2v_sin_sse4
CPPMPAS_FORCEINLINE double mvec_cos(double x) { return ::cppwrf::glibc::mvec::cos(arg(x)); }  // _ZGVbN2v_cos_sse4
CPPMPAS_FORCEINLINE double mvec_pow(double x, double y) {  // _ZGVbN2vv_pow_sse4
  return ::cppwrf::glibc::mvec::pow(arg(x), arg(y));
}
CPPMPAS_FORCEINLINE double mvec_exp(double x) { return ::cppwrf::glibc::mvec::exp(arg(x)); }  // _ZGVbN2v_exp_sse4

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

// WSM6 (physics/mp_wsm6.hpp): log (__log_fma), atan (__atan_fma, mp_wsm6_init's pi) and the 2-lane log
// _ZGVbN2v_log_sse4 (CPPWRF's core/libmvec_x86.hpp, shared). The vector log takes an approximate reciprocal with
// rcpps, which the hardware defines: libmvec_tables.hpp holds Intel's (dumped on an i7-10750H, the workstation);
// lanes outside the normal range call the scalar log, as libmvec does.
CPPMPAS_FORCEINLINE double log(double x) { return ::cppwrf::glibc::log(arg(x)); }    // __log_fma
CPPMPAS_FORCEINLINE double atan(double x) { return ::cppwrf::glibc::atan(arg(x)); }  // __atan_fma
CPPMPAS_FORCEINLINE double mvec_log(double x) { return ::cppwrf::mvec::log(arg(x)); }  // _ZGVbN2v_log_sse4 (Intel rcpps)

// YSU (physics/bl_ysu.hpp): tanh, a plain FUNC of libm (no ifunc; libm 0x5e560, s_tanh.c calling __expm1),
// CPPWRF's port. The revised surface layer's tables (physics/sf_sfclayrev.hpp, the MPAS -O3 variant of
// sf_sfclayrev_init): the 2-lane atan _ZGVbN2v_atan_sse4.
CPPMPAS_FORCEINLINE double tanh(double x) { return ::cppwrf::glibc::tanh(arg(x)); }        // tanh (s_tanh.c)
CPPMPAS_FORCEINLINE double mvec_atan(double x) { return ::cppwrf::glibc::mvec::atan(arg(x)); }  // _ZGVbN2v_atan_sse4
// Noah LSM (physics/sf_noah.hpp, TDFCND): log10 is __ieee754_log10 (sysdeps/ieee754/dbl-64/e_log10.c, no FMA
// variant; its log is __log_fma). The reference calls it through log10@plt; no libmvec lane in Noah.
// log10@plt is the x86-64 libm's default log10, the wrapper of math/w_log10_compat.c:25-41 (LIBM_SVID_COMPAT), as
// asin / acos above: x < 0 returns __kernel_standard's NAN (+qNaN, k_standard.c type 19), not __ieee754_log10's
// (x - x) / (x - x); +-0 gives -HUGE_VAL (type 18), which __ieee754_log10 also returns. Any other x is
// __ieee754_log10's.
CPPMPAS_FORCEINLINE double log10(double x) {  // log10 (w_log10_compat.c) -> __ieee754_log10
  x = arg(x);
  return x < 0.0 ? __builtin_nan("") : ::cppwrf::glibc::log10(x);
}

// WSM6 reflectivity (physics/mp_radar.hpp, mp_wsm6_refl.hpp): the 2-lane log10 _ZGVbN2v_log10_sse4 (CPPWRF's
// core/libmvec_x86.hpp, Intel rcpps tables as the reference runs preload them), and the complex functions gfortran
// calls for COMPLEX*16 ABS, LOG and SQRT: cabs (__cabs = __hypot, no ifunc), clog (__clog), csqrt (__csqrt)
// (CPPWRF's core/glibc_math.hpp ports).
// A lane outside the normal range calls log10@plt (svml_d_log102_core_sse4.S): the wrapper of log10 above, so a
// negative lane is +qNaN; CPPWRF's lane takes __ieee754_log10 (core/glibc_math.hpp, CPPWRF_GLIBC_MATH) for the rest.
CPPMPAS_FORCEINLINE double mvec_log10(double x) {  // _ZGVbN2v_log10_sse4
  x = arg(x);
  return x < 0.0 ? __builtin_nan("") : ::cppwrf::mvec::log10(x);
}
CPPMPAS_FORCEINLINE double cabs(double re, double im) { return ::cppwrf::glibc::hypot(arg(re), arg(im)); }
CPPMPAS_FORCEINLINE void clog(double re, double im, double* r, double* i) {
  ::cppwrf::glibc::clog(arg(re), arg(im), r, i);
}
CPPMPAS_FORCEINLINE void csqrt(double re, double im, double* r, double* i) {
  ::cppwrf::glibc::csqrt(arg(re), arg(im), r, i);
}

}  // namespace cppmpas::libm
