// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// Working precision. MPAS: RKIND = selected_real_kind(12) unless SINGLE_PRECISION (src/framework/mpas_kind_types.F:24-28);
// the reference is the PRECISION=double build. Pattern from CPPWRF 7fac5bf src/core/precision.hpp (SHARED_WITH_CPPWRF.md).
namespace cppmpas {
#ifdef CPPMPAS_SINGLE_PRECISION
#error "CPPMPAS_SINGLE_PRECISION: the FP32 (MPAS real*4) build has no verified port yet"
#else
using real = double;
inline constexpr const char* precision_name = "float64";
#endif
}  // namespace cppmpas

// Every double operation on the host is rounded to its own type (SSE2 x86-64, aarch64): no x87 excess precision.
#if defined(__FLT_EVAL_METHOD__) && !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__)
static_assert(__FLT_EVAL_METHOD__ == 0, "CPPMPAS needs float and double arithmetic without excess precision");
#endif
