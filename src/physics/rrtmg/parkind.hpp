// SPDX-License-Identifier: GPL-3.0-or-later AND LicenseRef-AER-RRTMG
// Shared with CPPWRF (SHARED_WITH_CPPWRF.md): CPPWRF 7fac5bf src/phys/ra_rrtmg_common/parkind.hpp, copied unchanged except the namespace
// (cppmpas::phys::rrtmg), the CPPMPAS_ macro prefix, the include paths and the real type.
#pragma once
// RRTMG's real kind (module parkind of WRF's module_ra_rrtmg_lw.F, which the LW and SW modules use as
// "use parkind, only : im => kind_im, rb => kind_rb"). WRF builds it with kind_rb = kind(1.0), the default REAL:
// double in WRF -r8, float in WRF real*4. RRTMG's arithmetic here (LW, SW, the shared code in ra_rrtmg_common, the
// radiation driver's staging, Dudhia SW) is in rb, the model's real: float in the FP32 builds, double in the FP64
// builds, on the host and on the device. A constant written rb(c) is the double c rounded to rb (c itself in the
// FP64 builds).
#include "core/libm.hpp"  // CPPMPAS: libm:: of the kernels is cppmpas::libm
#include "core/precision.hpp"

// CPPMPAS addition: the load and warp-vote macros the RRTMG views and stages use, from CPPWRF 7fac5bf
// src/core/gpu_runtime.hpp:156-171 (renamed CPPWRF_ -> CPPMPAS_, logic unchanged)
#if defined(__CUDA_ARCH__)
#define CPPMPAS_LDG(ptr) __ldg(ptr)
#define CPPMPAS_WARP_ANY(pred) __any_sync(__activemask(), (pred))
#define CPPMPAS_WARP_ALL(pred) __all_sync(__activemask(), (pred))
#else
#define CPPMPAS_LDG(ptr) (*(ptr))
#define CPPMPAS_WARP_ANY(pred) (pred)
#define CPPMPAS_WARP_ALL(pred) (pred)
#endif

namespace cppmpas::phys::rrtmg {

using rb = ::cppmpas::real;  // kind_rb

}  // namespace cppmpas::phys::rrtmg
