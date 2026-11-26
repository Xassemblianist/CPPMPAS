// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// Kessler warm-rain microphysics, the column kernel of CPPWRF (commit 7fac5bf, src/phys/mp_kessler/kessler.hpp,
// last changed in e6892bb), a statement-by-statement port of WRF's phys/module_mp_kessler.F. Shared with CPPWRF
// (CONSTITUTION.md 2.6, SHARED_WITH_CPPWRF.md): the struct KesslerKernel below is CPPWRF's, copied unchanged except
//  * the MPAS-specific parameter es_cap: MPAS v8.4.2's copy of the scheme
//    (src/core_atmosphere/physics/physics_wrf/module_mp_kessler.F, 236 lines) is WRF's plus one line, 219
//    "es = min(es, 0.99*pressure)" (diff against WRF 4.8.0 phys/module_mp_kessler.F: only that line); MPAS calls
//    the kernel with es_cap = true;
//  * CPPWRF_HD spelled CPPMPAS_HD.
// The shims above it map CPPWRF's names to CPPMPAS's: View3 / View2 (the same strided views), rpow (scalar pow),
// vm::pow / vm::exp (lane(i, n): libmvec's 2-lane function for the pairs of iterations from the start of the loop,
// the scalar libm for the odd last iteration and for a 1-trip loop, CPPWRF core/backend.hpp:197).
//
// The MPAS reference executable (bin/atmosphere_model of the mpas-build asset 153c7206...,
// __module_mp_kessler_MOD_kessler at 0x355040) was read for the lanes, as MPAS builds the scheme with -O3 where
// WRF uses -O2 -ftree-vectorize:
//  * terminal velocity (module_mp_kessler.F:102-110): the vectorised copy (_ZGVbN2vv_pow at 0x355ab1) is guarded
//    by the overlap test of rhok(1) with the stores of rhok(k:k+1) (0x355901: min(kts+1, kte) offset > kts
//    offset), which always fails for kte >= 2, so every level calls pow@plt (0x357224): rpow;
//  * split sedimentation steps (:166-171): k pairs from kts through _ZGVbN2vv_pow (0x356227), an odd last level
//    pow@plt (0x356334): vm::pow(.., k, nz);
//  * production / evaporation (:203-236): the i loop in pairs from its with _ZGVbN2vv_pow (0x356787, 0x356871,
//    0x3569f2, 0x356a20) and _ZGVbN2v_exp (0x3568dc), an odd last i with pow@plt / exp@plt: vm::pow / vm::exp
//    (.., li, ni).
// No FMA, no ymm register. NINT is lround@plt (0x355e36, 0x3563d9).
#include <algorithm>
#include <cmath>

#include "core/backend.hpp"
#include "core/libm.hpp"
#include "core/precision.hpp"

namespace cppmpas::phys::kessler {

template <class T> struct View3 {
  T* p = nullptr;
  int sj = 0, sk = 0;
  CPPMPAS_HD T& operator()(int i, int j, int k) const { return p[i + j * sj + k * sk]; }
};
template <class T> struct View2 {
  T* p = nullptr;
  int sj = 0;
  CPPMPAS_HD T& operator()(int i, int j) const { return p[i + j * sj]; }
};
CPPMPAS_HD inline double rpow(double a, double b) { return ::cppmpas::libm::pow(a, b); }
namespace vm {
CPPMPAS_HD inline bool lane(int i, int n) { return n > 1 && !((n & 1) && i == n - 1); }
CPPMPAS_HD inline double pow(double a, double b, int i, int n) {
  return lane(i, n) ? ::cppmpas::libm::mvec_pow(a, b) : ::cppmpas::libm::pow(a, b);
}
CPPMPAS_HD inline double exp(double a, int i, int n) {
  return lane(i, n) ? ::cppmpas::libm::mvec_exp(a) : ::cppmpas::libm::exp(a);
}
}  // namespace vm

static constexpr int kMaxLevels = 128;

struct KesslerKernel {
  View3<real> t;
  View3<real> qv;
  View3<real> qc;
  View3<real> qr;
  View3<const real> rho;
  View3<const real> pii;
  View3<const real> dz8w;
  View3<const real> z;
  View2<real> rainnc;
  View2<real> rainncv;
  real dt;
  real xlv;
  real cp;
  real ep2;
  real svp1;
  real svp2;
  real svp3;
  real svpt0;
  real rhowater;
  int nz;
  // WRF's i loop of the microphysics part (its..ite = i0..i0+ni-1): gfortran vectorises it through
  // libmvec (cppwrf::vm lanes). Of the fall-speed k loops, the first stays scalar in WRF's object (its
  // vectorised version is guarded by an overlap check of rhok(1) with the rhok(k) stores, which always
  // fails) and the one of the split sedimentation steps (nfall > 1) is vectorised over k = kts..kte
  int i0 = 0, ni = 1;
  // MPAS v8.4.2 physics_wrf/module_mp_kessler.F:219 "es = min(es, 0.99*pressure)", a line WRF 4.8.0's
  // phys/module_mp_kessler.F does not have: true for MPAS, false gives the CPPWRF kernel
  bool es_cap = false;

  CPPMPAS_HD void operator()(int i, int j) const {
    using std::max;
    using std::min;
    using std::sqrt;
    using std::round;

    // Fixed-size local arrays per column (Fortran module_mp_kessler.F lines 65)
    real vt[kMaxLevels];
    real prodk[kMaxLevels];
    real vtden[kMaxLevels];
    real rdzk[kMaxLevels];
    real rhok[kMaxLevels];
    real factor[kMaxLevels];
    real rdzw[kMaxLevels];
    real prod[kMaxLevels];

    // Fortran line 81: f5 = svp2*(svpt0-svp3)*xlv/cp
    const real f5 = svp2 * (svpt0 - svp3) * xlv / cp;

    // Fortran line 23-28 parameters:
    constexpr real c1 = real(0.001);
    constexpr real c2 = real(0.001);
    constexpr real c3 = real(2.2);
    constexpr real c4 = real(0.875);
    constexpr real max_cr_sedimentation = real(0.75);

    // Fortran lines 102-124: Terminal velocity calculation and advection setup
    real crmax = real(0);
    for (int k = 0; k < nz; ++k) {
      prodk[k] = qr(i, j, k);
      rhok[k] = rho(i, j, k);
      real qrr = max(real(0), qr(i, j, k) * real(0.001) * rhok[k]);
      vtden[k] = sqrt(rhok[0] / rhok[k]);
      vt[k] = real(36.34) * rpow(qrr, real(0.1364)) * vtden[k];
      rdzw[k] = real(1) / dz8w(i, j, k);
      crmax = max(vt[k] * dt * rdzw[k], crmax);
    }
    for (int k = 0; k < nz - 1; ++k) {
      rdzk[k] = real(1) / (z(i, j, k + 1) - z(i, j, k));
    }
    rdzk[nz - 1] = real(1) / (z(i, j, nz - 1) - z(i, j, nz - 2));

    // Fortran lines 125-128: Courant number and sedimentation timestep splitting
    int nfall = max(1, int(round(real(0.5) + crmax / max_cr_sedimentation)));
    real dtfall = dt / real(nfall);
    real time_sediment = dt;

    // Fortran lines 134-193: Column sedimentation loop
    while (nfall > 0) {
      time_sediment = time_sediment - dtfall;
      for (int k = 0; k < nz - 1; ++k) {
        factor[k] = dtfall * rdzk[k] / rhok[k];
      }
      factor[nz - 1] = dtfall * rdzk[nz - 1];

      // Surface rain fallout
      real ppt = rhok[0] * prodk[0] * vt[0] * dtfall / rhowater;
      rainncv(i, j) = ppt * real(1000.0);
      rainnc(i, j) = rainnc(i, j) + ppt * real(1000.0);

      // Fallout flux upstream
      for (int k = 0; k < nz - 1; ++k) {
        prodk[k] = prodk[k] - factor[k] * (rhok[k] * prodk[k] * vt[k] - rhok[k + 1] * prodk[k + 1] * vt[k + 1]);
      }
      prodk[nz - 1] = prodk[nz - 1] - factor[nz - 1] * prodk[nz - 1] * vt[nz - 1];

      if (nfall > 1) {
        nfall = nfall - 1;
        crmax = real(0);
        for (int k = 0; k < nz; ++k) {
          real qrr = max(real(0), prodk[k] * real(0.001) * rhok[k]);
          vt[k] = real(36.34) * vm::pow(qrr, real(0.1364), k, nz) * vtden[k];
          crmax = max(vt[k] * time_sediment * rdzw[k], crmax);
        }
        int nfall_new = max(1, int(round(real(0.5) + crmax / max_cr_sedimentation)));
        if (nfall_new != nfall) {
          nfall = nfall_new;
          dtfall = time_sediment / real(nfall);
        }
      } else {
        for (int k = 0; k < nz; ++k) {
          prod[k] = prodk[k];
        }
        nfall = 0;
      }
    }

    // Fortran lines 203-236: Microphysics production & evaporation
    const int li = i - i0;
    for (int k = 0; k < nz; ++k) {
      real factorn = real(1.0) / (real(1.0) + c3 * dt * vm::pow(max(real(0.0), qr(i, j, k)), c4, li, ni));
      real qrprod = qc(i, j, k) * (real(1.0) - factorn) + factorn * c1 * dt * max(qc(i, j, k) - c2, real(0.0));
      real rcgs = real(0.001) * rho(i, j, k);

      qc(i, j, k) = max(qc(i, j, k) - qrprod, real(0.0));
      qr(i, j, k) = (qr(i, j, k) + prod[k] - qr(i, j, k));
      qr(i, j, k) = max(qr(i, j, k) + qrprod, real(0.0));

      real temp = pii(i, j, k) * t(i, j, k);
      real pressure = real(1.000e+05) * vm::pow(pii(i, j, k), real(1004.0 / 287.0), li, ni);
      real gam = real(2.5e+06) / (real(1004.0) * pii(i, j, k));
      real es = real(1000.0) * svp1 * vm::exp(svp2 * (temp - svpt0) / (temp - svp3), li, ni);
      if (es_cap) es = min(es, real(0.99) * pressure);  // MPAS-specific parameter (see the header)
      real qvs = ep2 * es / (pressure - es);
      prod[k] = (qv(i, j, k) - qvs) / (real(1.0) + pressure / (pressure - es) * qvs * f5 / ((temp - svp3) * (temp - svp3)));
      real dim_qvs_qv = max(real(0.0), qvs - qv(i, j, k));
      real ern_num = (real(1.6) + real(124.9) * vm::pow(rcgs * qr(i, j, k), real(0.2046), li, ni)) * vm::pow(rcgs * qr(i, j, k), real(0.525), li, ni);
      real ern_den = real(2.55e8) / (pressure * qvs) + real(5.4e5);
      real ern_term1 = dt * (ern_num / ern_den) * (dim_qvs_qv / (rcgs * qvs));
      real ern_term2 = max(-prod[k] - qc(i, j, k), real(0.0));
      real ern = min(ern_term1, min(ern_term2, qr(i, j, k)));

      real product = max(prod[k], -qc(i, j, k));
      t(i, j, k) = t(i, j, k) + gam * (product - ern);
      qv(i, j, k) = max(qv(i, j, k) - product + ern, real(0.0));
      qc(i, j, k) = qc(i, j, k) + product;
      qr(i, j, k) = qr(i, j, k) - ern;
    }
  }
};

}  // namespace cppmpas::phys::kessler
