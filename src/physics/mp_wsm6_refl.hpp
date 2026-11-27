// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// refl10cm_wsm6 of MMM-physics a4baf7f mp_wsm6.F90 (MW6:2275-2444), the 10 cm radar reflectivity of one WSM6
// column, which MPAS's compute_radar_reflectivity calls (mpas_atmphys_driver_microphysics.F:713). Rayleigh
// scattering at 10 cm; melting snow / graupel below the melting level as water-coated ice spheres
// (rayleigh_soak_wetgraupel, physics/mp_radar.hpp).
// Taken from CPPWRF 7fac5bf src/phys/mp_wsm6/wsm6.hpp:2525-2649 (Wsm6ReflKernel, the same mp_wsm6.F90): the
// statements are CPPWRF's; the column comes in as arrays 1..kte (kts = 1) and the result is dBZ(k), not yet
// limited to -35 (MPAS's driver does that). REAL is double (-fdefault-real-8): CPPWRF's Real(...) casts vanish.
// Lanes as in the reference executable (__mp_wsm6_MOD_refl10cm_wsm6 0x498620): the first loop's pow / exp, the
// ze loop's pow and the melting loops' pow / exp are scalar (pow@plt, exp@plt); the final dBZ loop
// (MW6:2439-2441, k = kte down to kts) is vectorised: _ZGVbN2v_log10 on the pairs (kte, kte-1), (kte-2, kte-3),
// ..., an odd last k = kts (and a 1-trip loop) with log10@plt, the sum (ze_rain + ze_snow) + ze_graupel, times
// 1.d18, log10, times 10. No FMA, no ymm register, no single-precision conversion (sngl is the default REAL
// kind, real*8, in this build).
#include "core/error.hpp"
#include "core/libm.hpp"
#include "physics/mp_radar.hpp"

namespace cppmpas::phys::wsm6 {

// mp_wsm6's module values refl10cm_wsm6 reads (MW6:16-43 parameters; pi, n0g from mp_wsm6_init) and the radar
// tables mp_wsm6_init's radar_init call leaves (MW6:174-184)
struct Wsm6Radar {
  double n0r = 8.e6, n0s = 2.e6, n0smax = 1.e11, alpha = .12, n0g = 4.e6, pi = 0.;
  radar::RadarTables t{};
};

// mp_wsm6_init MW6:102 (pi) and 174-184: xam_x = PI * den_x / 6., xbm = 3., xmu = 0., radar_init
inline Wsm6Radar wsm6_radar_init(double pi, double denr, double dens, double deng, double n0g) {
  Wsm6Radar r;
  r.pi = pi;
  r.n0g = n0g;
  r.t = radar::radar_init(pi * denr / 6., 3., 0., pi * dens / 6., 3., 0., pi * deng / 6., 3., 0.);
  return r;
}

inline constexpr int kReflMaxLevels = 128;

// MW6:2275-2444 on the column qv1d ... p1d (1..kte), dBZ out (1..kte)
CPPMPAS_HD inline void refl10cm_wsm6(const Wsm6Radar& w, const double* qv1d, const double* qr1d, const double* qs1d,
                          const double* qg1d, const double* t1d, const double* p1d, double* dBZ, int kte) {
  using radar::detail::c_max;
  using radar::detail::c_min;
  if (kte < 1 || kte > kReflMaxLevels) ::cppmpas::kernel_fail("refl10cm_wsm6: nVertLevels outside 1..128");
  const radar::RadarTables& t = w.t;
  const int nz = kte;
  constexpr double R = 287.;
  bool l_qr[kReflMaxLevels], l_qs[kReflMaxLevels], l_qg[kReflMaxLevels];
  double temp[kReflMaxLevels], rr[kReflMaxLevels], rs[kReflMaxLevels], rg[kReflMaxLevels];
  double ze_rain[kReflMaxLevels], ze_snow[kReflMaxLevels], ze_graupel[kReflMaxLevels];
  double ilamr[kReflMaxLevels], ilams[kReflMaxLevels], ilamg[kReflMaxLevels];
  double n0_r[kReflMaxLevels], n0_s[kReflMaxLevels], n0_g[kReflMaxLevels];
  for (int k = 0; k < nz; ++k) dBZ[k] = -35.0;  // MW6:2304-2306
  for (int k = 0; k < nz; ++k) {               // MW6:2311-2350
    temp[k] = t1d[k];
    const double temp_c = c_min(-0.001, temp[k] - 273.15);
    const double qvk = c_max(1.e-10, qv1d[k]);
    const double pres = p1d[k];
    const double rho = 0.622 * pres / (R * temp[k] * (qvk + 0.622));
    n0_r[k] = n0_s[k] = n0_g[k] = ilamr[k] = ilams[k] = ilamg[k] = 0.0;  // read only where l_q* is set
    if (qr1d[k] > 1.e-9) {
      rr[k] = qr1d[k] * rho;
      n0_r[k] = w.n0r;
      const double lamr = libm::pow(t.xam_r * t.xcrg[2] * n0_r[k] / rr[k], 1. / t.xcre[0]);
      ilamr[k] = 1. / lamr;
      l_qr[k] = true;
    } else {
      rr[k] = 1.e-12;
      l_qr[k] = false;
    }
    if (qs1d[k] > 1.e-9) {
      rs[k] = qs1d[k] * rho;
      n0_s[k] = c_min(w.n0smax, w.n0s * libm::exp(-w.alpha * temp_c));
      const double lams = libm::pow(t.xam_s * t.xcsg[2] * n0_s[k] / rs[k], 1. / t.xcse[0]);
      ilams[k] = 1. / lams;
      l_qs[k] = true;
    } else {
      rs[k] = 1.e-12;
      l_qs[k] = false;
    }
    if (qg1d[k] > 1.e-9) {
      rg[k] = qg1d[k] * rho;
      n0_g[k] = w.n0g;
      const double lamg = libm::pow(t.xam_g * t.xcgg[2] * n0_g[k] / rg[k], 1. / t.xcge[0]);
      ilamg[k] = 1. / lamg;
      l_qg[k] = true;
    } else {
      rg[k] = 1.e-12;
      l_qg[k] = false;
    }
  }
  // MW6:2355-2365: the level of the start of melting (k_0 the level above)
  bool melti = false;
  int k_0 = 0;
  for (int k = nz - 2; k >= 0; --k)
    if (temp[k] > 273.15 && l_qr[k] && (l_qs[k + 1] || l_qg[k + 1])) {
      k_0 = k + 1 > k_0 ? k + 1 : k_0;
      melti = true;
      break;
    }
  // MW6:2373-2384: Rayleigh: rain at all temperatures, dry snow and graupel
  const double cs = 0.176 / 0.93 * (6.0 / w.pi) * (6.0 / w.pi) * (t.xam_s / 900.0) * (t.xam_s / 900.0);
  const double cg = 0.176 / 0.93 * (6.0 / w.pi) * (6.0 / w.pi) * (t.xam_g / 900.0) * (t.xam_g / 900.0);
  for (int k = 0; k < nz; ++k) {
    ze_rain[k] = 1.e-22;
    ze_snow[k] = 1.e-22;
    ze_graupel[k] = 1.e-22;
    if (l_qr[k]) ze_rain[k] = n0_r[k] * t.xcrg[3] * libm::pow(ilamr[k], t.xcre[3]);
    if (l_qs[k]) ze_snow[k] = cs * n0_s[k] * t.xcsg[3] * libm::pow(ilams[k], t.xcse[3]);
    if (l_qg[k]) ze_graupel[k] = cg * n0_g[k] * t.xcgg[3] * libm::pow(ilamg[k], t.xcge[3]);
  }
  // MW6:2395-2437: melting snow and graupel, the meltwater fraction from the amount at the melting level
  if (melti && k_0 >= 1) {
    for (int k = k_0 - 1; k >= 0; --k) {
      if (l_qs[k] && l_qs[k_0]) {
        const double fmelt_s = c_max(0.005, c_min(1.0 - rs[k] / rs[k_0], 0.99));
        double eta = 0.0;
        const double lams = 1. / ilams[k];
        for (int n = 0; n < radar::nrbins; ++n) {
          const double x = t.xam_s * libm::pow(t.xxDs[n], t.xbm_s);
          const double cback = radar::rayleigh_soak_wetgraupel(x, t.xocms, t.xobms, fmelt_s, t.m_w_0, t.m_i_0, t.pi5,
                                                               t.lamda4);
          const double f_d = n0_s[k] * libm::pow(t.xxDs[n], t.xmu_s) * libm::exp(-lams * t.xxDs[n]);
          eta = eta + f_d * cback * t.simpson[n] * t.xdts[n];
        }
        ze_snow[k] = t.lamda4 / (t.pi5 * t.k_w) * eta;
      }
      if (l_qg[k] && l_qg[k_0]) {
        const double fmelt_g = c_max(0.005, c_min(1.0 - rg[k] / rg[k_0], 0.99));
        double eta = 0.0;
        const double lamg = 1. / ilamg[k];
        for (int n = 0; n < radar::nrbins; ++n) {
          const double x = t.xam_g * libm::pow(t.xxDg[n], t.xbm_g);
          const double cback = radar::rayleigh_soak_wetgraupel(x, t.xocmg, t.xobmg, fmelt_g, t.m_w_0, t.m_i_0, t.pi5,
                                                               t.lamda4);
          const double f_d = n0_g[k] * libm::pow(t.xxDg[n], t.xmu_g) * libm::exp(-lamg * t.xxDg[n]);
          eta = eta + f_d * cback * t.simpson[n] * t.xdtg[n];
        }
        ze_graupel[k] = t.lamda4 / (t.pi5 * t.k_w) * eta;
      }
    }
  }
  // MW6:2439-2441, vectorised from k = kte: lane j = kte - k of nz
  for (int k = nz - 1; k >= 0; --k) {
    const int j = nz - 1 - k;
    const bool lane = nz > 1 && !((nz & 1) && j == nz - 1);
    const double ze = (ze_rain[k] + ze_snow[k] + ze_graupel[k]) * 1.e18;
    dBZ[k] = 10. * (lane ? libm::mvec_log10(ze) : libm::log10(ze));
  }
}

}  // namespace cppmpas::phys::wsm6
