// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// The MYNN PBL's module constants (MMM-physics a4baf7f bl_mynn_subroutines.F90:2-72, module bl_mynn_common, BS), the
// parameters of bl_mynn_subroutines (BS:360-404), bl_mynn_init (bl_mynn.F90:24-122, BM) and the bl_mynn_* options of
// MPAS's mynn_bl_driver (physics_wrf/module_bl_mynn.F, MB). The scheme is physics/bl_mynn.hpp.
#include "core/precision.hpp"

namespace cppmpas::phys::mynn {

// bl_mynn_common (BS:2-72): the host constants bl_mynn_init sets (BM:24-122) and the derived ones
struct Const {
  real cp, cpv, cice, cliq, ep_1, ep_2, grav, karman, p1000mb, r_d, r_v, rvovrd, svp1, svp2, svp3, svpt0, xlf, xls,
      xlv;
  real ep_3, gtr, p608, rcp, t0c, tv0, xlscp, xlvcp;
};
// parameters of bl_mynn_common (BS:53-70)
inline constexpr real tice = 240.0, tkmin = 253.0, tref = 300.0, onethird = 1. / 3., twothirds = 2. / 3.;
inline constexpr real b1 = 24.0, b2 = 15.0, cphh_st = 5.0, cphm_st = 5.0, cphh_unst = 16.0, cphm_unst = 16.0;
// parameters of bl_mynn_subroutines (BS:360-404); gfortran folds them one correctly rounded operation per node
inline constexpr real pr = 0.74, g1 = 0.235, c2 = 0.729, c3 = 0.340, c4 = 0.0, c5 = 0.2;
inline constexpr real a1 = b1 * (1.0 - 3.0 * g1) / 6.0;
inline constexpr real c1 = g1 - 1.0 / (3.0 * a1 * 2.88449914061481660);
inline constexpr real a2 = a1 * (g1 - c1) / (g1 * pr);
inline constexpr real g2 = b2 / b1 * (1.0 - c3) + 2.0 * a1 / b1 * (3.0 - 2.0 * c2);
inline constexpr real cc2 = 1.0 - c2, cc3 = 1.0 - c3;
inline constexpr real e1c = 3.0 * a2 * b2 * cc3, e2c = 9.0 * a1 * a2 * cc2, e3c = 9.0 * a2 * a2 * cc2 * (1.0 - c5);
inline constexpr real e4c = 12.0 * a1 * a2 * cc2, e5c = 6.0 * a1 * a1;
inline constexpr real qmin = 0.0, zmax = 1.0, Sqfac = 3.0;
inline constexpr real gpw = 5. / 3., qcgmin = 1.e-8, qkemin = 1.e-12;
inline constexpr real tliq = 269.;
inline constexpr real rr2 = 0.7071068, rrp = 0.3989423;
inline constexpr real CKmod = 1.;

// BM:24-122 bl_mynn_init (MPAS passes mpas_atmphys_constants: mpas_atmphys_driver_pbl.F:746-747)
inline Const bl_mynn_init(real con_cp, real con_cpv, real con_cice, real con_cliq, real con_ep1, real con_ep2,
                          real con_grav, real con_karman, real con_p0, real con_rd, real con_rv, real con_svp1,
                          real con_svp2, real con_svp3, real con_svpt0, real con_xlf, real con_xls, real con_xlv) {
  Const c{};
  c.cp = con_cp;  // BM:93-111
  c.cpv = con_cpv;
  c.cliq = con_cliq;
  c.cice = con_cice;
  c.ep_1 = con_ep1;
  c.ep_2 = con_ep2;
  c.grav = con_grav;
  c.karman = con_karman;
  c.p1000mb = con_p0;
  c.r_d = con_rd;
  c.r_v = con_rv;
  c.rvovrd = c.r_v / c.r_d;
  c.svp1 = con_svp1;
  c.svp2 = con_svp2;
  c.svp3 = con_svp3;
  c.svpt0 = con_svpt0;
  c.xlf = con_xlf;
  c.xls = con_xls;
  c.xlv = con_xlv;
  c.ep_3 = 1. - c.ep_2;  // BM:114-121
  c.gtr = c.grav / tref;
  c.p608 = c.ep_1;
  c.rcp = c.r_d / c.cp;
  c.t0c = c.svpt0;
  c.tv0 = c.p608 * tref;
  c.xlscp = (c.xlv + c.xlf) / c.cp;
  c.xlvcp = c.xlv / c.cp;
  return c;
}

// The bl_mynn_* options of mynn_bl_driver (physics_wrf/module_bl_mynn.F, MB:286-315: the integer namelist options
// become logicals) and the flags MPAS passes (mpas_atmphys_vars.F: spp_pbl = icloud_bl = 0)
struct Options {
  bool tkeadvect = false;
  int tkebudget = 0, cloudpdf = 2, mixlength = 2, stfunc = 1, topdown = 0, scaleaware = 1, dheat_opt = 1, edmf = 1,
      edmf_dd = 0, edmf_mom = 0, edmf_tke = 0, output = 0, mixscalars = 1, cloudmix = 1, mixqt = 0;
  real closure = 2.5;
  int icloud_bl = 0, spp_pbl = 0;
  bool do_DAcycling = false;
  // the option values without a verified port, or nullptr
  const char* refused() const {
    if (tkeadvect) return "config_mynn_tkeadvect = .true.";
    if (tkebudget != 0) return "config_mynn_tkebudget /= 0";
    if (cloudpdf != 2) return "config_mynn_cloudpdf /= 2";
    if (mixlength != 2) return "config_mynn_mixlength /= 2";
    if (stfunc != 1) return "config_mynn_stfunc /= 1";
    if (topdown != 0) return "config_mynn_topdown /= 0";
    if (scaleaware != 1 && scaleaware != 0) return "config_mynn_scaleaware not 0 or 1";
    if (dheat_opt != 1 && dheat_opt != 0) return "config_mynn_dheat_opt not 0 or 1";
    if (edmf != 1) return "config_mynn_edmf /= 1";
    if (edmf_dd != 0) return "config_mynn_edmf_dd /= 0";
    if (edmf_mom != 0) return "config_mynn_edmf_mom /= 0";
    if (edmf_tke != 0) return "config_mynn_edmf_tke /= 0";
    if (output != 0) return "config_mynn_edmf_output /= 0";
    if (mixscalars != 1) return "config_mynn_mixscalars /= 1";
    if (cloudmix != 1) return "config_mynn_mixclouds /= 1";
    if (mixqt != 0) return "config_mynn_mixqt /= 0";
    if (!(closure <= 2.5)) return "config_mynn_closure > 2.5";
    if (icloud_bl != 0) return "icloud_bl /= 0";
    if (spp_pbl != 0) return "spp_pbl /= 0";
    return nullptr;
  }
};

}  // namespace cppmpas::phys::mynn
