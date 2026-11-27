// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// MPAS-Atmosphere's Grell-Freitas convection (physics/physics_wrf/module_cu_gf.mpas.F, GF, built with -Dmpas),
// statement by statement, for one column: the column part of cu_grell_freitas (GF:54-453), CUP_gf (GF:458-1986),
// CUP_gf_sh (GF:4151-5085) and the routines they call: cup_dd_edt (GF:1989-2124), cup_dd_moisture_new
// (GF:2127-2263), cup_env (GF:2265-2397), cup_env_clev (GF:2400-2502), cup_forcing_ens_3d (GF:2505-2984), cup_kbcon
// (GF:2987-3103), cup_MAXIMI (GF:3176-3229), cup_minimi (GF:3232-3280), cup_up_aa0 (GF:3283-3350), neg_check
// (GF:3475-3566), cup_output_ens_3d (GF:3569-3786), cup_up_moisture (GF:3788-4058) and satvap (GF:4126-4149).
// (g3init, gdinit and cup_ktop are not called by MPAS.)
//
// MPAS's own scheme: WRF 4.8.0's Grell-Freitas is the GFL rewrite (phys/module_cu_gfl_*.F, CPPWRF src/phys/cu_gfl),
// a different code; CPPWRF's port served only as a guide to layout (SHARED_WITH_CPPWRF.md).
//
// Every column is independent: each loop "do i" of the Fortran touches column i only (no cross-column sum; the
// random-number branches are off: irandom = 0), so the columns run one at a time, each with its arrays (1-based,
// levels 1..ktf) in the order of the Fortran's k loops. ktf = kte (kde = kte + 1), itf = ite.
// The module's parameters as MPAS sets them (GF:68-90, 733-737, 4356-4373): autoconv = 1, aeroevap = 1,
// training = 0, use_excess = use_excess_sh = 0, maxiens = maxens = maxens2 = 1, maxens3 = ensdim = 16, iens = 1,
// ens4 = 9, ccnclean = 250, beta = 0.02, tcrit = 258; cup_env is always called with itest = -1. Only the branches
// these values reach are ported (the Berry autoconversion, aerosol evaporation, use_excess >= 1, maxens > 1 and
// itest >= 0 code is unreachable in MPAS). ens4 = 9 copies of omeg / mconv are filled with the same values in
// cu_grell_freitas (GF:315-331), so one column value stands for the nine; the loops over the nine (a_ave,
// xff_ens3(5), (9), (14), the omeg minimum / maximum) keep their order and repetitions. ierrc (the error text) and the
// totmas / neg della diagnostics written to unit 0 (GF:1686-1692, 4768-4775, 4824-4825) are not reproduced: they set
// no value.
//
// Reference object code (module_cu_gf.mpas.F built standalone with the reference flags, -Dmpas): no FMA, no ymm
// register, no libmvec call. libm: satvap's log (2) and pow(10., x) (4), scalar (log@plt / pow@plt; satvap is not
// inlined into cup_env); cup_env's log calls (itest = 0, 1) and cup_dd_edt's pow (aeroevap > 1) are unreachable.
// log(10.), log(6.1071) / log(10.) and log(1013.246) / log(10.) are folded constants (kLn10, kLog6_1071,
// kLog1013_246). x**2, x**3 are multiplies; CUP_gf_sh's float(k)**power_entr (power_entr = 2.) is folded to a
// multiply (no pow call in cup_gf_sh). gfortran's MAX(a, b) / MIN(a, b) are maxsd / minsd: maxr / minr below
// (a > b ? a : b with a constant second, as gfortran orders them); the order matters only for a NaN or -0.0, which
// the sums that follow absorb (each starts from +0.).
//
// Undefined in MPAS (not reproduced, reported instead):
//  * CUP_gf_sh reads qco(i, ktop + 1), never set, at the cloud top (GF:4784) where it is multiplied by detup = 0;
//    the C++ column holds 0 there (the product is +-0 either way, added to the detupk term).
//  * CUP_gf_sh's loop k = kbcon..kbcon + tun_lim (GF:4538) runs past ktf when kbcon > ktf - 8 (writes past the end of
//    cd): refused (kernel_fail).
#include "core/backend.hpp"
#include "core/libm.hpp"
#include "core/precision.hpp"

namespace cppmpas::phys::cu_gf {

inline constexpr int kMaxLevels = 100;
inline constexpr int KM = kMaxLevels + 2;  // 1-based levels 1..ktf, ktf <= kMaxLevels

// gfortran MAX / MIN as maxsd / minsd (header)
CPPMPAS_HD inline real maxr(real a, real b) { return a > b ? a : b; }
CPPMPAS_HD inline real minr(real a, real b) { return a < b ? a : b; }
CPPMPAS_HD inline int maxi(int a, int b) { return a > b ? a : b; }
CPPMPAS_HD inline int mini(int a, int b) { return a < b ? a : b; }

// folded constants of satvap (the reference object's .rodata)
inline constexpr real kLn10 = 0x1.26bb1bbb55516p+1;        // log(10.)
inline constexpr real kLog6_1071 = 0x1.9258f81f79246p-1;   // log(6.1071) / log(10.)
inline constexpr real kLog1013_246 = 0x1.80bb4409e541ap+1;  // log(1013.246) / log(10.)

inline constexpr int maxens3 = 16, ens4 = 9;
inline constexpr real ccnclean = 250., beta = 0.02, tcrit = 258.;  // GF:87-90

// GF:4126-4149
CPPMPAS_HD inline real satvap(real temp2) {
  const real temp = temp2 - 273.155;
  if (temp < -20.) {  // ice saturation
    const real toot = 273.16 / temp2;
    const real toto = 1 / toot;
    const real eilog = -9.09718 * (toot - 1) - 3.56654 * (libm::log(toot) / kLn10) + .876793 * (1 - toto) + kLog6_1071;
    return libm::pow(10., eilog);
  }
  const real tsot = 373.16 / temp2;
  const real ewlog = -7.90298 * (tsot - 1) + 5.02808 * (libm::log(tsot) / kLn10);
  const real ewlog2 = ewlog - 1.3816e-07 * (libm::pow(10., 11.344 * (1 - (1 / tsot))) - 1);
  const real ewlog3 = ewlog2 + .0081328 * (libm::pow(10., -3.49149 * (tsot - 1)) - 1);
  const real ewlog4 = ewlog3 + kLog1013_246;
  return libm::pow(10., ewlog4);
}

// cup_env with itest = -1 (GF:2265-2397): z is not changed; tv and the AE / BE / HT constants feed nothing here
CPPMPAS_HD inline void cup_env(const real* z, real* qes, real* he, real* hes, const real* t, const real* q,
                               const real* p, int ierr, int ktf) {
  for (int k = 1; k <= ktf; ++k) {  // GF:2331-2349
    if (ierr == 0) {
      const real e = satvap(t[k]);
      qes[k] = 0.622 * e / maxr(p[k] - e, 1.e-8);
      if (qes[k] <= 1.E-08) qes[k] = 1.E-08;
      if (qes[k] < q[k]) qes[k] = q[k];
    }
  }
  for (int k = 1; k <= ktf; ++k) {  // GF:2387-2395
    if (ierr == 0) {
      he[k] = 9.81 * z[k] + 1004. * t[k] + 2.5E06 * q[k];  // itest <= 0
      hes[k] = 9.81 * z[k] + 1004. * t[k] + 2.5E06 * qes[k];
      if (he[k] >= hes[k]) he[k] = hes[k];
    }
  }
}

// GF:2400-2502
CPPMPAS_HD inline void cup_env_clev(const real* t, const real* qes, const real* q, const real* he, const real* hes,
                                    const real* z, const real* p, real* qes_cup, real* q_cup, real* he_cup,
                                    real* hes_cup, real* z_cup, real* p_cup, real* gamma_cup, real* t_cup, real psur,
                                    int ierr, real z1, real xl, real rv, real cp, int ktf) {
  for (int k = 1; k <= ktf; ++k) {  // GF:2456-2467
    qes_cup[k] = 0.;
    q_cup[k] = 0.;
    hes_cup[k] = 0.;
    he_cup[k] = 0.;
    z_cup[k] = 0.;
    p_cup[k] = 0.;
    t_cup[k] = 0.;
    gamma_cup[k] = 0.;
  }
  for (int k = 2; k <= ktf; ++k) {  // GF:2468-2483
    if (ierr == 0) {
      qes_cup[k] = .5 * (qes[k - 1] + qes[k]);
      q_cup[k] = .5 * (q[k - 1] + q[k]);
      hes_cup[k] = .5 * (hes[k - 1] + hes[k]);
      he_cup[k] = .5 * (he[k - 1] + he[k]);
      if (he_cup[k] > hes_cup[k]) he_cup[k] = hes_cup[k];
      z_cup[k] = .5 * (z[k - 1] + z[k]);
      p_cup[k] = .5 * (p[k - 1] + p[k]);
      t_cup[k] = .5 * (t[k - 1] + t[k]);
      gamma_cup[k] = (xl / cp) * (xl / (rv * t_cup[k] * t_cup[k])) * qes_cup[k];
    }
  }
  if (ierr == 0) {  // GF:2484-2500 (z_cup(1) and p_cup(1) are set twice; the first values are overwritten)
    qes_cup[1] = qes[1];
    q_cup[1] = q[1];
    hes_cup[1] = 9.81 * z1 + 1004. * t[1] + 2.5e6 * qes[1];
    he_cup[1] = 9.81 * z1 + 1004. * t[1] + 2.5e6 * q[1];
    z_cup[1] = z1;
    p_cup[1] = psur;
    t_cup[1] = t[1];
    gamma_cup[1] = xl / cp * (xl / (rv * t_cup[1] * t_cup[1])) * qes_cup[1];
  }
}

// GF:3176-3229
CPPMPAS_HD inline int cup_maximi(const real* array, int ks, int ke, int ierr) {
  int maxx = ks;
  if (ierr == 0) {
    real x = array[ks];
    for (int k = ks; k <= ke; ++k) {
      const real xar = array[k];
      if (xar >= x) {
        x = xar;
        maxx = k;
      }
    }
  }
  return maxx;
}

// GF:3232-3280. kte: the upper bound of MPAS's array(its:ite, kts:kte); a search past it is undefined in MPAS
CPPMPAS_HD inline int cup_minimi(const real* array, int ks, int kend, int ierr, int kte) {
  int kt = ks;
  if (ierr == 0) {
    real x = array[ks];
    const int kstop = maxi(ks + 1, kend);
    if (kstop > kte) kernel_fail("cup_minimi: MPAS's search GF:3269-3272 would read past kte (undefined)");
    for (int k = ks + 1; k <= kstop; ++k) {
      if (array[k] < x) {
        x = array[k];
        kt = k;
      }
    }
  }
  return kt;
}

// GF:2987-3103 with use_excess = 0 (zqexec = ztexec = 0)
CPPMPAS_HD inline void cup_kbcon(real cap_inc, int iloop, int& k22, int& kbcon, const real* he_cup,
                                 const real* hes_cup, real& hkb, int& ierr, int kbmax, const real* p_cup, real cap_max,
                                 real xl, real cp, real ztexec, real zqexec, int kte) {
  // kte: the upper bound of MPAS's hes_cup / he_cup(its:ite, kts:kte); a level past it is undefined in MPAS
  constexpr int use_excess = 0;
  kbcon = 1;
  if (ierr != 0) return;
  kbcon = k22 + 1;
  if (iloop == 5) kbcon = k22;
  for (;;) {                      // label 32
    if (kbcon > kte || k22 > kte)
      kernel_fail("cup_kbcon: MPAS's search GF:3045-3100 would read past kte (undefined; kbmax + 2 > kte)");
    const real hetest = hkb;      // GF:3055-3061
    if (hetest < hes_cup[kbcon]) {  // go to 31 (GF:3045-3053)
      kbcon = kbcon + 1;
      if (kbcon > kbmax + 2) {
        if (iloop != 4) ierr = 3;
        return;
      }
      continue;
    }
    if (kbcon - k22 == 1) return;  // GF:3069-3070
    if (iloop == 5 && (kbcon - k22) == 0) return;
    real pbcdif = -p_cup[kbcon] + p_cup[k22];  // GF:3071-3077
    real plus = maxr(cap_max - real(iloop - 1) * cap_inc, 25.);
    if (iloop == 4) plus = cap_max;
    if (iloop == 5) plus = 25.;
    if (iloop == 5 && cap_max > 25) pbcdif = -p_cup[kbcon] + cap_max;
    if (pbcdif > plus) {  // GF:3078-3100
      k22 = k22 + 1;
      kbcon = k22 + 1;
      hkb = he_cup[k22] + real(use_excess) * (xl * zqexec + cp * ztexec);
      if (iloop == 5) kbcon = k22;
      if (kbcon > kbmax + 2) {
        if (iloop != 4) ierr = 3;
        return;
      }
      continue;  // go to 32
    }
    return;
  }
}

// GF:3283-3350
CPPMPAS_HD inline real cup_up_aa0(const real* z, const real* zu, const real* dby, const real* gamma_cup,
                                  const real* t_cup, int kbcon, int ktop, int ierr, int ktf) {
  real aa0 = 0.;
  for (int k = 2; k <= ktf; ++k) {
    if (ierr != 0) continue;
    if (k <= kbcon) continue;
    if (k > ktop) continue;
    const real dz = z[k] - z[k - 1];
    const real da = zu[k] * dz * (9.81 / (1004. * ((t_cup[k])))) * dby[k - 1] / (1. + gamma_cup[k]);
    if (k == ktop && da <= 0.) continue;
    aa0 = aa0 + da;
    if (aa0 < 0.) aa0 = 0.;
  }
  return aa0;
}

// GF:2127-2263 (iloop = 1)
CPPMPAS_HD inline void cup_dd_moisture_new(const real* zd, const real* hcd, const real* hes_cup, real* qcd,
                                           const real* qes_cup, real* pwd, const real* q_cup, const real* z_cup,
                                           const real* dd_massentr, const real* dd_massdetr, int jmin, int& ierr,
                                           const real* gamma_cup, real& pwev, real& bu, real* qrcd, const real* q,
                                           int iloop, real xl, int ktf) {
  bu = 0.;
  pwev = 0.;
  for (int k = 1; k <= ktf; ++k) {
    qcd[k] = 0.;
    qrcd[k] = 0.;
    pwd[k] = 0.;
  }
  if (ierr != 0) return;
  int k = jmin;  // GF:2204-2218
  real dz = z_cup[k + 1] - z_cup[k];
  qcd[k] = q_cup[k];
  real dh = hcd[k] - hes_cup[k];
  if (dh < 0)
    qrcd[k] = (qes_cup[k] + (1. / xl) * (gamma_cup[k] / (1. + gamma_cup[k])) * dh);
  else
    qrcd[k] = qes_cup[k];
  pwd[jmin] = zd[jmin] * minr(qcd[k] - qrcd[k], 0.);
  qcd[k] = qrcd[k];
  pwev = pwev + pwd[jmin];
  bu = dz * dh;
  for (int ki = jmin - 1; ki >= 1; --ki) {  // GF:2219-2247
    dz = z_cup[ki + 1] - z_cup[ki];
    qcd[ki] = (qcd[ki + 1] * zd[ki + 1] - .5 * dd_massdetr[ki] * qcd[ki + 1] + dd_massentr[ki] * q[ki]) /
              (zd[ki + 1] - .5 * dd_massdetr[ki] + dd_massentr[ki]);
    dh = hcd[ki] - hes_cup[ki];
    bu = bu + dz * dh;
    qrcd[ki] = qes_cup[ki] + (1. / xl) * (gamma_cup[ki] / (1. + gamma_cup[ki])) * dh;
    real dqeva = qcd[ki] - qrcd[ki];
    if (dqeva > 0.) {
      dqeva = 0.;
      qrcd[ki] = qcd[ki];
    }
    pwd[ki] = zd[ki] * dqeva;
    qcd[ki] = qrcd[ki];
    pwev = pwev + pwd[ki];
  }
  if (pwev == 0 && iloop == 1) ierr = 7;  // GF:2250-2259
  if (bu >= 0 && iloop == 1) ierr = 7;
}

// GF:3788-4058 for name = 'deep', use_excess = 0, autoconv = 1 (c0 = .002, iall = 0). qch, qrcb, clw_allh, pwh,
// pwavh, psumh and prop_b feed nothing for autoconv = 1 (Berry) and are left out; clw_all is computed, not used by
// CUP_gf.
CPPMPAS_HD inline void cup_up_moisture(int ierr, const real* z_cup, real* qc, real* qrc, real* pw, real& pwav,
                                       int kbcon, int ktop, const real* dby, real* clw_all, const real* q,
                                       const real* gamma_cup, const real* zu, const real* qes_cup, int k22,
                                       const real* qe_cup, real xl, real zqexec, const real* up_massentr,
                                       const real* up_massdetr, real& psum, int ktf) {
  constexpr int use_excess = 0;
  const real c0 = .002;
  pwav = 0.;
  psum = 0.;
  for (int k = 1; k <= ktf; ++k) {  // GF:3889-3901
    pw[k] = 0.;
    qc[k] = 0.;
    if (ierr == 0) qc[k] = qes_cup[k];
    clw_all[k] = 0.;
    qrc[k] = 0.;
  }
  if (ierr == 0) {  // GF:3902-3918
    for (int k = 2; k <= kbcon - 1; ++k) {
      const real dz = z_cup[k] - z_cup[k - 1];
      qc[k] = qe_cup[k22] + real(use_excess) * zqexec;
      if (qc[k] > qes_cup[kbcon - 1]) {
        pw[k] = zu[k] * (qc[k] - qes_cup[kbcon - 1]);
        qc[k] = qes_cup[kbcon - 1];
        pwav = pwav + pw[k];
        psum = psum + pw[k] * dz;
      }
    }
  }
  for (int k = 2; k <= ktf; ++k) {  // GF:3940-4047
    if (ierr != 0) continue;
    if (k < kbcon) continue;
    if (k > ktop) continue;
    const real dz = z_cup[k] - z_cup[k - 1];
    const real qrch = qes_cup[k] + (1. / xl) * (gamma_cup[k] / (1. + gamma_cup[k])) * dby[k];
    qc[k] = (qc[k - 1] * zu[k - 1] - .5 * up_massdetr[k - 1] * qc[k - 1] + up_massentr[k - 1] * q[k - 1]) /
            (zu[k - 1] - .5 * up_massdetr[k - 1] + up_massentr[k - 1]);
    if (qc[k] <= qrch) qc[k] = qrch;
    clw_all[k] = qc[k] - qrch;
    qrc[k] = (qc[k] - qrch);
    // autoconv = 1 (GF:4025-4041)
    qrc[k] = ((qc[k] - qrch) * zu[k] - qrc[k - 1] * (.5 * up_massdetr[k - 1])) /
             (zu[k] + .5 * up_massdetr[k - 1] + c0 * dz * zu[k]);
    pw[k] = c0 * dz * qrc[k] * zu[k];
    if (qrc[k] < 0) {
      qrc[k] = 0.;
      pw[k] = 0.;
    }
    qc[k] = qrc[k] + qrch;
    pwav = pwav + pw[k];  // GF:4045-4046
    psum = psum + clw_all[k] * zu[k] * dz;
  }
}

// GF:1989-2124 with maxens2 = 1, aeroevap = 1
CPPMPAS_HD inline void cup_dd_edt(int ierr, const real* us, const real* vs, const real* z, int ktop, int kbcon,
                                  real& edt, const real* p, real pwav, real pwev, real edtmax, real edtmin, real& edtc,
                                  int ktf) {
  edt = 0.;
  real vws = 0., sdp = 0., vshear = 0.;
  edtc = 0.;
  for (int kk = 1; kk <= ktf - 1; ++kk) {  // GF:2056-2068
    if (ierr != 0) continue;
    if (kk <= mini(ktop, ktf) && kk >= kbcon) {
      vws = vws + (std::fabs((us[kk + 1] - us[kk]) / (z[kk + 1] - z[kk])) +
                   std::fabs((vs[kk + 1] - vs[kk]) / (z[kk + 1] - z[kk]))) *
                      (p[kk] - p[kk + 1]);
      sdp = sdp + p[kk] - p[kk + 1];
    }
    if (kk == ktf - 1) vshear = 1.e3 * vws / sdp;
  }
  if (ierr == 0) {  // GF:2069-2113
    real pef = (1.591 - .639 * vshear + .0953 * (vshear * vshear) - .00496 * (vshear * vshear * vshear));
    if (pef > 0.9) pef = 0.9;
    if (pef < 0.1) pef = 0.1;
    const real zkbc = z[kbcon] * 3.281e-3;
    real prezk = .02;
    if (zkbc > 3.)
      prezk = .96729352 + zkbc * (-.70034167 + zkbc * (.162179896 + zkbc * (-1.2569798E-2 + zkbc * (4.2772E-4 - zkbc * 5.44E-6))));
    if (zkbc > 25) prezk = 2.4;
    real pefb = 1. / (1. + prezk);
    if (pefb > 0.9) pefb = 0.9;
    if (pefb < 0.1) pefb = 0.1;
    edt = 1. - .5 * (pefb + pef);
    const real einc = .2 * edt;
    edtc = edt + real(1 - 2) * einc;  // k = 1
  }
  if (ierr == 0) {  // GF:2114-2122
    edtc = -edtc * pwav / pwev;
    if (edtc > edtmax) edtc = edtmax;
    if (edtc < edtmin) edtc = edtmin;
  }
}

// GF:2505-2984 for name = 'deeps', maxens = maxens2 = 1, iens = iedt = 1 (nall = 0), irandom = 0. xf_ens and pr_ens
// are 1..16; omeg is the column (one of the nine identical ens4 copies), omeg2 = omeg(i, 2, 1).
CPPMPAS_HD inline void cup_forcing_ens_3d(real& closure_n, real xland, real& aa0, real aa1, real xaa0, real mbdt,
                                          real dtime, int& ierr, int ierr2, int ierr3, real* xf_ens, real axx,
                                          real mconv, const real* omeg, int k22, const real* pr_ens, int kbcon,
                                          int icoic) {
  const real fens4 = real(ens4);
  real ens_adj = 1.;
  if (ierr > 995) {  // GF:2630-2633
    aa0 = 0.;
    ierr = 0;
  }
  if (ierr == 0) {
    ens_adj = 1.;  // GF:2635-2637
    if (ierr2 > 0 && ierr3 == 0) ens_adj = 0.;
    if (ierr2 > 0 && ierr3 > 0) ens_adj = 0.;
    real a_ave = 0.;  // GF:2643-2650
    for (int ne = 1; ne <= ens4; ++ne) a_ave = a_ave + axx;
    a_ave = maxr(a_ave / fens4, 0.);
    a_ave = minr(a_ave, aa1);
    a_ave = maxr(a_ave, 0.);
    real xff_ens3[maxens3 + 1];
    for (int ne = 1; ne <= 16; ++ne) xff_ens3[ne] = 0.;
    const real xff0 = (aa1 - aa0) / dtime;  // GF:2654-2668
    xff_ens3[1] = maxr((aa1 - aa0) / dtime, 0.);
    xff_ens3[2] = maxr((a_ave - aa0) / dtime, 0.);
    xff_ens3[3] = maxr((aa1 - aa0) / dtime, 0.);
    xff_ens3[13] = maxr((aa1 - aa0) / dtime, 0.);
    xff_ens3[14] = 0.;  // GF:2677-2686
    for (int ne = 1; ne <= ens4; ++ne) xff_ens3[14] = xff_ens3[14] - omeg[k22] / (fens4 * 9.81);
    if (xff_ens3[14] < 0.) xff_ens3[14] = 0.;
    xff_ens3[5] = 0.;
    for (int ne = 1; ne <= ens4; ++ne) xff_ens3[5] = xff_ens3[5] - omeg[kbcon] / (fens4 * 9.81);
    if (xff_ens3[5] < 0.) xff_ens3[5] = 0.;
    xff_ens3[4] = -omeg[2] / 9.81;  // GF:2690-2709
    for (int k = 2; k <= kbcon - 1; ++k)
      for (int ne = 1; ne <= ens4; ++ne) {
        const real xomg = -omeg[k] / 9.81;
        if (xomg < xff_ens3[4]) xff_ens3[4] = xomg;
      }
    if (xff_ens3[4] < 0.) xff_ens3[4] = 0.;
    xff_ens3[6] = -omeg[2] / 9.81;
    for (int k = 2; k <= kbcon - 1; ++k)
      for (int ne = 1; ne <= ens4; ++ne) {
        const real xomg = -omeg[k] / 9.81;
        if (xomg > xff_ens3[6]) xff_ens3[6] = xomg;
      }
    if (xff_ens3[6] < 0.) xff_ens3[6] = 0.;
    xff_ens3[5] = xff_ens3[6];
    xff_ens3[4] = xff_ens3[6];
    xff_ens3[7] = mconv;  // GF:2714-2728
    xff_ens3[8] = mconv;
    xff_ens3[9] = mconv;
    for (int ne = 2; ne <= ens4; ++ne)
      if (mconv > xff_ens3[7]) xff_ens3[7] = mconv;
    for (int ne = 2; ne <= ens4; ++ne)
      if (mconv < xff_ens3[8]) xff_ens3[8] = mconv;
    for (int ne = 2; ne <= ens4; ++ne) xff_ens3[9] = xff_ens3[9] + mconv;
    xff_ens3[9] = xff_ens3[9] / fens4;
    xff_ens3[15] = mconv;  // GF:2736
    xff_ens3[10] = aa0 / (60. * 20.);  // GF:2741-2750
    xff_ens3[11] = aa0 / (60. * 20.);
    xff_ens3[16] = aa0 / (60. * 20.);
    xff_ens3[12] = aa0 / (60. * 20.);
    if (icoic == 0) {  // GF:2755-2768
      if (xff0 < 0.) {
        xff_ens3[1] = 0.;
        xff_ens3[2] = 0.;
        xff_ens3[3] = 0.;
        xff_ens3[13] = 0.;
        xff_ens3[10] = 0.;
        xff_ens3[11] = 0.;
        xff_ens3[12] = 0.;
      }
      if (xff0 < 0 && xland < 0.1)
        for (int ne = 1; ne <= 16; ++ne) xff_ens3[ne] = 0.;
    }
    real xk = (xaa0 - aa1) / mbdt;  // GF:2773-2780
    if (xk <= 0 && xk > -1.e-2) xk = -1.e-2;
    if (xk > 0 && xk < 1.e-2) xk = 1.e-2;
    // ne = 1 (GF:2784-2949), nall = 0
    if (xland < 0.1) {  // GF:2806-2824 (maxens > 0)
      if (ierr2 > 0 || ierr3 > 0) {
        xff_ens3[1] = ens_adj * xff_ens3[1];
        xff_ens3[2] = ens_adj * xff_ens3[2];
        xff_ens3[3] = ens_adj * xff_ens3[3];
        xff_ens3[13] = ens_adj * xff_ens3[13];
        xff_ens3[10] = ens_adj * xff_ens3[10];
        xff_ens3[11] = ens_adj * xff_ens3[11];
        xff_ens3[12] = ens_adj * xff_ens3[12];
        xff_ens3[16] = ens_adj * xff_ens3[16];
        xff_ens3[7] = ens_adj * xff_ens3[7];
        xff_ens3[8] = ens_adj * xff_ens3[8];
        xff_ens3[9] = ens_adj * xff_ens3[9];
        xff_ens3[15] = ens_adj * xff_ens3[15];
      }
    }
    // iresulte = 1 (GF:2833-2836)
    if (xff0 >= 0.) {  // GF:2842-2848
      if (xff_ens3[1] > 0) xf_ens[1] = maxr(-xff_ens3[1] / xk, 0.);
      if (xff_ens3[2] > 0) xf_ens[2] = maxr(-xff_ens3[2] / xk, 0.);
      if (xff_ens3[3] > 0) xf_ens[3] = maxr(-xff_ens3[3] / xk, 0.);
      if (xff_ens3[13] > 0) xf_ens[13] = maxr(-xff_ens3[13] / xk, 0.);
    }
    xf_ens[4] = maxr(xff_ens3[4], 0.);  // GF:2852-2870
    xf_ens[5] = maxr(xff_ens3[5], 0.);
    xf_ens[6] = maxr(xff_ens3[6], 0.);
    xf_ens[14] = maxr(xff_ens3[14], 0.);
    real a1 = maxr(pr_ens[7], 1.e-3);
    xf_ens[7] = maxr(xff_ens3[7] / a1, 0.);
    a1 = maxr(pr_ens[8], 1.e-3);
    xf_ens[8] = maxr(xff_ens3[8] / a1, 0.);
    a1 = maxr(pr_ens[9], 1.e-3);
    xf_ens[9] = maxr(xff_ens3[9] / a1, 0.);
    a1 = maxr(pr_ens[15], 1.e-3);
    xf_ens[15] = maxr(xff_ens3[15] / a1, 0.);
    if (xk < 0.) {
      xf_ens[10] = maxr(-xff_ens3[10] / xk, 0.);
      xf_ens[11] = maxr(-xff_ens3[11] / xk, 0.);
      xf_ens[12] = maxr(-xff_ens3[12] / xk, 0.);
      xf_ens[16] = maxr(-xff_ens3[16] / xk, 0.);
    }
    if (icoic >= 1) {  // GF:2871-2889, in order
      closure_n = 0.;
      for (int n = 1; n <= 16; ++n) xf_ens[n] = xf_ens[icoic];
    }
  } else if (ierr != 20 && ierr != 0) {  // GF:2976-2980
    for (int n = 1; n <= 16; ++n) xf_ens[n] = 0.;
  }
}

// GF:3569-3786 for name = 'deep', nx = nx2 = 1, iens = 1 (n = 1..16). xf_ens is updated as MPAS does (its scaled
// values are not used afterwards).
CPPMPAS_HD inline void cup_output_ens_3d(real* xf_ens, int& ierr, const real* dellat, const real* dellaq,
                                         const real* dellaqc, const real* subt_ens, const real* subq_ens, real* subt,
                                         real* subq, real* outtem, real* outq, real* outqc, const real* zu,
                                         real* sub_mas, real& pre, const real* pw, real& xmb, int ktop,
                                         const real* pr_ens, real sig, int ktf) {
  constexpr int nx = 1;
  for (int k = 1; k <= ktf; ++k) {  // GF:3675-3684
    outtem[k] = 0.;
    outq[k] = 0.;
    outqc[k] = 0.;
    subt[k] = 0.;
    subq[k] = 0.;
    sub_mas[k] = 0.;
  }
  pre = 0.;
  xmb = 0.;
  real xfac1 = 0.;
  if (ierr == 0)  // GF:3692-3701
    for (int n = 1; n <= 16; ++n)
      if (pr_ens[n] <= 0.) xf_ens[n] = 0.;
  if (ierr == 0) {  // GF:3708-3737
    int k = 0;
    real xmb_ave = 0.;
    for (int n = 1; n <= 16; ++n) {
      k = k + 1;
      xmb_ave = xmb_ave + xf_ens[n];
    }
    xmb_ave = xmb_ave / real(k);
    if (xmb_ave <= 0.) {
      ierr = 13;
      xmb_ave = 0.;
    }
    xmb = sig * xmb_ave;
    if (xmb == 0.) ierr = 19;
    if (xmb > 100.) ierr = 19;
    xfac1 = xmb;
  }
  for (int k = 1; k <= ktf; ++k) {  // GF:3738-3765
    real dtt = 0., dtts = 0., dtq = 0., dtqs = 0., dtqc = 0., dtpw = 0.;
    if (ierr == 0 && k <= ktop) {
      dtt = dtt + dellat[k];  // n = 1..nx
      dtts = dtts + subt_ens[k];
      dtq = dtq + dellaq[k];
      dtqs = dtqs + subq_ens[k];
      dtqc = dtqc + dellaqc[k];
      dtpw = dtpw + pw[k];
      outtem[k] = xmb * dtt / real(nx);
      subt[k] = xmb * dtts / real(nx);
      outq[k] = xmb * dtq / real(nx);
      subq[k] = xmb * dtqs / real(nx);
      outqc[k] = xmb * dtqc / real(nx);
      pre = pre + xmb * dtpw / real(nx);
      sub_mas[k] = zu[k] * xmb;
    }
  }
  if (ierr == 0)  // GF:3767-3773
    for (int k = 1; k <= 16; ++k) xf_ens[k] = sig * xf_ens[k] * xfac1;
}

// GF:3475-3566
CPPMPAS_HD inline void neg_check(real* subt, real* subq, real dt, const real* q, real* outq, real* outt, real* outqc,
                                 real& pret, int ktf) {
  real thresh = 300.01;
  {
    real qmemf = 1.;
    for (int k = 1; k <= ktf; ++k) {
      const real qmem = (subt[k] + outt[k]) * 86400.;
      if (qmem > 2. * thresh) {
        const real qmem2 = 2. * thresh / qmem;
        qmemf = minr(qmemf, qmem2);
      }
      if (qmem < -thresh) {
        const real qmem2 = -thresh / qmem;
        qmemf = minr(qmemf, qmem2);
      }
    }
    for (int k = 1; k <= ktf; ++k) {
      subq[k] = subq[k] * qmemf;
      subt[k] = subt[k] * qmemf;
      outq[k] = outq[k] * qmemf;
      outt[k] = outt[k] * qmemf;
      outqc[k] = outqc[k] * qmemf;
    }
    pret = pret * qmemf;
  }
  thresh = 1.e-10;
  {
    real qmemf = 1.;
    for (int k = 1; k <= ktf - 1; ++k) {
      const real qmem = subq[k] + outq[k];
      if (std::fabs(qmem) > 0.) {
        const real qtest = q[k] + (subq[k] + outq[k]) * dt;
        if (qtest < thresh) {
          const real qmem1 = outq[k] + subq[k];
          const real qmem2 = (thresh - q[k]) / dt;
          qmemf = minr(qmemf, qmem2 / qmem1);
          qmemf = maxr(qmemf, 0.);
        }
      }
    }
    for (int k = 1; k <= ktf; ++k) {
      subq[k] = subq[k] * qmemf;
      subt[k] = subt[k] * qmemf;
      outq[k] = outq[k] * qmemf;
      outt[k] = outt[k] * qmemf;
      outqc[k] = outqc[k] * qmemf;
    }
    pret = pret * qmemf;
  }
}

// The column inputs of CUP_gf / CUP_gf_sh as cu_grell_freitas builds them (GF:252-385)
struct GfColumn {
  int ktf = 0;
  real zo[KM], t2d[KM], q2d[KM], tn[KM], qo[KM], po[KM], p2d[KM], us[KM], vs[KM], rhoi[KM];
  real tshall[KM], qshall[KM], dhdt[KM];
  real omeg[KM];  // omeg(i, k, n), n = 1..ens4 (identical)
  real mconv = 0;  // mconv(i, n), n = 1..ens4 (identical)
  real ter11 = 0, psur = 0, aaeq = 0, ztexec = 0, zqexec = 0, ccn = 0, dx = 0, xland = 0, tscl_kf = 0;
  int kpbl = 0;
};

// CUP_gf's outputs (the cu_grell_freitas arrays it writes for column i)
struct GfDeepOut {
  real outqc[KM], sub_mas[KM], outt[KM], outq[KM], cupclw[KM], subt[KM], subq[KM];
  real pre = 0, xmb = 0;
  int kbcon = 0, ktop = 0;
};

// GF:458-1986 for column c (name 'deep'). On entry o holds cu_grell_freitas' initial values of the intent(inout)
// arrays (GF:235-312): ktop = kbcon = 0, xmb = 0, pre (pret) = 0 and outqc, outt, outq, subt, subq, sub_mas = 0.
CPPMPAS_HD inline void cup_gf(const GfColumn& c, real dtime, real xl, real rv, real cp, real g, int ichoice,
                              GfDeepOut& o) {
  constexpr int use_excess = 0;
  const int ktf = c.ktf;
  const real* zo = c.zo;
  const real* t = c.t2d;
  const real* q = c.q2d;
  const real* tn = c.tn;
  const real* qo = c.qo;
  const real* po = c.po;
  const real* p = c.p2d;
  const real z1 = c.ter11, psur = c.psur;
  const real* omeg = c.omeg;
  int& kbcon = o.kbcon;
  int& ktop = o.ktop;

  real entr_rate_2d[KM] = {}, mentrd_rate_2d[KM] = {}, he[KM] = {}, hes[KM] = {}, qes[KM] = {}, z[KM] = {};
  real heo[KM] = {}, heso[KM] = {}, qeso[KM] = {};
  real xhe[KM] = {}, xhes[KM] = {}, xqes[KM] = {}, xz[KM] = {}, xt[KM] = {}, xq[KM] = {};
  real qes_cup[KM] = {}, q_cup[KM] = {}, he_cup[KM] = {}, hes_cup[KM] = {}, z_cup[KM] = {}, p_cup[KM] = {},
       gamma_cup[KM] = {}, t_cup[KM] = {};
  real qeso_cup[KM] = {}, qo_cup[KM] = {}, heo_cup[KM] = {}, heso_cup[KM] = {}, zo_cup[KM] = {}, po_cup[KM] = {},
       gammao_cup[KM] = {}, tn_cup[KM] = {};
  real xqes_cup[KM] = {}, xq_cup[KM] = {}, xhe_cup[KM] = {}, xhes_cup[KM] = {}, xz_cup[KM] = {}, xt_cup[KM] = {};
  real dby[KM] = {}, hc[KM] = {}, zu[KM] = {}, zd[KM] = {}, clw_all[KM] = {};
  real dbyo[KM] = {}, qco[KM] = {}, qrcdo[KM] = {}, pwdo[KM] = {}, pwo[KM] = {}, hcdo[KM] = {}, qcdo[KM] = {},
       dbydo[KM] = {}, hco[KM] = {}, qrco[KM] = {}, zuo[KM] = {}, zdo[KM] = {};
  real xdby[KM] = {}, xhc[KM] = {}, xzu[KM] = {};
  real cd[KM] = {}, cdd[KM] = {}, dellah[KM] = {}, dellaq[KM] = {}, dellat[KM] = {}, dellaqc[KM] = {},
       dsubt[KM] = {}, dsubh[KM] = {}, dsubq[KM] = {};
  real up_massentr[KM] = {}, up_massdetr[KM] = {}, up_massentro[KM] = {}, up_massdetro[KM] = {},
       dd_massentro[KM] = {}, dd_massdetro[KM] = {};
  real dellat_ens[KM] = {}, dellaqc_ens[KM] = {}, dellaq_ens[KM] = {}, pwo_ens[KM] = {}, subt_ens[KM] = {},
       subq_ens[KM] = {};
  real smth[KM] = {};
  real xf_ens[maxens3 + 1] = {}, pr_ens[maxens3 + 1] = {};

  const int levadj = 5;  // GF:733-737
  const real zustart = .1, zufinal = 1.;
  real closure_n = 16.;  // GF:738-745
  real xland1 = 1.;
  if (c.xland > 1.5) xland1 = 0.;
  const real cap_max_increment = 25.;

  // GF:752-768
  real entr_rate = 7.e-5;
  real radius = .2 / entr_rate;
  real frh = 3.14 * (radius * radius) / c.dx / c.dx;
  if (frh > 0.7) {
    frh = .7;
    radius = std::sqrt(frh * c.dx * c.dx / 3.14);
    entr_rate = .2 / radius;
  }
  const real sig = (1. - frh) * (1. - frh);
  real mentrd_rate = entr_rate;
  const real xlamdd = mentrd_rate;

  for (int k = 1; k <= ktf; ++k) {  // GF:773-784
    z[k] = zo[k];
    xz[k] = zo[k];
    o.cupclw[k] = 0.;
    cd[k] = 1. * entr_rate;
    cdd[k] = xlamdd;
    hcdo[k] = 0.;
    qrcdo[k] = 0.;
    dellaqc[k] = 0.;
  }
  const real edtmax = 1., edtmin = .1;  // GF:828-838
  const real depth_min = 1000.;
  const real cap_maxs = 75.;
  int kbmax = 1;  // GF:839-853
  real aa0 = 0., aa1 = 0., edt = 0.;
  const int kstabm = ktf - 1;
  int ierr = 0, ierr2 = 0, ierr3 = 0;
  const real cap_max = cap_maxs;
  const real zkbmax = 4000., zcutdown = 3000., z_detr = 1250.;  // GF:857-865
  const real mbdt_ens1 = (real(1) - 3.) * dtime * 1.e-3 + dtime * 5.E-03;  // GF:867-869
  if (ierr != 20)  // GF:876-883
    for (int k = 1; k <= 16; ++k) {
      xf_ens[k] = 0.;
      pr_ens[k] = 0.;
    }

  // GF:887-908
  cup_env(z, qes, he, hes, t, q, p, ierr, ktf);
  cup_env(zo, qeso, heo, heso, tn, qo, po, ierr, ktf);
  cup_env_clev(t, qes, q, he, hes, z, p, qes_cup, q_cup, he_cup, hes_cup, z_cup, p_cup, gamma_cup, t_cup, psur, ierr,
               z1, xl, rv, cp, ktf);
  cup_env_clev(tn, qeso, qo, heo, heso, zo, po, qeso_cup, qo_cup, heo_cup, heso_cup, zo_cup, po_cup, gammao_cup,
               tn_cup, psur, ierr, z1, xl, rv, cp, ktf);
  int kdet = 0;  // set below wherever it is read
  if (ierr == 0) {  // GF:909-934
    if (c.aaeq < -0.1) ierr = 20;
    for (int k = 1; k <= ktf; ++k)
      if (zo_cup[k] > zkbmax + z1) {
        kbmax = k;
        break;
      }
    for (int k = 1; k <= ktf; ++k)
      if (zo_cup[k] > z_detr + z1) {
        kdet = k;
        break;
      }
  }

  // GF:941-966
  int k22 = cup_maximi(heo_cup, 3, kbmax, ierr);
  if (ierr == 0) {
    frh = q_cup[k22] / qes_cup[k22];
    if (omeg[k22] < 0. && frh >= 0.99 && sig < 0.091) ierr = 1200;
    if (k22 >= kbmax) ierr = 2;
  }
  // GF:970-982
  real hkb = 0., hkbo = 0.;
  if (ierr == 0) {
    hkb = he_cup[k22];
    hkbo = heo_cup[k22] + real(use_excess) * (xl * c.zqexec + cp * c.ztexec);
  }
  cup_kbcon(cap_max_increment, 1, k22, kbcon, heo_cup, heso_cup, hkbo, ierr, kbmax, po_cup, cap_max, xl, cp,
            c.ztexec, c.zqexec, ktf);  // GF:985-989
  const int kstabi = cup_minimi(heso_cup, kbcon, kstabm, ierr, ktf);  // GF:993-995
  if (ierr == 0)  // GF:997-1004
    for (int k = k22; k <= kbcon; ++k) {
      frh = q_cup[k] / qes_cup[k];
      if (omeg[k] < -1.e-6 && frh >= 0.99 && sig < 0.091) ierr = 1200;
    }

  // GF:1021-1100: smooth normalized mass flux profile
  int kbegzu = 0;
  real zuhe = 0., dh = 0., dz = 0.;
  if (ierr == 0) {
    for (int k = 1; k <= ktf; ++k) {
      frh = minr(qo_cup[k] / qeso_cup[k], 1.);
      entr_rate_2d[k] = entr_rate * (1.3 - frh);
    }
    zuhe = zustart;
    int kstart = 1;
    frh = (zufinal - zustart) / ((real(kbcon * kbcon)) - (real(kstart * kstart)));
    dh = zuhe - frh * (real(kstart * kstart));
    for (int k = kstart; k <= kbcon - 1; ++k) {
      dz = z_cup[k + 1] - z_cup[k];
      if (p_cup[k] > p_cup[kstabi]) cd[k] = 1.e-6;
      entr_rate_2d[k] = ((frh * (real((k + 1) * (k + 1))) + dh) / zuhe - 1. + cd[k] * dz) / dz;
      zuhe = zuhe + entr_rate_2d[k] * dz * zuhe - cd[k] * dz * zuhe;
    }
    kbegzu = kstabi + 4;
    kbegzu = mini(kbegzu, ktf - 1);
    int kfinalzu = kbegzu + 1;
    for (int k = 1; k <= ktf; ++k) cd[k] = entr_rate_2d[kbcon];
    for (int k = kbcon; k <= kbegzu; ++k) {
      cd[k] = entr_rate_2d[kbcon];
      if (p_cup[k] > p_cup[kstabi]) cd[k] = 1.e-6;
      dz = z_cup[k + 1] - z_cup[k];
      zuhe = zuhe + entr_rate_2d[k] * dz * zuhe - cd[k] * dz * zuhe;
    }
    for (int k = kstabi; k <= ktf - 2; ++k)
      if ((hkb - hes_cup[k]) < 0) {
        kfinalzu = k - 3;
        break;  // go to 411
      }
    kfinalzu = maxi(kfinalzu, kbegzu + 1);
    kfinalzu = mini(kfinalzu, ktf - 1);
    frh = -(0.2 - zuhe) / ((real(kfinalzu * kfinalzu)) - (real(kbegzu * kbegzu)));
    dh = zuhe + frh * (real(kbegzu * kbegzu));
    for (int k = kbegzu + 1; k <= kfinalzu; ++k) {
      dz = z_cup[k + 1] - z_cup[k];
      cd[k] = -((-frh * (real((k + 1) * (k + 1))) + dh) / zuhe - 1. - entr_rate_2d[k] * dz) / dz;
      zuhe = zuhe + entr_rate_2d[k] * dz * zuhe - cd[k] * dz * zuhe;
    }
    for (int k = kfinalzu + 1; k <= ktf; ++k) cd[k] = entr_rate_2d[k];
    for (int k = 2; k <= ktf - 2; ++k) {
      const real dzm1 = z_cup[k] - z_cup[k - 1];
      dz = z_cup[k + 1] - z_cup[k];
      const real dzp1 = z_cup[k + 2] - z_cup[k + 1];
      smth[k] = .25 * (dzm1 * cd[k - 1] + 2. * dz * cd[k] + dzp1 * cd[k + 1]);
    }
    for (int k = 2; k <= ktf - 2; ++k) {
      dz = z_cup[k + 1] - z_cup[k];
      cd[k] = smth[k] / dz;
    }
    for (int k = 1; k <= ktf; ++k) smth[k] = 0.;
    for (int k = 2; k <= ktf - 2; ++k) {
      const real dzm1 = z_cup[k] - z_cup[k - 1];
      dz = z_cup[k + 1] - z_cup[k];
      const real dzp1 = z_cup[k + 2] - z_cup[k + 1];
      smth[k] = .25 * (dzm1 * entr_rate_2d[k - 1] + 2. * dz * entr_rate_2d[k] + dzp1 * entr_rate_2d[k + 1]);
    }
    for (int k = 2; k <= ktf - 2; ++k) {
      dz = z_cup[k + 1] - z_cup[k];
      entr_rate_2d[k] = smth[k] / dz;
    }
    zuhe = zustart;
    for (int k = 2; k <= kbegzu; ++k) {
      dz = z_cup[k + 1] - z_cup[k];
      frh = zuhe;
      zuhe = zuhe + entr_rate_2d[k] * dz * zuhe - cd[k] * dz * zuhe;
    }
  }

  // GF:1190-1265: mass entrainment and detrainment
  for (int k = 1; k <= ktf; ++k) {
    hc[k] = 0.;
    dby[k] = 0.;
    hco[k] = 0.;
    dbyo[k] = 0.;
  }
  if (ierr == 0) {
    for (int k = 1; k <= kbcon - 1; ++k) {
      hc[k] = hkb;
      hco[k] = hkbo;
    }
    const int k = kbcon;
    hc[k] = hkb;
    dby[kbcon] = hkb - hes_cup[k];
    hco[k] = hkbo;
    dbyo[kbcon] = hkbo - heso_cup[k];
  }
  if (ierr == 0) {
    zu[1] = zustart;
    zuo[1] = zustart;
    for (int k = 2; k <= ktf - 1; ++k) {
      dz = zo_cup[k] - zo_cup[k - 1];
      up_massentro[k - 1] = entr_rate_2d[k - 1] * dz * zuo[k - 1];
      up_massdetro[k - 1] = cd[k - 1] * dz * zuo[k - 1];
      zuo[k] = zuo[k - 1] + up_massentro[k - 1] - up_massdetro[k - 1];
      if (zuo[k] < 0.05) {
        zuo[k] = .05;
        up_massdetro[k - 1] = zuo[k - 1] - .05 + up_massentro[k - 1];
        cd[k - 1] = up_massdetro[k - 1] / dz / zuo[k - 1];
      }
      zu[k] = zuo[k];
      up_massentr[k - 1] = up_massentro[k - 1];
      up_massdetr[k - 1] = up_massdetro[k - 1];
    }
    for (int k = kbcon + 1; k <= ktf - 1; ++k) {
      hc[k] = (hc[k - 1] * zu[k - 1] - .5 * up_massdetr[k - 1] * hc[k - 1] + up_massentr[k - 1] * he[k - 1]) /
              (zu[k - 1] - .5 * up_massdetr[k - 1] + up_massentr[k - 1]);
      dby[k] = hc[k] - hes_cup[k];
      hco[k] = (hco[k - 1] * zuo[k - 1] - .5 * up_massdetro[k - 1] * hco[k - 1] + up_massentro[k - 1] * heo[k - 1]) /
               (zuo[k - 1] - .5 * up_massdetro[k - 1] + up_massentro[k - 1]);
      dbyo[k] = hco[k] - heso_cup[k];
    }
    for (int k = kbcon + 1; k <= ktf; ++k)
      if (dbyo[k] < 0) {
        ktop = k - 1;
        break;  // go to 41
      }
    if (ktop < kbcon + 2) ierr = 5;
    for (int k = ktop + 1; k <= ktf; ++k) {
      hc[k] = hes_cup[k];
      hco[k] = heso_cup[k];
      dby[k] = 0.;
      dbyo[k] = 0.;
      zu[k] = 0.;
      zuo[k] = 0.;
      cd[k] = 0.;
      entr_rate_2d[k] = 0.;
      up_massentr[k] = 0.;
      up_massdetr[k] = 0.;
      up_massentro[k] = 0.;
      up_massdetro[k] = 0.;
    }
  }

  // GF:1267-1279
  int kzdown = 0;
  if (ierr == 0) {
    real zktop = (zo_cup[ktop] - z1) * .6;
    zktop = minr(zktop + z1, zcutdown + z1);
    for (int k = 1; k <= ktf; ++k)
      if (zo_cup[k] > zktop) {
        kzdown = k;
        break;  // go to 37
      }
  }
  // GF:1283-1324: downdraft originating level
  int jmin = cup_minimi(heso_cup, k22, kzdown, ierr, ktf);
  if (ierr == 0) {
    int jmini = jmin;
    bool keep_going = true;
    while (keep_going) {
      keep_going = false;
      if (jmini - 1 < kdet) kdet = jmini - 1;
      if (jmini >= ktop - 1) jmini = ktop - 2;
      const int ki = jmini;
      hcdo[ki] = heso_cup[ki];
      dz = zo_cup[ki + 1] - zo_cup[ki];
      dh = 0.;
      for (int k = ki - 1; k >= 1; --k) {
        hcdo[k] = heso_cup[jmini];
        dz = zo_cup[k + 1] - zo_cup[k];
        dh = dh + dz * (hcdo[k] - heso_cup[k]);
        if (dh > 0.) {
          jmini = jmini - 1;
          if (jmini > 5) {
            keep_going = true;
          } else {
            ierr = 9;
            break;  // exit
          }
        }
      }
    }
    jmin = jmini;
    if (jmini <= 5) ierr = 4;
  }
  // GF:1329-1336
  if (ierr == 0)
    if (-zo_cup[kbcon] + zo_cup[ktop] < depth_min) ierr = 6;

  // GF:1342-1354
  for (int k = 1; k <= ktf; ++k) {
    zd[k] = 0.;
    zdo[k] = 0.;
    cdd[k] = 0.;
    dd_massentro[k] = 0.;
    dd_massdetro[k] = 0.;
    hcdo[k] = heso_cup[k];
    dbydo[k] = 0.;
  }
  // GF:1358-1429: normalized downdraft mass flux
  real bud = 0.;
  if (ierr == 0) {
    for (int k = 1; k <= ktf; ++k) mentrd_rate_2d[k] = mentrd_rate;  // mentrd_rate_2d(i,:), kts..kte
    for (int k = 1; k <= jmin; ++k) cdd[k] = xlamdd;
    cdd[jmin] = 0.;
    zd[jmin] = 0.2;
    zdo[jmin] = 0.2;
    frh = (zdo[jmin] - 1.) / (-real((jmin - levadj) * (jmin - levadj)) + real(jmin * jmin));
    dh = zdo[jmin] - frh * real(jmin * jmin);
    zuhe = zdo[jmin];
    for (int ki = jmin - 1; ki >= jmin - levadj; --ki) {
      cdd[ki] = 0.;
      dz = z_cup[ki + 1] - z_cup[ki];
      mentrd_rate_2d[ki] = ((frh * real(ki * ki) + dh) / zuhe - 1.) / dz;
      zuhe = zuhe + mentrd_rate_2d[ki] * dz * zuhe;
    }
    int kstart = maxi(kbcon, kdet) - 1;
    kstart = mini(jmin - levadj, kstart);
    kstart = maxi(2, kstart);
    if (kstart < jmin - levadj - 1) {
      for (int ki = jmin - levadj - 1; ki >= kstart; --ki) {
        dz = z_cup[ki + 1] - z_cup[ki];
        mentrd_rate_2d[ki] = mentrd_rate;
        cdd[ki] = xlamdd;
        zuhe = zuhe - cdd[ki] * dz * zuhe + mentrd_rate_2d[ki] * dz * zuhe;
      }
    }
    frh = (zuhe - beta) / (real(kstart * kstart) - 1.);
    dh = beta - frh;
    mentrd_rate_2d[kstart] = 0.;
    for (int ki = kstart + 1; ki >= 1; --ki) {
      mentrd_rate_2d[ki] = 0.;
      dz = z_cup[ki + 1] - z_cup[ki];
      cdd[ki] = maxr((1. - (frh * real(ki * ki) + dh) / zuhe) / dz, 0.);
      zuhe = zuhe - cdd[ki] * dz * zuhe;
    }
    for (int ki = jmin - 1; ki >= 1; --ki) {
      mentrd_rate = mentrd_rate_2d[ki];
      const real dzo = zo_cup[ki + 1] - zo_cup[ki];
      dd_massentro[ki] = mentrd_rate * dzo * zdo[ki + 1];
      dd_massdetro[ki] = cdd[ki] * dzo * zdo[ki + 1];
      zdo[ki] = zdo[ki + 1] + dd_massentro[ki] - dd_massdetro[ki];
    }
    dbydo[jmin] = hcdo[jmin] - heso_cup[jmin];
    bud = dbydo[jmin] * (zo_cup[jmin + 1] - zo_cup[jmin]);
    for (int ki = jmin - 1; ki >= 1; --ki) {
      const real dzo = zo_cup[ki + 1] - zo_cup[ki];
      hcdo[ki] = (hcdo[ki + 1] * zdo[ki + 1] - .5 * dd_massdetro[ki] * hcdo[ki + 1] + dd_massentro[ki] * heo[ki]) /
                 (zdo[ki + 1] - .5 * dd_massdetro[ki] + dd_massentro[ki]);
      dbydo[ki] = hcdo[ki] - heso_cup[ki];
      bud = bud + dbydo[ki] * dzo;
    }
  }
  if (bud > 0) ierr = 7;

  // GF:1510-1528
  real pwevo = 0., bu = 0.;
  cup_dd_moisture_new(zdo, hcdo, heso_cup, qcdo, qeso_cup, pwdo, qo_cup, zo_cup, dd_massentro, dd_massdetro, jmin,
                      ierr, gammao_cup, pwevo, bu, qrcdo, qo, 1, xl, ktf);
  real pwavo = 0., psum = 0.;
  cup_up_moisture(ierr, zo_cup, qco, qrco, pwo, pwavo, kbcon, ktop, dbyo, clw_all, qo, gammao_cup, zuo, qeso_cup,
                  k22, qo_cup, xl, c.zqexec, up_massentr, up_massdetr, psum, ktf);
  for (int k = 1; k <= ktf; ++k) o.cupclw[k] = qrco[k];

  // GF:1532-1554
  aa0 = cup_up_aa0(z, zu, dby, gamma_cup, t_cup, kbcon, ktop, ierr, ktf);
  aa1 = cup_up_aa0(zo, zuo, dbyo, gammao_cup, tn_cup, kbcon, ktop, ierr, ktf);
  if (ierr == 0)
    if (aa1 == 0.) ierr = 17;
  const real axx = aa1;

  // GF:1559-1562
  real edtc = 0.;
  cup_dd_edt(ierr, c.us, c.vs, zo, ktop, kbcon, edt, po, pwavo, pwevo, edtmax, edtmin, edtc, ktf);
  // iedt = 1 (GF:1563-1964)
  real edto = 0.;
  if (ierr == 0) {
    edt = edtc;
    edto = edtc;
  }
  for (int k = 1; k <= ktf; ++k) {  // GF:1576-1600
    subt_ens[k] = 0.;
    subq_ens[k] = 0.;
    dellat_ens[k] = 0.;
    dellaq_ens[k] = 0.;
    dellaqc_ens[k] = 0.;
    pwo_ens[k] = 0.;
  }
  for (int k = 1; k <= ktf; ++k) {
    dellah[k] = 0.;
    dsubt[k] = 0.;
    dsubh[k] = 0.;
    dellaq[k] = 0.;
    dsubq[k] = 0.;
  }
  if (ierr == 0) {  // GF:1641-1727
    real dp = 100. * (po_cup[1] - po_cup[2]);
    dellah[1] = (edto * zdo[2] * hcdo[2] - edto * zdo[2] * heo_cup[2]) * g / dp;
    dellaq[1] = (edto * zdo[2] * qrcdo[2] - edto * zdo[2] * qo_cup[2]) * g / dp;
    dsubt[1] = 0.;
    dsubq[1] = 0.;
    for (int k = 2; k <= ktop; ++k) {
      real entupk = 0., detupk = 0., entdoj = 0.;
      real detdo = edto * dd_massdetro[k];
      real entdo = edto * dd_massentro[k];
      real entup = up_massentro[k];
      real detup = up_massdetro[k];
      real subin = -zdo[k + 1] * edto;
      real subdown = -zdo[k] * edto;
      if (k == jmin) entdoj = edto * zdo[k];
      if (k == ktop) {
        detupk = zuo[ktop];
        subin = 0.;
        subdown = 0.;
        detdo = 0.;
        entdo = 0.;
        entup = 0.;
        detup = 0.;
      }
      // totmas (GF:1680-1692): written to unit 0 when |totmas| > 1e-6, no value set
      dp = 100. * (po_cup[k] - po_cup[k + 1]);
      dellah[k] = (detup * .5 * (hco[k + 1] + hco[k]) + detdo * .5 * (hcdo[k + 1] + hcdo[k]) - entup * heo[k] -
                   entdo * heo[k] + subin * heo_cup[k + 1] - subdown * heo_cup[k] +
                   detupk * (hco[ktop] - heo_cup[ktop]) - entupk * heo_cup[k22] - entdoj * heo_cup[jmin]) *
                  g / dp;
      dellaq[k] = (detup * .5 * (qco[k + 1] + qco[k] - qrco[k + 1] - qrco[k]) + detdo * .5 * (qrcdo[k + 1] + qrcdo[k]) -
                   entup * qo[k] - entdo * qo[k] + subin * qo_cup[k + 1] - subdown * qo_cup[k] +
                   detupk * (qco[ktop] - qrco[ktop] - qo_cup[ktop]) - entupk * qo_cup[k22] - entdoj * qo_cup[jmin]) *
                  g / dp;
      if (k < ktop) {
        dsubt[k] = (zuo[k + 1] * heo_cup[k + 1] - zuo[k] * heo_cup[k]) * g / dp;
        dsubq[k] = (zuo[k + 1] * qo_cup[k + 1] - zuo[k] * qo_cup[k]) * g / dp;
      }
    }
  }
  // GF:1731-1746
  for (int k = 1; k <= ktf - 1; ++k) {
    dellaqc[k] = 0.;
    if (ierr == 0) {
      if (k == ktop - 0) dellaqc[k] = .01 * zuo[ktop] * qrco[ktop] * 9.81 / (po_cup[k] - po_cup[k + 1]);
      if (k < ktop && k > kbcon) {
        dz = zo_cup[k + 1] - zo_cup[k];
        dellaqc[k] = .01 * 9.81 * up_massdetro[k] * .5 * (qrco[k] + qrco[k + 1]) / (po_cup[k] - po_cup[k + 1]);
      }
      dellaqc[k] = maxr(dellaqc[k], 0.);
    }
  }
  // GF:1750-1778
  const real mbdt = mbdt_ens1;
  real xaa0_ens = 0.;
  for (int k = 1; k <= ktf; ++k) {
    dellat[k] = 0.;
    if (ierr == 0) {
      dsubh[k] = dsubt[k];
      xhe[k] = (dsubt[k] + dellah[k]) * mbdt + heo[k];
      xq[k] = (dsubq[k] + dellaq[k] + dellaqc[k]) * mbdt + qo[k];
      dellat[k] = (1. / cp) * (dellah[k] - xl * dellaq[k]);
      dsubt[k] = (1. / cp) * (dsubt[k] - xl * dsubq[k]);
      xt[k] = (dellat[k] + dsubt[k] - dellaqc[k] * xl / cp) * mbdt + tn[k];
      if (xq[k] <= 0.) xq[k] = 1.E-08;
    }
  }
  real xhkb = 0.;
  if (ierr == 0) {
    xhkb = hkbo + (dsubh[k22] + dellah[k22]) * mbdt;
    xhe[ktf] = heo[ktf];
    xq[ktf] = qo[ktf];
    xt[ktf] = tn[ktf];
    if (xq[ktf] <= 0.) xq[ktf] = 1.E-08;
  }
  // GF:1782-1793 (cup_env_clev writes po_cup and gamma_cup again, as MPAS passes them)
  cup_env(xz, xqes, xhe, xhes, xt, xq, po, ierr, ktf);
  cup_env_clev(xt, xqes, xq, xhe, xhes, xz, po, xqes_cup, xq_cup, xhe_cup, xhes_cup, xz_cup, po_cup, gamma_cup,
               xt_cup, psur, ierr, z1, xl, rv, cp, ktf);
  // GF:1805-1848
  for (int k = 1; k <= ktf; ++k) {
    xhc[k] = 0.;
    xdby[k] = 0.;
  }
  if (ierr == 0) {
    for (int k = 1; k <= kbcon - 1; ++k) xhc[k] = xhkb;
    const int k = kbcon;
    xhc[k] = xhkb;
    xdby[kbcon] = xhkb - xhes_cup[k];
  }
  if (ierr == 0) {
    for (int k = 1; k <= ktf; ++k) xzu[k] = zuo[k];  // xzu(i,:) = zuo(i,:), kts..kte
    for (int k = kbcon + 1; k <= ktop; ++k) {
      xhc[k] = (xhc[k - 1] * xzu[k - 1] - .5 * up_massdetro[k - 1] * xhc[k - 1] + up_massentro[k - 1] * xhe[k - 1]) /
               (xzu[k - 1] - .5 * up_massdetro[k - 1] + up_massentro[k - 1]);
      xdby[k] = xhc[k] - xhes_cup[k];
    }
    for (int k = ktop + 1; k <= ktf; ++k) {
      xhc[k] = xhes_cup[k];
      xdby[k] = 0.;
      xzu[k] = 0.;
    }
  }
  // GF:1853-1904
  const real xaa0 = cup_up_aa0(xz, xzu, xdby, gamma_cup, xt_cup, kbcon, ktop, ierr, ktf);
  if (ierr == 0) {  // nens = 1, nall = 0
    xaa0_ens = xaa0;
    for (int k = 1; k <= ktf; ++k)
      if (k <= ktop)
        for (int nens3 = 1; nens3 <= maxens3; ++nens3) pr_ens[nens3] = pr_ens[nens3] + pwo[k];  // all branches
    if (pr_ens[7] < 1.e-6) {
      ierr = 18;
      for (int nens3 = 1; nens3 <= maxens3; ++nens3) pr_ens[nens3] = 0.;
    }
    for (int nens3 = 1; nens3 <= maxens3; ++nens3)
      if (pr_ens[nens3] < 1.e-4) pr_ens[nens3] = 0.;
  }
  // GF:1913-1932
  ierr2 = ierr;
  ierr3 = ierr;
  int k22x = k22, kbconx = 0;
  cup_kbcon(cap_max_increment, 2, k22x, kbconx, heo_cup, heso_cup, hkbo, ierr2, kbmax, po_cup, cap_max, xl, cp,
            c.ztexec, c.zqexec, ktf);
  cup_kbcon(cap_max_increment, 3, k22x, kbconx, heo_cup, heso_cup, hkbo, ierr3, kbmax, po_cup, cap_max, xl, cp,
            c.ztexec, c.zqexec, ktf);
  // GF:1937-1943
  cup_forcing_ens_3d(closure_n, xland1, aa0, aa1, xaa0_ens, mbdt_ens1, dtime, ierr, ierr2, ierr3, xf_ens, axx,
                     c.mconv, omeg, k22, pr_ens, kbcon, ichoice);
  // GF:1945-1963
  for (int k = 1; k <= ktf; ++k) {
    if (ierr == 0) {
      subt_ens[k] = dsubt[k];
      subq_ens[k] = dsubq[k];
      dellat_ens[k] = dellat[k];
      dellaq_ens[k] = dellaq[k];
      dellaqc_ens[k] = dellaqc[k];
      pwo_ens[k] = pwo[k] + edt * pwdo[k];
    } else {
      subt_ens[k] = 0.;
      subq_ens[k] = 0.;
      dellat_ens[k] = 0.;
      dellaq_ens[k] = 0.;
      dellaqc_ens[k] = 0.;
      pwo_ens[k] = 0.;
    }
  }
  // GF:1968-1981
  cup_output_ens_3d(xf_ens, ierr, dellat_ens, dellaq_ens, dellaqc_ens, subt_ens, subq_ens, o.subt, o.subq, o.outt,
                    o.outq, o.outqc, zuo, o.sub_mas, o.pre, pwo_ens, o.xmb, ktop, pr_ens, sig, ktf);
  if (ierr == 0) o.pre = maxr(o.pre, 0.);
  (void)closure_n;
  (void)bu;
  (void)zd;  // zd(i, jmin) = 0.2 (GF:1365): set, never read
}

// CUP_gf_sh's outputs (GF:4151-4160)
struct GfShallowOut {
  real outqc[KM], outt[KM], outq[KM], cupclws[KM];
  real xmb_out = 0;
  int kbcon = 0, ktop = 0, k22 = 0, ierr = 0;
};

// GF:4151-5085 for column c (tn = tshall, qo = qshall, po = p = p2d). On entry o holds cu_grell_freitas' initial
// values (outqcs, outts, outqs, cupclws = 0; k22s, kbcons, ktops, ierrs = 0); pre (pret) is not written.
CPPMPAS_HD inline void cup_gf_sh(const GfColumn& c, real dtime, real xl, real rv, real cp, real g, int ichoice,
                                 GfShallowOut& o) {
  constexpr int use_excess = 0;
  const int ktf = c.ktf;
  const real* zo = c.zo;
  const real* t = c.t2d;
  const real* q = c.q2d;
  const real* tn = c.tshall;
  const real* qo = c.qshall;
  const real* po = c.p2d;
  const real* p = c.p2d;
  const real z1 = c.ter11, psur = c.psur;
  int& ierr = o.ierr;
  int& kbcon = o.kbcon;
  int& ktop = o.ktop;
  int& k22 = o.k22;
  const int kpbl = c.kpbl;

  real entr_rate_2d[KM] = {}, he[KM] = {}, hes[KM] = {}, qes[KM] = {}, z[KM] = {};
  real heo[KM] = {}, heso[KM] = {}, qeso[KM] = {};
  real xhe[KM] = {}, xhes[KM] = {}, xqes[KM] = {}, xz[KM] = {}, xt[KM] = {}, xq[KM] = {};
  real qes_cup[KM] = {}, q_cup[KM] = {}, he_cup[KM] = {}, hes_cup[KM] = {}, z_cup[KM] = {}, p_cup[KM] = {},
       gamma_cup[KM] = {}, t_cup[KM] = {};
  real qeso_cup[KM] = {}, qo_cup[KM] = {}, heo_cup[KM] = {}, heso_cup[KM] = {}, zo_cup[KM] = {}, po_cup[KM] = {},
       gammao_cup[KM] = {}, tn_cup[KM] = {};
  real xqes_cup[KM] = {}, xq_cup[KM] = {}, xhe_cup[KM] = {}, xhes_cup[KM] = {}, xz_cup[KM] = {}, xt_cup[KM] = {};
  real dby[KM] = {}, hc[KM] = {}, zu[KM] = {};
  real dbyo[KM] = {}, qco[KM] = {}, hco[KM] = {}, qrco[KM] = {}, zuo[KM] = {};
  real xdby[KM] = {}, xhc[KM] = {}, xzu[KM] = {};
  real cd[KM] = {}, dellah[KM] = {}, dellaq[KM] = {}, dellat[KM] = {}, dellaqc[KM] = {}, dsubt[KM] = {},
       dsubh[KM] = {}, dsubq[KM] = {};
  real up_massentr[KM] = {}, up_massdetr[KM] = {}, up_massentro[KM] = {}, up_massdetro[KM] = {};

  const real zustart = .1, zufinal = 1.;  // GF:4356-4360 (levadj, power_entr = 2., day: see the header)
  o.xmb_out = 0.;  // GF:4361-4367
  real xland1 = 1.;
  if (c.xland > 1.5) xland1 = 0.;
  (void)xland1;
  const real cap_max_increment = 25.;
  const real entr_rate = .2 / 200.;  // GF:4372-4373
  const int tun_lim = 7;
  for (int k = 1; k <= ktf; ++k) {  // GF:4378-4388
    up_massentro[k] = 0.;
    up_massdetro[k] = 0.;
    z[k] = zo[k];
    xz[k] = zo[k];
    qrco[k] = 0.;
    cd[k] = 1. * entr_rate;
    dellaqc[k] = 0.;
  }
  const real cap_maxs = 25.;  // GF:4394-4412 (depth_min = 50 feeds nothing)
  int kbmax = 1;
  real aa0 = 0., aa1 = 0.;
  real cap_max = cap_maxs;
  const real zkbmax = 4000.;
  // GF:4416-4437
  cup_env(z, qes, he, hes, t, q, p, ierr, ktf);
  cup_env(zo, qeso, heo, heso, tn, qo, po, ierr, ktf);
  cup_env_clev(t, qes, q, he, hes, z, p, qes_cup, q_cup, he_cup, hes_cup, z_cup, p_cup, gamma_cup, t_cup, psur, ierr,
               z1, xl, rv, cp, ktf);
  cup_env_clev(tn, qeso, qo, heo, heso, zo, po, qeso_cup, qo_cup, heo_cup, heso_cup, zo_cup, po_cup, gammao_cup,
               tn_cup, psur, ierr, z1, xl, rv, cp, ktf);
  if (ierr == 0) {  // GF:4438-4454
    if (c.aaeq < -0.1) ierr = 20;
    for (int k = 1; k <= ktf; ++k)
      if (zo_cup[k] > zkbmax + z1) {
        kbmax = k;
        break;
      }
    kbmax = mini(kbmax, ktf - 4);
  }
  // GF:4461-4480
  k22 = cup_maximi(heo_cup, 3, kbmax, ierr);
  if (kpbl > 5) cap_max = po_cup[kpbl];
  if (ierr == 0) {
    if (k22 > kbmax) ierr = 2;
    if (kpbl > 5) {
      k22 = kpbl;
      ierr = 0;
    }
  }
  // GF:4485-4505
  real hkb = 0., hkbo = 0., qkbo = 0.;
  if (ierr == 0) {
    hkb = he_cup[k22] + real(use_excess) * (xl * c.zqexec + cp * c.ztexec);
    hkbo = heo_cup[k22] + real(use_excess) * (xl * c.zqexec + cp * c.ztexec);
    qkbo = qo_cup[k22] + real(use_excess) * (xl * c.zqexec);
    for (int k = 1; k <= k22; ++k) {
      hkb = maxr(hkb, he_cup[k]);
      hkbo = maxr(hkbo, heo_cup[k]);
      qkbo = maxr(qkbo, qo_cup[k]);
    }
  }
  cup_kbcon(cap_max_increment, 5, k22, kbcon, heo_cup, heso_cup, hkbo, ierr, kbmax, po_cup, cap_max, xl, cp,
            c.ztexec, c.zqexec, ktf);  // GF:4506-4510
  // GF:4514-4551
  real frh = 0., dh = 0., dz = 0., zuhe = 0.;
  if (ierr == 0) {
    if (kbcon > ktf - 4) {
      ierr = 231;
    } else {
      if (kbcon + tun_lim + 1 > ktf)
        kernel_fail("cup_gf_sh: kbcon + 7 >= ktf, MPAS's loop GF:4538 would run past the column (undefined)");
      for (int k = 1; k <= ktf; ++k) {
        frh = minr(qo_cup[k] / qeso_cup[k], 1.);
        entr_rate_2d[k] = entr_rate * (1.3 - frh);
        cd[k] = entr_rate_2d[k];
      }
      zuhe = zustart;
      const int kstart = 1;
      frh = (zufinal - zustart) / ((real(kbcon) * real(kbcon)) - (real(kstart) * real(kstart)));
      dh = zuhe - frh * (real(kstart) * real(kstart));
      for (int k = kstart; k <= kbcon - 1; ++k) {
        dz = z_cup[k + 1] - z_cup[k];
        cd[k] = 0.;
        entr_rate_2d[k] = ((frh * (real((k + 1)) * real((k + 1))) + dh) / zuhe - 1. + cd[k] * dz) / dz;
        zuhe = zuhe + entr_rate_2d[k] * dz * zuhe - cd[k] * dz * zuhe;
      }
      frh = -(0.1 - zuhe) / ((real(kbcon + tun_lim) * real(kbcon + tun_lim)) - (real(kbcon - 1) * real(kbcon - 1)));
      dh = zuhe + frh * (real(kbcon) * real(kbcon));
      for (int k = kbcon; k <= kbcon + tun_lim; ++k) {
        dz = z_cup[k + 1] - z_cup[k];
        cd[k] = -((-frh * (real((k + 1)) * real((k + 1))) + dh) / zuhe - 1. - entr_rate_2d[k] * dz) / dz;
        zuhe = zuhe + entr_rate_2d[k] * dz * zuhe - cd[k] * dz * zuhe;
      }
      for (int k = kbcon + tun_lim + 1; k <= ktf; ++k) {
        entr_rate_2d[k] = 0.;
        cd[k] = 0.;
      }
    }
  }
  // GF:4555-4580
  for (int k = 1; k <= ktf; ++k) {
    hc[k] = 0.;
    dby[k] = 0.;
    hco[k] = 0.;
    dbyo[k] = 0.;
  }
  if (ierr == 0) {
    for (int k = 1; k <= kbcon - 1; ++k) {
      hc[k] = hkb;
      hco[k] = hkbo;
      qco[k] = qkbo;
    }
    const int k = kbcon;
    hc[k] = hkb;
    qco[k] = qkbo;
    dby[kbcon] = hkb - hes_cup[k];
    hco[k] = hkbo;
    dbyo[kbcon] = hkbo - heso_cup[k];
    const real trash = qeso_cup[k] + (1. / xl) * (gammao_cup[k] / (1. + gammao_cup[k])) * dbyo[k];
    qrco[k] = maxr(qco[k] - trash, 0.);
  }
  // GF:4583-4660
  if (ierr == 0) {
    zu[1] = zustart;
    zuo[1] = zustart;
    for (int k = 2; k <= ktf - 1; ++k) {
      dz = zo_cup[k] - zo_cup[k - 1];
      up_massentro[k - 1] = entr_rate_2d[k - 1] * dz * zuo[k - 1];
      up_massdetro[k - 1] = cd[k - 1] * dz * zuo[k - 1];
      zuo[k] = zuo[k - 1] + up_massentro[k - 1] - up_massdetro[k - 1];
      if (zuo[k] < 0.05) {
        zuo[k] = .05;
        up_massdetro[k - 1] = zuo[k - 1] - .05 + up_massentro[k - 1];
        cd[k - 1] = up_massdetro[k - 1] / dz / zuo[k - 1];
      }
      zu[k] = zuo[k];
      up_massentr[k - 1] = up_massentro[k - 1];
      up_massdetr[k - 1] = up_massdetro[k - 1];
    }
    for (int k = kbcon + 1; k <= ktf - 1; ++k) {
      hc[k] = (hc[k - 1] * zu[k - 1] - .5 * up_massdetr[k - 1] * hc[k - 1] + up_massentr[k - 1] * he[k - 1]) /
              (zu[k - 1] - .5 * up_massdetr[k - 1] + up_massentr[k - 1]);
      dby[k] = hc[k] - hes_cup[k];
      hco[k] = (hco[k - 1] * zuo[k - 1] - .5 * up_massdetro[k - 1] * hco[k - 1] + up_massentro[k - 1] * heo[k - 1]) /
               (zuo[k - 1] - .5 * up_massdetro[k - 1] + up_massentro[k - 1]);
      dbyo[k] = hco[k] - heso_cup[k];
    }
    for (int k = kbcon + 1; k <= ktf; ++k)
      if (dbyo[k] < 0) {
        ktop = k - 1;
        break;  // go to 41
      }
    if (ktop < kbcon + 1) {
      ierr = 5;
    } else if (ktop > ktf - 2) {
      ierr = 5;
    } else {
      for (int k = kbcon + 1; k <= ktop; ++k) {
        const real trash = qeso_cup[k] + (1. / xl) * (gammao_cup[k] / (1. + gammao_cup[k])) * dbyo[k];
        qco[k] = (qco[k - 1] * zuo[k - 1] - .5 * up_massdetr[k - 1] * qco[k - 1] + up_massentr[k - 1] * qo[k - 1]) /
                 (zuo[k - 1] - .5 * up_massdetr[k - 1] + up_massentr[k - 1]);
        qrco[k] = maxr(qco[k] - trash, 0.);
        o.cupclws[k] = qrco[k] * .5;
      }
      for (int k = ktop + 1; k <= ktf; ++k) {
        hc[k] = hes_cup[k];
        hco[k] = heso_cup[k];
        dby[k] = 0.;
        dbyo[k] = 0.;
        zu[k] = 0.;
        zuo[k] = 0.;
        cd[k] = 0.;
        entr_rate_2d[k] = 0.;
        up_massentr[k] = 0.;
        up_massdetr[k] = 0.;
        up_massentro[k] = 0.;
        up_massdetro[k] = 0.;
      }
    }
  }
  // GF:4665-4680
  aa0 = cup_up_aa0(z, zu, dby, gamma_cup, t_cup, kbcon, ktop, ierr, ktf);
  aa1 = cup_up_aa0(zo, zuo, dbyo, gammao_cup, tn_cup, kbcon, ktop, ierr, ktf);
  if (ierr == 0)
    if (aa1 == 0.) ierr = 17;
  // GF:4689-4808
  for (int k = 1; k <= ktf; ++k) {
    dellah[k] = 0.;
    dsubt[k] = 0.;
    dsubh[k] = 0.;
    dellaq[k] = 0.;
    dsubq[k] = 0.;
  }
  if (ierr == 0) {
    real dp = 100. * (po_cup[1] - po_cup[2]);
    dsubt[1] = 0.;
    dsubq[1] = 0.;
    for (int k = 2; k <= ktop; ++k) {
      real subin = 0., subdown = 0.;
      const real entupk = 0.;
      real detupk = 0.;
      real entup = up_massentro[k];
      real detup = up_massdetro[k];
      if (k == ktop) {
        detupk = zuo[ktop];
        subin = 0.;
        subdown = 0.;
        entup = 0.;
        detup = 0.;
      }
      // totmas (GF:4762-4775): written to unit 0 when |totmas| > 1e-6, no value set
      dp = 100. * (po_cup[k] - po_cup[k + 1]);
      dellah[k] = (detup * .5 * (hco[k + 1] + hco[k]) - entup * heo[k] + subin * heo_cup[k + 1] -
                   subdown * heo_cup[k] + detupk * (hco[ktop] - heo_cup[ktop]) - entupk * heo_cup[k22]) *
                  g / dp;
      dellaq[k] = (detup * .5 * (qco[k + 1] + qco[k] - qrco[k + 1] - qrco[k]) - entup * qo[k] +
                   subin * qo_cup[k + 1] - subdown * qo_cup[k] + detupk * (qco[ktop] - qrco[ktop] - qo_cup[ktop]) -
                   entupk * qo_cup[k22]) *
                  g / dp;
      if (k < ktop) {
        dsubt[k] = (zuo[k + 1] * heo_cup[k + 1] - zuo[k] * heo_cup[k]) * g / dp;
        dsubq[k] = (zuo[k + 1] * qo_cup[k + 1] - zuo[k] * qo_cup[k]) * g / dp;
      }
    }
  }
  // GF:4812-4829
  for (int k = 1; k <= ktf - 1; ++k) {
    dellaqc[k] = 0.;
    if (ierr == 0) {
      if (k == ktop - 0) dellaqc[k] = .01 * zuo[ktop] * qrco[ktop] * 9.81 / (po_cup[k] - po_cup[k + 1]);
      if (k < ktop && k > kbcon) {
        dz = zo_cup[k + 1] - zo_cup[k];
        dellaqc[k] = .01 * 9.81 * up_massdetro[k] * .5 * (qrco[k] + qrco[k + 1]) / (po_cup[k] - po_cup[k + 1]);
      }
      dellaqc[k] = maxr(dellaqc[k], 0.);
    }
  }
  // GF:4833-4859
  const real mbdt = 3.e-4;
  for (int k = 1; k <= ktf; ++k) {
    dellat[k] = 0.;
    if (ierr == 0) {
      dsubh[k] = dsubt[k];
      dellaq[k] = dellaq[k] + dellaqc[k];
      dellaqc[k] = 0.;
      xhe[k] = (dsubt[k] + dellah[k]) * mbdt + heo[k];
      xq[k] = (dsubq[k] + dellaq[k]) * mbdt + qo[k];
      dellat[k] = (1. / cp) * (dellah[k] - xl * dellaq[k]);
      dsubt[k] = (1. / cp) * (dsubt[k] - xl * dsubq[k]);
      xt[k] = (dellat[k] + dsubt[k]) * mbdt + tn[k];
      if (xq[k] <= 0.) xq[k] = 1.E-08;
    }
  }
  real xhkb = 0.;
  if (ierr == 0) {
    xhkb = hkbo + (dsubh[k22] + dellah[k22]) * mbdt;
    xhe[ktf] = heo[ktf];
    xq[ktf] = qo[ktf];
    xt[ktf] = tn[ktf];
    if (xq[ktf] <= 0.) xq[ktf] = 1.E-08;
  }
  // GF:4863-4874
  cup_env(xz, xqes, xhe, xhes, xt, xq, po, ierr, ktf);
  cup_env_clev(xt, xqes, xq, xhe, xhes, xz, po, xqes_cup, xq_cup, xhe_cup, xhes_cup, xz_cup, po_cup, gamma_cup,
               xt_cup, psur, ierr, z1, xl, rv, cp, ktf);
  // GF:4886-4927
  for (int k = 1; k <= ktf; ++k) {
    xhc[k] = 0.;
    xdby[k] = 0.;
  }
  if (ierr == 0) {
    for (int k = 1; k <= kbcon - 1; ++k) xhc[k] = xhkb;
    const int k = kbcon;
    xhc[k] = xhkb;
    xdby[kbcon] = xhkb - xhes_cup[k];
  }
  if (ierr == 0) {
    for (int k = 1; k <= ktf; ++k) xzu[k] = zuo[k];
    for (int k = kbcon + 1; k <= ktop; ++k) {
      xhc[k] = (xhc[k - 1] * xzu[k - 1] - .5 * up_massdetro[k - 1] * xhc[k - 1] + up_massentro[k - 1] * xhe[k - 1]) /
               (xzu[k - 1] - .5 * up_massdetro[k - 1] + up_massentro[k - 1]);
      xdby[k] = xhc[k] - xhes_cup[k];
    }
    for (int k = ktop + 1; k <= ktf; ++k) {
      xhc[k] = xhes_cup[k];
      xdby[k] = 0.;
      xzu[k] = 0.;
    }
  }
  // GF:4932-4935
  const real xaa0 = cup_up_aa0(xz, xzu, xdby, gamma_cup, xt_cup, kbcon, ktop, ierr, ktf);
  // GF:4939-5079: shallow forcing
  real xmb = 0.;
  real xff_shal[10];
  for (int n = 1; n <= 9; ++n) xff_shal[n] = 0.;
  if (ierr == 0) {
    const real xmbmax = 0.1;
    real xkshal = (xaa0 - aa1) / mbdt;
    if (xkshal >= 0.) xkshal = +1.e6;
    if (xkshal > -1.e-4 && xkshal < 0.) xkshal = -1.e-4;
    xff_shal[1] = maxr(-(aa1 - aa0) / (xkshal * dtime), 0.);
    xff_shal[1] = minr(xmbmax, xff_shal[1]);
    xff_shal[2] = maxr(-(aa1 - aa0) / (xkshal * dtime), 0.);
    xff_shal[2] = minr(xmbmax, xff_shal[2]);
    xff_shal[3] = maxr(-(aa1 - aa0) / (xkshal * dtime), 0.);
    xff_shal[3] = minr(xmbmax, xff_shal[3]);
    if (aa1 <= 0) {
      xff_shal[1] = 0.;
      xff_shal[2] = 0.;
      xff_shal[3] = 0.;
    }
    if (aa1 - aa0 <= 0.) {
      xff_shal[1] = 0.;
      xff_shal[2] = 0.;
      xff_shal[3] = 0.;
    }
    real blqe = 0.;  // GF:4964-4975
    real trash = 0.;
    if (k22 < kpbl + 1) {
      for (int k = 1; k <= kbcon - 1; ++k) blqe = blqe + 100. * c.dhdt[k] * (po_cup[k] - po_cup[k + 1]) / g;
      trash = maxr((hc[kbcon] - he_cup[kbcon]), 1.e1);
      xff_shal[7] = maxr(blqe / trash, 0.);
      xff_shal[7] = minr(xmbmax, xff_shal[7]);
    } else {
      xff_shal[7] = 0.;
    }
    real xff_max = 0;  // GF:4977-4988 (mpas)
    if (xkshal < -1.1e-04) {
      xff_shal[4] = maxr(-aa0 / (xkshal * c.tscl_kf), 0.);
      xff_shal[4] = minr(xmbmax, xff_shal[4]);
      xff_shal[5] = xff_shal[4];
      xff_shal[6] = xff_shal[4];
    } else {
      xff_shal[4] = 0.;
      xff_shal[5] = 0.;
      xff_shal[6] = 0.;
    }
    xff_shal[8] = xff_shal[7];
    xff_shal[9] = xff_shal[7];
    real fsum = 0.;  // GF:5006-5018
    for (int k = 1; k <= 9; ++k) {
      if (ichoice > 0) {
        xmb = xmb + xff_shal[ichoice];
      } else if (ichoice == 0) {
        xmb = xmb + xff_shal[k];
      } else if (ichoice < 0) {
        xff_max = xff_shal[1];  // maxval(xff_shal(1:9))
        for (int n = 2; n <= 9; ++n)
          if (xff_shal[n] > xff_max) xff_max = xff_shal[n];
        xmb = xmb + xff_max;
      }
      fsum = fsum + 1.;
    }
    xmb = minr(xmbmax, xmb / fsum);
    if (xmb == 0.) ierr = 22;
    if (xmb < 0.) ierr = 21;
  }
  if (ierr != 0) {  // GF:5028-5037
    k22 = 0;
    kbcon = 0;
    ktop = 0;
    xmb = 0;
    for (int k = 1; k <= ktf; ++k) {
      o.outt[k] = 0.;
      o.outq[k] = 0.;
      o.outqc[k] = 0.;
    }
  } else {  // GF:5038-5078
    real trash = 0.;
    for (int k = 2; k <= ktop; ++k) trash = maxr(trash, 86400. * (dsubt[k] + dellat[k]) * xmb);
    if (trash > 100.) xmb = xmb * 100. / trash;
    trash = 0.;
    for (int k = 2; k <= ktop; ++k) trash = minr(trash, 86400. * (dsubt[k] + dellat[k]) * xmb);
    if (trash < -100.) xmb = -xmb * 100. / trash;
    for (int k = 2; k <= ktop; ++k) {
      trash = q[k] + (dsubq[k] + dellaq[k]) * xmb * dtime;
      if (trash < 1.e-12) {
        trash = (1.e-12 - q[k]) / ((dsubq[k] + dellaq[k]) * dtime);
        xmb = (1.e-12 - q[k]) / ((dsubq[k] + dellaq[k]) * dtime);
      }
    }
    o.xmb_out = xmb;
    for (int k = 2; k <= ktop; ++k) {
      o.outt[k] = (dsubt[k] + dellat[k]) * xmb;
      o.outq[k] = (dsubq[k] + dellaq[k]) * xmb;
    }
  }
}

// One column of cu_grell_freitas (GF:54-453), the arrays of column i gathered into levels 1..kte
struct GfDriverIn {
  int kte = 0, itimestep = 0;
  real dt = 0, dxCell = 0, areaCell = 0, ht = 0, xland = 0, gsw = 0, hfx = 0, qfx = 0, p8w1 = 0;
  int kpbl = 0;
  real u[KM], v[KM], t[KM], q[KM], p[KM], pi[KM], rho[KM], dz8w[KM], w[KM];  // w(i,k,j) at k = 1..kte
  real rqvften[KM], rthften[KM], rthblten[KM], rqvblten[KM], rthraten[KM];
  real xlv = 0, cp = 0, g = 0, r_v = 0;
  int ichoice_deep = 0, ichoice_shallow = 0, ishallow_g3 = 0;
};
struct GfDriverOut {
  real hbot = 0, htop = 0, xmb_total = 0, raincv = 0, pratec = 0, xmb_shallow = 0;
  int k22_shallow = 0, kbcon_shallow = 0, ktop_shallow = 0;
  real gdc[KM], gdc2[KM], rthcuten[KM], rqvcuten[KM], rqccuten[KM], rqicuten[KM];
};

CPPMPAS_HD inline void cu_grell_freitas_column(const GfDriverIn& in, GfDriverOut& out) {
  const int kts = 1, kte = in.kte, ktf = kte;  // ktf = min(kte, kde - 1)
  // GF:165-199
  out.hbot = real(kte);
  out.htop = real(kts);
  out.xmb_total = 0.;
  out.raincv = 0.;
  out.pratec = 0.;
  out.k22_shallow = 0;
  out.kbcon_shallow = 0;
  out.ktop_shallow = 0;
  out.xmb_shallow = 0.;
  for (int k = kts; k <= kte; ++k) {  // do k = jts, kte
    out.gdc[k] = 0.;
    out.gdc2[k] = 0.;
  }
  // GF:204-274 (apri_*, edti_out, gswi, pri_ens, xfi_ens: not read by MPAS's CUP_gf)
  GfColumn c;
  c.ktf = ktf;
  c.tscl_kf = in.dxCell / 25.;
  c.dx = in.dxCell;
  c.ter11 = maxr(in.ht, 0.);
  c.zo[kts] = c.ter11 + 0.5 * in.dz8w[1];
  for (int k = kts + 1; k <= ktf; ++k) c.zo[k] = c.zo[k - 1] + 0.5 * (in.dz8w[k - 1] + in.dz8w[k]);
  c.psur = in.p8w1 * .01;
  c.kpbl = in.kpbl;
  c.xland = in.xland;
  c.ccn = 1500.;
  real cuten = 0., pret = 0.;
  c.aaeq = 0.;
  c.zqexec = 0.;
  c.ztexec = 0.;
  GfDeepOut d;
  GfShallowOut s;
  d.kbcon = 0;
  d.ktop = 0;
  d.xmb = 0.;
  s.ierr = 0;
  s.kbcon = 0;
  s.ktop = 0;
  s.k22 = 0;
  s.xmb_out = 0.;
  // GF:276-312
  for (int k = kts; k <= ktf; ++k) {
    c.us[k] = in.u[k];
    c.vs[k] = in.v[k];
    c.rhoi[k] = in.rho[k];
    c.t2d[k] = in.t[k];
    c.q2d[k] = in.q[k];
    if (c.q2d[k] < 1.e-08) c.q2d[k] = 1.e-08;
    c.tn[k] = c.t2d[k] + (in.rthften[k] + in.rthraten[k] + in.rthblten[k]) * in.pi[k] * in.dt;
    c.qo[k] = c.q2d[k] + (in.rqvften[k] + in.rqvblten[k]) * in.dt;
    if (c.tn[k] < 200.) c.tn[k] = c.t2d[k];
    if (c.qo[k] < 1.e-08) c.qo[k] = 1.e-08;
    const real phh = in.p[k];
    c.po[k] = phh * .01;
    c.p2d[k] = c.po[k];
    d.cupclw[k] = 0.;
    d.outq[k] = 0.;
    d.outqc[k] = 0.;
    d.outt[k] = 0.;
    d.sub_mas[k] = 0.;
    d.subq[k] = 0.;
    d.subt[k] = 0.;
    c.tshall[k] = c.t2d[k] + in.rthblten[k] * in.pi[k] * in.dt;
    c.qshall[k] = c.q2d[k] + in.rqvblten[k] * in.dt;
    c.dhdt[k] = in.cp * in.rthblten[k] * in.pi[k] + in.xlv * in.rqvblten[k];
    s.cupclws[k] = 0.;
    s.outqc[k] = 0.;
    s.outq[k] = 0.;
    s.outt[k] = 0.;
  }
  // GF:224-233, 315-331: omeg and mconv (the ens4 copies are identical)
  for (int k = kts; k <= ktf; ++k) c.omeg[k] = 0.;
  c.mconv = 0.;
  for (int k = kts + 1; k <= ktf; ++k) c.omeg[k] = -(in.g * 0.5 * (in.rho[k] + in.rho[k - 1]) * in.w[k]);
  for (int k = kts + 1; k <= ktf; ++k) {
    const real dq = (c.q2d[k] - c.q2d[k - 1]);
    c.mconv = c.mconv + c.omeg[k] * dq / in.g;
  }
  if (c.mconv < 0.) c.mconv = 0.;

  // GF:390-401
  d.pre = pret;
  cup_gf(c, in.dt, in.xlv, in.r_v, in.cp, in.g, in.ichoice_deep, d);
  pret = d.pre;
  neg_check(d.subt, d.subq, in.dt, c.q2d, d.outq, d.outt, d.outqc, pret, ktf);

  // GF:404-419
  if (in.ishallow_g3 == 1) {
    cup_gf_sh(c, in.dt, in.xlv, in.r_v, in.cp, in.g, in.ichoice_shallow, s);
    out.xmb_shallow = s.xmb_out;
    out.k22_shallow = s.k22;
    out.kbcon_shallow = s.kbcon;
    out.ktop_shallow = s.ktop;
  }

  // GF:421-431
  out.xmb_total = d.xmb;
  if (pret > 0.) {
    cuten = 1.;
    out.pratec = pret;
    out.raincv = pret * in.dt;
    if (real(d.ktop) > out.htop) out.htop = real(d.ktop) + .001;
    if (real(d.kbcon) < out.hbot) out.hbot = real(d.kbcon) + .001;
  }
  // GF:434-449
  for (int k = kts; k <= kte; ++k) {
    out.rthcuten[k] = (s.outt[k] + (d.subt[k] + d.outt[k]) * cuten) / in.pi[k];
    out.rqvcuten[k] = s.outq[k] + (d.subq[k] + d.outq[k]) * cuten;
    if (c.t2d[k] < tcrit) {
      out.rqccuten[k] = 0.;
      out.rqicuten[k] = s.outqc[k] + d.outqc[k] * cuten;
      out.gdc2[k] = s.cupclws[k] + d.cupclw[k] * cuten;
    } else {
      out.rqicuten[k] = 0.;
      out.rqccuten[k] = s.outqc[k] + d.outqc[k] * cuten;
      out.gdc[k] = s.cupclws[k] + d.cupclw[k] * cuten;
    }
  }
}

}  // namespace cppmpas::phys::cu_gf
