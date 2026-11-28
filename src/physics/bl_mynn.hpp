// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// MYNN planetary boundary layer as MPAS-Atmosphere v8.4.2 runs it: MMM-physics a4baf7f (the reference executable's
// physics_mmm; cppmpas-data `inputs` mmm-physics-a4baf7f.tar.gz) bl_mynn.F90 (BM), bl_mynn_subroutines.F90 (BS) and
// mynn_shared.F90 (MS), statement by statement, one column per call (CPPMPAS_HD, the k loops inside the thread).
// Written from MPAS's Fortran: CPPWRF ports MYNN-EDMF v6 (WRF 4.8), whose routines differ from this v5 line, so no
// CPPWRF code is shared.
//
// Ported (the path of MPAS's defaults, Registry.xml config_mynn_*): bl_mynn_init (BM:24-122), bl_mynn_run for one
// column (BM:130-1240), and from BS: mym_initialize (473-633), mym_level2 (677-788), mym_length (808-1202, mixlength 2),
// mym_turbulence (1577-2108, closure < 3), mym_predict (2156-2526, closure <= 2.5), mym_condensation (2562-2981,
// cloudpdf 2), mynn_tendencies (2987-4092, mixqt 0), moisture_check (4095-4182), retrieve_exchange_coeffs
// (4316-4340), tridiag2 (4380-4412), get_pblh (4475-4615), dmp_mf (4637-5785, momentum and TKE transport off, no
// chemistry), condensation_edmf (5789-5846), scale_aware (6273-6344), phim (6346-6395), phih (6398-6444); MS
// esat_blend, qsat_blend (waterice absent), xl_blend. Refused (BlMynnOptions::check): bl_mynn_mixlength 0 / 1 (BS
// 880-1057 and boulac_length), bl_mynn_closure > 2.5 (the prognostic tsq / qsq / cov of level 3), bl_mynn_cloudpdf
// other than 2 (BS:2649-2772, erf), bl_mynn_stfunc 0 (BM:831-837), bl_mynn_topdown (topdown_cloudrad), bl_mynn_edmf
// off and bl_mynn_edmf_mom / _tke / _dd on (DMP_mf's momentum and TKE plumes, ddmf_jpl), bl_mynn_mixqt,
// bl_mynn_cloudmix off, bl_mynn_mixscalars off, bl_mynn_tkebudget, bl_mynn_tkeadvect, bl_mynn_output, spp_pbl > 0,
// icloud_bl > 0, and the number concentrations of the aerosol-aware Thompson scheme (flag_qnc / qnwfa / qnifa /
// qnbca) and ozone (flag_qoz).
//
// Reference object code (bin/atmosphere_model sha256 524a1cf4...): bl_mynn.F90, bl_mynn_subroutines.F90 and
// mynn_shared.F90 compiled standalone with the reference flags (gfortran 14.2 -O3 -ffree-line-length-none
// -fconvert=big-endian -ffree-form -fdefault-real-8 -fdefault-double-8, no -Dmpas, as the build log) are
// instruction-identical to the executable's (normalised objdump, 0 differing lines in all 24 routines). No FMA, no
// ymm register. Every untyped literal is a double (-fdefault-real-8). The libm calls (-g build, objdump -l):
//  * libmvec lanes (pairs of iterations from the loop start; the odd last iteration and a 1-trip loop call the scalar
//    function: lane(i, n)):
//      BM:975 qke1(k)**1.5 _ZGVbN2vv_pow and BM:978 exp(-10000./max(p1(k),1.)) _ZGVbN2v_exp, the dissipative heating
//        loop k = kts..kte-1;
//      BS:5033 l**d in DMP_mf's second plume loop i = 1..8 (four pairs, _ZGVbN2vv_pow); the first loop (BS:5007,
//        the sum cn) calls the scalar pow, so the two loops' l**d may differ.
//  * scalar libm (FMA variants, core/libm.hpp): pow (mym_initialize 586, 589, 601; mym_length 1110, 1132, 1178;
//    mym_condensation 2635-2636; DMP_mf 5046, 5071; condensation_edmf 5812; scale_aware 6298, 6328; phim / phih),
//    tanh (mym_length 1135 = 1162 = 1182, computed once; get_pblh 4594; DMP_mf 4977-5015), exp (mym_condensation
//    2863, 2919-2923; DMP_mf 5716-5720), atan (2856, 5661, phim / phih), log10 (2833, 2841), log (phim / phih).
//    Powers with a real exponent stay pow calls (x**0.5 in phih included); x**2. and integer powers are multiplies.
//  * folded constants (no call in the object; the .rodata value is used): mym_initialize's (b1*pmz)**(2.0/3.0)
//    (BS:529, 0x4020a402fcc79298) and phh*b2/(b1*pmz)**(1.0/3.0) (BS:537, 0x4014cd03bbf9773e; pmz = phh = 1. are
//    constants there), and exp(-1.) of BS:2869-2870 (0x3fd78b56362cef38).
//  * MAX(a, b) and MIN(a, b) are maxsd / minsd with a as the destination: a > b ? a : b, a < b ? a : b.
#include <cstdio>

#include "core/backend.hpp"
#include "core/libm.hpp"
#include "core/precision.hpp"
#include "physics/bl_mynn_common.hpp"

namespace cppmpas::phys::mynn {

inline constexpr int kMaxLevels = 100;      // nVertLevels <= kMaxLevels
inline constexpr int KM = kMaxLevels + 3;   // 1-based levels 1..kte+1
inline constexpr int NUP = 8;               // BS:4718

CPPMPAS_HD inline real fmx(real a, real b) { return a > b ? a : b; }  // MAX(a, b): maxsd, a the destination
CPPMPAS_HD inline real fmn(real a, real b) { return a < b ? a : b; }  // MIN(a, b): minsd
CPPMPAS_HD inline real fabs_(real a) { return __builtin_fabs(a); }  // ABS: andpd (sign cleared)
// the libmvec rule (CPPWRF 7fac5bf src/core/backend.hpp:197): lanes for the pairs from the loop start, the scalar
// libm for an odd last iteration and a 1-trip loop; i is the 0-based iteration, n the trip count
CPPMPAS_HD inline bool lane(int i, int n) { return n > 1 && !((n & 1) && i == n - 1); }

// ---- mynn_shared.F90 ----
inline constexpr real j0 = .611583699e03, j1 = .444606896e02, j2 = .143177157e01, j3 = .264224321e-1,
                      j4 = .299291081e-3, j5 = .203154182e-5, j6 = .702620698e-8, j7 = .379534310e-11,
                      j8 = -.321582393e-13;
inline constexpr real k0 = .609868993e03, k1 = .499320233e02, k2 = .184672631e01, k3 = .402737184e-1,
                      k4 = .565392987e-3, k5 = .521693933e-5, k6 = .307839583e-7, k7 = .105785160e-9,
                      k8 = .161444444e-12;
CPPMPAS_HD inline real esat_l(real xc) {
  return j0 + xc * (j1 + xc * (j2 + xc * (j3 + xc * (j4 + xc * (j5 + xc * (j6 + xc * (j7 + xc * j8)))))));
}
CPPMPAS_HD inline real esat_i(real xc) {
  return k0 + xc * (k1 + xc * (k2 + xc * (k3 + xc * (k4 + xc * (k5 + xc * (k6 + xc * (k7 + xc * k8)))))));
}
// MS:46-63
CPPMPAS_HD inline real esat_blend(real t, real t0c, real tic) {
  const real xc = fmx(t - t0c, -80.);
  if (t >= t0c) return esat_l(xc);
  if (t <= tic) return esat_i(xc);
  const real esl = esat_l(xc);
  const real esi = esat_i(xc);
  const real chi = (273.16 - t) / 20.16;
  return (1. - chi) * esl + chi * esi;
}
// MS:67-97 (waterice absent: wrt = 'b')
CPPMPAS_HD inline real qsat_blend(real t, real t0c, real tic, real p) {
  const real xc = fmx(t - t0c, -80.);
  if (t >= t0c) {
    const real esl = esat_l(xc);
    return 0.622 * esl / (p - esl);
  }
  if (t <= tic) {
    const real esi = esat_i(xc);
    return 0.622 * esi / (p - esi);
  }
  const real esl = esat_l(xc);
  const real esi = esat_i(xc);
  const real rslf = 0.622 * esl / (p - esl);
  const real rsif = 0.622 * esi / (p - esi);
  const real chi = (t0c - t) / (t0c - tic);
  return (1. - chi) * rslf + chi * rsif;
}
// MS:101-118
CPPMPAS_HD inline real xl_blend(real t, real t0c, real tic, real cice, real cliq, real cpv, real xls, real xlv) {
  if (t >= t0c) return xlv + (cpv - cliq) * (t - t0c);
  if (t <= tic) return xls + (cpv - cice) * (t - t0c);
  const real xlvt = xlv + (cpv - cliq) * (t - t0c);
  const real xlst = xls + (cpv - cice) * (t - t0c);
  const real chi = (t0c - t) / (t0c - tic);
  return (1. - chi) * xlvt + chi * xlst;
}

// ---- bl_mynn_subroutines.F90 ----
// Column arrays are 1-based (index k = Fortran k), kts = 1.

// BS:677-788
CPPMPAS_HD inline void mym_level2(const Const& C, int kte, const real* dz, const real* u, const real* v,
                                  const real* thl, const real* thetav, const real* qw, const real* ql,
                                  const real* vt, const real* vq, real* dtl, real* dqw, real* dtv, real* gm,
                                  real* gh, real* sm, real* sh) {
  (void)thetav;
  (void)ql;
  const int kts = 1;
  for (int k = kts + 1; k <= kte; ++k) {  // BS:725
    const real dzk = 0.5 * (dz[k] + dz[k - 1]);
    const real afk = dz[k] / (dz[k] + dz[k - 1]);
    const real abk = 1.0 - afk;
    real duz = (u[k] - u[k - 1]) * (u[k] - u[k - 1]) + (v[k] - v[k - 1]) * (v[k] - v[k - 1]);
    duz = duz / (dzk * dzk);
    const real dtz = (thl[k] - thl[k - 1]) / (dzk);
    const real dqz = (qw[k] - qw[k - 1]) / (dzk);
    const real vtt = 1.0 + vt[k] * abk + vt[k - 1] * afk;  // BS:734
    const real vqq = C.tv0 + vq[k] * abk + vq[k - 1] * afk;
    const real dtq = vtt * dtz + vqq * dqz;
    dtl[k] = dtz;  // BS:740
    dqw[k] = dqz;
    dtv[k] = dtq;
    gm[k] = duz;  // BS:747
    gh[k] = -dtq * C.gtr;
    const real ri = -gh[k] / fmx(duz, 1.0e-10);  // BS:751
    real a2fac;
    if (CKmod == 1) {  // BS:754
      a2fac = 1. / (1. + fmx(ri, 0.0));
    } else {
      a2fac = 1.;
    }
    const real rfc = g1 / (g1 + g2);  // BS:760
    const real f1 = b1 * (g1 - c1) + 3.0 * a2 * a2fac * (1.0 - c2) * (1.0 - c5) + 2.0 * a1 * (3.0 - 2.0 * c2);
    const real f2 = b1 * (g1 + g2) - 3.0 * a1 * (1.0 - c2);
    const real rf1 = b1 * (g1 - c1) / f1;
    const real rf2 = b1 * g1 / f2;
    const real smc = a1 / (a2 * a2fac) * f1 / f2;
    const real shc = 3.0 * (a2 * a2fac) * (g1 + g2);
    const real ri1 = 0.5 / smc;  // BS:769
    const real ri2 = rf1 * smc;
    const real ri3 = 4.0 * rf2 * smc - 2.0 * ri2;
    const real ri4 = ri2 * ri2;
    const real rf = fmn(ri1 * (ri + ri2 - __builtin_sqrt(ri * ri - ri3 * ri + ri4)), rfc);  // BS:775
    sh[k] = shc * (rfc - rf) / (1.0 - rf);
    sm[k] = smc * (rf1 - rf) / (rf2 - rf) * sh[k];
  }
}

// BS:808-1202, bl_mynn_mixlength = 2 (CASE (2), BS:1059-1192)
CPPMPAS_HD inline void mym_length(const Const& C, int kte, real xland, const real* dz, real dx, const real* zw,
                                  real rmo, real flt, real fltv, real flq, const real* vt, const real* vq,
                                  const real* u1, const real* v1, const real* qke, const real* dtv, real* el,
                                  real zi, const real* theta, real* qkw, real Psig_bl, const real* cldfra_bl1D,
                                  const real* edmf_w1, const real* edmf_a1) {
  (void)xland;
  (void)dx;
  (void)flt;
  (void)flq;
  (void)vt;
  (void)vq;
  (void)theta;
  const int kts = 1;
  const real ctau = 1000.;  // BS:873
  real qtke[KM];
  // BS:1061-1069 (Uonset and Ugrid are not used)
  (void)u1;
  (void)v1;
  const real cns = 3.5, alp1 = 0.22, alp2 = 0.30, alp3 = 2.0, alp4 = 5.0, alp5 = alp2, alp6 = 50.0;
  const real zi2 = fmx(zi, 300.);  // BS:1073
  real h1 = fmx(0.3 * zi2, 300.);  // BS:1076
  h1 = fmn(h1, 600.);
  const real h2 = h1 * 0.5;
  qtke[kts] = fmx(0.5 * qke[kts], 0.01);  // BS:1080
  qkw[kts] = __builtin_sqrt(fmx(qke[kts], 1.0e-4));
  for (int k = kts + 1; k <= kte; ++k) {  // BS:1083
    const real afk = dz[k] / (dz[k] + dz[k - 1]);
    const real abk = 1.0 - afk;
    qkw[k] = __builtin_sqrt(fmx(qke[k] * abk + qke[k - 1] * afk, 1.0e-3));
    qtke[k] = 0.5 * (qkw[k] * qkw[k]);
  }
  real elt = 1.0e-5;  // BS:1090
  real vsc = 1.0e-5;
  const real PBLH_PLUS_ENT = fmx(zi + h1, 100.);  // BS:1094
  int k = kts + 1;
  real zwk = zw[k];
  while (zwk <= PBLH_PLUS_ENT) {
    // BS:1094-1104 has no bound on k: with the PBL top above the column MPAS reads dz(k), qkw(k) past kte (undefined)
    if (k > kte) {
      kernel_fail("mym_length: zi + h1 is above the model top (MPAS reads dz and qkw past kte there)");
      return;
    }
    const real dzk = 0.5 * (dz[k] + dz[k - 1]);
    const real qdz = fmn(fmx(qkw[k] - qmin, 0.03), 30.0) * dzk;
    elt = elt + qdz * zwk;
    vsc = vsc + qdz;
    k = k + 1;
    zwk = zw[k];
  }
  elt = fmn(fmx(alp1 * elt / vsc, 10.), 400.);  // BS:1106
  const real vflx = fltv;  // BS:1109
  vsc = libm::pow(C.gtr * elt * fmx(vflx, 0.0), onethird);
  el[kts] = 0.0;  // BS:1113
  // the wstar of BS:1132 = 1159 (loop invariant, one pow in the object)
  const real wstar = 1.25 * libm::pow(C.gtr * zi * fmx(vflx, 1.0e-4), onethird);
  for (k = kts + 1; k <= kte; ++k) {  // BS:1116
    zwk = zw[k];
    const real dzk = 0.5 * (dz[k] + dz[k - 1]);
    const real cldavg = 0.5 * (cldfra_bl1D[k - 1] + cldfra_bl1D[k]);
    (void)cldavg;
    // tanh of BS:1135 = 1162 = 1182: one call per level
    const real wt = .5 * libm::tanh((zwk - (zi2 + h1)) / h2) + .5;
    real elb, elf, elb_mf;
    if (dtv[k] > 0.0) {  // BS:1122
      const real bv = fmx(__builtin_sqrt(C.gtr * dtv[k]), 0.001);
      elb_mf = fmx(alp2 * qkw[k], alp6 * edmf_a1[k - 1] * edmf_w1[k - 1]) / bv *
               (1.0 + alp3 * __builtin_sqrt(vsc / (bv * elt)));
      elb = fmn(fmx(alp5 * qkw[k], alp6 * edmf_a1[k] * edmf_w1[k]) / bv, zwk);
      real tau_cloud = fmn(fmx(ctau * wstar / C.grav, 30.), 150.);  // BS:1133
      tau_cloud = tau_cloud * (1. - wt) + 50. * wt;                 // BS:1136
      elf = fmn(fmx(tau_cloud * __builtin_sqrt(fmn(qtke[k], 40.)), alp6 * edmf_a1[k] * edmf_w1[k] / bv), zwk);
    } else {
      real tau_cloud = fmn(fmx(ctau * wstar / C.grav, 50.), 200.);  // BS:1160
      tau_cloud = tau_cloud * (1. - wt) + fmx(100., dzk * 0.25) * wt;  // BS:1164
      elb = fmn(tau_cloud * __builtin_sqrt(fmn(qtke[k], 40.)), zwk);
      elf = elb;
      elb_mf = elb;
    }
    elf = elf / (1. + (elf / 800.));  // BS:1171
    elb_mf = fmx(elb_mf, 0.01);
    real els;
    if (rmo > 0.0) {  // BS:1175
      els = C.karman * zwk / (1.0 + cns * fmn(zwk * rmo, zmax));
    } else {
      els = C.karman * zwk * libm::pow(1.0 - alp4 * zwk * rmo, 0.2);
    }
    el[k] = __builtin_sqrt(els * els / (1. + (els * els / (elt * elt)) + (els * els / (elb_mf * elb_mf))));
    el[k] = el[k] * (1. - wt) + elf * wt;  // BS:1186
    const real el_les = fmn(els / (1. + (els / 12.)), elb_mf);
    el[k] = el[k] * Psig_bl + (1. - Psig_bl) * el_les;
  }
}

// BS:473-633
CPPMPAS_HD inline void mym_initialize(const Const& C, int kte, real xland, const real* dz, real dx, const real* zw,
                                      const real* u, const real* v, const real* thl, const real* qw, real zi,
                                      const real* theta, const real* thetav, real* sh, real* sm, real ust, real rmo,
                                      real* el, real* qke, real* tsq, real* qsq, real* cov, real Psig_bl,
                                      const real* cldfra_bl1D, const real* edmf_w1, const real* edmf_a1,
                                      bool INITIALIZE_QKE) {
  const int kts = 1;
  real ql[KM], pdk[KM], pdt[KM], pdq[KM], pdc[KM], dtl[KM], dqw[KM], dtv[KM], gm[KM], gh[KM], qkw[KM], vt[KM],
      vq[KM];
  // pmz = phh = 1. and flt = fltv = flq = 0. (BS:504-505, saved locals never changed)
  const real flt = 0., fltv = 0., flq = 0.;
  for (int k = kts; k <= kte; ++k) {  // BS:511
    ql[k] = 0.0;
    vt[k] = 0.0;
    vq[k] = 0.0;
  }
  mym_level2(C, kte, dz, u, v, thl, thetav, qw, ql, vt, vq, dtl, dqw, dtv, gm, gh, sm, sh);  // BS:518
  el[kts] = 0.0;  // BS:526
  // (b1*pmz)**(2.0/3.0) and phh*b2/(b1*pmz)**(1.0/3.0): folded by gfortran (no pow in the object), the .rodata values
  constexpr real b1pmz_2_3 = __builtin_bit_cast(real, 0x4020a402fcc79298ull);  // 8.320335292207616
  constexpr real phm0 = __builtin_bit_cast(real, 0x4014cd03bbf9773eull);       // 5.20020955762976
  if (INITIALIZE_QKE) {
    qke[kts] = 1.5 * (ust * ust) * b1pmz_2_3;  // BS:529
    for (int k = kts + 1; k <= kte; ++k)       // BS:530
      qke[k] = qke[kts] * fmx((ust * 700. - zw[k]) / (fmx(ust, 0.01) * 700.), 0.01);
  }
  real phm = phm0;  // BS:537
  tsq[kts] = phm * ((flt / ust) * (flt / ust));
  qsq[kts] = phm * ((flq / ust) * (flq / ust));
  cov[kts] = phm * (flt / ust) * (flq / ust);
  for (int k = kts + 1; k <= kte; ++k) {  // BS:542
    const real vkz = C.karman * zw[k];
    el[k] = vkz / (1.0 + vkz / 100.0);
    tsq[k] = 0.0;
    qsq[k] = 0.0;
    cov[k] = 0.0;
  }
  const int lmax = 5;  // BS:554
  for (int l = 1; l <= lmax; ++l) {
    mym_length(C, kte, xland, dz, dx, zw, rmo, flt, fltv, flq, vt, vq, u, v, qke, dtv, el, zi, theta, qkw, Psig_bl,
               cldfra_bl1D, edmf_w1, edmf_a1);
    for (int k = kts + 1; k <= kte; ++k) {  // BS:572
      const real elq = el[k] * qkw[k];
      pdk[k] = elq * (sm[k] * gm[k] + sh[k] * gh[k]);
      pdt[k] = elq * sh[k] * (dtl[k] * dtl[k]);
      pdq[k] = elq * sh[k] * (dqw[k] * dqw[k]);
      pdc[k] = elq * sh[k] * dtl[k] * dqw[k];
    }
    const real vkz = C.karman * 0.5 * dz[kts];  // BS:582
    const real elv = 0.5 * (el[kts + 1] + el[kts]) / vkz;
    if (INITIALIZE_QKE) {  // BS:586 (b1*pmz = 24.)
      const real m = fmx(ust, 0.02);
      qke[kts] = 1.0 * (m * m) * libm::pow(b1 * 1. * elv, 2.0 / 3.0);
    }
    phm = 1. * b2 / libm::pow(b1 * 1. / (elv * elv), 1.0 / 3.0);  // BS:589
    tsq[kts] = phm * ((flt / ust) * (flt / ust));
    qsq[kts] = phm * ((flq / ust) * (flq / ust));
    cov[kts] = phm * (flt / ust) * (flq / ust);
    for (int k = kts + 1; k <= kte - 1; ++k) {  // BS:594
      const real b1l = b1 * 0.25 * (el[k + 1] + el[k]);
      const real tmpq = fmn(fmx(b1l * (pdk[k + 1] + pdk[k]), qkemin), 125.);
      if (INITIALIZE_QKE) qke[k] = libm::pow(tmpq, twothirds);
      real b2l;
      if (qke[k] <= 0.0) {
        b2l = 0.0;
      } else {
        b2l = b2 * (b1l / b1) / __builtin_sqrt(qke[k]);
      }
      tsq[k] = b2l * (pdt[k + 1] + pdt[k]);
      qsq[k] = b2l * (pdq[k + 1] + pdq[k]);
      cov[k] = b2l * (pdc[k + 1] + pdc[k]);
    }
  }
  if (INITIALIZE_QKE) {  // BS:622
    qke[kts] = 0.5 * (qke[kts] + qke[kts + 1]);
    qke[kte] = qke[kte - 1];
  }
  tsq[kte] = tsq[kte - 1];
  qsq[kte] = qsq[kte - 1];
  cov[kte] = cov[kte - 1];
}

// BS:1577-2108 (closure < 3: the level-2.5 branch; bl_mynn_tkebudget = .false.; spp_pbl = 0)
CPPMPAS_HD inline void mym_turbulence(const Const& C, int kte, real xland, real closure, const real* dz, real dx,
                                      const real* zw, const real* u, const real* v, const real* thl,
                                      const real* thetav, const real* ql, const real* qw, const real* qke,
                                      const real* tsq, const real* qsq, const real* cov, const real* vt,
                                      const real* vq, real rmo, real flt, real fltv, real flq, real zi,
                                      const real* theta, real* sh, real* sm, real* el, real* dfm, real* dfh,
                                      real* dfq, real* tcd, real* qcd, real* pdk, real* pdt, real* pdq, real* pdc,
                                      real Psig_bl, real Psig_shcu, const real* cldfra_bl1D, const real* edmf_w1,
                                      const real* edmf_a1, const real* TKEprodTD) {
  (void)closure;
  (void)tsq;
  (void)qsq;
  (void)cov;
  (void)Psig_shcu;
  const int kts = 1;
  real qkw[KM], dtl[KM], dqw[KM], dtv[KM], gm[KM], gh[KM];
  const real Prlimit = 5.0;  // BS:1649
  mym_level2(C, kte, dz, u, v, thl, thetav, qw, ql, vt, vq, dtl, dqw, dtv, gm, gh, sm, sh);  // BS:1664
  mym_length(C, kte, xland, dz, dx, zw, rmo, flt, fltv, flq, vt, vq, u, v, qke, dtv, el, zi, theta, qkw, Psig_bl,
             cldfra_bl1D, edmf_w1, edmf_a1);  // BS:1670
  for (int k = kts + 1; k <= kte; ++k) {  // BS:1684
    const real dzk = 0.5 * (dz[k] + dz[k - 1]);
    const real afk = dz[k] / (dz[k] + dz[k - 1]);
    (void)afk;
    const real elsq = el[k] * el[k];
    real q3sq = qkw[k] * qkw[k];
    const real q2sq = b1 * elsq * (sm[k] * gm[k] + sh[k] * gh[k]);
    sh[k] = fmx(sh[k], 1e-5);  // BS:1694 (sh20, sm20 are not used)
    real duz = (u[k] - u[k - 1]) * (u[k] - u[k - 1]) + (v[k] - v[k - 1]) * (v[k] - v[k - 1]);  // BS:1697
    duz = duz / (dzk * dzk);
    const real ri = -gh[k] / fmx(duz, 1.0e-10);
    real a2fac;
    if (CKmod == 1) {
      a2fac = 1. / (1. + fmx(ri, 0.0));
    } else {
      a2fac = 1.;
    }
    const real gmel = gm[k] * elsq;  // BS:1717
    const real ghel = gh[k] * elsq;
    real dlsq = elsq;  // BS:1736
    if (q3sq / dlsq < -gh[k]) q3sq = -dlsq * gh[k];
    real qdiv, e1, e2, e3, e4, eden;
    if (q3sq < q2sq) {  // BS:1739
      qdiv = __builtin_sqrt(q3sq / q2sq);
      sh[k] = sh[k] * qdiv;
      sm[k] = sm[k] * qdiv;
      const real qd2 = qdiv * qdiv;
      e1 = q3sq - e1c * ghel * a2fac * qd2;  // BS:1776
      e2 = q3sq - e2c * ghel * a2fac * qd2;
      e3 = e1 + e3c * ghel * (a2fac * a2fac) * qd2;
      e4 = e1 - e4c * ghel * a2fac * qd2;
      eden = e2 * e4 + e3 * e5c * gmel * qd2;
      eden = fmx(eden, 1.0e-20);
    } else {
      e1 = q3sq - e1c * ghel * a2fac;  // BS:1792
      e2 = q3sq - e2c * ghel * a2fac;
      e3 = e1 + e3c * ghel * (a2fac * a2fac);
      e4 = e1 - e4c * ghel * a2fac;
      eden = e2 * e4 + e3 * e5c * gmel;
      eden = fmx(eden, 1.0e-20);
      qdiv = 1.0;
      sm[k] = q3sq * a1 * (e3 - 3.0 * c1 * e4) / eden;  // BS:1801
      sh[k] = q3sq * (a2 * a2fac) * (e2 + 3.0 * c1 * e5c * gmel) / eden;
    }
    const real sh25max = 4., sh25min = 0.0;  // BS:1817
    if (sh[k] > sh25max) sh[k] = sh25max;  // BS:1840
    if (sh[k] < sh25min) sh[k] = sh25min;
    const real shb = fmx(sh[k], 0.002);  // BS:1853
    sm[k] = fmn(sm[k], Prlimit * shb);
    const real gamt = 0.0, gamq = 0.0, gamv = 0.0;  // BS:2000 (closure < 3)
    const real cldavg = 0.5 * (cldfra_bl1D[k - 1] + cldfra_bl1D[k]);  // BS:2007
    if (edmf_a1[k] > 0.001 || cldavg > 0.02) {
      sm[k] = fmx(sm[k], 0.03 * fmn(10. * edmf_a1[k] * edmf_w1[k], 1.0));
      sh[k] = fmx(sh[k], 0.03 * fmn(10. * edmf_a1[k] * edmf_w1[k], 1.0));
      sm[k] = fmx(sm[k], 0.05 * fmn(cldavg, 1.0));
      sh[k] = fmx(sh[k], 0.05 * fmn(cldavg, 1.0));
    }
    const real elq = el[k] * qkw[k];  // BS:2017
    const real elh = elq * qdiv;
    pdk[k] = elq * (sm[k] * gm[k] + sh[k] * gh[k] + gamv) + 0.5 * TKEprodTD[k];  // BS:2022
    pdt[k] = elh * (sh[k] * dtl[k] + gamt) * dtl[k];
    pdq[k] = elh * (sh[k] * dqw[k] + gamq) * dqw[k];
    pdc[k] = elh * (sh[k] * dtl[k] + gamt) * dqw[k] * 0.5 + elh * (sh[k] * dqw[k] + gamq) * dtl[k] * 0.5;
    tcd[k] = elq * gamt;  // BS:2032
    qcd[k] = elq * gamq;
    dfm[k] = elq * sm[k] / dzk;  // BS:2036
    dfh[k] = elq * sh[k] / dzk;
    dfq[k] = dfm[k];
  }
  dfm[kts] = 0.0;  // BS:2079
  dfh[kts] = 0.0;
  dfq[kts] = 0.0;
  tcd[kts] = 0.0;
  qcd[kts] = 0.0;
  tcd[kte] = 0.0;
  qcd[kte] = 0.0;
  for (int k = kts; k <= kte - 1; ++k) {  // BS:2089
    const real dzk = dz[k];
    tcd[k] = (tcd[k + 1] - tcd[k]) / (dzk);
    qcd[k] = (qcd[k + 1] - qcd[k]) / (dzk);
  }
}

// BS:4380-4412
CPPMPAS_HD inline void tridiag2(int n, const real* a, const real* b, const real* c, const real* d, real* x) {
  real cp[KM], dp[KM];
  cp[1] = c[1] / b[1];
  dp[1] = d[1] / b[1];
  for (int i = 2; i <= n; ++i) {
    const real m = b[i] - cp[i - 1] * a[i];
    cp[i] = c[i] / m;
    dp[i] = (d[i] - dp[i - 1] * a[i]) / m;
  }
  x[n] = dp[n];
  for (int i = n - 1; i >= 1; --i) x[i] = dp[i] - cp[i] * x[i + 1];
}

// BS:2156-2526 (closure = 2.5: no prognostic qsq / tsq / cov; bl_mynn_tkebudget = .false.)
CPPMPAS_HD inline void mym_predict(const Const& C, int kte, real closure, real delt, const real* dz, real ust,
                                   real flt, real flq, real pmz, real phh, const real* el, const real* dfq,
                                   const real* rho, real* pdk, real* pdt, real* pdq, real* pdc, real* qke, real* tsq,
                                   real* qsq, real* cov, const real* s_aw, const real* s_awqke,
                                   bool bl_mynn_edmf_tke) {
  (void)closure;
  const int kts = 1;
  real qkw[KM], bp[KM], rp[KM], df3q[KM], dtz[KM], a[KM], b[KM], c[KM], d[KM], x[KM], rhoinv[KM], rhoz[KM],
      kqdz[KM], kmdz[KM];
  real onoff;
  if (!bl_mynn_edmf_tke) {  // BS:2200
    onoff = 0.0;
  } else {
    onoff = 1.0;
  }
  const real vkz = C.karman * 0.5 * dz[kts];  // BS:2207
  for (int k = kts; k <= kte; ++k) {  // BS:2211
    qkw[k] = __builtin_sqrt(fmx(qke[k], 0.0));
    df3q[k] = Sqfac * dfq[k];
    dtz[k] = delt / dz[k];
  }
  rhoz[kts] = rho[kts];  // BS:2221
  rhoinv[kts] = 1. / rho[kts];
  kqdz[kts] = rhoz[kts] * df3q[kts];
  kmdz[kts] = rhoz[kts] * dfq[kts];
  for (int k = kts + 1; k <= kte; ++k) {
    rhoz[k] = (rho[k] * dz[k - 1] + rho[k - 1] * dz[k]) / (dz[k - 1] + dz[k]);
    rhoz[k] = fmx(rhoz[k], 1E-4);
    rhoinv[k] = 1. / fmx(rho[k], 1E-4);
    kqdz[k] = rhoz[k] * df3q[k];
    kmdz[k] = rhoz[k] * dfq[k];
  }
  rhoz[kte + 1] = rhoz[kte];  // BS:2232
  kqdz[kte + 1] = rhoz[kte + 1] * df3q[kte];
  kmdz[kte + 1] = rhoz[kte + 1] * dfq[kte];
  for (int k = kts + 1; k <= kte - 1; ++k) {  // BS:2237
    kqdz[k] = fmx(kqdz[k], 0.5 * s_aw[k]);
    kqdz[k] = fmx(kqdz[k], -0.5 * (s_aw[k] - s_aw[k + 1]));
    kmdz[k] = fmx(kmdz[k], 0.5 * s_aw[k]);
    kmdz[k] = fmx(kmdz[k], -0.5 * (s_aw[k] - s_aw[k + 1]));
  }
  const real pdk1 = 2.0 * (ust * ust * ust) * pmz / (vkz);  // BS:2245
  const real phm = 2.0 / ust * phh / (vkz);
  const real pdt1 = phm * (flt * flt);
  const real pdq1 = phm * (flq * flq);
  const real pdc1 = phm * flt * flq;
  (void)pdt1;
  (void)pdq1;
  (void)pdc1;
  pdk[kts] = pdk1 - pdk[kts + 1];  // BS:2252
  pdt[kts] = pdt[kts + 1];
  pdq[kts] = pdq[kts + 1];
  pdc[kts] = pdc[kts + 1];
  for (int k = kts; k <= kte - 1; ++k) {  // BS:2263
    const real b1l = b1 * 0.5 * (el[k + 1] + el[k]);
    bp[k] = 2. * qkw[k] / b1l;
    rp[k] = pdk[k + 1] + pdk[k];
  }
  for (int k = kts; k <= kte - 1; ++k) {  // BS:2275
    a[k] = -dtz[k] * kqdz[k] * rhoinv[k] + 0.5 * dtz[k] * rhoinv[k] * s_aw[k] * onoff;
    b[k] = 1. + dtz[k] * (kqdz[k] + kqdz[k + 1]) * rhoinv[k] +
           0.5 * dtz[k] * rhoinv[k] * (s_aw[k] - s_aw[k + 1]) * onoff + bp[k] * delt;
    c[k] = -dtz[k] * kqdz[k + 1] * rhoinv[k] - 0.5 * dtz[k] * rhoinv[k] * s_aw[k + 1] * onoff;
    d[k] = rp[k] * delt + qke[k] + dtz[k] * rhoinv[k] * (s_awqke[k] - s_awqke[k + 1]) * onoff;
  }
  a[kte] = 0.;  // BS:2311
  b[kte] = 1.;
  c[kte] = 0.;
  d[kte] = qke[kte];
  tridiag2(kte, a, b, c, d, x);
  for (int k = kts; k <= kte; ++k) {  // BS:2319
    qke[k] = fmx(x[k], 1.e-4);
    qke[k] = fmn(qke[k], 150.);
  }
  for (int k = kts; k <= kte - 1; ++k) {  // BS:2390 (closure <= 2.5: level 2 diagnostic)
    real b2l;
    if (qkw[k] <= 0.0) {
      b2l = 0.0;
    } else {
      b2l = b2 * 0.25 * (el[k + 1] + el[k]) / qkw[k];
    }
    qsq[k] = b2l * (pdq[k + 1] + pdq[k]);
  }
  qsq[kte] = qsq[kte - 1];
  for (int k = kts; k <= kte - 1; ++k) {  // BS:2505 (closure < 3)
    real b2l;
    if (qkw[k] <= 0.0) {
      b2l = 0.0;
    } else {
      b2l = b2 * 0.25 * (el[k + 1] + el[k]) / qkw[k];
    }
    tsq[k] = b2l * (pdt[k + 1] + pdt[k]);
    cov[k] = b2l * (pdc[k + 1] + pdc[k]);
  }
  tsq[kte] = tsq[kte - 1];
  cov[kte] = cov[kte - 1];
}

// BS:2562-2981, bl_mynn_cloudpdf = 2 (CASE (2, -2), BS:2774-2955), spp_pbl as passed (0 in MPAS)
CPPMPAS_HD inline void mym_condensation(const Const& C, int kte, real dx, const real* dz, const real* zw, real xland,
                                        const real* thl, const real* qw, const real* qv, const real* qc,
                                        const real* qi, const real* qs, const real* p, const real* exner,
                                        const real* tsq, const real* qsq, const real* cov, const real* Sh,
                                        const real* el, real* qc_bl1D, real* qi_bl1D, real* cldfra_bl1D, real PBLH1,
                                        real HFX1, real* vt, real* vq, const real* th, real* sgm, real rmo,
                                        int spp_pbl, const real* rstoch_col) {
  (void)dx;
  (void)zw;
  (void)thl;
  (void)qv;
  (void)tsq;
  (void)cov;
  (void)Sh;
  (void)el;
  (void)HFX1;
  (void)rmo;
  const int kts = 1;
  const real qpct_sfc = 0.025, qpct_pbl = 0.030, qpct_trp = 0.040, rhcrit = 0.83, rhmax = 1.02;  // BS:2601-2605
  real alp[KM], bet[KM], a[KM], b[KM], q1[KM], RH[KM];
  // BS:2632-2643: the tropopause level
  int k;
  for (k = kte - 3; k >= kts; --k) {
    const real theta1 = th[k];
    const real theta2 = th[k + 2];
    const real ht1 = 44307.692 * (1.0 - libm::pow(p[k] / 101325., 0.190));
    const real ht2 = 44307.692 * (1.0 - libm::pow(p[k + 2] / 101325., 0.190));
    if ((((theta2 - theta1) / (ht2 - ht1)) < 10. / 1500.) && (ht1 < 19000.) && (ht1 > 4000.)) break;  // goto 86
  }
  const int k_tropo = k + 2 > kts + 2 ? k + 2 : kts + 2;  // MAX(kts+2, k+2)
  const real pblh2 = fmx(10., PBLH1);  // BS:2778 MAX(10._kind_phys,pblh1)
  real zagl = 0.;
  real dzm1 = 0.;
  // exp(-1.): folded by gfortran (BS:2869-2870, no call), the .rodata value
  constexpr real em1 = __builtin_bit_cast(real, 0x3fd78b56362cef38ull);
  for (k = kts; k <= kte - 1; ++k) {  // BS:2781
    zagl = zagl + 0.5 * (dz[k] + dzm1);
    dzm1 = dz[k];
    const real t = th[k] * exner[k];  // BS:2785
    const real xl = xl_blend(t, C.t0c, tice, C.cice, C.cliq, C.cpv, C.xls, C.xlv);
    const real qsat_tk = qsat_blend(t, C.t0c, tice, p[k]);
    RH[k] = fmx(fmn(rhmax, qw[k] / fmx(1.E-10, qsat_tk)), 0.001);
    const real dqsl = qsat_tk * C.ep_2 * C.xlv / (C.r_d * (t * t));  // BS:2791
    alp[k] = 1.0 / (1.0 + dqsl * C.xlvcp);
    bet[k] = dqsl * exner[k];
    const real rsl = xl * qsat_tk / (C.r_v * (t * t));  // BS:2795
    const real cpm = C.cp + qw[k] * C.cpv;
    a[k] = 1. / (1. + xl * rsl / cpm);
    b[k] = a[k] * rsl;
    const real qw_pert = qw[k] + qw[k] * 0.5 * rstoch_col[k] * real(spp_pbl);  // BS:2802
    const real qmq = qw_pert - qsat_tk;
    const real r3sq = fmx(qsq[k], 0.0);  // BS:2809
    sgm[k] = __builtin_sqrt(r3sq);
    sgm[k] = fmn(sgm[k], qsat_tk * 0.666);
    real wt = fmx(500. - fmx(dz[k] - 100., 0.0), 0.0) / 500.;  // BS:2817
    sgm[k] = sgm[k] + sgm[k] * 0.2 * (1.0 - wt);
    real qpct = qpct_pbl * wt + qpct_trp * (1.0 - wt);  // BS:2821
    qpct = fmn(qpct, fmx(qpct_sfc, qpct_pbl * zagl / 500.));
    sgm[k] = fmx(sgm[k], qsat_tk * qpct);
    q1[k] = qmq / sgm[k];  // BS:2825
    real rh_hack = RH[k];
    (void)rh_hack;
    const real wt2 = fmn(fmx(zagl - pblh2, 0.0) / 300., 1.0);
    if ((qi[k] + qs[k]) > 1.e-9 && (zagl > pblh2)) {  // BS:2832
      rh_hack = fmn(rhmax, rhcrit + wt2 * 0.045 * (9.0 + libm::log10(qi[k] + qs[k])));
      RH[k] = fmx(RH[k], rh_hack);
      const real q1_rh = -3. + 3. * (RH[k] - rhcrit) / (1. - rhcrit);
      q1[k] = fmx(q1_rh, q1[k]);
    }
    if (qc[k] > 1.e-6 && (zagl > pblh2)) {  // BS:2840
      rh_hack = fmn(rhmax, rhcrit + wt2 * 0.08 * (6.0 + libm::log10(qc[k])));
      RH[k] = fmx(RH[k], rh_hack);
      const real q1_rh = -3. + 3. * (RH[k] - rhcrit) / (1. - rhcrit);
      q1[k] = fmx(q1_rh, q1[k]);
    }
    real q1k = q1[k];  // BS:2848
    cldfra_bl1D[k] = fmx(0., fmn(1., 0.5 + 0.36 * libm::atan(1.8 * (q1[k] + 0.2))));  // BS:2856
    const real maxqc = fmx(qw[k] - qsat_tk, 0.0);  // BS:2861
    real ql_water, ql_ice;
    if (q1k < 0.) {
      ql_water = sgm[k] * libm::exp(1.2 * q1k - 1.);
      ql_ice = sgm[k] * libm::exp(1.2 * q1k - 1.);
    } else if (q1k > 2.) {
      ql_water = fmn(sgm[k] * q1k, maxqc);
      ql_ice = sgm[k] * q1k;
    } else {
      ql_water = fmn(sgm[k] * (em1 + 0.66 * q1k + 0.086 * (q1k * q1k)), maxqc);
      ql_ice = sgm[k] * (em1 + 0.66 * q1k + 0.086 * (q1k * q1k));
    }
    if (cldfra_bl1D[k] < 0.001) {  // BS:2878
      ql_ice = 0.0;
      ql_water = 0.0;
      cldfra_bl1D[k] = 0.0;
    }
    const real liq_frac = fmn(1.0, fmx(0.0, (t - tice) / (tliq - tice)));  // BS:2884
    qc_bl1D[k] = liq_frac * ql_water;
    qi_bl1D[k] = (1.0 - liq_frac) * ql_ice;
    if (k >= k_tropo) {  // BS:2890
      cldfra_bl1D[k] = 0.;
      qc_bl1D[k] = 0.;
      qi_bl1D[k] = 0.;
    }
    if ((xland - 1.5) >= 0) {  // BS:2899
      q1k = fmx(q1[k], -2.5);
    } else {
      q1k = fmx(q1[k], -2.0);
    }
    real Fng;
    if (q1k >= 1.0) {  // BS:2916
      Fng = 1.0;
    } else if (q1k >= -1.7 && q1k < 1.0) {
      Fng = libm::exp(-0.4 * (q1k - 1.0));
    } else if (q1k >= -2.5 && q1k < -1.7) {
      Fng = 3.0 + libm::exp(-3.8 * (q1k + 1.7));
    } else {
      Fng = fmn(23.9 + libm::exp(-1.6 * (q1k + 2.5)), 60.);
    }
    real cfmax = fmn(cldfra_bl1D[k], 0.6);  // BS:2926
    const real zsl = fmn(fmx(25., 0.1 * pblh2), 100.);
    wt = fmn(zagl / zsl, 1.0);
    cfmax = cfmax * wt;
    const real bb = b[k] * t / th[k];  // BS:2932
    const real qww = 1. + 0.61 * qw[k];  // BS:2938
    const real alpha = 0.61 * th[k];
    const real beta = (th[k] / t) * (xl / C.cp) - 1.61 * th[k];
    vt[k] = qww - cfmax * beta * bb * Fng - 1.;
    vq[k] = alpha + cfmax * beta * a[k] * Fng - C.tv0;
    const real fac_damp = fmn(zagl * 0.0025, 1.0);  // BS:2950
    const real rd = fmx(0.0, (RH[k] - 0.92)) / 0.145;
    const real cld_factor = 1.0 + fac_damp * fmn(rd * rd, 0.37);
    cldfra_bl1D[k] = fmn(1., cld_factor * cldfra_bl1D[k]);
  }
  (void)alp;  // alp and bet feed only cloudpdf 0 / 1 (BS:2792-2793 sets them)
  (void)bet;
  // BS:2968-2973 (ql is not used: the cloudpdf 2 branch does not set it)
  vt[kte] = vt[kte - 1];
  vq[kte] = vq[kte - 1];
  qc_bl1D[kte] = 0.;
  qi_bl1D[kte] = 0.;
  cldfra_bl1D[kte] = 0.;
}

// BS:4095-4182
CPPMPAS_HD inline void moisture_check(const Const& C, int kte, real delt, const real* dp, const real* exner,
                                      real* qv, real* qc, real* qi, real* qs, real* th, real* dqv, real* dqc,
                                      real* dqi, real* dqs, real* dth) {
  const real qvmin = 1e-20, qcmin = 0.0, qimin = 0.0;
  real dqv2 = 0.;
  for (int k = kte; k >= 1; --k) {  // BS:4124
    const real dqc2 = fmx(0.0, qcmin - qc[k]);
    const real dqi2 = fmx(0.0, qimin - qi[k]);
    const real dqs2 = fmx(0.0, qimin - qs[k]);
    dqc[k] = dqc[k] + dqc2 / delt;  // BS:4130
    dqi[k] = dqi[k] + dqi2 / delt;
    dqs[k] = dqs[k] + dqs2 / delt;
    dqv[k] = dqv[k] - (dqc2 + dqi2 + dqs2) / delt;
    dth[k] = dth[k] + C.xlvcp / exner[k] * (dqc2 / delt) + C.xlscp / exner[k] * ((dqi2 + dqs2) / delt);
    qc[k] = qc[k] + dqc2;  // BS:4137
    qi[k] = qi[k] + dqi2;
    qs[k] = qs[k] + dqs2;
    qv[k] = qv[k] - dqc2 - dqi2 - dqs2;
    th[k] = th[k] + C.xlvcp / exner[k] * dqc2 + C.xlscp / exner[k] * (dqi2 + dqs2);
    dqv2 = fmx(0.0, qvmin - qv[k]);  // BS:4145
    dqv[k] = dqv[k] + dqv2 / delt;
    qv[k] = qv[k] + dqv2;
    if (k != 1) {
      qv[k - 1] = qv[k - 1] - dqv2 * dp[k] / dp[k - 1];
      dqv[k - 1] = dqv[k - 1] - dqv2 * dp[k] / dp[k - 1] / delt;
    }
    qv[k] = fmx(qv[k], qvmin);  // BS:4152
    qc[k] = fmx(qc[k], qcmin);
    qi[k] = fmx(qi[k], qimin);
    qs[k] = fmx(qs[k], qimin);
  }
  if (dqv2 > 1.e-20) {  // BS:4160
    real sum = 0.0;
    for (int k = 1; k <= kte; ++k)
      if (qv[k] > 2.0 * qvmin) sum = sum + qv[k] * dp[k];
    const real aa = dqv2 * dp[1] / fmx(1.e-20, sum);
    if (aa < 0.5) {
      for (int k = 1; k <= kte; ++k) {
        if (qv[k] > 2.0 * qvmin) {
          const real dum = aa * qv[k];
          qv[k] = qv[k] - dum;
          dqv[k] = dqv[k] - dum / delt;
        }
      }
    }
  }
}

// BS:2987-4092 for bl_mynn_mixqt = .false., bl_mynn_cloudmix = .true., FLAG_QNC = FLAG_QNWFA = FLAG_QNIFA =
// FLAG_QNBCA = FLAG_OZONE = .false. (checked by the caller); the mass-flux terms of the solvers as in the source
struct TendIn {
  const real *dz, *rho, *u, *v, *th, *tk, *qv, *qc, *qi, *qs, *qnc, *qni, *p, *exner;
  const real *sqv, *sqc, *sqi, *sqs;
  const real *qnwfa, *qnifa, *qnbca, *ozone;
  const real *tcd, *qcd, *dfm, *dfh, *diss_heat;
  const real *s_aw, *s_awthl, *s_awqt, *s_awqv, *s_awqc, *s_awu, *s_awv, *s_awqnc, *s_awqni, *s_awqnwfa,
      *s_awqnifa, *s_awqnbca, *sd_aw, *sd_awthl, *sd_awqt, *sd_awqv, *sd_awqc, *sd_awu, *sd_awv;
  const real *sub_thl, *sub_sqv, *sub_u, *sub_v, *det_thl, *det_sqv, *det_sqc, *det_u, *det_v;
  real delt, psfc, ust, flt, flq, flqv, flqc, wspd, uoce, voce;
  bool FLAG_QC, FLAG_QI, FLAG_QNI, FLAG_QS, bl_mynn_edmf_mom, bl_mynn_mixscalars;
};
struct TendOut {
  real *thl, *sqw;  // inout (thl is updated)
  real *du, *dv, *dth, *dqv, *dqc, *dqi, *dqs, *dqni, *dqnc, *dqnwfa, *dqnifa, *dqnbca, *dozone;
};
CPPMPAS_HD inline void mynn_tendencies(const Const& C, int kte, const TendIn& I, const TendOut& O) {
  const int kts = 1;
  const real* dz = I.dz;
  const real* rho = I.rho;
  const real* p = I.p;
  const real* s_aw = I.s_aw;
  const real* sd_aw = I.sd_aw;
  const real delt = I.delt;
  real* thl = O.thl;
  real dtz[KM], delp[KM], sqv2[KM], sqc2[KM], sqi2[KM], sqs2[KM], qni2[KM], rhoinv[KM], a[KM], b[KM], c[KM],
      d[KM], x[KM], rhoz[KM], khdz[KM], kmdz[KM];
  const real nonloc = 1.0;  // BS:3081
  const real dztop = .5 * (dz[kte] + dz[kte - 1]);  // BS:3083
  (void)dztop;
  real onoff;
  if (!I.bl_mynn_edmf_mom) {  // BS:3088
    onoff = 0.0;
  } else {
    onoff = 1.0;
  }
  const real rhosfc = I.psfc / (C.r_d * (I.tk[kts] + C.p608 * I.qv[kts]));  // BS:3096
  dtz[kts] = delt / dz[kts];
  rhoz[kts] = rho[kts];
  rhoinv[kts] = 1. / rho[kts];
  khdz[kts] = rhoz[kts] * I.dfh[kts];
  kmdz[kts] = rhoz[kts] * I.dfm[kts];
  delp[kts] = I.psfc - (p[kts + 1] * dz[kts] + p[kts] * dz[kts + 1]) / (dz[kts] + dz[kts + 1]);
  for (int k = kts + 1; k <= kte; ++k) {  // BS:3103
    dtz[k] = delt / dz[k];
    rhoz[k] = (rho[k] * dz[k - 1] + rho[k - 1] * dz[k]) / (dz[k - 1] + dz[k]);
    rhoz[k] = fmx(rhoz[k], 1E-4);
    rhoinv[k] = 1. / fmx(rho[k], 1E-4);
    khdz[k] = rhoz[k] * I.dfh[k];
    kmdz[k] = rhoz[k] * I.dfm[k];
  }
  for (int k = kts + 1; k <= kte - 1; ++k)  // BS:3112
    delp[k] = (p[k] * dz[k - 1] + p[k - 1] * dz[k]) / (dz[k] + dz[k - 1]) -
              (p[k + 1] * dz[k] + p[k] * dz[k + 1]) / (dz[k] + dz[k + 1]);
  delp[kte] = delp[kte - 1];
  rhoz[kte + 1] = rhoz[kte];
  khdz[kte + 1] = rhoz[kte + 1] * I.dfh[kte];
  kmdz[kte + 1] = rhoz[kte + 1] * I.dfm[kte];
  for (int k = kts + 1; k <= kte - 1; ++k) {  // BS:3122
    khdz[k] = fmx(khdz[k], 0.5 * s_aw[k]);
    khdz[k] = fmx(khdz[k], -0.5 * (s_aw[k] - s_aw[k + 1]));
    kmdz[k] = fmx(kmdz[k], 0.5 * s_aw[k]);
    kmdz[k] = fmx(kmdz[k], -0.5 * (s_aw[k] - s_aw[k + 1]));
  }
  const real ust = I.ust, wspd = I.wspd;
  for (int k = kts; k <= kte; ++k) O.dth[k] = 0.0;  // BS:3131

  // u (BS:3137-3192)
  {
    const int k = kts;
    a[k] = -dtz[k] * kmdz[k] * rhoinv[k];
    b[k] = 1. + dtz[k] * (kmdz[k + 1] + rhosfc * (ust * ust) / wspd) * rhoinv[k] -
           0.5 * dtz[k] * rhoinv[k] * s_aw[k + 1] * onoff - 0.5 * dtz[k] * rhoinv[k] * sd_aw[k + 1] * onoff;
    c[k] = -dtz[k] * kmdz[k + 1] * rhoinv[k] - 0.5 * dtz[k] * rhoinv[k] * s_aw[k + 1] * onoff -
           0.5 * dtz[k] * rhoinv[k] * sd_aw[k + 1] * onoff;
    d[k] = I.u[k] + dtz[k] * I.uoce * (ust * ust) / wspd - dtz[k] * rhoinv[k] * I.s_awu[k + 1] * onoff +
           dtz[k] * rhoinv[k] * I.sd_awu[k + 1] * onoff + I.sub_u[k] * delt + I.det_u[k] * delt;
  }
  for (int k = kts + 1; k <= kte - 1; ++k) {
    a[k] = -dtz[k] * kmdz[k] * rhoinv[k] + 0.5 * dtz[k] * rhoinv[k] * s_aw[k] * onoff +
           0.5 * dtz[k] * rhoinv[k] * sd_aw[k] * onoff;
    b[k] = 1. + dtz[k] * (kmdz[k] + kmdz[k + 1]) * rhoinv[k] +
           0.5 * dtz[k] * rhoinv[k] * (s_aw[k] - s_aw[k + 1]) * onoff +
           0.5 * dtz[k] * rhoinv[k] * (sd_aw[k] - sd_aw[k + 1]) * onoff;
    c[k] = -dtz[k] * kmdz[k + 1] * rhoinv[k] - 0.5 * dtz[k] * rhoinv[k] * s_aw[k + 1] * onoff -
           0.5 * dtz[k] * rhoinv[k] * sd_aw[k + 1] * onoff;
    d[k] = I.u[k] + dtz[k] * rhoinv[k] * (I.s_awu[k] - I.s_awu[k + 1]) * onoff -
           dtz[k] * rhoinv[k] * (I.sd_awu[k] - I.sd_awu[k + 1]) * onoff + I.sub_u[k] * delt + I.det_u[k] * delt;
  }
  a[kte] = 0;
  b[kte] = 1.;
  c[kte] = 0.;
  d[kte] = I.u[kte];
  tridiag2(kte, a, b, c, d, x);
  for (int k = kts; k <= kte; ++k) O.du[k] = (x[k] - I.u[k]) / delt;

  // v (BS:3198-3253)
  {
    const int k = kts;
    a[k] = -dtz[k] * kmdz[k] * rhoinv[k];
    b[k] = 1. + dtz[k] * (kmdz[k + 1] + rhosfc * (ust * ust) / wspd) * rhoinv[k] -
           0.5 * dtz[k] * rhoinv[k] * s_aw[k + 1] * onoff - 0.5 * dtz[k] * rhoinv[k] * sd_aw[k + 1] * onoff;
    c[k] = -dtz[k] * kmdz[k + 1] * rhoinv[k] - 0.5 * dtz[k] * rhoinv[k] * s_aw[k + 1] * onoff -
           0.5 * dtz[k] * rhoinv[k] * sd_aw[k + 1] * onoff;
    d[k] = I.v[k] + dtz[k] * I.voce * (ust * ust) / wspd - dtz[k] * rhoinv[k] * I.s_awv[k + 1] * onoff +
           dtz[k] * rhoinv[k] * I.sd_awv[k + 1] * onoff + I.sub_v[k] * delt + I.det_v[k] * delt;
  }
  for (int k = kts + 1; k <= kte - 1; ++k) {
    a[k] = -dtz[k] * kmdz[k] * rhoinv[k] + 0.5 * dtz[k] * rhoinv[k] * s_aw[k] * onoff +
           0.5 * dtz[k] * rhoinv[k] * sd_aw[k] * onoff;
    b[k] = 1. + dtz[k] * (kmdz[k] + kmdz[k + 1]) * rhoinv[k] +
           0.5 * dtz[k] * rhoinv[k] * (s_aw[k] - s_aw[k + 1]) * onoff +
           0.5 * dtz[k] * rhoinv[k] * (sd_aw[k] - sd_aw[k + 1]) * onoff;
    c[k] = -dtz[k] * kmdz[k + 1] * rhoinv[k] - 0.5 * dtz[k] * rhoinv[k] * s_aw[k + 1] * onoff -
           0.5 * dtz[k] * rhoinv[k] * sd_aw[k + 1] * onoff;
    d[k] = I.v[k] + dtz[k] * rhoinv[k] * (I.s_awv[k] - I.s_awv[k + 1]) * onoff -
           dtz[k] * rhoinv[k] * (I.sd_awv[k] - I.sd_awv[k + 1]) * onoff + I.sub_v[k] * delt + I.det_v[k] * delt;
  }
  a[kte] = 0;
  b[kte] = 1.;
  c[kte] = 0.;
  d[kte] = I.v[kte];
  tridiag2(kte, a, b, c, d, x);
  for (int k = kts; k <= kte; ++k) O.dv[k] = (x[k] - I.v[k]) / delt;

  // thl (BS:3258-3321)
  {
    const int k = kts;
    a[k] = -dtz[k] * khdz[k] * rhoinv[k];
    b[k] = 1. + dtz[k] * (khdz[k + 1] + khdz[k]) * rhoinv[k] - 0.5 * dtz[k] * rhoinv[k] * s_aw[k + 1] -
           0.5 * dtz[k] * rhoinv[k] * sd_aw[k + 1];
    c[k] = -dtz[k] * khdz[k + 1] * rhoinv[k] - 0.5 * dtz[k] * rhoinv[k] * s_aw[k + 1] -
           0.5 * dtz[k] * rhoinv[k] * sd_aw[k + 1];
    d[k] = thl[k] + dtz[k] * rhosfc * I.flt * rhoinv[k] + I.tcd[k] * delt - dtz[k] * rhoinv[k] * I.s_awthl[k + 1] -
           dtz[k] * rhoinv[k] * I.sd_awthl[k + 1] + I.diss_heat[k] * delt + I.sub_thl[k] * delt +
           I.det_thl[k] * delt;
  }
  for (int k = kts + 1; k <= kte - 1; ++k) {
    a[k] = -dtz[k] * khdz[k] * rhoinv[k] + 0.5 * dtz[k] * rhoinv[k] * s_aw[k] + 0.5 * dtz[k] * rhoinv[k] * sd_aw[k];
    b[k] = 1. + dtz[k] * (khdz[k] + khdz[k + 1]) * rhoinv[k] + 0.5 * dtz[k] * rhoinv[k] * (s_aw[k] - s_aw[k + 1]) +
           0.5 * dtz[k] * rhoinv[k] * (sd_aw[k] - sd_aw[k + 1]);
    c[k] = -dtz[k] * khdz[k + 1] * rhoinv[k] - 0.5 * dtz[k] * rhoinv[k] * s_aw[k + 1] -
           0.5 * dtz[k] * rhoinv[k] * sd_aw[k + 1];
    d[k] = thl[k] + I.tcd[k] * delt + dtz[k] * rhoinv[k] * (I.s_awthl[k] - I.s_awthl[k + 1]) +
           dtz[k] * rhoinv[k] * (I.sd_awthl[k] - I.sd_awthl[k + 1]) + I.diss_heat[k] * delt +
           I.sub_thl[k] * delt + I.det_thl[k] * delt;
  }
  a[kte] = 0.;
  b[kte] = 1.;
  c[kte] = 0.;
  d[kte] = thl[kte];
  tridiag2(kte, a, b, c, d, x);
  for (int k = kts; k <= kte; ++k) thl[k] = x[k];

  // bl_mynn_mixqt = .false.: sqw2 = sqw (BS:3384-3386), not used further

  // sqc (BS:3393-3444)
  if (I.FLAG_QC) {
    const int k = kts;
    a[k] = -dtz[k] * khdz[k] * rhoinv[k];
    b[k] = 1. + dtz[k] * (khdz[k + 1] + khdz[k]) * rhoinv[k] - 0.5 * dtz[k] * rhoinv[k] * s_aw[k + 1] -
           0.5 * dtz[k] * rhoinv[k] * sd_aw[k + 1];
    c[k] = -dtz[k] * khdz[k + 1] * rhoinv[k] - 0.5 * dtz[k] * rhoinv[k] * s_aw[k + 1] -
           0.5 * dtz[k] * rhoinv[k] * sd_aw[k + 1];
    d[k] = I.sqc[k] + dtz[k] * rhosfc * I.flqc * rhoinv[k] + I.qcd[k] * delt - dtz[k] * rhoinv[k] * I.s_awqc[k + 1] -
           dtz[k] * rhoinv[k] * I.sd_awqc[k + 1] + I.det_sqc[k] * delt;
    for (int kk = kts + 1; kk <= kte - 1; ++kk) {
      const int q = kk;
      a[q] = -dtz[q] * khdz[q] * rhoinv[q] + 0.5 * dtz[q] * rhoinv[q] * s_aw[q] + 0.5 * dtz[q] * rhoinv[q] * sd_aw[q];
      b[q] = 1. + dtz[q] * (khdz[q] + khdz[q + 1]) * rhoinv[q] +
             0.5 * dtz[q] * rhoinv[q] * (s_aw[q] - s_aw[q + 1]) + 0.5 * dtz[q] * rhoinv[q] * (sd_aw[q] - sd_aw[q + 1]);
      c[q] = -dtz[q] * khdz[q + 1] * rhoinv[q] - 0.5 * dtz[q] * rhoinv[q] * s_aw[q + 1] -
             0.5 * dtz[q] * rhoinv[q] * sd_aw[q + 1];
      d[q] = I.sqc[q] + I.qcd[q] * delt + dtz[q] * rhoinv[q] * (I.s_awqc[q] - I.s_awqc[q + 1]) +
             dtz[q] * rhoinv[q] * (I.sd_awqc[q] - I.sd_awqc[q + 1]) + I.det_sqc[q] * delt;
    }
    a[kte] = 0.;
    b[kte] = 1.;
    c[kte] = 0.;
    d[kte] = I.sqc[kte];
    tridiag2(kte, a, b, c, d, sqc2);
  } else {
    for (int k = kts; k <= kte; ++k) sqc2[k] = I.sqc[k];
  }

  // sqv (BS:3447-3521)
  {
    const int k = kts;
    real qvflux = I.flqv;  // BS:3470
    if (qvflux < 0.0) qvflux = fmx(qvflux, (fmn(0.9 * I.sqv[kts] - 1e-8, 0.0) / dtz[kts]));
    a[k] = -dtz[k] * khdz[k] * rhoinv[k];
    b[k] = 1. + dtz[k] * (khdz[k + 1] + khdz[k]) * rhoinv[k] - 0.5 * dtz[k] * rhoinv[k] * s_aw[k + 1] -
           0.5 * dtz[k] * rhoinv[k] * sd_aw[k + 1];
    c[k] = -dtz[k] * khdz[k + 1] * rhoinv[k] - 0.5 * dtz[k] * rhoinv[k] * s_aw[k + 1] -
           0.5 * dtz[k] * rhoinv[k] * sd_aw[k + 1];
    d[k] = I.sqv[k] + dtz[k] * rhosfc * qvflux * rhoinv[k] + I.qcd[k] * delt - dtz[k] * rhoinv[k] * I.s_awqv[k + 1] -
           dtz[k] * rhoinv[k] * I.sd_awqv[k + 1] + I.sub_sqv[k] * delt + I.det_sqv[k] * delt;
  }
  for (int k = kts + 1; k <= kte - 1; ++k) {
    a[k] = -dtz[k] * khdz[k] * rhoinv[k] + 0.5 * dtz[k] * rhoinv[k] * s_aw[k] + 0.5 * dtz[k] * rhoinv[k] * sd_aw[k];
    b[k] = 1. + dtz[k] * (khdz[k] + khdz[k + 1]) * rhoinv[k] + 0.5 * dtz[k] * rhoinv[k] * (s_aw[k] - s_aw[k + 1]) +
           0.5 * dtz[k] * rhoinv[k] * (sd_aw[k] - sd_aw[k + 1]);
    c[k] = -dtz[k] * khdz[k + 1] * rhoinv[k] - 0.5 * dtz[k] * rhoinv[k] * s_aw[k + 1] -
           0.5 * dtz[k] * rhoinv[k] * sd_aw[k + 1];
    d[k] = I.sqv[k] + I.qcd[k] * delt + dtz[k] * rhoinv[k] * (I.s_awqv[k] - I.s_awqv[k + 1]) +
           dtz[k] * rhoinv[k] * (I.sd_awqv[k] - I.sd_awqv[k + 1]) + I.sub_sqv[k] * delt + I.det_sqv[k] * delt;
  }
  a[kte] = 0.;
  b[kte] = 1.;
  c[kte] = 0.;
  d[kte] = I.sqv[kte];
  tridiag2(kte, a, b, c, d, sqv2);

  // sqi (BS:3526-3570)
  if (I.FLAG_QI) {
    {
      const int k = kts;
      a[k] = -dtz[k] * khdz[k] * rhoinv[k];
      b[k] = 1. + dtz[k] * (khdz[k + 1] + khdz[k]) * rhoinv[k];
      c[k] = -dtz[k] * khdz[k + 1] * rhoinv[k];
      d[k] = I.sqi[k];
    }
    for (int k = kts + 1; k <= kte - 1; ++k) {
      a[k] = -dtz[k] * khdz[k] * rhoinv[k];
      b[k] = 1. + dtz[k] * (khdz[k] + khdz[k + 1]) * rhoinv[k];
      c[k] = -dtz[k] * khdz[k + 1] * rhoinv[k];
      d[k] = I.sqi[k];
    }
    a[kte] = 0.;
    b[kte] = 1.;
    c[kte] = 0.;
    d[kte] = I.sqi[kte];
    tridiag2(kte, a, b, c, d, sqi2);
  } else {
    for (int k = kts; k <= kte; ++k) sqi2[k] = I.sqi[k];
  }

  for (int k = kts; k <= kte; ++k) sqs2[k] = I.sqs[k];  // BS:3576-3607 (snow mixing disabled: .false.)

  // qni (BS:3612-3647)
  if (I.FLAG_QNI && I.bl_mynn_mixscalars) {
    {
      const int k = kts;
      a[k] = -dtz[k] * khdz[k] * rhoinv[k];
      b[k] = 1. + dtz[k] * (khdz[k + 1] + khdz[k]) * rhoinv[k] - 0.5 * dtz[k] * rhoinv[k] * s_aw[k + 1] * nonloc;
      c[k] = -dtz[k] * khdz[k + 1] * rhoinv[k] - 0.5 * dtz[k] * rhoinv[k] * s_aw[k + 1] * nonloc;
      d[k] = I.qni[k] - dtz[k] * rhoinv[k] * I.s_awqni[k + 1] * nonloc;
    }
    for (int k = kts + 1; k <= kte - 1; ++k) {
      a[k] = -dtz[k] * khdz[k] * rhoinv[k] + 0.5 * dtz[k] * rhoinv[k] * s_aw[k] * nonloc;
      b[k] = 1. + dtz[k] * (khdz[k] + khdz[k + 1]) * rhoinv[k] +
             0.5 * dtz[k] * rhoinv[k] * (s_aw[k] - s_aw[k + 1]) * nonloc;
      c[k] = -dtz[k] * khdz[k + 1] * rhoinv[k] - 0.5 * dtz[k] * rhoinv[k] * s_aw[k + 1] * nonloc;
      d[k] = I.qni[k] + dtz[k] * rhoinv[k] * (I.s_awqni[k] - I.s_awqni[k + 1]) * nonloc;
    }
    a[kte] = 0.;
    b[kte] = 1.;
    c[kte] = 0.;
    d[kte] = I.qni[kte];
    tridiag2(kte, a, b, c, d, x);
    for (int k = kts; k <= kte; ++k) qni2[k] = x[k];
  } else {
    for (int k = kts; k <= kte; ++k) qni2[k] = I.qni[k];
  }
  // qnc, qnwfa, qnifa, qnbca (BS:3653-3814) and ozone (BS:3819-3851): their flags are .false.
  for (int k = kts; k <= kte; ++k) O.dozone[k] = 0.0;

  for (int k = kts; k <= kte; ++k) O.dqv[k] = (sqv2[k] - I.sqv[k]) / delt;  // BS:3891
  // bl_mynn_cloudmix = .true. (BS:3896-3965)
  if (I.FLAG_QC) {
    for (int k = kts; k <= kte; ++k) O.dqc[k] = (sqc2[k] - I.sqc[k]) / delt;
  } else {
    for (int k = kts; k <= kte; ++k) O.dqc[k] = 0.;
  }
  for (int k = kts; k <= kte; ++k) O.dqnc[k] = 0.;  // FLAG_QNC = .false.
  if (I.FLAG_QI) {
    for (int k = kts; k <= kte; ++k) O.dqi[k] = (sqi2[k] - I.sqi[k]) / delt;
  } else {
    for (int k = kts; k <= kte; ++k) O.dqi[k] = 0.;
  }
  for (int k = kts; k <= kte; ++k) O.dqs[k] = 0.;  // BS:3943 (disabled)
  if (I.FLAG_QNI && I.bl_mynn_mixscalars) {
    for (int k = kts; k <= kte; ++k) O.dqni[k] = (qni2[k] - I.qni[k]) / delt;
  } else {
    for (int k = kts; k <= kte; ++k) O.dqni[k] = 0.;
  }
  moisture_check(C, kte, delt, delp, I.exner, sqv2, sqc2, sqi2, sqs2, thl, O.dqv, O.dqc, O.dqi, O.dqs, O.dth);
  for (int k = kts; k <= kte; ++k)  // BS:3985
    if (O.dozone[k] * delt + I.ozone[k] < 0.) O.dozone[k] = -I.ozone[k] * 0.99 / delt;
  if (I.FLAG_QI) {  // BS:3994
    for (int k = kts; k <= kte; ++k)
      O.dth[k] = (thl[k] + C.xlvcp / I.exner[k] * sqc2[k] + C.xlscp / I.exner[k] * (sqi2[k] + I.sqs[k]) - I.th[k]) /
                 delt;
  } else {
    for (int k = kts; k <= kte; ++k) O.dth[k] = (thl[k] + C.xlvcp / I.exner[k] * sqc2[k] - I.th[k]) / delt;
  }
  for (int k = kts; k <= kte; ++k) {  // BS:4018-4048 (flags .false.)
    O.dqnwfa[k] = 0.;
    O.dqnifa[k] = 0.;
    O.dqnbca[k] = 0.;
  }
}

// BS:4316-4340
CPPMPAS_HD inline void retrieve_exchange_coeffs(int kte, const real* dfm, const real* dfh, const real* dz, real* K_m,
                                                real* K_h) {
  const int kts = 1;
  K_m[kts] = 0.;
  K_h[kts] = 0.;
  for (int k = kts + 1; k <= kte; ++k) {
    const real dzk = 0.5 * (dz[k] + dz[k - 1]);
    K_m[k] = dfm[k] * dzk;
    K_h[k] = dfh[k] * dzk;
  }
}

// BS:4475-4615
CPPMPAS_HD inline void get_pblh(int kte, real& zi, const real* thetav1D, const real* qke1D, const real* zw1D,
                                const real* dz1D, real landsea, int& kzi) {
  const int kts = 1;
  const real sbl_lim = 200., sbl_damp = 400.;
  kzi = 2;  // BS:4513
  int k = kts + 1;
  int kthv = 1;
  real minthv = 9.E9;
  while (zw1D[k] <= 200.) {  // BS:4519
    // no bound on k in MPAS: a column top at or below 200 m reads thetav1D(k) past kte (undefined)
    if (k > kte) {
      kernel_fail("get_pblh: the model top is at or below 200 m (MPAS reads thetav1D past kte there)");
      return;
    }
    if (minthv > thetav1D[k]) {
      minthv = thetav1D[k];
      kthv = k;
    }
    k = k + 1;
  }
  (void)kthv;
  zi = 0.;
  real delt_thv;
  if ((landsea - 1.5) >= 0) {  // BS:4532
    delt_thv = 1.0;
  } else {
    delt_thv = 1.25;
  }
  zi = 0.;
  for (k = kts + 1; k <= kte - 1; ++k) {  // BS:4543
    if (thetav1D[k] >= (minthv + delt_thv)) {
      zi = zw1D[k] - dz1D[k - 1] *
                         fmn((thetav1D[k] - (minthv + delt_thv)) / fmx(thetav1D[k] - thetav1D[k - 1], 1E-6), 1.0);
    }
    if (k == kte - 1) zi = zw1D[kts + 1];  // EXIT SAFEGUARD
    if (zi != 0.0) break;
  }
  const real maxqke = fmx(qke1D[kts], 0.);  // BS:4560
  real TKEeps = maxqke / 40.;
  TKEeps = fmx(TKEeps, 0.02);
  real PBLH_TKE = 0.;
  for (k = kts + 1; k <= kte - 1; ++k) {  // BS:4569
    const real qtke = fmx(qke1D[k] / 2., 0.);
    const real qtkem1 = fmx(qke1D[k - 1] / 2., 0.);
    if (qtke <= TKEeps) {
      PBLH_TKE = zw1D[k] - dz1D[k - 1] * fmn((TKEeps - qtke) / fmx(qtkem1 - qtke, 1E-6), 1.0);
      PBLH_TKE = fmx(PBLH_TKE, zw1D[kts + 1]);
    }
    if (k == kte - 1) PBLH_TKE = zw1D[kts + 1];  // EXIT SAFEGUARD
    if (PBLH_TKE != 0.) break;
  }
  PBLH_TKE = fmn(PBLH_TKE, zi + 350.);  // BS:4591
  PBLH_TKE = fmx(PBLH_TKE, fmx(zi - 350., 10.));
  const real wt = .5 * libm::tanh((zi - sbl_lim) / sbl_damp) + .5;  // BS:4594
  if (maxqke <= 0.05) {
  } else {
    zi = PBLH_TKE * (1. - wt) + zi * wt;
  }
  for (k = kts + 1; k <= kte - 1; ++k) {  // BS:4603
    if (zw1D[k] >= zi) {
      kzi = k - 1;
      break;
    }
  }
}

// BS:5789-5846 (an out-of-line function in the object, condensation_edmf.isra.0)
CPPMPAS_HD inline void condensation_edmf(const Const& C, real QT, real THL, real P, real zagl, real& THV, real& QC) {
  const int niter = 50;  // BS:5808
  const real diff = 1.e-6;
  const real EXN = libm::pow(P / C.p1000mb, C.rcp);  // BS:5812
  for (int i = 1; i <= niter; ++i) {
    const real T = EXN * THL + C.xlvcp * QC;
    const real QS = qsat_blend(T, C.t0c, tice, P);
    const real QCOLD = QC;
    QC = 0.5 * QC + 0.5 * fmx((QT - QS), 0.);
    if (fabs_(QC - QCOLD) < diff) break;
  }
  const real T = EXN * THL + C.xlvcp * QC;  // BS:5822
  const real QS = qsat_blend(T, C.t0c, tice, P);
  QC = fmx(QT - QS, 0.);
  if (zagl < 100.) QC = 0.;  // BS:5827
  THV = (THL + C.xlvcp * QC) * (1. + QT * (C.rvovrd - 1.) - C.rvovrd * QC);  // BS:5830
}

// BS:4637-5785 DMP_mf with momentum_opt = tke_opt = .false. (bl_mynn_edmf_mom / _tke), scalar_opt =
// bl_mynn_mixscalars, mix_chem = .false., env_subs = .false. (BS:407)
struct MfIn {
  const real *zw, *dz, *p, *rho, *u, *v, *w, *th, *thl, *thv, *tk, *qt, *qv, *qc, *qke, *qnc, *qni, *qnwfa, *qnifa,
      *qnbca, *exner;
  real dt, ust, flt, fltv, flq, flqv, pblh, dx, landsea, ts, Psig_shcu;
  int kpbl;
  bool scalar_opt;
  const real* rstoch_col;
};
struct MfOut {
  real *vt, *vq, *sgm;                                                     // inout
  real *edmf_a, *edmf_w, *edmf_qt, *edmf_thl, *edmf_ent, *edmf_qc;          // out
  real *s_aw, *s_awthl, *s_awqt, *s_awqv, *s_awqc, *s_awu, *s_awv, *s_awqke, *s_awqnc, *s_awqni, *s_awqnwfa,
      *s_awqnifa, *s_awqnbca;                                              // kts:kte+1
  real *sub_thl, *sub_sqv, *sub_u, *sub_v, *det_thl, *det_sqv, *det_sqc, *det_u, *det_v;
  real *qc_bl1d, *cldfra_bl1d;                                             // inout
  real *maxwidth, *maxmf, *ztop;
  int* ktop;
};
CPPMPAS_HD inline void dmp_mf(const Const& C, int kte, const MfIn& I, const MfOut& O) {
  const int kts = 1;
  const real* zw = I.zw;
  const real* dz = I.dz;
  const real pblh = I.pblh, dx = I.dx, landsea = I.landsea, flt = I.flt, fltv = I.fltv, flq = I.flq;
  // BS:4724-4730 (kts:kte+1, 1:NUP); 1-based [k][i]
  real UPW[KM][NUP + 1], UPTHL[KM][NUP + 1], UPQT[KM][NUP + 1], UPQC[KM][NUP + 1], UPA[KM][NUP + 1],
      UPU[KM][NUP + 1], UPV[KM][NUP + 1], UPTHV[KM][NUP + 1], UPQKE[KM][NUP + 1], UPQNC[KM][NUP + 1],
      UPQNI[KM][NUP + 1], UPQNWFA[KM][NUP + 1], UPQNIFA[KM][NUP + 1], UPQNBCA[KM][NUP + 1], ENT[KM][NUP + 1];
  real exneri[KM], dzi[KM], rhoz[KM], edmf_th[KM], envm_thl[KM], envm_sqv[KM], envm_sqc[KM], envm_u[KM],
      envm_v[KM];
  const real Atot = 0.10, lmax = 1000., lmin = 300., dlmin = 0., dcut = 1.2;  // BS:4752-4758
  const real cf_thresh = 0.5, fluxportion = 0.75, dzpmax = 300., pgfac = 0.00;
  for (int k = kts; k <= kte + 1; ++k)  // BS:4836-4855
    for (int i = 1; i <= NUP; ++i) {
      UPW[k][i] = 0.;
      UPTHL[k][i] = 0.;
      UPTHV[k][i] = 0.;
      UPQT[k][i] = 0.;
      UPA[k][i] = 0.;
      UPU[k][i] = 0.;
      UPV[k][i] = 0.;
      UPQC[k][i] = 0.;
      UPQKE[k][i] = 0.;
      UPQNC[k][i] = 0.;
      UPQNI[k][i] = 0.;
      UPQNWFA[k][i] = 0.;
      UPQNIFA[k][i] = 0.;
      UPQNBCA[k][i] = 0.;
      if (k <= kte) ENT[k][i] = 0.001;
    }
  for (int k = kts; k <= kte; ++k) {  // BS:4857-4862
    O.edmf_a[k] = 0.;
    O.edmf_w[k] = 0.;
    O.edmf_qt[k] = 0.;
    O.edmf_thl[k] = 0.;
    O.edmf_ent[k] = 0.;
    O.edmf_qc[k] = 0.;
  }
  for (int k = kts; k <= kte + 1; ++k) {  // BS:4868-4880
    O.s_aw[k] = 0.;
    O.s_awthl[k] = 0.;
    O.s_awqt[k] = 0.;
    O.s_awqv[k] = 0.;
    O.s_awqc[k] = 0.;
    O.s_awu[k] = 0.;
    O.s_awv[k] = 0.;
    O.s_awqke[k] = 0.;
    O.s_awqnc[k] = 0.;
    O.s_awqni[k] = 0.;
    O.s_awqnwfa[k] = 0.;
    O.s_awqnifa[k] = 0.;
    O.s_awqnbca[k] = 0.;
  }
  for (int k = kts; k <= kte; ++k) {  // BS:4886-4894
    O.sub_thl[k] = 0.;
    O.sub_sqv[k] = 0.;
    O.sub_u[k] = 0.;
    O.sub_v[k] = 0.;
    O.det_thl[k] = 0.;
    O.det_sqv[k] = 0.;
    O.det_sqc[k] = 0.;
    O.det_u[k] = 0.;
    O.det_v[k] = 0.;
  }
  real nup2 = NUP;  // BS:4895
  real maxw = 0.0;  // BS:4899
  real cloud_base = 9000.0;
  int k50 = 0;
  for (int k = 1; k <= kte - 1; ++k) {  // BS:4901
    if (zw[k] > pblh + 500.) break;
    real wpbl = I.w[k];
    if (I.w[k] < 0.) wpbl = 2. * I.w[k];
    maxw = fmx(maxw, fabs_(wpbl));
    if (zw[k] <= 50.) k50 = k;
    const real qc_sgs = fmx(I.qc[k], O.qc_bl1d[k]);
    if (qc_sgs > 1E-5 && (O.cldfra_bl1d[k] >= 0.5) && cloud_base == 9000.0) cloud_base = 0.5 * (zw[k] + zw[k + 1]);
  }
  maxw = fmx(0., maxw - 1.0);  // BS:4919
  real Psig_w = fmx(0.0, 1.0 - maxw);
  Psig_w = fmn(Psig_w, I.Psig_shcu);
  real fltv2 = fltv;  // BS:4924
  if (Psig_w == 0.0 && fltv > 0.0) fltv2 = -1. * fltv;
  bool superadiabatic = false;  // BS:4929
  real hux;
  if ((landsea - 1.5) >= 0) {
    hux = -0.001;
  } else {
    hux = -0.005;
  }
  const real tvs = I.ts * (1.0 + C.p608 * I.qv[kts]);  // BS:4935
  for (int k = 1; k <= (k50 - 1 > 1 ? k50 - 1 : 1); ++k) {
    if (k == 1) {
      if ((I.thv[k] - tvs) / (0.5 * dz[k]) < hux) {
        superadiabatic = true;
      } else {
        superadiabatic = false;
        break;
      }
    } else {
      if ((I.thv[k] - I.thv[k - 1]) / (0.5 * (dz[k] + dz[k - 1])) < hux) {
        superadiabatic = true;
      } else {
        superadiabatic = false;
        break;
      }
    }
  }
  real maxwidth = fmn(dx * dcut, lmax);  // BS:4963
  maxwidth = fmn(maxwidth, 1.1 * pblh);
  if ((landsea - 1.5) < 0) {  // BS:4967
    maxwidth = fmn(maxwidth, 0.5 * cloud_base);
  } else {
    maxwidth = fmn(maxwidth, 0.9 * cloud_base);
  }
  const real wspd_pbl = __builtin_sqrt(fmx(I.u[kts] * I.u[kts] + I.v[kts] * I.v[kts], 0.01));  // BS:4973
  (void)wspd_pbl;
  real width_flx;
  if ((landsea - 1.5) < 0) {  // BS:4976
    width_flx = fmx(fmn(1000. * (0.6 * libm::tanh((fltv - 0.040) / 0.04) + .5), 1000.), 0.);
  } else {
    width_flx = fmx(fmn(1000. * (0.6 * libm::tanh((fltv - 0.007) / 0.02) + .5), 1000.), 0.);
  }
  maxwidth = fmn(maxwidth, width_flx);  // BS:4981
  real minwidth = lmin;
  if (maxwidth >= (lmax - 1.0) && fltv > 0.2) minwidth = lmin + dlmin * fmn((fltv - 0.2) / 0.3, 1.);
  if (maxwidth <= minwidth) {  // BS:4987
    nup2 = 0;
    maxwidth = 0.0;
  }
  int ktop = 0;  // BS:4993
  real ztop = 0.0;
  real maxmf = 0.0;
  real dl = 0.;
  const real* thl = I.thl;
  if (fltv2 > 0.002 && (maxwidth > minwidth) && superadiabatic) {  // BS:4998
    real cn = 0.;
    const real d = -1.9;
    dl = (maxwidth - minwidth) / real(NUP - 1);
    for (int i = 1; i <= NUP; ++i) {  // BS:5004 (scalar pow)
      const real l = minwidth + dl * real(i - 1);
      cn = cn + libm::pow(l, d) * (l * l) / (dx * dx) * dl;
    }
    const real Cc = Atot / cn;  // BS:5009
    real acfac;
    if ((landsea - 1.5) < 0) {  // BS:5012
      acfac = .5 * libm::tanh((fltv2 - 0.02) / 0.05) + .5;
    } else {
      acfac = .5 * libm::tanh((fltv2 - 0.01) / 0.03) + .5;
    }
    real ac_wsp;
    if (wspd_pbl <= 10.) {  // BS:5021
      ac_wsp = 1.0;
    } else {
      ac_wsp = 1.0 - fmn((wspd_pbl - 10.0) / 15., 1.0);
    }
    acfac = acfac * ac_wsp;
    for (int i = 1; i <= NUP; ++i) {  // BS:5030: _ZGVbN2vv_pow, four pairs (i = 1,2 / 3,4 / 5,6 / 7,8)
      const real l = minwidth + dl * real(i - 1);
      const real N = Cc * libm::mvec_pow(l, d);
      UPA[1][i] = N * l * l / (dx * dx) * dl;
      UPA[1][i] = UPA[1][i] * acfac;
    }
    const real z0 = 50.;  // BS:5042
    const real pwmin = 0.1;
    const real pwmax = 0.4;
    const real wstar = fmx(1.E-2, libm::pow(C.gtr * fltv2 * pblh, onethird));  // BS:5046
    const real qstar = fmx(flq, 1.0E-5) / wstar;
    const real thstar = flt / wstar;
    const real csigma = 1.34;  // BS:5050-5054 (water and land)
    real exc_fac;
    if ((landsea - 1.5) >= 0) {  // BS:5059 (env_subs = .false.)
      exc_fac = 0.58 * 4.0;
    } else {
      exc_fac = 0.58;
    }
    exc_fac = exc_fac * ac_wsp;
    const real zp = libm::pow(z0 / pblh, onethird);  // BS:5071-5073: one pow
    const real sigmaW = csigma * wstar * zp * (1 - 0.8 * z0 / pblh);
    const real sigmaQT = csigma * qstar * zp;
    const real sigmaTH = csigma * thstar * zp;
    const real wmin = fmn(sigmaW * pwmin, 0.1);  // BS:5077
    const real wmax = fmn(sigmaW * pwmax, 0.5);
    for (int i = 1; i <= NUP; ++i) {  // BS:5081
      UPW[1][i] = wmin + real(i) / real(NUP) * (wmax - wmin);
      UPU[1][i] = (I.u[kts] * dz[kts + 1] + I.u[kts + 1] * dz[kts]) / (dz[kts] + dz[kts + 1]);
      UPV[1][i] = (I.v[kts] * dz[kts + 1] + I.v[kts + 1] * dz[kts]) / (dz[kts] + dz[kts + 1]);
      UPQC[1][i] = 0.0;
      const real exc_heat = exc_fac * UPW[1][i] * sigmaTH / sigmaW;
      UPTHV[1][i] = (I.thv[kts] * dz[kts + 1] + I.thv[kts + 1] * dz[kts]) / (dz[kts] + dz[kts + 1]) + exc_heat;
      UPTHL[1][i] = (thl[kts] * dz[kts + 1] + thl[kts + 1] * dz[kts]) / (dz[kts] + dz[kts + 1]) + exc_heat;
      const real exc_moist = exc_fac * UPW[1][i] * sigmaQT / sigmaW;
      UPQT[1][i] = (I.qt[kts] * dz[kts + 1] + I.qt[kts + 1] * dz[kts]) / (dz[kts] + dz[kts + 1]) + exc_moist;
      UPQKE[1][i] = (I.qke[kts] * dz[kts + 1] + I.qke[kts + 1] * dz[kts]) / (dz[kts] + dz[kts + 1]);
      UPQNC[1][i] = (I.qnc[kts] * dz[kts + 1] + I.qnc[kts + 1] * dz[kts]) / (dz[kts] + dz[kts + 1]);
      UPQNI[1][i] = (I.qni[kts] * dz[kts + 1] + I.qni[kts + 1] * dz[kts]) / (dz[kts] + dz[kts + 1]);
      UPQNWFA[1][i] = (I.qnwfa[kts] * dz[kts + 1] + I.qnwfa[kts + 1] * dz[kts]) / (dz[kts] + dz[kts + 1]);
      UPQNIFA[1][i] = (I.qnifa[kts] * dz[kts + 1] + I.qnifa[kts + 1] * dz[kts]) / (dz[kts] + dz[kts + 1]);
      UPQNBCA[1][i] = (I.qnbca[kts] * dz[kts + 1] + I.qnbca[kts + 1] * dz[kts]) / (dz[kts] + dz[kts + 1]);
    }
    for (int k = kts; k <= kte; ++k) {  // BS:5119
      envm_thl[k] = thl[k];
      envm_sqv[k] = I.qv[k];
      envm_sqc[k] = I.qc[k];
      envm_u[k] = I.u[k];
      envm_v[k] = I.v[k];
    }
    for (int k = kts; k <= kte - 1; ++k)  // BS:5124
      rhoz[k] = (I.rho[k] * dz[k + 1] + I.rho[k + 1] * dz[k]) / (dz[k + 1] + dz[k]);
    rhoz[kte] = I.rho[kte];
    const real dxsa = 1. - fmn(fmx((12000.0 - dx) / (12000.0 - 3000.0), 0.), 1.);  // BS:5130
    real dzp = 0.;
    for (int i = 1; i <= NUP; ++i) {  // BS:5133
      real QCn = 0.;
      int overshoot = 0;
      const real l = minwidth + dl * real(i - 1);
      for (int k = kts + 1; k <= kte - 1; ++k) {  // BS:5137
        const real wmin_k = 0.3 + l * 0.0005;  // BS:5140
        ENT[k][i] = 0.33 / (fmn(fmx(UPW[k - 1][i], wmin_k), 0.9) * l);
        ENT[k][i] = fmx(ENT[k][i], 0.0003);
        if (zw[k] >= fmn(pblh + 1500., 4000.)) ENT[k][i] = ENT[k][i] + (zw[k] - fmn(pblh + 1500., 4000.)) * 5.0E-6;
        ENT[k][i] = ENT[k][i] * (1.0 - I.rstoch_col[k]);  // BS:5158
        ENT[k][i] = fmn(ENT[k][i], 0.9 / (zw[k + 1] - zw[k]));
        const real Uk = (I.u[k] * dz[k + 1] + I.u[k + 1] * dz[k]) / (dz[k + 1] + dz[k]);  // BS:5163
        const real Ukm1 = (I.u[k - 1] * dz[k] + I.u[k] * dz[k - 1]) / (dz[k - 1] + dz[k]);
        const real Vk = (I.v[k] * dz[k + 1] + I.v[k + 1] * dz[k]) / (dz[k + 1] + dz[k]);
        const real Vkm1 = (I.v[k - 1] * dz[k] + I.v[k] * dz[k - 1]) / (dz[k - 1] + dz[k]);
        const real EntExp = ENT[k][i] * (zw[k + 1] - zw[k]);  // BS:5169
        const real EntExm = EntExp * 0.3333;
        const real QTn = UPQT[k - 1][i] * (1. - EntExp) + I.qt[k] * EntExp;
        const real THLn = UPTHL[k - 1][i] * (1. - EntExp) + thl[k] * EntExp;
        const real Un = UPU[k - 1][i] * (1. - EntExm) + I.u[k] * EntExm + dxsa * pgfac * (Uk - Ukm1);
        const real Vn = UPV[k - 1][i] * (1. - EntExm) + I.v[k] * EntExm + dxsa * pgfac * (Vk - Vkm1);
        const real QKEn = UPQKE[k - 1][i] * (1. - EntExp) + I.qke[k] * EntExp;
        const real QNCn = UPQNC[k - 1][i] * (1. - EntExp) + I.qnc[k] * EntExp;
        const real QNIn = UPQNI[k - 1][i] * (1. - EntExp) + I.qni[k] * EntExp;
        const real QNWFAn = UPQNWFA[k - 1][i] * (1. - EntExp) + I.qnwfa[k] * EntExp;
        const real QNIFAn = UPQNIFA[k - 1][i] * (1. - EntExp) + I.qnifa[k] * EntExp;
        const real QNBCAn = UPQNBCA[k - 1][i] * (1. - EntExp) + I.qnbca[k] * EntExp;
        const real qc_ent = QCn;  // BS:5184
        const real qt_ent = QTn;
        const real thl_ent = THLn;
        const real Pk = (I.p[k] * dz[k + 1] + I.p[k + 1] * dz[k]) / (dz[k + 1] + dz[k]);  // BS:5206
        real THVn;
        condensation_edmf(C, QTn, THLn, Pk, zw[k + 1], THVn, QCn);
        const real THVk = (I.thv[k] * dz[k + 1] + I.thv[k + 1] * dz[k]) / (dz[k + 1] + dz[k]);  // BS:5211
        const real THVkm1 = (I.thv[k - 1] * dz[k] + I.thv[k] * dz[k - 1]) / (dz[k - 1] + dz[k]);
        const real B = C.grav * (THVn / THVk - 1.0);  // BS:5215
        real BCOEFF;
        if (B > 0.) {
          BCOEFF = 0.15;
        } else {
          BCOEFF = 0.2;
        }
        real Wn;
        if (UPW[k - 1][i] < 0.2) {  // BS:5230
          Wn = UPW[k - 1][i] + (-2. * ENT[k][i] * UPW[k - 1][i] + BCOEFF * B / fmx(UPW[k - 1][i], 0.2)) *
                                   fmn(zw[k] - zw[k - 1], 250.);
        } else {
          Wn = UPW[k - 1][i] + (-2. * ENT[k][i] * UPW[k - 1][i] + BCOEFF * B / UPW[k - 1][i]) *
                                   fmn(zw[k] - zw[k - 1], 250.);
        }
        if (Wn > UPW[k - 1][i] + fmn(1.25 * (zw[k] - zw[k - 1]) / 200., 2.0))  // BS:5237
          Wn = UPW[k - 1][i] + fmn(1.25 * (zw[k] - zw[k - 1]) / 200., 2.0);
        if (Wn < UPW[k - 1][i] - fmn(1.25 * (zw[k] - zw[k - 1]) / 200., 2.0))
          Wn = UPW[k - 1][i] - fmn(1.25 * (zw[k] - zw[k - 1]) / 200., 2.0);
        Wn = fmn(fmx(Wn, 0.0), 3.0);
        if (k == kts + 1 && Wn == 0.) {  // BS:5248
          nup2 = 0;
          break;
        }
        if (Wn <= 0.0 && overshoot == 0) {  // BS:5264
          overshoot = 1;
          if (THVk - THVkm1 > 0.0) {
            const real bvf = __builtin_sqrt(C.gtr * (THVk - THVkm1) / dz[k]);
            const real Frz = UPW[k - 1][i] / (bvf * dz[k]);
            dzp = dz[k] * fmx(fmn(Frz, 1.0), 0.0);
          }
        } else {
          dzp = dz[k];
        }
        const real aratio = fmn(UPA[k - 1][i] / (1. - UPA[k - 1][i]), 0.5);  // BS:5285
        const real detturb = 0.00008;
        const real oow = -0.060 / fmx(1.0, (0.5 * (Wn + UPW[k - 1][i])));
        const real detrate = fmn(fmx(oow * (Wn - UPW[k - 1][i]) / dz[k], detturb), .0002);
        const real detrateUV = fmn(fmx(oow * (Wn - UPW[k - 1][i]) / dz[k], detturb), .0001);
        envm_thl[k] = envm_thl[k] + (0.5 * (thl_ent + UPTHL[k - 1][i]) - thl[k]) * detrate * aratio * fmn(dzp, dzpmax);
        const real qv_ent = 0.5 * (fmx(qt_ent - qc_ent, 0.) + fmx(UPQT[k - 1][i] - UPQC[k - 1][i], 0.));
        envm_sqv[k] = envm_sqv[k] + (qv_ent - I.qv[k]) * detrate * aratio * fmn(dzp, dzpmax);
        if (UPQC[k - 1][i] > 1E-8) {  // BS:5293
          real qc_grid;
          if (I.qc[k] > 1E-6) {
            qc_grid = I.qc[k];
          } else {
            qc_grid = O.cldfra_bl1d[k] * O.qc_bl1d[k];
          }
          envm_sqc[k] = envm_sqc[k] + fmx(UPA[k - 1][i] * 0.5 * (QCn + UPQC[k - 1][i]) - qc_grid, 0.0) * detrate *
                                          aratio * fmn(dzp, dzpmax);
        }
        envm_u[k] = envm_u[k] + (0.5 * (Un + UPU[k - 1][i]) - I.u[k]) * detrateUV * aratio * fmn(dzp, dzpmax);
        envm_v[k] = envm_v[k] + (0.5 * (Vn + UPV[k - 1][i]) - I.v[k]) * detrateUV * aratio * fmn(dzp, dzpmax);
        if (Wn > 0.) {  // BS:5304
          UPW[k][i] = Wn;
          UPTHV[k][i] = THVn;
          UPTHL[k][i] = THLn;
          UPQT[k][i] = QTn;
          UPQC[k][i] = QCn;
          UPU[k][i] = Un;
          UPV[k][i] = Vn;
          UPQKE[k][i] = QKEn;
          UPQNC[k][i] = QNCn;
          UPQNI[k][i] = QNIn;
          UPQNWFA[k][i] = QNWFAn;
          UPQNIFA[k][i] = QNIFAn;
          UPQNBCA[k][i] = QNBCAn;
          UPA[k][i] = UPA[k - 1][i];
          ktop = ktop > k ? ktop : k;
        } else {
          break;
        }
      }
    }
  } else {
    nup2 = 0.;  // BS:5352
  }
  ktop = ktop < kte - 1 ? ktop : kte - 1;  // BS:5355
  if (ktop == 0) {
    ztop = 0.0;
  } else {
    ztop = zw[ktop];
  }
  if (nup2 > 0) {  // BS:5362
    for (int i = 1; i <= NUP; ++i)
      for (int k = kts; k <= kte - 1; ++k) {
        O.s_aw[k + 1] = O.s_aw[k + 1] + rhoz[k] * UPA[k][i] * UPW[k][i] * Psig_w;
        O.s_awthl[k + 1] = O.s_awthl[k + 1] + rhoz[k] * UPA[k][i] * UPW[k][i] * UPTHL[k][i] * Psig_w;
        O.s_awqt[k + 1] = O.s_awqt[k + 1] + rhoz[k] * UPA[k][i] * UPW[k][i] * UPQT[k][i] * Psig_w;
        const real qc_plume = UPQC[k][i];
        O.s_awqc[k + 1] = O.s_awqc[k + 1] + rhoz[k] * UPA[k][i] * UPW[k][i] * qc_plume * Psig_w;
        O.s_awqv[k + 1] = O.s_awqt[k + 1] - O.s_awqc[k + 1];
      }
    if (I.scalar_opt) {  // BS:5410 (k up to kte: rhoz(kte) = rho(kte))
      for (int k = kts; k <= kte; ++k)
        for (int i = 1; i <= NUP; ++i) {
          O.s_awqnc[k + 1] = O.s_awqnc[k + 1] + rhoz[k] * UPA[k][i] * UPW[k][i] * UPQNC[k][i] * Psig_w;
          O.s_awqni[k + 1] = O.s_awqni[k + 1] + rhoz[k] * UPA[k][i] * UPW[k][i] * UPQNI[k][i] * Psig_w;
          O.s_awqnwfa[k + 1] = O.s_awqnwfa[k + 1] + rhoz[k] * UPA[k][i] * UPW[k][i] * UPQNWFA[k][i] * Psig_w;
          O.s_awqnifa[k + 1] = O.s_awqnifa[k + 1] + rhoz[k] * UPA[k][i] * UPW[k][i] * UPQNIFA[k][i] * Psig_w;
          O.s_awqnbca[k + 1] = O.s_awqnbca[k + 1] + rhoz[k] * UPA[k][i] * UPW[k][i] * UPQNBCA[k][i] * Psig_w;
        }
    }
    real flx1;
    if (O.s_aw[kts + 1] != 0.) {  // BS:5425
      dzi[kts] = 0.5 * (dz[kts] + dz[kts + 1]);
      flx1 = fmx(O.s_aw[kts + 1] * (I.th[kts] - I.th[kts + 1]) / dzi[kts], 1.0e-5);
    } else {
      flx1 = 0.0;
    }
    real adjustment = 1.0;  // BS:5433
    if (flx1 > fluxportion * flt / dz[kts] && flx1 > 0.0) {
      adjustment = fluxportion * flt / dz[kts] / flx1;
      for (int k = kts; k <= kte + 1; ++k) {
        O.s_aw[k] = O.s_aw[k] * adjustment;
        O.s_awthl[k] = O.s_awthl[k] * adjustment;
        O.s_awqt[k] = O.s_awqt[k] * adjustment;
        O.s_awqc[k] = O.s_awqc[k] * adjustment;
        O.s_awqv[k] = O.s_awqv[k] * adjustment;
        O.s_awqnc[k] = O.s_awqnc[k] * adjustment;
        O.s_awqni[k] = O.s_awqni[k] * adjustment;
        O.s_awqnwfa[k] = O.s_awqnwfa[k] * adjustment;
        O.s_awqnifa[k] = O.s_awqnifa[k] * adjustment;
        O.s_awqnbca[k] = O.s_awqnbca[k] * adjustment;
        for (int i = 1; i <= NUP; ++i) UPA[k][i] = UPA[k][i] * adjustment;
      }
    }
    for (int k = kts; k <= kte - 1; ++k)  // BS:5464
      for (int i = 1; i <= NUP; ++i) {
        O.edmf_a[k] = O.edmf_a[k] + UPA[k][i];
        O.edmf_w[k] = O.edmf_w[k] + rhoz[k] * UPA[k][i] * UPW[k][i];
        O.edmf_qt[k] = O.edmf_qt[k] + rhoz[k] * UPA[k][i] * UPQT[k][i];
        O.edmf_thl[k] = O.edmf_thl[k] + rhoz[k] * UPA[k][i] * UPTHL[k][i];
        O.edmf_ent[k] = O.edmf_ent[k] + rhoz[k] * UPA[k][i] * ENT[k][i];
        O.edmf_qc[k] = O.edmf_qc[k] + rhoz[k] * UPA[k][i] * UPQC[k][i];
      }
    for (int k = kts; k <= kte - 1; ++k) {  // BS:5474
      if (O.edmf_a[k] > 0.) {
        O.edmf_w[k] = O.edmf_w[k] / O.edmf_a[k];
        O.edmf_qt[k] = O.edmf_qt[k] / O.edmf_a[k];
        O.edmf_thl[k] = O.edmf_thl[k] / O.edmf_a[k];
        O.edmf_ent[k] = O.edmf_ent[k] / O.edmf_a[k];
        O.edmf_qc[k] = O.edmf_qc[k] / O.edmf_a[k];
        O.edmf_a[k] = O.edmf_a[k] * Psig_w;
        if (O.edmf_a[k] * O.edmf_w[k] > maxmf) maxmf = O.edmf_a[k] * O.edmf_w[k];
      }
    }
    for (int k = kts; k <= kte - 1; ++k) {  // BS:5584
      exneri[k] = (I.exner[k] * dz[k + 1] + I.exner[k + 1] * dz[k]) / (dz[k + 1] + dz[k]);
      edmf_th[k] = O.edmf_thl[k] + C.xlvcp / exneri[k] * O.edmf_qc[k];
      dzi[k] = 0.5 * (dz[k] + dz[k + 1]);
    }
    for (int k = kts + 1; k <= kte - 2; ++k) {  // BS:5593
      if (k > ktop) break;
      if (0.5 * (O.edmf_qc[k] + O.edmf_qc[k - 1]) > 0.0 && (O.cldfra_bl1d[k] < cf_thresh)) {
        const real Aup = (O.edmf_a[k] * dzi[k - 1] + O.edmf_a[k - 1] * dzi[k]) / (dzi[k - 1] + dzi[k]);
        const real QTp = (O.edmf_qt[k] * dzi[k - 1] + O.edmf_qt[k - 1] * dzi[k]) / (dzi[k - 1] + dzi[k]);
        // BS:5598 THp and BS:5603-5605 esat_blend / qsl are not used
        real QCp;
        if (O.edmf_qc[k] > 0.0 && O.edmf_qc[k - 1] > 0.0) {  // BS:5608
          QCp = (O.edmf_qc[k] * dzi[k - 1] + O.edmf_qc[k - 1] * dzi[k]) / (dzi[k - 1] + dzi[k]);
        } else {
          QCp = fmx(O.edmf_qc[k], O.edmf_qc[k - 1]);
        }
        const real tk = I.tk[k];
        const real xl = xl_blend(tk, C.t0c, tice, C.cice, C.cliq, C.cpv, C.xls, C.xlv);  // BS:5615
        const real qsat_tk = qsat_blend(tk, C.t0c, tice, I.p[k]);
        const real rsl = xl * qsat_tk / (C.r_v * (tk * tk));
        const real cpm = C.cp + I.qt[k] * C.cpv;
        const real a = 1. / (1. + xl * rsl / cpm);
        const real b9 = a * rsl;
        const real q2p = C.xlvcp / I.exner[k];  // BS:5624
        const real pt = thl[k] + q2p * QCp * Aup;
        const real bb = b9 * tk / pt;
        const real qww = 1. + 0.61 * I.qt[k];  // BS:5632
        const real alpha = 0.61 * pt;
        const real beta = pt * xl / (tk * C.cp) - 1.61 * pt;
        real sigq = 10. * Aup * (QTp - I.qt[k]);  // BS:5649
        sigq = fmx(sigq, qsat_tk * 0.02);
        sigq = fmn(sigq, qsat_tk * 0.25);
        const real qmq = a * (I.qt[k] - qsat_tk);  // BS:5654
        real Q1 = qmq / sigq;
        real mf_cf;
        if ((landsea - 1.5) >= 0) {  // BS:5657
          mf_cf = fmn(fmx(0.5 + 0.36 * libm::atan(1.55 * Q1), 0.01), 0.6);
          mf_cf = fmx(mf_cf, 1.2 * Aup);
          mf_cf = fmn(mf_cf, 5.0 * Aup);
        } else {
          mf_cf = fmn(fmx(0.5 + 0.36 * libm::atan(1.55 * Q1), 0.01), 0.6);
          mf_cf = fmx(mf_cf, 1.8 * Aup);
          mf_cf = fmn(mf_cf, 5.0 * Aup);
        }
        if (QCp * Aup > 5e-5) {  // BS:5684-5700 (water and land alike)
          O.qc_bl1d[k] = 1.86 * (QCp * Aup) - 2.2e-5;
        } else {
          O.qc_bl1d[k] = 1.18 * (QCp * Aup);
        }
        O.cldfra_bl1d[k] = mf_cf;
        Q1 = fmx(Q1, -2.25);  // BS:5708
        real Fng;
        if (Q1 >= 1.0) {  // BS:5713
          Fng = 1.0;
        } else if (Q1 >= -1.7 && Q1 < 1.0) {
          Fng = libm::exp(-0.4 * (Q1 - 1.0));
        } else if (Q1 >= -2.5 && Q1 < -1.7) {
          Fng = 3.0 + libm::exp(-3.8 * (Q1 + 1.7));
        } else {
          Fng = fmn(23.9 + libm::exp(-1.6 * (Q1 + 2.5)), 60.);
        }
        O.vt[k] = qww - (1.5 * Aup) * beta * bb * Fng - 1.;  // BS:5724
        O.vq[k] = alpha + (1.5 * Aup) * beta * a * Fng - C.tv0;
      }
    }
  }
  (void)UPTHV;    // BS:5307: not read (momentum / scalar options)
  (void)edmf_th;  // BS:5586: read by no statement in the options ported
  if (ktop > 0) {  // BS:5732
    real maxqc = O.edmf_qc[1];
    for (int k = 2; k <= ktop; ++k) maxqc = fmx(maxqc, O.edmf_qc[k]);
    if (maxqc < 1.E-8) maxmf = -1.0 * maxmf;
  }
  if (O.edmf_w[1] > 4.0) {  // BS:5740 (debugging print; edmf_w(1) <= rhoz(1)*wmax <= 0.5*rho)
    printf(" flq: %g fltv: %g\n pblh: %g\n", flq, fltv2, pblh);
  }
  *O.maxwidth = maxwidth;
  *O.maxmf = maxmf;
  *O.ztop = ztop;
  *O.ktop = ktop;
}

// BS:6273-6344
CPPMPAS_HD inline void scale_aware(real dx, real PBL1, real& Psig_bl, real& Psig_shcu) {
  Psig_bl = 1.0;
  Psig_shcu = 1.0;
  real dxdh = fmx(2.5 * dx, 10.) / fmn(PBL1, 3000.);  // BS:6290
  real pw = libm::pow(dxdh, 0.667);  // one pow for the two dxdh**0.667 of BS:6298
  Psig_bl = ((dxdh * dxdh) + 0.106 * pw) / ((dxdh * dxdh) + 0.066 * pw + 0.071);
  dxdh = fmx(2.5 * dx, 10.) / fmn(PBL1 + 500., 3500.);  // BS:6301
  pw = libm::pow(dxdh, 0.667);
  Psig_shcu = ((dxdh * dxdh) + 0.145 * pw) / ((dxdh * dxdh) + 0.172 * pw + 0.170);  // BS:6328
  if (Psig_bl > 1.0) Psig_bl = 1.0;
  if (Psig_bl < 0.0) Psig_bl = 0.0;
  if (Psig_shcu > 1.0) Psig_shcu = 1.0;
  if (Psig_shcu < 0.0) Psig_shcu = 0.0;
}

// BS:6346-6395
CPPMPAS_HD inline real phim(real zet) {
  const real am_st = 6.1, bm_st = 2.5, rbm_st = 1. / bm_st, am_unst = 10.;
  real phi_m;
  if (zet >= 0.0) {
    const real dummy_0 = 1 + libm::pow(zet, bm_st);
    const real dummy_1 = zet + libm::pow(dummy_0, rbm_st);
    const real dummy_11 = 1 + libm::pow(dummy_0, rbm_st - 1) * libm::pow(zet, bm_st - 1);
    const real dummy_2 = (-am_st / dummy_1) * dummy_11;
    phi_m = 1 - zet * dummy_2;
  } else {
    real dummy_0 = libm::pow(1.0 - cphm_unst * zet, 0.25);
    phi_m = 1. / dummy_0;
    const real dummy_psi = 2. * libm::log(0.5 * (1. + dummy_0)) + libm::log(0.5 * (1. + dummy_0 * dummy_0)) -
                           2. * libm::atan(dummy_0) + 1.570796;
    dummy_0 = (1. - am_unst * zet);  // BS:6373
    real dummy_1 = libm::pow(dummy_0, 0.333333);
    real dummy_11 = -0.33333 * am_unst * libm::pow(dummy_0, -0.6666667);
    real dummy_2 = 0.33333 * (dummy_1 * dummy_1 + dummy_1 + 1.);
    real dummy_22 = 0.3333 * dummy_11 * (2. * dummy_1 + 1.);
    const real dummy_3 = 0.57735 * (2. * dummy_1 + 1.);
    const real dummy_33 = 1.1547 * dummy_11;
    const real dummy_4 = 1.5 * libm::log(dummy_2) - 1.73205 * libm::atan(dummy_3) + 1.813799364;
    const real dummy_44 = (1.5 / dummy_2) * dummy_22 - 1.73205 * dummy_33 / (1. + dummy_3 * dummy_3);
    dummy_0 = zet * zet;  // BS:6383
    dummy_1 = 1. / (1. + dummy_0);
    dummy_11 = 2. * zet;
    dummy_2 = ((1 - phi_m) / zet + dummy_11 * dummy_4 + dummy_0 * dummy_44) * dummy_1;
    dummy_22 = -dummy_11 * (dummy_psi + dummy_0 * dummy_4) * (dummy_1 * dummy_1);
    phi_m = 1. - zet * (dummy_2 + dummy_22);
  }
  return phi_m;
}

// BS:6398-6444
CPPMPAS_HD inline real phih(real zet) {
  const real ah_st = 5.3, bh_st = 1.1, rbh_st = 1. / bh_st, ah_unst = 34.;
  real ph;
  if (zet >= 0.0) {
    const real dummy_0 = 1 + libm::pow(zet, bh_st);
    const real dummy_1 = zet + libm::pow(dummy_0, rbh_st);
    const real dummy_11 = 1 + libm::pow(dummy_0, rbh_st - 1) * libm::pow(zet, bh_st - 1);
    const real dummy_2 = (-ah_st / dummy_1) * dummy_11;
    ph = 1 - zet * dummy_2;
  } else {
    real dummy_0 = libm::pow(1.0 - cphh_unst * zet, 0.5);  // pow, not sqrt (object)
    const real phh = 1. / dummy_0;
    const real dummy_psi = 2. * libm::log(0.5 * (1. + dummy_0));
    dummy_0 = (1. - ah_unst * zet);  // BS:6425
    real dummy_1 = libm::pow(dummy_0, 0.333333);
    real dummy_11 = -0.33333 * ah_unst * libm::pow(dummy_0, -0.6666667);
    real dummy_2 = 0.33333 * (dummy_1 * dummy_1 + dummy_1 + 1.);
    real dummy_22 = 0.3333 * dummy_11 * (2. * dummy_1 + 1.);
    const real dummy_3 = 0.57735 * (2. * dummy_1 + 1.);
    const real dummy_33 = 1.1547 * dummy_11;
    const real dummy_4 = 1.5 * libm::log(dummy_2) - 1.73205 * libm::atan(dummy_3) + 1.813799364;
    const real dummy_44 = (1.5 / dummy_2) * dummy_22 - 1.73205 * dummy_33 / (1. + dummy_3 * dummy_3);
    dummy_0 = zet * zet;  // BS:6435
    dummy_1 = 1. / (1. + dummy_0);
    dummy_11 = 2. * zet;
    dummy_2 = ((1 - phh) / zet + dummy_11 * dummy_4 + dummy_0 * dummy_44) * dummy_1;
    dummy_22 = -dummy_11 * (dummy_psi + dummy_0 * dummy_4) * (dummy_1 * dummy_1);
    ph = 1. - zet * (dummy_2 + dummy_22);
  }
  return ph;
}

// a (ims:ime, kms:kme) array of mpas_atmphys_vars at column i, level k (1-based), or an (ims:ime) one at i
struct V3 {
  real* p = nullptr;
  int ldi = 0;
  CPPMPAS_HD real& operator()(int i, int k) const { return p[static_cast<long>(k - 1) * ldi + (i - 1)]; }
};
struct V2 {
  real* p = nullptr;
  CPPMPAS_HD real& operator()(int i) const { return p[i - 1]; }
};
struct V2i {
  int* p = nullptr;
  CPPMPAS_HD int& operator()(int i) const { return p[i - 1]; }
};

// physics_wrf/module_bl_mynn.F mynn_bl_driver (MB:19-746) for one column i (j = 1): the copies into the its:ite
// arrays, bl_mynn_pre_run (physics_wrf/bl_mynn_pre.F:63-144, BP), bl_mynn_run (BM) and bl_mynn_post_run
// (physics_wrf/bl_mynn_post.F:63-155, BO), and the copies back, with the arguments mpas_atmphys_driver_pbl.F:904-971
// passes (nbca, qoz, rnbcablten and rqozblten absent; every other optional argument present).
struct BlMynnKernel {
  Const C{};
  Options o{};
  int kte = 0;
  real delt = 0;
  int initflag = 0;
  bool do_restart = false;
  bool f_qc = false, f_qi = false, f_qs = false, f_qoz = false, f_nc = false, f_ni = false, f_nifa = false,
       f_nwfa = false, f_nbca = false;
  V2 dx, xland, ps, ts, qsfc, ust, ch, hfx, qfx, rmol, wspd, uoce, voce, znt;
  V3 dz, u, v, w, th, tt, p, exner, rho, qv, qc, qi, qs, nc, ni, nifa, nwfa, rthraten;
  V2 pblh, maxwidth, maxmf, ztop_plume;
  V2i kpbl, ktop_plume;
  V3 cldfra_bl, qc_bl, qi_bl, el_pbl, qke, qke_adv, cov, qsq, tsq, sh3d, sm3d;
  V3 rublten, rvblten, rthblten, rqvblten, rqcblten, rqiblten, rqsblten, rncblten, rniblten, rnifablten, rnwfablten;
  V3 edmf_a, edmf_w, edmf_qt, edmf_thl, edmf_ent, edmf_qc, sub_thl, sub_sqv, det_thl, det_sqv;
  V3 exch_h, exch_m, dqke, qwt, qshear, qbuoy, qdiss;

  CPPMPAS_HD void operator()(int i) const {
    const int kts = 1;
    const Const& c = C;
    // ---- MB:324-458: the its:ite copies, cloud species and number concentrations (nbca, qoz absent: 0) ----
    real qv_hv[KM], qc_hv[KM], qi_hv[KM], qs_hv[KM], nc_hv[KM], ni_hv[KM], nifa_hv[KM], nwfa_hv[KM], nbca_hv[KM],
        qoz_hv[KM];
    for (int k = kts; k <= kte; ++k) {
      qv_hv[k] = qv(i, k);
      qc_hv[k] = 0.;
      qi_hv[k] = 0.;
      qs_hv[k] = 0.;
      qoz_hv[k] = 0.;
      nc_hv[k] = 0.;
      ni_hv[k] = 0.;
      nifa_hv[k] = 0.;
      nwfa_hv[k] = 0.;
      nbca_hv[k] = 0.;
    }
    for (int k = kts; k <= kte; ++k) {  // MB:397-453
      if (f_qc) qc_hv[k] = qc(i, k);
      if (f_qi) qi_hv[k] = qi(i, k);
      if (f_qs) qs_hv[k] = qs(i, k);
      if (f_nc) nc_hv[k] = nc(i, k);
      if (f_ni) ni_hv[k] = ni(i, k);
      if (f_nifa) nifa_hv[k] = nifa(i, k);
      if (f_nwfa) nwfa_hv[k] = nwfa(i, k);
    }
    // ---- BP:96-141 bl_mynn_pre_run ----
    real sqv1[KM], sqc1[KM], sqi1[KM], sqs1[KM];
    for (int k = kts; k <= kte; ++k) {
      sqc1[k] = 0.;
      sqi1[k] = 0.;
      sqs1[k] = 0.;  // intent(out), not set without f_qs (f_qs is required)
    }
    for (int k = kts; k <= kte; ++k) sqv1[k] = qv_hv[k] / (1. + qv_hv[k]);
    if (f_qc)
      for (int k = kts; k <= kte; ++k) sqc1[k] = qc_hv[k] / (1. + qv_hv[k]);
    if (f_qi)
      for (int k = kts; k <= kte; ++k) sqi1[k] = qi_hv[k] / (1. + qv_hv[k]);
    if (f_qs)  // BP:135: qs/(1.+qs), as MPAS (not qs/(1.+qv))
      for (int k = kts; k <= kte; ++k) sqs1[k] = qs_hv[k] / (1. + qs_hv[k]);

    // ---- BM:446-1236 bl_mynn_run for column i ----
    real dqke_c[KM];
    for (int k = kts; k <= kte; ++k) dqke_c[k] = qke(i, k);  // BM:450 (present(dqke))
    const real dx1 = dx(i), xland1 = xland(i), ps1 = ps(i), ts1 = ts(i), ust1 = ust(i), hfx1 = hfx(i),
               qfx1 = qfx(i), wspd1 = wspd(i), uoce1 = uoce(i), voce1 = voce(i);
    real rmol1 = rmol(i);
    (void)qsfc;
    (void)ch;
    (void)znt;
    real dz1[KM], u1[KM], v1[KM], w1[KM], th1[KM], tk1[KM], p1[KM], ex1[KM], rho1[KM], sh1[KM], sm1[KM],
        rthraten1[KM], rstoch_col[KM], kzero[KM];
    for (int k = kts; k <= kte; ++k) {  // BM:489-513
      dz1[k] = dz(i, k);
      u1[k] = u(i, k);
      v1[k] = v(i, k);
      w1[k] = w(i, k);
      th1[k] = th(i, k);
      tk1[k] = tt(i, k);
      p1[k] = p(i, k);
      ex1[k] = exner(i, k);
      rho1[k] = rho(i, k);
      rthraten1[k] = rthraten(i, k);
      kzero[k] = 0.;
    }
    w1[kte + 1] = w(i, kte + 1);  // BM:514-516
    (void)rthraten1;
    for (int k = kts; k <= kte; ++k) rstoch_col[k] = 0.;  // BM:519-523 (spp_pbl /= 1)
    int kpbl1 = kpbl(i);  // BM:526
    real pblh1 = pblh(i);
    real cldfra_bl1[KM], qc_bl1[KM], qi_bl1[KM], el_pbl1[KM], qke1[KM], qke_adv1[KM], cov1[KM], qsq1[KM],
        tsq1[KM];
    for (int k = kts; k <= kte; ++k) {  // BM:529-544
      cldfra_bl1[k] = cldfra_bl(i, k);
      qc_bl1[k] = qc_bl(i, k);
      qi_bl1[k] = qi_bl(i, k);
      el_pbl1[k] = el_pbl(i, k);
      qke1[k] = qke(i, k);
      qke_adv1[k] = qke_adv(i, k);
      cov1[k] = cov(i, k);
      qsq1[k] = qsq(i, k);
      tsq1[k] = tsq(i, k);
      sh1[k] = sh3d(i, k);
      sm1[k] = sm3d(i, k);
    }
    int ktop_plume1 = 0;  // BM:547-551
    real maxwidth1 = 0., maxmf1 = 0., ztop_plume1 = 0.;
    real du1[KM], dv1[KM], dth1[KM], dqv1[KM], dqc1[KM], dqi1[KM], dqs1[KM], dqnc1[KM], dqni1[KM], dqnifa1[KM],
        dqnwfa1[KM], dqnbca1[KM], dqozone1[KM];
    real edmf_a1[KM], edmf_w1[KM], edmf_qt1[KM], edmf_thl1[KM], edmf_ent1[KM], edmf_qc1[KM], sub_thl1[KM],
        sub_sqv1[KM], det_thl1[KM], det_sqv1[KM];
    real exch_h1[KM], exch_m1[KM], sub_u[KM], sub_v[KM], det_sqc[KM], det_u[KM], det_v[KM];
    real s_aw1[KM], s_awthl1[KM], s_awqt1[KM], s_awqv1[KM], s_awqc1[KM], s_awu1[KM], s_awv1[KM], s_awqke1[KM],
        s_awqnc1[KM], s_awqni1[KM], s_awqnwfa1[KM], s_awqnifa1[KM], s_awqnbca1[KM];
    real sd_aw1[KM], sd_awthl1[KM], sd_awqt1[KM], sd_awqv1[KM], sd_awqc1[KM], sd_awu1[KM], sd_awv1[KM];
    for (int k = kts; k <= kte; ++k) {  // BM:553-609
      du1[k] = 0.;
      dv1[k] = 0.;
      dth1[k] = 0.;
      dqv1[k] = 0.;
      dqc1[k] = 0.;
      dqi1[k] = 0.;
      dqs1[k] = 0.;
      dqnc1[k] = 0.;
      dqni1[k] = 0.;
      dqnifa1[k] = 0.;
      dqnwfa1[k] = 0.;
      dqnbca1[k] = 0.;
      dqozone1[k] = 0.;
      edmf_a1[k] = 0.;
      edmf_w1[k] = 0.;
      edmf_qc1[k] = 0.;
      edmf_ent1[k] = 0.;
      edmf_qt1[k] = 0.;
      edmf_thl1[k] = 0.;
      sub_thl1[k] = 0.;
      sub_sqv1[k] = 0.;
      det_thl1[k] = 0.;
      det_sqv1[k] = 0.;
      exch_h1[k] = 0.;
      exch_m1[k] = 0.;
      sub_u[k] = 0.;
      sub_v[k] = 0.;
      det_sqc[k] = 0.;
      det_u[k] = 0.;
      det_v[k] = 0.;
    }
    for (int k = kts; k <= kte + 1; ++k) {  // BM:610-633
      s_aw1[k] = 0.;
      s_awthl1[k] = 0.;
      s_awqt1[k] = 0.;
      s_awqv1[k] = 0.;
      s_awqc1[k] = 0.;
      s_awu1[k] = 0.;
      s_awv1[k] = 0.;
      s_awqke1[k] = 0.;
      s_awqnc1[k] = 0.;
      s_awqni1[k] = 0.;
      s_awqnwfa1[k] = 0.;
      s_awqnifa1[k] = 0.;
      s_awqnbca1[k] = 0.;
      sd_aw1[k] = 0.;
      sd_awthl1[k] = 0.;
      sd_awqt1[k] = 0.;
      sd_awqv1[k] = 0.;
      sd_awqc1[k] = 0.;
      sd_awu1[k] = 0.;
      sd_awv1[k] = 0.;
    }
    real qv1[KM], qc1[KM], qi1[KM], qs1[KM];
    for (int k = kts; k <= kte; ++k) {  // BM:640
      qv1[k] = sqv1[k] / (1. - sqv1[k]);
      qc1[k] = sqc1[k] / (1. - sqv1[k]);
      qi1[k] = sqi1[k] / (1. - sqv1[k]);
      qs1[k] = sqs1[k] / (1. - sqv1[k]);
    }
    (void)qs1;
    real zw[KM];
    zw[kts] = 0.;  // BM:647
    for (int k = kts + 1; k <= kte + 1; ++k) zw[k] = zw[k - 1] + dz1[k - 1];
    real sqw[KM], thl[KM], thetav[KM];
    for (int k = kts; k <= kte; ++k) {  // BM:688
      sqw[k] = sqv1[k] + sqc1[k] + sqi1[k];
      thl[k] = th1[k] - c.xlvcp / ex1[k] * sqc1[k] - c.xlscp / ex1[k] * (sqi1[k]);
      thetav[k] = th1[k] * (1. + 0.608 * sqv1[k]);
    }
    real psig_bl, psig_shcu;
    const bool scaleaware = o.scaleaware == 1;
    if (initflag > 0 && !do_restart) {  // BM:703
      const bool cycling = o.do_DAcycling;
      bool initialize_qke;
      if ((do_restart || cycling)) {
        initialize_qke = qke1[kts] < 0.0002;
      } else {
        initialize_qke = true;
      }
      if (!do_restart || !cycling) {  // BM:716
        for (int k = kts; k <= kte; ++k) {
          sh1[k] = 0.;
          sm1[k] = 0.;
          el_pbl1[k] = 0.;
          tsq1[k] = 0.;
          qsq1[k] = 0.;
          cov1[k] = 0.;
          cldfra_bl1[k] = 0.;
          qc_bl1[k] = 0.;
          qi_bl1[k] = 0.;
          qke1[k] = 0.;
        }
      }
      if (initialize_qke)  // BM:736
        for (int k = kts; k <= kte; ++k)
          qke1[k] = 5. * ust1 * fmx((ust1 * 700. - zw[k]) / (fmx(ust1, 0.01) * 700.), 0.01);
      get_pblh(kte, pblh1, thetav, qke1, zw, dz1, xland1, kpbl1);  // BM:743
      if (scaleaware) {
        scale_aware(dx1, pblh1, psig_bl, psig_shcu);
      } else {
        psig_bl = 1.;
        psig_shcu = 1.;
      }
      mym_initialize(c, kte, xland1, dz1, dx1, zw, u1, v1, thl, sqv1, pblh1, th1, thetav, sh1, sm1, ust1, rmol1,
                     el_pbl1, qke1, tsq1, qsq1, cov1, psig_bl, cldfra_bl1, edmf_w1, edmf_a1, initialize_qke);
    }
    get_pblh(kte, pblh1, thetav, qke1, zw, dz1, xland1, kpbl1);  // BM:792 (tkeadvect, tkebudget, icloud_bl off)
    if (scaleaware) {
      scale_aware(dx1, pblh1, psig_bl, psig_shcu);
    } else {
      psig_bl = 1.;
      psig_shcu = 1.;
    }
    const real cpm = c.cp * (1. + 0.84 * qv1[kts]);  // BM:803 (exnerg is not used)
    const real flqv = qfx1 / rho1[kts];  // BM:813
    const real flqc = 0.0;
    const real th_sfc = ts1 / ex1[kts];
    const real flq = flqv + flqc;  // BM:818
    const real flt = hfx1 / (rho1[kts] * cpm) - c.xlvcp * flqc / ex1[kts];
    const real fltv = flt + flqv * c.p608 * th_sfc;
    rmol1 = -c.karman * c.gtr * fltv / fmx(ust1 * ust1 * ust1, 1.0e-6);  // BM:823
    real zet = 0.5 * dz1[kts] * rmol1;
    zet = fmx(zet, -20.);
    zet = fmn(zet, 20.);
    const real phi_m = phim(zet);  // BM:840 (bl_mynn_stfunc = 1)
    const real pmz = phi_m - zet;
    const real phh = phih(zet);
    (void)phi_m;
    real vt[KM], vq[KM], sgm[KM];
    for (int k = kts; k <= kte; ++k) {  // BM:849
      vt[k] = 0.;
      vq[k] = 0.;
      sgm[k] = 0.;
    }
    mym_condensation(c, kte, dx1, dz1, zw, xland1, thl, sqw, sqv1, sqc1, sqi1, sqs1, p1, ex1, tsq1, qsq1, cov1, sh1,
                     el_pbl1, qc_bl1, qi_bl1, cldfra_bl1, pblh1, hfx1, vt, vq, th1, sgm, rmol1, o.spp_pbl,
                     rstoch_col);  // BM:855
    real tkeprodtd[KM];
    for (int k = kts; k <= kte; ++k) tkeprodtd[k] = 0.;  // BM:877-881 (bl_mynn_topdown off)
    {  // BM:886-917 dmp_mf (bl_mynn_edmf on)
      MfIn mi{zw, dz1, p1, rho1, u1, v1, w1, th1, thl, thetav, tk1, sqw, sqv1, sqc1, qke1, nc_hv, ni_hv, nwfa_hv,
              nifa_hv, nbca_hv, ex1, delt, ust1, flt, fltv, flq, flqv, pblh1, dx1, xland1, th_sfc, psig_shcu, kpbl1,
              o.mixscalars == 1, rstoch_col};
      MfOut mo{vt, vq, sgm, edmf_a1, edmf_w1, edmf_qt1, edmf_thl1, edmf_ent1, edmf_qc1, s_aw1, s_awthl1, s_awqt1,
               s_awqv1, s_awqc1, s_awu1, s_awv1, s_awqke1, s_awqnc1, s_awqni1, s_awqnwfa1, s_awqnifa1, s_awqnbca1,
               sub_thl1, sub_sqv1, sub_u, sub_v, det_thl1, det_sqv1, det_sqc, det_u, det_v, qc_bl1, cldfra_bl1,
               &maxwidth1, &maxmf1, &ztop_plume1, &ktop_plume1};
      dmp_mf(c, kte, mi, mo);
    }
    const real delt2 = delt;  // BM:937
    real dfm[KM], dfh[KM], dfq[KM], tcd[KM], qcd[KM], pdk[KM], pdt[KM], pdq[KM], pdc[KM];
    mym_turbulence(c, kte, xland1, o.closure, dz1, dx1, zw, u1, v1, thl, thetav, sqc1, sqw, qke1, tsq1, qsq1, cov1, vt,
                   vq, rmol1, flt, fltv, flq, pblh1, th1, sh1, sm1, el_pbl1, dfm, dfh, dfq, tcd, qcd, pdk, pdt, pdq,
                   pdc, psig_bl, psig_shcu, cldfra_bl1, edmf_w1, edmf_a1, tkeprodtd);  // BM:939
    mym_predict(c, kte, o.closure, delt2, dz1, ust1, flt, flq, pmz, phh, el_pbl1, dfq, rho1, pdk, pdt, pdq, pdc,
                qke1, tsq1, qsq1, cov1, s_aw1, s_awqke1, o.edmf_tke == 1);  // BM:963
    real diss_heat[KM];
    if (o.dheat_opt == 1) {  // BM:972: the loop k = kts..kte-1 runs in libmvec pairs (lane)
      const int n = kte - 1;
      for (int k = kts; k <= kte - 1; ++k) {
        const bool ln = lane(k - kts, n);
        const real q15 = ln ? libm::mvec_pow(qke1[k], 1.5) : libm::pow(qke1[k], 1.5);
        diss_heat[k] = fmn(fmx(1.0 * (q15) / (b1 * fmx(0.5 * (el_pbl1[k] + el_pbl1[k + 1]), 1.)) / c.cp, 0.0), 0.002);
        const real ea = -10000. / fmx(p1[k], 1.);
        diss_heat[k] = diss_heat[k] * (ln ? libm::mvec_exp(ea) : libm::exp(ea));
      }
      diss_heat[kte] = 0.;
    } else {
      for (int k = kts; k <= kte; ++k) diss_heat[k] = 0.;
    }
    {  // BM:986 mynn_tendencies (kzero replaces qs and sqs)
      TendIn ti{dz1, rho1, u1, v1, th1, tk1, qv1, qc1, qi1, kzero, nc_hv, ni_hv, p1, ex1, sqv1, sqc1, sqi1, kzero,
                nwfa_hv, nifa_hv, nbca_hv, qoz_hv, tcd, qcd, dfm, dfh, diss_heat, s_aw1, s_awthl1, s_awqt1, s_awqv1,
                s_awqc1, s_awu1, s_awv1, s_awqnc1, s_awqni1, s_awqnwfa1, s_awqnifa1, s_awqnbca1, sd_aw1, sd_awthl1,
                sd_awqt1, sd_awqv1, sd_awqc1, sd_awu1, sd_awv1, sub_thl1, sub_sqv1, sub_u, sub_v, det_thl1,
                det_sqv1, det_sqc, det_u, det_v, delt, ps1, ust1, flt, flq, flqv, flqc, wspd1, uoce1, voce1, f_qc,
                f_qi, f_ni, f_qs, o.edmf_mom == 1, o.mixscalars == 1};
      TendOut to{thl, sqw, du1, dv1, dth1, dqv1, dqc1, dqi1, dqs1, dqni1, dqnc1, dqnwfa1, dqnifa1, dqnbca1, dqozone1};
      mynn_tendencies(c, kte, ti, to);
    }
    retrieve_exchange_coeffs(kte, dfm, dfh, dz1, exch_m1, exch_h1);  // BM:1069
    // BM:1075-1129 tendencies (bl_mynn_cloudmix, bl_mynn_mixscalars on; qnc / qnifa / qnwfa / qnbca off)
    real rqv_hv[KM], rqc_hv[KM], rqi_hv[KM], rqs_hv[KM], rni_hv[KM];
    for (int k = kts; k <= kte; ++k) {
      rqc_hv[k] = 0.;
      rqi_hv[k] = 0.;
      rqs_hv[k] = 0.;
      rni_hv[k] = 0.;
      rublten(i, k) = du1[k];
      rvblten(i, k) = dv1[k];
      rthblten(i, k) = dth1[k];
      rqv_hv[k] = dqv1[k];
    }
    for (int k = kts; k <= kte; ++k) {
      if (f_qc) rqc_hv[k] = dqc1[k];
      if (f_qi) rqi_hv[k] = dqi1[k];
      if (f_qs) rqs_hv[k] = dqs1[k];
      if (f_ni) rni_hv[k] = dqni1[k];
    }
    kpbl(i) = kpbl1;  // BM:1132
    ktop_plume(i) = ktop_plume1;
    pblh(i) = pblh1;
    for (int k = kts; k <= kte; ++k) {  // BM:1137-1153
      cldfra_bl(i, k) = cldfra_bl1[k];
      qc_bl(i, k) = qc_bl1[k];
      qi_bl(i, k) = qi_bl1[k];
      el_pbl(i, k) = el_pbl1[k];
      qke(i, k) = qke1[k];
      qke_adv(i, k) = qke_adv1[k];
      cov(i, k) = cov1[k];
      qsq(i, k) = qsq1[k];
      tsq(i, k) = tsq1[k];
      sh3d(i, k) = sh1[k];
      sm3d(i, k) = sm1[k];
    }
    for (int k = kts; k <= kte; ++k) {  // BM:1157-1192 (bl_mynn_tkebudget off)
      (void)dqke_c[k];
      qshear(i, k) = 0.;
      qbuoy(i, k) = 0.;
      qwt(i, k) = 0.;
      qdiss(i, k) = 0.;
      dqke(i, k) = 0.;
    }
    for (int k = kts; k <= kte; ++k) {  // BM:1196-1222 (bl_mynn_output off)
      edmf_a(i, k) = 0.;
      edmf_w(i, k) = 0.;
      edmf_qt(i, k) = 0.;
      edmf_thl(i, k) = 0.;
      edmf_ent(i, k) = 0.;
      edmf_qc(i, k) = 0.;
      sub_thl(i, k) = 0.;
      sub_sqv(i, k) = 0.;
      det_thl(i, k) = 0.;
      det_sqv(i, k) = 0.;
    }
    maxwidth(i) = maxwidth1;  // BM:1228
    maxmf(i) = maxmf1;
    ztop_plume(i) = ztop_plume1;
    for (int k = kts; k <= kte; ++k) {  // BM:1232
      exch_h(i, k) = exch_h1[k];
      exch_m(i, k) = exch_m1[k];
    }

    // ---- BO:96-151 bl_mynn_post_run ----
    real sqv_n[KM];
    for (int k = kts; k <= kte; ++k) {
      const real sq = qv_hv[k] / (1. + qv_hv[k]);
      sqv_n[k] = sq + rqv_hv[k] * delt;
      const real rq = sqv_n[k] / (1. - sqv_n[k]);
      rqv_hv[k] = (rq - qv_hv[k]) / delt;
    }
    if (f_qc)
      for (int k = kts; k <= kte; ++k) {
        const real sq = qc_hv[k] / (1. + qv_hv[k]);
        const real sqx = sq + rqc_hv[k] * delt;
        const real rq = sqx * (1. + sqv_n[k]);
        rqc_hv[k] = (rq - qc_hv[k]) / delt;
      }
    if (f_qi)
      for (int k = kts; k <= kte; ++k) {
        const real sq = qi_hv[k] / (1. + qv_hv[k]);
        const real sqx = sq + rqi_hv[k] * delt;
        const real rq = sqx * (1. + sqv_n[k]);
        rqi_hv[k] = (rq - qi_hv[k]) / delt;
      }
    if (f_qs)  // BO:139-147: MPAS converts rqsblten too
      for (int k = kts; k <= kte; ++k) {
        const real sq = qs_hv[k] / (1. + qv_hv[k]);
        const real sqx = sq + rqs_hv[k] * delt;
        const real rq = sqx * (1. + sqv_n[k]);
        rqs_hv[k] = (rq - qs_hv[k]) / delt;
      }
    // ---- MB:599-660: the tendencies back (f_nc, f_nifa, f_nwfa: refused by the caller) ----
    for (int k = kts; k <= kte; ++k) {
      rqvblten(i, k) = rqv_hv[k];
      if (f_qc) rqcblten(i, k) = rqc_hv[k];
      if (f_qi) rqiblten(i, k) = rqi_hv[k];
      if (f_qs) rqsblten(i, k) = rqs_hv[k];
      if (f_ni) rniblten(i, k) = rni_hv[k];
    }
  }
};
}  // namespace cppmpas::phys::mynn
