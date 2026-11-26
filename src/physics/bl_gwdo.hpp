// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// Orographic gravity wave drag of the YSU suite (MMM-physics bl_gwdo.F90, pinned commit a4baf7f): bl_gwdo_run for one
// column, statement by statement. MPAS's variant of CPPWRF's GWDO kernel (SHARED_WITH_CPPWRF.md).
//
// Source comparison: MPAS v8.4.2 builds GWDO from MMM-physics a4baf7f (cppmpas-data `inputs`,
// mmm-physics-a4baf7f.tar.gz) through physics_wrf/module_bl_gwdo.F. Its bl_gwdo.F90 is the KIM GWDO of WRF 4.6.1
// (byte-identical to WRF 4.6.1's phys/physics_mmm/bl_gwdo.F90). CPPWRF ported WRF 4.8.0's bl_gwdo.F90 (CPPWRF 7fac5bf
// src/phys/bl_gwdo/gwdo.hpp, last changed in a985afe), a later revision of the scheme that differs in most of the
// drag computation (diff of the two files without comments: the reference level kref = level of 2 var capped at
// kte / 2 with omax and kbomax for the blocking layer, against MPAS's kbl = level of 2 var capped at kte; the
// oc / ol clamps; the ldrag tests (omax, bnv2 < 0); the stress formula gfac / dx_eff / nonhydrostatic factor against
// MPAS's gfobnv / cleff = dxmeter; MPAS's "oa <= 0 or kp1 >= kpblmin" test in the saturation branch; MPAS's one
// dtfac per column below kbl against WRF 4.8.0's dtfac per level with the mesosphere limiter; the blocking-drag area
// dxmeter**2 against delx*dely). The two revisions share the column layout, the initialisation, the u1 / v1
// rotation, the wind-direction lookup, the Richardson number and Brunt-Vaisala loops, the projected wind, the base
// stress efact / coefm powers, the wave-breaking loop and the final rotation; the delta is too wide for parameters of
// CPPWRF's kernel, so this file is the MPAS revision written in CPPWRF's form: CPPWRF's GwdoColumn layout and
// conventions (level arrays 1-based, k = 1 the lowest level, interfaces 1 .. kte + 1; x**2 and x**2. are x*x,
// x**real is pow, nint is lround; c_min / c_max of CPPWRF src/core/backend.hpp:147-149), with the statements of
// MMM-physics a4baf7f bl_gwdo.F90:55-643 (cited by line, GW).
//
// Object code: bl_gwdo.F90 compiled standalone with the reference flags (gfortran 14.2.0 -O3 -ffree-line-length-none
// -fconvert=big-endian -ffree-form -fdefault-real-8 -fdefault-double-8, no -Dmpas, as the build log) is
// instruction-identical to the reference executable's (bin/atmosphere_model, sha256 524a1cf4...: bl_gwdo_run 6476
// instructions, bl_gwdo_init, bl_gwdo_finalize). No FMA, no ymm register, no libmvec call: its libm calls are
// atan2@plt (GW:338), pow@plt x2 (GW:464, 466: __ieee754_atan2_fma, __ieee754_pow_fma) and lround@plt (nint, GW:339),
// all scalar, so the column kernel is the MPAS one.
#include <cmath>

#include "core/backend.hpp"
#include "core/libm.hpp"

namespace cppmpas::phys::bl_gwdo {

inline constexpr int kMaxLevels = 128;
inline constexpr int KM = kMaxLevels + 2;

// CPPWRF 7fac5bf src/core/backend.hpp:147-149
template <class T>
CPPMPAS_HD inline T c_min(T a, T b) { return (a < b) ? a : b; }
template <class T>
CPPMPAS_HD inline T c_max(T a, T b) { return (a > b) ? a : b; }

// One column of bl_gwdo_run's arguments (CPPWRF's GwdoColumn): level arrays 1..kte (prsi 1..kte + 1).
// rublten, rvblten in: the PBL tendencies; out: + the drag
template <typename T>
struct GwdoColumn {
  T uproj[KM], vproj[KM], t1[KM], q1[KM], prslk[KM], zl[KM], prsl[KM], prsi[KM], rublten[KM], rvblten[KM];
  T dtaux3d[KM], dtauy3d[KM];
  T var, oc1, oa2d1, oa2d2, oa2d3, oa2d4, ol2d1, ol2d2, ol2d3, ol2d4, sina, cosa, dxmeter;
  T dusfcg, dvsfcg;
};

// GW:55-643 for one column (its = ite): kte levels, kme = kte + 1 (module_bl_gwdo.F passes kme = kte + 1)
template <typename T>
CPPMPAS_HD inline void bl_gwdo_run(GwdoColumn<T>& c, int kte, T g_, T cp_, T rd_, T fv_, T pi_, T deltim) {
  constexpr int kts = 1;  // GW:120
  // GW:141-154, 188-189, 192-195
  const T ric = T(0.25), dw2min = T(1.), rimin = T(-100.), bnv2min = T(1.0e-5), efmin = T(0.0), efmax = T(10.0);
  const T xl = T(4.0e4), gmax = T(1.), veleps = T(1.0), frc = T(1.0), ce = T(0.8), cg = T(0.5);
  constexpr int kpblmin = 2, mdir = 8;
  constexpr int nwdir[mdir] = {6, 7, 5, 8, 2, 3, 1, 4};
  const T frmax = T(10.), olmin = T(1.0e-5), odmin = T(0.1), odmax = T(10.);
  // critac = 1.0e-5 (GW:148) is not used

  T taup[KM], velco[KM], bnv2[KM], usqj[KM], taud[KM], rho[KM], vtk[KM], vtj[KM], del[KM], u1[KM], v1[KM];
  T dudt[KM], dvdt[KM], dtaux2d[KM], dtauy2d[KM], taufb[KM];

  // GW:212-214
  const int lcap = kte;
  const int lcapp1 = lcap + 1;
  const T fdir = T(mdir) / (T(2.0) * pi_);
  // GW:223-234
  const T delx = c.dxmeter, dely = c.dxmeter;
  T dxy4[5], dxy4p[5];
  dxy4[1] = delx;
  dxy4[2] = dely;
  dxy4[3] = std::sqrt(delx * delx + dely * dely);
  dxy4[4] = dxy4[3];
  dxy4p[1] = dxy4[2];
  dxy4p[2] = dxy4[1];
  dxy4p[3] = dxy4[4];
  dxy4p[4] = dxy4[3];
  const T cleff = c.dxmeter;
  // GW:238-253
  bool ldrag = false, icrilv = false;
  int klowtop = 0, kbl = 0;
  T dtaux = T(0.), dtauy = T(0.), xn = T(0.), yn = T(0.);
  T ubar = T(0.), vbar = T(0.), rhobar = T(0.), ulow = T(0.);
  T oa = T(0.), ol = T(0.), taub = T(0.);
  for (int k = 0; k < KM; ++k) {
    usqj[k] = T(0.); bnv2[k] = T(0.); vtj[k] = T(0.); vtk[k] = T(0.);
    taup[k] = T(0.); taud[k] = T(0.); dtaux2d[k] = T(0.); dtauy2d[k] = T(0.);
    velco[k] = T(0.);  // not set by the Fortran; read only where GW:386-396 set it
    taufb[k] = T(0.0);  // GW:252
  }
  T dtfac = T(1.0);
  T xlinv = T(1.0) / xl;
  // locals the Fortran leaves undefined where they are not used (ldrag columns)
  T coefm = T(0.), bnv = T(0.), fr = T(0.), brvf = T(0.), rulow = T(0.), od = T(0.), olp = T(0.), dxy = T(0.),
    dxyp = T(0.);
  // GW:254-263
  for (int k = kts; k <= kte; ++k) {
    vtj[k] = c.t1[k] * (T(1.) + fv_ * c.q1[k]);
    vtk[k] = vtj[k] / c.prslk[k];
    rho[k] = T(1.) / rd_ * c.prsl[k] / vtj[k];
    del[k] = c.prsi[k] - c.prsi[k + 1];
    u1[k] = c.uproj[k] * c.cosa - c.vproj[k] * c.sina;
    v1[k] = c.uproj[k] * c.sina + c.vproj[k] * c.cosa;
  }
  // GW:276-293: the level of 2 var above the lowest level
  const T zlowtop = T(2.) * c.var;
  bool kloop1 = true;
  for (int k = kts + 1; k <= kte; ++k) {
    if (zlowtop > T(0.)) {
      if (kloop1 && c.zl[k] - c.zl[1] >= zlowtop) {
        klowtop = k + 1;
        kloop1 = false;
      }
    }
  }
  // GW:295-303
  const int kpblmax = kte;
  kbl = klowtop;
  kbl = c_max(c_min(kbl, kpblmax), kpblmin);
  // komax = kbl (GW:303): not used further
  // GW:305-308
  const T delks = T(1.0) / (c.prsi[1] - c.prsi[kbl]);
  const T delks1 = T(1.0) / (c.prsl[1] - c.prsl[kbl]);
  // GW:312-321: mean below kbl
  for (int k = kts; k <= kpblmax; ++k) {
    if (k < kbl) {
      const T rcsks = del[k] * delks;
      const T rdelks = del[k] * delks;
      ubar = ubar + rcsks * u1[k];
      vbar = vbar + rcsks * v1[k];
      rhobar = rhobar + rdelks * rho[k];
    }
  }
  // GW:325-358: the low-level wind direction and the orography in it
  {
    T oa4[5], ol4[5], ol4p[5];
    oa4[1] = c.oa2d1; oa4[2] = c.oa2d2; oa4[3] = c.oa2d3; oa4[4] = c.oa2d4;
    ol4[1] = c.ol2d1; ol4[2] = c.ol2d2; ol4[3] = c.ol2d3; ol4[4] = c.ol2d4;
    const T wdir = libm::atan2(ubar, vbar) + pi_;
    const int idir = int(std::lround(fdir * wdir) % mdir) + 1;
    const int nwd = nwdir[idir - 1];
    oa = T(1 - 2 * ((nwd - 1) / 4)) * oa4[(nwd - 1) % 4 + 1];
    ol = ol4[(nwd - 1) % 4 + 1];
    ol4p[1] = ol4[2]; ol4p[2] = ol4[1]; ol4p[3] = ol4[4]; ol4p[4] = ol4[3];
    olp = ol4p[(nwd - 1) % 4 + 1];
    od = olp / c_max(ol, olmin);
    od = c_min(od, odmax);
    od = c_max(od, odmin);
    dxy = dxy4[(nwd - 1) % 4 + 1];
    dxyp = dxy4p[(nwd - 1) % 4 + 1];
  }
  // GW:362-377: Richardson number (usqj) and N**2 (bnv2) at the interfaces
  for (int k = kts; k <= kte - 1; ++k) {
    const T ti = T(2.0) / (c.t1[k] + c.t1[k + 1]);
    const T rdz = T(1.) / (c.zl[k + 1] - c.zl[k]);
    const T tem1 = u1[k] - u1[k + 1];
    const T tem2 = v1[k] - v1[k + 1];
    const T dw2 = tem1 * tem1 + tem2 * tem2;
    const T shr2 = c_max(dw2, dw2min) * rdz * rdz;
    const T bvf2 = g_ * (g_ / cp_ + rdz * (vtj[k + 1] - vtj[k])) * ti;
    usqj[k] = c_max(bvf2 / shr2, rimin);
    bnv2[k] = T(2.0) * g_ * rdz * (vtk[k + 1] - vtk[k]) / (vtk[k + 1] + vtk[k]);
  }
  // GW:381-384
  ulow = c_max(std::sqrt(ubar * ubar + vbar * vbar), T(1.0));
  rulow = T(1.) / ulow;
  // GW:386-396: the wind projected on the low-level wind
  for (int k = kts; k <= kte - 1; ++k) {
    velco[k] = T(0.5) * ((u1[k] + u1[k + 1]) * ubar + (v1[k] + v1[k + 1]) * vbar);
    velco[k] = velco[k] * rulow;
    if ((velco[k] < veleps) && (velco[k] > T(0.))) velco[k] = veleps;
  }
  // GW:400-410: no drag with a critical level below kbl
  ldrag = velco[1] <= T(0.);
  for (int k = kpblmin; k <= kpblmax; ++k)
    if (k < kbl) ldrag = ldrag || velco[k] <= T(0.);
  // GW:414-430: the mean usqj and bnv2 below kbl
  {
    const T wtkbj = (c.prsl[1] - c.prsl[2]) * delks1;
    bnv2[1] = wtkbj * bnv2[1];
    usqj[1] = wtkbj * usqj[1];
  }
  for (int k = kpblmin; k <= kpblmax; ++k) {
    if (k < kbl) {
      const T rdelks = (c.prsl[k] - c.prsl[k + 1]) * delks1;
      bnv2[1] = bnv2[1] + bnv2[k] * rdelks;
      usqj[1] = usqj[1] + usqj[k] * rdelks;
    }
  }
  // GW:433-437
  ldrag = ldrag || bnv2[1] <= T(0.0);
  ldrag = ldrag || ulow == T(1.0);
  ldrag = ldrag || c.var <= T(0.0);
  // GW:441-445
  for (int k = kpblmin; k <= kpblmax; ++k)
    if (k < kbl) usqj[k] = usqj[1];
  // GW:449-457
  if (!ldrag) {
    bnv = std::sqrt(bnv2[1]);
    fr = bnv * rulow * c.var * od;
    fr = c_min(fr, frmax);
    xn = ubar * rulow;
    yn = vbar * rulow;
  }
  // GW:461-477: the base-level stress
  if (!ldrag) {
    T efact = libm::pow(oa + T(2.), ce * fr / frc);
    efact = c_min(c_max(efact, efmin), efmax);
    coefm = libm::pow(T(1.) + ol, oa + T(1.));
    xlinv = coefm / cleff;
    const T tem = fr * fr * c.oc1;
    const T gfobnv = gmax * tem / ((tem + cg) * bnv);
    taub = xlinv * rhobar * ulow * ulow * ulow * gfobnv * efact;
  } else {
    taub = T(0.0);
    xn = T(0.0);
    yn = T(0.0);
  }
  // GW:481-485
  for (int k = kts; k <= kpblmax; ++k)
    if (k <= kbl) taup[k] = taub;
  // GW:487-531: the stress above kbl (wave breaking)
  for (int k = kpblmin; k <= kte - 1; ++k) {
    const int kp1 = k + 1;
    if (k >= kbl) {
      icrilv = icrilv || (usqj[k] < ric) || (velco[k] <= T(0.0));
      brvf = c_max(bnv2[k], bnv2min);
      brvf = std::sqrt(brvf);
    }
    if (k >= kbl && (!ldrag)) {
      if (!icrilv && taup[k] > T(0.0)) {
        const T temv = T(1.0) / velco[k];
        const T tem1 = coefm / dxy * (rho[kp1] + rho[k]) * brvf * velco[k] * T(0.5);
        T hd = std::sqrt(taup[k] / tem1);
        const T fro = brvf * hd * temv;
        const T tem2 = std::sqrt(usqj[k]);
        const T tem = T(1.) + tem2 * fro;
        const T rim = usqj[k] * (T(1.) - fro) / (tem * tem);
        if (rim <= ric) {
          if ((oa <= T(0.)) || (kp1 >= kpblmin)) {
            const T temc = T(2.0) + T(1.0) / tem2;
            hd = velco[k] * (T(2.) * std::sqrt(temc) - temc) / brvf;
            taup[kp1] = tem1 * hd * hd;
          }
        } else {
          taup[kp1] = taup[k];
        }
      }
    }
  }
  // GW:534-540
  if (lcap < kte) {
    for (int klcap = lcapp1; klcap <= kte; ++klcap) taup[klcap] = c.prsi[klcap] / c.prsi[lcap] * taup[lcap];
  }
  // GW:543-578: the flow-blocking drag
  if (!ldrag) {
    int kblk = 0;
    T fbdpe = T(0.0), fbdke = T(0.0), zblk = T(0.);
    for (int k = kte; k >= kpblmin; --k) {
      if (kblk == 0 && k <= kbl) {
        fbdpe = fbdpe + bnv2[k] * (c.zl[kbl] - c.zl[k]) * del[k] / g_ / rho[k];
        fbdke = T(0.5) * (u1[k] * u1[k] + v1[k] * v1[k]);
        if (fbdpe >= fbdke) {
          kblk = k;
          kblk = c_min(kblk, kbl);
          zblk = c.zl[kblk] - c.zl[kts];
        }
      }
    }
    if (kblk != 0) {
      const T fbdcd = c_max(T(2.0) - T(1.0) / od, T(0.0));
      taufb[kts] = T(0.5) * rhobar * coefm / (c.dxmeter * c.dxmeter) * fbdcd * dxyp * olp * zblk * (ulow * ulow);
      const T tautem = taufb[kts] / T(kblk - kts);
      for (int k = kts + 1; k <= kblk; ++k) taufb[k] = taufb[k - 1] - tautem;
      for (int k = kts; k <= kte + 1; ++k) taup[k] = taup[k] + taufb[k];
    }
  }
  // GW:582-586: - g d(tau) / dp
  for (int k = kts; k <= kte; ++k) taud[k] = T(1.) * (taup[k + 1] - taup[k]) * g_ / del[k];
  // GW:595-602: one damping factor per column against wind reversal below kbl
  for (int k = kts; k <= kpblmax - 1; ++k) {
    if (k <= kbl) {
      if (taud[k] != T(0.)) dtfac = c_min(dtfac, std::fabs(velco[k] / (deltim * taud[k])));
    }
  }
  // GW:604-622
  T dusfc = T(0.), dvsfc = T(0.);
  for (int k = kts; k <= kte; ++k) {
    taud[k] = taud[k] * dtfac;
    dtaux = taud[k] * xn;
    dtauy = taud[k] * yn;
    dtaux2d[k] = dtaux;
    dtauy2d[k] = dtauy;
    dudt[k] = dtaux;
    dvdt[k] = dtauy;
    dusfc = dusfc + dtaux * del[k];
    dvsfc = dvsfc + dtauy * del[k];
  }
  // GW:624-627
  dusfc = (-(T(1.) / g_)) * dusfc;
  dvsfc = (-(T(1.) / g_)) * dvsfc;
  // GW:631-643: back to the model grid
  for (int k = kts; k <= kte; ++k) {
    c.rublten[k] = c.rublten[k] + dudt[k] * c.cosa + dvdt[k] * c.sina;
    c.rvblten[k] = c.rvblten[k] - dudt[k] * c.sina + dvdt[k] * c.cosa;
    c.dtaux3d[k] = dtaux2d[k] * c.cosa + dtauy2d[k] * c.sina;
    c.dtauy3d[k] = -dtaux2d[k] * c.sina + dtauy2d[k] * c.cosa;
  }
  c.dusfcg = dusfc * c.cosa + dvsfc * c.sina;
  c.dvsfcg = -dusfc * c.sina + dvsfc * c.cosa;
}

}  // namespace cppmpas::phys::bl_gwdo
