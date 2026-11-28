// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// The MYNN surface layer as MPAS v8.4.2 runs it: MMM-physics a4baf7f sf_mynn.F90 (SM) sf_mynn_init and sf_mynn_run
// with the helpers on MPAS's path, and the wrapper physics_wrf/sf_mynn_pre.F (SP) sf_mynn_pre_run, for one column
// (one thread per cell). module_sf_mynn.F's sfclay_mynn (the copies around these two) is in mpas_atmphys_sfclayer.cpp.
//
// MPAS's path (every other option is a module parameter of MPAS, not a namelist option):
//  * psi_opt = 0 (SM:72): the Cheng-Brutsaert / blended COARE functions, never the GFS ones;
//  * coare_opt = 3.0 (SM:78): charnock_1955 and fairall_etal_2003, never edson_etal_2013 / fairall_etal_2014;
//  * isftcflx = 0, iz0tlnd = 0, isfflx = 1 (mpas_atmphys_driver_sfclayer.F:30-32, both optional arguments present):
//    water: charnock_1955 and fairall_etal_2003; land: andreas_2002 over snow (snowh >= 0.1) else zilitinkevich_1995
//    with czil = 0.085; never davis_etal_2008, taylor_yelland_2001, garratt_1992 or yang_2008; the flux branch, never
//    the isfflx = 0 zeros; no dissipative heating (isftcflx = 0);
//  * spp_pbl = 0 (mpas_atmphys_vars.F:417): f_spp = .false., rstoch1d = 0, zntstoch = znt;
//  * itimestep is driver_sfclayer's initflag, 0 or 1 (mpas_atmphys_driver_sfclayer.F:885-886): SM:579 and SM:664
//    (itimestep .le. 1) always take li_etal_2010's first guess; the mol-based first guess (SM:582-584, 667-669) is
//    dead code in MPAS and is not ported (sf_mynn_column refuses another itimestep);
//  * ustm, ck, cka, cd and cda are present (SM:740, 879).
// The diagnostics wstar and qstar (SM:356-365, 775) are outputs of sf_mynn_run that sfclay_mynn keeps in its local
// (ims:ime) arrays and drops; wstar enters wspd. qstar is read only by yang_2008 (iz0tlnd = 2), so the uninitialised
// local qstar of module_sf_mynn.F:215 never reaches a result.
//
// Object code: sf_mynn.F90 compiled standalone with the reference flags (gfortran 14.2.0 -O3 -ffree-line-length-none
// -fconvert=big-endian -ffree-form -fdefault-real-8 -fdefault-double-8, no -Dmpas, as the reference build log) is
// instruction-identical to the reference executable's (bin/atmosphere_model sha256 524a1cf4...), see
// tests/sfmynn_column. No FMA and no ymm register. With -g:
//  * SM:285-306 is vectorised: thgb (SM:289) and thcon (SM:292) call _ZGVbN2vv_pow for the pairs of i from its and the
//    scalar pow for an odd last i and for its = ite (the 1-trip loop): CPPWRF's lane(i, n) rule (sfclayrev::vm::pow).
//    Nothing else of the loop calls libm.
//  * every other libm call of sf_mynn_run and of the helpers is scalar (exp, log, pow, atan; the helpers charnock,
//    fairall, andreas are inlined; zilitinkevich_1995, li_etal_2010, zolrib and the psi*_stable / _unstable lookups
//    are calls with the *_full functions inlined into them).
//  * psi_init (SM:1967-2004) is inlined into sf_mynn_init and vectorised exactly as sf_sfclayrev_init is (the same
//    libmvec calls, _ZGVbN2vv_pow x14, _ZGVbN2v_log x14, _ZGVbN2v_atan x6, the same arithmetic opcodes and counts):
//    init_sfmynn_tables is init_sfclayrev_tables_o3, and sfmynn_column compares the tables of both codes.
//  * folded constants (GCC, correctly rounded; each equals what glibc's log / exp give at run time for the same
//    argument, checked bit for bit, so the C++ computes them): 1./273.15 (SM:311); log(10./1e-4) of charnock_1955
//    (SM:1165); in li_etal_2010 log(100000.) once zaz0 > 1e5 and log(100.) once z0zt > 100; in andreas_2002 the whole
//    exp(bt0_r + bt1_r*log(1000.) + bt2_r*log(1000.)**2) (and the bq one) once ren2 > 1000; 2.*ATAN(1.) and
//    4.*ATAN(1.)/sqrt(3.) of the unstable functions; e**2. (garratt, not on the path).
//  * powers: x**2 and x**2. are products (rib**2, alfa**2, beta**2, log(ren2)**2, psix2**2, zolf**2, ym**2.); every
//    other real power is a pow call, (1.-16.*zolf)**.5 included (scalar: pow, not sqrt).
//  * MAX / MIN are maxsd / minsd; no ±0 can reach them with equal operands on this path (c_max / c_min).
//  * x/8. and x/2. are multiplications by 0.125 and 0.5 (exact).
// Untyped literals are REAL(8) (-fdefault-real-8), so every constant is the double of its decimal text; the local
// parameters g = 9.81 of charnock_1955 and 9.8 of andreas_2002 are not MPAS's gravity.
//
// The MMM-physics version is MPAS's (the older single-module MYNN surface layer). CPPWRF ports WRF 4.8's split
// MYNN-SFC (src/phys/sf_mynn). Of WRF 4.8's routines (comment-stripped diff) only li_etal_2010 and zolrib are MPAS's
// statements (but for _kind_phys suffixes and the table arguments); its table lookups end at nzol + 1 < 1000 where
// MPAS's end at nzol + 1 <= 1000 (SM:2147, 2172, 2197, 2222); fairall_etal_2003 (WRF clamps zq, MPAS copies the
// clamped zt), zilitinkevich_1995 (other options) and andreas_2002 (* g_inv against MPAS's / 9.8, SPP) differ;
// charnock_1955 is not in WRF 4.8's compiled files. So no CPPWRF code is used here; sf_sfclayrev.hpp's psi*_full
// functions (the same Fortran text, SM:2009-2060) are.
#include <cmath>

#include "core/backend.hpp"
#include "core/libm.hpp"
#include "physics/sf_sfclayrev.hpp"

namespace cppmpas::phys::sfmynn {

using sfclayrev::c_max;
using sfclayrev::c_min;
CPPMPAS_HD inline real c_abs(real x) { return std::fabs(x); }  // ABS: andpd with the sign mask
using PsiTables = sfclayrev::SFClayRevTables<real>;  // psim_stab, psim_unstab, psih_stab, psih_unstab (SM:81-82)

// sf_mynn_init -> psi_init(psi_opt = 0) (SM:92-106, 1967-1989); see the header
inline void init_sfmynn_tables(PsiTables& t) { sfclayrev::init_sfclayrev_tables_o3(t); }

// SM:2135-2232 (psi_opt = 0): the tables for nzol + 1 <= 1000, the full functions beyond
CPPMPAS_HD inline real psim_stable(real zolf, const PsiTables& t) {
  const int nzol = static_cast<int>(zolf * real(100.));
  const real rzol = zolf * real(100.) - static_cast<real>(nzol);
  if (nzol + 1 <= 1000) {
    if (nzol < 0) kernel_fail("sf_mynn: psim_stable reads psim_stab below 0 (undefined in MPAS)");
    return t.psim_stab[nzol] + rzol * (t.psim_stab[nzol + 1] - t.psim_stab[nzol]);
  }
  return sfclayrev::psim_stable_full(zolf);
}
CPPMPAS_HD inline real psih_stable(real zolf, const PsiTables& t) {
  const int nzol = static_cast<int>(zolf * real(100.));
  const real rzol = zolf * real(100.) - static_cast<real>(nzol);
  if (nzol + 1 <= 1000) {
    if (nzol < 0) kernel_fail("sf_mynn: psih_stable reads psih_stab below 0 (undefined in MPAS)");
    return t.psih_stab[nzol] + rzol * (t.psih_stab[nzol + 1] - t.psih_stab[nzol]);
  }
  return sfclayrev::psih_stable_full(zolf);
}
CPPMPAS_HD inline real psim_unstable(real zolf, const PsiTables& t) {
  const int nzol = static_cast<int>(-zolf * real(100.));
  const real rzol = -zolf * real(100.) - static_cast<real>(nzol);
  if (nzol + 1 <= 1000) {
    if (nzol < 0) kernel_fail("sf_mynn: psim_unstable reads psim_unstab below 0 (undefined in MPAS)");
    return t.psim_unstab[nzol] + rzol * (t.psim_unstab[nzol + 1] - t.psim_unstab[nzol]);
  }
  return sfclayrev::psim_unstable_full(zolf);
}
CPPMPAS_HD inline real psih_unstable(real zolf, const PsiTables& t) {
  const int nzol = static_cast<int>(-zolf * real(100.));
  const real rzol = -zolf * real(100.) - static_cast<real>(nzol);
  if (nzol + 1 <= 1000) {
    if (nzol < 0) kernel_fail("sf_mynn: psih_unstable reads psih_unstab below 0 (undefined in MPAS)");
    return t.psih_unstab[nzol] + rzol * (t.psih_unstab[nzol + 1] - t.psih_unstab[nzol]);
  }
  return sfclayrev::psih_unstable_full(zolf);
}

// li_etal_2010, SM:1721-1788
CPPMPAS_HD inline void li_etal_2010(real& zl, real rib, real zaz0, real z0zt) {
  constexpr real au11 = 0.045, bu11 = 0.003, bu12 = 0.0059, bu21 = -0.0828, bu22 = 0.8845, bu31 = 0.1739,
                 bu32 = -0.9213, bu33 = -0.1057;
  constexpr real aw11 = 0.5738, aw12 = -0.4399, aw21 = -4.901, aw22 = 52.50, bw11 = -0.0539, bw12 = 1.540,
                 bw21 = -0.669, bw22 = -3.282;
  constexpr real as11 = 0.7529, as21 = 14.94, bs11 = 0.1569, bs21 = -0.3091, bs22 = -1.303;
  real zaz02 = zaz0;  // :1749-1751
  if (zaz0 < real(100.0)) zaz02 = real(100.);
  if (zaz0 > real(100000.0)) zaz02 = real(100000.);
  real z0zt2 = z0zt;  // :1754-1756
  if (z0zt < real(0.5)) z0zt2 = real(0.5);
  if (z0zt > real(100.0)) z0zt2 = real(100.);
  const real alfa = libm::log(zaz02);  // :1758-1759 (log(1e5), log(100.) folded on the clamped paths: same bits)
  const real beta = libm::log(z0zt2);
  if (rib <= real(0.0)) {  // :1761-1769
    zl = au11 * alfa * (rib * rib) +
         ((bu11 * beta + bu12) * (alfa * alfa) + (bu21 * beta + bu22) * alfa + (bu31 * (beta * beta) + bu32 * beta + bu33)) *
             rib;
    zl = c_max(zl, real(-15.));
    zl = c_min(zl, real(0.));
  } else if (rib > real(0.0) && rib <= real(0.2)) {  // :1770-1778
    zl = ((aw11 * beta + aw12) * alfa + (aw21 * beta + aw22)) * (rib * rib) +
         ((bw11 * beta + bw12) * alfa + (bw21 * beta + bw22)) * rib;
    zl = c_min(zl, real(4.));
    zl = c_max(zl, real(0.));
  } else {  // :1779-1786
    zl = (as11 * alfa + as21) * rib + bs11 * alfa + bs21 * beta + bs22;
    zl = c_min(zl, real(20.));
    zl = c_max(zl, real(1.));
  }
}

// zolrib, SM:1894-1964 (psi_opt = 0); zol1 is inout as in MPAS
CPPMPAS_HD inline real zolrib(real ri, real za, real z0, real zt, real logz0, real logzt, real& zol1,
                              const PsiTables& t) {
  constexpr int nmax = 20;
  if (zol1 * ri < real(0.)) zol1 = real(0.);  // :1915-1918
  real zolold, res;
  if (ri < real(0.)) {  // :1920-1926
    zolold = real(-99999.);
    res = real(-66666.);
  } else {
    zolold = real(99999.);
    res = real(66666.);
  }
  int n = 1;
  while (c_abs(zolold - res) > real(0.01) && n < nmax) {  // :1929-1951
    if (n == 1)
      zolold = zol1;
    else
      zolold = res;
    const real zol20 = zolold * z0 / za;
    const real zol3 = zolold + zol20;
    const real zolt = zolold * zt / za;
    real psit2, psix2;
    if (ri < real(0)) {
      psit2 = c_max(logzt - (psih_unstable(zol3, t) - psih_unstable(zolt, t)), real(1.0));
      psix2 = c_max(logz0 - (psim_unstable(zol3, t) - psim_unstable(zol20, t)), real(1.0));
    } else {
      psit2 = c_max(logzt - (psih_stable(zol3, t) - psih_stable(zolt, t)), real(1.0));
      psix2 = c_max(logz0 - (psim_stable(zol3, t) - psim_stable(zol20, t)), real(1.0));
    }
    res = ri * (psix2 * psix2) / psit2;
    n = n + 1;
  }
  if (n == nmax && c_abs(zolold - res) > real(0.01))  // :1953-1956
    li_etal_2010(res, ri, za / z0, z0 / zt);
  return res;
}

// charnock_1955, SM:1144-1172 (inlined in MPAS's object)
CPPMPAS_HD inline void charnock_1955(real& z_0, real ustar, real wsp10, real visc, real zu) {
  constexpr real g = 9.81, czo2 = 0.011;
  const real wsp10m = wsp10 * libm::log(real(10.) / real(1e-4)) / libm::log(zu / real(1e-4));  // log(1e5) folded
  const real czc = czo2 + real(0.007) * c_min(c_max((wsp10m - real(10.)) / real(8.), real(0.)), real(1.0));
  z_0 = czc * ustar * ustar / g + (real(0.11) * visc / c_max(ustar, real(0.05)));
  z_0 = c_max(z_0, real(1.27e-7));
  z_0 = c_min(z_0, real(2.85e-3));
}

// fairall_etal_2003, SM:1253-1302 (spp_pbl false); zq = max(zt, ...) as written (MPAS: zq = zt)
CPPMPAS_HD inline void fairall_etal_2003(real& zt, real& zq, real ren) {
  if (ren <= real(2.)) {
    zt = real(5.5e-5) * libm::pow(ren, real(-0.60));
    zq = zt;
  } else {
    zt = real(5.5e-5) * libm::pow(ren, real(-0.60));
    zq = zt;
  }
  zt = c_min(zt, real(1.0e-4));
  zt = c_max(zt, real(2.0e-9));
  zq = c_min(zt, real(1.0e-4));
  zq = c_max(zt, real(2.0e-9));
}

// zilitinkevich_1995, SM:999-1073, called over land only (landsea = xland < 1.5) with iz0tlnd2 = 0 and spp_pbl false:
// the land branch with czil = 0.085
CPPMPAS_HD inline void zilitinkevich_1995(real z_0, real& zt, real& zq, real restar, real karman, real landsea) {
  if (landsea - real(1.5) > real(0)) kernel_fail("sf_mynn: zilitinkevich_1995's water branch has no port");
  const real czil = 0.085;  // :1054
  zt = z_0 * libm::exp(-karman * czil * std::sqrt(restar));
  zt = c_min(zt, real(0.75) * z_0);
  zq = z_0 * libm::exp(-karman * czil * std::sqrt(restar));
  zq = c_min(zq, real(0.75) * z_0);
}

// andreas_2002, SM:1403-1457 (inlined in MPAS's object; with ren2 clamped to 1000 GCC folds the exp of the rough
// branch, which equals the run-time value bit for bit)
CPPMPAS_HD inline void andreas_2002(real z_0, real bvisc, real ustar, real& zt, real& zq) {
  (void)z_0;
  constexpr real bt0_s = 1.25, bt0_t = 0.149, bt0_r = 0.317, bt1_s = 0.0, bt1_t = -0.55, bt1_r = -0.565,
                 bt2_s = 0.0, bt2_t = 0.0, bt2_r = -0.183;
  constexpr real bq0_s = 1.61, bq0_t = 0.351, bq0_r = 0.396, bq1_s = 0.0, bq1_t = -0.628, bq1_r = -0.512,
                 bq2_s = 0.0, bq2_t = 0.0, bq2_r = -0.180;
  const real zntsno = real(0.135) * bvisc / ustar +
                      (real(0.035) * (ustar * ustar) / real(9.8)) *
                          (real(5.) * libm::exp(-(((ustar - real(0.18)) / real(0.1)) * ((ustar - real(0.18)) / real(0.1)))) +
                           real(1.));
  real ren2 = ustar * zntsno / bvisc;
  if (ren2 > real(1000.)) ren2 = real(1000.);
  const real l = libm::log(ren2);  // log(ren2) of both lines (one call in MPAS's object)
  if (ren2 <= real(0.135)) {
    zt = zntsno * libm::exp(bt0_s + bt1_s * l + bt2_s * (l * l));
    zq = zntsno * libm::exp(bq0_s + bq1_s * l + bq2_s * (l * l));
  } else if (ren2 > real(0.135) && ren2 < real(2.5)) {
    zt = zntsno * libm::exp(bt0_t + bt1_t * l + bt2_t * (l * l));
    zq = zntsno * libm::exp(bq0_t + bq1_t * l + bq2_t * (l * l));
  } else {
    zt = zntsno * libm::exp(bt0_r + bt1_r * l + bt2_r * (l * l));
    zq = zntsno * libm::exp(bq0_r + bq1_r * l + bq2_r * (l * l));
  }
}

// The column of one sfclay_mynn call: sf_mynn_pre_run's level-1/2 slices (SP:127-159) of cell i, then
// sf_mynn_run (SM:130-996) for that cell. ii = i - its and n = ite - its + 1 place the cell in the vectorised loop
// SM:285 (lanes). Arguments in MPAS's names; the inout and out scalars are the cell's elements of the *_hv arrays.
struct SfMynnColumnIn {
  // SP: the slices of dz8w, u3d, v3d, p3d (pres_hyd_p), t3d, rho3d, qv3d at kts and kts + 1
  real dz1, dz2, u1, u2, v1, v2, p1, t1, rho1, qv1;
  // sfclay_mynn's 2D inputs (module_sf_mynn.F:252-261)
  real mavail, pblh, xland, tsk, psfcpa, snowh, dx;
};
struct SfMynnColumnInOut {
  real regime, hfx, qfx, lh, mol, rmol, qgh, qsfc, znt, zol, ust, cpm, chs2, cqs2, chs, ch, flhc, flqc, gz1oz0,
      wspd, br, psim, psih, ustm;                   // inout (ustm: in through module_sf_mynn.F:306-310)
  real u10, v10, th2, t2, q2, ck, cka, cd, cda;  // out
};
struct SfMynnConsts {
  real cp, g, rovcp, r, xlv, svp1, svp2, svp3, svpt0, ep1, ep2, karman;
};

CPPMPAS_HD inline void sf_mynn_column(const SfMynnColumnIn& in, SfMynnColumnInOut& o, const SfMynnConsts& c,
                                      int itimestep, int ii, int n, const PsiTables& tb) {
  constexpr real wmin = 0.1, vconvc = 1.25, prt = 1.;  // SM:75-76, 155
  if (itimestep != 0 && itimestep != 1) kernel_fail("sf_mynn: itimestep (initflag) other than 0 or 1 has no port");
  // sf_mynn_pre_run SP:125-159
  const real dz8w1d = in.dz1, u1d = in.u1, v1d = in.v1, qv1d = in.qv1, p1d = in.p1, t1d = in.t1, rho1d = in.rho1;
  const real dz2w1d = in.dz2, u1d2 = in.u2, v1d2 = in.v2;
  if (itimestep == 1) {  // SP:152-159 (qstar(i) = 0 is dropped: see the header)
    o.ust = c_max(real(0.04) * std::sqrt(u1d * u1d + v1d * v1d), real(0.001));
    o.mol = real(0.);
    o.qsfc = qv1d / (real(1.) + qv1d);
  }
  const real mavail = in.mavail, pblh = in.pblh, xland = in.xland, tsk = in.tsk, psfcpa = in.psfcpa,
             snowh = in.snowh, dx = in.dx;
  const real ep1 = c.ep1, ep2 = c.ep2, karman = c.karman, g = c.g, cp = c.cp, rovcp = c.rovcp;

  // sf_mynn_run
  const real ep3 = real(1.) - ep2;  // :283
  // :285-306 (vectorised: the two pows are lanes)
  const real psfc = psfcpa / real(1000.);
  const real thgb = tsk * sfclayrev::vm::pow(real(100.) / psfc, rovcp, ii, n);
  real pl = p1d / real(1000.);
  const real thcon = sfclayrev::vm::pow(real(100.) / pl, rovcp, ii, n);
  const real th1d = t1d * thcon;
  const real tc1d = t1d - real(273.15);
  const real qvsh = qv1d / (real(1.) + qv1d);
  const real tvcon = (real(1.) + ep1 * qvsh);
  const real thv1d = th1d * tvcon;
  const real za = real(0.5) * dz8w1d;
  const real za2 = dz8w1d + real(0.5) * dz2w1d;
  const real govrth = g / th1d;

  // :308-339
  real e1, qsfcmr;
  if (tsk < real(273.15))
    e1 = c.svp1 * libm::exp(real(4648) * (real(1.) / real(273.15) - real(1.) / tsk) -
                            real(11.64) * libm::log(real(273.15) / tsk) + real(0.02265) * (real(273.15) - tsk));
  else
    e1 = c.svp1 * libm::exp(c.svp2 * (tsk - c.svpt0) / (tsk - c.svp3));
  if (xland > real(1.5) || o.qsfc <= real(0.0)) {
    o.qsfc = ep2 * e1 / (psfc - ep3 * e1);
    qsfcmr = ep2 * e1 / (psfc - e1);
  } else {
    qsfcmr = o.qsfc / (real(1.) - o.qsfc);
  }
  if (tsk < real(273.15))  // :327-334 (the ice formula follows tsk, the argument t1d)
    e1 = c.svp1 * libm::exp(real(4648) * (real(1.) / real(273.15) - real(1.) / t1d) -
                            real(11.64) * libm::log(real(273.15) / t1d) + real(0.02265) * (real(273.15) - t1d));
  else
    e1 = c.svp1 * libm::exp(c.svp2 * (t1d - c.svpt0) / (t1d - c.svp3));
  pl = p1d / real(1000.);
  o.qgh = ep2 * e1 / (pl - e1);
  o.cpm = cp * (real(1.) + real(0.84) * qv1d);

  // :341-393
  o.wspd = std::sqrt(u1d * u1d + v1d * v1d);
  const real thvgb = thgb * (real(1.) + ep1 * o.qsfc);
  const real dthvdz = (thv1d - thvgb);
  const real fluxc = c_max(o.hfx / rho1d / cp + ep1 * thvgb * o.qfx / rho1d, real(0.));
  real wstar;
  if (xland > real(1.5) || o.qsfc <= real(0.0))
    wstar = vconvc * libm::pow(g / tsk * pblh * fluxc, real(.33));
  else
    wstar = vconvc * libm::pow(g / tsk * c_min(real(1.5) * pblh, real(4000.)) * fluxc, real(.33));
  const real vsgd = real(0.32) * libm::pow(c_max(dx / real(5000.) - real(1.), real(0.)), real(.33));
  o.wspd = std::sqrt(o.wspd * o.wspd + wstar * wstar + vsgd * vsgd);
  o.wspd = c_max(o.wspd, wmin);
  o.br = govrth * za * dthvdz / (o.wspd * o.wspd);
  if (itimestep == 1) {
    o.br = c_max(o.br, real(-2.0));
    o.br = c_min(o.br, real(2.0));
  } else {
    o.br = c_max(o.br, real(-4.0));
    o.br = c_min(o.br, real(4.0));
  }

  // :404-787
  const real visc = real(1.326e-5) * (real(1.) + real(6.542e-3) * tc1d + real(8.301e-6) * tc1d * tc1d -
                                      real(4.84e-9) * tc1d * tc1d * tc1d);
  real zntstoch, restar, z_t, z_q;
  if ((xland - real(1.5)) >= real(0)) {  // water: isftcflx = 0, coare_opt = 3.0
    charnock_1955(o.znt, o.ust, o.wspd, visc, za);
    zntstoch = o.znt;  // spp_pbl false
    restar = c_max(o.ust * zntstoch / visc, real(0.1));
    fairall_etal_2003(z_t, z_q, restar);
  } else {  // land
    zntstoch = o.znt;
    restar = c_max(o.ust * zntstoch / visc, real(0.1));
    if (snowh >= real(0.1))
      andreas_2002(zntstoch, visc, o.ust, z_t, z_q);
    else
      zilitinkevich_1995(zntstoch, z_t, z_q, restar, karman, xland);  // iz0tlnd = 0
  }
  const real zratio = zntstoch / z_t;  // :539

  o.gz1oz0 = libm::log((za + zntstoch) / zntstoch);  // :541-546
  real gz1ozt = libm::log((za + zntstoch) / z_t);
  const real gz2oz0 = libm::log((real(2.0) + zntstoch) / zntstoch);
  real gz2ozt = libm::log((real(2.0) + zntstoch) / z_t);
  const real gz10oz0 = libm::log((real(10.) + zntstoch) / zntstoch);
  const real gz10ozt = libm::log((real(10.) + zntstoch) / z_t);
  (void)gz2oz0;

  real psim10 = 0, psih10 = 0, psih2 = 0;  // set by every branch MPAS takes (br is not NaN)
  if (o.br > real(0.0)) {  // :569-635
    if (o.br > real(0.2))
      o.regime = real(1.);
    else
      o.regime = real(2.);
    li_etal_2010(o.zol, o.br, za / zntstoch, zratio);  // itimestep <= 1
    o.zol = zolrib(o.br, za, zntstoch, z_t, o.gz1oz0, gz1ozt, o.zol, tb);
    // SM:591 max(zol, 0.0): in the reference object code (both copies of the loop) maxsd has +0.0 as the
    // destination and zol as the source, so zol = -0.0 (and NaN) is kept: 0 > zol ? 0 : zol
    o.zol = c_max(real(0.0), o.zol);
    o.zol = c_min(o.zol, real(20.));
    const real zolzt = o.zol * z_t / za;
    const real zolz0 = o.zol * zntstoch / za;
    const real zolza = o.zol * (za + zntstoch) / za;
    const real zol10 = o.zol * (real(10.) + zntstoch) / za;
    const real zol2 = o.zol * (real(2.) + zntstoch) / za;
    // water and land: the same table functions (:601-627)
    o.psim = psim_stable(zolza, tb) - psim_stable(zolz0, tb);
    o.psih = psih_stable(zolza, tb) - psih_stable(zolzt, tb);
    psim10 = psim_stable(zol10, tb) - psim_stable(zolz0, tb);
    psih10 = psih_stable(zol10, tb) - psih_stable(zolz0, tb);
    psih2 = psih_stable(zol2, tb) - psih_stable(zolz0, tb);
    o.rmol = o.zol / za;
  } else if (o.br == real(0.)) {  // :637-655
    o.regime = real(3.);
    o.psim = real(0.0);
    o.psih = o.psim;
    psim10 = real(0.);
    psih10 = real(0.);
    psih2 = real(0.);
    o.zol = real(0.);
    o.rmol = o.zol / za;
  } else if (o.br < real(0.)) {  // :657-723
    o.regime = real(4.);
    li_etal_2010(o.zol, o.br, za / zntstoch, zratio);  // itimestep <= 1
    o.zol = zolrib(o.br, za, zntstoch, z_t, o.gz1oz0, gz1ozt, o.zol, tb);
    o.zol = c_max(o.zol, real(-20.0));
    o.zol = c_min(o.zol, real(0.0));
    const real zolzt = o.zol * z_t / za;
    const real zolz0 = o.zol * zntstoch / za;
    const real zolza = o.zol * (za + zntstoch) / za;
    const real zol10 = o.zol * (real(10.) + zntstoch) / za;
    const real zol2 = o.zol * (real(2.) + zntstoch) / za;
    o.psim = psim_unstable(zolza, tb) - psim_unstable(zolz0, tb);
    o.psih = psih_unstable(zolza, tb) - psih_unstable(zolzt, tb);
    psim10 = psim_unstable(zol10, tb) - psim_unstable(zolz0, tb);
    psih10 = psih_unstable(zol10, tb) - psih_unstable(zolz0, tb);
    psih2 = psih_unstable(zol2, tb) - psih_unstable(zolz0, tb);
    o.psih = c_min(o.psih, real(0.9) * gz1ozt);
    o.psim = c_min(o.psim, real(0.9) * o.gz1oz0);
    psih2 = c_min(psih2, real(0.9) * gz2ozt);
    psim10 = c_min(psim10, real(0.9) * gz10oz0);
    psih10 = c_min(psih10, real(0.9) * gz10ozt);
    o.rmol = o.zol / za;
  } else {
    kernel_fail("sf_mynn: br is NaN (MPAS keeps psim10 / psih10 / psih2 undefined)");
  }

  // :731-775
  real psix = o.gz1oz0 - o.psim;
  o.ust = real(0.5) * o.ust + real(0.5) * karman * o.wspd / psix;
  const real wspdi = c_max(std::sqrt(u1d * u1d + v1d * v1d), wmin);
  o.ustm = real(0.5) * o.ustm + real(0.5) * karman * wspdi / psix;
  if ((xland - real(1.5)) < real(0.)) {
    o.ust = c_max(o.ust, real(0.005));
    o.ustm = o.ust;
  }
  gz1ozt = libm::log((za + zntstoch) / z_t);
  gz2ozt = libm::log((real(2.0) + zntstoch) / z_t);
  real dtg = thv1d - thvgb;
  o.mol = karman * dtg / c_max(gz1ozt - o.psih, real(1.)) / prt;
  // qstar = karman*dqg/psiq/prt (:774-775) is dropped (see the header); psiq2 / psiq10 there are unused

  // :799-933
  psix = o.gz1oz0 - o.psim;
  const real psix10 = gz10oz0 - psim10;
  const real psit = c_max(gz1ozt - o.psih, real(1.0));
  const real psit2 = c_max(gz2ozt - psih2, real(1.0));
  const real psiq = c_max(libm::log((za + z_q) / z_q) - o.psih, real(1.0));
  const real psiq2 = c_max(libm::log((real(2.0) + z_q) / z_q) - psih2, real(1.0));
  const real psiq10 = c_max(libm::log((real(10.0) + z_q) / z_q) - psih10, real(1.0));
  // isfflx = 1 (:831-887)
  o.flqc = rho1d * mavail * o.ust * karman / psiq;
  o.flhc = rho1d * o.cpm * o.ust * karman / psit;
  o.qfx = o.flqc * (qsfcmr - qv1d);
  o.qfx = c_max(o.qfx, real(-0.02));
  o.lh = c.xlv * o.qfx;
  if (xland - real(1.5) > real(0.)) {  // isftcflx = 0: no dissipative heating
    o.hfx = o.flhc * (thgb - th1d);
  } else if (xland - real(1.5) < real(0.)) {
    o.hfx = o.flhc * (thgb - th1d);
    o.hfx = c_max(o.hfx, real(-250.));
  }
  o.chs = o.ust * karman / psit;
  o.ch = o.flhc / (o.cpm * rho1d);
  o.cqs2 = o.ust * karman / psiq2;
  o.chs2 = o.ust * karman / psit2;
  o.ck = (karman / psix10) * (karman / psiq10);
  o.cd = (karman / psix10) * (karman / psix10);
  o.cka = (karman / psix) * (karman / psiq);
  o.cda = (karman / psix) * (karman / psix);
  // :896-917
  if (za <= real(7.0)) {
    if (za2 > real(7.0) && za2 < real(13.0)) {
      o.u10 = u1d2;
      o.v10 = v1d2;
    } else {
      o.u10 = u1d * libm::log(real(10.) / zntstoch) / libm::log(za / zntstoch);
      o.v10 = v1d * libm::log(real(10.) / zntstoch) / libm::log(za / zntstoch);
    }
  } else if (za > real(7.0) && za < real(13.0)) {
    o.u10 = u1d * libm::log(real(10.) / zntstoch) / libm::log(za / zntstoch);
    o.v10 = v1d * libm::log(real(10.) / zntstoch) / libm::log(za / zntstoch);
  } else {
    o.u10 = u1d * psix10 / psix;
    o.v10 = v1d * psix10 / psix;
  }
  // :923-935
  dtg = th1d - thgb;
  o.th2 = thgb + dtg * psit2 / psit;
  if ((th1d > thgb && (o.th2 < thgb || o.th2 > th1d)) || (th1d < thgb && (o.th2 > thgb || o.th2 < th1d)))
    o.th2 = thgb + real(2.) * (th1d - thgb) / za;
  o.t2 = o.th2 * libm::pow(psfc / real(100.), rovcp);
  o.q2 = qsfcmr + (qv1d - qsfcmr) * psiq2 / psiq;
  o.q2 = c_max(o.q2, c_min(qsfcmr, qv1d));
  o.q2 = c_min(o.q2, real(1.05) * qv1d);
}

}  // namespace cppmpas::phys::sfmynn
