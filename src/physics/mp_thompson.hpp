// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// MPAS v8.4.2's Thompson microphysics, physics/physics_wrf/module_mp_thompson.F (MT), the "reference: WRF-v4.1.4"
// copy of the scheme MPAS builds with -Dmpas (Nt_c / mu_c per column from mp_gt_driver's ntc / muc arguments, the
// lookup tables read from MPAS's MP_THOMPSON_*_DATA.DBL, GAMMLN / WGAMMA / RSLF / RSIF of
// mpas_atmphys_functions.F). Statement by statement, for one column (kts = 1 .. kte), on the host and the device
// (CPPMPAS_HD): mp_thompson (MT:1400-3642) and calc_refl10cm (MT:4924-5235). The module values thompson_init sets
// (MT:433-993) come in a ThompsonView (physics/mp_thompson_init.hpp builds them).
//
// MPAS-own code: this is not CPPWRF's src/phys/mp_thompson (WRF 4.8.0's scheme: hail-aware graupel densities, 37
// graupel and snow table bins, other process rates; `diff` of the two module_mp_thompson.F: 2233 differing lines).
// CPPWRF's kernel was the guide to the layout (one column per call, the column arrays on the stack, the libmvec
// lane rule, the table index blocks), nothing of it is shared.
//
// The reference executable (bin/atmosphere_model sha256 524a1cf4..., __module_mp_thompson_MOD_mp_thompson, 11814
// instructions; MT compiled alone with the reference compiler and flags gives the same instructions) was read for
// the libm calls: 113 pow@plt, 11 __powidf2 (Fortran's REAL**INTEGER, `powi` below), 17 log10, 6 log, 5 exp, 24
// lround (NINT), 4 RSLF and 1 RSIF calls, and 5 _ZGVbN2vv_pow. A -g build of the same object puts the vector pows at
// MT:1941, 1944 (rain slopes, the k loop kte .. kts, both pows), MT:3021, 3024 (the same loop after the tendency
// update) and MT:3251 (rain sedimentation: the two (lamr+fv_r)**(-cre) of one k paired into one call, SLP). Every
// other libm call is scalar. No FMA, no ymm register. With -fdefault-real-8 / -fdefault-double-8 every REAL, DOUBLE
// PRECISION, literal and SNGL / DBLE is a double. x**2.0 (REAL constant 2) is x*x (GCC folds pow(x, 2.0)); x**0.0 is
// 1; the other REAL powers stay pow calls. MAX / MIN: c_max / c_min with the Fortran operand order.
//
// Not ported (dead or refused): the aerosol-aware branches (is_aerosol_aware = .false. for mp_thompson; the
// aerosol scavenging rates pna_* / pnd_* with Eff_aero and the nwfa / nifa tendencies feed only nwfa1d / nifa1d,
// which mp_gt_driver writes back only when aerosol-aware), calc_effectRad (has_reqc = 0: config_microp_re is
// refused), iiwarm = .false. (a PARAMETER).
#include <cmath>
#include <cstdint>

#include "core/backend.hpp"
#include "core/libm.hpp"
#include "physics/mp_radar.hpp"

namespace cppmpas::phys {
// physics/mpas_atmphys_functions.F:156-186 / 188-218 (mpas_atmphys_microphysics.cpp)
double rslf(double p, double t);
double rsif(double p, double t);
}  // namespace cppmpas::phys

namespace cppmpas::phys::thompson {

inline constexpr int kMaxLevels = 128;

// ---- module parameters (MT:84-363) ----
inline constexpr double T_0 = 273.15;
inline constexpr double PI = 3.1415926536;
inline constexpr double rho_w = 1000.0, rho_s = 100.0, rho_g = 500.0, rho_i = 890.0;
inline constexpr double Nt_c_max = 1999.E6;
inline constexpr double naIN0 = 1.5E6, naIN1 = 0.5E6, naCCN0 = 300.0E6, naCCN1 = 50.0E6;
inline constexpr double mu_r = 0.0, mu_g = 0.0, mu_i = 0.0;
inline constexpr double mu_s = 0.6357, Kap0 = 490.6, Kap1 = 17.46, Lam0 = 20.78, Lam1 = 3.29;
inline constexpr double gonv_min = 1.E4, gonv_max = 3.E6;
inline constexpr double am_r = PI * rho_w / 6.0, bm_r = 3.0;
inline constexpr double am_s = 0.069, bm_s = 2.0;
inline constexpr double am_g = PI * rho_g / 6.0, bm_g = 3.0;
inline constexpr double am_i = PI * rho_i / 6.0, bm_i = 3.0;
inline constexpr double av_r = 4854.0, bv_r = 1.0, fv_r = 195.0;
inline constexpr double av_s = 40.0, bv_s = 0.55, fv_s = 100.0;
inline constexpr double av_g = 442.0, bv_g = 0.89;
inline constexpr double av_i = 1847.5, bv_i = 1.0;
inline constexpr double av_c = 0.316946E8, bv_c = 2.0;
inline constexpr double C_cube = 0.5, C_sqrd = 0.15;
inline constexpr double Ef_si = 0.05, Ef_rs = 0.95, Ef_rg = 0.75, Ef_ri = 0.95;
inline constexpr double R1 = 1.E-12, R2 = 1.E-6, eps = 1.E-15;
inline constexpr double TNO = 5.0, ATO = 0.304;
inline constexpr double rho_not = 101325.0 / (287.05 * 298.0);
inline constexpr double Sc = 0.632;
inline constexpr double HGFR = 235.16;
inline constexpr double Rv = 461.5, oRv = 1. / Rv, R = 287.04, Cp = 1004.0;
inline constexpr double lsub = 2.834E6, lvap0 = 2.5E6, lfus = lsub - lvap0, olfus = 1. / lfus;
inline constexpr double xm0i = 1.E-12, D0c = 1.E-6, D0r = 50.E-6, D0s = 200.E-6, D0g = 250.E-6;
inline constexpr int nbins = 100, nbc = nbins, nbi = nbins, nbr = nbins, nbs = nbins, nbg = nbins;
inline constexpr int ntb_c = 37, ntb_i = 64, ntb_r = 37, ntb_s = 28, ntb_g = 28, ntb_g1 = 28, ntb_r1 = 37,
                     ntb_i1 = 55, ntb_t = 9, ntb_IN = 55;
inline constexpr int IFDRY = 0;
// sa / sb (MT:332-337): the snow moment coefficients (Field et al. 2005), 1-based as in MT
inline constexpr double sa[11] = {0., 5.065339, -0.062659, -3.032362, 0.029469, -0.000285,
                                  0.31255, 0.000204, 0.003199, 0.0, -0.015952};
inline constexpr double sb[11] = {0., 0.476221, -0.015896, 0.165977, 0.007468, -0.000141,
                                  0.060366, 0.000079, 0.000594, 0.0, -0.003577};
// r_c(1), r_i(1), r_r(1), r_s(1), r_g(1), Nt_i(1), Nt_IN(1) (MT:252-325): the lower ends of the table ranges
inline constexpr double r_c1 = 1.e-6, r_i1 = 1.e-10, r_r1 = 1.e-6, r_s1 = 1.e-5, r_g1 = 1.e-5, Nt_i1 = 1.0,
                        Nt_IN1 = 1.0;

// ---- the module variables thompson_init sets (MT:387-417, 491-778) and the lookup tables (MT:362-384) ----
// Fortran layouts, 1-based accessors. Trivially copyable: the tables are pointers (host or device memory).
struct ThompsonView {
  double cce[5 * 15], ccg[5 * 15], ocg1[15], ocg2[15];
  double cie[7], cig[7], oig1, oig2, obmi;
  double cre[13], crg[13], ore1, org1, org2, org3, obmr;
  double cse[18], csg[18], oams, obms, ocms;
  double cge[12], cgg[12], oge1, ogg1, ogg2, ogg3, oamg, obmg, ocmg;
  double t1_qr_qc, t1_qr_qi, t2_qr_qi, t1_qg_qc, t1_qs_qc, t1_qs_qi, t1_qr_ev, t2_qr_ev;
  double t1_qs_sd, t2_qs_sd, t1_qg_sd, t2_qg_sd, t1_qs_me, t2_qs_me, t1_qg_me, t2_qg_me;
  double Sc3, D0i, xm0s, xm0g;
  int nic1, nic2, nii2, nii3, nir2, nir3, nis2, nig2, nig3, niIN2;
  double Dc[nbc], dtc[nbc], Di[nbi], dti[nbi], Dr[nbr], dtr[nbr], Ds[nbs], dts[nbs], Dg[nbg], dtg[nbg];
  double t_Nc[nbc];
  // tables (REAL(KIND=R8SIZE)): (ntb_g1, ntb_g, ntb_r1, ntb_r), (ntb_s, ntb_t, ntb_r1, ntb_r),
  // (ntb_r, ntb_r1, 45, ntb_IN), (ntb_c, nbc, 45, ntb_IN), (ntb_i, ntb_i1), (nbr, nbc), (nbs, nbc)
  const double *tcg_racg, *tmr_racg, *tcr_gacr, *tmg_gacr, *tnr_racg, *tnr_gacr;
  const double *tcs_racs1, *tmr_racs1, *tcs_racs2, *tmr_racs2, *tcr_sacr1, *tms_sacr1, *tcr_sacr2, *tms_sacr2,
      *tnr_racs1, *tnr_racs2, *tnr_sacr1, *tnr_sacr2;
  const double *tpi_qrfz, *tni_qrfz, *tpg_qrfz, *tnr_qrfz, *tpi_qcfz, *tni_qcfz;
  const double *tpi_ide, *tps_iaus, *tni_iaus;
  const double *t_Efrw, *t_Efsw;
  // the radar values thompson_init's radar_init call leaves (MT:978-988, mp_radar's module variables)
  radar::RadarTables radar;

  CPPMPAS_HD double CCE(int i, int n) const { return cce[(n - 1) * 5 + (i - 1)]; }
  CPPMPAS_HD double CCG(int i, int n) const { return ccg[(n - 1) * 5 + (i - 1)]; }
  static CPPMPAS_HD std::size_t i4(int i, int j, int k, int m, int ni, int nj, int nk) {
    return std::size_t(i - 1) + std::size_t(ni) * (std::size_t(j - 1) + std::size_t(nj) * (std::size_t(k - 1) +
                                                                                       std::size_t(nk) * (m - 1)));
  }
  CPPMPAS_HD std::size_t gr(int i, int j, int k, int m) const { return i4(i, j, k, m, ntb_g1, ntb_g, ntb_r1); }
  CPPMPAS_HD std::size_t sr(int i, int j, int k, int m) const { return i4(i, j, k, m, ntb_s, ntb_t, ntb_r1); }
  CPPMPAS_HD std::size_t rfz(int i, int j, int k, int m) const { return i4(i, j, k, m, ntb_r, ntb_r1, 45); }
  CPPMPAS_HD std::size_t cfz(int i, int j, int k, int m) const { return i4(i, j, k, m, ntb_c, nbc, 45); }
  CPPMPAS_HD std::size_t ia(int i, int j) const { return std::size_t(i - 1) + std::size_t(ntb_i) * (j - 1); }
};

namespace detail {
CPPMPAS_HD inline double c_min(double a, double b) { return (a < b) ? a : b; }
CPPMPAS_HD inline double c_max(double a, double b) { return (a > b) ? a : b; }
CPPMPAS_HD inline int i_min(int a, int b) { return (a < b) ? a : b; }
CPPMPAS_HD inline int i_max(int a, int b) { return (a > b) ? a : b; }
// NINT: lround@plt
CPPMPAS_HD inline int nint(double x) { return static_cast<int>(std::lround(x)); }
// REAL**INTEGER with a variable exponent: libgcc's __powidf2 (square-and-multiply, 1/y for m < 0); CPPWRF 7fac5bf
// src/core/backend.hpp:183-195 (also physics/sf_noah.hpp)
CPPMPAS_HD inline double powi(double x, int m) {
  unsigned int n = m < 0 ? 0u - unsigned(m) : unsigned(m);
  double y = (n % 2) ? x : 1.0;
  while (n >>= 1) {
    x = x * x;
    if (n % 2) y = y * x;
  }
  return m < 0 ? 1.0 / y : y;
}
// libmvec's 2-lane pow for the pairs of a vectorised loop from its first iteration i = 0, the scalar pow for an odd
// last iteration and a 1-trip loop (CPPWRF's lane(i, n), physics/mp_kessler.hpp vm::)
CPPMPAS_HD inline double vpow(double a, double b, int i, int n) {
  const bool lane = n > 1 && !((n & 1) && i == n - 1);
  return lane ? ::cppmpas::libm::mvec_pow(a, b) : ::cppmpas::libm::pow(a, b);
}
// The lookup-table index blocks of MT (e.g. MT:2034-2045): n = the exponent with 1 <= x / 10**n < 10 among
// NINT(ALOG10(x)) - 1 .. + 1 (the last tried when none fits), then INT(x / 10**n) + 10 * (n - n2) - (n - n2)
CPPMPAS_HD inline int table_index(double x, int n2) {
  const int nic = nint(::cppmpas::libm::log10(x));
  int n = nic - 1;
  for (int nn = nic - 1; nn <= nic + 1; ++nn) {
    n = nn;
    if ((x / powi(10., nn)) >= 1.0 && (x / powi(10., nn)) < 10.0) break;
  }
  return static_cast<int>(x / powi(10., n)) + 10 * (n - n2) - (n - n2);
}
}  // namespace detail

// The snow moment of order `cse` (MT:1855-1867 form): a_ * smo2**b_
CPPMPAS_HD inline double snow_moment(double tc0, double c, double smo2) {
  const double loga_ = sa[1] + sa[2] * tc0 + sa[3] * c + sa[4] * tc0 * c + sa[5] * tc0 * tc0 + sa[6] * c * c +
                       sa[7] * tc0 * tc0 * c + sa[8] * tc0 * c * c + sa[9] * tc0 * tc0 * tc0 + sa[10] * c * c * c;
  const double a_ = ::cppmpas::libm::pow(10.0, loga_);
  const double b_ = sb[1] + sb[2] * tc0 + sb[3] * c + sb[4] * tc0 * c + sb[5] * tc0 * tc0 + sb[6] * c * c +
                    sb[7] * tc0 * tc0 * c + sb[8] * tc0 * c * c + sb[9] * tc0 * tc0 * tc0 + sb[10] * c * c * c;
  return a_ * ::cppmpas::libm::pow(smo2, b_);
}

// mp_thompson (MT:1400-3642) for one column, arrays 1-based (index 0 unused, 1..kte). Nt_c / mu_c: the module
// variables mp_gt_driver sets from ntc / muc before the call (MT:1144-1147; mu_c is not read by mp_thompson).
// Returns the precipitation pptrain .. pptice (added to, as MT's INOUT arguments).
CPPMPAS_HD inline void mp_thompson(const ThompsonView& P, double Nt_c, double* qv1d, double* qc1d, double* qi1d,
                                   double* qr1d, double* qs1d, double* qg1d, double* ni1d, double* nr1d,
                                   double* nc1d, double* t1d, const double* p1d,
                                   const double* w1d, const double* dzq, double& pptrain, double& pptsnow,
                                   double& pptgraul, double& pptice, double* rainprod, double* evapprod, int kts,
                                   int kte, double dt) {
  using namespace detail;
  using ::cppmpas::libm::exp;
  using ::cppmpas::libm::log;
  using ::cppmpas::libm::log10;
  using ::cppmpas::libm::pow;
  using std::fabs;
  using std::sqrt;
  constexpr int NK = kMaxLevels + 2;
  double tten[NK], qvten[NK], qcten[NK], qiten[NK], qrten[NK], qsten[NK], qgten[NK], niten[NK], nrten[NK],
      ncten[NK];
  double prw_vcd[NK], pnc_wcd[NK], pnc_wau[NK], pnc_rcw[NK], pnc_scw[NK], pnc_gcw[NK];
  double prr_wau[NK], prr_rcw[NK], prr_rcs[NK], prr_rcg[NK], prr_sml[NK], prr_gml[NK], prr_rci[NK], prv_rev[NK],
      pnr_wau[NK], pnr_rcs[NK], pnr_rcg[NK], pnr_rci[NK], pnr_sml[NK], pnr_gml[NK], pnr_rev[NK], pnr_rcr[NK],
      pnr_rfz[NK];
  double pri_inu[NK], pni_inu[NK], pri_ihm[NK], pni_ihm[NK], pri_wfz[NK], pni_wfz[NK], pri_rfz[NK], pni_rfz[NK],
      pri_ide[NK], pni_ide[NK], pri_rci[NK], pni_rci[NK], pni_sci[NK], pni_iau[NK], pri_iha[NK], pni_iha[NK];
  double prs_iau[NK], prs_sci[NK], prs_rcs[NK], prs_scw[NK], prs_sde[NK], prs_ihm[NK], prs_ide[NK];
  double prg_scw[NK], prg_rfz[NK], prg_gde[NK], prg_gcw[NK], prg_rci[NK], prg_rcs[NK], prg_rcg[NK], prg_ihm[NK];
  constexpr double zeroD0 = 0.0;
  double temp[NK], pres[NK], qv[NK];
  double rc[NK], ri[NK], rr[NK], rs[NK], rg[NK], ni[NK], nr[NK], nc[NK];
  double rho[NK], rhof[NK], rhof2[NK];
  double qvs[NK], qvsi[NK], delQvs[NK];
  double satw[NK], sati[NK], ssatw[NK], ssati[NK];
  double diffu[NK], visco[NK], vsc2[NK], tcond[NK], lvap[NK], ocp[NK], lvt2[NK];
  double ilamr[NK], ilamg[NK], N0_r[NK], N0_g[NK];
  double mvd_r[NK], mvd_c[NK];
  double smob[NK], smo2[NK], smo1[NK], smo0[NK], smoc[NK], smod[NK], smoe[NK], smof[NK];
  double sed_r[NK], sed_s[NK], sed_g[NK], sed_i[NK], sed_n[NK], sed_c[NK];
  double vtik[NK], vtnik[NK], vtrk[NK], vtnrk[NK], vtsk[NK], vtgk[NK], vtck[NK], vtnck[NK];
  double vts_boost[NK];
  bool L_qc[NK], L_qi[NK], L_qr[NK], L_qs[NK], L_qg[NK];
  double onstep[6];
  int ksed1[6];
  double rgvm, delta_tp, orho, lfus2;
  double N0_exp, N0_min, lam_exp, lamc = 0, lamr, lamg;
  double lami, ilami, ilamc;
  double xDc = 0, Dc_b, Dc_g, xDi, xDs, xDg;
  double zeta1, zeta, taud, tau;
  double stoke_g;
  double vti, vtr, vts, vtg, vtc;
  double Mrat, ils1, ils2, t1_vts, t2_vts, t3_vts, t4_vts, C_snow;
  double tf;
  double tempc, tc0, xnc, xri, xni, xmi, oxmi, xrc, xrr, xnr;
  double xsat, rate_max, sump, ratio;
  double clap, fcd, dfcd;
  double otemp, rvs, rvs_p, rvs_pp, gamsc, alphsc, t1_evap, t1_subl;
  double r_frac, g_frac;
  double Ef_rw, Ef_sw, Ef_gw, Ef_rr;
  double dtsave, odts, odt, odzq, hgt_agl, SR;
  double xslw1, ygra1, zans1, eva_factor;
  int k, n, nstep, k_0;
  int idx_tc, idx_t, idx_s, idx_g1, idx_g, idx_r1, idx_r, idx_i1, idx_i, idx_c, idx, idx_n, idx_IN;
  bool no_micro;
  int nu_c = 0;
  (void)rgvm;
  (void)idx_c;
  (void)idx_n;
  (void)Nt_IN1;

  no_micro = true;  // MT:1519-1523
  dtsave = dt;
  odt = 1. / dt;
  odts = 1. / dtsave;

  for (k = kts; k <= kte; ++k) {  // MT:1540-1617
    tten[k] = 0.; qvten[k] = 0.; qcten[k] = 0.; qiten[k] = 0.; qrten[k] = 0.; qsten[k] = 0.; qgten[k] = 0.;
    niten[k] = 0.; nrten[k] = 0.; ncten[k] = 0.;
    prw_vcd[k] = 0.;
    pnc_wcd[k] = 0.; pnc_wau[k] = 0.; pnc_rcw[k] = 0.; pnc_scw[k] = 0.; pnc_gcw[k] = 0.;
    prv_rev[k] = 0.; prr_wau[k] = 0.; prr_rcw[k] = 0.; prr_rcs[k] = 0.; prr_rcg[k] = 0.; prr_sml[k] = 0.;
    prr_gml[k] = 0.; prr_rci[k] = 0.; pnr_wau[k] = 0.; pnr_rcs[k] = 0.; pnr_rcg[k] = 0.; pnr_rci[k] = 0.;
    pnr_sml[k] = 0.; pnr_gml[k] = 0.; pnr_rev[k] = 0.; pnr_rcr[k] = 0.; pnr_rfz[k] = 0.;
    pri_inu[k] = 0.; pni_inu[k] = 0.; pri_ihm[k] = 0.; pni_ihm[k] = 0.; pri_wfz[k] = 0.; pni_wfz[k] = 0.;
    pri_rfz[k] = 0.; pni_rfz[k] = 0.; pri_ide[k] = 0.; pni_ide[k] = 0.; pri_rci[k] = 0.; pni_rci[k] = 0.;
    pni_sci[k] = 0.; pni_iau[k] = 0.; pri_iha[k] = 0.; pni_iha[k] = 0.;
    prs_iau[k] = 0.; prs_sci[k] = 0.; prs_rcs[k] = 0.; prs_scw[k] = 0.; prs_sde[k] = 0.; prs_ihm[k] = 0.;
    prs_ide[k] = 0.;
    prg_scw[k] = 0.; prg_rfz[k] = 0.; prg_gde[k] = 0.; prg_gcw[k] = 0.; prg_rci[k] = 0.; prg_rcs[k] = 0.;
    prg_rcg[k] = 0.; prg_ihm[k] = 0.;
  }
  for (k = kts; k <= kte; ++k) {  // MT:1618-1621
    rainprod[k] = 0.;
    evapprod[k] = 0.;
  }
  for (k = kts; k <= kte; ++k) {  // MT:1633-1644
    smo0[k] = 0.; smo1[k] = 0.; smo2[k] = 0.; smob[k] = 0.; smoc[k] = 0.; smod[k] = 0.; smoe[k] = 0.;
    smof[k] = 0.; mvd_r[k] = 0.; mvd_c[k] = 0.;
  }

  // Put column of data into local arrays (MT:1649-1749)
  for (k = kts; k <= kte; ++k) {
    temp[k] = t1d[k];
    qv[k] = c_max(1.E-10, qv1d[k]);
    pres[k] = p1d[k];
    rho[k] = 0.622 * pres[k] / (R * temp[k] * (qv[k] + 0.622));

    if (qc1d[k] > R1) {
      no_micro = false;
      rc[k] = qc1d[k] * rho[k];
      nc[k] = c_max(2., c_min(nc1d[k] * rho[k], Nt_c_max));
      L_qc[k] = true;
      nu_c = i_min(15, nint(1000.E6 / nc[k]) + 2);
      lamc = pow(nc[k] * am_r * P.CCG(2, nu_c) * P.ocg1[nu_c - 1] / rc[k], P.obmr);
      xDc = (bm_r + nu_c + 1.) / lamc;
      if (xDc < D0c) {
        lamc = P.CCE(2, nu_c) / D0c;
      } else if (xDc > D0r * 2.) {
        lamc = P.CCE(2, nu_c) / (D0r * 2.);
      }
      nc[k] = c_min(Nt_c_max, P.CCG(1, nu_c) * P.ocg2[nu_c - 1] * rc[k] / am_r * pow(lamc, bm_r));
      nc[k] = Nt_c;  // .NOT. is_aerosol_aware
    } else {
      qc1d[k] = 0.0;
      nc1d[k] = 0.0;
      rc[k] = R1;
      nc[k] = 2.;
      L_qc[k] = false;
    }

    if (qi1d[k] > R1) {
      no_micro = false;
      ri[k] = qi1d[k] * rho[k];
      ni[k] = c_max(R2, ni1d[k] * rho[k]);
      L_qi[k] = true;
      if (ni[k] <= R2) {
        lami = P.cie[1] / 5.E-6;
        ni[k] = c_min(9999.E3, P.cig[0] * P.oig2 * ri[k] / am_i * pow(lami, bm_i));
      }
      lami = pow(am_i * P.cig[1] * P.oig1 * ni[k] / ri[k], P.obmi);
      ilami = 1. / lami;
      xDi = (bm_i + mu_i + 1.) * ilami;
      if (xDi < 5.E-6) {
        lami = P.cie[1] / 5.E-6;
        ni[k] = c_min(9999.E3, P.cig[0] * P.oig2 * ri[k] / am_i * pow(lami, bm_i));
      } else if (xDi > 300.E-6) {
        lami = P.cie[1] / 300.E-6;
        ni[k] = P.cig[0] * P.oig2 * ri[k] / am_i * pow(lami, bm_i);
      }
    } else {
      qi1d[k] = 0.0;
      ni1d[k] = 0.0;
      ri[k] = R1;
      ni[k] = R2;
      L_qi[k] = false;
    }

    if (qr1d[k] > R1) {
      no_micro = false;
      rr[k] = qr1d[k] * rho[k];
      nr[k] = c_max(R2, nr1d[k] * rho[k]);
      if (nr[k] <= R2) {
        mvd_r[k] = 1.0E-3;
        lamr = (3.0 + mu_r + 0.672) / mvd_r[k];
        nr[k] = P.crg[1] * P.org3 * rr[k] * pow(lamr, bm_r) / am_r;
      }
      L_qr[k] = true;
      lamr = pow(am_r * P.crg[2] * P.org2 * nr[k] / rr[k], P.obmr);
      mvd_r[k] = (3.0 + mu_r + 0.672) / lamr;
      if (mvd_r[k] > 2.5E-3) {
        mvd_r[k] = 2.5E-3;
        lamr = (3.0 + mu_r + 0.672) / mvd_r[k];
        nr[k] = P.crg[1] * P.org3 * rr[k] * pow(lamr, bm_r) / am_r;
      } else if (mvd_r[k] < D0r * 0.75) {
        mvd_r[k] = D0r * 0.75;
        lamr = (3.0 + mu_r + 0.672) / mvd_r[k];
        nr[k] = P.crg[1] * P.org3 * rr[k] * pow(lamr, bm_r) / am_r;
      }
    } else {
      qr1d[k] = 0.0;
      nr1d[k] = 0.0;
      rr[k] = R1;
      nr[k] = R2;
      L_qr[k] = false;
    }
    if (qs1d[k] > R1) {
      no_micro = false;
      rs[k] = qs1d[k] * rho[k];
      L_qs[k] = true;
    } else {
      qs1d[k] = 0.0;
      rs[k] = R1;
      L_qs[k] = false;
    }
    if (qg1d[k] > R1) {
      no_micro = false;
      rg[k] = qg1d[k] * rho[k];
      L_qg[k] = true;
    } else {
      qg1d[k] = 0.0;
      rg[k] = R1;
      L_qg[k] = false;
    }
  }

  // Derive various thermodynamic variables frequently used (MT:1769-1811)
  for (k = kts; k <= kte; ++k) {
    tempc = temp[k] - 273.15;
    rhof[k] = sqrt(rho_not / rho[k]);
    rhof2[k] = sqrt(rhof[k]);
    qvs[k] = rslf(pres[k], temp[k]);
    delQvs[k] = c_max(0.0, rslf(pres[k], 273.15) - qv[k]);
    if (tempc <= 0.0) {
      qvsi[k] = rsif(pres[k], temp[k]);
    } else {
      qvsi[k] = qvs[k];
    }
    satw[k] = qv[k] / qvs[k];
    sati[k] = qv[k] / qvsi[k];
    ssatw[k] = satw[k] - 1.;
    ssati[k] = sati[k] - 1.;
    if (fabs(ssatw[k]) < eps) ssatw[k] = 0.0;
    if (fabs(ssati[k]) < eps) ssati[k] = 0.0;
    if (no_micro && ssati[k] > 0.0) no_micro = false;
    diffu[k] = 2.11E-5 * pow(temp[k] / 273.15, 1.94) * (101325. / pres[k]);
    if (tempc >= 0.0) {
      visco[k] = (1.718 + 0.0049 * tempc) * 1.0E-5;
    } else {
      visco[k] = (1.718 + 0.0049 * tempc - 1.2E-5 * tempc * tempc) * 1.0E-5;
    }
    ocp[k] = 1. / (Cp * (1. + 0.887 * qv[k]));
    vsc2[k] = sqrt(rho[k] / visco[k]);
    lvap[k] = lvap0 + (2106.0 - 4218.0) * tempc;
    tcond[k] = (5.69 + 0.0168 * tempc) * 1.0E-5 * 418.936;
  }

  if (no_micro) return;  // MT:1818

  // Snow moments (MT:1823-1908); bm_s = 2: smo2 = smob
  for (k = kts; k <= kte; ++k) {
    if (!L_qs[k]) continue;
    tc0 = c_min(-0.1, temp[k] - 273.15);
    smob[k] = rs[k] * P.oams;
    smo2[k] = smob[k];
    {  // 0th moment (MT:1849-1852)
      const double loga_ = sa[1] + sa[2] * tc0 + sa[5] * tc0 * tc0 + sa[9] * tc0 * tc0 * tc0;
      const double a_ = pow(10.0, loga_);
      const double b_ = sb[1] + sb[2] * tc0 + sb[5] * tc0 * tc0 + sb[9] * tc0 * tc0 * tc0;
      smo0[k] = a_ * pow(smo2[k], b_);
    }
    {  // 1st moment (MT:1855-1865)
      const double loga_ = sa[1] + sa[2] * tc0 + sa[3] + sa[4] * tc0 + sa[5] * tc0 * tc0 + sa[6] +
                           sa[7] * tc0 * tc0 + sa[8] * tc0 + sa[9] * tc0 * tc0 * tc0 + sa[10];
      const double a_ = pow(10.0, loga_);
      const double b_ = sb[1] + sb[2] * tc0 + sb[3] + sb[4] * tc0 + sb[5] * tc0 * tc0 + sb[6] + sb[7] * tc0 * tc0 +
                        sb[8] * tc0 + sb[9] * tc0 * tc0 * tc0 + sb[10];
      smo1[k] = a_ * pow(smo2[k], b_);
    }
    smoc[k] = snow_moment(tc0, P.cse[0], smo2[k]);   // bm_s+1 (MT:1868-1878)
    smoe[k] = snow_moment(tc0, P.cse[12], smo2[k]);  // bv_s+2 (MT:1881-1891)
    smof[k] = snow_moment(tc0, P.cse[15], smo2[k]);  // 1+(bv_s+1)/2 (MT:1894-1904)
  }

  // Graupel y-intercept and slope (MT:1913-1933)
  N0_min = gonv_max;
  k_0 = kts;
  for (k = kte; k >= kts; --k)
    if (temp[k] >= 270.65) k_0 = i_max(k_0, k);
  for (k = kte; k >= kts; --k) {
    if (k > k_0 && L_qr[k] && mvd_r[k] > 100.E-6) {
      xslw1 = 4.01 + log10(mvd_r[k]);
    } else {
      xslw1 = 0.01;
    }
    ygra1 = 4.31 + log10(c_max(5.E-5, rg[k]));
    zans1 = 3.1 + (100. / (300. * xslw1 * ygra1 / (10. / xslw1 + 1. + 0.25 * ygra1) + 30. + 10. * ygra1));
    N0_exp = pow(10., zans1);
    N0_exp = c_max(gonv_min, c_min(N0_exp, gonv_max));
    N0_min = c_min(N0_exp, N0_min);
    N0_exp = N0_min;
    lam_exp = pow(N0_exp * am_g * P.cgg[0] / rg[k], P.oge1);
    lamg = lam_exp * pow(P.cgg[2] * P.ogg2 * P.ogg1, P.obmg);
    ilamg[k] = 1. / lamg;
    N0_g[k] = N0_exp / (P.cgg[1] * lam_exp) * pow(lamg, P.cge[1]);
  }

  // Rain y-intercept and slope (MT:1939-1945): vectorised from kte, _ZGVbN2vv_pow (MT:1941, 1944)
  for (k = kte; k >= kts; --k) {
    const int li = kte - k, ln = kte - kts + 1;
    lamr = vpow(am_r * P.crg[2] * P.org2 * nr[k] / rr[k], P.obmr, li, ln);
    ilamr[k] = 1. / lamr;
    mvd_r[k] = (3.0 + mu_r + 0.672) / lamr;
    N0_r[k] = nr[k] * P.org2 * vpow(lamr, P.cre[1], li, ln);
  }

  // Warm-rain process terms (MT:1950-2016)
  for (k = kts; k <= kte; ++k) {
    if (L_qr[k] && mvd_r[k] > D0r) {
      Ef_rr = 1.0 - exp(2300.0 * (mvd_r[k] - 1950.0E-6));
      pnr_rcr[k] = Ef_rr * 2.0 * nr[k] * rr[k];
    }

    mvd_c[k] = D0c;
    if (L_qc[k]) {
      nu_c = i_min(15, nint(1000.E6 / nc[k]) + 2);
      xDc = c_max(D0c * 1.E6, (pow(rc[k] / (am_r * nc[k]), P.obmr)) * 1.E6);
      lamc = pow(nc[k] * am_r * P.CCG(2, nu_c) * P.ocg1[nu_c - 1] / rc[k], P.obmr);
      mvd_c[k] = (3.0 + nu_c + 0.672) / lamc;
    }

    if (rc[k] > 0.01e-3) {
      Dc_g = (pow(P.CCG(3, nu_c) * P.ocg2[nu_c - 1], P.obmr) / lamc) * 1.E6;
      Dc_b = pow(xDc * xDc * xDc * Dc_g * Dc_g * Dc_g - xDc * xDc * xDc * xDc * xDc * xDc, 1. / 6.);
      zeta1 = 0.5 * ((6.25E-6 * xDc * Dc_b * Dc_b * Dc_b - 0.4) + fabs(6.25E-6 * xDc * Dc_b * Dc_b * Dc_b - 0.4));
      zeta = 0.027 * rc[k] * zeta1;
      taud = 0.5 * ((0.5 * Dc_b - 7.5) + fabs(0.5 * Dc_b - 7.5)) + R1;
      tau = 3.72 / (rc[k] * taud);
      prr_wau[k] = zeta / tau;
      prr_wau[k] = c_min(rc[k] * odts, prr_wau[k]);
      pnr_wau[k] = prr_wau[k] / (am_r * nu_c * 200. * D0r * D0r * D0r);
      pnc_wau[k] = c_min(nc[k] * odts, prr_wau[k] / (am_r * mvd_c[k] * mvd_c[k] * mvd_c[k]));
    }

    if (L_qr[k] && mvd_r[k] > D0r && mvd_c[k] > D0c) {
      lamr = 1. / ilamr[k];
      idx = 1 + static_cast<int>(nbr * log(mvd_r[k] / P.Dr[0]) / log(P.Dr[nbr - 1] / P.Dr[0]));
      idx = i_min(idx, nbr);
      Ef_rw = P.t_Efrw[(idx - 1) + std::size_t(nbr) * (static_cast<int>(mvd_c[k] * 1.E6) - 1)];
      prr_rcw[k] = rhof[k] * P.t1_qr_qc * Ef_rw * rc[k] * N0_r[k] * (pow(lamr + fv_r, -P.cre[8]));
      prr_rcw[k] = c_min(rc[k] * odts, prr_rcw[k]);
      pnc_rcw[k] = rhof[k] * P.t1_qr_qc * Ef_rw * nc[k] * N0_r[k] * (pow(lamr + fv_r, -P.cre[8]));
      pnc_rcw[k] = c_min(nc[k] * odts, pnc_rcw[k]);
    }
    // MT:2003-2014: rain collecting aerosols (pna_rca, pnd_rcd): aerosol-aware only
  }

  // Frozen species process terms (MT:2021-2551)
  for (k = kts; k <= kte; ++k) {
    vts_boost[k] = 1.0;
    xDs = 0.0;
    if (L_qs[k]) xDs = smoc[k] / smob[k];

    tempc = temp[k] - 273.15;
    idx_tc = i_max(1, i_min(nint(-tempc), 45));
    idx_t = static_cast<int>((tempc - 2.5) / 5.) - 1;
    idx_t = i_max(1, -idx_t);
    idx_t = i_min(idx_t, ntb_t);

    if (rc[k] > r_c1) {
      idx_c = i_max(1, i_min(table_index(rc[k], P.nic2), ntb_c));
    } else {
      idx_c = 1;
    }

    idx_n = nint(1.0 + double(nbc) * log(nc[k] / P.t_Nc[0]) / P.nic1);
    idx_n = i_max(1, i_min(idx_n, nbc));

    if (ri[k] > r_i1) {
      idx_i = i_max(1, i_min(table_index(ri[k], P.nii2), ntb_i));
    } else {
      idx_i = 1;
    }

    if (ni[k] > Nt_i1) {
      idx_i1 = i_max(1, i_min(table_index(ni[k], P.nii3), ntb_i1));
    } else {
      idx_i1 = 1;
    }

    if (rr[k] > r_r1) {
      idx_r = i_max(1, i_min(table_index(rr[k], P.nir2), ntb_r));
      lamr = 1. / ilamr[k];
      lam_exp = lamr * pow(P.crg[2] * P.org2 * P.org1, bm_r);
      N0_exp = P.org1 * rr[k] / am_r * pow(lam_exp, P.cre[0]);
      idx_r1 = i_max(1, i_min(table_index(N0_exp, P.nir3), ntb_r1));
    } else {
      idx_r = 1;
      idx_r1 = ntb_r1;
    }

    if (rs[k] > r_s1) {
      idx_s = i_max(1, i_min(table_index(rs[k], P.nis2), ntb_s));
    } else {
      idx_s = 1;
    }

    if (rg[k] > r_g1) {
      idx_g = i_max(1, i_min(table_index(rg[k], P.nig2), ntb_g));
      lamg = 1. / ilamg[k];
      lam_exp = lamg * pow(P.cgg[2] * P.ogg2 * P.ogg1, bm_g);
      N0_exp = P.ogg1 * rg[k] / am_g * pow(lam_exp, P.cge[0]);
      idx_g1 = i_max(1, i_min(table_index(N0_exp, P.nig3), ntb_g1));
    } else {
      idx_g = 1;
      idx_g1 = ntb_g1;
    }

    // Deposition/sublimation prefactor (MT:2156-2173)
    otemp = 1. / temp[k];
    rvs = rho[k] * qvsi[k];
    rvs_p = rvs * otemp * (lsub * otemp * oRv - 1.);
    rvs_pp = rvs * (otemp * (lsub * otemp * oRv - 1.) * otemp * (lsub * otemp * oRv - 1.) +
                    (-2. * lsub * otemp * otemp * otemp * oRv) + otemp * otemp);
    gamsc = lsub * diffu[k] / tcond[k] * rvs_p;
    alphsc = 0.5 * (gamsc / (1. + gamsc)) * (gamsc / (1. + gamsc)) * rvs_pp / rvs_p * rvs / rvs_p;
    alphsc = c_max(1.E-9, alphsc);
    xsat = ssati[k];
    if (fabs(xsat) < 1.E-9) xsat = 0.;
    t1_subl = 4. * PI *
              (1.0 - alphsc * xsat + 2. * alphsc * alphsc * xsat * xsat -
               5. * alphsc * alphsc * alphsc * xsat * xsat * xsat) /
              (1. + gamsc);

    // Snow collecting cloud water (MT:2176-2186)
    if (L_qc[k] && mvd_c[k] > D0c) {
      if (xDs > D0s) {
        idx = 1 + static_cast<int>(nbs * log(xDs / P.Ds[0]) / log(P.Ds[nbs - 1] / P.Ds[0]));
        idx = i_min(idx, nbs);
        Ef_sw = P.t_Efsw[(idx - 1) + std::size_t(nbs) * (static_cast<int>(mvd_c[k] * 1.E6) - 1)];
        prs_scw[k] = rhof[k] * P.t1_qs_qc * Ef_sw * rc[k] * smoe[k];
        prs_scw[k] = c_min(rc[k] * odts, prs_scw[k]);
        pnc_scw[k] = rhof[k] * P.t1_qs_qc * Ef_sw * nc[k] * smoe[k];
        pnc_scw[k] = c_min(nc[k] * odts, pnc_scw[k]);
      }

      // Graupel collecting cloud water (MT:2189-2209)
      if (rg[k] >= r_g1 && mvd_c[k] > D0c) {
        xDg = (bm_g + mu_g + 1.) * ilamg[k];
        vtg = rhof[k] * av_g * P.cgg[5] * P.ogg3 * pow(ilamg[k], bv_g);
        stoke_g = mvd_c[k] * mvd_c[k] * vtg * rho_w / (9. * visco[k] * xDg);
        if (xDg > D0g) {
          if (stoke_g >= 0.4 && stoke_g <= 10.) {
            Ef_gw = 0.55 * log10(2.51 * stoke_g);
          } else if (stoke_g < 0.4) {
            Ef_gw = 0.0;
          } else if (stoke_g > 10) {
            Ef_gw = 0.77;
          }
          prg_gcw[k] = rhof[k] * P.t1_qg_qc * Ef_gw * rc[k] * N0_g[k] * pow(ilamg[k], P.cge[8]);
          pnc_gcw[k] = rhof[k] * P.t1_qg_qc * Ef_gw * nc[k] * N0_g[k] * pow(ilamg[k], P.cge[8]);
          pnc_gcw[k] = c_min(nc[k] * odts, pnc_gcw[k]);
        }
      }
    }
    // MT:2212-2233: snow and graupel collecting aerosols: aerosol-aware only

    // Rain collecting snow (MT:2238-2270)
    const std::size_t is = P.sr(idx_s, idx_t, idx_r1, idx_r);
    if (rr[k] >= r_r1) {
      if (rs[k] >= r_s1) {
        if (temp[k] < T_0) {
          prr_rcs[k] = -(P.tmr_racs2[is] + P.tcr_sacr2[is] + P.tmr_racs1[is] + P.tcr_sacr1[is]);
          prs_rcs[k] = P.tmr_racs2[is] + P.tcr_sacr2[is] - P.tcs_racs1[is] - P.tms_sacr1[is];
          prg_rcs[k] = P.tmr_racs1[is] + P.tcr_sacr1[is] + P.tcs_racs1[is] + P.tms_sacr1[is];
          prr_rcs[k] = c_max(-rr[k] * odts, prr_rcs[k]);
          prs_rcs[k] = c_max(-rs[k] * odts, prs_rcs[k]);
          prg_rcs[k] = c_min((rr[k] + rs[k]) * odts, prg_rcs[k]);
          pnr_rcs[k] = P.tnr_racs1[is] + P.tnr_racs2[is] + P.tnr_sacr1[is] + P.tnr_sacr2[is];
          pnr_rcs[k] = c_min(nr[k] * odts, pnr_rcs[k]);
        } else {
          prs_rcs[k] = -P.tcs_racs1[is] - P.tms_sacr1[is] + P.tmr_racs2[is] + P.tcr_sacr2[is];
          prs_rcs[k] = c_max(-rs[k] * odts, prs_rcs[k]);
          prr_rcs[k] = -prs_rcs[k];
        }
      }

      // Rain collecting graupel (MT:2275-2298)
      if (rg[k] >= r_g1) {
        const std::size_t ig = P.gr(idx_g1, idx_g, idx_r1, idx_r);
        if (temp[k] < T_0) {
          prg_rcg[k] = P.tmr_racg[ig] + P.tcr_gacr[ig];
          prg_rcg[k] = c_min(rr[k] * odts, prg_rcg[k]);
          prr_rcg[k] = -prg_rcg[k];
          pnr_rcg[k] = P.tnr_racg[ig] + P.tnr_gacr[ig];
          pnr_rcg[k] = c_min(nr[k] * odts, pnr_rcg[k]);
        } else {
          prr_rcg[k] = P.tcg_racg[ig];
          prr_rcg[k] = c_min(rg[k] * odts, prr_rcg[k]);
          prg_rcg[k] = -prr_rcg[k];
          pnr_rcg[k] = -5.0 * P.tnr_gacr[ig];
        }
      }
    }

    // Processes below 0C (MT:2306-2507)
    if (temp[k] < T_0) {
      vts_boost[k] = 1.0;
      rate_max = (qv[k] - qvsi[k]) * rho[k] * odts * 0.999;

      // xni = 1.0 * 1000. (MT:2329, dustyIce without aerosols): the IN index of 1000 per m3
      xni = 1.0 * 1000.;
      if (xni > Nt_IN1) {
        idx_IN = i_max(1, i_min(table_index(xni, P.niIN2), ntb_IN));
      } else {
        idx_IN = 1;
      }

      // Freezing of water drops (MT:2348-2357)
      if (rr[k] > r_r1) {
        const std::size_t ir = P.rfz(idx_r, idx_r1, idx_tc, idx_IN);
        prg_rfz[k] = P.tpg_qrfz[ir] * odts;
        pri_rfz[k] = P.tpi_qrfz[ir] * odts;
        pni_rfz[k] = P.tni_qrfz[ir] * odts;
        pnr_rfz[k] = P.tnr_qrfz[ir] * odts;
        pnr_rfz[k] = c_min(nr[k] * odts, pnr_rfz[k]);
      } else if (rr[k] > R1 && temp[k] < HGFR) {
        pri_rfz[k] = rr[k] * odts;
        pni_rfz[k] = nr[k] * odts;
      }

      if (rc[k] > r_c1) {  // MT:2359-2368
        const std::size_t ic = P.cfz(idx_c, idx_n, idx_tc, idx_IN);
        pri_wfz[k] = P.tpi_qcfz[ic] * odts;
        pri_wfz[k] = c_min(rc[k] * odts, pri_wfz[k]);
        pni_wfz[k] = P.tni_qcfz[ic] * odts;
        pni_wfz[k] = c_min(c_min(nc[k] * odts, pri_wfz[k] / (2. * xm0i)), pni_wfz[k]);
      } else if (rc[k] > R1 && temp[k] < HGFR) {
        pri_wfz[k] = rc[k] * odts;
        pni_wfz[k] = nc[k] * odts;
      }

      // Deposition nucleation (MT:2372-2383)
      if ((ssati[k] >= 0.25) || (ssatw[k] > eps && temp[k] < 253.15)) {
        xnc = c_min(250.E3, TNO * exp(ATO * (T_0 - temp[k])));
        xni = ni[k] + (pni_rfz[k] + pni_wfz[k]) * dtsave;
        pni_inu[k] = 0.5 * (xnc - xni + fabs(xnc - xni)) * odts;
        pri_inu[k] = c_min(rate_max, xm0i * pni_inu[k]);
        pni_inu[k] = pri_inu[k] / xm0i;
      }
      // MT:2386-2393: Koop freezing of aqueous aerosols: aerosol-aware only

      // Deposition/sublimation of cloud ice (MT:2398-2433)
      if (L_qi[k]) {
        lami = pow(am_i * P.cig[1] * P.oig1 * ni[k] / ri[k], P.obmi);
        ilami = 1. / lami;
        xDi = c_max(P.D0i, (bm_i + mu_i + 1.) * ilami);
        xmi = am_i * pow(xDi, bm_i);
        oxmi = 1. / xmi;
        pri_ide[k] = C_cube * t1_subl * diffu[k] * ssati[k] * rvs * P.oig1 * P.cig[4] * ni[k] * ilami;

        if (pri_ide[k] < 0.0) {
          pri_ide[k] = c_max(c_max(-ri[k] * odts, pri_ide[k]), rate_max);
          pni_ide[k] = pri_ide[k] * oxmi;
          pni_ide[k] = c_max(-ni[k] * odts, pni_ide[k]);
        } else {
          pri_ide[k] = c_min(pri_ide[k], rate_max);
          prs_ide[k] = (1.0 - P.tpi_ide[P.ia(idx_i, idx_i1)]) * pri_ide[k];
          pri_ide[k] = P.tpi_ide[P.ia(idx_i, idx_i1)] * pri_ide[k];
        }

        if ((idx_i == ntb_i) || (xDi > 5.0 * D0s)) {
          prs_iau[k] = ri[k] * .99 * odts;
          pni_iau[k] = ni[k] * .95 * odts;
        } else if (xDi < 0.1 * D0s) {
          prs_iau[k] = 0.;
          pni_iau[k] = 0.;
        } else {
          prs_iau[k] = P.tps_iaus[P.ia(idx_i, idx_i1)] * odts;
          prs_iau[k] = c_min(ri[k] * .99 * odts, prs_iau[k]);
          pni_iau[k] = P.tni_iaus[P.ia(idx_i, idx_i1)] * odts;
          pni_iau[k] = c_min(ni[k] * .95 * odts, pni_iau[k]);
        }
      }

      // Deposition/sublimation of snow / graupel (MT:2437-2460)
      if (L_qs[k]) {
        C_snow = C_sqrd + (tempc + 1.5) * (C_cube - C_sqrd) / (-30. + 1.5);
        C_snow = c_max(C_sqrd, c_min(C_snow, C_cube));
        prs_sde[k] = C_snow * t1_subl * diffu[k] * ssati[k] * rvs *
                     (P.t1_qs_sd * smo1[k] + P.t2_qs_sd * rhof2[k] * vsc2[k] * smof[k]);
        if (prs_sde[k] < 0.) {
          prs_sde[k] = c_max(c_max(-rs[k] * odts, prs_sde[k]), rate_max);
        } else {
          prs_sde[k] = c_min(prs_sde[k], rate_max);
        }
      }

      if (L_qg[k] && ssati[k] < -eps) {
        prg_gde[k] = C_cube * t1_subl * diffu[k] * ssati[k] * rvs * N0_g[k] *
                     (P.t1_qg_sd * pow(ilamg[k], P.cge[9]) + P.t2_qg_sd * vsc2[k] * rhof2[k] * pow(ilamg[k], P.cge[10]));
        if (prg_gde[k] < 0.) {
          prg_gde[k] = c_max(c_max(-rg[k] * odts, prg_gde[k]), rate_max);
        } else {
          prg_gde[k] = c_min(prg_gde[k], rate_max);
        }
      }

      // Snow collecting cloud ice, rain collecting cloud ice (MT:2463-2482)
      if (L_qi[k]) {
        lami = pow(am_i * P.cig[1] * P.oig1 * ni[k] / ri[k], P.obmi);
        ilami = 1. / lami;
        xDi = c_max(P.D0i, (bm_i + mu_i + 1.) * ilami);
        xmi = am_i * pow(xDi, bm_i);
        oxmi = 1. / xmi;
        if (rs[k] >= r_s1) {
          prs_sci[k] = P.t1_qs_qi * rhof[k] * Ef_si * ri[k] * smoe[k];
          pni_sci[k] = prs_sci[k] * oxmi;
        }

        if (rr[k] >= r_r1 && mvd_r[k] > 4. * xDi) {
          lamr = 1. / ilamr[k];
          pri_rci[k] = rhof[k] * P.t1_qr_qi * Ef_ri * ri[k] * N0_r[k] * (pow(lamr + fv_r, -P.cre[8]));
          pnr_rci[k] = rhof[k] * P.t1_qr_qi * Ef_ri * ni[k] * N0_r[k] * (pow(lamr + fv_r, -P.cre[8]));
          pni_rci[k] = pri_rci[k] * oxmi;
          prr_rci[k] = rhof[k] * P.t2_qr_qi * Ef_ri * ni[k] * N0_r[k] * (pow(lamr + fv_r, -P.cre[7]));
          prr_rci[k] = c_min(rr[k] * odts, prr_rci[k]);
          prg_rci[k] = pri_rci[k] + prr_rci[k];
        }
      }

      // Hallett-Mossop (MT:2485-2498)
      if (prg_gcw[k] > eps && tempc > -8.0) {
        tf = 0.;
        if (tempc >= -5.0 && tempc < -3.0) {
          tf = 0.5 * (-3.0 - tempc);
        } else if (tempc > -8.0 && tempc < -5.0) {
          tf = 0.33333333 * (8.0 + tempc);
        }
        pni_ihm[k] = 3.5E8 * tf * prg_gcw[k];
        pri_ihm[k] = xm0i * pni_ihm[k];
        prs_ihm[k] = prs_scw[k] / (prs_scw[k] + prg_gcw[k]) * pri_ihm[k];
        prg_ihm[k] = prg_gcw[k] / (prs_scw[k] + prg_gcw[k]) * pri_ihm[k];
      }

      // Rimed snow into graupel (MT:2504-2511)
      if (prs_scw[k] > 2.0 * prs_sde[k] && prs_sde[k] > eps) {
        r_frac = c_min(30.0, prs_scw[k] / prs_sde[k]);
        g_frac = c_min(0.95, 0.15 + (r_frac - 2.) * .028);
        vts_boost[k] = c_min(1.5, 1.1 + (r_frac - 2.) * .014);
        prg_scw[k] = g_frac * prs_scw[k];
        prs_scw[k] = (1. - g_frac) * prs_scw[k];
      }

    } else {
      // Melting and sublimation of snow and graupel (MT:2517-2560)
      if (L_qs[k]) {
        prr_sml[k] = (tempc * tcond[k] - lvap0 * diffu[k] * delQvs[k]) *
                     (P.t1_qs_me * smo1[k] + P.t2_qs_me * rhof2[k] * vsc2[k] * smof[k]);
        if (prr_sml[k] > 0.) {
          prr_sml[k] = prr_sml[k] + 4218. * olfus * tempc * (prr_rcs[k] + prs_scw[k]);
        }
        prr_sml[k] = c_min(rs[k] * odts, c_max(0.0, prr_sml[k]));
        pnr_sml[k] = smo0[k] / rs[k] * prr_sml[k] * pow(10.0, -0.25 * tempc);
        pnr_sml[k] = c_min(smo0[k] * odts, pnr_sml[k]);

        if (ssati[k] < 0.) {
          prs_sde[k] = C_cube * t1_subl * diffu[k] * ssati[k] * rvs *
                       (P.t1_qs_sd * smo1[k] + P.t2_qs_sd * rhof2[k] * vsc2[k] * smof[k]);
          prs_sde[k] = c_max(-rs[k] * odts, prs_sde[k]);
        }
      }

      if (L_qg[k]) {
        prr_gml[k] = (tempc * tcond[k] - lvap0 * diffu[k] * delQvs[k]) * N0_g[k] *
                     (P.t1_qg_me * pow(ilamg[k], P.cge[9]) + P.t2_qg_me * rhof2[k] * vsc2[k] * pow(ilamg[k], P.cge[10]));
        prr_gml[k] = c_min(rg[k] * odts, c_max(0.0, prr_gml[k]));
        pnr_gml[k] = N0_g[k] * P.cgg[1] * pow(ilamg[k], P.cge[1]) / rg[k] * prr_gml[k] * pow(10.0, -0.5 * tempc);

        if (ssati[k] < 0.) {
          prg_gde[k] = C_cube * t1_subl * diffu[k] * ssati[k] * rvs * N0_g[k] *
                       (P.t1_qg_sd * pow(ilamg[k], P.cge[9]) + P.t2_qg_sd * vsc2[k] * rhof2[k] * pow(ilamg[k], P.cge[10]));
          prg_gde[k] = c_max(-rg[k] * odts, prg_gde[k]);
        }
      }

      if (dt > 120.) {  // MT:2555-2559
        prr_rcw[k] = prr_rcw[k] + prs_scw[k] + prg_gcw[k];
        prs_scw[k] = 0.;
        prg_gcw[k] = 0.;
      }
    }
  }

  // Ensure we do not deplete more hydrometeor species than exists (MT:2569-2662)
  for (k = kts; k <= kte; ++k) {
    sump = pri_inu[k] + pri_ide[k] + prs_ide[k] + prs_sde[k] + prg_gde[k] + pri_iha[k];
    rate_max = (qv[k] - qvsi[k]) * rho[k] * odts * 0.999;
    if ((sump > eps && sump > rate_max) || (sump < -eps && sump < rate_max)) {
      ratio = rate_max / sump;
      pri_inu[k] = pri_inu[k] * ratio;
      pri_ide[k] = pri_ide[k] * ratio;
      pni_ide[k] = pni_ide[k] * ratio;
      prs_ide[k] = prs_ide[k] * ratio;
      prs_sde[k] = prs_sde[k] * ratio;
      prg_gde[k] = prg_gde[k] * ratio;
      pri_iha[k] = pri_iha[k] * ratio;
    }

    sump = -prr_wau[k] - pri_wfz[k] - prr_rcw[k] - prs_scw[k] - prg_scw[k] - prg_gcw[k];
    rate_max = -rc[k] * odts;
    if (sump < rate_max && L_qc[k]) {
      ratio = rate_max / sump;
      prr_wau[k] = prr_wau[k] * ratio;
      pri_wfz[k] = pri_wfz[k] * ratio;
      prr_rcw[k] = prr_rcw[k] * ratio;
      prs_scw[k] = prs_scw[k] * ratio;
      prg_scw[k] = prg_scw[k] * ratio;
      prg_gcw[k] = prg_gcw[k] * ratio;
    }

    sump = pri_ide[k] - prs_iau[k] - prs_sci[k] - pri_rci[k];
    rate_max = -ri[k] * odts;
    if (sump < rate_max && L_qi[k]) {
      ratio = rate_max / sump;
      pri_ide[k] = pri_ide[k] * ratio;
      prs_iau[k] = prs_iau[k] * ratio;
      prs_sci[k] = prs_sci[k] * ratio;
      pri_rci[k] = pri_rci[k] * ratio;
    }

    sump = -prg_rfz[k] - pri_rfz[k] - prr_rci[k] + prr_rcs[k] + prr_rcg[k];
    rate_max = -rr[k] * odts;
    if (sump < rate_max && L_qr[k]) {
      ratio = rate_max / sump;
      prg_rfz[k] = prg_rfz[k] * ratio;
      pri_rfz[k] = pri_rfz[k] * ratio;
      prr_rci[k] = prr_rci[k] * ratio;
      prr_rcs[k] = prr_rcs[k] * ratio;
      prr_rcg[k] = prr_rcg[k] * ratio;
    }

    sump = prs_sde[k] - prs_ihm[k] - prr_sml[k] + prs_rcs[k];
    rate_max = -rs[k] * odts;
    if (sump < rate_max && L_qs[k]) {
      ratio = rate_max / sump;
      prs_sde[k] = prs_sde[k] * ratio;
      prs_ihm[k] = prs_ihm[k] * ratio;
      prr_sml[k] = prr_sml[k] * ratio;
      prs_rcs[k] = prs_rcs[k] * ratio;
    }

    sump = prg_gde[k] - prg_ihm[k] - prr_gml[k] + prg_rcg[k];
    rate_max = -rg[k] * odts;
    if (sump < rate_max && L_qg[k]) {
      ratio = rate_max / sump;
      prg_gde[k] = prg_gde[k] * ratio;
      prg_ihm[k] = prg_ihm[k] * ratio;
      prr_gml[k] = prr_gml[k] * ratio;
      prg_rcg[k] = prg_rcg[k] * ratio;
    }

    pri_ihm[k] = prs_ihm[k] + prg_ihm[k];
    ratio = c_min(fabs(prr_rcg[k]), fabs(prg_rcg[k]));
    prr_rcg[k] = ratio * std::copysign(1.0, prr_rcg[k]);
    prg_rcg[k] = -prr_rcg[k];
    if (temp[k] > T_0) {
      ratio = c_min(fabs(prr_rcs[k]), fabs(prs_rcs[k]));
      prr_rcs[k] = ratio * std::copysign(1.0, prr_rcs[k]);
      prs_rcs[k] = -prr_rcs[k];
    }
  }

  // Tendencies (MT:2668-2841)
  for (k = kts; k <= kte; ++k) {
    orho = 1. / rho[k];
    lfus2 = lsub - lvap[k];

    qvten[k] = qvten[k] + (-pri_inu[k] - pri_iha[k] - pri_ide[k] - prs_ide[k] - prs_sde[k] - prg_gde[k]) * orho;

    qcten[k] = qcten[k] + (-prr_wau[k] - pri_wfz[k] - prr_rcw[k] - prs_scw[k] - prg_scw[k] - prg_gcw[k]) * orho;

    ncten[k] = ncten[k] + (-pnc_wau[k] - pnc_rcw[k] - pni_wfz[k] - pnc_scw[k] - pnc_gcw[k]) * orho;

    xrc = c_max(R1, (qc1d[k] + qcten[k] * dtsave) * rho[k]);
    xnc = c_max(2., (nc1d[k] + ncten[k] * dtsave) * rho[k]);
    if (xrc > R1) {
      nu_c = i_min(15, nint(1000.E6 / xnc) + 2);
      lamc = pow(xnc * am_r * P.CCG(2, nu_c) * P.ocg1[nu_c - 1] / rc[k], P.obmr);
      xDc = (bm_r + nu_c + 1.) / lamc;
      if (xDc < D0c) {
        lamc = P.CCE(2, nu_c) / D0c;
        xnc = P.CCG(1, nu_c) * P.ocg2[nu_c - 1] * xrc / am_r * pow(lamc, bm_r);
        ncten[k] = (xnc - nc1d[k] * rho[k]) * odts * orho;
      } else if (xDc > D0r * 2.) {
        lamc = P.CCE(2, nu_c) / (D0r * 2.);
        xnc = P.CCG(1, nu_c) * P.ocg2[nu_c - 1] * xrc / am_r * pow(lamc, bm_r);
        ncten[k] = (xnc - nc1d[k] * rho[k]) * odts * orho;
      }
    } else {
      ncten[k] = -nc1d[k] * odts;
    }
    xnc = c_max(0., (nc1d[k] + ncten[k] * dtsave) * rho[k]);
    if (xnc > Nt_c_max) ncten[k] = (Nt_c_max - nc1d[k] * rho[k]) * odts * orho;

    qiten[k] = qiten[k] + (pri_inu[k] + pri_iha[k] + pri_ihm[k] + pri_wfz[k] + pri_rfz[k] + pri_ide[k] -
                           prs_iau[k] - prs_sci[k] - pri_rci[k]) *
                              orho;

    niten[k] = niten[k] + (pni_inu[k] + pni_iha[k] + pni_ihm[k] + pni_wfz[k] + pni_rfz[k] + pni_ide[k] -
                           pni_iau[k] - pni_sci[k] - pni_rci[k]) *
                              orho;

    xri = c_max(R1, (qi1d[k] + qiten[k] * dtsave) * rho[k]);
    xni = c_max(R2, (ni1d[k] + niten[k] * dtsave) * rho[k]);
    if (xri > R1) {
      lami = pow(am_i * P.cig[1] * P.oig1 * xni / xri, P.obmi);
      ilami = 1. / lami;
      xDi = (bm_i + mu_i + 1.) * ilami;
      if (xDi < 5.E-6) {
        lami = P.cie[1] / 5.E-6;
        xni = c_min(9999.E3, P.cig[0] * P.oig2 * xri / am_i * pow(lami, bm_i));
        niten[k] = (xni - ni1d[k] * rho[k]) * odts * orho;
      } else if (xDi > 300.E-6) {
        lami = P.cie[1] / 300.E-6;
        xni = P.cig[0] * P.oig2 * xri / am_i * pow(lami, bm_i);
        niten[k] = (xni - ni1d[k] * rho[k]) * odts * orho;
      }
    } else {
      niten[k] = -ni1d[k] * odts;
    }
    xni = c_max(0., (ni1d[k] + niten[k] * dtsave) * rho[k]);
    if (xni > 9999.E3) niten[k] = (9999.E3 - ni1d[k] * rho[k]) * odts * orho;

    qrten[k] = qrten[k] + (prr_wau[k] + prr_rcw[k] + prr_sml[k] + prr_gml[k] + prr_rcs[k] + prr_rcg[k] -
                           prg_rfz[k] - pri_rfz[k] - prr_rci[k]) *
                              orho;

    nrten[k] = nrten[k] + (pnr_wau[k] + pnr_sml[k] + pnr_gml[k] -
                           (pnr_rfz[k] + pnr_rcr[k] + pnr_rcg[k] + pnr_rcs[k] + pnr_rci[k] + pni_rfz[k])) *
                              orho;

    xrr = c_max(R1, (qr1d[k] + qrten[k] * dtsave) * rho[k]);
    xnr = c_max(R2, (nr1d[k] + nrten[k] * dtsave) * rho[k]);
    if (xrr > R1) {
      lamr = pow(am_r * P.crg[2] * P.org2 * xnr / xrr, P.obmr);
      mvd_r[k] = (3.0 + mu_r + 0.672) / lamr;
      if (mvd_r[k] > 2.5E-3) {
        mvd_r[k] = 2.5E-3;
        lamr = (3.0 + mu_r + 0.672) / mvd_r[k];
        xnr = P.crg[1] * P.org3 * xrr * pow(lamr, bm_r) / am_r;
        nrten[k] = (xnr - nr1d[k] * rho[k]) * odts * orho;
      } else if (mvd_r[k] < D0r * 0.75) {
        mvd_r[k] = D0r * 0.75;
        lamr = (3.0 + mu_r + 0.672) / mvd_r[k];
        xnr = P.crg[1] * P.org3 * xrr * pow(lamr, bm_r) / am_r;
        nrten[k] = (xnr - nr1d[k] * rho[k]) * odts * orho;
      }
    } else {
      qrten[k] = -qr1d[k] * odts;
      nrten[k] = -nr1d[k] * odts;
    }

    qsten[k] = qsten[k] + (prs_iau[k] + prs_sde[k] + prs_sci[k] + prs_scw[k] + prs_rcs[k] + prs_ide[k] -
                           prs_ihm[k] - prr_sml[k]) *
                              orho;

    qgten[k] = qgten[k] + (prg_scw[k] + prg_rfz[k] + prg_gde[k] + prg_rcg[k] + prg_gcw[k] + prg_rci[k] +
                           prg_rcs[k] - prg_ihm[k] - prr_gml[k]) *
                              orho;

    if (temp[k] < T_0) {
      tten[k] = tten[k] + (lsub * ocp[k] *
                               (pri_inu[k] + pri_ide[k] + prs_ide[k] + prs_sde[k] + prg_gde[k] + pri_iha[k]) +
                           lfus2 * ocp[k] *
                               (pri_wfz[k] + pri_rfz[k] + prg_rfz[k] + prs_scw[k] + prg_scw[k] + prg_gcw[k] +
                                prg_rcs[k] + prs_rcs[k] + prr_rci[k] + prg_rcg[k])) *
                              orho * (1 - IFDRY);
    } else {
      tten[k] = tten[k] + (lfus * ocp[k] * (-prr_sml[k] - prr_gml[k] - prr_rcg[k] - prr_rcs[k]) +
                           lsub * ocp[k] * (prs_sde[k] + prg_gde[k])) *
                              orho * (1 - IFDRY);
    }
  }

  // Update variables for TAU+1 before condensation & sedimention (MT:2847-2871)
  for (k = kts; k <= kte; ++k) {
    temp[k] = t1d[k] + dt * tten[k];
    otemp = 1. / temp[k];
    tempc = temp[k] - 273.15;
    qv[k] = c_max(1.E-10, qv1d[k] + dt * qvten[k]);
    rho[k] = 0.622 * pres[k] / (R * temp[k] * (qv[k] + 0.622));
    rhof[k] = sqrt(rho_not / rho[k]);
    rhof2[k] = sqrt(rhof[k]);
    qvs[k] = rslf(pres[k], temp[k]);
    ssatw[k] = qv[k] / qvs[k] - 1.;
    if (fabs(ssatw[k]) < eps) ssatw[k] = 0.0;
    diffu[k] = 2.11E-5 * pow(temp[k] / 273.15, 1.94) * (101325. / pres[k]);
    if (tempc >= 0.0) {
      visco[k] = (1.718 + 0.0049 * tempc) * 1.0E-5;
    } else {
      visco[k] = (1.718 + 0.0049 * tempc - 1.2E-5 * tempc * tempc) * 1.0E-5;
    }
    vsc2[k] = sqrt(rho[k] / visco[k]);
    lvap[k] = lvap0 + (2106.0 - 4218.0) * tempc;
    tcond[k] = (5.69 + 0.0168 * tempc) * 1.0E-5 * 418.936;
    ocp[k] = 1. / (Cp * (1. + 0.887 * qv[k]));
    lvt2[k] = lvap[k] * lvap[k] * ocp[k] * oRv * otemp * otemp;
    // MT:2870: nwfa(k), aerosol-aware only
  }

  for (k = kts; k <= kte; ++k) {  // MT:2873-2935
    if ((qc1d[k] + qcten[k] * dt) > R1) {
      rc[k] = (qc1d[k] + qcten[k] * dt) * rho[k];
      nc[k] = c_max(2., c_min((nc1d[k] + ncten[k] * dt) * rho[k], Nt_c_max));
      nc[k] = Nt_c;
      L_qc[k] = true;
    } else {
      rc[k] = R1;
      nc[k] = 2.;
      L_qc[k] = false;
    }

    if ((qi1d[k] + qiten[k] * dt) > R1) {
      ri[k] = (qi1d[k] + qiten[k] * dt) * rho[k];
      ni[k] = c_max(R2, (ni1d[k] + niten[k] * dt) * rho[k]);
      L_qi[k] = true;
    } else {
      ri[k] = R1;
      ni[k] = R2;
      L_qi[k] = false;
    }

    if ((qr1d[k] + qrten[k] * dt) > R1) {
      rr[k] = (qr1d[k] + qrten[k] * dt) * rho[k];
      nr[k] = c_max(R2, (nr1d[k] + nrten[k] * dt) * rho[k]);
      L_qr[k] = true;
      lamr = pow(am_r * P.crg[2] * P.org2 * nr[k] / rr[k], P.obmr);
      mvd_r[k] = (3.0 + mu_r + 0.672) / lamr;
      if (mvd_r[k] > 2.5E-3) {
        mvd_r[k] = 2.5E-3;
        lamr = (3.0 + mu_r + 0.672) / mvd_r[k];
        nr[k] = P.crg[1] * P.org3 * rr[k] * pow(lamr, bm_r) / am_r;
      } else if (mvd_r[k] < D0r * 0.75) {
        mvd_r[k] = D0r * 0.75;
        lamr = (3.0 + mu_r + 0.672) / mvd_r[k];
        nr[k] = P.crg[1] * P.org3 * rr[k] * pow(lamr, bm_r) / am_r;
      }
    } else {
      rr[k] = R1;
      nr[k] = R2;
      L_qr[k] = false;
    }

    if ((qs1d[k] + qsten[k] * dt) > R1) {
      rs[k] = (qs1d[k] + qsten[k] * dt) * rho[k];
      L_qs[k] = true;
    } else {
      rs[k] = R1;
      L_qs[k] = false;
    }

    if ((qg1d[k] + qgten[k] * dt) > R1) {
      rg[k] = (qg1d[k] + qgten[k] * dt) * rho[k];
      L_qg[k] = true;
    } else {
      rg[k] = R1;
      L_qg[k] = false;
    }
  }

  // Snow moments and graupel / rain intercepts again (MT:2941-3025)
  for (k = kts; k <= kte; ++k) {
    smo2[k] = 0.;
    smob[k] = 0.;
    smoc[k] = 0.;
    smod[k] = 0.;
  }
  for (k = kts; k <= kte; ++k) {
    if (!L_qs[k]) continue;
    tc0 = c_min(-0.1, temp[k] - 273.15);
    smob[k] = rs[k] * P.oams;
    smo2[k] = smob[k];
    smoc[k] = snow_moment(tc0, P.cse[0], smo2[k]);   // bm_s+1 (MT:2971-2981)
    smod[k] = snow_moment(tc0, P.cse[13], smo2[k]);  // bm_s+bv_s (MT:2984-2994)
  }

  N0_min = gonv_max;  // MT:3001-3017
  k_0 = kts;
  for (k = kte; k >= kts; --k)
    if (temp[k] >= 270.65) k_0 = i_max(k_0, k);
  for (k = kte; k >= kts; --k) {
    if (k > k_0 && L_qr[k] && mvd_r[k] > 100.E-6) {
      xslw1 = 4.01 + log10(mvd_r[k]);
    } else {
      xslw1 = 0.01;
    }
    ygra1 = 4.31 + log10(c_max(5.E-5, rg[k]));
    zans1 = 3.1 + (100. / (300. * xslw1 * ygra1 / (10. / xslw1 + 1. + 0.25 * ygra1) + 30. + 10. * ygra1));
    N0_exp = pow(10., zans1);
    N0_exp = c_max(gonv_min, c_min(N0_exp, gonv_max));
    N0_min = c_min(N0_exp, N0_min);
    N0_exp = N0_min;
    lam_exp = pow(N0_exp * am_g * P.cgg[0] / rg[k], P.oge1);
    lamg = lam_exp * pow(P.cgg[2] * P.ogg2 * P.ogg1, P.obmg);
    ilamg[k] = 1. / lamg;
    N0_g[k] = N0_exp / (P.cgg[1] * lam_exp) * pow(lamg, P.cge[1]);
  }

  for (k = kte; k >= kts; --k) {  // MT:3020-3025, _ZGVbN2vv_pow (MT:3021, 3024)
    const int li = kte - k, ln = kte - kts + 1;
    lamr = vpow(am_r * P.crg[2] * P.org2 * nr[k] / rr[k], P.obmr, li, ln);
    ilamr[k] = 1. / lamr;
    mvd_r[k] = (3.0 + mu_r + 0.672) / lamr;
    N0_r[k] = nr[k] * P.org2 * vpow(lamr, P.cre[1], li, ln);
  }

  // Cloud water condensation and evaporation (MT:3035-3130)
  for (k = kts; k <= kte; ++k) {
    orho = 1. / rho[k];
    if ((ssatw[k] > eps) || (ssatw[k] < -eps && L_qc[k])) {
      clap = (qv[k] - qvs[k]) / (1. + lvt2[k] * qvs[k]);
      for (n = 1; n <= 3; ++n) {
        fcd = qvs[k] * exp(lvt2[k] * clap) - qv[k] + clap;
        dfcd = qvs[k] * lvt2[k] * exp(lvt2[k] * clap) + 1.;
        clap = clap - fcd / dfcd;
      }
      xrc = rc[k] + clap * rho[k];
      xnc = 0.;
      if (xrc > R1) {
        prw_vcd[k] = clap * odt;
        if (clap > eps) {
          xnc = Nt_c;
          pnc_wcd[k] = 0.5 * (xnc - nc[k] + fabs(xnc - nc[k])) * odts * orho;
        }
        // MT:3057-3109: the evaporation branch needs is_aerosol_aware
      } else {
        prw_vcd[k] = -rc[k] * orho * odt;
        pnc_wcd[k] = -nc[k] * orho * odt;
      }

      qvten[k] = qvten[k] - prw_vcd[k];
      qcten[k] = qcten[k] + prw_vcd[k];
      ncten[k] = ncten[k] + pnc_wcd[k];
      tten[k] = tten[k] + lvap[k] * ocp[k] * prw_vcd[k] * (1 - IFDRY);
      rc[k] = c_max(R1, (qc1d[k] + dt * qcten[k]) * rho[k]);
      if (rc[k] == R1) L_qc[k] = false;
      nc[k] = c_max(2., c_min((nc1d[k] + ncten[k] * dt) * rho[k], Nt_c_max));
      nc[k] = Nt_c;
      qv[k] = c_max(1.E-10, qv1d[k] + dt * qvten[k]);
      temp[k] = t1d[k] + dt * tten[k];
      rho[k] = 0.622 * pres[k] / (R * temp[k] * (qv[k] + 0.622));
      qvs[k] = rslf(pres[k], temp[k]);
      ssatw[k] = qv[k] / qvs[k] - 1.;
    }
  }

  // Rain evaporation (MT:3136-3209)
  for (k = kts; k <= kte; ++k) {
    if ((ssatw[k] < -eps) && L_qr[k] && (!(prw_vcd[k] > 0.))) {
      tempc = temp[k] - 273.15;
      otemp = 1. / temp[k];
      orho = 1. / rho[k];
      rhof[k] = sqrt(rho_not * orho);
      rhof2[k] = sqrt(rhof[k]);
      diffu[k] = 2.11E-5 * pow(temp[k] / 273.15, 1.94) * (101325. / pres[k]);
      if (tempc >= 0.0) {
        visco[k] = (1.718 + 0.0049 * tempc) * 1.0E-5;
      } else {
        visco[k] = (1.718 + 0.0049 * tempc - 1.2E-5 * tempc * tempc) * 1.0E-5;
      }
      vsc2[k] = sqrt(rho[k] / visco[k]);
      lvap[k] = lvap0 + (2106.0 - 4218.0) * tempc;
      tcond[k] = (5.69 + 0.0168 * tempc) * 1.0E-5 * 418.936;
      ocp[k] = 1. / (Cp * (1. + 0.887 * qv[k]));

      rvs = rho[k] * qvs[k];
      rvs_p = rvs * otemp * (lvap[k] * otemp * oRv - 1.);
      rvs_pp = rvs * (otemp * (lvap[k] * otemp * oRv - 1.) * otemp * (lvap[k] * otemp * oRv - 1.) +
                      (-2. * lvap[k] * otemp * otemp * otemp * oRv) + otemp * otemp);
      gamsc = lvap[k] * diffu[k] / tcond[k] * rvs_p;
      alphsc = 0.5 * (gamsc / (1. + gamsc)) * (gamsc / (1. + gamsc)) * rvs_pp / rvs_p * rvs / rvs_p;
      alphsc = c_max(1.E-9, alphsc);
      xsat = c_min(-1.E-9, ssatw[k]);
      t1_evap = 2. * PI *
                (1.0 - alphsc * xsat + 2. * alphsc * alphsc * xsat * xsat -
                 5. * alphsc * alphsc * alphsc * xsat * xsat * xsat) /
                (1. + gamsc);

      lamr = 1. / ilamr[k];
      if (qv[k] / qvs[k] < 0.95 && rr[k] * orho <= 1.E-8) {
        prv_rev[k] = rr[k] * orho * odts;
      } else {
        prv_rev[k] = t1_evap * diffu[k] * (-ssatw[k]) * N0_r[k] * rvs *
                     (P.t1_qr_ev * pow(ilamr[k], P.cre[9]) +
                      P.t2_qr_ev * vsc2[k] * rhof2[k] * (pow(lamr + 0.5 * fv_r, -P.cre[10])));
        rate_max = c_min((rr[k] * orho * odts), (qvs[k] - qv[k]) * odts);
        prv_rev[k] = c_min(rate_max, prv_rev[k] * orho);

        if (prr_gml[k] > 0.0) {
          eva_factor = c_min(1.0, 0.01 + (0.99 - 0.01) * (tempc / 20.0));
          prv_rev[k] = prv_rev[k] * eva_factor;
        }
      }

      pnr_rev[k] = c_min(nr[k] * 0.99 * orho * odts, prv_rev[k] * nr[k] / rr[k]);

      qrten[k] = qrten[k] - prv_rev[k];
      qvten[k] = qvten[k] + prv_rev[k];
      nrten[k] = nrten[k] - pnr_rev[k];
      tten[k] = tten[k] - lvap[k] * ocp[k] * prv_rev[k] * (1 - IFDRY);

      rr[k] = c_max(R1, (qr1d[k] + dt * qrten[k]) * rho[k]);
      qv[k] = c_max(1.E-10, qv1d[k] + dt * qvten[k]);
      nr[k] = c_max(R2, (nr1d[k] + dt * nrten[k]) * rho[k]);
      temp[k] = t1d[k] + dt * tten[k];
      rho[k] = 0.622 * pres[k] / (R * temp[k] * (qv[k] + 0.622));
    }
  }
  for (k = kts; k <= kte; ++k) {  // MT:3210-3217
    evapprod[k] = prv_rev[k] - (c_min(zeroD0, prs_sde[k]) + c_min(zeroD0, prg_gde[k]));
    rainprod[k] = prr_wau[k] + prr_rcw[k] + prs_scw[k] + prg_scw[k] + prs_iau[k] + prg_gcw[k] + prs_sci[k] +
                  pri_rci[k];
  }

  // Fall speeds and the number of sedimentation sub-steps (MT:3227-3424)
  nstep = 0;
  for (int i = 1; i <= 5; ++i) {
    onstep[i] = 1.0;
    ksed1[i] = 1;
  }
  for (k = kte + 1; k >= kts; --k) {
    vtrk[k] = 0.; vtnrk[k] = 0.; vtik[k] = 0.; vtnik[k] = 0.; vtsk[k] = 0.; vtgk[k] = 0.; vtck[k] = 0.;
    vtnck[k] = 0.;
  }

  auto any = [&](const bool* L) {
    for (int kk = kts; kk <= kte; ++kk)
      if (L[kk]) return true;
    return false;
  };

  if (any(L_qr)) {  // MT:3243-3274; MT:3251's two pows paired in one _ZGVbN2vv_pow
    for (k = kte; k >= kts; --k) {
      vtr = 0.;
      rhof[k] = sqrt(rho_not / rho[k]);
      if (rr[k] > R1) {
        lamr = pow(am_r * P.crg[2] * P.org2 * nr[k] / rr[k], P.obmr);
        vtr = rhof[k] * av_r * P.crg[5] * P.org3 * pow(lamr, P.cre[2]) *
              (::cppmpas::libm::mvec_pow(lamr + fv_r, -P.cre[5]));
        vtrk[k] = vtr;
        vtr = rhof[k] * av_r * P.crg[6] / P.crg[11] * pow(lamr, P.cre[11]) *
              (::cppmpas::libm::mvec_pow(lamr + fv_r, -P.cre[6]));
        vtnrk[k] = vtr;
      } else {
        vtrk[k] = vtrk[k + 1];
        vtnrk[k] = vtnrk[k + 1];
      }

      if (c_max(vtrk[k], vtnrk[k]) > 1.E-3) {
        ksed1[1] = i_max(ksed1[1], k);
        delta_tp = dzq[k] / (c_max(vtrk[k], vtnrk[k]));
        nstep = i_max(nstep, static_cast<int>(dt / delta_tp + 1.));
      }
    }
    if (ksed1[1] == kte) ksed1[1] = kte - 1;
    if (nstep > 0) onstep[1] = 1. / double(nstep);
  }

  if (any(L_qc)) {  // MT:3277-3299
    hgt_agl = 0.;
    for (k = kts; k <= kte - 1; ++k) {
      if (rc[k] > R2) ksed1[5] = k;
      hgt_agl = hgt_agl + dzq[k];
      if (hgt_agl > 500.0) break;
    }

    for (k = ksed1[5]; k >= kts; --k) {
      vtc = 0.;
      if (rc[k] > R1 && w1d[k] < 1.E-1) {
        nu_c = i_min(15, nint(1000.E6 / nc[k]) + 2);
        lamc = pow(nc[k] * am_r * P.CCG(2, nu_c) * P.ocg1[nu_c - 1] / rc[k], P.obmr);
        ilamc = 1. / lamc;
        vtc = rhof[k] * av_c * P.CCG(5, nu_c) * P.ocg2[nu_c - 1] * (ilamc * ilamc);  // ilamc**bv_c, bv_c = 2.0
        vtck[k] = vtc;
        vtc = rhof[k] * av_c * P.CCG(4, nu_c) * P.ocg1[nu_c - 1] * (ilamc * ilamc);
        vtnck[k] = vtc;
      }
    }
  }

  if (any(L_qi)) {  // MT:3305-3334
    nstep = 0;
    for (k = kte; k >= kts; --k) {
      vti = 0.;
      if (ri[k] > R1) {
        lami = pow(am_i * P.cig[1] * P.oig1 * ni[k] / ri[k], P.obmi);
        ilami = 1. / lami;
        vti = rhof[k] * av_i * P.cig[2] * P.oig2 * ilami;  // ilami**bv_i, bv_i = 1.0
        vtik[k] = vti;
        vti = rhof[k] * av_i * P.cig[5] / P.cig[6] * ilami;
        vtnik[k] = vti;
      } else {
        vtik[k] = vtik[k + 1];
        vtnik[k] = vtnik[k + 1];
      }

      if (vtik[k] > 1.E-3) {
        ksed1[2] = i_max(ksed1[2], k);
        delta_tp = dzq[k] / vtik[k];
        nstep = i_max(nstep, static_cast<int>(dt / delta_tp + 1.));
      }
    }
    if (ksed1[2] == kte) ksed1[2] = kte - 1;
    if (nstep > 0) onstep[2] = 1. / double(nstep);
  }

  if (any(L_qs)) {  // MT:3338-3376
    nstep = 0;
    for (k = kte; k >= kts; --k) {
      vts = 0.;
      if (rs[k] > R1) {
        xDs = smoc[k] / smob[k];
        Mrat = 1. / xDs;
        ils1 = 1. / (Mrat * Lam0 + fv_s);
        ils2 = 1. / (Mrat * Lam1 + fv_s);
        t1_vts = Kap0 * P.csg[3] * pow(ils1, P.cse[3]);
        t2_vts = Kap1 * pow(Mrat, mu_s) * P.csg[9] * pow(ils2, P.cse[9]);
        ils1 = 1. / (Mrat * Lam0);
        ils2 = 1. / (Mrat * Lam1);
        t3_vts = Kap0 * P.csg[0] * pow(ils1, P.cse[0]);
        t4_vts = Kap1 * pow(Mrat, mu_s) * P.csg[6] * pow(ils2, P.cse[6]);
        vts = rhof[k] * av_s * (t1_vts + t2_vts) / (t3_vts + t4_vts);
        if (temp[k] > (T_0 + 0.1)) {
          SR = rs[k] / (rs[k] + rr[k]);
          vtsk[k] = vts * SR + (1. - SR) * vtrk[k];
        } else {
          vtsk[k] = vts * vts_boost[k];
        }
      } else {
        vtsk[k] = vtsk[k + 1];
      }

      if (vtsk[k] > 1.E-3) {
        ksed1[3] = i_max(ksed1[3], k);
        delta_tp = dzq[k] / vtsk[k];
        nstep = i_max(nstep, static_cast<int>(dt / delta_tp + 1.));
      }
    }
    if (ksed1[3] == kte) ksed1[3] = kte - 1;
    if (nstep > 0) onstep[3] = 1. / double(nstep);
  }

  if (any(L_qg)) {  // MT:3380-3405
    nstep = 0;
    for (k = kte; k >= kts; --k) {
      vtg = 0.;
      if (rg[k] > R1) {
        vtg = rhof[k] * av_g * P.cgg[5] * P.ogg3 * pow(ilamg[k], bv_g);
        if (temp[k] > T_0) {
          vtgk[k] = c_max(vtg, vtrk[k]);
        } else {
          vtgk[k] = vtg;
        }
      } else {
        vtgk[k] = vtgk[k + 1];
      }

      if (vtgk[k] > 1.E-3) {
        ksed1[4] = i_max(ksed1[4], k);
        delta_tp = dzq[k] / vtgk[k];
        nstep = i_max(nstep, static_cast<int>(dt / delta_tp + 1.));
      }
    }
    if (ksed1[4] == kte) ksed1[4] = kte - 1;
    if (nstep > 0) onstep[4] = 1. / double(nstep);
  }

  // Sedimentation (MT:3413-3541)
  if (any(L_qr)) {
    nstep = nint(1. / onstep[1]);
    for (n = 1; n <= nstep; ++n) {
      for (k = kte; k >= kts; --k) {
        sed_r[k] = vtrk[k] * rr[k];
        sed_n[k] = vtnrk[k] * nr[k];
      }
      k = kte;
      odzq = 1. / dzq[k];
      orho = 1. / rho[k];
      qrten[k] = qrten[k] - sed_r[k] * odzq * onstep[1] * orho;
      nrten[k] = nrten[k] - sed_n[k] * odzq * onstep[1] * orho;
      rr[k] = c_max(R1, rr[k] - sed_r[k] * odzq * dt * onstep[1]);
      nr[k] = c_max(R2, nr[k] - sed_n[k] * odzq * dt * onstep[1]);
      for (k = ksed1[1]; k >= kts; --k) {
        odzq = 1. / dzq[k];
        orho = 1. / rho[k];
        qrten[k] = qrten[k] + (sed_r[k + 1] - sed_r[k]) * odzq * onstep[1] * orho;
        nrten[k] = nrten[k] + (sed_n[k + 1] - sed_n[k]) * odzq * onstep[1] * orho;
        rr[k] = c_max(R1, rr[k] + (sed_r[k + 1] - sed_r[k]) * odzq * dt * onstep[1]);
        nr[k] = c_max(R2, nr[k] + (sed_n[k + 1] - sed_n[k]) * odzq * dt * onstep[1]);
      }

      if (rr[kts] > R1 * 1000.) pptrain = pptrain + sed_r[kts] * dt * onstep[1];
    }
  }

  if (any(L_qc)) {
    for (k = kte; k >= kts; --k) {
      sed_c[k] = vtck[k] * rc[k];
      sed_n[k] = vtnck[k] * nc[k];
    }
    for (k = ksed1[5]; k >= kts; --k) {
      odzq = 1. / dzq[k];
      orho = 1. / rho[k];
      qcten[k] = qcten[k] + (sed_c[k + 1] - sed_c[k]) * odzq * orho;
      ncten[k] = ncten[k] + (sed_n[k + 1] - sed_n[k]) * odzq * orho;
      rc[k] = c_max(R1, rc[k] + (sed_c[k + 1] - sed_c[k]) * odzq * dt);
      nc[k] = c_max(10., nc[k] + (sed_n[k + 1] - sed_n[k]) * odzq * dt);
    }
  }

  if (any(L_qi)) {
    nstep = nint(1. / onstep[2]);
    for (n = 1; n <= nstep; ++n) {
      for (k = kte; k >= kts; --k) {
        sed_i[k] = vtik[k] * ri[k];
        sed_n[k] = vtnik[k] * ni[k];
      }
      k = kte;
      odzq = 1. / dzq[k];
      orho = 1. / rho[k];
      qiten[k] = qiten[k] - sed_i[k] * odzq * onstep[2] * orho;
      niten[k] = niten[k] - sed_n[k] * odzq * onstep[2] * orho;
      ri[k] = c_max(R1, ri[k] - sed_i[k] * odzq * dt * onstep[2]);
      ni[k] = c_max(R2, ni[k] - sed_n[k] * odzq * dt * onstep[2]);
      for (k = ksed1[2]; k >= kts; --k) {
        odzq = 1. / dzq[k];
        orho = 1. / rho[k];
        qiten[k] = qiten[k] + (sed_i[k + 1] - sed_i[k]) * odzq * onstep[2] * orho;
        niten[k] = niten[k] + (sed_n[k + 1] - sed_n[k]) * odzq * onstep[2] * orho;
        ri[k] = c_max(R1, ri[k] + (sed_i[k + 1] - sed_i[k]) * odzq * dt * onstep[2]);
        ni[k] = c_max(R2, ni[k] + (sed_n[k + 1] - sed_n[k]) * odzq * dt * onstep[2]);
      }

      if (ri[kts] > R1 * 1000.) pptice = pptice + sed_i[kts] * dt * onstep[2];
    }
  }

  if (any(L_qs)) {
    nstep = nint(1. / onstep[3]);
    for (n = 1; n <= nstep; ++n) {
      for (k = kte; k >= kts; --k) sed_s[k] = vtsk[k] * rs[k];
      k = kte;
      odzq = 1. / dzq[k];
      orho = 1. / rho[k];
      qsten[k] = qsten[k] - sed_s[k] * odzq * onstep[3] * orho;
      rs[k] = c_max(R1, rs[k] - sed_s[k] * odzq * dt * onstep[3]);
      for (k = ksed1[3]; k >= kts; --k) {
        odzq = 1. / dzq[k];
        orho = 1. / rho[k];
        qsten[k] = qsten[k] + (sed_s[k + 1] - sed_s[k]) * odzq * onstep[3] * orho;
        rs[k] = c_max(R1, rs[k] + (sed_s[k + 1] - sed_s[k]) * odzq * dt * onstep[3]);
      }

      if (rs[kts] > R1 * 1000.) pptsnow = pptsnow + sed_s[kts] * dt * onstep[3];
    }
  }

  if (any(L_qg)) {
    nstep = nint(1. / onstep[4]);
    for (n = 1; n <= nstep; ++n) {
      for (k = kte; k >= kts; --k) sed_g[k] = vtgk[k] * rg[k];
      k = kte;
      odzq = 1. / dzq[k];
      orho = 1. / rho[k];
      qgten[k] = qgten[k] - sed_g[k] * odzq * onstep[4] * orho;
      rg[k] = c_max(R1, rg[k] - sed_g[k] * odzq * dt * onstep[4]);
      for (k = ksed1[4]; k >= kts; --k) {
        odzq = 1. / dzq[k];
        orho = 1. / rho[k];
        qgten[k] = qgten[k] + (sed_g[k + 1] - sed_g[k]) * odzq * onstep[4] * orho;
        rg[k] = c_max(R1, rg[k] + (sed_g[k + 1] - sed_g[k]) * odzq * dt * onstep[4]);
      }

      if (rg[kts] > R1 * 1000.) pptgraul = pptgraul + sed_g[kts] * dt * onstep[4];
    }
  }

  // Instant melting / freezing (MT:3547-3567)
  for (k = kts; k <= kte; ++k) {
    xri = c_max(0.0, qi1d[k] + qiten[k] * dt);
    if ((temp[k] > T_0) && (xri > 0.0)) {
      qcten[k] = qcten[k] + xri * odt;
      ncten[k] = ncten[k] + ni1d[k] * odt;
      qiten[k] = qiten[k] - xri * odt;
      niten[k] = -ni1d[k] * odt;
      tten[k] = tten[k] - lfus * ocp[k] * xri * odt * (1 - IFDRY);
    }

    xrc = c_max(0.0, qc1d[k] + qcten[k] * dt);
    if ((temp[k] < HGFR) && (xrc > 0.0)) {
      lfus2 = lsub - lvap[k];
      xnc = nc1d[k] + ncten[k] * dt;
      qiten[k] = qiten[k] + xrc * odt;
      niten[k] = niten[k] + xnc * odt;
      qcten[k] = qcten[k] - xrc * odt;
      ncten[k] = ncten[k] - xnc * odt;
      tten[k] = tten[k] + lfus2 * ocp[k] * xrc * odt * (1 - IFDRY);
    }
  }

  // Apply the tendencies (MT:3572-3640)
  for (k = kts; k <= kte; ++k) {
    t1d[k] = t1d[k] + tten[k] * dt;
    qv1d[k] = c_max(1.E-10, qv1d[k] + qvten[k] * dt);
    qc1d[k] = qc1d[k] + qcten[k] * dt;
    nc1d[k] = c_max(2. / rho[k], c_min(nc1d[k] + ncten[k] * dt, Nt_c_max));
    // MT:3577-3580: nwfa1d / nifa1d, written back by mp_gt_driver only when aerosol-aware

    if (qc1d[k] <= R1) {
      qc1d[k] = 0.0;
      nc1d[k] = 0.0;
    } else {
      nu_c = i_min(15, nint(1000.E6 / (nc1d[k] * rho[k])) + 2);
      lamc = pow(am_r * P.CCG(2, nu_c) * P.ocg1[nu_c - 1] * nc1d[k] / qc1d[k], P.obmr);
      xDc = (bm_r + nu_c + 1.) / lamc;
      if (xDc < D0c) {
        lamc = P.CCE(2, nu_c) / D0c;
      } else if (xDc > D0r * 2.) {
        lamc = P.CCE(2, nu_c) / (D0r * 2.);
      }
      nc1d[k] = c_min(P.CCG(1, nu_c) * P.ocg2[nu_c - 1] * qc1d[k] / am_r * pow(lamc, bm_r), Nt_c_max / rho[k]);
    }

    qi1d[k] = qi1d[k] + qiten[k] * dt;
    ni1d[k] = c_max(R2 / rho[k], ni1d[k] + niten[k] * dt);
    if (qi1d[k] <= R1) {
      qi1d[k] = 0.0;
      ni1d[k] = 0.0;
    } else {
      lami = pow(am_i * P.cig[1] * P.oig1 * ni1d[k] / qi1d[k], P.obmi);
      ilami = 1. / lami;
      xDi = (bm_i + mu_i + 1.) * ilami;
      if (xDi < 5.E-6) {
        lami = P.cie[1] / 5.E-6;
      } else if (xDi > 300.E-6) {
        lami = P.cie[1] / 300.E-6;
      }
      ni1d[k] = c_min(P.cig[0] * P.oig2 * qi1d[k] / am_i * pow(lami, bm_i), 9999.E3 / rho[k]);
    }
    qr1d[k] = qr1d[k] + qrten[k] * dt;
    nr1d[k] = c_max(R2 / rho[k], nr1d[k] + nrten[k] * dt);
    if (qr1d[k] <= R1) {
      qr1d[k] = 0.0;
      nr1d[k] = 0.0;
    } else {
      lamr = pow(am_r * P.crg[2] * P.org2 * nr1d[k] / qr1d[k], P.obmr);
      mvd_r[k] = (3.0 + mu_r + 0.672) / lamr;
      if (mvd_r[k] > 2.5E-3) {
        mvd_r[k] = 2.5E-3;
      } else if (mvd_r[k] < D0r * 0.75) {
        mvd_r[k] = D0r * 0.75;
      }
      lamr = (3.0 + mu_r + 0.672) / mvd_r[k];
      nr1d[k] = P.crg[1] * P.org3 * qr1d[k] * pow(lamr, bm_r) / am_r;
    }
    qs1d[k] = qs1d[k] + qsten[k] * dt;
    if (qs1d[k] <= R1) qs1d[k] = 0.0;
    qg1d[k] = qg1d[k] + qgten[k] * dt;
    if (qg1d[k] <= R1) qg1d[k] = 0.0;
  }
}

// calc_refl10cm (MT:4924-5235) for one column, arrays 1-based (1..kte); dBZ out. Lanes: the dBZ loop (MT:5186-5188,
// k = kte .. kts) is vectorised, _ZGVbN2v_log10 on the pairs from kte, an odd last k = kts and a 1-trip loop with
// log10@plt (the reference's one _ZGVbN2v_log10 is at MT:5187). Every other libm call is scalar.
CPPMPAS_HD inline void calc_refl10cm(const ThompsonView& P, const double* qv1d, const double* qc1d, const double* qr1d,
                                     const double* nr1d, const double* qs1d, const double* qg1d, const double* t1d,
                                     const double* p1d, double* dBZ, int kts, int kte) {
  using namespace detail;
  using ::cppmpas::libm::exp;
  using ::cppmpas::libm::log10;
  using ::cppmpas::libm::pow;
  using std::sqrt;
  (void)qc1d;
  const radar::RadarTables& T = P.radar;
  constexpr int NK = kMaxLevels + 2;
  double temp[NK], pres[NK], qv[NK], rho[NK], rhof[NK];
  double rc[NK], rr[NK], nr[NK], rs[NK], rg[NK];
  double ilamr[NK], ilamg[NK], N0_r[NK], N0_g[NK];
  double mvd_r[NK];
  double smob[NK], smo2[NK], smoc[NK], smoz[NK];
  double oM3, M0, Mrat, slam1, slam2;
  double ze_rain[NK], ze_snow[NK], ze_graupel[NK];
  double N0_exp, N0_min, lam_exp, lamr, lamg;
  double tc0;
  double fmelt_s, fmelt_g;
  int k, k_0, n;
  bool melti;
  bool L_qr[NK], L_qs[NK], L_qg[NK];
  double cback, x, eta, f_d;
  double xslw1, ygra1, zans1;

  for (k = kts; k <= kte; ++k) dBZ[k] = -35.0;

  for (k = kts; k <= kte; ++k) {  // MT:4970-5003
    temp[k] = t1d[k];
    qv[k] = c_max(1.E-10, qv1d[k]);
    pres[k] = p1d[k];
    rho[k] = 0.622 * pres[k] / (R * temp[k] * (qv[k] + 0.622));
    rhof[k] = sqrt(rho_not / rho[k]);
    rc[k] = c_max(R1, qc1d[k] * rho[k]);
    if (qr1d[k] > R1) {
      rr[k] = qr1d[k] * rho[k];
      nr[k] = c_max(R2, nr1d[k] * rho[k]);
      lamr = pow(am_r * P.crg[2] * P.org2 * nr[k] / rr[k], P.obmr);
      ilamr[k] = 1. / lamr;
      N0_r[k] = nr[k] * P.org2 * pow(lamr, P.cre[1]);
      mvd_r[k] = (3.0 + mu_r + 0.672) * ilamr[k];
      L_qr[k] = true;
    } else {
      rr[k] = R1;
      nr[k] = R1;
      mvd_r[k] = 50.E-6;
      L_qr[k] = false;
    }
    if (qs1d[k] > R2) {
      rs[k] = qs1d[k] * rho[k];
      L_qs[k] = true;
    } else {
      rs[k] = R1;
      L_qs[k] = false;
    }
    if (qg1d[k] > R2) {
      rg[k] = qg1d[k] * rho[k];
      L_qg[k] = true;
    } else {
      rg[k] = R1;
      L_qg[k] = false;
    }
  }

  auto any = [&](const bool* L) {
    for (int kk = kts; kk <= kte; ++kk)
      if (L[kk]) return true;
    return false;
  };

  for (k = kts; k <= kte; ++k) {  // MT:5008-5013
    smo2[k] = 0.;
    smob[k] = 0.;
    smoc[k] = 0.;
    smoz[k] = 0.;
  }
  if (any(L_qs)) {  // MT:5014-5058, bm_s = 2: smo2 = smob
    for (k = kts; k <= kte; ++k) {
      if (!L_qs[k]) continue;
      tc0 = c_min(-0.1, temp[k] - 273.15);
      smob[k] = rs[k] * P.oams;
      smo2[k] = smob[k];
      smoc[k] = snow_moment(tc0, P.cse[0], smo2[k]);  // bm_s+1
      smoz[k] = snow_moment(tc0, P.cse[2], smo2[k]);  // bm_s*2
    }
  }

  if (any(L_qg)) {  // MT:5064-5089
    N0_min = gonv_max;
    k_0 = kts;
    for (k = kte; k >= kts; --k)
      if (temp[k] >= 270.65) k_0 = i_max(k_0, k);
    for (k = kte; k >= kts; --k) {
      if (k > k_0 && L_qr[k] && mvd_r[k] > 100.E-6) {
        xslw1 = 4.01 + log10(mvd_r[k]);
      } else {
        xslw1 = 0.01;
      }
      ygra1 = 4.31 + log10(c_max(5.E-5, rg[k]));
      zans1 = 3.1 + (100. / (300. * xslw1 * ygra1 / (10. / xslw1 + 1. + 0.25 * ygra1) + 30. + 10. * ygra1));
      N0_exp = pow(10., zans1);
      N0_exp = c_max(gonv_min, c_min(N0_exp, gonv_max));
      N0_min = c_min(N0_exp, N0_min);
      N0_exp = N0_min;
      lam_exp = pow(N0_exp * am_g * P.cgg[0] / rg[k], P.oge1);
      lamg = lam_exp * pow(P.cgg[2] * P.ogg2 * P.ogg1, P.obmg);
      ilamg[k] = 1. / lamg;
      N0_g[k] = N0_exp / (P.cgg[1] * lam_exp) * pow(lamg, P.cge[1]);
    }
  }

  melti = false;  // MT:5094-5104
  k_0 = kts;
  for (k = kte - 1; k >= kts; --k) {
    if ((temp[k] > 273.15) && L_qr[k] && (L_qs[k + 1] || L_qg[k + 1])) {
      k_0 = i_max(k + 1, k_0);
      melti = true;
      break;
    }
  }

  for (k = kts; k <= kte; ++k) {  // MT:5111-5121
    ze_rain[k] = 1.e-22;
    ze_snow[k] = 1.e-22;
    ze_graupel[k] = 1.e-22;
    if (L_qr[k]) ze_rain[k] = N0_r[k] * P.crg[3] * pow(ilamr[k], P.cre[3]);
    if (L_qs[k])
      ze_snow[k] = (0.176 / 0.93) * (6.0 / PI) * (6.0 / PI) * (am_s / 900.0) * (am_s / 900.0) * smoz[k];
    if (L_qg[k])
      ze_graupel[k] = (0.176 / 0.93) * (6.0 / PI) * (6.0 / PI) * (am_g / 900.0) * (am_g / 900.0) * N0_g[k] *
                      P.cgg[3] * pow(ilamg[k], P.cge[3]);
  }

  if (melti && k_0 >= 2) {  // MT:5131-5179 (iiwarm = .false.)
    for (k = k_0 - 1; k >= kts; --k) {
      if (L_qs[k] && L_qs[k_0]) {
        fmelt_s = c_max(0.05, c_min(1.0 - rs[k] / rs[k_0], 0.99));
        eta = 0.;
        oM3 = 1. / smoc[k];
        M0 = (smob[k] * oM3);
        Mrat = smob[k] * M0 * M0 * M0;
        slam1 = M0 * Lam0;
        slam2 = M0 * Lam1;
        for (n = 1; n <= radar::nrbins; ++n) {
          x = am_s * (T.xxDs[n - 1] * T.xxDs[n - 1]);  // xxDs(n)**bm_s, bm_s = 2.0
          cback = radar::rayleigh_soak_wetgraupel(x, P.ocms, P.obms, fmelt_s, T.m_w_0, T.m_i_0, T.pi5, T.lamda4);
          f_d = Mrat * (Kap0 * exp(-slam1 * T.xxDs[n - 1]) +
                        Kap1 * pow(M0 * T.xxDs[n - 1], mu_s) * exp(-slam2 * T.xxDs[n - 1]));
          eta = eta + f_d * cback * T.simpson[n - 1] * T.xdts[n - 1];
        }
        ze_snow[k] = T.lamda4 / (T.pi5 * T.k_w) * eta;
      }

      if (L_qg[k] && L_qg[k_0]) {
        fmelt_g = c_max(0.05, c_min(1.0 - rg[k] / rg[k_0], 0.99));
        eta = 0.;
        lamg = 1. / ilamg[k];
        for (n = 1; n <= radar::nrbins; ++n) {
          x = am_g * pow(T.xxDg[n - 1], bm_g);
          cback = radar::rayleigh_soak_wetgraupel(x, P.ocmg, P.obmg, fmelt_g, T.m_w_0, T.m_i_0, T.pi5, T.lamda4);
          f_d = N0_g[k] * exp(-lamg * T.xxDg[n - 1]);  // xxDg(n)**mu_g, mu_g = 0.0: 1
          eta = eta + f_d * cback * T.simpson[n - 1] * T.xdtg[n - 1];
        }
        ze_graupel[k] = T.lamda4 / (T.pi5 * T.k_w) * eta;
      }
    }
  }

  for (k = kte; k >= kts; --k) {  // MT:5186-5188
    const int li = kte - k, ln = kte - kts + 1;
    const bool lane = ln > 1 && !((ln & 1) && li == ln - 1);
    const double z = (ze_rain[k] + ze_snow[k] + ze_graupel[k]) * 1.e18;
    dBZ[k] = 10. * (lane ? ::cppmpas::libm::mvec_log10(z) : log10(z));
  }
}

}  // namespace cppmpas::phys::thompson
