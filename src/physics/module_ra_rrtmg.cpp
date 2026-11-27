// SPDX-License-Identifier: GPL-3.0-or-later AND LicenseRef-AER-RRTMG
// See module_ra_rrtmg.hpp. LW = physics/physics_wrf/module_ra_rrtmg_lw.F, SW = module_ra_rrtmg_sw.F, VI =
// module_ra_rrtmg_vinterp.F (MPAS v8.4.2). With -fdefault-real-8 every REAL and literal is a double.
#include "physics/module_ra_rrtmg.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <string>

#include "core/error.hpp"
#include "core/field.hpp"
#include "core/libm.hpp"
#include "physics/mpas_atmphys_constants.hpp"
#include "physics/rrtmg/kissvec.hpp"
#include "physics/rrtmg/rrtmg_lw_stages.hpp"
#include "physics/rrtmg/rrtmg_sw_stages.hpp"

namespace cppmpas::phys {

namespace {

using rrtmg::rb;
using rrtmg::SD;

// retab (LW 11547-11566): re(T) of hexagonal columns, 180 K .. 274 K
constexpr real retab_data[95] = {
    5.92779, 6.26422, 6.61973, 6.99539, 7.39234, 7.81177, 8.25496, 8.72323, 9.21800, 9.74075, 10.2930,
    10.8765, 11.4929, 12.1440, 12.8317, 13.5581, 14.2319, 15.0351, 15.8799, 16.7674, 17.6986, 18.6744,
    19.6955, 20.7623, 21.8757, 23.0364, 24.2452, 25.5034, 26.8125, 27.7895, 28.6450, 29.4167, 30.1088,
    30.7306, 31.2943, 31.8151, 32.3077, 32.7870, 33.2657, 33.7540, 34.2601, 34.7892, 35.3442, 35.9255,
    36.5316, 37.1602, 37.8078, 38.4720, 39.1508, 39.8442, 40.5552, 41.2912, 42.0635, 42.8876, 43.7863,
    44.7853, 45.9170, 47.2165, 48.7221, 50.4710, 52.4980, 54.8315, 57.4898, 60.4785, 63.7898, 65.5604,
    71.2885, 75.4113, 79.7368, 84.2351, 88.8833, 93.6658, 98.5739, 103.603, 108.752, 114.025, 119.424,
    124.954, 130.630, 136.457, 142.446, 148.608, 154.956, 161.503, 168.262, 175.248, 182.473, 189.952,
    197.699, 205.728, 214.055, 222.694, 231.661, 240.971, 250.639};
// the standard profile above the model top (LW 11726-11749)
constexpr real pprof_data[60] = {
    1000.00, 855.47, 731.82, 626.05, 535.57, 458.16, 391.94, 335.29, 286.83, 245.38, 209.91, 179.57,
    153.62,  131.41, 112.42, 96.17,  82.27,  70.38,  60.21,  51.51,  44.06,  37.69,  32.25,  27.59,
    23.60,   20.19,  17.27,  14.77,  12.64,  10.81,  9.25,   7.91,   6.77,   5.79,   4.95,   4.24,
    3.63,    3.10,   2.65,   2.27,   1.94,   1.66,   1.42,   1.22,   1.04,   0.89,   0.76,   0.65,
    0.56,    0.48,   0.41,   0.35,   0.30,   0.26,   0.22,   0.19,   0.16,   0.14,   0.12,   0.10};
constexpr real tprof_data[60] = {
    286.96, 281.07, 275.16, 268.11, 260.56, 253.02, 245.62, 238.41, 231.57, 225.91, 221.72, 217.79,
    215.06, 212.74, 210.25, 210.16, 210.69, 212.14, 213.74, 215.37, 216.82, 217.94, 219.03, 220.18,
    221.37, 222.64, 224.16, 225.88, 227.63, 229.51, 231.50, 233.73, 236.18, 238.78, 241.60, 244.44,
    247.35, 250.33, 253.32, 256.30, 259.22, 262.12, 264.80, 266.50, 267.59, 268.44, 268.69, 267.76,
    266.13, 263.96, 261.54, 258.93, 256.15, 253.23, 249.89, 246.67, 243.48, 240.25, 236.66, 233.86};
// O3DATA (LW 12345-12384)
constexpr real o3sum_data[31] = {
    5.297E-8, 5.852E-8, 6.579E-8, 7.505E-8, 8.577E-8, 9.895E-8, 1.175E-7, 1.399E-7, 1.677E-7, 2.003E-7, 2.571E-7,
    3.325E-7, 4.438E-7, 6.255E-7, 8.168E-7, 1.036E-6, 1.366E-6, 1.855E-6, 2.514E-6, 3.240E-6, 4.033E-6, 4.854E-6,
    5.517E-6, 6.089E-6, 6.689E-6, 1.106E-5, 1.462E-5, 1.321E-5, 9.856E-6, 5.960E-6, 5.960E-6};
constexpr real ppsum_data[31] = {955.890, 850.532, 754.599, 667.742, 589.841, 519.421, 455.480, 398.085,
                                 347.171, 301.735, 261.310, 225.360, 193.419, 165.490, 141.032, 120.125,
                                 102.689, 87.829,  75.123,  64.306,  55.086,  47.209,  40.535,  34.795,
                                 29.865,  19.122,  9.277,   4.660,   2.421,   1.294,   0.647};
constexpr real o3win_data[31] = {
    4.629E-8, 4.686E-8, 5.017E-8, 5.613E-8, 6.871E-8, 8.751E-8, 1.138E-7, 1.516E-7, 2.161E-7, 3.264E-7, 4.968E-7,
    7.338E-7, 1.017E-6, 1.308E-6, 1.625E-6, 2.011E-6, 2.516E-6, 3.130E-6, 3.840E-6, 4.703E-6, 5.486E-6, 6.289E-6,
    6.993E-6, 7.494E-6, 8.197E-6, 9.632E-6, 1.113E-5, 1.146E-5, 9.389E-6, 6.135E-6, 6.135E-6};
constexpr real ppwin_data[31] = {955.747, 841.783, 740.199, 649.538, 568.404, 495.815, 431.069, 373.464,
                                 322.354, 277.190, 237.635, 203.433, 174.070, 148.949, 127.408, 108.915,
                                 93.114,  79.551,  67.940,  58.072,  49.593,  42.318,  36.138,  30.907,
                                 26.362,  16.423,  7.583,   3.620,   1.807,   0.938,   0.469};

// NINT of a REAL (LW 11822): the reference calls lround
CPPMPAS_HD inline long nint(double x) {
#if defined(CPPMPAS_DEVICE_CODE)
  return ::lround(x);
#else
  return std::lround(x);
#endif
}

// Columns per launch of the CUDA build (the work arrays of a chunk: about 1.3 GB for the LW, 1.9 GB for the SW at
// 60 layers); the CPU build runs one column at a time
#if defined(CPPMPAS_USE_CUDA)
constexpr int kChunk = 4096;
#else
constexpr int kChunk = 1;
#endif

// The wrappers' column arrays that are not rrtmg_lw's / rrtmg_sw's inputs (LW 11613-11663, SW 10086-10130), in the
// chunk's work buffer with the column index fastest, as the stages' views (rrtmg_lw_stages.hpp LwLayout); element
// k of MPAS's 1-based array at [k - 1]. The other column arrays of the wrappers are the stages' views themselves:
// plev = pz, tlev = tz, play = pavel, tlay = tavel, h2ovmr / co2vmr / o3vmr / n2ovmr / ch4vmr / o2vmr = wkl(0, 1,
// 2, 3, 5, 6), ccl4vmr / cfc11vmr / cfc12vmr / cfc22vmr = wx(0..3) (LW only), cldfrac, clwpth, ciwpth, cswpth,
// rel = relqmc, rei = reicmc, res = resnmc, o3mmr (LW: the LwViews' o3mmr).
struct RadWork {
  SD pw1d, tw1d, p1d, dz1d, t1d, qv1d, qc1d, qi1d, qs1d, cldfra1d, pdel, reliq, reice, resnow, cliqwp, cicewp, csnowp;
  SD hlay, varint, o3mmr, o3clim1d;
};
struct RadWorkLayout {
  int nlev, nl, noz;  // kte - kts + 1, the layer dimension, noznlevels
  std::size_t nd = 0;
  CPPMPAS_HD RadWorkLayout(int nlev_, int nl_, int noz_) : nlev(nlev_), nl(nl_), noz(noz_) {
    views(nullptr, 0, 1, &nd);
  }
  CPPMPAS_HD RadWork views(rb* d, int col, int ncol, std::size_t* nd_out = nullptr) const {
    RadWork w{};
    std::size_t od = 0;
    auto v1 = [&](int n) {
      rb* q = d ? d + od * std::size_t(ncol) + col : nullptr;
      od += std::size_t(n);
      return SD{q, ncol};
    };
    const int L1 = nlev + 1;
    w.pw1d = v1(L1); w.tw1d = v1(L1); w.p1d = v1(L1); w.dz1d = v1(L1); w.t1d = v1(L1); w.qv1d = v1(L1);
    w.qc1d = v1(L1); w.qi1d = v1(L1); w.qs1d = v1(L1); w.cldfra1d = v1(L1); w.pdel = v1(L1); w.reliq = v1(L1);
    w.reice = v1(L1); w.resnow = v1(L1); w.cliqwp = v1(L1); w.cicewp = v1(L1); w.csnowp = v1(L1);
    w.hlay = v1(nl + 1); w.varint = v1(nl + 1); w.o3mmr = v1(nl + 1); w.o3clim1d = v1(noz);
    if (nd_out) *nd_out = od;
    return w;
  }
};

// relcalc (LW 12447-12497) for one column (ncol = pcols = 1): rel(1:pver) from t(1:pver)
CPPMPAS_HD inline void relcalc(int pver, SD t, real landfrac, real landm, real icefrac, SD rel, real snowh) {
  (void)landfrac;
  const real tmelt = 273.16, rliqocean = 14.0, rliqice = 14.0, rliqland = 8.0;  // LW 12483-12486
  for (int k = 1; k <= pver; ++k) {
    real r = rel[k - 1];
    const real tk = t[k - 1];
    r = rliqland + (rliqocean - rliqland) * std::min(1.0, std::max(0.0, (tmelt - tk) * 0.05));  // LW 12491
    r = r + (rliqocean - r) * std::min(1.0, std::max(0.0, snowh * 10.));                       // LW 12493
    r = r + (rliqocean - r) * std::min(1.0, std::max(0.0, 1.0 - landm));                       // LW 12495
    r = r + (rliqice - r) * std::min(1.0, std::max(0.0, icefrac));                             // LW 12497
    rel[k - 1] = r;
  }
}

// reicalc (LW 12499-12527) for one column
CPPMPAS_HD inline void reicalc(int pver, SD t, SD re, const real* retab) {
  for (int k = 1; k <= pver; ++k) {
    const real tk = t[k - 1];
    int index = static_cast<int>(tk - 179.);       // LW 12519 (INT truncates)
    index = std::min(std::max(index, 1), 94);      // LW 12520
    const real corr = tk - static_cast<int>(tk);   // LW 12521
    re[k - 1] = retab[index - 1] * (1. - corr) + retab[index] * corr;  // LW 12522-12523
  }
}

// inirad (LW 12308-12331) with o3data (LW 12334-12444): o3prof(kts:kte+1) from plev(kts:kte+2); element k at
// [k - kts]. MPAS's copy of plev into prlevh (LW 12412-12414) is read in place.
CPPMPAS_HD inline void inirad(SD o3prof, SD plev, int kts, int kte, const RrtmgMpasConsts& kc) {
  auto O3 = [&](int k) -> rb& { return o3prof[k - kts]; };
  auto PR = [&](int k) { return plev[k - kts]; };
  for (int k = kts; k <= kte + 1; ++k) O3(k) = 0.;  // LW 12326-12328
  real ppann[31], o3ann[31], o3wrk[31], ppwrk[31], ppwrkh[32];
  for (int k = 0; k < 31; ++k) ppann[k] = kc.ppsum[k];  // LW 12386-12388
  o3ann[0] = 0.5 * (kc.o3sum[0] + kc.o3win[0]);          // LW 12390
  for (int k = 1; k < 31; ++k)                            // LW 12392-12395
    o3ann[k] = kc.o3win[k - 1] +
               (kc.o3win[k] - kc.o3win[k - 1]) / (kc.ppwin[k] - kc.ppwin[k - 1]) * (kc.ppsum[k] - kc.ppwin[k - 1]);
  for (int k = 1; k < 31; ++k) o3ann[k] = 0.5 * (o3ann[k] + kc.o3sum[k]);  // LW 12397-12399
  for (int k = 0; k < 31; ++k) {                                            // LW 12401-12404
    o3wrk[k] = o3ann[k];
    ppwrk[k] = ppann[k];
  }
  ppwrkh[0] = 1100.;                                                       // LW 12416
  for (int k = 1; k < 31; ++k) ppwrkh[k] = (ppwrk[k] + ppwrk[k - 1]) / 2.;  // LW 12417-12419
  ppwrkh[31] = 0.;                                                         // LW 12420
  for (int k = kts; k <= kte + 1; ++k) {                                   // LW 12421-12442
    for (int jj = 1; jj <= 31; ++jj) {
      const real pb1 = (-(PR(k) - ppwrkh[jj - 1])) >= 0. ? 0. : PR(k) - ppwrkh[jj - 1];
      const real pb2 = (-(PR(k) - ppwrkh[jj])) >= 0. ? 0. : PR(k) - ppwrkh[jj];
      const real pt1 = (-(PR(k + 1) - ppwrkh[jj - 1])) >= 0. ? 0. : PR(k + 1) - ppwrkh[jj - 1];
      const real pt2 = (-(PR(k + 1) - ppwrkh[jj])) >= 0. ? 0. : PR(k + 1) - ppwrkh[jj];
      O3(k) = O3(k) + (pb2 - pb1 - pt2 + pt1) * o3wrk[jj - 1];
    }
    O3(k) = O3(k) / (PR(k) - PR(k + 1));
  }
}

// mcica_subcol_lw / _sw (LW 2099-2248, SW 1430-1555) with generate_stochastic_clouds (irng = 0): CPPWRF's
// KissRng and mcica_subcolumns (physics/rrtmg/kissvec.hpp). play in mb (pmid = play*1.e2, LW 2201); cldfrac(1:nlay);
// cldy(isubcol, ilev, cloudy) receives the subcolumns of the layers with cldf > 0 (for the others the Fortran's
// cldfmcl stays the 0 rrtmg_lwrad / rrtmg_swrad set). Returns false (no call of cldy) for icld = 0 (LW 2194).
template <class OUT>
CPPMPAS_HD bool mcica_subcol(int icld, int permuteseed, int nsubcol, int nlay, SD play, SD cldfrac, SD hgt,
                             int idcor, int juldat, real lat, OUT cldy) {
  if (icld == 0) return false;  // LW 2194
  if (icld < 0 || icld > 5) kernel_fail("mcica_subcol: MCICA_SUBCOL: INVALID ICLD");
  if (icld == 4)
    kernel_fail("config_radt_cld_overlap = 'exponential' has no verified port (MPAS's -O3 runs its alpha loop in "
                "_ZGVbN2v_exp lanes, module_ra_rrtmg_lw.F:2552)");
  if (nlay > rrtmg::kMcicaMaxLayers) kernel_fail("mcica_subcol: more than 128 layers");
  // generate_stochastic_clouds: KISS seed from the bottom four layers (LW 2441-2452)
  if (play[0] * 1.e2 < play[1] * 1.e2)
    kernel_fail("MCICA_SUBCOL: KISSVEC SEED GENERATOR REQUIRES PMID FROM BOTTOM FOUR LAYERS.");
  rrtmg::KissRng rng;
  rng.init_from_pressure(play[0] * 1.e2, play[1] * 1.e2, play[2] * 1.e2, play[3] * 1.e2);
  rng.permute(permuteseed);
  // Zo_inv = 1 / decorr_len (LW 2225-2240, 2420); alpha of the exponential-random overlap (LW 2590-2603)
  const real zo_inv = 1. / rrtmg::mcica_decorr_len(idcor, juldat, lat);
  auto cldf_of = [&](int l) { const real cf = cldfrac[l]; return cf < 1.0e-20 ? 0.0 : cf; };  // LW 2425-2432
  auto alpha_of = [&](int l) {
    if (icld == 5 && cldf_of(l) == 0.0 && cldf_of(l - 1) > 0.0) return 0.0;
    return libm::exp(-(hgt[l] - hgt[l - 1]) * zo_inv);
  };
  rrtmg::mcica_subcolumns(icld, rng, nsubcol, nlay, cldf_of, alpha_of, [&](int isub, int ilev, bool cld) {
    if (cldf_of(ilev) > 0.0) cldy(isub, ilev, cld);
  });
  return true;
}

// ---- rrtmg_lwrad (LW 11574-12305) -------------------------------------------------------------------------------
// trace gases (LW 11681-11704), molecular weight ratios (LW 11718-11720), buffer layer (LW 11569)
constexpr real lw_co2 = 379.e-6, lw_ch4 = 1774.e-9, lw_n2o = 319.e-9, lw_cfc11 = 0.251e-9, lw_cfc12 = 0.538e-9,
               lw_cfc22 = 0.169e-9, lw_ccl4 = 0.093e-9, lw_o2 = 0.209488;
constexpr real amdw = 1.607793, amdo = 0.603461;
constexpr real deltap = 4.;

// The number of layers of column i (LW 11782-11784, 11822): NINT is lround
CPPMPAS_HD inline int lw_nlayers(const RrtmgRadIn& in, int i) {
  const real pw1d_top = in.p8w(i, in.kte + 1) / 100.;  // pw1d(kte+1) (LW 11784)
  return in.kte + std::max(static_cast<int>(nint(pw1d_top / deltap)), 1);
}

// Column i of rrtmg_lwrad's column loop (LW 11766-11768) up to rrtmg_lw's inatm and the McICA subcolumns
// (LW 11771-12191), then inatm (LW 11197-11533) and the diffusivity angles; setcoef is left to the per-layer stage
// (iscal[0] = laytrop counted there).
CPPMPAS_HD void rrtmg_lwrad_column(const RrtmgRadIn& in, PCView emiss, int i, int nlayers, const rrtmg::LwViews& v,
                                   const RadWork& w, const RrtmgMpasConsts& kc) {
  using physics_constants::gravity;
  const int kts = 1, kte = in.kte;
  const int icld = in.cldovrlp;  // LW 11754
  int inflglw = 2, iceflglw = 3;  // LW 11757-11759
  const int liqflglw = 1;
  auto K = [](int k) { return k - 1; };  // 1-based column arrays

  const real lat = in.xlat(i);  // LW 11771
  // column sounding (LW 11782-11800)
  for (int k = kts; k <= kte + 1; ++k) {
    w.pw1d[K(k)] = in.p8w(i, k) / 100.;
    w.tw1d[K(k)] = in.t8w(i, k);
  }
  for (int k = kts; k <= kte; ++k) {
    w.p1d[K(k)] = in.p3d(i, k) / 100.;
    w.dz1d[K(k)] = in.dz8w(i, k);
    w.t1d[K(k)] = in.t3d(i, k);
    w.qv1d[K(k)] = std::max(in.qv3d(i, k), 1.e-12);
    w.qc1d[K(k)] = 0.;
    w.qi1d[K(k)] = 0.;
    w.qs1d[K(k)] = 0.;
    w.cldfra1d[K(k)] = 0.;
  }
  if (in.icloud > 0) {  // LW 11803-11812 (qr3d is not passed: qr1d stays 0 and is not used)
    for (int k = kts; k <= kte; ++k) {
      w.qc1d[K(k)] = std::max(in.qc3d(i, k), 0.);
      w.qi1d[K(k)] = std::max(in.qi3d(i, k), 0.);
      w.qs1d[K(k)] = std::max(in.qs3d(i, k), 0.);
      w.cldfra1d[K(k)] = in.cldfra3d(i, k);
    }
  }
  // buffer layers (LW 11822): lw_nlayers, passed in
  // LW 11872-11914: local arrays (reliq.. on kts:kte, the rest 0)
  for (int k = kts; k <= kte; ++k) {
    w.reliq[K(k)] = 10.;
    w.reice[K(k)] = 10.;
    w.resnow[K(k)] = 10.;
    w.cliqwp[K(k)] = 0.;
    w.cicewp[K(k)] = 0.;
    w.csnowp[K(k)] = 0.;
  }
  // input sounding (LW 11917-11937)
  v.pz[K(1)] = w.pw1d[K(1)];
  v.tz[K(1)] = w.tw1d[K(1)];
  for (int k = kts; k <= kte; ++k) {
    v.pz[K(k + 1)] = w.pw1d[K(k + 1)];
    v.tz[K(k + 1)] = w.tw1d[K(k + 1)];
    v.pavel[K(k)] = w.p1d[K(k)];
    w.pdel[K(k)] = v.pz[K(k)] - v.pz[K(k + 1)];
    v.tavel[K(k)] = w.t1d[K(k)];
    v.wkl[0][K(k)] = w.qv1d[K(k)] * amdw;  // h2ovmr
    v.wkl[1][K(k)] = lw_co2;
    v.wkl[6][K(k)] = lw_o2;
    v.wkl[5][K(k)] = lw_ch4;
    v.wkl[3][K(k)] = lw_n2o;
    v.wx[1][K(k)] = lw_cfc11;
    v.wx[2][K(k)] = lw_cfc12;
    v.wx[3][K(k)] = lw_cfc22;
    v.wx[0][K(k)] = lw_ccl4;
  }
  // layer mid-point heights (LW 11942-11950)
  real dzsum = 0.;
  for (int k = kts; k <= kte; ++k) {
    w.hlay[K(k)] = dzsum + 0.5 * w.dz1d[K(k)];
    dzsum = dzsum + w.dz1d[K(k)];
  }
  for (int k = kte + 1; k <= nlayers; ++k) {
    w.hlay[K(k)] = dzsum + 0.5 * w.dz1d[K(kte)];
    dzsum = dzsum + w.dz1d[K(kte)];
  }
  // buffer pressures (LW 11955-11964)
  for (int k = kte + 1; k <= nlayers; ++k) {
    v.pz[K(k + 1)] = std::max(v.pz[K(k)] - deltap, 0.00);
    v.pavel[K(k)] = 0.5 * (v.pz[K(k)] + v.pz[K(k + 1)]);
  }
  v.pz[K(nlayers + 1)] = 0.00;
  v.pavel[K(nlayers)] = 0.5 * (v.pz[K(nlayers)] + v.pz[K(nlayers + 1)]);
  // table temperatures at the column pressures (LW 11967-11989)
  constexpr int nproflevs = 60;
  for (int k = 1; k <= nlayers + 1; ++k) {
    int klev = nproflevs;
    if (kc.pprof[nproflevs - 1] < v.pz[K(k)]) {
      for (int kk = 2; kk <= nproflevs; ++kk) {
        if (kc.pprof[kk - 1] < v.pz[K(k)]) {
          klev = kk - 1;
          break;
        }
      }
    } else {
      klev = nproflevs;
    }
    real vark, vark1, wght;
    if (klev != nproflevs) {
      vark = kc.tprof[klev - 1];
      vark1 = kc.tprof[klev];
      wght = (v.pz[K(k)] - kc.pprof[klev - 1]) / (kc.pprof[klev] - kc.pprof[klev - 1]);
    } else {
      vark = kc.tprof[klev - 1];
      vark1 = kc.tprof[klev - 1];
      wght = 0.0;
    }
    w.varint[K(k)] = wght * (vark1 - vark) + vark;
  }
  for (int k = kte + 1; k <= nlayers + 1; ++k) {  // LW 11992-11997
    v.tz[K(k)] = w.varint[K(k)] + (v.tz[K(kte)] - w.varint[K(kte)]);
    v.tavel[K(k - 1)] = 0.5 * (v.tz[K(k)] + v.tz[K(k - 1)]);
  }
  for (int k = kte + 1; k <= nlayers; ++k) {  // LW 12000-12010
    v.wkl[0][K(k)] = v.wkl[0][K(kte)];
    v.wkl[1][K(k)] = v.wkl[1][K(kte)];
    v.wkl[6][K(k)] = v.wkl[6][K(kte)];
    v.wkl[5][K(k)] = v.wkl[5][K(kte)];
    v.wkl[3][K(k)] = v.wkl[3][K(kte)];
    v.wx[1][K(k)] = v.wx[1][K(kte)];
    v.wx[2][K(k)] = v.wx[2][K(kte)];
    v.wx[3][K(k)] = v.wx[3][K(kte)];
    v.wx[0][K(k)] = v.wx[0][K(kte)];
  }
  // ozone (LW 12013-12027)
  inirad(v.o3mmr, v.pz, kts, nlayers - 1, kc);
  if (in.o3input == 2) {
    for (int k = 1; k <= in.noznlevels; ++k) w.o3clim1d[K(k)] = in.o3clim(i, k);
    vinterp_ozn_column(
        nlayers, [&](int k) { return v.pavel[K(k)]; }, [&](int k) { return in.pin(k); }, in.noznlevels,
        [&](int k) { return w.o3clim1d[K(k)]; }, [&](int k) -> rb& { return v.o3mmr[K(k)]; });
    for (int k = kts; k <= nlayers; ++k) v.wkl[2][K(k)] = v.o3mmr[K(k)];
  } else {
    for (int k = kts; k <= nlayers; ++k) v.wkl[2][K(k)] = v.o3mmr[K(k)] * amdo;
  }
  // cloud optical properties (LW 12032-12176); inflglw = 2 > 0
  for (int k = kts; k <= kte; ++k) v.cldfrac[K(k)] = w.cldfra1d[K(k)];
  const real gravmks = gravity;  // LW 12048
  {  // has_reqc = 0: relcalc below the model top (LW 12069-12074)
    const real landfrac = 2. - in.xland(i);
    const real landm = landfrac;
    const real snowh = 0.001 * in.snow(i);
    const real icefrac = in.xice(i);
    relcalc(kte - kts + 1, v.tavel, landfrac, landm, icefrac, w.reliq, snowh);
  }
  reicalc(kte - kts + 1, v.tavel, w.reice, kc.retab);  // has_reqi = 0 (LW 12099-12105)
  if (iceflglw == 3) {
    for (int k = kts; k <= kte; ++k) {
      w.reice[K(k)] = w.reice[K(k)] * 1.0315;
      w.reice[K(k)] = std::min(140.0, w.reice[K(k)]);
    }
  }
  for (int k = kts; k <= kte; ++k) w.resnow[K(k)] = 10.;  // has_reqs = 0 (LW 12117-12121)
  for (int k = kts; k <= kte; ++k) {                       // LW 12124-12127
    const real gliqwp = w.qc1d[K(k)] * w.pdel[K(k)] * 100.0 / gravmks * 1000.0;
    w.cliqwp[K(k)] = gliqwp / std::max(0.01, v.cldfrac[K(k)]);
  }
  if (iceflglw == 3) {  // LW 12130-12135
    for (int k = kts; k <= kte; ++k) {
      const real gicewp = (w.qi1d[K(k)] + w.qs1d[K(k)]) * w.pdel[K(k)] * 100.0 / gravmks * 1000.0;
      w.cicewp[K(k)] = gicewp / std::max(0.01, v.cldfrac[K(k)]);
    }
  }
  for (int k = kts; k <= kte; ++k) {  // LW 12157-12164
    v.clwpth[K(k)] = w.cliqwp[K(k)];
    v.ciwpth[K(k)] = w.cicewp[K(k)];
    v.cswpth[K(k)] = w.csnowp[K(k)];
    v.relqmc[K(k)] = w.reliq[K(k)];
    v.reicmc[K(k)] = w.reice[K(k)];
    v.resnmc[K(k)] = w.resnow[K(k)];
  }
  for (int k = kte + 1; k <= nlayers; ++k) {  // LW 12165-12175
    v.clwpth[K(k)] = 0.;
    v.ciwpth[K(k)] = 0.;
    v.cswpth[K(k)] = 0.;
    v.relqmc[K(k)] = 10.;
    v.reicmc[K(k)] = 10.;
    v.resnmc[K(k)] = 10.;
    v.cldfrac[K(k)] = 0.;
  }

  // rrtmg_lw (LW 10705-11194): inatm's inputs (LW 11197-11533) besides the wrapper's arrays above
  for (int l = 0; l < nlayers; ++l) {
    v.wkl[4][l] = rb(0.0);           // co
    for (int g = 0; g < 140; ++g) {  // cldfmcl, taucmcl = 0 (LW 11899-11903)
      v.cldfmc[g][l] = rb(0.0);
      v.taucmc[g][l] = rb(0.0);
    }
  }
  for (int nb = 0; nb < 16; ++nb) v.semiss[nb] = emiss(i);  // emis (LW 11774-11776), semiss (inatm)
  // McICA (LW 12181-12191): iplon = 1, irng = 0, permuteseed = 150
  std::uint64_t mask_lo = 0, mask_hi = 0;
  mcica_subcol(icld, 150, 140, nlayers, v.pavel, v.cldfrac, w.hlay, in.idcor, in.julday, lat,
               [&](int isub, int ilev, bool cld) {
                 v.cldfmc[isub][ilev] = cld ? rb(1.0) : rb(0.0);
                 if (cld) {
                   if (ilev < 64) mask_lo |= (std::uint64_t(1) << ilev);
                   else mask_hi |= (std::uint64_t(1) << (ilev - 64));
                 }
               });
  int top_cld_plus1 = 0;
  for (int ilev = 0; ilev < nlayers; ++ilev) {
    const int bit = ilev < 64 ? int((mask_lo >> ilev) & 1u) : int((mask_hi >> (ilev - 64)) & 1u);
    v.icldlyr[ilev] = bit;
    if (bit) top_cld_plus1 = ilev + 1;
  }
  v.iscal[1] = inflglw;
  v.iscal[2] = iceflglw;
  v.iscal[3] = liqflglw;
  v.iscal[4] = top_cld_plus1;
  rrtmg::rrtmg_lw_inatm(nlayers, v);
  rrtmg::rrtmg_lw_secdiff(v.dscal[0], v.secdiff);
  v.iscal[0] = 0;  // laytrop, counted by rrtmg_lw_setcoef_stage
}

// rrtmg_lw's heating rates (LW 11160-11192: uflx(1, k+1) = totuflux(k), hr(1, k+1) = htr(k)) and rrtmg_lwrad's outputs
// (LW 12212-12256) for column i
CPPMPAS_HD void rrtmg_lwrad_finish(const RrtmgRadIn& in, const RrtmgLwradOut& out, int i, int nlayers,
                                   const rrtmg::LwViews& v, rb heatfac) {
  const int kts = 1, kte = in.kte;
  rrtmg::rrtmg_lw_heating(nlayers, v.pz, v.fnet, v.fnetc, v.htr, v.htrc, heatfac);
  out.glw(i) = v.totdflux[0];
  out.olr(i) = v.totuflux[nlayers];
  out.lwcf(i) = v.totuclfl[nlayers] - v.totuflux[nlayers];
  out.lwupt(i) = v.totuflux[nlayers];
  out.lwuptc(i) = v.totuclfl[nlayers];
  out.lwdnt(i) = v.totdflux[nlayers];
  out.lwdntc(i) = v.totdclfl[nlayers];
  out.lwupb(i) = v.totuflux[0];
  out.lwupbc(i) = v.totuclfl[0];
  out.lwdnb(i) = v.totdflux[0];
  out.lwdnbc(i) = v.totdclfl[0];
  for (int k = kts; k <= kte; ++k) {  // LW 12241-12244
    const real tten1d = v.htr[k - 1] / 86400.;
    out.rthratenlw(i, k) = tten1d / in.pi3d(i, k);
  }
  if (in.icloud > 0) {  // LW 12247-12255
    for (int k = kts; k <= kte; ++k) {
      if (in.cldfra3d(i, k) > 0.) {
        out.rre_cloud(i, k) = v.relqmc[k - 1];
        out.rre_ice(i, k) = v.reicmc[k - 1];
        out.rre_snow(i, k) = v.resnmc[k - 1];
      }
    }
  }
}

// ---- rrtmg_swrad (SW 10025-10667) -------------------------------------------------------------------------------
constexpr real sw_co2 = 379.e-6, sw_ch4 = 1774.e-9, sw_n2o = 319.e-9, sw_o2 = 0.209488;  // SW 10137-10148

// Column i of rrtmg_swrad's column loop (SW 10194-10549) up to rrtmg_sw's McICA subcolumns; a night column
// (coszrs <= 0) gets its zero outputs (SW 10649-10661) and iscal[0] = 0, which every later stage skips. inatm_sw,
// setcoef_sw and taumol_sw are left to the per-layer stages (iscal[1] = laytrop counted there).
CPPMPAS_HD void rrtmg_swrad_column(const RrtmgRadIn& in, PCView albedo, PCView xlong, real xtime, real gmt, real radt,
                                   real degrad, real declin, const RrtmgSwradOut& out, int i,
                                   const rrtmg::SwViews& v, const RadWork& w, const RrtmgMpasConsts& kc) {
  using physics_constants::gravity;
  const int kts = 1, kte = in.kte;
  const int icld = in.cldovrlp;  // SW 10173
  int inflgsw = 2, iceflgsw = 3;  // SW 10176-10178
  const int liqflgsw = 1;
  const int nlay = (kte - kts + 1) + 1;  // SW 10258
  auto K = [](int k) { return k - 1; };

  // zenith angle (SW 10196-10206): sincos of xxlat and of declin, cos of hrang (the reference's calls)
  bool dorrsw = true;
  const real xt24 = std::fmod(xtime + radt * 0.5, 1440.);
  const real tloctm = gmt + xt24 / 60. + xlong(i) / 15.;
  const real hrang = 15. * (tloctm - 12.) * degrad;
  const real xxlat = in.xlat(i) * degrad;
  real sinlat, coslat, sindec, cosdec;
  libm::sincos(xxlat, &sinlat, &coslat);
  libm::sincos(declin, &sindec, &cosdec);
  const real coszrs = sinlat * sindec + coslat * cosdec * libm::cos(hrang);
  out.coszr(i) = coszrs;
  if (coszrs <= 0.0) dorrsw = false;
  if (!dorrsw) {  // SW 10649-10661
    out.swupt(i) = 0.;
    out.swuptc(i) = 0.;
    out.swdnt(i) = 0.;
    out.swdntc(i) = 0.;
    out.swupb(i) = 0.;
    out.swupbc(i) = 0.;
    out.swdnb(i) = 0.;
    out.swdnbc(i) = 0.;
    v.iscal[0] = 0;
    return;
  }
  const real lat = in.xlat(i);            // SW 10212
  for (int k = kts; k <= kte + 1; ++k) {  // SW 10216-10219
    w.pw1d[K(k)] = in.p8w(i, k) / 100.;
    w.tw1d[K(k)] = in.t8w(i, k);
  }
  for (int k = kts; k <= kte; ++k) {  // SW 10221-10234
    w.p1d[K(k)] = in.p3d(i, k) / 100.;
    w.dz1d[K(k)] = in.dz8w(i, k);
    w.t1d[K(k)] = in.t3d(i, k);
    w.qv1d[K(k)] = std::max(in.qv3d(i, k), 1.e-12);
    w.qc1d[K(k)] = 0.;
    w.qi1d[K(k)] = 0.;
    w.qs1d[K(k)] = 0.;
    w.cldfra1d[K(k)] = 0.;
  }
  if (in.icloud > 0) {  // SW 10237-10246
    for (int k = kts; k <= kte; ++k) {
      w.qc1d[K(k)] = std::max(in.qc3d(i, k), 0.);
      w.qi1d[K(k)] = std::max(in.qi3d(i, k), 0.);
      w.qs1d[K(k)] = std::max(in.qs3d(i, k), 0.);
      w.cldfra1d[K(k)] = in.cldfra3d(i, k);
    }
  }
  for (int k = kts; k <= kte; ++k) {  // SW 10262-10270
    w.reliq[K(k)] = 10.;
    w.reice[K(k)] = 10.;
    w.resnow[K(k)] = 10.;
    w.cliqwp[K(k)] = 0.;
    w.cicewp[K(k)] = 0.;
    w.csnowp[K(k)] = 0.;
  }
  for (int k = 1; k <= nlay; ++k) {
    v.clwpth[K(k)] = v.ciwpth[K(k)] = v.cswpth[K(k)] = 0.;
    v.relqmc[K(k)] = v.reicmc[K(k)] = v.resnmc[K(k)] = v.cldfrac[K(k)] = 0.;
  }
  // input sounding (SW 10325-10354)
  v.pz[K(1)] = w.pw1d[K(1)];
  v.tz[K(1)] = w.tw1d[K(1)];
  for (int k = kts; k <= kte; ++k) {
    v.pz[K(k + 1)] = w.pw1d[K(k + 1)];
    v.tz[K(k + 1)] = w.tw1d[K(k + 1)];
    v.pavel[K(k)] = w.p1d[K(k)];
    w.pdel[K(k)] = v.pz[K(k)] - v.pz[K(k + 1)];
    v.tavel[K(k)] = w.t1d[K(k)];
    v.wkl[0][K(k)] = w.qv1d[K(k)] * amdw;  // h2ovmr
    v.wkl[1][K(k)] = sw_co2;
    v.wkl[6][K(k)] = sw_o2;
    v.wkl[5][K(k)] = sw_ch4;
    v.wkl[3][K(k)] = sw_n2o;
  }
  v.pz[K(kte + 2)] = 1.0e-5;
  v.tz[K(kte + 2)] = v.tz[K(kte + 1)] + 0.0;
  v.pavel[K(kte + 1)] = 0.5 * v.pz[K(kte + 1)];
  v.tavel[K(kte + 1)] = v.tz[K(kte + 1)] + 0.0;
  v.wkl[0][K(kte + 1)] = v.wkl[0][K(kte)];
  v.wkl[1][K(kte + 1)] = v.wkl[1][K(kte)];
  v.wkl[6][K(kte + 1)] = v.wkl[6][K(kte)];
  v.wkl[5][K(kte + 1)] = v.wkl[5][K(kte)];
  v.wkl[3][K(kte + 1)] = v.wkl[3][K(kte)];
  real dzsum = 0.;  // SW 10359-10364
  for (int k = kts; k <= kte; ++k) {
    w.hlay[K(k)] = dzsum + 0.5 * w.dz1d[K(k)];
    dzsum = dzsum + w.dz1d[K(k)];
  }
  w.hlay[K(kte + 1)] = dzsum + 0.5 * w.dz1d[K(kte)];
  // ozone (SW 10367-10381)
  inirad(w.o3mmr, v.pz, kts, kte, kc);
  if (in.o3input == 2) {
    for (int k = 1; k <= in.noznlevels; ++k) w.o3clim1d[K(k)] = in.o3clim(i, k);
    vinterp_ozn_column(
        kte + 1, [&](int k) { return v.pavel[K(k)]; }, [&](int k) { return in.pin(k); }, in.noznlevels,
        [&](int k) { return w.o3clim1d[K(k)]; }, [&](int k) -> rb& { return w.o3mmr[K(k)]; });
    for (int k = kts; k <= kte + 1; ++k) v.wkl[2][K(k)] = w.o3mmr[K(k)];
  } else {
    for (int k = kts; k <= kte + 1; ++k) v.wkl[2][K(k)] = w.o3mmr[K(k)] * amdo;
  }
  // surface albedo, solar constant, zenith angle (SW 10384-10397): dyofyr = 0, adjes = 1 (the solar constant and
  // adjes enter through RrtmgSwCallConfig::adjflux)
  const real asdir = albedo(i), asdif = albedo(i), aldir = albedo(i), aldif = albedo(i);
  const real coszen = coszrs;
  // cloud optical properties (SW 10400-10535)
  for (int k = kts; k <= kte; ++k) v.cldfrac[K(k)] = w.cldfra1d[K(k)];
  const real gravmks = gravity;
  {  // has_reqc = 0 (SW 10434-10443)
    const real landfrac = 2. - in.xland(i);
    const real landm = landfrac;
    const real snowh = 0.001 * in.snow(i);
    const real icefrac = in.xice(i);
    relcalc(kte - kts + 1, v.tavel, landfrac, landm, icefrac, w.reliq, snowh);
  }
  reicalc(kte - kts + 1, v.tavel, w.reice, kc.retab);  // has_reqi = 0 (SW 10463-10471)
  if (iceflgsw == 3) {
    for (int k = kts; k <= kte; ++k) {
      w.reice[K(k)] = w.reice[K(k)] * 1.0315;
      w.reice[K(k)] = std::min(140.0, w.reice[K(k)]);
    }
  }
  for (int k = kts; k <= kte; ++k) w.resnow[K(k)] = 10.;  // has_reqs = 0 (SW 10481-10485)
  for (int k = kts; k <= kte; ++k) {                       // SW 10488-10491
    const real gliqwp = w.qc1d[K(k)] * w.pdel[K(k)] * 100.0 / gravmks * 1000.0;
    w.cliqwp[K(k)] = gliqwp / std::max(0.01, v.cldfrac[K(k)]);
  }
  if (iceflgsw == 3) {  // SW 10494-10499
    for (int k = kts; k <= kte; ++k) {
      const real gicewp = (w.qi1d[K(k)] + w.qs1d[K(k)]) * w.pdel[K(k)] * 100.0 / gravmks * 1000.0;
      w.cicewp[K(k)] = gicewp / std::max(0.01, v.cldfrac[K(k)]);
    }
  }
  for (int k = kts; k <= kte; ++k) {  // SW 10521-10528
    v.clwpth[K(k)] = w.cliqwp[K(k)];
    v.ciwpth[K(k)] = w.cicewp[K(k)];
    v.cswpth[K(k)] = w.csnowp[K(k)];
    v.relqmc[K(k)] = w.reliq[K(k)];
    v.reicmc[K(k)] = w.reice[K(k)];
    v.resnmc[K(k)] = w.resnow[K(k)];
  }
  v.clwpth[K(kte + 1)] = 0.;  // SW 10531-10537
  v.ciwpth[K(kte + 1)] = 0.;
  v.cswpth[K(kte + 1)] = 0.;
  v.relqmc[K(kte + 1)] = 10.;
  v.reicmc[K(kte + 1)] = 10.;
  v.resnmc[K(kte + 1)] = 10.;
  v.cldfrac[K(kte + 1)] = 0.;

  // rrtmg_sw (SW 8857-9665) on CPPWRF's stages: the inputs besides the wrapper's arrays above
  for (int l = 0; l < nlay; ++l) {
    v.wkl[4][l] = rb(0.0);                                   // co
    for (int g = 0; g < 112; ++g) v.cldfmc[l][g] = rb(0.0);  // cldfmcl = 0 (SW 10292)
  }
  v.iscal[2] = inflgsw;
  v.iscal[3] = iceflgsw;
  v.iscal[4] = liqflgsw;
  // McICA (SW 10540-10549): iplon = 1, irng = 0, permuteseed = 1
  int any_cld_subcol = 0;
  mcica_subcol(icld, 1, 112, nlay, v.pavel, v.cldfrac, w.hlay, in.idcor, in.julday, lat,
               [&](int isub, int ilev, bool cld) {
                 v.cldfmc[ilev][isub] = cld ? rb(1.0) : rb(0.0);
                 if (cld) any_cld_subcol = 1;
               });
  v.iscal[5] = any_cld_subcol;
  v.iscal[0] = 1;  // a day column
  v.iscal[1] = 0;  // laytrop, counted by rrtmg_sw_setcoef_stage
  rrtmg::rrtmg_sw_core_scalars(v, coszen, asdir, asdif, aldir, aldif);
}

// rrtmg_sw's heating rates (SW 9558-9566) and rrtmg_swrad's outputs (SW 10598-10647: swuflx(1, k+1) = fu(k),
// swdkdir(1, 1) = fddir(0), swdkdif(1, 1) = fd(0) - fddir(0)) for a day column i
CPPMPAS_HD void rrtmg_swrad_finish(const RrtmgRadIn& in, const RrtmgSwradOut& out, int i, const rrtmg::SwViews& v,
                                   rb heatfac) {
  const int kts = 1, kte = in.kte, nlay = kte + 1;
  for (int l = 0; l < nlay; ++l) {
    const rb zdpgcp = heatfac / v.pdp[l];
    v.swhrc[l] = ((v.cd[l + 1] - v.cu[l + 1]) - (v.cd[l] - v.cu[l])) * zdpgcp;
    v.swhr[l] = ((v.fd[l + 1] - v.fu[l + 1]) - (v.fd[l] - v.fu[l])) * zdpgcp;
  }
  v.swhrc[nlay - 1] = rb(0.);
  v.swhr[nlay - 1] = rb(0.);
  const real coszrs = out.coszr(i);  // the zenith angle of rrtmg_swrad_column (SW 10206)
  out.gsw(i) = v.fd[0] - v.fu[0];
  out.swcf(i) = (v.fd[kte + 1] - v.fu[kte + 1]) - (v.cd[kte + 1] - v.cu[kte + 1]);
  out.swupt(i) = v.fu[kte + 1];
  out.swuptc(i) = v.cu[kte + 1];
  out.swdnt(i) = v.fd[kte + 1];
  out.swdntc(i) = v.cd[kte + 1];
  out.swupb(i) = v.fu[0];
  out.swupbc(i) = v.cu[0];
  out.swdnb(i) = v.fd[0];
  out.swdnbc(i) = v.cd[0];
  out.swddir(i) = v.fddir[0];
  out.swddni(i) = out.swddir(i) / coszrs;
  out.swddif(i) = v.fd[0] - v.fddir[0];
  for (int k = kts; k <= kte; ++k) {  // SW 10640-10644
    const real tten1d = v.swhr[k - 1] / 86400.;
    out.rthratensw(i, k) = tten1d / in.pi3d(i, k);
  }
}

void check_effective_radii(const RrtmgRadIn& in, const char* who) {
  if (in.has_reqc != 0 || in.has_reqi != 0 || in.has_reqs != 0)
    throw Error(std::string(who) + ": the microphysics' effective radii (has_reqc / has_reqi / has_reqs = 1, "
                                   "config_microp_re = true) have no verified port");
  if (in.cldovrlp == 4)  // mcica_subcol, called for every column (LW 12181, SW 10540)
    throw Error("config_radt_cld_overlap = 'exponential' has no verified port (MPAS's -O3 runs its alpha loop in "
                "_ZGVbN2v_exp lanes, module_ra_rrtmg_lw.F:2552)");
  if (in.o3input == 2 && in.noznlevels > kMaxOznLevels)
    throw Error(std::string(who) + ": more than 64 ozone levels");
}

}  // namespace

void rrtmg_mpas_fill_consts(RrtmgMpasConsts& k) {
  for (int n = 0; n < 95; ++n) k.retab[n] = retab_data[n];
  for (int n = 0; n < 60; ++n) {
    k.pprof[n] = pprof_data[n];
    k.tprof[n] = tprof_data[n];
  }
  for (int n = 0; n < 31; ++n) {
    k.o3sum[n] = o3sum_data[n];
    k.ppsum[n] = ppsum_data[n];
    k.o3win[n] = o3win_data[n];
    k.ppwin[n] = ppwin_data[n];
  }
}

// ---- rrtmg_lwrad: the column loop as CPPWRF's LW batch (rrtmg_lw_batch.cpp RrtmgLwBatch::run) --------------------
void rrtmg_lwrad(const RrtmgRadIn& in, PCView emiss, const RrtmgLwradOut& out, const rrtmg::RrtmgLwTables* tables,
                 const RrtmgMpasConsts* kc, int its, int ite) {
  using physics_constants::cp;
  if (its != 1) throw Error("rrtmg_lwrad: its /= 1");
  check_effective_radii(in, "rrtmg_lwrad");
  if (ite < its) return;
  const int ncol = ite - its + 1;
  const rb heatfac = rb(9.8066) * rb(8.6400e4) / (cp * rb(1.e2));  // rrtmg_lw_ini: lwdatinit (LW 8187, 8207, 8244)
  const RrtmgRadIn a = in;
  // the layers of each column (LW 11822)
  FStorage<int> nlay_s(static_cast<std::size_t>(ncol), 0);
  int* const nlay = nlay_s.data();
  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {
    const int nl = lw_nlayers(a, i);
    if (nl > rrtmg::RRTMG_LW_SWEEP_MAX) kernel_fail("rrtmg_lwrad: more than 128 layers");
    nlay[i - its] = nl;
  });
  int nlmax = 0;
  for (int c = 0; c < ncol; ++c) nlmax = std::max(nlmax, nlay[c]);
  const int chunk = std::min(ncol, kChunk);
  const int noz = in.o3input == 2 ? in.noznlevels : 0;
  const std::size_t nd = rrtmg::LwLayout(nlmax, false).nd, ni = rrtmg::LwLayout(nlmax, false).ni;
  const std::size_t nw = RadWorkLayout(in.kte, nlmax, noz).nd;
  FStorage<rb> dbuf(nd * std::size_t(chunk)), wbuf(nw * std::size_t(chunk));
  FStorage<int> ibuf(ni * std::size_t(chunk));
  rb* const d = dbuf.data();
  rb* const wd = wbuf.data();
  int* const ip = ibuf.data();
  for (int c0 = 0; c0 < ncol; c0 += chunk) {
    const int n = std::min(chunk, ncol - c0);
    int nlc = 0;  // the chunk's largest number of layers: its layer dimension
    for (int c = c0; c < c0 + n; ++c) nlc = std::max(nlc, nlay[c]);
    const rrtmg::LwLayout lay(nlc, false);
    const RadWorkLayout wl(in.kte, nlc, noz);
    const std::size_t sn = std::size_t(n);
    // per column: the wrapper (LW 11766-12191), inatm, secdiff
    parallel_for_n(sn, [=] CPPMPAS_LAMBDA(std::size_t x) {
      const int ci = int(x), c = c0 + ci;
      rrtmg_lwrad_column(a, emiss, its + c, nlay[c], lay.views(d, ip, ci, n), wl.views(wd, ci, n), *kc);
    });
    // per (column, layer): setcoef (laytrop counted), then taumol and cldprmc
    parallel_for_n(sn * std::size_t(nlc), [=] CPPMPAS_LAMBDA(std::size_t x) {
      const int ci = int(x % sn), l = int(x / sn), c = c0 + ci;
      if (l >= nlay[c]) return;
      rrtmg::rrtmg_lw_setcoef_stage(l, a.tsk(its + c), lay.views(d, ip, ci, n), *tables);
    });
    parallel_for_n(sn * std::size_t(nlc), [=] CPPMPAS_LAMBDA(std::size_t x) {
      const int ci = int(x % sn), l = int(x / sn), c = c0 + ci;
      if (l >= nlay[c]) return;
      rrtmg::rrtmg_lw_layer(nlay[c], l, lay.views(d, ip, ci, n), *tables);
    });
    // per (column, g-point): rtrnmc's sweeps; per (column, level): its band sums and fluxes
    if (nlc <= 60)
      parallel_for_n(sn * 140, [=] CPPMPAS_LAMBDA(std::size_t x) {
        const int ci = int(x % sn), g = int(x / sn), c = c0 + ci;
        rrtmg::rrtmg_lw_gpointstage<60>(nlay[c], g, lay.views(d, ip, ci, n), *tables);
      });
    else
      parallel_for_n(sn * 140, [=] CPPMPAS_LAMBDA(std::size_t x) {
        const int ci = int(x % sn), g = int(x / sn), c = c0 + ci;
        rrtmg::rrtmg_lw_gpointstage<rrtmg::RRTMG_LW_SWEEP_MAX>(nlay[c], g, lay.views(d, ip, ci, n), *tables);
      });
    parallel_for_n(sn * std::size_t(nlc + 1), [=] CPPMPAS_LAMBDA(std::size_t x) {
      const int ci = int(x % sn), lev = int(x / sn), c = c0 + ci;
      if (lev > nlay[c]) return;
      rrtmg::rrtmg_lw_bandsum_level(nlay[c], lev, lay.views(d, ip, ci, n), *tables);
    });
    // per column: heating rates and the outputs
    parallel_for_n(sn, [=] CPPMPAS_LAMBDA(std::size_t x) {
      const int ci = int(x), c = c0 + ci;
      rrtmg_lwrad_finish(a, out, its + c, nlay[c], lay.views(d, ip, ci, n), heatfac);
    });
  }
}

// ---- rrtmg_swrad: the column loop as CPPWRF's SW batch (rrtmg_sw_batch.cpp RrtmgSwBatch::run) --------------------
void rrtmg_swrad(const RrtmgRadIn& in, PCView albedo, PCView xlong, real xtime, real gmt, real radt, real degrad,
                 real declin, real solcon, const RrtmgSwradOut& out, const rrtmg::RrtmgSwTables* tables,
                 const RrtmgMpasConsts* kc, int its, int ite) {
  using physics_constants::cp;
  if (its != 1) throw Error("rrtmg_swrad: its /= 1");
  check_effective_radii(in, "rrtmg_swrad");
  const int kts = 1, kte = in.kte;
  const int nlay = (kte - kts + 1) + 1;  // SW 10258
  if (nlay > rrtmg::RRTMG_SW_COL_MAX) throw Error("rrtmg_swrad: more than 128 layers");
  if (ite < its) return;
  const int ncol = ite - its + 1;
  // rrsw_con heatfac (rrtmg_sw_ini: swdatinit, SW 4856, 4876, 4912); rrsw_scon (SW 123)
  const rb heatfac = rb(9.8066) * rb(8.6400e4) / (cp * rb(1.e2));
  const rb rrsw_scon = rb(1.36822e+03);
  const real scon = solcon;  // SW 10396-10397
  const real adjes = 1.0;
  rrtmg::RrtmgSwCallConfig cfg{};
  cfg.num_levels = nlay - 1;
  cfg.nlayers = nlay;
  cfg.heatfac = heatfac;
  cfg.adjflux = adjes * (scon / rrsw_scon);  // inatm_sw (SW 9909-9916)
  const RrtmgRadIn a = in;
  const int chunk = std::min(ncol, kChunk);
  const int noz = in.o3input == 2 ? in.noznlevels : 0;
  const rrtmg::SwLayout lay(nlay);
  const RadWorkLayout wl(kte, nlay, noz);
  FStorage<rb> dbuf(lay.nd * std::size_t(chunk)), wbuf(wl.nd * std::size_t(chunk));
  FStorage<int> ibuf(lay.ni * std::size_t(chunk));
  rb* const d = dbuf.data();
  rb* const wd = wbuf.data();
  int* const ip = ibuf.data();
  for (int c0 = 0; c0 < ncol; c0 += chunk) {
    const int n = std::min(chunk, ncol - c0);
    const std::size_t sn = std::size_t(n);
    // per column: the zenith angle, night columns, the wrapper (SW 10194-10549)
    parallel_for_n(sn, [=] CPPMPAS_LAMBDA(std::size_t x) {
      const int ci = int(x);
      rrtmg_swrad_column(a, albedo, xlong, xtime, gmt, radt, degrad, declin, out, its + c0 + ci,
                         lay.views(d, ip, ci, n), wl.views(wd, ci, n), *kc);
    });
    // per (column, layer): inatm_sw and setcoef_sw; per column: taumol_sw's solar source; per (column, layer): taumol_sw
    parallel_for_n(sn * std::size_t(nlay), [=] CPPMPAS_LAMBDA(std::size_t x) {
      rrtmg::rrtmg_sw_setcoef_stage(int(x / sn), lay.views(d, ip, int(x % sn), n), *tables);
    });
    parallel_for_n(sn, [=] CPPMPAS_LAMBDA(std::size_t x) {
      rrtmg::rrtmg_sw_solfr_stage(nlay, lay.views(d, ip, int(x), n), *tables);
    });
    parallel_for_n(sn * std::size_t(nlay), [=] CPPMPAS_LAMBDA(std::size_t x) {
      rrtmg::rrtmg_sw_layer(nlay, int(x / sn), lay.views(d, ip, int(x % sn), n), *tables);
    });
    // per (column, g-point): spcvmc_sw; per (column, level): its sums
    if (nlay <= 48)
      parallel_for_n(sn * 112, [=] CPPMPAS_LAMBDA(std::size_t x) {
        rrtmg::rrtmg_sw_gpoint<48>(cfg, int(x / sn), lay.views(d, ip, int(x % sn), n), *tables);
      });
    else
      parallel_for_n(sn * 112, [=] CPPMPAS_LAMBDA(std::size_t x) {
        rrtmg::rrtmg_sw_gpoint<rrtmg::RRTMG_SW_COL_MAX>(cfg, int(x / sn), lay.views(d, ip, int(x % sn), n), *tables);
      });
    parallel_for_n(sn * std::size_t(nlay + 1), [=] CPPMPAS_LAMBDA(std::size_t x) {
      rrtmg::rrtmg_sw_level(int(x / sn), lay.views(d, ip, int(x % sn), n), *tables);
    });
    // per day column: heating rates and the outputs
    parallel_for_n(sn, [=] CPPMPAS_LAMBDA(std::size_t x) {
      const int ci = int(x);
      const rrtmg::SwViews v = lay.views(d, ip, ci, n);
      if (v.iscal[0] == 0) return;
      rrtmg_swrad_finish(a, out, its + c0 + ci, v, heatfac);
    });
  }
}

}  // namespace cppmpas::phys
