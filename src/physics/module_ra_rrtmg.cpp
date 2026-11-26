// SPDX-License-Identifier: GPL-3.0-or-later AND LicenseRef-AER-RRTMG
// See module_ra_rrtmg.hpp. LW = physics/physics_wrf/module_ra_rrtmg_lw.F, SW = module_ra_rrtmg_sw.F, VI =
// module_ra_rrtmg_vinterp.F (MPAS v8.4.2). With -fdefault-real-8 every REAL and literal is a double.
#include "physics/module_ra_rrtmg.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

#include "core/error.hpp"
#include "core/libm.hpp"
#include "physics/mpas_atmphys_constants.hpp"
#include "physics/rrtmg/kissvec.hpp"
#include "physics/rrtmg/rrtmg_lw_stages.hpp"
#include "physics/rrtmg/rrtmg_sw_stages.hpp"

namespace cppmpas::phys {

namespace {

using rrtmg::rb;

// retab (LW 11547-11566): re(T) of hexagonal columns, 180 K .. 274 K
constexpr real retab[95] = {
    5.92779, 6.26422, 6.61973, 6.99539, 7.39234, 7.81177, 8.25496, 8.72323, 9.21800, 9.74075, 10.2930,
    10.8765, 11.4929, 12.1440, 12.8317, 13.5581, 14.2319, 15.0351, 15.8799, 16.7674, 17.6986, 18.6744,
    19.6955, 20.7623, 21.8757, 23.0364, 24.2452, 25.5034, 26.8125, 27.7895, 28.6450, 29.4167, 30.1088,
    30.7306, 31.2943, 31.8151, 32.3077, 32.7870, 33.2657, 33.7540, 34.2601, 34.7892, 35.3442, 35.9255,
    36.5316, 37.1602, 37.8078, 38.4720, 39.1508, 39.8442, 40.5552, 41.2912, 42.0635, 42.8876, 43.7863,
    44.7853, 45.9170, 47.2165, 48.7221, 50.4710, 52.4980, 54.8315, 57.4898, 60.4785, 63.7898, 65.5604,
    71.2885, 75.4113, 79.7368, 84.2351, 88.8833, 93.6658, 98.5739, 103.603, 108.752, 114.025, 119.424,
    124.954, 130.630, 136.457, 142.446, 148.608, 154.956, 161.503, 168.262, 175.248, 182.473, 189.952,
    197.699, 205.728, 214.055, 222.694, 231.661, 240.971, 250.639};

// relcalc (LW 12447-12497) for one column (ncol = pcols = 1): rel(1:pver) from t(1:pver)
void relcalc(int pver, const real* t, real landfrac, real landm, real icefrac, real* rel, real snowh) {
  (void)landfrac;
  const real tmelt = 273.16, rliqocean = 14.0, rliqice = 14.0, rliqland = 8.0;  // LW 12483-12486
  for (int k = 1; k <= pver; ++k) {
    real& r = rel[k - 1];
    const real tk = t[k - 1];
    r = rliqland + (rliqocean - rliqland) * std::min(1.0, std::max(0.0, (tmelt - tk) * 0.05));  // LW 12491
    r = r + (rliqocean - r) * std::min(1.0, std::max(0.0, snowh * 10.));                       // LW 12493
    r = r + (rliqocean - r) * std::min(1.0, std::max(0.0, 1.0 - landm));                       // LW 12495
    r = r + (rliqice - r) * std::min(1.0, std::max(0.0, icefrac));                             // LW 12497
  }
}

// reicalc (LW 12499-12527) for one column
void reicalc(int pver, const real* t, real* re) {
  for (int k = 1; k <= pver; ++k) {
    const real tk = t[k - 1];
    int index = static_cast<int>(tk - 179.);       // LW 12519 (INT truncates)
    index = std::min(std::max(index, 1), 94);      // LW 12520
    const real corr = tk - static_cast<int>(tk);   // LW 12521
    re[k - 1] = retab[index - 1] * (1. - corr) + retab[index] * corr;  // LW 12522-12523
  }
}

// mcica_subcol_lw / _sw (LW 2099-2248, SW 1430-1555) with generate_stochastic_clouds (irng = 0): CPPWRF's
// KissRng and mcica_subcolumns (physics/rrtmg/kissvec.hpp). play in mb (pmid = play*1.e2, LW 2201); cldfrac(1:nlay);
// cldy(isubcol, ilev, cloudy) receives the subcolumns of the layers with cldf > 0 (for the others the Fortran's
// cldfmcl stays the 0 rrtmg_lwrad / rrtmg_swrad set). Returns false (no call of cldy) for icld = 0 (LW 2194).
template <class OUT>
bool mcica_subcol(int icld, int permuteseed, int nsubcol, int nlay, const real* play, const real* cldfrac,
                  const real* hgt, int idcor, int juldat, real lat, OUT cldy) {
  if (icld == 0) return false;  // LW 2194
  if (icld < 0 || icld > 5) throw Error("mcica_subcol: MCICA_SUBCOL: INVALID ICLD");
  if (icld == 4)
    throw Error("config_radt_cld_overlap = 'exponential' has no verified port (MPAS's -O3 runs its alpha loop "
                "in _ZGVbN2v_exp lanes, module_ra_rrtmg_lw.F:2552)");
  if (nlay > rrtmg::kMcicaMaxLayers) throw Error("mcica_subcol: more than 128 layers");
  // generate_stochastic_clouds: KISS seed from the bottom four layers (LW 2441-2452)
  if (play[0] * 1.e2 < play[1] * 1.e2)
    throw Error("MCICA_SUBCOL: KISSVEC SEED GENERATOR REQUIRES PMID FROM BOTTOM FOUR LAYERS.");
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

}  // namespace

// ---- vinterp_ozn (VI:16-176) ------------------------------------------------------------------------------------
void vinterp_ozn(int lchnk, int ncol, int pcols, int pver, const real* pmid_in, const real* pin_in, int levsiz,
                 const real* ozmix_in, real* o3vmr_out) {
  (void)lchnk;
  using Z = std::size_t;
  std::vector<real> pin(static_cast<Z>(levsiz)), ozmix(Z(pcols) * Z(levsiz)), pmid(Z(pcols) * Z(pver)),
      o3vmr(Z(pcols) * Z(pver));
  std::vector<int> kupper(static_cast<Z>(pcols));
  auto O = [&](int i, int k) -> real& { return ozmix[Z(i - 1) + Z(pcols) * Z(k - 1)]; };
  auto P = [&](int i, int k) -> real& { return pmid[Z(i - 1) + Z(pcols) * Z(k - 1)]; };
  auto V = [&](int i, int k) -> real& { return o3vmr[Z(i - 1) + Z(pcols) * Z(k - 1)]; };
  auto PIN = [&](int k) -> real& { return pin[Z(k - 1)]; };
  for (int k = 1; k <= levsiz; ++k) PIN(k) = pin_in[k - 1];  // VI:73-75
  for (int i = 1; i <= pcols; ++i)                           // VI:76-80
    for (int k = 1; k <= levsiz; ++k) O(i, k) = ozmix_in[Z(i - 1) + Z(pcols) * Z(k - 1)];
  for (int i = 1; i <= pcols; ++i)  // VI:81-86
    for (int k = 1; k <= pver; ++k) {
      const int kk = pver - k + 1;
      P(i, kk) = pmid_in[Z(i - 1) + Z(pcols) * Z(k - 1)];
    }
  for (int i = 1; i <= ncol; ++i) kupper[Z(i - 1)] = 1;  // VI:91-93
  for (int k = 1; k <= pver; ++k) {                       // VI:95
    int kkstart = levsiz;                                 // VI:100-103
    for (int i = 1; i <= ncol; ++i) kkstart = std::min(kkstart, kupper[Z(i - 1)]);
    int kount = 0;
    bool found = false;
    for (int kk = kkstart; kk <= levsiz - 1; ++kk) {  // VI:109-131
      for (int i = 1; i <= ncol; ++i) {
        if (PIN(kk) < P(i, k) && P(i, k) <= PIN(kk + 1)) {
          kupper[Z(i - 1)] = kk;
          kount = kount + 1;
        }
      }
      if (kount == ncol) {  // VI:120-130
        for (int i = 1; i <= ncol; ++i) {
          const int ku = kupper[Z(i - 1)];
          const real dpu = P(i, k) - PIN(ku);
          const real dpl = PIN(ku + 1) - P(i, k);
          V(i, k) = (O(i, ku) * dpl + O(i, ku + 1) * dpu) / (dpl + dpu);
        }
        found = true;
        break;  // goto 35
      }
    }
    if (found) continue;
    for (int i = 1; i <= ncol; ++i) {  // VI:137-148
      if (P(i, k) < PIN(1)) {
        V(i, k) = O(i, 1) * P(i, k) / PIN(1);
      } else if (P(i, k) > PIN(levsiz)) {
        V(i, k) = O(i, levsiz);
      } else {
        const int ku = kupper[Z(i - 1)];
        const real dpu = P(i, k) - PIN(ku);
        const real dpl = PIN(ku + 1) - P(i, k);
        V(i, k) = (O(i, ku) * dpl + O(i, ku + 1) * dpu) / (dpl + dpu);
      }
    }
  }
  for (int i = 1; i <= pcols; ++i)  // VI:160-165
    for (int k = 1; k <= pver; ++k) {
      const int kk = pver - k + 1;
      o3vmr_out[Z(i - 1) + Z(pcols) * Z(kk - 1)] = V(i, k);
    }
}

// ---- inirad / o3data (LW 12308-12444) -------------------------------------------------------------------------
void inirad(real* o3prof, const real* plev, int kts, int kte) {
  // o3prof(kts:kte+1), plev(kts:kte+2): element k at [k - kts]
  auto O3 = [&](int k) -> real& { return o3prof[k - kts]; };
  auto PL = [&](int k) { return plev[k - kts]; };
  for (int k = kts; k <= kte + 1; ++k) O3(k) = 0.;  // LW 12326-12328
  // O3DATA (LW 12334-12444)
  static constexpr real o3sum[31] = {
      5.297E-8, 5.852E-8, 6.579E-8, 7.505E-8, 8.577E-8, 9.895E-8, 1.175E-7, 1.399E-7, 1.677E-7, 2.003E-7, 2.571E-7,
      3.325E-7, 4.438E-7, 6.255E-7, 8.168E-7, 1.036E-6, 1.366E-6, 1.855E-6, 2.514E-6, 3.240E-6, 4.033E-6, 4.854E-6,
      5.517E-6, 6.089E-6, 6.689E-6, 1.106E-5, 1.462E-5, 1.321E-5, 9.856E-6, 5.960E-6, 5.960E-6};
  static constexpr real ppsum[31] = {955.890, 850.532, 754.599, 667.742, 589.841, 519.421, 455.480, 398.085,
                                     347.171, 301.735, 261.310, 225.360, 193.419, 165.490, 141.032, 120.125,
                                     102.689, 87.829,  75.123,  64.306,  55.086,  47.209,  40.535,  34.795,
                                     29.865,  19.122,  9.277,   4.660,   2.421,   1.294,   0.647};
  static constexpr real o3win[31] = {
      4.629E-8, 4.686E-8, 5.017E-8, 5.613E-8, 6.871E-8, 8.751E-8, 1.138E-7, 1.516E-7, 2.161E-7, 3.264E-7, 4.968E-7,
      7.338E-7, 1.017E-6, 1.308E-6, 1.625E-6, 2.011E-6, 2.516E-6, 3.130E-6, 3.840E-6, 4.703E-6, 5.486E-6, 6.289E-6,
      6.993E-6, 7.494E-6, 8.197E-6, 9.632E-6, 1.113E-5, 1.146E-5, 9.389E-6, 6.135E-6, 6.135E-6};
  static constexpr real ppwin[31] = {955.747, 841.783, 740.199, 649.538, 568.404, 495.815, 431.069, 373.464,
                                     322.354, 277.190, 237.635, 203.433, 174.070, 148.949, 127.408, 108.915,
                                     93.114,  79.551,  67.940,  58.072,  49.593,  42.318,  36.138,  30.907,
                                     26.362,  16.423,  7.583,   3.620,   1.807,   0.938,   0.469};
  real ppann[31], o3ann[31], o3wrk[31], ppwrk[31], ppwrkh[32];
  for (int k = 0; k < 31; ++k) ppann[k] = ppsum[k];  // LW 12386-12388
  o3ann[0] = 0.5 * (o3sum[0] + o3win[0]);           // LW 12390
  for (int k = 1; k < 31; ++k)                       // LW 12392-12395
    o3ann[k] = o3win[k - 1] + (o3win[k] - o3win[k - 1]) / (ppwin[k] - ppwin[k - 1]) * (ppsum[k] - ppwin[k - 1]);
  for (int k = 1; k < 31; ++k) o3ann[k] = 0.5 * (o3ann[k] + o3sum[k]);  // LW 12397-12399
  for (int k = 0; k < 31; ++k) {                                         // LW 12401-12404
    o3wrk[k] = o3ann[k];
    ppwrk[k] = ppann[k];
  }
  std::vector<real> prlevh(static_cast<std::size_t>(kte + 2 - kts + 1));
  for (int k = kts; k <= kte + 2; ++k) prlevh[k - kts] = PL(k);  // LW 12412-12414
  auto PR = [&](int k) { return prlevh[k - kts]; };
  ppwrkh[0] = 1100.;                                                   // LW 12416
  for (int k = 1; k < 31; ++k) ppwrkh[k] = (ppwrk[k] + ppwrk[k - 1]) / 2.;  // LW 12417-12419
  ppwrkh[31] = 0.;                                                     // LW 12420
  for (int k = kts; k <= kte + 1; ++k) {                               // LW 12421-12442
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

// ---- rrtmg_lwrad (LW 11574-12305) -------------------------------------------------------------------------------
void rrtmg_lwrad(const RrtmgRadIn& in, const std::vector<real>& emiss, RrtmgLwradOut& out,
                 const rrtmg::RrtmgLwTables& tables, int its, int ite) {
  using physics_constants::cp;
  using physics_constants::gravity;
  if (its != 1) throw Error("rrtmg_lwrad: its /= 1");
  const int kts = 1, kte = in.kte;
  if (in.has_reqc != 0 || in.has_reqi != 0 || in.has_reqs != 0)
    throw Error("rrtmg_lwrad: the microphysics' effective radii (has_reqc / has_reqi / has_reqs = 1, "
                "config_microp_re = true) have no verified port");
  // trace gases (LW 11681-11704), molecular weight ratios (LW 11718-11720), buffer layer (LW 11569)
  const real co2 = 379.e-6, ch4 = 1774.e-9, n2o = 319.e-9, cfc11 = 0.251e-9, cfc12 = 0.538e-9, cfc22 = 0.169e-9,
             ccl4 = 0.093e-9, o2 = 0.209488;
  const real amdw = 1.607793, amdo = 0.603461;
  const real deltap = 4.;
  // standard profile above the model top (LW 11726-11749)
  constexpr int nproflevs = 60;
  static constexpr real pprof[nproflevs] = {
      1000.00, 855.47, 731.82, 626.05, 535.57, 458.16, 391.94, 335.29, 286.83, 245.38, 209.91, 179.57,
      153.62,  131.41, 112.42, 96.17,  82.27,  70.38,  60.21,  51.51,  44.06,  37.69,  32.25,  27.59,
      23.60,   20.19,  17.27,  14.77,  12.64,  10.81,  9.25,   7.91,   6.77,   5.79,   4.95,   4.24,
      3.63,    3.10,   2.65,   2.27,   1.94,   1.66,   1.42,   1.22,   1.04,   0.89,   0.76,   0.65,
      0.56,    0.48,   0.41,   0.35,   0.30,   0.26,   0.22,   0.19,   0.16,   0.14,   0.12,   0.10};
  static constexpr real tprof[nproflevs] = {
      286.96, 281.07, 275.16, 268.11, 260.56, 253.02, 245.62, 238.41, 231.57, 225.91, 221.72, 217.79,
      215.06, 212.74, 210.25, 210.16, 210.69, 212.14, 213.74, 215.37, 216.82, 217.94, 219.03, 220.18,
      221.37, 222.64, 224.16, 225.88, 227.63, 229.51, 231.50, 233.73, 236.18, 238.78, 241.60, 244.44,
      247.35, 250.33, 253.32, 256.30, 259.22, 262.12, 264.80, 266.50, 267.59, 268.44, 268.69, 267.76,
      266.13, 263.96, 261.54, 258.93, 256.15, 253.23, 249.89, 246.67, 243.48, 240.25, 236.66, 233.86};

  const int icld = in.cldovrlp;  // LW 11754
  int inflglw = 2, iceflglw = 3;  // LW 11757-11759
  const int liqflglw = 1;
  const rb heatfac = rb(9.8066) * rb(8.6400e4) / (cp * rb(1.e2));  // rrtmg_lw_ini: lwdatinit (LW 8187, 8207, 8244)

  const int nlev = kte - kts + 1;
  std::vector<real> pw1d(nlev + 1), tw1d(nlev + 1), p1d(nlev), dz1d(nlev), t1d(nlev), qv1d(nlev), qc1d(nlev),
      qi1d(nlev), qs1d(nlev), cldfra1d(nlev), pdel(nlev), reliq(nlev), reice(nlev), resnow(nlev), cliqwp(nlev),
      cicewp(nlev), csnowp(nlev);
  std::vector<rb> d;
  std::vector<int> iw;
  auto K = [](int k) { return static_cast<std::size_t>(k - 1); };  // 1-based local arrays

  for (int i = its; i <= ite; ++i) {  // LW 11766-11768
    const real lat = (*in.xlat)[K(i)];  // LW 11771
    // column sounding (LW 11782-11800)
    for (int k = kts; k <= kte + 1; ++k) {
      pw1d[K(k)] = (*in.p8w)(i, k) / 100.;
      tw1d[K(k)] = (*in.t8w)(i, k);
    }
    for (int k = kts; k <= kte; ++k) {
      p1d[K(k)] = (*in.p3d)(i, k) / 100.;
      dz1d[K(k)] = (*in.dz8w)(i, k);
      t1d[K(k)] = (*in.t3d)(i, k);
      qv1d[K(k)] = std::max((*in.qv3d)(i, k), 1.e-12);
      qc1d[K(k)] = 0.;
      qi1d[K(k)] = 0.;
      qs1d[K(k)] = 0.;
      cldfra1d[K(k)] = 0.;
    }
    if (in.icloud > 0) {  // LW 11803-11812 (qr3d is not passed: qr1d stays 0 and is not used)
      for (int k = kts; k <= kte; ++k) {
        qc1d[K(k)] = std::max((*in.qc3d)(i, k), 0.);
        qi1d[K(k)] = std::max((*in.qi3d)(i, k), 0.);
        qs1d[K(k)] = std::max((*in.qs3d)(i, k), 0.);
        cldfra1d[K(k)] = (*in.cldfra3d)(i, k);
      }
    }
    // buffer layers (LW 11822): NINT is lround
    const int nlayers = kte + std::max(static_cast<int>(std::lround(pw1d[K(kte + 1)] / deltap)), 1);
    if (nlayers > rrtmg::RRTMG_LW_SWEEP_MAX) throw Error("rrtmg_lwrad: more than 128 layers");
    // column arrays (LW 11827-11868), 1-based: plev, tlev (nlayers + 1), the layer arrays (nlayers)
    std::vector<real> plev(nlayers + 1), tlev(nlayers + 1), play(nlayers), hlay(nlayers), tlay(nlayers),
        h2ovmr(nlayers), o3vmr(nlayers), co2vmr(nlayers), o2vmr(nlayers), ch4vmr(nlayers), n2ovmr(nlayers),
        cfc11vmr(nlayers), cfc12vmr(nlayers), cfc22vmr(nlayers), ccl4vmr(nlayers), clwpth(nlayers),
        ciwpth(nlayers), cswpth(nlayers), rel(nlayers), rei(nlayers), res(nlayers), cldfrac(nlayers),
        o3mmr(nlayers), varint(nlayers + 1);
    // LW 11872-11914: local arrays (reliq.. on kts:kte, the rest 0)
    for (int k = kts; k <= kte; ++k) {
      reliq[K(k)] = 10.;
      reice[K(k)] = 10.;
      resnow[K(k)] = 10.;
      cliqwp[K(k)] = 0.;
      cicewp[K(k)] = 0.;
      csnowp[K(k)] = 0.;
    }
    // input sounding (LW 11917-11937)
    plev[K(1)] = pw1d[K(1)];
    tlev[K(1)] = tw1d[K(1)];
    const real tsfc = (*in.tsk)[K(i)];
    for (int k = kts; k <= kte; ++k) {
      plev[K(k + 1)] = pw1d[K(k + 1)];
      tlev[K(k + 1)] = tw1d[K(k + 1)];
      play[K(k)] = p1d[K(k)];
      pdel[K(k)] = plev[K(k)] - plev[K(k + 1)];
      tlay[K(k)] = t1d[K(k)];
      h2ovmr[K(k)] = qv1d[K(k)] * amdw;
      co2vmr[K(k)] = co2;
      o2vmr[K(k)] = o2;
      ch4vmr[K(k)] = ch4;
      n2ovmr[K(k)] = n2o;
      cfc11vmr[K(k)] = cfc11;
      cfc12vmr[K(k)] = cfc12;
      cfc22vmr[K(k)] = cfc22;
      ccl4vmr[K(k)] = ccl4;
    }
    // layer mid-point heights (LW 11942-11950)
    real dzsum = 0.;
    for (int k = kts; k <= kte; ++k) {
      hlay[K(k)] = dzsum + 0.5 * dz1d[K(k)];
      dzsum = dzsum + dz1d[K(k)];
    }
    for (int k = kte + 1; k <= nlayers; ++k) {
      hlay[K(k)] = dzsum + 0.5 * dz1d[K(kte)];
      dzsum = dzsum + dz1d[K(kte)];
    }
    // buffer pressures (LW 11955-11964)
    for (int k = kte + 1; k <= nlayers; ++k) {
      plev[K(k + 1)] = std::max(plev[K(k)] - deltap, 0.00);
      play[K(k)] = 0.5 * (plev[K(k)] + plev[K(k + 1)]);
    }
    plev[K(nlayers + 1)] = 0.00;
    play[K(nlayers)] = 0.5 * (plev[K(nlayers)] + plev[K(nlayers + 1)]);
    // table temperatures at the column pressures (LW 11967-11989)
    for (int k = 1; k <= nlayers + 1; ++k) {
      int klev = nproflevs;
      if (pprof[nproflevs - 1] < plev[K(k)]) {
        for (int kk = 2; kk <= nproflevs; ++kk) {
          if (pprof[kk - 1] < plev[K(k)]) {
            klev = kk - 1;
            break;
          }
        }
      } else {
        klev = nproflevs;
      }
      real vark, vark1, wght;
      if (klev != nproflevs) {
        vark = tprof[klev - 1];
        vark1 = tprof[klev];
        wght = (plev[K(k)] - pprof[klev - 1]) / (pprof[klev] - pprof[klev - 1]);
      } else {
        vark = tprof[klev - 1];
        vark1 = tprof[klev - 1];
        wght = 0.0;
      }
      varint[K(k)] = wght * (vark1 - vark) + vark;
    }
    for (int k = kte + 1; k <= nlayers + 1; ++k) {  // LW 11992-11997
      tlev[K(k)] = varint[K(k)] + (tlev[K(kte)] - varint[K(kte)]);
      tlay[K(k - 1)] = 0.5 * (tlev[K(k)] + tlev[K(k - 1)]);
    }
    for (int k = kte + 1; k <= nlayers; ++k) {  // LW 12000-12010
      h2ovmr[K(k)] = h2ovmr[K(kte)];
      co2vmr[K(k)] = co2vmr[K(kte)];
      o2vmr[K(k)] = o2vmr[K(kte)];
      ch4vmr[K(k)] = ch4vmr[K(kte)];
      n2ovmr[K(k)] = n2ovmr[K(kte)];
      cfc11vmr[K(k)] = cfc11vmr[K(kte)];
      cfc12vmr[K(k)] = cfc12vmr[K(kte)];
      cfc22vmr[K(k)] = cfc22vmr[K(kte)];
      ccl4vmr[K(k)] = ccl4vmr[K(kte)];
    }
    // ozone (LW 12013-12027)
    inirad(o3mmr.data(), plev.data(), kts, nlayers - 1);
    if (in.o3input == 2) {
      std::vector<real> o3clim1d(static_cast<std::size_t>(in.noznlevels));
      for (int k = 1; k <= in.noznlevels; ++k) o3clim1d[K(k)] = (*in.o3clim)(i, k);
      vinterp_ozn(1, 1, 1, nlayers, play.data(), in.pin, in.noznlevels, o3clim1d.data(), o3mmr.data());
      for (int k = kts; k <= nlayers; ++k) o3vmr[K(k)] = o3mmr[K(k)];
    } else {
      for (int k = kts; k <= nlayers; ++k) o3vmr[K(k)] = o3mmr[K(k)] * amdo;
    }
    // cloud optical properties (LW 12032-12176); inflglw = 2 > 0
    for (int k = kts; k <= kte; ++k) cldfrac[K(k)] = cldfra1d[K(k)];
    const real gravmks = gravity;  // LW 12048
    {  // has_reqc = 0: relcalc below the model top (LW 12069-12074)
      const real landfrac = 2. - (*in.xland)[K(i)];
      const real landm = landfrac;
      const real snowh = 0.001 * (*in.snow)[K(i)];
      const real icefrac = (*in.xice)[K(i)];
      relcalc(kte - kts + 1, tlay.data(), landfrac, landm, icefrac, reliq.data(), snowh);
    }
    reicalc(kte - kts + 1, tlay.data(), reice.data());  // has_reqi = 0 (LW 12099-12105)
    if (iceflglw == 3) {
      for (int k = kts; k <= kte; ++k) {
        reice[K(k)] = reice[K(k)] * 1.0315;
        reice[K(k)] = std::min(140.0, reice[K(k)]);
      }
    }
    for (int k = kts; k <= kte; ++k) resnow[K(k)] = 10.;  // has_reqs = 0 (LW 12117-12121)
    for (int k = kts; k <= kte; ++k) {                     // LW 12124-12127
      const real gliqwp = qc1d[K(k)] * pdel[K(k)] * 100.0 / gravmks * 1000.0;
      cliqwp[K(k)] = gliqwp / std::max(0.01, cldfrac[K(k)]);
    }
    if (iceflglw == 3) {  // LW 12130-12135
      for (int k = kts; k <= kte; ++k) {
        const real gicewp = (qi1d[K(k)] + qs1d[K(k)]) * pdel[K(k)] * 100.0 / gravmks * 1000.0;
        cicewp[K(k)] = gicewp / std::max(0.01, cldfrac[K(k)]);
      }
    }
    for (int k = kts; k <= kte; ++k) {  // LW 12157-12164
      clwpth[K(k)] = cliqwp[K(k)];
      ciwpth[K(k)] = cicewp[K(k)];
      cswpth[K(k)] = csnowp[K(k)];
      rel[K(k)] = reliq[K(k)];
      rei[K(k)] = reice[K(k)];
      res[K(k)] = resnow[K(k)];
    }
    for (int k = kte + 1; k <= nlayers; ++k) {  // LW 12165-12175
      clwpth[K(k)] = 0.;
      ciwpth[K(k)] = 0.;
      cswpth[K(k)] = 0.;
      rel[K(k)] = 10.;
      rei[K(k)] = 10.;
      res[K(k)] = 10.;
      cldfrac[K(k)] = 0.;
    }

    // rrtmg_lw (LW 10705-11194) on CPPWRF's stages: inatm's inputs (LW 11197-11533)
    const rrtmg::LwLayout lay(nlayers);
    if (d.size() < lay.nd) d.assign(lay.nd, rb(0.0));
    if (iw.size() < lay.ni) iw.assign(lay.ni, 0);
    const rrtmg::LwViews v = lay.views(d.data(), iw.data(), 0, 1);
    v.pz[0] = plev[0];
    v.tz[0] = tlev[0];
    for (int l = 0; l < nlayers; ++l) {
      v.pavel[l] = play[l];
      v.tavel[l] = tlay[l];
      v.pz[l + 1] = plev[l + 1];
      v.tz[l + 1] = tlev[l + 1];
      v.wkl[0][l] = h2ovmr[l];
      v.wkl[1][l] = co2vmr[l];
      v.wkl[2][l] = o3vmr[l];
      v.wkl[3][l] = n2ovmr[l];
      v.wkl[4][l] = rb(0.0);
      v.wkl[5][l] = ch4vmr[l];
      v.wkl[6][l] = o2vmr[l];
      v.wx[0][l] = ccl4vmr[l];
      v.wx[1][l] = cfc11vmr[l];
      v.wx[2][l] = cfc12vmr[l];
      v.wx[3][l] = cfc22vmr[l];
      v.reicmc[l] = rei[l];
      v.relqmc[l] = rel[l];
      v.resnmc[l] = res[l];
      v.ciwpth[l] = ciwpth[l];
      v.clwpth[l] = clwpth[l];
      v.cswpth[l] = cswpth[l];
      v.cldfrac[l] = cldfrac[l];
      for (int g = 0; g < 140; ++g) {  // cldfmcl, taucmcl = 0 (LW 11899-11903)
        v.cldfmc[g][l] = rb(0.0);
        v.taucmc[g][l] = rb(0.0);
      }
    }
    for (int nb = 0; nb < 16; ++nb) v.semiss[nb] = emiss[K(i)];  // emis (LW 11774-11776), semiss (inatm)
    // McICA (LW 12181-12191): iplon = 1, irng = 0, permuteseed = 150
    std::uint64_t mask_lo = 0, mask_hi = 0;
    mcica_subcol(icld, 150, 140, nlayers, play.data(), cldfrac.data(), hlay.data(), in.idcor, in.julday, lat,
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
    int laytrop = 0;
    rrtmg::rrtmg_lw_setcoef(nlayers, 1, v.pavel, v.tavel, v.tz, tsfc, v.semiss, v.coldry, v.wkl, v.wbrodl, laytrop,
                            v.jp, v.jt, v.jt1, v.planklay, v.planklev, v.plankbnd, v.colh2o, v.colco2, v.colo3,
                            v.coln2o, v.colco, v.colch4, v.colo2, v.colbrd, v.fac00, v.fac01, v.fac10, v.fac11,
                            v.rat_h2oco2, v.rat_h2oco2_1, v.rat_h2oo3, v.rat_h2oo3_1, v.rat_h2on2o, v.rat_h2on2o_1,
                            v.rat_h2och4, v.rat_h2och4_1, v.rat_n2oco2, v.rat_n2oco2_1, v.rat_o3co2, v.rat_o3co2_1,
                            v.selffac, v.selffrac, v.indself, v.forfac, v.forfrac, v.indfor, v.minorfrac,
                            v.scaleminor, v.scaleminorn2, v.indminor, tables);
    v.iscal[0] = laytrop;
    for (int l = 0; l < nlayers; ++l) rrtmg::rrtmg_lw_layer(nlayers, l, v, tables);
    for (int b = 0; b < 16; ++b) rrtmg::rrtmg_lw_bandstage(nlayers, b, v, tables);
    for (int lev = 0; lev <= nlayers; ++lev) rrtmg::rrtmg_lw_level(lev, v, tables);
    rrtmg::rrtmg_lw_heating(nlayers, v.pz, v.fnet, v.fnetc, v.htr, v.htrc, heatfac);
    // uflx(1, k+1) = totuflux(k), hr(1, k+1) = htr(k) (LW 11182-11192)

    // outputs (LW 12212-12256)
    (*out.glw)[K(i)] = v.totdflux[0];
    (*out.olr)[K(i)] = v.totuflux[nlayers];
    (*out.lwcf)[K(i)] = v.totuclfl[nlayers] - v.totuflux[nlayers];
    (*out.lwupt)[K(i)] = v.totuflux[nlayers];
    (*out.lwuptc)[K(i)] = v.totuclfl[nlayers];
    (*out.lwdnt)[K(i)] = v.totdflux[nlayers];
    (*out.lwdntc)[K(i)] = v.totdclfl[nlayers];
    (*out.lwupb)[K(i)] = v.totuflux[0];
    (*out.lwupbc)[K(i)] = v.totuclfl[0];
    (*out.lwdnb)[K(i)] = v.totdflux[0];
    (*out.lwdnbc)[K(i)] = v.totdclfl[0];
    for (int k = kts; k <= kte; ++k) {  // LW 12241-12244
      const real tten1d = v.htr[k - 1] / 86400.;
      (*out.rthratenlw)(i, k) = tten1d / (*in.pi3d)(i, k);
    }
    if (in.icloud > 0) {  // LW 12247-12255
      for (int k = kts; k <= kte; ++k) {
        if ((*in.cldfra3d)(i, k) > 0.) {
          (*out.rre_cloud)(i, k) = rel[K(k)];
          (*out.rre_ice)(i, k) = rei[K(k)];
          (*out.rre_snow)(i, k) = res[K(k)];
        }
      }
    }
  }
}

// ---- rrtmg_swrad (SW 10025-10667) -------------------------------------------------------------------------------
void rrtmg_swrad(const RrtmgRadIn& in, const std::vector<real>& albedo, const std::vector<real>& xlong, real xtime,
                 real gmt, real radt, real degrad, real declin, real solcon, RrtmgSwradOut& out,
                 const rrtmg::RrtmgSwTables& tables, int its, int ite) {
  using physics_constants::cp;
  using physics_constants::gravity;
  if (its != 1) throw Error("rrtmg_swrad: its /= 1");
  const int kts = 1, kte = in.kte;
  if (in.has_reqc != 0 || in.has_reqi != 0 || in.has_reqs != 0)
    throw Error("rrtmg_swrad: the microphysics' effective radii (has_reqc / has_reqi / has_reqs = 1, "
                "config_microp_re = true) have no verified port");
  const real co2 = 379.e-6, ch4 = 1774.e-9, n2o = 319.e-9, o2 = 0.209488;  // SW 10137-10148
  const real amdw = 1.607793, amdo = 0.603461;                              // SW 10162-10164
  const int icld = in.cldovrlp;                                             // SW 10173
  int inflgsw = 2, iceflgsw = 3;                                            // SW 10176-10178
  const int liqflgsw = 1;
  const int nlay = (kte - kts + 1) + 1;  // SW 10258
  if (nlay > rrtmg::RRTMG_SW_COL_MAX) throw Error("rrtmg_swrad: more than 128 layers");
  // rrsw_con heatfac (rrtmg_sw_ini: swdatinit, SW 4856, 4876, 4912); rrsw_scon (SW 123)
  const rb heatfac = rb(9.8066) * rb(8.6400e4) / (cp * rb(1.e2));
  const rb rrsw_scon = rb(1.36822e+03);

  const int nlev = kte - kts + 1;
  std::vector<real> pw1d(nlev + 1), tw1d(nlev + 1), p1d(nlev), dz1d(nlev), t1d(nlev), qv1d(nlev), qc1d(nlev),
      qi1d(nlev), qs1d(nlev), cldfra1d(nlev), pdel(nlev), reliq(nlev), reice(nlev), resnow(nlev), cliqwp(nlev),
      cicewp(nlev), csnowp(nlev);
  std::vector<real> plev(nlay + 1), tlev(nlay + 1), play(nlay), hlay(nlay), tlay(nlay), h2ovmr(nlay), o3vmr(nlay),
      co2vmr(nlay), o2vmr(nlay), ch4vmr(nlay), n2ovmr(nlay), clwpth(nlay), ciwpth(nlay), cswpth(nlay), rel(nlay),
      rei(nlay), res(nlay), cldfrac(nlay), o3mmr(nlay);
  std::vector<rb> d;
  std::vector<int> iw;
  auto K = [](int k) { return static_cast<std::size_t>(k - 1); };

  for (int i = its; i <= ite; ++i) {
    // zenith angle (SW 10196-10206): sincos of xxlat and of declin, cos of hrang (the reference's calls)
    bool dorrsw = true;
    const real xt24 = std::fmod(xtime + radt * 0.5, 1440.);
    const real tloctm = gmt + xt24 / 60. + xlong[K(i)] / 15.;
    const real hrang = 15. * (tloctm - 12.) * degrad;
    const real xxlat = (*in.xlat)[K(i)] * degrad;
    real sinlat, coslat, sindec, cosdec;
    libm::sincos(xxlat, &sinlat, &coslat);
    libm::sincos(declin, &sindec, &cosdec);
    const real coszrs = sinlat * sindec + coslat * cosdec * libm::cos(hrang);
    (*out.coszr)[K(i)] = coszrs;
    if (coszrs <= 0.0) dorrsw = false;
    if (!dorrsw) {  // SW 10649-10661
      (*out.swupt)[K(i)] = 0.;
      (*out.swuptc)[K(i)] = 0.;
      (*out.swdnt)[K(i)] = 0.;
      (*out.swdntc)[K(i)] = 0.;
      (*out.swupb)[K(i)] = 0.;
      (*out.swupbc)[K(i)] = 0.;
      (*out.swdnb)[K(i)] = 0.;
      (*out.swdnbc)[K(i)] = 0.;
      continue;
    }
    const real lat = (*in.xlat)[K(i)];  // SW 10212
    for (int k = kts; k <= kte + 1; ++k) {  // SW 10216-10219
      pw1d[K(k)] = (*in.p8w)(i, k) / 100.;
      tw1d[K(k)] = (*in.t8w)(i, k);
    }
    for (int k = kts; k <= kte; ++k) {  // SW 10221-10234
      p1d[K(k)] = (*in.p3d)(i, k) / 100.;
      dz1d[K(k)] = (*in.dz8w)(i, k);
      t1d[K(k)] = (*in.t3d)(i, k);
      qv1d[K(k)] = std::max((*in.qv3d)(i, k), 1.e-12);
      qc1d[K(k)] = 0.;
      qi1d[K(k)] = 0.;
      qs1d[K(k)] = 0.;
      cldfra1d[K(k)] = 0.;
    }
    if (in.icloud > 0) {  // SW 10237-10246
      for (int k = kts; k <= kte; ++k) {
        qc1d[K(k)] = std::max((*in.qc3d)(i, k), 0.);
        qi1d[K(k)] = std::max((*in.qi3d)(i, k), 0.);
        qs1d[K(k)] = std::max((*in.qs3d)(i, k), 0.);
        cldfra1d[K(k)] = (*in.cldfra3d)(i, k);
      }
    }
    for (int k = kts; k <= kte; ++k) {  // SW 10262-10270
      reliq[K(k)] = 10.;
      reice[K(k)] = 10.;
      resnow[K(k)] = 10.;
      cliqwp[K(k)] = 0.;
      cicewp[K(k)] = 0.;
      csnowp[K(k)] = 0.;
    }
    for (int k = 1; k <= nlay; ++k) {
      clwpth[K(k)] = ciwpth[K(k)] = cswpth[K(k)] = 0.;
      rel[K(k)] = rei[K(k)] = res[K(k)] = cldfrac[K(k)] = 0.;
    }
    // input sounding (SW 10325-10354)
    plev[K(1)] = pw1d[K(1)];
    tlev[K(1)] = tw1d[K(1)];
    for (int k = kts; k <= kte; ++k) {
      plev[K(k + 1)] = pw1d[K(k + 1)];
      tlev[K(k + 1)] = tw1d[K(k + 1)];
      play[K(k)] = p1d[K(k)];
      pdel[K(k)] = plev[K(k)] - plev[K(k + 1)];
      tlay[K(k)] = t1d[K(k)];
      h2ovmr[K(k)] = qv1d[K(k)] * amdw;
      co2vmr[K(k)] = co2;
      o2vmr[K(k)] = o2;
      ch4vmr[K(k)] = ch4;
      n2ovmr[K(k)] = n2o;
    }
    plev[K(kte + 2)] = 1.0e-5;
    tlev[K(kte + 2)] = tlev[K(kte + 1)] + 0.0;
    play[K(kte + 1)] = 0.5 * plev[K(kte + 1)];
    tlay[K(kte + 1)] = tlev[K(kte + 1)] + 0.0;
    h2ovmr[K(kte + 1)] = h2ovmr[K(kte)];
    co2vmr[K(kte + 1)] = co2vmr[K(kte)];
    o2vmr[K(kte + 1)] = o2vmr[K(kte)];
    ch4vmr[K(kte + 1)] = ch4vmr[K(kte)];
    n2ovmr[K(kte + 1)] = n2ovmr[K(kte)];
    real dzsum = 0.;  // SW 10359-10364
    for (int k = kts; k <= kte; ++k) {
      hlay[K(k)] = dzsum + 0.5 * dz1d[K(k)];
      dzsum = dzsum + dz1d[K(k)];
    }
    hlay[K(kte + 1)] = dzsum + 0.5 * dz1d[K(kte)];
    // ozone (SW 10367-10381)
    inirad(o3mmr.data(), plev.data(), kts, kte);
    if (in.o3input == 2) {
      std::vector<real> o3clim1d(static_cast<std::size_t>(in.noznlevels));
      for (int k = 1; k <= in.noznlevels; ++k) o3clim1d[K(k)] = (*in.o3clim)(i, k);
      vinterp_ozn(1, 1, 1, kte + 1, play.data(), in.pin, in.noznlevels, o3clim1d.data(), o3mmr.data());
      for (int k = kts; k <= kte + 1; ++k) o3vmr[K(k)] = o3mmr[K(k)];
    } else {
      for (int k = kts; k <= kte + 1; ++k) o3vmr[K(k)] = o3mmr[K(k)] * amdo;
    }
    // surface albedo, solar constant, zenith angle (SW 10384-10397): dyofyr = 0, adjes = 1
    const real asdir = albedo[K(i)], asdif = albedo[K(i)], aldir = albedo[K(i)], aldif = albedo[K(i)];
    const real scon = solcon;
    const real coszen = coszrs;
    const real adjes = 1.0;
    // cloud optical properties (SW 10400-10535)
    for (int k = kts; k <= kte; ++k) cldfrac[K(k)] = cldfra1d[K(k)];
    const real gravmks = gravity;
    {  // has_reqc = 0 (SW 10434-10443)
      const real landfrac = 2. - (*in.xland)[K(i)];
      const real landm = landfrac;
      const real snowh = 0.001 * (*in.snow)[K(i)];
      const real icefrac = (*in.xice)[K(i)];
      relcalc(kte - kts + 1, tlay.data(), landfrac, landm, icefrac, reliq.data(), snowh);
    }
    reicalc(kte - kts + 1, tlay.data(), reice.data());  // has_reqi = 0 (SW 10463-10471)
    if (iceflgsw == 3) {
      for (int k = kts; k <= kte; ++k) {
        reice[K(k)] = reice[K(k)] * 1.0315;
        reice[K(k)] = std::min(140.0, reice[K(k)]);
      }
    }
    for (int k = kts; k <= kte; ++k) resnow[K(k)] = 10.;  // has_reqs = 0 (SW 10481-10485)
    for (int k = kts; k <= kte; ++k) {                     // SW 10488-10491
      const real gliqwp = qc1d[K(k)] * pdel[K(k)] * 100.0 / gravmks * 1000.0;
      cliqwp[K(k)] = gliqwp / std::max(0.01, cldfrac[K(k)]);
    }
    if (iceflgsw == 3) {  // SW 10494-10499
      for (int k = kts; k <= kte; ++k) {
        const real gicewp = (qi1d[K(k)] + qs1d[K(k)]) * pdel[K(k)] * 100.0 / gravmks * 1000.0;
        cicewp[K(k)] = gicewp / std::max(0.01, cldfrac[K(k)]);
      }
    }
    for (int k = kts; k <= kte; ++k) {  // SW 10521-10528
      clwpth[K(k)] = cliqwp[K(k)];
      ciwpth[K(k)] = cicewp[K(k)];
      cswpth[K(k)] = csnowp[K(k)];
      rel[K(k)] = reliq[K(k)];
      rei[K(k)] = reice[K(k)];
      res[K(k)] = resnow[K(k)];
    }
    clwpth[K(kte + 1)] = 0.;  // SW 10531-10537
    ciwpth[K(kte + 1)] = 0.;
    cswpth[K(kte + 1)] = 0.;
    rel[K(kte + 1)] = 10.;
    rei[K(kte + 1)] = 10.;
    res[K(kte + 1)] = 10.;
    cldfrac[K(kte + 1)] = 0.;

    // rrtmg_sw (SW 8857-9665) on CPPWRF's stages
    const rrtmg::SwLayout lay(nlay);
    if (d.size() < lay.nd) d.assign(lay.nd, rb(0.0));
    if (iw.size() < lay.ni) iw.assign(lay.ni, 0);
    const rrtmg::SwViews v = lay.views(d.data(), iw.data(), 0, 1);
    v.pz[0] = plev[0];
    v.tz[0] = tlev[0];
    for (int l = 0; l < nlay; ++l) {
      v.pavel[l] = play[l];
      v.tavel[l] = tlay[l];
      v.pz[l + 1] = plev[l + 1];
      v.tz[l + 1] = tlev[l + 1];
      v.wkl[0][l] = h2ovmr[l];
      v.wkl[1][l] = co2vmr[l];
      v.wkl[2][l] = o3vmr[l];
      v.wkl[3][l] = n2ovmr[l];
      v.wkl[4][l] = rb(0.0);
      v.wkl[5][l] = ch4vmr[l];
      v.wkl[6][l] = o2vmr[l];
      v.cldfrac[l] = cldfrac[l];
      v.clwpth[l] = clwpth[l];
      v.ciwpth[l] = ciwpth[l];
      v.cswpth[l] = cswpth[l];
      v.relqmc[l] = rel[l];
      v.reicmc[l] = rei[l];
      v.resnmc[l] = res[l];
      for (int g = 0; g < 112; ++g) v.cldfmc[l][g] = rb(0.0);  // cldfmcl = 0 (SW 10292)
    }
    v.iscal[2] = inflgsw;
    v.iscal[3] = iceflgsw;
    v.iscal[4] = liqflgsw;
    // McICA (SW 10540-10549): iplon = 1, irng = 0, permuteseed = 1
    int any_cld_subcol = 0;
    mcica_subcol(icld, 1, 112, nlay, play.data(), cldfrac.data(), hlay.data(), in.idcor, in.julday, lat,
                 [&](int isub, int ilev, bool cld) {
                   v.cldfmc[ilev][isub] = cld ? rb(1.0) : rb(0.0);
                   if (cld) any_cld_subcol = 1;
                 });
    v.iscal[5] = any_cld_subcol;
    rrtmg::rrtmg_sw_core_prep(nlay, v, tables, coszen, asdir, asdif, aldir, aldif);
    rrtmg::RrtmgSwCallConfig c{};
    c.num_levels = nlay - 1;
    c.nlayers = nlay;
    c.heatfac = heatfac;
    c.adjflux = adjes * (scon / rrsw_scon);  // inatm_sw (SW 9909-9916)
    for (int g = 0; g < 112; ++g) rrtmg::rrtmg_sw_gpoint(c, g, v, tables);
    for (int lev = 0; lev <= nlay; ++lev) rrtmg::rrtmg_sw_level(lev, v, tables);
    // swhr (rrtmg_sw, SW 9558-9566)
    for (int l = 0; l < nlay; ++l) {
      const rb zdpgcp = heatfac / v.pdp[l];
      v.swhrc[l] = ((v.cd[l + 1] - v.cu[l + 1]) - (v.cd[l] - v.cu[l])) * zdpgcp;
      v.swhr[l] = ((v.fd[l + 1] - v.fu[l + 1]) - (v.fd[l] - v.fu[l])) * zdpgcp;
    }
    v.swhrc[nlay - 1] = rb(0.);
    v.swhr[nlay - 1] = rb(0.);

    // outputs (SW 10598-10647): swuflx(1, k+1) = fu(k), swdkdir(1, 1) = fddir(0), swdkdif(1, 1) = fd(0) - fddir(0)
    (*out.gsw)[K(i)] = v.fd[0] - v.fu[0];
    (*out.swcf)[K(i)] = (v.fd[kte + 1] - v.fu[kte + 1]) - (v.cd[kte + 1] - v.cu[kte + 1]);
    (*out.swupt)[K(i)] = v.fu[kte + 1];
    (*out.swuptc)[K(i)] = v.cu[kte + 1];
    (*out.swdnt)[K(i)] = v.fd[kte + 1];
    (*out.swdntc)[K(i)] = v.cd[kte + 1];
    (*out.swupb)[K(i)] = v.fu[0];
    (*out.swupbc)[K(i)] = v.cu[0];
    (*out.swdnb)[K(i)] = v.fd[0];
    (*out.swdnbc)[K(i)] = v.cd[0];
    (*out.swddir)[K(i)] = v.fddir[0];
    (*out.swddni)[K(i)] = (*out.swddir)[K(i)] / coszrs;
    (*out.swddif)[K(i)] = v.fd[0] - v.fddir[0];
    for (int k = kts; k <= kte; ++k) {  // SW 10640-10644
      const real tten1d = v.swhr[k - 1] / 86400.;
      (*out.rthratensw)(i, k) = tten1d / (*in.pi3d)(i, k);
    }
  }
}

}  // namespace cppmpas::phys
