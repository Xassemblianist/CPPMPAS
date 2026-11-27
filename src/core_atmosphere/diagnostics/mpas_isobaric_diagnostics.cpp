// SPDX-License-Identifier: GPL-3.0-or-later
// mpas_isobaric_diagnostics.F (IS), statement by statement: isobaric_diagnostics_compute (IS:62-306),
// interp_diagnostics (IS:310-1023), interp_tofixed_pressure (IS:1027-1110), compute_slp (IS:1113-1243) and
// compute_layer_mean (IS:1262-1351). One rank, one block: nCells = nCellsSolve, the halo exchange of pressure_p
// (IS:452-453) is a no-op.
//
// The libm calls are the reference executable's (objdump -d of atmosphere_model,
// __mpas_isobaric_diagnostics_MOD_interp_diagnostics.constprop.0 and _compute_slp):
//  * the three pressure2 loops (IS:631-664) are vectorised over iCell, two cells per iteration from iCell = 1:
//    _ZGVbN2v_log, _ZGVbN2v_log, _ZGVbN2v_exp per pair (0x241a9c..0x241acb, 0x241e93..0x241ec2,
//    0x2422f6..0x242324), the last cell of an odd nCells with log / log / exp (0x241bf8..0x241c43, ...); with
//    nCells = 1 the scalar loop alone (0x2468d0);
//  * the temperature / dewpoint loop (IS:680-693) is vectorised over k, two levels per iteration from k = 1:
//    one _ZGVbN2v_log per pair (0x2448d0; log(evp/6.112) is evaluated once), the last level of an odd nVertLevels
//    with log (0x244a44);
//  * compute_slp: the t_at_pconst / z_at_pconst loop (IS:1197-1220) calls log twice and pow (scalar,
//    0x23ece8..0x23edc4); the slp loop (IS:1238-1242) is vectorised over icol, two columns per iteration from
//    icol = 1, with _ZGVbN2v_exp (0x23f16b) and exp for the last column of an odd ncol (0x23f22f).
#include <algorithm>
#include <cstddef>
#include <string>

#include "core/backend.hpp"
#include "core/error.hpp"
#include "core/field.hpp"
#include "core/libm.hpp"
#include "core_atmosphere/diagnostics/mpas_atm_diagnostics.hpp"
#include "physics/mpas_atmphys_constants.hpp"

namespace cppmpas::diag {

namespace {

using physics_constants::rvord;

// IS:1027-1110 interp_tofixed_pressure; arrays (ncol, nlev) as MPAS's
// The column loops are joined by kkstart and kount, so the routine runs as one sequence (one GPU thread, serial).
void interp_tofixed_pressure(int ncol, int nlev_in, int nlev_out, const FArray<real>& pres_in_a,
                             const FArray<real>& field_in_a, const FArray<real>& pres_out_a, FArray<real>& field_out_a) {
  FStorage<int> kupper_a(static_cast<std::size_t>(ncol) + 1);
  int* kupper = kupper_a.data();
  const auto pres_in = pres_in_a.cview(), field_in = field_in_a.cview(), pres_out = pres_out_a.cview();
  const auto field_out = field_out_a.view();
  serial([=] CPPMPAS_LAMBDA() {
  for (int icol = 1; icol <= ncol; ++icol) kupper[icol] = 1;
  for (int k = 1; k <= nlev_out; ++k) {
    int kkstart = nlev_in;
    for (int icol = 1; icol <= ncol; ++icol) kkstart = kkstart < kupper[icol] ? kkstart : kupper[icol];  // min
    int kount = 0;
    bool done = false;
    for (int kk = kkstart; kk <= nlev_in - 1 && !done; ++kk) {
      for (int icol = 1; icol <= ncol; ++icol) {
        if (pres_out(icol, k) > pres_in(icol, kk) && pres_out(icol, k) <= pres_in(icol, kk + 1)) {
          kupper[icol] = kk;
          kount = kount + 1;
        }
      }
      if (kount == ncol) {
        for (int icol = 1; icol <= ncol; ++icol) {
          const real dpu = pres_out(icol, k) - pres_in(icol, kupper[icol]);
          const real dpl = pres_in(icol, kupper[icol] + 1) - pres_out(icol, k);
          field_out(icol, k) =
              (field_in(icol, kupper[icol]) * dpl + field_in(icol, kupper[icol] + 1) * dpu) / (dpl + dpu);
        }
        done = true;  // goto 35
      }
    }
    if (done) continue;
    for (int icol = 1; icol <= ncol; ++icol) {
      if (pres_out(icol, k) < pres_in(icol, 1)) {
        field_out(icol, k) = field_in(icol, 1) * pres_out(icol, k) / pres_in(icol, 1);
      } else if (pres_out(icol, k) > pres_in(icol, nlev_in)) {
        field_out(icol, k) = field_in(icol, nlev_in);
      } else {
        const real dpu = pres_out(icol, k) - pres_in(icol, kupper[icol]);
        const real dpl = pres_in(icol, kupper[icol] + 1) - pres_out(icol, k);
        field_out(icol, k) =
            (field_in(icol, kupper[icol]) * dpl + field_in(icol, kupper[icol] + 1) * dpu) / (dpl + dpu);
      }
    }
  }
  });
}

// IS:1113-1243 compute_slp; p in hPa, t in K, height in m; slp(1:ncol)
// A column without a trapping level ends the routine for all columns, so it runs as one sequence (serial).
void compute_slp(int ncol, int nlev_in, const FArray<real>& t_a, const FArray<real>& height_a, const FArray<real>& p_a,
                 int index_qv, const FArray<real>& scalars_a, FArray<real>& slp_a) {
  constexpr real rr = 287.0, grav = 9.80616, gamma = 0.0065;  // IS:1134
  constexpr real tc = 273.16 + 17.5, pconst = 100.;           // IS:1137
  constexpr bool mm5_test = true;                             // IS:1140
  constexpr real gamma_rr_grav = gamma * rr / grav;           // folded by gfortran
  FStorage<int> level_a(static_cast<std::size_t>(ncol) + 1);
  FStorage<real> t_surf_a(static_cast<std::size_t>(ncol) + 1), t_msl_a(static_cast<std::size_t>(ncol) + 1);
  int* level = level_a.data();
  real *t_surf = t_surf_a.data(), *t_msl = t_msl_a.data();
  const auto t = t_a.cview(), height = height_a.cview(), p = p_a.cview(), scalars = scalars_a.cview();
  const auto slp = slp_a.view();
  serial([=] CPPMPAS_LAMBDA() {
  auto zero_slp = [&]() {  // slp(:) = 0.0 (the dummy's extent: ncol)
    for (int icol = 1; icol <= ncol; ++icol) slp(icol) = 0.0;
  };
  for (int icol = 1; icol <= ncol; ++icol) {  // IS:1159-1181
    level[icol] = -1;
    int k = 1;
    bool found = false;
    while (!found && k <= nlev_in) {
      if (p(k, icol) < p(1, icol) - pconst) {
        level[icol] = k;
        found = true;
      }
      k = k + 1;
    }
    if (level[icol] == -1) {  // 'Troubles finding level ... MSLP field will not be computed' (log only)
      zero_slp();
      return;
    }
  }
  for (int icol = 1; icol <= ncol; ++icol) {  // IS:1186-1221
    const int klo = level[icol] - 1 > 1 ? level[icol] - 1 : 1;  // max(level - 1, 1)
    const int khi = klo + 1 < nlev_in - 1 ? klo + 1 : nlev_in - 1;  // min(klo + 1, nlev_in - 1)
    if (klo == khi) {  // 'Trapping levels are weird.' (log only)
      zero_slp();
      return;
    }
    const real plo = p(klo, icol);
    const real phi = p(khi, icol);
    const real tlo = t(klo, icol) * (1. + 0.608 * scalars(index_qv, klo, icol));
    const real thi = t(khi, icol) * (1. + 0.608 * scalars(index_qv, khi, icol));
    const real zlo = 0.5 * (height(klo, icol) + height(klo + 1, icol));
    const real zhi = 0.5 * (height(khi, icol) + height(khi + 1, icol));
    const real p_at_pconst = p(1, icol) - pconst;
    const real l1 = libm::log(p_at_pconst / phi), l2 = libm::log(plo / phi);  // log@plt (0x23ece8, 0x23ed03)
    const real t_at_pconst = thi - (thi - tlo) * l1 * l2;
    const real z_at_pconst = zhi - (zhi - zlo) * l1 * l2;
    t_surf[icol] = t_at_pconst * libm::pow(p(1, icol) / p_at_pconst, gamma_rr_grav);  // pow@plt (0x23edc4)
    t_msl[icol] = t_at_pconst + gamma * z_at_pconst;
  }
  if (mm5_test) {  // IS:1226-1236
    for (int icol = 1; icol <= ncol; ++icol) {
      const bool l1 = t_msl[icol] < tc;
      const bool l2 = t_surf[icol] <= tc;
      const bool l3 = !l1;
      if (l2 && l3) {
        t_msl[icol] = tc;
      } else {
        const real d = t_surf[icol] - tc;
        t_msl[icol] = tc - 0.005 * (d * d);
      }
    }
  }
  for (int icol = 1; icol <= ncol; ++icol) {  // IS:1238-1242: _ZGVbN2v_exp in pairs from icol = 1, exp for an odd tail
    const real z_half_lowest = 0.5 * (height(1, icol) + height(2, icol));
    const real x = (2. * grav * z_half_lowest) / (rr * (t_msl[icol] + t_surf[icol]));
    const bool lane = icol <= ncol - ncol % 2;
    slp(icol) = p(1, icol) * (lane ? libm::mvec_exp(x) : libm::exp(x));
  }
  });
}

// IS:1262-1351 compute_layer_mean; field_in / press_in (nCells, nVertLevels), pressure increasing with index
// wtop_* / wbot_* live across the cells in MPAS: the cell loop runs as one sequence (serial).
void compute_layer_mean(FArray<real>& layerMean_a, real p1, real p2, const FArray<real>& field_in_a,
                        const FArray<real>& press_in_a) {
  const int nCells = field_in_a.extent(0);
  const int nVertLevels = field_in_a.extent(1);
  if (press_in_a(1, 1) > press_in_a(1, nVertLevels)) {  // 'Error in compute_layer_mean' (log only)
    for (std::size_t n = 0; n < layerMean_a.size(); ++n) layerMean_a.data()[n] = 0.0;
    return;
  }
  const auto field_in = field_in_a.cview(), press_in = press_in_a.cview();
  const auto layerMean = layerMean_a.view();
  real p_top, p_bot;
  if (p1 < p2) {
    p_top = p1;
    p_bot = p2;
  } else {
    p_top = p2;
    p_bot = p1;
  }
  serial([=] CPPMPAS_LAMBDA() {
  real wtop_p = 0, wtop_m = 0, wbot_p = 0, wbot_m = 0;  // unset in MPAS until a trapping level is found
  for (int iCell = 1; iCell <= nCells; ++iCell) {
    int k_bot = -1, k_top = -1;
    for (int k = 1; k <= nVertLevels - 1; ++k) {
      if (press_in(iCell, k) <= p_top && press_in(iCell, k + 1) > p_top) {
        k_top = k;
        wtop_p = (p_top - press_in(iCell, k)) / (press_in(iCell, k + 1) - press_in(iCell, k));
        wtop_m = (press_in(iCell, k + 1) - p_top) / (press_in(iCell, k + 1) - press_in(iCell, k));
      }
      if (press_in(iCell, k) <= p_bot && press_in(iCell, k + 1) > p_bot) {
        k_bot = k;
        wbot_m = (p_bot - press_in(iCell, k)) / (press_in(iCell, k + 1) - press_in(iCell, k));
        wbot_p = (press_in(iCell, k + 1) - p_bot) / (press_in(iCell, k + 1) - press_in(iCell, k));
      }
    }
    if (k_top == -1 || k_bot == -1) {
      layerMean(iCell) = 0.0;
    } else if (k_top == k_bot) {
      layerMean(iCell) = wtop_m * field_in(iCell, k_top) + wtop_p * field_in(iCell, k_top + 1);
      layerMean(iCell) = layerMean(iCell) + wbot_m * field_in(iCell, k_bot) + wbot_p * field_in(iCell, k_bot + 1);
      layerMean(iCell) = 0.5 * layerMean(iCell);
    } else {
      real wtotal = press_in(iCell, k_top + 1) - p_top;
      real temp = wtop_m * field_in(iCell, k_top) + wtop_p * field_in(iCell, k_top + 1);
      layerMean(iCell) = wtotal * 0.5 * (field_in(iCell, k_top + 1) + temp);
      for (int k = k_top + 1; k <= k_bot - 1; ++k) {
        const real w = press_in(iCell, k + 1) - press_in(iCell, k);
        wtotal = wtotal + w;
        layerMean(iCell) = layerMean(iCell) + w * 0.5 * (field_in(iCell, k) + field_in(iCell, k + 1));
      }
      const real w = p_bot - press_in(iCell, k_bot);
      wtotal = wtotal + w;
      temp = wbot_m * field_in(iCell, k_bot) + wbot_p * field_in(iCell, k_bot + 1);
      layerMean(iCell) = layerMean(iCell) + w * 0.5 * (field_in(iCell, k_bot) + temp);
      layerMean(iCell) = layerMean(iCell) / wtotal;
    }
  }
  });
}

const char* const kLevels[8] = {"50", "100", "200", "250", "500", "700", "850", "925"};

struct Needs {
  bool mslp = false, temp = false, relhum = false, dewpoint = false, w = false, uzonal = false, umeridional = false,
       vorticity = false, height = false, t_isobaric = false, z_isobaric = false, meanT_500_300 = false;
};

// The copies into and out of interp_tofixed_pressure's (ncol, nlev) arrays (IS:719-1046): one thread per column.
void flip_in(FArray<real>& dst_a, const FArray<real>& src_a, int ncol, int nlev, real scale) {
  const auto dst = dst_a.view();
  const auto src = src_a.cview();
  if (scale == 1.0)
    parallel_for(1, ncol, [=] CPPMPAS_LAMBDA(int iCell) {
      for (int k = 1; k <= nlev; ++k) dst(iCell, nlev + 1 - k) = src(k, iCell);
    });
  else
    parallel_for(1, ncol, [=] CPPMPAS_LAMBDA(int iCell) {
      for (int k = 1; k <= nlev; ++k) dst(iCell, nlev + 1 - k) = src(k, iCell) * scale;
    });
}
void copy_level(FArray<real>& out_a, const FArray<real>& fi_a, int ncol, int l) {
  const auto out = out_a.view();
  const auto fi = fi_a.cview();
  parallel_for(1, ncol, [=] CPPMPAS_LAMBDA(int iCell) { out(iCell) = fi(iCell, l); });
}
void fill_levels(FArray<real>& pi_a, int ncol, int nl, const real* lv) {
  const auto pi = pi_a.view();
  FStorage<real> lv_a(lv, lv + nl);
  const real* lvp = lv_a.data();
  parallel_for(1, ncol, [=] CPPMPAS_LAMBDA(int iCell) {
    for (int l = 1; l <= nl; ++l) pi(iCell, l) = lvp[l - 1];
  });
}

// IS:310-1023 interp_diagnostics (time_lev = 1). Every loop over cells or vertices is one kernel (one thread per cell
// or vertex, the level loops inside the thread in MPAS's order).
void interp_diagnostics(DiagnosticsContext& c, const Needs& need) {
  const Mesh& m = *c.mesh;
  State& s = *c.state;
  Diag& d = *c.diag;
  RegistryPools& p = *c.pools;
  const int nCells = m.nCells, nVertLevels = m.nVertLevels, nVertices = m.nVertices;
  const int vertexDegree = m.vertexDegree;
  const int nIsoLevelsT = p.dim("nIsoLevelsT", false), nIsoLevelsZ = p.dim("nIsoLevelsZ", false);
  const int nVertLevelsP1 = nVertLevels + 1;
  int index_qv = 0;
  for (std::size_t k = 0; k < s.scalar_names.size(); ++k)
    if (s.scalar_names[k] == "qv") index_qv = static_cast<int>(k) + 1;
  if (index_qv == 0) throw Error("isobaric diagnostics: qv is not an active scalar");
  const auto cellsOnVertex = m.cellsOnVertex.cview();
  const auto areaTriangle = m.areaTriangle.cview();
  const auto kiteAreasOnVertex = m.kiteAreasOnVertex.cview();
  const auto height = m.zgrid.cview();
  const auto theta_m = s.theta_m.cview();
  const auto scalars = s.scalars.cview();
  const auto exner = d.exner.cview();
  const auto pressure_b = d.pressure_base.cview();
  const auto pressure_p = d.pressure_p.cview();
  auto& relhum = p.R("relhum");
  auto& t_iso_levels_a = p.R("t_iso_levels");
  auto& z_iso_levels_a = p.R("z_iso_levels");
  const auto t_iso_levels = t_iso_levels_a.view(), z_iso_levels = z_iso_levels_a.view();
  const auto t_isobaric = p.R("t_isobaric").view();
  const auto z_isobaric = p.R("z_isobaric").view();
  auto& meanT_500_300 = p.R("meanT_500_300");
  auto lev = [&](const char* f, int l) -> FArray<real>& { return p.R(std::string(f) + "_" + kLevels[l] + "hPa"); };
  auto& mslp_a = p.R("mslp");

  // IS:609-614 (allocated, not initialised: temperature and dewpoint are only read where they were computed)
  FArray<real> pressure_a({nVertLevels, nCells}), pressureCp1_a({nVertLevels, nCells + 1});
  FArray<real> pressure2_a({nVertLevelsP1, nCells}), pressure_v_a({nVertLevels, nVertices});
  FArray<real> temperature_a({nVertLevels, nCells}), dewpoint_a({nVertLevels, nCells});
  const auto pressure = pressure_a.view(), pressureCp1 = pressureCp1_a.view(), pressure2 = pressure2_a.view(),
             pressure_v = pressure_v_a.view(), temperature = temperature_a.view(), dewpoint = dewpoint_a.view();

  if (need.t_isobaric) {  // IS:616-622
    t_iso_levels_a(1) = 30000.0;
    t_iso_levels_a(2) = 35000.0;
    t_iso_levels_a(3) = 40000.0;
    t_iso_levels_a(4) = 45000.0;
    t_iso_levels_a(5) = 50000.0;
  }
  if (need.z_isobaric) {  // IS:624-638
    for (int k = 1; k <= 13; ++k) z_iso_levels_a(k) = 30000.0 + 5000.0 * (k - 1);  // 30000.0 ... 90000.0, exact
  }

  // IS:641-651 total pressure at cell centres (mass points), in hPa
  parallel_for(1, nCells, [=] CPPMPAS_LAMBDA(int iCell) {
    for (int k = 1; k <= nVertLevels; ++k) {
      pressure(k, iCell) = (pressure_p(k, iCell) + pressure_b(k, iCell)) / 100.;
      pressureCp1(k, iCell) = pressure(k, iCell);
    }
  });
  parallel_for(nCells + 1, nCells + 1, [=] CPPMPAS_LAMBDA(int iCell) {
    for (int k = 1; k <= nVertLevels; ++k)
      pressureCp1(k, iCell) = (pressure_p(k, iCell) + pressure_b(k, iCell)) / 100.;
  });

  // IS:654-688 total pressure at w points; vectorised over iCell (pairs from iCell = 1, odd tail scalar)
  const int nPair = nCells - nCells % 2;
#define logexp(iCell, w1, a, w2, b)                                                                     \
  ((iCell) <= nPair ? libm::mvec_exp((w1) * libm::mvec_log(a) + (w2) * libm::mvec_log(b))              \
                    : libm::exp((w1) * libm::log(a) + (w2) * libm::log(b)))
  parallel_for(1, nCells, [=] CPPMPAS_LAMBDA(int iCell) {
    const int k = nVertLevelsP1;
    const real z0 = height(k, iCell);
    const real z1 = 0.5 * (height(k, iCell) + height(k - 1, iCell));
    const real z2 = 0.5 * (height(k - 1, iCell) + height(k - 2, iCell));
    const real w1 = (z0 - z2) / (z1 - z2);
    const real w2 = 1. - w1;
    pressure2(k, iCell) = logexp(iCell, w1, pressure(k - 1, iCell), w2, pressure(k - 2, iCell));
  });
  parallel_for(1, nCells, [=] CPPMPAS_LAMBDA(int iCell) {
    for (int k = 2; k <= nVertLevels; ++k) {
      const real w1 = (height(k, iCell) - height(k - 1, iCell)) / (height(k + 1, iCell) - height(k - 1, iCell));
      const real w2 = (height(k + 1, iCell) - height(k, iCell)) / (height(k + 1, iCell) - height(k - 1, iCell));
      pressure2(k, iCell) = logexp(iCell, w1, pressure(k, iCell), w2, pressure(k - 1, iCell));
    }
  });
  parallel_for(1, nCells, [=] CPPMPAS_LAMBDA(int iCell) {
    const int k = 1;
    const real z0 = height(k, iCell);
    const real z1 = 0.5 * (height(k, iCell) + height(k + 1, iCell));
    const real z2 = 0.5 * (height(k + 1, iCell) + height(k + 2, iCell));
    const real w1 = (z0 - z2) / (z1 - z2);
    const real w2 = 1. - w1;
    pressure2(k, iCell) = logexp(iCell, w1, pressure(k, iCell), w2, pressure(k + 1, iCell));
  });
#undef logexp

  // IS:691-700 total pressure at vertices (mass points)
  parallel_for(1, nVertices, [=] CPPMPAS_LAMBDA(int iVert) {
    for (int k = 1; k <= nVertLevels; ++k) pressure_v(k, iVert) = 0.;
    for (int k = 1; k <= nVertLevels; ++k) {
      for (int iVertD = 1; iVertD <= vertexDegree; ++iVertD)
        pressure_v(k, iVert) = pressure_v(k, iVert) + kiteAreasOnVertex(iVertD, iVert) *
                                                          pressureCp1(k, cellsOnVertex(iVertD, iVert));
      pressure_v(k, iVert) = pressure_v(k, iVert) / areaTriangle(iVert);
    }
  });

  // IS:702-716 temperature and dewpoint; vectorised over k (pairs from k = 1, an odd last level scalar)
  if (need.temp || need.relhum || need.dewpoint || need.mslp) {
    const int kPair = nVertLevels - nVertLevels % 2;
    parallel_for(1, nCells, [=] CPPMPAS_LAMBDA(int iCell) {
      for (int k = 1; k <= nVertLevels; ++k) {
        temperature(k, iCell) =
            (theta_m(k, iCell) / (1. + rvord * scalars(index_qv, k, iCell))) * exner(k, iCell);
        real evp = pressure(k, iCell) * scalars(index_qv, k, iCell) / (scalars(index_qv, k, iCell) + 0.622);
        evp = evp < 1.0e-8 ? 1.0e-8 : evp;  // std::max(evp, 1.0e-8)
        const real le = k <= kPair ? libm::mvec_log(evp / 6.112) : libm::log(evp / 6.112);
        dewpoint(k, iCell) = (243.5 * le) / (17.67 - le);
        dewpoint(k, iCell) = dewpoint(k, iCell) + 273.15;
      }
    });
  }

  // IS:719-743 interpolation to fixed pressure levels of fields at cell centres and mass points
  int nIntP = 8;
  FArray<real> field_interp({nCells, nIntP}), press_interp({nCells, nIntP});
  static const real kPress[8] = {50.0, 100.0, 200.0, 250.0, 500.0, 700.0, 850.0, 925.0};
  fill_levels(press_interp, nCells, 8, kPress);
  {
    FArray<real> press_in({nCells, nVertLevels}), field_in({nCells, nVertLevels});
    flip_in(press_in, pressure_a, nCells, nVertLevels, 1.0);
    auto mass_level = [&](bool needed, const FArray<real>& src, const char* name) {
      if (!needed) return;
      flip_in(field_in, src, nCells, nVertLevels, 1.0);
      interp_tofixed_pressure(nCells, nVertLevels, nIntP, press_in, field_in, press_interp, field_interp);
      for (int l = 0; l < 8; ++l) copy_level(lev(name, l), field_interp, nCells, l + 1);
    };
    mass_level(need.temp, temperature_a, "temperature");     // IS:745-764
    mass_level(need.relhum, relhum, "relhum");               // IS:767-785
    mass_level(need.dewpoint, dewpoint_a, "dewpoint");       // IS:787-805
    mass_level(need.uzonal, d.uReconstructZonal, "uzonal");  // IS:807-825
    mass_level(need.umeridional, d.uReconstructMeridional, "umeridional");  // IS:827-845
  }

  if (need.w || need.height) {  // IS:850-902: fields at cell centres and w points
    FArray<real> press_in({nCells, nVertLevelsP1}), field_in({nCells, nVertLevelsP1});
    flip_in(press_in, pressure2_a, nCells, nVertLevelsP1, 1.0);
    auto w_level = [&](const FArray<real>& src, const char* name) {
      flip_in(field_in, src, nCells, nVertLevelsP1, 1.0);
      interp_tofixed_pressure(nCells, nVertLevelsP1, nIntP, press_in, field_in, press_interp, field_interp);
      for (int l = 0; l < 8; ++l) copy_level(lev(name, l), field_interp, nCells, l + 1);
    };
    w_level(m.zgrid, "height");
    w_level(s.w, "w");
  }

  if (need.vorticity) {  // IS:907-956: fields at vertices and mass points
    nIntP = 8;
    FArray<real> fi({nVertices, nIntP}), pi({nVertices, nIntP});
    fill_levels(pi, nVertices, 8, kPress);
    FArray<real> press_in({nVertices, nVertLevels}), field_in({nVertices, nVertLevels});
    flip_in(press_in, pressure_v_a, nVertices, nVertLevels, 1.0);
    flip_in(field_in, d.vorticity, nVertices, nVertLevels, 1.0);
    interp_tofixed_pressure(nVertices, nVertLevels, nIntP, press_in, field_in, pi, fi);
    for (int l = 0; l < 8; ++l) copy_level(lev("vorticity", l), fi, nVertices, l + 1);
  }

  if (need.mslp) {  // IS:961-964: mslp(:) = mslp(:) * 100.0 over the whole array (the garbage element too)
    compute_slp(nCells, nVertLevels, temperature_a, m.zgrid, pressure_a, index_qv, s.scalars, mslp_a);
    const auto mslp = mslp_a.view();
    parallel_for_n(mslp_a.size(), [=] CPPMPAS_LAMBDA(std::size_t n) { mslp.data()[n] = mslp.data()[n] * 100.0; });
  }

  if (need.t_isobaric || need.meanT_500_300) {  // IS:973-1013
    if (!(need.temp || need.relhum || need.dewpoint || need.mslp))
      throw Error("isobaric diagnostics: t_isobaric / meanT_500_300 without temperature, relhum, dewpoint or mslp "
                  "read MPAS's uninitialised temperature array (IS:973-990); no port");
    FArray<real> field_in({nCells, nVertLevels}), press_in({nCells, nVertLevels});
    FArray<real> fi({nCells, nIsoLevelsT}), pi({nCells, nIsoLevelsT});
    {
      const auto piv = pi.view();
      parallel_for(1, nCells, [=] CPPMPAS_LAMBDA(int iCell) {
        for (int k = 1; k <= nIsoLevelsT; ++k) piv(iCell, k) = t_iso_levels(k);
      });
    }
    flip_in(field_in, temperature_a, nCells, nVertLevels, 1.0);
    flip_in(press_in, pressure_a, nCells, nVertLevels, 100.0);
    if (need.t_isobaric) {
      interp_tofixed_pressure(nCells, nVertLevels, nIsoLevelsT, press_in, field_in, pi, fi);
      const auto fiv = fi.cview();
      parallel_for(1, nCells, [=] CPPMPAS_LAMBDA(int iCell) {
        for (int k = 1; k <= nIsoLevelsT; ++k) t_isobaric(k, iCell) = fiv(iCell, k);
      });
    }
    if (need.meanT_500_300) compute_layer_mean(meanT_500_300, 50000.0, 30000.0, field_in, press_in);
  }

  if (need.z_isobaric) {  // IS:1017-1046
    FArray<real> field_in({nCells, nVertLevelsP1}), press_in({nCells, nVertLevelsP1});
    FArray<real> fi({nCells, nIsoLevelsZ}), pi({nCells, nIsoLevelsZ});
    {
      const auto piv = pi.view();
      parallel_for(1, nCells, [=] CPPMPAS_LAMBDA(int iCell) {
        for (int k = 1; k <= nIsoLevelsZ; ++k) piv(iCell, k) = z_iso_levels(k);
      });
    }
    flip_in(field_in, m.zgrid, nCells, nVertLevelsP1, 1.0);
    flip_in(press_in, pressure2_a, nCells, nVertLevelsP1, 100.0);
    interp_tofixed_pressure(nCells, nVertLevelsP1, nIsoLevelsZ, press_in, field_in, pi, fi);
    const auto fiv = fi.cview();
    parallel_for(1, nCells, [=] CPPMPAS_LAMBDA(int iCell) {
      for (int k = 1; k <= nIsoLevelsZ; ++k) z_isobaric(k, iCell) = fiv(iCell, k);
    });
  }
}

}  // namespace

// IS:62-306 isobaric_diagnostics_compute
void isobaric_diagnostics_compute(DiagnosticsContext& c) {
  Needs need;
  bool need_any_diags = false;
  need.mslp = c.field_will_be_written("mslp");
  need_any_diags = need_any_diags || need.mslp;
  // IS:106-297: the per-level flags, in MPAS's order; each level sets its group
  auto group = [&](const char* name, bool& g) {
    for (int l = 0; l < 8; ++l) {
      const bool n = c.field_will_be_written(std::string(name) + "_" + kLevels[l] + "hPa");
      g = g || n;
      need_any_diags = need_any_diags || n;
    }
  };
  group("relhum", need.relhum);
  group("dewpoint", need.dewpoint);
  group("temperature", need.temp);
  group("height", need.height);
  group("uzonal", need.uzonal);
  group("umeridional", need.umeridional);
  group("w", need.w);
  group("vorticity", need.vorticity);
  need.t_isobaric = c.field_will_be_written("t_isobaric");
  need_any_diags = need_any_diags || need.t_isobaric;
  need.z_isobaric = c.field_will_be_written("z_isobaric");
  need_any_diags = need_any_diags || need.z_isobaric;
  need.meanT_500_300 = c.field_will_be_written("meanT_500_300");
  need_any_diags = need_any_diags || need.meanT_500_300;
  if (need_any_diags) interp_diagnostics(c, need);
}

}  // namespace cppmpas::diag
