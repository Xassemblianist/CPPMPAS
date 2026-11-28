// SPDX-License-Identifier: GPL-3.0-or-later
// advance_scalars, atm_advance_scalars_work and atm_advance_scalars_mono_work of
// src/core_atmosphere/dynamics/mpas_atm_time_integration.F (TI); see atm_advance_scalars.hpp. Statement by
// statement, in MPAS's loop and operation order (a*b*c is (a*b)*c, - a*b is -(a*b)).
//
// Reference object code (objdump -d of bin/atmosphere_model, mpas-build asset 153c7206..., gfortran 14.2 -O3, no
// -march): atm_advance_scalars_work (0x120510) and atm_advance_scalars_mono_work (0x1145d0) call no libm or
// libmvec function and have no FMA and no ymm register. gfortran lowers MAX / MIN to MAX_EXPR / MIN_EXPR, which
// x86-64 emits as maxsd / minsd (maxpd / minpd) with the variable as the destination: max(0.0_RKIND, x) is
// "x > 0 ? x : 0" and min(0.0_RKIND, x) is "x < 0 ? x : 0" (e.g. the limiter TI:5606-5611 at 0x11b145-0x11b163,
// maxsd %xmm2(0.0),%xmm1(sf) then minsd %xmm6(1.0),%xmm1; scale_arr TI:5516-5519 at 0x11affa / 0x11b074), so a
// -0.0 or a NaN argument gives +0.0. fmax0 / fmin0 / fmin1 below are those instructions. The maxima and minima
// of two or more variables (s_max / s_min, the rescale factors) take their first argument as the running value;
// they only differ from the instructions in the sign of a zero between equal zeros, which no result depends on
// (the s_max / s_min products are clipped by the limiter, the rescale factors are >= +0).
//
// Every loop over cells or edges is a kernel (core/backend.hpp) with one thread per element; the level, scalar and
// neighbour loops stay inside the thread in MPAS's order, and MPAS's per-routine work arrays (scalar_weight2,
// scalar_tend_column, wdtn, flux_upwind_arr, ica, swa) are the thread's own. The loop over the scalars of the
// monotonic transport (TI:5282) is the host's, each pass over the mesh one launch.
#include "core_atmosphere/atm_advance_scalars.hpp"

#include <cmath>
#include <string>
#include <vector>

#include "core/error.hpp"
#include "core_atmosphere/atm_views.hpp"

namespace cppmpas {

namespace {

// flux4 / flux3 statement functions (TI:4730-4735, 5177-5182)
CPPMPAS_HD inline real flux4(real q_im2, real q_im1, real q_i, real q_ip1, real ua) {
  return ua * (7. * (q_i + q_im1) - (q_ip1 + q_im2)) / 12.0;
}
CPPMPAS_HD inline real flux3(real q_im2, real q_im1, real q_i, real q_ip1, real ua, real coef3) {
  return flux4(q_im2, q_im1, q_i, q_ip1, ua) + coef3 * std::fabs(ua) * ((q_ip1 - q_im2) - 3. * (q_i - q_im1)) / 12.0;
}
CPPMPAS_HD inline real sign1(real x) { return std::copysign(1.0, x); }  // sign(1.0_RKIND, x)

// maxsd / minsd with the variable as destination (see the header comment)
CPPMPAS_HD inline real fmax0(real x) { return x > 0.0 ? x : 0.0; }  // max(0.0_RKIND, x)
CPPMPAS_HD inline real fmin0(real x) { return x < 0.0 ? x : 0.0; }  // min(0.0_RKIND, x)
CPPMPAS_HD inline real fmin1(real x) { return x < 1.0 ? x : 1.0; }  // min(1.0_RKIND, x)
// Fortran MAX / MIN of variables, left to right
CPPMPAS_HD inline real vmax(real a, real b) { return b > a ? b : a; }
CPPMPAS_HD inline real vmin(real a, real b) { return b < a ? b : a; }

constexpr real eps = 1.e-20;  // TI:5175
constexpr int SCALE_IN = 1, SCALE_OUT = 2;  // TI:5174
// the capacity of the GPU's per-thread scalar work arrays (scalar_tend_column, wdtn of atm_advance_scalars_work):
// the scalar sets of the suites are {qv}, {qv, qc, qr}, mesoscale_reference's {qv, qc, qr, qi, qs, qg} and
// convection_permitting's {qv, qc, qr, qi, qs, qg, ni, nr}. A larger count is refused on both builds (the CPU's work
// arrays would take it; one range for both builds)
constexpr int kMaxScalars = 8;

// TI:4683-4939; edges 1..nEdges, cells 1..nCells. advance_density: config_split_dynamics_transport (TI:3049)
void atm_advance_scalars_work(const DycoreConfig& cfg, const Mesh& m, const State& s, DycoreFields& f, real dt,
                              FArray<real>& horiz_flux_arr_a, int rk_step, bool advance_density) {
  const int nVertLevels = m.nVertLevels;
  const int num_scalars = s.scalars.extent(0);
  const auto scalar_old = s.scalars.cview(), rho_zz_old = s.rho_zz.cview(), rho_zz_new = f.rho_zz2.cview();
  const auto scalar_new = f.scalars2.view();
  const auto uhAvg = f.ruAvg.cview(), wwAvg = f.wwAvg.cview();
  const auto scalar_tend_save = f.tend_scalars.cview();
  const auto fnm = m.fzm.cview(), fnp = m.fzp.cview(), rdnw = m.rdzw.cview();
  const auto horiz_flux_arr = horiz_flux_arr_a.view();
  const real coef_3rd_order = cfg.config_coef_3rd_order;

  real weight_time_new = 1.;  // :4749-4756
  if (!advance_density) {
    weight_time_new = 1.;                     // :4750
  } else {
    if (rk_step == 1 && cfg.config_time_integration_order == 3) weight_time_new = 1. / 3;
    if (rk_step == 1 && cfg.config_time_integration_order == 2) weight_time_new = 1. / 2;
    if (rk_step == 2) weight_time_new = 1. / 2;
    if (rk_step == 3) weight_time_new = 1.;
  }
  const real weight_time_old = 1. - weight_time_new;  // :4757

  parallel_for(1, m.nEdges, [=, m = view(m)] CPPMPAS_LAMBDA(int iEdge) {  // :4762-4844 (no LBCs: the full flux of :4764-4823)
    LocalArray<real, kMaxColumn * 10> scalar_weight2(nVertLevels * 10 + 1);
    auto sw2 = [&](int k, int j) -> real& { return scalar_weight2[(j - 1) * nVertLevels + k]; };
    int ica[11];
    if (m.nAdvCellsForEdge(iEdge) != 10) {  // :4801-4821 (case default: the edges next to a pentagon)
      for (int k = 1; k <= nVertLevels; ++k)
        for (int iScalar = 1; iScalar <= num_scalars; ++iScalar) horiz_flux_arr(iScalar, k, iEdge) = 0.0;
      for (int j = 1; j <= m.nAdvCellsForEdge(iEdge); ++j) {
        const int iAdvCell = m.advCellsForEdge(j, iEdge);
        for (int k = 1; k <= nVertLevels; ++k)
          for (int iScalar = 1; iScalar <= num_scalars; ++iScalar) {
            const real scalar_weight = m.adv_coefs(j, iEdge) + sign1(uhAvg(k, iEdge)) * m.adv_coefs_3rd(j, iEdge);
            horiz_flux_arr(iScalar, k, iEdge) =
                horiz_flux_arr(iScalar, k, iEdge) + scalar_weight * scalar_new(iScalar, k, iAdvCell);
          }
      }
      return;
    }
    for (int j = 1; j <= 10; ++j)  // :4772-4776
      for (int k = 1; k <= nVertLevels; ++k)
        sw2(k, j) = m.adv_coefs(j, iEdge) + sign1(uhAvg(k, iEdge)) * m.adv_coefs_3rd(j, iEdge);
    for (int j = 1; j <= 10; ++j) ica[j] = m.advCellsForEdge(j, iEdge);  // :4778-4780
    for (int k = 1; k <= nVertLevels; ++k)  // :4784-4798
      for (int iScalar = 1; iScalar <= num_scalars; ++iScalar)
        horiz_flux_arr(iScalar, k, iEdge) =
            sw2(k, 1) * scalar_new(iScalar, k, ica[1]) + sw2(k, 2) * scalar_new(iScalar, k, ica[2]) +
            sw2(k, 3) * scalar_new(iScalar, k, ica[3]) + sw2(k, 4) * scalar_new(iScalar, k, ica[4]) +
            sw2(k, 5) * scalar_new(iScalar, k, ica[5]) + sw2(k, 6) * scalar_new(iScalar, k, ica[6]) +
            sw2(k, 7) * scalar_new(iScalar, k, ica[7]) + sw2(k, 8) * scalar_new(iScalar, k, ica[8]) +
            sw2(k, 9) * scalar_new(iScalar, k, ica[9]) + sw2(k, 10) * scalar_new(iScalar, k, ica[10]);
  });

  parallel_for(1, m.nCellsSolve(), [=, m = view(m)] CPPMPAS_LAMBDA(int iCell) {  // :4859-4932 (bdyMaskCell = 0 <= nRelaxZone)
    LocalArray<real, (kMaxScalars + 1) * (kMaxColumn - 1)> stc((num_scalars + 1) * (nVertLevels + 1));
    LocalArray<real, (kMaxScalars + 1) * kMaxColumn> wdtn((num_scalars + 1) * (nVertLevels + 2));
    auto STC = [&](int q, int k) -> real& { return stc[k * (num_scalars + 1) + q]; };
    auto WDTN = [&](int q, int k) -> real& { return wdtn[k * (num_scalars + 1) + q]; };
    for (int k = 1; k <= nVertLevels; ++k)  // :4863-4871 (DO_PHYSICS: scalar_tend_save is kept)
      for (int iScalar = 1; iScalar <= num_scalars; ++iScalar) STC(iScalar, k) = 0.0;
    for (int i = 1; i <= m.nEdgesOnCell(iCell); ++i) {  // :4873-4885
      const int iEdge = m.edgesOnCell(i, iCell);
      for (int k = 1; k <= nVertLevels; ++k)
        for (int iScalar = 1; iScalar <= num_scalars; ++iScalar)
          STC(iScalar, k) = STC(iScalar, k) - m.edgesOnCell_sign(i, iCell) * uhAvg(k, iEdge) * horiz_flux_arr(iScalar, k, iEdge);
    }
    for (int k = 1; k <= nVertLevels; ++k)  // :4887-4892
      for (int iScalar = 1; iScalar <= num_scalars; ++iScalar)
        STC(iScalar, k) = STC(iScalar, k) * m.invAreaCell(iCell) + scalar_tend_save(iScalar, k, iCell);
    for (int iScalar = 1; iScalar <= num_scalars; ++iScalar) {  // :4899-4906
      WDTN(iScalar, 1) = 0.0;
      WDTN(iScalar, 2) = wwAvg(2, iCell) * (fnm(2) * scalar_new(iScalar, 2, iCell) + fnp(2) * scalar_new(iScalar, 2 - 1, iCell));
      WDTN(iScalar, nVertLevels) =
          wwAvg(nVertLevels, iCell) * (fnm(nVertLevels) * scalar_new(iScalar, nVertLevels, iCell) +
                                       fnp(nVertLevels) * scalar_new(iScalar, nVertLevels - 1, iCell));
      WDTN(iScalar, nVertLevels + 1) = 0.0;
    }
    for (int k = 3; k <= nVertLevels - 1; ++k)  // :4908-4914
      for (int iScalar = 1; iScalar <= num_scalars; ++iScalar)
        WDTN(iScalar, k) = flux3(scalar_new(iScalar, k - 2, iCell), scalar_new(iScalar, k - 1, iCell),
                                 scalar_new(iScalar, k, iCell), scalar_new(iScalar, k + 1, iCell), wwAvg(k, iCell),
                                 coef_3rd_order);
    for (int k = 1; k <= nVertLevels; ++k)  // :4920-4928
      for (int iScalar = 1; iScalar <= num_scalars; ++iScalar) {
        const real rho_zz_new_inv =
            1.0 / (weight_time_old * rho_zz_old(k, iCell) + weight_time_new * rho_zz_new(k, iCell));
        scalar_new(iScalar, k, iCell) =
            (scalar_old(iScalar, k, iCell) * rho_zz_old(k, iCell) +
             dt * (STC(iScalar, k) - rdnw(k) * (WDTN(iScalar, k + 1) - WDTN(iScalar, k)))) *
            rho_zz_new_inv;
      }
  });
}

struct MonoWork {
  FArray<real> scalar_old, scalar_new, s_max, s_min, flux_arr, wdtn, flux_upwind_tmp, flux_tmp, scale_arr;
  FArray<real> rho_zz_int;  // module array of TI:71, allocated by advance_scalars (TI:3009-3014)
};

// TI:5106-5783; cells / edges are the whole ranges. advance_density: config_split_dynamics_transport (TI:3059)
void atm_advance_scalars_mono_work(const DycoreConfig& cfg, const Mesh& m, State& s, DycoreFields& f, real dt,
                                   MonoWork& w, bool advance_density) {
  const int nVertLevels = m.nVertLevels, nCells = m.nCells, nCellsSolve = m.nCellsSolve();
  const int num_scalars = s.scalars.extent(0);
  const auto scalars_old = s.scalars.view(), scalars_new = f.scalars2.view(), scalar_tend = f.tend_scalars.view();
  const auto rho_zz_old = s.rho_zz.cview(), uhAvg = f.ruAvg.cview(),
             wwAvg = f.wwAvg.cview();
  const auto fnm = m.fzm.cview(), fnp = m.fzp.cview(), rdnw = m.rdzw.cview();
  const real coef_3rd_order = cfg.config_coef_3rd_order;
  const auto scalar_old = w.scalar_old.view(), scalar_new = w.scalar_new.view(), s_max = w.s_max.view(),
             s_min = w.s_min.view();
  const auto flux_arr = w.flux_arr.view(), wdtn = w.wdtn.view(), flux_upwind_tmp = w.flux_upwind_tmp.view(),
             flux_tmp = w.flux_tmp.view();
  const auto scale_arr = w.scale_arr.view();

  parallel_for(1, nCellsSolve, [=] CPPMPAS_LAMBDA(int iCell) {  // :5199-5213 (DO_PHYSICS: scalar_tend is not zeroed first)
    for (int k = 1; k <= nVertLevels; ++k)
      for (int iScalar = 1; iScalar <= num_scalars; ++iScalar) {
        scalars_old(iScalar, k, iCell) =
            scalars_old(iScalar, k, iCell) + dt * scalar_tend(iScalar, k, iCell) / rho_zz_old(k, iCell);
        scalar_tend(iScalar, k, iCell) = 0.0;
      }
  });
  // :5223 halo exchange 'dynamics:scalars_old': no-op on one rank
  const auto rho_zz_int = w.rho_zz_int.view();
  if (advance_density) {  // :5235-5278
    parallel_for(1, nCellsSolve, [=, m = view(m)] CPPMPAS_LAMBDA(int iCell) {
      for (int k = 1; k <= nVertLevels; ++k) rho_zz_int(k, iCell) = 0.0;
      for (int i = 1; i <= m.nEdgesOnCell(iCell); ++i) {
        const int iEdge = m.edgesOnCell(i, iCell);
        for (int k = 1; k <= nVertLevels; ++k)
          rho_zz_int(k, iCell) = rho_zz_int(k, iCell) - m.edgesOnCell_sign(i, iCell) * uhAvg(k, iEdge) *
                                                            m.dvEdge(iEdge) * m.invAreaCell(iCell);
      }
    });
    parallel_for(1, nCellsSolve, [=] CPPMPAS_LAMBDA(int iCell) {
      for (int k = 1; k <= nVertLevels; ++k)
        rho_zz_int(k, iCell) =
            rho_zz_old(k, iCell) + dt * (rho_zz_int(k, iCell) - rdnw(k) * (wwAvg(k + 1, iCell) - wwAvg(k, iCell)));
    });
  }
  // :5585-5615, 5725-5735: rho_zz_int with advance_density, else rho_zz_new
  const auto rho_lim = advance_density ? w.rho_zz_int.cview() : f.rho_zz2.cview();

  for (int iScalar = 1; iScalar <= num_scalars; ++iScalar) {  // :5282-5781
    parallel_for(1, nCells, [=] CPPMPAS_LAMBDA(int iCell) {  // :5284-5289
      for (int k = 1; k <= nVertLevels; ++k) {
        scalar_old(k, iCell) = scalars_old(iScalar, k, iCell);
        scalar_new(k, iCell) = scalars_new(iScalar, k, iCell);
      }
    });
    parallel_for(1, nVertLevels, [=] CPPMPAS_LAMBDA(int k) {  // :5297-5300
      scalar_old(k, nCells + 1) = 0.0;
      scalar_new(k, nCells + 1) = 0.0;
    });

    parallel_for(1, nCellsSolve, [=, m = view(m)] CPPMPAS_LAMBDA(int iCell) {  // :5332-5406
      wdtn(1, iCell) = 0.0;
      wdtn(nVertLevels + 1, iCell) = 0.0;
      int k = 1;
      s_max(k, iCell) = vmax(scalar_old(1, iCell), scalar_old(2, iCell));
      s_min(k, iCell) = vmin(scalar_old(1, iCell), scalar_old(2, iCell));
      k = 2;
      wdtn(k, iCell) = wwAvg(k, iCell) * (fnm(k) * scalar_new(k, iCell) + fnp(k) * scalar_new(k - 1, iCell));
      s_max(k, iCell) = vmax(vmax(scalar_old(k - 1, iCell), scalar_old(k, iCell)), scalar_old(k + 1, iCell));
      s_min(k, iCell) = vmin(vmin(scalar_old(k - 1, iCell), scalar_old(k, iCell)), scalar_old(k + 1, iCell));
      for (k = 3; k <= nVertLevels - 1; ++k) {
        wdtn(k, iCell) = flux3(scalar_new(k - 2, iCell), scalar_new(k - 1, iCell), scalar_new(k, iCell),
                               scalar_new(k + 1, iCell), wwAvg(k, iCell), coef_3rd_order);
        s_max(k, iCell) = vmax(vmax(scalar_old(k - 1, iCell), scalar_old(k, iCell)), scalar_old(k + 1, iCell));
        s_min(k, iCell) = vmin(vmin(scalar_old(k - 1, iCell), scalar_old(k, iCell)), scalar_old(k + 1, iCell));
      }
      k = nVertLevels;
      wdtn(k, iCell) = wwAvg(k, iCell) * (fnm(k) * scalar_new(k, iCell) + fnp(k) * scalar_new(k - 1, iCell));
      s_max(k, iCell) = vmax(scalar_old(k, iCell), scalar_old(k - 1, iCell));
      s_min(k, iCell) = vmin(scalar_old(k, iCell), scalar_old(k - 1, iCell));
      if (m.nEdgesOnCell(iCell) != 6) {  // :5394-5404 (case default: the pentagons)
        for (int i = 1; i <= m.nEdgesOnCell(iCell); ++i)
          for (k = 1; k <= nVertLevels; ++k) {
            s_max(k, iCell) = vmax(s_max(k, iCell), scalar_old(k, m.cellsOnCell(i, iCell)));
            s_min(k, iCell) = vmin(s_min(k, iCell), scalar_old(k, m.cellsOnCell(i, iCell)));
          }
        return;
      }
      for (k = 1; k <= nVertLevels; ++k) {  // :5374-5392
        real mx = s_max(k, iCell), mn = s_min(k, iCell);
        for (int j = 1; j <= 6; ++j) {
          mx = vmax(mx, scalar_old(k, m.cellsOnCell(j, iCell)));
          mn = vmin(mn, scalar_old(k, m.cellsOnCell(j, iCell)));
        }
        s_max(k, iCell) = mx;
        s_min(k, iCell) = mn;
      }
    });

    parallel_for(1, m.nEdges, [=, m = view(m)] CPPMPAS_LAMBDA(int iEdge) {  // :5414-5471
      int ica[11];
      real swa[11][3];
      const int cell1 = m.cellsOnEdge(1, iEdge);
      const int cell2 = m.cellsOnEdge(2, iEdge);
      if (cell1 <= nCellsSolve || cell2 <= nCellsSolve) {
        if (m.nAdvCellsForEdge(iEdge) != 10) {  // :5450-5463 (case default)
          for (int k = 1; k <= nVertLevels; ++k) flux_arr(k, iEdge) = 0.0;
          for (int i = 1; i <= m.nAdvCellsForEdge(iEdge); ++i) {
            const int iCell = m.advCellsForEdge(i, iEdge);
            for (int k = 1; k <= nVertLevels; ++k) {
              const real scalar_weight =
                  uhAvg(k, iEdge) * (m.adv_coefs(i, iEdge) + sign1(uhAvg(k, iEdge)) * m.adv_coefs_3rd(i, iEdge));
              flux_arr(k, iEdge) = flux_arr(k, iEdge) + scalar_weight * scalar_new(k, iCell);
            }
          }
          return;
        }
        for (int jj = 1; jj <= 10; ++jj) {  // :5429-5433
          ica[jj] = m.advCellsForEdge(jj, iEdge);
          swa[jj][1] = m.adv_coefs(jj, iEdge) + m.adv_coefs_3rd(jj, iEdge);
          swa[jj][2] = m.adv_coefs(jj, iEdge) - m.adv_coefs_3rd(jj, iEdge);
        }
        for (int k = 1; k <= nVertLevels; ++k) {  // :5435-5446
          const int ii = uhAvg(k, iEdge) > 0 ? 1 : 2;  // merge(1, 2, uhAvg(k,iEdge) > 0)
          flux_arr(k, iEdge) =
              uhAvg(k, iEdge) *
              (swa[1][ii] * scalar_new(k, ica[1]) + swa[2][ii] * scalar_new(k, ica[2]) +
               swa[3][ii] * scalar_new(k, ica[3]) + swa[4][ii] * scalar_new(k, ica[4]) +
               swa[5][ii] * scalar_new(k, ica[5]) + swa[6][ii] * scalar_new(k, ica[6]) +
               swa[7][ii] * scalar_new(k, ica[7]) + swa[8][ii] * scalar_new(k, ica[8]) +
               swa[9][ii] * scalar_new(k, ica[9]) + swa[10][ii] * scalar_new(k, ica[10]));
        }
      } else {
        for (int k = 1; k <= nVertLevels; ++k) flux_arr(k, iEdge) = 0.0;
      }
    });

    parallel_for(1, nCellsSolve, [=] CPPMPAS_LAMBDA(int iCell) {  // :5482-5521
      Column<real> flux_upwind_arr(nVertLevels + 2);
      int k = 1;
      scalar_new(k, iCell) = scalar_old(k, iCell) * rho_zz_old(k, iCell);
      for (k = 2; k <= nVertLevels; ++k) {
        scalar_new(k, iCell) = scalar_old(k, iCell) * rho_zz_old(k, iCell);
        flux_upwind_arr[k] = dt * (fmax0(wwAvg(k, iCell)) * scalar_old(k - 1, iCell) +
                                   fmin0(wwAvg(k, iCell)) * scalar_old(k, iCell));
      }
      for (k = 1; k <= nVertLevels - 1; ++k)
        scalar_new(k, iCell) = scalar_new(k, iCell) - flux_upwind_arr[k + 1] * rdnw(k);
      for (k = 2; k <= nVertLevels; ++k) {
        scalar_new(k, iCell) = scalar_new(k, iCell) + flux_upwind_arr[k] * rdnw(k);
        wdtn(k, iCell) = dt * wdtn(k, iCell) - flux_upwind_arr[k];
      }
      for (k = 1; k <= nVertLevels; ++k) {  // :5516-5519
        scale_arr(k, SCALE_IN, iCell) = -(rdnw(k) * (fmin0(wdtn(k + 1, iCell)) - fmax0(wdtn(k, iCell))));
        scale_arr(k, SCALE_OUT, iCell) = -(rdnw(k) * (fmax0(wdtn(k + 1, iCell)) - fmin0(wdtn(k, iCell))));
      }
    });

    parallel_for(1, m.nEdges, [=, m = view(m)] CPPMPAS_LAMBDA(int iEdge) {  // :5531-5548
      const int cell1 = m.cellsOnEdge(1, iEdge);
      const int cell2 = m.cellsOnEdge(2, iEdge);
      for (int k = 1; k <= nVertLevels; ++k) {
        flux_upwind_tmp(k, iEdge) = m.dvEdge(iEdge) * dt *
                                    (fmax0(uhAvg(k, iEdge)) * scalar_old(k, cell1) +
                                     fmin0(uhAvg(k, iEdge)) * scalar_old(k, cell2));
        flux_tmp(k, iEdge) = dt * flux_arr(k, iEdge) - flux_upwind_tmp(k, iEdge);
      }
      // :5541 (config_apply_lbcs .and. mask == nRelaxZone) .or. (mask == nRelaxZone-1): bdyMaskEdge = 0, nRelaxZone = 5
    });

    parallel_for(1, nCellsSolve, [=, m = view(m)] CPPMPAS_LAMBDA(int iCell) {  // :5561-5574
      for (int i = 1; i <= m.nEdgesOnCell(iCell); ++i) {
        const int iEdge = m.edgesOnCell(i, iCell);
        for (int k = 1; k <= nVertLevels; ++k) {
          scalar_new(k, iCell) =
              scalar_new(k, iCell) - m.edgesOnCell_sign(i, iCell) * flux_upwind_tmp(k, iEdge) * m.invAreaCell(iCell);
          scale_arr(k, SCALE_OUT, iCell) =
              scale_arr(k, SCALE_OUT, iCell) - fmax0(m.edgesOnCell_sign(i, iCell) * flux_tmp(k, iEdge)) * m.invAreaCell(iCell);
          scale_arr(k, SCALE_IN, iCell) =
              scale_arr(k, SCALE_IN, iCell) - fmin0(m.edgesOnCell_sign(i, iCell) * flux_tmp(k, iEdge)) * m.invAreaCell(iCell);
        }
      }
    });

    parallel_for(1, nCellsSolve, [=] CPPMPAS_LAMBDA(int iCell) {  // :5585-5615
      for (int k = 1; k <= nVertLevels; ++k) {
        real scale_factor = (s_max(k, iCell) * rho_lim(k, iCell) - scalar_new(k, iCell)) /
                            (scale_arr(k, SCALE_IN, iCell) + eps);
        scale_arr(k, SCALE_IN, iCell) = fmin1(fmax0(scale_factor));
        scale_factor = (s_min(k, iCell) * rho_lim(k, iCell) - scalar_new(k, iCell)) /
                       (scale_arr(k, SCALE_OUT, iCell) - eps);
        scale_arr(k, SCALE_OUT, iCell) = fmin1(fmax0(scale_factor));
      }
    });
    // :5630 halo exchange 'dynamics:scale': no-op on one rank

    parallel_for(1, m.nEdges, [=, m = view(m)] CPPMPAS_LAMBDA(int iEdge) {  // :5640-5664
      const int cell1 = m.cellsOnEdge(1, iEdge);
      const int cell2 = m.cellsOnEdge(2, iEdge);
      if (cell1 <= nCellsSolve || cell2 <= nCellsSolve)
        for (int k = 1; k <= nVertLevels; ++k) {
          const real flux_upwind = m.dvEdge(iEdge) * dt *
                                   (fmax0(uhAvg(k, iEdge)) * scalar_old(k, cell1) +
                                    fmin0(uhAvg(k, iEdge)) * scalar_old(k, cell2));
          flux_arr(k, iEdge) = dt * flux_arr(k, iEdge) - flux_upwind;
        }
      // :5654 the LBC test, false (as :5541)
    });

    parallel_for(1, m.nEdges, [=, m = view(m)] CPPMPAS_LAMBDA(int iEdge) {  // :5672-5686
      const int cell1 = m.cellsOnEdge(1, iEdge);
      const int cell2 = m.cellsOnEdge(2, iEdge);
      if (cell1 <= nCellsSolve || cell2 <= nCellsSolve)
        for (int k = 1; k <= nVertLevels; ++k) {
          real flux = flux_arr(k, iEdge);
          flux = fmax0(flux) * vmin(scale_arr(k, SCALE_OUT, cell1), scale_arr(k, SCALE_IN, cell2)) +
                 fmin0(flux) * vmin(scale_arr(k, SCALE_IN, cell1), scale_arr(k, SCALE_OUT, cell2));
          flux_arr(k, iEdge) = flux;
        }
    });

    parallel_for(1, nCellsSolve, [=] CPPMPAS_LAMBDA(int iCell) {  // :5700-5707
      for (int k = 2; k <= nVertLevels; ++k) {
        real flux = wdtn(k, iCell);
        flux = fmax0(flux) * vmin(scale_arr(k - 1, SCALE_OUT, iCell), scale_arr(k, SCALE_IN, iCell)) +
               fmin0(flux) * vmin(scale_arr(k, SCALE_OUT, iCell), scale_arr(k - 1, SCALE_IN, iCell));
        wdtn(k, iCell) = flux;
      }
    });

    parallel_for(1, nCellsSolve, [=, m = view(m)] CPPMPAS_LAMBDA(int iCell) {  // :5714-5739
      for (int i = 1; i <= m.nEdgesOnCell(iCell); ++i) {
        const int iEdge = m.edgesOnCell(i, iCell);
        for (int k = 1; k <= nVertLevels; ++k)
          scalar_new(k, iCell) = scalar_new(k, iCell) - m.edgesOnCell_sign(i, iCell) * flux_arr(k, iEdge) * m.invAreaCell(iCell);
      }
      for (int k = 1; k <= nVertLevels; ++k)  // :5725-5735 (rho_zz_int with advance_density)
        scalar_new(k, iCell) =
            (scalar_new(k, iCell) + (-(rdnw(k) * (wdtn(k + 1, iCell) - wdtn(k, iCell))))) / rho_lim(k, iCell);
    });

    parallel_for(1, nCells, [=] CPPMPAS_LAMBDA(int iCell) {  // :5770-5777 (bdyMaskCell = 0 <= nSpecZone)
      for (int k = 1; k <= nVertLevels; ++k) scalars_new(iScalar, k, iCell) = fmax0(scalar_new(k, iCell));
    });
  }
}

}  // namespace

void check_scalar_capacity(int num_scalars) {
  if (num_scalars > kMaxScalars)  // both builds (kMaxScalars)
    throw Error("advance_scalars: more than " + std::to_string(kMaxScalars) +
                " scalars (the capacity of the per-thread work arrays)");
}

void check_advance_scalars_supported(const DycoreConfig& cfg, const Mesh& m) {
  auto no = [](const std::string& what) { throw Error("advance_scalars: " + what + " has no verified port"); };
  if (!cfg.config_monotonic) no("config_monotonic = false (the non-monotonic rk_step 3, TI:3044)");
  if (cfg.config_positive_definite) no("config_positive_definite = true");
  if (cfg.config_apply_lbcs) no("config_apply_lbcs = true");
  if (m.nVertLevels < 4) no("nVertLevels < 4");
}

// TI:2911-3085 (rk_step < 3: atm_advance_scalars; rk_step 3 and config_monotonic: atm_advance_scalars_mono)
void advance_scalars(const DycoreConfig& cfg, const Mesh& m, State& s, DycoreFields& f, int rk_step,
                     const real* rk_timestep) {
  check_advance_scalars_supported(cfg, m);
  const int nCells = m.nCells, nEdges = m.nEdges, nVertLevels = m.nVertLevels;
  const int num_scalars = s.scalars.extent(0);
  fits_column(nVertLevels, "advance_scalars");
  check_scalar_capacity(num_scalars);
  if (rk_step < 3) {  // :3021-3029, 3044-3049
    FArray<real> horiz_flux_array({num_scalars, nVertLevels, nEdges + 1}, 0.0);
    atm_advance_scalars_work(cfg, m, s, f, rk_timestep[rk_step], horiz_flux_array, rk_step,
                             cfg.config_split_dynamics_transport);
  } else {  // :2993-3005, 3031-3039, 3050-3059 (the scratch arrays' garbage elements zeroed)
    MonoWork w;
    w.scalar_old.allocate({nVertLevels, nCells + 1}, 0.0);
    w.scalar_new.allocate({nVertLevels, nCells + 1}, 0.0);
    w.s_max.allocate({nVertLevels, nCells + 1}, 0.0);
    w.s_min.allocate({nVertLevels, nCells + 1}, 0.0);
    w.flux_arr.allocate({nVertLevels, nEdges + 1}, 0.0);
    w.wdtn.allocate({nVertLevels + 1, nCells + 1}, 0.0);
    w.flux_upwind_tmp.allocate({nVertLevels, nEdges + 1}, 0.0);
    w.flux_tmp.allocate({nVertLevels, nEdges + 1}, 0.0);
    w.scale_arr.allocate({nVertLevels, 2, nCells + 1}, 0.0);  // halo_scratch 'scale' (mpas_allocate_scratch_field)
    if (cfg.config_split_dynamics_transport)
      w.rho_zz_int.allocate({nVertLevels, nCells + 1}, 0.0);  // TI:3010-3011 (garbage column zeroed)
    else
      w.rho_zz_int.allocate({1, 1}, 0.0);
    atm_advance_scalars_mono_work(cfg, m, s, f, rk_timestep[rk_step], w, cfg.config_split_dynamics_transport);
  }
}

}  // namespace cppmpas
