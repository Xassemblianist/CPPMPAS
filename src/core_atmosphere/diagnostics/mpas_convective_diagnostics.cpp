// SPDX-License-Identifier: GPL-3.0-or-later
// mpas_convective_diagnostics.F (CV), statement by statement: see mpas_atm_diagnostics.hpp.
// With -fdefault-real-8 (the reference's FFLAGS) every literal without a kind (0.5, 0.0001, 611.2, 56.0, ...) and every
// parameter is a double; gfortran folds the parameter expressions (lv1, rp00, rddcp, ...) one correctly rounded
// operation per node, as the constexpr expressions below. Integer powers (**2, **(-1)) are a multiply and 1/x.
// The libm calls, as the reference executable makes them (objdump of bin/atmosphere_model,
// __mpas_convective_diagnostics_MOD_convective_diagnostics_compute at 0x236dc0, getcape and getqvs / getqvi / getthe
// inlined into it):
//  * CV:363-374, the temperature / dewpoint loop: the k loop is vectorised in pairs (k = 1,2 / 3,4 / ...) with one
//    _ZGVbN2v_log of evp/6.112 (both logs of the statement are one call); an odd last level calls log@plt;
//  * CV:709-717, getcape's first loop: pairs of k with _ZGVbN2vv_pow ((p*rp00)**rddcp) and _ZGVbN2v_exp (getqvs),
//    an odd last level pow@plt / exp@plt;
//  * everything else of getcape (getthe's log / pow / exp, the parcel ascent's pow, log, and the exp of getqvs,
//    getqvi and th2) is scalar: log@plt, pow@plt, exp@plt.
// convective_diagnostics_update (0x23c940) calls no libm: sqrtsd / sqrtpd; MAX(old, new) is maxsd with the new value
// as the destination (new > old ? new : old); MAX(0., x) is x > 0 ? x : 0.
#include <cmath>
#include <string>
#include <vector>

#include "core/error.hpp"
#include "core/libm.hpp"
#include "core_atmosphere/diagnostics/mpas_atm_diagnostics.hpp"
#include "physics/mpas_atmphys_constants.hpp"

namespace cppmpas::diag {

namespace {

inline real maxsd(real a, real b) { return a > b ? a : b; }  // gfortran MAX as maxsd, a the destination
inline real minsd(real a, real b) { return a < b ? a : b; }

// CV:528-544 (z, column_values 1-based through the pointers' [-1])
real column_height_value(const real* column_values, const real* z, real z_interp, int n) {
  int kz = 1;
  for (int k = 1; k <= n; ++k)
    if (z[k - 1] <= z_interp) kz = k;
  kz = kz < n - 1 ? kz : n - 1;
  const real wz = (z[kz] - z_interp) / (z[kz] - z[kz - 1]);
  const real wzp1 = 1. - wz;
  return wz * column_values[kz - 1] + wzp1 * column_values[kz];
}

// CV:548-565
real integral_zstaggered(const real* column_values, const real* z, real zbot, real ztop, int n) {
  real s = 0.;
  for (int k = 1; k <= n; ++k) {
    const real zb = maxsd(z[k - 1], zbot);
    const real zt = minsd(z[k], ztop);
    s = s + column_values[k - 1] * maxsd(zt - zb, 0.);
  }
  return s;
}

// CV:569-589
real integral_zpoint(const real* column_values, const real* z, real zbot, real ztop, int n) {
  real s = 0.;
  for (int k = 1; k <= n - 1; ++k) {
    const real zb = maxsd(z[k - 1], zbot);
    const real zt = minsd(z[k], ztop);
    const real dz = maxsd(zt - zb, 0.);
    const real zr_midpoint = (0.5 * (zt + zb) - z[k - 1]) / (z[k] - z[k - 1]);
    const real midpoint_value = column_values[k - 1] + (column_values[k] - column_values[k - 1]) * zr_midpoint;
    s = s + dz * midpoint_value;
  }
  return s;
}

// CV:1043-1054 getqvs (eps = 287.04/461.5 folded); vector: the exp lane of the first getcape loop
inline real getqvs_es(real t, bool lane) {
  const real a = 17.67 * (t - 273.15) / (t - 29.65);
  return 611.2 * (lane ? libm::mvec_exp(a) : libm::exp(a));
}
inline real getqvs(real p, real t, bool lane = false) {
  constexpr real eps = 287.04 / 461.5;
  const real es = getqvs_es(t, lane);
  return eps * es / (p - es);
}

// CV:1060-1071 getqvi
inline real getqvi(real p, real t) {
  constexpr real eps = 287.04 / 461.5;
  const real es = 611.2 * libm::exp(21.8745584 * (t - 273.15) / (t - 7.66));
  return eps * es / (p - es);
}

// CV:1077-1093 getthe
inline real getthe(real p, real t, real td, real q) {
  real tlcl;
  if ((td - t) >= -0.1) {
    tlcl = t;
  } else {
    tlcl = 56.0 + 1. / (1. / (td - 56.0) + 0.00125 * libm::log(t / td));
  }
  return t * libm::pow(100000.0 / p, 0.2854 * (1.0 - 0.28 * q)) *
         libm::exp(((3376.0 / tlcl) - 2.54) * q * (1.0 + 0.81 * q));
}

// CV:595-1037 getcape with source = 2 (most unstable parcel), adiabat = 1 (pseudoadiabatic, liquid only),
// debug_level = 0 (no log writes)
void getcape(int nk, const real* p_in, const real* t_in, const real* td_in, real& cape, real& cin) {
  constexpr real pinc = 100.0;
  constexpr real g = 9.81;
  constexpr real p00 = 100000.0;
  constexpr real cp = 1005.7;
  constexpr real rd = 287.04;
  constexpr real rv = 461.5;
  constexpr real xlv = 2501000.0;
  constexpr real xls = 2836017.0;
  constexpr real t0 = 273.15;
  constexpr real cpv = 1875.0;
  constexpr real cpl = 4190.0;
  constexpr real cpi = 2118.636;
  constexpr real lv1 = xlv + (cpl - cpv) * t0;
  constexpr real lv2 = cpl - cpv;
  constexpr real ls1 = xls + (cpi - cpv) * t0;
  constexpr real ls2 = cpi - cpv;
  constexpr real rp00 = 1.0 / p00;
  constexpr real reps = rv / rd;
  constexpr real rddcp = rd / cp;
  constexpr real cpdg = cp / g;
  constexpr real converge = 0.0002;
  constexpr real rd_cp = rd / cp;  // the constant rd/cp of CV:932

  std::vector<real> p(nk + 1), t(nk + 1), td(nk + 1), pi(nk + 1), q(nk + 1), th(nk + 1), thv(nk + 1);
  // CV:708-717; pairs of k use the libmvec lanes, an odd last k the scalar calls
  const int nvec = (nk / 2) * 2;
  for (int k = 1; k <= nk; ++k) {
    const bool lane = k <= nvec;
    p[k] = 100.0 * p_in[k - 1];
    t[k] = 273.15 + t_in[k - 1];
    td[k] = 273.15 + td_in[k - 1];
    pi[k] = lane ? libm::mvec_pow(p[k] * rp00, rddcp) : libm::pow(p[k] * rp00, rddcp);
    q[k] = getqvs(p[k], td[k], lane);
    th[k] = t[k] / pi[k];
    thv[k] = th[k] * (1.0 + reps * q[k]) / (1.0 + q[k]);
  }
  // CV:721-725: z (used by source = 3 only)

  // CV:733-754 source = 2
  int kmax = 1;
  real maxthe;
  if (p[1] < 50000.0) {
    kmax = 1;
    maxthe = getthe(p[1], t[1], td[1], q[1]);
  } else {
    maxthe = 0.0;
    for (int k = 1; k <= nk; ++k) {
      if (p[k] >= 50000.0) {
        const real the = getthe(p[k], t[k], td[k], q[k]);
        if (the > maxthe) {
          maxthe = the;
          kmax = k;
        }
      }
    }
  }
  (void)maxthe;

  // CV:820-855
  real narea = 0.0;
  int k = kmax;
  real th2 = th[kmax];
  real pi2 = pi[kmax];
  real p2 = p[kmax];
  real t2 = t[kmax];
  real thv2 = thv[kmax];
  real qv2 = q[kmax];
  real b2 = 0.0;
  real ql2 = 0.0;
  real qi2 = 0.0;
  real qt = qv2;
  cape = 0.0;
  cin = 0.0;
  bool doit = true;
  const bool ice = false;  // adiabat = 1

  // CV:866-1031 the ascent (getthe of the parcel, CV:857 and :961, only feeds debug output)
  while (doit && (k < nk)) {
    k = k + 1;
    const real b1 = b2;
    real dp = p[k - 1] - p[k];
    int nloop;
    if (dp < pinc) {
      nloop = 1;
    } else {
      nloop = 1 + static_cast<int>(dp / pinc);
      dp = dp / static_cast<real>(nloop);
    }
    for (int n = 1; n <= nloop; ++n) {
      const real p1 = p2;
      const real t1 = t2;
      const real th1 = th2;
      const real qv1 = qv2;
      const real ql1 = ql2;
      const real qi1 = qi2;
      p2 = p2 - dp;
      pi2 = libm::pow(p2 * rp00, rddcp);
      real thlast = th1;
      int i = 0;
      bool not_converged = true;
      while (not_converged) {
        i = i + 1;
        t2 = thlast * pi2;
        real fliq, fice;
        if (ice) {
          fliq = maxsd(minsd((t2 - 233.15) / (273.15 - 233.15), 1.0), 0.0);
          fice = 1.0 - fliq;
        } else {
          fliq = 1.0;
          fice = 0.0;
        }
        qv2 = minsd(qt, fliq * getqvs(p2, t2) + fice * getqvi(p2, t2));
        qi2 = maxsd(fice * (qt - qv2), 0.0);
        ql2 = maxsd(qt - qv2 - qi2, 0.0);
        const real tbar = 0.5 * (t1 + t2);
        const real qvbar = 0.5 * (qv1 + qv2);
        const real qlbar = 0.5 * (ql1 + ql2);
        const real qibar = 0.5 * (qi1 + qi2);
        const real lhv = lv1 - lv2 * tbar;
        const real lhs = ls1 - ls2 * tbar;
        const real rm = rd + rv * qvbar;
        const real cpm = cp + cpv * qvbar + cpl * qlbar + cpi * qibar;
        th2 = th1 * libm::exp(lhv * (ql2 - ql1) / (cpm * tbar) + lhs * (qi2 - qi1) / (cpm * tbar) +
                              (rm / cpm - rd_cp) * libm::log(p2 / p1));
        if (i > 100) return;  // CV:937-947 lack of convergence (debug_level = 0: no message)
        if (std::fabs(th2 - thlast) > converge) {
          thlast = thlast + 0.3 * (th2 - thlast);
        } else {
          not_converged = false;
        }
      }
      // CV:957 (cloud: unused); CV:959-963 pseudoadiabat
      qt = qv2;
      ql2 = 0.0;
      qi2 = 0.0;
    }
    thv2 = th2 * (1.0 + reps * qv2) / (1.0 + qv2 + ql2 + qi2);
    b2 = g * (thv2 - thv[k]) / thv[k];
    const real dz = -cpdg * 0.5 * (thv[k] + thv[k - 1]) * (pi[k] - pi[k - 1]);
    // CV:983-1018
    real parea;
    if ((b2 >= 0.0) && (b1 < 0.0)) {
      const real frac = b2 / (b2 - b1);
      parea = 0.5 * b2 * dz * frac;
      narea = narea - 0.5 * b1 * dz * (1.0 - frac);
      cin = cin + narea;
      narea = 0.0;
    } else if ((b2 < 0.0) && (b1 > 0.0)) {
      const real frac = b1 / (b1 - b2);
      parea = 0.5 * b1 * dz * frac;
      narea = -0.5 * b2 * dz * (1.0 - frac);
    } else if (b2 < 0.0) {
      parea = 0.0;
      narea = narea - 0.5 * dz * (b1 + b2);
    } else {
      parea = 0.5 * dz * (b1 + b2);
      narea = 0.0;
    }
    cape = cape + maxsd(parea, 0.0);
    if ((p[k] <= 10000.0) && (b2 < 0.0)) doit = false;
  }
}

}  // namespace

void convective_diagnostics_setup(DiagnosticsContext& c) {
  // CV:73-82
  c.is_needed_updraft_helicity = c.output_stream_inclusion_count("updraft_helicity_max") > 0;
  c.is_needed_w_max = c.output_stream_inclusion_count("w_velocity_max") > 0;
  c.is_needed_lml_wsp_max = c.output_stream_inclusion_count("wind_speed_level1_max") > 0;
}

// CV:99-198
void convective_diagnostics_update(DiagnosticsContext& c) {
  if (!(c.is_needed_updraft_helicity || c.is_needed_w_max || c.is_needed_lml_wsp_max)) return;
  const Mesh& m = *c.mesh;
  const State& s = *c.state;
  const Diag& d = *c.diag;
  const int nCells = m.nCells, nVertLevels = m.nVertLevels, nCellsSolve = m.nCellsSolve();
  const int nVertices = m.nVertices, nVertexDegree = m.vertexDegree;
  const int nVertLevelsP1 = nVertLevels + 1;
  auto& updraft_helicity_max = c.pools->R("updraft_helicity_max");
  auto& w_velocity_max = c.pools->R("w_velocity_max");
  auto& wind_speed_level1_max = c.pools->R("wind_speed_level1_max");
  const auto& vorticity = d.vorticity;
  const auto& uzonal = d.uReconstructZonal;
  const auto& umeridional = d.uReconstructMeridional;
  const auto& w = s.w;
  std::vector<real> z_agl(nVertLevelsP1);

  if (c.is_needed_w_max) {  // CV:156-160
    for (int iCell = 1; iCell <= nCellsSolve; ++iCell) {
      real mx = w(1, iCell);  // maxval(w(1:nVertLevels,iCell))
      for (int k = 2; k <= nVertLevels; ++k) mx = maxsd(w(k, iCell), mx);
      w_velocity_max(iCell) = maxsd(mx, w_velocity_max(iCell));
    }
  }
  if (c.is_needed_lml_wsp_max) {  // CV:162-166
    for (int iCell = 1; iCell <= nCellsSolve; ++iCell)
      wind_speed_level1_max(iCell) =
          maxsd(std::sqrt(uzonal(1, iCell) * uzonal(1, iCell) + umeridional(1, iCell) * umeridional(1, iCell)),
                wind_speed_level1_max(iCell));
  }
  if (c.is_needed_updraft_helicity) {  // CV:168-193
    FArray<real> updraft_helicity({nVertLevels, nCells + 1}, 0.);
    for (int i = 1; i <= nVertices; ++i)
      for (int j = 1; j <= nVertexDegree; ++j) {
        const int iCell = m.cellsOnVertex(j, i);
        for (int k = 1; k <= nVertLevels; ++k)
          updraft_helicity(k, iCell) = updraft_helicity(k, iCell) + m.kiteAreasOnVertex(j, i) * vorticity(k, i);
      }
    for (int iCell = 1; iCell <= nCellsSolve; ++iCell)
      for (int k = 1; k <= nVertLevels; ++k)
        updraft_helicity(k, iCell) = maxsd(0.5 * (w(k, iCell) + w(k + 1, iCell)), 0.) *
                                     maxsd(updraft_helicity(k, iCell) / m.areaCell(iCell), 0.);
    for (int iCell = 1; iCell <= nCellsSolve; ++iCell) {
      for (int k = 1; k <= nVertLevelsP1; ++k) z_agl[k - 1] = m.zgrid(k, iCell) - m.zgrid(1, iCell);
      const real uph = integral_zstaggered(&updraft_helicity(1, iCell), z_agl.data(), 2000.0, 5000.0, nVertLevels);
      updraft_helicity_max(iCell) = maxsd(uph, updraft_helicity_max(iCell));
    }
  }
}

// CV:227-487
void convective_diagnostics_compute(DiagnosticsContext& c) {
  using physics_constants::rvord;
  bool need_any_diags = false;
  auto need = [&](const char* f) {
    const bool b = c.field_will_be_written(f);
    need_any_diags = need_any_diags || b;
    return b;
  };
  const bool need_cape = need("cape");
  const bool need_cin = need("cin");
  (void)need("lcl");  // CV:289-292: lcl and lfc are asked for but never computed (they keep their values)
  (void)need("lfc");
  const bool need_srh_01km = need("srh_0_1km");
  const bool need_srh_03km = need("srh_0_3km");
  (void)need("uzonal_surface");  // CV:297-312: the surface values are set whenever any diagnostic is needed
  const bool need_uzonal_1km = need("uzonal_1km");
  const bool need_uzonal_6km = need("uzonal_6km");
  (void)need("umeridional_surface");
  const bool need_umerid_1km = need("umeridional_1km");
  const bool need_umerid_6km = need("umeridional_6km");
  (void)need("temperature_surface");
  (void)need("dewpoint_surface");
  if (!need_any_diags) return;

  const Mesh& m = *c.mesh;
  const State& s = *c.state;
  const Diag& d = *c.diag;
  RegistryPools& p = *c.pools;
  const int nCells = m.nCells, nCellsSolve = m.nCellsSolve(), nVertLevels = m.nVertLevels;
  const int nVertLevelsP1 = nVertLevels + 1;
  auto& cape = p.R("cape");
  auto& cin = p.R("cin");
  (void)p.R("lcl");
  (void)p.R("lfc");
  auto& srh_0_1km = p.R("srh_0_1km");
  auto& srh_0_3km = p.R("srh_0_3km");
  auto& uzonal_surface = p.R("uzonal_surface");
  auto& uzonal_1km = p.R("uzonal_1km");
  auto& uzonal_6km = p.R("uzonal_6km");
  auto& umeridional_surface = p.R("umeridional_surface");
  auto& umeridional_1km = p.R("umeridional_1km");
  auto& umeridional_6km = p.R("umeridional_6km");
  auto& temperature_surface = p.R("temperature_surface");
  auto& dewpoint_surface = p.R("dewpoint_surface");
  const auto& height = m.zgrid;
  const auto& umeridional = d.uReconstructMeridional;
  const auto& uzonal = d.uReconstructZonal;
  const auto& theta_m = s.theta_m;
  const auto& scalars = s.scalars;
  const auto& exner = d.exner;
  const auto& pressure_base = d.pressure_base;
  const auto& pressure_p = d.pressure_p;
  int index_qv = 0;
  for (std::size_t k = 0; k < s.scalar_names.size(); ++k)
    if (s.scalar_names[k] == "qv") index_qv = static_cast<int>(k) + 1;
  if (index_qv == 0) throw Error("convective_diagnostics_compute: qv is not an active scalar");

  constexpr real dev_motion = 7.5, z_bunker_bot = 0., z_bunker_top = 6000.;
  FArray<real> temperature({nVertLevels, nCells}), dewpoint({nVertLevels, nCells});
  std::vector<real> dudz(nVertLevels + 1), dvdz(nVertLevels + 1), zp(nVertLevels), zrel(nVertLevels + 1),
      srh(nVertLevels + 1), p_in(nVertLevels), t_in(nVertLevels), td_in(nVertLevels), ucol(nVertLevels),
      vcol(nVertLevels);

  // CV:363-374: k in pairs through _ZGVbN2v_log, an odd last level through log@plt
  const int nvec = (nVertLevels / 2) * 2;
  for (int iCell = 1; iCell <= nCells; ++iCell)
    for (int k = 1; k <= nVertLevels; ++k) {
      temperature(k, iCell) =
          (theta_m(k, iCell) / (1. + rvord * scalars(index_qv, k, iCell))) * exner(k, iCell);
      real evp = 0.01 * (pressure_base(k, iCell) + pressure_p(k, iCell)) * scalars(index_qv, k, iCell) /
                 (scalars(index_qv, k, iCell) + 0.622);
      evp = maxsd(evp, 1.0e-8);
      const real lg = k <= nvec ? libm::mvec_log(evp / 6.112) : libm::log(evp / 6.112);
      dewpoint(k, iCell) = (243.5 * lg) / (17.67 - lg);
      dewpoint(k, iCell) = dewpoint(k, iCell) + 273.15;
    }

  // CV:378-448
  for (int iCell = 1; iCell <= nCellsSolve; ++iCell) {
    for (int k = 1; k <= nVertLevels; ++k) {
      zp[k - 1] = 0.5 * (height(k, iCell) + height(k + 1, iCell)) - height(1, iCell);
      ucol[k - 1] = uzonal(k, iCell);
      vcol[k - 1] = umeridional(k, iCell);
    }
    for (int k = 1; k <= nVertLevelsP1; ++k) zrel[k - 1] = height(k, iCell) - height(1, iCell);

    uzonal_surface(iCell) = uzonal(1, iCell);
    umeridional_surface(iCell) = umeridional(1, iCell);
    temperature_surface(iCell) = temperature(1, iCell);
    dewpoint_surface(iCell) = dewpoint(1, iCell);
    if (need_uzonal_1km) uzonal_1km(iCell) = column_height_value(ucol.data(), zp.data(), 1000.0, nVertLevels);
    if (need_umerid_1km) umeridional_1km(iCell) = column_height_value(vcol.data(), zp.data(), 1000.0, nVertLevels);
    if (need_uzonal_6km) uzonal_6km(iCell) = column_height_value(ucol.data(), zp.data(), 6000.0, nVertLevels);
    if (need_umerid_6km) umeridional_6km(iCell) = column_height_value(vcol.data(), zp.data(), 6000.0, nVertLevels);

    if (need_srh_01km || need_srh_03km) {
      real u_srh_bot = uzonal(1, iCell);
      real v_srh_bot = umeridional(1, iCell);
      if (z_bunker_bot > zp[0]) {
        u_srh_bot = column_height_value(ucol.data(), zp.data(), z_bunker_bot, nVertLevels);
        v_srh_bot = column_height_value(vcol.data(), zp.data(), z_bunker_bot, nVertLevels);
      }
      const real u_srh_top = column_height_value(ucol.data(), zp.data(), z_bunker_top, nVertLevels);
      const real v_srh_top = column_height_value(vcol.data(), zp.data(), z_bunker_top, nVertLevels);
      const real u_shear = u_srh_top - u_srh_bot;
      const real v_shear = v_srh_top - v_srh_bot;
      const real u_mean = integral_zstaggered(ucol.data(), zrel.data(), z_bunker_bot, z_bunker_top, nVertLevels) /
                          (z_bunker_top - z_bunker_bot);
      const real v_mean = integral_zstaggered(vcol.data(), zrel.data(), z_bunker_bot, z_bunker_top, nVertLevels) /
                          (z_bunker_top - z_bunker_bot);
      const real shear_magnitude = maxsd(std::sqrt(u_shear * u_shear + v_shear * v_shear), 0.0001);
      const real u_storm = u_mean + dev_motion * v_shear / shear_magnitude;
      const real v_storm = v_mean - dev_motion * u_shear / shear_magnitude;

      for (int k = 2; k <= nVertLevels - 1; ++k) {
        dudz[k] = (uzonal(k, iCell) - uzonal(k - 1, iCell)) / (0.5 * (height(k + 1, iCell) - height(k - 1, iCell)));
        dvdz[k] = (umeridional(k, iCell) - umeridional(k - 1, iCell)) /
                  (0.5 * (height(k + 1, iCell) - height(k - 1, iCell)));
      }
      dudz[1] = dudz[2];
      dvdz[1] = dvdz[2];
      dudz[nVertLevels] = dudz[nVertLevels - 1];
      dvdz[nVertLevels] = dvdz[nVertLevels - 1];

      for (int k = 2; k <= nVertLevels; ++k)
        srh[k - 1] = -((0.5 * (uzonal(k, iCell) + uzonal(k - 1, iCell)) - u_storm) * dvdz[k]) +
                     (0.5 * (umeridional(k, iCell) + umeridional(k - 1, iCell)) - v_storm) * dudz[k];
      srh[0] = -((uzonal(1, iCell) - u_storm) * dvdz[1]) + (umeridional(1, iCell) - v_storm) * dudz[1];
      srh[nVertLevelsP1 - 1] = srh[nVertLevels - 1];
      for (int k = 1; k <= nVertLevels + 1; ++k) srh[k - 1] = maxsd(srh[k - 1], 0.);

      if (need_srh_01km) srh_0_1km(iCell) = integral_zpoint(srh.data(), zrel.data(), 0.0, 1000.0, nVertLevelsP1);
      if (need_srh_03km) srh_0_3km(iCell) = integral_zpoint(srh.data(), zrel.data(), 0.0, 3000.0, nVertLevelsP1);
    }
  }

  // CV:451-474
  if (need_cape || need_cin) {
    for (int iCell = 1; iCell <= nCellsSolve; ++iCell) {
      for (int k = 1; k <= nVertLevels; ++k) {
        p_in[k - 1] = (pressure_p(k, iCell) + pressure_base(k, iCell)) / 100.0;
        t_in[k - 1] = temperature(k, iCell) - 273.15;
        td_in[k - 1] = dewpoint(k, iCell) - 273.15;
      }
      real cape_out = 0., cin_out = 0.;
      getcape(nVertLevels, p_in.data(), t_in.data(), td_in.data(), cape_out, cin_out);
      cape(iCell) = cape_out;
      cin(iCell) = cin_out;
    }
  }
}

void convective_diagnostics_reset(DiagnosticsContext& c) {
  // CV:515-523: the whole arrays (the garbage element included)
  for (const char* f : {"updraft_helicity_max", "w_velocity_max", "wind_speed_level1_max"})
    if (c.field_will_be_written(f)) {
      auto& a = c.pools->R(f);
      for (std::size_t n = 0; n < a.size(); ++n) a.data()[n] = 0.;
    }
}

}  // namespace cppmpas::diag
