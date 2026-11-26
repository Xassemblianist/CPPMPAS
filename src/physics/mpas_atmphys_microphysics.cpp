// SPDX-License-Identifier: GPL-3.0-or-later
// See mpas_atmphys_microphysics.hpp. PI = src/core_atmosphere/physics/mpas_atmphys_interface.F,
// DRV = physics/mpas_atmphys_driver_microphysics.F (MPAS v8.4.2). Statement by statement.
//
// Reference object code (bin/atmosphere_model, sha256 524a1cf4...): microphysics_from_mpas and
// microphysics_to_mpas have no FMA, no ymm register and no libmvec call; the only libm call is the scalar pow@plt
// of the exner function (PI:926). precip_from_mpas / precip_to_mpas call no libm. max(0._RKIND, qv) is a maxsd with
// qv as the destination ("qv > 0 ? qv : 0").
#include "physics/mpas_atmphys_microphysics.hpp"

#include <vector>

#include "core/error.hpp"
#include "core/libm.hpp"
#include "physics/module_mp_wsm6.hpp"
#include "physics/mpas_atmphys_constants.hpp"

namespace cppmpas::phys {

namespace {
inline real fmax0(real x) { return x > 0.0 ? x : 0.0; }  // max(0._RKIND, x)
}  // namespace

void driver_microphysics_wsm6(const wsm6::Wsm6Params<real>& params, MicrophysicsPools& m, int its, int ite) {
  using namespace physics_constants;
  if (m.l_diags)  // DRV:466-482
    throw Error("driver_microphysics: l_diags (compute_radar_reflectivity, compute_relhum) has no verified port");
  if (m.l_acrain) throw Error("driver_microphysics: l_acrain (config_bucket_update) has no verified port");
  if (its != 1) throw Error("driver_microphysics: its /= 1");
  if (!m.index_qv || !m.index_qc || !m.index_qr || !m.index_qi || !m.index_qs || !m.index_qg)
    throw Error("driver_microphysics: WSM6 needs qv, qc, qr, qi, qs and qg");
  const int kts = 1, kte = m.nVertLevels;
  if (m.mp_top_level < 1 || m.mp_top_level > kte) throw Error("driver_microphysics: mp_top_level outside 1..nVertLevels");

  // allocate_microphysics (DRV:98-176): (ims:ime, kms:kme, jms:jme), ims = its = 1, ime = ite, kme = kte + 1
  const int ldi = ite;
  const std::size_t n3 = static_cast<std::size_t>(ldi) * (kte + 1);
  std::vector<real> rho_p(n3), th_p(n3), pi_p(n3), pres_p(n3), z_p(n3), dz_p(n3), w_p(n3), qv_p(n3), qc_p(n3),
      qr_p(n3), qi_p(n3), qs_p(n3), qg_p(n3), rainprod_p(n3), evapprod_p(n3), recloud_p(n3), reice_p(n3),
      resnow_p(n3);
  std::vector<real> rainnc_p(ldi), rainncv_p(ldi), sr_p(ldi), snownc_p(ldi), snowncv_p(ldi), graupelnc_p(ldi),
      graupelncv_p(ldi);
  auto P = [&](std::vector<real>& a, int i, int k) -> real& { return a[static_cast<std::size_t>(k - 1) * ldi + (i - 1)]; };
  const auto& qv = [&](int k, int i) -> real& { return m.scalars(m.index_qv, k, i); };
  const auto& qc = [&](int k, int i) -> real& { return m.scalars(m.index_qc, k, i); };
  const auto& qr = [&](int k, int i) -> real& { return m.scalars(m.index_qr, k, i); };
  const auto& qi = [&](int k, int i) -> real& { return m.scalars(m.index_qi, k, i); };
  const auto& qs = [&](int k, int i) -> real& { return m.scalars(m.index_qs, k, i); };
  const auto& qg = [&](int k, int i) -> real& { return m.scalars(m.index_qg, k, i); };

  // precip_from_MPAS (DRV:530-562)
  for (int i = its; i <= ite; ++i) {
    rainncv_p[i - 1] = 0.;
    rainnc_p[i - 1] = 0.;
  }
  for (int i = its; i <= ite; ++i) m.rainncv(i) = 0.;
  for (int i = its; i <= ite; ++i) {  // mp_wsm6
    snowncv_p[i - 1] = 0.;
    graupelncv_p[i - 1] = 0.;
    snownc_p[i - 1] = 0.;
    graupelnc_p[i - 1] = 0.;
    sr_p[i - 1] = 0.;
  }
  for (int i = its; i <= ite; ++i) {
    m.snowncv(i) = 0.;
    m.graupelncv(i) = 0.;
    m.sr(i) = 0.;
  }

  // microphysics_from_MPAS (PI:598-823)
  for (int k = kts; k <= kte; ++k)  // PI:656-675
    for (int i = its; i <= ite; ++i) {
      P(qv_p, i, k) = qv(k, i);
      P(qc_p, i, k) = qc(k, i);
      P(qr_p, i, k) = qr(k, i);
      P(rho_p, i, k) = m.zz(k, i) * m.rho_zz(k, i);
      P(th_p, i, k) = m.theta_m(k, i) / (1. + R_v / R_d * fmax0(P(qv_p, i, k)));
      P(pi_p, i, k) = m.exner(k, i);
      P(pres_p, i, k) = m.pressure_base(k, i) + m.pressure_p(k, i);
      P(z_p, i, k) = m.zgrid(k, i);
      P(dz_p, i, k) = m.zgrid(k + 1, i) - m.zgrid(k, i);
      P(w_p, i, k) = m.w(k, i);
    }
  for (int k = kts; k <= kte; ++k)  // PI:693-708
    for (int i = its; i <= ite; ++i) {
      P(qi_p, i, k) = qi(k, i);
      P(qs_p, i, k) = qs(k, i);
      P(qg_p, i, k) = qg(k, i);
      P(rainprod_p, i, k) = m.rainprod(k, i);
      P(evapprod_p, i, k) = m.evapprod(k, i);
      P(recloud_p, i, k) = m.re_cloud(k, i);
      P(reice_p, i, k) = m.re_ice(k, i);
      P(resnow_p, i, k) = m.re_snow(k, i);
    }
  for (int k = kts; k <= kte; ++k)  // PI:776-786
    for (int i = its; i <= ite; ++i) {
      m.rthmpten(k, i) = m.theta_m(k, i) / (1. + R_v / R_d * fmax0(qv(k, i)));
      m.rqvmpten(k, i) = qv(k, i);
      m.rqcmpten(k, i) = qc(k, i);
      m.rqrmpten(k, i) = qr(k, i);
      m.rqimpten(k, i) = qi(k, i);
      m.rqsmpten(k, i) = qs(k, i);
      m.rqgmpten(k, i) = qg(k, i);
    }

  // DRV:437-461: wsm6 (rainprod_p / evapprod_p are untouched without WRF_CHEM; re_* are untouched with
  // has_reqc = 0, module_mp_wsm6.hpp)
  {
    module_mp_wsm6::Wsm6Columns a;
    a.ldi = ldi;
    a.th = th_p.data();
    a.q = qv_p.data();
    a.qc = qc_p.data();
    a.qr = qr_p.data();
    a.qi = qi_p.data();
    a.qs = qs_p.data();
    a.qg = qg_p.data();
    a.den = rho_p.data();
    a.pii = pi_p.data();
    a.p = pres_p.data();
    a.delz = dz_p.data();
    a.rain = rainnc_p.data();
    a.rainncv = rainncv_p.data();
    a.snow = snownc_p.data();
    a.snowncv = snowncv_p.data();
    a.graupel = graupelnc_p.data();
    a.graupelncv = graupelncv_p.data();
    a.sr = sr_p.data();
    module_mp_wsm6::wsm6(params, m.dt_microp, m.has_reqc, m.has_reqi, m.has_reqs, its, ite, kts, kte, a);
  }

  // precip_to_MPAS (DRV:567-653), l_acrain = .false.
  for (int i = its; i <= ite; ++i) m.precipw(i) = 0.;
  for (int i = its; i <= ite; ++i) {
    for (int k = kts; k <= kte; ++k) {
      const real rho_a = P(rho_p, i, k) / (1. + P(qv_p, i, k));
      m.precipw(i) = m.precipw(i) + P(qv_p, i, k) * rho_a * P(dz_p, i, k);
    }
    m.rainncv(i) = rainnc_p[i - 1];
    m.rainnc(i) = m.rainnc(i) + m.rainncv(i);
  }
  for (int i = its; i <= ite; ++i) {  // DRV:636-648
    m.snowncv(i) = snownc_p[i - 1];
    m.graupelncv(i) = graupelnc_p[i - 1];
    m.sr(i) = (snownc_p[i - 1] + graupelnc_p[i - 1]) / (rainnc_p[i - 1] + 1.e-12);
    m.snownc(i) = m.snownc(i) + m.snowncv(i);
    m.graupelnc(i) = m.graupelnc(i) + m.graupelncv(i);
  }

  // microphysics_to_MPAS (PI:826-1138)
  const int mp_top_level = m.mp_top_level;
  const real dt_dyn = m.dt_dyn;
  for (int k = kts; k <= mp_top_level; ++k)  // PI:903-934
    for (int i = its; i <= ite; ++i) {
      m.rt_diabatic_tend(k, i) = m.theta_m(k, i);
      m.dtheta_dt_mp(k, i) = m.theta_m(k, i) / (1. + rvord * qv(k, i));
      qv(k, i) = P(qv_p, i, k);
      qc(k, i) = P(qc_p, i, k);
      qr(k, i) = P(qr_p, i, k);
      m.theta_m(k, i) = P(th_p, i, k) * (1. + rvord * P(qv_p, i, k));
      m.rt_diabatic_tend(k, i) = (m.theta_m(k, i) - m.rt_diabatic_tend(k, i)) / dt_dyn;
      m.dtheta_dt_mp(k, i) = (m.theta_m(k, i) / (1. + rvord * qv(k, i)) - m.dtheta_dt_mp(k, i)) / (dt_dyn);
      m.rtheta_p(k, i) = m.rho_zz(k, i) * m.theta_m(k, i) - m.rtheta_base(k, i);
      m.exner(k, i) = libm::pow(m.zz(k, i) * (R_d / P0) * (m.rtheta_p(k, i) + m.rtheta_base(k, i)), rcv);
      m.pressure_p(k, i) = m.zz(k, i) * R_d * (m.exner(k, i) * m.rtheta_p(k, i)
                         + (m.exner(k, i) - m.exner_base(k, i)) * m.rtheta_base(k, i));
    }
  if (kte > mp_top_level)  // PI:938-946
    for (int i = its; i <= ite; ++i)
      for (int k = mp_top_level + 1; k <= kte; ++k) {
        qc(k, i) = 0.;
        qr(k, i) = 0.;
        m.rt_diabatic_tend(k, i) = 0.;
      }
  for (int i = its; i <= ite; ++i) {  // PI:950-964
    const real tem1 = m.zgrid(2, i) - m.zgrid(1, i);
    const real tem2 = m.zgrid(3, i) - m.zgrid(2, i);
    const real rho1 = m.rho_zz(1, i) * m.zz(1, i) * (1. + P(qv_p, i, 1));
    const real rho2 = m.rho_zz(2, i) * m.zz(2, i) * (1. + P(qv_p, i, 2));
    m.tend_sfc_pressure(i) = m.surface_pressure(i);
    m.surface_pressure(i) = 0.5 * gravity * (m.zgrid(2, i) - m.zgrid(1, i)) * (rho1 - 0.5 * (rho2 - rho1) * tem1 / (tem1 + tem2));
    m.surface_pressure(i) = m.surface_pressure(i) + m.pressure_p(1, i) + m.pressure_base(1, i);
    m.tend_sfc_pressure(i) = (m.surface_pressure(i) - m.tend_sfc_pressure(i)) / dt_dyn;
  }
  for (int k = kts; k <= mp_top_level; ++k)  // PI:983-997
    for (int i = its; i <= ite; ++i) {
      qi(k, i) = P(qi_p, i, k);
      qs(k, i) = P(qs_p, i, k);
      qg(k, i) = P(qg_p, i, k);
      m.rainprod(k, i) = P(rainprod_p, i, k);
      m.evapprod(k, i) = P(evapprod_p, i, k);
      m.re_cloud(k, i) = P(recloud_p, i, k);
      m.re_ice(k, i) = P(reice_p, i, k);
      m.re_snow(k, i) = P(resnow_p, i, k);
    }
  if (kte > mp_top_level)  // PI:1000-1014
    for (int i = its; i <= ite; ++i)
      for (int k = mp_top_level + 1; k <= kte; ++k) {
        qi(k, i) = 0.;
        qs(k, i) = 0.;
        qg(k, i) = 0.;
        m.rainprod(k, i) = 0.;
        m.evapprod(k, i) = 0.;
        m.re_cloud(k, i) = re_qc_bg;
        m.re_ice(k, i) = re_qi_bg;
        m.re_snow(k, i) = re_qs_bg;
      }
  for (int i = its; i <= ite; ++i)  // PI:1076-1086
    for (int k = kts; k <= mp_top_level; ++k) {
      m.rthmpten(k, i) = (m.theta_m(k, i) / (1. + R_v / R_d * fmax0(qv(k, i))) - m.rthmpten(k, i)) / dt_dyn;
      m.rqvmpten(k, i) = (qv(k, i) - m.rqvmpten(k, i)) / dt_dyn;
      m.rqcmpten(k, i) = (qc(k, i) - m.rqcmpten(k, i)) / dt_dyn;
      m.rqrmpten(k, i) = (qr(k, i) - m.rqrmpten(k, i)) / dt_dyn;
      m.rqimpten(k, i) = (qi(k, i) - m.rqimpten(k, i)) / dt_dyn;
      m.rqsmpten(k, i) = (qs(k, i) - m.rqsmpten(k, i)) / dt_dyn;
      m.rqgmpten(k, i) = (qg(k, i) - m.rqgmpten(k, i)) / dt_dyn;
    }
  if (kte > mp_top_level)  // PI:1089-1101
    for (int i = its; i <= ite; ++i)
      for (int k = mp_top_level + 1; k <= kte; ++k) {
        m.rthmpten(k, i) = 0.;
        m.rqvmpten(k, i) = 0.;
        m.rqcmpten(k, i) = 0.;
        m.rqrmpten(k, i) = 0.;
        m.rqimpten(k, i) = 0.;
        m.rqsmpten(k, i) = 0.;
        m.rqgmpten(k, i) = 0.;
      }
}

}  // namespace cppmpas::phys
