// SPDX-License-Identifier: GPL-3.0-or-later
// See mpas_atmphys_microphysics.hpp. PI = src/core_atmosphere/physics/mpas_atmphys_interface.F,
// DRV = physics/mpas_atmphys_driver_microphysics.F (MPAS v8.4.2). Statement by statement.
//
// Reference object code (bin/atmosphere_model, sha256 524a1cf4...): microphysics_from_mpas and
// microphysics_to_mpas have no FMA, no ymm register and no libmvec call; the only libm call is the scalar pow@plt
// of the exner function (PI:926). precip_from_mpas / precip_to_mpas call no libm. max(0._RKIND, qv) is a maxsd with
// qv as the destination ("qv > 0 ? qv : 0").
#include "physics/mpas_atmphys_microphysics.hpp"

#include <cstddef>

#include "core/backend.hpp"
#include "physics/mpas_atmphys_interface.hpp"

#include "core/error.hpp"
#include "core/libm.hpp"
#include "physics/module_mp_wsm6.hpp"
#include "physics/mpas_atmphys_constants.hpp"

namespace cppmpas::phys {

namespace {
CPPMPAS_HD inline real fmax0(real x) { return x > 0.0 ? x : 0.0; }  // max(0._RKIND, x)
}  // namespace

CPPMPAS_HD real rslf(real p, real t) {  // mpas_atmphys_functions.F:156-186
  constexpr real C0 = .611583699E03, C1 = .444606896E02, C2 = .143177157E01, C3 = .264224321E-1,
                 C4 = .299291081E-3, C5 = .203154182E-5, C6 = .702620698E-8, C7 = .379534310E-11,
                 C8 = -.321582393E-13;
  const real xt = t - 273.16;
  const real X = xt > -80. ? xt : -80.;
  real ESL = C0 + X * (C1 + X * (C2 + X * (C3 + X * (C4 + X * (C5 + X * (C6 + X * (C7 + X * C8)))))));
  const real pc = p * 0.15;
  ESL = ESL < pc ? ESL : pc;
  return .622 * ESL / (p - ESL);
}

CPPMPAS_HD real rsif(real p, real t) {  // mpas_atmphys_functions.F:191-217
  constexpr real C0 = .609868993E03, C1 = .499320233E02, C2 = .184672631E01, C3 = .402737184E-1,
                 C4 = .565392987E-3, C5 = .521693933E-5, C6 = .307839583E-7, C7 = .105785160E-9,
                 C8 = .161444444E-12;
  const real xt = t - 273.16;
  const real X = xt > -80. ? xt : -80.;
  real ESI = C0 + X * (C1 + X * (C2 + X * (C3 + X * (C4 + X * (C5 + X * (C6 + X * (C7 + X * C8)))))));
  const real pc = p * 0.15;
  ESI = ESI < pc ? ESI : pc;
  return .622 * ESI / (p - ESI);
}

// One kernel per loop of DRV / PI (one thread per cell i, the k loops inside the thread in MPAS's order).
void driver_microphysics_wsm6(const wsm6::Wsm6Params<real>& params, MicrophysicsPools& mp, int its, int ite) {
  using namespace physics_constants;
  if (mp.l_diags && !mp.radar) throw Error("driver_microphysics: l_diags without the WSM6 radar tables");
  if (mp.l_acrain) throw Error("driver_microphysics: l_acrain (config_bucket_update) has no verified port");
  if (its != 1) throw Error("driver_microphysics: its /= 1");
  if (mp.index_qv <= 0 || mp.index_qc <= 0 || mp.index_qr <= 0 || mp.index_qi <= 0 || mp.index_qs <= 0 ||
      mp.index_qg <= 0)  // an inactive constituent has MPAS's index -1 (scalar_index)
    throw Error("driver_microphysics: WSM6 needs qv, qc, qr, qi, qs and qg");
  const int kts = 1, kte = mp.nVertLevels;
  if (mp.mp_top_level < 1 || mp.mp_top_level > kte) throw Error("driver_microphysics: mp_top_level outside 1..nVertLevels");
  // the surface pressure of PI:950-964 reads atmospheric level 2 and interface 3: a 1-level column is past it
  if (kte < 2) throw Error("driver_microphysics: nVertLevels < 2 (the surface pressure of PI:950-964)");
  if (mp.l_diags && kte > wsm6::kReflMaxLevels) throw Error("refl10cm_wsm6: nVertLevels outside 1..128");
  KernelCopy<wsm6::Wsm6Radar> radar_k;
  MicrophysicsPools m = mp;
  if (m.l_diags) m.radar = kernel_ptr(*mp.radar, radar_k);

  // allocate_microphysics (DRV:98-176): (ims:ime, kms:kme, jms:jme), ims = its = 1, ime = ite, kme = kte + 1
  const int ldi = ite;
  const std::size_t n3 = static_cast<std::size_t>(ldi) * (kte + 1);
  WorkArray rho_a(n3), th_a(n3), pi_a(n3), pres_a(n3), z_a(n3), dz_a(n3), w_a(n3), qv_a(n3), qc_a(n3), qr_a(n3),
      qi_a(n3), qs_a(n3), qg_a(n3), rainprod_a(n3), evapprod_a(n3), recloud_a(n3), reice_a(n3), resnow_a(n3);
  WorkArray rainnc_a(ldi), rainncv_a(ldi), sr_a(ldi), snownc_a(ldi), snowncv_a(ldi), graupelnc_a(ldi),
      graupelncv_a(ldi);
  real *rho_p = rho_a.data(), *th_p = th_a.data(), *pi_p = pi_a.data(), *pres_p = pres_a.data(), *z_p = z_a.data(),
       *dz_p = dz_a.data(), *w_p = w_a.data(), *qv_p = qv_a.data(), *qc_p = qc_a.data(), *qr_p = qr_a.data(),
       *qi_p = qi_a.data(), *qs_p = qs_a.data(), *qg_p = qg_a.data(), *rainprod_p = rainprod_a.data(),
       *evapprod_p = evapprod_a.data(), *recloud_p = recloud_a.data(), *reice_p = reice_a.data(),
       *resnow_p = resnow_a.data();
  real *rainnc_p = rainnc_a.data(), *rainncv_p = rainncv_a.data(), *sr_p = sr_a.data(), *snownc_p = snownc_a.data(),
       *snowncv_p = snowncv_a.data(), *graupelnc_p = graupelnc_a.data(), *graupelncv_p = graupelncv_a.data();
#define P(a, i, k) (a)[static_cast<std::size_t>((k) - 1) * ldi + ((i) - 1)]
#define qv(k, i) m.scalars(m.index_qv, k, i)
#define qc(k, i) m.scalars(m.index_qc, k, i)
#define qr(k, i) m.scalars(m.index_qr, k, i)
#define qi(k, i) m.scalars(m.index_qi, k, i)
#define qs(k, i) m.scalars(m.index_qs, k, i)
#define qg(k, i) m.scalars(m.index_qg, k, i)

  // precip_from_MPAS (DRV:530-562)
  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {
    rainncv_p[i - 1] = 0.;
    rainnc_p[i - 1] = 0.;
  });
  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) { m.rainncv(i) = 0.; });
  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {  // mp_wsm6
    snowncv_p[i - 1] = 0.;
    graupelncv_p[i - 1] = 0.;
    snownc_p[i - 1] = 0.;
    graupelnc_p[i - 1] = 0.;
    sr_p[i - 1] = 0.;
  });
  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {
    m.snowncv(i) = 0.;
    m.graupelncv(i) = 0.;
    m.sr(i) = 0.;
  });

  // microphysics_from_MPAS (PI:598-823)
  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {  // PI:656-675
    for (int k = kts; k <= kte; ++k) {
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
  });
  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {  // PI:693-708
    for (int k = kts; k <= kte; ++k) {
      P(qi_p, i, k) = qi(k, i);
      P(qs_p, i, k) = qs(k, i);
      P(qg_p, i, k) = qg(k, i);
      P(rainprod_p, i, k) = m.rainprod(k, i);
      P(evapprod_p, i, k) = m.evapprod(k, i);
      P(recloud_p, i, k) = m.re_cloud(k, i);
      P(reice_p, i, k) = m.re_ice(k, i);
      P(resnow_p, i, k) = m.re_snow(k, i);
    }
  });
  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {  // PI:776-786
    for (int k = kts; k <= kte; ++k) {
      m.rthmpten(k, i) = m.theta_m(k, i) / (1. + R_v / R_d * fmax0(qv(k, i)));
      m.rqvmpten(k, i) = qv(k, i);
      m.rqcmpten(k, i) = qc(k, i);
      m.rqrmpten(k, i) = qr(k, i);
      m.rqimpten(k, i) = qi(k, i);
      m.rqsmpten(k, i) = qs(k, i);
      m.rqgmpten(k, i) = qg(k, i);
    }
  });

  // DRV:437-461: wsm6 (rainprod_p / evapprod_p are untouched without WRF_CHEM; re_* are untouched with
  // has_reqc = 0, module_mp_wsm6.hpp)
  {
    module_mp_wsm6::Wsm6Columns a;
    a.ldi = ldi;
    a.th = th_p;
    a.q = qv_p;
    a.qc = qc_p;
    a.qr = qr_p;
    a.qi = qi_p;
    a.qs = qs_p;
    a.qg = qg_p;
    a.den = rho_p;
    a.pii = pi_p;
    a.p = pres_p;
    a.delz = dz_p;
    a.rain = rainnc_p;
    a.rainncv = rainncv_p;
    a.snow = snownc_p;
    a.snowncv = snowncv_p;
    a.graupel = graupelnc_p;
    a.graupelncv = graupelncv_p;
    a.sr = sr_p;
    module_mp_wsm6::wsm6(params, m.dt_microp, m.has_reqc, m.has_reqi, m.has_reqs, its, ite, kts, kte, a);
  }

  if (m.l_diags) {  // DRV:467-482: mp_wsm6 computes the radar reflectivity, then the relative humidity
    // compute_radar_reflectivity, mp_wsm6 (DRV:690-730); the column arrays are the thread's own
    parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {
      constexpr int KR = wsm6::kReflMaxLevels;
      double p1d[KR], t1d[KR], qv1d[KR], qr1d[KR], qs1d[KR], qg1d[KR], dBZ1d[KR], zp[KR];
      for (int k = kts; k <= kte; ++k) {
        p1d[k - 1] = P(pres_p, i, k);
        t1d[k - 1] = P(th_p, i, k) * P(pi_p, i, k);
        qv1d[k - 1] = P(qv_p, i, k);
        qr1d[k - 1] = P(qr_p, i, k);
        qs1d[k - 1] = P(qs_p, i, k);
        qg1d[k - 1] = P(qg_p, i, k);
        dBZ1d[k - 1] = -35.;
        zp[k - 1] = P(z_p, i, k) - P(z_p, i, 1) + 0.5 * P(dz_p, i, k);  // height AGL
      }
      wsm6::refl10cm_wsm6(*m.radar, qv1d, qr1d, qs1d, qg1d, t1d, p1d, dBZ1d, kte);
      int kp = 1;
      for (int k = kts; k <= kte; ++k) {
        dBZ1d[k - 1] = -35. > dBZ1d[k - 1] ? -35. : dBZ1d[k - 1];  // max(-35._RKIND, dBZ1d(k))
        m.refl10cm(k, i) = dBZ1d[k - 1];
        if (zp[k - 1] < 1000.) kp = k;
      }
      real mx = dBZ1d[0];  // maxval(dBZ1d(:))
      for (int k = kts + 1; k <= kte; ++k)
        if (dBZ1d[k - 1] > mx) mx = dBZ1d[k - 1];
      m.refl10cm_max(i) = mx;
      // DRV:722-724 read zp(kp+1) and dBZ1d(kp+1): with no midpoint at or above 1 km AGL (kp = kte, also a 1-level
      // column) that is past the arrays MPAS allocated (kts:kte), undefined in MPAS: refused
      if (kp >= kte) kernel_fail("refl10cm_1km: no level midpoint at or above 1 km AGL (MPAS reads zp(kte+1))");
      const real w1 = (zp[kp] - 1000.) / (zp[kp] - zp[kp - 1]);
      const real w2 = 1.0 - w1;
      m.refl10cm_1km(i) = w1 * dBZ1d[kp - 1] + w2 * dBZ1d[kp];
      // max(refl10cm_1km_max(i), refl10cm_1km(i)): maxsd keeps the first argument unless the second is greater (a
      // tie of -0.0 / +0.0 and a NaN keep refl10cm_1km_max)
      m.refl10cm_1km_max(i) = m.refl10cm_1km(i) > m.refl10cm_1km_max(i) ? m.refl10cm_1km(i) : m.refl10cm_1km_max(i);
    });
    // compute_relhum (DRV:798-850)
    parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {
      for (int k = kts; k <= kte; ++k) m.relhum(k, i) = 0.;
    });
    parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {
      for (int k = kts; k <= kte; ++k) {
        const real p1 = P(pres_p, i, k);
        const real t1 = P(th_p, i, k) * P(pi_p, i, k);
        const real tempc = t1 - 273.16;
        real qvs1 = rslf(p1, t1);
        if (tempc <= 0.) qvs1 = rsif(p1, t1);
        const real qv1 = P(qv_p, i, k);
        m.relhum(k, i) = qv1 / qvs1;
        m.relhum(k, i) = m.relhum(k, i) * 100.;
      }
    });
  }

  // precip_to_MPAS (DRV:567-653), l_acrain = .false.
  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) { m.precipw(i) = 0.; });
  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {
    for (int k = kts; k <= kte; ++k) {
      const real rho_a = P(rho_p, i, k) / (1. + P(qv_p, i, k));
      m.precipw(i) = m.precipw(i) + P(qv_p, i, k) * rho_a * P(dz_p, i, k);
    }
    m.rainncv(i) = rainnc_p[i - 1];
    m.rainnc(i) = m.rainnc(i) + m.rainncv(i);
  });
  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {  // DRV:636-648
    m.snowncv(i) = snownc_p[i - 1];
    m.graupelncv(i) = graupelnc_p[i - 1];
    m.sr(i) = (snownc_p[i - 1] + graupelnc_p[i - 1]) / (rainnc_p[i - 1] + 1.e-12);
    m.snownc(i) = m.snownc(i) + m.snowncv(i);
    m.graupelnc(i) = m.graupelnc(i) + m.graupelncv(i);
  });

  // microphysics_to_MPAS (PI:826-1138)
  const int mp_top_level = m.mp_top_level;
  const real dt_dyn = m.dt_dyn;
  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {  // PI:903-934
    for (int k = kts; k <= mp_top_level; ++k) {
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
  });
  if (kte > mp_top_level)  // PI:938-946
    parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {
      for (int k = mp_top_level + 1; k <= kte; ++k) {
        qc(k, i) = 0.;
        qr(k, i) = 0.;
        m.rt_diabatic_tend(k, i) = 0.;
      }
    });
  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {  // PI:950-964
    const real tem1 = m.zgrid(2, i) - m.zgrid(1, i);
    const real tem2 = m.zgrid(3, i) - m.zgrid(2, i);
    const real rho1 = m.rho_zz(1, i) * m.zz(1, i) * (1. + P(qv_p, i, 1));
    const real rho2 = m.rho_zz(2, i) * m.zz(2, i) * (1. + P(qv_p, i, 2));
    m.tend_sfc_pressure(i) = m.surface_pressure(i);
    m.surface_pressure(i) = 0.5 * gravity * (m.zgrid(2, i) - m.zgrid(1, i)) * (rho1 - 0.5 * (rho2 - rho1) * tem1 / (tem1 + tem2));
    m.surface_pressure(i) = m.surface_pressure(i) + m.pressure_p(1, i) + m.pressure_base(1, i);
    m.tend_sfc_pressure(i) = (m.surface_pressure(i) - m.tend_sfc_pressure(i)) / dt_dyn;
  });
  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {  // PI:983-997
    for (int k = kts; k <= mp_top_level; ++k) {
      qi(k, i) = P(qi_p, i, k);
      qs(k, i) = P(qs_p, i, k);
      qg(k, i) = P(qg_p, i, k);
      m.rainprod(k, i) = P(rainprod_p, i, k);
      m.evapprod(k, i) = P(evapprod_p, i, k);
      m.re_cloud(k, i) = P(recloud_p, i, k);
      m.re_ice(k, i) = P(reice_p, i, k);
      m.re_snow(k, i) = P(resnow_p, i, k);
    }
  });
  if (kte > mp_top_level)  // PI:1000-1014
    parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {
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
    });
  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {  // PI:1076-1086
    for (int k = kts; k <= mp_top_level; ++k) {
      m.rthmpten(k, i) = (m.theta_m(k, i) / (1. + R_v / R_d * fmax0(qv(k, i))) - m.rthmpten(k, i)) / dt_dyn;
      m.rqvmpten(k, i) = (qv(k, i) - m.rqvmpten(k, i)) / dt_dyn;
      m.rqcmpten(k, i) = (qc(k, i) - m.rqcmpten(k, i)) / dt_dyn;
      m.rqrmpten(k, i) = (qr(k, i) - m.rqrmpten(k, i)) / dt_dyn;
      m.rqimpten(k, i) = (qi(k, i) - m.rqimpten(k, i)) / dt_dyn;
      m.rqsmpten(k, i) = (qs(k, i) - m.rqsmpten(k, i)) / dt_dyn;
      m.rqgmpten(k, i) = (qg(k, i) - m.rqgmpten(k, i)) / dt_dyn;
    }
  });
  if (kte > mp_top_level)  // PI:1089-1101
    parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {
      for (int k = mp_top_level + 1; k <= kte; ++k) {
        m.rthmpten(k, i) = 0.;
        m.rqvmpten(k, i) = 0.;
        m.rqcmpten(k, i) = 0.;
        m.rqrmpten(k, i) = 0.;
        m.rqimpten(k, i) = 0.;
        m.rqsmpten(k, i) = 0.;
        m.rqgmpten(k, i) = 0.;
      }
    });
#undef P
#undef qv
#undef qc
#undef qr
#undef qi
#undef qs
#undef qg
}

}  // namespace cppmpas::phys
