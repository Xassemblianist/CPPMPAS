// SPDX-License-Identifier: GPL-3.0-or-later
// driver_microphysics for config_microp_scheme = 'mp_thompson' (see mpas_atmphys_microphysics.hpp). MPAS v8.4.2:
// DRV = physics/mpas_atmphys_driver_microphysics.F, PI = physics/mpas_atmphys_interface.F, MT =
// physics/physics_wrf/module_mp_thompson.F, IM = physics/mpas_atmphys_init_microphysics.F. Statement by statement.
//
// Every loop over the cells is a kernel (core/backend.hpp parallel_for, one thread per cell i with the k loops of MPAS
// inside it, in MPAS's order; the (i, k) statements of these loops are independent, so the loop order MPAS writes,
// k outside i, gives the same values). The column scheme (physics/mp_thompson.hpp) is CPPMPAS_HD. The driver's
// work arrays are WorkArrays and the ThompsonView a KernelCopy (managed memory in the CUDA build), and the tables it
// points to are FStorage (physics/mp_thompson_init.hpp), as for WSM6.
//
// Reference object code: driver_microphysics, precip_from/to_mpas and microphysics_from/to_mpas as for WSM6
// (mpas_atmphys_microphysics.cpp); mp_gt_driver (__module_mp_thompson_MOD_mp_gt_driver: no libm call, no FMA;
// its 22 comisd / 6 divsd / 18 addsd ... are the copies, rho, nc1d = Nt_c/rho, the precipitation sums and th = t1d/pii;
// the standalone build differs from the reference only in the mpas_log_write call set-up, same arithmetic counts).
#include <cmath>
#include <cstddef>
#include <string>

#include "core/backend.hpp"
#include "core/error.hpp"
#include "core/libm.hpp"
#include "physics/mp_thompson.hpp"
#include "physics/mpas_atmphys_constants.hpp"
#include "physics/mpas_atmphys_interface.hpp"
#include "physics/mpas_atmphys_microphysics.hpp"

namespace cppmpas::phys {

namespace {
CPPMPAS_HD inline real fmax0(real x) { return x > 0.0 ? x : 0.0; }  // max(0._RKIND, x)
}  // namespace

void microp_timestep(const std::string& microp_scheme, real dt_dyn, real& dt_microp, int& n_microp) {
  dt_microp = dt_dyn;  // mpas_atmphys_manager.F:686-694
  n_microp = 1;
  if (microp_scheme == "mp_thompson" || microp_scheme == "mp_thompson_aerosols") {
    dt_microp = 90.;
    n_microp = static_cast<int>(std::lround(dt_dyn / dt_microp));  // nint
    n_microp = n_microp > 1 ? n_microp : 1;
    dt_microp = dt_dyn / n_microp;
    if (dt_dyn <= dt_microp) dt_microp = dt_dyn;
  }
}

void init_thompson_clouddroplets_forMPAS(int nCellsSolve, FView<const int> landmask, FView<real> nt_c,
                                         FView<real> mu_c) {
  parallel_for(1, nCellsSolve, [=] CPPMPAS_LAMBDA(int iCell) {  // IM:76-83
    if (landmask(iCell) == 1) {
      nt_c(iCell) = 300.e6;
    } else if (landmask(iCell) == 0) {
      nt_c(iCell) = 100.e6;
    }
    const real a = 1000.e6 / nt_c(iCell) + 2.;
    mu_c(iCell) = 15. < a ? 15. : a;  // MIN(15., (1000.e6/nt_c(iCell) + 2.))
  });
}

void driver_microphysics_thompson(const thompson::ThompsonView& tv, MicrophysicsPools& mp, int its, int ite) {
  using namespace physics_constants;
  MicrophysicsPools m = mp;  // the kernels capture the pool views by value
  if (m.l_acrain) throw Error("driver_microphysics: l_acrain (config_bucket_update) has no verified port");
  if (m.has_reqc != 0 || m.has_reqi != 0 || m.has_reqs != 0)
    throw Error("mp_gt_driver: has_reqc / has_reqi / has_reqs (config_microp_re, calc_effectRad) has no verified port");
  if (its != 1) throw Error("driver_microphysics: its /= 1");
  // an inactive constituent has MPAS's index -1 (scalar_index)
  if (m.index_qv <= 0 || m.index_qc <= 0 || m.index_qr <= 0 || m.index_qi <= 0 || m.index_qs <= 0 ||
      m.index_qg <= 0 || m.index_ni <= 0 || m.index_nr <= 0)
    throw Error("driver_microphysics: mp_thompson needs qv, qc, qr, qi, qs, qg, ni and nr");
  const int kts = 1, kte = m.nVertLevels;
  if (kte < 3 || kte > thompson::kMaxLevels) throw Error("mp_thompson: nVertLevels outside 3..128");
  if (m.mp_top_level < 1 || m.mp_top_level > kte) throw Error("driver_microphysics: mp_top_level outside 1..nVertLevels");
  if (m.n_microp < 1) throw Error("driver_microphysics: n_microp < 1");

  // allocate_microphysics (DRV:98-176): (ims:ime, kms:kme, jms:jme), ims = its = 1, ime = ite, kme = kte + 1
  const int ldi = ite;
  const std::size_t n3 = static_cast<std::size_t>(ldi) * (kte + 1);
  WorkArray rho_a(n3), th_a(n3), pi_a(n3), pres_a(n3), z_a(n3), dz_a(n3), w_a(n3), qv_a(n3), qc_a(n3),
      qr_a(n3), qi_a(n3), qs_a(n3), qg_a(n3), ni_a(n3), nr_a(n3), rainprod_a(n3), evapprod_a(n3), recloud_a(n3),
      reice_a(n3), resnow_a(n3);
  WorkArray rainnc_a(ldi), rainncv_a(ldi), sr_a(ldi), snownc_a(ldi), snowncv_a(ldi), graupelnc_a(ldi),
      graupelncv_a(ldi), ntc_a(ldi), muc_a(ldi);
  real *rho_p = rho_a.data(), *th_p = th_a.data(), *pi_p = pi_a.data(), *pres_p = pres_a.data(), *z_p = z_a.data(),
       *dz_p = dz_a.data(), *w_p = w_a.data(), *qv_p = qv_a.data(), *qc_p = qc_a.data(), *qr_p = qr_a.data(),
       *qi_p = qi_a.data(), *qs_p = qs_a.data(), *qg_p = qg_a.data(), *ni_p = ni_a.data(), *nr_p = nr_a.data(),
       *rainprod_p = rainprod_a.data(), *evapprod_p = evapprod_a.data(), *recloud_p = recloud_a.data(),
       *reice_p = reice_a.data(), *resnow_p = resnow_a.data();
  real *rainnc_p = rainnc_a.data(), *rainncv_p = rainncv_a.data(), *sr_p = sr_a.data(), *snownc_p = snownc_a.data(),
       *snowncv_p = snowncv_a.data(), *graupelnc_p = graupelnc_a.data(), *graupelncv_p = graupelncv_a.data(),
       *ntc_p = ntc_a.data(), *muc_p = muc_a.data();
  KernelCopy<thompson::ThompsonView> tv_k;
  const thompson::ThompsonView* T = kernel_ptr(tv, tv_k);
#define P(a, i, k) (a)[static_cast<std::size_t>((k) - 1) * ldi + ((i) - 1)]
#define qv(k, i) m.scalars(m.index_qv, k, i)
#define qc(k, i) m.scalars(m.index_qc, k, i)
#define qr(k, i) m.scalars(m.index_qr, k, i)
#define qi(k, i) m.scalars(m.index_qi, k, i)
#define qs(k, i) m.scalars(m.index_qs, k, i)
#define qg(k, i) m.scalars(m.index_qg, k, i)
#define ni(k, i) m.scalars(m.index_ni, k, i)
#define nr(k, i) m.scalars(m.index_nr, k, i)

  // precip_from_MPAS (DRV:530-562)
  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {
    rainncv_p[i - 1] = 0.;
    rainnc_p[i - 1] = 0.;
    m.rainncv(i) = 0.;
    snowncv_p[i - 1] = 0.;  // mp_thompson
    graupelncv_p[i - 1] = 0.;
    snownc_p[i - 1] = 0.;
    graupelnc_p[i - 1] = 0.;
    sr_p[i - 1] = 0.;
    m.snowncv(i) = 0.;
    m.graupelncv(i) = 0.;
    m.sr(i) = 0.;
  });

  // microphysics_from_MPAS (PI:598-823)
  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {
    for (int k = kts; k <= kte; ++k) {  // PI:656-675
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
    for (int k = kts; k <= kte; ++k) {  // PI:693-708
      P(qi_p, i, k) = qi(k, i);
      P(qs_p, i, k) = qs(k, i);
      P(qg_p, i, k) = qg(k, i);
      P(rainprod_p, i, k) = m.rainprod(k, i);
      P(evapprod_p, i, k) = m.evapprod(k, i);
      P(recloud_p, i, k) = m.re_cloud(k, i);
      P(reice_p, i, k) = m.re_ice(k, i);
      P(resnow_p, i, k) = m.re_snow(k, i);
    }
    muc_p[i - 1] = m.mu_c(i);  // PI:710-729
    ntc_p[i - 1] = m.nt_c(i);
    for (int k = kts; k <= kte; ++k) {
      P(ni_p, i, k) = ni(k, i);
      P(nr_p, i, k) = nr(k, i);
    }
    for (int k = kts; k <= kte; ++k) {  // PI:776-786
      m.rthmpten(k, i) = m.theta_m(k, i) / (1. + R_v / R_d * fmax0(qv(k, i)));
      m.rqvmpten(k, i) = qv(k, i);
      m.rqcmpten(k, i) = qc(k, i);
      m.rqrmpten(k, i) = qr(k, i);
      m.rqimpten(k, i) = qi(k, i);
      m.rqsmpten(k, i) = qs(k, i);
      m.rqgmpten(k, i) = qg(k, i);
    }
    for (int k = kts; k <= kte; ++k) {  // PI:792-797
      m.rnimpten(k, i) = ni(k, i);
      m.rnrmpten(k, i) = nr(k, i);
    }
  });

  // DRV:378-398: n_microp calls of mp_gt_driver (MT:1001-1387) with dt_in = dt_microp
  const real dt_microp = m.dt_microp;
  for (int istep = 1; istep <= m.n_microp; ++istep) {
    parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {
      constexpr int NK = thompson::kMaxLevels + 2;
      double qv1d[NK], qc1d[NK], qi1d[NK], qr1d[NK], qs1d[NK], qg1d[NK], ni1d[NK], nr1d[NK], nc1d[NK], t1d[NK],
          p1d[NK], w1d[NK], dz1d[NK], rho[NK], rainprod1d[NK], evapprod1d[NK];
      double pptrain = 0., pptsnow = 0., pptgraul = 0., pptice = 0.;  // MT:1133-1143
      rainncv_p[i - 1] = 0.;
      snowncv_p[i - 1] = 0.;
      graupelncv_p[i - 1] = 0.;
      sr_p[i - 1] = 0.;
      const double Nt_c = ntc_p[i - 1];  // MT:1144-1147 (mu_c is not read by mp_thompson)
      for (int k = kts; k <= kte; ++k) {
        rainprod1d[k] = 0.;
        evapprod1d[k] = 0.;
      }
      for (int k = kts; k <= kte; ++k) {  // MT:1153-1167
        t1d[k] = P(th_p, i, k) * P(pi_p, i, k);
        p1d[k] = P(pres_p, i, k);
        w1d[k] = P(w_p, i, k);
        dz1d[k] = P(dz_p, i, k);
        qv1d[k] = P(qv_p, i, k);
        qc1d[k] = P(qc_p, i, k);
        qi1d[k] = P(qi_p, i, k);
        qr1d[k] = P(qr_p, i, k);
        qs1d[k] = P(qs_p, i, k);
        qg1d[k] = P(qg_p, i, k);
        ni1d[k] = P(ni_p, i, k);
        nr1d[k] = P(nr_p, i, k);
        rho[k] = 0.622 * p1d[k] / (thompson::R * t1d[k] * (qv1d[k] + 0.622));
      }
      for (int k = kts; k <= kte; ++k) nc1d[k] = Nt_c / rho[k];  // MT:1175-1179 (nwfa1d / nifa1d unused)

      thompson::mp_thompson(*T, Nt_c, qv1d, qc1d, qi1d, qr1d, qs1d, qg1d, ni1d, nr1d, nc1d, t1d, p1d, w1d, dz1d,
                            pptrain, pptsnow, pptgraul, pptice, rainprod1d, evapprod1d, kts, kte, dt_microp);

      rainncv_p[i - 1] = pptrain + pptsnow + pptgraul + pptice;  // MT:1192-1206
      rainnc_p[i - 1] = rainnc_p[i - 1] + pptrain + pptsnow + pptgraul + pptice;
      snowncv_p[i - 1] = pptsnow + pptice;
      snownc_p[i - 1] = snownc_p[i - 1] + pptsnow + pptice;
      graupelncv_p[i - 1] = pptgraul;
      graupelnc_p[i - 1] = graupelnc_p[i - 1] + pptgraul;
      sr_p[i - 1] = (pptsnow + pptgraul + pptice) / (rainncv_p[i - 1] + 1.e-12);

      for (int k = kts; k <= kte; ++k) {  // MT:1226-1353 (the max / min trackers and the log messages are output only)
        P(qv_p, i, k) = qv1d[k];
        P(qc_p, i, k) = qc1d[k];
        P(qi_p, i, k) = qi1d[k];
        P(qr_p, i, k) = qr1d[k];
        P(qs_p, i, k) = qs1d[k];
        P(qg_p, i, k) = qg1d[k];
        P(ni_p, i, k) = ni1d[k];
        P(nr_p, i, k) = nr1d[k];
        P(th_p, i, k) = t1d[k] / P(pi_p, i, k);
        P(rainprod_p, i, k) = rainprod1d[k];
        P(evapprod_p, i, k) = evapprod1d[k];
        if (qv1d[k] < 0.0) {  // MT:1325-1333
          if (k < kte - 2 && k > kts + 1) {
            const double a = 0.5 * (P(qv_p, i, k - 1) + P(qv_p, i, k + 1));
            P(qv_p, i, k) = 1.E-7 > a ? 1.E-7 : a;  // MAX(1.E-7, ...)
          } else {
            P(qv_p, i, k) = 1.E-7;
          }
        }
      }
    });
  }

  if (m.l_diags) {  // DRV:467-482
    // compute_radar_reflectivity, mp_thompson (DRV:738-791)
    parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {
      constexpr int NK = thompson::kMaxLevels + 2;
      double p1d[NK], t1d[NK], qv1d[NK], qc1d[NK], qr1d[NK], qs1d[NK], qg1d[NK], nr1d[NK], dBZ1d[NK], zp[NK];
      for (int k = kts; k <= kte; ++k) {
        p1d[k] = P(pres_p, i, k);
        t1d[k] = P(th_p, i, k) * P(pi_p, i, k);
        qv1d[k] = P(qv_p, i, k);
        qc1d[k] = P(qc_p, i, k);
        qr1d[k] = P(qr_p, i, k);
        qs1d[k] = P(qs_p, i, k);
        qg1d[k] = P(qg_p, i, k);
        nr1d[k] = P(nr_p, i, k);
        dBZ1d[k] = -35.;
        zp[k] = P(z_p, i, k) - P(z_p, i, 1) + 0.5 * P(dz_p, i, k);  // height AGL
      }
      thompson::calc_refl10cm(*T, qv1d, qc1d, qr1d, nr1d, qs1d, qg1d, t1d, p1d, dBZ1d, kts, kte);
      int kp = 1;
      for (int k = kts; k <= kte; ++k) {
        dBZ1d[k] = -35. > dBZ1d[k] ? -35. : dBZ1d[k];  // max(-35._RKIND, dBZ1d(k))
        m.refl10cm(k, i) = dBZ1d[k];
        if (zp[k] < 1000.) kp = k;
      }
      real mx = dBZ1d[1];  // maxval(dBZ1d(:))
      for (int k = kts + 1; k <= kte; ++k)
        if (dBZ1d[k] > mx) mx = dBZ1d[k];
      m.refl10cm_max(i) = mx;
      const real w1 = (zp[kp + 1] - 1000.) / (zp[kp + 1] - zp[kp]);
      const real w2 = 1.0 - w1;
      m.refl10cm_1km(i) = w1 * dBZ1d[kp] + w2 * dBZ1d[kp + 1];
      // max(refl10cm_1km_max(i), refl10cm_1km(i)): maxsd keeps the first argument unless the second is greater (a
      // tie of -0.0 / +0.0 and a NaN keep refl10cm_1km_max)
      m.refl10cm_1km_max(i) = m.refl10cm_1km(i) > m.refl10cm_1km_max(i) ? m.refl10cm_1km(i) : m.refl10cm_1km_max(i);
    });
    // compute_relhum (DRV:798-850)
    parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {
      for (int k = kts; k <= kte; ++k) m.relhum(k, i) = 0.;
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
  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {
    m.precipw(i) = 0.;
    for (int k = kts; k <= kte; ++k) {
      const real rho_a = P(rho_p, i, k) / (1. + P(qv_p, i, k));
      m.precipw(i) = m.precipw(i) + P(qv_p, i, k) * rho_a * P(dz_p, i, k);
    }
    m.rainncv(i) = rainnc_p[i - 1];
    m.rainnc(i) = m.rainnc(i) + m.rainncv(i);
    m.snowncv(i) = snownc_p[i - 1];  // DRV:636-648
    m.graupelncv(i) = graupelnc_p[i - 1];
    m.sr(i) = (snownc_p[i - 1] + graupelnc_p[i - 1]) / (rainnc_p[i - 1] + 1.e-12);
    m.snownc(i) = m.snownc(i) + m.snowncv(i);
    m.graupelnc(i) = m.graupelnc(i) + m.graupelncv(i);
  });

  // microphysics_to_MPAS (PI:826-1138)
  const int mp_top_level = m.mp_top_level;
  const real dt_dyn = m.dt_dyn;
  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {
    for (int k = kts; k <= mp_top_level; ++k) {  // PI:903-934
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
      for (int k = mp_top_level + 1; k <= kte; ++k) {
        qc(k, i) = 0.;
        qr(k, i) = 0.;
        m.rt_diabatic_tend(k, i) = 0.;
      }
    {  // PI:950-964
      const real tem1 = m.zgrid(2, i) - m.zgrid(1, i);
      const real tem2 = m.zgrid(3, i) - m.zgrid(2, i);
      const real rho1 = m.rho_zz(1, i) * m.zz(1, i) * (1. + P(qv_p, i, 1));
      const real rho2 = m.rho_zz(2, i) * m.zz(2, i) * (1. + P(qv_p, i, 2));
      m.tend_sfc_pressure(i) = m.surface_pressure(i);
      m.surface_pressure(i) = 0.5 * gravity * (m.zgrid(2, i) - m.zgrid(1, i)) * (rho1 - 0.5 * (rho2 - rho1) * tem1 / (tem1 + tem2));
      m.surface_pressure(i) = m.surface_pressure(i) + m.pressure_p(1, i) + m.pressure_base(1, i);
      m.tend_sfc_pressure(i) = (m.surface_pressure(i) - m.tend_sfc_pressure(i)) / dt_dyn;
    }
    for (int k = kts; k <= mp_top_level; ++k) {  // PI:983-997
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
    for (int k = kts; k <= kte; ++k) {  // PI:1016-1030: ni, nr on every level
      ni(k, i) = P(ni_p, i, k);
      nr(k, i) = P(nr_p, i, k);
    }
    for (int k = kts; k <= mp_top_level; ++k) {  // PI:1076-1086
      m.rthmpten(k, i) = (m.theta_m(k, i) / (1. + R_v / R_d * fmax0(qv(k, i))) - m.rthmpten(k, i)) / dt_dyn;
      m.rqvmpten(k, i) = (qv(k, i) - m.rqvmpten(k, i)) / dt_dyn;
      m.rqcmpten(k, i) = (qc(k, i) - m.rqcmpten(k, i)) / dt_dyn;
      m.rqrmpten(k, i) = (qr(k, i) - m.rqrmpten(k, i)) / dt_dyn;
      m.rqimpten(k, i) = (qi(k, i) - m.rqimpten(k, i)) / dt_dyn;
      m.rqsmpten(k, i) = (qs(k, i) - m.rqsmpten(k, i)) / dt_dyn;
      m.rqgmpten(k, i) = (qg(k, i) - m.rqgmpten(k, i)) / dt_dyn;
    }
    if (kte > mp_top_level)  // PI:1089-1101
      for (int k = mp_top_level + 1; k <= kte; ++k) {
        m.rthmpten(k, i) = 0.;
        m.rqvmpten(k, i) = 0.;
        m.rqcmpten(k, i) = 0.;
        m.rqrmpten(k, i) = 0.;
        m.rqimpten(k, i) = 0.;
        m.rqsmpten(k, i) = 0.;
        m.rqgmpten(k, i) = 0.;
      }
    for (int k = kts; k <= kte; ++k) {  // PI:1104-1113
      m.rnimpten(k, i) = (ni(k, i) - m.rnimpten(k, i)) / dt_dyn;
      m.rnrmpten(k, i) = (nr(k, i) - m.rnrmpten(k, i)) / dt_dyn;
    }
  });
#undef P
#undef qv
#undef qc
#undef qr
#undef qi
#undef qs
#undef qg
#undef ni
#undef nr
}

}  // namespace cppmpas::phys
