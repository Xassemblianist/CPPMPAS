// SPDX-License-Identifier: GPL-3.0-or-later
// See mpas_atmphys_pbl.hpp. PB = src/core_atmosphere/physics/mpas_atmphys_driver_pbl.F, BY =
// physics/physics_wrf/module_bl_ysu.F (MPAS v8.4.2).
#include "physics/mpas_atmphys_pbl.hpp"

#include <vector>

#include "core/error.hpp"
#include "physics/bl_ysu.hpp"
#include "physics/mpas_atmphys_constants.hpp"

namespace cppmpas::phys {

void driver_pbl(int itimestep, PblPools& m, const PhysicsVars& v, int its, int ite) {
  using namespace physics_constants;
  (void)itimestep;  // initflag (PB:820-821) is MYNN's
  if (m.pbl_scheme != "bl_ysu")
    throw Error("driver_pbl: config_pbl_scheme = '" + m.pbl_scheme + "' has no verified port");
  if (its != 1) throw Error("driver_pbl: its /= 1");
  if (!m.f_qc || !m.f_qi)
    throw Error("driver_pbl: YSU without cloud water or cloud ice (f_qc / f_qi = .false.) has no verified port");
  const int kts = 1, kte = m.nVertLevels;
  if (kte > ysu::kMaxLevels) throw Error("driver_pbl: more than 128 levels");
  if (v.ime < ite || v.kme != kte + 1) throw Error("driver_pbl: MPAS_to_physics has not filled its..ite");
  const int n = ite, kme = kte + 1;

  // allocate_pbl (PB:87-193): (ims:ime, jms:jme) and (ims:ime, kms:kme, jms:jme)
  std::vector<real> hfx_p(n), qfx_p(n), ust_p(n), wspd_p(n), xland_p(n), hpbl_p(n), znt_p(n), uoce_p(n), voce_p(n);
  std::vector<int> kpbl_p(n);
  PArray rublten_p, rvblten_p, rthblten_p, rqvblten_p, rqcblten_p, rqiblten_p, rthraten_p, kzh_p, kzm_p, kzq_p,
      exch_p;
  for (PArray* a : {&rublten_p, &rvblten_p, &rthblten_p, &rqvblten_p, &rqcblten_p, &rqiblten_p, &rthraten_p, &kzh_p,
                    &kzm_p, &kzq_p, &exch_p})
    a->allocate(n, kme);
  std::vector<real> br_p(n), ctopo_p(n), ctopo2_p(n), delta_p(n), psih_p(n), psim_p(n), u10_p(n), v10_p(n),
      wstar_p(n);
  auto I = [](int i) { return static_cast<std::size_t>(i - 1); };

  // pbl_from_MPAS (PB:306-529)
  for (int i = its; i <= ite; ++i) {  // PB:357-371
    hfx_p[I(i)] = m.hfx(i);
    hpbl_p[I(i)] = m.hpbl(i);
    qfx_p[I(i)] = m.qfx(i);
    ust_p[I(i)] = m.ust(i);
    wspd_p[I(i)] = m.wspd(i);
    xland_p[I(i)] = m.xland(i);
    kpbl_p[I(i)] = 1;
    znt_p[I(i)] = m.znt(i);
    uoce_p[I(i)] = 0.;
    voce_p[I(i)] = 0.;
  }
  for (int k = kts; k <= kte; ++k)  // PB:372-376
    for (int i = its; i <= ite; ++i) rthraten_p(i, k) = m.rthratenlw(k, i) + m.rthratensw(k, i);
  const int ysu_pblmix = m.config_ysu_pblmix ? 1 : 0;  // PB:392-393
  for (int i = its; i <= ite; ++i) {  // PB:395-409
    br_p[I(i)] = m.br(i);
    psim_p[I(i)] = m.fm(i);
    psih_p[I(i)] = m.fh(i);
    u10_p[I(i)] = m.u10(i);
    v10_p[I(i)] = m.v10(i);
    delta_p[I(i)] = m.delta(i);
    wstar_p[I(i)] = m.wstar(i);
    ctopo_p[I(i)] = 1.;
    ctopo2_p[I(i)] = 1.;
  }
  for (int k = kts; k <= kte; ++k)  // PB:411-417
    for (int i = its; i <= ite; ++i) exch_p(i, k) = 0.;
  for (int k = kts; k <= kte; ++k)  // PB:512-527
    for (int i = its; i <= ite; ++i) {
      rublten_p(i, k) = 0.;
      rvblten_p(i, k) = 0.;
      rthblten_p(i, k) = 0.;
      rqvblten_p(i, k) = 0.;
      rqcblten_p(i, k) = 0.;
      rqiblten_p(i, k) = 0.;
      kzh_p(i, k) = 0.;
      kzm_p(i, k) = 0.;
      kzq_p(i, k) = 0.;
    }

  // PB:827-851: ysu (BY:12-470). The wrapper's _hv copies (BY:318-396) hold the same values; bl_ysu_run computes column
  // by column (bl_ysu.hpp), with the views on the (ims:ime, kms:kme) arrays: (i - 1, 0, k - 1), sk = ime.
  {
    using ysu::View2;
    using ysu::View3;
    const int ld = v.ime;
    auto c3 = [&](const PArray& a) { return View3<const real>{a.data(), 0, ld}; };
    auto w3 = [&](PArray& a) { return View3<real>{a.data(), 0, a.ldi}; };
    auto c2 = [](const std::vector<real>& a) { return View2<const real>{a.data(), 0}; };
    auto w2 = [](std::vector<real>& a) { return View2<real>{a.data(), 0}; };
    ysu::YSUKernel k{};
    k.ux = c3(v.u_p);
    k.vx = c3(v.v_p);
    k.tx = c3(v.t_p);
    k.qvx = c3(v.qv_p);
    k.qcx = c3(v.qc_p);
    k.qix = c3(v.qi_p);
    k.p2d = c3(v.pres_hyd_p);
    k.p2di = c3(v.pres2_hyd_p);
    k.pi2d = c3(v.pi_p);
    k.dz8w2d = c3(v.dz_p);
    k.rthraten = View3<const real>{rthraten_p.data(), 0, rthraten_p.ldi};
    k.psfcpa = View2<const real>{v.psfc_p.data(), 0};
    k.znt = c2(znt_p);
    k.ust = c2(ust_p);
    k.psim = c2(psim_p);
    k.psih = c2(psih_p);
    k.xland = c2(xland_p);
    k.hfx = c2(hfx_p);
    k.qfx = c2(qfx_p);
    k.wspd = c2(wspd_p);
    k.br = c2(br_p);
    k.u10 = w2(u10_p);
    k.v10 = w2(v10_p);
    k.uox = c2(uoce_p);
    k.vox = c2(voce_p);
    k.ctopo = c2(ctopo_p);
    k.ctopo2 = c2(ctopo2_p);
    k.utnp = w3(rublten_p);
    k.vtnp = w3(rvblten_p);
    k.ttnp = w3(rthblten_p);  // ttend / pi2d: BY:452
    k.qvtnp = w3(rqvblten_p);
    k.qctnp = w3(rqcblten_p);
    k.qitnp = w3(rqiblten_p);
    k.exch_hx = w3(kzh_p);
    k.exch_mx = w3(kzm_p);
    k.hpbl = w2(hpbl_p);
    k.kpbl1d = View2<int>{kpbl_p.data(), 0};
    k.wstar = w2(wstar_p);
    k.delta = w2(delta_p);
    k.dt = m.dt_pbl;
    k.nz = kte;
    k.f_qc = m.f_qc;
    k.f_qi = m.f_qi;
    k.ysu_topdown_pblmix = ysu_pblmix == 1;  // BY:315-316
    k.have_surface_diag = false;             // dusfc, dvsfc, dtsfc, dqsfc are not passed
    k.cp = cp;
    k.g = gravity;
    k.rovcp = rcp;
    k.rd = R_d;
    k.rovg = rdg;
    k.ep1 = ep_1;
    k.ep2 = ep_2;
    k.karman = karman;
    k.xlv = xlv;
    k.rv = R_v;
    k.ent_eff_mpas = true;  // MMM-physics a4baf7f bl_ysu.F90:855
    for (int i = its; i <= ite; ++i) k(i - 1, 0);
  }

  // pbl_to_MPAS (PB:532-724)
  for (int i = its; i <= ite; ++i) {  // PB:585-590
    m.hpbl(i) = hpbl_p[I(i)];
    m.kpbl(i) = kpbl_p[I(i)];
  }
  for (int k = kts; k <= kte; ++k)  // PB:592-607
    for (int i = its; i <= ite; ++i) {
      m.rublten(k, i) = rublten_p(i, k);
      m.rvblten(k, i) = rvblten_p(i, k);
      m.rthblten(k, i) = rthblten_p(i, k);
      m.rqvblten(k, i) = rqvblten_p(i, k);
      m.rqcblten(k, i) = rqcblten_p(i, k);
      m.rqiblten(k, i) = rqiblten_p(i, k);
      m.kzh(k, i) = kzh_p(i, k);
      m.kzm(k, i) = kzm_p(i, k);
      m.kzq(k, i) = kzh_p(i, k);  // PB:604: kzh_p, as MPAS
    }
  for (int i = its; i <= ite; ++i) {  // PB:616-626
    m.delta(i) = delta_p[I(i)];
    m.wstar(i) = wstar_p[I(i)];
  }
  for (int k = kts; k <= kte; ++k)
    for (int i = its; i <= ite; ++i) m.exch_h(k, i) = exch_p(i, k);
}

}  // namespace cppmpas::phys
