// SPDX-License-Identifier: GPL-3.0-or-later
// See mpas_atmphys_pbl.hpp. PB = src/core_atmosphere/physics/mpas_atmphys_driver_pbl.F, BY =
// physics/physics_wrf/module_bl_ysu.F (MPAS v8.4.2).
#include "physics/mpas_atmphys_pbl.hpp"

#include <vector>

#include "core/error.hpp"
#include "physics/bl_ysu.hpp"
#include "physics/mpas_atmphys_constants.hpp"

namespace cppmpas::phys {

namespace {
CPPMPAS_HD inline std::size_t I(int i) { return static_cast<std::size_t>(i - 1); }
}  // namespace

// One kernel per loop of PB (one thread per cell i, the k loops inside the thread).
void driver_pbl(int itimestep, PblPools& mp, const PhysicsVars& vp, int its, int ite) {
  using namespace physics_constants;
  (void)itimestep;  // initflag (PB:820-821) is MYNN's
  if (mp.pbl_scheme != "bl_ysu")
    throw Error("driver_pbl: config_pbl_scheme = '" + mp.pbl_scheme + "' has no verified port");
  if (its != 1) throw Error("driver_pbl: its /= 1");
  if (!mp.f_qc || !mp.f_qi)
    throw Error("driver_pbl: YSU without cloud water or cloud ice (f_qc / f_qi = .false.) has no verified port");
  const int kts = 1, kte = mp.nVertLevels;
  if (kte > ysu::kMaxLevels) throw Error("driver_pbl: more than 128 levels");
  if (vp.ime < ite || vp.kme != kte + 1) throw Error("driver_pbl: MPAS_to_physics has not filled its..ite");
  const PblArrays m = mp;
  const int n = ite, kme = kte + 1;

  // allocate_pbl (PB:87-193): (ims:ime, jms:jme) and (ims:ime, kms:kme, jms:jme)
  WorkArray hfx_a(n), qfx_a(n), ust_a(n), wspd_a(n), xland_a(n), hpbl_a(n), znt_a(n), uoce_a(n), voce_a(n);
  WorkArrayInt kpbl_a(n);
  PArray rublten_a, rvblten_a, rthblten_a, rqvblten_a, rqcblten_a, rqiblten_a, rthraten_a, kzh_a, kzm_a, kzq_a,
      exch_a;
  for (PArray* a : {&rublten_a, &rvblten_a, &rthblten_a, &rqvblten_a, &rqcblten_a, &rqiblten_a, &rthraten_a, &kzh_a,
                    &kzm_a, &kzq_a, &exch_a})
    a->allocate(n, kme);
  WorkArray br_a(n), ctopo_a(n), ctopo2_a(n), delta_a(n), psih_a(n), psim_a(n), u10_a(n), v10_a(n), wstar_a(n);
  real *hfx_p = hfx_a.data(), *qfx_p = qfx_a.data(), *ust_p = ust_a.data(), *wspd_p = wspd_a.data(),
       *xland_p = xland_a.data(), *hpbl_p = hpbl_a.data(), *znt_p = znt_a.data(), *uoce_p = uoce_a.data(),
       *voce_p = voce_a.data(), *br_p = br_a.data(), *ctopo_p = ctopo_a.data(), *ctopo2_p = ctopo2_a.data(),
       *delta_p = delta_a.data(), *psih_p = psih_a.data(), *psim_p = psim_a.data(), *u10_p = u10_a.data(),
       *v10_p = v10_a.data(), *wstar_p = wstar_a.data();
  int* kpbl_p = kpbl_a.data();
  const PView rublten_p = rublten_a.view(), rvblten_p = rvblten_a.view(), rthblten_p = rthblten_a.view(),
              rqvblten_p = rqvblten_a.view(), rqcblten_p = rqcblten_a.view(), rqiblten_p = rqiblten_a.view(),
              rthraten_p = rthraten_a.view(), kzh_p = kzh_a.view(), kzm_p = kzm_a.view(), kzq_p = kzq_a.view(),
              exch_p = exch_a.view();

  // pbl_from_MPAS (PB:306-529)
  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {  // PB:357-371
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
  });
  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {  // PB:372-376
    for (int k = kts; k <= kte; ++k) rthraten_p(i, k) = m.rthratenlw(k, i) + m.rthratensw(k, i);
  });
  const int ysu_pblmix = m.config_ysu_pblmix ? 1 : 0;  // PB:392-393
  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {  // PB:395-409
    br_p[I(i)] = m.br(i);
    psim_p[I(i)] = m.fm(i);
    psih_p[I(i)] = m.fh(i);
    u10_p[I(i)] = m.u10(i);
    v10_p[I(i)] = m.v10(i);
    delta_p[I(i)] = m.delta(i);
    wstar_p[I(i)] = m.wstar(i);
    ctopo_p[I(i)] = 1.;
    ctopo2_p[I(i)] = 1.;
  });
  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {  // PB:411-417
    for (int k = kts; k <= kte; ++k) exch_p(i, k) = 0.;
  });
  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {  // PB:512-527
    for (int k = kts; k <= kte; ++k) {
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
  });

  // PB:827-851: ysu (BY:12-470). The wrapper's _hv copies (BY:318-396) hold the same values; bl_ysu_run computes column
  // by column (bl_ysu.hpp), with the views on the (ims:ime, kms:kme) arrays: (i - 1, 0, k - 1), sk = ime.
  {
    using ysu::View2;
    using ysu::View3;
    const PhysicsVars& v = vp;
    const int ld = v.ime;
    auto c3 = [&](const PArray& a) { return View3<const real>{a.data(), 0, ld}; };
    auto w3 = [&](PArray& a) { return View3<real>{a.data(), 0, a.ldi}; };
    auto c2 = [](const WorkArray& a) { return View2<const real>{a.data(), 0}; };
    auto w2 = [](WorkArray& a) { return View2<real>{a.data(), 0}; };
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
    k.rthraten = View3<const real>{rthraten_a.data(), 0, rthraten_a.ldi};
    k.psfcpa = View2<const real>{v.psfc_p.data(), 0};
    k.znt = c2(znt_a);
    k.ust = c2(ust_a);
    k.psim = c2(psim_a);
    k.psih = c2(psih_a);
    k.xland = c2(xland_a);
    k.hfx = c2(hfx_a);
    k.qfx = c2(qfx_a);
    k.wspd = c2(wspd_a);
    k.br = c2(br_a);
    k.u10 = w2(u10_a);
    k.v10 = w2(v10_a);
    k.uox = c2(uoce_a);
    k.vox = c2(voce_a);
    k.ctopo = c2(ctopo_a);
    k.ctopo2 = c2(ctopo2_a);
    k.utnp = w3(rublten_a);
    k.vtnp = w3(rvblten_a);
    k.ttnp = w3(rthblten_a);  // ttend / pi2d: BY:452
    k.qvtnp = w3(rqvblten_a);
    k.qctnp = w3(rqcblten_a);
    k.qitnp = w3(rqiblten_a);
    k.exch_hx = w3(kzh_a);
    k.exch_mx = w3(kzm_a);
    k.hpbl = w2(hpbl_a);
    k.kpbl1d = View2<int>{kpbl_a.data(), 0};
    k.wstar = w2(wstar_a);
    k.delta = w2(delta_a);
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
    parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) { k(i - 1, 0); });
  }

  // pbl_to_MPAS (PB:532-724)
  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {  // PB:585-590
    m.hpbl(i) = hpbl_p[I(i)];
    m.kpbl(i) = kpbl_p[I(i)];
  });
  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {  // PB:592-607
    for (int k = kts; k <= kte; ++k) {
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
  });
  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {  // PB:616-626
    m.delta(i) = delta_p[I(i)];
    m.wstar(i) = wstar_p[I(i)];
  });
  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {
    for (int k = kts; k <= kte; ++k) m.exch_h(k, i) = exch_p(i, k);
  });
}

}  // namespace cppmpas::phys
