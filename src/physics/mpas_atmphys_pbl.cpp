// SPDX-License-Identifier: GPL-3.0-or-later
// See mpas_atmphys_pbl.hpp. PB = src/core_atmosphere/physics/mpas_atmphys_driver_pbl.F, BY =
// physics/physics_wrf/module_bl_ysu.F, MB = physics/physics_wrf/module_bl_mynn.F (MPAS v8.4.2).
#include "physics/mpas_atmphys_pbl.hpp"

#include <vector>

#include "core/error.hpp"
#include "core/libm.hpp"
#include "physics/bl_mynn.hpp"
#include "physics/bl_ysu.hpp"
#include "physics/mpas_atmphys_constants.hpp"

namespace cppmpas::phys {

namespace {
CPPMPAS_HD inline std::size_t I(int i) { return static_cast<std::size_t>(i - 1); }

mynn::V3 v3(const PArray& a) { return mynn::V3{const_cast<real*>(a.data()), a.ldi}; }
mynn::V2 v2(WorkArray& a) { return mynn::V2{a.data()}; }

// PB:757-973 with config_pbl_scheme = 'bl_mynn'
void driver_pbl_mynn(int itimestep, PblPools& mp, const PhysicsVars& vp, int its, int ite) {
  using namespace physics_constants;
  if (its != 1) throw Error("driver_pbl: its /= 1");
  if (const char* r = mp.mynn.refused())
    throw Error(std::string("driver_pbl: bl_mynn with ") + r + " has no verified port");
  if (!mp.f_qc || !mp.f_qi || !mp.f_qs || !mp.f_ni)
    throw Error("driver_pbl: bl_mynn without qc, qi, qs or ni (f_qc / f_qi / f_qs / f_ni = .false.) has no verified port");
  if (mp.f_nc || mp.f_nifa || mp.f_nwfa || mp.f_qoz || mp.f_nbca)
    throw Error("driver_pbl: bl_mynn with nc, nifa, nwfa (mp_thompson_aerosols), qoz or nbca has no verified port");
  const int kts = 1, kte = mp.nVertLevels;
  if (kte > mynn::kMaxLevels) throw Error("driver_pbl: bl_mynn with more than mynn::kMaxLevels levels");
  if (kte < 4) throw Error("driver_pbl: bl_mynn with fewer than 4 levels");
  if (vp.ime < ite || vp.kme != kte + 1) throw Error("driver_pbl: MPAS_to_physics has not filled its..ite");
  if (!vp.ni_p.data() || !vp.nc_p.data()) throw Error("driver_pbl: bl_mynn needs MPAS_to_physics' nc_p .. nwfa_p");
  const PblArrays m = mp;
  const int n = ite, kme = kte + 1;

  // allocate_pbl (PB:87-193): the arrays of every scheme and the bl_mynn ones
  WorkArray hfx_a(n), qfx_a(n), ust_a(n), wspd_a(n), xland_a(n), hpbl_a(n), znt_a(n), uoce_a(n), voce_a(n);
  WorkArrayInt kpbl_a(n), kbl_plume_a(n);
  WorkArray dx_a(n), ch_a(n), qsfc_a(n), rmol_a(n), tsk_a(n), maxwidthbl_a(n), maxmfbl_a(n), zbl_plume_a(n);
  PArray rublten_a, rvblten_a, rthblten_a, rqvblten_a, rqcblten_a, rqiblten_a, rthraten_a, kzh_a, kzm_a, kzq_a;
  PArray cov_a, qke_a, qsq_a, tsq_a, qkeadv_a, elpbl_a, tkepbl_a, sh3d_a, sm3d_a, dqke_a, qbuoy_a, qdiss_a, qshear_a,
      qwt_a, qcbl_a, qibl_a, cldfrabl_a, edmfa_a, edmfw_a, edmfqt_a, edmfthl_a, edmfent_a, edmfqc_a, subthl_a,
      subqv_a, detthl_a, detqv_a, rqsblten_a, rncblten_a, rniblten_a, rnifablten_a, rnwfablten_a, pattern_spp_a;
  for (PArray* a : {&rublten_a, &rvblten_a, &rthblten_a, &rqvblten_a, &rqcblten_a, &rqiblten_a, &rthraten_a, &kzh_a,
                    &kzm_a, &kzq_a, &cov_a, &qke_a, &qsq_a, &tsq_a, &qkeadv_a, &elpbl_a, &tkepbl_a, &sh3d_a,
                    &sm3d_a, &dqke_a, &qbuoy_a, &qdiss_a, &qshear_a, &qwt_a, &qcbl_a, &qibl_a, &cldfrabl_a,
                    &edmfa_a, &edmfw_a, &edmfqt_a, &edmfthl_a, &edmfent_a, &edmfqc_a, &subthl_a, &subqv_a,
                    &detthl_a, &detqv_a, &rqsblten_a, &rncblten_a, &rniblten_a, &rnifablten_a, &rnwfablten_a,
                    &pattern_spp_a})
    a->allocate(n, kme);
  real *hfx_p = hfx_a.data(), *qfx_p = qfx_a.data(), *ust_p = ust_a.data(), *wspd_p = wspd_a.data(),
       *xland_p = xland_a.data(), *hpbl_p = hpbl_a.data(), *znt_p = znt_a.data(), *uoce_p = uoce_a.data(),
       *voce_p = voce_a.data(), *dx_p = dx_a.data(), *ch_p = ch_a.data(), *qsfc_p = qsfc_a.data(),
       *rmol_p = rmol_a.data(), *tsk_p = tsk_a.data(), *maxwidthbl_p = maxwidthbl_a.data(),
       *maxmfbl_p = maxmfbl_a.data(), *zbl_plume_p = zbl_plume_a.data();
  int *kpbl_p = kpbl_a.data(), *kbl_plume_p = kbl_plume_a.data();
  const PView rthraten_p = rthraten_a.view(), elpbl_p = elpbl_a.view(), cov_p = cov_a.view(), qke_p = qke_a.view(),
              qsq_p = qsq_a.view(), tsq_p = tsq_a.view(), tkepbl_p = tkepbl_a.view(), qkeadv_p = qkeadv_a.view(),
              sh3d_p = sh3d_a.view(), sm3d_p = sm3d_a.view(), cldfrabl_p = cldfrabl_a.view(),
              qcbl_p = qcbl_a.view(), qibl_p = qibl_a.view(), edmfa_p = edmfa_a.view(), edmfent_p = edmfent_a.view(),
              edmfqc_p = edmfqc_a.view(), edmfqt_p = edmfqt_a.view(), edmfthl_p = edmfthl_a.view(),
              edmfw_p = edmfw_a.view(), subthl_p = subthl_a.view(), subqv_p = subqv_a.view(),
              detthl_p = detthl_a.view(), detqv_p = detqv_a.view(), dqke_p = dqke_a.view(), qbuoy_p = qbuoy_a.view(),
              qdiss_p = qdiss_a.view(), qshear_p = qshear_a.view(), qwt_p = qwt_a.view(),
              rqsblten_p = rqsblten_a.view(), rncblten_p = rncblten_a.view(), rniblten_p = rniblten_a.view(),
              rnifablten_p = rnifablten_a.view(), rnwfablten_p = rnwfablten_a.view(),
              pattern_spp_pbl = pattern_spp_a.view(), rublten_p = rublten_a.view(), rvblten_p = rvblten_a.view(),
              rthblten_p = rthblten_a.view(), rqvblten_p = rqvblten_a.view(), rqcblten_p = rqcblten_a.view(),
              rqiblten_p = rqiblten_a.view(), kzh_p = kzh_a.view(), kzm_p = kzm_a.view(), kzq_p = kzq_a.view();

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
  const real len_disp = m.config_len_disp;
  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {  // PB:451-458 (meshDensity**0.25: scalar pow@plt)
    dx_p[I(i)] = len_disp / libm::pow(m.meshDensity(i), 0.25);
    ch_p[I(i)] = m.ch(i);
    qsfc_p[I(i)] = m.qsfc(i);
    rmol_p[I(i)] = m.rmol(i);
    tsk_p[I(i)] = m.skintemp(i);
  });
  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {  // PB:460-505
    for (int k = kts; k <= kte; ++k) {
      elpbl_p(i, k) = m.el_pbl(k, i);
      cov_p(i, k) = m.cov(k, i);
      qke_p(i, k) = m.qke(k, i);
      qsq_p(i, k) = m.qsq(k, i);
      tsq_p(i, k) = m.tsq(k, i);
      tkepbl_p(i, k) = m.tke_pbl(k, i);
      qkeadv_p(i, k) = m.qke_adv(k, i);
      sh3d_p(i, k) = m.sh3d(k, i);
      sm3d_p(i, k) = m.sm3d(k, i);
      cldfrabl_p(i, k) = m.cldfrac_bl(k, i);
      qcbl_p(i, k) = m.qc_bl(k, i);
      qibl_p(i, k) = m.qi_bl(k, i);
      edmfa_p(i, k) = m.edmf_a(k, i);
      edmfent_p(i, k) = m.edmf_ent(k, i);
      edmfqc_p(i, k) = m.edmf_qc(k, i);
      edmfqt_p(i, k) = m.edmf_qt(k, i);
      edmfthl_p(i, k) = m.edmf_thl(k, i);
      edmfw_p(i, k) = m.edmf_w(k, i);
      subthl_p(i, k) = m.sub_thl(k, i);
      subqv_p(i, k) = m.sub_qv(k, i);
      detthl_p(i, k) = m.det_thl(k, i);
      detqv_p(i, k) = m.det_qv(k, i);
      dqke_p(i, k) = 0.;
      qbuoy_p(i, k) = 0.;
      qdiss_p(i, k) = 0.;
      qshear_p(i, k) = 0.;
      qwt_p(i, k) = 0.;
      rqsblten_p(i, k) = 0.;
      rncblten_p(i, k) = 0.;
      rniblten_p(i, k) = 0.;
      rnifablten_p(i, k) = 0.;
      rnwfablten_p(i, k) = 0.;
      pattern_spp_pbl(i, k) = 0.;
    }
    kbl_plume_p[I(i)] = 0;
    maxwidthbl_p[I(i)] = 0.;
    maxmfbl_p[I(i)] = 0.;
    zbl_plume_p[I(i)] = 0;
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

  // PB:820-821
  int initflag = 1;
  if (m.config_do_restart || itimestep > 1) initflag = 0;

  // PB:904-966 mynn_bl_driver (MB:19-746; physics/bl_mynn.hpp, one thread per column)
  {
    const PhysicsVars& v = vp;
    mynn::BlMynnKernel K{};
    K.C = m.mynn_const;
    K.o = m.mynn;
    K.kte = kte;
    K.delt = m.dt_pbl;
    K.initflag = initflag;
    K.do_restart = m.config_do_restart;
    K.f_qc = m.f_qc;
    K.f_qi = m.f_qi;
    K.f_qs = m.f_qs;
    K.f_qoz = m.f_qoz;
    K.f_nc = m.f_nc;
    K.f_ni = m.f_ni;
    K.f_nifa = m.f_nifa;
    K.f_nwfa = m.f_nwfa;
    K.f_nbca = m.f_nbca;
    K.dx = v2(dx_a);
    K.xland = v2(xland_a);
    K.ps = mynn::V2{const_cast<real*>(v.psfc_p.data())};
    K.ts = v2(tsk_a);
    K.qsfc = v2(qsfc_a);
    K.ust = v2(ust_a);
    K.ch = v2(ch_a);
    K.hfx = v2(hfx_a);
    K.qfx = v2(qfx_a);
    K.rmol = v2(rmol_a);
    K.wspd = v2(wspd_a);
    K.uoce = v2(uoce_a);
    K.voce = v2(voce_a);
    K.znt = v2(znt_a);
    K.dz = v3(v.dz_p);
    K.u = v3(v.u_p);
    K.v = v3(v.v_p);
    K.w = v3(v.w_p);
    K.th = v3(v.th_p);
    K.tt = v3(v.t_p);
    K.p = v3(v.pres_hyd_p);
    K.exner = v3(v.pi_p);
    K.rho = v3(v.rho_p);
    K.qv = v3(v.qv_p);
    K.qc = v3(v.qc_p);
    K.qi = v3(v.qi_p);
    K.qs = v3(v.qs_p);
    K.nc = v3(v.nc_p);
    K.ni = v3(v.ni_p);
    K.nifa = v3(v.nifa_p);
    K.nwfa = v3(v.nwfa_p);
    K.rthraten = v3(rthraten_a);
    K.pblh = v2(hpbl_a);
    K.kpbl = mynn::V2i{kpbl_a.data()};
    K.maxwidth = v2(maxwidthbl_a);
    K.maxmf = v2(maxmfbl_a);
    K.ktop_plume = mynn::V2i{kbl_plume_a.data()};
    K.ztop_plume = v2(zbl_plume_a);
    K.cldfra_bl = v3(cldfrabl_a);
    K.qc_bl = v3(qcbl_a);
    K.qi_bl = v3(qibl_a);
    K.dqke = v3(dqke_a);
    K.qke_adv = v3(qkeadv_a);
    K.tsq = v3(tsq_a);
    K.qsq = v3(qsq_a);
    K.cov = v3(cov_a);
    K.el_pbl = v3(elpbl_a);
    K.rublten = v3(rublten_a);
    K.rvblten = v3(rvblten_a);
    K.rthblten = v3(rthblten_a);
    K.rqvblten = v3(rqvblten_a);
    K.rqcblten = v3(rqcblten_a);
    K.rqiblten = v3(rqiblten_a);
    K.rqsblten = v3(rqsblten_a);
    K.rncblten = v3(rncblten_a);
    K.rniblten = v3(rniblten_a);
    K.rnifablten = v3(rnifablten_a);
    K.rnwfablten = v3(rnwfablten_a);
    K.edmf_a = v3(edmfa_a);
    K.edmf_w = v3(edmfw_a);
    K.edmf_qt = v3(edmfqt_a);
    K.edmf_thl = v3(edmfthl_a);
    K.edmf_ent = v3(edmfent_a);
    K.edmf_qc = v3(edmfqc_a);
    K.sub_thl = v3(subthl_a);
    K.sub_sqv = v3(subqv_a);
    K.det_thl = v3(detthl_a);
    K.det_sqv = v3(detqv_a);
    K.exch_h = v3(kzh_a);
    K.exch_m = v3(kzm_a);
    K.qke = v3(qke_a);
    K.qwt = v3(qwt_a);
    K.qshear = v3(qshear_a);
    K.qbuoy = v3(qbuoy_a);
    K.qdiss = v3(qdiss_a);
    K.sh3d = v3(sh3d_a);
    K.sm3d = v3(sm3d_a);
    parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) { K(i); });
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
  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {  // PB:659-690
    for (int k = kts; k <= kte; ++k) {
      m.el_pbl(k, i) = elpbl_p(i, k);
      m.cov(k, i) = cov_p(i, k);
      m.qke(k, i) = qke_p(i, k);
      m.qsq(k, i) = qsq_p(i, k);
      m.tsq(k, i) = tsq_p(i, k);
      m.sh3d(k, i) = sh3d_p(i, k);
      m.sm3d(k, i) = sm3d_p(i, k);
      m.tke_pbl(k, i) = tkepbl_p(i, k);
      m.qke_adv(k, i) = qkeadv_p(i, k);
      m.cldfrac_bl(k, i) = cldfrabl_p(i, k);
      m.qc_bl(k, i) = qcbl_p(i, k);
      m.qi_bl(k, i) = qibl_p(i, k);
      m.edmf_a(k, i) = edmfa_p(i, k);
      m.edmf_ent(k, i) = edmfent_p(i, k);
      m.edmf_qc(k, i) = edmfqc_p(i, k);
      m.edmf_qt(k, i) = edmfqt_p(i, k);
      m.edmf_thl(k, i) = edmfthl_p(i, k);
      m.edmf_w(k, i) = edmfw_p(i, k);
      m.sub_thl(k, i) = subthl_p(i, k);
      m.sub_qv(k, i) = subqv_p(i, k);
      m.det_thl(k, i) = detthl_p(i, k);
      m.det_qv(k, i) = detqv_p(i, k);
      m.dqke(k, i) = dqke_p(i, k);
      m.qbuoy(k, i) = qbuoy_p(i, k);
      m.qdiss(k, i) = qdiss_p(i, k);
      m.qshear(k, i) = qshear_p(i, k);
      m.qwt(k, i) = qwt_p(i, k);
      m.rqsblten(k, i) = rqsblten_p(i, k);
    }
  });
  if (m.f_ni)  // PB:692-700
    parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {
      for (int k = kts; k <= kte; ++k) m.rniblten(k, i) = rniblten_p(i, k);
    });
  // PB:701-711 (f_nc .and. f_nifa .and. f_nwfa): refused above
  if (m.diag_kbl_plume.p)  // the column check's copy of the plume diagnostics (not in MPAS)
    for (int i = its; i <= ite; ++i) {
      m.diag_kbl_plume(i) = kbl_plume_a[I(i)];
      m.diag_maxwidthbl(i) = maxwidthbl_a[I(i)];
      m.diag_maxmfbl(i) = maxmfbl_a[I(i)];
      m.diag_zbl_plume(i) = zbl_plume_a[I(i)];
    }
}
}  // namespace

mynn::Const mynn_init_pbl() {  // PB:744-747
  using namespace physics_constants;
  return mynn::bl_mynn_init(cp, cpv, cice, cliq, ep_1, ep_2, gravity, karman, P0, R_d, R_v, svp1, svp2, svp3, svpt0,
                            xlf, xls, xlv);
}

// One kernel per loop of PB (one thread per cell i, the k loops inside the thread).
void driver_pbl(int itimestep, PblPools& mp, const PhysicsVars& vp, int its, int ite) {
  using namespace physics_constants;
  if (mp.pbl_scheme == "bl_mynn") {
    driver_pbl_mynn(itimestep, mp, vp, its, ite);
    return;
  }
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
