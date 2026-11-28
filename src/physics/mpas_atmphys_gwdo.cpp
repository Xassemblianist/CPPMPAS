// SPDX-License-Identifier: GPL-3.0-or-later
// See mpas_atmphys_gwdo.hpp. GD = src/core_atmosphere/physics/mpas_atmphys_driver_gwdo.F, BG =
// physics/physics_wrf/module_bl_gwdo.F (MPAS v8.4.2).
#include "physics/mpas_atmphys_gwdo.hpp"


#include "core/error.hpp"
#include "core/libm.hpp"
#include "physics/bl_gwdo.hpp"
#include "physics/mpas_atmphys_constants.hpp"

namespace cppmpas::phys {

// One kernel per loop of GD and one thread per column for gwdo (one thread per cell i, the k loops inside).
void driver_gwdo(int itimestep, GwdoPools& mp, const PhysicsVars& vp, int its, int ite) {
  using namespace physics_constants;
  (void)itimestep;  // gwdo's itimestep argument is not used by bl_gwdo_run
  if (mp.gwdo_scheme != "bl_ysu_gwdo")
    throw Error("driver_gwdo: config_gwdo_scheme = '" + mp.gwdo_scheme + "' has no verified port");
  if (mp.config_ngw_scheme) throw Error("driver_gwdo: config_ngw_scheme = true has no verified port");
  if (its != 1) throw Error("driver_gwdo: its /= 1");
  const int kts = 1, kte = mp.nVertLevels;
  if (kte > bl_gwdo::kMaxLevels) throw Error("driver_gwdo: more than 128 levels");
  // bl_gwdo's kbl >= kpblmin = 2 reads prsl(2) (bl_gwdo.F90:298-307): one level is past MPAS's column
  if (kte < 2) throw Error("driver_gwdo: fewer than 2 levels (bl_gwdo reads level kpblmin = 2)");
  if (vp.ime < ite || vp.kme != kte + 1) throw Error("driver_gwdo: MPAS_to_physics has not filled its..ite");
  const GwdoArrays m = mp;
  const PhysicsVarsView v = vp.view();
  const int n = ite, kme = kte + 1;

  // allocate_gwdo (GD:75-115)
  PArray cosa_a, sina_a, dx_a, dusfcg_a, dvsfcg_a, var2d_a, con_a, oa1_a, oa2_a, oa3_a, oa4_a, ol1_a, ol2_a, ol3_a,
      ol4_a;
  WorkArrayInt kpbl_a(n);
  for (PArray* a : {&cosa_a, &sina_a, &dx_a, &dusfcg_a, &dvsfcg_a, &var2d_a, &con_a, &oa1_a, &oa2_a, &oa3_a, &oa4_a,
                    &ol1_a, &ol2_a, &ol3_a, &ol4_a})
    a->allocate(n, 1);
  PArray dtaux3d_a, dtauy3d_a, rublten_a, rvblten_a, rthblten_a;
  for (PArray* a : {&dtaux3d_a, &dtauy3d_a, &rublten_a, &rvblten_a, &rthblten_a}) a->allocate(n, kme);
  const PView cosa_p = cosa_a.view(), sina_p = sina_a.view(), dx_p = dx_a.view(), dusfcg_p = dusfcg_a.view(),
              dvsfcg_p = dvsfcg_a.view(), var2d_p = var2d_a.view(), con_p = con_a.view(), oa1_p = oa1_a.view(),
              oa2_p = oa2_a.view(), oa3_p = oa3_a.view(), oa4_p = oa4_a.view(), ol1_p = ol1_a.view(),
              ol2_p = ol2_a.view(), ol3_p = ol3_a.view(), ol4_p = ol4_a.view();
  const PView dtaux3d_p = dtaux3d_a.view(), dtauy3d_p = dtauy3d_a.view(), rublten_p = rublten_a.view(),
              rvblten_p = rvblten_a.view(), rthblten_p = rthblten_a.view();
  int* kpbl_p = kpbl_a.data();

  // gwdo_from_MPAS (GD:285-552)
  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {  // GD:350-363
    var2d_p(i) = m.var2d(i);
    con_p(i) = m.con(i);
    oa1_p(i) = m.oa1(i);
    oa2_p(i) = m.oa2(i);
    oa3_p(i) = m.oa3(i);
    oa4_p(i) = m.oa4(i);
    ol1_p(i) = m.ol1(i);
    ol2_p(i) = m.ol2(i);
    ol3_p(i) = m.ol3(i);
    ol4_p(i) = m.ol4(i);
  });
  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {  // GD:529-538 (scalar pow@plt, GD:533)
    sina_p(i) = 0.;
    cosa_p(i) = 1.;
    dx_p(i) = m.config_len_disp / libm::pow(m.meshDensity(i), 0.25);
    kpbl_p[i - 1] = m.kpbl(i);
    dusfcg_p(i) = m.dusfcg(i);
    dvsfcg_p(i) = m.dvsfcg(i);
  });
  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {  // GD:540-550
    for (int k = kts; k <= kte; ++k) {
      dtaux3d_p(i, k) = m.dtaux3d(k, i);
      dtauy3d_p(i, k) = m.dtauy3d(k, i);
      rublten_p(i, k) = m.rublten(k, i);
      rvblten_p(i, k) = m.rvblten(k, i);
      rthblten_p(i, k) = m.rthblten(k, i);
    }
  });

  // GD:757-778: gwdo (BG:16-238) for j = 1: the _hv copies of its..ite (BG:158-196) and bl_gwdo_run (BG:198-214) on
  // each column (the columns of bl_gwdo_run are independent; no vectorised libm call), then BG:218-234
  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {
    bl_gwdo::GwdoColumn<real> c;
    for (int k = kts; k <= kte + 1; ++k) c.prsi[k] = v.pres2_hydd_p(i, k);  // p3di = pres2_hydd_p
    for (int k = kts; k <= kte; ++k) {
      c.rublten[k] = rublten_p(i, k);
      c.rvblten[k] = rvblten_p(i, k);
      c.uproj[k] = v.u_p(i, k);
      c.vproj[k] = v.v_p(i, k);
      c.t1[k] = v.t_p(i, k);
      c.q1[k] = v.qv_p(i, k);
      c.prsl[k] = v.pres_hydd_p(i, k);  // p3d = pres_hydd_p
      c.prslk[k] = v.pi_p(i, k);
      c.zl[k] = v.zmid_p(i, k);  // z = zmid_p
    }
    c.dxmeter = dx_p(i);
    c.var = var2d_p(i);
    c.oc1 = con_p(i);
    c.sina = sina_p(i);
    c.cosa = cosa_p(i);
    c.oa2d1 = oa1_p(i);
    c.oa2d2 = oa2_p(i);
    c.oa2d3 = oa3_p(i);
    c.oa2d4 = oa4_p(i);
    c.ol2d1 = ol1_p(i);
    c.ol2d2 = ol2_p(i);
    c.ol2d3 = ol3_p(i);
    c.ol2d4 = ol4_p(i);
    // g_ = gravity, cp_ = cp, rd_ = R_d, fv_ = ep_1, pi_ = pii, deltim = dt_pbl (GD:760-765, BG:212-213); rv_ is
    // not used by bl_gwdo_run
    bl_gwdo::bl_gwdo_run<real>(c, kte, gravity, cp, R_d, ep_1, pii, m.dt_pbl);
    for (int k = kts; k <= kte; ++k) {
      rublten_p(i, k) = c.rublten[k];
      rvblten_p(i, k) = c.rvblten[k];
      dtaux3d_p(i, k) = c.dtaux3d[k];
      dtauy3d_p(i, k) = c.dtauy3d[k];
    }
    dusfcg_p(i) = c.dusfcg;
    dvsfcg_p(i) = c.dvsfcg;
  });
  (void)kpbl_p;  // kpbl2d is not used by gwdo (BG)

  // gwdo_to_MPAS (GD:555-689)
  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {  // GD:668-673
    m.dusfcg(i) = dusfcg_p(i);
    m.dvsfcg(i) = dvsfcg_p(i);
  });
  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {  // GD:675-687
    for (int k = kts; k <= kte; ++k) {
      m.dtaux3d(k, i) = dtaux3d_p(i, k);
      m.dtauy3d(k, i) = dtauy3d_p(i, k);
      m.rubldiff(k, i) = rublten_p(i, k) - m.rublten(k, i);
      m.rvbldiff(k, i) = rvblten_p(i, k) - m.rvblten(k, i);
      m.rublten(k, i) = rublten_p(i, k);
      m.rvblten(k, i) = rvblten_p(i, k);
      m.rthblten(k, i) = rthblten_p(i, k);
    }
  });
}

}  // namespace cppmpas::phys
