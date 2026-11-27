// SPDX-License-Identifier: GPL-3.0-or-later
// See mpas_atmphys_gwdo.hpp. GD = src/core_atmosphere/physics/mpas_atmphys_driver_gwdo.F, BG =
// physics/physics_wrf/module_bl_gwdo.F (MPAS v8.4.2).
#include "physics/mpas_atmphys_gwdo.hpp"

#include <memory>

#include "core/error.hpp"
#include "core/libm.hpp"
#include "physics/bl_gwdo.hpp"
#include "physics/mpas_atmphys_constants.hpp"

namespace cppmpas::phys {

void driver_gwdo(int itimestep, GwdoPools& m, const PhysicsVars& v, int its, int ite) {
  using namespace physics_constants;
  (void)itimestep;  // gwdo's itimestep argument is not used by bl_gwdo_run
  if (m.gwdo_scheme != "bl_ysu_gwdo")
    throw Error("driver_gwdo: config_gwdo_scheme = '" + m.gwdo_scheme + "' has no verified port");
  if (m.config_ngw_scheme) throw Error("driver_gwdo: config_ngw_scheme = true has no verified port");
  if (its != 1) throw Error("driver_gwdo: its /= 1");
  const int kts = 1, kte = m.nVertLevels;
  if (kte > bl_gwdo::kMaxLevels) throw Error("driver_gwdo: more than 128 levels");
  if (v.ime < ite || v.kme != kte + 1) throw Error("driver_gwdo: MPAS_to_physics has not filled its..ite");
  const int n = ite, kme = kte + 1;

  // allocate_gwdo (GD:75-115)
  PArray cosa_p, sina_p, dx_p, dusfcg_p, dvsfcg_p, var2d_p, con_p, oa1_p, oa2_p, oa3_p, oa4_p, ol1_p, ol2_p, ol3_p,
      ol4_p;
  std::unique_ptr<int[]> kpbl_p(new int[n]);
  for (PArray* a : {&cosa_p, &sina_p, &dx_p, &dusfcg_p, &dvsfcg_p, &var2d_p, &con_p, &oa1_p, &oa2_p, &oa3_p, &oa4_p,
                    &ol1_p, &ol2_p, &ol3_p, &ol4_p})
    a->allocate(n, 1);
  PArray dtaux3d_p, dtauy3d_p, rublten_p, rvblten_p, rthblten_p;
  for (PArray* a : {&dtaux3d_p, &dtauy3d_p, &rublten_p, &rvblten_p, &rthblten_p}) a->allocate(n, kme);

  // gwdo_from_MPAS (GD:285-552)
  for (int i = its; i <= ite; ++i) {  // GD:350-363
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
  }
  for (int i = its; i <= ite; ++i) {  // GD:529-538 (scalar pow@plt, GD:533)
    sina_p(i) = 0.;
    cosa_p(i) = 1.;
    dx_p(i) = m.config_len_disp / libm::pow(m.meshDensity(i), 0.25);
    kpbl_p[i - 1] = m.kpbl(i);
    dusfcg_p(i) = m.dusfcg(i);
    dvsfcg_p(i) = m.dvsfcg(i);
  }
  for (int k = kts; k <= kte; ++k)  // GD:540-550
    for (int i = its; i <= ite; ++i) {
      dtaux3d_p(i, k) = m.dtaux3d(k, i);
      dtauy3d_p(i, k) = m.dtauy3d(k, i);
      rublten_p(i, k) = m.rublten(k, i);
      rvblten_p(i, k) = m.rvblten(k, i);
      rthblten_p(i, k) = m.rthblten(k, i);
    }

  // GD:757-778: gwdo (BG:16-238) for j = 1: the _hv copies of its..ite (BG:158-196) and bl_gwdo_run (BG:198-214) on
  // each column (the columns of bl_gwdo_run are independent; no vectorised libm call), then BG:218-234
  for (int i = its; i <= ite; ++i) {
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
  }
  (void)kpbl_p;  // kpbl2d is not used by gwdo (BG)

  // gwdo_to_MPAS (GD:555-689)
  for (int i = its; i <= ite; ++i) {  // GD:668-673
    m.dusfcg(i) = dusfcg_p(i);
    m.dvsfcg(i) = dvsfcg_p(i);
  }
  for (int k = kts; k <= kte; ++k)  // GD:675-687
    for (int i = its; i <= ite; ++i) {
      m.dtaux3d(k, i) = dtaux3d_p(i, k);
      m.dtauy3d(k, i) = dtauy3d_p(i, k);
      m.rubldiff(k, i) = rublten_p(i, k) - m.rublten(k, i);
      m.rvbldiff(k, i) = rvblten_p(i, k) - m.rvblten(k, i);
      m.rublten(k, i) = rublten_p(i, k);
      m.rvblten(k, i) = rvblten_p(i, k);
      m.rthblten(k, i) = rthblten_p(i, k);
    }
}

}  // namespace cppmpas::phys
