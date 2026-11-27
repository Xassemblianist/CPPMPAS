// SPDX-License-Identifier: GPL-3.0-or-later
// See mpas_atmphys_convection.hpp. CV = src/core_atmosphere/physics/mpas_atmphys_driver_convection.F (MPAS v8.4.2).
#include "physics/mpas_atmphys_convection.hpp"

#include "core/error.hpp"
#include "core/libm.hpp"
#include "physics/cu_ntiedtke.hpp"
#include "physics/module_cu_gf.hpp"
#include "physics/module_cu_ntiedtke.hpp"
#include "physics/mpas_atmphys_constants.hpp"

namespace cppmpas::phys {

namespace {
// CPPWRF's lane rule (CPPWRF 7fac5bf src/core/backend.hpp:197)
CPPMPAS_HD inline bool lane(int i, int n) { return n > 1 && !((n & 1) && i == n - 1); }

void check_scheme(const std::string& s, const char* where) {
  if (s != "cu_ntiedtke" && s != "cu_grell_freitas")
    throw Error(std::string(where) + ": config_convection_scheme = '" + s + "' has no verified port");
}
}  // namespace

void update_convection_step1(const ConvectionPools& m, int its, int ite) {
  (void)its;
  (void)ite;
  check_scheme(m.convection_scheme, "update_convection_step1");  // CV:1014-1056: only cu_kain_fritsch acts
}

namespace {
// CV:405-621 for config_convection_scheme = 'cu_grell_freitas'. One kernel per loop of CV (one thread per cell i,
// the k loops inside the thread).
void driver_convection_gf(int itimestep, ConvectionPools& mp, const PhysicsVars& vp, int its, int ite) {
  using namespace physics_constants;
  const int kts = 1, kte = mp.nVertLevels;
  const int n = ite, kme = kte + 1;
  const ConvectionArrays m = mp;
  const PhysicsVarsView v = vp.view();
  // allocate_convection (CV:103-176): the common arrays, zeroed, and the GF arrays (CV:146-176)
  PArray pratec_a, raincv_a, rthcuten_a, rqvcuten_a, rqccuten_a, rqicuten_a;
  PArray ht_a, cubot_a, cutop_a, dx_a, area_a, gsw_a, hfx_a, qfx_a, xland_a, xmb_total_a, xmb_shallow_a;
  PArray rthblten_a, rthdynten_a, rthraten_a, rqvblten_a, rqvdynten_a, rucuten_a, rvcuten_a, qccu_a, qicu_a;
  IArray kpbl_a, k22_shallow_a, kbcon_shallow_a, ktop_shallow_a, kbot_shallow_a, ktop_deep_a;
  for (PArray* x : {&pratec_a, &raincv_a, &ht_a, &cubot_a, &cutop_a, &dx_a, &area_a, &gsw_a, &hfx_a, &qfx_a,
                    &xland_a, &xmb_total_a, &xmb_shallow_a})
    x->allocate(n, 1);
  for (PArray* x : {&rthcuten_a, &rqvcuten_a, &rqccuten_a, &rqicuten_a, &rthblten_a, &rthdynten_a, &rthraten_a,
                    &rqvblten_a, &rqvdynten_a, &rucuten_a, &rvcuten_a, &qccu_a, &qicu_a})
    x->allocate(n, kme);
  for (IArray* x : {&kpbl_a, &k22_shallow_a, &kbcon_shallow_a, &ktop_shallow_a, &kbot_shallow_a, &ktop_deep_a})
    x->allocate(n);
  const PView pratec_p = pratec_a.view(), raincv_p = raincv_a.view(), rthcuten_p = rthcuten_a.view(),
              rqvcuten_p = rqvcuten_a.view(), rqccuten_p = rqccuten_a.view(), rqicuten_p = rqicuten_a.view();
  const PView ht_p = ht_a.view(), cubot_p = cubot_a.view(), cutop_p = cutop_a.view(), dx_p = dx_a.view(),
              area_p = area_a.view(), gsw_p = gsw_a.view(), hfx_p = hfx_a.view(), qfx_p = qfx_a.view(),
              xland_p = xland_a.view(), xmb_total_p = xmb_total_a.view(), xmb_shallow_p = xmb_shallow_a.view();
  const PView rthblten_p = rthblten_a.view(), rthdynten_p = rthdynten_a.view(), rthraten_p = rthraten_a.view(),
              rqvblten_p = rqvblten_a.view(), rqvdynten_p = rqvdynten_a.view(), rucuten_p = rucuten_a.view(),
              rvcuten_p = rvcuten_a.view(), qccu_p = qccu_a.view(), qicu_p = qicu_a.view();
  const IView kpbl_p = kpbl_a.view(), k22_shallow_p = k22_shallow_a.view(), kbcon_shallow_p = kbcon_shallow_a.view(),
              ktop_shallow_p = ktop_shallow_a.view(), ktop_deep_p = ktop_deep_a.view();

  // convection_from_MPAS (CV:624-860), called with dt_dyn (CV:461)
  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {  // CV:668-679
    raincv_p(i) = m.raincv(i);
    pratec_p(i) = m.cuprec(i);
    for (int k = kts; k <= kte; ++k) {
      rthcuten_p(i, k) = m.rthcuten(k, i);
      rqvcuten_p(i, k) = m.rqvcuten(k, i);
      rqccuten_p(i, k) = m.rqccuten(k, i);
      rqicuten_p(i, k) = m.rqicuten(k, i);
    }
  });
  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {  // CV:716-750
    dx_p(i) = m.config_len_disp / libm::pow(m.meshDensity(i), 0.25);  // scalar pow@plt (header)
    area_p(i) = m.areaCell(i);
    ht_p(i) = m.zgrid(1, i);
    xland_p(i) = m.xland(i);
    gsw_p(i) = m.gsw(i);
    hfx_p(i) = m.hfx(i);
    qfx_p(i) = m.qfx(i);
    kpbl_p(i) = m.kpbl(i);
    cubot_p(i) = m.cubot(i);
    cutop_p(i) = m.cutop(i);
    xmb_total_p(i) = m.xmb_total(i);
    xmb_shallow_p(i) = m.xmb_shallow(i);
    k22_shallow_p(i) = m.k22_shallow(i);
    kbcon_shallow_p(i) = m.kbcon_shallow(i);
    ktop_shallow_p(i) = m.ktop_shallow(i);
    ktop_deep_p(i) = m.ktop_deep(i);
    for (int k = kts; k <= kte; ++k) {
      qccu_p(i, k) = m.qc_cu(k, i);
      qicu_p(i, k) = m.qi_cu(k, i);
      rthblten_p(i, k) = m.rthblten(k, i);
      rthdynten_p(i, k) = m.rthdynten(k, i);
      rthraten_p(i, k) = m.rthratenlw(k, i) + m.rthratensw(k, i);
      rqvblten_p(i, k) = m.rqvblten(k, i);
      rqvdynten_p(i, k) = m.rqvdynten(k, i);
      rucuten_p(i, k) = m.rucuten(k, i);
      rvcuten_p(i, k) = m.rvcuten(k, i);
    }
  });

  // CV:482-514
  GfArgs a;
  a.itimestep = itimestep;
  a.dt = m.dt_dyn;
  a.dxCell = dx_p.cview();
  a.areaCell = area_p.cview();
  a.u = v.u_p.cview();
  a.v = v.v_p.cview();
  a.w = v.w_p.cview();
  a.t = v.t_p.cview();
  a.q = v.qv_p.cview();
  a.rho = v.rho_p.cview();
  a.p = v.pres_hyd_p.cview();
  a.pi = v.pi_p.cview();
  a.p8w = v.pres2_hyd_p.cview();
  a.dz8w = v.dz_p.cview();
  a.ht = ht_p.cview();
  a.xland = xland_p.cview();
  a.gsw = gsw_p.cview();
  a.xlv = xlv;
  a.cp = cp;
  a.g = gravity;
  a.r_v = R_v;
  a.hfx = hfx_p.cview();
  a.qfx = qfx_p.cview();
  a.rthblten = rthblten_p.cview();
  a.rqvblten = rqvblten_p.cview();
  a.rthften = rthdynten_p.cview();
  a.rqvften = rqvdynten_p.cview();
  a.rthraten = rthraten_p.cview();
  a.kpbl = kpbl_p;
  a.raincv = raincv_p;
  a.pratec = pratec_p;
  a.htop = cutop_p;
  a.hbot = cubot_p;
  a.k22_shallow = k22_shallow_p;
  a.kbcon_shallow = kbcon_shallow_p;
  a.ktop_shallow = ktop_shallow_p;
  a.xmb_total = xmb_total_p;
  a.xmb_shallow = xmb_shallow_p;
  a.gdc = qccu_p;
  a.gdc2 = qicu_p;
  a.rthcuten = rthcuten_p;
  a.rqvcuten = rqvcuten_p;
  a.rqccuten = rqccuten_p;
  a.rqicuten = rqicuten_p;
  a.ichoice_deep = m.config_gfconv_closure_deep;
  a.ichoice_shallow = m.config_gfconv_closure_shallow;
  a.ishallow_g3 = 1;  // ishallow (mpas_atmphys_vars.F:340)
  a.kts = kts;
  a.kte = kte;
  cu_grell_freitas(a, its, ite);

  // convection_to_MPAS (CV:863-988)
  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {  // CV:899-910
    m.raincv(i) = raincv_p(i);
    m.cuprec(i) = pratec_p(i);
    for (int k = kts; k <= kte; ++k) {
      m.rthcuten(k, i) = rthcuten_p(i, k);
      m.rqvcuten(k, i) = rqvcuten_p(i, k);
      m.rqccuten(k, i) = rqccuten_p(i, k);
      m.rqicuten(k, i) = rqicuten_p(i, k);
    }
  });
  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {  // CV:929-947
    m.cubot(i) = cubot_p(i);
    m.cutop(i) = cutop_p(i);
    m.xmb_total(i) = xmb_total_p(i);
    m.xmb_shallow(i) = xmb_shallow_p(i);
    m.k22_shallow(i) = k22_shallow_p(i);
    m.kbcon_shallow(i) = kbcon_shallow_p(i);
    m.ktop_shallow(i) = ktop_shallow_p(i);
    m.ktop_deep(i) = ktop_deep_p(i);
    for (int k = kts; k <= kte; ++k) {
      m.qc_cu(k, i) = qccu_p(i, k);
      m.qi_cu(k, i) = qicu_p(i, k);
      m.rucuten(k, i) = rucuten_p(i, k);
      m.rvcuten(k, i) = rvcuten_p(i, k);
    }
  });
}
}  // namespace

// One kernel per loop of CV (one thread per cell i, the k loops inside the thread).
void driver_convection(int itimestep, ConvectionPools& mp, const PhysicsVars& vp, int its, int ite) {
  using namespace physics_constants;
  check_scheme(mp.convection_scheme, "driver_convection");
  if (mp.convection_scheme == "cu_grell_freitas") {
    if (its != 1) throw Error("driver_convection: its /= 1");
    if (mp.nVertLevels > cu_gf::kMaxLevels) throw Error("driver_convection: more than 100 levels");
    if (vp.ime < ite || vp.kme != mp.nVertLevels + 1)
      throw Error("driver_convection: MPAS_to_physics has not filled its..ite");
    if (!mp.areaCell.p || !mp.zgrid.p || !mp.gsw.p || !mp.kpbl.p || !mp.cubot.p || !mp.cutop.p || !mp.xmb_total.p ||
        !mp.xmb_shallow.p || !mp.qc_cu.p || !mp.qi_cu.p || !mp.k22_shallow.p || !mp.kbcon_shallow.p ||
        !mp.ktop_shallow.p || !mp.ktop_deep.p)
      throw Error("driver_convection: the cu_grell_freitas_in package arrays are not set");
    driver_convection_gf(itimestep, mp, vp, its, ite);
    return;
  }
  if (its != 1) throw Error("driver_convection: its /= 1");
  const int kts = 1, kte = mp.nVertLevels;
  if (kte > cu_ntiedtke::kMaxLevels) throw Error("driver_convection: more than 100 levels");
  if (vp.ime < ite || vp.kme != kte + 1) throw Error("driver_convection: MPAS_to_physics has not filled its..ite");
  const ConvectionArrays m = mp;
  const int n = ite, kme = kte + 1;
  const PhysicsVarsView v = vp.view();

  // allocate_convection (CV:103-265) for cu_ntiedtke: (ims:ime, jms:jme) and (ims:ime, kms:kme, jms:jme), zeroed
  // where CV zeroes them; every element of its..ite is set again below
  PArray pratec_a, raincv_a, hfx_a, qfx_a, xland_a, dx_a;
  PArray rthcuten_a, rqvcuten_a, rqccuten_a, rqicuten_a, rucuten_a, rvcuten_a, rqvften_a, rthften_a;
  for (PArray* a : {&pratec_a, &raincv_a, &hfx_a, &qfx_a, &xland_a, &dx_a}) a->allocate(n, 1);
  for (PArray* a : {&rthcuten_a, &rqvcuten_a, &rqccuten_a, &rqicuten_a, &rucuten_a, &rvcuten_a, &rqvften_a,
                    &rthften_a})
    a->allocate(n, kme);
  const PView pratec_p = pratec_a.view(), raincv_p = raincv_a.view(), hfx_p = hfx_a.view(), qfx_p = qfx_a.view(),
              xland_p = xland_a.view(), dx_p = dx_a.view();
  const PView rthcuten_p = rthcuten_a.view(), rqvcuten_p = rqvcuten_a.view(), rqccuten_p = rqccuten_a.view(),
              rqicuten_p = rqicuten_a.view(), rucuten_p = rucuten_a.view(), rvcuten_p = rvcuten_a.view(),
              rqvften_p = rqvften_a.view(), rthften_p = rthften_a.view();

  // convection_from_MPAS (CV:624-860), called with dt_dyn (CV:461)
  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {  // CV:668-679
    raincv_p(i) = m.raincv(i);
    pratec_p(i) = m.cuprec(i);
    for (int k = kts; k <= kte; ++k) {
      rthcuten_p(i, k) = m.rthcuten(k, i);
      rqvcuten_p(i, k) = m.rqvcuten(k, i);
      rqccuten_p(i, k) = m.rqccuten(k, i);
      rqicuten_p(i, k) = m.rqicuten(k, i);
    }
  });
  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {  // CV:800-805
    xland_p(i) = m.xland(i);
    hfx_p(i) = m.hfx(i);
    qfx_p(i) = m.qfx(i);
  });
  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {  // CV:807-812
    for (int k = kts; k <= kte; ++k) {
      rucuten_p(i, k) = m.rucuten(k, i);
      rvcuten_p(i, k) = m.rvcuten(k, i);
    }
  });
  // CV:838-841: vectorised over i (_ZGVbN2vv_pow lanes, mpas_atmphys_convection.hpp)
  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {
    dx_p(i) = m.config_len_disp / (lane(i - its, ite - its + 1) ? libm::mvec_pow(m.meshDensity(i), 0.25)
                                                                 : libm::pow(m.meshDensity(i), 0.25));
  });
  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {  // CV:843-849
    for (int k = kts; k <= kte; ++k) {
      rqvften_p(i, k) = m.rqvdynten(k, i) + m.rqvblten(k, i);
      rthften_p(i, k) =
          (m.rthdynten(k, i) + m.rthblten(k, i) + m.rthratenlw(k, i) + m.rthratensw(k, i)) * v.pi_p(i, k);
    }
  });

  // CV:476-478
  int initflag = 1;
  if (m.config_do_restart || itimestep > 1) initflag = 0;

  // CV:580-608
  NtiedtkeArgs a;
  a.dt = m.dt_dyn;
  a.itimestep = initflag;
  a.stepcu = m.n_cu;
  a.raincv = raincv_p;
  a.pratec = pratec_p;
  a.qfx = qfx_a.cview();
  a.hfx = hfx_a.cview();
  a.xland = xland_a.cview();
  a.dx = dx_a.cview();
  a.u3d = v.u_p.cview();
  a.v3d = v.v_p.cview();
  a.w = v.w_p.cview();
  a.t3d = v.t_p.cview();
  a.qv3d = v.qv_p.cview();
  a.qc3d = v.qc_p.cview();
  a.qi3d = v.qi_p.cview();
  a.pi3d = v.pi_p.cview();
  a.rho3d = v.rho_p.cview();
  a.qvften = rqvften_a.cview();
  a.thften = rthften_a.cview();
  a.dz8w = v.dz_p.cview();
  a.pcps = v.pres_hyd_p.cview();
  a.p8w = v.pres2_hyd_p.cview();
  a.f_qc = m.f_qc;
  a.f_qi = m.f_qi;
  a.rthcuten = rthcuten_p;
  a.rqvcuten = rqvcuten_p;
  a.rqccuten = rqccuten_p;
  a.rqicuten = rqicuten_p;
  a.rucuten = rucuten_p;
  a.rvcuten = rvcuten_p;
  a.grav = gravity;
  a.xlf = xlf;
  a.xls = xls;
  a.xlv = xlv;
  a.rd = R_d;
  a.rv = R_v;
  a.cp = cp;
  a.kts = kts;
  a.kte = kte;
  cu_ntiedtke_driver(a, its, ite);

  // convection_to_MPAS (CV:863-988)
  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {  // CV:899-910
    m.raincv(i) = raincv_p(i);
    m.cuprec(i) = pratec_p(i);
    for (int k = kts; k <= kte; ++k) {
      m.rthcuten(k, i) = rthcuten_p(i, k);
      m.rqvcuten(k, i) = rqvcuten_p(i, k);
      m.rqccuten(k, i) = rqccuten_p(i, k);
      m.rqicuten(k, i) = rqicuten_p(i, k);
    }
  });
  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {  // CV:975-982
    for (int k = kts; k <= kte; ++k) {
      m.rucuten(k, i) = rucuten_p(i, k);
      m.rvcuten(k, i) = rvcuten_p(i, k);
    }
  });
}

void update_convection_step2(ConvectionPools& mp, int its, int ite) {
  const ConvectionArrays m = mp;
  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {  // CV:1089-1097
    m.rainc(i) = m.rainc(i) + m.cuprec(i) * m.dt_dyn;
    if (m.l_acrain && m.config_bucket_rainc > 0. && m.rainc(i) > m.config_bucket_rainc) {
      m.i_rainc(i) = m.i_rainc(i) + 1;
      m.rainc(i) = m.rainc(i) - m.config_bucket_rainc;
    }
  });
}

}  // namespace cppmpas::phys
