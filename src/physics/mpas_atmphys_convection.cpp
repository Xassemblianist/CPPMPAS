// SPDX-License-Identifier: GPL-3.0-or-later
// See mpas_atmphys_convection.hpp. CV = src/core_atmosphere/physics/mpas_atmphys_driver_convection.F (MPAS v8.4.2).
#include "physics/mpas_atmphys_convection.hpp"

#include "core/error.hpp"
#include "core/libm.hpp"
#include "physics/cu_ntiedtke.hpp"
#include "physics/module_cu_ntiedtke.hpp"
#include "physics/mpas_atmphys_constants.hpp"

namespace cppmpas::phys {

namespace {
// CPPWRF's lane rule (CPPWRF 7fac5bf src/core/backend.hpp:197)
inline bool lane(int i, int n) { return n > 1 && !((n & 1) && i == n - 1); }

void check_scheme(const std::string& s, const char* where) {
  if (s != "cu_ntiedtke")
    throw Error(std::string(where) + ": config_convection_scheme = '" + s + "' has no verified port");
}
}  // namespace

void update_convection_step1(const ConvectionPools& m, int its, int ite) {
  (void)its;
  (void)ite;
  check_scheme(m.convection_scheme, "update_convection_step1");  // CV:1014-1056: only cu_kain_fritsch acts
}

void driver_convection(int itimestep, ConvectionPools& m, const PhysicsVars& v, int its, int ite) {
  using namespace physics_constants;
  check_scheme(m.convection_scheme, "driver_convection");
  if (its != 1) throw Error("driver_convection: its /= 1");
  const int kts = 1, kte = m.nVertLevels;
  if (kte > cu_ntiedtke::kMaxLevels) throw Error("driver_convection: more than 100 levels");
  if (v.ime < ite || v.kme != kte + 1) throw Error("driver_convection: MPAS_to_physics has not filled its..ite");
  const int n = ite, kme = kte + 1;

  // allocate_convection (CV:103-265) for cu_ntiedtke: (ims:ime, jms:jme) and (ims:ime, kms:kme, jms:jme), zeroed
  // where CV zeroes them; every element of its..ite is set again below
  PArray pratec_p, raincv_p, hfx_p, qfx_p, xland_p, dx_p;
  PArray rthcuten_p, rqvcuten_p, rqccuten_p, rqicuten_p, rucuten_p, rvcuten_p, rqvften_p, rthften_p;
  for (PArray* a : {&pratec_p, &raincv_p, &hfx_p, &qfx_p, &xland_p, &dx_p}) a->allocate(n, 1);
  for (PArray* a : {&rthcuten_p, &rqvcuten_p, &rqccuten_p, &rqicuten_p, &rucuten_p, &rvcuten_p, &rqvften_p,
                    &rthften_p})
    a->allocate(n, kme);

  // convection_from_MPAS (CV:624-860), called with dt_dyn (CV:461)
  for (int i = its; i <= ite; ++i) {  // CV:668-679
    raincv_p(i) = m.raincv(i);
    pratec_p(i) = m.cuprec(i);
    for (int k = kts; k <= kte; ++k) {
      rthcuten_p(i, k) = m.rthcuten(k, i);
      rqvcuten_p(i, k) = m.rqvcuten(k, i);
      rqccuten_p(i, k) = m.rqccuten(k, i);
      rqicuten_p(i, k) = m.rqicuten(k, i);
    }
  }
  for (int i = its; i <= ite; ++i) {  // CV:800-805
    xland_p(i) = m.xland(i);
    hfx_p(i) = m.hfx(i);
    qfx_p(i) = m.qfx(i);
  }
  for (int k = kts; k <= kte; ++k)  // CV:807-812
    for (int i = its; i <= ite; ++i) {
      rucuten_p(i, k) = m.rucuten(k, i);
      rvcuten_p(i, k) = m.rvcuten(k, i);
    }
  // CV:838-841: vectorised over i (_ZGVbN2vv_pow lanes, mpas_atmphys_convection.hpp)
  for (int i = its; i <= ite; ++i)
    dx_p(i) = m.config_len_disp / (lane(i - its, ite - its + 1) ? libm::mvec_pow(m.meshDensity(i), 0.25)
                                                                  : libm::pow(m.meshDensity(i), 0.25));
  for (int k = kts; k <= kte; ++k)  // CV:843-849
    for (int i = its; i <= ite; ++i) {
      rqvften_p(i, k) = m.rqvdynten(k, i) + m.rqvblten(k, i);
      rthften_p(i, k) =
          (m.rthdynten(k, i) + m.rthblten(k, i) + m.rthratenlw(k, i) + m.rthratensw(k, i)) * v.pi_p(i, k);
    }

  // CV:476-478
  int initflag = 1;
  if (m.config_do_restart || itimestep > 1) initflag = 0;

  // CV:580-608
  NtiedtkeArgs a;
  a.dt = m.dt_dyn;
  a.itimestep = initflag;
  a.stepcu = m.n_cu;
  a.raincv = &raincv_p;
  a.pratec = &pratec_p;
  a.qfx = &qfx_p;
  a.hfx = &hfx_p;
  a.xland = &xland_p;
  a.dx = &dx_p;
  a.u3d = &v.u_p;
  a.v3d = &v.v_p;
  a.w = &v.w_p;
  a.t3d = &v.t_p;
  a.qv3d = &v.qv_p;
  a.qc3d = &v.qc_p;
  a.qi3d = &v.qi_p;
  a.pi3d = &v.pi_p;
  a.rho3d = &v.rho_p;
  a.qvften = &rqvften_p;
  a.thften = &rthften_p;
  a.dz8w = &v.dz_p;
  a.pcps = &v.pres_hyd_p;
  a.p8w = &v.pres2_hyd_p;
  a.f_qc = m.f_qc;
  a.f_qi = m.f_qi;
  a.rthcuten = &rthcuten_p;
  a.rqvcuten = &rqvcuten_p;
  a.rqccuten = &rqccuten_p;
  a.rqicuten = &rqicuten_p;
  a.rucuten = &rucuten_p;
  a.rvcuten = &rvcuten_p;
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
  for (int i = its; i <= ite; ++i) {  // CV:899-910
    m.raincv(i) = raincv_p(i);
    m.cuprec(i) = pratec_p(i);
    for (int k = kts; k <= kte; ++k) {
      m.rthcuten(k, i) = rthcuten_p(i, k);
      m.rqvcuten(k, i) = rqvcuten_p(i, k);
      m.rqccuten(k, i) = rqccuten_p(i, k);
      m.rqicuten(k, i) = rqicuten_p(i, k);
    }
  }
  for (int k = kts; k <= kte; ++k)  // CV:975-982
    for (int i = its; i <= ite; ++i) {
      m.rucuten(k, i) = rucuten_p(i, k);
      m.rvcuten(k, i) = rvcuten_p(i, k);
    }
}

void update_convection_step2(ConvectionPools& m, int its, int ite) {
  for (int i = its; i <= ite; ++i) {  // CV:1089-1097
    m.rainc(i) = m.rainc(i) + m.cuprec(i) * m.dt_dyn;
    if (m.l_acrain && m.config_bucket_rainc > 0. && m.rainc(i) > m.config_bucket_rainc) {
      m.i_rainc(i) = m.i_rainc(i) + 1;
      m.rainc(i) = m.rainc(i) - m.config_bucket_rainc;
    }
  }
}

}  // namespace cppmpas::phys
