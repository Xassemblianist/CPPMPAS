// SPDX-License-Identifier: GPL-3.0-or-later
// See mpas_atmphys_interface.hpp. PI = src/core_atmosphere/physics/mpas_atmphys_interface.F (MPAS v8.4.2).
#include "physics/mpas_atmphys_interface.hpp"

#include <cstdio>

#include "core/error.hpp"
#include "core/libm.hpp"
#include "physics/mpas_atmphys_constants.hpp"

namespace cppmpas::phys {

namespace {
CPPMPAS_HD inline real fmax0(real x) { return x > 0.0 ? x : 0.0; }  // max(0., x): maxsd with x as the destination
// CPPWRF's lane rule (CPPWRF 7fac5bf src/core/backend.hpp:197): libmvec for the pairs of i from its, the scalar
// libm for an odd last i and a 1-trip loop
CPPMPAS_HD inline bool lane(int i, int n) { return n > 1 && !((n & 1) && i == n - 1); }
}  // namespace

void PhysicsVars::allocate_forall_physics(const std::string& microp_scheme, const std::string& pbl_scheme, int ime_,
                                          int kme_) {
  if (microp_scheme == "mp_thompson_aerosols")  // PI:131-137 (nifa_p, nwfa_p)
    throw Error("allocate_forall_physics: config_microp_scheme = 'mp_thompson_aerosols' has no verified port");
  if (pbl_scheme == "bl_mynn")  // PI:139-147 (nc_p, ni_p, nifa_p, nwfa_p)
    throw Error("allocate_forall_physics: config_pbl_scheme = 'bl_mynn' has no verified port");
  ime = ime_;
  kme = kme_;
  for (PArray* a : {&psfc_p, &ptop_p, &psfc_hyd_p, &psfc_hydd_p}) a->allocate(ime, 1);  // PI:100-101, 150-151
  for (PArray* a : {&u_p, &v_p, &fzm_p, &fzp_p, &zz_p, &pres_p, &pi_p, &z_p, &zmid_p, &dz_p, &t_p, &th_p, &al_p,
                    &rho_p, &rh_p, &znu_p, &w_p, &pres2_p, &t2_p, &qv_p, &qc_p, &qr_p, &qi_p, &qs_p, &qg_p,
                    &pres_hyd_p, &pres_hydd_p, &pres2_hyd_p, &pres2_hydd_p, &znu_hyd_p})  // PI:103-129, 152-156
    a->allocate(ime, kme);
}

// One kernel per loop of PI (one thread per cell i, the k loops inside the thread in MPAS's order).
void MPAS_to_physics(const ToPhysicsPools& mp, PhysicsVars& vp, int its, int ite) {
  using namespace physics_constants;
  if (its != 1) throw Error("MPAS_to_physics: its /= 1");
  if (mp.microp_scheme == "mp_thompson_aerosols")  // PI:357-375
    throw Error("MPAS_to_physics: config_microp_scheme = 'mp_thompson_aerosols' has no verified port");
  if (mp.pbl_scheme == "bl_mynn")  // PI:377-425
    throw Error("MPAS_to_physics: config_pbl_scheme = 'bl_mynn' has no verified port");
  if (mp.index_qv <= 0 || mp.index_qc <= 0 || mp.index_qr <= 0 || mp.index_qi <= 0 || mp.index_qs <= 0 ||
      mp.index_qg <= 0)  // an inactive constituent has MPAS's index -1 (scalar_index)
    throw Error("MPAS_to_physics: the scalars qv, qc, qr, qi, qs and qg are needed (PI:316-321)");
  const int kts = 1, kte = mp.nVertLevels;
  if (vp.ime < ite || vp.kme != kte + 1) throw Error("MPAS_to_physics: allocate_forall_physics was not called");
  const ToPhysicsArrays m = mp;
  const PhysicsVarsView v = vp.view();

  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {
    const auto q = [&](int idx, int k, int ii) { return m.scalars(idx, k, ii); };
    for (int k = kts; k <= kte; ++k) {  // PI:325-355
      v.qv_p(i, k) = fmax0(q(m.index_qv, k, i));
      v.qc_p(i, k) = fmax0(q(m.index_qc, k, i));
      v.qr_p(i, k) = fmax0(q(m.index_qr, k, i));
      v.qi_p(i, k) = fmax0(q(m.index_qi, k, i));
      v.qs_p(i, k) = fmax0(q(m.index_qs, k, i));
      v.qg_p(i, k) = fmax0(q(m.index_qg, k, i));

      v.u_p(i, k) = m.uReconstructZonal(k, i);
      v.v_p(i, k) = m.uReconstructMeridional(k, i);

      v.zz_p(i, k) = m.zz(k, i);
      v.rho_p(i, k) = m.zz(k, i) * m.rho_zz(k, i);
      v.rho_p(i, k) = v.rho_p(i, k) * (1. + v.qv_p(i, k));
      v.th_p(i, k) = m.theta_m(k, i) / (1. + R_v / R_d * v.qv_p(i, k));
      v.t_p(i, k) = v.th_p(i, k) * m.exner(k, i);

      v.pi_p(i, k) = m.exner(k, i);
      v.pres_p(i, k) = m.pressure_p(k, i) + m.pressure_base(k, i);

      v.zmid_p(i, k) = 0.5 * (m.zgrid(k + 1, i) + m.zgrid(k, i));
      v.dz_p(i, k) = m.zgrid(k + 1, i) - m.zgrid(k, i);
    }
  });

  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {  // PI:428-439
    const real tem1 = m.zgrid(2, i) - m.zgrid(1, i);
    const real tem2 = m.zgrid(3, i) - m.zgrid(2, i);
    const real rho1 = m.rho_zz(1, i) * m.zz(1, i) * (1. + v.qv_p(i, 1));
    const real rho2 = m.rho_zz(2, i) * m.zz(2, i) * (1. + v.qv_p(i, 2));
    m.surface_pressure(i) = 0.5 * gravity * (m.zgrid(2, i) - m.zgrid(1, i)) * (rho1 - 0.5 * (rho2 - rho1) * tem1 / (tem1 + tem2));
    m.surface_pressure(i) = m.surface_pressure(i) + m.pressure_p(1, i) + m.pressure_base(1, i);
  });
  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {  // PI:441-445
    for (int k = kts; k <= kte; ++k) v.znu_p(i, k) = v.pres_p(i, k) / m.surface_pressure(i);
  });

  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {  // PI:449-456
    for (int k = kts; k <= kte + 1; ++k) {
      v.w_p(i, k) = m.w(k, i);
      v.z_p(i, k) = m.zgrid(k, i);
    }
  });

  {  // PI:460-476: the log only, in cell order (the kernel marks the cells)
    WorkArrayInt low(static_cast<std::size_t>(ite), 0);
    int* lowp = low.data();
    parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {
      if (v.pres_p(i, 1) <= v.pres_p(i, 2)) lowp[i - 1] = 1;
    });
    for (int i = its; i <= ite; ++i)
      if (low[i - 1]) std::fprintf(stderr, "--- subroutine MPAS_to_phys - pressure(1) < pressure(2): i = %d\n", i);
  }

  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {  // PI:479-489
    for (int k = kts + 1; k <= kte; ++k) {
      const real tem1 = 1. / (m.zgrid(k + 1, i) - m.zgrid(k - 1, i));
      v.fzm_p(i, k) = (m.zgrid(k, i) - m.zgrid(k - 1, i)) * tem1;
      v.fzp_p(i, k) = (m.zgrid(k + 1, i) - m.zgrid(k, i)) * tem1;
      v.t2_p(i, k) = v.fzm_p(i, k) * v.t_p(i, k) + v.fzp_p(i, k) * v.t_p(i, k - 1);
      v.pres2_p(i, k) = v.fzm_p(i, k) * v.pres_p(i, k) + v.fzp_p(i, k) * v.pres_p(i, k - 1);
    }
  });

  {  // PI:494-506, k = kte+1; the i loop is vectorised (libmvec log / exp lanes)
    const int n = ite - its + 1;
    parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {
      const int k = kte + 1;
      const real z0 = m.zgrid(k, i);
      const real z1 = 0.5 * (m.zgrid(k, i) + m.zgrid(k - 1, i));
      const real z2 = 0.5 * (m.zgrid(k - 1, i) + m.zgrid(k - 2, i));
      const real w1 = (z0 - z2) / (z1 - z2);
      const real w2 = 1. - w1;
      v.t2_p(i, k) = w1 * v.t_p(i, k - 1) + w2 * v.t_p(i, k - 2);
      if (lane(i - its, n))
        v.pres2_p(i, k) = libm::mvec_exp(w1 * libm::mvec_log(v.pres_p(i, k - 1)) + w2 * libm::mvec_log(v.pres_p(i, k - 2)));
      else
        v.pres2_p(i, k) = libm::exp(w1 * libm::log(v.pres_p(i, k - 1)) + w2 * libm::log(v.pres_p(i, k - 2)));
    });
  }

  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {  // PI:510-523, k = kts
    const int k = kts;
    const real z0 = m.zgrid(k, i);
    const real z1 = 0.5 * (m.zgrid(k, i) + m.zgrid(k + 1, i));
    const real z2 = 0.5 * (m.zgrid(k + 1, i) + m.zgrid(k + 2, i));
    const real w1 = (z0 - z2) / (z1 - z2);
    const real w2 = 1. - w1;
    v.t2_p(i, k) = w1 * v.t_p(i, k) + w2 * v.t_p(i, k + 1);
    v.pres2_p(i, k) = w1 * v.pres_p(i, k) + w2 * v.pres_p(i, k + 1);
    v.psfc_p(i) = m.surface_pressure(i);
  });

  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {  // PI:526-550: hydrostatic pressures
    int k = kte + 1;
    v.pres2_hyd_p(i, k) = v.pres2_p(i, k);
    v.pres2_hydd_p(i, k) = v.pres2_p(i, k);
    for (k = kte; k >= 1; --k) {
      const real rho_a = v.rho_p(i, k) / (1. + v.qv_p(i, k));
      v.pres2_hyd_p(i, k) = v.pres2_hyd_p(i, k + 1) + gravity * v.rho_p(i, k) * v.dz_p(i, k);
      v.pres2_hydd_p(i, k) = v.pres2_hydd_p(i, k + 1) + gravity * rho_a * v.dz_p(i, k);
    }
    for (k = kte; k >= 1; --k) {
      v.pres_hyd_p(i, k) = 0.5 * (v.pres2_hyd_p(i, k + 1) + v.pres2_hyd_p(i, k));
      v.pres_hydd_p(i, k) = 0.5 * (v.pres2_hydd_p(i, k + 1) + v.pres2_hydd_p(i, k));
    }
    v.psfc_hyd_p(i) = v.pres2_hyd_p(i, 1);
    v.psfc_hydd_p(i) = v.pres2_hydd_p(i, 1);
    for (k = kte; k >= 1; --k) v.znu_hyd_p(i, k) = v.pres_hyd_p(i, k) / v.psfc_hyd_p(i);
  });

  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) { m.plrad(i) = v.pres2_p(i, kte + 1); });  // PI:553-557
}

}  // namespace cppmpas::phys
