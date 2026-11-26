// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// MPAS-Atmosphere's interface from the dynamics to the physics (physics/mpas_atmphys_interface.F, PI): the arrays
// allocate_forall_physics allocates in mpas_atmphys_vars (PI:86-158) and MPAS_to_physics (PI:235-559), statement by
// statement. One rank, one thread (the reference runs 1 MPI rank, OMP_NUM_THREADS = 1): its..ite = 1..nCellsSolve,
// ims = its, ime = ite, kts..kte = 1..nVertLevels, kme = kte + 1, jts = jte = 1 (mpas_atmphys_manager.F:661-669).
//
// Reference object code (bin/atmosphere_model, sha256 524a1cf4...): __mpas_atmphys_interface_MOD_mpas_to_physics
// has no FMA and no ymm register; its libm calls are the top-of-the-model pressure of PI:504,
// exp(w1*log(pres_p(k-1)) + w2*log(pres_p(k-2))), vectorised over i: _ZGVbN2v_log x2 and _ZGVbN2v_exp, with the
// scalar log@plt / exp@plt copies for an odd last i and a 1-trip loop (CPPWRF's lane(i, n) rule, pairs of i from its).
// CPU vendor: _ZGVbN2v_log rounds an rcpps result (Intel's table in libm.hpp mvec_log; the AMD reference runners'
// may differ); pres2_p(kte+1) seeds every hydrostatic pressure (pres2_hyd_p, pres_hyd_p) the surface layer and YSU
// read. _ZGVbN2v_exp is table-driven and vendor-independent.
// MPAS's file compiled standalone with the reference flags has the reference's arithmetic instruction counts (only the
// set-up of the stub pool calls differs).
#include <cstddef>
#include <string>
#include <vector>

#include "core/field.hpp"
#include "core/precision.hpp"

namespace cppmpas::phys {

// A (ims:ime, kms:kme, jms:jme) or (ims:ime, jms:jme) array of mpas_atmphys_vars, 1-based, jms = jme = 1.
struct PArray {
  std::vector<real> a;
  int ldi = 0;
  void allocate(int ni, int nk) {
    ldi = ni;
    a.assign(static_cast<std::size_t>(ni) * nk, real(0));
  }
  real& operator()(int i, int k = 1) { return a[static_cast<std::size_t>(k - 1) * ldi + (i - 1)]; }
  const real& operator()(int i, int k = 1) const { return a[static_cast<std::size_t>(k - 1) * ldi + (i - 1)]; }
  real* data() { return a.data(); }
  const real* data() const { return a.data(); }
};

// mpas_atmphys_vars' arrays shared by all physics drivers (allocate_forall_physics, PI:100-156)
struct PhysicsVars {
  int ime = 0, kme = 0;  // ims = 1, kms = 1
  PArray psfc_p, ptop_p;
  PArray u_p, v_p, fzm_p, fzp_p, zz_p, pres_p, pi_p, z_p, zmid_p, dz_p, t_p, th_p, al_p, rho_p, rh_p, znu_p;
  PArray w_p, pres2_p, t2_p;
  PArray qv_p, qc_p, qr_p, qi_p, qs_p, qg_p;
  PArray psfc_hyd_p, psfc_hydd_p, pres_hyd_p, pres_hydd_p, pres2_hyd_p, pres2_hydd_p, znu_hyd_p;
  // PI:86-158 for config_microp_scheme / config_pbl_scheme (mp_thompson_aerosols and bl_mynn are refused)
  void allocate_forall_physics(const std::string& microp_scheme, const std::string& pbl_scheme, int ime_, int kme_);
};

// The pool arrays MPAS_to_physics reads or writes (Registry names; Fortran layout (k, iCell), 1-based).
struct ToPhysicsPools {
  int nVertLevels = 0;
  std::string microp_scheme, pbl_scheme;
  // mesh
  FView<const real> zgrid, zz;
  // state (time_lev)
  FView<const real> rho_zz, theta_m, w, scalars;
  int index_qv = 0, index_qc = 0, index_qr = 0, index_qi = 0, index_qs = 0, index_qg = 0;
  // diag
  FView<real> surface_pressure;
  FView<const real> exner, pressure_base, pressure_p, uReconstructZonal, uReconstructMeridional;
  // diag_physics
  FView<real> plrad;
};

// PI:235-559 on cells its..ite (its = 1)
void MPAS_to_physics(const ToPhysicsPools& m, PhysicsVars& v, int its, int ite);

}  // namespace cppmpas::phys
