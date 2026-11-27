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

#include "core/backend.hpp"
#include "core/error.hpp"
#include "core/field.hpp"
#include "core/precision.hpp"

namespace cppmpas::phys {

// Kernel view of a PArray: trivially copyable, CPPMPAS_HD access with the PArray's 1-based (i, k) indexing.
struct PCView {
  const real* p = nullptr;
  int ldi = 0;
  CPPMPAS_HD const real& operator()(int i, int k = 1) const {
    return p[static_cast<std::size_t>(k - 1) * ldi + (i - 1)];
  }
};
struct PView {
  real* p = nullptr;
  int ldi = 0;
  CPPMPAS_HD real& operator()(int i, int k = 1) const { return p[static_cast<std::size_t>(k - 1) * ldi + (i - 1)]; }
  CPPMPAS_HD PCView cview() const { return PCView{p, ldi}; }
};

// A (ims:ime, kms:kme, jms:jme) or (ims:ime, jms:jme) array of mpas_atmphys_vars, 1-based, jms = jme = 1.
// Storage: core/field.hpp FStorage (host memory in the CPU build, CUDA managed memory in the CUDA build), so the
// physics kernels and any host code read the same array.
struct PArray {
  FStorage<real> a;
  int ldi = 0;
  void allocate(int ni, int nk) {
    ldi = ni;
    a.assign(static_cast<std::size_t>(ni) * nk, real(0));
  }
  real& operator()(int i, int k = 1) { return a[static_cast<std::size_t>(k - 1) * ldi + (i - 1)]; }
  const real& operator()(int i, int k = 1) const { return a[static_cast<std::size_t>(k - 1) * ldi + (i - 1)]; }
  real* data() { return a.data(); }
  const real* data() const { return a.data(); }
  PView view() { return PView{a.data(), ldi}; }
  PCView view() const { return PCView{a.data(), ldi}; }
  PCView cview() const { return PCView{a.data(), ldi}; }
};

// the arrays of PhysicsVars (one list: the struct, its kernel view and allocate_forall_physics)
#define CPPMPAS_PHYSICS_VARS(X)                                                                                  \
  X(psfc_p) X(ptop_p) X(u_p) X(v_p) X(fzm_p) X(fzp_p) X(zz_p) X(pres_p) X(pi_p) X(z_p) X(zmid_p) X(dz_p) X(t_p)   \
  X(th_p) X(al_p) X(rho_p) X(rh_p) X(znu_p) X(w_p) X(pres2_p) X(t2_p) X(qv_p) X(qc_p) X(qr_p) X(qi_p) X(qs_p)     \
  X(qg_p) X(psfc_hyd_p) X(psfc_hydd_p) X(pres_hyd_p) X(pres_hydd_p) X(pres2_hyd_p) X(pres2_hydd_p) X(znu_hyd_p)

// Kernel view of PhysicsVars: the same member names, as PViews
struct PhysicsVarsView {
  int ime = 0, kme = 0;
#define CPPMPAS_X(a) PView a;
  CPPMPAS_PHYSICS_VARS(CPPMPAS_X)
#undef CPPMPAS_X
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
  // the kernel view (the arrays stay where they are; a const PhysicsVars is read through it by convention)
  PhysicsVarsView view() const {
    PhysicsVarsView w;
    w.ime = ime;
    w.kme = kme;
#define CPPMPAS_X(a) w.a = PView{const_cast<real*>(a.data()), a.ldi};
    CPPMPAS_PHYSICS_VARS(CPPMPAS_X)
#undef CPPMPAS_X
    return w;
  }
};

// A trivially copyable object (a scheme's tables) where kernels can read it: in managed memory in the CUDA build, in
// host memory in the CPU build. get() is the address a kernel reads.
template <class T> struct KernelCopy {
  FStorage<T> s;
  void set(const T& x) { s.assign(1, x); }
  const T* get() const {
    if (s.empty()) throw Error("KernelCopy: the object was not set");
    return s.data();
  }
};

// The address of x where kernels can read it: a managed copy in keep (CUDA build), x itself (CPU build)
template <class T> const T* kernel_ptr(const T& x, KernelCopy<T>& keep) {
#if defined(CPPMPAS_USE_CUDA)
  keep.set(x);
  return keep.get();
#else
  (void)keep;
  return &x;
#endif
}

// A 0-based work array of a driver (a (ims:ime) std::vector of MPAS's allocate_* before) where kernels can write it
using WorkArray = FStorage<real>;
using WorkArrayInt = FStorage<int>;

// The pool arrays MPAS_to_physics reads or writes (Registry names; Fortran layout (k, iCell), 1-based).
// ToPhysicsArrays is the part a kernel captures (trivially copyable); ToPhysicsPools adds the configs.
struct ToPhysicsArrays {
  int nVertLevels = 0;
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
struct ToPhysicsPools : ToPhysicsArrays {
  std::string microp_scheme, pbl_scheme;
};

// PI:235-559 on cells its..ite (its = 1)
void MPAS_to_physics(const ToPhysicsPools& m, PhysicsVars& v, int its, int ite);

}  // namespace cppmpas::phys
