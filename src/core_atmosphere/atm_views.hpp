// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// Kernel views (core/views.hpp) of the dycore's fields (atm_dycore.hpp DycoreFields) and of the solve-diagnostics
// work arrays (atm_solve_diagnostics.hpp SolveDiagWork), under the members' own names.
#include "core/views.hpp"
#include "core_atmosphere/atm_dycore.hpp"
#include "core_atmosphere/atm_solve_diagnostics.hpp"

namespace cppmpas {

template <bool C> struct SolveDiagWorkViewT {
  FView<VElem<C, real>> ke_vertex, ke_edge;
};
template <class W> auto view_sdw(W& w) {
  SolveDiagWorkViewT<std::is_const_v<W>> v;
  v.ke_vertex = w.ke_vertex.view();
  v.ke_edge = w.ke_edge.view();
  return v;
}
inline SolveDiagWorkViewT<true> view(const SolveDiagWork& w) { return view_sdw(w); }
inline SolveDiagWorkViewT<false> view(SolveDiagWork& w) { return view_sdw(w); }

#define CPPMPAS_DYCORE_REAL_ARRAYS(X)                                                                                 \
  X(u2) X(w2) X(theta_m2) X(rho_zz2) X(scalars2) X(ru_save) X(rw_save) X(rtheta_p_save) X(rho_p_save)                \
  X(rho_zz_old_split) X(cqw) X(cqu) X(cofrz) X(cofwr) X(cofwz) X(coftz) X(cofwt) X(a_tri) X(alpha_tri) X(gamma_tri)  \
  X(ru_p) X(rw_p) X(rtheta_pp) X(rtheta_pp_old) X(rho_pp) X(ruAvg) X(wwAvg) X(ruAvg_split) X(wwAvg_split)            \
  X(h_divergence) X(eddy_visc_horz) X(tend_u) X(tend_w) X(tend_theta) X(tend_rho) X(tend_u_euler) X(tend_w_euler)    \
  X(tend_theta_euler) X(rt_diabatic_tend) X(tend_scalars) X(rthdynten) X(tend_ru_physics) X(tend_rtheta_physics)     \
  X(tend_rho_physics) X(qtot) X(dpdz) X(delsq_theta) X(delsq_w) X(delsq_divergence) X(delsq_u) X(delsq_vorticity)

template <bool C> struct DycoreFieldsViewT {
#define CPPMPAS_X(n) FView<VElem<C, real>> n;
  CPPMPAS_DYCORE_REAL_ARRAYS(CPPMPAS_X)
#undef CPPMPAS_X
  real Time1, Time2;
  SolveDiagWorkViewT<C> sdw;
  int moist_start, moist_end, index_qv;
};
template <class F> auto view_dycore(F& f) {
  DycoreFieldsViewT<std::is_const_v<F>> v;
#define CPPMPAS_X(n) v.n = f.n.view();
  CPPMPAS_DYCORE_REAL_ARRAYS(CPPMPAS_X)
#undef CPPMPAS_X
  v.Time1 = f.Time1;
  v.Time2 = f.Time2;
  v.sdw = view(f.sdw);
  v.moist_start = f.moist_start;
  v.moist_end = f.moist_end;
  v.index_qv = f.index_qv;
  return v;
}
inline DycoreFieldsViewT<true> view(const DycoreFields& f) { return view_dycore(f); }
inline DycoreFieldsViewT<false> view(DycoreFields& f) { return view_dycore(f); }

// the numeric options of DycoreConfig (its strings are host-only), under their names
#define CPPMPAS_DYCORE_CONFIG_NUMBERS(X)                                                                              \
  X(int, config_time_integration_order) X(real, config_dt) X(bool, config_split_dynamics_transport)                  \
  X(int, config_number_of_sub_steps) X(int, config_dynamics_split_steps) X(real, config_h_mom_eddy_visc2)            \
  X(real, config_h_mom_eddy_visc4) X(real, config_v_mom_eddy_visc2) X(real, config_h_theta_eddy_visc2)               \
  X(real, config_h_theta_eddy_visc4) X(real, config_v_theta_eddy_visc2) X(real, config_len_disp)                     \
  X(bool, config_mix_scalars) X(real, config_del4u_div_factor) X(bool, config_scalar_advection)                      \
  X(bool, config_positive_definite) X(bool, config_monotonic) X(real, config_coef_3rd_order) X(bool, config_mix_full) \
  X(real, config_smdiv) X(real, config_apvm_upwinding) X(real, config_mpas_cam_coef) X(bool, config_rayleigh_damp_u)  \
  X(bool, config_apply_lbcs)
struct DycoreConfigView {
#define CPPMPAS_X(t, n) t n;
  CPPMPAS_DYCORE_CONFIG_NUMBERS(CPPMPAS_X)
#undef CPPMPAS_X
};
inline DycoreConfigView view(const DycoreConfig& c) {
  DycoreConfigView v;
#define CPPMPAS_X(t, n) v.n = c.n;
  CPPMPAS_DYCORE_CONFIG_NUMBERS(CPPMPAS_X)
#undef CPPMPAS_X
  return v;
}

}  // namespace cppmpas
