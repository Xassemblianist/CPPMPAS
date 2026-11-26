// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// The dissipation routines of src/core_atmosphere/dynamics/mpas_atm_dissipation_models.F (DM) that
// atm_compute_dyn_tend calls at rk_step 1, on the paths of the mountain-wave (MW) and JW baroclinic-wave references:
// les_model_opt = none; the eddy viscosity of smagorinsky_2d (JW) or 2d_fixed (MW, set in atm_compute_dyn_tend);
// 2nd-order horizontal mixing with eddy_visc_horz; the 4th-order filters when h_mom_eddy_visc4 / h_theta_eddy_visc4
// > 0 (JW, from smagorinsky_2d); vertical mixing of the full state (config_mix_full = .true.) when v_*_eddy_visc2 > 0
// (MW), none when v_*_eddy_visc2 = 0 (JW); no scalar mixing. Any other branch is an Error (atm_dycore.hpp
// check_dycore_supported refuses the configurations first).
#include "core/mesh.hpp"

namespace cppmpas {

struct DissipationArgs {
  bool config_mix_full = true;
  real config_del4u_div_factor = 10;  // Registry.xml:189
  real h_mom_eddy_visc4 = 0, v_mom_eddy_visc2 = 0;
  real h_theta_eddy_visc4 = 0, v_theta_eddy_visc2 = 0, prandtl_inv = 1;
  const State* state_init = nullptr;  // u_init, v_init, t_init of the perturbation mixing (config_mix_full = .false.)
};

// DM:119-204, called at TI:6348-6352 for config_horiz_mixing = '2d_smagorinsky': kdiff (eddy_visc_horz) of every
// cell from u and v of the current state and diagnostics, h_mom_eddy_visc4 = h_theta_eddy_visc4 =
// config_visc4_2dsmag * config_len_disp**3.
void smagorinsky_2d(const Mesh& m, FArray<real>& kdiff, const FArray<real>& u, const FArray<real>& v, real c_s,
                    real config_len_disp, real invDt, real config_visc4_2dsmag, real& h_mom_eddy_visc4,
                    real& h_theta_eddy_visc4);
// DM:579-947. u, divergence, vorticity, rho_edge of the current solve diagnostics; tend_u_euler accumulates.
void u_dissipation_3d(const Mesh& m, const DissipationArgs& a, const FArray<real>& eddy_visc_horz,
                      FArray<real>& delsq_u, FArray<real>& delsq_vorticity, FArray<real>& delsq_divergence,
                      const FArray<real>& u, const FArray<real>& divergence,
                      const FArray<real>& vorticity, const FArray<real>& rho_edge, FArray<real>& tend_u_euler);
// DM:951-1153. tend_w_euler is zeroed and set.
void w_dissipation_3d(const Mesh& m, const DissipationArgs& a, FArray<real>& delsq_w, const FArray<real>& w,
                      const FArray<real>& rho_edge, const FArray<real>& rho_zz, const FArray<real>& eddy_visc_horz,
                      FArray<real>& tend_w_euler);
// DM:1157-1606. tend_theta_euler is zeroed and set.
void scalar_dissipation_3d_les(const Mesh& m, const DissipationArgs& a, FArray<real>& delsq_theta,
                               const FArray<real>& theta_m, const FArray<real>& rho_edge, const FArray<real>& rho_zz,
                               const FArray<real>& eddy_visc_horz, FArray<real>& tend_theta_euler);

// The perturbation branches of the vertical mixing (config_mix_full = .false., the supercell reference), in
// atm_dissipation_perturbation.cpp: DM:850-882 (u_dissipation_3d) and DM:1442-1466 (scalar_dissipation_3d_les).
void u_vertical_mixing_perturbation(const Mesh& m, const DissipationArgs& a, const FArray<real>& u,
                                    const FArray<real>& rho_edge, FArray<real>& tend_u_euler);
void theta_vertical_mixing_perturbation(const Mesh& m, const DissipationArgs& a, const FArray<real>& theta_m,
                                        const FArray<real>& rho_zz, FArray<real>& tend_theta_euler);

}  // namespace cppmpas
