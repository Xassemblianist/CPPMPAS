// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// The dissipation routines of src/core_atmosphere/dynamics/mpas_atm_dissipation_models.F (DM) that
// atm_compute_dyn_tend calls at rk_step 1, on the path of the mountain-wave reference: les_model_opt = none,
// 2nd-order horizontal mixing with eddy_visc_horz, no 4th-order filter (h_mom_eddy_visc4 = h_theta_eddy_visc4 = 0),
// vertical mixing of the full state (config_mix_full = .true.) with v_*_eddy_visc2 > 0, no scalar mixing. Any
// other branch is an Error (atm_dycore.hpp check_dycore_supported refuses the configurations first).
#include "core/mesh.hpp"

namespace cppmpas {

struct DissipationArgs {
  bool config_mix_full = true;
  real h_mom_eddy_visc4 = 0, v_mom_eddy_visc2 = 0;
  real h_theta_eddy_visc4 = 0, v_theta_eddy_visc2 = 0, prandtl_inv = 1;
  const State* state_init = nullptr;  // u_init, v_init, t_init of the perturbation mixing (config_mix_full = .false.)
};

// DM:579-947. u, divergence, vorticity, rho_edge of the current solve diagnostics; tend_u_euler accumulates.
void u_dissipation_3d(const Mesh& m, const DissipationArgs& a, const FArray<real>& eddy_visc_horz,
                      FArray<real>& delsq_u, const FArray<real>& u, const FArray<real>& divergence,
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
