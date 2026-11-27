// SPDX-License-Identifier: GPL-3.0-or-later
// The vertical mixing of the perturbation from the initial 1-D state (config_mix_full = .false.) of
// src/core_atmosphere/dynamics/mpas_atm_dissipation_models.F (DM), as the supercell reference takes it:
// u_dissipation_3d DM:850-882 and scalar_dissipation_3d_les DM:1442-1466 (called from atm_dissipation.cpp).
// Reference object code (bin/atmosphere_model of the mpas-build asset 153c7206...): u_dissipation_3d (0x25bf50)
// calls sincos@plt once per edge (0x25e5b9, angleEdge(iEdge) hoisted out of the k loop; cos in 0x260(%rsp), sin in
// 0x268(%rsp)) and evaluates u_mix(k) = (u(k) - u_init(k)*cos) + v_init(k)*sin (0x25e640-0x25e661: mulpd, subpd,
// mulpd, addpd); no FMA, no other libm call on these branches.
// Both loops are kernels (core/backend.hpp), one thread per edge / cell with the k loops inside; u_mix is the
// thread's own column.
#include "core/error.hpp"
#include "core/libm.hpp"
#include "core/views.hpp"
#include "core_atmosphere/atm_dissipation.hpp"

namespace cppmpas {

void u_vertical_mixing_perturbation(const Mesh& m, const DissipationArgs& a, const FArray<real>& u,
                                    const FArray<real>& rho_edge, FArray<real>& tend_u_euler) {
  if (a.state_init == nullptr) throw Error("u_dissipation_3d: no u_init / v_init for the perturbation mixing");
  const int nVertLevels = m.nVertLevels;
  fits_column(nVertLevels, "u_dissipation_3d");
  const auto u_init = a.state_init->u_init.cview(), v_init = a.state_init->v_init.cview();
  parallel_for(1, m.nEdgesSolve(), [=, m = view(m), u = view(u), rho_edge = view(rho_edge),
                                    tend_u_euler = view(tend_u_euler)] CPPMPAS_LAMBDA(int iEdge) {  // DM:853-880
    Column<real> u_mix(nVertLevels + 1);
    const int cell1 = m.cellsOnEdge(1, iEdge);
    const int cell2 = m.cellsOnEdge(2, iEdge);
    real sina = 0, cosa = 0;
    libm::sincos(m.angleEdge(iEdge), &sina, &cosa);
    for (int k = 1; k <= nVertLevels; ++k)  // DM:859-862
      u_mix[k] = u(k, iEdge) - u_init(k) * cosa + v_init(k) * sina;
    for (int k = 2; k <= nVertLevels - 1; ++k) {  // DM:864-879
      const real z1 = 0.5 * (m.zgrid(k - 1, cell1) + m.zgrid(k - 1, cell2));
      const real z2 = 0.5 * (m.zgrid(k, cell1) + m.zgrid(k, cell2));
      const real z3 = 0.5 * (m.zgrid(k + 1, cell1) + m.zgrid(k + 1, cell2));
      const real z4 = 0.5 * (m.zgrid(k + 2, cell1) + m.zgrid(k + 2, cell2));
      const real zm = 0.5 * (z1 + z2);
      const real z0 = 0.5 * (z2 + z3);
      const real zp = 0.5 * (z3 + z4);
      tend_u_euler(k, iEdge) = tend_u_euler(k, iEdge) + rho_edge(k, iEdge) * a.v_mom_eddy_visc2 *
                                                            ((u_mix[k + 1] - u_mix[k]) / (zp - z0) -
                                                             (u_mix[k] - u_mix[k - 1]) / (z0 - zm)) /
                                                            (0.5 * (zp - zm));  // DM:876-878
    }
  });
}

void theta_vertical_mixing_perturbation(const Mesh& m, const DissipationArgs& a, const FArray<real>& theta_m,
                                        const FArray<real>& rho_zz, FArray<real>& tend_theta_euler) {
  if (a.state_init == nullptr) throw Error("scalar_dissipation_3d_les: no t_init for the perturbation mixing");
  const int nVertLevels = m.nVertLevels;
  const auto t_init = a.state_init->t_init.cview();
  parallel_for(1, m.nCellsSolve(), [=, m = view(m), theta_m = view(theta_m), rho_zz = view(rho_zz),
                                    tend_theta_euler = view(tend_theta_euler)] CPPMPAS_LAMBDA(int iCell) {  // DM:1447-1464
    for (int k = 2; k <= nVertLevels - 1; ++k) {
      const real z1 = m.zgrid(k - 1, iCell);
      const real z2 = m.zgrid(k, iCell);
      const real z3 = m.zgrid(k + 1, iCell);
      const real z4 = m.zgrid(k + 2, iCell);
      const real zm = 0.5 * (z1 + z2);
      const real z0 = 0.5 * (z2 + z3);
      const real zp = 0.5 * (z3 + z4);
      tend_theta_euler(k, iCell) =
          tend_theta_euler(k, iCell) +
          a.v_theta_eddy_visc2 * a.prandtl_inv * rho_zz(k, iCell) *
              (((theta_m(k + 1, iCell) - t_init(k + 1, iCell)) - (theta_m(k, iCell) - t_init(k, iCell))) / (zp - z0) -
               ((theta_m(k, iCell) - t_init(k, iCell)) - (theta_m(k - 1, iCell) - t_init(k - 1, iCell))) / (z0 - zm)) /
              (0.5 * (zp - zm));  // DM:1460-1462
    }
  });
}

}  // namespace cppmpas
