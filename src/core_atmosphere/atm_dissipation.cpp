// SPDX-License-Identifier: GPL-3.0-or-later
// src/core_atmosphere/dynamics/mpas_atm_dissipation_models.F (DM), the branches the mountain-wave reference takes.
// Reference object code (bin/atmosphere_model of the mpas-build asset 153c7206...): u_dissipation_3d (0x25bf50),
// w_dissipation_3d (0x25a1d0) and scalar_dissipation_3d_les (0x256140) have no FMA and no libm call on these paths
// (u_dissipation_3d's one sincos@plt is DM:860-861, the perturbation branch, not taken); every statement is the
// Fortran one evaluated left to right.
#include "core_atmosphere/atm_dissipation.hpp"

#include <algorithm>

#include "core/error.hpp"

namespace cppmpas {

void u_dissipation_3d(const Mesh& m, const DissipationArgs& a, const FArray<real>& eddy_visc_horz,
                      FArray<real>& delsq_u, const FArray<real>& u, const FArray<real>& divergence,
                      const FArray<real>& vorticity, const FArray<real>& rho_edge, FArray<real>& tend_u_euler) {
  const int nVertLevels = m.nVertLevels;
  const real tau_12_factor = 0.0;  // DM:686-687, les_model_opt = LES_MODEL_NONE
  for (int iEdge = 1; iEdge <= m.nEdges; ++iEdge) {  // DM:690-724
    const int cell1 = m.cellsOnEdge(1, iEdge);
    const int cell2 = m.cellsOnEdge(2, iEdge);
    const int vertex1 = m.verticesOnEdge(1, iEdge);
    const int vertex2 = m.verticesOnEdge(2, iEdge);
    const real r_dc = m.invDcEdge(iEdge);
    const real r_dv = std::min(m.invDvEdge(iEdge), 4 * m.invDcEdge(iEdge));
    for (int k = 1; k <= nVertLevels; ++k) delsq_u(k, iEdge) = 0.0;
    for (int k = 1; k <= nVertLevels; ++k) {
      const real u_diffusion = (divergence(k, cell2) - divergence(k, cell1)) * r_dc -
                               (vorticity(k, vertex2) - vorticity(k, vertex1)) * r_dv;  // DM:709-710
      const real u_diffusion_les =
          u_diffusion + tau_12_factor * (divergence(k, cell2) - divergence(k, cell1)) * r_dc;  // DM:713
      delsq_u(k, iEdge) = delsq_u(k, iEdge) + u_diffusion;
      const real kdiffu = 0.5 * (eddy_visc_horz(k, cell1) + eddy_visc_horz(k, cell2));
      tend_u_euler(k, iEdge) =
          tend_u_euler(k, iEdge) + rho_edge(k, iEdge) * kdiffu * u_diffusion_les * m.meshScalingDel2(iEdge);
    }
  }
  if (a.h_mom_eddy_visc4 > 0.0) throw Error("u_dissipation_3d: the 4th-order filter (DM:728-811) is not ported");
  if (a.v_mom_eddy_visc2 > 0.0) {  // DM:816
    if (!a.config_mix_full)
      throw Error("u_dissipation_3d: vertical mixing of the perturbation state (DM:848-884) is not ported");
    for (int iEdge = 1; iEdge <= m.nEdges; ++iEdge) {  // DM:823-844
      const int cell1 = m.cellsOnEdge(1, iEdge);
      const int cell2 = m.cellsOnEdge(2, iEdge);
      for (int k = 2; k <= nVertLevels - 1; ++k) {
        const real z1 = 0.5 * (m.zgrid(k - 1, cell1) + m.zgrid(k - 1, cell2));
        const real z2 = 0.5 * (m.zgrid(k, cell1) + m.zgrid(k, cell2));
        const real z3 = 0.5 * (m.zgrid(k + 1, cell1) + m.zgrid(k + 1, cell2));
        const real z4 = 0.5 * (m.zgrid(k + 2, cell1) + m.zgrid(k + 2, cell2));
        const real zm = 0.5 * (z1 + z2);
        const real z0 = 0.5 * (z2 + z3);
        const real zp = 0.5 * (z3 + z4);
        tend_u_euler(k, iEdge) = tend_u_euler(k, iEdge) +
                                 rho_edge(k, iEdge) * a.v_mom_eddy_visc2 *
                                     ((u(k + 1, iEdge) - u(k, iEdge)) / (zp - z0) -
                                      (u(k, iEdge) - u(k - 1, iEdge)) / (z0 - zm)) /
                                     (0.5 * (zp - zm));  // DM:840-842
      }
    }
  } else {
    throw Error("u_dissipation_3d: v_mom_eddy_visc2 <= 0 (no vertical mixing, DM:816 not taken) is not verified");
  }
}

void w_dissipation_3d(const Mesh& m, const DissipationArgs& a, FArray<real>& delsq_w, const FArray<real>& w,
                      const FArray<real>& rho_edge, const FArray<real>& rho_zz, const FArray<real>& eddy_visc_horz,
                      FArray<real>& tend_w_euler) {
  const int nVertLevels = m.nVertLevels;
  for (int iCell = 1; iCell <= m.nCells; ++iCell) {  // DM:1032-1067
    for (int k = 1; k <= nVertLevels; ++k) delsq_w(k, iCell) = 0.0;
    for (int k = 1; k <= nVertLevels + 1; ++k) tend_w_euler(k, iCell) = 0.0;
    const real r_areaCell = m.invAreaCell(iCell);
    for (int i = 1; i <= m.nEdgesOnCell(iCell); ++i) {
      const int iEdge = m.edgesOnCell(i, iCell);
      const real edge_sign = 0.5 * r_areaCell * m.edgesOnCell_sign(i, iCell) * m.dvEdge(iEdge) * m.invDcEdge(iEdge);
      const int cell1 = m.cellsOnEdge(1, iEdge);
      const int cell2 = m.cellsOnEdge(2, iEdge);
      for (int k = 2; k <= nVertLevels; ++k) {
        real w_turb_flux = edge_sign * (rho_edge(k, iEdge) + rho_edge(k - 1, iEdge)) * (w(k, cell2) - w(k, cell1));
        delsq_w(k, iCell) = delsq_w(k, iCell) + w_turb_flux;
        w_turb_flux = w_turb_flux * m.meshScalingDel2(iEdge) * 0.25 *
                      (eddy_visc_horz(k, cell1) + eddy_visc_horz(k, cell2) + eddy_visc_horz(k - 1, cell1) +
                       eddy_visc_horz(k - 1, cell2));  // DM:1061-1063
        tend_w_euler(k, iCell) = tend_w_euler(k, iCell) + w_turb_flux;
      }
    }
  }
  if (a.h_mom_eddy_visc4 > 0.0) throw Error("w_dissipation_3d: the 4th-order filter (DM:1073-1100) is not ported");
  if (a.v_mom_eddy_visc2 > 0.0) {  // DM:1102-1119
    for (int iCell = 1; iCell <= m.nCells; ++iCell) {
      for (int k = 2; k <= nVertLevels; ++k) {
        tend_w_euler(k, iCell) = tend_w_euler(k, iCell) + a.v_mom_eddy_visc2 * 0.5 *
                                                              (rho_zz(k, iCell) + rho_zz(k - 1, iCell)) *
                                                              ((w(k + 1, iCell) - w(k, iCell)) * m.rdzw(k) -
                                                               (w(k, iCell) - w(k - 1, iCell)) * m.rdzw(k - 1)) *
                                                              m.rdzu(k);  // DM:1111-1113
      }
    }
  } else {
    throw Error("w_dissipation_3d: v_mom_eddy_visc2 <= 0 (no vertical mixing, DM:1102 not taken) is not verified");
  }
}

void scalar_dissipation_3d_les(const Mesh& m, const DissipationArgs& a, FArray<real>& delsq_theta,
                               const FArray<real>& theta_m, const FArray<real>& rho_edge, const FArray<real>& rho_zz,
                               const FArray<real>& eddy_visc_horz, FArray<real>& tend_theta_euler) {
  const int nVertLevels = m.nVertLevels;
  for (int iCell = 1; iCell <= m.nCells; ++iCell) {  // DM:1268-1300
    for (int k = 1; k <= nVertLevels; ++k) {
      delsq_theta(k, iCell) = 0.0;
      tend_theta_euler(k, iCell) = 0.0;
    }
    const real r_areaCell = m.invAreaCell(iCell);
    for (int i = 1; i <= m.nEdgesOnCell(iCell); ++i) {
      const int iEdge = m.edgesOnCell(i, iCell);
      const real edge_sign = r_areaCell * m.edgesOnCell_sign(i, iCell) * m.dvEdge(iEdge) * m.invDcEdge(iEdge);
      const real pr_scale = a.prandtl_inv * m.meshScalingDel2(iEdge);
      const int cell1 = m.cellsOnEdge(1, iEdge);
      const int cell2 = m.cellsOnEdge(2, iEdge);
      for (int k = 1; k <= nVertLevels; ++k) {
        real theta_turb_flux = edge_sign * (theta_m(k, cell2) - theta_m(k, cell1)) * rho_edge(k, iEdge);
        delsq_theta(k, iCell) = delsq_theta(k, iCell) + theta_turb_flux;
        theta_turb_flux =
            theta_turb_flux * 0.5 * (eddy_visc_horz(k, cell1) + eddy_visc_horz(k, cell2)) * pr_scale;  // DM:1295
        tend_theta_euler(k, iCell) = tend_theta_euler(k, iCell) + theta_turb_flux;
      }
    }
  }
  if (a.h_theta_eddy_visc4 > 0.0)
    throw Error("scalar_dissipation_3d_les: the 4th-order filter (DM:1306-1332) is not ported");
  if (a.v_theta_eddy_visc2 > 0.0) {  // DM:1414
    if (!a.config_mix_full)
      throw Error("scalar_dissipation_3d_les: mixing of the perturbation from t_init (DM:1442-1466) is not ported");
    for (int iCell = 1; iCell <= m.nCells; ++iCell) {  // DM:1421-1438
      for (int k = 2; k <= nVertLevels - 1; ++k) {
        const real z1 = m.zgrid(k - 1, iCell);
        const real z2 = m.zgrid(k, iCell);
        const real z3 = m.zgrid(k + 1, iCell);
        const real z4 = m.zgrid(k + 2, iCell);
        const real zm = 0.5 * (z1 + z2);
        const real z0 = 0.5 * (z2 + z3);
        const real zp = 0.5 * (z3 + z4);
        tend_theta_euler(k, iCell) = tend_theta_euler(k, iCell) +
                                     a.v_theta_eddy_visc2 * a.prandtl_inv * rho_zz(k, iCell) *
                                         ((theta_m(k + 1, iCell) - theta_m(k, iCell)) / (zp - z0) -
                                          (theta_m(k, iCell) - theta_m(k - 1, iCell)) / (z0 - zm)) /
                                         (0.5 * (zp - zm));  // DM:1434-1436
      }
    }
  } else {
    throw Error("scalar_dissipation_3d_les: v_theta_eddy_visc2 <= 0 (DM:1414 not taken) is not verified");
  }
}

}  // namespace cppmpas
