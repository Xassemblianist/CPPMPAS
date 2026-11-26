// SPDX-License-Identifier: GPL-3.0-or-later
// src/core_atmosphere/dynamics/mpas_atm_dissipation_models.F (DM), the branches the mountain-wave and JW references
// take. Reference object code (bin/atmosphere_model of the mpas-build asset 153c7206...): u_dissipation_3d
// (0x25bf50-0x25eff0), w_dissipation_3d (0x25a1d0) and scalar_dissipation_3d_les (0x256140) have no FMA, no ymm
// register and no libm call on these paths, the 4th-order filters included (u_dissipation_3d's one sincos@plt is
// DM:860-861, the perturbation branch, not taken); every statement is the Fortran one evaluated left to right.
// smagorinsky_2d (0x263480-0x264730) calls only malloc / memset / free; its k loop of DM:183-191 is 2-lane SSE2
// (sqrtpd, then minpd with the bound as the source operand, 0x264447-0x26444f) with a scalar epilogue for an odd
// nVertLevels (sqrtsd / minsd, 0x264507-0x26450f): the same per-element operations.
#include "core_atmosphere/atm_dissipation.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

#include "core/error.hpp"

namespace cppmpas {

// DM:119-204 (cellStart..cellEnd = every cell on one thread)
void smagorinsky_2d(const Mesh& m, FArray<real>& kdiff, const FArray<real>& u, const FArray<real>& v, real c_s,
                    real config_len_disp, real invDt, real config_visc4_2dsmag, real& h_mom_eddy_visc4,
                    real& h_theta_eddy_visc4) {
  const int nVertLevels = m.nVertLevels;
  std::vector<real> dudx(nVertLevels + 1), dudy(nVertLevels + 1), dvdx(nVertLevels + 1), dvdy(nVertLevels + 1);
  // (c_s*config_len_disp)**2 and (0.01*config_len_disp**2)*invDt: integer powers as multiplies (0x26437a-0x26438e;
  // the reference forms (len*len)*0.01, the same product)
  const real cl = c_s * config_len_disp;
  const real kfac = cl * cl;
  const real kmax = (0.01 * (config_len_disp * config_len_disp)) * invDt;
  for (int iCell = 1; iCell <= m.nCells; ++iCell) {
    for (int k = 1; k <= nVertLevels; ++k) {  // DM:159-164
      dudx[k] = 0.0;
      dudy[k] = 0.0;
      dvdx[k] = 0.0;
      dvdy[k] = 0.0;
    }
    for (int iEdge = 1; iEdge <= m.nEdgesOnCell(iCell); ++iEdge) {  // DM:167-179 (iEdge: the local edge number)
      const int e = m.edgesOnCell(iEdge, iCell);
      const real c2 = m.deformation_coef_c2(iEdge, iCell), s2 = m.deformation_coef_s2(iEdge, iCell);
      const real cs = m.deformation_coef_cs(iEdge, iCell);
      for (int k = 1; k <= nVertLevels; ++k) {
        dudx[k] = dudx[k] + c2 * u(k, e) - cs * v(k, e);
        dudy[k] = dudy[k] + cs * u(k, e) - s2 * v(k, e);
        dvdx[k] = dvdx[k] + cs * u(k, e) + c2 * v(k, e);
        dvdy[k] = dvdy[k] + s2 * u(k, e) + cs * v(k, e);
      }
    }
    for (int k = 1; k <= nVertLevels; ++k) {  // DM:183-191
      const real d_11 = 2 * dudx[k];
      const real d_22 = 2 * dvdy[k];
      const real d_12 = dudy[k] + dvdx[k];
      const real kd = kfac * std::sqrt(0.25 * ((d_11 - d_22) * (d_11 - d_22)) + d_12 * d_12);
      kdiff(k, iCell) = kd < kmax ? kd : kmax;  // min(kdiff, bound): minsd / minpd keep the bound for a NaN
    }
  }
  // DM:199-200: config_len_disp**3 as (len*len)*len (0x264575-0x26458c)
  h_mom_eddy_visc4 = config_visc4_2dsmag * (config_len_disp * config_len_disp * config_len_disp);
  h_theta_eddy_visc4 = h_mom_eddy_visc4;
}

void u_dissipation_3d(const Mesh& m, const DissipationArgs& a, const FArray<real>& eddy_visc_horz,
                      FArray<real>& delsq_u, FArray<real>& delsq_vorticity, FArray<real>& delsq_divergence,
                      const FArray<real>& u, const FArray<real>& divergence,
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
  if (a.h_mom_eddy_visc4 > 0.0) {  // DM:728-811
    for (int iVertex = 1; iVertex <= m.nVertices; ++iVertex) {  // DM:735-751
      for (int k = 1; k <= nVertLevels; ++k) delsq_vorticity(k, iVertex) = 0.0;
      for (int i = 1; i <= m.vertexDegree; ++i) {
        const int iEdge = m.edgesOnVertex(i, iVertex);
        const real edge_sign = m.invAreaTriangle(iVertex) * m.dcEdge(iEdge) * m.edgesOnVertex_sign(i, iVertex);
        for (int k = 1; k <= nVertLevels; ++k)
          delsq_vorticity(k, iVertex) = delsq_vorticity(k, iVertex) + edge_sign * delsq_u(k, iEdge);
      }
    }
    for (int iCell = 1; iCell <= m.nCells; ++iCell) {  // DM:754-772
      for (int k = 1; k <= nVertLevels; ++k) delsq_divergence(k, iCell) = 0.0;
      const real r = m.invAreaCell(iCell);
      for (int i = 1; i <= m.nEdgesOnCell(iCell); ++i) {
        const int iEdge = m.edgesOnCell(i, iCell);
        const real edge_sign = r * m.dvEdge(iEdge) * m.edgesOnCell_sign(i, iCell);
        for (int k = 1; k <= nVertLevels; ++k)
          delsq_divergence(k, iCell) = delsq_divergence(k, iCell) + edge_sign * delsq_u(k, iEdge);
      }
    }
    for (int iEdge = 1; iEdge <= m.nEdgesSolve(); ++iEdge) {  // DM:781-807
      const int cell1 = m.cellsOnEdge(1, iEdge);
      const int cell2 = m.cellsOnEdge(2, iEdge);
      const int vertex1 = m.verticesOnEdge(1, iEdge);
      const int vertex2 = m.verticesOnEdge(2, iEdge);
      const real u_mix_scale = m.meshScalingDel4(iEdge) * a.h_mom_eddy_visc4;
      const real r_dc = u_mix_scale * a.config_del4u_div_factor * m.invDcEdge(iEdge);
      const real r_dv = u_mix_scale * std::min(m.invDvEdge(iEdge), 4 * m.invDcEdge(iEdge));
      for (int k = 1; k <= nVertLevels; ++k) {
        const real u_diffusion = rho_edge(k, iEdge) * ((delsq_divergence(k, cell2) - delsq_divergence(k, cell1)) * r_dc -
                                                       (delsq_vorticity(k, vertex2) - delsq_vorticity(k, vertex1)) * r_dv);
        tend_u_euler(k, iEdge) = tend_u_euler(k, iEdge) - u_diffusion;  // DM:804
      }
    }
  }
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
  }  // v_mom_eddy_visc2 <= 0: no vertical mixing (JW)
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
  if (a.h_mom_eddy_visc4 > 0.0) {  // DM:1073-1100 (cellSolveStart..cellSolveEnd)
    for (int iCell = 1; iCell <= m.nCellsSolve(); ++iCell) {
      const real r_areaCell = a.h_mom_eddy_visc4 * m.invAreaCell(iCell);
      for (int i = 1; i <= m.nEdgesOnCell(iCell); ++i) {
        const int iEdge = m.edgesOnCell(i, iCell);
        const int cell1 = m.cellsOnEdge(1, iEdge);
        const int cell2 = m.cellsOnEdge(2, iEdge);
        const real edge_sign =
            m.meshScalingDel4(iEdge) * r_areaCell * m.dvEdge(iEdge) * m.edgesOnCell_sign(i, iCell) * m.invDcEdge(iEdge);
        for (int k = 2; k <= nVertLevels; ++k)
          tend_w_euler(k, iCell) = tend_w_euler(k, iCell) - edge_sign * (delsq_w(k, cell2) - delsq_w(k, cell1));
      }
    }
  }
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
  }  // v_mom_eddy_visc2 <= 0: no vertical mixing (JW)
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
  if (a.h_theta_eddy_visc4 > 0.0) {  // DM:1306-1332 (cellSolveStart..cellSolveEnd)
    for (int iCell = 1; iCell <= m.nCellsSolve(); ++iCell) {
      const real r_areaCell = a.h_theta_eddy_visc4 * a.prandtl_inv * m.invAreaCell(iCell);
      for (int i = 1; i <= m.nEdgesOnCell(iCell); ++i) {
        const int iEdge = m.edgesOnCell(i, iCell);
        const real edge_sign =
            m.meshScalingDel4(iEdge) * r_areaCell * m.dvEdge(iEdge) * m.edgesOnCell_sign(i, iCell) * m.invDcEdge(iEdge);
        const int cell1 = m.cellsOnEdge(1, iEdge);
        const int cell2 = m.cellsOnEdge(2, iEdge);
        for (int k = 1; k <= nVertLevels; ++k)
          tend_theta_euler(k, iCell) =
              tend_theta_euler(k, iCell) - edge_sign * (delsq_theta(k, cell2) - delsq_theta(k, cell1));
      }
    }
  }
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
  }  // v_theta_eddy_visc2 <= 0: no vertical mixing (JW)
}

}  // namespace cppmpas
