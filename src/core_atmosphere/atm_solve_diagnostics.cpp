// SPDX-License-Identifier: GPL-3.0-or-later
// atm_compute_solve_diagnostics, src/core_atmosphere/dynamics/mpas_atm_time_integration.F:6967-7437 (TI below).
// Floating point as the reference object code (gfortran 14.2 -O3, x86-64 without -march): no FMA, no libm call;
// the k loops are vectorised two lanes wide with mulpd / addpd / subpd, which round each lane exactly as the scalar
// statement does, so every statement below is the Fortran statement evaluated left to right. Real literals are
// default real promoted to double (-fdefault-real-8): 0.5, 0.25, 1.0 - .375 = 0.625 and 1. - ke_fact = 0.375 are
// exact.
#include "core_atmosphere/atm_solve_diagnostics.hpp"

#include <algorithm>

#include "core/error.hpp"

namespace cppmpas {

void allocate_solve_diag_work(const Mesh& m, SolveDiagWork& w) {
  const int nL = m.nVertLevels;
  w.ke_vertex.allocate({nL, m.nVertices + 1}, 0.0);  // mpas_atm_core.F:492, :497
  w.ke_edge.allocate({nL, m.nEdges + 1}, 0.0);       // mpas_atm_core.F:493, :498
}

void allocate_solve_diag_fields(const Mesh& m, Diag& d) {
  const int nL = m.nVertLevels, nC = m.nCells + 1, nE = m.nEdges + 1, nV = m.nVertices + 1;
  for (auto* a : {&d.v, &d.rho_edge, &d.pv_edge, &d.gradPVt, &d.gradPVn}) a->allocate({nL, nE}, 0.0);
  for (auto* a : {&d.vorticity, &d.pv_vertex}) a->allocate({nL, nV}, 0.0);
  for (auto* a : {&d.divergence, &d.ke, &d.pv_cell, &d.uReconstructX, &d.uReconstructY, &d.uReconstructZ,
                  &d.uReconstructZonal, &d.uReconstructMeridional})
    a->allocate({nL, nC}, 0.0);
}

// TI:6967-7058 gathers the arrays and calls the work routine TI:7061-7437; the work routine is inlined here.
void atm_compute_solve_diagnostics(real dt, real config_apvm_upwinding, const Mesh& m, const FArray<real>& h,
                                   const FArray<real>& u, Diag& d, SolveDiagWork& w, int cellStart, int cellEnd,
                                   int vertexStart, int vertexEnd, int edgeStart, int edgeEnd, int rk_step) {
  const int nVertLevels = m.nVertLevels;
  const int vertexDegree = m.vertexDegree;
  if (vertexDegree != 3)
    throw Error("atm_compute_solve_diagnostics: vertexDegree /= 3 (TI:7263 sums EdgesOnVertex(1:3) only)");
  if (!(config_apvm_upwinding > 0.0))
    throw Error("atm_compute_solve_diagnostics: config_apvm_upwinding <= 0 (no APVM upwinding, TI:7372) is not ported");
  auto& h_edge = d.rho_edge;
  auto& vorticity = d.vorticity;
  auto& divergence = d.divergence;
  auto& ke = d.ke;
  auto& v = d.v;
  auto& pv_vertex = d.pv_vertex;
  auto& pv_edge = d.pv_edge;
  auto& pv_cell = d.pv_cell;
  auto& gradPVt = d.gradPVt;
  auto& gradPVn = d.gradPVn;
  auto& ke_edge = w.ke_edge;
  auto& ke_vertex = w.ke_vertex;

  // TI:7136-7154 height (rho_zz) and kinetic energy times area at edges
  for (int iEdge = edgeStart; iEdge <= edgeEnd; ++iEdge) {
    const int cell1 = m.cellsOnEdge(1, iEdge);
    const int cell2 = m.cellsOnEdge(2, iEdge);
    for (int k = 1; k <= nVertLevels; ++k) h_edge(k, iEdge) = 0.5 * (h(k, cell1) + h(k, cell2));
    const real efac = m.dcEdge(iEdge) * m.dvEdge(iEdge);
    for (int k = 1; k <= nVertLevels; ++k) ke_edge(k, iEdge) = efac * (u(k, iEdge) * u(k, iEdge));  // u**2 = u*u
  }

  // TI:7160-7180 circulation and relative vorticity at vertices
  for (int iVertex = vertexStart; iVertex <= vertexEnd; ++iVertex) {
    for (int k = 1; k <= nVertLevels; ++k) vorticity(k, iVertex) = 0.0;
    for (int i = 1; i <= vertexDegree; ++i) {
      const int iEdge = m.edgesOnVertex(i, iVertex);
      const real s = m.edgesOnVertex_sign(i, iVertex) * m.dcEdge(iEdge);
      for (int k = 1; k <= nVertLevels; ++k) vorticity(k, iVertex) = vorticity(k, iVertex) + s * u(k, iEdge);
    }
    for (int k = 1; k <= nVertLevels; ++k) vorticity(k, iVertex) = vorticity(k, iVertex) * m.invAreaTriangle(iVertex);
  }

  // TI:7187-7207 divergence at cell centres
  for (int iCell = cellStart; iCell <= cellEnd; ++iCell) {
    for (int k = 1; k <= nVertLevels; ++k) divergence(k, iCell) = 0.0;
    for (int i = 1; i <= m.nEdgesOnCell(iCell); ++i) {
      const int iEdge = m.edgesOnCell(i, iCell);
      const real s = m.edgesOnCell_sign(i, iCell) * m.dvEdge(iEdge);
      for (int k = 1; k <= nVertLevels; ++k) divergence(k, iCell) = divergence(k, iCell) + s * u(k, iEdge);
    }
    const real r = m.invAreaCell(iCell);
    for (int k = 1; k <= nVertLevels; ++k) divergence(k, iCell) = divergence(k, iCell) * r;
  }

  // TI:7220-7239 kinetic energy at cell centres
  for (int iCell = cellStart; iCell <= cellEnd; ++iCell) {
    for (int k = 1; k <= nVertLevels; ++k) ke(k, iCell) = 0.0;
    for (int i = 1; i <= m.nEdgesOnCell(iCell); ++i) {
      const int iEdge = m.edgesOnCell(i, iCell);
      for (int k = 1; k <= nVertLevels; ++k) ke(k, iCell) = ke(k, iCell) + 0.25 * ke_edge(k, iEdge);
    }
    for (int k = 1; k <= nVertLevels; ++k) ke(k, iCell) = ke(k, iCell) * m.invAreaCell(iCell);
  }

  // TI:7243-7301 hollingsworth = .true. (a parameter, TI:7125): AG's KE construction
  for (int iVertex = vertexStart; iVertex <= vertexEnd; ++iVertex) {  // TI:7253-7266
    const real r = 0.25 * m.invAreaTriangle(iVertex);
    const int e1 = m.edgesOnVertex(1, iVertex), e2 = m.edgesOnVertex(2, iVertex), e3 = m.edgesOnVertex(3, iVertex);
    for (int k = 1; k <= nVertLevels; ++k)
      ke_vertex(k, iVertex) = (ke_edge(k, e1) + ke_edge(k, e2) + ke_edge(k, e3)) * r;
  }
  const real ke_fact = 1.0 - .375;  // TI:7274
  for (int iCell = cellStart; iCell <= cellEnd; ++iCell)  // TI:7278-7282
    for (int k = 1; k <= nVertLevels; ++k) ke(k, iCell) = ke_fact * ke(k, iCell);
  for (int iCell = cellStart; iCell <= cellEnd; ++iCell) {  // TI:7286-7298
    const real r = m.invAreaCell(iCell);
    for (int i = 1; i <= m.nEdgesOnCell(iCell); ++i) {
      const int iVertex = m.verticesOnCell(i, iCell);
      const int j = m.kiteForCell(i, iCell);
      for (int k = 1; k <= nVertLevels; ++k)
        ke(k, iCell) = ke(k, iCell) + (1. - ke_fact) * m.kiteAreasOnVertex(j, iVertex) * ke_vertex(k, iVertex) * r;
    }
  }

  // TI:7307-7331 tangential velocity (Thuburn et al. 2009), at init and at rk_step 3 only
  bool reconstruct_v = true;
  if (rk_step != 0 && rk_step != 3) reconstruct_v = false;
  if (reconstruct_v) {
    for (int iEdge = edgeStart; iEdge <= edgeEnd; ++iEdge) {
      for (int k = 1; k <= nVertLevels; ++k) v(k, iEdge) = 0.0;
      for (int i = 1; i <= m.nEdgesOnEdge(iEdge); ++i) {
        const int eoe = m.edgesOnEdge(i, iEdge);
        for (int k = 1; k <= nVertLevels; ++k) v(k, iEdge) = v(k, iEdge) + m.weightsOnEdge(i, iEdge) * u(k, eoe);
      }
    }
  }

  // TI:7341-7353 absolute vorticity at vertices
  for (int iVertex = vertexStart; iVertex <= vertexEnd; ++iVertex)
    for (int k = 1; k <= nVertLevels; ++k) pv_vertex(k, iVertex) = (m.fVertex(iVertex) + vorticity(k, iVertex));

  // TI:7364-7369 pv at edges
  for (int iEdge = edgeStart; iEdge <= edgeEnd; ++iEdge)
    for (int k = 1; k <= nVertLevels; ++k)
      pv_edge(k, iEdge) = 0.5 * (pv_vertex(k, m.verticesOnEdge(1, iEdge)) + pv_vertex(k, m.verticesOnEdge(2, iEdge)));

  // TI:7373-7432 APVM upwinding (config_apvm_upwinding > 0, checked above)
  for (int iCell = cellStart; iCell <= cellEnd; ++iCell) {  // TI:7381-7397
    for (int k = 1; k <= nVertLevels; ++k) pv_cell(k, iCell) = 0.0;
    const real r = m.invAreaCell(iCell);
    for (int i = 1; i <= m.nEdgesOnCell(iCell); ++i) {
      const int iVertex = m.verticesOnCell(i, iCell);
      const int j = m.kiteForCell(i, iCell);
      for (int k = 1; k <= nVertLevels; ++k)
        pv_cell(k, iCell) = pv_cell(k, iCell) + m.kiteAreasOnVertex(j, iVertex) * pv_vertex(k, iVertex) * r;
    }
  }
  const real r = config_apvm_upwinding * dt;  // TI:7417
  for (int iEdge = edgeStart; iEdge <= edgeEnd; ++iEdge) {  // TI:7420-7430
    const real r1 = 1.0 * m.invDvEdge(iEdge);
    const real r2 = 1.0 * m.invDcEdge(iEdge);
    const int v1 = m.verticesOnEdge(1, iEdge), v2 = m.verticesOnEdge(2, iEdge);
    const int c1 = m.cellsOnEdge(1, iEdge), c2 = m.cellsOnEdge(2, iEdge);
    for (int k = 1; k <= nVertLevels; ++k) {
      gradPVt(k, iEdge) = (pv_vertex(k, v2) - pv_vertex(k, v1)) * r1;
      gradPVn(k, iEdge) = (pv_cell(k, c2) - pv_cell(k, c1)) * r2;
      pv_edge(k, iEdge) = pv_edge(k, iEdge) - r * (v(k, iEdge) * gradPVt(k, iEdge) + u(k, iEdge) * gradPVn(k, iEdge));
    }
  }
}

}  // namespace cppmpas
