// SPDX-License-Identifier: GPL-3.0-or-later
// atm_compute_solve_diagnostics, src/core_atmosphere/dynamics/mpas_atm_time_integration.F:6967-7437 (TI below).
// Floating point as the reference executable (bin/atmosphere_model of the mpas-build asset 153c7206..., objdump -d:
// __atm_time_integration_MOD_atm_compute_solve_diagnostics_work at 0x104ed0, called from the wrapper at 0x108ecd):
// no FMA, no libm or libmvec call (its only calls are memset@plt at 0x1050dc, 0x1057d4, 0x105e09, 0x106d0d,
// 0x1073bf), no negation, SSE2 only; the k loops are vectorised two lanes wide with mulpd / addpd / subpd, which round
// each lane exactly as the scalar statement does, so every statement below is the Fortran statement evaluated left
// to right. Real literals are default real promoted to double (-fdefault-real-8): 0.5, 0.25, 1.0 - .375 = 0.625
// and 1. - ke_fact = 0.375 are exact.
#include "core_atmosphere/atm_solve_diagnostics.hpp"

#include <string>
#include <utility>

#include "core/error.hpp"
#include "core_atmosphere/atm_views.hpp"

namespace cppmpas {

void allocate_solve_diag_work(const Mesh& m, SolveDiagWork& w) {
  const int nL = m.nVertLevels;
  if (nL < 1) throw Error("allocate_solve_diag_work: nVertLevels < 1 (no vertical mesh)");
  w.ke_vertex.allocate({nL, m.nVertices + 1}, 0.0);  // mpas_atm_core.F:492, :497
  w.ke_edge.allocate({nL, m.nEdges + 1}, 0.0);       // mpas_atm_core.F:493, :498
}

void allocate_solve_diag_fields(const Mesh& m, Diag& d) {
  const int nL = m.nVertLevels, nC = m.nCells + 1, nE = m.nEdges + 1, nV = m.nVertices + 1;
  if (nL < 1) throw Error("allocate_solve_diag_fields: nVertLevels < 1 (no vertical mesh)");
  for (auto* a : {&d.v, &d.rho_edge, &d.pv_edge, &d.gradPVt, &d.gradPVn}) a->allocate({nL, nE}, 0.0);
  for (auto* a : {&d.vorticity, &d.pv_vertex}) a->allocate({nL, nV}, 0.0);
  for (auto* a : {&d.divergence, &d.ke, &d.pv_cell, &d.uReconstructX, &d.uReconstructY, &d.uReconstructZ,
                  &d.uReconstructZonal, &d.uReconstructMeridional})
    a->allocate({nL, nC}, 0.0);
}

namespace {
[[noreturn]] void fail(const std::string& what) { throw Error("atm_compute_solve_diagnostics: " + what); }
// an array the routine reads or writes: exactly MPAS's (nVertLevels, n+1)
void need_exact(const FArray<real>& a, int nL, int n1, const char* name) {
  if (a.rank() != 2 || a.extent(0) != nL || a.extent(1) != n1)
    fail(std::string(name) + " is not allocated (nVertLevels, " + std::to_string(n1) + ")");
}
// a mesh array: rank r and at least the extents the routine indexes (d0, then 1..n)
template <class T> void need_mesh(const FArray<T>& a, int d0, int n, const char* name) {
  const bool ok = d0 == 0 ? (a.rank() == 1 && a.extent(0) >= n) : (a.rank() == 2 && a.extent(0) >= d0 && a.extent(1) >= n);
  if (!ok) fail(std::string("mesh array ") + name + " is not allocated with the extents it reads");
}
}  // namespace

// TI:6967-7058 gathers the arrays and calls the work routine TI:7061-7437; the work routine is inlined here.
void atm_compute_solve_diagnostics(real dt, real config_apvm_upwinding, const Mesh& m, const FArray<real>& h,
                                   const FArray<real>& u, Diag& d, SolveDiagWork& w, int cellStart, int cellEnd,
                                   int vertexStart, int vertexEnd, int edgeStart, int edgeEnd,
                                   std::optional<int> rk_step) {
  const int nVertLevels = m.nVertLevels;
  const int vertexDegree = m.vertexDegree;
  if (vertexDegree != 3) fail("vertexDegree /= 3 (TI:7263 sums EdgesOnVertex(1:3) only)");
  if (!(config_apvm_upwinding > 0.0)) fail("config_apvm_upwinding <= 0 (no APVM upwinding, TI:7372) is not ported");
  if (nVertLevels < 1) fail("nVertLevels < 1 (no vertical mesh)");
  // one thread (mpas_atm_core.F:502-517, mpas_atm_threading.F:114-125): the whole ranges
  if (cellStart != 1 || cellEnd != m.nCells || vertexStart != 1 || vertexEnd != m.nVertices || edgeStart != 1 ||
      edgeEnd != m.nEdges)
    fail("the cell / vertex / edge ranges must be 1..nCells, 1..nVertices, 1..nEdges (one thread; ranges " +
         std::to_string(cellStart) + ".." + std::to_string(cellEnd) + ", " + std::to_string(vertexStart) + ".." +
         std::to_string(vertexEnd) + ", " + std::to_string(edgeStart) + ".." + std::to_string(edgeEnd) + ")");
  if (rk_step && (*rk_step < 1 || *rk_step > 3))
    fail("rk_step = " + std::to_string(*rk_step) + " (atm_srk3 passes 1, 2, 3: TI:2556)");
  const int nC1 = m.nCells + 1, nE1 = m.nEdges + 1, nV1 = m.nVertices + 1;
  need_exact(h, nVertLevels, nC1, "rho_zz");
  need_exact(u, nVertLevels, nE1, "u");
  for (auto [a, n] : {std::pair{&d.v, "v"}, {&d.rho_edge, "rho_edge"}, {&d.pv_edge, "pv_edge"},
                      {&d.gradPVt, "gradPVt"}, {&d.gradPVn, "gradPVn"}, {&w.ke_edge, "ke_edge"}})
    need_exact(*a, nVertLevels, nE1, n);
  for (auto [a, n] : {std::pair{&d.vorticity, "vorticity"}, {&d.pv_vertex, "pv_vertex"}, {&w.ke_vertex, "ke_vertex"}})
    need_exact(*a, nVertLevels, nV1, n);
  for (auto [a, n] : {std::pair{&d.divergence, "divergence"}, {&d.ke, "ke"}, {&d.pv_cell, "pv_cell"}})
    need_exact(*a, nVertLevels, nC1, n);
  need_mesh(m.nEdgesOnCell, 0, m.nCells, "nEdgesOnCell");
  need_mesh(m.nEdgesOnEdge, 0, m.nEdges, "nEdgesOnEdge");
  for (auto [a, n] : {std::pair{&m.invAreaCell, "invAreaCell"}}) need_mesh(*a, 0, m.nCells, n);
  for (auto [a, n] : {std::pair{&m.dcEdge, "dcEdge"}, {&m.dvEdge, "dvEdge"}, {&m.invDvEdge, "invDvEdge"},
                      {&m.invDcEdge, "invDcEdge"}})
    need_mesh(*a, 0, m.nEdges, n);
  for (auto [a, n] : {std::pair{&m.invAreaTriangle, "invAreaTriangle"}, {&m.fVertex, "fVertex"}})
    need_mesh(*a, 0, m.nVertices, n);
  need_mesh(m.cellsOnEdge, 2, m.nEdges, "cellsOnEdge");
  need_mesh(m.verticesOnEdge, 2, m.nEdges, "verticesOnEdge");
  need_mesh(m.edgesOnEdge, m.maxEdges2, m.nEdges, "edgesOnEdge");
  need_mesh(m.weightsOnEdge, m.maxEdges2, m.nEdges, "weightsOnEdge");
  need_mesh(m.edgesOnVertex, 3, m.nVertices, "edgesOnVertex");
  need_mesh(m.edgesOnVertex_sign, 3, m.nVertices, "edgesOnVertex_sign");
  need_mesh(m.kiteAreasOnVertex, 3, m.nVertices, "kiteAreasOnVertex");
  need_mesh(m.edgesOnCell, m.maxEdges, m.nCells, "edgesOnCell");
  need_mesh(m.edgesOnCell_sign, m.maxEdges, m.nCells, "edgesOnCell_sign");
  need_mesh(m.verticesOnCell, m.maxEdges, m.nCells, "verticesOnCell");
  need_mesh(m.kiteForCell, m.maxEdges, m.nCells, "kiteForCell");
  // kiteForCell(i, iCell) indexes kiteAreasOnVertex(j, iVertex) (TI:7291, 7390); atm_compute_signs leaves it 0 when
  // it logs "Unexpected error while identifying kiteForCell" (mpas_atm_core.F:1262-1264), and MPAS would then read
  // kiteAreasOnVertex(0, iVertex), outside the array
  for (int iCell = 1; iCell <= m.nCells; ++iCell)
    for (int i = 1; i <= m.nEdgesOnCell(iCell); ++i)
      if (m.kiteForCell(i, iCell) < 1 || m.kiteForCell(i, iCell) > vertexDegree)
        fail("kiteForCell(" + std::to_string(i) + ", " + std::to_string(iCell) + ") = " +
             std::to_string(m.kiteForCell(i, iCell)) + " outside 1..vertexDegree (mpas_atm_core.F:1262-1264)");
  // kernel views of the arrays (core/views.hpp); each loop over elements is one kernel, one thread per element, with
  // the k loop and the neighbour loops in the thread in MPAS's order (core/backend.hpp)
  auto h_edge = d.rho_edge.view();
  auto vorticity = d.vorticity.view();
  auto divergence = d.divergence.view();
  auto ke = d.ke.view();
  auto v = d.v.view();
  auto pv_vertex = d.pv_vertex.view();
  auto pv_edge = d.pv_edge.view();
  auto pv_cell = d.pv_cell.view();
  auto gradPVt = d.gradPVt.view();
  auto gradPVn = d.gradPVn.view();
  auto ke_edge = w.ke_edge.view();
  auto ke_vertex = w.ke_vertex.view();

  // TI:7136-7154 height (rho_zz) and kinetic energy times area at edges
  parallel_for(edgeStart, edgeEnd, [=, m = view(m), h = view(h), u = view(u)] CPPMPAS_LAMBDA(int iEdge) {
    const int cell1 = m.cellsOnEdge(1, iEdge);
    const int cell2 = m.cellsOnEdge(2, iEdge);
    for (int k = 1; k <= nVertLevels; ++k) h_edge(k, iEdge) = 0.5 * (h(k, cell1) + h(k, cell2));
    const real efac = m.dcEdge(iEdge) * m.dvEdge(iEdge);
    for (int k = 1; k <= nVertLevels; ++k) ke_edge(k, iEdge) = efac * (u(k, iEdge) * u(k, iEdge));  // u**2 = u*u
  });

  // TI:7160-7180 circulation and relative vorticity at vertices
  parallel_for(vertexStart, vertexEnd, [=, m = view(m), u = view(u)] CPPMPAS_LAMBDA(int iVertex) {
    for (int k = 1; k <= nVertLevels; ++k) vorticity(k, iVertex) = 0.0;
    for (int i = 1; i <= vertexDegree; ++i) {
      const int iEdge = m.edgesOnVertex(i, iVertex);
      const real s = m.edgesOnVertex_sign(i, iVertex) * m.dcEdge(iEdge);
      for (int k = 1; k <= nVertLevels; ++k) vorticity(k, iVertex) = vorticity(k, iVertex) + s * u(k, iEdge);
    }
    for (int k = 1; k <= nVertLevels; ++k) vorticity(k, iVertex) = vorticity(k, iVertex) * m.invAreaTriangle(iVertex);
  });

  // TI:7187-7207 divergence at cell centres
  parallel_for(cellStart, cellEnd, [=, m = view(m), u = view(u)] CPPMPAS_LAMBDA(int iCell) {
    for (int k = 1; k <= nVertLevels; ++k) divergence(k, iCell) = 0.0;
    for (int i = 1; i <= m.nEdgesOnCell(iCell); ++i) {
      const int iEdge = m.edgesOnCell(i, iCell);
      const real s = m.edgesOnCell_sign(i, iCell) * m.dvEdge(iEdge);
      for (int k = 1; k <= nVertLevels; ++k) divergence(k, iCell) = divergence(k, iCell) + s * u(k, iEdge);
    }
    const real r = m.invAreaCell(iCell);
    for (int k = 1; k <= nVertLevels; ++k) divergence(k, iCell) = divergence(k, iCell) * r;
  });

  // TI:7220-7239 kinetic energy at cell centres
  parallel_for(cellStart, cellEnd, [=, m = view(m)] CPPMPAS_LAMBDA(int iCell) {
    for (int k = 1; k <= nVertLevels; ++k) ke(k, iCell) = 0.0;
    for (int i = 1; i <= m.nEdgesOnCell(iCell); ++i) {
      const int iEdge = m.edgesOnCell(i, iCell);
      for (int k = 1; k <= nVertLevels; ++k) ke(k, iCell) = ke(k, iCell) + 0.25 * ke_edge(k, iEdge);
    }
    for (int k = 1; k <= nVertLevels; ++k) ke(k, iCell) = ke(k, iCell) * m.invAreaCell(iCell);
  });

  // TI:7243-7301 hollingsworth = .true. (a parameter, TI:7125): AG's KE construction
  parallel_for(vertexStart, vertexEnd, [=, m = view(m)] CPPMPAS_LAMBDA(int iVertex) {  // TI:7253-7266
    const real r = 0.25 * m.invAreaTriangle(iVertex);
    const int e1 = m.edgesOnVertex(1, iVertex), e2 = m.edgesOnVertex(2, iVertex), e3 = m.edgesOnVertex(3, iVertex);
    for (int k = 1; k <= nVertLevels; ++k)
      ke_vertex(k, iVertex) = (ke_edge(k, e1) + ke_edge(k, e2) + ke_edge(k, e3)) * r;
  });
  const real ke_fact = 1.0 - .375;  // TI:7274
  parallel_for(cellStart, cellEnd, [=] CPPMPAS_LAMBDA(int iCell) {  // TI:7278-7282
    for (int k = 1; k <= nVertLevels; ++k) ke(k, iCell) = ke_fact * ke(k, iCell);
  });
  parallel_for(cellStart, cellEnd, [=, m = view(m)] CPPMPAS_LAMBDA(int iCell) {  // TI:7286-7298
    const real r = m.invAreaCell(iCell);
    for (int i = 1; i <= m.nEdgesOnCell(iCell); ++i) {
      const int iVertex = m.verticesOnCell(i, iCell);
      const int j = m.kiteForCell(i, iCell);
      for (int k = 1; k <= nVertLevels; ++k)
        ke(k, iCell) = ke(k, iCell) + (1. - ke_fact) * m.kiteAreasOnVertex(j, iVertex) * ke_vertex(k, iVertex) * r;
    }
  });

  // TI:7307-7331 tangential velocity (Thuburn et al. 2009): rk_step absent (init) or 3; at rk_steps 1 and 2 v keeps
  // the values it has (those of the last rk_step 3, or of the init), and TI:7428 below uses them
  bool reconstruct_v = true;
  if (rk_step.has_value()) {
    if (*rk_step != 3) reconstruct_v = false;
  }
  if (reconstruct_v) {
    parallel_for(edgeStart, edgeEnd, [=, m = view(m), u = view(u)] CPPMPAS_LAMBDA(int iEdge) {
      for (int k = 1; k <= nVertLevels; ++k) v(k, iEdge) = 0.0;
      for (int i = 1; i <= m.nEdgesOnEdge(iEdge); ++i) {
        const int eoe = m.edgesOnEdge(i, iEdge);
        for (int k = 1; k <= nVertLevels; ++k) v(k, iEdge) = v(k, iEdge) + m.weightsOnEdge(i, iEdge) * u(k, eoe);
      }
    });
  }

  // TI:7341-7353 absolute vorticity at vertices
  parallel_for(vertexStart, vertexEnd, [=, m = view(m)] CPPMPAS_LAMBDA(int iVertex) {
    for (int k = 1; k <= nVertLevels; ++k) pv_vertex(k, iVertex) = (m.fVertex(iVertex) + vorticity(k, iVertex));
  });

  // TI:7364-7369 pv at edges
  parallel_for(edgeStart, edgeEnd, [=, m = view(m)] CPPMPAS_LAMBDA(int iEdge) {
    for (int k = 1; k <= nVertLevels; ++k)
      pv_edge(k, iEdge) = 0.5 * (pv_vertex(k, m.verticesOnEdge(1, iEdge)) + pv_vertex(k, m.verticesOnEdge(2, iEdge)));
  });

  // TI:7372-7434 APVM upwinding (config_apvm_upwinding > 0, checked above)
  parallel_for(cellStart, cellEnd, [=, m = view(m)] CPPMPAS_LAMBDA(int iCell) {  // TI:7381-7397
    for (int k = 1; k <= nVertLevels; ++k) pv_cell(k, iCell) = 0.0;
    const real r = m.invAreaCell(iCell);
    for (int i = 1; i <= m.nEdgesOnCell(iCell); ++i) {
      const int iVertex = m.verticesOnCell(i, iCell);
      const int j = m.kiteForCell(i, iCell);
      for (int k = 1; k <= nVertLevels; ++k)
        pv_cell(k, iCell) = pv_cell(k, iCell) + m.kiteAreasOnVertex(j, iVertex) * pv_vertex(k, iVertex) * r;
    }
  });
  const real r = config_apvm_upwinding * dt;  // TI:7417
  parallel_for(edgeStart, edgeEnd, [=, m = view(m), u = view(u)] CPPMPAS_LAMBDA(int iEdge) {  // TI:7420-7430
    const real r1 = 1.0 * m.invDvEdge(iEdge);
    const real r2 = 1.0 * m.invDcEdge(iEdge);
    const int v1 = m.verticesOnEdge(1, iEdge), v2 = m.verticesOnEdge(2, iEdge);
    const int c1 = m.cellsOnEdge(1, iEdge), c2 = m.cellsOnEdge(2, iEdge);
    for (int k = 1; k <= nVertLevels; ++k) {
      gradPVt(k, iEdge) = (pv_vertex(k, v2) - pv_vertex(k, v1)) * r1;
      gradPVn(k, iEdge) = (pv_cell(k, c2) - pv_cell(k, c1)) * r2;
      pv_edge(k, iEdge) = pv_edge(k, iEdge) - r * (v(k, iEdge) * gradPVt(k, iEdge) + u(k, iEdge) * gradPVn(k, iEdge));
    }
  });
}

}  // namespace cppmpas
