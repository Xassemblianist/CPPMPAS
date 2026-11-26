// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// atm_compute_solve_diagnostics (src/core_atmosphere/dynamics/mpas_atm_time_integration.F:6967-7437): the diagnostics
// of the horizontal flow that the dycore computes once before the first step (mpas_atm_core.F:532) and after every
// acoustic / RK stage (mpas_atm_time_integration.F:2556). One rank, one thread.
#include "core/mesh.hpp"

namespace cppmpas {

// Module work arrays of atm_time_integration (mpas_atm_time_integration.F:98-99), allocated (nVertLevels, n+1) by
// the caller (mpas_atm_core.F:492-499, mpas_atm_time_integration.F:2543-2550) with the garbage column zeroed.
struct SolveDiagWork {
  FArray<real> ke_vertex;  // (nVertLevels, nVertices+1)
  FArray<real> ke_edge;    // (nVertLevels, nEdges+1)
};
// mpas_atm_core.F:492-499: allocate(ke_vertex(nVertLevels,nVertices+1)), allocate(ke_edge(nVertLevels,nEdges+1)),
// ke_vertex(k,nVertices+1) = 0, ke_edge(k,nEdges+1) = 0. The other elements are uninitialised in MPAS and are all
// written by atm_compute_solve_diagnostics before they are read (edge range 1..nEdges, vertex range 1..nVertices).
void allocate_solve_diag_work(const Mesh& mesh, SolveDiagWork& w);

// The diag-pool fields atm_compute_solve_diagnostics and mpas_reconstruct write, with MPAS's extents and the
// Registry default (0: Registry.xml:1846-1908 give no default_value).
void allocate_solve_diag_fields(const Mesh& mesh, Diag& diag);

// mpas_atm_time_integration.F:6967-7437 (atm_compute_solve_diagnostics + _work). h = rho_zz and u of the time level
// in use. rk_step = 0 means the optional argument is absent (v is then computed, :7307-7310).
// config_apvm_upwinding <= 0 (no APVM, :7372 not taken) is an error: no configuration tested takes it.
void atm_compute_solve_diagnostics(real dt, real config_apvm_upwinding, const Mesh& mesh, const FArray<real>& h,
                                   const FArray<real>& u, Diag& diag, SolveDiagWork& w, int cellStart, int cellEnd,
                                   int vertexStart, int vertexEnd, int edgeStart, int edgeEnd, int rk_step = 0);

}  // namespace cppmpas
