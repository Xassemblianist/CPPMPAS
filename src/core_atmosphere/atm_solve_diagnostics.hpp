// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// atm_compute_solve_diagnostics (src/core_atmosphere/dynamics/mpas_atm_time_integration.F:6967-7437): the diagnostics
// of the horizontal flow that the dycore computes once before the first step (mpas_atm_core.F:532) and after every
// acoustic / RK stage (mpas_atm_time_integration.F:2556). One rank, one thread.
#include <optional>

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
// Registry default (0: Registry.xml:1846-1909 give no default_value).
void allocate_solve_diag_fields(const Mesh& mesh, Diag& diag);

// mpas_atm_time_integration.F:6967-7437 (atm_compute_solve_diagnostics + _work). h = rho_zz and u of the time level
// in use. rk_step is Fortran's optional argument (TI:6983, 7116): std::nullopt is "not present" (the call at
// mpas_atm_core.F:532-535); atm_srk3 passes 1, 2, 3 (TI:2556). v is recomputed when rk_step is absent or 3 and left
// as it is for 1 and 2 (TI:7307-7310); any other present value is an error (no MPAS caller passes one).
// One rank, one thread: the ranges must be 1..nCells, 1..nVertices, 1..nEdges (mpas_atm_threading.F:114-125 with
// nThreads = 1; the reference is built without OpenMP, so the !$OMP BARRIERs at TI:7211, 7269, 7356, 7401 are
// comments there); anything else is an error, as are arrays of other extents and kiteForCell outside 1..3.
// config_apvm_upwinding <= 0 or NaN (no APVM, :7372 not taken) is an error: no configuration tested takes it.
void atm_compute_solve_diagnostics(real dt, real config_apvm_upwinding, const Mesh& mesh, const FArray<real>& h,
                                   const FArray<real>& u, Diag& diag, SolveDiagWork& w, int cellStart, int cellEnd,
                                   int vertexStart, int vertexEnd, int edgeStart, int edgeEnd,
                                   std::optional<int> rk_step = std::nullopt);

}  // namespace cppmpas
