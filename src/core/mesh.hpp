// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// MPAS-Atmosphere's mesh pool (var_struct "mesh", src/core_atmosphere/Registry.xml:1356-1670) and the initial state
// of the input stream (Registry.xml:602-623), held as MPAS holds them on one MPI rank with one block:
//  * cells in file order (src/framework/mpas_block_decomp.F:303-315);
//  * edges and vertices in MPAS's local order: first seen while walking edgesOnCell / verticesOnCell cell by cell
//    (mpas_bootstrapping.F:261-266 -> mpas_block_creator.F:318 -> mpas_block_decomp.F:633-685, then
//    mpas_block_decomp.F:565-620 keeps the owned order);
//  * connectivity as local 1-based indices, unmatched IDs -> n+1 (mpas_stream_manager.F:5619-5636);
//  * one garbage element on every nCells / nEdges / nVertices dimension, holding the Registry default value
//    (0 for indexTo*ID, mpas_block_creator.F:1069, 1095, 1122).
#include <string>
#include <vector>

#include "core/field.hpp"
#include "core/precision.hpp"

namespace cppmpas {

enum class Element { Cell, Edge, Vertex };
const char* element_dim_name(Element e);  // "nCells", "nEdges", "nVertices"

// Local <-> file (global ID) numbering of one element kind. Index 0 is unused (1-based on both sides).
struct ElementOrder {
  std::vector<int> local_to_file;  // local index i (1..n) -> global ID = 1-based record in the file
  std::vector<int> file_to_local;  // global ID g (1..n) -> local index
  int n() const { return local_to_file.empty() ? 0 : static_cast<int>(local_to_file.size()) - 1; }
  bool identity() const;
};

struct Mesh {
  // dimensions (nCellsSolve == nCells etc.: one rank, one block, no halo elements)
  int nCells = 0, nEdges = 0, nVertices = 0, maxEdges = 0, maxEdges2 = 0, vertexDegree = 0;
  int nVertLevels = 0;  // 0 when only a grid (no vertical mesh) was read
  int nCellsSolve() const { return nCells; }
  int nEdgesSolve() const { return nEdges; }
  int nVerticesSolve() const { return nVertices; }

  // global attributes (mpas_bootstrapping.F:287-379)
  bool on_a_sphere = true;
  double sphere_radius = 1.0;
  bool is_periodic = false;
  double x_period = 0.0, y_period = 0.0;
  std::string file_id, parent_id, mesh_spec;

  // index maps, (n+1) with garbage element 0
  FArray<int> indexToCellID, indexToEdgeID, indexToVertexID;
  // connectivity, local indices
  FArray<int> nEdgesOnCell, nEdgesOnEdge;
  FArray<int> cellsOnCell, edgesOnCell, verticesOnCell;   // (maxEdges, nCells+1)
  FArray<int> cellsOnEdge, verticesOnEdge;                 // (TWO, nEdges+1)
  FArray<int> edgesOnEdge;                                 // (maxEdges2, nEdges+1)
  FArray<int> cellsOnVertex, edgesOnVertex;                // (vertexDegree, nVertices+1)
  FArray<int> bdyMaskCell, bdyMaskEdge, bdyMaskVertex;
  // horizontal geometry
  FArray<real> latCell, lonCell, xCell, yCell, zCell, latEdge, lonEdge, xEdge, yEdge, zEdge;
  FArray<real> latVertex, lonVertex, xVertex, yVertex, zVertex;
  FArray<real> weightsOnEdge, dvEdge, dcEdge, angleEdge, areaCell, areaTriangle, kiteAreasOnVertex;
  FArray<real> fEdge, fVertex, meshDensity, nominalMinDc;
  FArray<real> edgeNormalVectors, localVerticalUnitVectors, cellTangentPlane;
  // vertical mesh (init file)
  FArray<real> cf1, cf2, cf3, zgrid, rdzw, dzu, rdzu, fzm, fzp, zxu, zz, zb, zb3, dss, deriv_two, coeffs_reconstruct;
  // computed by atm_mpas_init_block (src/core_atmosphere/mpas_atm_core.F:366-617), Registry.xml:1434-1669;
  // allocated with the Registry default (0) by allocate_init_fields()
  FArray<real> invAreaCell, invDvEdge, invDcEdge, invAreaTriangle;
  FArray<real> edgesOnCell_sign, edgesOnVertex_sign, zb_cell, zb3_cell;
  FArray<int> kiteForCell;
  FArray<int> advCellsForEdge, nAdvCellsForEdge;
  FArray<real> adv_coefs, adv_coefs_3rd;
  FArray<real> deformation_coef_c2, deformation_coef_s2, deformation_coef_cs, deformation_coef_c, deformation_coef_s;
  FArray<real> meshScalingDel2, meshScalingDel4, meshScalingRegionalCell, meshScalingRegionalEdge;
  FArray<real> etp, etm, ewp, ewm;
  FArray<real> specZoneMaskCell, specZoneMaskEdge, specZoneMaskVertex;
  FArray<int> nearestRelaxationCell;

  ElementOrder cell_order, edge_order, vertex_order;
  const ElementOrder& order(Element e) const;
};

// The initial-state arrays of the atmosphere input stream the dycore reads (Registry.xml:602-623).
struct State {
  std::string xtime, initial_time;          // raw StrLen (64) characters as stored in the file
  FArray<real> u_init, v_init, t_init;
  FArray<real> u, w, rho, theta, rho_base, theta_base, surface_pressure;
  std::vector<std::string> scalar_names;    // constituents of var_array scalars, in Registry order
  FArray<real> scalars;                     // (num_scalars, nVertLevels, nCells+1)
  // time level 1 of the prognostic state that atm_init_coupled_diagnostics sets (mpas_atm_time_integration.F:7498,
  // 7500); u, w and scalars above are time level 1 as read
  FArray<real> theta_m, rho_zz;             // (nVertLevels, nCells+1)
};

// The diag pool fields atm_init_coupled_diagnostics sets (mpas_atm_time_integration.F:7501-7523); rho and theta of
// the diag pool are State::rho / State::theta (read from the input stream).
struct Diag {
  FArray<real> ru;                                       // (nVertLevels, nEdges+1)
  FArray<real> rw;                                       // (nVertLevelsP1, nCells+1)
  FArray<real> rho_p, rtheta_base, rtheta_p, exner, exner_base, pressure_p, pressure_base;  // (nVertLevels, nCells+1)
  // set by atm_compute_solve_diagnostics (mpas_atm_time_integration.F:6967-7437, Registry.xml:1846-1908) and
  // mpas_reconstruct (mpas_vector_reconstruction.F:205-339, Registry.xml:1882-1894); allocated by
  // allocate_solve_diag_fields (core_atmosphere/atm_solve_diagnostics.hpp)
  FArray<real> v, rho_edge, pv_edge, gradPVt, gradPVn;                // (nVertLevels, nEdges+1)
  FArray<real> vorticity, pv_vertex;                                  // (nVertLevels, nVertices+1)
  FArray<real> divergence, ke, pv_cell;                               // (nVertLevels, nCells+1)
  FArray<real> uReconstructX, uReconstructY, uReconstructZ, uReconstructZonal, uReconstructMeridional;
};

// MPAS's local numbering on one rank from the file's cell connectivity (global IDs, Fortran order):
// nEdgesOnCell(nCells), edgesOnCell / verticesOnCell (maxEdges, nCells), cellsOnEdge (2, nEdges),
// cellsOnVertex (vertexDegree, nVertices). Fills mesh.cell_order / edge_order / vertex_order and indexTo*ID.
void build_one_rank_order(Mesh& mesh, const std::vector<int>& nEdgesOnCell, const std::vector<int>& edgesOnCell,
                          const std::vector<int>& verticesOnCell, const std::vector<int>& cellsOnEdge,
                          const std::vector<int>& cellsOnVertex);

// Copy an array whose slowest (last Fortran) dimension is the element dimension between MPAS local order and file
// order. inner = number of values per element. The garbage element is not part of the file.
template <class T>
void local_to_file(const ElementOrder& o, std::size_t inner, const T* local, T* file) {
  for (int g = 1; g <= o.n(); ++g) {
    const T* src = local + static_cast<std::size_t>(o.file_to_local[g] - 1) * inner;
    T* dst = file + static_cast<std::size_t>(g - 1) * inner;
    for (std::size_t r = 0; r < inner; ++r) dst[r] = src[r];
  }
}
template <class T>
void file_to_local(const ElementOrder& o, std::size_t inner, const T* file, T* local) {
  for (int i = 1; i <= o.n(); ++i) {
    const T* src = file + static_cast<std::size_t>(o.local_to_file[i] - 1) * inner;
    T* dst = local + static_cast<std::size_t>(i - 1) * inner;
    for (std::size_t r = 0; r < inner; ++r) dst[r] = src[r];
  }
}

}  // namespace cppmpas
