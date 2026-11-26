// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// Kernel views of the mesh pool, the state and the diag fields (core/mesh.hpp): trivially copyable structs with the
// members of Mesh / State / Diag under the same names, each array an FView of the same storage. A kernel captures
// view(m) under the name m (core/backend.hpp), so MPAS's statements (m.cellsOnEdge(1, iEdge), s.u(k, iEdge), ...)
// read the same in the kernel as on the host. view of a const object gives read-only arrays.
// The lists name every array of the structs; a new member is added to its list here.
#include <type_traits>

#include "core/mesh.hpp"

namespace cppmpas {

template <class T> FView<T> view(FArray<T>& a) { return a.view(); }
template <class T> FView<const T> view(const FArray<T>& a) { return a.view(); }
template <class T> FView<T> view(FView<T> v) { return v; }

// a = v over the whole array, garbage element included (one thread per element)
template <class T> void fill_all(FArray<T>& a, T v) {
  const FView<T> p = a.view();
  parallel_for_n(a.size(), [=] CPPMPAS_LAMBDA(std::size_t i) { p[i] = v; });
}

#define CPPMPAS_MESH_INT_ARRAYS(X)                                                                                    \
  X(indexToCellID) X(indexToEdgeID) X(indexToVertexID) X(nEdgesOnCell) X(nEdgesOnEdge) X(cellsOnCell) X(edgesOnCell) \
  X(verticesOnCell) X(cellsOnEdge) X(verticesOnEdge) X(edgesOnEdge) X(cellsOnVertex) X(edgesOnVertex) X(bdyMaskCell)  \
  X(bdyMaskEdge) X(bdyMaskVertex) X(kiteForCell) X(advCellsForEdge) X(nAdvCellsForEdge) X(nearestRelaxationCell)
#define CPPMPAS_MESH_REAL_ARRAYS(X)                                                                                   \
  X(latCell) X(lonCell) X(xCell) X(yCell) X(zCell) X(latEdge) X(lonEdge) X(xEdge) X(yEdge) X(zEdge) X(latVertex)      \
  X(lonVertex) X(xVertex) X(yVertex) X(zVertex) X(weightsOnEdge) X(dvEdge) X(dcEdge) X(angleEdge) X(areaCell)        \
  X(areaTriangle) X(kiteAreasOnVertex) X(fEdge) X(fVertex) X(meshDensity) X(nominalMinDc) X(edgeNormalVectors)       \
  X(localVerticalUnitVectors) X(cellTangentPlane) X(cf1) X(cf2) X(cf3) X(zgrid) X(rdzw) X(dzu) X(rdzu) X(fzm) X(fzp) \
  X(zxu) X(zz) X(zb) X(zb3) X(dss) X(deriv_two) X(coeffs_reconstruct) X(invAreaCell) X(invDvEdge) X(invDcEdge)       \
  X(invAreaTriangle) X(edgesOnCell_sign) X(edgesOnVertex_sign) X(zb_cell) X(zb3_cell) X(adv_coefs) X(adv_coefs_3rd)  \
  X(deformation_coef_c2) X(deformation_coef_s2) X(deformation_coef_cs) X(deformation_coef_c) X(deformation_coef_s)   \
  X(meshScalingDel2) X(meshScalingDel4) X(meshScalingRegionalCell) X(meshScalingRegionalEdge) X(etp) X(etm) X(ewp)   \
  X(ewm) X(specZoneMaskCell) X(specZoneMaskEdge) X(specZoneMaskVertex)
#define CPPMPAS_STATE_REAL_ARRAYS(X)                                                                                  \
  X(u_init) X(v_init) X(t_init) X(u) X(w) X(rho) X(theta) X(rho_base) X(theta_base) X(surface_pressure) X(scalars)    \
  X(theta_m) X(rho_zz) X(Time)
#define CPPMPAS_DIAG_REAL_ARRAYS(X)                                                                                   \
  X(ru) X(rw) X(rho_p) X(rtheta_base) X(rtheta_p) X(exner) X(exner_base) X(pressure_p) X(pressure_base) X(v)         \
  X(rho_edge) X(pv_edge) X(gradPVt) X(gradPVn) X(vorticity) X(pv_vertex) X(divergence) X(ke) X(pv_cell)              \
  X(uReconstructX) X(uReconstructY) X(uReconstructZ) X(uReconstructZonal) X(uReconstructMeridional)

// view element type: const T for a read-only view
template <bool C, class T> using VElem = std::conditional_t<C, const T, T>;

template <bool C> struct MeshViewT {
  int nCells, nEdges, nVertices, maxEdges, maxEdges2, vertexDegree, nVertLevels;
  bool on_a_sphere, is_periodic;
  double sphere_radius, x_period, y_period;
  CPPMPAS_HD int nCellsSolve() const { return nCells; }
  CPPMPAS_HD int nEdgesSolve() const { return nEdges; }
  CPPMPAS_HD int nVerticesSolve() const { return nVertices; }
#define CPPMPAS_X(n) FView<VElem<C, int>> n;
  CPPMPAS_MESH_INT_ARRAYS(CPPMPAS_X)
#undef CPPMPAS_X
#define CPPMPAS_X(n) FView<VElem<C, real>> n;
  CPPMPAS_MESH_REAL_ARRAYS(CPPMPAS_X)
#undef CPPMPAS_X
};
template <bool C> struct StateViewT {
#define CPPMPAS_X(n) FView<VElem<C, real>> n;
  CPPMPAS_STATE_REAL_ARRAYS(CPPMPAS_X)
#undef CPPMPAS_X
};
template <bool C> struct DiagViewT {
#define CPPMPAS_X(n) FView<VElem<C, real>> n;
  CPPMPAS_DIAG_REAL_ARRAYS(CPPMPAS_X)
#undef CPPMPAS_X
};

template <class M> auto view_mesh(M& m) {
  MeshViewT<std::is_const_v<M>> v;
  v.nCells = m.nCells;
  v.nEdges = m.nEdges;
  v.nVertices = m.nVertices;
  v.maxEdges = m.maxEdges;
  v.maxEdges2 = m.maxEdges2;
  v.vertexDegree = m.vertexDegree;
  v.nVertLevels = m.nVertLevels;
  v.on_a_sphere = m.on_a_sphere;
  v.is_periodic = m.is_periodic;
  v.sphere_radius = m.sphere_radius;
  v.x_period = m.x_period;
  v.y_period = m.y_period;
#define CPPMPAS_X(n) v.n = m.n.view();
  CPPMPAS_MESH_INT_ARRAYS(CPPMPAS_X)
  CPPMPAS_MESH_REAL_ARRAYS(CPPMPAS_X)
#undef CPPMPAS_X
  return v;
}
inline MeshViewT<true> view(const Mesh& m) { return view_mesh(m); }
inline MeshViewT<false> view(Mesh& m) { return view_mesh(m); }

template <class S> auto view_state(S& s) {
  StateViewT<std::is_const_v<S>> v;
#define CPPMPAS_X(n) v.n = s.n.view();
  CPPMPAS_STATE_REAL_ARRAYS(CPPMPAS_X)
#undef CPPMPAS_X
  return v;
}
inline StateViewT<true> view(const State& s) { return view_state(s); }
inline StateViewT<false> view(State& s) { return view_state(s); }

template <class D> auto view_diag(D& d) {
  DiagViewT<std::is_const_v<D>> v;
#define CPPMPAS_X(n) v.n = d.n.view();
  CPPMPAS_DIAG_REAL_ARRAYS(CPPMPAS_X)
#undef CPPMPAS_X
  return v;
}
inline DiagViewT<true> view(const Diag& d) { return view_diag(d); }
inline DiagViewT<false> view(Diag& d) { return view_diag(d); }

}  // namespace cppmpas
