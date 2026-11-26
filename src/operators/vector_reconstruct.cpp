// SPDX-License-Identifier: GPL-3.0-or-later
// src/operators/mpas_vector_reconstruction.F:205-339 (VR below). Each loop over cells is a kernel, one thread per cell
// (core/backend.hpp), the k loop and the edgesOnCell loop inside the thread in MPAS's order.
#include "operators/vector_reconstruct.hpp"

#include <string>
#include <utility>

#include "core/error.hpp"
#include "core/libm.hpp"
#include "core/views.hpp"

namespace cppmpas {

namespace {
// u is (nVertLevels, nEdges+1) and every output (nVertLevels, nCells+1), as MPAS allocates them
void need_exact(const FArray<real>& a, int nL, int n1, const char* name) {
  if (a.rank() != 2 || a.extent(0) != nL || a.extent(1) != n1)
    throw Error(std::string("mpas_reconstruct_2d: ") + name + " is not allocated (nVertLevels, " + std::to_string(n1) +
                ")");
}
}  // namespace

void mpas_reconstruct_2d(const Mesh& m, const FArray<real>& u, FArray<real>& uReconstructX,
                         FArray<real>& uReconstructY, FArray<real>& uReconstructZ, FArray<real>& uReconstructZonal,
                         FArray<real>& uReconstructMeridional) {
  const int nCells = m.nCellsSolve();  // VR:236-240, 255-259: includeHalos absent (every MPAS-Atmosphere call)
  const int nVertLevels = m.nVertLevels;
  const auto& coeffs_reconstruct = m.coeffs_reconstruct;
  if (nVertLevels < 1) throw Error("mpas_reconstruct_2d: nVertLevels < 1 (no vertical mesh)");
  need_exact(u, nVertLevels, m.nEdges + 1, "u");
  for (auto [a, n] : {std::pair{&uReconstructX, "uReconstructX"}, {&uReconstructY, "uReconstructY"},
                      {&uReconstructZ, "uReconstructZ"}, {&uReconstructZonal, "uReconstructZonal"},
                      {&uReconstructMeridional, "uReconstructMeridional"}})
    need_exact(*a, nVertLevels, nCells + 1, n);
  if (coeffs_reconstruct.rank() != 3 || coeffs_reconstruct.extent(0) != 3 ||
      coeffs_reconstruct.extent(1) < m.maxEdges || coeffs_reconstruct.extent(2) < nCells)
    throw Error("mpas_reconstruct_2d: coeffs_reconstruct is not (3, maxEdges, nCells) (mpas_init_reconstruct first)");

  parallel_for(1, nCells, [=, m = view(m), coeffs_reconstruct = view(coeffs_reconstruct), u = view(u),
                           uReconstructX = view(uReconstructX), uReconstructY = view(uReconstructY),
                           uReconstructZ = view(uReconstructZ)] CPPMPAS_LAMBDA(int iCell) {  // VR:274-299
    for (int k = 1; k <= nVertLevels; ++k) {        // VR:277-281 (0.0: default real, exact)
      uReconstructX(k, iCell) = 0.0;
      uReconstructY(k, iCell) = 0.0;
      uReconstructZ(k, iCell) = 0.0;
    }
    for (int i = 1; i <= m.nEdgesOnCell(iCell); ++i) {  // VR:286-298
      const int iEdge = m.edgesOnCell(i, iCell);
      for (int k = 1; k <= nVertLevels; ++k) {
        uReconstructX(k, iCell) = uReconstructX(k, iCell) + coeffs_reconstruct(1, i, iCell) * u(k, iEdge);
        uReconstructY(k, iCell) = uReconstructY(k, iCell) + coeffs_reconstruct(2, i, iCell) * u(k, iEdge);
        uReconstructZ(k, iCell) = uReconstructZ(k, iCell) + coeffs_reconstruct(3, i, iCell) * u(k, iEdge);
      }
    }
  });

  if (m.on_a_sphere) {
    parallel_for(1, nCells, [=, m = view(m), uReconstructX = view(uReconstructX), uReconstructY = view(uReconstructY),
                             uReconstructZ = view(uReconstructZ), uReconstructZonal = view(uReconstructZonal),
                             uReconstructMeridional = view(uReconstructMeridional)] CPPMPAS_LAMBDA(int iCell) {  // VR:309-322
      // VR:310-313: gfortran fuses cos / sin of the same argument into one sincos@plt each, latCell then lonCell
      // (reference executable, __mpas_vector_reconstruction_MOD_mpas_reconstruct_2d at 0x588510: sincos@plt at
      // 0x589a57 and 0x589a8b, the only libm calls; no sin / cos / libmvec call, no FMA; the k loop is vectorised two
      // lanes wide with mulpd / addpd / subpd)
      real slat, clat, slon, clon;
      libm::sincos(m.latCell(iCell), &slat, &clat);
      libm::sincos(m.lonCell(iCell), &slon, &clon);
      for (int k = 1; k <= nVertLevels; ++k) {
        // VR:316-320. Fortran's unary minus applies to the product: -a*b + c is -(a*b) + c. The reference evaluates
        // these as c - a*b (subpd / subsd at 0x589b67, 0x589b82, 0x589bf0, 0x589c18, and 0x589f01, 0x589f2d,
        // 0x58a001, 0x58a042 in the second copy of the loop): no negation, so a NaN keeps its sign. Equal to
        // (-a)*b + c for every non-NaN value, signed zeros included; g++ makes the same c - a*b of this form.
        uReconstructZonal(k, iCell) = -(uReconstructX(k, iCell) * slon) + uReconstructY(k, iCell) * clon;
        uReconstructMeridional(k, iCell) =
            -((uReconstructX(k, iCell) * clon + uReconstructY(k, iCell) * slon) * slat) + uReconstructZ(k, iCell) * clat;
      }
    });
  } else {
    parallel_for(1, nCells, [=, uReconstructX = view(uReconstructX), uReconstructY = view(uReconstructY),
                             uReconstructZonal = view(uReconstructZonal),
                             uReconstructMeridional = view(uReconstructMeridional)] CPPMPAS_LAMBDA(int iCell) {  // VR:329-334
      for (int k = 1; k <= nVertLevels; ++k) {
        uReconstructZonal(k, iCell) = uReconstructX(k, iCell);
        uReconstructMeridional(k, iCell) = uReconstructY(k, iCell);
      }
    });
  }
}

}  // namespace cppmpas
