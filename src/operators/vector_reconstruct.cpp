// SPDX-License-Identifier: GPL-3.0-or-later
// src/operators/mpas_vector_reconstruction.F:205-339 (VR below).
#include "operators/vector_reconstruct.hpp"

#include "core/libm.hpp"

namespace cppmpas {

void mpas_reconstruct_2d(const Mesh& m, const FArray<real>& u, FArray<real>& uReconstructX,
                         FArray<real>& uReconstructY, FArray<real>& uReconstructZ, FArray<real>& uReconstructZonal,
                         FArray<real>& uReconstructMeridional) {
  const int nCells = m.nCellsSolve();  // VR:259-263, includeHalos absent
  const int nVertLevels = m.nVertLevels;
  const auto& coeffs_reconstruct = m.coeffs_reconstruct;

  for (int iCell = 1; iCell <= nCells; ++iCell) {  // VR:274-299
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
  }

  if (m.on_a_sphere) {
    for (int iCell = 1; iCell <= nCells; ++iCell) {  // VR:309-322
      // VR:310-313: gfortran fuses cos / sin of the same argument into one sincos@plt each (latCell, then lonCell);
      // the iCell loop is not vectorised (gfortran 14.2 -O3 object code of this routine: two sincos calls, no libmvec)
      real slat, clat, slon, clon;
      libm::sincos(m.latCell(iCell), &slat, &clat);
      libm::sincos(m.lonCell(iCell), &slon, &clon);
      for (int k = 1; k <= nVertLevels; ++k) {
        uReconstructZonal(k, iCell) = -uReconstructX(k, iCell) * slon + uReconstructY(k, iCell) * clon;
        uReconstructMeridional(k, iCell) =
            -(uReconstructX(k, iCell) * clon + uReconstructY(k, iCell) * slon) * slat + uReconstructZ(k, iCell) * clat;
      }
    }
  } else {
    for (int iCell = 1; iCell <= nCells; ++iCell) {  // VR:329-334
      for (int k = 1; k <= nVertLevels; ++k) {
        uReconstructZonal(k, iCell) = uReconstructX(k, iCell);
        uReconstructMeridional(k, iCell) = uReconstructY(k, iCell);
      }
    }
  }
}

}  // namespace cppmpas
