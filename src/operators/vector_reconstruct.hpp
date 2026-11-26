// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// mpas_reconstruct (src/operators/mpas_vector_reconstruction.F:41-44 -> mpas_reconstruct_2d, :205-339): the cell-centre
// reconstruction of an edge-normal vector field from coeffs_reconstruct.
#include "core/mesh.hpp"

namespace cppmpas {

// mpas_vector_reconstruction.F:205-339, includeHalos absent (cells 1..nCellsSolve). u (nVertLevels, nEdges+1);
// outputs (nVertLevels, nCells+1).
void mpas_reconstruct_2d(const Mesh& mesh, const FArray<real>& u, FArray<real>& uReconstructX,
                         FArray<real>& uReconstructY, FArray<real>& uReconstructZ, FArray<real>& uReconstructZonal,
                         FArray<real>& uReconstructMeridional);

}  // namespace cppmpas
