// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// src/operators/mpas_vector_operations.F and mpas_rbf_interpolation.F / mpas_vector_reconstruction.F: the init-time
// set-up of the cell tangent planes and of the reconstruction coefficients.
#include "core/mesh.hpp"

namespace cppmpas {

// mpas_vector_operations.F:95-101
CPPMPAS_HD void mpas_unit_vec_in_r3(real xin[3]);
// mpas_vector_operations.F:113-121
CPPMPAS_HD void mpas_cross_product_in_r3(const real p_1[3], const real p_2[3], real p_out[3]);
// mpas_vector_operations.F:652-771 (via mpas_rbf_interp_initialize, mpas_rbf_interpolation.F:131-139): recomputes
// edgeNormalVectors, localVerticalUnitVectors and cellTangentPlane
void mpas_initialize_vectors(Mesh& mesh);

// mpas_rbf_interpolation.F:1079-1145. Arrays as Fortran (pointCount, 3) column-major, 0-based storage here:
// sourcePoints[(i-1) + pointCount*(d-1)]; planeBasisVectors[(r-1) + 2*(d-1)]. MPAS's automatic arrays are the
// caller's thread's own here (kernel work arrays): in the CUDA build pointCount <= kMaxRbfPoints (the CPU build has
// no limit).
inline constexpr int kMaxRbfPoints = 16;
CPPMPAS_HD void mpas_rbf_interp_func_3D_plane_vec_const_dir_comp_coeffs(int pointCount, const real* sourcePoints,
                                                                       const real* unitVectors,
                                                                       const real destinationPoint[3], real alpha,
                                                                       const real* planeBasisVectors,
                                                                       real* coefficients);
// mpas_rbf_interpolation.F:1670-1698 with elgs 1782-1846; A (N,N) column-major; CUDA build: N <= kMaxRbfPoints + 2
CPPMPAS_HD void mpas_legs(real* A, int N, real* B, real* X, int* INDX);

// mpas_vector_reconstruction.F:61-191 (includeHalos absent: nCellsSolve)
void mpas_init_reconstruct(Mesh& mesh);

}  // namespace cppmpas
