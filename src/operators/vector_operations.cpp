// SPDX-License-Identifier: GPL-3.0-or-later
#include "operators/vector_operations.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

#include "core/error.hpp"
#include "operators/geometry_utils.hpp"

namespace cppmpas {

// gfortran expands SUM(a(1:n)*b(1:n)) with an accumulator that starts at +0.0 and adds the products in index order
// (reference executable: pxor then addsd per term), so a sum of -0.0 terms is +0.0.

// mpas_vector_operations.F:95-101
void mpas_unit_vec_in_r3(real xin[3]) {
  const real mag = std::sqrt(xin[0] * xin[0] + xin[1] * xin[1] + xin[2] * xin[2]);  // :99
  xin[0] = xin[0] / mag;                                                              // :100
  xin[1] = xin[1] / mag;
  xin[2] = xin[2] / mag;
}

// mpas_vector_operations.F:113-121
void mpas_cross_product_in_r3(const real p_1[3], const real p_2[3], real p_out[3]) {
  p_out[0] = p_1[1] * p_2[2] - p_1[2] * p_2[1];
  p_out[1] = p_1[2] * p_2[0] - p_1[0] * p_2[2];
  p_out[2] = p_1[0] * p_2[1] - p_1[1] * p_2[0];
}

// mpas_vector_operations.F:652-771
void mpas_initialize_vectors(Mesh& m) {
  if (m.is_periodic)
    throw Error("mpas_initialize_vectors: the is_periodic branches (mpas_vector_operations.F:717, 729, 742) are not ported");
  const int nCells = m.nCells, nEdges = m.nEdges;
  // the boundary-edge branches (cell1 or cell2 == nCells+1, :722-724, :734-736) are not verified
  for (int iEdge = 1; iEdge <= nEdges; ++iEdge)
    if (m.cellsOnEdge(1, iEdge) == nCells + 1 || m.cellsOnEdge(2, iEdge) == nCells + 1)
      throw Error("mpas_initialize_vectors: edge " + std::to_string(iEdge) +
                  " is a boundary edge (cellsOnEdge = nCells+1): that branch is not verified");
  auto& lvu = m.localVerticalUnitVectors;
  auto& env = m.edgeNormalVectors;
  auto& ctp = m.cellTangentPlane;
  std::fill(env.data(), env.data() + env.size(), 0.0);  // :693-694, whole arrays
  std::fill(lvu.data(), lvu.data() + lvu.size(), 0.0);

  for (int iCell = 1; iCell <= nCells; ++iCell) {  // :697-706
    if (m.on_a_sphere) {
      real v[3] = {m.xCell(iCell), m.yCell(iCell), m.zCell(iCell)};
      mpas_unit_vec_in_r3(v);
      for (int d = 1; d <= 3; ++d) lvu(d, iCell) = v[d - 1];
    } else {
      lvu(1, iCell) = 0.;
      lvu(2, iCell) = 0.;
      lvu(3, iCell) = 1.;
    }
  }

  for (int iEdge = 1; iEdge <= nEdges; ++iEdge) {  // :711-754, is_periodic = .false. branches
    const int cell1 = m.cellsOnEdge(1, iEdge);
    const int cell2 = m.cellsOnEdge(2, iEdge);
    real v[3];
    if (cell1 == nCells + 1) {  // :722-724
      v[0] = m.xCell(cell2) - m.xEdge(iEdge);
      v[1] = m.yCell(cell2) - m.yEdge(iEdge);
      v[2] = m.zCell(cell2) - m.zEdge(iEdge);
    } else if (cell2 == nCells + 1) {  // :734-736
      v[0] = m.xEdge(iEdge) - m.xCell(cell1);
      v[1] = m.yEdge(iEdge) - m.yCell(cell1);
      v[2] = m.zEdge(iEdge) - m.zCell(cell1);
    } else {  // :747-749
      v[0] = m.xCell(cell2) - m.xCell(cell1);
      v[1] = m.yCell(cell2) - m.yCell(cell1);
      v[2] = m.zCell(cell2) - m.zCell(cell1);
    }
    mpas_unit_vec_in_r3(v);  // :753
    for (int d = 1; d <= 3; ++d) env(d, iEdge) = v[d - 1];
  }

  for (int iCell = 1; iCell <= nCells; ++iCell) {  // :756-769
    const int iEdge = m.edgesOnCell(1, iCell);
    const real rHat[3] = {lvu(1, iCell), lvu(2, iCell), lvu(3, iCell)};
    real normalDotRHat = 0.0;  // :761 SUM
    for (int d = 1; d <= 3; ++d) normalDotRHat = normalDotRHat + env(d, iEdge) * rHat[d - 1];
    real xHatPlane[3];
    for (int d = 1; d <= 3; ++d) xHatPlane[d - 1] = env(d, iEdge) - normalDotRHat * rHat[d - 1];  // :762
    mpas_unit_vec_in_r3(xHatPlane);
    real yHatPlane[3];
    mpas_cross_product_in_r3(rHat, xHatPlane, yHatPlane);  // :765
    mpas_unit_vec_in_r3(yHatPlane);
    for (int d = 1; d <= 3; ++d) {  // :767-768
      ctp(d, 1, iCell) = xHatPlane[d - 1];
      ctp(d, 2, iCell) = yHatPlane[d - 1];
    }
  }
}

namespace {

// evaluate_rbf, mpas_rbf_interpolation.F:1369-1376 (integer 1 -> 1.0)
real evaluate_rbf(real rSquared) { return 1.0 / std::sqrt(1.0 + rSquared); }

// mpas_set_up_vector_dirichlet_rbf_matrix_and_rhs, mpas_rbf_interpolation.F:1527-1559; arrays (pointCount, dims)
void set_up_vector_dirichlet_rbf_matrix_and_rhs(int pointCount, int dimensions, const real* sourcePoints,
                                                const real* unitVectors, const real* destinationPoint, real alpha,
                                                real* matrix, real* rhs) {
  auto S = [&](int i, int d) { return sourcePoints[(i - 1) + pointCount * (d - 1)]; };
  auto U = [&](int i, int d) { return unitVectors[(i - 1) + pointCount * (d - 1)]; };
  auto M = [&](int i, int j) -> real& { return matrix[(i - 1) + pointCount * (j - 1)]; };
  for (int j = 1; j <= pointCount; ++j) {
    for (int i = j; i <= pointCount; ++i) {
      real s = 0.0;  // :1546 SUM((p_i - p_j)**2) / alpha**2
      for (int d = 1; d <= dimensions; ++d) {
        const real t = S(i, d) - S(j, d);
        s = s + t * t;
      }
      const real rSquared = s / (alpha * alpha);
      const real rbfValue = evaluate_rbf(rSquared);
      real dot = 0.0;  // :1548
      for (int d = 1; d <= dimensions; ++d) dot = dot + U(i, d) * U(j, d);
      M(i, j) = rbfValue * dot;  // :1549-1550
      M(j, i) = M(i, j);
    }
  }
  for (int j = 1; j <= pointCount; ++j) {  // :1554-1557
    real s = 0.0;
    for (int d = 1; d <= dimensions; ++d) {
      const real t = destinationPoint[d - 1] - S(j, d);
      s = s + t * t;
    }
    const real rSquared = s / (alpha * alpha);
    for (int d = 1; d <= dimensions; ++d) rhs[(j - 1) + pointCount * (d - 1)] = evaluate_rbf(rSquared) * U(j, d);
  }
}

// elgs, mpas_rbf_interpolation.F:1782-1846 (inlined into mpas_legs in the reference)
void elgs(real* A, int N, int* INDX) {
  auto a = [&](int i, int j) -> real& { return A[(i - 1) + N * (j - 1)]; };
  std::vector<real> C(static_cast<std::size_t>(N) + 1);
  for (int I = 1; I <= N; ++I) INDX[I - 1] = I;
  for (int I = 1; I <= N; ++I) {  // :1805-1812
    real C1 = 0.0;
    for (int J = 1; J <= N; ++J) C1 = maxsd(C1, std::fabs(a(I, J)));
    C[I] = C1;
  }
  for (int J = 1; J <= N - 1; ++J) {  // :1816-1844
    real PI1 = 0.0;
    int K = 0;
    for (int I = J; I <= N; ++I) {
      const real PI = std::fabs(a(INDX[I - 1], J)) / C[INDX[I - 1]];
      if (PI > PI1) {
        PI1 = PI;
        K = I;
      }
    }
    // Without a pivot, MPAS reuses K from the previous column's DO K loop (N+1: out of bounds); not reproduced.
    if (K == 0) throw Error("mpas_legs/elgs: no pivot in column " + std::to_string(J) + " (singular RBF matrix)");
    const int ITMP = INDX[J - 1];
    INDX[J - 1] = INDX[K - 1];
    INDX[K - 1] = ITMP;
    for (int I = J + 1; I <= N; ++I) {
      const real PJ = a(INDX[I - 1], J) / a(INDX[J - 1], J);
      a(INDX[I - 1], J) = PJ;
      for (int KK = J + 1; KK <= N; ++KK) a(INDX[I - 1], KK) = a(INDX[I - 1], KK) - PJ * a(INDX[J - 1], KK);
    }
  }
}

}  // namespace

// mpas_rbf_interpolation.F:1670-1698
void mpas_legs(real* A, int N, real* B, real* X, int* INDX) {
  // MPAS's only caller passes pointCount + 2 >= 3; N < 1 would index X(0) (public routine: checked)
  if (N < 1) throw Error("mpas_legs: N < 1");
  auto a = [&](int i, int j) { return A[(i - 1) + N * (j - 1)]; };
  elgs(A, N, INDX);
  for (int I = 1; I <= N - 1; ++I)
    for (int J = I + 1; J <= N; ++J) B[INDX[J - 1] - 1] = B[INDX[J - 1] - 1] - a(INDX[J - 1], I) * B[INDX[I - 1] - 1];
  X[N - 1] = B[INDX[N - 1] - 1] / a(INDX[N - 1], N);
  for (int I = N - 1; I >= 1; --I) {
    X[I - 1] = B[INDX[I - 1] - 1];
    for (int J = I + 1; J <= N; ++J) X[I - 1] = X[I - 1] - a(INDX[I - 1], J) * X[J - 1];
    X[I - 1] = X[I - 1] / a(INDX[I - 1], I);
  }
}

// mpas_rbf_interpolation.F:1079-1145
void mpas_rbf_interp_func_3D_plane_vec_const_dir_comp_coeffs(int pointCount, const real* sourcePoints,
                                                            const real* unitVectors, const real destinationPoint[3],
                                                            real alpha, const real* planeBasisVectors,
                                                            real* coefficients) {
  const int pc = pointCount;
  const int matrixSize = pc + 2;  // :1104
  auto SP = [&](int i, int d) { return sourcePoints[(i - 1) + pc * (d - 1)]; };
  auto UV = [&](int i, int d) { return unitVectors[(i - 1) + pc * (d - 1)]; };
  auto PB = [&](int r, int d) { return planeBasisVectors[(r - 1) + 2 * (d - 1)]; };
  std::vector<real> planarSourcePoints(static_cast<std::size_t>(pc) * 2), planarUnitVectors(static_cast<std::size_t>(pc) * 2);
  auto PSP = [&](int i, int d) -> real& { return planarSourcePoints[(i - 1) + pc * (d - 1)]; };
  auto PUV = [&](int i, int d) -> real& { return planarUnitVectors[(i - 1) + pc * (d - 1)]; };
  const std::size_t ms = static_cast<std::size_t>(matrixSize);
  std::vector<real> matrix(ms * ms, 0.0), rhs(ms * 2, 0.0), coeffs(ms * 2, 0.0);  // :1106-1108
  auto MAT = [&](int i, int j) -> real& { return matrix[(i - 1) + ms * (j - 1)]; };
  auto RHS = [&](int i, int d) -> real& { return rhs[(i - 1) + ms * (d - 1)]; };

  for (int i = 1; i <= pc; ++i) {  // :1110-1115, four 3-term SUMs
    for (int r = 1; r <= 2; ++r) {
      real s = 0.0;
      for (int d = 1; d <= 3; ++d) s = s + SP(i, d) * PB(r, d);
      PSP(i, r) = s;
      real u = 0.0;
      for (int d = 1; d <= 3; ++d) u = u + UV(i, d) * PB(r, d);
      PUV(i, r) = u;
    }
  }
  real planarDestinationPoint[2];  // :1116-1117
  for (int r = 1; r <= 2; ++r) {
    real s = 0.0;
    for (int d = 1; d <= 3; ++d) s = s + destinationPoint[d - 1] * PB(r, d);
    planarDestinationPoint[r - 1] = s;
  }

  std::vector<real> matrixWork(static_cast<std::size_t>(pc) * pc), rhsWork(static_cast<std::size_t>(pc) * 2);  // :1119-1121
  set_up_vector_dirichlet_rbf_matrix_and_rhs(pc, 2, planarSourcePoints.data(), planarUnitVectors.data(),
                                             planarDestinationPoint, alpha, matrixWork.data(), rhsWork.data());
  for (int j = 1; j <= pc; ++j)  // :1123-1124
    for (int i = 1; i <= pc; ++i) MAT(i, j) = matrixWork[(i - 1) + pc * (j - 1)];
  for (int d = 1; d <= 2; ++d)
    for (int i = 1; i <= pc; ++i) RHS(i, d) = rhsWork[(i - 1) + pc * (d - 1)];
  for (int i = 1; i <= pc; ++i) {  // :1126-1129
    MAT(i, pc + 1) = PUV(i, 1);
    MAT(i, pc + 2) = PUV(i, 2);
    MAT(pc + 1, i) = MAT(i, pc + 1);
    MAT(pc + 2, i) = MAT(i, pc + 2);
  }
  for (int i = 1; i <= 2; ++i) RHS(pc + i, i) = 1.0;  // :1130-1132

  std::vector<real> matrixCopy = matrix;  // :1135-1137
  std::vector<int> pivotIndices(ms);
  mpas_legs(matrix.data(), matrixSize, &rhs[0], &coeffs[0], pivotIndices.data());
  mpas_legs(matrixCopy.data(), matrixSize, &rhs[ms], &coeffs[ms], pivotIndices.data());

  for (int i = 1; i <= 3; ++i)  // :1140-1143
    for (int p = 1; p <= pc; ++p)
      coefficients[(p - 1) + pc * (i - 1)] = PB(1, i) * coeffs[p - 1] + PB(2, i) * coeffs[ms + (p - 1)];
}

// mpas_vector_reconstruction.F:61-191
void mpas_init_reconstruct(Mesh& m) {
  if (m.is_periodic)
    throw Error("mpas_init_reconstruct: the is_periodic branch (mpas_vector_reconstruction.F:141-144) is not ported");
  auto& cr = m.coeffs_reconstruct;
  std::fill(cr.data(), cr.data() + cr.size(), 0.0);  // :124
  const int nCells = m.nCellsSolve();                 // :117-121
  for (int iCell = 1; iCell <= nCells; ++iCell) {
    const int pointCount = m.nEdgesOnCell(iCell);
    const real cellCenter[3] = {m.xCell(iCell), m.yCell(iCell), m.zCell(iCell)};
    const std::size_t pc = static_cast<std::size_t>(pointCount);
    std::vector<real> loc(pc * 3), nrm(pc * 3), coeffs(pc * 3);
    for (int i = 1; i <= pointCount; ++i) {  // :139-151
      const int iEdge = m.edgesOnCell(i, iCell);
      loc[(i - 1) + pc * 0] = m.xEdge(iEdge);
      loc[(i - 1) + pc * 1] = m.yEdge(iEdge);
      loc[(i - 1) + pc * 2] = m.zEdge(iEdge);
      for (int d = 1; d <= 3; ++d) nrm[(i - 1) + pc * (d - 1)] = m.edgeNormalVectors(d, iEdge);
    }
    real alpha = 0.0;  // :153-158
    for (int i = 1; i <= pointCount; ++i) {
      real s = 0.0;
      for (int d = 1; d <= 3; ++d) {
        const real t = cellCenter[d - 1] - loc[(i - 1) + pc * (d - 1)];
        s = s + t * t;
      }
      const real r = std::sqrt(s);
      alpha = alpha + r;
    }
    alpha = alpha / static_cast<real>(pointCount);
    real tangentPlane[6];  // (2,3): :160-161
    for (int d = 1; d <= 3; ++d) {
      tangentPlane[0 + 2 * (d - 1)] = m.cellTangentPlane(d, 1, iCell);
      tangentPlane[1 + 2 * (d - 1)] = m.cellTangentPlane(d, 2, iCell);
    }
    mpas_rbf_interp_func_3D_plane_vec_const_dir_comp_coeffs(pointCount, loc.data(), nrm.data(), cellCenter, alpha,
                                                            tangentPlane, coeffs.data());
    for (int i = 1; i <= pointCount; ++i)  // :181-183
      for (int d = 1; d <= 3; ++d) cr(d, i, iCell) = coeffs[(i - 1) + pc * (d - 1)];
  }
}

}  // namespace cppmpas
