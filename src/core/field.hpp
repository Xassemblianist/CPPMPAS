// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// MPAS arrays as MPAS holds them in memory: Fortran (column-major) order, 1-based indices, and the extents MPAS
// allocates, i.e. one garbage element on every nCells / nEdges / nVertices dimension
// (src/framework/mpas_domain_routines.F:251-265, src/framework/mpas_block_creator.F:1634-1635), so a(k, nCells+1)
// exists exactly as in MPAS. Fresh arrays hold the Registry default value everywhere, garbage element included
// (mpas_block_creator.F:1649, 1680).
#include <array>
#include <cstddef>
#include <string>
#include <vector>

#include "core/backend.hpp"
#include "core/error.hpp"

namespace cppmpas {

// Kernel view: trivially copyable, CPPMPAS_HD access, 1-based Fortran indexing. The offset is computed in 64 bits:
// on the NCAR x1.65536002 mesh (nEdges ~ 1.97e8) u(k, iEdge) with nVertLevels = 55 is ~1.1e10 elements in.
template <class T> struct FView {
  T* p = nullptr;
  int e0 = 1, e1 = 1, e2 = 1;  // extents of the first three dimensions
  CPPMPAS_HD std::ptrdiff_t offset(int i, int j = 1, int k = 1, int l = 1) const {
    using D = std::ptrdiff_t;
    return D(i - 1) + D(e0) * (D(j - 1) + D(e1) * (D(k - 1) + D(e2) * D(l - 1)));
  }
  CPPMPAS_HD T& operator()(int i) const { return p[i - 1]; }
  CPPMPAS_HD T& operator()(int i, int j) const { return p[offset(i, j)]; }
  CPPMPAS_HD T& operator()(int i, int j, int k) const { return p[offset(i, j, k)]; }
  CPPMPAS_HD T& operator()(int i, int j, int k, int l) const { return p[offset(i, j, k, l)]; }
};

template <class T> class FArray {
 public:
  FArray() = default;
  // extents in Fortran order (fastest first), garbage elements included by the caller
  explicit FArray(std::vector<int> extents, T fill = T{}) { allocate(std::move(extents), fill); }
  void allocate(std::vector<int> extents, T fill = T{}) {
    if (extents.empty() || extents.size() > 4) throw Error("FArray: rank must be 1..4");
    std::size_t n = 1;
    for (int e : extents) {
      if (e < 0) throw Error("FArray: negative extent");
      n *= static_cast<std::size_t>(e);
    }
    ext_ = std::move(extents);
    data_.assign(n, fill);
  }
  int rank() const { return static_cast<int>(ext_.size()); }
  int extent(int d) const { return d < rank() ? ext_[d] : 1; }  // d is 0-based
  const std::vector<int>& extents() const { return ext_; }
  std::size_t size() const { return data_.size(); }
  bool allocated() const { return !ext_.empty(); }
  T* data() { return data_.data(); }
  const T* data() const { return data_.data(); }
  T& operator()(int i) { return data_[i - 1]; }
  const T& operator()(int i) const { return data_[i - 1]; }
  T& operator()(int i, int j) { return data_[off(i, j, 1, 1)]; }
  const T& operator()(int i, int j) const { return data_[off(i, j, 1, 1)]; }
  T& operator()(int i, int j, int k) { return data_[off(i, j, k, 1)]; }
  const T& operator()(int i, int j, int k) const { return data_[off(i, j, k, 1)]; }
  T& operator()(int i, int j, int k, int l) { return data_[off(i, j, k, l)]; }
  const T& operator()(int i, int j, int k, int l) const { return data_[off(i, j, k, l)]; }
  FView<T> view() { return FView<T>{data_.data(), extent(0), extent(1), extent(2)}; }
  FView<const T> view() const { return FView<const T>{data_.data(), extent(0), extent(1), extent(2)}; }

 private:
  std::size_t off(int i, int j, int k, int l) const {
    const std::size_t e0 = extent(0), e1 = extent(1), e2 = extent(2);
    return static_cast<std::size_t>(i - 1) +
           e0 * (static_cast<std::size_t>(j - 1) + e1 * (static_cast<std::size_t>(k - 1) + e2 * (l - 1)));
  }
  std::vector<int> ext_;
  std::vector<T> data_;
};

}  // namespace cppmpas
