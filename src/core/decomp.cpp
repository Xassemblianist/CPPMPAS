// SPDX-License-Identifier: GPL-3.0-or-later
#include "core/decomp.hpp"

#include <algorithm>
#include <string>
#include <unordered_set>

namespace cppmpas {

// Arguments are checked (negative counts, a count beyond its dimension, a search past the list): these routines are
// public, MPAS's take them as valid. Offsets into the (maxEdges, nCells) and (maxCells, nEdges) arrays are computed
// in 64 bits, as gfortran computes the offset of a 2-D array element: maxEdges * nCells may exceed 2^31.
namespace {
std::size_t at2(int j, int i, int dim1) {
  return static_cast<std::size_t>(j - 1) + static_cast<std::size_t>(dim1) * static_cast<std::size_t>(i - 1);
}
}  // namespace

std::vector<int> block_decomp_one_block_cell_list(int nCellsGlobal) {
  // mpas_block_decomp.F:304-315
  if (nCellsGlobal < 0) throw Error("block_decomp_one_block_cell_list: nCellsGlobal < 0");
  std::vector<int> local_cell_list(nCellsGlobal);
  for (int i = 1; i <= nCellsGlobal; ++i) local_cell_list[i - 1] = i;
  return local_cell_list;
}

std::vector<int> block_decomp_all_edges_in_block(int maxEdges, int nCells, const int* nEdgesOnCell,
                                                 const int* edgesOnCell) {
  // mpas_block_decomp.F:647-655: count the distinct IDs (the hash only answers membership)
  if (maxEdges < 0 || nCells < 0) throw Error("block_decomp_all_edges_in_block: negative dimension");
  for (int i = 1; i <= nCells; ++i)
    if (nEdgesOnCell[i - 1] < 0 || nEdgesOnCell[i - 1] > maxEdges)
      throw Error("block_decomp_all_edges_in_block: nEdgesOnCell(" + std::to_string(i) + ") outside 0..maxEdges");
  std::unordered_set<int> h;
  for (int i = 1; i <= nCells; ++i)
    for (int j = 1; j <= nEdgesOnCell[i - 1]; ++j) h.insert(edgesOnCell[at2(j, i, maxEdges)]);
  const int nEdges = static_cast<int>(h.size());
  std::vector<int> edgeList(nEdges);
  h.clear();
  // mpas_block_decomp.F:662-676: first-seen order
  int k = 0;
  for (int i = 1; i <= nCells; ++i) {
    for (int j = 1; j <= nEdgesOnCell[i - 1]; ++j) {
      const int e = edgesOnCell[at2(j, i, maxEdges)];
      if (!h.count(e)) {
        k = k + 1;
        if (k > nEdges) throw Error("block_decomp_all_edges_in_block: Trying to add more edges than expected.");
        edgeList[k - 1] = e;
        h.insert(e);
      }
    }
  }
  if (k < nEdges) throw Error("block_decomp_all_edges_in_block: Listed fewer edges than expected.");
  return edgeList;
}

int block_decomp_partitioned_edge_list(int nCells, const int* cellIDList, int maxCells, int nEdges,
                                       const int* cellsOnEdge, int* edgeIDList) {
  if (nCells < 0 || maxCells < 0 || nEdges < 0) throw Error("block_decomp_partitioned_edge_list: negative dimension");
  std::unordered_set<int> h;  // mpas_block_decomp.F:581-585
  for (int i = 1; i <= nCells; ++i) h.insert(cellIDList[i - 1]);
  int lastEdge = 0;
  int ghostEdgeStart = nEdges + 1;
  const std::vector<int> edgeIDListLocal(edgeIDList, edgeIDList + nEdges);  // :591
  for (int i = 1; i <= nEdges; ++i) {  // mpas_block_decomp.F:593-612
    int j = 1;
    for (; j <= maxCells; ++j)
      if (cellsOnEdge[at2(j, i, maxCells)] != 0) break;
    if (j > maxCells)
      throw Error("Error in block_decomp_partitioned_edge_list: edge/vertex is not adjacent to any valid cells");
    if (h.count(cellsOnEdge[at2(j, i, maxCells)])) {
      lastEdge = lastEdge + 1;
      edgeIDList[lastEdge - 1] = edgeIDListLocal[i - 1];
    } else {
      ghostEdgeStart = ghostEdgeStart - 1;
      edgeIDList[ghostEdgeStart - 1] = edgeIDListLocal[i - 1];
    }
    if (ghostEdgeStart <= lastEdge)
      throw Error("block_decomp_partitioned_edge_list: Somehow we have more edges than we thought we should.");
  }
  if (ghostEdgeStart != lastEdge + 1)
    throw Error("block_decomp_partitioned_edge_list: Somehow we didn't have enough edges to fill edgeIDList.");
  return ghostEdgeStart;
}

int mpas_binary_search(const std::vector<std::pair<int, int>>& sorted, int n2, int key) {
  // mpas_sort.F:430-446 with n1 = 1; sorted[k-1].first is array(1,k)
  if (n2 < 0 || static_cast<std::size_t>(n2) > sorted.size())
    throw Error("mpas_binary_search: n2 outside 0..size of the list");
  // (l + u) / 2 in 64 bits: the same k as MPAS's default-integer sum whenever that does not overflow
  const auto mid = [](int l, int u) { return static_cast<int>((static_cast<long long>(l) + u) / 2); };
  int result = n2 + 1;
  int l = 1, u = n2;
  int k = mid(l, u);
  while (u >= l) {
    if (sorted[k - 1].first == key) {
      result = k;
      break;
    } else if (sorted[k - 1].first < key) {
      l = k + 1;
      k = mid(l, u);
    } else {
      u = k - 1;
      k = mid(l, u);
    }
  }
  return result;
}

void postread_reindex(FArray<int>& field, int innerDim, int outerDim, const FArray<int>& indexToID,
                      int indexSpaceDim) {
  // mpas_stream_manager.F:5619-5624: sortedID(1,i) = indexToID(i), sortedID(2,i) = i, mpas_quicksort (IDs are unique,
  // so the sort order is fully determined)
  if (innerDim < 0 || outerDim < 0 || indexSpaceDim < 0 || indexSpaceDim > indexToID.extent(0) ||
      innerDim > field.extent(0) || outerDim > field.extent(1))
    throw Error("postread_reindex: dimension outside its array");
  std::vector<std::pair<int, int>> sortedID(indexSpaceDim);
  for (int i = 1; i <= indexSpaceDim; ++i) sortedID[i - 1] = {indexToID(i), i};
  std::sort(sortedID.begin(), sortedID.end());
  for (int i = 1; i < indexSpaceDim; ++i)
    if (sortedID[i].first == sortedID[i - 1].first)
      throw Error("postread_reindex: duplicate global ID " + std::to_string(sortedID[i].first));
  // mpas_stream_manager.F:5627-5636
  for (int i = 1; i <= outerDim; ++i) {
    for (int j = 1; j <= innerDim; ++j) {
      const int k = mpas_binary_search(sortedID, indexSpaceDim, field(j, i));
      if (k <= indexSpaceDim) {
        field(j, i) = sortedID[k - 1].second;
      } else {
        field(j, i) = indexSpaceDim + 1;
      }
    }
  }
}

}  // namespace cppmpas
