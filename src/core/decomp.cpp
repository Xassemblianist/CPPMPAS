// SPDX-License-Identifier: GPL-3.0-or-later
#include "core/decomp.hpp"

#include <algorithm>
#include <unordered_set>

namespace cppmpas {

std::vector<int> block_decomp_one_block_cell_list(int nCellsGlobal) {
  // mpas_block_decomp.F:304-315
  std::vector<int> local_cell_list(nCellsGlobal);
  for (int i = 1; i <= nCellsGlobal; ++i) local_cell_list[i - 1] = i;
  return local_cell_list;
}

std::vector<int> block_decomp_all_edges_in_block(int maxEdges, int nCells, const int* nEdgesOnCell,
                                                 const int* edgesOnCell) {
  // mpas_block_decomp.F:647-655: count the distinct IDs (the hash only answers membership)
  std::unordered_set<int> h;
  for (int i = 1; i <= nCells; ++i)
    for (int j = 1; j <= nEdgesOnCell[i - 1]; ++j) h.insert(edgesOnCell[(j - 1) + maxEdges * (i - 1)]);
  const int nEdges = static_cast<int>(h.size());
  std::vector<int> edgeList(nEdges);
  h.clear();
  // mpas_block_decomp.F:662-676: first-seen order
  int k = 0;
  for (int i = 1; i <= nCells; ++i) {
    for (int j = 1; j <= nEdgesOnCell[i - 1]; ++j) {
      const int e = edgesOnCell[(j - 1) + maxEdges * (i - 1)];
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
  std::unordered_set<int> h;  // mpas_block_decomp.F:581-585
  for (int i = 1; i <= nCells; ++i) h.insert(cellIDList[i - 1]);
  int lastEdge = 0;
  int ghostEdgeStart = nEdges + 1;
  const std::vector<int> edgeIDListLocal(edgeIDList, edgeIDList + nEdges);  // :591
  for (int i = 1; i <= nEdges; ++i) {  // mpas_block_decomp.F:593-612
    int j = 1;
    for (; j <= maxCells; ++j)
      if (cellsOnEdge[(j - 1) + maxCells * (i - 1)] != 0) break;
    if (j > maxCells)
      throw Error("Error in block_decomp_partitioned_edge_list: edge/vertex is not adjacent to any valid cells");
    if (h.count(cellsOnEdge[(j - 1) + maxCells * (i - 1)])) {
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
  int result = n2 + 1;
  int l = 1, u = n2;
  int k = (l + u) / 2;
  while (u >= l) {
    if (sorted[k - 1].first == key) {
      result = k;
      break;
    } else if (sorted[k - 1].first < key) {
      l = k + 1;
      k = (l + u) / 2;
    } else {
      u = k - 1;
      k = (l + u) / 2;
    }
  }
  return result;
}

void postread_reindex(FArray<int>& field, int innerDim, int outerDim, const FArray<int>& indexToID,
                      int indexSpaceDim) {
  // mpas_stream_manager.F:5619-5624: sortedID(1,i) = indexToID(i), sortedID(2,i) = i, mpas_quicksort (IDs are unique,
  // so the sort order is fully determined)
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
