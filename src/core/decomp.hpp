// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// MPAS's local element numbering on one MPI rank with one block, ported routine by routine. All lists hold MPAS
// global IDs (1-based); arrays are Fortran order as in MPAS.
#include <vector>

#include "core/field.hpp"

namespace cppmpas {

// src/framework/mpas_block_decomp.F:303-315 (total_blocks == 1): local_cell_list(i) = i.
std::vector<int> block_decomp_one_block_cell_list(int nCellsGlobal);

// src/framework/mpas_block_decomp.F:633-685 mpas_block_decomp_all_edges_in_block: the edges (or, called with
// verticesOnCell, the vertices) of the block in first-seen order over do i=1,nCells; do j=1,nEdgesOnCell(i).
// edgesOnCell: Fortran (maxEdges, nCells) of global IDs.
std::vector<int> block_decomp_all_edges_in_block(int maxEdges, int nCells, const int* nEdgesOnCell,
                                                 const int* edgesOnCell);

// src/framework/mpas_block_decomp.F:565-620 mpas_block_decomp_partitioned_edge_list: owned edges (first nonzero
// cell of cellsOnEdge is an owned cell) keep their order at the front, ghost edges are placed from the end
// backwards. cellsOnEdge: Fortran (maxCells, nEdges) of global cell IDs, edge i of the current edgeIDList.
// Returns ghostEdgeStart (nEdges+1 when every edge is owned).
int block_decomp_partitioned_edge_list(int nCells, const int* cellIDList, int maxCells, int nEdges,
                                       const int* cellsOnEdge, int* edgeIDList);

// src/framework/mpas_sort.F:421-448 mpas_binary_search over sorted (id, index) pairs: the position of key, or
// n2+1 when absent (n1 = 1).
int mpas_binary_search(const std::vector<std::pair<int, int>>& sorted, int n2, int key);

// src/framework/mpas_stream_manager.F:5413-5645 (postread_reindex, loop 5619-5636): every entry (j=1..innerDim, i=1..outerDim) of a
// connectivity field read with global IDs becomes the local index of that ID in indexToID(1:indexSpaceDim), or
// indexSpaceDim+1 when the ID is not local (0 included). The garbage column outerDim+1 is not touched.
void postread_reindex(FArray<int>& field, int innerDim, int outerDim, const FArray<int>& indexToID,
                      int indexSpaceDim);

}  // namespace cppmpas
