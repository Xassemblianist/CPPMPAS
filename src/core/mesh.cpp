// SPDX-License-Identifier: GPL-3.0-or-later
#include "core/mesh.hpp"

#include "core/decomp.hpp"

namespace cppmpas {

const char* element_dim_name(Element e) {
  switch (e) {
    case Element::Cell: return "nCells";
    case Element::Edge: return "nEdges";
    case Element::Vertex: return "nVertices";
  }
  return "?";
}

bool ElementOrder::identity() const {
  for (int i = 1; i <= n(); ++i)
    if (local_to_file[i] != i) return false;
  return true;
}

const ElementOrder& Mesh::order(Element e) const {
  switch (e) {
    case Element::Cell: return cell_order;
    case Element::Edge: return edge_order;
    case Element::Vertex: return vertex_order;
  }
  throw Error("Mesh::order: bad element");
}

namespace {

ElementOrder make_order(const std::vector<int>& list, int nFile, const char* what) {
  const int n = static_cast<int>(list.size());
  if (n != nFile)
    throw Error(std::string("mesh: ") + std::to_string(nFile - n) + " of the " + std::to_string(nFile) + " " + what +
                " in the file are not on any cell; MPAS would leave them out of the block (mesh not handled)");
  ElementOrder o;
  o.local_to_file.assign(n + 1, 0);
  o.file_to_local.assign(n + 1, 0);
  for (int i = 1; i <= n; ++i) {
    const int g = list[i - 1];
    if (g < 1 || g > nFile || o.file_to_local[g] != 0)
      throw Error(std::string("mesh: invalid or repeated ") + what + " ID " + std::to_string(g));
    o.local_to_file[i] = g;
    o.file_to_local[g] = i;
  }
  return o;
}

void set_index_array(FArray<int>& a, const ElementOrder& o) {
  // indexTo*ID(1:n) = local list, garbage element 0 (mpas_block_creator.F:1068-1069, 1094-1095, 1121-1122)
  a.allocate({o.n() + 1}, 0);
  for (int i = 1; i <= o.n(); ++i) a(i) = o.local_to_file[i];
}

}  // namespace

void build_one_rank_order(Mesh& m, const std::vector<int>& nEdgesOnCell, const std::vector<int>& edgesOnCell,
                          const std::vector<int>& verticesOnCell, const std::vector<int>& cellsOnEdge,
                          const std::vector<int>& cellsOnVertex) {
  // cells: one block on one rank -> identity (mpas_block_decomp.F:303-315)
  const std::vector<int> cellList = block_decomp_one_block_cell_list(m.nCells);
  m.cell_order = make_order(cellList, m.nCells, "cells");

  // mpas_bootstrapping.F:261-266: mpas_block_creator_build_0_and_1halo_edge_fields once for edges (edgesOnCell,
  // cellsOnEdge) and once for vertices (verticesOnCell, cellsOnVertex). The cell fields of the block are in local
  // cell order, which is file order here.
  auto build = [&](const std::vector<int>& xOnCell, const std::vector<int>& cellsOnX, int degree, int nFile,
                   const char* what) {
    // mpas_block_creator.F:318 -> mpas_block_decomp_all_edges_in_block
    std::vector<int> list =
        block_decomp_all_edges_in_block(m.maxEdges, m.nCells, nEdgesOnCell.data(), xOnCell.data());
    const int nLocal = static_cast<int>(list.size());
    // cellsOnEdge of the block, edge i = list(i) (mpas_dmpar_alltoall_field, mpas_block_creator.F:366)
    std::vector<int> cellsOnLocal(static_cast<std::size_t>(degree) * nLocal);
    for (int i = 1; i <= nLocal; ++i) {
      const int g = list[i - 1];
      if (g < 1 || g > nFile)
        throw Error(std::string("mesh: ") + what + " ID " + std::to_string(g) + " out of range 1.." +
                    std::to_string(nFile));
      for (int j = 1; j <= degree; ++j)
        cellsOnLocal[static_cast<std::size_t>(j - 1) + static_cast<std::size_t>(degree) * (i - 1)] =
            cellsOnX[static_cast<std::size_t>(j - 1) + static_cast<std::size_t>(degree) * (g - 1)];
    }
    // mpas_block_creator.F:391-393 -> mpas_block_decomp_partitioned_edge_list
    const int ghostStart = block_decomp_partitioned_edge_list(m.nCells, cellList.data(), degree, nLocal,
                                                              cellsOnLocal.data(), list.data());
    if (ghostStart != nLocal + 1)
      throw Error(std::string("mesh: ") + what + " owned by no local cell on one rank (mesh not handled)");
    return make_order(list, nFile, what);
  };
  m.edge_order = build(edgesOnCell, cellsOnEdge, 2, m.nEdges, "edges");
  m.vertex_order = build(verticesOnCell, cellsOnVertex, m.vertexDegree, m.nVertices, "vertices");

  set_index_array(m.indexToCellID, m.cell_order);
  set_index_array(m.indexToEdgeID, m.edge_order);
  set_index_array(m.indexToVertexID, m.vertex_order);
}

}  // namespace cppmpas
