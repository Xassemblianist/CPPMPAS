// SPDX-License-Identifier: GPL-3.0-or-later
// Unit test of the ports in src/core/decomp.cpp on a hand-made 3-cell strip with known answers.
//
//   cell 1 | cell 2 | cell 3      edgesOnCell (global IDs): 1: (4, 1, 2)  2: (2, 5, 3)  3: (3, 6, 7)
//
// first-seen edge order: 4 1 2 5 3 6 7 (mpas_block_decomp.F:662-676). With owned cells {2}, the owned edges keep
// their order at the front and the others go to the end backwards (mpas_block_decomp.F:596-609).
#include <cstdio>
#include <vector>

#include "core/decomp.hpp"

using namespace cppmpas;

static int fails = 0;
static void expect(bool ok, const char* what) {
  if (!ok) {
    std::printf("FAIL %s\n", what);
    ++fails;
  }
}

int main() {
  const int maxEdges = 3, nCells = 3;
  const std::vector<int> nEdgesOnCell = {3, 3, 3};
  const std::vector<int> edgesOnCell = {4, 1, 2, 2, 5, 3, 3, 6, 7};
  std::vector<int> list = block_decomp_all_edges_in_block(maxEdges, nCells, nEdgesOnCell.data(), edgesOnCell.data());
  expect(list == std::vector<int>({4, 1, 2, 5, 3, 6, 7}), "first-seen edge order");

  // cellsOnEdge (2, nEdges) for the edges of list, in list order: edge 4:(1,0) 1:(1,0) 2:(1,2) 5:(2,0) 3:(2,3)
  // 6:(3,0) 7:(3,0)
  const std::vector<int> cellsOnEdge = {1, 0, 1, 0, 1, 2, 2, 0, 2, 3, 3, 0, 3, 0};
  std::vector<int> one_rank = list;
  const std::vector<int> all = block_decomp_one_block_cell_list(3);
  expect(all == std::vector<int>({1, 2, 3}), "one-block cell list");
  int ghost = block_decomp_partitioned_edge_list(3, all.data(), 2, 7, cellsOnEdge.data(), one_rank.data());
  expect(ghost == 8 && one_rank == list, "one rank: every edge owned, order kept");
  const std::vector<int> owned = {2};
  std::vector<int> part = list;
  ghost = block_decomp_partitioned_edge_list(1, owned.data(), 2, 7, cellsOnEdge.data(), part.data());
  // owned (first nonzero cell == 2): 5, 3; ghosts from the end: 4 -> slot 7, 1 -> 6, 2 -> 5, 6 -> 4, 7 -> 3
  expect(ghost == 3 && part == std::vector<int>({5, 3, 7, 6, 2, 1, 4}), "partitioned edge list with ghosts");

  // postread_reindex: global IDs -> local indices, absent (0, 99) -> n+1, garbage column untouched
  FArray<int> idx({8}, 0);
  for (int i = 1; i <= 7; ++i) idx(i) = list[i - 1];
  FArray<int> conn({3, 4}, -5);  // (maxEdges, nCells+1)
  const int vals[9] = {4, 1, 2, 2, 5, 3, 0, 99, 7};
  for (int k = 0; k < 9; ++k) conn(k % 3 + 1, k / 3 + 1) = vals[k];
  postread_reindex(conn, 3, 3, idx, 7);
  const int want[9] = {1, 2, 3, 3, 4, 5, 8, 8, 7};
  bool ok = true;
  for (int k = 0; k < 9; ++k) ok &= conn(k % 3 + 1, k / 3 + 1) == want[k];
  ok &= conn(1, 4) == -5 && conn(2, 4) == -5 && conn(3, 4) == -5;
  expect(ok, "postread_reindex");

  // only slots 1..nEdgesOnCell are walked (mpas_block_decomp.F:649, 664): cell 2 has 2 edges and junk in slot 3
  const std::vector<int> nEdgesShort = {3, 2};
  const std::vector<int> edgesJunk = {4, 1, 2, 2, 5, 99};
  expect(block_decomp_all_edges_in_block(3, 2, nEdgesShort.data(), edgesJunk.data()) == std::vector<int>({4, 1, 2, 5}),
         "unused slots are not walked");

  std::vector<std::pair<int, int>> sorted = {{2, 3}, {5, 1}, {9, 2}};
  expect(mpas_binary_search(sorted, 3, 5) == 2 && mpas_binary_search(sorted, 3, 4) == 4, "mpas_binary_search");
  expect(mpas_binary_search(sorted, 3, 2) == 1 && mpas_binary_search(sorted, 3, 9) == 3, "binary search: first, last");
  expect(mpas_binary_search(sorted, 3, 1) == 4 && mpas_binary_search(sorted, 3, 10) == 4, "binary search: below, above");
  const std::vector<std::pair<int, int>> one = {{5, 1}}, none;
  expect(mpas_binary_search(one, 1, 5) == 1 && mpas_binary_search(one, 1, 4) == 2 && mpas_binary_search(none, 0, 4) == 1,
         "binary search: n = 1 and n = 0");
  if (fails == 0) std::printf("OK decomp unit test\n");
  return fails ? 1 : 0;
}
