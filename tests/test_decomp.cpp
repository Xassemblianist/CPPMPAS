// SPDX-License-Identifier: GPL-3.0-or-later
// Unit test of the ports in src/core/decomp.cpp on a hand-made 3-cell strip with known answers.
//
//   cell 1 | cell 2 | cell 3      edgesOnCell (global IDs): 1: (4, 1, 2)  2: (2, 5, 3)  3: (3, 6, 7)
//
// first-seen edge order: 4 1 2 5 3 6 7 (mpas_block_decomp.F:662-676). With owned cells {2}, the owned edges keep
// their order at the front and the others go to the end backwards (mpas_block_decomp.F:596-609).
#include <sys/mman.h>

#include <cstdio>
#include <string>
#include <vector>

#include "core/decomp.hpp"
#include "core/error.hpp"

using namespace cppmpas;

static int fails = 0;
static void expect(bool ok, const char* what) {
  if (!ok) {
    std::printf("FAIL %s\n", what);
    ++fails;
  }
}
template <class F> static bool throws_with(F f, const std::string& text) {
  try {
    f();
  } catch (const Error& e) {
    return std::string(e.what()).find(text) != std::string::npos;
  }
  return false;
}
// n ints of address space, untouched pages not backed by memory (a 2^31-element array costs a few pages here)
static int* sparse_ints(std::size_t n) {
  void* p = mmap(nullptr, n * sizeof(int), PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS | MAP_NORESERVE, -1, 0);
  return p == MAP_FAILED ? nullptr : static_cast<int*>(p);
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

  // offsets past 2^31 elements (maxEdges * nCells > INT_MAX): gfortran indexes edgesOnCell(j, i) in 64 bits
  {
    const int big = 1 << 30;  // maxEdges (or maxCells) = 2^30: element (1, 3) is at offset 2^31
    const std::size_t n = (std::size_t(1) << 31) + 1;
    int* a = sparse_ints(n);
    expect(a != nullptr, "mmap of 2^31 ints");
    if (a) {
      a[std::size_t(1) << 31] = 42;
      const std::vector<int> nE = {0, 0, 1};
      expect(block_decomp_all_edges_in_block(big, 3, nE.data(), a) == std::vector<int>({42}),
             "edgesOnCell(1, 3) at offset 2^31");
      a[0] = 1;  // cellsOnEdge (maxCells = 2^30, 3 edges): the first cell of every edge is cell 1
      a[std::size_t(1) << 30] = 1;
      a[std::size_t(1) << 31] = 1;
      std::vector<int> ids = {7, 8, 9};
      const std::vector<int> cells = {1};
      expect(block_decomp_partitioned_edge_list(1, cells.data(), big, 3, a, ids.data()) == 4 &&
                 ids == std::vector<int>({7, 8, 9}),
             "cellsOnEdge(1, 3) at offset 2^31");
      munmap(a, n * sizeof(int));
    }
  }

  // arguments the public routines refuse
  expect(throws_with([] { (void)block_decomp_one_block_cell_list(-1); }, "nCellsGlobal < 0"), "cell list of -1 cells");
  {
    const std::vector<int> nE = {3, 4};
    expect(throws_with([&] { (void)block_decomp_all_edges_in_block(3, 2, nE.data(), edgesOnCell.data()); },
                       "nEdgesOnCell(2) outside 0..maxEdges"),
           "nEdgesOnCell > maxEdges");
    const std::vector<int> nN = {3, -1};
    expect(throws_with([&] { (void)block_decomp_all_edges_in_block(3, 2, nN.data(), edgesOnCell.data()); },
                       "nEdgesOnCell(2) outside 0..maxEdges"),
           "nEdgesOnCell < 0");
    std::vector<int> l2 = list;
    expect(throws_with([&] { (void)block_decomp_partitioned_edge_list(3, all.data(), 2, -1, cellsOnEdge.data(), l2.data()); },
                       "negative dimension"),
           "partitioned edge list of -1 edges");
  }
  expect(throws_with([&] { (void)mpas_binary_search(sorted, 4, 5); }, "n2 outside 0..size"), "binary search past the list");
  expect(throws_with([&] { (void)mpas_binary_search(sorted, -1, 5); }, "n2 outside 0..size"), "binary search, n2 < 0");
  {
    FArray<int> c2({3, 4}, 0);
    expect(throws_with([&] { postread_reindex(c2, 3, 3, idx, 9); }, "dimension outside its array"),
           "postread_reindex: indexSpaceDim past indexToID");
    expect(throws_with([&] { postread_reindex(c2, 4, 3, idx, 7); }, "dimension outside its array"),
           "postread_reindex: innerDim past the field");
  }

  // FArray: an element count that overflows size_t, a negative dimension index
  expect(throws_with([] { FArray<char> x({1 << 30, 1 << 30, 1 << 30, 1 << 30}); }, "number of elements overflows"),
         "FArray: 2^120 elements");
  {
    FArray<int> x({2, 3});
    expect(throws_with([&] { (void)x.extent(-1); }, "negative dimension index"), "FArray::extent(-1)");
  }

  if (fails == 0) std::printf("OK decomp unit test\n");
  return fails ? 1 : 0;
}
