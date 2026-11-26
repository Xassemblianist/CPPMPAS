// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// Writers of MPAS-format files. Elements are written in file (global ID) order and connectivity as global IDs, as
// MPAS writes history and restart streams (blockWrite = .false., mpas_stream_manager.F:2841, 2861):
//  * the I/O decomposition of a field on nCells / nEdges / nVertices is indexTo*ID(1:n*Solve)
//    (src/framework/add_field_indices.inc:19-29), so local element i goes to file record indexTo*ID(i);
//  * connectivity is converted back to global IDs by prewrite_reindex (mpas_stream_manager.F:4894-5210), which
//    writes 0 in the unused slots j > nEdgesOnCell(i) / nEdgesOnEdge(i) (5053, 5074, 5095, 5152) and maps a garbage
//    index n+1 to indexTo*ID(n+1) = 0.
#include <string>
#include <vector>

#include "io/mpas_reader.hpp"

namespace cppmpas {

// Writes the given fields (empty: every field that was read) with the source file's format, dimensions, global
// attributes and variable attributes: an MPAS file read and written back unchanged, except for the rules above.
void write_round_trip(MpasInput& in, const std::string& path, std::vector<std::string> fields = {});

struct HistoryOptions {
  // core % git_version: -DMPAS_GIT_VERSION=$(VER) (Makefile:1021-1026); the reference recipe builds with VER=v8.4.2
  std::string git_version = "v8.4.2";
  int nprocs = 1;
  // empty: a random 10-character [0-9a-z] id, as MPAS writes (mpas_stream_manager.F:4263, 4306-4310; random_id.c)
  std::string file_id = "";
};

// Writes an MPAS history file for the given fields: CDF-5 (SMIOL's default format, smiolf.F90:426-430, with no
// fformat passed by MPAS_io_open, mpas_io.F:420-421), Registry dimensions, types and attributes
// (gen_inc.c:43-58), xtime as char(Time, StrLen = 64), global attributes of mpas_subdriver.F:513-530 and
// mpas_stream_manager.F:4306-4310. NOT yet written: the namelist options MPAS adds as global attributes
// (mpas_subdriver.F:532-557); they need the namelist port.
void write_history(MpasInput& in, const std::string& path, const std::vector<std::string>& fields,
                   const HistoryOptions& opt = {});

}  // namespace cppmpas
