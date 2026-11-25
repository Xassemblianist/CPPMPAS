// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// Reader of MPAS grid and init files into the mesh pool and initial state as MPAS holds them in memory on one rank
// (src/core/mesh.hpp). What MPAS would read but CPPMPAS has no port for is an error.
#include <map>
#include <string>
#include <vector>

#include "core/mesh.hpp"
#include "io/netcdf_file.hpp"

namespace cppmpas {

// What the source file looked like, so a file can be written back with the same metadata (round trip).
struct SourceInfo {
  std::string path;
  int format = 0;
  std::vector<RawAttribute> global_atts;
  std::map<std::string, std::vector<RawAttribute>> var_atts;
  std::map<std::string, std::pair<std::size_t, bool>> dims;  // length, unlimited
  std::vector<std::string> vars_read;                         // in read order
};

struct MpasInput {
  Mesh mesh;
  State state;
  bool has_state = false;
  SourceInfo source;
  std::vector<std::string> warnings;  // MPAS's own warnings (mpas_bootstrapping.F), e.g. a missing attribute
};

// A mesh ("grid") file as the init_atmosphere core reads it: the horizontal mesh only.
MpasInput read_mpas_grid(const std::string& path);
// An init file as the atmosphere core's "input" stream reads it (Registry.xml:602-648). scalars: the active
// constituents of var_array scalars in Registry order (from the packages the namelist selects: qv for JW and MW,
// qv, qc, qr for SC with mp_kessler, DYCORE_MAP.md §0.3).
MpasInput read_mpas_init(const std::string& path, const std::vector<std::string>& scalars);

// Field access by MPAS name (for the writer and tests). Exactly one of r / i / c is set.
struct FieldRef {
  std::string name;
  FArray<real>* r = nullptr;
  FArray<int>* i = nullptr;
  std::string* c = nullptr;
  int scalar_slot = 0;  // > 0: constituent of State::scalars
};
FieldRef field_ref(MpasInput& in, const std::string& name);  // throws Error for a field CPPMPAS does not handle
std::vector<std::string> handled_mesh_fields();
std::vector<std::string> handled_init_fields();             // without scalars constituents

}  // namespace cppmpas
