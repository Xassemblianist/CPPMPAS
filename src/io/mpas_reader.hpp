// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// Reader of MPAS grid and init files into the mesh pool and initial state as MPAS holds them in memory on one rank
// (src/core/mesh.hpp). What MPAS would read but CPPMPAS has no port for is an error.
#include <map>
#include <string>
#include <vector>

#include "core/mesh.hpp"
#include "io/netcdf_file.hpp"
#include "io/registry_pools.hpp"

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
  Diag diag;
  bool has_state = false;
  RegistryPools pools;  // the fields of the other var_structs (read_input_pools), by name
  // fields a caller holds elsewhere (the dycore's diag / tend / tend_physics arrays of a run), by Registry name
  std::map<std::string, FArray<real>*> held;
  SourceInfo source;
  std::vector<std::string> warnings;  // MPAS's own warnings (mpas_bootstrapping.F), e.g. a missing attribute
};

// A mesh ("grid") file as the init_atmosphere core reads it: the horizontal mesh only.
MpasInput read_mpas_grid(const std::string& path);
// An init file as the atmosphere core's "input" stream reads it (Registry.xml:602-648). scalars: the active
// constituents of var_array scalars in Registry order (from the packages the namelist selects: qv for JW and MW,
// qv, qc, qr for SC with mp_kessler, Registry.xml:1701-1712). qv is required; a name given twice is an error.
MpasInput read_mpas_init(const std::string& path, const std::vector<std::string>& scalars);

// The fields of the input stream that belong to the Registry pools (sfc_input, diag_physics, diag's relhum, ...;
// Registry.xml:567-648) and to the active packages, read into in.pools with MPAS's memory extents and local order.
// A field of an inactive package is skipped ("Requested field ... is deactivated due to packages",
// mpas_bootstrapping.F); a field the file lacks keeps its Registry default with MPAS's warning
// "Variable <name> not in input file." (mpas_io_streams.F). dims: the namelist-defined dimensions (nLags,
// nOznLevels, nAerLevels, cam_dim1); nSoilLevels and nMonths come from the file when it has them.
void read_input_pools(const std::string& path, MpasInput& in, const ActivePackages& packages,
                      const std::map<std::string, int>& dims);
// The same for a restart run (mpas_atm_core.F:129-149): the pool fields of the restart stream (Registry.xml:650-
// 925: the physics tendencies, the surface and soil state, relhum, the mesh field qv_init, ...), read
// after read_mpas_restart. held: restart-stream fields the caller reads into its own arrays (the dycore's ru_p, rw_p,
// rt_diabatic_tend, rthdynten); the fields read_mpas_restart read are skipped too. A pool field the restart stream
// does not hold keeps its Registry default, as in MPAS, which reads only the restart stream on a restart.
void read_restart_pools(const std::string& path, MpasInput& in, const ActivePackages& packages,
                        const std::map<std::string, int>& dims, const std::vector<std::string>& held);
// Field access by MPAS name (for the writer and tests). Exactly one of r / i / c is set.
struct FieldRef {
  std::string name;
  FArray<real>* r = nullptr;
  FArray<int>* i = nullptr;
  std::string* c = nullptr;
  int scalar_slot = 0;  // > 0: constituent of State::scalars
};
FieldRef field_ref(MpasInput& in, const std::string& name);  // throws Error for a field CPPMPAS does not handle
// Reads the Registry field `name` of a file opened by the reader into a with MPAS's memory extents and local order
// (the Registry default in the garbage element), checking its type and dimensions; the record is the file's only one.
void read_registry_field(const NcFile& f, const Mesh& m, const std::string& name, FArray<real>& a);
// read_registry_field into the field field_ref(in, name) names, recorded as read (in.source).
void read_field_ref(const NcFile& f, MpasInput& in, const std::string& name);
std::vector<std::string> handled_mesh_fields();
std::vector<std::string> handled_init_fields();             // without scalars constituents

}  // namespace cppmpas
