// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// The Registry pools beyond the mesh, state and diag fields CPPMPAS holds in its own structs (core/mesh.hpp): the
// var_structs of the physics (sfc_input, diag_physics, tend_physics, atm_input, Registry.xml:2400-3990) and the diag /
// tend fields the physics coupling reads (relhum, tend_u_phys, ...). MPAS allocates every field of a var_struct whose
// packages are active, with the Registry dimensions plus one garbage element on nCells / nEdges / nVertices, and fills
// it with the Registry default value (mpas_block_creator.F:1634-1680; gen_inc.c:615-630 for an empty default: 0);
// a field of an inactive package is not allocated (mpas_pool_get_array returns a null pointer). Here a field is
// allocated the first time it is asked for, which gives the same values as long as nothing else has touched it.
// Asking for a field whose packages are inactive is an Error (MPAS would dereference a null pointer).
#include <map>
#include <set>
#include <string>
#include <vector>

#include "core/field.hpp"
#include "core/mesh.hpp"
#include "core/precision.hpp"
#include "io/registry.hpp"

namespace cppmpas {

// The packages of Registry.xml:466-491 that a run activates (mpas_atm_core_interface.F:atm_setup_packages and
// physics/mpas_atmphys_packages.F); every other package is inactive.
struct ActivePackages {
  std::set<std::string> on;
  bool active(const RegistryVar& rv) const;  // no packages attribute, or one of them active
};

struct RegistryPools {
  ActivePackages packages;
  // dimension lengths: the mesh's (set_mesh_dims), those read from the input file (nSoilLevels, nMonths) and those
  // defined from the namelist (nLags, nOznLevels, nAerLevels, cam_dim1; Registry.xml:37-47)
  std::map<std::string, int> dims;
  void set_mesh_dims(const Mesh& m);
  std::map<std::string, FArray<real>> r;
  std::map<std::string, FArray<int>> i;
  std::map<std::string, std::string> c;  // text fields: StrLen (64) characters, blank-padded
  std::vector<std::string> read_from_input;  // the fields the input stream read, in read order

  int dim(const std::string& d, bool garbage) const;  // memory extent (garbage: +1 on element dimensions)
  std::vector<int> extents(const RegistryVar& rv) const;
  bool has(const std::string& name) const { return r.count(name) || i.count(name) || c.count(name); }
  FArray<real>& R(const std::string& name);
  FArray<int>& I(const std::string& name);
  std::string& C(const std::string& name);
  FView<real> v(const std::string& name) { return R(name).view(); }
  FView<int> vi(const std::string& name) { return I(name).view(); }
};

}  // namespace cppmpas
