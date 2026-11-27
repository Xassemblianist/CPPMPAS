// SPDX-License-Identifier: GPL-3.0-or-later
#include "io/registry_pools.hpp"

#include <sstream>

#include "core/error.hpp"

namespace cppmpas {

bool ActivePackages::active(const RegistryVar& rv) const {
  if (rv.packages.empty()) return true;
  std::stringstream ss(rv.packages);
  std::string p;
  while (std::getline(ss, p, ';'))
    if (on.count(p)) return true;
  return false;
}

void RegistryPools::set_mesh_dims(const Mesh& m) {
  dims["nCells"] = m.nCells;
  dims["nEdges"] = m.nEdges;
  dims["nVertices"] = m.nVertices;
  dims["maxEdges"] = m.maxEdges;
  dims["maxEdges2"] = m.maxEdges2;
  dims["vertexDegree"] = m.vertexDegree;
  dims["TWO"] = 2;
  dims["THREE"] = 3;
  dims["R3"] = 3;
  dims["FIFTEEN"] = 15;
  dims["TWENTYONE"] = 21;
  if (m.nVertLevels > 0) {
    dims["nVertLevels"] = m.nVertLevels;
    dims["nVertLevelsP1"] = m.nVertLevels + 1;
    dims["nVertLevelsP2"] = m.nVertLevels + 2;
  }
}

int RegistryPools::dim(const std::string& d, bool garbage) const {
  const auto it = dims.find(d);
  if (it == dims.end()) throw Error("RegistryPools: dimension '" + d + "' is not defined");
  const bool element = d == "nCells" || d == "nEdges" || d == "nVertices";
  return it->second + (garbage && element ? 1 : 0);
}

std::vector<int> RegistryPools::extents(const RegistryVar& rv) const {
  std::vector<int> ext;
  for (const auto& d : rv.space_dims()) ext.push_back(dim(d, true));
  if (ext.empty()) ext.push_back(1);
  return ext;
}

namespace {
const RegistryVar& active_var(const RegistryPools& p, const std::string& name, const char* type) {
  const RegistryVar& rv = registry_var(name);
  if (rv.type != type) throw Error("RegistryPools: '" + name + "' is " + rv.type + ", not " + type);
  if (!p.packages.active(rv))
    throw Error("RegistryPools: '" + name + "' belongs to the inactive package(s) " + rv.packages +
                " (MPAS does not allocate it)");
  return rv;
}
}  // namespace

FArray<real>& RegistryPools::R(const std::string& name) {
  auto it = r.find(name);
  if (it != r.end()) return it->second;
  const RegistryVar& rv = active_var(*this, name, "real");
  // src/tools/registry/gen_inc.c:615-630: an empty default is 0
  const real def = rv.default_value.empty() ? 0.0 : std::stod(rv.default_value);
  FArray<real>& a = r[name];
  a.allocate(extents(rv), def);
  return a;
}

FArray<int>& RegistryPools::I(const std::string& name) {
  auto it = i.find(name);
  if (it != i.end()) return it->second;
  const RegistryVar& rv = active_var(*this, name, "integer");
  const int def = rv.default_value.empty() ? 0 : std::stoi(rv.default_value);
  FArray<int>& a = i[name];
  a.allocate(extents(rv), def);
  return a;
}

std::string& RegistryPools::C(const std::string& name) {
  auto it = c.find(name);
  if (it != c.end()) return it->second;
  const RegistryVar& rv = active_var(*this, name, "text");
  std::string& s = c[name];
  s = rv.default_value;
  s.resize(64, ' ');
  return s;
}

}  // namespace cppmpas
