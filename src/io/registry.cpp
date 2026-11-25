// SPDX-License-Identifier: GPL-3.0-or-later
#include "io/registry.hpp"

#include <map>
#include <sstream>

#include "core/error.hpp"

namespace cppmpas {

namespace {

std::vector<std::string> split_ws(const std::string& s) {
  std::istringstream is(s);
  std::vector<std::string> out;
  std::string w;
  while (is >> w) out.push_back(w);
  return out;
}

std::vector<RegistryVar> build() {
  std::vector<RegistryVar> v;
  auto add = [&v](const char* name, const char* vs, const char* type, const char* dims, const char* array,
                  const char* def, const char* pkgs, int tl, const char* attrs) {
    RegistryVar r;
    r.name = name;
    r.var_struct = vs;
    r.type = type;
    r.dims = split_ws(dims);
    r.array = array;
    r.default_value = def;
    r.packages = pkgs;
    r.time_levs = tl;
    std::istringstream is(attrs);
    std::string line;
    while (std::getline(is, line)) {
      const auto eq = line.find('=');
      if (eq == std::string::npos) continue;
      r.attrs.emplace_back(line.substr(0, eq), line.substr(eq + 1));
    }
    v.push_back(std::move(r));
  };
#define CPPMPAS_REGISTRY_VAR(n, s, t, d, a, def, p, tl, at) add(n, s, t, d, a, def, p, tl, at);
#include "io/registry_atm_table.inc"
#undef CPPMPAS_REGISTRY_VAR
  return v;
}

std::vector<std::string> build_input() {
  std::vector<std::string> v;
#define CPPMPAS_INPUT_STREAM_VAR(n) v.emplace_back(n);
#include "io/registry_atm_table.inc"
#undef CPPMPAS_INPUT_STREAM_VAR
  return v;
}

const std::map<std::string, std::size_t>& index() {
  static const std::map<std::string, std::size_t> idx = [] {
    std::map<std::string, std::size_t> m;
    const auto& v = registry_vars();
    for (std::size_t i = 0; i < v.size(); ++i) m[v[i].name] = i;
    return m;
  }();
  return idx;
}

}  // namespace

std::vector<std::string> RegistryVar::space_dims() const {
  std::vector<std::string> d = dims;
  if (has_time()) d.pop_back();
  return d;
}

const std::vector<RegistryVar>& registry_vars() {
  static const std::vector<RegistryVar> v = build();
  return v;
}

const std::vector<std::string>& input_stream_vars() {
  static const std::vector<std::string> v = build_input();
  return v;
}

const RegistryVar* find_registry_var(const std::string& name) {
  const auto& idx = index();
  const auto it = idx.find(name);
  return it == idx.end() ? nullptr : &registry_vars()[it->second];
}

const RegistryVar& registry_var(const std::string& name) {
  const RegistryVar* r = find_registry_var(name);
  if (!r) throw Error("variable '" + name + "' is not defined in MPAS-Atmosphere's Registry");
  return *r;
}

}  // namespace cppmpas
