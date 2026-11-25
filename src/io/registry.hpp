// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// MPAS-Atmosphere v8.4.2 Registry metadata (src/core_atmosphere/Registry.xml), generated into
// registry_atm_table.inc by tools/gen_registry_table.py.
#include <string>
#include <utility>
#include <vector>

namespace cppmpas {

struct RegistryVar {
  std::string name, var_struct, type;            // type: "real", "integer", "text"
  std::vector<std::string> dims;                 // Fortran order, as in the Registry (may end with "Time")
  std::string array;                             // var_array name for constituents ("scalars"), else empty
  std::string default_value, packages;
  int time_levs = 1;
  std::vector<std::pair<std::string, std::string>> attrs;  // NetCDF attributes MPAS writes, in order
  bool has_time() const { return !dims.empty() && dims.back() == "Time"; }
  std::vector<std::string> space_dims() const;   // dims without Time
};

const RegistryVar* find_registry_var(const std::string& name);
const RegistryVar& registry_var(const std::string& name);  // throws Error for a name the Registry does not define
const std::vector<RegistryVar>& registry_vars();
// variables of the atmosphere "input" stream (Registry.xml:602-648, nested "invariant" stream expanded)
const std::vector<std::string>& input_stream_vars();

}  // namespace cppmpas
