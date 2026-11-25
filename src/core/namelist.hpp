// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// The namelist options of namelist.atmosphere that the ported init-time code reads, with the Registry defaults
// (src/core_atmosphere/Registry.xml, nml_option default_value) for options a namelist leaves out.
// Only the options listed in kInitOptions are taken; the reading of the complete namelist (every option either
// ported or rejected, CONSTITUTION.md 2.4) is a separate unit. Options whose value selects a code path that is not
// ported are errors here (config_do_restart = true).
#include <map>
#include <string>

#include "core/precision.hpp"

namespace cppmpas {

struct InitConfig {
  real config_coef_3rd_order = 0.25;              // Registry.xml:244
  bool config_h_ScaleWithMesh = true;              // Registry.xml:274
  real config_zd = 22000.0;                        // Registry.xml:291
  real config_xnutr = 0.2;                         // Registry.xml:296
  real config_epssm_minimum = 0.1;                 // Registry.xml:326
  real config_epssm_maximum = 0.5;                 // Registry.xml:331
  real config_epssm_transition_bottom_z = 30000.;  // Registry.xml:336
  real config_epssm_transition_top_z = 50000.;     // Registry.xml:341
  bool config_do_restart = false;                  // Registry.xml:395
  bool config_do_DAcycling = false;                // Registry.xml:400
};

// Raw "key = value" pairs of a Fortran namelist file (keys lower-case; values as written, quotes removed).
std::map<std::string, std::string> read_namelist_file(const std::string& path);
// InitConfig from a namelist.atmosphere file (empty path: the Registry defaults).
InitConfig read_init_config(const std::string& path);

}  // namespace cppmpas
