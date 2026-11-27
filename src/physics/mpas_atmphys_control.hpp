// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// MPAS-Atmosphere's physics configuration: the namelist options of the &physics record the physics reads
// (Registry.xml:2200-2380), physics_namelist_check (physics/mpas_atmphys_control.F:92-332: the schemes of
// config_physics_suite, the checks, moist_physics), the physics packages (physics/mpas_atmphys_packages.F:84-176, and
// the 'sfclayer' package of Registry.xml:488-490) and the refusal of every configuration whose path is not ported.
#include <string>
#include <vector>

#include "core/namelist.hpp"
#include "core/precision.hpp"
#include "io/registry_pools.hpp"

namespace cppmpas::phys {

struct PhysicsConfig {
  // &nhyd_model / &restart / &decomposition values the physics reads
  real config_dt = 0, config_len_disp = 0;
  bool config_do_restart = false;
  std::string config_start_time;
  std::string config_les_model;
  // &physics
  std::string input_soil_data;
  int input_soil_temperature_lag = 0, num_soil_layers = 0, months = 0, noznlev = 0, naerlev = 0, camdim1 = 0;
  bool config_frac_seaice = true, config_sfc_albedo = true, config_sfc_snowalbedo = true;
  bool config_sst_update = false, config_sstdiurn_update = false, config_deepsoiltemp_update = false;
  bool config_o3climatology = true, config_microp_re = false, config_ysu_pblmix = false, config_urban_physics = false;
  int config_n_microp = 1;
  real config_microphysics_top = 45000.0;
  std::string config_radtlw_interval, config_radtsw_interval, config_conv_interval, config_pbl_interval;
  std::string config_camrad_abs_update, config_greeness_update, config_bucket_update;
  std::string config_physics_suite, config_microp_scheme, config_convection_scheme, config_lsm_scheme,
      config_pbl_scheme, config_gwdo_scheme, config_radt_cld_scheme, config_radt_lw_scheme, config_radt_sw_scheme,
      config_sfclayer_scheme, config_radt_cld_overlap, config_radt_cld_dcorrlen;
  bool config_ngw_scheme = false, config_ugwp_diags = false, config_oml1d = false;
  real config_bucket_radt = 0, config_bucket_rainc = 0, config_bucket_rainnc = 0;
  std::string config_sounding_interval;
  // set by physics_namelist_check; the namelist's config_sfclayer_scheme before it (the 'sfclayer' package)
  bool moist_physics = false;
  std::string config_sfclayer_scheme_namelist;
};

// The options above from namelist.atmosphere (effective_namelist: every value already checked against the verified
// references), config_len_disp as atm_core_init left it, then physics_namelist_check. log: MPAS's log lines.
PhysicsConfig read_physics_config(const std::string& namelist_path, const InitConfig& init_cfg,
                                  std::vector<std::string>& log);
// mpas_atmphys_control.F:92-332
void physics_namelist_check(PhysicsConfig& c, std::vector<std::string>& log);
// mpas_atm_core_interface.F:atm_setup_packages for the physics and the packages of the run
ActivePackages physics_packages(const PhysicsConfig& c);
// Refuses (Error naming the option) every physics configuration whose path is not ported and verified: CPPMPAS
// runs config_physics_suite = 'none' (no physics driver) and 'mesoscale_reference' with RG's options.
void check_physics_supported(const PhysicsConfig& c);

}  // namespace cppmpas::phys
