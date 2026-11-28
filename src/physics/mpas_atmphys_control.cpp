// SPDX-License-Identifier: GPL-3.0-or-later
#include "physics/mpas_atmphys_control.hpp"

#include <cstdlib>
#include <map>

#include "core/error.hpp"

namespace cppmpas::phys {

PhysicsConfig read_physics_config(const std::string& namelist_path, const InitConfig& init_cfg,
                                  std::vector<std::string>& log) {
  const auto kv = effective_namelist(namelist_path);
  auto get = [&](const char* k) -> const std::string& {
    const auto it = kv.find(k);
    if (it == kv.end()) throw Error(std::string("read_physics_config: option ") + k + " is not in the Registry table");
    return it->second;
  };
  auto r = [&](const char* k) { return std::strtod(get(k).c_str(), nullptr); };
  auto i = [&](const char* k) { return static_cast<int>(std::strtol(get(k).c_str(), nullptr, 10)); };
  auto l = [&](const char* k) { return get(k) == "true"; };
  PhysicsConfig c;
  c.config_dt = init_cfg.config_dt;
  c.config_len_disp = init_cfg.config_len_disp;
  c.config_do_restart = init_cfg.config_do_restart;
  c.config_start_time = init_cfg.config_start_time;
  c.config_les_model = get("config_les_model");
  c.input_soil_data = get("input_soil_data");
  c.input_soil_temperature_lag = i("input_soil_temperature_lag");
  c.num_soil_layers = i("num_soil_layers");
  c.months = i("months");
  c.noznlev = i("noznlev");
  c.naerlev = i("naerlev");
  c.camdim1 = i("camdim1");
  c.config_frac_seaice = l("config_frac_seaice");
  c.config_sfc_albedo = l("config_sfc_albedo");
  c.config_sfc_snowalbedo = l("config_sfc_snowalbedo");
  c.config_sst_update = l("config_sst_update");
  c.config_sstdiurn_update = l("config_sstdiurn_update");
  c.config_deepsoiltemp_update = l("config_deepsoiltemp_update");
  c.config_o3climatology = l("config_o3climatology");
  c.config_microp_re = l("config_microp_re");
  c.config_ysu_pblmix = l("config_ysu_pblmix");
  c.config_urban_physics = l("config_urban_physics");
  c.config_n_microp = i("config_n_microp");
  c.config_microphysics_top = r("config_microphysics_top");
  c.config_radtlw_interval = get("config_radtlw_interval");
  c.config_radtsw_interval = get("config_radtsw_interval");
  c.config_conv_interval = get("config_conv_interval");
  c.config_pbl_interval = get("config_pbl_interval");
  c.config_camrad_abs_update = get("config_camrad_abs_update");
  c.config_greeness_update = get("config_greeness_update");
  c.config_bucket_update = get("config_bucket_update");
  c.config_physics_suite = get("config_physics_suite");
  c.config_microp_scheme = get("config_microp_scheme");
  c.config_convection_scheme = get("config_convection_scheme");
  c.config_lsm_scheme = get("config_lsm_scheme");
  c.config_pbl_scheme = get("config_pbl_scheme");
  c.config_gwdo_scheme = get("config_gwdo_scheme");
  c.config_radt_cld_scheme = get("config_radt_cld_scheme");
  c.config_radt_lw_scheme = get("config_radt_lw_scheme");
  c.config_radt_sw_scheme = get("config_radt_sw_scheme");
  c.config_sfclayer_scheme = get("config_sfclayer_scheme");
  c.config_radt_cld_overlap = get("config_radt_cld_overlap");
  c.config_radt_cld_dcorrlen = get("config_radt_cld_dcorrlen");
  c.config_ngw_scheme = l("config_ngw_scheme");
  c.config_ugwp_diags = l("config_ugwp_diags");
  c.config_oml1d = l("config_oml1d");
  c.config_bucket_radt = r("config_bucket_radt");
  c.config_bucket_rainc = r("config_bucket_rainc");
  c.config_bucket_rainnc = r("config_bucket_rainnc");
  c.config_sounding_interval = get("config_sounding_interval");
  c.config_sfclayer_scheme_namelist = c.config_sfclayer_scheme;
  physics_namelist_check(c, log);
  return c;
}

// mpas_atmphys_control.F:92-332
void physics_namelist_check(PhysicsConfig& c, std::vector<std::string>& log) {
  log.push_back("----- Setting up physics suite '" + c.config_physics_suite + "' -----");
  auto suite = [](std::string& x, const char* v) {
    if (x == "suite") x = v;
  };
  if (c.config_physics_suite == "mesoscale_reference") {  // :157-166
    suite(c.config_microp_scheme, "mp_wsm6");
    suite(c.config_convection_scheme, "cu_ntiedtke");
    suite(c.config_pbl_scheme, "bl_ysu");
    suite(c.config_gwdo_scheme, "bl_ysu_gwdo");
    suite(c.config_radt_lw_scheme, "rrtmg_lw");
    suite(c.config_radt_sw_scheme, "rrtmg_sw");
    suite(c.config_radt_cld_scheme, "cld_fraction");
    suite(c.config_sfclayer_scheme, "sf_monin_obukhov_rev");
    suite(c.config_lsm_scheme, "sf_noah");
  } else if (c.config_physics_suite == "convection_permitting") {  // :168-178
    suite(c.config_microp_scheme, "mp_thompson");
    suite(c.config_convection_scheme, "cu_grell_freitas");
    suite(c.config_pbl_scheme, "bl_mynn");
    suite(c.config_gwdo_scheme, "bl_ysu_gwdo");
    suite(c.config_radt_lw_scheme, "rrtmg_lw");
    suite(c.config_radt_sw_scheme, "rrtmg_sw");
    suite(c.config_radt_cld_scheme, "cld_fraction");
    suite(c.config_sfclayer_scheme, "sf_mynn");
    suite(c.config_lsm_scheme, "sf_noah");
  } else if (c.config_physics_suite == "none") {  // :180-190
    for (std::string* x : {&c.config_microp_scheme, &c.config_convection_scheme, &c.config_pbl_scheme,
                           &c.config_gwdo_scheme, &c.config_radt_lw_scheme, &c.config_radt_sw_scheme,
                           &c.config_radt_cld_scheme, &c.config_sfclayer_scheme, &c.config_lsm_scheme})
      suite(*x, "off");
  } else {
    throw Error("Unrecognized choice of physics suite: config_physics_suite = '" + c.config_physics_suite + "'");
  }
  auto one_of = [](const std::string& x, std::initializer_list<const char*> v) {
    for (const char* s : v)
      if (x == s) return true;
    return false;
  };
  if (!one_of(c.config_microp_scheme, {"off", "mp_kessler", "mp_thompson", "mp_thompson_aerosols", "mp_wsm6"}))
    throw Error("illegal value for config_microp_scheme: " + c.config_microp_scheme);  // :200-209
  if (!one_of(c.config_convection_scheme, {"off", "cu_grell_freitas", "cu_kain_fritsch", "cu_tiedtke", "cu_ntiedtke"}))
    throw Error("illegal value for config_convection_scheme: " + c.config_convection_scheme);  // :212-222
  if (!one_of(c.config_pbl_scheme, {"off", "bl_mynn", "bl_ysu"}))
    throw Error("illegal value for pbl_scheme: " + c.config_pbl_scheme);  // :225-233
  if (!one_of(c.config_gwdo_scheme, {"off", "bl_ysu_gwdo", "bl_ugwp_gwdo"}))
    throw Error("illegal value for gwdo_scheme: " + c.config_gwdo_scheme);  // :236-244
  if (!one_of(c.config_radt_lw_scheme, {"off", "cam_lw", "rrtmg_lw"}))
    throw Error("illegal value for longwave radiation scheme: " + c.config_radt_lw_scheme);  // :247-255
  if (!one_of(c.config_radt_sw_scheme, {"off", "cam_sw", "rrtmg_sw"}))
    throw Error("illegal value for shortwave radiation _scheme: " + c.config_radt_sw_scheme);  // :258-266
  if (!one_of(c.config_radt_cld_scheme, {"off", "cld_incidence", "cld_fraction", "cld_fraction_thompson"}))
    throw Error("illegal value for calculation of cloud fraction: " + c.config_radt_cld_scheme);  // :269-278
  if ((c.config_radt_lw_scheme != "off" && c.config_radt_cld_scheme == "off") ||
      (c.config_radt_sw_scheme != "off" && c.config_radt_cld_scheme == "off")) {  // :279-289
    log.push_back("    config_radt_cld_scheme is not set for radiation calculation");
    log.push_back("    switch calculation of cloud fraction to config_radt_cld_scheme = cld_incidence");
    c.config_radt_cld_scheme = "cld_incidence";
  }
  if (!one_of(c.config_sfclayer_scheme, {"off", "sf_mynn", "sf_monin_obukhov", "sf_monin_obukhov_rev"})) {  // :292-310
    throw Error("illegal value for surface layer scheme: " + c.config_sfclayer_scheme);
  } else {
    if (c.config_pbl_scheme == "bl_mynn") {
      c.config_sfclayer_scheme = "sf_mynn";
    } else if (c.config_pbl_scheme == "bl_ysu") {
      if (c.config_sfclayer_scheme != "sf_monin_obukhov" && c.config_sfclayer_scheme != "sf_monin_obukhov_rev")
        throw Error("wrong choice for surface layer scheme with YSU PBL: " + c.config_sfclayer_scheme);
    }
  }
  if (c.config_lsm_scheme != "off" && c.config_sfclayer_scheme == "off")  // :313-321 ('off ' compares as 'off')
    throw Error("land surface scheme: set config_sfclayer_scheme different than off");
  else if (!one_of(c.config_lsm_scheme, {"off", "sf_noah", "sf_noahmp"}))
    throw Error("illegal value for land surface scheme: " + c.config_lsm_scheme);
  c.moist_physics = true;  // :324-331
  if (c.config_microp_scheme == "off" && c.config_convection_scheme == "off" && c.config_lsm_scheme == "off" &&
      c.config_pbl_scheme == "off" && c.config_radt_lw_scheme == "off" && c.config_radt_sw_scheme == "off" &&
      c.config_sfclayer_scheme == "off")
    c.moist_physics = false;
  log.push_back("    config_microp_scheme       = " + c.config_microp_scheme);
  log.push_back("    config_convection_scheme   = " + c.config_convection_scheme);
  log.push_back("    config_pbl_scheme          = " + c.config_pbl_scheme);
  log.push_back("    config_gwdo_scheme         = " + c.config_gwdo_scheme);
  log.push_back("    config_radt_cld_scheme     = " + c.config_radt_cld_scheme);
  log.push_back("    config_radt_lw_scheme      = " + c.config_radt_lw_scheme);
  log.push_back("    config_radt_sw_scheme      = " + c.config_radt_sw_scheme);
  log.push_back("    config_sfclayer_scheme     = " + c.config_sfclayer_scheme);
  log.push_back("    config_lsm_scheme          = " + c.config_lsm_scheme);
}

ActivePackages physics_packages(const PhysicsConfig& c) {
  ActivePackages p;
  // mpas_atmphys_packages.F:84-176 (after physics_namelist_check)
  if (c.config_microp_scheme == "mp_kessler") p.on.insert("mp_kessler_in");
  if (c.config_microp_scheme == "mp_thompson") p.on.insert("mp_thompson_in");
  if (c.config_microp_scheme == "mp_thompson_aerosols") p.on.insert("mp_thompson_aers_in");
  if (c.config_microp_scheme == "mp_wsm6") p.on.insert("mp_wsm6_in");
  if (c.config_convection_scheme == "cu_grell_freitas") p.on.insert("cu_grell_freitas_in");
  if (c.config_convection_scheme == "cu_kain_fritsch") p.on.insert("cu_kain_fritsch_in");
  if (c.config_convection_scheme == "cu_tiedtke" || c.config_convection_scheme == "cu_ntiedtke")
    p.on.insert("cu_ntiedtke_in");
  if (c.config_pbl_scheme == "bl_mynn") p.on.insert("bl_mynn_in");
  if (c.config_pbl_scheme == "bl_ysu") p.on.insert("bl_ysu_in");
  if (c.config_lsm_scheme == "sf_noahmp") p.on.insert("sf_noahmp_in");
  // Registry.xml:485-490 (active_when, evaluated before physics_namelist_check, on the namelist's value)
  if (c.config_les_model != "none") p.on.insert("les");
  const std::string& sl = c.config_sfclayer_scheme_namelist;
  if ((c.config_physics_suite != "none" && sl == "suite") || sl == "sf_monin_obukhov" ||
      sl == "sf_monin_obukhov_rev" || sl == "sf_mynn")
    p.on.insert("sfclayer");
  // mpas_atm_core_interface.F:145-255: iau, limited_area, jedi_da off; no_invariant_stream (no invariant stream with
  // input_interval initial_only in the run's streams); the UGWP streams with bl_ugwp_gwdo only (refused)
  p.on.insert("no_invariant_stream");
  return p;
}

void check_physics_supported(const PhysicsConfig& c) {
  auto no = [](const std::string& what) { throw Error("physics: " + what + " has no verified port"); };
  if (c.config_physics_suite == "none") return;  // the dycore's own checks refuse any scheme it does not run
  if (c.config_physics_suite != "mesoscale_reference") no("config_physics_suite = '" + c.config_physics_suite + "'");
  auto want = [&](const std::string& v, const char* name, const char* value) {
    if (v != value) no(std::string(name) + " = '" + v + "' (mesoscale_reference: '" + value + "')");
  };
  want(c.config_microp_scheme, "config_microp_scheme", "mp_wsm6");
  want(c.config_convection_scheme, "config_convection_scheme", "cu_ntiedtke");
  want(c.config_pbl_scheme, "config_pbl_scheme", "bl_ysu");
  want(c.config_gwdo_scheme, "config_gwdo_scheme", "bl_ysu_gwdo");
  want(c.config_radt_lw_scheme, "config_radt_lw_scheme", "rrtmg_lw");
  want(c.config_radt_sw_scheme, "config_radt_sw_scheme", "rrtmg_sw");
  want(c.config_radt_cld_scheme, "config_radt_cld_scheme", "cld_fraction");
  want(c.config_sfclayer_scheme, "config_sfclayer_scheme", "sf_monin_obukhov_rev");
  want(c.config_lsm_scheme, "config_lsm_scheme", "sf_noah");
  if (c.config_les_model != "none") no("config_les_model = '" + c.config_les_model + "'");
  // the time-dependent surface updates and the bucket limits: not ported
  if (c.config_sst_update) no("config_sst_update = true (physics_update_sst)");
  if (c.config_sstdiurn_update) no("config_sstdiurn_update = true (physics_update_sstskin)");
  if (c.config_deepsoiltemp_update) no("config_deepsoiltemp_update = true (physics_update_deepsoiltemp)");
  if (c.config_bucket_update != "none") no("config_bucket_update = '" + c.config_bucket_update + "'");
  if (c.config_conv_interval != "none") no("config_conv_interval = '" + c.config_conv_interval + "'");
  if (c.config_pbl_interval != "none") no("config_pbl_interval = '" + c.config_pbl_interval + "'");
  if (c.config_microp_re) no("config_microp_re = true (the WSM6 effective radii)");
  if (c.config_oml1d) no("config_oml1d = true");
  if (c.config_ngw_scheme) no("config_ngw_scheme = true");
  if (c.config_urban_physics) no("config_urban_physics = true");
  if (c.config_n_microp != 1) no("config_n_microp != 1");
  // false (xice_threshold 0.5, MG:740-744; the surface layer without its open-water call, whose values the sea-ice
  // driver blends, mpas_atmphys_driver_seaice.F:275-281, 387-405) is not ported; namelist_table.inc refuses it on
  // the command line too
  if (!c.config_frac_seaice) no("config_frac_seaice = false");
  if (c.config_sounding_interval != "none") no("config_sounding_interval = '" + c.config_sounding_interval + "'");
}

}  // namespace cppmpas::phys
