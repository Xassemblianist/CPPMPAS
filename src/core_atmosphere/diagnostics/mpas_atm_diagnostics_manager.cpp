// SPDX-License-Identifier: GPL-3.0-or-later
// mpas_atm_diagnostics_manager.F (MG) and the diagnostics without a port (cloud, PV, soundings): see
// mpas_atm_diagnostics.hpp.
#include <initializer_list>

#include "core/error.hpp"
#include "core_atmosphere/diagnostics/mpas_atm_diagnostics.hpp"

namespace cppmpas::diag {

namespace {

// The diag fields of Registry_cloud.xml and Registry_pv.xml: a diagnostic without a port must not be written.
void refuse_unported(DiagnosticsContext& c, const char* what, std::initializer_list<const char*> fields) {
  for (const char* f : fields)
    if (c.field_will_be_written(f))
      throw Error(std::string(what) + ": the field " + f + " would be written; this diagnostic has no verified port");
}

}  // namespace

void mpas_atm_diag_setup(DiagnosticsContext& c, const std::string& config_sounding_interval) {
  // MG:51-61: diagnostic_template_setup and isobaric / cloud / pv_diagnostics_setup only keep pool pointers
  convective_diagnostics_setup(c);
  // soundings_setup (mpas_soundings.F:86-...): soundings only for config_sounding_interval /= 'none'
  if (config_sounding_interval != "none")
    throw Error("config_sounding_interval = '" + config_sounding_interval + "': the soundings have no verified port");
}

void mpas_atm_diag_update(DiagnosticsContext& c) {
  // MG:84-85: diagnostic_template_update is empty
  convective_diagnostics_update(c);
}

void mpas_atm_diag_compute(DiagnosticsContext& c) {
  // MG:114-119 in this order; diagnostic_template_compute is empty, soundings_compute returns at once when the
  // soundings are off (refused in setup otherwise)
  isobaric_diagnostics_compute(c);
  refuse_unported(c, "cloud_diagnostics_compute (mpas_cloud_diagnostics.F)",
                  {"cldfrac_low_UPP", "cldfrac_mid_UPP", "cldfrac_high_UPP", "cldfrac_tot_UPP"});
  convective_diagnostics_compute(c);
  refuse_unported(c, "pv_diagnostics_compute (mpas_pv_diagnostics.F)",
                  {"ertel_pv", "u_pv", "v_pv", "theta_pv", "vort_pv", "iLev_DT", "depv_dt_lw", "depv_dt_sw",
                   "depv_dt_bl", "depv_dt_cu", "depv_dt_mp", "dtheta_dt_mix", "depv_dt_mix", "depv_dt_diab",
                   "depv_dt_fric", "tend_u_phys", "depv_dt_diab_pv", "depv_dt_fric_pv", "inTropo"});
}

void mpas_atm_diag_reset(DiagnosticsContext& c) {
  // MG:143-144: diagnostic_template_reset is empty
  convective_diagnostics_reset(c);
}

}  // namespace cppmpas::diag
