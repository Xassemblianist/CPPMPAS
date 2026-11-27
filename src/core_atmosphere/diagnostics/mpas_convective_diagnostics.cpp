// SPDX-License-Identifier: GPL-3.0-or-later
// mpas_convective_diagnostics.F (CV): see mpas_atm_diagnostics.hpp.
#include <string>

#include "core/error.hpp"
#include "core_atmosphere/diagnostics/mpas_atm_diagnostics.hpp"

namespace cppmpas::diag {

void convective_diagnostics_setup(DiagnosticsContext& c) {
  // CV:73-82
  c.is_needed_updraft_helicity = c.output_stream_inclusion_count("updraft_helicity_max") > 0;
  c.is_needed_w_max = c.output_stream_inclusion_count("w_velocity_max") > 0;
  c.is_needed_lml_wsp_max = c.output_stream_inclusion_count("wind_speed_level1_max") > 0;
}

void convective_diagnostics_update(DiagnosticsContext& c) {
  if (c.is_needed_updraft_helicity || c.is_needed_w_max || c.is_needed_lml_wsp_max)
    throw Error("convective_diagnostics_update: updraft_helicity_max / w_velocity_max / wind_speed_level1_max have "
                "no verified port");
}

void convective_diagnostics_compute(DiagnosticsContext& c) {
  static const char* const fields[] = {"cape", "cin", "lcl", "lfc", "srh_0_1km", "srh_0_3km", "uzonal_surface",
                                       "uzonal_1km", "uzonal_6km", "umeridional_surface", "umeridional_1km",
                                       "umeridional_6km", "temperature_surface", "dewpoint_surface"};
  for (const char* f : fields)
    if (c.field_will_be_written(f))
      throw Error(std::string("convective diagnostics: ") + f + " has no verified port");
}

void convective_diagnostics_reset(DiagnosticsContext& c) {
  // CV:515-523: the whole arrays (the garbage element included)
  for (const char* f : {"updraft_helicity_max", "w_velocity_max", "wind_speed_level1_max"})
    if (c.field_will_be_written(f)) {
      auto& a = c.pools->R(f);
      for (std::size_t n = 0; n < a.size(); ++n) a.data()[n] = 0.;
    }
}

}  // namespace cppmpas::diag
