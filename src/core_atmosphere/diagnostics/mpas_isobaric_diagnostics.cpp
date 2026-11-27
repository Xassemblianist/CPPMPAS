// SPDX-License-Identifier: GPL-3.0-or-later
// mpas_isobaric_diagnostics.F (IS): see mpas_atm_diagnostics.hpp.
#include <string>

#include "core/error.hpp"
#include "core_atmosphere/diagnostics/mpas_atm_diagnostics.hpp"

namespace cppmpas::diag {

void isobaric_diagnostics_compute(DiagnosticsContext& c) {
  static const char* const fields[] = {"mslp", "t_isobaric", "z_isobaric", "meanT_500_300"};
  static const char* const lev[] = {"relhum", "dewpoint", "temperature", "height", "uzonal", "umeridional", "w",
                                    "vorticity"};
  static const char* const hpa[] = {"50", "100", "200", "250", "500", "700", "850", "925"};
  for (const char* f : fields)
    if (c.field_will_be_written(f)) throw Error(std::string("isobaric diagnostics: ") + f + " has no verified port");
  for (const char* l : lev)
    for (const char* p : hpa)
      if (c.field_will_be_written(std::string(l) + "_" + p + "hPa"))
        throw Error(std::string("isobaric diagnostics: ") + l + "_" + p + "hPa has no verified port");
}

}  // namespace cppmpas::diag
