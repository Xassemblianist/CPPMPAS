// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// MPAS-Atmosphere's diagnostics framework (src/core_atmosphere/diagnostics/), statement by statement:
//  * the manager, mpas_atm_diagnostics_manager.F (MG): mpas_atm_diag_setup (MG:31-62), mpas_atm_diag_update
//    (MG:75-87), mpas_atm_diag_compute (MG:100-121), mpas_atm_diag_reset (MG:134-146);
//  * mpas_atm_diagnostics_utils.F: MPAS_field_will_be_written (UT:79-121) and mpas_stream_inclusion_count
//    (UT:130-213), answered by the caller's stream manager (DiagnosticsContext);
//  * the diagnostics the manager calls: diagnostic_template (no-ops), mpas_isobaric_diagnostics.F,
//    mpas_cloud_diagnostics.F, mpas_convective_diagnostics.F, mpas_pv_diagnostics.F, mpas_soundings.F.
// The diagnostics read and write MPAS's state / diag pools: the arrays CPPMPAS holds in State and Diag (core/mesh.hpp)
// and every other diag field in the Registry pools (io/registry_pools.hpp; mslp, cape, relhum, ...).
// Ported: isobaric and convective diagnostics. Cloud and PV diagnostics and soundings have no verified port: a run in
// which one of their fields would be written (or config_sounding_interval /= 'none') is an error.
#include <functional>
#include <string>

#include "core/mesh.hpp"
#include "io/registry_pools.hpp"

namespace cppmpas::diag {

struct DiagnosticsContext {
  const Mesh* mesh = nullptr;
  State* state = nullptr;  // time level 1 (the state after the time-level shift, as atm_core_run computes them)
  Diag* diag = nullptr;
  RegistryPools* pools = nullptr;
  // UT:79-121 MPAS_field_will_be_written: the field is active in an active output stream whose output alarm rings
  std::function<bool(const std::string&)> field_will_be_written;
  // UT:130-213 mpas_stream_inclusion_count(name, direction=MPAS_STREAM_OUTPUT): the active output streams with an
  // output interval other than 'none' that hold the field
  std::function<int(const std::string&)> output_stream_inclusion_count;
  // module variables of mpas_convective_diagnostics.F (CV:41-43), set by convective_diagnostics_setup
  bool is_needed_updraft_helicity = false, is_needed_w_max = false, is_needed_lml_wsp_max = false;
};

// MG:31-62 (diagnostic_template_setup ... soundings_setup); config_sounding_interval from the namelist
void mpas_atm_diag_setup(DiagnosticsContext& c, const std::string& config_sounding_interval);
void mpas_atm_diag_update(DiagnosticsContext& c);   // MG:75-87
void mpas_atm_diag_compute(DiagnosticsContext& c);  // MG:100-121
void mpas_atm_diag_reset(DiagnosticsContext& c);    // MG:134-146

// mpas_isobaric_diagnostics.F isobaric_diagnostics_compute (IS:62-...)
void isobaric_diagnostics_compute(DiagnosticsContext& c);

// mpas_convective_diagnostics.F (CV)
void convective_diagnostics_setup(DiagnosticsContext& c);    // CV:53-85
void convective_diagnostics_update(DiagnosticsContext& c);   // CV:99-198
void convective_diagnostics_compute(DiagnosticsContext& c);  // CV:227-487
void convective_diagnostics_reset(DiagnosticsContext& c);    // CV:501-525

}  // namespace cppmpas::diag
