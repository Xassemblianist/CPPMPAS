// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// The initial state of a restart run (config_do_restart = .true.): what MPAS-Atmosphere's atm_core_init reads from
// the "restart" stream instead of the "input" stream (src/core_atmosphere/mpas_atm_core.F:129-149; the restart stream
// is Registry.xml:650-925), on one rank. The mesh comes from the same file: atm_get_mesh_stream returns 'restart'
// for a restart run (mpas_atm_core_interface.F:394-402).
//
// The stream is read at the start time of the clock with MPAS_STREAM_EXACT_TIME (mpas_stream_manager.F:3529-3531):
// a file without a record whose xtime is the start time is an error (:3783-3786). Every variable of the restart
// stream the file holds is either read into the array the port uses or named in a list of fields no ported routine
// reads (with the reason); any other one is an error.
#include <string>
#include <vector>

#include "core/mesh.hpp"
#include "core_atmosphere/atm_dycore.hpp"
#include "io/mpas_reader.hpp"

namespace cppmpas {

struct RestartInput {
  // mesh, the state of time level 1 (u, w, rho_zz, theta_m, scalars, xtime, Time, initial_time, ...), and the diag
  // fields of the restart stream that atm_init_coupled_diagnostics would set on a cold start (ru, rw, rho_p,
  // rtheta_base, rtheta_p, exner, exner_base, pressure_p, pressure_base)
  MpasInput in;
  // restart-stream fields the port holds in DycoreFields (diag ru_p / rw_p, tend rt_diabatic_tend, tend_physics
  // rthdynten), copied there by restart_dycore_fields once the dycore fields are allocated
  FArray<real> ru_p, rw_p, rt_diabatic_tend, rthdynten;
};

// Reads a restart file as MPAS's restart stream does; start_stamp is the clock's start time (atm_start_time_stamp).
// scalars: the active constituents (as read_mpas_init).
RestartInput read_mpas_restart(const std::string& path, const std::vector<std::string>& scalars,
                               const std::string& start_stamp);

// The variables of the restart stream that MPAS reads but no ported routine uses on the verified path, each with the
// reason (for the tests and the error message); physics pools are covered by a rule, see mpas_restart.cpp.
bool restart_field_not_used(const std::string& name);

// The restart-stream values that live in DycoreFields (after allocate_dycore_fields): ru_p, rw_p, rt_diabatic_tend,
// rthdynten, and the state Time of time level 1 (f.Time1 = Time read, Registry.xml:668).
void restart_dycore_fields(const RestartInput& r, DycoreFields& f);

// MPAS reads the initial-conditions stream at the clock's start time with MPAS_STREAM_EXACT_TIME
// (mpas_stream_manager.F:3776-3786): the file's xtime must be the start time, or the read fails.
void check_stream_time(const std::string& path, const std::string& xtime, const std::string& start_stamp);

}  // namespace cppmpas
