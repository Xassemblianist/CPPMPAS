// SPDX-License-Identifier: GPL-3.0-or-later
// The restart stream of MPAS-Atmosphere v8.4.2 as atm_core_init reads it (see mpas_restart.hpp).
#include "io/mpas_restart.hpp"

#include <algorithm>
#include <set>

#include "core/error.hpp"
#include "io/netcdf_file.hpp"
#include "io/registry.hpp"

namespace cppmpas {

namespace {

// Fields of the restart stream beyond the input stream's state that the restart run reads instead of computing them:
// atm_init_coupled_diagnostics (mpas_atm_time_integration.F:7440-7665) is skipped for a restart run
// (mpas_atm_core.F:522), so theta_m, rho_zz (state, time level 1) and the coupled diag fields come from the file.
const std::vector<std::string> kRestartCoupled = {"theta_m",       "rho_zz", "pressure_p", "exner",    "exner_base",
                                                  "rtheta_base",   "pressure_base", "ru", "rw", "rtheta_p", "rho_p"};

// Restart-stream fields no ported routine reads on the verified path (physics suite 'none', no LES):
//  * relhum: diagnosed for output only (mpas_atm_diagnostics), read by no routine of the dycore;
//  * circulation: no routine of MPAS v8.4.2 sets or reads it (TI:53, 2232, 2266 remove its last user);
//  * qv_init: mesh pool; atm_compute_dyn_tend_work receives it (TI:6083, 6118, 6224) and never references it;
//  * tend_sfc_pressure: set and read by the physics interface only (mpas_atmphys_interface.F:956-962);
//  * rthratenlw, rthratensw: added by physics_get_tend only with radiation on (mpas_atmphys_todynamics.F:468-484);
// and, by rule, every field of the sfc_input and diag_physics pools: read only by the physics drivers (and ustm / hfx
// / qfx by the LES surface branch of atm_compute_dyn_tend, TI:5969-5971), all of which check_dycore_supported refuses.
const std::set<std::string> kRestartNotUsed = {"relhum",     "circulation", "qv_init", "tend_sfc_pressure",
                                               "rthratenlw", "rthratensw"};

std::string rtrim(std::string s) {
  while (!s.empty() && (s.back() == ' ' || s.back() == '\0')) s.pop_back();
  return s;
}

}  // namespace

bool restart_field_not_used(const std::string& name) {
  if (kRestartNotUsed.count(name)) return true;
  const RegistryVar* rv = find_registry_var(name);
  return rv && (rv->var_struct == "sfc_input" || rv->var_struct == "diag_physics");
}

void check_stream_time(const std::string& path, const std::string& xtime, const std::string& start_stamp) {
  if (rtrim(xtime) != start_stamp)
    throw Error("File " + path + " does not contain the time " + start_stamp + " (its xtime is '" + rtrim(xtime) +
                "'; MPAS reads the initial state at the start time, mpas_stream_manager.F:3776-3786)");
}

RestartInput read_mpas_restart(const std::string& path, const std::vector<std::string>& scalars,
                               const std::string& start_stamp) {
  RestartInput r;
  // the mesh (invariant stream, Registry.xml:657) and the restart-stream fields the input stream has as well:
  // initial_time, xtime, u_init, v_init, t_init, u, w, rho, theta, rho_base, theta_base, surface_pressure, scalars
  r.in = read_mpas_init(path, scalars);
  MpasInput& in = r.in;
  // the one record of the file (the reader refuses more) must be the start time (MPAS_STREAM_EXACT_TIME)
  check_stream_time(path, in.state.xtime, start_stamp);
  const NcFile f = NcFile::open_read(path);
  for (const auto& n : kRestartCoupled) {
    if (!f.has_var(n)) throw Error(path + ": variable " + n + " of the restart stream is missing");
    read_field_ref(f, in, n);
  }
  if (!f.has_var("Time")) throw Error(path + ": variable Time of the restart stream is missing");
  read_field_ref(f, in, "Time");  // state Time, time level 1 (Registry.xml:668)
  const std::pair<const char*, FArray<real>*> extra[] = {
      {"ru_p", &r.ru_p}, {"rw_p", &r.rw_p}, {"rt_diabatic_tend", &r.rt_diabatic_tend}, {"rthdynten", &r.rthdynten}};
  std::vector<std::string> read_extra;
  for (const auto& [n, a] : extra) {
    if (!f.has_var(n)) throw Error(path + ": variable " + std::string(n) + " of the restart stream is missing");
    read_registry_field(f, in.mesh, n, *a);
    read_extra.push_back(n);
  }
  // every other variable of the restart stream in the file: one the port knows no ported routine reads
  for (const auto& n : restart_stream_vars()) {
    if (n == "scalars") continue;  // var_array: the active constituents, read by read_mpas_init
    if (!f.has_var(n)) continue;   // a field of an inactive package: MPAS does not read it either
    const auto& vr = in.source.vars_read;
    if (std::find(vr.begin(), vr.end(), n) != vr.end()) continue;
    if (std::find(read_extra.begin(), read_extra.end(), n) != read_extra.end()) continue;
    if (restart_field_not_used(n)) continue;
    throw Error(path + ": MPAS reads variable " + n + " from the restart stream; CPPMPAS has no port for it");
  }
  return r;
}

void restart_dycore_fields(const RestartInput& r, DycoreFields& f) {
  const auto copy = [](const FArray<real>& from, FArray<real>& to, const char* name) {
    if (from.extents() != to.extents())
      throw Error(std::string("restart_dycore_fields: ") + name + " read with other extents than the dycore holds");
    std::copy(from.data(), from.data() + from.size(), to.data());
  };
  copy(r.ru_p, f.ru_p, "ru_p");
  copy(r.rw_p, f.rw_p, "rw_p");
  copy(r.rt_diabatic_tend, f.rt_diabatic_tend, "rt_diabatic_tend");
  copy(r.rthdynten, f.rthdynten, "rthdynten");
  if (!r.in.state.Time.allocated() || r.in.state.Time.size() != 1) throw Error("restart_dycore_fields: Time was not read");
  f.Time1 = r.in.state.Time(1);
}

}  // namespace cppmpas
