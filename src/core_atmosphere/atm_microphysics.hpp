// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// The cloud microphysics of MPAS-Atmosphere as the supercell reference runs it (config_physics_suite = 'none',
// config_microp_scheme = 'mp_kessler', every other scheme off): the parts of the physics initialisation of
// atm_mpas_init_block that reach the microphysics (src/core_atmosphere/mpas_atm_core.F:565-597 ->
// physics/mpas_atmphys_init.F physics_init, physics/mpas_atmphys_manager.F physics_run_init), and the call at the
// bottom of atm_srk3 (dynamics/mpas_atm_time_integration.F:2806-2817): pre_microphysics and post_microphysics
// (empty without OpenACC, physics/mpas_atmphys_interface.F:562-595, 1141-1182) around driver_microphysics
// (physics/mpas_atmphys_driver_microphysics.F:317-500): precip_from_MPAS, microphysics_from_MPAS
// (mpas_atmphys_interface.F:598-823), kessler (physics/mp_kessler.hpp, shared with CPPWRF), precip_to_MPAS and
// microphysics_to_MPAS (mpas_atmphys_interface.F:826-1138). One rank, one thread (its..ite = 1..nCellsSolve).
// Every other scheme or branch is refused.
#include <string>

#include "core/mesh.hpp"

namespace cppmpas {

// Registry fields of the physics pools the Kessler path reads or writes, and the module variables of
// mpas_atmphys_vars that physics_run_init sets.
struct MicrophysicsFields {
  bool active = false;  // config_microp_scheme /= 'off' (moist_physics, mpas_atmphys_control.F:310-318)
  int mp_top_level = 0;          // tend_physics (mpas_atmphys_init.F:276-282)
  real dt_dyn = 0, dt_microp = 0;  // mpas_atmphys_manager.F:707-708
  FArray<int> i_rainnc;          // diag_physics (nCells+1)
  FArray<real> rainnc, rainncv, precipw;  // diag_physics (nCells+1)
  FArray<real> dtheta_dt_mp;     // diag_physics (nVertLevels, nCells+1)
  FArray<real> tend_sfc_pressure;  // diag (nCells+1)
};

// Refuses (Error) the physics configurations without a verified port. microp_scheme is the effective
// config_microp_scheme (with suite 'none': 'suite' and 'off' are off).
void check_microphysics_supported(const std::string& namelist_path, const Mesh& m);

// The Kessler part of mpas_atm_core.F:565-597: physics_run_init's dt_dyn / dt_microp (mpas_atmphys_manager.F:
// 707-708) and physics_init's i_rainnc = 0 and mp_top_level (mpas_atmphys_init.F:268-282). dt is the clock's
// time step.
void microphysics_init(const std::string& namelist_path, const Mesh& m, real dt, MicrophysicsFields& mp);

// driver_microphysics on time level 2 (theta_m2, rho_zz2, w2, scalars2 of the state pool) and the diag / tend
// fields it updates: exner, pressure_p, rtheta_p, rt_diabatic_tend, State::surface_pressure.
void driver_microphysics(const Mesh& m, State& s, FArray<real>& theta_m, FArray<real>& rho_zz, FArray<real>& scalars,
                         Diag& d, FArray<real>& rt_diabatic_tend, MicrophysicsFields& mp);

}  // namespace cppmpas
