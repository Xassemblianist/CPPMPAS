// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// MPAS-Atmosphere's physics in the model run, statement by statement:
//  * mpas_atm_core.F:565-597 (atm_core_init: physics_tables_init, physics_registry_init, physics_run_init,
//    physics_init) and :1031-1037 (atm_do_timestep: physics_timetracker, physics_driver);
//  * physics/mpas_atmphys_control.F:335-390 (physics_registry_init);
//  * physics/mpas_atmphys_manager.F (MG): physics_run_init (MG:378-766) with the clock alarms of the radiation, the
//    greenness update and the physics diagnostics, and physics_timetracker (MG:137-375);
//  * physics/mpas_atmphys_init.F (IN): physics_init (IN:83-481), init_dirs_forphys (IN:556-580), r3_normalize
//    (IN:583-596), init_physics_flags (IN:484-545);
//  * physics/mpas_atmphys_driver.F (DR): physics_driver (DR:91-399);
//  * the coupling of mpas_atm_time_integration.F (TI): physics_get_tend at TI:2141-2151, the dynamics tendencies of
//    the convection (TI:2758-2788) and the microphysics at TI:2806-2817 (driver_microphysics, DRV:317-500).
// One rank, one block, one thread: its..ite = 1..nCellsSolve, the halo exchanges are no-ops.
//
// Ported for config_physics_suite = 'mesoscale_reference' with RG's options (physics/mpas_atmphys_control.hpp
// check_physics_supported): WSM6, New Tiedtke, YSU with GWDO, the revised MM5 surface layer, Noah with the sea-ice
// driver and the surface diagnostics, RRTMG LW and SW with 'cld_fraction' and the ozone climatology. A diagnostics
// alarm that rings (l_diags: an output or diagnostics stream with an interval) is refused: the WSM6 diagnostics
// (reflectivity, relhum) are not ported.
#include <functional>
#include <string>
#include <vector>

#include "core/mesh.hpp"
#include "io/registry_pools.hpp"
#include "physics/mpas_atmphys_control.hpp"
#include "physics/mpas_atmphys_interface.hpp"
#include "physics/mpas_atmphys_radiation.hpp"
#include "physics/mpas_atmphys_sfclayer.hpp"
#include "physics/mp_wsm6.hpp"
#include "physics/sf_noah.hpp"

namespace cppmpas::phys {

// A recurring alarm of the MPAS clock (mpas_timekeeping.F mpas_add_clock_alarm / mpas_is_alarm_ringing /
// mpas_reset_clock_alarm), in whole seconds since the start time: it rings once the clock reaches its ring time and a
// reset moves the ring time past the clock by whole intervals.
struct ClockAlarm {
  bool defined = false;
  long long ring = 0, interval = 0;
  bool ringing(long long now) const { return defined && now >= ring; }
  void reset(long long now) {
    while (ring <= now) ring += interval;
  }
};

struct PhysicsRun {
  PhysicsConfig cfg;
  std::string table_dir;  // LANDUSE.TBL, VEGPARM.TBL, SOILPARM.TBL, GENPARM.TBL, OZONE_*.TBL
  std::string rrtmg_lw_tables, rrtmg_sw_tables;  // the RRTMG table images (physics/mpas_atmphys_radiation.hpp)
  // tables of the schemes
  sfclayrev::SFClayRevTables<real> sfclayrev_tables;
  noah::NoahTables noah_tables;
  RrtmgTables rrtmg;
  wsm6::Wsm6Params<real> wsm6;
  // mpas_atmphys_vars / mpas_atmphys_manager
  real dt_dyn = 0, dt_microp = 0, dt_cu = 0, dt_pbl = 0, xice_threshold = 0;
  int n_cu = 1, n_microp = 1, mp_top_level = 0, isurban = 0;
  int has_reqc = 0, has_reqi = 0, has_reqs = 0;
  bool l_radtlw = false, l_radtsw = false, l_conv = false, l_diags = false, l_acrain = false, l_acradt = false;
  bool f_qc = false, f_qr = false, f_qi = false, f_qs = false, f_qg = false;
  int index_qv = 0, index_qc = 0, index_qr = 0, index_qi = 0, index_qs = 0, index_qg = 0;
  RadiationVars rad;
  ClockAlarm radtlw_alarm, radtsw_alarm, green_alarm, diag_alarm;
  std::string start_time;  // YYYY-MM-DD_hh:mm:ss of the clock's start
};

// The arrays of the dycore the physics reads or writes besides the pools (core_atmosphere/atm_dycore.hpp
// DycoreFields): they are MPAS's tend_physics / tend / diag fields of the same names.
struct DycorePhysicsArrays {
  FArray<real>* rthdynten = nullptr;  // tend_physics (written by atm_compute_dyn_tend, TI:6912)
  FArray<real>* tend_ru_physics = nullptr;
  FArray<real>* tend_rtheta_physics = nullptr;
  FArray<real>* tend_rho_physics = nullptr;
  FArray<real>* tend_scalars = nullptr;  // tend 'scalars_tend'
  FArray<real>* rt_diabatic_tend = nullptr;
};

// atm_core_init's physics (mpas_atm_core.F:565-597) for a cold start: mminlu, physics_registry_init (dzs),
// physics_run_init and physics_init. clock_now: the clock's current (= start) time stamp. log: MPAS's messages.
void physics_core_init(PhysicsRun& r, const Mesh& m, State& s, Diag& d, RegistryPools& p,
                       const std::string& clock_now, std::vector<std::string>& log);

// physics_timetracker (MG:137-375) at the clock time now_stamp, now_s seconds after the start; dt the time step.
void physics_timetracker(PhysicsRun& r, const Mesh& m, RegistryPools& p, const std::string& now_stamp,
                         long long now_s, real dt, std::vector<std::string>& log);

// physics_driver (DR:91-399) on time level 1 of the state.
void physics_driver(PhysicsRun& r, const Mesh& m, State& s, Diag& d, RegistryPools& p, int itimestep,
                    real xtime_s, const DycorePhysicsArrays& a);

// TI:2141-2151: physics_get_tend (TD:98-279) with rk_step = 1, dynamics_substep = 1; rho_zz2: rho_zz of time level 2
void physics_get_tend_for_dynamics(PhysicsRun& r, const Mesh& m, const State& s, Diag& d, RegistryPools& p,
                                   const FArray<real>& rho_zz2, const DycorePhysicsArrays& a);

// TI:2758-2788: rqvdynten and rthdynten for cu_ntiedtke from the scalars of both time levels and theta_m of level 2
void dynamics_tendencies_for_convection(PhysicsRun& r, const Mesh& m, const State& s, RegistryPools& p,
                                        const FArray<real>& scalars2, const FArray<real>& theta_m2, real dt,
                                        bool config_monotonic, const DycorePhysicsArrays& a);

// TI:2806-2817: driver_microphysics (DRV:317-500, mp_wsm6) on time level 2
void microphysics_for_dynamics(PhysicsRun& r, const Mesh& m, State& s, Diag& d, RegistryPools& p,
                               FArray<real>& theta_m2, FArray<real>& rho_zz2, FArray<real>& scalars2,
                               const FArray<real>& w2, int itimestep, const DycorePhysicsArrays& a);

// CT:459-533 (mpas_atm_core.F:1597-1603): YSU GWDO needs the subgrid orography fields
void physics_compatibility_check(const PhysicsConfig& c, const Mesh& m, RegistryPools& p);

}  // namespace cppmpas::phys
