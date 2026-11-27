// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// One time step of MPAS-Atmosphere's dynamical core as the mountain-wave (MW) and JW baroclinic-wave references run it: atm_timestep ->
// atm_srk3 (src/core_atmosphere/dynamics/mpas_atm_time_integration.F:1789-1850, 1855-2882; tag TI below), with the
// dissipation of mpas_atm_dissipation_models.F (DM, core_atmosphere/atm_dissipation.hpp), the time-level shift and
// the output diagnostics of atm_core_run (mpas_atm_core.F:816-848, 916-965; CORE). One rank, one block, one thread
// (the reference is built without OpenMP: every thread range is the whole range, every halo exchange a no-op,
// DYCORE_MAP.md section 3).
//
// Ported and verified (MW: mpas842-dbl-mw-001, planar, steps 1-10 and 300 / 301 / 600; JW: mpas842-dbl-jw-001,
// x1.40962 sphere, steps 1-10): SRK3 of order 2, coupled dynamics-transport (config_split_dynamics_transport =
// false), no scalar advection, horizontal mixing 2d_fixed (MW) or 2d_smagorinsky with the 4th-order filters it
// switches on (JW), vertical mixing of the full state (MW) or none (JW, v_*_eddy_visc2 = 0), no LES, no Rayleigh
// damping of u, no CAM top filter, physics suite 'none', no LBCs, no IAU, fEdge = 0 (MW) or u_init = v_init = 0 (JW),
// an even nVertLevels >= 4. Supercell (SC, mpas842-dbl-sc-001 steps 1-10, 600, 601) adds: coupled monotonic scalar
// transport of qv, qc, qr (atm_advance_scalars.hpp), vertical mixing of the perturbation from the initial state
// (config_mix_full = .false., atm_dissipation_perturbation.cpp) and Kessler microphysics at the end of the step
// (atm_microphysics.hpp). Every other configuration is refused by check_dycore_supported.
#include <functional>
#include <string>
#include <vector>

#include "core/mesh.hpp"
#include "core/namelist.hpp"
#include "core_atmosphere/atm_microphysics.hpp"
#include "core_atmosphere/atm_solve_diagnostics.hpp"

namespace cppmpas {

// The namelist options the dycore reads (Registry.xml line of each option in atm_dycore.cpp), from
// effective_namelist (every value already checked against the verified references) and the InitConfig that
// atm_core_init updated (config_len_disp).
struct DycoreConfig {
  std::string config_time_integration;      // Registry.xml:64
  int config_time_integration_order = 2;    // :69
  real config_dt = 0;                       // :74
  bool config_split_dynamics_transport = true;  // :99
  int config_number_of_sub_steps = 4;       // :104
  int config_dynamics_split_steps = 3;      // :109
  real config_h_mom_eddy_visc2 = 0, config_h_mom_eddy_visc4 = 0, config_v_mom_eddy_visc2 = 0;      // :114-124
  real config_h_theta_eddy_visc2 = 0, config_h_theta_eddy_visc4 = 0, config_v_theta_eddy_visc2 = 0;  // :129-139
  std::string config_horiz_mixing;          // :144
  std::string config_les_model;             // :149
  real config_len_disp = 0;                 // :174 (after atm_core_init, mpas_atm_core.F:155-201)
  real config_visc4_2dsmag = 0.05;          // :179
  bool config_mix_scalars = false;          // :184
  real config_del4u_div_factor = 10;        // :189
  bool config_scalar_advection = true;      // :229
  bool config_positive_definite = false;    // :234
  bool config_monotonic = true;             // :239
  real config_coef_3rd_order = 0.25;        // :244
  real config_smagorinsky_coef = 0.125;     // :249
  bool config_mix_full = true;              // :254
  bool config_h_ScaleWithMesh = true;       // :274 (from the InitConfig)
  real config_smdiv = 0.1;                  // :264
  real config_apvm_upwinding = 0.5;         // :269
  real config_mpas_cam_coef = 0;            // :301
  bool config_rayleigh_damp_u = false;      // :311
  bool config_apply_lbcs = false;           // :349
  std::string config_IAU_option;            // :424
  std::string config_physics_suite, config_microp_scheme, config_convection_scheme, config_pbl_scheme;
  std::string config_radt_lw_scheme, config_radt_sw_scheme, config_radt_cld_scheme;
  std::string config_lsm_scheme, config_sfclayer_scheme, config_gwdo_scheme;
  bool config_oml1d = false;
  std::string config_calendar_type;
};
DycoreConfig read_dycore_config(const std::string& namelist_path, const InitConfig& init_cfg);

// Everything the dycore holds beyond the mesh pool and the fields of atm_core_init: time level 2 of the state pool,
// the diag / tend / tend_physics fields of the dycore and the module work arrays of atm_time_integration, with
// MPAS's extents (one garbage element) and the Registry default 0.
struct DycoreFields {
  // state, time level 2 (Registry.xml:1672-1711); time level 1 is State::u, w, theta_m, rho_zz, scalars
  FArray<real> u2, w2, theta_m2, rho_zz2, scalars2;
  std::string xtime2;
  real Time2 = 0.0;  // state Time (CF time, s) of time level 2; time level 1 is State::Time(1)
  // diag
  FArray<real> ru_save, rw_save, rtheta_p_save, rho_p_save, rho_zz_old_split, cqw, cqu;
  FArray<real> cofrz, cofwr, cofwz, coftz, cofwt, a_tri, alpha_tri, gamma_tri;
  FArray<real> ru_p, rw_p, rtheta_pp, rtheta_pp_old, rho_pp, ruAvg, wwAvg, ruAvg_split, wwAvg_split;
  FArray<real> h_divergence, eddy_visc_horz;
  // tend
  FArray<real> tend_u, tend_w, tend_theta, tend_rho, tend_u_euler, tend_w_euler, tend_theta_euler;
  FArray<real> rt_diabatic_tend, tend_scalars;
  // tend_physics
  FArray<real> rthdynten, tend_ru_physics, tend_rtheta_physics, tend_rho_physics;
  // module arrays of atm_time_integration (TI:93-99): qtot, dpdz and the delsq_* scratch of atm_srk3
  FArray<real> qtot, dpdz, delsq_theta, delsq_w, delsq_divergence, delsq_u, delsq_vorticity;
  SolveDiagWork sdw;
  int moist_start = 1, moist_end = 1, index_qv = 1;
  MicrophysicsFields mp;  // the Kessler path (atm_microphysics.hpp); inactive unless microphysics_init set it
  // The physics suite's coupling (physics/mpas_atmphys_driver.hpp), set by the caller with a physics suite:
  // TI:2141-2151 physics_get_tend, TI:2758-2788 the dynamics tendencies of the convection, TI:2806-2817
  // driver_microphysics. Unset: the suite 'none' (zero physics tendencies, the Kessler path of mp above).
  struct PhysicsHooks {
    std::function<void(const State&, Diag&, DycoreFields&)> get_tend;
    std::function<void(const State&, DycoreFields&, real)> dyn_tend_for_convection;
    std::function<void(State&, Diag&, DycoreFields&)> microphysics;
  } phys;
};
void allocate_dycore_fields(const Mesh& mesh, const State& state, DycoreFields& f);

// Refuses (Error naming the option or the property) every configuration whose code path is not ported and verified.
void check_dycore_supported(const DycoreConfig& cfg, const Mesh& mesh, const State& state);
// Refuses mesh, state, diag and dycore arrays without the extents the routines below index (every exported routine
// calls it: no routine reads or writes outside an array), an nVertLevels the dycore is not verified for, and a mesh
// with boundary or specified-zone cells (bdyMask* / specZoneMask* != 0: the LBC branches are not ported).
void check_dycore_arrays(const Mesh& m, const State& s, const Diag& d, const DycoreFields& f);

// The routines of atm_srk3 in its order (TI line ranges in atm_dycore.cpp). time level 1 = state, 2 = f.*2.
void atm_rk_integration_setup(const Mesh& m, State& s, Diag& d, DycoreFields& f);
void atm_compute_moist_coefficients(const Mesh& m, DycoreFields& f);
void physics_get_tend(const Mesh& m, const State& s, DycoreFields& f);
void atm_compute_vert_imp_coefs(const Mesh& m, const State& s, const Diag& d, DycoreFields& f, real dts);
void atm_compute_dyn_tend(const DycoreConfig& cfg, const Mesh& m, const State& s, Diag& d, DycoreFields& f,
                          int rk_step, real dt);
void atm_set_smlstep_pert_variables(const Mesh& m, DycoreFields& f);
void atm_advance_acoustic_step(const Mesh& m, const State& s, const Diag& d, DycoreFields& f, real dts,
                               int small_step);
void atm_divergence_damping_3d(const DycoreConfig& cfg, const Mesh& m, const State& s, DycoreFields& f, real dts);
void atm_recover_large_step_variables(const Mesh& m, State& s, Diag& d, DycoreFields& f, real dt, int ns,
                                      int rk_step);
void atm_rk_dynamics_substep_finish(const Mesh& m, State& s, Diag& d, DycoreFields& f, int dynamics_substep,
                                    int dynamics_split);

// where (scalars_2 < 0.0) scalars_2 = 0.0 (TI:2799-2800), the whole array including the garbage column
void clamp_negative_scalars(FArray<real>& scalars);

// atm_srk3 (TI:1855-2882) for one step of length dt (the clock's time step).
void atm_srk3(const DycoreConfig& cfg, const Mesh& m, State& s, Diag& d, DycoreFields& f, real dt);
// atm_timestep (TI:1789-1850): atm_srk3, then xtime(2) = now + dt and Time(2) = Time(1) + dt.
void atm_timestep(const DycoreConfig& cfg, const Mesh& m, State& s, Diag& d, DycoreFields& f, real dt,
                  const std::string& xtime_new);
// mpas_pool_shift_time_levels(state) (mpas_atm_core.F:823): time level 2 becomes time level 1 and back.
void shift_time_levels(State& s, DycoreFields& f);
// atm_compute_output_diagnostics(state, 1, diag, mesh) (mpas_atm_core.F:916-965): theta, rho of the diag pool
// (State::theta, State::rho) from time level 1, and pressure = pressure_base + pressure_p (Diag::pressure, allocated
// here on the first call). MPAS calls it only when an output stream writes (mpas_atm_core.F:686, 848).
void atm_compute_output_diagnostics(const Mesh& m, State& s, Diag& d);

// "YYYY-MM-DD_hh:mm:ss" + seconds (proleptic Gregorian calendar, year >= 0; mpas_timekeeping's gregorian clock).
std::string add_seconds(const std::string& timestamp, long long seconds);

}  // namespace cppmpas
