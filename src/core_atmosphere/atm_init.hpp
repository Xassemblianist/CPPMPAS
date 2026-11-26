// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// The computations MPAS-Atmosphere's atm_mpas_init_block (src/core_atmosphere/mpas_atm_core.F:366-617) makes before
// the first time step, in MPAS's order, on one rank: a cold start (config_do_restart = .false.), or a restart
// (config_do_restart = .true., the state read by io/mpas_restart.hpp), which skips atm_init_coupled_diagnostics
// (CORE:522). config_do_DAcycling = .true. is not ported.
// Includes the initial atm_compute_solve_diagnostics (CORE:532-535, core_atmosphere/atm_solve_diagnostics.hpp) and the
// first mpas_reconstruct (CORE:546-563, operators/vector_reconstruct.hpp).
// Not included (later unit): for moist physics (SC) physics_tables_init / physics_registry_init / physics_run_init /
// physics_init (CORE:565-597). These change no field atm_mpas_init_block computes (SC initdump test).
#include <string>
#include <vector>

#include "core/mesh.hpp"
#include "core/namelist.hpp"

namespace cppmpas {

// Allocates the fields atm_mpas_init_block sets with MPAS's extents and Registry defaults (mesh pool:
// Registry.xml:1434-1669; state theta_m / rho_zz; diag ru ... pressure_base). restart: the state and diag ones were
// read from the restart stream and are kept.
void allocate_init_fields(Mesh& mesh, State& state, Diag& diag, bool restart = false);

// mpas_atm_core.F:1184-1271
void atm_compute_signs(Mesh& mesh, std::vector<std::string>& log);
// mpas_atm_core.F:459-473
void atm_compute_inverses(Mesh& mesh);
// mpas_atm_core.F:1620-1850
void atm_initialize_deformation_weights(Mesh& mesh, int nCells, bool on_a_sphere, real sphere_radius);
// mpas_atm_core.F:1391-1536
void atm_adv_coef_compression(Mesh& mesh);
// mpas_atm_core.F:1539-1558
void atm_couple_coef_3rd_order(Mesh& mesh, const InitConfig& cfg);
// mpas_atm_core.F:1124-1181
void atm_compute_mesh_scaling(Mesh& mesh, const InitConfig& cfg);
// mpas_atm_core.F:1274-1388
void atm_compute_damping_coefs(Mesh& mesh, const InitConfig& cfg);
// src/core_atmosphere/dynamics/mpas_atm_boundaries.F:659-753
void mpas_atm_setup_bdy_masks(Mesh& mesh);
// src/core_atmosphere/dynamics/mpas_atm_time_integration.F:7440-7665 (time level 1, whole cell / edge ranges)
void atm_init_coupled_diagnostics(const Mesh& mesh, State& state, Diag& diag);

// atm_mpas_init_block (mpas_atm_core.F:366-617) without the parts listed at the top. log receives the messages MPAS
// writes to its log (e.g. CORE:1263); allocates the fields first.
void atm_mpas_init_block(Mesh& mesh, State& state, Diag& diag, const InitConfig& cfg, std::vector<std::string>& log);

// atm_core_init (mpas_atm_core.F:155-201): nominalMinDc = config_len_disp when that is > 0, else config_len_disp =
// nominalMinDc; both <= 0 is an error. cfg.config_len_disp is updated as MPAS updates the config.
void atm_core_init_len_disp(Mesh& mesh, InitConfig& cfg, std::vector<std::string>& log);
// startTimeStamp of the clock (mpas_atm_core.F:238-239) for a config_start_time in YYYY-MM-DD_hh:mm:ss form.
std::string atm_start_time_stamp(const InitConfig& cfg);
// mpas_atm_run_compatibility (mpas_atm_core.F:1578-1617): limited-area and dynamics checks.
void mpas_atm_run_compatibility(const Mesh& mesh, const InitConfig& cfg, std::vector<std::string>& log);
// What atm_core_init (mpas_atm_core.F:155-262) does to the fields of a cold start: len_disp, the compatibility
// checks, atm_mpas_init_block, then xtime = startTimeStamp (:261-262). cfg.config_len_disp is updated as in
// atm_core_init_len_disp.
void atm_core_init(Mesh& mesh, State& state, Diag& diag, InitConfig& cfg, std::vector<std::string>& log);

// Arrays atm_core_init sets that are proven bit-identical against the MPAS references (tests init_restart_*,
// init_step0_*, initdump_*_*): all of them for an even nVertLevels; for an odd nVertLevels exner, exner_base,
// pressure_p and pressure_base go through the scalar pow of the last level, which no reference exercises. With
// scaled_mesh (mesh_scaling_unverified) the four meshScaling* arrays are not verified either.
std::vector<std::string> init_block_verified_fields(int nVertLevels, bool scaled_mesh = false);
// config_h_ScaleWithMesh with a meshDensity(1:nCells) other than 1: the mesh scaling (mpas_atm_core.F:1124-1181)
// takes pow of values no reference has (JW, the only one with config_h_ScaleWithMesh = true, has meshDensity = 1).
bool mesh_scaling_unverified(const Mesh& mesh, const InitConfig& cfg);
// Every array atm_core_init sets (or overwrites).
const std::vector<std::string>& init_block_computed_fields();

}  // namespace cppmpas
