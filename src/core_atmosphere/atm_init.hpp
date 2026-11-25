// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// The computations MPAS-Atmosphere's atm_mpas_init_block (src/core_atmosphere/mpas_atm_core.F:366-617) makes before
// the first time step, in MPAS's order, for a cold start (config_do_restart = .false.) on one rank.
// Not included (later units): atm_compute_solve_diagnostics and the first mpas_reconstruct (CORE:524-541, 545-563),
// and for moist physics (SC) physics_tables_init / physics_registry_init / physics_run_init / physics_init
// (CORE:568-596). These do not change any field computed here.
#include <string>
#include <vector>

#include "core/mesh.hpp"
#include "core/namelist.hpp"

namespace cppmpas {

// Allocates the fields atm_mpas_init_block sets with MPAS's extents and Registry defaults (mesh pool:
// Registry.xml:1434-1669; state theta_m / rho_zz; diag ru ... pressure_base).
void allocate_init_fields(Mesh& mesh, State& state, Diag& diag);

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

}  // namespace cppmpas
