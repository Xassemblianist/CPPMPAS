// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// Scalar transport of MPAS-Atmosphere's atm_srk3 (src/core_atmosphere/dynamics/mpas_atm_time_integration.F, TI):
// advance_scalars (TI:2911-3085) with atm_advance_scalars(_work) for rk_steps 1 and 2 (TI:4541-4939) and
// atm_advance_scalars_mono(_work) for rk_step 3 (TI:4953-5785), on the path of the supercell reference:
// coupled dynamics-transport (config_split_dynamics_transport = .false., so advance_density = .false.),
// config_monotonic = .true., config_positive_definite = .false., no LBCs (bdyMaskCell = bdyMaskEdge = 0), every
// edge with nAdvCellsForEdge = 10 and every cell with nEdgesOnCell = 6 (the hexagon cases of TI:4770, 5372, 5425).
// One rank, one thread: the thread ranges are the whole ranges and the halo exchanges (TI:5223, 5630) are no-ops.
// Any other branch is an Error.
#include "core/mesh.hpp"
#include "core_atmosphere/atm_dycore.hpp"

namespace cppmpas {

// Refuses (Error) the scalar-transport configurations and meshes whose path is not ported.
void check_advance_scalars_supported(const DycoreConfig& cfg, const Mesh& m);

// advance_scalars('scalars', ...) at TI:2504-2507 for rk_step 1..3: scalars of time level 1 (State::scalars) and 2
// (f.scalars2), rho_zz of both levels, diag ruAvg / wwAvg, tend scalars_tend (f.tend_scalars).
// rk_timestep[1..3] as atm_srk3 holds it.
void advance_scalars(const DycoreConfig& cfg, const Mesh& m, State& s, DycoreFields& f, int rk_step,
                     const real* rk_timestep);

}  // namespace cppmpas
