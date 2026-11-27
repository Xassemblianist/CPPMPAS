// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// MPAS-Atmosphere's microphysics driver for WSM6 (physics/mpas_atmphys_driver_microphysics.F, DRV) and its
// interface to the dynamics (physics/mpas_atmphys_interface.F, PI), statement by statement, on the arrays of the
// MPAS pools (Fortran layout: (k, iCell), 1-based, k fastest; scalars (index, k, iCell)). One rank, one thread:
// its..ite = 1..nCellsSolve, kts..kte = 1..nVertLevels, jts = jte = 1 (mpas_atmphys_manager.F:653-663;
// dt_dyn, dt_microp: :684-686).
//
// driver_microphysics_wsm6 = DRV:317-500 with config_microp_scheme = 'mp_wsm6': allocate_microphysics
// (DRV:98-176), precip_from_MPAS (DRV:503-564), microphysics_from_MPAS (PI:598-823), wsm6 (DRV:437-461, then
// physics/module_mp_wsm6.hpp), precip_to_MPAS (DRV:567-653), microphysics_to_MPAS (PI:826-1138).
// With l_diags (the step before a diagnostics alarm, physics_timetracker MG:365-372), after wsm6: 
// compute_radar_reflectivity (DRV:656-795, mp_wsm6: refl10cm_wsm6 of physics/mp_wsm6_refl.hpp) and compute_relhum
// (DRV:798-850, with RSLF / RSIF of physics/mpas_atmphys_functions.F:156-217, the functions the reference
// executable calls: __mpas_atmphys_functions_MOD_rslf / _rsif from driver_microphysics), DRV:467-482.
// Refused: l_acrain (config_bucket_update), has_reqc / has_reqi / has_reqs = 1 (config_microp_re).
#include "core/field.hpp"
#include "core/precision.hpp"
#include "physics/mp_wsm6.hpp"
#include "physics/mp_wsm6_refl.hpp"

namespace cppmpas::phys {

// The pool arrays driver_microphysics reads or writes for WSM6 (Registry names), and the module variables of
// mpas_atmphys_vars it uses.
struct MicrophysicsPools {
  int nVertLevels = 0;
  // mesh
  FView<const real> zgrid, zz;
  // state, time level 2 (mpas_atm_time_integration.F:2806-2817 calls the driver with time_lev = 2)
  FView<real> rho_zz, theta_m, scalars;  // scalars (num_scalars, nVertLevels, nCells+1)
  FView<const real> w;
  int index_qv = 0, index_qc = 0, index_qr = 0, index_qi = 0, index_qs = 0, index_qg = 0;
  // diag
  FView<real> exner, pressure_p, rtheta_p, surface_pressure;
  FView<const real> pressure_base, exner_base, rtheta_base;
  // tend
  FView<real> rt_diabatic_tend, tend_sfc_pressure;
  // diag_physics
  FView<real> dtheta_dt_mp, precipw, rainnc, rainncv, snownc, snowncv, graupelnc, graupelncv, sr;
  FView<real> rainprod, evapprod, re_cloud, re_ice, re_snow;
  // tend_physics
  FView<real> rthmpten, rqvmpten, rqcmpten, rqrmpten, rqimpten, rqsmpten, rqgmpten;
  int mp_top_level = 0;
  // mpas_atmphys_vars
  real dt_dyn = 0, dt_microp = 0;
  int has_reqc = 0, has_reqi = 0, has_reqs = 0;
  bool l_diags = false, l_acrain = false;
  // with l_diags: diag relhum (nVertLevels, nCells+1), diag_physics refl10cm (nVertLevels, nCells+1),
  // refl10cm_max, refl10cm_1km, refl10cm_1km_max (nCells+1), and mp_wsm6_init's radar tables
  FView<real> relhum, refl10cm, refl10cm_max, refl10cm_1km, refl10cm_1km_max;
  const wsm6::Wsm6Radar* radar = nullptr;
};

// physics/mpas_atmphys_functions.F:156-186 RSLF and :188-218 RSIF (liquid / ice saturation mixing ratio); the
// object code's MAX(-80., T-273.16) is maxsd (x > -80 ? x : -80) and MIN(ESL, P*0.15) minsd (ESL < P*0.15 ? ESL :
// P*0.15), no FMA
CPPMPAS_HD real rslf(real p, real t);
CPPMPAS_HD real rsif(real p, real t);

// DRV:317-500 for mp_wsm6 on cells its..ite. params: module_mp_wsm6::mpas_wsm6_params() (init_microphysics).
void driver_microphysics_wsm6(const wsm6::Wsm6Params<real>& params, MicrophysicsPools& m, int its, int ite);

}  // namespace cppmpas::phys
