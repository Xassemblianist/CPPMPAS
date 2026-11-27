// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// MPAS-Atmosphere's gravity wave drag over orography for the YSU suite (physics/mpas_atmphys_driver_gwdo.F, GD, and
// physics/physics_wrf/module_bl_gwdo.F, BG), statement by statement, on the arrays of the MPAS pools (Fortran layout,
// 1-based) and the module arrays of mpas_atmphys_interface.hpp. One rank, one thread: its..ite = 1..nCellsSolve,
// kts..kte = 1..nVertLevels, jts = jte = 1.
//
// driver_gwdo = GD:692-834 with config_gwdo_scheme = 'bl_ysu_gwdo': allocate_gwdo (GD:75-177), gwdo_from_MPAS
// (GD:285-552: the static subgrid orography var2d, con, oa1-4, ol1-4 of sfc_input, sina = 0, cosa = 1,
// dx_p = len_disp / meshDensity**0.25, kpbl, dusfcg / dvsfcg, dtaux3d / dtauy3d and the PBL tendencies), the gwdo call
// (GD:757-777, then BG:16-238: the _hv copies of its..ite around bl_gwdo_run, physics/bl_gwdo.hpp, with
// dxmeter = dx_p, deltim = dt_pbl, kme = kte + 1), gwdo_to_MPAS (GD:555-689: dusfcg, dvsfcg, dtaux3d, dtauy3d,
// rubldiff / rvbldiff = the drag, and rublten / rvblten / rthblten). Refused: 'bl_ugwp_gwdo' and config_ngw_scheme.
//
// Reference object code (bin/atmosphere_model, sha256 524a1cf4...): driver_gwdo (gwdo_from_MPAS and gwdo_to_MPAS
// inlined) has no FMA, no ymm register and no libmvec call; its one libm call is the scalar pow@plt of dx_p (GD:533;
// a -g build maps it). MPAS's file compiled standalone with the reference flags has the reference's arithmetic opcode
// and libm call counts (only the set-up of the stub pool calls differs). module_bl_gwdo.F's gwdo is
// instruction-identical to the standalone build with the reference flags (-Dmpas), no arithmetic.
#include <string>

#include "core/field.hpp"
#include "core/precision.hpp"
#include "physics/mpas_atmphys_interface.hpp"

namespace cppmpas::phys {

// GwdoArrays is the part a kernel captures (trivially copyable); GwdoPools adds the scheme's name.
struct GwdoArrays {
  int nVertLevels = 0;
  // configs
  bool config_ngw_scheme = false;
  real config_len_disp = 0;
  // mesh
  FView<const real> meshDensity;
  // sfc_input
  FView<const real> var2d, con, oa1, oa2, oa3, oa4, ol1, ol2, ol3, ol4;
  // diag_physics
  FView<const int> kpbl;
  FView<real> dusfcg, dvsfcg, dtaux3d, dtauy3d, rubldiff, rvbldiff;
  // tend_physics
  FView<real> rublten, rvblten, rthblten;
  // mpas_atmphys_vars
  real dt_pbl = 0;  // mpas_atmphys_manager.F:542 (config_pbl_interval = 'none': dt)
};
struct GwdoPools : GwdoArrays {
  std::string gwdo_scheme;
};

// GD:692-834 on cells its..ite (its = 1). v: MPAS_to_physics' arrays.
void driver_gwdo(int itimestep, GwdoPools& m, const PhysicsVars& v, int its, int ite);

}  // namespace cppmpas::phys
