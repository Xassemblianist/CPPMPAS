// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// MPAS-Atmosphere's PBL driver for YSU (physics/mpas_atmphys_driver_pbl.F, PB), statement by statement, on the arrays
// of the MPAS pools (Fortran layout, 1-based) and the module arrays of mpas_atmphys_interface.hpp. One rank, one
// thread: its..ite = 1..nCellsSolve, kts..kte = 1..nVertLevels, jts = jte = 1.
//
// driver_pbl = PB:757-973 with config_pbl_scheme = 'bl_ysu': allocate_pbl (PB:87-193), pbl_from_MPAS (PB:306-529), the
// ysu call (PB:827-851, then physics_wrf/module_bl_ysu.F: the copies to its..ite columns around bl_ysu_run and
// rthblten = rthblten_hv / pi3d_hv, module_bl_ysu.F:452; physics/bl_ysu.hpp), pbl_to_MPAS (PB:532-724).
// MPAS passes ctopo = ctopo2 = 1, uoce = voce = 0, flag_bep = .false., idiff = 0, no BEP arrays
// (mpas_atmphys_vars.F:377-382), exch_h = kzh_p; pbl_to_MPAS writes kzq from kzh_p (PB:604) and exch_h from exch_p,
// which stays 0 (PB:414, 623). Refused: 'bl_mynn', and flag_qc / flag_qi = .false. (bl_ysu_run leaves qctnp / qitnp
// unset then).
//
// Reference object code (bin/atmosphere_model, sha256 524a1cf4...): driver_pbl and pbl_from_mpas.isra.0 (pbl_to_MPAS
// is inlined) have no FMA, no ymm register and no libmvec call; their one libm call is MYNN's dx_p pow (PB:452).
// module_bl_ysu.F's ysu is instruction-identical to the standalone build with the reference flags; MPAS's driver file
// has the reference's arithmetic instruction counts (only the set-up of the stub pool calls differs).
#include <string>

#include "core/field.hpp"
#include "core/precision.hpp"
#include "physics/mpas_atmphys_interface.hpp"

namespace cppmpas::phys {

// PblArrays is the part a kernel captures (trivially copyable); PblPools adds the scheme's name.
struct PblArrays {
  int nVertLevels = 0;
  // configs
  bool config_ysu_pblmix = false;
  // sfc_input
  FView<const real> xland;
  // diag_physics: in
  FView<const real> hfx, qfx, ust, wspd, znt, br, fm, fh, u10, v10;
  // diag_physics: inout and out
  FView<real> hpbl, delta, wstar, kzh, kzm, kzq, exch_h;
  FView<int> kpbl;
  // tend_physics
  FView<const real> rthratenlw, rthratensw;
  FView<real> rublten, rvblten, rthblten, rqvblten, rqcblten, rqiblten;
  // mpas_atmphys_vars
  real dt_pbl = 0;                  // mpas_atmphys_manager.F:542 (config_pbl_interval = 'none': dt)
  bool f_qc = true, f_qi = true;    // mpas_atmphys_init.F:501-515
};
struct PblPools : PblArrays {
  std::string pbl_scheme;
};

// PB:757-973 on cells its..ite (its = 1). v: MPAS_to_physics' arrays.
void driver_pbl(int itimestep, PblPools& m, const PhysicsVars& v, int its, int ite);

}  // namespace cppmpas::phys
