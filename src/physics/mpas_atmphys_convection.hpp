// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// MPAS-Atmosphere's convection driver for New Tiedtke (physics/mpas_atmphys_driver_convection.F, CV), statement by
// statement, on the arrays of the MPAS pools (Fortran layout, 1-based) and the module arrays of
// mpas_atmphys_interface.hpp. One rank, one thread: its..ite = 1..nCellsSolve, kts..kte = 1..nVertLevels, jts = jte = 1.
//
// driver_convection = CV:405-621 with config_convection_scheme = 'cu_ntiedtke': allocate_convection (CV:103-265),
// convection_from_MPAS (CV:624-860), initflag (CV:476-478), the cu_ntiedtke_driver call (CV:580-608;
// module_cu_ntiedtke.hpp), convection_to_MPAS (CV:863-988). update_convection_step1 (CV:991-1058) does nothing for
// cu_ntiedtke; update_convection_step2 (CV:1061-1099) accumulates rainc. physics_driver (mpas_atmphys_driver.F:344-369)
// calls update_convection_step1, then driver_convection when l_conv, then update_convection_step2.
// Refused: every other config_convection_scheme ('cu_kain_fritsch', 'cu_tiedtke', 'cu_grell_freitas'), its /= 1,
// more than 100 levels.
//
// Reference object code (bin/atmosphere_model, sha256 524a1cf4...): driver_convection (convection_from_MPAS and
// convection_to_MPAS inlined), update_convection_step1 and update_convection_step2 have no FMA and no ymm register;
// MPAS's file compiled standalone with the reference flags has the reference's arithmetic opcode and libm call
// counts (only the set-up of the stub pool calls differs). Its libm calls: dx_p = len_disp / meshDensity**0.25 of
// the Grell-Freitas (CV:718) and Kain-Fritsch (CV:769) branches, scalar pow@plt, and of the cu_ntiedtke branch
// (CV:840): vectorised over i, _ZGVbN2vv_pow for the pairs of i from its, the scalar pow@plt copies for an odd last i
// and a 1-trip loop (a -g build maps the calls; the loop is entered after overlap and stride tests of dx_p and
// meshDensity that MPAS's arrays pass).
#include <string>

#include "core/field.hpp"
#include "core/precision.hpp"
#include "physics/mpas_atmphys_interface.hpp"

namespace cppmpas::phys {

// ConvectionArrays is the part a kernel captures (trivially copyable); ConvectionPools adds the scheme's name.
struct ConvectionArrays {
  int nVertLevels = 0;
  // configs
  bool config_do_restart = false;
  real config_len_disp = 0, config_bucket_rainc = 0;
  // mesh
  FView<const real> meshDensity;
  // sfc_input
  FView<const real> xland;
  // diag_physics
  FView<const real> hfx, qfx;
  FView<real> cuprec, raincv, rainc;
  FView<int> i_rainc;
  // tend_physics
  FView<const real> rqvdynten, rqvblten, rthdynten, rthblten, rthratenlw, rthratensw;
  FView<real> rthcuten, rqvcuten, rqccuten, rqicuten, rucuten, rvcuten;
  // mpas_atmphys_vars
  real dt_dyn = 0;               // mpas_atmphys_vars.F:161 (dt); cudt = dt_cu / 60 (CV:464) is Kain-Fritsch's
  int n_cu = 1;                  // mpas_atmphys_manager.F:701-702: nint(dt_cu / dt_dyn), at least 1
  bool f_qc = true, f_qi = true; // mpas_atmphys_init.F
  bool l_acrain = false;         // config_bucket_update /= 'none' (mpas_atmphys_manager.F)
  // config_bucket_rainc above: the bucket of update_convection_step2 (CV:1092-1096)
};
struct ConvectionPools : ConvectionArrays {
  std::string convection_scheme;
};

// CV:991-1058 (nothing for cu_ntiedtke; other schemes are refused)
void update_convection_step1(const ConvectionPools& m, int its, int ite);
// CV:405-621 on cells its..ite (its = 1). v: MPAS_to_physics' arrays.
void driver_convection(int itimestep, ConvectionPools& m, const PhysicsVars& v, int its, int ite);
// CV:1061-1099
void update_convection_step2(ConvectionPools& m, int its, int ite);

}  // namespace cppmpas::phys
