// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// MPAS-Atmosphere's surface layer driver for the revised MM5 scheme (physics/mpas_atmphys_driver_sfclayer.F, SL),
// statement by statement, on the arrays of the MPAS pools (Fortran layout, 1-based) and the module arrays of
// mpas_atmphys_interface.hpp. One rank, one thread: its..ite = 1..nCellsSolve, jts = jte = 1.
//
// driver_sfclayer = SL:839-1088 with config_sfclayer_scheme = 'sf_monin_obukhov_rev' or 'sf_mynn': allocate_sfclayer
// (SL:97-227), sfclayer_from_MPAS (SL:352-592), the one or two (config_frac_seaice) scheme calls, sfclayer_to_MPAS
// (SL:595-796).
//  * sf_monin_obukhov_rev: the sfclayrev calls (SL:951-1013, then physics_wrf/module_sf_sfclayrev.F and
//    sf_sfclayrev_pre.F, the copies around sf_sfclayrev_run, and physics/sf_sfclayrev.hpp);
//  * sf_mynn: the sfclay_mynn calls (SL:1014-1077, then physics_wrf/module_sf_mynn.F, the copies around
//    sf_mynn_pre_run and sf_mynn_run, and physics/sf_mynn.hpp) with itimestep = initflag (SL:885-886).
// init_sfclayer (SL:799-836) is init_sfclayrev_tables_o3 for both (sf_mynn_init's psi_init builds the same tables
// with the same vector code, sf_mynn.hpp). Refused: 'sf_monin_obukhov' (module_sf_sfclay).
//
// Reference object code (bin/atmosphere_model, sha256 524a1cf4...): driver_sfclayer and sfclayer_from_mpas.isra.0
// (sfclayer_to_MPAS is inlined) have no FMA, no ymm register and no libmvec call; the one libm call is the scalar
// pow@plt of dx_p = len_disp / meshDensity**0.25 (SL:425). MPAS's file compiled standalone with the reference flags has
// the reference's arithmetic instruction counts (only the set-up of the stub pool calls differs);
// module_sf_sfclayrev.F's sfclayrev and sf_sfclayrev_pre_run are instruction-identical, and so are module_sf_mynn.F's
// sfclay_mynn and sf_mynn_pre_run (tests/sfmynn_column): copies only, no libm call; sf_mynn_pre_run's
// itimestep = 1 loop (sf_mynn_pre.F:153-158) is vectorised with sqrtpd / maxpd / divpd, which round as the scalar
// statements do.
#include <string>
#include <vector>

#include "core/field.hpp"
#include "core/precision.hpp"
#include "physics/mpas_atmphys_interface.hpp"
#include "physics/sf_sfclayrev.hpp"

namespace cppmpas::phys {

// SfclayerArrays is the part a kernel captures (trivially copyable); SfclayerPools adds the scheme's name.
struct SfclayerArrays {
  // configs
  bool config_frac_seaice = true, config_do_restart = false;
  real config_len_disp = 0;
  // mesh
  FView<const real> meshDensity, areaCell;
  // sfc_input
  FView<const real> skintemp, xland, sst, xice;
  // diag_physics: in
  FView<const real> hpbl, mavail;
  // diag_physics: inout and out
  FView<real> br, cpm, chs, chs2, cqs2, flhc, flqc, gz1oz0, hfx, qfx, qgh, qsfc, lh, mol, psim, psih, regime, rmol,
      ust, ustm, wspd, znt, zol, q2, t2m, th2m, u10, v10, cd, cda, ck, cka, fh, fm;
  // sf_mynn only (SL:568-587, 773-790): diag_physics qcg (in) and ch (inout), sfc_input snowh (in)
  FView<const real> qcg, snowh;
  FView<real> ch;
  // mpas_atmphys_vars
  real xice_threshold = 0;  // mpas_atmphys_manager.F:740-744
};
struct SfclayerPools : SfclayerArrays {
  std::string sfclayer_scheme;
};

// init_sfclayer (SL:799-836) for 'sf_monin_obukhov_rev' (sf_sfclayrev_init) and 'sf_mynn' (sf_mynn_init), as the MPAS
// executable computes them: the same four psi tables
void init_sfclayer(const std::string& sfclayer_scheme, sfclayrev::SFClayRevTables<real>& tables);

// The *_sea arrays of mpas_atmphys_vars (mpas_atmphys_vars.F:934-937) that outlive driver_sfclayer: allocate_sfclayer
// allocates them, deallocate_sfclayer leaves them, driver_seaice reads them and deallocate_seaice frees them
// (mpas_atmphys_driver_seaice.F:152-160, 336-460). Filled only with config_frac_seaice (1-based: element i - 1).
struct SfclayerSeaArrays {
  WorkArray chs, chs2, cqs2, cpm, hfx, qfx, qgh, qsfc, lh;
};

// SL:839-1088 on cells its..ite (its = 1). v: MPAS_to_physics' arrays (pres_hyd_p, psfc_p, t_p, u_p, v_p, qv_p, dz_p;
// sf_mynn also rho_p and qc_p, and levels 1 and 2 of u_p, v_p and dz_p). tables: init_sfclayer's.
// sea_out (optional): receives the *_sea arrays the sea-ice driver reads.
void driver_sfclayer(int itimestep, const sfclayrev::SFClayRevTables<real>& tables, SfclayerPools& m,
                     const PhysicsVars& v, int its, int ite, SfclayerSeaArrays* sea_out = nullptr);

}  // namespace cppmpas::phys
