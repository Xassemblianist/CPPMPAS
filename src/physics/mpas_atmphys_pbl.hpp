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
// which stays 0 (PB:414, 623). Refused: flag_qc / flag_qi = .false. (bl_ysu_run leaves qctnp / qitnp unset then).
//
// driver_pbl with config_pbl_scheme = 'bl_mynn' (PB:853-969): allocate_pbl's bl_mynn arrays (PB:136-182),
// pbl_from_MPAS (PB:306-529 with its bl_mynn branch PB:423-506: dx_p = config_len_disp / meshDensity**0.25, the
// scalar pow@plt), initflag (PB:820-821: 1 at the first step of a run that is not a restart), mynn_bl_driver
// (physics_wrf/module_bl_mynn.F with bl_mynn_pre / bl_mynn_post and MMM-physics bl_mynn_run: physics/bl_mynn.hpp,
// one thread per column) and pbl_to_MPAS (PB:532-724 with its bl_mynn branch PB:627-706: kzq from kzh_p, rqsblten,
// rniblten under f_ni). tke_pbl goes to tkepbl_p and back unchanged (mynn_bl_driver does not take it); delta and
// wstar are YSU's. init_pbl (PB:727-754) is mynn_init_pbl: bl_mynn_init with mpas_atmphys_constants. Refused: the
// config_mynn_* values without a verified port (mynn::Options::refused), f_qc / f_qi / f_qs / f_ni = .false. (MPAS
// has qc, qi, qs and ni with bl_mynn: Registry.xml packages bl_mynn_in), f_nc / f_nifa / f_nwfa (the
// aerosol-aware Thompson scheme), f_qoz and f_nbca (not in the Registry), more than mynn::kMaxLevels levels.
//
// Reference object code (bin/atmosphere_model, sha256 524a1cf4...): driver_pbl and pbl_from_mpas.isra.0 (pbl_to_MPAS
// is inlined) have no FMA, no ymm register and no libmvec call; their one libm call is MYNN's dx_p pow (PB:452).
// module_bl_ysu.F's ysu is instruction-identical to the standalone build with the reference flags; MPAS's driver file
// has the reference's arithmetic instruction counts (only the set-up of the stub pool calls differs).
#include <string>

#include "core/field.hpp"
#include "core/precision.hpp"
#include "physics/bl_mynn_common.hpp"
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
  // bl_mynn (PB:316-338, 551-562): configs, mesh, sfc_input, diag_physics and tend_physics
  bool config_do_restart = false;
  real config_len_disp = 0;
  mynn::Options mynn;                         // config_mynn_* and config_do_DAcycling (PB:854-870)
  mynn::Const mynn_const{};                   // bl_mynn_common after init_pbl (mynn_init_pbl)
  bool f_qs = true, f_qoz = false, f_nc = false, f_ni = true, f_nifa = false, f_nwfa = false, f_nbca = false;
  FView<const real> meshDensity, skintemp, ch, qsfc, rmol;
  FView<real> el_pbl, cov, qke, qke_adv, qsq, tsq, tke_pbl, sh3d, sm3d, cldfrac_bl, qc_bl, qi_bl, edmf_a, edmf_ent,
      edmf_qc, edmf_qt, edmf_thl, edmf_w, sub_thl, sub_qv, det_thl, det_qv, dqke, qbuoy, qdiss, qshear, qwt;
  FView<real> rqsblten, rniblten;
  // not pool arrays: when set, driver_pbl copies mpas_atmphys_vars' kbl_plume_p, maxwidthbl_p, maxmfbl_p and
  // zbl_plume_p (allocate_pbl's, written by mynn_bl_driver, never read by MPAS) there, for the column check
  FView<int> diag_kbl_plume;
  FView<real> diag_maxwidthbl, diag_maxmfbl, diag_zbl_plume;
};
struct PblPools : PblArrays {
  std::string pbl_scheme;
};

// PB:757-973 on cells its..ite (its = 1). v: MPAS_to_physics' arrays.
void driver_pbl(int itimestep, PblPools& m, const PhysicsVars& v, int its, int ite);

// PB:727-754 init_pbl for 'bl_mynn': bl_mynn_init(cp, cpv, cice, cliq, ep_1, ep_2, gravity, karman, P0, R_d, R_v,
// svp1, svp2, svp3, svpt0, xlf, xls, xlv) of mpas_atmphys_constants
mynn::Const mynn_init_pbl();

}  // namespace cppmpas::phys
