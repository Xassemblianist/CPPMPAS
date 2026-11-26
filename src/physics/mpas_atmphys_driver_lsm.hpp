// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// MPAS-Atmosphere's Noah land-surface path (config_lsm_scheme = 'sf_noah'), statement by statement, on the arrays of
// the MPAS pools (Fortran layout: (iCell) and (nSoilLevels, iCell), 1-based). One rank, one thread: its..ite =
// 1..nCellsSolve, jts = jte = 1, kts = 1 (mpas_atmphys_manager.F:653-663).
//
//  * init_lsm (mpas_atmphys_driver_lsm.F:636-664) -> noah_init_forMPAS -> lsminit and soil_veg_gen_parm
//    (mpas_atmphys_lsm_noahinit.F, physics/mpas_atmphys_lsm_noahinit.cpp);
//  * driver_lsm (mpas_atmphys_driver_lsm.F:667-761): allocate_lsm (:122-192), lsm_from_MPAS (:268-476), lsm of
//    physics_wrf/module_sf_noahdrv.F with MPAS's arguments (:705-749; physics/module_sf_noahdrv.cpp),
//    lsm_to_MPAS (:479-633);
//  * driver_seaice (mpas_atmphys_driver_seaice.F:463-557): allocate_seaice (:36-101), seaice_from_MPAS
//    (:178-333) with correct_tsk_over_seaice (mpas_atmphys_lsm_shared.F:22-57), seaice_noah
//    (module_sf_noah_seaice_drv.F, :488-516), seaice_to_MPAS (:336-460);
//  * atmphys_sfc_diagnostics (mpas_atmphys_sfc_diagnostics.F:27-121) for sf_noah,
// in the order of physics_driver (mpas_atmphys_driver.F:283-315).
// Refused (hard errors): config_lsm_scheme other than 'sf_noah' (sf_noahmp is not ported), more than
// noah::kMaxSoilLevels soil levels. With config_do_restart lsminit only reads the tables (:124); the Noah state then
// comes from the restart stream, whose reading is not part of this port.
#include <string>

#include "core/field.hpp"
#include "core/precision.hpp"
#include "physics/sf_noah.hpp"

namespace cppmpas::phys {

// The configs the Noah path reads (Registry defaults: config_sfc_albedo, config_sfc_snowalbedo = true,
// input_soil_data = 'STAS')
struct LsmConfig {
  std::string config_lsm_scheme = "sf_noah";
  std::string config_microp_scheme = "mp_wsm6";
  std::string config_convection_scheme = "cu_ntiedtke";
  std::string input_soil_data = "STAS";
  bool config_sfc_albedo = true;      // usemonalb
  bool config_sfc_snowalbedo = true;  // lsminit's input_sfc_albedo
  bool config_do_restart = false;
};

// The pool arrays the Noah path reads or writes (Registry names; (iCell) or (nSoilLevels, iCell))
struct LsmPools {
  int nCells = 0, nCellsSolve = 0, nSoilLevels = 0;
  // sfc_input
  std::string mminlu;
  int isice = 0;
  FView<int> isltyp, ivgtyp;
  FView<real> shdmin, shdmax, snoalb, sfc_albbck, snow, snowc, snowh, tmn, skintemp, vegfra, xice, xland;
  FView<real> dzs, sh2o, smcrel, smois, tslb;
  // diag_physics
  FView<real> acsnom, acsnow, canwat, chs, chs2, chklowq, cpm, cqs2, glw, grdflx, gsw, hfx, lai, lh, noahres, potevp,
      qfx, qgh, qsfc, br, sfc_albedo, sfc_emibck, sfc_emiss, sfcrunoff, smstav, smstot, snotime, snopcx, swddif, swddir,
      udrunoff, z0, znt, sr, rainncv, raincv, q2, t2m, th2m;
  // diag
  FView<real> surface_pressure;
};

// The module variables of mpas_atmphys_vars the Noah path uses: the dimensions (physics_run_init), the physics
// arrays of MPAS_to_physics (mpas_atmphys_interface.F; (ims:ime, kms:kme), i fastest, jms = jme = 1) and the
// surface layer's open-water values of the sea-ice cells (*_sea, mpas_atmphys_driver_sfclayer.F)
struct LsmVars {
  int ims = 1, ime = 0, kms = 1, kme = 0;
  int num_soils = 0;
  int isurban = 0;  // mpas_atmphys_landuse.F (landuse_init_forMPAS)
  real dt_pbl = 0;
  real xice_threshold = 0;  // mpas_atmphys_manager.F:740-744
  FView<const real> dz_p, pres2_hyd_p, t_p, qv_p;
  FView<const real> cpm_sea, chs_sea, chs2_sea, cqs2_sea, hfx_sea, qfx_sea, qgh_sea, qsfc_sea, lh_sea;
};

// mpas_atmphys_driver_lsm.F:636-664 (init_lsm); reads VEGPARM.TBL, SOILPARM.TBL and GENPARM.TBL in table_dir
void init_lsm(const LsmConfig& cfg, LsmPools& p, const std::string& table_dir, noah::NoahTables& tables);
// mpas_atmphys_driver_lsm.F:667-761 (with allocate_lsm / deallocate_lsm around it, mpas_atmphys_driver.F:285-292)
void driver_lsm(int itimestep, const LsmConfig& cfg, const noah::NoahTables& tables, LsmPools& p, const LsmVars& v,
                int its, int ite);
// mpas_atmphys_driver_seaice.F:463-557 (with allocate_seaice / deallocate_seaice, mpas_atmphys_driver.F:303-310)
void driver_seaice(const LsmConfig& cfg, LsmPools& p, const LsmVars& v, int its, int ite);
// mpas_atmphys_sfc_diagnostics.F:27-121
void atmphys_sfc_diagnostics(const LsmConfig& cfg, LsmPools& p, int its, int ite);

}  // namespace cppmpas::phys
