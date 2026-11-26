// SPDX-License-Identifier: GPL-3.0-or-later
// MPAS-Atmosphere's Noah path; see mpas_atmphys_driver_lsm.hpp.
#include "physics/mpas_atmphys_driver_lsm.hpp"

#include <algorithm>
#include <string>
#include <vector>

#include "core/error.hpp"
#include "core/libm.hpp"
#include "physics/module_sf_noahdrv.hpp"
#include "physics/mpas_atmphys_constants.hpp"
#include "physics/mpas_atmphys_lsm_noahinit.hpp"

namespace cppmpas::phys {

namespace pc = physics_constants;

namespace {

void check_lsm_config(const LsmConfig& cfg, const LsmPools& p) {
  if (cfg.config_lsm_scheme != "sf_noah")
    throw Error("config_lsm_scheme = '" + cfg.config_lsm_scheme + "': only 'sf_noah' is ported");
  if (p.nSoilLevels < 1 || p.nSoilLevels > noah::kMaxSoilLevels)
    throw Error("nSoilLevels = " + std::to_string(p.nSoilLevels) + ": the Noah port takes 1.." +
                std::to_string(noah::kMaxSoilLevels));
}

// The local arrays of mpas_atmphys_vars that allocate_lsm / allocate_seaice allocate (ims:ime, jms:jme) and
// (ims:ime, 1:num_soils, jms:jme), jms = jme = 1. MPAS's allocate leaves them undefined until *_from_MPAS fills its..ite.
struct Local {
  int ims, ni, ns;
  std::vector<real> dzs_p, smcrel_p, sh2o_p, smois_p, tslb_p;
  std::vector<real> acsnom_p, acsnow_p, canwat_p, chs_p, chs2_p, chklowq_p, cpm_p, cqs2_p, glw_p, grdflx_p, gsw_p,
      hfx_p, lai_p, lh_p, noahres_p, potevp_p, qfx_p, qgh_p, qsfc_p, qz0_p, rainbl_p, br_p, sfc_albbck_p,
      sfc_albedo_p, sfc_emibck_p, sfc_emiss_p, sfcrunoff_p, shdmin_p, shdmax_p, smstav_p, smstot_p, snoalb_p,
      snotime_p, snopcx_p, snow_p, snowc_p, snowh_p, sr_p, swddir_p, swddif_p, swdown_p, tmn_p, tsk_p, udrunoff_p,
      vegfra_p, xice_p, xland_p, z0_p, znt_p, flxsnow_p, fvbsnow_p, fbursnow_p, fgsnsnow_p, frc_urb_p, ust_urb_p;
  std::vector<real> albsi_p, icedepth_p, snowsi_p, tsk_sea, tsk_ice;
  std::vector<int> isltyp_p, ivgtyp_p, utype_urb_p;
  Local(const LsmVars& v) : ims(v.ims), ni(v.ime - v.ims + 1), ns(v.num_soils) {
    dzs_p.assign(ns, 0);
    for (auto* f : {&smcrel_p, &sh2o_p, &smois_p, &tslb_p}) f->assign(std::size_t(ni) * ns, 0);
    for (auto* f : {&acsnom_p, &acsnow_p, &canwat_p, &chs_p, &chs2_p, &chklowq_p, &cpm_p, &cqs2_p, &glw_p, &grdflx_p,
                    &gsw_p, &hfx_p, &lai_p, &lh_p, &noahres_p, &potevp_p, &qfx_p, &qgh_p, &qsfc_p, &qz0_p, &rainbl_p,
                    &br_p, &sfc_albbck_p, &sfc_albedo_p, &sfc_emibck_p, &sfc_emiss_p, &sfcrunoff_p, &shdmin_p,
                    &shdmax_p, &smstav_p, &smstot_p, &snoalb_p, &snotime_p, &snopcx_p, &snow_p, &snowc_p, &snowh_p,
                    &sr_p, &swddir_p, &swddif_p, &swdown_p, &tmn_p, &tsk_p, &udrunoff_p, &vegfra_p, &xice_p, &xland_p,
                    &z0_p, &znt_p, &flxsnow_p, &fvbsnow_p, &fbursnow_p, &fgsnsnow_p, &frc_urb_p, &ust_urb_p,
                    &albsi_p, &icedepth_p, &snowsi_p, &tsk_sea, &tsk_ice})
      f->assign(ni, 0);
    for (auto* f : {&isltyp_p, &ivgtyp_p, &utype_urb_p}) f->assign(ni, 0);
  }
  // (i, j = 1) and (i, n, j = 1), 1-based i in ims..ime
  real& at(std::vector<real>& a, int i) { return a[std::size_t(i - ims)]; }
  int& at(std::vector<int>& a, int i) { return a[std::size_t(i - ims)]; }
  real& at(std::vector<real>& a, int i, int n) { return a[std::size_t(i - ims) + std::size_t(ni) * (n - 1)]; }
};

// mpas_atmphys_driver_lsm.F:268-476
void lsm_from_MPAS(const LsmConfig& cfg, const LsmPools& p, Local& L, int its, int ite) {
  const int num_soils = L.ns;
  for (int n = 1; n <= num_soils; ++n) L.dzs_p[n - 1] = p.dzs(n, its);  // :360-362
  for (int n = 1; n <= num_soils; ++n)                                    // :364-373
    for (int i = its; i <= ite; ++i) {
      L.at(L.sh2o_p, i, n) = p.sh2o(n, i);
      L.at(L.smcrel_p, i, n) = p.smcrel(n, i);
      L.at(L.smois_p, i, n) = p.smois(n, i);
      L.at(L.tslb_p, i, n) = p.tslb(n, i);
    }
  for (int i = its; i <= ite; ++i) {  // :375-440
    L.at(L.acsnom_p, i) = p.acsnom(i);
    L.at(L.acsnow_p, i) = p.acsnow(i);
    L.at(L.canwat_p, i) = p.canwat(i);
    L.at(L.chs_p, i) = p.chs(i);
    L.at(L.chs2_p, i) = p.chs2(i);
    L.at(L.chklowq_p, i) = p.chklowq(i);
    L.at(L.cpm_p, i) = p.cpm(i);
    L.at(L.cqs2_p, i) = p.cqs2(i);
    L.at(L.glw_p, i) = p.glw(i);
    L.at(L.grdflx_p, i) = p.grdflx(i);
    L.at(L.gsw_p, i) = p.gsw(i);
    L.at(L.hfx_p, i) = p.hfx(i);
    L.at(L.lai_p, i) = p.lai(i);
    L.at(L.lh_p, i) = p.lh(i);
    L.at(L.noahres_p, i) = p.noahres(i);
    L.at(L.potevp_p, i) = p.potevp(i);
    L.at(L.qfx_p, i) = p.qfx(i);
    L.at(L.qgh_p, i) = p.qgh(i);
    L.at(L.qsfc_p, i) = p.qsfc(i);
    L.at(L.br_p, i) = p.br(i);
    L.at(L.sfc_albedo_p, i) = p.sfc_albedo(i);
    L.at(L.sfc_emibck_p, i) = p.sfc_emibck(i);
    L.at(L.sfc_emiss_p, i) = p.sfc_emiss(i);
    L.at(L.sfcrunoff_p, i) = p.sfcrunoff(i);
    L.at(L.smstav_p, i) = p.smstav(i);
    L.at(L.smstot_p, i) = p.smstot(i);
    L.at(L.snotime_p, i) = p.snotime(i);
    L.at(L.snopcx_p, i) = p.snopcx(i);
    L.at(L.swddif_p, i) = p.swddif(i);
    L.at(L.swddir_p, i) = p.swddir(i);
    L.at(L.udrunoff_p, i) = p.udrunoff(i);
    L.at(L.z0_p, i) = p.z0(i);
    L.at(L.znt_p, i) = p.znt(i);
    L.at(L.isltyp_p, i) = p.isltyp(i);
    L.at(L.ivgtyp_p, i) = p.ivgtyp(i);
    L.at(L.shdmin_p, i) = p.shdmin(i);
    L.at(L.shdmax_p, i) = p.shdmax(i);
    L.at(L.snoalb_p, i) = p.snoalb(i);
    L.at(L.sfc_albbck_p, i) = p.sfc_albbck(i);
    L.at(L.snow_p, i) = p.snow(i);
    L.at(L.snowc_p, i) = p.snowc(i);
    L.at(L.snowh_p, i) = p.snowh(i);
    L.at(L.tmn_p, i) = p.tmn(i);
    L.at(L.tsk_p, i) = p.skintemp(i);
    L.at(L.vegfra_p, i) = p.vegfra(i);
    L.at(L.xice_p, i) = p.xice(i);
    L.at(L.xland_p, i) = p.xland(i);
    L.at(L.qz0_p, i) = 0.;
    L.at(L.flxsnow_p, i) = 0.;
    L.at(L.fvbsnow_p, i) = 0.;
    L.at(L.fbursnow_p, i) = 0.;
    L.at(L.fgsnsnow_p, i) = 0.;
    L.at(L.frc_urb_p, i) = 0.;
    L.at(L.ust_urb_p, i) = 0.;
    L.at(L.utype_urb_p, i) = 0;
  }
  for (int i = its; i <= ite; ++i) {  // :442-447
    L.at(L.sr_p, i) = 0.;
    L.at(L.rainbl_p, i) = 0.;
  }
  if (cfg.config_microp_scheme != "off")  // :449-459
    for (int i = its; i <= ite; ++i) {
      L.at(L.sr_p, i) = p.sr(i);
      L.at(L.rainbl_p, i) = L.at(L.rainbl_p, i) + p.rainncv(i);
    }
  if (cfg.config_convection_scheme != "off")  // :460-468
    for (int i = its; i <= ite; ++i) L.at(L.rainbl_p, i) = L.at(L.rainbl_p, i) + p.raincv(i);
  for (int i = its; i <= ite; ++i) L.at(L.swdown_p, i) = p.gsw(i) / (1. - p.sfc_albedo(i));  // :470-474
}

// mpas_atmphys_driver_lsm.F:479-633
void lsm_to_MPAS(const LsmConfig& cfg, LsmPools& p, Local& L, int its, int ite) {
  const int num_soils = L.ns;
  for (int n = 1; n <= num_soils; ++n)  // :565-574
    for (int i = its; i <= ite; ++i) {
      p.sh2o(n, i) = L.at(L.sh2o_p, i, n);
      p.smcrel(n, i) = L.at(L.smcrel_p, i, n);
      p.smois(n, i) = L.at(L.smois_p, i, n);
      p.tslb(n, i) = L.at(L.tslb_p, i, n);
    }
  for (int i = its; i <= ite; ++i) {  // :576-621
    p.acsnom(i) = L.at(L.acsnom_p, i);
    p.acsnow(i) = L.at(L.acsnow_p, i);
    p.canwat(i) = L.at(L.canwat_p, i);
    p.chs(i) = L.at(L.chs_p, i);
    p.chs2(i) = L.at(L.chs2_p, i);
    p.chklowq(i) = L.at(L.chklowq_p, i);
    p.cpm(i) = L.at(L.cpm_p, i);
    p.cqs2(i) = L.at(L.cqs2_p, i);
    p.glw(i) = L.at(L.glw_p, i);
    p.grdflx(i) = L.at(L.grdflx_p, i);
    p.gsw(i) = L.at(L.gsw_p, i);
    p.hfx(i) = L.at(L.hfx_p, i);
    p.lai(i) = L.at(L.lai_p, i);
    p.lh(i) = L.at(L.lh_p, i);
    p.noahres(i) = L.at(L.noahres_p, i);
    p.potevp(i) = L.at(L.potevp_p, i);
    p.qfx(i) = L.at(L.qfx_p, i);
    p.qgh(i) = L.at(L.qgh_p, i);
    p.qsfc(i) = L.at(L.qsfc_p, i);
    p.br(i) = L.at(L.br_p, i);
    p.sfc_albedo(i) = L.at(L.sfc_albedo_p, i);
    p.sfc_emibck(i) = L.at(L.sfc_emibck_p, i);
    p.sfc_emiss(i) = L.at(L.sfc_emiss_p, i);
    p.sfcrunoff(i) = L.at(L.sfcrunoff_p, i);
    p.smstav(i) = L.at(L.smstav_p, i);
    p.smstot(i) = L.at(L.smstot_p, i);
    p.snotime(i) = L.at(L.snotime_p, i);
    p.snopcx(i) = L.at(L.snopcx_p, i);
    p.udrunoff(i) = L.at(L.udrunoff_p, i);
    p.z0(i) = L.at(L.z0_p, i);
    p.znt(i) = L.at(L.znt_p, i);
    p.snoalb(i) = L.at(L.snoalb_p, i);
    p.sfc_albbck(i) = L.at(L.sfc_albbck_p, i);
    p.snow(i) = L.at(L.snow_p, i);
    p.snowc(i) = L.at(L.snowc_p, i);
    p.snowh(i) = L.at(L.snowh_p, i);
    p.skintemp(i) = L.at(L.tsk_p, i);
    p.tmn(i) = L.at(L.tmn_p, i);
    p.vegfra(i) = L.at(L.vegfra_p, i);
    p.xice(i) = L.at(L.xice_p, i);
    p.xland(i) = L.at(L.xland_p, i);
  }
  if (cfg.config_microp_scheme != "off")  // :623-631
    for (int i = its; i <= ite; ++i) p.sr(i) = L.at(L.sr_p, i);
}

// mpas_atmphys_lsm_shared.F:22-57
void correct_tsk_over_seaice(int ims, int its, int ite, real xice_thresh, const std::vector<real>& xice,
                             const std::vector<real>& tsk, std::vector<real>& tsk_sea, std::vector<real>& tsk_ice) {
  for (int i = its; i <= ite; ++i) {
    const std::size_t ii = std::size_t(i - ims);
    tsk_sea[ii] = tsk[ii];  // :39-40
    tsk_ice[ii] = tsk[ii];
    if (xice[ii] >= xice_thresh && xice[ii] <= 1.) {  // :42
      tsk_sea[ii] = std::max(tsk_sea[ii], 271.4);  // :44 (maxsd: the second operand when either is NaN)
      if (xice[ii] < 0.2 && tsk_ice[ii] < 253.15) tsk_ice[ii] = 253.15;  // :48
      if (xice[ii] < 0.1 && tsk_ice[ii] < 263.15) tsk_ice[ii] = 263.15;  // :49
    }
  }
}

// mpas_atmphys_driver_seaice.F:178-333
void seaice_from_MPAS(const LsmConfig& cfg, const LsmPools& p, const LsmVars& v, Local& L, int its, int ite) {
  const int num_soils = L.ns;
  for (int i = its; i <= ite; ++i) {  // :238-282
    L.at(L.xice_p, i) = p.xice(i);
    L.at(L.glw_p, i) = p.glw(i);
    L.at(L.qgh_p, i) = p.qgh(i);
    L.at(L.snoalb_p, i) = p.snoalb(i);
    L.at(L.br_p, i) = p.br(i);
    L.at(L.chs_p, i) = p.chs(i);
    L.at(L.swdown_p, i) = p.gsw(i) / (1. - p.sfc_albedo(i));
    for (int n = 1; n <= num_soils; ++n) L.at(L.tslb_p, i, n) = p.tslb(n, i);
    L.at(L.z0_p, i) = p.z0(i);
    L.at(L.snow_p, i) = p.snow(i);
    L.at(L.snowc_p, i) = p.snowc(i);
    L.at(L.snowh_p, i) = p.snowh(i);
    L.at(L.tsk_p, i) = p.skintemp(i);
    L.at(L.cqs2_p, i) = p.cqs2(i);
    L.at(L.acsnom_p, i) = p.acsnom(i);
    L.at(L.acsnow_p, i) = p.acsnow(i);
    L.at(L.sfcrunoff_p, i) = p.sfcrunoff(i);
    L.at(L.albsi_p, i) = 0.65;   // seaice_albedo_default (mpas_atmphys_vars.F:651)
    L.at(L.snowsi_p, i) = 0.001;  // seaice_snowdepth_min (:654)
    L.at(L.icedepth_p, i) = 3.0;  // seaice_thickness_default (:652)
    L.at(L.hfx_p, i) = p.hfx(i);
    L.at(L.lh_p, i) = p.lh(i);
    L.at(L.qfx_p, i) = p.qfx(i);
    L.at(L.znt_p, i) = p.znt(i);
    L.at(L.grdflx_p, i) = p.grdflx(i);
    L.at(L.qsfc_p, i) = p.qsfc(i);
    L.at(L.chs2_p, i) = p.chs2(i);
    if (p.xice(i) >= v.xice_threshold && p.xice(i) <= 1.) {  // :275-281
      L.at(L.sfc_albedo_p, i) = (p.sfc_albedo(i) - 0.08 * (1. - p.xice(i))) / p.xice(i);
      L.at(L.sfc_emiss_p, i) = (p.sfc_emiss(i) - 0.98 * (1. - p.xice(i))) / p.xice(i);
    } else {
      L.at(L.sfc_emiss_p, i) = p.sfc_emiss(i);
      L.at(L.sfc_albedo_p, i) = p.sfc_albedo(i);
    }
  }
  correct_tsk_over_seaice(L.ims, its, ite, v.xice_threshold, L.xice_p, L.tsk_p, L.tsk_sea, L.tsk_ice);  // :285-286
  for (int i = its; i <= ite; ++i) L.at(L.tsk_p, i) = L.at(L.tsk_ice, i);  // :287-289
  for (int i = its; i <= ite; ++i) {  // :293-296
    L.at(L.sr_p, i) = 0.;
    L.at(L.rainbl_p, i) = 0.;
  }
  if (cfg.config_microp_scheme != "off")  // :297-304
    for (int i = its; i <= ite; ++i) {
      L.at(L.sr_p, i) = p.sr(i);
      L.at(L.rainbl_p, i) = L.at(L.rainbl_p, i) + p.rainncv(i);
    }
  if (cfg.config_convection_scheme != "off")  // :305-310
    for (int i = its; i <= ite; ++i) L.at(L.rainbl_p, i) = L.at(L.rainbl_p, i) + p.raincv(i);
  for (int i = its; i <= ite; ++i) {  // :313-329 (sf_noah)
    L.at(L.noahres_p, i) = p.noahres(i);
    L.at(L.potevp_p, i) = p.potevp(i);
    L.at(L.snopcx_p, i) = p.snopcx(i);
  }
}

// mpas_atmphys_driver_seaice.F:336-460
void seaice_to_MPAS(LsmPools& p, const LsmVars& v, Local& L, int its, int ite) {
  const int num_soils = L.ns;
  auto sea = [&](const FView<const real>& f, int i) { return f(i - v.ims + 1); };  // (ims:ime, jms:jme) of the sfclayer
  for (int i = its; i <= ite; ++i) {  // :387-405
    const real x = L.at(L.xice_p, i);
    if (x >= v.xice_threshold && x <= 1.) {
      p.cpm(i) = x * p.cpm(i) + (1. - x) * sea(v.cpm_sea, i);
      L.at(L.chs_p, i) = x * L.at(L.chs_p, i) + (1. - x) * sea(v.chs_sea, i);
      L.at(L.chs2_p, i) = x * L.at(L.chs2_p, i) + (1. - x) * sea(v.chs2_sea, i);
      L.at(L.cqs2_p, i) = x * L.at(L.cqs2_p, i) + (1. - x) * sea(v.cqs2_sea, i);
      L.at(L.hfx_p, i) = x * L.at(L.hfx_p, i) + (1. - x) * sea(v.hfx_sea, i);
      L.at(L.lh_p, i) = x * L.at(L.lh_p, i) + (1. - x) * sea(v.lh_sea, i);
      L.at(L.qfx_p, i) = x * L.at(L.qfx_p, i) + (1. - x) * sea(v.qfx_sea, i);
      L.at(L.qgh_p, i) = x * L.at(L.qgh_p, i) + (1. - x) * sea(v.qgh_sea, i);
      L.at(L.qsfc_p, i) = x * L.at(L.qsfc_p, i) + (1. - x) * sea(v.qsfc_sea, i);
      L.at(L.tsk_p, i) = x * L.at(L.tsk_p, i) + (1. - x) * L.at(L.tsk_sea, i);
      L.at(L.sfc_albedo_p, i) = x * L.at(L.sfc_albedo_p, i) + (1. - x) * 0.08;
      L.at(L.sfc_emiss_p, i) = x * L.at(L.sfc_emiss_p, i) + (1. - x) * 0.98;
    }
  }
  for (int i = its; i <= ite; ++i) {  // :408-438
    for (int n = 1; n <= num_soils; ++n) p.tslb(n, i) = L.at(L.tslb_p, i, n);
    p.z0(i) = L.at(L.z0_p, i);
    p.snow(i) = L.at(L.snow_p, i);
    p.snowc(i) = L.at(L.snowc_p, i);
    p.snowh(i) = L.at(L.snowh_p, i);
    p.skintemp(i) = L.at(L.tsk_p, i);
    p.acsnom(i) = L.at(L.acsnom_p, i);
    p.acsnow(i) = L.at(L.acsnow_p, i);
    p.sfcrunoff(i) = L.at(L.sfcrunoff_p, i);
    p.znt(i) = L.at(L.znt_p, i);
    p.grdflx(i) = L.at(L.grdflx_p, i);
    p.chs(i) = L.at(L.chs_p, i);
    p.chs2(i) = L.at(L.chs2_p, i);
    p.cqs2(i) = L.at(L.cqs2_p, i);
    p.qsfc(i) = L.at(L.qsfc_p, i);
    p.qgh(i) = L.at(L.qgh_p, i);
    p.hfx(i) = L.at(L.hfx_p, i);
    p.qfx(i) = L.at(L.qfx_p, i);
    p.lh(i) = L.at(L.lh_p, i);
    p.sfc_albedo(i) = L.at(L.sfc_albedo_p, i);
    p.sfc_emiss(i) = L.at(L.sfc_emiss_p, i);
  }
  for (int i = its; i <= ite; ++i) {  // :440-456 (sf_noah)
    p.noahres(i) = L.at(L.noahres_p, i);
    p.potevp(i) = L.at(L.potevp_p, i);
    p.snopcx(i) = L.at(L.snopcx_p, i);
  }
}

}  // namespace

// mpas_atmphys_driver_lsm.F:636-664 -> mpas_atmphys_lsm_noahinit.F noah_init_forMPAS / lsminit (:76-207)
void init_lsm(const LsmConfig& cfg, LsmPools& p, const std::string& table_dir, noah::NoahTables& tables) {
  check_lsm_config(cfg, p);
  tables = noah::NoahTables{};
  soil_veg_gen_parm(table_dir, p.mminlu, cfg.input_soil_data, tables);  // :120
  if (!cfg.config_do_restart) {                                            // :124
    int errflag = 0;
    for (int iCell = 1; iCell <= p.nCells; ++iCell) {  // :127-135
      if (p.isltyp(iCell) < 1) errflag = 1;
      if (!cfg.config_sfc_snowalbedo) p.snoalb(iCell) = tables.maxalb[p.ivgtyp(iCell) - 1] * 0.01;
    }
    if (errflag == 1)  // :136-138
      throw Error("module_sf_noahlsm.F: lsminit: out of range value of ISLTYP. Is this field in the input?");
    constexpr real blim = 5.5, hlice = 3.335e5, t0 = 273.15;  // :101-103
    for (int iCell = 1; iCell <= p.nCells; ++iCell) {  // :140-170
      const int st = p.isltyp(iCell) - 1;
      real bx = tables.bb[st];
      real smcmax = tables.maxsmc[st];
      real psisat = tables.satpsi[st];
      if ((bx > 0.0) && (smcmax > 0.0) && (psisat > 0.0)) {
        for (int ns = 1; ns <= p.nSoilLevels; ++ns) {
          if (p.tslb(ns, iCell) < 273.149) {
            bx = tables.bb[st];
            smcmax = tables.maxsmc[st];
            psisat = tables.satpsi[st];
            if (bx > blim) bx = blim;
            // fk = (((hlice/(gravity*(-psisat))) * ((tslb-t0)/tslb))**(-1/bx))*smcmax (:153-154)
            real fk = libm::pow((hlice / (pc::gravity * (-psisat))) * ((p.tslb(ns, iCell) - t0) / p.tslb(ns, iCell)),
                                -(1. / bx)) *
                      smcmax;
            if (fk < 0.02) fk = 0.02;
            p.sh2o(ns, iCell) = std::min(fk, p.smois(ns, iCell));  // :156 (minsd operand order as MIN(fk, smois))
            real free = 0;
            noah::NoahKernel::frh2o(free, p.tslb(ns, iCell), p.smois(ns, iCell), p.sh2o(ns, iCell), smcmax, bx,
                                    psisat);  // :157-158
            p.sh2o(ns, iCell) = free;
          } else {
            p.sh2o(ns, iCell) = p.smois(ns, iCell);
          }
        }
      } else {
        for (int ns = 1; ns <= p.nSoilLevels; ++ns) p.sh2o(ns, iCell) = p.smois(ns, iCell);
      }
    }
    // fndsnowh = .true. (:96): the snow-height branch (:172-178) is not taken
  }
}

// mpas_atmphys_driver_lsm.F:667-761
void driver_lsm(int itimestep, const LsmConfig& cfg, const noah::NoahTables& tables, LsmPools& p, const LsmVars& v,
                int its, int ite) {
  check_lsm_config(cfg, p);
  if (v.num_soils != p.nSoilLevels) throw Error("driver_lsm: num_soils /= nSoilLevels");
  Local L(v);  // allocate_lsm (:122-192)
  lsm_from_MPAS(cfg, p, L, its, ite);  // :697
  noahdrv::LsmArgs a;                  // :705-749
  a.dz8w = v.dz_p.p; a.p8w3d = v.pres2_hyd_p.p; a.t3d = v.t_p.p; a.qv3d = v.qv_p.p;
  a.xland = L.xland_p.data(); a.xice = L.xice_p.data(); a.ivgtyp = L.ivgtyp_p.data(); a.isltyp = L.isltyp_p.data();
  a.tmn = L.tmn_p.data(); a.vegfra = L.vegfra_p.data(); a.shdmin = L.shdmin_p.data(); a.shdmax = L.shdmax_p.data();
  a.snoalb = L.snoalb_p.data(); a.glw = L.glw_p.data(); a.swdown = L.swdown_p.data(); a.rainbl = L.rainbl_p.data();
  a.embck = L.sfc_emibck_p.data(); a.sr = L.sr_p.data(); a.qgh = L.qgh_p.data(); a.qz0 = L.qz0_p.data();
  a.tsk = L.tsk_p.data(); a.hfx = L.hfx_p.data(); a.qfx = L.qfx_p.data(); a.lh = L.lh_p.data();
  a.grdflx = L.grdflx_p.data(); a.qsfc = L.qsfc_p.data(); a.cqs2 = L.cqs2_p.data(); a.chs = L.chs_p.data();
  a.chs2 = L.chs2_p.data(); a.snow = L.snow_p.data(); a.snowc = L.snowc_p.data(); a.snowh = L.snowh_p.data();
  a.canwat = L.canwat_p.data(); a.smstav = L.smstav_p.data(); a.smstot = L.smstot_p.data();
  a.sfcrunoff = L.sfcrunoff_p.data(); a.udrunoff = L.udrunoff_p.data(); a.acsnom = L.acsnom_p.data();
  a.acsnow = L.acsnow_p.data(); a.snotime = L.snotime_p.data(); a.snopcx = L.snopcx_p.data();
  a.emiss = L.sfc_emiss_p.data(); a.rib = L.br_p.data(); a.potevp = L.potevp_p.data();
  a.albedo = L.sfc_albedo_p.data(); a.albbck = L.sfc_albbck_p.data(); a.z0 = L.z0_p.data(); a.znt = L.znt_p.data();
  a.lai = L.lai_p.data(); a.noahres = L.noahres_p.data(); a.chklowq = L.chklowq_p.data();
  a.sh2o = L.sh2o_p.data(); a.smois = L.smois_p.data(); a.tslb = L.tslb_p.data(); a.smcrel = L.smcrel_p.data();
  a.dzs = L.dzs_p.data(); a.isurban = v.isurban; a.isice = p.isice;
  a.dt = v.dt_pbl; a.myj = false; a.itimestep = itimestep; a.frpcpn = false; a.rdlai2d = false;  // :31-34
  a.opt_thcnd = 1; a.ua_phys = false; a.fasdas = 0;  // mpas_atmphys_vars.F:821-830
  a.num_soil_layers = v.num_soils; a.xice_threshold = v.xice_threshold; a.usemonalb = cfg.config_sfc_albedo;
  a.sf_urban_physics = 0;  // :37
  a.ims = v.ims; a.ime = v.ime; a.its = its; a.ite = ite; a.kms = v.kms; a.kts = v.kms;
  noahdrv::lsm(tables, a);
  lsm_to_MPAS(cfg, p, L, its, ite);  // :757
}

// mpas_atmphys_driver_seaice.F:463-557
void driver_seaice(const LsmConfig& cfg, LsmPools& p, const LsmVars& v, int its, int ite) {
  check_lsm_config(cfg, p);
  if (v.num_soils != p.nSoilLevels) throw Error("driver_seaice: num_soils /= nSoilLevels");
  Local L(v);  // allocate_seaice (:36-101)
  seaice_from_MPAS(cfg, p, v, L, its, ite);  // :484
  noahdrv::SeaiceArgs a;                     // :488-516
  a.dz8w = v.dz_p.p; a.p8w3d = v.pres2_hyd_p.p; a.t3d = v.t_p.p; a.qv3d = v.qv_p.p;
  a.xice = L.xice_p.data(); a.snoalb2d = L.snoalb_p.data(); a.glw = L.glw_p.data(); a.swdown = L.swdown_p.data();
  a.rainbl = L.rainbl_p.data(); a.sr = L.sr_p.data(); a.qgh = L.qgh_p.data(); a.tsk = L.tsk_p.data();
  a.hfx = L.hfx_p.data(); a.qfx = L.qfx_p.data(); a.lh = L.lh_p.data(); a.grdflx = L.grdflx_p.data();
  a.qsfc = L.qsfc_p.data(); a.emiss = L.sfc_emiss_p.data(); a.albedo = L.sfc_albedo_p.data(); a.rib = L.br_p.data();
  a.cqs2 = L.cqs2_p.data(); a.chs = L.chs_p.data(); a.chs2 = L.chs2_p.data(); a.z02d = L.z0_p.data();
  a.znt = L.znt_p.data(); a.tslb = L.tslb_p.data(); a.snow = L.snow_p.data(); a.snowc = L.snowc_p.data();
  a.snowh2d = L.snowh_p.data(); a.acsnow = L.acsnow_p.data(); a.acsnom = L.acsnom_p.data();
  a.sfcrunoff = L.sfcrunoff_p.data(); a.dt = v.dt_pbl; a.frpcpn = false;  // driver_seaice.F:20
  a.noahres = L.noahres_p.data(); a.potevp = L.potevp_p.data(); a.snopcx = L.snopcx_p.data();
  // mpas_atmphys_vars.F:637-654 (parameters)
  a.seaice_albedo_opt = 0; a.seaice_albedo_default = 0.65; a.seaice_thickness_opt = 0;
  a.seaice_thickness_default = 3.0; a.seaice_snowdepth_opt = 0; a.seaice_snowdepth_max = 1.e10;
  a.seaice_snowdepth_min = 0.001;
  a.xice_threshold = v.xice_threshold; a.num_soil_layers = v.num_soils; a.sf_urban_physics = 0;
  a.ims = v.ims; a.ime = v.ime; a.its = its; a.ite = ite; a.kms = v.kms; a.kts = v.kms;
  noahdrv::seaice_noah(a);
  seaice_to_MPAS(p, v, L, its, ite);  // :553
}

// mpas_atmphys_sfc_diagnostics.F:27-121, config_lsm_scheme = 'sf_noah' (:63-79)
void atmphys_sfc_diagnostics(const LsmConfig& cfg, LsmPools& p, int its, int ite) {
  (void)its;
  (void)ite;
  if (cfg.config_lsm_scheme != "sf_noah")
    throw Error("atmphys_sfc_diagnostics: config_lsm_scheme = '" + cfg.config_lsm_scheme + "' is not ported");
  for (int i = 1; i <= p.nCellsSolve; ++i) {
    const real rho = p.surface_pressure(i) / (pc::R_d * p.skintemp(i));
    if (p.cqs2(i) < 1.e-5) {
      p.q2(i) = p.qsfc(i);
    } else {
      p.q2(i) = p.qsfc(i) - p.qfx(i) / (rho * p.cqs2(i));
    }
    if (p.chs2(i) < 1.e-5) {
      p.t2m(i) = p.skintemp(i);
    } else {
      p.t2m(i) = p.skintemp(i) - p.hfx(i) / (rho * pc::cp * p.chs2(i));
    }
    p.th2m(i) = p.t2m(i) * libm::pow(pc::P0 / p.surface_pressure(i), pc::rcp);
  }
}

}  // namespace cppmpas::phys
