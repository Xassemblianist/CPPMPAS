// SPDX-License-Identifier: GPL-3.0-or-later
// MPAS-Atmosphere's Noah path; see mpas_atmphys_driver_lsm.hpp.
#include "physics/mpas_atmphys_driver_lsm.hpp"

#include <algorithm>
#include <string>
#include <vector>

#include "core/error.hpp"
#include "core/libm.hpp"
#include "physics/module_sf_noahdrv.hpp"
#include "physics/mpas_atmphys_interface.hpp"
#include "physics/mpas_atmphys_constants.hpp"
#include "physics/mpas_atmphys_lsm_noahinit.hpp"

namespace cppmpas::phys {

namespace pc = physics_constants;

namespace {

void check_lsm_config(const LsmConfig& cfg, const LsmPools& p) {
  if (cfg.config_lsm_scheme != "sf_noah")
    throw Error("config_lsm_scheme = '" + cfg.config_lsm_scheme + "': only 'sf_noah' is ported");
  // one layer: SHFLX / SMFLX / the ice and glacial HRT read layer 2 of zsoil, stc, smc (undefined in MPAS)
  if (p.nSoilLevels < 2 || p.nSoilLevels > noah::kMaxSoilLevels)
    throw Error("nSoilLevels = " + std::to_string(p.nSoilLevels) + ": the Noah port takes 2.." +
                std::to_string(noah::kMaxSoilLevels));
}

// The local arrays of mpas_atmphys_vars that allocate_lsm / allocate_seaice allocate (ims:ime, jms:jme) and
// (ims:ime, 1:num_soils, jms:jme), jms = jme = 1. MPAS's allocate leaves them undefined until *_from_MPAS fills its..ite.
// Work arrays a kernel can write (WorkArray); LocalView is their kernel view with the same names and at().
#define CPPMPAS_LSM_SOIL(X) X(smcrel_p) X(sh2o_p) X(smois_p) X(tslb_p)
#define CPPMPAS_LSM_REAL(X) \
  X(acsnom_p) X(acsnow_p) X(canwat_p) X(chs_p) X(chs2_p) X(chklowq_p) X(cpm_p) X(cqs2_p) X(glw_p) X(grdflx_p) \
  X(gsw_p) X(hfx_p) X(lai_p) X(lh_p) X(noahres_p) X(potevp_p) X(qfx_p) X(qgh_p) X(qsfc_p) X(qz0_p) X(rainbl_p) \
  X(br_p) X(sfc_albbck_p) X(sfc_albedo_p) X(sfc_emibck_p) X(sfc_emiss_p) X(sfcrunoff_p) X(shdmin_p) X(shdmax_p) \
  X(smstav_p) X(smstot_p) X(snoalb_p) X(snotime_p) X(snopcx_p) X(snow_p) X(snowc_p) X(snowh_p) X(sr_p) \
  X(swddir_p) X(swddif_p) X(swdown_p) X(tmn_p) X(tsk_p) X(udrunoff_p) X(vegfra_p) X(xice_p) X(xland_p) X(z0_p) \
  X(znt_p) X(flxsnow_p) X(fvbsnow_p) X(fbursnow_p) X(fgsnsnow_p) X(frc_urb_p) X(ust_urb_p) X(albsi_p) \
  X(icedepth_p) X(snowsi_p) X(tsk_sea) X(tsk_ice) 
#define CPPMPAS_LSM_INT(X) X(isltyp_p) X(ivgtyp_p) X(utype_urb_p)
struct LocalView {
  int ims, ni, ns;
  real* dzs_p;
#define CPPMPAS_X(a) real* a;
  CPPMPAS_LSM_SOIL(CPPMPAS_X) CPPMPAS_LSM_REAL(CPPMPAS_X)
#undef CPPMPAS_X
#define CPPMPAS_X(a) int* a;
  CPPMPAS_LSM_INT(CPPMPAS_X)
#undef CPPMPAS_X
  // (i, j = 1) and (i, n, j = 1), 1-based i in ims..ime
  CPPMPAS_HD real& at(real* a, int i) const { return a[std::size_t(i - ims)]; }
  CPPMPAS_HD int& at(int* a, int i) const { return a[std::size_t(i - ims)]; }
  CPPMPAS_HD real& at(real* a, int i, int n) const { return a[std::size_t(i - ims) + std::size_t(ni) * (n - 1)]; }
};
struct Local {
  int ims, ni, ns;
  WorkArray dzs_p;
#define CPPMPAS_X(a) WorkArray a;
  CPPMPAS_LSM_SOIL(CPPMPAS_X) CPPMPAS_LSM_REAL(CPPMPAS_X)
#undef CPPMPAS_X
#define CPPMPAS_X(a) WorkArrayInt a;
  CPPMPAS_LSM_INT(CPPMPAS_X)
#undef CPPMPAS_X
  Local(const LsmVars& v) : ims(v.ims), ni(v.ime - v.ims + 1), ns(v.num_soils) {
    dzs_p.assign(ns, 0);
#define CPPMPAS_X(a) a.assign(std::size_t(ni) * ns, 0);
    CPPMPAS_LSM_SOIL(CPPMPAS_X)
#undef CPPMPAS_X
#define CPPMPAS_X(a) a.assign(ni, 0);
    CPPMPAS_LSM_REAL(CPPMPAS_X) CPPMPAS_LSM_INT(CPPMPAS_X)
#undef CPPMPAS_X
  }
  LocalView view() {
    LocalView w;
    w.ims = ims;
    w.ni = ni;
    w.ns = ns;
    w.dzs_p = dzs_p.data();
#define CPPMPAS_X(a) w.a = a.data();
    CPPMPAS_LSM_SOIL(CPPMPAS_X) CPPMPAS_LSM_REAL(CPPMPAS_X) CPPMPAS_LSM_INT(CPPMPAS_X)
#undef CPPMPAS_X
    return w;
  }
};

// mpas_atmphys_driver_lsm.F:268-476
// One kernel per loop (one thread per cell i; the n loops inside the thread).
void lsm_from_MPAS(const LsmConfig& cfg, const LsmPools& pools, Local& loc, int its, int ite) {
  const LsmArrays p = pools;
  const LocalView L = loc.view();
  const int num_soils = L.ns;
  const bool microp = cfg.config_microp_scheme != "off", convection = cfg.config_convection_scheme != "off";
  serial([=] CPPMPAS_LAMBDA() {
    for (int n = 1; n <= num_soils; ++n) L.dzs_p[n - 1] = p.dzs(n, its);  // :360-362
  });
  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {  // :364-373
    for (int n = 1; n <= num_soils; ++n) {
      L.at(L.sh2o_p, i, n) = p.sh2o(n, i);
      L.at(L.smcrel_p, i, n) = p.smcrel(n, i);
      L.at(L.smois_p, i, n) = p.smois(n, i);
      L.at(L.tslb_p, i, n) = p.tslb(n, i);
    }
  });
  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {  // :375-440
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
  });
  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {  // :442-447
    L.at(L.sr_p, i) = 0.;
    L.at(L.rainbl_p, i) = 0.;
  });
  if (microp)  // :449-459
    parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {
      L.at(L.sr_p, i) = p.sr(i);
      L.at(L.rainbl_p, i) = L.at(L.rainbl_p, i) + p.rainncv(i);
    });
  if (convection)  // :460-468
    parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) { L.at(L.rainbl_p, i) = L.at(L.rainbl_p, i) + p.raincv(i); });
  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {  // :470-474
    L.at(L.swdown_p, i) = p.gsw(i) / (1. - p.sfc_albedo(i));
  });
}

// mpas_atmphys_driver_lsm.F:479-633
void lsm_to_MPAS(const LsmConfig& cfg, LsmPools& pools, Local& loc, int its, int ite) {
  const LsmArrays p = pools;
  const LocalView L = loc.view();
  const int num_soils = L.ns;
  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {  // :565-574
    for (int n = 1; n <= num_soils; ++n) {
      p.sh2o(n, i) = L.at(L.sh2o_p, i, n);
      p.smcrel(n, i) = L.at(L.smcrel_p, i, n);
      p.smois(n, i) = L.at(L.smois_p, i, n);
      p.tslb(n, i) = L.at(L.tslb_p, i, n);
    }
  });
  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {  // :576-621
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
  });
  if (cfg.config_microp_scheme != "off")  // :623-631
    parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) { p.sr(i) = L.at(L.sr_p, i); });
}

// mpas_atmphys_lsm_shared.F:22-57
void correct_tsk_over_seaice(int ims, int its, int ite, real xice_thresh, const real* xice, const real* tsk,
                             real* tsk_sea, real* tsk_ice) {
  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {
    const std::size_t ii = std::size_t(i - ims);
    tsk_sea[ii] = tsk[ii];  // :39-40
    tsk_ice[ii] = tsk[ii];
    if (xice[ii] >= xice_thresh && xice[ii] <= 1.) {  // :42
      tsk_sea[ii] = std::max(tsk_sea[ii], 271.4);  // :44 (maxsd: the second operand when either is NaN)
      if (xice[ii] < 0.2 && tsk_ice[ii] < 253.15) tsk_ice[ii] = 253.15;  // :48
      if (xice[ii] < 0.1 && tsk_ice[ii] < 263.15) tsk_ice[ii] = 263.15;  // :49
    }
  });
}

// mpas_atmphys_driver_seaice.F:178-333
void seaice_from_MPAS(const LsmConfig& cfg, const LsmPools& pools, const LsmVars& v, Local& loc, int its, int ite) {
  const LsmArrays p = pools;
  const LocalView L = loc.view();
  const int num_soils = L.ns;
  const bool microp = cfg.config_microp_scheme != "off", convection = cfg.config_convection_scheme != "off";
  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {  // :238-282
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
  });
  correct_tsk_over_seaice(L.ims, its, ite, v.xice_threshold, L.xice_p, L.tsk_p, L.tsk_sea, L.tsk_ice);  // :285-286
  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) { L.at(L.tsk_p, i) = L.at(L.tsk_ice, i); });  // :287-289
  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {  // :293-296
    L.at(L.sr_p, i) = 0.;
    L.at(L.rainbl_p, i) = 0.;
  });
  if (microp)  // :297-304
    parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {
      L.at(L.sr_p, i) = p.sr(i);
      L.at(L.rainbl_p, i) = L.at(L.rainbl_p, i) + p.rainncv(i);
    });
  if (convection)  // :305-310
    parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) { L.at(L.rainbl_p, i) = L.at(L.rainbl_p, i) + p.raincv(i); });
  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {  // :313-329 (sf_noah)
    L.at(L.noahres_p, i) = p.noahres(i);
    L.at(L.potevp_p, i) = p.potevp(i);
    L.at(L.snopcx_p, i) = p.snopcx(i);
  });
}

// mpas_atmphys_driver_seaice.F:336-460
void seaice_to_MPAS(LsmPools& pools, const LsmVars& v, Local& loc, int its, int ite) {
  const LsmArrays p = pools;
  const LocalView L = loc.view();
  const int num_soils = L.ns;
  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {  // :387-405
    auto sea = [&](const FView<const real>& f, int ii) { return f(ii - v.ims + 1); };  // (ims:ime) of the sfclayer
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
  });
  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {  // :408-438
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
  });
  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {  // :440-456 (sf_noah)
    p.noahres(i) = L.at(L.noahres_p, i);
    p.potevp(i) = L.at(L.potevp_p, i);
    p.snopcx(i) = L.at(L.snopcx_p, i);
  });
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
      // :129 reads MAXALB(IVGTYP) with no check of the category: outside 1..LUCATS MPAS reads out of bounds
      if (!cfg.config_sfc_snowalbedo && (p.ivgtyp(iCell) < 1 || p.ivgtyp(iCell) > tables.lucats))
        throw Error("module_sf_noahlsm.F: lsminit: IVGTYP outside the vegetation table (MAXALB(IVGTYP) out of bounds)");
      if (!cfg.config_sfc_snowalbedo) p.snoalb(iCell) = tables.maxalb[p.ivgtyp(iCell) - 1] * 0.01;
    }
    if (errflag == 1)  // :136-138
      throw Error("module_sf_noahlsm.F: lsminit: out of range value of ISLTYP. Is this field in the input?");
    // :140-170 read BB, MAXSMC, SATPSI(ISLTYP) with no upper check: a category above SLCATS is out of bounds
    for (int iCell = 1; iCell <= p.nCells; ++iCell)
      if (p.isltyp(iCell) > tables.slcats)
        throw Error("module_sf_noahlsm.F: lsminit: ISLTYP above the soil table (BB(ISLTYP) out of bounds)");
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
  // seaice_from_MPAS / seaice_to_MPAS read the open-water values of the surface layer (*_sea) of every sea-ice cell
  // (:275-281, 387-405); physics_driver fills them with config_frac_seaice = true only (false, xice_threshold 0.5,
  // is not ported: check_physics_supported refuses it)
  for (const FView<const real>* f : {&v.cpm_sea, &v.chs_sea, &v.chs2_sea, &v.cqs2_sea, &v.hfx_sea, &v.qfx_sea,
                                     &v.qgh_sea, &v.qsfc_sea, &v.lh_sea})
    if (!f->p) throw Error("driver_seaice: the surface layer's open-water arrays (*_sea) are not set");
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
void atmphys_sfc_diagnostics(const LsmConfig& cfg, LsmPools& pools, int its, int ite) {
  (void)its;
  (void)ite;
  if (cfg.config_lsm_scheme != "sf_noah")
    throw Error("atmphys_sfc_diagnostics: config_lsm_scheme = '" + cfg.config_lsm_scheme + "' is not ported");
  const LsmArrays p = pools;
  parallel_for(1, p.nCellsSolve, [=] CPPMPAS_LAMBDA(int i) {
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
  });
}

}  // namespace cppmpas::phys
