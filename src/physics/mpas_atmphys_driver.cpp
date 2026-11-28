// SPDX-License-Identifier: GPL-3.0-or-later
// See mpas_atmphys_driver.hpp. CORE = mpas_atm_core.F, CT = physics/mpas_atmphys_control.F, MG =
// physics/mpas_atmphys_manager.F, IN = physics/mpas_atmphys_init.F, DR = physics/mpas_atmphys_driver.F,
// TI = dynamics/mpas_atm_time_integration.F (MPAS v8.4.2).
#include "physics/mpas_atmphys_driver.hpp"

#include <cmath>
#include <cstdlib>
#include <limits>

#include "core/error.hpp"
#include "core/libm.hpp"
#include "physics/module_mp_wsm6.hpp"
#include "physics/mpas_atmphys_constants.hpp"
#include "physics/mpas_atmphys_convection.hpp"
#include "physics/mpas_atmphys_driver_lsm.hpp"
#include "physics/mpas_atmphys_gwdo.hpp"
#include "physics/mpas_atmphys_landuse.hpp"
#include "physics/mpas_atmphys_microphysics.hpp"
#include "physics/mpas_atmphys_pbl.hpp"
#include "physics/mpas_atmphys_todynamics.hpp"

namespace cppmpas::phys {

namespace {

template <class T> FView<const T> cv(const FArray<T>& a) { return a.view(); }

// A time interval string of the namelist ("DD_hh:mm:ss", "hh:mm:ss", "mm:ss" or "ss"; mpas_set_timeInterval) in
// whole seconds. Fractions are refused.
long long interval_seconds(const std::string& t) {
  std::string s = t;
  long long days = 0;
  const auto us = s.find('_');
  if (us != std::string::npos) {
    days = std::atoll(s.substr(0, us).c_str());
    s = s.substr(us + 1);
  }
  long long v[3] = {0, 0, 0};
  int n = 0;
  std::size_t pos = 0;
  while (pos <= s.size() && n < 3) {
    const auto c = s.find(':', pos);
    const std::string part = s.substr(pos, c == std::string::npos ? std::string::npos : c - pos);
    if (part.empty() || part.find_first_not_of("0123456789") != std::string::npos)
      throw Error("time interval '" + t + "' is not [DD_]hh:mm:ss in whole seconds");
    v[n++] = std::atoll(part.c_str());
    if (c == std::string::npos) break;
    pos = c + 1;
  }
  long long sec = 0;
  if (n == 3) sec = v[0] * 3600 + v[1] * 60 + v[2];
  else if (n == 2) sec = v[0] * 60 + v[1];
  else sec = v[0];
  sec += days * 86400;
  if (sec <= 0) throw Error("time interval '" + t + "' is not positive");
  return sec;
}

int num(const std::string& ts, std::size_t p, std::size_t n) { return std::atoi(ts.substr(p, n).c_str()); }

// mpas_get_time(..., DoY=): the day of the year of a YYYY-MM-DD_hh:mm:ss stamp (gregorian)
int day_of_year(const std::string& ts) {
  const int y = num(ts, 0, 4), mo = num(ts, 5, 2), d = num(ts, 8, 2);
  static const int md[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
  const bool leap = (y % 4 == 0 && y % 100 != 0) || y % 400 == 0;
  int doy = d;
  for (int k = 1; k < mo; ++k) doy += (k == 2 && leap) ? 29 : md[k - 1];
  return doy;
}

int scalar_index(const State& s, const char* n) {
  for (std::size_t k = 0; k < s.scalar_names.size(); ++k)
    if (s.scalar_names[k] == n) return static_cast<int>(k) + 1;
  return -1;  // MPAS's index of an inactive constituent (mpas_atm_core_interface.F)
}

RadiationPools radiation_pools(PhysicsRun& r, const Mesh& m, RegistryPools& p) {
  RadiationPools q;
  q.nVertLevels = m.nVertLevels;
  q.nOznLevels = p.dim("nOznLevels", false);
  q.nMonths = p.dim("nMonths", false);
  q.radt_lw_scheme = r.cfg.config_radt_lw_scheme;
  q.radt_sw_scheme = r.cfg.config_radt_sw_scheme;
  q.radt_cld_scheme = r.cfg.config_radt_cld_scheme;
  q.radt_cld_overlap = r.cfg.config_radt_cld_overlap;
  q.radt_cld_dcorrlen = r.cfg.config_radt_cld_dcorrlen;
  q.microp_scheme = r.cfg.config_microp_scheme;
  q.config_o3climatology = r.cfg.config_o3climatology;
  q.config_microp_re = r.cfg.config_microp_re;
  q.latCell = cv(m.latCell);
  q.lonCell = cv(m.lonCell);
  q.skintemp = cv(p.R("skintemp"));
  q.snow = cv(p.R("snow"));
  q.xice = cv(p.R("xice"));
  q.xland = cv(p.R("xland"));
  q.pin = p.v("pin");
  q.ozmixm = p.v("ozmixm");
  q.sfc_albedo = cv(p.R("sfc_albedo"));
  q.sfc_emiss = cv(p.R("sfc_emiss"));
#define CPPMPAS_P(x) q.x = p.v(#x);
  CPPMPAS_P(cldfrac) CPPMPAS_P(o3clim) CPPMPAS_P(o3vmr) CPPMPAS_P(rre_cloud) CPPMPAS_P(rre_ice) CPPMPAS_P(rre_snow)
  CPPMPAS_P(glw) CPPMPAS_P(lwcf) CPPMPAS_P(lwdnb) CPPMPAS_P(lwdnbc) CPPMPAS_P(lwdnt) CPPMPAS_P(lwdntc)
  CPPMPAS_P(lwupb) CPPMPAS_P(lwupbc) CPPMPAS_P(lwupt) CPPMPAS_P(lwuptc) CPPMPAS_P(olrtoa)
  CPPMPAS_P(coszr) CPPMPAS_P(gsw) CPPMPAS_P(swcf) CPPMPAS_P(swdnb) CPPMPAS_P(swdnbc) CPPMPAS_P(swdnt)
  CPPMPAS_P(swdntc) CPPMPAS_P(swupb) CPPMPAS_P(swupbc) CPPMPAS_P(swupt) CPPMPAS_P(swuptc) CPPMPAS_P(swddir)
  CPPMPAS_P(swddni) CPPMPAS_P(swddif) CPPMPAS_P(rthratenlw) CPPMPAS_P(rthratensw)
#undef CPPMPAS_P
  return q;
}

LsmConfig lsm_config(const PhysicsRun& r) {
  LsmConfig c;
  c.config_lsm_scheme = r.cfg.config_lsm_scheme;
  c.config_microp_scheme = r.cfg.config_microp_scheme;
  c.config_convection_scheme = r.cfg.config_convection_scheme;
  c.input_soil_data = r.cfg.input_soil_data;
  c.config_sfc_albedo = r.cfg.config_sfc_albedo;
  c.config_sfc_snowalbedo = r.cfg.config_sfc_snowalbedo;
  c.config_do_restart = r.cfg.config_do_restart;
  return c;
}

LsmPools lsm_pools(const Mesh& m, State& s, RegistryPools& p) {
  LsmPools q;
  q.nCells = m.nCells;
  q.nCellsSolve = m.nCellsSolve();
  q.nSoilLevels = p.dim("nSoilLevels", false);
  std::string lu = p.C("mminlu");
  while (!lu.empty() && lu.back() == ' ') lu.pop_back();
  q.mminlu = lu;
  q.isice = p.I("isice_lu")(1);
  q.isltyp = p.vi("isltyp");
  q.ivgtyp = p.vi("ivgtyp");
#define CPPMPAS_P(x) q.x = p.v(#x);
  CPPMPAS_P(shdmin) CPPMPAS_P(shdmax) CPPMPAS_P(snoalb) CPPMPAS_P(sfc_albbck) CPPMPAS_P(snow) CPPMPAS_P(snowc)
  CPPMPAS_P(snowh) CPPMPAS_P(tmn) CPPMPAS_P(skintemp) CPPMPAS_P(vegfra) CPPMPAS_P(xice) CPPMPAS_P(xland)
  CPPMPAS_P(dzs) CPPMPAS_P(sh2o) CPPMPAS_P(smcrel) CPPMPAS_P(smois) CPPMPAS_P(tslb)
  CPPMPAS_P(acsnom) CPPMPAS_P(acsnow) CPPMPAS_P(canwat) CPPMPAS_P(chs) CPPMPAS_P(chs2) CPPMPAS_P(chklowq)
  CPPMPAS_P(cpm) CPPMPAS_P(cqs2) CPPMPAS_P(glw) CPPMPAS_P(grdflx) CPPMPAS_P(gsw) CPPMPAS_P(hfx) CPPMPAS_P(lai)
  CPPMPAS_P(lh) CPPMPAS_P(noahres) CPPMPAS_P(potevp) CPPMPAS_P(qfx) CPPMPAS_P(qgh) CPPMPAS_P(qsfc) CPPMPAS_P(br)
  CPPMPAS_P(sfc_albedo) CPPMPAS_P(sfc_emibck) CPPMPAS_P(sfc_emiss) CPPMPAS_P(sfcrunoff) CPPMPAS_P(smstav)
  CPPMPAS_P(smstot) CPPMPAS_P(snotime) CPPMPAS_P(snopcx) CPPMPAS_P(swddif) CPPMPAS_P(swddir) CPPMPAS_P(udrunoff)
  CPPMPAS_P(z0) CPPMPAS_P(znt) CPPMPAS_P(sr) CPPMPAS_P(rainncv) CPPMPAS_P(raincv) CPPMPAS_P(q2) CPPMPAS_P(t2m)
  CPPMPAS_P(th2m)
#undef CPPMPAS_P
  q.surface_pressure = s.surface_pressure.view();  // diag surface_pressure (State in CPPMPAS)
  return q;
}

FView<const real> pv(const PArray& a, int kme) { return FView<const real>{a.data(), a.ldi, kme, 1}; }
FView<const real> sv(const WorkArray& a) {
  return FView<const real>{a.data(), static_cast<int>(a.size()), 1, 1};
}

ConvectionPools convection_pools(PhysicsRun& r, const Mesh& m, RegistryPools& p, const DycorePhysicsArrays& a) {
  ConvectionPools q;
  q.nVertLevels = m.nVertLevels;
  q.convection_scheme = r.cfg.config_convection_scheme;
  q.config_do_restart = r.cfg.config_do_restart;
  q.config_len_disp = r.cfg.config_len_disp;
  q.config_bucket_rainc = r.cfg.config_bucket_rainc;
  q.meshDensity = cv(m.meshDensity);
  q.xland = cv(p.R("xland"));
  q.hfx = cv(p.R("hfx"));
  q.qfx = cv(p.R("qfx"));
  q.cuprec = p.v("cuprec");
  q.raincv = p.v("raincv");
  q.rainc = p.v("rainc");
  q.i_rainc = p.vi("i_rainc");
  q.rqvdynten = cv(p.R("rqvdynten"));
  q.rqvblten = cv(p.R("rqvblten"));
  q.rthdynten = cv(*a.rthdynten);
  q.rthblten = cv(p.R("rthblten"));
  q.rthratenlw = cv(p.R("rthratenlw"));
  q.rthratensw = cv(p.R("rthratensw"));
  q.rthcuten = p.v("rthcuten");
  q.rqvcuten = p.v("rqvcuten");
  q.rqccuten = p.v("rqccuten");
  q.rqicuten = p.v("rqicuten");
  q.rucuten = p.v("rucuten");
  q.rvcuten = p.v("rvcuten");
  q.dt_dyn = r.dt_dyn;
  q.n_cu = r.n_cu;
  q.f_qc = r.f_qc;
  q.f_qi = r.f_qi;
  q.l_acrain = r.l_acrain;
  if (r.cfg.config_convection_scheme == "cu_grell_freitas") {  // CV:427-433, 655-720: cu_grell_freitas_in
    q.config_gfconv_closure_deep = r.cfg.config_gfconv_closure_deep;
    q.config_gfconv_closure_shallow = r.cfg.config_gfconv_closure_shallow;
    q.areaCell = cv(m.areaCell);
    q.zgrid = cv(m.zgrid);
    q.gsw = cv(p.R("gsw"));
    q.kpbl = cv(p.I("kpbl"));
#define CPPMPAS_P(x) q.x = p.v(#x);
    CPPMPAS_P(cubot) CPPMPAS_P(cutop) CPPMPAS_P(xmb_total) CPPMPAS_P(xmb_shallow) CPPMPAS_P(qc_cu) CPPMPAS_P(qi_cu)
#undef CPPMPAS_P
    q.k22_shallow = p.vi("k22_shallow");
    q.kbcon_shallow = p.vi("kbcon_shallow");
    q.ktop_shallow = p.vi("ktop_shallow");
    q.ktop_deep = p.vi("ktop_deep");
  }
  return q;
}

// mynn::Options of the namelist (PB:854-870: mynn_bl_driver's bl_mynn_* arguments) and config_do_DAcycling (PB:813)
mynn::Options mynn_options(const PhysicsConfig& c) {
  mynn::Options o;
  o.tkeadvect = c.config_mynn_tkeadvect;
  o.tkebudget = c.config_mynn_tkebudget;
  o.cloudpdf = c.config_mynn_cloudpdf;
  o.mixlength = c.config_mynn_mixlength;
  o.stfunc = c.config_mynn_stfunc;
  o.topdown = c.config_mynn_topdown;
  o.scaleaware = c.config_mynn_scaleaware;
  o.dheat_opt = c.config_mynn_dheat_opt;
  o.edmf = c.config_mynn_edmf;
  o.edmf_dd = c.config_mynn_edmf_dd;
  o.edmf_mom = c.config_mynn_edmf_mom;
  o.edmf_tke = c.config_mynn_edmf_tke;
  o.output = c.config_mynn_edmf_output;
  o.mixscalars = c.config_mynn_mixscalars;
  o.cloudmix = c.config_mynn_mixclouds;
  o.mixqt = c.config_mynn_mixqt;
  o.closure = c.config_mynn_closure;
  o.do_DAcycling = c.config_do_DAcycling;
  return o;  // icloud_bl = spp_pbl = 0 (mpas_atmphys_vars.F:418-419)
}

}  // namespace

// ---------------------------------------------------------------------------------------------------------------
void physics_core_init(PhysicsRun& r, const Mesh& m, State& s, Diag& d, RegistryPools& p,
                       const std::string& clock_now, std::vector<std::string>& log) {
  (void)d;
  const PhysicsConfig& c = r.cfg;
  check_physics_supported(c);
  if (!c.moist_physics) return;  // CORE:581 (the suite 'none': no physics at all)
  {  // CORE:569-578
    std::string& mminlu = p.C("mminlu");
    if (mminlu.find_first_not_of(' ') == std::string::npos) {
      log.push_back("No information on land use dataset is available.");
      log.push_back("Assume that we are using 'USGS'.");
      mminlu = "USGS";
      mminlu.resize(64, ' ');
    }
  }
  // CT:393-440 physics_tables_init: its test (/= mp_thompson .or. /= mp_thompson_aerosols) returns at once
  // CT:335-390 physics_registry_init
  if (!c.config_do_restart && (c.config_lsm_scheme == "sf_noah" || c.config_lsm_scheme == "sf_noahmp")) {
    auto& dzs = p.R("dzs");
    if (dzs.extent(0) < 4) throw Error("physics_registry_init: nSoilLevels < 4");
    for (int iCell = 1; iCell <= m.nCells; ++iCell) {
      dzs(1, iCell) = 0.10;
      dzs(2, iCell) = 0.30;
      dzs(3, iCell) = 0.60;
      dzs(4, iCell) = 1.00;
    }
  }
  // ---- MG:378-766 physics_run_init
  const real dt = c.config_dt;  // the clock's time step
  r.start_time = c.config_start_time.substr(0, 19);
  {
    const std::string& st = r.start_time;
    physics_radiation_time_init(r.rad, num(st, 0, 4), day_of_year(st), num(st, 11, 2), num(st, 14, 2),
                                num(st, 17, 2), 0, 1);  // MG:455-461
  }
  auto alarm = [&](const std::string& interval, ClockAlarm& al, real& dt_out) {  // MG:467-546
    if (interval != "none") {
      al.defined = true;
      al.ring = 0;  // alarmStartTime = startTime
      al.interval = interval_seconds(interval);
      dt_out = static_cast<real>(al.interval);
    } else {
      dt_out = dt;
    }
  };
  alarm(c.config_radtlw_interval, r.radtlw_alarm, r.rad.dt_radtlw);
  alarm(c.config_radtsw_interval, r.radtsw_alarm, r.rad.dt_radtsw);
  {
    ClockAlarm conv, pbl;  // config_conv_interval / config_pbl_interval = 'none' (check_physics_supported)
    alarm(c.config_conv_interval, conv, r.dt_cu);
    alarm(c.config_pbl_interval, pbl, r.dt_pbl);
  }
  r.green_alarm.defined = true;  // MG:549-555
  r.green_alarm.ring = 0;
  r.green_alarm.interval = interval_seconds(c.config_greeness_update);
  // MG:620-645 the diagnostics alarm: defined by the caller from the output streams (r.diag_alarm)
  log.push_back("DT_RADTLW = " + std::to_string(r.rad.dt_radtlw));
  log.push_back("DT_RADTSW = " + std::to_string(r.rad.dt_radtsw));
  log.push_back("DT_CU     = " + std::to_string(r.dt_cu));
  log.push_back("DT_PBL    = " + std::to_string(r.dt_pbl));
  r.dt_dyn = dt;  // MG:680-694
  microp_timestep(c.config_microp_scheme, r.dt_dyn, r.dt_microp, r.n_microp);
  r.n_cu = std::max(static_cast<int>(std::lround(r.dt_cu / r.dt_dyn)), 1);  // MG:701-702 (nint)
  r.l_conv = r.l_radtlw = r.l_radtsw = r.l_diags = r.l_acrain = r.l_acradt = false;
  r.xice_threshold = c.config_frac_seaice ? 0.02 : 0.5;  // MG:740-744
  r.has_reqc = r.has_reqi = r.has_reqs = 0;            // MG:746-760 (config_microp_re = false)

  // ---- IN:83-481 physics_init
  const int julday = day_of_year(clock_now);  // IN:180-181
  {  // IN:556-596 init_dirs_forphys
    auto &east = p.R("east"), &north = p.R("north");
    for (int iCell = 1; iCell <= m.nCells; ++iCell) {
      const real lon = m.lonCell(iCell), lat = m.latCell(iCell);
      east(1, iCell) = -libm::sin(lon);
      east(2, iCell) = libm::cos(lon);
      east(3, iCell) = 0.0;
      {  // r3_normalize: mi = 1.0 / sqrt(ax**2 + ay**2 + az**2)
        const real mi = 1.0 / std::sqrt(east(1, iCell) * east(1, iCell) + east(2, iCell) * east(2, iCell) +
                                         east(3, iCell) * east(3, iCell));
        east(1, iCell) = east(1, iCell) * mi;
        east(2, iCell) = east(2, iCell) * mi;
        east(3, iCell) = east(3, iCell) * mi;
      }
      north(1, iCell) = -(libm::cos(lon) * libm::sin(lat));
      north(2, iCell) = -(libm::sin(lon) * libm::sin(lat));
      north(3, iCell) = libm::cos(lat);
      {
        const real mi = 1.0 / std::sqrt(north(1, iCell) * north(1, iCell) + north(2, iCell) * north(2, iCell) +
                                         north(3, iCell) * north(3, iCell));
        north(1, iCell) = north(1, iCell) * mi;
        north(2, iCell) = north(2, iCell) * mi;
        north(3, iCell) = north(3, iCell) * mi;
      }
    }
  }
  // IN:484-545 init_physics_flags
  r.index_qv = scalar_index(s, "qv");
  r.index_qc = scalar_index(s, "qc");
  r.index_qr = scalar_index(s, "qr");
  r.index_qi = scalar_index(s, "qi");
  r.index_qs = scalar_index(s, "qs");
  r.index_qg = scalar_index(s, "qg");
  r.f_qc = r.index_qc > -1;
  r.f_qr = r.index_qr > -1;
  r.f_qi = r.index_qi > -1;
  r.f_qs = r.index_qs > -1;
  r.f_qg = r.index_qg > -1;
  r.index_nc = scalar_index(s, "nc");
  r.index_ni = scalar_index(s, "ni");
  r.index_nifa = scalar_index(s, "nifa");
  r.index_nwfa = scalar_index(s, "nwfa");
  r.index_nr = scalar_index(s, "nr");  // the state's index_nr (Thompson's microphysics_from_MPAS, PI:710-714)
  r.f_nc = r.index_nc > -1;
  r.f_ni = r.index_ni > -1;
  r.f_nifa = r.index_nifa > -1;
  r.f_nwfa = r.index_nwfa > -1;
  const int nCellsSolve = m.nCellsSolve();
  if (!c.config_do_restart && c.config_convection_scheme != "off") {  // IN:189-193
    auto& i_rainc = p.I("i_rainc");
    for (int iCell = 1; iCell <= nCellsSolve; ++iCell) i_rainc(iCell) = 0;
  }
  if (!c.config_do_restart && c.config_microp_scheme != "off") {  // IN:194-198
    auto& i_rainnc = p.I("i_rainnc");
    for (int iCell = 1; iCell <= nCellsSolve; ++iCell) i_rainnc(iCell) = 0;
  }
  {  // IN:202-212 mp_top_level
    r.mp_top_level = 0;
    real layer_height = 0.5 / m.rdzw(1);
    if (layer_height <= c.config_microphysics_top) r.mp_top_level = 1;
    for (int k = 2; k <= m.nVertLevels; ++k) {
      layer_height = layer_height + m.dzu(k);
      if (layer_height <= c.config_microphysics_top) r.mp_top_level = k;
    }
    log.push_back("physics init: config_microphysics_top and level " + std::to_string(c.config_microphysics_top) +
                  " " + std::to_string(r.mp_top_level));
  }
  if (!c.config_do_restart) {  // IN:216-251
    for (const char* n : {"i_acswdnb", "i_acswdnbc", "i_acswdnt", "i_acswdntc", "i_acswupb", "i_acswupbc", "i_acswupt",
                          "i_acswuptc", "i_aclwdnb", "i_aclwdnbc", "i_aclwdnt", "i_aclwdntc", "i_aclwupb",
                          "i_aclwupbc", "i_aclwupt", "i_aclwuptc"}) {
      auto& a = p.I(n);
      for (int iCell = 1; iCell <= nCellsSolve; ++iCell) a(iCell) = 0;
    }
    for (const char* n : {"acswdnb", "acswdnbc", "acswdnt", "acswdntc", "acswupb", "acswupbc", "acswupt", "acswuptc",
                          "aclwdnb", "aclwdnbc", "aclwdnt", "aclwdntc", "aclwupb", "aclwupbc", "aclwupt",
                          "aclwuptc"}) {
      auto& a = p.R(n);
      for (int iCell = 1; iCell <= nCellsSolve; ++iCell) a(iCell) = 0.;
    }
  }
  if (!c.config_do_restart) {  // IN:253-262
    auto &xicem = p.R("xicem"), &xice = p.R("xice"), &sstsk = p.R("sstsk"), &sst = p.R("sst");
    for (int iCell = 1; iCell <= nCellsSolve; ++iCell) xicem(iCell) = xice(iCell);
    for (int iCell = 1; iCell <= nCellsSolve; ++iCell) sstsk(iCell) = sst(iCell);
  }
  if (!c.config_do_restart) {  // IN:301-314
    auto &nsteps_accum = p.R("nsteps_accum"), &ndays_accum = p.R("ndays_accum"), &tday_accum = p.R("tday_accum");
    auto &tyear_accum = p.R("tyear_accum"), &tyear_mean = p.R("tyear_mean"), &tmn = p.R("tmn");
    auto& tlag = p.R("tlag");
    const int nLags = p.dim("nLags", false);
    for (int iCell = 1; iCell <= nCellsSolve; ++iCell) {
      nsteps_accum(iCell) = 0.;
      ndays_accum(iCell) = 0.;
      tday_accum(iCell) = 0.;
      tyear_accum(iCell) = 0.;
      tyear_mean(iCell) = tmn(iCell);
      for (int iLag = 1; iLag <= nLags; ++iLag) tlag(iLag, iCell) = tmn(iCell);
    }
  }
  if (c.config_o3climatology) {  // IN:317
    RadiationPools q = radiation_pools(r, m, p);
    init_o3climatology(r.table_dir, q, m.nCells);
  }
  {  // IN:320
    LanduseConfig lc;
    lc.config_do_restart = c.config_do_restart;
    lc.config_frac_seaice = c.config_frac_seaice;
    lc.config_sfc_albedo = c.config_sfc_albedo;
    lc.xice_threshold = r.xice_threshold;
    r.isurban = landuse_init_forMPAS(r.table_dir, julday, lc, m, p, log);
  }
  // IN:323-324 init_convection: nothing for cu_ntiedtke and cu_grell_freitas (CV:387-400: Kain-Fritsch only)
  // IN:327-328 init_microphysics (DRV:260-314)
  if (c.config_microp_scheme == "mp_wsm6") {  // mp_wsm6_init MW6:174-184 (denr = rho_r, dens = rho_s, DRV:305)
    r.wsm6 = module_mp_wsm6::mpas_wsm6_params();
    r.wsm6_radar = wsm6::wsm6_radar_init(r.wsm6.pi, physics_constants::rho_r, physics_constants::rho_s, r.wsm6.deng,
                                         r.wsm6.n0g);
  } else if (c.config_microp_scheme == "mp_thompson") {  // DRV:294-296
    if (r.thompson_dir.empty()) throw Error("physics_init: mp_thompson needs the directory of its MP_THOMPSON_*_DATA.DBL tables");
    thompson::thompson_init(r.thompson_dir, r.thompson);  // l_mp_tables = .true.: the tables are read
    init_thompson_clouddroplets_forMPAS(nCellsSolve, cv(p.I("landmask")), p.v("nt_c"), p.v("mu_c"));
  } else if (c.config_microp_scheme != "off") {
    throw Error("physics_init: config_microp_scheme = '" + c.config_microp_scheme + "' has no verified port");
  }
  // IN:331 init_pbl (PB:727-754): nothing for bl_ysu; bl_mynn_init for bl_mynn
  if (c.config_pbl_scheme == "bl_mynn") r.mynn_const = mynn_init_pbl();
  if (c.config_sfclayer_scheme != "off") init_sfclayer(c.config_sfclayer_scheme, r.sfclayrev_tables);  // IN:334
  if (c.config_lsm_scheme == "sf_noah") {  // IN:337-345
    log.push_back("--- initialize NOAH LSM tables");
    LsmConfig lc = lsm_config(r);
    LsmPools lp = lsm_pools(m, s, p);
    init_lsm(lc, lp, r.table_dir, r.noah_tables);
    log.push_back("--- end initialize NOAH LSM tables");
  }
  // IN:353-380 init_radiation_sw / init_radiation_lw: rrtmg_initsw_forMPAS / rrtmg_initlw_forMPAS
  if (c.config_radt_sw_scheme != "off" || c.config_radt_lw_scheme != "off")
    init_radiation_rrtmg(r.rrtmg_lw_tables, r.rrtmg_sw_tables, r.rrtmg);
}

// ---------------------------------------------------------------------------------------------------------------
void physics_timetracker(PhysicsRun& r, const Mesh& m, RegistryPools& p, const std::string& now_stamp,
                         long long now_s, real dt, std::vector<std::string>& log) {
  const PhysicsConfig& c = r.cfg;
  // MG:199-206
  physics_radiation_time_update(r.rad, num(now_stamp, 0, 4), day_of_year(now_stamp), num(now_stamp, 11, 2),
                                num(now_stamp, 14, 2), num(now_stamp, 17, 2), 0, 1);
  if (r.green_alarm.ringing(now_s)) {  // MG:216-220
    r.green_alarm.reset(now_s);
    log.push_back("--- time to update background surface albedo, greeness fraction.");
    physics_update_surface(now_stamp, c.config_sfc_albedo, m.nCellsSolve(), p);
  }
  // MG:222-235: config_sst_update, config_sstdiurn_update, config_deepsoiltemp_update are false
  auto radt = [&](const std::string& scheme, const std::string& interval, ClockAlarm& al, bool& l, const char* what) {
    if (scheme == "off") return;  // MG:241-270
    l = false;
    if (interval != "none") {
      if (al.ringing(now_s)) {
        al.reset(now_s);
        l = true;
      }
    } else {
      l = true;
    }
    log.push_back(std::string("--- time to run the ") + what + (l ? "T" : "F"));
  };
  radt(c.config_radt_lw_scheme, c.config_radtlw_interval, r.radtlw_alarm, r.l_radtlw, "LW radiation scheme L_RADLW = ");
  radt(c.config_radt_sw_scheme, c.config_radtsw_interval, r.radtsw_alarm, r.l_radtsw, "SW radiation scheme L_RADSW = ");
  if (c.config_convection_scheme != "off") {  // MG:272-283 (config_conv_interval = 'none')
    r.l_conv = true;
    log.push_back("--- time to run the convection scheme L_CONV    = T");
  }
  if (c.config_o3climatology) {  // MG:286-296
    if ((r.l_radtlw && c.config_radt_lw_scheme == "rrtmg_lw") || (r.l_radtsw && c.config_radt_sw_scheme == "rrtmg_sw")) {
      RadiationPools q = radiation_pools(r, m, p);
      o3climatology_from_MPAS(r.rad.curr_julday, q, m.nCellsSolve());
      log.push_back("--- time to update the ozone climatology for RRTMG radiation codes");
    }
  }
  // MG:320-345: the bucket alarms are not created (config_bucket_update = 'none'): never ringing
  if (c.config_convection_scheme != "off") {
    r.l_acrain = false;
    log.push_back("--- time to apply limit to accumulated rainc and rainnc L_ACRAIN   = F");
  }
  if (c.config_radt_lw_scheme != "off" || c.config_radt_sw_scheme != "off") {
    r.l_acradt = false;
    log.push_back("--- time to apply limit to accumulated radiation diags. L_ACRADT   = F");
  }
  // MG:356-371
  r.l_diags = false;
  if (r.diag_alarm.ringing(now_s)) r.diag_alarm.reset(now_s);
  if (r.diag_alarm.defined && now_s + static_cast<long long>(dt) >= r.diag_alarm.ring) r.l_diags = true;
  log.push_back(std::string("--- time to calculate additional physics_diagnostics               = ") +
                (r.l_diags ? "T" : "F"));
}

// ---------------------------------------------------------------------------------------------------------------
void physics_driver(PhysicsRun& r, const Mesh& m, State& s, Diag& d, RegistryPools& p, int itimestep,
                    real xtime_s, const DycorePhysicsArrays& a) {
  const PhysicsConfig& c = r.cfg;
  if (c.config_convection_scheme == "off" && c.config_lsm_scheme == "off" && c.config_pbl_scheme == "off" &&
      c.config_radt_lw_scheme == "off" && c.config_radt_sw_scheme == "off" && c.config_sfclayer_scheme == "off")
    return;  // DR:170-176
  const int its = 1, ite = m.nCellsSolve(), nV = m.nVertLevels;
  // DR:200 allocate_forall_physics, DR:203-207 MPAS_to_physics (time_lev = 1)
  PhysicsVars v;
  v.allocate_forall_physics(c.config_microp_scheme, c.config_pbl_scheme, ite, nV + 1);
  {
    ToPhysicsPools tp;
    tp.nVertLevels = nV;
    tp.microp_scheme = c.config_microp_scheme;
    tp.pbl_scheme = c.config_pbl_scheme;
    tp.zgrid = cv(m.zgrid);
    tp.zz = cv(m.zz);
    tp.rho_zz = cv(s.rho_zz);
    tp.theta_m = cv(s.theta_m);
    tp.w = cv(s.w);
    tp.scalars = cv(s.scalars);
    tp.index_qv = r.index_qv;
    tp.index_qc = r.index_qc;
    tp.index_qr = r.index_qr;
    tp.index_qi = r.index_qi;
    tp.index_qs = r.index_qs;
    tp.index_qg = r.index_qg;
    tp.index_ni = r.index_ni;  // bl_mynn (PI:390-397)
    tp.f_ni = r.f_ni;
    tp.f_nc = r.f_nc;
    tp.f_nifa = r.f_nifa;
    tp.f_nwfa = r.f_nwfa;
    tp.surface_pressure = s.surface_pressure.view();
    tp.exner = cv(d.exner);
    tp.pressure_base = cv(d.pressure_base);
    tp.pressure_p = cv(d.pressure_p);
    tp.uReconstructZonal = cv(d.uReconstructZonal);
    tp.uReconstructMeridional = cv(d.uReconstructMeridional);
    tp.plrad = p.v("plrad");
    MPAS_to_physics(tp, v, its, ite);
  }
  if (r.l_radtlw || r.l_radtsw) {  // DR:210-246
    // in the CUDA build the cloudiness and the RRTMG drivers run as kernels on the arrays MPAS_to_physics wrote and
    // return synchronised (GPU_RADIATION_INTERFACE: physics/mpas_atmphys_radiation.hpp)
    RadiationPools q = radiation_pools(r, m, p);
    r.rad.sw_ran = false;
    driver_cloudiness(q, v, r.rad, its, ite);
    if (r.l_radtsw) driver_radiation_sw(itimestep, q, v, r.rad, r.rrtmg, xtime_s, its, ite);
    if (r.l_radtlw) driver_radiation_lw(xtime_s, q, v, r.rad, r.rrtmg, its, ite);
  }
  // DR:229-233: config_bucket_update = 'none'
  SfclayerSeaArrays sea;
  if (c.config_sfclayer_scheme != "off") {  // DR:249-256
    SfclayerPools q;
    q.sfclayer_scheme = c.config_sfclayer_scheme;
    q.config_frac_seaice = c.config_frac_seaice;
    q.config_do_restart = c.config_do_restart;
    q.config_len_disp = c.config_len_disp;
    q.meshDensity = cv(m.meshDensity);
    q.areaCell = cv(m.areaCell);
    q.skintemp = cv(p.R("skintemp"));
    q.xland = cv(p.R("xland"));
    q.sst = cv(p.R("sst"));
    q.xice = cv(p.R("xice"));
    q.hpbl = cv(p.R("hpbl"));
    q.mavail = cv(p.R("mavail"));
#define CPPMPAS_P(x) q.x = p.v(#x);
    CPPMPAS_P(br) CPPMPAS_P(cpm) CPPMPAS_P(chs) CPPMPAS_P(chs2) CPPMPAS_P(cqs2) CPPMPAS_P(flhc) CPPMPAS_P(flqc)
    CPPMPAS_P(gz1oz0) CPPMPAS_P(hfx) CPPMPAS_P(qfx) CPPMPAS_P(qgh) CPPMPAS_P(qsfc) CPPMPAS_P(lh) CPPMPAS_P(mol)
    CPPMPAS_P(psim) CPPMPAS_P(psih) CPPMPAS_P(regime) CPPMPAS_P(rmol) CPPMPAS_P(ust) CPPMPAS_P(ustm) CPPMPAS_P(wspd)
    CPPMPAS_P(znt) CPPMPAS_P(zol) CPPMPAS_P(q2) CPPMPAS_P(t2m) CPPMPAS_P(th2m) CPPMPAS_P(u10) CPPMPAS_P(v10)
    CPPMPAS_P(cd) CPPMPAS_P(cda) CPPMPAS_P(ck) CPPMPAS_P(cka) CPPMPAS_P(fh) CPPMPAS_P(fm)
#undef CPPMPAS_P
    if (c.config_sfclayer_scheme == "sf_mynn") {  // SL:568-587, 773-790
      q.qcg = cv(p.R("qcg"));
      q.snowh = cv(p.R("snowh"));
      q.ch = p.v("ch");
    }
    q.xice_threshold = r.xice_threshold;
    driver_sfclayer(itimestep, r.sfclayrev_tables, q, v, its, ite, &sea);
  }
  // DR:259: config_oml1d = false
  if (c.config_lsm_scheme != "off") {  // DR:262-294
    LsmConfig lc = lsm_config(r);
    LsmPools lp = lsm_pools(m, s, p);
    LsmVars lv;
    lv.ims = 1;
    lv.ime = ite;
    lv.kms = 1;
    lv.kme = nV + 1;
    lv.num_soils = p.dim("nSoilLevels", false);
    lv.isurban = r.isurban;
    lv.dt_pbl = r.dt_pbl;
    lv.xice_threshold = r.xice_threshold;
    lv.dz_p = pv(v.dz_p, nV + 1);
    lv.pres2_hyd_p = pv(v.pres2_hyd_p, nV + 1);
    lv.t_p = pv(v.t_p, nV + 1);
    lv.qv_p = pv(v.qv_p, nV + 1);
    if (c.config_frac_seaice) {
      lv.cpm_sea = sv(sea.cpm);
      lv.chs_sea = sv(sea.chs);
      lv.chs2_sea = sv(sea.chs2);
      lv.cqs2_sea = sv(sea.cqs2);
      lv.hfx_sea = sv(sea.hfx);
      lv.qfx_sea = sv(sea.qfx);
      lv.qgh_sea = sv(sea.qgh);
      lv.qsfc_sea = sv(sea.qsfc);
      lv.lh_sea = sv(sea.lh);
    }
    if (c.config_lsm_scheme != "sf_noah") throw Error("physics_driver: config_lsm_scheme = '" + c.config_lsm_scheme + "'");
    driver_lsm(itimestep, lc, r.noah_tables, lp, lv, its, ite);
    driver_seaice(lc, lp, lv, its, ite);
    atmphys_sfc_diagnostics(lc, lp, its, ite);
  }
  if (c.config_pbl_scheme != "off" && c.config_sfclayer_scheme != "off") {  // DR:297-305
    PblPools q;
    q.nVertLevels = nV;
    q.pbl_scheme = c.config_pbl_scheme;
    q.config_ysu_pblmix = c.config_ysu_pblmix;
    q.xland = cv(p.R("xland"));
#define CPPMPAS_C(x) q.x = cv(p.R(#x));
#define CPPMPAS_P(x) q.x = p.v(#x);
    CPPMPAS_C(hfx) CPPMPAS_C(qfx) CPPMPAS_C(ust) CPPMPAS_C(wspd) CPPMPAS_C(znt) CPPMPAS_C(br) CPPMPAS_C(fm)
    CPPMPAS_C(fh) CPPMPAS_C(u10) CPPMPAS_C(v10)
    CPPMPAS_P(hpbl) CPPMPAS_P(kzh) CPPMPAS_P(kzm) CPPMPAS_P(kzq)
    if (c.config_pbl_scheme == "bl_ysu") {  // delta, wstar, exch_h: bl_ysu_in (Registry.xml), YSU's (PB:316-338)
      CPPMPAS_P(delta) CPPMPAS_P(wstar) CPPMPAS_P(exch_h)
    }
    CPPMPAS_C(rthratenlw) CPPMPAS_C(rthratensw)
    CPPMPAS_P(rublten) CPPMPAS_P(rvblten) CPPMPAS_P(rthblten) CPPMPAS_P(rqvblten) CPPMPAS_P(rqcblten)
    CPPMPAS_P(rqiblten)
#undef CPPMPAS_C
#undef CPPMPAS_P
    q.kpbl = p.vi("kpbl");
    q.dt_pbl = r.dt_pbl;
    q.f_qc = r.f_qc;
    q.f_qi = r.f_qi;
    if (c.config_pbl_scheme == "bl_mynn") {  // PB:316-338, 423-506, 551-562, 627-706
      q.config_do_restart = c.config_do_restart;
      q.config_len_disp = c.config_len_disp;
      q.mynn = mynn_options(c);
      q.mynn_const = r.mynn_const;
      q.f_qs = r.f_qs;
      q.f_qoz = false;  // IN:505: qoz is not in the Registry
      q.f_nc = r.f_nc;
      q.f_ni = r.f_ni;
      q.f_nifa = r.f_nifa;
      q.f_nwfa = r.f_nwfa;
      q.f_nbca = false;  // IN:524: nbca is not in the Registry
      q.meshDensity = cv(m.meshDensity);
#define CPPMPAS_C(x) q.x = cv(p.R(#x));
#define CPPMPAS_P(x) q.x = p.v(#x);
      CPPMPAS_C(skintemp) CPPMPAS_C(ch) CPPMPAS_C(qsfc) CPPMPAS_C(rmol)
      CPPMPAS_P(el_pbl) CPPMPAS_P(cov) CPPMPAS_P(qke) CPPMPAS_P(qke_adv) CPPMPAS_P(qsq) CPPMPAS_P(tsq) CPPMPAS_P(tke_pbl)
      CPPMPAS_P(sh3d) CPPMPAS_P(sm3d) CPPMPAS_P(cldfrac_bl) CPPMPAS_P(qc_bl) CPPMPAS_P(qi_bl) CPPMPAS_P(edmf_a)
      CPPMPAS_P(edmf_ent) CPPMPAS_P(edmf_qc) CPPMPAS_P(edmf_qt) CPPMPAS_P(edmf_thl) CPPMPAS_P(edmf_w) CPPMPAS_P(sub_thl)
      CPPMPAS_P(sub_qv) CPPMPAS_P(det_thl) CPPMPAS_P(det_qv) CPPMPAS_P(dqke) CPPMPAS_P(qbuoy) CPPMPAS_P(qdiss)
      CPPMPAS_P(qshear) CPPMPAS_P(qwt) CPPMPAS_P(rqsblten) CPPMPAS_P(rniblten)
#undef CPPMPAS_C
#undef CPPMPAS_P
    }
    driver_pbl(itimestep, q, v, its, ite);
  }
  if (c.config_gwdo_scheme != "off") {  // DR:308-316
    GwdoPools q;
    q.nVertLevels = nV;
    q.gwdo_scheme = c.config_gwdo_scheme;
    q.config_ngw_scheme = c.config_ngw_scheme;
    q.config_len_disp = c.config_len_disp;
    q.meshDensity = cv(m.meshDensity);
#define CPPMPAS_C(x) q.x = cv(p.R(#x));
#define CPPMPAS_P(x) q.x = p.v(#x);
    CPPMPAS_C(var2d) CPPMPAS_C(con) CPPMPAS_C(oa1) CPPMPAS_C(oa2) CPPMPAS_C(oa3) CPPMPAS_C(oa4) CPPMPAS_C(ol1)
    CPPMPAS_C(ol2) CPPMPAS_C(ol3) CPPMPAS_C(ol4)
    CPPMPAS_P(dusfcg) CPPMPAS_P(dvsfcg) CPPMPAS_P(dtaux3d) CPPMPAS_P(dtauy3d) CPPMPAS_P(rubldiff) CPPMPAS_P(rvbldiff)
    CPPMPAS_P(rublten) CPPMPAS_P(rvblten) CPPMPAS_P(rthblten)
#undef CPPMPAS_C
#undef CPPMPAS_P
    q.kpbl = cv(p.I("kpbl"));
    q.dt_pbl = r.dt_pbl;
    driver_gwdo(itimestep, q, v, its, ite);
  }
  if (c.config_convection_scheme != "off") {  // DR:319-341
    ConvectionPools q = convection_pools(r, m, p, a);
    update_convection_step1(q, its, ite);
    if (r.l_conv) driver_convection(itimestep, q, v, its, ite);
    update_convection_step2(q, its, ite);
  }
}

// ---------------------------------------------------------------------------------------------------------------
void physics_get_tend_for_dynamics(PhysicsRun& r, const Mesh& m, const State& s, Diag& d, RegistryPools& p,
                                   const FArray<real>& rho_zz2, const DycorePhysicsArrays& a) {
  const PhysicsConfig& c = r.cfg;
  GetTendPools q;
  q.nCells = m.nCells;
  q.nEdges = m.nEdges;
  q.nCellsSolve = m.nCellsSolve();
  q.nEdgesSolve = m.nEdgesSolve();
  q.nVertLevels = m.nVertLevels;
  q.pbl_scheme = c.config_pbl_scheme;
  q.convection_scheme = c.config_convection_scheme;
  q.microp_scheme = c.config_microp_scheme;
  q.radt_lw_scheme = c.config_radt_lw_scheme;
  q.radt_sw_scheme = c.config_radt_sw_scheme;
  q.cellsOnEdge = cv(m.cellsOnEdge);
  q.edgeNormalVectors = cv(m.edgeNormalVectors);
  q.east = cv(p.R("east"));
  q.north = cv(p.R("north"));
  q.theta_m = cv(s.theta_m);  // TD:163-166: theta_m, scalars of time level 1, rho_zz of time level 2
  q.scalars = cv(s.scalars);
  q.mass = cv(rho_zz2);
  q.mass_edge = cv(d.rho_edge);
  q.index_qv = r.index_qv;
  q.index_qc = r.index_qc;
  q.index_qr = r.index_qr;
  q.index_qi = r.index_qi;
  q.index_qs = r.index_qs;
  q.index_ni = r.index_ni;
  q.tend_u_phys = p.v("tend_u_phys");
  if (c.config_pbl_scheme != "off") {
#define CPPMPAS_C(x) q.x = cv(p.R(#x));
    CPPMPAS_C(rublten) CPPMPAS_C(rvblten) CPPMPAS_C(rthblten) CPPMPAS_C(rqvblten) CPPMPAS_C(rqcblten) CPPMPAS_C(rqiblten)
    q.rublten_Edge = p.v("rublten_Edge");
    if (c.config_pbl_scheme == "bl_mynn") {  // TD:405-412
      CPPMPAS_C(rqsblten) CPPMPAS_C(rniblten)
    }
  }
  if (c.config_convection_scheme != "off") {
    CPPMPAS_C(rucuten) CPPMPAS_C(rvcuten) CPPMPAS_C(rthcuten) CPPMPAS_C(rqvcuten) CPPMPAS_C(rqccuten) CPPMPAS_C(rqicuten)
    q.rucuten_Edge = p.v("rucuten_Edge");
  }
  if (c.config_radt_lw_scheme != "off") CPPMPAS_C(rthratenlw)
  if (c.config_radt_sw_scheme != "off") CPPMPAS_C(rthratensw)
#undef CPPMPAS_C
  q.tend_scalars = a.tend_scalars->view();
  physics_get_tend(q, 1, 1, a.tend_ru_physics->view(), a.tend_rtheta_physics->view(), a.tend_rho_physics->view(),
                   [](const char*) {});  // one rank: the halo exchanges are no-ops
}

// ---------------------------------------------------------------------------------------------------------------
void dynamics_tendencies_for_convection(PhysicsRun& r, const Mesh& m, const State& s, RegistryPools& p,
                                        const FArray<real>& scalars2, const FArray<real>& theta_m2, real dt,
                                        bool config_monotonic, const DycorePhysicsArrays& a) {
  const PhysicsConfig& c = r.cfg;
  if (c.config_convection_scheme != "cu_grell_freitas" && c.config_convection_scheme != "cu_ntiedtke") return;
  using physics_constants::rgas;
  using physics_constants::rv;
  const auto rqvdynten = p.v("rqvdynten");
  const auto rthdynten = a.rthdynten->view();
  const auto scalars2_v = scalars2.cview(), scalars1_v = s.scalars.cview(), theta_m2_v = theta_m2.cview();
  const int qv = r.index_qv, nV = m.nVertLevels;
  // TI:2769-2773: whole arrays (the garbage column included); one thread per cell, the k loop inside
  parallel_for(1, m.nCells + 1, [=] CPPMPAS_LAMBDA(int i) {
    const auto& scalars2 = scalars2_v;
    for (int k = 1; k <= nV; ++k)
      rqvdynten(k, i) = config_monotonic ? (scalars2(qv, k, i) - scalars1_v(qv, k, i)) / dt : 0.;
  });
  parallel_for(1, m.nCellsSolve(), [=] CPPMPAS_LAMBDA(int iCell) {  // TI:2777-2783
    const auto& scalars2 = scalars2_v;
    const auto& theta_m2 = theta_m2_v;
    for (int k = 1; k <= nV; ++k) {
      const real fac_m = 1. / (1. + rv / rgas * scalars2(qv, k, iCell));
      const real theta_local = theta_m2(k, iCell) * fac_m;
      rthdynten(k, iCell) = fac_m * (rthdynten(k, iCell) - theta_local * rv / rgas * rqvdynten(k, iCell));
    }
  });
}

// ---------------------------------------------------------------------------------------------------------------
void microphysics_for_dynamics(PhysicsRun& r, const Mesh& m, State& s, Diag& d, RegistryPools& p,
                               FArray<real>& theta_m2, FArray<real>& rho_zz2, FArray<real>& scalars2,
                               const FArray<real>& w2, int itimestep, const DycorePhysicsArrays& a) {
  (void)itimestep;
  const PhysicsConfig& c = r.cfg;
  if (c.config_microp_scheme == "off") return;
  if (c.config_microp_scheme != "mp_wsm6" && c.config_microp_scheme != "mp_thompson")
    throw Error("microphysics: config_microp_scheme = '" + c.config_microp_scheme + "' with the physics suite");
  const bool thompson = c.config_microp_scheme == "mp_thompson";
  MicrophysicsPools q;
  q.nVertLevels = m.nVertLevels;
  q.zgrid = cv(m.zgrid);
  q.zz = cv(m.zz);
  q.rho_zz = rho_zz2.view();
  q.theta_m = theta_m2.view();
  q.scalars = scalars2.view();
  q.w = cv(w2);
  q.index_qv = r.index_qv;
  q.index_qc = r.index_qc;
  q.index_qr = r.index_qr;
  q.index_qi = r.index_qi;
  q.index_qs = r.index_qs;
  q.index_qg = r.index_qg;
  q.exner = d.exner.view();
  q.pressure_p = d.pressure_p.view();
  q.rtheta_p = d.rtheta_p.view();
  q.surface_pressure = s.surface_pressure.view();
  q.pressure_base = cv(d.pressure_base);
  q.exner_base = cv(d.exner_base);
  q.rtheta_base = cv(d.rtheta_base);
  q.rt_diabatic_tend = a.rt_diabatic_tend->view();
  q.tend_sfc_pressure = p.v("tend_sfc_pressure");
#define CPPMPAS_P(x) q.x = p.v(#x);
  CPPMPAS_P(dtheta_dt_mp) CPPMPAS_P(precipw) CPPMPAS_P(rainnc) CPPMPAS_P(rainncv) CPPMPAS_P(snownc)
  CPPMPAS_P(snowncv) CPPMPAS_P(graupelnc) CPPMPAS_P(graupelncv) CPPMPAS_P(sr) CPPMPAS_P(rainprod) CPPMPAS_P(evapprod)
  CPPMPAS_P(re_cloud) CPPMPAS_P(re_ice) CPPMPAS_P(re_snow) CPPMPAS_P(rthmpten) CPPMPAS_P(rqvmpten) CPPMPAS_P(rqcmpten)
  CPPMPAS_P(rqrmpten) CPPMPAS_P(rqimpten) CPPMPAS_P(rqsmpten) CPPMPAS_P(rqgmpten)
#undef CPPMPAS_P
  q.mp_top_level = r.mp_top_level;
  q.dt_dyn = r.dt_dyn;
  q.dt_microp = r.dt_microp;
  q.has_reqc = r.has_reqc;
  q.has_reqi = r.has_reqi;
  q.has_reqs = r.has_reqs;
  q.l_diags = r.l_diags;
  if (r.l_diags) {  // DRV:467-482: compute_radar_reflectivity (diag_physics) and compute_relhum (diag)
    q.relhum = p.v("relhum");
    q.refl10cm = p.v("refl10cm");
    q.refl10cm_max = p.v("refl10cm_max");
    q.refl10cm_1km = p.v("refl10cm_1km");
    q.refl10cm_1km_max = p.v("refl10cm_1km_max");
    if (!thompson) q.radar = &r.wsm6_radar;
  }
  q.l_acrain = r.l_acrain;
  if (thompson) {  // DRV:378-398: ni, nr, nt_c, mu_c, rnimpten, rnrmpten and n_microp sub-steps
    q.index_ni = r.index_ni;
    q.index_nr = r.index_nr;
    q.nt_c = cv(p.R("nt_c"));
    q.mu_c = cv(p.R("mu_c"));
    q.rnimpten = p.v("rnimpten");
    q.rnrmpten = p.v("rnrmpten");
    q.n_microp = r.n_microp;
    driver_microphysics_thompson(r.thompson.v, q, 1, m.nCellsSolve());
  } else {
    driver_microphysics_wsm6(r.wsm6, q, 1, m.nCellsSolve());
  }
}

// ---------------------------------------------------------------------------------------------------------------
// CT:459-533 physics_compatibility_check
void physics_compatibility_check(const PhysicsConfig& c, const Mesh& m, RegistryPools& p) {
  if (c.config_gwdo_scheme != "bl_ysu_gwdo") return;
  const auto& var2d = p.R("var2d");
  real maxvar2d = -std::numeric_limits<real>::max();
  for (int iCell = 1; iCell <= m.nCellsSolve(); ++iCell) maxvar2d = std::max(maxvar2d, var2d(iCell));
  const int iswater = p.I("iswater_lu")(1);
  const auto& ivgtyp = p.I("ivgtyp");
  bool all_water = true;
  for (int iCell = 1; iCell <= m.nCellsSolve(); ++iCell)
    if (ivgtyp(iCell) != iswater) all_water = false;
  if (maxvar2d <= 0.0 && !all_water)
    throw Error("The YSU GWDO scheme requires valid var2d, con, oa{1,2,3,4}, and ol{1,2,3,4} fields, but these fields "
                "appear to be zero everywhere in the model input.");
}

}  // namespace cppmpas::phys
