// SPDX-License-Identifier: GPL-3.0-or-later
// See mpas_atmphys_radiation.hpp. CL = physics/mpas_atmphys_driver_cloudiness.F, RS = mpas_atmphys_driver_radiation_sw.F,
// RL = mpas_atmphys_driver_radiation_lw.F, OZ = mpas_atmphys_o3climatology.F, MG = mpas_atmphys_manager.F,
// CS = physics_wrf/module_ra_cam_support.F (MPAS v8.4.2).
#include "physics/mpas_atmphys_radiation.hpp"

#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>

#include "core/error.hpp"
#include "core/libm.hpp"
#include "physics/mpas_atmphys_constants.hpp"

namespace cppmpas::phys {

namespace {

// mpas_atmphys_constants.F:88-90
constexpr real solcon_0 = 1370.;
constexpr real degrad = 3.1415926 / 180.;
constexpr real dpd = 360. / 365.;

std::size_t I(int i) { return static_cast<std::size_t>(i - 1); }

// the cloud overlap and decorrelation length options of RL:850-861 / RS:904-915 (cldovrlp is undefined for any other
// value of config_radt_cld_overlap: refused)
void cloud_overlap(const RadiationPools& m, int& cldovrlp, int& idcor) {
  const std::string& o = m.radt_cld_overlap;
  if (o == "none") cldovrlp = 0;
  else if (o == "random") cldovrlp = 1;
  else if (o == "maximum_random") cldovrlp = 2;
  else if (o == "maximum") cldovrlp = 3;
  else if (o == "exponential") cldovrlp = 4;
  else if (o == "exponential_random") cldovrlp = 5;
  else throw Error("config_radt_cld_overlap = '" + o + "' is not an MPAS option");
  idcor = 0;
  if (o == "exponential" || o == "exponential_random") {
    if (m.radt_cld_dcorrlen == "constant") idcor = 0;
    else if (m.radt_cld_dcorrlen == "latitude_varying") idcor = 1;
    else throw Error("config_radt_cld_dcorrlen = '" + m.radt_cld_dcorrlen + "' is not an MPAS option");
  }
}

void check_microp(const RadiationPools& m) {
  if (m.config_microp_re)
    throw Error("config_microp_re = true (the microphysics' effective radii in RRTMG) has no verified port");
  if (m.microp_scheme == "mp_thompson_aerosols")
    throw Error("radiation with 'mp_thompson_aerosols' (aer_opt = 3, calc_aerosol_rrtmg_sw) has no verified port");
}

// list-directed READ of one real per record (gfortran's list-directed real input is correctly rounded, as strtod)
real read_real(std::ifstream& f, const std::string& name) {
  std::string line;
  while (std::getline(f, line)) {
    const char* s = line.c_str();
    char* end = nullptr;
    errno = 0;
    const double x = std::strtod(s, &end);
    if (end == s) continue;  // a blank record: list-directed input reads on
    return x;
  }
  throw Error("init_o3climatology: " + name + " ends early");
}

}  // namespace

void init_radiation_rrtmg(const std::string& lw_tables_bin, const std::string& sw_tables_bin, RrtmgTables& t) {
  if (!rrtmg::load_rrtmg_lw_tables(lw_tables_bin, t.lw))
    throw Error("init_radiation_lw: cannot read the RRTMG LW tables " + lw_tables_bin);
  if (!rrtmg::load_rrtmg_sw_tables(sw_tables_bin, t.sw))
    throw Error("init_radiation_sw: cannot read the RRTMG SW tables " + sw_tables_bin);
}

// MG:457-463
void physics_radiation_time_init(RadiationVars& r, int yr, int doy, int h, int m, int s, long s_n, long s_d) {
  if (s_d == 0) throw Error("physics_run_init: S_d = 0");
  const real utc_h = real(h) + real(m) / 60.0 + real(s + int(s_n / s_d)) / 3600.0;
  r.year = yr;
  r.gmt = utc_h;
  r.julday = doy;
  r.curr_julday = real(r.julday - 1);
}

// MG:205-212
void physics_radiation_time_update(RadiationVars& r, int yr, int doy, int h, int m, int s, long s_n, long s_d) {
  if (s_d == 0) throw Error("physics_timetracker: S_d = 0");
  const real utc_h = real(h) + real(m) / 60.0 + real(s + int(s_n / s_d)) / 3600.0;
  r.year = yr;
  r.julday = doy;
  r.curr_julday = real(r.julday - 1) + utc_h / 24.0;
}

// RS:1009-1050
void radconst(real& declin, real& solcon, real julian, real degrad_, real dpd_) {
  declin = 0.;
  solcon = 0.;
  const real obecl = 23.5 * degrad_;
  const real sinob = libm::sin(obecl);
  real sxlong = 0.;
  if (julian >= 80.) sxlong = dpd_ * (julian - 80.);
  if (julian < 80.) sxlong = dpd_ * (julian + 285.);
  sxlong = sxlong * degrad_;
  const real arg = sinob * libm::sin(sxlong);
  declin = libm::asin(arg);
  const real djul = julian * 360. / 365.;
  const real rjul = djul * degrad_;
  // the reference takes sincos(rjul) and sincos(2*rjul) (2*rjul: the integer 2 converted, an exact doubling)
  real s1, c1, s2, c2;
  libm::sincos(rjul, &s1, &c1);
  libm::sincos(2 * rjul, &s2, &c2);
  const real eccfac = 1.000110 + 0.034221 * c1 + 0.001280 * s1 + 0.000719 * c2 + 0.000077 * s2;
  solcon = solcon_0 * eccfac;
}

// OZ:62-180
void init_o3climatology(const std::string& dir, RadiationPools& m, int nCells) {
  constexpr int latsiz = 64, lonsiz = 1;  // OZ:23-24
  const int levsiz = m.nOznLevels, num_months = m.nMonths;
  auto open = [&](const std::string& name) {
    std::ifstream f(dir + "/" + name);
    if (!f) throw Error("subroutine oznini: failure opening " + name);
    return f;
  };
  {
    std::ifstream f = open("OZONE_PLEV.TBL");
    for (int k = 1; k <= levsiz; ++k) m.pin(k) = read_real(f, "OZONE_PLEV.TBL");  // OZ:107-109
  }
  real lat_ozone[latsiz];
  {
    std::ifstream f = open("OZONE_LAT.TBL");
    for (int j = 1; j <= latsiz; ++j) lat_ozone[j - 1] = read_real(f, "OZONE_LAT.TBL");  // OZ:119-122
  }
  // ozmixin(lonsiz, levsiz, latsiz, num_months) (OZ:134-143)
  std::vector<real> ozmixin(static_cast<std::size_t>(lonsiz) * levsiz * latsiz * num_months);
  auto OZI = [&](int i, int k, int j, int mm) -> real& {
    return ozmixin[I(i) + std::size_t(lonsiz) * (I(k) + std::size_t(levsiz) * (I(j) + std::size_t(latsiz) * I(mm)))];
  };
  {
    std::ifstream f = open("OZONE_DAT.TBL");
    for (int mm = 1; mm <= num_months; ++mm)
      for (int j = 1; j <= latsiz; ++j)
        for (int k = 1; k <= levsiz; ++k)
          for (int i = 1; i <= lonsiz; ++i) OZI(i, k, j, mm) = read_real(f, "OZONE_DAT.TBL");
  }
  for (int iCell = 1; iCell <= nCells; ++iCell) {  // OZ:154-184
    const real lat = m.latCell(iCell) / degrad;
    int i1, i2;
    if (lat > lat_ozone[latsiz - 1]) {
      i1 = latsiz;
      i2 = latsiz;
    } else if (lat < lat_ozone[0]) {
      i1 = 1;
      i2 = 1;
    } else {
      int i = 1;
      for (; i <= latsiz; ++i) {
        if (i == latsiz)  // lat_ozone(latsiz + 1) is past the array (lat == lat_ozone(latsiz))
          throw Error("init_o3climatology: a cell latitude equal to the last ozone latitude");
        if (lat >= lat_ozone[i - 1] && lat < lat_ozone[i]) break;
      }
      i1 = i;
      i2 = i + 1;
    }
    for (int mm = 1; mm <= num_months; ++mm)
      for (int k = 1; k <= levsiz; ++k)
        for (int j = 1; j <= lonsiz; ++j) {
          const real dlat = lat_ozone[i2 - 1] - lat_ozone[i1 - 1];
          const real dlatCell = lat - lat_ozone[i1 - 1];
          if (dlat == 0.) {
            m.ozmixm(mm, k, iCell) = OZI(j, k, i1, mm);
          } else {
            m.ozmixm(mm, k, iCell) = OZI(j, k, i1, mm) + (OZI(j, k, i2, mm) - OZI(j, k, i1, mm)) * dlatCell / dlat;
          }
        }
  }
}

// OZ:228-325 with getfactors (CS:2528-2611; r8 = double, ozncyc = .true.)
void o3climatology_from_MPAS(real julian, RadiationPools& m, int nCellsSolve) {
  static constexpr int date_oz[12] = {16, 45, 75, 105, 136, 166, 197, 228, 258, 289, 319, 350};
  const int nMonths = m.nMonths;
  if (nMonths != 12) throw Error("o3climatology_from_MPAS: nMonths /= 12");
  double intjulian_r8 = double(julian + 1.0);  // OZ:279
  int ijul = int(intjulian_r8);
  intjulian_r8 = intjulian_r8 - double(ijul);
  ijul = ijul % 365;
  if (ijul == 0) ijul = 365;
  intjulian_r8 = intjulian_r8 + double(ijul);
  int np1 = 1;  // OZ:289-296
  bool finddate = false;
  for (int mm = 1; mm <= nMonths; ++mm) {
    if (date_oz[mm - 1] > intjulian_r8 && !finddate) {
      np1 = mm;
      finddate = true;
    }
  }
  const double cdayozp_r8 = date_oz[np1 - 1];
  double cdayozm_r8;
  int np, nm;
  if (np1 > 1) {
    cdayozm_r8 = date_oz[np1 - 2];
    np = np1;
    nm = np - 1;
  } else {
    cdayozm_r8 = date_oz[11];
    np = np1;
    nm = nMonths;
  }
  // getfactors (CS:2576-2603): the sanity prints (CS:2575-2587, 2605-2608) write to unit 6 only
  constexpr double daysperyear = 365.;
  const double cday = intjulian_r8;
  double fact1_r8, fact2_r8;
  if (np1 == 1) {
    const double deltat = cdayozp_r8 + daysperyear - cdayozm_r8;
    if (cday > cdayozp_r8) {
      fact1_r8 = (cdayozp_r8 + daysperyear - cday) / deltat;
      fact2_r8 = (cday - cdayozm_r8) / deltat;
    } else {
      fact1_r8 = (cdayozp_r8 - cday) / deltat;
      fact2_r8 = (cday + daysperyear - cdayozm_r8) / deltat;
    }
  } else {
    const double deltat = cdayozp_r8 - cdayozm_r8;
    fact1_r8 = (cdayozp_r8 - cday) / deltat;
    fact2_r8 = (cday - cdayozm_r8) / deltat;
  }
  const real fact1 = real(fact1_r8), fact2 = real(fact2_r8);
  for (int k = 1; k <= m.nOznLevels; ++k)  // OZ:317-321
    for (int iCell = 1; iCell <= nCellsSolve; ++iCell)
      m.o3clim(k, iCell) = fact1 * m.ozmixm(nm, k, iCell) + fact2 * m.ozmixm(np, k, iCell);
}

// ---- cloudiness ---------------------------------------------------------------------------------------------------
void driver_cloudiness(RadiationPools& m, const PhysicsVars& v, RadiationVars& r, int its, int ite) {
  using physics_constants::ep_2;
  if (its != 1) throw Error("driver_cloudiness: its /= 1");
  const int kts = 1, kte = m.nVertLevels, n = ite, kme = kte + 1;
  if (m.radt_cld_scheme != "cld_fraction")
    throw Error("driver_cloudiness: config_radt_cld_scheme = '" + m.radt_cld_scheme + "' has no verified port");
  // allocate_cloudiness (CL:64-75)
  for (PArray* a : {&r.cldfrac_p, &r.qvrad_p, &r.qcrad_p, &r.qirad_p, &r.qsrad_p}) a->allocate(n, kme);
  r.xland_cl.assign(static_cast<std::size_t>(n), 0.);
  // cloudiness_from_MPAS (CL:92-143); dx_p (CL:124-126) left out (mpas_atmphys_radiation.hpp)
  for (int i = its; i <= ite; ++i) r.xland_cl[I(i)] = m.xland(i);
  for (int k = kts; k <= kte; ++k)
    for (int i = its; i <= ite; ++i) {
      r.qvrad_p(i, k) = v.qv_p(i, k);
      r.qcrad_p(i, k) = v.qc_p(i, k);
      r.qirad_p(i, k) = v.qi_p(i, k);
      r.qsrad_p(i, k) = v.qs_p(i, k);
      r.cldfrac_p(i, k) = 0.;
    }
  // calc_cldfraction (CL:284-370): cldfrac_p, t_p, pres_p, qvrad_p, qcrad_p, qirad_p, qs_p
  const real alpha0 = 100., gamma = 0.49, qcldmin = 1.e-12, pexp = 0.25, rhgrid = 1.0;
  const real svp1 = 0.61078, svp2 = 17.2693882, svpi2 = 21.8745584, svp3 = 35.86, svpi3 = 7.66, svpt0 = 273.15;
  for (int k = kts; k <= kte; ++k)
    for (int i = its; i <= ite; ++i) r.cldfrac_p(i, k) = 0.;
  for (int k = kts; k <= kte; ++k)
    for (int i = its; i <= ite; ++i) {
      const real t = v.t_p(i, k), p = v.pres_p(i, k);
      const real qv = r.qvrad_p(i, k), qc = r.qcrad_p(i, k), qi = r.qirad_p(i, k), qs = v.qs_p(i, k);
      const real esw = 1000. * svp1 * libm::exp(svp2 * (t - svpt0) / (t - svp3));   // CL:329
      const real esi = 1000. * svp1 * libm::exp(svpi2 * (t - svpt0) / (t - svpi3));  // CL:330
      const real qvsw = ep_2 * esw / (p - esw);
      const real qvsi = ep_2 * esi / (p - esi);
      const real qcld = qc + qi + qs;
      real weight;
      if (qcld < qcldmin) weight = 0.;
      else weight = (qi + qs) / qcld;
      const real qvs = (1 - weight) * qvsw + weight * qvsi;
      real rhum = qv / qvs;
      if (qcld < qcldmin) {
        r.cldfrac_p(i, k) = 0.;
      } else if (rhum >= rhgrid) {
        r.cldfrac_p(i, k) = 1.;
      } else {
        const real subsat = std::max(1.e-10, rhgrid * qvs - qv);
        const real denom = libm::pow(subsat, gamma);
        const real arg = std::max(-6.9, -alpha0 * qcld / denom);
        rhum = std::max(1.e-10, rhum);
        r.cldfrac_p(i, k) = libm::pow(rhum / rhgrid, pexp) * (1. - libm::exp(arg));  // CL:361
        if (r.cldfrac_p(i, k) < 0.01) r.cldfrac_p(i, k) = 0.;
      }
    }
  // cloudiness_to_MPAS (CL:146-172)
  for (int k = kts; k <= kte; ++k)
    for (int i = its; i <= ite; ++i) m.cldfrac(k, i) = r.cldfrac_p(i, k);
}

// ---- shortwave ----------------------------------------------------------------------------------------------------
void driver_radiation_sw(int itimestep, RadiationPools& m, const PhysicsVars& v, RadiationVars& r,
                         const RrtmgTables& t, real xtime_s, int its, int ite) {
  (void)itimestep;
  if (its != 1) throw Error("driver_radiation_sw: its /= 1");
  if (m.radt_sw_scheme != "rrtmg_sw")
    throw Error("driver_radiation_sw: config_radt_sw_scheme = '" + m.radt_sw_scheme + "' has no verified port");
  check_microp(m);
  if (r.cldfrac_p.ldi != ite) throw Error("driver_radiation_sw: driver_cloudiness has not run");
  const int kts = 1, kte = m.nVertLevels, n = ite, kme = kte + 1;
  const int nozn = m.nOznLevels;
  // allocate_radiation_sw (RS:106-235, the arrays rrtmg_sw uses)
  for (std::vector<real>* a : {&r.xlat_p, &r.xlon_p, &r.sfc_albedo_p, &r.snow_p, &r.tsk_p, &r.xice_p, &r.xland_p,
                               &r.coszr_p, &r.gsw_p, &r.swcf_p, &r.swdnb_p, &r.swdnbc_p, &r.swdnt_p, &r.swdntc_p,
                               &r.swupb_p, &r.swupbc_p, &r.swupt_p, &r.swuptc_p, &r.swddir_p, &r.swddni_p,
                               &r.swddif_p})
    a->assign(static_cast<std::size_t>(n), 0.);
  for (PArray* a : {&r.f_ice, &r.f_rain, &r.rthratensw_p, &r.recloud_p, &r.reice_p, &r.resnow_p}) a->allocate(n, kme);
  r.pin_p.assign(static_cast<std::size_t>(nozn), 0.);
  r.o3clim_p.allocate(n, nozn);

  const real xtime_m = xtime_s / 60.;  // RS:874
  // radiation_sw_from_MPAS (RS:355-711)
  for (int i = its; i <= ite; ++i) {  // RS:419-429
    r.xlat_p[I(i)] = m.latCell(i) / degrad;
    r.xlon_p[I(i)] = m.lonCell(i) / degrad;
    r.sfc_albedo_p[I(i)] = m.sfc_albedo(i);
    r.snow_p[I(i)] = m.snow(i);
    r.tsk_p[I(i)] = m.skintemp(i);
    r.xice_p[I(i)] = m.xice(i);
    r.xland_p[I(i)] = m.xland(i);
  }
  r.sw_ran = true;
  for (int k = kts; k <= kte; ++k)  // RS:430-436
    for (int i = its; i <= ite; ++i) r.cldfrac_p(i, k) = m.cldfrac(k, i);
  for (int k = kts; k <= kte; ++k)  // RS:439-445
    for (int i = its; i <= ite; ++i) {
      r.f_ice(i, k) = 0.;
      r.f_rain(i, k) = 0.;
    }
  for (int i = its; i <= ite; ++i) {  // RS:448-462 (assigned 0 above)
    r.coszr_p[I(i)] = r.gsw_p[I(i)] = r.swcf_p[I(i)] = 0.;
    r.swdnb_p[I(i)] = r.swdnbc_p[I(i)] = r.swdnt_p[I(i)] = r.swdntc_p[I(i)] = 0.;
    r.swupb_p[I(i)] = r.swupbc_p[I(i)] = r.swupt_p[I(i)] = r.swuptc_p[I(i)] = 0.;
  }
  for (int k = kts; k <= kte; ++k)
    for (int i = its; i <= ite; ++i) r.rthratensw_p(i, k) = 0.;
  // aer_opt = 0 (RS:471); swddir_p, swddni_p, swddif_p = 0 (RS:476-482, assigned above); tauaer_p = 0,
  // ssaaer_p = 1, asyaer_p = 0 (RS:484-492): the aerosol values of CPPWRF's spcvmc_sw (module_ra_rrtmg.hpp)
  if (m.microp_scheme == "mp_thompson" || m.microp_scheme == "mp_wsm6") {  // RS:494-531, config_microp_re false
    for (int k = kts; k <= kte; ++k)
      for (int i = its; i <= ite; ++i) {
        r.recloud_p(i, k) = 0.;
        r.reice_p(i, k) = 0.;
        r.resnow_p(i, k) = 0.;
      }
  }
  // swdnflx_p ... = 0 (RS:588-597): not passed to rrtmg_swrad, never copied to MPAS
  if (m.config_o3climatology) {  // RS:600-623
    for (int k = 1; k <= nozn; ++k) r.pin_p[I(k)] = m.pin(k);
    for (int k = 1; k <= nozn; ++k)
      for (int i = its; i <= ite; ++i) r.o3clim_p(i, k) = m.o3clim(k, i);
  } else {
    for (int k = 1; k <= nozn; ++k) r.pin_p[I(k)] = 0.0;
    for (int k = 1; k <= nozn; ++k)
      for (int i = its; i <= ite; ++i) r.o3clim_p(i, k) = 0.0;
  }

  radconst(r.declin, r.solcon, r.curr_julday, degrad, dpd);  // RS:884
  const real radt = r.dt_radtsw / 60.;                         // RS:896
  RrtmgRadIn in;
  in.kte = kte;
  in.o3input = 0;
  if (m.config_o3climatology) in.o3input = 2;  // RS:901-902
  cloud_overlap(m, in.cldovrlp, in.idcor);     // RS:904-915
  in.p3d = &v.pres_hyd_p; in.p8w = &v.pres2_hyd_p; in.pi3d = &v.pi_p; in.t3d = &v.t_p; in.t8w = &v.t2_p;
  in.dz8w = &v.dz_p; in.qv3d = &r.qvrad_p; in.qc3d = &r.qcrad_p; in.qi3d = &r.qirad_p; in.qs3d = &r.qsrad_p;
  in.cldfra3d = &r.cldfrac_p; in.tsk = &r.tsk_p; in.xland = &r.xland_p; in.xice = &r.xice_p; in.snow = &r.snow_p;
  in.xlat = &r.xlat_p; in.julday = r.julday; in.icloud = r.icloud; in.noznlevels = nozn; in.pin = r.pin_p.data();
  in.o3clim = &r.o3clim_p; in.has_reqc = r.has_reqc; in.has_reqi = r.has_reqi; in.has_reqs = r.has_reqs;
  in.re_cloud = &r.recloud_p; in.re_ice = &r.reice_p; in.re_snow = &r.resnow_p;
  RrtmgSwradOut out{&r.coszr_p, &r.gsw_p,    &r.swcf_p,   &r.swupt_p,  &r.swuptc_p, &r.swdnt_p,
                    &r.swdntc_p, &r.swupb_p, &r.swupbc_p, &r.swdnb_p,  &r.swdnbc_p, &r.swddir_p,
                    &r.swddni_p, &r.swddif_p, &r.rthratensw_p};
  rrtmg_swrad(in, r.sfc_albedo_p, r.xlon_p, xtime_m, r.gmt, radt, degrad, r.declin, r.solcon, out, t.sw, its, ite);

  // radiation_sw_to_MPAS (RS:714-794)
  for (int i = its; i <= ite; ++i) {
    m.coszr(i) = r.coszr_p[I(i)];
    m.gsw(i) = r.gsw_p[I(i)];
    m.swcf(i) = r.swcf_p[I(i)];
    m.swdnb(i) = r.swdnb_p[I(i)];
    m.swdnbc(i) = r.swdnbc_p[I(i)];
    m.swdnt(i) = r.swdnt_p[I(i)];
    m.swdntc(i) = r.swdntc_p[I(i)];
    m.swupb(i) = r.swupb_p[I(i)];
    m.swupbc(i) = r.swupbc_p[I(i)];
    m.swupt(i) = r.swupt_p[I(i)];
    m.swuptc(i) = r.swuptc_p[I(i)];
  }
  for (int i = its; i <= ite; ++i) {  // RS:774-780
    m.swddir(i) = r.swddir_p[I(i)];
    m.swddni(i) = r.swddni_p[I(i)];
    m.swddif(i) = r.swddif_p[I(i)];
  }
  for (int k = kts; k <= kte; ++k)
    for (int i = its; i <= ite; ++i) m.rthratensw(k, i) = r.rthratensw_p(i, k);
}

// ---- longwave -----------------------------------------------------------------------------------------------------
void driver_radiation_lw(real xtime_s, RadiationPools& m, const PhysicsVars& v, RadiationVars& r,
                         const RrtmgTables& t, int its, int ite) {
  (void)xtime_s;
  if (its != 1) throw Error("driver_radiation_lw: its /= 1");
  if (m.radt_lw_scheme != "rrtmg_lw")
    throw Error("driver_radiation_lw: config_radt_lw_scheme = '" + m.radt_lw_scheme + "' has no verified port");
  check_microp(m);
  if (r.cldfrac_p.ldi != ite) throw Error("driver_radiation_lw: driver_cloudiness has not run");
  const int kts = 1, kte = m.nVertLevels, n = ite, kme = kte + 1;
  const int nozn = m.nOznLevels;
  // allocate_radiation_lw (RL:103-206)
  for (std::vector<real>* a : {&r.sfc_emiss_p, &r.snow_p, &r.tsk_p, &r.xice_p, &r.xland_p, &r.glw_p, &r.lwcf_p,
                               &r.lwdnb_p, &r.lwdnbc_p, &r.lwdnt_p, &r.lwdntc_p, &r.lwupb_p, &r.lwupbc_p, &r.lwupt_p,
                               &r.lwuptc_p, &r.olrtoa_p})
    a->assign(static_cast<std::size_t>(n), 0.);
  for (PArray* a : {&r.f_ice, &r.f_rain, &r.rthratenlw_p, &r.recloud_p, &r.reice_p, &r.resnow_p, &r.rrecloud_p,
                    &r.rreice_p, &r.rresnow_p})
    a->allocate(n, kme);
  r.pin_p.assign(static_cast<std::size_t>(nozn), 0.);
  r.o3clim_p.allocate(n, nozn);

  // radiation_lw_from_MPAS (RL:297-610)
  for (int i = its; i <= ite; ++i) {  // RL:364-372
    r.sfc_emiss_p[I(i)] = m.sfc_emiss(i);
    r.tsk_p[I(i)] = m.skintemp(i);
    r.snow_p[I(i)] = m.snow(i);
    r.xice_p[I(i)] = m.xice(i);
    r.xland_p[I(i)] = m.xland(i);
  }
  for (int k = kts; k <= kte; ++k)  // RL:373-379
    for (int i = its; i <= ite; ++i) r.cldfrac_p(i, k) = m.cldfrac(k, i);
  for (int k = kts; k <= kte; ++k)  // RL:382-389
    for (int i = its; i <= ite; ++i) {
      r.f_ice(i, k) = 0.;
      r.f_rain(i, k) = 0.;
    }
  for (int i = its; i <= ite; ++i) {  // RL:391-404
    r.glw_p[I(i)] = r.lwcf_p[I(i)] = r.lwdnb_p[I(i)] = r.lwdnbc_p[I(i)] = r.lwdnt_p[I(i)] = 0.;
    r.lwdntc_p[I(i)] = r.lwupb_p[I(i)] = r.lwupbc_p[I(i)] = r.lwupt_p[I(i)] = r.lwuptc_p[I(i)] = 0.;
    r.olrtoa_p[I(i)] = 0.;
  }
  for (int k = kts; k <= kte; ++k)  // RL:405-409
    for (int i = its; i <= ite; ++i) r.rthratenlw_p(i, k) = 0.;
  const bool microp_re_scheme =
      m.microp_scheme == "mp_thompson" || m.microp_scheme == "mp_thompson_aerosols" || m.microp_scheme == "mp_wsm6";
  if (microp_re_scheme) {  // RL:416-458 (config_microp_re false)
    for (int k = kts; k <= kte; ++k)
      for (int i = its; i <= ite; ++i) {
        r.recloud_p(i, k) = 0.;
        r.reice_p(i, k) = 0.;
        r.resnow_p(i, k) = 0.;
      }
    for (int k = kts; k <= kte; ++k)
      for (int i = its; i <= ite; ++i) {
        r.rrecloud_p(i, k) = 0.;
        r.rreice_p(i, k) = 0.;
        r.rresnow_p(i, k) = 0.;
      }
  }
  // lwdnflx_p ... = 0 (RL:459-467): not passed to rrtmg_lwrad, never copied to MPAS
  if (m.config_o3climatology) {  // RL:470-516
    for (int k = 1; k <= nozn; ++k) r.pin_p[I(k)] = m.pin(k);
    for (int k = 1; k <= nozn; ++k)
      for (int i = its; i <= ite; ++i) r.o3clim_p(i, k) = m.o3clim(k, i);
    const int nlevs = kte - kts + 1, ncols = ite - its + 1;
    PArray p2d, o32d;  // (its:ite, kts:kte)
    p2d.allocate(ncols, nlevs);
    o32d.allocate(ncols, nlevs);
    for (int i = its; i <= ite; ++i)
      for (int k = kts; k <= kte; ++k) {
        o32d(i, k) = 0.;
        p2d(i, k) = v.pres_hyd_p(i, k) / 100.;
      }
    vinterp_ozn(1, ncols, ncols, nlevs, p2d.data(), r.pin_p.data(), nozn, r.o3clim_p.data(), o32d.data());
    for (int i = its; i <= ite; ++i)
      for (int k = kts; k <= kte; ++k) m.o3vmr(k, i) = o32d(i, k);
  } else {
    for (int k = 1; k <= nozn; ++k) r.pin_p[I(k)] = 0.0;
    for (int k = 1; k <= nozn; ++k)
      for (int i = its; i <= ite; ++i) r.o3clim_p(i, k) = 0.0;
  }

  RrtmgRadIn in;
  in.kte = kte;
  in.o3input = 0;
  if (m.config_o3climatology) in.o3input = 2;  // RL:847-848
  cloud_overlap(m, in.cldovrlp, in.idcor);     // RL:850-861
  if (in.idcor == 1 && !r.sw_ran)
    throw Error("driver_radiation_lw: 'latitude_varying' needs xlat_p, which only this step's shortwave sets");
  if (r.xlat_p.size() < static_cast<std::size_t>(n)) r.xlat_p.assign(static_cast<std::size_t>(n), 0.);
  in.p3d = &v.pres_hyd_p; in.p8w = &v.pres2_hyd_p; in.pi3d = &v.pi_p; in.t3d = &v.t_p; in.t8w = &v.t2_p;
  in.dz8w = &v.dz_p; in.qv3d = &r.qvrad_p; in.qc3d = &r.qcrad_p; in.qi3d = &r.qirad_p; in.qs3d = &r.qsrad_p;
  in.cldfra3d = &r.cldfrac_p; in.tsk = &r.tsk_p; in.xland = &r.xland_p; in.xice = &r.xice_p; in.snow = &r.snow_p;
  in.xlat = &r.xlat_p; in.julday = r.julday; in.icloud = r.icloud; in.noznlevels = nozn; in.pin = r.pin_p.data();
  in.o3clim = &r.o3clim_p; in.has_reqc = r.has_reqc; in.has_reqi = r.has_reqi; in.has_reqs = r.has_reqs;
  in.re_cloud = &r.recloud_p; in.re_ice = &r.reice_p; in.re_snow = &r.resnow_p;
  RrtmgLwradOut out{&r.glw_p,   &r.olrtoa_p, &r.lwcf_p,      &r.lwupt_p,  &r.lwuptc_p,  &r.lwdnt_p,
                    &r.lwdntc_p, &r.lwupb_p, &r.lwupbc_p,    &r.lwdnb_p,  &r.lwdnbc_p,  &r.rthratenlw_p,
                    &r.rrecloud_p, &r.rreice_p, &r.rresnow_p};
  rrtmg_lwrad(in, r.sfc_emiss_p, out, t.lw, its, ite);

  // radiation_lw_to_MPAS (RL:613-718)
  for (int i = its; i <= ite; ++i) {
    m.glw(i) = r.glw_p[I(i)];
    m.lwcf(i) = r.lwcf_p[I(i)];
    m.lwdnb(i) = r.lwdnb_p[I(i)];
    m.lwdnbc(i) = r.lwdnbc_p[I(i)];
    m.lwdnt(i) = r.lwdnt_p[I(i)];
    m.lwdntc(i) = r.lwdntc_p[I(i)];
    m.lwupb(i) = r.lwupb_p[I(i)];
    m.lwupbc(i) = r.lwupbc_p[I(i)];
    m.lwupt(i) = r.lwupt_p[I(i)];
    m.lwuptc(i) = r.lwuptc_p[I(i)];
    m.olrtoa(i) = r.olrtoa_p[I(i)];
  }
  for (int k = kts; k <= kte; ++k)
    for (int i = its; i <= ite; ++i) m.rthratenlw(k, i) = r.rthratenlw_p(i, k);
  if (microp_re_scheme) {  // RL:684-714: config_microp_re false
    for (int k = kts; k <= kte; ++k)
      for (int i = its; i <= ite; ++i) {
        m.rre_cloud(k, i) = 0.;
        m.rre_ice(k, i) = 0.;
        m.rre_snow(k, i) = 0.;
      }
  }
}

}  // namespace cppmpas::phys
