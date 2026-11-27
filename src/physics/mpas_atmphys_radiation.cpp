// SPDX-License-Identifier: GPL-3.0-or-later
// See mpas_atmphys_radiation.hpp. CL = physics/mpas_atmphys_driver_cloudiness.F, RS = mpas_atmphys_driver_radiation_sw.F,
// RL = mpas_atmphys_driver_radiation_lw.F, OZ = mpas_atmphys_o3climatology.F, MG = mpas_atmphys_manager.F,
// CS = physics_wrf/module_ra_cam_support.F (MPAS v8.4.2).
#include "physics/mpas_atmphys_radiation.hpp"

#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
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
  t.lw_store.resize(0);  // value-initialised in place (the table structs are large)
  t.lw_store.resize(1);
  t.sw_store.resize(0);
  t.sw_store.resize(1);
  t.k_store.resize(0);
  t.k_store.resize(1);
  if (!rrtmg::load_rrtmg_lw_tables(lw_tables_bin, t.lw_store[0]))
    throw Error("init_radiation_lw: cannot read the RRTMG LW tables " + lw_tables_bin);
  if (!rrtmg::load_rrtmg_sw_tables(sw_tables_bin, t.sw_store[0]))
    throw Error("init_radiation_sw: cannot read the RRTMG SW tables " + sw_tables_bin);
  rrtmg_mpas_fill_consts(t.k_store[0]);
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


// OZ:62-180: the tables are read on the host, the latitude interpolation of every cell is a kernel
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
  FStorage<real> lat_s(latsiz);
  {
    std::ifstream f = open("OZONE_LAT.TBL");
    for (int j = 1; j <= latsiz; ++j) lat_s[std::size_t(j - 1)] = read_real(f, "OZONE_LAT.TBL");  // OZ:119-122
  }
  // ozmixin(lonsiz, levsiz, latsiz, num_months) (OZ:134-143)
  FStorage<real> ozmixin(static_cast<std::size_t>(lonsiz) * levsiz * latsiz * num_months);
  const FView<real> OZI{ozmixin.data(), lonsiz, levsiz, latsiz};
  {
    std::ifstream f = open("OZONE_DAT.TBL");
    for (int mm = 1; mm <= num_months; ++mm)
      for (int j = 1; j <= latsiz; ++j)
        for (int k = 1; k <= levsiz; ++k)
          for (int i = 1; i <= lonsiz; ++i) OZI(i, k, j, mm) = read_real(f, "OZONE_DAT.TBL");
  }
  const real* const lat_ozone = lat_s.data();  // lat_ozone(j) at [j - 1]
  const FView<const real> latCell = m.latCell;
  const FView<real> ozmixm = m.ozmixm;
  const real dg = degrad;
  parallel_for(1, nCells, [=] CPPMPAS_LAMBDA(int iCell) {  // OZ:154-184
    const real lat = latCell(iCell) / dg;
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
          kernel_fail("init_o3climatology: a cell latitude equal to the last ozone latitude");
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
            ozmixm(mm, k, iCell) = OZI(j, k, i1, mm);
          } else {
            ozmixm(mm, k, iCell) = OZI(j, k, i1, mm) + (OZI(j, k, i2, mm) - OZI(j, k, i1, mm)) * dlatCell / dlat;
          }
        }
  });
}

// OZ:279-312 with getfactors (CS:2528-2611; r8 = double, ozncyc = .true.)
void o3climatology_months(real julian, int nMonths, int& nm, int& np, real& fact1, real& fact2) {
  static constexpr int date_oz[12] = {16, 45, 75, 105, 136, 166, 197, 228, 258, 289, 319, 350};
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
  fact1 = real(fact1_r8);
  fact2 = real(fact2_r8);
}

// OZ:228-325: the months and factors on the host, the interpolation of every cell a kernel (OZ:317-321)
void o3climatology_from_MPAS(real julian, RadiationPools& m, int nCellsSolve) {
  int nm = 0, np = 0;
  real fact1 = 0, fact2 = 0;
  o3climatology_months(julian, m.nMonths, nm, np, fact1, fact2);
  const int noz = m.nOznLevels;
  const FView<real> o3clim = m.o3clim, ozmixm = m.ozmixm;
  parallel_for(1, nCellsSolve, [=] CPPMPAS_LAMBDA(int iCell) {
    for (int k = 1; k <= noz; ++k) o3clim(k, iCell) = fact1 * ozmixm(nm, k, iCell) + fact2 * ozmixm(np, k, iCell);
  });
}

// ---- cloudiness ---------------------------------------------------------------------------------------------------
void driver_cloudiness(RadiationPools& m, const PhysicsVars& v, RadiationVars& r, int its, int ite) {
  if (its != 1) throw Error("driver_cloudiness: its /= 1");
  const int kts = 1, kte = m.nVertLevels, n = ite, kme = kte + 1;
  if (m.radt_cld_scheme != "cld_fraction")
    throw Error("driver_cloudiness: config_radt_cld_scheme = '" + m.radt_cld_scheme + "' has no verified port");
  // allocate_cloudiness (CL:64-75)
  for (PArray* a : {&r.cldfrac_p, &r.qvrad_p, &r.qcrad_p, &r.qirad_p, &r.qsrad_p}) a->allocate(n, kme);
  r.xland_cl.allocate(n, 1);
  const PView cldfrac_p = r.cldfrac_p.view(), qvrad_p = r.qvrad_p.view(), qcrad_p = r.qcrad_p.view(),
              qirad_p = r.qirad_p.view(), qsrad_p = r.qsrad_p.view(), xland_p = r.xland_cl.view();
  const PCView qv_p = v.qv_p.cview(), qc_p = v.qc_p.cview(), qi_p = v.qi_p.cview(), qs_p = v.qs_p.cview(),
               t_p = v.t_p.cview(), pres_p = v.pres_p.cview();
  const FView<const real> xland = m.xland;
  const FView<real> cldfrac = m.cldfrac;
  // cloudiness_from_MPAS (CL:92-143); dx_p (CL:124-126) left out (mpas_atmphys_radiation.hpp)
  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {
    xland_p(i) = xland(i);
    for (int k = kts; k <= kte; ++k) {
      qvrad_p(i, k) = qv_p(i, k);
      qcrad_p(i, k) = qc_p(i, k);
      qirad_p(i, k) = qi_p(i, k);
      qsrad_p(i, k) = qs_p(i, k);
      cldfrac_p(i, k) = 0.;
    }
  });
  // calc_cldfraction (CL:284-370): cldfrac_p, t_p, pres_p, qvrad_p, qcrad_p, qirad_p, qs_p
  const real ep_2 = physics_constants::ep_2;
  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {
    const real alpha0 = 100., gamma = 0.49, qcldmin = 1.e-12, pexp = 0.25, rhgrid = 1.0;
    const real svp1 = 0.61078, svp2 = 17.2693882, svpi2 = 21.8745584, svp3 = 35.86, svpi3 = 7.66, svpt0 = 273.15;
    for (int k = kts; k <= kte; ++k) cldfrac_p(i, k) = 0.;
    for (int k = kts; k <= kte; ++k) {
      const real t = t_p(i, k), p = pres_p(i, k);
      const real qv = qvrad_p(i, k), qc = qcrad_p(i, k), qi = qirad_p(i, k), qs = qs_p(i, k);
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
        cldfrac_p(i, k) = 0.;
      } else if (rhum >= rhgrid) {
        cldfrac_p(i, k) = 1.;
      } else {
        const real subsat = std::max(1.e-10, rhgrid * qvs - qv);
        const real denom = libm::pow(subsat, gamma);
        const real arg = std::max(-6.9, -alpha0 * qcld / denom);
        rhum = std::max(1.e-10, rhum);
        cldfrac_p(i, k) = libm::pow(rhum / rhgrid, pexp) * (1. - libm::exp(arg));  // CL:361
        if (cldfrac_p(i, k) < 0.01) cldfrac_p(i, k) = 0.;
      }
    }
  });
  // cloudiness_to_MPAS (CL:146-172)
  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {
    for (int k = kts; k <= kte; ++k) cldfrac(k, i) = cldfrac_p(i, k);
  });
}

namespace {

// the rrtmg_lwrad / rrtmg_swrad arguments of RL:858-878 / RS:922-946
RrtmgRadIn rrtmg_arguments(const RadiationPools& m, const PhysicsVars& v, const RadiationVars& r) {
  RrtmgRadIn in;
  in.kte = m.nVertLevels;
  in.o3input = 0;
  if (m.config_o3climatology) in.o3input = 2;  // RL:847-848, RS:901-902
  cloud_overlap(m, in.cldovrlp, in.idcor);     // RL:850-861, RS:904-915
  in.p3d = v.pres_hyd_p.cview(); in.p8w = v.pres2_hyd_p.cview(); in.pi3d = v.pi_p.cview(); in.t3d = v.t_p.cview();
  in.t8w = v.t2_p.cview(); in.dz8w = v.dz_p.cview(); in.qv3d = r.qvrad_p.cview(); in.qc3d = r.qcrad_p.cview();
  in.qi3d = r.qirad_p.cview(); in.qs3d = r.qsrad_p.cview(); in.cldfra3d = r.cldfrac_p.cview();
  in.tsk = r.tsk_p.cview(); in.xland = r.xland_p.cview(); in.xice = r.xice_p.cview(); in.snow = r.snow_p.cview();
  in.xlat = r.xlat_p.cview(); in.julday = r.julday; in.icloud = r.icloud; in.noznlevels = m.nOznLevels;
  in.pin = r.pin_p.cview(); in.o3clim = r.o3clim_p.cview();
  in.has_reqc = r.has_reqc; in.has_reqi = r.has_reqi; in.has_reqs = r.has_reqs;
  return in;
}

}  // namespace

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
  for (PArray* a : {&r.xlat_p, &r.xlon_p, &r.sfc_albedo_p, &r.snow_p, &r.tsk_p, &r.xice_p, &r.xland_p, &r.coszr_p,
                    &r.gsw_p, &r.swcf_p, &r.swdnb_p, &r.swdnbc_p, &r.swdnt_p, &r.swdntc_p, &r.swupb_p, &r.swupbc_p,
                    &r.swupt_p, &r.swuptc_p, &r.swddir_p, &r.swddni_p, &r.swddif_p})
    a->allocate(n, 1);
  for (PArray* a : {&r.f_ice, &r.f_rain, &r.rthratensw_p, &r.recloud_p, &r.reice_p, &r.resnow_p}) a->allocate(n, kme);
  r.pin_p.allocate(nozn, 1);
  r.o3clim_p.allocate(n, nozn);

  const real xtime_m = xtime_s / 60.;  // RS:874
  // radiation_sw_from_MPAS (RS:355-711)
  {
    const FView<const real> latCell = m.latCell, lonCell = m.lonCell, sfc_albedo = m.sfc_albedo, snow = m.snow,
                            skintemp = m.skintemp, xice = m.xice, xland = m.xland;
    const FView<real> cldfrac = m.cldfrac, o3clim = m.o3clim;
    const PView xlat_p = r.xlat_p.view(), xlon_p = r.xlon_p.view(), sfc_albedo_p = r.sfc_albedo_p.view(),
                snow_p = r.snow_p.view(), tsk_p = r.tsk_p.view(), xice_p = r.xice_p.view(),
                xland_p = r.xland_p.view(), cldfrac_p = r.cldfrac_p.view(), f_ice = r.f_ice.view(),
                f_rain = r.f_rain.view(), coszr_p = r.coszr_p.view(), gsw_p = r.gsw_p.view(),
                swcf_p = r.swcf_p.view(), swdnb_p = r.swdnb_p.view(), swdnbc_p = r.swdnbc_p.view(),
                swdnt_p = r.swdnt_p.view(), swdntc_p = r.swdntc_p.view(), swupb_p = r.swupb_p.view(),
                swupbc_p = r.swupbc_p.view(), swupt_p = r.swupt_p.view(), swuptc_p = r.swuptc_p.view(),
                rthratensw_p = r.rthratensw_p.view(), recloud_p = r.recloud_p.view(), reice_p = r.reice_p.view(),
                resnow_p = r.resnow_p.view(), o3clim_p = r.o3clim_p.view();
    const bool re_scheme = m.microp_scheme == "mp_thompson" || m.microp_scheme == "mp_wsm6";
    const bool o3c = m.config_o3climatology;
    const real dg = degrad;
    parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {
      xlat_p(i) = latCell(i) / dg;  // RS:419-429
      xlon_p(i) = lonCell(i) / dg;
      sfc_albedo_p(i) = sfc_albedo(i);
      snow_p(i) = snow(i);
      tsk_p(i) = skintemp(i);
      xice_p(i) = xice(i);
      xland_p(i) = xland(i);
      for (int k = kts; k <= kte; ++k) cldfrac_p(i, k) = cldfrac(k, i);  // RS:430-436
      for (int k = kts; k <= kte; ++k) {                                  // RS:439-445
        f_ice(i, k) = 0.;
        f_rain(i, k) = 0.;
      }
      coszr_p(i) = gsw_p(i) = swcf_p(i) = 0.;  // RS:448-462 (assigned 0 above)
      swdnb_p(i) = swdnbc_p(i) = swdnt_p(i) = swdntc_p(i) = 0.;
      swupb_p(i) = swupbc_p(i) = swupt_p(i) = swuptc_p(i) = 0.;
      for (int k = kts; k <= kte; ++k) rthratensw_p(i, k) = 0.;
      // aer_opt = 0 (RS:471); swddir_p, swddni_p, swddif_p = 0 (RS:476-482, assigned above); tauaer_p = 0,
      // ssaaer_p = 1, asyaer_p = 0 (RS:484-492): the aerosol values of CPPWRF's spcvmc_sw (module_ra_rrtmg.hpp)
      if (re_scheme) {  // RS:494-531, config_microp_re false
        for (int k = kts; k <= kte; ++k) {
          recloud_p(i, k) = 0.;
          reice_p(i, k) = 0.;
          resnow_p(i, k) = 0.;
        }
      }
      // swdnflx_p ... = 0 (RS:588-597): not passed to rrtmg_swrad, never copied to MPAS
      if (o3c) {  // RS:600-623
        for (int k = 1; k <= nozn; ++k) o3clim_p(i, k) = o3clim(k, i);
      } else {
        for (int k = 1; k <= nozn; ++k) o3clim_p(i, k) = 0.0;
      }
    });
    const FView<real> pin = m.pin;
    const PView pin_p = r.pin_p.view();
    parallel_for(1, nozn, [=] CPPMPAS_LAMBDA(int k) { pin_p(k) = o3c ? pin(k) : 0.0; });  // RS:600-623
  }
  r.sw_ran = true;

  radconst(r.declin, r.solcon, r.curr_julday, degrad, dpd);  // RS:884
  const real radt = r.dt_radtsw / 60.;                         // RS:896
  const RrtmgRadIn in = rrtmg_arguments(m, v, r);
  RrtmgSwradOut out;
  out.coszr = r.coszr_p.view(); out.gsw = r.gsw_p.view(); out.swcf = r.swcf_p.view(); out.swupt = r.swupt_p.view();
  out.swuptc = r.swuptc_p.view(); out.swdnt = r.swdnt_p.view(); out.swdntc = r.swdntc_p.view();
  out.swupb = r.swupb_p.view(); out.swupbc = r.swupbc_p.view(); out.swdnb = r.swdnb_p.view();
  out.swdnbc = r.swdnbc_p.view(); out.swddir = r.swddir_p.view(); out.swddni = r.swddni_p.view();
  out.swddif = r.swddif_p.view(); out.rthratensw = r.rthratensw_p.view();
  rrtmg_swrad(in, r.sfc_albedo_p.cview(), r.xlon_p.cview(), xtime_m, r.gmt, radt, degrad, r.declin, r.solcon, out,
              t.sw(), t.consts(), its, ite);

  // radiation_sw_to_MPAS (RS:714-794)
  {
    const PCView coszr_p = r.coszr_p.cview(), gsw_p = r.gsw_p.cview(), swcf_p = r.swcf_p.cview(),
                 swdnb_p = r.swdnb_p.cview(), swdnbc_p = r.swdnbc_p.cview(), swdnt_p = r.swdnt_p.cview(),
                 swdntc_p = r.swdntc_p.cview(), swupb_p = r.swupb_p.cview(), swupbc_p = r.swupbc_p.cview(),
                 swupt_p = r.swupt_p.cview(), swuptc_p = r.swuptc_p.cview(), swddir_p = r.swddir_p.cview(),
                 swddni_p = r.swddni_p.cview(), swddif_p = r.swddif_p.cview(),
                 rthratensw_p = r.rthratensw_p.cview();
    const FView<real> coszr = m.coszr, gsw = m.gsw, swcf = m.swcf, swdnb = m.swdnb, swdnbc = m.swdnbc,
                      swdnt = m.swdnt, swdntc = m.swdntc, swupb = m.swupb, swupbc = m.swupbc, swupt = m.swupt,
                      swuptc = m.swuptc, swddir = m.swddir, swddni = m.swddni, swddif = m.swddif,
                      rthratensw = m.rthratensw;
    parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {
      coszr(i) = coszr_p(i);
      gsw(i) = gsw_p(i);
      swcf(i) = swcf_p(i);
      swdnb(i) = swdnb_p(i);
      swdnbc(i) = swdnbc_p(i);
      swdnt(i) = swdnt_p(i);
      swdntc(i) = swdntc_p(i);
      swupb(i) = swupb_p(i);
      swupbc(i) = swupbc_p(i);
      swupt(i) = swupt_p(i);
      swuptc(i) = swuptc_p(i);
      swddir(i) = swddir_p(i);  // RS:774-780
      swddni(i) = swddni_p(i);
      swddif(i) = swddif_p(i);
      for (int k = kts; k <= kte; ++k) rthratensw(k, i) = rthratensw_p(i, k);
    });
  }
}

// ---- longwave -----------------------------------------------------------------------------------------------------
void driver_radiation_lw(real xtime_s, RadiationPools& m, const PhysicsVars& v, RadiationVars& r,
                         const RrtmgTables& t, int its, int ite) {
  if (its != 1) throw Error("driver_radiation_lw: its /= 1");
  if (m.radt_lw_scheme != "rrtmg_lw")
    throw Error("driver_radiation_lw: config_radt_lw_scheme = '" + m.radt_lw_scheme + "' has no verified port");
  check_microp(m);
  if (r.cldfrac_p.ldi != ite) throw Error("driver_radiation_lw: driver_cloudiness has not run");
  const int kts = 1, kte = m.nVertLevels, n = ite, kme = kte + 1;
  const int nozn = m.nOznLevels;
  // allocate_radiation_lw (RL:103-206)
  for (PArray* a : {&r.sfc_emiss_p, &r.snow_p, &r.tsk_p, &r.xice_p, &r.xland_p, &r.glw_p, &r.lwcf_p, &r.lwdnb_p,
                    &r.lwdnbc_p, &r.lwdnt_p, &r.lwdntc_p, &r.lwupb_p, &r.lwupbc_p, &r.lwupt_p, &r.lwuptc_p,
                    &r.olrtoa_p})
    a->allocate(n, 1);
  for (PArray* a : {&r.f_ice, &r.f_rain, &r.rthratenlw_p, &r.recloud_p, &r.reice_p, &r.resnow_p, &r.rrecloud_p,
                    &r.rreice_p, &r.rresnow_p})
    a->allocate(n, kme);
  r.pin_p.allocate(nozn, 1);
  r.o3clim_p.allocate(n, nozn);

  const bool microp_re_scheme =
      m.microp_scheme == "mp_thompson" || m.microp_scheme == "mp_thompson_aerosols" || m.microp_scheme == "mp_wsm6";
  const bool o3c = m.config_o3climatology;
  // radiation_lw_from_MPAS (RL:297-610)
  {
    const FView<const real> sfc_emiss = m.sfc_emiss, skintemp = m.skintemp, snow = m.snow, xice = m.xice,
                            xland = m.xland;
    const FView<real> cldfrac = m.cldfrac, o3clim = m.o3clim;
    const PView sfc_emiss_p = r.sfc_emiss_p.view(), tsk_p = r.tsk_p.view(), snow_p = r.snow_p.view(),
                xice_p = r.xice_p.view(), xland_p = r.xland_p.view(), cldfrac_p = r.cldfrac_p.view(),
                f_ice = r.f_ice.view(), f_rain = r.f_rain.view(), glw_p = r.glw_p.view(), lwcf_p = r.lwcf_p.view(),
                lwdnb_p = r.lwdnb_p.view(), lwdnbc_p = r.lwdnbc_p.view(), lwdnt_p = r.lwdnt_p.view(),
                lwdntc_p = r.lwdntc_p.view(), lwupb_p = r.lwupb_p.view(), lwupbc_p = r.lwupbc_p.view(),
                lwupt_p = r.lwupt_p.view(), lwuptc_p = r.lwuptc_p.view(), olrtoa_p = r.olrtoa_p.view(),
                rthratenlw_p = r.rthratenlw_p.view(), recloud_p = r.recloud_p.view(), reice_p = r.reice_p.view(),
                resnow_p = r.resnow_p.view(), rrecloud_p = r.rrecloud_p.view(), rreice_p = r.rreice_p.view(),
                rresnow_p = r.rresnow_p.view(), o3clim_p = r.o3clim_p.view();
    parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {
      sfc_emiss_p(i) = sfc_emiss(i);  // RL:364-372
      tsk_p(i) = skintemp(i);
      snow_p(i) = snow(i);
      xice_p(i) = xice(i);
      xland_p(i) = xland(i);
      for (int k = kts; k <= kte; ++k) cldfrac_p(i, k) = cldfrac(k, i);  // RL:373-379
      for (int k = kts; k <= kte; ++k) {                                  // RL:382-389
        f_ice(i, k) = 0.;
        f_rain(i, k) = 0.;
      }
      glw_p(i) = lwcf_p(i) = lwdnb_p(i) = lwdnbc_p(i) = lwdnt_p(i) = 0.;  // RL:391-404
      lwdntc_p(i) = lwupb_p(i) = lwupbc_p(i) = lwupt_p(i) = lwuptc_p(i) = 0.;
      olrtoa_p(i) = 0.;
      for (int k = kts; k <= kte; ++k) rthratenlw_p(i, k) = 0.;  // RL:405-409
      if (microp_re_scheme) {  // RL:416-458 (config_microp_re false)
        for (int k = kts; k <= kte; ++k) {
          recloud_p(i, k) = 0.;
          reice_p(i, k) = 0.;
          resnow_p(i, k) = 0.;
        }
        for (int k = kts; k <= kte; ++k) {
          rrecloud_p(i, k) = 0.;
          rreice_p(i, k) = 0.;
          rresnow_p(i, k) = 0.;
        }
      }
      // lwdnflx_p ... = 0 (RL:459-467): not passed to rrtmg_lwrad, never copied to MPAS
      if (o3c) {  // RL:470-516
        for (int k = 1; k <= nozn; ++k) o3clim_p(i, k) = o3clim(k, i);
      } else {
        for (int k = 1; k <= nozn; ++k) o3clim_p(i, k) = 0.0;
      }
    });
    const FView<real> pin = m.pin;
    const PView pin_p = r.pin_p.view();
    parallel_for(1, nozn, [=] CPPMPAS_LAMBDA(int k) { pin_p(k) = o3c ? pin(k) : 0.0; });
  }
  if (o3c) {  // the o3vmr diagnostic (RL:476-500)
    const int nlevs = kte - kts + 1, ncols = ite - its + 1;
    if (nozn > kMaxOznLevels) throw Error("driver_radiation_lw: more than 64 ozone levels");
    // vinterp_ozn runs column by column: its multi-column search (VI:95-153) gives each column the values of its
    // own search when OZONE_PLEV's pressures are non-decreasing (one interval per pressure) and the column's
    // pressures are non-decreasing downwards (no column's interval below the last one it found)
    for (int k = 1; k < nozn; ++k)
      if (r.pin_p(k + 1) < r.pin_p(k))
        throw Error("driver_radiation_lw: OZONE_PLEV.TBL's pressures are not ordered (vinterp_ozn's multi-column "
                    "search has no verified port then)");
    PArray p2d, o32d;  // (its:ite, kts:kte)
    p2d.allocate(ncols, nlevs);
    o32d.allocate(ncols, nlevs);
    const PView p2dv = p2d.view(), o32dv = o32d.view(), o3clim_p = r.o3clim_p.view();
    const PCView pres_hyd_p = v.pres_hyd_p.cview(), pin_p = r.pin_p.cview();
    const FView<real> o3vmr = m.o3vmr;
    parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {
      for (int k = kts; k <= kte; ++k) {
        o32dv(i, k) = 0.;
        p2dv(i, k) = pres_hyd_p(i, k) / 100.;
      }
      for (int k = kts; k < kte; ++k)
        if (p2dv(i, k + 1) > p2dv(i, k))
          kernel_fail("driver_radiation_lw: a column whose hydrostatic pressure grows upwards (the o3vmr diagnostic's "
                      "multi-column vinterp_ozn has no verified port then)");
      vinterp_ozn_column(
          nlevs, [&](int k) { return p2dv(i, k); }, [&](int k) { return pin_p(k); }, nozn,
          [&](int k) { return o3clim_p(i, k); }, [&](int k) -> real& { return o32dv(i, k); });
      for (int k = kts; k <= kte; ++k) o3vmr(k, i) = o32dv(i, k);
    });
  }

  RrtmgRadIn in = rrtmg_arguments(m, v, r);
  if (in.idcor == 1 && !r.sw_ran)
    throw Error("driver_radiation_lw: 'latitude_varying' needs xlat_p, which only this step's shortwave sets");
  if (r.xlat_p.ldi < n) {
    r.xlat_p.allocate(n, 1);
    in.xlat = r.xlat_p.cview();
  }
  RrtmgLwradOut out;
  out.glw = r.glw_p.view(); out.olr = r.olrtoa_p.view(); out.lwcf = r.lwcf_p.view(); out.lwupt = r.lwupt_p.view();
  out.lwuptc = r.lwuptc_p.view(); out.lwdnt = r.lwdnt_p.view(); out.lwdntc = r.lwdntc_p.view();
  out.lwupb = r.lwupb_p.view(); out.lwupbc = r.lwupbc_p.view(); out.lwdnb = r.lwdnb_p.view();
  out.lwdnbc = r.lwdnbc_p.view(); out.rthratenlw = r.rthratenlw_p.view(); out.rre_cloud = r.rrecloud_p.view();
  out.rre_ice = r.rreice_p.view(); out.rre_snow = r.rresnow_p.view();
  rrtmg_lwrad(in, r.sfc_emiss_p.cview(), out, t.lw(), t.consts(), its, ite);

  // radiation_lw_to_MPAS (RL:613-718)
  {
    const PCView glw_p = r.glw_p.cview(), lwcf_p = r.lwcf_p.cview(), lwdnb_p = r.lwdnb_p.cview(),
                 lwdnbc_p = r.lwdnbc_p.cview(), lwdnt_p = r.lwdnt_p.cview(), lwdntc_p = r.lwdntc_p.cview(),
                 lwupb_p = r.lwupb_p.cview(), lwupbc_p = r.lwupbc_p.cview(), lwupt_p = r.lwupt_p.cview(),
                 lwuptc_p = r.lwuptc_p.cview(), olrtoa_p = r.olrtoa_p.cview(),
                 rthratenlw_p = r.rthratenlw_p.cview();
    const FView<real> glw = m.glw, lwcf = m.lwcf, lwdnb = m.lwdnb, lwdnbc = m.lwdnbc, lwdnt = m.lwdnt,
                      lwdntc = m.lwdntc, lwupb = m.lwupb, lwupbc = m.lwupbc, lwupt = m.lwupt, lwuptc = m.lwuptc,
                      olrtoa = m.olrtoa, rthratenlw = m.rthratenlw, rre_cloud = m.rre_cloud, rre_ice = m.rre_ice,
                      rre_snow = m.rre_snow;
    parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {
      glw(i) = glw_p(i);
      lwcf(i) = lwcf_p(i);
      lwdnb(i) = lwdnb_p(i);
      lwdnbc(i) = lwdnbc_p(i);
      lwdnt(i) = lwdnt_p(i);
      lwdntc(i) = lwdntc_p(i);
      lwupb(i) = lwupb_p(i);
      lwupbc(i) = lwupbc_p(i);
      lwupt(i) = lwupt_p(i);
      lwuptc(i) = lwuptc_p(i);
      olrtoa(i) = olrtoa_p(i);
      for (int k = kts; k <= kte; ++k) rthratenlw(k, i) = rthratenlw_p(i, k);
      if (microp_re_scheme) {  // RL:684-714: config_microp_re false
        for (int k = kts; k <= kte; ++k) {
          rre_cloud(k, i) = 0.;
          rre_ice(k, i) = 0.;
          rre_snow(k, i) = 0.;
        }
      }
    });
  }

  if (const char* dir = std::getenv("CPPMPAS_RADIATION_DUMP")) {  // tests/radiation_replay
    char name[64];
    std::snprintf(name, sizeof name, "/radiation.%lld.bin", static_cast<long long>(std::llround(xtime_s)));
    write_radiation_dump(std::string(dir) + name, xtime_s, m, v, r, ite);
  }
}

// ---- the radiation-replay dump ------------------------------------------------------------------------------------
void write_radiation_dump(const std::string& path, real xtime_s, const RadiationPools& m, const PhysicsVars& v,
                          const RadiationVars& r, int ite) {
  std::FILE* f = std::fopen(path.c_str(), "wb");
  if (!f) throw Error("write_radiation_dump: cannot write " + path);
  auto rec = [&](const char* name, int kind, const void* p, std::int64_t count, std::size_t elem) {
    char nm[32] = {};
    std::strncpy(nm, name, sizeof nm - 1);
    const std::int32_t k = kind;
    if (std::fwrite(nm, 1, 32, f) != 32 || std::fwrite(&k, 4, 1, f) != 1 || std::fwrite(&count, 8, 1, f) != 1 ||
        std::fwrite(p, elem, std::size_t(count), f) != std::size_t(count))
      throw Error("write_radiation_dump: short write to " + path);
  };
  auto ints = [&](const char* name, const std::vector<std::int32_t>& x) { rec(name, 1, x.data(), std::int64_t(x.size()), 4); };
  auto reals = [&](const char* name, const std::vector<double>& x) { rec(name, 0, x.data(), std::int64_t(x.size()), 8); };
  auto str = [&](const char* name, const std::string& s) { rec(name, 2, s.data(), std::int64_t(s.size()), 1); };
  const int nV = m.nVertLevels, noz = m.nOznLevels;
  auto cells = [&](const FView<const real>& a) {  // a(i), i = 1..ite
    std::vector<double> x(static_cast<std::size_t>(ite));
    for (int i = 1; i <= ite; ++i) x[std::size_t(i - 1)] = a(i);
    return x;
  };
  auto levels = [&](const FView<const real>& a, int nk) {  // a(k, i), k = 1..nk, i = 1..ite
    std::vector<double> x(static_cast<std::size_t>(nk) * ite);
    for (int i = 1; i <= ite; ++i)
      for (int k = 1; k <= nk; ++k) x[std::size_t(k - 1) + std::size_t(nk) * (i - 1)] = a(k, i);
    return x;
  };
  auto cv = [](FView<real> a) { return FView<const real>{a.p, a.e0, a.e1, a.e2}; };
  int nm = 0, np = 0;
  real fact1 = 0, fact2 = 0;
  if (m.config_o3climatology) o3climatology_months(r.curr_julday, m.nMonths, nm, np, fact1, fact2);
  ints("ints", {ite, nV, noz, m.nMonths, r.julday, r.year, r.icloud, r.has_reqc, r.has_reqi, r.has_reqs,
                r.sw_ran ? 1 : 0, m.config_o3climatology ? 1 : 0, m.config_microp_re ? 1 : 0, v.ime, v.kme, nm, np});
  reals("reals", {xtime_s, r.gmt, r.curr_julday, r.dt_radtsw, r.dt_radtlw});
  str("radt_lw_scheme", m.radt_lw_scheme);
  str("radt_sw_scheme", m.radt_sw_scheme);
  str("radt_cld_scheme", m.radt_cld_scheme);
  str("radt_cld_overlap", m.radt_cld_overlap);
  str("radt_cld_dcorrlen", m.radt_cld_dcorrlen);
  str("microp_scheme", m.microp_scheme);
  reals("latCell", cells(m.latCell));
  reals("lonCell", cells(m.lonCell));
  reals("skintemp", cells(m.skintemp));
  reals("snow", cells(m.snow));
  reals("xice", cells(m.xice));
  reals("xland", cells(m.xland));
  reals("sfc_albedo", cells(m.sfc_albedo));
  reals("sfc_emiss", cells(m.sfc_emiss));
  if (m.config_o3climatology) {
    std::vector<double> pin(static_cast<std::size_t>(noz)), om(std::size_t(noz) * ite), op(std::size_t(noz) * ite);
    for (int k = 1; k <= noz; ++k) pin[std::size_t(k - 1)] = m.pin(k);
    for (int i = 1; i <= ite; ++i)
      for (int k = 1; k <= noz; ++k) {
        om[std::size_t(k - 1) + std::size_t(noz) * (i - 1)] = m.ozmixm(nm, k, i);
        op[std::size_t(k - 1) + std::size_t(noz) * (i - 1)] = m.ozmixm(np, k, i);
      }
    reals("pin", pin);
    reals("ozmixm_nm", om);
    reals("ozmixm_np", op);
  }
  auto parr = [&](const char* name, const PArray& a) { rec(name, 0, a.data(), std::int64_t(a.a.size()), 8); };
  parr("@pres_hyd_p", v.pres_hyd_p); parr("@pres2_hyd_p", v.pres2_hyd_p); parr("@pi_p", v.pi_p);
  parr("@t_p", v.t_p); parr("@t2_p", v.t2_p); parr("@dz_p", v.dz_p); parr("@pres_p", v.pres_p);
  parr("@qv_p", v.qv_p); parr("@qc_p", v.qc_p); parr("@qi_p", v.qi_p); parr("@qs_p", v.qs_p);
  // the outputs
  if (m.config_o3climatology) {
    reals("o3clim", levels(cv(m.o3clim), noz));
    reals("o3vmr", levels(cv(m.o3vmr), nV));
  }
  reals("cldfrac", levels(cv(m.cldfrac), nV));
  reals("rre_cloud", levels(cv(m.rre_cloud), nV));
  reals("rre_ice", levels(cv(m.rre_ice), nV));
  reals("rre_snow", levels(cv(m.rre_snow), nV));
  reals("rthratenlw", levels(cv(m.rthratenlw), nV));
#define CPPMPAS_D(x) reals(#x, cells(cv(m.x)));
  CPPMPAS_D(glw) CPPMPAS_D(lwcf) CPPMPAS_D(lwdnb) CPPMPAS_D(lwdnbc) CPPMPAS_D(lwdnt) CPPMPAS_D(lwdntc)
  CPPMPAS_D(lwupb) CPPMPAS_D(lwupbc) CPPMPAS_D(lwupt) CPPMPAS_D(lwuptc) CPPMPAS_D(olrtoa)
  if (r.sw_ran) {
    reals("rthratensw", levels(cv(m.rthratensw), nV));
    CPPMPAS_D(coszr) CPPMPAS_D(gsw) CPPMPAS_D(swcf) CPPMPAS_D(swdnb) CPPMPAS_D(swdnbc) CPPMPAS_D(swdnt)
    CPPMPAS_D(swdntc) CPPMPAS_D(swupb) CPPMPAS_D(swupbc) CPPMPAS_D(swupt) CPPMPAS_D(swuptc) CPPMPAS_D(swddir)
    CPPMPAS_D(swddni) CPPMPAS_D(swddif)
  }
#undef CPPMPAS_D
  if (std::fclose(f) != 0) throw Error("write_radiation_dump: cannot close " + path);
}

}  // namespace cppmpas::phys
