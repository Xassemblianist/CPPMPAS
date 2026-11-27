// SPDX-License-Identifier: GPL-3.0-or-later
// See mpas_atmphys_landuse.hpp. LU = physics/mpas_atmphys_landuse.F, US = physics/mpas_atmphys_update_surface.F,
// DT = physics/mpas_atmphys_date_time.F (MPAS v8.4.2).
#include "physics/mpas_atmphys_landuse.hpp"

#include <cstdio>
#include <sstream>

#include "core/error.hpp"
#include "physics/fortran_list_read.hpp"

namespace cppmpas::phys {

int landuse_init_forMPAS(const std::string& table_dir, int julday, const LanduseConfig& cfg, const Mesh& mesh,
                         RegistryPools& p, std::vector<std::string>& log) {
  constexpr int max_cats = 100, max_seas = 12;  // LU:110-111
  // albd(max_cats, max_seas), ...: (ic, is) 1-based
  std::vector<real> albd(max_cats * max_seas), slmo(albd.size()), sfem(albd.size()), sfz0(albd.size()),
      therin(albd.size()), scfx(albd.size()), sfhc(albd.size());
  auto A = [](std::vector<real>& a, int ic, int is) -> real& {
    return a[static_cast<std::size_t>(is - 1) * max_cats + (ic - 1)];
  };
  const std::string mminlu = p.C("mminlu");
  std::string lutype;
  int lucats = 0, luseas = 0;
  {  // LU:142-189
    FortranListFile u((table_dir.empty() ? std::string(".") : table_dir) + "/LANDUSE.TBL", "LANDUSE.TBL",
                      "subroutine landuse_init_forMPAS");
    int lumatch = 0;
    while (lumatch == 0) {
      if (!u.read_an(lutype, 35)) throw Error("landuse_init_forMPAS: end of file reading LANDUSE.TBL");  // :156
      if (!u.read_list(&lucats, &luseas)) throw Error("landuse_init_forMPAS: end of file reading LANDUSE.TBL");
      if (feq(lutype, mminlu)) {
        std::ostringstream m;
        m << "    landuse type = " << lutype.substr(0, lutype.find_last_not_of(' ') + 1) << " found " << lucats
          << " categories " << luseas << " seasons";
        log.push_back(m.str());
        lumatch = 1;
      } else {
        log.push_back("    skipping over lutype = " + lutype.substr(0, lutype.find_last_not_of(' ') + 1));
        for (int is = 1; is <= luseas; ++is) {  // :167-172 (iostat=ierr: an end of file is not an error)
          if (!u.eof()) u.skip();
          for (int ic = 1; ic <= lucats; ++ic)
            if (!u.eof()) u.skip();
        }
      }
    }
    if (lucats > max_cats || luseas > max_seas) throw Error("landuse_init_forMPAS: more categories than max_cats");
    for (int is = 1; is <= luseas; ++is) {  // :176-183
      if (!u.eof()) u.skip();
      for (int ic = 1; ic <= lucats; ++ic) {
        real li = 0;
        if (!u.read_list(&li, &A(albd, ic, is), &A(slmo, ic, is), &A(sfem, ic, is), &A(sfz0, ic, is),
                         &A(therin, ic, is), &A(scfx, ic, is), &A(sfhc, ic, is)))
          throw Error("landuse_init_forMPAS: end of file reading LANDUSE.TBL");
      }
    }
  }
  int isurban = 0;  // the module variable keeps its value for an unknown lutype (LU:189-205 case default)
  const std::string lt = lutype.substr(0, lutype.find_last_not_of(' ') + 1);
  if (lt == "OLD" || lt == "USGS" || lt == "LW12")
    isurban = 1;
  else if (lt == "MODIFIED_IGBP_MODIS_NOAH")
    isurban = 13;
  else if (lt == "SiB")
    isurban = 11;
  if (cfg.config_do_restart) {  // LU:226-230
    log.push_back("--- config_do_restart =T");
    log.push_back("--- skip the end of landuse_init_forMPAS");
    return isurban;
  }
  const int isice = p.I("isice_lu")(1), iswater = p.I("iswater_lu")(1);
  auto& ivgtyp = p.I("ivgtyp");
  auto &snoalb = p.R("snoalb"), &snowc = p.R("snowc"), &xice = p.R("xice"), &albbck = p.R("sfc_albbck");
  auto &embck = p.R("sfc_emibck"), &sfc_albedo = p.R("sfc_albedo"), &sfc_emiss = p.R("sfc_emiss");
  auto &thc = p.R("thc"), &z0 = p.R("z0"), &znt = p.R("znt");
  // mavail and ust belong to the 'sfclayer' package: associated only when it is active (LU:127-128, 133-136)
  FArray<real>* mavail = p.packages.active(registry_var("mavail")) ? &p.R("mavail") : nullptr;
  FArray<real>* ust = p.packages.active(registry_var("ust")) ? &p.R("ust") : nullptr;
  for (int iCell = 1; iCell <= mesh.nCells; ++iCell) {  // LU:233-280
    int isn = 1;
    if (julday < 105 || julday >= 288) isn = 2;
    if (mesh.latCell(iCell) < 0.) isn = 3 - isn;
    int is = ivgtyp(iCell);
    if (is == 0) is = iswater;
    if (is < 1 || is > lucats || isice < 1 || isice > lucats)
      throw Error("landuse_init_forMPAS: a land-use category outside the table (MPAS would read out of bounds)");
    if (!cfg.config_sfc_albedo) albbck(iCell) = A(albd, is, isn) / 100.;
    sfc_albedo(iCell) = albbck(iCell);
    if (snowc(iCell) > 0.5) {
      albbck(iCell) = A(albd, isice, isn) / 100.;
      if (cfg.config_sfc_albedo)
        sfc_albedo(iCell) = snoalb(iCell);
      else
        sfc_albedo(iCell) = albbck(iCell) / (1 + A(scfx, is, isn));
    }
    thc(iCell) = A(therin, is, isn) / 100.;
    z0(iCell) = A(sfz0, is, isn) / 100.;
    znt(iCell) = z0(iCell);
    embck(iCell) = A(sfem, is, isn);
    sfc_emiss(iCell) = embck(iCell);
    if (mavail) (*mavail)(iCell) = A(slmo, is, isn);
    if (ust) (*ust)(iCell) = 0.0001;
    if (xice(iCell) >= cfg.xice_threshold) {
      albbck(iCell) = A(albd, isice, isn) / 100.;
      embck(iCell) = A(sfem, isice, isn);
      if (cfg.config_frac_seaice) {
        sfc_albedo(iCell) = xice(iCell) * albbck(iCell) + (1 - xice(iCell)) * 0.08;
        sfc_emiss(iCell) = xice(iCell) * embck(iCell) + (1 - xice(iCell)) * 0.98;
      } else {
        sfc_albedo(iCell) = albbck(iCell);
        sfc_emiss(iCell) = embck(iCell);
      }
      thc(iCell) = A(therin, isice, isn) / 100.;
      z0(iCell) = A(sfz0, isice, isn) / 100.;
      znt(iCell) = z0(iCell);
      if (mavail) (*mavail)(iCell) = A(slmo, isice, isn);
    }
  }
  return isurban;
}

namespace {

// DT:131-149: read(date, '(nX, In)'); blanks read as zeros (an internal unit, BLANK='NULL' with nothing but blanks)
int field(const std::string& d, std::size_t pos, std::size_t n) {
  int v = 0;
  bool neg = false, any = false;
  for (std::size_t i = pos; i < pos + n; ++i) {
    const char c = i < d.size() ? d[i] : ' ';
    if (c == ' ') continue;
    if (c == '-' && !any) {
      neg = true;
      continue;
    }
    if (c < '0' || c > '9') throw Error("split_date_char: '" + d + "' is not a date");
    v = v * 10 + (c - '0');
    any = true;
  }
  return neg ? -v : v;
}

// DT:61-62: data mmd /.../ makes mmd a saved variable: a leap year once seen keeps February at 29 days
int mmd[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

}  // namespace

void get_julgmt(const std::string& date_str, int& julyr, int& julday, real& gmt) {
  const int ny = field(date_str, 0, 4), nm = field(date_str, 5, 2), nd = field(date_str, 8, 2);
  const int nh = field(date_str, 11, 2), ni = field(date_str, 14, 2), ns = field(date_str, 17, 2);
  gmt = nh + real(ni) / 60. + real(ns) / 3600.;  // DT:68 (float: default real, 8 bytes)
  const int my1 = ny % 4, my2 = ny % 100, my3 = ny % 400;
  if ((my1 == 0 && my2 != 0) || my3 == 0) mmd[1] = 29;
  julday = nd;
  julyr = ny;
  for (int monss = 1; monss <= nm - 1; ++monss) julday = julday + mmd[monss - 1];
}

namespace {

// DT:152-225; field_in(12, npoints)
void monthly_interp_to_date(int npoints, const std::string& date_str, const FArray<real>& field_in,
                            FArray<real>& field_out) {
  int middle[14];
  int julyr = 0, julday = 0;
  real gmt = 0;
  for (int l = 1; l <= 12; ++l) {
    char mon[3];
    std::snprintf(mon, sizeof mon, "%02d", l);
    get_julgmt(date_str.substr(0, 4) + "-" + mon + "-15_00:00:00.0000", julyr, julday, gmt);
    middle[l] = julyr * 1000 + julday;
  }
  middle[0] = middle[1] - 31;
  middle[13] = middle[12] + 31;
  int target_julyr = 0, target_julday = 0;
  get_julgmt(date_str, target_julyr, target_julday, gmt);
  const int target_date = target_julyr * 1000 + target_julday;
  for (int l = 0; l <= 12; ++l) {
    if (middle[l] < target_date && middle[l + 1] >= target_date) {
      for (int n = 1; n <= npoints; ++n) {
        const int int_month = l;
        int month1, month2;
        if (int_month == 0 || int_month == 12) {
          month1 = 12;
          month2 = 1;
        } else {
          month1 = int_month;
          month2 = month1 + 1;
        }
        field_out(n) = (field_in(month2, n) * real(target_date - middle[l]) +
                        field_in(month1, n) * real(middle[l + 1] - target_date)) /
                       real(middle[l + 1] - middle[l]);
      }
      break;
    }
  }
}

// DT:228-258
void monthly_min_max(int npoints, const FArray<real>& field_in, FArray<real>& field_min, FArray<real>& field_max) {
  for (int n = 1; n <= npoints; ++n) {
    real minner = field_in(1, n), maxxer = field_in(1, n);
    for (int nn = 2; nn <= 12; ++nn) {
      if (field_in(nn, n) < minner) minner = field_in(nn, n);
      if (field_in(nn, n) > maxxer) maxxer = field_in(nn, n);
      field_min(n) = minner;
      field_max(n) = maxxer;
    }
  }
}

}  // namespace

void physics_update_surface(const std::string& current_date, bool config_sfc_albedo, int nCellsSolve,
                            RegistryPools& p) {
  auto& landmask = p.I("landmask");
  auto &albedo12m = p.R("albedo12m"), &sfc_albbck = p.R("sfc_albbck"), &greenfrac = p.R("greenfrac");
  auto &vegfra = p.R("vegfra"), &shdmin = p.R("shdmin"), &shdmax = p.R("shdmax");
  if (albedo12m.extent(0) != 12 || greenfrac.extent(0) != 12)
    throw Error("physics_update_surface: nMonths != 12 (monthly_interp_to_date's dimension(12, npoints))");
  if (config_sfc_albedo) {  // US:85-93
    monthly_interp_to_date(nCellsSolve, current_date, albedo12m, sfc_albbck);
    for (int iCell = 1; iCell <= nCellsSolve; ++iCell) {
      sfc_albbck(iCell) = sfc_albbck(iCell) / 100.;
      if (landmask(iCell) == 0) sfc_albbck(iCell) = 0.08;
    }
  }
  monthly_interp_to_date(nCellsSolve, current_date, greenfrac, vegfra);  // US:96
  monthly_min_max(nCellsSolve, greenfrac, shdmin, shdmax);             // US:99
}

}  // namespace cppmpas::phys
