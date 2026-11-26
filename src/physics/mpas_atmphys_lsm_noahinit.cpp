// SPDX-License-Identifier: GPL-3.0-or-later
// soil_veg_gen_parm of mpas_atmphys_lsm_noahinit.F:210-589; see mpas_atmphys_lsm_noahinit.hpp. lsminit is in
// mpas_atmphys_driver_lsm.cpp (init_lsm).
#include "physics/mpas_atmphys_lsm_noahinit.hpp"

#include <cerrno>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>

#include "core/error.hpp"

namespace cppmpas::phys {

namespace {

// A formatted sequential file as gfortran reads it: records are lines; a list-directed READ (fmt = *) starts at
// the next record, takes its items from as many records as it needs and skips the rest of its last record.
// Value separators: a comma or a slash with optional blanks around it, or blanks. A slash ends the READ (the
// remaining items keep their values). An undelimited character value ends at a blank, comma or slash; a delimited
// one is between apostrophes or quotes. A null value (two commas) is refused: no Noah table has one and MPAS's
// READs would leave the item unchanged.
class Unit {
 public:
  Unit(const std::string& path, const std::string& what) : name_(what) {
    std::ifstream f(path);
    if (!f) throw Error("mpas_atmphys_lsm_noahinit.F: soil_veg_gen_parm: failure opening " + what);
    std::string l;
    while (std::getline(f, l)) {
      if (!l.empty() && l.back() == '\r') l.pop_back();
      lines_.push_back(l);
    }
  }
  bool eof() const { return rec_ >= lines_.size(); }
  // READ (u, '(A)') s: the whole record. false at the end of the file (END=).
  bool read_a(std::string& s) {
    if (eof()) return false;
    s = lines_[rec_++];
    return true;
  }
  // READ (u, '(A4)') s: the first 4 characters, blank padded
  bool read_a4(std::string& s) {
    if (eof()) return false;
    s = lines_[rec_++].substr(0, 4);
    return true;
  }
  // READ (u, *) with no items: skips one record
  void skip() {
    if (eof()) throw Error("soil_veg_gen_parm: end of file reading " + name_);
    ++rec_;
  }
  // READ (u, *) items: each item a pointer to int, double or std::string. false when the READ meets the end of
  // the file (the END= branch; the callers without END= make it an error, as the Fortran stops)
  template <class... T> bool read_list(T*... items) {
    if (eof()) return false;
    line_ = lines_[rec_++];
    pos_ = 0;
    slash_ = false;
    comma_ = true;
    try {
      (item(items), ...);
    } catch (const EndOfFile&) {
      return false;
    }
    return true;
  }

 private:
  struct EndOfFile {};
  [[noreturn]] void end_error() const { throw EndOfFile{}; }
  // the next value as a string; false after a slash. comma_ tells whether a comma has been taken since the last
  // value (true at the start of the READ): a comma then is a null value.
  bool next(std::string& v) {
    if (slash_) return false;
    for (;;) {
      while (pos_ < line_.size() && (line_[pos_] == ' ' || line_[pos_] == '\t')) ++pos_;
      if (pos_ >= line_.size()) {  // the end of a record is a blank: the READ continues in the next record
        if (eof()) end_error();
        line_ = lines_[rec_++];
        pos_ = 0;
        continue;
      }
      if (line_[pos_] != ',') break;
      if (comma_) throw Error("soil_veg_gen_parm: null value in " + name_);
      comma_ = true;
      ++pos_;
    }
    if (line_[pos_] == '/') {
      slash_ = true;
      return false;
    }
    v.clear();
    const char q = line_[pos_];
    if (q == '\'' || q == '"') {
      ++pos_;
      while (pos_ < line_.size()) {
        if (line_[pos_] == q) {
          if (pos_ + 1 < line_.size() && line_[pos_ + 1] == q) {
            v += q;
            pos_ += 2;
            continue;
          }
          ++pos_;
          break;
        }
        v += line_[pos_++];
      }
    } else {
      while (pos_ < line_.size() && line_[pos_] != ' ' && line_[pos_] != '\t' && line_[pos_] != ',' &&
             line_[pos_] != '/')
        v += line_[pos_++];
    }
    comma_ = false;
    return true;
  }
  void item(int* x) {
    std::string v;
    if (!next(v)) return;
    errno = 0;
    char* e = nullptr;
    const long r = std::strtol(v.c_str(), &e, 10);
    if (v.empty() || *e != '\0' || errno) throw Error("soil_veg_gen_parm: bad integer '" + v + "' in " + name_);
    *x = int(r);
  }
  void item(double* x) {
    std::string v;
    if (!next(v)) return;
    // gfortran's list-directed REAL(8) input: the decimal string correctly rounded, as strtod (no Fortran D
    // exponent or signed exponent without a letter appears in the tables; anything strtod does not take whole is
    // refused)
    errno = 0;
    char* e = nullptr;
    const double r = std::strtod(v.c_str(), &e);
    if (v.empty() || *e != '\0' || errno == ERANGE) throw Error("soil_veg_gen_parm: bad real '" + v + "' in " + name_);
    *x = r;
  }
  void item(std::string* x) {
    std::string v;
    if (!next(v)) return;
    *x = v;
  }
  std::string name_;
  std::vector<std::string> lines_;
  std::size_t rec_ = 0;
  std::string line_;
  std::size_t pos_ = 0;
  bool slash_ = false, comma_ = true;
};

// Fortran's character comparison (.eq.): the shorter operand padded with blanks
bool feq(const std::string& a, const std::string& b) {
  auto rt = [](std::string s) {
    while (!s.empty() && s.back() == ' ') s.pop_back();
    return s;
  };
  return rt(a) == rt(b);
}

}  // namespace

void soil_veg_gen_parm(const std::string& table_dir, const std::string& mminlu, const std::string& mminsl,
                       noah::NoahTables& t) {
  const std::string dir = table_dir.empty() ? std::string(".") : table_dir;
  std::string astring, lutype, sltype;
  int iindex = 0;
  // ---- VEGPARM.TBL (mpas_atmphys_lsm_noahinit.F:255-363)
  {
    Unit u(dir + "/VEGPARM.TBL", "VEGPARM.TBL");
    constexpr int loop_max = 10;  // :232
    int lumatch = 0;              // :265
    int loop_count = 0;           // :267
    bool end = !u.read_a(astring);  // :268 (end=2002)
    while (!end && lumatch == 0) {  // :269 find_lutype
      if (!u.read_list(&lutype)) {  // :270 (end=2002)
        end = true;
        break;
      }
      if (!u.read_list(&t.lucats, &iindex)) throw Error("soil_veg_gen_parm: end of file reading VEGPARM.TBL");
      if (feq(lutype, mminlu)) {  // :273
        lumatch = 1;
      } else {
        loop_count = loop_count + 1;  // :279
        for (;;) {                    // :282 find_vegetation_parameter_flag
          if (!u.read_a(astring)) {   // :283 (end=2002)
            end = true;
            break;
          }
          if (astring.substr(0, 21) == "Vegetation Parameters") {
            break;
          } else if (loop_count >= loop_max) {
            throw Error("mpas_atmphys_lsm_noahinit.F: too many loops in VEGPARM.TBL");  // :287
          }
        }
      }
    }
    if (!end) {
      if (noah::kMaxVegCats < t.lucats)  // :293-312 (size(shdtbl) ... < lucats)
        throw Error("table sizes too small for value of lucats in module_sf_noahdrv.f");
      if (feq(lutype, mminlu)) {  // :314
        for (int lc = 0; lc < t.lucats; ++lc)  // :315-320
          if (!u.read_list(&iindex, &t.shdtbl[lc], &t.nrotbl[lc], &t.rstbl[lc], &t.rgltbl[lc], &t.hstbl[lc],
                           &t.snuptbl[lc], &t.maxalb[lc], &t.laimintbl[lc], &t.laimaxtbl[lc], &t.emissmintbl[lc],
                           &t.emissmaxtbl[lc], &t.albedomintbl[lc], &t.albedomaxtbl[lc], &t.z0mintbl[lc],
                           &t.z0maxtbl[lc], &t.ztopvtbl[lc], &t.zbotvtbl[lc]))
            throw Error("soil_veg_gen_parm: end of file reading VEGPARM.TBL");
        auto rd = [&](auto* x) {
          if (!u.read_list(x)) throw Error("soil_veg_gen_parm: end of file reading VEGPARM.TBL");
        };
        u.skip(); rd(&t.topt_data);    // :322-323
        u.skip(); rd(&t.cmcmax_data);  // :324-325
        u.skip(); rd(&t.cfactr_data);  // :326-327
        u.skip(); rd(&t.rsmax_data);   // :328-329
        u.skip(); rd(&t.bare);         // :330-331
        u.skip(); rd(&t.natural);      // :332-333
        u.skip(); u.skip(); u.skip(); rd(&t.lcz[0]);  // :334-337
        for (int n = 1; n < 11; ++n) {  // :338-357
          u.skip();
          rd(&t.lcz[n]);
        }
      }
    }
    if (lumatch == 0)  // :360-363 (2002 continue)
      throw Error("land use dataset '" + mminlu + "' not found in VEGPARM.TBL.");
  }
  // ---- SOILPARM.TBL (:441-488)
  {
    Unit u(dir + "/SOILPARM.TBL", "SOILPARM.TBL");
    int lumatch = 0;  // :450
    u.skip();         // :451
    if (u.read_a4(sltype)) {  // :452 (format a4, end=2003)
      if (!u.read_list(&t.slcats, &iindex)) throw Error("soil_veg_gen_parm: end of file reading SOILPARM.TBL");
      if (feq(sltype, mminsl)) lumatch = 1;  // :455-460
      // :462-474: the size checks of MPAS do nothing (an empty IF); the arrays hold NSLTYPE = 30
      if (feq(sltype, mminsl)) {  // :475
        if (t.slcats > noah::kMaxSoilCats)
          throw Error("soil_veg_gen_parm: SOILPARM.TBL has more than NSLTYPE = 30 categories (MPAS writes out of "
                      "bounds)");
        for (int lc = 0; lc < t.slcats; ++lc)  // :476-480
          if (!u.read_list(&iindex, &t.bb[lc], &t.drysmc[lc], &t.f11[lc], &t.maxsmc[lc], &t.refsmc[lc],
                           &t.satpsi[lc], &t.satdk[lc], &t.satdw[lc], &t.wltsmc[lc], &t.qtz[lc]))
            throw Error("soil_veg_gen_parm: end of file reading SOILPARM.TBL");
      }
    }
    if (lumatch == 0)  // :482-488 (2003 continue)
      throw Error("soil texture in input file does not match soilparm table: inconsistent or missing soilparm file");
  }
  // ---- GENPARM.TBL (:524-571)
  {
    Unit u(dir + "/GENPARM.TBL", "GENPARM.TBL");
    auto rd = [&](auto* x) {
      if (!u.read_list(x)) throw Error("soil_veg_gen_parm: end of file reading GENPARM.TBL");
    };
    int num_slope = 0;
    u.skip();         // :529
    u.skip();         // :530
    rd(&num_slope);   // :531
    t.slpcats = num_slope;  // :533
    if (noah::kMaxSlopeCats < num_slope)  // :535-537
      throw Error("NUM_SLOPE too large for slope_data arrayin module_sf_noahdrv");
    for (int lc = 0; lc < t.slpcats; ++lc) rd(&t.slope_data[lc]);  // :539-541
    u.skip(); rd(&t.sbeta_data);   // :542-543
    u.skip(); rd(&t.fxexp_data);   // :544-545
    u.skip(); rd(&t.csoil_data);   // :546-547
    u.skip(); rd(&t.salp_data);    // :548-549
    u.skip(); rd(&t.refdk_data);   // :550-551
    u.skip(); rd(&t.refkdt_data);  // :552-553
    u.skip(); rd(&t.frzk_data);    // :554-555
    u.skip(); rd(&t.zbot_data);    // :556-557
    u.skip(); rd(&t.czil_data);    // :558-559
    u.skip(); rd(&t.smlow_data);   // :560-561
    u.skip(); rd(&t.smhigh_data);  // :562-563
    u.skip(); rd(&t.lvcoef_data);  // :564-565
  }
}

}  // namespace cppmpas::phys
