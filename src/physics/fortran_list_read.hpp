// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// gfortran's formatted sequential input of the physics tables (VEGPARM.TBL, SOILPARM.TBL, GENPARM.TBL, LANDUSE.TBL):
// list-directed READs and the A edit descriptor, as the table readers of mpas_atmphys_lsm_noahinit.F and
// mpas_atmphys_landuse.F use them.
#include <cerrno>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>

#include "core/error.hpp"

namespace cppmpas::phys {

// A formatted sequential file as gfortran reads it: records are lines; a list-directed READ (fmt = *) starts at
// the next record, takes its items from as many records as it needs and skips the rest of its last record.
// Value separators: a comma or a slash with optional blanks around it, or blanks. A slash ends the READ (the
// remaining items keep their values). An undelimited character value ends at a blank, comma or slash; a delimited
// one is between apostrophes or quotes. A null value (two commas) is refused: no Noah table has one and MPAS's
// READs would leave the item unchanged.
class FortranListFile {
 public:
  // who: the reading routine, for the messages ("soil_veg_gen_parm", "landuse_init_forMPAS")
  FortranListFile(const std::string& path, const std::string& what, const std::string& who = "soil_veg_gen_parm")
      : name_(what), who_(who) {
    std::ifstream f(path);
    if (!f) {
      if (who == "soil_veg_gen_parm") throw Error("mpas_atmphys_lsm_noahinit.F: soil_veg_gen_parm: failure opening " + what);
      throw Error(who + ": failure opening " + what);
    }
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
  // READ (u, '(An)') s: the first n characters of the record (shorter records: what there is)
  bool read_an(std::string& s, std::size_t n) {
    if (eof()) return false;
    s = lines_[rec_++].substr(0, n);
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
    if (eof()) throw Error(who_ + ": end of file reading " + name_);
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
      if (comma_) throw Error(who_ + ": null value in " + name_);
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
    if (v.empty() || *e != '\0' || errno) throw Error(who_ + ": bad integer '" + v + "' in " + name_);
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
    if (v.empty() || *e != '\0' || errno == ERANGE) throw Error(who_ + ": bad real '" + v + "' in " + name_);
    *x = r;
  }
  void item(std::string* x) {
    std::string v;
    if (!next(v)) return;
    *x = v;
  }
  std::string name_, who_;
  std::vector<std::string> lines_;
  std::size_t rec_ = 0;
  std::string line_;
  std::size_t pos_ = 0;
  bool slash_ = false, comma_ = true;
};

// Fortran's character comparison (.eq.): the shorter operand padded with blanks
inline bool feq(const std::string& a, const std::string& b) {
  auto rt = [](std::string s) {
    while (!s.empty() && s.back() == ' ') s.pop_back();
    return s;
  };
  return rt(a) == rt(b);
}

}  // namespace cppmpas::phys
