// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// Thin RAII wrapper over the NetCDF-C API: every failing call throws cppmpas::Error naming the file and variable.
#include <netcdf.h>

#include <cstddef>
#include <string>
#include <vector>

#include "core/error.hpp"

namespace cppmpas {

// An attribute as stored in a file: type, element count and raw bytes, copied without conversion.
struct RawAttribute {
  std::string name;
  nc_type type = NC_NAT;
  std::size_t len = 0;
  std::vector<unsigned char> bytes;
  std::string text() const { return std::string(bytes.begin(), bytes.end()); }  // for NC_CHAR
};

struct VarInfo {
  int id = -1;
  nc_type type = NC_NAT;
  std::vector<std::string> dims;     // C order (slowest first), as in the file
  std::vector<std::size_t> shape;
};

class NcFile {
 public:
  static NcFile open_read(const std::string& path);
  // cmode: NC_CLOBBER | NC_64BIT_DATA (CDF-5), NC_64BIT_OFFSET (CDF-2) or 0 (classic)
  static NcFile create(const std::string& path, int cmode);
  NcFile(NcFile&& o) noexcept;
  NcFile& operator=(NcFile&&) = delete;
  NcFile(const NcFile&) = delete;
  ~NcFile();
  void close();

  const std::string& path() const { return path_; }
  int format() const;  // NC_FORMAT_CLASSIC, NC_FORMAT_64BIT_OFFSET, NC_FORMAT_CDF5, NC_FORMAT_NETCDF4, ...
  static const char* format_name(int format);

  bool has_dim(const std::string& name) const;
  std::size_t dim_len(const std::string& name) const;
  bool dim_unlimited(const std::string& name) const;
  bool has_var(const std::string& name) const;
  VarInfo var(const std::string& name) const;
  std::vector<std::string> var_names() const;

  // attributes; varname "" is the global attribute set
  bool has_att(const std::string& varname, const std::string& att) const;
  RawAttribute att(const std::string& varname, const std::string& att) const;
  std::vector<RawAttribute> atts(const std::string& varname) const;

  // whole variable, or one record (index 0 of the first, unlimited dimension) when record >= 0
  std::vector<double> get_double(const std::string& name, long record = -1) const;
  std::vector<int> get_int(const std::string& name, long record = -1) const;
  std::vector<char> get_char(const std::string& name, long record = -1) const;

  // definitions (define mode) and data
  void def_dim(const std::string& name, std::size_t len);  // len = NC_UNLIMITED for the record dimension
  int def_var(const std::string& name, nc_type type, const std::vector<std::string>& dims);
  void put_att(const std::string& varname, const RawAttribute& a);
  void put_att_text(const std::string& varname, const std::string& name, const std::string& value);
  void put_att_double(const std::string& varname, const std::string& name, double value);
  void put_att_int(const std::string& varname, const std::string& name, int value);
  void enddef();
  void put_double(const std::string& name, const std::vector<double>& v, long record = -1);
  void put_int(const std::string& name, const std::vector<int>& v, long record = -1);
  void put_char(const std::string& name, const std::vector<char>& v, long record = -1);

 private:
  NcFile(int ncid, std::string path) : ncid_(ncid), path_(std::move(path)) {}
  int varid(const std::string& name) const;
  void check(int status, const std::string& what) const;
  void record_range(const std::string& name, long record, std::vector<std::size_t>& start,
                    std::vector<std::size_t>& count) const;
  int ncid_ = -1;
  std::string path_;
};

}  // namespace cppmpas
