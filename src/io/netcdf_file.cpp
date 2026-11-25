// SPDX-License-Identifier: GPL-3.0-or-later
#include "io/netcdf_file.hpp"

namespace cppmpas {

void NcFile::check(int status, const std::string& what) const {
  if (status != NC_NOERR) throw Error(path_ + ": " + what + ": " + nc_strerror(status));
}

NcFile NcFile::open_read(const std::string& path) {
  int id = -1;
  const int st = nc_open(path.c_str(), NC_NOWRITE, &id);
  if (st != NC_NOERR) throw Error(path + ": cannot open: " + nc_strerror(st));
  return NcFile(id, path);
}

NcFile NcFile::create(const std::string& path, int cmode) {
  int id = -1;
  const int st = nc_create(path.c_str(), NC_CLOBBER | cmode, &id);
  if (st != NC_NOERR) throw Error(path + ": cannot create: " + nc_strerror(st));
  // values are always written in full: no fill pass (the bytes written are the bytes of the data)
  int old = 0;
  nc_set_fill(id, NC_NOFILL, &old);
  return NcFile(id, path);
}

NcFile::NcFile(NcFile&& o) noexcept : ncid_(o.ncid_), path_(std::move(o.path_)) { o.ncid_ = -1; }
NcFile::~NcFile() {
  if (ncid_ >= 0) nc_close(ncid_);
}
void NcFile::close() {
  if (ncid_ >= 0) {
    const int id = ncid_;
    ncid_ = -1;
    check(nc_close(id), "close");
  }
}

int NcFile::format() const {
  int f = 0;
  check(nc_inq_format(ncid_, &f), "inq_format");
  return f;
}

const char* NcFile::format_name(int f) {
  switch (f) {
    case NC_FORMAT_CLASSIC: return "classic (CDF-1)";
    case NC_FORMAT_64BIT_OFFSET: return "64-bit offset (CDF-2)";
    case NC_FORMAT_CDF5: return "64-bit data (CDF-5)";
    case NC_FORMAT_NETCDF4: return "netCDF-4/HDF5";
    case NC_FORMAT_NETCDF4_CLASSIC: return "netCDF-4 classic model";
  }
  return "unknown";
}

bool NcFile::has_dim(const std::string& name) const {
  int id = -1;
  return nc_inq_dimid(ncid_, name.c_str(), &id) == NC_NOERR;
}

std::size_t NcFile::dim_len(const std::string& name) const {
  int id = -1;
  check(nc_inq_dimid(ncid_, name.c_str(), &id), "dimension " + name + " not found");
  std::size_t len = 0;
  check(nc_inq_dimlen(ncid_, id, &len), "dimension " + name);
  return len;
}

bool NcFile::dim_unlimited(const std::string& name) const {
  int id = -1;
  check(nc_inq_dimid(ncid_, name.c_str(), &id), "dimension " + name + " not found");
  int n = 0;
  check(nc_inq_unlimdims(ncid_, &n, nullptr), "inq_unlimdims");
  std::vector<int> ids(n);
  if (n) check(nc_inq_unlimdims(ncid_, &n, ids.data()), "inq_unlimdims");
  for (int u : ids)
    if (u == id) return true;
  return false;
}

bool NcFile::has_var(const std::string& name) const {
  int id = -1;
  return nc_inq_varid(ncid_, name.c_str(), &id) == NC_NOERR;
}

int NcFile::varid(const std::string& name) const {
  if (name.empty()) return NC_GLOBAL;
  int id = -1;
  check(nc_inq_varid(ncid_, name.c_str(), &id), "variable " + name + " not found");
  return id;
}

VarInfo NcFile::var(const std::string& name) const {
  VarInfo v;
  v.id = varid(name);
  int ndims = 0;
  check(nc_inq_var(ncid_, v.id, nullptr, &v.type, &ndims, nullptr, nullptr), "variable " + name);
  std::vector<int> dimids(ndims);
  if (ndims) check(nc_inq_vardimid(ncid_, v.id, dimids.data()), "variable " + name);
  for (int d : dimids) {
    char buf[NC_MAX_NAME + 1];
    std::size_t len = 0;
    check(nc_inq_dim(ncid_, d, buf, &len), "variable " + name);
    v.dims.emplace_back(buf);
    v.shape.push_back(len);
  }
  return v;
}

std::vector<std::string> NcFile::var_names() const {
  int n = 0;
  check(nc_inq_nvars(ncid_, &n), "inq_nvars");
  std::vector<std::string> out;
  for (int i = 0; i < n; ++i) {
    char buf[NC_MAX_NAME + 1];
    check(nc_inq_varname(ncid_, i, buf), "inq_varname");
    out.emplace_back(buf);
  }
  return out;
}

bool NcFile::has_att(const std::string& varname, const std::string& att) const {
  int vid = NC_GLOBAL;
  if (!varname.empty() && nc_inq_varid(ncid_, varname.c_str(), &vid) != NC_NOERR) return false;
  int aid = -1;
  return nc_inq_attid(ncid_, vid, att.c_str(), &aid) == NC_NOERR;
}

RawAttribute NcFile::att(const std::string& varname, const std::string& name) const {
  const int vid = varid(varname);
  RawAttribute a;
  a.name = name;
  const std::string where = (varname.empty() ? "global" : varname) + " attribute " + name;
  check(nc_inq_att(ncid_, vid, name.c_str(), &a.type, &a.len), where);
  std::size_t size = 0;
  check(nc_inq_type(ncid_, a.type, nullptr, &size), where);
  if (a.type == NC_STRING || a.type > NC_MAX_ATOMIC_TYPE) throw Error(path_ + ": " + where + ": type not handled");
  a.bytes.resize(a.len * size);
  if (a.len) check(nc_get_att(ncid_, vid, name.c_str(), a.bytes.data()), where);
  return a;
}

std::vector<RawAttribute> NcFile::atts(const std::string& varname) const {
  const int vid = varid(varname);
  int n = 0;
  if (vid == NC_GLOBAL)
    check(nc_inq_natts(ncid_, &n), "inq_natts");
  else
    check(nc_inq_varnatts(ncid_, vid, &n), "inq_varnatts " + varname);
  std::vector<RawAttribute> out;
  for (int i = 0; i < n; ++i) {
    char buf[NC_MAX_NAME + 1];
    check(nc_inq_attname(ncid_, vid, i, buf), "inq_attname");
    out.push_back(att(varname, buf));
  }
  return out;
}

void NcFile::record_range(const std::string& name, long record, std::vector<std::size_t>& start,
                          std::vector<std::size_t>& count) const {
  const VarInfo v = var(name);
  start.assign(v.dims.size(), 0);
  count = v.shape;
  if (record >= 0) {
    if (v.dims.empty() || !dim_unlimited(v.dims[0])) throw Error(path_ + ": " + name + " has no record dimension");
    start[0] = static_cast<std::size_t>(record);
    count[0] = 1;
  }
}

namespace {
std::size_t product(const std::vector<std::size_t>& c) {
  std::size_t n = 1;
  for (auto x : c) n *= x;
  return n;
}
}  // namespace

std::vector<double> NcFile::get_double(const std::string& name, long record) const {
  std::vector<std::size_t> s, c;
  record_range(name, record, s, c);
  std::vector<double> out(product(c));
  if (!out.empty()) check(nc_get_vara_double(ncid_, varid(name), s.data(), c.data(), out.data()), "read " + name);
  return out;
}
std::vector<int> NcFile::get_int(const std::string& name, long record) const {
  std::vector<std::size_t> s, c;
  record_range(name, record, s, c);
  std::vector<int> out(product(c));
  if (!out.empty()) check(nc_get_vara_int(ncid_, varid(name), s.data(), c.data(), out.data()), "read " + name);
  return out;
}
std::vector<char> NcFile::get_char(const std::string& name, long record) const {
  std::vector<std::size_t> s, c;
  record_range(name, record, s, c);
  std::vector<char> out(product(c));
  if (!out.empty()) check(nc_get_vara_text(ncid_, varid(name), s.data(), c.data(), out.data()), "read " + name);
  return out;
}

void NcFile::def_dim(const std::string& name, std::size_t len) {
  int id = -1;
  check(nc_def_dim(ncid_, name.c_str(), len, &id), "def_dim " + name);
}

int NcFile::def_var(const std::string& name, nc_type type, const std::vector<std::string>& dims) {
  std::vector<int> ids;
  for (const auto& d : dims) {
    int id = -1;
    check(nc_inq_dimid(ncid_, d.c_str(), &id), "def_var " + name + ": dimension " + d);
    ids.push_back(id);
  }
  int vid = -1;
  check(nc_def_var(ncid_, name.c_str(), type, static_cast<int>(ids.size()), ids.data(), &vid), "def_var " + name);
  return vid;
}

void NcFile::put_att(const std::string& varname, const RawAttribute& a) {
  check(nc_put_att(ncid_, varid(varname), a.name.c_str(), a.type, a.len, a.bytes.data()),
        "put_att " + varname + ":" + a.name);
}
void NcFile::put_att_text(const std::string& varname, const std::string& name, const std::string& value) {
  check(nc_put_att_text(ncid_, varid(varname), name.c_str(), value.size(), value.data()),
        "put_att " + varname + ":" + name);
}
void NcFile::put_att_double(const std::string& varname, const std::string& name, double value) {
  check(nc_put_att_double(ncid_, varid(varname), name.c_str(), NC_DOUBLE, 1, &value),
        "put_att " + varname + ":" + name);
}
void NcFile::put_att_int(const std::string& varname, const std::string& name, int value) {
  check(nc_put_att_int(ncid_, varid(varname), name.c_str(), NC_INT, 1, &value), "put_att " + varname + ":" + name);
}
void NcFile::enddef() { check(nc_enddef(ncid_), "enddef"); }

void NcFile::put_double(const std::string& name, const std::vector<double>& v, long record) {
  std::vector<std::size_t> s, c;
  record_range(name, -1, s, c);
  if (record >= 0) {
    s[0] = static_cast<std::size_t>(record);
    c[0] = 1;
  }
  if (product(c) != v.size()) throw Error(path_ + ": write " + name + ": size mismatch");
  check(nc_put_vara_double(ncid_, varid(name), s.data(), c.data(), v.data()), "write " + name);
}
void NcFile::put_int(const std::string& name, const std::vector<int>& v, long record) {
  std::vector<std::size_t> s, c;
  record_range(name, -1, s, c);
  if (record >= 0) {
    s[0] = static_cast<std::size_t>(record);
    c[0] = 1;
  }
  if (product(c) != v.size()) throw Error(path_ + ": write " + name + ": size mismatch");
  check(nc_put_vara_int(ncid_, varid(name), s.data(), c.data(), v.data()), "write " + name);
}
void NcFile::put_char(const std::string& name, const std::vector<char>& v, long record) {
  std::vector<std::size_t> s, c;
  record_range(name, -1, s, c);
  if (record >= 0) {
    s[0] = static_cast<std::size_t>(record);
    c[0] = 1;
  }
  if (product(c) != v.size()) throw Error(path_ + ": write " + name + ": size mismatch");
  check(nc_put_vara_text(ncid_, varid(name), s.data(), c.data(), v.data()), "write " + name);
}

}  // namespace cppmpas
