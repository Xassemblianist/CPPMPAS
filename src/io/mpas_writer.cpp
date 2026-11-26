// SPDX-License-Identifier: GPL-3.0-or-later
#include "io/mpas_writer.hpp"

#include <algorithm>
#include <cstdio>
#include <map>
#include <random>
#include <set>

#include "io/registry.hpp"

namespace cppmpas {

namespace {

constexpr std::size_t kStrLen = 64;  // ShortStrKIND (mpas_io_streams.F:1518-1519)

bool is_element_dim(const std::string& d) { return d == "nCells" || d == "nEdges" || d == "nVertices"; }
Element element_of(const std::string& d) {
  if (d == "nCells") return Element::Cell;
  if (d == "nEdges") return Element::Edge;
  return Element::Vertex;
}

std::size_t file_dim(const Mesh& m, const std::string& d) {
  if (d == "nCells") return m.nCells;
  if (d == "nEdges") return m.nEdges;
  if (d == "nVertices") return m.nVertices;
  if (d == "maxEdges") return m.maxEdges;
  if (d == "maxEdges2") return m.maxEdges2;
  if (d == "vertexDegree") return m.vertexDegree;
  if (d == "TWO") return 2;
  if (d == "R3") return 3;
  if (d == "FIFTEEN") return 15;
  if (d == "nVertLevels" && m.nVertLevels > 0) return m.nVertLevels;
  if (d == "nVertLevelsP1" && m.nVertLevels > 0) return m.nVertLevels + 1;
  if (d == "StrLen") return kStrLen;
  throw Error("dimension '" + d + "' is not handled by CPPMPAS");
}

std::vector<std::string> file_dims(const RegistryVar& rv) {
  std::vector<std::string> d;
  if (rv.has_time()) d.push_back("Time");
  auto sp = rv.space_dims();
  for (auto it = sp.rbegin(); it != sp.rend(); ++it) d.push_back(*it);
  if (rv.type == "text") d.push_back("StrLen");
  return d;
}

// Defines the dimensions of a field not yet defined, in the order MPAS defines them: Fortran order (StrLen first for
// character fields, mpas_io_streams.F:1518-1519, 1610-1611), the inner dimensions, the outer dimension, then Time
// (mpas_io_streams.F:1700-1745). dims is the C-order list of file_dims(). len(d) gives the length of d.
template <class Len>
void define_dims(NcFile& f, const std::vector<std::string>& dims, std::vector<std::string>& defined, Len len) {
  for (auto it = dims.rbegin(); it != dims.rend(); ++it) {
    if (std::find(defined.begin(), defined.end(), *it) != defined.end()) continue;
    f.def_dim(*it, len(*it));
    defined.push_back(*it);
  }
}

// gen_random (random_id.c:35-46) with idLength = 10 (mpas_stream_manager.F:4263, 4306-4310): ten characters of
// [0-9a-z]. MPAS seeds rand() with time(NULL) (random_id.c:21), so the id is not reproducible from run to run.
std::string random_file_id() {
  static const char alphanum[] = "0123456789abcdefghijklmnopqrstuvwxyz";
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<int> pick(0, static_cast<int>(sizeof(alphanum)) - 2);
  std::string id;
  for (int i = 0; i < 10; ++i) id += alphanum[pick(gen)];
  return id;
}

// prewrite_reindex (mpas_stream_manager.F:4894-5210): local connectivity -> global IDs, still in local element order.
std::vector<int> prewrite_reindex(const Mesh& m, const std::string& name, const FArray<int>& a) {
  const int inner = a.extent(0);
  const int outer = a.extent(1) - 1;
  std::vector<int> out(static_cast<std::size_t>(inner) * outer, 0);
  auto put = [&](int j, int i, int v) { out[(j - 1) + static_cast<std::size_t>(inner) * (i - 1)] = v; };
  if (name == "cellsOnCell" || name == "edgesOnCell" || name == "verticesOnCell") {
    const FArray<int>& id = name == "cellsOnCell" ? m.indexToCellID
                            : name == "edgesOnCell" ? m.indexToEdgeID
                                                    : m.indexToVertexID;
    for (int i = 1; i <= m.nCellsSolve(); ++i) {  // 5047-5054, 5068-5075, 5089-5096
      for (int j = 1; j <= m.nEdgesOnCell(i); ++j) put(j, i, id(a(j, i)));
      for (int j = m.nEdgesOnCell(i) + 1; j <= m.maxEdges; ++j) put(j, i, 0);  // UNUSED_CELL/EDGE/VERTEX = 0
    }
  } else if (name == "cellsOnEdge" || name == "verticesOnEdge") {
    const FArray<int>& id = name == "cellsOnEdge" ? m.indexToCellID : m.indexToVertexID;
    for (int i = 1; i <= m.nEdgesSolve(); ++i) {  // 5110-5113, 5128-5131
      put(1, i, id(a(1, i)));
      put(2, i, id(a(2, i)));
    }
  } else if (name == "edgesOnEdge") {
    for (int i = 1; i <= m.nEdgesSolve(); ++i) {  // 5146-5153
      for (int j = 1; j <= m.nEdgesOnEdge(i); ++j) put(j, i, m.indexToEdgeID(a(j, i)));
      for (int j = m.nEdgesOnEdge(i) + 1; j <= m.maxEdges2; ++j) put(j, i, 0);
    }
  } else if (name == "cellsOnVertex" || name == "edgesOnVertex") {
    const FArray<int>& id = name == "cellsOnVertex" ? m.indexToCellID : m.indexToEdgeID;
    for (int i = 1; i <= m.nVerticesSolve(); ++i)  // 5167-5171, 5185-5189
      for (int j = 1; j <= m.vertexDegree; ++j) put(j, i, id(a(j, i)));
  } else {
    for (int i = 1; i <= outer; ++i)
      for (int j = 1; j <= inner; ++j) put(j, i, a(j, i));
  }
  return out;
}

const std::set<std::string> kConnectivity = {"cellsOnCell", "edgesOnCell", "verticesOnCell", "cellsOnEdge",
                                             "verticesOnEdge", "edgesOnEdge", "cellsOnVertex", "edgesOnVertex"};

// The values of one field in file order (element records by global ID), without the garbage element.
template <class T> std::vector<T> file_values(const Mesh& m, const RegistryVar& rv, const FArray<T>& a,
                                              const std::vector<T>& local_values) {
  const auto sp = rv.space_dims();
  std::size_t inner = 1;
  for (std::size_t k = 0; k + 1 < sp.size(); ++k) inner *= file_dim(m, sp[k]);
  if (!sp.empty() && is_element_dim(sp.back())) {
    const ElementOrder& o = m.order(element_of(sp.back()));
    std::vector<T> out(inner * o.n());
    local_to_file(o, inner, local_values.data(), out.data());
    return out;
  }
  (void)a;
  std::size_t n = inner * (sp.empty() ? 1 : file_dim(m, sp.back()));
  return std::vector<T>(local_values.begin(), local_values.begin() + n);
}

void write_values(NcFile& f, MpasInput& in, const std::string& name) {
  const RegistryVar& rv = registry_var(name);
  const FieldRef ref = field_ref(in, name);
  const Mesh& m = in.mesh;
  const long record = rv.has_time() ? 0 : -1;
  if (ref.c) {
    if (ref.c->size() > kStrLen) throw Error("write " + name + ": value longer than StrLen = 64");
    std::vector<char> v(ref.c->begin(), ref.c->end());
    v.resize(kStrLen, ' ');
    f.put_char(name, v, record);
  } else if (ref.i) {
    if (!ref.i->allocated()) throw Error("write " + name + ": the field was neither read nor computed");
    std::vector<int> local = kConnectivity.count(name)
                                 ? prewrite_reindex(m, name, *ref.i)
                                 : std::vector<int>(ref.i->data(), ref.i->data() + ref.i->size());
    f.put_int(name, file_values(m, rv, *ref.i, local), record);
  } else if (ref.scalar_slot > 0) {
    const FArray<real>& s = *ref.r;
    std::vector<double> local(static_cast<std::size_t>(m.nVertLevels) * m.nCells);
    for (int i = 1; i <= m.nCells; ++i)
      for (int k = 1; k <= m.nVertLevels; ++k)
        local[static_cast<std::size_t>(i - 1) * m.nVertLevels + (k - 1)] = s(ref.scalar_slot, k, i);
    std::vector<double> out(local.size());
    local_to_file(m.cell_order, m.nVertLevels, local.data(), out.data());
    f.put_double(name, out, record);
  } else {
    if (!ref.r->allocated()) throw Error("write " + name + ": the field was neither read nor computed");
    std::vector<double> local(ref.r->data(), ref.r->data() + ref.r->size());
    f.put_double(name, file_values(m, rv, *ref.r, local), record);
  }
}

nc_type nc_type_of(const RegistryVar& rv) {
  return rv.type == "real" ? NC_DOUBLE : rv.type == "integer" ? NC_INT : NC_CHAR;
}

int cmode_of(int format) {
  switch (format) {
    case NC_FORMAT_CLASSIC: return 0;
    case NC_FORMAT_64BIT_OFFSET: return NC_64BIT_OFFSET;
    case NC_FORMAT_CDF5: return NC_64BIT_DATA;
  }
  throw Error("output format not handled");
}

std::string yes_no(bool b) { return b ? "YES" : "NO"; }

}  // namespace

void write_round_trip(MpasInput& in, const std::string& path, std::vector<std::string> fields) {
  if (fields.empty()) fields = in.source.vars_read;
  for (const auto& n : fields) {
    if (std::find(in.source.vars_read.begin(), in.source.vars_read.end(), n) == in.source.vars_read.end())
      throw Error("round trip: field " + n + " was not read from " + in.source.path);
  }
  NcFile f = NcFile::create(path, cmode_of(in.source.format));
  std::vector<std::string> defined;
  for (const auto& n : fields) {
    const RegistryVar& rv = registry_var(n);
    const auto dims = file_dims(rv);
    define_dims(f, dims, defined, [&](const std::string& d) {
      const auto it = in.source.dims.find(d);
      if (it == in.source.dims.end()) throw Error("round trip: dimension " + d + " not in the source");
      return it->second.second ? NC_UNLIMITED : it->second.first;
    });
    f.def_var(n, nc_type_of(rv), dims);
    for (const auto& a : in.source.var_atts[n]) f.put_att(n, a);
  }
  for (const auto& a : in.source.global_atts) f.put_att("", a);
  f.enddef();
  for (const auto& n : fields) write_values(f, in, n);
  f.close();
}

void write_history(MpasInput& in, const std::string& path, const std::vector<std::string>& fields,
                   const HistoryOptions& opt) {
  if (fields.empty()) throw Error("history: no fields to write");
  const Mesh& m = in.mesh;
  // every field is checked before the file is created: no partial file is left behind, and a file already at path
  // is not lost to a refused request
  // MPAS_io_put_att_text (mpas_io.F:6036-6045, 6208-6210): trim, cut at StrKIND, trim; SMIOL then writes
  // strlen of the C string (smiol.c:1731), which ends at the first NUL
  auto trimmed = [](std::string s) {
    while (!s.empty() && s.back() == ' ') s.pop_back();
    s = s.substr(0, 512);
    while (!s.empty() && s.back() == ' ') s.pop_back();
    return s.substr(0, s.find('\0'));
  };
  // mpas_atm_core.F:271-274: Time's units are 'seconds since ' // initial_time with '_' -> ' ' (written trimmed)
  std::string time_units;
  if (std::find(fields.begin(), fields.end(), "Time") != fields.end()) {
    std::string t = in.state.initial_time;
    while (!t.empty() && t.back() == ' ') t.pop_back();
    if (t.empty()) throw Error("write Time: initial_time was not read");
    std::replace(t.begin(), t.end(), '_', ' ');
    time_units = trimmed("seconds since " + t);
  }
  std::set<std::string> names;
  for (const auto& n : fields) {
    if (!names.insert(n).second) throw Error("history: field " + n + " is listed twice");
    const FieldRef ref = field_ref(in, n);  // a name the Registry does not define, or one without a port, is an error
    const RegistryVar& rv = registry_var(n);
    // a character field MPAS writes always holds a value (xtime from the clock, initial_time from the init file)
    const bool held = ref.c ? !ref.c->empty() : ref.i ? ref.i->allocated() : ref.r->allocated();
    if (!held) throw Error("write " + n + ": the field was neither read nor computed");
    if (ref.c && ref.c->size() > kStrLen) throw Error("write " + n + ": value longer than StrLen = 64");
    for (const auto& d : file_dims(rv))
      if (d != "Time") (void)file_dim(m, d);  // a dimension CPPMPAS does not handle is an error
    for (const auto& kv : rv.attrs)
      if (kv.first == "_FillValue") throw Error("history: numeric _FillValue attributes are not handled (" + n + ")");
  }
  if (opt.nprocs < 1 || opt.nprocs > 9) throw Error("history: nprocs outside 1..9 is not handled");
  NcFile f = NcFile::create(path, NC_64BIT_DATA);
  try {  // an error the checks above do not foresee (a NetCDF error, ...) removes the file
    std::vector<std::string> defined;
    for (const auto& n : fields) {
      const RegistryVar& rv = registry_var(n);
      const auto dims = file_dims(rv);
      define_dims(f, dims, defined,
                  [&](const std::string& d) { return d == "Time" ? NC_UNLIMITED : file_dim(m, d); });
      f.def_var(n, nc_type_of(rv), dims);
      // MPAS trims these too (mpas_io_streams.F:4341); no attribute of a field CPPMPAS writes ends in a blank
      for (const auto& [k, v] : rv.attrs) {
        // mpas_atm_core.F:271-274: Time's units are 'seconds since ' // initial_time with '_' -> ' ' (written trimmed)
        if (n == "Time" && k == "units") {
          f.put_att_text(n, k, time_units);
        } else {
          f.put_att_text(n, k, v);
        }
      }
    }
    // mpas_subdriver.F:497-530 (character values are written trimmed, mpas_io.F:6036, 6208-6210)
    f.put_att_text("", "model_name", "mpas");
    f.put_att_text("", "core_name", "atmosphere");
    f.put_att_text("", "version", "8.4.2");
    f.put_att_text("", "source", "MPAS");        // mpas_atm_core_interface.F:56
    f.put_att_text("", "Conventions", "MPAS");   // mpas_atm_core_interface.F:55
    f.put_att_text("", "git_version", trimmed(opt.git_version));
    f.put_att_text("", "on_a_sphere", yes_no(m.on_a_sphere));
    f.put_att_double("", "sphere_radius", m.sphere_radius);
    f.put_att_text("", "is_periodic", yes_no(m.is_periodic));
    f.put_att_double("", "x_period", m.x_period);
    f.put_att_double("", "y_period", m.y_period);
    f.put_att_text("", "history", "mpirun -n " + std::to_string(opt.nprocs) + " ./atmosphere_model");
    f.put_att_text("", "parent_id", trimmed(m.parent_id));
    f.put_att_text("", "mesh_spec", trimmed(m.mesh_spec));
    f.put_att_text("", "file_id", trimmed(opt.file_id.empty() ? random_file_id() : opt.file_id));
    f.enddef();
    for (const auto& n : fields) write_values(f, in, n);
    f.close();
  } catch (...) {
    std::remove(path.c_str());
    throw;
  }
}

}  // namespace cppmpas
