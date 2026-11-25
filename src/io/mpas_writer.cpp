// SPDX-License-Identifier: GPL-3.0-or-later
#include "io/mpas_writer.hpp"

#include <algorithm>
#include <map>
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
    for (const auto& d : dims) {
      if (std::find(defined.begin(), defined.end(), d) != defined.end()) continue;
      const auto it = in.source.dims.find(d);
      if (it == in.source.dims.end()) throw Error("round trip: dimension " + d + " not in the source");
      f.def_dim(d, it->second.second ? NC_UNLIMITED : it->second.first);
      defined.push_back(d);
    }
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
  NcFile f = NcFile::create(path, NC_64BIT_DATA);
  std::vector<std::string> defined;
  for (const auto& n : fields) {
    const RegistryVar& rv = registry_var(n);
    (void)field_ref(in, n);  // a field without a port is an error here, before anything is written
    const auto dims = file_dims(rv);
    for (const auto& d : dims) {
      if (std::find(defined.begin(), defined.end(), d) != defined.end()) continue;
      f.def_dim(d, d == "Time" ? NC_UNLIMITED : file_dim(m, d));
      defined.push_back(d);
    }
    f.def_var(n, nc_type_of(rv), dims);
    for (const auto& [k, v] : rv.attrs) {
      if (k == "_FillValue") throw Error("history: numeric _FillValue attributes are not handled (" + n + ")");
      f.put_att_text(n, k, v);
    }
  }
  // mpas_subdriver.F:497-530 (character values are written trimmed, mpas_io.F:6036, 6208-6210)
  f.put_att_text("", "model_name", "mpas");
  f.put_att_text("", "core_name", "atmosphere");
  f.put_att_text("", "version", "8.4.2");
  f.put_att_text("", "source", "MPAS");        // mpas_atm_core_interface.F:56
  f.put_att_text("", "Conventions", "MPAS");   // mpas_atm_core_interface.F:55
  f.put_att_text("", "git_version", opt.git_version);
  f.put_att_text("", "on_a_sphere", yes_no(m.on_a_sphere));
  f.put_att_double("", "sphere_radius", m.sphere_radius);
  f.put_att_text("", "is_periodic", yes_no(m.is_periodic));
  f.put_att_double("", "x_period", m.x_period);
  f.put_att_double("", "y_period", m.y_period);
  if (opt.nprocs < 1 || opt.nprocs > 9) throw Error("history: nprocs outside 1..9 is not handled");
  f.put_att_text("", "history", "mpirun -n " + std::to_string(opt.nprocs) + " ./atmosphere_model");
  auto trimmed = [](std::string s) {
    while (!s.empty() && (s.back() == ' ' || s.back() == '\0')) s.pop_back();
    return s;
  };
  f.put_att_text("", "parent_id", trimmed(m.parent_id));
  f.put_att_text("", "mesh_spec", trimmed(m.mesh_spec));
  f.put_att_text("", "file_id", opt.file_id);
  f.enddef();
  for (const auto& n : fields) write_values(f, in, n);
  f.close();
}

}  // namespace cppmpas
