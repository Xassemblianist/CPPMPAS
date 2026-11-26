// SPDX-License-Identifier: GPL-3.0-or-later
#include "io/mpas_reader.hpp"

#include <algorithm>
#include <cstring>
#include <limits>
#include <set>

#include "core/decomp.hpp"
#include "io/registry.hpp"

namespace cppmpas {

namespace {

// ShortStrKIND: the length of the StrLen dimension of character fields (mpas_io_streams.F:1518-1519, 1610-1611)
constexpr std::size_t kStrLen = 64;
// StrKIND (mpas_kind_types.F:32): the length of the character variables MPAS reads the global attributes into
// (mpas_bootstrapping.F:103) and of domain % parent_id / mesh_spec (mpas_domain_types.inc:32-33)
constexpr std::size_t kStrKIND = 512;

const std::vector<std::string> kMeshRequired = {
    "latCell", "lonCell", "xCell", "yCell", "zCell", "indexToCellID", "latEdge", "lonEdge", "xEdge", "yEdge",
    "zEdge", "indexToEdgeID", "latVertex", "lonVertex", "xVertex", "yVertex", "zVertex", "indexToVertexID",
    "cellsOnEdge", "nEdgesOnCell", "nEdgesOnEdge", "edgesOnCell", "edgesOnEdge", "weightsOnEdge", "dvEdge",
    "dcEdge", "angleEdge", "areaCell", "areaTriangle", "cellsOnCell", "verticesOnCell", "verticesOnEdge",
    "edgesOnVertex", "cellsOnVertex", "kiteAreasOnVertex", "meshDensity"};
// in the invariant stream (Registry.xml:499-566) but not always in a grid file (the x1.40962 grid has no fEdge /
// fVertex: init_atmosphere computes them); required in an init file, read from a grid file when present
const std::vector<std::string> kMeshInitOnly = {
    "fEdge", "fVertex", "edgeNormalVectors", "localVerticalUnitVectors", "cellTangentPlane", "nominalMinDc", "bdyMaskCell",
    "bdyMaskEdge", "bdyMaskVertex", "cf1", "cf2", "cf3", "zgrid", "rdzw", "dzu", "rdzu", "fzm", "fzp", "zxu", "zz",
    "zb", "zb3", "dss", "deriv_two", "coeffs_reconstruct"};
const std::vector<std::string> kState = {"initial_time", "xtime", "u_init", "v_init", "t_init", "u", "w", "rho",
                                         "theta", "rho_base", "theta_base", "surface_pressure", "Time"};
// Variables of the atmosphere input stream that MPAS reads but the dycore does not use in the verified cases
// (physics suite none, SC with Kessler microphysics only): the DO_PHYSICS fields of the invariant and
// input streams (Registry.xml:567-597, 625-646) and the diagnostic "relhum". (Time, Registry.xml:615, is read:
// it is time level 1 of the state's CF time, mpas_atm_core.F:137.)
const std::set<std::string> kInputNotUsed = {
    "relhum", "isltyp", "ivgtyp", "mminlu", "isice_lu", "iswater_lu", "landmask", "shdmin", "shdmax",
    "snoalb", "albedo12m", "greenfrac", "var2d", "con", "oa1", "oa2", "oa3", "oa4", "ol1", "ol2", "ol3", "ol4",
    "soilcomp", "soilcl1", "soilcl2", "soilcl3", "soilcl4", "sfc_albbck", "skintemp", "snow", "snowc", "snowh",
    "sst", "tmn", "vegfra", "seaice", "xice", "xland", "u10", "v10", "q2", "t2m", "precipw", "dzs", "sh2o", "smois",
    "tslb", "h_oml_initial"};

bool is_element_dim(const std::string& d) { return d == "nCells" || d == "nEdges" || d == "nVertices"; }
Element element_of(const std::string& d) {
  if (d == "nCells") return Element::Cell;
  if (d == "nEdges") return Element::Edge;
  return Element::Vertex;
}

// Dimension length in the file (garbage == false) or in memory (garbage == true: +1 on element dimensions).
std::size_t dim_size(const Mesh& m, const std::string& d, bool garbage) {
  const int g = garbage ? 1 : 0;
  if (d == "nCells") return m.nCells + g;
  if (d == "nEdges") return m.nEdges + g;
  if (d == "nVertices") return m.nVertices + g;
  if (d == "maxEdges") return m.maxEdges;
  if (d == "maxEdges2") return m.maxEdges2;
  if (d == "vertexDegree") return m.vertexDegree;
  if (d == "TWO") return 2;
  if (d == "R3") return 3;
  if (d == "FIFTEEN") return 15;
  if (d == "nVertLevels" && m.nVertLevels > 0) return m.nVertLevels;
  if (d == "nVertLevelsP1" && m.nVertLevels > 0) return m.nVertLevels + 1;
  throw Error("dimension '" + d + "' is not handled by CPPMPAS");
}

std::string join(const std::vector<std::string>& v) {
  std::string s;
  for (const auto& x : v) s += (s.empty() ? "" : ", ") + x;
  return s;
}

// The file layout MPAS defines for a Registry variable: Time first, then the Fortran dimensions reversed, then StrLen
// for character variables.
std::vector<std::string> file_dims(const RegistryVar& rv) {
  std::vector<std::string> d;
  if (rv.has_time()) d.push_back("Time");
  auto sp = rv.space_dims();
  for (auto it = sp.rbegin(); it != sp.rend(); ++it) d.push_back(*it);
  if (rv.type == "text") d.push_back("StrLen");
  return d;
}

// Checks name, type and dimensions of a variable against the Registry and returns the record to read (-1: none).
long check_var(const NcFile& f, const Mesh& m, const RegistryVar& rv, const std::string& fileName) {
  const VarInfo vi = f.var(fileName);
  const nc_type want = rv.type == "real" ? NC_DOUBLE : rv.type == "integer" ? NC_INT : NC_CHAR;
  if (vi.type != want) {
    throw Error(f.path() + ": variable " + fileName + " has NetCDF type " + std::to_string(vi.type) +
                ", MPAS-Atmosphere (" + rv.type + ", PRECISION=double) needs " +
                (want == NC_DOUBLE ? "double" : want == NC_INT ? "int" : "char") +
                " (conversion on read is not handled)");
  }
  const std::vector<std::string> want_dims = file_dims(rv);
  if (vi.dims != want_dims)
    throw Error(f.path() + ": variable " + fileName + " has dimensions (" + join(vi.dims) + "), MPAS expects (" +
                join(want_dims) + ")");
  long record = -1;
  for (std::size_t k = 0; k < vi.dims.size(); ++k) {
    const std::string& d = vi.dims[k];
    if (d == "Time") {
      if (vi.shape[k] != 1)
        throw Error(f.path() + ": variable " + fileName + " has " + std::to_string(vi.shape[k]) +
                    " time records; reading an init file with more than one record is not handled");
      record = 0;
    } else if (d == "StrLen") {
      if (vi.shape[k] != kStrLen)
        throw Error(f.path() + ": StrLen = " + std::to_string(vi.shape[k]) + ", MPAS writes 64");
    } else if (vi.shape[k] != dim_size(m, d, false)) {
      throw Error(f.path() + ": variable " + fileName + ": dimension " + d + " = " + std::to_string(vi.shape[k]) +
                  ", expected " + std::to_string(dim_size(m, d, false)));
    }
  }
  return record;
}

template <class T> T default_of(const RegistryVar& rv) {
  if (rv.default_value.empty()) return T{0};  // src/tools/registry/gen_inc.c:615-630
  if constexpr (std::is_same_v<T, int>)
    return std::stoi(rv.default_value);
  else
    return static_cast<T>(std::stod(rv.default_value));
}

// Allocates a with MPAS's memory extents (garbage element included), fills it with the Registry default
// (mpas_block_creator.F:1649, 1680, ...) and reads the file values into it in MPAS local order.
template <class T>
void read_into(const NcFile& f, const Mesh& m, const RegistryVar& rv, FArray<T>& a) {
  const long record = check_var(f, m, rv, rv.name);
  const auto sp = rv.space_dims();
  std::vector<int> ext;
  for (const auto& d : sp) ext.push_back(static_cast<int>(dim_size(m, d, true)));
  if (ext.empty()) ext.push_back(1);
  a.allocate(ext, default_of<T>(rv));
  std::vector<T> buf;
  if constexpr (std::is_same_v<T, int>)
    buf = f.get_int(rv.name, record);
  else
    buf = f.get_double(rv.name, record);
  if (!sp.empty() && is_element_dim(sp.back())) {
    std::size_t inner = 1;
    for (std::size_t k = 0; k + 1 < sp.size(); ++k) inner *= dim_size(m, sp[k], false);
    file_to_local(m.order(element_of(sp.back())), inner, buf.data(), a.data());
  } else {
    std::copy(buf.begin(), buf.end(), a.data());
  }
}

void read_text(const NcFile& f, const Mesh& m, const RegistryVar& rv, std::string& s) {
  const long record = check_var(f, m, rv, rv.name);
  const std::vector<char> buf = f.get_char(rv.name, record);
  s.assign(buf.begin(), buf.end());
}

std::vector<int> raw_int(const NcFile& f, const Mesh& m, const std::string& name) {
  const RegistryVar& rv = registry_var(name);
  const long record = check_var(f, m, rv, name);
  return f.get_int(name, record);
}

std::string text_att(const NcFile& f, const std::string& name) {
  const RawAttribute a = f.att("", name);
  if (a.type != NC_CHAR) throw Error(f.path() + ": global attribute " + name + " is not text");
  if (a.len > kStrKIND)
    throw Error(f.path() + ": global attribute " + name + " has " + std::to_string(a.len) +
                " characters; MPAS reads it into a character(len=StrKIND = 512) variable (mpas_bootstrapping.F:103): "
                "longer values are not handled");
  return a.text();
}

double real_att(const NcFile& f, const std::string& name) {
  const RawAttribute a = f.att("", name);
  if (a.len != 1) throw Error(f.path() + ": global attribute " + name + " is not a scalar");
  if (a.type == NC_DOUBLE) {
    double v;
    std::memcpy(&v, a.bytes.data(), sizeof v);
    return v;
  }
  if (a.type == NC_FLOAT) {
    float v;
    std::memcpy(&v, a.bytes.data(), sizeof v);
    return v;
  }
  throw Error(f.path() + ": global attribute " + name + " is not a real number");
}

// Fortran's trim: trailing blanks only (SMIOLf_inquire_att copies every byte of the attribute, NULs included,
// smiolf.F90:1630-1631)
std::string rtrim(std::string s) {
  while (!s.empty() && s.back() == ' ') s.pop_back();
  return s;
}

// A mesh dimension: MPAS holds it as a default (32-bit) integer and allocates one garbage element more
// (mpas_domain_routines.F:251-265), and an empty mesh is not a mesh.
int mesh_dim(const NcFile& f, const char* d) {
  const std::size_t n = f.dim_len(d);
  const std::size_t max = static_cast<std::size_t>(std::numeric_limits<int>::max()) - 1;
  if (n < 1 || n > max)
    throw Error(f.path() + ": dimension " + d + " = " + std::to_string(n) + " outside 1.." + std::to_string(max));
  return static_cast<int>(n);
}

void read_dims_and_attributes(const NcFile& f, MpasInput& in) {
  Mesh& m = in.mesh;
  const int fmt = f.format();
  // The reference reads through SMIOL (PnetCDF): CDF-1, CDF-2 and CDF-5 files only.
  if (fmt != NC_FORMAT_CLASSIC && fmt != NC_FORMAT_64BIT_OFFSET && fmt != NC_FORMAT_CDF5)
    throw Error(f.path() + ": NetCDF format " + NcFile::format_name(fmt) +
                " is not handled (the reference MPAS build reads through SMIOL/PnetCDF: CDF-1, CDF-2, CDF-5)");
  in.source.path = f.path();
  in.source.format = fmt;
  for (const char* d : {"nCells", "nEdges", "nVertices", "maxEdges", "maxEdges2", "vertexDegree"})
    if (!f.has_dim(d)) throw Error(f.path() + ": required dimension " + std::string(d) + " is missing");
  m.nCells = mesh_dim(f, "nCells");
  m.nEdges = mesh_dim(f, "nEdges");
  m.nVertices = mesh_dim(f, "nVertices");
  m.maxEdges = mesh_dim(f, "maxEdges");
  m.maxEdges2 = mesh_dim(f, "maxEdges2");
  m.vertexDegree = mesh_dim(f, "vertexDegree");
  if (f.has_dim("TWO") && f.dim_len("TWO") != 2) throw Error(f.path() + ": dimension TWO is not 2");

  // mpas_bootstrapping.F:287-299
  if (!f.has_att("", "on_a_sphere")) {
    in.warnings.push_back("Attribute on_a_sphere not found; setting on_a_sphere to 'YES' (mpas_bootstrapping.F:289-291)");
    m.on_a_sphere = true;
  } else {
    m.on_a_sphere = text_att(f, "on_a_sphere").find("YES") != std::string::npos;
  }
  // mpas_bootstrapping.F:301-308
  if (!f.has_att("", "sphere_radius")) {
    in.warnings.push_back("Attribute sphere_radius not found; setting sphere_radius to 1.0 (mpas_bootstrapping.F:303-305)");
    m.sphere_radius = 1.0;
  } else {
    m.sphere_radius = real_att(f, "sphere_radius");
  }
  // mpas_bootstrapping.F:313-339: file_id, parent_id, mesh_spec (a missing one is a warning)
  std::string file_id, parent_id;
  if (f.has_att("", "file_id")) {
    file_id = text_att(f, "file_id");
  } else {
    in.warnings.push_back("Attribute file_id not found; using '' for the previous file_id (mpas_bootstrapping.F:313-318)");
  }
  if (f.has_att("", "parent_id")) {
    parent_id = text_att(f, "parent_id");
  } else {
    in.warnings.push_back("Attribute parent_id not found; setting parent_id to '' (mpas_bootstrapping.F:320-325)");
  }
  m.file_id = file_id;
  // domain % parent_id is character(len=StrKIND): the assignment cuts the joined value at 512 characters
  m.parent_id = rtrim(file_id).empty() ? parent_id : (rtrim(file_id) + "\n" + rtrim(parent_id)).substr(0, kStrKIND);
  if (f.has_att("", "mesh_spec")) {
    m.mesh_spec = text_att(f, "mesh_spec");
  } else {
    in.warnings.push_back("Attribute mesh_spec not found; setting mesh_spec to '0.0' (mpas_bootstrapping.F:333-337)");
    m.mesh_spec = "0.0";
  }
  // mpas_bootstrapping.F:346-379
  if (!m.on_a_sphere) {
    if (!f.has_att("", "is_periodic")) {
      in.warnings.push_back("is_periodic attribute not found; setting is_periodic to .false. (mpas_bootstrapping.F:349-351)");
      m.is_periodic = false;
    } else {
      m.is_periodic = text_att(f, "is_periodic").find("YES") != std::string::npos;
    }
  } else {
    m.is_periodic = false;
  }
  if (m.is_periodic)
    throw Error(f.path() + ": is_periodic = YES: the periodic-plane paths (mpas_fix_periodicity in "
                "mpas_initialize_vectors and mpas_init_reconstruct, x_period / y_period) are not ported "
                "(only the is_periodic = .false. branches are)");
  m.x_period = 0.0;
  m.y_period = 0.0;

  for (const auto& a : f.atts("")) in.source.global_atts.push_back(a);
}

void check_identity_ids(const NcFile& f, const Mesh& m) {
  // The I/O decomposition reads local element i from file record indexTo*ID(i) (add_field_indices.inc), while the
  // bootstrapping builds those IDs from the connectivity values; both agree only when the file's indexTo*ID is
  // 1..n.
  for (const char* name : {"indexToCellID", "indexToEdgeID", "indexToVertexID"}) {
    const std::vector<int> v = raw_int(f, m, name);
    for (std::size_t i = 0; i < v.size(); ++i)
      if (v[i] != static_cast<int>(i) + 1)
        throw Error(f.path() + ": " + name + "(" + std::to_string(i + 1) + ") = " + std::to_string(v[i]) +
                    ": files whose global IDs are not 1..n are not handled");
  }
}

// Every used slot of the connectivity must hold a global ID 1..n. postread_reindex maps anything else to the garbage
// element n+1 (mpas_stream_manager.F:5627-5636), which only the boundary branches (cell == nCells+1) of MPAS handle,
// and those are not ported. MPAS marks a missing neighbour across a mesh boundary with 0. Unused slots
// (j > nEdgesOnCell(i), j > nEdgesOnEdge(i)) may hold anything: they are reindexed and written back as 0.
void check_slots(const NcFile& f, const char* name, const std::vector<int>& a, int inner, int outer,
                 const std::vector<int>* used, int n, const char* boundary_element) {
  for (int i = 1; i <= outer; ++i) {
    const int nu = used ? (*used)[i - 1] : inner;
    for (int j = 1; j <= nu; ++j) {
      const int v = a[static_cast<std::size_t>(j - 1) + static_cast<std::size_t>(inner) * (i - 1)];
      if (v == 0 && boundary_element)
        throw Error(f.path() + ": " + boundary_element + " " + std::to_string(i) + " lies on the mesh boundary (" +
                    name + " = 0): meshes with boundaries (limited-area or non-periodic planar) are not handled; "
                    "their boundary branches (cell == nCells+1) are not ported");
      if (v < 1 || v > n)
        throw Error(f.path() + ": " + name + "(" + std::to_string(j) + ", " + std::to_string(i) + ") = " +
                    std::to_string(v) + " outside 1.." + std::to_string(n) +
                    ": MPAS would turn it into the garbage element (mesh not handled)");
    }
  }
}

void check_count(const NcFile& f, const char* name, const std::vector<int>& c, int lo, int hi, const char* bound) {
  for (std::size_t i = 0; i < c.size(); ++i)
    if (c[i] < lo || c[i] > hi)
      throw Error(f.path() + ": " + name + "(" + std::to_string(i + 1) + ") = " + std::to_string(c[i]) + " outside " +
                  std::to_string(lo) + ".." + bound);
}

// Connectivity fields and the element kind their values index (mpas_stream_manager.F:5455-5610)
struct ConnInfo {
  const char* name;
  Element space;
};
const ConnInfo kConn[] = {{"cellsOnCell", Element::Cell},    {"edgesOnCell", Element::Edge},
                          {"verticesOnCell", Element::Vertex}, {"cellsOnEdge", Element::Cell},
                          {"verticesOnEdge", Element::Vertex}, {"edgesOnEdge", Element::Edge},
                          {"cellsOnVertex", Element::Cell},   {"edgesOnVertex", Element::Edge}};

const FArray<int>& index_array(const Mesh& m, Element e) {
  return e == Element::Cell ? m.indexToCellID : e == Element::Edge ? m.indexToEdgeID : m.indexToVertexID;
}

void read_field(const NcFile& f, MpasInput& in, const std::string& name) {
  FieldRef ref = field_ref(in, name);
  const RegistryVar& rv = registry_var(name);
  if (name == "indexToCellID" || name == "indexToEdgeID" || name == "indexToVertexID") {
    // read through the decomposition: equal to the bootstrapped list (identity checked), garbage element 0
    check_var(f, in.mesh, rv, name);
  } else if (ref.r) {
    read_into(f, in.mesh, rv, *ref.r);
  } else if (ref.i) {
    read_into(f, in.mesh, rv, *ref.i);
    for (const auto& c : kConn) {
      if (name == c.name) {
        // postread_reindex (mpas_stream_manager.F:5413-5645), all innerDim slots of elements 1..outerDim
        const Mesh& m = in.mesh;
        const int outer = static_cast<int>(ref.i->extent(1)) - 1;
        const int space = c.space == Element::Cell ? m.nCells : c.space == Element::Edge ? m.nEdges : m.nVertices;
        postread_reindex(*ref.i, ref.i->extent(0), outer, index_array(m, c.space), space);
      }
    }
  } else {
    read_text(f, in.mesh, rv, *ref.c);
  }
  in.source.vars_read.push_back(name);
  in.source.var_atts[name] = f.atts(name);
}

void read_mesh(const NcFile& f, MpasInput& in, bool init) {
  read_dims_and_attributes(f, in);
  Mesh& m = in.mesh;
  for (const auto& n : kMeshRequired)
    if (!f.has_var(n)) throw Error(f.path() + ": required mesh variable " + n + " is missing");
  check_identity_ids(f, m);
  const std::vector<int> nEdgesOnCell = raw_int(f, m, "nEdgesOnCell");
  const std::vector<int> cellsOnCell = raw_int(f, m, "cellsOnCell");
  const std::vector<int> edgesOnCell = raw_int(f, m, "edgesOnCell");
  const std::vector<int> verticesOnCell = raw_int(f, m, "verticesOnCell");
  const std::vector<int> cellsOnEdge = raw_int(f, m, "cellsOnEdge");
  const std::vector<int> cellsOnVertex = raw_int(f, m, "cellsOnVertex");
  {
    const std::vector<int> nEdgesOnEdge = raw_int(f, m, "nEdgesOnEdge");
    check_count(f, "nEdgesOnCell", nEdgesOnCell, 1, m.maxEdges, "maxEdges");
    check_count(f, "nEdgesOnEdge", nEdgesOnEdge, 0, m.maxEdges2, "maxEdges2");
    check_slots(f, "cellsOnCell", cellsOnCell, m.maxEdges, m.nCells, &nEdgesOnCell, m.nCells, "cell");
    check_slots(f, "cellsOnEdge", cellsOnEdge, 2, m.nEdges, nullptr, m.nCells, "edge");
    check_slots(f, "cellsOnVertex", cellsOnVertex, m.vertexDegree, m.nVertices, nullptr, m.nCells, "vertex");
    check_slots(f, "edgesOnCell", edgesOnCell, m.maxEdges, m.nCells, &nEdgesOnCell, m.nEdges, nullptr);
    check_slots(f, "verticesOnCell", verticesOnCell, m.maxEdges, m.nCells, &nEdgesOnCell, m.nVertices, nullptr);
    check_slots(f, "verticesOnEdge", raw_int(f, m, "verticesOnEdge"), 2, m.nEdges, nullptr, m.nVertices, nullptr);
    check_slots(f, "edgesOnEdge", raw_int(f, m, "edgesOnEdge"), m.maxEdges2, m.nEdges, &nEdgesOnEdge, m.nEdges,
                nullptr);
    check_slots(f, "edgesOnVertex", raw_int(f, m, "edgesOnVertex"), m.vertexDegree, m.nVertices, nullptr, m.nEdges,
                nullptr);
  }
  build_one_rank_order(m, nEdgesOnCell, edgesOnCell, verticesOnCell, cellsOnEdge, cellsOnVertex);

  if (init) {
    if (!f.has_dim("nVertLevels")) throw Error(f.path() + ": dimension nVertLevels is missing (not an init file?)");
    m.nVertLevels = mesh_dim(f, "nVertLevels");
    if (f.has_dim("nVertLevelsP1") && f.dim_len("nVertLevelsP1") != static_cast<std::size_t>(m.nVertLevels) + 1)
      throw Error(f.path() + ": nVertLevelsP1 != nVertLevels + 1");
  }
  for (const auto& n : kMeshRequired) read_field(f, in, n);
  for (const auto& n : kMeshInitOnly) {
    const RegistryVar& rv = registry_var(n);
    bool vertical = false;
    for (const auto& d : rv.dims) vertical |= (d == "nVertLevels" || d == "nVertLevelsP1");
    if (f.has_var(n)) {
      if (!init && vertical) continue;  // vertical fields, read as a grid without a vertical mesh: not read
      read_field(f, in, n);
    } else if (init) {
      throw Error(f.path() + ": required variable " + n + " of the input stream is missing");
    } else {
      // allocate with the Registry default, as MPAS allocates every mesh-pool field
      if (vertical) continue;
      FieldRef ref = field_ref(in, n);
      std::vector<int> ext;
      for (const auto& d : rv.space_dims()) ext.push_back(static_cast<int>(dim_size(m, d, true)));
      if (ext.empty()) ext.push_back(1);
      if (ref.r) ref.r->allocate(ext, default_of<real>(rv));
      if (ref.i) ref.i->allocate(ext, default_of<int>(rv));
    }
  }
  // limited-area meshes (bdyMask* /= 0) are not handled: config_apply_lbcs = false and the boundary code
  // (mpas_atm_boundaries.F) is not ported
  for (auto* a : {&m.bdyMaskCell, &m.bdyMaskEdge, &m.bdyMaskVertex}) {
    if (!a->allocated()) continue;
    for (std::size_t k = 0; k < a->size(); ++k)
      if (a->data()[k] != 0)
        throw Error(f.path() + ": bdyMask* has nonzero entries: limited-area meshes are not handled");
  }
}

void record_dims(const NcFile& f, SourceInfo& s) {
  for (const auto& v : s.vars_read)
    for (const auto& d : f.var(v).dims) s.dims[d] = {f.dim_len(d), f.dim_unlimited(d)};
}

}  // namespace

std::vector<std::string> handled_mesh_fields() {
  std::vector<std::string> v = kMeshRequired;
  v.insert(v.end(), kMeshInitOnly.begin(), kMeshInitOnly.end());
  return v;
}
std::vector<std::string> handled_init_fields() {
  std::vector<std::string> v = handled_mesh_fields();
  v.insert(v.end(), kState.begin(), kState.end());
  return v;
}

FieldRef field_ref(MpasInput& in, const std::string& name) {
  Mesh& m = in.mesh;
  State& s = in.state;
  FieldRef r;
  r.name = name;
#define CPPMPAS_R(x) \
  if (name == #x) { r.r = &m.x; return r; }
#define CPPMPAS_I(x) \
  if (name == #x) { r.i = &m.x; return r; }
#define CPPMPAS_S(x) \
  if (name == #x) { r.r = &s.x; return r; }
#define CPPMPAS_D(x) \
  if (name == #x) { r.r = &in.diag.x; return r; }
  CPPMPAS_I(indexToCellID) CPPMPAS_I(indexToEdgeID) CPPMPAS_I(indexToVertexID)
  CPPMPAS_I(nEdgesOnCell) CPPMPAS_I(nEdgesOnEdge) CPPMPAS_I(cellsOnCell) CPPMPAS_I(edgesOnCell)
  CPPMPAS_I(verticesOnCell) CPPMPAS_I(cellsOnEdge) CPPMPAS_I(verticesOnEdge) CPPMPAS_I(edgesOnEdge)
  CPPMPAS_I(cellsOnVertex) CPPMPAS_I(edgesOnVertex) CPPMPAS_I(bdyMaskCell) CPPMPAS_I(bdyMaskEdge)
  CPPMPAS_I(bdyMaskVertex)
  CPPMPAS_R(latCell) CPPMPAS_R(lonCell) CPPMPAS_R(xCell) CPPMPAS_R(yCell) CPPMPAS_R(zCell)
  CPPMPAS_R(latEdge) CPPMPAS_R(lonEdge) CPPMPAS_R(xEdge) CPPMPAS_R(yEdge) CPPMPAS_R(zEdge)
  CPPMPAS_R(latVertex) CPPMPAS_R(lonVertex) CPPMPAS_R(xVertex) CPPMPAS_R(yVertex) CPPMPAS_R(zVertex)
  CPPMPAS_R(weightsOnEdge) CPPMPAS_R(dvEdge) CPPMPAS_R(dcEdge) CPPMPAS_R(angleEdge) CPPMPAS_R(areaCell)
  CPPMPAS_R(areaTriangle) CPPMPAS_R(kiteAreasOnVertex) CPPMPAS_R(fEdge) CPPMPAS_R(fVertex) CPPMPAS_R(meshDensity)
  CPPMPAS_R(nominalMinDc) CPPMPAS_R(edgeNormalVectors) CPPMPAS_R(localVerticalUnitVectors)
  CPPMPAS_R(cellTangentPlane) CPPMPAS_R(cf1) CPPMPAS_R(cf2) CPPMPAS_R(cf3) CPPMPAS_R(zgrid) CPPMPAS_R(rdzw)
  CPPMPAS_R(dzu) CPPMPAS_R(rdzu) CPPMPAS_R(fzm) CPPMPAS_R(fzp) CPPMPAS_R(zxu) CPPMPAS_R(zz) CPPMPAS_R(zb)
  CPPMPAS_R(zb3) CPPMPAS_R(dss) CPPMPAS_R(deriv_two) CPPMPAS_R(coeffs_reconstruct)
  CPPMPAS_S(u_init) CPPMPAS_S(v_init) CPPMPAS_S(t_init) CPPMPAS_S(u) CPPMPAS_S(w) CPPMPAS_S(rho)
  CPPMPAS_S(theta) CPPMPAS_S(rho_base) CPPMPAS_S(theta_base) CPPMPAS_S(surface_pressure)
  // computed by atm_mpas_init_block (core_atmosphere/atm_init.hpp); allocated only after it ran
  CPPMPAS_R(invAreaCell) CPPMPAS_R(invDvEdge) CPPMPAS_R(invDcEdge) CPPMPAS_R(invAreaTriangle)
  CPPMPAS_R(edgesOnCell_sign) CPPMPAS_R(edgesOnVertex_sign) CPPMPAS_R(zb_cell) CPPMPAS_R(zb3_cell)
  CPPMPAS_I(kiteForCell) CPPMPAS_I(advCellsForEdge) CPPMPAS_I(nAdvCellsForEdge) CPPMPAS_R(adv_coefs)
  CPPMPAS_R(adv_coefs_3rd) CPPMPAS_R(deformation_coef_c2) CPPMPAS_R(deformation_coef_s2)
  CPPMPAS_R(deformation_coef_cs) CPPMPAS_R(deformation_coef_c) CPPMPAS_R(deformation_coef_s)
  CPPMPAS_R(meshScalingDel2) CPPMPAS_R(meshScalingDel4) CPPMPAS_R(meshScalingRegionalCell)
  CPPMPAS_R(meshScalingRegionalEdge) CPPMPAS_R(etp) CPPMPAS_R(etm) CPPMPAS_R(ewp) CPPMPAS_R(ewm)
  CPPMPAS_R(specZoneMaskCell) CPPMPAS_R(specZoneMaskEdge) CPPMPAS_R(specZoneMaskVertex)
  CPPMPAS_I(nearestRelaxationCell)
  CPPMPAS_S(theta_m) CPPMPAS_S(rho_zz) CPPMPAS_S(Time)
  CPPMPAS_D(ru) CPPMPAS_D(rw) CPPMPAS_D(rho_p) CPPMPAS_D(rtheta_base) CPPMPAS_D(rtheta_p) CPPMPAS_D(exner)
  CPPMPAS_D(exner_base) CPPMPAS_D(pressure_p) CPPMPAS_D(pressure_base)
  CPPMPAS_D(v) CPPMPAS_D(rho_edge) CPPMPAS_D(vorticity) CPPMPAS_D(divergence) CPPMPAS_D(ke) CPPMPAS_D(pv_edge)
  CPPMPAS_D(pv_vertex) CPPMPAS_D(pv_cell) CPPMPAS_D(gradPVt) CPPMPAS_D(gradPVn) CPPMPAS_D(uReconstructX)
  CPPMPAS_D(uReconstructY) CPPMPAS_D(uReconstructZ) CPPMPAS_D(uReconstructZonal) CPPMPAS_D(uReconstructMeridional)
#undef CPPMPAS_D
#undef CPPMPAS_R
#undef CPPMPAS_I
#undef CPPMPAS_S
  if (name == "xtime") { r.c = &s.xtime; return r; }
  if (name == "initial_time") { r.c = &s.initial_time; return r; }
  for (std::size_t k = 0; k < s.scalar_names.size(); ++k) {
    if (name == s.scalar_names[k]) {
      r.r = &s.scalars;
      r.scalar_slot = static_cast<int>(k) + 1;
      return r;
    }
  }
  if (!find_registry_var(name)) throw Error("variable '" + name + "' is not defined in MPAS-Atmosphere's Registry");
  throw Error("variable '" + name + "' has no CPPMPAS port (not read or written)");
}

MpasInput read_mpas_grid(const std::string& path) {
  const NcFile f = NcFile::open_read(path);
  MpasInput in;
  read_mesh(f, in, false);
  record_dims(f, in.source);
  return in;
}

MpasInput read_mpas_init(const std::string& path, const std::vector<std::string>& scalars) {
  const NcFile f = NcFile::open_read(path);
  MpasInput in;
  read_mesh(f, in, true);
  Mesh& m = in.mesh;
  for (const auto& n : kState) {
    if (!f.has_var(n)) throw Error(f.path() + ": required variable " + n + " of the input stream is missing");
    read_field(f, in, n);
  }
  // var_array scalars: (num_scalars, nVertLevels, nCells+1), constituents in Registry order
  // (Registry.xml:1700-1749); each active constituent is read from its own NetCDF variable
  std::vector<std::string> order;
  for (const auto& rv : registry_vars())
    if (rv.array == "scalars") order.push_back(rv.name);
  for (auto it = scalars.begin(); it != scalars.end(); ++it) {
    if (std::find(order.begin(), order.end(), *it) == order.end())
      throw Error("scalar '" + *it + "' is not a constituent of MPAS's var_array scalars");
    if (std::find(scalars.begin(), it, *it) != it) throw Error("scalar '" + *it + "' is given twice");
  }
  for (const auto& s : order)
    if (std::find(scalars.begin(), scalars.end(), s) != scalars.end()) in.state.scalar_names.push_back(s);
  // qv belongs to no package (Registry.xml:1701-1703): every MPAS-Atmosphere run has it
  if (std::find(scalars.begin(), scalars.end(), "qv") == scalars.end())
    throw Error("qv is not among the active scalars; MPAS-Atmosphere always has qv (Registry.xml:1701-1703)");
  const int ns = static_cast<int>(in.state.scalar_names.size());
  in.state.scalars.allocate({ns, m.nVertLevels, m.nCells + 1}, 0.0);
  for (int k = 1; k <= ns; ++k) {
    const std::string& n = in.state.scalar_names[k - 1];
    if (!f.has_var(n)) throw Error(f.path() + ": active scalar " + n + " is missing from the input file");
    const RegistryVar& rv = registry_var(n);
    const long record = check_var(f, m, rv, n);
    const std::vector<double> buf = f.get_double(n, record);  // (nCells, nVertLevels) in file order
    for (int i = 1; i <= m.nCells; ++i) {
      const int g = m.cell_order.local_to_file[i];
      for (int lev = 1; lev <= m.nVertLevels; ++lev)
        in.state.scalars(k, lev, i) = buf[static_cast<std::size_t>(g - 1) * m.nVertLevels + (lev - 1)];
    }
    in.source.vars_read.push_back(n);
    in.source.var_atts[n] = f.atts(n);
  }
  // every other variable of the input stream present in the file must be one the port knows it does not need
  const std::vector<std::string> handled = handled_init_fields();
  for (const auto& n : input_stream_vars()) {
    if (n == "scalars") continue;
    if (std::find(handled.begin(), handled.end(), n) != handled.end()) continue;
    if (!f.has_var(n)) continue;
    if (kInputNotUsed.count(n)) continue;
    throw Error(f.path() + ": MPAS reads variable " + n + " from the input stream; CPPMPAS has no port for it");
  }
  in.has_state = true;
  record_dims(f, in.source);
  return in;
}

}  // namespace cppmpas
