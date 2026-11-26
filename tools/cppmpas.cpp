// SPDX-License-Identifier: GPL-3.0-or-later
// cppmpas: command-line access to CPPMPAS's mesh / I-O layer.
//   cppmpas info      <file> [--init] [--scalars qv,qc,qr]
//   cppmpas roundtrip <in> <out> [--init] [--scalars ...] [--fields a,b,...]
//   cppmpas history   <init> <out> --fields a,b,... [--scalars ...] [--file-id ID]
//   cppmpas dump      <in> <out.nc> [--init] [--scalars ...]
//   cppmpas init      <init> <out> --fields a,b,... [--namelist namelist.atmosphere] [--scalars ...] [--file-id ID]
//                     [--unverified]
// init runs what atm_core_init does to the fields of a cold start (core_atmosphere/atm_init.hpp: len_disp, the
// compatibility checks, atm_mpas_init_block, xtime) on an init file and writes the
// given fields as MPAS writes them to a history / restart file (global order, Registry layout).
// dump writes every array as MPAS holds it in memory on one rank: local order, local 1-based indices, the garbage
// element included (dimensions <dim>_mem); compare with MPAS memory dumps after nothing more than a type change.
// The var_array scalars is written whole, as "scalars" (num_scalars, nVertLevels, nCells+1) with the attribute
// constituents; character fields as their StrLen characters.
#include <algorithm>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "core/namelist.hpp"
#include "core_atmosphere/atm_init.hpp"
#include "io/mpas_reader.hpp"
#include "io/mpas_writer.hpp"
#include "io/registry.hpp"

using namespace cppmpas;

namespace {

std::vector<std::string> split(const std::string& s) {
  std::vector<std::string> v;
  std::stringstream ss(s);
  std::string x;
  while (std::getline(ss, x, ',')) {
    if (x.empty()) throw Error("empty name in the list '" + s + "'");
    v.push_back(x);
  }
  if (v.empty() || s.back() == ',') throw Error("empty name in the list '" + s + "'");
  return v;
}

struct Args {
  std::vector<std::string> pos;
  bool init = false, unverified = false;
  std::vector<std::string> scalars{"qv"}, fields;
  std::string file_id, namelist;
};

Args parse(int argc, char** argv) {
  Args a;
  for (int i = 2; i < argc; ++i) {
    const std::string s = argv[i];
    auto next = [&]() -> std::string {
      if (i + 1 >= argc) throw Error("option " + s + " needs a value");
      return argv[++i];
    };
    if (s == "--init") a.init = true;
    else if (s == "--unverified") a.unverified = true;
    else if (s == "--scalars") a.scalars = split(next());
    else if (s == "--fields") a.fields = split(next());
    else if (s == "--file-id") a.file_id = next();
    else if (s == "--namelist") a.namelist = next();
    else if (s.rfind("--", 0) == 0) throw Error("unknown option " + s);
    else a.pos.push_back(s);
  }
  return a;
}

MpasInput load(const Args& a, const std::string& path) {
  MpasInput in = a.init ? read_mpas_init(path, a.scalars) : read_mpas_grid(path);
  for (const auto& w : in.warnings) std::cerr << "MPAS warning: " << w << "\n";
  return in;
}

void dump(MpasInput& in, const std::string& out) {
  NcFile f = NcFile::create(out, NC_64BIT_DATA);
  const Mesh& m = in.mesh;
  std::vector<std::string> names = in.source.vars_read;
  std::vector<std::string> defined;
  auto mem_dim = [&](const std::string& d) -> std::pair<std::string, std::size_t> {
    if (d == "nCells") return {"nCells_mem", static_cast<std::size_t>(m.nCells) + 1};
    if (d == "nEdges") return {"nEdges_mem", static_cast<std::size_t>(m.nEdges) + 1};
    if (d == "nVertices") return {"nVertices_mem", static_cast<std::size_t>(m.nVertices) + 1};
    return {d + "_mem", 0};
  };
  struct Item { std::string name; FieldRef ref; };
  std::vector<Item> items;
  for (const auto& n : names) {
    FieldRef r = field_ref(in, n);
    if (r.scalar_slot > 0) continue;  // written once below, as the whole var_array
    if (r.c) {
      if (std::find(defined.begin(), defined.end(), "StrLen") == defined.end()) {
        f.def_dim("StrLen", 64);
        defined.push_back("StrLen");
      }
      f.def_var(n, NC_CHAR, {"StrLen"});
      items.push_back({n, r});
      continue;
    }
    const FArray<real>* ra = r.r;
    const FArray<int>* ia = r.i;
    const std::vector<int>& ext = ra ? ra->extents() : ia->extents();
    const RegistryVar& rv = registry_var(n);
    std::vector<std::string> sp = rv.space_dims();
    if (sp.empty()) sp.push_back("ONE");
    std::vector<std::string> dims;
    for (std::size_t k = sp.size(); k-- > 0;) {
      auto [dn, len] = mem_dim(sp[k]);
      if (len == 0) len = static_cast<std::size_t>(ext[k]);
      if (std::find(defined.begin(), defined.end(), dn) == defined.end()) {
        f.def_dim(dn, len);
        defined.push_back(dn);
      }
      dims.push_back(dn);
    }
    f.def_var(n, ra ? NC_DOUBLE : NC_INT, dims);
    items.push_back({n, r});
  }
  const FArray<real>& sc = in.state.scalars;
  if (sc.allocated()) {
    f.def_dim("num_scalars_mem", static_cast<std::size_t>(sc.extent(0)));
    if (std::find(defined.begin(), defined.end(), "nVertLevels_mem") == defined.end())
      f.def_dim("nVertLevels_mem", static_cast<std::size_t>(sc.extent(1)));
    if (std::find(defined.begin(), defined.end(), "nCells_mem") == defined.end())
      f.def_dim("nCells_mem", static_cast<std::size_t>(sc.extent(2)));
    f.def_var("scalars", NC_DOUBLE, {"nCells_mem", "nVertLevels_mem", "num_scalars_mem"});
    std::string names_att;
    for (const auto& s : in.state.scalar_names) names_att += (names_att.empty() ? "" : " ") + s;
    f.put_att_text("scalars", "constituents", names_att);
  }
  f.enddef();
  if (sc.allocated()) f.put_double("scalars", std::vector<double>(sc.data(), sc.data() + sc.size()));
  for (const auto& it : items) {
    if (it.ref.c) {
      std::vector<char> v(it.ref.c->begin(), it.ref.c->end());
      v.resize(64, ' ');
      f.put_char(it.name, v);
    } else if (it.ref.r)
      f.put_double(it.name, std::vector<double>(it.ref.r->data(), it.ref.r->data() + it.ref.r->size()));
    else
      f.put_int(it.name, std::vector<int>(it.ref.i->data(), it.ref.i->data() + it.ref.i->size()));
  }
  f.close();
}

int run(int argc, char** argv) {
  if (argc < 2) throw Error("usage: cppmpas info|roundtrip|history|init|dump ...");
  const std::string cmd = argv[1];
  Args a = parse(argc, argv);
  if (cmd == "info") {
    if (a.pos.size() != 1) throw Error("usage: cppmpas info <file> [--init] [--scalars ...]");
    MpasInput in = load(a, a.pos[0]);
    const Mesh& m = in.mesh;
    std::cout << "file " << in.source.path << " (" << NcFile::format_name(in.source.format) << ")\n"
              << "nCells " << m.nCells << " nEdges " << m.nEdges << " nVertices " << m.nVertices << " maxEdges "
              << m.maxEdges << " maxEdges2 " << m.maxEdges2 << " vertexDegree " << m.vertexDegree << " nVertLevels "
              << m.nVertLevels << "\n"
              << "on_a_sphere " << m.on_a_sphere << " sphere_radius " << m.sphere_radius << " is_periodic "
              << m.is_periodic << "\n";
    for (auto [label, o] : {std::pair{"cells", &m.cell_order}, {"edges", &m.edge_order}, {"vertices", &m.vertex_order}}) {
      std::cout << label << " local order " << (o->identity() ? "= file order" : "!= file order") << ", first:";
      for (int i = 1; i <= std::min(12, o->n()); ++i) std::cout << " " << o->local_to_file[i];
      std::cout << "\n";
    }
    std::cout << "fields read " << in.source.vars_read.size() << "\n";
  } else if (cmd == "roundtrip") {
    if (a.pos.size() != 2) throw Error("usage: cppmpas roundtrip <in> <out> [--init] [--scalars ...] [--fields ...]");
    MpasInput in = load(a, a.pos[0]);
    write_round_trip(in, a.pos[1], a.fields);
  } else if (cmd == "history") {
    if (a.pos.size() != 2 || a.fields.empty())
      throw Error("usage: cppmpas history <init> <out> --fields a,b,... [--scalars ...] [--file-id ID]");
    a.init = true;
    MpasInput in = load(a, a.pos[0]);
    HistoryOptions opt;
    opt.file_id = a.file_id;
    write_history(in, a.pos[1], a.fields, opt);
    std::cerr << "NOTE: the namelist options MPAS writes as global attributes (mpas_subdriver.F:532-557) are not "
                 "written yet\n";
  } else if (cmd == "init") {
    if (a.pos.size() != 2 || a.fields.empty())
      throw Error("usage: cppmpas init <init> <out> --fields a,b,... [--namelist file] [--scalars ...] [--file-id ID] "
                  "[--unverified]");
    a.init = true;
    InitConfig cfg = read_init_config(a.namelist);
    // the scalars MPAS allocates are those the namelist's packages activate (Registry.xml:1701-1748)
    const std::vector<std::string> want = active_scalars(cfg);
    if (a.scalars != want) {
      std::string w;
      for (const auto& x : want) w += (w.empty() ? "" : ",") + x;
      throw Error("--scalars must be " + w + ": the scalars the namelist activates (MPAS allocates exactly those)");
    }
    MpasInput in = load(a, a.pos[0]);
    // an array this port computes that is not yet proven against MPAS needs --unverified (CONSTITUTION.md 3.1)
    const auto& computed = init_block_computed_fields();
    const auto verified = init_block_verified_fields(in.mesh.nVertLevels, mesh_scaling_unverified(in.mesh, cfg));
    for (const auto& f : a.fields) {
      const bool is_computed = std::find(computed.begin(), computed.end(), f) != computed.end();
      if (is_computed && !a.unverified && std::find(verified.begin(), verified.end(), f) == verified.end())
        throw Error("field " + f + " is computed by a port not yet verified against MPAS output; --unverified "
                    "writes it anyway");
    }
    std::vector<std::string> log;
    atm_core_init(in.mesh, in.state, in.diag, cfg, log);
    for (const auto& l : log) std::cerr << "MPAS log: " << l << "\n";
    HistoryOptions opt;
    opt.file_id = a.file_id;
    write_history(in, a.pos[1], a.fields, opt);
  } else if (cmd == "dump") {
    if (a.pos.size() != 2) throw Error("usage: cppmpas dump <in> <out.nc> [--init] [--scalars ...]");
    MpasInput in = load(a, a.pos[0]);
    dump(in, a.pos[1]);
  } else {
    throw Error("unknown command " + cmd);
  }
  return 0;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    return run(argc, argv);
  } catch (const std::exception& e) {
    std::cerr << "ERROR: " << e.what() << "\n";
    return 1;
  }
}
