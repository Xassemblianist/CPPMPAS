// SPDX-License-Identifier: GPL-3.0-or-later
// cppmpas: command-line access to CPPMPAS's mesh / I-O layer.
//   cppmpas info      <file> [--init] [--scalars qv,qc,qr]
//   cppmpas roundtrip <in> <out> [--init] [--scalars ...] [--fields a,b,...]
//   cppmpas history   <init> <out> --fields a,b,... [--scalars ...] [--file-id ID]
//   cppmpas dump      <in> <out.nc> [--init] [--scalars ...]
//   cppmpas init      <init> <out> --fields a,b,... [--namelist namelist.atmosphere] [--scalars ...] [--file-id ID]
//                     [--unverified]
//   cppmpas run       <init> <outdir> --namelist namelist.atmosphere --steps N [--scalars ...] [--file-id ID]
//                     [--every K] [--write-at N1,N2,... --write-fields a,b,...]
// init runs what atm_core_init does to the fields of a cold start (core_atmosphere/atm_init.hpp: len_disp, the
// compatibility checks, atm_mpas_init_block, xtime) on an init file and writes the
// given fields as MPAS writes them to a history / restart file (global order, Registry layout).
// With config_do_restart = true in the namelist, run's <init> is a restart file, read as MPAS reads its restart stream
// at config_start_time (io/mpas_restart.hpp), and the run continues from it.
// run integrates N time steps of the dycore (core_atmosphere/atm_dycore.hpp: atm_core_init, then per step
// atm_timestep, the time-level shift and atm_compute_output_diagnostics, as atm_core_run does) and writes, at the
// start and after every step, the fields of the references' "bisection" stream (xtime, Time, u, w, rho_zz, theta_m,
// rho, theta, pressure_p, scalars) to <outdir>/step.$Y-$M-$D_$h.$m.$s.nc, as MPAS names them; --every K writes
// them at the start and after every K-th step only (0: never). --write-at writes, after each listed step, the
// --write-fields to <outdir>/state.$Y-$M-$D_$h.$m.$s.nc (to compare with MPAS's history / restart streams).
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
#include "core_atmosphere/atm_dycore.hpp"
#include "core_atmosphere/atm_init.hpp"
#include "io/mpas_reader.hpp"
#include "io/mpas_restart.hpp"
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

long whole(const std::string& opt, const std::string& v) {
  std::size_t used = 0;
  long n = -1;
  try {
    n = std::stol(v, &used);
  } catch (const std::exception&) {
    used = 0;
  }
  if (used != v.size() || n < 0 || n > 1000000) throw Error(opt + " needs a whole number 0..1000000, not '" + v + "'");
  return n;
}

struct Args {
  std::vector<std::string> pos;
  bool init = false, unverified = false;
  std::vector<std::string> scalars{"qv"}, fields;
  std::string file_id, namelist;
  long steps = -1, every = 1;
  std::vector<long> write_at;
  std::vector<std::string> write_fields;
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
    else if (s == "--steps") {
      a.steps = whole(s, next());
    } else if (s == "--every") {
      a.every = whole(s, next());
    } else if (s == "--write-at") {
      for (const auto& x : split(next())) a.write_at.push_back(whole(s, x));
    } else if (s == "--write-fields") {
      a.write_fields = split(next());
    }
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
  if (argc < 2) throw Error("usage: cppmpas info|roundtrip|history|init|run|dump ...");
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
    if (cfg.config_do_restart)
      throw Error("config_do_restart = true: init reads an input-stream file (a cold start); a restart run is "
                  "'cppmpas run' on the restart file");
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
  } else if (cmd == "run") {
    if (a.pos.size() != 2 || a.steps < 0 || a.namelist.empty() || !a.fields.empty())
      throw Error("usage: cppmpas run <init> <outdir> --namelist file --steps N [--scalars ...] [--file-id ID] "
                  "[--every K] [--write-at N1,... --write-fields a,...]");
    a.init = true;
    InitConfig cfg = read_init_config(a.namelist);
    const std::vector<std::string> want = active_scalars(cfg);
    if (a.scalars != want) throw Error("--scalars must be the scalars the namelist activates");
    // config_do_restart: <init> is the restart file, read as the restart stream (io/mpas_restart.hpp)
    RestartInput rst;
    if (cfg.config_do_restart) {
      rst = read_mpas_restart(a.pos[0], a.scalars, atm_start_time_stamp(cfg));
      for (const auto& w : rst.in.warnings) std::cerr << "MPAS warning: " << w << "\n";
    } else {
      rst.in = load(a, a.pos[0]);
    }
    MpasInput& in = rst.in;
    const std::string xtime_read = in.state.xtime;
    std::vector<std::string> log;
    atm_core_init(in.mesh, in.state, in.diag, cfg, log);
    for (const auto& l : log) std::cerr << "MPAS log: " << l << "\n";
    const DycoreConfig dcfg = read_dycore_config(a.namelist, cfg);
    check_dycore_supported(dcfg, in.mesh, in.state);  // before any file is written
    // the input stream, too, is read at the start time (the restart reader checks its own)
    if (!cfg.config_do_restart) check_stream_time(a.pos[0], xtime_read, atm_start_time_stamp(cfg));
    DycoreFields f;
    allocate_dycore_fields(in.mesh, in.state, f);
    if (cfg.config_do_restart) restart_dycore_fields(rst, f);
    microphysics_init(a.namelist, in.mesh, dcfg.config_dt, f.mp);  // mpas_atm_core.F:565-597 (SC: Kessler)
    const std::string start = in.state.xtime.substr(0, 19);
    std::vector<std::string> fields = {"xtime", "Time", "u", "w", "rho_zz", "theta_m", "rho", "theta", "pressure_p"};
    for (const auto& sn : in.state.scalar_names) fields.push_back(sn);
    if (a.write_at.empty() != a.write_fields.empty()) throw Error("--write-at and --write-fields go together");
    for (long n : a.write_at)
      if (n > a.steps) throw Error("--write-at " + std::to_string(n) + " is after the last step");
    for (const auto& w : a.write_fields) (void)field_ref(in, w);  // a field without a port: refused before the run
    auto write_step = [&](long n) {
      atm_compute_output_diagnostics(in.mesh, in.state);  // mpas_atm_core.F:686 (t = 0) and :848
      in.state.Time.allocate({1}, f.Time1);
      std::string ts = in.state.xtime.substr(0, 19);
      for (auto& c : ts)
        if (c == ':') c = '.';
      HistoryOptions opt;
      opt.file_id = a.file_id;
      if (a.every > 0 && n % a.every == 0) write_history(in, a.pos[1] + "/step." + ts + ".nc", fields, opt);
      if (std::find(a.write_at.begin(), a.write_at.end(), n) != a.write_at.end())
        write_history(in, a.pos[1] + "/state." + ts + ".nc", a.write_fields, opt);
    };
    write_step(0);
    const long long dt_s = static_cast<long long>(dcfg.config_dt);
    for (long n = 1; n <= a.steps; ++n) {
      atm_timestep(dcfg, in.mesh, in.state, in.diag, f, dcfg.config_dt, add_seconds(start, n * dt_s));
      shift_time_levels(in.state, f);  // mpas_atm_core.F:823
      write_step(n);
    }
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
