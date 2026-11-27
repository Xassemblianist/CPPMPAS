// SPDX-License-Identifier: GPL-3.0-or-later
// cppmpas: command-line access to CPPMPAS's mesh / I-O layer.
//   cppmpas info      <file> [--init] [--scalars qv,qc,qr]
//   cppmpas roundtrip <in> <out> [--init] [--scalars ...] [--fields a,b,...]
//   cppmpas history   <init> <out> --fields a,b,... [--scalars ...] [--file-id ID]
//   cppmpas dump      <in> <out.nc> [--init] [--scalars ...]
//   cppmpas init      <init> <out> --fields a,b,... [--namelist namelist.atmosphere] [--scalars ...] [--file-id ID]
//                     [--unverified]
//   cppmpas run       <init> <outdir> --namelist namelist.atmosphere --steps N [--scalars ...] [--file-id ID]
//                     [--every K] [--write-at N1,N2,... --write-fields a,b,...] [--tables DIR --rrtmg DIR]
//                     [--streams streams.atmosphere]
// With config_physics_suite = 'mesoscale_reference', run reads the input stream's physics fields and runs MPAS's
// physics (physics/mpas_atmphys_driver.hpp): --tables names the directory with LANDUSE.TBL, VEGPARM.TBL, SOILPARM.TBL,
// GENPARM.TBL and OZONE_*.TBL (the reference build's tables), --rrtmg the one with rrtmg_lw_tables.bin and
// rrtmg_sw_tables.bin (the RRTMG tables MPAS's initialisation computes, physics/mpas_atmphys_radiation.hpp). The
// physics runs on the host between the dycore's kernels.
//   cppmpas device    (the backend: "cpu", or the CUDA device and a check that kernels run on it)
// The CUDA build of this file is cppmpas_gpu: the same commands with every kernel on the GPU (core/backend.hpp); it
// exits with status 77 (CTest's skip code) and a SKIP message when there is no usable CUDA device.
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
// --streams streams.atmosphere (with the physics suite) writes MPAS's output streams instead, as atm_core_run does
// (io/mpas_streams.hpp: every active output stream whose alarm rings, its fields, file names and the namelist options
// as global attributes; the physics diagnostics alarm from the output / diagnostics intervals, MG:615-645; the
// diagnostics manager; atm_reset_diagnostics; restart_timestamp). The run lasts config_run_duration (--steps may be
// left out). Without --streams a --write-at file stands for a ringing output stream of the diagnostics manager.
// dump writes every array as MPAS holds it in memory on one rank: local order, local 1-based indices, the garbage
// element included (dimensions <dim>_mem); compare with MPAS memory dumps after nothing more than a type change.
// The var_array scalars is written whole, as "scalars" (num_scalars, nVertLevels, nCells+1) with the attribute
// constituents; character fields as their StrLen characters.
#include <algorithm>
#include <cstdio>
#include <initializer_list>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "core/backend.hpp"
#include "core/namelist.hpp"
#include "core/views.hpp"
#include "core_atmosphere/atm_dycore.hpp"
#include "core_atmosphere/atm_init.hpp"
#include "io/mpas_reader.hpp"
#include "io/mpas_restart.hpp"
#include "io/mpas_streams.hpp"
#include "io/mpas_writer.hpp"
#include "io/registry.hpp"
#include "physics/mpas_atmphys_driver.hpp"
#include "core_atmosphere/diagnostics/mpas_atm_diagnostics.hpp"

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
  std::string file_id, namelist, tables, rrtmg, streams;
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
    else if (s == "--tables") a.tables = next();
    else if (s == "--rrtmg") a.rrtmg = next();
    else if (s == "--streams") a.streams = next();
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
    // explicit element type: nvcc 12 cannot deduce a braced list of std::pair
    using Order = std::pair<const char*, decltype(&m.cell_order)>;
    for (auto [label, o] : std::initializer_list<Order>{{"cells", &m.cell_order}, {"edges", &m.edge_order},
                                                        {"vertices", &m.vertex_order}}) {
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
    if (a.pos.size() != 2 || (a.steps < 0 && a.streams.empty()) || a.namelist.empty() || !a.fields.empty())
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
    // the physics configuration: physics_namelist_check runs in atm_setup_packages, before the fields are read
    phys::PhysicsRun pr;
    pr.cfg = phys::read_physics_config(a.namelist, cfg, log);
    const bool physics = pr.cfg.config_physics_suite != "none";
    if (physics) {
      if (cfg.config_do_restart) throw Error("run: a restart with the physics suite has no verified port");
      if (a.tables.empty() || a.rrtmg.empty()) throw Error("run: the physics suite needs --tables and --rrtmg");
      phys::check_physics_supported(pr.cfg);
      const std::map<std::string, int> dims = {{"nLags", pr.cfg.input_soil_temperature_lag},
                                               {"nOznLevels", pr.cfg.noznlev},
                                               {"nAerLevels", pr.cfg.naerlev},
                                               {"cam_dim1", pr.cfg.camdim1},
                                               {"nMonths", pr.cfg.months},
                                               {"nSoilLevels", pr.cfg.num_soil_layers},
                                               {"nIsoLevelsT", 5},    // Registry_isobaric.xml:6-10
                                               {"nIsoLevelsZ", 13}};
      read_input_pools(a.pos[0], in, phys::physics_packages(pr.cfg), dims);
      for (const auto& w : in.warnings) std::cerr << "MPAS warning: " << w << "\n";
      pr.table_dir = a.tables;
      pr.rrtmg_lw_tables = a.rrtmg + "/rrtmg_lw_tables.bin";
      pr.rrtmg_sw_tables = a.rrtmg + "/rrtmg_sw_tables.bin";
    }
    atm_core_init(in.mesh, in.state, in.diag, cfg, log);
    pr.cfg.config_len_disp = cfg.config_len_disp;  // as atm_core_init left it (mpas_atm_core.F:155-201)
    const DycoreConfig dcfg = read_dycore_config(a.namelist, cfg);
    check_dycore_supported(dcfg, in.mesh, in.state);  // before any file is written
    // the input stream, too, is read at the start time (the restart reader checks its own)
    if (!cfg.config_do_restart) check_stream_time(a.pos[0], xtime_read, atm_start_time_stamp(cfg));
    DycoreFields f;
    allocate_dycore_fields(in.mesh, in.state, f);
    if (cfg.config_do_restart) restart_dycore_fields(rst, f);
    const std::string start = in.state.xtime.substr(0, 19);
    int itimestep = 1;
    phys::DycorePhysicsArrays pa;
    if (physics) {
      phys::physics_compatibility_check(pr.cfg, in.mesh, in.pools);  // mpas_atm_core.F:1597-1603
      phys::physics_core_init(pr, in.mesh, in.state, in.diag, in.pools, start, log);  // mpas_atm_core.F:565-597
      pa.rthdynten = &f.rthdynten;
      pa.tend_ru_physics = &f.tend_ru_physics;
      pa.tend_rtheta_physics = &f.tend_rtheta_physics;
      pa.tend_rho_physics = &f.tend_rho_physics;
      pa.tend_scalars = &f.tend_scalars;
      pa.rt_diabatic_tend = &f.rt_diabatic_tend;
      Mesh& mesh = in.mesh;
      RegistryPools& pools = in.pools;
      f.phys.get_tend = [&](const State& s, Diag& d, DycoreFields& ff) {
        phys::physics_get_tend_for_dynamics(pr, mesh, s, d, pools, ff.rho_zz2, pa);
      };
      f.phys.dyn_tend_for_convection = [&](const State& s, DycoreFields& ff, real dt) {
        phys::dynamics_tendencies_for_convection(pr, mesh, s, pools, ff.scalars2, ff.theta_m2, dt,
                                                 dcfg.config_monotonic, pa);
      };
      f.phys.microphysics = [&](State& s, Diag& d, DycoreFields& ff) {
        phys::microphysics_for_dynamics(pr, mesh, s, d, pools, ff.theta_m2, ff.rho_zz2, ff.scalars2, ff.w2, itimestep,
                                        pa);
      };
    } else {
      microphysics_init(a.namelist, in.mesh, dcfg.config_dt, f.mp);  // mpas_atm_core.F:565-597 (SC: Kessler)
    }
    for (const auto& l : log) std::cerr << "MPAS log: " << l << "\n";
    log.clear();
    std::vector<std::string> fields = {"xtime", "Time", "u", "w", "rho_zz", "theta_m", "rho", "theta", "pressure_p"};
    for (const auto& sn : in.state.scalar_names) fields.push_back(sn);
    if (a.write_at.empty() != a.write_fields.empty()) throw Error("--write-at and --write-fields go together");
    for (long n : a.write_at)
      if (n > a.steps) throw Error("--write-at " + std::to_string(n) + " is after the last step");
    // the dycore's arrays of the diag / tend / tend_physics pools (Registry names)
    in.held = {{"ru_p", &f.ru_p}, {"rw_p", &f.rw_p}, {"rt_diabatic_tend", &f.rt_diabatic_tend},
               {"rthdynten", &f.rthdynten}, {"rtheta_pp", &f.rtheta_pp}, {"rho_pp", &f.rho_pp}};
    for (const auto& w : a.write_fields) (void)field_ref(in, w);  // a field without a port: refused before the run
    // --streams: MPAS's output streams (io/mpas_streams.hpp) written as atm_core_run writes them
    // (mpas_atm_core.F:625-919), with the namelist options as global attributes (mpas_subdriver.F:532-557)
    const bool streams_mode = !a.streams.empty();
    StreamManager sm;
    HistoryOptions stream_opt;
    stream_opt.file_id = a.file_id;
    if (streams_mode) {
      if (!physics) throw Error("run --streams: only the physics suite's runs have a verified port of the streams");
      sm = make_stream_manager(read_streams_file(a.streams), in.pools.packages, in.state.scalar_names,
                               cfg.config_do_restart);
      for (const auto& na : namelist_attributes(a.namelist)) {
        std::string v = na.value;
        // physics_namelist_check (atm_setup_packages) runs before add_stream_attributes: the resolved schemes
        const std::map<std::string, const std::string*> resolved = {
            {"config_microp_scheme", &pr.cfg.config_microp_scheme},
            {"config_convection_scheme", &pr.cfg.config_convection_scheme},
            {"config_lsm_scheme", &pr.cfg.config_lsm_scheme},
            {"config_pbl_scheme", &pr.cfg.config_pbl_scheme},
            {"config_gwdo_scheme", &pr.cfg.config_gwdo_scheme},
            {"config_radt_cld_scheme", &pr.cfg.config_radt_cld_scheme},
            {"config_radt_lw_scheme", &pr.cfg.config_radt_lw_scheme},
            {"config_radt_sw_scheme", &pr.cfg.config_radt_sw_scheme},
            {"config_sfclayer_scheme", &pr.cfg.config_sfclayer_scheme}};
        if (auto it = resolved.find(na.name); it != resolved.end()) v = *it->second;
        stream_opt.config_atts.push_back({na.name, {na.type, v}});
      }
      // MG:615-645: the physics diagnostics alarm from the 'output' and 'diagnostics' streams' output intervals
      const OutputStream* out = sm.find("output");
      const OutputStream* dia = sm.find("diagnostics");
      if (!out) throw Error("run --streams: MPAS's physics_run_init reads the 'output' stream's interval (MG:615)");
      long long step = 0;
      if (out->spec.output_interval != "none") {
        step = interval_seconds(out->spec.output_interval);
        if (dia && dia->spec.output_interval != "none")
          step = std::min(step, interval_seconds(dia->spec.output_interval));
      } else if (dia && dia->spec.output_interval != "none") {
        step = interval_seconds(dia->spec.output_interval);
      }
      if (step > 0) {
        pr.diag_alarm.defined = true;
        pr.diag_alarm.ring = 0;
        pr.diag_alarm.interval = step;
      }
      // the clock runs from the start to the stop time (mpas_atm_core.F:338-357, 723)
      const auto eff = effective_namelist(a.namelist);
      if (eff.at("config_stop_time") != "none")
        throw Error("run --streams: config_stop_time has no verified port (config_run_duration does)");
      const long long dur = interval_seconds(eff.at("config_run_duration"));
      const long long dts = static_cast<long long>(dcfg.config_dt);
      if (static_cast<real>(dts) != dcfg.config_dt || dur % dts != 0)
        throw Error("run --streams: config_run_duration is not a whole number of time steps");
      if (a.steps >= 0 && a.steps != dur / dts)
        throw Error("run --streams: --steps differs from config_run_duration / config_dt");
      a.steps = dur / dts;
    }
    // the diagnostics manager (core_atmosphere/diagnostics/mpas_atm_diagnostics.hpp). Without --streams a --write-at
    // file stands for an output stream whose alarm rings, its --write-fields for the stream's fields.
    diag::DiagnosticsContext dc;
    bool writing_state = false;
    long long clock_s = 0;  // the clock's time, seconds after the start
    dc.mesh = &in.mesh;
    dc.state = &in.state;
    dc.diag = &in.diag;
    dc.pools = &in.pools;
    dc.field_will_be_written = [&](const std::string& n) {
      if (streams_mode) return sm.field_will_be_written(n, clock_s);
      return writing_state && std::find(a.write_fields.begin(), a.write_fields.end(), n) != a.write_fields.end();
    };
    dc.output_stream_inclusion_count = [&](const std::string& n) {
      if (streams_mode) return sm.output_inclusion_count(n);
      return static_cast<int>(std::find(a.write_fields.begin(), a.write_fields.end(), n) != a.write_fields.end());
    };
    if (physics) diag::mpas_atm_diag_setup(dc, pr.cfg.config_sounding_interval);
    const long long dt_s = static_cast<long long>(dcfg.config_dt);
    auto write_step = [&](long n) {
      clock_s = n * dt_s;
      const bool step_file = a.every > 0 && n % a.every == 0;
      const bool state_file = std::find(a.write_at.begin(), a.write_at.end(), n) != a.write_at.end();
      const bool ringing = streams_mode && sm.any_ringing(clock_s);
      writing_state = state_file;
      // mpas_atm_core.F:686 (t = 0) and :848: only when an output stream writes
      if (step_file || state_file || ringing) atm_compute_output_diagnostics(in.mesh, in.state, in.diag);
      if (physics) {  // mpas_atm_core.F:692-694 (t = 0: reset first) and :855-856
        if (n == 0) diag::mpas_atm_diag_reset(dc);
        diag::mpas_atm_diag_update(dc);
        diag::mpas_atm_diag_compute(dc);
      }
      std::string ts = in.state.xtime.substr(0, 19);
      const std::string stamp = ts;
      for (auto& c : ts)
        if (c == ':') c = '.';
      HistoryOptions opt;
      opt.file_id = a.file_id;
      if (step_file) write_history(in, a.pos[1] + "/step." + ts + ".nc", fields, opt);
      if (state_file) write_history(in, a.pos[1] + "/state." + ts + ".nc", a.write_fields, opt);
      if (ringing) {
        // mpas_stream_mgr_write (mpas_atm_core.F:701, 862): every active output stream whose alarm rings
        for (const auto& o : sm.outputs)
          if (o.ringing(clock_s))
            write_history(in, a.pos[1] + "/" + expand_filename(o.spec.filename_template, stamp), o.fields, stream_opt);
        if (n > 0) {
          // mpas_atm_core.F:884-895 atm_reset_diagnostics (DO_PHYSICS): refl10cm_1km_max(:) = 0 when allocated
          const OutputStream* dia = sm.find("diagnostics");
          const RegistryVar& rv = registry_var("refl10cm_1km_max");
          if (dia && dia->ringing(clock_s) && in.pools.packages.active(rv)) {
            auto& r = in.pools.R("refl10cm_1km_max");
            for (std::size_t k = 0; k < r.size(); ++k) r.data()[k] = 0.;
          }
          // mpas_atm_core.F:899-908: the restart_timestamp file after a restart write (list-directed: a blank first)
          const OutputStream* rs = sm.find("restart");
          if (rs && rs->ringing(clock_s)) {
            std::FILE* tf = std::fopen((a.pos[1] + "/restart_timestamp").c_str(), "w");
            if (!tf) throw Error("run: cannot write restart_timestamp");
            std::fprintf(tf, " %s\n", stamp.c_str());
            std::fclose(tf);
          }
        }
      }
      if (physics) diag::mpas_atm_diag_reset(dc);  // mpas_atm_core.F:717, 914
      if (ringing)  // mpas_atm_core.F:719, 917
        for (auto& o : sm.outputs)
          if (o.ringing(clock_s)) o.reset(clock_s);
    };
    write_step(0);
    for (long n = 1; n <= a.steps; ++n) {
      itimestep = static_cast<int>(n);  // atm_core_run's itimestep (mpas_atm_core.F:722, 832)
      if (physics) {  // mpas_atm_core.F:1031-1037 (xtime_s = s + s_n / s_d of now - start)
        const long long now_s = (n - 1) * dt_s;
        const std::string now = add_seconds(start, now_s);
        phys::physics_timetracker(pr, in.mesh, in.pools, now, now_s, dcfg.config_dt, log);
        phys::physics_driver(pr, in.mesh, in.state, in.diag, in.pools, itimestep, static_cast<real>(now_s), pa);
        for (const auto& l : log) std::cerr << "MPAS log: " << l << "\n";
        log.clear();
      }
      atm_timestep(dcfg, in.mesh, in.state, in.diag, f, dcfg.config_dt, add_seconds(start, n * dt_s));
      shift_time_levels(in.state, f);  // mpas_atm_core.F:823
      write_step(n);
    }
  } else if (cmd == "device") {
    if (!a.pos.empty()) throw Error("usage: cppmpas device");
#if defined(CPPMPAS_USE_CUDA)
    std::cout << "backend cuda: " << gpu::describe() << "\n";
#else
    std::cout << "backend cpu\n";
#endif
    // a kernel writes the architecture it was compiled for (0 for the host): the backend that runs the kernels
    FArray<int> arch({64}, -1);
    const auto av = arch.view();
    parallel_for(1, 64, [=] CPPMPAS_LAMBDA(int i) {
#if defined(__CUDA_ARCH__)
      av(i) = __CUDA_ARCH__;
#else
      av(i) = 0;
#endif
    });
    for (int i = 1; i <= 64; ++i)
      if (arch(i) != arch(1)) throw Error("device check: the kernel's threads disagree");
    if (arch(1) < 0) throw Error("device check: the kernel did not run");
    std::cout << "kernels run on " << (arch(1) > 0 ? "the GPU, sm_" + std::to_string(arch(1) / 10) : std::string("the CPU"))
              << "\n";
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
#if defined(CPPMPAS_USE_CUDA)
  {
    std::string why;
    if (!gpu::init(why)) {
      std::cerr << "cppmpas_gpu: SKIP: no usable CUDA device (" << why << ")\n";
      return gpu::kNoDeviceExit;
    }
  }
#endif
  try {
    return run(argc, argv);
  } catch (const std::exception& e) {
    std::cerr << "ERROR: " << e.what() << "\n";
    return 1;
  }
}
