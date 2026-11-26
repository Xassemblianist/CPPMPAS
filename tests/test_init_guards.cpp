// SPDX-License-Identifier: GPL-3.0-or-later
// Unit tests of the guards and atm_core_init steps of src/core_atmosphere/atm_init.cpp and
// src/operators/vector_operations.cpp on a reference init file:
//  * atm_initialize_deformation_weights refuses nEdgesOnCell > 24 (MPAS's cell_list is dimension(25)) before it
//    copies, and the periodic plane;
//  * mpas_initialize_vectors refuses a boundary edge (cellsOnEdge = nCells+1);
//  * mpas_atm_setup_bdy_masks refuses a limited-area mask;
//  * atm_core_init_len_disp (mpas_atm_core.F:155-201): both branches and the error;
//  * atm_start_time_stamp: YYYY-MM-DD_hh:mm:ss only;
//  * mpas_atm_run_compatibility (mpas_atm_core.F:1578-1617): the two limited-area errors and the warnings;
//  * atm_core_init sets xtime = config_start_time, blank-padded to StrLen (mpas_atm_core.F:261-262), whatever the
//    init file held.
//
// usage: test_init_guards <init file> <its namelist.atmosphere>
#include <cstdio>
#include <string>
#include <vector>

#include "core/error.hpp"
#include "core_atmosphere/atm_init.hpp"
#include "io/mpas_reader.hpp"
#include "operators/vector_operations.hpp"

using namespace cppmpas;

static int fails = 0;
static void expect(bool ok, const char* what) {
  if (!ok) {
    std::printf("FAIL %s\n", what);
    ++fails;
  }
}
template <class F> static bool throws_with(F f, const std::string& text) {
  try {
    f();
  } catch (const Error& e) {
    if (std::string(e.what()).find(text) != std::string::npos) return true;
    std::printf("  unexpected error: %s\n", e.what());
    return false;
  }
  return false;
}
static bool has(const std::vector<std::string>& log, const std::string& text) {
  for (const auto& l : log)
    if (l.find(text) != std::string::npos) return true;
  return false;
}

int main(int argc, char** argv) {
  if (argc != 3) {
    std::printf("usage: test_init_guards <init file> <namelist.atmosphere>\n");
    return 2;
  }
  const std::string path = argv[1];
  auto fresh = [&] {
    MpasInput in = read_mpas_init(path, {"qv"});
    allocate_init_fields(in.mesh, in.state, in.diag);
    return in;
  };

  // --- atm_initialize_deformation_weights
  {
    MpasInput in = fresh();
    Mesh& m = in.mesh;
    m.nEdgesOnCell(1) = 25;  // n = 26 > 25: the copy would write cell_list(26)
    expect(throws_with([&] { atm_initialize_deformation_weights(m, m.nCells, m.on_a_sphere, m.sphere_radius); },
                       "nEdgesOnCell > 24 overflows"),
           "nEdgesOnCell = 25 is refused");
    m.nEdgesOnCell(1) = 24;  // n = 25 passes the check (and then reads cellsOnCell past maxEdges: not run)
    MpasInput p = fresh();
    p.mesh.is_periodic = true;
    expect(throws_with([&] { atm_initialize_deformation_weights(p.mesh, p.mesh.nCells, p.mesh.on_a_sphere,
                                                                p.mesh.sphere_radius); },
                       "periodic-plane path"),
           "periodic mesh refused by atm_initialize_deformation_weights");
  }
  // --- mpas_initialize_vectors
  {
    MpasInput in = fresh();
    in.mesh.cellsOnEdge(2, 7) = in.mesh.nCells + 1;
    expect(throws_with([&] { mpas_initialize_vectors(in.mesh); }, "edge 7 is a boundary edge"),
           "boundary edge refused by mpas_initialize_vectors");
    MpasInput in1 = fresh();
    in1.mesh.cellsOnEdge(1, 9) = in1.mesh.nCells + 1;
    expect(throws_with([&] { mpas_initialize_vectors(in1.mesh); }, "edge 9 is a boundary edge"),
           "boundary edge (first cell) refused by mpas_initialize_vectors");
  }
  // --- mpas_atm_setup_bdy_masks
  for (int which = 0; which < 3; ++which) {
    MpasInput in = fresh();
    Mesh& m = in.mesh;
    (which == 0 ? m.bdyMaskCell(5) : which == 1 ? m.bdyMaskEdge(5) : m.bdyMaskVertex(5)) = 7;
    expect(throws_with([&] { mpas_atm_setup_bdy_masks(m); }, "limited-area mesh"),
           which == 0 ? "bdyMaskCell > 0 refused" : which == 1 ? "bdyMaskEdge > 0 refused" : "bdyMaskVertex > 0 refused");
  }
  // --- atm_core_init_len_disp
  {
    MpasInput in = fresh();
    InitConfig cfg;
    std::vector<std::string> log;
    in.mesh.nominalMinDc(1) = 0.0;
    cfg.config_len_disp = 577.35;
    atm_core_init_len_disp(in.mesh, cfg, log);
    expect(in.mesh.nominalMinDc(1) == 577.35, "nominalMinDc = config_len_disp when config_len_disp > 0");
    expect(has(log, "based on namelist option config_len_disp"), "log: nominalMinDc from config_len_disp");
    expect(!has(log, "WARNING"), "no warning when nominalMinDc was 0");
    log.clear();
    in.mesh.nominalMinDc(1) = 600.0;
    atm_core_init_len_disp(in.mesh, cfg, log);
    expect(in.mesh.nominalMinDc(1) == 577.35 && has(log, "WARNING: nominalMinDc was read from input file"),
           "a different positive nominalMinDc is overwritten with a warning");
    log.clear();
    in.mesh.nominalMinDc(1) = 1000.0;
    cfg.config_len_disp = 0.0;
    atm_core_init_len_disp(in.mesh, cfg, log);
    expect(in.mesh.nominalMinDc(1) == 1000.0 && cfg.config_len_disp == 1000.0,
           "config_len_disp = nominalMinDc when config_len_disp <= 0");
    in.mesh.nominalMinDc(1) = 0.0;
    cfg.config_len_disp = -1.0;
    expect(throws_with([&] { atm_core_init_len_disp(in.mesh, cfg, log); }, "Both config_len_disp and nominalMinDc"),
           "config_len_disp and nominalMinDc both <= 0 is an error");
  }
  // --- atm_start_time_stamp
  {
    InitConfig cfg;
    cfg.config_start_time = "0000-01-01_00:00:00";
    expect(atm_start_time_stamp(cfg) == "0000-01-01_00:00:00", "canonical start time");
    for (const char* bad : {"file", "0000-01-01", "0000-01-01_00:00:00.5", "2010-10-23 00:00:00"}) {
      cfg.config_start_time = bad;
      expect(throws_with([&] { (void)atm_start_time_stamp(cfg); }, "only YYYY-MM-DD_hh:mm:ss is handled"), bad);
    }
  }
  // --- mpas_atm_run_compatibility
  {
    MpasInput in = fresh();
    InitConfig cfg;
    std::vector<std::string> log;
    mpas_atm_run_compatibility(in.mesh, cfg, log);
    expect(has(log, " - Maximum value in bdyMaskCell = 0") && !has(log, "WARNING"), "global mesh, defaults: no warning");
    cfg.config_apply_lbcs = true;
    expect(throws_with([&] { mpas_atm_run_compatibility(in.mesh, cfg, log); }, "config_apply_lbcs = true, but no"),
           "config_apply_lbcs on a global mesh is an error");
    cfg.config_apply_lbcs = false;
    in.mesh.bdyMaskCell(in.mesh.nCells) = 1;
    expect(throws_with([&] { mpas_atm_run_compatibility(in.mesh, cfg, log); }, "Boundary cells found"),
           "boundary cells without config_apply_lbcs is an error");
    in.mesh.bdyMaskCell(in.mesh.nCells) = 0;
    log.clear();
    cfg.config_epssm = 0.1;
    cfg.config_positive_definite = true;
    mpas_atm_run_compatibility(in.mesh, cfg, log);
    expect(has(log, "config_epssm in the &nhyd_model namelist group is no longer supported"), "config_epssm warning");
    expect(has(log, "config_positive_definite = true will enable monotonic advection"), "positive definite warning");
  }

  // --- atm_core_init: xtime
  {
    MpasInput in = read_mpas_init(path, {"qv"});
    InitConfig cfg = read_init_config(argv[2]);
    in.state.xtime = "1999-12-31_23:59:59";
    std::vector<std::string> log;
    atm_core_init(in.mesh, in.state, in.diag, cfg, log);
    std::string want = cfg.config_start_time;
    want.resize(64, ' ');
    expect(in.state.xtime == want, "atm_core_init: xtime = config_start_time, blank-padded");
  }

  if (fails == 0) std::printf("OK init guards unit test\n");
  return fails ? 1 : 0;
}
