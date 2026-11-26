// SPDX-License-Identifier: GPL-3.0-or-later
// Unit tests of the guards and atm_core_init steps of src/core_atmosphere/atm_init.cpp and
// src/operators/vector_operations.cpp on a reference init file:
//  * atm_initialize_deformation_weights refuses nEdgesOnCell outside 1..maxEdges, and > 24 (MPAS's cell_list is
//    dimension(25)) before it copies, and the periodic plane;
//  * mpas_initialize_vectors refuses a boundary edge (cellsOnEdge = nCells+1);
//  * mpas_atm_setup_bdy_masks refuses a limited-area mask;
//  * atm_core_init_len_disp (mpas_atm_core.F:155-201): both branches and the error;
//  * atm_start_time_stamp: YYYY-MM-DD_hh:mm:ss only;
//  * mpas_atm_run_compatibility (mpas_atm_core.F:1578-1617): the two limited-area errors and the warnings;
//  * atm_compute_solve_diagnostics refuses config_apvm_upwinding <= 0 (no APVM, mpas_atm_time_integration.F:7372);
//  * init_block_verified_fields: everything for an even nVertLevels >= 4, not the pow-dependent fields of an odd
//    one or of nVertLevels <= 2;
//  * atm_core_init sets xtime = config_start_time, blank-padded to StrLen (mpas_atm_core.F:261-262), whatever the
//    init file held;
//  * atm_mpas_init_block with config_do_restart refuses state / diag fields the restart stream did not fill.
//
// usage: test_init_guards <init file> <its namelist.atmosphere>
#include <algorithm>
#include <cstdio>
#include <exception>
#include <string>
#include <vector>

#include "core/error.hpp"
#include "core/namelist.hpp"
#include "core_atmosphere/atm_init.hpp"
#include "core_atmosphere/atm_solve_diagnostics.hpp"
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
    // nEdgesOnCell outside 1..maxEdges would index past a column of the (maxEdges, nCells+1) arrays
    for (int bad : {0, -1, m.maxEdges + 1}) {
      m.nEdgesOnCell(1) = bad;
      expect(throws_with([&] { atm_initialize_deformation_weights(m, m.nCells, m.on_a_sphere, m.sphere_radius); },
                         "nEdgesOnCell(1) outside 1..maxEdges"),
             "nEdgesOnCell outside 1..maxEdges is refused");
    }
    // MPAS's dimension(25): a mesh with maxEdges = 25. cellsOnCell gets 25 rows; cell 1's added neighbours are
    // nCells+1, so the cell is skipped (:1710-1716) and no other array is read past its old maxEdges rows
    const int rows = m.maxEdges;
    m.maxEdges = 25;
    FArray<int> wide({25, m.nCells + 1}, m.nCells + 1);
    for (int i = 1; i <= m.nCells + 1; ++i)
      for (int j = 1; j <= rows; ++j) wide(j, i) = m.cellsOnCell(j, i);
    m.cellsOnCell = wide;
    m.nEdgesOnCell(1) = 25;  // n = 26 > 25: the copy would write cell_list(26)
    expect(throws_with([&] { atm_initialize_deformation_weights(m, m.nCells, m.on_a_sphere, m.sphere_radius); },
                       "nEdgesOnCell > 24 overflows"),
           "nEdgesOnCell = 25 is refused");
    m.nEdgesOnCell(1) = 24;  // n = 25: accepted, no error of any kind
    bool ok = true;
    try {
      atm_initialize_deformation_weights(m, m.nCells, m.on_a_sphere, m.sphere_radius);
    } catch (const std::exception& e) {
      std::printf("  unexpected error: %s\n", e.what());
      ok = false;
    }
    expect(ok, "nEdgesOnCell = 24 (n = 25) is accepted");
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
    cfg.config_start_time = "2010-12-31_23:59:59";
    expect(atm_start_time_stamp(cfg) == "2010-12-31_23:59:59", "last second of a year");
    cfg.config_start_time = "2010-02-28_00:00:00";
    expect(atm_start_time_stamp(cfg) == "2010-02-28_00:00:00", "February 28");
    for (const char* bad : {"file", "0000-01-01", "0000-01-01_00:00:00.5", "2010-10-23 00:00:00"}) {
      cfg.config_start_time = bad;
      expect(throws_with([&] { (void)atm_start_time_stamp(cfg); }, "only YYYY-MM-DD_hh:mm:ss is handled"), bad);
    }
    for (const char* bad : {"0000-13-01_00:00:00", "0000-00-01_00:00:00", "0000-01-00_00:00:00", "0000-01-32_00:00:00",
                            "0000-04-31_00:00:00", "2000-02-29_00:00:00", "0000-01-01_24:00:00", "0000-01-01_00:60:00",
                            "0000-01-01_00:00:60"}) {
      cfg.config_start_time = bad;
      expect(throws_with([&] { (void)atm_start_time_stamp(cfg); }, "a field out of its range is not handled"), bad);
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
    in.mesh.bdyMaskCell(in.mesh.nCells) = 1;
    expect(throws_with([&] { mpas_atm_run_compatibility(in.mesh, cfg, log); }, "lbc_in stream) is not ported"),
           "config_apply_lbcs with boundary cells: the lbc_in check is not ported");
    in.mesh.bdyMaskCell(in.mesh.nCells) = 0;
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
    MpasInput g = fresh();
    g.mesh.bdyMaskCell = FArray<int>();  // never allocated (a Mesh filled by a library caller)
    InitConfig c0;
    expect(throws_with([&] { mpas_atm_run_compatibility(g.mesh, c0, log); }, "bdyMaskCell was not read"),
           "mpas_atm_run_compatibility without bdyMaskCell is an error");
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

  // --- atm_core_init: config_len_disp = nominalMinDc reaches the caller's config (mpas_atm_core.F:190-192 sets the
  // pool's config_len_disp); a bad start time is refused before anything changes
  {
    MpasInput in = read_mpas_init(path, {"qv"});
    InitConfig cfg = read_init_config(argv[2]);
    in.mesh.nominalMinDc(1) = 600.25;
    cfg.config_len_disp = 0.0;
    std::vector<std::string> log;
    atm_core_init(in.mesh, in.state, in.diag, cfg, log);
    expect(cfg.config_len_disp == 600.25, "atm_core_init: the caller's config_len_disp = nominalMinDc");
    MpasInput b = read_mpas_init(path, {"qv"});
    InitConfig bad = read_init_config(argv[2]);
    bad.config_start_time = "0000-13-01_00:00:00";
    b.mesh.nominalMinDc(1) = 0.0;
    expect(throws_with([&] { atm_core_init(b.mesh, b.state, b.diag, bad, log); }, "a field out of its range"),
           "atm_core_init: a bad start time is an error");
    expect(b.mesh.nominalMinDc(1) == 0.0, "atm_core_init: nominalMinDc untouched when the start time is refused");
  }

  // --- the mesh scaling is unverified for a meshDensity other than 1 with config_h_ScaleWithMesh
  {
    MpasInput in = read_mpas_init(path, {"qv"});
    InitConfig cfg;
    for (int i = 1; i <= in.mesh.nCells; ++i) in.mesh.meshDensity(i) = 1.0;
    cfg.config_h_ScaleWithMesh = true;
    expect(!mesh_scaling_unverified(in.mesh, cfg), "meshDensity = 1: the mesh scaling is verified");
    in.mesh.meshDensity(in.mesh.nCells) = 2.0;
    expect(mesh_scaling_unverified(in.mesh, cfg), "meshDensity(nCells) = 2: the mesh scaling is not verified");
    cfg.config_h_ScaleWithMesh = false;
    expect(!mesh_scaling_unverified(in.mesh, cfg), "config_h_ScaleWithMesh = false: nothing is scaled");
    cfg.config_h_ScaleWithMesh = true;
    in.mesh.meshDensity = FArray<real>();
    expect(throws_with([&] { (void)mesh_scaling_unverified(in.mesh, cfg); }, "meshDensity was not read"),
           "mesh_scaling_unverified without meshDensity");
    const auto v = init_block_verified_fields(70, true);
    bool none = true;
    for (const auto& n : v) none = none && n.rfind("meshScaling", 0) != 0;
    expect(none && std::find(v.begin(), v.end(), "pv_edge") != v.end(),
           "verified fields of a scaled mesh: no meshScaling*, the others");
  }

  // --- mpas_legs: N < 1 is an error (the routine is public)
  {
    std::vector<real> A(4, 1.0), B(4, 1.0), X(4, 0.0);
    std::vector<int> I(4, 1);
    expect(throws_with([&] { mpas_legs(A.data() + 2, 0, B.data() + 2, X.data() + 2, I.data() + 2); }, "mpas_legs: N < 1"),
           "mpas_legs with N = 0 is an error");
  }

  // --- atm_compute_solve_diagnostics: the non-APVM branch is not ported
  {
    MpasInput in = fresh();
    InitConfig cfg = read_init_config(argv[2]);
    std::vector<std::string> log;
    atm_core_init(in.mesh, in.state, in.diag, cfg, log);
    SolveDiagWork w;
    allocate_solve_diag_work(in.mesh, w);
    const Mesh& m = in.mesh;
    expect(throws_with([&] { atm_compute_solve_diagnostics(cfg.config_dt, 0.0, m, in.state.rho_zz, in.state.u, in.diag,
                                                           w, 1, m.nCells, 1, m.nVertices, 1, m.nEdges); },
                       "config_apvm_upwinding <= 0"),
           "solve diagnostics: config_apvm_upwinding = 0 refused");
  }

  // --- verified fields
  {
    const auto& all = init_block_computed_fields();
    auto in_list = [](const std::vector<std::string>& l, const char* n) {
      for (const auto& x : l)
        if (x == n) return true;
      return false;
    };
    const auto even = init_block_verified_fields(70), odd = init_block_verified_fields(71);
    const auto scaled = init_block_verified_fields(70, true);
    const auto two = init_block_verified_fields(2);  // even, but the reference's scalar pow loop (nVertLevels <= 2)
    // exactly the computed fields, each once; minus the four pow-tail fields for an odd nVertLevels, minus the four
    // meshScaling* fields for a scaled mesh
    auto exactly = [&](const std::vector<std::string>& l, const std::vector<std::string>& minus) {
      if (l.size() + minus.size() != all.size()) return false;
      for (const auto& n : all) {
        const bool out = std::find(minus.begin(), minus.end(), n) != minus.end();
        if (in_list(l, n.c_str()) == out || std::count(l.begin(), l.end(), n) > 1) return false;
      }
      return true;
    };
    expect(exactly(even, {}), "verified fields: all computed fields for an even nVertLevels");
    expect(exactly(odd, {"exner", "exner_base", "pressure_p", "pressure_base"}),
           "verified fields: all but the pow tail (exner, exner_base, pressure_p, pressure_base) for an odd nVertLevels");
    expect(exactly(two, {"exner", "exner_base", "pressure_p", "pressure_base"}),
           "verified fields: all but the pow-dependent fields for nVertLevels = 2 (scalar pow loop, no reference)");
    expect(exactly(scaled, {"meshScalingDel2", "meshScalingDel4", "meshScalingRegionalCell", "meshScalingRegionalEdge"}),
           "verified fields: all but meshScaling* for a scaled mesh");
  }

  // --- atm_mpas_init_block with config_do_restart = true on an input-stream file: the coupled fields were never read
  {
    MpasInput in = read_mpas_init(path, {"qv"});
    InitConfig cfg = read_init_config(argv[2]);
    cfg.config_do_restart = true;
    std::vector<std::string> log;
    expect(throws_with([&] { atm_mpas_init_block(in.mesh, in.state, in.diag, cfg, log); },
                       "config_do_restart = true, but theta_m was not read from the restart stream"),
           "a restart init block refuses fields the restart stream did not fill");
  }

  if (fails == 0) std::printf("OK init guards unit test\n");
  return fails ? 1 : 0;
}
