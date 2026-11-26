// SPDX-License-Identifier: GPL-3.0-or-later
// Unit tests of src/core_atmosphere/atm_dycore.cpp that the MW step files cannot show:
//  * check_dycore_supported accepts the MW reference configuration and refuses every configuration whose path is
//    not ported and verified (one option or mesh property changed at a time, each with its own message), including
//    values the namelist reader would already refuse (the dycore check must not rely on it);
//  * add_seconds: the clock's YYYY-MM-DD_hh:mm:ss arithmetic (gregorian, year 0 a leap year), and its refusals;
//  * read_dycore_config takes the options of the namelist file and config_len_disp / config_coef_3rd_order /
//    config_apvm_upwinding / config_dt from the InitConfig that atm_core_init left.
//
// usage: test_dycore <MW init file> <MW namelist.atmosphere>
#include <cmath>
#include <cstdio>
#include <functional>
#include <string>
#include <vector>

#include "core/error.hpp"
#include "core/namelist.hpp"
#include "core_atmosphere/atm_dycore.hpp"
#include "core_atmosphere/atm_init.hpp"
#include "io/mpas_reader.hpp"

using namespace cppmpas;

static int fails = 0;
static void expect(bool ok, const std::string& what) {
  if (!ok) {
    std::printf("FAIL %s\n", what.c_str());
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

int main(int argc, char** argv) {
  if (argc != 3) {
    std::printf("usage: test_dycore <MW init file> <MW namelist.atmosphere>\n");
    return 2;
  }
  InitConfig icfg = read_init_config(argv[2]);
  MpasInput in = read_mpas_init(argv[1], {"qv"});
  std::vector<std::string> log;
  atm_core_init(in.mesh, in.state, in.diag, icfg, log);
  const DycoreConfig mw = read_dycore_config(argv[2], icfg);

  // the options as the MW namelist (and the Registry defaults) give them
  expect(mw.config_dt == 6.0 && mw.config_number_of_sub_steps == 6 && mw.config_time_integration == "SRK3" &&
             mw.config_time_integration_order == 2 && !mw.config_split_dynamics_transport &&
             mw.config_horiz_mixing == "2d_fixed" && mw.config_h_theta_eddy_visc2 == 10.0 &&
             mw.config_v_mom_eddy_visc2 == 10.0 && mw.config_v_theta_eddy_visc2 == 10.0 &&
             mw.config_h_mom_eddy_visc4 == 0.0 && mw.config_mix_full && !mw.config_scalar_advection &&
             mw.config_smdiv == 0.1 && mw.config_len_disp == 577.35 && mw.config_coef_3rd_order == 0.25 &&
             mw.config_apvm_upwinding == 0.5 && mw.config_IAU_option == "off" && mw.config_physics_suite == "none" &&
             mw.config_les_model == "none" && !mw.config_apply_lbcs && !mw.config_rayleigh_damp_u &&
             mw.config_mpas_cam_coef == 0.0 && !mw.config_mix_scalars && mw.config_smagorinsky_coef == 0.125 &&
             mw.config_visc4_2dsmag == 0.05 && mw.config_del4u_div_factor == 10.0,
         "read_dycore_config: MW values");
  {  // config_len_disp comes from the InitConfig atm_core_init updated, not from the file
    InitConfig c2 = icfg;
    c2.config_len_disp = 123.0;
    expect(read_dycore_config(argv[2], c2).config_len_disp == 123.0, "read_dycore_config: len_disp from InitConfig");
  }
  try {
    check_dycore_supported(mw, in.mesh, in.state);
  } catch (const Error& e) {
    expect(false, std::string("MW configuration refused: ") + e.what());
  }

  struct Case {
    const char* text;
    std::function<void(DycoreConfig&)> change;
  };
  const std::vector<Case> cases = {
      {"config_time_integration = 'RK4'", [](DycoreConfig& c) { c.config_time_integration = "RK4"; }},
      {"config_time_integration_order = 3", [](DycoreConfig& c) { c.config_time_integration_order = 3; }},
      {"config_split_dynamics_transport = true", [](DycoreConfig& c) { c.config_split_dynamics_transport = true; }},
      {"config_scalar_advection = true", [](DycoreConfig& c) { c.config_scalar_advection = true; }},
      {"config_positive_definite = true", [](DycoreConfig& c) { c.config_positive_definite = true; }},
      {"config_horiz_mixing = '2d_smag'", [](DycoreConfig& c) { c.config_horiz_mixing = "2d_smag"; }},
      {"config_horiz_mixing = '2d_fixed '", [](DycoreConfig& c) { c.config_horiz_mixing = "2d_fixed "; }},
      {"config_les_model = '3d_smagorinsky'", [](DycoreConfig& c) { c.config_les_model = "3d_smagorinsky"; }},
      {"config_mpas_cam_coef > 0", [](DycoreConfig& c) { c.config_mpas_cam_coef = 0.2; }},
      {"config_rayleigh_damp_u = true", [](DycoreConfig& c) { c.config_rayleigh_damp_u = true; }},
      {"config_mix_scalars = true", [](DycoreConfig& c) { c.config_mix_scalars = true; }},
      {"config_mix_full = false", [](DycoreConfig& c) { c.config_mix_full = false; }},
      {"config_h_mom_eddy_visc4 > 0", [](DycoreConfig& c) { c.config_h_mom_eddy_visc4 = 1.0; }},
      {"config_h_theta_eddy_visc4 > 0", [](DycoreConfig& c) { c.config_h_theta_eddy_visc4 = 1.0; }},
      {"config_v_mom_eddy_visc2 < 0", [](DycoreConfig& c) { c.config_v_mom_eddy_visc2 = -1.0; }},
      {"config_v_mom_eddy_visc2 < 0", [](DycoreConfig& c) { c.config_v_mom_eddy_visc2 = std::nan(""); }},
      {"config_v_theta_eddy_visc2 < 0", [](DycoreConfig& c) { c.config_v_theta_eddy_visc2 = -1.0; }},
      {"config_apply_lbcs = true", [](DycoreConfig& c) { c.config_apply_lbcs = true; }},
      {"config_IAU_option = 'on'", [](DycoreConfig& c) { c.config_IAU_option = "on"; }},
      {"config_physics_suite = 'convection_permitting'",
       [](DycoreConfig& c) { c.config_physics_suite = "convection_permitting"; }},
      {"config_microp_scheme = 'mp_kessler'", [](DycoreConfig& c) { c.config_microp_scheme = "mp_kessler"; }},
      {"config_convection_scheme = 'cu_tiedtke'", [](DycoreConfig& c) { c.config_convection_scheme = "cu_tiedtke"; }},
      {"config_pbl_scheme = 'bl_ysu'", [](DycoreConfig& c) { c.config_pbl_scheme = "bl_ysu"; }},
      {"radiation schemes", [](DycoreConfig& c) { c.config_radt_lw_scheme = "rrtmg_lw"; }},
      {"radiation schemes", [](DycoreConfig& c) { c.config_radt_sw_scheme = "rrtmg_sw"; }},
      {"config_dt that is not a whole number", [](DycoreConfig& c) { c.config_dt = 6.5; }},
      {"config_dt that is not a whole number", [](DycoreConfig& c) { c.config_dt = 0.0; }},
      {"config_dt that is not a whole number", [](DycoreConfig& c) { c.config_dt = 86401.0; }},
      {"config_number_of_sub_steps < 2", [](DycoreConfig& c) { c.config_number_of_sub_steps = 1; }},
  };
  for (const auto& cs : cases) {
    DycoreConfig c = mw;
    cs.change(c);
    expect(throws_with([&] { check_dycore_supported(c, in.mesh, in.state); }, cs.text), std::string("refuses ") + cs.text);
  }
  // the JW options on the MW mesh: 2d_smagorinsky, no vertical mixing (v_*_eddy_visc2 = 0)
  {
    DycoreConfig c = mw;
    c.config_horiz_mixing = "2d_smagorinsky";
    c.config_v_mom_eddy_visc2 = c.config_v_theta_eddy_visc2 = 0.0;
    try {
      check_dycore_supported(c, in.mesh, in.state);
    } catch (const Error& e) {
      expect(false, std::string("JW mixing options refused: ") + e.what());
    }
  }
  // 'suite' means 'off' with the suite 'none' (mpas_atmphys_control.F:157-165)
  {
    DycoreConfig c = mw;
    c.config_microp_scheme = c.config_convection_scheme = c.config_pbl_scheme = "suite";
    c.config_radt_lw_scheme = c.config_radt_sw_scheme = "off";
    try {
      check_dycore_supported(c, in.mesh, in.state);
    } catch (const Error& e) {
      expect(false, std::string("'suite' / 'off' schemes refused: ") + e.what());
    }
  }
  // mesh and state properties
  {
    Mesh m = in.mesh;
    m.fEdge(17) = 1.0e-4;  // MW's u_init is not 0
    expect(throws_with([&] { check_dycore_supported(mw, m, in.state); }, "fEdge != 0 with u_init or v_init != 0"),
           "refuses fEdge != 0 with u_init != 0");
    {  // u_init = v_init = 0 (JW): the reference_u term is +-0, fEdge != 0 is taken
      State s0 = in.state;
      for (int k = 1; k <= m.nVertLevels; ++k) s0.u_init(k) = s0.v_init(k) = 0.0;
      try {
        check_dycore_supported(mw, m, s0);
      } catch (const Error& e) {
        expect(false, std::string("fEdge != 0 with u_init = v_init = 0 refused: ") + e.what());
      }
      s0.v_init(m.nVertLevels) = 1.0;  // the last level of v_init alone
      expect(throws_with([&] { check_dycore_supported(mw, m, s0); }, "fEdge != 0 with u_init or v_init != 0"),
             "refuses fEdge != 0 with v_init != 0");
    }
    m = in.mesh;
    m.nVertLevels = 69;
    expect(throws_with([&] { check_dycore_supported(mw, m, in.state); }, "nVertLevels = 69"), "refuses odd nVertLevels");
    m.nVertLevels = 2;
    expect(throws_with([&] { check_dycore_supported(mw, m, in.state); }, "nVertLevels = 2"), "refuses nVertLevels < 4");
    m = in.mesh;
    m.vertexDegree = 4;
    expect(throws_with([&] { check_dycore_supported(mw, m, in.state); }, "vertexDegree"), "refuses vertexDegree 4");
    State s = in.state;
    s.scalar_names = {"qv", "qc", "qr"};
    expect(throws_with([&] { check_dycore_supported(mw, in.mesh, s); }, "scalars other than qv"),
           "refuses qv, qc, qr");
  }

  // Time (F3 of the dycore review): time level 1 is the Time of the input stream (Registry.xml:615), time level 2
  // starts as its copy (mpas_atm_core.F:213-218), atm_timestep sets Time(2) = Time(1) + dt (TI:1844-1849) and the
  // shift makes it time level 1
  {
    MpasInput t = in;
    t.state.Time(1) = 1800.0;
    DycoreFields tf;
    allocate_dycore_fields(t.mesh, t.state, tf);
    expect(tf.Time2 == 1800.0, "Time: time level 2 starts as time level 1");
    atm_timestep(mw, t.mesh, t.state, t.diag, tf, 6.0, "0000-01-01_00:30:06");
    shift_time_levels(t.state, tf);
    expect(t.state.Time(1) == 1806.0 && tf.Time2 == 1800.0, "Time: 1800 + 6 after one step");
    State s0 = in.state;
    s0.Time = FArray<real>();
    DycoreFields f0;
    expect(throws_with([&] { allocate_dycore_fields(in.mesh, s0, f0); }, "the state Time was not read"),
           "allocate_dycore_fields refuses a state without Time");
    expect(in.state.Time.rank() == 1 && in.state.Time(1) == 0.0, "Time is read from the MW init file (0)");
  }

  // add_seconds
  expect(add_seconds("0000-01-01_00:00:00", 6) == "0000-01-01_00:00:06", "add_seconds 6");
  expect(add_seconds("0000-01-01_00:00:54", 6) == "0000-01-01_00:01:00", "add_seconds minute");
  expect(add_seconds("0000-01-01_00:59:59", 1) == "0000-01-01_01:00:00", "add_seconds hour");
  expect(add_seconds("0000-01-01_23:59:59", 1) == "0000-01-02_00:00:00", "add_seconds day");
  expect(add_seconds("0000-02-28_23:59:59", 1) == "0000-02-29_00:00:00", "add_seconds year 0 is a leap year");
  expect(add_seconds("1900-02-28_23:59:59", 1) == "1900-03-01_00:00:00", "add_seconds 1900 is not");
  expect(add_seconds("2000-02-28_12:00:00", 86400) == "2000-02-29_12:00:00", "add_seconds 2000 is");
  expect(add_seconds("2001-02-28_23:59:59", 1) == "2001-03-01_00:00:00", "add_seconds 2001 is not");
  expect(add_seconds("2010-10-23_00:00:00", 9 * 86400 + 3600) == "2010-11-01_01:00:00", "add_seconds 31-day month");
  expect(add_seconds("1999-12-31_23:59:59", 1) == "2000-01-01_00:00:00", "add_seconds year");
  expect(add_seconds("2024-04-30_00:00:00", 86400) == "2024-05-01_00:00:00", "add_seconds 30-day month");
  expect(add_seconds("0000-01-01_00:00:00", 0) == "0000-01-01_00:00:00", "add_seconds 0");
  expect(throws_with([] { add_seconds("0000-01-01_00:00", 6); }, "is not YYYY-MM-DD_hh:mm:ss"), "add_seconds short");
  expect(throws_with([] { add_seconds("0000-01-01_00:00:00x", 6); }, "is not YYYY-MM-DD_hh:mm:ss"), "add_seconds long");
  expect(throws_with([] { add_seconds("2001-02-29_00:00:00", 6); }, "is not a valid date"), "add_seconds Feb 29");
  expect(throws_with([] { add_seconds("2001-13-01_00:00:00", 6); }, "is not a valid date"), "add_seconds month 13");
  expect(throws_with([] { add_seconds("2001-01-01_24:00:00", 6); }, "is not a valid date"), "add_seconds hour 24");
  expect(throws_with([] { add_seconds("0000-01-01_00:00:00", -1); }, "negative"), "add_seconds negative");

  std::printf("%s (%d failures)\n", fails ? "FAILED" : "OK", fails);
  return fails ? 1 : 0;
}
