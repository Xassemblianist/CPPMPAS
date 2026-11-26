// SPDX-License-Identifier: GPL-3.0-or-later
// Unit test of src/core/namelist.cpp on variants of a reference namelist.atmosphere: the values read_init_config
// gives (bit for bit, as gfortran's namelist read of MPAS gives them), and the forms that are errors.
//  * reals: d exponents, a value that underflows to +0.0 is the accepted 0.0 (gfortran reads it as +0.0), -0.0 and
//    one that underflows to -0.0 are not (the sign of zero reaches the results: dss = config_xnutr * ...);
//  * integers: a value outside a 32-bit integer is an error (gfortran's read fails, MPAS stops);
//  * logicals and names in any case; character values in ' or ", '' inside; at most StrKIND = 512 characters;
//  * a second group of the same name is an error (MPAS reads only the first, gen_inc.c:855-857);
//  * the line of an error is the line of the option's name;
//  * active_scalars of the accepted physics options.
//
// usage: test_namelist <reference namelist.atmosphere (MW)> <scratch dir>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iterator>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

#include "core/error.hpp"
#include "core/namelist.hpp"

using namespace cppmpas;

static int fails = 0;
static void expect(bool ok, const std::string& what) {
  if (!ok) {
    std::printf("FAIL %s\n", what.c_str());
    ++fails;
  }
}
static std::uint64_t bits(double x) {
  std::uint64_t b;
  std::memcpy(&b, &x, 8);
  return b;
}

int main(int argc, char** argv) {
  if (argc != 3) {
    std::printf("usage: test_namelist <namelist.atmosphere> <scratch dir>\n");
    return 2;
  }
  std::string ref;
  {
    std::ifstream in(argv[1]);
    std::stringstream ss;
    ss << in.rdbuf();
    ref = ss.str();
  }
  const std::string path = std::string(argv[2]) + "/test_namelist.nml";
  // the reference with one substitution (the pattern must match once)
  auto variant = [&](const std::string& pattern, const std::string& repl) {
    const std::regex re(pattern);
    const auto n = std::distance(std::sregex_iterator(ref.begin(), ref.end(), re), std::sregex_iterator());
    if (n != 1) {
      std::printf("test setup: %s matched %ld times\n", pattern.c_str(), static_cast<long>(n));
      std::exit(2);
    }
    std::ofstream(path) << std::regex_replace(ref, re, repl);
  };
  auto read = [&](const std::string& pattern, const std::string& repl) {
    variant(pattern, repl);
    return read_init_config(path);
  };
  auto error_of = [&](const std::string& pattern, const std::string& repl) -> std::string {
    variant(pattern, repl);
    try {
      (void)read_init_config(path);
    } catch (const Error& e) {
      return e.what();
    } catch (const std::exception& e) {
      return std::string("not a cppmpas::Error: ") + e.what();
    }
    return "(no error)";
  };
  auto refused = [&](const char* name, const std::string& pattern, const std::string& repl, const std::string& text) {
    const std::string e = error_of(pattern, repl);
    expect(e.find(text) != std::string::npos && e.rfind("not a cppmpas::Error", 0) != 0,
           std::string(name) + ": want '" + text + "', got: " + e);
  };

  // --- the reference itself
  {
    const InitConfig c = read("config_zd = 10500.0", "config_zd = 10500.0");
    expect(c.config_zd == 10500.0 && c.config_xnutr == 0.1 && c.config_dt == 6.0 && c.config_len_disp == 577.35 &&
               !c.config_h_ScaleWithMesh && c.config_start_time == "0000-01-01_00:00:00" &&
               c.config_microp_scheme == "suite" && c.config_physics_suite == "none",
           "reference values");
  }
  // --- reals
  expect(bits(read("config_zd = 10500.0", "config_zd = 1.05d4").config_zd) == bits(10500.0), "1.05d4 is 10500.0");
  expect(bits(read("config_zd = 10500.0", "config_zd = 105D2").config_zd) == bits(10500.0), "105D2 is 10500.0");
  {
    const InitConfig c = read("config_xnutr = 0.1", "config_xnutr = 1e-400");
    expect(bits(c.config_xnutr) == 0, "1e-400 underflows to +0.0, the accepted 0.0");
  }
  expect(bits(read("config_xnutr = 0.1", "config_xnutr = 0.0e0").config_xnutr) == 0, "0.0e0 is +0.0");
  refused("negative zero", "config_xnutr = 0.1", "config_xnutr = -0.0", "config_xnutr = -0.0 has no verified port");
  refused("underflow to negative zero", "config_xnutr = 0.1", "config_xnutr = -1e-400",
          "config_xnutr = -1e-400 has no verified port");
  refused("real overflow", "config_xnutr = 0.1", "config_xnutr = 1e999", "config_xnutr = 1e999 has no verified port");
  // the value on the next line (gfortran reads it); a comment between '=' and the value is refused
  expect(bits(read("config_zd = 10500.0", "config_zd =\n      10500.0").config_zd) == bits(10500.0),
         "value on the line after '='");
  refused("comment after '='", "config_zd = 10500.0", "config_zd = ! c\n 10500.0", "config_zd has no value");
  // --- integers
  refused("integer above 2^31-1", "config_number_of_sub_steps = 6", "config_number_of_sub_steps = 2147483648",
          "config_number_of_sub_steps = 2147483648 is not a 32-bit integer");
  refused("integer below -2^31", "config_number_of_sub_steps = 6", "config_number_of_sub_steps = -2147483649",
          "is not a 32-bit integer");
  refused("integer beyond 64 bits", "config_number_of_sub_steps = 6",
          "config_number_of_sub_steps = 99999999999999999999", "is not a 32-bit integer");
  refused("integer in range, not a reference value", "config_number_of_sub_steps = 6",
          "config_number_of_sub_steps = 2147483647", "config_number_of_sub_steps = 2147483647 has no verified port");
  // --- logicals, names, characters
  expect(read("config_h_ScaleWithMesh = false", "CONFIG_H_SCALEWITHMESH = .t.").config_h_ScaleWithMesh,
         "CONFIG_H_SCALEWITHMESH = .t. is true");
  expect(read("config_h_ScaleWithMesh = false", "config_h_scalewithmesh = T").config_h_ScaleWithMesh,
         "config_h_scalewithmesh = T is true");
  expect(!read("config_h_ScaleWithMesh = false", "config_h_ScaleWithMesh = .FALSE.").config_h_ScaleWithMesh,
         ".FALSE. is false");
  expect(read("config_start_time = '0000-01-01_00:00:00'", "config_start_time = \"0000-01-01_00:00:00\"")
                 .config_start_time == "0000-01-01_00:00:00",
         "a double-quoted character value");
  refused("'' inside quotes", "config_physics_suite = 'none'", "config_physics_suite = 'no''ne'",
          "config_physics_suite = 'no'ne' has no verified port");
  refused("longer than StrKIND", "config_physics_suite = 'none'",
          "config_physics_suite = '" + std::string(513, 'x') + "'", "is not at most StrKIND = 512 characters long");
  // --- groups
  {
    variant("&damping", "&damping\n    config_xnutr = 0.1\n/\n&damping");
    std::string e = "(no error)";
    try {
      (void)read_namelist_file(path);
    } catch (const Error& x) {
      e = x.what();
    }
    expect(e.find("group &damping appears twice") != std::string::npos, "a second &damping group: " + e);
  }
  // --- the line of an error is the line of the name (here 2 lines above its '=')
  {
    variant("    config_xnutr = 0.1\n", "    config_xnutr = 0.1\n    config_not_an_option\n\n = 1\n");
    std::string e;
    try {
      (void)read_namelist_file(path);
    } catch (const Error& x) {
      e = x.what();
    }
    int line = 0;
    {
      std::ifstream in(path);
      std::string l;
      for (int n = 1; std::getline(in, l); ++n)
        if (l.find("config_not_an_option") != std::string::npos) line = n;
    }
    expect(e.find(path + ":" + std::to_string(line) + ": config_not_an_option is not an option") != std::string::npos,
           "error line: " + e);
  }
  // --- active scalars
  {
    InitConfig c;
    c.config_physics_suite = "none";
    expect(active_scalars(c) == std::vector<std::string>{"qv"}, "microphysics 'suite' of the suite 'none': qv");
    c.config_microp_scheme = "mp_kessler";
    expect(active_scalars(c) == std::vector<std::string>{"qv", "qc", "qr"}, "mp_kessler: qv, qc, qr");
    c.config_microp_scheme = "mp_wsm6";
    bool threw = false;
    try {
      (void)active_scalars(c);
    } catch (const Error&) {
      threw = true;
    }
    expect(threw, "mp_wsm6 is not handled");
  }
  std::remove(path.c_str());
  if (fails == 0) std::printf("OK namelist unit test\n");
  return fails ? 1 : 0;
}
