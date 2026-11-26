// SPDX-License-Identifier: GPL-3.0-or-later
// Unit tests of src/core/field.hpp, src/core/mesh.hpp, src/io/netcdf_file.cpp and the Registry table:
//  * FView offsets in 64 bits at the sizes of the NCAR x1.65536002 mesh (nCells 65,536,002; nEdges 196,608,000);
//  * an empty ElementOrder has 0 elements;
//  * NcFile writes one record only of a variable whose first dimension is the record dimension, creates files
//    without the fill pass, and reports a failing query instead of answering "absent";
//  * Registry defaults: a <var> with missing_value and no default_value takes the missing value
//    (src/tools/registry/gen_inc.c:1709-1715);
//  * the history writer refuses a character field that was not read (a grid file has no xtime).
//
// usage: test_io_core <scratch dir> <grid file>
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

#include "core/field.hpp"
#include "core/mesh.hpp"
#include "io/mpas_reader.hpp"
#include "io/mpas_writer.hpp"
#include "io/netcdf_file.hpp"
#include "io/registry.hpp"

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
    return std::string(e.what()).find(text) != std::string::npos;
  }
  return false;
}
template <class F> static bool throws(F f) {
  try {
    f();
  } catch (const Error&) {
    return true;
  }
  return false;
}

int main(int argc, char** argv) {
  if (argc != 3) {
    std::printf("usage: test_io_core <scratch dir> <grid file>\n");
    return 2;
  }
  const std::string dir = argv[1];

  // --- FView: u(nVertLevels = 55, nEdges + 1), zb(nVertLevelsP1 = 56, 2, nEdges + 1) on x1.65536002
  {
    volatile int nVertLevels = 55, nEdges = 196608000;  // runtime values: no constant folding of the products
    FView<double> u{nullptr, nVertLevels, nEdges + 1, 1};
    const std::int64_t want_u = 54 + std::int64_t{55} * (std::int64_t{nEdges} - 1);
    expect(u.offset(55, nEdges) == want_u, "FView offset u(55, nEdges) in 64 bits");
    FView<double> zb{nullptr, nVertLevels + 1, 2, nEdges + 1};
    const std::int64_t want_zb = 55 + std::int64_t{56} * (1 + std::int64_t{2} * (std::int64_t{nEdges} - 1));
    expect(zb.offset(56, 2, nEdges) == want_zb, "FView offset zb(56, 2, nEdges) in 64 bits");
    FView<double> s{nullptr, 3, nVertLevels, 65536003};
    const std::int64_t want_s = 2 + std::int64_t{3} * (54 + std::int64_t{55} * (std::int64_t{65536002} - 1));
    expect(s.offset(3, 55, 65536002) == want_s, "FView offset scalars(3, 55, nCells) in 64 bits");
    // small arrays: same element as FArray
    FArray<double> a({4, 3, 2}, 0.0);
    for (std::size_t k = 0; k < a.size(); ++k) a.data()[k] = static_cast<double>(k);
    expect(&a.view()(3, 2, 2) == &a(3, 2, 2), "FView and FArray address the same element");
  }

  // --- ElementOrder
  expect(ElementOrder{}.n() == 0, "empty ElementOrder has 0 elements");

  // --- NcFile
  {
    const std::string path = dir + "/test_io_core.nc";
    NcFile f = NcFile::create(path, NC_64BIT_DATA);
    f.def_dim("Time", NC_UNLIMITED);
    f.def_dim("nCells", 3);
    f.def_var("scalar", NC_DOUBLE, {});
    f.def_var("fixed", NC_DOUBLE, {"nCells"});
    f.def_var("rec", NC_DOUBLE, {"Time", "nCells"});
    f.def_var("unwritten", NC_DOUBLE, {"nCells"});
    f.enddef();
    expect(throws([&] { f.put_double("scalar", {1.0}, 0); }), "record write of a scalar variable is an error");
    expect(throws([&] { f.put_double("fixed", {1.0, 2.0, 3.0}, 0); }),
           "record write of a variable without record dimension is an error");
    expect(!throws([&] { f.put_double("fixed", {1.0, 2.0, 3.0}); }), "whole write of a fixed variable");
    expect(!throws([&] { f.put_double("rec", {1.0, 2.0, 3.0}, 0); }), "record 0 of a record variable");
    expect(!throws([&] { f.put_double("rec", {4.0, 5.0, 6.0}, 1); }), "record 1 of a record variable");
    expect(!throws([&] { f.put_double("scalar", {7.0}); }), "whole write of a scalar variable");
    f.close();
    expect(throws([&] { (void)f.has_dim("nCells"); }), "has_dim on a closed file is an error, not 'absent'");
    expect(throws([&] { (void)f.has_var("rec"); }), "has_var on a closed file is an error, not 'absent'");
    expect(throws([&] { (void)f.has_att("", "x"); }), "has_att on a closed file is an error, not 'absent'");

    const NcFile r = NcFile::open_read(path);
    expect(r.get_double("fixed") == std::vector<double>({1.0, 2.0, 3.0}), "fixed values");
    expect(r.get_double("rec", 1) == std::vector<double>({4.0, 5.0, 6.0}), "record 1 values");
    expect(r.get_double("scalar") == std::vector<double>({7.0}), "scalar value");
    expect(r.has_dim("nCells") && !r.has_dim("nope") && r.has_var("rec") && !r.has_var("nope"), "has_dim / has_var");
    expect(!r.has_att("", "nope") && !r.has_att("nope", "x") && !r.has_att("rec", "x"), "has_att absent");
    // no fill pass: a variable never written does not hold NC_FILL_DOUBLE
    bool filled = false;
    for (double v : r.get_double("unwritten")) filled |= (v == NC_FILL_DOUBLE);
    expect(!filled, "files are created with NC_NOFILL");
    std::remove(path.c_str());
  }

  // --- Registry defaults (gen_inc.c:1709-1715)
  expect(registry_var("tvxy").default_value == "-9999.0", "tvxy default = its missing_value -9999.0");
  expect(registry_var("isnowxy").default_value == "-9999", "isnowxy default = its missing_value -9999");
  expect(registry_var("latCell").default_value.empty(), "latCell has no default (0)");
  expect(registry_var("specZoneMaskCell").default_value == "0.0", "specZoneMaskCell default_value 0.0");

  // --- history writer: xtime of a grid file was never read
  {
    MpasInput in = read_mpas_grid(argv[2]);
    const std::string out = dir + "/test_io_core_xtime.nc";
    expect(throws([&] { write_history(in, out, {"latCell", "xtime"}, {"v8.4.2", 1, "0123456789"}); }),
           "writing an unread character field is an error");
    expect(!throws([&] { write_history(in, out, {"latCell"}, {"v8.4.2", 1, "0123456789"}); }),
           "writing a field that was read");
    std::remove(out.c_str());
    // no partial file: nprocs is checked before the file is created; an error after it (here NetCDF's, a name
    // defined twice) removes the file
    auto exists = [](const std::string& p) {
      std::FILE* f = std::fopen(p.c_str(), "r");
      if (f) std::fclose(f);
      return f != nullptr;
    };
    expect(throws_with([&] { write_history(in, out, {"latCell"}, {"v8.4.2", 10, "0123456789"}); }, "nprocs outside 1..9") &&
               !exists(out),
           "nprocs = 10 is refused before the file is created");
    expect(throws([&] { write_history(in, out, {"latCell", "latCell"}, {"v8.4.2", 1, "0123456789"}); }) && !exists(out),
           "a NetCDF error while writing leaves no file");
    std::remove(out.c_str());
    // each of the two guards against a name the Registry does not define (write_history meets registry_var first)
    const std::string undefined = "is not defined in MPAS-Atmosphere's Registry";
    expect(throws_with([&] { (void)registry_var("notavar"); }, undefined), "registry_var of an undefined name");
    expect(throws_with([&] { (void)field_ref(in, "notavar"); }, undefined), "field_ref of an undefined name");
    expect(throws_with([&] { (void)field_ref(in, "relhum"); }, "has no CPPMPAS port"), "field_ref of an unported name");
  }

  if (fails == 0) std::printf("OK io/core unit test\n");
  return fails ? 1 : 0;
}
