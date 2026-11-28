// SPDX-License-Identifier: GPL-3.0-or-later
// The C++ side of the MYNN surface layer column check: test_sfmynn_column <in.bin> <out.bin> runs init_sfclayer, then
// per step MPAS_to_physics and driver_sfclayer with config_sfclayer_scheme = 'sf_mynn' (src/physics) on the columns of
// in.bin in the order harness.F90 runs MPAS's own code, and writes the same records (tests/sfmynn_column/columns.py).
// config_pbl_scheme 'bl_ysu' as in harness.F90 (it only selects MPAS_to_physics' bl_mynn number concentrations).
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <map>
#include <string>
#include <vector>

#include "core/error.hpp"
#include "core/field.hpp"
#include "physics/mpas_atmphys_interface.hpp"
#include "physics/mpas_atmphys_sfclayer.hpp"

using cppmpas::FView;
using cppmpas::real;
namespace ph = cppmpas::phys;

namespace {

struct Rec {
  std::string name;
  std::int32_t kind = 0, nd = 0, dims[4] = {1, 1, 1, 1};
  std::vector<double> r;
  std::vector<std::int32_t> iv;
  std::size_t n() const {
    std::size_t s = 1;
    for (int d = 0; d < nd; ++d) s *= static_cast<std::size_t>(dims[d]);
    return s;
  }
};
std::vector<Rec> recs;
std::map<std::string, std::size_t> by_name;

Rec& rec(const std::string& n) {
  auto it = by_name.find(n);
  if (it == by_name.end()) throw cppmpas::Error("in.bin has no record " + n);
  return recs[it->second];
}
FView<real> R(const std::string& n) {
  Rec& r = rec(n);
  return FView<real>{r.r.data(), r.dims[0], r.nd > 1 ? r.dims[1] : 1, r.nd > 2 ? r.dims[2] : 1};
}
FView<const real> C(const std::string& n) {
  FView<real> v = R(n);
  return FView<const real>{v.p, v.e0, v.e1, v.e2};
}

void write_header(std::FILE* o, const std::string& name, std::int32_t kind, std::int32_t nd, const std::int32_t* dims) {
  char nm[32];
  std::memset(nm, ' ', 32);
  std::memcpy(nm, name.data(), name.size());
  std::fwrite(nm, 1, 32, o);
  std::int32_t meta[6] = {kind, nd, dims[0], dims[1], dims[2], dims[3]};
  std::fwrite(meta, 4, 6, o);
}
void put_r(std::FILE* o, const std::string& name, std::int32_t nd, std::vector<std::int32_t> dims, const double* a,
           std::size_t n) {
  dims.resize(4, 1);
  write_header(o, name, 0, nd, dims.data());
  std::fwrite(a, 8, n, o);
}
void marker(std::FILE* o, const char* name, int step) {
  const std::int32_t dims[4] = {1, 1, 1, 1};
  write_header(o, name, 1, 1, dims);
  const std::int32_t s = step;
  std::fwrite(&s, 4, 1, o);
}
void dump_all(std::FILE* o) {
  for (const Rec& r : recs) {
    if (r.name == "hdr" || r.name == "dtl") continue;
    write_header(o, r.name, r.kind, r.nd, r.dims);
    if (r.kind == 0)
      std::fwrite(r.r.data(), 8, r.r.size(), o);
    else
      std::fwrite(r.iv.data(), 4, r.iv.size(), o);
  }
}

}  // namespace

// test_sfmynn_column guards: the hard errors of the sf_mynn path. Each must throw an Error naming the reason.
int guards() {
  int fails = 0;
  auto expect_error = [&](const char* what, const char* needle, auto run) {
    bool ok = false;
    try {
      run();
    } catch (const cppmpas::Error& e) {
      ok = std::strstr(e.what(), needle) != nullptr;
      if (!ok) std::printf("%s: unexpected error: %s\n", what, e.what());
    }
    std::printf("%s %s\n", ok ? "ok  " : "FAIL", what);
    fails += !ok;
  };
  static ph::sfclayrev::SFClayRevTables<real> tables;
  expect_error("sf_mynn with one level", "fewer than 2 levels", [] {
    ph::SfclayerPools q;
    q.sfclayer_scheme = "sf_mynn";
    ph::PhysicsVars v;
    v.ime = 1;
    v.kme = 2;
    ph::driver_sfclayer(1, tables, q, v, 1, 1);
  });
  expect_error("sf_mynn without qcg, snowh, ch", "are not set", [] {
    ph::SfclayerPools q;
    q.sfclayer_scheme = "sf_mynn";
    ph::PhysicsVars v;
    v.ime = 1;
    v.kme = 56;
    ph::driver_sfclayer(1, tables, q, v, 1, 1);
  });
  return fails ? 1 : 0;
}

int main(int argc, char** argv) {
  if (argc == 2 && std::strcmp(argv[1], "guards") == 0) return guards();
  if (argc != 3) {
    std::fprintf(stderr, "usage: test_sfmynn_column <in.bin> <out.bin> | guards\n");
    return 2;
  }
  try {
    std::FILE* in = std::fopen(argv[1], "rb");
    if (!in) return 2;
    for (;;) {
      char nm[32];
      if (std::fread(nm, 1, 32, in) != 32) break;
      Rec r;
      r.name.assign(nm, 32);
      r.name.erase(r.name.find_last_not_of(' ') + 1);
      std::int32_t meta[6];
      if (std::fread(meta, 4, 6, in) != 6) throw cppmpas::Error("short in.bin");
      r.kind = meta[0];
      r.nd = meta[1];
      for (int d = 0; d < 4; ++d) r.dims[d] = meta[2 + d];
      if (r.kind == 0) {
        r.r.resize(r.n());
        if (std::fread(r.r.data(), 8, r.r.size(), in) != r.r.size()) throw cppmpas::Error("short in.bin");
      } else {
        r.iv.resize(r.n());
        if (std::fread(r.iv.data(), 4, r.iv.size(), in) != r.iv.size()) throw cppmpas::Error("short in.bin");
      }
      by_name[r.name] = recs.size();
      recs.push_back(std::move(r));
    }
    std::fclose(in);
    const auto& hdr = rec("hdr").iv;
    const int nC = hdr[0], nV = hdr[2], nsteps = hdr[3];
    const bool frac_seaice = hdr[4] == 1, do_restart = hdr[5] == 1;
    const real len_disp = rec("dtl").r[1];
    const std::string microp = "mp_wsm6", pbl = "bl_ysu", sfclayer = "sf_mynn";

    ph::ToPhysicsPools tp;
    tp.nVertLevels = nV;
    tp.microp_scheme = microp;
    tp.pbl_scheme = pbl;
    tp.zgrid = C("zgrid");
    tp.zz = C("zz");
    tp.rho_zz = C("rho_zz");
    tp.theta_m = C("theta_m");
    tp.w = C("w");
    tp.scalars = C("scalars");
    tp.index_qv = 1; tp.index_qc = 2; tp.index_qr = 3; tp.index_qi = 4; tp.index_qs = 5; tp.index_qg = 6;
    tp.surface_pressure = R("surface_pressure");
    tp.exner = C("exner");
    tp.pressure_base = C("pressure_base");
    tp.pressure_p = C("pressure_p");
    tp.uReconstructZonal = C("uReconstructZonal");
    tp.uReconstructMeridional = C("uReconstructMeridional");
    tp.plrad = R("plrad");

    ph::SfclayerPools sp;
    sp.sfclayer_scheme = sfclayer;
    sp.config_frac_seaice = frac_seaice;
    sp.config_do_restart = do_restart;
    sp.config_len_disp = len_disp;
    sp.meshDensity = C("meshDensity");
    sp.areaCell = C("areaCell");
    sp.skintemp = C("skintemp");
    sp.xland = C("xland");
    sp.sst = C("sst");
    sp.xice = C("xice");
    sp.snowh = C("snowh");
    sp.hpbl = C("hpbl");
    sp.mavail = C("mavail");
    sp.qcg = C("qcg");
    sp.br = R("br"); sp.cpm = R("cpm"); sp.chs = R("chs"); sp.chs2 = R("chs2"); sp.cqs2 = R("cqs2");
    sp.flhc = R("flhc"); sp.flqc = R("flqc"); sp.gz1oz0 = R("gz1oz0"); sp.hfx = R("hfx"); sp.qfx = R("qfx");
    sp.qgh = R("qgh"); sp.qsfc = R("qsfc"); sp.lh = R("lh"); sp.mol = R("mol"); sp.psim = R("psim");
    sp.psih = R("psih"); sp.regime = R("regime"); sp.rmol = R("rmol"); sp.ust = R("ust"); sp.ustm = R("ustm");
    sp.wspd = R("wspd"); sp.znt = R("znt"); sp.zol = R("zol"); sp.q2 = R("q2"); sp.t2m = R("t2m");
    sp.th2m = R("th2m"); sp.u10 = R("u10"); sp.v10 = R("v10"); sp.cd = R("cd"); sp.cda = R("cda");
    sp.ck = R("ck"); sp.cka = R("cka"); sp.fh = R("fh"); sp.fm = R("fm"); sp.ch = R("ch");
    sp.xice_threshold = frac_seaice ? 0.02 : 0.5;  // mpas_atmphys_manager.F:740-744

    ph::sfclayrev::SFClayRevTables<real> tables;
    ph::init_sfclayer(sfclayer, tables);

    std::FILE* out = std::fopen(argv[2], "wb");
    if (!out) return 2;
    put_r(out, "psim_stab", 1, {1001}, tables.psim_stab, 1001);
    put_r(out, "psim_unstab", 1, {1001}, tables.psim_unstab, 1001);
    put_r(out, "psih_stab", 1, {1001}, tables.psih_stab, 1001);
    put_r(out, "psih_unstab", 1, {1001}, tables.psih_unstab, 1001);
    for (int step = 1; step <= nsteps; ++step) {
      ph::PhysicsVars v;
      v.allocate_forall_physics(microp, pbl, nC, nV + 1);
      ph::MPAS_to_physics(tp, v, 1, nC);
      marker(out, "#A", step);
      dump_all(out);
      auto dumpv = [&](const char* name, const ph::PArray& a, int k0, int k1) {
        std::vector<double> buf;
        for (int k = k0; k <= k1; ++k)
          for (int i = 1; i <= nC; ++i) buf.push_back(a(i, k));
        put_r(out, name, 2, {nC, k1 - k0 + 1}, buf.data(), buf.size());
      };
      dumpv("@qv_p", v.qv_p, 1, nV); dumpv("@qc_p", v.qc_p, 1, nV);
      dumpv("@u_p", v.u_p, 1, nV); dumpv("@v_p", v.v_p, 1, nV);
      dumpv("@rho_p", v.rho_p, 1, nV); dumpv("@th_p", v.th_p, 1, nV); dumpv("@t_p", v.t_p, 1, nV);
      dumpv("@pi_p", v.pi_p, 1, nV); dumpv("@dz_p", v.dz_p, 1, nV);
      dumpv("@pres_hyd_p", v.pres_hyd_p, 1, nV);
      put_r(out, "@psfc_p", 1, {nC}, v.psfc_p.data(), nC);
      ph::SfclayerSeaArrays sea;
      ph::driver_sfclayer(step, tables, sp, v, 1, nC, &sea);
      marker(out, "#B", step);
      dump_all(out);
      if (frac_seaice) {
        put_r(out, "@chs_sea", 1, {nC}, sea.chs.data(), nC);
        put_r(out, "@chs2_sea", 1, {nC}, sea.chs2.data(), nC);
        put_r(out, "@cqs2_sea", 1, {nC}, sea.cqs2.data(), nC);
        put_r(out, "@cpm_sea", 1, {nC}, sea.cpm.data(), nC);
        put_r(out, "@hfx_sea", 1, {nC}, sea.hfx.data(), nC);
        put_r(out, "@qfx_sea", 1, {nC}, sea.qfx.data(), nC);
        put_r(out, "@qgh_sea", 1, {nC}, sea.qgh.data(), nC);
        put_r(out, "@qsfc_sea", 1, {nC}, sea.qsfc.data(), nC);
        put_r(out, "@lh_sea", 1, {nC}, sea.lh.data(), nC);
      }
    }
    // a short write (a full disk) is a failed run, never an output that looks complete
    const bool write_error = std::ferror(out) != 0;
    if (std::fclose(out) != 0 || write_error) {
      std::fprintf(stderr, "cannot write %s\n", argv[2]);
      return 1;
    }
  } catch (const std::exception& e) {
    std::fprintf(stderr, "%s\n", e.what());
    return 1;
  }
  return 0;
}
