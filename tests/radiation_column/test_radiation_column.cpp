// SPDX-License-Identifier: GPL-3.0-or-later
// The C++ side of the radiation column check: test_radiation_column <in.bin> <out.bin> <tables dir> <lw.bin> <sw.bin>
// runs init_o3climatology (OZONE_*.TBL of <tables dir>), the RRTMG tables (<lw.bin>, <sw.bin>: memory images),
// then per step o3climatology_from_MPAS, driver_cloudiness, driver_radiation_sw and driver_radiation_lw
// (src/physics) on the columns of in.bin in the order harness.F90 runs MPAS's own code, and writes the same records
// (tests/radiation_column/columns.py).
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <map>
#include <string>
#include <vector>

#include "core/error.hpp"
#include "core/field.hpp"
#include "physics/mpas_atmphys_interface.hpp"
#include "physics/mpas_atmphys_radiation.hpp"

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
void marker(std::FILE* o, const char* name, int step) {
  const std::int32_t dims[4] = {1, 1, 1, 1};
  write_header(o, name, 1, 1, dims);
  const std::int32_t s = step;
  std::fwrite(&s, 4, 1, o);
}
void dump_all(std::FILE* o) {
  for (const Rec& r : recs) {
    if (r.name == "hdr" || r.name == "dtl" || r.name == "time" || r.name[0] == '@') continue;
    write_header(o, r.name, r.kind, r.nd, r.dims);
    if (r.kind == 0)
      std::fwrite(r.r.data(), 8, r.r.size(), o);
    else
      std::fwrite(r.iv.data(), 4, r.iv.size(), o);
  }
}

}  // namespace

int main(int argc, char** argv) {
  if (argc != 6) {
    std::fprintf(stderr, "usage: test_radiation_column <in.bin> <out.bin> <tables dir> <lw.bin> <sw.bin>\n");
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
    const int nC = hdr[0], nV = hdr[1], nsteps = hdr[2];
    const char* overlaps[] = {"none", "random", "maximum_random", "maximum", "exponential", "exponential_random"};

    ph::RadiationPools m;
    m.nVertLevels = nV;
    m.nOznLevels = hdr[3];
    m.nMonths = hdr[4];
    m.radt_lw_scheme = "rrtmg_lw";
    m.radt_sw_scheme = "rrtmg_sw";
    m.radt_cld_scheme = "cld_fraction";
    m.radt_cld_overlap = overlaps[hdr[5]];
    m.radt_cld_dcorrlen = hdr[6] == 1 ? "latitude_varying" : "constant";
    m.microp_scheme = "mp_wsm6";
    m.config_o3climatology = hdr[7] == 1;
    m.config_microp_re = false;
    m.latCell = C("latCell"); m.lonCell = C("lonCell");
    m.skintemp = C("skintemp"); m.snow = C("snow"); m.xice = C("xice"); m.xland = C("xland");
    m.pin = R("pin"); m.ozmixm = R("ozmixm");
    m.sfc_albedo = C("sfc_albedo"); m.sfc_emiss = C("sfc_emiss");
    m.cldfrac = R("cldfrac"); m.o3clim = R("o3clim"); m.o3vmr = R("o3vmr");
    m.rre_cloud = R("rre_cloud"); m.rre_ice = R("rre_ice"); m.rre_snow = R("rre_snow");
    m.glw = R("glw"); m.lwcf = R("lwcf"); m.lwdnb = R("lwdnb"); m.lwdnbc = R("lwdnbc"); m.lwdnt = R("lwdnt");
    m.lwdntc = R("lwdntc"); m.lwupb = R("lwupb"); m.lwupbc = R("lwupbc"); m.lwupt = R("lwupt");
    m.lwuptc = R("lwuptc"); m.olrtoa = R("olrtoa");
    m.coszr = R("coszr"); m.gsw = R("gsw"); m.swcf = R("swcf"); m.swdnb = R("swdnb"); m.swdnbc = R("swdnbc");
    m.swdnt = R("swdnt"); m.swdntc = R("swdntc"); m.swupb = R("swupb"); m.swupbc = R("swupbc");
    m.swupt = R("swupt"); m.swuptc = R("swuptc"); m.swddir = R("swddir"); m.swddni = R("swddni");
    m.swddif = R("swddif");
    m.rthratenlw = R("rthratenlw"); m.rthratensw = R("rthratensw");

    // MPAS_to_physics' arrays (the "@" records)
    ph::PhysicsVars v;
    v.ime = nC;
    v.kme = nV + 1;
    auto setv = [&](ph::PArray& a, const char* name, int nk) {
      a.allocate(nC, nV + 1);
      const Rec& r = rec(name);
      for (int k = 1; k <= nk; ++k)
        for (int i = 1; i <= nC; ++i) a(i, k) = r.r[static_cast<std::size_t>(i - 1) + std::size_t(nC) * (k - 1)];
    };
    setv(v.pres_hyd_p, "@pres_hyd_p", nV); setv(v.pres2_hyd_p, "@pres2_hyd_p", nV + 1);
    setv(v.pi_p, "@pi_p", nV); setv(v.t_p, "@t_p", nV); setv(v.t2_p, "@t2_p", nV + 1);
    setv(v.dz_p, "@dz_p", nV); setv(v.pres_p, "@pres_p", nV);
    setv(v.qv_p, "@qv_p", nV); setv(v.qc_p, "@qc_p", nV); setv(v.qi_p, "@qi_p", nV); setv(v.qs_p, "@qs_p", nV);

    ph::RadiationVars r;
    r.dt_radtsw = rec("dtl").r[0];
    r.dt_radtlw = rec("dtl").r[1];
    ph::RrtmgTables t;
    ph::init_radiation_rrtmg(argv[4], argv[5], t);
    if (m.config_o3climatology) ph::init_o3climatology(argv[3], m, nC + 1);

    std::FILE* o = std::fopen(argv[2], "wb");
    if (!o) return 2;
    marker(o, "#I", 0);
    dump_all(o);
    const auto& tm = rec("time").r;
    for (int step = 1; step <= nsteps; ++step) {
      const real xtime_s = tm[std::size_t(4 * (step - 1))];
      r.julday = static_cast<int>(std::lround(tm[std::size_t(4 * (step - 1) + 1)]));
      r.curr_julday = tm[std::size_t(4 * (step - 1) + 2)];
      r.gmt = tm[std::size_t(4 * (step - 1) + 3)];
      r.year = 2026;
      r.sw_ran = false;
      if (m.config_o3climatology) ph::o3climatology_from_MPAS(r.curr_julday, m, nC);
      ph::driver_cloudiness(m, v, r, 1, nC);
      ph::driver_radiation_sw(step, m, v, r, t, xtime_s, 1, nC);
      ph::driver_radiation_lw(xtime_s, m, v, r, t, 1, nC);
      marker(o, "#S", step);
      dump_all(o);
    }
    std::fclose(o);
  } catch (const std::exception& e) {
    std::fprintf(stderr, "test_radiation_column: %s\n", e.what());
    return 1;
  }
  return 0;
}
