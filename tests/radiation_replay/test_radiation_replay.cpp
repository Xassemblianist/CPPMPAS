// SPDX-License-Identifier: GPL-3.0-or-later
// The radiation step of a model run, replayed: reads a dump of driver_radiation_lw (physics/mpas_atmphys_radiation.hpp
// write_radiation_dump: the step's inputs and the outputs the model computed), runs the step again on the dump's
// inputs (o3climatology_from_MPAS, driver_cloudiness, driver_radiation_sw if the step ran it, driver_radiation_lw)
// and compares every output pool array with the dump's, bit for bit (tolerance 0, -0.0 /= +0.0, NaN payloads).
// This compares CPPMPAS with CPPMPAS (the CPU run's own outputs): it proves the replay reproduces the run, and the
// CUDA build of it that the device reproduces the CPU. It is not a comparison with MPAS: the MPAS proof of the
// radiation is the RG tests (dycore_RG_*, dycore_RG_streams_*, against mpas842-dbl-rg-004) and the radiation_column
// check (MPAS's own radiation code built as a column program).
// The output arrays start as a NaN of their own payload, so an element the step leaves unwritten shows.
// Built twice from this source: test_radiation_replay (the CPU library: the replay reproduces the CPU model run) and,
// with CPPMPAS_ENABLE_CUDA, test_gpu_radiation_replay (the CUDA library: every loop of the radiation drivers and
// every RRTMG stage a kernel on the device, against the CPU model's outputs). Without a usable device the CUDA
// executable exits 77 (skipped).
//
// usage: test_radiation_replay <dump> <rrtmg_lw_tables.bin> <rrtmg_sw_tables.bin>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <string>
#include <vector>

#include "core/backend.hpp"
#include "core/error.hpp"
#include "core/field.hpp"
#include "physics/mpas_atmphys_interface.hpp"
#include "physics/mpas_atmphys_radiation.hpp"

namespace ph = cppmpas::phys;
using cppmpas::FArray;
using cppmpas::real;

namespace {

struct Rec {
  int kind = 0;
  std::vector<double> r;
  std::vector<std::int32_t> iv;
  std::string s;
};

std::map<std::string, Rec> read_dump(const char* path) {
  std::FILE* f = std::fopen(path, "rb");
  if (!f) throw cppmpas::Error(std::string("cannot open ") + path);
  std::map<std::string, Rec> recs;
  for (;;) {
    char nm[32];
    const std::size_t got = std::fread(nm, 1, 32, f);
    if (got == 0 && std::feof(f)) break;
    if (got != 32) throw cppmpas::Error("short record name (a truncated dump)");
    nm[31] = '\0';
    std::int32_t kind = 0;
    std::int64_t n = 0;
    if (std::fread(&kind, 4, 1, f) != 1 || std::fread(&n, 8, 1, f) != 1 || n < 0)
      throw cppmpas::Error("short record header");
    Rec r;
    r.kind = kind;
    const std::size_t sn = static_cast<std::size_t>(n);
    bool ok = true;
    if (kind == 0) {
      r.r.resize(sn);
      ok = std::fread(r.r.data(), 8, sn, f) == sn;
    } else if (kind == 1) {
      r.iv.resize(sn);
      ok = std::fread(r.iv.data(), 4, sn, f) == sn;
    } else {
      r.s.resize(sn);
      ok = std::fread(r.s.data(), 1, sn, f) == sn;
    }
    if (!ok) throw cppmpas::Error(std::string("short record ") + nm);
    recs[nm] = std::move(r);
  }
  std::fclose(f);
  return recs;
}

// a quiet NaN of a payload no computation produces
double sentinel() {
  const std::uint64_t b = 0x7ff8dead0000beefULL;
  double x;
  std::memcpy(&x, &b, 8);
  return x;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc != 4) {
    std::fprintf(stderr, "usage: %s <dump> <rrtmg_lw_tables.bin> <rrtmg_sw_tables.bin>\n", argv[0]);
    return 2;
  }
  unsetenv("CPPMPAS_RADIATION_DUMP");  // the replay writes no dump of its own
  try {
#if defined(CPPMPAS_USE_CUDA)
    std::string why;
    if (!cppmpas::gpu::init(why)) {
      std::printf("SKIP: no usable CUDA device (%s)\n", why.c_str());
      return cppmpas::gpu::kNoDeviceExit;
    }
    std::printf("device: %s\n", cppmpas::gpu::describe().c_str());
#endif
    std::map<std::string, Rec> d = read_dump(argv[1]);
    auto rec = [&](const char* name) -> const Rec& {
      auto it = d.find(name);
      if (it == d.end()) throw cppmpas::Error(std::string("the dump has no record ") + name);
      return it->second;
    };
    const auto& iv = rec("ints").iv;
    const auto& rv = rec("reals").r;
    const int ite = iv[0], nV = iv[1], noz = iv[2], nMonths = iv[3];
    if (ite < 1 || nV < 1) throw cppmpas::Error("the dump holds no cell or no level: nothing to compare");
    const bool sw_ran = iv[10] != 0, o3c = iv[11] != 0;
    const int ime = iv[13], kme = iv[14], nm = iv[15], np = iv[16];
    const real xtime_s = rv[0];
    const double nan = sentinel();

    // the pools: cells 1..ite and MPAS's garbage element, outputs filled with the sentinel
    std::map<std::string, FArray<real>> pool;
    auto cells_in = [&](const char* name) {
      FArray<real>& a = pool[name];
      a.allocate({ite + 1}, nan);
      const auto& x = rec(name).r;
      for (int i = 1; i <= ite; ++i) a(i) = x[std::size_t(i - 1)];
      return a.cview();
    };
    auto out1 = [&](const char* name) {
      FArray<real>& a = pool[name];
      a.allocate({ite + 1}, nan);
      return a.view();
    };
    auto out2 = [&](const char* name, int nk) {
      FArray<real>& a = pool[name];
      a.allocate({nk, ite + 1}, nan);
      return a.view();
    };
    ph::RadiationPools m;
    m.nVertLevels = nV;
    m.nOznLevels = noz;
    m.nMonths = nMonths;
    m.radt_lw_scheme = rec("radt_lw_scheme").s;
    m.radt_sw_scheme = rec("radt_sw_scheme").s;
    m.radt_cld_scheme = rec("radt_cld_scheme").s;
    m.radt_cld_overlap = rec("radt_cld_overlap").s;
    m.radt_cld_dcorrlen = rec("radt_cld_dcorrlen").s;
    m.microp_scheme = rec("microp_scheme").s;
    m.config_o3climatology = o3c;
    m.config_microp_re = iv[12] != 0;
    m.latCell = cells_in("latCell");
    m.lonCell = cells_in("lonCell");
    m.skintemp = cells_in("skintemp");
    m.snow = cells_in("snow");
    m.xice = cells_in("xice");
    m.xland = cells_in("xland");
    m.sfc_albedo = cells_in("sfc_albedo");
    m.sfc_emiss = cells_in("sfc_emiss");
    {
      FArray<real>& pin = pool["pin"];
      pin.allocate({noz}, nan);
      FArray<real>& oz = pool["ozmixm"];
      oz.allocate({nMonths, noz, ite + 1}, nan);  // the months the step does not read stay NaN
      if (o3c) {
        const auto &p = rec("pin").r, &om = rec("ozmixm_nm").r, &op = rec("ozmixm_np").r;
        for (int k = 1; k <= noz; ++k) pin(k) = p[std::size_t(k - 1)];
        for (int i = 1; i <= ite; ++i)
          for (int k = 1; k <= noz; ++k) {
            oz(nm, k, i) = om[std::size_t(k - 1) + std::size_t(noz) * (i - 1)];
            oz(np, k, i) = op[std::size_t(k - 1) + std::size_t(noz) * (i - 1)];
          }
      }
      m.pin = pin.view();
      m.ozmixm = oz.view();
    }
    m.o3clim = out2("o3clim", noz);
    m.o3vmr = out2("o3vmr", nV);
    m.cldfrac = out2("cldfrac", nV);
    m.rre_cloud = out2("rre_cloud", nV);
    m.rre_ice = out2("rre_ice", nV);
    m.rre_snow = out2("rre_snow", nV);
    m.rthratenlw = out2("rthratenlw", nV);
    m.rthratensw = out2("rthratensw", nV);
#define CPPMPAS_O(x) m.x = out1(#x);
    CPPMPAS_O(glw) CPPMPAS_O(lwcf) CPPMPAS_O(lwdnb) CPPMPAS_O(lwdnbc) CPPMPAS_O(lwdnt) CPPMPAS_O(lwdntc)
    CPPMPAS_O(lwupb) CPPMPAS_O(lwupbc) CPPMPAS_O(lwupt) CPPMPAS_O(lwuptc) CPPMPAS_O(olrtoa)
    CPPMPAS_O(coszr) CPPMPAS_O(gsw) CPPMPAS_O(swcf) CPPMPAS_O(swdnb) CPPMPAS_O(swdnbc) CPPMPAS_O(swdnt)
    CPPMPAS_O(swdntc) CPPMPAS_O(swupb) CPPMPAS_O(swupbc) CPPMPAS_O(swupt) CPPMPAS_O(swuptc) CPPMPAS_O(swddir)
    CPPMPAS_O(swddni) CPPMPAS_O(swddif)
#undef CPPMPAS_O

    // MPAS_to_physics' arrays the radiation reads
    ph::PhysicsVars v;
    v.ime = ime;
    v.kme = kme;
    auto setv = [&](ph::PArray& a, const char* name) {
      a.allocate(ime, kme);
      const auto& x = rec(name).r;
      if (x.size() != a.a.size()) throw cppmpas::Error(std::string("record ") + name + " has the wrong size");
      for (std::size_t e = 0; e < x.size(); ++e) a.a[e] = x[e];
    };
    setv(v.pres_hyd_p, "@pres_hyd_p"); setv(v.pres2_hyd_p, "@pres2_hyd_p"); setv(v.pi_p, "@pi_p");
    setv(v.t_p, "@t_p"); setv(v.t2_p, "@t2_p"); setv(v.dz_p, "@dz_p"); setv(v.pres_p, "@pres_p");
    setv(v.qv_p, "@qv_p"); setv(v.qc_p, "@qc_p"); setv(v.qi_p, "@qi_p"); setv(v.qs_p, "@qs_p");

    ph::RadiationVars r;
    r.julday = iv[4];
    r.year = iv[5];
    r.icloud = iv[6];
    r.has_reqc = iv[7];
    r.has_reqi = iv[8];
    r.has_reqs = iv[9];
    r.gmt = rv[1];
    r.curr_julday = rv[2];
    r.dt_radtsw = rv[3];
    r.dt_radtlw = rv[4];
    ph::RrtmgTables t;
    ph::init_radiation_rrtmg(argv[2], argv[3], t);

    // the radiation step as physics_timetracker (MG:286-296) and physics_driver (DR:210-246) run it
    if (o3c) ph::o3climatology_from_MPAS(r.curr_julday, m, ite);
    r.sw_ran = false;
    ph::driver_cloudiness(m, v, r, 1, ite);
    if (sw_ran) ph::driver_radiation_sw(1, m, v, r, t, xtime_s, 1, ite);
    ph::driver_radiation_lw(xtime_s, m, v, r, t, 1, ite);

    // the comparison
    std::vector<std::string> names = {"cldfrac", "rre_cloud", "rre_ice", "rre_snow", "rthratenlw", "glw", "lwcf",
                                      "lwdnb", "lwdnbc", "lwdnt", "lwdntc", "lwupb", "lwupbc", "lwupt", "lwuptc",
                                      "olrtoa"};
    if (o3c) names.insert(names.begin(), {"o3clim", "o3vmr"});
    if (sw_ran)
      for (const char* x : {"rthratensw", "coszr", "gsw", "swcf", "swdnb", "swdnbc", "swdnt", "swdntc", "swupb",
                            "swupbc", "swupt", "swuptc", "swddir", "swddni", "swddif"})
        names.push_back(x);
    std::printf("replay of %s: %d cells, %d levels, xtime %.0f s, %s\n", argv[1], ite, nV, xtime_s,
                sw_ran ? "LW and SW" : "LW only");
    int bad = 0;
    for (const std::string& name : names) {
      const auto& ref = rec(name.c_str()).r;
      const FArray<real>& a = pool.at(name);
      const int nk = a.rank() == 1 ? 1 : a.extent(0);
      if (ref.size() != std::size_t(nk) * std::size_t(ite)) {
        std::printf("  %-11s the dump has %zu values, expected %d\n", name.c_str(), ref.size(), nk * ite);
        ++bad;
        continue;
      }
      std::size_t ndiff = 0, first = 0;
      for (std::size_t e = 0; e < ref.size(); ++e) {
        if (std::memcmp(&ref[e], a.data() + e, 8) != 0) {
          if (ndiff == 0) first = e;
          ++ndiff;
        }
      }
      if (ndiff == 0) {
        std::printf("  %-11s %8zu values  IDENTICAL\n", name.c_str(), ref.size());
      } else {
        const std::size_t i = first / std::size_t(nk) + 1, k = first % std::size_t(nk) + 1;
        std::printf("  %-11s %8zu values  DIFFER: %zu, first at cell %zu level %zu: dump %.17g, replay %.17g\n",
                    name.c_str(), ref.size(), ndiff, i, k, ref[first], a.data()[first]);
        ++bad;
      }
    }
    std::printf("%zu fields compared, %d differ\n", names.size(), bad);
    return bad == 0 ? 0 : 1;
  } catch (const std::exception& e) {
    std::fprintf(stderr, "test_radiation_replay: %s\n", e.what());
    return 1;
  }
}
