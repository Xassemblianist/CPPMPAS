// SPDX-License-Identifier: GPL-3.0-or-later
// Unit tests of src/core_atmosphere/atm_microphysics.cpp and src/physics/mp_kessler.hpp that the SC step files cannot
// show:
//  * microphysics_init on the SC reference: Kessler active, dt_dyn = dt_microp = 3 s, mp_top_level = nVertLevels
//    (mpas_atmphys_init.F:276-282); a model top above config_microphysics_top (mp_top_level < nVertLevels, the
//    untaken PI:938-946) is refused; the MW namelist (microphysics off) leaves it inactive;
//  * the MPAS parameter es_cap of the Kessler kernel (module_mp_kessler.F:219 es = min(es, 0.99*pressure)): it
//    changes a column where es > 0.99*pressure and leaves an ordinary column unchanged. No SC column reaches the cap
//    (es is a few kPa at most against >= 5 kPa of pressure), so the step files cannot show it;
//  * driver_microphysics refuses a scheme's constituent whose index is MPAS's -1 (inactive, scalar_index), for WSM6
//    and Thompson, before it reads any array.
//
// usage: test_microphysics <SC init file> <SC namelist.atmosphere> <MW namelist.atmosphere>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

#include "core/error.hpp"
#include "core_atmosphere/atm_init.hpp"
#include "core_atmosphere/atm_microphysics.hpp"
#include "io/mpas_reader.hpp"
#include "physics/mp_kessler.hpp"
#include "physics/mp_thompson.hpp"
#include "physics/mpas_atmphys_microphysics.hpp"

using namespace cppmpas;

static int fails = 0;
static void expect(bool ok, const std::string& what) {
  if (!ok) {
    std::printf("FAIL %s\n", what.c_str());
    ++fails;
  }
}

// one column of nz levels through the kernel; returns qv, qc, qr, t after the call
// bit for bit (-0.0 differs from +0.0; a NaN equals only the same NaN)
static bool same_bits(const std::vector<real>& a, const std::vector<real>& b) {
  return a.size() == b.size() && (a.empty() || std::memcmp(a.data(), b.data(), a.size() * sizeof(real)) == 0);
}
static bool all_finite(const std::vector<real>& a) {
  for (real x : a)
    if (!std::isfinite(x)) return false;
  return !a.empty();
}

static std::vector<real> column(bool es_cap, real pii, real theta, real qv0) {
  const int nz = 4;
  std::vector<real> t(nz, theta), qv(nz, qv0), qc(nz, 1.e-4), qr(nz, 1.e-4), rho(nz, 0.1), pi(nz, pii), dz(nz, 500.0),
      z(nz), rainnc(1, 0.0), rainncv(1, 0.0);
  for (int k = 0; k < nz; ++k) z[k] = 250.0 + 500.0 * k;
  using phys::kessler::View2;
  using phys::kessler::View3;
  phys::kessler::KesslerKernel kk{};
  kk.t = View3<real>{t.data(), 0, 1};
  kk.qv = View3<real>{qv.data(), 0, 1};
  kk.qc = View3<real>{qc.data(), 0, 1};
  kk.qr = View3<real>{qr.data(), 0, 1};
  kk.rho = View3<const real>{rho.data(), 0, 1};
  kk.pii = View3<const real>{pi.data(), 0, 1};
  kk.dz8w = View3<const real>{dz.data(), 0, 1};
  kk.z = View3<const real>{z.data(), 0, 1};
  kk.rainnc = View2<real>{rainnc.data(), 0};
  kk.rainncv = View2<real>{rainncv.data(), 0};
  kk.dt = 3.0;
  kk.xlv = 2.50e6;
  kk.cp = 7.0 * 287.0 / 2.0;
  kk.ep2 = 287.0 / 461.6;
  kk.svp1 = 0.6112;
  kk.svp2 = 17.67;
  kk.svp3 = 29.65;
  kk.svpt0 = 273.15;
  kk.rhowater = 1000.0;
  kk.nz = nz;
  kk.i0 = 0;
  kk.ni = 1;
  kk.es_cap = es_cap;
  kk(0, 0);
  std::vector<real> out;
  for (int k = 0; k < nz; ++k) {
    out.push_back(qv[k]);
    out.push_back(qc[k]);
    out.push_back(qr[k]);
    out.push_back(t[k]);
  }
  return out;
}

int main(int argc, char** argv) {
  if (argc != 4) {
    std::printf("usage: test_microphysics <SC init file> <SC namelist.atmosphere> <MW namelist.atmosphere>\n");
    return 2;
  }
  InitConfig icfg = read_init_config(argv[2]);
  MpasInput in = read_mpas_init(argv[1], {"qv", "qc", "qr"});
  std::vector<std::string> log;
  atm_core_init(in.mesh, in.state, in.diag, icfg, log);

  MicrophysicsFields mp;
  microphysics_init(argv[2], in.mesh, 3.0, mp);
  expect(mp.active, "SC: Kessler active");
  expect(mp.mp_top_level == in.mesh.nVertLevels, "SC: mp_top_level = nVertLevels (" + std::to_string(mp.mp_top_level) + ")");
  expect(mp.dt_dyn == 3.0 && mp.dt_microp == 3.0, "SC: dt_dyn = dt_microp = 3");
  expect(mp.i_rainnc.allocated() && mp.i_rainnc(1) == 0 && mp.rainnc(1) == 0.0 && !std::signbit(mp.rainnc(1)),
         "SC: i_rainnc, rainnc +0.0");
  {  // a model top above config_microphysics_top = 45 km: refused
    Mesh m = in.mesh;
    for (int k = 1; k <= m.nVertLevels; ++k) m.dzu(k) = 3.0 * m.dzu(k);
    MicrophysicsFields mp2;
    bool threw = false;
    try {
      microphysics_init(argv[2], m, 3.0, mp2);
    } catch (const Error& e) {
      threw = std::strstr(e.what(), "mp_top_level") != nullptr;
    }
    expect(threw, "refuses mp_top_level < nVertLevels");
  }
  {
    MicrophysicsFields mp3;
    microphysics_init(argv[3], in.mesh, 6.0, mp3);
    expect(!mp3.active, "MW: microphysics off");
  }
  {  // a constituent the scheme needs is inactive (index -1): refused
    auto refused = [](auto run) {
      try {
        run();
      } catch (const Error& e) {
        return std::strstr(e.what(), "needs") != nullptr;
      }
      return false;
    };
    phys::MicrophysicsPools p;
    p.nVertLevels = 3;
    p.mp_top_level = 3;
    p.index_qv = 1; p.index_qc = 2; p.index_qr = 3; p.index_qi = 4; p.index_qs = 5; p.index_qg = 6;
    p.index_ni = 7; p.index_nr = 8;
    static phys::thompson::ThompsonView tv{};
    for (int* idx : {&p.index_qg, &p.index_ni}) {
      phys::MicrophysicsPools q = p;
      *(idx == &p.index_qg ? &q.index_qg : &q.index_ni) = -1;
      expect(refused([&] { phys::driver_microphysics_thompson(tv, q, 1, 1); }),
             std::string("Thompson refuses index ") + (idx == &p.index_qg ? "qg" : "ni") + " = -1");
    }
    phys::MicrophysicsPools w = p;
    w.index_qg = -1;
    expect(refused([&] { phys::driver_microphysics_wsm6(phys::wsm6::Wsm6Params<real>{}, w, 1, 1); }),
           "WSM6 refuses index qg = -1");
  }
  // es_cap: pii = 0.1 (p = 1e5 * 0.1^(1004/287) = 32 Pa) and theta = 3000 K (T = 300 K, es = 3.5 kPa): capped
  {
    const std::vector<real> capped = column(true, 0.1, 3000.0, 1.e-3), uncapped = column(false, 0.1, 3000.0, 1.e-3);
    expect(all_finite(capped) && all_finite(uncapped) && !same_bits(capped, uncapped),
           "es_cap changes a capped column (finite values)");
  }
  expect(same_bits(column(true, 0.95, 300.0, 1.e-2), column(false, 0.95, 300.0, 1.e-2)),
         "es_cap leaves an ordinary column (bit for bit)");
  if (fails) std::printf("%d failures\n", fails);
  else std::printf("all passed\n");
  return fails ? 1 : 0;
}
