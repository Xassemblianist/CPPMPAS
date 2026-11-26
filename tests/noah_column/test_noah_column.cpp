// SPDX-License-Identifier: GPL-3.0-or-later
// The C++ side of the Noah column check: test_noah_column <in.bin> <out.bin> <table dir> runs init_lsm, then per step
// driver_lsm, driver_seaice and atmphys_sfc_diagnostics (src/physics/mpas_atmphys_driver_lsm.cpp) on the points of
// in.bin as harness.F90 runs MPAS's own code, and writes the same out.bin layout (tests/noah_column/columns.py).
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

#include "core/field.hpp"
#include "physics/mpas_atmphys_driver_lsm.hpp"

using cppmpas::FView;
using cppmpas::real;

namespace {
constexpr int NR1 = 53, NS2 = 5, NA = 4, NSEA = 9;
const char* r1names[NR1] = {
    "acsnom", "acsnow", "canwat", "chs", "chs2", "chklowq", "cpm", "cqs2", "glw", "grdflx", "gsw", "hfx", "lai",
    "lh", "noahres", "potevp", "qfx", "qgh", "qsfc", "br", "sfc_albedo", "sfc_emibck", "sfc_emiss", "sfcrunoff",
    "smstav", "smstot", "snotime", "snopcx", "swddif", "swddir", "udrunoff", "z0", "znt", "sr", "rainncv", "raincv",
    "q2", "t2m", "th2m", "shdmin", "shdmax", "snoalb", "sfc_albbck", "snow", "snowc", "snowh", "tmn", "skintemp",
    "vegfra", "xice", "xland", "surface_pressure", "spare"};
}  // namespace

int main(int argc, char** argv) {
  if (argc != 4) {
    std::fprintf(stderr, "usage: test_noah_column <in.bin> <out.bin> <table dir>\n");
    return 2;
  }
  std::FILE* in = std::fopen(argv[1], "rb");
  if (!in) return 2;
  std::int32_t hdr[10];
  double dt;
  char lu[64];
  if (std::fread(hdr, 4, 10, in) != 10 || std::fread(&dt, 8, 1, in) != 1 || std::fread(lu, 1, 64, in) != 64) return 2;
  const int nC = hdr[0], nsoil = hdr[1], nsteps = hdr[2], nV = hdr[3];
  std::vector<int> ivgtyp(nC), isltyp(nC);
  std::vector<real> r1(std::size_t(nC) * NR1), s2(std::size_t(nsoil) * nC * NS2), a3(std::size_t(nC) * (nV + 1) * NA),
      sea(std::size_t(nC) * NSEA);
  auto rd = [&](void* p, std::size_t sz, std::size_t n) { return std::fread(p, sz, n, in) == n; };
  if (!rd(ivgtyp.data(), 4, nC) || !rd(isltyp.data(), 4, nC) || !rd(r1.data(), 8, r1.size()) ||
      !rd(s2.data(), 8, s2.size()) || !rd(a3.data(), 8, a3.size()) || !rd(sea.data(), 8, sea.size()))
    return 2;
  std::fclose(in);

  std::string mminlu(lu, 64);
  while (!mminlu.empty() && mminlu.back() == ' ') mminlu.pop_back();
  cppmpas::phys::LsmConfig cfg;
  cfg.config_sfc_albedo = hdr[6] != 0;
  cfg.config_sfc_snowalbedo = hdr[7] != 0;
  cfg.config_do_restart = hdr[8] != 0;

  cppmpas::phys::LsmPools p;
  p.nCells = nC; p.nCellsSolve = nC; p.nSoilLevels = nsoil;
  p.mminlu = mminlu; p.isice = hdr[5];
  p.ivgtyp = FView<int>{ivgtyp.data(), nC, 1, 1};
  p.isltyp = FView<int>{isltyp.data(), nC, 1, 1};
  auto f1 = [&](const char* name) {
    for (int k = 0; k < NR1; ++k)
      if (std::string(name) == r1names[k]) return FView<real>{r1.data() + std::size_t(k) * nC, nC, 1, 1};
    throw cppmpas::Error(std::string("no field ") + name);
  };
  auto f2 = [&](int k) { return FView<real>{s2.data() + std::size_t(k) * nsoil * nC, nsoil, nC, 1}; };
#define F1(x) p.x = f1(#x)
  F1(shdmin); F1(shdmax); F1(snoalb); F1(sfc_albbck); F1(snow); F1(snowc); F1(snowh); F1(tmn); F1(skintemp);
  F1(vegfra); F1(xice); F1(xland);
  F1(acsnom); F1(acsnow); F1(canwat); F1(chs); F1(chs2); F1(chklowq); F1(cpm); F1(cqs2); F1(glw); F1(grdflx);
  F1(gsw); F1(hfx); F1(lai); F1(lh); F1(noahres); F1(potevp); F1(qfx); F1(qgh); F1(qsfc); F1(br); F1(sfc_albedo);
  F1(sfc_emibck); F1(sfc_emiss); F1(sfcrunoff); F1(smstav); F1(smstot); F1(snotime); F1(snopcx); F1(swddif);
  F1(swddir); F1(udrunoff); F1(z0); F1(znt); F1(sr); F1(rainncv); F1(raincv); F1(q2); F1(t2m); F1(th2m);
  F1(surface_pressure);
#undef F1
  p.dzs = f2(0); p.sh2o = f2(1); p.smcrel = f2(2); p.smois = f2(3); p.tslb = f2(4);

  // physics_run_init (as harness.F90)
  cppmpas::phys::LsmVars v;
  v.ims = 1; v.ime = nC; v.kms = 1; v.kme = nV + 1;
  v.num_soils = nsoil; v.isurban = hdr[4]; v.dt_pbl = dt;
  v.xice_threshold = hdr[9] == 0 ? 0.5 : 0.02;
  auto fa = [&](int k) { return FView<const real>{a3.data() + std::size_t(k) * nC * (nV + 1), nC, nV + 1, 1}; };
  v.dz_p = fa(0); v.t_p = fa(1); v.qv_p = fa(2); v.pres2_hyd_p = fa(3);
  auto fs = [&](int k) { return FView<const real>{sea.data() + std::size_t(k) * nC, nC, 1, 1}; };
  v.cpm_sea = fs(0); v.chs_sea = fs(1); v.chs2_sea = fs(2); v.cqs2_sea = fs(3); v.hfx_sea = fs(4);
  v.qfx_sea = fs(5); v.qgh_sea = fs(6); v.qsfc_sea = fs(7); v.lh_sea = fs(8);

  std::FILE* out = std::fopen(argv[2], "wb");
  if (!out) return 2;
  try {
    cppmpas::phys::noah::NoahTables t;
    cppmpas::phys::init_lsm(cfg, p, argv[3], t);
    std::int32_t ti[66] = {t.lucats, t.bare, t.natural};
    for (int n = 0; n < 11; ++n) ti[3 + n] = t.lcz[n];
    ti[14] = t.slcats;
    ti[15] = t.slpcats;
    for (int n = 0; n < 50; ++n) ti[16 + n] = t.nrotbl[n];
    std::fwrite(ti, 4, 66, out);
    auto wa = [&](const real* a, int n) { std::fwrite(a, 8, n, out); };
    for (const real* a : {t.shdtbl, t.rstbl, t.rgltbl, t.hstbl, t.snuptbl, t.maxalb, t.laimintbl, t.laimaxtbl,
                          t.emissmintbl, t.emissmaxtbl, t.albedomintbl, t.albedomaxtbl, t.z0mintbl, t.z0maxtbl,
                          t.ztopvtbl, t.zbotvtbl})
      wa(a, 50);
    for (real x : {t.topt_data, t.cmcmax_data, t.cfactr_data, t.rsmax_data}) wa(&x, 1);
    for (const real* a : {t.bb, t.drysmc, t.f11, t.maxsmc, t.refsmc, t.satpsi, t.satdk, t.satdw, t.wltsmc, t.qtz})
      wa(a, 30);
    wa(t.slope_data, 30);
    for (real x : {t.sbeta_data, t.fxexp_data, t.csoil_data, t.salp_data, t.refdk_data, t.refkdt_data, t.frzk_data,
                   t.zbot_data, t.czil_data, t.smlow_data, t.smhigh_data, t.lvcoef_data})
      wa(&x, 1);
    wa(r1.data(), int(r1.size()));
    wa(s2.data(), int(s2.size()));
    for (int step = 1; step <= nsteps; ++step) {
      cppmpas::phys::driver_lsm(step, cfg, t, p, v, 1, nC);
      cppmpas::phys::driver_seaice(cfg, p, v, 1, nC);
      cppmpas::phys::atmphys_sfc_diagnostics(cfg, p, 1, nC);
      wa(r1.data(), int(r1.size()));
      wa(s2.data(), int(s2.size()));
    }
  } catch (const std::exception& e) {
    std::fprintf(stderr, "%s\n", e.what());
    return 1;
  }
  std::fclose(out);
  return 0;
}
