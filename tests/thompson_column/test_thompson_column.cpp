// SPDX-License-Identifier: GPL-3.0-or-later
// The C++ side of the Thompson column check: test_thompson_column <in.bin> <out.bin> <init.bin> <tables dir>
// [l_diags] runs init_thompson_clouddroplets_forMPAS, thompson_init (dump to init.bin) and driver_microphysics_thompson
// (src/physics/mpas_atmphys_microphysics_thompson.cpp) on the columns of in.bin as harness.F90 runs MPAS's own code,
// and writes the same out.bin layout (tests/thompson_column/columns.py).
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <utility>
#include <vector>

#include "core/field.hpp"
#include "physics/mp_thompson_init.hpp"
#include "physics/mpas_atmphys_constants.hpp"
#include "physics/mpas_atmphys_microphysics.hpp"

using cppmpas::FArray;
using cppmpas::real;

int main(int argc, char** argv) {
  if (argc != 5 && !(argc == 6 && std::strcmp(argv[5], "l_diags") == 0)) {
    std::fprintf(stderr, "usage: test_thompson_column <in.bin> <out.bin> <init.bin> <tables dir> [l_diags]\n");
    return 2;
  }
  const bool diags = argc == 6;
  std::FILE* in = std::fopen(argv[1], "rb");
  if (!in) return 2;
  std::int32_t hdr[5];
  double dt;
  if (std::fread(hdr, 4, 5, in) != 5 || std::fread(&dt, 8, 1, in) != 1) return 2;
  const int nC = hdr[0], nV = hdr[1], nsteps = hdr[2], nscal = hdr[4];
  auto rd = [&](FArray<real>& a, std::vector<int> ext) {
    a.allocate(ext);
    if (std::fread(a.data(), 8, a.size(), in) != a.size()) throw cppmpas::Error("short in.bin");
  };
  FArray<real> zgrid, zz, rho_zz, theta_m, w, scalars, exner, pb, pp, exb, rtb, rtp, sfc, rainprod, evapprod, rec, rei,
      res, rainnc, snownc, graupelnc, rland, nt_c, mu_c;
  rd(zgrid, {nV + 1, nC}); rd(zz, {nV, nC}); rd(rho_zz, {nV, nC}); rd(theta_m, {nV, nC}); rd(w, {nV + 1, nC});
  rd(scalars, {nscal, nV, nC});
  rd(exner, {nV, nC}); rd(pb, {nV, nC}); rd(pp, {nV, nC}); rd(exb, {nV, nC}); rd(rtb, {nV, nC}); rd(rtp, {nV, nC});
  rd(sfc, {nC});
  rd(rainprod, {nV, nC}); rd(evapprod, {nV, nC}); rd(rec, {nV, nC}); rd(rei, {nV, nC}); rd(res, {nV, nC});
  rd(rainnc, {nC}); rd(snownc, {nC}); rd(graupelnc, {nC}); rd(rland, {nC}); rd(nt_c, {nC}); rd(mu_c, {nC});
  std::fclose(in);
  FArray<int> landmask({nC});
  for (int i = 1; i <= nC; ++i) landmask(i) = static_cast<int>(std::lround(rland(i)));
  // harness.F90 sets these to -999 (fields the driver writes before it reads them)
  FArray<real> rainncv({nC}, -999.), snowncv({nC}, -999.), graupelncv({nC}, -999.), sr({nC}, -999.),
      precipw({nC}, -999.), tsp({nC}, -999.), rtd({nV, nC}, -999.), dth({nV, nC}, -999.);
  FArray<real> ten[9];
  for (auto& t : ten) t.allocate({nV, nC}, -999.);

  try {
    // mpas_atmphys_manager.F:686-694
    real dt_microp = 0;
    int n_microp = 0;
    cppmpas::phys::microp_timestep("mp_thompson", dt, dt_microp, n_microp);

    // init_microphysics (DRV:293-295)
    cppmpas::phys::thompson::ThompsonTables tables;
    cppmpas::phys::thompson::thompson_init(argv[4], tables);
    cppmpas::phys::thompson::dump_thompson_init(tables, argv[3]);
    cppmpas::phys::init_thompson_clouddroplets_forMPAS(nC, std::as_const(landmask).view(), nt_c.view(), mu_c.view());

    cppmpas::phys::MicrophysicsPools m;
    m.nVertLevels = nV;
    m.zgrid = std::as_const(zgrid).view(); m.zz = std::as_const(zz).view(); m.rho_zz = rho_zz.view();
    m.theta_m = theta_m.view(); m.scalars = scalars.view(); m.w = std::as_const(w).view();
    m.index_qv = 1; m.index_qc = 2; m.index_qr = 3; m.index_qi = 4; m.index_qs = 5; m.index_qg = 6;
    m.index_ni = 7; m.index_nr = 8;
    m.exner = exner.view(); m.pressure_p = pp.view(); m.rtheta_p = rtp.view(); m.surface_pressure = sfc.view();
    m.pressure_base = std::as_const(pb).view(); m.exner_base = std::as_const(exb).view();
    m.rtheta_base = std::as_const(rtb).view();
    m.rt_diabatic_tend = rtd.view(); m.tend_sfc_pressure = tsp.view(); m.dtheta_dt_mp = dth.view();
    m.precipw = precipw.view(); m.rainnc = rainnc.view(); m.rainncv = rainncv.view(); m.snownc = snownc.view();
    m.snowncv = snowncv.view(); m.graupelnc = graupelnc.view(); m.graupelncv = graupelncv.view(); m.sr = sr.view();
    m.rainprod = rainprod.view(); m.evapprod = evapprod.view(); m.re_cloud = rec.view(); m.re_ice = rei.view();
    m.re_snow = res.view();
    m.rthmpten = ten[0].view(); m.rqvmpten = ten[1].view(); m.rqcmpten = ten[2].view(); m.rqrmpten = ten[3].view();
    m.rqimpten = ten[4].view(); m.rqsmpten = ten[5].view(); m.rqgmpten = ten[6].view();
    m.rnimpten = ten[7].view(); m.rnrmpten = ten[8].view();
    m.nt_c = std::as_const(nt_c).view(); m.mu_c = std::as_const(mu_c).view();
    m.mp_top_level = hdr[3];
    m.dt_dyn = dt;
    m.dt_microp = dt_microp;
    m.n_microp = n_microp;
    FArray<real> refl10cm({nV, nC}, -999.), refl_max({nC}, -999.), refl_1km({nC}, -999.), refl_1km_max({nC}, 3.),
        relhum({nV, nC}, -999.);
    m.l_diags = diags;
    m.relhum = relhum.view(); m.refl10cm = refl10cm.view(); m.refl10cm_max = refl_max.view();
    m.refl10cm_1km = refl_1km.view(); m.refl10cm_1km_max = refl_1km_max.view();

    std::FILE* out = std::fopen(argv[2], "wb");
    if (!out) return 2;
    auto wr = [&](const FArray<real>& a) { std::fwrite(a.data(), 8, a.size(), out); };
    wr(nt_c); wr(mu_c);
    const double nm = n_microp;
    std::fwrite(&dt_microp, 8, 1, out);
    std::fwrite(&nm, 8, 1, out);
    for (int step = 1; step <= nsteps; ++step) {
      cppmpas::phys::driver_microphysics_thompson(tables.v, m, 1, nC);
      wr(scalars); wr(theta_m); wr(rtp); wr(exner); wr(pp); wr(sfc); wr(tsp); wr(rtd); wr(dth); wr(precipw);
      wr(rainnc); wr(rainncv); wr(snownc); wr(snowncv); wr(graupelnc); wr(graupelncv); wr(sr); wr(rainprod);
      wr(evapprod); wr(rec); wr(rei); wr(res);
      for (auto& t : ten) wr(t);
      if (diags) {
        wr(refl10cm); wr(refl_max); wr(refl_1km); wr(refl_1km_max); wr(relhum);
      }
    }
    std::fclose(out);
  } catch (const std::exception& e) {
    std::fprintf(stderr, "%s\n", e.what());
    return 1;
  }
  return 0;
}
