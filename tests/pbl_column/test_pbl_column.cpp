// SPDX-License-Identifier: GPL-3.0-or-later
// The C++ side of the surface layer / YSU PBL column check: test_pbl_column <in.bin> <out.bin> runs MPAS_to_physics,
// driver_sfclayer, driver_pbl and physics_get_tend (src/physics) on the columns of in.bin in the order harness.F90 runs
// MPAS's own code, and writes the same records (tests/pbl_column/columns.py).
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <map>
#include <string>
#include <vector>

#include "core/error.hpp"
#include "core/field.hpp"
#include "physics/mpas_atmphys_interface.hpp"
#include "physics/mpas_atmphys_pbl.hpp"
#include "physics/mpas_atmphys_sfclayer.hpp"
#include "physics/mpas_atmphys_todynamics.hpp"

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
FView<int> RI(const std::string& n) {
  Rec& r = rec(n);
  return FView<int>{r.iv.data(), r.dims[0], r.nd > 1 ? r.dims[1] : 1, 1};
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

int main(int argc, char** argv) {
  if (argc != 3) {
    std::fprintf(stderr, "usage: test_pbl_column <in.bin> <out.bin>\n");
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
    const int nC = hdr[0], nCt = hdr[1], nEs = hdr[2], nE = hdr[3], nV = hdr[4], nsteps = hdr[5];
    const bool frac_seaice = hdr[6] == 1, ysu_pblmix = hdr[7] == 1;
    const char* cu_names[] = {"off", "cu_ntiedtke", "cu_kain_fritsch", "cu_grell_freitas"};
    const std::string convection = cu_names[hdr[8]];
    const std::string radt_lw = hdr[9] == 1 ? "rrtmg_lw" : "off", radt_sw = hdr[9] == 1 ? "rrtmg_sw" : "off";
    const real dt = rec("dtl").r[0], len_disp = rec("dtl").r[1];
    const std::string microp = "mp_wsm6", pbl = "bl_ysu", sfclayer = "sf_monin_obukhov_rev";

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
    sp.config_len_disp = len_disp;
    sp.meshDensity = C("meshDensity");
    sp.areaCell = C("areaCell");
    sp.skintemp = C("skintemp");
    sp.xland = C("xland");
    sp.sst = C("sst");
    sp.xice = C("xice");
    sp.hpbl = C("hpbl");
    sp.mavail = C("mavail");
    sp.br = R("br"); sp.cpm = R("cpm"); sp.chs = R("chs"); sp.chs2 = R("chs2"); sp.cqs2 = R("cqs2");
    sp.flhc = R("flhc"); sp.flqc = R("flqc"); sp.gz1oz0 = R("gz1oz0"); sp.hfx = R("hfx"); sp.qfx = R("qfx");
    sp.qgh = R("qgh"); sp.qsfc = R("qsfc"); sp.lh = R("lh"); sp.mol = R("mol"); sp.psim = R("psim");
    sp.psih = R("psih"); sp.regime = R("regime"); sp.rmol = R("rmol"); sp.ust = R("ust"); sp.ustm = R("ustm");
    sp.wspd = R("wspd"); sp.znt = R("znt"); sp.zol = R("zol"); sp.q2 = R("q2"); sp.t2m = R("t2m");
    sp.th2m = R("th2m"); sp.u10 = R("u10"); sp.v10 = R("v10"); sp.cd = R("cd"); sp.cda = R("cda");
    sp.ck = R("ck"); sp.cka = R("cka"); sp.fh = R("fh"); sp.fm = R("fm");
    sp.xice_threshold = frac_seaice ? 0.02 : 0.5;  // mpas_atmphys_manager.F:740-744

    ph::PblPools bp;
    bp.nVertLevels = nV;
    bp.pbl_scheme = pbl;
    bp.config_ysu_pblmix = ysu_pblmix;
    bp.xland = C("xland");
    bp.hfx = C("hfx"); bp.qfx = C("qfx"); bp.ust = C("ust"); bp.wspd = C("wspd"); bp.znt = C("znt");
    bp.br = C("br"); bp.fm = C("fm"); bp.fh = C("fh"); bp.u10 = C("u10"); bp.v10 = C("v10");
    bp.hpbl = R("hpbl"); bp.delta = R("delta"); bp.wstar = R("wstar"); bp.kzh = R("kzh"); bp.kzm = R("kzm");
    bp.kzq = R("kzq"); bp.exch_h = R("exch_h");
    bp.kpbl = RI("kpbl");
    bp.rthratenlw = C("rthratenlw"); bp.rthratensw = C("rthratensw");
    bp.rublten = R("rublten"); bp.rvblten = R("rvblten"); bp.rthblten = R("rthblten");
    bp.rqvblten = R("rqvblten"); bp.rqcblten = R("rqcblten"); bp.rqiblten = R("rqiblten");
    bp.dt_pbl = dt;
    bp.f_qc = true;
    bp.f_qi = true;

    ph::GetTendPools gp;
    gp.nCells = nCt; gp.nEdges = nE; gp.nCellsSolve = nC; gp.nEdgesSolve = nEs; gp.nVertLevels = nV;
    gp.pbl_scheme = pbl; gp.convection_scheme = convection; gp.microp_scheme = microp;
    gp.radt_lw_scheme = radt_lw; gp.radt_sw_scheme = radt_sw;
    {
      FView<int> coe = RI("cellsOnEdge");
      gp.cellsOnEdge = FView<const int>{coe.p, coe.e0, coe.e1, coe.e2};
    }
    gp.edgeNormalVectors = C("edgeNormalVectors"); gp.east = C("east"); gp.north = C("north");
    gp.theta_m = C("theta_m"); gp.scalars = C("scalars"); gp.mass = C("rho_zz"); gp.mass_edge = C("rho_edge");
    gp.index_qv = 1; gp.index_qc = 2; gp.index_qr = 3; gp.index_qi = 4; gp.index_qs = 5;
    gp.tend_u_phys = R("tend_u_phys");
    gp.rublten = C("rublten"); gp.rvblten = C("rvblten"); gp.rthblten = C("rthblten"); gp.rqvblten = C("rqvblten");
    gp.rqcblten = C("rqcblten"); gp.rqiblten = C("rqiblten");
    gp.rucuten = C("rucuten"); gp.rvcuten = C("rvcuten"); gp.rthcuten = C("rthcuten"); gp.rqvcuten = C("rqvcuten");
    gp.rqccuten = C("rqccuten"); gp.rqrcuten = C("rqrcuten"); gp.rqicuten = C("rqicuten"); gp.rqscuten = C("rqscuten");
    gp.rthratenlw = C("rthratenlw"); gp.rthratensw = C("rthratensw");
    gp.rublten_Edge = R("rublten_Edge"); gp.rucuten_Edge = R("rucuten_Edge");
    gp.tend_scalars = R("scalars_tend");

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
      dumpv("@qv_p", v.qv_p, 1, nV); dumpv("@qc_p", v.qc_p, 1, nV); dumpv("@qr_p", v.qr_p, 1, nV);
      dumpv("@qi_p", v.qi_p, 1, nV); dumpv("@qs_p", v.qs_p, 1, nV); dumpv("@qg_p", v.qg_p, 1, nV);
      dumpv("@u_p", v.u_p, 1, nV); dumpv("@v_p", v.v_p, 1, nV); dumpv("@zz_p", v.zz_p, 1, nV);
      dumpv("@rho_p", v.rho_p, 1, nV); dumpv("@th_p", v.th_p, 1, nV); dumpv("@t_p", v.t_p, 1, nV);
      dumpv("@pi_p", v.pi_p, 1, nV); dumpv("@pres_p", v.pres_p, 1, nV); dumpv("@zmid_p", v.zmid_p, 1, nV);
      dumpv("@dz_p", v.dz_p, 1, nV); dumpv("@znu_p", v.znu_p, 1, nV);
      dumpv("@w_p", v.w_p, 1, nV + 1); dumpv("@z_p", v.z_p, 1, nV + 1);
      dumpv("@fzm_p", v.fzm_p, 2, nV); dumpv("@fzp_p", v.fzp_p, 2, nV);
      dumpv("@t2_p", v.t2_p, 1, nV + 1); dumpv("@pres2_p", v.pres2_p, 1, nV + 1);
      dumpv("@pres2_hyd_p", v.pres2_hyd_p, 1, nV + 1); dumpv("@pres2_hydd_p", v.pres2_hydd_p, 1, nV + 1);
      dumpv("@pres_hyd_p", v.pres_hyd_p, 1, nV); dumpv("@pres_hydd_p", v.pres_hydd_p, 1, nV);
      dumpv("@znu_hyd_p", v.znu_hyd_p, 1, nV);
      put_r(out, "@psfc_p", 1, {nC}, v.psfc_p.data(), nC);
      put_r(out, "@psfc_hyd_p", 1, {nC}, v.psfc_hyd_p.data(), nC);
      put_r(out, "@psfc_hydd_p", 1, {nC}, v.psfc_hydd_p.data(), nC);
      ph::driver_sfclayer(step, tables, sp, v, 1, nC);
      marker(out, "#B", step);
      dump_all(out);
      ph::driver_pbl(step, bp, v, 1, nC);
      marker(out, "#C", step);
      dump_all(out);
      for (int rk = 1; rk <= 2; ++rk) {
        ph::physics_get_tend(gp, rk, 1, R("tend_ru_physics"), R("tend_rtheta_physics"), R("tend_rho_physics"),
                             [](const char*) {});  // one rank: the halo cells hold what in.bin gives them
        marker(out, rk == 1 ? "#D" : "#E", step);
        dump_all(out);
      }
    }
    std::fclose(out);
  } catch (const std::exception& e) {
    std::fprintf(stderr, "%s\n", e.what());
    return 1;
  }
  return 0;
}
