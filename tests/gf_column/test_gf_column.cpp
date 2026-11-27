// SPDX-License-Identifier: GPL-3.0-or-later
// The C++ side of the Grell-Freitas column check: test_gf_column <in.bin> <out.bin> runs MPAS_to_physics,
// update_convection_step1, driver_convection (cu_grell_freitas), update_convection_step2 and physics_get_tend
// (src/physics) on the columns of in.bin in the order harness.F90 runs MPAS's own code, and writes the same records
// (tests/gf_column/columns.py). The record reader and writer are tests/cu_gwdo_column/test_cu_gwdo_column.cpp's.
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <map>
#include <string>
#include <vector>

#include "core/error.hpp"
#include "core/field.hpp"
#include "physics/mpas_atmphys_convection.hpp"
#include "physics/mpas_atmphys_interface.hpp"
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
    std::fprintf(stderr, "usage: test_gf_column <in.bin> <out.bin>\n");
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
    const bool acrain = hdr[6] == 1, do_restart = hdr[7] == 1;
    const int closure_deep = hdr[8], closure_shallow = hdr[9];
    const real dt = rec("dtl").r[0], len_disp = rec("dtl").r[1], bucket_rainc = rec("dtl").r[2];
    const std::string microp = "mp_wsm6", pbl = "bl_ysu";

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

    ph::ConvectionPools cv;
    cv.nVertLevels = nV;
    cv.convection_scheme = "cu_grell_freitas";
    cv.config_do_restart = do_restart;
    cv.config_len_disp = len_disp;
    cv.config_bucket_rainc = bucket_rainc;
    cv.meshDensity = C("meshDensity");
    cv.xland = C("xland");
    cv.hfx = C("hfx"); cv.qfx = C("qfx");
    cv.cuprec = R("cuprec"); cv.raincv = R("raincv"); cv.rainc = R("rainc");
    cv.i_rainc = RI("i_rainc");
    cv.rqvdynten = C("rqvdynten"); cv.rqvblten = C("rqvblten"); cv.rthdynten = C("rthdynten");
    cv.rthblten = C("rthblten"); cv.rthratenlw = C("rthratenlw"); cv.rthratensw = C("rthratensw");
    cv.rthcuten = R("rthcuten"); cv.rqvcuten = R("rqvcuten"); cv.rqccuten = R("rqccuten");
    cv.rqicuten = R("rqicuten"); cv.rucuten = R("rucuten"); cv.rvcuten = R("rvcuten");
    cv.config_gfconv_closure_deep = closure_deep;
    cv.config_gfconv_closure_shallow = closure_shallow;
    cv.areaCell = C("areaCell");
    cv.zgrid = C("zgrid");
    cv.gsw = C("gsw");
    {
      FView<int> kp = RI("kpbl");
      cv.kpbl = FView<const int>{kp.p, kp.e0, kp.e1, kp.e2};
    }
    cv.cubot = R("cubot"); cv.cutop = R("cutop"); cv.xmb_total = R("xmb_total"); cv.xmb_shallow = R("xmb_shallow");
    cv.qc_cu = R("qc_cu"); cv.qi_cu = R("qi_cu");
    cv.k22_shallow = RI("k22_shallow"); cv.kbcon_shallow = RI("kbcon_shallow"); cv.ktop_shallow = RI("ktop_shallow");
    cv.ktop_deep = RI("ktop_deep");
    cv.dt_dyn = dt;
    cv.n_cu = 1;
    cv.f_qc = true;
    cv.f_qi = true;
    cv.l_acrain = acrain;

    ph::GetTendPools gp;
    gp.nCells = nCt; gp.nEdges = nE; gp.nCellsSolve = nC; gp.nEdgesSolve = nEs; gp.nVertLevels = nV;
    gp.pbl_scheme = pbl; gp.convection_scheme = "cu_grell_freitas"; gp.microp_scheme = microp;
    gp.radt_lw_scheme = "rrtmg_lw"; gp.radt_sw_scheme = "rrtmg_sw";
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

    std::FILE* out = std::fopen(argv[2], "wb");
    if (!out) return 2;
    for (int step = 1; step <= nsteps; ++step) {
      ph::PhysicsVars v;
      v.allocate_forall_physics(microp, pbl, nC, nV + 1);
      ph::MPAS_to_physics(tp, v, 1, nC);
      marker(out, "#A", step);
      dump_all(out);
      ph::update_convection_step1(cv, 1, nC);
      ph::driver_convection(step, cv, v, 1, nC);  // l_conv: config_conv_interval = 'none'
      ph::update_convection_step2(cv, 1, nC);
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
