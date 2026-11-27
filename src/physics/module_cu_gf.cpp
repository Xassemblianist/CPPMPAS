// SPDX-License-Identifier: GPL-3.0-or-later
// See module_cu_gf.hpp. GF = src/core_atmosphere/physics/physics_wrf/module_cu_gf.mpas.F (MPAS v8.4.2).
#include "physics/module_cu_gf.hpp"

#include "core/backend.hpp"
#include "core/error.hpp"

namespace cppmpas::phys {

void cu_grell_freitas(const GfArgs& a, int its, int ite) {
  namespace gf = cu_gf;
  const int kts = a.kts, kte = a.kte;
  if (kts != 1) throw Error("cu_grell_freitas: kts /= 1");
  if (kte > gf::kMaxLevels) throw Error("cu_grell_freitas: more than 100 levels");
  // cup_MAXIMI with ks = 3 (GF:941, 4461) reads level 3 and CUP_gf_sh uses ktf - 4: fewer levels index past MPAS's
  // arrays
  if (kte < 3) throw Error("cu_grell_freitas: fewer than 3 levels");
  if (its != 1) throw Error("cu_grell_freitas: its /= 1");
  if (a.ishallow_g3 != 1) throw Error("cu_grell_freitas: ishallow_g3 /= 1 has no verified port");
  if (a.ichoice_deep > 16)
    throw Error("cu_grell_freitas: config_gfconv_closure_deep > 16 (MPAS indexes past xf_ens)");
  if (a.ichoice_shallow > 9)
    throw Error("cu_grell_freitas: config_gfconv_closure_shallow > 9 (MPAS indexes past xff_shal)");

  // one thread per column i: its GfDriverIn / GfDriverOut and CUP_gf's work arrays are the thread's own
  const GfArgs A = a;
  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {
    gf::GfDriverIn c;
    gf::GfDriverOut o;
    c.kte = kte;
    c.itimestep = A.itimestep;
    c.dt = A.dt;
    c.dxCell = A.dxCell(i);
    c.areaCell = A.areaCell(i);
    c.ht = A.ht(i);
    c.xland = A.xland(i);
    c.gsw = A.gsw(i);
    c.hfx = A.hfx(i);
    c.qfx = A.qfx(i);
    c.p8w1 = A.p8w(i, 1);
    c.kpbl = A.kpbl(i);
    for (int k = kts; k <= kte; ++k) {
      c.u[k] = A.u(i, k);
      c.v[k] = A.v(i, k);
      c.t[k] = A.t(i, k);
      c.q[k] = A.q(i, k);
      c.p[k] = A.p(i, k);
      c.pi[k] = A.pi(i, k);
      c.rho[k] = A.rho(i, k);
      c.dz8w[k] = A.dz8w(i, k);
      c.w[k] = A.w(i, k);
      c.rqvften[k] = A.rqvften(i, k);
      c.rthften[k] = A.rthften(i, k);
      c.rthblten[k] = A.rthblten(i, k);
      c.rqvblten[k] = A.rqvblten(i, k);
      c.rthraten[k] = A.rthraten(i, k);
    }
    c.xlv = A.xlv;
    c.cp = A.cp;
    c.g = A.g;
    c.r_v = A.r_v;
    c.ichoice_deep = A.ichoice_deep;
    c.ichoice_shallow = A.ichoice_shallow;
    c.ishallow_g3 = A.ishallow_g3;
    if (c.kpbl > kte) kernel_fail("cu_grell_freitas: kpbl > kte");

    gf::cu_grell_freitas_column(c, o);

    A.hbot(i) = o.hbot;
    A.htop(i) = o.htop;
    A.xmb_total(i) = o.xmb_total;
    A.raincv(i) = o.raincv;
    A.pratec(i) = o.pratec;
    A.k22_shallow(i) = o.k22_shallow;
    A.kbcon_shallow(i) = o.kbcon_shallow;
    A.ktop_shallow(i) = o.ktop_shallow;
    A.xmb_shallow(i) = o.xmb_shallow;
    for (int k = kts; k <= kte; ++k) {
      A.gdc(i, k) = o.gdc[k];
      A.gdc2(i, k) = o.gdc2[k];
      A.rthcuten(i, k) = o.rthcuten[k];
      A.rqvcuten(i, k) = o.rqvcuten[k];
      A.rqccuten(i, k) = o.rqccuten[k];
      A.rqicuten(i, k) = o.rqicuten[k];
    }
  });
}

}  // namespace cppmpas::phys
