// SPDX-License-Identifier: GPL-3.0-or-later
// See module_cu_gf.hpp. GF = src/core_atmosphere/physics/physics_wrf/module_cu_gf.mpas.F (MPAS v8.4.2).
#include "physics/module_cu_gf.hpp"

#include <memory>

#include "core/error.hpp"

namespace cppmpas::phys {

void cu_grell_freitas(const GfArgs& a, int its, int ite) {
  namespace gf = cu_gf;
  const int kts = a.kts, kte = a.kte;
  if (kts != 1) throw Error("cu_grell_freitas: kts /= 1");
  if (kte > gf::kMaxLevels) throw Error("cu_grell_freitas: more than 100 levels");
  if (its != 1) throw Error("cu_grell_freitas: its /= 1");
  if (a.ishallow_g3 != 1) throw Error("cu_grell_freitas: ishallow_g3 /= 1 has no verified port");
  if (a.ichoice_deep > 16)
    throw Error("cu_grell_freitas: config_gfconv_closure_deep > 16 (MPAS indexes past xf_ens)");
  if (a.ichoice_shallow > 9)
    throw Error("cu_grell_freitas: config_gfconv_closure_shallow > 9 (MPAS indexes past xff_shal)");

  auto in = std::make_unique<gf::GfDriverIn>();
  auto out = std::make_unique<gf::GfDriverOut>();
  for (int i = its; i <= ite; ++i) {
    gf::GfDriverIn& c = *in;
    c.kte = kte;
    c.itimestep = a.itimestep;
    c.dt = a.dt;
    c.dxCell = (*a.dxCell)(i);
    c.areaCell = (*a.areaCell)(i);
    c.ht = (*a.ht)(i);
    c.xland = (*a.xland)(i);
    c.gsw = (*a.gsw)(i);
    c.hfx = (*a.hfx)(i);
    c.qfx = (*a.qfx)(i);
    c.p8w1 = (*a.p8w)(i, 1);
    c.kpbl = (*a.kpbl)(i);
    for (int k = kts; k <= kte; ++k) {
      c.u[k] = (*a.u)(i, k);
      c.v[k] = (*a.v)(i, k);
      c.t[k] = (*a.t)(i, k);
      c.q[k] = (*a.q)(i, k);
      c.p[k] = (*a.p)(i, k);
      c.pi[k] = (*a.pi)(i, k);
      c.rho[k] = (*a.rho)(i, k);
      c.dz8w[k] = (*a.dz8w)(i, k);
      c.w[k] = (*a.w)(i, k);
      c.rqvften[k] = (*a.rqvften)(i, k);
      c.rthften[k] = (*a.rthften)(i, k);
      c.rthblten[k] = (*a.rthblten)(i, k);
      c.rqvblten[k] = (*a.rqvblten)(i, k);
      c.rthraten[k] = (*a.rthraten)(i, k);
    }
    c.xlv = a.xlv;
    c.cp = a.cp;
    c.g = a.g;
    c.r_v = a.r_v;
    c.ichoice_deep = a.ichoice_deep;
    c.ichoice_shallow = a.ichoice_shallow;
    c.ishallow_g3 = a.ishallow_g3;
    if (c.kpbl > kte) throw Error("cu_grell_freitas: kpbl > kte");

    gf::cu_grell_freitas_column(c, *out);

    const gf::GfDriverOut& o = *out;
    (*a.hbot)(i) = o.hbot;
    (*a.htop)(i) = o.htop;
    (*a.xmb_total)(i) = o.xmb_total;
    (*a.raincv)(i) = o.raincv;
    (*a.pratec)(i) = o.pratec;
    (*a.k22_shallow)(i) = o.k22_shallow;
    (*a.kbcon_shallow)(i) = o.kbcon_shallow;
    (*a.ktop_shallow)(i) = o.ktop_shallow;
    (*a.xmb_shallow)(i) = o.xmb_shallow;
    for (int k = kts; k <= kte; ++k) {
      (*a.gdc)(i, k) = o.gdc[k];
      (*a.gdc2)(i, k) = o.gdc2[k];
      (*a.rthcuten)(i, k) = o.rthcuten[k];
      (*a.rqvcuten)(i, k) = o.rqvcuten[k];
      (*a.rqccuten)(i, k) = o.rqccuten[k];
      (*a.rqicuten)(i, k) = o.rqicuten[k];
    }
  }
}

}  // namespace cppmpas::phys
