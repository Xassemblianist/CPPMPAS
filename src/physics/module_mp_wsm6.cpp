// SPDX-License-Identifier: GPL-3.0-or-later
// See module_mp_wsm6.hpp. MW = physics/physics_wrf/module_mp_wsm6.F, DRV = physics/mpas_atmphys_driver_microphysics.F.
#include "physics/module_mp_wsm6.hpp"

#include <vector>

#include "core/error.hpp"
#include "physics/mpas_atmphys_constants.hpp"

namespace cppmpas::phys::module_mp_wsm6 {

phys::wsm6::Wsm6Params<real> mpas_wsm6_params() {
  namespace c = physics_constants;
  phys::wsm6::Wsm6Params<real> w;
  w.init(c::rho_a, c::rho_r, c::rho_s, c::cliq, c::cpv, 0);  // DRV:305-307, hail_opt = 0 (DRV:91)
  // DRV:439-460: the arguments of mp_wsm6_run (MW:148-154)
  w.g = c::gravity;
  w.cpd = c::cp;
  w.cpv = c::cpv;
  w.rd = c::R_d;
  w.rv = c::R_v;
  w.t0c = c::svpt0;
  w.ep1 = c::ep_1;
  w.ep2 = c::ep_2;
  w.qmin = c::epsilon;
  w.xls = c::xls;
  w.xlv0 = c::xlv;
  w.xlf0 = c::xlf;
  w.den0 = c::rho_a;
  w.denr = c::rho_w;
  w.cliq = c::cliq;
  w.cice = c::cice;
  w.psat = c::psat;
  return w;
}

void wsm6(const phys::wsm6::Wsm6Params<real>& params, real delt, int has_reqc, int has_reqi, int has_reqs, int its, int ite,
          int kts, int kte, const Wsm6Columns& a) {
  if (kts != 1) throw Error("wsm6: kts /= 1");
  if (kte > phys::wsm6::kMaxLevels) throw Error("wsm6: more than 128 levels");
  if (has_reqc == 1 && has_reqi == 1 && has_reqs == 1)  // MW:209-210: do_microp_re
    throw Error("wsm6: the effective radii of mp_wsm6_effectRad_run (config_microp_re) have no verified port");
  const int nz = kte;
  std::vector<real> th(nz), q(nz), qc(nz), qr(nz), qi(nz), qs(nz), qg(nz), den(nz), pii(nz), p(nz), delz(nz);
  auto at = [&](int i, int k) { return static_cast<std::size_t>(k - 1) * a.ldi + (i - 1); };
  for (int i = its; i <= ite; ++i) {
    // MW:114-145: the inputs of column i (t_hv = th*pii is formed inside the kernel, operator() of mp_wsm6.hpp)
    for (int k = 1; k <= nz; ++k) {
      th[k - 1] = a.th[at(i, k)];
      q[k - 1] = a.q[at(i, k)];
      qc[k - 1] = a.qc[at(i, k)];
      qr[k - 1] = a.qr[at(i, k)];
      qi[k - 1] = a.qi[at(i, k)];
      qs[k - 1] = a.qs[at(i, k)];
      qg[k - 1] = a.qg[at(i, k)];
      den[k - 1] = a.den[at(i, k)];
      pii[k - 1] = a.pii[at(i, k)];
      p[k - 1] = a.p[at(i, k)];
      delz[k - 1] = a.delz[at(i, k)];
    }
    real rain = a.rain[i - 1], snow = a.snow[i - 1], graupel = a.graupel[i - 1];
    real rainncv = 0, snowncv = 0, graupelncv = 0, sr = 0;
    // MW:148-158: mp_wsm6_run over its..ite; column i is iteration i - its of the vectorised i loops
    params(nz, delt, th.data(), q.data(), qc.data(), qr.data(), qi.data(), qs.data(), qg.data(), den.data(),
           pii.data(), p.data(), delz.data(), rain, rainncv, snow, snowncv, graupel, graupelncv, sr, i - its,
           ite - its + 1);
    // MW:160-188 (th = t_hv/pii in the kernel)
    a.rain[i - 1] = rain;
    a.rainncv[i - 1] = rainncv;
    a.sr[i - 1] = sr;
    for (int k = 1; k <= nz; ++k) {
      a.th[at(i, k)] = th[k - 1];
      a.q[at(i, k)] = q[k - 1];
      a.qc[at(i, k)] = qc[k - 1];
      a.qi[at(i, k)] = qi[k - 1];
      a.qr[at(i, k)] = qr[k - 1];
      a.qs[at(i, k)] = qs[k - 1];
      a.qg[at(i, k)] = qg[k - 1];
    }
    a.snow[i - 1] = snow;
    a.snowncv[i - 1] = snowncv;
    a.graupel[i - 1] = graupel;
    a.graupelncv[i - 1] = graupelncv;
  }
  // MW:208-231: mp_wsm6_effectRad_run with do_microp_re = .false. returns at once (refused above otherwise)
}

}  // namespace cppmpas::phys::module_mp_wsm6
