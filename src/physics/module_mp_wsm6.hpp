// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// MPAS's WSM6 set-up and wrapper: mp_wsm6_init as init_microphysics calls it
// (physics/mpas_atmphys_driver_microphysics.F:305-307) and subroutine wsm6 of physics/physics_wrf/module_mp_wsm6.F
// (MPAS v8.4.2), which copies each column into the arguments of mp_wsm6_run (physics/mp_wsm6.hpp, shared with
// CPPWRF) and calls mp_wsm6_effectRad_run. Arrays are the driver's (ims:ime, kms:kme, jms:jme) with jms = jme = 1,
// i fastest.
#include "core/precision.hpp"
#include "physics/mp_wsm6.hpp"

namespace cppmpas::phys::module_mp_wsm6 {

// mp_wsm6_init(den0=rho_a, denr=rho_r, dens=rho_s, cl=cliq, cpv=cpv, hail_opt=0) (DRV:91, 305-307), and the
// constants driver_microphysics passes to wsm6 (DRV:439-460: g = gravity, cpd = cp, rd = R_d, rv = R_v,
// t0c = svpt0, ep1 = ep_1, ep2 = ep_2, qmin = epsilon, xls, xlv0 = xlv, xlf0 = xlf, den0 = rho_a, denr = rho_w,
// cliq, cice, psat), which mp_wsm6_run reads as its arguments
phys::wsm6::Wsm6Params<real> mpas_wsm6_params();

struct Wsm6Columns {
  int ldi = 0;  // ime - ims + 1: the stride of k
  real *th = nullptr, *q = nullptr, *qc = nullptr, *qr = nullptr, *qi = nullptr, *qs = nullptr, *qg = nullptr;
  const real *den = nullptr, *pii = nullptr, *p = nullptr, *delz = nullptr;
  real *rain = nullptr, *rainncv = nullptr, *snow = nullptr, *snowncv = nullptr, *graupel = nullptr,
       *graupelncv = nullptr, *sr = nullptr;
};

// module_mp_wsm6.F:18-235 for its..ite (1-based, ims = 1) and kts = 1..kte. has_reqc / has_reqi / has_reqs = 1 (the
// effective radii of mp_wsm6_effectRad_run, config_microp_re) are refused: with 0 mp_wsm6_effectRad_run returns at
// once (mp_wsm6_effectRad.F90:113) and re_cloud / re_ice / re_snow are left as they are.
void wsm6(const phys::wsm6::Wsm6Params<real>& params, real delt, int has_reqc, int has_reqi, int has_reqs, int its, int ite,
          int kts, int kte, const Wsm6Columns& a);

}  // namespace cppmpas::phys::module_mp_wsm6
