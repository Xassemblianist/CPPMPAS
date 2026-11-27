// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// MPAS-Atmosphere's New Tiedtke wrapper (physics/physics_wrf/module_cu_ntiedtke.F, NT, with cu_ntiedtke_pre.F, PR, and
// cu_ntiedtke_post.F, PO), statement by statement: cu_ntiedtke_driver (NT:20-276) calls cu_ntiedtke_init (NT:161-165)
// and, for the row j = 1, cu_ntiedtke_pre_run (PR:63-183), cu_ntiedtke_run (MMM-physics cu_ntiedtke.F90, the shared
// column kernel of physics/cu_ntiedtke.hpp) and cu_ntiedtke_post_run (PO:63-116) on its..ite, then copies the
// tendencies back (NT:240-272).
//
// MPAS calls cu_ntiedtke_run with all columns its..ite at once (lq = im = ite - its + 1); the columns are independent
// except for the one vectorised libm loop of the kernel (cu_ntiedtke.hpp: cuinin's cuadjtqn, _ZGVbN2v_exp for the
// pairs of columns from its), so the kernel runs column by column with NtConst::cuinin_mvec set by CPPWRF's lane rule.
//
// Reference object code (bin/atmosphere_model, sha256 524a1cf4...): cu_ntiedtke_driver, ntiedtkeinit,
// cu_ntiedtke_pre_run and cu_ntiedtke_post_run are instruction-identical to MPAS's files compiled standalone with the
// reference flags (-Dmpas): no FMA, no ymm register, no libm call.
#include "core/precision.hpp"
#include "physics/mpas_atmphys_interface.hpp"

namespace cppmpas::phys {

// The arguments of cu_ntiedtke_driver (NT:20-33): (ims:ime, kms:kme) arrays of mpas_atmphys_vars (PArray, 1-based)
struct NtiedtkeArgs {
  real dt = 0;
  int itimestep = 0, stepcu = 0;  // itimestep: driver_convection's initflag
  PArray *raincv = nullptr, *pratec = nullptr;
  const PArray *qfx = nullptr, *hfx = nullptr, *xland = nullptr, *dx = nullptr;
  const PArray *u3d = nullptr, *v3d = nullptr, *w = nullptr, *t3d = nullptr, *qv3d = nullptr, *qc3d = nullptr,
               *qi3d = nullptr, *pi3d = nullptr, *rho3d = nullptr, *qvften = nullptr, *thften = nullptr,
               *dz8w = nullptr, *pcps = nullptr, *p8w = nullptr;
  bool f_qc = true, f_qi = true;
  real grav = 0, xlf = 0, xls = 0, xlv = 0, rd = 0, rv = 0, cp = 0;
  PArray *rthcuten = nullptr, *rqvcuten = nullptr, *rqccuten = nullptr, *rqicuten = nullptr, *rucuten = nullptr,
         *rvcuten = nullptr;
  int kts = 1, kte = 0;
};

// NT:20-276 on its..ite, j = 1
void cu_ntiedtke_driver(const NtiedtkeArgs& a, int its, int ite);

}  // namespace cppmpas::phys
