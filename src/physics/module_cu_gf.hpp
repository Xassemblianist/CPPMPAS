// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// MPAS-Atmosphere's cu_grell_freitas (physics/physics_wrf/module_cu_gf.mpas.F:54-453, GF), statement by statement, for
// the row j = 1 of cells its..ite: the initialisations (GF:165-312), the moisture convergence (GF:315-331), CUP_gf,
// neg_check, CUP_gf_sh (ishallow_g3 = 1) and the tendencies (GF:421-449). The columns are independent (cu_gf.hpp), so
// each runs as one kernel call (cu_gf_column) on its own arrays. use_excess = use_excess_sh = 0: the w* excess block
// (GF:351-374) is not reached; umean / vmean / pmean (GF:376-385) feed nothing.
//
// Refused: kts /= 1, more than cu_gf::kMaxLevels levels, ichoice_deep > 16 or ichoice_shallow > 9 (indices past
// xf_ens / xff_shal), ishallow_g3 /= 1.
#include <vector>

#include "core/precision.hpp"
#include "physics/cu_gf.hpp"
#include "physics/mpas_atmphys_interface.hpp"

namespace cppmpas::phys {

// An integer (ims:ime, jms:jme) array of mpas_atmphys_vars, 1-based, jms = jme = 1
struct IArray {
  std::vector<int> a;
  void allocate(int ni) { a.assign(static_cast<std::size_t>(ni), 0); }
  int& operator()(int i) { return a[static_cast<std::size_t>(i - 1)]; }
  const int& operator()(int i) const { return a[static_cast<std::size_t>(i - 1)]; }
};

// The arguments of cu_grell_freitas (GF:54-65), (ims:ime, kms:kme) arrays of mpas_atmphys_vars
struct GfArgs {
  int itimestep = 0;
  real dt = 0;
  const PArray *dxCell = nullptr, *areaCell = nullptr;
  const PArray *u = nullptr, *v = nullptr, *t = nullptr, *w = nullptr, *q = nullptr, *p = nullptr, *pi = nullptr,
               *rho = nullptr, *dz8w = nullptr, *p8w = nullptr;
  const PArray *xland = nullptr, *ht = nullptr, *hfx = nullptr, *qfx = nullptr, *gsw = nullptr;
  const PArray *rqvften = nullptr, *rthften = nullptr, *rthblten = nullptr, *rqvblten = nullptr, *rthraten = nullptr;
  const IArray* kpbl = nullptr;
  real xlv = 0, cp = 0, g = 0, r_v = 0;
  int ichoice_deep = 0, ichoice_shallow = 0, ishallow_g3 = 0;
  PArray *htop = nullptr, *hbot = nullptr;
  IArray *k22_shallow = nullptr, *kbcon_shallow = nullptr, *ktop_shallow = nullptr;
  PArray *xmb_total = nullptr, *xmb_shallow = nullptr, *raincv = nullptr, *pratec = nullptr, *gdc = nullptr,
         *gdc2 = nullptr;
  PArray *rthcuten = nullptr, *rqvcuten = nullptr, *rqccuten = nullptr, *rqicuten = nullptr;
  int kts = 1, kte = 0;
};

// GF:54-453 on its..ite, j = 1
void cu_grell_freitas(const GfArgs& a, int its, int ite);

}  // namespace cppmpas::phys
