// SPDX-License-Identifier: GPL-3.0-or-later
// See module_cu_ntiedtke.hpp. NT = src/core_atmosphere/physics/physics_wrf/module_cu_ntiedtke.F, PR = physics_wrf/
// cu_ntiedtke_pre.F, PO = physics_wrf/cu_ntiedtke_post.F (MPAS v8.4.2).
#include "physics/module_cu_ntiedtke.hpp"

#include <cmath>

#include "core/backend.hpp"

#include "core/error.hpp"
#include "physics/cu_ntiedtke.hpp"

namespace cppmpas::phys {

namespace {
// CPPWRF's lane rule (CPPWRF 7fac5bf src/core/backend.hpp:197): column i - its of n in a libmvec pair (the pairs from
// its; the odd last column and a 1-column loop are scalar)
CPPMPAS_HD inline bool lane(int i, int n) { return n > 1 && !((n & 1) && i == n - 1); }
}  // namespace

void cu_ntiedtke_driver(const NtiedtkeArgs& a, int its, int ite) {
  namespace nt = cu_ntiedtke;
  using nt::KM;
  const int kts = a.kts, kte = a.kte;
  if (kts != 1) throw Error("cu_ntiedtke_driver: kts /= 1");
  if (kte > nt::kMaxLevels) throw Error("cu_ntiedtke_driver: more than 100 levels");
  if (its != 1) throw Error("cu_ntiedtke_driver: its /= 1");

  // NT:161-165: cu_ntiedtke_init (cu_ntiedtke_common)
  const nt::NtConst<real> c0 = nt::make_ntconst<real>(a.cp, a.rd, a.rv, a.xlv, a.xls, a.xlf, a.grav);
  // NT:167-171: cu_act_flag = .true. (driver_convection's local, not used after the call for cu_ntiedtke)

  // NT:173-274, j = jts = jte = 1. PR:101-103
  const int im = ite - its + 1, kx = kte - kts + 1;
  const real delt = a.dt * real(a.stepcu);  // PR:105
  // one thread per column i (the columns of cu_ntiedtke_run are independent)
  parallel_for(its, ite, [=] CPPMPAS_LAMBDA(int i) {
    nt::NtConst<real> c = c0;  // the column's copy: cuinin_mvec below is the column's lane
    // the _hv copies (NT:175-204) are the (i, k) values; cu_ntiedtke_pre_run (PR:63-183) for column i
    const int slimsk = int(std::fabs(a.xland(i) - 2.));  // PR:107-109 (integer = real: truncation)
    real zi[KM], zl[KM], dot[KM];
    real ghti[KM], prsi[KM], ghtl[KM], omg[KM], prsl[KM];
    real tf[KM], qvf[KM], qcf[KM], qif[KM], uf[KM], vf[KM], qvftenz[KM], thftenz[KM];
    zi[kts] = 0.;  // PR:111-114
    for (int k = kts; k <= kte; ++k) zi[k + 1] = zi[k] + a.dz8w(i, k);  // PR:115-119
    for (int k = kts; k <= kte; ++k) {                                       // PR:120-125
      zl[k] = 0.5 * (zi[k] + zi[k + 1]);
      dot[k] = -0.5 * a.grav * a.rho3d(i, k) * (a.w(i, k) + a.w(i, k + 1));
    }
    {
      int pp = 0;  // PR:127-135
      for (int k = kts; k <= kte + 1; ++k) {
        const int zz = kte + 1 - pp;
        ghti[zz] = zi[k];
        prsi[zz] = a.p8w(i, k);
        pp = pp + 1;
      }
      pp = 0;  // PR:136-145
      for (int k = kts; k <= kte; ++k) {
        const int zz = kte - pp;
        ghtl[zz] = zl[k];
        omg[zz] = dot[k];
        prsl[zz] = a.pcps(i, k);
        pp = pp + 1;
      }
      pp = 0;  // PR:147-159
      for (int k = kts; k <= kte; ++k) {
        const int zz = kte - pp;
        tf[zz] = a.t3d(i, k);
        qvf[zz] = a.qv3d(i, k);
        qcf[zz] = a.qc3d(i, k);
        qif[zz] = a.qi3d(i, k);
        uf[zz] = a.u3d(i, k);
        vf[zz] = a.v3d(i, k);
        pp = pp + 1;
      }
      if (a.itimestep == 1) {  // PR:161-178
        for (int k = kts; k <= kte; ++k) {
          qvftenz[k] = 0.;
          thftenz[k] = 0.;
        }
      } else {
        pp = 0;
        for (int k = kts; k <= kte; ++k) {
          const int zz = kte - pp;
          qvftenz[zz] = a.qvften(i, k);
          thftenz[zz] = a.thften(i, k);
          pp = pp + 1;
        }
      }
    }

    // NT:220-227: cu_ntiedtke_run (lq = im columns; column i here). The -O3 lanes of cuinin's cuadjtqn: column
    // i - its of im (cu_ntiedtke.hpp)
    c.cuinin_mvec = lane(i - its, im);
    real rn = 0.;
    nt::cu_ntiedtke_run<real>(c, kx, uf, vf, tf, qvf, qcf, qif, qvftenz, thftenz, ghtl, ghti, omg, prsl, prsi,
                              a.qfx(i), a.hfx(i), rn, slimsk, delt, a.dx(i), nt::ColumnSlots{});

    // NT:229-238: cu_ntiedtke_post_run (PO:63-116) into the _hv arrays, then NT:240-272
    const real rdelt = 1. / delt;           // PO:91-92
    a.raincv(i) = rn / real(a.stepcu);   // PO:94-97, NT:240-243
    a.pratec(i) = rn / (real(a.stepcu) * a.dt);
    int pp = 0;
    for (int k = kts; k <= kte; ++k) {  // PO:99-111, NT:245-272
      const int zz = kte - pp;
      a.rthcuten(i, k) = (tf[zz] - a.t3d(i, k)) / a.pi3d(i, k) * rdelt;
      a.rqvcuten(i, k) = (qvf[zz] - a.qv3d(i, k)) * rdelt;
      if (a.f_qc) a.rqccuten(i, k) = (qcf[zz] - a.qc3d(i, k)) * rdelt;  // NT:254-262
      if (a.f_qi) a.rqicuten(i, k) = (qif[zz] - a.qi3d(i, k)) * rdelt;  // NT:264-272
      a.rucuten(i, k) = (uf[zz] - a.u3d(i, k)) * rdelt;
      a.rvcuten(i, k) = (vf[zz] - a.v3d(i, k)) * rdelt;
      pp = pp + 1;
    }
  });
}

}  // namespace cppmpas::phys
