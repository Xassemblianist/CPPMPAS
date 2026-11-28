// SPDX-License-Identifier: GPL-3.0-or-later
// The Kessler path of MPAS-Atmosphere's cloud microphysics; see atm_microphysics.hpp. PI =
// src/core_atmosphere/physics/mpas_atmphys_interface.F, DRV = physics/mpas_atmphys_driver_microphysics.F,
// INIT = physics/mpas_atmphys_init.F, MGR = physics/mpas_atmphys_manager.F. Statement by statement.
//
// Reference object code (bin/atmosphere_model of the mpas-build asset 153c7206...): microphysics_from_mpas
// (0x1e8d80) and microphysics_to_mpas (0x1e39d0) have no FMA, no ymm register and no libmvec call; the only libm
// call is one pow@plt in microphysics_to_mpas (the exner of PI:926, in the i loop innermost of PI:903-905: scalar
// pow for every cell). The max(0._RKIND, qv) of PI:723 is a maxsd with qv as the destination ("qv > 0 ? qv : 0").
//
// The loops of driver_microphysics over the cells (i) are kernels (core/backend.hpp), one thread per column with the
// k loops inside, and kessler is CPPWRF's column kernel launched one thread per column as CPPWRF launches it
// (CPPWRF 7fac5bf src/phys/mp_kessler/kessler.cpp: column_for over the (i, j) columns); its libmvec lanes depend
// only on the column's index (vm::lane(i - i0, ni)), so every column takes the reference's lane in either build.
// The (ims:ime, kms:kme) work arrays are FArrays, in the managed memory of the CUDA build.
#include "core_atmosphere/atm_microphysics.hpp"

#include <cstdlib>

#include "core/error.hpp"
#include "core/libm.hpp"
#include "core/namelist.hpp"
#include "core/views.hpp"
#include "physics/mp_kessler.hpp"

namespace cppmpas {

namespace {

// src/framework/mpas_constants.F:43-56 and physics/mpas_atmphys_constants.F (R_d => rgas, folded as gfortran folds
// the parameter expressions, one correctly rounded operation per node)
constexpr real rgas = 287.0;
constexpr real rv = 461.6;
constexpr real cp = 7.0 * rgas / 2.0;
constexpr real rvord = rv / rgas;
constexpr real gravity = 9.80616;
constexpr real R_d = rgas;
constexpr real R_v = 461.6;          // mpas_atmphys_constants.F:42
constexpr real P0 = 100000.;         // :40
constexpr real ep_2 = R_d / R_v;     // :44
constexpr real rcv = R_d / (cp - R_d);  // :49
constexpr real rho_w = 1000.0;       // :54
constexpr real svp1 = 0.6112, svp2 = 17.67, svp3 = 29.65, svpt0 = 273.15;  // :55-58
constexpr real xlv = 2.50e6;         // :59

CPPMPAS_HD inline real fmax0(real x) { return x > 0.0 ? x : 0.0; }  // max(0._RKIND, x): maxsd, x the destination

std::string option(const std::string& namelist_path, const char* key) {
  const auto kv = effective_namelist(namelist_path);
  const auto it = kv.find(key);
  if (it == kv.end()) throw Error(std::string("microphysics: option ") + key + " is not in the Registry table");
  return it->second;
}

}  // namespace

void check_microphysics_supported(const std::string& namelist_path, const Mesh& m) {
  auto no = [](const std::string& what) { throw Error("microphysics: " + what + " has no verified port"); };
  const std::string mp = option(namelist_path, "config_microp_scheme");
  if (option(namelist_path, "config_physics_suite") != "none")
    no("config_physics_suite = '" + option(namelist_path, "config_physics_suite") + "'");
  if (mp != "mp_kessler" && mp != "off" && mp != "suite") no("config_microp_scheme = '" + mp + "'");
  if (option(namelist_path, "config_bucket_update") != "none") no("config_bucket_update (l_acrain, DRV:625)");
  // a restart of the Kessler path (i_rainnc, rainnc from the restart file, INIT:262-273) is not verified
  if (mp == "mp_kessler" && option(namelist_path, "config_do_restart") != "false") no("config_do_restart = true with Kessler");
  if (m.nVertLevels < 3) no("nVertLevels < 3 (the surface pressure of PI:950-964)");
}

void microphysics_init(const std::string& namelist_path, const Mesh& m, real dt, MicrophysicsFields& mp) {
  check_microphysics_supported(namelist_path, m);
  const std::string scheme = option(namelist_path, "config_microp_scheme");
  mp.active = scheme == "mp_kessler";
  if (!mp.active) return;
  const int nC = m.nCells + 1, nVertLevels = m.nVertLevels;
  mp.i_rainnc.allocate({nC}, 0);
  mp.rainnc.allocate({nC}, 0.0);
  mp.rainncv.allocate({nC}, 0.0);
  mp.precipw.allocate({nC}, 0.0);
  mp.dtheta_dt_mp.allocate({nVertLevels, nC}, 0.0);
  mp.tend_sfc_pressure.allocate({nC}, 0.0);
  mp.dt_dyn = dt;               // MGR:707
  mp.dt_microp = mp.dt_dyn;     // MGR:708 (n_microp = 1: not Thompson, MGR:710-716)
  const auto i_rainnc = mp.i_rainnc.view();
  parallel_for(1, m.nCellsSolve(), [=] CPPMPAS_LAMBDA(int iCell) { i_rainnc(iCell) = 0; });  // INIT:268-273
  const real config_microphysics_top = std::strtod(option(namelist_path, "config_microphysics_top").c_str(), nullptr);
  mp.mp_top_level = 0;  // INIT:276-282
  real layer_height = 0.5 / m.rdzw(1);
  if (layer_height <= config_microphysics_top) mp.mp_top_level = 1;
  for (int k = 2; k <= nVertLevels; ++k) {
    layer_height = layer_height + m.dzu(k);
    if (layer_height <= config_microphysics_top) mp.mp_top_level = k;
  }
  if (mp.mp_top_level != nVertLevels)  // PI:938-946 (kte > mp_top_level) is not taken by the reference
    throw Error("microphysics: mp_top_level = " + std::to_string(mp.mp_top_level) + " < nVertLevels (PI:938-946) has "
                "no verified port");
}

void driver_microphysics(const Mesh& m, State& s, FArray<real>& theta_m_a, FArray<real>& rho_zz_a,
                         FArray<real>& scalars_a, Diag& d, FArray<real>& rt_diabatic_tend_a, MicrophysicsFields& mp) {
  if (!mp.active) return;
  // module dims of physics_run_init (MGR:661-669): its..ite = 1..nCellsSolve, kts..kte = 1..nVertLevels,
  // kms..kme = 1..nVertLevels+1, jts = jte = 1
  const int ite = m.nCellsSolve(), kte = m.nVertLevels, kme = m.nVertLevels + 1;
  // as microphysics_init (check_microphysics_supported): the surface pressure of PI:950-964 reads levels 1-2 and
  // interface 3; a direct call with a thinner column is refused before any array is read
  if (kte < 3) throw Error("driver_microphysics: nVertLevels < 3 (the surface pressure of PI:950-964)");
  int index_qv = 0, index_qc = 0, index_qr = 0;
  for (std::size_t q = 0; q < s.scalar_names.size(); ++q) {
    if (s.scalar_names[q] == "qv") index_qv = static_cast<int>(q) + 1;
    if (s.scalar_names[q] == "qc") index_qc = static_cast<int>(q) + 1;
    if (s.scalar_names[q] == "qr") index_qr = static_cast<int>(q) + 1;
  }
  if (!index_qv || !index_qc || !index_qr) throw Error("driver_microphysics: Kessler needs qv, qc and qr");
  const auto zz = m.zz.cview(), zgrid = m.zgrid.cview();
  const auto exner = d.exner.view(), pressure_p = d.pressure_p.view(), rtheta_p = d.rtheta_p.view();
  const auto pressure_b = d.pressure_base.cview(), exner_b = d.exner_base.cview(), rtheta_b = d.rtheta_base.cview();
  const auto scalars = scalars_a.view(), theta_m = theta_m_a.view();
  const auto rho_zz = rho_zz_a.cview();
  const auto rt_diabatic_tend = rt_diabatic_tend_a.view();
  const auto mp_rainncv = mp.rainncv.view(), mp_rainnc = mp.rainnc.view(), mp_precipw = mp.precipw.view();
  const auto mp_dtheta_dt_mp = mp.dtheta_dt_mp.view(), mp_tend_sfc_pressure = mp.tend_sfc_pressure.view();
  const auto surface_pressure = s.surface_pressure.view();

  // allocate_microphysics (DRV:98-127): (ims:ime, kms:kme, jms:jme), i fastest
  FArray<real> rho_p_a({ite, kme}), th_p_a({ite, kme}), pi_p_a({ite, kme}), pres_p_a({ite, kme}), z_p_a({ite, kme}),
      dz_p_a({ite, kme}), qv_p_a({ite, kme}), qc_p_a({ite, kme}), qr_p_a({ite, kme}), rainnc_p_a({ite}),
      rainncv_p_a({ite});
  const auto rho_p = rho_p_a.view(), th_p = th_p_a.view(), pi_p = pi_p_a.view(), pres_p = pres_p_a.view(),
             z_p = z_p_a.view(), dz_p = dz_p_a.view(), qv_p = qv_p_a.view(), qc_p = qc_p_a.view(),
             qr_p = qr_p_a.view(), rainnc_p = rainnc_p_a.view(), rainncv_p = rainncv_p_a.view();
#define P(a, i, k) a(i, k)

  parallel_for(1, ite, [=] CPPMPAS_LAMBDA(int i) {  // precip_from_MPAS (DRV:530-541)
    rainncv_p[i - 1] = 0.0;
    rainnc_p[i - 1] = 0.0;
  });
  parallel_for(1, ite, [=] CPPMPAS_LAMBDA(int i) { mp_rainncv(i) = 0.0; });

  // microphysics_from_MPAS (PI:711-731); w_p (PI:731) is read only by Thompson and WSM6
  parallel_for(1, ite, [=] CPPMPAS_LAMBDA(int i) {
    for (int k = 1; k <= kte; ++k) {
      P(qv_p, i, k) = scalars(index_qv, k, i);
      P(qc_p, i, k) = scalars(index_qc, k, i);
      P(qr_p, i, k) = scalars(index_qr, k, i);
      P(rho_p, i, k) = zz(k, i) * rho_zz(k, i);
      P(th_p, i, k) = theta_m(k, i) / (1.0 + R_v / R_d * fmax0(P(qv_p, i, k)));
      P(pi_p, i, k) = exner(k, i);
      P(pres_p, i, k) = pressure_b(k, i) + pressure_p(k, i);
      P(z_p, i, k) = zgrid(k, i);
      P(dz_p, i, k) = zgrid(k + 1, i) - zgrid(k, i);
    }
  });

  {  // DRV:371-385
    using phys::kessler::View2;
    using phys::kessler::View3;
    if (kte > phys::kessler::kMaxLevels) throw Error("kessler: more than 128 levels");
    phys::kessler::KesslerKernel kk{};
    kk.t = View3<real>{th_p.data(), 0, ite};
    kk.qv = View3<real>{qv_p.data(), 0, ite};
    kk.qc = View3<real>{qc_p.data(), 0, ite};
    kk.qr = View3<real>{qr_p.data(), 0, ite};
    kk.rho = View3<const real>{rho_p.data(), 0, ite};
    kk.pii = View3<const real>{pi_p.data(), 0, ite};
    kk.dz8w = View3<const real>{dz_p.data(), 0, ite};
    kk.z = View3<const real>{z_p.data(), 0, ite};
    kk.rainnc = View2<real>{rainnc_p.data(), 0};
    kk.rainncv = View2<real>{rainncv_p.data(), 0};
    kk.dt = mp.dt_microp;
    kk.xlv = xlv;
    kk.cp = cp;
    kk.ep2 = ep_2;
    kk.svp1 = svp1;
    kk.svp2 = svp2;
    kk.svp3 = svp3;
    kk.svpt0 = svpt0;
    kk.rhowater = rho_w;
    kk.nz = kte;
    kk.i0 = 0;
    kk.ni = ite;
    kk.es_cap = true;  // MPAS's module_mp_kessler.F:219
    parallel_for(0, ite - 1, [=] CPPMPAS_LAMBDA(int i) { kk(i, 0); });  // one thread per column
  }

  // l_diags (compute_relhum, DRV:467-482) only at output times: relhum is not ported (no compared stream has it)
  parallel_for(1, ite, [=] CPPMPAS_LAMBDA(int i) {  // precip_to_MPAS (DRV:605-633), l_acrain = .false.
    mp_precipw(i) = 0.0;
  });
  parallel_for(1, ite, [=] CPPMPAS_LAMBDA(int i) {
    for (int k = 1; k <= kte; ++k) {
      const real rho_a = P(rho_p, i, k) / (1.0 + P(qv_p, i, k));
      mp_precipw(i) = mp_precipw(i) + P(qv_p, i, k) * rho_a * P(dz_p, i, k);
    }
    mp_rainncv(i) = rainnc_p[i - 1];
    mp_rainnc(i) = mp_rainnc(i) + mp_rainncv(i);
  });

  const int mp_top_level = mp.mp_top_level;  // microphysics_to_MPAS (PI:903-934)
  const real dt_dyn = mp.dt_dyn;
  parallel_for(1, ite, [=] CPPMPAS_LAMBDA(int i) {
    for (int k = 1; k <= mp_top_level; ++k) {
      rt_diabatic_tend(k, i) = theta_m(k, i);
      mp_dtheta_dt_mp(k, i) = theta_m(k, i) / (1.0 + rvord * scalars(index_qv, k, i));
      scalars(index_qv, k, i) = P(qv_p, i, k);
      scalars(index_qc, k, i) = P(qc_p, i, k);
      scalars(index_qr, k, i) = P(qr_p, i, k);
      theta_m(k, i) = P(th_p, i, k) * (1.0 + rvord * P(qv_p, i, k));
      rt_diabatic_tend(k, i) = (theta_m(k, i) - rt_diabatic_tend(k, i)) / dt_dyn;
      mp_dtheta_dt_mp(k, i) = (theta_m(k, i) / (1.0 + rvord * scalars(index_qv, k, i)) - mp_dtheta_dt_mp(k, i)) / (dt_dyn);
      rtheta_p(k, i) = rho_zz(k, i) * theta_m(k, i) - rtheta_b(k, i);
      exner(k, i) = libm::pow(zz(k, i) * (R_d / P0) * (rtheta_p(k, i) + rtheta_b(k, i)), rcv);
      pressure_p(k, i) = zz(k, i) * R_d * (exner(k, i) * rtheta_p(k, i) + (exner(k, i) - exner_b(k, i)) * rtheta_b(k, i));
    }
  });
  if (kte > mp_top_level) throw Error("microphysics_to_MPAS: kte > mp_top_level (PI:938-946) is not ported");

  parallel_for(1, ite, [=] CPPMPAS_LAMBDA(int i) {  // PI:950-964
    const real tem1 = zgrid(2, i) - zgrid(1, i);
    const real tem2 = zgrid(3, i) - zgrid(2, i);
    const real rho1 = rho_zz(1, i) * zz(1, i) * (1. + P(qv_p, i, 1));
    const real rho2 = rho_zz(2, i) * zz(2, i) * (1. + P(qv_p, i, 2));
    mp_tend_sfc_pressure(i) = surface_pressure(i);
    surface_pressure(i) = 0.5 * gravity * (zgrid(2, i) - zgrid(1, i)) * (rho1 - 0.5 * (rho2 - rho1) * tem1 / (tem1 + tem2));
    surface_pressure(i) = surface_pressure(i) + pressure_p(1, i) + pressure_b(1, i);
    mp_tend_sfc_pressure(i) = (surface_pressure(i) - mp_tend_sfc_pressure(i)) / dt_dyn;
  });
#undef P
}

}  // namespace cppmpas
