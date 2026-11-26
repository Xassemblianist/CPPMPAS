// SPDX-License-Identifier: GPL-3.0-or-later
// atm_srk3 and the routines it calls (src/core_atmosphere/dynamics/mpas_atm_time_integration.F, TI) on the path of
// the mountain-wave reference; see atm_dycore.hpp. Statement by statement: MPAS's loop order, operation order and
// parenthesisation (Fortran evaluates a*b*c/d as ((a*b)*c)/d, -a*b as -(a*b)), integer powers as multiplies.
// Every loop over cells / edges / vertices is a kernel (core/backend.hpp parallel_for): one thread per element, the k
// loop and the neighbour loops inside the thread in MPAS's order; the work arrays MPAS allocates per routine are the
// thread's own (Column), so the CPU and the GPU build run the same statements on the same values.
//
// Reference object code (objdump -d of bin/atmosphere_model, mpas-build asset 153c7206..., gfortran 14.2 -O3, no
// -march): no FMA and no ymm register in any of these routines; their only libm / libmvec calls are
//  * atm_recover_large_step_variables_work (0x125d40): the rk_step 3 k loop is vectorised, levels (1,2), (3,4), ...
//    call _ZGVbN2vv_pow (0x126ed3), an odd last level (or nVertLevels = 1) pow@plt (0x127016);
//  * atm_compute_dyn_tend_work (0x109430): the reference_u of TI:6547, in the j loop over edgesOnEdge unrolled by 2
//    (see coriolis_reference below): sincos@plt (0x10b2cf, 0x10b69d, 0x10b905), _ZGVbN2v_cos / _ZGVbN2v_sin
//    (0x10b543-0x10b575).
#include "core_atmosphere/atm_dycore.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "core/error.hpp"
#include "core/libm.hpp"
#include "core_atmosphere/atm_dissipation.hpp"
#include "core_atmosphere/atm_views.hpp"
#include "operators/vector_reconstruct.hpp"

namespace cppmpas {

namespace {

// src/framework/mpas_constants.F:43-56 (folded with one correctly rounded operation per node, as gfortran does)
constexpr real rgas = 287.0;
constexpr real rv = 461.6;
constexpr real cp = 7.0 * rgas / 2.0;
constexpr real rvord = rv / rgas;
constexpr real gravity = 9.80616;
constexpr real prandtl = 1.0;

// flux4 / flux3 statement functions (TI:6303-6308): ua*(...)/12.0 is (ua*(...))/12.0
CPPMPAS_HD inline real flux4(real q_im2, real q_im1, real q_i, real q_ip1, real ua) {
  return ua * (7. * (q_i + q_im1) - (q_ip1 + q_im2)) / 12.0;
}
CPPMPAS_HD inline real flux3(real q_im2, real q_im1, real q_i, real q_ip1, real ua, real coef3) {
  return flux4(q_im2, q_im1, q_i, q_ip1, ua) + coef3 * std::fabs(ua) * ((q_ip1 - q_im2) - 3. * (q_i - q_im1)) / 12.0;
}

// sign(1.0_RKIND, x) is copysign(1.0, x) in gfortran (-0.0 gives -1.0)
CPPMPAS_HD inline real sign1(real x) { return std::copysign(1.0, x); }

}  // namespace

// ---------------------------------------------------------------------------------------------------------------
DycoreConfig read_dycore_config(const std::string& namelist_path, const InitConfig& init_cfg) {
  const auto kv = effective_namelist(namelist_path);
  auto get = [&](const char* k) -> const std::string& {
    const auto it = kv.find(k);
    if (it == kv.end()) throw Error(std::string("read_dycore_config: option ") + k + " is not in the Registry table");
    return it->second;
  };
  auto r = [&](const char* k) { return std::strtod(get(k).c_str(), nullptr); };
  auto i = [&](const char* k) { return static_cast<int>(std::strtol(get(k).c_str(), nullptr, 10)); };
  auto l = [&](const char* k) { return get(k) == "true"; };
  DycoreConfig c;
  c.config_time_integration = get("config_time_integration");
  c.config_time_integration_order = i("config_time_integration_order");
  c.config_dt = init_cfg.config_dt;
  c.config_split_dynamics_transport = l("config_split_dynamics_transport");
  c.config_number_of_sub_steps = i("config_number_of_sub_steps");
  c.config_dynamics_split_steps = i("config_dynamics_split_steps");
  c.config_h_mom_eddy_visc2 = r("config_h_mom_eddy_visc2");
  c.config_h_mom_eddy_visc4 = r("config_h_mom_eddy_visc4");
  c.config_v_mom_eddy_visc2 = r("config_v_mom_eddy_visc2");
  c.config_h_theta_eddy_visc2 = r("config_h_theta_eddy_visc2");
  c.config_h_theta_eddy_visc4 = r("config_h_theta_eddy_visc4");
  c.config_v_theta_eddy_visc2 = r("config_v_theta_eddy_visc2");
  c.config_horiz_mixing = get("config_horiz_mixing");
  c.config_les_model = get("config_les_model");
  c.config_len_disp = init_cfg.config_len_disp;  // as atm_core_init left it (mpas_atm_core.F:155-201)
  c.config_mix_scalars = l("config_mix_scalars");
  c.config_del4u_div_factor = r("config_del4u_div_factor");
  c.config_scalar_advection = l("config_scalar_advection");
  c.config_positive_definite = l("config_positive_definite");
  c.config_monotonic = l("config_monotonic");
  c.config_coef_3rd_order = init_cfg.config_coef_3rd_order;
  c.config_mix_full = l("config_mix_full");
  c.config_smdiv = r("config_smdiv");
  c.config_apvm_upwinding = init_cfg.config_apvm_upwinding;
  c.config_mpas_cam_coef = r("config_mpas_cam_coef");
  c.config_rayleigh_damp_u = l("config_rayleigh_damp_u");
  c.config_apply_lbcs = l("config_apply_lbcs");
  c.config_IAU_option = get("config_iau_option");
  c.config_physics_suite = get("config_physics_suite");
  c.config_microp_scheme = get("config_microp_scheme");
  c.config_convection_scheme = get("config_convection_scheme");
  c.config_pbl_scheme = get("config_pbl_scheme");
  c.config_radt_lw_scheme = get("config_radt_lw_scheme");
  c.config_radt_sw_scheme = get("config_radt_sw_scheme");
  return c;
}

void check_dycore_supported(const DycoreConfig& c, const Mesh& m, const State& s) {
  auto no = [](const std::string& what) { throw Error("dycore: " + what + " has no verified port"); };
  if (c.config_time_integration != "SRK3") no("config_time_integration = '" + c.config_time_integration + "'");  // TI:1828
  if (c.config_time_integration_order != 2)  // TI:2066, 2210 (order 3 not verified)
    no("config_time_integration_order = " + std::to_string(c.config_time_integration_order));
  if (c.config_split_dynamics_transport) no("config_split_dynamics_transport = true");  // TI:2054, 2674
  if (c.config_scalar_advection) no("config_scalar_advection = true (atm_advance_scalars, TI:2504-2507)");
  if (c.config_positive_definite) no("config_positive_definite = true");
  if (c.config_horiz_mixing != "2d_fixed") no("config_horiz_mixing = '" + c.config_horiz_mixing + "'");  // TI:6346
  if (c.config_les_model != "none") no("config_les_model = '" + c.config_les_model + "'");
  if (c.config_mpas_cam_coef > 0.0) no("config_mpas_cam_coef > 0 (TI:6390)");
  if (c.config_rayleigh_damp_u) no("config_rayleigh_damp_u = true (TI:6608)");
  if (c.config_mix_scalars) no("config_mix_scalars = true (DM:1334)");
  if (!c.config_mix_full) no("config_mix_full = false (DM:848, 1442)");
  if (c.config_h_mom_eddy_visc4 > 0.0) no("config_h_mom_eddy_visc4 > 0 (DM:728, 1073)");
  if (c.config_h_theta_eddy_visc4 > 0.0) no("config_h_theta_eddy_visc4 > 0 (DM:1306)");
  if (!(c.config_v_mom_eddy_visc2 > 0.0)) no("config_v_mom_eddy_visc2 <= 0 (DM:816, 1102 not taken)");
  if (!(c.config_v_theta_eddy_visc2 > 0.0)) no("config_v_theta_eddy_visc2 <= 0 (DM:1414 not taken)");
  if (c.config_apply_lbcs) no("config_apply_lbcs = true");
  if (c.config_IAU_option != "off") no("config_IAU_option = '" + c.config_IAU_option + "' (TI:2168)");
  if (c.config_physics_suite != "none") no("config_physics_suite = '" + c.config_physics_suite + "'");
  // mpas_atmphys_control.F:157-165: with the suite 'none', 'suite' is 'off'
  const auto off = [](const std::string& x) { return x == "suite" || x == "off"; };
  if (!off(c.config_microp_scheme)) no("config_microp_scheme = '" + c.config_microp_scheme + "' (TI:2806)");
  if (!off(c.config_convection_scheme)) no("config_convection_scheme = '" + c.config_convection_scheme + "'");
  if (!off(c.config_pbl_scheme)) no("config_pbl_scheme = '" + c.config_pbl_scheme + "'");
  if (!off(c.config_radt_lw_scheme) || !off(c.config_radt_sw_scheme)) no("radiation schemes");
  // the clock's dt is config_dt through a 9-decimal string (mpas_timekeeping.F:1449-1498): exact for whole seconds
  if (!(c.config_dt > 0.0) || c.config_dt != std::floor(c.config_dt) || c.config_dt > 86400.0)
    no("config_dt that is not a whole number of seconds in 1..86400");
  if (c.config_number_of_sub_steps < 2) no("config_number_of_sub_steps < 2");
  if (s.scalar_names.size() != 1 || s.scalar_names[0] != "qv") no("scalars other than qv alone");
  // odd nVertLevels: the last level of the rk_step-3 Exner pow is scalar pow@plt (TI:4439), never exercised
  if (m.nVertLevels < 4 || m.nVertLevels % 2 != 0) no("nVertLevels = " + std::to_string(m.nVertLevels) +
                                                      " (only an even nVertLevels >= 4 is verified)");
  // the reference_u term (TI:6546-6550) calls sincos / libmvec lanes whose values only enter through fEdge: with
  // fEdge = 0 (MW) the term is +-0; any other fEdge is not verified
  for (int iEdge = 1; iEdge <= m.nEdges; ++iEdge)
    if (m.fEdge(iEdge) != 0.0) no("fEdge != 0 (the reference_u term of TI:6547-6548)");
  if (m.vertexDegree != 3) no("vertexDegree != 3");
}

// ---------------------------------------------------------------------------------------------------------------
void allocate_dycore_fields(const Mesh& m, const State& s, DycoreFields& f) {
  const int nC = m.nCells + 1, nE = m.nEdges + 1, nV = m.nVertices + 1, nL = m.nVertLevels, nL1 = nL + 1;
  if (nL < 1) throw Error("allocate_dycore_fields: no vertical levels");
  const int ns = s.scalars.extent(0);
  f.u2.allocate({nL, nE}, 0.0);
  f.w2.allocate({nL1, nC}, 0.0);
  f.theta_m2.allocate({nL, nC}, 0.0);
  f.rho_zz2.allocate({nL, nC}, 0.0);
  f.scalars2.allocate({ns, nL, nC}, 0.0);
  for (auto* a : {&f.ru_save, &f.cqu, &f.ru_p, &f.ruAvg, &f.ruAvg_split, &f.tend_u, &f.tend_u_euler,
                  &f.tend_ru_physics})
    a->allocate({nL, nE}, 0.0);
  for (auto* a : {&f.rtheta_p_save, &f.rho_p_save, &f.rho_zz_old_split, &f.cqw, &f.cofwr, &f.cofwz, &f.cofwt,
                  &f.a_tri, &f.alpha_tri, &f.gamma_tri, &f.rtheta_pp, &f.rtheta_pp_old, &f.rho_pp, &f.h_divergence,
                  &f.eddy_visc_horz, &f.tend_theta, &f.tend_rho, &f.tend_theta_euler, &f.rt_diabatic_tend,
                  &f.rthdynten, &f.tend_rtheta_physics, &f.tend_rho_physics})
    a->allocate({nL, nC}, 0.0);
  for (auto* a : {&f.rw_save, &f.coftz, &f.rw_p, &f.wwAvg, &f.wwAvg_split, &f.tend_w, &f.tend_w_euler})
    a->allocate({nL1, nC}, 0.0);
  f.cofrz.allocate({nL}, 0.0);
  f.tend_scalars.allocate({ns, nL, nC}, 0.0);
  // module arrays (TI:2023, 2227-2234): allocated (uninitialised in MPAS, every element read is written first)
  for (auto* a : {&f.qtot, &f.dpdz, &f.delsq_theta, &f.delsq_w, &f.delsq_divergence}) a->allocate({nL, nC}, 0.0);
  f.delsq_u.allocate({nL, nE}, 0.0);
  f.delsq_vorticity.allocate({nL, nV}, 0.0);
  allocate_solve_diag_work(m, f.sdw);
  f.index_qv = 0;
  for (std::size_t k = 0; k < s.scalar_names.size(); ++k)
    if (s.scalar_names[k] == "qv") f.index_qv = static_cast<int>(k) + 1;
  if (f.index_qv == 0) throw Error("allocate_dycore_fields: qv is not an active scalar");
  // array_group "moist" (Registry.xml:1701-1748): qv, qc, qr, ... in Registry order, all of the active scalars
  f.moist_start = 1;
  f.moist_end = ns;
}

// ---------------------------------------------------------------------------------------------------------------
// TI:3088-3185
void atm_rk_integration_setup(const Mesh& m, State& s, Diag& d, DycoreFields& f) {
  const int nVertLevels = m.nVertLevels, num_scalars = s.scalars.extent(0);
  parallel_for(1, nVertLevels, [=, m = view(m), f = view(f)] CPPMPAS_LAMBDA(int k) {
    f.theta_m2(k, m.nCells + 1) = 0.0;  // :3139
  });
  parallel_for(1, m.nEdges, [=, s = view(s), d = view(d), f = view(f)] CPPMPAS_LAMBDA(int iEdge) {  // :3144-3150
    for (int k = 1; k <= nVertLevels; ++k) {
      f.ru_save(k, iEdge) = d.ru(k, iEdge);
      f.u2(k, iEdge) = s.u(k, iEdge);
    }
  });
  parallel_for(1, m.nCells, [=, s = view(s), d = view(d), f = view(f)] CPPMPAS_LAMBDA(int iCell) {  // :3153-3162
    for (int k = 1; k <= nVertLevels; ++k) {
      f.rtheta_p_save(k, iCell) = d.rtheta_p(k, iCell);
      f.rho_p_save(k, iCell) = d.rho_p(k, iCell);
      f.theta_m2(k, iCell) = s.theta_m(k, iCell);
      f.rho_zz2(k, iCell) = s.rho_zz(k, iCell);
      f.rho_zz_old_split(k, iCell) = s.rho_zz(k, iCell);
    }
  });
  parallel_for(1, m.nCells, [=, s = view(s), d = view(d), f = view(f)] CPPMPAS_LAMBDA(int iCell) {  // :3165-3171
    for (int k = 1; k <= nVertLevels + 1; ++k) {
      f.rw_save(k, iCell) = d.rw(k, iCell);
      f.w2(k, iCell) = s.w(k, iCell);
    }
  });
  parallel_for(1, m.nCells, [=, s = view(s), f = view(f)] CPPMPAS_LAMBDA(int iCell) {  // :3174-3181
    for (int k = 1; k <= nVertLevels; ++k)
      for (int j = 1; j <= num_scalars; ++j) f.scalars2(j, k, iCell) = s.scalars(j, k, iCell);
  });
}

// TI:3188-3283 (scalars of time level 2)
void atm_compute_moist_coefficients(const Mesh& m, DycoreFields& f) {
  const int nVertLevels = m.nVertLevels;
  const auto scalars = f.scalars2.cview();
  parallel_for(1, m.nCells, [=, f = view(f)] CPPMPAS_LAMBDA(int iCell) {  // :3237-3246
    for (int k = 1; k <= nVertLevels; ++k) {
      f.qtot(k, iCell) = 0.0;
      for (int iq = f.moist_start; iq <= f.moist_end; ++iq) f.qtot(k, iCell) = f.qtot(k, iCell) + scalars(iq, k, iCell);
    }
  });
  parallel_for(1, m.nCells, [=, f = view(f)] CPPMPAS_LAMBDA(int iCell) {  // :3252-3258
    for (int k = 2; k <= nVertLevels; ++k) {
      const real qtotal = 0.5 * (f.qtot(k, iCell) + f.qtot(k - 1, iCell));
      f.cqw(k, iCell) = 1.0 / (1.0 + qtotal);
    }
  });
  parallel_for(1, m.nEdges, [=, m = view(m), f = view(f)] CPPMPAS_LAMBDA(int iEdge) {  // :3265-3279 (cell1 <= nCellsSolve always true, one rank)
    const int cell1 = m.cellsOnEdge(1, iEdge);
    const int cell2 = m.cellsOnEdge(2, iEdge);
    for (int k = 1; k <= nVertLevels; ++k) {
      real qtotal = 0.0;
      for (int iq = f.moist_start; iq <= f.moist_end; ++iq)
        qtotal = qtotal + 0.5 * (scalars(iq, k, cell1) + scalars(iq, k, cell2));
      f.cqu(k, iEdge) = 1.0 / (1.0 + qtotal);
    }
  });
}

// physics_get_tend with every scheme off (src/core_atmosphere/physics/mpas_atmphys_todynamics.F:98-279, work routine
// :304-498): the tendencies are zeroed (:211-216), no scheme block runs, then :489-495 (theta_m and scalars of time
// level 1). R_v/R_d is the folded quotient (mpas_atmphys_constants.F:11, 42: R_d = rgas, R_v = 461.6).
void physics_get_tend(const Mesh& m, const State& s, DycoreFields& f) {
  const int nVertLevels = m.nVertLevels;
  constexpr real R_v_R_d = rv / rgas;
  fill_all(f.tend_scalars, 0.0);         // :213
  fill_all(f.tend_ru_physics, 0.0);      // :214
  fill_all(f.tend_rtheta_physics, 0.0);  // :215
  fill_all(f.tend_rho_physics, 0.0);     // :216
  parallel_for(1, m.nCellsSolve(), [=, s = view(s), f = view(f)] CPPMPAS_LAMBDA(int i) {  // :489-495, tend_th = 0 (:210-211)
    for (int k = 1; k <= nVertLevels; ++k) {
      const real coeff = (1. + R_v_R_d * s.scalars(f.index_qv, k, i));
      real tend_th = 0.0;
      tend_th = coeff * tend_th + R_v_R_d * s.theta_m(k, i) * f.tend_scalars(f.index_qv, k, i) / coeff;
      f.tend_rtheta_physics(k, i) = f.tend_rtheta_physics(k, i) + tend_th;
    }
  });
}

// TI:3286-3507 (t = theta_m of time level 2, dts = rk_sub_timestep(1))
void atm_compute_vert_imp_coefs(const Mesh& m, const State& s, const Diag& d, DycoreFields& f, real dts) {
  const int nVertLevels = m.nVertLevels;
  fits_column(nVertLevels, "atm_compute_vert_imp_coefs");
  const auto zz = m.zz.cview(), fzm = m.fzm.cview(), fzp = m.fzp.cview(), rdzu = m.rdzu.cview(), rdzw = m.rdzw.cview();
  const auto etp = m.etp.cview(), ewp = m.ewp.cview();
  const auto p = d.exner.cview(), pb = d.exner_base.cview(), rt = d.rtheta_p.cview(), rtb = d.rtheta_base.cview(),
             rb = s.rho_base.cview();
  const auto t = f.theta_m2.cview();
  auto cofwr = f.cofwr.view(), cofwz = f.cofwz.view(), coftz = f.coftz.view(), cofwt = f.cofwt.view(),
       cofrz = f.cofrz.view();
  auto a_tri = f.a_tri.view(), alpha_tri = f.alpha_tri.view(), gamma_tri = f.gamma_tri.view(), cqw = f.cqw.view();
  constexpr real rcv = rgas / (cp - rgas);  // :3422
  constexpr real c2 = cp * rcv;             // :3423
  parallel_for(1, nVertLevels, [=] CPPMPAS_LAMBDA(int k) {
    cofrz(k) = rdzw(k);  // :3428-3430
  });
  parallel_for(1, m.nCellsSolve(), [=, f = view(f)] CPPMPAS_LAMBDA(int iCell) {    // :3435-3502
    Column<real> b_tri(nVertLevels + 1), c_tri(nVertLevels + 1);
    for (int k = 2; k <= nVertLevels; ++k) cofwr(k, iCell) = .5 * gravity * (fzm(k) * zz(k, iCell) + fzp(k) * zz(k - 1, iCell));
    coftz(1, iCell) = 0.0;
    for (int k = 2; k <= nVertLevels; ++k) {  // :3445-3449
      cofwz(k, iCell) = c2 * (fzm(k) * zz(k, iCell) + fzp(k) * zz(k - 1, iCell)) * rdzu(k) * cqw(k, iCell) *
                        (fzm(k) * p(k, iCell) + fzp(k) * p(k - 1, iCell));
      coftz(k, iCell) = (fzm(k) * t(k, iCell) + fzp(k) * t(k - 1, iCell));
    }
    coftz(nVertLevels + 1, iCell) = 0.0;
    for (int k = 1; k <= nVertLevels; ++k) {  // :3453-3465
      const real qtotal = f.qtot(k, iCell);
      cofwt(k, iCell) = .5 * rcv * zz(k, iCell) * gravity * rb(k, iCell) / (1. + qtotal) * p(k, iCell) /
                        ((rtb(k, iCell) + rt(k, iCell)) * pb(k, iCell));
    }
    a_tri(1, iCell) = 0.;  // :3467-3471
    b_tri[1] = 1.;
    c_tri[1] = 0.;
    gamma_tri(1, iCell) = 0.;
    alpha_tri(1, iCell) = 0.;
    for (int k = 2; k <= nVertLevels; ++k) {  // :3475-3493
      a_tri(k, iCell) = -(cofwz(k, iCell) * coftz(k - 1, iCell) * rdzw(k - 1) * zz(k - 1, iCell)) +
                        cofwr(k, iCell) * cofrz(k - 1) - cofwt(k - 1, iCell) * coftz(k - 1, iCell) * rdzw(k - 1);
      a_tri(k, iCell) = a_tri(k, iCell) * etp(k - 1) * ewp(k - 1);
      b_tri[k] = +cofwz(k, iCell) * coftz(k, iCell) *
                     (etp(k) * rdzw(k) * zz(k, iCell) + etp(k - 1) * rdzw(k - 1) * zz(k - 1, iCell)) -
                 coftz(k, iCell) * (etp(k) * cofwt(k, iCell) * rdzw(k) - etp(k - 1) * cofwt(k - 1, iCell) * rdzw(k - 1)) +
                 cofwr(k, iCell) * (etp(k) * cofrz(k) - etp(k - 1) * cofrz(k - 1));
      b_tri[k] = b_tri[k] * ewp(k);
      c_tri[k] = -(cofwz(k, iCell) * coftz(k + 1, iCell) * rdzw(k) * zz(k, iCell)) - cofwr(k, iCell) * cofrz(k) +
                 cofwt(k, iCell) * coftz(k + 1, iCell) * rdzw(k);
      c_tri[k] = c_tri[k] * etp(k) * ewp(k + 1);
    }
    c_tri[nVertLevels] = 0.0;  // :3494
    for (int k = 2; k <= nVertLevels; ++k) {  // :3497-3500
      alpha_tri(k, iCell) = 1. / (1.0 + (dts * dts) * (b_tri[k] - a_tri(k, iCell) * gamma_tri(k - 1, iCell)));
      gamma_tri(k, iCell) = (dts * dts) * c_tri[k] * alpha_tri(k, iCell);
    }
  });
}

// ---------------------------------------------------------------------------------------------------------------
namespace {

// reference_u = u_init(k)*cos(angleEdge(eoe)) - v_init(k)*sin(angleEdge(eoe)) (TI:6547) with the sin / cos the
// reference calls. atm_compute_dyn_tend_work's j loop (TI:6542) is unrolled by 2: for nEdgesOnEdge >= 3 the pairs
// (1,2), (3,4), ..., ((n-3)/2*2+1, (n-3)/2*2+2) (loop 0x10b2a0-0x10b80b), then the rest one j at a time
// (0x10b8d7-0x10ba01); nEdgesOnEdge <= 2 takes only the second loop (0x10b1df -> 0x10f218). The first j of a pair
// and every j of the second loop use one sincos@plt per j (0x10b2cf, 0x10b905). The second j of a pair: for
// nVertLevels >= 4 (0x10b2e9) the k loop runs 4 levels per iteration over k = 1 .. 4*(nVertLevels/4) calling
// _ZGVbN2v_cos and _ZGVbN2v_sin with angleEdge(eoe) in both lanes (0x10b543-0x10b575), the remaining levels (and
// all levels for nVertLevels < 4) use one sincos@plt (0x10b69d).
struct CoriolisAngles {
  real c_scalar = 0, s_scalar = 0, c_vec = 0, s_vec = 0;
  int k_vec_end = 0;  // levels 1..k_vec_end use c_vec / s_vec
};

}  // namespace

void atm_compute_dyn_tend(const DycoreConfig& cfg, const Mesh& m, const State& s, Diag& d, DycoreFields& f,
                          int rk_step, real dt) {
  const int nVertLevels = m.nVertLevels;
  fits_column(nVertLevels, "atm_compute_dyn_tend");
  const auto rho_zz = f.rho_zz2.cview(), u = f.u2.cview(), w = f.w2.cview(), theta_m = f.theta_m2.cview();
  const auto theta_m_save = s.theta_m.cview();  // time level 1 (TI:5934)
  const auto rb = s.rho_base.cview(), rr_save = f.rho_p_save.cview(), ru = d.ru.cview(), ru_save = f.ru_save.cview(),
             rw = d.rw.cview();
  const auto rw_save = f.rw_save.cview(), ke = d.ke.cview(), pv_edge = d.pv_edge.cview(), pp = d.pressure_p.cview(),
             rho_edge = d.rho_edge.cview();
  const auto zz = m.zz.cview(), fzm = m.fzm.cview(), fzp = m.fzp.cview(), rdzw = m.rdzw.cview(), rdzu = m.rdzu.cview();
  auto tend_u = f.tend_u.view(), tend_w = f.tend_w.view(), tend_theta = f.tend_theta.view(),
       tend_rho = f.tend_rho.view();
  auto tend_u_euler = f.tend_u_euler.view(), tend_w_euler = f.tend_w_euler.view(),
       tend_theta_euler = f.tend_theta_euler.view();
  auto h_divergence = f.h_divergence.view(), dpdz = f.dpdz.view(), eddy_visc_horz = f.eddy_visc_horz.view();
  const real coef_3rd_order = cfg.config_coef_3rd_order;

  const real prandtl_inv = 1.0 / prandtl;  // :6322
  const real invDt = 1.0 / dt;              // :6323 (read by smagorinsky_2d only)
  (void)invDt;
  const real v_mom_eddy_visc2 = cfg.config_v_mom_eddy_visc2;      // :6326
  const real v_theta_eddy_visc2 = cfg.config_v_theta_eddy_visc2;  // :6327
  real h_mom_eddy_visc4 = cfg.config_h_mom_eddy_visc4, h_theta_eddy_visc4 = cfg.config_h_theta_eddy_visc4;

  if (rk_step == 1) {  // :6329-6410
    parallel_for(1, m.nEdges, [=] CPPMPAS_LAMBDA(int iEdge) {
      for (int k = 1; k <= nVertLevels; ++k) tend_u_euler(k, iEdge) = 0.0;
    });
    if (cfg.config_les_model != "none") throw Error("atm_compute_dyn_tend: LES models (TI:6371-6388) are not ported");
    if (cfg.config_horiz_mixing == "2d_fixed") {  // :6354-6367
      parallel_for(1, m.nCells, [=, cfg = view(cfg)] CPPMPAS_LAMBDA(int iCell) {
        for (int k = 1; k <= nVertLevels; ++k) eddy_visc_horz(k, iCell) = cfg.config_h_theta_eddy_visc2;
      });
      h_mom_eddy_visc4 = cfg.config_h_mom_eddy_visc4;
      h_theta_eddy_visc4 = cfg.config_h_theta_eddy_visc4;
    } else {
      throw Error("atm_compute_dyn_tend: config_horiz_mixing = '" + cfg.config_horiz_mixing + "' is not ported");
    }
    if (cfg.config_mpas_cam_coef > 0.0) throw Error("atm_compute_dyn_tend: the CAM top filter (TI:6390) is not ported");
  }

  parallel_for(1, m.nCells, [=, m = view(m)] CPPMPAS_LAMBDA(int iCell) {  // :6417-6434
    for (int k = 1; k <= nVertLevels; ++k) h_divergence(k, iCell) = 0.0;
    for (int i = 1; i <= m.nEdgesOnCell(iCell); ++i) {
      const int iEdge = m.edgesOnCell(i, iCell);
      const real edge_sign = m.edgesOnCell_sign(i, iCell) * m.dvEdge(iEdge);
      for (int k = 1; k <= nVertLevels; ++k) h_divergence(k, iCell) = h_divergence(k, iCell) + edge_sign * ru(k, iEdge);
    }
  });
  parallel_for(1, m.nCells, [=, m = view(m)] CPPMPAS_LAMBDA(int iCell) {  // :6441-6447
    const real r = m.invAreaCell(iCell);
    for (int k = 1; k <= nVertLevels; ++k) h_divergence(k, iCell) = h_divergence(k, iCell) * r;
  });

  if (rk_step == 1) {  // :6455-6472
    parallel_for(1, m.nCells, [=, f = view(f)] CPPMPAS_LAMBDA(int iCell) {
      for (int k = 1; k <= nVertLevels; ++k) {
        tend_rho(k, iCell) = -h_divergence(k, iCell) - rdzw(k) * (rw(k + 1, iCell) - rw(k, iCell)) +
                             f.tend_rho_physics(k, iCell);
        dpdz(k, iCell) = -(gravity * (rb(k, iCell) * (f.qtot(k, iCell)) + rr_save(k, iCell) * (1. + f.qtot(k, iCell))));
      }
    });
  }

  parallel_for(1, m.nEdgesSolve(), [=, m = view(m), s = view(s), f = view(f)] CPPMPAS_LAMBDA(int iEdge) {  // :6482-6573
    Column<real> wduz(nVertLevels + 2), q(nVertLevels + 1);
    const int cell1 = m.cellsOnEdge(1, iEdge);
    const int cell2 = m.cellsOnEdge(2, iEdge);
    if (rk_step == 1) {  // :6489-6497
      for (int k = 1; k <= nVertLevels; ++k)
        tend_u_euler(k, iEdge) =
            -(f.cqu(k, iEdge) * ((pp(k, cell2) - pp(k, cell1)) * m.invDcEdge(iEdge) / (.5 * (zz(k, cell2) + zz(k, cell1))) -
                                 0.5 * m.zxu(k, iEdge) * (dpdz(k, cell1) + dpdz(k, cell2))));
    }
    wduz[1] = 0.;  // :6501-6512
    int k = 2;
    wduz[k] = 0.5 * (rw(k, cell1) + rw(k, cell2)) * (fzm(k) * u(k, iEdge) + fzp(k) * u(k - 1, iEdge));
    for (k = 3; k <= nVertLevels - 1; ++k)
      wduz[k] = flux3(u(k - 2, iEdge), u(k - 1, iEdge), u(k, iEdge), u(k + 1, iEdge), 0.5 * (rw(k, cell1) + rw(k, cell2)),
                      1.0);
    k = nVertLevels;
    wduz[k] = 0.5 * (rw(k, cell1) + rw(k, cell2)) * (fzm(k) * u(k, iEdge) + fzp(k) * u(k - 1, iEdge));
    wduz[nVertLevels + 1] = 0.;
    for (k = 1; k <= nVertLevels; ++k) tend_u(k, iEdge) = -(rdzw(k) * (wduz[k + 1] - wduz[k]));  // :6516-6518

    for (k = 1; k <= nVertLevels; ++k) q[k] = 0.0;  // :6523-6538
    for (int j = 1; j <= m.nEdgesOnEdge(iEdge); ++j) {
      const int eoe = m.edgesOnEdge(j, iEdge);
      for (k = 1; k <= nVertLevels; ++k) {
        const real workpv = 0.5 * (pv_edge(k, iEdge) + pv_edge(k, eoe));
        q[k] = q[k] + m.weightsOnEdge(j, iEdge) * u(k, eoe) * workpv;
      }
    }
    {  // perturbation_coriolis = .true. (TI:6300), :6540-6552
      const int n = m.nEdgesOnEdge(iEdge);
      const int pairs = n >= 3 ? (n - 3) / 2 + 1 : 0;
      for (int j = 1; j <= n; ++j) {
        const int eoe = m.edgesOnEdge(j, iEdge);
        const real angle = m.angleEdge(eoe);
        CoriolisAngles ang;
        const bool second_of_pair = j <= 2 * pairs && j % 2 == 0;
        libm::sincos(angle, &ang.s_scalar, &ang.c_scalar);
        if (second_of_pair && nVertLevels >= 4) {
          ang.k_vec_end = 4 * (nVertLevels / 4);
          ang.c_vec = libm::mvec_cos(angle);
          ang.s_vec = libm::mvec_sin(angle);
        }
        for (k = 1; k <= nVertLevels; ++k) {
          const bool vec = k <= ang.k_vec_end;
          const real cosa = vec ? ang.c_vec : ang.c_scalar, sina = vec ? ang.s_vec : ang.s_scalar;
          const real reference_u = s.u_init(k) * cosa - s.v_init(k) * sina;
          q[k] = q[k] - m.weightsOnEdge(j, iEdge) * reference_u * m.fEdge(iEdge);
        }
      }
    }
    for (k = 1; k <= nVertLevels; ++k)  // :6556-6571 (CURVATURE is not defined in the reference build)
      tend_u(k, iEdge) = tend_u(k, iEdge) + rho_edge(k, iEdge) * (q[k] - (ke(k, cell2) - ke(k, cell1)) * m.invDcEdge(iEdge)) -
                         u(k, iEdge) * 0.5 * (h_divergence(k, cell1) + h_divergence(k, cell2));
  });

  DissipationArgs da;
  da.config_mix_full = cfg.config_mix_full;
  da.h_mom_eddy_visc4 = h_mom_eddy_visc4;
  da.v_mom_eddy_visc2 = v_mom_eddy_visc2;
  da.h_theta_eddy_visc4 = h_theta_eddy_visc4;
  da.v_theta_eddy_visc2 = v_theta_eddy_visc2;
  da.prandtl_inv = prandtl_inv;
  if (rk_step == 1)  // :6582-6601
    u_dissipation_3d(m, da, f.eddy_visc_horz, f.delsq_u, f.u2, d.divergence, d.vorticity, d.rho_edge, f.tend_u_euler);

  if (cfg.config_rayleigh_damp_u) throw Error("atm_compute_dyn_tend: Rayleigh damping of u (TI:6608) is not ported");
  parallel_for(1, m.nEdgesSolve(), [=, f = view(f)] CPPMPAS_LAMBDA(int iEdge) {  // :6635-6642
    for (int k = 1; k <= nVertLevels; ++k)
      tend_u(k, iEdge) = tend_u(k, iEdge) + tend_u_euler(k, iEdge) + f.tend_ru_physics(k, iEdge);
  });

  // ----------- rhs for w
  parallel_for(1, m.nCellsSolve(), [=, m = view(m)] CPPMPAS_LAMBDA(int iCell) {  // :6654-6697
    Column<real> ru_edge_w(nVertLevels + 1), flux_arr(nVertLevels + 1);
    for (int k = 1; k <= nVertLevels + 1; ++k) tend_w(k, iCell) = 0.0;
    for (int i = 1; i <= m.nEdgesOnCell(iCell); ++i) {
      const int iEdge = m.edgesOnCell(i, iCell);
      for (int k = 2; k <= nVertLevels; ++k) ru_edge_w[k] = fzm(k) * ru(k, iEdge) + fzp(k) * ru(k - 1, iEdge);
      for (int k = 1; k <= nVertLevels; ++k) flux_arr[k] = 0.0;
      for (int j = 1; j <= m.nAdvCellsForEdge(iEdge); ++j) {
        const int iAdvCell = m.advCellsForEdge(j, iEdge);
        for (int k = 2; k <= nVertLevels; ++k) {
          const real scalar_weight = m.adv_coefs(j, iEdge) + sign1(ru_edge_w[k]) * m.adv_coefs_3rd(j, iEdge);
          flux_arr[k] = flux_arr[k] + scalar_weight * w(k, iAdvCell);
        }
      }
      for (int k = 2; k <= nVertLevels; ++k)
        tend_w(k, iCell) = tend_w(k, iCell) - m.edgesOnCell_sign(i, iCell) * ru_edge_w[k] * flux_arr[k];
    }
  });
  if (rk_step == 1)  // :6725-6740
    w_dissipation_3d(m, da, f.delsq_w, f.w2, d.rho_edge, f.rho_zz2, f.eddy_visc_horz, f.tend_w_euler);

  parallel_for(1, m.nCellsSolve(), [=, m = view(m), f = view(f)] CPPMPAS_LAMBDA(int iCell) {  // :6755-6790
    Column<real> wdwz(nVertLevels + 2);
    wdwz[1] = 0.0;
    int k = 2;
    wdwz[k] = 0.25 * (rw(k, iCell) + rw(k - 1, iCell)) * (w(k, iCell) + w(k - 1, iCell));
    for (k = 3; k <= nVertLevels - 1; ++k)
      wdwz[k] = flux3(w(k - 2, iCell), w(k - 1, iCell), w(k, iCell), w(k + 1, iCell), 0.5 * (rw(k, iCell) + rw(k - 1, iCell)),
                      1.0);
    k = nVertLevels;
    wdwz[k] = 0.25 * (rw(k, iCell) + rw(k - 1, iCell)) * (w(k, iCell) + w(k - 1, iCell));
    wdwz[nVertLevels + 1] = 0.0;
    for (k = 2; k <= nVertLevels; ++k)
      tend_w(k, iCell) = tend_w(k, iCell) * m.invAreaCell(iCell) - rdzu(k) * (wdwz[k + 1] - wdwz[k]);
    if (rk_step == 1) {
      for (k = 2; k <= nVertLevels; ++k)
        tend_w_euler(k, iCell) =
            tend_w_euler(k, iCell) - f.cqw(k, iCell) * (rdzu(k) * (pp(k, iCell) - pp(k - 1, iCell)) -
                                                        (fzm(k) * dpdz(k, iCell) + fzp(k) * dpdz(k - 1, iCell)));
    }
  });
  parallel_for(1, m.nCellsSolve(), [=] CPPMPAS_LAMBDA(int iCell) {  // :6797-6803
    for (int k = 2; k <= nVertLevels; ++k) tend_w(k, iCell) = tend_w(k, iCell) + tend_w_euler(k, iCell);
  });

  // ----------- rhs for theta
  parallel_for(1, m.nCellsSolve(), [=, m = view(m)] CPPMPAS_LAMBDA(int iCell) {  // :6814-6848
    Column<real> flux_arr(nVertLevels + 1);
    for (int k = 1; k <= nVertLevels; ++k) tend_theta(k, iCell) = 0.0;
    for (int i = 1; i <= m.nEdgesOnCell(iCell); ++i) {
      const int iEdge = m.edgesOnCell(i, iCell);
      for (int k = 1; k <= nVertLevels; ++k) flux_arr[k] = 0.0;
      for (int j = 1; j <= m.nAdvCellsForEdge(iEdge); ++j) {
        const int iAdvCell = m.advCellsForEdge(j, iEdge);
        for (int k = 1; k <= nVertLevels; ++k) {
          const real scalar_weight = m.adv_coefs(j, iEdge) + sign1(ru(k, iEdge)) * m.adv_coefs_3rd(j, iEdge);
          flux_arr[k] = flux_arr[k] + scalar_weight * theta_m(k, iAdvCell);
        }
      }
      for (int k = 1; k <= nVertLevels; ++k)
        tend_theta(k, iCell) = tend_theta(k, iCell) - m.edgesOnCell_sign(i, iCell) * ru(k, iEdge) * flux_arr[k];
    }
  });
  if (rk_step > 1) {  // :6853-6874
    parallel_for(1, m.nCellsSolve(), [=, m = view(m)] CPPMPAS_LAMBDA(int iCell) {
      for (int i = 1; i <= m.nEdgesOnCell(iCell); ++i) {
        const int iEdge = m.edgesOnCell(i, iCell);
        const int cell1 = m.cellsOnEdge(1, iEdge);
        const int cell2 = m.cellsOnEdge(2, iEdge);
        for (int k = 1; k <= nVertLevels; ++k) {
          const real flux = m.edgesOnCell_sign(i, iCell) * m.dvEdge(iEdge) * (ru_save(k, iEdge) - ru(k, iEdge)) * 0.5 *
                            (theta_m_save(k, cell2) + theta_m_save(k, cell1));
          tend_theta(k, iCell) = tend_theta(k, iCell) - flux;
        }
      }
    });
  }
  if (rk_step == 1)  // :6881-6903
    scalar_dissipation_3d_les(m, da, f.delsq_theta, f.theta_m2, d.rho_edge, f.rho_zz2, f.eddy_visc_horz,
                              f.tend_theta_euler);

  parallel_for(1, m.nCellsSolve(), [=, m = view(m), f = view(f)] CPPMPAS_LAMBDA(int iCell) {  // :6913-6939
    Column<real> wdtz(nVertLevels + 2);
    wdtz[1] = 0.0;
    int k = 2;
    wdtz[k] = rw(k, iCell) * (fzm(k) * theta_m(k, iCell) + fzp(k) * theta_m(k - 1, iCell));
    wdtz[k] = wdtz[k] + (rw_save(k, iCell) - rw(k, iCell)) * (fzm(k) * theta_m_save(k, iCell) + fzp(k) * theta_m_save(k - 1, iCell));
    for (k = 3; k <= nVertLevels - 1; ++k) {
      wdtz[k] = flux3(theta_m(k - 2, iCell), theta_m(k - 1, iCell), theta_m(k, iCell), theta_m(k + 1, iCell), rw(k, iCell),
                      coef_3rd_order);
      wdtz[k] = wdtz[k] + (rw_save(k, iCell) - rw(k, iCell)) *
                              (fzm(k) * theta_m_save(k, iCell) + fzp(k) * theta_m_save(k - 1, iCell));
    }
    k = nVertLevels;
    wdtz[k] = rw_save(k, iCell) * (fzm(k) * theta_m(k, iCell) + fzp(k) * theta_m(k - 1, iCell));
    wdtz[nVertLevels + 1] = 0.0;
    for (k = 1; k <= nVertLevels; ++k) {
      tend_theta(k, iCell) = tend_theta(k, iCell) * m.invAreaCell(iCell) - rdzw(k) * (wdtz[k + 1] - wdtz[k]);
      f.rthdynten(k, iCell) = (tend_theta(k, iCell) - tend_rho(k, iCell) * theta_m(k, iCell)) / rho_zz(k, iCell);
      tend_theta(k, iCell) = tend_theta(k, iCell) + rho_zz(k, iCell) * f.rt_diabatic_tend(k, iCell);
    }
  });
  parallel_for(1, m.nCellsSolve(), [=, f = view(f)] CPPMPAS_LAMBDA(int iCell) {  // :6944-6950
    for (int k = 1; k <= nVertLevels; ++k)
      tend_theta(k, iCell) = tend_theta(k, iCell) + tend_theta_euler(k, iCell) + f.tend_rtheta_physics(k, iCell);
  });
}

// ---------------------------------------------------------------------------------------------------------------
// TI:3510-3642 (bdyMaskCell = 0 <= nRelaxZone for every cell)
void atm_set_smlstep_pert_variables(const Mesh& m, DycoreFields& f) {
  const int nVertLevels = m.nVertLevels;
  const auto u_tend = f.tend_u.cview();
  auto w_tend = f.tend_w.view();
  parallel_for(1, m.nCellsSolve(), [=, m = view(m)] CPPMPAS_LAMBDA(int iCell) {
    for (int i = 1; i <= m.nEdgesOnCell(iCell); ++i) {  // :3622-3631
      const int iEdge = m.edgesOnCell(i, iCell);
      for (int k = 2; k <= nVertLevels; ++k) {
        const real flux = m.edgesOnCell_sign(i, iCell) * (m.fzm(k) * u_tend(k, iEdge) + m.fzp(k) * u_tend(k - 1, iEdge));
        w_tend(k, iCell) = w_tend(k, iCell) -
                           (m.zb_cell(k, i, iCell) + sign1(u_tend(k, iEdge)) * m.zb3_cell(k, i, iCell)) * flux;
      }
    }
    for (int k = 2; k <= nVertLevels; ++k)  // :3634-3636
      w_tend(k, iCell) = (m.fzm(k) * m.zz(k, iCell) + m.fzp(k) * m.zz(k - 1, iCell)) * w_tend(k, iCell);
  });
}

// TI:3645-4109: rho_zz and w of time level 2, theta_m of time level 1 (TI:3713-3718); specZoneMask* = 0
void atm_advance_acoustic_step(const Mesh& m, const State& s, const Diag& d, DycoreFields& f, real dts,
                               int small_step) {
  const int nVertLevels = m.nVertLevels;
  fits_column(nVertLevels, "atm_advance_acoustic_step");
  const auto rho_zz = f.rho_zz2.cview(), w = f.w2.cview(), theta_m = s.theta_m.cview();
  const auto zz = m.zz.cview(), fzm = m.fzm.cview(), fzp = m.fzp.cview(), rdzw = m.rdzw.cview(), dss = m.dss.cview();
  const auto etp = m.etp.cview(), etm = m.etm.cview(), ewp = m.ewp.cview(), ewm = m.ewm.cview();
  const auto exner = d.exner.cview(), rw = d.rw.cview(), rw_save = f.rw_save.cview();
  const auto cqu = f.cqu.cview(), cofwt = f.cofwt.cview(), coftz = f.coftz.cview(), cofrz = f.cofrz.cview(),
             cofwr = f.cofwr.cview(), cofwz = f.cofwz.cview();
  const auto a_tri = f.a_tri.cview(), alpha_tri = f.alpha_tri.cview(), gamma_tri = f.gamma_tri.cview();
  const auto tend_ru = f.tend_u.cview(), tend_rho = f.tend_rho.cview(), tend_rt = f.tend_theta.cview(),
             tend_rw = f.tend_w.cview();
  auto ru_p = f.ru_p.view(), rw_p = f.rw_p.view(), rtheta_pp = f.rtheta_pp.view(),
       rtheta_pp_old = f.rtheta_pp_old.view();
  auto ruAvg = f.ruAvg.view(), wwAvg = f.wwAvg.view(), rho_pp = f.rho_pp.view();
  constexpr real rcv = rgas / (cp - rgas);  // :3876
  constexpr real c2 = cp * rcv;             // :3877

  if (small_step != 1) {  // :3881-3921
    parallel_for(1, m.nEdges, [=, m = view(m)] CPPMPAS_LAMBDA(int iEdge) {
      const int cell1 = m.cellsOnEdge(1, iEdge);
      const int cell2 = m.cellsOnEdge(2, iEdge);
      for (int k = 1; k <= nVertLevels; ++k) {
        real pgrad = ((rtheta_pp(k, cell2) - rtheta_pp(k, cell1)) * m.invDcEdge(iEdge)) / (.5 * (zz(k, cell2) + zz(k, cell1)));
        pgrad = cqu(k, iEdge) * 0.5 * c2 * (exner(k, cell1) + exner(k, cell2)) * pgrad;
        pgrad = pgrad + 0.5 * m.zxu(k, iEdge) * gravity * (rho_pp(k, cell1) + rho_pp(k, cell2));
        ru_p(k, iEdge) = ru_p(k, iEdge) + dts * (tend_ru(k, iEdge) - (1.0 - m.specZoneMaskEdge(iEdge)) * pgrad);
      }
      for (int k = 1; k <= nVertLevels; ++k) ruAvg(k, iEdge) = ruAvg(k, iEdge) + ru_p(k, iEdge);
    });
  } else {  // :3924-3950
    parallel_for(1, m.nEdges, [=] CPPMPAS_LAMBDA(int iEdge) {
      for (int k = 1; k <= nVertLevels; ++k) ru_p(k, iEdge) = dts * tend_ru(k, iEdge);
      for (int k = 1; k <= nVertLevels; ++k) ruAvg(k, iEdge) = ru_p(k, iEdge);
    });
  }
  if (small_step == 1) {  // :3954-3972
    parallel_for(1, m.nCells, [=] CPPMPAS_LAMBDA(int iCell) {
      for (int k = 1; k <= nVertLevels; ++k) rtheta_pp_old(k, iCell) = 0.0;
    });
  } else {
    parallel_for(1, m.nCells, [=] CPPMPAS_LAMBDA(int iCell) {
      for (int k = 1; k <= nVertLevels; ++k) rtheta_pp_old(k, iCell) = rtheta_pp(k, iCell);
    });
  }

  parallel_for(1, m.nCellsSolve(), [=, m = view(m)] CPPMPAS_LAMBDA(int iCell) {  // :3978-4105
    Column<real> ts(nVertLevels + 1), rs(nVertLevels + 1);
    if (small_step == 1) {
      for (int k = 1; k <= nVertLevels; ++k) {
        wwAvg(k, iCell) = 0.0;
        rho_pp(k, iCell) = 0.0;
        rtheta_pp(k, iCell) = 0.0;
        rw_p(k, iCell) = 0.0;
      }
      wwAvg(nVertLevels + 1, iCell) = 0.0;
      rw_p(nVertLevels + 1, iCell) = 0.0;
    }
    if (m.specZoneMaskCell(iCell) != 0.0)
      return kernel_fail("atm_advance_acoustic_step: specified-zone cells (TI:4093-4103) are not ported");
    for (int k = 1; k <= nVertLevels; ++k) {
      ts[k] = 0.0;
      rs[k] = 0.0;
    }
    for (int i = 1; i <= m.nEdgesOnCell(iCell); ++i) {  // :4001-4012
      const int iEdge = m.edgesOnCell(i, iCell);
      const int cell1 = m.cellsOnEdge(1, iEdge);
      const int cell2 = m.cellsOnEdge(2, iEdge);
      for (int k = 1; k <= nVertLevels; ++k) {
        const real flux = m.edgesOnCell_sign(i, iCell) * dts * m.dvEdge(iEdge) * ru_p(k, iEdge) * m.invAreaCell(iCell);
        rs[k] = rs[k] - flux;
        ts[k] = ts[k] - flux * 0.5 * (theta_m(k, cell2) + theta_m(k, cell1));
      }
    }
    for (int k = 1; k <= nVertLevels; ++k) {  // :4020-4026
      rs[k] = rho_pp(k, iCell) + dts * tend_rho(k, iCell) + rs[k] -
              dts * cofrz(k) * (ewm(k + 1) * rw_p(k + 1, iCell) - ewm(k) * rw_p(k, iCell));
      ts[k] = rtheta_pp(k, iCell) + dts * tend_rt(k, iCell) + ts[k] -
              dts * rdzw(k) * (ewm(k + 1) * coftz(k + 1, iCell) * rw_p(k + 1, iCell) - ewm(k) * coftz(k, iCell) * rw_p(k, iCell));
    }
    for (int k = 2; k <= nVertLevels; ++k) wwAvg(k, iCell) = wwAvg(k, iCell) + ewm(k) * rw_p(k, iCell);  // :4030-4032
    for (int k = 2; k <= nVertLevels; ++k) {  // :4036-4047
      rw_p(k, iCell) =
          rw_p(k, iCell) +
          dts * (tend_rw(k, iCell) -
                 cofwz(k, iCell) * ((etp(k) * zz(k, iCell) * ts[k] - etp(k - 1) * zz(k - 1, iCell) * ts[k - 1]) +
                                    (etm(k) * zz(k, iCell) * rtheta_pp(k, iCell) - etm(k - 1) * zz(k - 1, iCell) * rtheta_pp(k - 1, iCell))) -
                 cofwr(k, iCell) * ((etp(k) * rs[k] + etp(k - 1) * rs[k - 1]) +
                                    (etm(k) * rho_pp(k, iCell) + etm(k - 1) * rho_pp(k - 1, iCell))) +
                 cofwt(k, iCell) * (etp(k) * ts[k] + etm(k) * rtheta_pp(k, iCell)) +
                 cofwt(k - 1, iCell) * (etp(k - 1) * ts[k - 1] + etm(k - 1) * rtheta_pp(k - 1, iCell)));
    }
    for (int k = 2; k <= nVertLevels; ++k)  // :4053-4055
      rw_p(k, iCell) = (rw_p(k, iCell) - (dts * dts) * a_tri(k, iCell) * rw_p(k - 1, iCell)) * alpha_tri(k, iCell);
    for (int k = nVertLevels; k >= 1; --k)  // :4059-4061
      rw_p(k, iCell) = rw_p(k, iCell) - gamma_tri(k, iCell) * rw_p(k + 1, iCell);
    for (int k = 2; k <= nVertLevels; ++k) {  // :4067-4073
      rw_p(k, iCell) = (rw_p(k, iCell) + (rw_save(k, iCell) - rw(k, iCell)) -
                        dts * dss(k, iCell) * (fzm(k) * zz(k, iCell) + fzp(k) * zz(k - 1, iCell)) *
                            (fzm(k) * rho_zz(k, iCell) + fzp(k) * rho_zz(k - 1, iCell)) * w(k, iCell)) /
                           (1.0 + dts * dss(k, iCell)) -
                       (rw_save(k, iCell) - rw(k, iCell));
    }
    for (int k = 2; k <= nVertLevels; ++k) wwAvg(k, iCell) = wwAvg(k, iCell) + ewp(k) * rw_p(k, iCell);  // :4078-4080
    for (int k = 1; k <= nVertLevels; ++k) {  // :4086-4091
      rho_pp(k, iCell) = rs[k] - dts * cofrz(k) * (ewp(k + 1) * rw_p(k + 1, iCell) - ewp(k) * rw_p(k, iCell));
      rtheta_pp(k, iCell) = ts[k] - dts * rdzw(k) * (ewp(k + 1) * coftz(k + 1, iCell) * rw_p(k + 1, iCell) -
                                                     ewp(k) * coftz(k, iCell) * rw_p(k, iCell));
    }
  });
}

// TI:4112-4193 (theta_m of time level 1)
void atm_divergence_damping_3d(const DycoreConfig& cfg, const Mesh& m, const State& s, DycoreFields& f, real dts) {
  const int nVertLevels = m.nVertLevels;
  const auto theta_m = s.theta_m.cview();
  const real rdts = 1.0 / dts;                                                   // :4153
  const real coef_divdamp = 2.0 * cfg.config_smdiv * cfg.config_len_disp * rdts;  // :4154
  parallel_for(1, m.nEdges, [=, m = view(m), f = view(f)] CPPMPAS_LAMBDA(int iEdge) {  // :4162-4188
    const int cell1 = m.cellsOnEdge(1, iEdge);
    const int cell2 = m.cellsOnEdge(2, iEdge);
    for (int k = 1; k <= nVertLevels; ++k) {
      const real divCell1 = -(f.rtheta_pp(k, cell1) - f.rtheta_pp_old(k, cell1));
      const real divCell2 = -(f.rtheta_pp(k, cell2) - f.rtheta_pp_old(k, cell2));
      f.ru_p(k, iEdge) = f.ru_p(k, iEdge) + coef_divdamp * (divCell2 - divCell1) * (1.0 - m.specZoneMaskEdge(iEdge)) /
                                                (theta_m(k, cell1) + theta_m(k, cell2));
    }
  });
}

// TI:4196-4527 (state of time level 2; bdyMaskCell = 0 <= nRelaxZone for every cell)
void atm_recover_large_step_variables(const Mesh& m, State& s, Diag& d, DycoreFields& f, real dt, int ns,
                                      int rk_step) {
  (void)s;
  const int nVertLevels = m.nVertLevels;
  auto rho_zz = f.rho_zz2.view(), w = f.w2.view(), theta_m = f.theta_m2.view(), u = f.u2.view();
  auto rw = d.rw.view(), rtheta_p = d.rtheta_p.view(), rho_p = d.rho_p.view(), ru = d.ru.view(),
       exner = d.exner.view(), pressure_p = d.pressure_p.view();
  auto wwAvg = f.wwAvg.view(), ruAvg = f.ruAvg.view();
  const auto zz = m.zz.cview(), fzm = m.fzm.cview(), fzp = m.fzp.cview();
  const auto rtheta_base = d.rtheta_base.cview(), exner_base = d.exner_base.cview(), rho_base = s.rho_base.cview();
  constexpr real rcv = rgas / (cp - rgas);  // :4382
  constexpr real p0 = 1.0e+05;              // :4383
  constexpr real rgas_p0 = rgas / p0;       // (rgas/p0) of :4439, folded
  const real cf1 = m.cf1(1), cf2 = m.cf2(1), cf3 = m.cf3(1);

  parallel_for(1, nVertLevels, [=, m = view(m)] CPPMPAS_LAMBDA(int k) {
    rho_zz(k, m.nCells + 1) = 1.0;  // :4389-4391
  });
  const real invNs = 1 / static_cast<real>(ns);                         // :4397

  parallel_for(1, m.nCells, [=, f = view(f)] CPPMPAS_LAMBDA(int iCell) {  // :4401-4427
    for (int k = 1; k <= nVertLevels; ++k) {
      rho_p(k, iCell) = f.rho_p_save(k, iCell) + f.rho_pp(k, iCell);
      rho_zz(k, iCell) = rho_p(k, iCell) + rho_base(k, iCell);
    }
    rw(1, iCell) = 0.0;
    w(1, iCell) = 0.0;
    for (int k = 2; k <= nVertLevels; ++k) {
      wwAvg(k, iCell) = f.rw_save(k, iCell) + (wwAvg(k, iCell) * invNs);
      rw(k, iCell) = f.rw_save(k, iCell) + f.rw_p(k, iCell);
      w(k, iCell) = rw(k, iCell) / (fzm(k) * zz(k, iCell) + fzp(k) * zz(k - 1, iCell));
    }
    rw(nVertLevels + 1, iCell) = 0.0;
    w(nVertLevels + 1, iCell) = 0.0;
  });

  if (rk_step == 3) {  // :4430-4445
    parallel_for(1, m.nCells, [=, f = view(f)] CPPMPAS_LAMBDA(int iCell) {
      auto lane = [&](int k, bool vector) {
        rtheta_p(k, iCell) = f.rtheta_p_save(k, iCell) + f.rtheta_pp(k, iCell) -
                             dt * rho_zz(k, iCell) * f.rt_diabatic_tend(k, iCell);
        theta_m(k, iCell) = (rtheta_p(k, iCell) + rtheta_base(k, iCell)) / rho_zz(k, iCell);
        const real x = zz(k, iCell) * rgas_p0 * (rtheta_p(k, iCell) + rtheta_base(k, iCell));
        exner(k, iCell) = vector ? libm::mvec_pow(x, rcv) : libm::pow(x, rcv);
        pressure_p(k, iCell) = zz(k, iCell) * rgas *
                               (exner(k, iCell) * rtheta_p(k, iCell) + rtheta_base(k, iCell) * (exner(k, iCell) - exner_base(k, iCell)));
      };
      if (nVertLevels == 1) {
        lane(1, false);
        return;
      }
      int k = 1;
      for (; k + 1 <= nVertLevels; k += 2) {
        lane(k, true);
        lane(k + 1, true);
      }
      if (k <= nVertLevels) lane(k, false);
    });
  } else {  // :4447-4456
    parallel_for(1, m.nCells, [=, f = view(f)] CPPMPAS_LAMBDA(int iCell) {
      for (int k = 1; k <= nVertLevels; ++k) {
        rtheta_p(k, iCell) = f.rtheta_p_save(k, iCell) + f.rtheta_pp(k, iCell);
        theta_m(k, iCell) = (rtheta_p(k, iCell) + rtheta_base(k, iCell)) / rho_zz(k, iCell);
      }
    });
  }

  parallel_for(1, m.nEdges, [=, m = view(m), f = view(f)] CPPMPAS_LAMBDA(int iEdge) {  // :4467-4479
    const int cell1 = m.cellsOnEdge(1, iEdge);
    const int cell2 = m.cellsOnEdge(2, iEdge);
    for (int k = 1; k <= nVertLevels; ++k) {
      ruAvg(k, iEdge) = f.ru_save(k, iEdge) + (ruAvg(k, iEdge) * invNs);
      ru(k, iEdge) = f.ru_save(k, iEdge) + f.ru_p(k, iEdge);
      u(k, iEdge) = 2. * ru(k, iEdge) / (rho_zz(k, cell1) + rho_zz(k, cell2));
    }
  });

  parallel_for(1, m.nCells, [=, m = view(m)] CPPMPAS_LAMBDA(int iCell) {  // :4486-4523
    for (int i = 1; i <= m.nEdgesOnCell(iCell); ++i) {
      const int iEdge = m.edgesOnCell(i, iCell);
      real flux = (cf1 * ru(1, iEdge) + cf2 * ru(2, iEdge) + cf3 * ru(3, iEdge));
      w(1, iCell) = w(1, iCell) + m.edgesOnCell_sign(i, iCell) *
                                      (m.zb_cell(1, i, iCell) + sign1(flux) * m.zb3_cell(1, i, iCell)) * flux;
      for (int k = 2; k <= nVertLevels; ++k) {
        flux = (fzm(k) * ru(k, iEdge) + fzp(k) * ru(k - 1, iEdge));
        w(k, iCell) = w(k, iCell) + m.edgesOnCell_sign(i, iCell) *
                                        (m.zb_cell(k, i, iCell) + sign1(flux) * m.zb3_cell(k, i, iCell)) * flux;
      }
    }
    w(1, iCell) = w(1, iCell) / (cf1 * rho_zz(1, iCell) + cf2 * rho_zz(2, iCell) + cf3 * rho_zz(3, iCell));
    for (int k = 2; k <= nVertLevels; ++k) w(k, iCell) = w(k, iCell) / (fzm(k) * rho_zz(k, iCell) + fzp(k) * rho_zz(k - 1, iCell));
  });
}

// TI:7668-7832 (dynamics_split = 1 when config_split_dynamics_transport = .false.)
void atm_rk_dynamics_substep_finish(const Mesh& m, State& s, Diag& d, DycoreFields& f, int dynamics_substep,
                                    int dynamics_split) {
  (void)d;
  const int nVertLevels = m.nVertLevels;
  parallel_for(1, nVertLevels, [=, m = view(m), s = view(s)] CPPMPAS_LAMBDA(int k) {
    s.theta_m(k, m.nCells + 1) = 0.0;  // :7730
  });
  const real inv_dynamics_split = 1.0 / static_cast<real>(dynamics_split);  // :7733
  if (dynamics_substep < dynamics_split)
    throw Error("atm_rk_dynamics_substep_finish: dynamics substeps (TI:7735-7767) are not ported");
  if (dynamics_substep == 1) {  // :7769-7785
    parallel_for(1, m.nEdges, [=, f = view(f)] CPPMPAS_LAMBDA(int iEdge) {
      for (int k = 1; k <= nVertLevels; ++k) f.ruAvg_split(k, iEdge) = f.ruAvg(k, iEdge);
    });
    parallel_for(1, m.nCells, [=, f = view(f)] CPPMPAS_LAMBDA(int iCell) {
      for (int k = 1; k <= nVertLevels + 1; ++k) f.wwAvg_split(k, iCell) = f.wwAvg(k, iCell);
    });
  } else {
    throw Error("atm_rk_dynamics_substep_finish: dynamics substeps (TI:7786-7803) are not ported");
  }
  if (dynamics_substep == dynamics_split) {  // :7805-7829
    parallel_for(1, m.nEdges, [=, f = view(f)] CPPMPAS_LAMBDA(int iEdge) {
      for (int k = 1; k <= nVertLevels; ++k) f.ruAvg(k, iEdge) = f.ruAvg_split(k, iEdge) * inv_dynamics_split;
    });
    parallel_for(1, m.nCells, [=, f = view(f)] CPPMPAS_LAMBDA(int iCell) {
      for (int k = 1; k <= nVertLevels + 1; ++k) f.wwAvg(k, iCell) = f.wwAvg_split(k, iCell) * inv_dynamics_split;
    });
    parallel_for(1, m.nCells, [=, s = view(s), f = view(f)] CPPMPAS_LAMBDA(int iCell) {
      for (int k = 1; k <= nVertLevels; ++k) s.rho_zz(k, iCell) = f.rho_zz_old_split(k, iCell);
    });
  }
}

// ---------------------------------------------------------------------------------------------------------------
// TI:1855-2882
void atm_srk3(const DycoreConfig& cfg, const Mesh& m, State& s, Diag& d, DycoreFields& f, real dt) {
  check_dycore_supported(cfg, m, s);
  const int nVertLevels = m.nVertLevels;
  parallel_for(1, nVertLevels, [=, m = view(m), f = view(f)] CPPMPAS_LAMBDA(int k) {
    f.qtot(k, m.nCells + 1) = 0.0;  // :2023-2028
  });
  // :2031-2047 scratch physics tendencies, garbage column zeroed
  parallel_for(1, nVertLevels, [=, m = view(m), f = view(f)] CPPMPAS_LAMBDA(int k) {
    f.tend_rtheta_physics(k, m.nCells + 1) = 0.0;
    f.tend_rho_physics(k, m.nCells + 1) = 0.0;
    f.tend_ru_physics(k, m.nEdges + 1) = 0.0;
  });
  // :2053-2061 (config_split_dynamics_transport = .false.)
  const int dynamics_split = 1;
  const real dt_dynamics = dt;
  const int number_of_sub_steps = cfg.config_number_of_sub_steps;  // :2064
  real rk_timestep[4], rk_sub_timestep[4];
  int number_sub_steps[4];
  // :2080-2092 (config_time_integration_order = 2)
  rk_timestep[1] = dt_dynamics / 2.;
  rk_timestep[2] = dt_dynamics / 2.;
  rk_timestep[3] = dt_dynamics;
  rk_sub_timestep[1] = dt_dynamics / static_cast<real>(number_of_sub_steps);
  rk_sub_timestep[2] = dt_dynamics / static_cast<real>(number_of_sub_steps);
  rk_sub_timestep[3] = dt_dynamics / static_cast<real>(number_of_sub_steps);
  number_sub_steps[1] = std::max(1, number_of_sub_steps / 2);
  number_sub_steps[2] = std::max(1, number_of_sub_steps / 2);
  number_sub_steps[3] = number_of_sub_steps;
  // :2105 halo exchange: no-op on one rank
  atm_rk_integration_setup(m, s, d, f);  // :2113
  atm_compute_moist_coefficients(m, f);  // :2129
  physics_get_tend(m, s, f);             // :2141-2151 (pre/post_physics_get_tend are OpenACC-only)
  // :2168 IAU off
  for (int dynamics_substep = 1; dynamics_substep <= dynamics_split; ++dynamics_substep) {  // :2178
    atm_compute_vert_imp_coefs(m, s, d, f, rk_sub_timestep[1]);  // :2184-2193
    for (int rk_step = 1; rk_step <= 3; ++rk_step) {              // :2207
      // :2227-2245 delsq_* and dpdz are re-allocated per rk_step with the garbage column zeroed
      parallel_for(1, nVertLevels, [=, m = view(m), f = view(f)] CPPMPAS_LAMBDA(int k) {
        f.delsq_theta(k, m.nCells + 1) = 0.0;
        f.delsq_w(k, m.nCells + 1) = 0.0;
        f.delsq_divergence(k, m.nCells + 1) = 0.0;
        f.delsq_u(k, m.nEdges + 1) = 0.0;
        f.delsq_vorticity(k, m.nVertices + 1) = 0.0;
        f.dpdz(k, m.nCells + 1) = 0.0;
      });
      atm_compute_dyn_tend(cfg, m, s, d, f, rk_step, dt);  // :2250 (dt: the full time step)
      atm_set_smlstep_pert_variables(m, f);                // :2291
      for (int small_step = 1; small_step <= number_sub_steps[rk_step]; ++small_step) {  // :2359-2408
        atm_advance_acoustic_step(m, s, d, f, rk_sub_timestep[rk_step], small_step);    // :2372
        atm_divergence_damping_3d(cfg, m, s, f, rk_sub_timestep[rk_step]);              // :2401
      }
      atm_recover_large_step_variables(m, s, d, f, rk_timestep[rk_step], number_sub_steps[rk_step], rk_step);  // :2427
      // :2504 scalar advection off (check_dycore_supported)
      // :2543-2551 ke_vertex / ke_edge re-allocated, garbage column zeroed
      parallel_for(1, nVertLevels, [=, m = view(m), f = view(f)] CPPMPAS_LAMBDA(int k) {
        f.sdw.ke_vertex(k, m.nVertices + 1) = 0.0;
        f.sdw.ke_edge(k, m.nEdges + 1) = 0.0;
      });
      atm_compute_solve_diagnostics(dt, cfg.config_apvm_upwinding, m, f.rho_zz2, f.u2, d, f.sdw, 1, m.nCells, 1,
                                    m.nVertices, 1, m.nEdges, rk_step);  // :2556
    }
    atm_rk_dynamics_substep_finish(m, s, d, f, dynamics_substep, dynamics_split);  // :2645
  }
  // :2744 mpas_reconstruct of u at time level 2
  mpas_reconstruct_2d(m, f.u2, d.uReconstructX, d.uReconstructY, d.uReconstructZ, d.uReconstructZonal,
                      d.uReconstructMeridional);
  // :2799-2800 where (scalars_2 < 0.0) scalars_2 = 0.0, the whole array (garbage column included)
  const auto scalars_2 = f.scalars2.view();
  parallel_for_n(f.scalars2.size(), [=] CPPMPAS_LAMBDA(std::size_t i) {
    if (scalars_2[i] < 0.0) scalars_2[i] = 0.0;
  });
  // :2806 microphysics off; :2880 summarize_timestep writes only to the log
}

void atm_timestep(const DycoreConfig& cfg, const Mesh& m, State& s, Diag& d, DycoreFields& f, real dt,
                  const std::string& xtime_new) {
  atm_srk3(cfg, m, s, d, f, dt);  // TI:1829
  f.xtime2 = xtime_new;           // TI:1836-1842
  f.Time2 = f.Time1 + dt;         // TI:1844-1849
}

void shift_time_levels(State& s, DycoreFields& f) {
  std::swap(s.u, f.u2);
  std::swap(s.w, f.w2);
  std::swap(s.theta_m, f.theta_m2);
  std::swap(s.rho_zz, f.rho_zz2);
  std::swap(s.scalars, f.scalars2);
  std::swap(s.xtime, f.xtime2);
  std::swap(f.Time1, f.Time2);
}

void atm_compute_output_diagnostics(const Mesh& m, State& s) {
  int index_qv = 0;
  for (std::size_t k = 0; k < s.scalar_names.size(); ++k)
    if (s.scalar_names[k] == "qv") index_qv = static_cast<int>(k) + 1;
  if (index_qv == 0) throw Error("atm_compute_output_diagnostics: qv is not an active scalar");
  parallel_for(1, m.nCells, [=, m = view(m), s = view(s)] CPPMPAS_LAMBDA(int iCell) {  // mpas_atm_core.F:957-963
    for (int k = 1; k <= m.nVertLevels; ++k) {
      s.theta(k, iCell) = s.theta_m(k, iCell) / (1. + rvord * s.scalars(index_qv, k, iCell));
      s.rho(k, iCell) = s.rho_zz(k, iCell) * m.zz(k, iCell);
    }
  });
}

std::string add_seconds(const std::string& ts, long long seconds) {
  int Y, M, D, h, mi, sec;
  char tail;
  if (ts.size() != 19 || std::sscanf(ts.c_str(), "%4d-%2d-%2d_%2d:%2d:%2d%c", &Y, &M, &D, &h, &mi, &sec, &tail) != 6)
    throw Error("time stamp '" + ts + "' is not YYYY-MM-DD_hh:mm:ss");
  if (seconds < 0) throw Error("add_seconds: negative interval");
  auto leap = [](int y) { return (y % 4 == 0 && y % 100 != 0) || y % 400 == 0; };
  auto mdays = [&](int y, int mo) {
    static const int n[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    return mo == 2 && leap(y) ? 29 : n[mo - 1];
  };
  if (Y < 0 || M < 1 || M > 12 || D < 1 || D > mdays(Y, M) || h > 23 || mi > 59 || sec > 59 || h < 0 || mi < 0 || sec < 0)
    throw Error("time stamp '" + ts + "' is not a valid date");
  long long t = static_cast<long long>(h) * 3600 + mi * 60 + sec + seconds;
  long long days = t / 86400;
  t %= 86400;
  while (days-- > 0) {
    if (++D > mdays(Y, M)) {
      D = 1;
      if (++M > 12) {
        M = 1;
        ++Y;
      }
    }
  }
  char buf[32];
  std::snprintf(buf, sizeof buf, "%04d-%02d-%02d_%02d:%02d:%02d", Y, M, D, static_cast<int>(t / 3600),
                static_cast<int>(t / 60 % 60), static_cast<int>(t % 60));
  return buf;
}

}  // namespace cppmpas
